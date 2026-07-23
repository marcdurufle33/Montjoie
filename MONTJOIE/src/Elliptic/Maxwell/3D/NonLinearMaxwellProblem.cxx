#ifndef MONTJOIE_FILE_NON_LINEAR_MAXWELL_PROBLEM_CXX

namespace Montjoie
{

  //! default constructor
  template<class TypeEquation>
  NonLinearMaxwellProblem<TypeEquation>::NonLinearMaxwellProblem()
    : NonLinearOpticsProblem(), Glob_mat_Kh(var_laplace)
  {
    var_laplace.SetFrequency(0.0);
    nb_components_sismo = 1;
  }
  

  //! returns the final time
  template<class TypeEquation>
  Real_wp NonLinearMaxwellProblem<TypeEquation>::GetTimeStep() const
  {
    return dt; 
  }
  
  
  //! returns the initial time
  template<class TypeEquation>
  Real_wp NonLinearMaxwellProblem<TypeEquation>::GetInitialTime() const
  {
    return initial_time; 
  }
  
  
  //! returns the final time
  template<class TypeEquation>
  Real_wp NonLinearMaxwellProblem<TypeEquation>::GetFinalTime() const
  {
    return final_time; 
  }
  
  
  //! returns the number of degrees of freedom
  template<class TypeEquation>
  int NonLinearMaxwellProblem<TypeEquation>::GetNbDof() const
  {
    return var_laplace.GetNbDof()*(1+index.omega_polarization.GetM()); 
  }
  
  
  //! modifies parameters of the problem with a line of the data file
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    NonLinearOpticsProblem::SetInputData(description_field, parameters);
    
    if ((description_field != "Frequency") && (description_field != "SismoPoint") 
        && (description_field != "FileOutputSismoPoint") && (description_field != "Polarization"))
      var_laplace.SetInputData(description_field, parameters);

    if (description_field == "SismoPointGrid")
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
            cout << "SismoPointGrid needs one parameter, for instance :" << endl;
            cout << "SismoPointGrid = xmin xmax ymin ymax nbx nby" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        Real_wp xmin(0), xmax(0), ymin(0), ymax(0); int nbx(1), nby(1);
        if (parameters(0) == "RADIAL")
          {
            xmin = 0; xmax = to_num<Real_wp>(parameters(1))/z0_adim;
            nbx = to_num<int>(parameters(2));
            ymin = 0; ymax = 0; nby = 1;
          }
        else
          {
            if (parameters.GetM() <= 5)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "SismoPointGrid needs one parameter, for instance :" << endl;
                cout << "SismoPointGrid = xmin xmax ymin ymax nbx nby" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
            xmin = to_num<Real_wp>(parameters(0))/z0_adim;
            xmax = to_num<Real_wp>(parameters(1))/z0_adim;
            ymin = to_num<Real_wp>(parameters(2))/z0_adim;
            ymax = to_num<Real_wp>(parameters(3))/z0_adim;
            nbx = to_num<int>(parameters(4));
            nby = to_num<int>(parameters(5));
          }    
        
        int nbz = grid_sismo.GetNbPointsGrid();
        GenerateGridPoints(xmin, xmax, ymin, ymax, nbx, nby, nbz,
                           grid_sismoND.GetGlobalCoordinate(), grid_sismoND.GetTheta());
        
      }
    
  }
  
  
  //! computes points of seismograms
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::GenerateGridPoints(const Real_wp& xmin, const Real_wp& xmax, const Real_wp& ymin,
                       const Real_wp& ymax, int nbx, int nby, int nbz,
                       VectR2& GlobalCoord, VectReal_wp& TetaInterp)
  {
    GlobalCoord.Reallocate(nbx*nby*nbz);
    TetaInterp.Reallocate(nbx*nby*nbz);
    
    Real_wp dx(0), dy(0);
    if (nbx > 1)
      dx = (xmax-xmin) / (nbx-1);

    if (nby > 1)
      dy = (ymax-ymin) / (nby-1);
    
    Real_wp r(0), theta(0);
    for (int i = 0; i < nbx; i++)
      for (int j = 0; j < nby; j++)
        {
          Real_wp x = xmin + i*dx;
          Real_wp y = ymin + j*dy;
          CartesianToPolar(x, y, r, theta);
          
          for (int k = 0; k < nbz; k++)
            {
              Real_wp z = grid_sismo.GetGlobalCoordinate(k);
              int num = k*(nbx*nby) + j*nbx + i;
              GlobalCoord(num).Init(r, z);
              TetaInterp(num) = theta;
            }
        }
  }
  
  
  //! avoids duplicates in seismograms
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::ForceUniquenessSismo(GridInterpolation<Dimension>& grid, IVect& proc)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
    if (nb_proc <= 1)
      return;
    
    int nb_pts = grid.GetNbGlobalPoints();
    proc.Reallocate(nb_pts);
    proc.Fill(nb_proc+1);    
    for (int i = 0; i < nb_pts; i++)
      if (grid.GetElementNumber(i) >= 0)
        proc(i) = rank_proc;
    
    IVect proc_loc(proc);
    MPI_Allreduce(proc_loc.GetData(), proc.GetData(), proc.GetM(),
                  MPI_INTEGER, MPI_MIN, var_laplace.comm_group_mode);
    
    for (int i = 0; i < nb_pts; i++)
      if (grid.GetElementNumber(i) >= 0)
        {
          if (proc(i) != rank_proc)
            grid.SetElementNumber(i, -1);
        }
#endif
  }
  
  
  //! computes right hand side and stores Dirichlet dofs
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::ComputeRightHandSide()
  {
    VectReal_wp space_source;
    // computation of F = space_source as a full Vector
    var_laplace.ComputeRightHandSide(space_source, false);
    
    // isolating source associated with dirichlet condition
    int nb_dir = var_laplace.GetNbDirichletDof();
    const IVect& DirDof = var_laplace.GetDirichletDofNumber();
    dirichlet_source.Reallocate(nb_dir);
    for (int i = 0; i < nb_dir; i++)
      {
	dirichlet_source.Index(i) = DirDof(i);
	dirichlet_source.Value(i) = space_source(DirDof(i));
	space_source(DirDof(i)) = 0;
      }
    
    // removing Dirichlet dofs
    dirichlet_dof.Reallocate(nb_dir);
    for (int i = 0; i < nb_dir; i++)
      dirichlet_dof(i) = DirDof(i);
    
    IVect num;
    var_laplace.SetDirichletDofs(0, num);
  }
  
  
  //! constructing arrays needed to complete simulation
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::ConstructAll(const string& input_file, const string& name_element, int num)
  {
    var_laplace.InitIndices(PhysicalConstant::nb_max_indices);
    
    // the data file is read
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
    ReadLinesFile(input_file, lines_data_file, var_laplace.comm_group_mode);
#else
    int nb_proc(1), rank_proc(0);
    ReadLinesFile(input_file, lines_data_file);
#endif
    
    ReadInputFile(lines_data_file, *this);
    
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (rank_proc == 0)
      if (var_laplace.DOSSIER_output.size() == 0)
        {
          var_laplace.DOSSIER_output = "[STIFFOUT]/"; 
          EcritDossier(var_laplace.DOSSIER_output, input_file, num);
        }
    
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      MPI_Bcast_string(var_laplace.DOSSIER_output, 0, var_laplace.comm_group_mode);
#endif    

    // adimensionalization of parameters contained in var_laplace
    var_laplace.PerformAdimensionalization(t0_adim, z0_adim);

    Real_wp eps = 1e-6*dt;
    for (int i = 0; i < var_laplace.output_grid_param.GetM(); i++)
      var_laplace.output_grid_param(i).SetThresholdTime(eps);
    
    // constructing the mesh and finite elements
    var_laplace.mesh_num.SetSameNumberPeriodicDofs();
    var_laplace.ComputeMeshAndFiniteElement(name_element);
    
    if (nb_proc == 1)
      var_laplace.mesh.Write("test.mesh");
    
    // damping coefficient
    //var_laplace.ref_sigma.Fill(index.alpha_damping);
    
    // computing geometrical quantities in order to compute stiffness matrix
    var_laplace.ComputeMassMatrix();
    var_laplace.ComputeQuasiPeriodicPhase();
    
    var_laplace.PerformOtherInitializations();
    
    // calcul du second membre en espace et on enleve les noeuds de Dirichlet
    this->ComputeRightHandSide();

    // getting z0, zmax and the number of points along z-axis
    Real_wp z0 = var_laplace.GetYmin();
    Real_wp zmax = var_laplace.GetYmax();
    int Ne = -1;
    for (int i = 0; i < var_laplace.mesh.GetNbVertices(); i++)
      if (var_laplace.IsVertexOnAxis(i))
        Ne++;
    
    dz_translation = (zmax-z0) / Ne;
    nb_iterations_translation = 0;
    DISP(dz_translation);
    
    // storing original positions of vertices
    OriginalVertex.Reallocate(var_laplace.mesh.GetNbVertices());
    for (int i = 0; i < var_laplace.mesh.GetNbVertices(); i++)
      OriginalVertex(i) = var_laplace.mesh.Vertex(i);
    
    // on ne traite que le premier mode
    var_laplace.SetCurrentModeNumber(0);
    
    // computing stiffness matrix
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(0.0);
    nat_mat.SetCoefDamping(0.0);
    nat_mat.SetCoefStiffness(1.0);

    abort();
    //var_laplace.AddMatrixWithBC(Glob_mat_Kh, nat_mat);
    Glob_mat_Kh.CompressMatrix();
    
    // computing mass matrix mat_mass and damping matrix
    var_laplace.GetMassMatrix(mat_mass);
    //var_laplace.GetDampingMatrix(mat_sigma);
    
    for (int i = 0; i < dirichlet_dof.GetM(); i++)
      mat_mass(dirichlet_dof(i)) = 1.0;
    
    // inversion of mass matrix
    invMat_mass.Reallocate(mat_mass.GetM());
    for (int i = 0; i < mat_mass.GetM(); i++)
      invMat_mass(i) = 1.0/mat_mass(i);
    
    mat_mass.Clear();
    
    // velocity of translation
    if (velocity_translate == 0.0)
      velocity_translate = real(index.GetGroupVelocity(omega_L));
  }
  
  
  //! allocation of arrays needed for the time scheme
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>::InitTimeIterations()
  {
    // neighboring processors ?
    neighboring_proc.Clear();
    for (int i = 0; i < var_laplace.MatchingNumber_Subdomain.GetM(); i++)
      if (var_laplace.MatchingDofOrig_Subdomain(i).GetM() > 0)
        neighboring_proc.PushBack(var_laplace.MatchingNumber_Subdomain(i));
    
    // initializing seismogramms
    grid_sismoND.InitInterpolationGrid(var_laplace.mesh);
    grid_sismoND.LocalizePoints(var_laplace.mesh);
    ForceUniquenessSismo(grid_sismoND, previous_proc_sismo);
    
    if (nb_points_sismo <= 0)
      {
	nb_points_sismo = toInteger(round((2.0*Tmax)/dt_sismo))+2;
	// searching the closest power of 2
	int k = 2;
	while (k < nb_points_sismo)
	  k *= 2;
	
	nb_points_sismo = k;
      }
    
    int nb_pts = grid_sismoND.GetNbGlobalPoints();
    time_begin_sismo.Reallocate(nb_pts);
    En_sismo.Reallocate(nb_pts);
    last_sismo_point_number.Reallocate(nb_pts);
    nb_components_sismo = 1;
    last_sismo_point_number.Fill(-1);
    time_begin_sismo.Fill(0);
    
    // initialization of snapshots
    var_laplace.InitOutput(Real_wp(0));
  }  
  
  
  //! matrix vector product with stiffness matrixKh, V = Kh U
  template<class TypeEquation> template<class Vector1>
  void NonLinearMaxwellProblem<TypeEquation>
  ::MltStiffness(const Vector1& U, Vector1& V)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::STIFFNESS);
    
    this->Glob_mat_Kh.MltVector(U, V);
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::STIFFNESS);
  }  
  
  
  //! computes Px = f(tn, X) where the evolution system is equal to dU/dt = f(t, U)
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Px,
                     bool invert, bool source)
  {
    if (!var_laplace.FirstOrderFormulation())
      {
        cout << "Scheme not compatible with first-order formulation" << endl;
        abort();
      }
    
    // X contains (E, P_k, H, Q_k)
    abort();
    
    // getting E from displacement D
    //GetElectricFieldFromDisplacement(X, E_current);
    
    // enforcing Dirichlet condition on the left side
    //if (var_laplace.condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
    //E_current(0) = GetPulseTime(tn);
    
    // extracting H
    /*int offset = nb_points_z*(1+index.omega_polarization.GetM());
    for (int i = 0; i < H_current.GetM(); i++)
      H_current(i) = X(offset + i);

    Real_wp coef = index.c0;
    MltStiffness(coef, SeldonNoTrans, H_current, ProdTmp);
    MltStiffness(coef, SeldonTrans, E_current, H_current);

    // forming Px = | -c0 \int H d phi / dz - \sigma \int E phi
    //              | omega_k Q_k
    //              | c0 \int dE/dz psi 
    //              | omega_k ( \alpha_k E - P_k)
    for (int i = 0; i < nb_points_z; i++)
      Px(i) = (-ProdTmp(i) - mat_sigma(i)*E_current(i))*invMat_mass(i);
    
    offset = nb_points_z*(1+index.omega_polarization.GetM());
    for (int i = 0; i < H_current.GetM(); i++)
      Px(offset+i) = H_current(i)*invMat_mass(nb_points_z+i);
    
    int offsetQ = offset + H_current.GetM();
    offset = nb_points_z;
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp wk = index.omega_polarization(k);
        Real_wp sk = 0;
        if (index.sigma_polarization.GetM() > 0)
          sk = index.sigma_polarization(k)*wk;
        
        Real_wp alpha = index.alpha_polarization(k);
        for (int i = 0; i < nb_points_z; i++)
          {
            Px(offset+i) = wk*X(offsetQ+i);
            Px(offsetQ+i) = wk*(alpha*E_current(i) - X(offset+i) - sk*X(offsetQ+i));
          }
        
        offset += nb_points_z;
        offsetQ += nb_points_z;
      }

    if (var_laplace.condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
    Px(0) = 0; */
  }
    
  
  //! computes ProdX = f(tn, X, Xp) where the evolution system is equal to 
  //! d^2 U/dt^2 = f(t, U, dU/dt)
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::EvaluateFunctionS(const Real_wp& tn, const VectReal_wp& X, const VectReal_wp& Xp,
		      VectReal_wp& ProdX, bool invert, bool source)
  {
    if (var_laplace.FirstOrderFormulation())
      {
        cout << "Scheme not compatible with second-order formulation" << endl;
        abort();
      }
    
    // X contains (E, P_k)
    
    // getting E from displacement D
    VectReal_wp E_current(var_laplace.GetNbDof()), ProdTmp(var_laplace.GetNbDof());
    GetElectricFieldFromDisplacement(tn, X, E_current);
    
    // enforcing Dirichlet condition
    Real_wp pulse_time;
    GetPulseTime(tn, pulse_time); //DISP(omega_L); DISP(tn); DISP(pulse_time);
    for (int i = 0; i < dirichlet_source.GetM(); i++)
      E_current(dirichlet_source.Index(i)) = dirichlet_source.Value(i)*pulse_time;
    
    // product with stiffness matrix
    Real_wp coef = index.c0*index.c0; //DISP(coef);
    MltStiffness(E_current, ProdTmp);
    //E_current.Write("E.dat"); ProdTmp.Write("ProdE.dat");
    //Glob_mat_Kh.WriteText("KhTest.dat");
    //invMat_mass.WriteText("invMhTest.dat");
    
    // forming result = | - \int d E/dz d \varphi_i dz
    //                  | \omega_k^2 (alpha_k E - P_k)
    #pragma omp parallel for
    for (int i = 0; i < var_laplace.GetNbDof(); i++)
      ProdX(i) = -coef*ProdTmp(i)*invMat_mass(i);
    
    // damping terms
    if (mat_sigma.GetM() > 0)
      {
	VectReal_wp Eprime(var_laplace.GetNbDof());
	abort();
        //GetDeriveElectricFieldFromDisplacement(tn, Xp, E_current, Eprime);
        
        #pragma omp parallel for
	for (int i = 0; i < var_laplace.GetNbDof(); i++)
	  ProdX(i) -= mat_sigma(i)*invMat_mass(i)*Eprime(i);
      }
    
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::PROD);
    
    // linear polarization
    int offset = var_laplace.GetNbDof();
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp wk2 = square(index.omega_polarization(k));
        Real_wp alpha = index.alpha_polarization(k);
        if (index.sigma_polarization.GetM() > 0)
          {
            Real_wp sk = index.sigma_polarization(k);
            #pragma omp parallel for
            for (int i = 0; i < var_laplace.GetNbDof(); i++)
              ProdX(offset + i) = wk2*(alpha*E_current(i) - X(offset+i) - sk*Xp(offset+i));
          }
        else
          #pragma omp parallel for
          for (int i = 0; i < var_laplace.GetNbDof(); i++)
            ProdX(offset + i) = wk2*(alpha*E_current(i) - X(offset+i));
        
        offset += var_laplace.GetNbDof();
      }

    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::PROD);

    // canceling Dirichlet dofs
    for (int i = 0; i < dirichlet_source.GetM(); i++)
      ProdX(dirichlet_source.Index(i)) = 0.0;
  }
  
  
  //! returns the two-norm of a vector
  template<class TypeEquation> template<class T>
  Real_wp NonLinearMaxwellProblem<TypeEquation>
  ::GetNorm2Vector(Vector<T>& x)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
    if (nb_proc > 1)
      {
        DistributedVector<T> xdis(var_laplace.GetOverlappedDofNumber(), var_laplace.comm_group_mode);
        
        xdis.SetData(x);
        Real_wp res = Norm2(xdis);
        xdis.Nullify();
        
        return res;
      }
    else
      return Norm2(x);
#else
    return Norm2(x);
#endif
  }

  
  //! advancing time-scheme, computes E^n+1 from E^n and E^n-1
  template<class TypeEquation> template<class T>
  void NonLinearMaxwellProblem<TypeEquation>
  ::AdvanceScheme(int nt, const Real_wp& t, All_TimeScheme<T>& time_scheme,
                  Vector<T>& E_current)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    // replaces E^n by E^n+1
    time_scheme.Advance(t, nt, *this);
        
    // translates mesh if needed
    UpdateMesh(nt, t, time_scheme);
    
    // retrieving E from displacement
    GetElectricFieldFromDisplacement(t, time_scheme.GetIterate(), E_current);
    
    // enforcing Dirichlet condition
    T pulse_time;
    GetPulseTime(t+dt, pulse_time);
    for (int i = 0; i < dirichlet_source.GetM(); i++)
      E_current(dirichlet_source.Index(i)) = dirichlet_source.Value(i)*pulse_time;
    
    if (print_level >= 4)
      if (nt%100 == 0)
        {
          cout.setf(ios::scientific);
          Real_wp normE = GetNorm2Vector(E_current)*E0_adim;
          if (rank_proc == 0)
            {
              cout << "At time t = " << t*t0_adim << endl;
              cout << " ||En || = " << normE << endl;          
            }
        }
    
    if (print_level >= 2)
      {
        if (nt%1000 == 0)
          if (rank_proc == 0)
            {
              cout.setf(ios::fixed);
              cout << "At time t = " << t*t0_adim << endl;
              cout << "Time spent in stiffness product : "
                   << glob_chrono.GetSeconds(VirtualTimer::STIFFNESS) << endl;
              cout << "Time spent in inversion of displacement : " 
                   << glob_chrono.GetSeconds(VirtualTimer::FLUX) << endl;
              //cout << "Time spent to translate back solution : "
              //     << glob_chrono.GetSeconds(VirtualTimer::EXTRAPOL) << endl;
              cout << "Time spent in polarization scheme : "
                   << glob_chrono.GetSeconds(VirtualTimer::PROD) << endl;
              //cout << "Time spent in other parts of the scheme : "
              //     << glob_chrono.GetSeconds(VirtualTimer::SCHEME) << endl;
              cout << "Time spent in outputs : "
                   << glob_chrono.GetSeconds(VirtualTimer::OUTPUT) << endl;
              cout << "Time spent in the global simulation : "
                   << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
              cout.unsetf(ios::fixed);
              cout.setf(ios::scientific);
            }
      }
    else
      {
        cout.setf(ios::fixed);
        if (rank_proc == 0)
          if (nt%10000 == 0)
            cout << "At time t = " << t*t0_adim << endl;
        
        cout.unsetf(ios::fixed);
        cout.setf(ios::scientific);
      }    
    
    
    //int test_input; cout << "we wait" << endl; cin >> test_input;
  }
  
  
  //! translates the mesh if needed
  /*!
    The mesh is translated if needed in order to follow the solution
  */  
  template<class TypeEquation> template<class T>
  void NonLinearMaxwellProblem<TypeEquation>
  ::UpdateMesh(int nt, const Real_wp& t,
               All_TimeScheme<T>& time_scheme)
  {
  }  
  
  
  //! translating mesh of N elements and Ndof dofs
  template<class TypeEquation> template<class T>
  void NonLinearMaxwellProblem<TypeEquation>
  ::TranslateMesh(int N, int Ndof, int Ndof_H,
                  int nb_proc, int offset_elt, int nb_elt_left,
                  All_TimeScheme<T>& time_scheme)
  {
  }
  
  
  //! translates the field En of n dofs
  template<class TypeEquation> template<class T>
  void NonLinearMaxwellProblem<TypeEquation>
  ::TranslateIterate(Vector<T>& En, int Ndof, int Ndof_H,
                     int nb_proc, int offset_elt, int nb_elt_left)
  {
  }
  

  //! finalizing output of seismogramms if needed
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>::CloseBuffers()
  {
    for (int i = 0; i < En_sismo.GetM(); i++)
      if (En_sismo(i).GetM() > 0)
	WriteSismo(i);
  }
  
  
  //! writing seismogramm i
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>::WriteSismo(int i)
  {
    string name_file = GetBaseString(file_output_sismo) + "_P" + to_str(i);
    
    int N = last_sismo_point_number(i)+1;
    Real_wp t0 = -Tmax + time_begin_sismo(i);
    Real_wp t1 = -Tmax + time_begin_sismo(i) + (N-1)*dt_sismo;
    VectReal_wp u(N);
    string DOSSIER = var_laplace.DOSSIER_output;
    int nb_u = nb_components_sismo;
    for (int m = 0; m < nb_u; m++)
      {
        string ext;
        if (m == 0)
          ext = "_timeX.dat";
        else if (m == 1)
          ext = "_timeY.dat";
        else
          ext = "_timeZ.dat";
        
        // extracting component m
        for (int k = 0; k < N; k++)
          u(k) = E0_adim*En_sismo(i)(nb_u*k+m);
        
        // En is written by using the "loadND" format
        GridInterpolationFull<Dimension1> grid;
        int type = GridInterpolationFull<Dimension1>::LINE;
        grid.Init(type, t0*t0_adim, t1*t0_adim, N);
        
        int double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
        bool ascii = false;
        WriteMatlab(u, grid, DOSSIER+name_file+ext, double_prec, ascii);
    
        // computing Fourier transform of the signal
        if (omega_fourier.GetM() > 0)
          {    
            if (m == 0)
              ext = "_harmonicX.dat";
            else if (m == 1)
              ext = "_harmonicY.dat";
            else
              ext = "_harmonicZ.dat";
            
            N = En_sismo(i).GetM()/nb_u;
            u.Reallocate(N);
            for (int k = 0; k < N; k++)
              u(k) = E0_adim*t0_adim*En_sismo(i)(nb_u*k+m);
            
            if (fft_time.GetNbPoints() != N)
              fft_time.Init(N);

            // step for pulsation
            int facteur = N/nb_points_sismo;
            if (N%nb_points_sismo != 0)
              {
                DISP(N); DISP(nb_points_sismo); abort();
              }
            
            Real_wp domega = 2.0*pi_wp/(dt*nb_points_sismo);
            
            int n0 = toInteger(floor(omega_fourier(0)/domega));
            int n1 = toInteger(ceil(omega_fourier(omega_fourier.GetM()-1)/domega));
            if (n0 < 0)
              n0 = 0;
            
            Real_wp omega0 = domega*n0/t0_adim;
            Real_wp omega1 = domega*n1/t0_adim;
            
            if ((n1 > N/2) || (n0 > n1))
              continue;
    
            int istep = 1;
            n0 *= facteur;
            n1 *= facteur;
            istep *= facteur;
        
            VectComplex_wp uhat(N/2+1);
            fft_time.ApplyForward(u, uhat);
    
            Real_wp coef_fft_u = 2.0*dt/sqrt(2.0*pi_wp);
            VectComplex_wp usub((n1-n0)/istep+1);
            int num = 0;
            for (int k = n0; k <= n1; k += istep)
              {
                usub(num) = coef_fft_u*conj(uhat(k));
                num++;
              }
            
            grid.Init(type, omega0, omega1, usub.GetM());
            WriteMatlab(usub, grid, DOSSIER+name_file+ext, double_prec, ascii);
          }
      }
  }
  
  
  //! interpolation with basis functions
  template<class TypeEquation> template<class T>
  TinyVector<T, 1> NonLinearMaxwellProblem<TypeEquation>
  ::GetInterpolate(const Vector<T>& E, int ne, const R_N& pt_loc, const R_N& pt_glob, const MatrixN_N& dfjm1)
  {
    Vector<TinyVector<Real_wp, 1> > phi;
    var_laplace.GetReferenceElement(ne).ComputeValuesPhi(pt_loc, phi, dfjm1,
							 var_laplace.mesh_num, ne);
    
    TinyVector<T, 1> val_En;
    IVect Nodle = var_laplace.GetDofNumberOnElement(ne);
    for (int j = 0; j < phi.GetM(); j++)
      Add(E(Nodle(j)), phi(j), val_En);        
    
    return val_En;
  }
  
  
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::CopyValueSismo(const TinyVector<Real_wp, 1>& val, VectReal_wp& En, int offset)
  {
    En(offset) = val(0);
  }
  
  
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::CopyValueSismo(const TinyVector<Complex_wp, 1>& val, VectReal_wp& En, int offset)
  {
    En(offset) = real(val(0));
    En(offset+1) = imag(val(0));
  }
  
  
  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::CopyValueSismo(const TinyVector<Real_wp, 3>& val, VectReal_wp& En, int offset)
  {
    En(offset) = val(0);
    En(offset+1) = val(1);
    En(offset+2) = val(2);
  }
    

  template<class TypeEquation>
  void NonLinearMaxwellProblem<TypeEquation>
  ::SetDirichletCondition(const Real_wp& t, int n, Vector<Real_wp>& Y, Real_wp alpha)
  {
  }

  
  //! snapshot is computed and written in the output file
  template<class TypeEquation> template<class T>
  void NonLinearMaxwellProblem<TypeEquation>
  ::WriteSnapshot(int nt, const Real_wp& t, Vector<T>& En)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::OUTPUT);
    
    // snapshot on interpolation grids
    var_laplace.WriteSnapshot(nt, t, En, 1);
    
    // seismogramms
    int nb_u = sizeof(T)/sizeof(Real_wp);
    nb_components_sismo = nb_u;
    TinyVector<T, 1> val_En;
    for (int i = 0; i < En_sismo.GetM(); i++)
      {
	if (grid_sismoND.GetElementNumber(i) >= 0)	
	  {
	    // evaluating En at the given point
	    val_En = GetInterpolate(En, grid_sismoND.GetElementNumber(i),
                                    grid_sismoND.GetLocalCoordinate(i),
                                    grid_sismoND.GetGlobalCoordinate(i),
                                    grid_sismoND.GetDFjm1(i));
	    
	    bool write_sismo = true;
	    if (En_sismo(i).GetM() <= 0)
	      {
		/* if (Norm2(val_En) < amplitude_impulse*threshold_sismo)
		  write_sismo = false;
		else
                { */
		    time_begin_sismo(i) = t;
		    En_sismo(i).Reallocate(nb_u*nb_points_sismo);
		    En_sismo(i).Fill(0);
                    //}
	      }
	    
	    if (write_sismo)
	      {
		Real_wp coef = (t-time_begin_sismo(i))/dt_sismo;
		int n = toInteger(round(coef));
		if (abs(coef-n) < 1e-5)
		  {
		    if (n*nb_u >= En_sismo(i).GetM())
		      {
			int nb_old = En_sismo(i).GetM();
			En_sismo(i).Resize(2*n*nb_u);
			for (int k = nb_old; k < 2*n*nb_u; k++)
			  En_sismo(i)(k) = 0.0;
		      }
		    
		    last_sismo_point_number(i) = n;
                    CopyValueSismo(val_En, En_sismo(i), nb_u*n);
		  }
	      }
	  }
	else
	  {
	    if (En_sismo(i).GetM() > 0)
	      {
		WriteSismo(i);
		En_sismo(i).Clear();
	      }	    
	  }
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::OUTPUT);
  }
    
};

#define MONTJOIE_FILE_NON_LINEAR_MAXWELL_PROBLEM_CXX
#endif
