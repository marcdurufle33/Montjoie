#ifndef MONTJOIE_FILE_NON_LINEAR_MAXWELL_1D_CXX

#include "NonLinearMaxwell1D_Inline.cxx"

namespace Montjoie
{

  StiffnessOperatorMaxwell* StiffnessOperatorMaxwell::GetNewObject(int r, MaxwellProblem1D& var)
  {
    switch (r)
      {
      case 2 : return new TinyStiffnessOperatorMaxwell<2>(var); break;
      case 4 : return new TinyStiffnessOperatorMaxwell<4>(var); break;
      case 6 : return new TinyStiffnessOperatorMaxwell<6>(var); break;
      case 8 : return new TinyStiffnessOperatorMaxwell<8>(var); break;
      case 10 : return new TinyStiffnessOperatorMaxwell<10>(var); break;
	//case 12 : return new TinyStiffnessOperatorMaxwell<12>(var); break;
	//case 14 : return new TinyStiffnessOperatorMaxwell<14>(var); break;
	//case 16 : return new TinyStiffnessOperatorMaxwell<16>(var); break;
      default :
	{
	  cout << "Order not implemented" << endl;
	  abort();
	}
      }
  }
  
  
  /********************************
   * TinyStiffnessOperatorMaxwell *
   ********************************/

  
  template<int r>
  TinyStiffnessOperatorMaxwell<r>::TinyStiffnessOperatorMaxwell(MaxwellProblem1D& var_)
    : var(var_), var_laplace(var_.var_laplace)
  {
    solver_real.Reallocate(10);
    solver_cplx.Reallocate(10);
  }

  
  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::GetGradientMatrix()
  {
    var_laplace.GetGradientMatrix(mat_tiny_stiff);
  }
    

  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::
  AddMatrixBand(int N, GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "AddMatrixBand not available" << endl;
    //mat_tiny_band_stiff.Reallocate(N, N);
    //var_laplace.AddMatrixFEM(mat_tiny_band_stiff, nat_mat);
  }
  

  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::
  SetFreeStorage(int N, GlobalGenericMatrix<Real_wp>& nat_mat, const Real_wp& theta)
  {
    if (theta > 0)
      {
        cout << "theta not available" << endl;
	//mat_tiny_band_stiff.Reallocate(N, N);
	//var_laplace.AddMatrixFEM(mat_tiny_band_stiff, nat_mat);
      }
    
    var_laplace.GetStiffnessMatrix(mat_tiny_stiff);
  }
  

  template<int r>
  void TinyStiffnessOperatorMaxwell<r>
  ::MltTiny(const Real_wp& alpha, const SeldonTranspose& trans, const VectComplex_wp& U, VectComplex_wp& V)
  {
    cout << "MltTiny not available" << endl;
    //Real_wp zero(0); MltAdd(alpha, trans, mat_tiny_band_stiff, U, zero, V);
  }

  
  template<int r>
  void TinyStiffnessOperatorMaxwell<r>
  ::MltFree(const Real_wp& alpha, const SeldonTranspose& trans, const VectComplex_wp& U, VectComplex_wp& V)
  {
    int offset_E = 0, offset_H = 0;
    TinyVector<Complex_wp, r+1> uloc, vloc; V.Fill(0);
    if (var_laplace.FirstOrderFormulation())
      {
	if (trans.NoTrans())
	  for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
	    {
	      TinyVectorLoop<r+1>::ExtractVector(U, offset_H, uloc);
              
	      MltTrans(mat_tiny_stiff, uloc, vloc);
	      Mlt(alpha, vloc);
              
	      TinyVectorLoop<r+1>::AddVector(vloc, offset_E, V);
	      offset_E += r; offset_H += r+1;
	    }
	else
	  for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
	    {
	      TinyVectorLoop<r+1>::ExtractVector(U, offset_E, uloc);
              
	      Mlt(mat_tiny_stiff, uloc, vloc);
	      Mlt(alpha, vloc);
              
	      TinyVectorLoop<r+1>::AddVector(vloc, offset_H, V);
	      offset_E += r; offset_H += r+1;
	    }
      }
    else
      {
	for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
	  {
	    TinyVectorLoop<r+1>::ExtractVector(U, offset_E, uloc);
            
	    Mlt(mat_tiny_stiff, uloc, vloc);
	    Mlt(alpha, vloc);
            
	    TinyVectorLoop<r+1>::AddVector(vloc, offset_E, V);
	    offset_E += r;
	  }
      }
  }
  

  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::ComputeJacobian(const VectComplex_wp& En,
							const VectReal_wp& scale)
  {
    cout << "ComputteJacobian not available" << endl;
    //mat_tiny_band_DF = mat_tiny_band_stiff;
    abort();
    // decommenter apres correction
    // ComputeJacobianMatrix(En, scale, mat_tiny_band_DF);
    //mat_tiny_band_DF.Factorize();
  }
  

  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::SolveTiny(VectComplex_wp& x_sol)
  {
    cout << "SolveTiny not available" << endl;
    //mat_tiny_band_DF.Solve(x_sol);
  }
  
  
  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::MltFirstBlock(const Complex_wp& E0, const Real_wp& alpha, int offset,
						      VectComplex_wp& b_src)
  {
    TinyVector<Complex_wp, r+1> uloc, vloc;
    uloc(0) = E0;
    
    Mlt(mat_tiny_stiff, uloc, vloc);
    Mlt(alpha, vloc);
    
    for (int i = 0; i < r+1; i++)
      b_src(offset+i) += vloc(i);
  }
  

  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::FactorizeFast(Matrix<Real_wp, General, BandedCol>& mat_band_DF, bool hg)
  {
    solver_opt.Factorize(mat_band_DF, hg);
  }
  
  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::SolveFast(VectComplex_wp& Fu)
  {
    solver_opt.Solve(Fu);
  }


  template<int r>
  void TinyStiffnessOperatorMaxwell<r>
  ::FactorizeReal(Matrix<Real_wp, General, BandedCol>& mat_band_DF, bool hg, int i)
  {
    solver_real(i).Factorize(mat_band_DF, hg);
  }


  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::SolveReal(int i, VectComplex_wp& Fu)
  {
    solver_real(i).Solve(Fu);
  }


  template<int r>
  void TinyStiffnessOperatorMaxwell<r>
  ::FactorizeComplex(Matrix<Complex_wp, General, BandedCol>& mat_band_DF, bool hg, int i)
  {
    solver_cplx(i).Factorize(mat_band_DF, hg);
  }


  template<int r>
  void TinyStiffnessOperatorMaxwell<r>::SolveComplex(int i, VectComplex_wp& Fu)
  {
    solver_cplx(i).Solve(Fu);
  }

  
  /********************
   * MaxwellProblem1D *
   ********************/

  
  //! default constructor
  MaxwellProblem1D::MaxwellProblem1D() : NonLinearOpticsProblem()
  {
    nb_points_z = 1;
    
    type_storage_stiff = FREE_STORAGE;
    type_storage_jacobian = TINY_BAND_STORAGE;

    var_laplace.SetOmega(0.0);
    
    print_level = 0;
    
    // default directory is the current directory
    DOSSIER = "./";
    var_stiff = NULL;
  }
  
  
  //! modifies parameters of the problem with a line of the data file
  void MaxwellProblem1D
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    NonLinearOpticsProblem::SetInputData(description_field, parameters);
    var_laplace.SetInputData(description_field, parameters);
    var_helmholtz.SetInputData(description_field, parameters);
    newton_solver.SetInputData(description_field, parameters);

    var_laplace.SetOmega(0.0);
    
    if (description_field == "OrderDiscretization")
      {
        // OrderDiscretization = r
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of MaxwellProblem1D" << endl;
            cout << "OrderDiscretization needs one parameter, for instance :" << endl;
            cout << "OrderDiscretization = r" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        int order = to_num<int>(parameters(0));
        var_laplace.mesh.SetOrder(order);
	var_helmholtz.mesh.SetOrder(order);
	var_stiff = StiffnessOperatorMaxwell::GetNewObject(order, *this);
      }
    else if (description_field == "StorageStiffness")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of MaxwellProblem1D" << endl;
            cout << "StorageStiffness needs one parameter, for instance :" << endl;
            cout << "StorageStiffness = Free" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        if (parameters(0) == "Band")
          type_storage_stiff = BAND_STORAGE;
        else if (parameters(0) == "TinyBand")
          type_storage_stiff = TINY_BAND_STORAGE;
        else if (parameters(0) == "Sparse")
          type_storage_stiff = SPARSE_STORAGE;
        else
          type_storage_stiff = FREE_STORAGE;
      }
    else if (description_field == "DirectoryOutput")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MaxwellProblem1D" << endl;
	    cout << "DirectoryOutput needs one parameter, for instance :" << endl;
	    cout << "DirectoryOutput = chemin" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (parameters(0) == "AUTO")
          DOSSIER.clear();
        else
          DOSSIER = parameters(0);
      }
  }
  
  
  //! constructing arrays needed to complete simulation
  void MaxwellProblem1D::ConstructAll(const string& input_file, int num)
  {
    // the data file is read
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
    ReadLinesFile(input_file, lines_data_file, var_laplace.comm_group_mode);
#else
    int rank_proc(0), nb_proc(1);
    ReadLinesFile(input_file, lines_data_file);
#endif
    
    ReadInputFile(lines_data_file, *this);
    
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (rank_proc == 0)
      if (DOSSIER.size() == 0)
        {
          DOSSIER = "[STIFFOUT]/"; 
          EcritDossier(DOSSIER, input_file, num);
        }
    
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      MPI_Bcast_string(DOSSIER, 0, var_laplace.comm_group_mode);
#endif    
    
    if (harmonic_resolution)
      {
        if (nb_proc > 1)
          {
            cout << "not implemented" << endl; abort();
          }
        
        ConstructHarmonic();
        return;
      }
    
    // adimensionalization of parameters contained in var_laplace
    VectString& mesh_data = var_laplace.GetMeshData(0);
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        mesh_data(2) = to_str(to_num<Real_wp>(mesh_data(2))/z0_adim);
        mesh_data(3) = to_str(to_num<Real_wp>(mesh_data(3))/z0_adim);
        
        int Ngrid = var_laplace.GetSectionGrid().GetNbPointsGrid();
        if (Ngrid > 1)
          {
            Real_wp xmin = var_laplace.GetSectionGrid().GetXmin()/z0_adim;
            Real_wp xmax = var_laplace.GetSectionGrid().GetXmax()/z0_adim;
            var_laplace.GetSectionGrid().Init(xmin, xmax, Ngrid);
          }
        
        var_laplace.GetOutputParameters().AdimTime(t0_adim);
        
        Real_wp delta = var_laplace.mesh.GetThicknessPML()/z0_adim;
        var_laplace.mesh.SetThicknessPML(delta);
      }

    Real_wp z0 = to_num<Real_wp>(mesh_data(2));
    Real_wp zmax = to_num<Real_wp>(mesh_data(3));
    int Ne = to_num<int>(mesh_data(1))-1;
    dz_translation = (zmax-z0) / Ne;
    nb_iterations_translation = 0;
    
    Complex_wp vg;
    vg = index.GetGroupVelocity(omega_L);
    
    // constructing the mesh and finite elements
    var_laplace.ComputeMeshAndFiniteElement(string("EDGE_LOBATTO"));
    
    OriginalVertex.Reallocate(var_laplace.mesh.GetNbElt()+1);
    if (var_laplace.mesh.GetNbElt() == 0)
      OriginalVertex(0) = zmax;
    
    for (int i = 0; i <= var_laplace.mesh.GetNbElt(); i++)
      OriginalVertex(i) = var_laplace.mesh.Vertex(i);
    
    // setting rho = 1 and mu = 1 for Laplace equation
    /*for (int ref = 0; ref < var_laplace.ref_rho.GetM(); ref++)
      {
	var_laplace.ref_rho(ref).SetConstant(Real_wp(1));
	var_laplace.ref_mu(ref).SetConstant(Real_wp(1));
	var_laplace.ref_sigma(ref).SetConstant(index.alpha_damping);
      }
    */
    
    // computing geometrical quantities in order to compute stiffness matrix
    var_laplace.ComputeMassMatrix();

    // computing stiffness matrix in the required storage
    // for second-order formulation, the stiffness matrix is equal
    // to \int d/dx (\varphi_i) d/dx (\varphi_j) dx 
    // for the fist order formulation, it is equal to
    // \int d/dx(\varphi_j) \varphi_i dx
    GlobalGenericMatrix<Real_wp> nat_mat;
    int Nscal = var_laplace.mesh.GetNbDof();
    int N = var_laplace.GetNbDof();
    int Nvec = N - Nscal;
    int m = var_laplace.mesh.GetOrder();
    nat_mat.SetCoefMass(0.0);
    nat_mat.SetCoefStiffness(1.0);
    nat_mat.SetCoefDamping(0.0);
    
    type_storage_jacobian = type_storage_stiff;

    if (var_laplace.FirstOrderFormulation())
      {
        Matrix<Real_wp, General, ArrayRowSparse> A(N, N);
        mat_sparse_stiff.Reallocate(Nscal, Nvec);
        var_laplace.AddMatrixFEM(A, nat_mat);
        //A.WriteText("Kh.dat");
        
        if (type_storage_stiff == SPARSE_STORAGE)
          {
            IVect col_interac(2*(m+1));
            VectReal_wp val_interac(2*(m+1));
            for (int i = 0; i < Nscal; i++)
              {
                int size_row = A.GetRowSize(i);
                if (size_row > col_interac.GetM())
                  {
                    abort();
                  }
                
                for (int j = 0; j < size_row; j++)
                  {
                    col_interac(j) = A.Index(i, j) - Nscal;
                    val_interac(j) = A.Value(i, j);
                  }
                
                mat_sparse_stiff.AddInteractionRow(i, size_row, col_interac, val_interac);
              }
          }
        else
          {
            type_storage_stiff = FREE_STORAGE;
	    var_stiff->GetGradientMatrix();
          }
      }
    else
      {
        if (type_storage_stiff == SPARSE_STORAGE)
          {
            mat_sparse_stiff.Reallocate(N, N);
            var_laplace.AddMatrixFEM(mat_sparse_stiff, nat_mat);
          }
        else if (type_storage_stiff == BAND_STORAGE)
          {
            mat_band_stiff.Reallocate(N, N, m, m);
            var_laplace.AddMatrixFEM(mat_band_stiff, nat_mat);
          }
        else if (type_storage_stiff == TINY_BAND_STORAGE)
          {
	    var_stiff->AddMatrixBand(N, nat_mat);
          }
        else if (type_storage_stiff == FREE_STORAGE)
          {
            type_storage_jacobian = TINY_BAND_STORAGE;
	    var_stiff->SetFreeStorage(N, nat_mat, theta);
          }
        else
          {
            cout << "not implemented " << endl;
            abort();
          }
        
      }
    
    // computing mass matrix mat_mass and damping matrix
    var_laplace.GetMassMatrix(mat_mass);
    var_laplace.GetDampingMatrix(mat_sigma);
    
    invMat_mass.Reallocate(mat_mass.GetM());
    for (int i = 0; i < mat_mass.GetM(); i++)
      invMat_mass(i) = 1.0/mat_mass(i);

    //invMat_mass.WriteText("invMh.dat");
    // vector is not cleared for ApplyMassMatrix
    // mat_mass.Clear();
    
    nb_points_z = Nscal;

    // velocity of translation
    if (velocity_translate == 0.0)
      velocity_translate = real(index.GetGroupVelocity(omega_L));
  }
  
  
  //! allocation of arrays needed for the time scheme
  int MaxwellProblem1D::SetInitialCondition(bool explicit_scheme)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif

    int Nscal = var_laplace.mesh.GetNbDof();
    int N = var_laplace.GetNbDof();
    int Nvec = N - Nscal;

    if ((nb_proc > 1) && (!explicit_scheme))
      {
        cout << "not implemented" << endl;
        abort();
      }    
    
    // initializing seismogramms
    grid_sismo.LocalizePoints(var_laplace.mesh);
    ForceUniquenessSismo(grid_sismo, previous_proc_sismo);
    if (nb_points_sismo <= 0)
      {
	nb_points_sismo = toInteger(round((2.0*Tmax)/dt_sismo))+2;
	// searching the closest power of 2
	int k = 2;
	while (k < nb_points_sismo)
	  k *= 2;
	
	nb_points_sismo = k;
      }
    
    time_begin_sismo.Reallocate(grid_sismo.GetNbPointsGrid());
    En_sismo.Reallocate(grid_sismo.GetNbPointsGrid());
    last_sismo_point_number.Reallocate(grid_sismo.GetNbPointsGrid());
    last_sismo_point_number.Fill(-1);
    time_begin_sismo.Fill(0);
    
    if (var_laplace.GetSectionGrid().GetGlobalNumberPoints() > 0)
      {
        zmin_display = var_laplace.GetSectionGrid().GetXmin();
        zmax_display = var_laplace.GetSectionGrid().GetXmax();        
      }
    else
      {
        zmin_display = -1.0;
        zmax_display = -1.0;
      }    
    
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      {
        MPI_Comm& comm = var_laplace.comm_group_mode;
        Real_wp z0(zmin_display), z1(zmax_display);
        Vector<int64_t> ztmp;
        MpiAllreduce(comm, &z0, ztmp, &zmin_display, 1, MPI_MIN);
        MpiAllreduce(comm, &z1, ztmp, &zmax_display, 1, MPI_MAX);
        
        int p(nb_points_sismo);
        MPI_Allreduce(&p, &nb_points_sismo, 1, MPI_INTEGER, MPI_MAX, comm);
      }
#endif
    
    var_laplace.GetOutputParameters().InitTime(Real_wp(0));

    int size_sys = Nscal;
    
    if (!explicit_scheme)
      {
        Real_wp coef = 0.5*this->dt*this->index.c0;
        Mlt(coef, mat_sigma);
      }
    else
      Mlt(index.c0, mat_sigma);
        
    // allocating arrays for time iterates
    if (explicit_scheme)
      {
        E_current.Reallocate(Nscal);
        E_current.Fill(0);
        ProdTmp.Reallocate(Nscal);
        ProdTmp.Fill(0);        
        
        size_sys += this->index.omega_polarization.GetM()*nb_points_z;
        if (var_laplace.FirstOrderFormulation())
          {
            size_sys += this->index.omega_polarization.GetM()*nb_points_z + Nvec;
            H_current.Reallocate(Nvec);
            H_current.Fill(0);
          }                
	
	if (this->var_raman.GetAlpha() != Real_wp(0))
	  size_sys += Nscal;
      }
    else
      {
        // home-made theta scheme
        E_current.Reallocate(Nscal); E_prev.Reallocate(Nscal); E_next.Reallocate(Nscal);
        E_current.Fill(0); E_prev.Fill(0); E_next.Fill(0);
        ProdTmp.Reallocate(Nscal); source_newton.Reallocate(Nscal);
        ProdTmp.Fill(0); source_newton.Fill(0);    
        
        mass_termEn.Reallocate(Nscal);
        mass_termEn_prev.Reallocate(Nscal);
        mass_termEn.Fill(0); mass_termEn_prev.Fill(0);
        
        int Np = index.omega_polarization.GetM();
        P_current.Reallocate(Nscal*Np); P_prev.Reallocate(Nscal*Np);
        P_next.Reallocate(Nscal*Np); rhs_polar.Reallocate(Nscal*Np);
        P_current.Fill(0); P_prev.Fill(0); P_next.Fill(0); rhs_polar.Fill(0);
        
        if (var_laplace.FirstOrderFormulation())
          {
            H_current.Reallocate(Nvec); H_next.Reallocate(Nvec);
            H_current.Fill(0); H_next.Fill(0);
          }
        
        // initializing jacobian matrix
        Vector<Real_wp> scale;
        ComputeAndFactoriseDiff(E_current, scale);
        
        // and Newton solver
        this->newton_solver.Init(*this, E_current);
      }
  
    return size_sys;
  }
  
  
  //! constructs arrays needed for resolution in time-harmonic domain
  void MaxwellProblem1D::ConstructHarmonic()
  {  
    // adimensionalization of parameters contained in var_helmholtz
    VectString& mesh_data = var_helmholtz.GetMeshData(0);
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        mesh_data(2)
          = to_str(to_num<Real_wp>(mesh_data(2))/z0_adim);
        
        mesh_data(3)
          = to_str(to_num<Real_wp>(mesh_data(3))/z0_adim);
        
        int Ngrid = var_helmholtz.GetSectionGrid().GetNbPointsGrid();
        if (Ngrid > 1)
          {
            Real_wp xmin = var_helmholtz.GetSectionGrid().GetXmin()/z0_adim;
            Real_wp xmax = var_helmholtz.GetSectionGrid().GetXmax()/z0_adim;
            var_helmholtz.GetSectionGrid().Init(xmin, xmax, Ngrid);
          }
        
        var_helmholtz.GetOutputParameters().AdimTime(t0_adim);
        
        Real_wp delta = var_helmholtz.mesh.GetThicknessPML()/z0_adim;
        var_helmholtz.mesh.SetThicknessPML(delta);
      }
    
    // constructing the mesh and finite elements
    var_helmholtz.ComputeMeshAndFiniteElement(string("EDGE_LOBATTO"));
    
    // we set rho and mu to 1
    for (int ref = 0; ref < var_helmholtz.ref_rho.GetM(); ref++)
      {
	var_helmholtz.ref_rho(ref).SetConstant(Complex_wp(1, 0));
	var_helmholtz.ref_mu(ref).SetConstant(Complex_wp(1, 0));    
      }
    
    nb_points_z = var_helmholtz.mesh.GetNbDof();
    
    // initialising spectrum to compute
    this->ConstructSpectrum(this->DOSSIER);

    // calcul des matrices de masse
    mat_mass_harmonic.Reallocate(omega.GetM());
    for (int i = 0; i < omega.GetM(); i++)
      {
	var_helmholtz.SetOmega(omega(i));
	var_helmholtz.kwave = 0.0;
	var_helmholtz.pml_damping.SetMaximumVelocity(abs(omega(i) / kwave(i)));
	if (omega(i) == 0.0)
	  var_helmholtz.pml_damping.SetMaximumVelocity(1.0);
	
	var_helmholtz.ComputeMassMatrix();
	var_helmholtz.GetMassMatrix(mat_mass_harmonic(i));
      }

    // initializing seismogramms
    grid_sismo.LocalizePoints(var_helmholtz.mesh);
    if (nb_points_sismo <= 0)
      {
	nb_points_sismo = toInteger(round((2.0*Tmax)/dt_sismo))+2;
	// searching the closest power of 2
	int k = 2;
	while (k < nb_points_sismo)
	  k *= 2;
	
	nb_points_sismo = k;
      }
  }
  
  
  //! matrix vector product with stiffness matrixKh, V = Kh U
  void MaxwellProblem1D::MltStiffness(const Real_wp& alpha, const SeldonTranspose& trans,
				      const Vector<Complex_wp>& U, Vector<Complex_wp>& V)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::STIFFNESS);
        
    Real_wp zero(0.0);
    switch (type_storage_stiff)
      {
      case SPARSE_STORAGE : MltAdd(alpha, trans, mat_sparse_stiff, U, zero, V); break;
      case BAND_STORAGE : MltAdd(alpha, trans, mat_band_stiff, U, zero, V); break;
      case TINY_BAND_STORAGE : var_stiff->MltTiny(alpha, trans, U, V); break;
      case FREE_STORAGE : var_stiff->MltFree(alpha, trans, U, V); break;
      }
    
    if (var_laplace.FirstOrderFormulation())
      {
        if (trans.NoTrans())
          var_laplace.AddDomains(V);
      }
    else
      var_laplace.AddDomains(V);
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::STIFFNESS);
  }  
  
    
  //! evaluates "right hand side" involved in Newton iteration
  /*!
    Somehow the non-linear system f(E^n+1, E^n, E^n-1)
    is written as f(E^n+1) + right_hand_side(E^n, E^n-1)
  */
  void MaxwellProblem1D::EvaluateRightHandSide(const Real_wp& t)
  {
    if (!var_laplace.FirstOrderFormulation())
      {
        if (print_level >= 2)
          glob_chrono.Start(VirtualTimer::SCHEME);
        
        // part Kh [(1-2 theta)U^n + theta U^(n-1)]
        Real_wp coef = square(index.c0*dt);
        for (int i = 0; i < nb_points_z; i++)
          ProdTmp(i) = (1.0-2.0*theta)*E_current(i) + theta*E_prev(i);
        
        if (print_level >= 2)
          glob_chrono.Stop(VirtualTimer::SCHEME);
        
        MltStiffness(coef, SeldonNoTrans, ProdTmp, source_newton);
        
        if (print_level >= 2)
          glob_chrono.Start(VirtualTimer::SCHEME);
        
        // we compute :
        // source = eps_inf Mh (-2 E^n + E^{n-1}) +gamma Mh(-2 |E^n|^2 E^n + |E^{n-1}|^2 E^(n-1) )
        //          + c^2 dt^2 Kh [(1-2 theta)E^n + theta E^(n-1)]
        
        // +damping term  = -c0 dt / 2 S E^{n-1} 
        for (int i = 0; i < nb_points_z; i++)
          {
            mass_termEn(i) = E_current(i)*(index.eps_inf + gammaNL*absSquare(E_current(i)));
            source_newton(i) += mat_mass(i)
              *(-2.0*mass_termEn(i) + mass_termEn_prev(i)) - mat_sigma(i)*E_prev(i);
          }
        
        if (print_level >= 2)
          glob_chrono.Stop(VirtualTimer::SCHEME);
      }
    else
      {
        Real_wp coef = index.c0*dt;
        // part Rh H^{n+1/2}
        MltStiffness(coef, SeldonNoTrans, H_current, ProdTmp);
        
        if (print_level >= 2)
          glob_chrono.Start(VirtualTimer::SCHEME);
        
        // we compute :
        // source = -eps_inf Mh E^n + c dt/2 Sh E^n + c dt Rh H^{n+1/2} - gamma Mh |E^n|^2 E^n	
        for (int i = 0; i < nb_points_z; i++)
          source_newton(i) = ProdTmp(i) + mat_sigma(i)*E_current(i)
            - mat_mass(i)*E_current(i)*(index.eps_inf + gammaNL*absSquare(E_current(i)));
        
        if (print_level >= 2)
          glob_chrono.Stop(VirtualTimer::SCHEME);
      }
    
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::PROD);
    
    // polarization right hand side :
    // 1 / (1/omega_k^2 + sigma_k dt/2 + theta dt^2 ) [ (2/omega_k^2 - dt^2 (1-2 theta) ) P^n 
    //  + (-1/omega_k^2 + sigma_k dt/2 - dt^2 theta) P^{n-1}
    //                                  + alpha_k dt^2 ( theta E^{n-1} + (1- 2 theta) E^n)
    //
    // we have the relation P^{n+1} = beta E^{n+1} + rhs_polar
    int offset = 0;
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp wk = index.omega_polarization(k);
        Real_wp alpha = index.alpha_polarization(k);
        Real_wp dt2 = dt*dt, invWk2 = 1.0/square(wk);
        Real_wp sk(0);
        if (index.sigma_polarization.GetM() > 0)
          sk = index.sigma_polarization(k);
        
	Real_wp coef = 1.0 / (invWk2 + 0.5*dt*sk + theta_polar*dt2);
        Real_wp coefPn = (2.0*invWk2 - dt2*(1.0-2.0*theta_polar))*coef;
        Real_wp coefPn_prev = (-invWk2 + sk*0.5*dt - dt2*theta_polar)*coef;
        Real_wp coefEn = alpha*dt2*(1.0 - 2.0*theta_polar)*coef;
        Real_wp coefEn_prev = alpha*dt2*theta_polar*coef;
        for (int i = 0; i < nb_points_z; i++)
          rhs_polar(offset + i) = coefPn*P_current(offset+i) + coefPn_prev*P_prev(offset+i)
            + coefEn*E_current(i) + coefEn_prev*E_prev(i);
        
        offset += nb_points_z;
      }
    
    // updating source_newton with polarization
    offset = 0;
    if (!var_laplace.FirstOrderFormulation())
      {
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          {
            for (int i = 0; i < nb_points_z; i++)
              source_newton(i) += mat_mass(i)*(rhs_polar(offset+i)
                                               - 2.0*P_current(offset+i) + P_prev(offset+i));
            
            offset += nb_points_z;
          }
      }
    else
      {
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          {
            for (int i = 0; i < nb_points_z; i++)
              source_newton(i) += mat_mass(i)*(rhs_polar(offset+i) - P_current(offset+i));
            
            offset += nb_points_z;
          }
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::PROD);
    
    // enforcing Dirichlet condition on the left side
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      source_newton(0) = GetPulseTime(t);
  }
  
  
  //! computation of non-linear system f to solve, scheme = f(En)
  void MaxwellProblem1D::ComputeScheme(const Vector<Complex_wp>& En, Vector<Complex_wp>& scheme)
  {
    if (harmonic_resolution)
      {
        Complex_wp zero; SetComplexZero(zero);
        
        // loop over frequencies
        Vector<Complex_wp> Ex(nb_points_z), Ey(nb_points_z),
          ProdEx(nb_points_z), ProdEy(nb_points_z);
        
        ProdEx.Fill(0); ProdEy.Fill(0);
	scheme.Fill(0);
	
	// linear part
	VectComplex_wp phase(omega.GetM());
        for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
          {
            for (int i = 0; i < nb_points_z; i++)
              {
                Ex(i) = En(nb_points_z*2*n_freq + i);
                Ey(i) = En(nb_points_z*(2*n_freq+1) + i);
              }

	    abort();
	    // a decommenter apres correction
            //Mlt(fem_harmonic_matrix(n_freq), Ex, ProdEx);
            //Mlt(fem_harmonic_matrix(n_freq), Ey, ProdEy);
            
            for (int i = 0; i < nb_points_z; i++)
              {
                scheme(nb_points_z*2*n_freq + i) = ProdEx(i);
                scheme(nb_points_z*(2*n_freq+1) + i) = ProdEy(i);
              }
	  }
        
        ProdEx.Clear(); ProdEy.Clear();
	
	// non-linear part
	VectReal_wp Ex_time(nb_points_time), Ey_time(nb_points_time); 
	Ex.Reallocate(omega.GetM()); Ey.Reallocate(omega.GetM());
	Ex.Fill(0); Ey.Fill(0); Ex_time.Fill(0); Ey_time.Fill(0);
        Vector<VectComplex_wp> VecEx, VecEy; VectComplex_wp PolarEx, PolarEy;
	for (int i = 0; i < nb_points_z; i++)
	  {
	    if (var_helmholtz.EnvelopeToCompute())
	      for (int n = 0; n < omega.GetM(); n++)
		phase(n) = exp(Iwp*kwave(n)*var_helmholtz.GetCoordinateDof(i));
	    else
	      phase.Fill(1);
	    
	    if (type_spectrum == CENTERED)
	      {		
		for (int n = 0; n < omega.GetM(); n++)
		  {
		    Ex(n) = conj(En(nb_points_z*2*n + i)*phase(n));
		    Ey(n) = conj(En(nb_points_z*(2*n+1) + i)*phase(n));
		  }
	    
		// returning back in time
		fft_time.ApplyInverse(Ex, Ex_time);
		fft_time.ApplyInverse(Ey, Ey_time);
		
		Mlt(1.0/coef_fft, Ex_time);
		Mlt(1.0/coef_fft, Ey_time);
		
		// forming |E|^2 E
		Real_wp moduleE;
		for (int j = 0; j < nb_points_time; j++)
		  {
		    moduleE = gammaNL*(square(Ex_time(j)) + square(Ey_time(j)));
		    Ex_time(j) *= moduleE;
		    Ey_time(j) *= moduleE;
		  }
		
		// going to Fourier domain
		fft_time.ApplyForward(Ex_time, Ex);
		fft_time.ApplyForward(Ey_time, Ey);
		
		Mlt(coef_fft, Ex); Mlt(coef_fft, Ey);
		
		for (int n = 0; n < omega.GetM(); n++)
		  {
		    scheme(nb_points_z*2*n + i) += mat_mass_harmonic(n)(i)*conj(Ex(n)*phase(n));
		    scheme(nb_points_z*(2*n+1) + i)
                      += mat_mass_harmonic(n)(i)*conj(Ey(n)*phase(n));
		  }
	      }
            else if (type_spectrum == ODD_FREQUENCIES)
              {
                Ex.Clear(); Ey.Clear(); Ex_time.Clear(); Ey_time.Clear();
                int N = number_odd_frequencies;
                if (i == 0)
                  {
                    VecEx.Reallocate(N); VecEy.Reallocate(N);
                    PolarEx.Reallocate(N); PolarEy.Reallocate(N);
                    for (int n = 0; n < N; n++)
                      {
                        VecEx(n).Reallocate(nb_points_time);
                        VecEy(n).Reallocate(nb_points_time);
                      }
                  }
                
                int offset = 0;
                for (int n = 0; n < N; n++)
                  {
                    for (int k = 0; k < nb_points_time; k++)
                      {
                        VecEx(n)(k) = En(nb_points_z*2*(offset+k) + i)*phase(offset+k);
                        VecEy(n)(k) = En(nb_points_z*(2*(offset+k)+1) + i)*phase(offset+k);
                      }
                    
                    fft_envelope.ApplyForward(VecEx(n));
                    fft_envelope.ApplyForward(VecEy(n));
                    
                    Mlt(1.0/coef_fft, VecEx(n));
                    Mlt(1.0/coef_fft, VecEy(n));
                    
                    offset += nb_points_time;
                  }
                
                Real_wp coef_gamma = gammaNL;
                Complex_wp Ek_dot_El, Ex_k, Ey_k, Ex_l, Ey_l;
                int Nmax = 2*N-1;
                for (int j = 0; j < nb_points_time; j++)
                  {
                    PolarEx.Fill(zero); PolarEy.Fill(zero);
                    for (int k2 = -N; k2 < N; k2++)
                      for (int l2 = k2; l2 < N; l2++)
                        {
                          int k = 2*k2+1; int l = 2*l2+1;
                          if (k2 == l2)
                            coef_gamma = 0.25*gammaNL;
                          else
                            coef_gamma = 0.5*gammaNL;
                          
                          // forming coef_gamma Ek . El
                          if (k < 0)
                            {
                              Ex_k = conj(VecEx(-k2-1)(j));
                              Ey_k = conj(VecEy(-k2-1)(j));
                            }
                          else
                            {
                              Ex_k = VecEx(k2)(j);
                              Ey_k = VecEy(k2)(j);
                            }

                          if (l < 0)
                            {
                              Ex_l = conj(VecEx(-l2-1)(j));
                              Ey_l = conj(VecEy(-l2-1)(j));
                            }
                          else
                            {
                              Ex_l = VecEx(l2)(j);
                              Ey_l = VecEy(l2)(j);
                            }
                          
                          Ek_dot_El = coef_gamma*(Ex_k*Ex_l + Ey_k*Ey_l);
                          
                          // then adding terms E_k \cdot E_l E_{n-k-l}
                          for (int n2 = 0; n2 < N; n2++)
                            {
                              int n = 2*n2+1;
                              int p = n-k-l;
                              if (p < 0)
                                {
                                  if (p >= -Nmax)
                                    {
                                      int p2 = -(p+1)/2;
                                      PolarEx(n2) += Ek_dot_El*conj(VecEx(p2)(j));
                                      PolarEy(n2) += Ek_dot_El*conj(VecEy(p2)(j));
                                    }
                                }
                              else
                                {
                                  if (p <= Nmax)
                                    {
                                      int p2 = (p-1)/2;
                                      PolarEx(n2) += Ek_dot_El*VecEx(p2)(j);
                                      PolarEy(n2) += Ek_dot_El*VecEy(p2)(j);
                                    }
                                }
                            }
                        }
                    
                    // putting non-linear polarization in VecEx, VecEy
                    for (int n = 0; n < N; n++)
                      {
                        VecEx(n)(j) = PolarEx(n);
                        VecEy(n)(j) = PolarEy(n);
                      }
                  }
                
                // going to Fourier domain
                offset = 0;
                for (int n = 0; n < N; n++)
                  {
                    fft_envelope.ApplyInverse(VecEx(n));
                    fft_envelope.ApplyInverse(VecEy(n));
		
                    Mlt(coef_fft, VecEx(n)); Mlt(coef_fft, VecEy(n));
		
                    for (int k = 0; k < nb_points_time; k++)
                      {
                        int p = offset+k;
                        scheme(nb_points_z*2*p + i) 
                          += mat_mass_harmonic(p)(i)*VecEx(n)(k)*conj(phase(p));
                        
                        scheme(nb_points_z*(2*p+1) + i)
                          += mat_mass_harmonic(p)(i)*VecEy(n)(k)*conj(phase(p));
                      }
                    
                    offset += nb_points_time;
                  }
              }
	    else if (type_spectrum == SHIFTED_SINGLE)
	      {	
		for (int n = 0; n < omega.GetM(); n++)
		  {
		    Ex(n) = En(nb_points_z*2*n + i)*phase(n);
		    Ey(n) = En(nb_points_z*(2*n+1) + i)*phase(n);
		  }
		
		// returning back in time
		fft_envelope.ApplyForward(Ex);
		fft_envelope.ApplyForward(Ey);
		
		Mlt(1.0/coef_fft, Ex);
		Mlt(1.0/coef_fft, Ey);
		
		// shifted case => 1/4 (2 |E|^2 E + E.E E*)
		Real_wp moduleE; Complex_wp EdotE;
		for (int j = 0; j < nb_points_time; j++)
		  {
		    moduleE = 0.5*gammaNL*(absSquare(Ex(j)) + absSquare(Ey(j)));
		    EdotE = 0.25*gammaNL*(Ex(j)*Ex(j) + Ey(j)*Ey(j));
		    Ex(j) = moduleE*Ex(j) + EdotE*conj(Ex(j));
		    Ey(j) = moduleE*Ey(j) + EdotE*conj(Ey(j));
		  }
		
		// going to Fourier domain
		fft_envelope.ApplyInverse(Ex);
		fft_envelope.ApplyInverse(Ey);
		
		Mlt(coef_fft, Ex); Mlt(coef_fft, Ey);
		
		for (int n = 0; n < omega.GetM(); n++)
		  {
		    scheme(nb_points_z*2*n + i) += mat_mass_harmonic(n)(i)*Ex(n)*conj(phase(n));
		    scheme(nb_points_z*(2*n+1) + i) 
                      += mat_mass_harmonic(n)(i)*Ey(n)*conj(phase(n));
		  }
	      }
	  }
	
	for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
	  {
	    // Dirichlet condition
            if (var_helmholtz.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
              {
                scheme(nb_points_z*2*n_freq) 
                  = En(nb_points_z*2*n_freq) - source_newton(nb_points_z*2*n_freq);
                
                scheme(nb_points_z*(2*n_freq+1))
                  = En(nb_points_z*(2*n_freq+1)) - source_newton(nb_points_z*(2*n_freq+1));
              }
	  }

	
        return;
      }
    
    if (!var_laplace.FirstOrderFormulation())
      {
        if (theta > 0)
          {
            Real_wp coef = square(index.c0*dt)*theta;
            MltStiffness(coef, SeldonNoTrans, En, ProdTmp);
            for (int i = 0; i < nb_points_z; i++)
              scheme(i) = source_newton(i) + ProdTmp(i);
          }
        else
          for (int i = 0; i < nb_points_z; i++)
            scheme(i) = source_newton(i);
        
        if (print_level >= 2)
          glob_chrono.Start(VirtualTimer::SCHEME);
        
        for (int i = 0; i < nb_points_z; i++)
          scheme(i) += mat_mass(i)*En(i)*(index.eps_inf + gammaNL*absSquare(En(i)))
            + mat_sigma(i)*En(i);
        
        if (print_level >= 2)
          glob_chrono.Stop(VirtualTimer::SCHEME);
      }
    else
      {
        if (print_level >= 2)
          glob_chrono.Start(VirtualTimer::SCHEME);
        
        for (int i = 0; i < nb_points_z; i++)
          scheme(i) = source_newton(i) + mat_sigma(i)*En(i)
            + mat_mass(i)*En(i)*(index.eps_inf + gammaNL*absSquare(En(i)));
        
        if (print_level >= 2)
          glob_chrono.Stop(VirtualTimer::SCHEME);
      }
    
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::PROD);
    
    // adding contribution due to linear polarization
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp dt2 = square(dt), invWk2 = 1.0/square(index.omega_polarization(k));
        Real_wp sk(0);
        if (index.sigma_polarization.GetM() > 0)
          sk = index.sigma_polarization(k);
        
        Real_wp coef = index.alpha_polarization(k)*dt2*theta_polar
          / (invWk2 + 0.5*dt*sk + theta_polar*dt2);
        
        for (int i = 0; i < nb_points_z; i++)
          scheme(i) += coef*mat_mass(i)*En(i);
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::PROD);
    
    // dirichlet condition
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      scheme(0) = En(0) - source_newton(0);
  }

  
  //! computes the jacobian matrix from the stiffness matrix
  template<class MatrixSparse>
  void MaxwellProblem1D
  ::ComputeJacobianMatrix(const Vector<Complex_wp>& En, const Vector<Real_wp>& scale,
                          MatrixSparse& mat_DF)
  {
    Real_wp coef = square(index.c0*dt)*theta;
    Mlt(coef, mat_DF);
    
    // mass part and non linear part
    for (int i = 0; i < nb_points_z; i++)
      mat_DF.AddInteraction(i, i, mat_mass(i)*(index.eps_inf + gammaNL*absSquare(En(i))));
    
    // term due to linear polarization
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp dt2 = square(dt), invWk2 = 1.0/square(index.omega_polarization(k));
        Real_wp sk(0);
        if (index.sigma_polarization.GetM() > 0)
          sk = index.sigma_polarization(k);
        
        Real_wp coef = index.alpha_polarization(k)*dt2*theta_polar
          / (invWk2 + 0.5*dt*sk + theta_polar*dt2);
        
        for (int i = 0; i < nb_points_z; i++)
          mat_DF.AddInteraction(i, i, coef*mat_mass(i));
      }
    
    // damping term
    for (int i = 0; i < nb_points_z; i++)
      mat_DF.AddInteraction(i, i, mat_sigma(i));
    
    // enforcing Dirichlet condition for first dof
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        mat_DF.ClearRow(0);
        mat_DF.Set(0, 0, 1.0);
      }
    
    /*
      uncomment the following lines if you want to check the jacobian matrix
      
      mat_DF.WriteText("mat_df.dat");
      
      int N = nb_points_z;
      Matrix<Real_wp, General, ArrayRowSparse> DF_num(N, N);
      Real_wp h = 1e-6;
      Vector<Complex_wp> Etmp(N), f_plus(N), f_minus(N);
      f_plus.Fill(0); f_minus.Fill(0);
      for (int j = 0; j < N; j++)
      {
      Copy(En, Etmp); Etmp(j) += h;
      ComputeScheme(Etmp, f_plus);
      
      Copy(En, Etmp); Etmp(j) -= h;
      ComputeScheme(Etmp, f_minus);
      
      for (int i = 0; i < N; i++)
      {
      Real_wp val = real((f_plus(i) - f_minus(i))/(2.0*h));
      if (abs(val) > 1e-12)
      DF_num.AddInteraction(i, j, val);
      }
      }
      
      DF_num.WriteText("df_num.dat"); */
    
    // factorisation of mat_DF    
  }
  
    
  //! computing jacobian matrix and factorising it
  void MaxwellProblem1D::ComputeAndFactoriseDiff(const Vector<Complex_wp>& En,
                                                    const Vector<Real_wp>& scale)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::FACTO);
    
    if (scale.GetM() > 0)
      {
        cout << "not implemented" << endl;
        abort();
      }
    
    if (harmonic_resolution)
      {
	abort();
	// a decommenter apres correction
        //if (fem_harmonic_matrix.GetM() != omega.GetM())
	//{
	//  fem_harmonic_matrix.Reallocate(omega.GetM());
	//  fem_harmonic_mat_lu.Reallocate(omega.GetM());
	//}
        
        GlobalGenericMatrix<Complex_wp> nat_mat;

        // loop over frequencies
        for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
          {
            if (imag(kwave(n_freq)) != Real_wp(0))
              {
                cout << "Case not handled" << endl;
                abort();
              }
            
	    var_helmholtz.SetOmega(omega(n_freq));	    
	    for (int ref = 0; ref < var_helmholtz.ref_rho.GetM(); ref++)
	      {
		if (omega(n_freq) == 0.0)
		  var_helmholtz.ref_rho(ref).SetConstant(Complex_wp(1, 0));
		else
		  var_helmholtz.ref_rho(ref).SetConstant(square(kwave(n_freq)/omega(n_freq)));
	      }
	    
	    var_helmholtz.kwave = kwave(n_freq);
	    var_helmholtz.pml_damping.SetMaximumVelocity(var_helmholtz.GetVelocityOfMedia(1));
            
            // computing geometrical quantities in order to compute global matrix
            var_helmholtz.ComputeMassMatrix();

	    abort();

	    // a decommenter apres correction
            //if (fem_harmonic_matrix(n_freq).GetM() == nb_points_z)
	    //fem_harmonic_matrix(n_freq).Fill(0);
            //else
	    //{
	    //  fem_harmonic_matrix(n_freq).Reallocate(nb_points_z, nb_points_z);
	    //  fem_harmonic_matrix(n_freq).Fill(0);
	    // }
            
            //var_helmholtz.AddMatrixFEM(fem_harmonic_matrix(n_freq), nat_mat);
            //var_helmholtz.AddBoundaryTerms(fem_harmonic_matrix(n_freq), nat_mat);
	    
            //fem_harmonic_mat_lu(n_freq) = fem_harmonic_matrix(n_freq);
            //fem_harmonic_mat_lu(n_freq).Factorize();
          }
        
        return;
      }
    
    if (var_laplace.FirstOrderFormulation())
      {
        diag_inv_DF.Reallocate(nb_points_z);
        for (int i = 0; i < nb_points_z; i++)
          diag_inv_DF(i) = mat_sigma(i) + mat_mass(i)*(index.eps_inf + gammaNL*absSquare(En(i)));
        
        // term due to linear polarization
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          {
            Real_wp dt2 = square(dt), invWk2 = 1.0/square(index.omega_polarization(k));
            Real_wp sk(0);
            if (index.sigma_polarization.GetM() > 0)
              sk = index.sigma_polarization(k);
            
            Real_wp coef = index.alpha_polarization(k)*dt2*theta_polar 
              / (invWk2 + 0.5*dt*sk + theta_polar*dt2);
            
            for (int i = 0; i < nb_points_z; i++)
              diag_inv_DF(i) += coef*mat_mass(i);
          }
        
        // Dirichlet condition
        if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
          diag_inv_DF(0) = 1.0;
	
        // inverting diagonal
        for (int i = 0; i < nb_points_z; i++)
          diag_inv_DF(i) = 1.0/diag_inv_DF(i);
        
        return;	
      }
    
    if (theta == 0)
      {
	// for explicit schemes the jacobian matrix is diagonal (because of mass lumping)
        diag_inv_DF.Reallocate(nb_points_z);
        for (int i = 0; i < nb_points_z; i++)
          diag_inv_DF(i) = mat_sigma(i) + mat_mass(i)*(index.eps_inf + gammaNL*absSquare(En(i)));
        
	// term due to linear polarization
	for (int k = 0; k < index.omega_polarization.GetM(); k++)
	  {
            Real_wp dt2 = square(dt), invWk2 = 1.0/square(index.omega_polarization(k));
            Real_wp sk(0);
            if (index.sigma_polarization.GetM() > 0)
              sk = index.sigma_polarization(k);
            
            Real_wp coef = index.alpha_polarization(k)*dt2*theta_polar
              / (invWk2 + 0.5*dt*sk + theta_polar*dt2);
	    
            for (int i = 0; i < nb_points_z; i++)
	      diag_inv_DF(i) += coef*mat_mass(i);
	  }
	
        // Dirichlet condition
        if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
	  diag_inv_DF(0) = 1.0;
	
        // inverting diagonal
        for (int i = 0; i < nb_points_z; i++)
          diag_inv_DF(i) = 1.0/diag_inv_DF(i);
        
        return;
      }
    
    switch (type_storage_jacobian)
      {
      case SPARSE_STORAGE : 
        mat_sparse_DF = mat_sparse_stiff;
        ComputeJacobianMatrix(En, scale, mat_sparse_DF);
        mat_sparse_lu.Factorize(mat_sparse_DF, true, true);
        break;
      case BAND_STORAGE :
        mat_band_DF = mat_band_stiff;
        ComputeJacobianMatrix(En, scale, mat_band_DF);
        mat_band_DF.Factorize();
        break;
      case TINY_BAND_STORAGE :
	var_stiff->ComputeJacobian(En, scale);
        break;
      default :
        cout << "Case not handled" << endl;
        abort();
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::FACTO);
  }
  
  
  //! solves by jacobian matrix : x_sol solves DF x_sol = rhs
  void MaxwellProblem1D::SolveDifferential(const Vector<Complex_wp>& rhs,
                                              Vector<Complex_wp>& x_sol)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::SOLVE);
    
    if (harmonic_resolution)
      {
        // loop over frequencies
        Vector<Complex_wp> Ex(nb_points_z), Ey(nb_points_z);	
	x_sol.Fill(0);
        for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
          {
            for (int i = 0; i < nb_points_z; i++)
              {
                Ex(i) = rhs(nb_points_z*2*n_freq + i);
                Ey(i) = rhs(nb_points_z*(2*n_freq+1) + i);
              }

	    abort();
	    // a decommenter apres correction
            //fem_harmonic_mat_lu(n_freq).Solve(Ex);
            //fem_harmonic_mat_lu(n_freq).Solve(Ey);
            
            for (int i = 0; i < nb_points_z; i++)
              {
                x_sol(nb_points_z*2*n_freq + i) = Ex(i);
                x_sol(nb_points_z*(2*n_freq+1) + i) = Ey(i);
              }
          }
        
        return;
      }

    if ((theta == 0) || var_laplace.FirstOrderFormulation())
      {
        for (int i = 0; i < nb_points_z; i++)
          x_sol(i) = rhs(i)*diag_inv_DF(i);
      }
    else
      {
        Copy(rhs, x_sol);
        switch (type_storage_jacobian)
          {
          case SPARSE_STORAGE : mat_sparse_lu.Solve(x_sol); break;
          case BAND_STORAGE : mat_band_DF.Solve(x_sol); break;
          case TINY_BAND_STORAGE : var_stiff->SolveTiny(x_sol); break;
          default :
            cout << "Case not handled" << endl;
            abort();
          }
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::SOLVE);
  }
  

  //! computes Px = f(tn, X) where the evolution system is equal to dU/dt = f(t, U)
  void MaxwellProblem1D
  ::EvaluateFunction(const Real_wp& tn, const VectComplex_wp& X, VectComplex_wp& Px,
                     bool invert, bool source)
  {
    if (!var_laplace.FirstOrderFormulation())
      {
        cout << "Scheme not compatible with second-order formulation" << endl;
        abort();
      }
    
    // X contains (E, P_k, H, Q_k)

    VectComplex_wp E_temp(nb_points_z);
    E_temp.Fill(0);
    
    // getting E from displacement D
    GetElectricFieldFromDisplacement(tn, X, E_temp);
    
    // enforcing Dirichlet condition on the left side
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	if (source)
	  E_temp(0) = GetPulseTime(tn);
	else
	  SetComplexZero(E_temp(0));
      }
    
    // extracting H
    int offset = nb_points_z*(1+index.omega_polarization.GetM());
    for (int i = 0; i < H_current.GetM(); i++)
      H_current(i) = X(offset + i);

    Real_wp coef = index.c0;
    MltStiffness(coef, SeldonNoTrans, H_current, ProdTmp);
    MltStiffness(coef, SeldonTrans, E_temp, H_current);

    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_NEUMANN)
      if (source)
	ProdTmp(0) -= GetPulseTime(tn);
    
    // forming Px = | -c0 \int H d phi / dz - \sigma \int E phi
    //              | omega_k Q_k
    //              | c0 \int dE/dz psi 
    //              | omega_k ( \alpha_k E - P_k)
    if (invert)
      for (int i = 0; i < nb_points_z; i++)
	Px(i) = (-ProdTmp(i) - mat_sigma(i)*E_temp(i))*invMat_mass(i);
    else
      for (int i = 0; i < nb_points_z; i++)
	Px(i) = (-ProdTmp(i) - mat_sigma(i)*E_temp(i));
      
    offset = nb_points_z*(1+index.omega_polarization.GetM());
    if (invert)
      for (int i = 0; i < H_current.GetM(); i++)
	Px(offset+i) = H_current(i)*invMat_mass(nb_points_z+i);
    else
      for (int i = 0; i < H_current.GetM(); i++)
	Px(offset+i) = H_current(i);
    
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
            Px(offsetQ+i) = wk*(alpha*E_temp(i) - X(offset+i) - sk*X(offsetQ+i));
          }
        
        offset += nb_points_z;
        offsetQ += nb_points_z;
      }

    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      Px(0) = 0;
  }
    

  //! computes ProdX = f(tn, X, Xp) where the evolution system is equal to 
  //! d^2 U/dt^2 = f(t, U, dU/dt)
  void MaxwellProblem1D
  ::EvaluateFunctionS(const Real_wp& tn, const VectComplex_wp& X, const VectComplex_wp& Xp,
		      VectComplex_wp& ProdX, bool invert, bool source)
  {
    // warning : only used for second order formulation
    if (var_laplace.FirstOrderFormulation())
      {
        cout << "Scheme not compatible with first-order formulation" << endl;
        abort();
      }
    
    // X contains (E, P_k)
    
    // getting E from displacement D
    GetElectricFieldFromDisplacement(tn, X, E_current);

    // enforcing Dirichlet condition on the left side
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      E_current(0) = GetPulseTime(tn);
    
    Real_wp coef = index.c0*index.c0;
    MltStiffness(coef, SeldonNoTrans, E_current, ProdTmp);
    
    // forming result = | - \int d E/dz d \varphi_i dz
    //                  | \omega_k^2 (alpha_k E - P_k)
    for (int i = 0; i < nb_points_z; i++)
      ProdX(i) = -ProdTmp(i)*invMat_mass(i);

    if ((var_laplace.GetRightBoundaryCondition() == BoundaryConditionEnum::LINE_ABSORBING)
        || (index.alpha_damping != Real_wp(0)))
      {
	VectComplex_wp Eprime(nb_points_z);
        GetDeriveElectricFieldFromDisplacement(tn, X, Xp, E_current, Eprime);
	
	for (int i = 0; i < nb_points_z; i++)
	  ProdX(i) -= mat_sigma(i)*invMat_mass(i)*Eprime(i);
      }
    
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::PROD);

    int offset = nb_points_z;
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp wk2 = square(index.omega_polarization(k));
        Real_wp alpha = index.alpha_polarization(k);
        if (index.sigma_polarization.GetM() > 0)
          {
            Real_wp sk = index.sigma_polarization(k);
            for (int i = 0; i < nb_points_z; i++)
              ProdX(offset + i) = wk2*(alpha*E_current(i) - X(offset+i) - sk*Xp(offset+i));
          }
        else
          for (int i = 0; i < nb_points_z; i++)
            ProdX(offset + i) = wk2*(alpha*E_current(i) - X(offset+i));
        
        offset += nb_points_z;
      }
    
    if (this->var_raman.GetAlpha() != 0)
      {
	int offset_r = nb_points_z*(index.omega_polarization.GetM()+1);
	Real_wp alpha_r, beta_r, gamma_r;
	this->var_raman.GetCoefPDE(alpha_r, beta_r, gamma_r);
	for (int i = 0; i < nb_points_z; i++)
	  {
	    Real_wp moduleE2 = absSquare(E_current(i));
	    ProdX(offset_r + i) = alpha_r*moduleE2 - beta_r*X(offset_r+i) - gamma_r*Xp(offset_r+i);
	  }
      }

    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::PROD);

    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      ProdX(0) = 0;
  }  


  void MaxwellProblem1D::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
						int nb_deriv, Vector<Real_wp>& b_src)
  {
    // implemented in complex only
    abort();
  }


  //! Adds alpha F to b_src
  void MaxwellProblem1D::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
						int nb_deriv, Vector<Complex_wp>& b_src)
  {
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      {	
	Complex_wp E0 = GetPulseTime(t);

	Real_wp coef = index.c0;
	int offset = nb_points_z*(1+index.omega_polarization.GetM());
	var_stiff->MltFirstBlock(E0, coef*alpha, offset, b_src);
	
	int offsetQ = offset + H_current.GetM();
	for (int k = 0; k < index.omega_polarization.GetM(); k++)
	  {
	    Real_wp wk = index.omega_polarization(k);
	    Real_wp alpha_ = alpha*index.alpha_polarization(k);
	    b_src(offsetQ) += wk*alpha_*E0;
	  }
      }

    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_NEUMANN)
      b_src(0) += alpha*GetPulseTime(t);
  }

  void MaxwellProblem1D::FillInitialCondition(VectComplex_wp& E0)
  {
    // condition initiale
    Real_wp z_center = 1.0; // centre physique
    Real_wp beta = 10.0;
    Real_wp k = omega_L*t0_adim;
    for (int i = 0; i < nb_points_z; i++)
      {
	Real_wp z = var_laplace.GetCoordinateDof(i)*z0_adim;
	E0(i) = exp(-beta*square(z - z_center))*exp(Iwp*k*z);	
      }

    for (int i = nb_points_z; i < var_laplace.GetNbDof(); i++)
      {
	Real_wp z = var_laplace.GetCoordinateDof(i)*z0_adim;
	Complex_wp term = exp(-beta*square(z-z_center))*exp(Iwp*k*z);
	E0(i) = -term;
      }
    
    GetElectricFieldFromDisplacement(0.0, E0, E_current);
  }

  
  //! the time schemes asks the user to factorize
  //! alpha M_h + gamma A_h
  //! gamma_implist = \gamma * \delta_t
  void MaxwellProblem1D::FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, 
							const Real_wp& beta,
							const Real_wp& gamma)
  {
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(alpha);
    nat_mat.SetCoefStiffness(gamma*gamma);
    nat_mat.SetCoefDamping(gamma);
    gamma_implicit = gamma;

    var_laplace.SetFirstOrderFormulation(false);
    var_laplace.ComputeMassMatrix();

    int m = var_laplace.mesh.GetOrder();
    int Nscal = var_laplace.mesh.GetNbDof();
    mat_band_DF.Reallocate(Nscal, Nscal, m, m);
    var_laplace.AddMatrixFEM(mat_band_DF, nat_mat);
    var_laplace.AddBoundaryTerms(mat_band_DF, nat_mat);

    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	IVect col_number(1); col_number(0) = 0;
	EraseCol(col_number, mat_band_DF);
	EraseRow(col_number, mat_band_DF);
	
	mat_band_DF.AddInteraction(0, 0, 1.0);
      }
    
    // mat_band_DF.WriteText("mat.dat") ;

    // methode plus lente a priori
    // mat_band_DF.Factorize();

    // methode rapide
    var_stiff->FactorizeFast(mat_band_DF, false);

    var_laplace.SetFirstOrderFormulation(true);
    
    /*int N = var_laplace.GetNbDof();
    mat_first_lu.Reallocate(N, N);

    VectComplex_wp Ones(N), MhOnes(N);
    for (int j = 0; j < N; j++)
      {
	Ones.Zero();
	Ones(j) = 1.0;
	
	EvaluateFunction(0.0, Ones, MhOnes, false);
	for (int i = 0; i < N; i++)
	  mat_first_lu(i, j) = real(MhOnes(i));
      }

    Mlt(-gamma, mat_first_lu);

    for (int i = 0; i < N; i++)
      mat_first_lu(i, i) += alpha/invMat_mass(i);

      GetLU(mat_first_lu, pivot_DF);*/
  }  

  void MaxwellProblem1D::SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectComplex_wp& X,
						    VectComplex_wp& Y)
  {
    /* Copy(X, Y);
    SolveLU(mat_first_lu, pivot_DF, Y);
    
    return; */
    
    int NvectU = var_laplace.mesh.GetNbDof();
    int Ndof_global = var_laplace.GetNbDof();
    int NvectV = Ndof_global - NvectU;

    // extracting H
    int offset = NvectU*(1+index.omega_polarization.GetM());
    Real_wp coef = index.c0;
    
    // first step : we compute F_U + gamma deltat R_h B_h^-1 Fv
    VectComplex_wp Fu(NvectU);
    Fu.Zero();

    // X = (Fu Fv)'
    // first we compute R_h B_h^-1 Fv
    for(int i = 0; i < NvectV; i++)
      H_current(i) = invMat_mass(NvectU+i)*X(NvectU + i);

    MltStiffness(-coef, SeldonNoTrans, H_current, ProdTmp);
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      SetComplexZero(ProdTmp(0));
    
    for(int i = 0; i < NvectU; i++)
      Fu(i) = X(i) + gamma_implicit*ProdTmp(i);

    // second step : we solve
    // (D_h + gamma_implicit*S_h - gamma_implicit² R_h B_h(-1) R_h(*)) U = Fu
    
    // methode plus lente de resolution
    // mat_band_DF.Solve(Fu);
    
    // resolution rapide
    var_stiff->SolveFast(Fu);

    // we fill the U part in Y
    for(int i = 0; i < NvectU; i++)
      Y(i) = Fu(i);

    // third step : we reconstruct V
    // V = B_h^{-1} ( F_v - gamma dt Rh' U)
    MltStiffness(-coef, SeldonTrans, Fu, H_current);
    
    // we fill the V part in Y
    for(int i = 0; i < NvectV; i++)
      Y(offset + i) = invMat_mass(NvectU+i)*(X(NvectU+i) - gamma_implicit*H_current(i));
    
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      Y(0) = 0;
  }


  void MaxwellProblem1D::SolveMassMatrix(VectComplex_wp& x)
  {
    int NvectU = var_laplace.mesh.GetNbDof();
    int Ndof_global = var_laplace.GetNbDof();
    int NvectV = Ndof_global - NvectU;

    Real_wp invEps = 1.0/index.eps_inf;
    for (int i = 0; i < nb_points_z; i++)
      x(i) *= invEps*invMat_mass(i);
    
    for (int i = 0; i < NvectV; i++)
      x(nb_points_z+i) *= invMat_mass(nb_points_z+i);
    
    if (index.omega_polarization.GetM() > 0)
      {
	cout << "Not implemented " << endl;
	abort();
      }
  }

  
  void MaxwellProblem1D
  ::ApplyMassMatrix(const Real_wp& alpha, const Real_wp& tn, const VectComplex_wp& Uh,
		    const Real_wp& beta, VectComplex_wp& Prod_Uh)
  {
    int NvectU = var_laplace.mesh.GetNbDof();
    int Ndof_global = var_laplace.GetNbDof();
    int NvectV = Ndof_global - NvectU;

    if (beta == Real_wp(0))
      Prod_Uh.Zero();
    else
      Mlt(beta, Prod_Uh);
    
    Real_wp Eps = index.eps_inf;
    for (int i = 0; i < nb_points_z; i++)
      Prod_Uh(i) += alpha*Eps*mat_mass(i)*Uh(i);
    
    for (int i = 0; i < NvectV; i++)
      Prod_Uh(nb_points_z+i) += alpha*mat_mass(nb_points_z+i)*Uh(nb_points_z+i);
    
    if (index.omega_polarization.GetM() > 0)
      {
	cout << "Not implemented " << endl;
	abort();
      }
  }
  
  
  //! returns the two-norm of a vector
  Real_wp MaxwellProblem1D::GetNorm2Vector(const Vector<Complex_wp>& x)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
    
    if (nb_proc > 1)
      {
        DistributedVector<Complex_wp> xdis(var_laplace.GetOverlappedDofNumber(), var_laplace.comm_group_mode);
        
        xdis.SetData(x.GetM(), x.GetData());
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


  void MaxwellProblem1D::FactorizeOperatorReal(const VectReal_wp& alpha, 
					       const VectReal_wp& beta,
					       const VectReal_wp& gamma)    
  {
    GlobalGenericMatrix<Real_wp> nat_mat;
    int nb_stages = gamma.GetM();

    var_laplace.SetFirstOrderFormulation(false);
    var_laplace.ComputeMassMatrix();

    vect_dt_gamma_implicit.Reallocate(nb_stages);
    
    for(int i = 0; i < nb_stages; i++)
      {
	nat_mat.SetCoefMass(alpha(i));
	nat_mat.SetCoefStiffness(gamma(i)*gamma(i));
	nat_mat.SetCoefDamping(gamma(i));
	vect_dt_gamma_implicit(i) = gamma(i);

	int m = var_laplace.mesh.GetOrder();
	int Nscal = var_laplace.mesh.GetNbDof();
	mat_band_DF.Reallocate(Nscal, Nscal, m, m);
	var_laplace.AddMatrixFEM(mat_band_DF, nat_mat);
	var_laplace.AddBoundaryTerms(mat_band_DF, nat_mat);

	if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
	  {
	    IVect col_number(1); col_number(0) = 0;
	    EraseRow(col_number, mat_band_DF);
	    EraseCol(col_number, mat_band_DF);
	    
	    mat_band_DF.AddInteraction(0, 0, 1.0);
	  }
    
	// methode plus lente a priori
	// vect_mat_band_DF(i).Factorize();

	// methode rapide
	var_stiff->FactorizeReal(mat_band_DF, false, i);
	
      }
    
    var_laplace.SetFirstOrderFormulation(true);
	
  }

  void MaxwellProblem1D::SolveOperatorReal(const Real_wp& t, const VectComplex_wp& X, 
					   VectComplex_wp& Y, int num_system)
  {
    int NvectU = var_laplace.mesh.GetNbDof();
    int Ndof_global = var_laplace.GetNbDof();
    int NvectV = Ndof_global - NvectU;

    // extracting H
    int offset = NvectU*(1+index.omega_polarization.GetM());
    Real_wp coef = index.c0;
    
    // first step : we compute F_U + gamma deltat R_h B_h^-1 Fv
    VectComplex_wp Fu(NvectU);
    Fu.Zero();

    // X = (Fu Fv)'
    // first we compute R_h B_h^-1 Fv
    for(int i = 0; i < NvectV; i++)
      H_current(i) = invMat_mass(NvectU+i)*X(NvectU + i);

    MltStiffness(-coef, SeldonNoTrans, H_current, ProdTmp);
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      SetComplexZero(ProdTmp(0));
    
    for(int i = 0; i < NvectU; i++)
      Fu(i) = X(i) + vect_dt_gamma_implicit(num_system)*ProdTmp(i);

    // second step : we solve
    // (D_h + gamma_implicit*S_h - gamma_implicit² R_h B_h(-1) R_h(*)) U = Fu
    
    // methode plus lente de resolution
    //vect_mat_band_DF(num_system).Solve(Fu);

    // resolution rapide
    var_stiff->SolveReal(num_system, Fu);
    
    // we fill the U part in Y
    for(int i = 0; i < NvectU; i++)
      Y(i) = Fu(i);

    // third step : we reconstruct V
    // V = B_h^{-1} ( F_v - gamma dt Rh' U)
    MltStiffness(-coef, SeldonTrans, Fu, H_current);
    
    // we fill the V part in Y
    for(int i = 0; i < NvectV; i++)
      Y(offset + i) = invMat_mass(NvectU+i)*(X(NvectU+i) - 
					     vect_dt_gamma_implicit(num_system)*H_current(i));
    
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      Y(0) = 0;

  }
 

  //! Factorizeoperator complex
  void MaxwellProblem1D::FactorizeOperatorComplex(const VectComplex_wp& alpha, 
						  const VectComplex_wp& beta,
						  const VectComplex_wp& gamma)    
  {
    GlobalGenericMatrix<Complex_wp> nat_mat;
    // number of complex conjugate roots
    int nb_cplxRoots = gamma.GetM();

    var_laplace.SetFirstOrderFormulation(false);
    var_laplace.ComputeMassMatrix();

    vect_dt_gamma_implicit_cplx.Reallocate(nb_cplxRoots);
    Matrix<Complex_wp, General, BandedCol> mat_band_DF_cplx;

    for(int i = 0; i < nb_cplxRoots; i++)
      {
	nat_mat.SetCoefMass(alpha(i));
	nat_mat.SetCoefStiffness(gamma(i)*gamma(i));
	nat_mat.SetCoefDamping(gamma(i));
	vect_dt_gamma_implicit_cplx(i) = gamma(i);

	int m = var_laplace.mesh.GetOrder();
	int Nscal = var_laplace.mesh.GetNbDof();
	mat_band_DF_cplx.Reallocate(Nscal, Nscal, m, m);
	var_laplace.AddMatrixFEM(mat_band_DF_cplx, nat_mat);
	var_laplace.AddBoundaryTerms(mat_band_DF_cplx, nat_mat);

	if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
	  {
	    IVect col_number(1); col_number(0) = 0;
	    EraseCol(col_number, mat_band_DF_cplx);
	    EraseRow(col_number, mat_band_DF_cplx);
	    
	    mat_band_DF_cplx.AddInteraction(0, 0, Complex_wp(1,0));
	  }
    
	var_stiff->FactorizeComplex(mat_band_DF_cplx, false, i);
      }
    
    var_laplace.SetFirstOrderFormulation(true);
	
  }


  void MaxwellProblem1D::SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, 
					      VectComplex_wp& Y, int num_system)
  {
    int NvectU = var_laplace.mesh.GetNbDof();
    int Ndof_global = var_laplace.GetNbDof();
    int NvectV = Ndof_global - NvectU;

    // extracting H
    int offset = NvectU*(1+index.omega_polarization.GetM());
    Real_wp coef = index.c0;
    
    // first step : we compute F_U + gamma deltat R_h B_h^-1 Fv
    VectComplex_wp Fu(NvectU);
    Fu.Zero();

    // X = (Fu Fv)'
    // first we compute R_h B_h^-1 Fv
    for(int i = 0; i < NvectV; i++)
      H_current(i) = invMat_mass(NvectU+i)*X(NvectU + i);

    MltStiffness(-coef, SeldonNoTrans, H_current, ProdTmp);
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      SetComplexZero(ProdTmp(0));
    
    for(int i = 0; i < NvectU; i++)
      Fu(i) = X(i) + vect_dt_gamma_implicit_cplx(num_system)*ProdTmp(i);

    // second step : we solve
    // (D_h + gamma_implicit*S_h - gamma_implicit² R_h B_h(-1) R_h(*)) U = Fu
    
    // methode plus lente de resolution
    // vect_mat_band_DF_cplx(num_system).Solve(Fu);

    // methode rapide de resolution
    var_stiff->SolveComplex(num_system, Fu);
    
    // we fill the U part in Y
    for(int i = 0; i < NvectU; i++)
      Y(i) = Fu(i);

    // third step : we reconstruct V
    // V = B_h^{-1} ( F_v - gamma dt Rh' U)
    MltStiffness(-coef, SeldonTrans, Fu, H_current);
    
    // we fill the V part in Y
    for(int i = 0; i < NvectV; i++)
      Y(offset + i) = invMat_mass(NvectU+i)*(X(NvectU+i) - 
					     vect_dt_gamma_implicit_cplx(num_system)*H_current(i));
    
    if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      Y(0) = 0;

  }

 
  //! advancing time-scheme, computes E^n+1 from E^n and E^n-1
  void MaxwellProblem1D::AdvanceScheme(int nt, const Real_wp& t,
				       All_TimeScheme<Complex_wp>& time_scheme)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::ALL);

    if (time_scheme.GetTimeSchemeType() == TimeSchemeEnum::THETA_SCHEME)
      {
        // right hand side of non-linear system is computed
        EvaluateRightHandSide(t+dt);
        
        // initializing E^n+1 = E^n
        E_next.Reallocate(nb_points_z);
        Copy(E_current, E_next);
        
        // on impose Dirichlet non-homogene sur E_next
        if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
          E_next(0) = source_newton(0);
        
        if (gammaNL == 0.0)
          {
            for (int i = 1; i < source_newton.GetM(); i++)
              source_newton(i) = -source_newton(i);
            
            // linear scheme
            SolveDifferential(source_newton, E_next);
            
            if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
              E_next(0) = source_newton(0);        
          }
        else
          {
            // print_level = 0 => no messages from Newton solver
            //newton_solver.SetPrintLevel(2);
            newton_solver.SetPrintLevel(0);
            
            // computes E^n+1 with Newton method
            newton_solver.Solve(*this, E_next);
          }
        
        // we consider the next time step
        UpdateScheme(nt, t);
      }
    else
      {
        time_scheme.Advance(t, nt, *this);
        //time_scheme.GetIterate().Write("En" + to_str(nt) + ".dat");
      }        
    
    UpdateMesh(nt, t, time_scheme);

    if (time_scheme.GetTimeSchemeType() != TimeSchemeEnum::THETA_SCHEME)
      {
        // retrieving E from displacement
        GetElectricFieldFromDisplacement(t+dt, time_scheme.GetIterate(), E_current);

	// enforcing Dirichlet condition on the left side
        if (var_laplace.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
          E_current(0) = GetPulseTime(t+dt);
      }

    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::ALL);

    if (print_level >= 4)
      if (nt%100 == 0)
        {
          //cout.setf(ios::scientific);
          cout << "At time t = " << t*t0_adim << endl;
          cout << " ||En || = " << GetNorm2Vector(E_current)*E0_adim << endl;
        }

    /* for (int i = 0; i < E_current.GetM(); i++)
      {
        Real_wp x = real(E_current(i));
        if ( (x != 0) && abs(x) < numeric_limits<double>::min() )
          {
            DISP(x);
            cout << "Denormalized value" << endl;
            // it's denormalized
          }   
          } */
 
    if (print_level >= 2)
      {
        if (nt%1000 == 0)
          {
            Real_wp normE = GetNorm2Vector(E_current);
            if (rank_proc == 0)
              {
                //cout.setf(ios::fixed);
                cout << "At time t = " << t*t0_adim << endl;
                /*cout << "Time spent in stiffness product : "
                  << glob_chrono.GetSeconds(VirtualTimer::STIFFNESS) << endl;
                  cout << "Time spent in inversion of displacement : " 
                  << glob_chrono.GetSeconds(VirtualTimer::FLUX) << endl;
                  cout << "Time spent in Solve function : "
                  << glob_chrono.GetSeconds(VirtualTimer::SOLVE) << endl;
                  cout << "Time spent to compute and factorize jacobian : "
                  << glob_chrono.GetSeconds(VirtualTimer::FACTO) << endl;
                  cout << "Time spent to translate back solution : "
                  << glob_chrono.GetSeconds(VirtualTimer::EXTRAPOL) << endl;
                  cout << "Time spent in polarization scheme : "
                  << glob_chrono.GetSeconds(VirtualTimer::PROD) << endl;
                  cout << "Time spent in other parts of the scheme : "
                  << glob_chrono.GetSeconds(VirtualTimer::SCHEME) << endl;
                  cout << "Time spent in outputs : "
                  << glob_chrono.GetSeconds(VirtualTimer::OUTPUT) << endl; */
                
                cout << " ||En || = " << normE*E0_adim << endl;
                if (normE > 1e100)
                  {
                    cout << "Scheme instable" << endl;
                    abort();
                  }
                
                cout << "Time spent in the global simulation : "
                     << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;
                //cout.unsetf(ios::fixed);
                //cout.setf(ios::scientific);
              }
          }
      }
    else
      {
        //cout.setf(ios::fixed);
        if (nt%10000 == 0)
          cout << "At time t = " << t*t0_adim << endl;

        //cout.unsetf(ios::fixed);
        //cout.setf(ios::scientific);
      }    
    
    
    //int test_input; cout << "we wait" << endl; cin >> test_input;
  }
  
  
  //! Updates time scheme (u^n+1 and u^n are written in u^n and u^n-1)
  void MaxwellProblem1D::UpdateScheme(int nt, const Real_wp& t)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::EXTRAPOL);
    
    // computes P^n+1 from E^n+1
    int offset = 0;
    for (int k = 0; k < index.omega_polarization.GetM(); k++)
      {
        Real_wp dt2 = square(dt), invWk2 = 1.0/square(index.omega_polarization(k));
        Real_wp sk(0);
        if (index.sigma_polarization.GetM() > 0)
          sk = index.sigma_polarization(k);
        
        Real_wp beta = index.alpha_polarization(k)*dt2*theta_polar
          / (invWk2 + 0.5*dt*sk + theta_polar*dt2);
	
        for (int i = 0; i < nb_points_z; i++)
	  P_next(offset + i) = beta*E_next(i) + rhs_polar(offset + i);
	
	offset += nb_points_z;
      }
    
    // for first-order formulation, H^{n+3/2} is computed
    if (var_laplace.FirstOrderFormulation())
      {
        Real_wp coef = index.c0*dt;
	MltStiffness(coef, SeldonTrans, E_next, H_next);
	int Nscal = E_next.GetM();	
	for (int i = 0; i < H_next.GetM(); i++)
	  {
	    H_next(i) = (mat_mass(Nscal+i) - mat_sigma(Nscal+i))*H_current(i) + H_next(i);
	    H_next(i) /= mat_mass(Nscal+i) + mat_sigma(Nscal+i);
	  }
	
	Copy(H_next, H_current);
      }
    
    // next time step
    Copy(P_current, P_prev);
    Copy(P_next, P_current);
    
    Copy(E_current, E_prev);
    Copy(E_next, E_current);
    Copy(mass_termEn, mass_termEn_prev);

    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::EXTRAPOL);
  }    
  
  
  //! translates the mesh if needed
  /*!
    The mesh is translated if needed in order to follow the solution
  */  
  void MaxwellProblem1D::UpdateMesh(int nt, const Real_wp& t,
				    All_TimeScheme<Complex_wp>& time_scheme)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
    int r = var_laplace.mesh.GetOrder();
#else
    int rank_proc(0), nb_proc(1);
#endif

    if (nb_iterations_check_mesh > 0)
      if ((nt%nb_iterations_check_mesh == 0) && (t > 2.0*Tmax))
	{
          // changing to Neumann boundary condition
          var_laplace.SetBoundaryCondition(BoundaryConditionEnum::LINE_NEUMANN,
					   var_laplace.GetRightBoundaryCondition());
	  
          // we compare with the amplitude of Dirichlet inhomogeneous condition
          Real_wp threshold = threshold_translate*amplitude_impulse;
          
          int nb_elt = 0, nb_dof = 0, nb_dof_L2 = 0;
	  if (velocity_translate > 0)
            {
              // elements located below z0 are skipped
              Real_wp z0 = (t-t0_begin_translate)*velocity_translate;
              // counting the number of elements and dofs to skip
              for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
                {
                  int ri = var_laplace.mesh.GetOrderElement(i);
                  bool drop_elt = false;
                  
                  if (var_laplace.mesh.Vertex(i+1) < z0)
                    drop_elt = true;
                  
                  if (!drop_elt)
                    break;
                  else
                    {
                      nb_elt++;
                      nb_dof += ri;
		      nb_dof_L2 += ri+1;
                    }
                }              
            }
          else
            { 
              // counting the number of elements and dofs to skip
              for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
                {
                  int ri = var_laplace.mesh.GetOrderElement(i);
                  bool drop_elt = true;
                  // last degree of freedom is not checked (since it belongs to the next element)
                  for (int j = 0; j < ri; j++)
                    if (abs(E_prev(var_laplace.mesh.GetNumberDof(i, j))) > threshold)
                      drop_elt = false;
                  
                  if (!drop_elt)
                    break;
                  else
                    {
                      nb_elt++;
                      nb_dof += ri;
		      nb_dof_L2 += ri+1;
                    }
                }
              
              // checking that value of E on last dof is close to 0
              if (rank_proc == nb_proc-1)
                if (abs(E_current(nb_points_z-1)) > threshold)
                  {
                    cout << "Value on last degree of freedom is too high" << endl;
                    cout << "Maybe the mesh should be increased" << endl;
                    abort();
                  }
            }
          
          int nb_proc_dropped = 0, nb_elt_dropped = 0, nb_elt_left = 0;
#ifdef SELDON_WITH_MPI
          if (nb_proc > 1)
            {
              IVect num(6);
              num(0) = nb_elt; num(1) = nb_dof; num(2) = nb_dof_L2;
              num(3) = var_laplace.mesh.GetNbElt(); num(4) = 0; num(5) = nb_elt;
              MPI_Comm& comm = var_laplace.comm_group_mode; MPI_Status status;
              
              IVect num_left(num);
              if (rank_proc > 0)
                {
                  MPI_Recv(num_left.GetData(), num.GetM(), MPI_INTEGER, rank_proc-1, 26, comm, &status);
                  
                  // if the previous processor has one element left, we cannot drop
                  // any element in the current processor
                  if (num_left(0) != num_left(3))
                    {
                      nb_elt = 0; nb_dof = 0; nb_dof_L2 = 0;
                      num(0) = 0; num(1) = 0; num(2) = 0;
                      num(4) = num_left(4); num(5) = num_left(5);
                    }
                  else
                    num(4) = num_left(4) + 1;
                }
              
              if (rank_proc < nb_proc-1)
                MPI_Send(num.GetData(), num.GetM(), MPI_INTEGER, rank_proc+1, 26, comm);
              
              // sending num to all the processors (especially num(4) and num(5))
              MPI_Bcast(num.GetData(), num.GetM(), MPI_INTEGER, nb_proc-1, comm);
              nb_proc_dropped = num(4); nb_elt_dropped = num(5);
              
              nb_elt_left = num_left(3);
              nb_elt = num_left(3)*nb_proc_dropped + nb_elt_dropped;
              nb_dof = nb_elt_dropped*r + 1; nb_dof_L2 = nb_elt_dropped*(r+1);
            }
#endif
          
	  if (nb_elt > 0)
	    TranslateMesh(nb_elt, nb_dof, nb_dof_L2,
                          nb_proc_dropped, nb_elt_dropped, nb_elt_left, time_scheme);
                    
	}

  }  
  
  
  //! translating mesh of N elements and Ndof dofs
  void MaxwellProblem1D::TranslateMesh(int N, int Ndof, int Ndof_H,
				       int nb_proc_d, int offset_elt, int nb_elt_left,
				       All_TimeScheme<Complex_wp>& time_scheme)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
#endif

    // points of the mesh are moved (regular mesh)
    int Ne = var_laplace.mesh.GetNbElt();
    nb_iterations_translation += N;
    for (int i = 0; i <= Ne; i++)
      var_laplace.mesh.Vertex(i) = OriginalVertex(i) + nb_iterations_translation*dz_translation;
    
    if (nb_iterations_translation > 1000000000)
      {
        for (int i = 0; i <= Ne; i++)
          OriginalVertex(i) = var_laplace.mesh.Vertex(i);
        
        nb_iterations_translation = 0;
      }
    
    if (time_scheme.GetTimeSchemeType() == TimeSchemeEnum::THETA_SCHEME)
      {
        // translating iterates
        TranslateIterate(E_prev, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
        TranslateIterate(E_current, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
        TranslateIterate(E_next, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
        TranslateIterate(mass_termEn_prev, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
        
        TranslateIterate(P_prev, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
        TranslateIterate(P_current, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
        TranslateIterate(P_next, Ndof, Ndof_H, nb_proc_d, offset_elt, nb_elt_left);
      }
    else
      {
        for (int k = 0; k < time_scheme.GetNumberOfIterates(); k++)
          TranslateIterate(time_scheme.GetIterate(k), Ndof, Ndof_H,
                           nb_proc_d, offset_elt, nb_elt_left);
      }
    
    // updating points on seismogramms
    IVect proc_sismo;
    grid_sismo.LocalizePoints(var_laplace.mesh);
    ForceUniquenessSismo(grid_sismo, proc_sismo);
    
    // exchanging seismogramms if needed
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      {
        MPI_Comm& comm = var_laplace.comm_group_mode;
        MPI_Status status; Vector<int64_t> xtmp;
        for (int i = 0; i < proc_sismo.GetM(); i++)
          if ((previous_proc_sismo(i) < nb_proc) && (proc_sismo(i) < nb_proc))
            {
              if (previous_proc_sismo(i) != proc_sismo(i))
                {
                  // the point is moved to another processor
                  if (rank_proc == previous_proc_sismo(i))
                    {
                      IVect m(2); m(0) = En_sismo(i).GetM(); m(1) = last_sismo_point_number(i);
                      MPI_Ssend(m.GetData(), 2, MPI_INTEGER, proc_sismo(i), 38, comm);
                      MpiSsend(comm, &time_begin_sismo(i), xtmp, 1,
                               proc_sismo(i), 39);
                      MpiSsend(comm, En_sismo(i), xtmp, En_sismo(i).GetM(),
                               proc_sismo(i), 40);
                      En_sismo(i).Clear();
                    }
                  else if (rank_proc == proc_sismo(i))
                    {
                      IVect m(2);
                      MPI_Recv(m.GetData(), 2, MPI_INTEGER, previous_proc_sismo(i), 38, comm, &status);
                      En_sismo(i).Reallocate(m(0));
                      last_sismo_point_number(i) = m(1);
                      MpiRecv(comm, &time_begin_sismo(i), xtmp, 1,
                              previous_proc_sismo(i), 39, status);
                      
                      MpiRecv(comm, En_sismo(i), xtmp, m(0),
                              previous_proc_sismo(i), 40, status);
                    }
                }
            }
      }
#endif
    
    previous_proc_sismo = proc_sismo;

    // updating extremities of the interpolation grid
    zmin_display += N*dz_translation;
    zmax_display += N*dz_translation;

    /*for (int i = 0; i < proc_sismo.GetM(); i++)
      if ((grid_sismo.GetGlobalCoordinate(i) <= zmax_display)
          && (grid_sismo.GetGlobalCoordinate(i) >= zmin_display))
        {
          if (proc_sismo(i) >= nb_proc)
            {
              sleep(rank_proc);
              DISP(grid_sismo.GetGlobalCoordinate(i));
              DISP(var_laplace.mesh.Vertex(0)); DISP(var_laplace.mesh.Vertex(Ne));
              sleep(nb_proc-rank_proc);
              
              abort();
            }
            }    */
  }
  
  
  //! avoids duplicates in sismogramms
  void MaxwellProblem1D::ForceUniquenessSismo(GridInterpolation<Dimension1>& sis, IVect& proc)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);

    if (nb_proc <= 1)
      return;
    
    proc.Reallocate(sis.GetNbPointsGrid());
    proc.Fill(nb_proc+1);
    for (int i = 0; i < sis.GetNbPointsGrid(); i++)
      if (sis.GetElementNumber(i) >= 0)
        proc(i) = rank_proc;
    
    IVect proc_loc(proc);
    MPI_Allreduce(proc_loc.GetData(), proc.GetData(), proc.GetM(),
                  MPI_INTEGER, MPI_MIN, var_laplace.comm_group_mode);
    
    for (int i = 0; i < sis.GetNbPointsGrid(); i++)
      if (sis.GetElementNumber(i) >= 0)
        {
          if (proc(i) != rank_proc)
            sis.SetElementNumber(i, -1);
        }
#endif
  }
  
  
  //! translates the field En of n dofs
  template<class T>
  void MaxwellProblem1D::TranslateIterate(Vector<T>& En, int Ndof, int Ndof_H,
					  int nb_proc_d, int offset_elt, int nb_elt_left)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_laplace.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    T last_value;
    if (var_laplace.FirstOrderFormulation())
      {
        if (nb_proc > 1)
          {
            cout << "not implemented" << endl; abort();
          }
        
	int Nscal = var_laplace.mesh.GetNbDof();
	int N = var_laplace.GetNbDof();
	int Nvec = N - Nscal;
	int size_sys = Nscal*(index.omega_polarization.GetM()*2+1) + Nvec;	
	if (En.GetM() != size_sys)
	  {
	    cout << "Case not implemented " << endl;
	    abort();
	  }
	
	int offset = 0;
	// translating unknowns E and P
	for (int k = 0; k <= index.omega_polarization.GetM(); k++)
	  {
	    for (int i = 0; i < nb_points_z-Ndof; i++)
	      En(offset + i) = En(offset + i + Ndof);
	    
            last_value = En(offset + nb_points_z - 1);
	    for (int i = nb_points_z-Ndof; i < nb_points_z; i++)
	      En(offset + i) = last_value;
	    
	    offset += nb_points_z;
	  }

	// translating H
	for (int i = 0; i < Nvec-Ndof_H; i++)
	  En(offset + i) = En(offset + i + Ndof_H);
	
        last_value = En(offset + Nvec - 1);
	for (int i = Nvec-Ndof_H; i < Nvec; i++)
	  En(offset + i) = last_value;
	
	offset += Nvec;
	
	// translating unknowns Q
	for (int k = 0; k < index.omega_polarization.GetM(); k++)
	  {
	    for (int i = 0; i < nb_points_z-Ndof; i++)
	      En(offset + i) = En(offset + i + Ndof);
	    
            last_value = En(offset + nb_points_z - 1);
	    for (int i = nb_points_z-Ndof; i < nb_points_z; i++)
	      En(offset + i) = last_value;
	    
	    offset += nb_points_z;
	  }
	
	return;
      }
    
    int nb_u = En.GetM() / nb_points_z;
    int offset = 0;
    if (nb_proc > 1)
      {
#ifdef SELDON_WITH_MPI
        MPI_Comm& comm = var_laplace.comm_group_mode;

        MPI_Request request_small, request_large; MPI_Status status;
	int r = var_laplace.mesh.GetOrder();
        int ns = offset_elt*r+1, nl = (nb_elt_left-offset_elt)*r+1; 
        // broadcasting last value to all processors
        Vector<T> last_val(nb_u);
        if (rank_proc == nb_proc-1)
          {            
            for (int k = 0; k < nb_u; k++)
              {
                last_val(k) = En(offset + nb_points_z - 1);
                offset += nb_points_z;
              }
          }
        
        Vector<int64_t> val_tmp;
        MpiBcast(comm, last_val, val_tmp, nb_u, nb_proc-1);
        
        // sending datas
        Vector<T> send_large, send_small, recv_small, recv_large;
        Vector<int64_t> send_tmp, recv_tmp;
        if (rank_proc >= nb_proc_d)
          {
            send_large.Reallocate(nl*nb_u);
            send_small.Reallocate(ns*nb_u);
            for (int i = 0; i < nb_u; i++)
              {
                for (int j = 0; j < ns; j++)
                  send_small(ns*i + j) = last_val(i);
                
                for (int j = 0; j < nl; j++)
                  send_large(nl*i + j) = last_val(i);                
              }
            
            if (rank_proc > nb_proc_d)
              {
                // sending values from 0 to offset_elt to the processor rank_proc-nb_proc_d-1
                for (int i = 0; i < nb_u; i++)
                  for (int j = 0; j < min(ns, nb_points_z); j++)
                    send_small(ns*i + j) = En(i*nb_points_z + j);
                
                request_small = MpiIsend(comm, send_small, val_tmp,
                                         ns*nb_u, rank_proc-nb_proc_d-1, 34);
              }
            
            // sending values from offset_elt until the end to the processor rank_proc-nb_proc_d
            for (int i = 0; i < nb_u; i++)
              for (int j = ns-1; j < nb_points_z; j++)
                send_large(nl*i + j-ns+1) = En(i*nb_points_z + j);
            
            if (nb_proc > 0)
              request_large = MpiIsend(comm, send_large, send_tmp, nl*nb_u,
                                       rank_proc-nb_proc_d, 33);
          }
        
        // receiving datas
        if (rank_proc < nb_proc-nb_proc_d)
          {
            if (rank_proc < nb_proc-nb_proc_d-1)
              {
                recv_small.Reallocate(ns*nb_u);
                MpiRecv(comm, recv_small, recv_tmp, ns*nb_u,
                        rank_proc+nb_proc_d+1, 34, status);
              }
            
            if (nb_proc_d > 0)
              {
                recv_large.Reallocate(nl*nb_u);
                MpiRecv(comm, recv_large, recv_tmp, nl*nb_u,
                        rank_proc+nb_proc_d, 33, status);
              }
            else
              recv_large = send_large;
          }
        
        // waiting end of sends
        MPI_Wait(&request_small, &status);
        MPI_Wait(&request_large, &status);
        
        // reconstructing En
        if (rank_proc < nb_proc-nb_proc_d)
          {
            for (int i = 0; i < nb_u; i++)
              for (int j = 0; j < min(nl, nb_points_z); j++)
                En(i*nb_points_z + j) = recv_large(i*nl + j);
          }
        else
          {
            for (int i = 0; i < nb_u; i++)
              for (int j = 0; j < min(nl, nb_points_z); j++)
                En(i*nb_points_z + j) = last_val(i);
          }
        
        if (rank_proc < nb_proc-nb_proc_d-1)
          {
            for (int i = 0; i < nb_u; i++)
              for (int j = nl-1; j < nb_points_z; j++)
                En(i*nb_points_z + j) = recv_small(i*ns + j-nl+1);
          }
        else
          {
            for (int i = 0; i < nb_u; i++)
              for (int j = nl-1; j < nb_points_z; j++)
                En(i*nb_points_z + j) = last_val(i);
          }
#endif
        return;
      }

    // sequential case
    offset = 0;
    for (int k = 0; k < nb_u; k++)
      {
	for (int i = 0; i < nb_points_z-Ndof; i++)
	  En(offset + i) = En(offset + i + Ndof);
	
        last_value = En(offset + nb_points_z - 1);
	for (int i = nb_points_z-Ndof; i < nb_points_z; i++)
	  En(offset + i) = last_value;
	
	offset += nb_points_z;
      }
  }
  
  
  //! calcul du second membre en harmonique
  void MaxwellProblem1D::ComputeRightHandSide(VectComplex_wp& source_rhs)
  {
    source_rhs.Reallocate(2*nb_points_z*omega.GetM());
    source_rhs.Fill(0);
    for (int i = 0; i < omega.GetM(); i++)
      {
        Complex_wp pulse = GetPulseHarmonic(omega(i));
	// Condition Ex(0) = ...
        source_rhs(2*nb_points_z*i) = pulse*polarization_init(0)*exp(-Iwp*phase_init(0));
        // Condition Ey(0) = ...
        source_rhs(2*nb_points_z*i+nb_points_z)
          = pulse*polarization_init(1)*exp(-Iwp*phase_init(1));
      }
  }
  
  
  //! writing seismogramm i
  void MaxwellProblem1D::WriteSismo(int i)
  {
    string name_file = GetBaseString(file_output_sismo) + "_P" + to_str(i) + ".dat";
    
    int N = last_sismo_point_number(i)+1;
    Real_wp t0 = -Tmax + time_begin_sismo(i);
    Real_wp t1 = -Tmax + time_begin_sismo(i) + (N-1)*dt_sismo;
    VectComplex_wp u(N);
    for (int k = 0; k < N; k++)
      u(k) = E0_adim*En_sismo(i)(k);
    
    // En is written by using the "loadND" format
    GridInterpolationFull<Dimension1> grid;
    int type = GridInterpolationFull<Dimension1>::LINE;
    grid.Init(type, t0*t0_adim, t1*t0_adim, N);
    
    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
    bool ascii = false;
    WriteMatlab(u, grid, DOSSIER+name_file, double_prec, ascii);
    
    // computing Fourier transform of the signal
    if (omega_fourier.GetM() <= 0)
      return;
    
    N = En_sismo(i).GetM();
    u.Reallocate(N);
    for (int k = 0; k < N; k++)
      u(k) = E0_adim*t0_adim*En_sismo(i)(k);
    
    if (fft_time.GetNbPoints() != N)
      fft_time.Init(N);

    // step for pulsation
    int facteur = N/nb_points_sismo;
    if (N%nb_points_sismo != 0)
      {
        DISP(N); DISP(nb_points_sismo); abort();
      }
    
    Real_wp domega = 2.0*pi_wp/(dt*nb_points_sismo);
    //DISP(N); DISP(dt); DISP(domega); DISP(omega_L);
    
    //DISP(omega_fourier(0)); DISP(omega_fourier(omega_fourier.GetM()-1));
    int n0 = toInteger(floor(omega_fourier(0)/domega));
    int n1 = toInteger(ceil(omega_fourier(omega_fourier.GetM()-1)/domega));
    if (n0 < 0)
      n0 = 0;
    
    Real_wp omega0 = domega*n0/t0_adim;
    Real_wp omega1 = domega*n1/t0_adim;
    
    if ((n1 > N/2) || (n0 > n1))
      return;
    
    int istep = 1;
    
    //DISP(n0); DISP(n1); DISP(omega0); DISP(omega1); DISP(istep);
    n0 *= facteur;
    n1 *= facteur;
    istep *= facteur;
    //DISP(n0); DISP(n1); DISP(istep);
    
    VectReal_wp ux(N), uy(N);
    for (int k = 0; k < N; k++)
      {
	ux(k) = real(u(k));
	uy(k) = imag(u(k));
      }
    
    VectComplex_wp ux_hat(N/2+1), uy_hat(N/2+1);
    fft_time.ApplyForward(ux, ux_hat);
    fft_time.ApplyForward(uy, uy_hat);
    
    Real_wp coef_fft_u = 2.0*dt/sqrt(2.0*pi_wp);
    VectComplex_wp usub_x((n1-n0)/istep+1), usub_y((n1-n0)/istep+1);
    int num = 0;
    for (int k = n0; k <= n1; k += istep)
      {
	usub_x(num) = coef_fft_u*conj(ux_hat(k));
	usub_y(num) = coef_fft_u*conj(uy_hat(k));
	num++;
      }
    
    string name_file_x = GetBaseString(file_output_sismo) + "_P" + to_str(i) + "_harmonicX.dat";
    string name_file_y = GetBaseString(file_output_sismo) + "_P" + to_str(i) + "_harmonicY.dat";
    
    grid.Init(type, omega0, omega1, usub_x.GetM());
    WriteMatlab(usub_x, grid, DOSSIER+name_file_x, double_prec, ascii);
    WriteMatlab(usub_y, grid, DOSSIER+name_file_y, double_prec, ascii);
  }
  
  
  //! snapshot is computed and written in the output file
  void MaxwellProblem1D::WriteSnapshot(int nt, const Real_wp& t, Vector<Complex_wp>& En)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::OUTPUT);
    
    // snapshot on interpolation grids
    if (var_laplace.GetOutputParameters().SnapshotToStore(t))
      {        
	// setting the name of the output file 
        int num = var_laplace.GetOutputParameters().GetSnapshotNumber();
        Real_wp norm_En = GetNorm2Vector(En);
        string name_file = GetBaseString(var_laplace.GetOutputParameters().GetFileName(1))
          + NumberToString(num) + ".dat";
        
        if (rank_proc == 0)
          {
            //cout.setf(ios::scientific);
            cout << "Writing at t = " << t*t0_adim << endl;
            cout << "|| En || = " << norm_En*E0_adim << endl;
            DISP(name_file);
          }
        
	// writing the output vector
        WriteOutputFile(var_laplace, En, name_file);
        var_laplace.GetOutputParameters().IncrementSnapshot();
      }
     
    // seismogramms
    VectReal_wp val_phi; Complex_wp val_En;
    for (int i = 0; i < En_sismo.GetM(); i++)
      {
	if (grid_sismo.GetElementNumber(i) >= 0)	
	  {
	    // evaluating En at the given point
	    val_En = var_laplace.
	      GetInterpolate(En, grid_sismo.GetElementNumber(i),
			     grid_sismo.GetLocalCoordinate(i), val_phi);
	    
	    bool write_sismo = true;
	    if (En_sismo(i).GetM() <= 0)
	      {
		if (abs(val_En) < amplitude_impulse*threshold_sismo)
		  write_sismo = false;
		else
		  {
		    time_begin_sismo(i) = t;
		    En_sismo(i).Reallocate(nb_points_sismo);
		    En_sismo(i).Fill(0);
		  }
	      }
	    
	    if (write_sismo)
	      {
		Real_wp coef = (t-time_begin_sismo(i))/dt_sismo;
		int n = toInteger(round(coef));
		if (abs(coef-n) < 1e-5)
		  {
		    if (n >= En_sismo(i).GetM())
		      {
			int nb_old = En_sismo(i).GetM();
			En_sismo(i).Resize(2*n);
			for (int k = nb_old; k < 2*n; k++)
			  En_sismo(i)(k) = 0.0;
		      }
		    
		    last_sismo_point_number(i) = n;
		    En_sismo(i)(n) = val_En;
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
  
  
  //! interpolating En and writing it on a file
  void MaxwellProblem1D::WriteOutputFile(const VarProblem_1D& var, Vector<Complex_wp>& En,
					 const string& name_file)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_laplace.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    // the interpolation of En on regular points is computed
    Vector<Complex_wp> En_interp;
    var.ComputeInterpolationU(En, var.GetSectionGrid(), En_interp);
    
    if (rank_proc == 0)
      {
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          Mlt(E0_adim, En_interp);
        
        // interpolation of En is written by using the "loadND" format
        GridInterpolationFull<Dimension1> grid;
        int type = GridInterpolationFull<Dimension1>::LINE;
        grid.Init(type, zmin_display*z0_adim, zmax_display*z0_adim,
                  var.GetSectionGrid().GetGlobalNumberPoints());
        
        int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
        bool ascii = false;
        WriteMatlab(En_interp, grid, DOSSIER+name_file, double_prec, ascii);
      }
  }

  
  //! solution of the harmonic problem
  void MaxwellProblem1D::SolveNewton(VectComplex_wp& source, VectComplex_wp& x_sol)
  {
    x_sol.Reallocate(source.GetM());
    x_sol.Fill(0);

    if (initialisation_with_schrodinger)
      {
	VectComplex_wp E0;
	KerrProblem var;
	var.Init(*this);
	
	int Nf = omega.GetM();
	if (type_polarization == ELLIPTIC)
	  E0.Reallocate(2*Nf);
	else
	  E0.Reallocate(Nf);
	
	for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
	  {
	    x_sol(nb_points_z*2*n_freq) = source(nb_points_z*2*n_freq);
	    x_sol(nb_points_z*(2*n_freq+1)) = source(nb_points_z*(2*n_freq+1));
	    E0(n_freq) = x_sol(nb_points_z*2*n_freq);
	    if (type_polarization == ELLIPTIC)
	      E0(n_freq + Nf) = x_sol(nb_points_z*(2*n_freq+1));
	  }

	// computing interval between dofs
	int r = var_helmholtz.mesh.GetOrder();
	const ElementReference<Dimension1, 1>& Fb = var_helmholtz.GetReferenceElementH1(0);
	VectReal_wp TimeStep(r);
	for (int i = 0; i < r; i++)
	  TimeStep(i) = Fb.Points(i+1) - Fb.Points(i);
	
	var.time_scheme.SetInitialCondition(0.0, var.dz, E0, var);
	
	Complex_wp coef = polarization_init(1)*exp(Iwp*(phase_init(0) - phase_init(1)));
	for (int i = 0; i < nb_points_z-1; i++)
	  {
	    int n = i/r;
	    var.dz = var_helmholtz.mesh.Vertex(n+1) - var_helmholtz.mesh.Vertex(n);
	    var.time_scheme.ChangeTimeStep(TimeStep(i%r)*var.dz);
	    var.time_scheme.Advance(var_helmholtz.GetCoordinateDof(i), i, var);
	    
	    const Vector<Complex_wp>& X = var.time_scheme.GetIterate();
	    if (type_polarization == ELLIPTIC)
	      for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
		{
		  x_sol(nb_points_z*(2*n_freq) + i+1) = X(n_freq);
		  x_sol(nb_points_z*(2*n_freq+1) + i+1) = X(n_freq + Nf);
		}
	    else
	      for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
		{
		  x_sol(nb_points_z*(2*n_freq) + i+1) = X(n_freq);
		  x_sol(nb_points_z*(2*n_freq+1) + i+1) = X(n_freq)*coef;
		}
	    
	    //X.Write("xtest.dat");
	    //int test_input; cout << "waiting..." << endl; cin >> test_input;
	  }
	
	//return;
      }
    
    VectReal_wp scale;
    ComputeAndFactoriseDiff(x_sol, scale);

    if (gammaNL == 0.0)
      {
        SolveDifferential(source, x_sol);
      }
    else
      {
        source_newton = source;
        
        newton_solver.Init(*this, x_sol);
        newton_solver.SetPrintLevel(2);
	newton_solver.EnableDecreasingResidue(false);
        
        newton_solver.Solve(*this, x_sol);
      }
  }
  
  
  //! writing seismogramm i
  void MaxwellProblem1D::WriteHarmonicSismo(VectComplex_wp& Ex, VectComplex_wp& Ey, int i)
  {
    string name_file = GetBaseString(file_output_sismo) + "_P" + to_str(i) + ".dat";

    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
    bool ascii = false;
    int type = GridInterpolationFull<Dimension1>::LINE;
    GridInterpolationFull<Dimension1> grid;
    
    Real_wp z = grid_sismo.GetGlobalCoordinate(i);
    Real_wp omega0(0), omega1(0);
    Complex_wp vg;
    vg = index.GetGroupVelocity(omega_L);
    
    if (type_spectrum == CENTERED)
      {
	omega0 = omega(0);
	omega1 = omega(omega.GetM()-1);
	
	// returning back in time
	Complex_wp phase(1.0);
	for (int j = 0; j < Ex.GetM(); j++)
	  {
	    if (var_helmholtz.EnvelopeToCompute())
	      phase = exp(Iwp*(kwave(j)-omega(j)/vg)*z);
	    else
	      phase = exp(-Iwp*omega(j)/vg*z);
	    
	    Ex(j) = conj(Ex(j)*phase);
	    Ey(j) = conj(Ey(j)*phase);
	  }
	
	VectReal_wp Ex_time(nb_points_time), Ey_time(nb_points_time);
	Ex_time.Fill(0); Ey_time.Fill(0);
	
	fft_time.ApplyInverse(Ex, Ex_time);
	fft_time.ApplyInverse(Ey, Ey_time);
	
	Mlt(E0_adim/coef_fft, Ex_time);
	Mlt(E0_adim/coef_fft, Ey_time);
    
	VectComplex_wp u(nb_points_time);
	for (int j = 0; j < nb_points_time/2; j++)
	  {
	    u(j) = Complex_wp(Ex_time(nb_points_time/2+j), Ey_time(nb_points_time/2+j));
	    u(nb_points_time/2+j) = Complex_wp(Ex_time(j), Ey_time(j));
	  }
	
	Real_wp t0 = -Tmax + z/real(vg), t1 = z/real(vg) + time_interval(nb_points_time/2-1);
	grid.Init(type, t0*t0_adim, t1*t0_adim, nb_points_time);
	
	WriteMatlab(u, grid, DOSSIER+name_file, double_prec, ascii);
      }
    else if (type_spectrum == ODD_FREQUENCIES)
      {
        // writing only harmonic seismogramms
        VectComplex_wp Ex_chap(nb_points_time), Ey_chap(nb_points_time);
        Real_wp coef = E0_adim*t0_adim;
        
        // multiplication by phase
        Complex_wp phase;
        for (int j = 0; j < Ex.GetM(); j++)
          {
	    if (var_helmholtz.EnvelopeToCompute())
	      phase = exp(Iwp*(kwave(j)-omega(j)/vg)*z);
	    else
	      phase = exp(-Iwp*omega(j)/vg*z);  
	    
            phase *= coef;
	    Ex(j) *= phase; Ey(j) *= phase;
          }
        
        // writing each frequency in a different file
        int offset = 0;
        for (int n = 0; n < number_odd_frequencies; n++)
          {
            for (int k = 0; k < nb_points_time/2; k++)
              {
                Ex_chap(k) = Ex(nb_points_time/2+offset+k);
                Ex_chap(nb_points_time/2+k) = Ex(offset+k);
                Ey_chap(k) = Ey(nb_points_time/2+offset+k);
                Ey_chap(nb_points_time/2+k) = Ey(offset+k);
              }
            
            omega0 = omega(offset+nb_points_time/2);
            omega1 = omega(offset+nb_points_time/2-1);     
            
            string name_file_x = GetBaseString(file_output_sismo) + "_P"
              + to_str(i) + "_harmonicX_Freq" + to_str(n) + ".dat";
            
            string name_file_y = GetBaseString(file_output_sismo)
              + "_P" + to_str(i) + "_harmonicY_Freq" + to_str(n) + ".dat";
            
            grid.Init(type, omega0/t0_adim, omega1/t0_adim, Ex_chap.GetM());
            WriteMatlab(Ex_chap, grid, DOSSIER+name_file_x, double_prec, ascii);
            WriteMatlab(Ey_chap, grid, DOSSIER+name_file_y, double_prec, ascii);
            offset += nb_points_time;
          }
        
        return;
      }
    else
      {
	omega0 = omega(nb_points_time/2);
	omega1 = omega(nb_points_time/2-1);

	// returning back in time
	Complex_wp phase(1.0);
	VectComplex_wp Ex_time(nb_points_time), Ey_time(nb_points_time);
	for (int j = 0; j < Ex.GetM(); j++)
	  {
	    if (var_helmholtz.EnvelopeToCompute())
	      phase = exp(Iwp*(kwave(j)-omega(j)/vg)*z);
	    else
	      phase = exp(-Iwp*omega(j)/vg*z);  
	    
	    Ex(j) *= phase; Ey(j) *= phase;
	    Ex_time(j) = Ex(j);
	    Ey_time(j) = Ey(j);
	  }
	
	fft_envelope.ApplyForward(Ex_time);
	fft_envelope.ApplyForward(Ey_time);
	
	Mlt(E0_adim/coef_fft, Ex_time);
	Mlt(E0_adim/coef_fft, Ey_time);
	
	for (int j = 0; j < nb_points_time; j++)
	  {
	    phase = exp(-Iwp*omega_L*time_interval(j));
	    Ex_time(j) *= phase;
	    Ey_time(j) *= phase;
	  }
	
	Complex_wp vloc;
	for (int j = 0; j < nb_points_time/2; j++)
	  {
	    vloc = Ex_time(nb_points_time/2+j);
	    Ex_time(nb_points_time/2+j) = Ex_time(j);
	    Ex_time(j) = vloc;
	    
	    vloc = Ey_time(nb_points_time/2+j);
	    Ey_time(nb_points_time/2+j) = Ey_time(j);
	    Ey_time(j) = vloc;
	  }
	
	Real_wp t0 = -Tmax + z/real(vg), t1 = z/real(vg) + time_interval(nb_points_time/2-1);
	grid.Init(type, t0*t0_adim, t1*t0_adim, nb_points_time);
	
	string name_file_tx = GetBaseString(file_output_sismo) + "_P" + to_str(i) + "_timeX.dat";
	string name_file_ty = GetBaseString(file_output_sismo) + "_P" + to_str(i) + "_timeY.dat";

	WriteMatlab(Ex_time, grid, DOSSIER+name_file_tx, double_prec, ascii);	
	WriteMatlab(Ey_time, grid, DOSSIER+name_file_ty, double_prec, ascii);	
      }
    
    // writing harmonic field as well
    Real_wp coef = E0_adim*t0_adim;
    if (type_spectrum == CENTERED)
      for (int j = 0; j < Ex.GetM(); j++)
	{
	  Ex(j) = coef*conj(Ex(j));
	  Ey(j) = coef*conj(Ey(j));
	}
    else
      {
	Complex_wp Eloc;
	for (int j = 0; j < nb_points_time/2; j++)
	  {
	    Eloc = Ex(j);
	    Ex(j) = coef*Ex(nb_points_time/2 + j);
	    Ex(nb_points_time/2 + j) = coef*Eloc;
	    
	    Eloc = Ey(j);
	    Ey(j) = coef*Ey(nb_points_time/2 + j);
	    Ey(nb_points_time/2 + j) = coef*Eloc;
	  }
      }
    
    string name_file_x = GetBaseString(file_output_sismo) + "_P" + to_str(i) + "_harmonicX.dat";
    string name_file_y = GetBaseString(file_output_sismo) + "_P" + to_str(i) + "_harmonicY.dat";
    
    grid.Init(type, omega0/t0_adim, omega1/t0_adim, Ex.GetM());
    WriteMatlab(Ex, grid, DOSSIER+name_file_x, double_prec, ascii);
    WriteMatlab(Ey, grid, DOSSIER+name_file_y, double_prec, ascii);
  }
  
  
  //! writing the harmonic solution
  void MaxwellProblem1D::WriteDatas(VectComplex_wp& sol)
  {
    // snapshot on interpolation grids
    VectComplex_wp Ex(nb_points_z), Ey(nb_points_z);
    for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
      {
	for (int i = 0; i < nb_points_z; i++)
	  {
	    Ex(i) = sol(nb_points_z*(2*n_freq) + i);
	    Ey(i) = sol(nb_points_z*(2*n_freq+1) + i);
	  }
	
	// setting the name of the output file 
        string name_file_x = GetBaseString(var_helmholtz.GetOutputParameters().GetFileName(1))
          + NumberToString(n_freq) + "_harmonicX.dat";
        
        string name_file_y = GetBaseString(var_helmholtz.GetOutputParameters().GetFileName(1))
          + NumberToString(n_freq) + "_harmonicY.dat";
        
	// writing the output vector
        WriteOutputFile(var_helmholtz, Ex, name_file_x);
        WriteOutputFile(var_helmholtz, Ey, name_file_y);
      }
    
    // snapshot on seismogramms
    VectComplex_wp En_x(omega.GetM()), En_y(omega.GetM());
    En_x.Fill(0); En_y.Fill(0);
    Complex_wp val_Ex, val_Ey; VectReal_wp val_phi;
    for (int i = 0; i < grid_sismo.GetNbPointsGrid(); i++)
      {
	if (grid_sismo.GetElementNumber(i) >= 0)	
	  {
	    for (int n_freq = 0; n_freq < omega.GetM(); n_freq++)
	      {
		Ex.SetData(nb_points_z, &sol(nb_points_z*(2*n_freq)));
		Ey.SetData(nb_points_z, &sol(nb_points_z*(2*n_freq+1)));
		
		// evaluating En at the given point
		val_Ex = var_helmholtz.
		  GetInterpolate(Ex, grid_sismo.GetElementNumber(i),
				 grid_sismo.GetLocalCoordinate(i), val_phi);

		val_Ey = var_helmholtz.
		  GetInterpolate(Ey, grid_sismo.GetElementNumber(i),
				 grid_sismo.GetLocalCoordinate(i), val_phi);
		
		En_x(n_freq) = val_Ex;
		En_y(n_freq) = val_Ey;
		
		Ex.Nullify(); Ey.Nullify();
	      }
	    
	    WriteHarmonicSismo(En_x, En_y, i);
	  }
      }

  }
  
};

#define MONTJOIE_FILE_NON_LINEAR_MAXWELL_1D_CXX
#endif
