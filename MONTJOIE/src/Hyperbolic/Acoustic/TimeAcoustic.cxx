#ifndef MONTJOIE_FILE_TIME_ACOUSTIC_CXX

namespace Montjoie
{
  
  
  /*********************
   * TimeAcoustic_Base *
   *********************/
  
  
  //! returns true if the media whose reference is ref contains damping
  template<class Dimension>
  bool TimeAcoustic_Base<Dimension>::IsDampedMedia(int ref)
  {
    return (!this->var_helm.ref_sigma(ref).IsZero());
  }
  
  
  //! initialization before time iterations
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>::InitTimeIterations()
  {        
#ifdef MONTJOIE_WITH_TIME_REVERSAL
    if (time_reversal.IsPresent())
      {
	time_reversal.InitComputation();
	if (time_reversal.GetSimulationType() == time_reversal.INVERSE)
	  var_time.SetFinalTimeSource(time_reversal.GetFinalTime() - time_reversal.GetInitialTime());
      }
#endif

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    
    // modification of source for HDG formulation
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	int level = -1;
	if (var_time.SourceDoesNotDependOnTime())
	  {
	    int nb_dof_L = var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
	    int Nvol = var_problem.GetNbMainUnknownDof();
	    int nb_dof_u = Nvol;
	    int nb_dof_v = Dimension::dim_N*Nvol;
	    
	    VectReal_wp Lambda(nb_dof_L);
	    VectReal_wp Uh(nb_dof_u), Vh(nb_dof_v); Uh.Zero(); Vh.Zero();
	    this->GetExtrapolationLambda(Uh, Vh, level, Lambda);
	    Lambda.Zero();
	    for (int i = 0; i < var_time.lambda_vector_source.GetM(); i++)
	      Lambda(var_time.lambda_vector_source.Index(i)) = var_time.lambda_vector_source.Value(i);

	    var_time.Glob_mat_Dh->SolveOperatorCh(Lambda);

	    VectReal_wp Y(nb_dof_u + nb_dof_v); Y.Zero();
	    for (int i = 0; i < var_time.sparse_vector_source.GetM(); i++)
	      Y(var_time.sparse_vector_source.Index(i)) = var_time.sparse_vector_source.Value(i);
	    
	    VectReal_wp Prod_Uh, Prod_Vh;
	    Prod_Uh.SetData(nb_dof_u, Y.GetData());
	    Prod_Vh.SetData(nb_dof_v, &Y(nb_dof_u));
	    if (var_time.Glob_mat_Ch_Lambda.GetM() > 0)
	      MltAdd(Real_wp(-1), var_time.Glob_mat_Ch_Lambda, Lambda, Real_wp(1), Y);
	    else
	      this->MltAddStiffness(Real_wp(1), level, level, Uh, Vh, Lambda, Real_wp(1), Prod_Uh, Prod_Vh);

	    Prod_Uh.Nullify(); Prod_Vh.Nullify();
	    
	    Copy(Y, var_time.sparse_vector_source);
	  }
      }
    else
      {
	int Nvol = mesh_num.GetNbDof();
	int Npml = mesh_num.GetNbDofPML();
	int N = var_problem.GetNbDof();
	int d = Dimension::dim_N;
	int nb_dof_u = Nvol + (d-1)*Npml;
	int nb_dof_v = N - nb_dof_u;

	// to allocate arrays
	if ( (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
             && !var_problem.FirstOrderFormulationDG())
          {
            // SIPG, nothing to do
          }
        else
          {
            VectReal_wp Uh(nb_dof_u), Vh(nb_dof_v);
            Uh.Zero(); Vh.Zero();
            this->GetExtrapolationAcoustic(Uh, Vh, -1, true, true);
          }
      }
    
    // checking that scheme is potentially stable
    switch (var_time.var_time_scheme.GetTimeSchemeType())
      {
      case TimeSchemeEnum::OPTIMAL_LEAP_FROG_PML:
        {
          // checking that there is no damping
          bool damping = false;
          for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
            {
              int ref = var_problem.mesh.BoundaryRef(i).GetReference();
              if (var_problem.mesh.GetBoundaryCondition(ref)
                  == BoundaryConditionEnum::LINE_ABSORBING)
                damping = true;
            }
          
          if (var_boundary.GetNbEltPML() > 0)
            damping = true;
          
          Vector<bool> RefUsed(var_helm.ref_sigma.GetM());
          RefUsed.Fill(false);
          for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
            {
              int ref = var_problem.mesh.Element(i).GetReference();
              RefUsed(ref) = true;
            }
          
          for (int ref = 0; ref < var_helm.ref_sigma.GetM(); ref++)
            if (RefUsed(ref))
              if (!var_helm.ref_sigma(ref).IsZero())
                damping = true;
          
          if (damping)
            {
              cout << "Scheme instable if there are damping terms in the equation" << endl;
              abort();
            }
        }
        break;
      default :
        {
          /*Real_wp cfl = var_time.var_time_scheme.GetCflScheme();
          if (cfl == Real_wp(0))
            {
              cout << "Time scheme unstable with wave equation" << endl;
              abort();
	      }*/
        }
      }
  }
  
    
  //! additional parameters for acoustic equation in data file
  /*!
    \param[in] description_field keyword
    \param[in] parameters list of values
    \param[in] nb_param number of parameters
   */
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    
#ifdef MONTJOIE_WITH_TIME_REVERSAL
    time_reversal.SetInputData(description_field, parameters);  
#endif
    
    if (description_field == "DiagrammeDirectivite")
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of TimeAcoustic_Base" << endl;
	    cout << "DiagrammeDirectivite needs more parameters, for instance :" << endl;
	    cout << "DiagrammeDirectivite = YES t0 file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        // on effectue un calcul de la norme L2 de la pression pour les points
        // entres dans SismoOutsidePoints
	
	VarComputationRCS_Laplace<Dimension>& output_rcs_param =
	  static_cast<VarComputationRCS_Laplace<Dimension>& >(var_output.GetParameterOutputRCS());
	
	bool compute_L2_norm = false;
	Real_wp t_begin_calculL2(0);
	string file_name_normL2;
        if (parameters(0) == "YES")
          {
            compute_L2_norm = true;
            t_begin_calculL2 = to_num<Real_wp>(parameters(1));
            file_name_normL2 = parameters(2);
          }
	
	output_rcs_param.SetL2NormOutput(compute_L2_norm, t_begin_calculL2, file_name_normL2);
      }
  }
  

  //! fills how memory is used in the class
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::GetMemoryUsed(map<string, size_t>& var) const
  {
    var["AcousticVector"] = Seldon::GetMemorySize(extrapolV) + Seldon::GetMemorySize(extrapolU)
      + Seldon::GetMemorySize(Vneighbor) + Seldon::GetMemorySize(Uneighbor)
      + Seldon::GetMemorySize(evalUn_quad) + Seldon::GetMemorySize(evalVn_quad)
      + evalLambda_n.GetMemorySize();
  }
  
    
  //! compute Y = G(tn,X) if we put wave equation with pml to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
                if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
			       const VectReal_wp& X, VectReal_wp& Y, bool invert, bool source)
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int nb_dof_u = 0, nb_dof_v = 0, N = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	int Nvol = var_problem.GetNbMainUnknownDof();
	int d = Dimension::dim_N;
	nb_dof_u = Nvol;
	nb_dof_v = Nvol*d;
	N = nb_dof_u + nb_dof_v;
      }
    else
      {
	int Nvol = mesh_num.GetNbDof();
	int Npml = mesh_num.GetNbDofPML();
	N = var_problem.GetNbDof();
	int d = Dimension::dim_N;
	nb_dof_u = Nvol + (d-1)*Npml;
	nb_dof_v = N - nb_dof_u;
      }
    
    Real_wp one(1), zero(0);
    
    VectReal_wp Uh, Vh;
    VectReal_wp Prod_Uh, Prod_Vh;
    
    if ((N > X.GetM()) || (N > Y.GetM()))
      {
        cout << "X or Y is not large enough" << endl;
        DISP(N);
        DISP(X.GetM()); DISP(Y.GetM());
        abort();
      }

    // we retrieve U and V
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    Uh.SetData(nb_dof_u, x_ptr); Prod_Uh.SetData(nb_dof_u, y_ptr);
    x_ptr += nb_dof_u; y_ptr += nb_dof_u;
    
    Vh.SetData(nb_dof_v, x_ptr); Prod_Vh.SetData(nb_dof_v, y_ptr);
    int level = MatrixVectorProductLevel::ALL_LEVELS;


    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {	
	// HDG case
	int nb_dof_L = mesh_num.GetNbDof();
	VectReal_wp& Lambda = evalLambda_n;
	Lambda.Reallocate(nb_dof_L);	
	Y.Zero();

	if (var_time.Glob_mat_RhS.GetM() > 0)
	  MltAdd(-Real_wp(1), var_time.Glob_mat_RhS, X, Real_wp(0), Lambda);
	else
	  this->GetExtrapolationLambda(Uh, Vh, level, Lambda);
        
	if (source)
	  this->AddTimeSourceHDG(Real_wp(1), tn, nb_deriv, Y, Lambda);

	var_time.Glob_mat_Dh->SolveOperatorCh(Lambda);

	// testing MltAddStiffnessScalar/Vectorial
	if (var_time.Glob_mat_Ch_Lambda.GetM() > 0)
	  {
	    MltAdd(Real_wp(-1), var_time.Glob_mat_Ch_Lambda, Lambda, one, Y);
	    MltAdd(Real_wp(-1), var_time.Glob_mat_RhV, X, one, Y);
	  }
	else
	  this->MltAddStiffness(one, level, level, Uh, Vh, Lambda, one, Prod_Uh, Prod_Vh);

	var_time.ApplyOperatorSh(-one, tn, Uh, one, Prod_Uh);
        
	if (invert)
	  var_time.SolveMassMatrix(Y);

        //int test_input; cout << "Waiting" << endl; cin >> test_input;
      }	
    else
      {
	// general case (LDG or continuous)
	this->GetExtrapolationAcoustic(Uh, Vh, level, true, true);
	
	// dU/dt = Dh^{-1} ( Rh^S V + Fh - Sh U)
	var_time.ApplyOperatorRhScalar(one, tn, Vh, zero, Prod_Uh, false);
	var_time.ApplyOperatorSh(-one, tn, Uh, one, Prod_Uh);

	if (source)
	  {
	    var_time.AddScalarSourceAtTime(one, tn, nb_deriv, Prod_Uh);
	    var_time.SetDirichletConditionSource(tn, nb_deriv+1, Prod_Uh);
	  }
	else
	  var_boundary.ImposeNullDirichletCondition(Prod_Uh);
	
	if (invert)
	  var_time.SolveOperatorDh(Prod_Uh);

	// dV/dt = Bh^{-1} ( Rh^V U - ShV V + FhV)
	var_time.ApplyOperatorRhVectorial(one, tn, Uh, zero, Prod_Vh, false);
	if (source)
	  var_time.AddVectorialSourceAtTime(one, tn, nb_deriv, Prod_Vh);
	
	var_time.ApplyOperatorShVectorial(-one, tn, Vh, one, Prod_Vh);
	
	if (invert)
	  var_time.SolveOperatorBh(Prod_Vh);    	
      }
    
    // nullify temporary vectors
    Uh.Nullify(); Vh.Nullify(); Prod_Uh.Nullify(); Prod_Vh.Nullify();
  }
  
  
  //! evaluation of Y'' = f(tn, Y, dY) (second-order formulation)
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::EvaluateFunctionS(const Real_wp& tn, const VectReal_wp& Y, const VectReal_wp& Yp,
		      VectReal_wp& ProdY, bool invert, bool source)
  {
    VectReal_wp Ydir(var_boundary.GetNbDirichletDof()), Ydir_p;
    VectReal_wp& Y_ = const_cast<VectReal_wp&>(Y);
    VectReal_wp& Yp_ = const_cast<VectReal_wp&>(Yp);
    
    // Dirichlet condition is set on Y
    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      Ydir(i) = Y(var_boundary.GetDirichletDofNumber(i));
    
    if (source)
      var_time.SetDirichletCondition(tn, 0, Y_);
    else
      var_boundary.ImposeNullDirichletCondition(Y_);
    
    // and on y'
    if (Yp.GetM() > 0)
      {
	Ydir_p.Reallocate(var_boundary.GetNbDirichletDof());
	for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	  Ydir_p(i) = Yp(var_boundary.GetDirichletDofNumber(i));
    
	if (source)
	  var_time.SetDirichletCondition(tn, 1, Yp_);
	else
	  var_boundary.ImposeNullDirichletCondition(Yp_);
      }
    
    // then we use the modified y and yp to compute y'' = Dh^{-1}(K_h y - S_h y' + F)
    Real_wp one(1), zero(0);
    var_time.ApplyOperatorKh(one, tn, Y_, zero, ProdY);
    if (source)
      var_time.AddScalarSourceAtTime(one, tn, 1, ProdY);
    
    // on ajoute la partie -Sh dY/dt
    if (Yp.GetM() > 0)
      var_time.ApplyOperatorSh(-one, tn, Yp_, one, ProdY);
    
    // on inverse par la matrice de masse si besoin
    var_boundary.ImposeNullDirichletCondition(ProdY);
    if (invert)
      var_time.SolveOperatorDh(ProdY);
    
    // stored values on dirichlet dofs are recovered
    if (Yp.GetM() > 0)
      for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	Yp_(var_boundary.GetDirichletDofNumber(i)) = Ydir_p(i);

    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      Y_(var_boundary.GetDirichletDofNumber(i)) = Ydir(i);

    // on Dirichlet dofs, we put 0
    var_boundary.ImposeNullDirichletCondition(ProdY);
    
    // alternative solution : y''_i = f'''_i on Dirichlet dofs
    //var_time.SetDirichletCondition(tn, 2, ProdY);
  }
  

  template<class Dimension>
  void TimeAcoustic_Base<Dimension>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
			     const VectReal_wp& X, VectReal_wp& Y, bool invert, bool source)
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    int N = var_problem.GetNbDof();
    int d = Dimension::dim_N;
    int nb_dof_u = Nvol + (d-1)*Npml;
    int nb_dof_v = N - nb_dof_u;
    
    VectReal_wp Uh, Vh;
    VectReal_wp Prod_Uh, Prod_Vh;
    
    Real_wp zero(0);
    
    // X is split in (Uh, Vh, Vh_pml, Vh_diamond, Vh_star)
    // idem for Y
    if ((N > X.GetM()) || (N > Y.GetM()))
      {
        cout << "X or Y is not large enough" << endl;
        DISP(N);
        DISP(X.GetM()); DISP(Y.GetM());
        abort();
      }
    
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    Uh.SetData(nb_dof_u, x_ptr); Prod_Uh.SetData(nb_dof_u, y_ptr);
    x_ptr += nb_dof_u; y_ptr += nb_dof_u;
    
    Vh.SetData(nb_dof_v, x_ptr); Prod_Vh.SetData(nb_dof_v, y_ptr);
    
    this->GetExtrapolationAcoustic(Uh, Vh, level, true, true);
    
    var_time.MltAddStiffnessScalar(alpha, level, Vh, zero, Prod_Uh);
    if (source)
      var_time.AddScalarSourceAtTime(alpha, tn, nb_deriv, Prod_Uh);
    
    if (invert)
      var_time.SolveOperatorDh(Prod_Uh);
    
    var_time.MltAddStiffnessVectorial(alpha, level, Uh, zero, Prod_Vh);
    if (source)
      var_time.AddVectorialSourceAtTime(alpha, tn, nb_deriv, Prod_Vh);
    
    if (invert)
      var_time.SolveOperatorBh(Prod_Vh);
    
    // nullify temporary vectors
    Uh.Nullify(); Vh.Nullify();
    Prod_Uh.Nullify(); Prod_Vh.Nullify();    
  }
  
  
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Vh,
			  const Real_wp& beta, VectReal_wp& Yh, bool extrapol)
  {
    //int level = MatrixVectorProductLevel::ALL_LEVELS;
    int level = var_time.GetOperatorLevel();
    if (extrapol)
      this->GetExtrapolationAcoustic(Yh, Vh, level, false, true);
    
    if (var_time.Glob_mat_RhS.GetM() > 0)
      MltAdd(alpha, var_time.Glob_mat_RhS, Vh, beta, Yh);
    else
      var_time.MltAddStiffnessScalar(alpha, level, Vh, beta, Yh);
  }
  
  
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>::
  ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			   const Real_wp& beta, VectReal_wp& Yh, bool extrapol)
  {
    int level = var_time.GetOperatorLevel();
    if (extrapol)
      this->GetExtrapolationAcoustic(Uh, Yh, level, true, false);
    
    if (var_time.Glob_mat_RhV.GetM() > 0)
      MltAdd(alpha, var_time.Glob_mat_RhV, Uh, beta, Yh);
    else
      var_time.MltAddStiffnessVectorial(alpha, level, Uh, beta, Yh);
  }
  
  
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>::
  GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y)
  {        
    // on appelle le calcul du rayonnement
    VarComputationRCS_Laplace<Dimension>& output_rcs_param =
      static_cast<VarComputationRCS_Laplace<Dimension>& >(var_output.GetParameterOutputRCS());
    
    if (output_rcs_param.GetNbPointsOutside() > 0)
      {
        // on evalue dY/dt
        VectReal_wp dY(Y.GetM()); dY.Fill(0);
        var_time.EvaluateFunction(tn, Y, dY);
        
        // on appelle le calcul du champ en des points exterieurs au domaine
        output_rcs_param.
          WriteOutput(nb_iter, tn, var_time.GetTimeStep(), Y, dY);
      }
  }
  
  
  //! writing other snapshots in file
  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::GiveVectorialIterate(int nb_iter, const Real_wp& t, VectReal_wp Vh)
  {
#ifdef MONTJOIE_WITH_TIME_REVERSAL
    if (time_reversal.IsPresent())
      time_reversal.WriteOutput(nb_iter, t, Vh);
#endif
  }


  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::MltAddStiffness(const Real_wp& alpha, int level, int lvl, const VectReal_wp& Uh, const VectReal_wp& Vh,
		    const VectReal_wp& Lambda, const Real_wp& beta, VectReal_wp& ProdUh, VectReal_wp& ProdVh)
  {
    cout << "Not available for continuous elements" << endl;
    abort();
  }


  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::GetExtrapolationLambda(const VectReal_wp& Uh, const VectReal_wp& Vh, int level, VectReal_wp& Lambda)
  {
    cout << "Not available for continuous elements" << endl;
    abort();
  }


  template<class Dimension>
  void TimeAcoustic_Base<Dimension>
  ::AddTimeSourceHDG(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& Y, VectReal_wp& Lambda)
  {
    if (t >= var_time.tlimit_source)
      return;
    
    if (var_time.SourceDoesNotDependOnTime())
      var_time.AddPrimitiveSourceAtTime(alpha, t, n, Y);
    else
      {
	VectReal_wp b_src(var_problem.GetNbDof());
	b_src.Zero();
	
	Vector<VectReal_wp> b_vec;
        b_vec.SetData(1, &b_src);
        
	VirtualSourceFEM<Real_wp, Dimension>* f;
	f = var_source.GetNewSourceEquationObject(0);
	
	f->Init(t, var_time.GetTimeStep(), var_problem.print_level, n, true);
	
	Vector<VirtualSourceFEM<Real_wp, Dimension>* > f_vec(1);
	f_vec(0) = f;
	var_source.AddSurfaceSource(Real_wp(1), b_vec, f_vec);
        b_vec.Nullify();
	
	int nb_dof_L = var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
	int nb_dof_uv = var_problem.GetNbDof() - nb_dof_L;
	
	for (int i = 0; i < nb_dof_L; i++)
	  Lambda(i) += alpha*b_src(i);

	for (int i = 0; i < nb_dof_uv; i++)
	  Y(i) += alpha*b_src(nb_dof_L + i);
      }
  }
  
  
  /***********************************
   * DiscontinuousMassMatrixAcoustic *
   ***********************************/
  
  
  //! computes the mass matrices
  template<class Dimension>
  void DiscontinuousMassMatrixAcoustic<Dimension>::Init(bool compute_time)
  {
    Vector<bool> diag_elt;
    int type = var_problem.GetMassMatrixType(diag_elt);
    bool diag_matrices = true;
    if (type != FemMassMatrix::DIAGONAL)
      diag_matrices = false;

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    int N = var_problem.GetNbDof();
    int d = Dimension::dim_N;
    int nb_dof_u = Nvol + (d-1)*Npml;
    int nb_dof_v = N - nb_dof_u;
    Matrix<Real_wp, Symmetric, DiagonalRow> Dh(N, N), Sh(N, N);    
    
    // evaluation of mass matrix (diagonal part)
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(1.0);
    nat_mat.SetCoefDamping(0.0);
    nat_mat.SetCoefStiffness(0.0);
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    var_problem.AddMatrixWithBC(Dh, nat_mat);
    
    // and damping matrix (diagonal part)
    nat_mat.SetCoefMass(0.0);
    nat_mat.SetCoefDamping(1.0);
    nat_mat.SetCoefStiffness(0.0);
    var_boundary.SetCoefficientDirichlet(Real_wp(0));
    var_problem.AddMatrixWithBC(Sh, nat_mat);
        
    bool compute_Dh, compute_invDh, compute_DhMinusdtSh,
      compute_invDhPlusdtSh, compute_Sh;
    
    this->FindMatricesToCompute(compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh);
    
    bool compute_Bh, compute_invBh, compute_BhMinusdtSh,
      compute_invBhPlusdtSh, compute_ShVec;
    
    this->FindMatricesToComputeVec(compute_Bh, compute_invBh, compute_BhMinusdtSh,
				   compute_invBhPlusdtSh, compute_ShVec);
    
    Real_wp deltat = var_time.GetTimeStep();
    int offset_u = 0, offset_v = nb_dof_u;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	offset_u = mesh_num.GetNbDof();
	nb_dof_u = var_problem.GetNbMainUnknownDof();
	nb_dof_v = nb_dof_u*Dimension::dim_N;
	offset_v = offset_u + nb_dof_u;
      }

    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    if (diag_matrices)
      {
	if (compute_Dh || compute_invDh)
	  {
	    if (this->Dh != NULL)
	      delete this->Dh;
	    
	    this->Dh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->Dh->SetCoefficient(1.0, 0.0);
	  }
	
	if (compute_Dh)
	  {
	    VectReal_wp& diagonal = this->Dh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_u);
	    for (int i = 0; i < nb_dof_u; i++)
	      diagonal(i) = Dh(offset_u + i, offset_u + i);
	  }
	
	if (compute_invDh)
	  {
	    VectReal_wp& invDiagonal = this->Dh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_u);
	    for (int i = 0; i < nb_dof_u; i++)
	      invDiagonal(i) = 1.0/Dh(offset_u + i, offset_u + i);
	  }
	
	if (compute_DhMinusdtSh)
	  {
	    if (this->DhMinusdtSh != NULL)
	      delete this->DhMinusdtSh;
	    
	    this->DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->DhMinusdtSh->SetCoefficient(1.0, -0.5*var_time.GetTimeStep());
	    
	    VectReal_wp& diagonal = this->DhMinusdtSh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_u);
	    for (int i = 0; i < nb_dof_u; i++)
	      diagonal(i) = Dh(offset_u + i, offset_u + i) - 0.5*deltat*Sh(offset_u+i, offset_u+i);
	  }
	
	if (compute_invDhPlusdtSh)
	  {
	    if (this->DhPlusdtSh != NULL)
	      delete this->DhPlusdtSh;
	    
	    this->DhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->DhPlusdtSh->SetCoefficient(1.0, 0.5*var_time.GetTimeStep());

	    VectReal_wp& invDiagonal = this->DhPlusdtSh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_u);
	    for (int i = 0; i < nb_dof_u; i++)
	      invDiagonal(i) = 1.0 / (Dh(offset_u+i, offset_u+i) + 0.5*deltat*Sh(offset_u+i, offset_u+i));
	  }
	
	if (compute_Sh)
	  {
	    if (this->Sh != NULL)
	      delete this->Sh;
	    
	    this->Sh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->Sh->SetCoefficient(0.0, 1.0);	
	    
	    VectReal_wp& diagonal = this->Sh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_u);
	    for (int i = 0; i < nb_dof_u; i++)
	      diagonal(i) = Sh(offset_u+i, offset_u+i);
	  }
	
	if ((compute_Bh) || (compute_invBh))
	  {
	    if (this->Bh != NULL)
	      delete this->Bh;

	    this->Bh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->Bh->SetCoefficient(1.0, 0.0);	
	  }

	if (compute_Bh)
	  {
	    VectReal_wp& diagonal = this->Bh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_v);
	    for (int i = 0; i < nb_dof_v; i++)
	      diagonal(i) = Dh(offset_v+i, offset_v+i);          
	  }
	
	if (compute_invBh)
	  {
	    VectReal_wp& invDiagonal = this->Bh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_v);
	    for (int i = 0; i < nb_dof_v; i++)
	      invDiagonal(i) = 1.0/Dh(offset_v+i, offset_v+i);          
	  }
	
	if (compute_ShVec)
	  {
	    if (this->ShVec != NULL)
	      delete this->ShVec;

	    this->ShVec = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->ShVec->SetCoefficient(0.0, 1.0);	
	    VectReal_wp& diagonal = this->ShVec->GetDiagonal();
	    diagonal.Reallocate(nb_dof_v);
	    for (int i = 0; i < nb_dof_v; i++)
	      diagonal(i) = Sh(offset_v+i, offset_v+i);          
	  }
	
	if (compute_BhMinusdtSh)
	  {
	    if (this->BhMinusdtSh != NULL)
	      delete this->BhMinusdtSh;

	    this->BhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->BhMinusdtSh->SetCoefficient(1.0, -0.5*var_time.GetTimeStep());

	    VectReal_wp& diagonal = this->BhMinusdtSh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_v);
	    for (int i = 0; i < nb_dof_v; i++)
	      diagonal(i) = Dh(offset_v+i, offset_v+i)
		- 0.5*deltat*Sh(offset_v+i, offset_v+i);          
	  }
	
	if (compute_invBhPlusdtSh)
	  {
	    if (this->BhPlusdtSh != NULL)
	      delete this->BhPlusdtSh;

	    this->BhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    this->BhPlusdtSh->SetCoefficient(1.0, 0.5*var_time.GetTimeStep());

	    VectReal_wp& invDiagonal = this->BhPlusdtSh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_v);
	    for (int i = 0; i < nb_dof_v; i++)
	      invDiagonal(i)
                = 1.0/(Dh(offset_v+i, offset_v+i) + 0.5*deltat*Sh(offset_v+i, offset_v+i));
	  }

	DiscontinuousUnsteadyMassMatrix<Dimension>::InitSurface();
      }
    else
      {
	// premiere etape : on calcule les matrices comme sans les PMLs
	DiscontinuousUnsteadyMassMatrix<Dimension>::Init(false);
	
	// et on rajoute les PMLs
	if (mesh_num.GetNbDofPML() > 0)
	  {
	    int d = Dimension::dim_N;
	    int nb_scal = 1;
	    int nb_vec = Dimension::dim_N+1;
	    int nodl = var_problem.offset_dof_unknown(1);
	    int nb_new = (d-1)*var_boundary.GetNbEltPML();

	    // on verifie qu'il y a condensation de masse sur les PML
	    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
	      if (var_problem.InsidePML(i))
		{
		  if (!var_problem.GetReferenceElement(i).LumpedMassMatrix())
		    {
		      cout << "not implemented" << endl;
		      abort();
		    }
		}
	    
	    if (compute_Dh || compute_invDh)
	      {
		VectReal_wp diag_Dh(nb_dof_u);
		for (int i = 0; i < nb_dof_u; i++)
		  diag_Dh(i) = Dh(i, i);
                
		this->Dh->SetDiagonalPML(nb_new, 0, nb_scal, nodl, diag_Dh);
	      }
	    
	    if (compute_DhMinusdtSh)
	      {
		VectReal_wp diag_Dh(nb_dof_u);
		for (int i = 0; i < nb_dof_u; i++)
		  diag_Dh(i) = Dh(i, i) - 0.5*deltat*Sh(i, i);
		
		this->DhMinusdtSh
		  ->SetDiagonalPML(nb_new, 0, nb_scal, nodl, diag_Dh);
	      }
	    
	    if (compute_invDhPlusdtSh)
	      {
		VectReal_wp diag_Dh(nb_dof_u);
		for (int i = 0; i < nb_dof_u; i++)
		  diag_Dh(i) = (Dh(i, i) + 0.5*deltat*Sh(i, i));
		
		this->DhPlusdtSh
		  ->SetDiagonalPML(nb_new, 0, nb_scal, nodl, diag_Dh);
	      }
	
	    if (compute_Sh)
	      {
		VectReal_wp diag_Dh(nb_dof_u);
		for (int i = 0; i < nb_dof_u; i++)
		  diag_Dh(i) = Sh(i, i);
		
	      this->Sh
		->SetDiagonalPML(nb_new, 0, nb_scal, nodl, diag_Dh);
	      }
	    
	    if (compute_Bh)
	      {
		VectReal_wp diag_Dh(nb_dof_v);
		for (int i = 0; i < nb_dof_v; i++)
		  diag_Dh(i) = Dh(nb_dof_u+i, nb_dof_u+i);
		
		this->Bh->SetDiagonalPML(0, nb_scal, nb_vec, Nvol, diag_Dh);
	      }
	    
	    if (compute_ShVec)
	      {
		VectReal_wp diag_Dh(nb_dof_v);
		for (int i = 0; i < nb_dof_v; i++)
		  diag_Dh(i) = Sh(nb_dof_u+i, nb_dof_u+i);
		  
		this->ShVec->SetDiagonalPML(0, nb_scal, nb_vec, Nvol, diag_Dh);
	      }
	    
	    if (compute_BhMinusdtSh)
	      {
		VectReal_wp diag_Dh(nb_dof_v);
		for (int i = 0; i < nb_dof_v; i++)
		  diag_Dh(i) = Dh(nb_dof_u+i, nb_dof_u+i) 
		    - 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i); 
	      
		this->BhMinusdtSh
		  ->SetDiagonalPML(0, nb_scal, nb_vec, Nvol, diag_Dh);
	      }
	    
	    if (compute_invBhPlusdtSh)
	      {
		VectReal_wp diag_Dh(nb_dof_v);
		for (int i = 0; i < nb_dof_v; i++)
		  diag_Dh(i) = (Dh(nb_dof_u+i, nb_dof_u+i)
                                + 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i));
		
		this->BhPlusdtSh
		  ->SetDiagonalPML(0, nb_scal, nb_vec, Nvol, diag_Dh);
	      }
	  }
      }
  }


  /********************************
   * ContinuousMassMatrixAcoustic *
   ********************************/
  

  //! Computes the mass matrices
  template<class Dimension>
  void ContinuousMassMatrixAcoustic<Dimension>::Init(bool compute_time)
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    int N = var_problem.GetNbDof();
    int d = Dimension::dim_N;
    int nb_dof_u = Nvol + (d-1)*Npml;
    int nb_dof_v = N - nb_dof_u;
    Matrix<Real_wp, Symmetric, DiagonalRow> Dh(N, N);    

    int Npml_all = Npml;
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&Npml, &Npml_all, 1, MPI_INTEGER, MPI_MAX, var_problem.comm_group_mode);
#endif

    Vector<bool> diag_elt;
    int type = var_problem.GetMassMatrixType(diag_elt);
    bool diag_matrices = true;
    if (type != FemMassMatrix::DIAGONAL)
      diag_matrices = false;

    bool compute_Dh, compute_invDh, compute_DhMinusdtSh,
      compute_invDhPlusdtSh, compute_Sh;
    
    this->FindMatricesToCompute(compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh);
    
    bool compute_Bh, compute_invBh, compute_BhMinusdtSh,
      compute_invBhPlusdtSh, compute_ShVec;
    
    this->FindMatricesToComputeVec(compute_Bh, compute_invBh, compute_BhMinusdtSh,
				   compute_invBhPlusdtSh, compute_ShVec);
    
    Real_wp deltat = var_time.GetTimeStep();
    if (!diag_matrices)
      {        
        ContinuousUnsteadyMassMatrix<Dimension>::Init(compute_time);
        return;
      }

    // evaluation of mass matrix
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(1.0);
    nat_mat.SetCoefDamping(0.0);
    nat_mat.SetCoefStiffness(0.0);
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    var_problem.AddMatrixWithBC(Dh, nat_mat);
    
    if ((compute_Dh) || (compute_invDh))
      {
	if (this->Dh != NULL)
	  delete this->Dh;
	
	this->Dh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
      }

    if (compute_Dh)
      {
	VectReal_wp& diagonal = this->Dh->GetDiagonal();
        diagonal.Reallocate(nb_dof_u);
        for (int i = 0; i < nb_dof_u; i++)
          diagonal(i) = Dh(i, i);
      }
    
    if (compute_invDh)
      {
	VectReal_wp& invDiagonal = this->Dh->GetInverseDiagonal();
        invDiagonal.Reallocate(nb_dof_u);
        for (int i = 0; i < nb_dof_u; i++)
          invDiagonal(i) = Dh(i, i);
        
        var_problem.AddDomains(invDiagonal);
        for (int i = 0; i < nb_dof_u; i++)
          invDiagonal(i) = 1.0/invDiagonal(i);
      }
    
    if (Npml_all == 0)
      {
        Matrix<Real_wp, Symmetric, DiagonalRow> Sh(N, N);            
        // in that case, mass matrix is diagonal
        nat_mat.SetCoefMass(0.0);
        nat_mat.SetCoefDamping(1.0);
        nat_mat.SetCoefStiffness(0.0);
	var_boundary.SetCoefficientDirichlet(Real_wp(0));
        var_problem.AddMatrixWithBC(Sh, nat_mat);
                
        if (compute_DhMinusdtSh)
          {
	    if (this->DhMinusdtSh != NULL)
	      delete this->DhMinusdtSh;

            this->DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    
	    VectReal_wp& diagonal = this->DhMinusdtSh->GetDiagonal();
            diagonal.Reallocate(nb_dof_u);
            for (int i = 0; i < nb_dof_u; i++)
              diagonal(i) = Dh(i, i) - 0.5*deltat*Sh(i, i);
          }
        
        if (compute_invDhPlusdtSh)
          {
	    if (this->DhPlusdtSh != NULL)
	      delete this->DhPlusdtSh;

            this->DhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& invDiagonal = this->DhPlusdtSh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_u);
            for (int i = 0; i < nb_dof_u; i++)
              invDiagonal(i) = Dh(i, i) + 0.5*deltat*Sh(i, i);
            
            var_problem.AddDomains(invDiagonal);
            for (int i = 0; i < nb_dof_u; i++)
              invDiagonal(i) = 1.0/invDiagonal(i);
          }
	
        if (compute_Sh)
          {
	    if (this->Sh != NULL)
	      delete this->Sh;

	    this->Sh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
            VectReal_wp& diagonal = this->Sh->GetDiagonal();
	    diagonal.Reallocate(nb_dof_u);
            for (int i = 0; i < nb_dof_u; i++)
              diagonal(i) = Sh(i, i);
          }
        
        if ((compute_Bh) || (compute_invBh))
	  {
	    if (this->Bh != NULL)
	      delete this->Bh;

	    this->Bh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	  }

	if (compute_Bh)
	  {
	    VectReal_wp& diagonal = this->Bh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i);          
          }
	
        if (compute_invBh)
          {
	    VectReal_wp& invDiagonal = this->Bh->GetInverseDiagonal();
            invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i) = 1.0/Dh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_ShVec)
          {
	    if (this->ShVec != NULL)
	      delete this->ShVec;

	    this->ShVec = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->ShVec->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Sh(nb_dof_u+i, nb_dof_u+i);          
          }
	
        if (compute_BhMinusdtSh)
          {
	    if (this->BhMinusdtSh != NULL)
	      delete this->BhMinusdtSh;

	    this->BhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->BhMinusdtSh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i) - 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i);          
          }
	
        if (compute_invBhPlusdtSh)
          {
	    if (this->BhPlusdtSh != NULL)
	      delete this->BhPlusdtSh;

	    this->BhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& invDiagonal = this->BhPlusdtSh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i)
                = 1.0/(Dh(nb_dof_u+i, nb_dof_u+i) + 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i));
          }
      }
    else
      {
        Vector<IVect> pattern(Nvol + nb_dof_v);
        DistributedMatrix<Real_wp, General, BlockDiagRow> Sh;
        
        for (int i = 0; i < Nvol; i++)
          {
            int num = mesh_num.GetDofPML(i);
            if (num >= 0)
              {
                pattern(i).Reallocate(d);
                pattern(i)(0) = i;
                pattern(i)(1) = Nvol + num;
                if (d == 3)
                  pattern(i)(2) = Nvol + Npml + num;
              }
            else
              {
                pattern(i).Reallocate(1);
                pattern(i)(0) = i;
              }
          }
        
        for (int i = 0; i < nb_dof_v; i++)
          {
            pattern(Nvol+i).Reallocate(1);
            pattern(Nvol+i)(0) = nb_dof_u+i;
          }
        
        var_problem.InitDistributedMatrix(Sh);
        Sh.SetPattern(pattern);
        pattern.Resize(Nvol);
        
        nat_mat.SetCoefMass(0.0);
        nat_mat.SetCoefDamping(1.0);
        nat_mat.SetCoefStiffness(0.0);
	var_boundary.SetCoefficientDirichlet(Real_wp(0));
        var_problem.AddMatrixWithBC(Sh, nat_mat);
	
	if (this->Sh != NULL)
	  delete this->Sh;

	this->Sh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_DIAGONAL_UNSYM);
	DistributedMatrix<Real_wp, General, BlockDiagRow>&
	  Sh_block_diagonal = this->Sh->GetUnsymmetricBlockDiagonal();
	
        var_problem.InitDistributedMatrix(Sh_block_diagonal);
        Sh_block_diagonal.SetPattern(pattern);
        for (int i = 0; i < Nvol; i++)
          {
            int n = Sh_block_diagonal.GetBlockSize(i);
            for (int j = 0; j < n; j++)
              for (int k = 0; k < n; k++)
                Sh_block_diagonal.Value(i, j, k) = Sh.Value(i, j, k);
          }
	
        if (compute_DhMinusdtSh)
          {
	    if (this->DhMinusdtSh != NULL)
	      delete this->DhMinusdtSh;

	    this->DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_DIAGONAL_UNSYM);
	    DistributedMatrix<Real_wp, General, BlockDiagRow>&
	      block_diagonal = this->DhMinusdtSh->GetUnsymmetricBlockDiagonal();

            block_diagonal = Sh_block_diagonal;            
            Mlt(-0.5*deltat, block_diagonal);
            for (int i = 0; i < nb_dof_u; i++)
	      block_diagonal.AddInteraction(i, i, Dh(i, i));          
          }
        
        if (compute_invDhPlusdtSh)
          {
	    if (this->DhPlusdtSh != NULL)
	      delete this->DhPlusdtSh;

	    this->DhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::BLOCK_DIAGONAL_UNSYM);
	    DistributedMatrix<Real_wp, General, BlockDiagRow>&
	      block_diagonal = this->DhPlusdtSh->GetUnsymmetricBlockDiagonal();
	    
            block_diagonal = Sh_block_diagonal;            
            Mlt(0.5*deltat, block_diagonal);
            for (int i = 0; i < nb_dof_u; i++)
	      block_diagonal.AddInteraction(i, i, Dh(i, i));          
            
            this->DhPlusdtSh->Invert(false);
          }
        
        if (!compute_Sh)
	  Sh_block_diagonal.Clear();
	
        if ((compute_Bh) || (compute_invBh))
	  {
	    if (this->Bh != NULL)
	      delete this->Bh;

	    this->Bh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	  }

        if (compute_Bh)
          {
	    VectReal_wp& diagonal = this->Bh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_invBh)
          {
	    VectReal_wp& invDiagonal = this->Bh->GetInverseDiagonal();
	    invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i) = 1.0/Dh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_ShVec)
          {
	    if (this->ShVec != NULL)
	      delete this->ShVec;

	    this->ShVec = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->ShVec->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Sh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_BhMinusdtSh)
          {
	    if (this->BhMinusdtSh != NULL)
	      delete this->BhMinusdtSh;

	    this->BhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& diagonal = this->BhMinusdtSh->GetDiagonal();
            diagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              diagonal(i) = Dh(nb_dof_u+i, nb_dof_u+i) - 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i);          
          }
        
        if (compute_invBhPlusdtSh)
          {
	    if (this->BhPlusdtSh != NULL)
	      delete this->BhPlusdtSh;
	    
	    this->BhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
	    VectReal_wp& invDiagonal = this->BhPlusdtSh->GetInverseDiagonal();
            invDiagonal.Reallocate(nb_dof_v);
            for (int i = 0; i < nb_dof_v; i++)
              invDiagonal(i) = 1.0/
		(Dh(nb_dof_u+i, nb_dof_u+i) + 0.5*deltat*Sh(nb_dof_u+i, nb_dof_u+i));
          }
      }

    var_boundary.SetCoefficientDirichlet(Real_wp(1));
  }
  
  
  /*******************
   * TimeAcoustic_Eq *
   *******************/
  
  
  template<class TypeEquation>
  TimeAcoustic_Eq<TypeEquation>::TimeAcoustic_Eq()
    : TimeAcoustic_Base<Dimension>(this->GetLeafClass())
  {
    this->var_harmonic.compute_dfjm1 = true;
  }


  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>::GetMemoryUsed(map<string, size_t>& var) const
  {
    VarInstationary<TypeEquation>::GetMemoryUsed(var);
    TimeAcoustic_Base<Dimension>::GetMemoryUsed(var);
  }

  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>::InitTimeIterations()
  {
    // use of split PML => jacobian matrices DFi must be diagonal
    this->CheckIdentityJacobianInPML();
    
    VarInstationary<TypeEquation>::InitTimeIterations();
    TimeAcoustic_Base<Dimension>::InitTimeIterations();
  }
  

  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
    ::AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y)
  {
#ifdef MONTJOIE_WITH_TIME_REVERSAL
    if (this->time_reversal.IsPresent() && 
	(this->time_reversal.GetSimulationType() == this->time_reversal.INVERSE))
      return this->time_reversal.AddScalarSourceAtTime(alpha, tn, n, Y);
#endif
    
    VarInstationary<TypeEquation>::AddScalarSourceAtTime(alpha, tn, n, Y);
  }
  

  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::AddVectorialSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y)
  {
#ifdef MONTJOIE_WITH_TIME_REVERSAL
    if (this->time_reversal.IsPresent() &&
	(this->time_reversal.GetSimulationType() == this->time_reversal.INVERSE))
      return;
#endif
    
    VarInstationary<TypeEquation>::AddVectorialSourceAtTime(alpha, tn, n, Y);
  }
  
  
  //! reads a line of the data file
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>::SetInputData(const string& description_field, const VectString& parameters)
  {    
    VarInstationary<TypeEquation>::SetInputData(description_field, parameters);  
    TimeAcoustic_Base<Dimension>::SetInputData(description_field, parameters);
  }
  
  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
			       const VectReal_wp& X, VectReal_wp& Y,
			       bool invert, bool source)
  {
    TimeAcoustic_Base<Dimension>::EvaluateDerivativeFunction(tn, nb_deriv, X, Y, invert, source);
  }
    
  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::EvaluateFunctionS(const Real_wp& tn,const VectReal_wp& Y, const VectReal_wp& Yp,
		      VectReal_wp& ProdY, bool invert, bool source)
  {
    TimeAcoustic_Base<Dimension>::EvaluateFunctionS(tn, Y, Yp, ProdY, invert, source);
  }
  

  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::EvaluateFunction(const Real_wp& tn,const VectReal_wp& Y,
                     VectReal_wp& ProdY, bool invert, bool source)
  {
    VarInstationary<TypeEquation>::EvaluateFunction(tn, Y, ProdY, invert, source);
  }


  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::EvaluateFunction(const Real_wp& tn, const Real_wp& alpha, int level,
                          const VectReal_wp& X, VectReal_wp& Y,
                          bool invert_mass, bool source)
  {
    VarInstationary<TypeEquation>::EvaluateFunction(tn, alpha, level, X, Y, invert_mass, source);
  }

  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
			       const VectReal_wp& X, VectReal_wp& Y,
			       bool invert, bool source)
  {
    TimeAcoustic_Base<Dimension>::EvaluateDerivativeFunction(tn, nb_deriv, alpha, level, X, Y, invert, source);
  }
  
  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Vh,
			  const Real_wp& beta, VectReal_wp& Yh, bool extrapol)
  {
    TimeAcoustic_Base<Dimension>::ApplyOperatorRhScalar(alpha, t, Vh, beta, Yh, extrapol);
  }
  
  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
			     const Real_wp& beta, VectReal_wp& Yh, bool extrapol)
  {
    TimeAcoustic_Base<Dimension>::ApplyOperatorRhVectorial(alpha, t, Uh, beta, Yh, extrapol);
  }
  
  
  //! treats a new iterate
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>::GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y)
  {
    VarInstationary_Base::GiveIterate(nb_iter, tn, Y);
    TimeAcoustic_Base<Dimension>::GiveIterate(nb_iter, tn, Y);
  }


  //! treats a new vectorial iterate
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>::GiveVectorialIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y)
  {
    VarInstationary_Base::GiveVectorialIterate(nb_iter, tn, Y);
    TimeAcoustic_Base<Dimension>::GiveVectorialIterate(nb_iter, tn, Y);
  }
  
  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>::
  ComputeExtrapolationVectorial(const VectReal_wp& Uh, const VectReal_wp& Vh)
  {
    int level = this->global_level_operator;
    this->GetExtrapolationAcoustic(Uh, Vh, level, false, true);
  }
  
  
  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::ComputeExtrapolationScalar(const VectReal_wp& Uh, const VectReal_wp& Vh)
  {
    int level = this->global_level_operator;
    this->GetExtrapolationAcoustic(Uh, Vh, level, true, false);    
  }    


  template<class TypeEquation>
  void TimeAcoustic_Eq<TypeEquation>
  ::GetExtrapolationAcoustic(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
			     bool extrapol_u, bool extrapol_v)
  {
    if (this->Glob_mat_RhV.GetM() > 0)
      return;
    
    Montjoie::GetExtrapolationAcoustic(this->GetLeafClass(), U, V, lvl, 
				       this->extrapolU, this->extrapolV, this->Uneighbor, this->Vneighbor,
				       this->evalUn_quad, this->evalVn_quad, extrapol_u, extrapol_v);
  }

  template<class TypeEquation>
  bool TimeAcoustic_Eq<TypeEquation>::SplitSystem() const 
  {
    return true; 
  }
  
  
  template<class TypeEquation>
  VirtualMassMatrix* TimeAcoustic_Eq<TypeEquation>::GetNewMassMatrix()
  {
    if (this->var_harmonic.FormulationDG() != ElementReference_Base::CONTINUOUS)
      return new DiscontinuousMassMatrixAcoustic<Dimension>(this->GetLeafClass());
    
    return new ContinuousMassMatrixAcoustic<Dimension>(this->GetLeafClass());
  }


  /*********************
   * HyperbolicProblem *
   *********************/
  
  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquation<Dimension> >
  ::AddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha_,
                              const VectReal_wp& Uh, VectReal_wp& Vh,
                              const Real_wp& coef_mu, const Real_wp& coef_sigma)
  {
    const ElementReference<Dimension, 1>& Fb = this->var_harmonic.GetReferenceElementH1(i);
    int ref = this->var_harmonic.mesh.Element(i).GetReference();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    const VectReal_wp& WeightsND = Fb.WeightsND();
    Real_wp coef; typename Dimension::R_N sigma, epsilon;
    bool affine = this->var_harmonic.mesh.IsElementAffine(i);
    Real_wp alpha = coef_mu*alpha_;
    if (coef_mu != Real_wp(0))
      for (int j = 0; j < nb_points_quad; j++)
        {
          if (affine)
            coef = this->var_harmonic.Glob_jacobian(i)(0)*WeightsND(j)*alpha;
          else
            coef = this->var_harmonic.Glob_jacobian(i)(j)*alpha;
          
          sigma(0) = Uh(offset);
          sigma(1) = Uh(offset+1);
          sigma(2) = Uh(offset+2);
          
          this->var_harmonic.ref_invMu(ref).MltMatrix(this->var_harmonic, i, j, sigma, epsilon);
          
          Vh(offset) += coef*epsilon(0);
          Vh(offset+1) += coef*epsilon(1);
          Vh(offset+2) += coef*epsilon(2);
          
          offset += 3;
        }
    else
      offset += 3*nb_points_quad;
    
    if (this->var_harmonic.InsidePML(i))
      {
        cout << "Not implemented" << endl;
        abort();
      }
  }

  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquation<Dimension> >
  ::SolveElementStiffnessTensorGen(int i, int& offset, const Real_wp& alpha,
                                   const VectReal_wp& Uh, VectReal_wp& Vh,
                                   const Real_wp& coef_sigma, bool add)
  {
    const ElementReference<Dimension, 1>& Fb = this->var_harmonic.GetReferenceElementH1(i);
    int ref = this->var_harmonic.mesh.Element(i).GetReference();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    const VectReal_wp& WeightsND = Fb.WeightsND();
    Real_wp coef; typename Dimension::R_N sigma, epsilon;
    bool affine = this->var_harmonic.mesh.IsElementAffine(i);
    if (alpha != Real_wp(0))
      for (int j = 0; j < nb_points_quad; j++)
        {
          if (affine)
            coef = alpha/(this->var_harmonic.Glob_jacobian(i)(0)*WeightsND(j));
          else
            coef = alpha / this->var_harmonic.Glob_jacobian(i)(j);
          
          sigma(0) = Uh(offset);
          sigma(1) = Uh(offset+1);
          sigma(2) = Uh(offset+2);
          
          this->var_harmonic.ref_mu(ref).MltMatrix(this->var_harmonic, i, j, sigma, epsilon);
          
          if (add)
            {
              Vh(offset) += coef*epsilon(0);
              Vh(offset+1) += coef*epsilon(1);
              Vh(offset+2) += coef*epsilon(2);
            }
          else
            {
              Vh(offset) = coef*epsilon(0);
              Vh(offset+1) = coef*epsilon(1);
              Vh(offset+2) = coef*epsilon(2);
            }

              
          offset += 3;
        }
    else
      offset += 3*nb_points_quad;
    
    if (this->var_harmonic.InsidePML(i))
      {
        cout << "Not implemented" << endl;
        abort();
      }
  }
  
  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquation<Dimension> >
  ::SolveElementStiffnessTensor(int i, int& offset, VectReal_wp& Vh,
                                const Real_wp& coef_sigma)
  {
    SolveElementStiffnessTensorGen(i, offset, Real_wp(1), Vh, Vh, coef_sigma, false);
  }
    
  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquation<Dimension> >
  ::SolveAddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
                                   const VectReal_wp& Uh, VectReal_wp& Vh,
                                   const Real_wp& coef_sigma)
  {
    SolveElementStiffnessTensorGen(i, offset, alpha, Uh, Vh, coef_sigma, true);
  }
  

  template<class Dimension>
  void HyperbolicProblem<AcousticEquation<Dimension> >
  ::MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			  const Real_wp& beta, VectReal_wp& C)
  {
    FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >&
      Kh = static_cast<FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >& >(*this->Glob_mat_Kh);

    MltAdd_SquareHex_ScalarH1(alpha, level, Kh, B, beta, C);
  }

  template<class Dimension>
  void HyperbolicProblem<AcousticEquation<Dimension> >
  ::MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
			     const Real_wp& beta, VectReal_wp& C)
  {
    FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >&
      Kh = static_cast<FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >& >(*this->Glob_mat_Kh);

    MltAdd_SquareHex_VectorialH1(alpha, level, Kh, B, beta, C);
  }

  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			  const Real_wp& beta, VectReal_wp& C)
  {
    MltAddStiffnessScalar_AcousticDG(alpha, level, *this, B, beta, C);
  }

  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
			     const Real_wp& beta, VectReal_wp& C)
  {
    MltAddStiffnessVectorial_AcousticDG(alpha, level, *this, B, beta, C);
  }


  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::MltAddStiffness(const Real_wp& alpha, int level, int lvl, const VectReal_wp& Uh, const VectReal_wp& Vh,
		    const VectReal_wp& Lambda, const Real_wp& beta, VectReal_wp& ProdUh, VectReal_wp& ProdVh)
  {
    MltAddStiffnessAcousticHDG(alpha, level, lvl, *this, Uh, Vh, Lambda,
			       this->evalUn_quad, this->evalVn_quad,
			       this->extrapolU, this->extrapolV, beta, ProdUh, ProdVh);
  }


  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::GetExtrapolationLambda(const VectReal_wp& Uh, const VectReal_wp& Vh, int level, VectReal_wp& Lambda)
  {
    GetExtrapolationAcousticHDG(*this, level, Uh, Vh, this->evalUn_quad, this->evalVn_quad,
				this->extrapolU, this->extrapolV, Lambda);
  }
  

  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma)
  {
    // isotropic case only
    if (num == 0)
      {
	rho = this->var_harmonic.ref_rho(ref).GetCoefficient(this->var_harmonic, i, j);
	sigma = this->var_harmonic.ref_sigma(ref).GetCoefficient(this->var_harmonic, i, j);
      }
    else
      {
	rho = this->var_harmonic.ref_invMu(ref).GetCoefficient(this->var_harmonic, i, j, num-1, num-1);
	sigma = 0.0;
      }
  }


  template<class Dimension>
  Real_wp HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::GetCoefficientTauHDG(int num_face) const
  {
    return this->var_harmonic.alpha_penalization*this->var_harmonic.Glob_CoefPenalDG(num_face);
  }


  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t, int nb_deriv, VectReal_wp& b_src)
  {
    if ((this->var_harmonic.FormulationDG() == ElementReference_Base::HDG)
	&& (!this->SourceDoesNotDependOnTime()))
      {
	if (t >= this->tlimit_source)
	  return;

	int nb_dof_L = this->var_harmonic.GetNbMeshDof();
	VectReal_wp Lambda(nb_dof_L);
	Lambda.Zero();
	this->AddTimeSourceHDG(alpha, t, nb_deriv, b_src, Lambda);
	
	this->Glob_mat_Dh->SolveOperatorCh(Lambda);

	if (this->Glob_mat_Ch_Lambda.GetM() > 0)
	  MltAdd(Real_wp(-1), this->Glob_mat_Ch_Lambda, Lambda, Real_wp(1), b_src);
	else
	  MltAddLambdaAcousticHDG(Real_wp(1), -1, *this, Lambda, b_src);
      }
    else
      TimeAcoustic_Eq<AcousticEquationDG<Dimension> >::AddPrimitiveTimeSource(alpha, t, nb_deriv, b_src);
  }


  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::GetExtrapolationLambdaHDG(int level, const VectReal_wp& Y, VectReal_wp& ProdL)
  {
    int Nvol = this->var_harmonic.GetNbMainUnknownDof();
    int d = Dimension::dim_N;
    int nb_dof_u = Nvol;
    int nb_dof_v = Nvol*d;

    VectReal_wp Uh, Vh;
    Real_wp* y_ptr = Y.GetData();

    Uh.SetData(nb_dof_u, y_ptr); y_ptr += nb_dof_u;
    Vh.SetData(nb_dof_v, y_ptr);

    this->GetExtrapolationLambda(Uh, Vh, level, ProdL);

    Uh.Nullify(); Vh.Nullify();
  }

  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::MltAddStiffnessHDG(const Real_wp& alpha, int level, int level_inside,
		       const VectReal_wp& Y, const VectReal_wp& L,
		       const Real_wp& beta, VectReal_wp& ProdY)
  {
    int Nvol = this->var_harmonic.GetNbMainUnknownDof();
    int d = Dimension::dim_N;
    int nb_dof_u = Nvol;
    int nb_dof_v = Nvol*d;

    VectReal_wp Uh, Vh, Prod_Uh, Prod_Vh;
    Real_wp* y_ptr = Y.GetData(); Real_wp* py_ptr = ProdY.GetData();

    Uh.SetData(nb_dof_u, y_ptr); y_ptr += nb_dof_u;
    Vh.SetData(nb_dof_v, y_ptr);

    Prod_Uh.SetData(nb_dof_u, py_ptr); py_ptr += nb_dof_u;
    Prod_Vh.SetData(nb_dof_v, py_ptr);
    
    this->MltAddStiffness(alpha, level, level_inside, Uh, Vh, L, beta, Prod_Uh, Prod_Vh);

    Uh.Nullify(); Vh.Nullify();
    Prod_Uh.Nullify(); Prod_Vh.Nullify();
  }
  
  
  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::MltAddLambdaHDG(const Real_wp& alpha, int level,
		    const VectReal_wp& L, VectReal_wp& ProdY)
  {
    MltAddLambdaAcousticHDG(alpha, level, *this, L, ProdY);
  }


  template<class Dimension>
  void HyperbolicProblem<AcousticEquationDG<Dimension> >
  ::MltAddLambdaHDG(const Complex_wp& alpha, int level,
		    const VectComplex_wp& L, VectComplex_wp& ProdY)
  {
    MltAddLambdaAcousticHDG(alpha, level, *this, L, ProdY);
  }

}

#define MONTJOIE_FILE_TIME_ACOUSTIC_CXX
#endif
