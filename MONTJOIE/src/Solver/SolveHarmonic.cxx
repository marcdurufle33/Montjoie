#ifndef MONTJOIE_FILE_SOLVE_HARMONIC_CXX

namespace Montjoie
{
  
  /***********************
   * All_LinearSolver *
   ***********************/

  
  //! Sets parameters to default values
  void All_LinearSolver::InitDefaultValues()
  {
    iterative_resolution = false;
    type_choice_solver_auto = AUTOMATIC_SOLVER;
    
    estimation_conditioning_number_asked = false;
    nb_iterations_powers_iterative = 10;
    
    mat_lu_real.SetPrintLevel(var_problem.print_level);
    mat_lu_cplx.SetPrintLevel(var_problem.print_level);
    mat_lu_real.SetThresholdMatrix(toDouble(var_computation.GetThresholdMatrix()));
    mat_lu_cplx.SetThresholdMatrix(toDouble(var_computation.GetThresholdMatrix()));
        
    static_condensation = false;
    compress_matrix = false;

    storage_matrix_asked = false;
    
    scaling_matrix = false;
    sym_matrix = false;
    complex_matrix = false;
    add_diagonal_part = false;
    use_old_pattern_Ah = false;
    new_pattern_Ah = true;
    nnz_old_Ah = 0;

    use_real_preconditioning = false;
    shallow_copy_real_precond = false;
    use_csr_matrices = false;
    facto_completed = false;

    use_symmetrization = true;
  }

  
  //! destructs object
  All_LinearSolver::~All_LinearSolver()
  {
    delete mat_iterative_real;
    delete mat_iterative_cplx;
    if (!shallow_copy_real_precond)
      delete precond_real;

    delete precond_cplx;
    delete solver_real;
    delete solver_cplx;
  }
  

  //! Selects the best solver
  template<class Dimension, class Prop>
  void All_LinearSolver
  ::SelectOptimalLinearSolver(int order, int N, Dimension dim, Prop property)
  {
    if (type_choice_solver_auto == USER_CHOICE)
      return;
    
    // best direct solver
    mat_lu_real.SelectOptimalDirectSolver(N, Real_wp(0), dim);
    mat_lu_cplx.SelectOptimalDirectSolver(N, Complex_wp(0), dim);
    
    // for direct solver
    int nb_max_direct, nb_max_no_condensation, order_min_condensation;
    mat_lu_real.GetMaximumSize(Real_wp(0), dim, nb_max_direct,
			       nb_max_no_condensation, order_min_condensation);
    
    // best iterative solver
    iterative_solver_real.SelectOptimalIterativeSolver(N, Real_wp(0), dim, property);
    iterative_solver_cplx.SelectOptimalIterativeSolver(N, Complex_wp(0), dim, property);
    
    // for iterative solver
    int nb_max_multigrid(0), nb_max_ilut(0), order_min_multigrid(1), nb_max_iterative(0);
    if (precond_real != NULL)
      precond_real->GetMaximumSize(Real_wp(0), dim, nb_max_ilut, nb_max_iterative,
				   nb_max_multigrid, order_min_multigrid);
 
    if (precond_cplx != NULL)
      precond_cplx->GetMaximumSize(Real_wp(0), dim, nb_max_ilut, nb_max_iterative,
				   nb_max_multigrid, order_min_multigrid);
   
    int type_preconditioner = PreconditioningEnum::IDENTITY;
    if (type_choice_solver_auto == AUTOMATIC_SOLVER)
      {
	if (N < nb_max_direct)
	  {
	    if ((order >= order_min_condensation)||(N < nb_max_no_condensation))
	      {
		// direct solver
		iterative_resolution = false;
		//if (N >= nb_max_no_condensation)
                //static_condensation = true;
	      }
	    else if (N < nb_max_ilut)
	      {
		// iterative solver with ilut preconditioner
		iterative_resolution = true;
		type_preconditioner = PreconditioningEnum::LOW_ORDER;
	      }
	    else if ((N < nb_max_multigrid)&&(order >= order_min_multigrid))
	      {
		// iterative solver with multigrid preconditioner
		iterative_resolution = true;
		type_preconditioner = PreconditioningEnum::MULTIGRID;
	      }
	    else if (N < nb_max_iterative)
	      {
		// iterative solver
		iterative_resolution = true;
		type_preconditioner = PreconditioningEnum::IDENTITY;
	      }
	    else
	      {
		cout<<"Very large sparse linear system, are you sure to solve it ? "<< endl;
		cout << "Type 1 for yes!" << endl;
		int reponse = 0; cin>>reponse;
		if (reponse != 1)
		  abort();
		
		iterative_resolution = true;
		type_preconditioner = PreconditioningEnum::IDENTITY;
	      }
	  }
	    
      }
    else if (type_choice_solver_auto == DIRECT_SOLVER)
      {
	iterative_resolution = false;
	//if (N >= nb_max_no_condensation)
        //static_condensation = true;
      }
    else if (type_choice_solver_auto == ITERATIVE_SOLVER)
      {
	// iterative solver
	iterative_resolution = true;
	type_preconditioner = PreconditioningEnum::IDENTITY;
      }
    else if (type_choice_solver_auto == MULTIGRID_SOLVER)
      {
	// iterative solver with multigrid preconditioner
	iterative_resolution = true;
	type_preconditioner = PreconditioningEnum::MULTIGRID;
      }
    else if (type_choice_solver_auto == ILUT_SOLVER)
      {
	// iterative solver with multigrid preconditioner
	iterative_resolution = true;
	type_preconditioner = PreconditioningEnum::LOW_ORDER;
      }
    
    precond_real->SetPreconditioningType(type_preconditioner);
    precond_cplx->SetPreconditioningType(type_preconditioner);
    precond_real->SetOptimalParameters(order, N, *this);
    precond_cplx->SetOptimalParameters(order, N, *this);
  }
  
  
  //! sets solver with a line of the datafile
  void All_LinearSolver
  ::SetInputData(const string& description_field, const Vector<string>& parameters)
  {
    precond_real->SetInputData(description_field, parameters);
    precond_cplx->SetInputData(description_field, parameters);

    iterative_solver_real.SetInputData(description_field, parameters);
    iterative_solver_cplx.SetInputData(description_field, parameters);
    
    if (!description_field.compare("TypeSolver"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "TypeSolver needs more parameters, for instance :" << endl;
	    cout << "TypeSolver = AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("AUTO"))
	  type_choice_solver_auto = AUTOMATIC_SOLVER;
	else if (!parameters(0).compare("DIRECT"))
	  type_choice_solver_auto = DIRECT_SOLVER;
	else if (!parameters(0).compare("ITERATIVE"))
	  type_choice_solver_auto = ITERATIVE_SOLVER;
	else if (!parameters(0).compare("MULTIGRID"))
	  type_choice_solver_auto = MULTIGRID_SOLVER;
	else if (!parameters(0).compare("ILUT"))
	  type_choice_solver_auto = ILUT_SOLVER;
      }    
    else if (!description_field.compare("TypeResolution")) 
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "TypeResolution needs more parameters, for instance :" << endl;
	    cout << "TypeResolution = MUMPS" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_choice_solver_auto = USER_CHOICE;
	iterative_resolution = mat_lu_real.SelectDirectSolver(parameters);
        iterative_resolution = mat_lu_cplx.SelectDirectSolver(parameters);
      }
    else if (!description_field.compare("RefineSolution"))
      {
        if (parameters(0) == "YES")
          {
            this->mat_lu_real.RefineSolution();
            this->mat_lu_cplx.RefineSolution();
          }
        else
          {
            this->mat_lu_real.DoNotRefineSolution();
            this->mat_lu_cplx.DoNotRefineSolution();
          }
      }
    else if (!description_field.compare("PivotThreshold"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "PivotThreshold needs more parameters, for instance :" << endl;
	    cout << "PivotThreshold = epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        mat_lu_real.SetPivotThreshold(to_num<double>(parameters(0)));
        mat_lu_cplx.SetPivotThreshold(to_num<double>(parameters(0)));
      }
    else if (!description_field.compare("RealPreconditioning"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "RealPreconditioning needs more parameters, for instance :" << endl;
	    cout << "RealPreconditioning = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (parameters(0) == "YES")
	  use_real_preconditioning = true;
	else
	  use_real_preconditioning = false;
      }
    else if (!description_field.compare("NbThreadsPerNode"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "NbThreadsPerNode needs more parameters, for instance :" << endl;
	    cout << "NbThreadsPerNode = n" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        mat_lu_real.SetNumberOfThreadPerNode(to_num<int>(parameters(0)));
        mat_lu_cplx.SetNumberOfThreadPerNode(to_num<int>(parameters(0)));
      }
    else if (!description_field.compare("PrintLevel")) 
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "PrintLevel needs more parameters, for instance :" << endl;
	    cout << "PrintLevel = level" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        int level = to_num<int>(parameters(0));
        SetPrintLevel(level);
      }
    else if (!description_field.compare("EstimationConditionNumber"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "EstimationConditionNumber needs more parameters, for instance :" << endl;
	    cout << "EstimationConditionNumber = YES n" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("YES"))
	  {
	    estimation_conditioning_number_asked = true;
	    nb_iterations_powers_iterative = to_num<int>(parameters(1));  
	  }  
      }
    else if (!description_field.compare("StaticCondensation"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "StaticCondensation needs more parameters, for instance :" << endl;
	    cout << "StaticCondensation = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	static_condensation = false;
	if (!parameters(0).compare("YES"))
	  {
	    static_condensation = true;
	    if (parameters.GetM() > 1)
	      {
		if (parameters(1) == "Uncompressed")
		  compress_matrix = false;
		else
		  compress_matrix = true;

		if (parameters(1) == "Store")
		  {
		    solver_real->SetStoreBlock();
		    solver_cplx->SetStoreBlock();
		  }
		else
		  {
		    solver_real->SetStoreBlock(false);
		    solver_cplx->SetStoreBlock(false);
		  }
	      }

	    if (parameters.GetM() > 2)
	      if (parameters(2) == "Optimal")
		{
		  solver_real->SetOptimalCondensation(true);
		  solver_cplx->SetOptimalCondensation(true);
		}
	  }
      }
    else if (!description_field.compare("MatrixCSR"))
      {
	if (parameters(0) == "YES")
	  use_csr_matrices = true;
	else
	  use_csr_matrices = false;
      }
    else if (!description_field.compare("ScalingMatrix"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "ScalingMatrix needs more parameters, for instance :" << endl;
	    cout << "ScalingMatrix = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          scaling_matrix = true;
	else
          scaling_matrix = false;
      }
    else if (!description_field.compare("StorageMatrix"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "StorageMatrix needs more parameters, for instance :" << endl;
	    cout << "StorageMatrix = YES nom_fichier" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	this->storage_matrix_asked = false;
	if (!parameters(0).compare("YES"))
	  {
	    this->storage_matrix_asked = true;
	    if (parameters(1).size() > 1)
	      this->file_name_matrix_stored = string(parameters(1));
	  }
      }    
    else if (!description_field.compare("MumpsMemoryCoefficient"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of All_LinearSolver" << endl;
	    cout << "MumpsMemoryCoefficient needs more parameters, for instance :" << endl;
	    cout << "MumpsMemoryCoefficient = coef_min coef_max increment" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// coefficients that the user can set for 
	// memory overcost that happens when calling Mumps
        double coef_init, coef_max, inc_coef;
        coef_init = to_num<double>(parameters(0));
        coef_max = to_num<double>(parameters(1));
        inc_coef = to_num<double>(parameters(2));
        mat_lu_real.SetCoefficientEstimationNeededMemory(coef_init);
        mat_lu_real.SetMaximumCoefficientEstimationNeededMemory(coef_max);
        mat_lu_real.SetIncreaseCoefficientEstimationNeededMemory(inc_coef);

        mat_lu_cplx.SetCoefficientEstimationNeededMemory(coef_init);
        mat_lu_cplx.SetMaximumCoefficientEstimationNeededMemory(coef_max);
        mat_lu_cplx.SetIncreaseCoefficientEstimationNeededMemory(inc_coef);
      }
  }

  
  //! parameters of preconditioning
  void All_LinearSolver
  ::SetInputPreconditioning(const string& description_field, const Vector<string>& parameters)
  {
    precond_cplx->SetInputPreconditioning(description_field, parameters);
  }


  //! enables static condensation
  void All_LinearSolver::EnableStaticCondensation(bool flag, bool store, bool compress)
  {
    if (flag)
      {
	static_condensation = true;
	if (store)
	  {
	    solver_real->SetStoreBlock();
	    solver_cplx->SetStoreBlock();
	  }
	else
	  {
	    solver_real->SetStoreBlock(false);
	    solver_cplx->SetStoreBlock(false);
	  }

	compress_matrix = compress;
      }
    else
      {
	static_condensation = false;
	compress_matrix = false;
      }
  }


  //! enables writing of finite element matrix
  void All_LinearSolver::EnableMatrixStorage(bool flag, const string& file_name)
  {
    storage_matrix_asked = flag;
    file_name_matrix_stored = file_name;
  }
  
  
  //! copies parameters contained in the linear solver given in argument
  void All_LinearSolver::CopyParameter(const All_LinearSolver& solver)
  {
    type_choice_solver_auto = solver.type_choice_solver_auto;
    iterative_resolution = solver.iterative_resolution;

    mat_lu_real.CopyParameter(solver.mat_lu_real);
    mat_lu_cplx.CopyParameter(solver.mat_lu_cplx);
    
    iterative_solver_real.CopyParameter(solver.iterative_solver_real);
    iterative_solver_cplx.CopyParameter(solver.iterative_solver_cplx);

    precond_real->CopyParameter(*solver.precond_real);
    precond_cplx->CopyParameter(*solver.precond_cplx);
    
    static_condensation = solver.static_condensation;
    compress_matrix = solver.compress_matrix;
    use_symmetrization = solver.use_symmetrization;
    
    solver_real->CopyParameter(*solver.solver_real);
    solver_cplx->CopyParameter(*solver.solver_cplx);

    scaling_matrix = solver.scaling_matrix;
    use_real_preconditioning = solver.use_real_preconditioning;
  }
  

  //! to handled overloaded virtual functions
  void All_LinearSolver::CopyParameter(const Preconditioner_Base<Complex_wp>&)
  {
    cout << "this method should not be called" << endl;
    abort();
  }
  
  
  //! adds memory used by solver
  void All_LinearSolver::GetMemoryUsed(map<string, size_t>& var, bool add) const
  {
    size_t size = 0;
    if (!iterative_resolution)
      {
        size = mat_lu_real.GetMemorySize() + mat_lu_cplx.GetMemorySize();
        if (add)
	  var["Factorisation"] += size;
	else
	  var["Factorisation"] = size;

	size = solver_real->GetMemorySize() + solver_cplx->GetMemorySize();
	if (add)
	  var["CondensationSolver"] += size;
	else
	  var["CondensationSolver"] = size;
      }
    else
      {
        if (complex_matrix)
          size = iterative_solver_cplx.GetNbIterativeVectors()
            *sizeof(Complex_wp)*var_problem.GetNbDof();
        else
          size = iterative_solver_real.GetNbIterativeVectors()
            *sizeof(Real_wp)*var_problem.GetNbDof();

	size += iterative_solver_real.GetMemorySize() + iterative_solver_cplx.GetMemorySize();
        if (add)
	  var["IterativeVectors"] += size;
	else
	  var["IterativeVectors"] = size;
        
        size = mat_iterative_real->GetMemorySize() + mat_iterative_cplx->GetMemorySize();
	size += diagonal_real_part_to_add.GetMemorySize();
	size += diagonal_complex_part_to_add.GetMemorySize();
        if (add)
	  var["IterativeMatrix"] += size;
	else
	  var["IterativeMatrix"] = size;
        
        size = precond_cplx->GetMemorySize();
	if (!shallow_copy_real_precond)
	  size += precond_real->GetMemorySize();
	
        if (add)
	  var["Preconditioning"] += size;
	else
	  var["Preconditioning"] = size;

	size = solver_real->GetMemorySize() + solver_cplx->GetMemorySize();
	if (add)
	  var["CondensationSolver"] += size;
	else
	  var["CondensationSolver"] = size;
      }
  }
  

  //! returns size of memory used by the object in bytes
  size_t All_LinearSolver::GetMemorySize() const
  {
    size_t taille = mat_iterative_real->GetMemorySize() + mat_iterative_cplx->GetMemorySize();
    taille += mat_lu_real.GetMemorySize() + mat_lu_cplx.GetMemorySize();
    taille += iterative_solver_real.GetMemorySize() + iterative_solver_cplx.GetMemorySize();
    taille += precond_real->GetMemorySize() + precond_cplx->GetMemorySize();
    taille += sizeof(*this);
    return taille;
  }
  

  //! direct factorisation and preconditioning are erased
  void All_LinearSolver::ClearFactorization()
  {
    // clearing LU factorisation (when direct solver is selected)
    mat_lu_real.Clear();
    mat_lu_cplx.Clear();
    
    // clearing preconditioning (when iterative solver is selected
    precond_real->Clear();
    precond_cplx->Clear();
  }
  
  
  //! releases memory used by the object
  void All_LinearSolver::Clear()
  {
    ClearFactorization();
    
    mat_iterative_real->Clear();
    mat_iterative_cplx->Clear();
    
    precond_real->Clear();
    precond_cplx->Clear();

    diagonal_real_part_to_add.Clear();
    diagonal_complex_part_to_add.Clear();
  }
  
  
  //! computation of finite element matrix and factorization
  /*!
    \param[in] nat_mat mass, damping and stiffness coefficient
    \param[in] precond if set to false, the preconditioning is not recomputed
    The finite element matrix is computed with coefficients of nat_mat.
    If a direct solver has been selected, the matrix will be factorized.
    If an iterative solver has been selected, the preconditioning is constructed
    if precond is true, otherwise the preconditioner previously computed will be used.
  */
  template<class T>
  void All_LinearSolver
  ::PerformFactorizationStep(const GlobalGenericMatrix<T>& nat_mat, bool precond)
  {
    facto_completed = true;
    if (this->static_condensation)
      var_computation.SetLeafStaticCondensation(true);

    use_symmetrization = var_boundary.GetInitialSymmetrization();
    
    T z;
    var_computation.SetSymmetrizationUse(use_symmetrization);
    GetCondensedSolver(z).SetSymmetrization(use_symmetrization);
    GetCondensedSolver(z).SetCompressionSystem(this->compress_matrix);
    
    complex_matrix = IsComplexNumber(z);
    sym_matrix = var_problem.IsSymmetricProblem();
    GetCondensedSolver(z).SetSymmetryElementaryMatrix(var_computation.IsSymmetricElementaryMatrix(nat_mat));
    GetNatureMatrix(z) = nat_mat;
 
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if (rank_proc == 0)
      if (var_problem.print_level >= 0)
        {
          cout << "The linear solver will be constructed " << endl;
          if (iterative_resolution)
            cout << "You have chosen an iterative solver, the preconditioning will be computed" << endl;
          else
            cout << "You have chosen a direct solver, the matrix will be computed and factorized" << endl;
        }
    
    if (iterative_resolution)
      {
        var_problem.var_chrono.SetMessage("PreconditioningConstruction",
                                          "compute the preconditioning");
        var_problem.var_chrono.Start("PreconditioningConstruction");
	// computation of the iterative matrix and the preconditioner
	ComputeMatrixAndPreconditioning(GetIterativeMatrix(z), GetPreconditioning(z),
					nat_mat, precond);
	
        var_problem.var_chrono.Stop("PreconditioningConstruction");
        
	// reducing space used by the matrix
        // lines commented because preconditioning may use stored matrices
        //var_problem.var_chrono.SetMessage("ConversionMatrix", "convert sparse matrices");
        //var_problem.var_chrono.Start("ConversionMatrix");
	//mat_iterative_real->CompressMatrix();	
	//mat_iterative_cplx->CompressMatrix();	
        //var_problem.var_chrono.Stop("ConversionMatrix");
      }
    else
      {
	Real_wp tmp = var_computation.GetThresholdMatrix();
        var_computation.SetThresholdMatrix(mat_lu_real.GetThresholdMatrix());    
	
        // computation and factorization of the direct matrix
	if (sym_matrix)
          {
	    ComputeMatrixAndFactorize(GetDirectSymmetricMatrix(z),
				      GetDirectFactorization(z), nat_mat);
	  }
        else
          {
	    ComputeMatrixAndFactorize(GetDirectUnsymmetricMatrix(z),
				      GetDirectFactorization(z), nat_mat);
	  }
	
	var_computation.SetThresholdMatrix(tmp);
      }

    
    if (this->static_condensation)
      var_computation.SetLeafStaticCondensation(false);

    var_computation.SetSymmetrizationUse(false);
  }


  //! solves A x = b with iterative solver in parallel
  template<class T>
  void All_LinearSolver::SolveIterativeSystem(Vector<T>& source_rhs, Vector<T>& x_sol)
  {
    T z;
    
#ifdef SELDON_WITH_MPI    
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    if (nb_proc > 1)
      {
	DistributedVector<T>* source, *solution;
	if (compress_matrix)
	  {
	    source = GetCondensedSolver(z).AllocateDistributedVector(source_rhs);
	    solution = GetCondensedSolver(z).AllocateDistributedVector(x_sol);
	  }
	else
	  {
	    source = var_comm.AllocateDistributedVector(source_rhs);
	    solution = var_comm.AllocateDistributedVector(x_sol);
	  }
	
	if (rank_proc != 0)
	  GetIterativeSolver(z).SetPrintLevel(0);
	
	GetIterativeSolver(z).Solve(GetIterativeMatrix(z), *solution, *source);
	    
	var_comm.NullifyDistributedVector(source);
	var_comm.NullifyDistributedVector(solution);
      }
    else
      {
#endif
	GetIterativeSolver(z).Solve(GetIterativeMatrix(z), x_sol, source_rhs);
	
#ifdef SELDON_WITH_MPI
      }
#endif
  }
  

  //! applies complex preconditioning
  void All_LinearSolver::Solve(const VirtualMatrix<Complex_wp>& A,
			       const Vector<Complex_wp>& r, Vector<Complex_wp>& z)
  {
    if (use_real_preconditioning)
      precond_real->Solve(SeldonNoTrans, *mat_iterative_real, r, z);
    else
      precond_cplx->Solve(A, r, z);
  }


  //! applies complex preconditioning
  void All_LinearSolver::TransSolve(const VirtualMatrix<Complex_wp>& A,
				    const Vector<Complex_wp>& r, Vector<Complex_wp>& z)
  {
    if (use_real_preconditioning)
      precond_real->Solve(SeldonTrans, *mat_iterative_real, r, z);
    else
      precond_cplx->TransSolve(A, r, z);
  }
  

  //! solves the linear system A x_sol = source_rhs
  /*!
    \param[inout] x_sol on input the right hand side, on output the solution    
  */
  template<class T>
  void All_LinearSolver
  ::ComputeSolution(Vector<T>& x_sol, bool dirichlet_inhg, bool assemble)
  {
    ComputeSolution(x_sol, GetNatureMatrix(T(0)), dirichlet_inhg, assemble);
  }

  
  //! solves the linear system A x_sol = source_rhs
  /*!
    \param[inout] x_sol on input the right hand side, on output the solution    
  */
  template<class T>
  void All_LinearSolver
  ::ComputeSolution(Vector<T>& x_sol,
		    const GlobalGenericMatrix<T>& nature_matrix,
		    bool dirichlet_inhg, bool assemble)
  {
    if (this->use_symmetrization)
      var_computation.ModifySourceSymmetry(x_sol);

    Vector<T> source_dirichlet(x_sol.GetM());
    source_dirichlet.Zero();
    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      {
        source_dirichlet(var_boundary.GetDirichletDofNumber(i))
          = x_sol(var_boundary.GetDirichletDofNumber(i));        
      }
    
    T z;
    // modification of the right hand side 
    if (static_condensation)
      {
	if (assemble)
	  {
	    var_comm.AssembleDirichlet(source_dirichlet);
	    GetCondensedSolver(z)
              .ModifyRhsStaticCondensation(SeldonNoTrans, x_sol, nature_matrix);
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
      }

    if (dirichlet_inhg)
      {
	if (!assemble)
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
	
	var_boundary.ApplyDirichletCondition(SeldonNoTrans, GetIterativeMatrix(z), x_sol, 0);
      }

    if (assemble)
      var_comm.AddDomains(x_sol);

    // Right hand side is modified for inhomogeneous Dirichlet conditions    
    // Then we solve the linear system A x_sol = source_rhs
    if (iterative_resolution)
      {	
	if (compress_matrix)	  
	  {
	    Vector<T> xtmp;
	    GetCondensedSolver(z).ExtractVector(x_sol, xtmp);

	    Vector<T> source_rhs(xtmp);
	    xtmp.Zero();
	    
	    SolveIterativeSystem(source_rhs, xtmp);
	    
	    GetCondensedSolver(z).ScatterVector(xtmp, x_sol);
	  }
	else
	  {
	    Vector<T> source_rhs(x_sol);
	    x_sol.Zero();

	    SolveIterativeSystem(source_rhs, x_sol);
	  }
      }
    else
      {        
	if (compress_matrix)	  
	  {
	    Vector<T> xtmp;
	    GetCondensedSolver(z).ExtractVector(x_sol, xtmp);
	    
	    GetDirectFactorization(z).Solve(xtmp);

	    GetCondensedSolver(z).ScatterVector(xtmp, x_sol);
	  }
	else
	  GetDirectFactorization(z).Solve(x_sol);
      }

    // inhomogeneous Dirichlet is set
    if (dirichlet_inhg)
      GetIterativeMatrix(z).ImposeDirichletCondition(SeldonNoTrans, x_sol, 0);

    // in the case of a static condensation, internal nodes are recomputed
    if (static_condensation)
      {
	GetCondensedSolver(z).RecomposeSolution(SeldonNoTrans, x_sol, nature_matrix);

	// inhomogeneous Dirichlet is repeated
	for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	  {
	    x_sol(var_boundary.GetDirichletDofNumber(i))
	      = source_dirichlet(var_boundary.GetDirichletDofNumber(i));
	  }
      }
  }


  //! solves the linear system A x_sol = source_rhs
  /*!
    \param[inout] x_sol on input the right hand side, on output the solution    
  */
  template<class T>
  void All_LinearSolver
  ::ComputeSolution(const SeldonTranspose& trans, Vector<T>& x_sol,
		    const GlobalGenericMatrix<T>& nature_matrix,
		    bool dirichlet_inhg, bool assemble)
  {
    if (this->use_symmetrization)
      var_computation.ModifySourceSymmetry(x_sol);

    T z;
    // modification of the right hand side 
    if (static_condensation)
      {
	if (assemble)
	  {
	    //var_comm.AssembleDirichlet(source_dirichlet);
	    GetCondensedSolver(z)
              .ModifyRhsStaticCondensation(trans, x_sol, nature_matrix);
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
      }

    if (dirichlet_inhg)
      {
	if (!assemble)
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
	
        var_boundary.ApplyDirichletCondition(trans, GetIterativeMatrix(z), x_sol, 0);
      }
    
    if (assemble)
      var_comm.AddDomains(x_sol);

    // Right hand side is modified for inhomogeneous Dirichlet conditions    
    // Then we solve the linear system A x_sol = source_rhs
    if (iterative_resolution)
      {	
        cout << "Not implemented" << endl;
        abort();
      }
    else
      {        
	if (compress_matrix)	  
	  {
	    Vector<T> xtmp;
	    GetCondensedSolver(z).ExtractVector(x_sol, xtmp);
	    
	    GetDirectFactorization(z).Solve(trans, xtmp);

	    GetCondensedSolver(z).ScatterVector(xtmp, x_sol);
	  }
	else
	  GetDirectFactorization(z).Solve(trans, x_sol);
      }

    // inhomogeneous Dirichlet is set
    if (dirichlet_inhg)
      GetIterativeMatrix(z).ImposeDirichletCondition(trans, x_sol, 0);
    
    // in the case of a static condensation, internal nodes are recomputed
    if (static_condensation)
      {
        Vector<T> source_dirichlet(x_sol.GetM());
        source_dirichlet.Zero();
        for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
          {
            source_dirichlet(var_boundary.GetDirichletDofNumber(i))
              = x_sol(var_boundary.GetDirichletDofNumber(i));        
          }
        
	GetCondensedSolver(z).RecomposeSolution(trans, x_sol, nature_matrix);
        
	// inhomogeneous Dirichlet is repeated
	for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	  {
	    x_sol(var_boundary.GetDirichletDofNumber(i))
	      = source_dirichlet(var_boundary.GetDirichletDofNumber(i));
	  }
      }
  }


  //! solves multiple right hand sides
  /*!
    \param[inout] x_sol on input right hand sides, on output solutions of the linear system
  */
  template<class T>
  void All_LinearSolver
  ::ComputeSolution(Matrix<T, General, ColMajor>& M_sol,
		    const GlobalGenericMatrix<T>& nat_mat, bool dirichlet_inhg, bool assemble)
  {
    ComputeSolution(SeldonNoTrans, M_sol, nat_mat, dirichlet_inhg, assemble);
    
  }


  //! solves multiple right hand sides
  /*!
    \param[inout] x_sol on input right hand sides, on output solutions of the linear system
  */
  template<class T>
  void All_LinearSolver
  ::ComputeSolution(const SeldonTranspose& trans, Matrix<T, General, ColMajor>& M_sol,
		    const GlobalGenericMatrix<T>& nat_mat, bool dirichlet_inhg, bool assemble)
  {
    typedef T Complexe; T z;
    // modification of the right hand side 
    GetIterativeMatrix(z).SetNbDirichletCondition(M_sol.GetN());

    if (this->use_symmetrization)
      var_computation.ModifySourceSymmetry(M_sol);
      
    // Right hand side is modified for inhomogeneous Dirichlet conditions
    Vector<Complexe> x_sol, source_rhs, xtmp;
    x_sol.Reallocate(M_sol.GetM());
    if (dirichlet_inhg)
      {
	if (!assemble)
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
	
	for (int k = 0; k < M_sol.GetN(); k++)
	  {
	    x_sol.SetData(M_sol.GetM(), &M_sol(0, k));
	    var_boundary.ApplyDirichletCondition(trans, GetIterativeMatrix(z), x_sol, k);
	    x_sol.Nullify();
	  }
      }

    // assembling values between processors if needed
    if (assemble)
      for (int k = 0; k < M_sol.GetN(); k++)
	{
	  x_sol.SetData(M_sol.GetM(), &M_sol(0, k));
	  var_comm.AddDomains(x_sol);
	  x_sol.Nullify();
	}

    // Then we solve the linear system A x_sol = source_rhs
    if (iterative_resolution)
      {
	for (int k = 0; k < M_sol.GetN(); k++)
	  {
	    GetCol(M_sol, k, x_sol);

	    if (compress_matrix)	  
	      {
		GetCondensedSolver(z).ExtractVector(x_sol, xtmp);

		source_rhs = xtmp;
		xtmp.Zero();
		
		SolveIterativeSystem(source_rhs, xtmp);
		
		GetCondensedSolver(z).ScatterVector(xtmp, x_sol);
	      }
	    else
	      {
		source_rhs = x_sol;
		x_sol.Zero();
		
		SolveIterativeSystem(source_rhs, x_sol);
	      }
	    
	    SetCol(x_sol, k, M_sol);
	  }
      }
    else
      {
	if (compress_matrix)
	  {
	    Matrix<T, General, ColMajor> M_tmp;
	    for (int k = 0; k < M_sol.GetN(); k++)
	      {
		GetCol(M_sol, k, x_sol);
		
		GetCondensedSolver(z).ExtractVector(x_sol, xtmp);
		if (k == 0)
		  M_tmp.Reallocate(xtmp.GetM(), M_sol.GetN());

		SetCol(xtmp, k, M_tmp);
	      }

	    GetDirectFactorization(z).Solve(trans, M_tmp);

	    for (int k = 0; k < M_sol.GetN(); k++)
	      {
		GetCol(M_tmp, k, xtmp);
		GetCondensedSolver(z).ScatterVector(xtmp, x_sol);
		SetCol(x_sol, k, M_sol);
	      }	    
	  }
	else
	  GetDirectFactorization(z).Solve(trans, M_sol);
      }
    
    // in the case of a static condensation, internal nodes are recomputed
    if (static_condensation)
      {
	cout << "Not implemented" << endl;
	abort();
      }

    // inhomogeneous Dirichlet is set
    if (dirichlet_inhg)
      for (int k = 0; k < M_sol.GetN(); k++)
	{
	  x_sol.SetData(M_sol.GetM(), &M_sol(0, k));
	  GetIterativeMatrix(z).ImposeDirichletCondition(trans, x_sol, k);
	  x_sol.Nullify();
	}
  }

  
  //! computation of the direct matrix and factorization
  /*!
    \param[out] mat_direct finite element matrix
    \param[out] mat_lu LU matrix
    \param[in] nature_matrix mass and stiffness coefficients
  */
  template<class T, class Prop, class Storage>
  void All_LinearSolver
  ::ComputeMatrixAndFactorize(DistributedMatrix<T, Prop, Storage>& Ah,
			      All_MatrixLU<T>& matLU,
			      const GlobalGenericMatrix<T>& nature_matrix)
  {
    if (this->static_condensation)
      var_computation.SetLeafStaticCondensation(true);

    T z;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    // the finite element matrix is computed
    Ah.Zero();
    
    var_problem.var_chrono.SetMessage("ComputationMatrix", "compute finite element matrix");
    var_problem.var_chrono.Start("ComputationMatrix");

    var_computation.AddMatrixWithBC(GetIterativeMatrix(T(0)), Ah, nature_matrix, 0, 0,
				    &GetCondensedSolver(z));
    
    if (add_diagonal_part)
      {
	T value;
	if (complex_matrix)
	  for (int i = 0; i < Ah.GetM(); i++)
	    {
	      to_complex(diagonal_complex_part_to_add(i), value);
	      Ah.AddInteraction(i, i, value);
	    }
	else
	  for (int i = 0; i < Ah.GetM(); i++)
	    {
	      to_complex(diagonal_real_part_to_add(i), value);
	      Ah.AddInteraction(i, i, value);
	    }
      }
    
    if (this->static_condensation)
      GetCondensedSolver(z).Init(nature_matrix, this->compress_matrix);

    if (this->compress_matrix)
      GetCondensedSolver(z).CompressMatrixGen(Ah);

    var_problem.var_chrono.Stop("ComputationMatrix");
    
    if (var_problem.print_level >= 2)
      {
        if (rank_proc == 0)
	  cout << rank_proc << " Matrix computed successfully " << endl;
	  
        var_problem.var_chrono.DisplayTime("ComputationMatrix");
      }
    
    // if the user wants the matrix to be stored in a file
    int mode_number = var_boundary.GetCurrentModeNumber();
    if (this->storage_matrix_asked)
      {
        if (var_comm.GetNbProcPerMode() == 1)
          {
            if (GetExtension(this->file_name_matrix_stored) == "mtx")
              WriteMatrixMarket(Ah, GetBaseString(this->file_name_matrix_stored)
                                +to_str(mode_number)+".mtx");
            else
              Ah.WriteText(GetBaseString(this->file_name_matrix_stored)
                           +to_str(mode_number)+".dat");
          }
        else
          Ah.WriteText(GetBaseString(this->file_name_matrix_stored)
                       +to_str(mode_number)+".dat");
      }
    
    var_problem.var_chrono.SetMessage("FactorizationMatrix", "factorize finite element matrix");
    var_problem.var_chrono.Start("FactorizationMatrix");
    
    if (estimation_conditioning_number_asked)
      {
	// extremals eigenvalues evaluated with iterative power method
 	T lambda_max,lambda_min, one, zero;
	SetComplexOne(one); SetComplexZero(zero);
	
	// random vectors
	int nodl = Ah.GetM();
	Vector<T> X(nodl), Y(nodl);
	X.FillRand();
	
	// iteration of X = A*X to get the biggest eigenvalue
	for (int i = 0; i < nb_iterations_powers_iterative; i++)
	  {
	    Mlt(1.0/Norm2(X), X);
	    Copy(X, Y);
	    MltAdd(one, SeldonNoTrans,
		   Ah, Y, zero, X, true);
	  }
	
	lambda_max = Norm2(X)/Norm2(Y);
	
	// the matrix is factorized
	matLU.Factorize(Ah, false, scaling_matrix);
	
	// iteration of X = A^{-1} X to get the smallest eigenvalue
	Y.FillRand();
	for (int i = 0; i < nb_iterations_powers_iterative; i++)
	  {
	    Mlt(1.0/Norm2(Y), Y);
	    Copy(Y, X);
	    matLU.Solve(X);
	  }
	
	lambda_min = Norm2(Y) / Norm2(X);
	
	if ((rank_proc == 0)||(var_problem.print_level >= 10))
	  cout << "Smallest eigenvalue " << lambda_min << "\n Biggest eigenvalue " << lambda_max
	       << "\n CONDITION NUMBER " << abs(lambda_max/lambda_min) << "     Estimated " << endl;
      }
    else
      {
	int nnz_new_Ah = Ah.GetNonZeros();
	if (use_old_pattern_Ah)
	  {
	    DistributedMatrix<T, Prop, typename GenericStorage<Storage>::StorageCSR> Ah_CSR;
	    Copy(Ah, Ah_CSR); Ah.Clear();
	    if (new_pattern_Ah)
	      matLU.PerformAnalysis(Ah_CSR);
	    else
	      {
		if (nnz_new_Ah != nnz_old_Ah)
		  {
		    cout << "The profile of Ah changed" << endl;
		    DISP(nnz_new_Ah); DISP(nnz_old_Ah);
		    abort();
		  }
	      }
	    
	    matLU.PerformFactorization(Ah_CSR, scaling_matrix);	
	    new_pattern_Ah = false;
	  }
	else
	  matLU.Factorize(Ah, false, scaling_matrix);

	nnz_old_Ah = nnz_new_Ah;
      }

    var_problem.var_chrono.Stop("FactorizationMatrix");
    
    // printing error messages
    if (rank_proc == 0)
      {
        int ierr = 0;
        int type = matLU.GetInfoFactorization(ierr);
        PrintFactorizationFailed(type, ierr);
      }

    if (this->static_condensation)
      var_computation.SetLeafStaticCondensation(false);

  } // end method ComputeMatrixAndFactorize
  

  //! computation of iterative matrix and preconditioning
  /*!
    \param[out] mat_iterative iterative matrix computed
    \param[out] prec preconditioning constructed
    \param[in] precond if true the preconditioning is constructed
  */
  template<class T>
  void All_LinearSolver
  ::ComputeMatrixAndPreconditioning(FemMatrixFreeClass_Base<T>& Ah,
				    All_Preconditioner_Base<T>& prec,
				    const GlobalGenericMatrix<T>& nat_mat, bool precond)
  {
    Ah.SetSymmetrizationUse(this->use_symmetrization);
    
    if (this->static_condensation)
      var_computation.SetLeafStaticCondensation(true);

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    // computation of finite element matrix
    Ah.Clear(); T z;
    var_computation.AddMatrixWithBC(Ah, nat_mat, 0, 0, &GetCondensedSolver(z));
    
    if (GetIterativeSolver(z).GetSolverType() == iterative_solver_real.JACOBI)
      {
        Vector<T> diag;
        var_computation.ComputeDiagonalMatrix(diag, Ah, nat_mat);
        for (int i = 0; i < diag.GetM(); i++)
	  diag(i) = 1.0/diag(i);
        
        GetIterativeSolver(z).SetInvDiagonalJacobi(diag);
      }
    
    if (this->use_csr_matrices)
      Ah.CompressMatrix();
    
    if (this->static_condensation)
      {
	if (this->compress_matrix)
	  Ah.SetCondensedSolver(&GetCondensedSolver(z));
	else
	  Ah.SetCondensedSolver(NULL);	
      }
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 2)
	cout << rank_proc << " Matrix successfully computed" << endl;

    int nodl = var_problem.GetNbDof();
    if (scaling_matrix)
      {
	VectReal_wp diagonal_scale_left, diagonal_scale_right;
	// symmetric scaling
	diagonal_scale_left.Reallocate(nodl);
	diagonal_scale_right.Reallocate(nodl);
	diagonal_scale_left.Fill(0);
	diagonal_scale_right.Fill(0);
	
	GetRowColSum(diagonal_scale_left, diagonal_scale_right, Ah);
	
	for (int i = 0; i < nodl; i++)
	  if (abs(diagonal_scale_left(i)) > 1e-20)
	    diagonal_scale_left(i) = 1.0/diagonal_scale_left(i);
	  else
	    diagonal_scale_left(i) = 1.0;
	
	for (int i = 0; i < nodl; i++)
	  if (abs(diagonal_scale_right(i)) > 1e-20)
	    diagonal_scale_right(i) = 1.0/diagonal_scale_right(i);
	  else
	    diagonal_scale_right(i) = 1.0;
	
	GetIterativeSolver(z).SetScaling(diagonal_scale_left, diagonal_scale_right);
	ScaleMatrix(Ah, diagonal_scale_left, diagonal_scale_right);
      }

    if (this->static_condensation)
      GetCondensedSolver(z).Init(nat_mat, this->compress_matrix);

    if (this->compress_matrix)
      GetCondensedSolver(z).CompressMatrixGen(Ah);
    
    
    // the preconditioning is constructed
    if (precond)
      {
	bool use_real_prec = use_real_preconditioning;
	if (!IsComplexNumber(z))
	  use_real_prec = false;
	
	if (use_real_prec)
	  {
	    GlobalGenericMatrix<Real_wp> nat_mat_real;
	    nat_mat_real.SetCoefMass(realpart(nat_mat.GetCoefMass()));
	    nat_mat_real.SetCoefDamping(realpart(nat_mat.GetCoefDamping()));
	    nat_mat_real.SetCoefStiffness(realpart(nat_mat.GetCoefStiffness()));

	    CopyReal(Ah, *mat_iterative_real);
	    if (this->compress_matrix)
	      solver_real->CopyPtr(GetCondensedSolver(z));

	    if (this->static_condensation)
	      {
		if (this->compress_matrix)
		  mat_iterative_real->SetCondensedSolver(solver_real);
		else
		  mat_iterative_real->SetCondensedSolver(NULL);	
	      }
	    
	    precond_real->ConstructPreconditioner(var_computation, *this,
						  *mat_iterative_real, nat_mat_real, *solver_real);
	  }
	else
	  {
	    prec.ConstructPreconditioner(var_computation, *this,
					 Ah, nat_mat, GetCondensedSolver(z));
	  }	  
      }
    
    // if the storage of the matrix is asked
    // Lots of product matrix vectors are realized (same number of
    // MV products than the size of the matrix)
    // very slow procedure, so store the matrix only for small cases
    if (this->storage_matrix_asked)
      {
	int mode_number = var_boundary.GetCurrentModeNumber();
	T one, zero;
	SetComplexZero(zero);
	SetComplexOne(one);

        if (var_comm.GetNbProcPerMode() == 1)
          {
            Matrix<T, General, ArrayRowSparse> Ah_full(nodl, nodl);
            Vector<T> Ones(nodl), Ah_Ones(nodl);
            for (int i = 0; i < nodl; i++)
              {
                Ones.Fill(0); Ones(i) = one;
                cout << "colonne " << i << endl;
                Ah.MltAddVector(one, Ones, zero, Ah_Ones);
                //Ah.MltVector(SeldonTrans, Ones, Ah_Ones);
                for (int j = 0; j < nodl; j++)
                  if (Ah_Ones(j) != zero)
                    Ah_full.AddInteraction(j, i, Ah_Ones(j));
              }
            
            Ah_full.WriteText(string("Ah_full")+to_str(mode_number)+".dat");
            
            if (precond)
              {
                Ah_full.Clear(); Ah_full.Reallocate(nodl, nodl);
                for (int i = 0; i < nodl; i++)
                  {
                    Ones.Fill(0); Ones(i) = one;
                    cout << "colonne precond " << i << endl;
                    prec.Solve(Ah, Ones, Ah_Ones);
                    for (int j = 0; j < nodl; j++)
                      if (Ah_Ones(j) != zero)
                        Ah_full.AddInteraction(j, i, Ah_Ones(j));
                  }
                
                Ah_full.WriteText(string("Precond_full")+to_str(mode_number)+".dat");
              }
          }
        else
          {
            Vector<T> Ones(nodl), Ah_Ones(nodl);
            Ones.Fill(0); Ah_Ones.Fill(0);
            string file_name = string("Ah_full")+to_str(mode_number)+"P"
              + to_str(rank_proc) + ".dat";
            
            ofstream file_out(file_name.data());
            file_out.precision(15);
            for (int i = 0; i < var_comm.GetNbGlobalDof(); i++)
              {
                for (int j = 0; j < nodl; j++)
                  {
                    if (var_comm.GetGlobalDofNumber(j) == i)
                      Ones(j) = one;
                    else
                      Ones(j) = zero;
                  }
                
                cout << "colonne " << i << endl;
                Ah.MltAddVector(one, Ones, zero, Ah_Ones);
                for (int j = 0; j < nodl; j++)
                  if (Ah_Ones(j) != zero)
                    file_out << var_comm.GetGlobalDofNumber(j)+1 
                             << " " << i+1 << " " << Ah_Ones(j) << '\n';
	      }
            
            file_out.close();
	  }
      }

    if (this->static_condensation)
      var_computation.SetLeafStaticCondensation(false);
    
  } // end ComputeMatrixAndPreconditioning


  /********************************
   * CondensationBlockSolver_Base *
   ********************************/


  //! returns the memory used by the object in bytes
  template<class T, class Prop, class Storage>
  size_t CondensedBlockSchur<T, Prop, Storage>::GetMemorySize() const
  {
    size_t taille = inv_a22.GetMemorySize() + num_ddl.GetMemorySize();
    taille += a12.GetMemorySize() + a21.GetMemorySize();
    return taille;
  }
  
  
  //! inits default values
  template<class T>
  void CondensationBlockSolver_Fem<T>::InitDefaultValues()
  {
    size_compressed_system = 0;
    store_block = false;
    nb_scalar_dof = 0;
    nb_global_dof = 0;
    optim_condensation = false;
    use_symmetrization = false;
    compress_system = true;
    symmetric_elem_matrix = false;
    
    IndexCondensedRows = &IndexCondensedRows_store;
    GlobalCondensedRows = &GlobalCondensedRows_store;
    OverlappedRows = &OverlappedRows_store;
    OverlappedProcs = &OverlappedProcs_store;
    SharingProcs = &SharingProcs_store;
    SharingRows = &SharingRows_store;
  }
  

  //! if store is true, the blocks issued from Schur complement are stored
  template<class T>
  void CondensationBlockSolver_Fem<T>::SetStoreBlock(bool store)
  {
    store_block = store;
  }
  

  //! copies parameters of another solver
  template<class T>
  void CondensationBlockSolver_Fem<T>::CopyParameter(const CondensationBlockSolver_Fem<T>& solver)
  {
    store_block = solver.store_block;
    optim_condensation = solver.optim_condensation;
  }

  
  //! the matrix A is compressed (internal dofs are removed)
  template<class T> template<class Prop, class Storage>
  void CondensationBlockSolver_Fem<T>::CompressMatrixGen(DistributedMatrix<T, Prop, Storage>& Ah)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#endif

    CompressMatrix(Ah, *IndexCondensedRows);
    size_compressed_system = Ah.GetM();
    
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      Ah.Init(nb_global_dof, GlobalCondensedRows,
	      OverlappedRows, OverlappedProcs, nb_scalar_dof, Ah.GetNbScalarUnknowns(),
	      SharingProcs, SharingRows, Ah.GetCommunicator());
#endif
    
  }


  //! the matrix A is compressed
  template<class T>
  void CondensationBlockSolver_Fem<T>::CompressMatrixGen(FemMatrixFreeClass_Base<T>& A)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#endif

    if (!compress_system)
      return;
    
    if (A.IsSymmetric())
      {
	if (A.mat_iterative_sym.GetM() > 0)
	  {
	    typename GenericStorage<T>::SparseSymMatrix& Ah = A.mat_iterative_sym;

	    CompressMatrix(Ah, *IndexCondensedRows);
	    size_compressed_system = Ah.GetM();
	    
#ifdef SELDON_WITH_MPI
	    if (nb_proc > 1)
	      Ah.Init(nb_global_dof, GlobalCondensedRows,
		      OverlappedRows, OverlappedProcs, nb_scalar_dof, Ah.GetNbScalarUnknowns(),
		      SharingProcs, SharingRows, Ah.GetCommunicator());
#endif
	  }

	if (A.matCSR_iterative_sym.GetM() > 0)
	  {
	    typename GenericStorage<T>::SparseSymMatrixCSR& Ah = A.matCSR_iterative_sym;

	    CompressMatrix(Ah, *IndexCondensedRows);
	    size_compressed_system = Ah.GetM();
	    
#ifdef SELDON_WITH_MPI
	    if (nb_proc > 1)
	      Ah.Init(nb_global_dof, GlobalCondensedRows,
		      OverlappedRows, OverlappedProcs, nb_scalar_dof, Ah.GetNbScalarUnknowns(),
		      SharingProcs, SharingRows, Ah.GetCommunicator());
#endif
	  }
      }
    else
      {
	if (A.mat_iterative_unsym.GetM() > 0)
	  {
	    typename GenericStorage<T>::SparseMatrix& Ah = A.mat_iterative_unsym;

	    CompressMatrix(Ah, *IndexCondensedRows);
	    size_compressed_system = Ah.GetM();
	    
#ifdef SELDON_WITH_MPI
	    if (nb_proc > 1)
	      Ah.Init(nb_global_dof, GlobalCondensedRows,
		      OverlappedRows, OverlappedProcs, nb_scalar_dof, Ah.GetNbScalarUnknowns(),
		      SharingProcs, SharingRows, Ah.GetCommunicator());
#endif
	  }

	if (A.matCSR_iterative_unsym.GetM() > 0)
	  {
	    typename GenericStorage<T>::SparseMatrixCSR& Ah = A.matCSR_iterative_unsym;

	    CompressMatrix(Ah, *IndexCondensedRows);
	    size_compressed_system = Ah.GetM();
	    
#ifdef SELDON_WITH_MPI
	    if (nb_proc > 1)
	      Ah.Init(nb_global_dof, GlobalCondensedRows,
		      OverlappedRows, OverlappedProcs, nb_scalar_dof, Ah.GetNbScalarUnknowns(),
		      SharingProcs, SharingRows, Ah.GetCommunicator());
#endif
	  }
      }
    
    // changes the size set in A
    A.Reallocate(size_compressed_system, size_compressed_system);
  }


  //! compress the vector x_sol : xtmp = x_sol(dof_compressed)
  template<class T> template<class T0>
  void CondensationBlockSolver_Fem<T>::ExtractVector(const Vector<T0>& x_sol,
						     Vector<T0>& xtmp) const
  {
    xtmp.Reallocate(size_compressed_system);
    for (int i = 0; i < IndexCondensedRows->GetM(); i++)
      if ((*IndexCondensedRows)(i) >= 0)
	xtmp((*IndexCondensedRows)(i)) = x_sol(i);
  }


  //! uncompress the vector xtmp : x_sol(dof_compressed) = xtmp
  template<class T> template<class T0>
  void CondensationBlockSolver_Fem<T>::ScatterVector(const Vector<T0>& xtmp, Vector<T0>& x_sol) const
  {
    for (int i = 0; i < IndexCondensedRows->GetM(); i++)
      if ((*IndexCondensedRows)(i) >= 0)
	x_sol(i) = xtmp((*IndexCondensedRows)(i));
  }


  template<class T> template<class T0>
  void CondensationBlockSolver_Fem<T>::AddDomains(Vector<T0>& X, int nb_u) const
  {
#ifdef SELDON_WITH_MPI
    // for DG formulation, no need to assemble solution
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;
    
    if (nb_u == -1)
      {
	if (var_problem.FormulationDG() == ElementReference_Base::HDG)
	  nb_u = var_problem.nb_unknowns_hdg;
	else
	  nb_u = var_problem.nb_unknowns;
      }
    
    if (var_comm.GetNbProcPerMode() <= 1)
      return;

    AssembleVector(X, MPI_SUM, *this->SharingProcs, *this->SharingRows,
		   var_comm.comm_group_mode, nb_scalar_dof, nb_u, 23);
#endif

  }
  

#ifdef SELDON_WITH_MPI
  template<class T> DistributedVector<T>*
  CondensationBlockSolver_Fem<T>::AllocateDistributedVector(Vector<T>& x)
  {
    DistributedVector<T>* d;
    d = new DistributedVector<T>(*this->OverlappedRows, var_comm.comm_group_mode);  
    d->SetData(x.GetM(), x.GetData());
    return d;
  }
#endif

  
  //! computes the Schur complement from elementary matrix mat_interac
  template<class T> template<class Prop0, class Storage0, class Prop, class Storage>
  void CondensationBlockSolver_Fem<T>
  ::GetSchurMatrixGen(int i, int nb_dof_loc, Matrix<T, Prop0, Storage0>& mat_interac,
		      Matrix<T>& a12, Matrix<T>& a21, Matrix<T, Prop, Storage>& a11,
		      Matrix<T, Prop, Storage>& a22, Vector<int>& num_ddl_final)
  {
    Vector<int> intern_node;
    
    // we compute schur complement a11 - a12 a22^{-1} a21
    int nb_dof_edges, nb_dof_int;
    var_computation.GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges, nb_dof_int, intern_node);

    if (nb_dof_int == 0)
      {
	a22.Clear(); a21.Clear(); a11.Clear(); a12.Clear();
	return;
      }

    a22.Reallocate(nb_dof_int, nb_dof_int);
    a21.Reallocate(nb_dof_int, nb_dof_edges);
    a12.Reallocate(nb_dof_edges, nb_dof_int);
    a11.Reallocate(nb_dof_edges, nb_dof_edges);
    a22.Zero(); a21.Zero(); a12.Zero(); a11.Zero();

    // setting zero for dofs not present
    T zero; SetComplexZero(zero);
    T one; SetComplexOne(one);
    for (int i = 0; i < nb_dof_loc; i++)
      if (num_ddl_final(i) < 0)
	{
	  for (int j = 0; j < nb_dof_loc; j++)
	    mat_interac(i, j) = zero;

	  mat_interac(i, i) = one;
	}
    
    // matrices for Schur complement (static condensation)
    Matrix<T> a_tmp;
    a_tmp.Reallocate(nb_dof_int, nb_dof_edges);
    
    for (int j = 0; j < nb_dof_loc; j++)
      if (intern_node(j) < 0)
	{
	  int num_j = -intern_node(j)-1;
	  for (int k = 0; k < nb_dof_loc; k++)
	    {
	      if (intern_node(k) < 0)
		a22(num_j, -intern_node(k)-1) = mat_interac(j, k);
	      else
		{
		  a21(num_j, intern_node(k)) = mat_interac(j, k);
		  a12(intern_node(k), num_j) = mat_interac(k, j);
		}
	    }
	}

    // we compute the schur complement
    Seldon::GetInverse(a22);
    Seldon::Mlt(a22, a21, a_tmp);
    Seldon::Mlt(a12, a_tmp, a11);
    for (int j = 0; j < nb_dof_loc; j++)
      if (intern_node(j) >= 0)
	{
	  int kmin = 0;
	  if (IsSymmetricMatrix(mat_interac))
	    kmin = j;
	  
	  for (int k = kmin; k < nb_dof_loc; k++)
	    if (intern_node(k) >= 0)
	      mat_interac(j, k) -= a11(intern_node(j), intern_node(k));
	}

    a11.Zero();
    for (int j = 0; j < nb_dof_loc; j++)
      for (int k = 0; k < nb_dof_loc; k++)
	if ((intern_node(j) >= 0) && (intern_node(k) >= 0))
	  a11(intern_node(j), intern_node(k)) = mat_interac(j, k);
    
    // for internal dofs, we set the identity matrix
    if (!compress_system)
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    int kmin = 0;
	    if (IsSymmetricMatrix(mat_interac))
	      kmin = j;
	    
	    for (int k = kmin; k < nb_dof_loc; k++)
	      {
		if ((intern_node(j) < 0) || (intern_node(k) < 0))
		  {
		    mat_interac(j, k) = 0;
		    if (j == k)
		      mat_interac(j, j) = 1;
		  }
	      }
	  }
      }
    else
      {
	mat_interac.Reallocate(nb_dof_edges, nb_dof_edges);
	for (int j = 0; j < nb_dof_edges; j++)
	  for (int k = 0; k < nb_dof_edges; k++)
	    mat_interac(j, k) = a11(j, k);
	
	Vector<int> num_ddl(num_ddl_final);
	num_ddl_final.Reallocate(nb_dof_edges);
	for (int k = 0; k < nb_dof_loc; k++)
	  if (intern_node(k) >= 0)
	    num_ddl_final(intern_node(k)) = num_ddl(k);
      }
    
    if (store_block)
      {
	// we replace a_12 by a_12 a_22^{-1}
	a_tmp = a12;
	Mlt(a_tmp, a22, a12);
      }
  }
  

  //! copies arrays needed to compress the matrix
  template<class T> template<class T0>
  void CondensationBlockSolver_Fem<T>::CopyPtr(const CondensationBlockSolver_Fem<T0>& solver)
  {
    nb_scalar_dof = solver.nb_scalar_dof;
    nb_global_dof = solver.nb_global_dof;
    size_compressed_system = solver.size_compressed_system;

    IndexCondensedRows = solver.IndexCondensedRows;
    GlobalCondensedRows = solver.GlobalCondensedRows;
    OverlappedRows = solver.OverlappedRows;
    OverlappedProcs = solver.OverlappedProcs;
    SharingProcs = solver.SharingProcs;
    SharingRows = solver.SharingRows;
  }
  

  //! overloaded method
  template<class T>
  void CondensationBlockSolver_Fem<T>::Init(const GlobalGenericMatrix<T>& nat_mat, bool compress_matrix)
  {
    this->level_volume.SetNbElt(var_computation.GetNbElt(), var_boundary.GetNbEltPML());

    if (compress_matrix)
      {
	nb_scalar_dof = 0;
	nb_global_dof = 0;
	var_computation.GetStaticCondensedRows(IndexCondensedRows_store, GlobalCondensedRows_store,
					       OverlappedRows_store, OverlappedProcs_store,
					       nb_scalar_dof, nb_global_dof,
					       SharingProcs_store, SharingRows_store);

	IndexCondensedRows = &IndexCondensedRows_store;
	GlobalCondensedRows = &GlobalCondensedRows_store;
	OverlappedRows = &OverlappedRows_store;
	OverlappedProcs = &OverlappedProcs_store;
	SharingProcs = &SharingProcs_store;
	SharingRows = &SharingRows_store;
      }
  }
  
  
  //! modifies the elementary matrix because of static condensation
  template<class T>
  void CondensationBlockSolver_Fem<T>
  ::ModifyElementaryMatrix(int ic, IVect& num_ddl, VirtualMatrix<T>& mat_interac,
			   const GlobalGenericMatrix<T>& nat_mat)
  {
    if (optim_condensation)
      return;
    
    int nb_dof_loc = num_ddl.GetM();
    int i = this->GetGlobalElementNumber();
    
    if (var_computation.GetLeafStaticCondensation())
      {
        Matrix<T>& A = dynamic_cast<Matrix<T>& >(mat_interac);
	Matrix<T> a12, a21; Vector<int> num_ddl_kept(num_ddl);
	if (symmetric_elem_matrix)
	  {
	    Matrix<T, Symmetric, RowSymPacked> a11, a22;
	    GetSchurMatrixGen(i, nb_dof_loc, A, a12, a21, a11, a22, num_ddl_kept);
	    if (store_block)
	      {
		int nb_elt = var_computation.GetNbElt();
		if (block_sym.GetM() != nb_elt)
		  block_sym.Reallocate(nb_elt);		    

		a11.Clear();
		block_sym(ic).inv_a22 = a22;
		block_sym(ic).a12 = a12;
		block_sym(ic).a21 = a21;
		block_sym(ic).num_ddl = num_ddl;
	      }
	  }
	else
	  {
	    Matrix<T, General, RowMajor> a11, a22;
	    GetSchurMatrixGen(i, nb_dof_loc, A, a12, a21, a11, a22, num_ddl_kept);
	    if (store_block)
	      {
		int nb_elt = var_computation.GetNbElt();
		if (block_unsym.GetM() != nb_elt)
		  block_unsym.Reallocate(nb_elt);		    
		
		a11.Clear();
		block_unsym(ic).inv_a22 = a22;
		block_unsym(ic).a12 = a12;
		block_unsym(ic).a21 = a21;
		block_unsym(ic).num_ddl = num_ddl;
	      }
	  }

	if (compress_system)
	  num_ddl = num_ddl_kept;
      }
  }


  //! modifies the original source of the non-condensed system 
  template<class T> template<class Prop0, class Storage0, class Prop, class Storage>
  void CondensationBlockSolver_Fem<T>
  ::ModifyRhsGen(const SeldonTranspose& trans, Vector<T>& b_source, int i, const IVect& num_ddl, Matrix<T, Prop0, Storage0>& mat_elt_unsym,
		 Matrix<T>& a12, Matrix<T>& a21, Matrix<T, Prop, Storage>& a22) const
  {
    if (optim_condensation)
      return;

    int nb_dof_loc = num_ddl.GetM();

    Vector<int> intern_node;
    int nb_dof_edges, nb_dof_int;
    var_computation.GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges, nb_dof_int, intern_node);

    if (nb_dof_int == 0)
      return;
    
    Vector<T> b, f;
    b.Reallocate(nb_dof_int); b.Zero();
    f.Reallocate(nb_dof_edges); f.Zero();
    
    if (!store_block)
      {
	// case where the elementary matrices are not stored
	Vector<int> ipivot;
	
	// A = [a11,a12;a21,a22] where a22 refer to inside-inside node matrix
	// schur complement provide (a11-a12 a22^{-1} a21) x = f_1 - a12 a22^{-1} f_2
	a22.Reallocate(nb_dof_int, nb_dof_int);
	
	a12.Reallocate(nb_dof_edges, nb_dof_int);
	for (int j = 0; j < nb_dof_loc; j++)
	  if (intern_node(j) < 0)
	    {
	      int num_j = -intern_node(j)-1;
	      for (int k = 0; k < nb_dof_loc; k++)
		{
		  if (intern_node(k) < 0)
		    a22(num_j, -intern_node(k)-1) = mat_elt_unsym(j, k);
		  else
		    a12(intern_node(k), num_j) = mat_elt_unsym(k, j);
		}
	    }
	
	// new right hand side : f_1 - a12 a22^{-1} f_2
	Seldon::GetLU(a22, ipivot);
	
	for (int j = 0; j < nb_dof_loc; j++)
	  if (intern_node(j) < 0)
	    {
	      int num_dof = num_ddl(j);
	      if (num_dof >= 0)
		b(-intern_node(j)-1) = b_source(num_dof);
	    }
	
	Seldon::SolveLU(a22, ipivot, b);
	
	Seldon::Mlt(a12, b, f);
	
	for (int k = 0; k < nb_dof_loc; k++)
	  if (intern_node(k) >= 0)
	    {
	      int num_dof = num_ddl(k);
	      if (num_dof >= 0)
		if ((var_problem.FormulationDG() == ElementReference_Base::HDG)
		    || (!var_boundary.IsDofDirichlet(num_dof)) )
		  b_source(num_dof) -= f(intern_node(k));
	    }
      }
    else
      {
	// case where the elementary blocks are stored
	for (int j = 0; j < nb_dof_loc; j++)
	  if (intern_node(j) < 0)
	    {
	      int num_dof = num_ddl(j);
	      if (num_dof >= 0)
		b(-intern_node(j)-1) = b_source(num_dof);
	    }

	Mlt(a12, b, f);
	for (int k = 0; k < nb_dof_loc; k++)
	  if (intern_node(k) >= 0)
	    {
	      int num_dof = num_ddl(k);
	      if (num_dof >= 0)
		if ((var_problem.FormulationDG() == ElementReference_Base::HDG)
		    || (!var_boundary.IsDofDirichlet(num_dof)) )
		  b_source(num_dof) -= f(intern_node(k));
	    }	
      }
  }
  
  
  //! modifying right hand side because of static condensation
  /*!
    \param[in,out] b_source right hand side to modify
  */
  template<class T>
  void CondensationBlockSolver_Fem<T>
  ::ModifyRhsStaticCondensation(const SeldonTranspose& trans,
                                Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const
  {
    if (optim_condensation)
      return;

    IVect num_ddl;
    Matrix<T> mat_elt_unsym;
    
    int nb_elt = var_computation.GetNbElt();
    Matrix<T, Symmetric, RowSymPacked> a22_sym;
    Matrix<T> a12, a21, a22_unsym;

    // static condensation, the rhs is modified so that, inside nodes are eliminated
    var_computation.SetLeafStaticCondensation(true);
    for (int i = 0; i < nb_elt; i++)
      {	
	if (!store_block)
	  var_computation.ComputeElementaryMatrix(i, num_ddl, mat_elt_unsym,
						  const_cast<CondensationBlockSolver_Fem<T>& >(*this), nat_mat);

	if (symmetric_elem_matrix)
	  {
	    if (store_block)
	      this->ModifyRhsGen(trans, b_source, i, block_sym(i).num_ddl, mat_elt_unsym,
				 block_sym(i).a12, block_sym(i).a21, block_sym(i).inv_a22);
	    else
	      this->ModifyRhsGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_sym);
	  }
	else
	  {
	    if (store_block)
	      this->ModifyRhsGen(trans, b_source, i, block_unsym(i).num_ddl, mat_elt_unsym,
				 block_unsym(i).a12, block_unsym(i).a21, block_unsym(i).inv_a22);
	    else
	      this->ModifyRhsGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_unsym);
	  }
      }

    var_computation.SetLeafStaticCondensation(false);
  }
  
  
  //! recomposition of the solution (on internal dofs) because of static condensation
  /*!
    \param[in,out] x_sol solution to recompose
   */
  template<class T> template<class Prop0, class Storage0, class Prop, class Storage>
  void CondensationBlockSolver_Fem<T>
  ::RecomposeSolGen(const SeldonTranspose&, Vector<T>& x_sol, int i, const IVect& num_ddl, Matrix<T, Prop0, Storage0>& mat_elt_unsym,
		    Matrix<T>& a12, Matrix<T>& a21, Matrix<T, Prop, Storage>& a22) const
  {
    if (optim_condensation)
      return;

    int nb_dof_loc = num_ddl.GetM();

    Vector<int> intern_node;
    int nb_dof_edges, nb_dof_int;
    var_computation.GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges, nb_dof_int, intern_node);

    if (nb_dof_int == 0)
      return;
	
    Vector<T> x, y, f;
    y.Reallocate(nb_dof_int); y.Zero();
    x.Reallocate(nb_dof_edges); x.Zero();

    T one; SetComplexOne(one);
    if (!store_block)
      {
	Vector<int> ipivot;
	a22.Reallocate(nb_dof_int, nb_dof_int);
	
	a21.Reallocate(nb_dof_int, nb_dof_edges);
	for (int j = 0; j < nb_dof_loc; j++)
	  if (intern_node(j) < 0)
	    {
	      int num_j = -intern_node(j)-1;
	      for (int k = 0; k < nb_dof_loc; k++)
		{
		  if (intern_node(k) < 0)
		    a22(num_j, -intern_node(k)-1) = mat_elt_unsym(j, k);
		  else
		    a21(num_j, intern_node(k)) = mat_elt_unsym(j, k);
		}
	    }
	
	Seldon::GetLU(a22, ipivot);
	
	// y = f_2 - a22^{-1} a_21 x
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    int num_dof = num_ddl(j);
	    if (num_dof >= 0)
	      {		
		if (intern_node(j) >= 0)
		  x(intern_node(j)) = x_sol(num_dof);
		else
		  y(-intern_node(j)-1) = x_sol(num_dof);
	      }
	  }
	
	MltAdd(-one, a21, x, one, y);
	Seldon::SolveLU(a22, ipivot, y);
	    
	for (int k = 0; k < nb_dof_loc; k++)
	  if (intern_node(k) < 0)
	    {
	      int num_dof = num_ddl(k);
	      if (num_dof >= 0)
		x_sol(num_dof) = y(-intern_node(k)-1);
	  }
      }
    else
      {
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    int num_dof = num_ddl(j);
	    if (num_dof >= 0)
	      {		
		if (intern_node(j) >= 0)
		  x(intern_node(j)) = x_sol(num_dof);
		else
		  y(-intern_node(j)-1) = x_sol(num_dof);
	      }
	  }
	
	MltAdd(-one, a21, x, one, y);

	f = y;
	Mlt(a22, f, y);
	
	for (int k = 0; k < nb_dof_loc; k++)
	  if (intern_node(k) < 0)
	    {
	      int num_dof = num_ddl(k);
	      if (num_dof >= 0)
		x_sol(num_dof) = y(-intern_node(k)-1);
	    }
      }
  }

  
  //! modifying right hand side because of static condensation
  /*!
    \param[in,out] b_source right hand side to modify
  */
  template<class T>
  void CondensationBlockSolver_Fem<T>
  ::RecomposeSolution(const SeldonTranspose& trans, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const
  {
    if (optim_condensation)
      return;

    IVect num_ddl;
    Matrix<T> mat_elt_unsym;
    
    int nb_elt = var_computation.GetNbElt();
    Matrix<T, Symmetric, RowSymPacked> a22_sym;
    Matrix<T> a12, a21, a22_unsym;
    
    // static condensation, the rhs is modified so that, inside nodes are eliminated
    var_computation.SetLeafStaticCondensation(true);
    for (int i = 0; i < nb_elt; i++)
      {	
	if (!store_block)
	  var_computation.ComputeElementaryMatrix(i, num_ddl, mat_elt_unsym,
						  const_cast<CondensationBlockSolver_Fem<T>& >(*this), nat_mat);

	if (symmetric_elem_matrix)
	  {
	    if (store_block)
	      this->RecomposeSolGen(trans, b_source, i, block_sym(i).num_ddl, mat_elt_unsym,
				    block_sym(i).a12, block_sym(i).a21, block_sym(i).inv_a22);
	    else
	      this->RecomposeSolGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_sym);
	  }
	else
	  {
	    if (store_block)
	      this->RecomposeSolGen(trans, b_source, i, block_unsym(i).num_ddl, mat_elt_unsym,
				    block_unsym(i).a12, block_unsym(i).a21, block_unsym(i).inv_a22);
	    else
	      this->RecomposeSolGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_unsym);
	  }
      }
    
    var_computation.SetLeafStaticCondensation(false);
  }


  //! returns the memory used by the object in bytes
  template<class T>
  size_t CondensationBlockSolver_Fem<T>::GetMemorySize() const
  {
    size_t taille = 48;
    taille += IndexCondensedRows_store.GetMemorySize() + GlobalCondensedRows_store.GetMemorySize()
      + OverlappedRows_store.GetMemorySize() + OverlappedProcs_store.GetMemorySize()
      + SharingProcs_store.GetMemorySize() + Seldon::GetMemorySize(SharingRows_store);

    taille += block_unsym.GetM()*sizeof(void*);
    for (int i = 0; i < block_unsym.GetM(); i++)
      taille += block_unsym(i).GetMemorySize();

    taille += block_sym.GetM()*sizeof(void*);
    for (int i = 0; i < block_sym.GetM(); i++)
      taille += block_sym(i).GetMemorySize();

    return taille;
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_SOLVE_HARMONIC_CXX
#endif
