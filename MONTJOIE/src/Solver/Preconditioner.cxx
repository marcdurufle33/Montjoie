#ifndef MONTJOIE_FILE_PRECONDITIONER_CXX

namespace Montjoie
{
    
  
#ifdef SELDON_WITH_PRECONDITIONING
  
  ///////////////////////////
  // JacobiPreconditioning //
  
  
  //! default constructor
  template<class T>
  JacobiPreconditioning<T>::JacobiPreconditioning()
  {
    omega = 1;
    size_block = 1;
    nb_iterations = 1;
  }
  
  
  //! returns memory used by the object in bytes
  template<class T>
  size_t JacobiPreconditioning<T>::GetMemorySize() const
  {
    size_t taille = sizeof(T)*diagonal.GetM() + sizeof(*this);
    taille += block_diag_unsym.GetMemorySize();
    taille += block_diag_sym.GetMemorySize();
    return taille;
  }
  
  
  //! no additional parameter
  template<class T>
  void JacobiPreconditioning<T>::SetInputData(const string& keyword, const VectString& param)
  {
  }
  
  
  //! line of a data file
  template<class T>
  void JacobiPreconditioning<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param)
  {
    if (param.GetM() <= 0)
      {
	cout << "In SetInputPreconditioning of JacobiPreconditioning" << endl;
	cout << "JACOBI needs more parameters, for instance :" << endl;
	cout << "JACOBI omega" << endl;
	cout << "Current parameters are : " << endl << param << endl;
	abort();
      }

    nb_iterations = 1;
    if (!keyword.compare("BLOCK_JACOBI"))
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputPreconditioning of JacobiPreconditioning" << endl;
	    cout << "BLOCK_JACOBI needs more parameters, for instance :" << endl;
	    cout << "BLOCK_JACOBI size_block omega" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
	
	size_block = to_num<int>(param(0));
	omega = to_num<Real_wp>(param(1));
	if (param.GetM() > 2)
	  nb_iterations = to_num<int>(param(2));
      }
    else
      {
	size_block = 1;
	omega = to_num<Real_wp>(param(0));
	if (param.GetM() > 1)
	  nb_iterations = to_num<int>(param(1));
      }
  }
  
  
  //! computation of the diagonal of the matrix (or block-diagonal)
  template<class T> template<class MatrixSparse>
  void JacobiPreconditioning<T>::
  ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
			  const MatrixSparse& mat, const GlobalGenericMatrix<T>& nat_mat,
			  CondensationBlockSolver_Fem<T>& cond_solver)
  {
    if (size_block == 1)
      {
	var.ComputeDiagonalMatrix(diagonal, mat, nat_mat);
	
        for (int i = 0; i < diagonal.GetM(); i++)
	  diagonal(i) = omega/diagonal(i);
      }
    else
      {
	block_diag_sym.Clear();
	block_diag_unsym.Clear();
	
	if (IsSymmetricMatrix(mat))
	  {
	    var.ComputeBlockDiagonalMatrix(block_diag_sym, size_block, mat, nat_mat);
	    GetInverse(block_diag_sym); Mlt(omega, block_diag_sym);
	  }
	else
	  {
	    var.ComputeBlockDiagonalMatrix(block_diag_unsym, size_block, mat, nat_mat);
	    GetInverse(block_diag_unsym); Mlt(omega, block_diag_unsym);
	  }
      }
  }
  
  
  //! application of preconditioning
  /*!
    If init_guess_null is null, x is assumed to be 0, such that the application
    of preconditioning is linear (x = M^-1 b where M is a matrix)
    If init_guess_null is false, x is different from 0, and we start from this value
    of x to compute the new value with Jacobi algorithm, in this case the preconditioning
    is non linear (and not adapted to be used directly with iterative solvers)
   */
  template<class T> template<class T0>
  void JacobiPreconditioning<T>::SolveGen(const SeldonTranspose& trans,
					  const VirtualMatrix<T>& A, const Vector<T0>& b,
					  Vector<T0>& x, bool init_guess_null)
  {
    T one; SetComplexOne(one);
    if ((nb_iterations == 1)&&(init_guess_null))
      {
	// no need of an intermediary vector
	if (size_block > 1)
	  {
	    if (block_diag_sym.GetM() > 0)
	      Mlt(trans, block_diag_sym, b, x);
	    else
	      Mlt(trans, block_diag_unsym, b, x);
	  }
	else
	  for (int i = 0; i < x.GetM(); i++)
	    x(i) = b(i)*diagonal(i);
      }
    else
      {
	Vector<T0> r(b);
	for (int n = 0; n < nb_iterations; n++)
	  {
	    if (init_guess_null)
	      {
		if (size_block > 1)
		  {
		    if (block_diag_sym.GetM() > 0)
		      Mlt(trans, block_diag_sym, b, x);
		    else
		      Mlt(trans, block_diag_unsym, b, x);
		  }
		else
		  for (int i = 0; i < x.GetM(); i++)
                    x(i) = b(i)*diagonal(i);
	      }
	    else
	      {
		Copy(b, r);
		A.MltAddVector(-one, trans, x, one, r);
		if (size_block > 1)
		  {
		    if (block_diag_sym.GetM() > 0)
		      MltAdd(one, trans, block_diag_sym, b, one, x);
		    else
		      MltAdd(one, trans, block_diag_unsym, b, one, x);
		  }
		else
		  for (int i = 0; i < x.GetM(); i++)
		    x(i) += r(i)*diagonal(i);
	      }
	    
	    init_guess_null = false;
	  }
      }
  }
  

  //! applies preconditioning (damped Jacobi iteration)
  template<class T>
  void JacobiPreconditioning<T>::Solve(const VirtualMatrix<T>& A, const Vector<T>& b,
				       Vector<T>& x, bool init_guess_null)
  {
    SolveGen(SeldonNoTrans, A, b, x, init_guess_null);
  }
 

  //! application of transpose preconditioning  
  template<class T>
  void JacobiPreconditioning<T>::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b,
					    Vector<T>& x, bool init_guess_null)
  {
    SolveGen(SeldonTrans, A, b, x, init_guess_null);
  }

  
  //! Applies Jacobi preconditioning x = M^-1 b
  template<class T>
  void JacobiPreconditioning<T>::Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)  
  {
    Solve(A, b, x, true);
  }


  //! Applies Jacobi preconditioning x = M^-T b
  template<class T>
  void JacobiPreconditioning<T>::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)  
  {
    TransSolve(A, b, x, true);
  }
  

  //! copies parameters of another jacobi preconditioning
  template<class T>
  void JacobiPreconditioning<T>::CopyParameter(const Preconditioner_Base<T>& precond)
  {
    const JacobiPreconditioning<T>& prec = dynamic_cast<const JacobiPreconditioning<T>& >(precond);
	
    omega = prec.omega;
    nb_iterations = prec.nb_iterations;
    size_block = prec.size_block;
  }
  
  
  //! clears memory used by preconditioning
  template<class T>
  void JacobiPreconditioning<T>::Clear()
  {
    diagonal.Clear();
    block_diag_sym.Clear();
    block_diag_unsym.Clear();
  }
  
    
  // JacobiPreconditioning //
  ///////////////////////////


  //////////////////////////////
  // SubdomainPreconditioning //


  template<class T>
  void SubdomainPreconditioning_Base<T>::InitDefaultValues()
  {
    SetComplexOne(alpha_laplace);

    // default direct solver -> MUMPS
    param_direct_solver.Reallocate(1);
    param_direct_solver(0) = All_MatrixLU<T>::GetAvailableSolver();
  }
  
  
  //! returns memory used by the object in bytes
  template<class T>
  size_t SubdomainPreconditioning_Base<T>::GetMemorySize() const
  {
    size_t taille = OperatorP.GetM()*sizeof(Real_wp) + mat_lu.GetMemorySize();
    taille += Seldon::GetMemorySize(param_direct_solver) + sizeof(*this);
    return taille;
  }
  
  
  //! additional parameters for preconditioning
  template<class T>
  void SubdomainPreconditioning_Base<T>
  ::SetInputData(const string& keyword, const VectString& param)
  {
    if (!keyword.compare("DampingParameters"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of SubdomainPreconditioning" << endl;
	    cout << "DampingParameters needs more parameters, for instance :" << endl;
	    cout << "DampingParameters = alpha" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	alpha_laplace = to_num<T>(param(0));
      }
    else if (!keyword.compare("DirectSolver"))
      {
	param_direct_solver = param;
      }

  }
  
  
  //! parameters present in the line TypeResolution = GMRES SUBDOMAIN
  template<class T>
  void SubdomainPreconditioning_Base<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param)
  {
  }
    

  //! copies parameters of a similar object
  template<class T>
  void SubdomainPreconditioning_Base<T>::CopyParameter(const Preconditioner_Base<T>& prec)
  {
    const SubdomainPreconditioning_Base<T>& prec_sub
      = dynamic_cast<const SubdomainPreconditioning_Base<T>& >(prec);

    param_direct_solver = prec_sub.param_direct_solver;
    alpha_laplace = prec_sub.alpha_laplace;    
  }

  
  //! applies preconditioning, x = M^{-1} b
  template<class T> template<class T0>
  void SubdomainPreconditioning_Base<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<T0>& b, Vector<T0>& x)
  {
#ifdef SELDON_WITH_MPI
    for (int i = 0; i < x.GetM(); i++)
      x(i) = OperatorP(i)*b(i);
    
    mat_lu.Solve(trans, x);
    
    var_subdomain.AddDomains(x);
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= OperatorP(i);
    
#else
    cout << "not implemented in sequential " << endl;
    abort();
#endif

  }


  //! applies preconditioning, x = M^{-1} b
  template<class T>
  void SubdomainPreconditioning_Base<T>::
  Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonNoTrans, A, b, x);
  }

  
  //! applies transpose of preconditioning, x = M^{-T} b
  template<class T>
  void SubdomainPreconditioning_Base<T>::
  TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonTrans, A, b, x);
  } 
  

  //! constructs preconditioning
  template<class T, class Dimension>
  void SubdomainPreconditioning_Dim<T, Dimension>
  ::ConstructPreconditioner(VarComputationProblem& var_computation, All_LinearSolver& solver,
			    const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nature,
			    CondensationBlockSolver_Fem<T>& cond_solver)
  {
#ifdef SELDON_WITH_MPI
    // operator P
    this->OperatorP.Reallocate(var_subdomain.GetNbDof());
    this->OperatorP.Fill(1.0);
    var_subdomain.AddDomains(this->OperatorP);
    for (int i = 0; i < var_subdomain.GetNbDof(); i++)
      this->OperatorP(i) = 1.0/this->OperatorP(i);
    
    // resolution of the local problem with absorbing boundary condition
    Vector<bool> ref_neighbor(var_subdomain.mesh.GetNbReferences()+1);
    ref_neighbor.Fill(false);
    
    for (int i = 1; i < ref_neighbor.GetM(); i++)
      if (var_subdomain.mesh.GetBoundaryCondition(i) == BoundaryConditionEnum::LINE_NEIGHBOR)
        {
          ref_neighbor(i) = true;
          var_subdomain.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_ABSORBING);
        }
    
    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(this->alpha_laplace);
    
    this->mat_lu.SelectDirectSolver(this->param_direct_solver);
    if (var_subdomain.IsSymmetricProblem())
      {
        DistributedMatrix<T, Symmetric, ArrayRowSymSparse> mat_direct;
        var_computation.AddMatrixWithBC(mat_direct, nat_mat);
        mat_direct.Init(mat_direct.GetM(), NULL, NULL, NULL,
                        mat_direct.GetM(), var_subdomain.nb_unknowns_scal,
                        NULL, NULL, MPI_COMM_SELF);
        
        this->mat_lu.Factorize(mat_direct);
      }
    else
      {
        DistributedMatrix<T, General, ArrayRowSparse> mat_direct;
        var_computation.AddMatrixWithBC(mat_direct, nat_mat);
        mat_direct.Init(mat_direct.GetM(), NULL, NULL, NULL,
                        mat_direct.GetM(), var_subdomain.nb_unknowns_scal,
                        NULL, NULL, MPI_COMM_SELF);
        
        this->mat_lu.Factorize(mat_direct);
      }
    
    for (int i = 1; i < ref_neighbor.GetM(); i++)
      if (ref_neighbor(i))
        var_subdomain.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_NEIGHBOR);
    
#else
    cout << "not implemented in sequential " << endl;
    abort();
#endif
    
  }
      
  
  // SubdomainPreconditioning //
  //////////////////////////////
  
  
  /////////////////////////////
  // LowOrderPreconditioning //

  
  //! default constructor
  template<class T>
  LowOrderPreconditioning_Base<T>::LowOrderPreconditioning_Base()
  {
  }
  
  
  //! returns size of memory used by the object in bytes
  template<class T>
  size_t LowOrderPreconditioning_Base<T>::GetMemorySize() const
  {
    size_t taille = NumberDof_LowOrder.GetM()*sizeof(int) + sizeof(*this);
    taille += ScalingDof_LowOrder.GetM()*sizeof(Real_wp);
    taille += mat_lu.GetMemorySize() + Seldon::GetMemorySize(param_direct_solver);
    return taille;
  }
  
  
  //! additional parameters
  template<class T>
  void LowOrderPreconditioning_Base<T>
  ::SetInputData(const string& keyword, const VectString& param)
  {
    if (!keyword.compare("DirectSolver"))
      {
	param_direct_solver = param;
      }
  }

  
  //! parameters in the line TypeResolution = GMRES LOW_ORDER
  template<class T>
  void LowOrderPreconditioning_Base<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param)
  {
  }
    

  //! applies preconditioning, x = M^{-1} b (or its transpose)
  template<class T> template<class T0>
  void LowOrderPreconditioning_Base<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	     const Vector<T0>& b, Vector<T0>& x)
  {
    Copy(b, x);
    mat_lu.Solve(trans, x);
  }
  

  //! applies preconditioning, x = M^{-1} b
  template<class T>
  void LowOrderPreconditioning_Base<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonNoTrans, A, b, x);
  }
  
  
  //! applies preconditioning, x = M^{-T} b
  template<class T>
  void LowOrderPreconditioning_Base<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonTrans, A, b, x);
  } 


  template<class T>
  void LowOrderPreconditioning_Base<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x)
  {
    SolveGen(trans, A, b, x);
  }


  template<class T>
  void LowOrderPreconditioning_Base<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x)
  {
    SolveGen(trans, A, b, x);
  }


  //! constructs preconditioning
  template<class T, class Dimension>
  void LowOrderPreconditioning_Dim<T, Dimension>
  ::ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver&,
			    const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nature,
			    CondensationBlockSolver_Fem<T>& cond_solver)
  {
    DistributedProblem<Dimension>* var_low_ptr = var_problem.GetNewEllipticProblem();
    DistributedProblem<Dimension>& var_low = *var_low_ptr;
    
    // subdividing the initial mesh
    var_low.CopyInputData(var_problem);
    var_low.mesh = var_problem.mesh;
    var_low.mesh_data.Clear();
    MeshNumbering<Dimension>& mesh_num_low = var_low.GetMeshNumbering(0);
    mesh_num_low.Clear();
    
    mesh_num_low.SetOrder(1);    

    MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    
    // using Gauss-Lobatto points as subdivision
    VectReal_wp points_div, omega;
    ComputeGaussLobatto(points_div, omega, mesh_num.GetOrder());
    var_low.mesh.SubdivideMesh(points_div);
    
    var_low.mesh.SetGeometryOrder(1);
    var_low.DoNotComputeGrid();
    var_low.ComputeMeshAndFiniteElement(var_problem.finite_element_name);
    var_low.ComputeMassMatrix();
    
    // for Maxwell's equations the dof numbers are different
    this->ComputeLocalProlongation(points_div, var_problem, var_low);
    
    // solving the problem with low order
    All_LinearSolver* solver_low;
    solver_low = var_low.GetComputationProblem().GetNewLinearSolver();
    
    this->mat_lu.SelectDirectSolver(this->param_direct_solver);
    var_low.GetComputationProblem().SetThresholdMatrix(this->mat_lu.GetThresholdMatrix());
    
    if (IsSymmetricMatrix(A))
      {
	DistributedMatrix<T, Symmetric, ArrayRowSymSparse> B;
	solver_low->ComputeMatrixAndFactorize(B, this->mat_lu, nature);
      }
    else
      {
	DistributedMatrix<T, General, ArrayRowSparse> B;
	solver_low->ComputeMatrixAndFactorize(B, this->mat_lu, nature);
      }
    
    delete solver_low;
    delete var_low_ptr;
  }
  
    
  // LowOrderPreconditioning //
  /////////////////////////////
  

  //////////////////////////////
  // MultigridPreconditioning //

  
  //! Sets values of attributes to default values
  template<class T>
  void MultigridPreconditioning_Base<T>::InitDefaultValues()
  {
    SetComplexOne(alpha_laplace);
    SetComplexOne(relaxation_laplace);
    coarse_solver = NULL;
    fine_solver = NULL;
    var_coarse = NULL;
    nb_smoothing_iterations = 1;
    nb_subcycles_multigrid = 2;
    use_f_cycle = false;
    minimal_order_multigrid = 1;
    type_smoother = JACOBI;
    mat_ssor_sym = NULL;
    mat_ssor_unsym = NULL;
    matCSR_ssor_sym = NULL;
    matCSR_ssor_unsym = NULL;
    mat_smoothing = NULL;
    cond_solver = NULL;
    mat_stored_smoothing = NULL;
    maximal_order_multigrid = 1;
    no_damping_on_first_smoother = false;
    nb_dof_coarse = 0;
    nb_iterations_multigrid = 1;
    
    // default direct solver -> MUMPS
    param_direct_solver.Reallocate(1);
    param_direct_solver(0) = All_MatrixLU<T>::GetAvailableSolver();
  }
  
  
  //! destructor
  template<class T>
  MultigridPreconditioning_Base<T>::~MultigridPreconditioning_Base()
  {
    if (coarse_solver != NULL)
      {
	delete coarse_solver;
	coarse_solver = NULL;
      }
    
    if (var_coarse != NULL)
      {
	delete var_coarse;
	var_coarse = NULL;
      }
    
    if (mat_stored_smoothing != NULL)
      delete mat_stored_smoothing;
  }
  
  
  //! returns size of memory used by the object in bytes
  template<class T>
  size_t MultigridPreconditioning_Base<T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += Vrestriction.GetMemorySize();
    taille += InvDh.GetM()*sizeof(Real_wp);
    taille += proj_coarse_to_fine.GetMemorySize();
    taille += Seldon::GetMemorySize(param_direct_solver);
    taille += mat_stored_sym.GetMemorySize() + mat_stored_unsym.GetMemorySize();
    taille += smoother_jacobi.GetMemorySize();
    if (mat_stored_smoothing != NULL)
      taille += mat_stored_smoothing->GetMemorySize();
    
    if (coarse_solver != NULL)
      taille += coarse_solver->GetMemorySize();

    if (var_coarse != NULL)
      taille += var_coarse->GetMemorySize();
        
    return taille;
  }
  
  
  //! parameters in the line TypeResolution = GMRES MULTIGRID parameters
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param)
  {
    if (!keyword.compare("MULTIGRID"))
      {
	if (param.GetM() <= 2)
	  {
	    cout << "In SetInputPreconditioning of MultigridPreconditioning" << endl;
	    cout << "MULTIGRID needs more parameters, for instance :" << endl;
	    cout << "MULTIGRID nb_subcycles rmin nb_smooth_iter" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	if (param(0) == "V")
	  nb_subcycles_multigrid = 1;
	else if (param(0) == "W")
	  nb_subcycles_multigrid = 2;
	else if (param(0) == "F")
	  {
	    nb_subcycles_multigrid = 2;
	    use_f_cycle = true;
	  }
	else
	  to_num(param(0), nb_subcycles_multigrid);
	
	to_num(param(1), minimal_order_multigrid);
	to_num(param(2), nb_smoothing_iterations);
      }
  }
    
  
  //! additional parameters
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::SetInputData(const string& keyword, const VectString& param)
  {
    if (!keyword.compare("Smoother"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of MultigridPreconditioning" << endl;
	    cout << "Smoother needs more parameters, for instance :" << endl;
	    cout << "Smoother = JACOBI" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	VectString other_param(param.GetM()-1);
	for (int i = 1; i < param.GetM(); i++)
	  other_param(i-1) = param(i);
	
	if (!param(0).compare("JACOBI"))
	  {
	    type_smoother = JACOBI;
	    smoother_jacobi.SetInputPreconditioning(param(0), other_param);
	  }
	else if (!param(0).compare("BLOCK_JACOBI"))
	  {
	    type_smoother = BLOCK_JACOBI;
	    smoother_jacobi.SetInputPreconditioning(param(0), other_param);
	  }
	else if (!param(0).compare("SSOR"))
	  {
	    if (param.GetM() <= 1)
	      {
		cout << "In SetInputData of MultigridPreconditioning" << endl;
		cout << "Smoother needs more parameters, for instance :" << endl;
		cout << "Smoother = SSOR omega" << endl;
		cout << "Current parameters are : " << endl << param << endl;
		abort();
	      }

	    type_smoother = SSOR;
	    smoother_ssor.SetParameterRelaxation(to_num<Real_wp>(param(1)));
	    if (param.GetM() > 2)
	      smoother_ssor.SetNumberIterations(to_num<int>(param(2)));
	  }
      }
    else if (!keyword.compare("DampingParameters"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of MultigridPreconditioning" << endl;
	    cout << "DampingParamters needs more parameters, for instance :" << endl;
	    cout << "DampingParameters = alpha" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	alpha_laplace = to_num<T>(param(0));
	to_num(param(1), relaxation_laplace);
      }
    else if (!keyword.compare("DirectSolver"))
      {
	param_direct_solver = param;
      }
  }
  
  
  //! sets coefficient alpha used to damp original equation
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::SetDampingCoefficient(const T& alpha)
  {
    alpha_laplace = alpha;
  }
  
    
  //! constructs multigrid preconditioning
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::FinalizePreconditioner(VarComputationProblem* var, int rf, int rc,
			   const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nature)
  {
    var_coarse = var;
    if (coarse_solver != NULL)
      delete coarse_solver;
           
    // resolution of coarse problem
    coarse_solver = var_coarse->GetNewLinearSolver();
    coarse_solver->CopyParameter(*fine_solver);
    coarse_solver->EnableRealPreconditioning(false);
    maximal_order_multigrid = max(maximal_order_multigrid, rf);
    
    All_Preconditioner_Base<T>& prec_coarse = coarse_solver->GetPreconditioning(T(0));
    MultigridPreconditioning_Base<T>& mg_coarse = prec_coarse.GetMultigridPreconditioning();

    // selecting solver for coarse problem
    if (rc <= minimal_order_multigrid)
      {
	// direct resolution
	coarse_solver->SetDirectSolver(param_direct_solver);
      }
    else
      {
	// multigrid solver
	coarse_solver->SetIterativeSolver(All_IterativeSolver<T>::PRECOND,
					  All_Preconditioner_Base<T>::MULTIGRID);
	
	mg_coarse.CopyParameter(*this);
      }
    
    // setting maximal order
    mg_coarse.SetMaximalOrder(maximal_order_multigrid);
    mg_coarse.nb_iterations_multigrid = this->nb_subcycles_multigrid;
    
    // effective resolution of coarse problem
    var_coarse->ComputeMassMatrix();
    
    {
      // setting coefficients of coarse problem
      T coef_mass, coef_stiff, coef_damp;
      Vector<string> param(2); param(0) = "1.0"; param(1) = to_str(relaxation_laplace);
      coarse_solver->SetInputData(string("DampingParameters"), param);
      to_complex(nature.GetCoefMass()*relaxation_laplace, coef_mass);
      to_complex(nature.GetCoefStiffness(), coef_stiff);
      to_complex(nature.GetCoefDamping(), coef_damp);
      nat_mat_coarse.SetCoefMass(coef_mass);
      nat_mat_coarse.SetCoefStiffness(coef_stiff);
      nat_mat_coarse.SetCoefDamping(coef_damp);
    }

    coarse_solver->PerformFactorizationStep(nat_mat_coarse);
    coarse_solver->SetPrintLevel(0);

    nb_dof_coarse = var_coarse->GetNbRows();
    if (fine_solver->StaticCondensation())
      nb_dof_coarse = coarse_solver->GetCondensedSolver(T(0)).GetNbRows();

    // smoother
    GlobalGenericMatrix<T> nat_mat = nature;
    if (maximal_order_multigrid == rf)
      {
        if (no_damping_on_first_smoother)
          nat_mat.SetCoefMass(nature.GetCoefMass());
      }

    if (nat_mat.GetCoefMass() == A.GetCoefMass())
      {
        FemMatrixFreeClass_Base<T>& Avar = const_cast<FemMatrixFreeClass_Base<T>&>(A);
        mat_smoothing = &Avar;
      }
    else
      {
	if (mat_stored_smoothing == NULL)
	  mat_stored_smoothing = var_fine.GetNewIterativeMatrix(T(0));

	CondensationBlockSolver_Fem<T>* other_cond
	  = static_cast<CondensationBlockSolver_Fem<T>* >(var_fine.GetNewCondensationSolver(T(0)));

	other_cond->CopyParameter(*cond_solver);
	if (fine_solver->StaticCondensation())
	  var_fine.SetLeafStaticCondensation(true);

	var_fine.AddMatrixWithBC(*mat_stored_smoothing, nat_mat, 0, 0, other_cond);
        mat_smoothing = mat_stored_smoothing;

	if (fine_solver->StaticCondensation())
	  var_fine.SetLeafStaticCondensation(false);

	delete other_cond;
	
	if (fine_solver->StaticCondensation())
	  {
	    cond_solver->CompressMatrixGen(mat_smoothing->mat_iterative_sym);
	    int N = mat_smoothing->mat_iterative_sym.GetM();
	    mat_smoothing->Reallocate(N, N);
	  }
      }

    /* int nb_dof_fine = fine_solver->GetCondensedSolver(T(0)).GetNbRows();
    Vector<T> r(nb_dof_fine), r_c(nb_dof_coarse);
    Matrix<Real_wp, General, ArrayRowSparse> mat(r_c.GetM(), r.GetM());

    for (int i = 0; i < r.GetM(); i++)
      {
	r.Zero();
	r(i) = Real_wp(1);
	
	this->GetRestriction(r, r_c);

	for (int j = 0; j < r_c.GetM(); j++)
	  if (realpart(r_c(j)) != Real_wp(0))
	    mat.AddInteraction(j, i, realpart(r_c(j)));
      }

    mat.WriteText("R.dat");

    mat.Clear();
    mat.Reallocate(r.GetM(), r_c.GetM());
    for (int i = 0; i < r_c.GetM(); i++)
      {
	r_c.Zero();
	r_c(i) = Real_wp(1);
	
	this->GetProlongation(r_c, r);

	for (int j = 0; j < r.GetM(); j++)
	  if (realpart(r(j)) != Real_wp(0))
	    mat.AddInteraction(j, i, realpart(r(j)));
      }

    mat.WriteText("P.dat");
    */
    this->ConstructSmoother(*mat_smoothing, nat_mat);
  }
  
  
  //! computes smoother
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::ConstructSmoother(const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nat_mat)
  {
    mat_stored_sym.Clear(); mat_stored_unsym.Clear();
    if ((type_smoother == JACOBI)||(type_smoother == BLOCK_JACOBI))
      {
	smoother_jacobi.ConstructPreconditioner(var_fine, *fine_solver,
						*mat_smoothing, nat_mat, *cond_solver);
      }
    else if (type_smoother == SSOR)
      {
	smoother_ssor.InitUnSymmetricPreconditioning();
	if (IsSymmetricMatrix(A))
          {
            if (!mat_smoothing->SucceedInAffectingPointer(mat_ssor_sym, matCSR_ssor_sym))
              {
                var_fine.AddMatrixWithBC(*mat_smoothing, mat_stored_sym, nat_mat);
                mat_ssor_sym = &mat_stored_sym;
              }
	    
	    mat_ssor_unsym = NULL;
	    matCSR_ssor_unsym = NULL;
          }
        else
          {
            if (!mat_smoothing->SucceedInAffectingPointer(mat_ssor_unsym, matCSR_ssor_unsym))
              {
                var_fine.AddMatrixWithBC(*mat_smoothing, mat_stored_unsym, nat_mat);
                mat_ssor_unsym = &mat_stored_unsym;
              }
	    
	    mat_ssor_sym = NULL;
	    matCSR_ssor_sym = NULL;
          }
      }
  }
  
  
  //! copies parameters of another multigrid preconditioning
  template<class T>
  void MultigridPreconditioning_Base<T>::
  CopyParameter(const Preconditioner_Base<T>& prec)
  {
    const MultigridPreconditioning_Base<T>& fine_mg = dynamic_cast<const MultigridPreconditioning_Base<T>& >(prec);
    
    alpha_laplace = fine_mg.alpha_laplace*fine_mg.relaxation_laplace;
    relaxation_laplace = fine_mg.relaxation_laplace;
    nb_smoothing_iterations = fine_mg.nb_smoothing_iterations;
    nb_subcycles_multigrid = fine_mg.nb_subcycles_multigrid;
    nb_iterations_multigrid = fine_mg.nb_iterations_multigrid;
    minimal_order_multigrid = fine_mg.minimal_order_multigrid;
    maximal_order_multigrid = fine_mg.maximal_order_multigrid;
    no_damping_on_first_smoother = fine_mg.no_damping_on_first_smoother;
    type_smoother = fine_mg.type_smoother;
    smoother_ssor = fine_mg.smoother_ssor;
    smoother_jacobi = fine_mg.smoother_jacobi;
    param_direct_solver = fine_mg.param_direct_solver;
    use_f_cycle = fine_mg.use_f_cycle;
  }
  
  
  //! applies preconditioning x = M^{-1} b
  template<class T> template<class T0>
  void MultigridPreconditioning_Base<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<T0>& b, Vector<T0>& x,
	     int nb_cycles, int nb_iter)
  {
    if (trans.Trans())
      {
	cout << "not implemented" << endl;
	abort();
      }

    if (nb_cycles == 0)
      nb_cycles = nb_subcycles_multigrid;

    if (nb_iter == 0)
      nb_iter = nb_iterations_multigrid;
    
    Vector<T0> r_coarse(nb_dof_coarse), r(b.GetM());
    r.Zero(); r_coarse.Zero(); x.Zero();
    
    MultiGrid_Vcycle(A, x, b, r_coarse, r, true, nb_cycles);

    if (use_f_cycle)
      nb_cycles = 1;
    
    for (int i = 1; i < nb_iter; i++)
      MultiGrid_Vcycle(A, x, b, r_coarse, r, false, nb_cycles);
  }
  

  //! applies preconditioning x = M^{-1} b
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonNoTrans, A, b, x, nb_subcycles_multigrid, nb_iterations_multigrid);
  }


  //! applies preconditioning x = M^{-T} b
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonTrans, A, b, x, nb_subcycles_multigrid, nb_iterations_multigrid);
  } 
  
  
  //! performs a cycle 
  template<class T> template<class T0>
  void MultigridPreconditioning_Base<T>
  ::MultiGrid_Vcycle(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
		     Vector<T0>& r_coarse, Vector<T0>& r, bool init_guess_null, int n)
  {
    T0 one; SetComplexOne(one);
    Vector<T0> x_coarse(r_coarse); x_coarse.Zero(); 

    //x.WriteText("x_init.dat");
    //b.WriteText("b.dat");
	
    // pre-smoothing step
    this->ApplyPreSmoother(*mat_smoothing, x, b, r, init_guess_null);

    //x.WriteText("x_pre.dat");
    
    // computes residual r = b - A x
    Copy(b, r);
    mat_smoothing->MltAddVector(-one, x, one, r);

    //r.WriteText("residu.dat");
	
    // projects residual to the coarse mesh
    this->GetRestriction(r, r_coarse);

    //r_coarse.WriteText("r_coarse.dat");
	
    if (this->GetCoarseOrder() <= this->minimal_order_multigrid)
      {
        // this is the coarsest grid : we use a direct solver
	x_coarse = r_coarse;
        coarse_solver->GetDirectFactorization(T(0)).Solve(x_coarse);
      }
    else
      {
        // there is a coarser grid
	MultigridPreconditioning_Base<T>& mg_coarse
	  = coarse_solver->GetPreconditioning(T(0)).GetMultigridPreconditioning();

	// solving the coarse problem with multigrid
	mg_coarse.SolveGen(SeldonNoTrans, *mg_coarse.mat_smoothing, r_coarse, x_coarse, n, n);
      }

    //x_coarse.WriteText("sol_coarse.dat");
	
    // projects solution to fine mesh
    this->GetProlongation(x_coarse, r);

    //r.WriteText("sol.dat");
	
    // updating solution
    Add(one, r, x);

    //x.WriteText("x_next.dat");
	
    // post-smoothing step
    this->ApplyPostSmoother(*mat_smoothing, x, b, r, false);

    //x.WriteText("x_post.dat");
    //cout << "waiting" << endl; int test_input; cin >> test_input;
  }
  
  
  //! applies pre-smoothing to vector x
  template<class T> template<class T0>
  void MultigridPreconditioning_Base<T>
  ::ApplyPreSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			Vector<T0>& r, bool init_guess_null)
  {
    if ((type_smoother == JACOBI)||(type_smoother == BLOCK_JACOBI))
      smoother_jacobi.SolveGen(SeldonNoTrans, A, b, x, init_guess_null);
    else if (type_smoother == SSOR)
      {
        if (mat_ssor_sym != NULL)
          smoother_ssor.SolveGen(SeldonNoTrans, *mat_ssor_sym, b, x, init_guess_null);
	else if (matCSR_ssor_sym != NULL)
          smoother_ssor.SolveGen(SeldonNoTrans, *matCSR_ssor_sym, b, x, init_guess_null);
        else if (mat_ssor_unsym != NULL)
          smoother_ssor.SolveGen(SeldonNoTrans, *mat_ssor_unsym, b, x, init_guess_null);
	else
	  smoother_ssor.SolveGen(SeldonNoTrans, *matCSR_ssor_unsym, b, x, init_guess_null);
      }
  }
  
  
  //! applies post-smoothing to vector x
  template<class T> template<class T0>
  void MultigridPreconditioning_Base<T>
  ::ApplyPostSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			 Vector<T0>& r, bool init_guess_null)
  {
    if ((type_smoother == JACOBI)||(type_smoother == BLOCK_JACOBI))
      smoother_jacobi.SolveGen(SeldonTrans, A, b, x, init_guess_null);
    else if (type_smoother == SSOR)
      {
        if (mat_ssor_sym != NULL)
          smoother_ssor.SolveGen(SeldonTrans, *mat_ssor_sym, b, x, init_guess_null);
	else if (matCSR_ssor_sym != NULL)
          smoother_ssor.SolveGen(SeldonTrans, *matCSR_ssor_sym, b, x, init_guess_null);
        else if (mat_ssor_unsym != NULL)
          smoother_ssor.SolveGen(SeldonTrans, *mat_ssor_unsym, b, x, init_guess_null);
	else
          smoother_ssor.SolveGen(SeldonTrans, *matCSR_ssor_unsym, b, x, init_guess_null);
      }
  }


  template<>
  void MultigridPreconditioning_Base<Complex_wp>
  ::ApplyPreSmoother(const VirtualMatrix<Complex_wp>& A, Vector<Real_wp>& x, const Vector<Real_wp>& b,
		     Vector<Real_wp>& r, bool init_guess_null)
  {
    cout << "Incompatible types" << endl;
    abort();
  }
  

  template<class T>
  void MultigridPreconditioning_Base<T>
  ::ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
		     Vector<Treal>& r, bool init_guess_null)
  {
    ApplyPreSmootherGen(A, x, b, r, init_guess_null);
  }
  

  template<class T>
  void MultigridPreconditioning_Base<T>
  ::ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
		     Vector<Tcplx>& r, bool init_guess_null)
  {
    ApplyPreSmootherGen(A, x, b, r, init_guess_null);
  }
  

  template<>
  void MultigridPreconditioning_Base<Complex_wp>
  ::ApplyPostSmoother(const VirtualMatrix<Complex_wp>& A, Vector<Real_wp>& x, const Vector<Real_wp>& b,
		     Vector<Real_wp>& r, bool init_guess_null)
  {
    cout << "Incompatible types" << endl;
    abort();
  }

  
  template<class T>
  void MultigridPreconditioning_Base<T>
  ::ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
		      Vector<Treal>& r, bool init_guess_null)
  {
    ApplyPostSmootherGen(A, x, b, r, init_guess_null);
  }
    

  template<class T>
  void MultigridPreconditioning_Base<T>
  ::ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
		      Vector<Tcplx>& r, bool init_guess_null)
  {
    ApplyPostSmootherGen(A, x, b, r, init_guess_null);
  }


  /********************************
   * MultigridPreconditioning_Dim *
   ********************************/

  
  //! adds local contribution to prolongation operator
  template<class T, class Dimension>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::ComputeProlongationElement(int i, int rf, int rc, Matrix<Real_wp>& LocalProlongation)
  {
    const MeshNumbering<Dimension>& mesh_fine = var_fine.GetMeshNumbering(0);
    const MeshNumbering<Dimension>& mesh_coarse = var_coarse->GetMeshNumbering(0);

    // incrementation of InvDh
    int nb_dof_fine = mesh_fine.GetNbLocalDof(i);
    int nb_dof_coarse = mesh_coarse.GetNbLocalDof(i);
    for (int j = 0; j < nb_dof_fine; j++)
      {
	int num_dof_fine = mesh_fine.Element(i).GetNumberDof(j);
        for (int p = 0; p < var_fine.nb_unknowns; p++)
	  this->InvDh(num_dof_fine + p*mesh_fine.GetNbDof()) += 1;
      }

    //if (rf <= fine_elt.GetMaximalOrderRestriction())
    if (this->Vrestriction.GetM() > 0)
      {
        // retrieving dofs with a negative sign
        Vector<bool> NegativeFine(nb_dof_fine), NegativeCoarse(nb_dof_coarse);
        NegativeFine.Fill(false); NegativeCoarse.Fill(false);
        const Vector<int>& neg_fine = mesh_fine.Element(i).GetNegativeDofNumber();
        const Vector<int>& neg_coarse = mesh_coarse.Element(i).GetNegativeDofNumber();
        for (int j = 0; j < neg_fine.GetM(); j++)
          NegativeFine(neg_fine(j)) = true;

        for (int j = 0; j < neg_coarse.GetM(); j++)
          NegativeCoarse(neg_coarse(j)) = true;
        
        // adding entries to Vrestriction
        Real_wp value;
        for (int j = 0; j < nb_dof_fine; j++)
          {
            for (int k = 0; k < nb_dof_coarse; k++)
              if (abs(LocalProlongation(j, k)) > 1e-12)
                {
                  int num_dof_fine = mesh_fine.Element(i).GetNumberDof(j);
                  int num_dof_coarse = mesh_coarse.Element(i).GetNumberDof(k);
                  value = LocalProlongation(j, k);
                                    
                  if (NegativeFine(j) != NegativeCoarse(k))
                    value = -value;
                  
                  for (int p = 0; p < var_fine.nb_unknowns; p++)
                    {
                      this->Vrestriction.
			AddInteraction(num_dof_coarse + p*mesh_coarse.GetNbDof(),
				       num_dof_fine + p*mesh_fine.GetNbDof(), value);
                    }
                }
          }
      }
  }


  //! constructs preconditioning
  template<class T, class Dimension>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
			    const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nat_mat,
			    CondensationBlockSolver_Fem<T>& cond_solver_)
  {
    this->fine_solver = &solver;
    this->cond_solver = &cond_solver_;
    var_coarse = var_fine.GetNewEllipticProblem();
    
    // computation of problem defined on coarse mesh
    var_coarse->CopyInputData(var_fine);
	
    int rf = var_fine.GetMeshNumbering(0).GetOrder();
    int rc = (rf+1)/2;
    var_coarse->GetMeshNumbering(0).SetOrder(rc);
    var_coarse->mesh.SetGeometryOrder(rc);
    
    var_coarse->DoNotComputeGrid();
    var_coarse->ComputeMeshAndFiniteElement(var_fine.finite_element_name);

    var_coarse->PerformOtherInitializations();
    
    // computation of restriction/prolongation operators
    TinyVector<Matrix<Real_wp>, 4> ProlongationElement;
    
    // on reference elements
    var_fine.ComputeLocalProlongation(this->proj_coarse_to_fine, *var_coarse, rc, ProlongationElement);
    
    int Nf = var_fine.offset_dof_unknown(var_fine.nb_unknowns_scal);
    int Nc = var_coarse->offset_dof_unknown(var_fine.nb_unknowns_scal);

    if (var_fine.FormulationDG() == ElementReference_Base::HDG)
      {
	Nf = var_fine.nb_unknowns_hdg*var_fine.GetNbMeshDof();
	Nc = var_coarse->nb_unknowns_hdg*var_coarse->GetNbMeshDof();
      }
    
    //if (rf <= TypeElement::ORDER_MAX_RESTRICTION_STD)
    if (false)
      {
	this->Vrestriction.Clear();
	this->Vrestriction.Reallocate(Nc, Nf);
      }
    
    // on all the mesh if needed
    this->InvDh.Reallocate(Nf); this->InvDh.Zero();
    if (var_fine.FormulationDG() == ElementReference_Base::HDG)
      {
	// not needed for HDG
	/*this->Vrestriction.Clear();
	this->Vrestriction.Reallocate(Nc, Nf);

	DISP(ProlongationElement(0));
	for (int i = 0; i < var_fine.mesh.GetNbBoundary(); i++)
	  {
	    
	  }*/
      }
    else
      for (int i = 0; i < var_fine.mesh.GetNbElt(); i++)
	{
	  int type_elt = var_fine.mesh.GetTypeElement(i);
	  ComputeProlongationElement(i, rf, rc, ProlongationElement(type_elt));
	}
    
    // assembling InvDh
    var_fine.AddDomains(this->InvDh);
    
    for (int i = 0; i < this->InvDh.GetM(); i++)
      this->InvDh(i) = Real_wp(1)/this->InvDh(i);
    
    //if (rf <= TypeElement::ORDER_MAX_RESTRICTION_STD)
    if (false)
      {	
	for (int i = 0; i < this->Vrestriction.GetM(); i++)
	  for (int j = 0; j < this->Vrestriction.GetRowSize(i); j++)
	    this->Vrestriction.Value(i, j) *= this->InvDh(this->Vrestriction.Index(i,j));
	
      }

    this->FinalizePreconditioner(&var_coarse->GetComputationProblem(), rf, rc, A, nat_mat);
  }
  

  //! computes restriction of r on the coarse mesh
  template<class T, class Dimension> template<class T0>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::GetRestrictionGen(const Vector<T0>& r0, Vector<T0>& r0_coarse) const
  {
    T0 zero, one; SetComplexZero(zero); SetComplexOne(one);
    r0_coarse.Zero();
    if (this->Vrestriction.GetM() > 0)
      {
	Mlt(this->Vrestriction, r0, r0_coarse);
      }
    else if (var_fine.FormulationDG() == ElementReference_Base::HDG)
      {
	Vector<T0> Ufine, Ucoarse;
	int offset_coarse = 0, offset_fine = 0;
	for (int i = 0; i < var_fine.mesh.GetNbBoundary(); i++)
	  {
	    int nb_dof_fine = var_fine.GetNbSurfaceDof(i);
	    int nb_dof_coarse = var_coarse->GetNbSurfaceDof(i);
	    int type_elt = var_fine.mesh.Boundary(i).GetHybridType();
	    Ufine.Reallocate(nb_dof_fine); Ucoarse.Reallocate(nb_dof_coarse);
	    for (int j = 0; j < nb_dof_fine; j++)
	      Ufine(j) = r0(offset_fine + j);
	    
	    this->proj_coarse_to_fine.TransposeProjectScalar(Ufine, Ucoarse, type_elt);

	    for (int j = 0; j < nb_dof_coarse; j++)
	      r0_coarse(offset_coarse + j) = Ucoarse(j);
	    
	    offset_fine += nb_dof_fine;
	    offset_coarse += nb_dof_coarse;	    
	  }
      }
    else
      {
	Vector<T0> r(var_fine.offset_dof_unknown(var_fine.nb_unknowns_scal));
	Vector<T0> r_coarse(var_coarse->offset_dof_unknown(var_fine.nb_unknowns_scal));
	r.Zero(); r_coarse.Zero();
	if ((this->fine_solver->StaticCondensation()) && (!var_fine.GetComputationProblem().LightStaticCondensation()))
	  this->cond_solver->ScatterVector(r0, r);
	else
	  r = r0;
	
	int nb_unknowns = var_fine.nb_unknowns;
	Vector<Vector<T0> > Ufine(nb_unknowns), Vcoarse(nb_unknowns);
	const MeshNumbering<Dimension>& mesh_coarse = var_coarse->GetMeshNumbering(0);
	const MeshNumbering<Dimension>& mesh_fine = var_fine.GetMeshNumbering(0);
	int Nvol_fine = mesh_fine.GetNbDof();
	int Nvol_coarse = mesh_coarse.GetNbDof();
	for (int e = 0; e < var_fine.mesh.GetNbElt(); e++)
	  {
	    int nb_dof_coarse = mesh_coarse.GetNbLocalDof(e);
	    int nb_dof_fine = mesh_fine.GetNbLocalDof(e);
            for (int m = 0; m < nb_unknowns; m++) 
	      {
		Vcoarse(m).Reallocate(nb_dof_coarse);
		Ufine(m).Reallocate(nb_dof_fine);
	      }
	    
	    // we get local variable U from global variable r
	    for (int i = 0; i < nb_dof_fine; i++)
	      {
		int num_dof = mesh_fine.Element(e).GetNumberDof(i);
		for (int m = 0; m < nb_unknowns; m++) 
		  Ufine(m)(i) = r(num_dof + m*Nvol_fine)*this->InvDh(num_dof);
	      }
   
            mesh_fine.number_map.ModifyLocalUnknownVector(mesh_fine, Ufine, e);
            
	    // local restriction
	    int type_elt = var_fine.mesh.GetTypeElement(e);
	    for (int m = 0; m < nb_unknowns; m++)
	      this->proj_coarse_to_fine.TransposeProjectScalar(Ufine(m), Vcoarse(m), type_elt);            
            
	    mesh_coarse.number_map.ModifyLocalUnknownVector(mesh_coarse, Vcoarse, e);
            
	    // we put local variable Vcoarse to global variable r_coarse
	    for (int i = 0; i < nb_dof_coarse; i++)
	      {
		int num_dof = mesh_coarse.Element(e).GetNumberDof(i);
		for (int m = 0; m < nb_unknowns; m++) 
		  r_coarse(num_dof+m*Nvol_coarse) += Vcoarse(m)(i);
	      }
	  }

	DistributedProblem_Base& var_comm = dynamic_cast<DistributedProblem_Base&>(*var_coarse);
	var_comm.AddDomains(r_coarse);
	
	if ((this->fine_solver->StaticCondensation()) && (!var_fine.GetComputationProblem().LightStaticCondensation()))
	  this->coarse_solver->GetCondensedSolver(T(0)).ExtractVector(r_coarse, r0_coarse);
	else
	  r0_coarse = r_coarse;
	
      }
	
    var_coarse->GetBoundaryConditionProblem().ImposeNullDirichletCondition(r0_coarse);
  }
  
  
  //! computes prolongation of r_coarse on a fine mesh
  template<class T, class Dimension> template<class T0>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::GetProlongationGen(const Vector<T0>& r0_coarse, Vector<T0>& r0) const
  {
    T0 zero, one; SetComplexZero(zero); SetComplexOne(one);
    r0.Zero();
    if (this->Vrestriction.GetM() > 0)
      {
	Mlt(SeldonTrans, this->Vrestriction, r0_coarse, r0);
      }
    else if (var_fine.FormulationDG() == ElementReference_Base::HDG)
      {
	Vector<T0> Ufine, Ucoarse;
	int offset_coarse = 0, offset_fine = 0;
	for (int i = 0; i < var_fine.mesh.GetNbBoundary(); i++)
	  {
	    int nb_dof_fine = var_fine.GetNbSurfaceDof(i);
	    int nb_dof_coarse = var_coarse->GetNbSurfaceDof(i);
	    int type_elt = var_fine.mesh.Boundary(i).GetHybridType();
	    Ufine.Reallocate(nb_dof_fine); Ucoarse.Reallocate(nb_dof_coarse);
	    for (int j = 0; j < nb_dof_coarse; j++)
	      Ucoarse(j) = r0_coarse(offset_coarse + j);

	    this->proj_coarse_to_fine.ProjectScalar(Ucoarse, Ufine, type_elt);
	    
	    for (int j = 0; j < nb_dof_fine; j++)
	      r0(offset_fine + j) = Ufine(j);

	    offset_fine += nb_dof_fine;
	    offset_coarse += nb_dof_coarse;	    
	  }
      }
    else
      {
	Vector<T0> r(var_fine.offset_dof_unknown(var_fine.nb_unknowns_scal));
	Vector<T0> r_coarse(var_coarse->offset_dof_unknown(var_fine.nb_unknowns_scal));
	r.Zero(); r_coarse.Zero();

	if ((this->fine_solver->StaticCondensation()) && (!var_fine.GetComputationProblem().LightStaticCondensation()))
	  this->coarse_solver->GetCondensedSolver(T(0)).ScatterVector(r0_coarse, r_coarse);
	else
	  r_coarse = r0_coarse;

	int nb_unknowns = var_coarse->nb_unknowns;
	Vector<Vector<T0> > Ucoarse(nb_unknowns), Vfine(nb_unknowns);
	const MeshNumbering<Dimension>& mesh_coarse = var_coarse->GetMeshNumbering(0);
	const MeshNumbering<Dimension>& mesh_fine = var_fine.GetMeshNumbering(0);
	int Nvol_coarse = mesh_coarse.GetNbDof();
	int Nvol_fine = mesh_fine.GetNbDof();
	
        for (int e = 0; e < var_fine.mesh.GetNbElt(); e++)
	  {
	    int nb_dof_coarse = mesh_coarse.GetNbLocalDof(e);
	    int nb_dof_fine = mesh_fine.GetNbLocalDof(e);
	    for (int m = 0; m < nb_unknowns; m++) 
	      {
		Ucoarse(m).Reallocate(nb_dof_coarse);
		Vfine(m).Reallocate(nb_dof_fine);
	      }
	    
	    // we get local variable Ucoarse from global variable r_coarse
	    for (int i = 0; i < nb_dof_coarse; i++)
	      {
		int num_dof = mesh_coarse.Element(e).GetNumberDof(i);
		for (int m = 0; m < nb_unknowns; m++) 
		  Ucoarse(m)(i) = r_coarse(num_dof+m*Nvol_coarse);
	      }
            
            mesh_coarse.number_map.ModifyLocalUnknownVector(mesh_coarse, Ucoarse, e);
            
	    // local prolongation
	    int type_elt = var_fine.mesh.GetTypeElement(e);
	    for (int m = 0; m < nb_unknowns; m++) 
	      this->proj_coarse_to_fine.ProjectScalar(Ucoarse(m), Vfine(m), type_elt);
	    
            mesh_fine.number_map.ModifyLocalUnknownVector(mesh_fine, Vfine, e);
            
	    // we put local variable Vfine to global variable r
	    for (int i = 0; i < nb_dof_fine; i++)
	      {
		int num_dof = mesh_fine.Element(e).GetNumberDof(i);
		for (int m = 0; m < nb_unknowns; m++) 
		  r(num_dof+m*Nvol_fine) += Vfine(m)(i);
	      }
	  }
	
	DistributedProblem_Base& var_comm = dynamic_cast<DistributedProblem_Base&>(var_fine);
	var_comm.AddDomains(r);
	
	for (int i = 0; i < Nvol_fine; i++)
	  for (int m = 0; m < nb_unknowns; m++) 
	    r(i+m*Nvol_fine) *= this->InvDh(i);

	if (this->fine_solver->StaticCondensation() && (!var_fine.GetComputationProblem().LightStaticCondensation()))
	  this->cond_solver->ExtractVector(r, r0);
	else
	  r0 = r;
      }
    
    var_fine.GetBoundaryConditionProblem().ImposeNullDirichletCondition(r0);
  }


  template<class T, class Dimension>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::GetRestriction(const Vector<Treal>& r, Vector<Treal>& r_coarse) const
  {
    GetRestrictionGen(r, r_coarse);
  }


  template<class T, class Dimension>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::GetRestriction(const Vector<Tcplx>& r, Vector<Tcplx>& r_coarse) const
  {
    GetRestrictionGen(r, r_coarse);
  }
    

  template<class T, class Dimension>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::GetProlongation(const Vector<Treal>& r_coarse, Vector<Treal>& r) const
  {
    GetProlongationGen(r_coarse, r);
  }


  template<class T, class Dimension>
  void MultigridPreconditioning_Dim<T, Dimension>
  ::GetProlongation(const Vector<Tcplx>& r_coarse, Vector<Tcplx>& r) const
  {
    GetProlongationGen(r_coarse, r);
  }

    
  // MultigridPreconditioning //
  //////////////////////////////
  
  
  ////////////////////////
  // All_Preconditioner //

  
  //! destructor
  template<class T>
  All_Preconditioner_Base<T>::~All_Preconditioner_Base()
  {
    delete prec_subdomain;
    delete prec_multigrid;
    delete prec_low_order;
    delete prec_local;
#ifdef SELDON_WITH_HYPRE
    delete prec_hypre;
#endif
  }


  //! Sets values of attributes to default values
  template<class T>
  void All_Preconditioner_Base<T>::InitDefaultValues()
  {
#ifdef SELDON_WITH_HYPRE
    prec_hypre = new HyprePreconditioner<T>();
#endif    

    type_precond = IDENTITY;
    threshold = 0;    
    SetComplexOne(alpha_laplace);
    mat_ssor_sym = NULL;
    mat_ssor_unsym = NULL;
    matCSR_ssor_sym = NULL;
    matCSR_ssor_unsym = NULL;
    mat_jacobi = NULL;
    mat_stored_jacobi = NULL;
    local_precond = false;
  }
  
  
  //! additional parameters
  template<class T>
  void All_Preconditioner_Base<T>
  ::SetInputData(const string& keyword, const VectString& param)
  {
    prec_low_order->SetInputData(keyword, param);
    prec_multigrid->SetInputData(keyword, param);
    prec_subdomain->SetInputData(keyword, param);
    
    if (keyword == "DampingParameters")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of All_Preconditioner" << endl;
	    cout << "DampingParameters needs more parameters, for instance :" << endl;
	    cout << "DampingParameters = alpha" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	alpha_laplace = to_num<T>(param(0));
      }
    else if (keyword == "LocalPreconditioning")
      {
	if (param(0) == "YES")
	  {
	    local_precond = true;
	    prec_local->SetInputPreconditioning(keyword, param);
	  }
	else
	  local_precond = false;
      }
  }
  
  
  //! parameters present in the line TypeResolution = GMRES JACOBI parameters
  template<class T>
  void All_Preconditioner_Base<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param)
  {
    if ((!keyword.compare("JACOBI"))||(!keyword.compare("BLOCK_JACOBI")))
      {
	prec_jacobi.SetInputPreconditioning(keyword, param);
	type_precond = JACOBI;
      }
    else if (!keyword.compare("SOR"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputPreconditioning of All_Preconditioner" << endl;
	    cout << "SOR omega needs more parameters, for instance :" << endl;
	    cout << "SOR omega" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	type_precond = SSOR;
	prec_ssor.SetParameterRelaxation(to_num<Real_wp>(param(0)));
	if (param.GetM() > 1)
	  prec_ssor.SetNumberIterations(to_num<int>(param(1)));
	
	// unsymmetric
	prec_ssor.InitUnSymmetricPreconditioning();
      }
    else if (!keyword.compare("SSOR"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputPreconditioning of All_Preconditioner" << endl;
	    cout << "SSOR omega needs more parameters, for instance :" << endl;
	    cout << "SSOR omega" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	type_precond = SSOR;
	prec_ssor.SetParameterRelaxation(to_num<Real_wp>(param(0)));
	if (param.GetM() > 1)
	  prec_ssor.SetNumberIterations(to_num<int>(param(1)));
	
	// symmetric
	prec_ssor.InitSymmetricPreconditioning();
      }
    else if (!keyword.compare("ILUT"))
      {
	VectString parameters(1); parameters(0) = string("ILUT");
	parameters.PushBack(param);
	mat_lu.SelectDirectSolver(parameters);
	type_precond = DIRECT;
      }
    else if (!keyword.compare("DIRECT"))
      {
	Dimension2 dim;
	mat_lu.SelectOptimalDirectSolver(10000, T(1), dim);
	type_precond = DIRECT;
      }
    else if (!keyword.compare("LOW_ORDER"))
      {
	prec_low_order->SetInputPreconditioning(keyword, param);
	type_precond = LOW_ORDER;
      }
    else if (!keyword.compare("MULTIGRID"))
      {
	prec_multigrid->SetInputPreconditioning(keyword, param);
	type_precond = MULTIGRID;
      }
    else if (!keyword.compare("SUBDOMAIN"))
      {
	prec_subdomain->SetInputPreconditioning(keyword, param);
	type_precond = SUBDOMAIN;
      }
    else if (!keyword.compare("IDENTITY"))
      {
	type_precond = IDENTITY;
      }
    else if (keyword == "HYPRE")
      {
	VectString other_param(param.GetM()-1);
	for (int i = 1; i < param.GetM(); i++)
	  other_param(i-1) = param(i);
	
#ifdef SELDON_WITH_HYPRE
	prec_hypre->SetInputPreconditioning(param(0), other_param);
	type_precond = HYPRE;
#else
        cout << "recompile Montjoie with Hypre" << endl;
        abort();
#endif
      }
    else
      {
        cout << "Unknown preconditioning " << endl;
        cout << keyword << endl;
        abort();
      }
  }
  
  
  //! sets damping coefficient alpha
  template<class T>
  void All_Preconditioner_Base<T>::SetDampingCoefficient(const T& alpha)
  {
    alpha_laplace = alpha;
    prec_multigrid->SetDampingCoefficient(alpha_laplace);
  }
  
  
  //! applies preconditioning, x = M^{-1} b
  template<class T>
  void All_Preconditioner_Base<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonNoTrans, A, b, x);
  }
  
  
  //! applies transpose preconditioning, x = M^{-T} b
  template<class T>
  void All_Preconditioner_Base<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonTrans, A, b, x);
  }


  //! applies preconditioning, x = M^{-1} b (or its transpose
  template<class T> template<class T0>
  void All_Preconditioner_Base<T>::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
					    const Vector<T0>& b, Vector<T0>& x)
  {
    switch(type_precond)
      {
      case SSOR :
	if (mat_ssor_sym != NULL)
          prec_ssor.SolveGen(trans, *mat_ssor_sym, b, x);
	else if (matCSR_ssor_sym != NULL)
          prec_ssor.SolveGen(trans, *matCSR_ssor_sym, b, x);
	else if (mat_ssor_unsym != NULL)
          prec_ssor.SolveGen(trans, *mat_ssor_unsym, b, x);
        else
          prec_ssor.SolveGen(trans, *matCSR_ssor_unsym, b, x);
	break;
      case DIRECT :
	Copy(b, x);
	mat_lu.Solve(trans, x);
	break;
      case JACOBI :
	prec_jacobi.SolveGen(trans, *mat_jacobi, b, x, true);
	break;
      case LOW_ORDER :
	prec_low_order->Solve(trans, A, b, x);
	break;
      case MULTIGRID :
	prec_multigrid->SolveGen(trans, A, b, x);
	break;
      case SUBDOMAIN :
	prec_subdomain->SolveGen(trans, A, b, x);
	break;
      case HYPRE :
	SolveHypreGen(trans, A, b, x);
	break;
      case LOCAL :
	prec_local->SolveGen(trans, A, b, x);
	break;
      default :
	Copy(b, x);
      }
  }

  
  template<>
  void All_Preconditioner_Base<Complex_wp>
  ::SolveHypreGen(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&,
		  const Vector<Real_wp>& b, Vector<Real_wp>& x)
  {
    cout << "Incompatibles types" << endl;
    abort();
  }

  
  template<>
  void All_Preconditioner_Base<Complex_wp>
  ::SolveHypreGen(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&,
		  const Vector<Complex_wp>& b, Vector<Complex_wp>& x)
  {
    cout << "Hypre compiled with real numbers" << endl;
    abort();
  }


  template<>
  void All_Preconditioner_Base<Real_wp>
  ::SolveHypreGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&,
		  const Vector<Complex_wp>& b, Vector<Complex_wp>& x)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<>
  void All_Preconditioner_Base<Real_wp>
  ::SolveHypreGen(const SeldonTranspose& trans, const VirtualMatrix<Real_wp>& A,
		  const Vector<Real_wp>& b, Vector<Real_wp>& x)
  {
#ifdef SELDON_WITH_HYPRE
    static_cast<HyprePreconditioner<Real_wp>* >(prec_hypre)->Solve(trans, A, b, x);
#else
    cout << "Recompile Montjoie with Hypre" << endl;
    abort();
#endif
  }
  
  
  template<>
  void All_Preconditioner_Base<Complex_wp>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<Complex_wp>& A,
	  const Vector<Real_wp>& b, Vector<Real_wp>& x)
  {
    cout << "Incompatibles types" << endl;
    abort();
  }
  

  template<class T>
  void All_Preconditioner_Base<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x)
  {
    SolveGen(trans, A, b, x);
  }


  template<class T>
  void All_Preconditioner_Base<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x)
  {
    SolveGen(trans, A, b, x);
  }

  
  //! returns memory used by the object in bytes
  template<class T>
  size_t All_Preconditioner_Base<T>::GetMemorySize() const
  {
    size_t taille = mat_lu.GetMemorySize() + sizeof(*this);
    taille += mat_stored_sym.GetMemorySize() + mat_stored_unsym.GetMemorySize();
    taille += prec_jacobi.GetMemorySize();
    if (mat_stored_jacobi != NULL)
      taille += mat_stored_jacobi->GetMemorySize();

    taille += prec_subdomain->GetMemorySize() + prec_multigrid->GetMemorySize();
    taille += prec_low_order->GetMemorySize();
    taille += prec_local->GetMemorySize();
    return taille;
  }
  
  
  //! sets maximal size of linear system for different preconditionings
  template<class T> template<class Dimension>
  void All_Preconditioner_Base<T>::
  GetMaximumSize(T a, Dimension dim, int& nb_max_ilut,
		 int& nb_max_iterative, int& nb_max_multigrid,
		 int& order_min_multigrid)
  {
  }
  
  
  //! sets optimal preconditioning parameters
  template<class T>
  void All_Preconditioner_Base<T>
  ::SetOptimalParameters(int order, int N, All_LinearSolver& solver)
  {
  }
  
  
  //! constructs preconditioning
  template<class T>
  void All_Preconditioner_Base<T>
  ::ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
			    const FemMatrixFreeClass_Base<T>& A, const GlobalGenericMatrix<T>& nature,
			    CondensationBlockSolver_Fem<T>& cond_solver)
  {
    GlobalGenericMatrix<T> nat_mat;
    FemMatrixFreeClass_Base<T>& Avar = const_cast<FemMatrixFreeClass_Base<T>&>(A);
    nat_mat.SetCoefMass(nature.GetCoefMass()*alpha_laplace);
    nat_mat.SetCoefStiffness(nature.GetCoefStiffness());
    nat_mat.SetCoefDamping(nature.GetCoefDamping());

    if (local_precond)
      {
	prec_local->ConstructPreconditioner(var, solver, A, nat_mat, cond_solver);
	type_precond = LOCAL;
	return;
      }
    
    if (var.GetPrintLevel() >= 4)
      mat_lu.SetPrintLevel(2);
    
    if (solver.StaticCondensation())
      var.SetLeafStaticCondensation(true);    

    CondensationBlockSolver_Fem<T>* other_cond
      = static_cast<CondensationBlockSolver_Fem<T>* >(var.GetNewCondensationSolver(T(0)));

    other_cond->CopyParameter(cond_solver);
    T one; SetComplexOne(one);
    switch (type_precond)
      {
      case SSOR :
	mat_stored_sym.Clear();
        mat_stored_unsym.Clear();
	if (IsSymmetricMatrix(A))
	  {
            prec_ssor.InitSymmetricPreconditioning();
            bool compute_mat = true;
            if (alpha_laplace == one)
              compute_mat = !Avar.SucceedInAffectingPointer(mat_ssor_sym, matCSR_ssor_sym);
            
            if (compute_mat)
              {
                var.AddMatrixWithBC(mat_stored_sym, nat_mat, 0, 0, other_cond);
		if (solver.StaticCondensation())
		  cond_solver.CompressMatrixGen(mat_stored_sym);
		
                mat_ssor_sym = &mat_stored_sym;
              }
	    
	    mat_ssor_unsym = NULL;
	    matCSR_ssor_unsym = NULL;
          }
	else
	  {
            prec_ssor.InitUnSymmetricPreconditioning();
            bool compute_mat = true;
            if (alpha_laplace == one)
              compute_mat = !Avar.SucceedInAffectingPointer(mat_ssor_unsym, matCSR_ssor_unsym);
            
            if (compute_mat)
              {
                var.AddMatrixWithBC(mat_stored_unsym, nat_mat, 0, 0, other_cond);
		if (solver.StaticCondensation())
		  cond_solver.CompressMatrixGen(mat_stored_unsym);

                mat_ssor_unsym = &mat_stored_unsym;
              }
	    
	    mat_ssor_sym = NULL;
	    matCSR_ssor_sym = NULL;
          }
	
	break;	
        
      case DIRECT :
	if (IsSymmetricMatrix(A))
          {
            DistributedMatrix<T, Symmetric, ArrayRowSymSparse> mat_sym;
            var.AddMatrixWithBC(mat_sym, nat_mat, 0, 0, other_cond);
	    if (solver.StaticCondensation())
	      cond_solver.CompressMatrixGen(mat_sym);

            //mat_sym.WriteText("mat_sym.dat");
	    mat_lu.Factorize(mat_sym);
          }
        else
          {
            DistributedMatrix<T, General, ArrayRowSparse> mat_unsym;
            var.AddMatrixWithBC(mat_unsym, nat_mat, 0, 0, other_cond);
	    if (solver.StaticCondensation())
	      cond_solver.CompressMatrixGen(mat_unsym);

            mat_lu.Factorize(mat_unsym);
          }
        
	break;	
        
      case JACOBI :
	if (alpha_laplace == one)
	  mat_jacobi = &Avar;
        else
          {
	    if (mat_stored_jacobi == NULL)
	      mat_stored_jacobi = var.GetNewIterativeMatrix(T(0));
	    
            var.AddMatrixWithBC(*mat_stored_jacobi, nat_mat, 0, 0, other_cond);
	    if (solver.StaticCondensation())
	      cond_solver.CompressMatrixGen(*mat_stored_jacobi);
	    
            mat_jacobi = mat_stored_jacobi;
          }
        
	prec_jacobi.ConstructPreconditioner(var, solver, *mat_jacobi, nat_mat, cond_solver);
	break;
	
      case LOW_ORDER :
	prec_low_order->ConstructPreconditioner(var, solver, A, nat_mat, cond_solver);
	break;
	
      case MULTIGRID :
	prec_multigrid->ConstructPreconditioner(var, solver, A, nat_mat, cond_solver);
	break;
	
      case SUBDOMAIN :
	prec_subdomain->ConstructPreconditioner(var, solver, A, nat_mat, cond_solver);
	break;

      case HYPRE :
#ifdef SELDON_WITH_HYPRE
        if (IsSymmetricMatrix(A))
          {
            DistributedMatrix<T, Symmetric, ArrayRowSymSparse> mat_sym;
            var.AddMatrixWithBC(mat_sym, nat_mat, 0, 0, other_cond);
            if (solver.StaticCondensation())
              cond_solver.CompressMatrixGen(mat_sym);
            
            static_cast<HyprePreconditioner<T>* >(prec_hypre)->ConstructPreconditioner(mat_sym, false);
          }
        else
          {
            DistributedMatrix<T, General, ArrayRowSparse> mat_unsym;
            var.AddMatrixWithBC(mat_unsym, nat_mat, 0, 0, other_cond);
            if (solver.StaticCondensation())
              cond_solver.CompressMatrixGen(mat_unsym);
            
            static_cast<HyprePreconditioner<T>* >(prec_hypre)->ConstructPreconditioner(mat_unsym, false);
          }
#else
        cout << "Recompile with Hypre" << endl;
        abort();
#endif
	break;	

      default:
	break;
      }

    delete other_cond;
    
    if (solver.StaticCondensation())
      var.SetLeafStaticCondensation(false);    

    mat_lu.SetPrintLevel(0);
  }


  template<class T>
  void All_Preconditioner_Base<T>::CopyParameter(const Preconditioner_Base<T>& precond)
  {
    const All_Preconditioner_Base<T>& prec = static_cast<const All_Preconditioner_Base<T>& >(precond);
    type_precond = prec.type_precond;
    threshold = prec.threshold;
    prec_jacobi.CopyParameter(prec.prec_jacobi);
    prec_ssor.CopyParameter(prec.prec_ssor);
    prec_subdomain->CopyParameter(*prec.prec_subdomain);
    prec_multigrid->CopyParameter(*prec.prec_multigrid);
    prec_low_order->CopyParameter(*prec.prec_low_order);
    alpha_laplace = prec.alpha_laplace;
    mat_lu.CopyParameter(prec.mat_lu);
    local_precond = prec.local_precond;
    prec_local->CopyParameter(*prec.prec_local);
#ifdef SELDON_WITH_HYPRE
    prec_hypre->CopyParameter(*prec.prec_hypre);
#endif
  }
  
  
  // All_Preconditioner //
  ////////////////////////


  /************************
   * LocalPreconditioning *
   ************************/


  //! default constructor
  template<class T>
  LocalPreconditioning_Base<T>::LocalPreconditioning_Base()
  {
    precond = NULL;
    iterative_matrix = NULL;
    local_cond = NULL;
  }


  //! destructor
  template<class T>
  LocalPreconditioning_Base<T>::~LocalPreconditioning_Base()
  {
    if (precond != NULL)
      delete precond;

    if (iterative_matrix != NULL)
      delete iterative_matrix;

    if (local_cond != NULL)
      delete local_cond;
  }


  //! modifies parameters of the object with a line of the data file
  template<class T>
  void LocalPreconditioning_Base<T>::SetInputPreconditioning(const string& keyword, const Vector<string>& param)
  {
    if (param(1) == "LOCAL_DT")
      {
	VectReal_wp local_dt;
	local_dt.ReadText(param(2));

	Real_wp dt_seuil = to_num<Real_wp>(param(3));
	int nb_elt = 0;
	for (int i = 0; i < local_dt.GetM(); i++)
	  if (local_dt(i) < dt_seuil)
	    nb_elt++;

	num_elt_precond.Reallocate(nb_elt);
	nb_elt = 0;
	for (int i = 0; i < local_dt.GetM(); i++)
	  if (local_dt(i) < dt_seuil)
	    num_elt_precond(nb_elt++) = i;

      }
    else if (param(1) == "ELEMENT")
      num_elt_precond.ReadText(param(2));
  }
  

  //! copies the parameters of prec
  template<class T>
  void LocalPreconditioning_Base<T>::CopyParameter(const LocalPreconditioning_Base<T>& prec)
  {
    num_elt_precond = prec.num_elt_precond;
  }


  //! returns the memory used by the object in bytes
  template<class T>
  size_t LocalPreconditioning_Base<T>::GetMemorySize() const
  {
    //return precond->GetMemorySize();
    return 0;
  }
  
  
  //! applies preconditioning
  template<class T> template<class T0>
  void LocalPreconditioning_Base<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	     const Vector<T0>& b, Vector<T0>& x)
  {
    Vector<T0> b_loc(ListeDof.GetM()), x_loc(ListeDof.GetM());
    x_loc.Zero();
    for (int i = 0; i < ListeDof.GetM(); i++)
      b_loc(i) = b(ListeDof(i));
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = b(i)*diagonal(i);

    precond->Solve(trans, *iterative_matrix, b_loc, x_loc);
    
    for (int i = 0; i < ListeDof.GetM(); i++)
      x(ListeDof(i)) = x_loc(i);
        

  }

  
  template<class T>
  void LocalPreconditioning_Base<T>::FinalizePreconditioner(VarComputationProblem& var,
							    VarComputationProblem& var_global, All_LinearSolver& solver,
							    const FemMatrixFreeClass_Base<T>& A,
							    const GlobalGenericMatrix<T>& nature,
							    CondensationBlockSolver_Fem<T>& cond_solver)
  {
    this->precond = var.GetNewPreconditioning(T(0));
    this->precond->CopyParameter(solver.GetPreconditioning(T(0)));
    this->precond->DisableLocalPreconditioning();
    solver.DisableLocalPreconditioning();

    this->iterative_matrix = var.GetNewIterativeMatrix(T(0));

    if (solver.StaticCondensation())
      var.SetLeafStaticCondensation(true);
    
    local_cond = static_cast<CondensationBlockSolver_Fem<T>* >(var.GetNewCondensationSolver(T(0)));
    local_cond->CopyParameter(cond_solver);
 
    var.AddMatrixWithBC(*this->iterative_matrix, nature, 0, 0, local_cond);

    if (solver.StaticCondensation())
      {
	if (solver.UseCompressedMatrix())
	  this->iterative_matrix->SetCondensedSolver(local_cond);
	else
	  this->iterative_matrix->SetCondensedSolver(NULL);
      }
    
    // this->iterative_matrix->mat_iterative_sym.WriteText("mat_local.dat");

    if (solver.StaticCondensation())
      local_cond->Init(nature, solver.UseCompressedMatrix());

    if (solver.UseCompressedMatrix())
      local_cond->CompressMatrixGen(*this->iterative_matrix);
    
    this->precond->ConstructPreconditioner(var, solver, *this->iterative_matrix, nature, *local_cond);

    Vector<int> NewListeDof(local_cond->GetNbRows());
    const Vector<int>& IndexLocal = local_cond->GetIndexCondensedRow();
    const Vector<int>& IndexGlobal = cond_solver.GetIndexCondensedRow();
    for (int i = 0; i < ListeDof.GetM(); i++)
      if (IndexLocal(i) >= 0)
	NewListeDof(IndexLocal(i)) = IndexGlobal(ListeDof(i));

    ListeDof = NewListeDof;

    if (solver.StaticCondensation())
      var.SetLeafStaticCondensation(false);
	
    var_global.ComputeDiagonalMatrix(this->diagonal, A, nature);
    for (int i = 0; i < this->diagonal.GetM(); i++)
      this->diagonal(i) = Real_wp(1) / this->diagonal(i);
    
    
  }

  
  /****************************
   * LocalPreconditioning_Dim *
   ****************************/
  
  
  //! default constructor
  template<class T, class Dimension>
  LocalPreconditioning_Dim<T, Dimension>::LocalPreconditioning_Dim()
  {
    var_local = NULL;
  }


  //! destructor
  template<class T, class Dimension>
  LocalPreconditioning_Dim<T, Dimension>::~LocalPreconditioning_Dim()
  {
    if (var_local != NULL)
      delete var_local;
    
    var_local = NULL;
  }


  //! constructs local preconditioning
  template<class T, class Dimension>
  void LocalPreconditioning_Dim<T, Dimension>::
  ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
			  const FemMatrixFreeClass_Base<T>& A,
			  const GlobalGenericMatrix<T>& nature,
			  CondensationBlockSolver_Fem<T>& cond_solver)
  {
    if (var_local != NULL)
      delete var_local;

    // we construct a new elliptic problem
    DistributedProblem<Dimension>& var_global = dynamic_cast<DistributedProblem<Dimension>& >(var);
    var_local = var_global.GetNewEllipticProblem();
    var_local->CopyInputData(var_global);

    int nb_elt_sub = this->num_elt_precond.GetM();
    int nb_vertex_sub = 0;
    Vector<bool> ElementOnSubdomain(var_global.mesh.GetNbElt());
    Vector<bool> VertexOnSubdomain(var_global.mesh.GetNbVertices());
    VertexOnSubdomain.Fill(false);
    ElementOnSubdomain.Fill(false);
    for (int i = 0; i < nb_elt_sub; i++)
      {
	int ne = this->num_elt_precond(i);
	ElementOnSubdomain(ne) = true;
	int nb_vert = var_global.mesh.Element(ne).GetNbVertices();
	for (int j = 0; j < nb_vert; j++)
	  {
	    int nv = var_global.mesh.Element(ne).numVertex(j);
	    if (!VertexOnSubdomain(nv))
	      {
		nb_vertex_sub++;
		VertexOnSubdomain(nv) = true;
	      }
	  }
      }

    Vector<int> RefUsed(var_global.mesh.GetNbReferences()+1);
    RefUsed.Fill(-1);

    for (int i = 0; i < var_global.mesh.GetNbBoundaryRef(); i++)
      RefUsed(var_global.mesh.BoundaryRef(i).GetReference()) = 0;
    
    // we change the mesh
    var_global.mesh.CreateSubmesh(var_local->mesh, nb_vertex_sub, nb_elt_sub,
				  VertexOnSubdomain, ElementOnSubdomain);

    var_local->mesh.AddBoundaryEdges();
    var_local->mesh.ProjectPointsOnCurves();
    
    for (int i = 0; i < var_local->mesh.GetNbBoundaryRef(); i++)
      if (RefUsed(var_local->mesh.BoundaryRef(i).GetReference()) == -1)
	RefUsed(var_local->mesh.BoundaryRef(i).GetReference()) = 1;

    for (int i = 0; i < RefUsed.GetM(); i++)
      if (RefUsed(i) == 1)
	var_local->mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_ABSORBING);    

    var_local->mesh.Write("sub.mesh");
    var_local->mesh_data.Clear();
    var_local->mesh.ClearInput();

    // we construct the local problem
    var_local->DoNotComputeGrid();
    var_local->ComputeMeshAndFiniteElement(var_global.finite_element_name);
    var_local->PerformOtherInitializations();
    var_local->ComputeMassMatrix();

    Vector<string> param(1);
    param(0) = string("sub.mesh");
    var_local->mesh_data.PushBack(param);
    var_local->mesh.SetPathName("./");
    
    // index arrays between local and global problem
    Vector<int> IndexDof(var_global.GetMeshNumbering(0).GetNbDof());
    IndexDof.Fill(-1);
    this->ListeDof.Reallocate(var_local->GetMeshNumbering(0).GetNbDof());
    for (int i = 0; i < nb_elt_sub; i++)
      {
	int ne = this->num_elt_precond(i);
	IVect num_dof = var_global.GetMeshNumbering(0).Element(ne).GetNodle();
	IVect num_dof_loc = var_local->GetMeshNumbering(0).Element(i).GetNodle();

	for (int j = 0; j < num_dof.GetM(); j++)
	  if (IndexDof(num_dof(j)) == -1)
	    {
	      IndexDof(num_dof(j)) = num_dof_loc(j);
	      this->ListeDof(num_dof_loc(j)) = num_dof(j);
	    }
      }

    this->FinalizePreconditioner(var_local->GetComputationProblem(), var, solver, A, nature, cond_solver);
  }
      
#else
  

  //! other parameters
  template<class T>
  void All_Preconditioner_Base<T>
  ::SetInputData(const string& keyword, const VectString& param)
  {
  }
  

  //! parameters present in the line TypeResolution = GMRES IDENTITY parameters
  template<class T>
  void All_Preconditioner_Base<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param)
  {
    if (!keyword.compare("IDENTITY"))
      {
        
      }
    else
      {
        cout << "Recompile with SELDON_WITH_PRECONDITIONING" << endl;
        abort();
      }
  }
  
#endif

} // namespace Montjoie

#define MONTJOIE_FILE_PRECONDITIONER_CXX
#endif

