#ifndef MONTJOIE_FILE_SOLVE_SYSTEM_CXX

namespace Montjoie
{
  
  /****************
   * All_MatrixLU *
   ****************/
  
  
  //! returns an available direct solver
  /*!
    The returned string can be used when calling SelectDirectSolver.
  */
  template <class T>
  string All_MatrixLU<T>::GetAvailableSolver()
  {

#ifdef SELDON_WITH_MUMPS
    return string("MUMPS");
#endif

#ifdef SELDON_WITH_PASTIX
    return string("PASTIX");
#endif

#ifdef SELDON_WITH_WSMP
    return string("WSMP");
#endif

#ifdef SELDON_WITH_PARDISO
    return string("PARDISO");
#endif

#ifdef SELDON_WITH_UMFPACK
    return string("UMFPACK");
#endif

#ifdef SELDON_WITH_SUPERLU
    return string("SUPERLU");
#endif
    
    return string("SELDON");
  }

    
  //! returns the best direct solver available
  template<class T> template<class Dimension>
  void All_MatrixLU<T>::SelectOptimalDirectSolver(int N, T a, Dimension& dim)
  {
    SelectOptimalDirectSolver();
  }
  
  
  //! returns the best direct solver available
  template<class T>
  void All_MatrixLU<T>::SelectOptimalDirectSolver()
  {
    this->type_solver = -1;
    
#ifdef SELDON_WITH_SUPERLU
    this->type_solver = this->SUPERLU;
#endif
#ifdef SELDON_WITH_UMFPACK
    this->type_solver = this->UMFPACK;
#endif
#ifdef SELDON_WITH_PARDISO
    this->type_solver = this->PARDISO;
#endif
#ifdef SELDON_WITH_WSMP
    this->type_solver = this->WSMP;
#endif
#ifdef SELDON_WITH_PASTIX
    this->type_solver = this->PASTIX;
#endif
#ifdef SELDON_WITH_MUMPS
    this->type_solver = this->MUMPS;
#endif
    
    // Default Seldon direct solver is not proposed because too slow
    // if no solver is found, the program is stopped
    if (this->type_solver == -1)
      {
	cout << "No available direct solver " << endl;
	abort();
      }

    this->InitSolver();
  }


  template<class T>
  void All_MatrixLU<T>::SelectDirectSolver(int type)
  {
    SparseDirectSolver<T>::SelectDirectSolver(type);
  }
   

  //! returns the best direct solver available
  //! sets direct solver if required
  /*!
    returns true if the asked resolution is iterative, 
    otherwise type_solver is changed depending on parameters
   */
  template<class T>
  bool All_MatrixLU<T>::SelectDirectSolver(const VectString& parameters)
  {
    if (parameters.GetM() < 1)
      {
        this->SelectOptimalDirectSolver();
        return false;
      }
    
    if (!parameters(0).compare("UMFPACK")) 
      this->SelectDirectSolver(this->UMFPACK);
    else if (!parameters(0).compare("SUPERLU"))
      this->SelectDirectSolver(this->SUPERLU);
    else if (!parameters(0).compare("MUMPS"))
      this->SelectDirectSolver(this->MUMPS);
    else if (!parameters(0).compare("PARDISO"))
      this->SelectDirectSolver(this->PARDISO);
    else if (!parameters(0).compare("PASTIX"))
      this->SelectDirectSolver(this->PASTIX);
    else if (!parameters(0).compare("WSMP"))
      this->SelectDirectSolver(this->WSMP);
    else if (!parameters(0).compare("SELDON"))
      this->SelectDirectSolver(this->SELDON_SOLVER);
    else if (!parameters(0).compare("ILUT"))
      {
#ifdef SELDON_WITH_PRECONDITIONING
	this->SelectDirectSolver(this->ILUT);
	this->threshold_matrix = to_num<double>(parameters(1));
        int type_ilu = IlutPreconditioning<T>::ILUT;
        if (!parameters(2).compare("ILU_0"))
          type_ilu = IlutPreconditioning<T>::ILU_0;
        else if (!parameters(2).compare("ILU_D"))
          type_ilu = IlutPreconditioning<T>::ILU_D;
        else if (!parameters(2).compare("ILU_K"))
          type_ilu = IlutPreconditioning<T>::ILU_K;
        else if (!parameters(2).compare("ILUT_K"))
          type_ilu = IlutPreconditioning<T>::ILUT_K;
        
	IlutPreconditioning<T>& mat_ilut =
	  dynamic_cast<IlutPreconditioning<T>& >(*this->solver);

        mat_ilut.SetFactorisationType(type_ilu);
        mat_ilut.SetFillLevel(to_num<int>(parameters(3)));
        mat_ilut.SetDroppingThreshold(to_num<double>(parameters(4)));
        mat_ilut.SetPivotThreshold(to_num<double>(parameters(5)));
        mat_ilut.SetDiagonalCoefficient(to_num<double>(parameters(6)));

        return false;
#else
	cout<<"Montjoie compiled without preconditioners "<<endl;
	abort();
#endif
      }
    else if (parameters(0) == "DIRECT")
      this->SelectOptimalDirectSolver();
    else
      return true;
  
    if (parameters.GetM() > 1)
      this->SelectOrdering(parameters(1));
    
    return false;
  }

    
  //! sets the type of ordering to use
  template<class T>
  void All_MatrixLU<T>::SelectOrdering(const string& keyword)
  {
    int type = SparseMatrixOrdering::AUTO;
    if (keyword == "Identity")
      type = SparseMatrixOrdering::IDENTITY;
    else if (keyword == "ReverseCuthillMcKee")
      type = SparseMatrixOrdering::REVERSE_CUTHILL_MCKEE;
    else if (keyword == "Pord")
      type = SparseMatrixOrdering::PORD;
    else if (keyword == "Scotch")
      type = SparseMatrixOrdering::SCOTCH;    
    else if (keyword == "Metis")
      type = SparseMatrixOrdering::METIS;
    else if (keyword == "AMD")
      type = SparseMatrixOrdering::AMD;        
    else if (keyword == "ColAMD")
      type = SparseMatrixOrdering::COLAMD;
    else if (keyword == "QAMD")
      type = SparseMatrixOrdering::QAMD;
    else if (keyword == "AMF")
      type = SparseMatrixOrdering::AMF;
    else if (keyword == "ParMetis")
      type = SparseMatrixOrdering::PARMETIS;
    else if (keyword == "PtScotch")
      type = SparseMatrixOrdering::PTSCOTCH;
    else if (keyword == "MMD_AT_PLUS_A")
      type = SparseMatrixOrdering::MMD_AT_PLUS_A;
    else if (keyword == "MMD_ATA")
      type = SparseMatrixOrdering::MMD_ATA;
    else
      type = to_num<int>(keyword);
    
    SparseDistributedSolver<T>::SelectOrdering(type);
  }
  
  
  //! Returns maximum number of rows of factorizable matrix 
  /*!
    \param[in] a dummy parameter to inform that the matrix is real
    \param[in] dim dummy parameter to inform that the problem to solve is 2-D
    \param[out] nb_max_direct maximum number of rows of a factorizable matrix
    \param[out] nb_max_no_condensation maximum number of rows of a solvable matrix
    without static condensation
    \param[out] order_min_condensation minimum order of approximation for
    static condensation
   */
  template<class T>
  void All_MatrixLU<T>::GetMaximumSize(Real_wp a, Dimension2 dim,
				       int& nb_max_direct,
				       int& nb_max_no_condensation,
				       int& order_min_condensation)
  {
    nb_max_direct = 1500000;
    nb_max_no_condensation = 1500000;
    order_min_condensation = 5;
  }
  
  
  //! Returns maximum number of rows of factorizable matrix 
  /*!
    \param[in] a dummy parameter to inform that the matrix is complex
    \param[in] dim dummy parameter to inform that the problem to solve is 2-D
    \param[out] nb_max_direct maximum number of rows of a factorizable matrix
    \param[out] nb_max_no_condensation maximum number of rows of a solvable matrix
    without static condensation
    \param[out] order_min_condensation minimum order of approximation for
    static condensation
   */
  template<class T>
  void All_MatrixLU<T>::GetMaximumSize(complex<Real_wp> a, Dimension2 dim,
				       int& nb_max_direct,
				       int& nb_max_no_condensation,
				       int& order_min_condensation)
  {
    nb_max_direct = 1200000;
    nb_max_no_condensation = 1200000;
    order_min_condensation = 5;
  }
  
  
  //! Returns maximum number of rows of factorizable matrix 
  /*!
    \param[in] a dummy parameter to inform that the matrix is real
    \param[in] dim dummy parameter to inform that the problem to solve is 3-D
    \param[out] nb_max_direct maximum number of rows of a factorizable matrix
    \param[out] nb_max_no_condensation maximum number of rows of a solvable matrix
    without static condensation
    \param[out] order_min_condensation minimum order of approximation for
    static condensation
   */
  template<class T>
  void All_MatrixLU<T>::GetMaximumSize(Real_wp a, Dimension3 dim,
				       int& nb_max_direct,
				       int& nb_max_no_condensation,
				       int& order_min_condensation)
  {
    nb_max_direct = 170000;
    nb_max_no_condensation = 140000;
    order_min_condensation = 5;
  }
  
  
  //! Returns maximum number of rows of factorizable matrix 
  /*!
    \param[in] a dummy parameter to inform that the matrix is complex
    \param[in] dim dummy parameter to inform that the problem to solve is 3-D
    \param[out] nb_max_direct maximum number of rows of a factorizable matrix
    \param[out] nb_max_no_condensation maximum number of rows of a solvable matrix
    without static condensation
    \param[out] order_min_condensation minimum order of approximation for
    static condensation
   */
  template<class T>
  void All_MatrixLU<T>::GetMaximumSize(complex<Real_wp> a, Dimension3 dim,
				       int& nb_max_direct,
				       int& nb_max_no_condensation,
				       int& order_min_condensation)
  {
    nb_max_direct = 150000;
    nb_max_no_condensation = 120000;
    order_min_condensation = 5;
  }


  //! copies input parameters of an another solver
  template<class T>
  void All_MatrixLU<T>::CopyParameter(const All_MatrixLU<T>& solver)
  {
    this->diagonal_scaling_left = solver.diagonal_scaling_left;
    this->diagonal_scaling_right = solver.diagonal_scaling_right;
    this->type_ordering = solver.type_ordering;
    this->type_solver = solver.type_solver;
    this->nb_threads_per_node = solver.nb_threads_per_node;
    this->threshold_matrix = solver.threshold_matrix;
    this->pivot_threshold = solver.pivot_threshold;
    this->refine_solution = solver.refine_solution;
    this->print_level = solver.print_level;
    this->enforce_unsym_ilut = solver.enforce_unsym_ilut;

    this->InitSolver();

    if (this->type_solver == this->ILUT)
      {
	IlutPreconditioning<T>& mat_ilut =
	  dynamic_cast<IlutPreconditioning<T>& >(*this->solver);

	IlutPreconditioning<T>& mat_ilut0 =
	  dynamic_cast<IlutPreconditioning<T>& >(*solver.solver);

	mat_ilut.SetFactorisationType(mat_ilut0.GetFactorisationType());
	mat_ilut.SetFillLevel(mat_ilut0.GetFillLevel());
	mat_ilut.SetDroppingThreshold(mat_ilut0.GetDroppingThreshold());
	mat_ilut.SetPivotThreshold(toDouble(mat_ilut0.GetPivotThreshold()));
	mat_ilut.SetDiagonalCoefficient(mat_ilut0.GetDiagonalCoefficient());
      }    
  }
  

  /**************************
   * VirtualMatrixTranspose *
   **************************/
  
  
  //! Default constructor
  template<class T>
  VirtualMatrixTranspose<T>::VirtualMatrixTranspose(const VirtualMatrix<T>& A, const SeldonTranspose& t)
    : VirtualMatrix<T>(A.GetM(), A.GetN()), A_ptr(&A), trans(t)
  {
  }
  
  
  //! y = beta y + alpha A x
  template<class T>
  void VirtualMatrixTranspose<T>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
                 const Treal& beta, Vector<Treal>& y) const
  {
    if (trans.Trans())
      A_ptr->MltAddVector(alpha, SeldonTrans, x, beta, y);
    else
      {
        cout << "not implemented" << endl;
        abort();
      }
  }
  

  //! y = beta y + alpha A x  
  template<class T>
  void VirtualMatrixTranspose<T>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
                 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    if (trans.Trans())
      A_ptr->MltAddVector(alpha, SeldonTrans, x, beta, y);
    else
      {
        cout << "not implemented" << endl;
        abort();
      }
  }

  
  //! y = beta y + alpha A^T x  
  template<class T>
  void VirtualMatrixTranspose<T>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& t,
                 const Vector<Treal>& x,
                 const Treal& beta, Vector<Treal>& y) const
  {
    if (trans.Trans())
      {
        if (t.NoTrans())
          A_ptr->MltAddVector(alpha, SeldonTrans, x, beta, y);
        else if (t.Trans())
          A_ptr->MltAddVector(alpha, SeldonNoTrans, x, beta, y);
        else
          {
            cout << "not implemented" << endl;
            abort();
          }
      }
    else
      {
        cout << "not implemented" << endl;
        abort();
      }    
  }
  
  
  //! y = beta y + alpha A^T x
  template<class T>
  void VirtualMatrixTranspose<T>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& t,
                 const Vector<Tcplx>& x,
                 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    if (trans.Trans())
      {
        if (t.NoTrans())
          A_ptr->MltAddVector(alpha, SeldonTrans, x, beta, y);
        else if (t.Trans())
          A_ptr->MltAddVector(alpha, SeldonNoTrans, x, beta, y);
        else
          {
            cout << "not implemented" << endl;
            abort();
          }
      }
    else
      {
        cout << "not implemented" << endl;
        abort();
      }    
  }

  
  //! y = A x
  template<class T>
  void VirtualMatrixTranspose<T>::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    if (trans.Trans())
      A_ptr->MltVector(SeldonTrans, x, y);
    else
      {
        cout << "not implemented" << endl;
        abort();
      }
  }
  
  
  //! y = A x
  template<class T>
  void VirtualMatrixTranspose<T>::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    if (trans.Trans())
      A_ptr->MltVector(SeldonTrans, x, y);
    else
      {
        cout << "not implemented" << endl;
        abort();
      }
  }
    
  
  //! y = A^T x
  template<class T>
  void VirtualMatrixTranspose<T>::MltVector(const SeldonTranspose& t,
                                            const Vector<Treal>& x, Vector<Treal>& y) const
  {
    if (trans.Trans())
      {
        if (t.NoTrans())
          A_ptr->MltVector(SeldonTrans, x, y);
        else if (t.Trans())
          A_ptr->MltVector(SeldonNoTrans, x, y);
        else
          {
            cout << "not implemented" << endl;
            abort();
          }
      }
    else
      {
        cout << "not implemented" << endl;
        abort();
      }    
  }
  
  
  //! y = A^T x
  template<class T>
  void VirtualMatrixTranspose<T>::MltVector(const SeldonTranspose& t,
                                            const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    if (trans.Trans())
      {
        if (t.NoTrans())
          A_ptr->MltVector(SeldonTrans, x, y);
        else if (t.Trans())
          A_ptr->MltVector(SeldonNoTrans, x, y);
        else
          {
            cout << "not implemented" << endl;
            abort();
          }
      }
    else
      {
        cout << "not implemented" << endl;
        abort();
      }    
  }
  

  /***********************
   * All_IterativeSolver *
   ***********************/


  //! Default constructor
  template <class T>
  All_IterativeSolver<T>::All_IterativeSolver()
  {
    print_level = 0;
    type_solver = COCG;
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    nb_max_iterations = 1000;
    tolerance = 1e-6;
    parameter_restart = 10;
    precond = &identity_precond;
    nb_iterations = 0;
  }
  
  
  //! constructor with a given preconditioning
  template <class T> 
  All_IterativeSolver<T>::All_IterativeSolver(Preconditioner_Base<T>& prec)
  {
    print_level = 0;
    type_solver = COCG;
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    nb_max_iterations = 1000;
    tolerance = 1e-6;
    parameter_restart = 10;
    precond = &prec;
  }


  //! sets iterative solver with a line of the data file
  template<class T>
  void All_IterativeSolver<T>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("NumberMaxIterations"))
      {
	// maximal number of iterations
	nb_max_iterations = to_num<int>(parameters(0));
	if (parameters.GetM() > 1)
          parameter_restart = to_num<int>(parameters(1));
      }
    else if (!description_field.compare("Tolerance"))
      {
	// stopping criterion for iterative resolution
	tolerance = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("TypeResolution"))
      {
        bool use_precond = true;
	if (!parameters(0).compare("QMR"))
	  type_solver = QMR;
	else if (!parameters(0).compare("TFQMR"))
	  type_solver = TFQMR;
	else if (!parameters(0).compare("GMRES"))
	  type_solver = GMRES;
	else if (!parameters(0).compare("CG"))
	  type_solver = CG;
 	else if (!parameters(0).compare("CGNE"))
	  type_solver = CGNE;
	else if (!parameters(0).compare("LSQR"))
	  type_solver = LSQR;
	else if (!parameters(0).compare("CGS"))
	  type_solver = CGS;
	else if (!parameters(0).compare("BICG"))
	  type_solver = BICG;
	else if (!parameters(0).compare("BICGSTAB"))
	  type_solver = BICGSTAB;
	else if (!parameters(0).compare("BICGCR"))
	  type_solver = BICGCR;
	else if (!parameters(0).compare("MINRES"))
	  type_solver = MINRES;
	else if (!parameters(0).compare("GCR"))
	  type_solver = GCR;
	else if (!parameters(0).compare("SYMMLQ"))
	  type_solver = SYMMLQ;
	else if (!parameters(0).compare("COCG"))
	  type_solver = COCG;
	else if (!parameters(0).compare("QMR_SYM"))
	  type_solver = QMR_SYM;
	else if (!parameters(0).compare("BICGSTABL"))
	  type_solver = BICGSTABL;
	else if (!parameters(0).compare("QCGS"))
	  type_solver = QCGS;
	else if (!parameters(0).compare("JACOBI"))
	  {
            type_solver = JACOBI;
            use_precond = false;
          }
	else if (!parameters(0).compare("GAUSS_SEIDEL"))
	  {
            type_solver = GAUSS_SEIDEL_RELAX;
            omega_relax = Real_wp(1); stage_relax = 0;
            use_precond = false;
          }
	else if (!parameters(0).compare("GAUSS_SEIDEL_RELAX"))
	  {
            type_solver = GAUSS_SEIDEL_RELAX;
            omega_relax = to_num<Real_wp>(parameters(1));
            if (parameters.GetM() > 1)
              stage_relax = to_num<int>(parameters(2));
            
            use_precond = false;
          }
	else
	  {
	    // direct solver
	    return;
	  }
	
	// second field -> preconditioner (if different from identity)
	if ((parameters.GetM() > 1) && use_precond)
	  {
	    string keyword = parameters(1);
	    VectString param;
	    if (parameters.GetM() > 2)
	      {
		param.Reallocate(parameters.GetM() - 2);
		for (int i = 2; i < parameters.GetM(); i++)
		  param(i-2) = parameters(i);
	      }
            
	    precond->SetInputPreconditioning(keyword, param);
	  }
      }
  }

  
  //! Iterative resolution of A x_sol = source_rhs
  /*!
    \param[in] mat iterative matrix
    \param[out] x_sol solution of the linear system
    \param[in] source_rhs right hand side
   */
  template<class T> template<class Vector1>
  int All_IterativeSolver<T>::Solve(const VirtualMatrix<T>& mat,
				    Vector1 & x_sol, const Vector1 & source_rhs)
  {
    Iteration<typename ClassComplexType<T>::Treal> iter(nb_max_iterations, tolerance);
    if (print_level < 1)
      iter.HideMessages();
    else
      iter.ShowMessages();
    
    iter.SetSolver(type_solver, parameter_restart, 0);
    if (print_level >= 4)
      iter.ShowFullHistory();
    
    Vector1 b_source(source_rhs);
    if (diagonal_scaling_left)
      for (int i = 0; i < x_sol.GetM(); i++)
	b_source(i) *= diagonal_scale_left(i);
    
    // null init guess ?
    if (Norm2(x_sol) <= epsilon_machine*Norm2(b_source))
      {
        x_sol.Fill(0);
        iter.SetInitGuess(true);
      }
    else
      iter.SetInitGuess(false);
    
    // now we can solve with the required iterative algorithm
    int success = Solve(mat, x_sol, b_source, *precond, iter);
    nb_iterations = iter.GetNumberIteration();
    
    if (print_level >= 2)
      {
        if (success == 0)
          cout << "Algorithm has converged in " << 
            iter.GetNumberIteration() << " iterations" << endl;
        else
          cout << "Iterative algorithm failed, Error code = " << success << endl;
      }
    
    if (diagonal_scaling_right)
      for (int i = 0; i < x_sol.GetM(); i++)
	x_sol(i) *= diagonal_scale_right(i);
    
    return success;
  }
  
  
  //! Iterative resolution of A^T x_sol = source_rhs
  /*!
    \param[in] mat iterative matrix
    \param[out] x_sol solution of the linear system
    \param[in] source_rhs right hand side
   */
  template<class T> template<class Vector1>
  int All_IterativeSolver<T>::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& mat,
				    Vector1 & x_sol, const Vector1 & source_rhs)
  {
    if (trans.NoTrans())
      {
        return Solve(mat, x_sol, source_rhs);
      }
    else if (trans.ConjTrans())
      {
        cout << "not implemented" << endl;
        abort();
      }
    
    VirtualMatrixTranspose<T> A(mat, SeldonTrans);
    return Solve(A, x_sol, source_rhs);
  }


  //! Iterative resolution of A x = b
  template<class T> template <class Vector1>
  int All_IterativeSolver<T>
  ::Solve(const VirtualMatrix<T>& A,  Vector1& x, const Vector1& b,
	  Preconditioner_Base<T>& M, Iteration<typename ClassComplexType<T>::Treal> & iter)
  {
    nb_iterations = 0;
    int ierr = -1;
    switch (type_solver)
      {
      case PRECOND: M.Solve(A, b, x); return 0;
      case GMRES: ierr = Gmres(A, x, b, M, iter); break;
      case QMR: ierr = Qmr(A, x, b, M, iter); break;
      case COCG: ierr = CoCg(A, x, b, M, iter); break;
      case BICGCR: ierr = BiCgcr(A, x, b, M, iter); break;
      case BICGSTAB: ierr = BiCgStab(A, x, b, M, iter); break;
      case TFQMR: ierr = TfQmr(A, x, b, M, iter); break;
      case CG: ierr = Cg(A, x, b, M, iter); break;
      case CGNE: ierr = Cgne(A, x, b, M, iter); break;
      case LSQR: ierr = Lsqr(A, x, b, M, iter); break;
      case CGS: ierr = Cgs(A, x, b, M, iter); break;
      case BICG: ierr = BiCg(A, x, b, M, iter); break;
      case MINRES: ierr = MinRes(A, x, b, M, iter); break;
      case SYMMLQ: ierr = Symmlq(A, x, b, M, iter); break;
      case GCR: ierr = Gcr(A, x, b, M, iter); break;
      case QMR_SYM: ierr = QmrSym(A, x, b, M, iter); break;
      case BICGSTABL: ierr = BiCgStabl(A, x, b, M, iter); break;
      case QCGS: ierr = QCgs(A, x, b, M, iter); break;
      case JACOBI: ierr = SolveJacobi(A, x, b, inv_diagonal, iter); break;
      case GAUSS_SEIDEL_RELAX: ierr = SolveSOR(A, x, b, omega_relax, 
                                               iter, stage_relax); break;
      default:
	{
	  cout << "Unavailable iterative solver " << type_solver << endl;
	  abort();
	}
      }

    nb_iterations = iter.GetNumberIteration();    
    return ierr;
  }


  //! return the number of iterative vectors needed by the solver
  template<class T> 
  int All_IterativeSolver<T>::GetNbIterativeVectors() const
  {
    int m = parameter_restart;
    switch (type_solver)
      {
      case GMRES: return m+4;
      case QMR: return 10;
      case COCG: return 4;
      case BICGCR: return 6;
      case BICGSTAB: return 8;
      case TFQMR: return 10;
      case CG: return 4;
      case CGNE: return 4;
      case LSQR: return 5;
      case CGS: return 9;
      case BICG: return 8;
      case MINRES: return 9;
      case SYMMLQ: return 8;
      case GCR: return 2*m+5;
      case QMR_SYM: return 7;
      case BICGSTABL: return 2*m+4;
      case QCGS: return 9;
      }
    
    return 0;
  }
  
  
  //! sets scaling (left and right)
  template<class T>
  void All_IterativeSolver<T>::SetScaling(VectReal_wp& diag_left, VectReal_wp& diag_right)
  {    
    diagonal_scale_left = true;
    diagonal_scale_right = true;
    diagonal_scale_left = diag_left;
    diagonal_scale_right = diag_right;
  }


  //! copies input parameters of an another solver
  template<class T>
  void All_IterativeSolver<T>::CopyParameter(const All_IterativeSolver<T>& solver)
  {
    this->print_level = solver.print_level;
    this->type_solver = solver.type_solver;
    this->diagonal_scaling_left = solver.diagonal_scaling_left;
    this->diagonal_scaling_right = solver.diagonal_scaling_right;
    this->nb_max_iterations = solver.nb_max_iterations;
    this->tolerance = solver.tolerance;
    this->parameter_restart = solver.parameter_restart;
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_SOLVE_SYSTEM_CXX
#endif
