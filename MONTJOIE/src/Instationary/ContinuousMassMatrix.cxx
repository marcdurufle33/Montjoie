#ifndef MONTJOIE_FILE_CONTINUOUS_MASS_MATRIX_CXX

namespace Montjoie
{
    
  //! default constructor
  template<class Dimension>
  ContinuousMassMatrix_Base
  ::ContinuousMassMatrix_Base(VarProblem<Dimension>& var, int type_mat,
			      IVect& non_lumped_elt, IVect& non_lumped_dof, IVect& inverse_dof,
			      IVect& lumped_dof)
    : VirtualMatrix<Real_wp>(non_lumped_dof.GetM(), non_lumped_dof.GetM()),
    NonLumpedElement(non_lumped_elt), NonLumpedDof(non_lumped_dof), InverseDof(inverse_dof),
    LumpedDof(lumped_dof)
  {    
    type_matrix = type_mat;
    free_mat = var.GetNewIterativeMatrix(Real_wp(0));
    nodl_mesh = var.GetMeshNumbering(0).GetNbDof();
    
    InitDefaultValues();
  }


  //! initialization of parameters with default values
  void ContinuousMassMatrix_Base::InitDefaultValues()
  {
    stopping_criterion = 1e-6;
    nb_max_iteration = 100;
    iterative_solver = true;
    print_level = 0;
    
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefDamping(0);
    nat_mat.SetCoefMass(1);
    nat_mat.SetCoefStiffness(0);
    
    free_mat->SetCoefficientMatrix(nat_mat);
    free_mat->IgnoreDirichletDof();
  }


  //! destructor
  ContinuousMassMatrix_Base::~ContinuousMassMatrix_Base()
  {
    if (free_mat != NULL)
      delete free_mat;
  }


  //! returns the number of rows of the mass matrix
  int ContinuousMassMatrix_Base::GetNumberOfRows() const
  {
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL : 
        {
          if (diagonal.GetM() > 0)
            return diagonal.GetM();
          
          return invDiagonal.GetM();
        }
        break;
      case FemMassMatrix::BLOCK_DIAGONAL : 
        {
          if (block_diagonal.GetM() > 0)
            return block_diagonal.GetM();
          
          return invBlock_diagonal.GetM();
        }
        break;
      case FemMassMatrix::MATRIX_SPARSE :
        {
          return mat_sparse.GetM();
        }
        break;
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
        {
          return mat_sparse_unsym.GetM();
        }
        break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
        {
          if (block_diagonal_unsym.GetM() > 0)
            return block_diagonal_unsym.GetM();
          
          return invBlock_diagonal_unsym.GetM();
        }
        break;
      case FemMassMatrix::DIAG_SPARSE:
	{
	  int N = invDiagonal.GetM();
	  if (diagonal.GetM() > 0)
	    N = diagonal.GetM();

	  N += mat_sparse.GetM();	  
	  return N;
	}
	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM:
	{
	  int N = invDiagonal.GetM();
	  if (diagonal.GetM() > 0)
	    N = diagonal.GetM();

	  N += mat_sparse_unsym.GetM();	  
	  return N;
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE :
        {
	  int N = invBlock_diagonal.GetM();
          if (block_diagonal.GetM() > 0)
            N = block_diagonal.GetM();
	  
	  N += mat_sparse.GetM();
          return N;
        }
        break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM :
        {
	  int N = invBlock_diagonal_unsym.GetM();
          if (block_diagonal_unsym.GetM() > 0)
            N = block_diagonal_unsym.GetM();

	  N += mat_sparse_unsym.GetM();
          return N;
        }
        break;
      default :
        return free_mat->GetM();        
      }
  }
  

  //! returns the size used by the object in bytes
  size_t ContinuousMassMatrix_Base::GetMemorySize() const
  {
    size_t taille = diagonal.GetMemorySize() + invDiagonal.GetMemorySize()
      + block_diagonal.GetMemorySize() + invBlock_diagonal.GetMemorySize()
      + block_diagonal_unsym.GetMemorySize() + invBlock_diagonal_unsym.GetMemorySize()
      + precond_diag.GetMemorySize() + mat_sparse.GetMemorySize() + mat_lu.GetMemorySize()
      + mat_sparse_unsym.GetMemorySize() + previous_iterate.GetMemorySize()
      + Xn_tmp.GetMemorySize() + Bn_tmp.GetMemorySize();

    if (free_mat != NULL)
      taille += free_mat->GetMemorySize();

    return taille;
  }
  

  //! sets iterative solver to use
  void ContinuousMassMatrix_Base::SetIterativeSolver(bool iterative, int n, const Real_wp& eps,
						     int lvl)
  {
    iterative_solver = iterative;
    nb_max_iteration = n;
    stopping_criterion = eps;
    print_level = lvl;
  }
  

  void ContinuousMassMatrix_Base::SetPrintLevel(int lvl)
  {
    print_level = lvl;
  }

  
  //! multiplication by mass matrix Vh = M Uh
  void ContinuousMassMatrix_Base
  ::MltMass(const VectReal_wp& Uh,  VectReal_wp& Vh)
  {
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (diagonal.GetM() <= 0)
            {
	      cout << "Diagonal not present" << endl;
              abort();
            }
          
	  for (int i = 0; i < diagonal.GetM(); i++)
	    Vh(i) = Uh(i)*diagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }

	  Seldon::MltAdd(Real_wp(1), block_diagonal, Uh, Real_wp(0), Vh, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal_unsym.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }

	  Seldon::MltAdd(Real_wp(1), block_diagonal_unsym, Uh, Real_wp(0), Vh, false);
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE :
	{
	  // exiting if the matrix has not been computed
	  if (mat_sparse.GetM() <= 0)
	    {
	      cout << "Sparse matrix not present" << endl;
              abort();
            }

	  Seldon::MltAdd(Real_wp(1), mat_sparse, Uh, Real_wp(0), Vh, false);
        }
	break;
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
	{
	  // exiting if the matrix has not been computed
	  if (mat_sparse_unsym.GetM() <= 0)
	    {
	      cout << "Sparse matrix not present" << endl;
              abort();
            }

	  Seldon::MltAdd(Real_wp(1), mat_sparse_unsym, Uh, Real_wp(0), Vh, false);
        }
	break;
      case FemMassMatrix::DIAG_SPARSE :
	{
	  int N = mat_sparse.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltSolution(Xn_tmp, Vh);
	    }
	  else
	    for (int i = 0; i < diagonal.GetM(); i++)
	      Vh(i) = Uh(i)*diagonal(i);
	}
	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM :
	{
	  int N = mat_sparse_unsym.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse_unsym, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltSolution(Xn_tmp, Vh);
	    }
	  else
	    for (int i = 0; i < diagonal.GetM(); i++)
	      Vh(i) = Uh(i)*diagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE :
	{
	  int N = mat_sparse.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltBlockSolution(block_diagonal, Xn_tmp, Uh, Vh);
	    }
	  else
	    Seldon::MltAdd(Real_wp(1), block_diagonal, Uh, Real_wp(0), Vh, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM :
	{
	  int N = mat_sparse_unsym.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse_unsym, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltBlockSolution(block_diagonal_unsym, Xn_tmp, Uh, Vh);
	    }
	  else
	    Seldon::MltAdd(Real_wp(1), block_diagonal_unsym, Uh, Real_wp(0), Vh, false);
	}
	break;
      default :
	{
          cout << "storage not handled " << endl;
          abort();
        }
	break;
      }
  }
  
  
  //! computation of B = B + alpha*M*X
  void ContinuousMassMatrix_Base
  ::MltMass(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& B)
  {
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (diagonal.GetM() <= 0)
	    {
	      cout << "Diagonal not present" << endl;
              abort();
            }
	  
	  if (alpha == Real_wp(1))
	    for (int i = 0; i < diagonal.GetM(); i++)
	      B(i) += diagonal(i)*X(i);
	  else if (alpha == Real_wp(-1))
	    for (int i = 0; i < diagonal.GetM(); i++)
              B(i) -= diagonal(i)*X(i);
	  else
	    for (int i = 0; i < diagonal.GetM(); i++)
	      B(i) += alpha*diagonal(i)*X(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }
	  
	  Seldon::MltAdd(alpha, block_diagonal, X, Real_wp(1), B, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal_unsym.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }
	  
	  Seldon::MltAdd(alpha, block_diagonal_unsym, X, Real_wp(1), B, false);
        }
	break;
      case FemMassMatrix::MATRIX_SPARSE :
	{
	  // exiting if the matrix has not been computed
	  if (mat_sparse.GetM() <= 0)
	    {
              cout << "Sparse matrix not present" << endl;
	      abort();
            }
          
	  Seldon::MltAdd(alpha, mat_sparse, X, Real_wp(1), B, false);
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
	{
	  // exiting if the matrix has not been computed
	  if (mat_sparse_unsym.GetM() <= 0)
	    {
              cout << "Sparse matrix not present" << endl;
	      abort();
            }
          
	  Seldon::MltAdd(alpha, mat_sparse_unsym, X, Real_wp(1), B, false);
	}
	break;
      case FemMassMatrix::DIAG_SPARSE:
	{
	  int N = mat_sparse.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(X, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltAddSolution(alpha, Xn_tmp, X, B);
	    }
	  else
	    {
	      if (alpha == Real_wp(1))
		for (int i = 0; i < diagonal.GetM(); i++)
		  B(i) += diagonal(i)*X(i);
	      else if (alpha == Real_wp(-1))
		for (int i = 0; i < diagonal.GetM(); i++)
		  B(i) -= diagonal(i)*X(i);
	      else
		for (int i = 0; i < diagonal.GetM(); i++)
		  B(i) += alpha*diagonal(i)*X(i);
	    }
	}
	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM:
	{
	  int N = mat_sparse_unsym.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(X, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse_unsym, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltAddSolution(alpha, Xn_tmp, X, B);
	    }
	  else
	    {
	      if (alpha == Real_wp(1))
		for (int i = 0; i < diagonal.GetM(); i++)
		  B(i) += diagonal(i)*X(i);
	      else if (alpha == Real_wp(-1))
		for (int i = 0; i < diagonal.GetM(); i++)
		  B(i) -= diagonal(i)*X(i);
	      else
		for (int i = 0; i < diagonal.GetM(); i++)
		  B(i) += alpha*diagonal(i)*X(i);
	    }
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE:
	{
	  int N = mat_sparse.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(X, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltAddBlockSolution(alpha, block_diagonal, Xn_tmp, X, B);
	    }
	  else
	    Seldon::MltAdd(alpha, block_diagonal, X, Real_wp(1), B, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM:
	{
	  int N = mat_sparse_unsym.GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(X, Bn_tmp);
	      Seldon::MltAdd(Real_wp(1), mat_sparse_unsym, Bn_tmp, Real_wp(0), Xn_tmp, false);
	      this->MltAddBlockSolution(alpha, block_diagonal_unsym, Xn_tmp, X, B);
	    }
	  else
	    Seldon::MltAdd(alpha, block_diagonal_unsym, X, Real_wp(1), B, false);
	}
	break;
      default :
        {
          cout << "storage not handled " << endl;
          abort();
        }
	break;
      }
  }
  

  //! solves symmetric system A X = B
  void ContinuousMassMatrix_Base::SolveSparseMass(const VectReal_wp& B, VectReal_wp& X)
  {
    if (iterative_solver)
      {
	Iteration<Real_wp> iter(nb_max_iteration, stopping_criterion);
	iter.HideMessages();
	int rank_proc = 0;
#ifdef SELDON_WITH_MPI
        int  nb_proc = 1;
	MPI_Comm_rank(mat_sparse.GetCommunicator(), &rank_proc);
	MPI_Comm_size(mat_sparse.GetCommunicator(), &nb_proc);
#endif
	if (rank_proc == 0)
	  {
	    if (print_level == 1)
	      iter.ShowMessages();
	    else if (print_level >= 2)
	      iter.ShowFullHistory();
	  }
	
	X.Zero();
	if (this->previous_iterate.GetM() > 0)
	  {
	    iter.SetInitGuess(false);
	    Copy(this->previous_iterate, X);
	  }
	else
	  iter.SetInitGuess(true);

#ifdef SELDON_WITH_MPI
	if (nb_proc > 1)
	  {
	    DistributedVector<Real_wp> Bpar(mat_sparse.GetOverlapRowNumber(), mat_sparse.GetCommunicator());
	    DistributedVector<Real_wp> Xpar(mat_sparse.GetOverlapRowNumber(), mat_sparse.GetCommunicator());
	    Bpar.SetData(B.GetM(), B.GetData());
	    Xpar.SetData(X.GetM(), X.GetData());
	    
	    Cg(mat_sparse, Xpar, Bpar, *this, iter);
	    
	    Bpar.Nullify();
	    Xpar.Nullify();
	  }
	else
#endif
	  Cg(mat_sparse, X, B, *this, iter);
      }
    else
      {
	X = B;
	mat_lu.Solve(X);
      }
  }


  //! solves unsymmetric system A X = B
  void ContinuousMassMatrix_Base::SolveSparseMassUnsym(const VectReal_wp& B, VectReal_wp& X)
  {
    if (iterative_solver)
      {
	Iteration<Real_wp> iter(nb_max_iteration, stopping_criterion);
	int rank_proc = 0;
	iter.HideMessages();
	
#ifdef SELDON_WITH_MPI
        int nb_proc = 1;
	MPI_Comm_rank(mat_sparse_unsym.GetCommunicator(), &rank_proc);
	MPI_Comm_size(mat_sparse_unsym.GetCommunicator(), &nb_proc);
#endif
	if (rank_proc == 0)
	  {
	    if (print_level == 1)
	      iter.ShowMessages();
	    else if (print_level >= 2)
	      iter.ShowFullHistory();
	  }
	
	X.Zero();
	if (this->previous_iterate.GetM() > 0)
	  {
	    iter.SetInitGuess(false);
	    Copy(this->previous_iterate, X);
	  }
	else
	  iter.SetInitGuess(true);

#ifdef SELDON_WITH_MPI
	if (nb_proc > 1)
	  {
	    DistributedVector<Real_wp> Bpar(mat_sparse_unsym.GetOverlapRowNumber(), mat_sparse_unsym.GetCommunicator());
	    DistributedVector<Real_wp> Xpar(mat_sparse_unsym.GetOverlapRowNumber(), mat_sparse_unsym.GetCommunicator());
	    Bpar.SetData(B.GetM(), B.GetData());
	    Xpar.SetData(X.GetM(), X.GetData());
	    
	    TfQmr(mat_sparse_unsym, Xpar, Bpar, *this, iter);
	    
	    Bpar.Nullify();
	    Xpar.Nullify();
	  }
	else
#endif
	  TfQmr(mat_sparse_unsym, X, B, *this, iter);
      }
    else
      {
	X = B;
	mat_lu.Solve(X);
      }
  }

  
  //! computes Vh = M^-1 Vh where M is the mass matrix
  void ContinuousMassMatrix_Base
  ::SolveMass(VectReal_wp& Vh)
  {  
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invDiagonal.GetM() <= 0)
            {
	      cout << "Diagonal not present" << endl;
              abort();
            }
          
	  for (int i = 0; i < invDiagonal.GetM(); i++)
	    Vh(i) *= invDiagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invBlock_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }
	  
          Seldon::Mlt(invBlock_diagonal, Vh, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
	{
	  // exiting if the matrix has not been computed
	  if (invBlock_diagonal_unsym.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }
	  
          Seldon::Mlt(invBlock_diagonal_unsym, Vh, false);
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE :
	{
	  int N = this->GetM();
	  if (Xn_tmp.GetM() != N)
	    Xn_tmp.Reallocate(N);

	  SolveSparseMass(Vh, Xn_tmp);
	  Vh = Xn_tmp;
	}	
	break;
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
	{
	  int N = this->GetM();
	  if (Xn_tmp.GetM() != N)
	    Xn_tmp.Reallocate(N);

	  SolveSparseMassUnsym(Vh, Xn_tmp);
	  Vh = Xn_tmp;
	}	
	break;	
      case FemMassMatrix::DIAG_SPARSE :
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Vh, Bn_tmp);
	      SolveSparseMass(Bn_tmp, Xn_tmp);
	      this->GetSolution(Xn_tmp, Vh);
	    }
	  else
	    for (int i = 0; i < invDiagonal.GetM(); i++)
	      Vh(i) *= invDiagonal(i);	  
	}
	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM:
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Vh, Bn_tmp);
	      SolveSparseMassUnsym(Bn_tmp, Xn_tmp);
	      this->GetSolution(Xn_tmp, Vh);
	    }
	  else
	    for (int i = 0; i < invDiagonal.GetM(); i++)
	      Vh(i) *= invDiagonal(i);	  
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE:
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Vh, Bn_tmp);
	      SolveSparseMass(Bn_tmp, Xn_tmp);
	      this->GetBlockSolution(invBlock_diagonal, Xn_tmp, Vh);
	    }
	  else
	    Seldon::Mlt(invBlock_diagonal, Vh, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM:
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Vh, Bn_tmp);
	      SolveSparseMassUnsym(Bn_tmp, Xn_tmp);
	      this->GetBlockSolution(invBlock_diagonal_unsym, Xn_tmp, Vh);
	    }
	  else
	    Seldon::Mlt(invBlock_diagonal_unsym, Vh, false);
	}
	break;
      default :
	{
          cout << "storage not handled " << endl;
          abort();
        }
	break;
      }
  }

  
  //! computes Vh = Vh + alpha*M^-1 Uh
  void ContinuousMassMatrix_Base
  ::SolveMass(const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {  
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invDiagonal.GetM() <= 0)
            {
	      cout << "Diagonal not present" << endl;
              abort();
            }
          
	  for (int i = 0; i < invDiagonal.GetM(); i++)
	    Vh(i) += alpha*Uh(i)*invDiagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invBlock_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }
	  
          Seldon::MltAdd(alpha, invBlock_diagonal, Uh, Real_wp(1), Vh, false);
          //VectReal_wp Uh2 = Uh;
	  //Seldon::SolveCholesky(SeldonNoTrans, invBlock_diagonal, Uh2);
          //Seldon::SolveCholesky(SeldonTrans, invBlock_diagonal, Uh2);
          //Add(alpha, Uh2, Vh);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
	{
	  // exiting if the matrix has not been computed
	  if (invBlock_diagonal_unsym.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
              abort();
            }
	  
          Seldon::MltAdd(alpha, invBlock_diagonal_unsym, Uh, Real_wp(1), Vh, false);
          //VectReal_wp Uh2 = Uh;
	  //Seldon::SolveCholesky(SeldonNoTrans, invBlock_diagonal, Uh2);
          //Seldon::SolveCholesky(SeldonTrans, invBlock_diagonal, Uh2);
          //Add(alpha, Uh2, Vh);
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE :
	{
	  int N = this->GetM();
	  if (Xn_tmp.GetM() != N)
	    Xn_tmp.Reallocate(N);

	  SolveSparseMass(Uh, Xn_tmp);
	  Vh += alpha*Xn_tmp;
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
	{
	  int N = this->GetM();
	  if (Xn_tmp.GetM() != N)
	    Xn_tmp.Reallocate(N);

	  SolveSparseMassUnsym(Uh, Xn_tmp);
	  Vh += alpha*Xn_tmp;
	}
	break;
	case FemMassMatrix::DIAG_SPARSE :
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      SolveSparseMass(Bn_tmp, Xn_tmp);
	      this->AddSolution(alpha, Xn_tmp, Uh, Vh);
	    }
	  else
	    for (int i = 0; i < invDiagonal.GetM(); i++)
	      Vh(i) += alpha*Uh(i)*invDiagonal(i);
	}
	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM:
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      SolveSparseMassUnsym(Bn_tmp, Xn_tmp);
	      this->AddSolution(alpha, Xn_tmp, Uh, Vh);
	    }
	  else
	    for (int i = 0; i < invDiagonal.GetM(); i++)
	      Vh(i) += alpha*Uh(i)*invDiagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE:
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      SolveSparseMass(Bn_tmp, Xn_tmp);
	      this->AddBlockSolution(alpha, invBlock_diagonal, Xn_tmp, Uh, Vh);
	    }
	  else
	    Seldon::MltAdd(alpha, invBlock_diagonal, Uh, Real_wp(1), Vh, false);
	}
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM:
	{
	  int N = this->GetM();
	  if (N > 0)
	    {
	      if (Bn_tmp.GetM() != N)
		Bn_tmp.Reallocate(N);
	      
	      if (Xn_tmp.GetM() != N)
		Xn_tmp.Reallocate(N);
	      
	      this->GetSource(Uh, Bn_tmp);
	      SolveSparseMassUnsym(Bn_tmp, Xn_tmp);
	      this->AddBlockSolution(alpha, invBlock_diagonal_unsym, Xn_tmp, Uh, Vh);
	    }
	  else
	    Seldon::MltAdd(alpha, invBlock_diagonal_unsym, Uh, Real_wp(1), Vh, false);
	}
	break;
      default :
	{
          cout << "storage not handled " << endl;
          abort();
        }
	break;
      }
  }


  //! computes diagonal preconditioning from sparse matrix A
  template<class Prop, class Storage>
  void ContinuousMassMatrix_Base::ComputePreconditioning(DistributedMatrix<Real_wp, Prop, Storage>& A)
  {
    precond_diag.Reallocate(A.GetM());
    if (A.GetM() > 0)
      {
	for (int i = 0; i < A.GetM(); i++)
	  precond_diag(i) = A(i, i);

#ifdef SELDON_WITH_MPI
	A.AssembleVec(precond_diag);
#endif
	
	for (int i = 0; i < A.GetM(); i++)
	  precond_diag(i) = 1.0/precond_diag(i);
      }
  }

  
  //! Inversion of mass matrix such that Solve can be called afterwards
  void ContinuousMassMatrix_Base::Invert(bool keep_mass)
  {
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
        invDiagonal.Reallocate(diagonal.GetM());
        for (int i = 0; i < diagonal.GetM(); i++)
          invDiagonal(i) = 1.0/diagonal(i);
        
        if (!keep_mass)
          diagonal.Clear();
        
        break;
      case FemMassMatrix::BLOCK_DIAGONAL :
        invBlock_diagonal = block_diagonal;
        GetInverse(invBlock_diagonal);
        
        if (!keep_mass)
          block_diagonal.Clear();
        
        break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
        invBlock_diagonal_unsym = block_diagonal_unsym;
        GetInverse(invBlock_diagonal_unsym);
        
        if (!keep_mass)
          block_diagonal_unsym.Clear();
        
        break;
      case FemMassMatrix::DIAG_SPARSE :
      case FemMassMatrix::BLOCK_DIAG_SPARSE :
      case FemMassMatrix::MATRIX_SPARSE :
	this->m_ = mat_sparse.GetM();
	this->n_ = mat_sparse.GetM();
	if (iterative_solver)
	  ComputePreconditioning(mat_sparse);
	else
	  {
	    if (this->m_ > 0)
	      {
                if (print_level == 1)
                  mat_lu.ShowMessages();
                else
                  mat_lu.HideMessages();
                
                mat_lu.Factorize(mat_sparse, keep_mass);

                int ierr;
                int type = mat_lu.GetInfoFactorization(ierr);
                PrintFactorizationFailed(type, ierr);
                mat_lu.HideMessages();
              }
	    else
	      mat_lu.Clear();
	  }

	if ((type_matrix == FemMassMatrix::BLOCK_DIAG_SPARSE)
	    && (block_diagonal.GetM() > 0))
	  {
	    invBlock_diagonal = block_diagonal;
	    GetInverse(invBlock_diagonal);
	    
	    if (!keep_mass)
	      block_diagonal.Clear();
	  }

	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM :
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM :
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
	this->m_ = mat_sparse_unsym.GetM();
	this->n_ = mat_sparse_unsym.GetM();
	if (iterative_solver)
	  ComputePreconditioning(mat_sparse_unsym);
	else
	  {
	    if (this->m_ > 0)
              {
                if (print_level == 1)
                  mat_lu.ShowMessages();
                else
                  mat_lu.HideMessages();
                
                mat_lu.Factorize(mat_sparse_unsym, keep_mass);

                int ierr;
                int type = mat_lu.GetInfoFactorization(ierr);
                PrintFactorizationFailed(type, ierr);
                mat_lu.HideMessages();
              }
	    else
	      mat_lu.Clear();
	  }

	if ((type_matrix == FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM)
	    && (block_diagonal_unsym.GetM() > 0))
	  {
	    invBlock_diagonal_unsym = block_diagonal_unsym;
	    GetInverse(invBlock_diagonal_unsym);
	    
	    if (!keep_mass)
	      block_diagonal_unsym.Clear();
	  }
	
	break;
      default :
	{
	  cout << "Storage not implemented" << endl;
	  abort();
	}
      }
  }
  
  
  //! extracts right hand side for non-lumped dofs
  void ContinuousMassMatrix_Base
  ::GetSource(const VectReal_wp& Uh, VectReal_wp& B) const
  {
    B.Reallocate(NonLumpedDof.GetM());
    for (int i = 0; i < NonLumpedDof.GetM(); i++)
      B(i) = Uh(NonLumpedDof(i));    
  }
  
  
  //! adds the solution of M*X = U, when U if partly filled with the solution for non-lumped dofs
  void ContinuousMassMatrix_Base
  ::AddSolution(const Real_wp& alpha, const VectReal_wp& Uh,
                const VectReal_wp& X, VectReal_wp& Y) const
  {
    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  Y(i) += alpha*invDiagonal(num++)*X(i);
	else
	  Y(i) += alpha*Uh(InverseDof(i));
      }
  }


  //! adds the solution of M*X = U, when U if partly filled with the solution for non-lumped dofs
  template<class Prop, class Storage>
  void ContinuousMassMatrix_Base
  ::AddBlockSolution(const Real_wp& alpha,
		     const DistributedMatrix<Real_wp, Prop, Storage>& M,
		     const VectReal_wp& Uh, const VectReal_wp& X, VectReal_wp& Y) const
  {
    VectReal_wp B;
    ExtractBlockSolution(M, X, B);

    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  Y(i) += alpha*B(num++);
	else
	  Y(i) += alpha*Uh(InverseDof(i));
      }
  }


  //! fills the solution of M*X = U, when U if partly filled with the solution for non-lumped dofs
  void ContinuousMassMatrix_Base
  ::GetSolution(const VectReal_wp& Uh, VectReal_wp& X) const
  {
    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  X(i) *= invDiagonal(num++);
	else
	  X(i) = Uh(InverseDof(i));
      }
  }


  //! fills the solution of M*X = U, when U if partly filled with the solution for non-lumped dofs
  template<class Prop, class Storage>
  void ContinuousMassMatrix_Base
  ::GetBlockSolution(const DistributedMatrix<Real_wp, Prop, Storage>& M,
		     const VectReal_wp& Uh, VectReal_wp& X) const
  {
    VectReal_wp B;
    ExtractBlockSolution(M, X, B);
  
    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  X(i) = B(num++);
	else
	  X(i) = Uh(InverseDof(i));
      }
  }


  //! fills the vector X = M*U
  void ContinuousMassMatrix_Base
  ::MltSolution(const VectReal_wp& Uh, VectReal_wp& X) const
  {
    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  X(i) *= diagonal(num++);
	else
	  X(i) = Uh(InverseDof(i));
      }
  }


  //! extracts block-diagonal from U and multiplies it with M, the result is stored in B
  template<class Prop, class Storage>
  void ContinuousMassMatrix_Base
  ::ExtractBlockSolution(const DistributedMatrix<Real_wp, Prop, Storage>& M,
			 const VectReal_wp& U, VectReal_wp& B) const
  {
    if (LumpedDof.GetM() > 0)
      {
	B.Reallocate(LumpedDof.GetM());
	for (int i = 0; i < LumpedDof.GetM(); i++)
	  B(i) = U(LumpedDof(i));    
	
	Mlt(M, B, false);
      }
  }
  

  //! fills the vector X = M*U
  template<class Prop, class Storage>
  void ContinuousMassMatrix_Base
  ::MltBlockSolution(const DistributedMatrix<Real_wp, Prop, Storage>& M,
		     const VectReal_wp& Uh, const VectReal_wp& U, VectReal_wp& X) const
  {
    VectReal_wp B;
    ExtractBlockSolution(M, U, B);
    
    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  X(i) = B(num++);
	else
	  X(i) = Uh(InverseDof(i));
      }
  }


  //! adds alpha*M*U to Y
  void ContinuousMassMatrix_Base
  ::MltAddSolution(const Real_wp& alpha, const VectReal_wp& Uh,
		   const VectReal_wp& X, VectReal_wp& Y) const
  {
    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  Y(i) += alpha*diagonal(num++)*X(i);
	else
	  Y(i) += alpha*Uh(InverseDof(i));
      }
  }


  //! adds alpha*M*U to Y
  template<class Prop, class Storage>
  void ContinuousMassMatrix_Base
  ::MltAddBlockSolution(const Real_wp& alpha, const DistributedMatrix<Real_wp, Prop, Storage>& M,
			const VectReal_wp& Uh, const VectReal_wp& X, VectReal_wp& Y) const
  {
    VectReal_wp B;
    ExtractBlockSolution(M, X, B);

    int num = 0;
    for (int i = 0; i < InverseDof.GetM(); i++)
      {
	if (InverseDof(i) == -1)
	  Y(i) += alpha*B(num++);
	else
	  Y(i) += alpha*Uh(InverseDof(i));
      }
  }

  
  //! application of diagonal preconditioning
  void ContinuousMassMatrix_Base
  ::Solve(const VirtualMatrix<Real_wp>& A, const VectReal_wp& r, VectReal_wp& z)
  {
    for (int i = 0; i < z.GetM(); i++)
      z(i) = r(i)*precond_diag(i);
  }

  
  //! application of diagonal preconditioning
  void ContinuousMassMatrix_Base
  ::TransSolve(const VirtualMatrix<Real_wp>& A, const VectReal_wp& r, VectReal_wp& z)
  {
    for (int i = 0; i < z.GetM(); i++)
      z(i) = r(i)*precond_diag(i);
  }

  
  /************************
   * ContinuousMassMatrix *
   ************************/


  template<class Dimension>
  ContinuousMassMatrix<Dimension>
  ::ContinuousMassMatrix(VarProblem<Dimension>& var, int type_mat,
			 IVect& non_lumped_elt, IVect& non_lumped_dof, IVect& inverse_dof,
			 IVect& lumped_dof)
    : ContinuousMassMatrix_Base(var, type_mat, non_lumped_elt, non_lumped_dof, inverse_dof, lumped_dof),
      mesh_num(var.GetMeshNumbering(0)), mesh(var.mesh), var_problem(var)
  {
  }
  

  //! matrix-vector product for a single element
  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltAddElement(const Real_wp& alpha, int i, int ref,
		  VectReal_wp& Uh, VectReal_wp& Vh,
		  IVect& num_ddl, const VectReal_wp& X, VectReal_wp& Y,
		  const ElementReference_Dim<Dimension>& Fb) const
  {
    /*
    int nb_dof = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    num_ddl.Reallocate(nb_dof); Uh.Reallocate(nb_dof);
    Vh.Reallocate(nb_quad);
    for (int j = 0; j < nb_dof; j++)
      {
	num_ddl(j) = InverseDof(mesh_num.Element(i).GetNumberDof(j));
	Uh(j) = X(num_ddl(j));
      }
    
    Fb.ApplyChTranspose(Uh, Vh);
    
    enum {nb_components_u = ElementReference<Dimension, type>::nb_components_u};
    TinyVector<Real_wp, nb_components_u> rho;
    bool variable = var_problem.UseNumericalIntegration(i);
    //bool affine = mesh.IsElementAffine(i);
    if (variable)
      {	
        if (affine)
          for (int m = 0; m < nb_components_u; m++)
            for (int j = 0; j < nb_quad; j++)
              {
                rho(m) = TypeEquation::GetMassCoefficient(var_time, i, j, m, ref);
                Vh(j*TypeElement::nb_components_u+m) 
                  *= var.Glob_jacobian(i)(0)*Fb.WeightsND(j)*rho(m);
              }
        else
          for (int m = 0; m < TypeElement::nb_components_u; m++)
            for (int j = 0; j < nb_quad; j++)
              {
                rho(m) = TypeEquation::GetMassCoefficient(var_time, i, j, m, ref);
                Vh(j*TypeElement::nb_components_u+m) *= var.Glob_jacobian(i)(j)*rho(m);
              }
      }
    else
      {
      for (int m = 0; m < TypeElement::nb_components_u; m++)
	  rho(m) = TypeEquation::GetMassCoefficient(var_time, i, 0, m, ref);
	
	for (int m = 0; m < TypeElement::nb_components_u; m++)
	  {
	    Real_wp mass = var.Glob_jacobian(i)(0)*rho(m);
	    for (int j = 0; j < nb_quad; j++)
	      Vh(j*TypeElement::nb_components_u+m) *= Fb.WeightsND(j)*mass;
	  }
      }
    
    Fb.ApplyCh(Vh, Uh);
    for (int j = 0; j < nb_dof; j++)
      Y(num_ddl(j)) += alpha*Uh(j);
    */

    cout << "Not implemented" << endl;
    abort();
  }
  
  
  //! matrix vector product with mass matrix Y = beta*Y + alpha*M*X
  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltAddVector(const Real_wp& alpha, const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Y) const
  {
    // diagonal part
    for (int i = 0; i < X.GetM(); i++)
      Y(i) = beta*Y(i) + alpha*diagonal(NonLumpedDof(i))*X(i);
    
    // and other elements
    VectReal_wp Uh, Vh; IVect num_ddl;
    for (int i1 = 0; i1 < this->NonLumpedElement.GetM(); i1++)
      {
	int i = this->NonLumpedElement(i1);
	int ref = mesh.Element(i).GetReference();
	MltAddElement(alpha, i, ref, Uh, Vh, num_ddl, X, Y,
		      var_problem.GetReferenceElement(i));
      }
  }


  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltAddVector(const Complex_wp& alpha, const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltAddVector(const Real_wp& alpha, const SeldonTranspose&, const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltAddVector(const Complex_wp& alpha, const SeldonTranspose&, const VectComplex_wp& X, const Complex_wp& beta, VectComplex_wp& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  //! matrix vector product with mass matrix Y = beta*Y + alpha*M*X
  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltVector(const VectReal_wp& X, VectReal_wp& Y) const
  {
    // diagonal part
    for (int i = 0; i < X.GetM(); i++)
      Y(i) = diagonal(this->NonLumpedDof(i))*X(i);
    
    // and other elements
    VectReal_wp Uh, Vh; IVect num_ddl;
    for (int i1 = 0; i1 < this->NonLumpedElement.GetM(); i1++)
      {
	int i = this->NonLumpedElement(i1);
	int ref = mesh.Element(i).GetReference();
	MltAddElement(Real_wp(1), i, ref, Uh, Vh, num_ddl, X, Y,
		      var_problem.GetReferenceElement(i));
      }
  }


  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltVector(const VectComplex_wp& X, VectComplex_wp& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltVector(const SeldonTranspose&, const VectReal_wp& X, VectReal_wp& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void ContinuousMassMatrix<Dimension>
  ::MltVector(const SeldonTranspose&, const VectComplex_wp& X, VectComplex_wp& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
}

#define MONTJOIE_FILE_CONTINUOUS_MASS_MATRIX_CXX
#endif
