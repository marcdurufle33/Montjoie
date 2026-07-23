#ifndef MONTJOIE_FILE_PRECONDITIONER_INLINE_CXX

namespace Seldon
{
#ifdef SELDON_WITH_HYPRE
  template<class T, class Prop, class Storage>
  inline void HyprePreconditioner<Complex_wp>::ConstructPreconditioner(DistributedMatrix<T, Prop, Storage>& A,
								       bool keep_matrix)
  {
    cout << "Hypre not compiled with complex" << endl;
    abort();
  }
    

  inline void HyprePreconditioner<Complex_wp>::
  Solve(const VirtualMatrix<Complex_wp>& A, const Vector<Complex_wp>& b, Vector<Complex_wp>& x)
  {
    cout << "Hypre not compiled with complex" << endl;
    abort();
  }


  inline void HyprePreconditioner<Complex_wp>
  ::TransSolve(const VirtualMatrix<Complex_wp>& A, const Vector<Complex_wp>& b, Vector<Complex_wp>& x)
  {
    cout << "Hypre not compiled with complex" << endl;
    abort();
  }
#endif
  
}

namespace Montjoie
{  
  template<class T>
  inline IdentityPreconditioner<T>::IdentityPreconditioner()
  {
  }
    
  template<class T>
  inline void IdentityPreconditioner<T>::SetInputData(const string& keyword, const VectString& param)
  {
  }
  
  template<class T>
  inline void IdentityPreconditioner<T>
  ::SetInputPreconditioning(const string& keyword, const VectString& param) 
  {
  }


#ifdef SELDON_WITH_PRECONDITIONING  
  template<class T>
  inline T JacobiPreconditioning<T>::GetRelaxationCoef() const 
  {
    return omega; 
  }


  //! default constructor
  template<class T> template<class TypeEquation>
  inline SubdomainPreconditioning_Base<T>
  ::SubdomainPreconditioning_Base(EllipticProblem<TypeEquation>& var)
    : var_subdomain(var)
  {
    InitDefaultValues();
  }
  
  
  template<class T>
  inline SubdomainPreconditioning_Base<T>::~SubdomainPreconditioning_Base() 
  {
  }

  template<class T>
  inline LowOrderPreconditioning_Base<T>::~LowOrderPreconditioning_Base()
  {
  }


  template<class T, class Dimension> template<class TypeEquation>
  inline SubdomainPreconditioning_Dim<T, Dimension>
  ::SubdomainPreconditioning_Dim(EllipticProblem<TypeEquation>& var)
    : SubdomainPreconditioning_Base<T>(var), var_subdomain(var), var_computation(var)
  {
  }


  //! constructor with the finite element object
  template<class T, class TypeEquation>
  inline SubdomainPreconditioning<T, TypeEquation>::
  SubdomainPreconditioning(EllipticProblem<TypeEquation>& var)
    : SubdomainPreconditioning_Dim<T, typename TypeEquation::Dimension>(var)
  {
  }



  //! default constructor
  template<class T, class Dimension> template<class TypeEquation>
  inline LowOrderPreconditioning_Dim<T, Dimension>
  ::LowOrderPreconditioning_Dim(EllipticProblem<TypeEquation>& var)
    : LowOrderPreconditioning_Base<T>(), var_problem(var)
  {
  }



  template<class T, class Dimension>
  inline void LowOrderPreconditioning_Dim<T, Dimension>
  ::ComputeLocalProlongation(const VectReal_wp&, DistributedProblem<Dimension>& var,
			     DistributedProblem<Dimension>& var_low)
  {
  }


  //! default constructor
  template<class T, class TypeEquation>
  inline LowOrderPreconditioning<T, TypeEquation>::LowOrderPreconditioning(EllipticProblem<TypeEquation>& var)
    : LowOrderPreconditioning_Dim<T, typename TypeEquation::Dimension>(var)
  {
  }
  
  

  //! constructor with a finite element object
  template<class T> template<class TypeEquation>
  inline MultigridPreconditioning_Base<T>
  ::MultigridPreconditioning_Base(EllipticProblem<TypeEquation>& var)
    : var_fine(var)
  {
    InitDefaultValues();
  }
  
  
  template<class T>
  inline void MultigridPreconditioning_Base<T>::SetMaximalOrder(int r)
  {
    maximal_order_multigrid = r;
  }


  //! constructor with a given problem
  template<class T, class Dimension> template<class TypeEquation>
  inline MultigridPreconditioning_Dim<T, Dimension>
  ::MultigridPreconditioning_Dim(EllipticProblem<TypeEquation>& var) :
    MultigridPreconditioning_Base<T>(var), var_fine(var) 
  {
    var_coarse = NULL;
  }

  
  template<class T, class Dimension>
  inline int MultigridPreconditioning_Dim<T, Dimension>
  ::GetCoarseOrder() const
  {
    return var_coarse->GetMeshNumbering(0).GetOrder(); 
  }  


  //! constructor with finite element object
  template<class T, class TypeEquation>
  inline MultigridPreconditioning<T, TypeEquation>::
  MultigridPreconditioning(EllipticProblem<TypeEquation>& var)
    : MultigridPreconditioning_Dim<T, typename TypeEquation::Dimension>(var)
  {
  }
  
  
  template<class T>
  inline MultigridPreconditioning_Base<T>& All_Preconditioner_Base<T>
  ::GetMultigridPreconditioning()
  {
    return *prec_multigrid; 
  }


  //! constructor with a finite element object
  template<class T> template<class TypeEquation>
  inline All_Preconditioner_Base<T>
  ::All_Preconditioner_Base(EllipticProblem<TypeEquation>& var)
  {
    prec_subdomain = new SubdomainPreconditioning<T, TypeEquation>(var);
    prec_multigrid = new MultigridPreconditioning<T, TypeEquation>(var);
    prec_low_order = new LowOrderPreconditioning<T, TypeEquation>(var);
    prec_local = new LocalPreconditioning<T, TypeEquation>(var);

    InitDefaultValues();
  }

  
  template<class T>
  inline void All_Preconditioner_Base<T>::SetPreconditioningType(int type)
  {
    type_precond = type;
  }


  template<class T> template<class Prop, class Storage, class Allocator>
  inline void All_Preconditioner_Base<T>::SetDirectPreconditioning(DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    type_precond = DIRECT;
    mat_lu.Factorize(A);
  }

  
  template<class T>
  inline void All_Preconditioner_Base<T>::Clear()
  {
  }


  //! constructor with finite element object
  template<class T, class TypeEquation>
  inline LocalPreconditioning<T, TypeEquation>::
  LocalPreconditioning(EllipticProblem<TypeEquation>& var)
    : LocalPreconditioning_Dim<T, typename TypeEquation::Dimension>()
  {
  }


#else

  //! constructor with a finite element object
  template<class T> template<class TypeEquation> inline
  All_Preconditioner_Base<T>::All_Preconditioner_Base(EllipticProblem<TypeEquation>& var)
  {
  }

  template<class T>
  inline All_Preconditioner_Base<T>::~All_Preconditioner_Base()
  {
  }

  template<class T>
  inline void All_Preconditioner_Base<T>::SetPreconditioningType(int type)
  {
    type_precond = type;
  }

  template<class T>
  inline size_t All_Preconditioner_Base<T>::GetMemorySize() const
  {
    return sizeof(*this); 
  }

  template<class T> template<class Dimension>
  inline void All_Preconditioner_Base<T>
  ::GetMaximumSize(T a, Dimension dim, int& nb_max_ilut,
		   int& nb_max_iterative, int& nb_max_multigrid,
		   int& order_min_multigrid) 
  {
  }
  
  template<class T>
  inline void All_Preconditioner_Base<T>::SetOptimalParameters(int order, int N, All_LinearSolver&) 
  {
  }
    
  template<class T>
  inline void All_Preconditioner_Base<T>
  ::ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
			    const FemMatrixFreeClass_Base<T>& A,
			    const GlobalGenericMatrix<T>& nature,
			    CondensationBlockSolver_Fem<T>& cond_solver) 
  {
  }
  
  template<class T>
  inline void All_Preconditioner_Base<T>::Clear() 
  {
  }

  template<class T, class TypeEquation>
  inline All_Preconditioner<T, TypeEquation>::All_Preconditioner(EllipticProblem<TypeEquation>& var)
    :  All_Preconditioner_Base<T>(var)
  {    
  }
#endif

}

#define MONTJOIE_FILE_PRECONDITIONER_INLINE_CXX
#endif
