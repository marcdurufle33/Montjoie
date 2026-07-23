#ifndef MONTJOIE_FILE_SOLVE_SYSTEM_INLINE_CXX

namespace Montjoie
{
  
  /****************
   * All_MatrixLU *
   ****************/
  
  
  //! Default constructor
  template <class T>
  inline All_MatrixLU<T>::All_MatrixLU() : SparseDistributedSolver<T>()
  {
  }
  

  /***********************
   * All_IterativeSolver *
   ***********************/

  
  //! Changes the level of displayed messages
  template <class T>
  inline void All_IterativeSolver<T>::SetPrintLevel(int print)
  {
    print_level = print;
  }
  

  //! returns the level of verbosity
  template <class T>
  inline int All_IterativeSolver<T>::GetPrintLevel() const
  {
    return print_level;
  }
  
    
  //! sets the iterative solver to use (Gmres, Cg, etc)
  template<class T> 
  inline void All_IterativeSolver<T>::SelectIterativeSolver(int type)
  {
    type_solver = type;
  }
  
  
  //! sets the maximum number of iterations
  template<class T> 
  inline void All_IterativeSolver<T>::SetMaxNumberIteration(int imax)
  {
    nb_max_iterations = imax;
  }
  
  
  //! returns the maximum number of iterations
  template<class T> 
  inline int All_IterativeSolver<T>::GetMaxNumberIteration() const
  {
    return nb_max_iterations;
  }

  
  //! returns the number of iterations used to converge
  template<class T> 
  inline int All_IterativeSolver<T>::GetNumberIteration() const
  {
    return nb_iterations;
  }

  
  //! sets the restart number (used for Gmres, Gcr, BicgStabl)
  template<class T> 
  inline void All_IterativeSolver<T>::SetRestart(int m)
  {
    parameter_restart = m;
  }


  //! returns the restart number (used for Gmres, Gcr, BicgStabl)
  template<class T> 
  inline int All_IterativeSolver<T>::GetRestart() const
  {
    return parameter_restart;
  }
  
  
  //! sets the stopping criterion
  template<class T> 
  inline void All_IterativeSolver<T>::SetStoppingCriterion(const Real_wp& tol)
  {
    tolerance = tol;
  }
  
  
  //! returns the stopping criterion
  template<class T> 
  inline Real_wp All_IterativeSolver<T>::GetStoppingCriterion() const
  {
    return tolerance;
  }
  
  
  //! sets preconditioning
  template<class T> 
  inline void All_IterativeSolver<T>::SetPreconditioning(Preconditioner_Base<T>& prec)
  {
    precond = &prec;
  }


  //! sets iterative algorithm used to solve the linear system
  template<class T> 
  inline void All_IterativeSolver<T>::SetSolverType(int type)
  {
    type_solver = type;
  }


  //! returns the solver used
  template<class T> 
  inline int All_IterativeSolver<T>::GetSolverType() const
  {
    return type_solver;
  }
  

  //! setting D^{-1} for jacobi solver
  template<class T> 
  inline void All_IterativeSolver<T>::SetInvDiagonalJacobi(const Vector<T>& invD)
  {
    inv_diagonal = invD;
  }
  
  
  //! returns memory used by the object in bytes
  template<class T> 
  inline size_t All_IterativeSolver<T>::GetMemorySize() const
  {
    size_t taille = sizeof(Real_wp)*(diagonal_scale_left.GetM()+diagonal_scale_right.GetM());
    taille += sizeof(*this);
    return taille;
  }
  
  
  //! Selects best iterative solver for unsymmetric matrix
  template<class T> template<class Dimension>
  inline void All_IterativeSolver<T>::
  SelectOptimalIterativeSolver(int N, T a, Dimension& dim, General prop)
  {
    type_solver = QMR;
  }
  
  
  //! Selects best iterative solver for symmetric matrix
  template<class T> template<class Dimension>
  inline void All_IterativeSolver<T>::
  SelectOptimalIterativeSolver(int N, T a, Dimension& dim, Symmetric prop)
  {
    type_solver = COCG;
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_SOLVE_SYSTEM_INLINE_CXX
#endif
