#ifndef MONTJOIE_FILE_SOLVE_SYSTEM_HXX

namespace Montjoie
{

  //! general class for direct solver
  template<class T>
  class All_MatrixLU : public SparseDistributedSolver<T>
  {
  public :
    All_MatrixLU();
    
    static string GetAvailableSolver();
    
    template<class Dimension>
    void SelectOptimalDirectSolver(int N, T a, Dimension& dim);
    
    void SelectOptimalDirectSolver();
    
    bool SelectDirectSolver(const VectString& parameters);
    void SelectDirectSolver(int);
    void SelectOrdering(const string&);
    
    void GetMaximumSize(Real_wp a, Dimension2 dim, int& nb_max_direct,
			int& nb_max_no_condensation, int& order_min_condensation);
    
    void GetMaximumSize(Real_wp a, Dimension3 dim, int& nb_max_direct,
			int& nb_max_no_condensation, int& order_min_condensation);
    
    void GetMaximumSize(complex<Real_wp> a, Dimension2 dim, int& nb_max_direct,
			int& nb_max_no_condensation, int& order_min_condensation);
    
    void GetMaximumSize(complex<Real_wp> a, Dimension3 dim, int& nb_max_direct,
			int& nb_max_no_condensation, int& order_min_condensation);

    void CopyParameter(const All_MatrixLU<T>&);
    
  };
  
  
  //! ghost class storing a matrix in order to handle solving the transpose system
  template<class T>
  class VirtualMatrixTranspose : public VirtualMatrix<T>
  {
    const VirtualMatrix<T>* A_ptr;
    SeldonTranspose trans;
    
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

  public:
    VirtualMatrixTranspose(const VirtualMatrix<T>& A, const SeldonTranspose& t);

    virtual void MltAddVector(const Treal& alpha, const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;

    virtual void MltAddVector(const Treal& alpha, const SeldonTranspose&,
			      const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const SeldonTranspose&,
			      const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const Vector<Treal>& x, Vector<Treal>& y) const;
    virtual void MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Treal>& x, Vector<Treal>& y) const;

    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Tcplx>& x, Vector<Tcplx>& y) const;

  };

  
  //! general class for iterative solver
  template<class T>
  class All_IterativeSolver
  {
  protected :
    int print_level; //!< level of display messages
    int type_solver; //!< iterative solver used (Gmres, Cg, etc)
    bool diagonal_scaling_left; //!< left scaling ?
    bool diagonal_scaling_right; //!< right scaling ?
    VectReal_wp diagonal_scale_left; //!< left scaling
    VectReal_wp diagonal_scale_right; //!< right scaling
    int nb_max_iterations; //!< maximum number of iterations
    int nb_iterations;
    Real_wp tolerance; //!< stopping criterion
    int parameter_restart; //!< restart parameter
    Real_wp omega_relax; int stage_relax;
    Vector<T> inv_diagonal;
    
    //! used preconditioning
    Preconditioner_Base<T>* precond;
    //! identity preconditioning
    Preconditioner_Base<T> identity_precond;
    
  public :
    // types of iterative solver
    enum { QMR, TFQMR, GMRES, CG, CGNE, LSQR, CGS,
	   BICG, BICGSTAB, BICGCR, MINRES, GCR, SYMMLQ,
	   COCG, QMR_SYM, BICGSTABL, QCGS, PRECOND,
           JACOBI, GAUSS_SEIDEL_RELAX};
    
    All_IterativeSolver();
    All_IterativeSolver(Preconditioner_Base<T>&);

    void SetInputData(const string& description_field, const VectString& parameters);    

    template<class Vector1>
    int Solve(const VirtualMatrix<T>& mat, Vector1 & x_sol, const Vector1 & source_rhs);

    template<class Vector1>
    int Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& mat,
              Vector1 & x_sol, const Vector1 & source_rhs);
    
    template <class Vector1>
    int Solve(const VirtualMatrix<T>& A,  Vector1& x, const Vector1& b,
	      Preconditioner_Base<T>& M, Iteration<typename ClassComplexType<T>::Treal> & iter);

    int GetNbIterativeVectors() const;
    void Clear();

    void CopyParameter(const All_IterativeSolver<T>&);
    
    // Inline methods
    void SetPrintLevel(int print);
    int GetPrintLevel() const;
    void SelectIterativeSolver(int);
    void SetMaxNumberIteration(int);
    int GetMaxNumberIteration() const;
    int GetNumberIteration() const;
    void SetRestart(int);
    int GetRestart() const;
    void SetStoppingCriterion(const Real_wp& );
    Real_wp GetStoppingCriterion() const;
    void SetPreconditioning(Preconditioner_Base<T>&);
    void SetSolverType(int);
    int GetSolverType() const;
    void SetInvDiagonalJacobi(const Vector<T>&);
    
    size_t GetMemorySize() const;    
    
    template<class Dimension>
    void SelectOptimalIterativeSolver(int N, T a, Dimension& dim, General prop);
    
    template<class Dimension>
    void SelectOptimalIterativeSolver(int N, T a, Dimension& dim, Symmetric prop);
    
    void SetScaling(VectReal_wp& diagonal_left, VectReal_wp& diagonal_right);
    
  };
  
}

#define MONTJOIE_FILE_SOLVE_SYSTEM_HXX
#endif
  
