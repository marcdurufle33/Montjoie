#ifndef MONTJOIE_FILE_PRECOND_MAXWELL_HXX

namespace Montjoie
{
  
  //! base class for Helmholtz decomposition of Maxwell's equations
  template<class T>
  class PreconditionerHelmholtzMaxwell
  {
  protected :
    //! Jacobi or SOR smoother
    int type_smoother;
    
    //! matrix used for SOR iteration (symmetric case)
    typename GenericStorage<T>::SparseSymMatrix* mat_ssor_sym, mat_stored_sym;
    //! matrix used for SOR iteration (non-symmetric case)
    typename GenericStorage<T>::SparseMatrix* mat_ssor_unsym, mat_stored_unsym;

    //! matrix used for SOR iteration (symmetric case)
    typename GenericStorage<T>::SparseSymMatrixCSR* matCSR_ssor_sym, matCSR_stored_sym;
    //! matrix used for SOR iteration (non-symmetric case)
    typename GenericStorage<T>::SparseMatrixCSR* matCSR_ssor_unsym, matCSR_stored_unsym;
    
    //! relaxation parameter
    Real_wp omega_ssor;
    
    //! projection from edge elements to nodal elements
    Matrix<Real_wp, General, RowSparse> Projection_NodalToEdge;
    //! finite element matrix associated with nodal elements
    typename GenericStorage<T>::SparseSymMatrix mat_galerkin_H1;
    //! diagonal of mat_galerkin_H1
    Vector<T> diagonal_galerkin_H1;
    //! diagonal for edge elements
    Vector<T> diagonal_jacobi;
    //! Dirichlet dofs for nodal elements and edge elements
    IVect DirichletHelmholtz_dof, DirichletMaxwell_dof;    

  public :
    enum{JACOBI, SSOR};

    PreconditionerHelmholtzMaxwell();
    
    size_t GetMemorySize() const;
    
    void SetParameterRelaxation(const Real_wp& w);
    
    void ConstructPreconditioner(int type, const VarComputationProblem&,
				 const FemMatrixFreeClass_Base<T>& A,
				 const GlobalGenericMatrix<T>& nature, bool same_nature);
    
    void ComputeHelmPrecond(const VarComputationProblem& var,
			    typename GenericStorage<T>::SparseSymMatrix& mat_gal,
			    Matrix<Real_wp, General, RowSparse>& mat_csr_proj,
			    Vector<T>& diag_gal, const GlobalGenericMatrix<T>& nature);

    template<class T0>
    void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<T0>& r, Vector<T0>& z);

    template<class T0>
    void SolveHelmholtzDecomposition(const VirtualMatrix<T>& A, const Vector<T0>& r,
                                     Vector<T0>& z, bool init_guess_null);

    template<class T0>
    void ApplyPreSmoothing(const VirtualMatrix<T>& A, const Vector<T0>& r,
                           Vector<T0>& z, bool init_guess_null);

    template<class T0>
    void ApplyPostSmoothing(const VirtualMatrix<T>& A, const Vector<T0>& r,
                            Vector<T0>& z, bool init_guess_null);
    
  };

  
  //! class for low-order preconditioning and Maxwell's equations
  template<class T>
  class LowOrderPreconditioning_Maxwell
    : public LowOrderPreconditioning_Dim<T, Dimension3>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    
  public :
    template<class TypeEquation>
    LowOrderPreconditioning_Maxwell(EllipticProblem<TypeEquation>& var);
    
    void GetLocalPermutation_LowOrder(const ElementReference<Dimension3, 2>& FaceFine,
				      const ElementReference<Dimension3, 2>& FaceCoarse,
                                      const VectReal_wp& Points,
                                      IVect& num_element, IVect& num_dof_local,
                                      VectReal_wp& LocalScaling);
    
    void ComputeLocalProlongation(const VectReal_wp& Points,
                                  DistributedProblem<Dimension3>& var,
                                  DistributedProblem<Dimension3>& var_low);

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x);
    
    void Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x);

    void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x);
    void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x);

  };
  

  //! low-order preconditioning for Maxwell's equations
  template<>
  class LowOrderPreconditioning<Complex_wp, HarmonicMaxwellEquation_3D>
    : public LowOrderPreconditioning_Maxwell<Complex_wp>
  {
  public :
    template<class TypeEquation>
    LowOrderPreconditioning(EllipticProblem<TypeEquation>& var);
    
  };
  
  
  //! class for multigrid preconditioning and Maxwell's equations
  template<class T>
  class MultigridPreconditioning_Maxwell
    : public MultigridPreconditioning_Dim<T, Dimension3>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    
  protected :
    //! smoother based on Helmholtz decomposition
    PreconditionerHelmholtzMaxwell<T> smoother_helm;

    // type of smoothers
    enum {SSOR, JACOBI, BLOCK_JACOBI, HELMHOLTZ_SSOR, HELMHOLTZ_JACOBI};

  private:
    VarProblem<Dimension3>& var_fine;
    
  public:
    template<class TypeEquation>
    MultigridPreconditioning_Maxwell(EllipticProblem<TypeEquation>& var);

    size_t GetMemorySize() const;

    void SetInputData(const string& keyword, const VectString& param);
    
    void ConstructSmoother(const FemMatrixFreeClass_Base<T>& A,
			   const GlobalGenericMatrix<T>& nature);

    template<class T0>
    void ApplyPreSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			     Vector<T0>& r, bool init_guess_null);

    template<class T0>
    void ApplyPostSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			      Vector<T0>& r, bool init_guess_null);
    
    void ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
                          Vector<Treal>& r, bool init_guess_null);

    void ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
                           Vector<Tcplx>& r, bool init_guess_null);

    void ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
			   Vector<Treal>& r, bool init_guess_null);
    
    void ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
                           Vector<Tcplx>& r, bool init_guess_null);

  };


  //! low-order preconditioning for Maxwell's equations
  template<>
  class MultigridPreconditioning<Complex_wp, HarmonicMaxwellEquation_3D>
    : public MultigridPreconditioning_Maxwell<Complex_wp>
  {
  public :
    template<class TypeEquation>
    MultigridPreconditioning(EllipticProblem<TypeEquation>& var);
    
  };
  
  
  //! base class for preconditioning for Maxwell's equations
  template<class T>
  class All_Preconditioner_Maxwell : public All_Preconditioner_Base<T>
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

  public :
    enum {IDENTITY, SSOR, DIRECT, JACOBI, LOW_ORDER, MULTIGRID, SUBDOMAIN,
          HELMHOLTZ_SSOR, HELMHOLTZ_JACOBI};
    
  protected :
    //! preconditioning based on Helmholtz decomposition
    PreconditionerHelmholtzMaxwell<T> prec_helm;
    
  public :    
    template<class TypeEquation>    
    All_Preconditioner_Maxwell(EllipticProblem<TypeEquation>& var);
    
    size_t GetMemorySize() const;

    void SetInputPreconditioning(const string& keyword, const Vector<string>& parameters);

    void ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
                                 const FemMatrixFreeClass_Base<T>& A,
				 const GlobalGenericMatrix<T>& nature,
				 CondensationBlockSolver_Fem<T>& cond_solver);
    
    void Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z);

    template<class T0>
    void SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
		  const Vector<T0>& b, Vector<T0>& x);

    void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x);
    void Solve(const SeldonTranspose&, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x);

  };


  //! preconditioning for Maxwell's equations
  template<class T>
  class All_Preconditioner<T, HarmonicMaxwellEquation_3D>
    : public All_Preconditioner_Maxwell<T>
  {
  public :
    template<class TypeEquation>    
    All_Preconditioner(EllipticProblem<TypeEquation>& var);
    
  };

}

#define MONTJOIE_FILE_PRECOND_MAXWELL_HXX
#endif
