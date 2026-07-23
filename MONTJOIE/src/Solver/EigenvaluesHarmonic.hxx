#ifndef MONTJOIE_FILE_EIGENVALUES_HARMONIC_HXX

namespace Montjoie
{
  // file used to interface eigenvalue solvers
  // with finite element matrices defined in Montjoie
  
  // for definition of methods, see file Eigenvalues_Harmonic.cxx
  
  ////////////////////////////
  // EIGENVALUE COMPUTATION //
  ////////////////////////////
  
  
  //! Class to solve an eigenvalue problem issued from finite element formulation
  /*!
    For finite element method, usually the eigenvalue problem reads as
    find \f$ (\omega, U_h) \mbox{ non-null so that } (-\omega^2 D_h + K_h)\, U_h = 0 \f$
    where D_h is the mass matrix and K_h the stiffness matrix
    For LDG formulation, the eigenvalue problem reads as a first order problem
    find \f$ (\omega, U_h) \mbox{ non-null so that } (-i \omega D_h + K_h)U_h = 0 \f$
    This class specifies some parameters associated to the seek of eigenvalues (omega)
    and eigenvectors (U_h), method used to find the eigenvalues, parameters of the method
  */
  template<class T>
  class EigenProblemMontjoie
    : public VirtualEigenProblem<T, T, T>,
      public InputDataProblem_Base
  {
  protected :
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;

    string file_name_eigenval;
    IVect local_col_numbers;
    
    //! linear solver
    All_LinearSolver& glob_solver;
    
    FemMatrixFreeClass_Base<T>* mat_mass;
    FemMatrixFreeClass_Base<T>* mat_stiff;
    GlobalGenericMatrix<Real_wp> nature_matrix_real;
    GlobalGenericMatrix<Complex_wp> nature_matrix_complex;

    //! Cholesky factorisation of mass matrix if required
    SparseCholeskySolver<Treal> chol_facto_mass_matrix;
    int type_solver_chol;
    
    //! temporary vectors for Cholesky
    Vector<Treal> Xchol_real, Xchol_imag;

    //! eigenvalue solver used
    int eigenvalue_solver_used;
    
  private:
    VarProblem_Base& var_problem;
    VarComputationProblem& var_computation;
    VarBoundaryCondition_Base& var_boundary;
    DistributedProblem_Base& var_comm;
    VarOutputProblem_Base& var_output;

    void InitDefaultValues();
    void InitComputation();

  public:
    template<class TypeEquation>
    EigenProblemMontjoie(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver);

    size_t GetMemorySize() const;
    const string& GetEigenvalFileName() const { return file_name_eigenval; }
    
    void InitMatrix(FemMatrixFreeClass_Base<T>& K, FemMatrixFreeClass_Base<T>& M);
    
    void SetInputData(const string& description_field,
                      const Vector<string>& parameters);
    
    const IVect& GetLocalColumnNumbers() const { return local_col_numbers; }

    
    // mass matrix stuff
    void ComputeDiagonalMass();    

    void ComputeMassForCholesky();
    void FactorizeCholeskyMass();
    void SetCholeskySolver(int solver);
    void GetFactoCholeskyMass(const Real_wp&);
    void GetFactoCholeskyMass(const Complex_wp&);
    
    template<class T0>
    void ExpandVector(const Vector<T0>& X, Vector<T0>& Y, bool assemble = false);

    template<class T0>
    void ContractVector(const Vector<T0>& X, Vector<T0>& Y);

    void SolveCholeskyMass(const SeldonTranspose&, Vector<Real_wp>&);
    void SolveCholeskyMass(const SeldonTranspose&, Vector<Complex_wp>&);

    void MltCholeskyMass(const SeldonTranspose&, Vector<Real_wp>&);
    void MltCholeskyMass(const SeldonTranspose&, Vector<Complex_wp>&);

    void ComputeMassMatrix();
    
    template<class T0>
    void MltMassGen(const SeldonTranspose&, const Vector<T0>& X, Vector<T0>& Y);

    void MltMass(const Vector<Treal>& X, Vector<Treal>& Y);
    void MltMass(const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void MltMass(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);
    void MltMass(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    // stiffness matrix
    void ComputeStiffnessMatrix();
    void ComputeStiffnessMatrix(const T& a, const T& b);    
    
    template<class T0>
    void ComputeStiffnessMatrix(const T0& a, const T0& b,
				FemMatrixFreeClass_Base<T0>& A, bool iterative);

    template<class T0>
    void MltStiffnessGen(const SeldonTranspose&, const Vector<T0>& X, Vector<T0>& Y);
    
    void MltStiffness(const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void MltStiffness(const T& a, const T& b,
                      const Vector<Treal>& X, Vector<Treal>& Y);

    void MltStiffness(const T& a, const T& b,
                      const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void MltStiffness(const SeldonTranspose&, const Vector<Treal>& X, Vector<Treal>& Y);
    void MltStiffness(const SeldonTranspose&, const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    // factorization of K + sigma M
    void ComputeAndFactorizeStiffnessMatrix(const Real_wp& a, const Real_wp& b,
					    int which = EigenProblem_Base<T>::COMPLEX_PART);

    void ComputeAndFactorizeStiffnessMatrix(const Complex_wp& a, const Complex_wp& b,
                                            int which = EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeAndFactoRealMatrix(const Real_wp&, const Real_wp& a, const Real_wp& b, int which);
    void ComputeAndFactoRealMatrix(const Complex_wp&, const Real_wp& a, const Real_wp& b, int which);

    void ComputeSolution(const Vector<Real_wp>& X, Vector<Real_wp>& Y);
    void ComputeSolution(const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);

    void ComputeSolution(const SeldonTranspose&,
			 const Vector<Real_wp>& X, Vector<Real_wp>& Y);
    
    void ComputeSolution(const SeldonTranspose&,
			 const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);


    template<class Complexe>
    void Solve(const SeldonTranspose& transA,
	       GlobalGenericMatrix<Complexe>& nat_mat, Vector<Complexe>& Y);
    
    void Clear();

    // main method called to retrieve all the eigenvalues and eigenvectors
    void ComputeEigenModes();
    
    bool IsSymmetricProblem() const;
    bool IsHermitianProblem() const;
    
  };

  template<class T>
  class PolynomialEigenProblemMontjoie
    : public PolynomialEigenProblem_Base<T>,
      public InputDataProblem_Base
  {
  protected:
    string file_name_eigenval;
    IVect local_col_numbers;
    
    //! linear solver
    All_LinearSolver& glob_solver;
    GlobalGenericMatrix<T> nat_mat_solver;
    
    // finite element matrices
    Vector<FemMatrixFreeClass_Base<T>* > mat_op;
    
  private:
    VarProblem_Base& var_problem;
    VarComputationProblem& var_computation;
    VarBoundaryCondition_Base& var_boundary;
    DistributedProblem_Base& var_comm;
    VarOutputProblem_Base& var_output;
    
  public:
    template<class TypeEquation>
    PolynomialEigenProblemMontjoie(EllipticProblem<TypeEquation>& var,
                                   All_LinearSolver& solver);
    
    void SetInputData(const string& description_field,
                      const Vector<string>& parameters);
    
    const IVect& GetLocalColumnNumbers() const;
    
    void InitComputation();
    
    void ExpandVector(const Vector<T>& X, Vector<T>& Y, bool assemble = false);
    void ContractVector(const Vector<T>& X, Vector<T>& Y);
    
    void ComputeOperator(int num, const Vector<T>& coef);
    void MltOperator(int num, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y);
    
    void FactorizeMass();
    void SolveMass(const SeldonTranspose& trans, const Vector<T>& x, Vector<T>& y);

    void FactorizeOperator(const Vector<T>& coef);
    void SolveOperator(const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y);
    
    // main method called to retrieve all the eigenvalues and eigenvectors
    void ComputeEigenModes();
    
  };

  template<class T>
  Real_wp GetSpectralRadius(VirtualMatrix<T>& mat_stiff, VirtualMatrix<Real_wp>& mat_mass);

  template<class T>
  Real_wp GetSpectralRadius(VirtualMatrix<T>& mat_stiff);
  
} // end namespace

#define MONTJOIE_FILE_EIGENVALUES_HARMONIC_HXX
#endif
