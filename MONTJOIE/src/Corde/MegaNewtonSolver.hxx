#ifndef MONTJOIE_FILE_MEGA_NEWTON_SOLVER_HXX

namespace Montjoie
{

  MontjoieTimer chrono_string;
  
  //! interface class for linear solvers
  class LinearSolverNewton_Base
  {
  public:
    virtual inline ~LinearSolverNewton_Base() {}
    virtual void Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A) = 0;
    virtual void Factorize(Matrix<Real_wp, General, ArrowCol>& A) = 0;
    virtual void Solve(VectReal_wp& b) = 0;
    
  };
  
  //! solver handling any sparse matrix
  class LinearSolverNewtonSparse : public LinearSolverNewton_Base
  {
    All_MatrixLU<Real_wp> mat_sparse_lu; // solveur direct (par exemple Mumps)
    
  public:
    virtual void Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A);
    virtual void Factorize(Matrix<Real_wp, General, ArrowCol>& A);
    virtual void Solve(VectReal_wp& b);

  };

  //! solver handling arrow matrices
  class LinearSolverNewtonArrow : public LinearSolverNewton_Base
  {
    Matrix<Real_wp, General, ArrowCol> mat_arrow_lu;
    
  public:
    virtual void Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A);
    virtual void Factorize(Matrix<Real_wp, General, ArrowCol>& A);
    virtual void Solve(VectReal_wp& b);

  };


  //! solver handling tiny band matrices
  template<int d>
  class LinearSolverNewtonTinyBand : public LinearSolverNewton_Base
  {
    TinyBandMatrix<Real_wp, d> mat_band_lu;
    
  public:
    virtual void Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A);
    virtual void Factorize(Matrix<Real_wp, General, ArrowCol>& A);
    virtual void Solve(VectReal_wp& b);

  };


  //! solver handling tiny arrow matrices
  template<int d, int m>
  class LinearSolverNewtonTinyArrow : public LinearSolverNewton_Base
  {
    TinyArrowMatrix<Real_wp, d, m> mat_arrow_lu;
    
  public:
    virtual void Factorize(Matrix<Real_wp, General, ArrayRowSparse>& A);
    virtual void Factorize(Matrix<Real_wp, General, ArrowCol>& A);
    virtual void Solve(VectReal_wp& b);

  };


  Matrix<Real_wp, General, ArrowCol> mat_arrow_lu; 
  //! classe implementant l'algorithme de Newton pour
  //! inverser un systeme non-lineaire
  class MegaNewtonSolver
  {
  protected :
    Real_wp tol; // test d'arret pour Newton
    Real_wp tol_convergence; // si |F(x)| < tol_convergence, on considere que l'algo a converge
    int nbitermax; // le nombre maximum d'iterations pour newton
    
    Vector<Real_wp> x_memory, x_prec; // solutions precedentes
// vecteurs intermediaires utilises lors des iterations de Newton
    VectReal_wp scheme, save_scheme, erreur_relative, b, rhs;
    
  public :
    Real_wp t_diff;
    
    // solveur lineaire utilise
    LinearSolverNewton_Base* linear_solver;
    
    // la matrice differentielle;
    Matrix<Real_wp, General, ArrayRowSparse> DiffMatrix;
    // et sa copie non scalée non factorisée
    Matrix<Real_wp, General, ArrayRowSparse> RowDiffMatrix;
    
    // matrice differentielle au format fleche
    Matrix<Real_wp, General, ArrowCol> mat_arrow;
    
    enum {LU_ARROW, LU_SPARSE, LU_TINY_ARROW}; // types de solveurs disponibles
    int type_matrix_lu; // type de solveur utilise (solveur fleche ou solveur direct generaliste)
    // pour un solveur fleche, largeur de bande et nombre de dernieres lignes/colonnes
    int size_band_lu, size_last_row_lu;
    
    VectReal_wp scale_newton; // scaling applique a chaque ligne de l'equation F(x) = 0
    Real_wp erreur; // norme | F(x) | 
    // on normalise le residu avec max(|F(x_0)|, threshold_norme_scheme)
    Real_wp threshold_norme_scheme; 
    bool scale_newton_test; // si oui, on applique un scaling aux equations
    bool force_reevaluation_jacobian; // si oui, la jacobienne est reevaluee a chaque iteration
    ofstream file_out_residual_newton;

    bool scheme_is_linear;
    
  public:
    MegaNewtonSolver();
    ~MegaNewtonSolver();

    VectReal_wp& GetScheme();

    template<class GenericPb>
    void RemplitTest(GenericPb & var);
        
    void SetInputData(const string & description_field, const Vector<string> &parameters);

    template<class GenericPb>
    void ComputeAndFactorizeDiff(GenericPb& var, Matrix<Real_wp>& x, bool compute_df = true);

    template<class GenericPb>
    void UpdateAndFactorizeDiff(GenericPb& var, Matrix<Real_wp>& x);

    void SolveDifferential(VectReal_wp& b, int print_level);

    void ApplyScaling(VectReal_wp& b);

    template<class GenericPb>
    void Init(GenericPb & var, Matrix<Real_wp> & x);

    template<class GenericPb>
    void Solve(GenericPb & var, Matrix<Real_wp> & x);

    template<class GenericPb>
    void Newton(GenericPb & var, Matrix<Real_wp> & x, int &nb,
                const Real_wp& norme_ref, bool reevaluate_jacobian = false);

    void EcritDiffMatrix(Real_wp t, string DOSSIER);
    void EcritScheme(Real_wp t, string DOSSIER);

    template<class GenericPb>
    void CompareDiffMatrix(GenericPb & var,Matrix<Real_wp> & val);
    
  };

}

#define MONTJOIE_FILE_MEGA_NEWTON_SOLVER_HXX
#endif
