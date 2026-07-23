#include "Montjoie.hxx"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;

class ModeEs_Solver;

#ifdef SELDON_WITH_SLEPC
// classe pour resoudre un pb aux valeurs propres non-lineaire avec Slepc
// forme splitte : T(L) = \sum_k f_k(L) A_k
// et f_k est une fraction rationnelle ou non

class ModeEsNonLinear2NDirichletEigenProblem : public NonLinearEigenProblem_Base<Complex_wp>
{
  ModeEs_Solver& var; int nb_rows, nb_dof_E3, nb_dof_H3, nb_dof_all;
  int nb_couches;
  Vector<int> nb_dof_couche, nb_dof_coucheE3;
  Vector<Vector<int> > IndexDofE3, NumGlobE3;
  Vector<Vector<int> > IndexDof, NumGlob;
  Vector<Matrix<Complexe, Prop, Storage>> vec_MhB, vec_KhB, vec_ChB, vec_RjhB, vec_RihB, vec_SjhB, vec_SihB, vec_DhB;
  // Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0, vec_Inv_D_tauBeta, vec_ChgmtVariables, vec_ktrSquared;
  SparseDistributedSolver<Complex_wp> mat_lu;
  Vector<int>& DofKeptDir; Vector<int>& IndexNonDirichlet;
  int countIteration;
  Complex_wp currentL;
  bool precond = false;

  // pour la formulation splittee
  // Matrix<Complexe, Prop, Storage> A1_sum, A1_sumH; // pour la f_i constante (= 1)
  // Vector<int> split_reference;
  Vector<bool> index_used;
//
// private:
//   void Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M);

public:
  ModeEsNonLinear2NDirichletEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<int>& DofDir, Vector<int>& IndexDir,
                              bool use_split = true);

  void ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X);
  void ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir);
  void RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>&,
                             Matrix<Complexe, Prop, Storage>&);

  void ComputeEsHs(const Complex_wp& L,
                   const Vector<Complex_wp>& eigenvector,
                   Vector<Complex_wp>& Es,
                   Vector<Complex_wp>& Hs,
                   Vector<Complex_wp>& E1,
                   Vector<Complex_wp>& E2,
                   Vector<Complex_wp>& H1,
                   Vector<Complex_wp>& H2);

  void ComputeOperator(const Complex_wp& L);

  // calcul de Y = T(L) X
  void MltOperator(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);

  void ComputeOperatorExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A);

  void ComputeJacobian(const Complex_wp& L);

  void MltJacobian(const Complex_wp& L, const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);

  void ComputeJacobianExplicit(const Complex_wp& L, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A);

  // factorisation de T(L)
  void ComputePreconditioning(const Complex_wp& L);
  void ComputePreconditioning(const VectComplex_wp& L, const VectComplex_wp& coef);
  void ComputeSplitPreconditioning(const Vector<int>& numL, const VectComplex_wp& coef);
  void ComputeExplicitPreconditioning(DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A);

  // application du preconditionneur Y = T(L)^{-1} X
  void ApplyPreconditioning(const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);

  void ComputeOperatorSplitExplicit(int i, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& A);
  void MltOperatorSplit(int i, const SeldonTranspose&, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y);

};

#endif
