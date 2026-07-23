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

class ModeEsPolynomialDirichletEigenProblem : public PolynomialSparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, General, ArrayRowSparse>, DistributedMatrix<Complex_wp, General, ArrayRowSparse>>
{

  friend class ModeEs_Solver;

  ModeEs_Solver& var; int nb_rows, nb_dof_Es, nb_dof_Hs, nb_dof_all;
  Vector<Matrix<Complexe, Prop, Storage> >& vec_Mh, &vec_Kh, &vec_Ch, &vec_Dh, &vec_Th, &vec_Uh;
  Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0, vec_Inv_D_tauBeta, vec_ktrSquared;
  // SparseDistributedSolver<Complex_wp> mat_lu;
  Vector<int>& DofKeptDir; Vector<int>& IndexDirichlet;

  // pour la formulation splittee
  DistributedMatrix<Complex_wp, General, ArrayRowSparse> A0, A1, A2, A3, A4;
  // Vector<int> split_reference;
  bool twisted; // vecteur sur les ref ensuite
  Vector<bool> vec_twisted;
  Complex_wp adim_coef;

private:
  void ComputeInvDBetaMatrix(const Complex_wp& L);
  // void Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M);
  void ComputeA0();
  void ComputeA1();
  void ComputeA2();
  void ComputeA3();
  void ComputeA4();

public:
  ModeEsPolynomialDirichletEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Th,
                              Vector<Matrix<Complexe, Prop, Storage> >& Uh,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& DiagTau0,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& ktr2,
                              Vector<int>& DofDir, Vector<int>& IndexDir);

  void ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X);
  void ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir);
  void RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>&,
                             Matrix<Complexe, Prop, Storage>&);

  void ComputeMatrix();

};

#endif
