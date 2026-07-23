#include "Montjoie.hxx"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;

class ModeEz_Solver;

#ifdef SELDON_WITH_SLEPC
// classe pour resoudre un pb aux valeurs propres non-lineaire avec Slepc
// forme splitte : T(L) = \sum_k f_k(L) A_k
// et f_k est une fraction rationnelle
class ModeEzNonLinearBGTEigenProblem : public NonLinearEigenProblem_Base<Complex_wp>
{
  ModeEz_Solver& var; int nb_rows, nb_dof;
  Vector<Matrix<Complexe, Prop, Storage> >& vec_Mh, &vec_Kh, &vec_Ch;
  Matrix<Complexe, Prop, Storage>& Kh_bord;
  VectReal_wp& Mh_bord;
  Matrix<Complexe, Prop, Storage>& Rh_bord;
  Matrix<Complexe, Prop, Storage> matMh_bord;
  SparseDistributedSolver<Complex_wp> mat_lu;
  // Vector<int>& DofKeptNonDir; Vector<int>& IndexNonDirichlet;

  // pour la formulation splittee
  Matrix<Complexe, Prop, Storage> Mh_sum, Mh_sumH;
  Vector<int> split_reference;

public:
  ModeEzNonLinearBGTEigenProblem(int N, ModeEz_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Matrix<Complexe, Prop, Storage>& Kh_surf,
                              VectReal_wp& Mh_surf,
                              Matrix<Complexe, Prop, Storage>& Rh_surf,
                              bool use_split = false);

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
