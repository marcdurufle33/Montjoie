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

class ModeEsNonLinear_DD_CLA_EigenProblem : public NonLinearEigenProblem_Base<Complex_wp>
{
  ModeEs_Solver& var; int nb_rows, nb_dof, nb_dof_all;
  int nb_couches;
  Vector<int> nb_dof_couche;
  Vector<Vector<int> > IndexDof, NumGlob;
  Vector<Matrix<Complexe, Prop, Storage>> &vec_Mh, &vec_Mr2h, &vec_Kh, &vec_Kr2h, &vec_Ch, &vec_Cr2h, &vec_Rjh, &vec_Rjr2h, &vec_Rih, &vec_Rir2h, &vec_Sjh, &vec_Sjr2h, &vec_Sih, &vec_Sir2h, &vec_Dh, &vec_Dr2h;
  Matrix<Complexe, Prop, Storage> R2;
  VectR2& vec_Pos;
  Matrix<Complexe, Prop, Storage>& Kh_bord;
  VectReal_wp& Mh_bord;
  Matrix<Complexe, Prop, Storage>& Rh_bord;
  Matrix<Complexe, Prop, Storage> matMh_bord;
  // Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0, vec_Inv_D_tauBeta, vec_ChgmtVariables, vec_ktrSquared;
  SparseDistributedSolver<Complex_wp> mat_lu;
  // Vector<int>& DofKeptDir; Vector<int>& IndexNonDirichlet;
  int countIteration;

  // pour la formulation splittee
  // Matrix<Complexe, Prop, Storage> A1_sum, A1_sumH; // pour la f_i constante (= 1)
  Vector<int> split_reference;
  Vector<bool> index_used;

private:
  void ComputeInvDBetaMatrix(const Complex_wp& L);
  void ComputeChgmtVariable(const Complex_wp& L);
  void Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M);

public:
  ModeEsNonLinear_DD_CLA_EigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& Mh,
                              Vector<Matrix<Complex_wp, Prop, Storage> >& Mr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<Matrix<Complexe, Prop, Storage> >& Cr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rjh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rjr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rih,
                              Vector<Matrix<Complexe, Prop, Storage> >& Rir2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sjh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sjr2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sih,
                              Vector<Matrix<Complexe, Prop, Storage> >& Sir2h,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dr2h,
                              VectReal_wp& Mh_surf,
                              Matrix<Complexe, Prop, Storage>& Kh_surf,
                              Matrix<Complexe, Prop, Storage>& Rh_surf,
                              VectR2& vec_Pos,
                              bool use_split = true);
  void ExtractAuxVariable(const VectComplex_wp& EHdir, const int& ref, VectComplex_wp& E3, VectComplex_wp& H3, VectComplex_wp& Fref, VectComplex_wp& Gref);
  void ExpandAuxVariable(const int& ref, const VectComplex_wp& E3, const VectComplex_wp& H3, const VectComplex_wp& Fref, const VectComplex_wp& Gref, VectComplex_wp& EHdir);
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
