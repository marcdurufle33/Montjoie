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

class ModeEsPoly2DirichletEigenProblem : public PolynomialSparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, General, ArrayRowSparse>, DistributedMatrix<Complex_wp, General, ArrayRowSparse>>
{

  friend class ModeEs_Solver;

  ModeEs_Solver& var; int nb_rows, nb_dof_E3, nb_dof_H3, nb_dof_all;
  int nb_couches;
  Vector<int> nb_dof_couche, nb_dof_coucheE3;
  Vector<Vector<int> > IndexDofE3, NumGlobE3;
  Vector<Vector<int> > IndexDof, NumGlob;
  Vector<Matrix<Complexe, Prop, Storage> > &vec_Mh, &vec_Mr2h, &vec_Kh, &vec_Ch, &vec_Rjh, &vec_Rih, &vec_Sjh, &vec_Sih, &vec_Dh;
  VectR2& vec_Pos;
  Vector<bool> index_used;

  Vector<int>& DofKeptDir; Vector<int>& IndexNonDirichlet;

  // pour la formulation splittee
  DistributedMatrix<Complex_wp, General, ArrayRowSparse> A0, A1, A2;

private:
  void ComputeA0();
  void ComputeA1();
  void ComputeA2();

public:
  ModeEsPoly2DirichletEigenProblem(int N, ModeEs_Solver& var0,
                                Vector<Matrix<Complex_wp, Prop, Storage> >& Mh,
                                Vector<Matrix<Complex_wp, Prop, Storage> >& Mr2h,
                                Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                                Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                                Vector<Matrix<Complexe, Prop, Storage> >& Rjh,
                                Vector<Matrix<Complexe, Prop, Storage> >& Rih,
                                Vector<Matrix<Complexe, Prop, Storage> >& Sjh,
                                Vector<Matrix<Complexe, Prop, Storage> >& Sih,
                                Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                                VectR2& vec_Pos,
                                Vector<int>& DofDir, Vector<int>& IndexDir);

  void ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X);
  void ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir);
  void RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>&,
                             Matrix<Complexe, Prop, Storage>&);

  void ComputeMatrix();

  void ExtractDirichletAuxVariable(const VectComplex_wp& EHdir, const int& ref, VectComplex_wp& E3, VectComplex_wp& H3, VectComplex_wp& V1, VectComplex_wp& V2, VectComplex_wp& W1, VectComplex_wp& W2);

  void ComputeEsHs(const Complex_wp& L,
                  const Vector<Complex_wp>& eigenvector,
                  Vector<Complex_wp>& Es,
                  Vector<Complex_wp>& Hs,
                  Vector<Complex_wp>& E1,
                  Vector<Complex_wp>& E2,
                  Vector<Complex_wp>& H1,
                  Vector<Complex_wp>& H2);

};

#endif
