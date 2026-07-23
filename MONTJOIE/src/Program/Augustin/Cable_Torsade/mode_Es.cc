#include "Montjoie.hxx"
#include "ModeEsNonLinear_DD_DirichletEigenProblem.hxx"
#include "ModeEsNonLinear_DD_CLA_EigenProblem.hxx"
#include "ModeEsNonLinear_Dk_DirichletEigenProblem.hxx"
#include "ModeEsNonLinear2NDirichletEigenProblem.hxx"
#include "ModeEsPoly1DirichletEigenProblem.hxx"
#include "ModeEsPoly2DirichletEigenProblem.hxx"
#define PATH_AUGUSTIN "./src/Program/Augustin"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;


class ModeEs_Solver : public InputDataProblem_Base,
                      public VarFiniteElement<Dimension2>
{
#ifdef SELDON_WITH_SLEPC
  friend class ModeEsNonLinear_DD_DirichletEigenProblem;
  friend class ModeEsNonLinear_DD_CLA_EigenProblem;
  friend class ModeEsNonLinear_Dk_DirichletEigenProblem;
  friend class ModeEsNonLinear2NDirichletEigenProblem;
  friend class ModeEsPolynomialDirichletEigenProblem;
  friend class ModeEsPoly1DirichletEigenProblem;
  friend class ModeEsPoly2DirichletEigenProblem;
#endif
protected:
  // pulsation omega
  Real_wp omega, torsion;

  // bords Dirichlet
  Vector<int> ref_bord;
  int ref_water;
  Real_wp ray;

  // pour les bords quasi-periodiques
  // TinyVector<int, 2> ref_periodic;
  // int num_mode_periodic;
  // Vector<Vector<bool> > is_dof_with_phase;

  // maillage
  Mesh<Dimension2> mesh; string name_mesh;

  // numerotation
  MeshNumbering<Dimension2> mesh_num;
  Vector<int> IndexDirichlet, DofKeptDir;
  Vector<int> IndexLocalBorder, NumGlobBorder;

  // indices physiques
  Vector<Complex_wp> ref_epsilon, ref_mu, ref_sigma;
  Complex_wp rhoInfini;
  // indices torsions
  Vector<Real_wp> ref_tau;

  // objet pour localiser les points sur le maillage
  GridInterpolation<Dimension2> all_points_display;

  // liste des grilles de sortie
  Vector<GridInterpolationFull<Dimension2> > var_grid;

  // type de solveur valeurs propres
  int type_eigensolver;
  enum {LINEAR_ARPACK, LINEAR_SLEPC, NL_RII, NL_SLP, NL_EIGS, POLYNOMIAL};
  bool use_split_formulation, use_explicit_matrix, use_default_petsc_solver;

  // type de résolution non LINÉAIRE
  int type_non_lin;
  enum {NL_2N, NL_4N_DD, NL_4N_DK};

  // type de condition au bord
  int type_CLA;
  enum {DIRICHLET, SOTOR};

  // type de coordonnées pour les champs transversaux
  int type_coord_other_fields;
  enum {CARTESIAN, POLAR, NONE};

  // shift pour le calcul des valeurs propres
  Complex_wp shift; int nb_asked_eigenval;
  Real_wp radRegion;

  // liste des conducteurs
  Vector<int> ref_conductor, ref_vol_conductor;
  Vector<Complex_wp> phase_conductor;
  Vector<Vector<int> > ddl_conductor;
  VectReal_wp Lk_conductor; Vector<int> num_conduc_ddl;
  VectReal_wp vec_Lk;

public:
  ModeEs_Solver() : mesh_num(mesh)
  { omega = 2.0*pi_wp;
    type_eigensolver = LINEAR_SLEPC; shift = Complex_wp(2.0,0.2);
    type_coord_other_fields = CARTESIAN;
    nb_asked_eigenval = 10;
    use_split_formulation = false; use_explicit_matrix = false;
    use_default_petsc_solver = false;
    torsion = 0.0;
    // type_CLA = DIRICHLET;
  }

  void InitIndices(int ref_max)
  {
    ref_epsilon.Reallocate(ref_max+1);
    ref_mu.Reallocate(ref_max+1);
    ref_sigma.Reallocate(ref_max+1);
    ref_tau.Reallocate(ref_max+1);
    ref_tau.Zero();
  }

  void SetInputData(const string& keyword, const VectString& parameters);

  // version simplifiee de ConstructFiniteElement (dans VarProblem.cxx)
  void ConstructFiniteElement(const string& name_elt);

  // traitement de Dirichlet (et quasi-periodique)
  void RenumberDofs();


  // calcul des matrices elements finis
  void ComputeFemMatrices(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Th,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Uh,
                          VectComplex_wp& vec_Vk,
                          Vector<VectComplex_wp >& vec_Vk2,
                          VectR2& vec_Pos,
                          Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_diag_D_tau0,
                          Vector<Matrix<Real_wp, Symmetric, DiagonalRow>>& vec_diag_D_tau0_Bis
                        );

  void ComputeFemMatricesSurf(int ref_surf,
                          VectReal_wp& Mh_surf,
                          Matrix<Complexe, Prop, Storage>& Kh_surf,
                          Matrix<Complexe, Prop, Storage>& Rh_surf
                        );


  void ComputeFemMatricesLin(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complex_wp, Prop, Storage> >& vec_M1h,
                          Vector<Matrix<Complex_wp, Prop, Storage> >& vec_M2h,
                          Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ti1h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ti2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Tj1h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Tj2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Di1h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Di2h,
                          VectR2& vec_Pos
                        );

  void ComputeFemMatricesPoly2(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                          VectR2& vec_Pos
                        );

  void ComputeFemMatricesNonLin(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                          VectR2& vec_Pos
                        );

  void ComputeFemMatricesNonLin2N(const Complex_wp& L,
                                  Vector<Matrix<Complex_wp, Prop, Storage> >& vec_MhB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_KhB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_ChB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_RjhB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_RihB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_SjhB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_SihB,
                                  Vector<Matrix<Complexe, Prop, Storage> >& vec_DhB
                                );

  void ComputeMassMatrix(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_MhB);

  void ComputeIntensity(const VectComplex_wp& Ez, const VectComplex_wp& vec_Vk, int k)
  {
    VectComplex_wp I0(Lk_conductor.GetM());
    I0.Zero();
    for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
      if (num_conduc_ddl(i) >= 0)
        {
          int n = num_conduc_ddl(i);
          int ref = ref_vol_conductor(n);
          Complex_wp sigma = ref_sigma(ref);
          I0(n) += vec_Vk(i)*Ez(i)*sigma*Lk_conductor(n);
        }

    cout << "Reference conducteurs = " << ref_vol_conductor << endl;
    cout << "Intensity for mode " << k << " = " << I0 << endl;
  }

  // void ComputeTwistedMaterialsDiagonalMatrix(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_diag_tau0,
  //                                             Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0,
  //                                           Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_ktrSquared)
  // {
  //   cout << "Construction des matrices de matériaux torsadés" << endl;
  //   // à refaire intégralement lors de la prise en compte des torsades multiples
  //   vec_D_tau0.Reallocate(this->ref_epsilon.GetM());
  //   vec_ktrSquared.Reallocate(this->ref_epsilon.GetM());
  //   Matrix<Complexe, Symmetric, DiagonalRow> diag_tau0(vec_diag_tau0(0)), D_tau0, ktr2;
  //   int N = diag_tau0.GetM();
  //   D_tau0.Reallocate(N,N);
  //   ktr2.Reallocate(N,N);
  //
  //   Vector<bool> index_used(this->ref_epsilon.GetM());
  //   index_used.Fill(false);
  //   for (int i = 0; i < this->mesh.GetNbElt(); i++)
  //     index_used(this->mesh.Element(i).GetReference()) = true;
  //
  //   for (int ref = 0; ref < index_used.GetM(); ref++)
  //     if (index_used(ref))
  //       {
  //         Complex_wp k2 = Iwp * this->omega * this->ref_mu(ref) * (-Iwp * this->omega * this->ref_epsilon(ref) + this->ref_sigma(ref));
  //         for (int i = 0 ; i < N ; i++)
  //           {
  //             D_tau0.Get(i,i) = -k2 * diag_tau0(i,i);
  //             ktr2.Get(i,i) = k2 * (diag_tau0(i,i) - 1.0);
  //           }
  //
  //         vec_D_tau0(ref) = D_tau0;
  //         vec_ktrSquared(ref) = ktr2;
  //       }
  // }



  void ComputeTwistedMaterialsDiagonalMatrix_Bis(const Vector<Matrix<Real_wp, Symmetric, DiagonalRow>>& vec_diag_tau0_Bis,
                                              Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                                              Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_ktrSquared,
                                              Vector<int>& nbRefNoeud);


  void ComputeMatrixV2ToV1(const Complex_wp& L,
                            const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                            const Vector<int>& nbRefNoeud,
                            Matrix<Complexe, Symmetric, DiagonalRow>& Dtau);

  void ComputeMatrixV2ToEs(const Complex_wp& L,
                            const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                            const Vector<int>& nbRefNoeud,
                            Matrix<Complexe, Symmetric, DiagonalRow>& Dtau2);


  void ComputeEsHs(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                  const Vector<int>& nbRefNoeud,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Uh,
                  const Vector<Matrix<Complexe, Prop, Storage>> &vec_Th,
                  const Vector<Matrix<Complexe, Symmetric, DiagonalRow>> &vec_ktrSquared,
                  const Complex_wp& L,
                  const Vector<Complex_wp>& V2, const Vector<Complex_wp>& W2,
                  Vector<Complex_wp>& Es, Vector<Complex_wp>& Hs);

  void ComputeMatrixF2ToF1(const Complex_wp& L,
                            const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                            const Vector<int>& nbRefNoeud,
                            Matrix<Complexe, Symmetric, DiagonalRow>& Dtau);

  void FromF2G2ToEsHs(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                  const Vector<int>& nbRefNoeud,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Uh,
                  const Vector<Matrix<Complexe, Prop, Storage>> &vec_Th,
                  const Vector<Matrix<Complexe, Symmetric, DiagonalRow>> &vec_ktrSquared,
                  const Complex_wp& L,
                  const Vector<Complex_wp>& F2, const Vector<Complex_wp>& G2,
                  Vector<Complex_wp>& Es, Vector<Complex_wp>& Hs);


  void AssemblageDirichletMatricesVerifDroit(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                          Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch,
                          const Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                          const Vector<Matrix<Complexe, Prop, Storage> >& vec_Th,
                          const Vector<Matrix<Complexe, Prop, Storage> >& vec_Uh,
                          const VectR2 vec_Pos);


  void ResolutionNonLinear_DD_Dirichlet(Vector<Complex_wp>& betas,
                                    Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3,
                                    Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                                    Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                                    VectR2& vec_Pos);

  void ResolutionNonLinear_Dk_Dirichlet(Vector<Complex_wp>& betas,
                                    Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3,
                                    Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                                    Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                                    VectR2& vec_Pos);

  void ResolutionNonLinear2NDirichlet(Vector<Complex_wp>& betas,
                                    Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3);

  void ResolutionNonLinear_DD_CLA(Vector<Complex_wp>& betas,
                                    Vector<VectComplex_wp>& E3, Vector<VectComplex_wp>& H3,
                                    Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                                    Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Kr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Cr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rjr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rih,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Rir2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sjr2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sih,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Sir2h,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Dh,
                                    Vector<Matrix<Complexe, Prop, Storage> >& vec_Dr2h,
                                    VectR2& vec_Pos,
                                    VectReal_wp& Mh_surf,
                                    Matrix<Complexe, Prop, Storage>& Kh_surf,
                                    Matrix<Complexe, Prop, Storage>& Rh_surf);


  // version simplifiee de InitGrid pour localiser les points du maillage sur les sorties
  void InitGrid();

  template<class T>
  void WriteOutputFile(const Vector<T>& x, const Vector<T>& y, int num_mode, Complex_wp beta);

  template<class T>
  void WriteOutputFileScalar(const Vector<T>& u, const string& name);


  // fonction principale pour lancer la simulation
  void RunAll(const string& input_file)
  {
    InitIndices(50);

    // Gauss-Lobatto finite elementFR

    string name_elt = "QUADRANGLE_LOBATTO";

    // on lit le fichier de donnees
    ReadInputFile(input_file, *this);

    // on construit le maillage
    this->mesh.Read(this->name_mesh);

    // on construit l'element fini
    this->ConstructFiniteElement(name_elt);

    // on numerote le maillage
    this->mesh_num.NumberMesh();
    cout << "Number of degrees of freedom = " << this->mesh_num.GetNbDof() << endl;

    // calcul des ddls de Dirichlet (pour les eliminer de la numerotation)
    this->mesh.Write("test.mesh");
    this->RenumberDofs();

    // calcul de la grille d'interpolation (pour sorties des vecteurs propres)
    this->InitGrid();

    // calcul des matrices de rigidite et de la matrice de masse
    Vector<Matrix<Complexe, Prop, Storage> > vec_Mh;
    Vector<Matrix<Complex_wp, Prop, Storage> > vec_M1h;
    Vector<Matrix<Complex_wp, Prop, Storage> > vec_M2h;
    Vector<Matrix<Complex_wp, Prop, Storage> > vec_Mr2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Rih;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Rir2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Sih;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Sir2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Rjh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Rjr2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Sjh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Sjr2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Ti1h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Ti2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Tj1h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Tj2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Di1h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Di2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Kh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Kr2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Ch;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Cr2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Dh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Dr2h;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Th;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Uh;
    VectReal_wp Mh_surf;
    Matrix<Complexe, Prop, Storage> Kh_surf;
    Matrix<Complexe, Prop, Storage> Rh_surf;
    DISP(ref_bord);


    VectComplex_wp vec_Vk;
    Vector<VectComplex_wp > vec_Vk2;
    VectR2 vec_Pos;
    Vector<Matrix<Complexe, Symmetric, DiagonalRow>> vec_diag_D_tau0, vec_D_tau0, vec_D_tau0_Bis, vec_ktrSquared;
    // Vector<Matrix<Complexe, Symmetric, DiagonalRow>> vec_diag_D_tau0, vec_D_tau0_Bis, vec_ktrSquared;
    Vector<Matrix<Real_wp, Symmetric, DiagonalRow>> vec_diag_D_tau0_Bis;
    Vector<int> nbRefNoeud;
    Matrix<Complexe, Symmetric, DiagonalRow> test, Dtau2, DTau;

    // Matrices calculées au début de chaque résolution (car parfois différentes)

    // this->ComputeFemMatrices(vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_Vk, vec_Vk2, vec_Pos, vec_diag_D_tau0, vec_diag_D_tau0_Bis);
    // this->ComputeFemMatricesPoly2(vec_Mh, vec_Mr2h, vec_Kh, vec_Ch, vec_Rjh, vec_Rih, vec_Sjh, vec_Sih, vec_Dh, vec_Pos);

    // DISP(vec_Rjh(1));
    // DISP(vec_Rih(1));

    // this->ComputeTwistedMaterialsDiagonalMatrix(vec_diag_D_tau0, vec_D_tau0, vec_ktrSquared);
    // this->ComputeTwistedMaterialsDiagonalMatrix_Bis(vec_diag_D_tau0_Bis, vec_D_tau0_Bis, vec_ktrSquared, nbRefNoeud);
    // Complex_wp zero(0.002,0.004);

    // for (int i = 0 ; i < 20 ; i++) {
    //   DISP(vec_ktrSquared(1).Get(i,i));
    // }
    int N = this->mesh_num.GetNbDof();
    DISP(N);


    // this->ComputeMatrixF2ToF1(zero, vec_D_tau0_Bis, nbRefNoeud, DTau);
    // this->ComputeMatrixV2ToV1(zero, vec_D_tau0_Bis, nbRefNoeud, DTau);
    // VectComplex_wp ones(N), res(N);
    // res.Zero();
    // for (int i = 0 ; i < N ; i++) {
    //   ones.Get(i) = 1.0;
    // }
    // Mlt(DTau, ones, res);
    // WriteOutputFileScalar(res, "F2ToF1");

    // VectComplex_wp ones(N), res(N);
    // res.Zero();
    // for (int i = 0 ; i < N ; i++) {
    //   ones.Get(i) = 1.0;
    // }
    //
    //
    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    // for (int ref = 0 ; ref < vec_diag_D_tau0_Bis.GetM() ; ref++) {
    //   if (index_used(ref)) {
    //     res.Zero();
    //     Mlt(vec_ktrSquared(ref), ones, res);
    //     WriteOutputFileScalar(res, "KTR2_"+ to_str(ref));
    //     res.Zero();
    //     Mlt(vec_D_tau0_Bis(ref), ones, res);
    //     WriteOutputFileScalar(res, "D_tau0_"+ to_str(ref));
    //   }
    // }
    // abort();


    // Vector<bool> index_used(this->ref_epsilon.GetM());
    // index_used.Fill(false);
    // for (int i = 0; i < this->mesh.GetNbElt(); i++)
    //   index_used(this->mesh.Element(i).GetReference()) = true;
    //
    // for (int i = 0 ; i < N ; i++) {
    //   cout << "Noeud " << i+1 << ", " << vec_Pos(i) << ", ";
    //   for (int ref = 0 ; ref < vec_diag_D_tau0_Bis.GetM() ; ref++) {
    //     if (index_used(ref)) {
    //       cout << "ref " << ref << " : " << vec_D_tau0_Bis(ref)(i,i) << ", ";
    //     }
    //   }
    //   cout << "nbRef " << nbRefNoeud(i);
    //   cout << endl;
    // }
    DISP(ref_water);
    DISP(ray);

    if ((type_eigensolver != LINEAR_ARPACK) && (type_eigensolver != LINEAR_SLEPC) && (type_eigensolver != POLYNOMIAL))
      {
        DISP(this->use_split_formulation);
        this->ComputeFemMatricesNonLin(vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h, vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h, vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos);
        if (type_CLA == SOTOR)
          ComputeFemMatricesSurf(ref_bord(0), Mh_surf, Kh_surf, Rh_surf);

#ifdef SELDON_WITH_SLEPC


        Vector<Complex_wp> betas;
        Vector<VectComplex_wp> E3, H3;

        if (type_non_lin == NL_4N_DD)
          {
            if (type_CLA == DIRICHLET)
              ResolutionNonLinear_DD_Dirichlet(betas, E3, H3, vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h, vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h, vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos);
            else if (type_CLA == SOTOR)
              ResolutionNonLinear_DD_CLA(betas, E3, H3, vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h, vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h, vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos, Mh_surf, Kh_surf, Rh_surf);
            else
              {
                cout << "Not implemented condition" << endl;
                abort();
              }
          }
        else if (type_non_lin == NL_4N_DK)
          {
            if (type_CLA == DIRICHLET)
              ResolutionNonLinear_Dk_Dirichlet(betas, E3, H3, vec_Mh, vec_Mr2h, vec_Kh, vec_Kr2h, vec_Ch, vec_Cr2h, vec_Rjh, vec_Rjr2h, vec_Rih, vec_Rir2h, vec_Sjh, vec_Sjr2h, vec_Sih, vec_Sir2h, vec_Dh, vec_Dr2h, vec_Pos);
            else
              {
                cout << "Not implemented condition" << endl;
                abort();
              }
          }
        else if (type_non_lin == NL_2N)
          {
            if (type_CLA == DIRICHLET)
              ResolutionNonLinear2NDirichlet(betas, E3, H3);
            else
              {
                cout << "Not implemented condition" << endl;
                abort();
              }
          }
        else
          {
            cout << "Not implemented formulation" << endl;
            abort();
          }



        // for (int k = 0; k < betas.GetM(); k++)
        //   {
        //     // this->FromF2G2ToEsHs(vec_D_tau0_Bis, nbRefNoeud, vec_Mh, vec_Uh, vec_Th, vec_ktrSquared, lambda.Get(k), V2, W2, Es, Hs);
        //
        //     // ComputeIntensity(Es, vec_Vk, k);
        //     // WriteOutputFileScalar(V2, "F2_"+ to_str(k));
        //     // WriteOutputFileScalar(W2, "G2_"+ to_str(k));
        //     cout << "VP " << k << ": " << betas(k) <<  endl;
        //     WriteOutputFileScalar(E3(k), "ModeE3_"+ to_str(k));
        //     WriteOutputFileScalar(H3(k), "ModeH3_"+ to_str(k));
        //     // WriteOutputFile(E3(k), H3(k), k, betas(k));
        //   }

#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
        return;
      }
    else if (type_eigensolver == POLYNOMIAL)
      {
        // this->ComputeFemMatricesLin(vec_Mh, vec_M1h, vec_M2h, vec_Mr2h, vec_Rih, vec_Sih, vec_Ti1h, vec_Ti2h, vec_Tj1h, vec_Tj2h, vec_Di1h, vec_Di2h, vec_Pos);
        this->ComputeFemMatricesPoly2(vec_Mh, vec_Mr2h, vec_Kh, vec_Ch, vec_Rjh, vec_Rih, vec_Sjh, vec_Sih, vec_Dh, vec_Pos);
#ifdef SELDON_WITH_SLEPC
        cout << "Résolution polynomiale" << endl;
        DistributedMatrix<Complex_wp, General, ArrayRowSparse> A2, A1, A0;

        // il faut les remplir
        // ModeEsPoly1DirichletEigenProblem var_eig(N, *this, vec_Mh, vec_M1h, vec_M2h, vec_Mr2h, vec_Rih, vec_Sih, vec_Ti1h, vec_Ti2h, vec_Tj1h, vec_Tj2h, vec_Di1h, vec_Di2h, this->DofKeptDir, this->IndexDirichlet);

        ModeEsPoly2DirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Mr2h, vec_Kh, vec_Ch, vec_Rjh, vec_Rih, vec_Sjh, vec_Sih, vec_Dh, vec_Pos, this->DofKeptDir, this->IndexDirichlet);

        // DISP(vec_Pos);

        A0 = var_eig.A0;
        A1 = var_eig.A1;
        A2 = var_eig.A2;
        // A3 = var_eig.A3;
        // A4 = var_eig.A4;

        Vector<DistributedMatrix<Complex_wp, General, ArrayRowSparse>* > list_op(2);
        list_op(0) = &A0;
        list_op(1) = &A1;
        // list_op(2) = &A2;
        // list_op(3) = &A3;


        var_eig.SetStoppingCriterion(1e-12);
        var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);
        var_eig.InitMatrix(list_op, A2); // celle de plus grand degré est à part
        var_eig.SetPrintLevel(6);

        var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
        var_eig.SetSpectralTransformation(true);

        Vector<Complex_wp> lambda, lambda_imag;
        Matrix<Complex_wp, General, ColMajor> eigen_vec;

        FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

        // DISP(lambda); DISP(lambda_imag);

        for (int k = 0; k < eigen_vec.GetN(); k++)
          {
            int nb_dof_E3 = DofKeptDir.GetM();
            int nb_dof_H3 = N;
            VectComplex_wp x(2*nb_dof_H3), X;
            VectComplex_wp y;
            GetCol(eigen_vec, k, X);
            int count13 = 0;
            int count14 = 0;
            int count15 = 0;
            for (int i = 0 ; i < X.GetM() ; i++)
              {
                if (abs(X(i)) <= 10e-13)
                  {
                    count13 ++;
                  }
                if (abs(X(i)) <= 10e-14)
                  {
                    count14 ++;
                  }
                if (abs(X(i)) <= 10e-15)
                  {
                    count15 ++;
                  }
              }
            DISP(X.GetM());
            DISP(count13);
            DISP(count14);
            DISP(count15);

            VectComplex_wp Es, Hs, E1, E2, H1, H2;
            var_eig.ComputeEsHs(lambda(k), X, Es, Hs, E1, E2, H1, H2);

            X.Resize(nb_dof_E3 + nb_dof_H3);
            VectComplex_wp Y(nb_dof_E3 + nb_dof_H3); Y.Zero();
            ModeEsNonLinear2NDirichletEigenProblem var_eig_model(N, *this,
                                                this->DofKeptDir, this->IndexDirichlet,
                                                false);

            var_eig_model.MltOperator(lambda(k), SeldonNoTrans, X, Y);
            cout << "VP " << k << ": " << lambda(k) << " ; ||T(VP) X|| : " << Norm2(Y) << " ; erreur relative : " << Norm2(Y) / Norm2(X) << endl;

            // this->ComputeMatrixV2ToEs(lambda.Get(k), vec_D_tau0_Bis, nbRefNoeud, Dtau2);
            var_eig.ExtractDirichlet(X, x);

            VectComplex_wp E3(nb_dof_H3); E3.Zero();
            VectComplex_wp H3(nb_dof_H3); H3.Zero();
            for (int i = 0; i < nb_dof_H3; i++)
              {
                E3(i) = x(i);
                H3(i) = x(nb_dof_H3 + i);
              }
            //x.Write("Ez_Hz0.dat");

            // this->ComputeEsHs(vec_D_tau0_Bis, nbRefNoeud, vec_Mh, vec_Uh, vec_Th, vec_ktrSquared, lambda.Get(k), V2, W2, Es, Hs);

            // ComputeIntensity(Es, vec_Vk, k);
            // WriteOutputFileScalar(V2, "V2_"+ to_str(k));
            // WriteOutputFileScalar(W2, "W2_"+ to_str(k));
            cout << "VP " << k << ": " << lambda(k) <<  endl;
            // WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
            // WriteOutputFileScalar(Hs, "ModeHs_"+ to_str(k));
            WriteOutputFileScalar(E3, "ModeE3_"+ to_str(k));
            WriteOutputFileScalar(H3, "ModeH3_"+ to_str(k));
            WriteOutputFileScalar(Es, "ModeEs_"+ to_str(k));
            WriteOutputFileScalar(Hs, "ModeHs_"+ to_str(k));
            // WriteOutputFile(E3, H3, k, lambda(k));
          }
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
        return;
      }
    else
      {
        // SYSTÈME LINÉAIRE

        int nb_dof_Es = DofKeptDir.GetM();
        int nb_dof_Hs = N;

        // Définitions des matrices
        Matrix<Complex_wp, General, ArrayRowSparse> K_tilde;
        Matrix<Complex_wp, General, ArrayRowSparse> M_tilde;

        // if (type_CLA == DIRICHLET)
        if (true)
          {
            AssemblageDirichletMatricesVerifDroit(K_tilde, M_tilde, vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_Pos);
          }

        cout << "Matrices créées" << endl;
        M_tilde.WriteText("Mt.dat");
        K_tilde.WriteText("Kt.dat");

        SparseEigenProblem<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>,
                           Matrix<Complex_wp, General, ArrayRowSparse> > var_eig;

        var_eig.SetStoppingCriterion(1e-12);
        var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);
        var_eig.SetComputationalMode(var_eig.INVERT_MODE);

        var_eig.InitMatrix(K_tilde, M_tilde);
        var_eig.SetPrintLevel(2);

        var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, this->shift, var_eig.SORTED_MODULUS);

        cout << "On calcule les valeurs propres/vecteurs propres" << endl;
        VectComplex_wp lambda, lambda_imag;
        Matrix<Complex_wp, General, ColMajor> eigen_vec;
        if (type_eigensolver == LINEAR_ARPACK)
          {
            GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);
          }
        else if (type_eigensolver == LINEAR_SLEPC)
         {
           GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec,
                                      TypeEigenvalueSolver::SLEPC);
         }


        VectComplex_wp betas;
        betas.Reallocate(lambda.GetM());
        Copy(lambda, betas);

        DISP(type_coord_other_fields);

        VectComplex_wp x;
        for (int i = 0; i < betas.GetM(); i++)
          {
            GetCol(eigen_vec, i, x);
            VectComplex_wp Ez;
            VectComplex_wp Hz; // on ne prend que les données de Ez et de Hz
            VectComplex_wp X = x; X.Resize(nb_dof_Es + nb_dof_Hs);
            VectComplex_wp Y = x; Y.Resize(nb_dof_Es + nb_dof_Hs);
            Ez.Reallocate(nb_dof_Hs); Ez.Zero();
            Hz.Reallocate(nb_dof_Hs); Hz.Zero();
    #ifdef SELDON_WITH_SLEPC
            // vérification (norme de T(L) X)
            // ModeEsNonLinear_DD_DirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_D_tau0_Bis, vec_ktrSquared,
            //                                     this->DofKeptDir, this->IndexDirichlet,
            //                                     this->use_split_formulation);
            //
            // // var_eig.MltOperator(betas(i), SeldonNoTrans, X, Y);
            // cout << "VP " << i << ": " << betas(i) << endl;
            //
            // x.Reallocate(2*nb_dof_Hs);
            // var_eig.ExtractDirichlet(X, x);
            // for (int i = 0; i < nb_dof_Hs; i++)
            //   {
            //     Ez(i) = x(i);
            //     Hz(i) = x(nb_dof_Hs + i);
            //   }
            // ComputeIntensity(Ez, vec_Vk, i);
            //
            // WriteOutputFile(Ez, Hz, i, betas(i));
    #endif
          }
      }
  }
};

#ifdef SELDON_WITH_SLEPC
// on inclut les fonctions definies dans le cxx
#include "ResolutionNonLinear.cxx"
#include "ModeEsNonLinear_DD_DirichletEigenProblem.cxx"
#include "ModeEsNonLinear_DD_CLA_EigenProblem.cxx"
#include "ModeEsNonLinear_Dk_DirichletEigenProblem.cxx"
#include "ModeEsNonLinear2NDirichletEigenProblem.cxx"
#include "ModeEsPoly1DirichletEigenProblem.cxx"
#include "ModeEsPoly2DirichletEigenProblem.cxx"
#endif
#include "IOSolution.cxx"
#include "Mesh.cxx"
#include "AssemblageMatricesDirichlet.cxx"
// #include "AssemblageMatrices.cxx"
#include "FiniteElement.cxx"
#include "Torsion.cxx"

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  if (argc < 2)
    {
      cout << "Data file not given" << endl;
      cout << "Usage : ./mode_Es.x data_file.ini" << endl;
      abort();
    }

  // on appelle RunAll avec le fichier de donnees
  ModeEs_Solver var;
  string input_file(argv[1]);
  var.RunAll(PATH_AUGUSTIN + to_str("/Datas/") + input_file);

  return FinalizeMontjoie();
}
