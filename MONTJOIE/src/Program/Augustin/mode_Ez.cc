#include "Montjoie.hxx"
#include "ModeEzNonLinearDirichletEigenProblem.hxx"
#include "ModeEzNonLinearBGTEigenProblem.hxx"
#define PATH_AUGUSTIN "./src/Program/Augustin"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;

class ModeEz_Solver : public InputDataProblem_Base,
                      public VarFiniteElement<Dimension2>
{
#ifdef SELDON_WITH_SLEPC
  friend class ModeEzNonLinearDirichletEigenProblem;
  friend class ModeEzNonLinearBGTEigenProblem;
#endif

protected:
  // pulsation omega
  Real_wp omega;

  // bords Dirichlet
  Vector<int> ref_bord;
  int ref_water;
  Real_wp ray;

  // pour les bords quasi-periodiques
  TinyVector<int, 2> ref_periodic;
  int num_mode_periodic;
  Vector<Vector<bool> > is_dof_with_phase;

  // maillage
  Mesh<Dimension2> mesh; string name_mesh;

  // numerotation
  MeshNumbering<Dimension2> mesh_num;
  Vector<int> IndexNonDirichlet, DofKeptNonDir;
  Vector<int> IndexLocalBorder, NumGlobBorder;

  // indices physiques
  Vector<Complex_wp> ref_epsilon, ref_mu, ref_sigma;
  Complex_wp rhoInfini;

  bool only_Ez = false;

  // type de condition au bord
  int type_CLA;
  enum {DIRICHLET, NEUMANN, BGT1, BGT2};
  bool linear;
  int linearization_step;

  // calcul valeurs propres laplacien
  bool compute_laplace_eigen;

  // objet pour localiser les points sur le maillage
  GridInterpolation<Dimension2> all_points_display;

  // liste des grilles de sortie
  Vector<GridInterpolationFull<Dimension2> > var_grid;

  // type de solveur valeurs propres
  int type_eigensolver;
  enum {LINEAR_ARPACK, LINEAR_SLEPC, NL_RII, NL_SLP, NL_EIGS, DIRECT};
  bool use_split_formulation, use_explicit_matrix, use_default_petsc_solver;

  // type de coordonnées pour les champs transversaux
  int type_coord_other_fields;
  enum {CARTESIAN, POLAR, NONE};

  // shift pour le calcul des valeurs propres
  Complex_wp shift; int nb_asked_eigenval;
  Real_wp minRe, maxRe, minIm, maxIm;

  // liste des conducteurs
  Vector<int> ref_conductor, ref_vol_conductor;
  Vector<int> ref_vol_isolant;
  Vector<Complex_wp> phase_conductor;
  Vector<Vector<int> > ddl_conductor;
  VectReal_wp Lk_conductor; Vector<int> num_conduc_ddl;
  VectReal_wp vec_Lk;

public:
  ModeEz_Solver() : mesh_num(mesh)
  { omega = 2.0*pi_wp; compute_laplace_eigen = false;
    type_eigensolver = LINEAR_ARPACK; shift = Complex_wp(1.6,1.2);
    type_coord_other_fields = NONE;
    num_mode_periodic = 1; ref_periodic.Init(-1, -1);
    nb_asked_eigenval = 10;
    use_split_formulation = false; use_explicit_matrix = false;
    use_default_petsc_solver = false;
  }

  void InitIndices(int ref_max)
  {
    ref_epsilon.Reallocate(ref_max+1);
    ref_mu.Reallocate(ref_max+1);
    ref_sigma.Reallocate(ref_max+1);
  }

  void SetInputData(const string& keyword, const VectString& parameters);

  // version simplifiee de ConstructFiniteElement (dans VarProblem.cxx)
  void ConstructFiniteElement(const string& name_elt);

  // traitement de Dirichlet, quasi-periodique et ddl conductors
  // À COMMENTER POUR BIEN COMPRENDRE
  void RenumberDofs();

  // calcul des matrices elements finis
  // À COMMENTER POUR BIEN COMPRENDRE
  void ComputeFemMatrices(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                          VectComplex_wp& vec_Vk,
                          Vector<VectComplex_wp >& vec_Vk2,
                          int ref_surf,
                          VectReal_wp& Mh_surf,
                          Matrix<Complexe, Prop, Storage>& Kh_surf,
                          Matrix<Complexe, Prop, Storage>& Rh_surf);

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

  // version simplifiee de InitGrid pour localiser les points du maillage sur les sorties
  void InitGrid();

  template<class T>
  void WriteOutputFile(const Vector<T>& x, const Vector<T>& y, int num_mode, Complex_wp beta);

  template<class T>
  void WriteOutputFileScalar(const Vector<T>& u, const string& name);

  void AssemblageDirichletMatricesOnlyEz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                          Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh);

  void AssemblageDirichletMatricesEzHz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                          Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch);

  void AssemblageMatricesOnlyEz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                          Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                          const VectReal_wp& Mh_surf,
                          const Matrix<Complexe, Prop, Storage>& Kh_surf);

  void AssemblageMatricesEzHz(Matrix<Complex_wp, General, ArrayRowSparse>& K_tilde,
                          Matrix<Complex_wp, General, ArrayRowSparse>& M_tilde,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                          const Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch,
                          const Matrix<Complexe, Prop, Storage>& Rh_surf,
                          const VectReal_wp& Mh_surf,
                          const Matrix<Complexe, Prop, Storage>& Kh_surf);

  void ResolutionNonLinearDirichlet(Vector<Complex_wp>& betas,
                                    Vector<VectComplex_wp>& Ez, Vector<VectComplex_wp>& Hz,
                                    Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                                    Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                                    Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch);



  void ResolutionNonLinearBGT(Vector<Complex_wp>& betas,
                              Vector<VectComplex_wp>& Ez, Vector<VectComplex_wp>& Hz,
                              Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                              Vector<Matrix<Complexe, Prop, Storage>>& vec_Kh,
                              Vector<Matrix<Complexe, Prop, Storage>>& vec_Ch,
                              VectReal_wp& Mh_surf,
                              Matrix<Complexe, Prop, Storage>& Kh_surf,
                              Matrix<Complexe, Prop, Storage>& Rh_surf);


  // fonction principale pour lancer la simulation
  void RunAll(const string& input_file)
  {
    InitIndices(50);

    // Gauss-Lobatto finite element
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
    DISP(ref_water);

    // calcul de la grille d'interpolation (pour sorties des vecteurs propres)
    this->InitGrid();

    // calcul des matrices de rigidite et de la matrice de masse
    Vector<Matrix<Complexe, Prop, Storage> > vec_Mh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Kh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Ch;
    VectComplex_wp vec_Vk;
    Vector<VectComplex_wp > vec_Vk2;
    VectReal_wp Mh_surf;
    Matrix<Complexe, Prop, Storage> Kh_surf;
    Matrix<Complexe, Prop, Storage> Rh_surf;
    DISP(ref_bord);
    this->ComputeFemMatrices(vec_Mh, vec_Kh, vec_Ch, vec_Vk, vec_Vk2, ref_bord(0), Mh_surf, Kh_surf, Rh_surf);

    DISP(this->omega);

    int N = this->mesh_num.GetNbDof();

    if (!linear)
      {
        DISP(this->use_split_formulation);
#ifdef SELDON_WITH_SLEPC

        Vector<Complex_wp> betas;
        Vector<VectComplex_wp> Ez, Hz;

        if (type_CLA == DIRICHLET)
          {
            ResolutionNonLinearDirichlet(betas, Ez, Hz, vec_Mh, vec_Kh, vec_Ch);
          }
        else if ((type_CLA == BGT1) || (type_CLA == BGT2))
          {
            ResolutionNonLinearBGT(betas, Ez, Hz, vec_Mh, vec_Kh, vec_Ch, Mh_surf, Kh_surf, Rh_surf);
          }
        else if (type_CLA == NEUMANN)
          {
            cout << "Not a linear problem" << endl;
            abort();
          }

        for (int k = 0; k < betas.GetM(); k++)
          {
            ComputeIntensity(Ez(k), vec_Vk, k);
            WriteOutputFile(Ez(k), Hz(k), k, betas(k));
          }
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
      }
    else
      {
        int nb_dof_Ez = DofKeptNonDir.GetM();
        int nb_dof_Hz = N;

        // Définitions des matrices
        Matrix<Complex_wp, General, ArrayRowSparse> K_tilde;
        Matrix<Complex_wp, General, ArrayRowSparse> M_tilde;

        if (type_CLA == DIRICHLET)
          {
            if (only_Ez)
              {
                AssemblageDirichletMatricesOnlyEz(K_tilde, M_tilde, vec_Mh, vec_Kh);
              }
            else
              {
                AssemblageDirichletMatricesEzHz(K_tilde, M_tilde, vec_Mh, vec_Kh, vec_Ch);
              }
          }
        else
          {
            if (only_Ez)
              {
                AssemblageMatricesOnlyEz(K_tilde, M_tilde, vec_Mh, vec_Kh, Mh_surf, Kh_surf);
              }
            else
              {
                AssemblageMatricesEzHz(K_tilde, M_tilde, vec_Mh, vec_Kh, vec_Ch, Rh_surf, Mh_surf, Kh_surf);
              }
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

        if (only_Ez)
            var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, square(this->shift), var_eig.SORTED_MODULUS);
        else
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
        if (only_Ez)
          {
            for (int i = 0; i < lambda.GetM(); i++) {
                betas.Get(i) = sqrt(lambda.Get(i));
            }
          }
        else
          {
            Copy(lambda, betas);
          }

        DISP(type_coord_other_fields);

        VectComplex_wp x;
        for (int i = 0; i < betas.GetM(); i++)
          {
            GetCol(eigen_vec, i, x);
            VectComplex_wp Ez;
            VectComplex_wp Hz; // on ne prend que les données de Ez et de Hz
            VectComplex_wp X = x; X.Resize(nb_dof_Ez + nb_dof_Hz);
            VectComplex_wp Y = x; Y.Resize(nb_dof_Ez + nb_dof_Hz);
            Ez.Reallocate(nb_dof_Hz); Ez.Zero();
            Hz.Reallocate(nb_dof_Hz); Hz.Zero();

    #ifdef SELDON_WITH_SLEPC
            // vérification (norme de T(L) X)
            ModeEzNonLinearDirichletEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
                                                this->DofKeptNonDir, this->IndexNonDirichlet);

            var_eig.MltOperator(betas(i), SeldonNoTrans, X, Y);
            cout << "VP " << i << ": " << betas(i) << " ; norme : " << Norm2(Y) / Norm2(X) << endl;

            x.Reallocate(2*nb_dof_Hz); x.Zero();
            if (type_CLA == DIRICHLET)
              var_eig.ExtractDirichlet(X, x);
            else
              x = X;
            for (int i = 0; i < nb_dof_Hz; i++)
              {
                Ez(i) = x(i);
                Hz(i) = x(nb_dof_Hz + i);
              }
            ComputeIntensity(Ez, vec_Vk, i);

            WriteOutputFile(Ez, Hz, i, betas(i));
    #endif
          }
      }

    return;
  }

};

#ifdef SELDON_WITH_SLEPC
// on inclut les fonctions definies dans le cxx
#include "ResolutionNonLinear.cxx"
#include "ModeEzNonLinearDirichletEigenProblem.cxx"
#include "ModeEzNonLinearBGTEigenProblem.cxx"
#endif
#include "IOSolution.cxx"
#include "Mesh.cxx"
#include "AssemblageMatricesDirichlet.cxx"
#include "AssemblageMatrices.cxx"
#include "FiniteElement.cxx"

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc < 2)
    {
      cout << "Data file not given" << endl;
      cout << "Usage : ./mode_Ez.x data_file.ini" << endl;
      abort();
    }

  // on appelle RunAll avec le fichier de donnees
  ModeEz_Solver var;
  string input_file(argv[1]);
  var.RunAll(PATH_AUGUSTIN + to_str("/Datas/") + input_file);

  return FinalizeMontjoie();
}
