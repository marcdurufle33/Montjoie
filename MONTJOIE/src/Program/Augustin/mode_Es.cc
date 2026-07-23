#include "Montjoie.hxx"
#define PATH_AUGUSTIN "./src/Program/Augustin"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;

class ModeEs_Solver;

#ifdef SELDON_WITH_SLEPC
// classe pour resoudre un pb aux valeurs propres non-lineaire avec Slepc
// forme splitte : T(L) = \sum_k f_k(L) A_k
// et f_k est une fraction rationnelle
class ModeEsNonLinearEigenProblem : public NonLinearEigenProblem_Base<Complex_wp>
{
  ModeEs_Solver& var; int nb_rows, nb_dof_Es, nb_dof_Hs, nb_dof_all;
  Vector<Matrix<Complexe, Prop, Storage> >& vec_Mh, &vec_Kh, &vec_Ch, &vec_Dh, &vec_Th, &vec_Uh;
  Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0, vec_Inv_D_tauBeta, vec_ChgmtVariables, vec_ktrSquared;
  SparseDistributedSolver<Complex_wp> mat_lu;
  Vector<int>& DofKeptDir; Vector<int>& IndexDirichlet;

  // pour la formulation splittee
  Matrix<Complexe, Prop, Storage> A1_sum, A1_sumH; // pour la f_i constante (= 1)
  Vector<int> split_reference;
  bool twisted; // vecteur sur les ref ensuite
  Vector<bool> vec_twisted;
  Vector<bool> index_used;

private:
  void ComputeInvDBetaMatrix(const Complex_wp& L);
  void ComputeChgmtVariable(const Complex_wp& L);
  void Mult(const Matrix<Complexe, Prop, Storage>& A, const  Matrix<Complexe, Symmetric, DiagonalRow>& D, Matrix<Complexe, Prop, Storage>& M);
  int GetMatrixIndex(int num);

public:
  ModeEsNonLinearEigenProblem(int N, ModeEs_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
                              Vector<Matrix<Complexe, Prop, Storage> >& Dh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Th,
                              Vector<Matrix<Complexe, Prop, Storage> >& Uh,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& DiagTau0,
                              Vector<Matrix<Complexe, Symmetric, DiagonalRow> >& ktr2,
                              Vector<int>& DofDir, Vector<int>& IndexDir,
                              bool use_split = false);

  void ExtractDirichlet(const VectComplex_wp& Xdir, VectComplex_wp& X);
  void ExpandDirichlet(const VectComplex_wp& X, VectComplex_wp& Xdir);
  void RemoveDirichletRowCol(const Matrix<Complexe, Prop, Storage>&,
                             Matrix<Complexe, Prop, Storage>&);

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

class ModeEsPolynomialProblem : public PolynomialSparseEigenProblem<Complex_wp, DistributedMatrix<Complex_wp, General, ArrayRowSparse>, DistributedMatrix<Complex_wp, General, ArrayRowSparse>>
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
  ModeEsPolynomialProblem(int N, ModeEs_Solver& var0,
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

class ModeEs_Solver : public InputDataProblem_Base,
                      public VarFiniteElement<Dimension2>
{
#ifdef SELDON_WITH_SLEPC
  friend class ModeEsNonLinearEigenProblem;
  friend class ModeEsPolynomialProblem;
#endif
protected:
  // pulsation omega
  Real_wp omega, torsion;

  // bords Dirichlet
  Vector<int> ref_dirichlet;

  // pour les bords quasi-periodiques
  // TinyVector<int, 2> ref_periodic;
  // int num_mode_periodic;
  // Vector<Vector<bool> > is_dof_with_phase;

  // maillage
  Mesh<Dimension2> mesh; string name_mesh;

  // numerotation
  MeshNumbering<Dimension2> mesh_num;
  Vector<int> IndexDirichlet, DofKeptDir;

  // indices physiques
  Vector<Complex_wp> ref_epsilon, ref_mu, ref_sigma;
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

  // type de coordonnées pour les champs transversaux
  int type_coord_other_fields;
  enum {CARTESIAN, POLAR, NONE};

  // shift pour le calcul des valeurs propres
  Complex_wp shift; int nb_asked_eigenval;

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
  }

  void InitIndices(int ref_max)
  {
    ref_epsilon.Reallocate(ref_max+1);
    ref_mu.Reallocate(ref_max+1);
    ref_sigma.Reallocate(ref_max+1);
    ref_tau.Reallocate(ref_max+1);
    ref_tau.Zero();
  }

  void SetInputData(const string& keyword, const VectString& parameters)
  {
    // on appelle SetInputData pour mesh et mesh_num
    mesh.SetInputData(keyword, parameters);
    mesh_num.SetInputData(keyword, parameters);

    // autres mots-cles
    if (keyword == "FileMesh")
      name_mesh = PATH_AUGUSTIN + to_str("/Datas/") + parameters(0);
    else if (keyword == "Frequency")
      this->omega = 2.0*pi_wp*to_num<Real_wp>(parameters(0)) + to_num<Real_wp>(parameters(1));
    else if (keyword == "Periode_Torsade")
      {
        if (parameters(0) == "Inf")
          this->torsion = 0.0;
        else
          this->torsion = 2.0*pi_wp/to_num<Real_wp>(parameters(0));
      }
    else if (keyword == "MateriauDielec")
      {
        // on garde le meme format que pour mode_maxwell (parameters(1) pas utilise)
        int ref = to_num<int>(parameters(0));
        ref_epsilon(ref) = to_num<Complex_wp>(parameters(2));
        ref_mu(ref) = to_num<Complex_wp>(parameters(3));
        ref_sigma(ref) = to_num<Complex_wp>(parameters(4));
        if (parameters(5) == "Inf")
          ref_tau(ref) = 0.0;
        else
          ref_tau(ref) = 2.0*pi_wp/to_num<Real_wp>(parameters(5));
      }
    else if (keyword == "ReferenceDirichlet")
      {
        // liste des bords de type Dirichlet
        ref_dirichlet.Reallocate(parameters.GetM());
        for (int i = 0; i < parameters.GetM(); i++)
          ref_dirichlet(i) = to_num<int>(parameters(i));
      }
    else if (keyword == "SismoPlane")
      {
      	// output on a plane
      	int i = var_grid.GetM();
      	var_grid.Resize(var_grid.GetM()+1);
      	var_grid(i).SetInputData(keyword, parameters);
      }
    else if ( keyword == "SismoPointsFile")
      {
      	// output on points given in a file
      	int i = var_grid.GetM();
      	var_grid.Resize(var_grid.GetM()+1);
      	var_grid(i).SetInputData(keyword, parameters);
      }
    else if (keyword == "SismoLine")
      {
      	// output on a plane
      	int i = var_grid.GetM();
      	var_grid.Resize(var_grid.GetM()+1);
      	var_grid(i).SetInputData(keyword, parameters);
      }
    else if ( keyword == "SismoCircle")
      {
      	// output on a circle
      	int i = var_grid.GetM();
      	var_grid.Resize(var_grid.GetM()+1);
      	var_grid(i).SetInputData(keyword, parameters);
      }
    else if (keyword == "Eigensolver")
      {
        cout << "Choix du solveur : " << parameters(0) << endl;
        if (parameters(0) == "Linear" || parameters(0) == "Linear_Arpack")
          type_eigensolver = LINEAR_ARPACK;
        else if (parameters(0) == "Linear_Slepc")
          type_eigensolver = LINEAR_SLEPC;
        else if (parameters(0) == "Rii")
          type_eigensolver = NL_RII;
        else if (parameters(0) == "Slp")
          type_eigensolver = NL_SLP;
        else if (parameters(0) == "Nleigs")
          type_eigensolver = NL_EIGS;
        else if (parameters(0) == "Polynomial")
          type_eigensolver = POLYNOMIAL;
        else
          {
            cout << "Unknown eigensolver " << parameters(0) << endl;
            abort();
          }

        for (int k = 1; k < parameters.GetM(); k++)
          {
            if (parameters(k) == "Split")
              use_split_formulation = true;
            if (parameters(k) == "Explicit")
              use_explicit_matrix = true;
            if (parameters(k) == "Default")
              use_default_petsc_solver = true;
          }
      }
    else if (keyword == "ShiftEigenvalue")
      {
        shift = to_num<Complex_wp>(parameters(0));
      }
    else if (keyword == "NumberEigenvalue")
      {
        nb_asked_eigenval = to_num<int>(parameters(0));
      }
    else if (keyword == "CoordOtherFields")
      {
        if (parameters(0) == "Cartesian")
          type_coord_other_fields = CARTESIAN;
        else if (parameters(0) == "Polar")
          type_coord_other_fields = POLAR;
        else if (parameters(0) == "None")
          type_coord_other_fields = NONE;
        else
          {
            cout << "Unknown type of cooordonates " << parameters(0) << ". The other fields won't be computed." << endl;
          }
      }
    else if (keyword == "ReferenceConductor")
      {
        int N = parameters.GetM()/2;
        ref_conductor.Reallocate(N);
        ref_vol_conductor.Reallocate(N);
        phase_conductor.Reallocate(N);
        Complex_wp alpha = exp(Iwp*2.0*pi_wp/N);
        Complex_wp phase(1, 0);
        for (int i = 0; i < N; i++)
          {
            ref_conductor(i) = to_num<int>(parameters(2*i));
            ref_vol_conductor(i) = to_num<int>(parameters(2*i+1));
            phase_conductor(i) = phase;
            phase *= alpha;
          }
      }
  }

  // version simplifiee de ConstructFiniteElement (dans VarProblem.cxx)
  void ConstructFiniteElement(const string& name_elt)
  {
    // previous finite elements are cleared if present
    this->ClearFiniteElement();

    // orders present in the mesh
    TinyVector<IVect, 4> order;
    this->mesh_num.GetOrder(order);

    Vector<bool> change_elt(this->mesh.GetNbElt());
    change_elt.Fill(true);

    this->AddFiniteElement(name_elt, order, change_elt, this->mesh, this->mesh_num,
                           0, false);

    this->mesh_num.GetOrderQuadrature(order);
    RemoveDuplicate(order(0));

    int type_integration_edge, type_integration_tri, type_integration_quad;
    Vector<MeshNumbering<Dimension2>* > other_mesh_num;
    this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, other_mesh_num, type_integration_edge,
				     type_integration_tri, type_integration_quad);

    this->mesh_num.number_map.ConstructQuadrature2D(order, type_integration_edge);
  }

  // traitement de Dirichlet (et quasi-periodique)
  void RenumberDofs()
  {
    Vector<bool> is_dirichlet(this->mesh.GetNbReferences()+1);
    is_dirichlet.Fill(false);
    for (int i = 0; i < this->ref_dirichlet.GetM(); i++)
      is_dirichlet(ref_dirichlet(i)) = true;

    // on met le ddl 0 en Dirichlet (dans le cas quasi-periodique avec m <> 0)
    Vector<int> DofDirichlet(this->mesh_num.GetNbDof());
    DofDirichlet.Fill(-1);

    int nb_dof_dirichlet = 0;

    cout << "Number of Dirichlet dofs = " << nb_dof_dirichlet << endl;

    // on renumerote
    if (nb_dof_dirichlet > 0)
      {
        Vector<int> IndexDof(this->mesh_num.GetNbDof());
        IndexDof.Fill(-1);
        int nodl = 0;
        for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
          {
            if (DofDirichlet(i) == -1)
              IndexDof(i) = nodl++;
          }

        //DISP(IndexDof);

        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          {
            //DISP(i); DISP(this->mesh_num.Element(i).GetNodle());
            for (int j = 0; j < this->mesh_num.Element(i).GetNbDof(); j++)
              {
                int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
                if (IndexDof(num_dof) <= -2)
                  {
                    this->mesh_num.Element(i).SetNumberDof(j, -IndexDof(num_dof)-2);
                  }
                else
                  this->mesh_num.Element(i).SetNumberDof(j, IndexDof(num_dof));
              }

            // DISP(this->is_dof_with_phase(i));
            //DISP(this->mesh_num.Element(i).GetNodle());
          }

        this->mesh_num.SetNbDof(nodl);
      }

    // boucle sur les aretes de bord pour lister les ddls Dirichlets
    nb_dof_dirichlet = 0;
    DofDirichlet.Fill(-1);
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
        int ref = this->mesh.BoundaryRef(i).GetReference();
        if (is_dirichlet(ref))
          {
            int num_elem = this->mesh.BoundaryRef(i).numElement(0);
            int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
            const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
            const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
            int nb_dof = Fb.GetNbDofBoundary(num_loc);
            for (int j = 0; j < nb_dof; j++)
              {
                int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
                int num_dof = Nodle(num_dof_loc);
                if ((num_dof >= 0) && (DofDirichlet(num_dof) == -1))
                  {
                    DofDirichlet(num_dof) = 0;
                    nb_dof_dirichlet++;
                  }
              }
          }
      }

    IndexDirichlet.Reallocate(this->mesh_num.GetNbDof());
    IndexDirichlet.Fill(-1);
    DofKeptDir.Reallocate(this->mesh_num.GetNbDof()-nb_dof_dirichlet);
    nb_dof_dirichlet = 0;
    for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
      if (DofDirichlet(i) == -1)
        {
          IndexDirichlet(i) = nb_dof_dirichlet;
          DofKeptDir(nb_dof_dirichlet) = i;
          nb_dof_dirichlet++;
        }
  }

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
                        )
  {
    Matrix<Complexe, Symmetric, DiagonalRow> diag_D_tau0;
    // on alloue les matrices
    int N = this->mesh_num.GetNbDof();
    vec_Mh.Reallocate(this->ref_epsilon.GetM());
    vec_Kh.Reallocate(this->ref_epsilon.GetM());
    vec_Ch.Reallocate(this->ref_epsilon.GetM());
    vec_Dh.Reallocate(this->ref_epsilon.GetM());
    vec_Th.Reallocate(this->ref_epsilon.GetM());
    vec_Uh.Reallocate(this->ref_epsilon.GetM());
    vec_diag_D_tau0.Reallocate(this->ref_epsilon.GetM()); // ce sera une ref sur les tau
    vec_diag_D_tau0_Bis.Reallocate(this->ref_epsilon.GetM());
    Lk_conductor.Reallocate(this->ref_epsilon.GetM());
    Lk_conductor.Zero();
    vec_Lk.Reallocate(this->ref_epsilon.GetM());
    vec_Lk.Zero();
    vec_Vk.Reallocate(N);
    vec_Vk.Zero();
    vec_Vk2.Reallocate(this->ref_epsilon.GetM());
    vec_Pos.Reallocate(N); vec_Pos.Fill(-21);
    diag_D_tau0.Reallocate(N,N);

    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int ref = 0; ref < index_used.GetM(); ref++)
      if (index_used(ref))
        {
          vec_Mh(ref).Reallocate(N, N);
          vec_Kh(ref).Reallocate(N, N);
          vec_Ch(ref).Reallocate(N, N);
          vec_Dh(ref).Reallocate(N, N);
          vec_Th(ref).Reallocate(N, N);
          vec_Uh(ref).Reallocate(N, N);
          vec_Vk2(ref).Reallocate(N);
          vec_diag_D_tau0(ref).Reallocate(N, N);
          vec_diag_D_tau0_Bis(ref).Reallocate(N, N);
        }

    VectR2 s;
    SetPoints<Dimension2> PointsElem;
    SetMatrices<Dimension2> MatricesElem;


    Real_wp threshold = 1e-15;

    // boucle sur les elements
    VectReal_wp Bzero, Bmass; Matrix2_2 dfjm1, Ctmp, Cref, Dtmp, Dref;
    VectR2 Dzero, Tadv, Uadv, Ezero; Vector<Matrix2_2> Astiff, Cstiff, Dstiff, Czero;
    TinyVector<bool, 4> null_stiff(true, false, true, true);
    TinyVector<bool, 4> null_adv(true, true, true, false);
    Matrix<Real_wp> mat_elem_Mh, mat_elem_Kh, mat_elem_Ch, mat_elem_Dh, mat_elem_Th, mat_elem_Uh;
    Cref(0, 1) = -1.0; Cref(1, 0) = 1.0; // pour "retourner" les dérivées
    DISP(this->mesh.GetNbElt());
    num_conduc_ddl.Reallocate(N); num_conduc_ddl.Fill(-1);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);
        /*if (!Fb.MassLumping())
          {
            cout << "Valid strategy for mass-lumped elements" << endl;
            abort();
            }*/

        int ref = this->mesh.Element(i).GetReference();
        // on recupere le numero du conducteur si present
        int num_conductor = -1;
        for (int k = 0; k < ref_vol_conductor.GetM(); k++)
          if (ref_vol_conductor(k) == ref)
            num_conductor = k;

        // calcul des matrices jacobiennes sur les points de quadrature
        this->mesh.GetVerticesElement(i, s);
        Fb.FjElem(s, PointsElem, this->mesh, i);
        Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

        int Nquad = Fb.GetNbPointsQuadratureInside();
        Bmass.Reallocate(Nquad); Astiff.Reallocate(Nquad);
        Cstiff.Reallocate(Nquad); Dstiff.Reallocate(Nquad);
        Uadv.Reallocate(Nquad); Tadv.Reallocate(Nquad);

        // on recupere le numero de conducteur pour chaque ddl
        VectReal_wp phi;
        int nb_dof = Fb.GetNbDof();
        IVect num_ddl = this->mesh_num.Element(i).GetNodle();
        if (num_conductor >= 0)
          for (int k = 0; k < nb_dof; k++)
            if (num_ddl(k) >= 0)
              num_conduc_ddl(num_ddl(k)) = num_conductor;

        // boucle sur les points de quadrature
        for (int j = 0; j < Nquad; j++)
          {
            // pour le calcul de Mh (matrice de masse)
            Real_wp weight = Fb.WeightsND(j);
            Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
            Real_wp x_quadGlob = PointsElem.GetPointQuadrature(j)(0);
            Real_wp y_quadGlob = PointsElem.GetPointQuadrature(j)(1);
            R2 TRef(x_quadGlob, y_quadGlob);
            R2 URef(y_quadGlob, - x_quadGlob);
            // R2 TRef(0.0, 0.0);
            // R2 URef(0.0, 0.0);
            Bmass(j) = jacob*weight;

            // calcul de Vk
            // pour avoir phi_i(xi_j) appeler Fb.GetValuePhiOnQuadraturePoint(j, phi);

            vec_Lk(ref) += jacob*weight;
            Fb.GetValuePhiOnQuadraturePoint(j, phi);
            for (int k = 0; k < nb_dof; k++)
              if (num_ddl(k) >= 0)
                vec_Vk2(ref)(num_ddl(k)) += jacob*weight*phi(k);

            if (num_conductor >= 0)
              {
                Lk_conductor(num_conductor) += jacob*weight;
                Fb.GetValuePhiOnQuadraturePoint(j, phi);
                for (int k = 0; k < nb_dof; k++)
                  if (num_ddl(k) >= 0)
                    vec_Vk(num_ddl(k)) += jacob*weight*phi(k);
              }

            // pour le calcul de Kh (matrice de rigidite)
            GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);
            MltTrans(dfjm1, dfjm1, Astiff(j));
            Mlt(jacob*weight, Astiff(j));

            // pour le calcul de Ch (matrice de couplage)
            MltTrans(Cref, dfjm1, Ctmp);
            Mlt(dfjm1, Ctmp, Cstiff(j));
            Mlt(jacob*weight, Cstiff(j));

            // pour le calcul de Dh
            // calcul de l'application sur les gradients
            Dref(0, 0) = y_quadGlob * y_quadGlob;
            Dref(1, 0) = x_quadGlob * y_quadGlob;
            Dref(0, 1) = -x_quadGlob * y_quadGlob;
            Dref(1, 1) = -x_quadGlob * x_quadGlob;

            // Dref(0, 0) = 0.0;
            // Dref(1, 0) = 0.0;
            // Dref(0, 1) = 0.0;
            // Dref(1, 1) = 0.0;


            MltTrans(Dref, dfjm1, Dtmp);
            Mlt(dfjm1, Dtmp, Dstiff(j));
            Mlt(jacob*weight, Dstiff(j));

            // Vadv(0) = y_quadGlob;
            // (grad(phi_j) vectoriel xi) phi_i
            Mlt(dfjm1, URef, Uadv(j));
            Mlt(jacob*weight, Uadv(j));

            // (grad(phi_j) scalaire xi) phi_i
            Mlt(dfjm1, TRef, Tadv(j));
            Mlt(jacob*weight, Tadv(j));
          }

        if (Fb.GetNbPointsDof() != Fb.GetNbDof())
          {
            cout << "Element fini pas nodal " << endl;
            abort();
          }

        // boucle sur les ddls (supposes nodaux)
        for (int k = 0; k < Fb.GetNbPointsDof(); k++)
          {
            Real_wp x_dofGlob = PointsElem.GetPointDof(k)(0);
            Real_wp y_dofGlob = PointsElem.GetPointDof(k)(1);

            if (num_ddl(k) >= 0)
              {
                vec_Pos(num_ddl(k))(0) = x_dofGlob;
                vec_Pos(num_ddl(k))(1) = y_dofGlob;

                Real_wp r2 = x_dofGlob * x_dofGlob + y_dofGlob * y_dofGlob;
                diag_D_tau0.Get(num_ddl(k), num_ddl(k)) = 1 + this->ref_tau(ref) * this->ref_tau(ref) * r2;
                vec_diag_D_tau0_Bis(ref).Get(num_ddl(k), num_ddl(k)) = 1 + this->ref_tau(ref) * this->ref_tau(ref) * r2;
              }
          }

        // calcul des matrices de masse et rigidite elementaires
        mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
        mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();
        mat_elem_Ch.Reallocate(nb_dof, nb_dof); mat_elem_Ch.Zero();
        mat_elem_Dh.Reallocate(nb_dof, nb_dof); mat_elem_Dh.Zero();
        mat_elem_Th.Reallocate(nb_dof, nb_dof); mat_elem_Th.Zero();
        mat_elem_Uh.Reallocate(nb_dof, nb_dof); mat_elem_Uh.Zero();

        Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Astiff, Dzero, Ezero, null_stiff, mat_elem_Kh);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Cstiff, Dzero, Ezero, null_stiff, mat_elem_Ch);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Dstiff, Dzero, Ezero, null_stiff, mat_elem_Dh);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Dzero, Tadv, null_adv, mat_elem_Th);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Czero, Dzero, Uadv, null_adv, mat_elem_Uh);

        // assemblage des matrices
        IVect permut(nb_dof); permut.Fill();
        Sort(num_ddl, permut);

        Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
        col.Fill(-1); val.Zero();
        for (int j = 0; j < nb_dof; j++)
          if (num_ddl(j) >= 0)
            {
              nb_val = 0;

              //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                  {
                    col(nb_val) = num_ddl(k);
                    val(nb_val) = mat_elem_Mh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              nb_val = 0;
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Kh(permut(j), permut(k))) >= threshold))
                  {

                    col(nb_val) = num_ddl(k);
                    val(nb_val) = mat_elem_Kh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Kh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              nb_val = 0;
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Ch(permut(j), permut(k))) >= threshold))
                  {
                    col(nb_val) = num_ddl(k);
                    val(nb_val) = mat_elem_Ch(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Ch(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              nb_val = 0;
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Dh(permut(j), permut(k))) >= threshold))
                  {
                    col(nb_val) = num_ddl(k);
                    val(nb_val) = mat_elem_Dh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Dh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              nb_val = 0;
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Th(permut(j), permut(k))) >= threshold))
                  {
                    col(nb_val) = num_ddl(k);
                    val(nb_val) = mat_elem_Th(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Th(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              nb_val = 0;
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Uh(permut(j), permut(k))) >= threshold))
                  {
                    col(nb_val) = num_ddl(k);
                    val(nb_val) = mat_elem_Uh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Uh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
            }
      }


    for (int ref = 0; ref < index_used.GetM(); ref++)
      if (index_used(ref))
        {
          vec_Ch(ref).RemoveSmallEntry(threshold);
          for (int i = 0; i < N; i++)
            {
              Real_wp Lk = vec_Lk(ref);
              vec_Vk2(ref)(i) /= Lk;
            }
        }

    DISP(Lk_conductor);
    DISP(vec_Lk);
    // on divise les Vk par l'aire pour avoir V_i = 1/|L_k| \int_{L_k} \varphi_i dx
    for (int i = 0; i < N; i++)
      if (num_conduc_ddl(i) >= 0)
        {
          Real_wp Lk = Lk_conductor(num_conduc_ddl(i));
          vec_Vk(i) /= Lk;
        }
    vec_diag_D_tau0(0) = diag_D_tau0;
  }

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
                                              Vector<int>& nbRefNoeud)
  {
    cout << "Construction des matrices de matériaux torsadés" << endl;
    // à refaire intégralement lors de la prise en compte des torsades multiples
    vec_D_tau0_Bis.Reallocate(this->ref_epsilon.GetM());
    vec_ktrSquared.Reallocate(this->ref_epsilon.GetM());
    int N = this->mesh_num.GetNbDof();
    nbRefNoeud.Reallocate(N);
    nbRefNoeud.Zero();

    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int ref = 0; ref < index_used.GetM(); ref++)
      if (index_used(ref))
        {
          DISP(ref);
          vec_D_tau0_Bis(ref).Reallocate(N,N);
          vec_ktrSquared(ref).Reallocate(N,N);
          Complex_wp k2 = Iwp * this->omega * this->ref_mu(ref) * (-Iwp * this->omega * this->ref_epsilon(ref) + this->ref_sigma(ref));
          for (int i = 0 ; i < N ; i++)
            {
              Real_wp val = vec_diag_tau0_Bis(ref)(i,i);
              if (val > 0.9) // tous les cas intéressants sont >= 1, on met 0.9 pour ne pas conter les éventuels 0 non comptabilisés numériquement
                {
                  vec_D_tau0_Bis(ref).Get(i,i) = -k2 * val;
                  vec_ktrSquared(ref).Get(i,i) = k2 * (val - 1.0);
                  nbRefNoeud.Get(i) += 1;
                }
              else
                {
                  vec_D_tau0_Bis(ref).Get(i,i) = 0.0;
                  vec_ktrSquared(ref).Get(i,i) = 0.0;
                }
            }
        }
  }


  void ComputeMatrixV2ToV1(const Complex_wp& L,
                            const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                            const Vector<int>& nbRefNoeud,
                            Matrix<Complexe, Symmetric, DiagonalRow>& Dtau)
  {
    int N = this->mesh_num.GetNbDof();
    Complex_wp beta = L * this->omega;
    Dtau.Reallocate(N,N);

    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int i = 0 ; i < N ; i++)
      {
        Dtau.Get(i,i) = 0.0;
        for (int ref = 0; ref < index_used.GetM(); ref++)
          if (index_used(ref))
            {
              Matrix<Complexe, Symmetric, DiagonalRow> D = vec_D_tau0_Bis(ref);
              if (D(i,i) != 0.0)
                Dtau.Get(i,i) += (D.Get(i,i) + beta*beta) / nbRefNoeud.Get(i);
                // Dtau.Get(i,i) += (D.Get(i,i) + beta*beta);
              // DISP(i+1);
              // cout << D(i,i) << ", beta : " << beta << endl;
              // DISP(Dtau2.Get(i,i));
              // DISP(nbRefNoeud.Get(i));
            }
      }
  }

  void ComputeMatrixV2ToEs(const Complex_wp& L,
                            const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                            const Vector<int>& nbRefNoeud,
                            Matrix<Complexe, Symmetric, DiagonalRow>& Dtau2)
  {
    int N = this->mesh_num.GetNbDof();
    Complex_wp beta = L * this->omega;
    Dtau2.Reallocate(N,N);

    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int i = 0 ; i < N ; i++)
      {
        Dtau2.Get(i,i) = 0.0;
        for (int ref = 0; ref < index_used.GetM(); ref++)
          if (index_used(ref))
            {
              Matrix<Complexe, Symmetric, DiagonalRow> D = vec_D_tau0_Bis(ref);
              if (D(i,i) != 0.0)
                Dtau2.Get(i,i) += (D.Get(i,i) + beta*beta)*(D.Get(i,i) + beta*beta);
                // Dtau2.Get(i,i) += (D.Get(i,i) + beta*beta)*(D.Get(i,i) + beta*beta) / nbRefNoeud.Get(i);
              // DISP(i+1);
              // cout << D(i,i) << ", beta : " << beta << endl;
              // DISP(Dtau2.Get(i,i));
              // DISP(nbRefNoeud.Get(i));
            }
      }
  }


  void ComputeEsHs(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                  const Vector<int>& nbRefNoeud,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Uh,
                  const Vector<Matrix<Complexe, Prop, Storage>> &vec_Th,
                  const Vector<Matrix<Complexe, Symmetric, DiagonalRow>> &vec_ktrSquared,
                  const Complex_wp& L,
                  const Vector<Complex_wp>& V2, const Vector<Complex_wp>& W2,
                  Vector<Complex_wp>& Es, Vector<Complex_wp>& Hs)
  {
    int N = this->mesh_num.GetNbDof();
    Vector<Complex_wp> V1(N);
    Vector<Complex_wp> W1(N);
    Matrix<Complexe, Symmetric, DiagonalRow> DTau;
    Matrix<Complexe, General, ArrayRowSparse> MM;
    Complex_wp beta = L * this->omega;

    ComputeMatrixV2ToV1(L, vec_D_tau0_Bis, nbRefNoeud, DTau);

    Es.Reallocate(N); Es.Zero();
    Hs.Reallocate(N); Hs.Zero();
    MM.Reallocate(N,N);

    Mlt(DTau, V2, V1);
    Mlt(DTau, W2, W1);

    for (int ref = 0; ref < this->ref_epsilon.GetM(); ref++)
      {
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
            Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
            Complex_wp k2 = -eps_sigma*mu_tilde;
            Complex_wp one(1,0);
            Vector<Complex_wp> V2t(N);
            Vector<Complex_wp> W2t(N);
            Mlt(vec_ktrSquared(ref), V2, V2t);
            Mlt(vec_ktrSquared(ref), W2, W2t);

            MltAdd(- k2 + beta * beta, vec_Mh(ref), V1, one, Es);
            MltAdd(- k2 + beta * beta, vec_Mh(ref), W1, one, Hs);
            MltAdd(- Iwp * beta * this->ref_tau(ref), vec_Uh(ref), V1, one, Es);
            MltAdd(- Iwp * beta * this->ref_tau(ref), vec_Uh(ref), W1, one, Hs);
            MltAdd(mu_tilde * this->ref_tau(ref), vec_Th(ref), W1, one, Es);
            MltAdd(- eps_sigma * this->ref_tau(ref), vec_Th(ref), V1, one, Hs);
            MltAdd(- 2 * mu_tilde * this->ref_tau(ref), vec_Mh(ref), W2t, one, Es);
            MltAdd(2 * eps_sigma * this->ref_tau(ref), vec_Mh(ref), V2t, one, Hs);

            Add(1.0, vec_Mh(ref), MM);
          }
      }
      SparseDistributedSolver<Complex_wp> mat_lu;
      mat_lu.Factorize(MM);
      mat_lu.Solve(Es);
      mat_lu.Solve(Hs);
  }

  void ComputeMatrixF2ToF1(const Complex_wp& L,
                            const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                            const Vector<int>& nbRefNoeud,
                            Matrix<Complexe, Symmetric, DiagonalRow>& Dtau)
  {
    int N = this->mesh_num.GetNbDof();
    Complex_wp beta = L * this->omega;
    Dtau.Reallocate(N,N);

    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int i = 0 ; i < N ; i++)
      {
        Dtau.Get(i,i) = 0.0;
        for (int ref = 0; ref < index_used.GetM(); ref++)
          if (index_used(ref))
            {
              Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
              Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
              Complex_wp k2 = -eps_sigma*mu_tilde;
              Matrix<Complexe, Symmetric, DiagonalRow> D = vec_D_tau0_Bis(ref);
              if (D(i,i) != 0.0)
                Dtau.Get(i,i) += ((D.Get(i,i) + beta*beta) / (-k2 + beta*beta)) / nbRefNoeud.Get(i);
                // Dtau.Get(i,i) += (D.Get(i,i) + beta*beta);
              // DISP(i+1);
              // cout << D(i,i) << ", beta : " << beta << endl;
              // DISP(Dtau2.Get(i,i));
              // DISP(nbRefNoeud.Get(i));
            }
      }
  }

  void FromF2G2ToEsHs(const Vector<Matrix<Complexe, Symmetric, DiagonalRow>>& vec_D_tau0_Bis,
                  const Vector<int>& nbRefNoeud,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Mh,
                  const Vector<Matrix<Complexe, Prop, Storage>>& vec_Uh,
                  const Vector<Matrix<Complexe, Prop, Storage>> &vec_Th,
                  const Vector<Matrix<Complexe, Symmetric, DiagonalRow>> &vec_ktrSquared,
                  const Complex_wp& L,
                  const Vector<Complex_wp>& F2, const Vector<Complex_wp>& G2,
                  Vector<Complex_wp>& Es, Vector<Complex_wp>& Hs)
  {
    int N = this->mesh_num.GetNbDof();
    Matrix<Complexe, General, ArrayRowSparse> MM;
    Complex_wp beta = L * this->omega;

    // ComputeMatrixF2ToF1(L, vec_D_tau0_Bis, nbRefNoeud, DTau);

    Es.Reallocate(N); Es.Zero();
    Hs.Reallocate(N); Hs.Zero();
    MM.Reallocate(N,N);


    for (int ref = 0; ref < this->ref_epsilon.GetM(); ref++)
      {
        if (vec_Mh(ref).GetM() > 0)
          {
            Complex_wp eps_sigma = -Iwp*this->omega*this->ref_epsilon(ref) + this->ref_sigma(ref);
            Complex_wp mu_tilde = -Iwp*this->omega*this->ref_mu(ref);
            Complex_wp k2 = -eps_sigma*mu_tilde;
            Complex_wp DBeta = - k2 + beta * beta;
            Complex_wp one(1,0);
            Vector<Complex_wp> F1(F2);
            Vector<Complex_wp> G1(G2);
            Vector<Complex_wp> temp(N); temp.Zero();
            Vector<Complex_wp> F2t(N);
            Vector<Complex_wp> G2t(N);

            Mlt(vec_ktrSquared(ref), F2, F2t);
            Mlt(vec_ktrSquared(ref), G2, G2t);
            Add(-1.0 / DBeta, F2t, F1);
            Add(-1.0 / DBeta, G2t, G1);
            //
            // for (int i = 0 ; i < N ; i++)
            //   {
            //     if (vec_ktrSquared(ref)(i,i) == 0)
            //      {
            //        F1.Get(i) = 0.0;
            //      }
            //   }

            WriteOutputFileScalar(F1, "F1_"+ to_str(ref));
            WriteOutputFileScalar(G1, "G1_"+ to_str(ref));


            MltAdd(one, vec_Mh(ref), F1, one, Es);
            MltAdd(one, vec_Mh(ref), G1, one, Hs);
            MltAdd(- Iwp * beta * this->ref_tau(ref) / DBeta, vec_Uh(ref), F1, one, Es);
            MltAdd(- Iwp * beta * this->ref_tau(ref) / DBeta, vec_Uh(ref), G1, one, Hs);
            MltAdd(mu_tilde * this->ref_tau(ref) / DBeta, vec_Th(ref), G1, one, Es);
            MltAdd(- eps_sigma * this->ref_tau(ref) / DBeta, vec_Th(ref), F1, one, Hs);
            MltAdd(- 2 * mu_tilde * this->ref_tau(ref) / (DBeta * DBeta), vec_Mh(ref), G2t, one, Es);
            MltAdd(2 * eps_sigma * this->ref_tau(ref) / (DBeta * DBeta), vec_Mh(ref), F2t, one, Hs);

            Add(1.0, vec_Mh(ref), MM);
          }
      }
      SparseDistributedSolver<Complex_wp> mat_lu;
      mat_lu.Factorize(MM);
      mat_lu.Solve(Es);
      mat_lu.Solve(Hs);
  }


  // version simplifiee de InitGrid pour localiser les points du maillage sur les sorties
  void InitGrid()
  {
    // boundaries of the grid are set to boundaries of the domain if required
    for (int i = 0; i < var_grid.GetM(); i++)
      {
	var_grid(i).SetXmin0(this->mesh.GetXmin());
	var_grid(i).SetXmax0(this->mesh.GetXmax());
	var_grid(i).SetYmin0(this->mesh.GetYmin());
	var_grid(i).SetYmax0(this->mesh.GetYmax());
	var_grid(i).SetZmin0(this->mesh.GetZmin());
	var_grid(i).SetZmax0(this->mesh.GetZmax());
      }

    all_points_display.SetXmin(this->mesh.GetXmin());
    all_points_display.SetXmax(this->mesh.GetXmax());
    all_points_display.SetYmin(this->mesh.GetYmin());
    all_points_display.SetYmax(this->mesh.GetYmax());
    all_points_display.SetZmin(this->mesh.GetZmin());
    all_points_display.SetZmax(this->mesh.GetZmax());

    // predefined grids are generated
    // all_points_display is an union of all the predefined grids
    // so that the localization step is done once for the grid all_points_display
    for (int i = 0; i < var_grid.GetM(); i++)
      var_grid(i).InitGrid(all_points_display);

    // pre-localization of the elements of the mesh so that
    // each new point to search can be fastly found
    // implementation of this method is in GridInterpolation.cxx
    all_points_display.InitInterpolationGrid(this->mesh);

    // implementation of this method in file GridInterpolation.cxx
    all_points_display.LocalizePoints(this->mesh);

    // releasing memory used to prelocalize points
    all_points_display.ClearPrelocalizationArrays();

    // the grids are compressed
    // all_points_display.CompressGrid(var_grid);
  }

  template<class T>
  void WriteOutputFile(const Vector<T>& x, const Vector<T>& y, int num_mode, Complex_wp beta)
  {
    GridInterpolation<Dimension2>& var_interp = all_points_display;


    bool compute_grad = false;
    if (type_coord_other_fields != NONE)
      compute_grad = true;

    Vector<TinyVector<Real_wp, 1> > phi; VectR2 grad_phi, s;
    // boucle sur les grilles
    for (int n = 0; n < this->var_grid.GetM(); n++)
      {
        Vector<T> trace_Ez;
        Vector<T> trace_Hz;
        // si on prend les coordonnées polaires
        Vector<Complex_wp> trace_Hteta;
        Vector<Complex_wp> trace_Er;
        // si on prend les coordonnées cartésiennes
        Vector<Complex_wp> trace_Ex;
        Vector<Complex_wp> trace_Ey;
        Vector<Complex_wp> trace_Hx;
        Vector<Complex_wp> trace_Hy;
        GridInterpolationFull<Dimension2>& var_gr = this->var_grid(n); // grille de rendu final (pas maillage)
        const IVect& list_points = var_gr.GetPointNumber(); // points de la grille
        int nnz = list_points.GetM(); // nombre de points
        trace_Ez.Reallocate(nnz); trace_Ez.Zero();
        trace_Hz.Reallocate(nnz); trace_Hz.Zero();
        trace_Hteta.Reallocate(nnz); trace_Hteta.Zero();
        trace_Er.Reallocate(nnz); trace_Er.Zero();
        trace_Ex.Reallocate(nnz); trace_Ex.Zero();
        trace_Ey.Reallocate(nnz); trace_Ey.Zero();
        trace_Hx.Reallocate(nnz); trace_Hx.Zero();
        trace_Hy.Reallocate(nnz); trace_Hy.Zero();
        // boucle sur les points de la grille
        for (int i1 = 0; i1 < nnz; i1++)
          {
            int i = list_points(i1); // point i
            int iquad = var_interp.GetElementNumber(i); // polygone correspondant dans le maillage
            R2 point_loc = var_interp.GetLocalCoordinate(i); // coordonnées locales
            if ((iquad >= 0) && (iquad < this->mesh.GetNbElt()))
              {
                const IVect& num_ddl = this->mesh_num.Element(iquad).GetNodle();
                this->mesh.GetVerticesElement(iquad, s);
                const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(iquad);

                // we compute the values of basis functions at point_loc
                Fb.ComputeValuesPhi(point_loc, phi, var_interp.GetDFjm1(i),
                                    this->mesh_num, iquad);

                if (compute_grad)
                  Fb.ComputeValuesGradientPhi(point_loc, grad_phi, var_interp.GetDFjm1(i),
                                              this->mesh_num, iquad);

                int ref = this->mesh.Element(iquad).GetReference();
                Complex_wp epsilon = this->ref_epsilon(ref);
                Complex_wp mu = this->ref_mu(ref);
                Complex_wp sigma = this->ref_sigma(ref);
                Complex_wp epsilon_tilde = epsilon + Iwp * sigma / this->omega;
                Complex_wp Delta = this->omega*this->omega*(- epsilon_tilde*mu + beta*beta);
                Complex_wp tempA = - Iwp * this->omega * epsilon_tilde / Delta;
                Complex_wp tempB = - Iwp * beta * this->omega / Delta;
                Complex_wp tempC = - Iwp * this->omega * mu / Delta;
                // cout << "Eps : " << epsilon << endl;
                // cout << "Mu : " << mu << endl;
                // cout << "Sigma : " << sigma << endl;
                // cout << "Eps_tild : " << epsilon_tilde << endl;
                // cout << "Delta : " << Delta << endl;
                // cout << "TempHt : " << tempHt << endl;
                // cout << "TempEr : " << tempEr << endl;

                // attention : pour un cas general il faut calculer Ex et Ey
                // et Hx, Hy
                // On trouve le vect r unitaire
                // DISP(var_interp.GetGlobalCoordinate(i));
                double x1 = var_interp.GetGlobalCoordinate(i)(0);
                double x2 = var_interp.GetGlobalCoordinate(i)(1);
                // DISP(point_loc); DISP(x1); DISP(x2);
                double norm = sqrt(x1*x1 + x2*x2);
                x1 /= norm;
                x2 /= norm;
                // DISP(norm); DISP(x1); DISP(x2);

                T valEz; SetComplexZero(valEz);
                T valHz; SetComplexZero(valHz);
                Complex_wp valEr; SetComplexZero(valEr);
                Complex_wp valHt; SetComplexZero(valHt);
                Complex_wp valEx; SetComplexZero(valEx);
                Complex_wp valEy; SetComplexZero(valEy);
                Complex_wp valHx; SetComplexZero(valHx);
                Complex_wp valHy; SetComplexZero(valHy);
                for (int i = 0; i < num_ddl.GetM(); i++)
                  if (num_ddl(i) >= 0)
                    {

                      // si on se trouve dans l'élément
                      valEz += phi(i)(0)*x(num_ddl(i));
                      valHz += phi(i)(0)*y(num_ddl(i));

                      if (type_coord_other_fields == CARTESIAN)
                        {
                          valEx += grad_phi(i)(0) * x(num_ddl(i)) * tempB;
                          valEy += grad_phi(i)(1)  * x(num_ddl(i)) * tempB;
                          valHx -= grad_phi(i)(1)  * x(num_ddl(i)) * tempA;
                          valHy += grad_phi(i)(0)  * x(num_ddl(i)) * tempA;

                          valEx += grad_phi(i)(1)  * y(num_ddl(i)) * tempC;
                          valEy -= grad_phi(i)(0)  * y(num_ddl(i)) * tempC;
                          valHx += grad_phi(i)(0) * y(num_ddl(i)) * tempB;
                          valHy += grad_phi(i)(1)  * y(num_ddl(i)) * tempB;
                        }
                    }

                trace_Ez(i1) = valEz;
                trace_Hz(i1) = valHz;
                if (type_coord_other_fields == CARTESIAN)
                  {
                    trace_Ex(i1) = valEx;
                    trace_Ey(i1) = valEy;
                    trace_Hx(i1) = valHx;
                    trace_Hy(i1) = valHy;
                  }
              }
          }

        // on ecrit au format lisible par Python (loadND)
        string file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEs" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
        string vtk_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEs" + to_str(num_mode) + "_G" + to_str(n) + ".vtk";
        WriteMatlab(trace_Ez, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        WriteVtk(trace_Ez, var_gr, vtk_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeHs" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
        WriteMatlab(trace_Hz, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);

        if (type_coord_other_fields == CARTESIAN)
          {
            file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEx" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Ex, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeEy" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Ey, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeHx" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Hx, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = PATH_AUGUSTIN + to_str("/Results/")+ "ModeHy" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Hy, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          }
      }
  }


  template<class T>
  void WriteOutputFileScalar(const Vector<T>& u, const string& name)
  {
    GridInterpolation<Dimension2>& var_interp = all_points_display;


    Vector<TinyVector<Real_wp, 1> > phi; VectR2 grad_phi, s;
    // boucle sur les grilles
    for (int n = 0; n < this->var_grid.GetM(); n++)
      {
        Vector<T> trace_u, trace_du_dx, trace_du_dy;
        GridInterpolationFull<Dimension2>& var_gr = this->var_grid(n); // grille de rendu final (pas maillage)
        const IVect& list_points = var_gr.GetPointNumber(); // points de la grille
        int nnz = list_points.GetM(); // nombre de points
        trace_u.Reallocate(nnz); trace_u.Zero();
        trace_du_dx.Reallocate(nnz); trace_du_dx.Zero();
        trace_du_dy.Reallocate(nnz); trace_du_dy.Zero();
        // boucle sur les points de la grille
        for (int i1 = 0; i1 < nnz; i1++)
          {
            int i = list_points(i1); // point i
            int iquad = var_interp.GetElementNumber(i); // polygone correspondant dans le maillage
            R2 point_loc = var_interp.GetLocalCoordinate(i); // coordonnées locales
            if ((iquad >= 0) && (iquad < this->mesh.GetNbElt()))
              {
                const IVect& num_ddl = this->mesh_num.Element(iquad).GetNodle();
                this->mesh.GetVerticesElement(iquad, s);
                const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(iquad);

                // we compute the values of basis functions at point_loc
                Fb.ComputeValuesPhi(point_loc, phi, var_interp.GetDFjm1(i),
                                    this->mesh_num, iquad);

                Fb.ComputeValuesGradientPhi(point_loc, grad_phi, var_interp.GetDFjm1(i),
                                            this->mesh_num, iquad);

                Complex_wp val_u(0, 0);
                R2_Complex_wp grad_u;
                for (int i = 0; i < num_ddl.GetM(); i++)
                  if (num_ddl(i) >= 0)
                    {

                      // si on se trouve dans l'élément
                      val_u += phi(i)(0)*u(num_ddl(i));
                      Add(u(num_ddl(i)), grad_phi(i), grad_u);
                    }

                trace_u(i1) = val_u;
                trace_du_dx(i1) = grad_u(0);
                trace_du_dy(i1) = grad_u(1);
              }
          }

        // on ecrit au format lisible par Python (loadND)
        string root = PATH_AUGUSTIN + to_str("/Results/") + GetBaseString(name);
        string file_name = root + "_G" + to_str(n) + ".dat";
        string file_name_dx = root + "_G" + to_str(n) + "_DX.dat";
        string file_name_dy = root + "_G" + to_str(n) + "_DY.dat";
        WriteMatlab(trace_u, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        WriteMatlab(trace_du_dx, var_gr, file_name_dx, OutputTypeEnum::DOUBLE_PRECISION, false);
        WriteMatlab(trace_du_dy, var_gr, file_name_dy, OutputTypeEnum::DOUBLE_PRECISION, false);
      }
  }

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
    Vector<Matrix<Complexe, Prop, Storage> > vec_Kh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Ch;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Dh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Th;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Uh;
    VectComplex_wp vec_Vk;
    Vector<VectComplex_wp > vec_Vk2;
    VectR2 vec_Pos;
    // Vector<Matrix<Complexe, Symmetric, DiagonalRow>> vec_diag_D_tau0, vec_D_tau0, vec_D_tau0_Bis, vec_ktrSquared;
    Vector<Matrix<Complexe, Symmetric, DiagonalRow>> vec_diag_D_tau0, vec_D_tau0_Bis, vec_ktrSquared;
    Vector<Matrix<Real_wp, Symmetric, DiagonalRow>> vec_diag_D_tau0_Bis;
    Vector<int> nbRefNoeud;
    Matrix<Complexe, Symmetric, DiagonalRow> test, Dtau2, DTau;
    this->ComputeFemMatrices(vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_Vk, vec_Vk2, vec_Pos, vec_diag_D_tau0, vec_diag_D_tau0_Bis);
    // this->ComputeTwistedMaterialsDiagonalMatrix(vec_diag_D_tau0, vec_D_tau0, vec_ktrSquared);
    this->ComputeTwistedMaterialsDiagonalMatrix_Bis(vec_diag_D_tau0_Bis, vec_D_tau0_Bis, vec_ktrSquared, nbRefNoeud);
    Complex_wp zero(0.002,0.004);

    // for (int i = 0 ; i < 20 ; i++) {
    //   DISP(vec_ktrSquared(1).Get(i,i));
    // }
    int N = this->mesh_num.GetNbDof();


    this->ComputeMatrixF2ToF1(zero, vec_D_tau0_Bis, nbRefNoeud, DTau);
    // this->ComputeMatrixV2ToV1(zero, vec_D_tau0_Bis, nbRefNoeud, DTau);
    VectComplex_wp ones(N), res(N);
    res.Zero();
    for (int i = 0 ; i < N ; i++) {
      ones.Get(i) = 1.0;
    }
    Mlt(DTau, ones, res);
    WriteOutputFileScalar(res, "F2ToF1");

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

    for (int ref = 0 ; ref < vec_diag_D_tau0_Bis.GetM() ; ref++) {
      if (index_used(ref)) {
        res.Zero();
        Mlt(vec_ktrSquared(ref), ones, res);
        WriteOutputFileScalar(res, "KTR2_"+ to_str(ref));
        res.Zero();
        Mlt(vec_D_tau0_Bis(ref), ones, res);
        WriteOutputFileScalar(res, "D_tau0_"+ to_str(ref));
      }
    }
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

    if ((type_eigensolver != LINEAR_ARPACK) && (type_eigensolver != LINEAR_SLEPC) && (type_eigensolver != POLYNOMIAL))
      {
        DISP(this->use_split_formulation);
#ifdef SELDON_WITH_SLEPC
        ModeEsNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_D_tau0_Bis, vec_ktrSquared,
                                            this->DofKeptDir, this->IndexDirichlet,
                                            this->use_split_formulation);

        if (this->use_explicit_matrix)
          var_eig.SetExplicitMatrix(); // en calculant les matrices

        var_eig.SetExactPreconditioning();
        var_eig.SetPrintLevel(4);

        var_eig.SetStoppingCriterion(1e-15);
        var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);

        var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
        SlepcParamNep& param = var_eig.GetSlepcParameters();
        if (type_eigensolver == NL_RII)
          param.SetEigensolverType(param.RII);
        else if (type_eigensolver == NL_SLP)
          param.SetEigensolverType(param.SLP);
        else
          param.SetEigensolverType(param.NLEIGS);

        param.SetIntervalRegion(0.1, 10.0, 0.001, 10.0);
        param.EnableCommandLineOptions();
        if (this->use_default_petsc_solver)
          param.SetDefaultPetscSolver();

        Vector<Complex_wp> lambda, lambda_imag;
        Matrix<Complex_wp, General, ColMajor> eigen_vec;

        DISP(var_eig.GetM());
        FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

        DISP(lambda); DISP(lambda_imag);

        for (int k = 0; k < eigen_vec.GetN(); k++)
          {
            int nb_dof_Hs = N;
            VectComplex_wp x(2*nb_dof_Hs), X;
            VectComplex_wp y;
            GetCol(eigen_vec, k, X);
            this->ComputeMatrixV2ToEs(lambda.Get(k), vec_D_tau0_Bis, nbRefNoeud, Dtau2);
            var_eig.ExtractDirichlet(X, x);

            VectComplex_wp E3(nb_dof_Hs); E3.Zero();
            VectComplex_wp H3(nb_dof_Hs); H3.Zero();
            VectComplex_wp V2(nb_dof_Hs);
            VectComplex_wp W2(nb_dof_Hs);
            VectComplex_wp Es(nb_dof_Hs);
            VectComplex_wp Hs(nb_dof_Hs);
            for (int i = 0; i < nb_dof_Hs; i++)
              {
                // Ez(i) = x(i);
                // Hz(i) = x(nb_dof_Hs + i);
                if (this->use_split_formulation)
                  {
                    E3(i) = x(i) * Dtau2.Get(i,i);
                    H3(i) = x(nb_dof_Hs + i) * Dtau2.Get(i,i);
                    // Ez(i) = x(i);
                    // Hz(i) = x(nb_dof_Hs + i);
                    V2(i) = x(i);
                    W2(i) = x(nb_dof_Hs + i);
                  }
                else
                  {
                    E3(i) = x(i);
                    H3(i) = x(nb_dof_Hs + i);
                  }

              }
            //x.Write("Ez_Hz0.dat");

            this->FromF2G2ToEsHs(vec_D_tau0_Bis, nbRefNoeud, vec_Mh, vec_Uh, vec_Th, vec_ktrSquared, lambda.Get(k), V2, W2, Es, Hs);

            ComputeIntensity(Es, vec_Vk, k);
            WriteOutputFileScalar(V2, "F2_"+ to_str(k));
            WriteOutputFileScalar(W2, "G2_"+ to_str(k));
            // WriteOutputFileScalar(E3, "E3_"+ to_str(k));
            // WriteOutputFileScalar(H3, "H3_"+ to_str(k));
            WriteOutputFile(Es, Hs, k, lambda(k));
          }
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
        return;
      }
    else if (type_eigensolver == POLYNOMIAL)
      {
#ifdef SELDON_WITH_SLEPC
        cout << "Résolution polynomiale" << endl;
        DistributedMatrix<Complex_wp, General, ArrayRowSparse> A4, A3, A2, A1, A0;

        // il faut les remplir
        ModeEsPolynomialProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_D_tau0_Bis, vec_ktrSquared,
                                            this->DofKeptDir, this->IndexDirichlet);

        A0 = var_eig.A0;
        A1 = var_eig.A1;
        A2 = var_eig.A2;
        A3 = var_eig.A3;
        A4 = var_eig.A4;

        Vector<DistributedMatrix<Complex_wp, General, ArrayRowSparse>* > list_op(4);
        list_op(0) = &A0;
        list_op(1) = &A1;
        list_op(2) = &A2;
        list_op(3) = &A3;


        var_eig.SetStoppingCriterion(1e-12);
        var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);
        var_eig.InitMatrix(list_op, A4);
        var_eig.SetPrintLevel(6);

        var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, shift);
        var_eig.SetSpectralTransformation(true);

        Vector<Complex_wp> lambda, lambda_imag;
        Matrix<Complex_wp, General, ColMajor> eigen_vec;

        FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

        DISP(lambda); DISP(lambda_imag);

        for (int k = 0; k < eigen_vec.GetN(); k++)
          {
            int nb_dof_Hs = N;
            VectComplex_wp x(2*nb_dof_Hs), X;
            VectComplex_wp y;
            GetCol(eigen_vec, k, X);
            this->ComputeMatrixV2ToEs(lambda.Get(k), vec_D_tau0_Bis, nbRefNoeud, Dtau2);
            var_eig.ExtractDirichlet(X, x);

            VectComplex_wp E3(nb_dof_Hs); E3.Zero();
            VectComplex_wp H3(nb_dof_Hs); H3.Zero();
            VectComplex_wp V2(nb_dof_Hs);
            VectComplex_wp W2(nb_dof_Hs);
            VectComplex_wp Es(nb_dof_Hs);
            VectComplex_wp Hs(nb_dof_Hs);
            for (int i = 0; i < nb_dof_Hs; i++)
              {
                E3(i) = x(i) * Dtau2.Get(i,i);
                H3(i) = x(nb_dof_Hs + i) * Dtau2.Get(i,i);
                // Ez(i) = x(i);
                // Hz(i) = x(nb_dof_Hs + i);
                V2(i) = x(i);
                W2(i) = x(nb_dof_Hs + i);

              }
            //x.Write("Ez_Hz0.dat");

            this->ComputeEsHs(vec_D_tau0_Bis, nbRefNoeud, vec_Mh, vec_Uh, vec_Th, vec_ktrSquared, lambda.Get(k), V2, W2, Es, Hs);

            ComputeIntensity(Es, vec_Vk, k);
            WriteOutputFileScalar(V2, "V2_"+ to_str(k));
            WriteOutputFileScalar(W2, "W2_"+ to_str(k));
            WriteOutputFileScalar(E3, "E3_"+ to_str(k));
            WriteOutputFileScalar(H3, "H3_"+ to_str(k));
            WriteOutputFile(Es, Hs, k, lambda(k));
          }
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
        return;
      }


    // SYSTÈME LINÉAIRE



  // construction des numerotations
  cout << "N = " << N << endl;
  DISP(vec_Pos.GetM());
  int nb_couches = vec_Mh.GetM();
  Vector<int> nb_dof_couche(nb_couches), nb_dof_coucheEs(nb_couches);
  nb_dof_couche.Zero(); nb_dof_coucheEs.Zero();
  Vector<Vector<int> > IndexDofEs(nb_couches), NumGlobEs(nb_couches);
  Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
  for (int ref = 0; ref < nb_couches; ref++)
    if (vec_Mh(ref).GetM() == N)
      {
        IndexDof(ref).Reallocate(N);
        IndexDof(ref).Fill(-1);

        IndexDofEs(ref).Reallocate(N);
        IndexDofEs(ref).Fill(-1);
      }

  // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
  // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
  // IndexDofEs, NumGlobEz pour l'inconnue Es (ddl Dirichlets enleves)
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      int ref = this->mesh.Element(i).GetReference();
      int nb_dof_loc = this->mesh_num.Element(i).GetNbDof();
      for (int j = 0; j < nb_dof_loc; j++)
        {
          int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
          if (IndexDof(ref)(num_dof) == -1) // si c'est un nouveau
            {
              IndexDof(ref)(num_dof) = nb_dof_couche(ref);
              NumGlob(ref).PushBack(num_dof);
              nb_dof_couche(ref)++;
            }

          if ((IndexDirichlet(num_dof) >= 0) && (IndexDofEs(ref)(num_dof) == -1)) // si c'est un nouveau
            {
              IndexDofEs(ref)(num_dof) = nb_dof_coucheEs(ref);
              NumGlobEs(ref).PushBack(num_dof);
              nb_dof_coucheEs(ref)++;
            }
        }
    }

    // Définitions des matrices
    Matrix<Complex_wp, General, ArrayRowSparse> K_tilde;
    Matrix<Complex_wp, General, ArrayRowSparse> M_tilde;

    int nb_dof_Es = DofKeptDir.GetM();
    int nb_dof_Hs = N;

    // Initialisation
    int ordreMatrices = nb_dof_Es + nb_dof_Hs; // lignes pour Es et Hs
    /* lignes pour V_ref et W_ref et leurs variables temporaires
    (pour outrepasser le caractère quadratique de l'équation) */
    for (int ref = 0 ; ref < nb_couches ; ref++)
      ordreMatrices += 4 * nb_dof_couche(ref) + 4 * nb_dof_coucheEs(ref);


    K_tilde.Reallocate(ordreMatrices, ordreMatrices);
    M_tilde.Reallocate(ordreMatrices, ordreMatrices);
    DISP(ordreMatrices);
    K_tilde.Zero(); M_tilde.Zero();

    DISP(this->ref_tau);
    cout << "Création des matrices" << endl;
    int inc = nb_dof_Es + nb_dof_Hs; // premier indice de la première couche dans la matrice
    for (int ref = 0 ; ref < nb_couches ; ref++)
      if (nb_dof_couche(ref) > 0)
        {
          Complex_wp epsilon = this->ref_epsilon(ref);
          Complex_wp mu = this->ref_mu(ref);
          Complex_wp sigma = this->ref_sigma(ref);
          Complex_wp eps_tilde = - Iwp * this->omega * epsilon + sigma;
          Complex_wp mu_tilde = - Iwp * this->omega * mu;
          Complex_wp k2 = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
          Real_wp tau = this->ref_tau(ref);
          // conversion en non-symetrique
          Matrix<Complexe, Prop, Storage> Kh, Mh, Ch, Dh, Th, Uh;
          // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
          Copy(vec_Kh(ref), Kh);
          Copy(vec_Mh(ref), Mh);
          Copy(vec_Ch(ref), Ch);
          Copy(vec_Dh(ref), Dh);
          Copy(vec_Th(ref), Th);
          Copy(vec_Uh(ref), Uh);


          int inc_tempV1 = 0;
          int inc_V1 = nb_dof_coucheEs(ref);
          int inc_tempV2 = 2*nb_dof_coucheEs(ref);
          int inc_V2 = 3*nb_dof_coucheEs(ref);
          int inc_tempW1 = 4 * nb_dof_coucheEs(ref);
          int inc_W1 = 4 * nb_dof_coucheEs(ref) + nb_dof_couche(ref);
          int inc_tempW2 = 4 * nb_dof_coucheEs(ref) + 2*nb_dof_couche(ref);
          int inc_W2 = 4 * nb_dof_coucheEs(ref) + 3*nb_dof_couche(ref);

          // int inc_tempV1 = nb_dof_coucheEs(ref) + nb_dof_couche(ref);
          // int inc_V1 = 0;
          // int inc_W1 = nb_dof_coucheEs(ref);
          // int inc_tempW1 = 2 * nb_dof_coucheEs(ref) + nb_dof_couche(ref);

          for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
            {
              int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet
              int ddl_loc_Es = IndexDofEs(ref)(m);
              R2 pos = vec_Pos(m);
              Real_wp r2 = pos(0)*pos(0) + pos(1)*pos(1);// à compléter
              if (IndexDirichlet(m) != -1) // si on n'est pas sur un noeud de Dirichlet
                {
                  // équation avec la variable temporaire : V_11 = beta omega V_1
                  K_tilde.Get(inc + inc_tempV1 + ddl_loc_Es, inc + inc_tempV1 + ddl_loc_Es) += 1.0;
                  M_tilde.Get(inc + inc_tempV1 + ddl_loc_Es, inc + inc_V1 + ddl_loc_Es) += this->omega;
                  // équation avec la variable temporaire : V_21 = beta omega V_2
                  K_tilde.Get(inc + inc_tempV2 + ddl_loc_Es, inc + inc_tempV2 + ddl_loc_Es) += 1.0;
                  M_tilde.Get(inc + inc_tempV2 + ddl_loc_Es, inc + inc_V2 + ddl_loc_Es) += this->omega;


                  // equation en Es : - Delta_tau,k V_1 + Es = 0
                  //                   k2(1+tau^2 r^2) V_1 + Es = beta^2 omega^2 V1 = beta omega V_11
                  // partie Ez dans l'equation en V_1
                  K_tilde.Get(inc + inc_V1 + ddl_loc_Es, IndexDirichlet(m)) += 1.0;
                  // partie k2(1+tau^2 r^2) V_1 dans l'equation en V_1
                  K_tilde.Get(inc + inc_V1 + ddl_loc_Es, inc + inc_V1 + ddl_loc_Es) += k2 * (1 + tau * tau * r2);
                  // partie beta omega V_11 dans l'equation en V1
                  M_tilde.Get(inc + inc_V1 + ddl_loc_Es, inc + inc_tempV1 + ddl_loc_Es) += this->omega;

                  // equation en V2 : - Delta_tau,k V_2 + V1 = 0
                  //                   k2(1+tau^2 r^2) V_2 + V1 = beta^2 omega^2 V2 = beta omega V_21
                  // partie Es dans l'equation en V_2
                  K_tilde.Get(inc + inc_V2 + ddl_loc_Es, inc + inc_V1 + ddl_loc_Es) += 1.0;
                  // partie k2(1+tau^2 r^2) V_2 dans l'equation en V_2
                  K_tilde.Get(inc + inc_V2 + ddl_loc_Es, inc + inc_V2 + ddl_loc_Es) += k2 * (1 + tau * tau * r2);
                  // partie beta omega V_21 dans l'equation en V2
                  M_tilde.Get(inc + inc_V2 + ddl_loc_Es, inc + inc_tempV2 + ddl_loc_Es) += this->omega;
                }

              // équation avec la variable temporaire : W_11 = beta omega W_1
              K_tilde.Get(inc + inc_tempW1 + ddl_loc, inc + inc_tempW1 + ddl_loc) += 1.0;
              M_tilde.Get(inc + inc_tempW1 + ddl_loc, inc + inc_W1 + ddl_loc) += this->omega;
              // équation avec la variable temporaire : W_21 = beta omega W_2
              K_tilde.Get(inc + inc_tempW2 + ddl_loc, inc + inc_tempW2 + ddl_loc) += 1.0;
              M_tilde.Get(inc + inc_tempW2 + ddl_loc, inc + inc_W2 + ddl_loc) += this->omega;


              // equation en Hs : - Delta_tau,k W_1 + Hs = 0
              //                   k2(1+tau^2 r^2) W_1 + Hs = beta^2 omega^2 W1 = beta omega W_11
              // partie Hs dans l'equation en W1
              K_tilde.Get(inc + inc_W1 + ddl_loc, nb_dof_Es + m) += 1.0;
              // partie k2(1+tau^2 r^2) W_1 dans l'equation en W1
              K_tilde.Get(inc + inc_W1 + ddl_loc, inc + inc_W1 + ddl_loc) += k2 * (1 + tau * tau * r2);
              // partie beta omega W_11 dans l'equation en W1
              M_tilde.Get(inc + inc_W1 + ddl_loc, inc + inc_tempW1 + ddl_loc) += this->omega;



              // equation en W2 : - Delta_tau,k W_2 + W1 = 0
              //                   k2(1+tau^2 r^2) W_2 + W1 = beta^2 omega^2 W2 = beta omega W_21
              // partie W1 dans l'equation en W2
              K_tilde.Get(inc + inc_W2 + ddl_loc, inc + inc_W1 + ddl_loc) += 1.0;
              // partie k2(1+tau^2 r^2) W_2 dans l'equation en W2
              K_tilde.Get(inc + inc_W2 + ddl_loc, inc + inc_W2 + ddl_loc) += k2 * (1.0 + tau * tau * r2);
              // partie beta omega W22 dans l'equation en W2
              M_tilde.Get(inc + inc_W2 + ddl_loc, inc + inc_tempW2 + ddl_loc) += this->omega;



              // équations principales (en V1 et V2)
              // (- i omega epsilon + sigma) M_h E_z + (-i omega epsilon + sigma) K_h V_k = i beta omega C_h W_k
              int size_row = Kh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Kh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_Es = IndexDofEs(ref)(n);

                  if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                    K_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += eps_tilde * Kh.Value(m, j);
                  K_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += mu_tilde * Kh.Value(m, j);
                  if (m == 7)
                  { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); DISP(Kh.Value(m, j)); }
                }

              size_row = Mh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Mh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_Es = IndexDofEs(ref)(n);
                  if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                    {
                      // - eps_tilde k^2 M_h V1 = - eps_tilde beta^2 omega^2 M_h V1 = - eps_tilde beta omega M_h V11
                      K_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += - eps_tilde * k2 * Mh.Value(m, j);
                      M_tilde.Get(IndexDirichlet(m), inc + inc_tempV1 + n_loc_Es) += - eps_tilde * this->omega * Mh.Value(m, j);
                      // - 2 k^4 tau M_h W2 = - 2 k^2 beta^2 omega^2 tau M_h W2 =  2 k^2 beta omega tau M_h W21
                      K_tilde.Get(IndexDirichlet(m), inc + inc_W2 + n_loc) += - 2.0 * k2 * k2 * tau * Mh.Value(m, j);
                      M_tilde.Get(IndexDirichlet(m), inc + inc_tempW2 + n_loc) += - 2.0 * k2 * tau * this->omega * Mh.Value(m, j);
                      // 4 k^4 tau^2 eps_tilde Mh V2 = 0
                      K_tilde.Get(IndexDirichlet(m), inc + inc_V2 + n_loc_Es) += 4.0 * k2 * tau * tau * eps_tilde * Mh.Value(m, j);

                      K_tilde.Get(nb_dof_Es + m, inc + inc_V2 + n_loc_Es) += 2.0 * k2 * k2 * tau * Mh.Value(m, j);
                      M_tilde.Get(nb_dof_Es + m, inc + inc_tempV2 + n_loc_Es) += 2.0 * k2 * tau * this->omega * Mh.Value(m, j);
                    }
                  K_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += - mu_tilde * k2 * Mh.Value(m, j);
                  M_tilde.Get(nb_dof_Es + m, inc + inc_tempW1 + n_loc) += - mu_tilde * this->omega * Mh.Value(m, j);
                  K_tilde.Get(nb_dof_Es + m, inc + inc_W2 + n_loc) += 4.0 * k2 * tau * tau * mu_tilde * Mh.Value(m, j);
                }

              size_row = Ch.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Ch.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_Es = IndexDofEs(ref)(n);

                  if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                    {
                      M_tilde.Get(IndexDirichlet(m), n_loc + inc_W1 + inc) += - Iwp * this->omega * Ch.Value(m, j);
                      M_tilde.Get(nb_dof_Es + m, n_loc_Es + inc_V1 + inc) += Iwp * this->omega * Ch.Value(m, j);
                    }
                }

              size_row = Uh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Uh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_Es = IndexDofEs(ref)(n);

                  if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                    {
                      // 0 = 2 i beta omega tau eps_tilde U_h V1
                      M_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += 2.0 * this->omega * Iwp * tau * eps_tilde * Uh.Value(m,j);
                      // 0 = 2 i beta omega k^2 tau^2 U_h W2
                      M_tilde.Get(IndexDirichlet(m), inc + inc_W2 + n_loc) += 2.0 * this->omega * Iwp * k2 * tau * tau * Uh.Value(m,j);

                      M_tilde.Get(nb_dof_Es + m, inc + inc_V2 + n_loc_Es) += -2.0 * this->omega * Iwp * k2 * tau * tau * Uh.Value(m,j);
                    }
                  M_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += 2.0 * this->omega * Iwp * tau * mu_tilde * Uh.Value(m,j);
                }

              size_row = Th.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Th.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_Es = IndexDofEs(ref)(n);

                  if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                    {
                      // 2 k^2 tau^2 eps_tilde T_h V2 = 0
                      K_tilde.Get(IndexDirichlet(m), inc + inc_V2 + n_loc_Es) += 2.0 * k2 * tau * tau * eps_tilde * Th.Value(m,j);
                    }
                  K_tilde.Get(nb_dof_Es + m, inc + inc_W2 + n_loc) += 2.0 * k2 * tau * tau * mu_tilde * Th.Value(m,j);
                }

              size_row = Dh.GetRowSize(m);
              for (int j = 0 ; j < size_row ; j++)
                {
                  int n = Dh.Index(m, j);
                  int n_loc = IndexDof(ref)(n);
                  int n_loc_Es = IndexDofEs(ref)(n);

                  if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                    {
                      // eps_tilde tau^2 D_h V1 = 0
                      K_tilde.Get(IndexDirichlet(m), inc + inc_V1 + n_loc_Es) += tau * tau * eps_tilde * Dh.Value(m,j);
                    }
                  K_tilde.Get(nb_dof_Es + m, inc + inc_W1 + n_loc) += tau * tau * mu_tilde * Dh.Value(m,j);
                }
            }


          inc += 4 * nb_dof_coucheEs(ref) + 4 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
          // inc += 2 * nb_dof_coucheEs(ref) + 2 * nb_dof_couche(ref);
          //DISP(inc);
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
        ModeEsNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch, vec_Dh, vec_Th, vec_Uh, vec_D_tau0_Bis, vec_ktrSquared,
                                            this->DofKeptDir, this->IndexDirichlet,
                                            this->use_split_formulation);

        // var_eig.MltOperator(betas(i), SeldonNoTrans, X, Y);
        cout << "VP " << i << ": " << betas(i) << endl;

        x.Reallocate(2*nb_dof_Hs);
        var_eig.ExtractDirichlet(X, x);
        for (int i = 0; i < nb_dof_Hs; i++)
          {
            Ez(i) = x(i);
            Hz(i) = x(nb_dof_Hs + i);
          }
        ComputeIntensity(Ez, vec_Vk, i);

        WriteOutputFile(Ez, Hz, i, betas(i));
#endif
      }

  }

};

#ifdef SELDON_WITH_SLEPC
// on inclut les fonctions definies dans le cxx
#include "mode_Es.cxx"
#include "mode_Es_poly.cxx"
#endif

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
