#include "Montjoie.hxx"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;

class ModeEz_Solver;
bool only_Ez = false;

#ifdef SELDON_WITH_SLEPC
// classe pour resoudre un pb aux valeurs propres non-lineaire avec Slepc
// forme splitte : T(L) = \sum_k f_k(L) A_k
// et f_k est une fraction rationnelle
class ModeEzNonLinearEigenProblem : public NonLinearEigenProblem_Base<Complex_wp>
{
  ModeEz_Solver& var; int nb_rows, nb_dof_Ez, nb_dof_Hz, nb_dof_all;
  Vector<Matrix<Complexe, Prop, Storage> >& vec_Mh, &vec_Kh, &vec_Ch;
  SparseDistributedSolver<Complex_wp> mat_lu;
  Vector<int>& DofKeptDir; Vector<int>& IndexDirichlet;

  // pour la formulation splittee
  Matrix<Complexe, Prop, Storage> Mh_sum, Mh_sumH;
  Vector<int> split_reference;

public:
  ModeEzNonLinearEigenProblem(int N, ModeEz_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Ch,
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
#endif

class ModeEz_Solver : public InputDataProblem_Base,
                      public VarFiniteElement<Dimension2>
{
#ifdef SELDON_WITH_SLEPC
  friend class ModeEzNonLinearEigenProblem;
#endif

protected:
  // pulsation omega
  Real_wp omega;

  // bords Dirichlet
  Vector<int> ref_dirichlet;

  // pour les bords quasi-periodiques
  TinyVector<int, 2> ref_periodic;
  int num_mode_periodic;
  Vector<Vector<bool> > is_dof_with_phase;

  // maillage
  Mesh<Dimension2> mesh; string name_mesh;

  // numerotation
  MeshNumbering<Dimension2> mesh_num;
  Vector<int> IndexDirichlet, DofKeptDir;

  // indices physiques
  Vector<Complex_wp> ref_epsilon, ref_mu, ref_sigma;

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

  void SetInputData(const string& keyword, const VectString& parameters)
  {
    // on appelle SetInputData pour mesh et mesh_num
    mesh.SetInputData(keyword, parameters);
    mesh_num.SetInputData(keyword, parameters);

    // autres mots-cles
    if (keyword == "FileMesh")
      name_mesh = parameters(0);
    else if (keyword == "Frequency")
      this->omega = 2.0*pi_wp*to_num<Real_wp>(parameters(0)) + to_num<Real_wp>(parameters(1));
    else if (keyword == "MateriauDielec")
      {
        // on garde le meme format que pour mode_maxwell (parameters(1) pas utilise)
        int ref = to_num<int>(parameters(0));
        ref_epsilon(ref) = to_num<Complex_wp>(parameters(2));
        ref_mu(ref) = to_num<Complex_wp>(parameters(3));
        ref_sigma(ref) = to_num<Complex_wp>(parameters(4));
      }
    else if (keyword == "ReferenceDirichlet")
      {
        // liste des bords de type Dirichlet
        ref_dirichlet.Reallocate(parameters.GetM());
        for (int i = 0; i < parameters.GetM(); i++)
          ref_dirichlet(i) = to_num<int>(parameters(i));
      }
    else if (keyword == "ReferencePeriodic")
      {
        // liste des bords de type periodic
        ref_periodic(0) = to_num<int>(parameters(0));
        ref_periodic(1) = to_num<int>(parameters(1));

        // on rajoute la condition de maillage periodique
        mesh.AddPeriodicCondition(ref_periodic, BoundaryConditionEnum::PERIODIC_THETA);
        mesh_num.SetFormulationForPeriodicCondition(mesh_num.STRONG_PERIODIC);
      }
    else if (keyword == "NumberModePeriodic")
      {
        num_mode_periodic = to_num<int>(parameters(0));
      }
    else if (keyword == "LaplaceEigenvalue")
      {
        if (parameters(0) == "YES")
          compute_laplace_eigen = true;
        else
          compute_laplace_eigen = false;
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
        else if (parameters(0) == "Direct")
          type_eigensolver = DIRECT;
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
    else if (keyword == "OnlyEz")
      {
        if (parameters(0) == "YES")
          only_Ez = true;
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
    else if (keyword == "ReferenceIsolant")
      {
        int N = parameters.GetM();
        ref_vol_isolant.Reallocate(N);
        for (int i = 0; i < N; i++)
            ref_vol_isolant(i) = to_num<int>(parameters(i));
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
    this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, other_mesh_num,
                                     type_integration_edge,
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

    int nb_dof_periodic = this->mesh_num.GetNbPeriodicDof();
    int nb_dof_dirichlet = 0;
    if ((nb_dof_periodic > 0) && (this->num_mode_periodic != 0))
      {
        DofDirichlet(0) = 0;
        nb_dof_dirichlet++;
      }

    for (int i = 0; i < nb_dof_periodic; i++)
      {
        int n = this->mesh_num.GetPeriodicDof(i);
        if (DofDirichlet(n) == -1)
          DofDirichlet(n) = 1;
      }

    cout << "Number of Dirichlet dofs = " << nb_dof_dirichlet << endl;
    cout << "Number of periodic dofs = " << nb_dof_periodic << endl;

    this->is_dof_with_phase.Reallocate(this->mesh.GetNbElt());
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        this->is_dof_with_phase(i).Reallocate(this->mesh_num.Element(i).GetNbDof());
        this->is_dof_with_phase(i).Fill(false);
      }

    // on renumerote
    if ((nb_dof_dirichlet > 0) || (nb_dof_periodic > 0))
      {
        Vector<int> IndexDof(this->mesh_num.GetNbDof());
        IndexDof.Fill(-1);
        int nodl = 0;
        for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
          {
            if (DofDirichlet(i) == -1)
              IndexDof(i) = nodl++;
          }

        for (int i = 0; i < nb_dof_periodic; i++)
          {
            int ndof = this->mesh_num.GetPeriodicDof(i);
            int n0 = this->mesh_num.GetOriginalPeriodicDof(i);
            if (DofDirichlet(ndof) == 1)
              {
                //DISP(i); DISP(n0+1); DISP(ndof+1);
                IndexDof(ndof) = -IndexDof(n0) - 2;
              }
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
                    this->is_dof_with_phase(i)(j) = true;
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

    // on recupere les ddls des conducteurs
    DofDirichlet.Fill(-1);
    ddl_conductor.Reallocate(ref_conductor.GetM());
    Vector<int> index_conductor(this->mesh.GetNbReferences()+1);
    index_conductor.Fill(-1);
    for (int i = 0; i < ref_conductor.GetM(); i++)
      {
        index_conductor(ref_conductor(i)) = i;
        ddl_conductor(i).Clear();
      }

    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
        int ref = this->mesh.BoundaryRef(i).GetReference();
        if (index_conductor(ref) >= 0)
          {
            int ic = index_conductor(ref);
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
                    ddl_conductor(ic).PushBack(num_dof);
                  }
              }
          }
      }

    for (int i = 0; i < ref_conductor.GetM(); i++)
      {
        DISP(ref_conductor(i));
        DISP(ref_vol_conductor(i));
        //DISP(ddl_conductor(i));
      }
  }

  // calcul des matrices elements finis
  void ComputeFemMatrices(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Ch,
                          VectComplex_wp& vec_Vk,
                          Vector<VectComplex_wp >& vec_Vk2)
  {
    // on alloue les matrices
    int N = this->mesh_num.GetNbDof();
    vec_Mh.Reallocate(this->ref_epsilon.GetM());
    vec_Kh.Reallocate(this->ref_epsilon.GetM());
    vec_Ch.Reallocate(this->ref_epsilon.GetM());
    Lk_conductor.Reallocate(this->ref_epsilon.GetM());
    Lk_conductor.Zero();
    vec_Lk.Reallocate(this->ref_epsilon.GetM());
    vec_Lk.Zero();
    vec_Vk.Reallocate(N);
    vec_Vk.Zero();
    vec_Vk2.Reallocate(this->ref_epsilon.GetM());

    Vector<bool> index_used(this->ref_epsilon.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int ref = 0; ref < index_used.GetM(); ref++)
      if (index_used(ref))
        {
          vec_Mh(ref).Reallocate(N, N);
          vec_Kh(ref).Reallocate(N, N);
          if (!only_Ez)
            vec_Ch(ref).Reallocate(N, N);
          vec_Vk2(ref).Reallocate(N);
        }

    VectR2 s;
    SetPoints<Dimension2> PointsElem;
    SetMatrices<Dimension2> MatricesElem;

    // calcul de la phase
    Complexe phase;
    SetComplexOne(phase);
    if (this->mesh_num.GetNbPeriodicDof() > 0)
      {
        Real_wp teta = this->mesh.GetPeriodicAlpha();
        to_complex(exp(Iwp*teta*this->num_mode_periodic), phase);
      }

    Real_wp threshold = 1e-15;

    // boucle sur les elements
    VectReal_wp Bzero, Bmass; Matrix2_2 dfjm1, Ctmp, Cref;
    VectR2 Dzero, Ezero; Vector<Matrix2_2> Astiff, Cstiff;
    TinyVector<bool, 4> null_term(true, false, true, true);
    Matrix<Real_wp> mat_elem_Mh, mat_elem_Kh, mat_elem_Ch;
    Cref(0, 1) = -1.0; Cref(1, 0) = 1.0;
    DISP(this->mesh.GetNbElt());
    num_conduc_ddl.Reallocate(N); num_conduc_ddl.Fill(-1);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);
        int ref = this->mesh.Element(i).GetReference();
        // on recupere le numero du conducteur si present
        int num_conductor = -1;
        for (int k = 0; k < ref_vol_conductor.GetM(); k++)
          if (ref_vol_conductor(k) == ref)
            num_conductor = k;

        // calcul des matrices jacobiennes sur les points de quadrature
        this->mesh.GetVerticesElement(i, s);
        Fb.FjElemQuadrature(s, PointsElem, this->mesh, i);
        Fb.DFjElemQuadrature(s, PointsElem, MatricesElem, this->mesh, i);

        int Nquad = Fb.GetNbPointsQuadratureInside();
        Bmass.Reallocate(Nquad); Astiff.Reallocate(Nquad);
        Cstiff.Reallocate(Nquad);

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
            if (!only_Ez)
              {
                MltTrans(Cref, dfjm1, Ctmp);
                Mlt(dfjm1, Ctmp, Cstiff(j));
                Mlt(jacob*weight, Cstiff(j));
              }
          }

        // calcul des matrices de masse et rigidite elementaires
        mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
        mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();
        if (!only_Ez)
          {
            mat_elem_Ch.Reallocate(nb_dof, nb_dof);
            mat_elem_Ch.Zero();
          }

        Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Astiff, Dzero, Ezero, null_term, mat_elem_Kh);
        if (!only_Ez)
          Fb.AddVariableElemMatrix(0, 0, Bzero, Cstiff, Dzero, Ezero, null_term, mat_elem_Ch);

        // assemblage des matrices
        IVect permut(nb_dof); permut.Fill();
        Sort(num_ddl, permut);

        Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
        col.Fill(-1); val.Zero();
        for (int j = 0; j < nb_dof; j++)
          if (num_ddl(j) >= 0)
            {
              nb_val = 0; Complexe coef1; SetComplexOne(coef1);
              if (this->is_dof_with_phase(i)(permut(j)))
                coef1 *= conjugate(phase);

              //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                  {
                    Complex_wp coef = coef1;
                    if (this->is_dof_with_phase(i)(permut(k)))
                      coef *= phase;

                    col(nb_val) = num_ddl(k);
                    val(nb_val) = coef*mat_elem_Mh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Mh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              nb_val = 0;
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Kh(permut(j), permut(k))) >= threshold))
                  {
                    Complexe coef = coef1;
                    if (this->is_dof_with_phase(i)(permut(k)))
                      coef *= phase;

                    col(nb_val) = num_ddl(k);
                    val(nb_val) = coef*mat_elem_Kh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Kh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

              if (!only_Ez)
                {
                  nb_val = 0;
                  for (int k = 0; k < nb_dof; k++)
                    if ((num_ddl(k) >= 0) && (abs(mat_elem_Ch(permut(j), permut(k))) >= threshold))
                      {
                        Complexe coef = coef1;
                        if (this->is_dof_with_phase(i)(permut(k)))
                          coef *= phase;

                        col(nb_val) = num_ddl(k);
                        val(nb_val) = coef*mat_elem_Ch(permut(j), permut(k));
                        nb_val++;
                      }

                  vec_Ch(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);
                }
            }
      }

    for (int ref = 0; ref < index_used.GetM(); ref++)
      if (index_used(ref))
        {
          if (!only_Ez)
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

    // calcul de la phase
    Complexe phase;
    SetComplexOne(phase);
    if (this->mesh_num.GetNbPeriodicDof() > 0)
      {
        Real_wp teta = this->mesh.GetPeriodicAlpha();
        to_complex(exp(Iwp*teta*this->num_mode_periodic), phase);
        //DISP(teta); DISP(phase);
      }

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
                      Complexe coef; SetComplexOne(coef);
                      if (this->is_dof_with_phase(iquad)(i))
                        coef = phase;

                      // si on se trouve dans l'élément
                      valEz += phi(i)(0)*coef*x(num_ddl(i));
                      if (!only_Ez)
                        valHz += phi(i)(0)*coef*y(num_ddl(i));

                      if (type_coord_other_fields == CARTESIAN)
                        {
                          if (only_Ez)
                            {
                              valEx += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempB;
                              valEy += grad_phi(i)(1) * coef * x(num_ddl(i)) * tempB;
                              valHx -= grad_phi(i)(1) * coef * x(num_ddl(i)) * tempA;
                              valHy += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempA;
                            }
                          else
                            {
                              valEx += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempB;
                              valEy += grad_phi(i)(1) * coef * x(num_ddl(i)) * tempB;
                              valHx -= grad_phi(i)(1) * coef * x(num_ddl(i)) * tempA;
                              valHy += grad_phi(i)(0) * coef * x(num_ddl(i)) * tempA;

                              valEx += grad_phi(i)(1) * coef * y(num_ddl(i)) * tempC;
                              valEy -= grad_phi(i)(0) * coef * y(num_ddl(i)) * tempC;
                              valHx += grad_phi(i)(0) * coef * y(num_ddl(i)) * tempB;
                              valHy += grad_phi(i)(1) * coef * y(num_ddl(i)) * tempB;
                            }
                        }
                      else if (type_coord_other_fields == POLAR)
                        {
                          double dPhi_dr = x1 * grad_phi(i)(0) + x2 * grad_phi(i)(1);
                          valHt += dPhi_dr * coef * x(num_ddl(i)) * tempA;
                          valEr += dPhi_dr * coef * x(num_ddl(i)) * tempB;
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
                else if (type_coord_other_fields == POLAR)
                  {
                    trace_Hteta(i1) = valHt;
                    trace_Er(i1) = valEr;
                  }
              }
          }

        // on ecrit au format lisible par Python (loadND)
        string file_name = "ModeEz" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
        WriteMatlab(trace_Ez, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
        if (!only_Ez)
          {
            file_name = "ModeHz" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Hz, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          }

        if (type_coord_other_fields == CARTESIAN)
          {
            file_name = "ModeEx" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Ex, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = "ModeEy" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Ey, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = "ModeHx" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Hx, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = "ModeHy" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Hy, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          }
        else if (type_coord_other_fields == POLAR)
          {
            file_name = "ModeEr" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Er, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
            file_name = "ModeHt" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
            WriteMatlab(trace_Hteta, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
          }
      }
  }


  template<class T>
  void WriteOutputFileScalar(const Vector<T>& u, const string& name)
  {
    GridInterpolation<Dimension2>& var_interp = all_points_display;

    // calcul de la phase
    Complexe phase;
    SetComplexOne(phase);
    if (this->mesh_num.GetNbPeriodicDof() > 0)
      {
        Real_wp teta = this->mesh.GetPeriodicAlpha();
        to_complex(exp(Iwp*teta*this->num_mode_periodic), phase);
        //DISP(teta); DISP(phase);
      }

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
                      Complexe coef; SetComplexOne(coef);
                      if (this->is_dof_with_phase(iquad)(i))
                        coef = phase;

                      // si on se trouve dans l'élément
                      val_u += phi(i)(0)*coef*u(num_ddl(i));
                      Add(coef*u(num_ddl(i)), grad_phi(i), grad_u);
                    }

                trace_u(i1) = val_u;
                trace_du_dx(i1) = grad_u(0);
                trace_du_dy(i1) = grad_u(1);
              }
          }

        // on ecrit au format lisible par Python (loadND)
        string root = GetBaseString(name);
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

    // calcul de la grille d'interpolation (pour sorties des vecteurs propres)
    this->InitGrid();

    // calcul des matrices de rigidite et de la matrice de masse
    Vector<Matrix<Complexe, Prop, Storage> > vec_Mh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Kh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Ch;
    VectComplex_wp vec_Vk;
    Vector<VectComplex_wp > vec_Vk2;
    this->ComputeFemMatrices(vec_Mh, vec_Kh, vec_Ch, vec_Vk, vec_Vk2);

    int N = this->mesh_num.GetNbDof();

    // calcul des petites valeurs propres du laplacian
    if (compute_laplace_eigen )
      {
        SparseEigenProblem<Complexe, Matrix<Complexe, Prop, Storage>,
                           Matrix<Complexe, Prop, Storage> > var_eig;

        var_eig.SetStoppingCriterion(1e-12);
        var_eig.SetNbAskedEigenvalues(3);
        var_eig.SetComputationalMode(var_eig.INVERT_MODE);

        var_eig.InitMatrix(vec_Kh(1), vec_Mh(1));

        var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.02, var_eig.SORTED_MODULUS);

        cout << "On calcule les valeurs propres/vecteurs propres" << endl;
        Vector<Complexe> lambda, lambda_imag;
        Matrix<Complexe, General, ColMajor> eigen_vec;
        GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);

        DISP(lambda);
        // on peut ecrire les vecteurs propres
        Vector<Complexe> x;
        VectComplex_wp y;
        for (int i = 0; i < lambda.GetM(); i++)
          {
            GetCol(eigen_vec, i, x);
            WriteOutputFile(x, y, i, 0);
          }

        return;
      }


    if ((type_eigensolver != LINEAR_ARPACK) && (type_eigensolver != LINEAR_SLEPC) && (type_eigensolver != DIRECT))
      {
        DISP(this->use_split_formulation);
#ifdef SELDON_WITH_SLEPC
        ModeEzNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
                                            this->DofKeptDir, this->IndexDirichlet,
                                            this->use_split_formulation);

        if (this->use_explicit_matrix)
          var_eig.SetExplicitMatrix(); // en calculant les matrices

        var_eig.SetExactPreconditioning();
        var_eig.SetPrintLevel(4);

        var_eig.SetStoppingCriterion(1e-12);
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
            int nb_dof_Hz = N;
            VectComplex_wp x(2*nb_dof_Hz), X;
            VectComplex_wp y;
            GetCol(eigen_vec, k, X);
            var_eig.ExtractDirichlet(X, x);
            VectComplex_wp Ez(nb_dof_Hz); Ez.Zero();
            VectComplex_wp Hz(nb_dof_Hz); Hz.Zero();
            for (int i = 0; i < nb_dof_Hz; i++)
              {
                Ez(i) = x(i);
                if (!only_Ez)
                  Hz(i) = x(nb_dof_Hz + i);
              }
            //x.Write("Ez_Hz0.dat");

            ComputeIntensity(Ez, vec_Vk, k);
            WriteOutputFile(Ez, Hz, k, lambda(k));
          }
#else
        cout << "Recompile with Slepc" << endl;
        abort();
#endif
        return;
      }

    // construction des numerotations
    cout << "N = " << N << endl;
    int nb_couches = vec_Mh.GetM();
    Vector<int> nb_dof_couche(nb_couches), nb_dof_coucheEz(nb_couches);
    nb_dof_couche.Zero(); nb_dof_coucheEz.Zero();
    Vector<Vector<int> > IndexDofEz(nb_couches), NumGlobEz(nb_couches);
    Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
    for (int ref = 0; ref < nb_couches; ref++)
      if (vec_Mh(ref).GetM() == N)
        {
          IndexDof(ref).Reallocate(N);
          IndexDof(ref).Fill(-1);

          IndexDofEz(ref).Reallocate(N);
          IndexDofEz(ref).Fill(-1);
        }

    // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
    // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
    // IndexDofEz, NumGlobEz pour l'inconnue Ez (ddl Dirichlets enleves)
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

            if ((IndexDirichlet(num_dof) >= 0) && (IndexDofEz(ref)(num_dof) == -1)) // si c'est un nouveau
              {
                IndexDofEz(ref)(num_dof) = nb_dof_coucheEz(ref);
                NumGlobEz(ref).PushBack(num_dof);
                nb_dof_coucheEz(ref)++;
              }
          }
      }


    if (type_eigensolver != DIRECT)
      {
        // Définitions des matrices
        Matrix<Complex_wp, General, ArrayRowSparse> K_tilde;
        Matrix<Complex_wp, General, ArrayRowSparse> M_tilde;

        int nb_dof_Ez = DofKeptDir.GetM();
        int nb_dof_Hz = N;

        if (only_Ez)
          {
            int ordreMatrices = nb_dof_Ez; // lignes pour Ez
            for (int ref = 0 ; ref < nb_couches ; ref++)
              ordreMatrices += nb_dof_coucheEz(ref); // lignes pour V_ref

            K_tilde.Reallocate(ordreMatrices, ordreMatrices);
            M_tilde.Reallocate(ordreMatrices, ordreMatrices);
            DISP(ordreMatrices);
            K_tilde.Zero(); M_tilde.Zero();

            cout << "Création des matrices" << endl;
            int inc = nb_dof_Ez; // premier indice de la première couche dans la matrice
            for (int ref = 0 ; ref < nb_couches ; ref++)
              if (nb_dof_couche(ref) > 0)
              {
                Complex_wp epsilon = this->ref_epsilon(ref);
                Complex_wp mu = this->ref_mu(ref);
                Complex_wp sigma = this->ref_sigma(ref);
                // conversion en non-symetrique
                Matrix<Complexe, Prop, Storage> Kh, Mh;
                // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
                Copy(vec_Kh(ref), Kh);
                Copy(vec_Mh(ref), Mh);
                for (int ddl_loc = 0 ; ddl_loc < nb_dof_coucheEz(ref) ; ddl_loc++)
                  {
                    int m = NumGlobEz(ref)(ddl_loc);
                    // partie -Ez dans l'equation en V_k
                    K_tilde.Get(inc + ddl_loc, IndexDirichlet(m)) -= 1.0;

                    // partie (-omega^2 epsilon - i omega sigma) mu V_k dans l'equation en V_k
                    K_tilde.Get(inc + ddl_loc, inc + ddl_loc) -= this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;

                    // partie -beta^2 omega^2 V_k dans l'equation en Vk
                    M_tilde.Get(inc + ddl_loc, inc + ddl_loc) -= this->omega * this->omega;
                    // cout << K_tilde.Get(inc + ddl_loc, m) << " ," << K_tilde.Get(inc + ddl_loc + 1, m) << endl;

                    // équation en Ez
                    int size_row = Kh.GetRowSize(m);
                    for (int j = 0 ; j < size_row ; j++)
                      {
                        int n = Kh.Index(m, j);
                        int n_loc = IndexDofEz(ref)(n);
                        if (IndexDirichlet(n) != -1)
                          K_tilde.Get(IndexDirichlet(m), n_loc + inc) += (- Iwp * this->omega * epsilon + sigma) * Kh.Value(m, j);
                        //if (m == 7)
                        // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
                      }

                    size_row = Mh.GetRowSize(m);
                    for (int j = 0 ; j < size_row ; j++)
                      {
                        int n = Mh.Index(m, j);
                        if (IndexDirichlet(n) != -1)
                          K_tilde.Get(IndexDirichlet(m), IndexDirichlet(n)) += (- Iwp * this->omega * epsilon + sigma) * Mh.Value(m, j);
                      }
                  }

                // cout << "Couche " << ref << endl;
                // for (int i = inc ; i <= inc + 1 ; i++) {
                //     for (int j = 0 ; j <= N ; j++) {
                //         cout << K_tilde.Get(i, j) << " ";
                //     }
                //     cout << endl;
                // }

                inc += nb_dof_coucheEz(ref); //on arrive au premier indice de la couche suivante dans la matrice
                //DISP(inc);
              }
          }
        else
          {
            // Initialisation
            int ordreMatrices = nb_dof_Ez + nb_dof_Hz; // lignes pour Ez et Hz
            /* lignes pour V_ref et W_ref et leurs variables temporaires
            (pour outrepasser le caractère quadratique de l'équation) */
            for (int ref = 0 ; ref < nb_couches ; ref++)
              ordreMatrices += 2 * nb_dof_couche(ref) + 2 * nb_dof_coucheEz(ref);


            K_tilde.Reallocate(ordreMatrices, ordreMatrices);
            M_tilde.Reallocate(ordreMatrices, ordreMatrices);
            DISP(ordreMatrices);
            K_tilde.Zero(); M_tilde.Zero();

            cout << "Création des matrices" << endl;
            int inc = nb_dof_Ez + nb_dof_Hz; // premier indice de la première couche dans la matrice
            for (int ref = 0 ; ref < nb_couches ; ref++)
              if (nb_dof_couche(ref) > 0)
                {
                  Complex_wp epsilon = this->ref_epsilon(ref);
                  Complex_wp mu = this->ref_mu(ref);
                  Complex_wp sigma = this->ref_sigma(ref);
                  // conversion en non-symetrique
                  Matrix<Complexe, Prop, Storage> Kh, Mh, Ch;
                  // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
                  Copy(vec_Kh(ref), Kh);
                  Copy(vec_Mh(ref), Mh);
                  Copy(vec_Ch(ref), Ch);

                  int inc_tempVk = 0;
                  int inc_Vk = nb_dof_coucheEz(ref);
                  int inc_tempWk = 2 * nb_dof_coucheEz(ref);
                  int inc_Wk = 2 * nb_dof_coucheEz(ref) + nb_dof_couche(ref);

                  for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
                    {
                      int m = NumGlob(ref)(ddl_loc); // numéro global avec les noeuds de Dirichlet
                      int ddl_loc_Ez = IndexDofEz(ref)(m);
                      if (IndexDirichlet(m) != -1) // si on n'est pas sur un noeud de Dirichlet
                        {
                          // équation avec la variable temporaire : TV_k = beta omega V_k
                          K_tilde.Get(inc + inc_tempVk + ddl_loc_Ez, inc + inc_tempVk + ddl_loc_Ez) = 1.0;
                          M_tilde.Get(inc + inc_tempVk + ddl_loc_Ez, inc + inc_Vk + ddl_loc_Ez) = this->omega;

                          // equation en Vk : - Delta_k V_k + E_z = 0
                          //                   (omega^2 epsilon mu + i omega sigma mu) V_k + Ez = beta^2 omega^2 Vk = beta omega TV_k
                          // partie Ez dans l'equation en V_k
                          K_tilde.Get(inc + inc_Vk + ddl_loc_Ez, IndexDirichlet(m)) = 1.0;
                          // partie (omega^2 epsilon mu + i omega sigma) mu V_k dans l'equation en V_k
                          K_tilde.Get(inc + inc_Vk + ddl_loc_Ez, inc + inc_Vk + ddl_loc_Ez) = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
                          // partie beta omega tempV_k dans l'equation en Vk
                          M_tilde.Get(inc + inc_Vk + ddl_loc_Ez, inc + inc_tempVk + ddl_loc_Ez) = this->omega;
                        }

                      // equation en TWk : TW_k = beta omega W_k
                      // équation avec la variable temporaire
                      K_tilde.Get(inc + inc_tempWk + ddl_loc, inc + inc_tempWk + ddl_loc) = 1.0;
                      M_tilde.Get(inc + inc_tempWk + ddl_loc, inc + inc_Wk + ddl_loc) = this->omega;

                      // equation en Wk : - Delta_k W_k + H_z = 0
                      //                   (omega^2 epsilon mu + i omega sigma mu) W_k + Ez = beta^2 omega^2 Wk = beta omega TW_k
                      // partie Hz dans l'equation en W_k
                      K_tilde.Get(inc + inc_Wk + ddl_loc, nb_dof_Ez + m) = 1.0;
                      // partie (omega^2 epsilon mu + i omega sigma) mu W_k dans l'equation en W_k
                      K_tilde.Get(inc + inc_Wk + ddl_loc, inc + inc_Wk + ddl_loc) = this->omega * this->omega * epsilon * mu + Iwp * this->omega * sigma * mu;
                      // partie beta omega tempW_k dans l'equation en Wk
                      M_tilde.Get(inc + inc_Wk + ddl_loc, inc + inc_tempWk + ddl_loc) = this->omega;


                      // équations en Ez et Hz
                      // (- i omega epsilon + sigma) M_h E_z + (-i omega epsilon + sigma) K_h V_k = i beta omega C_h W_k
                      int size_row = Kh.GetRowSize(m);
                      for (int j = 0 ; j < size_row ; j++)
                        {
                          int n = Kh.Index(m, j);
                          int n_loc = IndexDof(ref)(n);
                          int n_loc_Ez = IndexDofEz(ref)(n);

                          if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                            K_tilde.Get(IndexDirichlet(m), inc + inc_Vk + n_loc_Ez) = (- Iwp * this->omega * epsilon + sigma) * Kh.Value(m, j);
                          K_tilde.Get(nb_dof_Ez + m, inc + inc_Wk + n_loc) = (- Iwp * this->omega * mu) * Kh.Value(m, j);
                          //if (m == 7)
                          // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
                        }

                      size_row = Mh.GetRowSize(m);
                      for (int j = 0 ; j < size_row ; j++)
                        {
                          int n = Mh.Index(m, j);
                          if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                            K_tilde.Get(IndexDirichlet(m), IndexDirichlet(n)) += (- Iwp * this->omega * epsilon + sigma) * Mh.Value(m, j);
                          K_tilde.Get(nb_dof_Ez + m, nb_dof_Ez + n) += (- Iwp * this->omega * mu) * Mh.Value(m, j);
                        }

                      size_row = Ch.GetRowSize(m);
                      for (int j = 0 ; j < size_row ; j++)
                        {
                          int n = Ch.Index(m, j);
                          int n_loc = IndexDof(ref)(n);
                          int n_loc_Ez = IndexDofEz(ref)(n);

                          if ((IndexDirichlet(m) != -1) && (IndexDirichlet(n) != -1))
                            {
                              M_tilde.Get(IndexDirichlet(m), n_loc + inc_Wk + inc) = - Iwp * this->omega * Ch.Value(m, j);
                              M_tilde.Get(nb_dof_Ez + m, n_loc_Ez + inc_Vk + inc) = Iwp * this->omega * Ch.Value(m, j);
                            }
                        }
                    }


                  inc += 2 * nb_dof_coucheEz(ref) + 2 * nb_dof_couche(ref); //on arrive au premier indice de la couche suivante dans la matrice
                  //DISP(inc);
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
            ModeEzNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
                                                this->DofKeptDir, this->IndexDirichlet);

            var_eig.MltOperator(betas(i), SeldonNoTrans, X, Y);
            cout << "VP " << i << ": " << betas(i) << " ; norme : " << Norm2(Y) / Norm2(X) << endl;

            x.Reallocate(2*nb_dof_Hz);
            var_eig.ExtractDirichlet(X, x);
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
    else // si on est dans le cas direct
      {


        Matrix<Complex_wp, General, ArrayRowSparse> Ae, Ah;
        int ordreMatrices = N; // lignes pour Ez
        Ae.Reallocate(ordreMatrices, ordreMatrices);
        for (int ref = 0 ; ref < ref_conductor.GetM() ; ref++)
            ordreMatrices++;

        Ah.Reallocate(ordreMatrices, ordreMatrices);
        Ae.Zero(); Ah.Zero();

        cout << "Création des matrices" << endl;
        for (int ref = 0 ; ref < nb_couches ; ref++)
          if (nb_dof_couche(ref) > 0)
            {
              cout << "On est dans la couche " << ref << "." << endl;
              int num_conductor = -1;
              for (int k = 0; k < ref_vol_conductor.GetM(); k++)
                if (ref_vol_conductor(k) == ref)
                  num_conductor = k;

              Complex_wp epsilon = this->ref_epsilon(ref);
              Complex_wp mu = this->ref_mu(ref);
              Complex_wp sigma = this->ref_sigma(ref);
              // conversion en non-symetrique
              Matrix<Complexe, Prop, Storage> Kh, Mh;
              // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
              Copy(vec_Kh(ref), Kh);
              Copy(vec_Mh(ref), Mh);

              // Ae : -div((epsilon + i sigma/omega) grad phi) = 0
              //DISP(ref); DISP((epsilon + Iwp * sigma / this->omega));
              Add((epsilon + Iwp * sigma / this->omega), Kh, Ae);
              
              // pour Ah, deux termes
              // (1) : - div(mu^-1 grad(A)) - i omega sigma (A - Mk(A)) = Ik / |Lambda_k|
              // (2) : - sum Ai (Vk)i + wk = 0
              Complex_wp coef_sigma = this->omega*this->omega*epsilon + Iwp * this->omega * sigma;
              for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
                {
                  int m = NumGlob(ref)(ddl_loc);
                  int size_row = Kh.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Kh.Index(m, j);
                      Ah.Get(m, n) +=  Kh.Value(m, j) / mu;
                    }

                  size_row = Mh.GetRowSize(m);
                  for (int j = 0 ; j < size_row ; j++)
                    {
                      int n = Mh.Index(m, j);
                      Ah.Get(m, n) -=  coef_sigma * Mh.Value(m, j);
                    }

                  if (num_conductor >= 0)
                    {
                      Ah.Get(m, N + num_conductor) +=  Iwp * this->omega * sigma * vec_Vk(m) * Lk_conductor(num_conductor);
                      Ah.Get(N + num_conductor, m) -=  vec_Vk(m);
                    }
                }
            }

        // on impose Dirichlet pour Ae => phi = 0
        // Dirichlet sur le bord exterieur (pas conducteur)
        for (int i = 0; i < N; i++)
          if (IndexDirichlet(i) == -1)
            {
              Ae.ReallocateRow(i, 1); // cette ligne ne contient qu'une entrée non nulle
              Ae.Index(i, 0) = i; // la première entrée non nulle est à la colonne (réelle) i (sur la diagonale donc)
              Ae.Value(i, 0) = 1.0; // et elle vaut 1

              Ah.ReallocateRow(i, 1); // cette ligne ne contient qu'une entrée non nulle
              Ah.Index(i, 0) = i; // la première entrée non nulle est à la colonne (réelle) i (sur la diagonale donc)
              Ah.Value(i, 0) = 1.0; // et elle vaut 1
            }

        // pour Ae :
        // on met des zeros sur les lignes associees aux ddls conducteurs
        // sauf 1 sur la diagonale
        // et on met la phase sur le second membre => phi = V

        // pour Ah dans les conducteurs, terme w_k (deuxieme equation)
        VectComplex_wp rhsE(Ae.GetM());
        VectComplex_wp rhsH(Ah.GetM());
        rhsE.Zero(); rhsH.Zero();
        DISP(phase_conductor);
        for (int i = 0; i < ref_conductor.GetM(); i++)
          {
            DISP(i);
            int ref = ref_conductor(i);
            DISP(ref);
            cout << "On est dans la couche conducteur " << ref << "." << endl;
            Complex_wp sigma = this->ref_sigma(ref);

            // terme w_k pour A_h
            Ah.Get(N + i, N + i) += 1.0;

            // boucle sur les ddls surfaciques du conducteur (pour Ae)
            for (int j = 0; j < ddl_conductor(i).GetM(); j++)
              {
                // condition phi = V (pour Ae)
                int row = ddl_conductor(i)(j);
                Ae.ReallocateRow(row, 1);
                Ae.Index(row, 0) = row;
                Ae.Value(row, 0) = 1.0;
                rhsE(row) = phase_conductor(i);
              }
          }

        // calcul de B_i = V_{k, i} * phase (Intensite de 1A)
        for (int i = 0; i < N; i++)
          if (num_conduc_ddl(i) >= 0)
            {
              int j = num_conduc_ddl(i);
              rhsH(i) = phase_conductor(j) * vec_Vk(i);
            }

        SparseDistributedSolver<Complex_wp> mat_lu;
        mat_lu.Factorize(Ae);
        mat_lu.Solve(rhsE);
        mat_lu.Factorize(Ah);
        mat_lu.Solve(rhsH);
        cout << "Valeurs de w_k = " << endl;
        for (int i = 0; i < ref_conductor.GetM(); i++)
          cout << "Pour le conducteur " << i << " = " << rhsH(N+i) << endl;

        WriteOutputFileScalar(rhsE, "Phi");
        WriteOutputFileScalar(rhsH, "A");
      }

    return;
  }

};

#ifdef SELDON_WITH_SLEPC
// on inclut les fonctions definies dans le cxx
#include "mode_Ez.cxx"
#endif

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
  var.RunAll(input_file);

  return FinalizeMontjoie();
}
