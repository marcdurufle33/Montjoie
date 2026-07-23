#include "Montjoie.hxx"
#include <ctime>
#include <cmath>

#define PATH_AUGUSTIN "./src/Program/Augustin"

using namespace Montjoie;

// Complex_wp pour condition quasi-periodique
typedef Complex_wp Complexe;
typedef General Prop;
typedef ArrayRowSparse Storage;

class ModeU_Solver;

#ifdef SELDON_WITH_SLEPC
// classe pour resoudre un pb aux valeurs propres non-lineaire avec Slepc
// forme splitte : T(L) = \sum_k f_k(L) A_k
// et f_k est une fraction rationnelle
class ModeUNonLinearEigenProblem : public NonLinearEigenProblem_Base<Complex_wp>
{
  ModeU_Solver& var; int nb_dof;
  Vector<Matrix<Complexe, Prop, Storage> >& vec_Mh, &vec_Kh;
  Matrix<Complexe, Prop, Storage>& Kh_bord;
  VectReal_wp& Mh_bord;

  Matrix<Complexe, Prop, Storage> A0, A1, A2, A3;

  SparseDistributedSolver<Complex_wp> mat_lu;
  // Vector<int>& ptsBord;

  // pour la formulation splittee
  Matrix<Complexe, Prop, Storage> Mh_sum;
  Vector<int> split_reference;
  int type_CLA;
  enum {NEUMANN, BGT1, BGT2};

public:
  ModeUNonLinearEigenProblem(int N, ModeU_Solver& var0,
                              Vector<Matrix<Complexe, Prop, Storage> >& Mh,
                              Vector<Matrix<Complexe, Prop, Storage> >& Kh,
                              Matrix<Complexe, Prop, Storage>& Kh_surf,
                              VectReal_wp& Mh_surf,
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
  void MltOperatorSplit(int i, const SeldonTranspose&, const Vector<Complex_wp>& U, Vector<Complex_wp>& Y);

};
#endif

class ModeU_Solver : public InputDataProblem_Base,
                      public VarFiniteElement<Dimension2>
{
#ifdef SELDON_WITH_SLEPC
  friend class ModeUNonLinearEigenProblem;
#endif

protected:
  // pulsation omega
  Real_wp omega;

  // bords Dirichlet
  // CE NE SONT PLUS DES DIRICHLETS : IL FAUT LES INTÉGRER DANS DES MATRICES SURFACIQUES
  int ref_bord;
  Real_wp ray;

  // maillage
  Mesh<Dimension2> mesh; string name_mesh;

  // numerotation
  MeshNumbering<Dimension2> mesh_num;
  Vector<int> IndexDirichlet, DofKeptDir;

  // indices physiques
  Vector<Complex_wp> ref_rho;
  Complex_wp rhoInfini;

  // type de condition au bord
  int type_CLA;
  enum {NEUMANN, BGT1, BGT2};
  bool linear;
  int linearization_step;


  // objet pour localiser les points sur le maillage
  GridInterpolation<Dimension2> all_points_display;

  // liste des grilles de sortie
  Vector<GridInterpolationFull<Dimension2> > var_grid;

  // type de solveur valeurs propres
  int type_eigensolver;
  enum {LINEAR_ARPACK, LINEAR_SLEPC, NL_RII, NL_SLP, NL_EIGS};
  bool use_split_formulation, use_explicit_matrix, use_default_petsc_solver;

  // type de coordonnées pour les champs transversaux
  int type_coord_other_fields;
  enum {CARTESIAN, POLAR, NONE};

  // m maximal pour l'orthogonalisation
  int mMax = -1;
  int refOrtho = -1;
  bool orthogonalisation = false;

  // shift pour le calcul des valeurs propres
  Complex_wp shift; int nb_asked_eigenval;
  Real_wp minRe, maxRe, minIm, maxIm;


public:
  ModeU_Solver() : mesh_num(mesh)
  { omega = 2.0*pi_wp;
    type_CLA = NEUMANN; linear = true;
    type_eigensolver = LINEAR_ARPACK; shift = Complex_wp(1.6,0.0);
    type_coord_other_fields = NONE;
    nb_asked_eigenval = 10;
    use_split_formulation = false; use_explicit_matrix = false;
    use_default_petsc_solver = false;
    minRe = 0.1;
    maxRe = 10.0;
    minIm = -0.5;
    maxIm = 0.5;
  }

  void InitIndices(int ref_max)
  {
    ref_rho.Reallocate(ref_max+1);
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
      {
        // this->omega = 2.0*pi_wp*to_num<Real_wp>(parameters(0)) + to_num<Real_wp>(parameters(1));
        // cout << "Pulsation w = " << this->omega << endl;
        this->omega = to_num<Real_wp>(parameters(0)) + to_num<Real_wp>(parameters(1));
      }
    else if (keyword == "MateriauDielec")
      {
        // on garde le meme format que pour mode_maxwell (parameters(1) pas utilise)
        int ref = to_num<int>(parameters(0));
        ref_rho(ref) = to_num<Complex_wp>(parameters(2));
        if (parameters.GetM() == 3)
          {
            ref_rho(ref) = to_num<Complex_wp>(parameters(2));
          }
        else if (parameters.GetM() >= 4)
          {
            ref_rho(ref) = to_num<Complex_wp>(parameters(2)) + Iwp * to_num<Complex_wp>(parameters(3));
          }
      }
    else if (keyword == "Ray")
      {
        ray = to_num<Real_wp>(parameters(0));
      }
    else if (keyword == "ReferenceBord")
      {
        // liste des bords de type Dirichlet
        ref_bord = to_num<int>(parameters(0));
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
    else if (keyword == "Condition")
    {
      cout << "Condition : " << parameters(0);
      if (parameters(0) == "Neumann")
        cout << endl;
      else
        {
          cout << " " << parameters(1);
          if (parameters(1) == "Linear")
            cout << " " << parameters(2) << endl;
          else
            cout << endl;
        }

      if (parameters(0) == "Neumann")
        {
          type_CLA = NEUMANN;
          linear = true;
        }
      else if (parameters(0) == "BGT1")
        type_CLA = BGT1;
      else if (parameters(0) == "BGT2")
        type_CLA = BGT2;
      else
        {
          cout << "Unknown condition " << parameters(0) << endl;
          abort();
        }

      if (type_CLA != NEUMANN)
        {
          if (parameters(1) == "Linear")
            {
              linear = true;
              linearization_step = to_num<int>(parameters(2));
              type_eigensolver = LINEAR_ARPACK;
            }
          else if (parameters(1) == "Non_Linear")
            {
              linear = false;
              type_eigensolver = NL_EIGS;
              use_split_formulation = true;
            }
          else
            {
              cout << "Unknown condition resolution " << parameters(1) << endl;
              abort();
            }
        }
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
    else if (keyword == "Intervalle")
      {
        minRe = to_num<Real_wp>(parameters(0));
        maxRe = to_num<Real_wp>(parameters(1));
        minIm = to_num<Real_wp>(parameters(2));
        maxIm = to_num<Real_wp>(parameters(3));
        cout << "Intervalle de recherche : [" << minRe << ";" << maxRe << "] + i [" << minIm << ";" << maxIm << "]" << endl;
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
    else if (keyword == "Orthogonalization")
      {
        mMax = to_num<int>(parameters(0));
        orthogonalisation = true;

        if (parameters.GetM() > 1)
          {
            refOrtho = to_num<int>(parameters(1));
            cout << "Orthogonalisation sur la couche " << refOrtho << ", mMax = " << mMax << endl;
          }
        else
          {
            refOrtho = -1;
            cout << "Orthogonalisation sur le bord, mMax = " << mMax << endl;
          }
        cout << "Coefs max orthogonalization : " << mMax << endl;
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

  // calcule le theta correspondant à un point de coordonnées globales
  Real_wp atan2(R2 point)
    {
      Real_wp theta;
      if ((point(0) == 0) and (point(1) == 0)) // si on est au centre du domaine
        {
          return 0.0;
        }
      else if (point(0) >= 0)
        {
          return atan(point(1)/point(0));
        }
      else if (point(0) < 0)
        {
          return atan(point(1)/point(0)) + pi_wp;
        }
    }

  // calcul des matrices elements finis
  void ComputeFemMatrices(Vector<Matrix<Complex_wp, Prop, Storage> >& vec_Mh,
                          Vector<Matrix<Complexe, Prop, Storage> >& vec_Kh,
                          int ref_surf,
                          VectReal_wp& Mh_surf,
                          Matrix<Complexe, Prop, Storage>& Kh_surf,
                          // VectReal_wp& vec_theta,
                          VectReal_wp& vec_theta_Omega
                          // Matrix<Complexe, Prop, Storage>& full_Mh
                        )
  {
    // on alloue les matrices
    int N = this->mesh_num.GetNbDof();
    vec_Mh.Reallocate(this->ref_rho.GetM());
    vec_Kh.Reallocate(this->ref_rho.GetM());
    vec_theta_Omega.Reallocate(N);
    vec_theta_Omega.Zero();

    Vector<bool> index_used(this->ref_rho.GetM());
    index_used.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      index_used(this->mesh.Element(i).GetReference()) = true;

    for (int ref = 0; ref < index_used.GetM(); ref++)
      if (index_used(ref))
        {
          vec_Mh(ref).Reallocate(N, N);
          vec_Kh(ref).Reallocate(N, N);
        }

    VectR2 s;
    SetPoints<Dimension2> PointsElem;
    SetMatrices<Dimension2> MatricesElem;

    // calcul de la phase
    Complexe phase;
    SetComplexOne(phase);

    Real_wp threshold = 1e-15;

    // boucle sur les elements (calcul des matrices volumiques)
    VectReal_wp Bzero, Bmass; Matrix2_2 dfjm1, Ctmp, Cref;
    VectR2 Dzero, Ezero; Vector<Matrix2_2> Astiff, Cstiff;
    TinyVector<bool, 4> null_term(true, false, true, true);
    Matrix<Real_wp> mat_elem_Mh, mat_elem_Kh, mat_elem_Ch;
    Cref(0, 1) = -1.0; Cref(1, 0) = 1.0;
    DISP(this->mesh.GetNbElt());
    cout << "Calcul matrices volumiques" << endl;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(i);
        int ref = this->mesh.Element(i).GetReference();

        // calcul des matrices jacobiennes sur les points de quadrature
        this->mesh.GetVerticesElement(i, s);
        Fb.FjElem(s, PointsElem, this->mesh, i);
        Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);

        int nb_dof = Fb.GetNbDof();
        IVect num_ddl = this->mesh_num.Element(i).GetNodle();
        int Nquad = Fb.GetNbPointsQuadratureInside();
        Bmass.Reallocate(Nquad); Astiff.Reallocate(Nquad);
        Cstiff.Reallocate(Nquad);

        // boucle sur les points de quadrature
        for (int j = 0; j < Nquad; j++)
          {
            // pour le calcul de Mh (matrice de masse)
            Real_wp weight = Fb.WeightsND(j);
            Real_wp jacob = Det(MatricesElem.GetPointQuadrature(j));
            Bmass(j) = jacob*weight;

            // calcul de Vk
            // pour avoir phi_i(xi_j) appeler Fb.GetValuePhiOnQuadraturePoint(j, phi);


            // pour le calcul de Kh (matrice de rigidite)
            GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);
            MltTrans(dfjm1, dfjm1, Astiff(j));
            Mlt(jacob*weight, Astiff(j));

          }

        // boucle sur les ddls (supposes nodaux)
        for (int k = 0; k < Fb.GetNbPointsDof(); k++)
          {
            if (num_ddl(k) >= 0)
              {
                R2 pointQuadGlob = PointsElem.GetPointDof(k);
                vec_theta_Omega(num_ddl(k)) = atan2(pointQuadGlob);
              }
          }

        // calcul des matrices de masse et rigidite elementaires
        mat_elem_Mh.Reallocate(nb_dof, nb_dof); mat_elem_Mh.Zero();
        mat_elem_Kh.Reallocate(nb_dof, nb_dof); mat_elem_Kh.Zero();

        Fb.AddVariableMassMatrix(0, 0, Bmass, mat_elem_Mh);
        Fb.AddVariableElemMatrix(0, 0, Bzero, Astiff, Dzero, Ezero, null_term, mat_elem_Kh);

        // assemblage des matrices
        IVect permut(nb_dof); permut.Fill();
        Sort(num_ddl, permut);

        Vector<int> col(nb_dof); Vector<Complexe> val(nb_dof); int nb_val = 0;
        col.Fill(-1); val.Zero();
        for (int j = 0; j < nb_dof; j++)
          if (num_ddl(j) >= 0)
            {
              nb_val = 0; Complexe coef1; SetComplexOne(coef1);

              //DISP(i); DISP(j); DISP(coef1); DISP(num_ddl(j));
              for (int k = 0; k < nb_dof; k++)
                if ((num_ddl(k) >= 0) && (abs(mat_elem_Mh(permut(j), permut(k))) >= threshold))
                  {
                    Complex_wp coef = coef1;

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

                    col(nb_val) = num_ddl(k);
                    val(nb_val) = coef*mat_elem_Kh(permut(j), permut(k));
                    nb_val++;
                  }

              vec_Kh(ref).AddInteractionRow(num_ddl(j), nb_val, col, val, true);

            }
      }

    // boucle sur les aretes (matrices surfaciques)
    Mh_surf.Reallocate(N);
    Mh_surf.Zero();
    Kh_surf.Clear(); Kh_surf.Reallocate(N, N);
    Kh_surf.Zero();
    // vec_theta.Reallocate(N);
    // vec_theta.Zero();
    cout << "Calcul matrices surfaciques" << endl;
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
        int ref = this->mesh.BoundaryRef(i).GetReference();
        if (ref != ref_surf)
          continue;

        int num_elem = this->mesh.BoundaryRef(i).numElement(0);
        IVect num_ddl = this->mesh_num.Element(num_elem).GetNodle();
        int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
        const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);

        this->mesh.GetVerticesElement(num_elem, s);
        Fb.FjElemQuadrature(s, PointsElem, this->mesh, num_elem);
        Fb.DFjElemQuadrature(s, PointsElem, MatricesElem, this->mesh, num_elem);

        // restriction of DF_i on surface
        Fb.FjSurfaceElem(s, PointsElem, mesh, num_elem, num_loc);
        Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, mesh, num_elem, num_loc);

        // on calcule les gradients surfaciques des fonctions de base sur les points de quadrature
        // et la matrice de masse surfacique
        int nb_points_quadrature_edge = Fb.GetNbQuadBoundary(num_loc);
        int Ns = Fb.GetNbDofBoundary(num_loc);
        Matrix<R2> grad_phi_surf(Ns, nb_points_quadrature_edge);
        VectR2 grad_phi; R2 grad;
        for (int k = 0; k < nb_points_quadrature_edge; k++)
          {
            // jacobien surfacique
            Real_wp ds = MatricesElem.GetDsQuadratureBoundary(k);
            // poids d'integration sur l'intervalle unite
            Real_wp poids = Fb.WeightsQuadratureBoundary(k, num_loc);
            int num_dof_loc = Fb.GetLocalNumber(num_loc, k);
            int num_dof = num_ddl(num_dof_loc);
            Mh_surf(num_dof) += poids*ds;

            // // point de quadrature
            // // DISP(PointsElem.GetPointQuadratureBoundary(k));
            // R2 pointQuadGlob = PointsElem.GetPointQuadratureBoundary(k);
            // vec_theta(num_dof) = atan2(pointQuadGlob);

            // normale sur le point
            R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);

            // on calcule grad(phi) sur le point de quadrature
            GetInverse(MatricesElem.GetPointQuadratureBoundary(k), dfjm1);
            int num_point = Fb.GetQuadNumber(num_loc, k);
            Fb.GetGradientPhiOnQuadraturePoint(num_point, grad_phi);
            for (int i = 0; i < Ns; i++)
              {
                num_dof_loc = Fb.GetLocalNumber(num_loc, i);
                MltTrans(dfjm1, grad_phi(num_dof_loc), grad);
                Real_wp grad_dot_n = DotProd(grad, normale);
                grad -= grad_dot_n*normale;
                grad_phi_surf(i, k) = grad;
              }
          }

        // on calcule la matrice de rigidite surfacique
        for (int k = 0; k < nb_points_quadrature_edge; k++)
          {
            // jacobien surfacique
            Real_wp ds = MatricesElem.GetDsQuadratureBoundary(k);
            // poids d'integration sur l'intervalle unite
            Real_wp poids = Fb.WeightsQuadratureBoundary(k, num_loc);

            for (int i0 = 0; i0 < Ns; i0++)
              for (int j = 0; j < Ns; j++)
                {
                  int num_dof_loc_i0 = Fb.GetLocalNumber(num_loc, i0);
                  int num_dof_i0 = num_ddl(num_dof_loc_i0);
                  int num_dof_loc_j = Fb.GetLocalNumber(num_loc, j);
                  int num_dof_j = num_ddl(num_dof_loc_j);
                  Real_wp vloc = DotProd(grad_phi_surf(i0, k), grad_phi_surf(j, k));
                  Kh_surf.AddInteraction(num_dof_i0, num_dof_j, vloc*poids*ds);
                }
          }
      }

    // DISP(Mh_surf); Kh_surf.WriteText("KhSurf.dat");
    //exit(0);
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
  void WriteOutputFile(const Vector<T>& u, int num_mode, Complex_wp beta)
  {
    GridInterpolation<Dimension2>& var_interp = all_points_display;

    bool compute_grad = false;
    if (type_coord_other_fields != NONE)
      compute_grad = true;

    Vector<TinyVector<Real_wp, 1> > phi; VectR2 grad_phi, s;
    // boucle sur les grilles
    for (int n = 0; n < this->var_grid.GetM(); n++)
      {
        Vector<T> trace_u;
        GridInterpolationFull<Dimension2>& var_gr = this->var_grid(n); // grille de rendu final (pas maillage)
        const IVect& list_points = var_gr.GetPointNumber(); // points de la grille
        int nnz = list_points.GetM(); // nombre de points
        trace_u.Reallocate(nnz); trace_u.Zero();
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


                T valU; SetComplexZero(valU);
                for (int i = 0; i < num_ddl.GetM(); i++)
                  if (num_ddl(i) >= 0)
                    {
                      // si on se trouve dans l'élément
                      valU += phi(i)(0)*u(num_ddl(i));
                    }

                trace_u(i1) = valU;
              }
          }

        // on ecrit au format lisible par Python (loadND)
        string file_name = PATH_AUGUSTIN + to_str("/Results/") + "ModeU" + to_str(num_mode) + "_G" + to_str(n) + ".dat";
        cout << file_name << endl;
        WriteMatlab(trace_u, var_gr, file_name, OutputTypeEnum::DOUBLE_PRECISION, false);
      }
  }

  void ComputeCoefficientsDOmega(const VectComplex_wp& u, const VectReal_wp& Mh_surf, const VectReal_wp& theta, const int& k,
                            VectComplex_wp& coefsModes)
  {
    Real_wp threshold = 1e-15;
    VectComplex_wp coefs(2*mMax+1); coefs.Zero();
    Real_wp absMaxCoefs = 0.0;
    int mMode = - mMax - 1;
    for (int m = - mMax ; m <= mMax ; m++)
      {
        for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
          {
            coefs(m + mMax) += u(i) * exp(-Iwp * double(m) * theta(i)) * Mh_surf(i);
          }
        if (absMaxCoefs < abs(coefs(m + mMax)))
          {
            absMaxCoefs = abs(coefs(m + mMax));
            mMode = m;
          }
      }

    if (absMaxCoefs >= threshold)
      {
        if (mMode == 0)
          {
            coefsModes.Reallocate(1);
            coefsModes.Get(0) = coefs(mMax);
          }
        else
          {
            coefsModes.Reallocate(2);
            coefsModes.Get(0) = coefs(mMax - abs(mMode));
            coefsModes.Get(1) = coefs(mMax + abs(mMode));
          }
      }
    else
      {
        cout << "!!! Attention, mMax est sans doute trop petit !!!";
      }
    cout << "Coefficients of the mode " << k << " of m " << mMode << " : " << coefsModes << endl;
  }

  void ComputeCoefficientsFullCable(const VectComplex_wp& u, const Matrix<Complexe, Prop, Storage>& Mh, const VectReal_wp& theta, const int& k,
                            VectComplex_wp& coefsModes)
  {
    Real_wp threshold = 1e-15;
    VectComplex_wp coefs(2*mMax+1); coefs.Zero();
    Real_wp absMaxCoefs = 0.0;
    int mMode = - mMax - 1;
    for (int m = - mMax ; m <= mMax ; m++)
      {
        for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
          {
            int size_row = Mh.GetRowSize(i);
            if (size_row > 0) // on est dans le cable
              {
                for (int j = 0 ; j < size_row ; j++)
                  {
                    // c'est pas j !!!
                    int jGlob = Mh.Index(i, j);
                    coefs(m + mMax) += u(jGlob) * exp(-Iwp * double(m) * theta(jGlob)) * Mh.Value(i, j);
                  }
              }
          }
        // for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
        //   {
        //     int i = NumGlob(ref)(ddl_loc);
        //     int size_row = Mh.GetRowSize(i);
        //     for (int j = 0 ; j < size_row ; j++)
        //       {
        //         coefs(m + mMax) += u(j) * exp(-Iwp * double(m) * theta(j)) * Mh.Value(i, j);
        //       }
        //   }

        if (absMaxCoefs < abs(coefs(m + mMax)))
          {
            absMaxCoefs = abs(coefs(m + mMax));
            mMode = m;
          }
      }

    if (absMaxCoefs >= threshold)
      {
        if (mMode == 0)
          {
            coefsModes.Reallocate(1);
            coefsModes.Get(0) = coefs(mMax);
          }
        else
          {
            coefsModes.Reallocate(2);
            coefsModes.Get(0) = coefs(mMax - abs(mMode));
            coefsModes.Get(1) = coefs(mMax + abs(mMode));
          }
      }
    else
      {
        cout << "!!! Attention, mMax est sans doute trop petit !!!";
      }

    cout << "Coefficients of the mode " << k << " of m " << mMode << " : " << coefsModes << endl;
  }

  // fonction principale pour lancer la simulation
  void RunAll(const string& input_file)
  {
    InitIndices(50);

    // Gauss-Lobatto finite element
    string name_elt = "QUADRANGLE_LOBATTO";

    // on lit le fichier de donnees
    ReadInputFile(input_file, *this);

    if (linear && (type_eigensolver != LINEAR_ARPACK) && (type_eigensolver != LINEAR_SLEPC))
      {
        cout << "Illogical ini file : linear with " << type_eigensolver << endl;
        abort();
      }
    else if (!linear && ((type_eigensolver == LINEAR_ARPACK) || (type_eigensolver == LINEAR_SLEPC)))
      {
        cout << "Illogical ini file : non linear with " << type_eigensolver << endl;
        abort();
      }


    // on construit le maillage
    this->mesh.Read(this->name_mesh);

    // on construit l'element fini
    this->ConstructFiniteElement(name_elt);

    // on numerote le maillage
    this->mesh_num.NumberMesh();
    cout << "Number of degrees of freedom = " << this->mesh_num.GetNbDof() << endl;

    // calcul des ddls de Dirichlet (pour les eliminer de la numerotation)
    this->mesh.Write("test.mesh");
    MontjoieTimer chrono;
    clock_t debut, fin;


    // calcul de la grille d'interpolation (pour sorties des vecteurs propres)
    this->InitGrid();

    chrono.SetTimer(MontjoieTimer::ACCURATE_TIMER);
    chrono.Start("PVP");
    debut = clock();

    // calcul des matrices de rigidite et de la matrice de masse
    Vector<Matrix<Complexe, Prop, Storage> > vec_Mh;
    Vector<Matrix<Complexe, Prop, Storage> > vec_Kh;
    VectReal_wp Mh_surf;
    Matrix<Complexe, Prop, Storage> Kh_surf;
    // VectReal_wp vec_theta;
    VectReal_wp vec_theta_Omega;
    DISP(ref_bord);
    this->ComputeFemMatrices(vec_Mh, vec_Kh, ref_bord, Mh_surf, Kh_surf, vec_theta_Omega);



    int N = this->mesh_num.GetNbDof();


    int nb_dof_border = 0;
    Vector<int> IndexLocalBorder, NumGlobBorder;
    IndexLocalBorder.Reallocate(N);
    IndexLocalBorder.Fill(-1);
    rhoInfini = Complex_wp(0.0, 0.0);

    // IndexLocalBorder(j) renvoie le numero local du ddl global j dans la frontière
    // NumGlobBorder(i) renvoie le numero global j du ddl local i dans la frontière
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
        int ref = this->mesh.BoundaryRef(i).GetReference();
        if (ref == ref_bord)
          {
            int num_elem = this->mesh.BoundaryRef(i).numElement(0);
            int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
            const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
            const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
            if (realpart(rhoInfini) == 0.0)
              {
                int refElem = this->mesh.Element(num_elem).GetReference();
                rhoInfini = this->ref_rho(refElem);
              }
            int nb_dof = Fb.GetNbDofBoundary(num_loc);
            for (int j = 0; j < nb_dof; j++)
              {
                int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
                int num_dof = Nodle(num_dof_loc);
                // cout << "Arête " << i << ", ddl local " << j << endl;
                // DISP(num_dof_loc);
                // DISP(num_dof);
                if (IndexLocalBorder(num_dof) == -1)
                  {
                    IndexLocalBorder(num_dof) = nb_dof_border;
                    NumGlobBorder.PushBack(num_dof);
                    nb_dof_border++;
                  }
              }
          }
      }
    DISP(rhoInfini);


    if ((type_eigensolver != LINEAR_ARPACK) && (type_eigensolver != LINEAR_SLEPC))
      {
        cout << "Partie non linéaire" << endl;
        DISP(this->use_split_formulation);
#ifdef SELDON_WITH_SLEPC
        ModeUNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, Kh_surf, Mh_surf,
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

        param.SetIntervalRegion(minRe, maxRe, minIm, maxIm);
        param.EnableCommandLineOptions();
        if (this->use_default_petsc_solver)
          param.SetDefaultPetscSolver();

        Vector<Complex_wp> lambda, lambda_imag;
        Matrix<Complex_wp, General, ColMajor> eigen_vec;

        DISP(var_eig.GetM());
        FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);

        chrono.Stop("PVP");

        cout << "Temps de résolution : " << chrono.GetSeconds("PVP") << endl;

        DISP(lambda); DISP(lambda_imag);

        bool premierModeDoubleFait = false;
        Complex_wp a1, a2, b1, b2;
        VectComplex_wp u1, u2;

        for (int k = 0; k < eigen_vec.GetN(); k++)
          {
            int nb_dof = N;
            VectComplex_wp x;
            GetCol(eigen_vec, k, x);
            VectComplex_wp U1(nb_dof); U1.Zero();
            VectComplex_wp U2(nb_dof); U2.Zero();
            VectComplex_wp coefsModes;

            if (orthogonalisation)
              {
                if (this->refOrtho == -1)
                  {
                    ComputeCoefficientsDOmega(x, Mh_surf, vec_theta_Omega, k, coefsModes);
                  }
                else
                  {
                    ComputeCoefficientsFullCable(x, vec_Mh.Get(this->refOrtho), vec_theta_Omega, k, coefsModes);
                  }
                if (coefsModes.GetM() == 1)
                  {
                    for (int i = 0; i < nb_dof; i++)
                      {
                        U1(i) = x(i) / coefsModes.Get(0);
                      }

                    WriteOutputFile(U1, k, lambda(k));
                  }
                else if (coefsModes.GetM() == 2)
                  {
                    if (!premierModeDoubleFait)
                      {
                        a1 = coefsModes.Get(0);
                        a2 = coefsModes.Get(1);
                        u1 = x;
                        premierModeDoubleFait = true;
                      }
                    else
                      {
                        b1 = coefsModes.Get(0);
                        b2 = coefsModes.Get(1);
                        u2 = x;
                        premierModeDoubleFait = false;

                        Complex_wp Det = a1 * b2 - a2 * b1;

                        for (int i = 0; i < nb_dof; i++)
                          {
                            U1(i) = (b2 * u1(i) - a2 * u2(i)) / Det;
                            U2(i) = (-b1 * u1(i) + a1 * u2(i)) / Det;
                          }

                        WriteOutputFile(U1, k-1, lambda(k-1));
                        WriteOutputFile(U2, k, lambda(k));
                      }
                  }
                else
                  {
                    cout << "Problème orthogonalisation..." << endl;
                    WriteOutputFile(x, k, lambda(k));
                  }
              }
            else
              {
                WriteOutputFile(x, k, lambda(k));
              }
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
    Vector<int> nb_dof_couche(nb_couches);
    nb_dof_couche.Zero();
    Vector<Vector<int> > IndexDof(nb_couches), NumGlob(nb_couches);
    for (int ref = 0; ref < nb_couches; ref++)
      if (vec_Mh(ref).GetM() == N)
        {
          IndexDof(ref).Reallocate(N);
          IndexDof(ref).Fill(-1);
        }

    // IndexDof(ref)(j) renvoie le numero local du ddl global j dans la couche ref
    // NumGlob(ref)(i) renvoie le numero global j du ddl local i dans la couche ref
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
          }
      }


    // DISP(IndexLocalBorder);
    // DISP(NumGlobBorder);
    // DISP(nb_dof_border);

    // Définitions des matrices
    Matrix<Complex_wp, General, ArrayRowSparse> K_tilde;
    Matrix<Complex_wp, General, ArrayRowSparse> M_tilde;

    int nb_dof_u = N;
    int nbVarAuxiliaires = 0;

    // On cherche le nb de variables auxiliaires
    if ((type_CLA == BGT1 && linearization_step == 2) || (type_CLA == BGT2 && linearization_step == 1))
      nbVarAuxiliaires = 1;
    else if (type_CLA == BGT2 && linearization_step == 2)
      nbVarAuxiliaires = 2;
    else if (type_CLA == BGT2 && linearization_step == 3)
      nbVarAuxiliaires = 4;

    int ordreMatrices = nb_dof_u + nbVarAuxiliaires * nb_dof_border;
    // for (int ref = 0 ; ref < nb_couches ; ref++)
    //   ordreMatrices += nb_dof_coucheEz(ref); // lignes pour V_ref

    K_tilde.Reallocate(ordreMatrices, ordreMatrices);
    M_tilde.Reallocate(ordreMatrices, ordreMatrices);
    DISP(ordreMatrices);
    K_tilde.Zero(); M_tilde.Zero();

    cout << "Création des matrices" << endl;
    // int inc = nb_dof_u; // premier indice de la première couche dans la matrice

    // on écrit le problème sur Omega
    Complex_wp rhoPlus = 0;
    for (int ref = 0 ; ref < nb_couches ; ref++)
      if (nb_dof_couche(ref) > 0)
      {
        Complex_wp rho = this->ref_rho(ref);
        if (realpart(rho) > realpart(rhoPlus))
            rhoPlus = rho;
        cout << ref << " : " << rho << endl;
        // conversion en non-symetrique
        Matrix<Complexe, Prop, Storage> Kh, Mh;
        // Matrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
        Copy(vec_Kh(ref), Kh);
        Copy(vec_Mh(ref), Mh);
        // équation en Ez
        for (int ddl_loc = 0 ; ddl_loc < nb_dof_couche(ref) ; ddl_loc++)
          {
            int m = NumGlob(ref)(ddl_loc);
            int size_row = Kh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Kh.Index(m, j);
                K_tilde.Get(m, n) += Kh.Value(m, j);
                //if (m == 7)
                // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
              }

            size_row = Mh.GetRowSize(m);
            for (int j = 0 ; j < size_row ; j++)
              {
                int n = Mh.Index(m, j);
                K_tilde.Get(m, n) -= rho * this->omega * this->omega * Mh.Value(m, j);
                // M_tilde.Get(m, n) -= this->omega * this->omega * Mh.Value(m, j);
                M_tilde.Get(m, n) -= Mh.Value(m, j) * this->omega * this->omega;
              }


          }
        // cout << "Couche " << ref << endl;
        // for (int i = inc ; i <= inc + 1 ; i++) {
        //     for (int j = 0 ; j <= N ; j++) {
        //         cout << K_tilde.Get(i, j) << " ";
        //     }
        //     cout << endl;
        // }

        // inc += nb_dof_coucheEz(ref); //on arrive au premier indice de la couche suivante dans la matrice
        //DISP(inc);
      }

    // cout << "ATTENTION RHO_PLUS MANUEL !!!" << endl;
    // rhoPlus = 10.0;

    // partie sur le bord
    if (this->type_CLA != NEUMANN)
      {
        DISP(rhoInfini);
        Complex_wp alpha = this->omega * Sqrt(rhoPlus - rhoInfini);
        Real_wp R = this->ray;
        for (int ddl_loc_border = 0 ; ddl_loc_border < nb_dof_border ; ddl_loc_border++)
          {
            int m = NumGlobBorder(ddl_loc_border);

            if (this->type_CLA == BGT1)
              {
                if (this->linearization_step == 0)
                  K_tilde.Get(m, m) += (alpha + 1/(2*R)) * Mh_surf(m);
                else if (this->linearization_step == 1)
                  {
                    K_tilde.Get(m, m) += (alpha/2.0 - rhoInfini * this->omega * this->omega/(2.0*alpha) + 1/(2*R)) * Mh_surf(m);
                    M_tilde.Get(m, m) -= Mh_surf(m) / (2.0 * alpha) * this->omega * this->omega;
                  }
                else if (this->linearization_step == 2)
                  {
                    int mV = nb_dof_u + ddl_loc_border;
                    // partie heut gauche (eq de u en fonction de u)
                    K_tilde.Get(m, m) += (alpha/4.0 - rhoInfini * this->omega * this->omega/(4.0*alpha) + 1/(2*R)) * Mh_surf(m);
                    M_tilde.Get(m, m) -= Mh_surf(m) / (4.0 * alpha) * this->omega * this->omega;
                    // partie haut droite (eq de u en fonction de v)
                    K_tilde.Get(m, mV) -= rhoInfini * this->omega * this->omega * Mh_surf(m);
                    M_tilde.Get(m, mV) -= Mh_surf(m) * this->omega * this->omega;
                    // partie bas gauche (eq de v en fonction de u)
                    K_tilde.Get(mV, m) += 1.0;
                    // partie bas droite (eq de v en fonction de v)
                    K_tilde.Get(mV, mV) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                    M_tilde.Get(mV, mV) += 1.0 / alpha * this->omega * this->omega;
                  }
              }
            else if (this->type_CLA == BGT2)
              {
                if (this->linearization_step == 0)
                  {
                    K_tilde.Get(m, m) += (3/(2*R) - (9/(8*R) + R * rhoInfini * this->omega * this->omega)/(1.0 + R * alpha)) * Mh_surf(m);
                    M_tilde.Get(m, m) -= R * Mh_surf(m) / (1.0 + R * alpha) * this->omega * this->omega;

                    int size_row = Kh_surf.GetRowSize(m);
                    for (int j = 0 ; j < size_row ; j++)
                      {
                        int n = Kh_surf.Index(m, j);
                        K_tilde.Get(m, n) += R * Kh_surf.Value(m, j) / (2.0 * (1.0 + R * alpha));
                        //if (m == 7)
                        // { DISP(ref); DISP(vec_Kh(ref).Index(m, j)); DISP(n_loc+inc); DISP( K_tilde.Get(m, n_loc + inc)); }
                      }
                  }
                else if (this->linearization_step == 1)
                  {
                    int mV = nb_dof_u + ddl_loc_border;
                    // partie heut gauche (eq de u en fonction de u)
                    K_tilde.Get(m, m) += 3/(2*R) * Mh_surf(m);
                    // partie haut droite (eq de u en fonction de v)
                    K_tilde.Get(m, mV) -= (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                    int size_row = Kh_surf.GetRowSize(m);
                    for (int j = 0 ; j < size_row ; j++)
                      {
                        int n = Kh_surf.Index(m, j);
                        int nV = nb_dof_u + IndexLocalBorder(n);
                        K_tilde.Get(m, nV) += R * Kh_surf.Value(m, j) / 2.0;
                      }
                    M_tilde.Get(m, mV) -= R * Mh_surf(m) * this->omega * this->omega;
                    // partie bas gauche (eq de v en fonction de u)
                    K_tilde.Get(mV, m) += 1.0;
                    // partie bas droite (eq de v en fonction de v)
                    K_tilde.Get(mV, mV) -= 1.0 + alpha * R / 2.0 - R * rhoInfini * this->omega * this->omega / (2.0*alpha);
                    M_tilde.Get(mV, mV) += R / (2.0*alpha) * this->omega * this->omega;
                  }
                else if (this->linearization_step == 2)
                  {
                    int mV = nb_dof_u + ddl_loc_border;
                    int mW = nb_dof_u + nb_dof_border + ddl_loc_border;
                    // partie heut gauche (eq de u en fonction de u)
                    K_tilde.Get(m, m) += 3/(2*R) * Mh_surf(m);
                    // partie haut milieu (eq de u en fonction de v)
                    K_tilde.Get(m, mV) -= (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                    int size_row = Kh_surf.GetRowSize(m);
                    for (int j = 0 ; j < size_row ; j++)
                      {
                        int n = Kh_surf.Index(m, j);
                        int nV = nb_dof_u + IndexLocalBorder(n);
                        K_tilde.Get(m, nV) += R * Kh_surf.Value(m, j) / 2.0;
                      }
                    M_tilde.Get(m, mV) -= R * Mh_surf(m) * this->omega * this->omega;
                    // partie milieu gauche (eq de v en fonction de u)
                    K_tilde.Get(mV, m) += 1.0;
                    // partie milieu milieu (eq de v en fonction de v)
                    K_tilde.Get(mV, mV) -= 1.0 + alpha * R / 4.0 - R * rhoInfini * this->omega * this->omega / (4.0*alpha);
                    M_tilde.Get(mV, mV) += R / (4.0*alpha) * this->omega * this->omega;
                    // partie milieu droite (eq de v en fonction de w)
                    K_tilde.Get(mV, mW) += R * rhoInfini * this->omega * this->omega;
                    M_tilde.Get(mV, mW) += R * this->omega * this->omega;
                    // partie bas milieu (eq de w en fonction de v)
                    K_tilde.Get(mW, mV) += 1.0;
                    // partie bas droite (eq de v en fonction de w)
                    K_tilde.Get(mW, mW) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                    M_tilde.Get(mW, mW) += 1.0 / alpha * this->omega * this->omega;
                  }
                else if (this->linearization_step == 3)
                  {
                    int mU1 = nb_dof_u + ddl_loc_border;
                    int mV = nb_dof_u + nb_dof_border + ddl_loc_border;
                    int mW = nb_dof_u + 2 * nb_dof_border + ddl_loc_border;
                    int mX = nb_dof_u + 3 * nb_dof_border + ddl_loc_border;
                    // partie 1 1 (eq de u en fonction de u)
                    K_tilde.Get(m, m) += 3/(2*R) * Mh_surf(m);
                    // partie 1 2 (eq de u en fonction de u')
                    K_tilde.Get(m, mU1) -= (9/(8*R) + R * rhoInfini * this->omega * this->omega) * Mh_surf(m);
                    int size_row = Kh_surf.GetRowSize(m);
                    for (int j = 0 ; j < size_row ; j++)
                      {
                        int n = Kh_surf.Index(m, j);
                        int nU1 = nb_dof_u + IndexLocalBorder(n);
                        K_tilde.Get(m, nU1) += R * Kh_surf.Value(m, j) / 2.0;
                      }
                    M_tilde.Get(m, mU1) -= R * Mh_surf(m) * this->omega * this->omega;

                    // partie 2 1 (eq de u' en fonction de u)
                    K_tilde.Get(mU1, m) += 1.0;
                    // partie 2 2 (eq de u' en fonction de u')
                    K_tilde.Get(mU1, mU1) -= 1.0 + alpha * R / 8.0 - R * rhoInfini * this->omega * this->omega / (8.0*alpha);
                    M_tilde.Get(mU1, mU1) += R / (8.0*alpha) * this->omega * this->omega;
                    // partie 2 3 (eq de u' en fonction de v)
                    K_tilde.Get(mU1, mV) += R * rhoInfini * this->omega * this->omega / 2.0;
                    M_tilde.Get(mU1, mV) += R / 2.0 * this->omega * this->omega;
                    // partie 2 4 (eq de u' en fonction de v)
                    K_tilde.Get(mU1, mW) += R * rhoInfini * this->omega * this->omega;
                    M_tilde.Get(mU1, mW) += R * this->omega * this->omega;

                    // partie 3 2 (eq de v en fonction de u')
                    K_tilde.Get(mV, mU1) += 1.0;
                    // partie 3 3 (eq de v en fonction de v)
                    K_tilde.Get(mV, mV) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                    // M_tilde.Get(mV, mV) += R / alpha * this->omega * this->omega;
                    M_tilde.Get(mV, mV) += 1.0 / alpha * this->omega * this->omega;

                    // partie 4 2 (eq de w en fonction de u')
                    // K_tilde.Get(mW, m) += 1.0;
                    K_tilde.Get(mW, mU1) += 1.0;
                    // partie 4 4 (eq de w en fonction de w)
                    K_tilde.Get(mW, mW) -= alpha / 2.0 - rhoInfini * this->omega * this->omega / (2.0*alpha);
                    M_tilde.Get(mW, mW) += 1.0 / (2.0*alpha) * this->omega * this->omega;
                    // partie 4 5 (eq de w en fonction de x)
                    K_tilde.Get(mW, mX) += 2.0 * rhoInfini * this->omega * this->omega;
                    M_tilde.Get(mW, mX) += 2.0 * this->omega * this->omega;

                    // partie 5 4 (eq de x en fonction de w)
                    K_tilde.Get(mX, mW) += 1.0;
                    // partie 5 5 (eq de x en fonction de x)
                    K_tilde.Get(mX, mX) -= alpha - rhoInfini * this->omega * this->omega / alpha;
                    M_tilde.Get(mX, mX) += 1.0 / alpha * this->omega * this->omega;
                  }
              }
          }
      }
    M_tilde.WriteText("Mt.dat");
    K_tilde.WriteText("Kt.dat");

    SparseEigenProblem<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>,
                       Matrix<Complex_wp, General, ArrayRowSparse> > var_eig;

    var_eig.SetStoppingCriterion(1e-12);
    var_eig.SetNbAskedEigenvalues(nb_asked_eigenval);
    var_eig.SetComputationalMode(var_eig.INVERT_MODE);

    var_eig.InitMatrix(K_tilde, M_tilde);
    var_eig.SetPrintLevel(2);

    var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, square(this->shift), var_eig.SORTED_MODULUS);

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
    for (int i = 0; i < lambda.GetM(); i++) {
        betas.Get(i) = sqrt(lambda.Get(i));
    }
    chrono.Stop("PVP");
    fin = clock();

    DISP(betas);

    DISP(type_coord_other_fields);
    chrono.DisplayTime("PVP");
    cout << "Temps de résolution : " << (fin-debut)/1000 << endl;
    cout << debut << endl;
    cout << fin << endl;


    bool premierModeDoubleFait = false;
    Complex_wp a1, a2, b1, b2;
    VectComplex_wp u1, u2;

    // VectComplex_wp U;
    for (int i = 0; i < betas.GetM(); i++)
      {
        int k = betas.GetM() - i - 1;
        VectComplex_wp vecU;
        GetCol(eigen_vec, k, vecU);
        vecU.Resize(nb_dof_u);
        VectComplex_wp Y = vecU;
        DISP(vecU.GetM());
#ifdef SELDON_WITH_SLEPC
        // vérification (norme de T(L) X)
        type_CLA = BGT2;
        ModeUNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, Kh_surf, Mh_surf,
                                            this->use_split_formulation);
        var_eig.MltOperator(betas(k), SeldonNoTrans, vecU, Y);
        cout << "VP " << i << ": " << betas(k) << " ; norme : " << Norm2(Y) / Norm2(vecU) << endl;
        type_CLA = BGT1;
        ModeUNonLinearEigenProblem var_eig2(N, *this, vec_Mh, vec_Kh, Kh_surf, Mh_surf,
                                            this->use_split_formulation);
        var_eig2.MltOperator(betas(k), SeldonNoTrans, vecU, Y);
        cout << "VP " << i << ": " << betas(k) << " ; norme : " << Norm2(Y) / Norm2(vecU) << endl;

        // ComputeIntensity(Ez, vec_Vk, k);
#endif
        if (orthogonalisation)
          {
            int nb_dof = N;
            VectComplex_wp U1(nb_dof); U1.Zero();
            VectComplex_wp U2(nb_dof); U2.Zero();
            VectComplex_wp coefsModes;
            if (this->refOrtho == -1)
              {
                ComputeCoefficientsDOmega(vecU, Mh_surf, vec_theta_Omega, i, coefsModes);
              }
            else
              {
                ComputeCoefficientsFullCable(vecU, vec_Mh.Get(this->refOrtho), vec_theta_Omega, i, coefsModes);
              }
            if (coefsModes.GetM() == 1)
              {
                for (int k = 0; k < nb_dof; k++)
                  {
                    U1(k) = vecU(k) / coefsModes.Get(0);
                  }

                WriteOutputFile(U1, i, lambda(k));
              }
            else if (coefsModes.GetM() == 2)
              {
                if (!premierModeDoubleFait)
                  {
                    a1 = coefsModes.Get(0);
                    a2 = coefsModes.Get(1);
                    u1 = vecU;
                    premierModeDoubleFait = true;
                    cout << "Orthogonalisation avec le prochain mode." << endl;
                    if (i + 1 == betas.GetM()) // dernier mode
                      {
                        cout << "Le prochain mode n'existe pas :(" << endl;
                        WriteOutputFile(vecU, i, betas(k));
                      }
                  }
                else
                  {
                    b1 = coefsModes.Get(0);
                    b2 = coefsModes.Get(1);
                    u2 = vecU;
                    premierModeDoubleFait = false;

                    Complex_wp Det = a1 * b2 - a2 * b1;

                    for (int k = 0; k < nb_dof; k++)
                      {
                        U1(k) = (b2 * u1(k) - a2 * u2(k)) / Det;
                        U2(k) = (-b1 * u1(k) + a1 * u2(k)) / Det;
                      }

                    WriteOutputFile(U1, i-1, lambda(k+1));
                    WriteOutputFile(U2, i, lambda(k));
                  }
              }
            else
              {
                cout << "Problème orthogonalisation..." << endl;
                WriteOutputFile(vecU, i, betas(k));
              }
          }
        else
          {
            WriteOutputFile(vecU, i, betas(k));
          }
      }


//     VectComplex_wp x;
//     for (int i = 0; i < betas.GetM(); i++)
//       {
//         GetCol(eigen_vec, i, x);
//         VectComplex_wp Ez;
//         VectComplex_wp Hz; // on ne prend que les données de Ez et de Hz
//         VectComplex_wp X = x; X.Resize(nb_dof_Ez + nb_dof_Hz);
//         VectComplex_wp Y = x; Y.Resize(nb_dof_Ez + nb_dof_Hz);
//         Ez.Reallocate(nb_dof_Hz); Ez.Zero();
//         Hz.Reallocate(nb_dof_Hz); Hz.Zero();
//
// #ifdef SELDON_WITH_SLEPC
//         // vérification (norme de T(L) X)
//         ModeUNonLinearEigenProblem var_eig(N, *this, vec_Mh, vec_Kh, vec_Ch,
//                                             this->DofKeptDir, this->IndexDirichlet);
//
//         var_eig.MltOperator(betas(i), SeldonNoTrans, X, Y);
//         cout << "VP " << i << ": " << betas(i) << " ; norme : " << Norm2(Y) / Norm2(X) << endl;
//
//         x.Reallocate(2*nb_dof_Hz);
//         var_eig.ExtractDirichlet(X, x);
//         for (int i = 0; i < nb_dof_Hz; i++)
//           {
//             Ez(i) = x(i);
//             Hz(i) = x(nb_dof_Hz + i);
//           }
// #endif
//       }

    return;
  }

};

#ifdef SELDON_WITH_SLEPC
// on inclut les fonctions definies dans le cxx
#include "mode_u.cxx"
#endif

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc < 2)
    {
      cout << "Data file not given" << endl;
      cout << "Usage : ./mode_u.x data_file.ini" << endl;
      abort();
    }

  // on appelle RunAll avec le fichier de donnees
  ModeU_Solver var;
  string input_file(argv[1]);
  var.RunAll(PATH_AUGUSTIN + to_str("/Datas/") + input_file);

  return FinalizeMontjoie();
}
