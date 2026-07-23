#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#define MONTJOIE_WITH_TIME_REVERSAL

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

using namespace Montjoie;

namespace Montjoie
{
  template<class TypeElement, class TypeEquation>
  class VarMigration : public TimeAcoustic_Base<TypeElement, TypeEquation>
  {
  public :
    typedef typename TypeElement::Dimension Dimension;
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
    
    int nb_shots; // nombre total de shots
    string file_root; // racine des fichiers ou sont stockees les infos des recepteurs
    Vector<R_N> xs; // xs(i) position de l'emetteur i
    Vector<Vector<R_N> > xr; // xr(i) : positions des recepteurs associes au shot i
    int shot_number; // numero du shot en cours
    IVect list_points_sismo; // numero du point ou on veut connaitre la solution
    Vector<VectReal_wp> valU_time; // valU_time(i) valeur de la solution aux recepteurs et a l'iteration i en temps
    IVect ref_bord; // reference des bords du domaine (d'ou le signal sera retourne)
    Vector<VectReal_wp> Vn_time; // dU_dn_time(i) valeur de du/dn sur les bords a l'iteration i en temps
    MeshInterpolation<TypeElement> var_int; // integration sur le bord
    Matrix<int> NodleBoundary; // numerotation des ddls du bord
    int type_source;
    enum { SRC_DIRECTE, SRC_REVERSE, SRC_RECEIVERS};
    Real_wp radius_gaussian; // rayon de distribution pour la source spatiale
    Vector<Real_wp, VectSparse> src_emit; // source pour l'emission
    Vector<Vector<Real_wp, VectSparse> > src_recv; // source des recepteurs
    IVect OffsetQuadratureNumber;
    
    // constructeur
    VarMigration() : TimeAcoustic_Base<TypeElement, TypeEquation>()
    {
      shot_number = 0;
      nb_shots = 121;
      file_root = string("/home/mab/durufle/Migration/rcv/rcv_data_");    
      type_source = SRC_DIRECTE;
      radius_gaussian = 50.0;
    }
    
    void SetInputData(const string& keyword, const Vector<string>& parameters)
    {
      TimeAcoustic_Base<TypeElement, TypeEquation>::SetInputData(keyword, parameters);
      
      if (!keyword.compare("ReferenceBoundaryReverse"))
	{
	  ref_bord.Reallocate(parameters.GetM());
	  for (int i = 0; ref_bord.GetM(); i++)
	    ref_bord(i) = to_num<int>(parameters(i));
	}
      else if (!keyword.compare("EmitterRadius"))
	{
	  radius_gaussian = to_num<Real_wp>(parameters(0));
	}
      else if (!keyword.compare("ShotData"))
	{
	  if (!parameters(0).compare("REGULAR"))
	    {
	    }
	  else
	    {
	    }
	}
    }
    
    void InitTimeIterations()
    {
      TimeAcoustic_Base<TypeElement, TypeEquation>::InitTimeIterations();
      
      // nombre iterations en temps
      int nb_iteration_time = toInteger(ceil((this->final_time - this->initial_time)/this->deltat)) + 1;
      if (type_source == SRC_DIRECTE)
	{
	  valU_time.Reallocate(nb_iteration_time);
	  Vn_time.Reallocate(nb_iteration_time);
	}

    }
    
    void InitQuadrature()
    {
      int nb_elt = this->var_harmonic.mesh.GetNbElt();
      OffsetQuadratureNumber.Reallocate(nb_elt+1);
      OffsetQuadratureNumber(0) = 0;
      for (int i = 0; i < nb_elt; i++)
	{
	  OffsetQuadratureNumber(i+1) = OffsetQuadratureNumber(i)
	    + this->var_harmonic.GetNbPointsQuadratureInside(i);
	}
    }
    
    void SetParamSource(const R_N& pt)
    {
      for (int k = 0; k < Dimension::dim_N; k++)
	this->var_harmonic.source_space_param(0)(k) = pt(k);
      
      this->var_harmonic.source_space_param(0)(Dimension::dim_N) = radius_gaussian;      
    }
    
    
    void GetSparseSource(Vector<Real_wp, VectSparse>& x)
    {
      VectReal_wp space_source;
      this->var_harmonic.ComputeRightHandSide(space_source);
      
      // conversion to a sparse vector
      Real_wp threshold = epsilon_machine*Norm2(space_source);
      int nb_dof_source = 0;
      for (int i = 0; i < space_source.GetM(); i++)
	if (abs(space_source(i)) > threshold)
	  nb_dof_source++;
      
      x.Reallocate(nb_dof_source);
      for (int i = 0; i < space_source.GetM(); i++)
	if (abs(space_source(i)) > threshold)
	  {
	    x.Index(nb_dof_source) = i;
	    x.Value(nb_dof_source) = space_source(i);
	    nb_dof_source++;
	  }
      
    }

    
    void ComputeRightHandSide()
    {
      this->var_harmonic.source_space_param.Reallocate(1);
      this->var_harmonic.source_space_param(0).Reallocate(Dimension::dim_N+1);
      this->var_harmonic.source_space_param(0).Fill(0);
      
      if (type_source == SRC_DIRECTE)
	{
	  // source sur l'emetteur shot_number
	  SetParamSource(xs(shot_number));
	  
	  // on calcule la source et on stocke dans un vecteur creux
	  GetSparseSource(src_emit);
	}
      else if (type_source == SRC_REVERSE)
	{
	  // nothing to do in that case
	}
      else if (type_source == SRC_RECEIVERS)
	{
	  int nb_receivers = xr(shot_number).GetM();
	  src_recv.Reallocate(nb_receivers);
	  for (int i = 0; i < nb_receivers; i++)
	    {
	      // source sur le i-eme receveur
	      SetParamSource(xr(shot_number)(i));
	      
	      // on calcule la source et on stocke dans un vecteur creux
	      GetSparseSource(src_recv(i));
	    }
	}
    }
    
    // on lit la liste des receveurs dans les fichiers .H@
    void ReadReceiver()
    {
      xs.Reallocate(nb_shots);
      xr.Reallocate(nb_shots);
      for (int i = 0; i < nb_shots; i++)
	{
	  string file_name = file_root + to_str(i+1) + ".H@";
	  Vector<double> y;
	  y.ReadText(file_name);
	  xs(i).Init(y(3), y(4));
	  xr(i).Reallocate(int(y(0)));
	  for (int j = 0; j < xr(i).GetM(); j++)
	    xr(i)(j).Init(y(5*j+1), y(5*j+2)); 
	}

    }
    
    // stockage de la solution
    void WriteInstantane(int nb_iter, const Real_wp& t, const VectReal_wp& Un)
    {
      TimeAcoustic_Base<TypeElement, TypeEquation>::WriteInstantane(nb_iter, t, Un);
      
      if (type_source == SRC_DIRECTE)
	{
	  Vector<VectReal_wp> xsol(1);
	  xsol(0) = Un; VectReal_wp val_u; TinyVector<int, 1> ghost_param;
	  this->var_harmonic.ComputeInterpolationU(xsol, this->var_harmonic.all_points_display, val_u,
						   list_points_sismo, list_points_sismo.GetM(), ghost_param);
	  
	  // on stocke la valeur de u sur tous les recepteurs
	  valU_time(nb_iter) = val_u;	  
	}
      
    }

    void GiveVectorialIterate(int nb_iter, const Real_wp& t, const VectReal_wp& Vh)
    {
      if (type_source == SRC_DIRECTE)
	GetVn(Vh, Vn_time(nb_iter));
    }
    
    // on change le numero du shot a traiter
    void SetShotNumber(int i)
    {
      shot_number = i;
      
      if (type_source == SRC_DIRECTE)
	{
	  // in this case, we need to know solution on receivers
	  VectReal_wp TetaPoints(xr(i).GetM()); IVect NumPartie(xr(i).GetM());
	  TetaPoints.Zero(); NumPartie.Zero();
	  this->var_harmonic.all_points_display.Append(xr(i), TetaPoints, NumPartie, list_points_sismo);
	  
	  // localisation des recepteurs dans le maillage
	  abort();
	  //this->var_harmonic.all_points_display.LocalizePoints(this->var_harmonic.mesh, &this->var_harmonic);
	  
	  // 
	  // this->var_harmonic.source_space_param(0)
	}
    }
    
    void GetInterpolateTime(const Real_wp& tn, Vector<VectReal_wp>& all_Vn, VectReal_wp& Vn)
    {
      abort();
    }
    
    void AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int nb_deriv, VectReal_wp& b_src)
    {
      if (type_source == SRC_DIRECTE)
	{
	  // une seule source sur l'emetteur
	  Real_wp pulse;
	  SourceOnlyTime(tn, this->var_harmonic.GetFrequency(), nb_deriv, pulse);
	  pulse *= alpha;
	  
	  for (int i = 0; i < src_emit.GetM(); i++)
	    b_src(src_emit.Index(i)) += pulse*src_emit.Value(i);
	}
      else if (type_source == SRC_REVERSE)
	{
	  // source sur le bord du domaine
	  VectReal_wp Vn;
	  GetInterpolateTime(tn, Vn_time, Vn);
	  AddSurfacicSource(Real_wp(1), Vn, b_src);
	}
      else if (type_source == SRC_RECEIVERS)
	{
	  int nb_receivers = xr(shot_number).GetM();
	  // source sur les recepteurs
	  VectReal_wp pulse;
	  GetInterpolateTime(tn, valU_time, pulse);
	  for (int n = 0; n < nb_receivers; n++)
	    {
	      for (int i = 0; i < src_recv(n).GetM(); i++)
		b_src(src_recv(n).Index(i)) += pulse(n)*src_recv(n).Value(i);
	    }
	}
    }
    
    
    void AddSurfacicSource(const Real_wp& alpha, const VectReal_wp& Vn, VectReal_wp& b)
    {
      //SetPoints<Dimension> pts;
      //SetMatrices<Dimension> mat;
      //int nb_edges = var_int.mesh.GetNbBoundary();
      //int num_edge, num_elem, nb_pts_face, nb_dof_face, num_loc;
      // loop on edges of \Gamma
      //for (int i = 0; i < nb_edges; i++)
	{
	  abort();
	}
      
    }
    
    //! computation of v dot n
    void GetVn(const VectReal_wp& Vn, VectReal_wp& v_dot_n) const
    {
      int N = var_int.GetNbPointsQuadrature();
      int nb_edges = var_int.mesh.GetNbBoundary();
      int num_edge, num_elem, nb_pts_face, nb_dof_face, num_loc;
      v_dot_n.Reallocate(N); v_dot_n.Fill(0);
      // loop on edges of \Gamma
      for (int i = 0; i < nb_edges; i++)
	{
	  /*
	  
	  // we get values of v
	  for (int j = 0; j < nb_pts_quad; j++)
	    {
	      int num_dof_loc = this->var_harmonic.GetQuadNumber(num_elem, num_loc, j);
	      int offset = (OffsetQuadratureNumber(i) + num_dof_loc)*TypeEquation::nb_unknowns;
	      R_N normale = var_int.mesh.NormaleQuadrature(i, j);
	      
	      // computation of scalar product v. n
	      Real_wp vloc = 0;
	      for (int k = 0; k < Dimension::dim_N; k++)
		vloc += normale(k)*Vn(offset+k);
	      
	      // we store the value in output vector
	      v_dot_n(NodleBoundary(i, j)) = vloc;
	    }
	  */
	}
    }
    
  };
  
  
  class AcousticEquationMigration : public AcousticEquation<Dimension2>
  {
  };
  
  template<class TypeElement>
  class HyperbolicProblem<TypeElement, AcousticEquationMigration> : public VarMigration<TypeElement, AcousticEquationMigration>
  {
  public :
    HyperbolicProblem() : VarMigration<TypeElement, AcousticEquationMigration>()
    {}
    
  };
  
  /*
    template<class TypeElement, class Vector1>
    void MltAdd_StiffnessVectorial(const Real_wp& alpha, int level, const HyperbolicProblem<TypeElement,
    AcousticEquationMigration>& var,
    const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C, Vector1& extrapolU)
    {
    MltAdd_StiffnessVectorial_AcousticDG(alpha, level, var, B, beta, C, extrapolU);
    }
    
    template<class TypeElement, class Vector1>
    void MltAdd_StiffnessScalar(const Real_wp& alpha, int level, const HyperbolicProblem<TypeElement,
    AcousticEquationMigration>& var,
    const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C, Vector1& extrapolV)
    {
    MltAdd_StiffnessScalar_AcousticDG(alpha, level, var, B, beta, C, extrapolV);
    }
  

    template<class TypeElement, class Vector1>
    void MltAdd_AbsorbingScalar(const Real_wp& alpha, const HyperbolicProblem<TypeElement, 
    AcousticEquationMigration>& var_time,
    const Vector1& extrapolU, VectReal_wp& C)
    {
    MltAdd_AbsorbingScalar_AcousticDG(alpha, var_time, extrapolU, C);
    }

    template<class TypeElement, class Vector1>
    void MltAdd_AbsorbingVectorial(const Real_wp& alpha, const HyperbolicProblem<TypeElement, 
    AcousticEquationMigration>& var_time,
    const Vector1& extrapolV, VectReal_wp& C)
    {
    MltAdd_AbsorbingVectorial_AcousticDG(alpha, var_time, extrapolV, C);
    }
  */
}


template<class TypeElement, class TypeEquation>
void AddBoundaryReferenceDirect(HyperbolicProblem<TypeElement, TypeEquation>& var_direct,
				HyperbolicProblem<TypeElement, TypeEquation>& var_reverse)
{
  // on remplit ref_cond pour que ref_cond == 1 corresponde aux bords du domaine
  IVect ref_cond(var_direct.var_harmonic.mesh.GetNbReferences()+1), inv_ref_bord(ref_cond.GetM());
  int ref = 1; ref_cond.Fill(-1); inv_ref_bord.Fill(-1);
  for (int i = 0; i < var_direct.ref_bord.GetM(); i++)
    {
      ref_cond(var_direct.ref_bord(i)) = ref;
      inv_ref_bord(var_direct.ref_bord(i)) = i;
    }
  
  // pour le maillage "direct", il faut remettre les bords intermediaires (puisqu'on a mis des PMLs)
  // et changer les references des PML, ce qu'on fait en premier lieu
  IVect new_ref(var_direct.ref_bord.GetM());
  for (int i = 0; i < new_ref.GetM(); i++)
    new_ref(i) = var_direct.var_harmonic.mesh.GetNewReference();
  
  for (int i = 0; i < var_direct.var_harmonic.mesh.GetNbBoundaryRef(); i++)
    {
      int ref_b = var_direct.var_harmonic.mesh.BoundaryRef(i).GetReference();
      if (ref_cond(ref_b) == ref)
	var_direct.var_harmonic.mesh.BoundaryRef(i).SetReference(new_ref(inv_ref_bord(ref_b)));
    }
  
  // on compte toutes les aretes/faces a rajouter
  int nb_old_edges = var_direct.var_harmonic.mesh.GetNbBoundaryRef();
  int nb_new_edges = nb_old_edges;
  for (int i = 0; i < var_reverse.var_harmonic.mesh.GetNbBoundaryRef(); i++)
    if (ref_cond(var_direct.var_harmonic.mesh.BoundaryRef(i).GetReference()) == ref)
      nb_new_edges++;
  
  // puis on les rajoute
  var_direct.var_harmonic.mesh.ResizeBoundariesRef(nb_new_edges);
  nb_new_edges = nb_old_edges;
  for (int i = 0; i < var_reverse.var_harmonic.mesh.GetNbBoundaryRef(); i++)
    {
      int ref_b = var_reverse.var_harmonic.mesh.BoundaryRef(i).GetReference();
      if (ref_cond(ref_b) == ref)
	{
	  var_direct.var_harmonic.mesh.BoundaryRef(nb_new_edges) = 
	    var_reverse.var_harmonic.mesh.BoundaryRef(i);
	  
	  nb_new_edges++;
	}
    }
  
  // on refait la correspondance faces -> faces de reference
  var_direct.var_harmonic.mesh.ConstructCrossReferenceBoundaryRef();
  
  // maintenant on peut extraire le maillage surfacique
  //var_reverse.var_int.InitIntegrationSurfacicMesh(ref, ref_cond);
  //var_direct.var_int.InitIntegrationSurfacicMesh(ref, ref_cond);
  abort();
  
  // filling NodleBoundary
  /*nb_new_edges = var_direct.var_int.mesh.GetNbBoundaryRef();
  var_direct.NodleBoundary.Reallocate(nb_new_edges, var_direct.var_harmonic.GetNbMaxQuadBoundary());
  var_direct.NodleBoundary.Fill(-1); int nb = 0;
  for (int i = 0; i < nb_new_edges; i++)
    {
      int num_edge_ref = var_direct.var_int.mesh.ListeBoundaries(i);
      int num_elem = var_reverse.var_harmonic.mesh.BoundaryRef(num_edge_ref).numElement(0);
      int num_edge = var_reverse.var_harmonic.mesh.GetBoundary_FromBoundaryRef(num_edge_ref);
      int num_loc = var_reverse.var_harmonic.mesh.FindLocalBoundary(num_elem, num_edge);
      int nb_pts_quad = var_reverse.var_harmonic.GetNbPointsQuadBoundary(num_elem, num_loc);
      for (int j = 0; j < nb_pts_quad; j++)
	var_direct.NodleBoundary(i, j) = nb++;
    }
  
  var_direct.InitQuadrature();
  var_reverse.InitQuadrature();*/
}


template<class TypeElement, class TypeEquation>
void InitializeComputation(HyperbolicProblem<TypeElement, TypeEquation>& var_direct,
			   HyperbolicProblem<TypeElement, TypeEquation>& var_reverse, const string& input_file)
{
  typedef typename TypeEquation::TypeEquationStationary TypeEquation_Harmonic;
  var_direct.var_harmonic.InitIndices(50);
  
  // on lit le fichier de donnees
  ReadInputFile(input_file, var_direct);
  ReadInputFile(input_file, var_direct.glob_solver);
  ReadInputFile(input_file, var_reverse);
  ReadInputFile(input_file, var_reverse.glob_solver);
  
  // for the reverse simulation, we remove PML
  // var_reverse.var_harmonic.ResetValuesPML();
  
  // on lit le maillage et on numerote
  var_direct.var_harmonic.ComputeMeshAndFiniteElement();
  var_reverse.var_harmonic.ComputeMeshAndFiniteElement();
  AddBoundaryReferenceDirect(var_direct, var_reverse);
  // on calcule les donnees geometriques et indices variables
  var_direct.var_harmonic.ComputeMassMatrix();
  var_reverse.var_harmonic.ComputeMassMatrix();

  //  on calcule la masse et la rigidite
  typename TypeEquation_Harmonic::Nature_Matrix nat_mat;
  nat_mat.SetCoefMass(0.0); nat_mat.SetCoefStiffness(1.0);
  var_direct.Glob_mat_Kh.Clear(); var_reverse.Glob_mat_Kh.Clear();
  var_direct.var_harmonic.AddMatrixFEM(var_direct.Glob_mat_Kh, nat_mat);
  var_reverse.var_harmonic.AddMatrixFEM(var_reverse.Glob_mat_Kh, nat_mat);
    
  // condition Dirichlet et autres
  var_direct.var_harmonic.TreatDirichletCondition();
  var_reverse.var_harmonic.TreatDirichletCondition();

  var_direct.ComputeMassMatrix();
  var_reverse.ComputeMassMatrix();
  
  var_direct.InitTimeIterations();
  var_reverse.InitTimeIterations();
}


template<class TypeElement, class TypeEquation, class TypeScheme>
void GetInitialConditionReverse(const HyperbolicProblem<TypeElement, TypeEquation>& var_direct,
				const HyperbolicProblem<TypeElement, TypeEquation>& var_reverse,
				TypeScheme& scheme, VectReal_wp& U0, VectReal_wp& V0, VectReal_wp& V0_pml)
{
  typedef typename TypeElement::Dimension Dimension;
  const Mesh<Dimension>& meshr = var_reverse.var_harmonic.mesh;
  const Mesh<Dimension>& meshd = var_direct.var_harmonic.mesh;
  GhostIf<TypeEquation::TypeEquationStationary::FormulationDG> dg_form;
  
  // allocation des tableaux et mise a zero
  U0.Reallocate(meshr.GetNbDof()); U0.Fill(0);
  V0.Reallocate(var_reverse.GetNbVectorialUnknowns(dg_form)); V0.Fill(0);
  V0_pml.Reallocate(var_reverse.GetNbVectorialUnknowns_PML(dg_form)); V0_pml.Fill(0);

  // on boucle sur tous les elements du maillage 'reverse'
  int offset = 0;
  for (int i = 0; i < meshr.GetNbElt(); i++)
    {
      int nb_dof = var_direct.var_harmonic.GetNbLocalDof(i);
      int nb_quad = var_direct.var_harmonic.GetNbPointsQuadratureInside(i);
      for (int j = 0; j < nb_dof; j++)
	U0(meshr.GetNumberDof(i, j)) = scheme.Un(meshr.GetNumberDof(i, j));
      
      for (int j = 0; j < nb_quad; j++)
	{
	  for (int k = 0; k < Dimension::dim_N; k++)
	    {
	      int num = offset + Dimension::dim_N*j + k;
	      V0(num) = scheme.Vn(num);
	    }
	}
      
      offset += nb_quad*Dimension::dim_N;
    }
  
  scheme.Clear();
}


template<class TypeElement, class TypeEquation>
void AddConvolution(const Real_wp& alpha, 
		    const EllipticProblem<TypeElement, TypeEquation>& var_d,
		    const EllipticProblem<TypeElement, TypeEquation>& var_r,
		    const VectReal_wp& Un_dir, const VectReal_wp& Un_rev, VectReal_wp& image_rho)
{
  Vector<bool> DofUsed(var_d.mesh_num.GetNbDof()); DofUsed.Fill(false);
  for (int i = 0; i < var_r.mesh.GetNbElt(); i++)
    {
      int nb_quad = var_r.mesh_num.GetNbLocalDof(i);
      for (int j = 0; j < nb_quad; j++)
	{
	  int num_dof = var_d.mesh_num.Element(i).GetNumberDof(j);
	  int num_dofr = var_r.mesh_num.Element(i).GetNumberDof(j);
	  if (!DofUsed(num_dof))
	    {
	      image_rho(num_dof) += alpha*Un_rev(num_dofr)*Un_dir(num_dof);
	      DofUsed(num_dof) = true;
	    }
	}
    }
}


template<class TypeElement, class TypeEquation>
void RunAll(HyperbolicProblem<TypeElement, TypeEquation>& var_direct, const string& input_file)
{
  typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
  HyperbolicProblem<TypeElement, TypeEquation> var_reverse;
  InitializeComputation(var_direct, var_reverse, input_file);
  
  var_direct.ReadReceiver();
  VectReal_wp image_rho(var_direct.var_harmonic.mesh_num.GetNbDof());
  image_rho.Fill(0);
  int Ns = var_direct.GetNbScalarUnknowns();
  int Nv = var_direct.GetNbVectorialUnknowns();
  for (int i = 0; var_direct.nb_shots; i++)
    {
      var_direct.SetShotNumber(i);
      
      // calcul de la source en espace
      var_direct.ComputeRightHandSide();
      
      // on lance le calcul direct 
      // nombre d'iterations en temps
      Real_wp t0 = var_direct.initial_time;
      Real_wp tf = var_direct.final_time;
      Real_wp dt = var_direct.deltat, t = t0;
      int nb_max_iter = toInteger(ceil(abs(tf - t0)/dt));
      
      ModifiedEquationSystemIterator<Real_wp> direct_scheme;      
      VectReal_wp U0(Ns), V0(Nv);
      U0.Fill(0); V0.Fill(0);
      direct_scheme.SetInitialCondition(t0, dt, U0, V0, var_direct);
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
	{
	  t = t0 + (n_time+1)*dt;
	  direct_scheme.Advance(t, n_time, var_direct);
	}      
      
      // on initialise les deux schemas (receveurs et reverse)
      ModifiedEquationSystemIterator<Real_wp> recv_scheme, reverse_scheme;
      //GetInitialConditionReverse(var_direct, var_reverse, direct_scheme, U0, V0);
      reverse_scheme.SetInitialCondition(t0, dt, U0, V0, var_reverse);
      
      U0.Reallocate(Ns); V0.Reallocate(Nv);
      U0.Fill(0); V0.Fill(0);
      recv_scheme.SetInitialCondition(t0, dt, U0, V0, var_direct);
      
      // on boucle sur les schemas
      for (int n_time = 0; n_time < nb_max_iter; n_time++)
	{
	  t = t0 + (n_time+1)*dt;
	  // on avance les deux schemas
	  recv_scheme.Advance(t, n_time, var_reverse);
	  reverse_scheme.Advance(t, n_time, var_direct);
	  
	  // we add \int u(t, x) u'(t, x) dt
	  //AddConvolution(dt, var_direct.var_harmonic, var_reverse.var_harmonic,
	  //recv_scheme.Un, reverse_scheme.Un, image_rho);
	  
	}
    }      
  
  var_direct.var_harmonic.WriteOutputFile(image_rho, "ImageRho.dat");
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      HyperbolicProblem<AcousticEquationMigration> Vars;
      // HyperbolicProblem<Triangle_DG_Classical, AcousticEquationMigration> Vars;

      RunAll(Vars, file_name_data);

      cout<<" we destroy the variables "<<endl;      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"time2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
