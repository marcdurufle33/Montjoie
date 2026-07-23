#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG

#define MONTJOIE_WITH_TIME_REVERSAL

//#define MONTJOIE_WITH_TRANSMISSION

#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"
#include "Instationary/EvaluationCflWave.cxx"

using namespace Montjoie;

//bool inhg_case = false;

//! function f(t, y) = M^-1 K y
class VectorialFunctionFirst : public VirtualOdeSystem<Real_wp>
{
private:
  const Matrix<Real_wp, General, ArrayRowSparse>& K;
  int Nproche;
  Real_wp tn;
  All_MatrixLU<Real_wp> facto_real;
  Vector<All_MatrixLU<Complex_wp> > facto_cplx;
  const Vector<VectReal_wp>* wp_i_ptr;
  const Vector<VectReal_wp>* Fp_i_ptr;
  VectReal_wp F;
  Globatto<Real_wp> lob;
  bool presence_source;
  Real_wp dt;
  
public :
  VectorialFunctionFirst(Matrix<Real_wp, General, ArrayRowSparse>& Kh_sub, int N,
			 const Real_wp& deltat, const VectReal_wp& ci)
    : K(Kh_sub)
  {
    Nproche = N;
    lob.AffectPoints(ci);
    dt = deltat;
  }
  
  Real_wp GetSourceTime(const Real_wp& t, int n)
  {
    return exp(-t/4)*pow(-Real_wp(1)/4, n);
  }
  
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    Mlt(K, y, f);
    if(source)
      AddPrimitiveTimeSource(Real_wp(1), t, 0, f);
  }

  // f = beta f + alpha Kh y
  void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& y,
		       const Real_wp& beta, Vector<Real_wp>& f)
  {
    MltAdd(alpha, K, y, beta, f);
  }

  void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const Vector<Real_wp>& y,
		       const Real_wp& beta, Vector<Real_wp>& f)
  {
    f = beta*f + alpha*y;
  }

  void SolveOperatorDh(Vector<Real_wp>& f)
  {
  }

  void GiveNumberIterations(const Real_wp& dt, int N) {}
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha = 1.0) {}

  void FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c) 
  {
    Matrix<Real_wp, General, ArrayRowSparse> A(Nproche, Nproche);
    A.SetIdentity();
    Mlt(a, A);
    Add(-b, K, A);

    facto_real.Factorize(A);
  }

  Real_wp GetFinalTimeSource() const {return 1e300;}

  void Assemble(Vector<Real_wp>& Vh) const{}

  void SolveMassMatrix(Vector<Real_wp>& X){ }
  
  void ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t, 
		       const Vector<Real_wp>& x, 
		       const Real_wp& beta, 
		       Vector<Real_wp>& y)
  {
    y = beta*y + alpha*x;
  }

  // solving y - gamma f(t, y) = f
  void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& f, VectReal_wp& y)
  {
    y = f;
    facto_real.Solve(y);
  }


  void FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
				const VectComplex_wp& gamma)    
  {
    facto_cplx.Reallocate(alpha.GetM());
    for (int k = 0; k < alpha.GetM(); k++)
      {
	Matrix<Complex_wp, General, ArrayRowSparse> A(Nproche, Nproche);
	A.SetIdentity();
	
	Mlt(alpha(k), A);
	Add(-beta(k), K, A);
	facto_cplx(k).Factorize(A);
      }
  }
  
  void SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y,
			    int num_system)
  {
    Y = X;
    facto_cplx(num_system).Solve(Y);
  }

  void SetSource(const Real_wp& t, const Vector<VectReal_wp>& wp_i,
		 bool presence_source_, const Vector<VectReal_wp>& Fp_i)
  {
    tn = t;
    wp_i_ptr = &wp_i;
    Fp_i_ptr = &Fp_i;
    presence_source = presence_source_;
  }
    
  void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			      int nb_deriv, Vector<Real_wp>& b_src)
  {
    Real_wp tau = t - tn;
    const Vector<VectReal_wp>& wp_i = *wp_i_ptr;
    const Vector<VectReal_wp>& Fp_i = *Fp_i_ptr;

    F = wp_i(wp_i.GetM()-1);
    for (int j = wp_i.GetM()-2; j >= 0; j--)
      F = tau*F + wp_i(j);

    if (presence_source)
      {
	Real_wp xi = tau / dt;
	VectReal_wp phi;
	lob.ComputeValuesPhiRef(xi, phi);
	for (int i = 0; i < phi.GetM(); i++)
	  F += Fp_i(i)*phi(i);
      }
    
    b_src = b_src + alpha*F;
  }
  
};



namespace Montjoie
{

  class InputDataTimeScheme : public InputDataProblem_Base
  {
  public :
    bool local_implicit;
    int order_imp, s_imp;
    bool pade_scheme;
    int order_exp, s_exp;

    enum {FILE_SPLITTING, AUTO_SPLITTING};
    int type_splitting;
    Real_wp dt_splitting;
    string file_with_number_element;
    
    int p; // number of local time step 

    InputDataTimeScheme()
    {
      local_implicit = false;
      order_imp = 4; s_imp = 0;
      pade_scheme = true;
      order_exp = 4; s_exp = 0;
      type_splitting = FILE_SPLITTING;
      file_with_number_element = "num_implicit.dat";
      p = 2;
    }
    
    void SetInputData(const string& keyword, const Vector<string>& param)
    {
      if (keyword == "LocalImplicit")
	{
	  if (param(0) == "YES")
	    {
	      local_implicit = true;
	      if (param(1) == "Pade")
		pade_scheme = true;
	      else
		pade_scheme = false;

	      order_imp = to_num<int>(param(2));
	      s_imp = to_num<int>(param(3));
	    }
	  else
	    {
	      local_implicit = false;
	      p = to_num<int>(param(1));
	    }
	}
      else if (keyword == "ExplicitScheme")
	{
	  order_exp = to_num<int>(param(0));
	  s_exp = to_num<int>(param(1));
	}
      else if (keyword == "SplittingCoarseFineRegion")
	{
	  if (param(0) == "AUTO")
	    {
	      // case where we compute local time steps on each element
	      type_splitting = AUTO_SPLITTING;
	      dt_splitting = to_num<Real_wp>(param(1));
	    }
	  else if (param(0) == "FILE")
	    {
	      type_splitting = FILE_SPLITTING;
	      file_with_number_element = param(1);
	    }
	  else
	    {
	      cout << "Splitting not implemented" << endl;
	      abort();
	    }
	}
    }
    
  };      

  void ExtractMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>& K,
		     int nb_dof_L, int nb_dof_uv,
		     VectReal_wp& C, DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
		     DistributedMatrix<Real_wp, General, ArrayRowSparse>& Al,
		     DistributedMatrix<Real_wp, General, ArrayRowSparse>& B)
  {
    A.Reallocate(nb_dof_uv, nb_dof_uv);
    Al.Reallocate(nb_dof_uv, nb_dof_L);
    for (int i = 0; i < nb_dof_uv; i++)
      {
	int size_row = K.GetRowSize(nb_dof_L + i);
	int nb_val = 0, nb_val2 = 0;
	for (int j = 0; j < size_row; j++)
	  {
	    if (K.Index(nb_dof_L+i, j) < nb_dof_L)
	      nb_val++;
	    else
	      nb_val2++;
	  }
	
	Al.ReallocateRow(i, nb_val);
	A.ReallocateRow(i, nb_val2);
	nb_val = 0; nb_val2 = 0;
	for (int j = 0; j < size_row; j++)
	  {
	    if (K.Index(nb_dof_L+i, j) < nb_dof_L)
	      {
		Al.Index(i, nb_val) = K.Index(nb_dof_L + i, j);
		Al.Value(i, nb_val) = K.Value(nb_dof_L + i, j);
		nb_val++;
	      }
	    else
	      {
		A.Index(i, nb_val2) = K.Index(nb_dof_L + i, j) - nb_dof_L;
		A.Value(i, nb_val2) = K.Value(nb_dof_L + i, j);
		nb_val2++;
	      }
	      }
	
	K.ClearRow(nb_dof_L + i);
      }
    
    B.Reallocate(nb_dof_L, nb_dof_uv);
    C.Reallocate(nb_dof_L);
    C.Zero();
    for (int i = 0; i < nb_dof_L; i++)
      {
	int size_row = K.GetRowSize(i);
	int nb_val = 0;
	for (int j = 0; j < size_row; j++)
	  if (K.Index(i, j) >= nb_dof_L)
	    nb_val++;
	
	B.ReallocateRow(i, nb_val);
	nb_val = 0;
	for (int j = 0; j < size_row; j++)
	  {
	    if (K.Index(i, j) >= nb_dof_L)
	      {
		B.Index(i, nb_val) = K.Index(i, j) - nb_dof_L;
		B.Value(i, nb_val) = K.Value(i, j);
		nb_val++;
	      }
	    else
	      {
		if (K.Index(i, j) != i)
		  {
		    cout << "non-diagonal case " << endl;
		    abort();
		  }
		else
		  C(i) = K.Value(i, j);
	      }
	  }
	K.ClearRow(i);
      }
    
    K.Clear();
  }
  

  void ExtractSubMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
			const IVect& row_num, const IVect& index_row,
			const IVect& col_num, const IVect& index_col,
			DistributedMatrix<Real_wp, General, ArrayRowSparse>& Asub)
  {
    Asub.Clear();
    if ((row_num.GetM() <= 0) || (col_num.GetM() <= 0))
      return;
    
    Asub.Reallocate(row_num.GetM(), col_num.GetM());
    for (int i2 = 0; i2 < row_num.GetM(); i2++)
      {
	int i = row_num(i2);
	int size_row = 0;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index_col(A.Index(i, j)) >= 0)
	    size_row++;

	Asub.ReallocateRow(i2, size_row);
	size_row = 0;
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (index_col(A.Index(i, j)) >= 0)
	    {
	      Asub.Index(i2, size_row) = index_col(A.Index(i, j));
	      Asub.Value(i2, size_row) = A.Value(i, j);
	      size_row++;
	    }
      }
  }
  
  
  // on lance la simulation generale
  // var_time : donnees du probleme a simuler
  // input_file : fichier de donnees
  template<class TypeEquationTime>
  void RunAll(HyperbolicProblem<TypeEquationTime>& var_time, const string& input_file,
	      const string& name_element, const string& name_equation)
  {
    // on recupere le probleme stationnaire
    typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
    typedef typename TypeEquation::Dimension Dimension;
    EllipticProblem<TypeEquation>& var = var_time.var_harmonic;

    var.SetTypeEquation(name_equation);
    
    // on lit le fichier de donnees
    Vector<string> lines_data_file;
    ReadLinesFile(input_file, lines_data_file);

    InputDataTimeScheme data;
    
    var.InitIndices(10);
    ReadInputFile(lines_data_file, var_time);
    ReadInputFile(lines_data_file, data);
    
    // solveur utilise pour le theta schema
    //All_LinearSolver glob_solver(var);
    //ReadInputFile(input_file, glob_solver);
    
    // construction du maillage et element fini
    var.SetFirstOrderFormulation(true);
    var.ComputeMeshAndFiniteElement(name_element);
    var.mesh.Write("test.mesh");
    
    // calcul du jacobien, DF_i^{-1} 
    var.ComputeMassMatrix();
    
    // on repere noeuds de Dirichlet
    var.TreatDirichletCondition();
    var_time.ComputeRightHandSide();
    
    // autres initialisations pour le temporel
    var_time.Glob_mat_Dh = var_time.GetNewMassMatrix();
    var_time.Glob_mat_Dh->compute_Dh_default = true;
    var_time.Glob_mat_Dh->compute_Bh_default = true;
    var_time.ComputeMassMatrix();
    var_time.InitTimeIterations();
    
    // elements avec pas de temps implicite
    Vector<int> num_implicit;
    DISP(data.type_splitting);
    if (data.type_splitting == data.AUTO_SPLITTING)
      {
	VectReal_wp dt_elt;
	ComputeLocalTimeStep(var_time, dt_elt);
	int ne = 0;
	for (int i = 0; i < dt_elt.GetM(); i++)
	  if (dt_elt(i) < data.dt_splitting)
	    ne++;

	num_implicit.Reallocate(ne);
	ne = 0;
	for (int i = 0; i < dt_elt.GetM(); i++)
	  if (dt_elt(i) < data.dt_splitting)
	    num_implicit(ne++) = i;
	
      }
    else if (data.type_splitting == data.FILE_SPLITTING)
      num_implicit.ReadText(data.file_with_number_element);


    {
      num_implicit.WriteText("num.dat");
      Mesh<Dimension> mesh_order(var.mesh);

      for (int i = 0; i < mesh_order.GetNbElt(); i++)
	mesh_order.Element(i).SetReference(1);

      for (int i = 0; i < num_implicit.GetM(); i++)
	mesh_order.Element(num_implicit(i)).SetReference(2);

      mesh_order.Write("order.mesh");
    }
    
    // on recupere la diagonale D_h
    const VectReal_wp& Dh = var_time.Glob_mat_Dh->GetDiagonalDh();
    const VectReal_wp& Bh = var_time.Glob_mat_Dh->GetDiagonalBh();
    //DISP(Dh.GetM()); DISP(Bh.GetM());
    //Dh.WriteText("Dh.dat");
    //Bh.WriteText("Bh.dat");

    VectReal_wp invMass(Dh.GetM() + Bh.GetM());
    for (int i = 0; i < Dh.GetM(); i++)
      invMass(i) = Real_wp(1) / Dh(i);

    for (int i = 0; i < Bh.GetM(); i++)
      invMass(Dh.GetM() + i) = Real_wp(1) / Bh(i);

    invMass.WriteText("invMass_ref.dat");
    
    // calcul de la matrice de rigidite
    GlobalGenericMatrix<Real_wp> nat_mat;
    DistributedMatrix<Real_wp, General, ArrayRowSparse> Kh, A22, A21, A12;
    VectReal_wp A11_diag;
    nat_mat.SetCoefMass(0.0);
    nat_mat.SetCoefStiffness(1.0);
    
    var.AddMatrixWithBC(Kh, nat_mat);
    //Kh.WriteText("Kh.dat");
	
    int nb_dof_L = var.mesh_num.GetNbDof();
    int nb_dof_uv = var.GetNbDof() - nb_dof_L;
    ExtractMatrix(Kh, nb_dof_L, nb_dof_uv, A11_diag, A22, A21, A12);
    
    // on enleve DF_i^{-1}
    var.Glob_DFjm1.Clear();
    
    //A11_diag.WriteText("diag.dat");
    //A12.WriteText("A12.dat");
    //A21.WriteText("A21.dat");
    //A22.WriteText("A22.dat");
    
    // ddls implicites
    Real_wp dt = var_time.GetTimeStep();
    Vector<int> dof_implicit, ddl_vol_imp;
    int Nvol = var.GetNbMainUnknownDof();
    for (int i2 = 0; i2 < num_implicit.GetM(); i2++)
      {
	int i = num_implicit(i2);
	dof_implicit.PushBack(var.mesh_num.Element(i).GetNodle());
	
	int v1 = var.GetOffsetDofV(i);
	int v2 = var.GetOffsetDofV(i+1);
	int nb_dof_v = v2 - v1;
	Vector<int> dof_v((Dimension::dim_N+1)*nb_dof_v);
	for (int k = 0; k <= Dimension::dim_N; k++)
	  for (int i = v1; i < v2; i++)
	    dof_v(i - v1 + k*nb_dof_v) = k*Nvol + i;
	
	ddl_vol_imp.PushBack(dof_v);
      }
    
    RemoveDuplicate(dof_implicit);
    Sort(ddl_vol_imp);
    
    dof_implicit.WriteText("dof_implicit_ref.dat");
    ddl_vol_imp.WriteText("ddl_implicit_ref.dat");
    
    int Nimp = ddl_vol_imp.GetM();
    
    // IndexDof(i) = numero local du ddl global i parmi les ddls implicites
    Vector<int> IndexDof(nb_dof_uv);
    IndexDof.Fill(-1);
    for (int i = 0; i < ddl_vol_imp.GetM(); i++)
      IndexDof(ddl_vol_imp(i)) = i;

    // calcul des ddls lointains
    int Nlointain = nb_dof_uv - Nimp;
    Vector<int> IndexLointain(nb_dof_uv), dof_lointain;
    IndexLointain.Fill(-1);
    dof_lointain.Reallocate(Nlointain);
    Nlointain = 0;
    for (int i = 0; i < nb_dof_uv; i++)
      if (IndexDof(i) == -1)
	{
	  dof_lointain(Nlointain) = i;
	  IndexLointain(i) = Nlointain;
	  Nlointain++;
	}

    dof_lointain.WriteText("dof_lointain_ref.dat");
    
    // calcul des ddls proches
    Vector<int> IndexProche(nb_dof_uv);
    IndexProche.Fill(-1);
    VectReal_wp vec_aleatoire(nb_dof_uv), prod_aleatoire(nb_dof_L);
    vec_aleatoire.Zero();
    for (int i = 0; i < ddl_vol_imp.GetM(); i++)
      {
	int ig = ddl_vol_imp(i);
	IndexProche(ig) = 1;
	vec_aleatoire(ig) = Real_wp(rand()) / RAND_MAX;
      }
    
    Mlt(A12, vec_aleatoire, prod_aleatoire);
    MltAdd(Real_wp(1), A21, prod_aleatoire, Real_wp(1), vec_aleatoire); 
    
    for (int i = 0; i < vec_aleatoire.GetM(); i++)
      if (vec_aleatoire(i) != Real_wp(0))
	IndexProche(i) = 1;
    
    int Nproche = 0;
    for (int i = 0; i < nb_dof_uv; i++)
      if (IndexProche(i) == 1)
	Nproche++;

    IVect dof_proche(Nproche);
    Nproche = 0;
    for (int i = 0; i < nb_dof_uv; i++)
      if (IndexProche(i) == 1)
	{
	  dof_proche(Nproche) = i;
	  IndexProche(i) = Nproche;
	  Nproche++;
	}
    
    dof_proche.WriteText("dof_proche_ref.dat");

    // ddls pour Lambda proches
    int nb_dof_L_proche = 0;
    for (int i = 0; i < prod_aleatoire.GetM(); i++)
      if (prod_aleatoire(i) != Real_wp(0))
	nb_dof_L_proche++;

    IVect IndexLambda(nb_dof_L);
    IndexLambda.Fill(-1);
    IVect dofL_proche(nb_dof_L_proche);
    nb_dof_L_proche = 0;
    for (int i = 0; i < nb_dof_L; i++)
      if (prod_aleatoire(i) != Real_wp(0))
	{
	  IndexLambda(i) = nb_dof_L_proche;
	  dofL_proche(nb_dof_L_proche) = i;
	  nb_dof_L_proche++;
	}
    
    dofL_proche.WriteText("dofL_proche_ref.dat");

    // ddls pour Lambda lointains
    int nb_dof_L_lointain = 0;
    vec_aleatoire.Zero();
    for (int i = 0; i < dof_lointain.GetM(); i++)
      {
	int ig = dof_lointain(i);
	vec_aleatoire(ig) = Real_wp(rand()) / RAND_MAX;
      }
    
    Mlt(A12, vec_aleatoire, prod_aleatoire);
    
    for (int i = 0; i < prod_aleatoire.GetM(); i++)
      if (prod_aleatoire(i) != Real_wp(0))
	nb_dof_L_lointain++;

    IVect IndexLointainL(nb_dof_L);
    IndexLointainL.Fill(-1);
    IVect dofL_lointain(nb_dof_L_lointain);
    nb_dof_L_lointain = 0;
    for (int i = 0; i < prod_aleatoire.GetM(); i++)
      if (prod_aleatoire(i) != Real_wp(0))
	{
	  IndexLointainL(i) = nb_dof_L_lointain;
	  dofL_lointain(nb_dof_L_lointain) = i;
	  nb_dof_L_lointain++;
	}

    dofL_lointain.WriteText("dofL_lointain_ref.dat");

    // on extrait les matrices pour la partie proche
    DistributedMatrix<Real_wp, General, ArrayRowSparse> A12_sub, A21_sub, A22_sub;
    VectReal_wp A11_sub;
    A11_sub.Reallocate(nb_dof_L_proche);
    for (int i = 0; i < nb_dof_L_proche; i++)
      A11_sub(i) = A11_diag(dofL_proche(i));
    
    ExtractSubMatrix(A12, dofL_proche, IndexLambda, ddl_vol_imp, IndexDof, A12_sub);
    ExtractSubMatrix(A21, dof_proche, IndexProche, dofL_proche, IndexLambda, A21_sub);
    ExtractSubMatrix(A22, dof_proche, IndexProche, ddl_vol_imp, IndexDof, A22_sub);

    A11_sub.WriteText("A11_sub_ref.dat");
    A12_sub.WriteText("A12_sub_ref.dat");
    A21_sub.WriteText("A21_sub_ref.dat");
    A22_sub.WriteText("A22_sub_ref.dat");

    // on extrait les matrices pour la partie lointaine
    DistributedMatrix<Real_wp, General, ArrayRowSparse> A12_far, A21_far, A22_far;
    VectReal_wp A11_far;
    
    A11_far.Reallocate(nb_dof_L_lointain);
    for (int i = 0; i < nb_dof_L_lointain; i++)
      A11_far(i) = A11_diag(dofL_lointain(i));

    {
      IVect all_dof(nb_dof_uv); all_dof.Fill();
      IVect IndexAll(nb_dof_uv); IndexAll.Fill();
      
      ExtractSubMatrix(A12, dofL_lointain, IndexLointainL, dof_lointain, IndexLointain, A12_far);
      ExtractSubMatrix(A21, all_dof, IndexAll, dofL_lointain, IndexLointainL, A21_far);
      ExtractSubMatrix(A22, all_dof, IndexAll, dof_lointain, IndexLointain, A22_far);
    }

    A11_far.WriteText("A11_far_ref.dat");
    A12_far.WriteText("A12_far_ref.dat");
    A21_far.WriteText("A21_far_ref.dat");
    A22_far.WriteText("A22_far_ref.dat");

    // schema explicite : Runge-Kutta optimise
    LinearRungeKutta_Iterator<Real_wp> explicit_scheme;
    explicit_scheme.SetOrder(data.order_exp, data.s_exp, false);
    const UnivariatePolynomial<Real_wp>& P = explicit_scheme.GetStabilityFunction();
    int orderN = P.GetOrder();
    VectReal_wp points_ci = explicit_scheme.GetPoints();
    const Matrix<Real_wp>& bi = explicit_scheme.GetWeights();
    
    // variables pour le schema temporel
    int nb_iter = toInteger(ceil((var_time.GetFinalTime() - var_time.GetInitialTime())/dt));
    VectReal_wp F_next(nb_dof_uv), yn(nb_dof_uv), y_proche(Nproche), y_imp(Nimp), y_far(Nlointain);
    VectReal_wp w(nb_dof_uv), wp(Nproche), w_tild(nb_dof_uv), wp_tild(Nimp);
    Vector<VectReal_wp> w_n_i(orderN), EvalF(points_ci.GetM());
    for (int k = 0; k < orderN; k++)
      {
	w_n_i(k).Reallocate(nb_dof_uv);
	w_n_i(k).Zero();
      }
    
    for (int k = 0; k < points_ci.GetM(); k++)
      EvalF(k).Reallocate(nb_dof_uv);
    
    // vecteurs pour le schema zone fine    
    Vector<VectReal_wp> wp_i(orderN), EvalFproche(points_ci.GetM());
    for (int k = 0; k < orderN; k++)
      {
	wp_i(k).Reallocate(Nproche);
	wp_i(k).Zero();
      }
    
    for (int k = 0; k < points_ci.GetM(); k++)
      EvalFproche(k).Reallocate(Nproche);

    // condition initiale
    yn.Fill(0);
    var_time.SetInitialVector(var_time.GetInitialTime(), yn);
    

    /**************************************************
     * set initial condition for the implicit part
     **************************************************/

    
    for(int i = 0; i < Nproche; i++)
      y_proche(i) = yn(dof_proche(i));

    DistributedMatrix<Real_wp, General, ArrayRowSparse> KhP;
    
    {
      DistributedMatrix<Real_wp, General, ArrayRowSparse> C, D;
      VectReal_wp invA11_sub(A11_sub);
      for (int i = 0; i < invA11_sub.GetM(); i++)
	invA11_sub(i) = Real_wp(1) / A11_sub(i);
      
      C = A12_sub;
      ScaleLeftMatrix(C, invA11_sub);
      
      Mlt(A21_sub, C, D);
      Add(Real_wp(-1), A22_sub, D);
      
      VectReal_wp invMass_sub(Nproche);
      for (int i = 0; i < Nproche; i++)
	invMass_sub(i) = invMass(dof_proche(i));
      
      ScaleLeftMatrix(D, invMass_sub);

      KhP.Reallocate(Nproche, Nproche);
      for (int i = 0; i < Nproche; i++)
	{
	  KhP.ReallocateRow(i, D.GetRowSize(i));
	  for (int j = 0; j < D.GetRowSize(i); j++)
	    {
	      int dof_glob = ddl_vol_imp(D.Index(i, j));
	      KhP.Index(i, j) = IndexProche(dof_glob);
	      KhP.Value(i, j) = D.Value(i, j);
	    }
	}
    }

    //KhP.WriteText("KhP.dat");
    
    VectorialFunctionFirst sys(KhP, Nproche, dt, points_ci);
    //VectorialFunctionFirst sys(KhP, N);

    VirtualTimeScheme<Real_wp>* implicit_scheme = NULL;
    if (data.local_implicit)
      {
	if (data.pade_scheme)
	  implicit_scheme = new PadeScheme_Iterator<Real_wp>(data.order_imp, true);
	else
	  implicit_scheme =
	    new LinearSdirkScheme_Iterator<Real_wp>(data.order_imp, data.s_imp,
						    LinearSdirkScheme_Iterator<Real_wp>::STABLE_WEIGHTS);
	
	implicit_scheme->SetInitialCondition(var_time.GetInitialTime(), dt, y_proche, sys);
	
	//Reallocate y_proche it has been clearup in SetInitialCondition
	y_proche.Reallocate(Nproche);
	
	//VectReal_wp y0(yn);
	// scheme.SetInitialCondition(var_time.GetInitialTime(), dt, y0, sys);
      }
    else
      {
	/*****************************************************
	 * set initial condition for the local explicit case
	 *****************************************************/
	explicit_scheme.SetInitialCondition(var_time.GetInitialTime(), dt/data.p, y_proche, sys);
	//Reallocate y_proche it has been clearup in SetInitialCondition
	y_proche.Reallocate(Nproche);
      }
    
    VectReal_wp coef_alpha(orderN);
    
    // matrix used to recover q^{(l)}
    // we divide by 1 / dt^l
    Matrix<Real_wp> DerMat = explicit_scheme.GetDerivativeMatrix();
    for (int l = 1; l < DerMat.GetM(); l++)
      {
	Real_wp coef = Real_wp(1) / pow(dt, l);
	for (int j = 0; j < DerMat.GetN(); j++)
	  DerMat(l, j) *= coef;
      }

    VectReal_wp L_imp(nb_dof_L_proche), L_far(nb_dof_L_lointain);
    L_imp.Zero(); L_far.Zero();
    
    // main iteration
    // loop in time
    for (int nt = 0; nt < nb_iter; nt++)
      {
	// on ecrit l'instantane
	Real_wp t = var_time.GetInitialTime() + dt*nt;	
	var_time.WriteSnapshot(nt, t, yn);
	
	// on calcule la source aux points d'interpolation
	bool presence_source = false;
	if (t < var_time.GetFinalTimeSource())
	  {
	    presence_source = true;
	    for (int i = 0; i < points_ci.GetM(); i++)
	      {
		EvalF(i).Zero();
		var_time.AddPrimitiveTimeSource(Real_wp(1), t + points_ci(i)*dt, 0, EvalF(i));

		for (int j = 0; j < Dh.GetM(); j++)
		  EvalF(i)(j) /= Dh(j);
		
		// on extrait la partie implicite de F
		for (int j = 0; j < Nproche; j++)
		  EvalFproche(i)(j) = EvalF(i)(dof_proche(j));
	      }    
	  }
	
	// Calcul des w_n_j
	w_tild = yn;
	for (int j = 0; j < orderN; j++)
	  {
 	    // on extrait la partie implicite de w_tild
	    for (int i = 0; i < Nimp; i++)
	      y_imp(i) = w_tild(ddl_vol_imp(i));

	    // w = A (I-P) w_tild
	    for (int i = 0; i < Nlointain; i++)
	      y_far(i) = w_tild(dof_lointain(i));

	    Mlt(A12_far, y_far, L_far);
	    for (int k = 0; k < L_far.GetM(); k++)
	      L_far(k) /= A11_far(k);

	    Mlt(A21_far, L_far, w);
	    MltAdd(-Real_wp(1), A22_far, y_far, Real_wp(1), w);

	    for (int k = 0; k < w.GetM(); k++)
	      w(k) *= invMass(k);
	    
	    //  wp = AP w_tild
	    Mlt(A12_sub, y_imp, L_imp);
	    for (int k = 0; k < L_imp.GetM(); k++)
	      L_imp(k) /= A11_sub(k);

	    Mlt(A21_sub, L_imp, wp);
	    MltAdd(-Real_wp(1), A22_sub, y_imp, Real_wp(1), wp);

	    for (int k = 0; k < wp.GetM(); k++)
	      wp(k) *= invMass(dof_proche(k));

	    //w.WriteText("w_ref.dat"); wp.WriteText("wp_ref.dat");
	    //int test_input; cout << "Waiting" << endl; cin >> test_input;
	    
	    // on stocke w_n_j
	    w_n_i(j) = w;
	    for(int i = 0; i < Nproche; i++)
	      w(dof_proche(i)) += wp(i);
	    
	    Real_wp coef = (j+1)*P(j+1);
	    Mlt(coef, w_n_i(j));
	    coef_alpha(j) = Real_wp(1)/(j+1);
	    
	    // on rajoute la source
	    if (presence_source)
	      if (j < orderN-1)
		for (int i = 0; i < points_ci.GetM(); i++)
		  w += DerMat(j, i)*EvalF(i);
	    
	    // on passe a l'itere suivant
	    w_tild = w;
	  }
	
	// on extrait wp_0, wp_1, wp_2, wp_3, etc
	for(int i = 0; i < Nproche; i++)
	  {
	    int ip = dof_proche(i);
	    for (int j = 0; j < orderN; j++)
	      wp_i(j)(i) = w_n_i(j)(ip);
	  }

	if (data.local_implicit)
	  {
	    /********************************************************
	     * Fine part with implicit scheme
	     *******************************************************/

	    for(int i = 0; i < Nproche; i++)
	      y_proche(i) = yn(dof_proche(i));

	    sys.SetSource(t, wp_i, presence_source, EvalFproche);
	    implicit_scheme->GetIterate() = y_proche;
	    
	    implicit_scheme->Advance(t, nt, sys);
	    y_proche = implicit_scheme->GetIterate();

	    /*
	    sys.SetSource(t, w_n_0, w_n_1, w_n_2, w_n_3);
	    scheme.GetIterate() = yn;

	    scheme.Advance(t, nt, sys);

	    yn = scheme.GetIterate(); */
	  }
	else
	  {
	    /********************************************************
	     * Fine part with small time step explicit scheme
	     *******************************************************/
	    // solution de dy/dtau = patati dans la zone proche
	    int p = data.p;
	    Real_wp dtau = dt/p;
	    //DISP(p); DISP(dt); DISP(dtau);

	    for(int i = 0; i < Nproche; i++)
	      y_proche(i) = yn(dof_proche(i));
	    

	    //Using optimal linear Runge Kutta schemes
	    sys.SetSource(t, wp_i, presence_source, EvalFproche);
	    explicit_scheme.GetIterate() = y_proche;
	    
	    for(int m = 0; m < p; m++)
	      {
		explicit_scheme.Advance(t, nt, sys);
		t = t+dtau;
	      }
	    y_proche = explicit_scheme.GetIterate();
	    
	    /*/Using Classical RK4 like in Grote's paper With no source
	    VectReal_wp k1(Nproche), k2(Nproche), k3(Nproche), k4(Nproche);

	    for(int m = 0; m < p; m++)
	      {
		Mlt(Kh_sub, y_proche, k1);
		k1 += (((dtau*m/3)*wp_i(3) + wp_i(2)) * (dtau*m/2) + wp_i(1)) * (dtau*m) + wp_i(0);
		
		y_tmp = y_proche + (dtau/2)*k1;
		Mlt(Kh_sub, y_tmp, k2);
		Real_wp mh = m + 0.5;
		k2 += (((dtau*mh/3)*wp_i(3) + wp_i(2)) * (dtau*mh/2) + wp_i(1)) * (dtau*mh) + wp_i(0);
		
		y_tmp = y_proche + (dtau/2)*k2;
		Mlt(Kh_sub, y_tmp, k3);
		k3 += (((dtau*mh/3)*wp_i(3) + wp_i(2)) * (dtau*mh/2) + wp_i(1)) * (dtau*mh) + wp_i(0);
		
		y_tmp = y_proche + dtau*k3;
		Mlt(Kh_sub, y_tmp, k4);
		mh = m + 1.0;
		k4 += (((dtau*mh/3)*wp_i(3) + wp_i(2)) * (dtau*mh/2) + wp_i(1)) * (dtau*mh) + wp_i(0);
		
		y_proche += (dtau/6)*(k1 + Real_wp(2)*k2 + Real_wp(2)*k3 + k4);
		}*/
	    
	  }

	// on met y_proche dans yn
	for(int i = 0; i < Nproche; i++)
	  yn(dof_proche(i)) = y_proche(i); 
	
	// on calcule yn sur les ddls lointains
	for (int i = 0; i < nb_dof_uv; i++)
	  if (IndexProche(i) < 0)
	    {
	      Real_wp vloc = coef_alpha(orderN-1)*w_n_i(orderN-1)(i);
	      for (int j = orderN-2; j >= 0; j--)
		vloc = dt*vloc + coef_alpha(j)*w_n_i(j)(i);

	      for (int j = 0; j < bi.GetN(); j++)
		vloc += bi(0, j)*EvalF(j)(i);
	      
	      yn(i) += dt*vloc;
	    }

	//yn.WriteText("Unext_ref.dat");
	//y_proche.WriteText("yp_ref.dat");
	//int test_input; cout << "Waiting Final " << endl; cin >> test_input;
	
      }
    
  if (data.local_implicit)
    delete implicit_scheme;
  }
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc == 2)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      
      HyperbolicProblem<AcousticEquationDG<Dimension2> > Vars;
      
      cout<<" Acoustic Solver with Lobatto quadrangles Qk"<<endl; 
      
      RunAll(Vars, file_name_data, type_element, string("ACOUSTIC_HDG"));
      cout<<" we destroy the variables "<<endl;
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"test.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();

}
