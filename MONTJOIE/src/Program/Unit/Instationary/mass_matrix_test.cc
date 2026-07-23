#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG
#define MONTJOIE_WITH_LEGENDRE_DG

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"
#include "Hyperbolic/Maxwell/MontjoieTimeMaxwell3D.hxx"

using namespace Montjoie;

Real_wp threshold;
bool fast_check = true;

void Copy(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
	  DistributedMatrix<Complex_wp, General, ArrayRowSparse>& B)
{
  B.Clear();
  B.Reallocate(A.GetM(), A.GetN());
  for (int i = 0; i < B.GetM(); i++)
    {
      int s = A.GetRowSize(i);
      B.ReallocateRow(i, s);
      for (int j = 0; j < s; j++)
	{
	  B.Index(i, j) = A.Index(i, j);
	  B.Value(i, j) = A.Value(i, j);
	}
    }
}

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  x = complex<T>(rand(), rand())/RAND_MAX;
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Clear();
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRandNumber(x);
      A.Set(i, j, x);
    }
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    {
      DISP(x.GetM()); DISP(y.GetM());
      return false;
    }

  Real_wp yref = y.GetNormInf();
  if ((Norm2(x) == 0) && (Norm2(y) == 0))
    return true;
  
  if ((Norm2(x) <= eps) || isnan(yref))
    {
      DISP(Norm2(x)); DISP(Norm2(y));
      return false;
    }

  for (int i = 0; i < x.GetM(); i++)
    if (isnan(abs(x(i) - y(i))) || (abs(x(i) - y(i))/yref > eps))
      {
	DISP(i); DISP(x(i)); DISP(y(i));
	return false;
      }
  
  return true;
}

template<class TypeEquation>
void CheckBlockMass(HyperbolicProblem<TypeEquation>& var, GhostIf<false>& )
{
}

template<class TypeEquation>
void CheckBlockMass(HyperbolicProblem<TypeEquation>& var, GhostIf<true>&)
{
}

template<class T, class TypeEquation>
void GenerateRandomVectorPb(const EllipticProblem<TypeEquation>& var, Vector<T>& u, bool scalar, bool global_vec = false)
{
  int nb_dof_u = var.nb_unknowns_scal*var.GetNbScalarDof();
  int nb_dof_v = var.GetNbDof() - nb_dof_u;
  
  int Nu = var.GetNbGlobalScalarDof()*var.nb_unknowns_scal;
  int Nv = var.GetNbGlobalDof() - Nu;
  
  int nb_proc = var.GetNbProcPerMode();
  if (nb_proc == 1)
    {
      Nu = nb_dof_u;
      Nv = nb_dof_v;
    }
  
  Vector<T> u_glob;
  if (global_vec)
    GenerateRandomVector(u_glob, Nu + Nv);
  else
    {
      if (scalar)
	GenerateRandomVector(u_glob, Nu);
      else
	GenerateRandomVector(u_glob, Nv);
    }
  
  const Vector<int>& row_num = var.mesh_num.GlobDofNumber_Subdomain;
  
  Vector<Real_wp> coef_distrib(nb_dof_u);
  if (global_vec)
    coef_distrib.Reallocate(nb_dof_u + nb_dof_v);

  if (scalar || global_vec)
    {
      coef_distrib.Fill(1.0);
      var.AddDomains(coef_distrib);
      for (int i = 0; i < nb_dof_u; i++)
	coef_distrib(i) = 1.0/coef_distrib(i);
    }
  
  if (nb_proc == 1)
    u = u_glob;
  else
    {
      if (scalar || global_vec)
	for (int i = 0; i < coef_distrib.GetM(); i++)
	  u(i) = coef_distrib(i)*u_glob(row_num(i));
      else
	for (int i = 0; i < nb_dof_v; i++)
	  u(i) = u_glob(row_num(nb_dof_u + i)-Nu);
    }  
}

template<class TypeEquation>
void CheckProblem(const HyperbolicProblem<TypeEquation>& var_c, bool first_order, const string& name_equation,
		  const string& name_element, const string& input_file, bool split_system = true,
		  bool combined_mass = true)
{
  HyperbolicProblem<TypeEquation> var;
  
  var.var_harmonic.SetTypeEquation(name_equation);
  
  var.var_harmonic.InitIndices(100);
  ReadInputFile(input_file, var);
  
  var.glob_solver = var.var_harmonic.GetNewLinearSolver();
  ReadInputFile(input_file, *var.glob_solver);
  
  var.var_harmonic.SetFirstOrderFormulation(first_order);
  var.var_harmonic.ComputeMeshAndFiniteElement(name_element);

  int nb_proc = var.var_harmonic.GetNbProcPerMode();
  //int rank_proc = var.var_harmonic.GetRankProcMode();
  //string suffix = to_str(rank_proc) + ".dat";
  if (nb_proc == 1)  
    var.var_harmonic.mesh.Write("test.mesh");

  var.var_harmonic.PerformOtherInitializations();
  
  var.var_harmonic.ComputeMassMatrix();
  var.var_harmonic.ComputeQuasiPeriodicPhase();
  var.CheckIdentityJacobianInPML();

  VectReal_wp b_src;
  var.var_harmonic.ComputeRightHandSide(b_src);
  
  DistributedMatrix<Real_wp, General, ArrayRowSparse> Dh, Sh, Kh, invCh;
  GlobalGenericMatrix<Real_wp> nat_mat;
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(0.0);
  var.var_harmonic.SetCoefficientDirichlet(Real_wp(1));
  var.var_harmonic.AddMatrixWithBC(Dh, nat_mat);

  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(1.0);
  var.var_harmonic.SetCoefficientDirichlet(Real_wp(0));
  var.var_harmonic.AddMatrixWithBC(Sh, nat_mat);

  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(1.0);
  nat_mat.SetCoefDamping(0.0);
  var.var_harmonic.SetCoefficientDirichlet(Real_wp(0));
  var.var_harmonic.AddMatrixWithBC(Kh, nat_mat);

  int nb_dof_u = var.var_harmonic.nb_unknowns_scal*var.var_harmonic.GetNbScalarDof();
  int nb_dof_v = var.var_harmonic.GetNbDof() - nb_dof_u;
  int nb_dof_L = 0;

  if (var.var_harmonic.FormulationDG() == ElementReference_Base::HDG)
    {
      // Dh and Sh concerns only u and v part
      nb_dof_L = var.var_harmonic.nb_unknowns_hdg*var.var_harmonic.GetNbMeshDof();
      int nb_dof_uv = Dh.GetM() - nb_dof_L;
      nb_dof_u = var.var_harmonic.nb_unknowns_scal*var.var_harmonic.GetNbMainUnknownDof();
      nb_dof_v = var.var_harmonic.nb_unknowns_vec*var.var_harmonic.GetNbMainUnknownDof();

      Dh.WriteText("Dav.dat"); Sh.WriteText("Sav.dat"); Kh.WriteText("Kav.dat");
      
      DistributedMatrix<Real_wp, General, ArrayRowSparse> D(Dh), S(Sh), K(Kh), Cs;
      GetSubMatrix(D, nb_dof_L, nb_dof_L+nb_dof_uv, nb_dof_L, nb_dof_L+nb_dof_uv, Dh);
      GetSubMatrix(S, nb_dof_L, nb_dof_L+nb_dof_uv, nb_dof_L, nb_dof_L+nb_dof_uv, Sh);
      GetSubMatrix(K, nb_dof_L, nb_dof_L+nb_dof_uv, nb_dof_L, nb_dof_L+nb_dof_uv, Kh);

      // for invCh, we convert Ch to a block-diagonal matrix
      DistributedMatrix<Real_wp, General, BlockDiagRow> Ch;

      GetSubMatrix(K, 0, nb_dof_L, 0, nb_dof_L, Cs);
      ConvertToBlockDiagonal(Cs, Ch, 0, nb_dof_L);

      GetInverse(Ch);
      ConvertToSparse(Ch, invCh);

      // we compute Kh with a Schur complement
      DistributedMatrix<Real_wp, General, ArrayRowSparse> A, B, C;
      GetSubMatrix(K, 0, nb_dof_L, nb_dof_L, nb_dof_L+nb_dof_uv, A);
      GetSubMatrix(K, nb_dof_L, nb_dof_L+nb_dof_uv, 0, nb_dof_L, B);

      C.Reallocate(nb_dof_L, nb_dof_uv);
      Mlt(invCh, A, C);

      A.Reallocate(nb_dof_uv, nb_dof_uv);
      Mlt(B, C, A);

      Add(Real_wp(-1), A, Kh);

      // modification of right hand side
      VectReal_wp rhs_L(nb_dof_L), rhs_uv(nb_dof_uv);
      for (int i = 0; i < nb_dof_L; i++)
	rhs_L(i) = b_src(i);

      for (int i = 0; i < nb_dof_uv; i++)
	rhs_uv(i) = b_src(nb_dof_L + i);

      b_src = rhs_uv;
    }
  
  Dh.WriteText("Dh.dat"); Sh.WriteText("Sh.dat"); Kh.WriteText("Kh.dat");
  
  var.Glob_mat_Dh = var.GetNewMassMatrix();
  
  var.Glob_mat_Dh->compute_Dh_default = true;
  var.Glob_mat_Dh->compute_invDh_default = true;
  var.Glob_mat_Dh->compute_Sh_default = true;
  if (combined_mass)
    {
      var.Glob_mat_Dh->compute_DhMinusdtSh_default = true;
      var.Glob_mat_Dh->compute_invDhPlusdtSh_default = true;
    }

  var.Glob_mat_Dh->compute_Bh_default = true;
  var.Glob_mat_Dh->compute_invBh_default = true;
  if (combined_mass)
    {
      var.Glob_mat_Dh->compute_ShVec_default = true;
      var.Glob_mat_Dh->compute_BhMinusdtSh_default = true;
      var.Glob_mat_Dh->compute_invBhPlusdtSh_default = true;
    }

  var.ComputeRightHandSide();
  var.ComputeStiffnessMatrix();
  var.ComputeMassMatrix();
  var.InitTimeIterations();

  srand(Seed);
  
  Real_wp one(1), zero(0);
  Real_wp alpha; GetRandNumber(alpha);
  Real_wp beta; GetRandNumber(beta);
  Real_wp gamma; GetRandNumber(gamma);


  DISP(nb_dof_u); DISP(nb_dof_v);
  VectReal_wp Ones(nb_dof_v), ProdUh(nb_dof_u), ProdVh, Urand(nb_dof_u), Vrand(nb_dof_v);
  VectReal_wp Y(nb_dof_u);

  if (nb_proc > 1)
    {
      GenerateRandomVectorPb(var.var_harmonic, Urand, true);
      GenerateRandomVectorPb(var.var_harmonic, Vrand, false);
    }
  else
    {
      GenerateRandomVector(Urand, nb_dof_u);
      GenerateRandomVector(Vrand, nb_dof_v);
    }

  DistributedMatrix<Real_wp, General, ArrayRowSparse> DhScal, ShScal, RhS, RhV, ShVec, BhVec;
  if (fast_check)
    {
      if (first_order)
	{
	  GetSubMatrix(Dh, 0, nb_dof_u, 0, nb_dof_u, DhScal);
	  GetSubMatrix(Sh, 0, nb_dof_u, 0, nb_dof_u, ShScal);
	  GetSubMatrix(Dh, nb_dof_u, nb_dof_u+nb_dof_v, nb_dof_u, nb_dof_u+nb_dof_v, BhVec);
	  GetSubMatrix(Sh, nb_dof_u, nb_dof_u+nb_dof_v, nb_dof_u, nb_dof_u+nb_dof_v, ShVec);
	  GetSubMatrix(Kh, 0, nb_dof_u, nb_dof_u, nb_dof_u+nb_dof_v, RhS);
	  GetSubMatrix(Kh, nb_dof_u, nb_dof_u+nb_dof_v, 0, nb_dof_u, RhV);
	}
      else
	{
	  DhScal = Dh;
	  ShScal = Sh;
	}
    }
  else
    {
      if (nb_proc > 1)
	{
	  cout << "Not implemented" << endl;
	  abort();
	}
    }
  
  //
  // common functions (first or second-order formulation)
  //

  // testing ApplyOperatorSh
  ProdUh.Reallocate(nb_dof_u); Ones.Reallocate(nb_dof_u);
  if (fast_check)
    {
      if (nb_proc > 1)
	GenerateRandomVectorPb(var.var_harmonic, Ones, true);
      else
	GenerateRandomVector(Ones, nb_dof_u);
      
      var.ApplyOperatorSh(one, zero, Ones, zero, ProdUh);
      Mlt(ShScal, Ones, Y, false);
      if (!EqualVector(ProdUh, Y))
	{
	  cout << "ApplyOperatorSh incorrect" << endl;
	  abort();
	}

      ProdUh = Urand;
      var.ApplyOperatorSh(alpha, zero, Ones, beta, ProdUh);

      Y = Urand;
      MltAdd(alpha, ShScal, Ones, beta, Y, false);
      if (!EqualVector(ProdUh, Y))
	{
	  cout << "ApplyOperatorSh incorrect" << endl;
	  abort();
	}
    }
  else
    for (int i = 0; i < nb_dof_u; i++)
      {
	Ones.Fill(0); Ones(i) = 1.0;
	var.ApplyOperatorSh(one, zero, Ones, zero, ProdUh);
	for (int j = 0; j < nb_dof_u; j++)
	  {
	    Real_wp val = Sh(j, i);
	    if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
	      {
		cout << "ApplyOperatorSh incorrect" << endl;
		DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		abort();
	      }
	  }
	
	ProdUh = Urand;
	var.ApplyOperatorSh(alpha, zero, Ones, beta, ProdUh);
	for (int j = 0; j < nb_dof_u; j++)
	  {
	    Real_wp val = beta*Urand(j) + alpha*Sh(j, i);
	    if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
	      {
		cout << "ApplyOperatorSh incorrect" << endl;
		DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		abort();
	      }
	  }
      }
  
  // testing ApplyOperatorDh
  if (fast_check)
    {
      var.ApplyOperatorDh(one, zero, Ones, zero, ProdUh);
      Mlt(DhScal, Ones, Y, false);
      if (!EqualVector(ProdUh, Y))
	{
	  cout << "ApplyOperatorDh incorrect" << endl;
	  abort();
	}

      ProdUh = Urand;
      var.ApplyOperatorDh(alpha, zero, Ones, beta, ProdUh);

      Y = Urand;
      MltAdd(alpha, DhScal, Ones, beta, Y, false);
      if (!EqualVector(ProdUh, Y))
	{
	  cout << "ApplyOperatorDh incorrect" << endl;
	  abort();
	}
    }
  else
    for (int i = 0; i < nb_dof_u; i++)
      {
	Ones.Fill(0); Ones(i) = 1.0;
	var.ApplyOperatorDh(one, zero, Ones, zero, ProdUh);
	for (int j = 0; j < nb_dof_u; j++)
	  {
	    Real_wp val = Dh(j, i);
	    if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
	      {
		cout << "ApplyOperatorDh incorrect" << endl;
		DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		abort();
	      }
	  }
	
	ProdUh = Urand;
	var.ApplyOperatorDh(alpha, zero, Ones, beta, ProdUh);
	for (int j = 0; j < nb_dof_u; j++)
	  {
	    Real_wp val = beta*Urand(j) + alpha*Dh(j, i);
	    if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
	      {
		cout << "ApplyOperatorDh incorrect" << endl;
		DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		abort();
	      }
	  }
      }
    
  // testing SolveOperatorDh
  if (fast_check)
    {
      Y = Ones;
      var.ApplyOperatorDh(one, zero, Y, zero, ProdUh);
      var.SolveOperatorDh(ProdUh);
      if (!EqualVector(Ones, ProdUh))
	{
	  cout << "SolveOperatorDh incorrect" << endl;
	  abort();
	}
    }
  else
    for (int i = 0; i < nb_dof_u; i++)
      {
	Ones.Fill(0); Ones(i) = 1.0;
	Y = Ones;
	var.SolveOperatorDh(Y);
	var.ApplyOperatorDh(Real_wp(1), zero, Y, Real_wp(0), ProdUh);
	for (int j = 0; j < nb_dof_u; j++)
	  if (isnan(abs(ProdUh(j) - Ones(j))) || (abs(ProdUh(j) - Ones(j)) > 1e-10))
	    {
	      cout << "SolveOperatorDh incorrect" << endl;
	      DISP(i); DISP(j); DISP(ProdUh(j)); DISP(Ones(j));
	      abort();
	    }
      }
  
  // testing ApplyOperatorDhMinusdtSh
  if (combined_mass)
    {
      if (fast_check)
	{
	  var.ApplyOperatorDhMinusdtSh(one, zero, Ones, zero, ProdUh);
	  Mlt(DhScal, Ones, Y, false);
	  MltAdd(-0.5*var.GetTimeStep(), ShScal, Ones, one, Y, false);
	  if (!EqualVector(Y, ProdUh))
	    {
	      cout << "ApplyOperatorDhMinusdtSh incorrect" << endl;
	      abort();
	    }
	  
	  ProdUh = Urand;
	  var.ApplyOperatorDhMinusdtSh(alpha, zero, Ones, beta, ProdUh);
	  
	  Y = Urand;
	  MltAdd(alpha, DhScal, Ones, beta, Y, false);
	  MltAdd(-0.5*alpha*var.GetTimeStep(), ShScal, Ones, one, Y, false);
	  
	  if (!EqualVector(Y, ProdUh))
	    {
	      cout << "ApplyOperatorDhMinusdtSh incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_u; i++)
	  {
	    Ones.Zero(); Ones(i) = 1.0;      
	    var.ApplyOperatorDhMinusdtSh(one, zero, Ones, zero, ProdUh);
	    for (int j = 0; j < nb_dof_u; j++)
	      {
		Real_wp val = Dh(j, i) - 0.5*var.GetTimeStep()*Sh(j, i);
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "ApplyOperatorDhMinusdtSh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	    
	    ProdUh = Urand;
	    var.ApplyOperatorDhMinusdtSh(alpha, zero, Ones, beta, ProdUh);
	    for (int j = 0; j < nb_dof_u; j++)
	      {
		Real_wp val = beta*Urand(j) + alpha*(Dh(j, i) - 0.5*var.GetTimeStep()*Sh(j, i));
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "ApplyOperatorDhMinusdtSh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	  }
      
      
      // testing SolveOperatorDhPlusdtSh
      if (fast_check)
	{
	  Y = Ones;
	  var.ApplyOperatorDh(one, zero, Y, zero, ProdUh);
	  var.ApplyOperatorSh(Real_wp(0.5)*var.GetTimeStep(), zero, Y, one, ProdUh);
	  var.SolveOperatorDhPlusdtSh(ProdUh);
	  if (!EqualVector(ProdUh, Ones))
	    {
	      cout << "SolveOperatorDhPlusdtSh incorrect" << endl;
	      abort();
	    }
	  
	  var.ApplyOperatorDh(one/alpha, zero, Ones, zero, ProdUh);
	  var.ApplyOperatorSh(Real_wp(0.5)*var.GetTimeStep()/alpha, zero, Ones, one, ProdUh);
	  
	  Y = Urand;
	  var.SolveOperatorDhPlusdtSh(alpha, ProdUh, Y);
	  Y -= Urand;
	  
	  if (!EqualVector(Y, Ones))
	    {
	      cout << "SolveOperatorDhPlusdtSh incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_u; i++)
	  {
	    Ones.Zero(); Ones(i) = 1.0;  Y = Ones;
	    var.SolveOperatorDhPlusdtSh(Y);
	    var.ApplyOperatorDh(Real_wp(1), zero, Y, Real_wp(0), ProdUh);
	    var.ApplyOperatorSh(Real_wp(0.5)*var.GetTimeStep(), zero, Y, Real_wp(1), ProdUh);
	    for (int j = 0; j < nb_dof_u; j++)
	      if (isnan(abs(ProdUh(j) - Ones(j))) || abs(ProdUh(j) - Ones(j)) > 1e-10)
		{
		  cout << "SolveOperatorDhPlusdtSh incorrect" << endl;
		  DISP(i); DISP(j); DISP(ProdUh(j)); DISP(Ones(j));
		  abort();
		}
	    
	    Y = Urand;
	    var.SolveOperatorDhPlusdtSh(alpha, Ones, Y);
	    Y -= Urand;
	    var.ApplyOperatorDh(Real_wp(1), zero, Y, Real_wp(0), ProdUh);
	    var.ApplyOperatorSh(Real_wp(0.5)*var.GetTimeStep(), zero, Y, Real_wp(1), ProdUh);
	    for (int j = 0; j < nb_dof_u; j++)
	      {
		Real_wp val = alpha*Ones(j);
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "SolveOperatorDhPlusdtSh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	  }
    }
  
  if (!first_order)
    {
      // testing ApplyOperatorKh
      ProdUh.Reallocate(nb_dof_u); Ones.Reallocate(nb_dof_u);
      if (fast_check)
	{
	  var.ApplyOperatorKh(-one, zero, Ones, zero, ProdUh);
	  var.var_harmonic.ImposeNullDirichletCondition(ProdUh);
	  
	  Mlt(Kh, Ones, Y, false);
	  if (!EqualVector(ProdUh, Y))
	    {
	      cout << "ApplyOperatorKh incorrect" << endl;
	      abort();
	    }

	  ProdUh = Urand;
	  var.ApplyOperatorKh(-alpha, zero, Ones, beta, ProdUh);
	  for (int n2 = 0; n2 < var.var_harmonic.GetNbDirichletDof(); n2++)
	    {
	      int n = var.var_harmonic.GetDirichletDofNumber(n2);
	      ProdUh(n) = Urand(n)*beta;
	    }
	  
	  Y = Urand;
	  MltAdd(alpha, Kh, Ones, beta, Y, false);
	  if (!EqualVector(ProdUh, Y))
	    {
	      cout << "ApplyOperatorKh incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_u; i++)
	  {
	    Ones.Fill(0); Ones(i) = 1.0;
	    var.ApplyOperatorKh(-one, zero, Ones, zero, ProdUh);
	    var.var_harmonic.ImposeNullDirichletCondition(ProdUh);
	    for (int j = 0; j < nb_dof_u; j++)
	      {
		Real_wp val = Kh(j, i);
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "ApplyOperatorKh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	    
	    ProdUh = Urand;
	    var.ApplyOperatorKh(-alpha, zero, Ones, beta, ProdUh);
	    for (int n2 = 0; n2 < var.var_harmonic.GetNbDirichletDof(); n2++)
	      {
		int n = var.var_harmonic.GetDirichletDofNumber(n2);
		ProdUh(n) = Urand(n)*beta;
	      }
	    
	    for (int j = 0; j < nb_dof_u; j++)
	      {
		Real_wp val = beta*Urand(j) + alpha*Kh(j, i);
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "ApplyOperatorKh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	  }
    }
  else
    {	
      // testing ApplyOperatorRhScalar
      Ones.Reallocate(nb_dof_v); ProdUh.Reallocate(nb_dof_u);
      ProdUh.Zero();
      if (split_system)
	{
	  if (fast_check)
	    {
	      if (nb_proc > 1)
		GenerateRandomVectorPb(var.var_harmonic, Ones, false);
	      else
		GenerateRandomVector(Ones, nb_dof_v);
	      
	      ProdUh = Urand;
	      var.ApplyOperatorRhScalar(alpha, zero, Ones, beta, ProdUh);
	      for (int n2 = 0; n2 < var.var_harmonic.GetNbDirichletDof(); n2++)
		{
		  int n = var.var_harmonic.GetDirichletDofNumber(n2);
		  ProdUh(n) = Urand(n)*beta;
		}
	      
	      Y = Urand;
	      MltAdd(-alpha, RhS, Ones, beta, Y, false);
	      if (!EqualVector(Y, ProdUh))
		{
		  cout << "ApplyOperatorRhScalar incorrect" << endl;
		  abort();
		}
	    }
	  else
            for (int i = 0; i < nb_dof_v; i++)
	      {
		Ones.Zero(); Ones(i) = 1.0;
		ProdUh = Urand;
		var.ApplyOperatorRhScalar(alpha, zero, Ones, beta, ProdUh);
		for (int j = 0; j < nb_dof_u; j++)
		  if ((var.var_harmonic.GetNbDirichletDof() == 0) ||
		      (!var.var_harmonic.IsDofDirichlet(j)))
		    {
		      Real_wp val = beta*Urand(j) - alpha*Kh(j, nb_dof_u+i);
		      if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
			{
			  cout << "ApplyOperatorRhScalar incorrect" << endl;
			  DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
			  DISP(beta*Urand(j)); DISP(alpha);
			  abort();
			}
		    }      
	      }
	}
      
      // testing ApplyOperatorRhVectorial
      ProdVh.Reallocate(nb_dof_v);
      Ones.Reallocate(nb_dof_u);
      ProdVh.Zero();
      if (split_system)
	{
	  //SetComplexZero(beta); SetComplexOne(alpha);
	  if (fast_check)
	    {
	      if (nb_proc > 1)
		GenerateRandomVectorPb(var.var_harmonic, Ones, true);
	      else
		GenerateRandomVector(Ones, nb_dof_u);
	      
	      ProdVh = Vrand;
	      var.ApplyOperatorRhVectorial(alpha, zero, Ones, beta, ProdVh);
	      
	      Y = Vrand;
	      MltAdd(-alpha, RhV, Ones, beta, Y, false);
	      if (!EqualVector(Y, ProdVh))
		{
		  cout << "ApplyOperatorRhVectorial incorrect" << endl;
		  abort();
		}
	    }
	  else	    
	    for (int i = 0; i < nb_dof_u; i++)
	      {
		Ones.Zero(); Ones(i) = 1.0;
		ProdVh = Vrand;
		var.ApplyOperatorRhVectorial(alpha, zero, Ones, beta, ProdVh);
		for (int j = 0; j < nb_dof_v; j++)
		  {
		    Real_wp val = beta*Vrand(j) - alpha*Kh(nb_dof_u+j, i);
		    if (isnan(abs(ProdVh(j) - val)) || (abs(ProdVh(j) - val) > 1e-10))
		      {
			cout << "ApplyOperatorRhVectorial incorrect" << endl;
			DISP(i); DISP(j); DISP(ProdVh(j)); DISP(val);
			DISP(beta*Vrand(j)); DISP(alpha); DISP(Kh(nb_dof_u+j, i));
			abort();
		      }
		  }
	      }
	}
      
      // testing ApplyOperatorShVectorial
      ProdUh.Reallocate(nb_dof_v); Ones.Reallocate(nb_dof_v);
      cout << "Testing ApplyOperatorShVectorial" << endl;
      if (combined_mass)
	{
	  if (fast_check)
	    {
	      if (nb_proc > 1)
		GenerateRandomVectorPb(var.var_harmonic, Ones, false);
	      else
		GenerateRandomVector(Ones, nb_dof_v);

	      ProdUh = Vrand;
	      var.ApplyOperatorShVectorial(alpha, zero, Ones, beta, ProdUh);
	      
	      Y = Vrand;
	      MltAdd(alpha, ShVec, Ones, beta, Y, false);
	      if (!EqualVector(Y, ProdUh))
		{
		  cout << "ApplyOperatorShVectorial incorrect" << endl;
		  abort();
		}
	    }
	  else
	    for (int i = 0; i < nb_dof_v; i++)
	      {
		Ones.Zero(); Ones(i) = 1.0;
		ProdUh = Vrand;
		var.ApplyOperatorShVectorial(alpha, zero, Ones, beta, ProdUh);
		for (int j = 0; j < nb_dof_v; j++)
		  {
		    Real_wp val = beta*Vrand(j) + alpha*Sh(nb_dof_u+j, nb_dof_u+i);
		    if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		      {
			cout << "ApplyOperatorShVectorial incorrect" << endl;
			DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
			abort();
		      }
		  }
	      }
	}
      
      // testing ApplyOperatorBh/SolveOperatorBh
      cout << "Testing ApplyOperatorBh" << endl;
      if (fast_check)
	{	  
	  var.ApplyOperatorBh(one, zero, Ones, zero, ProdUh);
	  Mlt(BhVec, Ones, Y, false);

	  if (!EqualVector(Y, ProdUh))
	    {
	      cout << "ApplyOperatorBh incorrect" << endl;
	      abort();
	    }

	  var.SolveOperatorBh(ProdUh);
	  if (!EqualVector(Ones, ProdUh))
	    {
	      cout << "SolveOperatorBh incorrect" << endl;
	      abort();
	    }

	  ProdUh = Vrand;
	  var.ApplyOperatorBh(alpha, zero, Ones, beta, ProdUh);

	  Y = Vrand;
	  MltAdd(alpha, BhVec, Ones, beta, Y, false);

	  if (!EqualVector(Y, ProdUh))
	    {
	      cout << "ApplyOperatorBh incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_v; i++)
	  {
	    Ones.Zero(); Ones(i) = 1.0;
	    var.ApplyOperatorBh(one, zero, Ones, zero, ProdUh);
	    for (int j = 0; j < nb_dof_v; j++)
	      {
		Real_wp val = Dh(nb_dof_u+j, nb_dof_u+i);
		if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		  {
		    cout << "ApplyOperatorBh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	    
	    var.SolveOperatorBh(ProdUh);
	    for (int j = 0; j < nb_dof_v; j++)
	      if (isnan(abs(ProdUh(j) - Ones(j))) || (abs(ProdUh(j) - Ones(j)) > 1e-10))
		{
		  cout << "SolveOperatorBh incorrect" << endl;
		  DISP(i); DISP(j); DISP(ProdUh(j)); DISP(Ones(j));
		  abort();
		}
	    
	    ProdUh = Vrand;
	    var.ApplyOperatorBh(alpha, zero, Ones, beta, ProdUh);
	    for (int j = 0; j < nb_dof_v; j++)
	      {
		Real_wp val = beta*Vrand(j) + alpha*Dh(nb_dof_u+j, nb_dof_u+i);
		if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		  {
		    cout << "ApplyOperatorBh incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	  }
      
      
      // testing ApplyOperatorBhMinusdtSh / SolveOperatorBhPlusdtSh
      if (combined_mass)
	{
	  if (fast_check)
	    {
	      var.ApplyOperatorBhMinusdtSh(one, zero, Ones, zero, ProdUh);
	      Mlt(BhVec, Ones, Y, false);
	      MltAdd(-0.5*var.GetTimeStep(), ShVec, Ones, one, Y, false);
	      if (!EqualVector(Y, ProdUh))
		{
		  cout << "ApplyOperatorBhMinusdtSh incorrect" << endl;
		  abort();
		}
	      
	      var.ApplyOperatorShVectorial(var.GetTimeStep(), zero, Ones, one, ProdUh);
	      var.SolveOperatorBhPlusdtSh(ProdUh);
	      
	      if (!EqualVector(Ones, ProdUh))
		{
		  cout << "SolveOperatorBhPlusdtSh incorrect" << endl;
		  abort();
		}
	      
	      ProdUh = Vrand;
	      var.ApplyOperatorBhMinusdtSh(alpha, zero, Ones, beta, ProdUh);
	      
	      Y = Vrand;
	      MltAdd(alpha, BhVec, Ones, beta, Y, false);
	      MltAdd(-0.5*alpha*var.GetTimeStep(), ShVec, Ones, one, Y, false);
	      if (!EqualVector(Y, ProdUh))
		{
		  cout << "ApplyOperatorBhMinusdtSh incorrect" << endl;
		  abort();
		}
	    }
	  else
	    for (int i = 0; i < nb_dof_v; i++)
	      {
		Ones.Zero(); Ones(i) = Real_wp(1);
		var.ApplyOperatorBhMinusdtSh(one, zero, Ones, zero, ProdUh);
		for (int j = 0; j < nb_dof_v; j++)
		  {
		    Real_wp val = Dh(nb_dof_u+j, nb_dof_u+i) - 0.5*var.GetTimeStep()*Sh(nb_dof_u+j, nb_dof_u+i);
		    if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		      {
			cout << "ApplyOperatorBhMinusdtSh incorrect" << endl;
			DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
			abort();
		      }
		  }
		
		var.ApplyOperatorShVectorial(var.GetTimeStep(), zero, Ones, one, ProdUh);
		var.SolveOperatorBhPlusdtSh(ProdUh);
		for (int j = 0; j < nb_dof_v; j++)
		  if (isnan(abs(ProdUh(j) - Ones(j))) || (abs(ProdUh(j) - Ones(j)) > 1e-10))
		    {
		      cout << "SolveOperatorBhPlusdtSh incorrect" << endl;
		      DISP(i); DISP(j); DISP(ProdUh(j)); DISP(Ones(j));
		      abort();
		    }
		
		ProdUh = Vrand;
		var.ApplyOperatorBhMinusdtSh(alpha, zero, Ones, beta, ProdUh);
		for (int j = 0; j < nb_dof_v; j++)
		  {
		    Real_wp val = beta*Vrand(j) + alpha*(Dh(nb_dof_u+j, nb_dof_u+i) - 0.5*var.GetTimeStep()*Sh(nb_dof_u+j, nb_dof_u+i));
		    if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		      {
			cout << "ApplyOperatorBhMinusdtSh incorrect" << endl;
			DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
			abort();
		      }
		  }
	      }
	}
    }
  
  // testing AddPrimitiveSource
  Urand.Reallocate(nb_dof_u + nb_dof_v);
  ProdUh.Reallocate(nb_dof_u + nb_dof_v);
  
  if (nb_proc > 1)
    GenerateRandomVectorPb(var.var_harmonic, Urand, true, true);
  else
    GenerateRandomVector(Urand, nb_dof_u + nb_dof_v);

  Real_wp tn = Real_wp(0.9);
  Real_wp pulse; var.SourceOnlyTime(tn, 0, pulse);
  
  ProdUh = Urand;
  var.AddPrimitiveTimeSource(alpha, tn, 0, ProdUh);
  for (int i = 0; i < Urand.GetM(); i++)
    {
      Real_wp val = Urand(i) + alpha*pulse*b_src(i);
      if (isnan(ProdUh(i) - val) || (abs(ProdUh(i) - val) > 1e-10))
	{
	  cout << "AddPrimitiveTimeSource incorrect" << endl;
	  abort();
	}
    }
    
  if (!first_order)
    {
      Urand.Resize(nb_dof_u);
      ProdUh.Resize(nb_dof_u);
      Y.Reallocate(nb_dof_u);
      var.var_harmonic.ImposeNullDirichletCondition(Urand);
      
      // testing EvaluateFunctionS
      VectReal_wp yp;
      var.EvaluateFunctionS(tn, Urand, yp, ProdUh, false, false);

      Mlt(Kh, Urand, Y, false);
      var.var_harmonic.ImposeNullDirichletCondition(Y);
      for (int i = 0; i < nb_dof_u; i++)
	if (isnan(abs(Y(i) + ProdUh(i))) || (abs(Y(i) + ProdUh(i)) > 1e-10))
	  {
	    cout << "EvaluateFunctionS incorrect" << endl;
	    DISP(i); DISP(Y(i)); DISP(ProdUh(i));
	    abort();
	  }

      yp.Reallocate(nb_dof_u);
      GenerateRandomVectorPb(var.var_harmonic, yp, true);
      
      var.var_harmonic.ImposeNullDirichletCondition(yp);
      var.EvaluateFunctionS(tn, Urand, yp, ProdUh, true, false);
      
      Mlt(Kh, Urand, Y, false);
      MltAdd(Real_wp(1), Sh, yp, Real_wp(1), Y, false);
      var.var_harmonic.ImposeNullDirichletCondition(Y);
      var.SolveOperatorDh(Y);
      for (int i = 0; i < nb_dof_u; i++)
	if (isnan(abs(Y(i) + ProdUh(i))) || (abs(Y(i) + ProdUh(i)) > 1e-10))
	  {
	    cout << "EvaluateFunctionS incorrect" << endl;
	    DISP(i); DISP(Y(i)); DISP(ProdUh(i));
	    abort();
	  }
	    
      var.EvaluateFunctionS(tn, Urand, yp, ProdUh, true, true);

      var.SourceOnlyTime(tn, 1, pulse);
      Mlt(Kh, Urand, Y, false);
      MltAdd(Real_wp(1), Sh, yp, Real_wp(1), Y, false);
      Y -= pulse*b_src;
      var.var_harmonic.ImposeNullDirichletCondition(Y);
      var.SolveOperatorDh(Y);
      for (int i = 0; i < nb_dof_u; i++)
	if (isnan(abs(Y(i) + ProdUh(i))) || (abs(Y(i) + ProdUh(i)) > 1e-10))
	  {
	    cout << "EvaluateFunctionS incorrect" << endl;
	    DISP(i); DISP(Y(i)); DISP(ProdUh(i));
	    abort();
	  }
    }
  else
    {	
      Ones.Reallocate(nb_dof_u+nb_dof_v);
      ProdUh.Reallocate(nb_dof_u+nb_dof_v);
      ProdVh.Reallocate(nb_dof_u+nb_dof_v);
      Y.Reallocate(nb_dof_u+nb_dof_v);
      ProdUh.Zero();   ProdVh.Zero();
      
      cout << "Testing ApplyMassMatrix" << endl;
      
      if (fast_check)
	{
	  if (nb_proc > 1)
	    GenerateRandomVectorPb(var.var_harmonic, Ones, true, true);
	  else
	    GenerateRandomVector(Ones, nb_dof_u+nb_dof_v);
	  
	  var.ApplyMassMatrix(Real_wp(1), zero, Ones, zero, ProdUh);

	  Mlt(Dh, Ones, Y, false);
	  if (!EqualVector(Y, ProdUh))
	    {
	      cout << "ApplyMassMatrix incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_u+nb_dof_v; i++)
	  {
	    Ones.Zero();
	    Ones(i) = 1.0;
	    
	    var.ApplyMassMatrix(Real_wp(1), zero, Ones, zero, ProdUh);
	    for (int j = 0; j < nb_dof_u+nb_dof_v; j++)
	      {
		Real_wp val = Dh(j, i);
		if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		  {
		    cout << "ApplyMassMatrix incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }
	    
	    ProdUh = Urand;
	    var.ApplyMassMatrix(alpha, zero, Ones, beta, ProdUh);
	    for (int j = 0; j < nb_dof_u+nb_dof_v; j++)
	      {
		Real_wp val = beta*Urand(j) + alpha*Dh(j, i);
		if (isnan(abs(ProdUh(j) - val)) || (abs(ProdUh(j) - val) > 1e-10))
		  {
		    cout << "ApplyMassMatrix incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      } 
	  }
	    
      cout << "Testing SolveMassMatrix" << endl;
      
      if (fast_check)
	{
	  var.ApplyMassMatrix(Real_wp(1), zero, Ones, zero, ProdUh);
	  var.SolveMassMatrix(ProdUh);

	  if (!EqualVector(Ones, ProdUh))
	    {
	      cout << "SolveMassMatrix incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_u+nb_dof_v; i++)
	  {
	    Ones.Zero();
	    Ones(i) = 1.0;
	    var.ApplyMassMatrix(Real_wp(1), zero, Ones, zero, ProdUh);
	    var.SolveMassMatrix(ProdUh);
	    for (int j = 0; j < nb_dof_u+nb_dof_v; j++)
	      {
		if (isnan(abs(ProdUh(j) - Ones(j))) || (abs(ProdUh(j) - Ones(j)) > 1e-10))
		  {
		    cout << "SolveMassMatrix incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(ProdVh(j));
		    abort();
		  }
	      } 
	  }
	    
      // testing Ch operator
      if (var.var_harmonic.FormulationDG() == ElementReference_Base::HDG)
        {
          Ones.Reallocate(nb_dof_L);
          ProdUh.Reallocate(nb_dof_L);
          if (fast_check)
	    {
	      Ones.FillRand(); Ones *= 1e-9;

	      ProdUh = Ones;
	      var.Glob_mat_Dh->SolveOperatorCh(ProdUh);

	      Y.Reallocate(nb_dof_L);
	      Mlt(invCh, Ones, Y);

	      if (!EqualVector(Y, ProdUh))
		{
		  cout << "SolveOperatorCh incorrect" << endl;
		  abort();
		}
	    }
	  else
	    for (int i = 0; i < nb_dof_L; i++)
	      {
		Ones.Zero();
		Ones(i) = 1.0;
		
		ProdUh = Ones;
		var.Glob_mat_Dh->SolveOperatorCh(ProdUh);
		
		for (int j = 0; j < nb_dof_L; j++)
		  {
		    Real_wp val = invCh(j, i);
		    if (isnan(abs(ProdUh(j) - val))
			|| (abs(ProdUh(j) - val) > 1e-10))
		      {
			cout << "SolveOperatorCh incorrect" << endl;
			DISP(i); DISP(j); DISP(val); DISP(ProdUh(j));
			abort();
		      }
		  }
	      }
          
          Ones.Reallocate(nb_dof_u+nb_dof_v);
          ProdUh.Reallocate(nb_dof_u+nb_dof_v);
	  Y.Reallocate(nb_dof_u+nb_dof_v);
        }
      
      cout << "Testing EvaluateFunction" << endl;
      if (fast_check)
	{
	  if (nb_proc > 1)
	    GenerateRandomVectorPb(var.var_harmonic, Ones, true, true);
	  else
	    GenerateRandomVector(Ones, nb_dof_u+nb_dof_v);
	  
	  var.EvaluateDerivativeFunction(tn, 0, Ones, ProdUh, false, false);
	  var.EvaluateFunction(tn, Ones, ProdVh, false, false);
	  
	  Mlt(Kh, Ones, Y, false);
	  MltAdd(one, Sh, Ones, one, Y, false);
	  Y = -Y;

	  if (!EqualVector(Y, ProdUh))
	    {
	      cout << "EvaluateFunction incorrect" << endl;
	      abort();
	    }

	  if (!EqualVector(Y, ProdVh))
	    {
	      cout << "EvaluateFunction incorrect" << endl;
	      abort();
	    }
	}
      else
	for (int i = 0; i < nb_dof_u+nb_dof_v; i++)
	  {
	    Ones.Zero();
	    Ones(i) = 1.0;
	    
	    var.EvaluateDerivativeFunction(tn, 0, Ones, ProdUh, false, false);
	    var.EvaluateFunction(tn, Ones, ProdVh, false, false);
	    
	    for (int j = 0; j < nb_dof_u+nb_dof_v; j++)
	      {
		Real_wp val = Kh(j, i) + Sh(j, i);
		if (isnan(abs(ProdUh(j) + val + ProdVh(j))) ||
		    (abs(ProdUh(j) + val) > 1e-10) || (abs(ProdVh(j) + val) > 1e-10) )
		  {
		    cout << "EvaluateFunction incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(-val);
		    abort();
		  }
	      }
	  }
      
      var.EvaluateFunction(tn, Urand, ProdUh, true, false);
      Mlt(Kh, Urand, ProdVh, false);
      MltAdd(one, Sh, Urand, one, ProdVh, false);
      var.var_harmonic.ImposeNullDirichletCondition(ProdVh);
      var.SolveMassMatrix(ProdVh);

      for (int j = 0; j < nb_dof_u+nb_dof_v; j++)
	if (isnan(abs(ProdUh(j) + ProdVh(j))) ||
	    (abs(ProdUh(j) + ProdVh(j)) > 1e-10))
	  {
	    cout << "EvaluateFunction incorrect" << endl;
	    DISP(j); DISP(ProdUh(j)); DISP(-ProdVh(j));
	    abort();
	  }

      var.EvaluateFunction(tn, Urand, ProdUh, true, true);
      Mlt(Kh, Urand, ProdVh, false);
      MltAdd(Real_wp(1), Sh, Urand, Real_wp(1), ProdVh, false);
      
      var.SourceOnlyTime(tn, 0, pulse);
      ProdVh -= pulse*b_src;
      var.var_harmonic.ImposeNullDirichletCondition(ProdVh);
      var.SolveMassMatrix(ProdVh);
      for (int j = 0; j < nb_dof_u+nb_dof_v; j++)
	if (isnan(abs(ProdUh(j) + ProdVh(j))) ||
	    (abs(ProdUh(j) + ProdVh(j)) > 1e-10))
	  {
	    cout << "EvaluateFunction incorrect" << endl;
	    DISP(j); DISP(ProdUh(j)); DISP(-ProdVh(j));
	    abort();
	  }
      
    }

  cout << "Testing OperatorDhPlusGammaKh" << endl;

  var.FactorizeOperatorDhPlusGammaKh(Real_wp(1), beta, gamma);
  var.glob_solver->SetPrintLevel(-1);

  Y = Urand;
  
  MltAdd(Real_wp(1), Dh, Y, Real_wp(0), ProdUh, false);
  MltAdd(beta, Sh, Y, Real_wp(1), ProdUh, false);
  MltAdd(gamma, Kh, Y, Real_wp(1), ProdUh, false);

  var.SolveOperatorDhPlusGammaKh(zero, ProdUh, Y);
  
  for (int j = 0; j < Y.GetM(); j++)
    {
      if (isnan(abs(Y(j) - Urand(j))) || (abs(Y(j) - Urand(j)) > 1e-10))
	{
	  cout << "SolveOperatorDhPlusGammaKh incorrect" << endl;
	  DISP(j); DISP(ProdUh(j)); DISP(Urand(j));
	  abort();
	}
    }
    
  // testing FactorizeOperatorReal

  Ones.Reallocate(2);
  Ones.Fill(1);
  VectReal_wp beta_real(2), gamma_real(2);
  beta_real.FillRand(); beta_real *= 1e-10;
  gamma_real.FillRand(); gamma_real *= 1e-10;

  DISP(beta_real); DISP(gamma_real);
  var.FactorizeOperatorReal(Ones, beta_real, gamma_real);

  for (int num = 0; num < 2; num++)
    {
      MltAdd(Real_wp(1), Dh, Urand, Real_wp(0), Y, false);
      MltAdd(beta_real(num), Sh, Urand, Real_wp(1), Y, false);
      MltAdd(gamma_real(num), Kh, Urand, Real_wp(1), Y, false);

      var.SolveOperatorReal(zero, Y, ProdUh, num);
	    
      for (int j = 0; j < Y.GetM(); j++)
	{
	  if (isnan(abs(ProdUh(j) - Urand(j))) || (abs(ProdUh(j) - Urand(j)) > 1e-10))
	    {
	      cout << "SolveOperatorReal incorrect" << endl;
	      DISP(j); DISP(ProdUh(j)); DISP(Urand(j));
	      abort();
	    }
	}
    }

  return;
  
  // testing FactorizeOperatorComplex

  VectComplex_wp Ones_cplx(2);
  Ones_cplx.Reallocate(2);
  Ones_cplx.Fill(1);
  VectComplex_wp beta_cplx(2), gamma_cplx(2);
  GenerateRandomVector(beta_cplx, 2);
  GenerateRandomVector(gamma_cplx, 2);

  DISP(beta_cplx); DISP(gamma_cplx);
  var.FactorizeOperatorComplex(Ones_cplx, beta_cplx, gamma_cplx);

  VectComplex_wp Urand_cplx, Ycplx, ProdU_cplx;
  if (nb_proc > 1)
    GenerateRandomVectorPb(var.var_harmonic, Urand_cplx, true, true);
  else
    GenerateRandomVector(Urand_cplx, nb_dof_u+nb_dof_v);
  
  Ycplx.Reallocate(nb_dof_u+nb_dof_v);
  ProdU_cplx.Reallocate(nb_dof_u+nb_dof_v);
  for (int num = 0; num < 2; num++)
    {
      MltAdd(Complex_wp(1, 0), Dh, Urand_cplx, Complex_wp(0, 0), Ycplx, false);
      MltAdd(beta_cplx(num), Sh, Urand_cplx, Complex_wp(1, 0), Ycplx, false);
      MltAdd(gamma_cplx(num), Kh, Urand_cplx, Complex_wp(1, 0), Ycplx, false);

      var.SolveOperatorComplex(zero, Ycplx, ProdU_cplx, num);
	    
      for (int j = 0; j < Y.GetM(); j++)
	{
	  if (isnan(abs(ProdU_cplx(j) - Urand_cplx(j))) || (abs(ProdU_cplx(j) - Urand_cplx(j)) > 1e-10))
	    {
	      cout << "SolveOperatorComplex incorrect" << endl;
	      DISP(j); DISP(ProdU_cplx(j)); DISP(Urand_cplx(j));
	      abort();
	    }
	}
    }
      
}


void CheckMatrixProductLevel()
{
  int N = 100;
  int nb_elt_pml = rand()%(N/2);
  DISP(N); DISP(nb_elt_pml);
  
  MatrixVectorProductLevel lvl;
  lvl.SetNbElt(N, nb_elt_pml);

  lvl.SetLevel(MatrixVectorProductLevel::ALL_LEVELS);

  if (lvl.GetNbElt() != N)
    {
      cout << "SetNbElt / GetNbElt incorrect" << endl;
      abort();
    }

  int nb = 0;
  for (int i = 0; i < N; i++)
    if (lvl.TreatElement(i))
      nb++;
  
  if (nb != N)
    {
      cout << "TreatElement incorrect" << endl;
      abort();
    }

  for (int i = 0; i < N; i++)
    if (lvl.GetElementNumber(i) != i)
      {
	cout << "GetElementNumber incorrect" << endl;
	abort();
      }

  lvl.SetLevel(MatrixVectorProductLevel::LVL_NOPML);
  
  if (lvl.GetNbElt() != N-nb_elt_pml)
    {
      cout << "SetLevel incorrect" << endl;
      abort();
    }

  nb = 0;
  bool check_treat = true;
  for (int i = 0; i < N; i++)
    if (lvl.TreatElement(i))
      {
	nb++;
	if (i >= N-nb_elt_pml)
	  check_treat = false;
      }
  
  if ((nb != N-nb_elt_pml) || (!check_treat))
    {
      cout << "TreatElement incorrect" << endl;
      abort();
    }

  for (int i = 0; i < lvl.GetNbElt(); i++)
    if (lvl.GetElementNumber(i) != i)
      {
	cout << "GetElementNumber incorrect" << endl;
	abort();
      }

  lvl.SetLevel(MatrixVectorProductLevel::LVL_PML);
  
  if (lvl.GetNbElt() != nb_elt_pml)
    {
      cout << "SetLevel incorrect" << endl;
      abort();
    }

  nb = 0;
  check_treat = true;
  for (int i = 0; i < N; i++)
    if (lvl.TreatElement(i))
      {
	nb++;
	if (i < N-nb_elt_pml)
	  check_treat = false;
      }
  
  if ((nb != nb_elt_pml) || (!check_treat))
    {
      cout << "TreatElement incorrect" << endl;
      abort();
    }

  int offset = N-nb_elt_pml;
  for (int i = 0; i < lvl.GetNbElt(); i++)
    if (lvl.GetElementNumber(i) != offset + i)
      {
	cout << "GetElementNumber incorrect" << endl;
	abort();
      }

  // we test with three levels
  int nb_levels = 3;
  Vector<int> IndexElt(N);
  for (int i = 0; i < N; i++)
    IndexElt(i) = rand()%nb_levels;

  Vector<IVect> num(nb_levels);
  Vector<int> nb_elt_per_level(nb_levels);
  nb_elt_per_level.Zero();
  for (int i = 0; i < N; i++)
    nb_elt_per_level(IndexElt(i))++;

  for (int i = 0; i < nb_levels; i++)
    num(i).Reallocate(nb_elt_per_level(i));

  nb_elt_per_level.Zero();
  for (int i = 0; i < N; i++)
    {
      num(IndexElt(i))(nb_elt_per_level(IndexElt(i))) = i;
      nb_elt_per_level(IndexElt(i))++;
    }

  for (int i = 0; i < nb_levels; i++)
    {
      DISP(i); DISP(num(i));
    }

  lvl.SetLevelArray(num);
  
  for (int i = 0; i < nb_levels; i++)
    {
      lvl.SetLevel(i);
      if (lvl.GetNbElt() != nb_elt_per_level(i))
	{
	  cout << "SetLevel/GetLevelArray incorrect" << endl;
	  abort();
	}

      for (int k = 0; k < lvl.GetNbElt(); k++)
	if (lvl.GetElementNumber(k) != num(i)(k))
	  {
	    cout << "SetLevel/GetLevelArray incorrect" << endl;
	    abort();
	  }

      nb = 0;
      for (int k = 0; k < N; k++)
	if (lvl.TreatElement(k))
	  {
	    if (num(i)(nb) != k)
	      {
		cout << "SetLevel/GetLevelArray incorrect" << endl;
		abort();		
	      }
	    
	    nb++;
	  }
      
      if (nb != nb_elt_per_level(i))
	{
	  cout << "SetLevel/GetLevelArray incorrect" << endl;
	  abort();
	}
    }

  cout << "Memory used by lvl = " << lvl.GetMemorySize() << endl;
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e5*epsilon_machine;
  if (argc > 1)
    fast_check = false;

  
  /*{
    // hexahedral mesh (distorted)
    HyperbolicProblem<TimeMaxwellEquationHdg_3D> var;
    string input_file("src/Program/Unit/Instationary/cube_hexa.ini");
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_HCURL_LOBATTO", input_file, false, false);
  }

  {
    // hexahedral mesh (cubes)
    HyperbolicProblem<TimeMaxwellEquationHdg_3D> var;
    string input_file("src/Program/Unit/Instationary/cube_no_pml.ini");
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_HCURL_LOBATTO", input_file, false, false);
  }

  return FinalizeMontjoie();*/
    
  CheckMatrixProductLevel();
  
  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/time_disque.ini");
    CheckProblem(var, true, "ACOUSTIC_HDG", "TRIANGLE_LOBATTO", input_file, false);
  }
  
  
  {
    // hybrid mesh
    HyperbolicProblem<AcousticEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, false, "ACOUSTIC", "TETRAHEDRON_LOBATTO", input_file, true);
    CheckProblem(var, true, "ACOUSTIC", "TETRAHEDRON_LOBATTO", input_file, true);
  }
   
  {
    // case with PML
    HyperbolicProblem<AcousticEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube.ini");
    CheckProblem(var, true, "ACOUSTIC", "TETRAHEDRON_LOBATTO", input_file, true);
    CheckProblem(var, true, "ACOUSTIC", "TETRAHEDRON_LOBATTO", input_file, true);
    
  }
  
  {
    // hexahedral mesh
    HyperbolicProblem<AcousticEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_hexa.ini");
    CheckProblem(var, false, "ACOUSTIC", "TETRAHEDRON_LOBATTO", input_file, true);
    CheckProblem(var, true, "ACOUSTIC", "TETRAHEDRON_LOBATTO", input_file, true);
  }

  /* Testing Maxwell's equations */
  {
    // hybrid mesh
    HyperbolicProblem<TimeMaxwellEquation_3D> var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, false, "TIME_MAXWELL", "HEXAHEDRON_FIRST_FAMILY", input_file, true);
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_FIRST_FAMILY", input_file, true);
  }
  
  {
    // case with PML
    HyperbolicProblem<TimeMaxwellEquation_3D> var;
    string input_file("src/Program/Unit/Instationary/cube.ini");
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_HCURL_LOBATTO", input_file, true);
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_FIRST_FAMILY", input_file, true);
    
  }
  
  {
    // hexahedral mesh
    HyperbolicProblem<TimeMaxwellEquation_3D> var;
    string input_file("src/Program/Unit/Instationary/cube_hexa.ini");
    CheckProblem(var, false, "TIME_MAXWELL", "HEXAHEDRON_HCURL_LOBATTO", input_file, true);
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_HCURL_LOBATTO", input_file, true);
    CheckProblem(var, false, "TIME_MAXWELL", "HEXAHEDRON_FIRST_FAMILY", input_file, true);
    CheckProblem(var, true, "TIME_MAXWELL", "HEXAHEDRON_FIRST_FAMILY", input_file, true);    
  }

  {
    HyperbolicProblem<AcousticEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre.ini");
    CheckProblem(var, false, "ACOUSTIC", "TRIANGLE_LOBATTO", input_file);
  }

  {
    HyperbolicProblem<AcousticEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre.ini");
    CheckProblem(var, true, "ACOUSTIC", "TRIANGLE_LOBATTO", input_file);
  }

  {
    HyperbolicProblem<AcousticEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/disque.ini");
    CheckProblem(var, false, "ACOUSTIC", "TRIANGLE_LOBATTO", input_file);
  }

  {
    HyperbolicProblem<AcousticEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/disque.ini");
    CheckProblem(var, true, "ACOUSTIC", "TRIANGLE_LOBATTO", input_file);
  }
    
  {
    HyperbolicProblem<AcousticEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_pml_no_damp.ini");
    CheckProblem(var, true, "ACOUSTIC", "TRIANGLE_LOBATTO", input_file);
  }

  return FinalizeMontjoie();

  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg.ini");
    CheckProblem(var, true, "ACOUSTIC_DG", "TRIANGLE_LOBATTO", input_file, false);
  }

  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/disque_dg.ini");
    CheckProblem(var, true, "ACOUSTIC_DG", "TRIANGLE_LOBATTO", input_file, false);
  }
    
  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, true, "ACOUSTIC_DG", "TRIANGLE_LOBATTO", input_file, false);
  }

  {
    HyperbolicProblem<AcousticEquationDG<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, true, "ACOUSTIC_DG", "TETRAHEDRON_LOBATTO", input_file, false);
  }
  
  /* Testing AcousticEquation */

  {
    HyperbolicProblem<AcousticEquationDG<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, true, "ACOUSTIC_HDG", "TETRAHEDRON_LOBATTO", input_file, false);
  }

  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, true, "ACOUSTIC_HDG", "TRIANGLE_LOBATTO", input_file, false);
  }
  
  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg.ini");
    CheckProblem(var, true, "ACOUSTIC_HDG", "TRIANGLE_LOBATTO", input_file, false);
  }

  {
    HyperbolicProblem<AcousticEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/disque_dg.ini");
    CheckProblem(var, true, "ACOUSTIC_HDG", "TRIANGLE_LOBATTO", input_file, false);
  }

  cout << "All tests passed successfully" << endl;

  return FinalizeMontjoie();
}
