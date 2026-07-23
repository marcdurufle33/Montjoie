#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG
#define MONTJOIE_WITH_LEGENDRE_DG

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_HCURL

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"

using namespace Montjoie;

Real_wp threshold;
bool fast_check;

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
void GetRandNumber(T& x, bool reel = false)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x, bool reel = false)
{
  if (reel)
    x = complex<T>(rand(), 0)/RAND_MAX;
  else
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
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  Real_wp yref = y.GetNormInf();
  for (int i = 0; i < x.GetM(); i++)
    if (isnan(abs(x(i) - y(i)) / yref) || abs(x(i) - y(i))/yref > eps)
      {
	DISP(i); DISP(x(i)); DISP(y(i));
	return false;
      }
  
  return true;
}


template<class T, class TypeEquation>
void CheckProblem(const EllipticProblem<TypeEquation>& var_c, const T& z,
		  bool first_order, const string& name_equation,
		  const string& name_element, const string& input_file, bool force_sym,
		  bool check_trans = true, bool check_condensed = false)
{
  EllipticProblem<TypeEquation> var;
  
  var.SetTypeEquation(name_equation);
    
  var.InitIndices(100);
  ReadInputFile(input_file, var);
  
  var.SetFirstOrderFormulation(first_order);

  var.ComputeMeshAndFiniteElement(name_element);
    
  int nb_proc = var.GetNbProcPerMode();
  int rank_proc = var.GetRankProcMode();
  string suffix = to_str(rank_proc) + ".dat";
  if (nb_proc == 1)
    var.mesh.Write("test.mesh");
    
  var.PerformOtherInitializations();
    
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();

  if (check_condensed)
    {
      var.SetCoefficientDirichlet(Real_wp(1));
	
      All_LinearSolver solver(var);
      ReadInputFile(input_file, solver);

      solver.EnableMatrixStorage(true, "mat_glob.dat");
      solver.EnableStaticCondensation(false, false, false);
					    
      GlobalGenericMatrix<T> nat_mat;
      
      T alpha; GetRandNumber(alpha);
      T beta; GetRandNumber(beta);
      T gamma; GetRandNumber(gamma);
      //SetComplexOne(beta); SetComplexOne(gamma); SetComplexOne(alpha);
      
      //var.mesh_num.GlobDofNumber_Subdomain.WriteText("dof" + suffix);
	
      //SetComplexOne(alpha); //SetComplexOne(beta); SetComplexOne(gamma);
      //gamma = beta;

      nat_mat.SetCoefMass(alpha);
      nat_mat.SetCoefDamping(beta);
      nat_mat.SetCoefStiffness(gamma);

      solver.PerformFactorizationStep(nat_mat);
      srand(Seed);
      
      Vector<T> x_sol(var.GetNbDof()), b_source;
            
      Vector<T> x_sol_glob;
      GenerateRandomVector(x_sol_glob, max(var.GetNbDof(), var.GetNbGlobalDof()));

      Vector<int>& row_num = var.mesh_num.GlobDofNumber_Subdomain;

      Vector<Real_wp> coef_distrib(var.GetNbDof());
      coef_distrib.Fill(1.0);
      var.AddDomains(coef_distrib);
      for (int i = 0; i < x_sol.GetM(); i++)
	coef_distrib(i) = 1.0/coef_distrib(i);
      
      if (nb_proc == 1)
	x_sol = x_sol_glob;
      else
	for (int i = 0; i < x_sol.GetM(); i++)
	  x_sol(i) = coef_distrib(i)*x_sol_glob(row_num(i));
      
      //for (int i = 120; i < x_sol.GetM(); i++)
      //SetComplexZero(x_sol(i));
      
      x_sol.Write("source" + suffix);
      
      b_source = x_sol;
      solver.ComputeSolution(x_sol);

      x_sol.Write("sol" + suffix);

      solver.EnableMatrixStorage(true, "mat_cond.dat");
      solver.EnableStaticCondensation(true, true, true);

      solver.ClearFactorization();
      solver.PerformFactorizationStep(nat_mat);
      srand(Seed);
      
      solver.ComputeSolution(b_source);

      //b_source.Write("sol_cond.dat");
      if (!EqualVector(x_sol, b_source))
	{
	  cout << "Static condensation not working" << endl;
	  abort();
	}
    }

  var.SetSymmetrizationUse(force_sym);

  DistributedMatrix<T, General, ArrayRowSparse> Dh, Sh, Kh;
  GlobalGenericMatrix<T> nat_mat;
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(0.0);
  var.SetCoefficientDirichlet(Real_wp(1));
  var.AddMatrixWithBC(Dh, nat_mat);

  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(1.0);
  var.SetCoefficientDirichlet(Real_wp(0));
  var.AddMatrixWithBC(Sh, nat_mat);

  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(1.0);
  nat_mat.SetCoefDamping(0.0);
  var.SetCoefficientDirichlet(Real_wp(0));
  var.AddMatrixWithBC(Kh, nat_mat);
  Dh.WriteText("Dh.dat");
  Sh.WriteText("Sh.dat");
  Kh.WriteText("Kh.dat");
  
  T one; SetComplexOne(one);
  T alpha; GetRandNumber(alpha, true);
  T beta; GetRandNumber(beta, true);
  T gamma; GetRandNumber(gamma, true);
  //SetComplexZero(beta); SetComplexZero(gamma); SetComplexOne(alpha);
  //SetComplexZero(beta); SetComplexOne(gamma); SetComplexZero(alpha);
  //SetComplexOne(beta); SetComplexOne(gamma); SetComplexOne(alpha);
  
  int N = var.GetNbDof();
  DISP(N); DISP(var.mesh_num.GetNbDof());
  DISP(var.GetNbMainUnknownDof());
  
  Vector<T> Ones(N), ProdUh(N);
  Vector<T> Y(N), Yref(N);
  
  FemMatrixFreeClass<T, TypeEquation> Ah(var);

  nat_mat.SetCoefMass(alpha);
  nat_mat.SetCoefDamping(beta);
  nat_mat.SetCoefStiffness(gamma);
  
  var.AddMatrixWithBC(Ah, nat_mat);

  // testing matrix-free vector product (fast procedure with only one random vector)
  if (fast_check)
    {
      Vector<T> Ones_glob;
      GenerateRandomVector(Ones_glob, max(N, var.GetNbGlobalDof()));

      Vector<int>& row_num = var.mesh_num.GlobDofNumber_Subdomain;

      Vector<Real_wp> coef_distrib(var.GetNbDof());
      coef_distrib.Fill(1.0);
      var.AddDomains(coef_distrib);
      for (int i = 0; i < N; i++)
	coef_distrib(i) = 1.0/coef_distrib(i);
      
      if (nb_proc == 1)
	Ones = Ones_glob;
      else
	for (int i = 0; i < Ones.GetM(); i++)
	  Ones(i) = coef_distrib(i)*Ones_glob(row_num(i));

      var.ImposeNullDirichletCondition(Ones);
      
      Yref.Zero();
      MltAdd(alpha, Dh, Ones, one, Yref);
      MltAdd(beta, Sh, Ones, one, Yref);
      MltAdd(gamma, Kh, Ones, one, Yref);
      var.ImposeNullDirichletCondition(Yref);
 
      Y.Zero();
      Ah.MltVector(Ones, Y);

      if (!EqualVector(Yref, Y))
	{
	  cout << "MltVector incorrect" << endl;
	  DISP(Norm2(Y)); DISP(Norm2(Yref));
	  abort();
	}

      if (check_trans)
	{
	  Yref.Zero();
	  MltAdd(alpha, SeldonTrans, Dh, Ones, one, Yref);
	  MltAdd(beta, SeldonTrans, Sh, Ones, one, Yref);
	  MltAdd(gamma, SeldonTrans, Kh, Ones, one, Yref);
	  var.ImposeNullDirichletCondition(Yref);
	  
	  Y.FillRand();
	  Ah.MltVector(SeldonTrans, Ones, Y);
	  
	  if (!EqualVector(Yref, Y))
	    {
	      cout << "MltVector (transpose) incorrect" << endl;
	      DISP(Norm2(Y)); DISP(Norm2(Yref));
	      abort();
	    }

	  Yref.Zero();
	  MltAdd(conjugate(alpha), SeldonConjTrans, Dh, Ones, one, Yref);
	  MltAdd(conjugate(beta), SeldonConjTrans, Sh, Ones, one, Yref);
	  MltAdd(conjugate(gamma), SeldonConjTrans, Kh, Ones, one, Yref);
	  var.ImposeNullDirichletCondition(Yref);
	  
	  Y.FillRand();
	  Ah.MltVector(SeldonConjTrans, Ones, Y);
	  
	  if (!EqualVector(Yref, Y))
	    {
	      cout << "MltVector (conjugate transpose) incorrect" << endl;
	      DISP(Norm2(Y)); DISP(Norm2(Yref));
	      abort();
	    }
	}
      
      return;
    }

  if (nb_proc > 1)
    {
      cout << "Not implemented" << endl;
      abort();
    }
  
  // testing matrix-free vector product (slow procedure, all columns are tested)
  for (int i = 0; i < N; i++)
    if ((var.FormulationDG() != ElementReference_Base::CONTINUOUS)
	|| (!var.IsDofDirichlet(i)))
      
      {
	Ones.Fill(0); Ones(i) = 1.0;
	Ah.MltVector(Ones, ProdUh);
	
	for (int j = 0; j < N; j++)
	  {
	    T val = alpha*Dh(j, i) + beta*Sh(j, i) + gamma*Kh(j, i);
	    if ((var.FormulationDG() == ElementReference_Base::CONTINUOUS)
		&& var.IsDofDirichlet(j))
	      SetComplexZero(val);
	    
	    if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
	      {
		cout << "MltVector incorrect" << endl;
		DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		abort();
	      }
	  }
	
	if (check_trans)
	  {
	    Ah.MltVector(SeldonTrans, Ones, ProdUh);
	    
	    for (int j = 0; j < N; j++)
	      {
	        T val = alpha*Dh(i, j) + beta*Sh(i, j) + gamma*Kh(i, j);
		if ((var.FormulationDG() == ElementReference_Base::CONTINUOUS)
		    && var.IsDofDirichlet(j))
		  SetComplexZero(val);
		
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "MltVector (transpose) incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val);
		    abort();
		  }
	      }

	    Ah.MltVector(SeldonConjTrans, Ones, ProdUh);
	    
	    for (int j = 0; j < N; j++)
	      {
		T val = alpha*conjugate(Dh(i, j)) + beta*conjugate(Sh(i, j)) + gamma*conjugate(Kh(i, j));
		if ((var.FormulationDG() == ElementReference_Base::CONTINUOUS)
		    && var.IsDofDirichlet(j))
		  SetComplexZero(val);
		
		if (isnan(abs(ProdUh(j) - val)) || abs(ProdUh(j) - val) > 1e-10)
		  {
		    cout << "MltVector (transpose conjugate) incorrect" << endl;
		    DISP(i); DISP(j); DISP(ProdUh(j)); DISP(val); DISP(alpha); DISP(beta); DISP(gamma);
		    abort();
		  }
	      }
	  }
      }
}


template<class TypeEquation>
void CheckEigenProblem(const EllipticProblem<TypeEquation>& var_c, bool first_order, const string& name_equation,
		       const string& name_element, const string& input_file)
{
  EllipticProblem<TypeEquation> var;
    
  var.SetTypeEquation(name_equation);
    
  var.InitIndices(100);
  ReadInputFile(input_file, var);
  
  var.SetFirstOrderFormulation(first_order);
  var.SetSymmetrizationUse(false);
  var.EnableSymmetricDirichlet();

  int nb_proc = var.GetNbProcPerMode();
  if (nb_proc > 1)
    return;
  
  var.ComputeMeshAndFiniteElement(name_element);
  
  var.mesh.Write("test.mesh");
  var.PerformOtherInitializations();
  
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();

  DistributedMatrix<Real_wp, General, ArrayRowSparse> Dh, Kh;
  GlobalGenericMatrix<Real_wp> nat_mat;

  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefStiffness(0.0);
  nat_mat.SetCoefDamping(0.0);
  var.SetCoefficientDirichlet(Real_wp(1));
  var.AddMatrixWithBC(Dh, nat_mat);

  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(1.0);
  nat_mat.SetCoefDamping(1.0);
  Real_wp coef_dir(2.5);
  var.SetCoefficientDirichlet(coef_dir);
  var.AddMatrixWithBC(Kh, nat_mat);
  
  //Real_wp one(1), zero(0);
  Real_wp alpha; GetRandNumber(alpha);
  Real_wp beta; GetRandNumber(beta);
  Real_wp gamma; GetRandNumber(gamma);
  //SetComplexZero(beta); SetComplexZero(gamma); SetComplexOne(alpha);
  //SetComplexZero(beta); SetComplexOne(gamma); SetComplexZero(alpha);
  
  int N = var.GetNbDof();
  DISP(N); DISP(var.mesh_num.GetNbDof());
  DISP(var.GetNbMainUnknownDof()); DISP(alpha); DISP(beta);

  All_LinearSolver solver(var);
  EigenProblemMontjoie<Real_wp> var_eig(var, solver);
  //var_eig.SetDirichletCoef(coef_dir);
  
  FemMatrixFreeClass_Base<Real_wp>* Mh_free, *Kh_free;

  var.SetCoefficientDirichlet(Real_wp(1));
  Mh_free = var.GetNewIterativeMatrix(Real_wp(0));
  var.SetCoefficientDirichlet(coef_dir);
  Kh_free = var.GetNewIterativeMatrix(Real_wp(0));

  Mh_free->SetDirichletCondition(Dh, 0, 0);
  Kh_free->SetDirichletCondition(Kh, 0, 0);
  
  Dh.WriteText("Dh.dat");
  Kh.WriteText("Kh.dat");
  
  var_eig.InitMatrix(*Kh_free, *Mh_free);

  VectReal_wp X(N), Y(N), Yref(N);
  
  if (var_eig.DiagonalMass())
    {
      var_eig.ComputeDiagonalMass();
      var_eig.FactorizeDiagonalMass();

      VectReal_wp diag;
      var_eig.GetSqrtDiagonal(diag);

      diag.WriteText("sqrt_diag.dat");
      for (int i = 0; i < Dh.GetM(); i++)
	{
	  if ((Dh.GetRowSize(i) != 1) || (Dh.Index(i, 0) != i))
	    {
	      cout << "Dh is not diagonal" << endl;
	      abort();
	    }

	  if ((abs(diag(i)*diag(i) - Dh.Value(i, 0)) > threshold)
	      || isnan(abs(diag(i)*diag(i) - Dh.Value(i, 0))))
	    {
	      DISP(i); DISP(diag(i)); DISP(Dh.Value(i, 0));
	      cout << "Incorrect diagonal" << endl;
	      abort();
	    }
	}
    }

  // cas cholesky a faire

  // produit avec la matrice de masse
  var_eig.ComputeMassMatrix();
  if (fast_check)
    {
      GenerateRandomVector(X, N);            
      var_eig.MltMass(X, Y);
      
      Mlt(Dh, X, Yref);
      if (!EqualVector(Y, Yref))
	{
	  cout << "MltMass incorrect" << endl;
	  abort();
	}
    }
  else
    {
      for (int i = 0; i < N; i++)
	{
	  X.Zero(); X(i) = 1.0;

	  var_eig.MltMass(X, Y);
	  for (int j = 0; j < N; j++)
	    {
	      Real_wp val = Dh(j, i);
	      if (isnan(abs(Y(j) - val)) || abs(Y(j) - val) > 1e-10)
		{
		  cout << "MltMass incorrect" << endl;
		  abort();
		}
	    }
	}
    }

  // produit avec la rigidite
  var_eig.ComputeStiffnessMatrix();
  if (fast_check)
    {
      GenerateRandomVector(X, N);            
      var_eig.MltStiffness(X, Y);
      
      Mlt(Kh, X, Yref);
      if (!EqualVector(Y, Yref))
	{
	  cout << "MltStiffness incorrect" << endl;
	  abort();
	}
    }
  else
    {
      for (int i = 0; i < N; i++)
	{
	  X.Zero(); X(i) = 1.0;

	  var_eig.MltStiffness(X, Y);
	  for (int j = 0; j < N; j++)
	    {
	      Real_wp val = Kh(j, i);
	      if (isnan(abs(Y(j) - val)) || abs(Y(j) - val) > 1e-10)
		{
		  cout << "MltStiffness incorrect" << endl;
		  abort();
		}
	    }
	}
    }

  var_eig.ComputeStiffnessMatrix(alpha, beta);
  if (fast_check)
    {
      GenerateRandomVector(X, N);            
      var_eig.MltStiffness(alpha, beta, X, Y);

      Yref.Zero();
      MltAdd(alpha, Dh, X, Real_wp(1), Yref);
      MltAdd(beta, Kh, X, Real_wp(1), Yref);
      if (!EqualVector(Y, Yref))
	{
	  cout << "MltStiffness(a,b) incorrect" << endl;
	  abort();
	}
    }
  else
    {
      for (int i = 0; i < N; i++)
	{
	  X.Zero(); X(i) = 1.0;

	  var_eig.MltStiffness(alpha, beta, X, Y);
	  for (int j = 0; j < N; j++)
	    {
	      Real_wp val = alpha*Dh(j, i) + beta*Kh(j, i);
	      if (isnan(abs(Y(j) - val)) || abs(Y(j) - val) > 1e-10)
		{
		  DISP(i); DISP(j); DISP(val); DISP(Y(j));
		  cout << "MltStiffness(a,b) incorrect" << endl;
		  abort();
		}
	    }
	}
    }

  // testing ComputeSolution
  GenerateRandomVector(X, N);

  Yref.Zero();
  MltAdd(alpha, Dh, X, Real_wp(1), Yref);
  MltAdd(beta, Kh, X, Real_wp(1), Yref);

  Y.Zero();
  var_eig.ComputeAndFactorizeStiffnessMatrix(alpha, beta);
  var_eig.ComputeSolution(Yref, Y);

  if (!EqualVector(X, Y))
    {
      cout << "ComputeSolution incorrect" << endl;
      abort();
    }

  // not implemented
  /* Yref.Zero();
  MltAdd(alpha, SeldonTrans, Dh, X, Real_wp(1), Yref);
  MltAdd(beta, SeldonTrans, Kh, X, Real_wp(1), Yref);

  Y.Zero();
  var_eig.ComputeSolution(SeldonTrans, Yref, Y);

  if (!EqualVector(X, Y))
    {
      cout << "ComputeSolution incorrect" << endl;
      abort();
      } */

  // to correct
  /* Complex_wp alpha_c, beta_c;
  GetRandNumber(alpha_c); GetRandNumber(beta_c); DISP(alpha_c); DISP(beta_c);

  var_eig.ComputeAndFactorizeStiffnessMatrix(alpha_c, beta_c, EigenProblem_Base<Real_wp>::REAL_PART);

  Y.Zero();
  var_eig.ComputeSolution(Yref, Y);

  if (!EqualVector(X, Y))
    {
      cout << "ComputeSolution incorrect" << endl;
      abort();
    }
  */
  
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 10000*epsilon_machine;
  
  fast_check = true;
  if (argc == 2)
    fast_check = false;

  /*  {
    EllipticProblem<HelmholtzEquation<Dimension2> > var;
    string input_file("carre_dielec.ini");
    string type_element, type_equation;
    getElement_Equation(input_file, type_element, type_equation);
    
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, true, true, true);
    
  }
  
  

  return FinalizeMontjoie();
  */
  
  // checking HelmholtzEquation with TRIANGLE_LOBATTO / TETRAHEDRON_LOBATTO
  
  {
    // case with affine elements
    EllipticProblem<HelmholtzEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, true, true, true);
  }
  
  {
    // case with curved elements
    EllipticProblem<HelmholtzEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/disque.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, true, true, true);
  }
  
  {
    // case with PMLs
    EllipticProblem<HelmholtzEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_pml.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }
  
  {
    // case with PMLs (without damping for testing first order formulation)
    EllipticProblem<HelmholtzEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_pml_no_damp.ini");
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }
    
  {
    // case with affine elements and PML
    EllipticProblem<HelmholtzEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
  }

  {
    // case with affine elements and PML
    EllipticProblem<HelmholtzEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/bicouche.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
  }

  {
    // case with curved elements (hybrid mesh)
    EllipticProblem<HelmholtzEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TETRAHEDRON_LOBATTO", input_file, true, true, true);
  }

      
  // checking HelmholtzEquation with TRIANGLE_HIERARCHIC / TETRAHEDRON_HIERARCHIC
  {
    EllipticProblem<HelmholtzEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TRIANGLE_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TRIANGLE_HIERARCHIC", input_file, true, true, true);
  }

  {
    // currently not implemented
    // case with PMLs
    //EllipticProblem<HelmholtzEquation<Dimension2> > var;
    //string input_file("src/Program/Unit/Instationary/carre_pml.ini");
    //CheckProblem(var, true, "HELMHOLTZ", "TRIANGLE_HIERARCHIC", input_file, false, false);
  }

  {
    // case with curved elements (hybrid mesh)
    EllipticProblem<HelmholtzEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ", "TETRAHEDRON_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TETRAHEDRON_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ", "TETRAHEDRON_HIERARCHIC", input_file, true, true, true);
  }

  // checking LaplaceEquation with TRIANGLE_LOBATTO / TETRAHEDRON_LOBATTO
  {
    // case with affine elements
    EllipticProblem<LaplaceEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_LOBATTO", input_file, true, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TRIANGLE_LOBATTO", input_file);
  }
  
  {
    // case with curved elements
    EllipticProblem<LaplaceEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/disque.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_LOBATTO", input_file, true, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TRIANGLE_LOBATTO", input_file);
  }

  {
    // case with PMLs
    EllipticProblem<LaplaceEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_pml.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }
  
  {
    // case with PMLs (without damping for testing first order formulation)
    EllipticProblem<LaplaceEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_pml_no_damp.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_LOBATTO", input_file, false, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TRIANGLE_LOBATTO", input_file);
  }
  
  {
    // case with affine elements and PML
    EllipticProblem<LaplaceEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file);
  }

  {
    // case with curved elements (hybrid mesh)
    EllipticProblem<LaplaceEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file, true, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TETRAHEDRON_LOBATTO", input_file);
  }

  // checking LaplaceEquation with TRIANGLE_HIERARCHIC / TETRAHEDRON_HIERARCHIC
  {
    EllipticProblem<LaplaceEquation<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TRIANGLE_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TRIANGLE_HIERARCHIC", input_file, true, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TRIANGLE_HIERARCHIC", input_file);
  }

  {
    // currently not implemented
    // case with PMLs
    //EllipticProblem<LaplaceEquation<Dimension2> > var;
    //string input_file("src/Program/Unit/Instationary/carre_pml.ini");
    //CheckProblem(var, true, "LAPLACE", "TRIANGLE_HIERARCHIC", input_file, false, false);
  }

  {
    // case with curved elements (hybrid mesh)
    EllipticProblem<LaplaceEquation<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE", "TETRAHEDRON_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TETRAHEDRON_HIERARCHIC", input_file, false, true, true);
    CheckProblem(var, Complex_wp(0), true, "LAPLACE", "TETRAHEDRON_HIERARCHIC", input_file, true, true, true);
    CheckEigenProblem(var, true, "LAPLACE", "TETRAHEDRON_HIERARCHIC", input_file);
  }
  
  // checking HelmholtzEquationDG with TRIANGLE_LOBATTO / TETRAHEDRON_LOBATTO
  {
    EllipticProblem<HelmholtzEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, Complex_wp(0), false, "HELMHOLTZ_HDG", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }

  {
    EllipticProblem<HelmholtzEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, Complex_wp(0), true, "HELMHOLTZ_HDG", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }

  // checking LaplaceEquationDG with TRIANGLE_LOBATTO / TETRAHEDRON_LOBATTO
  
  {
    EllipticProblem<LaplaceEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE_HDG", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }

  {
    EllipticProblem<LaplaceEquationDG<Dimension2> > var;
    string input_file("src/Program/Unit/Instationary/carre_dg_hyb.ini");
    CheckProblem(var, Complex_wp(0), true, "LAPLACE_HDG", "TRIANGLE_LOBATTO", input_file, false, true, true);
  }

  {
    EllipticProblem<LaplaceEquationDG<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, Complex_wp(0), false, "LAPLACE_HDG", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
  }

  {
    EllipticProblem<LaplaceEquationDG<Dimension3> > var;
    string input_file("src/Program/Unit/Instationary/cube_max.ini");
    CheckProblem(var, Complex_wp(0), true, "LAPLACE_HDG", "TETRAHEDRON_LOBATTO", input_file, false, true, true);
  }
  
  cout << "All tests passed successfully" << endl;
  return FinalizeMontjoie();
}
