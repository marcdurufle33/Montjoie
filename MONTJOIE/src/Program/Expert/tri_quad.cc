#define MONTJOIE_WITH_TWO_DIM

// #include "lm.h"

//#ifndef LM_DBL_PREC
//#error Demo program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
//#endif


#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

class TriQuad
{
public :
  static int order, order_quad, nb_pts;
  static int nb_eval;
  static IVect OrderEta;
  static Matrix<Real_wp> ValExact;
  
  static int type_InitGuess;
  enum{Classic, Random};
  enum{GSL, LEVMAR, SA};
  
  TriQuad() { }
  
  static int GetNbUnknowns()
  {
    int n = order+1;
    for (int i = 0; i <= order; i++)
      n += 2*OrderEta(i);
    
    return n;
  }
  
  static int GetNbEquations()
  {
    int N = GetNbUnknowns();
    return max((order_quad+1)*(order_quad+2)/2,N);
  }
  
  template <class Vector1> static void FindInitGuess(Vector1& xsol)
  {
    Vector<double> x, y, w, wy;
    int nb = order+1;
    
    if (type_InitGuess == Classic)
      {
	ComputeGaussJacobi(y, wy, order, 1.0, 0.0);
	for (int i = 0; i <= order; i++)
	  {
	    xsol(i) = y(i);
	    ComputeGaussLegendre(x, w, OrderEta(i)-1);
	    for (int j = 0; j < OrderEta(i) ; j++)
	      xsol(nb++) = (1-y(i))*x(j);
	  }
	
	for (int i = 0; i <= order; i++)
	  {
	    ComputeGaussLegendre(x, w, OrderEta(i)-1);
	    for (int j = 0; j < OrderEta(i) ; j++)
	      xsol(nb++) = w(j)*wy(i);
	  }
      } 
    else
      {
	for (int i = 0; i <= order; i++)
	  {
	    xsol(i) = double(rand())/RAND_MAX;
	    for (int j = 0; j < OrderEta(i) ; j++)
	      xsol(nb++) = (1-xsol(i))*double(rand())/RAND_MAX;
	  }
	
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j < OrderEta(i) ; j++)
	    xsol(nb++) = 0.5/nb_pts;
      }
    
  }
  
  template <class Vector1> static void EvaluateObjective(const Vector1& pts_pds, Real_wp& feval)
  {
    feval = 0;
    
    for (int m = 0; m <= order_quad; m++)
      for (int n = 0; n <= order_quad-m ; n++)
	{ 
	  int nb = order+1;
	  Real_wp sum = 0.;
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j < OrderEta(i) ; j++)
	      {
		sum += pts_pds(nb_pts+nb)*pow(pts_pds(nb),m)*pow(pts_pds(i),n); 
		nb ++;
	      }
	  sum -= ValExact(m,n);
	  
	  feval += sum*sum;
	}
    
    nb_eval++;
  }
  
  template <class Vector1> static void EvaluateFunction(const Vector1& pts_pds, Vector1& grad)
  {
    
    grad.Fill(0);
    for (int m = 0; m <= order_quad; m++)
      for (int n = 0; n <= order_quad-m ; n++)
	{ 
	  int nb = order+1;
	  Real_wp sum = 0.;
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j < OrderEta(i) ; j++)
	      {
		sum += pts_pds(nb_pts+nb)*pow(pts_pds(nb),m)*pow(pts_pds(i),n); 
		nb ++;
	      }
	  sum -= ValExact(m,n);
	  nb = order+1;
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j < OrderEta(i) ; j++)
	      {
		if (n > 0)
		  grad(i) += 2*sum*n* pts_pds(nb_pts+nb)*pow(pts_pds(nb),m)*pow(pts_pds(i),n-1);
		if (m > 0)
		  grad(nb) +=  2*sum*m*pts_pds(nb_pts+nb)*pow(pts_pds(nb),m-1)*pow(pts_pds(i),n);
		grad(nb_pts+nb) +=  2*sum*pow(pts_pds(nb),m)*pow(pts_pds(i),n);
		nb ++;
	      }
	}
    
    nb_eval++;
  }
  
  
  static Real_wp GetStoppingCriterion()
  {
    return 1E-14;
  }
  
#ifdef MONTJOIE_WITH_GSL
  static double my_f(const gsl_vector* v, void* params)
  {
    int N = GetNbUnknowns();
    Vector<Real_wp> pts_pds(N); 
    Real_wp feval;
    for (int i = 0; i < N; i++)
      pts_pds(i) = gsl_vector_get(v, i);
    
    EvaluateObjective(pts_pds, feval);
    
    return toDouble(feval);
  }
  
  static void my_df(const gsl_vector* v, void* params, gsl_vector* df)
  {
    int N = GetNbUnknowns();
    Vector<Real_wp> pts_pds(N), grad(N);
    
    for (int i = 0; i < N; i++)
      pts_pds(i) = gsl_vector_get(v, i);
    
    EvaluateFunction(pts_pds, grad);
    
    for (int i = 0; i < N; i++)
      gsl_vector_set(df, i, grad(i));    
  }
  
  static void my_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df)
  {
    my_df(v, params, df);
    *f = my_f(v, params);
  }
#endif
  
  static void evalFunction(double* x, double* f, int nb_unknowns, int nb_equations, void* data)
  {
    int num_unknown = 0;
    for (int l = 0; l < GetNbEquations(); l++)
      f[l] = 0;
    
    for (int m = 0; m <= order_quad; m++)
      for (int n = 0; n <= order_quad-m ; n++)
	{ 
	  int nb = order+1;
	  Real_wp sum = 0.;
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j < OrderEta(i) ; j++)
	      {
		sum += x[nb_pts+nb]*pow(x[nb],m)*pow(x[i],n); 
		nb ++;
	      }
	  
	  sum -= ValExact(m,n);
	  f[num_unknown] = sum;
	  num_unknown++;
	}
    
  }
  
  static void evalGradient(double* x, double* df, int nb_unknowns, int nb_equations, void* data)
  {
    int num_component = 0;
    
    for (int l = 0; l < GetNbEquations()*GetNbUnknowns(); l++)
      df[l] = 0;
    
    for (int m = 0; m <= order_quad; m++)
      for (int n = 0; n <= order_quad-m ; n++)
	{ 
	  int nb = order+1;
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j < OrderEta(i) ; j++)
	      {
		if (n > 0)
		  df[num_component+i] = n*x[nb_pts+nb]*pow(x[nb],m)*pow(x[i],n-1);
		else
		  df[num_component+i] = 0;
		
		if (m > 0)
		  df[num_component+nb] =  m*x[nb_pts+nb]*pow(x[nb],m-1)*pow(x[i],n);
		else
		  df[num_component+nb] = 0;
		
		df[num_component+nb_pts+nb] =  pow(x[nb],m)*pow(x[i],n);
		// DISP(num_component); DISP(df[num_component+nb]);
		nb ++;
	      }
	  
	  num_component += nb_unknowns;
	}
  }
  
  static bool AcceptableRule(Vector<double>& xsol)
  {
    bool rule_ok = true;
    
    int nb = order+1;
    for (int i = 0; i < order+1; i++)
      {
	if ((xsol(i) >1) || (xsol(i)<0))
	  rule_ok = false;
	
	for (int j = 0; j < TriQuad::OrderEta(i) ; j++)
	  {
	    if ((xsol(nb) > 1-xsol(i)) || (xsol(nb)<0))
	      rule_ok = false;
	    
	    nb++;
	  }
      }
    
    for (int i = order+1+nb_pts; i < xsol.GetM(); i++)
      if (xsol(i) < 0)
	rule_ok = false;
    
    return rule_ok;
  }
  
};

int TriQuad::order, TriQuad::order_quad, TriQuad::nb_pts;
int TriQuad::nb_eval;
IVect TriQuad::OrderEta;
Matrix<Real_wp> TriQuad::ValExact;
int TriQuad::type_InitGuess;


template<class FunctionParam, class Vector1>
void FindQuadratureTriangle(FunctionParam& fct, Vector1& xsol, Real_wp& feval, int type_algo)
{
  TriQuad::nb_eval = 0;
  int N = fct.GetNbUnknowns();
  
  if (type_algo == fct.GSL)
    {
      // GSL class for definition of function f and derivative df
      gsl_multimin_function_fdf my_func;
      
      my_func.f = &fct.my_f;
      my_func.df = &fct.my_df;
      my_func.fdf = &fct.my_fdf;
      my_func.n = N;
      my_func.params = NULL;
      
      const gsl_multimin_fdfminimizer_type *T;
      gsl_multimin_fdfminimizer *s;
      
      // solution stored in a gsl vector
      gsl_vector *x;
      x = gsl_vector_alloc(my_func.n);
      
      // initial guess
      for (int i = 0; i < xsol.GetM(); i++)
	gsl_vector_set(x, i, xsol(i));
      
      T = gsl_multimin_fdfminimizer_vector_bfgs;
      // T = gsl_multimin_fdfminimizer_conjugate_pr;
      // T = gsl_multimin_fdfminimizer_conjugate_fr;
      // T = gsl_multimin_fdfminimizer_steepest_descent;
      s = gsl_multimin_fdfminimizer_alloc (T, my_func.n);
      
      gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.01, 1e-4);
      
      size_t iter = 0;
      int status;
      
      do
	{
	  iter++;
	  status = gsl_multimin_fdfminimizer_iterate(s);
	  
	  if (status)
	    break;
	  
	  status = gsl_multimin_test_gradient (s->gradient, fct.GetStoppingCriterion());
	  
	}
      while (status == GSL_CONTINUE && iter < 10000);
      
      for (int i = 0; i < xsol.GetM(); i++)
	xsol(i) = gsl_vector_get(s->x, i);
      
      // we return CFL
      fct.EvaluateObjective(xsol, feval);
      feval = sqrt(feval);
    }
  else
    if (type_algo == fct.LEVMAR)
      
      {
	int M = fct.GetNbEquations();
	Vector<Real_wp> x(N), f(M);
	x.Fill(0); f.Fill(0);
	x = xsol;
	/* double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing
	//opts[4]=-LM_DIFF_DELTA; // specifies central differencing to approximate Jacobian; more accurate but more expensive to compute!
	
	dlevmar_der(fct.evalFunction, fct.evalGradient, x.GetData(), f.GetData(), N, M, 1000, opts, info, NULL, NULL, NULL); // with analytic Jacobian
	xsol = x; 
	fct.evalFunction(x.GetData(), f.GetData(), N, M, NULL);
	fct.EvaluateObjective(xsol, feval); feval = sqrt(feval); */
	abort();
	// DISP(x); DISP(f);
      }
    else
      {
	int M = fct.GetNbEquations();
	int nbx = N-TriQuad::nb_pts;
	Vector<Real_wp> x(N), W(M), tau(M);
	IVect pivot(TriQuad::nb_pts);
	x.Fill(0); W.Fill(0);
	x = xsol;
	// evaluation initiale
	fct.EvaluateObjective(x, feval); feval = sqrt(feval);
	// température initiale
	Real_wp T=0.1, feval_t, p, Delta;

	int nbit = 1;
	// Matrix<Real_wp> A(M,TriQuad::nb_pts);
	// Matrix<Real_wp> A(TriQuad::nb_pts, TriQuad::nb_pts);
	Matrix<Real_wp> A(M, TriQuad::nb_pts);
	int nbr;
	
	Real_wp val_expo = 1.0;
	while ((T> 1E-6) && (nbit < 10000) && (val_expo > 1e-20))
	  {
	    int nb = TriQuad::order+1;
	    for (int i = 0; i < TriQuad::order+1; i++)
	      {
		bool pts_outside = true;
		while(pts_outside)
		  {
		    xsol(i) = x(i) + 0.1*(2*double(rand())/RAND_MAX-1);
		    pts_outside = ((xsol(i) >1) || (xsol(i)<0));
		  }
		
		for (int j = 0; j < TriQuad::OrderEta(i) ; j++)
		  {
		    bool pts_outside = true;
		    while(pts_outside)
		      {
			xsol(nb) = x(nb) + 0.1*(2*double(rand())/RAND_MAX-1);
			pts_outside = ((xsol(nb) > 1-xsol(i)) || (xsol(nb)<0));
		      }
		    		    
		    bool poids_positive = true;
		    while(poids_positive)
		      {
			xsol(nb+TriQuad::nb_pts) = x(nb+TriQuad::nb_pts)*(1.0 - 0.05*(2*double(rand())/RAND_MAX-1));
			poids_positive = (xsol(nb+TriQuad::nb_pts) < 0);
		      }
		    nb++;
		  }
		
	      }

	    /*
	    A.Fill(0);
	    nbr = 0;
	    
	    // W.Reallocate(TriQuad::nb_pts);
	    W.Reallocate(M);

	    for (int idiag = 0; idiag <= TriQuad::order_quad; idiag++)
	      for (int jdiag = 0; jdiag <= idiag; jdiag++)
		{
		  int n = jdiag;
		  int m = idiag-jdiag;
		  int nbc = 0;
		  
		  if (nbr < W.GetM())
		    {
		      DISP(m); DISP(n);
		      W(nbr) = TriQuad::ValExact(m,n);
		      for (int i = 0; i <= TriQuad::order; i++)
			for (int j = 0; j < TriQuad::OrderEta(i) ; j++)
			  {
			    A(nbr,nbc) = pow(xsol(nbc+TriQuad::order+1),m)*pow(xsol(i),n); 
			    nbc++;
			  }
		      DISP(nbr);
		    } 
		  nbr++;
		}
	    
	    DISP(M);
	    W.WriteText("Bh.dat");
	    A.WriteText("Mh.dat");
	    GetQR(A,tau);
	    for (int i = (TriQuad::order_quad+1)*(TriQuad::order_quad+2)/2; i < A.GetM(); i++)
	      A(i, i) = 1.0;
	    
	    DISP(A); SolveQR(A,tau,W); DISP(W);
	    // GetLU(A, pivot); SolveLU(A, pivot, W); 
	    
	    for (int i = nbx; i < N; i++)
	      xsol(i) = W(i-nbx);
	    */
 
	    fct.EvaluateObjective(xsol, feval_t);
	    Delta = feval_t-feval;

	    if (Delta<0)
	      {
		x = xsol;
		feval = feval_t;
	      }
	    else
	      {
		p = double(rand())/RAND_MAX;
		val_expo = exp(-1000*Delta/T);
		if (p < val_expo)
		  {
		    x = xsol;
		    feval = feval_t;
		  }
	      }
	    T = T*0.998;
	    nbit++;
	  }
	xsol = x;
	feval = sqrt(feval);
      }
  
}

int main(int argc, char** argv)
{
  if (argc != 3)
    {
      cout<<"entrez l'ordre d'approximation"<<endl;
      return -1;
    }
  
  InitMontjoie(argc, argv);
  
  cerr.precision(15);
  srand(time(NULL));

  TriQuad::order = to_num<int>(string(argv[1]));
  TriQuad::nb_pts = 0;
  TriQuad::OrderEta.Reallocate(TriQuad::order+1);
  for (int i = 0; i <= TriQuad::order; i++)
    {
      TriQuad::OrderEta(i) = TriQuad::order+1-i;
    }
  
 if (TriQuad::order == 3)
    {
      TriQuad::OrderEta(1)++;
      TriQuad::OrderEta(2)++;
      // TriQuad::OrderEta(3)++;
    }

  if (TriQuad::order == 4)
    {
      TriQuad::OrderEta(1)++;
      TriQuad::OrderEta(2)++;
      TriQuad::OrderEta(3)++;
      TriQuad::OrderEta(4)++;
    }
  
  for (int i = 0; i <= TriQuad::order; i++)
    {
      TriQuad::nb_pts += TriQuad::OrderEta(i);
    }

  if (TriQuad::order == 1)
    TriQuad::order_quad = 2;
  else if (TriQuad::order < 3)
    TriQuad::order_quad = 2*TriQuad::order-1;
  else
    TriQuad::order_quad = 2*TriQuad::order-2;

  TriQuad::type_InitGuess = TriQuad::Random;
  // TriQuad::type_InitGuess = TriQuad::Classic;

  TriQuad::ValExact.Reallocate(TriQuad::order_quad+1,TriQuad::order_quad+1);
  TriQuad::ValExact.Fill(0);
  for (int m = 0; m <= TriQuad::order_quad; m++)
    for (int n = 0; n <= TriQuad::order_quad-m; n++)
      TriQuad::ValExact(m,n) = tgamma(m+1)*tgamma(n+1)/tgamma(m+n+3);

  DISP(TriQuad::ValExact);

  int nb_iter_max =  to_num<int>(string(argv[2]));

  int N = TriQuad::GetNbUnknowns(); DISP(N); DISP(TriQuad::GetNbEquations());
  Vector<Real_wp> xmin(N);
  Real_wp fmin = 1;
  for (int i = 0; i< nb_iter_max; i++)
    {
      TriQuad fct;
      Vector<Real_wp> xsol(N);
      Real_wp feval;
      // solution initiale
      fct.FindInitGuess(xsol);
      fct.EvaluateObjective(xsol, feval);
      // computing quadrature rule by minimizing an objective function
      // FindQuadratureTriangle(fct, xsol, feval, TriQuad::SA);
      FindQuadratureTriangle(fct, xsol, feval, TriQuad::LEVMAR);
      //DISP(xsol); DISP(feval);
      if (TriQuad::AcceptableRule(xsol))
	if (feval < fmin)
	  {
	    fmin = feval;
	    xmin = xsol;
	  }
    }

  DISP(xmin); DISP(fmin);
  
  // retrieving points and weights
  Vector<R2> Points(TriQuad::nb_pts); Vector<Real_wp> Weights(TriQuad::nb_pts);
  
  int nb = TriQuad::order+1, nb_quad = 0;
  for (int i = 0; i <= TriQuad::order; i++)
    for (int j = 0; j < TriQuad::OrderEta(i); j++)
      {
	Points(nb_quad) = R2(xmin(nb), xmin(i)); 
	Weights(nb_quad) = xmin(nb+TriQuad::nb_pts);
	nb++; nb_quad++;
      }
  
  DISP(Points); DISP(Weights);
  
  // checking accuracy of quadrature rule
  // on teste la formule de quadrature
  for (int k = 0; k <= 2*TriQuad::order; k++)
    {
      bool integration_ok = true;
      for (int i = 0; i <= k; i++)
	{
	  int j = k-i;
	  // on integre x^i y^j
	  Real_wp val_exact = tgamma(i+1)*tgamma(j+1)/tgamma(k+3);
	  Real_wp val_app = 0;
	  for (int m = 0; m < TriQuad::nb_pts; m++)
	    val_app += Weights(m)*pow(Points(m)(0), i)*pow(Points(m)(1), j);

	  if (abs(val_app-val_exact) > 1e-10)
	    {
	      DISP(val_app); DISP(val_exact); DISP(i); DISP(j);
	      integration_ok = false;
	    }
	}
      
      if (!integration_ok)
	{
	  cout<<"L'ordre "<<k<<" n'est pas integre exactement"<<endl;
	  abort();
	}
    }
  
  cout<<"Quadrature rule perfect :) "<<endl;
  
  return FinalizeMontjoie();
}
