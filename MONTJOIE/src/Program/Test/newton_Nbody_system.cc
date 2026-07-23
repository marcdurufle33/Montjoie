#include "Instationary/MontjoieTime.hxx"
#include "Solver/MontjoieSolver.hxx"

using namespace Montjoie;


//! function f(t, y) for the n-body system
class AutonomousFunctionFirst  : public VirtualOdeSystem<Real_wp>,
				 public NonLinearEquations_Newton<Real_wp>
{
public :
  VectReal_wp masse;
  Real_wp G;
  Vector<ofstream> trajectoire;
  bool output;
  Matrix<Real_wp> DF; Vector<int> pivot;
  NewtonSolver<Real_wp> newton_solver;
  VectReal_wp previous_iterate, rhs_newton;
  Real_wp coef_mass, coef_stiff;
  VectReal_wp iterate_dt;
  Real_wp time_step;
  bool overwrite_y1;
  
  AutonomousFunctionFirst()
  {
    overwrite_y1 = false;
    output = false;
    int nb_body = 3;
    G = 1.0;
    masse.Reallocate(nb_body);
    masse(0) = 0.2;
    masse(1) = 1.0;
    masse(2) = 0.1;
    
    if (output)
      {
	trajectoire.Reallocate(nb_body);
	for (int i = 0; i < nb_body; i++)
	  {
	    string name = "P" + to_str(i) + ".txt";
	    trajectoire(i).precision(15);
	    trajectoire(i).open(name.data());        
	  }
      }
  }

  // size of the system 
  int GetM() { return 6*masse.GetM(); }
    
  // f(t, y) for the first-order formulation : y' = f(t, y)
  void EvaluateFunction(const Real_wp& t, const VectReal_wp& pos_vit, VectReal_wp& f,
			bool invert_mass = true, bool source = true)
  {
    int N = masse.GetM();
    Vector<R3> pos(N);
    
    // imposing dpos_i/dt = vi, and retrieving positions
    int offset = 0;
    for (int i = 0; i < N; i++)
      {
        // position of body i
        pos(i).Init(pos_vit(offset), pos_vit(offset+1), pos_vit(offset+2));
        
        // dpos_i/dt = vi
        f(offset) = pos_vit(offset+3);
        f(offset+1) = pos_vit(offset+4);
        f(offset+2) = pos_vit(offset+5);
        
        offset += 6;
      }
    
    // imposing dvi/dt = -G \sum_j m_j (xi - xj) / || xi - xj||^3
    R3 diff_x, acc_i; Real_wp r2, r3;
    offset = 0;
    for (int i = 0; i < N; i++)
      {
        acc_i.Fill(0);
        for (int j = 0; j < N; j++)
          if (j != i)
            {
              diff_x = pos(i) - pos(j);
              r2 = diff_x(0)*diff_x(0) + diff_x(1)*diff_x(1) + diff_x(2)*diff_x(2);
              r3 = r2*sqrt(r2);
              Add(-masse(j)/r3, diff_x, acc_i);
            }
        
        //if (i == N-1)
        //acc_i.Fill(0);
        
        f(offset+3) = G*acc_i(0);
        f(offset+4) = G*acc_i(1);
        f(offset+5) = G*acc_i(2);

        offset += 6;
      }
  }

  // f(t, y, y') for the second-order formulation : y'' = f(t, y, y')
  void EvaluateFunctionS(const Real_wp& t, const VectReal_wp& pos, const VectReal_wp& vit, VectReal_wp& f,
			 bool invert_mass = true, bool source = true)
  {
    int N = masse.GetM();
    Vector<R3> x(N);
    
    int offset = 0;
    for (int i = 0; i < N; i++)
      {
        // position of body i
        x(i).Init(pos(offset), pos(offset+1), pos(offset+2));
        offset += 3;
      }
    
    // f = -G \sum_j m_j (xi - xj) / || xi - xj||^3
    R3 diff_x, acc_i; Real_wp r2, r3;
    offset = 0;
    for (int i = 0; i < N; i++)
      {
        acc_i.Fill(0);
        for (int j = 0; j < N; j++)
          if (j != i)
            {
              diff_x = x(i) - x(j);
              r2 = diff_x(0)*diff_x(0) + diff_x(1)*diff_x(1) + diff_x(2)*diff_x(2);
              r3 = r2*sqrt(r2);
              Add(-masse(j)/r3, diff_x, acc_i);
            }
        
        //if (i == N-1)
        //acc_i.Fill(0);
        
        f(offset) = G*acc_i(0);
        f(offset+1) = G*acc_i(1);
        f(offset+2) = G*acc_i(2);

        offset += 3;
      }
  }

  void EvaluateJacobian(const VectReal_wp& pos)
  {
    Real_wp xi, yi, zi, xj, yj, zj, r, coef, coef2;
    Real_wp dcoef_dx, dcoef_dy, dcoef_dz;
    int N = masse.GetM();
    DF.Reallocate(3*N, 3*N);
    DF.Fill(0);
    for (int i = 0; i < N; i++)
      {
	int offset = 3*i;
	// position de la planete consideree
	xi = pos(offset);
	yi = pos(offset+1);
	zi = pos(offset+2);
       
	// equation dv/dt = -\sum_j G m_j/r^3 (x_i - x_j) 
	for (int j = 0; j < N; j++)
	  if (j != i)
	    {
	      // position de la planete attractrice
	      int off2 = 3*j;
	      xj = pos(off2);
	      yj = pos(off2 + 1);
	      zj = pos(off2 + 2);
	      
	      // distance au carre entre les deux planetes
	      r = sqrt(square(xi-xj) + square(yi-yj) + square(zi-zj));
             
	      // coefficient
	      coef = G*masse(j)/(r*r*r);
	      
	      // contribution de la planete j a l'acceleration de la planete i
	      coef2 = 3*G*masse(j)/pow(r, 5);
	      dcoef_dx = coef2*(xi-xj);
	      dcoef_dy = coef2*(yi-yj);
	      dcoef_dz = coef2*(zi-zj);
             
	      DF(offset, off2) += coef - dcoef_dx*(xi-xj);
	      DF(offset, off2+1) -= dcoef_dy*(xi-xj);
	      DF(offset, off2+2) -= dcoef_dz*(xi-xj);
             
	      DF(offset+1, off2) -= dcoef_dx*(yi-yj);
	      DF(offset+1, off2+1) += coef - dcoef_dy*(yi-yj);
	      DF(offset+1, off2+2) -= dcoef_dz*(yi-yj);
	      
	      DF(offset+2, off2) -= dcoef_dx*(zi-zj);
	      DF(offset+2, off2+1) -= dcoef_dy*(zi-zj);
	      DF(offset+2, off2+2) += coef - dcoef_dz*(zi-zj);

	      DF(offset, offset) += - coef + dcoef_dx*(xi-xj);
	      DF(offset, offset+1) += dcoef_dy*(xi-xj);
	      DF(offset, offset+2) += dcoef_dz*(xi-xj);
             
	      DF(offset+1, offset) += dcoef_dx*(yi-yj);
	      DF(offset+1, offset+1) += - coef + dcoef_dy*(yi-yj);
	      DF(offset+1, offset+2) += dcoef_dz*(yi-yj);
             
	      DF(offset+2, offset) += dcoef_dx*(zi-zj);
	      DF(offset+2, offset+1) += dcoef_dy*(zi-zj);
	      DF(offset+2, offset+2) += - coef + dcoef_dz*(zi-zj);
	    }
      }
    
    /* VectReal_wp pos_tmp(pos), vit, fm(pos), fp(pos);
    Real_wp h = 1e-6;
    for (int j = 0; j < 3*N; j++)
      {
	pos_tmp(j) += h;
	EvaluateFunction(Real_wp(0), pos_tmp, vit, fp);
	pos_tmp(j) -= 2.0*h;
	EvaluateFunction(Real_wp(0), pos_tmp, vit, fm);
	pos_tmp(j) += h;

	for (int i = 0; i < 3*N; i++)
	  {
	    Real_wp df_num = (fp(i) - fm(i))/(2.0*h);
	    DISP(df_num); DISP(DF(i, j));
	    if (abs(df_num - DF(i, j)) > 1e-8)
	      {
		abort();
	      }
	  }
	  } */
    
  }


  // methods called by RunSecondOrderScheme
  void GiveNumberIterations(const Real_wp& dt, int N) { time_step = dt; }
  
  void GiveIterate(int n, const Real_wp& t, VectReal_wp& Y) 
  {
    previous_iterate = Y;
    if (n == 1)
      iterate_dt = Y;
    
    /* if (n == 0)
      cout << "Initial positions = " << Y << endl;

    if (n == 1)
      cout << "position dt = " << Y << endl;

    if (n == 2)
    cout << "position 2 dt = " << Y << endl; */

    if (output)
      {
	int N = masse.GetM();
	int offset = 0;
	for (int i = 0; i < N; i++)
	  {
	    trajectoire(i) << Y(offset) << " " << Y(offset+1) << " " << Y(offset+2) << endl;
	    offset += 6;
	  }
      }
  }
  
  void GiveFinalIterate(int n, const Real_wp& t, VectReal_wp& Y)
  {
    if (output)
      for (int i = 0; i < masse.GetM(); i++)
	trajectoire(i).close();
  }
  
  // methods called by TetaScheme
  void ApplyOperatorKh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& U,
		       const Real_wp& beta, VectReal_wp& V)
  {
    VectReal_wp Fu(U), vit;
    EvaluateFunction(t, U, vit, Fu);
    
    if (beta == Real_wp(0))
      V.Fill(0);
    else
      Mlt(beta, V);
    
    Add(alpha, Fu, V);
  }

  void ApplyOperatorSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& U,
		       const Real_wp& beta, VectReal_wp& V)
  {
    if (beta == Real_wp(0))
      V.Fill(0);
    else
      Mlt(beta, V);
  }

  void ApplyOperatorDh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& U,
		       const Real_wp& beta, VectReal_wp& V)
  {
    if (beta == Real_wp(0))
      V.Fill(0);
    else
      Mlt(beta, V);

    Add(alpha, U, V);
  }

  void ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& U,
				const Real_wp& beta, VectReal_wp& V)
  {
    if (beta == Real_wp(0))
      V.Fill(0);
    else
      Mlt(beta, V);

    Add(alpha, U, V);
  }

  void AddScalarTimeSource(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& V)
  {}
  
  void SolveOperatorDh(VectReal_wp& V) {}
  
  void FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, const Real_wp& beta, const Real_wp& gamma)
  {
    coef_mass = alpha;
    coef_stiff = gamma;
    newton_solver.Init(*this, previous_iterate);
    // true => jacobian is reevaluated at each iteration
    newton_solver.ForceReevaluationJacobian(true);
    newton_solver.SetMaxNumberOfIterations(20);
  }

  void SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& rhs, VectReal_wp& U)
  {
    rhs_newton = rhs;
    U = previous_iterate;
    newton_solver.Solve(*this, U);
  }
  
  
  void SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha=1.0)
  {
    if (previous_iterate.GetM() == 0)
      previous_iterate = Y;
    
    if ((t == time_step) && overwrite_y1)
      Y = iterate_dt;
  }

  // methods called by NewtonSolver
  void ComputeScheme(const VectReal_wp& x, VectReal_wp& f)
  {
    // evaluation of f = -Y + (coef_mass Dh + coef_stiff Kh) X
    for (int i = 0; i < x.GetM(); i++)
      f(i) = -rhs_newton(i);
    
    ApplyOperatorDh(coef_mass, Real_wp(0), x, Real_wp(1), f);
    ApplyOperatorKh(-coef_stiff, Real_wp(0), x, Real_wp(1), f);    
  }

  void ComputeAndFactoriseDiff(const VectReal_wp& x, const VectReal_wp& scale)
  {
    EvaluateJacobian(x);
    
    Mlt(-coef_stiff, DF);
    for (int i = 0; i < x.GetM(); i++)
      DF(i, i) += coef_mass;
    
    GetLU(DF, pivot);
  }

  void SolveDifferential(const VectReal_wp& x, VectReal_wp& y)
  {
    y = x;
    SolveLU(DF, pivot, y);
  }

  Real_wp GetNorm2Vector(const VectReal_wp& x)
  {
    return Norm2(x);
  }

};


Real_wp GetRelativeError(const VectReal_wp& x, const VectReal_wp& y)
{
  Real_wp err = 0, sum = 0;
  for (int i = 0; i < x.GetM(); i++)
    {
      err += square(x(i) - y(i));
      sum += y(i)*y(i);
    }
  
  err = sqrt(err); sum = sqrt(sum);
  return err/sum;
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  Real_wp Tmax(1), dt, one(1), zero(0);
  dt = one/10;

  RungeKuttaNystrom_Iterator<Real_wp> scheme;
  scheme.SetOrder(8);
  
  AutonomousFunctionFirst fct_auto;
  int N = fct_auto.GetM();
  
  VectReal_wp y0(N/2), dy0(N/2), yinit(N/2), dyinit(N/2);

  y0(0) = 1.0; y0(1) = 0.0; y0(2) = 0.0;
  dy0(0) = 0.0; dy0(1) = 1.1; dy0(2) = 0.04;
  y0(3) = 0.0; y0(4) = 0.0; y0(5) = 0.0;
  dy0(3) = 0.0; dy0(4) = 0.0; dy0(5) = 0.0;
  y0(6) = 0.5; y0(7) = 0.8; y0(8) = 0.2;
  dy0(6) = -0.04; dy0(7) = 0.95; dy0(8) = 0.2;
  
  dyinit = dy0;
  yinit = y0;
  RunSecondOrderScheme(zero, Tmax, dt, y0, dy0, fct_auto, scheme);

  VectReal_wp yinit_dt(fct_auto.iterate_dt);
  VectReal_wp y_0(y0), dy_0(dy0);
  
  y0 = yinit;
  dy0 = dyinit;
  RunSecondOrderScheme(zero, Tmax, dt/2, y0, dy0, fct_auto, scheme);

  VectReal_wp yinit_dt2(fct_auto.iterate_dt);
  VectReal_wp y1(y0), dy1(dy0);

  y0 = yinit;
  dy0 = dyinit;
  RunSecondOrderScheme(zero, Tmax, dt/4, y0, dy0, fct_auto, scheme);

  VectReal_wp yinit_dt4(fct_auto.iterate_dt);
  VectReal_wp y2(y0), dy2(dy0);
  
  Real_wp err = GetRelativeError(y_0, y1);
  Real_wp err_half = GetRelativeError(y1, y2);
  
  int r = 8; Real_wp err_order = 0.2;
  Real_wp order = (log(err) - log(err_half))/log(2.0);
  DISP(err); DISP(err_half);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  /* if (order < r*(1.0- err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
      } */


  TetaScheme_Iterator<Real_wp> teta_scheme;
  teta_scheme.SetTheta(Real_wp(1)/12);
  fct_auto.overwrite_y1 = true;

  fct_auto.iterate_dt = yinit_dt;
  y0 = yinit;
  dy0 = dyinit;
  RunSecondOrderScheme(zero, Tmax, dt, y0, dy0, fct_auto, teta_scheme);
  
  VectReal_wp yteta_0(y0), dy_teta0(dy0);
  
  fct_auto.iterate_dt = yinit_dt2;
  y0 = yinit;
  dy0 = dyinit;
  RunSecondOrderScheme(zero, Tmax, dt/2, y0, dy0, fct_auto, teta_scheme);
  VectReal_wp yteta_1(y0), dy_teta1(dy0);

  fct_auto.iterate_dt = yinit_dt4;
  y0 = yinit;
  dy0 = dyinit;
  RunSecondOrderScheme(zero, Tmax, dt/4, y0, dy0, fct_auto, teta_scheme);
  VectReal_wp yteta_2(y0), dy_teta2(dy0);
  
  err = GetRelativeError(yteta_0, yteta_1);
  err_half = GetRelativeError(yteta_1, yteta_2);
  
  DISP(GetRelativeError(yteta_0, y_0));
  DISP(GetRelativeError(yteta_1, y1));
  DISP(GetRelativeError(yteta_2, y2));

  r = 2; err_order = 0.2;
  order = (log(err) - log(err_half))/log(2.0);
  DISP(err); DISP(err_half);
  cout << "Order expected : " << r << ", Order measured : " << order << endl;
  if (order < r*(1.0- err_order))
    {
      DISP(err); DISP(err_half); DISP(r); DISP(order);
      cout << "Scheme incorrect " << endl;
      abort();
    }

  return FinalizeMontjoie();
}
