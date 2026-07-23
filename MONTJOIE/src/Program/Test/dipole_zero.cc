#include "Montjoie.hxx"
#include "Output/MontjoieOutput.hxx"
#include "Output/OutputOpenCV.cxx"

using namespace Montjoie;

class DipoleEquation
  : public NonLinearEquations_Base<Real_wp>,
    public NonLinearEquations_Newton<Real_wp>
{
public:
  Real_wp gamma, Omega0, c0, Lambda0, z;
  Complex_wp Omega;
  Matrix<Real_wp> DF_store;
  bool mode_plus;

public :
  DipoleEquation()
  {
    c0 = PhysicalConstant::speed_light; DISP(c0);
    //gamma = 1.125e13;
    gamma = 1.125e8;
    Omega0 = 3.4586e15;
    Omega = Omega0 - Iwp*gamma;
    Lambda0 = 545e-9;
    z = 1.0/100*Lambda0;
  }

  DipoleEquation(Real_wp z0,bool mode)
  {
    c0 = 1.0;
    gamma = 0.1;
    Omega0 = 2.0*pi_wp;
    Omega = Omega0 - Iwp*gamma;
    Lambda0 = 1.0;
    z = z0;
    mode_plus = mode;
  }

  void GetGreenTensor(const R3& x, const R3& y, const Complex_wp& omega,
		      TinyMatrix<Complex_wp, Symmetric, 3, 3>& G)
  {
    Complex_wp omega2 = omega*omega;

    // in 3-D \phi(x,y) = exp(ik|x-y|) / (4 pi |x-y| )
    R3 xMinusy = x - y;

    // T = |x-y|   T2 = |x-y|^2   arg = ik |x-y|
    Real_wp T = Norm2(xMinusy); Complex_wp arg = Iwp*omega*T;
    Real_wp inv_T = Real_wp(1)/T, inv_T2 = inv_T*inv_T;
    
    // evaluation of green function
    Complex_wp phi = exp(arg) * inv_T / (4*pi_wp);
    Complex_wp phi_div_T = phi * inv_T;
        
    // evaluation of hessian matrix
    // initialization to (x_m - y_m) (x_l - y_l)
    GetNormalProjector(xMinusy, G);    

    Complex_wp alpha = (-omega2 - Real_wp(3)*Iwp*omega*inv_T + Real_wp(3) * inv_T2) * inv_T;
    Mlt(alpha, G);
    alpha = Iwp*omega - inv_T;
    G(0, 0) += alpha; G(1, 1) += alpha; G(2, 2) += alpha;
    
    Mlt(phi_div_T / omega2, G);

    // Green's tensor
    // G = \phi(x,y) I + 1/k^2 \nabla_y \nabla_y \phi(x,y)
    G(0, 0) += phi; G(1, 1) += phi; G(2, 2) += phi;
  }
  
  void ComputeGreenKernel(const R3& x, const R3& y, const Complex_wp& omega,
			  TinyMatrix<Complex_wp, Symmetric, 3, 3>& G)
  {
    TinyMatrix<Complex_wp, Symmetric, 3, 3> Gphi1, Gphi2, Gpq, Gs1, Gs2, Gxy;
    R3 s1(Real_wp(0), Real_wp(0), -z/2);
    R3 s2(Real_wp(0), Real_wp(0), z/2);
    
    // on evalue le champ "incident" sur les deux dipoles
    GetGreenTensor(y, s1, omega, Gphi1);
    GetGreenTensor(y, s2, omega, Gphi2);

    // coef sigma
    Complex_wp k = omega / c0;
    Complex_wp sigma = Real_wp(6)*pi_wp*gamma / (k*(Omega - omega));

    // on construit la matrice a inverser
    TinyMatrix<Complex_wp, Symmetric, 2, 2> mat;

    GetGreenTensor(s1, s2, omega, Gpq);
    GetGreenTensor(y, x, omega, Gxy);
	
    mat(0, 0) = Complex_wp(1, 0);
    mat(1, 1) = Complex_wp(1, 0);
    mat(0, 1) = -Gpq(0, 0)*sigma;
    DISP(omega); DISP(Gpq); DISP(sigma); DISP(k); DISP(mat);
    GetInverse(mat);
    DISP(mat);
    
    // on evalue G(x, s_i)
    GetGreenTensor(x, s1, omega, Gs1);
    GetGreenTensor(x, s2, omega, Gs2);
    
    // boucle sur les colonnes de G
    TinyVector<Complex_wp, 3> colG1, colG2, prod, vec_u1, vec_u2, Etot;
    TinyVector<Complex_wp, 2> rhs, sol;
    for (int j = 0; j < 3; j++)
      {
	GetCol(Gphi1, j, colG1);
	GetCol(Gphi2, j, colG2);

	rhs(0) = colG1(0); rhs(1) = colG2(0);
	Mlt(mat, rhs, sol);

	// la solution finale
	// s'ecrit Einc + \sum G(x, s_p) sigma_p w_p
	prod.Zero(); prod(0) = sigma * sol(0);
	Mlt(Gs1, prod, vec_u1);
	prod.Zero(); prod(0) = sigma * sol(1);
	Mlt(Gs2, prod, vec_u2);

	GetCol(Gxy, j, Etot);
	Etot = Etot + vec_u1 + vec_u2;
	SetCol(Etot, j, G);
      }
  }

  void EvaluateFunction(const VectReal_wp& x, VectReal_wp& y)
  {
    Complex_wp d_omega(x(0), x(1));
    Complex_wp k = (Omega + d_omega) / c0;
    Complex_wp kz = k*z;
    Complex_wp coef = -Real_wp(1.5)*gamma/ kz *
      (Real_wp(1) + Iwp/(k*z) - Real_wp(1)/(kz*kz))*exp(Iwp*k*z);

    DISP(d_omega);
    DISP(coef/d_omega);
    
    if (mode_plus)
      {
	y(0) = realpart(coef-d_omega);
	y(1) = imagpart(coef-d_omega);
      }
    else
      {
	y(0) = realpart(coef+d_omega);
	y(1) = imagpart(coef+d_omega);
      }
  }

  void ComputeScheme(const VectReal_wp& x, VectReal_wp& F)
  {
    EvaluateFunction(x, F);
  }

  // computes Jac = DF(x)
  void EvaluateJacobian(const VectReal_wp& x0, Matrix<Real_wp>& DF)
  {
    Real_wp h = pow(epsilon_machine, Real_wp(1)/3);
    VectReal_wp x(x0), fm(2), fp(2);
    for (int j = 0; j < 2; j++)
      {
	x = x0;
	x(j) -= h;
	EvaluateFunction(x, fm);

	x(j) += 2*h;
	EvaluateFunction(x, fp);

	//DISP(fm); DISP(fp);
	
	DF(0, j) = (fp(0) - fm(0)) / (2*h);
	DF(1, j) = (fp(1) - fm(1)) / (2*h);
      }

    //DISP(h); DISP(x0); DISP(DF);
  }
  
  void ComputeAndFactoriseDiff(const VectReal_wp& x, const VectReal_wp& scale)
  {
    //DISP(x); DISP(Omega);
    DF_store.Reallocate(2, 2);
    EvaluateJacobian(x, DF_store);

    //DISP(DF_store);
    
    GetInverse(DF_store);
  }

  void SolveDifferential(const VectReal_wp& b, VectReal_wp& x)
  {
    Mlt(DF_store, b, x);
    //DISP(b); DISP(x);
  }


  Real_wp GetNorm2Vector(const VectReal_wp& x)
  {
    return Norm2(x);
  }

};

void FollowRoots()
{
  Real_wp z0 = 1.0;
  Real_wp dz = 1.0/400;
  Real_wp zmin = 1.0;
  Real_wp zmax = 60.0;
  Real_wp z = z0;
  Real_wp z1 = 1;
  Real_wp log_dz = 0.01;
  DipoleEquation eq(z,true);
  DipoleEquation eq_m(z,false);

  VectReal_wp x_init(2),fvec(2), x_init_0, x_init_m_0;
  x_init(0)=0;x_init(1)=0;
  VectReal_wp x_init_m(2),fvec_m(2);
  x_init_m(0)=0;x_init_m(1)=0;


  Real_wp norme_ref(1); int nb_iter(1);
    
  NewtonSolver<Real_wp> solver;
  solver.ForceReevaluationJacobian();
  solver.SetPrintLevel(0);
  solver.SetMaxNumberOfIterations(100);
  solver.EnableDecreasingResidue(false);
  solver.SetStoppingCriterion(1e-11, 1e-10);

  if (false)
    {
      ofstream file_out("racine_double.dat");
      ofstream file_out_m("racine_m_double.dat");
      file_out.precision(16);
      file_out_m.precision(16);
      
      bool decreasing_step = true;
      
      Vector<TinyVector<Real_wp, 5> > tab_decrease;
      int n = 0; int nb_iteration = 0;
      while (z < zmax)
	{
	  nb_iteration++;
	  if (nb_iteration%1000 == 0)
	    cout << "z = " << z << endl;
	  
	  eq.z = z; eq_m.z = z;
	  solver.Init(eq, x_init); nb_iter = 0;
	  solver.Newton(eq, x_init, nb_iter, norme_ref, true, false);
	  
	  solver.Init(eq_m, x_init_m); nb_iter = 0;
	  solver.Newton(eq_m, x_init_m, nb_iter, norme_ref, true, false);
	  
	  eq.EvaluateFunction(x_init, fvec);
	  //DISP(x_init); DISP(fvec); DISP(z);
	  if (Norm2(fvec) > 1e-10)
	    {
	      cout << "Not converged for z = " << z << endl;
	      abort();
	    }
	  
	  eq_m.EvaluateFunction(x_init_m, fvec_m);
	  //DISP(x_init_m); DISP(fvec_m);
	  
	  if (Norm2(fvec_m) > 1e-10)
	    {
	      cout << "Minus not converged for z = " << z << endl;
	      abort();
	    }
	  
	  if (n == 0)
	    {
	      x_init_0 = x_init;
	      x_init_m_0 = x_init_m;
	    }
	  
	  if (decreasing_step)
	    {
	      TinyVector<Real_wp, 5> tmp(z, x_init(0), x_init(1), x_init_m(0), x_init_m(1));
	      tab_decrease.PushBack(tmp);
	      
	      z /= pow(10.0,log_dz);
	      if (z < zmin)
		{
		  decreasing_step = false;
		  for (int i = tab_decrease.GetM()-1; i >= 0; i--)
		    {
		      tmp = tab_decrease(i);
		      file_out << tmp(0) << " " << tmp(1) << " "  << tmp(2) << endl;
		      file_out_m << tmp(0) << " " << tmp(3) << " "  << tmp(4) << endl;
		    }
		  
		  z = z0*pow(10.0, log_dz);
		  x_init = x_init_0;
		  x_init_m = x_init_m_0;
		}
	    }
	  else
	    {
	      int Ndisplay = 1;
	      if (z > 1e2)
		Ndisplay = 3;
	      else if (z > 1e3)
		Ndisplay = 11;
	      else if (z > 1e4)
		Ndisplay = 97;
	      else if (z > 3e4)
		Ndisplay = 367;
	      else if (z > 1e5)
		Ndisplay = 951;
	      else if (z > 3e5)
		Ndisplay = 3221;
	      
	      if (nb_iteration%Ndisplay == 0)
		{
		  file_out << z << " " <<  x_init(0) << " " <<  x_init(1) << '\n';
		  file_out_m << z << " " <<  x_init_m(0) << " " <<  x_init_m(1) << '\n';
		}
	      	      
	      if (z < z1)
		z*=pow(10.0,log_dz);
	      else
		z+=dz;
	    }
	  
	  n++;
	}
      
      file_out.close();
      file_out_m.close();
    }

  //return;
  
  z0 = 28.3; z = z0;
    
  VectReal_wp x, y;
  Linspace(-pi_wp/2, pi_wp/2, 1000, x);
  Linspace(Real_wp(-0.4), Real_wp(0.05), 200, y);
  
  eq.z = z0;
  int m = x.GetM(), n = y.GetM();
  Matrix<Real_wp> Val(m, n);
  for (int i = 0; i < x.GetM(); i++)
    for (int j = 0; j < y.GetM(); j++)
      {
	VectReal_wp pt(2), f(2);
	pt(0) = x(i); pt(1) = y(j);
	eq.EvaluateFunction(pt, f);
	Val(i, j) = 1.0-log10(Norm2(f));
      }

  VectComplex_wp roots;
  for (int i = 1; i < x.GetM()-1; i++)
    for (int j = 1; j < y.GetM()-1; j++)
      {
	Real_wp max_neigh = -1e300;
	for (int p = -1; p <= 1; p++)
	  for (int q = -1; q <= 1; q++)
	    if ((p != 0) || (q != 0))
	      max_neigh = max(Val(i+p, j+q), max_neigh);
	
	if (Val(i, j) > max_neigh)
	  {
	    //DISP(i); DISP(j); DISP(x(i)); DISP(y(j));

	    x_init(0) = x(i);
	    x_init(1) = y(j);
      
	    solver.Init(eq, x_init); nb_iter = 0;
	    solver.Newton(eq, x_init, nb_iter, norme_ref, true, false);

	    eq.EvaluateFunction(x_init, fvec);
	    //DISP(x_init); DISP(fvec);
	    if (Norm2(fvec) < 1e-12)
	      {
		roots.PushBack(Complex_wp(x_init(0), x_init(1)));
	      }
	  }
    }

  DISP(roots); DISP(roots.GetM());
  VectComplex_wp roots0 = roots;
  
  ofstream file_out("racine.dat");
  file_out.precision(15);
  int nb_iteration = 0;
  bool decreasing_step = true;
  while (z < zmax)
    {
      nb_iteration++;
      if (nb_iteration%100 == 0)
	cout << "z = " << z << endl;
      
      eq.z = z; eq_m.z = z;
      file_out << z;
      for (int k = 0; k < roots.GetM(); k++)
	{
	  x_init(0) = realpart(roots(k));
	  x_init(1) = imagpart(roots(k));
	  
	  solver.Init(eq, x_init); nb_iter = 0;
	  solver.Newton(eq, x_init, nb_iter, norme_ref, true, false);

	  roots(k) = Complex_wp(x_init(0), x_init(1));
	  file_out << " " << x_init(0) << " " << x_init(1);
	}

      file_out << '\n';
      
      //solver.Init(eq_m, x_init_m); nb_iter = 0;
      //solver.Newton(eq_m, x_init_m, nb_iter, norme_ref, true, false);
      
      if (decreasing_step)
	z -= dz;
      else
	z += dz;
      
      if (z < zmin)
	{
	  decreasing_step = false;
	  z = z0 + dz;
	  roots = roots0;
	}
    }
  
  file_out.close();
  
  /*
  z0 = 100.0; z1 = 100.0;
  z = z0;
  eq.z = z; eq_m.z = z;
 
  for (int j = -5; j <= 5; j++)
    {
      x_init(0) = 2.0*pi_wp*j / z0;
      x_init(1) = 0.02;
      
      solver.Init(eq, x_init); nb_iter = 0;
      solver.Newton(eq, x_init, nb_iter, norme_ref, true, false);

      eq.EvaluateFunction(x_init, fvec);
      DISP(j); DISP(x_init); DISP(fvec); DISP(z);
    }
  */
}

void DisplayMovie()
{
  VectReal_wp x, y;
  Linspace(-pi_wp/8, pi_wp/8, 1600, x);
  Linspace(Real_wp(-0.4), Real_wp(0.2), 800, y);
  
  VectReal_wp z, log_z;
  //Linspace(10.0, 100.0, 2, z);
  Linspace(Real_wp(1.0), Real_wp(5.0), 21, log_z);
  z.Reallocate(log_z.GetM());
  for (int k = 0; k < log_z.GetM(); k++)
    z(k) = pow(Real_wp(10), log_z(k));
  
  VectReal_wp pt(2), f(2);
  int m = x.GetM(), n = y.GetM(); DISP(n); DISP(m);
  VectReal_wp Val(m*n); Vector<int> vgray;
  for (int p = 0; p < z.GetM(); p++)
    {      
      DipoleEquation eq(z(p),true);

      TinyMatrix<Complex_wp, Symmetric, 3, 3> G;
      Real_wp Gmin(300), Gmax(-300);
      for (int i = 0; i < x.GetM(); i++)
	for (int j = 0; j < y.GetM(); j++)
	  {
	    pt(0) = x(i); pt(1) = y(j);
	    eq.EvaluateFunction(pt, f);
	    R3 ptX(Real_wp(0), Real_wp(0), -Real_wp(0.25)*z(p));
	    R3 ptY(Real_wp(0), Real_wp(0), Real_wp(0.25)*z(p));
	    Complex_wp omega = eq.Omega + Complex_wp(x(i), y(j));
	    eq.ComputeGreenKernel(ptX, ptY, omega, G);
	    Val(m*j + i) = log10(abs(G(0, 0)));
	    if (Val(m*j + i) > Gmax)
	      Gmax = Val(m*j + i);

	    if (Val(m*j + i) < Gmin)
	      Gmin = Val(m*j + i);

	    exit(0);
	    //Val(m*j + i) = 1.0-log10(Norm2(f));
	  }

      DISP(p); DISP(z(p)); DISP(Gmin); DISP(Gmax);
      WriteJpeg(Val, "test" + NumberToString(p) + ".jpg", n, n, m, 0, 0, vgray,
		Gmin, Gmax, ColorMapEnum::HOT);
    }
}



int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  cout.precision(16);
  
  //FollowRoots();

  DisplayMovie();

  return FinalizeMontjoie();
}
