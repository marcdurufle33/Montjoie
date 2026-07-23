#define MONTJOIE_WITH_ONE_DIM

#include "Harmonic/MontjoieHarmonic.hxx"

#include "Elliptic/Maxwell/NonLinearOpticsProblem.hxx"
#include "Elliptic/Maxwell/NonLinearOpticsProblem.cxx"

using namespace Montjoie;


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  RamanEffect var;
  
  Real_wp tau1 = 12.2e-15;
  Real_wp tau2 = 32.0e-15;
  Real_wp Tmax_raman = 25.0*tau2;
  Real_wp alpha = 1.0;
  var.Init(tau1, tau2, Tmax_raman, alpha);
  
  Real_wp dt = 1.0e-16;
  var.Construct(2, dt);
  
  Real_wp tau = 70.0e-15/sqrt(2.0);
  Real_wp tau_impulse = tau / (2.0*sqrt(log(2.0)));
  DISP(tau_impulse);
  
  Real_wp Tmax = 800.0e-15;
  Real_wp t = 600.001e-15;
  int nb_iter = toInteger(round(t/dt));
  DISP(nb_iter);
  
  Real_wp lambda0 = 1.55e-6;
  Real_wp omega = 2*pi_wp*PhysicalConstant::speed_light / lambda0;

  Real_wp val_ref = 0;
  Globatto<Real_wp> gauss;
  gauss.ConstructQuadrature(3000);
  
  int nb_val = 16001; VectReal_wp time_interval(nb_val);
  Real_wp tmin = -8e-13;
  VectReal_wp val_anal(nb_val);
  for (int k = 0; k < nb_val; k++)
    val_anal(k) = var.EvaluateH(k*dt);
  
  val_anal.Write("raman_ref.dat");

  for (int k = 0; k < nb_val; k++)
    {
      t = Tmax + tmin + k*dt;
      val_ref = 0;
      for (int n = 0; n < gauss.GetNbPointsQuad(); n++)
	{
	  Real_wp s = gauss.Points(n)*t;
	  Real_wp t_n = t-s;
	  Real_wp pulse = 2e10*exp(-0.5*square((t_n-Tmax)/tau_impulse));
	  //pulse *= cos(omega*t_n);
	  //pulse *= abs(Complex_wp(0.8*cos(omega*t_n), 0.2*cos(omega*t_n+pi_wp/5)));
	  
	  Real_wp h = var.EvaluateH(s);
	  val_ref += h*pulse*pulse*gauss.Weights(n)*t*sin(2.0*omega*t_n);
	}
      
      if (abs(t - Tmax) < 1e-20)
	DISP(val_ref);

      val_anal(k) = val_ref;
    }
  
  val_anal.Write("conv_anal.dat");

  VectReal_wp Edeux(nb_val);
  VectComplex_wp E(2);
  for (int n = 0; n < nb_val; n++)
    {
      Real_wp t_n = n*dt; time_interval(n) = t_n;
      Real_wp pulse = 2e10*exp(-0.5*square((t_n-Tmax)/tau_impulse));
      E(0) = pulse*Complex_wp(cos(omega*t_n), sin(omega*t_n));
      //E(0) = pulse*Complex_wp(cos(omega*t_n), 0);
      E(1) = pulse*Complex_wp(0.8*cos(omega*t_n), 0.2*cos(omega*t_n+pi_wp/5));
      var.StoreNewField(t_n, E);
      
      var.InitCoefficient(t_n);
      Real_wp gamma = 1, eps = 0, invEps = 0;
      var.UpdateCoefficient(t_n, 0, gamma, eps, invEps);
      
      if (n%100 == 0)
	DISP(eps);
      
      Edeux(n) = pulse*pulse;
      val_anal(n) = eps;
    }
  
  Edeux.Write("Enum.dat");
  time_interval.Write("time_num.dat");
  val_anal.Write("conv_num.dat");
  exit(0);
  //var.moduleE(0).Write("moduleE0.dat");
  //var.moduleE(1).Write("moduleE1.dat");
  
  // now evaluating Raman effect for t slightly different
  var.InitCoefficient(t);
  //var.hEval.Write("evalH.dat");
  //var.level_romberg = 15;
  Real_wp gamma = 0.1, eps = 1.5, invEps = 1.0/eps;
  var.UpdateCoefficient(t, 0, gamma, eps, invEps);
  DISP(gamma); DISP(eps); DISP(invEps);

  gamma = 0.1; eps = 1.5; invEps = 1.0/eps;
  var.UpdateCoefficient(t, 1, gamma, eps, invEps);
  DISP(gamma); DISP(eps); DISP(invEps);
  
  return FinalizeMontjoie();
}

