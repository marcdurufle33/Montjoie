#include "Harmonic/MontjoieHarmonic.hxx"

using namespace Montjoie;

Real_wp threshold = 1e-12;

template<class T>
void GetRand(T& x)
{
  x = T(rand()) / RAND_MAX;
}

template<class T, int p>
void GetRand(TinyVector<T, p>& x)
{
  for (int i = 0; i < p; i++)
    x(i) = T(rand()) / RAND_MAX;
}

template<class T, class Dimension>
void CheckIncidentField(const IncidentWaveField<T, Dimension>& u)
{
  int nb_points = 100;
  TinyVector<Real_wp, Dimension::dim_N> x, xm, xp;
  TinyVector<T, Dimension::dim_N> df0;
  Real_wp h = pow(epsilon_machine, Real_wp(1)/3)/3;
  Real_wp epsilon = 100*h*h;
  T f, f0, fm, fp, df_num;
  // checking gradients of u
  for (int k = 0; k < nb_points; k++)
    {
      GetRand(x);
      
      u.EvaluateFunction(x, f);      
      u.EvaluateFunctionGradient(x, f0, df0);
      if ((abs(f - f0) > threshold) || isnan(abs(f - f0)))
	{
	  cout << "EvaluateFunctionGradient incorrect" << endl;
	  abort();
	}
      
      for (int i = 0; i < Dimension::dim_N; i++)
	{
	  xm = x; xm(i) -= h; xp = x; xp(i) += h;
	  u.EvaluateFunction(xm, fm);
	  u.EvaluateFunction(xp, fp);
	  
	  df_num = (fp - fm) / (2*h);
	  if ((abs(df_num - df0(i)) > epsilon) || isnan(abs(df_num - df0(i))))
	    {
	      DISP(i); DISP(df_num); DISP(df0); DISP(abs(df_num - df0(i)));
	      cout << "EvaluateFunctionGradient incorrect" << endl;
	      abort();
	    }
	}
    }
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  threshold = 1e4*epsilon_machine;
  
  R2 kwave; R3 kwave3D;
  GetRand(kwave); GetRand(kwave3D);

  {
    HankelIncidentField<Complex_wp, Dimension2> u(R2(1.4, 0.8), kwave);
    CheckIncidentField(u);
  }

  {
    HankelIncidentField<Complex_wp, Dimension3> u(R3(1.4, 0.8, 0.2), kwave3D);
    CheckIncidentField(u);
  }

  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
