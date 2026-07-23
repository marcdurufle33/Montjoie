#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

Real_wp GetRand()
{
  return Real_wp(rand())/RAND_MAX;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  SplineInterpolation<Complex_wp> spline;

  Matrix<Real_wp> Hz;
  Hz.ReadText("Hz_test.dat");
  int N = Hz.GetM();
  Real_wp xmin = Hz(0, 0), xmax = Hz(N-1, 0);
  for (int nb_points = 5; nb_points < 20; nb_points++)
    {
      VectReal_wp xdiv(nb_points);
      VectComplex_wp y(nb_points);
      xdiv(0) = xmin; xdiv(nb_points-1) = xmax;
      y(0) = Complex_wp(log(-Hz(0, 1)), log(-Hz(0, 2)));
      y(nb_points-1) = Complex_wp(log(-Hz(N-1, 1)), log(-Hz(N-1, 2)));

      Real_wp dx = (xmax-xmin)/(nb_points-1);
      for (int i = 1; i < nb_points-1; i++)
        {
          Real_wp x = xmin + i*dx;
          int j = 0;
          while (Hz(j, 0) < x)
            j++;
            
          xdiv(i) = Hz(j, 0);
          y(i) = Complex_wp(log(-Hz(j, 1)), log(-Hz(j, 2)));
        }

      DISP(xdiv); DISP(y);
      spline.Init(xdiv, y);
      Real_wp err_max = 0;
      Matrix<Real_wp> Hz_interp(N, 3);
      Hz_interp(0, 0) = exp(Hz(0, 0));
      Complex_wp vloc = exp(Hz(0, 0)) / Complex_wp(Hz(0, 1), Hz(0, 2));
      Hz_interp(0, 1) = real(vloc); Hz_interp(0, 2) = imag(vloc);
      vloc = exp(Hz(N-1, 0)) / Complex_wp(Hz(N-1, 1), Hz(N-1, 2));
      Hz_interp(N-1, 0) = exp(Hz(N-1, 0));
      Hz_interp(N-1, 1) = real(vloc); Hz_interp(N-1, 2) = imag(vloc);
            
      for (int i = 1; i < N-1; i++)
        {
          Complex_wp y_i = spline.Evaluate(Hz(i, 0));
          Complex_wp val = Complex_wp(-exp(real(y_i)), -exp(imag(y_i)));
          Real_wp ym = exp(Hz(i, 0));
          Complex_wp Hz_i = ym / val;
          Hz_interp(i, 0) = ym;
          Hz_interp(i, 1) = real(Hz_i); Hz_interp(i, 2) = imag(Hz_i);
          Complex_wp Hz_ref = ym / Complex_wp(Hz(i, 1), Hz(i, 2));
          Real_wp err = abs(Hz_i - Hz_ref) / abs(Hz_ref);
          if (err > err_max)
            err_max = err;
        }

      Hz_interp.WriteText("Hz_N" + to_str(nb_points) + ".dat");
      DISP(nb_points); DISP(err_max);
      
    }
 
  return FinalizeMontjoie();
}
