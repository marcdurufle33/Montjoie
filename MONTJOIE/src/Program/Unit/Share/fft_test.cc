#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

Real_wp threshold = 1e-11;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    {
      return false;
    }

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > threshold) || isnan(x(i)-y(i)))
      {
	DISP(i); DISP(x(i)); DISP(y(i));
	return false;
      }

  return true;
}

void CheckFft()
{
  int n = 200;
  VectComplex_wp x, y, x2, y2, z, xchap(n/2+1), xchap2(n/2+1);

  GenerateRandomVector(x, n);

  FftInterface<Complex_wp> fft;
  FftRealInterface fft_real;

  fft.SelectFftAlgorithm(fft.MANUAL);

  fft.Init(n);

  if (fft.GetNbPoints() != n)
    {
      cout << "GetNbPoints incorrect" << endl;
      abort();
    }

  y = x;
  fft.ApplyForward(y);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(x, x2))
    {
      cout << "ApplyForward/ApplyBackward incorrect" << endl;
      abort();
    }

  if (fft.GetMemorySize() < 3*sizeof(Complex_wp)*n)
    {
      cout << "GetMemorySize incorrect" << endl;
      abort();
    }

  fft_real.Init(n);

  if (fft_real.GetNbPoints() != n)
    {
      cout << "GetNbPoints incorrect" << endl;
      abort();
    }

  if (fft_real.GetMemorySize() < 3*sizeof(Real_wp)*n)
    {
      cout << "GetMemorySize incorrect" << endl;
      abort();
    }

  Vector<Real_wp> x_real, x_real2(n);

  GenerateRandomVector(x_real, n);

  fft_real.ApplyForward(x_real, xchap);
  z.Reallocate(n);
  for (int i = 0; i < n; i++)
    z(i) = Complex_wp(x_real(i), 0);
  
  fft.ApplyForward(z);

  for (int i = 0; i <= n/2; i++)
    if ((abs(xchap(i) - z(i)) > threshold) || isnan(abs(xchap(i) - z(i))))
      {
	cout << "ApplyForward (real) incorrect " << endl;
	abort();
      }

  fft_real.ApplyInverse(xchap, x_real2);
  if (!EqualVector(x_real, x_real2))
    {
      cout << "ApplyInverse (real) incorrect" << endl;
      abort();
    }

#ifdef MONTJOIE_WITH_GSL
  fft.SelectFftAlgorithm(fft.FFT_GSL);

  fft.Init(n);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (GSL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (GSL) incorrect" << endl;
      abort();
    }

  fft_real.SelectFftAlgorithm(fft.FFT_GSL);

  fft_real.Init(n);

  fft_real.ApplyForward(x_real, xchap2);
  fft_real.ApplyInverse(xchap, x_real2);

  if (!EqualVector(xchap, xchap2))
    {
      cout << "ApplyForward (GSL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x_real, x_real2))
    {
      cout << "ApplyBackward (GSL) incorrect" << endl;
      abort();
    }

#endif

#ifdef MONTJOIE_WITH_FFTW
  fft.SetNbThreads(1);
  fft.SelectFftAlgorithm(fft.FFTW);

  fft.Init(n);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (FFTW) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (FFTW) incorrect" << endl;
      abort();
    }

  fft_real.SelectFftAlgorithm(fft.FFTW);

  fft_real.Init(n);

  fft_real.ApplyForward(x_real, xchap2);
  fft_real.ApplyInverse(xchap, x_real2);

  if (!EqualVector(xchap, xchap2))
    {
      cout << "ApplyForward (FFTW) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x_real, x_real2))
    {
      cout << "ApplyBackward (FFTW) incorrect" << endl;
      abort();
    }

#endif

#ifdef SELDON_WITH_MKL
  fft.SelectFftAlgorithm(fft.FFT_MKL);

  fft.Init(n);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (MKL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (MKL) incorrect" << endl;
      abort();
    }

  fft_real.SelectFftAlgorithm(fft.FFT_MKL);

  fft_real.Init(n);

  fft_real.ApplyForward(x_real, xchap2);
  fft_real.ApplyInverse(xchap, x_real2);

  if (!EqualVector(xchap, xchap2))
    {
      cout << "ApplyForward (MKL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x_real, x_real2))
    {
      cout << "ApplyBackward (MKL) incorrect" << endl;
      abort();
    }

#endif

}

void CheckFft2D()
{
  int nx = 20, ny = 30;
  VectComplex_wp x, y, x2, y2;

  GenerateRandomVector(x, nx*ny);

  FftInterface<Complex_wp> fft;

  fft.SelectFftAlgorithm(fft.MANUAL);

  fft.Init(nx, ny);
  
  /*if (fft.GetNbPoints() != nx*ny)
    {
      cout << "GetNbPoints incorrect" << endl;
      abort();
      }*/

  y = x;
  fft.ApplyForward(y);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(x, x2))
    {
      cout << "ApplyForward/ApplyBackward incorrect" << endl;
      abort();
    }

  if (fft.GetMemorySize() < 3*sizeof(Complex_wp)*nx*ny)
    {
      cout << "GetMemorySize incorrect" << endl;
      abort();
    }

#ifdef MONTJOIE_WITH_GSL
  fft.SelectFftAlgorithm(fft.FFT_GSL);

  fft.Init(nx, ny);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (GSL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (GSL) incorrect" << endl;
      abort();
    }

#endif

#ifdef MONTJOIE_WITH_FFTW
  fft.SetNbThreads(1);
  fft.SelectFftAlgorithm(fft.FFTW);

  fft.Init(nx, ny);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (FFTW) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (FFTW) incorrect" << endl;
      abort();
    }

#endif

#ifdef SELDON_WITH_MKL
  fft.SelectFftAlgorithm(fft.FFT_MKL);

  fft.Init(nx, ny);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (MKL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (MKL) incorrect" << endl;
      abort();
    }

#endif

}


void CheckFft3D()
{
  int nx = 20, ny = 30, nz = 16;
  VectComplex_wp x, y, x2, y2;

  GenerateRandomVector(x, nx*ny*nz);

  FftInterface<Complex_wp> fft;

  fft.SelectFftAlgorithm(fft.MANUAL);

  fft.Init(nx, ny, nz);

  /*if (fft.GetNbPoints() != nx*ny*nz)
    {
      cout << "GetNbPoints incorrect" << endl;
      abort();
      }*/

  y = x;
  fft.ApplyForward(y);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(x, x2))
    {
      cout << "ApplyForward/ApplyBackward incorrect" << endl;
      abort();
    }

  if (fft.GetMemorySize() < 3*sizeof(Complex_wp)*nx*ny*nz)
    {
      cout << "GetMemorySize incorrect" << endl;
      abort();
    }

#ifdef MONTJOIE_WITH_GSL
  fft.SelectFftAlgorithm(fft.FFT_GSL);

  fft.Init(nx, ny, nz);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (GSL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (GSL) incorrect" << endl;
      abort();
    }

#endif

#ifdef MONTJOIE_WITH_FFTW
  fft.SetNbThreads(1);
  fft.SelectFftAlgorithm(fft.FFTW);

  fft.Init(nx, ny, nz);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (FFTW) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (FFTW) incorrect" << endl;
      abort();
    }

#endif

#ifdef SELDON_WITH_MKL
  fft.SelectFftAlgorithm(fft.FFT_MKL);

  fft.Init(nx, ny, nz);

  y2 = x;
  fft.ApplyForward(y2);

  x2 = y;
  fft.ApplyInverse(x2);

  if (!EqualVector(y, y2))
    {
      cout << "ApplyForward (MKL) incorrect" << endl;
      abort();
    }

  if (!EqualVector(x, x2))
    {
      cout << "ApplyBackward (MKL) incorrect" << endl;
      abort();
    }

#endif

}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  CheckFft();
  CheckFft2D();
  CheckFft3D();
  
  Real_wp a = -10.0, b = 10.0;
  int N = 200;
  
  VectReal_wp t;
  Linspace(a, b, N, t);
  
  VectReal_wp f(N), frev(N), foriginal;
  VectComplex_wp fchap(N/2+1);
  for (int i = 0; i < N; i++)
    f(i) = exp(-t(i)*t(i))*cos(2.0*pi_wp*t(i));

  foriginal = f;
  VectReal_wp omega(N/2+1);
  for (int i = 0; i < N/2+1; i++)
    omega(i) = 2.0*pi_wp*Real_wp(i)/(b-a);    
  
  FftRealInterface fft;
  fft.Init(N);
  
  fft.ApplyForward(f, fchap);
  
  //  fchap.Write("EvalFchap.dat");
  
  //for (int i = 0; i < N/2+1; i++)
  //fchap(i) *= -Iwp*omega(i);
  
  fft.ApplyInverse(fchap, frev);
  
  for (int i = 0; i < N; i++)
    if (abs(foriginal(i) - frev(i)) > threshold)
      {
        cout << "Fft incorrect" << endl;
        DISP(i); DISP(foriginal(i)); DISP(frev(i));
        abort();
      }
  
  // foriginal.Write("EvalF.dat");
  // frev.Write("EvalFder.dat");

  cout << "All tests passed successfully" << endl;
    
  return FinalizeMontjoie();
}
