#include "Montjoie.hxx"

using namespace Montjoie;

void cfft(Vector<Complex_wp>& x, VectComplex_wp& phase, int stride)
{
  /*const size_t N = x.GetM();
  if (N <= 1)
    return;
  
  // divide
  VectComplex_wp even(N/2), odd(N/2);
  for (int i = 0; i < N/2; i++)
    {
      even(i) = x(2*i);
      odd(i) = x(2*i+1);
    } 
 
  // conquer
  cfft(even, phase, 2*stride);
  cfft(odd, phase, 2*stride);
 
  // combine
  for (size_t k = 0; k < N/2; ++k)
    {
      Complex_wp t = phase(k*stride) * odd(k);
      x(k    ) = even(k) + t;
      x(k+N/2) = even(k) - t;
      } */

  // DFT
  int N = x.GetM(), k = N, n;
  double thetaT = pi_wp / N;
  Complex_wp phiT = Complex_wp(cos(thetaT), -sin(thetaT)), T;
  DISP(phiT);
  while (k > 1)
    {
      n = k;
      k >>= 1;
      phiT = phiT * phiT;
      T = 1.0;
      if (n == N)
        {
          for (int l = 0; l < k; l++)
            {
              int b = l + k;
              Complex_wp t = x(l) - x(b);
              x(l) += x(b);
              x(b) = t*T;
              T *= phiT;
            }                    
        }
       /*else if (n == N/2)
        {
          for (int l = 0; l < k; l++)
            {
              int a = l, b = l + k;
              Complex_wp t = x(a) - x(b);
              x(a) += x(b);
              x(b) = t*T;

              a += n; b += n;
              t = x(a) - x(b);
              x(a) += x(b);
              x(b) = t*T;
              T *= phiT;
            }
          
          } */
      else
        {
          for (int l = 0; l < k; l++)
            {
              for (int a = l; a < N; a += n)
                {
                  int b = a + k;
                  Complex_wp t = x(a) - x(b);
                  x(a) += x(b);
                  x(b) = t * T;
                }
              T *= phiT;
            }
        }
    }

  DISP(x);
  
  // Decimate
  int m = (int)log2(N);
  for (unsigned int a = 0; a < N; a++)
    {
      unsigned int b = a;
      // Reverse bits
      b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
      b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
      b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
      b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
      b = ((b >> 16) | (b << 16)) >> (32 - m);
      DISP(a); DISP(b);
      if (b > a)
        {
          Complex_wp t = x(a);
          x(a) = x(b);
          x(b) = t;
        }
    }
}

VectComplex_wp ComputePhaseFft(int N)
{
  VectComplex_wp phase(N); Complex_wp t;
  SetComplexOne(t);
  Complex_wp phi = std::polar(1.0, -2 * pi_wp * Real_wp(1) / N);
  for (int i = 0; i < N; i++)
    {
      phase(i) = t;
      t *= phi;
    }
  
  return phase;
}

void cfft2d(VectComplex_wp& u, int Nx, int Ny)
{
  VectComplex_wp phase_x = ComputePhaseFft(Nx);
  VectComplex_wp phase_y = ComputePhaseFft(Ny);
    
  VectComplex_wp ux(Nx);
  VectComplex_wp uy(Ny);
  for (int i = 0; i < Nx; i++)
    {
      int offset = i*Ny;
      for (int j = 0; j < Ny; j++)
        uy(j) = u(offset + j);
        
      cfft(uy, phase_y, 1);
        
      for (int j = 0; j < Ny; j++)
        u(offset + j) = uy(j);        
    }
  
  for (int j = 0; j < Ny; j++)
    {
      for (int i = 0; i < Nx; i++)
        ux(i) = u(i*Ny + j);
      
      cfft(ux, phase_x, 1);

      for (int i = 0; i < Nx; i++)
        u(i*Ny + j) = ux(i);        
    }
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  if (argc < 2)
    {
      cout << "Incorrect number of arguments" << endl;
      abort();
    }
  
  srand(0);

  int N = atoi(argv[1]); int Nx = N;
  int Ny = atoi(argv[2]), Nz = atoi(argv[3]);
  /*  VectComplex_wp z(4), phase;
  z(0) = Complex_wp(1.2, 2.5); z(1) = Complex_wp(1.3, -0.9);
  z(2) = Complex_wp(0.6, 1.4); z(3) = Complex_wp(-0.7, -1.1);

  cfft(z, phase, 1); DISP(z);

  z(0) = Complex_wp(1.2, 2.5); z(1) = Complex_wp(1.3, -0.9);
  z(2) = Complex_wp(0.6, 1.4); z(3) = Complex_wp(-0.7, -1.1); */

  /*
  VectComplex_wp z1, z2, z0(N), z1i, z2i;
  for (int i = 0; i < z0.GetM(); i++)
    z0(i) = Complex_wp(rand(), rand()) / Real_wp(RAND_MAX);

  z1 = z0;
  FftInterface<Complex_wp> fft;
  fft.SelectFftAlgorithm(fft.MANUAL);
  fft.Init(N);
  fft.ApplyForward(z1);

  z1i = z0;
  fft.ApplyInverse(z1i);
  
  DISP(z1);   DISP(z1i);

  z2 = z0;
  fft.SelectFftAlgorithm(fft.FFTW);
  fft.Init(N);
  fft.ApplyForward(z2);

  z2i = z0;
  fft.ApplyInverse(z2i);
      
  DISP(z2); DISP(z2i);
  DISP(Norm2(z1-z2));
  DISP(Norm2(z1i-z2i)); */

  VectComplex_wp z1, z2, z0(Nx*Ny*Nz), z1i, z2i;
  for (int i = 0; i < z0.GetM(); i++)
    z0(i) = Complex_wp(rand(), rand()) / Real_wp(RAND_MAX);

  DISP(z0);
  z1 = z0;
  FftInterface<Complex_wp> fft;
  fft.SelectFftAlgorithm(fft.MANUAL);
  fft.Init(Nx, Ny, Nz);
  fft.ApplyForward(z1);

  z1i = z0;
  fft.ApplyInverse(z1i);
  
  DISP(z1);   DISP(z1i);

  z2 = z0;
  fft.SelectFftAlgorithm(fft.FFTW);
  fft.Init(Nx, Ny, Nz);
  fft.ApplyForward(z2);

  z2i = z0;
  fft.ApplyInverse(z2i);
      
  DISP(z2); DISP(z2i);
  DISP(Norm2(z1-z2));
  DISP(Norm2(z1i-z2i));  
  
  return FinalizeMontjoie();
  
  
  //FftInterface<Complex_wp> fft;
  
  VectComplex_wp H(N*N);
  for (int i = 0; i < H.GetM(); i++)
    H(i) = Complex_wp(rand(), rand()) / Real_wp(RAND_MAX);
  
  glob_chrono.Start("FftComputation");
  
  //DISP(H);

  cfft2d(H, N, N);
  //fft.Init(N, N);
  //fft.ApplyForward(H);
  //DISP(H);

  glob_chrono.Stop("FftComputation");
  cout << "Time = " << glob_chrono.GetSeconds("FftComputation") << endl;
  
  FinalizeMontjoie();
}
