#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  if (argc != 5)
    {
      cout << "Entrez nx, ny, nz et le nombre d'iterations " << endl;
      return -1;
    }
  
  InitMontjoie(argc, argv);

  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int nz = atoi(argv[3]);
  int nb_iter = atoi(argv[4]);

  /* Vector<Real_wp> x_t(nx); Vector<Complex_wp> y_hat(nx/2+1);
  x_t.FillRand(); x_t *= 1e-9; x_t.Write("x_time.dat");
  FftRealInterface fft_time;
  fft_time.Init(nx);
  
  DISP(x_t);
  fft_time.ApplyForward(x_t, y_hat);
  DISP(y_hat);
  x_t.Fill(0);
  fft_time.ApplyInverse(y_hat, x_t);
  
  DISP(x_t); */

  FftInterface<Complex_wp> fft;
  fft.SelectFftAlgorithm(fft.MANUAL);
  if (nz == 1)
    {
      if (ny == 1)
        fft.Init(nx);
      else
        fft.Init(nx, ny);
    }
  else
    fft.Init(nx, ny, nz);
  
  
  srand(0);
  Vector<complex<double> > x(nx*ny*nz);
    
    
  for (int k = 0; k < x.GetM(); k++)
    x(k) = complex<double>(rand(), rand())/double(RAND_MAX);

  glob_chrono.Start(VirtualTimer::ALL);
  
  Complex_wp vloc = 0;
  for (int i = 0; i < nb_iter; i++)
    {  
      //x.Write("x0_t.dat");
      fft.ApplyForward(x); vloc += x(0);
      //x.Write("x1_t.dat");

      //Vector<TinyVector<complex<double>, 1> > x2(nx*ny*nz);
      //for (int j = 0; j < x2.GetM(); j++)
      //x2(j)(0) = x(j);

      fft.ApplyInverse(x);
      vloc += x(0);
      //x.Write("x2_t.dat");  
      
      /* TinyVector<Complex_wp, 1> val;
      if (nz == 1)
        {
          if (ny == 1)
            {
              for (int ix = 0; ix < nx; ix++)
                {
                  fft.ApplyInversePoint(ix, x2, val);
                  x(ix) = val(0);
                }
            }
          else
            {
              for (int ix = 0; ix < nx; ix++)
                for (int iy = 0; iy < ny; iy++)
                  {
                    fft.ApplyInversePoint(ix, iy, x2, val);
                    x(ix*ny+iy) = val(0);
                  } 
            }
        }
      else
        {
          for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
              for (int iz = 0; iz < nz; iz++)
                {
                  fft.ApplyInversePoint(ix, iy, iz, x2, val);
                  x(nz*(ix*ny+iy)+iz) = val(0);
                }
        }
      
        x.Write("x3.dat");            */
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  
  DISP(vloc);
  DISP(glob_chrono.GetSeconds(VirtualTimer::ALL));
  
  return FinalizeMontjoie();
}
