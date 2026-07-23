#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  int n = to_num<int>(argv[1]);
  int nb_iter = to_num<int>(argv[2]);
  
  Matrix<double> A(n, n);
  
  glob_chrono.Start(VirtualTimer::ALL);
  
  for (int k = 0; k < nb_iter; k++)
    {
      for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	  A.Get(i, j) = double(i+j);
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  DISP(glob_chrono.GetSeconds(VirtualTimer::ALL));
  
  return FinalizeMontjoie();
}
