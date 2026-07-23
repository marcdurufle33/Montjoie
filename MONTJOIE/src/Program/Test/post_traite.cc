#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Entrez un nom de fichier a post-traiter" << endl;
      abort();
    }
  
  string file_name(argv[1]);
  Matrix<double> val;
  val.ReadText(file_name);
  
  int N = val.GetM();
  ofstream file_rho("rho.don"), file_mu("mu.don"), file_sigma("sigma.don");
  ofstream file_radius("radius.don");
  file_rho.precision(15); file_mu.precision(15); file_sigma.precision(15);
  file_radius.precision(15);
  double omega = 2.1e6*2.0*M_PI;
  for (int i = 0; i < N; i++) 
    {
      double c = val(N-1-i, 1)/100;
      double r = val(N-1-i, 0);
      double rho = val(N-1-i, 2)*1e3;
      double sigma = (omega/100)/(rho*c*c);
      
      file_rho << r << " " << 1.0/(rho*c*c) << '\n';
      file_mu << r << " " << 1.0/rho << '\n';
      file_sigma << r << " " << sigma << '\n';
      file_radius << r << '\n';
    }
  
  file_radius.close();
  file_rho.close();
  file_mu.close();
  file_sigma.close();
}
