#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

// Legendre polynomes
class LegendrePolynome
{
public :
  LegendrePolynome(int n) {}
  
  
};

void ComputePn(int n, const LegendrePolynome& Pn, double x, DVect& Pn0)
{
  // allocate if needed
  Pn0.Reallocate(n);
  Pn0(0) = 1.0; Pn0(1) = x;
  for (int i = 2; i < n ; i++)
    Pn0(i) = ( x*Pn0(i-1)*(2*i-1) - Pn0(i-2)*(i-1))/i;
    
  for (int i = 0; i < n; i++)
    Pn0(i) *= sqrt((2.0*i+1)/2.0);
}

void ComputeSinus(int n, double x, DVect& Pn0)
{
  // allocate if needed
  Pn0.Reallocate(2*n+1);
  Pn0(0) = 0.5*sqrt(2);
  for (int i = 1; i <= n; i++)
    {
      Pn0(i) = cos(pi_wp*i*x);
      Pn0(i+n) = sin(pi_wp*i*x);
    }
}

double GetValue_ProfilM(double x)
{
  double M = 0;
  double xp = abs(x);
  M = sqrt(1+xp)-1;
  if (x < 0)
    M = -M;
  // return x;
  // return sin(pi_wp*x/2);
  return M;
}

int SINUS = 0, POLYNOME = 1;

int main(int argc,char** argv)
{
  InitMontjoie(argc, argv);
  
  cout.precision(16);
  if (argc < 3)
    {
      cout<<"Entrez deux arguments "<<endl;
      return -1;
    }

  int order_approx = atoi(argv[1]);
  int type_discretisation = atoi(argv[2]);
  int nb_elt = order_approx;
  if (type_discretisation == SINUS)
    nb_elt = 2*order_approx + 1;
  
  Matrix<double> Mh(nb_elt, nb_elt); Mh.Fill(0);
  LegendrePolynome Pn(nb_elt); DVect Pn0(nb_elt);
  int order = nb_elt+1;
  if (type_discretisation == SINUS)
    order = 5*order_approx;
  
  Globatto<Real_wp> gauss; gauss.ConstructQuadrature(order);
  cout<<"calcul de la matrice "<<endl;
  for (int k = 0; k < gauss.GetNbPointsQuad(); k++)
    {
      double x = 2.0*gauss.Points(k) - 1.0;
      double wn = 2.0*gauss.Weights(k);
      wn *= GetValue_ProfilM(x);
      if (type_discretisation == SINUS)
	ComputeSinus(order_approx, x, Pn0);
      else
	ComputePn(order_approx, Pn, x, Pn0);
      
      for (int i = 0; i < nb_elt; i++)
	for (int j = 0; j < nb_elt; j++)	
	  Mh(i,j) += wn*Pn0(i)*Pn0(j);
      
    }
  
  // DISP(Mh);
  
  Matrix<double> Ah(2*nb_elt, 2*nb_elt); Ah.Fill(0);
  for (int i = 0; i < nb_elt; i++)
    for (int j = 0; j < nb_elt; j++)
      {
	// part Mh
	Ah(i,j) = Mh(i,j);
	Ah(i+nb_elt,j+nb_elt) = Mh(i,j);
	
	// part Ih
	if (i == j)
	  Ah(i,j+nb_elt) = 1.0;
	
      }
  
  // part Eh
  Ah(nb_elt,0) = 1.0;
  
  Mh.Clear();
  int N = Ah.GetM();
  Vector<double> lambda_real(N), lambda_imag(N);
  Matrix<double> eigen_vector(N,N);
  cout<<"calcul des valeurs et vecteurs propres"<<endl;
  GetEigenvaluesEigenvectors(Ah, lambda_real, lambda_imag, eigen_vector);
  lambda_real.WriteText("Lh_r.dat"); lambda_imag.WriteText("Lh_i.dat");
  double Lmax = 0.0;
  for (int i = 0; i < lambda_real.GetM(); i++)
    Lmax = max(Lmax, lambda_real(i));
  
  DISP(Lmax); DISP(Lmax-sqrt(2));

  return FinalizeMontjoie();
}
