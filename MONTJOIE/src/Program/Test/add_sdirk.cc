#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  if (argc != 4)
    {
      cout << "Provide the order, zmax and gamma" << endl;
      cout << "Usage : ./test.x 10 2.5 0.25" << endl;
      abort();
    }

  InitMontjoie(argc, argv);
  
  int order = atoi(argv[1]);  
  Real_wp zmax = atoi(argv[2]);
  Real_wp gamma = to_num<Real_wp>(string(argv[3]));
  
  int nb_points = 10*order;
  VectReal_wp z;

  Linspace(Real_wp(0), zmax, nb_points, z);
  DISP(z);

  // si vrai, on fixe l'asymptote a un
  bool asymptot_one = true;

  // number of free coefficients
  int nb_terms = order+1;
  if (asymptot_one)
    nb_terms--;
  
  Matrix<Real_wp> Vdm; VectReal_wp y, tau;
  Vdm.Reallocate(2*nb_points, nb_terms);
  y.Reallocate(2*nb_points);

  Real_wp last_coef = pow(-gamma, order);
  for (int i = 0; i < nb_points; i++)
    {
      Complex_wp zi = Iwp*z(i);
      Complex_wp rhs = pow(Real_wp(1) - gamma*zi, order)*exp(zi);
      Real_wp poids = Real_wp(1) / abs(pow(Real_wp(1) - gamma*zi, order));
      rhs *= poids;
      
      // real part
      y(2*i) = realpart(rhs);
      for (int j = 0; j < nb_terms; j += 2)
	{
	  Vdm(2*i, j) = pow(z(i), j)*poids;
	  if (j%4 == 2)
	    Vdm(2*i, j) = -Vdm(2*i, j);
	}

      if ((asymptot_one) && (nb_terms%2 == 0))
	{
	  Real_wp coef_vdm = pow(z(i), nb_terms)*poids;
	  if (nb_terms%4 == 2)
	    coef_vdm = -coef_vdm;
	  
	  y(2*i) -= coef_vdm*last_coef;
	}

      // imaginary part
      y(2*i+1) = imagpart(rhs);
      for (int j = 1; j < nb_terms; j += 2)
	{
	  Vdm(2*i+1, j) = pow(z(i), j)*poids;
	  if (j%4 == 3)
	    Vdm(2*i+1, j) = -Vdm(2*i+1, j);
	}

      if ((asymptot_one) && (nb_terms%2 == 1))
	{
	  Real_wp coef_vdm = pow(z(i), nb_terms)*poids;
	  if (nb_terms%4 == 3)
	    coef_vdm = -coef_vdm;
	  
	  y(2*i+1) -= coef_vdm*last_coef;
	}

    }

  GetQR(Vdm, tau);
  SolveQR(Vdm, tau, y);

  if (asymptot_one)
    y.PushBack(last_coef);
  
  cout << "y = ";
  for (int j = order; j >= 0; j--)
    cout << y(j) << ", ";
  
  cout << endl;

  cout << "Asymptote = " << y(order) / pow(gamma, order) << endl;

  VectComplex_wp Y(nb_points), Yref(nb_points);
  for (int i = 0; i < nb_points; i++)
    {
      Complex_wp zi = Iwp*z(i);
      Complex_wp val_ref = exp(zi);

      // Horner algorithm
      Complex_wp val = y(order);
      Complex_wp denom; SetComplexOne(denom);
      for (int j = order-1; j >= 0; j--)
	{
	  val = val*zi + y(j);
	  denom *= Real_wp(1) - gamma*zi;
	}

      Y(i) = val / denom;
      Yref(i) = val_ref;
    }

  cout << "Erreur L^2 = " << Norm2(Y - Yref) / Norm2(Yref) << endl; 

  ofstream file_out("exp_app.dat"); file_out.precision(15);
  Linspace(Real_wp(0), Real_wp(20)*zmax, 20*nb_points, z);
  for (int i = 0; i < 20*nb_points; i++)
    {
      Complex_wp zi = Iwp*z(i);
      Complex_wp val_ref = exp(zi);

      // Horner algorithm
      Complex_wp val = y(order);
      Complex_wp denom; SetComplexOne(denom);
      for (int j = order-1; j >= 0; j--)
	{
	  val = val*zi + y(j);
	  denom *= Real_wp(1) - gamma*zi;
	}

      file_out << z(i) << " " <<  realpart(val / denom) << " " << imagpart(val/denom) << '\n';
    }
  
  file_out.close();

  return FinalizeMontjoie();
}

