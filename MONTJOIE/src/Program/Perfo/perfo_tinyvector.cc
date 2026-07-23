#define MONTJOIE_WITH_TWO_DIM

#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

template<int m>
void RunTestMlt(TinyVector<double, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre les vecteurs Seldon
  int jeton = glob_chrono.GetNumber();
  Vector<double> X0(m), Y0(m), X(m), Y(m);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      X0.FillRand(); Y0.FillRand();
      X.FillRand(); Y.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	{
	  Mlt(1.056, X0); Mlt(1.012, Y0); Mlt(1.033, X); Mlt(1.032, Y);
	  Mlt(0.556, X0); Mlt(0.512, Y0); Mlt(0.443, X); Mlt(0.432, Y);
	  Mlt(0.756, X0); Mlt(0.672, Y0); Mlt(0.783, X); Mlt(0.712, Y);
	  Mlt(1.156, X0); Mlt(1.042, Y0); Mlt(1.213, X); Mlt(1.092, Y);
	}
    }
  
  glob_chrono.Stop(jeton);
  cout<<"Big vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
    
  TinyVector<double, m> y, x0, y0;
  glob_chrono.Reset(jeton);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      x0.FillRand(); y0.FillRand();
      x.FillRand(); y.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	{
	  Mlt(1.056, x0); Mlt(1.012, y0); Mlt(1.033, x); Mlt(1.032, y);
	  Mlt(0.556, x0); Mlt(0.512, y0); Mlt(0.443, x); Mlt(0.432, y);
	  Mlt(0.756, x0); Mlt(0.672, y0); Mlt(0.783, x); Mlt(0.712, y);
	  Mlt(1.156, x0); Mlt(1.042, y0); Mlt(1.113, x); Mlt(1.092, y);
	}
    }

  glob_chrono.Stop(jeton);
  cout<<"Small vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
  
}


template<int m>
void RunTestAdd(TinyVector<double, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre les vecteurs Seldon
  int jeton = glob_chrono.GetNumber();
  Vector<double> X0(m), Y0(m), X(m), Y(m), Z(m);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      X0.FillRand(); Y0.FillRand();
      X.FillRand(); Y.FillRand(); Z.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	{
	  Add(2.356, Z, X0); Add(2.52, Z, Y0); Add(2.43, X, Z); Add(2.32, Y, Z);
	  Add(0.556, Z, X0); Add(0.512, Z, Y0); Add(0.443, X, Z); Add(0.432, Y, Z);
	  Add(0.756, Z, X0); Add(0.672, Z, Y0); Add(0.783, X, Z); Add(0.712, Y, Z);
	  Add(1.156, Z, X0); Add(1.042, Z, Y0); Add(1.213, X, Z); Add(1.092, Y, Z);
	}
    }
  
  glob_chrono.Stop(jeton);
  cout<<"Big vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
    
  TinyVector<double, m> y, x0, y0, z;
  glob_chrono.Reset(jeton);
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      x0.FillRand(); y0.FillRand();
      x.FillRand(); y.FillRand(); z.FillRand();
      for (int i = 0; i < nb_inner_iter; i++)
	{
	  Add(2.356, z, x0); Add(2.52, z, y0); Add(2.43, x, z); Add(2.32, y, z);
	  Add(0.556, z, x0); Add(0.512, z, y0); Add(0.443, x, z); Add(0.432, y, z);
	  Add(0.756, z, x0); Add(0.672, z, y0); Add(0.783, x, z); Add(0.712, y, z);
	  Add(1.156, z, x0); Add(1.042, z, y0); Add(1.113, x, z); Add(1.092, y, z);
	}
    }

  glob_chrono.Stop(jeton);
  cout<<"Small vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
  
}


template<int m>
void RunTestMltAdd(TinyVector<double, m>& x, int nb_inner_iter, int nb_outer_iter)
{
  // on chronometre le produit matrice vecteur de Seldon
  int jeton = glob_chrono.GetNumber();
  Matrix<double> A(m, m);
  Vector<double> X(m), Y(m); double alpha, beta(0), beta_;
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      A.FillRand(); Y.FillRand();
      Mlt(1e-9, A); Mlt(1e-9, Y);
      alpha = A(0,0)*1e-8; beta_ = Y(0);
      for (int i = 0; i < nb_inner_iter; i++)
	{
	  X.FillRand();
	  MltAdd(alpha, A, X, beta_, Y);
	}
      beta += Norm2(Y);
    }
  cout<<beta<<endl;
  
  glob_chrono.Stop(jeton);
  cout<<"Big vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
    
  TinyMatrix<double, General, m, m> a;
  TinyVector<double, m> y;
  glob_chrono.Reset(jeton); beta = 0;
  glob_chrono.Start(jeton);
  for (int n = 0; n < nb_outer_iter; n++)
    {
      a.FillRand(); y.FillRand();
      Mlt(1e-9, a); Mlt(1e-9, y);
      alpha = a(0,0)*1e-8; beta_ = Y(0);
      for (int i = 0; i < nb_inner_iter; i++)
	{
	  x.FillRand();
	  Mlt(beta_, y);
	  MltAdd(alpha, a, x, y);
	}
      beta += Norm2(y);
    }
  cout<<beta<<endl;

  glob_chrono.Stop(jeton);
  cout<<"Small vectors took "<<glob_chrono.GetSeconds(jeton)<<endl;
  
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  cout.precision(7);
  
  if (argc < 3)
    {
      cout<<"Enter more arguments"<<endl;
      return -1;
    }
  
  int nb_inner_iter = atoi(argv[1]);

  {
    TinyVector<complex<double>, 10> x, y, z;
    
    for (int i = 0; i < 10; i++)
      {
	//x(i, 0) = std::rand() / double(RAND_MAX);
	//y(i, 0) = std::rand() / double(RAND_MAX);
	
	x(i) = complex<double>(std::rand(), std::rand()) / double(RAND_MAX);
	y(i) = complex<double>(std::rand(), std::rand()) / double(RAND_MAX);
      }
    
    cout << "x = " << x << endl;
    
    complex<double> sum(0, 0);
    for (int i = 0; i < nb_inner_iter; i++)
      {
	z = complex<double>(2, -0.4)*x - 3.0*y;
	sum += z(std::rand()%5);
	z = complex<double>(2, -0.4)*x - 3.0*y + z;
	sum += z(std::rand()%5);
	z = complex<double>(2, -0.4)*x - 3.0*y + z;
	sum += z(std::rand()%5);
	z = complex<double>(2, -0.4)*x - 3.0*y + z;
	sum += z(std::rand()%5);
	z = complex<double>(2, -0.4)*x - 3.0*y + z;
	sum += z(std::rand()%5);
      }

    DISP(sum);
    
    return 0;
  }
  

  int nb_outer_iter = atoi(argv[2]);
  TinyVector<double, 10> x;
  RunTestMltAdd(x, nb_inner_iter, nb_outer_iter);
  
}
