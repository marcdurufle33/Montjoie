#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

void LaunchSimu(int nb_iter, int taille)
{
  Matrix<Real_wp> A(taille, taille);
  Real_wp a = 0; IVect pivot;
  for (int i = 0; i < nb_iter; i++)
    {
      A.FillRand();
      GetLU(A, pivot);
      a += MaxAbs(A);
    }
  
  DISP(a);
}

void CheckTimer(VirtualTimer& timer, int N, int nb_iter)
{
  timer.Start(VirtualTimer::ALL);
  
  LaunchSimu(nb_iter, N);

  timer.Stop(VirtualTimer::ALL);
  cout << "temps : " << timer.GetSeconds(VirtualTimer::ALL) << endl;
  
  timer.Reset(VirtualTimer::ALL);
  cout << "temps : " << timer.GetSeconds(VirtualTimer::ALL) << endl;
  
  int number = timer.GetNumber();
  timer.Start(number);

  LaunchSimu(nb_iter, 2*N);

  timer.Stop(number);
  cout << "temps : " << timer.GetSeconds(number) << endl;
  
  timer.ReleaseNumber(number);
  int n2 = timer.GetNumber();
  if (n2 != number)
    {
      cout << "Release number incorrect" << endl;
      abort();
    }

  timer.ReserveNumber(0);
  timer.Reset(0);
  timer.Start(0);

  LaunchSimu(nb_iter, int(1.5*N));

  timer.Stop(0);
  cout << "temps : " << timer.GetSeconds(0) << endl;
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc != 3)
    {
      cout << "Usage : ./timer.x n nb_iter" << endl;
      abort();
    }

  int nb_iter = atoi(argv[2]);
  int N = atoi(argv[1]);
  
  {
    cout << "Testing basic timer" << endl;
    BasicTimer timer;
    CheckTimer(timer, N, nb_iter);
  }  

  {
    cout << "Testing accurate timer" << endl;
    AccurateTimer timer;
    CheckTimer(timer, N, nb_iter);
  }  

  {
    cout << "Testing real timer" << endl;
    RealTimer timer;
    CheckTimer(timer, N, nb_iter);
  }  

  {
    MontjoieTimer timer;
    timer.SetTimer(MontjoieTimer::ACCURATE_TIMER);
    
    timer.SetMessage("Factorize", "factorize the matrix");
    timer.SetMessage("Mlt", "perform the matrix-matrix product");

    int taille = 2*N;
    Matrix<Real_wp> A(taille, taille), B(taille, taille), C(taille, taille);
    Real_wp a = 0; IVect pivot;
    for (int i = 0; i < nb_iter; i++)
      {
	timer.Start("Factorize");

	A.FillRand();
	GetLU(A, pivot);

	timer.Stop("Factorize");
	timer.Start("Mlt");
	
	Mlt(A, B, C);
	a += MaxAbs(C);

	timer.Stop("Mlt");
      }
    
    DISP(a);
    
    timer.DisplayAll();
  }
  
  return FinalizeMontjoie();
}
