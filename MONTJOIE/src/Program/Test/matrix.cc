#include "Montjoie.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  srand(2020);

  /*for (int k = 0; k < 6; k++)
    {
      int n = rand()%1281648;
      cout << "Numero matrice passage = " << n << endl;
      
      for (int i = 0; i < 3; i++)
        {
          Real_wp L = -4 + rand()%9;
          DISP(i); DISP(L);
        }
    }
  */
  
  TinyMatrix<double, General, 3, 3> J1, J2, J3, J4, J5, J6;
  J1(0, 0) = 4.0; J1(1, 1) = -2.0; J1(2, 2) = -2.0; J1(1, 2) = 1.0;
  J2(0, 0) = -1.0; J2(1, 1) = 1.0; J2(2, 2) = 2.0;
  J3(0, 0) = -2.0; J3(1, 1) = 3.0; J3(2, 2) = 4.0;
  J4(0, 0) = -2.0; J4(1, 1) = -1.0; J4(2, 2) = -2.0;
  J5(0, 0) = -3.0; J5(1, 1) = 1.0; J5(2, 2) = 1.0; J5(1, 2) = 1.0;
  J6(0, 0) = -2.0; J6(1, 1) = -2.0; J6(2, 2) = -2.0; J6(0, 1) = 1.0; J6(1, 2) = 1.0;
  
  TinyMatrix<double, General, 3, 3> A, B, C, R;

  int Nmax = 3; int p = 0;
  int Ncomb = pow(2*Nmax + 1, 9);
  ofstream file_out("totoR.dat");
  for (int i = 1; i <= Ncomb; i++)
    {
      int n = i; bool valid_P = true;
      for (int j = 0; j < 3; j++)
        {
          int nb_zero = 0;
          for (int k = 0; k < 3; k++)
            {
              A(k, j) = -Nmax + n%(2*Nmax+1);
              n = n / (2*Nmax+1);
              if (A(k, j) == Real_wp(0))
                nb_zero++;
            }

          if (nb_zero >= 2)
            valid_P = false;
        }

      Real_wp det = Det(A);
      if ((abs(det) == 1) && valid_P)
        {
          //file_out << "array([[" << A(0, 0) << ", " << A(0, 1) << ", " << A(0, 2) << "],["
          //         << A(1, 0) << ", " << A(1, 1) << ", " << A(1, 2) << "],["
          //         << A(2, 0) << ", " << A(2, 1) << ", " << A(2, 2) << "]])" << endl;
          //DISP(i); DISP(A);
          //exit(0);

          GetInverse(A, B);
          Mlt(A, J2, C); Mlt(C, B, R);
          bool lowR = true;
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
              if (abs(R(j, k)) >= 4)
                lowR = false;

          if (lowR)
            {
              p++;
              if (p == 405)
                {
                  DISP(A); DISP(J2); DISP(R);
                }
              //exit(0);
            }
        }
    }

  DISP(p);
  
  return FinalizeMontjoie();
}
