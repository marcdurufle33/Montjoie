#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout << "Provide two file names where matrices are stored" << endl;
      cout << "Usage : ./diff_matrix Ah.dat Bh.dat" << endl;
      abort();
    }
  
  string nameA(argv[1]), nameB(argv[2]);
  Matrix<Complex_wp, General, ArrayRowSparse> A, B;
  A.ReadText(nameA);
  B.ReadText(nameB);
  
  bool trans = false;
  if (argc == 4)
    if (string(argv[3]) == "t")
      trans = true;

  if (trans)
    Transpose(A);
  
  Real_wp normeA = Norm1(A);
  Real_wp normeB = Norm1(B);
  
  Complex_wp one(1, 0);
  Add(-one, A, B);
  
  Real_wp normeDiff = Norm1(B);
  Real_wp diff_max = 0; int imax = -1, jmax = -1;
  for (int i = 0; i < B.GetM(); i++)
    for (int j = 0; j < B.GetRowSize(i); j++)
      if (abs(B.Value(i, j)) > diff_max)
        {
          imax = i; jmax = B.Index(i, j);
          diff_max = abs(B.Value(i, j));
        }

  cout << "Maximum reached for i, j = " << imax << " " << jmax << endl;
  if ((imax >= 0) && (jmax >= 0))
    cout << " || A - B ||_inf = " << abs(B(imax, jmax)) << " A = " << A(imax, jmax) << " B = " << B(imax, jmax) + A(imax, jmax) << endl;
  
  cout << "norm(A) = " << normeA << endl;
  cout << "norm(B) = " << normeB << endl;
  cout << "norm(A-B) = " << normeDiff << endl;
  cout << "norm(A-B)/norm(A) = " << normeDiff/normeA << endl;
  
  return FinalizeMontjoie();
}
