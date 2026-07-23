#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout << "Fournissez deux noms de fichiers " << endl;
      abort();
    }
  
  bool ascii = true;
  //int double_prec = OutputTypeEnum::QUADRUPLE_PRECISION;
  VectReal_wp val, val_ref;
  GridInterpolationFull<Dimension2> grid;
    
  ReadMatlab(val_ref, grid, string(argv[1]), ascii);
  ReadMatlab(val, grid, string(argv[2]), ascii);
  
  Real_wp norm_val_ref = Norm2(val_ref);
  for (int i = 0; i < val.GetM(); i++)
    val(i) -= val_ref(i);
  
  Real_wp norm_diff = Norm2(val);
  cout << "Erreur L^2 = " << norm_diff/norm_val_ref << endl;
  
  return FinalizeMontjoie();
}
