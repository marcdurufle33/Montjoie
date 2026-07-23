#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Output/MontjoieOutput.hxx"

#include "Output/OutputOpenCV.cxx"

using namespace Montjoie;

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);

  // ici on fournit la transposee de u pour avoir une image correspondant a ce qu' on veut
  int m = 400, n = 400;
  VectReal_wp Val(m*n);

  VectReal_wp x, y;
  Linspace(-10.0, 350.0, m, x);
  Linspace(-10.0, 350.0, n, y);

  Mesh<Dimension2> mesh;
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      Val(m*j + i) = sin(0.02*x(i))*sin(0.02*y(j));

  mesh.Read("joly2d_view.mesh");
  WriteContour(mesh, x(0), x(m-1), y(0), y(n-1), m, n, Val);
  
  // example avec un petit cadre blanc autour
  WriteJpeg(Val, "toto.jpg", n, n, m, 0, 0, -1.0, 1.0, ColorMapEnum::JET);
  
  return FinalizeMontjoie();
}  

