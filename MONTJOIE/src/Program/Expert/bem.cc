#include "Montjoie.hxx"

#include "Harmonic/BoundaryIntegralEquation.cxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc != 2)
    {
      cout << "Entrez un fichier de donnees a utiliser" << endl;
      cout << "Utilisation : ./bem.x data.ini" << endl;
      abort();
    }

  string input_file(argv[1]);

  // we get the type of element selected by the user, and type of equation
  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  // running the simulation
  BoundaryIntegralEquation_Maxwell3D var;
  var.RunAll(input_file, type_element);

  return FinalizeMontjoie();
}
