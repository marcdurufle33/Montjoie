#include "FiniteElement/MontjoieFiniteElement.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  EdgeHierarchic edge;
  edge.ConstructFiniteElement(10, 10, 10, Globatto<Real_wp>::QUADRATURE_GAUSS,
			      EdgeHierarchic::SHEN);

  // testing derivatives
  Real_wp x = 0.23, h = 1e-6;
  VectReal_wp phi_p, phi_m, der_phi;
  edge.ComputeValuesPhiRef(x+h, phi_p);
  edge.ComputeValuesPhiRef(x-h, phi_m);
  edge.ComputeGradientPhiRef(x, der_phi);

  for (int i = 0; i < edge.GetNbDof(); i++)
    {
      Real_wp val_ex = der_phi(i);
      Real_wp val_num = (phi_p(i) - phi_m(i)) / (2*h);
      DISP(i); DISP(val_ex); DISP(val_num);
      if (abs(val_ex - val_num) > h)
	{
	  cout << "ComputeGradientPhiRef incorrect" << endl;
	  abort();
	}
    }
  
  return FinalizeMontjoie();
}
