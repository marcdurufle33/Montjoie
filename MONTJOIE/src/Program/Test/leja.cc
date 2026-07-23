#include "Montjoie.hxx"

using namespace Montjoie;

VectReal_wp ComputeLejaPoints(const Real_wp& z0, int N)
{
  Real_wp a(0), b(1);
  VectReal_wp pos(1);
  pos(0) = z0;
  VectReal_wp Lr, Li;
  for (int i = 0; i < N; i++)
    {
      UnivariatePolynomial<Real_wp> prod, dprod, Pol;
      Pol.SetOrder(1); Pol(1) = 1.0;
      prod.SetOrder(0); prod(0) = Real_wp(1);
      for (int j = 0; j < pos.GetM(); j++)
        {
          Pol(0) = -pos(j);
          prod = prod*Pol;
        }

      Real_wp new_pos = a;
      Real_wp val_min = abs(prod.Evaluate(a));

      Real_wp val = abs(prod.Evaluate(b));
      if (val > val_min)
        {
          val_min = val;
          new_pos = b;
        }
      
      DerivatePolynomial(prod, dprod);
      SolvePolynomialEquation(dprod, Lr, Li);
      for (int j = 0; j < Lr.GetM(); j++)
        if (abs(Li(j)) <= 1e-12)
          if ((Lr(j) <= b) && (Lr(j) >= a))
            {
              val = abs(prod.Evaluate(Lr(j)));
              if (val > val_min)
                {
                  val_min = val;
                  new_pos = Lr(j);
                }
            }

      pos.PushBack(new_pos);
      //DISP(i); DISP(new_pos); DISP(prod);
      //prod.WriteText("prod"+to_str(i)+".dat");
    }

  return pos;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  int N = atoi(argv[1]);
  VectReal_wp pos = ComputeLejaPoints(Real_wp(0.5), N);
  DISP(pos);
  pos.WriteText("leja.dat");
  
  return FinalizeMontjoie();
}
