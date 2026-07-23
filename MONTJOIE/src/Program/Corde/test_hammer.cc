#define MONTJOIE_WITH_ONE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Corde/MontjoieString.hxx"
#include "Instationary/MontjoieTime.hxx"

using namespace Montjoie;

template<int nb_base, int nb_quad, class TypeEquation>
class RkHammerScheme : public VirtualOdeSystem<Real_wp>
{
public :
  typedef Montjoie::Real_wp Real_wp;
  Shank<nb_base, nb_quad, TypeEquation>& var;
  Real_wp Ki, Ri, mH, tau;
  bool presence_tau;
  
  RkHammerScheme(Shank<nb_base, nb_quad, TypeEquation>& sys) : var(sys)
  {
    mH = var.hammer.masse_marteau;
    Ki = var.hammer.Ki(0);
    Ri = var.hammer.Ri(0);
    tau = 1e-5;
    presence_tau = true;
  }
  
  void EvaluateFunction(const Real_wp& tn, const VectReal_wp& xi, VectReal_wp& F)
  {
    F(0) = xi(1);
    Real_wp d = var.y_wall - xi(0);
    Real_wp phi = var.hammer.Phi(d);
    Real_wp phi_prime = var.hammer.PhiPrime(d);
    
    if (presence_tau)
      {
        F(0) = xi(1);
        F(1) = xi(2);
        F(2) = -1.0/(tau*mH)*(Ki*phi - Ri*xi(1)*phi_prime) - 1.0/tau*xi(2);
      }
    else
      F(1) = -(Ki*phi - Ri*xi(1)*phi_prime)/mH;
  }

  void SetDirichletCondition(const Real_wp& t, int nb_deriv, VectReal_wp& x, Real_wp alpha = 1)
  {}
  
  void GetForceDistance(const VectReal_wp& xi, Real_wp& Force, Real_wp& ecrase)
  {
    Real_wp d = var.y_wall - xi(0);
    ecrase = var.hammer.delta - d;
    if (ecrase < 0)
      ecrase = 0.0;
    
    Real_wp phi = var.hammer.Phi(d);
    Real_wp phi_prime = var.hammer.PhiPrime(d);
    Force = (Ki*phi - Ri*xi(1)*phi_prime);
  }
  
};

template<int nb_base, int nb_quad, class TypeEquation>
void RunAll(Shank<nb_base, nb_quad, TypeEquation>& var, const string& input_file)
{
  typedef Montjoie::Real_wp Real_wp;
  ReadInputFile(input_file, var);
  
  // var.ConstructAll(0);
  
  int nb_iterations = toInteger(var.t_end/var.Deltat);
  DISP(nb_iterations);
  
  RungeKutta_Iterator<Real_wp> rk;
  rk.SetOrder(4);
  
  RkHammerScheme<nb_base, nb_quad, TypeEquation> sys(var);
  VectReal_wp xsi(3); xsi.Fill(0);
  xsi(0) = var.hammer.eloignement_marteau; xsi(1) = var.hammer.vitesse_marteau;
  rk.SetInitialCondition(0.0, var.Deltat, xsi, sys);
  
  ofstream file_out("xsi.dat");
  file_out.precision(15);
  Real_wp Force, ecrase;
  for (int nt = 0; nt < nb_iterations; nt++)
    {
      Real_wp t = var.Deltat*nt;
      var.t_courant = t;

      rk.Advance(t, nt, sys);
      
      sys.GetForceDistance(rk.GetIterate(), Force, ecrase);
      
      file_out << t << " " << rk.GetIterate()(0) << " " << Force << " " << ecrase << endl;
    }
  
  file_out.close();
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  srand(0);
  
  if (argc < 2)
    {
      cout << "entrez le nom du fichier de donnees" << endl;
      abort();
    }
  
  Hammer hammer;
  Shank<5, 5, WaveEquationStiffString> var(hammer);
  string input_file(argv[1]);
  
  RunAll(var, input_file);
  
  return FinalizeMontjoie();
}
