#include "Montjoie.hxx"

using namespace Montjoie;

Real_wp fac(int m)
{
  return tgamma(Real_wp(m+1));
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 4)
    {
      cout << "Entrez l'ordre d'approximation et NN" << endl;
      abort();
    }

  int m = atoi(argv[1]); // nombre d'etapes
  int NN(atof(argv[2])); // calcul du dvt de Taylor de la dispersion
  int type_scheme = atoi(argv[3]); // -1 pour Pade, et nombre d'etapes additionnelles pour Linear-Sdirk

  UnivariatePolynomial<Real_wp> P, Q;
  if (type_scheme == -1)
    {
      // schemas de Pade
      P.SetOrder(m); Q.SetOrder(m); Real_wp one(1);
      for (int i = 0; i <= m; i ++)
	{
	  Real_wp coef = fac(m)*fac(2*m-i) / (fac(2*m) * fac(i) * fac(m-i));
	  P(i) = coef;
	}
      
      for(int i = 0; i <= m; i++)
	{
	  Q(i)=fac(2*m-i)*fac(m)/(fac(2*m)*fac(i)*fac(m-i));
	  if(i%2 != 0)
	    Q(i) *= -one; 
	}

      PadeScheme_Iterator<Real_wp> pade;
      
      pade.SetOrder(2*m);
      DISP(pade.GetWeights());
    }
  else
    {
      // schemas Sdirk lineaires
      LinearSdirkScheme_Iterator<Real_wp> sdirk;
      
      sdirk.SetOrder(m+1, type_scheme);
      P = sdirk.GetNumeratorStabilityFunction();
      Q = sdirk.GetDenominatorStabilityFunction();

      DISP(sdirk.GetWeights());
    }
  
  DISP(P); DISP(Q);

  int Np = 1000;
  int nb_sys = m; // ok pour Pade
  if (type_scheme > 0)
    nb_sys += type_scheme;

  // calcul de la dispersion sur Np points    
  ofstream file_out("dispersion.dat"); file_out.precision(15);
  Real_wp offset(0), argRz(0), argRz_prev(0);

  for (int i = 0; i <= Np; i++)
    {
      Real_wp z = Real_wp((i+1)*nb_sys)*pi_wp/(2*Np);
      Complex_wp Rz = P.Evaluate(Iwp*z) / Q.Evaluate(Iwp*z);
      argRz_prev = argRz; argRz = arg(Rz);
      if (argRz < 0)
	argRz += Real_wp(2)*pi_wp;
      
      if ((i > 0) && (abs(argRz -argRz_prev) > pi_wp))
	offset += Real_wp(2)*pi_wp;
      
      Real_wp err = (z - (offset + argRz))/z;
      //DISP(z); DISP(zi); DISP(arg(Rz)); DISP(argRz);
      file_out << z/nb_sys << " " << err << '\n';
    }
  file_out.close();

  // calcul de la dissipation sur Np points
  ofstream file_out2("dissipation.dat"); file_out2.precision(15);
  Real_wp modRz(0);
  
  for (int i = 0; i <= Np; i++)
    {
      Real_wp z = Real_wp((i+1)*nb_sys)*pi_wp/(2*Np);
      Complex_wp Rz = P.Evaluate(Iwp*z) / Q.Evaluate(Iwp*z);
      modRz = abs(Rz);
      
      Real_wp err = modRz - 1;
      //DISP(z); DISP(abs(Rz)); DISP(modRz);
      file_out2 << z/nb_sys << " " << err << '\n';
    }
  
  file_out2.close();



  bool search_taylor_dvt = true;
  if (search_taylor_dvt)
    {
      // searching Taylor expansion
      Np = 2*m+8;
      Matrix<Real_wp> VanDerMonde(Np+1, Np+1);
      VectReal_wp coef_relation(Np+1);
      offset = Real_wp(0); Real_wp z(0), z_prev(0);
      for (int i = 0; i <= Np; i++)
	{
	  z_prev = z; z = Real_wp((i+1)*m)*pi_wp/NN;
	  for (int j = 0; j <= Np; j++)
	    VanDerMonde(i, j) = pow(z, Real_wp(j));

	  Complex_wp Rz = P.Evaluate(Iwp*z) / Q.Evaluate(Iwp*z);
	  argRz_prev = argRz; argRz = arg(Rz);
	  if (argRz < 0)
	    argRz += Real_wp(2)*pi_wp;
	  
	  if ((i > 0) && (abs(argRz -argRz_prev) > pi_wp))
	    offset += Real_wp(2)*pi_wp;
	  
	  Real_wp err = (z - (offset + argRz))/z;	  
	  coef_relation(i) = err;

	  if (i == 1)
	    {
	      Real_wp order = (log(abs(coef_relation(1))) - log(abs(coef_relation(0)))) / (log(z) - log(z_prev));
	      DISP(order); DISP(coef_relation(0)); DISP(coef_relation(1)); DISP(z);
	      order = toInteger(round(order));
	      Real_wp cte = coef_relation(1) / pow(z, order);
	      DISP(cte);
	      Vector<int64_t> decomp_coef;
	      DecomposeContinuedFraction(cte, Real_wp(1e-15), decomp_coef);
	      int64_t num, denom;
	      GetNumeratorDenominator(decomp_coef, num, denom);
	      DISP(num); DISP(denom);
	    }
	}
      
      IVect ipivot_vdm(VanDerMonde.GetM()); GetLU(VanDerMonde, ipivot_vdm);
      SolveLU(VanDerMonde, ipivot_vdm, coef_relation);
      DISP(coef_relation);

      for (int i = 0; i <= Np; i++)
	{
	  if (abs(coef_relation(i)) > 1e-20)
	    {
	      DISP(i); DISP(coef_relation(i));
	      Vector<int64_t> decomp_coef;
	      DecomposeContinuedFraction(coef_relation(i), Real_wp(1e-20), decomp_coef);
	      int64_t num, denom;
	      GetNumeratorDenominator(decomp_coef, num, denom);
	      DISP(num); DISP(denom);
	    }
	}
    }
  
  
  return FinalizeMontjoie();
}
