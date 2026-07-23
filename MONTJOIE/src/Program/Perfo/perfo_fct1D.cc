#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 4)
    {
      cout << "Entrez l'ordre d'approximation le nombre de points et d'iteration" << endl;
      abort();
    }
  
  int r = atoi(argv[1]);
  int nb_pts = atoi(argv[2]);
  int nb_iter = atoi(argv[3]);
  
  {
    Globatto<Real_wp> lob;
    lob.ConstructQuadrature(r);
    
    MontjoieTimer chrono;
    chrono.SetMessage("EvaluateBasisFunctions", "evaluate basis functions");
    chrono.Start("EvaluateBasisFunctions");
    Vector<Real_wp> x(nb_pts);
    Vector<Real_wp> phi(r+1);
    for (int k = 0; k < nb_iter; k++)
      {
        x.FillRand();
        Mlt(1.0/RAND_MAX, x);
        
        for (int j = 0; j < x.GetM(); j++)
          lob.ComputeValuesPhiRef(x(j), phi);        
      }
    
    chrono.Stop("EvaluateBasisFunctions");

    chrono.SetMessage("EvaluateBasisFunctionsPhi1D", "evaluate basis functions with phi1D");
    chrono.Start("EvaluateBasisFunctionsPhi1D");
    for (int k = 0; k < nb_iter; k++)
      {
        x.FillRand();
        Mlt(1.0/RAND_MAX, x);
        
        for (int j = 0; j < x.GetM(); j++)
          for (int i = 0; i <= r; i++)
            phi(i) = lob.EvaluatePhi(i, x(j));
      }
    
    chrono.Stop("EvaluateBasisFunctionsPhi1D");
    
    chrono.DisplayAll();
  }


  
#ifdef MONTJOIE_WITH_MPFR
  if (false)
    {
      Globatto<double> lob;
      Globatto<Real_wp> lob_exact;
      
      lob.ConstructQuadrature(r, lob.QUADRATURE_GAUSS);
      lob_exact.ConstructQuadrature(r, lob.QUADRATURE_GAUSS);
      
      double err = 0;
      for (int i = 0; i <= r; i++)
        err = max(err, abs(lob.Points(i) - double(lob_exact.Points(i))));
      
      cout << "Erreur sur les points de Gauss : " << err << endl;
      
      // matrice de derivation
      Matrix<double> val;
      Matrix<Real_wp> val_exact;
      lob.ComputeGradPhi();  
      lob_exact.ComputeGradPhi();
      val_exact = lob_exact.GradPhi();
      
      err = 0;
      val = lob.GradPhi();  
      for (int i = 0; i <= r; i++)
        for (int j = 0; j <= r; j++)
          err = max(err, abs(val(i, j) - val_exact(i, j)));
      
      cout << "Erreur sur la matrice de derivation " << err << endl;
      
      int N = 101; Real_wp delta = 1e-13;
      Vector<Real_wp> x((r+2)*N);
      
      // interval [0,1]
      for (int i = 0; i < N; i++)
        x(i) = Real_wp(i)/(N-1);
      
      // intervals close to x_i
      Vector<Real_wp> xdiv;
      Linspace(-delta, delta, N, xdiv);
      for (int i = 0; i <= r; i++)
        for (int j = 0; j < N; j++)
          x((i+1)*N+j) = lob_exact.Points(i) + xdiv(j);
      
      Vector<Real_wp> f_exact(r+1);
      Vector<double> f(r+1);
      double errPhi(0), errFct(0);
      for (int i = 0; i < x.GetM(); i++)
        {
          for (int j = 0; j <= r; j++)
            f_exact(j) = lob_exact.EvaluatePhi(j, x(i));
          
          // error with dphi1D
          for (int j = 0; j <= r; j++)
            f(j) = lob.EvaluatePhi(j, double(x(i)));
          
          err = 0;
          for (int j = 0; j <= r; j++)
            err = max(err, abs(double(f_exact(j)) - f(j)));
          
          errPhi = max(errPhi, err);
          
          // error with ComputeValuesPhiRef
          lob.ComputeValuesPhiRef(double(x(i)), f);
          
          err = 0;
          for (int j = 0; j <= r; j++)
            err = max(err, abs(double(f_exact(j)) - f(j)));
          
          errFct = max(errFct, err);
        }
      
      cout << "Maximum error using phi1D = " << errPhi << endl;
      cout << "Maximum error using ComputeValuesPhiRef = " << errFct << endl;
    }
#endif
 
  return FinalizeMontjoie();
}
