
//////////////////////////////////////////////////////////////////////
// programme pour la résolution des problèmes cellulaires U10 et U11
//////////////////////////////////////////////////////////////////////



// APPEL DU PROGRAMME
// ./CellulairePeriodique.x FichierDonnee.ini






// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1


#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Elliptic/Helmholtz/CellulairePeriodique.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/CellulairePeriodique.cxx"
#endif

using namespace Montjoie;


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file, const string& name_element)
{
  typedef typename TypeEquation::Complexe Complexe;
  
  var.SetTypeEquation("none");
  // initialisation des indices
  var.InitIndices(PhysicalConstant::nb_max_indices);
  
  // on lit le fichier de donnees
  ReadInputFile(input_file, var);

  All_LinearSolver glob_solver(var);
  ReadInputFile(input_file, glob_solver);
  
  // on construit le maillage et les elements finis
  var.ComputeMeshAndFiniteElement(name_element);

  // autres initialisations
  var.PerformOtherInitializations();
  var.PeriodeX = var.GetXmax() - var.GetXmin();
  var.PeriodeY = var.GetYmax() - var.GetYmin();
  
  // on selectionne le meilleur solveur
  Dimension3 dim; Symmetric property;
  glob_solver.SelectOptimalLinearSolver(var.GetMeshNumbering(0).GetOrder(), var.GetNbDof(), dim, property);
  
  Vector<Complexe> x_sol(var.GetNbDof());
  
  Complexe D11_eps, D22_eps, A_eps, D11_mu, D22_mu, A_mu;
  
  /************************
   * Calculs pour epsilon *
   ************************/
  
  // penalisation, on resout -Penal
  // premier calcul en resolvant div(mu grad p) = - d mu/dx
  // avec mu = epsilon + i sigma/omega
  TinyMatrix<Complexe, Symmetric, 3, 3> mu;
  DISP(PhysicalConstant::adimensionalization);
  DISP(PhysicalConstant::ADIM_YES);
  for (int i = 0; i < var.ref_rho.GetM(); i++)
    {
      Complexe epsilon = var.ref_epsilon_max(i);
      // pour sigma, l'impedance apparait du fait de l'adimensionalisation des equations
      Complexe sigma(var.ref_sigma_max(i));
      if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
        sigma *= PhysicalConstant::impedance0;
      
      Complexe rho = epsilon + Iwp*sigma/var.GetOmega();
      mu.SetIdentity(); mu *= rho;
      var.ref_mu(i).SetConstant(mu);
    }
  
  // calcul de Ji, DFi, etc
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();
  
  // factorisation de la matrice
  GlobalGenericMatrix<Complexe> nat_mat;
  glob_solver.PerformFactorizationStep(nat_mat);
  
  // et on resout les trois problemes -dmu/dx, -dmu/dy, -dmu/dz
  var.type_source_cellule = var.SRC_MU_X;
  var.ComputeRightHandSide(x_sol);
  glob_solver.ComputeSolution(x_sol, nat_mat);
  var.WriteOutputFile(x_sol, 0);
  
  Complexe d11, d22, a3d;
  var.CalculD11D22A3D(x_sol, d11, d22, a3d);
  D11_eps = d11;
  
  // deuxieme probleme
  var.type_source_cellule = var.SRC_MU_Y;
  var.ComputeRightHandSide(x_sol);
  glob_solver.ComputeSolution(x_sol, nat_mat);
  var.WriteOutputFile(x_sol, 1);
  
  var.CalculD11D22A3D(x_sol, d11, d22, a3d);
  D22_eps = d22;
  
  // troisieme probleme
  var.type_source_cellule = var.SRC_ABSORBING;
  var.ComputeRightHandSide(x_sol);
  glob_solver.ComputeSolution(x_sol, nat_mat);
  
  var.CalculD11D22A3D(x_sol, d11, d22, a3d);
  A_eps = a3d;

  
  /*******************
   * Calculs pour mu *
   *******************/
  
  // penalisation, on resout -Penal
  // premier calcul en resolvant div(mu grad p) = - d mu/dx
  // avec mu = mu ...
  for (int i = 0; i < var.ref_rho.GetM(); i++)
    {
      Complexe mu_max = var.ref_mu_max(i);
      mu.SetIdentity(); mu *= mu_max;
      var.ref_mu(i).SetConstant(mu);
    }
  
  // calcul de Ji, DFi, etc
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();
  
  // factorisation de la matrice
  glob_solver.PerformFactorizationStep(nat_mat);
  
  // et on resout les trois problemes -dmu/dx, -dmu/dy, -dmu/dz
  var.type_source_cellule = var.SRC_MU_X;
  var.ComputeRightHandSide(x_sol);
  glob_solver.ComputeSolution(x_sol, nat_mat);
  
  var.CalculD11D22A3D(x_sol, d11, d22, a3d);
  D11_mu = d11;
  
  // deuxieme probleme
  var.type_source_cellule = var.SRC_MU_Y;
  var.ComputeRightHandSide(x_sol);
  glob_solver.ComputeSolution(x_sol, nat_mat);
  
  var.CalculD11D22A3D(x_sol, d11, d22, a3d);
  D22_mu = d22;
  
  // troisieme probleme
  var.type_source_cellule = var.SRC_ABSORBING;
  var.ComputeRightHandSide(x_sol);
  glob_solver.ComputeSolution(x_sol, nat_mat);
  
  var.CalculD11D22A3D(x_sol, d11, d22, a3d);
  A_mu = a3d;
  
  
  /*****************************
   * Impression des constantes *
   *****************************/
  
  
  // on affiche les coefficients obtenus
  cout << " " << endl ;
  cout << " " << endl ;
  cout << " " << endl ;
  cout << " " << endl ;
  cout << "/////////////////////////////////////////////// " << endl ;
  cout << " TRES IMPORTANT " << endl ; 
  DISP(D11_eps) ;
  DISP(D22_eps) ;
  DISP(A_eps) ;
  cout << "/////////////////////////////////////////////// " << endl ;
  cout << " " << endl ;
  cout << " " << endl ;
  cout << " " << endl ;
  cout << " " << endl ;
  // On calcule alpha
  
  double alpha ;
  alpha =0 ;
  if ( real(D11_eps) < 0 )
    alpha = -0.5 * real(D11_eps) ;
  if ( real(D22_eps) < 0 )
    alpha = max ( alpha, -0.5 * real(D22_eps)) ;
  if ( real(A_eps) < 0 )
    alpha = max ( alpha, -0.5 * real(A_eps)) ;
  DISP(alpha) ;
  alpha = max (0.5, alpha) ;
  DISP(alpha) ;
  
  ofstream fichier_constantes (var.GetConstantFileName().data());
  fichier_constantes.precision(12);
  
  fichier_constantes << "alpha = " << alpha << endl ;
  fichier_constantes << "D11 = " << D11_eps << endl ;
  fichier_constantes << "D22 = " << D22_eps << endl ;
  fichier_constantes << "A3D = " << A_eps << endl ;
  fichier_constantes << endl;
  fichier_constantes << "Ligne a inserer dans le fichier de donnees : " << endl;
  fichier_constantes << "ReferenceTransmission = 3 PARAM " << A_eps << " " << A_mu
                     << " " << D11_mu << " " << D22_mu << " " << D11_eps << " " << D22_eps << " " << " alpha delta" << endl;
  fichier_constantes << endl;
  fichier_constantes << "ReferenceTransmission = 3 PARAM " << real(A_eps) << " " << real(A_mu)
                     << " " << real(D11_mu) << " " << real(D22_mu) << " " << real(D11_eps) << " " << real(D22_eps) << " " << " alpha delta" << endl;
  
  fichier_constantes.close();
  
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc,argv);
  
  if (argc > 1)
    {
      // name of the data file
      string file_name_data;
      file_name_data = string(argv[1]);

      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      // instantiation avec le bon element fini
      EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension3> > Vars;
      RunAll(Vars, file_name_data, type_element);
      
      cout<<" we destroy the variables "<<endl;   
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"CellulairePeriodique3D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  cout<<"End of the program"<<endl; 
  
  return 0;
}
