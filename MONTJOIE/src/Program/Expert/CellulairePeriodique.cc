
//////////////////////////////////////////////////////////////////////
// programme pour la résolution des problèmes cellulaires U10 et U11
//////////////////////////////////////////////////////////////////////



// APPEL DU PROGRAMME
// ./CellulairePeriodique.x FichierDonnees.ini








// #define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Elliptic/Helmholtz/CellulairePeriodique.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/CellulairePeriodique.cxx"
#endif

using namespace Montjoie;

template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file,
	    const string& name_element)
{
  typedef typename TypeEquation::Complexe Complexe;
  
  // initialisation des indices
  var.InitIndices(PhysicalConstant::nb_max_indices);
  
  // on lit le fichier de donnees
  ReadInputFile(input_file, var);

  All_LinearSolver glob_solver(var);
  ReadInputFile(input_file, glob_solver);
  
  var.var_chrono.SetMessage("MeshGeneration", "construct and number the mesh");
  var.var_chrono.Start("MeshGeneration");
  
  // on construit le maillage et les elements finis
  var.ComputeMeshAndFiniteElement(name_element);

  // autres initialisations
  var.PerformOtherInitializations();
  var.PeriodeX = var.GetXmax() - var.GetXmin();
  
  var.var_chrono.Stop("MeshGeneration");
  
  // on selectionne le meilleur solveur
  Dimension2 dim; Symmetric property;
  glob_solver.SelectOptimalLinearSolver(var.GetMeshNumbering(0).GetOrder(), var.GetNbDof(), dim, property);
  
  Vector<Complexe> b_rhs(var.GetNbDof()), x_sol(var.GetNbDof());
  
  Complexe A_mu, B1_mu, B2_mu;
  
  /*******************
   * Calculs pour mu *
   *******************/
  
  TinyMatrix<Complexe, Symmetric, 2, 2> mu; Complexe b1, b2, a;
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
  GlobalGenericMatrix<Complexe> nat_mat;
  glob_solver.PerformFactorizationStep(nat_mat);
  
  if (!var.CasDirichlet)
    {
      // premier probleme
      var.type_source_cellule = var.SRC_MU;
      var.ComputeRightHandSide(b_rhs);
      x_sol = b_rhs;
      glob_solver.ComputeSolution(x_sol, nat_mat);
      //var.WriteOutputFile(x_sol, "SolCelluleP1.dat");
      var.WriteOutputFile(x_sol, 0);
      
      var.CalculB1B2(x_sol, b1, b2);
      var.CalculA(x_sol, a);
      B1_mu = b1; B2_mu = b2;
      
      DISP(var.SymmetricCase);
      // deuxieme probleme
      // exchanging Neumann and Dirichlet conditions
      if (var.SymmetricCase)
        {
          for (int i = 1; i <= var.mesh.GetNbReferences(); i++)
            {
              int cond = var.mesh.GetBoundaryCondition(i);
              if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
                var.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_NEUMANN);
              else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
                var.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_DIRICHLET);
            }
          
          var.TreatDirichletCondition();
          
          cout << "factorisation second systeme " << endl;
          glob_solver.PerformFactorizationStep(nat_mat);
        }
      
      var.type_source_cellule = var.SRC_ABSORBING;
      var.ComputeRightHandSide(b_rhs);
      x_sol = b_rhs;
      glob_solver.ComputeSolution(x_sol, nat_mat);
      //var.WriteOutputFile(x_sol, "SolCelluleP2.dat");
      var.WriteOutputFile(x_sol, 1);
      
      var.CalculB1B2(x_sol, b1, b2);
      var.CalculA(x_sol, a);
      A_mu = a;
    }
  else
    {
      // premier probleme
      var.source_valeur_absolue = false;
      var.type_source_cellule = var.SRC_ABSORBING;
      var.ComputeRightHandSide(b_rhs);
      x_sol = b_rhs;
      glob_solver.ComputeSolution(x_sol, nat_mat);
      var.WriteOutputFile(x_sol, 0);
      
      var.CalculA(x_sol, a, true);
      Complexe Sautldirichlet = a;
      DISP(Sautldirichlet) ;
      
      var.CalculA(x_sol, a);
      Complexe Moyenneldirichlet = a;
      DISP(Moyenneldirichlet) ;
      
      // deuxieme probleme
      var.source_valeur_absolue = true;
      var.type_source_cellule = var.SRC_ABSORBING;
      var.ComputeRightHandSide(b_rhs);
      x_sol = b_rhs;
      glob_solver.ComputeSolution(x_sol, nat_mat);
      var.WriteOutputFile(x_sol, 1);
      
      var.CalculA(x_sol, a, true);
      Sautldirichlet = a;
      DISP(Sautldirichlet) ;
      
      var.CalculA(x_sol, a);
      Moyenneldirichlet = a;
      DISP(Moyenneldirichlet) ;
    }
  
  ofstream fichier_constantes (var.GetConstantFileName().data());
  fichier_constantes.precision(12);
  
  fichier_constantes << "A = " << A_mu << endl ;
  fichier_constantes << "B1 = " << B1_mu << endl ;
  fichier_constantes << "B2 = " << B2_mu << endl ;

  fichier_constantes << "Ligne a inserer dans le fichier de donnees : " << endl;
  fichier_constantes << "ReferenceTransmission = 3 PARAM " << A_mu
                     << " " << B1_mu << " " << B2_mu << " alpha delta" << endl;
  
  fichier_constantes << endl;
  fichier_constantes << "ReferenceTransmission = 3 PARAM " << real(A_mu)
                     << " " << real(B1_mu) << " " << real(B2_mu) << " alpha delta" << endl;
  
  fichier_constantes.close();
}
	

int main(int argc, char **argv) 
{
  InitMontjoie(argc,argv);
  
  if (argc > 1)
    {
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      // instantiation avec le bon element fini
      EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension2> > Vars;
      RunAll(Vars, file_name_data, type_element);
      
      cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  cout<<"End of the program"<<endl; 
  
  return 0;
}
