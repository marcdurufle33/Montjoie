#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#define SELDON_WITH_PRECONDITIONING

#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

using namespace Montjoie;

// implementation du schema saute-mouton sur l'equation d'ordre 2
// D_h (U^{n+1} - 2 U^n + U^{n-1} ) / dt^2  + A_h U^n = F^n 
// var : definition du l'operateur Ah (laplacien)
// var_time : donnees du probleme temporel
template<class TypeEquationTime, class TypeEquation>
void RunLeapFrog(EllipticProblem<TypeEquation>& var,
		 HyperbolicProblem<TypeEquationTime>& var_time)
{
  // alias pour le maillage
  MeshNumbering<Dimension3>& mesh_num = var.mesh_num;
  int nb_iter = 0; Real_wp t = var_time.GetInitialTime();
  // Nvol : nombre ddls du maillage
  int Nvol = mesh_num.GetNbDof();
  // U0 -> U^{n-1},  U1 -> U^n, Ah_u1 -> A_h U^n
  VectReal_wp U0, U1, Ah_u1;
  U0.Reallocate(Nvol);
  U1.Reallocate(Nvol);
  Ah_u1.Reallocate(Nvol);
  U0.Fill(0); U1.Fill(0); Ah_u1.Fill(0);

  Real_wp dt = var_time.GetTimeStep();
  Real_wp unp1, dt2 = dt*dt;
  while (t <= (var_time.GetFinalTime() + var_time.epsilon_time))
    {
      // on ecrit les instantanes si demande
      var_time.WriteSnapshot(nb_iter, t, U1);
		
      Ah_u1.Zero();
      // on calcule le terme source F^n
      if (t < var_time.GetFinalTimeSource())
	var_time.AddPrimitiveSourceAtTime(1.0, t, 1, Ah_u1);
      
      // on calcule Ah_u1 = F^n - A_h U^n
      var_time.ApplyOperatorKh(Real_wp(1), t, U1, Real_wp(1), Ah_u1);
            
      // calcul de D_h^{-1} ( F^n - A_h U^n)
      var_time.SolveOperatorDh(Ah_u1);
     
      // Dh (Un+1 - 2 Un + Un-1) / dt2 + Ah Un =  F^n
      for (int i = 0; i < Nvol; i++)
	{
	  unp1 = 2.0*U1(i) - U0(i) + dt2*Ah_u1(i);
	  U0(i) = U1(i); U1(i) = unp1;
	}
      
      nb_iter++;
      t = var_time.GetInitialTime() + dt*nb_iter;
    }
  
}


// implementation du theta schema sur l'equation d'ordre 2
// D_h (U^{n+1} - 2 U^n + U^{n-1} ) / dt^2  + A_h (teta U^{n+1} + (1-2 teta) U^n + teta U^{n-1}) = F^n 
// teta >= 0.25 garantit la stabilite, et teta = 0.25 donne la meilleure precision
// var : definition du l'operateur Ah (laplacien)
// var_time : donnees du probleme temporel
// glob_solver : solveur utilise pour resoudre le systeme lineaire
template<class TypeEquationTime, class TypeEquation>
void RunThetaScheme(EllipticProblem<TypeEquation>& var,
		    HyperbolicProblem<TypeEquationTime>& var_time,
		    All_LinearSolver& glob_solver, GlobalGenericMatrix<Real_wp>& nat_mat)
{
  // on recupere le maillage
  MeshNumbering<Dimension3>& mesh_num = var.mesh_num;
  int nb_iter = 0; Real_wp t = var_time.GetInitialTime();
  // Nvol : nombre ddl maillages
  int Nvol = mesh_num.GetNbDof();
  // U0 -> U^{n-1}, U1 -> U^n, Ah_u0 -> dt^2 A_h U^n-1, Ah_u1 -> dt^2 A_h U^n
  VectReal_wp U0, U1, Ah_u1, Ah_u0, Prod_Uh;
  U0.Reallocate(Nvol);
  U1.Reallocate(Nvol); Prod_Uh.Reallocate(Nvol);
  Ah_u1.Reallocate(Nvol);   Ah_u0.Reallocate(Nvol);
  U0.Fill(0); U1.Fill(0);
  Ah_u1.Fill(0); Ah_u0.Fill(0); Prod_Uh.Fill(0);

  Real_wp dt = var_time.GetTimeStep();
  Real_wp dt2 = dt*dt;
  Real_wp teta = 0.25;
  DISP(var_time.GetInitialTime()); DISP(var_time.GetFinalTime());
  // a modifier si on veut verifier que le solveur converge
  // 0 -> aucune information n'est affichee
  // 3 -> toutes les iterations sont affichees
  glob_solver.SetPrintLevel(0);
  while (t <= (var_time.GetFinalTime() + var_time.epsilon_time))
    {
      // on ecrit les instantanes sur le disque
      var_time.WriteSnapshot(nb_iter, t, U1);
		
      Copy(Ah_u1, Ah_u0); Prod_Uh.Zero();
      // on calcule le terme source dt^2 F^n
      if (t < var_time.GetFinalTimeSource())
	{
	  var_time.AddPrimitiveSourceAtTime(1.0, t, 1, Prod_Uh);
	  Mlt(dt2, Prod_Uh);
	}
      
      // on calcule Ah_u1 = - (Delta t)^2 A_h U^n
      var_time.ApplyOperatorKh(Real_wp(dt2), t, U1, Real_wp(0), Ah_u1);
      // on ajoute D_h ( 2 U^n - U^{n-1})  a Prod_Uh
      Add(-2.0, U1, U0);
      var_time.ApplyOperatorDh(Real_wp(-1), t, U0, Real_wp(1), Prod_Uh);
      
      // maintenant on resout (D_h + teta dt^2 Ah) U^{n+1}
      //          = Dh(2 U^n - U^{n-1}) + F^n  + (1- 2 teta) dt^2 Ah U^n + teta Ah U^{n-1}
      Add(1.0-2.0*teta, Ah_u1, Prod_Uh);
      Add(teta, Ah_u0, Prod_Uh);
      
      Copy(U1, U0);
      var_time.SetDirichletCondition(t+dt, Prod_Uh);
      U1 = Prod_Uh;
      glob_solver.ComputeSolution(U1, nat_mat);
      
      nb_iter++;
      t = var_time.GetInitialTime() + dt*nb_iter;
    }
  
}


// on lance la simulation generale
// var_time : donnees du probleme a simuler
// input_file : fichier de donnees
template<class TypeEquationTime>
void RunAll(HyperbolicProblem<TypeEquationTime>& var_time, const string& input_file,
            const string& name_element, const string& name_equation)
{
  // on recupere le probleme stationnaire
  typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
  EllipticProblem<TypeEquation>& var = var_time.var_harmonic;
  
  // on lit le fichier de donnees
  var.InitIndices(PhysicalConstant::nb_max_indices);
  ReadInputFile(input_file, var_time);
  
  // solveur utilise pour le theta schema
  All_LinearSolver& glob_solver = *var_time.glob_solver;
  
  ReadInputFile(input_file, glob_solver);
  
  // construction du maillage et element fini
  var.ComputeMeshAndFiniteElement(name_element);

  // calcul du jacobien, DF_i^{-1} 
  var.ComputeMassMatrix();
  // on enleve DF_i^{-1}
  var.Glob_DFjm1.Clear();
  
  // on repere noeuds de Dirichlet
  var.TreatDirichletCondition();
  var_time.ComputeRightHandSide();
  
  // autres initialisations pour le temporel
  var_time.ComputeMassMatrix();
  var_time.InitTimeIterations();
  
  // calcul de la matrice de rigidite
  GlobalGenericMatrix<Real_wp> nat_mat;
  var_time.Glob_mat_Kh->Clear();
  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(1.0);
  var_time.Glob_mat_Kh->SetCoefficientMatrix(nat_mat);
  var.AddMatrixFEM(*var_time.Glob_mat_Kh, nat_mat);

  // type de schema ?
  if (var_time.GetTimeSchemeType() == TimeSchemeEnum::THETA_SCHEME)
    {
      // schema implicite
      // on doit resoudre (D_h + teta dt^2 K_h)
      Real_wp dt2 = square(var_time.GetTimeStep());
      nat_mat.SetCoefMass(1.0);
      nat_mat.SetCoefStiffness(dt2*0.25);
      glob_solver.SetPrintLevel(5);
      glob_solver.PerformFactorizationStep(nat_mat);
      
      // attention, on efface ici les faces et aretes du maillages
      // a commenter si necessaire
      var.mesh.ClearConnectivity();
      
      // teta-schema classique
      RunThetaScheme(var, var_time, glob_solver, nat_mat);
    }
  else
    {
      // attention, on efface ici les faces et aretes du maillages
      var.mesh.ClearConnectivity();
      
      // schema saute-mouton classique
      RunLeapFrog(var, var_time);
    }
}

// fonction principale
int main(int argc, char **argv) 
{
  // initialisation de Montjoie, ne pas enlever
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // TetrahedronClassical -> pas de condensation de masse
      // TetrahedronQuasiLumped -> quasi-condensation de masse
      // TetrahedronMassLumped -> condensation de masse      
      HyperbolicProblem<AcousticEquation<Dimension3> > var;

      string file_name_data(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      RunAll(var, file_name_data, type_element, type_equation);
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"lung.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
