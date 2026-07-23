#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Elliptic/Aeroacoustic/MontjoieAeroacousticAxisym.hxx"

using namespace Montjoie;

// probleme harmonique
template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& name_file,
	    const string& name_element, const string& name_equation)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;

  MontjoieTimer var_chrono;

#ifdef SELDON_WITH_MPI
  var.var_chrono.SetCommunicator(var.comm_group_mode);
  var_chrono.SetCommunicator(var.comm_group_mode);
  int rank_proc; MPI_Comm_rank(var.comm_group_mode, &rank_proc);
#else
  int rank_proc(0);
#endif
  
  var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var_chrono.Start("GlobalComputation");
  
  var_chrono.SetMessage("ConstructAll", "initialize the computation");
  var_chrono.Start("ConstructAll");

  Vector<string> lines_data_file;
  
  var.ConstructAll(name_file, name_element, lines_data_file);
  
  All_LinearSolver glob_solver(var);
  ReadInputFile(lines_data_file, glob_solver);
  
  // right hand side, solution and vector used to display solution
  VectComplex_wp source_rhs, x_sol, output_vector;
  
  int Nvol = var.GetNbDof();
  x_sol.Reallocate(Nvol);
  x_sol.Fill(0);
    
  GlobalGenericMatrix<Complex_wp> nature_matrix;
  
  // loop on the modes
  var.last_output_solution = false;
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();
  //var.WriteAllIndices();

  var_chrono.Stop("ConstructAll");

  var_chrono.SetMessage("ComputationRhs", "compute the right hand side");
  var_chrono.SetMessage("FactorisationMatrix", "compute and factorize finite element matrix");
  var_chrono.SetMessage("ComputationSolution", "solves linear systems");
  var_chrono.SetMessage("WriteSolution", "interpolate the solution on grid points");
 
  //int number_mode_prev = 1000000;
  int n = 0; Real_wp norme_max(0);
  while (!var.last_output_solution)
    {
      var_chrono.Start("FactorisationMatrix");
      
      // mode number, that we consider
      var.SetCurrentModeNumber(var.GetModeNumber(n));
      if (rank_proc == 0)
        cout << "We solve mode number " << var.GetCurrentModeNumber() << endl;
      
      // computation and factorisation of finite element matrix
      //if (abs(TypeEquation::number_mode) != abs(number_mode_prev))
      glob_solver.PerformFactorizationStep(nature_matrix);      
      
      var_chrono.Stop("FactorisationMatrix");
      var_chrono.Start("ComputationRhs");
      
      // computation of right hand side
      var.ComputeRightHandSide(source_rhs);

      var_chrono.Stop("ComputationRhs");
      var_chrono.Start("ComputationSolution");
      
      // we solve the linear system LU x = b
      x_sol = source_rhs;
      glob_solver.ComputeSolution(x_sol, nature_matrix);
      
      var_chrono.Stop("ComputationSolution");
      
      // we retrieve the infinite norm of x
      Real_wp x_inf = x_sol.GetNormInf();
#ifdef SELDON_WITH_MPI
      Real_wp max_x = x_inf; Vector<int64_t> xtmp;
      MpiAllreduce(var.comm_group_mode, &max_x, xtmp, &x_inf, 1, MPI_MAX);    
#endif
      
      if (rank_proc == 0)
        cout << "Infinite norm of x = " << x_inf << endl;
      
      norme_max = max(norme_max, x_inf);

      // then we update the number of modes if needed
      if ((var.NumberOfModesToBeComputed()) && (n == var.GetNbModes()-1))
        {
          if ((norme_max > 0) && (x_inf/norme_max > var.GetModeThreshold()))
            var.PushBackMode(var.GetCurrentModeNumber() + 1);
        }
      
      if (n == var.GetNbModes()-1)
	var.last_output_solution = true;
      
      if (rank_proc == 0)
        cout<<"Solution computed successfully"<<endl;

      var_chrono.Start("WriteSolution");
  
      var.WriteDatas(x_sol);
  
      var_chrono.Stop("WriteSolution");
      
      //number_mode_prev = TypeEquation::number_mode;
      n++;
      //int test_input; cout << "Waiting..." << endl; cin >> test_input;
    } // end loop over modes

  var_chrono.Stop("GlobalComputation");
  
  var_chrono.DisplayAll();
  var.var_chrono.DisplayTime("InterpolationGrid");  
} // end method RunAll


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
  int rank_proc(0);
#endif

  if (argc > 1)
    {
      
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
            
      if (rank_proc == 0)
        cout<<"Helmholtz Solver with Axi-symmetric geometry"<<endl; 
      
      // probleme harmonique
      EllipticProblem<HarmonicGalbrunEquationAxi> Vars;
      RunAll(Vars, file_name_data, type_element, type_equation);
      
      if (rank_proc == 0)
        cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz_axi.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  if (rank_proc == 0)
    cout<<"End of the program"<<endl; 

  return FinalizeMontjoie();
}
