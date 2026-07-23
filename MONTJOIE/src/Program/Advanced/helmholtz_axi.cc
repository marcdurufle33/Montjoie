#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"
#endif

//#include "Hyperbolic/Acoustic/AxiSymAcoustic.hxx"
//#include "Hyperbolic/Acoustic/AxiSymAcoustic.cxx"

using namespace Montjoie;

namespace Montjoie
{
  void GetTimeWave(const Real_wp& omega, const Real_wp& t,
                   const Real_wp& tau_impulse, const Real_wp& Tmax, Real_wp& pulse)
  {
    Real_wp phase_init = 0;
    pulse = exp(-0.5*square((t-Tmax)/tau_impulse));
    pulse *= cos(omega*t + phase_init);
  }
}  

// probleme temporel
/* template<class TypeElement, class TypeEquationTime, class Complexe>
void RunAll(HyperbolicProblem<TypeElement, TypeEquationTime>& var_time, const string& name_file, const Complexe& one)
{

typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;
  
  typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
  EllipticProblem<TypeElement, TypeEquation>& var = var_time.var_harmonic;
  
  // initialisation indices physiques
  var.InitIndices(PhysicalConstant::nb_max_indices);
  
  // lecture du fichier de donnees
  Vector<string> lines_data_file;
  
#ifdef SELDON_WITH_MPI
  ReadLinesFile(name_file, lines_data_file, MPI_COMM_WORLD);
#else
  ReadLinesFile(name_file, lines_data_file);
#endif

  ReadInputFile(lines_data_file, var_time);
  
  // construction du numero de dossier
  // si le dossier n'est pas connu
  //if (root_processor == rank_processor)
    //if (var.DOSSIER.size() == 0)
      //{
	//int num = -1;
	//var.DOSSIER = "[STIFFOUT]/"; 
	//EcritDossier(var.DOSSIER, name_file, num);
        //}
    
#ifdef SELDON_WITH_MPI
  //if (nb_processors > 1)
  //MPI_Bcast_string(var.DOSSIER, root_processor, MPI_COMM_WORLD);
#endif    

  All_LinearSolver<TypeElement, TypeEquation> glob_solver(var);
  ReadInputFile(lines_data_file, glob_solver);
  
  // construction maillage et elements finis
  var.mesh_num.SetSameNumberPeriodicDofs();
  var.ComputeMeshAndFiniteElement();
  
  // quantites geometriques et autres ddls
  var.ComputeMassMatrix();
  var.ComputeQuasiPeriodicPhase();

  var.PerformOtherInitializations();
  
  // on ne traite que le mode 0
  TypeEquation::number_mode = 0;
  var.UpdateDirichlet(0);
  
  // calcul de la source
  var_time.ComputeRightHandSide();

  IVect num;
  var.SetDirichletDofs(0, num);
  
  // calcul des matrices
  var.threshold_matrix = 0;
  FemMatrixFreeClass<TypeElement, TypeEquation> Kh(var), Mh(var), Sh(var);
  typename TypeEquation::Nature_Matrix nat_mat;
  Kh.IgnoreDirichletDof();   Mh.IgnoreDirichletDof();   Sh.IgnoreDirichletDof();
  Kh.SetCoefficientDirichlet(0.0);
  Mh.SetCoefficientDirichlet(0.0);
  Sh.SetCoefficientDirichlet(0.0);
 
  nat_mat.SetCoefMass(1.0);
  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(0.0);
  var.AddMatrixWithBC(Mh, nat_mat);
  
  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefDamping(1.0);
  var.AddMatrixWithBC(Sh, nat_mat);

  nat_mat.SetCoefDamping(0.0);
  nat_mat.SetCoefStiffness(1.0);
  var.AddMatrixWithBC(Kh, nat_mat);  

  //Real_wp fwhm_impulse = 1.0;
  Real_wp fwhm_impulse = 1.0/(pi_wp*7e8*1.5e-3);
  Real_wp tau_impulse = fwhm_impulse /(Real_wp(2) * sqrt(Real_wp(2)*log(Real_wp(2)))); 
  
  //Real_wp tau_impulse = 1.0/(2.0*pi_wp*7e8*1.0e-3);
  Real_wp Tmax = 8.0*tau_impulse;
  DISP(tau_impulse); DISP(Tmax);
  
  // initialisation des iterations en temps
  var_time.InitTimeIterations();
  
  // iterations en temps
  Real_wp dt = var_time.deltat;
  int nb_iterations_time = toInteger(ceil(var_time.final_time/dt));
  Real_wp theta = 0.25;

  int N = var.GetNbDof();
  VectReal_wp Mh_vec, Sh_vec;
  if (theta == Real_wp(0))
    {
      var.ComputeDiagonalMatrix(Mh_vec, Mh, nat_mat, false);
      var.ComputeDiagonalMatrix(Sh_vec, Sh, nat_mat, false);
      for (int i = 0; i < var_time.dirichlet_source.GetM(); i++)
        {
          Mh_vec(var_time.dirichlet_source.Index(i)) = 1.0;
          Sh_vec(var_time.dirichlet_source.Index(i)) = 0.0;
        }
      
      for (int i = 0; i < N; i++)
        Mh_vec(i) = Mh_vec(i) + 0.5*dt*Sh_vec(i);
      
      var_time.Assemble(Mh_vec); 
      
      for (int i = 0; i < N; i++)
        Mh_vec(i) = 1.0/Mh_vec(i);
    }
  else
    {
      nat_mat.SetCoefMass(1.0);
      nat_mat.SetCoefDamping(dt/2);
      nat_mat.SetCoefStiffness(theta*dt*dt);
      
      glob_solver.PerformFactorizationStep(nat_mat);
      glob_solver.SetPrintLevel(0);
    }
  
  Vector<Complexe> Un(N), Un_prev(N), ShUn(N), ShUn_prev(N), KhUn(N), KhUn_prev(N);
  Vector<Complexe> Un_next(N), rhs(N);
  Un.Fill(0); Un_prev.Fill(0);  ShUn.Fill(0); ShUn_prev.Fill(0);
  KhUn.Fill(0); KhUn_prev.Fill(0); Un_next.Fill(0);

  Real_wp eps = 1e-6*dt;
  for (int i = 0; i < var.output_grid_param.GetM(); i++)
    var.output_grid_param(i).SetThresholdTime(eps);

  //DISP(var_time.sparse_vector_source);
  //Mh_vec.WriteText("Mh.dat");
  // boucle en temps
  Complexe val, pulse;
  for (int nt = 0; nt < nb_iterations_time; nt++)
    {
      Real_wp t = nt*dt;
      if (root_processor == rank_processor)
        if (nt%100 == 0)
          cout << "t = " << t << ", || Un || = " << Norm2(Un) << endl;
      
      //for (int i = 0; i < N; i++)
      //to_complex(Un(i), Un_real(i));
      
      //var.WriteSnapshot(nt, t, Un_real);
      var.WriteSnapshot(nt, t, Un);
      rhs.Fill(0);

      if (t < (2.0*Tmax))
        {
          GetTimeWave(var.GetOmega(), t+dt, tau_impulse, Tmax, pulse);
          //DISP(pulse);
          for (int i = 0; i < var_time.sparse_vector_source.GetM(); i++)
            rhs(var_time.sparse_vector_source.Index(i)) = pulse*var_time.sparse_vector_source.Value(i);
        }
      
      // dt^2 (F^n - Kh U^n) - dt Sh (U^n - U^n-1)
      for (int i = 0; i < N; i++)
        {
          val = dt*dt*(-KhUn(i) + rhs(i)) - dt*(ShUn(i) - ShUn_prev(i));          
          rhs(i) = val;
        }
      
      var_time.Assemble(rhs);

      if (theta == Real_wp(0))
        for (int i = 0; i < N; i++)
          Un_next(i) = 2.0*Un(i) - Un_prev(i) + rhs(i)*Mh_vec(i);
      else
        {
          glob_solver.ComputeSolution(rhs, Un_next, nat_mat);
          for (int i = 0; i < N; i++)
            Un_next(i) += 2.0*Un(i) - Un_prev(i);
        }
      
      //Un_next.Write("Un"+suffix);
      
      //Un.WriteText("Un.dat");
      //Un_prev.WriteText("Un_prev.dat");
      //Un_next.WriteText("Un_next.dat");
      //int test_input; cout << "we wait" << endl; cin >> test_input;
      
      // on passe a l'itere suivant
      Copy(ShUn, ShUn_prev);
      Copy(KhUn, KhUn_prev); Copy(Un, Un_prev);
      MltAdd(1.0, SeldonNoTrans, Kh, Un_next, 0.0, KhUn, false);
      MltAdd(1.0, SeldonNoTrans, Sh, Un_next, 0.0, ShUn, false);
      
      Copy(Un_next, Un);
      //int test_input; cout << "Waiting..." << endl; cin >> test_input;
    }  
    } */

// probleme harmonique
template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& name_elt, const string& name_file)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;
  
  MontjoieTimer var_chrono;

#ifdef SELDON_WITH_MPI
  int rank_proc; MPI_Comm_rank(var.comm_group_mode, &rank_proc);
  var.var_chrono.SetCommunicator(var.comm_group_mode);
  var_chrono.SetCommunicator(var.comm_group_mode);
#else
  int rank_proc(0);
#endif
  
  var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var_chrono.Start("GlobalComputation");
  
  var_chrono.SetMessage("ConstructAll", "initialize the computation");
  var_chrono.Start("ConstructAll");

  Vector<string> lines_data_file;
  
  var.ConstructAll(name_file, name_elt, lines_data_file);
  
  All_LinearSolver glob_solver(var);
  ReadInputFile(lines_data_file, glob_solver);
  
  // right hand side, solution and vector used to display solution
  VectComplex_wp source_rhs, x_sol, output_vector;
  
  int Nvol = var.GetNbDof();
  x_sol.Reallocate(Nvol);
  x_sol.Fill(0);
    
  GlobalGenericMatrix<Complex_wp> nature_matrix;
  
  All_Preconditioner<Complex_wp, TypeEquation> prec(var);
  prec.SetPreconditioningType(prec.DIRECT);

  All_IterativeSolver<Complex_wp> iterative_solver;
  int nb_iter_max_eff = 30;
  iterative_solver.SetPrintLevel(6);
  iterative_solver.SetMaxNumberIteration(1000);
  iterative_solver.SetRestart(20);
  iterative_solver.SetStoppingCriterion(1e-6);
  
  if (var.include_flow_term)
    iterative_solver.SetSolverType(iterative_solver.GMRES);
  else
    iterative_solver.SetSolverType(iterative_solver.COCG);
  
  iterative_solver.SetPreconditioning(prec);
  
  FemMatrixFreeClass<Complex_wp, TypeEquation> mat_iterative(var);
	      
  // loop on the modes
  var.last_output_solution = false;
  var.ComputeMassMatrix();

  if (var.use_iterative_solver)
    var.AddMatrixWithBC(mat_iterative, nature_matrix);

  var_chrono.Stop("ConstructAll");

  var_chrono.SetMessage("ComputationRhs", "compute the right hand side");
  var_chrono.SetMessage("FactorisationMatrix", "compute and factorize finite element matrix");
  var_chrono.SetMessage("ComputationSolution", "solves linear systems");
  var_chrono.SetMessage("WriteSolution", "interpolate the solution on grid points");
 
  int number_mode_prev = 1000000;
  int n = 0; Real_wp norme_max(0);
  bool refacto_matrix = true;
  CondensationBlockSolver_Fem<Complex_wp>* other_cond
    = static_cast<CondensationBlockSolver_Fem<Complex_wp>* >(var.GetNewCondensationSolver(Complex_wp(0)));

  while (!var.last_output_solution)
    {
      var_chrono.Start("FactorisationMatrix");
      
      // mode number, that we consider
      var.SetCurrentModeNumber(var.GetModeNumber(n));
      var.UpdateDirichlet(n);
      if (rank_proc == 0)
        cout << "We solve mode number " << var.GetCurrentModeNumber() << endl;
      
      // computation and factorisation of finite element matrix
      if (!var.use_iterative_solver)
	{
	  if ((n == 1) && (var.same_pattern_LU))
	    glob_solver.UseOldPattern(true);
	  
	  if ((abs(var.GetCurrentModeNumber()) != abs(number_mode_prev)) || (var.include_flow_term))
	    glob_solver.PerformFactorizationStep(nature_matrix);    	  
	}
      else
	{
	  //var.GetDirichletDofNumber().WriteText("dir.dat");
	  //glob_solver.PerformFactorizationStep(nature_matrix);    	  
	  if ((n <= 1) || (refacto_matrix))
	    {
	      // direct solver for the two first modes
	      prec.ConstructPreconditioner(var, glob_solver, mat_iterative, nature_matrix, *other_cond);
	    }
	}
      
      var_chrono.Stop("FactorisationMatrix");
      var_chrono.Start("ComputationRhs");
      
      // computation of right hand side
      var.ComputeRightHandSide(source_rhs);

      var_chrono.Stop("ComputationRhs");
      var_chrono.Start("ComputationSolution");
      
      // we solve the linear system LU x = b
      if (!var.use_iterative_solver)
	{
	  x_sol = source_rhs;
	  glob_solver.ComputeSolution(x_sol, nature_matrix);
	}
      else
	{
	  var.ApplyDirichletCondition(SeldonNoTrans, mat_iterative, source_rhs);
	  
	  if ((n <= 1) || refacto_matrix)
	    {
	      prec.Solve(mat_iterative, source_rhs, x_sol);
	      refacto_matrix = false;
	    }
	  else
	    {
	      x_sol.Zero();
#ifdef SELDON_WITH_MPI
	      DistributedVector<Complex_wp>* source
		= var.AllocateDistributedVector(source_rhs);
            
	      DistributedVector<Complex_wp>* solution
		= var.AllocateDistributedVector(x_sol);

	      iterative_solver.Solve(mat_iterative, *solution, *source);

	      var.NullifyDistributedVector(source);
	      var.NullifyDistributedVector(solution);

#else
	      iterative_solver.Solve(mat_iterative, x_sol, source_rhs);
#endif

	      if (iterative_solver.GetNumberIteration() >= nb_iter_max_eff)
		refacto_matrix = true;	      
	    }
	  
	  mat_iterative.ImposeDirichletCondition(SeldonNoTrans, x_sol);	  
	}
      
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
            var.PushBackMode(var.GetCurrentModeNumber()+1);
        }
      
      if (n == var.GetNbModes()-1)
	var.last_output_solution = true;

      if (rank_proc == 0)
        cout<<"Solution computed successfully"<<endl;

      var_chrono.Start("WriteSolution");
  
      // si on veut recuperer la solution sur les points de la grille 0
      //int num_grid = 0; bool compute_grad = true;
      //Vector<VectComplex_wp> trace_u; Vector<VectComplex_wp> trace_grad_u;
      //var.ComputeInterpolationU_GradU(x_sol, num_grid, compute_grad, trace_u, trace_grad_u);
      //WriteMatlab(trace_u(0), var.var_grid(0), "toto.dat", true);
      
      // pour sortir la solution sur toutes les grilles et maillages
      var.WriteDatas(x_sol);
      
      var_chrono.Stop("WriteSolution");
      
      number_mode_prev = var.GetCurrentModeNumber();
      n++;
      //DISP(var.GetCurrentModeNumber());
      //DISP(var.NumberOfModesToBeComputed());
      //DISP(var.GetNbModes());
      //DISP(n);
      //int test_input; cout << "Waiting..." << endl; cin >> test_input;
    } // end loop over modes

  var_chrono.Stop("GlobalComputation");
  
  var_chrono.DisplayAll();

  cout << endl << endl;
  
  if (!var.use_iterative_solver)
    {
      var.var_chrono.DisplayTime("ComputationMatrix");      
      var.var_chrono.DisplayTime("FactorizationMatrix");
    }
  
  var.var_chrono.DisplayTime("InterpolationGrid");  
  delete other_cond;
    
} // end method RunAll


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
            
      EllipticProblem<HelmholtzEquationAxi> Vars;
      
      if (Vars.GetRankProcMode() == 0)
        cout<<"Helmholtz Solver with Axi-symmetric geometry"<<endl; 
      
      // probleme harmonique
      RunAll(Vars, type_element, file_name_data);
      
      // probleme temporel

      //HyperbolicProblem<TriangleRadau, AcousticEquationAxi> Vars;
      
      // calcul en reel ou complexe ?
      //Complex_wp one(1.0, 0.0);
      //Real_wp one(1.0);
      
      // on lance la simulation
      //RunAll(Vars, file_name_data, one);

      
      if (Vars.GetRankProcMode() == 0)
        cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz_axi.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
