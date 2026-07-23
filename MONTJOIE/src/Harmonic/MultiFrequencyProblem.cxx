#ifndef MONTJOIE_FILE_MULTI_FREQUENCY_PROBLEM_CXX

namespace Montjoie
{
  
  //! default constructor
  MultiFrequencyProblem::MultiFrequencyProblem()
  {
    wavelength_adim = Real_wp(1);
    nb_glob_parameters = 1;
    nb_parallel_proc_freq = 1;

    nb_parallel_proc_freq = 1;
    same_pattern_LU = true;
    use_iterative_solver = false;
    do_not_store_facto = false;
    nb_iter_max_eff_iterative_solver = 20;
    
    distribute_param = true;
    
    offset_omega = 0;
    offset_param = 0; nb_local_parameters = 1;
    
#ifdef SELDON_WITH_MPI
    comm_freq = MPI_COMM_WORLD;
#endif
    
  }

  
  //! reads the list of frequency from the input file
  void MultiFrequencyProblem::SetInputData(const string& keyword, const Vector<string>& param)
  { 
    if (keyword == "RangeFrequency")
      {
	Real_wp omega_min = 2*pi_wp*to_num<Real_wp>(param(0));
	Real_wp omega_max = 2*pi_wp*to_num<Real_wp>(param(1));
	int nb_omega = to_num<int>(param(2));
	Linspace(omega_min, omega_max, nb_omega, this->omega);
      }
    else if (keyword == "RangeWaveLength")
      {
        VectReal_wp L;
        Real_wp Lmin = to_num<Real_wp>(param(0));
        Real_wp Lmax = to_num<Real_wp>(param(1));
        int nb_lambda = to_num<int>(param(2));
        Linspace(Lmin, Lmax, nb_lambda, L);
        Mlt(1.0/this->wavelength_adim, L);
        this->omega.Reallocate(L.GetM());
        for (int i = 0; i < L.GetM(); i++)
          this->omega(i) = 2.0*pi_wp / L(i);

        //DISP(this->omega(0)); DISP(this->omega(L.GetM()-1));
      }
    else if (keyword == "WavelengthAdim")
      {
	wavelength_adim = to_num<Real_wp>(param(0));
      }    
    else if (keyword == "NbProcParallelFreq")
      {
	nb_parallel_proc_freq = to_num<int>(param(0));
      }
    else if (keyword == "UseIterativeSolverFreq")
      {
	if (param(0) == "YES")
          {
            use_iterative_solver = true;
            if (param.GetM() > 1)
              nb_iter_max_eff_iterative_solver = to_num<int>(param(1));
          }
        else
          use_iterative_solver = false;
      }
    else if (keyword == "UseSamePatternFreq")
      {
	if (param(0) == "YES")
	  {
	    same_pattern_LU = true;
	    do_not_store_facto = true;
	  }
	else
	  same_pattern_LU = false;	
      } 
    else if (keyword == "OutputImpedance")
      {
        int N = param.GetM();
        ref_impedance.Reallocate(N-1);
        for (int i = 0; i < ref_impedance.GetM(); i++)
          ref_impedance(i) = to_num<int>(param(i));

        file_output_impedance = param(N-1);
      }

  }     


  //! distributes the frequencies over the processors
  void MultiFrequencyProblem::DistributeFrequency(VectReal_wp& coef_input_frequency)
  {
    // then the computation of different frequencies is dispatched between processors    
#ifdef SELDON_WITH_MPI
    int nb_proc = 1, rank_proc = 0;
    MPI_Comm_size(comm_freq, &nb_proc);
    MPI_Comm_rank(comm_freq, &rank_proc);
    nb_parallel_proc_freq = min(nb_parallel_proc_freq, nb_proc);
    DISP(rank_proc); DISP(nb_proc); 

    // splitting processors into sub-groups
    int nb_effec_proc, first_proc(0), color_proc(0), key_proc;
    for (int rank = 0; rank < nb_parallel_proc_freq; rank++)
      {
	GetParallelDistributionPoints(nb_parallel_proc_freq, rank, nb_proc,
				      nb_effec_proc, first_proc, 1);
        
	if ((rank_proc >= first_proc) && (rank_proc < first_proc + nb_effec_proc))
	  {
	    color_proc = rank;
	    key_proc = rank_proc - first_proc;
	    break;
	  }
      }

    DISP(nb_effec_proc); DISP(color_proc); DISP(key_proc);
    MPI_Comm_split(comm_freq, color_proc, key_proc, &comm_group_mode);

    MPI_Comm_size(comm_group_mode, &nb_effec_proc);
    MPI_Comm_rank(comm_group_mode, &key_proc);

    // then frequencies and other parameters are dispatched in the sub-groups
    VectReal_wp all_omega(omega);
    int nb_omega_proc;

    if (nb_glob_parameters <= 1)
      distribute_param = false;
    
    int nb_parallel_omega(nb_parallel_proc_freq), nb_parallel_param(1);
    if (distribute_param)
      {
	//FindTwoFactors(nb_parallel_proc_freq, nb_parallel_param, nb_parallel_omega);
	nb_parallel_param = nb_proc / nb_parallel_omega;
      }
    
    GetParallelDistributionPoints(nb_parallel_omega, color_proc%nb_parallel_omega,
                                  all_omega.GetM(),
                                  nb_omega_proc, offset_omega, 1);

    if (distribute_param)
      {
        GetParallelDistributionPoints(nb_parallel_param, color_proc/nb_parallel_omega,
                                      nb_glob_parameters,
                                      nb_local_parameters, offset_param, 1);
      }
    else
      {
        offset_param = 0;
        nb_local_parameters = nb_glob_parameters;
      }

    omega.Reallocate(nb_omega_proc);
    for (int i = 0; i < nb_omega_proc; i++)
      omega(i) = all_omega(offset_omega + i);

    // we create also a communicator for each interval of frequencies
    int color_param = color_proc%nb_parallel_omega;
    int key_param = color_proc/nb_parallel_omega;
    MPI_Comm_split(comm_freq, color_param, key_param, &comm_param);

    color_param = 0;
    if (offset_param > 0)
      color_param = MPI_UNDEFINED;
    
    key_param = color_proc%nb_parallel_omega;
    MPI_Comm_split(comm_freq, color_param, key_param, &comm_first_param);
	
    // coef_input_frequency is also modified
    if (coef_input_frequency.GetM() > 0)
      {
        for (int i = 0; i < nb_omega_proc; i++)
          coef_input_frequency(i) = coef_input_frequency(offset_omega + i);
        
        coef_input_frequency.Resize(nb_omega_proc);
      }
#else
    nb_local_parameters = nb_glob_parameters;
#endif
  }

  
  //! computes all the frequencies of the problem
  void MultiFrequencyProblem::SolveAllFreq(VarComputationProblem& var_exp,
                                           All_LinearSolver* solver_exp, MultiFrequencyParameter& param)
  {
    int rank_proc(0); 
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(comm_freq, &rank_proc);
#endif
    DISP(rank_proc);
    
    VarBoundaryCondition_Base& var_boundary = var_exp.GetBoundaryConditionProblem();
    VarProblem_Base& var_problem = var_exp.GetVarProblemBase();
    DistributedProblem_Base& var_comm = var_exp.GetDistributedProblem();
    
    All_Preconditioner_Base<Complex_wp>* prec = var_exp.GetNewPreconditioning(Complex_wp());
    var_exp.SetPrintLevel(-1);
    prec->SetPreconditioningType(prec->DIRECT);
    
    All_IterativeSolver<Complex_wp> iterative_solver;
    if (var_exp.GetPrintLevel() >= 4)
      iterative_solver.SetPrintLevel(6);
    else if (var_exp.GetPrintLevel() >= 2)
      iterative_solver.SetPrintLevel(2);
    else
      iterative_solver.SetPrintLevel(0);
    
    if (var_comm.GetRankProcMode() > 0)
      iterative_solver.SetPrintLevel(0);
    
    iterative_solver.SetMaxNumberIteration(1000);
    iterative_solver.SetRestart(20);
    iterative_solver.SetStoppingCriterion(1e-8);  
    iterative_solver.SetSolverType(iterative_solver.COCG);
    
    iterative_solver.SetPreconditioning(*prec);

    // Finite element matrix is computed and factorized
    GlobalGenericMatrix<Complex_wp> nat_mat;
    FemMatrixFreeClass_Base<Complex_wp>* mat_iterative
      = var_exp.GetNewIterativeMatrix(Complex_wp(0));
    
    if (this->use_iterative_solver)
      {
        var_exp.SetStorageFiniteElementMatrix(var_exp.MATRIX_FREE);
        var_exp.AddMatrixWithBC(*mat_iterative, nat_mat);
      }

    Matrix<Complex_wp, General, ColMajor> sol_exp;    
    bool refacto_matrix = true;
    VectComplex_wp impedance(this->omega.GetM()); VectReal_wp freq(this->omega.GetM());
    for (int k = 0; k < this->omega.GetM(); k++)
      {
        int global_k = k;
#ifdef SELDON_WITH_MPI
        global_k = offset_omega + k;
#endif  
        
        param.SetPulsation(k, global_k, this->omega(k), var_exp, *this);
        var_problem.SetOmega(this->omega(k));

        if (var_comm.GetRankProcMode() == 0)
          cout << rank_proc << " Treating frequency " << global_k << " = " << this->omega(k) << endl;

        for (int m = 0; m < this->nb_local_parameters; m++)
          {
            int global_m = m;
#ifdef SELDON_WITH_MPI
            global_m = offset_param + m;
#endif
            param.SetParameters(k, global_k, m, global_m, var_exp, *this);
            
            // Computing geometry quantities
            var_exp.ComputeMassMatrix(false);
            if (!use_iterative_solver)
              {
                if ((k == 1) && (this->same_pattern_LU))
                  solver_exp->UseOldPattern(true);
                solver_exp->PerformFactorizationStep(nat_mat);
              }
            else
              {
                if ((k == 0) || (refacto_matrix))
                  prec->ConstructPreconditioner(var_exp, *solver_exp, *mat_iterative, nat_mat,
                                                solver_exp->GetCondensedSolver(Complex_wp(0,0)));
              }
            
            // the source for this frequency is retrieved
            param.SetSourceFrequency(k, global_k, m, global_m, var_exp, sol_exp);
            
            // then the solution is computed
            if (!use_iterative_solver)
              solver_exp->ComputeSolution(sol_exp, nat_mat);	
            else
              {
                VectComplex_wp source_exp(sol_exp.GetM()), sol_tmp(sol_exp.GetM());
                
                for (int ks = 0; ks < sol_exp.GetN(); ks++)
                  {
                    GetCol(sol_exp, ks, source_exp);
                    sol_tmp.Zero();
                    
                    var_boundary.ApplyDirichletCondition(SeldonNoTrans, *mat_iterative, source_exp);
                    
                    if ((k == 0) || refacto_matrix)
                      {
                        prec->Solve(*mat_iterative, source_exp, sol_tmp);
                        refacto_matrix = false;
                      }
                    else
                      {
                        sol_tmp.Zero();
#ifdef SELDON_WITH_MPI
                        DistributedVector<Complex_wp>* source
                          = var_comm.AllocateDistributedVector(source_exp);
                        
                        DistributedVector<Complex_wp>* solution
                          = var_comm.AllocateDistributedVector(sol_tmp);
                        
                        iterative_solver.Solve(*mat_iterative, *solution, *source);
                        
                        var_comm.NullifyDistributedVector(source);
                        var_comm.NullifyDistributedVector(solution);
                        
#else
                        iterative_solver.Solve(*mat_iterative, sol_tmp, source_exp);
#endif
                        
                        if (iterative_solver.GetNumberIteration() >= nb_iter_max_eff_iterative_solver)
                          refacto_matrix = true;	      
                      }
                    
                    mat_iterative->ImposeDirichletCondition(SeldonNoTrans, sol_tmp);
                    SetCol(sol_tmp, ks, sol_exp);
                  }
              }
            
            param.GiveSolutionFrequency(k, global_k, m, global_m, var_exp, sol_exp);
            
            //int test_input; cout << "waiting" << endl; cin >> test_input;
          }

        if (ref_impedance.GetM() > 0)
          {
            VectComplex_wp U0;
            GetCol(sol_exp, 0, U0);
            
            freq(k) = this->omega(k) / (2.0*pi_wp);
            impedance(k) = var_boundary.ComputeImpedanceCoefficient(ref_impedance, U0);
          }
      }
    
    if (ref_impedance.GetM() > 0)
      {
        GatherVectorFrequency(freq);
        GatherVectorFrequency(impedance);
        DISP(rank_proc);
        if (rank_proc == 0)
          {
            DISP(freq); DISP(impedance);
            ofstream file_out(file_output_impedance.data());
            file_out.precision(15);
            
            for (int k = 0; k < freq.GetM(); k++)
              file_out << freq(k) << " " << realpart(impedance(k)) << " " << imagpart(impedance(k)) << '\n';
            
            file_out.close();
          }
      }

    // on supprime le preconditioneur et la matrice iterative
    delete prec; delete mat_iterative;
  }

  template<class T>
  void MultiFrequencyProblem::GatherVectorFrequency(Vector<T>& x)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc = 1, rank_proc = 0;
    MPI_Comm_size(comm_freq, &nb_proc);
    MPI_Comm_rank(comm_freq, &rank_proc);
    if (nb_proc <= 1)
      return;
    
    int n0 = x.GetM(), nmax;
    MPI_Allreduce(&n0, &nmax, 1, MPI_INTEGER, MPI_MAX, comm_freq);
    
    Vector<int64_t> xtmp; 
    Vector<int> all_n(nb_proc);
    MpiGather(comm_freq, &n0, xtmp, all_n.GetData(), 1, 0);
    
    x.Resize(nmax);
    Vector<T> all_x(nmax*nb_proc);
    MpiGather(comm_freq, x, xtmp, all_x, nmax, 0);

    if (rank_proc == 0)
      {
        int nglob = 0;
        for (int i = 0; i < all_n.GetM(); i++)
          nglob += all_n(i);
        
        x.Resize(nglob); nglob = 0;
        for (int i = 0; i < nb_proc; i++)
          {
            for (int j = 0; j < all_n(i); j++)
              x(nglob+j) = all_x(nmax*i + j);
            
            nglob += all_n(i);
          }
      }
#endif
  }


  void MultiFrequencyParameter::SetPulsation(int m, int mg, const Real_wp& omega,
                                             VarComputationProblem&, MultiFrequencyProblem&)
  {
    // default case : nothing is done
  }
  

  void MultiFrequencyParameter::SetParameters(int k, int kg, int m, int mg,
                                              VarComputationProblem&, MultiFrequencyProblem&)
  {
    // default case : nothing is done
  }
  
}

#define MONTJOIE_FILE_MULTI_FREQUENCY_PROBLEM_CXX
#endif

