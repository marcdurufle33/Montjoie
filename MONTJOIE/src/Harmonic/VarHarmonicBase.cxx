#ifndef MONTJOIE_FILE_VAR_HARMONIC_BASE_CXX

namespace Montjoie
{

  //! modifies parameters of the problem with a line of the data file
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {    
    // we call SetInputData of parent classes
    VarBoundaryCondition<Complexe, Dimension>::SetInputData(description_field, parameters);
    VarOutputProblem<Dimension>::SetInputData(description_field, parameters);
    VarSourceProblem<Complexe, Dimension>::SetInputData(description_field, parameters);
    VarProblem<Dimension>::SetInputData(description_field, parameters);

#ifdef MONTJOIE_WITH_TRANSMISSION
    var_transmission_base.SetInputData(description_field, parameters);
#endif

    var_gibc_base.SetInputData(description_field, parameters);
  }
  

  //! Other initialization before the matrix computation
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>::PerformOtherInitializations()
  {
    // For DG formulation, in order to get same solution
    // than 2nd order formulation, source has to be divided by omega
    if (this->IsComplexProblem() && (this->FormulationDG() != ElementReference_Base::CONTINUOUS)
        && this->FirstOrderFormulationDG())
      this->coefficient_volumic_source /= this->GetOmega();
    
#ifdef MONTJOIE_WITH_WIRES
    var_wire_base.InitWires();
#endif

#ifdef MONTJOIE_WITH_TRANSMISSION
    var_transmission_base.InitTransmission();
#endif

    var_gibc_base.InitGIBC();

#ifdef SELDON_WITH_MPI
    // arrays needed to assemble Dirichlet dofs are constructed    
    this->ConstructDirichletComm();
#endif
  }


  //! restarting a computation with the same object
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>::Restart()
  {
    this->mesh_num.treat_periodic_condition_during_number = true;
    for (int n = 0; n < this->other_mesh_num.GetM(); n++)
      this->other_mesh_num(n)->treat_periodic_condition_during_number = true;
  }
  

  //!  constructs a finite element problem
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::ConstructAll(const string& input_file, const string& name_elt,
		 const string& name_equation, All_LinearSolver*& glob_solver,
		 bool compute_rho, bool delete_points, int num)
  {
    // we read the input file if input_file is different from NONE
    // if it is NONE, we suppose that parameters have been modified by hand
    Vector<string> lines_data_file;
    if (input_file.compare("NONE"))
      {    
#ifdef SELDON_WITH_MPI
        ReadLinesFile(input_file, lines_data_file, this->comm_group_mode);
#else
        ReadLinesFile(input_file, lines_data_file);
#endif
      }
    
    ConstructAll(lines_data_file, input_file, name_elt, name_equation, glob_solver,
		 compute_rho, delete_points, num);
  }

  
  //! constructs a finite element problem
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::ConstructAll(const Vector<string>& lines_data_file,
		 const string& input_file, const string& name_elt,
		 const string& name_equation, All_LinearSolver*& glob_solver, 
		 bool compute_rho, bool delete_points, int num)
  { 
    this->SetTypeEquation(name_equation);
    
    // initialization of physical properties
    this->InitIndices(PhysicalConstant::nb_max_indices);
    
    ReadInputFile(lines_data_file, *this);

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    glob_solver = this->GetNewLinearSolver();
    ReadInputFile(lines_data_file, *glob_solver);

    if (rank_proc == 0)
      if (this->print_level >= 0)
	cout << rank_proc << " The input file has been read"<<endl;
    
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (this->DOSSIER_output.size() == 0)
      {
	if (rank_proc == 0)
          {
            this->DOSSIER_output = "[STIFFOUT]/"; 
            EcritDossier(this->DOSSIER_output, input_file, num);
          }
        
#ifdef SELDON_WITH_MPI
        MPI_Bcast_string(this->DOSSIER_output, 0, this->comm_group_mode);
#endif
        
      }

#ifdef SELDON_WITH_MPI
    this->var_chrono.SetCommunicator(this->comm_group_mode);
#endif

    this->var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
    this->var_chrono.Start("GlobalComputation");
    
    this->var_chrono.SetMessage("MeshGeneration", "construct and number the mesh");
    this->var_chrono.Start("MeshGeneration");
    
    bool split_mesh = true;
    if (this->GetNbProcPerMode() == 1)
      split_mesh = false;

    // we read mesh and construct reference element (for finite element method)
    this->ComputeMeshAndFiniteElement(name_elt, split_mesh);
    
    // adding dofs for other variables
    this->PerformOtherInitializations();
    
    if (this->print_level >= 0)
      if (rank_proc == 0)
        {
	  if (this->dg_formulation == ElementReference_Base::DISCONTINUOUS)
	    cout << "Discontinuous Galerkin formulation is used " << endl;
	  else if (this->dg_formulation == ElementReference_Base::HDG)
	    cout << "HDG formulation is used " << endl;
	  else
	    cout << "Continuous Galerkin formulation is used " << endl;
	  
          if (this->GetNbProcPerMode() > 1)
            {
              cout << "Number of global scalar dofs " << this->nodl_mesh << endl;
              cout << "Number of global dofs for the overall problem : " << this->nodl_all << endl;
            }
          else
            {
              cout << "Number of scalar dofs " << this->mesh_num.GetNbDof() << endl;
              cout << "Number of dofs for the overall problem : " << this->GetNbDof() << endl;
            }
        }
    
    this->var_chrono.Stop("MeshGeneration");

    // selection of solver 
    Dimension dim; Symmetric property;
    glob_solver->SelectOptimalLinearSolver(this->mesh_num.GetOrder(), this->nodl, dim, property);
    
    if (this->GetNbProcPerMode() == 1)
      {
        if (this->print_level >= 2)
          this->mesh.Write(this->DOSSIER_output+string("test.mesh"));    
      }
    else
      {
        if (this->print_level >= 6)
	  this->mesh.Write(this->DOSSIER_output+string("test")
			   +to_str(rank_proc)+".mesh");    
      }

    this->var_chrono.SetMessage("JacobianMatricesComputation", "compute jacobian matrices");
    this->var_chrono.Start("JacobianMatricesComputation");
    
    // computation of Ji, DFi on the mesh
    this->ComputeMassMatrix(compute_rho, delete_points);

    this->var_chrono.Stop("JacobianMatricesComputation");
    if (this->print_level >= 1)
      this->var_chrono.DisplayTime("JacobianMatricesComputation");

    // computation of phase for quasi-periodic condition
    this->ComputeQuasiPeriodicPhase();
  }
  

  //! runs a complete simulation
  /*!
    The simulation is decomposed in several steps :
    - construction of the mesh, and finite element
    - computation of finite element matrix
    - factorization for direct solving or construction
    of preconditioner for iterative solving
    - right hand side computation
    - resolution of linear system (direct or iterative)
    - post-processing
  */
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::RunAll(const string& input_file, const string& name_elt,
	   const string& name_equation, int num)
  {
    // we read the input file if input_file is different from NONE
    // if it is NONE, we suppose that parameters have been modified by hand
    Vector<string> lines_data_file;
    if (input_file.compare("NONE"))
      {    
#ifdef SELDON_WITH_MPI
        ReadLinesFile(input_file, lines_data_file, this->comm_group_mode);
#else
        ReadLinesFile(input_file, lines_data_file);
#endif
      }
    
    RunAll(lines_data_file, input_file, name_elt, name_equation, num);
  }

  
  //! runs a complete simulation
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::RunAll(const Vector<string>& lines_data_file,
	   const string& input_file, const string& name_elt, const string& name_equation, int num)
  { 
    this->SetTypeEquation(name_equation);

    // initialization of physical properties
    this->InitIndices(PhysicalConstant::nb_max_indices);
    
    ReadInputFile(lines_data_file, *this);

#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(this->comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    All_LinearSolver* glob_solver;
    glob_solver = this->GetNewLinearSolver();
    ReadInputFile(lines_data_file, *glob_solver);

    TransparencySolver_Base* solver_transparent;
    solver_transparent = this->GetNewTransparentSolver(*glob_solver);
    ReadInputFile(lines_data_file, *solver_transparent);
    ReadInputFile(lines_data_file, output_rcs_base);
    
    if (rank_proc == 0)
      if (this->print_level >= 0)
	cout << rank_proc << " The input file has been read" << endl;
    
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (this->DOSSIER_output.size() == 0)
      {
	if (rank_proc == 0)
          {
            this->DOSSIER_output = "[STIFFOUT]/"; 
            EcritDossier(this->DOSSIER_output, input_file, num);
          }
        
#ifdef SELDON_WITH_MPI
        MPI_Bcast_string(this->DOSSIER_output, 0, this->comm_group_mode);
#endif
        
      }
    
    // variable storing size of different objects
    map<string, size_t> size_object;
    
#ifdef SELDON_WITH_MPI
    this->var_chrono.SetCommunicator(this->comm_group_mode);
#endif

    this->var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
    this->var_chrono.Start("GlobalComputation");
    
    this->var_chrono.SetMessage("MeshGeneration", "construct and number the mesh");
    this->var_chrono.Start("MeshGeneration");
    
    bool split_mesh = true;
    if (this->GetNbProcPerMode() == 1)
      split_mesh = false;

    // we read mesh and construct reference element (for finite element method)
    this->ComputeMeshAndFiniteElement(name_elt, split_mesh);
    
    // adding dofs for other variables
    this->PerformOtherInitializations();

    if (this->print_level >= 0)
      if (rank_proc == 0)
        {
	  if (this->dg_formulation == ElementReference_Base::DISCONTINUOUS)
	    cout << "Discontinuous Galerkin formulation is used " << endl;
	  else if (this->dg_formulation == ElementReference_Base::HDG)
	    cout << "HDG formulation is used " << endl;
	  else
	    cout << "Continuous Galerkin formulation is used " << endl;
	  
          if (this->GetNbProcPerMode() > 1)
            {
              cout << "Number of global scalar dofs " << this->nodl_mesh << endl;
              cout << "Number of global dofs for the overall problem : " << this->nodl_all << endl;
            }
          else
            {
              cout << "Number of scalar dofs " << this->mesh_num.GetNbDof() << endl;
              cout << "Number of dofs for the overall problem : " << this->GetNbDof() << endl;
            }
        }
    
    this->var_chrono.Stop("MeshGeneration");

    // selection of solver 
    Dimension dim; Symmetric property;
    glob_solver->SelectOptimalLinearSolver(this->mesh_num.GetOrder(), this->nodl, dim, property);
    
    if (this->GetNbProcPerMode() == 1)
      {
        if (this->print_level >= 2)
          this->mesh.Write(this->DOSSIER_output+string("test.mesh"));    
      }
    else
      {
        if (this->print_level >= 6)
	  this->mesh.Write(this->DOSSIER_output+string("test")
			   +to_str(rank_proc)+".mesh");    
      }
    
    if (false)
      {
      	this->ComputeMassMatrix();
    	DistributedMatrix<Complexe, General, ArrayRowSparse> Mh;
        DistributedMatrix<Complexe, General, ArrayRowSparse> Kh;
        DistributedMatrix<Complexe, General, ArrayRowSparse> Sh;
        
	// calcul matrice de masse
	GlobalGenericMatrix<Complexe> nat_mat;
	nat_mat.SetCoefMass(1.0);
	nat_mat.SetCoefStiffness(0.0);
	nat_mat.SetCoefDamping(0.0);
	
	Mh.Clear();
	this->AddMatrixWithBC(Mh, nat_mat);

        nat_mat.SetCoefMass(0.0);
	nat_mat.SetCoefStiffness(1.0);
	nat_mat.SetCoefDamping(0.0);
	
	Kh.Clear();
	this->AddMatrixWithBC(Kh, nat_mat);

        nat_mat.SetCoefMass(0.0);
	nat_mat.SetCoefStiffness(0.0);
	nat_mat.SetCoefDamping(1.0);
	
	Sh.Clear();
	this->AddMatrixWithBC(Sh, nat_mat);

        Mh.WriteText("Mh.dat");
        Sh.WriteText("Sh.dat");
        Kh.WriteText("Kh.dat");
        return;
    }
    

    EigenProblemMontjoie<Complexe>* eigen_solver;
    eigen_solver = GetNewEigenSolver(*glob_solver);
    ReadInputFile(lines_data_file, *eigen_solver);

    PolynomialEigenProblemMontjoie<Complexe>* poly_eigen_solver;
    poly_eigen_solver = GetNewPolynomialEigenSolver(*glob_solver);
    ReadInputFile(lines_data_file, *poly_eigen_solver);

    if (eigen_solver->GetNbAskedEigenvalues() > 0)
      {
	this->ComputeMassMatrix();
	eigen_solver->ComputeEigenModes();
	
	// computation successful, we exit the function
	delete glob_solver;
	delete eigen_solver;
	delete poly_eigen_solver;
        return;
      }

    if (poly_eigen_solver->GetNbAskedEigenvalues() > 0)
      {
	this->ComputeMassMatrix();
	poly_eigen_solver->ComputeEigenModes();
	
	// computation successful, we exit the function
	delete glob_solver;
	delete eigen_solver;
	delete poly_eigen_solver;
        return;
      }

    this->var_chrono.SetMessage("JacobianMatricesComputation", "compute jacobian matrices");
    this->var_chrono.Start("JacobianMatricesComputation");
    
    // computation of Ji, DFi on the mesh
    this->ComputeMassMatrix();

    this->var_chrono.Stop("JacobianMatricesComputation");
    if (this->print_level >= 1)
      this->var_chrono.DisplayTime("JacobianMatricesComputation");

    VectComplexe x_sol;
    // this vector will contain source and solution on modes (both)
    Vector<VectComplexe> glob_source(this->GetNbModes());
    
    // loop over modes (symmetric domain, cyclic domain, axisymmetric domain, etc)
    bool factorize_matrix = true;
    bool precond_matrix = false;
    VectReal_wp norme_rhs(this->GetNbModes());
    norme_rhs.Fill(0);
    int init_mode = 0, nb_modes_treated = 0;
    this->last_output_solution = false;
    if (this->GetNbProcPerMode() == 1)
      GetParallelDistributionPoints(nb_proc, rank_proc, this->GetNbModes(), nb_modes_treated, init_mode);
    else
      {
        init_mode = 0;
        nb_modes_treated = this->GetNbModes();
      }
    
    // initialisation of timers contained in the loop
    this->var_chrono.SetMessage("RightHandSideComputation", "compute the right hand side");
    this->var_chrono.SetMessage("SolutionComputation", "compute the solution of linear system");
    this->var_chrono.SetMessage("Output", "write the solution on asked outputs");
    for (int nmode = init_mode; nmode < init_mode+nb_modes_treated; nmode++)
      {
	// considered mode number
	this->number_mode = this->GetModeNumber(nmode);
        
        // computation of phase for quasi-periodic condition
        this->ComputeQuasiPeriodicPhase();
        
        // computation of source 
        this->var_chrono.Start("RightHandSideComputation");
        if (this->ModesNotStored())
          {
            if (rank_proc == 0)
              if (this->print_level >= 2)
                cout << "Right hand side is computed for mode " << this->number_mode << endl;

            x_sol.Zero();
            //if (nmode%this->nb_modes_periodic_xyz(2) == 0)
            //{
            this->ComputeRightHandSide(x_sol);
            factorize_matrix = true;

            //  }
            //else
            //factorize_matrix = false;
            
            Real_wp norme_max = Norm2(x_sol);
#ifdef SELDON_WITH_MPI
            if (nb_proc > 1)
              {
                Vector<int64_t> xtmp;
                Real_wp max_tmp = norme_max;
                MpiAllreduce(this->comm_group_mode, &max_tmp, xtmp, &norme_max, 1, MPI_MAX);
              }
#endif
            
            //DISP(norme_max); DISP(this->threshold_rhs);
            if (norme_max <= this->threshold_rhs)
              factorize_matrix = false;
            else
              factorize_matrix = true;
          }
        else
          factorize_matrix = true;

        this->var_chrono.Stop("RightHandSideComputation");        

        GlobalGenericMatrix<Complexe> nat_mat;
	//nat_mat.SetCoefStiffness(0.25); nat_mat.SetCoefMass(0.8); nat_mat.SetCoefDamping(0.4);
        if (factorize_matrix)
          {
            // computation and factorisation of finite element matrix	
            glob_solver->PerformFactorizationStep(nat_mat, !precond_matrix);
            precond_matrix = true;
          }
	
#ifdef SELDON_WITH_MPI
	// synchronization of all processors
	MPI_Barrier(this->comm_group_mode);
#endif

	// computation of the right hand side
	if (this->OnlyOneSource())
	  {
	    // only one right hand side
            this->var_chrono.Start("RightHandSideComputation");
            if (!this->ModesNotStored())
              {
                if (nmode == 0)
                  {
                    if (rank_proc == 0)
                      if (this->print_level >= 2)
                        cout << "Right hand sides are computed for all the modes" << endl;

                    this->ComputeRightHandSide(glob_source);
		    Real_wp norme_max = 0;
                    for (int k = 0; k < glob_source.GetM(); k++)
                      {
                        norme_rhs(k) = Norm2(glob_source(k));
                        norme_max = max(norme_max, norme_rhs(k));
                      }
                    
#ifdef SELDON_WITH_MPI
                    if (this->GetNbProcPerMode() > 1)
                      {
                        Vector<int64_t> xtmp;
                        VectReal_wp norme_tmp = norme_rhs;
                        Real_wp max_tmp = norme_max;
                        MpiAllreduce(this->comm_group_mode, norme_tmp, xtmp, norme_rhs,
                                     norme_tmp.GetM(), MPI_MAX);
                        
                        MpiAllreduce(this->comm_group_mode, &max_tmp, xtmp, &norme_max, 1, MPI_MAX);
                      }
#endif
                    
                    if (norme_max != Real_wp(0))
                      for (int k = 0; k < glob_source.GetM(); k++)
                        norme_rhs(k) /= norme_max;               
                    
                  }
                
                x_sol = glob_source(nmode);
              }

	    //srand(0); x_sol.FillRand();
            //for (int i = 81; i < x_sol.GetM(); i++)
            //  x_sol(i) = 0.0;
            
            //x_sol *= Real_wp(1e-9);
	    
            this->var_chrono.Stop("RightHandSideComputation");
            
            //Vector<VectComplexe> source_vec; source_vec.SetData(1, &source_rhs);
	    //this->WriteOutputFile(source_vec, string("source")+to_str(nmode));
            //source_vec.Nullify();
            if (this->print_level >= 6)
              {
                if (this->GetNbProcPerMode() > 1)
                  {
                    this->mesh_num.GlobDofNumber_Subdomain
                      .WriteText("numGlob"+to_str(rank_proc)+".dat");
                    
                    x_sol.Write(string("rhs"+to_str(nmode)+"P"
                                            + to_str(rank_proc)+".dat"));
                  }
                else
                  x_sol.Write(string("rhs"+to_str(nmode)+".dat"));
              }
	  }
        
        // displays memory used by different objects involved in the simulation
	size_t size_sol = x_sol.GetM();
        for (int k = 0; k < glob_source.GetM(); k++)
          size_sol += glob_source(k).GetM();
        
        if (this->print_level >= 1)
          {
            size_object["VectorRhsSolution"] = size_sol*sizeof(Complexe);
            glob_solver->GetMemoryUsed(size_object);
	    this->GetMemoryUsed(size_object);            
            this->DisplayMemoryUsed(size_object);
          }

        if (rank_proc == 0)
          if (this->print_level >= 2)
            cout << "The linear system will be solved" << endl;
	
        // computes solution
	if (factorize_matrix)
          {
            this->var_chrono.Start("SolutionComputation");
            if (this->OnlyOneSource())
              {
                //Vector<Complexe> x_sol_trans(x_sol);
                // computation of the solution
		if (this->ModesNotStored())
                  glob_solver->ComputeSolution(x_sol, nat_mat);
                else
                  {
                    if (norme_rhs(nmode) > this->threshold_rhs)
                      glob_solver->ComputeSolution(x_sol, nat_mat);
                  }
                
                //glob_solver->ComputeSolution(SeldonTrans, x_sol_trans, nat_mat);
                //x_sol_trans.Write("sol_t.dat");
                
		//Vector<VectComplexe> sol_vec; sol_vec.SetData(1, &x_sol);
                //this->WriteOutputFile(sol_vec, string("sol")+to_str(nmode));
                //sol_vec.Nullify();
                
                // for a transparent condition, we have to iterate to find the correct solution
                if (solver_transparent->UseTransparentCondition())
                  {
                    solver_transparent->Init();
                    // dans le cas d'une condition transparente,
                    // on doit iterer en plus pour avoir la solution
                    Vector<Complexe> source_rhs = x_sol;
                    solver_transparent->Solve(x_sol, source_rhs);
                  }
              }
            else
              {
                cout << "Multiple right hand sides not implemented "<<endl; // exit(0);
                
                // multiple right hand sides
		abort();
	      }
            
            this->var_chrono.Stop("SolutionComputation");
          }
        
        // writing solution in glob_source
        if ( (!this->ModesNotStored()) && (this->GetNbModes() > 1))
	  glob_source(nmode) = x_sol;
	
	// clearing the factorization of the matrix
	if (factorize_matrix)
          glob_solver->ClearFactorization();

        // next computation ?
        if (!this->ModesNotStored())
          {
            if (nmode < this->GetNbModes()-1)
              {
                if (norme_rhs(nmode+1) > this->threshold_rhs)
                  factorize_matrix = true;
                else
                  factorize_matrix = false;
              }
          }
                  
        this->var_chrono.Start("Output");
        
        if (this->print_level >= 6)
          {
            if (this->GetNbProcPerMode() > 1)
              x_sol.Write(string("sol")+to_str(nmode)+"P"+to_str(rank_proc)+".dat");
            else
              x_sol.Write(string("sol")+to_str(nmode)+".dat");
          }
        
        if (nmode == init_mode + nb_modes_treated-1)
          this->last_output_solution = true;
        
        if ((this->ModesNotStored())||(this->GetNbModes() == 1))
          {   
            
            if (output_rcs_base.RcsToBeComputed())
              output_rcs_base.ComputeRCS(x_sol);            

            // writing the solution on the asked outputs
            if ( (factorize_matrix) || (this->last_output_solution))
              this->WriteDatas(x_sol);
          }
        
        this->var_chrono.Stop("Output");
                	
      }

            
    if (rank_proc == 0)
      if (this->print_level >= 0)
        cout << rank_proc << "Solution computed, now post-processing is performed"<<endl;

    this->var_chrono.Start("Output");

    // on calcule aussi la ser / rcs (section equivalente radar / radar cross section)
    if ( (!this->ModesNotStored()) && (this->GetNbModes() > 1))
      {
        if (output_rcs_base.RcsToBeComputed())
          output_rcs_base.ComputeRCS(x_sol);
	
        this->WriteDatas(glob_source);
      }
    
    this->var_chrono.Stop("Output");
        
    // displaying time computations
    this->var_chrono.Stop("GlobalComputation");
    if (this->print_level >= 1)
      {
        this->var_chrono.DisplayAll();
        this->DisplayMemoryUsed(size_object);
      }
    
    if (rank_proc == 0)
      if (this->print_level >= 0)
	cout << rank_proc << " The variables are destructed"<<endl;
    
    delete glob_solver;

    delete solver_transparent;

  } // end method RunAll


  //! adds size used by different objects in VarHarmonic
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::GetMemoryUsed(map<string, size_t>& var) const
  {
    VarProblem<Dimension>::GetMemoryUsed(var);
    VarBoundaryCondition<Complexe, Dimension>::GetMemoryUsed(var);
    VarOutputProblem<Dimension>::GetMemoryUsed(var);
    VarSourceProblem<Complexe, Dimension>::GetMemoryUsed(var);    
  }
  
  
  //! returns size of memory used by the object in bytes
  template<class Complexe, class Dimension>
  size_t VarHarmonic_Base<Complexe, Dimension>::GetMemorySize() const
  {
    // using GetMemoryUsed of leaf class
    // such that GetMemorySize() does not need to be overloaded in derived classes
    map<string, size_t> var;
    this->GetMemoryUsed(var);
    map<string, size_t>::iterator it;
    size_t taille = 0;
    for (it = var.begin(); it != var.end(); ++it)
      taille += it->second;
    
    return taille;
  }
  
  
  //! computation of the inverse square root of the mass matrix
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::GetInverseSquareRootMassMatrix(VectComplexe& Dh)
  {
    this->GetMassMatrix(Dh);
    
    for (int i = 0; i < Dh.GetM(); i++)
      Dh(i) = 1.0/sqrt(abs(Dh(i)));
  }


  //! copies input parameters of another similar problem
  template<class Complexe, class Dimension>
  void VarHarmonic_Base<Complexe, Dimension>
  ::CopyInputData(const VarProblem_Base& var_base)
  {
    const VarHarmonic_Base<Complexe, Dimension>& var
      = static_cast<const VarHarmonic_Base<Complexe, Dimension>& >(var_base);
    
    VarProblem<Dimension>::CopyInputData(var);
    VarBoundaryCondition<Complexe, Dimension>::CopyInputData(var);
    VarOutputProblem<Dimension>::CopyInputData(var);
    VarSourceProblem<Complexe, Dimension>::CopyInputData(var);
  }
  
}

#define MONTJOIE_FILE_VAR_HARMONIC_BASE_CXX
#endif

