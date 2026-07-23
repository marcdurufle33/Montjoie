#ifndef MONTJOIE_FILE_VAR_INSTATIONARY_CXX

namespace Montjoie
{

  //! default constructor;
  MassMatrixParameter::MassMatrixParameter()
  {
    store_mass_matrix = false;
    iterative_solver_mass = false;
    stopping_criterion = 1e-10;
    nb_max_iteration = 100;
    print_level = 0;
  }
  

  //! Sets values of attributes to default values
  void VarInstationary_Base::InitDefaultValues()
  {
    glob_solver = NULL;
    Glob_mat_Dh = NULL;
    Glob_mat_Kh = NULL;

    time_step_to_be_computed = true;
    
    deltat = 0.01;
    initial_time = 0; final_time = 0;
    load_reprise = false; save_reprise = false;
    load_iter_reprise = 0; save_iter_reprise = 0;
    path_reprise = "./";
    bascule_reprise = -1;
    max_norme_sol = 1e30;
    
    file_with_local_dt = string();
    method_computation_cfl = CFL_EIGENVALUE;
    global_level_operator = level_time_scheme.ALL_LEVELS;

    tlimit_source = 0.0; 
    tinit_source = 0;
    default_threshold_cfl = 1e-4;
    epsilon_time = 1e-10;
    same_preconditioning = false;
    first_preconditioning_computed = false;
    ptr_real_precond = NULL;
    nb_iter_per_display = -1;
    display_computational_time = false;
    check_prod_stiffness_matrix = false;

    small_random_initial = false;
  }
  
  
  //! Destructor
  VarInstationary_Base::~VarInstationary_Base()
  {
    if (Glob_mat_Kh != NULL)
      delete Glob_mat_Kh;
    
    if (Glob_mat_Dh != NULL)
      delete Glob_mat_Dh;
    
    if (glob_solver != NULL)
      delete glob_solver;
    
    for (int i = 0; i < vec_solver.GetM(); i++)
      if (vec_solver(i) != NULL)
	delete vec_solver(i);
  }
  
  
  //! Reading of parameters in data file for instationnary problems only
  /*!
    \param[in] description_field keyword of the input file
    \param[in] parameters list of parameters associated to description_field
    \param[in] nb_param number of parameters
    \return 0 if the field has been found, -1 otherwise
  */
  void VarInstationary_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("TimeStep"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "TimeStep needs more parameters, for instance :" << endl;
	    cout << "TimeStep = dt" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("AUTO"))
	  {
            // in this case, the time step is computed automatically
            // through computation of maximal eigenvalues of the matrix
	    time_step_to_be_computed = true;

            // you can provide a guess about the initial time (used to compute the actual value)
            if (parameters.GetM() > 1)
              deltat = to_num<Real_wp>(parameters(1));
            else
              deltat = 1.0;
	  }
	else if (!parameters(0).compare("DIVISION"))
	  {
            if (parameters.GetM() <= 0)
              {
                cout << "In SetInputData of VarInstationary" << endl;
                cout << "TimeStep needs more parameters, for instance :" << endl;
                cout << "TimeStep = DIVISION Tfinal nb_iterations" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
	    time_step_to_be_computed = false;
            // nominal time step
	    deltat = to_num<Real_wp>(parameters(1))/to_num<Real_wp>(parameters(2));;
            if (parameters.GetM() > 4)
              {
                // for a local time step, you can provide
                // the time step for each element in a file
                if (!parameters(3).compare("FILE"))
                  {
                    file_with_local_dt = parameters(4);
                  }
              }
	  }
	else
	  {
	    time_step_to_be_computed = false;
            // nominal time step
	    deltat = to_num<Real_wp>(parameters(0));
            if (parameters.GetM() > 2)
              {
                // for a local time step, you can provide
                // the time step for each element in a file
                if (!parameters(1).compare("FILE"))
                  {
                    file_with_local_dt = parameters(2);
                  }
              }
	  }
      }
    else if (description_field == "RandomInitialCondition")
      {
        if (parameters(0) == "Small")
          small_random_initial = true;
        else
          small_random_initial = false;
      }
    else if (description_field == "DisplayRate")
      {
        nb_iter_per_display = to_num<int>(parameters(0));
        if (parameters(1) == "YES")
          display_computational_time = true;
        else 
          display_computational_time = false;
      }
    else if (!description_field.compare("EigenvalueTolerance"))
      {
        this->default_threshold_cfl = to_num<Real_wp>(parameters(0));
      }
    else if (description_field == "MethodEvaluationCFL")
      {
	if (parameters(0) == "MeshSize")
	  method_computation_cfl = CFL_MESH_SIZE_ELT;
	else if (parameters(0) == "MeshSizeFace")
	  method_computation_cfl = CFL_MESH_SIZE_FACE;
	else if (parameters(0) == "MeshSizeVertex")
	  method_computation_cfl = CFL_MESH_SIZE_VERTEX;
	else if (parameters(0) == "Eigenvalue")
	  method_computation_cfl = CFL_EIGENVALUE;
	else
	  {
	    cout << "Unknown method" << endl;
	    abort();
	  }
      }
    else if (!description_field.compare("TimeInterval"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "TimeInterval needs more parameters, for instance :" << endl;
	    cout << "TimeInterval = t0 tf" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	initial_time = to_num<Real_wp>(parameters(0));
	final_time = to_num<Real_wp>(parameters(1));
      }
    else if (!description_field.compare("OrderTimeScheme"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "OrderTimeScheme needs more parameters, for instance :" << endl;
	    cout << "OrderTimeScheme = r type_scheme" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        param_time_scheme = parameters;
        var_time_scheme.SetInputData(description_field, parameters);
      }
    else if (description_field == "SamePreconditioning")
      {
	if (parameters(0) == "YES")
	  same_preconditioning = true;
	else
	  same_preconditioning = false;
      }
    else if (description_field == "ExplicitMassMatrix")
      {
	if (parameters(0) == "YES")
	  param_mass.store_mass_matrix = true;
	else
	  param_mass.store_mass_matrix = false;
      }
    else if (description_field == "MassMatrixSolver")
      {
	if (parameters(0) == "Direct")
	  {
            param_mass.iterative_solver_mass = false;
	    if (parameters.GetM() >= 2)
	      param_mass.print_level = to_num<int>(parameters(1));
          }
	else if (parameters(0) == "Iterative")
	  {
	    param_mass.iterative_solver_mass = true;
	    if (parameters.GetM() >= 3)
	      {
		param_mass.stopping_criterion = to_num<Real_wp>(parameters(1));
		param_mass.nb_max_iteration = to_num<int>(parameters(2));
	      }

	    if (parameters.GetM() >= 4)
	      param_mass.print_level = to_num<int>(parameters(3));
	  }
	else
	  {
	    cout << "Unknown solver" << endl;
	    abort();
	  }
      }
    else if (!description_field.compare("FileOutputEnergy"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "FileOutputEnergy needs more parameters, for instance :" << endl;
	    cout << "FileOutputEnergy = file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        var_time_scheme.SetInputData(description_field, parameters);
      }
    else if (!description_field.compare("ParametersOutputEnergy"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "ParametersOutputEnergy needs more parameters, for instance :" << endl;
	    cout << "ParametersOutputEnergy = t_begin_ener t_end_ener deltat_ener (size_buffer)" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        var_time_scheme.SetInputData(description_field, parameters);
      }
    else if (!description_field.compare("InitialCondition"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "InitialCondition needs more parameters, for instance :" << endl;
	    cout << "InitialCondition = USER" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("USER"))
	  {
	    param_initial_condition.Reallocate(parameters.GetM()-1);
	    for (int i = 1; i < parameters.GetM(); i++)
	      param_initial_condition(i-1) = parameters(i);
	    
	  }
      }
    else if (!description_field.compare("TemporalSource"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "TemporalSource needs more parameters, for instance :" << endl;
	    cout << "TemporalSource = RICKER" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	param_time_source = parameters;
      }
    else if (!description_field.compare("LoadReprise"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "LoadReprise needs more parameters, for instance :" << endl;
	    cout << "LoadReprise = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (parameters(0) == "YES")
          load_reprise = true;
        else
          load_reprise = false;
      }
    else if (!description_field.compare("SaveReprise"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "SaveReprise needs more parameters, for instance :" << endl;
	    cout << "SaveReprise = num" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	save_reprise = true;
	save_iter_reprise = to_num<int>(parameters(0));
      }
    else if (description_field == "RegroupReprise")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "RegroupReprise needs one parameter, for instance :" << endl;
	    cout << "RegroupReprise = N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        int nb_proc_per_node = to_num<int>(parameters(0));
        output_reprise.RegroupWritingOnMaster(nb_proc_per_node);
      }
    else if (!description_field.compare("PathReprise"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "PathReprise needs more parameters, for instance :" << endl;
	    cout << "PathReprise = path" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
        path_reprise = parameters(0);
      }    
    else if (!description_field.compare("NormeMaxSolution"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarInstationary" << endl;
	    cout << "NormeMaxSolution needs more parameters, for instance :" << endl;
	    cout << "NormeMaxSolution = norme_max" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        max_norme_sol = to_num<Real_wp>(parameters(0));
      }
  }


  //! fills how memory is used by the current object
  void VarInstationary_Base::GetMemoryUsed(map<string, size_t>& var) const
  {
    var_problem.GetMemoryUsed(var);
    if ((glob_solver != NULL) && glob_solver->FactorizationCompleted())
      glob_solver->GetMemoryUsed(var);

    for (int i = 0; i < vec_solver.GetM(); i++)
      vec_solver(i)->GetMemoryUsed(var, true);
    
    size_t taille = 0;
    if (Glob_mat_Kh != NULL)
      taille += Glob_mat_Kh->GetMemorySize();

    var["StiffnessMatrix"] = taille + Glob_mat_RhS.GetMemorySize() +
      Glob_mat_RhV.GetMemorySize() + Glob_mat_Ch_Lambda.GetMemorySize();
    
    if (Glob_mat_Dh != NULL)
      var["MassMatrix"] = Glob_mat_Dh->GetMemorySize();

    var["SourceVector"] = sparse_vector_source.GetMemorySize() + lambda_vector_source.GetMemorySize()
      + dirichlet_source.GetMemorySize() + dirichlet_source_unassembled.GetMemorySize();

    var["VectorsTimeScheme"] = var_time_scheme.GetMemorySize();

    // if heavy objects are stored in var_time_scheme
    var_time_scheme.GetMemoryUsed(var);
  }
  
  
  //! allocates a VirtualTimeSource object with param_time_source
  VirtualTimeSource<Real_wp>* VarInstationary_Base::GetNewTimeSource()
  {
    const VectString& parameters = param_time_source;
    tinit_source = 0.0;
    tlimit_source = 0.0;
    if (parameters.GetM() <= 0)
      return NULL;
    
    Real_wp freq = var_problem.GetFrequency();
    if (!parameters(0).compare("RICKER"))
      {
	tlimit_source = 2.5/freq;
	return new TimeRickerSource(freq);
      }
    else if (!parameters(0).compare("DERIVATIVE_RICKER"))
      {
	tlimit_source = 2.5/freq;
	return new DerivativeTimeRickerSource(freq);
      }
    else if (!parameters(0).compare("MODIFIED_RICKER"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "TemporalSource needs more parameters, for instance :" << endl;
	    cout << "TemporalSource = MODIFIED_RICKER t0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	Real_wp center_time = to_num<Real_wp>(parameters(1));
	tinit_source = (-2.5+center_time)/(sqrt(Real_wp(2))*freq);
	tlimit_source = (2.5+center_time)/(sqrt(Real_wp(2))*freq);
	return new TimeModifiedRickerSource(freq, center_time);
      }
    else if (!parameters(0).compare("GAUSSIAN"))
      {
	Real_wp width_fct = 0;
	if (parameters.GetM() > 1)
	  width_fct = to_num<Real_wp>(parameters(1));
	    
	tlimit_source = 2.5/freq + width_fct;
	return new TimeGaussianSource(freq, 0.0, width_fct);
      }
    else if (!parameters(0).compare("HARMONIC"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "TemporalSource needs more parameters, for instance :" << endl;
	    cout << "TemporalSource = HARMONIC tf" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	tlimit_source = to_num<Real_wp>(parameters(1));
	return new TimeHarmonicSource(freq, tlimit_source);
      }
    else if (!parameters(0).compare("SINUS_GAUSSIAN"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "TemporalSource needs more parameters, for instance :" << endl;
	    cout << "TemporalSource = SINUS_GAUSSIAN alpha" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	Real_wp gaussian_time_troncature = to_num<Real_wp>(parameters(1));
	// toff = sqrt(log(1e6)/b) chosen such as the function begins at 1e-6 for t=0
	tlimit_source = 2.0*sqrt(log(1e6)/gaussian_time_troncature);
	return new TimeSinusGaussianSource(freq, 0.5*tlimit_source, gaussian_time_troncature);
      }
    else if (!parameters(0).compare("MODULATED_CARDINAL_SINUS"))
      {
	Real_wp bandwidth(0), center_time(0), correlation(0);
	if (parameters.GetM() == 4)
	  {
	    bandwidth = to_num<Real_wp>(parameters(1));
	    center_time = to_num<Real_wp>(parameters(2));
	    correlation = to_num<Real_wp>(parameters(3));
	  }
	else
	  {
	    bandwidth = freq;
	    correlation = 3. / freq;
	    center_time = 2. * correlation;
	  }
	
	tlimit_source = 1e300;
	return new TimeModulatedCardinalSinusSource(center_time, freq, bandwidth, correlation);
      }
    else if (!parameters(0).compare("RANDOM_MODULATED_CARDINAL_SINUS"))
      {
	tinit_source = to_num<Real_wp>(parameters(1));
	tlimit_source = to_num<Real_wp>(parameters(2));
	Real_wp tdelta_source = to_num<Real_wp>(parameters(3));
	Real_wp bandwidth(0), center_time(0), correlation(0);
	if (parameters.GetM() == 7)
	  {
	    bandwidth = to_num<Real_wp>(parameters(4));
	    center_time = to_num<Real_wp>(parameters(5));
	    correlation = to_num<Real_wp>(parameters(6));
	  }
	else
	  {
	    bandwidth = freq;
	    correlation = 3. / freq;
	    center_time = 2. * correlation;
	  }
	
	return new TimeRandomModulatedCardinalSource(freq, bandwidth, correlation, center_time,
						     tinit_source, tlimit_source, tdelta_source, false);
      }
    else if (!parameters(0).compare("MODULATED_CARDINAL_COSINUS"))
      {
	Real_wp bandwidth(0), center_time(0), correlation(0);
	if (parameters.GetM() == 4)
	  {
	    bandwidth = to_num<Real_wp>(parameters(1));
	    center_time = to_num<Real_wp>(parameters(2));
	    correlation = to_num<Real_wp>(parameters(3));
	  }
	else
	  {
	    bandwidth = freq;
	    correlation = 3. / freq;
	    center_time = 2. * correlation;
	  }

	tlimit_source = 1e300;
	return new TimeModulatedCardinalCosinusSource(center_time, freq, bandwidth, correlation);
      }
    else if (!parameters(0).compare("RANDOM_MODULATED_CARDINAL_COSINUS"))
      {
	tinit_source = to_num<Real_wp>(parameters(1));
	tlimit_source = to_num<Real_wp>(parameters(2));
	Real_wp tdelta_source = to_num<Real_wp>(parameters(3));
	Real_wp bandwidth(0), center_time(0), correlation(0);
	if (parameters.GetM() == 7)
	  {
	    bandwidth = to_num<Real_wp>(parameters(4));
	    center_time = to_num<Real_wp>(parameters(5));
	    correlation = to_num<Real_wp>(parameters(6));
	  }
	else
	  {
	    bandwidth = freq;
	    correlation = 3. / freq;
	    center_time = 2. * correlation;
	  }
	
	return new TimeRandomModulatedCardinalSource(freq, bandwidth, correlation, center_time,
						     tinit_source, tlimit_source, tdelta_source, true);
      }
    else if (!parameters(0).compare("FILE"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "TemporalSource needs more parameters, for instance :" << endl;
	    cout << "TemporalSource = FILE file_name t0 tf" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	string file_source_time = parameters(1);
	tinit_source = to_num<Real_wp>(parameters(2));
	tlimit_source = to_num<Real_wp>(parameters(3));
	return new TimeFileSource(freq, file_source_time, tinit_source, tlimit_source);
      }
    else if (!parameters(0).compare("USER"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarInstationary_Base" << endl;
	    cout << "TemporalSource needs more parameters, for instance :" << endl;
	    cout << "TemporalSource = USER tf" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	tlimit_source = 1e300;
	return new TimeUserSource(freq, tlimit_source, parameters);
      }
    
    return NULL;
  }
  
  
  //! returns true if the source is separated in time and in space g(x) h(t)
  bool VarInstationary_Base::SourceDoesNotDependOnTime() const
  {
    if (var_source.GetSourceType(0) == var_source.SRC_TOTAL_FIELD)
      return false;

    if (var_source.GetSourceType(0) == var_source.SRC_DIFFRACTED_FIELD)
      return false;
    
    return true;
  }
  
  
  //! Computation of the source terms
  /*!
    In the case of source in separated variables f(x) g(t)
    Vector \f$ (F)_i = \int f \varphi_i \f$ is computed before time iterations
    F is stored in sparse format in arrays dof_source/value_source_ddl
    so that F(dof_source(i)) = value_source_ddl(i)
  */
  void VarInstationary_Base::ComputeRightHandSide()
  {
    VectReal_wp space_source;
    // computation of F = space_source as a full Vector
    // implementation of this method in DefineSource_Elliptic.cxx
    if (!this->SourceDoesNotDependOnTime())
      {
	sparse_vector_source.Clear();
	return;
      }
    
    var_source.ComputeRightHandSide(space_source, false);
        
    //var_harmonic.WriteOutputFile(space_source, string("sourceH"));
    // isolating source associated with dirichlet condition
    int nb_dir = var_boundary.GetNbDirichletDof();
    const IVect& DirDof = var_boundary.GetDirichletDofNumber();
    dirichlet_source.Reallocate(nb_dir);
    VectReal_wp y(space_source); y.Zero();
    for (int i = 0; i < nb_dir; i++)
      {
	dirichlet_source.Index(i) = DirDof(i);
 	dirichlet_source.Value(i) = space_source(DirDof(i));
	y(DirDof(i)) = space_source(DirDof(i));
	space_source(DirDof(i)) = 0;
      }

    // unassembled and assembled Dirichlet are needed
    dirichlet_source_unassembled = dirichlet_source;
    var_comm.AssembleDirichlet(y);
    
    for (int i = 0; i < nb_dir; i++)
      dirichlet_source.Value(i) = y(DirDof(i));
    
    // removing small values
    Real_wp threshold = epsilon_machine*Norm2(space_source);
    for (int i = 0; i < space_source.GetM(); i++)
      if (abs(space_source(i)) <= threshold)
	space_source(i) = 0;

    // conversion to a sparse vector
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	VectReal_wp sourceL, sourceUV;
	int nb_dof_L = var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
	int nb_dof_uv = var_problem.GetNbDof() - nb_dof_L;
	sourceL.SetData(nb_dof_L, space_source.GetData());
	sourceUV.SetData(nb_dof_uv, &space_source(nb_dof_L));

	Copy(sourceUV, sparse_vector_source);
	Copy(sourceL, lambda_vector_source);
	sourceL.Nullify(); sourceUV.Nullify();
      }
    else
      Copy(space_source, sparse_vector_source);
  }

  
  //! fills the spatial source (without time dependence
  void VarInstationary_Base::FillSource(VectReal_wp& Bn) const
  {
    Bn.Fill(0);
    for (int i = 0; i < sparse_vector_source.GetM(); i++)
      Bn(sparse_vector_source.Index(i)) = sparse_vector_source.Value(i);
  }
  
  
  //! computation of temporal source f(t) 
  /*!
    \param[in] t0 time where source has to be computed
    \param[in] f0 frequency
    \param[in] nb_deriv we want to compute d^n f/dt^n
    \param[out] pulse result
  */
  void VarInstationary_Base
  ::SourceOnlyTime(const Real_wp& t0, int nb_deriv, Real_wp& pulse)
  {
    time_source.EvaluateDerivative(t0, nb_deriv, pulse);
  }
  
  
  //! Inhomogeneous Dirichlet condition
  void VarInstationary_Base
  ::SetDirichletCondition(const Real_wp& t, int nb_deriv, VectReal_wp& b_src, Real_wp alpha)
  {
    // first setting to 0 for homogeneous Dirichlet
    var_boundary.ImposeNullDirichletCondition(b_src);
    
    // then taking into account inhomogeneous Dirichlet terms
    if (SourceDoesNotDependOnTime())
      {
	int nb_dir = dirichlet_source.GetM();
	if (nb_dir <= 0)
	  return;

	Real_wp pulse(0); SourceOnlyTime(t, nb_deriv, pulse); 
	pulse *= alpha;
        for (int i = 0; i < dirichlet_source.GetM(); i++)
	  b_src(dirichlet_source.Index(i)) = pulse*dirichlet_source.Value(i);
      }
    else
      {
	SetDirichletConditionSource(t, nb_deriv, b_src, alpha);
	var_comm.AssembleDirichlet(b_src);
      }
  }


  //! adds the source :  b_src <- b_src + alpha*df/dt^nb_deriv
  /*!
    \param[in] t current time where we want to evaluate the source
    \param[in] alpha coefficient 
    \param[in] nb_deriv 0 for the primitive of the source f,
    1 for the source, 2 for the derivative, etc
    \param[in,out] b_src vector to be modified
  */
  void VarInstationary_Base
  ::AddPrimitiveSourceAtTime(const Real_wp& alpha, const Real_wp& t,
			     int nb_deriv, VectReal_wp& b_src)
  {
    if (t >= tlimit_source)
      return;

    if (SourceDoesNotDependOnTime())
      {
        Real_wp pulse(0);
        SourceOnlyTime(t, nb_deriv, pulse);
        //DISP(t); DISP(pulse); DISP(sparse_vector_source.GetM());
        Real_wp alpha_ = alpha*pulse;
        for (int i = 0; i < sparse_vector_source.GetM(); i++)
          b_src(sparse_vector_source.Index(i)) += alpha_*sparse_vector_source.Value(i);
      }
    else
      {
	VectReal_wp b_scal, b_vec;
	int Ns = this->GetNbScalarUnknowns();
	b_scal.SetData(Ns, &b_src(0));
	b_vec.SetData(b_src.GetM() - Ns, &b_src(Ns));
	
	AddScalarSourceAtTime(alpha, t, nb_deriv, b_scal);
	AddVectorialSourceAtTime(alpha, t, nb_deriv, b_vec);

	b_scal.Nullify();
	b_vec.Nullify();
      }
  }


  //! Adds the source for time-schemes
  void VarInstationary_Base
  ::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			   int nb_deriv, VectReal_wp& b_src)
  {
    VectReal_wp b_dir(var_boundary.GetNbDirichletDof());
    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      b_dir(i) = b_src(var_boundary.GetDirichletDofNumber(i));
    
    // usual space source
    AddPrimitiveSourceAtTime(alpha, t, nb_deriv, b_src);

    // and inhomogeneous Dirichlet
    SetDirichletConditionSource(t, nb_deriv+1, b_src);

    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      {
	Real_wp val = b_src(var_boundary.GetDirichletDofNumber(i));
	b_src(var_boundary.GetDirichletDofNumber(i)) = b_dir(i) + alpha*val;
      }
  }

  
  //! writing of Uh on a file,whose name changes according nb_iter (inst0001.dat inst0002.dat ...)
  /*!
    \param[in] Uh vector solution
    \param[in] t time
    \param[in,out] nb_iter number of iterations done before
    \remark this method increments the iteration number
  */
  void VarInstationary_Base
  ::WriteSnapshot(int nb_iter, const Real_wp& t, VectReal_wp& Uh, bool compute_time)
  {
    SaveDataReprise(nb_iter, Uh);
    
    if ((nb_iter == 0) && (var_problem.print_level >= 2))
      {
        if (compute_time)
          {
            glob_chrono.SetMessage("TimeIteration", "perform time iterations");
            glob_chrono.Start("TimeIteration");
          }

	// variable storing size of different objects
	map<string, size_t> size_object;
	
	// usage memory is displayed before the first call to Advance
	this->GetMemoryUsed(size_object);
	var_comm.DisplayMemoryUsed(size_object);
      }

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    var_output.WriteSnapshot(nb_iter, t, Uh, 1);
    
    if (nb_iter%nb_iter_per_display == 0)
      {
        Real_wp norme_sol = this->GetNormeSolution(Uh);
        if (rank_proc == 0)
          {
            cout<<" At iteration "<<nb_iter<<" t = "<<t<<endl;
            cout<<" Norme solution "<<norme_sol<<endl;
          }
        
        if (compute_time)
          {
            glob_chrono.Stop("TimeIteration");
            glob_chrono.DisplayAll();
	    double dt_loc, dt_sum, dt_min, dt_max;
	    glob_chrono.GetGlobalSeconds("TimeIteration", dt_loc, dt_sum, dt_min, dt_max);
	    double dt_comm;
	    glob_chrono.GetGlobalSeconds(VirtualTimer::COMM, dt_loc, dt_comm, dt_min, dt_max);
	    if(rank_proc == 0)
	      {
		cout << "Efficacite = " << 100.0*(dt_sum-dt_comm)/dt_sum << " %" << endl;
		cout << "Temps minimal d'attente = " << dt_min << " s" << endl;
		cout << "Temps maximal d'attente = " << dt_max << " s" << endl;
	      }
	    glob_chrono.Start("TimeIteration");
	  }
      }
  }
  
  
  //! loads files stored on the disk (reprise)
  void VarInstationary_Base::LoadDataReprise(VectReal_wp& Xh)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc, nb_proc;
    MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    // on lit les donnees de la reprise
    IVect num_iter;
    string prefix_load = this->path_reprise;
    string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
    if (rank_proc == 0)
      {
        num_iter.Read(prefix_load+"IterationAcousSave"+suffix);
        if (num_iter.GetM() != 3)
          {
            // on pense alors que l'ecriture critique a ete interrompu,
            // on regarde l'autre fichier dans ce cas
            num_iter.Read(prefix_load+"IterationAcousSaveNew"+suffix);
                
            if (num_iter.GetM() != 3)
              {
                cout << "Echec de la reprise" << endl;
                cout << "Veuillez corriger le probleme" << endl;
                abort();
              }
          }
      }
    else
      num_iter.Reallocate(3);
    
    // on broadcaste num_iter sur tous les procs
#ifdef SELDON_WITH_MPI
    MPI_Bcast(num_iter.GetData(), num_iter.GetM(), MPI_INTEGER, 0,
              var_comm.comm_group_mode);
#endif
    
    this->load_iter_reprise = num_iter(0);

    suffix = to_str(num_iter(2)) + "_P" + suffix;
    output_reprise.StartReading(path_reprise + "SolAcous" + suffix);
    
    output_reprise.Read(Xh);
    
    output_reprise.FinalizeReading();
    
    Real_wp t0 = this->initial_time + num_iter(0)*this->deltat;
    int nt = num_iter(0);
    for (int i = 0; i < var_output.output_grid_param.GetM(); i++)
      var_output.output_grid_param(i).ChangeTime(nt, t0);
    
    for (int i = 0; i < var_output.output_mesh_param.GetM(); i++)
      var_output.output_mesh_param(i).ChangeTime(nt, t0);
  }


  //! saves files stored on the disk (reprise)
  void VarInstationary_Base::SaveDataReprise(int nt, const VectReal_wp& Xh)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc, nb_proc;
    MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    if (save_reprise && (nt%save_iter_reprise == 0))
      {
        string suffix = to_str(rank_proc) + "N"+to_str(nb_proc)+".dat";
        string prefix = path_reprise;
        bool first_save = false;
        if (bascule_reprise == -1)
          {
            // premiere fois qu'on imprime
            bascule_reprise = 0;
            first_save = true;
          }
            
        IVect num_iter(3);
        int num_dossier = -1;
        num_iter(0) = nt;
        num_iter(1) = num_dossier;
        num_iter(2) = bascule_reprise;
        
        if (rank_proc == 0)
          {
            if (first_save)
              num_iter.Write(prefix+"IterationAcousSaveNew"+suffix);
            else
              {
                // ecriture critique ou on affirme que les 
                // anciens iteres sauves sont ceux qu'il faut considerer
                // pour une reprise
                IVect num_iter_old;
                num_iter_old.Read(prefix+"IterationAcousSaveNew"+suffix);
                num_iter_old.Write(prefix+"IterationAcousSave"+suffix);
                
                // fin de l'ecriture critique                  
                num_iter.Write(prefix+"IterationAcousSaveNew"+suffix);
              }
          }
        
        // on ecrit ici les nouveaux iteres (qui ne seront pas consideres pour la reprise)
        // le programme peut etre interrompu durant l'ecriture sans dommage collateral
        suffix = to_str(num_iter(2)) + "_P" + suffix;          
        output_reprise.StartWriting(prefix + "SolAcous" + suffix);
        
        output_reprise.Write(Xh);
       
        output_reprise.FinalizeWriting();

        bascule_reprise = 1 - bascule_reprise;

#ifdef SELDON_WITH_MPI          
        MPI_Barrier(var_comm.comm_group_mode);
#endif
        
      }
  }

  
  //! initialization of variables before running time iterations
  void VarInstationary_Base::InitTimeIterations()
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    // initialization of temporal sources
    this->epsilon_time = (final_time-initial_time)*1e-15;    
    this->InitSource();
    
    this->level_time_scheme.SetNbElt(var_computation.GetNbElt(), var_boundary.GetNbEltPML());
    if (!this->load_reprise)
      var_output.InitOutput(this->initial_time);
    
    // initialization of file names for vectorial unknowns
    if ( var_problem.FirstOrderFormulation()
	 && (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS) )
      {        
        for (int i = 0; i < var_output.output_grid_param.GetM(); i++)
          {
            var_output.output_grid_param(i).ResizeNbFile(2+dim_N);
            for (int j = 0; j < dim_N; j++)
              {
                string name_file = var_output.output_grid_param(i).GetFileName(1);
                name_file = GetBaseString(name_file) + "_V" 
                  + to_str(j) + "_." + GetExtension(name_file);
                
                var_output.output_grid_param(i).SetFileName(j+2, name_file);
              }
          }
        
        for (int i = 0; i < var_output.output_mesh_param.GetM(); i++)
          {
            var_output.output_grid_param(i).ResizeNbFile(2+dim_N);
            for (int j = 0; j < dim_N; j++)
              {
                string name_file = var_output.output_mesh_param(i).GetFileName(1);
                name_file = GetBaseString(name_file) + "_V" 
                  + to_str(j) + "_." + GetExtension(name_file);
                
                var_output.output_mesh_param(i).SetFileName(j+2, name_file);
              }
          }
      }
    
    // number of dofs 
    int nodl = this->GetNumberOfUnknowns();
    
    check_prod_stiffness_matrix = false;
    if (check_prod_stiffness_matrix && !this->IsImplicitScheme())
      {
        /*var_problem.WriteMesh("test.mesh");
        {          
          DistributedMatrix<Real_wp, General, ArrayRowSparse> A(nodl, nodl);
          GlobalGenericMatrix<Real_wp> nat_mat;
          nat_mat.SetCoefMass(1.0); nat_mat.SetCoefStiffness(0.0); nat_mat.SetCoefDamping(0.0);
          var_computation.AddMatrixWithBC(A, nat_mat);
          A.WriteText("Mh.dat");

          var_boundary.SetCoefficientDirichlet(0.0);
          A.Clear(); A.Reallocate(nodl, nodl);
          nat_mat.SetCoefMass(0.0); nat_mat.SetCoefStiffness(-1.0); nat_mat.SetCoefDamping(-1.0);
          var_computation.AddMatrixWithBC(A, nat_mat);
          A.WriteText("Kh.dat");
        }
        */
	Matrix<Real_wp, General, ArrayRowSparse> Ah_full(nodl, nodl);
	VectReal_wp Ones(nodl), Ah_Ones(nodl);
	
        int size_U = this->GetNbScalarUnknowns();
        int size_V = nodl - size_U; DISP(size_U); DISP(size_V);
            
	if (nb_proc == 1)
	  {
            //DISP(Glob_mat_Kh.nature_matrix.GetCoefMass());
            //DISP(Glob_mat_Kh.nature_matrix.GetCoefStiffness());
            if (Glob_mat_Kh->mat_iterative_unsym.GetM() > 0)
              Glob_mat_Kh->mat_iterative_unsym.WriteText("Ah_direct.dat");
            else if (Glob_mat_Kh->mat_iterative_sym.GetM() > 0)
              Glob_mat_Kh->mat_iterative_sym.WriteText("Ah_direct.dat");
            else
              {
                // computing the finite element matrix by computing
                // the matrix vector product with canonical vectors e_i
                // this procedure can be very slow
                if ( !var_problem.FirstOrderFormulation()) 
                  {
                    for (int i = 0; i < nodl; i++)
                      {
                        Ones.Zero(); Ones(i) = Real_wp(1);
                        cout << "colonne " << i << endl;
                        Glob_mat_Kh->MltVector(Ones, Ah_Ones);
                        for (int j = 0; j < nodl; j++)
                          if (Ah_Ones(j) != Real_wp(0))
                            Ah_full.AddInteraction(j, i, Ah_Ones(j));
                      }
                    Ah_full.WriteText("Ah_full.dat");	
                  }
              }
	    
            // Glob_mat_Dh.DhPlusdtSh.invDiagonal.Write("Dh.dat");
            Ah_full.Clear();
            //DISP(size_U); DISP(size_V);
            Ah_full.Reallocate(nodl, nodl);
            /* Ah_Ones.Reallocate(size_U);
            for (int i = 0; i < size_U; i++)
              {
                Ah_Ones.Zero(); Ah_Ones(i) = Real_wp(1);	    
                this->SolveOperatorDh(Ah_Ones);
                for (int j = 0; j < size_U; j++)
                  if (Ah_Ones(j) != Real_wp(0))
                    Ah_full.AddInteraction(j, i, Ah_Ones(j));
              }
            
            Ah_Ones.Reallocate(size_V);
            for (int i = 0; i < size_V; i++)
              {
                Ah_Ones.Zero(); Ah_Ones(i) = Real_wp(1);	    
                this->SolveOperatorBh(Ah_Ones);
                for (int j = 0; j < size_V; j++)
                  if (Ah_Ones(j) != Real_wp(0))
                    Ah_full.AddInteraction(j+size_U, i+size_U, Ah_Ones(j));
              }
	    */
            
            for (int i = 0; i < nodl; i++)
              {
                Ah_Ones.Zero(); Ah_Ones(i) = Real_wp(1);	    
                this->SolveMassMatrix(Ah_Ones);
                for (int j = 0; j < nodl; j++)
                  if (Ah_Ones(j) != Real_wp(0))
                    Ah_full.AddInteraction(j, i, Ah_Ones(j));
              }
	    
	    Ah_full.WriteText("invMh_full.dat");
            
	    /* Ah_full.Clear();
	    Ah_full.Reallocate(nodl, nodl);
            Ah_Ones.Reallocate(size_U);
            for (int i = 0; i < size_U; i++)
              {
                Ones.Zero(); Ones(i) = Real_wp(1);	    
                this->ApplyOperatorSh(Real_wp(1), Real_wp(0), Ones, Real_wp(0), Ah_Ones);
                for (int j = 0; j < size_U; j++)
                  if (Ah_Ones(j) != Real_wp(0))
                    Ah_full.AddInteraction(j, i, Ah_Ones(j));
              }
	    
	      Ah_full.WriteText("Sh_full.dat"); */
          }
        
        Ah_Ones.Reallocate(nodl);
	Ah_Ones.Fill(0);
	Ah_full.Clear(); Ah_full.Reallocate(nodl, nodl);
        
        // int size_V  = TypeEquation::nb_unknowns_vec*Nvol;
        if ((nb_proc > 1) && (var_problem.FirstOrderFormulation()))
          {
            const IVect& GlobDofNumber = var_comm.GetGlobalDofNumber();
            string file_name = string("Ah_second_P") + to_str(rank_proc) + ".dat";
            ofstream file_out(file_name.data()); file_out.precision(16);
            for (int i = 0; i < var_comm.GetNbGlobalDof(); i++)
	      {                
                Ones.Zero(); Ah_Ones.Zero();
                for (int j = 0; j < nodl; j++)
                  if (GlobDofNumber(j) == i)
                    Ones(j) = Real_wp(1);

                cout<<"Colonne " << i << endl;
                this->EvaluateDerivativeFunction(Real_wp(0), 0, Ones, Ah_Ones, false, false);
		
                //for (int j = size_U; j < nodl; j++)
		//Ah_Ones(j) *= -1;
                
                for (int j = 0; j < nodl; j++)
                  if (abs(Ah_Ones(j)) > 1e-15)
                    file_out << GlobDofNumber(j)+1 << " " << i+1 << " " << Ah_Ones(j) << '\n';
              }
            
            file_out << var_comm.GetNbGlobalDof() << " " 
                     << var_comm.GetNbGlobalDof() << " 0.0 " << endl;
            
            file_out.close();

	    cout << "coucou" << endl;
	    file_name = string("invMh_P") + to_str(rank_proc) + ".dat";
            file_out.open(file_name.data()); file_out.precision(16);
	    VectReal_wp Uh, Vh, ProdUh, ProdVh;
	    Uh.SetData(size_U, const_cast<Real_wp*>(&Ones(0)));
	    Vh.SetData(size_U, const_cast<Real_wp*>(&Ones(size_U)));
	    ProdUh.SetData(size_V, const_cast<Real_wp*>(&Ah_Ones(0)));
	    ProdVh.SetData(size_V, const_cast<Real_wp*>(&Ah_Ones(size_U)));
            for (int i = 0; i < var_comm.GetNbGlobalDof(); i++)
	      {
		DISP(i);
		Ah_Ones.Zero();
                for (int j = 0; j < nodl; j++)
                  if (GlobDofNumber(j) == i)
                    Ah_Ones(j) = Real_wp(1);
		
                if (i < size_U)
		  this->SolveOperatorDh(Ah_Ones);
		else
		  this->SolveOperatorBh(Ah_Ones);

		for (int j = 0; j < nodl; j++)
                  if (abs(Ah_Ones(j)) > 1e-15)
                    file_out << GlobDofNumber(j)+1 << " " << i+1 << " " << Ah_Ones(j) << '\n';
              }

	    file_out.close();
	    Uh.Nullify(); Vh.Nullify();
            ProdUh.Nullify(); ProdVh.Nullify();
#ifdef SELDON_WITH_MPI
	    MPI_Barrier(var_comm.comm_group_mode);
#endif
          }
        else if (var_problem.FirstOrderFormulation())
          {
            for (int i = 0; i < nodl; i++)
              {
                Ones.Zero(); Ones(i) = Real_wp(1);
                cout<<"Colonne "<<i<<endl;
		this->EvaluateDerivativeFunction(Real_wp(0), 0, Ones, Ah_Ones, false, false);
                
                //for (int j = size_U; j < this->nodl; j++)
		//Ah_Ones(j) *= -1;
                
                for (int j = 0; j < nodl; j++)
                  if (abs(Ah_Ones(j)) > 1e-15)
                    {
                      //DISP(j); DISP(i); DISP(Ah_Ones(j));
                      Ah_full.AddInteraction(j, i, Ah_Ones(j));
                    }
              }            
            
            Ah_full.WriteText("Ah_second.dat");
          }
        
        exit(0);
      }
    
    // evaluation of CFL if required
    if (time_step_to_be_computed)
      {
        Real_wp cfl = EvaluateCFL(false);
	if (rank_proc == 0)
          {
            cout << "Maximal time step = " << cfl << endl;
            if (!this->FirstOrderScheme())
              cout << "This time step is only valid when there is "
                   << "no PML and/or absorbing boundary conditions" << endl;
          }
        
        this->deltat = 0.999*cfl;
      }
  }
  

  //! returns the size of the vectorial solution
  int VarInstationary_Base::GetNbVectorialUnknowns() const
  {
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return var_problem.GetNbMeshDof()*var_problem.nb_unknowns_vec;
    else if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      return var_comm.GetNbMainUnknownDof()*var_problem.nb_unknowns_vec;
    else
      {
	if (var_comm.GetSizeOffsetDofV() <= 0)
	  return 0;
	
	return var_comm.GetOffsetDofV(var_computation.GetNbElt());;
      }
  }
  

  //! returns the size of the vectorial solution in PML layers
  int VarInstationary_Base::GetNbVectorialUnknownsPML() const
  {
    int nb_vec = var_problem.nb_unknowns_vec;
    int nb_dofV = 0;
    
    if (var_problem.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
      {
	for (int i = 0; i < var_computation.GetNbElt(); i++)
	  if (var_problem.ElementInsidePML(i))
	    {
	      int nb_pts_quad = var_problem.GetNbPointsQuadratureInside(i);
	      nb_dofV += nb_pts_quad;
	    }
      }
    else
      {
	for (int i = 0; i < var_computation.GetNbElt(); i++)
	  if (var_problem.ElementInsidePML(i))
	    {
	      int nb_pts_quad = var_problem.GetNbLocalDof(i);
	      nb_dofV += nb_pts_quad;
	    }
      }
    
    nb_dofV *= nb_vec;
    return nb_dofV;
  }
  
  
  //! computation of mass matrix
  void VarInstationary_Base::ComputeMassMatrix()
  {
    if (Glob_mat_Dh == NULL)
      Glob_mat_Dh = GetNewMassMatrix();
    
    this->Glob_mat_Dh->Init();    
  }


  //! computation of stiffness matrix K (and RhS, RhV if required)
  void VarInstationary_Base::ComputeStiffnessMatrix()
  {
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(0.0); nat_mat.SetCoefStiffness(1.0);
    nat_mat.SetCoefDamping(0.0);
    
    Glob_mat_Kh = var_computation.GetNewIterativeMatrix(Real_wp(0));
    Glob_mat_Kh->SetCoefficientMatrix(nat_mat);
    Glob_mat_Kh->IgnoreDirichletDof();

    var_computation.AddMatrixWithBC(*this->Glob_mat_Kh, nat_mat);

    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	int nb_dof_uv = var_comm.GetNbMainUnknownDof()*var_problem.nb_unknowns;
	int nb_dof_L = var_problem.GetNbMeshDof()*var_problem.nb_unknowns_hdg;
	
	DistributedMatrix<Real_wp, General, ArrayRowSparse>& K = Glob_mat_Kh->mat_iterative_unsym;

	// matrices are not extracted if matrix-free storage is selected
	if (K.GetM() <= 0)
	  return;

	// case where matrices are stored
	DistributedMatrix<Real_wp, General, ArrayRowSparse>& A = Glob_mat_RhV;
	DistributedMatrix<Real_wp, General, ArrayRowSparse>& Al = Glob_mat_Ch_Lambda;
	DistributedMatrix<Real_wp, General, ArrayRowSparse>& B = Glob_mat_RhS;
	
	A.Reallocate(nb_dof_uv, nb_dof_uv);
	Al.Reallocate(nb_dof_uv, nb_dof_L);
	for (int i = 0; i < nb_dof_uv; i++)
	  {
	    int size_row = K.GetRowSize(nb_dof_L + i);
	    int nb_val = 0, nb_val2 = 0;
	    for (int j = 0; j < size_row; j++)
	      {
		if (K.Index(nb_dof_L+i, j) < nb_dof_L)
		  nb_val++;
		else
		  nb_val2++;
	      }

	    Al.ReallocateRow(i, nb_val);
	    A.ReallocateRow(i, nb_val2);
	    nb_val = 0; nb_val2 = 0;
	    for (int j = 0; j < size_row; j++)
	      {
		if (K.Index(nb_dof_L+i, j) < nb_dof_L)
		  {
		    Al.Index(i, nb_val) = K.Index(nb_dof_L + i, j);
		    Al.Value(i, nb_val) = K.Value(nb_dof_L + i, j);
		    nb_val++;
		  }
		else
		  {
		    A.Index(i, nb_val2) = K.Index(nb_dof_L + i, j) - nb_dof_L;
		    A.Value(i, nb_val2) = K.Value(nb_dof_L + i, j);
		    nb_val2++;
		  }
	      }
	    
	    K.ClearRow(nb_dof_L + i);
	  }
	
	B.Reallocate(nb_dof_L, nb_dof_uv);
	//C.Reallocate(nb_dof_L);
	//C.Zero();
	for (int i = 0; i < nb_dof_L; i++)
	  {
	    int size_row = K.GetRowSize(i);
	    int nb_val = 0;
	    for (int j = 0; j < size_row; j++)
	      if (K.Index(i, j) >= nb_dof_L)
		nb_val++;
	    
	    B.ReallocateRow(i, nb_val);
	    nb_val = 0;
	    for (int j = 0; j < size_row; j++)
	      {
		if (K.Index(i, j) >= nb_dof_L)
		  {
		    B.Index(i, nb_val) = K.Index(i, j) - nb_dof_L;
		    B.Value(i, nb_val) = K.Value(i, j);
		    nb_val++;
		  }
		else
		  {
		    /* if (K.Index(i, j) != i)
		      {
			cout << "non-diagonal case " << endl;
			abort();
		      }
		    else
		    C(i) = K.Value(i, j); */
		  }
	      }
	    
	    K.ClearRow(i);
	  }
    
	K.Clear();
	
	return;
      }
    
    if ( (this->SplitSystem())
         && (var_problem.FirstOrderFormulation()))
      {
        if (this->Glob_mat_Kh->mat_iterative_sym.GetM() > 0)
          {
            //this->ExtractStiffnessMatrix(this->Glob_mat_Kh.mat_iterative_sym_real);
            cout << "Symmetric case not handled" << endl;
            abort();
          }
        else if (this->Glob_mat_Kh->mat_iterative_unsym.GetM() > 0)
          this->ExtractStiffnessMatrix(this->Glob_mat_Kh->mat_iterative_unsym);
      }
  }
  

  //! extracting scalar and vectorial stiffness matrices RhS and RhV
  template<class MatrixSparse>
  void VarInstationary_Base::ExtractStiffnessMatrix(MatrixSparse& K)
  {
    int nb_dof_u, nb_dof_v;
    nb_dof_u = this->GetNbScalarUnknowns();
    nb_dof_v = this->GetNbVectorialUnknowns();
    Glob_mat_RhS.Reallocate(nb_dof_u, nb_dof_v);
    Glob_mat_RhV.Reallocate(nb_dof_v, nb_dof_u);
    for (int i = 0; i < nb_dof_u; i++)
      {
        int size_row = K.GetRowSize(i);
        int nb_val = 0;
        for (int j = 0; j < size_row; j++)
          if (K.Index(i, j) >= nb_dof_u)
            nb_val++;
        
        Glob_mat_RhS.ReallocateRow(i, nb_val);
        nb_val = 0;
        for (int j = 0; j < size_row; j++)
          if (K.Index(i, j) >= nb_dof_u)
            {
              Glob_mat_RhS.Index(i, nb_val) = K.Index(i, j) - nb_dof_u;
              Glob_mat_RhS.Value(i, nb_val) = -K.Value(i, j);
              nb_val++;
            }
        
        K.ClearRow(i);
      }
    
    for (int i = 0; i < nb_dof_v; i++)
      {
        int size_row = K.GetRowSize(i+nb_dof_u);
        int nb_val = 0;
        for (int j = 0; j < size_row; j++)
          if (K.Index(i+nb_dof_u, j) < nb_dof_u)
            nb_val++;
        
        Glob_mat_RhV.ReallocateRow(i, nb_val);
        nb_val = 0;
        for (int j = 0; j < size_row; j++)
          if (K.Index(i+nb_dof_u, j) < nb_dof_u)
            {
              Glob_mat_RhV.Index(i, nb_val) = K.Index(i+nb_dof_u, j);
              Glob_mat_RhV.Value(i, nb_val) = -K.Value(i+nb_dof_u, j);
              nb_val++;
            }
        
        K.ClearRow(i+nb_dof_u);
      }

    K.Clear();
  }


  //! applying second-order time schemes
  void VarInstationary_Base::ApplySecondOrderScheme()
  {  
    // checking that no PML are present
    if (var_boundary.GetNbEltPML())
      {
        cout << "PML layers are implemented only for first-order formulations" << endl;
        abort();
      }
    
    if (var_problem.FirstOrderFormulationDG())
      {
	// no second order scheme for LDG formulation
	cout << "LDG formulation requires first-order schemes" << endl;
	abort();
      }
    else
      {
	int N = var_problem.GetNbDof();
	VectReal_wp Xh(2*N);
	
	// initial conditions for U(0) and dU/dt(0)
	Xh.Fill(0);
	this->SetInitialVector(this->initial_time, Xh);
	VectReal_wp U0(N), dU0_dt(N);
	for (int i = 0; i < N; i++)
	  {
	    U0(i) = Xh(i);
	    dU0_dt(i) = Xh(i+N);
	  }
    
	Xh.Clear();
	
	RunSecondOrderScheme(this->initial_time, this->final_time, this->deltat,
			     U0, dU0_dt, *this, var_time_scheme.GetTimeScheme());
      }
  }

  
  //! Running time iterations for "generic" time schemes
  /*!
    Time schemes like Runge-Kutta or Tal-Ezer, solve the evolution system
    dU/dt = A*U + F(t)   where A is generally a linear operator 
    all the time evolution problems solved by Montjoie 
    (maxwell, elastodynamic, aeroacoustic ...) can be written
    in this form, so these time schemes are "generic" to any equation
    solved by Montjoie
    
    the operation U -> A*U + F(t) is realized by the method EvaluateDerivativeFunction(t, n, X, Y)
    This method may be overloaded depending the solved equation
  */
  void VarInstationary_Base::RunTimeIterations()
  {    
    if (!this->FirstOrderScheme())
      {
        this->ApplySecondOrderScheme();
        return;
      }
    
    int N = this->GetNumberOfUnknowns();
    int Nscal = this->GetNbScalarUnknowns();
    int Nvec = this->GetNbVectorialUnknowns();
    
    VectReal_wp Xh(N);
    
    // initial conditions
    Xh.Zero();
    if (small_random_initial)
      {
        Xh.FillRand();
        Mlt(Real_wp(1e-270), Xh);
      }
    
    int num_iter0 = 0;
    if (this->load_reprise)
      {
        this->LoadDataReprise(Xh);
        num_iter0 = this->load_iter_reprise;
      }
    else
      this->SetInitialVector(this->initial_time, Xh);
    
    if (this->var_time_scheme.GetTimeSchemeType() == TimeSchemeEnum::LOCAL_IMPERIALE_SCHEME)
      {
	// running iterations
	//RunTimeScheme(this->initial_time, this->final_time, this->deltat,
	//this->local_deltat, Xh, var_leaf, this->var_time_scheme. local_time_scheme);
	abort();
      }
    else if (this->var_time_scheme.GetTimeSchemeType() == TimeSchemeEnum::LOCAL_PIPERNO_SCHEME)
      {
        VectReal_wp Uh(Nscal), Vh(Nvec);
        for (int i = 0; i < Nscal; i++)
          Uh(i) = Xh(i);
        
        for (int i = 0; i < Nvec; i++)
          Vh(i) = Xh(Nscal+i);
        
        Xh.Clear();
        
	abort();
	// running iterations
	//RunTimeScheme(this->initial_time, this->final_time, this->deltat,
	//this->local_deltat, Uh, Vh, var_leaf,
	//this->var_time_scheme.local_piperno_scheme);
      }
    else if (var_time_scheme.StaggeredScheme())
      {
        VectReal_wp Uh(Nscal), Vh(Nvec);
	for (int i = 0; i < Nscal; i++)
	  Uh(i) = Xh(i);
	
	for (int i = 0; i < Nvec; i++)
	  Vh(i) = Xh(i+Nscal);
	
	Xh.Clear();
	
	// run iterations
	RunFirstOrderScheme(this->initial_time, this->final_time,
                            this->deltat, Uh, Vh, *this, var_time_scheme.GetTimeScheme());
      }
    else
      {
        // run iterations
	RunTimeScheme(this->initial_time, this->final_time,
                      this->deltat, Xh, *this, var_time_scheme.GetTimeScheme(), num_iter0);
      }
    
    /*number of terms to use for Talezer scheme
      IVect nb_modes(16); nb_modes.Fill(0);
      for (int power = 2; power < 16; power++)
      {
      Real_wp threshold_rhs = pow(10.0,-Real_wp(power))*this->deltat/final_time;
      Real_wp kr = talezer_scheme.Rmax*this->deltat;
      Real_wp lb = GetLambertW0(2.0/(3.0*pi_wp*threshold_rhs*threshold_rhs)); // DISP(lb);
      nb_modes(power) = int(kr + 0.5*pow(1.5*lb, 2.0/3.0)*pow(kr,1.0/3.0))+1;
      }
    */
  }
  
  
  //! generic run all
  /*! 
    \param[in] input_file name of the input file used to launch the simulation
    Launching the complete simulation
    - Construction of the mesh and finite element
    - Computation of the source terms
    - Computations of matrices if necessary (implicit scheme)
    - Precomputations for mass matrices and boundary conditions
    - Finally, time iterations 
  */
  void VarInstationary_Base::RunAll(const string& input_file, const string& name_element,
				    const string& name_equation, int num)
  {

    var_problem.SetTypeEquation(name_equation);
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    if (rank_proc == 0)
      cout << rank_proc << " Beginning of simulation" << endl;
    
    var_problem.InitIndices(PhysicalConstant::nb_max_indices);
        
    // The input file is read, see file Reading_InputFiles.cxx
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(input_file, lines_data_file, var_comm.comm_group_mode);
#else
    ReadLinesFile(input_file, lines_data_file);
#endif

    ReadInputFile(lines_data_file, *this);

    glob_solver = var_computation.GetNewLinearSolver();
    ReadInputFile(lines_data_file, *glob_solver);
    
    int print_level = var_problem.print_level;
    if ((rank_proc == 0)||(print_level >= 10))
      if (print_level >= 0)
	{
	  cout << "Proc " << rank_proc;
	  cout << " out of " << nb_proc << " input file read" << endl;
	}

    if (nb_iter_per_display == -1)
      {
        if (print_level <= 0)
          nb_iter_per_display = 10000;
        if (print_level == 1)
          nb_iter_per_display = 1000;
        else if (print_level == 2)
          nb_iter_per_display = 200;
        else if (print_level == 3)
          nb_iter_per_display = 100;
        else if (print_level == 4)
          nb_iter_per_display = 20;
        else if (print_level == 5)
          nb_iter_per_display = 10;
        else if (print_level == 6)
          nb_iter_per_display = 5;
        else if (print_level == 7)
          nb_iter_per_display = 2;
        else
          nb_iter_per_display = 1;
      }
    
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (var_output.DOSSIER_output.size() == 0)
      {
	var_output.DOSSIER_output = "[STIFFOUT]/"; 
	EcritDossier(var_output.DOSSIER_output, input_file, num);
      }
    
    if (this->FirstOrderScheme())
      var_problem.SetFirstOrderFormulation(true);
    
    MontjoieTimer& var_chrono = var_problem.var_chrono;

#ifdef SELDON_WITH_MPI
    glob_chrono.SetCommunicator(var_comm.comm_group_mode);
    var_chrono.SetCommunicator(var_comm.comm_group_mode);
#endif
    
    var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
    var_chrono.Start("GlobalComputation");

    var_chrono.SetMessage("MeshGeneration", "construct and number the mesh");
    var_chrono.Start("MeshGeneration");
    
    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Constructing the mesh... " << endl;

    // Mesh is read, finite element constructed
    var_problem.SetSameNumberPeriodicDofs();
    
    // implementation of this method is done in VarProblem.cxx
    var_problem.ComputeMeshAndFiniteElement(name_element);

    //MPI_Barrier(var_comm.comm_group_mode);
    if ((print_level >= 6) && (rank_proc == 0))
      cout << rank_proc << "Mesh and finite element constructed " << endl;
    
    // Dirichlet condition (finding dofs associated to Dirichlet condition)
    // method implemented in TreatBoundaryConditions.cxx
    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Dirichlet condition ... " << endl;
    
    var_problem.PerformOtherInitializations();
    //MPI_Barrier(var_comm.comm_group_mode);

    var_chrono.Stop("MeshGeneration");

    var_chrono.SetMessage("JacobianMatricesComputation", "compute jacobian matrices");
    var_chrono.Start("JacobianMatricesComputation");

    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Computation of geometric quantities " << endl;

    // computation of DF_i^{-1} and other variables (J_i, normales)
    // method is implemented in Harmonic_MassMatrix.cxx
    // this->var_harmonic.ompute = true;
    var_computation.ComputeMassMatrix();
    var_boundary.ComputeQuasiPeriodicPhase();
    //MPI_Barrier(var_comm.comm_group_mode);    
    if ((rank_proc == 0)||(print_level >= 10))
      if (print_level >= 5)
	cout << rank_proc << " Geometric quantities OK " << endl;
    
    var_chrono.Stop("JacobianMatricesComputation");

    var_chrono.SetMessage("RightHandSideComputation", "compute the right hand side");
    var_chrono.Start("RightHandSideComputation");
    
    // computation of the right-hand side (source terms)
    // this method is implemented in this file VarInstationary.cxx
    this->ComputeRightHandSide();
    
    //MPI_Barrier(var_comm.comm_group_mode);    
    if ((rank_proc == 0)||(print_level >= 10))
      if (print_level >= 5)
	cout << rank_proc << " Right hand side computed " << endl;

    var_chrono.Stop("RightHandSideComputation");

    // for local time schemes, we compute deltat on each element
    if (this->var_time_scheme.LocalTimeStepping())
      this->ComputeCoarseFineRegion();
    
    var_chrono.SetMessage("ComputationMatrices", "compute mass and stiffness matrices");
    var_chrono.Start("ComputationMatrices");

    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Computing stiffness matrix (unsteady)... " << endl;

    // computation of stiffness matrix
    this->ComputeStiffnessMatrix();
    
    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Mass matrix for unsteady problems... " << endl;

    this->ComputeMassMatrix();

    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Mass matrix computed " << endl;
    
    var_chrono.Stop("ComputationMatrices");

#ifdef SELDON_WITH_MPI
    MPI_Barrier(var_comm.comm_group_mode);
#endif
    
    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Initialization iterations ... " << endl;
    
    // last initialization before time iterations
    this->InitTimeIterations();
    
    // mesh is written on the disk
    if (nb_proc == 1)
      var_problem.WriteMesh(var_output.DOSSIER_output+string("test.mesh"));
    
#ifdef SELDON_WITH_MPI
    MPI_Barrier(var_comm.comm_group_mode);
#endif
    
    if ((print_level >= 5) && (rank_proc == 0))
      cout << rank_proc << " Launching iterations ... " << endl;
    
    var_chrono.SetMessage("TimeIteration", "complete time iterations");
    var_chrono.Start("TimeIteration");
    
    // Running time iterations
    this->RunTimeIterations();
    
    var_chrono.Stop("TimeIteration");
    var_chrono.Stop("GlobalComputation");
    
    if (print_level >= 1)
      var_chrono.DisplayAll();

#ifdef SELDON_WITH_MPI
    MPI_Barrier(var_comm.comm_group_mode);
#endif

  }

  
  //! evaluation of Y'' = G(tn, Y, Yp) (second-order formulation) 
  void VarInstationary_Base::
  EvaluateFunctionS(const Real_wp& tn, const VectReal_wp& Y, const VectReal_wp& Yp,
		    VectReal_wp& ProdY, bool invert_mass, bool source)
  {
    VectReal_wp Ydir(var_boundary.GetNbDirichletDof()), Ydir_p;
    VectReal_wp& Y_ = const_cast<VectReal_wp&>(Y);
    VectReal_wp& Yp_ = const_cast<VectReal_wp&>(Yp);

    // Dirichlet condition is set on Y
    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      Ydir(i) = Y(var_boundary.GetDirichletDofNumber(i));
    
    if (source)
      this->SetDirichletCondition(tn, 0, Y_);
    else
      var_boundary.ImposeNullDirichletCondition(Y_);
    
    // and on y'
    if (Yp.GetM() > 0)
      {
	Ydir_p.Reallocate(var_boundary.GetNbDirichletDof());
	for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	  Ydir_p(i) = Yp(var_boundary.GetDirichletDofNumber(i));
    
	if (source)
	  this->SetDirichletCondition(tn, 1, Yp_);
	else
	  var_boundary.ImposeNullDirichletCondition(Yp_);
      }

    ProdY.Zero();
    if (source)
      this->AddPrimitiveSourceAtTime(Real_wp(1), tn, 1, ProdY); // DISP(Fh.GetM());
    
    Glob_mat_Kh->MltAddVector(Real_wp(-1), Y, Real_wp(1), ProdY, false);
    if (Yp.GetM() > 0)
      this->ApplyOperatorSh(Real_wp(-1), tn, Yp, Real_wp(1), ProdY);

    // on Dirichlet dofs, we put 0
    var_boundary.ImposeNullDirichletCondition(ProdY);

    if (invert_mass)
      this->SolveOperatorDh(ProdY);

    // stored values on dirichlet dofs are recovered
    if (Yp.GetM() > 0)
      for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	Yp_(var_boundary.GetDirichletDofNumber(i)) = Ydir_p(i);

    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      Y_(var_boundary.GetDirichletDofNumber(i)) = Ydir(i);

    // on Dirichlet dofs, we put 0 (again)
    var_boundary.ImposeNullDirichletCondition(ProdY);
    
    // alternative solution : y''_i = f'''_i on Dirichlet dofs
    // this->SetDirichletCondition(tn, 2, ProdY);
  }
  
  
  //! compute Y = G(tn,X) if we write the evolution problem as dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
    if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
    \param[in] invert_mass if true we invert G(tn, X) by mass matrix (for explicit schemes)
  */
  void VarInstationary_Base::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X,
			     VectReal_wp& Y, bool invert_mass, bool source)
  {
    Y.Fill(0);
    if (source)
      this->AddPrimitiveSourceAtTime(Real_wp(1), tn, nb_deriv, Y); // DISP(Fh.GetM());
    
    //DISP(Norm2(X)); DISP(Norm2(Y));
    Glob_mat_Kh->MltAddVector(Real_wp(-1), X, Real_wp(1), Y);
    
    //DISP(Norm2(Y));
    if (invert_mass)
      Glob_mat_Dh->SolveMassMatrix(Y);
    
    //DISP(Norm2(Y));
  }
  

  //! computation and factorisation of matrix Dh + beta Sh + gamma Kh
  void VarInstationary_Base::
  FactorizeOperatorDhPlusGammaKh(const Real_wp& alpha, const Real_wp& beta, const Real_wp& gamma)
  {
    GlobalGenericMatrix<Real_wp> nat_mat(alpha, beta, gamma);
    this->glob_solver->SetPrintLevel(var_problem.print_level);

    //this->var_computation.storage_matrix_asked = true;
    //this->var_computation.file_name_matrix_stored = "DhPlusKh.dat";

    var_boundary.SetCoefficientDirichlet(Real_wp(1));

    bool compute_prec = true;
    if (same_preconditioning && first_preconditioning_computed)
      compute_prec = false;
    
    this->glob_solver->PerformFactorizationStep(nat_mat, compute_prec);

    if (!compute_prec)
      this->glob_solver->SetRealPreconditioning(ptr_real_precond);
    else
      {
	ptr_real_precond = &this->glob_solver->GetPreconditioning(Real_wp(0));
	first_preconditioning_computed = true;
      }
    
    //this->var_computation.storage_matrix_asked = false;

    /*
    DISP(nat_mat.GetCoefStiffness());     DISP(nat_mat.GetCoefMass());

    int N = this->GetNumberOfUnknowns();
    VectReal_wp x(N), y(N);
    x.FillRand(); x *= 1e-9;

    x.WriteText("x.dat");
    
    y = x;
    this->SolveOperatorDhPlusGammaKh(Real_wp(0), x, y);
    
    y.WriteText("y.dat");
    int test_input; cout << "waiting" << endl; cin >> test_input;
    */
    /*int Nvol_glob = var_comm.GetNbGlobalDof();
    if (var_comm.GetNbProcPerMode() == 1)
      Nvol_glob = var_problem.GetNbDof();

    srand(0); DISP(var_problem.GetNbScalarDof());

    VectReal_wp xglob(Nvol_glob), x(var_problem.GetNbDof()), y(var_problem.GetNbDof());
    xglob.FillRand(); xglob *= 1e-9;
    const IVect& global_row = var_comm.GetGlobalDofNumber();
    if (var_comm.GetNbProcPerMode() > 1)
      {
	for (int i = 0; i < global_row.GetM(); i++)
	  x(i) = xglob(global_row(i));

	for (int i = 0; i < var_comm.GetNbOverlappedDof(); i++)
	  x(var_comm.GetOverlappedDofNumber(i)) = 0;
      }
    else
      x = xglob;
    
    int rank_proc = var_comm.GetRankProcMode();
    global_row.WriteText("numGlob"+to_str(rank_proc) + ".dat");
    x.Write("x_P" + to_str(rank_proc) + ".dat");

    y = x;
    this->glob_solver->ComputeSolution(y, nat_mat);
    
    y.Write("y_P" + to_str(rank_proc) + ".dat");
    abort(); */
    
    if (var_problem.print_level < 8)
      this->glob_solver->SetPrintLevel(0);        
  } 


  //! computation and factorisation of matrix Dh + beta Sh + gamma Kh
  void VarInstationary_Base::
  FactorizeOperatorReal(const VectReal_wp& alpha,
			const VectReal_wp& beta, const VectReal_wp& gamma)
  {
    for (int i = 0; i < vec_solver.GetM(); i++)
      if (vec_solver(i) != NULL)
	delete vec_solver(i);

    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    vec_solver.Reallocate(alpha.GetM());
    for (int n = 0; n < alpha.GetM(); n++)
      {
	vec_solver(n) = var_computation.GetNewLinearSolver();
	vec_solver(n)->CopyParameter(*this->glob_solver);
	
	GlobalGenericMatrix<Real_wp> nat_mat(alpha(n), beta(n), gamma(n));
	vec_solver(n)->SetPrintLevel(var_problem.print_level);

	bool compute_prec = true;
	if (same_preconditioning && first_preconditioning_computed)
	  compute_prec = false;
	
	vec_solver(n)->PerformFactorizationStep(nat_mat, compute_prec);
	
	if (!compute_prec)
	  this->vec_solver(n)->SetRealPreconditioning(ptr_real_precond);
	else
	  {
	    ptr_real_precond = &this->vec_solver(n)->GetPreconditioning(Real_wp(0));
	    first_preconditioning_computed = true;
	  }
	
	if (var_problem.print_level < 8)
	  vec_solver(n)->SetPrintLevel(0);
      }
  }


  //! computation and factorisation of matrix Dh + beta Sh + gamma Kh
  void VarInstationary_Base::
  FactorizeOperatorComplex(const VectComplex_wp& alpha,
			   const VectComplex_wp& beta, const VectComplex_wp& gamma)
  {
    for (int i = 0; i < vec_solver.GetM(); i++)
      if (vec_solver(i) != NULL)
	delete vec_solver(i);

    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    vec_solver.Reallocate(alpha.GetM());
    for (int n = 0; n < alpha.GetM(); n++)
      {
	vec_solver(n) = var_computation.GetNewLinearSolver();
	vec_solver(n)->CopyParameter(*this->glob_solver);
	
	GlobalGenericMatrix<Complex_wp> nat_mat(alpha(n), beta(n), gamma(n));
	vec_solver(n)->SetPrintLevel(var_problem.print_level);
	
	//this->var_computation.storage_matrix_asked = true;
	//this->var_computation.file_name_matrix_stored = "DhPlusKh.dat";

	bool compute_prec = true;
	if (same_preconditioning && first_preconditioning_computed)
	  compute_prec = false;
	
	vec_solver(n)->PerformFactorizationStep(nat_mat, compute_prec);

	if (!compute_prec)
	  this->vec_solver(n)->SetRealPreconditioning(ptr_real_precond);
	else
	  {
	    ptr_real_precond = &this->vec_solver(n)->GetPreconditioning(Real_wp(0));
	    first_preconditioning_computed = true;
	  }
	
	if (var_problem.print_level < 8)
	  vec_solver(n)->SetPrintLevel(0);
      }
  }
  

  //! solves the linear system  (Dh + beta Sh + gamma Kh) Y = B
  void VarInstationary_Base::
  SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& B, VectReal_wp& Y)
  {
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	VectReal_wp Yall(var_problem.GetNbDof());
	Yall.Zero();
	int offset_L = var_problem.GetNbMeshDof()*var_problem.nb_unknowns_hdg; 
	for (int i = 0; i < B.GetM(); i++)
	  Yall(offset_L + i) = B(i);
	
	glob_solver->ComputeSolution(Yall);

	for (int i = 0; i < B.GetM(); i++)
	  Y(i) = Yall(offset_L + i);

	return;
      }
    
    Y = B;
    glob_solver->ComputeSolution(Y);
  }

  
  //! solves the real system Dh + beta_n Sh + gamma_n Kh Y = B
  //! n is the number of the linear system
  void VarInstationary_Base::
  SolveOperatorReal(const Real_wp& t, const VectReal_wp& B, VectReal_wp& Y, int n)
  {
    if (n >= vec_solver.GetM())
      {
	cout << "Solver not allocated" << endl;
	abort();
      }

    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	VectReal_wp Yall(var_problem.GetNbDof());
	Yall.Zero();
	int offset_L = var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
	for (int i = 0; i < B.GetM(); i++)
	  Yall(offset_L + i) = B(i);
	
	vec_solver(n)->ComputeSolution(Yall);

	for (int i = 0; i < B.GetM(); i++)
	  Y(i) = Yall(offset_L + i);

	return;
      }

    Y = B;
    vec_solver(n)->ComputeSolution(Y);
  }  

  
  //! solves the complex system Dh + beta_n Sh + gamma_n Kh Y = B
  //! n is the number of the linear system
  void VarInstationary_Base::
  SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& B, VectComplex_wp& Y, int n)
  {
    if (n >= vec_solver.GetM())
      {
	cout << "Solver not allocated" << endl;
	abort();
      }

    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	VectComplex_wp Yall(var_problem.GetNbDof());
	Yall.Zero();
	int offset_L = var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
	for (int i = 0; i < B.GetM(); i++)
	  Yall(offset_L + i) = B(i);
	
	vec_solver(n)->ComputeSolution(Yall);

	for (int i = 0; i < B.GetM(); i++)
	  Y(i) = Yall(offset_L + i);
	
	return;
      }

    Y = B;
    vec_solver(n)->ComputeSolution(Y);
  }  
  
  
  //! solves mass matrix, X = M^{-1} X
  void VarInstationary_Base::SolveMassMatrix(VectReal_wp& X)
  {
    if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
      this->Assemble(X);

    Glob_mat_Dh->SolveMassMatrix(X);
  }


  //! Applies mass matrix, y = beta y + alpha M x
  void VarInstationary_Base::ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
					     const VectReal_wp& x, const Real_wp& beta, VectReal_wp& y)
  {
    Glob_mat_Dh->ApplyMassMatrix(alpha, t, x, beta, y);
  }
  

  //! compute Y = G(tn,X) if we write the evolution problem as dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
    if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
    \param[in] invert_mass if true we invert G(tn, X) by mass matrix (for explicit schemes)
  */
  void VarInstationary_Base::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha,
                             int level, const VectReal_wp& X,
			     VectReal_wp& Y, bool invert_mass, bool source)
  {
    Glob_mat_Kh->MltAddVector(-alpha, X, Real_wp(0), Y, false);
    if (source)
      this->AddPrimitiveSourceAtTime(alpha, tn, nb_deriv, Y);
        
    if (invert_mass)
      Glob_mat_Dh->SolveMassMatrix(Y);
  }
  

  //! returns mass and damping coefficient of element i (point j), unknown num
  void VarInstationary_Base
  ::GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp&, Real_wp&)
  {
    cout  << "GetMassDampingCoefficient not implemented for this equation" << endl;
    abort(); 
  }


  //! returns the tau parameter of HDG method
  Real_wp VarInstationary_Base::GetCoefficientTauHDG(int ref) const
  {
    cout  << "GetCoefficientTauHDG not implemented for this equation" << endl;
    abort();
    return Real_wp(0);
  }
  

  //! adds a part of operator Bh
  void VarInstationary_Base::AddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
                                                       const VectReal_wp& Uh, VectReal_wp& Vh,
                                                       const Real_wp&, const Real_wp&)
  {
    cout  << "AddElementStiffnessTensor not implemented for this equation" << endl;
    abort();
  }
  

  //! solves a part of operator Bh
  void VarInstationary_Base
  ::SolveElementStiffnessTensor(int i, int& offset, VectReal_wp& Vh, const Real_wp&)
  {
    cout  << "SolveElementStiffnessTensor not implemented for this equation" << endl;
    abort();
  }
    

  //! solves and adds a part of operator Bh
  void VarInstationary_Base
  ::SolveAddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
                                   const VectReal_wp& Uh, VectReal_wp& Vh, const Real_wp&)
  {
    cout  << "SolveAddElementStiffnessTensor not implemented for this equation" << endl;
    abort();
  }


  //! multiplies by the scalar part of stiffness matrix
  void VarInstationary_Base
  ::MltAddStiffnessScalar(const Real_wp&, int lvl, const VectReal_wp&,
                          const Real_wp&, VectReal_wp&)
  {
    cout  << "MltAddStiffnessScalar not implemented for this equation" << endl;
    abort(); 
  } 


  //! multiplies by the vectorial part of stiffness matrix
  void VarInstationary_Base
  ::MltAddStiffnessVectorial(const Real_wp&, int lvl, const VectReal_wp&,
                             const Real_wp&, VectReal_wp&)
  {
    cout  << "MltAddStiffnessVectorial not implemented for this equation" << endl;
    abort();
  }
  

  //! computes coarse and fine region
  void VarInstationary_Base::ComputeCoarseFineRegion()
  {
    try
      {
	LocalTimeScheme_Base& scheme
	  = dynamic_cast<LocalTimeScheme_Base& >(var_time_scheme.GetTimeScheme());
	
	IVect num_implicit;
	if (scheme.GetSplittingType() == LocalTimeScheme_Base::AUTO_SPLITTING)
	  {
	    VectReal_wp dt_elt;
	    this->ComputeLocalTimeStep(dt_elt);
	    
	    int ne = 0;
	    for (int i = 0; i < dt_elt.GetM(); i++)
	      if (dt_elt(i) < scheme.GetSplittingDt())
		ne++;
	    
	    num_implicit.Reallocate(ne);
	    ne = 0;
	    for (int i = 0; i < dt_elt.GetM(); i++)
	      if (dt_elt(i) < scheme.GetSplittingDt())
		num_implicit(ne++) = i;	
	  }
	else if (scheme.GetSplittingType() == LocalTimeScheme_Base::FILE_SPLITTING)
	  num_implicit.ReadText(scheme.GetNumberElementFileName());

	scheme.SetImplicitElement(num_implicit);
	this->DisplayCoarseFineRegion(num_implicit);
      }
    catch (const std::bad_cast&)
      {
	cout << "Time scheme not inherited from LocalTimeScheme_Base" << endl;
	abort();
      }
  }
  

  //! methods that computes L from (u, v)
  void VarInstationary_Base::GetExtrapolationLambdaHDG(int level, const VectReal_wp& Y, VectReal_wp& ProdL)
  {
    cout << "Only available for HDG formulation" << endl;
    abort();
  }


  //! methods that computes the (u,v) part of f(y)
  void VarInstationary_Base
  ::MltAddStiffnessHDG(const Real_wp& alpha, int lvl, int lvl2, const VectReal_wp& Y, const VectReal_wp& L,
		       const Real_wp& beta, VectReal_wp& ProdY)
  {
    cout << "Only available for HDG formulation" << endl;
    abort();
  }


  //! methods that updates ProdY with contribution of Lambda unknown only
  void VarInstationary_Base
  ::MltAddLambdaHDG(const Real_wp& alpha, int level, const VectReal_wp& L, VectReal_wp& ProdY)
  {
    cout << "Only available for HDG formulation" << endl;
    abort();
  }


  //! methods that updates ProdY with contribution of Lambda unknown only
  void VarInstationary_Base
  ::MltAddLambdaHDG(const Complex_wp& alpha, int level, const VectComplex_wp& L, VectComplex_wp& ProdY)
  {
    cout << "Only available for HDG formulation" << endl;
    abort();
  }

  
  /***********************
   * VarInstationary_Dim *
   ***********************/
  

  //! initialization of source
  template<class Dimension>
  void VarInstationary_Dim<Dimension>::InitSource()
  {    
    VirtualTimeSource<Real_wp>* fsrc = this->GetNewTimeSource();
    time_source.Init(fsrc, var_problem.GetFrequency(),
		     epsilon_time, this->initial_time);
    
    if (!this->SourceDoesNotDependOnTime())
      {
	var_source.InitIncidentField();
	
	IncidentWaveField<Real_wp, Dimension>* incident_wave;
	incident_wave = var_source.GetIncidentField(0, Real_wp(0));
	incident_wave->SetTimeSource(fsrc);

	IncidentWaveProjector<Real_wp, Dimension>* incident_proj;
	incident_proj = var_source.GetIncidentWaveProjector(0, Real_wp(0));
	incident_proj->Init(Real_wp(0));
	
	VirtualSourceFEM<Real_wp, Dimension>* f;
	f = var_source.GetNewSourceEquationObject(0);
	
	tlimit_source +=
	  f->Init(Real_wp(0), this->GetTimeStep(), var_problem.print_level, 0, true);

	delete f;
      }
  }

    
  //! adds source for the scalar unknowns only
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& t, int nb_deriv, VectReal_wp& b_src)
  {
    if (t >= tlimit_source)
      return;
    
    if (SourceDoesNotDependOnTime())
      {
        Real_wp pulse(0);
        SourceOnlyTime(t, nb_deriv, pulse);
        
        Real_wp alpha_ = alpha*pulse;
	int Nu = this->GetNbScalarUnknowns();
        for (int i = 0; i < sparse_vector_source.GetM(); i++)
          if (sparse_vector_source.Index(i) < Nu)
            b_src(sparse_vector_source.Index(i)) += alpha_*sparse_vector_source.Value(i);
      }
    else
      {
        Vector<VectReal_wp> b_vec;
        b_vec.SetData(1, &b_src);
        
	VirtualSourceFEM<Real_wp, Dimension>* f;
	f = var_source.GetNewSourceEquationObject(0);
	
	f->Init(t, this->GetTimeStep(), var_problem.print_level, nb_deriv, true);
	
	Vector<VirtualSourceFEM<Real_wp, Dimension>* > f_vec(1);
	f_vec(0) = f;
	var_source.AddSurfaceSource(alpha, b_vec, f_vec);
        b_vec.Nullify();
      }
  }

  
  //! adds source for the vectorial unknowns only
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::AddVectorialSourceAtTime(const Real_wp& alpha, const Real_wp& t, int nb_deriv,
			     VectReal_wp& b_src)
  {
    if (t >= tlimit_source)
      return;
    
    if (SourceDoesNotDependOnTime())
      {
	Real_wp pulse(0);
	SourceOnlyTime(t, nb_deriv, pulse);
	
	Real_wp alpha_ = alpha*pulse;
	int Nu = this->GetNbScalarUnknowns();
	for (int i = 0; i < sparse_vector_source.GetM(); i++)
	  if (sparse_vector_source.Index(i) >= Nu)
	    b_src(sparse_vector_source.Index(i)-Nu) += alpha_*sparse_vector_source.Value(i);
      }
    else
      {
        Vector<VectReal_wp> b_vec;
        b_vec.SetData(1, &b_src);
        
	VirtualSourceFEM<Real_wp, Dimension>* f;
	f = var_source.GetNewSourceEquationObject(0);
	
	Vector<VirtualSourceFEM<Real_wp, Dimension>* > f_vec(1);
	f_vec(0) = f;
	
	f->Init(t, this->GetTimeStep(), var_problem.print_level, nb_deriv, false);
	
	var_source.AddSurfaceSource(alpha, b_vec, f_vec);
        b_vec.Nullify();        	
      }
  }
      

  //! Inhomogeneous Dirichlet condition (unassembled)
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::SetDirichletConditionSource(const Real_wp& t, int nb_deriv, VectReal_wp& b_src,
				Real_wp alpha)
  {
    // first setting to 0 for homogeneous Dirichlet
    var_boundary.ImposeNullDirichletCondition(b_src);
    
    // then taking into account inhomogeneous Dirichlet terms
    if (SourceDoesNotDependOnTime())
      {
	int nb_dir = dirichlet_source.GetM();
	if (nb_dir <= 0)
	  return;
	
	Real_wp pulse(0); SourceOnlyTime(t, nb_deriv, pulse); 
	pulse *= alpha;
        for (int i = 0; i < dirichlet_source_unassembled.GetM(); i++)
	  b_src(dirichlet_source_unassembled.Index(i)) = pulse*dirichlet_source_unassembled.Value(i);
      }
    else
      {
	if (var_problem.FirstOrderFormulationDG() || (var_boundary.GetNbGlobalDirichletDof() <= 0))
	  return;
	
	Vector<VectReal_wp> b_vec;
        b_vec.SetData(1, &b_src);
        
	VirtualSourceFEM<Real_wp, Dimension>* f;
	f = var_source.GetNewSourceEquationObject(0);

	f->Init(t, this->GetTimeStep(), var_problem.print_level, nb_deriv, true);
	
	Vector<VirtualSourceFEM<Real_wp, Dimension>* > f_vec(1);
	f_vec(0) = f;
	var_source.SetSurfaceProjection(alpha, b_vec, f_vec);
        b_vec.Nullify();
      }
  }
  
  
  //! adds scalar source (ode interface)
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::AddScalarTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y)
  {
    this->AddScalarSourceAtTime(alpha, tn, n, Y);
  }


  //! adds vectorial source (ode interface)
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::AddVectorialTimeSource(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y)
  {
    this->AddVectorialSourceAtTime(alpha, tn, n, Y);
  }
  

  //! writes snapshots for the vectorial unknown 
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::WriteVectorialSnapshot(int& nb_iter, const Real_wp& t, VectReal_wp& Vh)
  {
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;
    
    bool snapshot = false;
    for (int i = 0; i < var_output.output_grid_param.GetM(); i++)
      if (var_output.output_grid_param(i).SnapshotToStore(t))
        snapshot = true;

    for (int i = 0; i < var_output.output_mesh_param.GetM(); i++)
      if (var_output.output_mesh_param(i).SnapshotToStore(t))
        snapshot = true;
    
    if (!snapshot)
      return;
    
    int nb_vec = var_problem.nb_unknowns_vec;
    VectReal_wp U(var_problem.offset_dof_unknown(var_problem.nb_unknowns_scal));
    //int Nvol = var_problem.mesh_num.GetNbDof();
    U.Fill(0);
    for (int k = 0; k < Dimension::dim_N; k++)
      {
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          {
            int offset = var_problem.GetOffsetDofV(i);
            for (int m = 0; m < var_problem.nb_unknowns_scal; m++)
              {
                int nm = var_problem.mesh_num_unknown(m);
                const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(nm);
                for (int j = 0; j < mesh_num.Element(i).GetNbDof(); j++)
                  {
                    int num_dof = mesh_num.Element(i).GetNumberDof(j);
                    U(num_dof + var_problem.offset_dof_unknown(m))
                      = Vh(offset + nb_vec*j + m*Dimension::dim_N + k);
                  }
              }
          }
        
        var_output.WriteSnapshot(nb_iter, t, U, k+2);
      }
  }


  //! checking that jacobian matrices are not variable in PML layers
  template<class Dimension>
  void VarInstationary_Dim<Dimension>::CheckIdentityJacobianInPML()
  {
    if (var_boundary.GetNbEltPML() > 0)
      if (var_problem.Glob_DFjm1.GetM() <= 0)
        {
          cout << "Jacobian matrices must be stored for checking them" << endl;
          abort();
        }
    
    int nb_elt = var_problem.mesh.GetNbElt();
    for (int i = 0; i < nb_elt; i++)
      if (var_problem.InsidePML(i))
        {
          if (Dimension::dim_N == 2)
            if (var_problem.Glob_DFjm1(i).GetM() > 1)
              {
                cout << "Expecting affine elements in PML layers" << endl;          
                abort();
              }
        }
  }
  

  //! what is the initial guess ?
  template<class Dimension>
  void VarInstationary_Dim<Dimension>
  ::SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y)
  {
    // only type of initial function => USER
    if (this->param_initial_condition.GetM() > 0)
      {
	VirtualProjectorFEM<Real_wp, Dimension>* initial_cond;
	initial_cond = GetNewInitialCondition();
	
	Vector<VectReal_wp> Yvec(1);
	Yvec(0).Reallocate(var_problem.GetNbDof());
	Yvec(0).Zero();
	
	Vector<VirtualProjectorFEM<Real_wp, Dimension>* > f_cond(1);
	f_cond(0) = initial_cond;
	var_source.AddVolumeProjection(1.0, Yvec, f_cond);	
	
	if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
	  {
	    int nb_scal = var_problem.nb_unknowns_scal;
	    int nb_vec = var_problem.nb_unknowns_vec;
	    int offset = var_problem.offset_dof_unknown(nb_scal);
            const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
	    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
	      {
		int nb_dof_elt = mesh_num.GetNbLocalDof(i);
		for (int j = 0; j < nb_dof_elt; j++)
		  for (int m = nb_scal; m < var_problem.nb_unknowns; m++)
		    Y(offset + j*nb_vec + m-nb_scal)
		      = Yvec(0)(mesh_num.Element(i).GetNumberDof(j) 
				+ m*mesh_num.GetNbDof());
		
		offset += nb_dof_elt*nb_vec;
	      }
	    
	  }
	else if (var_problem.FormulationDG() == ElementReference_Base::HDG)
	  {
	    int nb_dof_L = var_problem.nb_unknowns_hdg*var_problem.GetNbMeshDof();
	    int nb_dof_uv = var_problem.GetNbDof() - nb_dof_L;
	    for (int i = 0; i < nb_dof_uv; i++)
	      Y(i) = Yvec(0)(nb_dof_L + i);
	  }
	else
	  Y = Yvec(0);
	
	delete initial_cond;
      }
  }

  
  //! computes the local time steps for all elements of the mesh
  template<class Dimension>
  void VarInstationary_Dim<Dimension>::ComputeLocalTimeStep(VectReal_wp& local_dt)
  {
    Montjoie::ComputeLocalTimeStep(*this, local_dt);
  }
  

  //! writes the mesh with a difference reference for elements of num_implicit
  template<class Dimension>
  void VarInstationary_Dim<Dimension>::DisplayCoarseFineRegion(const IVect& num_implicit)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    string suffix;
    if (nb_proc > 1)
      suffix = "_" + to_str(rank_proc) + "N" + to_str(nb_proc);
    
    num_implicit.WriteText("num" + suffix + ".dat");
    Mesh<Dimension> mesh_order(var_problem.mesh);
    
    for (int i = 0; i < mesh_order.GetNbElt(); i++)
      mesh_order.Element(i).SetReference(1);
    
    for (int i = 0; i < num_implicit.GetM(); i++)
      mesh_order.Element(num_implicit(i)).SetReference(2);
    
    mesh_order.Write("order" + suffix + ".mesh");    
  }
  
  
  /***********************
   * VarInstationary_Fem *
   ***********************/

  
  //! returns norm of the solution
  /*!
    \param[in] Uh vector solution
   */
  template<class Dimension>
  Real_wp VarInstationary_Fem<Dimension>::GetNormeSolution(const VectReal_wp& Uh) const
  {
    Real_wp norme_sol(0);
    const Mesh<Dimension>& mesh = var_problem.mesh;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb = var_problem.GetReferenceElement(i);
	int nb_dof_elt = Fb.GetNbDof();
	int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
	IVect Nodle = var_problem.GetDofNumberOnElement(i);
	int nb_unknowns = var_problem.nb_unknowns;
	int Nvol = var_problem.offset_dof_unknown(1);
	nb_unknowns = min(nb_unknowns, int(Uh.GetM()/Nvol));
	// this is not actually the true L^2 norm
	// but an approximation (which is worthy if quadrature weights are positive...)
	// which is faster to compute (no interpolation is needed)
	if (var_problem.Glob_jacobian(i).GetM() > 1)
	  {
	    for (int m = 0; m < nb_unknowns; m++)
	      for (int j = 0; j < min(nb_dof_elt,nb_pts_quad); j++)
		{
		  int num_dof = Nodle(j) + var_problem.offset_dof_unknown(m);
		  norme_sol += var_problem.Glob_jacobian(i)(j)*Uh(num_dof)*Uh(num_dof);
		}
	  }
	else
	  {
	    Real_wp jacob = var_problem.Glob_jacobian(i)(0);
	    for (int m = 0; m < nb_unknowns; m++)
	      for (int j = 0; j < min(nb_dof_elt,nb_pts_quad); j++)
		{
		  int num_dof = Nodle(j) + var_problem.offset_dof_unknown(m);
		  norme_sol += jacob*Fb.WeightsND(j)*Uh(num_dof)*Uh(num_dof);
		}
	  }
      }
    
    Real_wp sum(0);
#ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    MpiReduce(var_problem.comm_group_mode, &norme_sol, xtmp, &sum, 1, MPI_SUM, 0);
#else
    sum = norme_sol;
#endif
    
    if (abs(norme_sol) > this->max_norme_sol)
      {
        cout << "Solution instable ?" << endl;
        abort();
      }

    return sum;
  }
  
}

#define MONTJOIE_FILE_VAR_INSTATIONARY_CXX
#endif
