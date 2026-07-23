#ifndef MONTJOIE_FILE_VAR_PROBLEM_BASE_CXX

namespace Montjoie
{
  
  /*******************
   * VarProblem_Base *
   *******************/
  
  bool VarProblem_Base::exit_if_no_boundary_condition(true);
  
  //! Sets attributes to default values
  void VarProblem_Base::InitDefaultValues()
  {    
    print_level = 0;
    
    dg_formulation = ElementReference_Base::CONTINUOUS;    
    if (first_order_formulation)
      dg_formulation = ElementReference_Base::DISCONTINUOUS;    

    sipg_formulation = true;
    compute_dfjm1 = false;

    alpha_penalization = 0.0;
    delta_penalization = 0.0;
    automatic_choice_penalization = true;
    upwind_fluxes = true;

    if (first_order_formulation)
      {
	alpha_penalization = -Real_wp(1);
	delta_penalization = -Real_wp(1);
      }
    
    // by default a unit brick
    xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0;
    zmin = 0.0; zmax = 1.0;
 
    nodl = 0;
    default_order = 1;

    order_over_integration = 0;
    dg_exact_integration = true;
    mixed_formulation = false;

    frequency = 1.0;
    omega = 2.0*pi_wp;
    omega2 = omega*omega;    
    wavelength_adim = 1.0;
    
    if (this->nb_unknowns <= 0)
      {
	cout<<"The number of unknowns is null "<<endl;
	cout<<"Fix this problem, then restart "<<endl;
	abort();
      }
  }
  

  //! Sets the pulsation of the problem
  void VarProblem_Base::SetOmega(const Real_wp& omega_)
  {
    omega = omega_;
    omega2 = omega*omega;
    frequency = omega/(2.0*pi_wp);
    this->UpdateWaveVector();
  }
  

  //! sets the frequency of the problem
  void VarProblem_Base::SetFrequency(const Real_wp& freq_)
  {
    frequency = freq_;
    omega = frequency*2.0*pi_wp;
    omega2 = omega*omega;
    this->UpdateWaveVector();
  }

  
  //! Modifies the object with a line of the data file
  void VarProblem_Base::SetInputData(const string& description_field,
				     const VectString& parameters)
  {    
    if (!description_field.compare("ExactIntegration"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "ExactIntegration needs more parameters, for instance :" << endl;
	    cout << "ExactIntegration = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	dg_exact_integration = false;
	if (!parameters(0).compare("YES"))
	  dg_exact_integration = true;
	
	if (parameters.GetM() > 1)
	  order_over_integration = to_num<int>(parameters(1));
      }
    else if (description_field == "MixedFormulation")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "MixedFormulation needs more parameters, for instance :" << endl;
	    cout << "MixedFormulation = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (parameters(0) == "YES")
          this->SetFirstOrderFormulation(true);
        else
          this->SetFirstOrderFormulation(false);
      }
    else if (!description_field.compare("Exit_IfNo_BoundaryCondition"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "Exit_IfNo_BoundaryCondition needs more parameters, for instance :" << endl;
	    cout << "Exit_IfNo_BoundaryCondition = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("YES"))
	  this->exit_if_no_boundary_condition = true;
	else
	  this->exit_if_no_boundary_condition = false;
      }
    else if (description_field == "StoreDFjm1")
      {
	if (parameters(0) == "YES")
	  this->compute_dfjm1 = true;
	else
	  this->compute_dfjm1 = false;
      }
    else if (!description_field.compare("FileMesh")) 
      {	
	mesh_data.Clear();
	mesh_data.PushBack(parameters);
      }
    else if (!description_field.compare("AdditionalMesh")) 
      {
	mesh_data.PushBack(parameters);
      }
    else if (!description_field.compare("PenalizationDG"))
      {
        if (parameters.GetM() <= 0)
          {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "PenalizationDG needs more parameters, for instance :" << endl;
	    cout << "Penalization = Upwind" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }
        
        if (parameters(0) == "Upwind")
          {
            upwind_fluxes = true;

	    if (this->FormulationDG() == ElementReference_Base::CONTINUOUS)
	      {
		alpha_penalization = Real_wp(0);
                delta_penalization = Real_wp(0);
	      }
            else if (this->first_order_formulation)
              {
                alpha_penalization = Real_wp(-1);
                delta_penalization = Real_wp(-1);
              }
            else
              {
                // setting alpha and delta to an high value for SIPG
                alpha_penalization = Real_wp(-10);
                delta_penalization = Real_wp(-10);
              }
          }
        else
          {
            if (parameters.GetM() <= 1)
              {
                cout << "In SetInputData of VarProblem" << endl;
                cout << "PenalizationDG needs more parameters, for instance :" << endl;
                cout << "Penalization = alpha delta" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }        
            
            upwind_fluxes = false;
            alpha_penalization = to_num<Real_wp>(parameters(0));
            delta_penalization = to_num<Real_wp>(parameters(1));
          }
      }
    else if (!description_field.compare("CoefficientPenalization"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarProblem" << endl;
	    cout << "CoefficientPenalization needs more parameters, for instance :" << endl;
	    cout << "CoefficientPenalization = AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }
        
        if (parameters(0) == "AUTO")
          automatic_choice_penalization = true;
        else
          automatic_choice_penalization = false;        
      }
    else if (!description_field.compare("MateriauDielec"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "MateriauDielec needs more parameters, for instance :" << endl;
	    cout << "MateriauDielec = ref ISOTROPE value" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int ref = to_num<int>(parameters(0));
	Vector<string> parameters_(parameters.GetM()-1);
	for (int i = 1; i < parameters.GetM(); i++)
	  parameters_(i-1) = parameters(i);
	
	this->SetIndices(ref, parameters_);
      }
    else if (!description_field.compare("PhysicalMedia"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "PhysicalMedia needs more parameters, for instance :" << endl;
	    cout << "PhysicalMedia = nom_media ref value" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        string name_media = parameters(0);
        int ref = to_num<int>(parameters(1));
	Vector<string> parameters_(parameters.GetM()-2);
	for (int i = 2; i < parameters.GetM(); i++)
	  parameters_(i-2) = parameters(i);
	
	this->SetPhysicalIndex(name_media, ref, parameters_);        
      }
    else if (description_field == "UnknownMeshNumbering")
      {
        if (parameters.GetM() > this->nb_unknowns)
          {
            cout << "Number of parameters " << parameters.GetM() << " is greater than the number of unknowns : " << this->nb_unknowns << endl;
            abort();
          }
        
        for (int i = 0; i < parameters.GetM(); i++)
          mesh_num_unknown(i) = to_num<int>(parameters(i));
      }
    else if (!description_field.compare("OrderGeometry"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MeshNumbering_Base" << endl;
	    cout << "OrderGeometry needs more parameters, for instance :" << endl;
	    cout << "OrderGeometry = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int r = to_num<int>(parameters(0));
	default_order = r;
      }
    else if (!description_field.compare("OrderDiscretization"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MeshNumbering_Base" << endl;
	    cout << "OrderDiscretization needs more parameters, for instance :" << endl;
	    cout << "OrderDiscretization = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if ((!parameters(0).compare("MAX_EDGE")) || (!parameters(0).compare("MEAN_EDGE")))
	  {
	  }
	else
	  {
	    int r = to_num<int>(parameters(0));
	    default_order = r;
	  }
      }
    else if (!description_field.compare("NonLinearSolver"))
      {
#ifdef MONTJOIE_WITH_TWO_DIM
        FjInverseProblem<Dimension2>::SetInputData(description_field, parameters);
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
        FjInverseProblem<Dimension3>::SetInputData(description_field, parameters);
#endif
      }
    else if (!description_field.compare("PrintLevel"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "PrintLevel needs more parameters, for instance :" << endl;
	    cout << "PrintLevel = level" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	print_level = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("Frequency"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "Frequency needs more parameters, for instance :" << endl;
	    cout << "Frequency = a b" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	this->omega = to_num<Real_wp>(parameters(0))*pi_wp*2.0
	  + to_num<Real_wp>(parameters(1));

	if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
	  this->omega *= this->wavelength_adim / PhysicalConstant::speed_light;
	
	// omega is the pulsation
	// omega2 is the square of omega
	this->omega2 = this->omega*this->omega;
	this->frequency = this->omega/(2.0*pi_wp);
      }
    else if ((description_field == "Wavelength") || (description_field == "WaveLength"))
      {
        // wavelength in the vacuum
	// Wavelength = lambda
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarProblemBase" << endl;
	    cout << "Wavelength needs one parameter, for instance :" << endl;
	    cout << "Wavelength = lambda" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

	Real_wp lambda = to_num<Real_wp>(parameters(0));
	Real_wp f;
	if (PhysicalConstant::adimensionalization != PhysicalConstant::ADIM_YES)
	  f = Real_wp(1)/lambda;
	else
	  f = wavelength_adim / lambda;

	SetFrequency(f);
      }
    else if (!description_field.compare("PhysicalFrequency"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "PhysicalFrequency needs more parameters, for instance :" << endl;
	    cout << "PhysicalFrequency = f0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	to_num(parameters(0), this->frequency);
	this->frequency *= Real_wp(1)/PhysicalConstant::speed_light; 
	this->omega = 2.0*pi_wp*this->frequency; 
	
	// omega is the pulsation
	// omega2 is the square of omega
	this->omega2 = this->omega*this->omega;
      }
    else if (description_field == "WavelengthAdim")
      {
	wavelength_adim = to_num<Real_wp>(parameters(0));
      }    
  }
  

  //! adds memory used by the current object
  void VarProblem_Base::GetMemoryUsed(map<string, size_t>& var) const
  {
  }


  //! computes the velocity of waves for each element of the mesh
  template<class Dimension>
  void VarProblem_Base::GetVelocityOnElements(VectReal_wp& velocity,
					      const Mesh<Dimension>& mesh)
  {
    // we assume that media are quasi-uniform on each reference
    // therefore we can associate a "mean" velocity for each reference
    int nb_domains = this->GetNbPhysicalIndices();
    VectReal_wp velocity_media(nb_domains);
    velocity_media.Fill(1.0);
    for (int ref = 0; ref < nb_domains; ref++)
      velocity_media(ref) = this->GetVelocityOfMedia(ref);
    
    velocity.Reallocate(mesh.GetNbElt());
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int ref = mesh.Element(i).GetReference();
        velocity(i) = velocity_media(ref);
      }
  }

  
  //! returns the velocity of waves for the physical domain whose reference is ref
  /*!
    By default, we consider that waves are travelling at a velocity equal to 1
    This method is usually overloaded in the derived classes
   */
  Real_wp VarProblem_Base::GetVelocityOfMedia(int ref) const
  {
    return Real_wp(1);
  }
    

  //! returns the velocity of waves at infinity
  /*!
    By default, we consider that waves are travelling at a velocity equal to 1
    This method is usually overloaded in the derived classes
   */
  Real_wp VarProblem_Base::GetVelocityOfInfinity() const
  {
    return Real_wp(1);
  }

  
  //! computes coefficients tau for HDG formulation
  void VarProblem_Base::ComputeTauCoefficient()
  {
    // default case : no HDG
  }

  
  //! copies input parameters of another similar problem
  void VarProblem_Base::CopyInputData(const VarProblem_Base& var)
  {
    print_level = var.print_level;
    
    dg_formulation = var.dg_formulation;
    sipg_formulation = var.sipg_formulation;
    alpha_penalization = var.alpha_penalization;
    delta_penalization = var.delta_penalization;
    upwind_fluxes = var.upwind_fluxes;
    mesh_num_unknown = var.mesh_num_unknown;
    name_other_elements = var.name_other_elements;
    
    automatic_choice_penalization = var.automatic_choice_penalization;
    mesh_data = var.mesh_data;
    finite_element_name = var.finite_element_name;
    
    default_order = var.default_order;
    dg_exact_integration = var.dg_exact_integration;
    mixed_formulation = var.mixed_formulation;
    
    omega = var.omega;
    omega2 = var.omega2;
    frequency = var.frequency;
  }


  //! sets the name of the equation that will be solved (can be overloaded)
  void VarProblem_Base::SetTypeEquation(const string&)
  {
    if (this->dg_formulation == ElementReference_Base::HDG)
      mesh_num_unknown.Reallocate(nb_unknowns_hdg + nb_unknowns);
    else
      mesh_num_unknown.Reallocate(nb_unknowns);
    
    mesh_num_unknown.Zero();
    this->InitPolarization();
  }


  /*************************
   * VarComputationProblem *
   *************************/
  
  
  //! Default constructor
  VarComputationProblem_Base::VarComputationProblem_Base()
  {
    threshold_matrix = 1e-30;
  }

  
  //! Sets values of attributes to default values
  void VarComputationProblem::InitDefaultValues()
  {
    symmetric_elementary_matrix = false;
    symmetric_global_matrix = false;
    sparse_elementary_matrix = false;
    
    storage_finite_element_matrix = MATRIX_AUTO_STORAGE;

    leaf_static_condensation = false;
    light_static_condensation = false;
    use_symmetrization_when_possible = false;
    erase_dirichlet_columns = false;
  }

  
  //! modifies the object with a line of the data file
  void VarComputationProblem
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "ExplicitMatrixFEM")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "ExplicitMatrixFEM needs more parameters, for instance :" << endl;
	    cout << "ExplicitMatrixFEM = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (parameters(0) == "YES")
          this->storage_finite_element_matrix = MATRIX_STORED;
        else if (parameters(0) == "NO")
          this->storage_finite_element_matrix = MATRIX_FREE;
        else
          this->storage_finite_element_matrix = MATRIX_AUTO_STORAGE;
      }
    else if (description_field == "TypeCondensation")
      {
	if (parameters(0) == "Light")
	  this->light_static_condensation = true;
	else
	  this->light_static_condensation = false;
      }
    else if (!description_field.compare("ThresholdMatrix"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "ThresholdMatrix needs more parameters, for instance :" << endl;
	    cout << "ThresholdMatrix = epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	this->threshold_matrix = to_num<Real_wp>(parameters(0));
      }
  }
  

  //! sets dofs for which a static condensation can be performed
  /*!
    \param[in] i element number
    \param[in] nb_dof_loc number of degrees of freedom associated with element i
    \param[out] nb_dof_edges number of degrees of freedom that can not be eliminated
    \param[out] nb_dof_int number of degrees of freedom that can be eliminated
    \param[out] intern_node intern_node(i) >= 0 for dofs on the boundary, in that
                case intern_node(i) is the local number among dofs on the boundary.
		intern_node(i) < 0 for intern dofs, in that case -intern_node(i)-1
		is the local number among intern dofs
   */
  void VarComputationProblem::GetInternalNodesElement(int i, int nb_dof_loc,
						      int& nb_dof_edges, int& nb_dof_int,
						      Vector<int>& intern_node) const
  {
    int nb_u = var_problem.nb_unknowns;
    int nb_mesh_num = var_problem.GetNbMeshNumberings();
    IVect nb_dof_elt(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      nb_dof_elt(nm) = var_problem.GetNbLocalDof(i, nm);
    
    intern_node.Reallocate(nb_dof_loc);
    intern_node.Fill(-1);

    if (this->light_static_condensation)
      {
	nb_dof_edges = 0;
        for (int n = 0; n < nb_u; n++)
          nb_dof_edges += nb_dof_elt(var_problem.mesh_num_unknown(n));
        
	nb_dof_int = nb_dof_loc - nb_dof_edges;
        
	for (int i = 0; i < nb_dof_edges; i++)
	  intern_node(i) = i;

	for (int i = nb_dof_edges; i < nb_dof_loc; i++)
	  intern_node(i) = -(i-nb_dof_edges)-1;
	
	return;
      }

    IVect nb_dof_bound(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      nb_dof_bound(nm) = var_problem.GetNbDofBoundaries(i, nm);

    nb_dof_edges = 0;
    for (int n = 0; n < nb_u; n++)
      nb_dof_edges += nb_dof_bound(var_problem.mesh_num_unknown(n));
    
    nb_dof_int = nb_dof_loc - nb_dof_edges;
    
    int nb = 0, offset = 0;
    for (int m = 0; m < nb_u; m++)
      {
        int nm = var_problem.mesh_num_unknown(m);
        for (int j = 0; j < nb_dof_bound(nm); j++)
          intern_node(offset + j) = nb++;
        
        offset += nb_dof_elt(nm);
      }
    
    nb = 1;
    for (int j = 0; j < nb_dof_loc; j++)
      if (intern_node(j) < 0)
	intern_node(j) = -nb++;
  }


  //! computation of the direct matrix
  /*!
    \param[out] mat_iterative objet needed to handle inhomogeneous Dirichlet condition
    \param[out] mat_direct the finite element matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods computes the finite element matrix with boundary conditions
  */
  template<class T, class MatrixSparse>
  void VarComputationProblem
  ::AddMatrixWithBC(FemMatrixFreeClass_Base<T>& mat_iterative,
		    MatrixSparse& mat_sp, const GlobalGenericMatrix<T>& nat_mat,
		    int offset_row, int offset_col,
		    CondensationBlockSolver_Fem<T>* solver, bool diag_matrix)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    mat_iterative.SetSymmetrizationUse(this->use_symmetrization_when_possible);
    int nodl = var_problem.GetNbDof();
    
    // we change size of the matrix if necessary
    mat_iterative.Reallocate(nodl, nodl);    
    if (mat_sp.GetM() < nodl)
      mat_sp.Resize(nodl, nodl);
    
    // boundary conditions and other extra-terms
    var_boundary.AddBoundaryConditionTerms(mat_sp, nat_mat, offset_row, offset_col);
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 4)
	cout << rank_proc << " Boundary Conditions OK" << endl;

    bool sym = var_problem.IsSymmetricProblem();
    
    if (sym)
      mat_iterative.InitSymmetricMatrix();
    else
      mat_iterative.InitUnsymmetricMatrix();
    
    // main interactions of the matrix (volumic integrals)
    this->AddMatrixFEM(mat_iterative, mat_sp, nat_mat, offset_row, offset_col,
                       solver, diag_matrix);
  }
  

  //! computation of the direct matrix
  /*!
    \param[out] mat_sp the finite element matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods computes the finite element matrix with boundary conditions
  */
  template<class T, class Prop, class Storage>
  void VarComputationProblem
  ::AddMatrixWithBC(DistributedMatrix<T, Prop, Storage>& mat_sp,
                    const GlobalGenericMatrix<T>& nat_mat,
		    int offset_row, int offset_col,
		    CondensationBlockSolver_Fem<T>* solver,
                    bool diag_matrix)
  {
    FemMatrixFreeClass_Base<T>* free_mat = GetNewIterativeMatrix(T(0));
    AddMatrixWithBC(*free_mat, mat_sp, nat_mat, offset_row, offset_col,
                    solver, diag_matrix);
    
    delete free_mat;
  }


  //! computation of the direct matrix
  /*!
    \param[out] mat_sp the finite element matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods computes the finite element matrix with boundary conditions
  */
  template<class T, class Prop>
  void VarComputationProblem
  ::AddMatrixWithBC(Matrix<T, Prop, DiagonalRow>& mat_sp,
                    const GlobalGenericMatrix<T>& nat_mat,
		    int offset_row, int offset_col,
		    CondensationBlockSolver_Fem<T>* solver,
                    bool diag_matrix)
  {
    FemMatrixFreeClass_Base<T>* free_mat = GetNewIterativeMatrix(T(0));
    AddMatrixWithBC(*free_mat, mat_sp, nat_mat, offset_row, offset_col,
                    solver, diag_matrix);
    
    delete free_mat;
  }
  
  
  //! computation of iterative matrix with associated boundary conditions
  /*!
    \param[out] mat_iterative iterative matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    This methods computes the matrix necessary for the iterative solver
    (it can be a "matrix-free" form)
   */
  template<class T>
  void VarComputationProblem
  ::AddMatrixWithBC(FemMatrixFreeClass_Base<T>& mat_iterative,
		    const GlobalGenericMatrix<T>& nat_mat,
		    int offset_row, int offset_col,
		    CondensationBlockSolver_Fem<T>* solver,
                    bool diag_matrix)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    int nodl = var_problem.GetNbDof();    
    bool sym = var_problem.IsSymmetricProblem();
    mat_iterative.SetSymmetrizationUse(this->use_symmetrization_when_possible);
    mat_iterative.Reallocate(nodl, nodl);

    bool presence_boundary_term = true;
    if (var_problem.FormulationDG() != ElementReference_Base::CONTINUOUS)
      {
        // for LDG or HDG formulation, boundary conditions are set elsewhere
        if (var_problem.FirstOrderFormulationDG())
          presence_boundary_term = false;
      }

    if (sym)
      {
        mat_iterative.InitSymmetricMatrix();
	
	// initialization of the iterative matrix and boundary terms
	if ((mat_iterative.mat_boundary_sym.GetM() < nodl) && (presence_boundary_term))
	  mat_iterative.mat_boundary_sym.Resize(nodl, nodl);
	
	// boundary conditions and extra-terms are added
	var_boundary.AddBoundaryConditionTerms(mat_iterative.mat_boundary_sym,
					       nat_mat, offset_row, offset_col);
      }
    else
      {
        mat_iterative.InitUnsymmetricMatrix();
        
	// initialization of the iterative matrix and boundary terms
	if ((mat_iterative.mat_boundary_unsym.GetM() < nodl) && (presence_boundary_term))
	  mat_iterative.mat_boundary_unsym.Resize(nodl, nodl);
	
	// boundary conditions and extra-terms are added
	var_boundary.AddBoundaryConditionTerms(mat_iterative.mat_boundary_unsym,
					       nat_mat, offset_row, offset_col);
      }

    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 4)
        cout << rank_proc << " Boundary Conditions OK" << endl;

    // volumetric terms are added
    this->AddMatrixFEM(mat_iterative, nat_mat, offset_row, offset_col,
                       solver, diag_matrix);    
  }
  
    
  //! computation of the finite element matrix (only interior terms)
  /*!
    \param[in,out] mat_iterative finite element matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods adds "volumic" terms to the matrix mat_direct
    The previous entries stored in the matrix are kept
   */
  template<class T>
  void VarComputationProblem
  ::AddMatrixFEM(FemMatrixFreeClass_Base<T>& mat_iterative,
		 const GlobalGenericMatrix<T>& nat_mat, int offset_row, int offset_col,
		 CondensationBlockSolver_Fem<T>* solver, bool diag_matrix)
  {
    mat_iterative.SetCoefficientMatrix(nat_mat);
    
    bool store_matrix = false;
    if (leaf_static_condensation)
      store_matrix = true;
    
    if (this->storage_finite_element_matrix == this->MATRIX_STORED)
      store_matrix = true;
    
    int nodl = var_problem.GetNbDof();
    if (store_matrix)
      {
        if (IsSymmetricMatrix(mat_iterative))
          {
	    if (mat_iterative.mat_boundary_sym.GetM() > 0)
	      {
		mat_iterative.mat_iterative_sym = mat_iterative.mat_boundary_sym;
		mat_iterative.mat_boundary_sym.Clear();
	      }
	    else
	      mat_iterative.mat_iterative_sym.Resize(nodl, nodl);
	    
	    this->AddMatrixFEM(mat_iterative, mat_iterative.mat_iterative_sym,
			       nat_mat, offset_row, offset_col, solver);
	  }
	else
	  {
	    if (mat_iterative.mat_boundary_unsym.GetM() > 0)
	      {
		mat_iterative.mat_iterative_unsym = mat_iterative.mat_boundary_unsym;
		mat_iterative.mat_boundary_unsym.Clear();
	      }
	    else
	      mat_iterative.mat_iterative_unsym.Resize(nodl, nodl);
		
	    this->AddMatrixFEM(mat_iterative, mat_iterative.mat_iterative_unsym,
			       nat_mat, offset_row, offset_col, solver, diag_matrix);
	  }
      }
    else
      {
	mat_iterative.SetDirichletCondition(offset_row, offset_col, erase_dirichlet_columns);
      }
  }
  
  
  //! computation of the finite element matrix (only interior terms)
  /*!
    \param[out] mat_iterative objet needed to handle inhomogeneous Dirichlet condition
    \param[in,out] mat_direct sparse matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods adds "volumic" terms to the matrix mat_direct
    The previous entries stored in the matrix are kept
   */
  template<class T, class MatrixSparse>
  void VarComputationProblem
  ::AddMatrixFEM(FemMatrixFreeClass_Base<T>& mat_iterative,
		 MatrixSparse& mat_direct,
		 const GlobalGenericMatrix<T>& nat_mat,
		 int offset_row, int offset_col,
		 CondensationBlockSolver_Fem<T>* solver, bool diag_matrix)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    var_comm.InitDistributedMatrix(mat_direct);
    
    if ((rank_proc == 0) || (var_problem.print_level >= 10))
      if (var_problem.print_level >= 2)
	cout << rank_proc << " We compute direct matrix " << endl;
    
    // we add boundary integrals (e.g. numerical fluxes in DG)
    if (var_problem.FormulationDG() != ElementReference_Base::HDG)
      this->AddElementaryFluxesDG(mat_direct, nat_mat, offset_row, offset_col);
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 1)
	cout << rank_proc << " generating matrix data..." << endl;

    CondensationBlockSolver_Base<T>* condensed_solver;
    if (solver != NULL)
      condensed_solver = solver;
    else
      condensed_solver = new CondensationBlockSolver_Base<T>;
    
    if (diag_matrix || this->IsDiagonalElementaryMatrix(nat_mat))
      {
        if ((rank_proc == 0) && (var_problem.print_level >= 2))
          cout << "Diagonal elementary matrix " << endl;
        
        // diagonal elementary matrix
        Matrix<T, Symmetric, DiagonalRow> mat_elt;        

        // assembling the matrix, that will call ComputeElementaryMatrix for each element
	AssembleMatrix(mat_direct, mat_elt, nat_mat, *this, *condensed_solver, offset_row, offset_col);
      }
    else if (this->IsSymmetricElementaryMatrix(nat_mat))
      {
        if (this->IsSparseElementaryMatrix(nat_mat))
          {
            if ((rank_proc == 0) && (var_problem.print_level >= 2))
              cout << "Symmetric sparse elementary matrix " << endl;
            
            // symmetric elementary matrix
            Matrix<T, Symmetric, ArrayRowSymSparse> mat_elt;
            
            // assembling the matrix, that will call ComputeElementaryMatrix for each element
            AssembleMatrix(mat_direct, mat_elt, nat_mat, *this, *condensed_solver, offset_row, offset_col);
          }
        else
          {
            if ((rank_proc == 0) && (var_problem.print_level >= 2))
              cout << "Symmetric dense elementary matrix " << endl;
            
            // symmetric elementary matrix
            Matrix<T, Symmetric, RowSymPacked> mat_elt;
            
            // assembling the matrix, that will call ComputeElementaryMatrix for each element
            AssembleMatrix(mat_direct, mat_elt, nat_mat, *this, *condensed_solver, offset_row, offset_col);
          }
      }
    else
      {
	if (this->IsSparseElementaryMatrix(nat_mat))
          {
            if ((rank_proc == 0) && (var_problem.print_level >= 2))
              cout << "Unsymmetric sparse elementary matrix " << endl;
            
            // unsymmetric elementary matrix
            Matrix<T, General, ArrayRowSparse> mat_elt;
            
            // assembling the matrix, that will call ComputeElementaryMatrix for each element
            AssembleMatrix(mat_direct, mat_elt, nat_mat, *this, *condensed_solver, offset_row, offset_col);

          }
        else
          {
            if ((rank_proc == 0) && (var_problem.print_level >= 2))
              cout << "Unsymmetric dense elementary matrix " << endl;
            
            // unsymmetric elementary matrix
            Matrix<T> mat_elt;
            
            // assembling the matrix, that will call ComputeElementaryMatrix for each element
            AssembleMatrix(mat_direct, mat_elt, nat_mat, *this, *condensed_solver, offset_row, offset_col);
          }
      }
    
    if (solver == NULL)
      delete condensed_solver;
    
    // quasi-periodic conditions are enforced
    var_boundary.SetPeriodicCondition(mat_direct);
    
    // and finally Dirichlet condition    
    mat_iterative.SetDirichletCondition(mat_direct, offset_row, offset_col,
                                        erase_dirichlet_columns);
    
    // small entries are removed
    mat_direct.RemoveSmallEntry(this->GetThresholdMatrix());
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 1)
	cout << rank_proc << " The matrix takes " <<
          GetHumanReadableMemory(mat_direct.GetMemorySize()) << endl;
  }
  

  //! computation of the finite element matrix (only interior terms)
  /*!
    \param[in,out] mat_direct sparse matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods adds "volumic" terms to the matrix mat_direct
    The previous entries stored in the matrix are kept
   */
  template<class T, class Prop, class Storage>
  void VarComputationProblem
  ::AddMatrixFEM(DistributedMatrix<T, Prop, Storage>& mat_direct,
                 const GlobalGenericMatrix<T>& nat_mat,
		 int offset_row, int offset_col,
		 CondensationBlockSolver_Fem<T>* solver,
                 bool diag_matrix)
  {
    FemMatrixFreeClass_Base<T>* free_mat = GetNewIterativeMatrix(T(0));
    AddMatrixFEM(*free_mat, mat_direct, nat_mat, offset_row, offset_col,
                 solver, diag_matrix);
    
    delete free_mat;
  }


  //! computation of the finite element matrix (only interior terms)
  /*!
    \param[in,out] mat_direct sparse matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] offset_row offset to add to row numbers
    \param[in] offset_col offset to add to column numbers
    this methods adds "volumic" terms to the matrix mat_direct
    The previous entries stored in the matrix are kept
   */
  template<class T, class Prop>
  void VarComputationProblem
  ::AddMatrixFEM(Matrix<T, Prop, DiagonalRow>& mat_direct,
                 const GlobalGenericMatrix<T>& nat_mat,
		 int offset_row, int offset_col,
		 CondensationBlockSolver_Fem<T>* solver,
                 bool diag_matrix)
  {
    FemMatrixFreeClass_Base<T>* free_mat = GetNewIterativeMatrix(T(0));
    AddMatrixFEM(*free_mat, mat_direct, nat_mat, offset_row, offset_col,
                 solver, diag_matrix);
    
    delete free_mat;
  }
  

  //! computation of the diagonal of iterative matrix
  template<class T>
  void VarComputationProblem
  ::ComputeDiagonalMatrix(Vector<T>& diagonal,
			  const GlobalGenericMatrix<T>& nat_mat, bool assemble)
  { 
    int nodl = var_problem.GetNbDof();
    diagonal.Reallocate(nodl); diagonal.Fill(0);

    // computing the matrix by considering it as diagonal
    Matrix<T, Symmetric, DiagonalRow> diag_mat;
    diag_mat.SetData(diagonal.GetM(), diagonal.GetData());

    CondensationBlockSolver_Fem<T>* null_ptr_c = NULL;
    this->AddMatrixWithBC(diag_mat, nat_mat, 0, 0, null_ptr_c, true);
    
    diag_mat.Nullify();
    
    // assembling diagonal
    if (assemble)
      var_comm.AddDomains(diagonal);    
  }
  
  
  //! computation of the diagonal of iterative matrix
  /*!
    \param[out] diagonal diagonal of the iterative matrix
    \param[in] mat iterative matrix (eventually a matrix-free form)
    \param[in] nat_mat mass and stiffness coefficient
   */
  template<class T>
  void VarComputationProblem
  ::ComputeDiagonalMatrix(Vector<T>& diagonal,
			  const FemMatrixFreeClass_Base<T>& mat,
			  const GlobalGenericMatrix<T>& nat_mat, bool assemble)
  {
    int nodl = mat.GetM();
    diagonal.Reallocate(nodl); diagonal.Zero();
    
    bool matrix_not_stored = false;
    if (mat.IsSymmetric())
      {
	if (mat.mat_iterative_sym.GetM() > 0)
	  for (int j = 0; j < nodl; j++)
	    diagonal(j) = mat.mat_iterative_sym(j, j);
	else if (mat.matCSR_iterative_sym.GetM() > 0)
	  for (int j = 0; j < nodl; j++)
	    diagonal(j) = mat.matCSR_iterative_sym(j, j);
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (mat.mat_iterative_unsym.GetM() > 0)
	  for (int j = 0; j < nodl; j++)
	    diagonal(j) = mat.mat_iterative_unsym(j, j);
	else if (mat.matCSR_iterative_unsym.GetM() > 0)
	  for (int j = 0; j < nodl; j++)
	    diagonal(j) = mat.matCSR_iterative_unsym(j, j);
	else
	  matrix_not_stored = true;
      }

    if (matrix_not_stored)
      {
	// computing the matrix by considering it as diagonal
	Matrix<T, Symmetric, DiagonalRow> diag_mat;
	diag_mat.SetData(diagonal.GetM(), diagonal.GetData());

        CondensationBlockSolver_Fem<T>* null_ptr_c = NULL;
	this->AddMatrixWithBC(diag_mat, nat_mat, 0, 0, null_ptr_c, true);
	
	diag_mat.Nullify();
      }
    
    // assembling diagonal
    if (assemble)
      var_comm.AddDomains(diagonal);
  }
  
  
  //! computation of diagonal matrix for a standard matrix
  template<class T, class Prop, class Storage, class Allocator>
  void VarComputationProblem::
  ComputeDiagonalMatrix(Vector<T>& diagonal,
			const Matrix<T, Prop, Storage, Allocator>& mat,
			const GlobalGenericMatrix<T>& nat_mat)
  {
    int nodl = mat.GetM();
    diagonal.Reallocate(nodl);
    for (int i = 0; i < nodl; i++)
      diagonal(i) = mat(i, i);
  }
  
  
  //! computing a pattern for a block-diagonal matrix
  /*!
    The pattern is computed so that dofs from a same block
    are not geometrically far from each other
    \param[in] nblock wished size for each block
    \param[out] num_ddl_blocks dof numbers grouped by blocks
   */
  void VarComputationProblem
  ::GetProfilBlockDiagonal(int nblock, Vector<IVect>& num_ddl_blocks)
  {
    abort();
    /*
    IVect liste_elem(1), pos_elem(1); liste_elem(0) = 0; pos_elem(0) = 0;
    IVect new_liste_elem, new_pos_elem;
    VectBool DofUsed(this->mesh_num.GetNbDof()); DofUsed.Fill(false);
    VectBool ElementUsed(this->mesh.GetNbElt()); ElementUsed.Fill(false); ElementUsed(0) = true;
    VectBool LocalDofUsed; IVect ListeLocalDof;
    IVect num_ddl; num_ddl_blocks.Reallocate(this->nodl);
    int nb_blocks = 0, size_all = 0;
    while (size_all < this->nb_unknowns*this->mesh_num.GetNbDof())
      {
	int size = 0; IVect ListeDofBlock;
	while (size < nblock)
	  {
	    new_liste_elem.Clear(); new_pos_elem.Clear();
	    // loop on all elements to obtain the wished size block
	    for (int n1 = 0; n1 < liste_elem.GetM(); n1++)
	      if (size < nblock)
		{
		  int num_elem = liste_elem(n1);
		  int nb_dof_loc = this->mesh_num.GetNbLocalDof(num_elem);
		  LocalDofUsed.Reallocate(nb_dof_loc); LocalDofUsed.Fill(false);
		  int nb_dof_available = 0;
		  // dof numbers that are still free
		  for (int j = 0; j < nb_dof_loc; j++)
		    {
		      num_ddl(j) = this->mesh_num.Element(num_elem).GetNumberDof(j);
		      if (num_ddl(j) >= 0)
			{
			  LocalDofUsed(j) = DofUsed(num_ddl(j));
			  if (!LocalDofUsed(j))
			    nb_dof_available++;
			}
		    }
		  
		  if (nb_dof_available > 0)
		    {
		      nb_dof_available = min(nblock-size, nb_dof_available);
		      // if some dofs are free, we pick them
                      this->PickNearDofs(pos_elem(n1), LocalDofUsed,
                                         ListeLocalDof, nb_dof_available, num_elem);
		      // DISP(ListeLocalDof);
		      
		      // updating the list of dofs and DofUsed
		      ListeDofBlock.Resize(size+nb_dof_available);
		      for (int j = 0; j < nb_dof_available; j++)
			{
			  ListeDofBlock(size+j) = num_ddl(ListeLocalDof(j));
			  DofUsed(ListeDofBlock(size+j)) = true;
			} 
		      size += nb_dof_available;
		      
		      // updating positioning
		      pos_elem(n1) = ListeLocalDof(nb_dof_available-1);
		    }
		}
	    
	    // DISP(liste_elem); DISP(pos_elem); DISP(ListeDofBlock); DISP(size);
	    
	    if (size < nblock)
	      {
		new_liste_elem.Clear();
		// not enough dofs have been found, we look for neighboring elements
		for (int n1 = 0; n1 < liste_elem.GetM(); n1++)
		  {
		    int num_elem = liste_elem(n1);
		    for (int j = 0; j < this->mesh.Element(num_elem).GetNbBoundary(); j++)
		      {
			int num_boundary = this->mesh.Element(num_elem).numBoundary(j);
			if (this->mesh.Boundary(num_boundary).GetNbElements() >= 2)
			  {
			    int num_elem2 = this->mesh.Boundary(num_boundary).numElement(0);
			    if (num_elem2 == num_elem)
			      num_elem2 = this->mesh.Boundary(num_boundary).numElement(1);
			    
			    new_liste_elem.PushBack(num_elem2);
			  }
		      }
		  }
		
		// sorting the list
		Assemble(new_liste_elem);
		
		// adding only unused elements
		liste_elem.Clear();
		for (int n1 = 0; n1 < new_liste_elem.GetM(); n1++)
		  if (!ElementUsed(new_liste_elem(n1)))
		    {
		      liste_elem.PushBack(new_liste_elem(n1));
		      ElementUsed(new_liste_elem(n1)) = true;
		    }
		
		// position set to zero
		pos_elem.Reallocate(liste_elem.GetM()); pos_elem.Fill(0);
	      }
	    
	    if ((size < nblock)&&(liste_elem.GetM() == 0))
	      {
		// no more elements to scan
		size = nblock;
	      }
	    
	  }
	
	// adding the block
	size = ListeDofBlock.GetM()*this->nb_unknowns;
        num_ddl_blocks(nb_blocks).Reallocate(size);
        int nb = 0;
	for (int j = 0; j < ListeDofBlock.GetM(); j++)
	  for (int m = 0; m < this->nb_unknowns; m++)
	    num_ddl_blocks(nb_blocks)(nb++) = ListeDofBlock(j)+m*this->mesh_num.GetNbDof();
	
	nb_blocks++;
        size_all += size;
	// if no more element to scan, we quit
	if (liste_elem.GetM() == 0)
	  size_all = this->nb_unknowns*this->mesh_num.GetNbDof();
      }
    
    num_ddl_blocks.Resize(nb_blocks);
    // DISP(num_ddl_blocks); DISP(size_blocks);
    */
  }
  
  
  //! computation of a block-diagonal submatrix extracted from the iterative matrix
  /*!
    \param[out] diagonal block-diagonal of the iterative matrix
    \param[in] mat iterative matrix (eventually a matrix-free form)
    \param[in] nat_mat mass and stiffness coefficients
  */
  template<class MatrixDiag, class T>
  void VarComputationProblem
  ::ComputeBlockDiagonalMatrix(MatrixDiag& diagonal, int size_block,
			       const FemMatrixFreeClass_Base<T>& mat,
			       const GlobalGenericMatrix<T>& nat_mat)
  {
    Vector<IVect> num_ddl_blocks;
    this->GetProfilBlockDiagonal(size_block, num_ddl_blocks);
    diagonal.SetPattern(num_ddl_blocks);
    
    bool matrix_not_stored = false;
    if (mat.IsSymmetric())
      {
	if (mat.mat_iterative_sym.GetM() > 0)
	  FillBlockDiagonal(mat.mat_iterative_sym, diagonal);
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (mat.mat_iterative_unsym.GetM() > 0)
	  FillBlockDiagonal(mat.mat_iterative_unsym, diagonal);
	else
	  matrix_not_stored = true;
      }

    if (matrix_not_stored)
      {
	AddMatrixWithBC(diagonal, nat_mat);
      }
  }
  
  
  //! computation of a block-diagonal matrix extracted from a standard matrix
  template<class Matrix1, class MatrixDiag, class T>
  void VarComputationProblem
  ::ComputeBlockDiagonalMatrix(MatrixDiag& diagonal, int size_block,
			       const Matrix1& mat,
			       const GlobalGenericMatrix<T>& nat_mat)
  {
    Vector<IVect> num_ddl_blocks;
    this->GetProfilBlockDiagonal(size_block, num_ddl_blocks);
    diagonal.SetPattern(num_ddl_blocks);
   
    FillBlockDiagonal(mat, diagonal);
  }
  
  
  //! copies input parameters of another similar problem
  void VarComputationProblem::CopyInputData(const VarComputationProblem& var)
  {
    this->threshold_matrix = var.threshold_matrix;
    this->storage_finite_element_matrix = var.storage_finite_element_matrix;
    this->leaf_static_condensation = var.leaf_static_condensation;
    this->light_static_condensation = var.light_static_condensation;
    symmetric_elementary_matrix = var.symmetric_elementary_matrix;
    symmetric_global_matrix = var.symmetric_global_matrix;
    use_symmetrization_when_possible = var.use_symmetrization_when_possible;    
  }


  //! returns true if the elementary matrix is symmetric
  bool VarComputationProblem::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const
  {
    return symmetric_elementary_matrix;
  }


  //! returns true if the elementary matrix is symmetric
  bool VarComputationProblem::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const
  {
    return symmetric_elementary_matrix;
  }

  
  //! returns true if the elementary mass matrix is diagonal
  bool VarComputationProblem::IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const
  {
    return false;
  }
  

  //! returns true if the elementary mass matrix is sparse
  bool VarComputationProblem::IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const
  {
    return false;
  }


  //! returns true if the elementary mass matrix is diagonal
  bool VarComputationProblem::IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const
  {
    return sparse_elementary_matrix;
  }
  

  //! returns true if the elementary mass matrix is sparse
  bool VarComputationProblem::IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const
  {    
    return sparse_elementary_matrix;
  }


  //! Adds terms due to numerical fluxes in the finite element matrix
  void VarComputationProblem
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col) 
  {
    if (var_problem.dg_formulation == ElementReference_Base::DISCONTINUOUS)
      {
	cout << "The function AddElementaryFluxesDG should be overloaded in the derived class." << endl
	     << "No generic function is available to compute the numerical fluxes" << endl;
	
	abort();
      }
  }


  //! Adds terms due to numerical fluxes in the finite element matrix
  void VarComputationProblem
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col) 
  {
    if (var_problem.dg_formulation == ElementReference_Base::DISCONTINUOUS)
      {
	cout << "The function AddElementaryFluxesDG should be overloaded in the derived class." << endl
	     << "No generic function is available to compute the numerical fluxes" << endl;
	
	abort();
      }
  }

  
  void VarComputationProblem::UpdateShiftAdimensionalization(Real_wp&, Real_wp&) 
  {
  }
  
  
  void VarComputationProblem::UpdateShiftAdimensionalization(Complex_wp&, Complex_wp&)
  {
  }
  
  
  void VarComputationProblem::UpdateEigenvaluesAdimensionalization(VectReal_wp&, VectReal_wp&,
								   Matrix<Real_wp, General, ColMajor>&)
  {
  }
  
  
  void VarComputationProblem::UpdateEigenvaluesAdimensionalization(VectComplex_wp&, VectComplex_wp&,
								   Matrix<Complex_wp, General, ColMajor>&)
  {
  }


  //! gives which rows are multiplied by -1 when a symmetrization is performed
  void VarComputationProblem::FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const
  {
    i0 = 0; i1 = 0; j0 = 0; j1 = 0;
  }
  
  
  //! modifies the source if the system is symmetrized (usually some rows are multiplied by -1)
  void VarComputationProblem::ModifySourceSymmetry(VectReal_wp& x) const
  {
    if (!var_problem.FirstOrderFormulation())
      return;

    int i0, i1, j0, j1;
    this->FindIntervalDofSignSymmetry(i0, i1, j0, j1);
    
    for (int i = i0; i < i1; i++)
      x(i) = -x(i);
    
    for (int i = j0; i < j1; i++)
      x(i) = -x(i);
  }


  //! modifies the source if the system is symmetrized (usually some rows are multiplied by -1
  void VarComputationProblem::ModifySourceSymmetry(VectComplex_wp& x) const
  {
    if (!var_problem.FirstOrderFormulation())
      return;

    int i0, i1, j0, j1;
    this->FindIntervalDofSignSymmetry(i0, i1, j0, j1);
    
    for (int i = i0; i < i1; i++)
      x(i) = -x(i);
    
    for (int i = j0; i < j1; i++)
      x(i) = -x(i);
  }


  //! modifies the source if the system is symmetrized (usually some rows are multiplied by -1
  void VarComputationProblem::ModifySourceSymmetry(Matrix<Real_wp, General, ColMajor>& A) const
  {
    if (!var_problem.FirstOrderFormulation())
      return;

    int i0, i1, j0, j1;
    this->FindIntervalDofSignSymmetry(i0, i1, j0, j1);
    
    for (int j = 0; j < A.GetN(); j++)
      for (int i = i0; i < i1; i++)
	A(i, j) = -A(i, j);

    for (int j = 0; j < A.GetN(); j++)
      for (int i = j0; i < j1; i++)
	A(i, j) = -A(i, j);
  }


  //! modifies the source if the system is symmetrized (usually some rows are multiplied by -1
  void VarComputationProblem::ModifySourceSymmetry(Matrix<Complex_wp, General, ColMajor>& A) const
  {
    if (!var_problem.FirstOrderFormulation())
      return;
    
    int i0, i1, j0, j1;
    this->FindIntervalDofSignSymmetry(i0, i1, j0, j1);
    
    for (int j = 0; j < A.GetN(); j++)
      for (int i = i0; i < i1; i++)
	A(i, j) = -A(i, j);

    for (int j = 0; j < A.GetN(); j++)
      for (int i = j0; i < j1; i++)
	A(i, j) = -A(i, j);
  }
  
}
  
#define MONTJOIE_FILE_VAR_PROBLEM_BASE_CXX
#endif
