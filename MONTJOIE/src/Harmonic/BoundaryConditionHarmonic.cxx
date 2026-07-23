#ifndef MONTJOIE_FILE_BOUNDARY_CONDITION_HARMONIC_CXX

namespace Montjoie
{

  /*****************************
   * VarBoundaryCondition_Base *
   *****************************/
  
  
  //! setting attributes to default values
  void VarBoundaryCondition_Base::InitDefaultValues()
  {
    nb_dof_dirichlet = 0;
    nb_dof_dirichlet_all = 0;
    
    matrix_symmetric_dir = false;    
    coef_dirichlet_matrix = Real_wp(1);
    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    order_high_conductivity = 0;
#endif

    type_symmetry = NO_SYMMETRY;
    list_number_mode.Reallocate(1);
    list_number_mode(0) = 0;
    plane_wave_quasi_periodic = false;
    number_mode = 0;
    nb_modes_periodic_xyz.Fill(1);
    do_not_store_modes = true;
    
    take_into_account_curvature_for_abc = false;
    gamma_cla_coef = 0.25;
    theta_cla_coef = 0.5;
    zeta_cla_coef = 0.25;

    grazing_abc = false;
   
    order_ABC = 1;

    nb_eltPML = 0;
    nb_eltPML_all = 0;    
    vsigma = 1.0;
    max_velocity_pml = 1.0;

    function_damping_pml = PML_PARABOLE;
    offset_damping_pml = 0.0;
    coef_mixed_damping_pml = 0.0;

    supported_components_BC.Reallocate(PhysicalConstant::nb_max_indices);
    //dispersive_pml = true;
    //omega_pml = Real_wp(1);
  }

    
  //! modification of boundary conditions with a line of the data file
  void VarBoundaryCondition_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("OrderAbsorbingBoundaryCondition"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "OrderAbsorbingBoundaryCondition needs more parameters, for instance :"
                 << endl;
	    cout << "OrderAbsorbingBoundaryCondition = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	order_ABC = atoi(parameters(0).data());
        if (parameters.GetM() > 1)
          {
            if (parameters(1) == "GRAZING")
              grazing_abc = true;
            else
              grazing_abc = false;

            if (parameters(1) == "PARAMETERS")
              {
                if (parameters.GetM() < 5)
                  {
                    cout << "In SetInputData of VarBoundaryCondition" << endl;
                    cout << "OrderAbsorbingBoundaryCondition with PARAMETERS needs more parameters, for instance :"
                         << endl;
                    cout << "OrderAbsorbingBoundaryCondition =  N PARAMETERS  gamma_cla_coef theta_cla_coef zeta_cla_coef" << endl;
                    cout << "Current parameters are : " << endl << parameters << endl;
                    abort();
                  }
                
                gamma_cla_coef = to_num<Real_wp>(parameters(2));
                theta_cla_coef = to_num<Real_wp>(parameters(3));
                zeta_cla_coef = to_num<Real_wp>(parameters(4));
                
                if ( (parameters.GetM() > 5) && (parameters(5) == "GRAZING"))
                  grazing_abc = true;
              }
          }
      }
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    else if (!description_field.compare("OrderHighConductivityBoundaryCondition"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "OrderHighConductivityBoundaryCondition needs more parameters,"
                 << " for instance :" << endl;
	    cout << "OrderHighConductivityBoundaryCondition = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	order_high_conductivity = atoi(parameters(0).data());
      }
#endif
    else if (!description_field.compare("TransverseDampingPML"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "TransverseDampingPML needs more parameters, for instance :" << endl;
	    cout << "TransverseDampingPML = vsigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }
        
        coef_mixed_damping_pml = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("DampingPML"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "DampingPML needs more parameters, for instance :" << endl;
	    cout << "DampingPML = vsigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (parameters(0) == "PARABOLE")
          {
	    if (parameters.GetM() <= 1)
	      {
		cout << "In SetInputData of VarBoundaryCondition" << endl;
		cout << "DampingPML needs more parameters, for instance :" << endl;
		cout << "DampingPML = PARABOLE vsigma" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            function_damping_pml = PML_PARABOLE;
            vsigma = to_num<Real_wp>(parameters(1));
          }
	else if (parameters(0) == "SHIFTED_PARABOLE")
          {
	    if (parameters.GetM() <= 2)
	      {
		cout << "In SetInputData of VarBoundaryCondition" << endl;
		cout << "DampingPML needs more parameters, for instance :" << endl;
		cout << "DampingPML = SHIFTED_PARABOLE vsigma offset" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            function_damping_pml = PML_SHIFTED_PARABOLE;
            vsigma = to_num<Real_wp>(parameters(1));
            offset_damping_pml = to_num<Real_wp>(parameters(2));
          }
	else if (parameters(0) == "LINEAR")
          {
	    if (parameters.GetM() <= 1)
	      {
		cout << "In SetInputData of VarBoundaryCondition" << endl;
		cout << "DampingPML needs more parameters, for instance :" << endl;
		cout << "DampingPML = LINEAR vsigma" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }
            
            function_damping_pml = PML_LINEAR;
            vsigma = to_num<Real_wp>(parameters(1));
          }
	else if (parameters(0) == "CONSTANT")
          {
	    if (parameters.GetM() <= 1)
	      {
		cout << "In SetInputData of VarBoundaryCondition" << endl;
		cout << "DampingPML needs more parameters, for instance :" << endl;
		cout << "DampingPML = CONSTANT vsigma" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            function_damping_pml = PML_CONSTANT;
            vsigma = to_num<Real_wp>(parameters(1));
          }
        else
          vsigma = to_num<Real_wp>(parameters(0));
      }
    /*    else if (!description_field.compare("DispersivePML"))
      {
	if (parameters(0) == "NO")
	  {
	    dispersive_pml = false;
	    this->omega_pml = 2.0*pi_wp*to_num<Real_wp>(parameters(1));
	  }
	else
	  dispersive_pml = true;
          }*/
    else if (!description_field.compare("ConditionReference"))
      {
	int position_word = -1;
	// we look for the first parameter, with begins by a letter
	// This parameter is the description of the boundary condition (DIRICHLET, NEUMANN, etc)
	for (int i = 0; i < parameters.GetM(); i++)
	  if (isalpha(parameters(i)[0]))
	    {
	      position_word = i;
	      break;
	    }
	
	if (position_word > 0)
	  {
	    // All the parameters at the left of the description of the boundary condition,
	    // are references. The references are stored
	    IVect ref(position_word);
	    for (int i = 0; i < position_word; i++)
	      {
		ref(i) = to_num<int>(parameters(i));
		
		this->CheckAndReallocateReferences(ref(i));	       				
	      }
	    
	    // if some values appear after the description of the boundary condition,
	    // we consider that these values are parameters of the BC
	    if (parameters.GetM() > (position_word+1))
	      this->FillParameterCondition(parameters, position_word, ref);
	    
	    bool periodic_condition = false;
	    int boundary_condition
	      = GetBoundaryConditionId(ref, position_word, parameters, periodic_condition);
            
	    if (periodic_condition)
	      {
		// We inform the mesh that a periodic condition is set
		if (position_word != 2)
		  {
		    cout << "Sorry we need two reference numbers for periodicity condition " << endl;
		    abort();
		  }
		
		this->AddPeriodicConditionMesh(TinyVector<int, 2>(ref(0), ref(1)), boundary_condition);		
	      }
	    else
	      {
		for (int i = 0; i < position_word; i++)
		  this->SetBoundaryConditionMesh(ref(i), boundary_condition);
	      }

            if (boundary_condition == BoundaryConditionEnum::LINE_SUPPORTED)
              {
                for (int j = 0; j < ref.GetM(); j++)
                  this->supported_components_BC(ref(j)).Clear();
                
                for (int i = (position_word+1); i < parameters.GetM(); i++)
                  {
                    int num = to_num<int>(parameters(i));
                    for (int j = 0; j < ref.GetM(); j++)
                      this->supported_components_BC(ref(j)).PushBack(num);                      
                  } 
              }
	  }
      } 
    else if (description_field == "ForceDirichletSymmetry")
      {
       	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "StorageModes needs more parameters, for instance :" << endl;
	    cout << "ForceDirichletSymmetry = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (parameters(0) == "YES")
          matrix_symmetric_dir = true;
        else
          matrix_symmetric_dir = false; 
      }
    else if (!description_field.compare("DirichletCoefMatrix"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition_Base" << endl;
	    cout << "DirichletCoefMatrix needs more parameters, for instance :" << endl;
	    cout << "DirichletCoefMatrix = coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
        coef_dirichlet_matrix = to_num<Real_wp>(parameters(0));
      }
    else if (!description_field.compare("StorageModes"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "StorageModes needs more parameters, for instance :" << endl;
	    cout << "StorageModes = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          do_not_store_modes = false;
        else
          do_not_store_modes = true;
      }
    else if (description_field == "PlaneWaveQuasiPeriodic")
      {
        if (parameters(0) == "YES")
          plane_wave_quasi_periodic = true;
        else
          plane_wave_quasi_periodic = false;
      }
  }


  //! returns the type of symmetrization to use when calling PerformFactorizationStep
  bool VarBoundaryCondition_Base::GetInitialSymmetrization() const
  {
    bool use_symmetrization = true;
    if (this->GetNbGlobalEltPML() > 0)
      use_symmetrization = false;

    return use_symmetrization;
  }
  
  
  //! returns the boundary condition id number associated with a keyword
  int VarBoundaryCondition_Base
  ::GetBoundaryConditionId(const IVect& ref, int pos, const VectString& parameters, bool& periodic)
  {
    if (pos >= parameters.GetM() )
      {
	cout << "Problem while reading the boundary condition" << endl;
	DISP(pos); DISP(parameters);
	abort();
      }
    
    string keyword = parameters(pos);
    periodic = false;
    
    // matching between description of BC and number of BC
    // we could use a map ...
    if (keyword == "DIRICHLET")
      return BoundaryConditionEnum::LINE_DIRICHLET;
    else if (keyword == "SUPPORTED")
      return BoundaryConditionEnum::LINE_SUPPORTED;
    else if (keyword == "NEUMANN")
      return BoundaryConditionEnum::LINE_NEUMANN;
    else if (keyword == "ABSORBING")
      return BoundaryConditionEnum::LINE_ABSORBING;
    else if (keyword == "DTN")
      return BoundaryConditionEnum::LINE_DTN;
    else if (keyword == "HIGH_CONDUCTIVITY")
      return BoundaryConditionEnum::LINE_HIGH_CONDUCTIVITY;
    else if (keyword == "IMPEDANCE")
      return BoundaryConditionEnum::LINE_IMPEDANCE;
    else if (keyword == "PERIODICITY")
      {
	periodic = true;
	return BoundaryConditionEnum::PERIODIC_CTE;
      }
    else if (keyword == "PERIODIC_X")
      {
	periodic = true;
	return BoundaryConditionEnum::PERIODIC_X;
      }
    else if (keyword == "PERIODIC_Y")
      {
	periodic = true;
	return BoundaryConditionEnum::PERIODIC_Y;
      }
    else if (keyword == "PERIODIC_Z")
      {
	periodic = true;
	return BoundaryConditionEnum::PERIODIC_Z;
      }
    else if (keyword == "CYCLIC")
      {
	periodic = true;
	return BoundaryConditionEnum::PERIODIC_THETA;
      }
    
    return BoundaryConditionEnum::LINE_INSIDE;
  }
  

  //! retrieves dof numbers of dofs associated with a given reference
  template<class Dimension>
  void VarBoundaryCondition_Base
  ::FindDofsOnReference(const VarProblem<Dimension>& var,
			const Vector<int>& ref_cond, int ref_target, Vector<int>& Dofs)
  {
    // similar function used for Dirichlet/Supported dofs
    Vector<bool> DofTaken(var.offset_dof_unknown(1));
    DofTaken.Fill(false);
    const Mesh<Dimension>& mesh = var.mesh;    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (ref_cond(ref) == ref_target)
          {
            int num_face = i;
            int ne = mesh.Boundary(num_face).numElement(0);
            if (ne >= 0)
              {
                int num_loc = mesh.Element(ne).GetPositionBoundary(num_face);
		const ElementReference_Dim<Dimension>& Fb = var.GetReferenceElement(ne);
		IVect Nodle = var.GetDofNumberOnElement(ne);
		int nb_dof = Fb.GetNbDofBoundary(num_loc);
		for (int j = 0; j < nb_dof; j++)
		  {
		    int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
		    int num_dof = Nodle(num_dof_loc);
		    if (num_dof >= 0)
		      DofTaken(num_dof) = true;
		  }
	      }
	  }
      }

#ifdef SELDON_WITH_MPI
    // assembling DofTaken
    var.ReduceDistributedVector(DofTaken, MPI_MAX, 1);
#endif
    
    // dofs are counted
    int nb_dof = 0;
    for (int i = 0; i < DofTaken.GetM(); i++)
      if (DofTaken(i))
	nb_dof++;
    
    Dofs.Reallocate(nb_dof);
    nb_dof = 0;
    for (int i = 0; i < DofTaken.GetM(); i++)
      if (DofTaken(i))
	Dofs(nb_dof++) = i;    
  }
  
  
  /////////////////////////
  // DIRICHLET CONDITION //
  
  
  //! internal function used to retrieve Dirichlet dofs on a boundary
  template<class Dimension>
  void VarBoundaryCondition_Base
  ::FindLocalDirichletDof(int i, int num_face, int ne, int num_loc, int ref_cond,
			  const VarProblem<Dimension>& var,
			  Vector<int>& IndexDirichlet,
			  const ElementReference_Dim<Dimension>& Fb, int num)
  {
    IVect Nodle = var.GetDofNumberOnElement(ne, num);
    int nb_dof = Fb.GetNbDofBoundary(num_loc);
    for (int j = 0; j < nb_dof; j++)
      { 
        int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
        int num_dof = Nodle(num_dof_loc);
        if (num_dof >= 0)
	  IndexDirichlet(num_dof) = ref_cond;
      }
  }
  
  
  //! treatment for dirichlet condition for any problem
  /*!
    Dirichlet_dof will contain the list of dofs associated
    to a Dirichlet condition
    nb_dof_dirichlet will be the number of Dirichlet dofs
  */
  template<class Dimension>
  void VarBoundaryCondition_Base
  ::TreatDirichletCondition(const VarProblem<Dimension>& var)
  {
    if ( (var.FormulationDG() != ElementReference_Base::CONTINUOUS) && var.FirstOrderFormulationDG())
      {
        // if local discontinuous galerkin or HDG => no dof with dirichlet
	nb_dof_dirichlet = 0;
	Dirichlet_dof.Clear();
	return;
      }
    
    const Mesh<Dimension>& mesh = var.mesh;
    
    // first stage : constructing array IndexDirichlet
    int nb_mesh_num = var.GetNbMeshNumberings();
    Vector<Vector<int> > IndexDirichlet(nb_mesh_num), IndexDirichletRef(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      {
        int n0 = -1;
        for (int j = 0; j < var_problem.nb_unknowns; j++)
          if (var_problem.mesh_num_unknown(j) == i)
            n0 = j;

        int N = var.GetMeshNumbering(i).GetNbDof();
        if (n0 >= 0)
          N = var.offset_dof_unknown(n0+1) - var.offset_dof_unknown(n0);
        
        IndexDirichlet(i).Reallocate(N);        
        IndexDirichlet(i).Fill(-1);

        IndexDirichletRef(i).Reallocate(N);        
        IndexDirichletRef(i).Fill(-1);
      }
    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_DIRICHLET
            || mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_SUPPORTED) 
          {
	    int ref_cond = mesh.GetBoundaryCondition(ref);
            int num_face = i;
            int ne = mesh.Boundary(num_face).numElement(0);
            if (ne >= 0)
              {
                int num_loc = mesh.Element(ne).GetPositionBoundary(num_face);
                for (int nm = 0; nm < nb_mesh_num; nm++)
                  {
                    FindLocalDirichletDof(i, num_face, ne, num_loc, ref_cond,
                                          var, IndexDirichlet(nm), var.GetReferenceElement(ne, nm), nm);

                    FindLocalDirichletDof(i, num_face, ne, num_loc, ref,
                                          var, IndexDirichletRef(nm), var.GetReferenceElement(ne, nm), nm);
                  }
              }
          }                    
      }
    
    AddDirichletPointConditions(IndexDirichlet);
          
#ifdef SELDON_WITH_MPI
    // assembling IndexDirichlet in order to retrieve Dirichlet dofs of neighboring domains
    for (int k = 0; k < IndexDirichlet.GetM(); k++)
      {
        var.ReduceDistributedVector(IndexDirichlet(k), MPI_MAX, k, true);
        var.ReduceDistributedVector(IndexDirichletRef(k), MPI_MAX, k, true);
      }
#endif
    
    // IsDofDirichlet similar to IndexDirichlet but for all unknowns
    is_dof_dirichlet.Reallocate(var.GetNbDof());
    is_dof_dirichlet.Fill(false);
    nb_dof_dirichlet = 0;
    int nb_dof_scalar = var.offset_dof_unknown(1);
    if (nb_mesh_num == 1)
      for (int i = 0; i < IndexDirichlet(0).GetM(); i++)
        {
          if (IndexDirichlet(0)(i) == BoundaryConditionEnum::LINE_SUPPORTED)
            {
              // supported dof
              int ref = IndexDirichletRef(0)(i);
              int nb_u = this->GetNbSupportedComponents(ref);
              for (int p1 = 0; p1 < nb_u; p1++)
                {
                  int p = this->GetSupportedComponent(ref, p1);
                  is_dof_dirichlet(i + p*nb_dof_scalar) = true;
                }
              
              nb_dof_dirichlet += nb_u;
            }
          else if (IndexDirichlet(0)(i) == BoundaryConditionEnum::LINE_DIRICHLET)
            {
              // Dirichlet dof
              for (int p = 0; p < var.nb_unknowns; p++)
                is_dof_dirichlet(i + p*nb_dof_scalar) = true;
              
              nb_dof_dirichlet += var.nb_unknowns;
            }
        }
    else
      {
        for (int num_u = 0; num_u < var.nb_unknowns_scal; num_u++)
          {
            int num = var.mesh_num_unknown(num_u);
            int offset = var.offset_dof_unknown(num_u);
            int Ndof = var.GetMeshNumbering(num).GetNbDof();
                        
            for (int i = 0; i < Ndof; i++)
              {
                int ref = IndexDirichletRef(num)(i);
                if (ref < 0)
                  continue;
                
                int nb_supp = this->GetNbSupportedComponents(ref);        
                bool supp = false;
                for (int p1 = 0; p1 < nb_supp; p1++)
                  if (this->GetSupportedComponent(ref, p1) == num_u)
                    supp = true;

                int p = IndexDirichlet(num)(i);
                bool dir_dof = false;
                if (p == BoundaryConditionEnum::LINE_DIRICHLET)
                  dir_dof = true;
                
                if ((p == BoundaryConditionEnum::LINE_SUPPORTED) && supp)
                  dir_dof = true;

                if (dir_dof)
                  {
                    is_dof_dirichlet(offset + i) = true;
                    nb_dof_dirichlet++;
                  }
              }
          }
      }
      
    // second stage : constructing array Dirichlet_dof
    if (nb_dof_dirichlet > 0)
      {
        this->Dirichlet_dof.Reallocate(nb_dof_dirichlet);
        nb_dof_dirichlet = 0;
        for (int i = 0; i < is_dof_dirichlet.GetM(); i++)
          if (is_dof_dirichlet(i))
            {
              this->Dirichlet_dof(nb_dof_dirichlet) = i;
              nb_dof_dirichlet++;
            }
      }

    nb_dof_dirichlet_all = nb_dof_dirichlet;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    MPI_Allreduce(&nb_dof_dirichlet, &nb_dof_dirichlet_all, 1, 
                  MPI_INTEGER, MPI_SUM, var.comm_group_mode);
#else
    int rank_proc(0);
#endif
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 4)
        cout << rank_proc << " Dirichlet done" << endl;
  }
  

  //! modifying dofs for which Dirichlet condition is applied
  void VarBoundaryCondition_Base
  ::SetDirichletDofs(int N, const IVect& num_dof)
  {
    if (N <= 0)
      {
        Dirichlet_dof.Clear();
        is_dof_dirichlet.Fill(false);
        nb_dof_dirichlet = 0;
        return;
      }
    
    nb_dof_dirichlet = N;
    Dirichlet_dof.Reallocate(N);
    for (int i = 0; i < N; i++)
      Dirichlet_dof(i) = num_dof(i);
    
    is_dof_dirichlet.Fill(false);
    for (int i = 0; i < N; i++)
      is_dof_dirichlet(num_dof(i)) = true;

    nb_dof_dirichlet_all = nb_dof_dirichlet;
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    MPI_Allreduce(&nb_dof_dirichlet, &nb_dof_dirichlet_all, 1, 
                  MPI_INTEGER, MPI_SUM, var_comm.comm_group_mode);
#else
    int rank_proc(0);
#endif

  }
  

  //! modifies the size of the problem
  void VarBoundaryCondition_Base::ResizeNbDof(int n)
  {
    if (this->is_dof_dirichlet.GetM() < n)
      {
        int ndof = this->is_dof_dirichlet.GetM();
        this->is_dof_dirichlet.Resize(n);
        for (int j = ndof; j < n; j++)
          this->is_dof_dirichlet(j) = false;
      }
    
    var_problem.SetNbDof(n);
  }


  //! estimation of higher eigenvalues and modification of Dirichlet coefficient
  /*!
    \param[in,out] mat_direct_real sparse matrix to be modified
    This methode estimates the higher eigenvalue of the sparse matrix.
    The dirichlet coefficient is set to twice this estimation
  */
  template<class T>
  Real_wp VarBoundaryCondition_Base::ComputeDirichletCoef(VirtualMatrix<T>& A)
  {
    int nodl = A.GetM();
    // estimation of the maximal eigenvalue if we want small eigenvalues,
    // to reject dirichlet condition at infinity
    Vector<T> X(nodl), Y(nodl);
    X.FillRand();
    X *= 1e-9;
    
    // on itère X = A*X
    for (int i = 0; i < 5; i++)
      {
        Mlt(1.0/Norm2(X), X);
	Seldon::Copy(X, Y);
	A.MltVector(X, Y);
      }

    Real_wp coef = 2.0*Norm2(X)/Norm2(Y);
    return coef;
  }

  
  //! modification of right hand side to take into account inhomogeneous Dirichlet condition
  void VarBoundaryCondition_Base
  ::ApplyDirichletCondition(const SeldonTranspose& trans, FemMatrixFreeClass_Base<Real_wp>& mat_sp,
			    Vector<Real_wp>& b_rhs, int k) const
  {
    if (mat_sp.FormulationDG() == ElementReference_Base::CONTINUOUS)
      this->ApplyPeriodicCondition(b_rhs);
    
    mat_sp.ApplyDirichletCondition(trans, b_rhs, k);
  }  


  //! modification of right hand side to take into account inhomogeneous Dirichlet condition
  void VarBoundaryCondition_Base
  ::ApplyDirichletCondition(const SeldonTranspose& trans, FemMatrixFreeClass_Base<Real_wp>& mat_sp,
			    Vector<Complex_wp>& b_rhs, int k) const
  {
    if (mat_sp.FormulationDG() == ElementReference_Base::CONTINUOUS)
      this->ApplyPeriodicCondition(b_rhs);
    
    mat_sp.ApplyDirichletCondition(trans, b_rhs, k);
  }  
  

  //! modification of right hand side to take into account inhomogeneous Dirichlet condition
  void VarBoundaryCondition_Base
  ::ApplyDirichletCondition(const SeldonTranspose& trans, FemMatrixFreeClass_Base<Complex_wp>& mat_sp,
			    Vector<Real_wp>& b_rhs, int k) const
  {
    cout << "Forbidden case" << endl;
    abort();
  }


  //! modification of right hand side to take into account inhomogeneous Dirichlet condition
  void VarBoundaryCondition_Base
  ::ApplyDirichletCondition(const SeldonTranspose& trans, FemMatrixFreeClass_Base<Complex_wp>& mat_sp,
			    Vector<Complex_wp>& b_rhs, int k) const
  {
    if (mat_sp.FormulationDG() == ElementReference_Base::CONTINUOUS)
      this->ApplyPeriodicCondition(b_rhs);
    
    mat_sp.ApplyDirichletCondition(trans, b_rhs, k);
  }  

  
  //! Updates Dirichlet dofs
  void VarBoundaryCondition_Base::UpdateDirichletDofs()
  {
    this->nb_dof_dirichlet = 0;
    for (int i = 0; i < this->is_dof_dirichlet.GetM(); i++)
      if (this->is_dof_dirichlet(i))
        this->nb_dof_dirichlet++;
    
    this->Dirichlet_dof.Reallocate(this->nb_dof_dirichlet);
    this->nb_dof_dirichlet = 0;
    for (int i = 0; i < this->is_dof_dirichlet.GetM(); i++)
      if (this->is_dof_dirichlet(i))
        this->Dirichlet_dof(this->nb_dof_dirichlet++) = i;    
    
    this->nb_dof_dirichlet_all = this->nb_dof_dirichlet;
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&this->nb_dof_dirichlet, &this->nb_dof_dirichlet_all, 1, 
                  MPI_INTEGER, MPI_SUM, var_comm.comm_group_mode);
#endif
  }

  
  // DIRICHLET CONDITION //
  /////////////////////////
  
  
  ////////////////////////
  // PERIODIC CONDITION //
  
  
  //! gives the number of modes along each direction
  /*!
    \param[out] nx number of modes in x
    \param[out] ny number of modes in y
    \param[out] nz number of modes in z
    \param[out] teta_sym if true the domain is also cyclic
  */
  void VarBoundaryCondition_Base
  ::GetPeriodicNumberModes(int& nx, int& ny, int& nz, bool& teta_sym) const
  {  
    teta_sym = false;
    nx = 1; ny = 1; nz = 1;
    if (plane_wave_quasi_periodic)
      return;
    
    switch (type_symmetry)
      {
      case PERIODIC_THETA :
        nx = nb_modes_periodic_xyz(0);
        teta_sym = true;
        break;
      case PERIODIC_X :
        nx = nb_modes_periodic_xyz(0);
        break;
      case PERIODIC_Y :
        nx = nb_modes_periodic_xyz(1);
        break;
      case PERIODIC_Z :
        nx = nb_modes_periodic_xyz(2);
        break;
      case PERIODIC_XY :
        nx = nb_modes_periodic_xyz(0);
        ny = nb_modes_periodic_xyz(1);
        break;
      case PERIODIC_XZ :
        nx = nb_modes_periodic_xyz(0);
        ny = nb_modes_periodic_xyz(2);
        break;
      case PERIODIC_YZ :
        nx = nb_modes_periodic_xyz(1);
        ny = nb_modes_periodic_xyz(2);
        break;
      case PERIODIC_XYZ :
        nx = nb_modes_periodic_xyz(0);
        ny = nb_modes_periodic_xyz(1);
        nz = nb_modes_periodic_xyz(2);
        break;
      case PERIODIC_ZTHETA :
        nx = nb_modes_periodic_xyz(0);
        ny = nb_modes_periodic_xyz(2);
        teta_sym = true;
        break;
      }
  }


  //! gives the mode numbers along each direction
  /*!
    \param[in] n mode number
    \param[out] ix mode number in x
    \param[out] iy mode number in y
    \param[out] iz mode number in z
    \param[out] teta_sym if true the domain is also cyclic
  */  
  void VarBoundaryCondition_Base
  ::GetPeriodicModes(int n, int& ix, int& iy, int& iz, bool& teta_sym) const
  {  
    teta_sym = false;
    ix = -1; iy = -1; iz = -1;
    if (plane_wave_quasi_periodic)
      return;
    
    switch (type_symmetry)
      {
      case PERIODIC_THETA :
        ix = n;
        teta_sym = true;
        break;
      case PERIODIC_X :
      case PERIODIC_Y :
      case PERIODIC_Z :
        ix = n;
        break;
      case PERIODIC_XY :
        ix = n/nb_modes_periodic_xyz(1);
        iy = n%nb_modes_periodic_xyz(1);
        break;
      case PERIODIC_XZ :
      case PERIODIC_YZ :
        ix = n/nb_modes_periodic_xyz(2);
        iy = n%nb_modes_periodic_xyz(2);
        break;
      case PERIODIC_XYZ :
        ix = n/(nb_modes_periodic_xyz(2)*nb_modes_periodic_xyz(1));
        iy = (n - ix*nb_modes_periodic_xyz(2)*nb_modes_periodic_xyz(1))/nb_modes_periodic_xyz(1);
        iz = n%nb_modes_periodic_xyz(1);
        break;
      case PERIODIC_ZTHETA :
        ix = n/nb_modes_periodic_xyz(2);
        iy = n%nb_modes_periodic_xyz(2);
        teta_sym = true;
        break;
      }
  }
  
  
  //! gives the phases along each direction
  /*!
    \param[in] n mode number
    \param[out] kx phase in x
    \param[out] ky phase in y
    \param[out] kz phase in z
  */  
  void VarBoundaryCondition_Base
  ::GetPeriodicModes(int n, Complex_wp& kx, Complex_wp& ky, Complex_wp& kz) const
  {
    kx = 0; ky = 0; kz = 0;
    if (plane_wave_quasi_periodic)
      {
        if (var_source.GetIncidentFieldType(0) == VarSourceProblem_Base::INCIDENT_PLANE_WAVE_CPLX)
          {
            const Vector<string>& param = var_source.GetParameterSource(0, 0);
            Real_wp omega = var_problem.GetOmega();
            for (int k = 0; k < param.GetM(); k++)
              {
                if (param(k) == "Kx")
                  kx = to_num<Complex_wp>(param(k+1));
                else if (param(k) == "Ky")
                  ky = to_num<Complex_wp>(param(k+1));
              }
            
            if (var_problem.GetDimension() == 2)
              ky = sqrt(omega*omega - kx*kx);
            else
              kz = sqrt(omega*omega - kx*kx - ky*ky);
          }
        else
          {
            Real_wp kx_r, ky_r, kz_r;
            var_problem.FillWaveVectorComponents(kx_r, ky_r, kz_r);
            kx = kx_r; ky = ky_r; kz = kz_r;
          }

        return;
      }
    
    switch (type_symmetry)
      {
      case PERIODIC_THETA :
        kx = 2.0*pi_wp*Real_wp(n)/nb_modes_periodic_xyz(0);
        break;
      case PERIODIC_X :
        kx = 2.0*pi_wp*Real_wp(n)/nb_modes_periodic_xyz(0);
        break;
      case PERIODIC_Y :
        ky = 2.0*pi_wp*Real_wp(n)/nb_modes_periodic_xyz(1);
        break;
      case PERIODIC_Z :
        kz = 2.0*pi_wp*Real_wp(n)/nb_modes_periodic_xyz(2);
        break;
      case PERIODIC_XY :
        {
          int ix = n/nb_modes_periodic_xyz(1);
          int iy = n%nb_modes_periodic_xyz(1);
          kx = 2.0*pi_wp*Real_wp(ix)/nb_modes_periodic_xyz(0);
          ky = 2.0*pi_wp*Real_wp(iy)/nb_modes_periodic_xyz(1);
        }
        break;
      case PERIODIC_XZ :
        {
          int ix = n/nb_modes_periodic_xyz(2);
          int iy = n%nb_modes_periodic_xyz(2);
          kx = 2.0*pi_wp*Real_wp(ix)/nb_modes_periodic_xyz(0);
          kz = 2.0*pi_wp*Real_wp(iy)/nb_modes_periodic_xyz(2);
        }
        break;
      case PERIODIC_YZ :
        {
          int ix = n/nb_modes_periodic_xyz(2);
          int iy = n%nb_modes_periodic_xyz(2);
          ky = 2.0*pi_wp*Real_wp(ix)/nb_modes_periodic_xyz(1);
          kz = 2.0*pi_wp*Real_wp(iy)/nb_modes_periodic_xyz(2);
        }
        break;
      case PERIODIC_XYZ :
        {
          int ix = n/(nb_modes_periodic_xyz(2)*nb_modes_periodic_xyz(1));
          int iy = (n - ix*nb_modes_periodic_xyz(2)*nb_modes_periodic_xyz(1))
            /nb_modes_periodic_xyz(1);
          
          int iz = n%nb_modes_periodic_xyz(1);
          kx = 2.0*pi_wp*Real_wp(ix)/nb_modes_periodic_xyz(0);
          ky = 2.0*pi_wp*Real_wp(iy)/nb_modes_periodic_xyz(1);
          kz = 2.0*pi_wp*Real_wp(iz)/nb_modes_periodic_xyz(2);
        }
        break;
      case PERIODIC_ZTHETA :
        {          
          int ix = n/nb_modes_periodic_xyz(2);
          int iy = n%nb_modes_periodic_xyz(2);
          kx = 2.0*pi_wp*Real_wp(ix)/nb_modes_periodic_xyz(0);
          kz = 2.0*pi_wp*Real_wp(iy)/nb_modes_periodic_xyz(2);
        }
        break;
      }
  }
  

  //! modifies matrix in order to take into account quasi-periodic conditions  
  template<class T, class Storage>
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(Matrix<T, Symmetric, Storage>&)
  {
    if (this->GetNbPeriodicDof() > 0)
      {
        cout << "not implemented "  << endl;
	abort();
      }
  }
    

  //! modifies matrix in order to take into account quasi-periodic conditions  
  template<class T>
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(Matrix<T, Symmetric, DiagonalRow>&)
  {
    if (this->GetNbPeriodicDof() > 0)
      {
        cout << "not implemented "  << endl;
	abort();
      }
  }
  
  
  //! modifies matrix in order to take into account quasi-periodic conditions
  template<class T>
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(Matrix<T, General, BlockDiagRow>&)
  {    
    if (this->GetNbPeriodicDof() > 0)
      {
        cout << "not implemented "  << endl;
	abort();
      }
  }
  
  
  //! modifies matrix in order to take into account quasi-periodic conditions
  template<class T>
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(FemMatrixFreeClass_Base<T>&)
  {    
    // to do
    if (this->GetNbPeriodicDof() > 0)
      {
        cout << "not implemented "  << endl;
	abort();
      }
  }
  
  
  //! modifies matrix in order to take into account quasi-periodic conditions
  template<class T>
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(Matrix<T, General, ArrayRowSparse>& mat_sp)
  {
    int nb_periodic_dof = this->GetNbPeriodicDof();
    if ((nb_periodic_dof <= 0)||(var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS))
      return;
    
    T value, value_conj; int nb = 0;
    IVect col_interac; Vector<T> val_interac;
    int nb_unknowns = var_problem.nb_unknowns;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nb_unknowns = var_problem.nb_unknowns_hdg;
    
    int Nvol = var_problem.offset_dof_unknown(1);
    // DISP(this->mesh.kwave_periodicity); 
    for (int i = 0; i < nb_periodic_dof; i++)
      {
	int num_per_, num_orig_;
	this->GetPeriodicDofNumbers(i, num_per_, num_orig_);
	
	// lines of periodic dofs are added to original dofs	
	this->GetPeriodicPhase(i, value);
        value_conj = Real_wp(1) / value;
        
        if ((var_problem.FormulationDG() == ElementReference_Base::HDG)
	    || (!this->IsDofDirichlet(num_per_)))
          for (int m = 0; m < nb_unknowns; m++)
            {
              int num_per = num_per_ + m*Nvol;
              int num_orig = num_orig_ + m*Nvol;
              
              // number of elements on the row
              nb = mat_sp.GetRowSize(num_per);
              // (col_interac,val_interac) -> row sparsely stored 
              col_interac.Reallocate(nb); val_interac.Reallocate(nb);
              for (int p = 0; p < nb; p++)
                {
                  col_interac(p) = mat_sp.Index(num_per, p);
                  val_interac(p) = value_conj*mat_sp.Value(num_per, p);
                }
              
              // we add to the line L_{num_orig), the line L_{num_per}
              mat_sp.AddInteractionRow(num_orig, nb, col_interac, val_interac);
              
              // we set a new equation u_{num_per} = phase*u_{num_orig}
              mat_sp.ReallocateRow(num_per, 2); 
              if (num_per > num_orig)
                {  
                  mat_sp.Index(num_per, 0) = num_orig; 
                  mat_sp.Value(num_per, 0) = -value;		  
                  mat_sp.Index(num_per, 1) = num_per; 
                  mat_sp.Value(num_per, 1) = 1.0;
                }
              else
                {
                  mat_sp.Index(num_per, 1) = num_orig; 
                  mat_sp.Value(num_per, 1) = -value;
		  mat_sp.Index(num_per, 0) = num_per; 
                  mat_sp.Value(num_per, 0) = 1.0;
                }
            }
      }
  }
  
  
  //! modifies matrix in order to take into account quasi-periodic conditions
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(Matrix<Complex_wp, General, ArrayRowComplexSparse>& mat_sp)
  {
    int nb_periodic_dof = this->GetNbPeriodicDof();
    if ((nb_periodic_dof <= 0)||(var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS))
      return;
    
    Complex_wp value, value_conj; int nb = 0;
    IVect col_interac; Vector<Complex_wp> val_interac;
    int nb_unknowns = var_problem.nb_unknowns;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nb_unknowns = var_problem.nb_unknowns_hdg;

    int Nvol = var_problem.offset_dof_unknown(1);
    // DISP(this->mesh.kwave_periodicity); 
    for (int i = 0; i < nb_periodic_dof; i++)
      {
	int num_per_, num_orig_;
	this->GetPeriodicDofNumbers(i, num_per_, num_orig_);
	
	// lines of periodic dofs are added to original dofs	
	this->GetPeriodicPhase(i, value);
        value_conj = Real_wp(1) / value;

	if ((var_problem.FormulationDG() == ElementReference_Base::HDG)
	    || (!this->IsDofDirichlet(num_per_)))
          for (int m = 0; m < nb_unknowns; m++)
            {
              int num_per = num_per_ + m*Nvol;
              int num_orig = num_orig_ + m*Nvol;
              
              // number of elements on the row
              nb = mat_sp.GetRealRowSize(num_per);
              int ni = mat_sp.GetImagRowSize(num_per);
              // (col_interac,val_interac) -> row sparsely stored 
              col_interac.Reallocate(nb+ni); val_interac.Reallocate(nb+ni);
              for (int p = 0; p < nb; p++)
                {
                  col_interac(p) = mat_sp.IndexReal(num_per, p);
                  val_interac(p) = value_conj*mat_sp.ValueReal(num_per, p);
                }
              
              for (int p = 0; p < ni; p++)
                {
                  col_interac(nb+p) = mat_sp.IndexImag(num_per, p);
                  val_interac(nb+p) = Iwp*value_conj*mat_sp.ValueImag(num_per, p);
                }
              
              // we add to the line L_{num_orig), the line phase*L_{num_per}
              mat_sp.AddInteractionRow(num_orig, nb, col_interac, val_interac);
              
              // we set a new equation u_{num_per} = phase*u_{num_orig}
              mat_sp.ReallocateRealRow(num_per, 2);
              mat_sp.ReallocateImagRow(num_per, 1); 
              if (num_per > num_orig)
                {
                  mat_sp.IndexReal(num_per, 0) = num_orig; 
                  mat_sp.ValueReal(num_per, 0) = -real(value);
                  mat_sp.IndexReal(num_per, 1) = num_per; 
                  mat_sp.ValueReal(num_per, 1) = 1.0;
                }
              else
                {
                  mat_sp.IndexReal(num_per, 1) = num_orig; 
                  mat_sp.ValueReal(num_per, 1) = -real(value);
                  mat_sp.IndexReal(num_per, 0) = num_per; 
                  mat_sp.ValueReal(num_per, 0) = 1.0;
                }
              
              mat_sp.IndexImag(num_per, 0) = num_orig; 
              mat_sp.ValueImag(num_per, 0) = -imag(value);
              
            }
      }
  }
  

#ifdef SELDON_WITH_MPI  
  //! modifies matrix in order to take into account quasi-periodic conditions
  template<class T>
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(DistributedMatrix<T, General, ArrayRowSparse>& mat_sp)
  {
    int nb_periodic_dof = this->GetNbPeriodicDof();
    if ((nb_periodic_dof <= 0)||(var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS))
      return;
        
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    T value, value_conj; int nb = 0;
    IVect col_interac; Vector<T> val_interac;
    int nb_unknowns = var_problem.nb_unknowns;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nb_unknowns = var_problem.nb_unknowns_hdg;
    
    int Nvol = var_problem.offset_dof_unknown(1);
    //int Nvol_glob = var_comm.GetNbGlobalScalarDof();
    // DISP(this->mesh.kwave_periodicity); 
    for (int i = 0; i < nb_periodic_dof; i++)
      {
	int num_per_, num_orig_;
	this->GetPeriodicDofNumbers(i, num_per_, num_orig_);
        
	// lines of periodic dofs are added to original dofs	
	this->GetPeriodicPhase(i, value);
        value_conj = Real_wp(1) / value;

	if ((var_problem.FormulationDG() == ElementReference_Base::HDG)
	    || (!this->IsDofDirichlet(num_per_)))
          for (int m = 0; m < nb_unknowns; m++)
            {
              int num_per = num_per_ + m*Nvol;
              
              // number of elements on the row
              nb = mat_sp.GetRowSize(num_per);
              // (col_interac,val_interac) -> row sparsely stored 
              col_interac.Reallocate(nb);
              val_interac.Reallocate(nb);
              for (int p = 0; p < nb; p++)
                {
                  col_interac(p) = mat_sp.Index(num_per, p);
                  val_interac(p) = value_conj*mat_sp.Value(num_per, p);
                }
              
              // we add to the line L_{num_orig), the line L_{num_per}
	      int num_orig = num_orig_ + m*Nvol;
	      mat_sp.AddInteractionRow(num_orig, nb, col_interac, val_interac);
              
	      // we set a new equation u_{num_per} = phase*u_{num_orig}
	      mat_sp.ReallocateRow(num_per, 2); 
	      if (num_per > num_orig)
		{  
		  mat_sp.Index(num_per, 0) = num_orig; 
		  mat_sp.Value(num_per, 0) = -value; 
		  mat_sp.Index(num_per, 1) = num_per; 
		  mat_sp.Value(num_per, 1) = 1.0;
		}
	      else
		{
		  mat_sp.Index(num_per, 1) = num_orig; 
		  mat_sp.Value(num_per, 1) = -value;
		  mat_sp.Index(num_per, 0) = num_per; 
		  mat_sp.Value(num_per, 0) = 1.0;
		}
            }
      }
  }
  
  
  //! modifies matrix in order to take into account quasi-periodic conditions
  void VarBoundaryCondition_Base
  ::SetPeriodicCondition(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>& mat_sp)
  {
    int nb_periodic_dof = this->GetNbPeriodicDof();
    if ((nb_periodic_dof <= 0)||(var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS))
      return;
    
    abort();
    Complex_wp value, value_conj; int nb = 0;
    IVect col_interac; Vector<Complex_wp> val_interac;
    int nb_unknowns = var_problem.nb_unknowns;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nb_unknowns = var_problem.nb_unknowns_hdg;
    
    int Nvol = var_problem.offset_dof_unknown(1);
    // DISP(this->mesh.kwave_periodicity); 
    for (int i = 0; i < nb_periodic_dof; i++)
      {
	int num_per_, num_orig_;
	this->GetPeriodicDofNumbers(i, num_per_, num_orig_);
	
	// lines of periodic dofs are added to original dofs	
	this->GetPeriodicPhase(i, value);
        value_conj = Real_wp(1) / value;
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    int num_per = num_per_ + m*Nvol;
	    int num_orig = num_orig_ + m*Nvol;
	    	    
	    // number of elements on the row
	    nb = mat_sp.GetRealRowSize(num_per);
	    int ni = mat_sp.GetImagRowSize(num_per);
	    // (col_interac,val_interac) -> row sparsely stored 
	    col_interac.Reallocate(nb+ni); val_interac.Reallocate(nb+ni);
	    for (int p = 0; p < nb; p++)
	      {
		col_interac(p) = mat_sp.IndexReal(num_per, p);
		val_interac(p) = value_conj*mat_sp.ValueReal(num_per, p);
	      }
	    
	    for (int p = 0; p < ni; p++)
	      {
		col_interac(nb+p) = mat_sp.IndexImag(num_per, p);
		val_interac(nb+p) = Iwp*value_conj*mat_sp.ValueImag(num_per, p);
	      }
	    
	    // we add to the line L_{num_orig), the line phase*L_{num_per}
	    mat_sp.AddInteractionRow(num_orig, nb, col_interac, val_interac);
	    
	    // we set a new equation u_{num_per} = phase*u_{num_orig}
	    mat_sp.ReallocateRealRow(num_per, 2);
	    mat_sp.ReallocateImagRow(num_per, 1); 
	    if (num_per > num_orig)
              {
                mat_sp.IndexReal(num_per, 0) = num_orig; 
                mat_sp.ValueReal(num_per, 0) = -real(value);
                mat_sp.IndexReal(num_per, 1) = num_per; 
                mat_sp.ValueReal(num_per, 1) = 1.0;
              }
            else
              {
                mat_sp.IndexReal(num_per, 1) = num_orig; 
                mat_sp.ValueReal(num_per, 1) = -real(value);
                mat_sp.IndexReal(num_per, 0) = num_per; 
                mat_sp.ValueReal(num_per, 0) = 1.0;
              }
            
            mat_sp.IndexImag(num_per, 0) = num_orig; 
	    mat_sp.ValueImag(num_per, 0) = -imag(value);
	    
	  }
      }
  }
#endif


  //! modifies right hand side because of quasi-periodic conditions
  template<class T>
  void VarBoundaryCondition_Base
  ::ApplyPeriodicCondition(Vector<T>& brhs) const
  {
    if (var_problem.GetMeshNumberingBase(0).GetFormulationForPeriodicCondition()
	!= MeshNumbering_Base<Real_wp>::STRONG_PERIODIC)
      return;
    
    int nb_periodic_dof = this->GetNbPeriodicDof();
    int nb_proc = var_comm.GetNbProcPerMode();
    //int rank_proc = var_comm.GetRankProcMode();
    int nb_unknowns = var_problem.nb_unknowns;
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      nb_unknowns = var_problem.nb_unknowns_hdg;
    
    int Nvol = var_problem.offset_dof_unknown(1);
    if (nb_proc <= 1)
      {
        for (int i = 0; i < nb_periodic_dof; i++)
          {
	    int num_per, num_orig;
	    this->GetPeriodicDofNumbers(i, num_per, num_orig);
            
            bool continuer = !this->IsDofDirichlet(num_orig);
            
            if ((num_per != num_orig) && continuer)
              {
                // lines of periodic dofs are added to original dofs
                T value, value_conj;
                this->GetPeriodicPhase(i, value);
		value_conj = Real_wp(1) / value;
            
                for (int m = 0; m < nb_unknowns; m++)
                  {
                    int offset = m*Nvol;
                    brhs(num_orig + offset) += brhs(num_per+offset) * value_conj;
                    brhs(num_per+offset) = 0 ;	    
                  }
              }
          }
      }
    else
      {
	// not implemented
      }
    
  }

  
  // PERIODIC CONDITION //
  ////////////////////////
  


  ////////////////
  // PML LAYERS //
  
    
  //! Computation of boundaries, and elements belongin to PML layers are marked
  /*!
    \param[in] alpha dummy parameter (not used)
  */
  template<class Dimension>
  void VarBoundaryCondition_Base::FindElementsInsidePML(DistributedProblem<Dimension>& var)
  {
    Mesh<Dimension>& mesh = var.mesh;
    
    // computation of boundaries for each domain, and taking the min/max of all elements
    // except elements located inside the PML
    Real_wp x0 = 1e300; Real_wp x1 = -1e300;
    Real_wp y0 = 1e300; Real_wp y1 = -1e300;
    Real_wp z0 = 1e300; Real_wp z1 = -1e300;

    Vector<bool> VertexUsed(mesh.GetNbVertices());
    VertexUsed.Fill(false);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (!mesh.Element(i).IsPML())
	for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
	  VertexUsed(mesh.Element(i).numVertex(j)) = true;

    for (int i = 0; i < mesh.GetNbVertices(); i++)
      if (VertexUsed(i))
	{
	  x0 = min(x0, mesh.Vertex(i)(0));
	  x1 = max(x1, mesh.Vertex(i)(0));
	  y0 = min(y0, mesh.Vertex(i)(1));
	  y1 = max(y1, mesh.Vertex(i)(1));
	  if (Dimension::dim_N == 3)
	    {
	      z0 = min(z0, mesh.Vertex(i)(2));
	      z1 = max(z1, mesh.Vertex(i)(2));
	    }
	}
    
    Real_wp xmin(x0), xmax(x1), ymin(y0), ymax(y1), zmin(z0), zmax(z1);
#ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    MpiAllreduce(var.comm_group_mode, &x0, xtmp, &xmin, 1, MPI_MIN);
    MpiAllreduce(var.comm_group_mode, &x1, xtmp, &xmax, 1, MPI_MAX);
    MpiAllreduce(var.comm_group_mode, &y0, xtmp, &ymin, 1, MPI_MIN);
    MpiAllreduce(var.comm_group_mode, &y1, xtmp, &ymax, 1, MPI_MAX);
    MpiAllreduce(var.comm_group_mode, &z0, xtmp, &zmin, 1, MPI_MIN);
    MpiAllreduce(var.comm_group_mode, &z1, xtmp, &zmax, 1, MPI_MAX);
#endif

    var.SetComputationalDomain(xmin, xmax, ymin, ymax, zmin, zmax);
    
    // now all elements in PML layers are marked
    nb_eltPML = 0;
    int ref_max = 0;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        if (mesh.Element(i).IsPML())         
          nb_eltPML++;
        
        ref_max = max(ref_max, mesh.Element(i).GetReference());
      }
    
    nb_eltPML_all = nb_eltPML;
    
#ifdef SELDON_WITH_MPI
    // counting number of PML elements in all the mesh
    MPI_Allreduce(&nb_eltPML, &nb_eltPML_all, 1, MPI_INTEGER, MPI_SUM, var.comm_group_mode);
#endif
    
    // computing maximal speed in PML layers
    max_velocity_pml = 0.0;
    Vector<bool> RefUsed(ref_max+1);
    RefUsed.Fill(false);
    if (nb_eltPML_all == 0)
      {          
        // max_velocity_pml will be the velocity of media close to absorbing boundary condition
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          if (mesh.GetBoundaryCondition(mesh.BoundaryRef(i).GetReference())
              == BoundaryConditionEnum::LINE_ABSORBING)
            {
              int num_elem = mesh.BoundaryRef(i).numElement(0);
              int ref = mesh.Element(num_elem).GetReference();
              RefUsed(ref) = true;
            }
      }
    else
      {
        for (int i = 0; i < mesh.GetNbElt(); i++)
          if (mesh.Element(i).IsPML())
            {
              int ref = mesh.Element(i).GetReference();
              RefUsed(ref) = true;
            }
      }
    
    // for physical indexes (epsilon, mu, sigma, etc)
    // we may want to know the value at infinity
    this->SetPhysicalIndexAtInfinity(RefUsed);

    // taking the maximum among elements of the current processor
    for (int ref = 0; ref <= ref_max; ref++)
      if (RefUsed(ref))
        {
          Real_wp velocity = var.GetVelocityOfMedia(ref);
          max_velocity_pml = max(max_velocity_pml, velocity);
        }
    
    // then max_velocity_pml is the maximum between processors
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    Real_wp max_speed = max_velocity_pml;
    MpiAllreduce(var.comm_group_mode, &max_speed, xtmp, &max_velocity_pml, 1, MPI_MAX);    
#else
    int rank_proc(0);
#endif
    
    if ((rank_proc == 0)||(var.print_level >= 10))
      if (var.print_level >= 1)
        cout << rank_proc << " Number elements inside PML " << nb_eltPML << endl;
    
  } // end method
  

  //! Evaluating PML damping function sigma(dx)
  /*!
    \param[in] dx distance between the interface and the current point of the PML
    \param[in] dsig coefficient of multiplication
    \param[in] a thickness of the PML (somehow maximal value of dx)
    \param[out] zeta damping coefficient sigma(dx)
    \param[out] zeta primitive of damping coefficient \int sigma(dx)
    Usually the damping function is a parabola with sigma(0) = 0
   */
  void VarBoundaryCondition_Base
  ::EvaluateFunctionTauPML(const Real_wp& dx, const Real_wp& dsig, const Real_wp& a,
			   Real_wp& zeta, Real_wp& zeta_primitive) const
  {
    switch (function_damping_pml)
      {
      case PML_PARABOLE :
        zeta = dsig*dx*dx;
        zeta_primitive = zeta*dx/3.0;
        break;
      case PML_CONSTANT :
        zeta = dsig*a*a;
        zeta_primitive = zeta*dx;
        break;
      case PML_LINEAR :
        zeta = dsig*abs(dx)*a;
        zeta_primitive = zeta*dx/2.0;
        break;
      case PML_SHIFTED_PARABOLE :
        zeta = dsig*(dx*dx + offset_damping_pml*a*a);
        zeta_primitive = dsig*dx*(dx*dx/3.0 + offset_damping_pml*a*a);
        break;
      }
  }
  
      
  //! PML damping coefficients
  /*!
    complex change of variable \tilde{x} = x + i/\omega \int_a^x sigma(s) ds
    This function computes the new variable tilde{x} and the
    coefficient zeta(i) = 1 / ( 1 + i sigma(i)/omega) 
    For mixed formulation (i.e. formulation that can be directly used for time-domain simulation)
    we set zeta(i) = sigma(i)
    In the equation, any derivative d/dx_i should be multiplied by zeta(i)
    \param[in] point point where coefficients need to be evaluated
    \param[out] zeta damping coefficient
    \param[out] point_prime modified point \tilde{x}
   */
  template<class Dimension>
  void VarBoundaryCondition_Base
  ::GetDampingFactorPML(Mesh<Dimension>& mesh, int num_pml, int type_pml,
			const typename Dimension::R_N& point,
			typename Dimension::R_N_Complex_wp& zeta,
			typename Dimension::R_N_Complex_wp& point_prime)
  {
    typename Dimension::R_N sigma, sigma_p; Real_wp one(1); Complex_wp coef_om;
    GetDampingTauPML(mesh, num_pml, type_pml, point, sigma, sigma_p);
    
    //if (!this->dispersive_pml)
    //  coef_om = Iwp / this->omega_pml;
    //else
    coef_om = Iwp / var_problem.GetOmega();
    
    if (mesh.GetPmlArea(num_pml).GetRadiusPML() > 0)
      {
        // radial pml
        typename Dimension::R_N x0 = mesh.GetPmlArea(num_pml).GetOriginRadialPML();
        Real_wp r = point.Distance(x0);
        
        zeta(0) = one / (one + coef_om*sigma(0));
        point_prime(0) = r + coef_om*sigma_p(0);
        
        // zeta(1) = (cos teta, sin teta)
        zeta(1) = Complex_wp(sigma(1), sigma_p(1));
        point_prime(1) = point_prime(0)/r;
        return;
      }
    
    if (!var_problem.FirstOrderFormulation())
      {
        if (abs(var_problem.GetOmega()) != Real_wp(0))
          for (int i = 0; i < zeta.GetM(); i++)
            zeta(i) = one/(one + coef_om*sigma(i));
      }
    else
      for (int i = 0; i < zeta.GetM(); i++)
        zeta(i) = sigma(i);
    
    for (int i = 0; i < zeta.GetM(); i++)
      point_prime(i) = point(i) + coef_om*sigma_p(i);
  }
  
  
  //! PML damping coefficients
  /*!
    \param[in] point point where coefficients need to be evaluated
    \param[out] zeta damping coefficient
    \param[out] zeta_primitive primitive of damping coefficient 
   */
  void VarBoundaryCondition_Base
  ::GetDampingTauPML(const Mesh<Dimension2>& mesh, int num_pml, int type_pml,
		     const R2& point, R2& zeta, R2& zeta_primitive) const
  {    
    const PmlRegionParameter<Dimension2>& pml_area = mesh.GetPmlArea(num_pml);
    Real_wp a = pml_area.GetThicknessPML();
    Real_wp dsig = log(1000e0)*3e0*this->vsigma*this->max_velocity_pml/(2*pow(a, 3));
    
    zeta.Fill(0); zeta_primitive.Fill(0);
    // for radial pml
    Real_wp Rmax = pml_area.GetRadiusPML();
    if (Rmax > 0)
      {
        R2 x0 = pml_area.GetOriginRadialPML();
        Real_wp r = point.Distance(x0);
        EvaluateFunctionTauPML(r-Rmax, dsig, a, zeta(0), zeta_primitive(0));
        
        // we store cos(theta) and sin(teta) in second coordinate of zeta/zeta_primitive
        zeta(1) = (point(0) - x0(0))/r;
        zeta_primitive(1) = (point(1) - x0(1))/r;
        
        return;
      }

    // damping in x-direction
    if (type_pml & PmlRegionParameter<Dimension2>::PML_Xm)
      {
        EvaluateFunctionTauPML(point(0) - var_problem.GetXmin(), dsig, a,
                               zeta(0), zeta_primitive(0));
      } 
    else if (type_pml & PmlRegionParameter<Dimension2>::PML_Xp)
      {
        EvaluateFunctionTauPML(point(0) - var_problem.GetXmax(), dsig, a,
                               zeta(0), zeta_primitive(0));
      }
    
    // damping in y-direction
    if (type_pml & PmlRegionParameter<Dimension2>::PML_Ym)
      {
        EvaluateFunctionTauPML(point(1) - var_problem.GetYmin(), dsig, a,
                               zeta(1), zeta_primitive(1));
      }
    else if (type_pml & PmlRegionParameter<Dimension2>::PML_Yp)
      {
        EvaluateFunctionTauPML(point(1) - var_problem.GetYmax(), dsig, a,
                               zeta(1), zeta_primitive(1));
      }

    // we can add some damping in x or y to ensure stability
    // (e.g. in elastodynamics)
    if (coef_mixed_damping_pml != 0)
      {
        Real_wp dx = zeta(0), dy = zeta(1);
        zeta(0) = max(dx, coef_mixed_damping_pml*dy);
        zeta(1) = max(coef_mixed_damping_pml*dx, dy);
      }
  }
  
  
  //! PML damping coefficients
  /*!
    \param[in] point point where coefficients need to be evaluated
    \param[out] zeta damping coefficient
    \param[out] zeta_prime primitive of damping coefficient 
   */
  void VarBoundaryCondition_Base
  ::GetDampingTauPML(const Mesh<Dimension3>& mesh, int num_pml, int type_pml,
		     const R3& point, R3& zeta, R3& zeta_primitive) const
  {
#ifdef MONTJOIE_WITH_THREE_DIM
    Real_wp a = mesh.GetPmlArea(num_pml).GetThicknessPML();
    Real_wp dsig = log(1000e0)*3e0*this->vsigma*this->max_velocity_pml/(2*pow(a, 3));
    
    zeta.Zero(); zeta_primitive.Zero();
    // damping in x-direction
    if (type_pml & PmlRegionParameter<Dimension3>::PML_Xm)
      {
        EvaluateFunctionTauPML(point(0) - var_problem.GetXmin(), dsig, a,
                               zeta(0), zeta_primitive(0));
      }
    else if (type_pml & PmlRegionParameter<Dimension3>::PML_Xp)
      {
        EvaluateFunctionTauPML(point(0) - var_problem.GetXmax(), dsig, a,
                               zeta(0), zeta_primitive(0));
      }
    
    // damping in y-direction
    if (type_pml & PmlRegionParameter<Dimension3>::PML_Ym)
      {
        EvaluateFunctionTauPML(point(1) - var_problem.GetYmin(), dsig, a,
                               zeta(1), zeta_primitive(1));
      }
    else if (type_pml & PmlRegionParameter<Dimension3>::PML_Yp)
      {
        EvaluateFunctionTauPML(point(1) - var_problem.GetYmax(), dsig, a,
                               zeta(1), zeta_primitive(1));
      }
    
    // damping in z-direction
    if (type_pml & PmlRegionParameter<Dimension3>::PML_Zm)
      {
        EvaluateFunctionTauPML(point(2) - var_problem.GetZmin(), dsig, a,
                               zeta(2), zeta_primitive(2));
      }
    else if (type_pml & PmlRegionParameter<Dimension3>::PML_Zp)
      {
        EvaluateFunctionTauPML(point(2) - var_problem.GetZmax(), dsig, a,
                               zeta(2), zeta_primitive(2));
      }
    
    // we can add some damping in x/y/z to ensure stability
    // (e.g. in elastodynamics)
    if (coef_mixed_damping_pml != 0)
      {
        Real_wp dx = zeta(0), dy = zeta(1), dz = zeta(2);
        zeta(0) = max(dx, coef_mixed_damping_pml*max(dy, dz));
        zeta(1) = max(coef_mixed_damping_pml*max(dx, dz), dy);
        zeta(2) = max(coef_mixed_damping_pml*max(dx, dy), dz);
      }
#endif
  }

  
  // PML LAYERS //  
  ////////////////
  

  /*****************
   * Other Methods *
   *****************/
  
  
  //! computation of surfacic term of the type
  //! adding \f$ alpha \int_\Gamma \, f(\varphi_i, \nabla \varphi_i) \varphi_j
  //! + fgrad(\varphi_i, \nabla \varphi_i) \nabla \varphi_j dx \f$
  /*!
    \param[in] ref_cond boundary condition associated to \f$ \Gamma \f$
    \param[in,out] mat_sp sparse matrix where interactions will be added
    \param[in] fimpedance description of the impedance functions f and fgrad
    f(varphi, nabla varphi) is evaluated when calling f.ApplyImpedancePhi
    fgrad(varphi, nabla varphi) is evaluated when calling f.ApplyImpedanceGrad
  */
  template<class T, class Dimension>
  void VarBoundaryCondition_Base
  ::AddMatrixImpedanceBoundary(const T& alpha, const Vector<int>& ref_cond, int ref_target,
			       const GlobalGenericMatrix<T>& nat_mat,
			       VirtualMatrix<T>& mat_sp, int offset_row, int offset_col,
			       ImpedanceFunction_Base<T, Dimension>& fimpedance,
			       bool changeColonne, bool changeLigne,
			       const VarProblem<Dimension>& var)
  {
    if (var.FormulationDG() != ElementReference_Base::CONTINUOUS)
      {
        // for LDG or HDG formulation, boundary conditions are set elsewhere
        if (var.FirstOrderFormulationDG())
          return;
      }
    
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;

    int nb_mesh_num = var.GetNbMeshNumberings();
    Vector<Matrix<TinyVector<Real_wp, 1> > > ValuePhi_H1(nb_mesh_num);
    Vector<Matrix<TinyVector<Real_wp, Dimension::dim_N> > > GradientPhi_H1(nb_mesh_num);

    Vector<Matrix<TinyVector<Real_wp, Dimension::dim_N> > > ValuePhi_Hdiv(nb_mesh_num);
    Vector<Matrix<TinyVector<Real_wp, 1> > > DivPhi_Hdiv(nb_mesh_num);

    enum {nb_u_hcurl = FiniteElementHcurl<Dimension>::nb_components_u,
          nb_du_hcurl = FiniteElementHcurl<Dimension>::nb_components_grad};
    
    Vector<Matrix<TinyVector<Real_wp, nb_u_hcurl> > > ValuePhi_Hcurl(nb_mesh_num);
    Vector<Matrix<TinyVector<Real_wp, nb_du_hcurl> > > CurlPhi_Hcurl(nb_mesh_num);
    
    // local vectors
    int size_phi = 0;
    int size_dphi = 0;
    Vector<int> off_phi(var.nb_unknowns+1), off_dphi(var.nb_unknowns+1);
    off_phi(0) = 0; off_dphi(0) = 0;
    Vector<int> type_u(nb_mesh_num); Vector<int> IndexMesh(var.nb_unknowns);
    Vector<Vector<int> > IndexUnknown(nb_mesh_num);
    for (int i = 0; i < var.nb_unknowns; i++)
      {
        int n = var.mesh_num_unknown(i);
        if (n == 0)
          type_u(n) = var.type_element;
        else
          type_u(n) = var.other_type_element(n-1);

        // arrays to go back and forth between mesh numbering and unknown numbering
        IndexMesh(i) = IndexUnknown(n).GetM();
        IndexUnknown(n).PushBack(i);
        
        int nb_u = var.GetNbComponentsUnknown(n);
        int nb_du = var.GetNbComponentsGradient(n);
        size_phi += nb_u; off_phi(i+1) = off_phi(i) + nb_u;
        size_dphi += nb_du; off_dphi(i+1) = off_dphi(i) + nb_du;
      }
    
    Vector<T> f_phi(size_phi), f_dphi(size_dphi);
    Vector<Vector<Vector<T> > > feval_phi(nb_mesh_num), feval_dphi(nb_mesh_num);
    Vector<Vector<Vector<T> > > contrib(nb_mesh_num), contrib_grad(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      {
        feval_phi(i).Reallocate(IndexUnknown(i).GetM());
        feval_dphi(i).Reallocate(IndexUnknown(i).GetM());
        contrib(i).Reallocate(IndexUnknown(i).GetM());
        contrib_grad(i).Reallocate(IndexUnknown(i).GetM());
      }

    // vertices of the elements
    typename Dimension::VectR_N s;
    
    // values to add to the matrix
    IVect col_interac;
    Vector<T> val_interac;
    T val;
    int c_interac = 0;
    int num_dof, num_col;

    const Mesh<Dimension>& mesh = var.mesh;
    Vector<IVect> num_dof_glob(nb_mesh_num);
    
    int rank_proc = var_comm.GetRankProcMode();

    // loop on referenced edges/faces
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
       	// reference of the edge/face
	int ref_line = mesh.BoundaryRef(i).GetReference();
	
        // the boundary condition is what the user asked
    	if ((ref_line > 0) && (ref_line < ref_cond.GetM()))
	  if (ref_target == ref_cond(ref_line))
	    {
	      int iquad = mesh.BoundaryRef(i).numElement(0);
	      int num_face = i;
	      int num_loc = mesh.Element(iquad).GetPositionBoundary(num_face);

	      const ElementReference_Dim<Dimension>& Fb = var.GetReferenceElement(iquad);
    	      int nb_points_quadrature_edge = Fb.GetNbQuadBoundary(num_loc);
	      int ref_domain = mesh.Element(iquad).GetReference();
              int nb_dof_all = 0;

	      // allocating these arrays to the maximal value
	      bool presence_grad = fimpedance.PresenceGradient();
	      for (int k = 0; k < var.nb_unknowns; k++)
		{
                  int n = var.mesh_num_unknown(k); int p = IndexMesh(k);
                  int nb_u = off_phi(n+1)-off_phi(n);
                  int nb_du = off_dphi(n+1) - off_dphi(n);
		  feval_phi(n)(p).Reallocate(nb_points_quadrature_edge*nb_u);
		  if (presence_grad)
		    feval_dphi(n)(p).Reallocate(nb_points_quadrature_edge*nb_du);

                  int nb_dof_loc = var.GetReferenceElement(iquad, n).GetNbDof();
		  contrib(n)(p).Reallocate(nb_dof_loc);
		  if (presence_grad)
		    contrib_grad(n)(p).Reallocate(nb_dof_loc);

                  nb_dof_all += nb_dof_loc;
		}
              
	      col_interac.Reallocate(nb_dof_all); col_interac.Fill(-1);
	      val_interac.Reallocate(nb_dof_all); val_interac.Zero();
	      
	      // we get the vertices of the element iquad
	      mesh.GetVerticesElement(iquad, s);
	      
	      // transformation F_i and DF_i are applied for quadrature points
	      // PointsElem => coordinates of quadrature points
	      // MatriceElem => jacobian matrices DFi on quadrature points
	      Fb.FjElemQuadrature(s, PointsElem, mesh, iquad);
	      Fb.DFjElemQuadrature(s, PointsElem, MatricesElem, mesh, iquad);
	      
	      // restriction of DF_i on surface
	      Fb.FjSurfaceElem(s, PointsElem, mesh, iquad, num_loc);    
	      Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, mesh, iquad, num_loc);
	      
	      // loop on quadrature points of the edge
	      int num_edge = num_face;
	      for (int k = 0; k < nb_points_quadrature_edge; k++)
		{
		  // the impedances here are assumed to be computed and stored in the 
		  // class derived from ImpedanceFunction_Base	
		  // evaluation of the impedance operator
		  fimpedance.EvaluateImpedancePhi(i, iquad, num_edge, num_loc, k, nat_mat, ref_domain,
						  PointsElem, MatricesElem);
		  
		  // and for gradient        
		  if (presence_grad)
		    fimpedance.EvaluateImpedanceGrad(i, iquad, num_edge, num_loc, k, nat_mat, ref_domain,
						     PointsElem, MatricesElem);
		}   

	      for (int n = 0; n < nb_mesh_num; n++)
                {
                  const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(n);
                  
                  switch(type_u(n))
                    {
                    case 1:
                      // evaluation of basis functions and gradient of basis functions on quadrature points
                      var.GetReferenceElementH1(iquad, n)
                        .ComputeValuesPhiBoundary(ValuePhi_H1(n),
                                                  MatricesElem, mesh_num, iquad, num_loc);
                  
                      // meme chose pour les gradients
                      var.GetReferenceElementH1(iquad, n)
                        .ComputeValuesGradientPhiBoundary(GradientPhi_H1(n), MatricesElem,
                                                          mesh_num, iquad, num_loc);

                      break;
                    case 2:
                      // evaluation of basis functions and gradient of basis functions on quadrature points
                      var.GetReferenceElementHcurl(iquad, n)
                        .ComputeValuesPhiBoundary(ValuePhi_Hcurl(n),
                                                  MatricesElem, mesh_num, iquad, num_loc);
                  
                      // meme chose pour les gradients
                      var.GetReferenceElementHcurl(iquad, n)
                        .ComputeValuesGradientPhiBoundary(CurlPhi_Hcurl(n), MatricesElem,
                                                          mesh_num, iquad, num_loc);
                      break;
                    case 3:
                      // evaluation of basis functions and gradient of basis functions on quadrature points
                      var.GetReferenceElementHdiv(iquad, n)
                        .ComputeValuesPhiBoundary(ValuePhi_Hdiv(n),
                                                  MatricesElem, mesh_num, iquad, num_loc);
                  
                      // meme chose pour les gradients
                      var.GetReferenceElementHdiv(iquad, n)
                        .ComputeValuesGradientPhiBoundary(DivPhi_Hdiv(n), MatricesElem,
                                                          mesh_num, iquad, num_loc);
                    }
                  
                  // we retrieve degrees of freedom
                  num_dof_glob(n) = var.GetDofNumberOnElement(iquad, n);
                }
              
	      // loop on unknowns (rows)
	      for (int m = 0; m < var.nb_unknowns; m++)
		{
                  int nm = var.mesh_num_unknown(m);
                  const ElementReference_Dim<Dimension>& Fb_m = var.GetReferenceElement(iquad, nm);
                  int type_m = type_u(nm);
		  FillZero(f_phi); FillZero(f_dphi);
                  
		  // loop on degrees of freedom involved in the computation (rows)
		  // a "tangential dof" if the impedance function involves only these dofs
		  // any degree of freedom if the impedance function involves other dofs
		  for (int k = 0; k < Fb_m.GetNbDof(); k++)
		    if ((!fimpedance.InvolveOnlyTangentialDofs() || Fb_m.IsTangentialDof(k, num_loc) ) 
			&& (num_dof_glob(nm)(k) >= 0) )
		      {
			num_dof = num_dof_glob(nm)(k) + var.offset_dof_unknown(m);
			
			// loop on quadrature points
			// we apply impedance operator on basis functions phi
			for (int j = 0; j < nb_points_quadrature_edge; j++)
			  {
			    switch (type_m)
                              {
                              case 1: fimpedance.ApplyImpedancePhi_H1(m, j, off_phi(m), ValuePhi_H1(nm)(k, j),
                                                                      GradientPhi_H1(nm)(k, j), f_phi);
                                break;
                              case 2: fimpedance.ApplyImpedancePhi_Hcurl(m, j, off_phi(m), ValuePhi_Hcurl(nm)(k, j),
                                                                         CurlPhi_Hcurl(nm)(k, j), f_phi);
                                break;
                              case 3: fimpedance.ApplyImpedancePhi_Hdiv(m, j, off_phi(m), ValuePhi_Hdiv(nm)(k, j),
                                                                        DivPhi_Hdiv(nm)(k, j), f_phi);
                                break;
                              }
                            
			    for (int p = 0; p < var.nb_unknowns; p++)
			      {
                                int np = var.mesh_num_unknown(p);
                                int ip = IndexMesh(p);
				// feval_phi(n)(j*nb_u+k) = f_phi(n*nb_u + k)
                                switch (type_u(np))
                                  {
                                  case 1: feval_phi(np)(ip)(j) = f_phi(off_phi(p)); break;
                                  case 2:
                                    TinyVectorLoop<nb_u_hcurl>::ExtractVector(f_phi, off_phi(p), feval_phi(np)(ip), j*nb_u_hcurl);
                                    break;
                                  case 3:
                                    TinyVectorLoop<Dimension::dim_N>::ExtractVector(f_phi, off_phi(p), feval_phi(np)(ip), j*Dimension::dim_N);
                                  }
			      }
			  }
                        
			// computation of the integral against basis functions phi
			for (int n = 0; n < nb_mesh_num; n++)
                          {
                            const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(n);
                            
                            switch (type_u(n))
                              {
                              case 1:
                                var.GetReferenceElementH1(iquad, n)
                                  .ComputeIntegralSurface(MatricesElem, feval_phi(n),
                                                          contrib(n), mesh_num, iquad, num_loc);
                                break;
                              case 2:
                                var.GetReferenceElementHcurl(iquad, n)
                                  .ComputeIntegralSurface(MatricesElem, feval_phi(n),
                                                          contrib(n), mesh_num, iquad, num_loc);
                                break;
                              case 3:
                                var.GetReferenceElementHdiv(iquad, n)
                                  .ComputeIntegralSurface(MatricesElem, feval_phi(n),
                                                          contrib(n), mesh_num, iquad, num_loc);
                                break;
                              }
                          }
                        
			if (presence_grad)
			  {
			    // same operations for gradient
			    for (int j = 0; j < nb_points_quadrature_edge; j++)
			      {
                                switch (type_m)
                                  {
                                  case 1:
                                    fimpedance.ApplyImpedanceGrad(m, j, off_dphi(m), ValuePhi_H1(nm)(k,j),
                                                                  GradientPhi_H1(nm)(k, j), f_dphi);
                                    break;
                                  case 2:
                                    fimpedance.ApplyImpedanceCurl(m, j, off_dphi(m), ValuePhi_Hcurl(nm)(k,j),
                                                                  CurlPhi_Hcurl(nm)(k, j), f_dphi);
                                    break;
                                  case 3:
                                    fimpedance.ApplyImpedanceDiv(m, j, off_dphi(m), ValuePhi_Hdiv(nm)(k,j),
                                                                 DivPhi_Hdiv(nm)(k, j), f_dphi);
                                    break;
                                  }
                                
                                for (int p = 0; p < var.nb_unknowns; p++)
                                  {
                                    int np = var.mesh_num_unknown(p);
                                    int ip = IndexMesh(p);
                                    // feval_phi(n)(j*nb_u+k) = f_phi(n*nb_u + k)
                                    switch (type_u(np))
                                      {
                                      case 1:
                                        TinyVectorLoop<Dimension::dim_N>::ExtractVector(f_dphi, off_dphi(p), feval_dphi(np)(ip), j*Dimension::dim_N);
                                        break;
                                      case 2:
                                        TinyVectorLoop<nb_du_hcurl>::ExtractVector(f_dphi, off_dphi(p), feval_dphi(np)(ip), j*nb_du_hcurl);
                                        break;
                                      case 3:
                                        feval_phi(np)(ip)(j) = f_dphi(off_dphi(p)); break;
                                      }
                                  }
                              }
			    
			    // computation of the integral against basis functions grad phi
                            for (int n = 0; n < nb_mesh_num; n++)
                              {
                                const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(n);
                                
                                switch (type_u(n))
                                  {
                                  case 1:
                                    var.GetReferenceElementH1(iquad, n)
                                      .ComputeIntegralSurfaceGradient(MatricesElem, feval_dphi(n),
                                                          contrib_grad(n), mesh_num, iquad, num_loc);
                                    break;
                                  case 2:
                                    var.GetReferenceElementHcurl(iquad, n)
                                      .ComputeIntegralSurfaceGradient(MatricesElem, feval_dphi(n),
                                                          contrib_grad(n), mesh_num, iquad, num_loc);
                                    break;
                                  case 3:
                                    var.GetReferenceElementHdiv(iquad, n)
                                      .ComputeIntegralSurfaceGradient(MatricesElem, feval_dphi(n),
                                                          contrib_grad(n), mesh_num, iquad, num_loc);
                                    break;
                                  }
                              }
			  }
			
			c_interac = 0;
			if (changeColonne)
			  {
			    // now loops on columns
			    for (int n = 0; n < var.nb_unknowns; n++)
                              {
                                int nn = var.mesh_num_unknown(n); int in = IndexMesh(n);
                                const ElementReference_Dim<Dimension>& Fb_n = var.GetReferenceElement(iquad, nn);
                                for (int j = 0; j < Fb_n.GetNbDof(); j++)
                                  if (!fimpedance.InvolveOnlyTangentialDofs() || Fb.IsTangentialDof(j, num_loc))
                                    {
                                      if (presence_grad)
                                        val = alpha*(contrib(nn)(in)(j) + contrib_grad(nn)(in)(j));
                                      else
                                        val = alpha*contrib(nn)(in)(j);
                                      
                                      // we add value to the matrix if not neglectible
                                      if (( abs(val) > var.GetThresholdMatrix()) && (num_dof_glob(nn)(j) >= 0))
                                        {
                                          // numbers of the columns
                                          num_col = num_dof_glob(nn)(j) + var.offset_dof_unknown(n);
                                          int proc = ProcColumnNumbers_Impedance(num_col);
                                          if (proc == rank_proc)
                                            {
                                              col_interac(c_interac) = NewColumnNumbers_Impedance(num_col);
                                              
                                              if (col_interac(c_interac) < 0)
                                                {
                                                  cout << "You need to specify in "
                                                       << " NewColumnNumbers_Impedance positive "
                                                       << " dof numbers " << endl;
                                                  DISP(num_col); DISP(col_interac(c_interac));
                                                  DISP(contrib(nn)(in)(j));
                                                  abort();
                                                }
                                              else
                                                {                      
                                                  col_interac(c_interac) += offset_col;
                                                  
                                                  // values
                                                  val_interac(c_interac) = val;
                                                  c_interac++;
                                                }
                                            }
                                          else
                                            {
#ifdef SELDON_WITH_MPI
                                              int p_row = rank_proc;
                                              if (changeLigne)
                                                {
                                                  num_dof = num_dof_glob(nm)(k) + var.offset_dof_unknown(m);
                                                  if (num_dof_glob(nm)(k) >= 0)
                                                    {
                                                      p_row = ProcRowNumbers_Impedance(num_dof);
                                                      num_dof = NewRowNumbers_Impedance(num_dof);
                                                    }
                                                }
                                              
                                              if ((num_dof >= 0) && (num_dof_glob(nn)(j) >= 0))
                                                {
                                                  if (p_row == rank_proc)
                                                    {
                                                      mat_sp
                                                        .AddDistantInteraction
                                                        (offset_row + num_dof,
                                                         NewColumnNumbers_Impedance(num_col), proc, val);
                                                    }
                                                  else
                                                    {
                                                      cout << "Forbidden case" << endl;
                                                      abort();
                                                    }
                                                }
#endif
                                            }
                                        }
                                    }
                              }
			  }
			else
			  {			    
			    for (int n = 0; n < var.nb_unknowns; n++)
                              {
                                int nn = var.mesh_num_unknown(n); int in = IndexMesh(n);
                                const ElementReference_Dim<Dimension>& Fb_n = var.GetReferenceElement(iquad, nn);
                                for (int j = 0; j < Fb_n.GetNbDof(); j++)
                                  if (!fimpedance.InvolveOnlyTangentialDofs() || Fb_n.IsTangentialDof(j, num_loc))
                                    {
                                      if (presence_grad)
                                        val = alpha*(contrib(nn)(in)(j) + contrib_grad(nn)(in)(j));
                                      else
                                        val = alpha*contrib(nn)(in)(j);
                                      
                                      // we add value to the matrix if not neglectible
                                      if (( abs(val) > var.GetThresholdMatrix()) && (num_dof_glob(nn)(j) >= 0))
                                        {
                                          // numbers of the columns
                                          col_interac(c_interac) = offset_col + num_dof_glob(nn)(j) + var.offset_dof_unknown(n);
                                          
                                          // values
                                          val_interac(c_interac) = val;
                                          c_interac++;
                                        }
                                    }
                              }
			  }
                        
			if (c_interac > 0)
			  {                      
			    int p_row = rank_proc;      
			    if (changeLigne) 
			      {
				num_dof = num_dof_glob(nm)(k) + var.offset_dof_unknown(m);
				p_row = ProcRowNumbers_Impedance(num_dof);
				num_dof = NewRowNumbers_Impedance(num_dof);
				if (num_dof < 0)
				  {
				    cout << "You need to specify in NewRowNumbers_Impedance "
					 << "positive dof numbers " << endl;
				    abort();
				  }
			      }
			    
			    // we add interations in the matrix mat_sp
			    if (p_row == rank_proc)
			      mat_sp.AddInteractionRow(offset_row + num_dof, c_interac,
						       col_interac, val_interac);
			    else
			      {
#ifdef SELDON_WITH_MPI
				for (int ni = 0; ni < c_interac; ni++)
				  mat_sp.AddRowDistantInteraction(num_dof, col_interac(ni),
								  p_row, val_interac(ni));
#endif
			      }
			  }
		      }
		}    	     
	    }
      }
  }
  
  
  //! copies input parameters of another similar problem
  void VarBoundaryCondition_Base::CopyInputData(const VarBoundaryCondition_Base& var)
  {
    supported_components_BC = var.supported_components_BC;
    order_high_conductivity = var.order_high_conductivity;
    
    list_number_mode = var.list_number_mode;
    number_mode = var.number_mode;    
    nb_modes_periodic_xyz = var.nb_modes_periodic_xyz;
    type_symmetry = var.type_symmetry;
    do_not_store_modes = var.do_not_store_modes;

    take_into_account_curvature_for_abc = var.take_into_account_curvature_for_abc;
    grazing_abc = var.grazing_abc;    
    order_ABC = var.order_ABC;
    vsigma = var.vsigma;
    function_damping_pml = var.function_damping_pml;
    offset_damping_pml = var.offset_damping_pml;
    coef_mixed_damping_pml = var.coef_mixed_damping_pml;
  }    


  //! returns the impedance for surfaces defined by ref
  Complex_wp VarBoundaryCondition_Base::ComputeImpedanceCoefficient(const IVect& ref, VectComplex_wp& U)
  {
    Complex_wp vloc; SetComplexZero(vloc);
    return vloc;
  }

  
  /****************************
   * VarBoundaryCondition_Dim *
   ****************************/


  //! default constructor
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>::InitDefaultValues()
  {
    SetComplexOne(Coef_Imped_ABC);
  }


  //! returns size of memory used by this object in bytes
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>::GetMemoryUsed(map<string, size_t>& var) const
  {
    size_t taille = sizeof(*this);
    taille += sizeof(int)*(NewColumnNumbers_Impedance.GetM() + ProcColumnNumbers_Impedance.GetM()
                           +NewRowNumbers_Impedance.GetM() + ProcRowNumbers_Impedance.GetM()
                           + Dirichlet_dof.GetM() + supported_components_BC.GetM() + list_number_mode.GetM());
    
    taille += sizeof(bool)*is_dof_dirichlet.GetM();
    taille += point_dirichlet_condition.GetMemorySize() +
      Seldon::GetMemorySize(Param_condition) + Seldon::GetMemorySize(Glob_TauPML)
      + Seldon::GetMemorySize(Glob_TauPrimePML);
    
    var["BoundaryCondition"] = taille;
  }


  //! retrieves parameters of a edge/face of the mesh
  template<class Complexe, class Dimension>
  const Matrix<int>& VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetGeometryPhaseData(int iquad, int num_loc,
			 int& num_face, int& ref_boundary, int& rf, bool& new_face, int& rot,
			 Complexe& phase, Complexe& phase_conj) const
  {
    const Mesh<Dimension>& mesh = mesh_reference;
    const MeshNumbering<Dimension>& mesh_num = var_comm.GetMeshNumbering(0);
    
    SetComplexOne(phase); SetComplexOne(phase_conj);
    num_face = mesh.Element(iquad).numBoundary(num_loc);
    ref_boundary = mesh.Boundary(num_face).GetReference();
    rf = mesh_num.GetOrderQuadrature(num_face);
    new_face = var_comm.IsNewFace(iquad)(num_loc);
 
    int num_elem2 = mesh.Boundary(num_face).numElement(0);
    if ((num_elem2 == iquad)&&(mesh.Boundary(num_face).GetNbElements()==2))
      num_elem2 = mesh.Boundary(num_face).numElement(1);

    // we retrieve the phase and conjugated phase for quasi-periodic conditions
#ifdef SELDON_WITH_MPI
    bool neighbor_face = (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_NEIGHBOR);
#endif
    int rot1 = mesh.Element(iquad).GetOrientationBoundary(num_loc), rot2 = 0;
    int nv = mesh.Boundary(num_face).GetNbVertices();
    rot = 0;
    if (num_elem2 != iquad)
      {
	int num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(num_face);
	if (num_pos2_face < 0)
	  {
	    int nf2 = mesh_num.GetPeriodicBoundary(num_face);
	    if (num_face > nf2)
	      if (mesh_num.GetFormulationForPeriodicCondition()
		  != MeshNumbering_Base<Real_wp>::STRONG_PERIODIC)
		{
		  this->GetPeriodicPhase(num_face, phase);
		  this->GetPeriodicPhase(nf2, phase_conj);
		}
	    
	    num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(nf2);
	  }
	
	rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
	rot = mesh_num.GetRotationFace(rot1, rot2, nv);
      }
    
#ifdef SELDON_WITH_MPI
    if (neighbor_face)
      {
	int type_per = mesh_num.GetPeriodicityTypeForBoundary(num_face);
	if (type_per >= 0)
	  {
	    int num_face_glob = mesh.GetGlobalBoundaryNumber(num_face);
	    int nf2_glob = mesh_num.GetPeriodicBoundary(num_face);
	    if (num_face_glob > nf2_glob)
	      {
		this->GetPeriodicPhase(num_face, phase);
		phase_conj = Real_wp(1) / phase;
	      }
	  }		  
      }
#endif
    
    return mesh_num.number_map.
      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
  }

  
  //! modification of boundary conditions with a line of the data file
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarBoundaryCondition_Base::SetInputData(description_field, parameters);
    
    if (!description_field.compare("PointDirichlet"))
      {
        R_N ptA;
	if (parameters.GetM() < ptA.GetM())
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "PointDirichlet needs more parameters, for instance :" << endl;
	    cout << "PointDirichlet = xA yA zA" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	for (int k = 0; k < ptA.GetM(); k++)
	  to_num(parameters(k), ptA(k));

        point_dirichlet_condition.PushBack(ptA);
      }
    else if (!description_field.compare("ModifiedImpedance"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "ModifiedImpedance needs more parameters, for instance :" << endl;
	    cout << "ModifiedImpedance = coef" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (parameters(0) == "CURVE")
          take_into_account_curvature_for_abc = true;
        else
          {
            take_into_account_curvature_for_abc = false;
            Coef_Imped_ABC = to_num<Complexe>(parameters(0));
          }
      }
    else if (!description_field.compare("NbModesPeriodic"))
      {
	if (parameters.GetM() < Dimension::dim_N)
	  {
	    cout << "In SetInputData of VarBoundaryCondition" << endl;
	    cout << "NbModesPeriodic needs more parameters, for instance :" << endl;
	    cout << "NbModesPeriodic = nX nY nZ" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        for (int i = 0; i < Dimension::dim_N; i++)
          nb_modes_periodic_xyz(i) = to_num<int>(parameters(i));
      }
  }


  //! checks if ref is not too large, resizing otherwise
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension> ::CheckAndReallocateReferences(int ref)
  {
    if (ref > mesh_reference.GetNbReferences())
      mesh_reference.ResizeNbReferences(ref);

    if (ref > Param_condition.GetM()-1)
      Param_condition.Resize(ref+1);
  }
  
  
  //! collects parameters in line ConditionReference of the data file
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::FillParameterCondition(const Vector<string>& parameters, int position_word, Vector<int>& ref)
  {
    // we store these values in parameters_condition
    Vector<Complexe> parameters_condition(parameters.GetM()-position_word-1);
    for (int i = (position_word+1); i < parameters.GetM(); i++)
      parameters_condition(i-position_word-1) = to_num<Complexe>(parameters(i));
    
    for (int i = 0; i < ref.GetM(); i++)
      Param_condition(ref(i)) = parameters_condition;		
  }
  
  
  //! sets dirichlet for asked points
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::AddDirichletPointConditions(Vector<Vector<int> >& IndexOtherDir)
  {
    Vector<bool> ddl_vertex_n(var_comm.GetNbMeshNumberings());
    ddl_vertex_n.Fill(false); bool ddl_vertex = false;
    for (int n = 0; n < ddl_vertex_n.GetM(); n++)
      if (var_comm.GetMeshNumbering(n).number_map.GetNbDofVertex(0) > 0)
        {
          ddl_vertex = true;
          ddl_vertex_n(n) = true;
        }
    
    for (int m = 0; m < this->point_dirichlet_condition.GetM(); m++)
      {
        // you may require Dirichlet condition on given points
        if (ddl_vertex)
	  {
	    for (int i = 0; i < mesh_reference.GetNbVertices(); i++)
	      if (mesh_reference.Vertex(i) == this->point_dirichlet_condition(m))
		{
                  for (int n = 0; n < ddl_vertex_n.GetM(); n++)
                    if (ddl_vertex_n(n))
                      IndexOtherDir(n)(i) = BoundaryConditionEnum::LINE_DIRICHLET;
                }
          }
      }
  }
  

  //! allocates array for PML coefficients
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>::AllocateTauPML()
  {
    if (this->nb_eltPML > 0)
      {
	this->Glob_TauPML.Reallocate(this->nb_eltPML);
	this->Glob_TauPrimePML.Reallocate(this->nb_eltPML);
      }
  }
  

  //! Computes tau and primitive of tau inside PML
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::ComputeStoreCoefficientsPML(int i1, int num_elem, const VectR_N& AllPoints)
  {
    int nb_points_all = AllPoints.GetM();
    this->Glob_TauPML(i1).Reallocate(nb_points_all);
    this->Glob_TauPrimePML(i1).Reallocate(nb_points_all);
    int num_pml = var_comm.mesh.Element(num_elem).GetNumberPML();
    int type_pml = var_comm.mesh.Element(num_elem).GetTypePML();
    for (int j = 0; j < nb_points_all; j++)
      this->GetDampingFactorPML(AllPoints(j), num_pml, type_pml,
                                this->Glob_TauPML(i1)(j), this->Glob_TauPrimePML(i1)(j));	
  }

  
  //! initialisation of computation for symmetric domains (for example cyclic domains)
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>::InitCyclicDomain()
  {
    Mesh<Dimension>& mesh = mesh_reference;
    
    bool per_theta = false, per_x = false, per_y = false, per_z = false;
    for (int n = 0; n < mesh.GetNbPeriodicReferences(); n++)
      {
        int type = mesh.GetPeriodicityTypeReference(n);
        if (type == BoundaryConditionEnum::PERIODIC_THETA)
          per_theta = true;
        
        if (type == BoundaryConditionEnum::PERIODIC_X)
          per_x = true;

        if (type == BoundaryConditionEnum::PERIODIC_Y)
          per_y = true;

        if (type == BoundaryConditionEnum::PERIODIC_Z)
          per_z = true;
      }
    
    type_symmetry = VarBoundaryCondition_Base::NO_SYMMETRY;
    int Nx = nb_modes_periodic_xyz(0);
    int Ny = nb_modes_periodic_xyz(1);
    int Nz = 1;
    if (Dimension::dim_N == 3)
      Nz = nb_modes_periodic_xyz(2);
    
    if (per_theta)
      {
	Real_wp alpha = mesh.GetPeriodicAlpha(); //DISP(alpha);
	int N = toInteger(round(2.0*pi_wp/abs(alpha))); //DISP(N);
        nb_modes_periodic_xyz(0) = N;
        nb_modes_periodic_xyz(1) = 1;
        if (per_z)
          {
            type_symmetry = VarBoundaryCondition_Base::PERIODIC_ZTHETA;            
            list_number_mode.Reallocate(N*Nz);
            list_number_mode.Fill();
          }
        else
          {
            if (Dimension::dim_N == 3)
              nb_modes_periodic_xyz(2) = 1;
            
            type_symmetry = VarBoundaryCondition_Base::PERIODIC_THETA;
            list_number_mode.Reallocate(N);
            list_number_mode.Fill();
          }
      }    
    else if (per_x)
      {
        if (per_y)
          {
            if (per_z)
              {
                type_symmetry = VarBoundaryCondition_Base::PERIODIC_XYZ;
                list_number_mode.Reallocate(Nx*Ny*Nz);
                list_number_mode.Fill();
              }
            else
              {
                type_symmetry = VarBoundaryCondition_Base::PERIODIC_XY;
                if (Dimension::dim_N == 3)
                  nb_modes_periodic_xyz(2) = 1;
                
                list_number_mode.Reallocate(Nx*Ny);
                list_number_mode.Fill();
              }
          }
        else if (per_z)
          {
            type_symmetry = VarBoundaryCondition_Base::PERIODIC_XZ;
            nb_modes_periodic_xyz(1) = 1;
            list_number_mode.Reallocate(Nx*Nz);
            list_number_mode.Fill();
          }
        else
          {
            type_symmetry = VarBoundaryCondition_Base::PERIODIC_X;
            nb_modes_periodic_xyz(1) = 1;
            if (Dimension::dim_N == 3)
              nb_modes_periodic_xyz(2) = 1;
            
            list_number_mode.Reallocate(Nx);
            list_number_mode.Fill();
          }
      }
    else if (per_y)
      {
        if (per_z)
          {
            type_symmetry = VarBoundaryCondition_Base::PERIODIC_YZ;
            nb_modes_periodic_xyz(0) = 1;
            list_number_mode.Reallocate(Ny*Nz);
            list_number_mode.Fill();
          }
        else
          {
            type_symmetry = VarBoundaryCondition_Base::PERIODIC_Y;
            nb_modes_periodic_xyz(0) = 1;
            if (Dimension::dim_N == 3)
              nb_modes_periodic_xyz(2) = 1;
            
            list_number_mode.Reallocate(Ny);
            list_number_mode.Fill();
          }
      }
    else if (per_z)
      {
        type_symmetry = VarBoundaryCondition_Base::PERIODIC_Z;
        nb_modes_periodic_xyz(0) = 1;
        nb_modes_periodic_xyz(1) = 1;
        list_number_mode.Reallocate(Nz);
        list_number_mode.Fill();
      }
    else
      {
        // no periodicity
        nb_modes_periodic_xyz.Fill(1);
      }
	
    if ((per_x) || (per_y) || (per_z) || (per_theta))
      {
	if (!var_comm.IsComplexProblem())
	  {
	    cout << "Cyclic domain implemented only for complex problems " << endl;
	    abort();
	  }
	
	if (var_comm.IsSymmetricProblem())
	  {
	    cout << "Cyclic domain induces non-symmetric matrix " << endl;
	    abort();
	  }
      }        
    
#ifdef SELDON_WITH_MPI
    // broadcasting modes to other processors
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
    if (nb_proc > 1)
      {
        MPI_Bcast(&type_symmetry, 1, MPI_INTEGER, 0, var_comm.comm_group_mode);
        MPI_Bcast(&nb_modes_periodic_xyz(0),
                  Dimension::dim_N, MPI_INTEGER, 0, var_comm.comm_group_mode);
        
        int nb = list_number_mode.GetM();
        MPI_Bcast(&nb, 1, MPI_INTEGER, 0, var_comm.comm_group_mode);
        if (rank_proc != 0)
          list_number_mode.Reallocate(nb);
        
        MPI_Bcast(list_number_mode.GetData(), nb, MPI_INTEGER, 0, var_comm.comm_group_mode);
      }
#endif
    
  }


  //! computation of phases for quasi-periodic conditions due to symmetry
  template<>
  void VarBoundaryCondition_Dim<Real_wp, Dimension2>
  ::ComputeQuasiPeriodicPhase()
  {
    // real numbers => nothing to do
  }
  

  //! computation of phases for quasi-periodic conditions due to symmetry  
  template<>
  void VarBoundaryCondition_Dim<Complex_wp, Dimension2>::ComputeQuasiPeriodicPhase()
  {
    Complex_wp phase, one(1, 0);
    Complex_wp kx, ky, kz; R_N vh;
    int n = this->number_mode;
    this->GetPeriodicModes(n, kx, ky, kz);
    MeshNumbering<Dimension2>& mesh_num = var_comm.GetMeshNumbering(0);
    if ((var_comm.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        && (var_comm.GetMeshNumbering(0).GetFormulationForPeriodicCondition()
            == MeshNumbering<Dimension2>::STRONG_PERIODIC))
      {
        periodic_dof_phase.Reallocate(mesh_num.GetNbPeriodicDof());
	periodic_dof_phase.Fill(one);
	mesh_num.InitPeriodicBoundary();
        for (int i = 0; i < mesh_num.GetNbPeriodicDof(); i++)
          {
            // translation vector (or angle for cyclic domains)
            vh = mesh_num.GetTranslationPeriodicDof(i);
            int type_dof = mesh_num.GetPeriodicityTypeForDof(i);
            Complex_wp arg(0, 0);
            switch (type_dof)
              {
              case BoundaryConditionEnum::PERIODIC_THETA :
                arg = kx*vh(0);
                break;
              case BoundaryConditionEnum::PERIODIC_CTE :
                arg = 0;
                break;
              case BoundaryConditionEnum::PERIODIC_X :
                arg = kx*vh(0);
                break;
              case BoundaryConditionEnum::PERIODIC_Y :
                arg = ky*vh(1);
                break;
              case BoundaryConditionEnum::PERIODIC_XY :
                arg = kx*vh(0) + ky*vh(1);
                break;
              }

            phase = exp(Iwp*arg);
            periodic_dof_phase(i) = phase;
          }        
      }
    else
      {
        periodic_dof_phase.Reallocate(var_comm.mesh.GetNbBoundaryRef());
	periodic_dof_phase.Fill(one);
        for (int i = 0; i < var_comm.mesh.GetNbBoundaryRef(); i++)
	  if (mesh_num.GetPeriodicBoundary(i) >= 0)
            {
              // translation vector (or angle for cyclic domains)
              vh = mesh_num.GetTranslationPeriodicBoundary(i);
              int type_dof = mesh_num.GetPeriodicityTypeForBoundary(i);
              Complex_wp arg(0, 0);
	      switch (type_dof)
                {
                case BoundaryConditionEnum::PERIODIC_THETA :
                  arg = kx*vh(0);
                  break;
                case BoundaryConditionEnum::PERIODIC_CTE :
                  arg = 0;
                  break;
                case BoundaryConditionEnum::PERIODIC_X :
                  arg = kx*vh(0);
                  break;
                case BoundaryConditionEnum::PERIODIC_Y :
                  arg = ky*vh(1);
                  break;
                case BoundaryConditionEnum::PERIODIC_XY :
                  arg = kx*vh(0) + ky*vh(1);
                  break;
                }
              
	      phase = exp(Iwp*arg);
              periodic_dof_phase(i) = phase;
	    }
      }
  }
  

  //! computation of phases for quasi-periodic conditions due to symmetry
  template<>
  void VarBoundaryCondition_Dim<Real_wp, Dimension3>::ComputeQuasiPeriodicPhase()
  {
    // real numbers => nothing to do
  }


  //! computation of phases for quasi-periodic conditions due to symmetry  
  template<>
  void VarBoundaryCondition_Dim<Complex_wp, Dimension3>::ComputeQuasiPeriodicPhase()
  {
#ifdef MONTJOIE_WITH_THREE_DIM
    Complex_wp phase, one(1, 0);
    Complex_wp kx, ky, kz; R_N vh;
    int n = this->number_mode;
    this->GetPeriodicModes(n, kx, ky, kz);
    MeshNumbering<Dimension3>& mesh_num = var_comm.GetMeshNumbering(0);
    if ((var_comm.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        && (mesh_num.GetFormulationForPeriodicCondition()
            == mesh_num.STRONG_PERIODIC))
      {
        periodic_dof_phase.Reallocate(mesh_num.GetNbPeriodicDof());
	periodic_dof_phase.Fill(one);
	mesh_num.InitPeriodicBoundary();
        for (int i = 0; i < mesh_num.GetNbPeriodicDof(); i++)
          {
            // translation vector (or angle for cyclic domains)
            vh = mesh_num.GetTranslationPeriodicDof(i);
            int type_dof = mesh_num.GetPeriodicityTypeForDof(i);
            Complex_wp arg(0, 0);
            switch (type_dof)
              {
              case BoundaryConditionEnum::PERIODIC_THETA :
                arg = kx*vh(0);
                break;
              case BoundaryConditionEnum::PERIODIC_CTE :
                arg = 0;
                break;
              case BoundaryConditionEnum::PERIODIC_X :
                arg = kx*vh(0);
                break;
              case BoundaryConditionEnum::PERIODIC_Y :
                arg = ky*vh(1);
                break;
              case BoundaryConditionEnum::PERIODIC_Z :
                arg = kz*vh(2);
                break;
              case BoundaryConditionEnum::PERIODIC_XY :
                arg = kx*vh(0) + ky*vh(1);
                break;
              case BoundaryConditionEnum::PERIODIC_XZ :
                arg = kx*vh(0) + kz*vh(2);
                break;
              case BoundaryConditionEnum::PERIODIC_YZ :
                arg = ky*vh(1) + kz*vh(2);
                break;
              case BoundaryConditionEnum::PERIODIC_XYZ :
                arg = kx*vh(0) + ky*vh(1) + kz*vh(2);
                break;
              case BoundaryConditionEnum::PERIODIC_ZTHETA :
                arg = kx*vh(0) + kz*vh(2);
                break;
              }
            
	    phase = exp(Iwp*arg);
	    periodic_dof_phase(i) = phase;
          }
        
      }
    else
      {
        periodic_dof_phase.Reallocate(var_comm.mesh.GetNbBoundaryRef());
	periodic_dof_phase.Fill(one);
        for (int i = 0; i < var_comm.mesh.GetNbBoundaryRef(); i++)
          if (mesh_num.GetPeriodicBoundary(i) >= 0)
            {
              // translation vector (or angle for cyclic domains)
              vh = mesh_num.GetTranslationPeriodicBoundary(i);
              int type_dof = mesh_num.GetPeriodicityTypeForBoundary(i);
              Complex_wp arg(0, 0);
              switch (type_dof)
                {
                case BoundaryConditionEnum::PERIODIC_THETA :
                  arg = kx*vh(0);
                  break;
                case BoundaryConditionEnum::PERIODIC_CTE :
                  arg = 0;
                  break;
                case BoundaryConditionEnum::PERIODIC_X :
                  arg = kx*vh(0);
                  break;
                case BoundaryConditionEnum::PERIODIC_Y :
                  arg = ky*vh(1);
                  break;
                case BoundaryConditionEnum::PERIODIC_Z :
                  arg = kz*vh(2);
                  break;
                case BoundaryConditionEnum::PERIODIC_XY :
                  arg = kx*vh(0) + ky*vh(1);
                  break;
                case BoundaryConditionEnum::PERIODIC_XZ :
                  arg = kx*vh(0) + kz*vh(2);
                  break;
                case BoundaryConditionEnum::PERIODIC_YZ :
                  arg = ky*vh(1) + kz*vh(2);
                  break;
                case BoundaryConditionEnum::PERIODIC_XYZ :
                  arg = kx*vh(0) + ky*vh(1) + kz*vh(2);
                  break;
                case BoundaryConditionEnum::PERIODIC_ZTHETA :
                  arg = kx*vh(0) + kz*vh(2);
                  break;
                }
              
	      phase = exp(Iwp*arg);
              periodic_dof_phase(i) = phase;
            }
      }
#endif
  }


  //! copies input parameters of another similar problem
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>::CopyInputData(const VarBoundaryCondition_Base& var_base)
  {
    VarBoundaryCondition_Base::CopyInputData(var_base);
    
    const VarBoundaryCondition_Dim<Complexe, Dimension>& var
      = static_cast<const VarBoundaryCondition_Dim<Complexe, Dimension>& >(var_base);

    point_dirichlet_condition = var.point_dirichlet_condition;
    Param_condition = var.Param_condition;
    Coef_Imped_ABC = var.Coef_Imped_ABC;
  }
  

  /************************
   * VarBoundaryCondition *
   ************************/
  

  //! treatment of boundary conditions for a wave-like equation
  /*!
    \param[in] mat_sp sparse matrix to modify
  */
  template<class Complexe, class Dimension> template<class T>
  void VarBoundaryCondition<Complexe, Dimension>
  ::AddBoundaryConditionGen(VirtualMatrix<T>& mat_sp, const GlobalGenericMatrix<T>& nat_mat,
			    int offset_row, int offset_col)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if ((rank_proc == 0)||(var_problem.print_level >= 10))
      if (var_problem.print_level >= 2)
	cout << rank_proc << " Treating impedance boundary conditions" << endl;
    
    T one; SetComplexOne(one);
    const IVect& ref_cond = var_problem.mesh.GetBoundaryCondition();
    
    // adding terms due to absorbing boundary condition
    ImpedanceFunction_Base<T, Dimension>* fct_imped_abc;
    fct_imped_abc = this->GetNewImpedanceABC(one);
    this->AddMatrixImpedanceBoundary(one, ref_cond, BoundaryConditionEnum::LINE_ABSORBING,
				     nat_mat, mat_sp, offset_row, offset_col, *fct_imped_abc);

    delete fct_imped_abc;
    
    // adding terms tue to impedance boundary condition
    ImpedanceFunction_Base<T, Dimension>* fct_imped_gen;
    fct_imped_gen = this->GetNewImpedanceGeneric(one);
    this->AddMatrixImpedanceBoundary(one, ref_cond, BoundaryConditionEnum::LINE_IMPEDANCE,
				     nat_mat, mat_sp, offset_row, offset_col, *fct_imped_gen);

    delete fct_imped_gen;
    
#ifdef MONTJOIE_WITH_TRANSMISSION
    // addings terms due to transmission conditions
    var_transmission_base.AddTransmissionTerms(one, nat_mat, mat_sp, 
					       offset_row, offset_col);
#endif

    // adding terms due to Generalized Impedance Boundary Conditions
    var_gibc_base.AddGibcTerms(one, nat_mat, mat_sp, offset_row, offset_col);
    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL    
    // adding terms due to high-conductivity boundary condition
    ImpedanceFunction_Base<T, Dimension>* fct_imped_conduc;
    fct_imped_conduc = this->GetNewImpedanceHighConductivity(one);
    this->AddMatrixImpedanceBoundary(one, ref_cond, BoundaryConditionEnum::LINE_HIGH_CONDUCTIVITY,
				     nat_mat, mat_sp, offset_row, offset_col, *fct_imped_conduc);
    
    delete fct_imped_conduc;
#endif

  }

  
  /****************************
   * VirtualImpedanceFunction *
   ****************************/
  
      
  //! sets attributes to default values
  template<class T, class Dimension>
  void VirtualImpedanceFunction<T, Dimension>::InitDefaultValues()
  {
    coef_phi.Reallocate(var_problem.nb_unknowns);
    coef_grad.Reallocate(var_problem.nb_unknowns);
    stored_coef_phi.Reallocate(var_problem.nb_unknowns);
    stored_coef_grad.Reallocate(var_problem.nb_unknowns);
    coef_phi.Zero(); coef_grad.Zero();
    read_param_condition = true;
    only_surface_dof_involved = true;
  }
  

  template<class T, class Dimension>
  bool VirtualImpedanceFunction<T, Dimension>::PresenceGradient() const
  {
    T zero; SetComplexZero(zero);
    for (int k = 0; k < coef_grad.GetM(); k++)
      if (coef_grad(k) != zero)
        return true;
    
    return false;
  }
  
 
  //! returns the impedance if it is a scalar coefficient
  template<class T, class Dimension>
  T VirtualImpedanceFunction<T, Dimension>
  ::GetCoefficient(int i, int num_elem,
		   int num_loc, int k, int ref_domain, int ref,
		   const SetPoints<Dimension>&, const SetMatrices<Dimension>&)
  {
    T imped;
    SetComplexZero(imped);
    
    if (read_param_condition)
      {
	SetComplexOne(imped);
	var_boundary.MltParamCondition(ref, 0, imped);
	var_boundary.MltMuIntegrationByParts(ref_domain, num_elem, num_loc, k, imped);
      }
    
    return imped;
  }
  

  //! Evaluates the impedance
  template<class T, class Dimension>
  void VirtualImpedanceFunction<T, Dimension>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			 const GlobalGenericMatrix<T>& nat_mat, int ref_domain,
			 const SetPoints<Dimension>& Pts,
			 const SetMatrices<Dimension>& Mat)
  {
    // uniform impedance, only the first quadrature point needs to be involved
    if (k > 0)
      return;
       
    // storing the normale
    vec_normale.Reallocate(Mat.GetNbPointsQuadratureBoundary());
    for (int kp = 0; kp < Mat.GetNbPointsQuadratureBoundary(); kp++)
      vec_normale(kp) = Mat.GetNormaleQuadratureBoundary(kp);
    
    int ref = var_problem.mesh.BoundaryRef(i).GetReference(); 
    T imped = this->GetCoefficient(i, num_elem, num_loc, k, ref_domain,
                                   ref, Pts, Mat);
    
    coef_phi.Fill(imped);
  }
  

  //! Evaluates the impedance for gradient terms
  template<class T, class Dimension>
  void VirtualImpedanceFunction<T, Dimension>
  ::EvaluateImpedanceGrad(int i, int num_elem, int num_edge, int num_loc, int k,
			  const GlobalGenericMatrix<T>& nat_mat, int ref_domain,
			  const SetPoints<Dimension>& Pts,
			  const SetMatrices<Dimension>& Mat)
  {
    // uniform impedance, only the first quadrature point needs to be involved
    if (k > 0)
      return;    
  }
  
  
  /**************************
   * ImpedanceFunction_Base *
   **************************/

  
  //! Applies impedance f(phi_boundary, grad_phi_boundary) for H^1 unknowns
  template<class T>
  void ImpedanceFunction_Base<T, Dimension2>::
  ApplyImpedancePhi_H1(int m, int j, int offset, const TinyVector<Real_wp, 1>& phi_boundary,
                       const R2& grad_phi_boundary, Vector<T>& f_phi)
  {
    T coef = this->coef_phi(m);
    if (this->stored_coef_phi.GetM() > 0)
      if (this->stored_coef_phi(m).GetM() > 0)
        coef = this->stored_coef_phi(m)(j);

    // identity matrix => only m component is modified
    f_phi(offset) = phi_boundary(0)*coef;
  }


  //! Applies impedance g(phi_boundary, grad_phi_boundary) for H^1 unknowns
  template<class T>
  void ImpedanceFunction_Base<T, Dimension2>::
  ApplyImpedanceGrad(int m, int j, int offset, const TinyVector<Real_wp, 1>& phi_boundary,
		     const R2& grad_phi_boundary, Vector<T>& g_phi)
  {
    if (this->read_param_condition)
      return;
    
    T coef = this->coef_grad(m);
    if (this->stored_coef_grad.GetM() > 0)
      if (this->stored_coef_grad(m).GetM() > 0)
        coef = this->stored_coef_grad(m)(j);
    
    // extracting tangential component of grad u
    T u_dot_n;
    TinyVector<T, 2> vec_u;
    u_dot_n = DotProd(grad_phi_boundary, this->vec_normale(j));
    vec_u = grad_phi_boundary;
    Add(-u_dot_n, this->vec_normale(j), vec_u);
    
    // then storing f = coef_grad grad_tangential phi on the m component
    vec_u *= coef;
    g_phi(offset) = vec_u(0);
    g_phi(offset+1) = vec_u(1);
  }


  //! Applies impedance f(phi_boundary, curl_phi_boundary) for H(curl) unknowns
  template<class T>
  void ImpedanceFunction_Base<T, Dimension2>::
  ApplyImpedancePhi_Hcurl(int m, int j, int offset, const TinyVector<Real_wp, 2>& phi_boundary,
                          const TinyVector<Real_wp, 1>& curl_phi_boundary, Vector<T>& f_phi)
  {
    T coef = this->coef_phi(m);
    if (this->stored_coef_phi.GetM() > 0)
      if (this->stored_coef_phi(m).GetM() > 0)
        coef = this->stored_coef_phi(m)(j);

    // extracting tangential component of u
    Real_wp u_dot_n;
    TinyVector<T, 2> vec_u;
    u_dot_n = DotProd(phi_boundary, this->vec_normale(j));
    vec_u(0) = phi_boundary(0) - u_dot_n*this->vec_normale(j)(0);
    vec_u(1) = phi_boundary(1) - u_dot_n*this->vec_normale(j)(1);
    
    // then storing f = coef_phi u_tangential phi on the m component
    vec_u *= coef;
    f_phi(offset) = vec_u(0);
    f_phi(offset+1) = vec_u(1);
  }


  //! Applies impedance g(phi_boundary, curl_phi_boundary) for H(curl) unknowns
  template<class T>
  void ImpedanceFunction_Base<T, Dimension2>::
  ApplyImpedanceCurl(int m, int j, int offset, const TinyVector<Real_wp, 2>& phi_boundary,
		     const TinyVector<Real_wp, 1>& curl_phi_boundary, Vector<T>& g_phi)
  {
    T coef = this->coef_grad(m);
    if (this->stored_coef_grad.GetM() > 0)
      if (this->stored_coef_grad(m).GetM() > 0)
        coef = this->stored_coef_grad(m)(j);

    // identity matrix => only m component is modified
    g_phi(offset) = curl_phi_boundary(0)*coef;
  }


  //! Applies impedance f(phi_boundary, div_phi_boundary) for H(div)
  template<class T>
  void ImpedanceFunction_Base<T, Dimension2>::
  ApplyImpedancePhi_Hdiv(int m, int j, int offset, const R2& phi_boundary,
                         const TinyVector<Real_wp, 1>& div_phi_boundary, Vector<T>& f_phi)
  {
    T coef = this->coef_phi(m);
    if (this->stored_coef_phi.GetM() > 0)
      if (this->stored_coef_phi(m).GetM() > 0)
        coef = this->stored_coef_phi(m)(j);
    
    // computing normal component of u
    T u_dot_n;
    TinyVector<T, 2> vec_u;
    u_dot_n = DotProd(phi_boundary, this->vec_normale(j));
    vec_u = this->vec_normale(j);
    
    // then storing f = coef_phi phi \cdot n  n  on the m component
    vec_u *= coef*u_dot_n;
    f_phi(offset) = vec_u(0);
    f_phi(offset+1) = vec_u(1);
  }


  //! Applies impedance g(phi_boundary, div_phi_boundary) for H(div)
  template<class T>
  void ImpedanceFunction_Base<T, Dimension2>::
  ApplyImpedanceDiv(int m, int j, int offset, const R2& phi_boundary,
                    const TinyVector<Real_wp, 1>& div_phi_boundary, Vector<T>& g_phi)
  {
    if (this->read_param_condition)
      return;

    T coef = this->coef_grad(m);
    if (this->stored_coef_grad.GetM() > 0)
      if (this->stored_coef_grad(m).GetM() > 0)
        coef = this->stored_coef_grad(m)(j);
    
    // identity matrix => only m component is modified
    g_phi(offset) = div_phi_boundary(0)*coef;    
  }


  /// 3D ///

  
  //! Applies impedance f(phi_boundary, grad_phi_boundary) for H^1
  template<class T>
  void ImpedanceFunction_Base<T, Dimension3>::
  ApplyImpedancePhi_H1(int m, int j, int offset, const TinyVector<Real_wp, 1>& phi_boundary,
                       const R3& grad_phi_boundary, Vector<T>& f_phi)
  {
    T coef = this->coef_phi(m);
    if (this->stored_coef_phi.GetM() > 0)
      if (this->stored_coef_phi(m).GetM() > 0)
        coef = this->stored_coef_phi(m)(j);

    // identity matrix => only m component is modified
    f_phi(offset) = phi_boundary(0)*coef;
  }


  //! Applies impedance g(phi_boundary, grad_phi_boundary) for H^1
  template<class T>
  void ImpedanceFunction_Base<T, Dimension3>::
  ApplyImpedanceGrad(int m, int j, int offset, const TinyVector<Real_wp, 1>& phi_boundary,
		     const R3& grad_phi_boundary, Vector<T>& g_phi)
  {
    if (this->read_param_condition)
      return;
    
    T coef = this->coef_grad(m);
    if (this->stored_coef_grad.GetM() > 0)
      if (this->stored_coef_grad(m).GetM() > 0)
        coef = this->stored_coef_grad(m)(j);

    // extracting tangential component of grad u
    T u_dot_n;
    TinyVector<T, 3> vec_u;
    u_dot_n = DotProd(grad_phi_boundary, this->vec_normale(j));
    vec_u = grad_phi_boundary;
    Add(-u_dot_n, this->vec_normale(j), vec_u);
    
    // then storing f = coef_grad grad_tangential phi on the m component
    vec_u *= coef;
    g_phi(offset) = vec_u(0);
    g_phi(offset+1) = vec_u(1);
    g_phi(offset+2) = vec_u(2);
  }

  
  //! Applies impedance f(phi_boundary, curl_phi_boundary) for H(curl)
  template<class T>
  void ImpedanceFunction_Base<T, Dimension3>::
  ApplyImpedancePhi_Hcurl(int m, int j, int offset, const TinyVector<Real_wp, 3>& phi_boundary,
                          const TinyVector<Real_wp, 3>& curl_phi_boundary, Vector<T>& f_phi)
  {
    T coef = this->coef_phi(m);
    if (this->stored_coef_phi.GetM() > 0)
      if (this->stored_coef_phi(m).GetM() > 0)
        coef = this->stored_coef_phi(m)(j);

    // extracting tangential component of u
    T u_dot_n;
    TinyVector<T, 3> vec_u;
    u_dot_n = DotProd(phi_boundary, this->vec_normale(j));
    vec_u(0) = phi_boundary(0) - u_dot_n*this->vec_normale(j)(0);
    vec_u(1) = phi_boundary(1) - u_dot_n*this->vec_normale(j)(1);
    vec_u(2) = phi_boundary(2) - u_dot_n*this->vec_normale(j)(2);
    
    // then storing f = coef_phi u_tangential phi on the m component
    vec_u *= coef;
    f_phi(offset) = vec_u(0);
    f_phi(offset+1) = vec_u(1);
    f_phi(offset+2) = vec_u(2);
  }


  //! Applies impedance g(phi_boundary, curl_phi_boundary) for H(curl)
  template<class T>
  void ImpedanceFunction_Base<T, Dimension3>::
  ApplyImpedanceCurl(int m, int j, int offset, const TinyVector<Real_wp, 3>& phi_boundary,
		     const TinyVector<Real_wp, 3>& curl_phi_boundary, Vector<T>& g_phi)
  {
    T coef = this->coef_grad(m);
    if (this->stored_coef_grad.GetM() > 0)
      if (this->stored_coef_grad(m).GetM() > 0)
        coef = this->stored_coef_grad(m)(j);

    // term n \cdot curl phi_i  n \cdot curl phi_j
    T curlU_dot_n = DotProd(curl_phi_boundary, this->vec_normale(j));
    TinyVector<T, 3> vec_u;
    vec_u(0) = curlU_dot_n*this->vec_normale(j)(0);
    vec_u(1) = curlU_dot_n*this->vec_normale(j)(1);
    vec_u(2) = curlU_dot_n*this->vec_normale(j)(2);
    
    // storing n \cdot curl phi_i  n on the m component
    vec_u *= coef;
    g_phi(offset) = vec_u(0);
    g_phi(offset+1) = vec_u(1);
    g_phi(offset+2) = vec_u(2);
  }


  //! Applies impedance f(phi_boundary, div_phi_boundary) for H(div)
  template<class T>
  void ImpedanceFunction_Base<T, Dimension3>::
  ApplyImpedancePhi_Hdiv(int m, int j, int offset, const R3& phi_boundary,
                         const TinyVector<Real_wp, 1>& div_phi_boundary, Vector<T>& f_phi)
  {
    T coef = this->coef_phi(m);
    if (this->stored_coef_phi.GetM() > 0)
      if (this->stored_coef_phi(m).GetM() > 0)
        coef = this->stored_coef_phi(m)(j);
    
    // computing normal component of u
    T u_dot_n;
    TinyVector<T, 3> vec_u;
    u_dot_n = DotProd(phi_boundary, this->vec_normale(j));
    vec_u = this->vec_normale(j);
    
    // then storing f = coef_phi phi \cdot n  n  on the m component
    vec_u *= coef*u_dot_n;
    f_phi(offset) = vec_u(0);
    f_phi(offset+1) = vec_u(1);
    f_phi(offset+2) = vec_u(2);
  }


  //! Applies impedance g(phi_boundary, div_phi_boundary) for H(div)
  template<class T>
  void ImpedanceFunction_Base<T, Dimension3>::
  ApplyImpedanceDiv(int m, int j, int offset, const R3& phi_boundary,
                    const TinyVector<Real_wp, 1>& div_phi_boundary, Vector<T>& g_phi)
  {
    if (this->read_param_condition)
      return;

    T coef = this->coef_grad(m);
    if (this->stored_coef_grad.GetM() > 0)
      if (this->stored_coef_grad(m).GetM() > 0)
        coef = this->stored_coef_grad(m)(j);
    
    // identity matrix => only m component is modified
    g_phi(offset) = div_phi_boundary(0)*coef;    
  }
  
}

#define MONTJOIE_FILE_BOUNDARY_CONDITION_HARMONIC_CXX
#endif
