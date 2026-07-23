#ifndef MONTJOIE_FILE_AXISYM_HELMHOLTZ_CXX

namespace Montjoie
{
  /************************
   * HelmholtzEquationAxi *
   ************************/


  void HelmholtzEquationAxi
  ::ComputeMassMatrix(EllipticProblem<HelmholtzEquationAxi>& var,
		      int i, const ElementReference_Dim<Dimension2>& Fb0)
  {
    if (var.use_iterative_solver)
      {
	const ElementReference<Dimension2, 1>& Fb = var.GetReferenceElementH1(i);
	int Ni = Fb.GetNbPointsQuadratureInside();
	var.Glob_mass_omega2(i).Reallocate(Ni);
	var.Glob_mass_iomega(i).Reallocate(Ni);
	var.Glob_mass_m2(i).Reallocate(Ni);

	if (var.include_flow_term)
	  {
	    var.Glob_grad_iomega(i).Reallocate(Ni);
	    var.Glob_mass_m_omega(i).Reallocate(Ni);
	  }

	if (var.include_double_gradM_flow)
	  var.Glob_grad_im(i).Reallocate(Ni);
	
	var.Glob_stiff(i).Reallocate(Ni);
	
	bool affine = var.mesh.IsElementAffine(i);
	int ref = var.mesh.Element(i).GetReference();
	Real_wp jacob; Matrix2_2 dfjm1;

	Complexe mu_teta;
	TinyMatrix<Complexe, Symmetric, 2, 2> mu_rz;
	TinyMatrix<Complexe, Symmetric, 3, 3> mu;
	TinyMatrix<Complexe, General, 2, 2> B, A;
	TinyVector<Complexe, 3> v0;
	TinyVector<Complexe, 2> vec_u, v_rz;
	
	for (int j = 0; j < Ni; j++)
	  {
	    if (affine)
	      {
		jacob = var.Glob_jacobian(i)(0)*Fb.WeightsND(j);
		dfjm1 = var.Glob_DFjm1(i)(0);                
		Mlt(1.0/var.Glob_jacobian(i)(0), dfjm1);
	      }
	    else
	      {
		jacob = var.Glob_jacobian(i)(j);
		dfjm1 = var.Glob_DFjm1(i)(j);
		Mlt(Fb.WeightsND(j)/jacob, dfjm1);
	      }

	    Complexe rho = var.ref_rho(ref).GetCoefficient(var, i, j);
	    mu = var.ref_mu(ref).GetCoefficient(var, i, j);	    
	    Complexe sigma = var.ref_sigma(ref).GetCoefficient(var, i, j);
	    if (var.include_flow_term)
	      {
		v0 = var.ref_v0(ref).GetCoefficient(var, i, j);        
		v_rz.Init(v0(0), v0(2));
	      }
	    
	    // mu_rz and mu_teta are extracted
	    mu_teta = mu(2, 2);
	    mu_rz(0, 0) = mu(0, 0);
	    mu_rz(0, 1) = mu(0, 1);
	    mu_rz(1, 1) = mu(1, 1);
	    
	    Complexe rtilde = var.Glob_rtilde(i)(j), beta; 
	    Complexe invRtilde; SetComplexZero(invRtilde);
	    if (abs(rtilde) > epsilon_machine*var.mesh.GetXmax())
	      invRtilde = 1.0/rtilde;
	    
	    if (var.include_double_gradM_flow)
	      {
		// to mu_rz, we subtract [m_r^2, m_r m_z; m_r m_z, m_z^2]
		// to mu_teta, we sutract -m_theta^2 		
		beta = var.ref_beta(ref).GetCoefficient(var, i, j);
		Complexe m_r = v0(0), m_z = v0(2);
		mu_rz(0, 0) -= beta*m_r*m_r;
		mu_rz(0, 1) -= beta*m_r*m_z;
		mu_rz(1, 1) -= beta*m_z*m_z;
		mu_teta -= beta*v0(1)*v0(1);

		// to v, we add m_theta m / omega [m_r;m_z]
		//v0(0) += beta*v0(1)*m_*m_r*invRtilde / var_problem.GetOmega();
		//v0(2) += beta*v0(1)*m_*m_z*invRtilde / var_problem.GetOmega();
	      }
	    
	    Mlt(dfjm1, mu_rz, B);
	    MltTrans(B, dfjm1, A);
	    
	    // term r \rho J_i \omega_k
	    var.Glob_mass_omega2(i)(j) = rtilde*jacob*rho;

	    // term r sigma J_i omega_k
	    var.Glob_mass_iomega(i)(j) = rtilde*jacob*sigma;

	    // term r DF_i^{-1} C DF_i^{-T} J_i \omega_k
	    // where C = mu_rz - \beta v0 v0^T 
	    var.Glob_stiff(i)(j) = A;
	    Mlt(rtilde*jacob, var.Glob_stiff(i)(j));
	    
	    // term (mu_teta - \beta v0_\theta^2) J_i \omega_k / r
	    var.Glob_mass_m2(i)(j) = mu_teta * jacob * invRtilde;

	    if (var.include_flow_term)
	      {
		Mlt(dfjm1, v_rz, vec_u);
		// term r DF_i^{-1} M_{rz} J_i \omega_k
		var.Glob_grad_iomega(i)(j) = rtilde*jacob*vec_u;
		// term - 2 M_\theta J_i \omega_k
		var.Glob_mass_m_omega(i)(j) = -2*v0(1)*jacob;
	      }

	    if (var.include_double_gradM_flow)
	      {
		var.Glob_grad_im(i)(j) = beta*v0(1)*jacob*vec_u;
	      }
	  }
      }
  }
  

  /********************
   * VarHelmholtz_Axi *
   ********************/

  
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::InitDefaultValues()
  {
    calcul_enveloppe = false;
    formulation_r3 = false;
    include_flow_term = false;
    include_double_gradM_flow = false;
    
    this->rho0 = 1.0;
    this->mu0 = 1.0;
    this->sigma0 = 0.0;
    
    use_iterative_solver = false;
    same_pattern_LU = false;
  }
  
  
  //! additionnal input datas for axisymmetric computation
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarAxisymProblem::SetInputData(description_field, parameters);
    
    if (!description_field.compare("CalculEnveloppe"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_Axi" << endl;
	    cout << "CalculEnveloppe needs 1 parameters, for instance :" << endl;
	    cout << "CalculEnveloppe = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        if (parameters(0) == "YES")
          calcul_enveloppe = true;
        else
          calcul_enveloppe = false;
      }
    else if (!description_field.compare("FormulationAxisymmetric"))
      {
        if (parameters(0) == "R3")
          formulation_r3 = true;
        else
          formulation_r3 = false;
      }
    else if (!description_field.compare("AddFlowTerm"))
      {
	if (parameters(0) == "YES")
	  include_flow_term = true;
	else if (parameters(0) == "GRAD")
	  {
	    include_flow_term = true;
	    include_double_gradM_flow = true;
	  }
	else
	  include_flow_term = false;
      }
    else if (!description_field.compare("UseIterativeSolverM"))
      {
	if (parameters(0) == "YES")
	  use_iterative_solver = true;
	else
	  use_iterative_solver = false;
      }
    else if (!description_field.compare("UseSamePatternM"))
      {
	if (parameters(0) == "YES")
	  same_pattern_LU = true;
	else
	  same_pattern_LU = false;	
      }
  }
  
  
  //! return true if additional unknown is required for absorbing boundary condition
  template<class Complexe>
  bool VarHelmholtz_Axi<Complexe>::UseAdditionalUnknownForABC() const
  {
    if (var_boundary.grazing_abc)
      {        
        if (var_boundary.GetOrderAbsorbingCondition() == 2)
          return true;
        
        if (var_boundary.GetOrderAbsorbingCondition() >= 4)
          return true;
      }
    
    if (var_boundary.GetOrderAbsorbingCondition() == 6)
      return true;

    //if (var_boundary.GetOrderAbsorbingCondition() == 23)
    //return true;
    
    return false;
  }

  
  //! la matrice est-elle symetrique ?
  template<class Complexe>
  bool VarHelmholtz_Axi<Complexe>::IsSymmetricProblem(bool eigen) const
  {
    if (calcul_enveloppe)
      return false;
    
    if (include_flow_term)
      return false;
  
    if (this->UseAdditionalUnknownForABC() && var_gibc_base.UseUnsymmetricImplementation())
      return false;
  
    return true;
  }
  

  //! we check if the mesh is on the left side of the revolution axis Oz or on the right
  /*!
    If it is on the left, the mesh is symmetrized so that x is always positive
    if the mesh has elements on the left and elements on the right, it is not a good mesh...
  */
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::ConstructFiniteElement(const string& name_elt)
  {
    // For Radau points, we need to construct Lobatto elements
    // and affect these elements to elements outside the axis
    string name_element = name_elt;
    if (name_elt == "TRIANGLE_RADAU")
      name_element = "TRIANGLE_LOBATTO";
    
    // orders present in the mesh
    TinyVector<IVect, 4> order;    
    var_problem.GetMeshNumbering(0).GetOrder(order);
    
    Vector<bool> change_elt(var_problem.mesh.GetNbElt());
    change_elt.Fill(true);
    if (this->force_diagonal_mass_invertible)
      for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
	if (this->Element_On_Axe(i))
	  change_elt(i) = false;
    
    var_problem.AddFiniteElement(name_element, order, change_elt, var_problem.mesh,
				 var_problem.GetMeshNumbering(0), var_problem.GetOverIntegration(),
				 var_problem.dg_formulation);
    
    // signaling that the geometry is axisymmetric
    var_problem.SetAxisymGeometry(true);
  }
  
    
  //! constructing arrays needed for the simulation
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::ConstructAll(const string& name_file, const string& name_elt,
						Vector<string>& lines_data_file)
  {            
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
    ReadLinesFile(name_file, lines_data_file, var_problem.comm_group_mode);
#else
    int rank_proc(0), nb_proc(1);
    ReadLinesFile(name_file, lines_data_file);
#endif
    
    All_LinearSolver* solver;
    ConstructAll(lines_data_file, name_file, name_elt, "none", solver);
  }


  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>
  ::ConstructAll(const Vector<string>& lines_data_file, const string& input_file,
                 const string& name_elt, const string& name_eq, All_LinearSolver*& glob_solver,
                 bool compute_rho, bool delete_points, int num)
  {
    this->InitIndices(50);
    var_problem.SetTypeEquation("none");
        
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    ReadInputFile(lines_data_file, var_output);
   
    glob_solver = this->GetNewLinearSolver();
    ReadInputFile(lines_data_file, *glob_solver);
 
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (rank_proc == 0)
      if (var_output.DOSSIER_output.size() == 0)
	{
	  int num = -1;
	  var_output.DOSSIER_output = "[STIFFOUT]/"; 
	  EcritDossier(var_output.DOSSIER_output, input_file, num);
        }
    
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      MPI_Bcast_string(var_output.DOSSIER_output, 0, var_problem.comm_group_mode);
#endif    
    
    // we read mesh and construct reference element (for finite element method)
    bool split_mesh = true;
    if (var_problem.GetNbProcPerMode() == 1)
      split_mesh = false;
    
    var_problem.ComputeMeshAndFiniteElement(name_elt, split_mesh);
    
    if (nb_proc == 1)
      var_problem.mesh.Write("test.mesh");
    
    var_problem.PerformOtherInitializations();

    // on ecrit tous les points nodaux pour lesquels il faut calculer les indices
    if (var_problem.print_level >= 8)
      var_problem.WriteNodalPointsMesh();
  }
  
  
  //! writes the index on the disk
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::WriteAllIndices()
  {
    // already done in VarHarmonic
  }
  
  
  //! performs adimensionalization of input parameters
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>
  ::PerformAdimensionalization(const Real_wp& t0_adim, const Real_wp& z0_adim)
  {
    // adimensionalisation of mesh and output grids
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
	// mesh
	int n0 = 2;
	if (var_problem.mesh_data(0)(0) == "REGULAR")
	  n0 = 2;
	else if (var_problem.mesh_data(0)(0) == "REGULAR_ANISO")
	  n0 = 3;
	else
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
	
	for (int i = n0; i < n0+4; i++)
	  var_problem.mesh_data(0)(i) = to_str(to_num<Real_wp>(var_problem.mesh_data(0)(i))/z0_adim);
	
        // right hand side
        for (int i = 0; i < var_source.GetNbParameterSource(0); i++)
          if (var_source.GetParameterSource(0, i).GetM() > 2)
            {
              Real_wp waist = to_num<Real_wp>(var_source.GetParameterSource(0, i)(2));
	      VectString& source_param = const_cast<VectString&>(var_source.GetParameterSource(0, i));
              source_param(2) = to_str(waist/z0_adim);
            }
        
	// output grids
        for (int i = 0; i < var_output.output_grid_param.GetM(); i++)
          var_output.output_grid_param(i).AdimTime(t0_adim);
        
        for (int i = 0; i < var_output.output_mesh_param.GetM(); i++)
	  var_output.output_mesh_param(i).AdimTime(t0_adim);
        
        for (int i = 0; i < var_output.var_grid.GetM(); i++)
	  var_output.var_grid(i).AdimSpace(z0_adim);
	
	// thickness of PML
        for (int num_pml = 0; num_pml < var_problem.mesh.GetNbPmlAreas(); num_pml++)
	  {
	    Real_wp delta = var_problem.mesh.GetPmlArea(num_pml).GetThicknessPML()/z0_adim;
	    var_problem.mesh.GetPmlArea(num_pml).SetThicknessPML(delta);
	  }
      }
  }
  
  
  //! allocation of arrays containing physical properties
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::InitIndices(int n)
  {
    Complexe rho0, mu0;
    SetComplexOne(rho0); SetComplexOne(mu0); 
    ref_rho.Reallocate(n+1);
    ref_mu.Reallocate(n+1);
    ref_sigma.Reallocate(n+1);
    ref_v0.Reallocate(n+1);
    ref_beta.Reallocate(n+1);
    for (int i = 0; i <= n; i++)
      {
	ref_rho(i).SetIdentity(); 
	ref_rho(i).Mlt(rho0);
	ref_mu(i).SetIdentity();
	ref_mu(i).Mlt(mu0);
	ref_sigma(i).Zero();
	ref_v0(i).Zero();
	ref_beta(i).SetIdentity(); 
      }
  }

  
  //! sets physical properties of media with reference i
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

    ref_rho(i).SetInputData(nb, parameters, parameters(0));
    ref_mu(i).SetInputData(nb, parameters, parameters(0));
    ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    if (include_flow_term)
      ref_v0(i).SetInputData(nb, parameters, parameters(0));

    if (include_double_gradM_flow)
      ref_beta(i).SetInputData(nb, parameters, parameters(0));
  }
  

  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

    if (name_media == "rho")
      ref_rho(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "mu")
      ref_mu(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "sigma")
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "M")
      {
        if (include_flow_term)
          ref_v0(i).SetInputData(nb, parameters, parameters(0));
      }
    else if (name_media == "beta")
      {
        if (include_double_gradM_flow)
          ref_beta(i).SetInputData(nb, parameters, parameters(0));
      }
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }
  

  //! returns the name associated with the physical index num
  template<class Complexe>
  string VarHelmholtz_Axi<Complexe>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("rho");
      case 1: return string("mu");
      case 2: return string("sigma");
      case 3: return string("M");
      case 4: return string("beta");
      }

    return string();
  }  
  
    
  //! fills arrays with the varying indices present in the current problem
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complexe>* >& rho_complex, IVect& num_ref,
		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    int nb = 0;
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        nb += ref_rho(i).GetNbVaryingMedia();
        nb += ref_mu(i).GetNbVaryingMedia();
        nb += ref_sigma(i).GetNbVaryingMedia();
	if (include_flow_term)
	  nb += ref_v0(i).GetNbVaryingMedia();

	if (include_double_gradM_flow)
	  nb += ref_beta(i).GetNbVaryingMedia();
      }
    
    rho_complex.Reallocate(nb);
    num_ref.Reallocate(nb);
    num_index.Reallocate(nb);
    num_component.Reallocate(nb);
    compute_grad.Reallocate(nb);
    compute_hess.Reallocate(nb);
    compute_grad.Fill(false);
    compute_hess.Fill(false);
    nb = 0;
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        int nb0 = nb;
        ref_rho(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        ref_mu(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 1;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        ref_sigma(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 2;
            num_ref(j) = i;
          }

	if (include_flow_term)
	  {
	    nb0 = nb;
	    ref_v0(i).GetVaryingMedia(nb, rho_complex, num_component);
	    for (int j = nb0; j < nb; j++)
	      {
		num_index(j) = 3;
		num_ref(j) = i;
	      }
	  }

	if (include_double_gradM_flow)
	  {
	    nb0 = nb;
	    ref_beta(i).GetVaryingMedia(nb, rho_complex, num_component);
	    for (int j = nb0; j < nb; j++)
	      {
		num_index(j) = 4;
		num_ref(j) = i;
	      }
	  }
      }
  }
  

  //! returns true if physical properties of media are varying
  template<class Complexe>
  bool VarHelmholtz_Axi<Complexe>::IsVaryingMedia(int i) const
  {
    if (ref_rho(i).IsVarying())
      return true;
    
    if (ref_sigma(i).IsVarying())
      return true;

    if (ref_mu(i).IsVarying())
      return true;

    if (include_flow_term)
      if (ref_v0(i).IsVarying())
        return true;

    if (include_double_gradM_flow)
      if (ref_beta(i).IsVarying())
	return true;
    
    return false;
  }
  
  
  //! finds values of rho0, mu0, sigma0
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>
  ::SetPhysicalIndexAtInfinity(const Vector<bool>& RefUsed)
  {
    // taking the maximum among elements of the current processor
    Vector<Real_wp> coef0(3);
    coef0.Fill(-1);
    for (int ref = 0; ref < RefUsed.GetM(); ref++)
      if (RefUsed(ref))
        {
          Real_wp rho = abs(ref_rho(ref).GetConstant());
          if (rho > coef0(0))
            coef0(0) = rho;

          Real_wp mu = abs(ref_mu(ref).GetConstant()(0, 0));
          if (mu > coef0(1))
            coef0(1) = mu;

          Real_wp sigma = abs(ref_sigma(ref).GetConstant());
          if (sigma > coef0(2))
            coef0(2) = sigma;
        }
    
    // then taking the maximum among processors
#ifdef SELDON_WITH_MPI
    Vector<Real_wp> max_coef0(3);
    Vector<int64_t> xtmp;
    MpiAllreduce(var_problem.comm_group_mode, coef0, xtmp, max_coef0, 3, MPI_MAX);    
    
    coef0 = max_coef0;
#endif
    
    this->rho0 = coef0(0); this->mu0 = coef0(1); this->sigma0 = coef0(2);
    
    if (coef0(0) == -1.0)
      this->rho0 = 1.0;

    if (coef0(1) == -1.0)
      this->mu0 = 1.0;

    if (coef0(2) == -1.0)
      this->sigma0 = 0.0;
    
    Complexe m_iomega; var_problem.GetMiomega(m_iomega); 
    Complexe rho0_tilde = this->rho0 + this->sigma0/m_iomega;
    kwave_envelope = var_problem.GetOmega()*sqrt(rho0_tilde/this->mu0);
  }
  
  
  //! affectation of Dirichlet dofs
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::UpdateDirichlet(int n)
  {
    if (this->formulation_r3)
      return;
    
    // for modes different from 0, we force u = 0 on the axis
    if (var_boundary.GetCurrentModeNumber() != 0)
      {
        for (int i = 0; i < Dof_On_Axe.GetM(); i++)
          var_boundary.SetDirichletDof(this->Dof_On_Axe(i), true);        
      }
    else
      for (int i = 0; i < Dof_On_Axe.GetM(); i++)
        var_boundary.SetDirichletDof(this->Dof_On_Axe(i), false);
    var_boundary.UpdateDirichletDofs();
    // mass lumping not needed for modes different from 0 and usual formulation
    if ((n > 0) && (this->force_diagonal_mass_invertible))
      {
        this->force_diagonal_mass_invertible = false;
	// extend TriangleLobatto for all the elements
	Vector<bool> change_elt(var_problem.mesh.GetNbElt());
	change_elt.Fill(true);
	var_problem.SetFiniteElement("TRIANGLE_LOBATTO", change_elt, 1);
        var_problem.ComputeMassMatrix();
      }
  }
  
  
  //! allocation of arrays needed for the computation of elementary matrices
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::AllocateMassMatrices()
  {
    Glob_rtilde.Reallocate(var_problem.mesh.GetNbElt());
    Glob_radius.Reallocate(var_problem.mesh.GetNbElt());
    Glob_DrPML.Reallocate(var_problem.mesh.GetNbElt());
    Glob_DzPML.Reallocate(var_problem.mesh.GetNbElt());

    this->coefficient_impedance_absorbing.Reallocate(this->GetNbPhysicalIndices());
    this->coefficient_impedance_absorbing.Fill(1.0);
    for (int i = 0; i < this->coefficient_impedance_absorbing.GetM(); i++)
      {
	Complexe rho = this->ref_rho(i).GetConstant(), mu = this->ref_mu(i).GetConstant()(0, 0);
        Complexe sigma = this->ref_sigma(i).GetConstant();
        Complexe rho_tilde;
	to_complex(rho + Iwp*sigma/var_problem.GetOmega(), rho_tilde);
        this->coefficient_impedance_absorbing(i) = sqrt(rho_tilde*mu);
      }

    if (use_iterative_solver)
      {
	if ((var_problem.GetBoundaryConditionProblem().GetNbGlobalEltPML() > 0)
	    || (this->formulation_r3) || (this->calcul_enveloppe))
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
	
	Glob_mass_omega2.Reallocate(var_problem.mesh.GetNbElt());
	Glob_mass_iomega.Reallocate(var_problem.mesh.GetNbElt());
	Glob_mass_m2.Reallocate(var_problem.mesh.GetNbElt());
	if (this->include_flow_term)
	  {
	    Glob_grad_iomega.Reallocate(var_problem.mesh.GetNbElt());
	    Glob_mass_m_omega.Reallocate(var_problem.mesh.GetNbElt());
	  }
	
	if (this->include_double_gradM_flow)
	  this->Glob_grad_im.Reallocate(var_problem.mesh.GetNbElt());
	
	Glob_stiff.Reallocate(var_problem.mesh.GetNbElt());
      }
  }
  
  
  //! computation of geometric quantities needed for the expression of elementary matrices
  template<class Complexe>
  void VarHelmholtz_Axi<Complexe>::ComputeLocalMassMatrix(int i)
  {
    int N = var_problem.Glob_PointsQuadrature(i).GetM();
    Glob_rtilde(i).Reallocate(N);
    if (!var_problem.InsidePML(i))
      {
        for (int j = 0; j < N; j++)
          Glob_rtilde(i)(j) = var_problem.Glob_PointsQuadrature(i)(j)(0);
      }
    else
      {
        Glob_radius(i).Reallocate(N);
        Glob_DrPML(i).Reallocate(N);
        Glob_DzPML(i).Reallocate(N);
        
        Complexe rtilde, dr, dz;
        int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
        for (int j = 0; j < N; j++)
          {
            dr = 1.0/var_boundary.GetTauPML(i1, j, 0);
            dz = 1.0/var_boundary.GetTauPML(i1, j, 1);
            
            rtilde = var_boundary.GetPrimitiveTauPML(i1, j, 0);
            Real_wp radius = var_problem.Glob_PointsQuadrature(i)(j)(0);
            
            Glob_radius(i)(j) = radius;
            Glob_rtilde(i)(j) = rtilde;
            Glob_DrPML(i)(j) = dr;
            Glob_DzPML(i)(j) = dz;
          }
      }        
  }


  template<>
  void VarHelmholtz_Axi<Real_wp>
  ::ModifyCoefficientPML(bool form_r3, TinyMatrix<Real_wp, Symmetric, 2, 2>& muPML,
			 Real_wp& rtilde, Real_wp& invRtilde, Real_wp& coef,
			 Real_wp& coef_cte, Real_wp& coef_grad, TinyVector<Real_wp, 3>& v0,
			 int iquad, int i1, int j)
  {
    muPML.SetIdentity();
  }
  

  template<>
  void VarHelmholtz_Axi<Complex_wp>
  ::ModifyCoefficientPML(bool form_r3, TinyMatrix<Complex_wp, Symmetric, 2, 2>& muPML,
			 Complex_wp& rtilde, Complex_wp& invRtilde, Complex_wp& coef,
			 Complex_wp& coef_cte, Complex_wp& coef_grad, TinyVector<Complex_wp, 3>& v0,
			 int iquad, int i1, int j)
  {
    int num_pml = var_problem.mesh.Element(iquad).GetNumberPML();
    if (var_problem.mesh.GetPmlArea(num_pml).GetRadiusPML() > 0)
      {
	if (this->calcul_enveloppe || form_r3)
	  {
	    cout << "not implemented" << endl;
	    abort();
	  }
	
	// circular PML
	Complex_wp d = 1.0/var_boundary.GetTauPML(i1, j, 0);
	Complex_wp dbar = var_boundary.GetPrimitiveTauPML(i1, j, 1);
        
	Complex_wp exp_iteta = var_boundary.GetTauPML(i1, j, 1);
	Real_wp cos_teta = real(exp_iteta), sin_teta = imag(exp_iteta);
        
	Complex_wp dbar_d = dbar/d, d_dbar = d/dbar;
	
	// mu is assumed to be isotropic here
	Complex_wp mu = muPML(0, 0);
	muPML(0, 0) = mu*(dbar_d*square(cos_teta) + d_dbar*square(sin_teta));
	muPML(1, 1) = mu*(dbar_d*square(sin_teta) + d_dbar*square(cos_teta));
	muPML(0, 1) = mu*cos_teta*sin_teta*(dbar_d - d_dbar);
        
	rtilde = var_boundary.GetPrimitiveTauPML(i1, j, 0)*cos_teta;
	if (abs(rtilde) > epsilon_machine*var_problem.mesh.GetXmax())
	  invRtilde = 1.0/rtilde;
	
	coef = d*dbar;
	coef_cte = d*dbar;
	
	if (this->include_flow_term)
	  {
	    v0(1) *= d*dbar;
	    Complex_wp mr = v0(0), mz = v0(2);
	    v0(0) = mr*(square(cos_teta)*dbar + square(sin_teta)*d) + mz*cos_teta*sin_teta*(dbar-d);
	    v0(2) = mz*(square(sin_teta)*dbar + square(cos_teta)*d) + mr*cos_teta*sin_teta*(dbar-d);
	  }
      }
    else
      {
	Complex_wp dr = this->Glob_DrPML(iquad)(j);
	Complex_wp dz = this->Glob_DzPML(iquad)(j);
	// mu is considered orthotropic here
	muPML(0, 0) *= dz/dr;
	muPML(1, 1) *= dr/dz;
	coef = dr*dz;
	coef_grad = dr/dz;
	coef_cte = coef;
	if (this->calcul_enveloppe)
	  coef_cte -= dr/dz;
	
	if (this->include_flow_term)
	  {
	    v0(1) *= dr*dz;
	    v0(0) *= dz;
	    v0(2) *= dr;
	  }
      }
  }
  
  
  //! computing elementary matrix
  template<class Complexe> template<class T>
  void VarHelmholtz_Axi<Complexe>
  ::ComputeElementaryMatrixHelmAxi(int iquad, IVect& num_dof, VirtualMatrix<T>& mat,
				   const GlobalGenericMatrix<T>& nat_mat,
				   const ElementReference<Dimension2, 1>& Fb)
  {
    // initializing matrix
    int nb_dof_elt = Fb.GetNbDof();
    mat.Reallocate(nb_dof_elt, nb_dof_elt);
    mat.Zero();
    
    // dof numbers
    num_dof = var_problem.GetDofNumberOnElement(iquad);
    
    // we compute CoefMass = (-omega^2 r + n^2 / r) Ji(xi_k) omega_k
    // this coefficient is involved in mass matrix \int CoefMass u phi dx
    
    // and CoefStiff = r Ji(xi_k) omega_k DFi(xi_k)^{-1} DFi(xi_k)^{*-1}
    // this coefficient is involved in stiffness matrix \int CoefStiff \nabla u \nabla phi dx
    int N = Fb.GetNbPointsQuadratureInside();
    
    Vector<T> CoefMass(N);
    Vector<TinyVector<T, 2> > CoefGrad(N), CoefGradMinus(N);
    Vector<TinyMatrix<T, General, 2, 2> > CoefStiff(N);

    bool form_r3 = this->formulation_r3;
    if (var_boundary.GetCurrentModeNumber() == 0)
      form_r3 = false;
    
    Complexe m_iomega, m_omega2; 
    var_problem.GetMiomega(m_iomega); var_problem.GetMomega2(m_omega2);
    
    Complexe mu_teta;
    TinyMatrix<Complexe, Symmetric, 2, 2> muPML, mu_rz;
    TinyMatrix<Complexe, Symmetric, 3, 3> mu;
    TinyMatrix<Complexe, General, 2, 2> B, A;
    TinyVector<Complexe, 3> v0;
    
    bool affine = var_problem.mesh.IsElementAffine(iquad);
    int ref = var_problem.mesh.Element(iquad).GetReference();
    
    Real_wp jacob; Complexe coef, coef_grad, coef_cte;
    SetComplexOne(coef); SetComplexOne(coef_grad); SetComplexOne(coef_cte);
    
    TinyMatrix<Real_wp, General, 2, 2> dfjm1;
    
    Real_wp n2 = square(var_boundary.GetCurrentModeNumber());
    Real_wp m_ = Real_wp(var_boundary.GetCurrentModeNumber());
    Complexe ikwave_envelope, radius(0), r2;
    int i1 = iquad - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();

    SetComplexOne(r2);
    to_complex(Iwp*this->kwave_envelope, ikwave_envelope);
    for (int j = 0; j < N; j++)
      {
	Complexe rho = this->ref_rho(ref).GetCoefficient(var_problem, iquad, j);
	mu = this->ref_mu(ref).GetCoefficient(var_problem, iquad, j);	    
        Complexe sigma = this->ref_sigma(ref).GetCoefficient(var_problem, iquad, j);
	if (this->include_flow_term)
          v0 = this->ref_v0(ref).GetCoefficient(var_problem, iquad, j);        
        
	// mu_rz and mu_teta are extracted
	mu_teta = mu(2, 2);
	mu_rz(0, 0) = mu(0, 0);
	mu_rz(0, 1) = mu(0, 1);
	mu_rz(1, 1) = mu(1, 1);
	
        Complexe rtilde = this->Glob_rtilde(iquad)(j); 
        Complexe invRtilde; SetComplexZero(invRtilde);
        if (abs(rtilde) > epsilon_machine*var_problem.mesh.GetXmax())
          invRtilde = 1.0/rtilde;
	
	if (this->include_double_gradM_flow)
	  {
	    // to mu_rz, we subtract [m_r^2, m_r m_z; m_r m_z, m_z^2]
	    // to mu_teta, we sutract -m_theta^2 
	    // to v, we add m_theta m / omega [m_r;m_z]
            Complexe beta = this->ref_beta(ref).GetCoefficient(var_problem, iquad, j);
	    Complexe m_r = v0(0), m_z = v0(2);
	    mu_rz(0, 0) -= beta*m_r*m_r;
	    mu_rz(0, 1) -= beta*m_r*m_z;
	    mu_rz(1, 1) -= beta*m_z*m_z;
	    mu_teta -= beta*v0(1)*v0(1);
	    v0(0) += beta*v0(1)*m_*m_r*invRtilde / var_problem.GetOmega();
	    v0(2) += beta*v0(1)*m_*m_z*invRtilde / var_problem.GetOmega();
	  }
	        
        if (form_r3)
          {
            if (var_problem.InsidePML(iquad))
              radius = this->Glob_radius(iquad)(j);
            else
              radius = rtilde;

            r2 = radius*radius;
          }
        
        if (affine)
          {
            jacob = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(j);
            dfjm1 = var_problem.Glob_DFjm1(iquad)(0);                
            Mlt(1.0/var_problem.Glob_jacobian(iquad)(0), dfjm1);
          }
        else
          {
            jacob = var_problem.Glob_jacobian(iquad)(j);
            dfjm1 = var_problem.Glob_DFjm1(iquad)(j);
            Mlt(Fb.WeightsND(j)/jacob, dfjm1);
          }
        
        if (var_problem.InsidePML(iquad))
          {
	    ModifyCoefficientPML(form_r3, mu_rz, rtilde, invRtilde, coef, coef_cte,
				 coef_grad, v0, iquad, i1, j);
	    
	    Mlt(dfjm1, mu_rz, B); 
            MltTrans(B, dfjm1, A);
          }
        else
          {
            coef_cte = 1.0;
            if (this->calcul_enveloppe)
              coef_cte = 0.0;
            
            Mlt(dfjm1, mu_rz, B);
	    MltTrans(B, dfjm1, A);
          }
	
        CoefStiff(j) = A;
        CoefStiff(j) *= jacob*rtilde*nat_mat.GetCoefStiffness();
        
        if (this->calcul_enveloppe)
          {
            T coef_grad_ = ikwave_envelope*coef_grad*mu(0,0)*rtilde*jacob*nat_mat.GetCoefStiffness();

            CoefGrad(j)(0) = dfjm1(0, 1);
            CoefGrad(j)(1) = dfjm1(1, 1);
            CoefGrad(j) *= -coef_grad_;

            CoefGradMinus(j)(0) = dfjm1(0, 1);
            CoefGradMinus(j)(1) = dfjm1(1, 1);
            CoefGradMinus(j) *= coef_grad_;
          }
        else if (form_r3)
          {
            T coef_grad_ = mu(0, 0)*rtilde*radius*jacob*nat_mat.GetCoefStiffness()/coef_grad;

            CoefGrad(j)(0) = dfjm1(0, 0);
            CoefGrad(j)(1) = dfjm1(1, 0);
            CoefGrad(j) *= coef_grad_;

            CoefGradMinus(j)(0) = dfjm1(0, 0);
            CoefGradMinus(j)(1) = dfjm1(1, 0);
            CoefGradMinus(j) *= coef_grad_;
            
            CoefStiff(j) *= r2;
          }
        
        T rho_tilde = m_omega2*rho*nat_mat.GetCoefMass() + m_iomega*sigma*nat_mat.GetCoefDamping();
        
        if (form_r3)
          {
            CoefMass(j) = rtilde*r2*jacob*rho_tilde*coef_cte;
            CoefMass(j) += n2*r2*invRtilde*mu(0, 0)*jacob*coef*nat_mat.GetCoefStiffness();
            CoefMass(j) += mu(0, 0)*rtilde*jacob/coef_grad*nat_mat.GetCoefStiffness();
          }
        else
          {
            CoefMass(j) = rtilde*jacob*rho_tilde*coef_cte;
            CoefMass(j) += n2*invRtilde*mu_teta*jacob*coef*nat_mat.GetCoefStiffness();
          }
	
	if (this->include_flow_term)
	  {
	    TinyVector<T, 2> vec_u(v0(0), v0(2)), vec_v;
	    Mlt(dfjm1, vec_u, vec_v);
	    
	    T coef_v = rtilde*r2*m_iomega*jacob;
	    CoefGrad(j)(0) += coef_v*vec_v(0)*nat_mat.GetCoefDamping();
	    CoefGrad(j)(1) += coef_v*vec_v(1)*nat_mat.GetCoefDamping();
	    
	    CoefGradMinus(j)(0) -= coef_v*vec_v(0)*nat_mat.GetCoefDamping();
	    CoefGradMinus(j)(1) -= coef_v*vec_v(1)*nat_mat.GetCoefDamping();
	    
	    if (form_r3)
	      CoefMass(j) += 2.0*rtilde*radius*m_iomega*jacob*vec_v(0)*nat_mat.GetCoefDamping();
	    
	    coef_v = -2.0*var_problem.GetOmega()*m_*jacob*r2*v0(1);
	    CoefMass(j) += coef_v*nat_mat.GetCoefMass();
	  }
      }
    
    TinyVector<bool, 4> null_term(false, false, true, true);
    
    // then computing the elementary matrix
    if ((this->calcul_enveloppe) || form_r3 || (this->include_flow_term))
      {
	null_term(2) = false;
	null_term(3) = false;
      }

    // on ajoute la matrice \int CoefMass phi_i phi_j + \int CoefStiff nabla phi_i nabla phi_j
    // + gradient matrices
    Fb.AddVariableElemMatrix(0, 0, CoefMass, CoefStiff, CoefGradMinus, CoefGrad, null_term, mat);

    // modification of signs if needed
    //const Mesh<Dimension2>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat, iquad);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat, iquad);    
  }


  //! computes the impedance coefficient
  template<class Complexe>
  Complexe VarHelmholtz_Axi<Complexe>::ComputeImpedanceCoefficient(const IVect& ref, Vector<Complexe>& U0)
  {
    IVect ref_cond(var_problem.mesh.GetNbReferences()+1);
    ref_cond.Fill(0);
    for (int i = 0; i < ref.GetM(); i++)
      ref_cond(ref(i)) = 1;
    
    int order = var_problem.GetMeshNumbering(0).GetOrder();
    Mesh<Dimension2> mesh_subdiv;
    MeshInterpolationFEM<Dimension2> interp_gamma;
    interp_gamma.SetGaussQuadrature(order);    
    interp_gamma.InitProjectionSurface(var_problem.mesh);
    interp_gamma.ComputeSurfaceMesh(ref_cond, var_problem.mesh, mesh_subdiv,
				    var_problem);
    
#ifdef SELDON_WITH_MPI
    interp_gamma.GatherQuadraturePoints(var_problem.comm_group_mode, false);
#else
    interp_gamma.GatherQuadraturePoints(false);
#endif

    Vector<Complexe> trace_En;
    Vector<Complexe> trace_Hn;
    var_problem.ComputeEnHnOnBoundary(interp_gamma, U0, trace_En, trace_Hn);
    
    int offset = 0;
    Complexe int_P, int_dP_dn; Real_wp int_one(0);
    SetComplexZero(int_P); SetComplexZero(int_dP_dn);
    for (int i = 0; i < interp_gamma.GetNbBoundary(); i++)
      {
        int Nquad = interp_gamma.GetNbPointsQuadrature(i);
        for (int k = 0; k < Nquad; k++)
          {
            Real_wp poids = interp_gamma.WeightsQuadrature(i, k);
            poids *= interp_gamma.PointsQuadrature(i, k)(0);
            int_P += trace_En(offset + k)*poids;
            int_one += poids;
            int_dP_dn += trace_Hn(offset + k)*poids;
          }

        offset += Nquad;
      }

#ifdef SELDON_WITH_MPI
    if (var_problem.GetNbProcPerMode() > 1)
      {
        Complexe int_P0(int_P), int_dP0(int_dP_dn); Vector<int64_t> xtmp;
        Real_wp int_one0(int_one);
        MpiReduce(var_problem.comm_group_mode, &int_P0, xtmp, &int_P, 1, MPI_SUM, 0);
        MpiReduce(var_problem.comm_group_mode, &int_one0, xtmp, &int_one, 1, MPI_SUM, 0);
        MpiReduce(var_problem.comm_group_mode, &int_dP0, xtmp, &int_dP_dn, 1, MPI_SUM, 0);
        
        if (var_problem.GetRankProcMode() == 0)
          {
            DISP(int_P); DISP(int_one); DISP(int_dP_dn);    
          }
        
      }
#endif
   
    Complexe imped = int_P / (2.0*pi_wp*int_one*int_dP_dn);
    //exit(0);
    
    return imped;
  }
  
  
  /**********************
   * VarHelmholtzAxi_Eq *
   **********************/
  
  
  //! constructor with a given problem
  template<class TypeEquation>
  VarHelmholtzAxi_Eq<TypeEquation>::VarHelmholtzAxi_Eq()
    : VarHelmholtz_Axi<typename TypeEquation::Complexe>(static_cast<EllipticProblem<TypeEquation>& >(*this))
  {
    // list of mode numbers involved in the computation
    // by default, we take only mode 0 (source is axisymmetric)
    this->list_number_mode.Reallocate(1);
    this->list_number_mode(0) = 0;
    
    this->threshold_rhs = epsilon_machine;
    this->do_not_store_modes = true;
 }   

  
  //! retrieve and treat referenced edges
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::PerformOtherInitializations()
  {        
    this->var_gibc.InitGIBC();
    
    this->ComputeDofOnAxe(*this);
  }

  
  template<class TypeEquation>
  int VarHelmholtzAxi_Eq<TypeEquation>
  ::GetBoundaryConditionId(const IVect& ref, int pos, const VectString& parameters, bool& periodic)
  {
    int cond = ImpedanceABC_HelmholtzBase::
      GetBoundaryConditionId(ref, pos, parameters, periodic,
                             this->order_ABC, this->gamma_cla_coef, 
                             this->take_into_account_curvature_for_abc);
    
    if (cond >= 1)
      return cond;
    
    return VarHarmonic<TypeEquation>::GetBoundaryConditionId(ref, pos, parameters, periodic);
  }
  

  //! adds a Dirac source
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::AddDiracSource(const Real_wp& alpha_, Vector<Vector<Real_wp> > & b_source,
		   Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const
  {
    bool dirac = false;
    for (int k = 0; k < f.GetM(); k++)
      if (f(k) != NULL)
	if (f(k)->IsDiracSource())
	  dirac = true;

    if (dirac)
      {
	cout << "not possible" << endl;
	abort();
      }
  }

  
  //! adds a Dirac source
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::AddDiracSource(const Complex_wp& alpha_, Vector<Vector<Complex_wp> > & b_source,
		   Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const
  {
    VectR3 point_source3D(1);
    point_source3D(0) = this->origine_phase3D;
    Vector<bool> vec_unknown(1);
    vec_unknown.Fill(false);
    this->AddDiracSourceAxisym(alpha_, vec_unknown, this->number_mode, b_source,
			       point_source3D, f);
  }


  //! adds a Dirac source
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeDiracSource(Vector<Complexe>& b_source, const VectR3& point_source3D,
                       const Vector<VectString>& param_src) const
  {
    b_source.Reallocate(this->GetNbDof()*point_source3D.GetM());
    b_source.Zero();
    
    Vector<Vector<Complexe> > b_vec;
    b_vec.SetData(1, &b_source);
        
    Vector<bool> vec_unknown(1);
    vec_unknown.Fill(false);
    
    Complexe one; SetComplexOne(one);
    
    Vector<VirtualSourceFEM<Complexe, Dimension2>* > f(point_source3D.GetM());
    Vector<VectString> params(1);
    for (int k = 0; k < f.GetM(); k++)
      {
        params(0) = param_src(k);
        f(k) = this->GetNewDiracSourceEquationObject(params);
      }
    
    this->AddDiracSourceAxisym(one, vec_unknown, this->number_mode, b_vec,
			       point_source3D, f);
    
    for (int k = 0; k < f.GetM(); k++)
      delete f(k);
    
    b_vec.Nullify();
  }
  
  
  //! computation of source vector
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeRightHandSide(Vector<Real_wp>& b_rhs, bool assemble)
  {
    VarHarmonic<TypeEquation>::ComputeRightHandSide(b_rhs, assemble);
    
    if ((this->number_mode != 0) && (!this->formulation_r3))
      for (int i = 0; i < this->Dof_On_Axe.GetM(); i++)
        b_rhs(this->Dof_On_Axe(i)) = 0;    
  }


  //! computation of source vector
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeRightHandSide(Vector<Complex_wp>& b_rhs, bool assemble)
  {
    VarHarmonic<TypeEquation>::ComputeRightHandSide(b_rhs, assemble);
    
    if ((this->number_mode != 0) && (!this->formulation_r3))
      for (int i = 0; i < this->Dof_On_Axe.GetM(); i++)
        b_rhs(this->Dof_On_Axe(i)) = 0;
  }
  
  
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeRightHandSide(Vector<Vector<Complexe> >& b_rhs, bool assemble)
  {
    VarHarmonic<TypeEquation>::ComputeRightHandSide(b_rhs, assemble);
  }

  
  //! multiplication by phase exp(- i m theta)
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<Complex_wp>& val_u,
			Vector<Complex_wp>& grad_u,
			int i, const GridInterpolation<Dimension2>& var_interp,
			int iquad, bool compute_grad) const
  {
    Real_wp teta = 0;
    if (var_interp.GetNbTheta() > 0)
      teta = var_interp.GetTheta(i);
    
    Complex_wp coef = exp(-Iwp*Complex_wp(this->number_mode)*teta);
    
    val_u *= coef;
    if (compute_grad)
      grad_u *= coef;
    
    if ((this->formulation_r3) && (this->number_mode != 0))
      val_u *= var_interp.GetGlobalCoordinate(i)(0);	
  } 
  
  
  //! multiplication by phase cos(m theta)
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<Real_wp>& val_u,
			Vector<Real_wp>& grad_u,
			int i, const GridInterpolation<Dimension2>& var_interp,
			int iquad, bool compute_grad) const
  {
    Real_wp teta = 0;
    if (var_interp.GetNbTheta() > 0)
      teta = var_interp.GetTheta(i);
    
    Real_wp coef = cos(Real_wp(this->number_mode)*teta);
    val_u *= coef;
    if (compute_grad)
      grad_u *= coef;
    
    if ((this->formulation_r3) && (this->number_mode != 0))
      val_u *= var_interp.GetGlobalCoordinate(i)(0);
  } 

    
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
			int, bool, bool) const
  {
  }
    
  
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
			int, bool, bool) const 
  {
  }
  
  
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& grad_quadrature,
			  int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const
  {
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(nb_points_quad);
    
    for (int j = 0; j < nb_points_quad; j++)
      {
	En_quad(j) = u_quadrature(0)(j);
	if (compute_H)
          {
            R2 grad_u(grad_quadrature(0)(j), grad_quadrature(1)(j));
            Hn_quad(j) = DotProd(grad_u, normale(j));
          }
      }
  }


  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& grad_quadrature,
			  int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const
  {
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(nb_points_quad);
    
    for (int j = 0; j < nb_points_quad; j++)
      {
	En_quad(j) = u_quadrature(0)(j);
	if (compute_H)
	  {
            R2_Complex_wp grad_u(grad_quadrature(0)(j), grad_quadrature(1)(j));
            Hn_quad(j) = DotProd(grad_u, normale(j));
          }
      }
  }

  
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeEnHnNodal(Vector<VectReal_wp>& u_nodal, Vector<VectReal_wp>& grad_nodal,
		     int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
		     Vector<VectReal_wp>& En_nodal, Vector<VectReal_wp>& Hn_nodal) const
  {
    int nb_points_nodal = u_nodal(0).GetM();
    En_nodal.Reallocate(1); Hn_nodal.Reallocate(1);
    En_nodal(0).Reallocate(nb_points_nodal);
    Hn_nodal(0).Reallocate(nb_points_nodal);
    
    TinyVector<Real_wp, 2> val_v;
    for (int j = 0; j < nb_points_nodal; j++)
      {
	En_nodal(0)(j) = u_nodal(0)(j);
	ExtractVector(grad_nodal, j, 0, val_v);
	Hn_nodal(0)(j) = DotProd(val_v, normale(j));
      }
  }

    
  template<class TypeEquation>
  void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeEnHnNodal(Vector<VectComplex_wp>& u_nodal, Vector<VectComplex_wp>& grad_nodal,
		     int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
		     Vector<VectComplex_wp>& En_nodal, Vector<VectComplex_wp>& Hn_nodal) const
  {
    int nb_points_nodal = u_nodal(0).GetM();
    En_nodal.Reallocate(1); Hn_nodal.Reallocate(1);
    En_nodal(0).Reallocate(nb_points_nodal);
    Hn_nodal(0).Reallocate(nb_points_nodal);
    
    TinyVector<Complex_wp, 2> val_v;
    for (int j = 0; j < nb_points_nodal; j++)
      {
	En_nodal(0)(j) = u_nodal(0)(j);
	ExtractVector(grad_nodal, j, 0, val_v);
	Hn_nodal(0)(j) = DotProd(val_v, normale(j));
      }
  }

  
  //! returns true if the finite element matrix is symmetric
  template<class TypeEquation>
  bool VarHelmholtzAxi_Eq<TypeEquation>::IsSymmetricProblem(bool eigen) const
  {
    if (!VarHelmholtz_Axi<Complexe>::IsSymmetricProblem(eigen))
      return false;
    
    return VarHarmonic<TypeEquation>::IsSymmetricProblem(eigen);
  }

  
  template<class TypeEquation>
  void VarHelmholtzAxi_Cplx<TypeEquation>
  ::ComputeRightHandSide(Vector<Complex_wp>& b_rhs, bool assemble)
  {
    VarHelmholtzAxi_Eq<TypeEquation>::ComputeRightHandSide(b_rhs, assemble);
    
    if (this->GetOrderAbsorbingCondition() == ImpedanceABC_HelmholtzBase::VISCOTHERMAL)
      {
        ImpedanceABC_Helm<Dimension2> imped(this->GetLeafClass());    
        GlobalGenericMatrix<Complex_wp> nat_mat;
        int i = -1, iquad = -1, num_edge = -1, num_loc = -1, ref_domain = -1, pos_corner;
        for (int i0 = 0; i0 < this->mesh.GetNbElt(); i0++)
          {
            bool presence_corner = false;
            for (int j = 0; j < this->mesh.Element(i0).GetNbVertices(); j++)
              {                
                int jp1 = (j+1)%this->mesh.Element(i0).GetNbVertices();
                int ne1 = this->mesh.Element(i0).numBoundary(j);
                int ne2 = this->mesh.Element(i0).numBoundary(jp1);
                int ref1 = this->mesh.Boundary(ne1).GetReference();
                int ref2 = this->mesh.Boundary(ne2).GetReference();
                if ((ref1 == 1) && (this->mesh.GetBoundaryCondition(ref2) == BoundaryConditionEnum::LINE_ABSORBING))
                  {
                    num_edge = ne2; num_loc = jp1; i = ne2; pos_corner = jp1;
                    presence_corner = true;
                  }

                if ((ref2 == 1) && (this->mesh.GetBoundaryCondition(ref1) == BoundaryConditionEnum::LINE_ABSORBING))
                  {
                    num_edge = ne1; num_loc = j; i = ne1; pos_corner = jp1;
                    presence_corner = true;
                  }
              }

            if (presence_corner)
              {
                iquad = i0;
                ref_domain = this->mesh.Element(i0).GetReference();
                break;
              }
          }

        if (num_edge >= 0)
          {
            SetPoints<Dimension2> pts; SetMatrices<Dimension2> mat;
            VectR2 s;
            const ElementReference_Dim<Dimension2>& Fb = this->GetReferenceElement(iquad);
            this->mesh.GetVerticesElement(iquad, s);
            Fb.FjElemQuadrature(s, pts, this->mesh, iquad);
            Fb.DFjElemQuadrature(s, pts, mat, this->mesh, iquad);
            Fb.FjSurfaceElem(s, pts, this->mesh, iquad, num_loc);    
            Fb.DFjSurfaceElem(s, pts, mat, this->mesh, iquad, num_loc);
            
            imped.EvaluateImpedancePhi(i, iquad, num_edge, num_loc, 0, nat_mat, ref_domain,
                                       pts, mat);
            
            int num_dof = this->mesh_num.Element(iquad).GetNumberDof(pos_corner);
            Complex_wp coef = imped.GetStoredGradCoefficient(0, 0);
            Real_wp radius = s(pos_corner)(0);
            b_rhs(num_dof) += coef*radius;
          }
      }
  }

  
  /****************************************
   * Impedance for axisymmetric Helmholtz *
   ****************************************/
  
  
  //! evaluation of impedance coefficient
  void ImpedanceABC<Complex_wp, HelmholtzEquationAxi>::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                       const SetPoints<Dimension2>& Pts,
                       const SetMatrices<Dimension2>& Mat)
  {
    ImpedanceABC_Helm<Dimension2>::
      EvaluateImpedancePhi(i, num_elem, num_edge, num_loc, k, nat_mat, ref, Pts, Mat);

    ImpedanceABC_Helm<Dimension2>::
      EvaluateImpedanceGrad(i, num_elem, num_edge, num_loc, k, nat_mat, ref, Pts, Mat);

    if (this->stored_coef_phi(0).GetM() <= 0)
      {
        abort();
      }
    
    //DISP(i); DISP(this->stored_coef_phi(0)(k));
    //DISP(this->stored_coef_grad(0)(k));
    Real_wp m2 = square(var_boundary.GetCurrentModeNumber());
    Real_wp radius = Pts.GetPointQuadratureBoundary(k)(0);
    this->stored_coef_phi(0)(k) *= radius;
    // cas particulier pour 36 car le terme dans le div grad surfacique est une matrice : anisotrope
    if( this->var_boundary.GetOrderAbsorbingCondition() == 36)
      {	 
	int num_point = this->var_problem.GetNbPointsQuadratureInside(num_elem);
	Dimension2::R_N normale = Mat.GetNormaleQuadratureBoundary(k);
	Complex_wp mu_n = this->var_helm.GetMuNormale(ref, num_elem, num_point+k, normale);
	Complex_wp rho = this->var_helm.GetRhoTilde(ref, num_elem, num_point+k);
	Complex_wp c = sqrt(mu_n/rho);
	Complex_wp ik = Iwp*this->var_problem.GetOmega()/c; // pas de signe - pour les RBC
	Real_wp k1 = Mat.GetK1QuadratureBoundary(k);
	Real_wp k2 = Mat.GetK2QuadratureBoundary(k);
	
	// correcting coefs of the diagonal of \mathcal{R}
	Complex_wp coef_cor_grad = Real_wp(1.0)+k1/(ik);
	Complex_wp coef_cor_phi  = Real_wp(1.0)+k2/(ik);
	
	if (radius > 0)
	  this->stored_coef_phi(0)(k) += coef_cor_phi*m2*this->stored_coef_grad(0)(k)/radius; 
	
	this->stored_coef_grad(0)(k) *= coef_cor_grad;
      }
    else
      {
	if (radius > 0)
	  this->stored_coef_phi(0)(k) += m2*this->stored_coef_grad(0)(k)/radius;		
      }
    
    this->stored_coef_grad(0)(k) *= radius;
    //this->stored_coef_grad(0)(k) = 0.0;
    
    if ((var_helm.UseFormulationR3()) && (var_boundary.GetCurrentModeNumber() != 0))
      {
	this->stored_coef_phi(0)(k) *= radius*radius;
	
	if (this->stored_coef_grad(0)(k) != Complex_wp(0, 0))
	  {
	    cout << "Not implemented for R^3 formulation" << endl;
	    abort();
	  }
      }   
  }
  
  
  //! evaluation of impedance coefficient
  template<class T>
  void ImpedanceABC<T, LaplaceEquationAxi>::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                       const GlobalGenericMatrix<T>& nat_mat, int ref,
                       const SetPoints<Dimension2>& Pts,
                       const SetMatrices<Dimension2>& Mat)
  {
    Real_wp radius = Pts.GetPointQuadratureBoundary(k)(0);
    Real_wp mu = var_helm.ref_mu(ref).GetConstant()(0, 0);
    Real_wp rho = var_helm.ref_rho(ref).GetConstant();
    Real_wp imped = sqrt(rho*mu);

    if (k == 0)
      {
        this->stored_coef_phi.Reallocate(1);
	int nb_points = Mat.GetNbPointsQuadratureBoundary();
	this->stored_coef_phi(0).Reallocate(nb_points);
	this->stored_coef_phi(0).Fill(0);
      }
   
    this->stored_coef_phi(0)(k) = radius*imped*nat_mat.GetCoefDamping();
    this->coef_grad = 0;
  }

  
  //! computes impedance
  template<class Complexe>
  void ImpedanceGenericAxiHelm<Complexe>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                         const GlobalGenericMatrix<Complexe>& nat_mat, int ref_d,
                         const SetPoints<Dimension2>& Pts,
                         const SetMatrices<Dimension2>& Mat)
  {
    ImpedanceFunction_Base<Complexe, Dimension2>::
      EvaluateImpedancePhi(i, num_elem, num_edge, num_loc, k, nat_mat, ref_d, Pts, Mat);

    if (k == 0)
      {
        this->stored_coef_phi.Reallocate(1);
        int nb_pts = Mat.GetNbPointsQuadratureBoundary();
        this->stored_coef_phi(0).Reallocate(nb_pts);
      }

    Real_wp radius = Pts.GetPointQuadratureBoundary(k)(0);
    this->stored_coef_phi(0)(k) = this->coef_phi(0)*radius;
    if (var_helm.UseFormulationR3() && (var_boundary.GetCurrentModeNumber() != 0))
      this->stored_coef_phi(0)(k) *= radius*radius;
  }

  
  /*****************************
   * PlaneWaveAxiIncidentField *
   *****************************/


  //! constructor
  PlaneWaveAxiIncidentField
  ::PlaneWaveAxiIncidentField(const R3& kwave3D, const Complex_wp& rho0, const Complex_wp& sigma0,
			      const Complex_wp& mu0, const Real_wp& omega, int m)
  {
    Complex_wp rho_tilde = rho0 + Iwp*sigma0/omega;
    VarAxisymProblem::Get_KwavePerp_Kz_Phase(rho_tilde, mu0, m, kwave3D, omega,
					     k_perp, kz, incidence_axial, phase);
    
    number_mode = m;    
  }

  
  //! computes component of incident field in Fourier expansion
  void PlaneWaveAxiIncidentField
  ::EvaluateFunction(const R2& x, Complex_wp& uinc) const
  {
    uinc = Complex_wp(cos(kz*x(1)), sin(kz*x(1)));    
    
    if (!incidence_axial)
      {
        Real_wp kr = k_perp*x(0);
        uinc *= phase*GetJn(-number_mode, kr);
      }
    else
      {
        if (number_mode != 0)
          uinc = 0;
      }
  }
  
  
  //! computes component of gradient of incident field in Fourier expansion
  void PlaneWaveAxiIncidentField
  ::EvaluateFunctionGradient(const R2& x, Complex_wp& u_inc,
			     TinyVector<Complex_wp, 2>& grad_u) const
  {
    u_inc = Complex_wp(cos(kz*x(1)), sin(kz*x(1)));
    Complex_wp du_dr = u_inc;
    Complex_wp du_dz = Iwp*kz*du_dr;
    
    if (incidence_axial)
      {
        du_dr = 0;
        if (number_mode != 0)
          {
	    du_dz = 0;
	    u_inc = 0;
	  }
      }
    else
      {
        Real_wp kr = k_perp*x(0);        
        Real_wp jn_, djn_;
        GetDeriveJn(-number_mode, kr, jn_, djn_);
        u_inc *= phase*jn_;
        du_dr *= phase*k_perp*djn_;
        du_dz *= phase*jn_;
      }
    
    grad_u.Init(du_dr, du_dz);
  }

  
  /********************************
   * DiffractedWaveSource_HelmAxi *
   ********************************/
  
  
  //! Sets attributes to default values
  void IncidentWaveProjector_HelmAxi::InitDefaultValues()
  {
    Complex_wp rho0 = var_helm.rho0 + Iwp*var_helm.sigma0/omega;
    Complex_wp mu0 = var_helm.mu0;
    
    VarAxisymProblem::Get_KwavePerp_Kz_Phase(rho0, mu0, number_mode,
					     var_helm.GetWaveVector(), omega,
					     k_perp, kz, incidence_axial, phase);
  }


  //! Evaluates incident field
  void IncidentWaveProjector_HelmAxi::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    Complex_wp u_inc;
    this->incident_wave.EvaluateFunction(x, u_inc);
    f(0) = u_inc;
    
    if (number_mode != 0)
      if (var_helm.UseFormulationR3() && (!incidence_axial))
	{
	  if (x(0) == 0)
	    {
	      f(0) = Complex_wp(cos(kz*x(1)), sin(kz*x(1)));
	      if (number_mode == 1)
		f(0) *= -0.5*k_perp*phase;
	      else if (number_mode == -1)
		f(0) *= 0.5*k_perp*phase;
	      else
		f(0) = Complex_wp(0, 0);              
	    }
	  else
	    f(0) /= x(0);
	}
  }
  
  
  //! Sets attributes to default values  
  void DiffractedWaveSource_HelmAxi::InitDefaultValues()
  {
    coef_vol = 0.0;
    coef_grad = 0.0;
    coef_m2 = 0.0;

    rho0 = var_helm.rho0 + Iwp*var_helm.sigma0/var_problem.GetOmega();
    mu0 = var_helm.mu0;
    number_mode = var_boundary.GetCurrentModeNumber();
    
    VarAxisymProblem::Get_KwavePerp_Kz_Phase(rho0, mu0, int(number_mode), var_helm.GetWaveVector(),
					     var_problem.GetOmega(),
					     k_perp, kz, incidence_axial, phase);
    
    formulation_r3 = var_helm.UseFormulationR3();
    if (var_boundary.GetCurrentModeNumber() == 0)
      formulation_r3 = false;
  }
  
  
  //! initialization for each element
  /*!
    \param[in] num_elem element number
    if physical properties are constant inside the element
    coefficients \f$ -\omega^2(\rho - \rho_0) \f$ and (mu-mu0) are computed
   */
  void DiffractedWaveSource_HelmAxi
  ::InitElement(int num_elem, const VectR2& s)
  {
    int ref = var_problem.mesh.Element(num_elem).GetReference();
    Complex_wp rho = var_helm.ref_rho(ref).GetConstant()
      + Iwp*var_helm.ref_sigma(ref).GetConstant()/var_problem.GetOmega();
    
    Complex_wp mu = var_helm.ref_mu(ref).GetConstant()(0, 0);
    
    if (rho != rho0)
      coef_vol = var_problem.GetSquareOmega()*(rho - rho0);
    else
      coef_vol = 0;
    
    if (mu != mu0)
      {
        coef_grad = mu0 - mu;
        coef_m2 = coef_grad*square(number_mode);
      }
    else
      {
        coef_grad = 0;
        coef_m2 = 0;
      }
    
    this->num_elem_ = num_elem;
    this->num_loc_ = -1;
    
    // for pml element, no source
    if (var_problem.InsidePML(num_elem))
      {
        coef_vol = 0;
        coef_grad = 0;
        coef_m2 = 0;
      }
  }
  
  
  //! volumetric source if rho <> rho0
  bool DiffractedWaveSource_HelmAxi
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    // if rho different from rho 0, we have a volumic source \int f \phi
    if ((coef_vol != Complex_wp(0, 0)) || (coef_m2 != Complex_wp(0, 0)))
      return true;
    
    return false;
  }
  
  
  //! evaluation of incident wave
  void DiffractedWaveSource_HelmAxi
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    // incident wave 
    incident_wave.EvaluateFunction(x, f(0));
    
    Complex_wp zero;
    SetComplexZero(zero);
    
    if (formulation_r3 && (!incidence_axial))
      {
        if (x(0) == 0)
          {
            f(0) = Complex_wp(cos(kz*x(1)), sin(kz*x(1)));
            if (var_boundary.GetCurrentModeNumber() == 1)
              f(0) *= -0.5*k_perp*phase;
            else if (var_boundary.GetCurrentModeNumber() == -1)
              f(0) *= 0.5*k_perp*phase;
            else
              f(0) = zero;              	    
          }
        else
          f(0) /= x(0);        
      }

    // for Dirichlet condition u = -u_inc
    if (this->num_loc_ >= 0)
      {
        f(0) = -f(0);
        // for PML, no source
        if (var_problem.InsidePML(i))
          f.Fill(zero);
      }
  }
  
  //! evaluation of volumetric source
  void DiffractedWaveSource_HelmAxi
  ::EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    Complex_wp u_inc;
    // volumetric source -> omega2 \int (rho-rho_0) u^inc \phi
    if (formulation_r3)
      {
        Complex_wp du_dr, du_dz; TinyVector<Complex_wp, 2> grad_u_inc;
	incident_wave.EvaluateFunctionGradient(x, u_inc, grad_u_inc);
	du_dr = grad_u_inc(0); du_dz = grad_u_inc(1);
        
        f(0) = u_inc*(coef_vol*x(0)*x(0) + coef_m2);
        f(0) += coef_grad*du_dr*x(0);
      }
    else
      {
	incident_wave.EvaluateFunction(x, u_inc);
        f(0) = u_inc*x(0)*coef_vol;
        if (x(0) > 1e-12)
          f(0) += u_inc*coef_m2/x(0);    
      }
  }
  
  
  //! volumetric source for \f$ \nabla phi \f$ if mu <> mu0
  bool DiffractedWaveSource_HelmAxi
  ::IsNonNullGradientSource(const VectR2& s)
  {
    // gradient source if mu different from mu_0
    if (coef_grad != Complex_wp(0))
      return true;
    
    return false;
  }
  
  
  //! evaluation of volumetric source g in \f$ \int g \nabla \varphi \f$
  void DiffractedWaveSource_HelmAxi
  ::EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    // gradient source -> \int (mu_0 - mu) \nabla u^{inc} \cdot \nabla \phi exp(ikx)
    Complex_wp du_dr, du_dz, u_inc; TinyVector<Complex_wp, 2> grad_u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_u_inc);
    du_dr = grad_u_inc(0); du_dz = grad_u_inc(1);
    
    f(0) = coef_grad*x(0)*du_dr;
    f(1) = coef_grad*x(0)*du_dz;
    if (formulation_r3)
      f *= x(0);
  }
  
  
  //! initialization before evaluation of surfacic source
  void DiffractedWaveSource_HelmAxi
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension2>::InitSurface(i, num_face, num_elem, num_loc);
    
    // for neumann condition and high conductivity condition
    //int ref = this->vars.mesh.Element(num_elem).GetReference();
    coef_grad = mu0;
    if (var_problem.InsidePML(num_elem))
      {
        coef_grad = 0;
      }
  }
  
  
  //! surfacic for Neumann and High-Conductivity boundary conditions
  bool DiffractedWaveSource_HelmAxi::IsNonNullSurfacicSource(int ref)
  {    
    // for neumann condition and high conductivity condition
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref;
    if ( (cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
         || (cond_ref == BoundaryConditionEnum::LINE_IMPEDANCE))
      return true;
    
    return false;
  }
  
  
  //! evaluation of surfacic source
  void DiffractedWaveSource_HelmAxi
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    Complex_wp du_dr, du_dz, u_inc;
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    TinyVector<Complex_wp, 2> grad_u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_u_inc);
    du_dr = grad_u_inc(0); du_dz = grad_u_inc(1);

    R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp du_inc_dn = du_dr*normale(0) + du_dz*normale(1);
    f.Fill(0);
    if (reference_condition == BoundaryConditionEnum::LINE_NEUMANN)
      {
	// neumann condition \int -\mu_0 du_inc/dn \phi 	
	f(0) = -coef_grad*du_inc_dn*x(0);
        if (formulation_r3)
          f(0) *= x(0);
      }
    else if (reference_condition == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
        // retrieving value of Z
	GlobalGenericMatrix<Complex_wp> nat_mat;
        int ref = var_problem.mesh.Element(this->num_elem_).GetReference();

	Complex_wp imped
	  = fct_imped.GetCoefficient(this->num_bound_ref_, this->num_elem_, this->num_loc_, k,
				     ref, this->ref_boundary_, PointsElem, MatricesElem);
	
        // neumann condition -\int \mu_0 du_inc/dn \phi 	
        f(0) = -coef_grad*du_inc_dn*x(0);        
        f(0) -= imped*u_inc*x(0);
        if (formulation_r3)
          f(0) *= x(0);        
      }
  }
  
  
  //! surfacic source \f$ \int f \frac{d\varphi}{ds} \f$  for high-conductivity boundary condition
  bool DiffractedWaveSource_HelmAxi
  ::IsNonNullSurfacicSourceGradient(int ref)
  {
    return false;
  }
  
  
  //! evaluation of surfacic source f in \f$ \int f \frac{d \varphi}{ds} \f$ 
  void DiffractedWaveSource_HelmAxi
  ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension2>& PointsElem,
				   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    f.Fill(0);
  }
  
  
  /***************************
   * TotalWaveSource_HelmAxi *
   ***************************/
  
  
  //! Sets values of attributes to zero
  void TotalWaveSource_HelmAxi::InitDefaultValues()
  {
    rho0 = var_helm.rho0 + Iwp*var_helm.sigma0/var_problem.GetOmega();
    mu0 = var_helm.mu0;
    Real_wp invC0 = sqrt(abs(rho0/mu0));
    k0 = var_problem.GetOmega()*invC0;
    number_mode = var_boundary.GetCurrentModeNumber();
    
    VarAxisymProblem::Get_KwavePerp_Kz_Phase(rho0, mu0, int(number_mode), var_helm.GetWaveVector(),
					     var_problem.GetOmega(),
					     k_perp, kz, incidence_axial, phase);

    formulation_r3 = var_helm.UseFormulationR3();
    if (var_boundary.GetCurrentModeNumber() == 0)
      formulation_r3 = false;
  }
  

  //! Is the element in the PML ?
  void TotalWaveSource_HelmAxi::InitElement(int num_elem, const VectR2& s)
  {
  }
  
  
  //! for Dirichlet condition
  void TotalWaveSource_HelmAxi
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    f.Fill(0);
  }
  
  
  //! volumetric source if PML element
  bool TotalWaveSource_HelmAxi
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    return false;
  }
  
  //! evaluation of volumetric source
  void TotalWaveSource_HelmAxi
  ::EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    f.Fill(0);
  }
  
  //! volumetric source if PML element
  bool TotalWaveSource_HelmAxi
  ::IsNonNullGradientSource(const VectR2& s)
  {
    return false;
  }
  
  //! evaluation of volumetric source f in \f$ \int f \nabla \varphi \f$
  void TotalWaveSource_HelmAxi
  ::EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    f.Fill(0);
  }
  
  //! initialization before computation of surfacic source
  void TotalWaveSource_HelmAxi
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension2>::InitSurface(i, num_face, num_elem, num_loc);
  }
  
  
  //! non null surfacic source for absorbing condition
  bool TotalWaveSource_HelmAxi
  ::IsNonNullSurfacicSource(int ref)
  {
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    if (cond_ref == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    return false;
  }
  
  //! evaluation of surfacic source
  void TotalWaveSource_HelmAxi
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    Complex_wp du_dr, du_dz;
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc;
    TinyVector<Complex_wp, 2> grad_u_inc;
    incident_wave.EvaluateFunctionGradient(x, u_inc, grad_u_inc);
    du_dr = grad_u_inc(0); du_dz = grad_u_inc(1);
    
    R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Complex_wp du_inc_dn = du_dr*normale(0) + du_dz*normale(1);
    
    f(0) = mu0*(du_inc_dn - Iwp*k0*u_inc)*x(0);
    if (formulation_r3)
      f(0) *= x(0);
  }


  /****************************
   * VolumetricSource_HelmAxi *
   ****************************/

  
  //! default constructor with given problem
  template<class T> template<class TypeEquation>
  VolumetricSource_AxiHelm<T>::
  VolumetricSource_AxiHelm(const EllipticProblem<TypeEquation>& var,
			   const Vector<VectString>& param)	   
    : VirtualSourceFEM<T, Dimension2>(var), var_boundary(var), var_helm(var), var_problem(var)
  { 
    coef_vol = var.GetCoefficientVolumeSource();
    fsrc = NULL;

    GaussianSourceField<T, Dimension3>* f_gauss = NULL;
    
    Vector<T> polar;
    var.GetPolarization(polar);

    unif_source_vol.Reallocate(var_problem.GetNbPhysicalIndices()+1);
    coef_source_vol.Reallocate(var_problem.GetNbPhysicalIndices()+1);
    unif_source_surf.Reallocate(var_problem.mesh.GetNbReferences()+1);
    coef_source_surf.Reallocate(var_problem.mesh.GetNbReferences()+1);
    unif_source_vol.Fill(false); unif_source_surf.Fill(false);
    coef_source_surf.Zero(); coef_source_vol.Zero();
    
    surface_source = false;
    bool toroidal_source = false;
    external_source = false;
    ref_surf = -1;
    for (int num = 0; num < param.GetM(); num++)
      {
	if (param(num).GetM() > 1)
	  if ((param(num)(0) == "SRC_VOLUME") || (param(num)(0) == "SRC_SURFACE"))
	  {
	    for (int k = 1; k < param(num).GetM(); k++)
	      {
		if (param(num)(k) == "Polarization")
		  {
		    for (int p = 0; p < polar.GetM(); p++)
		      polar(p) = to_num<T>(param(num)(k+p+1));
		    
		    break;
		  }
		
		if (param(num)(k) == "Torus")
		  toroidal_source = true;
	      }

	    int nb = 1;
            if (param(num)(2) == "UNIFORM")
              {
                int ref = to_num<int>(param(num)(1));
                if (param(num)(0) == "SRC_VOLUME")
                  {
                    unif_source_vol(ref) = true;
                    coef_source_vol(ref) = polar(0);
                  }

                if (param(num)(0) == "SRC_SURFACE")
                  {
                    unif_source_surf(ref) = true;
                    coef_source_surf(ref) = polar(0);
                  }
              }
            else
              {
                if (param(num)(0) == "SRC_SURFACE")
                  {
                    surface_source = true;	      
                    ref_surf = to_num<int>(param(num)(nb++));
                  }
                
                if (param(num)(nb) == "GAUSSIAN")
                  {
                    nb++;
                    f_gauss = new GaussianSourceField<T, Dimension3>();                
                    var.InitGaussianParameter(*f_gauss, param(num), nb);
                    f_gauss->SetPolarization(polar);
                    fsrc = f_gauss;
                  }
              }
	  }
      }

    if (f_gauss == NULL)
      return;
    
    Real_wp r = f_gauss->GetRadius();
    rmax = f_gauss->GetCutOffRadius();
    Real_wp x_c = sqrt(square(f_gauss->GetCenter()(0)) + square(f_gauss->GetCenter()(1)));
    Real_wp z_c = f_gauss->GetCenter()(2);
    center_source.Init(x_c, z_c);
    
    if ((x_c != Real_wp(0)) && (!toroidal_source))
      {
	int N = toInteger(ceil(2.0*pi_wp*(x_c+r)/r));
	int m = var_boundary.GetCurrentModeNumber();
	
	gauss_teta.ConstructQuadrature(max(4*N, 4*m));
	/* 
	T f_prev, f(0);
	for (int p = 0; p < 400; p += 10)
	  {
	    gauss_teta.ConstructQuadrature(p);
	    Real_wp cos_teta, sin_teta; R3 pt3D;
	    Vector<T> feval(1); f_prev = f;
	    f = T(0);
	    for (int k = 0; k < gauss_teta.GetNbPointsQuad(); k++)
	      {
		Real_wp teta = 2*pi_wp*gauss_teta.Points(k);
		cos_teta = cos(teta); sin_teta = sin(teta);
		pt3D(0) = x*cos_teta;
		pt3D(1) = x*sin_teta;
		pt3D(2) = z;
		this->fsrc->EvaluateFunction(pt3D, feval);
		feval *= exp(Iwp*m*teta);
		f += gauss_teta.Weights(k)*feval(0);
	      }
	    
	    DISP(p); DISP(f); DISP(abs(f_prev - f)/abs(f));
	  }
	*/
      }
    else
      gauss_teta.ConstructQuadrature(0);
  }


  template<class T>
  VolumetricSource_AxiHelm<T>::~VolumetricSource_AxiHelm()
  {
    if (fsrc != NULL)
      delete fsrc;
  }
  

  //! sets the volume source
  template<class T>
  void VolumetricSource_AxiHelm<T>::SetVolumeSourceFunction(VirtualSourceField<T, Dimension3>* f)
  {
    if (fsrc != NULL)
      delete fsrc;
    
    fsrc = f;
  }

  template<class T>
  void VolumetricSource_AxiHelm<T>::SetSurfaceSource(int ref, VirtualSourceField<T, Dimension3>* f)
  {
    if (fsrc != NULL)
      delete fsrc;
    
    surface_source = true;
    external_source = true;
    ref_surf = ref;
    fsrc = f;
  }

  
  //! volumetric source
  template<class T>
  bool VolumetricSource_AxiHelm<T>::
  IsNonNullVolumetricSource(const VectR2& s)
  {
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (unif_source_vol(ref))
      return true;

    if (fsrc == NULL)
      return false;
    
    return !surface_source;
  }

  template<class T>
  void VolumetricSource_AxiHelm<T>::
  EvaluateSourceMode(const R2& x, Vector<T>& f)
  {
    
    R3 pt3D;
    if (gauss_teta.GetOrder() != 0)
      {
	Vector<T> feval(f);
	Real_wp teta, cos_teta, sin_teta;
	T poids, mode;
	for (int k = 0; k < gauss_teta.GetNbPointsQuad(); k++)
	  {
	    teta = 2*pi_wp*gauss_teta.Points(k);
	    cos_teta = cos(teta); sin_teta = sin(teta);
	    pt3D(0) = x(0)*cos_teta;
	    pt3D(1) = x(0)*sin_teta;
	    pt3D(2) = x(1);
	    this->fsrc->EvaluateFunction(pt3D, feval);
	    this->var_helm.GetFourierMode(teta, mode);
	    poids = gauss_teta.Weights(k)*mode;
	    f(0) += poids*feval(0);
	  }
      }
    else
      {
	Vector<T> feval(f.GetM());
	pt3D.Init(x(0), 0, x(1));
        if ( (var_boundary.GetCurrentModeNumber() == 0) || external_source)
          {
            fsrc->EvaluateFunction(pt3D, feval);
            f(0) = feval(0);
          }
      }
  }

  
  //! Evaluation of volumetric source f (term \int f \varphi)
  template<class T>
  void VolumetricSource_AxiHelm<T>
  ::EvaluateVolumetricSource(int i, int j, const R2& x, Vector<T>& f)
  {
    // integration in theta
    f(0) = 0.0;

    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    if (unif_source_vol(ref))
      {
        f(0) = coef_source_vol(ref)*x(0);
        if ((var_helm.UseFormulationR3()) && (var_boundary.GetCurrentModeNumber() != 0))
          f(0) *= x(0);

        return;
      }
    
    if (fsrc == NULL)
      return;
      
    if (center_source.Distance(x) >= rmax)
      return;
    
    EvaluateSourceMode(x, f);
    f(0) *= coef_vol*x(0);
    if ((var_helm.UseFormulationR3()) && (var_boundary.GetCurrentModeNumber() != 0))
      Mlt(x(0), f);
  }
  
  
  //! returns true if there is a source term \int_\Gamma f \varphi 
  template<class T>
  bool VolumetricSource_AxiHelm<T>::IsNonNullSurfacicSource(int ref)
  {
    if (unif_source_surf(ref))
      return true;

    if (fsrc == NULL)
      return false;
    
    if (!surface_source)
      return false;
    
    if (ref == ref_surf)
      return true;
    
    return false;
  }
  
  
  //! evaluates f in source term \int_\Gamma f \varphi 
  template<class T>
  void VolumetricSource_AxiHelm<T>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, Vector<T>& f)
  {
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    
    // integration in theta
    f(0) = 0.0;

    int ref = this->ref_boundary_;
    if (unif_source_surf(ref))
      {
        f(0) = coef_source_surf(ref)*x(0);
        if ((var_helm.UseFormulationR3()) && (var_boundary.GetCurrentModeNumber() != 0))
          f(0) *= x(0);

        return;
      }

    if (fsrc == NULL)
      return;

    if (center_source.Distance(x) >= rmax)
      return;
    
    EvaluateSourceMode(x, f);
    
    f(0) *= x(0);
    if ((var_helm.UseFormulationR3()) && (var_boundary.GetCurrentModeNumber() != 0))
      Mlt(x(0), f);    
  }


  template<class T>
  void VolumetricSource_AxiHelm<T>
  ::EvaluateFunction(int i, int j, const R2& x, Vector<T>& f)
  {
    f(0) = 0.0;
    if (this->num_loc_ >= 0)
      {
        int ref = this->ref_boundary_;
        if (unif_source_surf(ref))
          f(0) = coef_source_surf(ref);

        if (fsrc == NULL)
          return;
        
        EvaluateSourceMode(x, f);
      }
  }


  void FemMatrixFreeClass<Complex_wp, HelmholtzEquationAxi>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "not possible" << endl;
    abort();
  }


  void FemMatrixFreeClass<Complex_wp, HelmholtzEquationAxi>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddScalarH1(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		   X, Complex_wp(1, 0), Y, false);
  }
  
}

#define MONTJOIE_FILE_AXISYM_HELMHOLTZ_CXX
#endif


