#ifndef MONTJOIE_FILE_VAR_PROBLEM_1D_CXX

namespace Montjoie
{
  
  /**********************
   * PmlDampingFunction *
   **********************/

  
  //! Default constructor
  PmlDampingFunction::PmlDampingFunction()
  {
    vsigma = 2.0;
    exponent_pml = 2.0;
    type_function = PML_PARABOLE;
    max_velocity = 1.0;
    offset_damping = 0.0;
    xmin = 0.0; xmax = 0.0;
    thickness = 1.0;
  }
  
  
  //! modifies parameters of the object with a line DampingPML = ...
  void PmlDampingFunction::SetInputData(const VectString& parameters)
  {
    if (parameters.GetM() <= 0)
      {
	cout << "In SetInputData of PmlDampingFunction" << endl;
	cout << "DampingPML needs more parameters, for instance :" << endl;
	cout << "DampingPML = vsigma" << endl;
	cout << "Current parameters are : " << endl << parameters << endl;
	abort();
      }
    
    if (parameters(0) == "PARABOLE")
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of PmlDampingFunction" << endl;
	    cout << "DampingPML needs more parameters, for instance :" << endl;
	    cout << "DampingPML = PARABOLE vsigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_function = PML_PARABOLE;
	vsigma = to_num<Real_wp>(parameters(1));
	exponent_pml = 2.0;
      }
    else if (parameters(0) == "SHIFTED_PARABOLE")
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of PmlDampingFunction" << endl;
	    cout << "DampingPML needs more parameters, for instance :" << endl;
	    cout << "DampingPML = SHIFTED_PARABOLE vsigma offset" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_function = PML_SHIFTED_PARABOLE;
	vsigma = to_num<Real_wp>(parameters(1));
	offset_damping = to_num<Real_wp>(parameters(2));
	exponent_pml = 2.0;
      }
    else if (parameters(0) == "CONSTANT")
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of PmlDampingFunction" << endl;
	    cout << "DampingPML needs more parameters, for instance :" << endl;
	    cout << "DampingPML = CONSTANT vsigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_function = PML_CONSTANT;
	vsigma = to_num<Real_wp>(parameters(1));
	exponent_pml = 0.0;
      }
    else if (parameters(0) == "EXPONENT")
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of PmlDampingFunction" << endl;
	    cout << "DampingPML needs more parameters, for instance :" << endl;
	    cout << "DampingPML = EXPONENT vsigma exponent" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_function = PML_EXPONENT;
	vsigma = to_num<Real_wp>(parameters(1));
	exponent_pml = to_num<Real_wp>(parameters(2));
      }
    else
      {
	type_function = PML_PARABOLE;
	vsigma = to_num<Real_wp>(parameters(0));
	exponent_pml = 2.0;
      }
  }
  
  
  //! Modifies the velocity associated with the PML
  void PmlDampingFunction::SetMaximumVelocity(const Real_wp& v)
  {
    max_velocity = v;
  }

  
  //! specifies the limit of the PML
  void PmlDampingFunction::SetPmlArea(const Real_wp& x0, const Real_wp& xN, const Real_wp& d)
  {
    xmin = x0; xmax = xN; thickness = d;
  }
  

  //! Evaluating PML damping function sigma(dx)
  /*!
    \param[in] dx distance between the interface and the current point of the PML
    \param[in] dsig coefficient of multiplication
    \param[in] a thickness of the PML (somehow maximal value of dx)
    \param[out] zeta damping coefficient sigma(dx)
    \param[out] zeta_primitive primitive of damping coefficient \int sigma(dx)
    Usually the damping function is a parabola with sigma(0) = 0
   */
  void PmlDampingFunction
  ::EvaluateFunctionTau(const Real_wp& dx, const Real_wp& dsig, const Real_wp& a,
			Real_wp& zeta, Real_wp& zeta_primitive) const
  {
    switch (type_function)
      {
      case PML_PARABOLE :
        zeta = dsig*dx*dx;
        zeta_primitive = zeta*dx/3.0;
        break;
      case PML_CONSTANT :
        zeta = dsig*a*a;
        zeta_primitive = zeta*dx;
        break;
      case PML_SHIFTED_PARABOLE :
        zeta = dsig*(dx*dx + offset_damping*a*a);
        zeta_primitive = dsig*dx*(dx*dx/3.0 + offset_damping*a*a);
        break;
      case PML_EXPONENT :
        zeta = dsig*pow(abs(dx), exponent_pml);
        zeta_primitive = zeta*dx / (exponent_pml+1);
        break;
      }
  }
  
  
  //! PML damping coefficients
  /*!
    \param[in] point point where coefficients need to be evaluated
    \param[out] zeta damping coefficient
    \param[out] zeta_primitive primitive of damping coefficient 
   */
  void PmlDampingFunction
  ::GetDampingTau(const Real_wp& point, Real_wp& zeta, Real_wp& zeta_primitive) const
  {
    Real_wp a = this->thickness;
    Real_wp dsig = log(1000e0)*3e0*this->vsigma*this->max_velocity
      /(2.0*pow(a, exponent_pml+1.0));
    
    zeta = 0.0; zeta_primitive = 0.0;
    
    if (point <= xmin)
      EvaluateFunctionTau(point - xmin, dsig, a, zeta, zeta_primitive);
    else if (point >= xmax)
      EvaluateFunctionTau(point - xmax, dsig, a, zeta, zeta_primitive);
  }


  /*****************
   * VarProblem_1D *
   *****************/
  
  
  //! default constructor
  VarProblem_1D::VarProblem_1D()
  {    
    type_coordinate = CARTESIAN;
    nb_unknowns = 1;
    
    condition_left_side = BoundaryConditionEnum::LINE_ABSORBING;
    condition_right_side = BoundaryConditionEnum::LINE_ABSORBING;
    
    nodl = 0;
    omega = 2.0*pi_wp; omega2 = omega*omega;
    
    print_level = 0;
    threshold_matrix = 0;
    mixed_formulation = false;
    
    calcul_enveloppe = false;
    xmin_ = 0.0; xmax_ = 0.0;
    
    nodl_mesh = 0;

#ifdef SELDON_WITH_MPI
    comm_group_mode = MPI_COMM_WORLD;
#endif    
  }


  //! reading the parameters contained in data file
  void VarProblem_1D
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("BoundaryCondition"))
      {
	if (parameters(0) == "PERIODICITY")
	  {
	    this->mesh.SetPeriodicExtremity();
	    condition_left_side = BoundaryConditionEnum::LINE_NEUMANN;
	    condition_right_side = BoundaryConditionEnum::LINE_NEUMANN;
	    return;
	  }
	
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarHelmholtz_1D" << endl;
	    cout << "BoundaryCondition needs more parameters, for instance :" << endl;
	    cout << "BoundaryCondition = DIRICHLET DIRICHLET" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        int nb = 0;
	condition_left_side = ReadBoundaryCondition(nb, parameters, 0);
	condition_right_side = ReadBoundaryCondition(nb, parameters, 1);
      }
    else if (!description_field.compare("OrderDiscretization"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "OrderDiscretization needs more parameters, for instance : " << endl;
            cout << "OrderDiscretization = order" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        int order = to_num<int>(parameters(0));
        this->mesh.SetOrder(order);
      }
    else if (!description_field.compare("InterpolationGrid"))
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "InterpolationGrid needs more parameters, for instance : " << endl;
            cout << "InterpolationGrid = xmin xmax nb_points" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        Real_wp xmin, xmax; int nbPoints;
        xmin = to_num<Real_wp>(parameters(0));
        xmax = to_num<Real_wp>(parameters(1));
        nbPoints = to_num<int>(parameters(2));
        var_section.Init(xmin, xmax, nbPoints);
      }
    else if (!description_field.compare("CoordinatesLaplacian"))
      {
        if (parameters(0) == "POLAR")
          type_coordinate = POLAR;
        else if (parameters(0) == "SPHERICAL")
          type_coordinate = SPHERICAL;
	else
	  type_coordinate = CARTESIAN;
      }
    else if (!description_field.compare("ParametersInterpolationGrid"))
      {
        output_section_param.SetInputData(parameters);
      }
    else if (!description_field.compare("FileInterpolationGrid"))
      {
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "FileInterpolationGrid needs more parameters, for instance : " << endl;
            cout << "FileInterpolationGrid = name_diffracted_field name_total_field" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
	
        output_section_param.SetDiffractedFieldFile(parameters(0));
        output_section_param.SetTotalFieldFile(parameters(1));
	if (parameters.GetM() > 2)
	  output_section_param.SetFileName(2, parameters(2));
      }
    else if (!description_field.compare("FileMesh")) 
      {	
	mesh_data.PushBack(parameters);
      }
    else if (description_field == "AdditionalLayers")
      {
	if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "AdditionalLayers needs more parameters, for instance : " << endl;
            cout << "AdditionalLayers = x nb_layers reference" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
	
	Real_wp pos = to_num<Real_wp>(parameters(0));
	int nb_layers = to_num<int>(parameters(1));
	int ref = to_num<int>(parameters(2));

        mesh.AddLayersInput(pos, nb_layers, ref);
      }
    else if (description_field == "RefinementVertex")
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "RefinementVertex needs more parameters, for instance : " << endl;
            cout << "RefinementVertex = x level ratio" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        Real_wp x = to_num<Real_wp>(parameters(0));
        int lvl = to_num<int>(parameters(1));
        Real_wp ratio = to_num<Real_wp>(parameters(2));
        mesh.AddRefinementVertex(x, lvl, ratio);
      }
    else if (description_field == "MixedFormulation")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_1D" << endl;
	    cout << "MixedFormulation needs more parameters, for instance :" << endl;
	    cout << "MixedFormulation = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (parameters(0) == "YES")
	  mixed_formulation = true;
	else
	  mixed_formulation = false;
      }
    else if (!description_field.compare("CalculEnveloppe"))
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
    else if (!description_field.compare("Frequency"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarHelmholtz_1D" << endl;
	    cout << "Frequency needs more parameters, for instance :" << endl;
	    cout << "Frequency = a b" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	omega = to_num<Real_wp>(parameters(0))*pi_wp*2.0 + to_num<Real_wp>(parameters(1));
	// omega is the pulsation
	// omega2 is the square of omega
	omega2 = omega*omega;
      }
    else if (!description_field.compare("PrintLevel"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_1D" << endl;
	    cout << "PrintLevel needs more parameters, for instance :" << endl;
	    cout << "PrintLevel = level" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
        print_level = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("DampingPML"))
      pml_damping.SetInputData(parameters);
    else if (!description_field.compare("AddPML"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of VarHelmholtz_1D" << endl;
	    cout << "AddPML needs more parameters, for instance :" << endl;
	    cout << "AddPML = YES XYZ delta nb_layers" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("YES"))
	  {
	    int pml_type = Mesh<Dimension1>::PML_NO;
	    
	    if (parameters(1).find("X-",0) != string::npos)
	      pml_type = Mesh<Dimension1>::PML_NEGATIVE_SIDE;
	    else if (parameters(1).find("X+",0) != string::npos)
	      pml_type = Mesh<Dimension1>::PML_POSITIVE_SIDE;
	    else if (parameters(1).find("X",0) != string::npos)
	      pml_type = Mesh<Dimension1>::PML_BOTH_SIDES;
	    
	    Real_wp delta_pml = to_num<Real_wp>(parameters(2));
	    int nb_layers_pml = 0;
	    if (!parameters(3).compare("AUTO"))
	      nb_layers_pml = 0;
	    else
	      to_num(parameters(3), nb_layers_pml);
	    
	    int ref_pml = 0;
	    if (parameters.GetM() > 4)
	      ref_pml = to_num<int>(parameters(4));
	    
	    mesh.SetAdditionPML(pml_type, nb_layers_pml, ref_pml);
	    mesh.SetThicknessPML(delta_pml);
	  }
      }    
  }
  

  int VarProblem_1D
  ::ReadBoundaryCondition(int& nb, const VectString& parameters, int side)
  {
    if (!parameters(nb).compare("DIRICHLET"))
      {
	nb++;
	return BoundaryConditionEnum::LINE_DIRICHLET;
      }
    else if (!parameters(nb).compare("NEUMANN"))
      {
	nb++;
	return BoundaryConditionEnum::LINE_NEUMANN;
      }
    else if ((parameters(nb) == "ABC")
	     || (parameters(nb) == "ABSORBING"))
      {
	nb++;
	return BoundaryConditionEnum::LINE_ABSORBING;
      }
    else if (!parameters(nb).compare("EXACT"))
      {
	nb++;
	return BoundaryConditionEnum::LINE_TRANSPARENT;
      }
    else if (!parameters(nb).compare("IMPEDANCE"))
      {		
	ReadImpedance(nb, parameters, side);     nb++;
	return BoundaryConditionEnum::LINE_IMPEDANCE;
      }
    else
      {
	cout << "Unknown boundary condition : " << parameters(nb) << endl;
	abort();
      }
    
    nb++;
    return BoundaryConditionEnum::LINE_INSIDE;
  }
  

  //! computes position of quadrature points as required by indexes
  void VarProblem_1D
  ::ComputeQuadraturePoints(VectReal_wp& pts)
  {
    VectReal_wp s(2);
    int nb = 0;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(i);
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
	nb += nb_points_quadrature + 2;
      }
    
    pts.Reallocate(nb);
    nb = 0;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(i);
	
	// number of integration points
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
	
	// We get the two extremities of the edge
	this->mesh.GetVerticesElement(i, s);
        
	for (int j = 0; j <= nb_points_quadrature+1; j++)
	  {
	    Real_wp x(0);
	    if (j == 0)
	      x = s(0);
	    else if (j == nb_points_quadrature+1)
	      x = s(1);
	    else
	      x = Fb.Points(j-1)*s(1) + (1.0-Fb.Points(j-1))*s(0);
	    
	    pts(nb) = x; nb++;
	  }
      }
  }


  //! retrieves the mass matrix without physical coefficient
  void VarProblem_1D::GetIntegralMass(Vector<Real_wp>& M)
  {
    VectReal_wp s(2);
    M.Reallocate(this->nodl);
    M.Zero();
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(i);
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();

	// We get the two extremities of the edge
	this->mesh.GetVerticesElement(i, s);

	// length of the edge
	Real_wp Ji = abs(s(1)-s(0));
	
	for (int j = 0; j < nb_points_quadrature; j++)
	  {
	    M(this->mesh.GetNumberDof(i, j)) += Ji*Fb.Weights(j);
	    if (this->mixed_formulation)
	      M(this->mesh.GetNbDof() + this->OffsetDofV(i) + j) = Ji*Fb.Weights(j);
	  }
      }
  }
  

  //! returns the stiffness matrix for an element of the mesh
  template<class MatrixB>
  void VarProblem_1D::GetStiffnessMatrix(MatrixB& mat_stiff)
  {
    if (mesh.GetNbElt() <= 0)
      return;
    
    // checking if the mesh is uniform
    int N = mesh.GetNbElt();
    Real_wp xmin = mesh.Vertex(0), xmax = mesh.Vertex(N);
    Real_wp h = (xmax - xmin) / N;
    Real_wp threshold = 1e6*h*epsilon_machine;
    int r = mesh.GetOrderElement(0);
    for (int i = 1; i < mesh.GetNbElt(); i++)
      {
        Real_wp hi = mesh.Vertex(i+1) - mesh.Vertex(i);
        if ( (abs(hi - h) > threshold) || (r != mesh.GetOrderElement(i)))
          {
            cout << "The mesh should be uniform in order to put 1/h in stiffness block" << endl;
            abort();
          }
      }
    
    const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(0);
    Real_wp coef = 1.0/h;
    const Matrix<Real_wp, Symmetric, RowSymPacked>& stiffness_matrix = Fb.GetStiffnessMatrix();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        mat_stiff(i, j) = coef*stiffness_matrix(i, j);
  }


  //! returns the gradient matrix for an element of the mesh
  template<class MatrixB>
  void VarProblem_1D::GetGradientMatrix(MatrixB& mat_stiff)
  {
    if (mesh.GetNbElt() <= 0)
      return;
    
    int r = mesh.GetOrderElement(0);
    const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(0);
    const Matrix<Real_wp>& gradient_matrix = Fb.GetGradientMatrix();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        mat_stiff(i, j) = gradient_matrix(i, j);
    
  }
  

  //! returns the local mass matrix for an element of the mesh
  template<class VectorB>
  void VarProblem_1D::GetLocalMassMatrix(VectorB& mat_mass)
  {
    if (mesh.GetNbElt() <= 0)
      return;

    int r = mesh.GetOrderElement(0);
    int N = mesh.GetNbElt();
    Real_wp xmin = mesh.Vertex(0), xmax = mesh.Vertex(N);
    Real_wp h = (xmax - xmin) / N;
    const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(0);
    Real_wp coef = h;
    for (int i = 0; i <= r; i++)
      mat_mass(i) = coef*Fb.Weights(i);    
  }
  

  //! construction of mesh and finite element
  void VarProblem_1D::ComputeMeshAndFiniteElement(const string& name_elt, bool split_mesh)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    // constructing mesh (for all the processors)
    if (this->mesh_data.GetM() > 0)
      this->mesh.ConstructMesh(this->mesh_data(0), xmin_, xmax_);
    
    // computing maximal speed in PML layers
    int ref_max = 0;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (mesh.Element(i).IsPML())
	ref_max = max(ref_max, mesh.Element(i).GetReference());
    
    Real_wp max_velocity_pml = 0.0;
    Vector<bool> RefUsed(ref_max+1);
    RefUsed.Fill(false);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (mesh.Element(i).IsPML())
        {
          int ref = mesh.Element(i).GetReference();
          RefUsed(ref) = true;
        }
    
    for (int ref = 0; ref <= ref_max; ref++)
      if (RefUsed(ref))
        {
          Real_wp velocity = this->GetVelocityOfMedia(ref);
          max_velocity_pml = max(max_velocity_pml, velocity);
        }
    
    //DISP(RefUsed); DISP(max_velocity_pml);
    pml_damping.SetPmlArea(this->GetXmin(), this->GetXmax(), this->mesh.GetThicknessPML());
    pml_damping.SetMaximumVelocity(max_velocity_pml);
    
    // constructing finite elements
    this->ClearFiniteElement();
    int order = this->mesh.GetOrder();
    // For Radau points, we need to construct Lobatto elements
    // and affect these elements to elements outside the axis
    string name_element = name_elt;
    if (name_elt == "EDGE_RADAU")
      name_element = "EDGE_LOBATTO";

    Vector<bool> change_elt(this->mesh.GetNbElt());
    change_elt.Fill(true);
    
    this->AddFiniteElement1D(name_element, order, change_elt, this->mesh);    

    change_elt.Fill(false); change_elt(0) = true;
    if (name_elt == "EDGE_RADAU")
      this->AddFiniteElement1D(name_elt, order, change_elt, this->mesh);
    
    // mesh is numbered
    this->mesh.NumberMesh();
    nodl = this->mesh.GetNbDof();
    
    if (mixed_formulation)
      {
	OffsetDofV.Reallocate(this->mesh.GetNbElt()+1);
	OffsetDofV(0) = 0;
	for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  OffsetDofV(i+1) = OffsetDofV(i) + this->mesh.GetOrderElement(i)+1;

	nodl += OffsetDofV(this->mesh.GetNbElt());
      }        
    
    // computation of interpolation grid (needed to display solution u)
    this->var_section.LocalizePoints(this->mesh);
    if (rank_proc == 0)
      if (this->print_level >= 1)
        cout<<"Interpolation Section OK"<<endl;
    
    // restricting computations for parallel computations
    nodl_mesh = 0;
    var_section.SetGlobalNumberPoints(var_section.GetNbPointsGrid());
    
#ifdef SELDON_WITH_MPI
    if ((nb_proc > 1) && (split_mesh))
      {
        int nb_elt = this->mesh.GetNbElt();
        IVect Epart(nb_elt);
        // last processor has less elements
        int nb_elt_per_proc = nb_elt / nb_proc;
        if (nb_proc*nb_elt_per_proc != nb_elt)
          nb_elt_per_proc++;
        
        for (int i = 0; i < nb_elt; i++)
          {
            int num_proc = i/nb_elt_per_proc;
            Epart(i) = num_proc;
            // LocalEpart(i) = i%nb_elt_per_proc;
          }
        
        // creating the mesh extracted from the global mesh
        int n0 = rank_proc*nb_elt_per_proc;
        int n1 = min(nb_elt, (rank_proc+1)*nb_elt_per_proc);
        Mesh<Dimension1> new_mesh;
        this->mesh.CreateSubmesh(new_mesh, n0, n1);
        
        nodl_mesh = nodl;
        nodl = new_mesh.GetNbDof();
        
        // changing boundary conditions
        if (rank_proc > 0)
          condition_left_side = BoundaryConditionEnum::LINE_NEIGHBOR;
        
        if (rank_proc < nb_proc-1)
          condition_right_side = BoundaryConditionEnum::LINE_NEIGHBOR;
        
        // updating mesh
        this->mesh = new_mesh;
        
        if (mixed_formulation)
          {
            OffsetDofV.Reallocate(this->mesh.GetNbElt()+1);
            OffsetDofV(0) = 0;
            for (int i = 0; i < this->mesh.GetNbElt(); i++)
              OffsetDofV(i+1) = OffsetDofV(i) + this->mesh.GetOrderElement(i)+1;
            
            nodl += OffsetDofV(this->mesh.GetNbElt());
          }

        // parallel connectivity
        if (rank_proc > 0)
          {
            OverlapDofNumber_Subdomain.Reallocate(1);
            OverlapDofNumber_Subdomain(0) = 0;
          }
        else
          OverlapDofNumber_Subdomain.Clear();
        
        // interpolation grid
        int nb_pts_glob = var_section.GetNbPointsGrid();
        Vector<IVect> ProcNumberPoints(nb_proc);
        IVect nb_points(nb_proc);
        nb_points.Fill(0);
        for (int i = 0; i < var_section.GetNbPointsGrid(); i++)
          {
            int num_elem = var_section.GetElementNumber(i);
            if (num_elem >= 0)
              nb_points(Epart(num_elem))++;
          }
        
        for (int i = 0; i < nb_proc; i++)
          ProcNumberPoints(i).Reallocate(nb_points(i));
        
        GridInterpolation<Dimension1> new_section;
        new_section.ReallocatePoints(nb_points(rank_proc));
        nb_points.Fill(0);
        for (int i = 0; i < var_section.GetNbPointsGrid(); i++)
          {
            int num_elem = var_section.GetElementNumber(i);
            if (num_elem >= 0)
              {
                int proc = Epart(num_elem);
                ProcNumberPoints(proc)(nb_points(proc)) = i;
                if (proc == rank_proc)
                  {
                    new_section.SetElementNumber(nb_points(proc), num_elem-n0);
                    new_section.SetLocalCoordinate(nb_points(proc),
                                                   var_section.GetLocalCoordinate(i));
                    new_section.SetGlobalCoordinate(nb_points(proc),
                                                    var_section.GetGlobalCoordinate(i));
                  }
                
                nb_points(proc)++;
              }
          }
        
        new_section.SetGlobalNumberPoints(nb_pts_glob);
        var_section = new_section;
        if (rank_proc == 0)
          var_section.ProcNumberPoints = ProcNumberPoints;
      }
#endif
  }
  
  
  //! assembles a vector 
  template<class T>
  void VarProblem_1D::AddDomains(Vector<T>& X, int nb_u) const
  {
    if (nodl_mesh <= 0) 
      return;
    
    if (nb_u == -1)
      nb_u = nb_unknowns;

#ifdef SELDON_WITH_MPI
    int tag = 15;
    const MPI_Comm& comm = comm_group_mode;
    int nb_proc; MPI_Comm_size(comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_group_mode, &rank_proc);

    Vector<T> xsend_left(nb_u), xsend_right(nb_u), xrecv_left(nb_u), xrecv_right(nb_u);
    MPI_Request request_send_left, request_send_right;
    Vector<int64_t> xsend_left_tmp, xsend_right_tmp;
    Vector<int64_t> xrecv_left_tmp, xrecv_right_tmp;
    // sending values on left and right
    if (rank_proc > 0)
      {
        for (int i = 0; i < nb_u; i++)
          xsend_left(i) = X(nodl*i);
        
        request_send_left = MpiIsend(comm, xsend_left, xsend_left_tmp,
                                     nb_u, rank_proc-1, tag);
      }
    
    if (rank_proc < nb_proc-1)
      {
        for (int i = 0; i < nb_u; i++)
          xsend_right(i) = X(nodl*(i+1)-1);
        
        request_send_right = MpiIsend(comm, xsend_right, xsend_right_tmp,
                                      nb_u, rank_proc+1, tag+1);
      }
    
    // receiving these values
    MPI_Status status;
    if (rank_proc > 0)
      MpiRecv(comm, xrecv_left, xrecv_left_tmp, nb_u,
              rank_proc-1, tag+1, status);
    
    if (rank_proc < nb_proc-1)
      MpiRecv(comm, xrecv_right, xrecv_right_tmp, nb_u,
              rank_proc+1, tag, status);
    
    // waiting end of communications
    if (rank_proc > 0)
      MPI_Wait(&request_send_left, &status);

    if (rank_proc < nb_proc-1)
      MPI_Wait(&request_send_right, &status);

    // then summing values
    if (rank_proc > 0)
      for (int i = 0; i < nb_u; i++)
        X(nodl*i) += xrecv_left(i);

    if (rank_proc < nb_proc-1)
      for (int i = 0; i < nb_u; i++)
        X(nodl*(i+1)-1) += xrecv_right(i);
#endif
  }
  

  //! we add to b_source(i)  \f$ \int f \varphi_i + \int g \varphi'_i \f$
  template<class T>
  void VarProblem_1D
  ::AddVolumetricSource(Vector<T> & b_source, VirtualSourceField<T, Dimension1>& f)
  { 
    VectReal_wp s(2), points;
    Vector<T> feval, contrib;
    
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	// we get vertices of the element i
	this->mesh.GetVerticesElement(i, s);
	const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(i);
	
	// number of integration points
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
        int nb_dof_elt = Fb.GetNbDof();
	Real_wp h = abs(s(1) - s(0));
        
	Fb.FjElem(s, points);
	
        feval.Reallocate(nb_points_quadrature);
        contrib.Reallocate(nb_dof_elt);
        feval.Fill(0);
        
	// if f is non null
	// we evaluate f at quadrature points
	for (int j = 0; j < nb_points_quadrature; j++)
	  {
	    f.EvaluateFunction(points(j), feval(j));
	    feval(j) *= h;
	  }
	
	// we compute \int f \phi_i
        Fb.ComputeIntegralRef(feval, contrib);
        
	// we add contributions to b_source
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = this->mesh.GetNumberDof(i, j);
	    b_source(num_dof) += contrib(j);
	  }
      }
  }
  
  
  //! surfacic source is added
  /* template<class T>
  void VarProblem_1D
  ::AddSurfacicSource(Vector<T>& b_source, VirtualSourceFEM<T, Dimension1>& f)
  { 
    T feval;
    if (f.IsNonNullLeftSource(condition_left_side))
      {
	f.EvaluateSurfacicSource(this->GetXmin(), feval);
	b_source(0) += feval;
      }
    
    if (f.IsNonNullRightSource(condition_right_side))
      {
	f.EvaluateSurfacicSource(this->GetXmax(), feval);
	b_source(this->mesh.GetNbDof()-1) += feval;
      }
  }
  
  
  //! projection of volumetric source
  template<class T>
  void VarProblem_1D
  ::AddVolumetricProjection(const T& alpha, Vector<T>& proj, VirtualProjectorFEM<T, Dimension1>& f)
  {
    VectReal_wp s(2), points;
    VectBool DofUsed(this->mesh.GetNbDof()); DofUsed.Fill(false);
    Complexe feval;
    
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	// we get the vertices of the face
	this->mesh.GetVerticesElement(i, s);
	
        // number of dofs in the face
        int order = this->mesh.GetOrderElement(i);
        const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(i);
	int nb_dof_elt = Fb.GetNbDof();
	
	// we make the transformation Fi on quadrature points
	Fb.FjElemDof(s, points);
	
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = this->mesh.GetNumberDof(i,j);
	    if (!DofUsed(num_dof))
	      {
		DofUsed(num_dof) = true;
		f.EvaluateFunction(points(j), feval);
                proj(num_dof) += alpha*feval;
	      }
	  }
      }
  }*/
  
  
  //! computation of right hand side
  template<class T>
  void VarProblem_1D::ComputeRightHandSide(Vector<T>& bsrc)
  {
    bsrc.Reallocate(this->mesh.GetNbDof()); bsrc.Fill(0);
    
    //GaussianSourceField<T, Dimension1> fsrc;
    //this->AddVolumetricSource(bsrc, fsrc);
    //this->AddVolumetricProjection(1.0, bsrc, fsrc);
    //this->AddSurfacicSource(bsrc, fsrc);
    
    // dirac
    SetComplexOne(bsrc(0));
    if (type_coordinate == POLAR)
      Mlt(1.0/(2*pi_wp), bsrc);    
    else if (type_coordinate == SPHERICAL)
      Mlt(1.0/(4*pi_wp), bsrc);
  }
  
  
  //! Dirichlet condition is treated
  void VarProblem_1D::TreatDirichletCondition()
  {
    if ((condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
        &&(condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET))
      {
	this->Dirichlet_dof.Reallocate(2);
	this->Dirichlet_dof(0) = 0;
	this->Dirichlet_dof(1) = this->mesh.GetNbDof() - 1;
      }
    else if (condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	this->Dirichlet_dof.Reallocate(1);
	this->Dirichlet_dof(0) = 0;
      }
    else if (condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	this->Dirichlet_dof.Reallocate(1);
	this->Dirichlet_dof(0) = this->mesh.GetNbDof() - 1;
      }    
    else
      this->Dirichlet_dof.Clear();
  }
  

  //! computation of finite element matrix (without boundary conditions)
  template<class T>
  void VarProblem_1D
  ::AddMatrixFEM(VirtualMatrix<T>& mat_direct, GlobalGenericMatrix<T>& nat_mat)
  {
    Matrix<T> mat_elt;
    CondensationBlockSolver_Base<T> solver;
    AssembleMatrix(mat_direct, mat_elt, nat_mat, *this, solver, 0, 0);
  }

  
  //! solution is interpolated and written in output file
  template<class T>
  void VarProblem_1D::WriteDatas(Vector<T>& U0)
  {    
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    int N = this->var_section.GetNbPointsGrid();
    Vector<T> output_vector(N);
    this->ComputeInterpolationU(U0, this->var_section, output_vector);
    if ((this->print_level >= 0) && (rank_proc == 0))
      cout<<"Writing"<<endl;
    
    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
    bool ascii = false;

    GridInterpolationFull<Dimension1> grid;
    int type = GridInterpolationFull<Dimension1>::LINE;
    grid.Init(type, this->var_section.GetXmin(), this->var_section.GetXmax(), N);

    string name_file = this->output_section_param.GetTotalFieldFile();
    WriteMatlab(output_vector, grid, name_file, double_prec, ascii);
    
    if (this->output_section_param.GetNbFile() > 2)
      {
	this->ComputeInterpolationU(U0, this->var_section, output_vector,
				    false, true);
	
	name_file = this->output_section_param.GetFileName(2);
	WriteMatlab(output_vector, grid, name_file, double_prec, ascii);
      }    
  }
  
  
  //! computation of the solution on a regular grid
  template<class T>
  void VarProblem_1D
  ::ComputeInterpolationU(const Vector<T> & U0, const GridInterpolation<Dimension1> & var_interp,
                          Vector<T>& trace_u, bool discontinuous, bool grad) const
  {
    trace_u.Reallocate(var_interp.GetNbPointsGrid());
    T scal_u; Real_wp point_loc;
    VectReal_wp val_phi;
    
    for (int i = 0; i < var_interp.GetNbPointsGrid(); i++)
      {	
	int iquad = var_interp.GetElementNumber(i);
	point_loc = var_interp.GetLocalCoordinate(i);
	scal_u = GetInterpolate(U0, iquad, point_loc, val_phi, discontinuous, grad);
	trace_u(i) = scal_u;
      }
    
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm_group_mode, &rank_proc);
    if (nb_proc <= 1)
      return;
    
    const MPI_Comm& comm = comm_group_mode;
    MPI_Status status;
    Vector<int64_t> trace_tmp;
    if (rank_proc != 0)
      MpiSsend(comm, trace_u, trace_tmp, trace_u.GetM(), 0, 24);
    else
      {
        Vector<T> trace_u_loc(trace_u);
        T zero; SetComplexZero(zero);
        trace_u.Reallocate(var_interp.GetGlobalNumberPoints());
        trace_u.Fill(zero);
        Vector<T> u_proc;
        for (int i = 0; i < nb_proc; i++)
          {
            if (i == 0)
              u_proc = trace_u_loc;
            else
              {
                int N = var_interp.ProcNumberPoints(i).GetM();
                u_proc.Reallocate(N);
                MpiRecv(comm, u_proc, trace_tmp, N, i, 24, status);
              }
            
            for (int j = 0; j < u_proc.GetM(); j++)
              trace_u(var_interp.ProcNumberPoints(i)(j)) = u_proc(j);
          }
      }
#endif
  }  
  

  //! performs interpolation on a single point of an element
  template<class T>
  T VarProblem_1D
  ::GetInterpolate(const Vector<T>& U0, int iquad, const Real_wp& point_loc,
		   VectReal_wp& val_phi, bool discontinuous, bool grad) const
  {    
    T scal_u;
    SetComplexZero(scal_u);
    if ((iquad>=0) && (iquad < this->mesh.GetNbElt()))
      {
	const ElementReference_Dim<Dimension1>& Fb = this->GetReferenceElement(iquad);
	int nb_dof_loc = Fb.GetNbDof();
	val_phi.Reallocate(nb_dof_loc);
	if (grad)
	  {
	    TinyVector<Real_wp, 2> s;
	    Fb.ComputeGradientPhiRef(point_loc, val_phi);
	    this->mesh.GetVerticesElement(iquad, s);
	    Real_wp h = abs(s(1) - s(0));
	    Mlt(1.0/h, val_phi);
	  }
	else
	  Fb.ComputeValuesPhiRef(point_loc, val_phi);
	
	if (discontinuous)
	  for (int j = 0; j < nb_dof_loc; j++)
	    scal_u += U0(OffsetDofV(iquad) + j)*val_phi(j);
	else
	  for (int j = 0; j < nb_dof_loc; j++)
	    scal_u += U0(this->mesh.GetNumberDof(iquad, j))*val_phi(j);
      }
    
    return scal_u;
  }

}

#define MONTJOIE_FILE_VAR_PROBLEM_1D_CXX
#endif
