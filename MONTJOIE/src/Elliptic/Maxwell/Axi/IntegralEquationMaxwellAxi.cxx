#ifndef MONTJOIE_FILE_INTEGRAL_EQUATION_MAXWELL_AXI_CXX

namespace Montjoie
{
  //! default constructor
  template<class TypeElement, class TypeEquation>
  VarIntegralEquationAxi<TypeElement, TypeEquation>
  ::VarIntegralEquationAxi(EllipticProblem<TypeElement, TypeEquation>& var)
    : var_harmonic(var)
  {
    name_file_Jb = "total_Jb"; name_file_Jt = "total_Jt"; name_file_Jnorme = "total_Jn";
    name_file_Kb = "total_Kb"; name_file_Kt = "total_Kt"; name_file_Knorme = "total_Kn";
    type_output_file_ie = 0;
    order_integration_phi = 0;
    order_integration_radiusX = 0; order_integration_radiusY = 0;
    automatic_integration_phi = true; automatic_integration_radius = true;
    nb_edges_integral = 0; nodl_integral = 0; boundary_condition = var.LINE_DIRICHLET;
    nb_dof_Jt = 0; nb_dof_Jb = 0; nodl_Jt = 0; nodl_Jb = 0;
    nbPoints_grid_boundary_ie = 200; space_step_IE = 0;
  }
  
  
  //! reading the line of a data file
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement, TypeEquation>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    int order = var_harmonic.mesh.GetOrder();
    if (!description_field.compare("OrderIntegrationPhi"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
	    cout << "OrderIntegrationPhi needs at least one parameter, for instance :" << endl;
	    cout << "OrderIntegrationPhi = AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("AUTO"))
	  automatic_integration_phi = true;
	else
	  {
	    automatic_integration_phi = false;
	    to_num(parameters(0), order_integration_phi);
	  }
      }
    else if (!description_field.compare("DisplayCurrents"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
	    cout << "DisplayCurrents needs at least three parameters, for instance :" << endl;
	    cout << "DisplayCurrents = 3D ASCII racine_fichier" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("2D"))
	  {	    
	    if (parameters.GetM() <= 0)
	      {
		cout << "In SetInputData of VarIntegralEquationAxi" << endl;
		cout << "DisplayCurrents needs 4 parameter, for instance :" << endl;
		cout << "DisplayCurrents = 2D ASCII N racine" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    if (!parameters(1).compare("ASCII"))
	      type_output_file_ie = 1;
	    else
	      type_output_file_ie = 0;
	    
	    to_num(parameters(2), nbPoints_grid_boundary_ie);
	    name_file_Jt = parameters(3) + string("_Jt.dat");
	    name_file_Jb = parameters(3) + string("_Jb.dat");
	    name_file_Kt = parameters(3) + string("_Kt.dat");
	    name_file_Kb = parameters(3) + string("_Kb.dat");
	    name_file_Jnorme = parameters(3) + string("_Jn.dat");
	    name_file_Knorme = parameters(3) + string("_Kn.dat");
	  }
	if (!parameters(0).compare("3D"))
	  {
	    if (!parameters(1).compare("ASCII"))
	      type_output_file_ie = 3;
	    else
	      type_output_file_ie = 2;
	    
	    name_file_Jt = parameters(2) + string("_Jt");
	    name_file_Jb = parameters(2) + string("_Jb");
	    name_file_Jnorme = parameters(2) + string("_Jn");
	    name_file_Knorme = parameters(2) + string("_Kn");
	    this->name_file_mesh3D_ie = parameters(2) + string("_J.mesh");
	  }
      }
    else if (!description_field.compare("OrderIntegrationAxis"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
	    cout << "OrderIntegrationAxis needs at least one parameter, for instance :" << endl;
	    cout << "OrderIntegrationAxis = AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("AUTO"))
	  {
	    order_integration_radiusX = 3*order+1;
	    order_integration_radiusY = 3*order+1;
	    automatic_integration_radius = true;
	  }
	else
	  {
	    if (parameters.GetM() <= 0)
	      {
		cout << "In SetInputData of VarIntegralEquationAxi" << endl;
		cout << "OrderIntegrationAxis needs two parameters, for instance :" << endl;
		cout << "OrderIntegrationAxis = rX rY" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }
	    
	    automatic_integration_radius = false;
	    to_num(parameters(0), order_integration_radiusX);
	    to_num(parameters(1), order_integration_radiusY);
	  }
      }
    else if (!description_field.compare("IntegralEquation"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
	    cout << "IntegralEquation needs at least two parameters, for instance :" << endl;
	    cout << "IntegralEquation = ONLY AUTO" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// first parameter is ONLY, YES or NO
	// ONLY -> a pure integral equation is used (no FEM)
	// YES  -> an integral equation is used, coupled with a finite element method
	// NO   -> no, integral equation, only the finite element method is used
	if (!parameters(0).compare("ONLY"))
	  {
	    this->coupling_integral_equation = this->BIE_ONLY;
	    if (!parameters(1).compare("AUTO"))
	      {
		this->order_integration_regular = order;
		this->order_integrationX = 2*order+1;
		this->order_integrationY = 2*order+1;
	      }
	    else
	      {
		if (parameters.GetM() <= 3)
		  {
		    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
		    cout << "IntegralEquation needs at 4 parameters, for instance :" << endl;
		    cout << "IntegralEquation = ONLY r0 rX rY" << endl;
		    cout << "Current parameters are : " << endl << parameters << endl;
		    abort();
		  }

		to_num(parameters(1), this->order_integration_regular);
		to_num(parameters(2), this->order_integrationX);
		to_num(parameters(3), this->order_integrationY);
	      }
	  }
	else if (!parameters(0).compare("YES"))
	  {
	    this->coupling_integral_equation = this->BIE_FEM;
	    if (!parameters(1).compare("AUTO"))
	      {
		this->order_integration_regular = order;
		this->order_integrationX = 2*order+1;
		this->order_integrationY = 2*order+1;
	      }
	    else
	      {
		if (parameters.GetM() <= 3)
		  {
		    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
		    cout << "IntegralEquation needs 4 parameters, for instance :" << endl;
		    cout << "IntegralEquation = ONLY r0 rX rY" << endl;
		    cout << "Current parameters are : " << endl << parameters << endl;
		    abort();
		  }
		
		to_num(parameters(1), this->order_integration_regular);
		to_num(parameters(2), this->order_integrationX);
		to_num(parameters(3), this->order_integrationY);
	      }
	  }
	else
	  this->coupling_integral_equation = this->FEM_ONLY;
      }
    else if (!description_field.compare("ConditionIntegralEquation"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
	    cout << "ConditionIntegralEquation needs at least one parameter, for instance :"
                 << endl;
	    cout << "ConditionIntegralEquation = DIRICHLET" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (!parameters(0).compare("DIRICHLET"))
	  this->condition_integral_equation = var_harmonic.LINE_DIRICHLET;
	else if (!parameters(0).compare("HIGH_CONDUCTIVITY"))
	  this->condition_integral_equation = var_harmonic.LINE_HIGH_CONDUCTIVITY;
	else if (!parameters(0).compare("IMPEDANCE"))
	  {
	    if (parameters.GetM() <= 1)
	      {
		cout << "In SetInputData of VarIntegralEquationAxi" << endl;
		cout << "ConditionIntegralEquation needs 2 parameters, for instance :" << endl;
		cout << "ConditionIntegralEquation = IMPEDANCE beta" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    this->condition_integral_equation = var_harmonic.LINE_IMPEDANCE;
	    to_num(parameters(1), this->beta_imped_integral);
	  }
      }
    else if (!description_field.compare("IntegralFormulation"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
	    cout << "IntegralFormulation needs 2 parameters, for instance :" << endl;
	    cout << "IntegralFormulation = type_singular_integration EFIE" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	to_num(parameters(0), this->singular_integration);
	if (!parameters(1).compare("EFIE"))
	  {
	    this->type_integral_formulation = this->EFIE_FORMULATION;
	    this->alpha_efie = 1.0; this->beta_mfie = 0.0;
	  } 
	else if (!parameters(1).compare("MFIE"))
	  {
	    this->type_integral_formulation = this->MFIE_FORMULATION;
	    this->alpha_efie = 0.0; this->beta_mfie = 1.0;
	  }
	else if (!parameters(1).compare("CFIE"))
	  {
	    this->type_integral_formulation = this->CFIE_FORMULATION;
	    if ( (parameters.GetM() < 3) || (!parameters(2).compare("AUTO")) )
	      {
		this->alpha_efie = 0.5;
		this->beta_mfie = 0.5;
	      }
	    else
	      {
		if (parameters.GetM() <= 3)
		  {
		    cout << "In SetInputData of VarIntegralEquationAxi" << endl;
		    cout << "IntegralFormulation needs 4 parameters, for instance :" << endl;
		    cout << "IntegralFormulation = type_singular_integration CFIE alpha beta"
                         << endl;
		    cout << "Current parameters are : " << endl << parameters << endl;
		    abort();
		  }
		
		this->alpha_efie = to_num<Complexe>(parameters(2));
		this->beta_mfie = to_num<Complexe>(parameters(3));
	      }
	  }
      }
  }
  
  
  //! constructs surfacics dofs associated with a boundary condition
  /*!
    Used to extract surfacic dofs for integral equation
    \param[in] ref boundary condition number 
    \param[out] mesh_surf surfacic mesh
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement, TypeEquation>
  ::GetDofsBoundary(int ref, SurfacicMesh<Dimension2>& mesh_surf)
  {
    boundary_condition = ref;
    Mesh<Dimension2>& mesh = var_harmonic.mesh;
    int r = mesh.GetOrder();
    nb_dof_Jt = r + 1;
    nb_dof_Jb = r;
    int nb_edges_ref = mesh.GetNbBoundaryRef();
    const TypeElement& Fb = var_harmonic.GetFaceBasis2(r);
    // first step
    // we get the list of the edges, which boundary condition is equal to ref
    // we get also the dofs
    // IndexEdgeInt_to_EdgeRef(i) will be the list of numbers of the edges 
    // where integral equation is set
    // IndexDofInt_to_DofVol(i) will be the list of numbers of volumics dofs
    // taken in account to couple volume/integral equation
    // PriorityNumDof  is equal to [0 1......N]
    // ListePoint : list of vertices where integral equation is set
    int num_dof, num_dof_loc; int ind = 0, nb = 0; nodl_integral = 0;
    IndexEdgeInt_to_EdgeRef.Reallocate(nb_edges_ref); 
    IndexDofInt_to_DofVol.Reallocate(nb_edges_ref*(2*r+1));
    IVect PriorityNumDof(nb_edges_ref*(2*r+1));
    IVect ListePoint(nb_edges_ref*2);
    // loop on all referenced edges of the volumic mesh
    for (int i = 0; i < nb_edges_ref; i++)
      {
	int ref_edge = mesh.BoundaryRef(i).GetReference();
	// we found an edge with the good boundary condition
	if ( ref == mesh.GetBoundaryCondition(ref_edge) )
	  {	   
	    // adding this edge number 
	    IndexEdgeInt_to_EdgeRef(ind) = i;
	    int num_edge = i;
	    int iquad = mesh.Boundary(num_edge).numElement(0);
	    int num_loc = mesh.Element(iquad).GetPositionBoundary(num_edge);
	    // adding the two vertices of the edge in array ListePoint
	    ListePoint(nb++) = mesh.GetEdge(num_edge).numVertex(0);
	    ListePoint(nb++) = mesh.GetEdge(num_edge).numVertex(1);
	      
	    // loop on scalar dofs (E_\theta)
	    for (int j = 0; j <= r; j++)
	      {
		num_dof_loc = Fb.GetLocalNumberScal(num_loc, j);
		IndexDofInt_to_DofVol(nodl_integral)
                  = mesh.Element(iquad).GetNumberDof(num_dof_loc);
		PriorityNumDof(nodl_integral) = nodl_integral;
		nodl_integral++;
	      }
	    ind++;
	  }
      }
    
    mesh.Write("decoup.mesh");
    // loop on all referenced edges of the volumic mesh
    // now we get vectorial dofs (E)
    for (int i = 0; i < nb_edges_ref; i++)
      {
	int ref_edge = mesh.BoundaryRef(i).GetReference();
	if ( ref == mesh.GetBoundaryCondition(ref_edge) )
	  {
	    int num_edge = i;
	    int iquad = mesh.Boundary(num_edge).numElement(0);
	    int num_loc = mesh.Element(iquad).GetPositionBoundary(num_edge);
	    
	    for (int j = 0; j < r; j++)
	      {
		num_dof_loc = Fb.GetLocalNumberVec(num_loc,j);
		IndexDofInt_to_DofVol(nodl_integral)
                  = mesh.Element(iquad).GetNumberDof(num_dof_loc);
		PriorityNumDof(nodl_integral) = nodl_integral;
		nodl_integral++;
	      }
	  }
      }

    // we sort dofs in order to have "Jt" dofs numbered before "Jb" dofs
    RemoveDuplicate(nodl_integral, IndexDofInt_to_DofVol, PriorityNumDof);
    Sort(nodl_integral, PriorityNumDof, IndexDofInt_to_DofVol);
    IndexDofInt_to_DofVol.Resize(nodl_integral);
    IndexEdgeInt_to_EdgeRef.Resize(ind);
    // sorting point
    Assemble(nb, ListePoint); ListePoint.Resize(nb);
    
    IVect Inverse_ListePoint(mesh.GetNbVertices()); Inverse_ListePoint.Fill(-1);
    for (int i = 0; i < nb; i++)
      Inverse_ListePoint(ListePoint(i)) = i;
    
    // we search a good initial point (the lowest vertex on the axis !)
    int num_init = -1;
    R2 point_init;
    for (int i = 0; i < nb; i++)
      {
	if (var_harmonic.Vertex_On_Axe(ListePoint(i))||(num_init==-1))
	  {
	    if (num_init == -1)
	      {
		point_init = mesh.Vertex(ListePoint(i));
		num_init = i;
	      }
	    else if (mesh.Vertex(ListePoint(i))(1) < point_init(1))
	      {
		point_init = mesh.Vertex(ListePoint(i));
		num_init = i;
	      }
	  }
      }
    
    
    if (var_harmonic.print_level >= 0)
      cout<<"Initial point of the surfacic mesh "<<num_init<<endl<<point_init<<endl;
    
    // if no initial point is found, we exit the program
    // the only luck to have this case is that no referenced edges have the right BC
    if (num_init == -1)
      {
	nodl_integral = 0;
	nb_edges_integral = 0;
	cout<<"Is there a referenced edge with the boundary condition "
	    <<ref<<" inside the mesh ?"<<endl;
	abort();
      }
    // DISP(mesh.edges_ref); DISP(mesh_surf.IndexEdgeInt_to_EdgeRef);
    
    // loop to recover all edges from this initial point
    // here we consider, that the integral equation is set on a curve (only one object)
    int iback = -1, ipoint, inext;
    int num_edge;
    ipoint = ListePoint(num_init); num_init = ipoint;
    bool test_loop = true; nb = 0;
    while (test_loop)
      {
	// next point on the boundary condition
	inext = mesh.FindFollowingVertex(iback, ipoint, ref,
                                         mesh.GetBoundaryCondition(), num_edge); 
	// DISP(iback); DISP(ipoint); DISP(inext); DISP(num_edge); DISP(ref);
	if ((inext != -1)&&((ipoint != num_init)||(nb==0)))
	  {
	    // we need referenced edge number and not global edge number
	    IndexEdgeInt_to_EdgeRef(nb) = num_edge;
	    nb++;
	  }
	else
	  test_loop = false;
	
	iback = ipoint;
	ipoint = inext;
	// DISP(iback); DISP(ipoint);
      }
    // DISP(mesh_surf.IndexEdgeInt_to_EdgeRef); DISP(nb);
    
    Inverse_DI_DV.Reallocate(mesh.GetNbDof()); Inverse_DI_DV.Fill(-1);
    for (int i = 0; i < nodl_integral; i++)
      Inverse_DI_DV(IndexDofInt_to_DofVol(i)) = i;
    
    Inverse_EI_ER.Reallocate(mesh.GetNbBoundaryRef()); Inverse_EI_ER.Fill(-1);
    for (int i = 0; i < ind; i++)
      Inverse_EI_ER(IndexEdgeInt_to_EdgeRef(i)) = i;
    
    // if we don't have all the edges with BC ref, we exit the program
    // one possible reason is that, two separated curves have the BC ref
    if (nb != ind)
      {
	cout<<"the mesh has a problem"<<endl;
	cout<<"Is there two separated objects ? Only one is allowed"<<endl;
	DISP(nb); DISP(ind); DISP(ListePoint);
	// mesh.Write("Volumetric.mesh");
	abort();
      }
    
    nb_edges_integral = nb; mesh_surf.ReallocateBoundariesRef(nb);
    // we construct now mesh associated to the boundary
    mesh_surf.ReallocateVertices(ListePoint.GetM());
    
    // vertices of integral equation mesh
    for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
      mesh_surf.Vertex(i) = mesh.Vertex(ListePoint(i));
    
    SetPoints<Dimension2> PointsElem;
    SetMatrices<Dimension2> MatricesElem;
    VectR2 s; 
    
    R2 point_loc, point_glob, normale_fj;
    Matrix2_2 mat_dfj, dfjm1;
    
    NodleInt.Reallocate(ind, 2*r+1); 
    SignVolDof_IE.Reallocate(nodl_integral); SignVolDof_IE.Fill(true);
    EdgeIE_On_Axe.Reallocate(nb_edges_integral); EdgeIE_On_Axe.Fill(false);
    
    this->lob_IE.ConstructQuadrature(r, this->lob_IE.QUADRATURE_LOBATTO);
    // Now, we want to compute r,z, ds, tx and tz at each nodal point
    int Ntot = (r+1)*nb_edges_integral;
    Ri_lob_IE.Reallocate(Ntot); Ri_lob_IE.Fill(0);
    Zi_lob_IE.Reallocate(Ntot); Zi_lob_IE.Fill(0);
    Ds_lob_IE.Reallocate(Ntot); Ds_lob_IE.Fill(0);
    Tx_lob_IE.Reallocate(Ntot); Tx_lob_IE.Fill(0);
    Tz_lob_IE.Reallocate(Ntot); Tz_lob_IE.Fill(0);
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    if (this->condition_integral_equation==var_harmonic.LINE_HIGH_CONDUCTIVITY)
      {
	K1_lob_IE.Reallocate(Ntot);
	K2_lob_IE.Reallocate(Ntot);
      }
#endif    
    R2 pt;
    // loop on edges of integral equation
    for (int i = 0; i < ind; i++)
      {
	int num_edge_ref = IndexEdgeInt_to_EdgeRef(i);
        int iquad = mesh.BoundaryRef(num_edge_ref).numElement(0);
        int num_loc = mesh.Element(iquad).GetPositionBoundary(num_edge_ref);
	
        int hyb = mesh.Element(iquad).GetHybridType();
	
	// edge of the mesh
	int n1 = mesh.Element(iquad).numVertex(Mesh<Dimension2>::FirstExtremityEdge(hyb, num_loc));
	int n2 = mesh.Element(iquad).numVertex(Mesh<Dimension2>::SecondExtremityEdge(hyb, num_loc));
	mesh_surf.BoundaryRef(i).Init(Inverse_ListePoint(n1), Inverse_ListePoint(n2),
				      mesh.BoundaryRef(num_edge_ref).GetReference() );
	
	// we test if edge is near from axis
	if (var_harmonic.Vertex_On_Axe(n1)||var_harmonic.Vertex_On_Axe(n2))
	  EdgeIE_On_Axe(i) = true;
	
	// vertices of the element
	mesh.GetVerticesElement(iquad, s);
	
	// computation of F_i(\xi_j) for all points, DF_i, surfacic elements ...
	Fb.FjElem(s, PointsElem, mesh, iquad);
	Fb.DFjElem(s, PointsElem, MatricesElem, mesh, iquad);
	Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, mesh, iquad, num_loc); 
	
	// updating NodleInt for scalar dofs (E_\theta for FEM, Jt for integral equation)
	for (int j = 0; j <= r; j++)
	  {
	    num_dof_loc = Fb.GetLocalNumberScal(num_loc,j);
	    num_dof = mesh.Element(iquad).GetNumberDof(num_dof_loc);
	    NodleInt(i, j) = Inverse_DI_DV(num_dof);
            abort();
	    //if (var_harmonic.IsNegativeSign(iquad, num_dof_loc))
            //SignVolDof_IE(NodleInt(i,j)) = false;
	  }
	
	// updating NodleInt for vectorial dofs (E \times n for FEM, Jb for IE)
	for (int j = 0; j < Fb.GetNbVectorialDofOnEdge(); j++)
	  {
	    num_dof_loc = Fb.GetLocalNumberVec(num_loc,j);
	    num_dof = mesh.Element(iquad).GetNumberDof(num_dof_loc);
	    NodleInt(i, r+1+j) = Inverse_DI_DV(num_dof);
	    abort();
            //if (var_harmonic.IsNegativeSign(iquad, num_dof_loc))
            //SignVolDof_IE(NodleInt(i,r+1+j)) = false;
	  }
	
	// loop on Lobatto points
	for (int k = 0; k < this->lob_IE.GetNbPointsQuad(); k++)
	  {
	    // we compute point_glob = F_i(\xi_k) where \xi_k is the Lobatto point on the edge
	    Fb.GetLocalCoordOnBoundary(num_loc, lob_IE.Points(k), point_loc);
	    Fb.Fj( s, PointsElem, point_loc, point_glob, mesh, iquad);
	    DISP(s); DISP(PointsElem); DISP(point_loc); DISP(point_glob);
	    // we get too DF_i
	    Fb.DFj(s, PointsElem, point_loc, mat_dfj, mesh, iquad);
	    GetInverse(mat_dfj, dfjm1);
	    // we get normale and surfacic element dsj
	    Real_wp dsj;
	    Fb.GetNormale(dfjm1, normale_fj, dsj, num_loc);
	    // DISP(k); DISP(point_loc); DISP(point_glob);
	    // DISP(mat_dfj); DISP(dfjm1); DISP(normale_fj); DISP(dsj); DISP(num_loc);
	    
	    
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
	    if (this->condition_integral_equation==var_harmonic.LINE_HIGH_CONDUCTIVITY)
	      {
		K1_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = MatricesElem.GetK1QuadratureBoundary(k);
		K2_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = MatricesElem.GetK2QuadratureBoundary(k);
	      }
#endif
            
	    Ri_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = point_glob(0);
	    Zi_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = point_glob(1);
	    Ds_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = dsj;
	    Tx_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = -normale_fj(1);
	    Tz_lob_IE(i*lob_IE.GetNbPointsQuad() + k) = normale_fj(0);
	    
	  }
      }
    DISP(Ri_lob_IE); DISP(Zi_lob_IE); DISP(Tx_lob_IE); DISP(Tz_lob_IE);
    DISP(nodl_integral); DISP(nb_edges_integral); DISP(NodleInt);
    DISP(Ds_lob_IE);
    // DISP(mesh_surf.NodleInt);
    
    int nb_dofT = 0, nb_dofB = 0;
    for (int i = 0; i < nb_edges_integral; i++)
      for (int j = 0; j <= r; j++)
	nb_dofT = max(NodleInt(i,j), nb_dofT);
    
    nb_dofT++;
    nb_dofB = nodl_integral - nb_dofT;
    
    nodl_Jt = nb_dofT;
    nodl_Jb = nb_dofB;
    // estimation of the mean edge
    space_step_IE = 0.0;
    for (int num = 0; num < nb_edges_integral; num++)
      {
	int i = IndexEdgeInt_to_EdgeRef(num);
	// two vertices of the edge
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	space_step_IE += mesh.Vertex(n1).Distance(mesh.Vertex(n2));
      }
    space_step_IE /= nb_edges_integral;
    // DISP(mesh_surf);
    
    VectReal_wp RiDof(nb_dofT), ZiDof(nb_dofT);
    Real_wp ri, zi, tx, tz, dsi, xi;
    for (int i = 0; i < nb_edges_integral; i++)
      {
	int offset = i*this->lob_IE.GetNbPointsQuad();
	for (int j = 0; j <= Fb.GetOrder(); j++)
	  {
	    xi = Fb.PointsNodal1D(j);
	    this->ComputeInterpolatedVal(ri, zi, tx, tz, dsi, this->lob_IE, Ri_lob_IE, 
					 Zi_lob_IE, Tx_lob_IE, Tz_lob_IE, Ds_lob_IE, offset, xi);
	    
	    int num_dof_Jt = this->GetNumberDofT(i,j);
	    RiDof(num_dof_Jt) = ri;
	    ZiDof(num_dof_Jt) = zi;
	  }
      }
    
    RiDof.WriteText("RiDof.dat");
    ZiDof.WriteText("ZiDof.dat");
  }
  
  
  //! method used to compute interpolation points in order to display currents (2-D)
  /*!
    nbPoints : number of interpolation points, which need to be computed
    EdgeInterp : edge number where the point is
    CoorInterp : local coordinate of this point on the edge, between [0,1]
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>
  ::ComputeSectionGrid_IE(SurfacicMesh<Dimension2>& mesh_surf )
  {
    if (nb_edges_integral <= 0)
      {
	nbPoints_grid_boundary_ie = 0;
	return;
      }
    // first step, we compute length of each edge on the boundary
    VectReal_wp LengthEdges_Cumul(nb_edges_integral);
    
    Real_wp dsj, cumul_ds = 0.0, ri, zi;
    
    // DISP(mesh_surf.DsDof_IE);
    // computation of lineic element
    for (int i = 0; i < nb_edges_integral; i++)
      {
	for (int k = 0; k < lob_IE.GetNbPointsQuad(); k++)
	  {
	    // dsj = DsDof_IE(mesh_surf.GetNumberDofT(i,k) );
	    dsj = Ds_lob_IE(lob_IE.GetNbPointsQuad()*i + k);
	    cumul_ds += dsj*lob_IE.Weights(k);
	    DISP(dsj); DISP(lob_IE.Weights(k));
	  }
	LengthEdges_Cumul(i) = cumul_ds;
      }
    
    CoorInterp_IE.Reallocate(nbPoints_grid_boundary_ie); 
    EdgeInterp_IE.Reallocate(nbPoints_grid_boundary_ie);
    RadiusInterp_IE.Reallocate(nbPoints_grid_boundary_ie); 
    ZiInterp_IE.Reallocate(nbPoints_grid_boundary_ie);
    // step for interpolation points
    Real_wp step_h = cumul_ds/(nbPoints_grid_boundary_ie-1), si;
    // loop on interpolation points
    int num_edge = 0, nb_iter;
    
    VectReal_wp vec_dsj(lob_IE.GetNbPointsQuad()); 
    Real_wp feval, diff_feval, test, test_prec, s_rhs, coor;
    // DISP(step_h); DISP(cumul_ds); DISP(LengthEdges_Cumul);
    for (int i = 0; i < nbPoints_grid_boundary_ie; i++)
      {
	// curvilinear abscisse of interpolation point
	si = step_h*i; // DISP(si);
	// we look for the edge where it is 
	while (si > LengthEdges_Cumul(num_edge))
	  {
	    num_edge++;
	    if (num_edge >= nb_edges_integral)
	      {
		num_edge--;
		break;
	      }
	  }
	
	EdgeInterp_IE(i) = num_edge; // DISP(num_edge);
	// evaluation of CoorInterp
	Real_wp sback, snext;
	if (num_edge == 0)
	  sback = 0.0;
	else
	  sback = LengthEdges_Cumul(num_edge-1);
	
	snext = LengthEdges_Cumul(num_edge);
	// isoparametric case 
	// CoorInterp(i) = (si-sback)/(snext-sback);
	
	for (int k = 0; k < lob_IE.GetNbPointsQuad(); k++)
	  vec_dsj(k) = Ds_lob_IE(lob_IE.GetNbPointsQuad()*num_edge + k);
	
	// other cases
	// use of newton algorithm to determine s
	s_rhs = si - sback;
	coor = (si-sback)/(snext-sback); // DISP(coor); DISP(s_rhs); DISP(vec_dsj);
	feval = -s_rhs;
	for (int k = 0; k < lob_IE.GetNbPointsQuad(); k++)
	  for (int j = 0; j < lob_IE.GetNbPointsQuad(); j++)
	    feval += vec_dsj(k)*coor*lob_IE.Weights(j)*lob_IE.EvaluatePhi(k,coor*lob_IE.Points(j));
	
	test = std::abs(feval - s_rhs); test_prec = test+1.0; nb_iter = 0;
	while ((test < test_prec)&&(nb_iter<=10))
	  {
	    diff_feval = 0.0;
	    for (int k = 0; k < lob_IE.GetNbPointsQuad(); k++)
	      diff_feval += vec_dsj(k)*lob_IE.EvaluatePhi(k, coor);
	    
	    coor -= feval/diff_feval;
	    
	    feval = -s_rhs;
	    for (int k = 0; k < lob_IE.GetNbPointsQuad(); k++)
	      for (int j = 0; j < lob_IE.GetNbPointsQuad(); j++)
		feval += vec_dsj(k)*coor*lob_IE.Weights(j)*lob_IE.EvaluatePhi(k,coor*lob_IE.Points(j));
	    
	    test_prec = test;
	    test = std::abs(feval - s_rhs); nb_iter++;
	  }
	// DISP(coor);
	if (boundary_condition == var_harmonic.LINE_ABSORBING)
	  CoorInterp_IE(i) = coor;
        else
	  CoorInterp_IE(i) = 1.0-coor;
	
	ri = 0.0; zi = 0.0;
	for (int k = 0; k <= lob_IE.GetOrder(); k++)
	  ri += Ri_lob_IE(lob_IE.GetNbPointsQuad()*num_edge + k)*lob_IE.EvaluatePhi(k, CoorInterp_IE(i));
	
	for (int k = 0; k <= lob_IE.GetOrder(); k++)
	  zi += Zi_lob_IE(lob_IE.GetNbPointsQuad()*num_edge + k)*lob_IE.EvaluatePhi(k, CoorInterp_IE(i));
	
	RadiusInterp_IE(i) = ri; ZiInterp_IE(i) = zi;
	DISP(i); DISP(ri); DISP(zi); DISP(coor); DISP(si); DISP(LengthEdges_Cumul);
      }
  }
  
  
  //! creates surfacic 3-D mesh of the object, by rotation of the boundary mesh
  /*!
    \param[in] ref reference number
    \param[in,out] mesh_surf 2-D surfacic mesh
    \param[out] mesh_refined 2-D subdivided surfacic mesh
   */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement, TypeEquation>
  ::ConstructSurfacicMesh_IE(int ref, SurfacicMesh<Dimension2>& mesh_surf,
                             Mesh<Dimension2>& mesh_refined)
  {
    const Mesh<Dimension2>& mesh = var_harmonic.mesh;
    const TypeElement& Fb = var_harmonic.GetFaceBasis2(mesh.GetOrder());
    // mesh.Write("volumic.mesh");
    // first step, we subdivide boundary mesh according to the order of approximation
    // creation of mesh_refined, subdivision of mesh_ie
    int nb_subdiv = mesh.GetOrder();
    // DISP(mesh_surf); DISP(mesh_surf.edges_ref);
    if (nb_subdiv > 1)
      {
	int nb_vertices = mesh_surf.GetNbVertices() + (nb_subdiv-1)*mesh_surf.GetNbBoundaryRef();
	mesh_refined.ReallocateVertices(nb_vertices);
	int nb_edges = mesh_surf.GetNbBoundaryRef()*nb_subdiv; // DISP(nb_subdiv); DISP(nb_edges);
	mesh_refined.ReallocateBoundariesRef(nb_edges);
	
	EdgeInterp_IE.Reallocate(nb_vertices);
	CoorInterp_IE.Reallocate(nb_vertices);
	
	// we copy old vertices
	nb_vertices = 0;
	for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
	  mesh_refined.Vertex(nb_vertices++) = mesh_surf.Vertex(i);
	
	SetPoints<Dimension2> PointsElem;
	SetMatrices<Dimension2> MatricesElem;
	VectR2 s; 
	
	R2 point_loc, point_glob;
	VectReal_wp pts_subdiv(nb_subdiv+1); IVect num(nb_subdiv+1);
	for (int i = 0; i <= nb_subdiv; i++)
	  pts_subdiv(i) = Real_wp(i)/nb_subdiv;
	
	nb_edges = 0;
	// new vertices and new edges
	for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
	  {
	    num(0) = mesh_surf.BoundaryRef(i).numVertex(0);
	    num(nb_subdiv) = mesh_surf.BoundaryRef(i).numVertex(1);
	    int num_edge_ref = IndexEdgeInt_to_EdgeRef(i);
            int num_elem = mesh.BoundaryRef(num_edge_ref).numElement(0);
            int num_loc = mesh.Element(num_elem).GetPositionBoundary(num_edge_ref);
	    
	    mesh.GetVerticesElement(num_elem, s);
	    Fb.FjElem(s, PointsElem, mesh, num_elem);
	    
	    for (int j = 1; j < nb_subdiv; j++)
	      {
		num(j) = nb_vertices;
		Fb.GetLocalCoordOnBoundary(num_loc, pts_subdiv(j), point_loc);
		Fb.Fj(s, PointsElem, point_loc, point_glob, mesh, num_elem);
		mesh_refined.Vertex(nb_vertices) = point_glob; 

		nb_vertices++;
	      }
	    
	    for (int j = 0; j < nb_subdiv; j++)
	      mesh_refined.BoundaryRef(nb_edges++).Init(num(j), num(j+1), ref);
	    
	    for (int j = 0; j <= nb_subdiv; j++)
	      {
		EdgeInterp_IE(num(j)) = i;
		CoorInterp_IE(num(j)) = pts_subdiv(j);
	      }
	  }
      }
    else
      mesh_refined = mesh_surf;
    
    // DISP(mesh_refined); DISP(EdgeInterp_IE); DISP(CoorInterp_IE);
    mesh_refined.Write("refine.mesh");
    // now we call method to generate 3-D mesh by rotation of the 2-D boundary mesh
    // mesh_boundary_inside is the obtained mesh
    IVect Index_FaceSurf_to_EdgeRef;
    IVect ref_cond(mesh_refined.GetNbReferences()); ref_cond.Fill();
    mesh_refined.GenerateSurfaceOfRevolution(ref, mesh_boundary_inside, Index_FaceSurf_to_EdgeRef, 
					     Index_VertexSurf_to_Vertex, AngleVertex, ref_cond);
    
  }
  
  
  //! writes currents J in files output_Jt, output_Jb
  /*!
    \param[in] U0 solution vector of the integral equation
    \param[in] output_Jt file name where the t-component of J will be stored
    \param[in] output_Jb file name where the b-component of J will be stored
    \param[in] output_Jnorme file name where the norm of J will be stored
    \param[in] IsUnknownJ if true J is stored, otherwise K is stored
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  WriteOutputFile_CurrentsJ(const VectComplex_wp& U0, const string& output_Jt,
			    const string& output_Jb, const string& output_Jnorme, bool IsUnknownJ)
  {
    if (type_output_file_ie/2 == 0)
      WriteOutputFile_Currents2D(U0, this->mesh_ie, output_Jt, output_Jb, IsUnknownJ);
    else
      WriteOutputFile_Currents3D(U0, mesh_refined_ie, mesh_boundary_inside,
                                 output_Jt, output_Jb, output_Jnorme, IsUnknownJ);
  }
  
  
  //! Add interpolation of Jt and Jb in vectors output_Jt and output_Jb
  /*!
    \param[in] U0 solution vector
    \param[in,out] output_Jt vector containing t-component of J 
    \param[in,out] output_Jb vector containing b-component of J 
    \param[in] mesh_surf surfacic mesh
    \param[in] IsUnknown_J true if the unknown is J, false if K
   */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  AddInterpolated_Currents2D(const VectComplex_wp& U0,
                             VectComplex_wp& output_Jt, VectComplex_wp& output_Jb,
			     const SurfacicMesh<Dimension2>& mesh_surf, bool IsUnknown_J) const
  {
    //int r = var_harmonic.mesh.GetOrder();
    int N = EdgeInterp_IE.GetM(); 
    //const TypeElement& Fb = var_harmonic.GetFaceBasis2(r);    
    Complex_wp Jt, Jb;
    
    for (int i = 0; i < N; i++)
      {
	//int num_edge = EdgeInterp_IE(i); // DISP(num_edge); DISP(mesh_ie.CoorInterp_IE(i));
	
	Jb = 0.0; Jt = 0.0;
        /*
	if (IsUnknown_J)
	  for (int j = 0; j < nb_dof_Jt; j++)
	    {
	      int num_dof_Jb = nodl_Jt + GetNumberDofT(num_edge, j);
	      Jb += Fb.GetValue_ScalBoundary(j, CoorInterp_IE(i))*U0(num_dof_Jb);
	    }
	else
	  for (int j = 0; j < nb_dof_Jb; j++)
	    {
	      int num_dof_Jb = GetNumberDofB(num_edge, j);
	      Jb += Fb.GetValue_VecBoundary(j, CoorInterp_IE(i))*U0(num_dof_Jb);
	    }
	
	for (int j = 0; j < nb_dof_Jt; j++)
	  {
	    int num_dof_Jt = GetNumberDofT(num_edge,j);
	    Jt += Fb.GetValue_ScalBoundary(j, CoorInterp_IE(i))*U0(num_dof_Jt);
	  }
        */
        abort();
	output_Jb(i) += Jb; output_Jt(i) += Jt;
      }
  }
  
  
  //! writes 2-D currents Jt and Jb on files name_Jt, name_Jb
  /*!
    \param[in] U0 solution vector of the integral equation
    \param[in] mesh_surf surfacic mesh
    \param[in] name_Jt file name where the t-component of J will be stored
    \param[in] name_Jb file name where the b-component of J will be stored
    \param[in] UnknownJ if true J is stored, otherwise K is stored
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement, TypeEquation>::
  WriteOutputFile_Currents2D(const VectComplex_wp& U0, const SurfacicMesh<Dimension2>& mesh_surf,
			     const string& name_Jt, const string& name_Jb, bool UnknownJ)
  {
    VectComplex_wp output_Jt, output_Jb;
    if (TypeEquation::number_mode == var_harmonic.number_minimum_mode)
      {
	output_Jt.Reallocate(EdgeInterp_IE.GetM()); output_Jt.Fill(0);
	output_Jb.Reallocate(EdgeInterp_IE.GetM()); output_Jb.Fill(0);
      }
    else
      {
	if (type_output_file_ie == 0)
	  {
	    output_Jb.Read(name_Jb);
	    output_Jt.Read(name_Jt);
	  }
	else
	  {
	    ifstream file_Jt(name_Jt.data()), file_Jb(name_Jb.data());
	    output_Jt.Reallocate(EdgeInterp_IE.GetM());
	    output_Jb.Reallocate(EdgeInterp_IE.GetM());
	    Real_wp tmp, val_real, val_imag;
	    for (int i = 0; i < nbPoints_grid_boundary_ie; i++)
	      {
		file_Jt>>tmp>>tmp>>val_real>>val_imag;
		output_Jt(i) = Complex_wp(val_real, val_imag);
		file_Jb>>tmp>>tmp>>val_real>>val_imag;
		output_Jb(i) = Complex_wp(val_real, val_imag);
	      }
	    file_Jt.close(); file_Jb.close();
	  }
      }
    
    AddInterpolated_Currents2D(U0, output_Jt, output_Jb, mesh_surf, UnknownJ);
    
    if (type_output_file_ie == 0)
      {
	output_Jb.Write(name_Jb);
	output_Jt.Write(name_Jt);
      }
    else
      {
	ofstream file_Jt(name_Jt.data()), file_Jb(name_Jb.data());
	// DISP(mesh_surf.RadiusInterp_IE); DISP(mesh_surf.ZiInterp_IE);
	// DISP(output_Jb); DISP(output_Jt);
	for (int i = 0; i < nbPoints_grid_boundary_ie; i++)
	  {
	    file_Jt<<RadiusInterp_IE(i)<<"  "<<ZiInterp_IE(i)<<"  "<<real(output_Jt(i))
                   <<"  "<<imag(output_Jt(i))<<endl;
	    
            file_Jb<<RadiusInterp_IE(i)<<"  "<<ZiInterp_IE(i)<<"  "<<real(output_Jb(i))
                   <<"  "<<imag(output_Jb(i))<<endl;
	  }
	file_Jt.close(); file_Jb.close();
      }
  }
  
  
  //! writes 3-D components of J in files name_Jt, name_Jb, name_Jnorme
  /*!
    \param[in] U0 solution vector of the integral equation
    \param[in] mesh2d surfacic 2-D mesh
    \param[in] mesh3d surfacic 3-D mesh
    \param[in] name_Jt file name where the t-component of J will be stored
    \param[in] name_Jb file name where the b-component of J will be stored
    \param[in] name_Jnorme file name where the norm of J will be stored
    \param[in] IsUnknown_J if true J is stored, otherwise K is stored
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  WriteOutputFile_Currents3D(const VectComplex_wp& U0, const Mesh<Dimension2>& mesh2d,
                             const Mesh<Dimension3>& mesh3d, 
			     const string& name_Jt, const string& name_Jb,
                             const string& name_Jnorme, bool IsUnknown_J)
  {
    //int r = var_harmonic.mesh.GetOrder();
    //const TypeElement& Fb = var_harmonic.GetFaceBasis2(r);    
    int N = mesh3d.GetNbVertices();
    VectComplex_wp output_Jt, output_Jb; VectReal_wp output_Jnorme;
    if (TypeEquation::number_mode == var_harmonic.number_minimum_mode)
      {
	output_Jt.Reallocate(N); output_Jt.Fill(0);
	output_Jb.Reallocate(N); output_Jb.Fill(0);
      }
    else
      {
	int success_read = ReadMedit(output_Jb, name_Jb);
	if (success_read != 0)
	  {
	    cout<<"we cannot write 3D currents"<<endl;
	    return;
	  }
	success_read = ReadMedit(output_Jt, name_Jt);
	if (success_read != 0)
	  {
	    cout<<"we cannot write 3D currents"<<endl;
	    return;
	  }
      }
    
    output_Jnorme.Reallocate(output_Jt.GetM());
    
    Complex_wp Jt, Jb, phase;
    Real_wp teta, x;
    for (int i = 0; i < N; i++)
      {
	// theta of the 3-D point
	teta = AngleVertex(i);
	int num_point = Index_VertexSurf_to_Vertex(i);
	//int num_edge = EdgeInterp_IE(num_point);
	x = CoorInterp_IE(num_point);
	
	// computation of Jt and Jb
	Jb = 0.0; Jt = 0.0;
        /*
	if (!IsUnknown_J)
	  for (int j = 0; j < nb_dof_Jb; j++)
	    {
	      int num_dof_Jb = GetNumberDofB(num_edge,j);
	      Jb += Fb.GetValue_VecBoundary(j, x)*U0(num_dof_Jb);
	    }
	else
	  for (int j = 0; j < nb_dof_Jt; j++)
	    {
	      int num_dof_Jb = nodl_Jt + GetNumberDofT(num_edge, j);
	      Jb += Fb.GetValue_ScalBoundary(j, x)*U0(num_dof_Jb);
	    }
	
	for (int j=0 ; j<nb_dof_Jt ; j++)
	  {
	    int num_dof_Jt = GetNumberDofT(num_edge, j);
	    Jt += Fb.GetValue_ScalBoundary(j, x)*U0(num_dof_Jt);
	  }
        */
        abort();
	phase = exp(-Iwp*Complex_wp(TypeEquation::number_mode)*teta);
	output_Jb(i) += Jb*phase; output_Jt(i) += Jt*phase;
	// output_Jnorme = sqrt(real(output_Jb(i)*conj(output_Jb(i))
        // + output_Jt(i)*conj(output_Jt(i))));
	output_Jnorme(i) = sqrt(real(output_Jb(i))*real(output_Jb(i)) 
                                + real(output_Jt(i))*real(output_Jt(i)));
      }
    
    // DISP(N)
    Dimension3 dim;
    WriteMedit(output_Jb, name_Jb, dim);
    WriteMedit(output_Jt, name_Jt, dim);
    WriteMedit(output_Jnorme, name_Jnorme, dim);
  }
  
  
  //! computation of ri, zi, txi, tzi and dsi
  /*!
    lob Gauss-Lobatto points used to interpolate function
    Ri_lob Ri_lob(offset1 : (offset1+N-1) ) evaluation of the function ri on Gauss-Lobatto points
    offset1 offset used in arrays Ri_lob, Zi_lob, Tx_lob, Tz_lob and Ds_lob
    x point where ri, zi, txi, tzi and dsi are evaluated
  */
  template<class TypeElement, class TypeEquation> inline
  void VarIntegralEquationAxi<TypeElement, TypeEquation>::
  ComputeInterpolatedVal(Real_wp& ri, Real_wp& zi, Real_wp& txi, Real_wp& tzi,
			 Real_wp& dsi, const Globatto<Real_wp>& lob,
			 const VectReal_wp& Ri_lob, const VectReal_wp& Zi_lob,
                         const VectReal_wp& Tx_lob, 
			 const VectReal_wp& Tz_lob, const VectReal_wp& Ds_lob,
                         int offset1, const Real_wp& x) const
  {
    // we compute ri, zi, txi, tzi and dsi, by interpolation with gauss-lobatto points
    ri = 0.0; zi = 0.0; txi = 0.0; tzi = 0.0; dsi = 0.0;
    Real_wp coef;
    for (int j = 0; j < lob.GetNbPointsQuad(); j++)
      {
	coef = lob.EvaluatePhi(j, x);
	// DISP(j); DISP(coef); DISP(Ri_lob(offset1+j));
	ri += coef*Ri_lob(offset1+j);
	zi += coef*Zi_lob(offset1+j);
	txi += coef*Tx_lob(offset1+j);
	tzi += coef*Tz_lob(offset1+j);
	dsi += coef*Ds_lob(offset1+j);
      }
  }
  
  
  //! computation of some intermediary variables for integral equation
  /*!
    G1 = int G(x,x') cos(m phi) dphi
    Gcos = int G(x,x') cos(phi) cos(m phi) dphi
    Gsin = I*int G(x,x') sin(phi) sin(m phi) dphi
    G(x,x') = exp(i k |x-x'|) / (4 pi |x-x'|)
    R = |x-x'| = sqrt{ dist + 2*prodr (1.0 - cos(phi) ) }
    where dist = (r'-r)^2 + (z'-z)^2 and prodr = r*r'
    integration from 0 to pi if first_subdiv is equal to 0 (for regular integrals)
    for singular integrals, we evaluate regular part,
    by integration from phi1 to pi (first_subdiv is equal to 1)
  */
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement, TypeEquation>::
  ComputeG1_Gcos_Gsin_EFIE(Globatto<Real_wp>& gauss,
			   const Real_wp& phi0, const Real_wp& phi1,
			   Real_wp& dist, Real_wp& prodr, Complex_wp& G1,
			   Complex_wp& Gcos, Complex_wp& Gsin) const
  {
    Real_wp phi_init, phi_end, cos_phi, sin_phi, phi, cos_mphi, sin_mphi, R, wn;
    Complex_wp G; G1 = 0.0; Gcos = 0.0; Gsin = 0.0;
    
    
    phi_init = phi0; phi_end = phi1;
    Complex_wp iw = Iwp*var_harmonic.GetOmega();
    for (int n = 0; n < gauss.GetNbPointsQuad(); n++)
      {
	phi = phi_init*(1.0-gauss.Points(n)) + phi_end*gauss.Points(n);
	cos_phi = cos(phi); sin_phi = sin(phi);
	cos_mphi = cos(TypeEquation::number_mode*phi);
        sin_mphi = sin(TypeEquation::number_mode*phi);
	// dist = \sqrt{ (r'- r)^2 + (z'-z)^2 + 2*r*r'(1.0-cos \phi) } 
	R = sqrt(dist + 2.0*prodr*(1.0-cos_phi));
	// green function
	G = exp(iw*R)/R;
	// weight of integration
	wn = gauss.Weights(n)*(phi_end-phi_init); 
	G1 += G*cos_mphi*wn;
	Gcos += G*cos_phi*cos_mphi*wn;
	Gsin += G*sin_phi*sin_mphi*wn;
      }
    
    G1 = -Iwp*G1; Gcos = -Iwp*Gcos;
  }
  
  
  //! computation of some intermediary variables for integral equation
  /*!
    same thing than for EFIE but with MFIE
  */
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeG1_Gcos_Gsin_MFIE(Globatto<Real_wp>& gauss,
			   const Real_wp& phi0, const Real_wp& phi1,
			   Real_wp& dist, Real_wp& prodr, Complex_wp& G1,
			   Complex_wp& Gcos, Complex_wp& Gsin) const
  {
    Real_wp phi_init, phi_end, cos_phi, sin_phi, phi, cos_mphi, sin_mphi, R, wn;
    Complex_wp G; G1 = 0.0; Gcos = 0.0; Gsin = 0.0;
    
    phi_init = phi0; phi_end = phi1;
    Complex_wp iw = Iwp*var_harmonic.GetOmega();
    for (int n = 0; n < gauss.GetNbPointsQuad(); n++)
      {
	phi = phi_init*(1.0-gauss.Points(n)) + phi_end*gauss.Points(n);
	cos_phi = cos(phi); sin_phi = sin(phi);
	cos_mphi = cos(TypeEquation::number_mode*phi);
        sin_mphi = sin(TypeEquation::number_mode*phi);
	// dist = \sqrt{ (r'- r)^2 + (z'-z)^2 + 2*r*r'(1.0-cos \phi) } 
	R = sqrt(dist + 2.0*prodr*(1.0-cos_phi));
	// gradient of green function
	G = exp(iw*R)/R * (iw/R - 1.0/(R*R));
	// weight of integration
	wn = gauss.Weights(n)*(phi_end-phi_init); 
	G1 += G*cos_mphi*wn;
	Gcos += G*cos_phi*cos_mphi*wn;
	Gsin += G*sin_phi*sin_mphi*wn;
      }
    phi_init = phi_end;
    Gsin = Iwp*Gsin;
  }
  
  
  //!  computation of some intermediary variables for integral equation
  /*!
    both EFIE and MFIE
  */
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeG1_Gcos_Gsin_CFIE(const Globatto<Real_wp>& gauss,
			   const Real_wp& phi0, const Real_wp& phi1, Real_wp& dist, Real_wp& prodr, 
			   Complex_wp& G1_efie,	Complex_wp& Gcos_efie, Complex_wp& Gsin_efie,
			   Complex_wp& G1_mfie,	Complex_wp& Gcos_mfie, Complex_wp& Gsin_mfie) const
  {
    Real_wp phi_init, phi_end, cos_phi, sin_phi, phi, cos_mphi, sin_mphi, R, wn;
    Complex_wp G,Gnabla; 
    G1_efie = 0.0; Gcos_efie = 0.0; Gsin_efie = 0.0;
    G1_mfie = 0.0; Gcos_mfie = 0.0; Gsin_mfie = 0.0;
    
    phi_init = phi0; phi_end = phi1;
    Complex_wp iw = Iwp*var_harmonic.GetOmega();
    for (int n = 0; n < gauss.GetNbPointsQuad(); n++)
      {
	phi = phi_init*(1.0-gauss.Points(n)) + phi_end*gauss.Points(n);
	cos_phi = cos(phi); sin_phi = sin(phi);
	cos_mphi = cos(TypeEquation::number_mode*phi);
        sin_mphi = sin(TypeEquation::number_mode*phi);
	// dist = \sqrt{ (r'- r)^2 + (z'-z)^2 + 2*r*r'(1.0-cos \phi) } 
	R = sqrt(dist + 2.0*prodr*(1.0-cos_phi));
	// green function and gradient
	G = exp(iw*R)/R;
	Gnabla = G * (iw/R - 1.0/(R*R));
	// weight of integration
	wn = gauss.Weights(n)*(phi_end-phi_init); 
	G1_efie += G*cos_mphi*wn;
	Gcos_efie += G*cos_phi*cos_mphi*wn;
	Gsin_efie += G*sin_phi*sin_mphi*wn;
	G1_mfie += Gnabla*cos_mphi*wn;
	Gcos_mfie += Gnabla*cos_phi*cos_mphi*wn;
	Gsin_mfie += Gnabla*sin_phi*sin_mphi*wn;
      }
    phi_init = phi_end;
    // }
    G1_efie = -Iwp*G1_efie; Gcos_efie = -Iwp*Gcos_efie;
    Gsin_mfie = Iwp*Gsin_mfie;
  }
  
  
  //! intermediary coefficients for EFIE
  /*!
    A1 = (omega r r' (t_x t_x' Gcos + t_z t_z' G1) - t_x t_x' G1 / omega) ds ds'
    A2 = - G1 t_x r' ds / omega
    A3 = - G1 t_x' r ds' / omega
    A4 = - G1 r r' / omega
  
    A5 = (omega Gcos r r' - m^2 / omega ) ds ds'
  
    A6 = (omega Gsin t_x r r' + i m G1 t_x / omega ) ds ds'
    A7 = i m G1 r ds' / omega

    A8 = -(omega Gsin t_x' r r' + i m G1 t_x' / omega ) ds ds'
    A9 = -i m G1 r' ds / omega
  */
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeCoefficients_EFIE(Real_wp& ri, Real_wp& zi, Real_wp& txi, Real_wp& tzi, Real_wp& dsi, 
			   Real_wp& rj, Real_wp& zj, Real_wp& txj, Real_wp& tzj, Real_wp& dsj,
			   Real_wp& prodr, Real_wp& prod_ds, Complex_wp& im, Complex_wp& m2, 
			   Complex_wp& G1, Complex_wp& Gcos, Complex_wp& Gsin,
                           Complex_wp& A1, Complex_wp& A2,
			   Complex_wp& A3, Complex_wp& A4,
                           Complex_wp& A5, Complex_wp& A6,Complex_wp& A7, 
			   Complex_wp& A8, Complex_wp& A9) const
  {
    Real_wp omega = var_harmonic.GetOmega();
    A1 = (omega*prodr*( txi*txj*Gcos + tzi*tzj*G1) - txi*txj*G1/omega)*prod_ds;
    A2 = -G1*txi*rj*dsi/omega;
    A3 = -G1*txj*ri*dsj/omega;
    A4 = -G1*prodr/omega;
    
    A5 = (omega*Gcos*prodr - m2/omega*G1)*prod_ds;
    
    A6 = (omega*Gsin*txi*prodr + im*G1*txi/omega)*prod_ds;
    A7 = im*G1*ri*dsj/omega;
    
    A8 = -(omega*Gsin*txj*prodr + im*G1*txj/omega)*prod_ds;
    A9 = -im*G1*rj*dsi/omega;
  }
  
  
  //! intermediary coefficients for MFIE
  /*!
    F11 = Gsin (r' t_z' t_x - r t_x' t_z - (z'-z) t_x' t_x)
    F12 = -G1 r' t_z + Gcos (r t_z + (z'-z) t_x )
    F21 = -G1 r t_z' + Gcos (r' t_z' - (z'-z) t_x')
    F22 = -Gsin (z'-z)
  */
  template<class TypeElement,class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeCoefficients_MFIE(Real_wp& ri, Real_wp& zi, Real_wp& txi, Real_wp& tzi, Real_wp& dsi, 
			   Real_wp& rj, Real_wp& zj, Real_wp& txj, Real_wp& tzj, Real_wp& dsj,
			   Real_wp& diffz, Real_wp& prodr, Complex_wp& G1,
                           Complex_wp& Gcos, Complex_wp& Gsin, 
			   Complex_wp& F11, Complex_wp& F12, Complex_wp& F21, Complex_wp& F22) const
  {
    F11 = -Gsin*(rj*tzj*txi - ri*txj*tzi - diffz*txj*txi)*prodr;
    F12 = G1*rj*tzi - Gcos*(ri*tzi+diffz*txi); F12 *= prodr;
    F21 = G1*ri*tzj - Gcos*(rj*tzj-diffz*txj); F21 *= prodr;
    F22 = Gsin*diffz*prodr;
  }
  
  
  //! add interactions to mat_tt, mat_bb, mat_tb, mat_bt (EFIE)
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  UpdateEFIEMatrices(int i, int j, const Complex_wp& A1, const Complex_wp& A2,
                     const Complex_wp& A3, const Complex_wp& A4,
		     const Complex_wp& A5, const Complex_wp& A6,
                     const Complex_wp& A7, const Complex_wp& A8, const Complex_wp& A9,
		     const Matrix<Real_wp>& Val_Phi1, const Matrix<Real_wp>& Grad_Phi1, 
                     const Matrix<Real_wp>& Val_Psi1,
		     const Matrix<Real_wp>& Val_Phi2, const Matrix<Real_wp>& Grad_Phi2,
                     const Matrix<Real_wp>& Val_Psi2,
		     Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
		     Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const
  {
    // interactions Jt <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_tt(k,l) += Val_Phi1(k,i)*Val_Phi2(l,j)*A1;
	  mat_tt(k,l) += Val_Phi1(k,i)*Grad_Phi2(l,j)*A2;
	  mat_tt(k,l) += Grad_Phi1(k,i)*Val_Phi2(l,j)*A3;
	  mat_tt(k,l) += Grad_Phi1(k,i)*Grad_Phi2(l,j)*A4;
	}
    
    // interactions Jb <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_bb.GetN(); l++)
	{
	  mat_bb(k,l) += Val_Phi1(k,i)*Val_Psi2(l,j)*A5;
	}
    
    // interactions Jt <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_tb.GetN(); l++)
	{
	  mat_tb(k,l) += Val_Phi1(k,i)*Val_Psi2(l,j)*A6;
	  mat_tb(k,l) += Grad_Phi1(k,i)*Val_Psi2(l,j)*A7;
	}
    
    // interactions Jb <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_bt(k,l) += Val_Phi1(k,i)*Val_Phi2(l,j)*A8;
	  mat_bt(k,l) += Val_Phi1(k,i)*Grad_Phi2(l,j)*A9;
	}
  }
  
  
  //! add interactions to mat_tt, mat_bb, mat_tb, mat_bt (MFIE)
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  UpdateMFIEMatrices(int i, int j, const Complex_wp& F11, const Complex_wp& F12,
                     const Complex_wp& F21, const Complex_wp& F22,
		     const Matrix<Real_wp>& Val_Phi1, const Matrix<Real_wp>& Val_Psi1,
		     const Matrix<Real_wp>& Val_Phi2, const Matrix<Real_wp>& Val_Psi2,
		     Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
		     Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const
  {
    // interactions Jt <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_tt(k,l) += Val_Phi1(k,i)*Val_Phi2(l,j)*F21;
	}
    
    // interactions Jb <-> Jb
    for (int k = 0; k < mat_bb.GetM(); k++)
      for (int l = 0; l < mat_bb.GetM(); l++)
	{
	  mat_bb(k,l) -= Val_Phi1(k,i)*Val_Psi2(l,j)*F12;
	}
    
    // interactions Jt <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_tb.GetN(); l++)
	{
	  mat_tb(k,l) += Val_Phi1(k,i)*Val_Psi2(l,j)*F22;
	}
    
    // interactions Jb <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_bt(k,l) -= Val_Phi1(k,i)*Val_Phi2(l,j)*F11;
	}
  }
  
  
  //! add interactions to mat_tt, mat_bb, mat_tb, mat_bt (EFIE)
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  UpdateEFIEMatrices(const Complex_wp& A1, const Complex_wp& A2,
                     const Complex_wp& A3, const Complex_wp& A4,
		     const Complex_wp& A5, const Complex_wp& A6,
                     const Complex_wp& A7, const Complex_wp& A8, const Complex_wp& A9,
		     const VectReal_wp& Val_Phi1, const VectReal_wp& Grad_Phi1,
                     const VectReal_wp& Val_Psi1,
		     const VectReal_wp& Val_Phi2, const VectReal_wp& Grad_Phi2,
                     const VectReal_wp& Val_Psi2,
		     Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
		     Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const
  {
    // interactions Jt <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_tt(k,l) += Val_Phi1(k)*Val_Phi2(l)*A1;
	  mat_tt(k,l) += Val_Phi1(k)*Grad_Phi2(l)*A2;
	  mat_tt(k,l) += Grad_Phi1(k)*Val_Phi2(l)*A3;
	  mat_tt(k,l) += Grad_Phi1(k)*Grad_Phi2(l)*A4;
	}
    
    // interactions Jb <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_bb.GetN(); l++)
	{
	  mat_bb(k,l) += Val_Phi1(k)*Val_Psi2(l)*A5;
	}
    
    // interactions Jt <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_tb.GetN(); l++)
	{
	  mat_tb(k,l) += Val_Phi1(k)*Val_Psi2(l)*A6;
	  mat_tb(k,l) += Grad_Phi1(k)*Val_Psi2(l)*A7;
	}
    
    // interactions Jb <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_bt(k,l) += Val_Phi1(k)*Val_Phi2(l)*A8;
	  mat_bt(k,l) += Val_Phi1(k)*Grad_Phi2(l)*A9;
	}
  }
  
  
  //! add interactions to mat_tt, mat_bb, mat_tb, mat_bt (MFIE)
  template<class TypeElement, class TypeEquation>
  inline void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  UpdateMFIEMatrices(const Complex_wp& F11, const Complex_wp& F12,
                     const Complex_wp& F21, const Complex_wp& F22,
		     const VectReal_wp& Val_Phi1, const VectReal_wp& Val_Psi1,
		     const VectReal_wp& Val_Phi2, const VectReal_wp& Val_Psi2,
		     Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
		     Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const
  {
    // interactions Jt <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_tt(k,l) += Val_Phi1(k)*Val_Phi2(l)*F21;
	}
    
    // interactions Jb <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_bb.GetN(); l++)
	{
	  mat_bb(k,l) -= Val_Phi1(k)*Val_Psi2(l)*F12;
	}
    
    // interactions Jt <-> Jb
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < mat_bb.GetN(); l++)
	{
	  mat_tb(k,l) += Val_Phi1(k)*Val_Psi2(l)*F22;
	}
    
    // interactions Jb <-> Jt
    for (int k = 0; k < nb_dof_Jt; k++)
      for (int l = 0; l < nb_dof_Jt; l++)
	{
	  mat_bt(k,l) -= Val_Phi1(k)*Val_Phi2(l)*F11;
	}
  }
  

  //! returns number of subdivisions for integration depending of phi 
  //! (integration on interval [0 2pi])
  template<class TypeElement, class TypeEquation>
  Real_wp VarIntegralEquationAxi<TypeElement,TypeEquation>::
  EvaluatePhi1(const Real_wp& r1, const Real_wp& z1, const Real_wp& r2,
               const Real_wp& z2, const Real_wp& Rmax) const
  {
    // initial angle so that a few oscillations are permitted
    Real_wp phi1 
      = min(pi_wp/(abs(TypeEquation::number_mode)
                   + 4.0*var_harmonic.GetXmax()*var_harmonic.GetFrequency()), 0.5*pi_wp);
    return phi1;
  }

  
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  PerformRegularIntegration(const Real_wp& phi_init, Globatto<Real_wp>& gauss_X,
                            Globatto<Real_wp>& gauss_Y, Globatto<Real_wp>& gauss_phi,
                            Matrix<Real_wp>& Val_PhiX, Matrix<Real_wp>& Val_PsiX,
                            Matrix<Real_wp>& Grad_PhiX, Matrix<Real_wp>& Val_PhiY,
                            Matrix<Real_wp>& Val_PsiY, Matrix<Real_wp>& Grad_PhiY,
			    int offset1, int offset2, Matrix<Complex_wp>& mat_tt,
                            Matrix<Complex_wp>& mat_bb, Matrix<Complex_wp>& mat_tb,
                            Matrix<Complex_wp>& mat_bt,
                            Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22,
			    Matrix<Complex_wp>& coup12, Matrix<Complex_wp>& coup21)
  {
    int Nx = gauss_X.GetNbPointsQuad();
    int Ny = gauss_Y.GetNbPointsQuad();
    Real_wp xi, xj, ri, zi, rj, zj, tzi, tzj,
      dsi, dsj, prod_ds, diffr, diffz, prodr, dist, Rmin, Rmax;
    
    Real_wp wn, txi, txj;
    Complex_wp m2 = TypeEquation::number_mode*TypeEquation::number_mode;
    Complex_wp Gcos_efie, Gsin_efie, G1_efie, Gcos_mfie, Gsin_mfie, G1_mfie;
    Complex_wp A1, A2, A3, A4, A5, A6, A7, A8, A9, F11, F12, F21, F22, zero(0);
    Complex_wp im(0, TypeEquation::number_mode), zc;
    
    bool efie_to_be_computed = false, mfie_to_be_computed = false;
    if (this->alpha_efie != zero)
      efie_to_be_computed = true;
    if (this->beta_mfie != zero)
      mfie_to_be_computed = true;

    if (phi_init != pi_wp)
      for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	  {
	    // xi is equal to local abscisse of s in the edge i1, xj for the edge i2
	    xi = gauss_X.Points(i); xj = gauss_Y.Points(j);
	    // for these abscisses s and s', evaluation of r, z, tx, tz and ds
	    ComputeInterpolatedVal(ri, zi, txi, tzi, dsi, lob_IE,
				   Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE,
                                   Ds_lob_IE, offset1, xi);
	    
	    ComputeInterpolatedVal(rj, zj, txj, tzj, dsj, lob_IE,
				   Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE,
                                   Ds_lob_IE, offset2, xj);
	    
	    prod_ds = dsi*dsj;
	    diffr = rj - ri;
	    diffz = zj - zi;
	    prodr = ri*rj;
	    dist = diffr*diffr + diffz*diffz;
	    // two extrema values of R for phi going from 0 to pi
	    Rmin = sqrt(dist); Rmax = sqrt(dist+4.0*prodr);
	    
	    // integration over phi
	    // computation of G1, Gcos and Gsin
	    ComputeG1_Gcos_Gsin_CFIE(gauss_phi, phi_init, pi_wp, 
				     dist, prodr, G1_efie, Gcos_efie, Gsin_efie,
				     G1_mfie, Gcos_mfie, Gsin_mfie);
	    
	    wn = gauss_X.Weights(i)*gauss_Y.Weights(j)/(2.0*pi_wp);
	    G1_efie *= wn; Gcos_efie *= wn; Gsin_efie *= wn;
	    G1_mfie *= wn; Gcos_mfie *= wn; Gsin_mfie *= wn;
	    
	    if (efie_to_be_computed)
	      ComputeCoefficients_EFIE(ri, zi, txi, tzi, dsi, rj, zj, txj, tzj, dsj,
                                       prodr, prod_ds, 
				       im, m2, G1_efie, Gcos_efie, Gsin_efie, 
				       A1, A2, A3, A4, A5, A6, A7, A8, A9);
	    
	    if (mfie_to_be_computed)
	      ComputeCoefficients_MFIE(ri, zi, txi, tzi, dsi, rj, zj, txj, tzj, dsj, diffz, prodr, 
				       G1_mfie, Gcos_mfie, Gsin_mfie, F11, F12, F21, F22);
	    
	    if (this->alpha_efie != zero)
	      {
		zc = this->alpha_efie;
		// add interactions on local matrices
		UpdateEFIEMatrices(i, j, zc*A1, zc*A2, zc*A3, zc*A4,
                                   zc*A5, zc*A6, zc*A7, zc*A8, zc*A9,
				   Val_PhiX, Grad_PhiX, Val_PhiX, Val_PhiY, Grad_PhiY, Val_PhiY,
				   mat_tt, mat_bb, mat_tb, mat_bt);
		
		zc *= prod_ds;
		if (this->coupling_integral_equation == this->BIE_FEM)
		  UpdateMFIEMatrices(i, j, zc*F11, zc*F12, zc*F21, zc*F22, Val_PhiX, Val_PsiX,
				     Val_PhiY, Val_PsiY, coup11, coup22, coup12, coup21);
		
	      } // end if this->alpha_efie different from 0
	    
	    if (this->beta_mfie != zero)
	      {
		zc = this->beta_mfie;
		if (this->coupling_integral_equation == this->BIE_FEM)
		  UpdateEFIEMatrices(i, j, zc*A1, zc*A2, zc*A3, zc*A4,
                                     zc*A5, zc*A6, zc*A7, zc*A8, zc*A9,
				     Val_PhiX, Grad_PhiX, Val_PsiX, Val_PhiY, Grad_PhiY, Val_PsiY,
				     coup11, coup22, coup12, coup21);
		
		zc *= prod_ds;
		// add interactions on local matrices
		UpdateMFIEMatrices(i, j, zc*F11, zc*F12, zc*F21, zc*F22, Val_PhiX, Val_PhiX,
				   Val_PhiY, Val_PhiY, mat_tt, mat_bb, mat_tb, mat_bt);
		
	      }
		  
	  } // end loop on quadrature points
  }
  
  
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  PerformGaussSquaredIntegration(const Real_wp& phi1, Globatto<Real_wp>& gauss_X,
                                 Globatto<Real_wp>& gauss_Y,
				 Globatto<Real_wp>& gauss_phi, Matrix<Real_wp>& Val_PhiX,
                                 Matrix<Real_wp>& Val_PsiX, Matrix<Real_wp>& Grad_PhiX,
				 int offset1, int offset2, Matrix<Complex_wp>& mat_tt,
                                 Matrix<Complex_wp>& mat_bb,
				 Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt,
                                 Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22,
				 Matrix<Complex_wp>& coup12, Matrix<Complex_wp>& coup21)
  {
    int Nx = gauss_X.GetNbPointsQuad();
    int Ny = gauss_Y.GetNbPointsQuad();
    Real_wp xi, xj, ri, zi, rj, zj, tzi, tzj, dsi, dsj, prod_ds, diffr, diffz, prodr, dist;
    Real_wp wn, zero(0), txi, txj;
    Complex_wp m2 = square(TypeEquation::number_mode);
    Complex_wp Gcos_efie, Gsin_efie, G1_efie, Gcos_mfie, Gsin_mfie, G1_mfie;
    Complex_wp A1, A2, A3, A4, A5, A6, A7, A8, A9, F11, F12, F21, F22;
    Complex_wp im(0, TypeEquation::number_mode), zc;
    
    bool efie_to_be_computed = false, mfie_to_be_computed = false;
    if (this->alpha_efie != zero)
      efie_to_be_computed = true;
    if (this->beta_mfie != zero)
      mfie_to_be_computed = true;

    //int r = var_harmonic.mesh.GetOrder();
    //const TypeElement& Fb = var_harmonic.GetFaceBasis(r);
    for (int i = 0; i < Nx; i++)
      {
	xi = gauss_X.Points(i);
	ComputeInterpolatedVal(ri, zi, txi, tzi, dsi, lob_IE,
			       Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE, Ds_lob_IE, offset1, xi);
		      
	VectReal_wp PointsY(Ny*4);
	VectReal_wp WeightsY(Ny*4);
		      
	Real_wp xs(0);
	if (offset1 == offset2)
	  xs = xi;
	else if (offset1 == (offset2-lob_IE.GetNbPointsQuad()))
	  {
	    xs = 0.0;
	  }
	else
	  xs = 1.0;
	
	Real_wp x1 = max(0.0, 2.0*xs-1.0);
	Real_wp x2 = min(1.0, 2.0*xs);
	
	int nb = 0;
	Real_wp xsi_square;
	for (int k = 0; k < Ny; k++)
	  {
	    xj = gauss_Y.Points(k)*x1;
	    PointsY(nb) = xj;
	    WeightsY(nb++) = x1*gauss_Y.Weights(k);
	    
	    xsi_square = gauss_Y.Points(k)*gauss_Y.Points(k);
	    xj = xsi_square*x1 + (1.0-xsi_square)*xs;
	    PointsY(nb) = xj;
	    WeightsY(nb++) = 2.0*(xs-x1)*gauss_Y.Weights(k)*gauss_Y.Points(k);
	    
	    xj = xsi_square*x2 + (1.0-xsi_square)*xs;
	    PointsY(nb) = xj;
	    WeightsY(nb++) = 2.0*(x2-xs)*gauss_Y.Weights(k)*gauss_Y.Points(k);
	    
	    xj = (1.0-gauss_Y.Points(k))*x2 + gauss_Y.Points(k);
	    PointsY(nb) = xj;
	    WeightsY(nb++) = (1.0-x2)*gauss_Y.Weights(k);
	  }
	
	Matrix<Real_wp> ValPhi_Sy_square(nb_dof_Jt, 4*Ny);
	Matrix<Real_wp> GradPhi_Sy_square(nb_dof_Jt, 4*Ny);
	Matrix<Real_wp> ValPsi_Sy_square(nb_dof_Jb, 4*Ny);
        /* 
	for (int k = 0; k < nb_dof_Jt; k++)
	  for (int j = 0; j < 4*Ny; j++)
	    {
	      ValPhi_Sy_square(k, j) = Fb.GetValue_ScalBoundary(k, PointsY(j));
	      GradPhi_Sy_square(k, j) = Fb.GetGradient_ScalBoundary(k, PointsY(j));
	    }
		      
	for (int k = 0; k < nb_dof_Jb; k++)
	  for (int j = 0; j < 4*Ny; j++)
	    ValPsi_Sy_square(k, j) = Fb.GetValue_VecBoundary(k, PointsY(j));
        */
        abort();
	for (int j = 0; j < 4*Ny; j++)
	  if (WeightsY(j) != 0)
	    {
	      xj = PointsY(j);
	      // for these abscisses s and s', evaluation of r, z, tx, tz and ds
	      ComputeInterpolatedVal(rj, zj, txj, tzj, dsj, lob_IE,
				     Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE,
                                     Ds_lob_IE, offset2, xj);
	      
	      prod_ds = dsi*dsj;
	      diffr = rj - ri;
	      diffz = zj - zi;
	      prodr = ri*rj;
	      dist = diffr*diffr + diffz*diffz;
	      
	      // integration over phi
	      // computation of G1, Gcos and Gsin
	      ComputeG1_Gcos_Gsin_CFIE(gauss_phi, Real_wp(0), phi1, 
				       dist, prodr, G1_efie, Gcos_efie, Gsin_efie,
				       G1_mfie, Gcos_mfie, Gsin_mfie);
			    
	      wn = gauss_X.Weights(i)*WeightsY(j)/(2.0*pi_wp);
	      G1_efie *= wn; Gcos_efie *= wn; Gsin_efie *= wn;
	      G1_mfie *= wn; Gcos_mfie *= wn; Gsin_mfie *= wn;
			    
	      if (efie_to_be_computed)
		ComputeCoefficients_EFIE(ri, zi, txi, tzi, dsi, rj, zj, txj, tzj,
                                         dsj, prodr, prod_ds, 
					 im, m2, G1_efie, Gcos_efie, Gsin_efie, 
					 A1, A2, A3, A4, A5, A6, A7, A8, A9);
	      
	      if (mfie_to_be_computed)
		ComputeCoefficients_MFIE(ri, zi, txi, tzi, dsi, rj, zj, txj, tzj, dsj,
                                         diffz, prodr, 
					 G1_mfie, Gcos_mfie, Gsin_mfie, F11, F12, F21, F22);
	      
	      if (this->alpha_efie != zero)
		{
		  zc = this->alpha_efie;
		  // add interactions on local matrices
		  UpdateEFIEMatrices(i, j, zc*A1, zc*A2, zc*A3, zc*A4, zc*A5,
                                     zc*A6, zc*A7, zc*A8, zc*A9,
				     Val_PhiX, Grad_PhiX, Val_PhiX, ValPhi_Sy_square,
                                     GradPhi_Sy_square, ValPhi_Sy_square,
				     mat_tt, mat_bb, mat_tb, mat_bt);
		  
		  zc *= prod_ds;
		  if (this->coupling_integral_equation == this->BIE_FEM)
		    UpdateMFIEMatrices(i, j, zc*F11, zc*F12, zc*F21, zc*F22, Val_PhiX, Val_PsiX,
				       ValPhi_Sy_square, ValPsi_Sy_square,
                                       coup11, coup22, coup12, coup21);
		  
		} // end if this->alpha_efie different from 0
	      
	      if (this->beta_mfie != zero)
		{
		  zc = this->beta_mfie;
		  if (this->coupling_integral_equation == this->BIE_FEM)
		    UpdateEFIEMatrices(i, j, zc*A1, zc*A2, zc*A3, zc*A4,
                                       zc*A5, zc*A6, zc*A7, zc*A8, zc*A9,
				       Val_PhiX, Grad_PhiX, Val_PsiX, ValPhi_Sy_square,
                                       GradPhi_Sy_square, ValPsi_Sy_square,
				       coup11, coup22, coup12, coup21);
		  
		  zc *= prod_ds;
		  // add interactions on local matrices
		  UpdateMFIEMatrices(i, j, zc*F11, zc*F12, zc*F21, zc*F22, Val_PhiX, Val_PhiX,
				     ValPhi_Sy_square, ValPhi_Sy_square,
                                     mat_tt, mat_bb, mat_tb, mat_bt);
		  
		}
	      
	      // DISP(mat_tt(0,0)); DISP(mat_bb(0,0));
	    }
      } // end loop on quadrature points
  }

  
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  PerformSingularIntegration(const Real_wp& phi1, Globatto<Real_wp>& gauss_X,
                             Globatto<Real_wp>& gauss_Y,
			     VectR2& CoorQuadY, VectReal_wp& WeightsQuadY,
			     int offset1, int offset2, Matrix<Complex_wp>& mat_tt,
                             Matrix<Complex_wp>& mat_bb,
			     Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt,
                             Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22,
			     Matrix<Complex_wp>& coup12, Matrix<Complex_wp>& coup21)
  {
    int nb_quad_sing = WeightsQuadY.GetM();
    // singular integration
    // adapted methods
    WeightsQuadY.Fill(0);
    int Nx = gauss_X.GetNbPointsQuad();
    Real_wp xi, xj, ri, zi, rj, zj, tzi, tzj, dsi, dsj, prod_ds, diffr, diffz, prodr, dist;
    Real_wp wn, txi, txj, phi, cos_phi, sin_phi, cos_mphi, sin_mphi, R;
    Complex_wp Gcos_efie, Gsin_efie, G1_efie, Gcos_mfie, Gsin_mfie, G1_mfie, G, zc;
    Complex_wp A1, A2, A3, A4, A5, A6, A7, A8, A9, F11, F12, F21, F22, zero(0);
    Complex_wp im(0, TypeEquation::number_mode);
    Complex_wp m2 = square(TypeEquation::number_mode);
    
    bool efie_to_be_computed = false, mfie_to_be_computed = false;
    if (this->alpha_efie != zero)
      efie_to_be_computed = true;
    if (this->beta_mfie != zero)
      mfie_to_be_computed = true;

    //int r = var_harmonic.mesh.GetOrder();
    //const TypeElement& Fb = var_harmonic.GetFaceBasis(r);
    VectReal_wp Value_PhiSx(nb_dof_Jt), Gradient_PhiSx(nb_dof_Jt), Value_PsiSx(nb_dof_Jt);
    VectReal_wp Value_PhiSy(nb_dof_Jt), Gradient_PhiSy(nb_dof_Jt), Value_PsiSy(nb_dof_Jt);
    
    Complex_wp iw = Iwp*var_harmonic.GetOmega();
    for (int i = 0; i < Nx; i++)
      {
	xi = gauss_X.Points(i);
	ComputeInterpolatedVal(ri, zi, txi, tzi, dsi, lob_IE,
			       Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE,
                               Ds_lob_IE, offset1, xi);
	
        /* 
	for (int k = 0; k < nb_dof_Jt; k++)
	  {
	    Value_PhiSx(k) = Fb.GetValue_ScalBoundary(k, xi);
	    Gradient_PhiSx(k) = Fb.GetGradient_ScalBoundary(k, xi);
	  }
	
	for (int k = 0; k < nb_dof_Jb; k++)
	  Value_PsiSx(k) = Fb.GetValue_VecBoundary(k, xi);
        */
        abort();
	Real_wp xs(0);
	// now we compute integration and weights points adapted to a singularity 
	// centered on point (ri,0,zi)
	// second method : duffy transformation
	if (offset1 == offset2)
	  xs = xi;
	else if (offset1 == (offset2-lob_IE.GetNbPointsQuad()))
	  {
	    xs = 0.0;
	  }
	else
	  xs = 1.0;
	
	GetQuadratureFormulas_Singularity(xs, phi1, mesh_sing, gauss_Y, CoorQuadY, WeightsQuadY);
	
	for (int n = 0; n < nb_quad_sing; n++)
	  if (abs(WeightsQuadY(n)) != 0)
	    {
	      xj = CoorQuadY(n)(0);
	      phi = CoorQuadY(n)(1);
	      ComputeInterpolatedVal(rj, zj, txj, tzj, dsj, lob_IE,
				     Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE,
                                     Ds_lob_IE, offset2, xj);
              
	      diffr = rj - ri;
	      diffz = zj - zi;
	      prodr = ri*rj; prod_ds = dsi*dsj;
	      dist = diffr*diffr + diffz*diffz;
	      cos_phi = cos(phi); sin_phi = sin(phi);
	      cos_mphi = cos(TypeEquation::number_mode*phi);
              sin_mphi = sin(TypeEquation::number_mode*phi);
	      
              // dist = \sqrt{ (r'- r)^2 + (z'-z)^2 + 2*r*r'(1.0-cos \phi) } 
	      R = sqrt(dist+2.0*prodr*(1.0-cos_phi));
	      // weight of integration
	      wn = WeightsQuadY(n);
	      
	      // green function
	      G = exp(iw*R)/R*wn;
	      G1_efie = -Iwp*G*cos_mphi;
	      Gcos_efie = -Iwp*G*cos_phi*cos_mphi;
	      Gsin_efie = G*sin_phi*sin_mphi;
	      
	      G *= (iw/R - 1.0/(R*R));
	      G1_mfie = G*cos_mphi;
	      Gcos_mfie = G*cos_phi*cos_mphi;
	      Gsin_mfie = Iwp*G*sin_phi*sin_mphi;
	      
	      wn = gauss_X.Weights(i)/(2.0*pi_wp);
	      
	      G1_efie *= wn; Gcos_efie *= wn; Gsin_efie *= wn;
	      G1_mfie *= wn; Gcos_mfie *= wn; Gsin_mfie *= wn;
	      
              /*
	      for (int k = 0; k < nb_dof_Jt; k++)
		{
		  Value_PhiSy(k) = Fb.GetValue_ScalBoundary(k, xj);
		  Gradient_PhiSy(k) = Fb.GetGradient_ScalBoundary(k, xj);
		}
	      
	      for (int k = 0;  k < nb_dof_Jb; k++)
		Value_PsiSy(k) = Fb.GetValue_VecBoundary(k, xj);
              */
              abort();
	      if (efie_to_be_computed)
		{
		  ComputeCoefficients_EFIE(ri, zi, txi, tzi, dsi, rj, zj, txj,
                                           tzj, dsj, prodr, prod_ds, 
					   im, m2, G1_efie, Gcos_efie, Gsin_efie, 
					   A1, A2, A3, A4, A5, A6, A7, A8, A9);
		}
	      
	      if (mfie_to_be_computed)
		{
		  ComputeCoefficients_MFIE(ri, zi, txi, tzi, dsi, rj, zj, txj,
                                           tzj, dsj, diffz, prodr, 
					   G1_mfie, Gcos_mfie, Gsin_mfie, F11, F12, F21, F22);
		  
		}
	      
	      if (this->alpha_efie != zero)
		{
		  zc = this->alpha_efie;
		  
		  UpdateEFIEMatrices(zc*A1, zc*A2, zc*A3, zc*A4, zc*A5, zc*A6, zc*A7, zc*A8, zc*A9,
				     Value_PhiSx, Gradient_PhiSx, Value_PhiSx, Value_PhiSy,
                                     Gradient_PhiSy, Value_PhiSy,
				     mat_tt, mat_bb, mat_tb, mat_bt);
		  zc *= prod_ds;
		  if (this->coupling_integral_equation == this->BIE_FEM)
		    UpdateMFIEMatrices(zc*F11, zc*F12, zc*F21, zc*F22, Value_PhiSx, Value_PsiSx,
				       Value_PhiSy, Value_PsiSy, coup11, coup22, coup12, coup21);
		  
		}
	      
	      if (this->beta_mfie != zero)
		{
		  zc = this->beta_mfie;
		  if (this->coupling_integral_equation == this->BIE_FEM)
		    UpdateEFIEMatrices(zc*A1, zc*A2, zc*A3, zc*A4, zc*A5,
                                       zc*A6, zc*A7, zc*A8, zc*A9,
				       Value_PhiSx, Gradient_PhiSx, Value_PsiSx,
                                       Value_PhiSy, Gradient_PhiSy, Value_PsiSy,
				       coup11, coup22, coup12, coup21);
		  
		  zc *= prod_ds;
		  UpdateMFIEMatrices(zc*F11, zc*F12, zc*F21, zc*F22, Value_PhiSx, Value_PhiSx,
				     Value_PhiSy, Value_PhiSy, mat_tt, mat_bb, mat_tb, mat_bt);
		  
		}
	      
	    } // end loop on quadrature points over j
      } // end loop on quadrature points over i
  }
  
  
  //! Initialization of some variables before computation of the matrix
  /*!
    mainly variables associated to geometry are computed :
    ri (radius), tx and tz (tangential vector, dsj (lineic element of integration) ...
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::InitComputation_IntegralEquation()
  {
    // different points of integration 
    gauss_reg_IE.ConstructQuadrature(this->order_integration_regular,
                                     gauss_reg_IE.QUADRATURE_GAUSS);
    gauss_sy_IE.ConstructQuadrature(this->order_integrationY, gauss_sy_IE.QUADRATURE_GAUSS);
    gauss_sx_IE.ConstructQuadrature(this->order_integrationX, gauss_sx_IE.QUADRATURE_GAUSS);
    gauss_syPhi_IE.ConstructQuadrature(this->order_integrationY,
                                       gauss_sy_IE.QUADRATURE_GAUSS_SQUARED);
    gauss_rx_IE.ConstructQuadrature(this->order_integration_radiusX,
                                    gauss_sy_IE.QUADRATURE_GAUSS_SQUARED);
    gauss_ry_IE.ConstructQuadrature(this->order_integration_radiusY,
                                    gauss_sx_IE.QUADRATURE_GAUSS_SQUARED);
    
    // DISP(gauss_reg_IE.Points); DISP(gauss_reg_IE.Weights); DISP(gauss_reg_IE.Cte_of_G);
    int nb = 0;
    Real_wp r, z, tx, tz, dsj;
    
    int Nb_pts = gauss_reg_IE.GetNbPointsQuad()*nb_edges_integral;
    Ri_gauss_IE.Reallocate(Nb_pts); Zi_gauss_IE.Reallocate(Nb_pts);
    Tx_gauss_IE.Reallocate(Nb_pts); Tz_gauss_IE.Reallocate(Nb_pts);
    Ds_gauss_IE.Reallocate(Nb_pts); Weight_gauss_IE.Reallocate(Nb_pts);

    for (int i = 0; i < nb_edges_integral; i++)
      for (int k = 0; k < gauss_reg_IE.GetNbPointsQuad(); k++)
	{
	  int offset1 = i*gauss_reg_IE.GetNbPointsQuad();
	  ComputeInterpolatedVal(r, z, tx, tz, dsj, lob_IE, Ri_lob_IE, 
				 Zi_lob_IE, Tx_lob_IE, Tz_lob_IE,
				 Ds_lob_IE, offset1, gauss_reg_IE.Points(k));
	  
	  this->Ri_gauss_IE(nb) = r;
	  this->Zi_gauss_IE(nb) = z;
	  this->Tx_gauss_IE(nb) = tx;
	  this->Tz_gauss_IE(nb) = tz;
	  this->Ds_gauss_IE(nb) = dsj;
	  this->Weight_gauss_IE(nb) = dsj*r*gauss_reg_IE.Weights(k);
	  nb++;
	}
    
    if (var_harmonic.print_level >= 0)
      {
	cout<<"Mean length of an edge "<<space_step_IE<<endl;
	cout<<"Number of dofs on integral equation "<<this->nodl_integral<<endl;
      }
    
    // now we compute \phi( \xi_j) d \phi/dx (\xi_j) and \psi(\xi_j)
    // where \phi are basis functions for Jt and psi basis functions for Jb
    // \phi are to be continuous (H^1) whereas \psi can be discontinuous (L^2)
    // R means regular (gauss_reg), Sx singular with integration over x (gauss_sx)
    // Sy singular with integration over y (gauss_sy), L Lobatto (lobs) 
    Val_PhiR_IE.Reallocate(nb_dof_Jt, gauss_reg_IE.GetNbPointsQuad());
    Grad_PhiR_IE.Reallocate(nb_dof_Jt, gauss_reg_IE.GetNbPointsQuad());
    Val_PsiR_IE.Reallocate(nb_dof_Jb, gauss_reg_IE.GetNbPointsQuad());
    
    Val_PhiSx_IE.Reallocate(nb_dof_Jt, gauss_sx_IE.GetNbPointsQuad());
    Grad_PhiSx_IE.Reallocate(nb_dof_Jt, gauss_sx_IE.GetNbPointsQuad());
    Val_PsiSx_IE.Reallocate(nb_dof_Jb, gauss_sx_IE.GetNbPointsQuad());
    
    Val_PhiSy_IE.Reallocate(nb_dof_Jt, gauss_sy_IE.GetNbPointsQuad());
    Grad_PhiSy_IE.Reallocate(nb_dof_Jt, gauss_sy_IE.GetNbPointsQuad());
    Val_PsiSy_IE.Reallocate(nb_dof_Jb, gauss_sy_IE.GetNbPointsQuad());
    
    Val_PhiRx_IE.Reallocate(nb_dof_Jt, gauss_rx_IE.GetNbPointsQuad());
    Grad_PhiRx_IE.Reallocate(nb_dof_Jt, gauss_rx_IE.GetNbPointsQuad());
    Val_PsiRx_IE.Reallocate(nb_dof_Jb, gauss_rx_IE.GetNbPointsQuad());
    
    Val_PhiRy_IE.Reallocate(nb_dof_Jt, gauss_ry_IE.GetNbPointsQuad());
    Grad_PhiRy_IE.Reallocate(nb_dof_Jt, gauss_ry_IE.GetNbPointsQuad());
    Val_PsiRy_IE.Reallocate(nb_dof_Jb, gauss_ry_IE.GetNbPointsQuad());
    
    //int order = var_harmonic.mesh.GetOrder();
    //const TypeElement& Fb = var_harmonic.GetFaceBasis(order);
    for (int i = 0; i < nb_dof_Jb; i++)
      {
	/* for (int j = 0; j < gauss_reg_IE.GetNbPointsQuad(); j++)
	  Val_PsiR_IE(i, j) = Fb.GetValue_VecBoundary(i, gauss_reg_IE.Points(j));
	
	for (int j = 0; j < gauss_sx_IE.GetNbPointsQuad(); j++)
	  Val_PsiSx_IE(i, j) = Fb.GetValue_VecBoundary(i, gauss_sx_IE.Points(j));
	
	for (int j = 0; j < gauss_sy_IE.GetNbPointsQuad(); j++)
	  Val_PsiSy_IE(i, j) = Fb.GetValue_VecBoundary(i, gauss_sy_IE.Points(j));
	
	for (int j = 0; j < gauss_rx_IE.GetNbPointsQuad(); j++)
	  Val_PsiRx_IE(i, j) = Fb.GetValue_VecBoundary(i, gauss_rx_IE.Points(j));
	
	for (int j = 0; j < gauss_ry_IE.GetNbPointsQuad(); j++)
        Val_PsiRy_IE(i,j) = Fb.GetValue_VecBoundary(i, gauss_ry_IE.Points(j)); */
        abort();
      }
    
    for (int i = 0; i < nb_dof_Jt; i++)
      {
	/* for (int j = 0; j < gauss_reg_IE.GetNbPointsQuad(); j++)
	  {
	    Val_PhiR_IE(i, j) = Fb.GetValue_ScalBoundary(i, gauss_reg_IE.Points(j));
	    Grad_PhiR_IE(i, j) = Fb.GetGradient_ScalBoundary(i, gauss_reg_IE.Points(j));
	  }
	
	for (int j = 0; j < gauss_sx_IE.GetNbPointsQuad(); j++)
	  {
	    Val_PhiSx_IE(i,j) = Fb.GetValue_ScalBoundary(i, gauss_sx_IE.Points(j));
	    Grad_PhiSx_IE(i,j) = Fb.GetGradient_ScalBoundary(i, gauss_sx_IE.Points(j));
	  }
	
	for (int j = 0; j < gauss_sy_IE.GetNbPointsQuad(); j++)
	  {
	    Val_PhiSy_IE(i,j) = Fb.GetValue_ScalBoundary(i, gauss_sy_IE.Points(j));
	    Grad_PhiSy_IE(i,j) = Fb.GetGradient_ScalBoundary(i, gauss_sy_IE.Points(j));
	  }
	
	for (int j = 0; j < gauss_rx_IE.GetNbPointsQuad(); j++)
	  {
	    Val_PhiRx_IE(i,j) = Fb.GetValue_ScalBoundary(i, gauss_rx_IE.Points(j));
	    Grad_PhiRx_IE(i,j) = Fb.GetGradient_ScalBoundary(i, gauss_rx_IE.Points(j));
	  }
	
	for (int j = 0; j < gauss_ry_IE.GetNbPointsQuad(); j++)
	  {
	    Val_PhiRy_IE(i,j) = Fb.GetValue_ScalBoundary(i, gauss_ry_IE.Points(j));
	    Grad_PhiRy_IE(i,j) = Fb.GetGradient_ScalBoundary(i, gauss_ry_IE.Points(j));
            } */
        abort();
	
      }
    
    if (var_harmonic.print_level >= 7)
      {
	DISP(Val_PhiR_IE); DISP(Grad_PhiR_IE); DISP(Val_PsiR_IE);
      }
    
    // if there is coupling with finite element, we need of mass matrices
    int Nt = nodl_Jt, Nb = nodl_Jb;
    // if (coupling_integral_equation == 1)
    ComputeMassMatrix_IntegralEquation(Nt, Nb, lob_IE, gauss_reg_IE, Ri_lob_IE, Zi_lob_IE,
				       Tx_lob_IE, Tz_lob_IE, Ds_lob_IE, Val_PhiR_IE, Val_PsiR_IE);
    
    // variables used for singularity integration
    if (this->singular_integration == QUADRATURE_HEXA_POLAR)
      ComputeQuadratureFormula_AxiSingularity(gauss_sy_IE, Points3D_sing, Weights3D_sing);
        
    // mesh adapted to the singularity (Duffy transformation with triangles)
    mesh_sing.ReallocateVertices(5);
    mesh_sing.Vertex(0).Init(0.0,0.0);
    mesh_sing.Vertex(1).Init(0.5,0.0);
    mesh_sing.Vertex(2).Init(1.0,0.0);	  
    mesh_sing.Vertex(3).Init(0.0,1.0);
    mesh_sing.Vertex(4).Init(1.0,1.0);
       
    // triangles
    mesh_sing.ReallocateElements(3);
    // first triangle
    mesh_sing.Element(0).InitTriangular(1, 3, 0, 1);
    // second triangle
    mesh_sing.Element(1).InitTriangular(1, 4, 3, 1);
    // third triangle
    mesh_sing.Element(2).InitTriangular(1, 2, 4, 1);
  
  } // end initialization for integral equation
  
  
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement, TypeEquation>
  ::AddMassTerms(Globatto<Real_wp>& gauss_reg, Matrix<Real_wp>& Val_PhiR,
                 Matrix<Real_wp>& Val_PsiR, int offset1,
                 Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
                 Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22)
  {
    Real_wp ri, zi, txi, tzi, dsi, wn, xi;
    Complex_wp zero(0);
    for (int i = 0; i < gauss_reg.GetNbPointsQuad(); i++)
      {
	xi = gauss_reg.Points(i);
	ComputeInterpolatedVal(ri, zi, txi, tzi, dsi, lob_IE,
			       Ri_lob_IE, Zi_lob_IE, Tx_lob_IE, Tz_lob_IE, Ds_lob_IE, offset1, xi);
	
	wn = gauss_reg.Weights(i)*ri*dsi;
	if (this->beta_mfie != zero)
	  {
	    // interactions Jt <-> Jt
	    for (int k = 0; k < nb_dof_Jt; k++)
	      for (int l = 0; l < nb_dof_Jt; l++)
		{
		  mat_tt(k,l) += 0.5*Val_PhiR(k,i)*Val_PhiR(l,i)*wn*this->beta_mfie;
		}
	    
	    // interactions Jb <-> Jb
	    for (int k = 0; k < nb_dof_Jt ; k++)
	      for (int l = 0; l < nb_dof_Jt; l++)
		{
		  mat_bb(k,l) += 0.5*Val_PhiR(k,i)*Val_PhiR(l,i)*wn*this->beta_mfie;
		} 
	  }
	
	if ((this->alpha_efie != zero)&&(this->coupling_integral_equation == this->BIE_FEM))
	  {
	    // interactions Jt <-> Jt
	    for (int k = 0; k < nb_dof_Jt; k++)
	      for (int l = 0; l < nb_dof_Jt; l++)
		{
		  coup11(k,l) += 0.5*Val_PhiR(k,i)*Val_PhiR(l,i)*wn*this->alpha_efie;
		}
	    
	    // interactions Jb <-> Jb
	    for (int k = 0; k < nb_dof_Jt; k++)
	      for (int l = 0; l < nb_dof_Jb; l++)
		{
		  coup22(k,l) += 0.5*Val_PhiR(k,i)*Val_PsiR(l,i)*wn*this->alpha_efie;
		} 
	      }
      }
  }
  
  
  //! computation of matrix coming from integral equation
  template<class TypeElement, class TypeEquation> template<class MatrixIntegral>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>
  ::ComputeMatrixIntegralEquation(MatrixIntegral& mat_integral)
  {
    mat_integral.Reallocate(2*nodl_Jt, 2*nodl_Jt); mat_integral.Fill(0);
    Matrix<Complex_wp> M11, M12, M21, M22;
    ComputeMatrixIntegralEquation(mat_integral, M11, M12, M21, M22);
  }
  
  
  //! computation of matrix coming from integral equation
  template<class TypeElement, class TypeEquation> template<class MatrixIntegral>
  void VarIntegralEquationAxi<TypeElement, TypeEquation>::
  ComputeMatrixIntegralEquation(MatrixIntegral& mat_integral, Matrix<Complex_wp>& M11,
                                Matrix<Complex_wp>& M12,
				 Matrix<Complex_wp>& M21, Matrix<Complex_wp>& M22)
  {    
    int order = var_harmonic.mesh.GetOrder();
    if (this->automatic_integration_phi)
      {
	// general rule for number of quadrature points for the integration over phi
	this->order_integration_phi = int(2.0*(0.5*Real_wp(abs(TypeEquation::number_mode))
                                               +4.0*var_harmonic.GetFrequency()
                                               *this->var_harmonic.GetXmax()));
	DISP(this->order_integration_phi);
      }
    
    if (this->automatic_integration_radius)
      {
	// general rule for quadrature points for integrals near radius
	this->order_integration_radiusY = max(2*abs(TypeEquation::number_mode), 3*order);
	if (gauss_ry_IE.GetOrder() != this->order_integration_radiusY)
	  gauss_ry_IE.ConstructQuadrature(this->order_integration_radiusY,
                                          gauss_sx_IE.QUADRATURE_GAUSS);
	
	DISP(this->order_integration_radiusY);
      }
    
    if (gauss_phi_IE.GetOrder() != this->order_integration_phi)
      gauss_phi_IE.ConstructQuadrature(this->order_integration_phi, gauss_sx_IE.QUADRATURE_GAUSS);
    
    // local matrices
    Matrix<Complex_wp> mat_tt(nb_dof_Jt, nb_dof_Jt), coup11(nb_dof_Jt, nb_dof_Jt);
    Matrix<Complex_wp> mat_bb(nb_dof_Jt, nb_dof_Jt), coup22(nb_dof_Jt, nb_dof_Jb);
    Matrix<Complex_wp> mat_tb(nb_dof_Jt, nb_dof_Jt), coup12(nb_dof_Jt, nb_dof_Jb);
    Matrix<Complex_wp> mat_bt(nb_dof_Jt, nb_dof_Jt), coup21(nb_dof_Jt, nb_dof_Jt);
        
    // allocation of coupling matrices if asked
    int Nt = nodl_Jt, Nb = nodl_Jb;
    int offset_Jt = 0; int offset_Jb = Nt;
    if (this->coupling_integral_equation == this->BIE_FEM)
      {
	M11.Reallocate(Nt, Nt); M12.Reallocate(Nt, Nb); M11.Fill(0); M12.Fill(0);
	M21.Reallocate(Nt, Nt); M22.Reallocate(Nt, Nb); M21.Fill(0); M22.Fill(0);
	offset_Jt += Nt+Nb; offset_Jb += Nt+Nb;
      }
    
    int nb_quad_sing(0);
    if (this->singular_integration != QUADRATURE_HEXA_POLAR)
      nb_quad_sing = 4*gauss_ry_IE.GetNbPointsQuad()*gauss_ry_IE.GetNbPointsQuad();
    
    VectReal_wp WeightsQuadY; VectR2 CoorQuadY;
    CoorQuadY.Reallocate(nb_quad_sing); CoorQuadY.Fill(R2(0,0));
    WeightsQuadY.Reallocate(nb_quad_sing); WeightsQuadY.Fill(0);
    
    Real_wp r1, r2, z1, z2, phi1;
    // Now we evaluate matrix
    // double loop on edges
    int old_percent = 0, new_percent = 0;    
    for (int i1 = 0; i1 < nb_edges_integral; i1++)
      for (int i2 = 0; i2 < nb_edges_integral; i2++)
	{
	  // displaying progression bar
	  new_percent = toInteger(round(Real_wp(i1)/(nb_edges_integral-1)*80));
	  for (int percent = old_percent; percent < new_percent; percent++)
	    { cout<<"#"; cout.flush(); }
	  
	  old_percent = new_percent;
	
	  // axial singularity ?
	  bool singularity_axis = false;
	  if ( abs(i1 - i2) <= 1)
	    if (EdgeIE_On_Axe(i1)||EdgeIE_On_Axe(i2))
	      singularity_axis = true;
	  
	  // initialization of local matrices
	  mat_tt.Fill(0); mat_bb.Fill(0); mat_bt.Fill(0); mat_tb.Fill(0);
	  if (this->coupling_integral_equation == this->BIE_FEM)
	    {
	      coup11.Fill(0); coup22.Fill(0); coup21.Fill(0); coup12.Fill(0);
	    }
	  
	  // offsets needed to access arrays Ri_lob, Zi_lob ...
	  int Nlob = lob_IE.GetNbPointsQuad();
	  int offset1 = Nlob*i1;
	  int offset2 = Nlob*i2;
	  
	  // extremities of the edge i1 (r1,z1) and of the edge i2 (r2,z2)
	  r1 = Ri_lob_IE(offset1+0); r2 = Ri_lob_IE(offset1+Nlob-1);
	  z1 = Zi_lob_IE(offset1+0); z2 = Zi_lob_IE(offset1+Nlob-1);
	  
	  // if i1 is near from i2, we compute number of subintervals for integration over phi
	  // singular integration -> we integrate regular part from phi1 to \pi 
          // (only one subdivision)
	  phi1 = EvaluatePhi1(r1, z1, r2, z2, var_harmonic.GetXmax());
	  
	  if ( abs(i1 - i2) > 1)
	    {
	      // no singularity	  
	      // regular integration from 0 to pi
	      PerformRegularIntegration(0.0, gauss_reg_IE, gauss_reg_IE, gauss_phi_IE,
                                        Val_PhiR_IE, Val_PsiR_IE, Grad_PhiR_IE,
					Val_PhiR_IE, Val_PsiR_IE, Grad_PhiR_IE,
                                        offset1, offset2, mat_tt, mat_bb, mat_tb, mat_bt,
					coup11, coup22, coup12, coup21);
	    }
	  else
	    {
	      // singularity
	      // regular integration from phi1 to pi
	      if ((singularity_axis)&&(i1==i2))
		{
		  // no regular part
		}
	      else
		PerformRegularIntegration(phi1, gauss_reg_IE, gauss_reg_IE, gauss_phi_IE,
                                          Val_PhiR_IE, Val_PsiR_IE, Grad_PhiR_IE,
					  Val_PhiR_IE, Val_PsiR_IE, Grad_PhiR_IE,
                                          offset1, offset2, mat_tt, mat_bb, mat_tb, mat_bt,
					  coup11, coup22, coup12, coup21);
	    }
	  
	  if ( abs(i1 - i2) <= 1)
	    {
	      if (this->singular_integration == QUADRATURE_GAUSS_SQUARE)
		{
		  if ((singularity_axis)&&(i1==i2))
		    PerformGaussSquaredIntegration(pi_wp, gauss_rx_IE, gauss_ry_IE,
                                                   gauss_syPhi_IE, Val_PhiRx_IE, Val_PsiRx_IE,
                                                   Grad_PhiRx_IE, offset1, offset2,
                                                   mat_tt, mat_bb, mat_tb, mat_bt,
                                                   coup11, coup22, coup12, coup21);
		  else
		    PerformGaussSquaredIntegration(phi1, gauss_sx_IE, gauss_sy_IE, gauss_ry_IE,
                                                   Val_PhiSx_IE, Val_PsiSx_IE, Grad_PhiSx_IE,
						   offset1, offset2, mat_tt, mat_bb,
                                                   mat_tb, mat_bt, coup11, coup22, coup12, coup21);
		}
	      else
		{
		  if ((singularity_axis)&&(i1==i2))
		    PerformSingularIntegration(pi_wp, gauss_rx_IE, gauss_ry_IE,
                                               CoorQuadY, WeightsQuadY, 
					       offset1, offset2, mat_tt, mat_bb, mat_tb, mat_bt,
					       coup11, coup22, coup12, coup21);
		  else
		    PerformSingularIntegration(phi1, gauss_sx_IE, gauss_sy_IE,
                                               CoorQuadY, WeightsQuadY, 
					       offset1, offset2, mat_tt, mat_bb, mat_tb, mat_bt,
					       coup11, coup22, coup12, coup21);
	      
		}
	    }
	  
	  if (i1 == i2)
	    {
	      // mass matrix for mfie_formulation
	      AddMassTerms(gauss_reg_IE, Val_PhiR_IE, Val_PsiR_IE, offset1,
                           mat_tt, mat_bb, coup11, coup22);
	      
	    } // end if i1==i2
	  
	  // we add contributions to the full matrix
	  // interactions Jt <-> Jt
	  for (int k = 0; k<nb_dof_Jt; k++)
	    for (int l = 0; l<nb_dof_Jt; l++)
	      {
		int num_dof_Jt = GetNumberDofT(i1,k);
		int num_dof2_Jt = GetNumberDofT(i2,l);
		mat_integral(offset_Jt + num_dof_Jt, offset_Jt + num_dof2_Jt) += mat_bt(k,l);
	      }
	  
	  // interactions Jb <-> Jb
	  for (int k = 0; k < nb_dof_Jt ; k++)
	    for (int l = 0; l < nb_dof_Jt; l++)
	      {
		int num_dof_Jb = GetNumberDofT(i1,k);
		int num_dof2_Jb = GetNumberDofT(i2,l);
		mat_integral(offset_Jb + num_dof_Jb, offset_Jb + num_dof2_Jb) -= mat_tb(k,l);
	      }
		    
	  // interactions Jt <-> Jb
	  for (int k = 0; k < nb_dof_Jt; k++)
	    for (int l = 0; l < nb_dof_Jt; l++)
	      {
		int num_dof_Jt = GetNumberDofT(i1,k);
		int num_dof2_Jb = GetNumberDofT(i2,l);
		mat_integral(offset_Jt  + num_dof_Jt, offset_Jb + num_dof2_Jb) += mat_bb(k,l);
	      }
	  
	  // interactions Jb <-> Jt
	  for (int k = 0; k < nb_dof_Jt; k++)
	    for (int l = 0; l < nb_dof_Jt; l++)
	      {
		int num_dof_Jb = GetNumberDofT(i1,k);
		int num_dof2_Jt = GetNumberDofT(i2,l);
		mat_integral(offset_Jb + num_dof_Jb, offset_Jt + num_dof2_Jt) -= mat_tt(k,l);
	      }
	  
	  // coupling matrices
	  if (this->coupling_integral_equation == this->BIE_FEM)
	    {
	      // interactions Jt <-> Jt
	      for (int k = 0; k < nb_dof_Jt; k++)
		for (int l = 0; l < nb_dof_Jt; l++)
		  {
		    int num_dof_Jt = GetNumberDofT(i1,k);
		    int num_dof2_Jt = GetNumberDofT(i2,l);
		    M11(num_dof_Jt, num_dof2_Jt) -= coup11(k,l);
		  }
	      
	      // interactions Jb <-> Jb
	      for (int k = 0; k < nb_dof_Jt; k++)
		for (int l = 0; l < nb_dof_Jb; l++)
		  {
		    int num_dof_Jb = GetNumberDofT(i1,k);
		    int num_dof2_Jb = GetNumberDofB(i2,l) - nodl_Jt;
		    M22(num_dof_Jb, num_dof2_Jb) -= coup22(k,l);
		  }
	      
	      // interactions Jt <-> Jb
	      for (int k = 0; k < nb_dof_Jt; k++)
		for (int l = 0; l < nb_dof_Jb; l++)
		  {
		    int num_dof_Jt = GetNumberDofT(i1,k);
		    int num_dof2_Jb = GetNumberDofB(i2,l) - nodl_Jt;
		    M12(num_dof_Jt, num_dof2_Jb) -= coup12(k,l);
		  }
	      
	      // interactions Jb <-> Jt
	      for (int k = 0; k < nb_dof_Jt; k++)
		for (int l = 0; l < nb_dof_Jt; l++)
		  {
		    int num_dof_Jb = GetNumberDofT(i1,k);
		    int num_dof2_Jt = GetNumberDofT(i2,l);
		    M21(num_dof_Jb, num_dof2_Jt) -= coup21(k,l);
		  }
	    } // end if coupling_integral_equation
	  
	} // end double loop on edges
    
    cout<<endl;

    // storage of matrix if asked
    if (var_harmonic.storage_matrix_asked)
      {
	if (TypeEquation::number_mode < 0)
	  mat_integral.Write(string("Mh_m"+to_str(-TypeEquation::number_mode)+".dat"));
	else
	  mat_integral.Write(string("Mh_p"+to_str(TypeEquation::number_mode)+".dat"));
      }
  }
  
  
  //! computation of quadrature formulas adapted to a singularity
  /*!
    \param[in] xs point of the singularity (xs,0)
    \param[in] phi1 domain of integration is [0,1]x[0,phi1]
    \param[out] mesh_sing internal mesh used
    \param[in] gauss_Y 1-D quadrature formula to use
    \param[out] CoorQuadY 2-D coordinates of integration points
    \param[out] WeightsQuadY 2-D integration weights
   */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  GetQuadratureFormulas_Singularity(const Real_wp& xs, const Real_wp& phi1,
                                    Mesh<Dimension2>& mesh_sing,
				    Globatto<Real_wp>& gauss_Y, VectR2& CoorQuadY,
                                    VectReal_wp& WeightsQuadY)
  {
    VectR2 s(4); Matrix2_2 mat_dfj; Real_wp coef;
    Real_wp x,y, xj, wn, phi;
    int Ny = gauss_Y.GetNbPointsQuad();
    int nb_quad_sing = WeightsQuadY.GetM();
    if (this->singular_integration == QUADRATURE_DUFFY)
      {
	Real_wp x1 = 0.0, x2 = 1.0;
	if (xs < 0.5)
	  x2 = 2.0*xs;
	else
	  x1 = 2.0*xs-1.0;
	// update of coordinates of vertices
	mesh_sing.Vertex(0).Init(x1,0.0);
	mesh_sing.Vertex(1).Init(xs,0.0);
	mesh_sing.Vertex(2).Init(x2,0.0);
	mesh_sing.Vertex(3).Init(x1,phi1);
	mesh_sing.Vertex(4).Init(x2,phi1);
	
	// DISP(xs); // DISP(mesh_sing);
	int nb = 0;
	// loop on all triangles
	for (int j = 0 ; j < mesh_sing.GetNbElt(); j++)
	  {
	    // we get vertices of the triangle j
	    mesh_sing.GetVerticesElement(j, s); 
	    
	    // loop on quadrature points
	    for (int k = 0; k < Ny; k++)
	      for (int l = 0; l < Ny; l++)
		{
		  x = gauss_Y.Points(k);
		  y = gauss_Y.Points(l);
		  for (int m = 0; m < 2; m++)
		    {
		      // coordinates of quadrature point
		      CoorQuadY(nb)(m) = (1.0-x)*s(0)(m) + x*(1.0-y)*s(1)(m) + x*y*s(2)(m);
		      mat_dfj(m,0) = s(1)(m)-s(0)(m) + y*(s(2)(m)-s(1)(m));
		      mat_dfj(m,1) = x*(s(2)(m)-s(1)(m));
		    }
		  // weight of integration
		  WeightsQuadY(nb) = gauss_Y.Weights(k)*gauss_Y.Weights(l)*abs(Det(mat_dfj));
		  nb++; // next point
		}
	  } // end loop on triangles
	
	if (xs < 0.5)
	  {
	    for (int k = 0; k < Ny; k++)
	      for (int l = 0; l < Ny; l++)
		{
		  CoorQuadY(nb)(0) = x2*(1.0-gauss_Y.Points(k)) + gauss_Y.Points(k);
		  CoorQuadY(nb)(1) = gauss_Y.Points(l)*phi1;
		  WeightsQuadY(nb) = gauss_Y.Weights(k)*gauss_Y.Weights(l)*(1.0-x2)*phi1;
		  nb++;
		}
	  }
	else
	  {
	    for (int k = 0; k < Ny; k++)
	      for (int l = 0; l < Ny; l++)
		{
		  CoorQuadY(nb)(0) = x1*gauss_Y.Points(k);
		  CoorQuadY(nb)(1) = gauss_Y.Points(l)*phi1;
		  WeightsQuadY(nb) = gauss_Y.Weights(k)*gauss_Y.Weights(l)*x1*phi1;
		  nb++;
		}
	  }
	
	coef = 0.0;
	for (int n = 0; n < nb_quad_sing; n++)
	  coef += WeightsQuadY(n);
	
	
	// DISP(phi1);
	// cout<<"Aire of rectangle (1,0) x (0,phi1) "<<coef<<endl;
      }
    else if (this->singular_integration == QUADRATURE_POLAR)
      {
	Real_wp x1 = 0.0, x2 = 1.0;
	if (xs < 0.5)
	  x2 = 2.0*xs;
	else
	  x1 = 2.0*xs-1.0;
	
	Real_wp teta, Rmax, rj;
	// teta1 angle between vector (xs,0) (x2,0) and (xs,0) (x2,phi1)
	Real_wp teta1 = atan(phi1/(x2-xs));
	// teta2 angle between vector (xs,0) (x2,phi1) and (xs,0) (x1,phi1)
	Real_wp teta2 = 0.5*pi_wp+atan((xs-x1)/phi1);
	int nb = 0; // DISP(teta1); DISP(teta2);
	// loop on quadrature points
	for (int k = 0; k < Ny; k++)
	  for (int l = 0; l < Ny; l++)
	    {
	      // first triangle (xs,0) (x2,0) (x2,phi1)
	      teta = gauss_Y.Points(k)*teta1; // DISP(teta);
	      Rmax = (x2-xs)/cos(teta); // DISP(rj):
	      rj = Rmax*gauss_Y.Points(l);
	      CoorQuadY(nb)(0) = xs + rj*cos(teta);
	      CoorQuadY(nb)(1) = rj*sin(teta);
	      WeightsQuadY(nb) = rj*Rmax*teta1 * gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      nb++;
	      
	      // second triangle (xs,0) (x2,phi1) (x1,phi1)
	      teta = gauss_Y.Points(k)*teta2 + (1.0-gauss_Y.Points(k))*teta1;
	      Rmax = phi1/sin(teta);
	      rj = Rmax*gauss_Y.Points(l);
	      CoorQuadY(nb)(0) = xs + rj*cos(teta);
	      CoorQuadY(nb)(1) = rj*sin(teta);
	      WeightsQuadY(nb) = rj*Rmax*(teta2-teta1) * gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      nb++;
	      
	      // third triangle (xs,0) (x1,phi1) (x1,0)
	      teta = (1.0-gauss_Y.Points(k))*teta2 + pi_wp*gauss_Y.Points(k);
	      Rmax = (x1-xs)/cos(teta);
	      rj = Rmax*gauss_Y.Points(l);
	      CoorQuadY(nb)(0) = xs + rj*cos(teta);
	      CoorQuadY(nb)(1) = rj*sin(teta);
	      WeightsQuadY(nb) = rj*Rmax*(pi_wp-teta2) * gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      nb++;
	    }
	
	if (xs < 0.5)
	  {
	    for (int k = 0; k < Ny; k++)
	      for (int l = 0; l < Ny; l++)
		{
		  CoorQuadY(nb)(0) = x2*(1.0-gauss_Y.Points(k)) + gauss_Y.Points(k);
		  CoorQuadY(nb)(1) = gauss_Y.Points(l)*phi1;
		  WeightsQuadY(nb) = gauss_Y.Weights(k)*gauss_Y.Weights(l)*(1.0-x2)*phi1;
		  nb++;
		}
	  }
	else
	  {
	    for (int k = 0; k < Ny; k++)
	      for (int l = 0; l < Ny; l++)
		{
		  CoorQuadY(nb)(0) = x1*gauss_Y.Points(k);
		  CoorQuadY(nb)(1) = gauss_Y.Points(l)*phi1;
		  WeightsQuadY(nb) = gauss_Y.Weights(k)*gauss_Y.Weights(l)*x1*phi1;
		  nb++;
		}
	  }
	coef = 0.0;
	for (int n = 0; n < nb_quad_sing; n++)
	  coef += WeightsQuadY(n);
	
	Mlt(phi1/coef, WeightsQuadY);
	// DISP(phi1);
	// cout<<"Aire of rectangle (1,0) x (0,phi1) "<<coef<<endl;
      }
    else if (this->singular_integration == QUADRATURE_GAUSS)
      {
	Real_wp x1 = max(0.0, 2.0*xs-1.0);
	Real_wp x2 = min(1.0, 2.0*xs);
	
	int nb = 0;
	// loop on quadrature points
	for (int k = 0; k < Ny; k++)
	  for (int l = 0; l < Ny; l++)
	    {
	      wn = phi1*gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      xj = gauss_Y.Points(k)*x1; phi = gauss_Y.Points(l)*phi1;
	      CoorQuadY(nb).Init(xj,phi);
	      WeightsQuadY(nb++) = x1*wn;
	      
	      xj = (1.0-gauss_Y.Points(k))*x1 + gauss_Y.Points(k)*xs;
	      CoorQuadY(nb).Init(xj,phi);
	      WeightsQuadY(nb++) = (xs-x1)*wn;
	      
	      xj = (1.0-gauss_Y.Points(k))*xs + gauss_Y.Points(k)*x2;
	      CoorQuadY(nb).Init(xj, phi);
	      WeightsQuadY(nb++) = (x2-xs)*wn;
	      
	      xj = (1.0-gauss_Y.Points(k))*x2 + gauss_Y.Points(k);
	      CoorQuadY(nb).Init(xj, phi);
	      WeightsQuadY(nb++) = (1.0-x2)*wn;
	    }
      }
    else if (this->singular_integration == QUADRATURE_GAUSS_SQUARE)
      {
	Real_wp x1 = max(0.0, 2.0*xs-1.0);
	Real_wp x2 = min(1.0, 2.0*xs);
	
	int nb = 0;
	Real_wp xsi_square, phi_square;
	// loop on quadrature points
	for (int k = 0; k < Ny; k++)
	  for (int l = 0; l < Ny; l++)
	    {
	      wn = phi1*gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      xj = gauss_Y.Points(k)*x1;
	      phi = gauss_Y.Points(l)*phi1;
	      CoorQuadY(nb).Init(xj,phi);
	      WeightsQuadY(nb++) = x1*wn;
	      
	      xsi_square = gauss_Y.Points(k)*gauss_Y.Points(k);
	      phi_square = gauss_Y.Points(l)*gauss_Y.Points(l);
	      xj = xsi_square*x1 + (1.0-xsi_square)*xs;
	      phi = phi_square*phi1;
	      wn = 4.0*gauss_Y.Points(l)*phi1*gauss_Y.Points(k)
                *gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      
              CoorQuadY(nb).Init(xj,phi);
	      WeightsQuadY(nb++) = (xs-x1)*wn;
	      
	      xj = (1.0-xsi_square)*xs + xsi_square*x2;
	      CoorQuadY(nb).Init(xj, phi);
	      WeightsQuadY(nb++) = (x2-xs)*wn;
	      
	      wn = phi1*gauss_Y.Weights(k)*gauss_Y.Weights(l);
	      xj = (1.0-gauss_Y.Points(k))*x2 + gauss_Y.Points(k);
	      phi = gauss_Y.Points(l)*phi1;
	      CoorQuadY(nb).Init(xj, phi);
	      WeightsQuadY(nb++) = (1.0-x2)*wn;
	    }
      }
    
    // mesh_sing.Write("singularity.mesh");
    // CoorQuadY.WriteText("Points.dat");
  }
  
  
  //! computation of right hand side for integral equation
  /*!
    written only for an incident plane wave
    The result is stored in vector rhs_integral
  */
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>
  ::ComputeRightHandSideIntegralEquation(VectComplex_wp& rhs_integral, int num_angle)
  {
    Real_wp tx, tz, r, z, dsj, wn;
    R3_Complex_wp valE, valH;
    R2_Complex_wp Erz, Hrz; Complex_wp Eteta, Hteta, Eb, Et,zero(0);
    R2 point_glob;
    
    // initialization
    int Nt = nodl_Jt, Nb = nodl_Jb;
    int offset_Jt = 0, offset_Jb = Nt;
    if (this->coupling_integral_equation == this->BIE_FEM)
      {
	offset_Jt += Nt+Nb; offset_Jb = 2*Nt+Nb;
	rhs_integral.Reallocate(this->nodl_integral+2*nodl_Jt);
      }
    else
      rhs_integral.Reallocate(2*nodl_Jt);
    
    DiffractedWaveSource<TypeElement, PlaneWaveSource<Dimension3>, TypeEquation> f(var_harmonic);
    rhs_integral.Fill(0);
    
    // loop on edges 
    int nb = 0;
    int Ntot_quad = nb_edges_integral*gauss_reg_IE.GetNbPointsQuad();
    for (int i = 0; i < nb_edges_integral; i++)
      {
	// loop on quadrature points of edge
	for (int k = 0; k < gauss_reg_IE.GetNbPointsQuad(); k++)
	  {
	    r = Ri_gauss_IE(nb); z = Zi_gauss_IE(nb);
	    tx = Tx_gauss_IE(nb); tz = Tz_gauss_IE(nb);
	    dsj = Ds_gauss_IE(nb); wn = Weight_gauss_IE(nb);
	    
	    // evalution of E^{inc}, Eteta^{inc}, H^{inc} and Hteta^{inc}
	    point_glob(0) = r; point_glob(1) = z;
	    
	    f.Evaluate_E_H(point_glob, valE, valH, num_angle*Ntot_quad+nb);
	    Erz(0) = valE(0); Eteta = valE(1); Erz(1) = valE(2);
	    Hrz(0) = valH(0); Hteta = valH(1); Hrz(1) = valE(2);
	    
	    Et = zero; Eb = zero; 
	    if (this->alpha_efie != zero)
	      {		
		// source is equal to \int_\Gamma E^{inc} J^t
                // = \int_\Gamma E_t t \cdot J^t + E_b b \cdot J^t 
		// multiplication by I
		Eb = -Iwp*Eteta; // E_b term
		Et = (tx*Erz(0) + tz*Erz(1)); // E_t term
		
		Eb *= this->alpha_efie; Et *= this->alpha_efie;
		
		// Eb = -(tx*Erz(0) + tz*Erz(1));
		// Et = -Iwp*Eteta;
		
	      }
	    
	    if (this->beta_mfie != zero)
	      {
		// source is \int_\Gamma n \times H^{inc} J^t dx
		//        \int_\Gamma -H_t b \cdot J^t + H_b t \cdot J^t
		Et -= (Iwp*Hteta) * this->beta_mfie; // term H_b
		Eb -= ( tx*Hrz(0) + tz*Hrz(1)) * this->beta_mfie; // term -H_t
		
	      }
	    // DISP(tx); DISP(tz);
	    // DISP(Et); DISP(Eb);
	    // we add contributions to rhs_integral
	    for (int j = 0; j < nb_dof_Jt; j++)
	      {
		int num_dof_Jt = GetNumberDofT(i,j);
		rhs_integral(offset_Jt + num_dof_Jt) += wn*Eb*Val_PhiR_IE(j,k);
	      }    
	    
	    for (int j = 0; j < nb_dof_Jt; j++)
	      {
		int num_dof_Jb = GetNumberDofT(i,j);
		rhs_integral(offset_Jb + num_dof_Jb) -= wn*Et*Val_PhiR_IE(j,k);
	      }    
	    nb++;
	  }
      }
    
    if (var_harmonic.storage_matrix_asked)
      {
	if (TypeEquation::number_mode < 0)
	  rhs_integral.Write(string("Bh_m"+to_str(-TypeEquation::number_mode)+".dat"));
	else
	  rhs_integral.Write(string("Bh_p"+to_str(TypeEquation::number_mode)+".dat"));
      }
  }
  
    
  //! 3-D quadrature formula to take into account singularity
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeQuadratureFormula_AxiSingularity(const Globatto<Real_wp>& gauss_reg,
                                          VectR3& Points3D, VectReal_wp& Weights3D) const
  {
    int Nr = gauss_reg.GetNbPointsQuad();
    int N = Nr*Nr*Nr*8;
    Points3D.Reallocate(N); Weights3D.Reallocate(N);
    
    Real_wp r, teta, z, teta_init, teta_end, rmax; int nb=0;
    for (int i=0 ; i<Nr ; i++)
      {
	z = 0.5*gauss_reg.Points(i);
	// first case
	teta_init = 0; teta_end = atan(1.0/z);
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = z/cos(teta);
	    for (int k=0 ; k<Nr ; k++)
	    {
	      r = gauss_reg.Points(k)*rmax;
	      Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
	      Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                *gauss_reg.Weights(j)*gauss_reg.Weights(i);
	      nb++;
	    }
	  }
	teta_init = teta_end; teta_end = 0.5*pi_wp;
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = 1.0/sin(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	      }
	  }
	
	// second case
	teta_init = 0.5*pi_wp; teta_end = pi_wp - atan(1.0/z);
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = 1.0/sin(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	      }
	  }
	teta_init = teta_end; teta_end = pi_wp;
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = -z/cos(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	      }
	  }
	
	z = 0.5 + 0.5*gauss_reg.Points(i);
	// third case 
	teta_init = 0; teta_end = atan(1.0/ (1.0-z) );
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = (1.0-z)/cos(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	    }
	  }
	teta_init = teta_end; teta_end = 0.5*pi_wp;
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = 1.0/sin(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	      }
	  }
	
	// fourth case
	teta_init = 0.5*pi_wp; teta_end = pi_wp - atan(1.0/(1.0-z));
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = 1.0/sin(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	      }
	  }
	teta_init = teta_end; teta_end = pi_wp;
	for (int j=0 ; j<Nr ; j++)
	  {
	    teta = teta_init*(1.0-gauss_reg.Points(j)) + gauss_reg.Points(j)*teta_end;
	    rmax = -(1.0-z)/cos(teta);
	    for (int k=0 ; k<Nr ; k++)
	      {
		r = gauss_reg.Points(k)*rmax;
		Points3D(nb).Init(z+r*cos(teta), z-r*cos(teta), r*sin(teta));
		Weights3D(nb) = r*rmax*(teta_end-teta_init)*gauss_reg.Weights(k)
                  *gauss_reg.Weights(j)*gauss_reg.Weights(i);
		nb++;
	      }
	  }
      }
    Mlt(Real_wp(1.0/(2.0*pi_wp)), Weights3D);
  }
  
  //! retrieve J from solution 
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  GetComponentsSolution_IntegralEquation(VectComplex_wp& full_sol, VectComplex_wp& J)
  {
    J.Copy(full_sol);
  }
  
  //! retrieve J and K from solution 
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  GetComponentsSolution_CouplingSystem(VectComplex_wp& full_sol, VectComplex_wp& J,
                                       VectComplex_wp& K)
  {
    J.Reallocate(2*nodl_Jt); K.Reallocate(this->nodl_integral);
    int Nt = nodl_Jt, Nb = nodl_Jb;
    for (int i = 0; i < Nt; i++)
      {
	K(i) = full_sol(i);
	J(i) = full_sol(Nt+Nb+i);
	J(Nt+i) = full_sol(2*Nt+Nb+i);
      }
    
    for (int i = 0; i < Nb; i++)
      K(Nt+i) = full_sol(Nt+i);
  }
  
  //! computation of mass matrices
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeMassMatrix_IntegralEquation(int Nt, int Nb, const Globatto<Real_wp>& lob,
                                     const Globatto<Real_wp>& gauss_reg,
				     const VectReal_wp& Ri_lob, const VectReal_wp& Zi_lob,
                                     const VectReal_wp& Tx_lob, const VectReal_wp& Tz_lob,
				     const VectReal_wp& Ds_lob, const Matrix<Real_wp>& Val_PhiR,
                                     const Matrix<Real_wp>& Val_PsiR)
  {
    this->matmass_Jt_IE.Clear(); this->matcrossed_Ch1.Clear();
    this->matmass_Jt_IE.Reallocate(Nt, Nt); 
    this->matcrossed_Ch1.Reallocate(Nt, Nb);
    
    Real_wp val, ri, zi, txi, tzi, dsi, xi, wn;
    for (int i1 = 0; i1 < nb_edges_integral; i1++)
      {
	int offset1 = lob.GetNbPointsQuad()*i1;
	for (int i = 0; i < gauss_reg.GetNbPointsQuad(); i++)
	  {
	    xi = gauss_reg.Points(i);
	    this->ComputeInterpolatedVal(ri, zi, txi, tzi, dsi, lob,
					 Ri_lob, Zi_lob, Tx_lob, Tz_lob, Ds_lob, offset1, xi);
	    
	    wn = ri*dsi*gauss_reg.Weights(i);	    
	    // interactions Jt <-> Jt
	    for (int k = 0; k < nb_dof_Jt; k++)
	      for (int l = 0; l < nb_dof_Jt; l++)
		{
		  int num_dof_Jt = this->GetNumberDofT(i1,k);
		  int num_dof2_Jt = this->GetNumberDofT(i1,l);
		  val = Val_PhiR(k,i)*Val_PhiR(l,i)*wn;
		  this->matmass_Jt_IE.AddInteraction(num_dof_Jt, num_dof2_Jt, val);
		}
		  
	    for (int k = 0; k < nb_dof_Jb; k++)
	      for (int l = 0; l < nb_dof_Jt; l++)
		{
		  val = Val_PsiR(k,i)*Val_PhiR(l,i)*wn;
		  int num_dof_Jb = this->GetNumberDofB(i1,k) - nodl_Jt;
		  int num_dof2_Jt = this->GetNumberDofT(i1,l);
		  this->matcrossed_Ch1.AddInteraction(num_dof2_Jt, num_dof_Jb, val);
		}
	  }
      } // end loop on edges
  }
  
  
  //! construction of global matrix from blocks M11, M21, M12 and M22
  template<class TypeElement, class TypeEquation> template<class MatrixIntegral>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeMatrixCoupling_IntegralEquation(MatrixIntegral& mat_integral,
                                         const Matrix<Complex_wp>& M11,
                                         const Matrix<Complex_wp>& M12,
					 const Matrix<Complex_wp>& M21,
                                         const Matrix<Complex_wp>& M22)
  {
    int Nt = nodl_Jt, Nb = nodl_Jb;
    for (int i = 0; i < Nt; i++)
      for (int j = 0; j < Nt; j++)
	mat_integral(Nt+Nb+i, j) = M11(i,j);
    
    for (int i = 0; i < Nt; i++)
      for (int j = 0; j < Nb; j++)
	mat_integral(Nt+Nb+i, Nt+j) = M12(i,j);
    
    for (int i = 0; i < Nt; i++)
      for (int j = 0; j < Nt; j++)
	mat_integral(2*Nt+Nb+i, j) = M21(i,j);
    
    for (int i = 0; i < Nt; i++)
      for (int j = 0; j < Nb; j++)
	mat_integral(2*Nt+Nb+i, Nt+j) = M22(i,j);
  }
  
  
  //! finite element matrix included in matrix coming from integral equation
  template<class TypeElement, class TypeEquation>
  void VarIntegralEquationAxi<TypeElement,TypeEquation>::
  ComputeMatrixCoupling_FiniteElement(Matrix<Complex_wp>& mat_integral,
                                      const Matrix<Complex_wp, Symmetric, RowSymPacked>& 
                                      schur_matrix)
  {
    int Nt = nodl_Jt, Nb = nodl_Jb;
    // DISP(Nt); DISP(Nb);
    Complex_wp coef;
    
    VectReal_wp DsDofB_IE(Nb), RiDofB_IE(Nb), ZiDofB_IE(Nb);
    VectReal_wp DsDofT_IE(Nt), RiDofT_IE(Nt), ZiDofT_IE(Nt);
    Real_wp ri,zi, tx, tz, dsi, xi;
    int r = var_harmonic.mesh.GetOrder();
    const TypeElement& Fb = var_harmonic.GetFaceBasis(r);
    for (int i = 0; i < nb_edges_integral; i++)
      {
	int offset = i*this->lob_IE.GetNbPointsQuad();
	for (int j = 0; j < nb_dof_Jb; j++)
	  {
	    xi = Fb.Points1D_dof(j);
	    this->ComputeInterpolatedVal(ri, zi, tx, tz, dsi, this->lob_IE, Ri_lob_IE, 
					 Zi_lob_IE, Tx_lob_IE, Tz_lob_IE, Ds_lob_IE, offset, xi);
	    
	    int num_dof_Jb = this->GetNumberDofB(i,j)-Nt;
	    DsDofB_IE(num_dof_Jb) = dsi;
	    RiDofB_IE(num_dof_Jb) = ri;
	    ZiDofB_IE(num_dof_Jb) = zi;
	  }
	
	for (int j = 0; j < nb_dof_Jt; j++)
	  {
	    xi = Fb.PointsNodal1D(j);
	    this->ComputeInterpolatedVal(ri, zi, tx, tz, dsi, this->lob_IE, Ri_lob_IE, 
					 Zi_lob_IE, Tx_lob_IE, Tz_lob_IE, Ds_lob_IE, offset, xi);
	    
	    int num_dof_Jt = this->GetNumberDofT(i,j);
	    DsDofT_IE(num_dof_Jt) = dsi;
	    RiDofT_IE(num_dof_Jt) = ri;
	    ZiDofT_IE(num_dof_Jt) = zi;
	  }
      }
    DsDofB_IE.Write("dsj_Jb.dat");
    RiDofB_IE.Write("ri_Jb.dat");
    ZiDofB_IE.Write("zi_Jb.dat");
    DsDofT_IE.Write("dsj_Jt.dat");
    RiDofT_IE.Write("ri_Jt.dat");
    ZiDofT_IE.Write("zi_Jt.dat");
    
    // part A11
    for (int i = 0; i < Nt; i++)
      for (int j = 0; j < Nt; j++)
	mat_integral(i,j) = schur_matrix(i, j);
    
    // part Iwp*A12*dsb
    for (int i = 0; i < Nt; i++)
      for (int j = 0; j < Nb; j++)
	{
	  coef = Iwp*DsDofB_IE(j);
	  
	  if (SignVolDof_IE(Nt+j))
	    mat_integral(i,Nt+j) = coef*schur_matrix(i, Nt+j);
	  else
	    mat_integral(i,Nt+j) = -coef*schur_matrix(i, Nt+j);
	}
    
    // part Iwp*dsb*A21
    for (int i = 0; i < Nb; i++)
      for (int j = 0; j < Nt; j++)
	{
	  coef = Iwp*DsDofB_IE(i);
	  
	  if (SignVolDof_IE(Nt+i))
	    mat_integral(Nt+i, j) = coef*schur_matrix(Nt+i, j);
	  else
	    mat_integral(Nt+i, j) = -coef*schur_matrix(Nt+i, j);
	}
    
    // part -dsb*A22*dsb
    for (int i = 0; i < Nb; i++)
      for (int j = 0; j < Nb; j++)
	{
	  coef = -DsDofB_IE(i)*DsDofB_IE(j);
	  
	  if (SignVolDof_IE(Nt+i) == SignVolDof_IE(Nt+j))
	    mat_integral(Nt+i, Nt+j) = coef*schur_matrix(Nt+i, Nt+j);
	  else
	    mat_integral(Nt+i, Nt+j) = -coef*schur_matrix(Nt+i, Nt+j);
	}
    
    // part -Iwp*C_h^1 P_b  -Iwp*C_h^2 P_t
    for (int i = 0; i < Nt; i++)
      for (int j1 = 0; j1 < this->matmass_Jt_IE.GetRowSize(i); j1++)
	{
	  int j = this->matmass_Jt_IE.Index(i,j1);
	  coef = -Iwp*this->matmass_Jt_IE.Value(i,j1);
	  
	  mat_integral(i, 2*Nt+Nb+j) = coef;
	}
    
    for (int i = 0; i < Nt; i++)
      for (int j1 = 0; j1 < this->matcrossed_Ch1.GetRowSize(i); j1++)
	{
	  int j = this->matcrossed_Ch1.Index(i,j1);
	  coef = -Iwp*this->matcrossed_Ch1.Value(i,j1);
	  
	  mat_integral(Nt+j, Nt+Nb+i) = coef;
	}

  }
  


} // namespace Montjoie

#define MONTJOIE_FILE_INTEGRAL_EQUATION_MAXWELL_AXI_CXX
#endif
