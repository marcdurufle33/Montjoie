#ifndef MONTJOIE_FILE_NUMBER_MESH_CXX

namespace Montjoie
{

  /**********************
   * MeshNumbering_Base *
   **********************/
  
  
  //! default constructor
  template<class T>
  MeshNumbering_Base<T>::MeshNumbering_Base()
  {
    order = 1;
    type_variable_order = CONSTANT_ORDER;
    use_max_velocity_order = false;
    
    nodl = 0;
    plane_wave_quasi_periodic = false;
    treat_periodic_condition_during_number = true;
    compute_dof_pml = false;
    drop_interface_pml_dof = false;
    nb_dof_pml = 0;
    type_formulation_periodic = SAME_PERIODIC_DOFS;
  }

  
  //! destructor
  template<class T>
  MeshNumbering_Base<T>::~MeshNumbering_Base()
  {
  }
  

  template<class T>
  void MeshNumbering_Base<T>::CopyInputData(const MeshNumbering_Base<T>& mesh_num)
  {
    order = mesh_num.order;
    type_variable_order = mesh_num.type_variable_order;
    compute_dof_pml = mesh_num.compute_dof_pml;
    type_formulation_periodic = mesh_num.type_formulation_periodic;
    plane_wave_quasi_periodic = mesh_num.plane_wave_quasi_periodic;
    drop_interface_pml_dof = mesh_num.drop_interface_pml_dof;
  }
  
  //! allocates n neighboring elements, previous values are lost
  /*!
    Neighboring elements are elements that are adjacent to the current mesh,
    but that are belonging to other processors
   */
  template<class T>
  void MeshNumbering_Base<T>::ReallocateNeighborElements(int n)
  {
    order_neighbor_elt.Reallocate(n);
    type_neighbor_elt.Reallocate(n);
    num_edge_neighbor_elt.Reallocate(n);
    order_neighbor_elt.Fill(0);
    type_neighbor_elt.Fill(0); num_edge_neighbor_elt.Fill(-1);
  }

  
  template<class T>
  void MeshNumbering_Base<T>::FinalizeNeighborElements(int N)
  {
    inv_num_edge_neighbor_elt.Reallocate(N);
    inv_num_edge_neighbor_elt.Fill(-1);
    for (int i = 0; i < num_edge_neighbor_elt.GetM(); i++)
      inv_num_edge_neighbor_elt(num_edge_neighbor_elt(i)) = i;
  }


  //! returns size of memory used by the object in bytes
  template<class T>
  size_t MeshNumbering_Base<T>::GetMemorySize() const
  {
    size_t taille = nb_dof_element.GetMemorySize() + order_inside.GetMemorySize()
      + order_edge.GetMemorySize() + order_face.GetMemorySize()
      + order_element.GetMemorySize() + order_boundary.GetMemorySize()
      + order_neighbor_elt.GetMemorySize() + type_neighbor_elt.GetMemorySize()
      + num_edge_neighbor_elt.GetMemorySize() + inv_num_edge_neighbor_elt.GetMemorySize();
    
    taille += mesh_size_variable_order.GetMemorySize() + coef_elt_variable_order.GetMemorySize();
    taille += nb_points_quadrature_on_boundary.GetMemorySize() + IndexDofPML.GetMemorySize();
    taille += PeriodicityBoundary.GetMemorySize();
    taille += type_periodicity_dof.GetMemorySize();
    taille += type_periodicity_boundary.GetMemorySize();
    taille += PeriodicDofs.GetMemorySize();
    taille += PeriodicDofs_Original.GetMemorySize();
    taille += sizeof(void*)*size_t(dof_element.GetM());
    for (int i = 0; i < dof_element.GetM(); i++)
      taille += dof_element(i).GetMemorySize();
    
    taille += number_map.GetMemorySize();
    taille += OffsetDofFaceNumber.GetMemorySize();
    taille += OffsetDofEdgeNumber.GetMemorySize();
    taille += OffsetDofElementNumber.GetMemorySize();
    taille += OffsetQuadElementNumber.GetMemorySize();
    taille += GlobDofNumber_Subdomain.GetMemorySize();
    taille += GlobDofPML_Subdomain.GetMemorySize();
    return taille;
  }

  
  //! sets input parameters from a line of the data file
  template<class T>
  void MeshNumbering_Base<T>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("OrderDiscretization"))
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
	    if (parameters.GetM() <= 2)
	      {
		cout << "In SetInputData of MeshNumbering_Base" << endl;
		cout << "OrderDiscretization needs more parameters, for instance :" << endl;
		cout << "OrderDiscretization = MAX_EDGE AUTO coef" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    if (!parameters(1).compare("AUTO"))
	      {
		T coef = to_num<T>(parameters(2));
		mesh_size_variable_order.Reallocate(20);
		mesh_size_variable_order(0) = 0;
		mesh_size_variable_order(1) = 0;
		mesh_size_variable_order(2) = 0.02;
		mesh_size_variable_order(3) = 0.15;
		mesh_size_variable_order(4) = 0.25;
		mesh_size_variable_order(5) = 0.41;
		mesh_size_variable_order(6) = 0.54;
		mesh_size_variable_order(7) = 0.7;
		mesh_size_variable_order(8) = 0.87;
		for (int i = 9; i < mesh_size_variable_order.GetM(); i++)
		  mesh_size_variable_order(i) = Real_wp(i)/8;
		
		Mlt(coef, mesh_size_variable_order);
	      }
	    else
	      {
		int n = parameters.GetM();
		mesh_size_variable_order.Reallocate(n);
		mesh_size_variable_order(0) = 0;
		mesh_size_variable_order(1) = 0;
                T coef = to_num<T>(parameters(1));
		for (int i = 1; i < parameters.GetM()-1; i++)
		  mesh_size_variable_order(i+1) = to_num<Real_wp>(parameters(i+1));
                
                Mlt(coef, mesh_size_variable_order);
	      }
	    
             if (!parameters(0).compare("MEAN_EDGE"))
               type_variable_order = MEAN_EDGE_ORDER;
             else
               type_variable_order = MAX_EDGE_ORDER;
	  }
	else if (parameters(0) == "REF")
	  {
	    type_variable_order = REF_ORDER;
	    order_for_each_ref.Reallocate(parameters.GetM());
	    order_for_each_ref.Zero();
	    for (int i = 1; i < parameters.GetM(); i++)
	      order_for_each_ref(i) = to_num<int>(parameters(i));

	    Vector<int> order_ref(order_for_each_ref.GetM());
	    for (int i = 0; i < order_ref.GetM(); i++)
	      order_ref(i) = order_for_each_ref(i);
	  }
	else
	  {
	    int r = to_num<int>(parameters(0));
	    this->SetOrder(r);
            
            if (parameters.GetM() > 1)
              {
                int rgeom = r;
                bool quad_mesh = false;
                if (parameters(1) == "QUAD")
                  {
                    quad_mesh = true;
                    if (parameters.GetM() > 2)
                      rgeom = to_num<int>(parameters(2));
                  }
                else
                  rgeom = to_num<int>(parameters(1));
                
                this->SetGeometryOrder(rgeom, quad_mesh);
              }
            else
              {
                // by default same order for geometry
                this->SetGeometryOrder(r);
              }
	  }
      }
    else if (description_field == "OrderStrategy")
      {
	if (parameters(0) == "Max")
	  use_max_velocity_order = true;
	else if (parameters(0) == "Min")
	  use_max_velocity_order = false;
	else
	  {
	    cout << "Unknown strategy" << endl;
	    abort();
	  }
      }
    else if (!description_field.compare("UseSameDofsForPeriodicCondition"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MeshNumbering_Base" << endl;
	    cout << "UseSameDofsForPeriodicCondition needs more parameters, for instance :" 
                 << endl;
	    cout << "UseSameDofsForPeriodicCondition = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	this->type_formulation_periodic = STRONG_PERIODIC;
	if (!parameters(0).compare("YES"))
	  this->type_formulation_periodic = SAME_PERIODIC_DOFS;
        else if (!parameters(0).compare("WEAK"))
	  this->type_formulation_periodic = WEAK_PERIODIC;
      }
    else if (description_field == "PlaneWaveQuasiPeriodic")
      {
        if (parameters(0) == "YES")
          plane_wave_quasi_periodic = true;
        else
          plane_wave_quasi_periodic = false;
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

	if (parameters.GetM() > 1)
	  {
	    bool quad_mesh = false;
	    if (parameters(1) == "QUAD")
	      quad_mesh = true;
	    
	    this->SetGeometryOrder(r, quad_mesh);
	  }
	else
	  this->SetGeometryOrder(r);
      }

  }


    //! dofs arrays are cleared
  template<class T>
  void MeshNumbering_Base<T>::ClearDofs()
  {
    nodl = 0;
    for (int i = 0; i < this->dof_element.GetM(); i++)
      this->dof_element(i).Clear();
    
    OffsetDofFaceNumber.Clear();
    OffsetDofEdgeNumber.Clear();
    OffsetDofElementNumber.Clear();
    OffsetDofVertexNumber.Clear();
    OffsetQuadElementNumber.Clear();

    PeriodicDofs_Original.Clear();
    type_periodicity_dof.Clear();
  }

  
  //! clears mesh 
  template<class T>
  void MeshNumbering_Base<T>::Clear()
  {
    nodl = 0; nb_dof_element.Clear();
    OffsetDofFaceNumber.Clear();
    OffsetDofEdgeNumber.Clear();
    OffsetDofVertexNumber.Clear();
    OffsetDofElementNumber.Clear();
    order_element.Clear(); order_inside.Clear();
    order_edge.Clear(); order_face.Clear();
    order_boundary.Clear();
    nb_points_quadrature_on_boundary.Clear();
    
  }
  

  //! numbering of the mesh
  template<class T>
  void MeshNumbering_Base<T>::NumberMesh(bool surface_mesh)
  {
    // previous dofs are removed
    this->ClearDofs();
    if (surface_mesh)
      this->PerformSurfaceMeshNumbering();
    else
      this->PerformMeshNumbering();
  }


  /*********************
   * MeshNumbering_Dim *
   *********************/
  
  
  //! default constructor
  template<class Dimension, class T>
  MeshNumbering_Dim<Dimension, T>::MeshNumbering_Dim(Mesh<Dimension, T>& m)
    : MeshNumbering_Base<T>(), mesh(m)
  {
  }


  //! copy constructor
  template<class Dimension, class T>
  MeshNumbering_Dim<Dimension, T>::MeshNumbering_Dim(const MeshNumbering_Dim<Dimension, T>& m)
    : MeshNumbering_Base<T>(m), mesh(m.mesh), TranslationPeriodicBoundary(m.TranslationPeriodicBoundary),
      TranslationPeriodicDof(m.TranslationPeriodicDof)
  {
  }
  

  //! operator = with another class of the same type
  template<class Dimension, class T>
  MeshNumbering_Dim<Dimension, T>& MeshNumbering_Dim<Dimension, T>
  ::operator=(const MeshNumbering_Dim<Dimension, T>& m)
  {
    static_cast<MeshNumbering_Base<T>& >(*this)
      = static_cast<const MeshNumbering_Base<T>& >(m);
    
    TranslationPeriodicBoundary = m.TranslationPeriodicBoundary;
    TranslationPeriodicDof = m.TranslationPeriodicDof;
        
    return *this;
  }
  

  //! allocates n elements (each element contains degrees of freedom), previous dofs are lost
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::ReallocateElements(int n)
  {
    this->dof_element.Reallocate(n);
    if (mesh.GetNbElt() > 0)
      for (int i = 0; i < n; i++)
        this->dof_element(i).ReallocateFaces(mesh.Element(i).GetNbBoundary());
  }
  

  //! informs that there are n periodic dofs
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::ReallocatePeriodicDof(int n)
  {
    this->type_periodicity_dof.Reallocate(n);
    this->PeriodicDofs.Reallocate(n);
    this->PeriodicDofs_Original.Reallocate(n);
    this->TranslationPeriodicDof.Reallocate(n);
  }


  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::ClearDofs()
  {
    MeshNumbering_Base<T>::ClearDofs();
    TranslationPeriodicDof.Clear();
  }


  //! returns size of memory used by the object in bytes
  template<class Dimension, class T>
  size_t MeshNumbering_Dim<Dimension, T>::GetMemorySize() const
  {
    size_t taille = MeshNumbering_Base<T>::GetMemorySize();

    taille += TranslationPeriodicBoundary.GetMemorySize();
    taille += TranslationPeriodicDof.GetMemorySize();
    
    return taille;
  }
  
  
  //! returns an array containing all the orders found in the mesh
  /*!
    \param[out] num num(0) will contain the orders of approximation used for tetrahedra
                    num(1) for pyramids, num(2) for prisms and num(3) for hexahedra
		    For 2-D meshes, num(0) will be associated with triangles and num(1)
		    with quadrilaterals
   */
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::GetOrder(TinyVector<IVect, 4>& num) const
  {
    if (this->type_variable_order != this->CONSTANT_ORDER)
      {
	char rmax = 0;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  rmax = max(rmax, this->order_element(i));
        
        for (int i = 0; i < mesh.GetNbEdges(); i++)
          rmax = max(rmax, this->order_edge(i));
        
	
        if (Dimension::dim_N == 3)
          for (int i = 0; i < mesh.GetNbBoundary(); i++)
            rmax = max(rmax, this->order_face(i));
        
	Vector<bool> OrderUsed_Tet(rmax+1), OrderUsed_Pyr(rmax+1);
	Vector<bool> OrderUsed_Wed(rmax+1), OrderUsed_Hex(rmax+1);
	OrderUsed_Tet.Fill(false); OrderUsed_Pyr.Fill(false);
	OrderUsed_Wed.Fill(false); OrderUsed_Hex.Fill(false);
        Vector<char> re(12), rf(6);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int type = mesh.Element(i).GetHybridType();
	    char r = this->order_element(i);
            char ri = this->order_inside(i);
            int nb_edges = mesh.Element(i).GetNbEdges();
            for (int j = 0; j < nb_edges; j++)
              re(j) = this->order_edge(mesh.Element(i).numEdge(j));

            int nb_faces = 0;
            if (Dimension::dim_N == 3)
              {
                nb_faces = mesh.Element(i).GetNbBoundary();
                for (int j = 0; j < nb_faces; j++)
                  rf(j) = this->order_face(mesh.Element(i).numBoundary(j));
              }
            
            switch (type)
              {
              case 0 :
                {
                  OrderUsed_Tet(r) = true;
                  OrderUsed_Tet(ri) = true;
                  for (int j = 0; j < nb_edges; j++)
                    OrderUsed_Tet(re(j)) = true;

                  for (int j = 0; j < nb_faces; j++)
                    OrderUsed_Tet(rf(j)) = true;
                }
                break;
              case 1 :
                {
                  OrderUsed_Pyr(r) = true;
                  OrderUsed_Pyr(ri) = true;
                  for (int j = 0; j < nb_edges; j++)
                    OrderUsed_Pyr(re(j)) = true;

                  for (int j = 0; j < nb_faces; j++)
                    OrderUsed_Pyr(rf(j)) = true;
                }
                break;
              case 2 :
                {
                  OrderUsed_Wed(r) = true;
                  OrderUsed_Wed(ri) = true;
                  for (int j = 0; j < nb_edges; j++)
                    OrderUsed_Wed(re(j)) = true;

                  for (int j = 0; j < nb_faces; j++)
                    OrderUsed_Wed(rf(j)) = true;
                }
                break;
              case 3 :
                {
                  OrderUsed_Hex(r) = true;
                  OrderUsed_Hex(ri) = true;
                  for (int j = 0; j < nb_edges; j++)
                    OrderUsed_Hex(re(j)) = true;

                  for (int j = 0; j < nb_faces; j++)
                    OrderUsed_Hex(rf(j)) = true;
                }
                break;
              }
	  }
	
	num(0).Clear(); num(1).Clear();
	num(2).Clear(); num(3).Clear();
	for (char r = 1; r <= rmax; r++)
	  {
	    if (OrderUsed_Tet(r))
	      num(0).PushBack(r);
	    
	    if (OrderUsed_Pyr(r))
	      num(1).PushBack(r);
	    
	    if (OrderUsed_Wed(r))
	      num(2).PushBack(r);
	    
	    if (OrderUsed_Hex(r))
	      num(3).PushBack(r);
	  }
      }
    else
      {
        int n0 = 0, n1 = 0, n2 = 0, n3 = 0;
        for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            int type = mesh.Element(i).GetHybridType();
            switch (type)
              {
              case 0 : n0++; break;
              case 1 : n1++; break;
              case 2 : n2++; break;
              case 3 : n3++; break;
              }
          }
        
        if (n0 > 0)
          {
            num(0).Reallocate(1);
            num(0)(0) = this->order;
          }

        if (n1 > 0)
          {
            num(1).Reallocate(1);
            num(1)(0) = this->order;
          }

        if (n2 > 0)
          {
            num(2).Reallocate(1);
            num(2)(0) = this->order;
          }
        
        if (n3 > 0)
          {
            num(3).Reallocate(1);
            num(3)(0) = this->order;
          }        
      }
  }


  //! changes order of approximation for geometry
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::SetGeometryOrder(int r, bool only_quad)
  {
    mesh.SetGeometryOrder(r, only_quad);
  }
  
  
  //! computation of the order associated with each edge, face, element of the mesh
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::ComputeVariableOrder(int dg_form)
  {
    if (this->type_variable_order == this->USER_ORDER)
      return;

    if (this->type_variable_order == this->REF_ORDER)
      {
	// order for each element is given through its reference
	this->order_inside.Reallocate(mesh.GetNbElt());
	this->order_element.Reallocate(mesh.GetNbElt());
        this->order = 1;
	this->order_edge.Reallocate(mesh.GetNbEdges());
	this->order_edge.Zero();
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = mesh.Element(i).GetReference();
	    char r = this->order_for_each_ref(ref);
	    this->order_element(i) = r;
	    this->order_inside(i) = r;
	    this->order = max(r, this->order);
	    for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	      {
		int ne = mesh.Element(i).numEdge(j);
		this->order_edge(ne) = max(this->order_edge(ne), r);
	      }	    
	  }
      }    
    else if (this->type_variable_order != this->CONSTANT_ORDER)
      {	
	if (this->coef_elt_variable_order.GetM() != mesh.GetNbElt())
	  {
	    this->coef_elt_variable_order.Reallocate(mesh.GetNbElt());
	    this->coef_elt_variable_order.Fill(1.0);
	  }
	
	this->order_edge.Reallocate(mesh.GetNbEdges());
	
	// we compute order for each edge
	for (int i = 0; i < mesh.GetNbEdges(); i++)
	  {
	    int n1 = mesh.GetEdge(i).numVertex(0);
	    int n2 = mesh.GetEdge(i).numVertex(1);
	    Real_wp length = mesh.Vertex(n1).Distance(mesh.Vertex(n2));
	    Real_wp coef = 1;
	    for (int j = 0; j < mesh.GetEdge(i).GetNbElements(); j++)
	      {
		Real_wp coef1 = this->coef_elt_variable_order(mesh.GetEdge(i).numElement(j));
		if (j == 0)
		  coef = coef1;
		else
		  {
		    if (this->use_max_velocity_order)
		      coef = max(coef, coef1);
		    else
		      coef = min(coef, coef1);
		  }
	      }
            
	    length *= coef;
	    int r = this->mesh_size_variable_order.GetM()-1;
	    for (int k = 1; k < this->mesh_size_variable_order.GetM(); k++) 
	      if (length >= this->mesh_size_variable_order(k))
		r = k;
            
	    this->order_edge(i) = r;
	  }	
	
      }
    
    if (this->type_variable_order == this->MAX_EDGE_ORDER)
      {
	// now elements
	this->order_inside.Reallocate(mesh.GetNbElt());
	this->order_element.Reallocate(mesh.GetNbElt());
        this->order = 1;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    char r = 1;
	    for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	      r = max(r, this->order_edge(mesh.Element(i).numEdge(j)));
	    
	    this->order_inside(i) = r;
	    this->order_element(i) = r;
            this->order = max(this->order, r);
	  }
        
      }
    else if (this->type_variable_order == this->MEAN_EDGE_ORDER)
      {	
	// now elements
	this->order_inside.Reallocate(mesh.GetNbElt());
	this->order_element.Reallocate(mesh.GetNbElt());
        this->order = 1;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    char r = 1;
            int nb_edges = mesh.Element(i).GetNbEdges();
	    for (int j = 0; j < nb_edges; j++)
	      r += this->order_edge(mesh.Element(i).numEdge(j));
	    
	    this->order_inside(i) = r/char(nb_edges);
	    this->order_element(i) = r/char(nb_edges);
            this->order = max(this->order, this->order_element(i));
	  }
      }
    
    // quadrature rule used for faces is the maximum
    // of order for the two adjoining elements to have exact integration
    if (this->type_variable_order != this->CONSTANT_ORDER)
      {
        this->order_boundary.Reallocate(mesh.GetNbBoundary());
        this->order_boundary.Fill(this->order);
        for (int i = 0; i < mesh.GetNbBoundary(); i++)
          {
            int num_elem = mesh.Boundary(i).numElement(0);
            if (num_elem >= 0)
              this->order_boundary(i) = this->order_element(num_elem);
               
            num_elem = mesh.Boundary(i).numElement(1); 
            if (num_elem >= 0)
              this->order_boundary(i)
                = max(this->order_element(num_elem), this->order_boundary(i));
          }
      }

  }
  
  
  //! numbering of the mesh
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::PerformMeshNumbering()
  {
    MeshNumbering<Dimension, T>& mesh_leaf = this->GetLeafClass();
    
    // retrieving informations about mesh
    int nb_vertices = mesh.GetNbVertices();
    int nb_edges = mesh.GetNbEdges();
    int nb_elt = mesh.GetNbElt();
    IVect NegativeDof(5000); int nb_negative_dof = 0;
    
    bool variable_order = this->IsOrderVariable();
    this->nb_dof_element.Reallocate(nb_elt);
    this->ReallocateElements(nb_elt);
    
    if (this->number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	// allocating dof numbers
	for (int i = 0; i < nb_elt; i++)
	  {
	    this->nb_dof_element(i) = mesh_leaf.GetNbDofElement(i, mesh.Element(i));
	    this->dof_element(i).ReallocateDof(this->nb_dof_element(i));
	    for (int j = 0; j < this->dof_element(i).GetNbDof(); j++)
	      this->dof_element(i).SetNumberDof(j, -1);
	  }
	
	// Degrees of freedom (Dofs) are numbered in the following order :
	// 1- dofs related to vertices
	// 2- dofs related to edges
	// 3- dofs related to faces (in 3-D)
	// 4- dofs related to elements
	
	// we init the total number of dofs
	this->nodl = 0;
	
	// First step : dofs on vertices
	if (!variable_order)
	  {
	    int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
            this->OffsetDofVertexNumber.Reallocate(mesh.GetNbVertices()+1);
	    this->OffsetDofVertexNumber.Fill();
	    Mlt(nb_dof_vertex, this->OffsetDofVertexNumber);
	    
	    for (int i = 0; i < nb_elt; i++)
	      {
		int nb_vertices_loc = mesh.Element(i).GetNbVertices();
		int offset_loc = 0;
		for (int j = 0; j < nb_vertices_loc; j++)
		  {
		    int nv = mesh.Element(i).numVertex(j);
		    for (int k = 0; k < nb_dof_vertex; k++)
		      this->dof_element(i).SetNumberDof(offset_loc + k, nv*nb_dof_vertex + k);
		    
		    offset_loc += nb_dof_vertex;
		  }
	      }
	    
	    this->nodl += nb_dof_vertex*nb_vertices;
	  }
	else
	  {
	    IVect offset_vertex(nb_vertices+1);
	    offset_vertex(0) = 0;

	    for (int i = 0; i < nb_vertices; i++)
              offset_vertex(i+1) = offset_vertex(i) + this->number_map.GetNbDofVertex(0);
            
            this->OffsetDofVertexNumber = offset_vertex;
	    
	    int nb_dof_vertex = this->number_map.GetNbDofVertex(0);
	    for (int i = 0; i < nb_elt; i++)
	      {
		int nb_vertices_loc = mesh.Element(i).GetNbVertices();		
		for (int j = 0; j < nb_vertices_loc; j++)
		  {
		    int nv = mesh.Element(i).numVertex(j);
		    for (int k = 0; k < nb_dof_vertex; k++)
		      this->Element(i).SetNumberDof(j*nb_dof_vertex + k,
						    offset_vertex(nv) + k);
		  }
	      }
	    
	    this->nodl += offset_vertex(nb_vertices);
	  }
	
	// Third step : dofs on edges
	if (!variable_order)
	  {
	    int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
	    int nb_dof_edge = this->number_map.GetNbDofEdge(this->order);
	    this->OffsetDofEdgeNumber.Reallocate(mesh.GetNbEdges()+1);
	    this->OffsetDofEdgeNumber.Fill();
	    Mlt(nb_dof_edge, this->OffsetDofEdgeNumber);
            for (int i = 0; i <= mesh.GetNbEdges(); i++)
              this->OffsetDofEdgeNumber(i) += this->nodl;
              
	    for (int i = 0; i < nb_elt; i++)
	      {
                int type_elt = mesh.Element(i).GetHybridType();
		int nb_vertices_loc = mesh.Element(i).GetNbVertices();
		int nb_edges_loc = mesh.Element(i).GetNbEdges();
		int offset_loc = nb_vertices_loc*nb_dof_vertex;
		nb_negative_dof = 0;
		for (int j = 0; j < nb_edges_loc; j++)
		  {
		    int ne = mesh.Element(i).numEdge(j);
		    
		    int ne_loc = mesh_leaf.FirstExtremityEdge(type_elt, j);
		    int num_vertex = mesh.Element(i).numVertex(ne_loc);
		    
		    if (mesh.GetEdge(ne).numVertex(0) == num_vertex)
		      {
			// same orientation between local and global edge
			for (int k = 0; k < nb_dof_edge; k++)
			  this->dof_element(i).SetNumberDof(offset_loc + k,
						      this->nodl + ne*nb_dof_edge + k);
		      }
		    else
		      {
			// opposite orientation between local and global edge
			for (int k = 0; k < nb_dof_edge; k++)
			  {
			    int num_dof_loc
                              = offset_loc + this->number_map.GetSymmetricEdgeDof(this->order, k);
			    
                            this->dof_element(i).SetNumberDof(num_dof_loc, this->nodl + ne*nb_dof_edge + k);
			    
			    if (this->number_map.IsSkewSymmetricEdgeDof(this->order, k))
			      {
				NegativeDof(nb_negative_dof) = num_dof_loc;
				nb_negative_dof++;
			      }
			  }
		      }
		    
		    offset_loc += nb_dof_edge;
		  }
		
		// if negative dofs, specifying it
		if (nb_negative_dof > 0)
		  this->dof_element(i).SetNegativeDofNumber(nb_negative_dof, NegativeDof);
	      }
	    
	    this->nodl += nb_dof_edge*nb_edges;
	  }
	else
	  {
	    IVect offset_edge(nb_edges+1);
	    offset_edge(0) = 0;
	    for (int i = 0; i < nb_edges; i++)
	      offset_edge(i+1) = offset_edge(i) + this->number_map.GetNbDofEdge(this->order_edge(i));
	    
	    this->OffsetDofEdgeNumber = offset_edge;
            for (int i = 0; i <= nb_edges; i++)
              this->OffsetDofEdgeNumber(i) += this->nodl;
            
	    for (int i = 0; i < nb_elt; i++)
	      {
                int type_elt = mesh.Element(i).GetHybridType();
		int nb_vertices_loc = mesh.Element(i).GetNbVertices();
		int nb_edges_loc = mesh.Element(i).GetNbEdges();
		int r = this->order_element(i);
		int offset_loc = nb_vertices_loc*this->number_map.GetNbDofVertex(r);
		int nb_dofE = this->number_map.GetNbDofEdge(r);
		
		nb_negative_dof = 0;
		for (int j = 0; j < nb_edges_loc; j++)
		  {
		    int ne = mesh.Element(i).numEdge(j);
		    int nb_dof_edge = this->number_map.GetNbDofEdge(this->order_edge(ne));
		    
		    int ne_loc = mesh_leaf.FirstExtremityEdge(type_elt, j);
		    int num_vertex = mesh.Element(i).numVertex(ne_loc);
		    
		    if (mesh.GetEdge(ne).numVertex(0) == num_vertex)
		      {
			// same orientation between local and global edge
			for (int k = 0; k < nb_dof_edge; k++)
			  this->dof_element(i).SetNumberDof(offset_loc + k,
							    this->nodl + offset_edge(ne) + k);
		      }
		    else
		      {
			// opposite orientation between local and global edge
			for (int k = 0; k < nb_dof_edge; k++)
			  {
			    int num_dof_loc = offset_loc + this->number_map.
			      GetSymmetricEdgeDof(this->order_edge(ne), k);
			    
			    this->dof_element(i).SetNumberDof(num_dof_loc,
                                                              this->nodl + offset_edge(ne) + k);
			    
			    if (this->number_map.IsSkewSymmetricEdgeDof(this->order_edge(ne), k))
			      {
				NegativeDof(nb_negative_dof) = num_dof_loc;
				nb_negative_dof++;
			      }
			  }
		      }
		    
		    offset_loc += nb_dofE;
		  }
		
		// if negative dofs, specifying it
		if (nb_negative_dof > 0)
		  this->dof_element(i).SetNegativeDofNumber(nb_negative_dof, NegativeDof);
		
	      }
	    
	    this->nodl += offset_edge(nb_edges);
	  }
		
	// Fourth step : dofs on faces
	mesh_leaf.NumberFaces();
	
	// Last step : dofs on elements
	for (int i = 0; i < nb_elt; i++)
	  {
	    int nb_dof_elt = 0, offset_loc = 0;
	    if (!variable_order)
	      {
		nb_dof_elt = this->number_map.GetNbDofElement(this->order, mesh.Element(i));
		offset_loc = this->dof_element(i).GetNbDof() - nb_dof_elt;
	      }
	    else
	      {
		nb_dof_elt = this->number_map.GetNbDofElement(this->order_inside(i), mesh.Element(i));
		offset_loc = this->dof_element(i).GetNbDof()
                  - this->number_map.GetNbDofElement(this->order_element(i), mesh.Element(i));
	      }
	    
	    for (int k = 0; k < nb_dof_elt; k++)
	      this->dof_element(i).SetNumberDof(offset_loc + k, this->nodl + k);
	    
	    this->nodl += nb_dof_elt;
	  }
      }
    else if (this->number_map.FormulationDG() == ElementReference_Base::HDG)
      {
	// counting dofs on each boundary (edge in 2-D, face in 3-D)
	this->nodl = 0;
	this->OffsetDofFaceNumber.Reallocate(mesh.GetNbBoundary()+1);
	this->OffsetDofFaceNumber(0) = 0;
	for (int i = 0; i < mesh.GetNbBoundary(); i++)
	  {
	    int r = mesh_leaf.GetOrderQuadrature(i);
	    int nb_dof = this->number_map.GetNbDofBoundary(r, mesh.Boundary(i));
	    this->nodl += nb_dof;
	    this->OffsetDofFaceNumber(i+1) = this->nodl;
	  }
	
	// then affecting dofs to elements
	for (int i = 0; i < nb_elt; i++)
	  {
	    int nb_dof = 0;
	    for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
	      {
		int nf = mesh.Element(i).numBoundary(j);
		nb_dof += this->OffsetDofFaceNumber(nf+1) - this->OffsetDofFaceNumber(nf);
	      }

	    this->nb_dof_element(i) = nb_dof;
	    this->dof_element(i).ReallocateDof(nb_dof);

	    int offset_loc = 0;
	    nb_negative_dof = 0;
	    for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
	      {
		int nf = mesh.Element(i).numBoundary(j);
		int rf = mesh_leaf.GetOrderQuadrature(nf);
		int nb_dof_boundary = this->OffsetDofFaceNumber(nf+1) - this->OffsetDofFaceNumber(nf);

		this->PerformLocalNumbering(i, j, nb_dof_boundary, offset_loc, nf, rf,
					    nb_negative_dof, NegativeDof);
		
		offset_loc += nb_dof_boundary;
	      }

	    // if negative dofs, specifying it
	    if (nb_negative_dof > 0)
	      this->dof_element(i).SetNegativeDofNumber(nb_negative_dof, NegativeDof);
	    
	  }
      }
    else
      {
	// dg formulation -> only dofs on elements
	
	// discontinuous galerkin, no number is stored
	// we store only an "offset number"
	// so that the global number is easily retrieved 
	// glob_number = local_number + offset
		
	this->nodl = 0;
	this->OffsetDofElementNumber.Reallocate(nb_elt+1);
	this->OffsetDofElementNumber(0) = 0;
	for (int i = 0; i < nb_elt; i++)
	  {
	    this->dof_element(i).ReallocateDof(1);
	    this->dof_element(i).SetNumberDof(0, this->nodl);
	    
	    int order = this->GetOrderElement(i);
	    int nb_dof_elt = this->number_map.GetNbDofElement(order, mesh.Element(i));
	    this->nb_dof_element(i) = nb_dof_elt;
	    this->nodl += nb_dof_elt;
            this->OffsetDofElementNumber(i+1) = this->nodl;
	  }

      }
    
    // offset for integration on boundaries
    this->nb_points_quadrature_on_boundary.Reallocate(mesh.GetNbBoundary());
    this->nb_points_quadrature_on_boundary.Fill(0);
    
    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      {
        int r = mesh_leaf.GetOrderQuadrature(i);
        int nb_points = this->number_map.GetNbPointsQuadBoundary(r, mesh.Boundary(i));;
        this->nb_points_quadrature_on_boundary(i) = nb_points;
      }
    
    this->OffsetQuadElementNumber.Reallocate(nb_elt+1);
    this->OffsetQuadElementNumber(0) = 0;
    for (int i = 0; i < nb_elt; i++)
      {
        int nb_points_faces = 0;
        for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
          nb_points_faces += this->nb_points_quadrature_on_boundary(mesh.Element(i).numBoundary(j));
	
        this->OffsetQuadElementNumber(i+1) = this->OffsetQuadElementNumber(i) + nb_points_faces;
      }
    
    // dofs can be modified because of periodic boundary conditions
    if (this->treat_periodic_condition_during_number)
      mesh_leaf.TreatPeriodicCondition();

    // dofs for pml may be numbered as well
    if (this->compute_dof_pml)
      {
        if (this->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          {
            Vector<bool> DofPML(this->nodl);
            if (this->drop_interface_pml_dof)
              {
                DofPML.Fill(true);
                for (int i = 0; i < mesh.GetNbElt(); i++)
                  if (!mesh.Element(i).IsPML())
                    {
                      for (int j = 0; j < this->dof_element(i).GetNbDof(); j++)
                        {
                          int num_dof = this->dof_element(i).GetNumberDof(j);
                          DofPML(num_dof) = false;
                        }
                    }
              }
            else
              {
                DofPML.Fill(false);
                for (int i = 0; i < mesh.GetNbElt(); i++)
                  if (mesh.Element(i).IsPML())
                    {
                      for (int j = 0; j < this->dof_element(i).GetNbDof(); j++)
                        {
                          int num_dof = this->dof_element(i).GetNumberDof(j);
                          DofPML(num_dof) = true;
                        }
                    }
              }
            
            this->nb_dof_pml = 0;
            this->IndexDofPML.Reallocate(this->nodl);
            this->IndexDofPML.Fill(-1);
            for (int i = 0; i < this->nodl; i++)
              if (DofPML(i))
                this->IndexDofPML(i) = this->nb_dof_pml++;
          }
        else
          {
            this->nb_dof_pml = 0;
            this->IndexDofPML.Reallocate(this->nodl);
            this->IndexDofPML.Fill(-1);
            for (int i = 0; i < mesh.GetNbElt(); i++)
              if (mesh.Element(i).IsPML())
                {
                  int nb_dof = this->OffsetDofElementNumber(i+1) - this->OffsetDofElementNumber(i);
                  for (int j = 0; j < nb_dof; j++)
                    this->IndexDofPML(this->OffsetDofElementNumber(i)+j) = this->nb_dof_pml + j;
                  
                  this->nb_dof_pml += nb_dof;
                }
          }
      }    
  }



  //! numbering of the mesh
  template<class Dimension, class T>
  void MeshNumbering_Dim<Dimension, T>::PerformSurfaceMeshNumbering()
  {
    MeshNumbering<Dimension, T>& mesh_leaf = this->GetLeafClass();
    
    // retrieving informations about mesh
    int nb_vertices = mesh.GetNbVertices();
    int nb_elt = mesh.GetNbBoundaryRef();
    
    bool variable_order = this->IsOrderVariable();
    this->nb_dof_element.Reallocate(nb_elt);
    this->ReallocateElements(nb_elt);

    if (this->number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	// allocating dof numbers
	for (int i = 0; i < nb_elt; i++)
	  {
	    this->nb_dof_element(i) = mesh_leaf.GetNbDofElement(i, mesh.BoundaryRef(i));
	    this->dof_element(i).ReallocateDof(this->nb_dof_element(i));
	    for (int j = 0; j < this->dof_element(i).GetNbDof(); j++)
	      this->dof_element(i).SetNumberDof(j, -1);
	  }
	
	// Degrees of freedom (Dofs) are numbered in the following order :
	// 1- dofs related to vertices
	// 2- dofs related to edges
	// 3- dofs related to faces (in 3-D)
	
	// we init the total number of dofs
	this->nodl = 0;
	
	// First step : dofs on vertices
	if (!variable_order)
	  {
	    int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
            this->OffsetDofVertexNumber.Reallocate(mesh.GetNbVertices()+1);
	    this->OffsetDofVertexNumber.Fill();
	    Mlt(nb_dof_vertex, this->OffsetDofVertexNumber);
	    
	    for (int i = 0; i < nb_elt; i++)
	      {
		int nb_vertices_loc = mesh.BoundaryRef(i).GetNbVertices();
		int offset_loc = 0;
		for (int j = 0; j < nb_vertices_loc; j++)
		  {
		    int nv = mesh.BoundaryRef(i).numVertex(j);
		    for (int k = 0; k < nb_dof_vertex; k++)
		      this->dof_element(i).SetNumberDof(offset_loc + k, nv*nb_dof_vertex + k);
		    
		    offset_loc += nb_dof_vertex;
		  }
	      }
	    
	    this->nodl += nb_dof_vertex*nb_vertices;
	  }
	else
	  {
	    IVect offset_vertex(nb_vertices+1);
	    offset_vertex(0) = 0;

	    for (int i = 0; i < nb_vertices; i++)
              offset_vertex(i+1) = offset_vertex(i) + this->number_map.GetNbDofVertex(0);
            
            this->OffsetDofVertexNumber = offset_vertex;
	    
	    int nb_dof_vertex = this->number_map.GetNbDofVertex(0);
	    for (int i = 0; i < nb_elt; i++)
	      {
		int nb_vertices_loc = mesh.BoundaryRef(i).GetNbVertices();		
		for (int j = 0; j < nb_vertices_loc; j++)
		  {
		    int nv = mesh.BoundaryRef(i).numVertex(j);
		    for (int k = 0; k < nb_dof_vertex; k++)
		      this->Element(i).SetNumberDof(j*nb_dof_vertex + k,
						    offset_vertex(nv) + k);
		  }
	      }
	    
	    this->nodl += offset_vertex(nb_vertices);
	  }
			
	// Third step : dofs on edges (in 3-D)
	mesh_leaf.NumberSurfaceEdges();
	
	// Last step : dofs on elements
	for (int i = 0; i < nb_elt; i++)
	  {
	    int nb_dof_elt = 0, offset_loc = 0;
	    if (!variable_order)
	      {
		nb_dof_elt = this->number_map.GetNbDofElement(this->order, mesh.BoundaryRef(i));
		offset_loc = this->dof_element(i).GetNbDof() - nb_dof_elt;
	      }
	    else
	      {
		nb_dof_elt = this->number_map.GetNbDofElement(this->order_inside(i), mesh.BoundaryRef(i));
		offset_loc = this->dof_element(i).GetNbDof()
                  - this->number_map.GetNbDofElement(this->order_element(i), mesh.BoundaryRef(i));
	      }
	    
	    for (int k = 0; k < nb_dof_elt; k++)
	      this->dof_element(i).SetNumberDof(offset_loc + k, this->nodl + k);
	    
	    this->nodl += nb_dof_elt;
	  }
      }
    else
      {
	cout << "not implemented " << endl;
	abort();
      }    
  }
  

#ifdef MONTJOIE_WITH_TWO_DIM
  template<>
  void MeshNumbering_Dim<Dimension2, Real_wp>
  ::PerformLocalNumbering(int i, int j, int nb_dof_boundary, int offset_loc, int nf, int rf,
			  int& nb_negative_dof, IVect& NegativeDof)    
  {
    // 2-D case
    int num_vertex = mesh.Element(i).numVertex(j);
    if (mesh.GetEdge(nf).numVertex(0) == num_vertex)
      {
	// same orientation between local and global edge
	for (int k = 0; k < nb_dof_boundary; k++)
	  this->dof_element(i).SetNumberDof(offset_loc + k,
					    this->OffsetDofFaceNumber(nf) + k);
      }
    else
      {
	// opposite orientation between local and global edge
	for (int k = 0; k < nb_dof_boundary; k++)
	  {
	    int num_dof_loc
	      = offset_loc + this->number_map.GetSymmetricEdgeDof(rf, k);
	    
	    this->dof_element(i).SetNumberDof(num_dof_loc, this->OffsetDofFaceNumber(nf) + k);
	    
	    if (this->number_map.IsSkewSymmetricEdgeDof(rf, k))
	      {
		NegativeDof(nb_negative_dof) = num_dof_loc;
		nb_negative_dof++;
	      }
	  }
      }
  }
#endif


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  void MeshNumbering_Dim<Dimension3, Real_wp>
  ::PerformLocalNumbering(int i, int j, int nb_dof_boundary, int offset_loc, int nf, int rf,
			  int& nb_negative_dof, IVect& NegativeDof)
  {
    // 3-D case
    if (this->number_map.DofInvariantByRotation(rf, this->mesh.Boundary(nf)))
      {
	int rot = this->mesh.Element(i).GetOrientationFace(j);
	// dofs invariation by rotation
	// we can find dof numbers obtained after rotation of the face
	for (int k = 0; k < nb_dof_boundary; k++)
	  {
	    int krot = this->number_map.
	      GetRotationFaceDof(rot, rf, k, this->mesh.Boundary(nf));
	    
	    int num_dof_loc = offset_loc + krot;
	    
	    this->Element(i).
	      SetNumberDof(num_dof_loc, this->OffsetDofFaceNumber(nf) + k);
	    
	    if (this->number_map.IsNegativeDof(rot, rf, k, this->mesh.Boundary(nf)))
	      {
		NegativeDof(nb_negative_dof) = num_dof_loc;
		nb_negative_dof++;
	      }
	  }
      }
    else
      {
	// dofs not invariant by rotation
	// local dofs are obtained by performing linear combination of global dofs
	// we only store global dofs (without any renumbering)
	for (int k = 0; k < nb_dof_boundary; k++)
	  {                    
	    int num_dof_loc = offset_loc + k;
	    
	    this->Element(i).
	      SetNumberDof(num_dof_loc, this->OffsetDofFaceNumber(nf) + k);
	    
	  }
	
	// and specifying an order for face
	this->Element(i).SetOrderFace(j, rf, offset_loc);
        
      }
  }
#endif
  
  
  /*****************************
   * MeshNumbering<Dimension2> *
   *****************************/
  
  
  //! initialization of static parameters
  template<class T>
  void InitStaticMeshData(const Dimension2& dim, const T& x)
  {
    // arrays containing the extremities of each edge of reference elements
    MeshNumbering<Dimension2, T>::FirstExtremityEdge.Fill(-1);
    MeshNumbering<Dimension2, T>::SecondExtremityEdge.Fill(-1);
    
    // edges of the square
    for (int i = 0; i < 4; i++)
      {
	MeshNumbering<Dimension2, T>::FirstExtremityEdge(1, i) = i;
	MeshNumbering<Dimension2, T>::SecondExtremityEdge(1, i) = (i+1)%4;
      }
    
    // edges of the triangle
    for (int i = 0; i < 3; i++)
      {
	MeshNumbering<Dimension2, T>::FirstExtremityEdge(0, i) = i;
	MeshNumbering<Dimension2, T>::SecondExtremityEdge(0, i) = (i+1)%3;
      }
    
    // entities for Gmsh
    NumberMap::NbNodesGmshEntity.Reallocate(120); NumberMap::NbNodesGmshEntity.Fill(-1);
    NumberMap::NbVerticesGmshEntity.Reallocate(120); NumberMap::NbVerticesGmshEntity.Fill(-1);
    NumberMap::DimensionGmshEntity.Reallocate(120); NumberMap::DimensionGmshEntity.Fill(-1);
    NumberMap::OrderGmshEntity.Reallocate(120); NumberMap::OrderGmshEntity.Fill(1);
    
    // 1 : 2-node line
    NumberMap::NbNodesGmshEntity(1) = 2; NumberMap::NbVerticesGmshEntity(1) = 2;
    NumberMap::DimensionGmshEntity(1) = 1; NumberMap::OrderGmshEntity(1) = 1;
    // 2 : 3-node triangle
    NumberMap::NbNodesGmshEntity(2) = 3; NumberMap::NbVerticesGmshEntity(2) = 3;
    NumberMap::DimensionGmshEntity(2) = 2; NumberMap::OrderGmshEntity(2) = 1;
    // 3 : 4-node quadrangle
    NumberMap::NbNodesGmshEntity(3) = 4; NumberMap::NbVerticesGmshEntity(3) = 4;
    NumberMap::DimensionGmshEntity(3) = 2; NumberMap::OrderGmshEntity(3) = 1;
    // 4 : 4-node tetrahedron
    NumberMap::NbNodesGmshEntity(4) = 4; NumberMap::NbVerticesGmshEntity(4) = 4;
    NumberMap::DimensionGmshEntity(4) = 3; NumberMap::OrderGmshEntity(4) = 1;
    // 5 : 8-node hexahedron
    NumberMap::NbNodesGmshEntity(5) = 8; NumberMap::NbVerticesGmshEntity(5) = 8;
    NumberMap::DimensionGmshEntity(5) = 3; NumberMap::OrderGmshEntity(5) = 1;
    // 6 : 6-node prism
    NumberMap::NbNodesGmshEntity(6) = 6; NumberMap::NbVerticesGmshEntity(6) = 6;
    NumberMap::DimensionGmshEntity(6) = 3; NumberMap::OrderGmshEntity(6) = 1;
    // 7 : 5-node pyramid
    NumberMap::NbNodesGmshEntity(7) = 5; NumberMap::NbVerticesGmshEntity(7) = 5;
    NumberMap::DimensionGmshEntity(7) = 3; NumberMap::OrderGmshEntity(7) = 1;
    // 8 : 3-node second order line
    NumberMap::NbNodesGmshEntity(8) = 3; NumberMap::NbVerticesGmshEntity(8) = 2;
    NumberMap::DimensionGmshEntity(8) = 1; NumberMap::OrderGmshEntity(8) = 2;
    // 9 : 6-node second order triangle
    NumberMap::NbNodesGmshEntity(9) = 6; NumberMap::NbVerticesGmshEntity(9) = 3;
    NumberMap::DimensionGmshEntity(9) = 2; NumberMap::OrderGmshEntity(9) = 2;
    // 10 : 9-node second order quadrangle
    NumberMap::NbNodesGmshEntity(10) = 9; NumberMap::NbVerticesGmshEntity(10) = 4;
    NumberMap::DimensionGmshEntity(10) = 2; NumberMap::OrderGmshEntity(10) = 2;
    // 11 : 10-node second order tetrahedron
    NumberMap::NbNodesGmshEntity(11) = 10; NumberMap::NbVerticesGmshEntity(11) = 4;
    NumberMap::DimensionGmshEntity(11) = 3; NumberMap::OrderGmshEntity(11) = 2;
    // 12 : 27-node second order hexahedron
    NumberMap::NbNodesGmshEntity(12) = 27; NumberMap::NbVerticesGmshEntity(12) = 8;
    NumberMap::DimensionGmshEntity(12) = 3; NumberMap::OrderGmshEntity(12) = 2;
    // 13 : 18-node second order prism
    NumberMap::NbNodesGmshEntity(13) = 18; NumberMap::NbVerticesGmshEntity(13) = 6;
    NumberMap::DimensionGmshEntity(13) = 3; NumberMap::OrderGmshEntity(18) = 2;
    // 14 : 14-node second order pyramid
    NumberMap::NbNodesGmshEntity(14) = 14; NumberMap::NbVerticesGmshEntity(14) = 5;
    NumberMap::DimensionGmshEntity(14) = 3; NumberMap::OrderGmshEntity(14) = 2;
    // 15 : 1-node point
    NumberMap::NbNodesGmshEntity(15) = 1; NumberMap::NbVerticesGmshEntity(15) = 1;
    NumberMap::DimensionGmshEntity(15) = 0; NumberMap::OrderGmshEntity(15) = 1;
    // 16 : 8-node second order quadrangle
    NumberMap::NbNodesGmshEntity(16) = 8; NumberMap::NbVerticesGmshEntity(16) = 4;
    NumberMap::DimensionGmshEntity(16) = 2; NumberMap::OrderGmshEntity(16) = 2;
    // 17 : 20-node second order hexahedron
    NumberMap::NbNodesGmshEntity(17) = 20; NumberMap::NbVerticesGmshEntity(17) = 8;
    NumberMap::DimensionGmshEntity(17) = 3; NumberMap::OrderGmshEntity(17) = 2;
    // 18 : 15-node second order prism
    NumberMap::NbNodesGmshEntity(18) = 15; NumberMap::NbVerticesGmshEntity(18) = 6;
    NumberMap::DimensionGmshEntity(18) = 3; NumberMap::OrderGmshEntity(18) = 2;
    // 19 : 13-node second order pyramid
    NumberMap::NbNodesGmshEntity(19) = 13; NumberMap::NbVerticesGmshEntity(19) = 5;
    NumberMap::DimensionGmshEntity(19) = 3; NumberMap::OrderGmshEntity(19) = 2;
    // 20 : 9-node third order triangle
    NumberMap::NbNodesGmshEntity(20) = 9; NumberMap::NbVerticesGmshEntity(20) = 3;
    NumberMap::DimensionGmshEntity(20) = 2; NumberMap::OrderGmshEntity(20) = 3;
    // 21 : 10-node third order triangle
    NumberMap::NbNodesGmshEntity(21) = 10; NumberMap::NbVerticesGmshEntity(21) = 3;
    NumberMap::DimensionGmshEntity(21) = 2; NumberMap::OrderGmshEntity(21) = 3;
    // 22 : 12-node fourth order triangle
    NumberMap::NbNodesGmshEntity(22) = 12; NumberMap::NbVerticesGmshEntity(22) = 3;
    NumberMap::DimensionGmshEntity(22) = 2; NumberMap::OrderGmshEntity(22) = 4;
    // 23 : 15-node fourth order triangle
    NumberMap::NbNodesGmshEntity(23) = 15; NumberMap::NbVerticesGmshEntity(23) = 3;
    NumberMap::DimensionGmshEntity(23) = 2; NumberMap::OrderGmshEntity(23) = 4;
    // 24 : 15-node fifth order triangle
    NumberMap::NbNodesGmshEntity(24) = 15; NumberMap::NbVerticesGmshEntity(24) = 3;
    NumberMap::DimensionGmshEntity(24) = 2; NumberMap::OrderGmshEntity(24) = 5;
    // 25 : 21-node fifth order triangle
    NumberMap::NbNodesGmshEntity(25) = 21; NumberMap::NbVerticesGmshEntity(25) = 3;
    NumberMap::DimensionGmshEntity(25) = 2; NumberMap::OrderGmshEntity(25) = 5;
    // 26 : 4-node third order line
    NumberMap::NbNodesGmshEntity(26) = 4; NumberMap::NbVerticesGmshEntity(26) = 2;
    NumberMap::DimensionGmshEntity(26) = 1; NumberMap::OrderGmshEntity(26) = 3;
    // 27 : 5-node fourth order line
    NumberMap::NbNodesGmshEntity(27) = 5; NumberMap::NbVerticesGmshEntity(27) = 2;
    NumberMap::DimensionGmshEntity(27) = 1; NumberMap::OrderGmshEntity(27) = 4;
    // 28 : 6-node fifth order line
    NumberMap::NbNodesGmshEntity(28) = 6; NumberMap::NbVerticesGmshEntity(28) = 2;
    NumberMap::DimensionGmshEntity(28) = 1; NumberMap::OrderGmshEntity(28) = 5;
    // 29 : 20-node third order tetrahedron
    NumberMap::NbNodesGmshEntity(29) = 20; NumberMap::NbVerticesGmshEntity(29) = 4;
    NumberMap::DimensionGmshEntity(29) = 3; NumberMap::OrderGmshEntity(29) = 3;
    // 30 : 35-node fourth order tetrahedron
    NumberMap::NbNodesGmshEntity(30) = 35; NumberMap::NbVerticesGmshEntity(30) = 4;
    NumberMap::DimensionGmshEntity(30) = 3; NumberMap::OrderGmshEntity(30) = 4;
    // 31 : 56-node fifth order tetrahedron
    NumberMap::NbNodesGmshEntity(31) = 56; NumberMap::NbVerticesGmshEntity(31) = 4;
    NumberMap::DimensionGmshEntity(31) = 3; NumberMap::OrderGmshEntity(31) = 5;
    // 32 : 34-node fourth order tetrahedron
    NumberMap::NbNodesGmshEntity(32) = 34; NumberMap::NbVerticesGmshEntity(32) = 4;
    NumberMap::DimensionGmshEntity(32) = 3; NumberMap::OrderGmshEntity(32) = 4;
    // 33 : 52-node fifth order tetrahedron
    NumberMap::NbNodesGmshEntity(33) = 52; NumberMap::NbVerticesGmshEntity(33) = 4;
    NumberMap::DimensionGmshEntity(33) = 3; NumberMap::OrderGmshEntity(33) = 5;


    // 36 : third order quadrangle
    NumberMap::NbNodesGmshEntity(36) = 16; NumberMap::NbVerticesGmshEntity(36) = 4;
    NumberMap::DimensionGmshEntity(36) = 2; NumberMap::OrderGmshEntity(36) = 3;
    // 37 : fourth order quadrangle
    NumberMap::NbNodesGmshEntity(37) = 25; NumberMap::NbVerticesGmshEntity(37) = 4;
    NumberMap::DimensionGmshEntity(37) = 2; NumberMap::OrderGmshEntity(37) = 4;
    // 38 : fifth order quadrangle
    NumberMap::NbNodesGmshEntity(38) = 36; NumberMap::NbVerticesGmshEntity(38) = 4;
    NumberMap::DimensionGmshEntity(38) = 2; NumberMap::OrderGmshEntity(38) = 5;
    // 39 : incomplete third order quadrangle
    NumberMap::NbNodesGmshEntity(39) = 12; NumberMap::NbVerticesGmshEntity(39) = 4;
    NumberMap::DimensionGmshEntity(39) = 2; NumberMap::OrderGmshEntity(39) = 3;
    // 40 : incomplete fourth order quadrangle
    NumberMap::NbNodesGmshEntity(40) = 16; NumberMap::NbVerticesGmshEntity(40) = 4;
    NumberMap::DimensionGmshEntity(40) = 2; NumberMap::OrderGmshEntity(40) = 4;
    // 41 : incomplete fifth order quadrangle
    NumberMap::NbNodesGmshEntity(41) = 20; NumberMap::NbVerticesGmshEntity(41) = 4;
    NumberMap::DimensionGmshEntity(41) = 2; NumberMap::OrderGmshEntity(41) = 5;
    // 42 : 28-node sixth order triangle
    NumberMap::NbNodesGmshEntity(42) = 28; NumberMap::NbVerticesGmshEntity(42) = 3;
    NumberMap::DimensionGmshEntity(42) = 2; NumberMap::OrderGmshEntity(42) = 6;
    // 43 : 36-node seventh order triangle
    NumberMap::NbNodesGmshEntity(43) = 36; NumberMap::NbVerticesGmshEntity(43) = 3;
    NumberMap::DimensionGmshEntity(43) = 2; NumberMap::OrderGmshEntity(43) = 7;
    // 44 : 45-node eigth order triangle
    NumberMap::NbNodesGmshEntity(44) = 45; NumberMap::NbVerticesGmshEntity(44) = 3;
    NumberMap::DimensionGmshEntity(44) = 2; NumberMap::OrderGmshEntity(44) = 8;
    // 45 : 55-node ninth order triangle
    NumberMap::NbNodesGmshEntity(45) = 55; NumberMap::NbVerticesGmshEntity(45) = 3;
    NumberMap::DimensionGmshEntity(45) = 2; NumberMap::OrderGmshEntity(45) = 9;
    // 46 : 66-node tenth order triangle
    NumberMap::NbNodesGmshEntity(46) = 66; NumberMap::NbVerticesGmshEntity(46) = 3;
    NumberMap::DimensionGmshEntity(46) = 2; NumberMap::OrderGmshEntity(46) = 10;    
    // 47 : sixth order quadrangle
    NumberMap::NbNodesGmshEntity(47) = 49; NumberMap::NbVerticesGmshEntity(47) = 4;
    NumberMap::DimensionGmshEntity(47) = 2; NumberMap::OrderGmshEntity(47) = 6;
    // 48 : seventh order quadrangle
    NumberMap::NbNodesGmshEntity(48) = 64; NumberMap::NbVerticesGmshEntity(48) = 4;
    NumberMap::DimensionGmshEntity(48) = 2; NumberMap::OrderGmshEntity(48) = 7;
    // 49 : eigth order quadrangle
    NumberMap::NbNodesGmshEntity(49) = 81; NumberMap::NbVerticesGmshEntity(49) = 4;
    NumberMap::DimensionGmshEntity(49) = 2; NumberMap::OrderGmshEntity(49) = 8;
    // 50 : ninth order quadrangle
    NumberMap::NbNodesGmshEntity(50) = 100; NumberMap::NbVerticesGmshEntity(50) = 4;
    NumberMap::DimensionGmshEntity(50) = 2; NumberMap::OrderGmshEntity(50) = 9;
    // 51 : tenth order quadrangle
    NumberMap::NbNodesGmshEntity(51) = 121; NumberMap::NbVerticesGmshEntity(51) = 4;
    NumberMap::DimensionGmshEntity(51) = 2; NumberMap::OrderGmshEntity(51) = 10;    
    // 52 : 18-node sixth order triangle
    NumberMap::NbNodesGmshEntity(52) = 18; NumberMap::NbVerticesGmshEntity(52) = 3;
    NumberMap::DimensionGmshEntity(52) = 2; NumberMap::OrderGmshEntity(52) = 6;
    // 53 : 21-node seventh order triangle
    NumberMap::NbNodesGmshEntity(53) = 21; NumberMap::NbVerticesGmshEntity(53) = 3;
    NumberMap::DimensionGmshEntity(53) = 2; NumberMap::OrderGmshEntity(53) = 7;
    // 54 : 24-node eigth order triangle
    NumberMap::NbNodesGmshEntity(54) = 24; NumberMap::NbVerticesGmshEntity(54) = 3;
    NumberMap::DimensionGmshEntity(54) = 2; NumberMap::OrderGmshEntity(54) = 8;
    // 55 : 27-node ninth order triangle
    NumberMap::NbNodesGmshEntity(55) = 27; NumberMap::NbVerticesGmshEntity(55) = 3;
    NumberMap::DimensionGmshEntity(55) = 2; NumberMap::OrderGmshEntity(55) = 9;
    // 56 : 30-node tenth order triangle
    NumberMap::NbNodesGmshEntity(56) = 30; NumberMap::NbVerticesGmshEntity(56) = 3;
    NumberMap::DimensionGmshEntity(56) = 2; NumberMap::OrderGmshEntity(56) = 10;    
    // 57 : incomplete sixth order quadrangle
    NumberMap::NbNodesGmshEntity(57) = 24; NumberMap::NbVerticesGmshEntity(47) = 4;
    NumberMap::DimensionGmshEntity(57) = 2; NumberMap::OrderGmshEntity(47) = 6;
    // 58 : incomplete seventh order quadrangle
    NumberMap::NbNodesGmshEntity(58) = 28; NumberMap::NbVerticesGmshEntity(48) = 4;
    NumberMap::DimensionGmshEntity(58) = 2; NumberMap::OrderGmshEntity(48) = 7;
    // 59 : incomplete eigth order quadrangle
    NumberMap::NbNodesGmshEntity(59) = 32; NumberMap::NbVerticesGmshEntity(49) = 4;
    NumberMap::DimensionGmshEntity(59) = 2; NumberMap::OrderGmshEntity(49) = 8;
    // 60 : incomplete ninth order quadrangle
    NumberMap::NbNodesGmshEntity(60) = 36; NumberMap::NbVerticesGmshEntity(50) = 4;
    NumberMap::DimensionGmshEntity(60) = 2; NumberMap::OrderGmshEntity(50) = 9;
    // 61 : incomplete tenth order quadrangle
    NumberMap::NbNodesGmshEntity(61) = 40; NumberMap::NbVerticesGmshEntity(51) = 4;
    NumberMap::DimensionGmshEntity(61) = 2; NumberMap::OrderGmshEntity(51) = 10;    
    // 62 : sixth-order edge
    NumberMap::NbNodesGmshEntity(62) = 7; NumberMap::NbVerticesGmshEntity(62) = 2;
    NumberMap::DimensionGmshEntity(62) = 1; NumberMap::OrderGmshEntity(62) = 6;    
    // 63 : seventh-order edge
    NumberMap::NbNodesGmshEntity(63) = 8; NumberMap::NbVerticesGmshEntity(63) = 2;
    NumberMap::DimensionGmshEntity(63) = 1; NumberMap::OrderGmshEntity(63) = 7;    
    // 64 : eighth-order edge
    NumberMap::NbNodesGmshEntity(64) = 9; NumberMap::NbVerticesGmshEntity(64) = 2;
    NumberMap::DimensionGmshEntity(64) = 1; NumberMap::OrderGmshEntity(64) = 8;    
    // 65 : ninth-order edge
    NumberMap::NbNodesGmshEntity(65) = 10; NumberMap::NbVerticesGmshEntity(65) = 2;
    NumberMap::DimensionGmshEntity(65) = 1; NumberMap::OrderGmshEntity(65) = 9;    
    // 66 : tenth-order edge
    NumberMap::NbNodesGmshEntity(66) = 11; NumberMap::NbVerticesGmshEntity(66) = 2;
    NumberMap::DimensionGmshEntity(66) = 1; NumberMap::OrderGmshEntity(66) = 10;    
    
    // 71 : 84-node sixth order tetrahedron
    NumberMap::NbNodesGmshEntity(71) = 84; NumberMap::NbVerticesGmshEntity(71) = 4;
    NumberMap::DimensionGmshEntity(71) = 3; NumberMap::OrderGmshEntity(71) = 6;
    // 72 : 120-node seventth order tetrahedron
    NumberMap::NbNodesGmshEntity(72) = 120; NumberMap::NbVerticesGmshEntity(72) = 4;
    NumberMap::DimensionGmshEntity(72) = 3; NumberMap::OrderGmshEntity(72) = 7;
    // 73 : 165-node eighth order tetrahedron
    NumberMap::NbNodesGmshEntity(73) = 165; NumberMap::NbVerticesGmshEntity(73) = 4;
    NumberMap::DimensionGmshEntity(73) = 3; NumberMap::OrderGmshEntity(73) = 8;
    // 74 : 220-node ninth order tetrahedron
    NumberMap::NbNodesGmshEntity(74) = 220; NumberMap::NbVerticesGmshEntity(74) = 4;
    NumberMap::DimensionGmshEntity(74) = 3; NumberMap::OrderGmshEntity(74) = 9;
    // 75 : 286-node tenth order tetrahedron
    NumberMap::NbNodesGmshEntity(75) = 286; NumberMap::NbVerticesGmshEntity(75) = 4;
    NumberMap::DimensionGmshEntity(75) = 3; NumberMap::OrderGmshEntity(75) = 10;
    
    // 76 : 74-node sixth order tetrahedron
    NumberMap::NbNodesGmshEntity(79) = 74; NumberMap::NbVerticesGmshEntity(79) = 4;
    NumberMap::DimensionGmshEntity(79) = 3; NumberMap::OrderGmshEntity(79) = 6;
    // 77 : 100-node seventh order tetrahedron
    NumberMap::NbNodesGmshEntity(80) = 100; NumberMap::NbVerticesGmshEntity(80) = 4;
    NumberMap::DimensionGmshEntity(80) = 3; NumberMap::OrderGmshEntity(80) = 7;
    // 78 : 135-node eighth order tetrahedron
    NumberMap::NbNodesGmshEntity(81) = 135; NumberMap::NbVerticesGmshEntity(81) = 4;
    NumberMap::DimensionGmshEntity(81) = 3; NumberMap::OrderGmshEntity(81) = 8;
    // 79 : 164-node ninth order tetrahedron
    NumberMap::NbNodesGmshEntity(82) = 164; NumberMap::NbVerticesGmshEntity(82) = 4;
    NumberMap::DimensionGmshEntity(82) = 3; NumberMap::OrderGmshEntity(82) = 9;
    // 80 : 202-node tenth order tetrahedron
    NumberMap::NbNodesGmshEntity(83) = 202; NumberMap::NbVerticesGmshEntity(83) = 4;
    NumberMap::DimensionGmshEntity(83) = 3; NumberMap::OrderGmshEntity(83) = 10;
    
    // 90 : third order prism
    NumberMap::NbNodesGmshEntity(90) = 40; NumberMap::NbVerticesGmshEntity(90) = 6;
    NumberMap::DimensionGmshEntity(90) = 3; NumberMap::OrderGmshEntity(90) = 3;
    // 91 : fourth order prism
    NumberMap::NbNodesGmshEntity(91) = 75; NumberMap::NbVerticesGmshEntity(91) = 6;
    NumberMap::DimensionGmshEntity(91) = 3; NumberMap::OrderGmshEntity(91) = 4;        
    // 92 : third order hexahedron
    NumberMap::NbNodesGmshEntity(92) = 64; NumberMap::NbVerticesGmshEntity(92) = 8;
    NumberMap::DimensionGmshEntity(92) = 3; NumberMap::OrderGmshEntity(92) = 3;
    // 93 : fourth order hexahedron
    NumberMap::NbNodesGmshEntity(93) = 125; NumberMap::NbVerticesGmshEntity(93) = 8;
    NumberMap::DimensionGmshEntity(93) = 3; NumberMap::OrderGmshEntity(93) = 4;
    // 94 : fifth order hexahedron
    NumberMap::NbNodesGmshEntity(94) = 216; NumberMap::NbVerticesGmshEntity(94) = 8;
    NumberMap::DimensionGmshEntity(94) = 3; NumberMap::OrderGmshEntity(94) = 5;
    // 95 : sixth order hexahedron
    NumberMap::NbNodesGmshEntity(95) = 343; NumberMap::NbVerticesGmshEntity(95) = 8;
    NumberMap::DimensionGmshEntity(95) = 3; NumberMap::OrderGmshEntity(95) = 6;
    // 96 : seventh order hexahedron
    NumberMap::NbNodesGmshEntity(96) = 512; NumberMap::NbVerticesGmshEntity(96) = 8;
    NumberMap::DimensionGmshEntity(96) = 3; NumberMap::OrderGmshEntity(96) = 7;
    // 97 : eighth order hexahedron
    NumberMap::NbNodesGmshEntity(97) = 729; NumberMap::NbVerticesGmshEntity(97) = 8;
    NumberMap::DimensionGmshEntity(97) = 3; NumberMap::OrderGmshEntity(97) = 8;
    // 98 : ninth order hexahedron
    NumberMap::NbNodesGmshEntity(98) = 1000; NumberMap::NbVerticesGmshEntity(98) = 8;
    NumberMap::DimensionGmshEntity(98) = 3; NumberMap::OrderGmshEntity(98) = 9;
    // 99 : incomplete third order hexahedron
    NumberMap::NbNodesGmshEntity(99) = 56; NumberMap::NbVerticesGmshEntity(99) = 8;
    NumberMap::DimensionGmshEntity(99) = 3; NumberMap::OrderGmshEntity(99) = 3;
    // 100 : incomplete fourth order hexahedron
    NumberMap::NbNodesGmshEntity(100) = 98; NumberMap::NbVerticesGmshEntity(100) = 8;
    NumberMap::DimensionGmshEntity(100) = 3; NumberMap::OrderGmshEntity(100) = 4;
    // 101 : incomplete fifth order hexahedron
    NumberMap::NbNodesGmshEntity(101) = 152; NumberMap::NbVerticesGmshEntity(101) = 8;
    NumberMap::DimensionGmshEntity(101) = 3; NumberMap::OrderGmshEntity(101) = 5;
    // 102 : incomplete sixth order hexahedron
    NumberMap::NbNodesGmshEntity(102) = 222; NumberMap::NbVerticesGmshEntity(102) = 8;
    NumberMap::DimensionGmshEntity(102) = 3; NumberMap::OrderGmshEntity(102) = 6;
    // 103 : incomplete seventh order hexahedron
    NumberMap::NbNodesGmshEntity(103) = 296; NumberMap::NbVerticesGmshEntity(103) = 8;
    NumberMap::DimensionGmshEntity(103) = 3; NumberMap::OrderGmshEntity(103) = 7;
    // 104 : incomplete eighth order hexahedron
    NumberMap::NbNodesGmshEntity(104) = 386; NumberMap::NbVerticesGmshEntity(104) = 8;
    NumberMap::DimensionGmshEntity(104) = 3; NumberMap::OrderGmshEntity(104) = 8;
    // 105 : incomplete ninth order hexahedron
    NumberMap::NbNodesGmshEntity(105) = 488; NumberMap::NbVerticesGmshEntity(105) = 8;
    NumberMap::DimensionGmshEntity(105) = 3; NumberMap::OrderGmshEntity(105) = 9;
    

    /**********************************
     *  Non standard entities in Gmsh *
     **********************************/
    
    // 109 : third order pyramid
    NumberMap::NbNodesGmshEntity(109) = 30; NumberMap::NbVerticesGmshEntity(109) = 5;
    NumberMap::DimensionGmshEntity(109) = 3; NumberMap::OrderGmshEntity(109) = 3;
    // 110 : fourth order pyramid
    NumberMap::NbNodesGmshEntity(110) = 55; NumberMap::NbVerticesGmshEntity(110) = 5;
    NumberMap::DimensionGmshEntity(110) = 3; NumberMap::OrderGmshEntity(110) = 4;
    // 111 : fifth order pyramid
    NumberMap::NbNodesGmshEntity(111) = 91; NumberMap::NbVerticesGmshEntity(111) = 5;
    NumberMap::DimensionGmshEntity(111) = 3; NumberMap::OrderGmshEntity(111) = 5;
    // 112 : fifth order prism
    NumberMap::NbNodesGmshEntity(112) = 126; NumberMap::NbVerticesGmshEntity(112) = 6;
    NumberMap::DimensionGmshEntity(112) = 3; NumberMap::OrderGmshEntity(112) = 5;

    // threshold for 2-D points
    R2::threshold = 1e-10;
  }

  
  //! default constructor
  template<class T>
  MeshNumbering<Dimension2, T>::MeshNumbering(Mesh<Dimension2, T>& mesh_)
    : MeshNumbering_Dim<Dimension2, T>(mesh_)
  {
  }
  

  //! returns maximal number of dofs on the edge i
  template<class T>
  int MeshNumbering<Dimension2, T>::GetNbDofElement(int i, const Edge<Dimension2>& elt) const
  {
    int r = this->order;
    if (this->type_variable_order != this->CONSTANT_ORDER)
      r = this->order_edge(i);
    
    int nb_dof =  this->number_map.GetNbDofVertex(r)*2;
    nb_dof += this->number_map.GetNbDofEdge(r);
    
    return nb_dof;
  }

  
  //! returns maximal number of dofs on the element i
  template<class T>
  int MeshNumbering<Dimension2, T>::GetNbDofElement(int i, const Face<Dimension2>& elt) const
  {
    int r = this->order;
    if (this->type_variable_order != this->CONSTANT_ORDER)
      r = this->order_element(i);
    
    int nb_dof =  this->number_map.GetNbDofVertex(r)*elt.GetNbVertices();
    nb_dof += this->number_map.GetNbDofEdge(r)*elt.GetNbEdges();
    
    if (elt.GetNbVertices() == 3)
      nb_dof += this->number_map.GetNbDofTriangle(r);
    else
      nb_dof += this->number_map.GetNbDofQuadrangle(r);
    
    return nb_dof;
  }


  //! returns the order of quadrature rules used in the mesh
  template<class T>
  void MeshNumbering<Dimension2, T>::GetOrderQuadrature(TinyVector<IVect, 4>& num) const
  {
    for (int t = 0; t < 4; t++)
      num(t).Clear();
    
    if (this->type_variable_order != this->CONSTANT_ORDER)
      {
	char rmax = 0;
	for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
	  rmax = max(rmax, this->order_boundary(i));
	
	Vector<bool> OrderUsed(rmax+1);
	OrderUsed.Fill(false);
	for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
	  OrderUsed(this->order_boundary(i)) = true;
        
	num(0).Clear();
	for (char r = 1; r <= rmax; r++)
          if (OrderUsed(r))
            num(0).PushBack(r);
      }
    else
      {
        num(0).Reallocate(1);
        num(0)(0) = this->order;
      }
  }
    
  
  //! has the face the same orientation in the two elements e1 and e2 ?
  template<class T>
  int MeshNumbering<Dimension2, T>
  ::GetBoundaryRotation(int f, int e1, int& e2, int pos1, int& pos2, int& rot) const
  {
    if (this->mesh.Boundary(f).GetNbElements() < 2)
      {
        e2 = -1;
        return -1;
      }

    e2 = this->mesh.Boundary(f).numElement(0);
    if (e1 == e2)
      e2 = this->mesh.Boundary(f).numElement(1);
    
    bool way_face1 = this->mesh.Element(e1).GetOrientationEdge(pos1);
    bool way_face2 = way_face1;
    pos2 = this->mesh.Element(e2).GetPositionBoundary(f); int f2 = -1;
    if (pos2 < 0)
      {
	f2 = this->PeriodicityBoundary(f);
        if (f2 >= 0)
          {
            pos2 = this->mesh.Element(e2).GetPositionBoundary(f2);
            way_face2 = this->mesh.Element(e2).GetOrientationEdge(pos2);
          }
      }
    else
      way_face2 = this->mesh.Element(e2).GetOrientationEdge(pos2);
    
    rot = 0;
    if (way_face1 != way_face2)
      rot = 1;
    
    return f2;
  }
  
  
  //! computes the orders associated for each edge, face, element
  template<class T>
  void MeshNumbering<Dimension2, T>::ComputeVariableOrder(int dg_form)
  {
    MeshNumbering_Dim<Dimension2, T>::ComputeVariableOrder(dg_form);
    
    if (dg_form == ElementReference_Base::CONTINUOUS
	&& (this->type_variable_order != this->CONSTANT_ORDER))
      {
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          {
            char r = this->order_inside(i);
	    for (int j = 0; j < this->mesh.Element(i).GetNbEdges(); j++)
	      r = max(r, this->order_edge(this->mesh.Element(i).numEdge(j)));
            
            this->order_element(i) = r;
          }
      }

    if (dg_form == ElementReference_Base::HDG
	&& (this->type_variable_order != this->CONSTANT_ORDER))
      {
	// order for each edge is equal to the order of quadrature
	for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
	  this->order_edge(i) = this->order_boundary(i);
      }
  }
  
  
  //! sets the order to use in the mesh (uniform, variable)
  template<class T>
  void MeshNumbering<Dimension2, T>::SetVariableOrder(int type)
  {
    this->type_variable_order = type;
    if (type == this->USER_ORDER)
      {
        this->order_element.Reallocate(this->mesh.GetNbElt());
        this->order_inside.Reallocate(this->mesh.GetNbElt());
	this->order_element.Fill(this->order);
        this->order_inside.Fill(this->order);
        
	this->order_edge.Reallocate(this->mesh.GetNbEdges());
        this->order_boundary.Reallocate(this->mesh.GetNbEdges());
        this->order_edge.Fill(this->order);
        this->order_boundary.Fill(this->order);
      }
    else if (type == this->CONSTANT_ORDER)
      {
        this->order_element.Clear(); this->order_inside.Clear();
        this->order_edge.Clear();
        this->order_boundary.Clear();
      }

  }
  
  
  //! numbering of faces
  template<class T>
  void MeshNumbering<Dimension2, T>::NumberFaces()
  {
    // no face in 2-D
  }


  //! numbering of edges
  template<class T>
  void MeshNumbering<Dimension2, T>::NumberSurfaceEdges()
  {
    // nothing to do
  }
  
  
  //! reconstructs arrays containing dof offsets for vertices, edges
  template<class T>
  void MeshNumbering<Dimension2, T>::ReconstructOffsetDofs()
  {
    if (this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    if (this->number_map.FormulationDG() == ElementReference_Base::HDG)
      {
	for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  {
	    int nb_edges_loc = this->mesh.Element(i).GetNbEdges();
	    int offset_loc = 0;
	    for (int j = 0; j < nb_edges_loc; j++)
	      {
		int ne = this->mesh.Element(i).numEdge(j);
		int re = this->GetOrderQuadrature(ne);
		int nb_dof_edge = this->number_map.GetNbDofEdge(re);
		this->OffsetDofFaceNumber(ne) = this->Element(i).GetNumberDof(offset_loc);
		for (int k = 1; k < nb_dof_edge; k++)
		  this->OffsetDofFaceNumber(ne)
		    = min(this->OffsetDofFaceNumber(ne), this->Element(i).GetNumberDof(offset_loc+k));

		offset_loc += nb_dof_edge;
	      }
	  }
	
	return;
      }
      
    int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
    int nb_max_vertex = 0, nb_max_edge = 0;
    
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int offset_loc = 0;
        if (nb_dof_vertex > 0)
          for (int j = 0; j < this->mesh.Element(i).GetNbVertices(); j++)
            {
              int nv = this->mesh.Element(i).numVertex(j);
              this->OffsetDofVertexNumber(nv) = this->Element(i).GetNumberDof(offset_loc);
              offset_loc += nb_dof_vertex;
              nb_max_vertex = max(nb_max_vertex,
				  this->OffsetDofVertexNumber(nv) + nb_dof_vertex);
            }
      }

    nb_max_edge = nb_max_vertex;
    this->OffsetDofEdgeNumber.Fill(nb_max_vertex);        
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int offset_loc = nb_dof_vertex*this->mesh.Element(i).GetNbVertices();
        int nb_edges_loc = this->mesh.Element(i).GetNbEdges();
        int r = this->GetOrderElement(i);
        int nb_dof_edge_nominal = this->number_map.GetNbDofEdge(r);
        for (int j = 0; j < nb_edges_loc; j++)
          {
            int ne = this->mesh.Element(i).numEdge(j);
            int re = this->GetOrderEdge(ne);
            int nb_dof_edge = this->number_map.GetNbDofEdge(re);

            if (nb_dof_edge > 0)
              {
                this->OffsetDofEdgeNumber(ne) = this->Element(i).GetNumberDof(offset_loc);
                for (int k = 1; k < nb_dof_edge; k++)
                  {
                    this->OffsetDofEdgeNumber(ne)
                      = min(this->OffsetDofEdgeNumber(ne),
                            this->Element(i).GetNumberDof(offset_loc+k));
                  }
                
                nb_max_edge = max(nb_max_edge, this->OffsetDofEdgeNumber(ne) + nb_dof_edge);
              }
            
            offset_loc += nb_dof_edge_nominal;
          }
      }
    
    this->OffsetDofVertexNumber(this->mesh.GetNbVertices()) = nb_max_vertex;
    this->OffsetDofEdgeNumber(this->mesh.GetNbEdges()) = nb_max_edge;
  }
  
  
  //! returns the number of dofs for an edge/face including vertices
  template<class T>
  int MeshNumbering<Dimension2, T>::GetNbDofBoundary(int i)
  {
    if (this->type_variable_order == this->CONSTANT_ORDER)
      return 2*this->number_map.GetNbDofVertex(this->order)
	+ this->number_map.GetNbDofEdge(this->order);
    else
      return 2*this->number_map.GetNbDofVertex(this->order)
	+ this->number_map.GetNbDofEdge(this->order_edge(i));
  }
  
  
  //! returns global dof numbers of boundary i
  template<class T>
  void MeshNumbering<Dimension2, T>::GetDofBoundary(int i, IVect& num)
  {   
    int nv = this->number_map.GetNbDofVertex(this->order);
    int ne = 0;
    if (this->type_variable_order == this->CONSTANT_ORDER)
      ne = this->number_map.GetNbDofEdge(this->order);
    else
      ne = this->number_map.GetNbDofEdge(this->order_edge(i));
    
    int nb_dof = 2*nv + ne;
    num.Reallocate(nb_dof);
    int n1 = this->mesh.BoundaryRef(i).numVertex(0);
    int n2 = this->mesh.BoundaryRef(i).numVertex(1);
    for (int k = 0; k < nv; k++)
      {
        num(k) = this->OffsetDofVertexNumber(n1) + k;
        num(k+nv) = this->OffsetDofVertexNumber(n2) + k;
      }
    
    for (int k = 0; k < ne; k++)
      num(k+2*nv) = this->OffsetDofEdgeNumber(i) + k;
    
  }
  
  
  //! constructing numbering used in Montjoie for nodes of a square
  template<class T>
  void MeshNumbering<Dimension2, T>::
  ConstructQuadrilateralNumbering(int order, Matrix<int>& NumNodes2D,
				  Matrix<int>& CoordinateNodes)
  {
    QuadrangleQuadrature::ConstructQuadrilateralNumbering(order, NumNodes2D, CoordinateNodes);
  }


  //! constructing numbering used in Montjoie for nodes of a triangle
  template<class T>
  void MeshNumbering<Dimension2, T>::
  ConstructTriangularNumbering(int order, Matrix<int>& NumNodes2D,
			       Matrix<int>& CoordinateNodes)
  {
    if (order < 0)
      {
	NumNodes2D.Clear();
	CoordinateNodes.Clear();
	return;
      }
      
    NumNodes2D.Reallocate(order+1, order+1);
    NumNodes2D.Fill(-1);
    // Num Nodes makes a bijection between a tensorial numbering (i,j)
    // three vertices of the triangle
    NumNodes2D(0, 0) = 0;
    NumNodes2D(order,0) = 1;
    NumNodes2D(0, order) = 2;
    
    // three edges
    int nb = 3;
    for (int i = 1; i < order; i++)
      NumNodes2D(i, 0) = nb++;

    for (int i = 1; i < order; i++)
      NumNodes2D(order-i, i) = nb++;

    for (int i = 1; i < order; i++)
      NumNodes2D(0, order-i) = nb++;
    
    // inside the triangle
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes2D(i, j) = nb++;
    
    CoordinateNodes.Reallocate((order+1)*(order+2)/2, 3);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	{
	  CoordinateNodes(NumNodes2D(i,j), 0) = order-i-j;
	  CoordinateNodes(NumNodes2D(i,j), 1) = i;
	  CoordinateNodes(NumNodes2D(i,j), 2) = j;
	}
    
  }


  //! treatment of periodic boundaries
  template<class T>
  void MeshNumbering<Dimension2, T>::TreatPeriodicCondition()
  {    
    int nb_boundaries_ref = this->mesh.GetNbBoundaryRef();
    
    // matching boundaries
    // PeriodicityBoundary(i) will be the global number
    // of the translated boundary to the global boundary i
    this->PeriodicityBoundary.Reallocate(nb_boundaries_ref);
    this->PeriodicityBoundary.Fill(-1);
    this->type_periodicity_boundary.Reallocate(nb_boundaries_ref);
    this->type_periodicity_boundary.Fill(0);
    this->TranslationPeriodicBoundary.Reallocate(nb_boundaries_ref);
    if (this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      this->type_formulation_periodic = this->WEAK_PERIODIC;

    if (this->number_map.FormulationDG() == ElementReference_Base::HDG)
      if (this->type_formulation_periodic == this->WEAK_PERIODIC)
	this->type_formulation_periodic = this->SAME_PERIODIC_DOFS;

    if (this->mesh.GetNbPeriodicReferences() <= 0)
      return;
    
    int nodl = this->nodl;
    int nb_vertices = this->mesh.GetNbVertices();
    
    IVect PeriodicityDof(nodl);
    IVect DofUsed(nodl); DofUsed.Fill(false);
    
    int nb_dof_periodic = 0;
    this->TranslationPeriodicDof.Reallocate(nodl);
    this->PeriodicDofs.Reallocate(nodl);
    this->PeriodicDofs_Original.Reallocate(nodl);
    this->type_periodicity_dof.Reallocate(nodl);
    
    // loop over all periodic boundaries specified by the user
    for (int num_per = 0; num_per < this->mesh.GetNbPeriodicReferences(); num_per++)
      {
	// references of the two corresponding boundaries
	int ref1 = this->mesh.GetPeriodicReference(num_per, 0);
	int ref2 = this->mesh.GetPeriodicReference(num_per, 1);
	
	// we get all the boundaries of reference 1, and of reference 2
	VectBool IsFace_Ref1(nb_boundaries_ref), IsFace_Ref2(nb_boundaries_ref);
	IsFace_Ref1.Fill(false);  IsFace_Ref2.Fill(false);
	int nb_faces_ref1 = 0, nb_faces_ref2 = 0;
	for (int i = 0; i < nb_boundaries_ref; i++)
	  {
	    if (this->mesh.BoundaryRef(i).GetReference() == ref1)
	      {
		nb_faces_ref1++;
		IsFace_Ref1(i) = true;
	      }
	    
	    if (this->mesh.BoundaryRef(i).GetReference() == ref2)
	      {
		nb_faces_ref2++;
		IsFace_Ref2(i) = true;
	      }
	  }

	//  computation of centers of faces
	VectR2 barycenter_faces1(nb_faces_ref1), barycenter_faces2(nb_faces_ref2);

	// index ref 1 :  face numbers related to first reference
	// index ref 2 :  face numbers related to second reference
	IVect Index_ref1(nb_faces_ref1), Index_ref2(nb_faces_ref2);
	R2 center_ref1, center_ref2;
	nb_faces_ref1 = 0; nb_faces_ref2 = 0;
	
	// we fill TypeVertex 
	// TypeVertex(i) = 1 if vertex i belongs to ref1
	// TypeVertex(i) = 2 if vertex i belongs to ref2
	// TypeVertex(i) = 3 if vertex i belongs to both boundaries
	// we get all the vertices of reference 1 and 2
	IVect TypeVertex(nb_vertices); TypeVertex.Fill(-1);
	for (int i = 0; i < nb_boundaries_ref; i++)
	  {
	    if (this->mesh.BoundaryRef(i).GetReference() == ref1)
	      {
		Index_ref1(nb_faces_ref1) = i;
	
		for (int j = 0; j < this->mesh.BoundaryRef(i).GetNbVertices(); j++)
		  {
		    barycenter_faces1(nb_faces_ref1) += this->mesh.Vertex(this->mesh.BoundaryRef(i).numVertex(j));
                    
		    // for circular periodicity, the center is both in ref1 and ref2
		    if (TypeVertex(this->mesh.BoundaryRef(i).numVertex(j)) == 2)
		      TypeVertex(this->mesh.BoundaryRef(i).numVertex(j)) = 3;
		    else
                      TypeVertex(this->mesh.BoundaryRef(i).numVertex(j)) = 1;
		  }
		
		Mlt(T(1.0/this->mesh.BoundaryRef(i).GetNbVertices()),
                    barycenter_faces1(nb_faces_ref1));
		center_ref1 += barycenter_faces1(nb_faces_ref1);
		nb_faces_ref1++;
	      }
	    
	    if (this->mesh.BoundaryRef(i).GetReference() == ref2)
	      {
		Index_ref2(nb_faces_ref2) = i;
	
		for (int j = 0; j < this->mesh.BoundaryRef(i).GetNbVertices(); j++)
		  {
		    barycenter_faces2(nb_faces_ref2) += this->mesh.Vertex(this->mesh.BoundaryRef(i).numVertex(j));
                        
		    // for circular periodicity, the center is both in ref1 and ref2
                    // -> TypeVertex = 3
		    if (TypeVertex(this->mesh.BoundaryRef(i).numVertex(j)) == 1)
		      TypeVertex(this->mesh.BoundaryRef(i).numVertex(j)) = 3;
		    else
		      TypeVertex(this->mesh.BoundaryRef(i).numVertex(j)) = 2;
		  }
		
		Mlt(T(1)/this->mesh.BoundaryRef(i).GetNbVertices(),
                    barycenter_faces2(nb_faces_ref2));
		
		center_ref2 += barycenter_faces2(nb_faces_ref2);
		nb_faces_ref2++;
	      }
	  }
                
	// center_ref1 and center_ref2 are the centers of the two boundaries ref1 and ref2
	Mlt(T(1.0/nb_faces_ref1), center_ref1);   Mlt(T(1.0/nb_faces_ref2), center_ref2);

	// number of vertices on each boundary
	int nb_vertices_ref1(0), nb_vertices_ref2(0);
	for (int i = 0; i < nb_vertices; i++)
	  {
	    if (TypeVertex(i) == 1)
	      nb_vertices_ref1++;
	    else if (TypeVertex(i) == 2)
	      nb_vertices_ref2++;
	    else if (TypeVertex(i) == 3)
	      {
		nb_vertices_ref1++;
		nb_vertices_ref2++;
	      }
	  }
	
	// vertices of reference 1 and 2, construction of indexing arrays
	VectR2 Vertices_Ref1(nb_vertices_ref1);
	VectR2 Vertices_Ref2(nb_vertices_ref2);
	
	// Index_Vert1(i) : vertex number of the vertex i belonging to ref 1
        // Index_Vert2(i) : vertex number of the vertex i belonging to ref 2
	// Vertices_Ref1(i)  : coordinates 
	// Vertices_Ref2(i) : coordinates
	IVect Index_Vert1(nb_vertices_ref1), Index_Vert2(nb_vertices_ref2);
	
	nb_vertices_ref1 = 0; nb_vertices_ref2 = 0;
	for (int i = 0; i < nb_vertices; i++)
	  {
	    // vertices belonging to boundary ref1 are retrieved
	    // vertex numbers are in Index_Vert1, coordinates in Vertices_Ref1
	    if ((TypeVertex(i) == 1)||(TypeVertex(i) == 3))
	      {
		Index_Vert1(nb_vertices_ref1) = i;
		Vertices_Ref1(nb_vertices_ref1) = this->mesh.Vertex(i);
	       
		nb_vertices_ref1++;
	      }
	    
	    // vertices belonging to boundary ref2 are retrieved
	    // vertex numbers are in Index_Vert2, coordinates in Vertices_Ref2
	    if ((TypeVertex(i) == 2)||(TypeVertex(i) == 3))
	      {
		Index_Vert2(nb_vertices_ref2) = i;
		Vertices_Ref2(nb_vertices_ref2) = this->mesh.Vertex(i);
		nb_vertices_ref2++;
	      }
	  }
	
	// translation vector between the two boundaries ?
	R2 vh_translat;
        int type_per = this->mesh.GetPeriodicityTypeReference(num_per);
        if (type_per == BoundaryConditionEnum::PERIODIC_THETA) 
	  {
            Real_wp theta = this->mesh.GetPeriodicAlpha();
            
	    // premiere composante du vecteur translation est vh
	    vh_translat(0) = theta;
	    
	    // centers of faces ref2 are rotated in order to compare them
	    // with centers of faces ref1
	    for (int i = 0; i < nb_faces_ref2; i++)
	      {
		// rotation de - theta
		R2 temp = barycenter_faces2(i) ; 
		barycenter_faces2(i)(0)
                  = cos( -vh_translat(0) ) * temp(0) - sin(- vh_translat(0)) *  temp(1);
		barycenter_faces2(i)(1)
                  = sin(-vh_translat(0) ) * temp(0) + cos(-vh_translat(0) ) *  temp(1); 
	      }                       
	  }
	else 
	  {
	    // cartesian coordinates
	    vh_translat = center_ref2 - center_ref1;
	    
	    // all centers of faces of ref2 are translated
	    for (int i = 0; i < nb_faces_ref2; i++)
	      barycenter_faces2(i) -= vh_translat;
          }
        
        if ((type_per != BoundaryConditionEnum::PERIODIC_CTE)
            && !this->plane_wave_quasi_periodic)
          {
            // storing 0 or 1 for translations
            if (abs(vh_translat(0)) < R2::threshold)
              vh_translat(0) = 0.0;
            else if (vh_translat(0) < 0)
              vh_translat(0) = -1.0;
            else if (vh_translat(0) > 0)
              vh_translat(0) = 1.0;
            
            if (abs(vh_translat(1)) < R2::threshold)
              vh_translat(1) = 0.0;
            else if (vh_translat(1) < 0)
              vh_translat(1) = -1.0;
            else if (vh_translat(1) > 0)
              vh_translat(1) = 1.0;
          }
        
	// all centers of faces are sorted, so that we can compare them
	Sort(nb_faces_ref1, barycenter_faces1, Index_ref1);
	Sort(nb_faces_ref2, barycenter_faces2, Index_ref2);
	if (type_per == BoundaryConditionEnum::PERIODIC_THETA)
	  {
	    // vertices of ref2 are rotated in order to compare them
	    // with vertices of ref1
	    for (int i = 0; i < nb_vertices_ref2; i++)
	      {
		R2 temp = Vertices_Ref2(i) ; 
		// rotation de - theta
		Vertices_Ref2(i)(0)
                  = cos( -vh_translat(0) ) * temp(0) - sin(- vh_translat(0)) *  temp(1);
		Vertices_Ref2(i)(1)
                  = sin(-vh_translat(0) ) * temp(0) + cos(-vh_translat(0) ) *  temp(1);
	      }
	  }
	else
	  {
	    // all vertices of ref2 are translated
	    for (int i = 0; i < nb_vertices_ref2; i++)
	      Vertices_Ref2(i) -= vh_translat;
	  }
	
	// vertices are sorted to have direct matching
	Sort(nb_vertices_ref1, Vertices_Ref1, Index_Vert1);
	Sort(nb_vertices_ref2, Vertices_Ref2, Index_Vert2);
	
	IVect InverseIndex_Vert(nb_vertices);
	// InverseIndex_Vert(i) -> numero local (parmi tous les sommets de ref1 )
        // du sommet global i
	for (int i = 0; i < nb_vertices_ref1; i++)
	  InverseIndex_Vert(Index_Vert1(i)) = i;
	
	// InverseIndex_Vert(i) -> numero local (parmi tous les sommets de ref2 )
        // du sommet global i
	for (int i = 0; i < nb_vertices_ref2; i++)
	  InverseIndex_Vert(Index_Vert2(i)) = i;
	
	// on construit la relation face <-> face_ref
	PeriodicityDof.Fill(-1);
	
	// degrees of freedom associated with vertices are treated
	int nb_dof_vertex = this->number_map.GetNbDofVertex(this->GetOrder());
	if (nb_dof_vertex > 0)
	  {
	    if (nb_vertices_ref1 != nb_vertices_ref2)
	      {
		cout << "the this->mesh should be the same between periodic boundaries" << endl;
		abort();
	      }
	    
	    if (this->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
	      {
		for (int i = 0; i < nb_vertices_ref1; i++)
		  {
		    int n1 = Index_Vert1(i);
		    int n2 = Index_Vert2(i);
		    for (int k = 0; k < nb_dof_vertex; k++)
		      PeriodicityDof(n2*nb_dof_vertex + k) = this->OffsetDofVertexNumber(n1) + k;
		  }
	      }
	  }
	
	// periodicite pour les aretes
	// on a la correspondance immediatement
	for (int i = 0; i < nb_faces_ref1; i++)
	  {
	    int num_face1_ref = Index_ref1(i), num_face2_ref = Index_ref2(i);
	    int num_face1 = Index_ref1(i);
	    int num_face2 = Index_ref2(i);
	    int num_elem1 = this->mesh.Boundary(num_face1).numElement(0);
	    int num_elem2 = this->mesh.Boundary(num_face2).numElement(0);
	    
	    if (this->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
	      {
		int r = this->GetOrderEdge(num_face1);
		int nb_dof_edge = this->number_map.GetNbDofEdge(r);
		
		int offset1 = this->OffsetDofEdgeNumber(num_face1);
		int offset2 = this->OffsetDofEdgeNumber(num_face2);
		
		for (int k = 0; k < nb_dof_edge; k++)
		  PeriodicityDof(offset2 + k) = offset1 + k;	    
	      }
	    
            this->mesh.Boundary(num_face1).AddElement(num_elem2);
            this->mesh.Boundary(num_face2).AddElement(num_elem1);
            
            this->PeriodicityBoundary(num_face1_ref) = num_face2_ref;
            this->PeriodicityBoundary(num_face2_ref) = num_face1_ref;
            this->type_periodicity_boundary(num_face1_ref) = type_per;
            this->type_periodicity_boundary(num_face2_ref) = type_per;
            this->TranslationPeriodicBoundary(num_face2_ref) = vh_translat;
            this->TranslationPeriodicBoundary(num_face1_ref) = -vh_translat;
          }
	
	// pour tous les ddl periodiques
	for (int i = 0; i < nodl; i++)
	  if (i != PeriodicityDof(i))
	    if ((PeriodicityDof(i) >= 0)&&(!DofUsed(i)))
	      {
                this->TranslationPeriodicDof(nb_dof_periodic) = vh_translat;
                this->PeriodicDofs(nb_dof_periodic) = i;
		this->PeriodicDofs_Original(nb_dof_periodic) = PeriodicityDof(i);
                this->type_periodicity_dof(nb_dof_periodic) = type_per;
		DofUsed(i) = true;
		nb_dof_periodic++;
	      }
      }
    
    this->PeriodicDofs.Resize(nb_dof_periodic);
    this->PeriodicDofs_Original.Resize(nb_dof_periodic);
    this->type_periodicity_dof.Resize(nb_dof_periodic);
    this->TranslationPeriodicDof.Resize(nb_dof_periodic);

    IVect IndexPeriodicDof(nodl);
    IndexPeriodicDof.Fill(-1);
    
    // searching original dofs
    for (int i = 0; i < nb_dof_periodic; i++)
      IndexPeriodicDof(this->PeriodicDofs(i)) = i;
    
    for (int i = 0; i < nb_dof_periodic; i++)
      {
        R2 v1 = this->TranslationPeriodicDof(i), v2, v3;
        int t1 = this->type_periodicity_dof(i), t2 = -1, t3 = -1;
        int n1 = this->PeriodicDofs_Original(i);
        int n2 = IndexPeriodicDof(n1);
        if (n2 >= 0)
          {
            t2 = this->type_periodicity_dof(n2);
            v2 = this->TranslationPeriodicDof(n2);
            n1 = this->PeriodicDofs_Original(n2);
            n2 = IndexPeriodicDof(n1);
            if (n2 >= 0)
              {
                abort();
              }
          }
        
        this->PeriodicDofs_Original(i) = n1;
        this->type_periodicity_dof(i) = BoundaryConditionEnum::GetCompositionPeriodicity(t1, t2, t3);
        this->TranslationPeriodicDof(i) = v1 + v2;
      }
    
    if ((this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS) 
        || this->type_formulation_periodic == this->WEAK_PERIODIC)
      {
        this->TranslationPeriodicDof.Clear();
        this->PeriodicDofs.Clear();
        this->PeriodicDofs_Original.Clear();
        this->type_periodicity_dof.Clear();        
      }
    else if (this->UseSameNumberForPeriodicDofs())
      {
	// we associate same dof numbers for periodic dofs
	// total number of dofs is thus decreased
	// this option can be used for a periodic condition (and not for a quasi-periodic)
	IVect NewNumber(nodl); NewNumber.Fill(-1); int nb = 0;	
	for (int i = 0; i < nodl; i++)
	  if (!DofUsed(i))
	    NewNumber(i) = nb++;
		
	for (int i = 0; i < nb_dof_periodic; i++)
	  {
            int n1 = this->PeriodicDofs_Original(i);
	    NewNumber(this->PeriodicDofs(i)) = NewNumber(n1);
	  }
        
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  {
            int nb_dof = this->Element(i).GetNbDof();
                            
            for (int j = 0; j < nb_dof; j++)
              {
                int num_dof = this->Element(i).GetNumberDof(j);                
                this->Element(i).SetNumberDof(j, NewNumber(num_dof));
              }
          }
        
	this->nodl = nb;
        this->ReconstructOffsetDofs();
	nb_dof_periodic = 0;
        this->PeriodicDofs.Clear();
        this->PeriodicDofs_Original.Clear();
        this->type_periodicity_dof.Clear();
	this->TranslationPeriodicDof.Clear();
      }
    else
      {
        //this->TranslationPeriodicBoundary.Clear();
        //this->PeriodicityBoundary.Fill(-1);
      }
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  template<class T>
  TinyMatrix<int, General, 2, 4> MeshNumbering<Dimension2, T>::FirstExtremityEdge;
  
  template<class T>
  TinyMatrix<int, General, 2, 4> MeshNumbering<Dimension2, T>::SecondExtremityEdge;
#endif
  
  
  /****************************
   * MeshNumbering<Dimension3> *
   ****************************/
  
  
#ifdef MONTJOIE_WITH_THREE_DIM


#ifndef SELDON_WITH_COMPILED_LIBRARY
  template<class T>
  TinyMatrix<int, General, 4, 12> MeshNumbering<Dimension3, T>::FirstExtremityEdge;
  
  template<class T>
  TinyMatrix<int, General, 4, 12> MeshNumbering<Dimension3, T>::SecondExtremityEdge;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::FirstExtremityFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::SecondExtremityFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::ThirdExtremityFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::FourthExtremityFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::FirstEdgeFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::SecondEdgeFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::ThirdEdgeFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 6> MeshNumbering<Dimension3, T>::FourthEdgeFace;
  
  template<class T>
  TinyMatrix<int, General, 4, 8> MeshNumbering<Dimension3, T>::FirstVertexTriedre;
  
  template<class T>
  TinyMatrix<int, General, 4, 8> MeshNumbering<Dimension3, T>::SecondVertexTriedre;
  
  template<class T>
  TinyMatrix<int, General, 4, 8> MeshNumbering<Dimension3, T>::ThirdVertexTriedre;

  template<class T>
  TinyMatrix<int, General, 4, 8> MeshNumbering<Dimension3, T>::FirstEdgeTriedre;
  
  template<class T>
  TinyMatrix<int, General, 4, 8> MeshNumbering<Dimension3, T>::SecondEdgeTriedre;
  
  template<class T>
  TinyMatrix<int, General, 4, 8> MeshNumbering<Dimension3, T>::ThirdEdgeTriedre;

  template<class T>
  TinyMatrix<bool, General, 4, 8> MeshNumbering<Dimension3, T>::BeginFirstEdgeTriedre;
  
  template<class T>
  TinyMatrix<bool, General, 4, 8> MeshNumbering<Dimension3, T>::BeginSecondEdgeTriedre;
  
  template<class T>
  TinyMatrix<bool, General, 4, 8> MeshNumbering<Dimension3, T>::BeginThirdEdgeTriedre;

  template<class T>
  TinyMatrix<int, General, 6, 6> MeshNumbering<Dimension3, T>::RotationFaceTri;
  
  template<class T>
  TinyMatrix<int, General, 8, 8> MeshNumbering<Dimension3, T>::RotationFaceQuad;
  
  template<class T>
  TinyMatrix<bool, General, 4, 6> MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral;
#endif
  
  
  //! initialization of static parameters
  template<class T>
  void InitStaticMeshData(const Dimension3& dim, const T& x)
  {
    // we use here hybrid type, ie 0 -> tetra, 1 -> pyramid, 2 -> wedge, 3 -> hexa
    MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral.Fill(false);
    // all faces of hexahedron are quadrilateral
    for (int j = 0; j < 6; j++)
      MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(3, j) = true;
    
    // first face of a pyramid is quadrilateral
    MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(1, 0) = true;
    // three faces of a wedge are quadrilateral
    MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(2, 1) = true;
    MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(2, 2) = true;
    MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(2, 3) = true;
    
    // filling arrays containing extremities of edges of reference elements
    MeshNumbering<Dimension3, T>::FirstExtremityEdge.Fill(-1);
    MeshNumbering<Dimension3, T>::SecondExtremityEdge.Fill(-1);
    
    // 6 edges for tetrahedron
    string data("0 1 \n 0 2 \n 0 3 \n 1 2 \n 1 3 \n 2 3 ");
    istringstream stream_data(data);
    
    Matrix<int> tet_edge; tet_edge.ReadText(stream_data);
    for (int i = 0; i < 6; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityEdge(0, i) = tet_edge(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityEdge(0, i) = tet_edge(i, 1);
      }
    
    // 8 edges for pyramid
    data = string("0 1 \n 1 2 \n 3 2 \n 0 3 \n 0 4 \n 1 4 \n 2 4 \n 3 4 ");
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> pyr_edge; pyr_edge.ReadText(stream_data);
    for (int i = 0; i < 8; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityEdge(1, i) = pyr_edge(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityEdge(1, i) = pyr_edge(i, 1);
      }
    
    // 9 edges for wedge
    data = string("0 1 \n 1 2 \n 0 2 \n 0 3 \n 1 4 \n 2 5 \n 3 4 \n 4 5 \n 3 5 ");
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> wed_edge; wed_edge.ReadText(stream_data);
    for (int i = 0; i < 9; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityEdge(2, i) = wed_edge(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityEdge(2, i) = wed_edge(i, 1);
      }

    // 12 edges for cube
    Matrix<int> hex_edge;
    HexahedronQuadrature::GetEdgeCube(hex_edge);
    
    for (int i = 0; i < 12; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityEdge(3, i) = hex_edge(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityEdge(3, i) = hex_edge(i, 1);
      }
    
    // filling arrays containing extremities of faces of reference elements
    MeshNumbering<Dimension3, T>::FirstExtremityFace.Fill(-1);
    MeshNumbering<Dimension3, T>::SecondExtremityFace.Fill(-1);
    MeshNumbering<Dimension3, T>::ThirdExtremityFace.Fill(-1);
    MeshNumbering<Dimension3, T>::FourthExtremityFace.Fill(-1);
    
    // four faces for tetrahedron
    data = string("0 1 2 \n 0 1 3 \n 0 2 3 \n 1 2 3");
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> tet_face; tet_face.ReadText(stream_data);
    for (int i = 0; i < 4; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityFace(0, i) = tet_face(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityFace(0, i) = tet_face(i, 1);
	MeshNumbering<Dimension3, T>::ThirdExtremityFace(0, i) = tet_face(i, 2);
      }

    // five faces for pyramid
    data = string("0 1 2 3\n 0 1 4 -1\n 1 2 4 -1\n 3 2 4 -1\n 0 3 4 -1");
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> pyr_face; pyr_face.ReadText(stream_data);
    for (int i = 0; i < 5; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityFace(1, i) = pyr_face(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityFace(1, i) = pyr_face(i, 1);
	MeshNumbering<Dimension3, T>::ThirdExtremityFace(1, i) = pyr_face(i, 2);
	MeshNumbering<Dimension3, T>::FourthExtremityFace(1, i) = pyr_face(i, 3);
      }
    
    // five faces for wedge
    data = string("0 1 2 -1\n 0 1 4 3\n 1 2 5 4\n 0 2 5 3\n 3 4 5 -1");
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> wed_face; wed_face.ReadText(stream_data);
    for (int i = 0; i < 5; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityFace(2, i) = wed_face(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityFace(2, i) = wed_face(i, 1);
	MeshNumbering<Dimension3, T>::ThirdExtremityFace(2, i) = wed_face(i, 2);
	MeshNumbering<Dimension3, T>::FourthExtremityFace(2, i) = wed_face(i, 3);
      }
    
    // six faces for cube
    data = string("0 3 7 4 \n 0 1 5 4 \n 0 1 2 3 \n 4 5 6 7\n 3 2 6 7\n 1 2 6 5");
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> hex_face; hex_face.ReadText(stream_data);
    for (int i = 0; i < 6; i++)
      {
	MeshNumbering<Dimension3, T>::FirstExtremityFace(3, i) = hex_face(i, 0);
	MeshNumbering<Dimension3, T>::SecondExtremityFace(3, i) = hex_face(i, 1);
	MeshNumbering<Dimension3, T>::ThirdExtremityFace(3, i) = hex_face(i, 2);
	MeshNumbering<Dimension3, T>::FourthExtremityFace(3, i) = hex_face(i, 3);
      }
    
    // now we find edges for each face
    MeshNumbering<Dimension3, T>::FirstEdgeFace.Fill(-1);
    MeshNumbering<Dimension3, T>::SecondEdgeFace.Fill(-1);
    MeshNumbering<Dimension3, T>::ThirdEdgeFace.Fill(-1);
    MeshNumbering<Dimension3, T>::FourthEdgeFace.Fill(-1);
    
    IVect n(4), ne(4);
    for (int num = 0; num < 4; num++)
      {
        for (int i = 0; i < 6; i++)
          {
            n(0) = MeshNumbering<Dimension3, T>::FirstExtremityFace(num, i);
            n(1) = MeshNumbering<Dimension3, T>::SecondExtremityFace(num, i);
            n(2) = MeshNumbering<Dimension3, T>::ThirdExtremityFace(num, i);
            n(3) = MeshNumbering<Dimension3, T>::FourthExtremityFace(num, i);
            
            int nv = 3;
            if (MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(num, i))
              nv = 4;
            
            ne.Fill(-1);
            // looking for each edge of the face
            for (int e = 0; e < nv; e++)
              {
                // two extremities of the edge
                int na = n(e), nb = n((e+1)%nv);
                
                // what is the edge number ?
                if ((na != -1)&&(nb != -1))
                  for (int k = 0; k < 12; k++)
                    {
                      int nc = MeshNumbering<Dimension3, T>::FirstExtremityEdge(num, k);
                      int nd = MeshNumbering<Dimension3, T>::SecondExtremityEdge(num, k);
                      if (((na == nc)&&(nb == nd))||((na == nd)&&(nb == nc)))
                        ne(e) = k;
                    }
              }
            
            MeshNumbering<Dimension3, T>::FirstEdgeFace(num, i) = ne(0);
            MeshNumbering<Dimension3, T>::SecondEdgeFace(num, i) = ne(1);
            MeshNumbering<Dimension3, T>::ThirdEdgeFace(num, i) = ne(2);
            MeshNumbering<Dimension3, T>::FourthEdgeFace(num, i) = ne(3);
          }
        
        // loop over points
        TinyVector<bool, 4> beg;
        for (int i = 0; i < 8; i++)
          {
            // loop over edges to find the extremities of the "triedre"
            // a triedre consists of four points, the center is the current vertex,
            // and the other points define the three directions
            int nb_ext = 0;
            for (int j = 0; j < 12; j++)
              {
                if (MeshNumbering<Dimension3, T>::FirstExtremityEdge(num, j) == i)
                  {
                    n(nb_ext) = MeshNumbering<Dimension3, T>::SecondExtremityEdge(num, j);
                    ne(nb_ext) = j; beg(nb_ext) = true;
                    nb_ext++;
                  }
                
                if (MeshNumbering<Dimension3, T>::SecondExtremityEdge(num, j) == i)
                  {
                    n(nb_ext) = MeshNumbering<Dimension3, T>::FirstExtremityEdge(num, j);
                    ne(nb_ext) = j; beg(nb_ext) = false;
                    nb_ext++;
                  }
              }
            
            // we store edges and vertices
            MeshNumbering<Dimension3, T>::FirstVertexTriedre(num, i) = n(0);
            MeshNumbering<Dimension3, T>::SecondVertexTriedre(num, i) = n(1);
            MeshNumbering<Dimension3, T>::ThirdVertexTriedre(num, i) = n(2);
            
            MeshNumbering<Dimension3, T>::FirstEdgeTriedre(num, i) = ne(0);
            MeshNumbering<Dimension3, T>::SecondEdgeTriedre(num, i) = ne(1);
            MeshNumbering<Dimension3, T>::ThirdEdgeTriedre(num, i) = ne(2);

            MeshNumbering<Dimension3, T>::BeginFirstEdgeTriedre(num, i) = beg(0);
            MeshNumbering<Dimension3, T>::BeginSecondEdgeTriedre(num, i) = beg(1);
            MeshNumbering<Dimension3, T>::BeginThirdEdgeTriedre(num, i) = beg(2);
          }
        
      }
    
    data = string("0 1 2 3 4 5 6 7 \n 3 0 1 2 7 4 5 6 \n 2 3 0 1 6 7 4 5 \n 1 2 3 0 5 6 7 4 \n")+
      string("4 7 6 5 0 3 2 1 \n 5 4 7 6 1 0 3 2 \n 6 5 4 7 2 1 0 3 \n 7 6 5 4 3 2 1 0 ");
    
    stream_data.clear(); stream_data.str(data);
    
    Matrix<int> rot_face; rot_face.ReadText(stream_data);
    for (int i = 0; i < rot_face.GetM(); i++)
      for (int j = 0; j < rot_face.GetN(); j++)
        MeshNumbering<Dimension3, T>::RotationFaceQuad(i, j) = rot_face(i, j);
    
    data = string("0 1 2 3 4 5 \n 2 0 1 5 3 4 \n 1 2 0 4 5 3 \n 3 5 4 0 2 1 \n")+
      string("4 3 5 1 0 2 \n 5 4 3 2 1 0");
    
    stream_data.clear(); stream_data.str(data);
    
    rot_face.ReadText(stream_data);
    for (int i = 0; i < rot_face.GetM(); i++)
      for (int j = 0; j < rot_face.GetN(); j++)
        MeshNumbering<Dimension3, T>::RotationFaceTri(i, j) = rot_face(i, j);
    
    // threshold for 3-D points
    R3::threshold = 1e-10;
  }
  
  
  //! default constructor
  template<class T>
  MeshNumbering<Dimension3, T>::MeshNumbering(Mesh<Dimension3, T>& mesh_)
    : MeshNumbering_Dim<Dimension3, T>(mesh_)
  {
  }
  

    //! returns the order of quadrature rules used in the mesh
  template<class T>
  void MeshNumbering<Dimension3, T>::GetOrderQuadrature(TinyVector<IVect, 4>& num) const
  {
    for (int t = 0; t < 4; t++)
      num(t).Clear();
    
    if (this->type_variable_order != this->CONSTANT_ORDER)
      {
	char rmax = 0;
	for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
	  rmax = max(rmax, this->order_boundary(i));

        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          rmax = max(rmax, this->order_element(i));
	
	Vector<bool> OrderUsed_Tri(rmax+1), OrderUsed_Quad(rmax+1);
	OrderUsed_Tri.Fill(false);
        OrderUsed_Quad.Fill(false);
	for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
	  {
            if (this->mesh.Boundary(i).GetHybridType() == 0)
              OrderUsed_Tri(this->order_boundary(i)) = true;
            else
              OrderUsed_Quad(this->order_boundary(i)) = true;
          }
        
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          for (int j = 0; j < this->mesh.Element(i).GetNbBoundary(); j++)
            {
              char r = this->order_element(i);
              int nf = this->mesh.Element(i).numBoundary(j);
              if (this->mesh.Boundary(nf).GetHybridType() == 0)
                OrderUsed_Tri(r) = true;
              else
                OrderUsed_Quad(r) = true;
            }
        
	num(0).Clear(); num(1).Clear();
	for (char r = 1; r <= rmax; r++)
          {
            if (OrderUsed_Tri(r))
              num(0).PushBack(r);

            if (OrderUsed_Quad(r))
              num(1).PushBack(r);
          }
      }
    else
      {
        bool presence_tri = false, presence_quad = false;
        for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
          {
            if (this->mesh.Boundary(i).GetHybridType() == 0)
              presence_tri = true;
            else
              presence_quad = true;
          }
        
        if (presence_tri)
          {
            num(0).Reallocate(1);
            num(0)(0) = this->order;
          }
        else
          num(0).Clear();
        
        if (presence_quad)
          {
            num(1).Reallocate(1);
            num(1)(0) = this->order;
          }
        else
          num(1).Clear();
      }
  }

  
  //! sets the order to use in the mesh (uniform, variable)
  template<class T>
  void MeshNumbering<Dimension3, T>::SetVariableOrder(int type)
  {
    this->type_variable_order = type;
    
    if (type == this->USER_ORDER)
      {
        this->order_element.Reallocate(this->mesh.GetNbElt());
        this->order_inside.Reallocate(this->mesh.GetNbElt());
	this->order_element.Fill(this->order);
        this->order_inside.Fill(this->order);
        
	this->order_edge.Reallocate(this->mesh.GetNbEdges());
        this->order_face.Reallocate(this->mesh.GetNbFaces());
        this->order_boundary.Reallocate(this->mesh.GetNbFaces());
        this->order_edge.Fill(this->order);
        this->order_face.Fill(this->order);
        this->order_boundary.Fill(this->order);
      }
    else if (type == this->CONSTANT_ORDER)
      {
        this->order_element.Clear(); this->order_inside.Clear();
        this->order_edge.Clear(); this->order_face.Clear();
        this->order_boundary.Clear();
      }
  }


  //! returns maximal number of dofs on the face i
  template<class T>
  int MeshNumbering<Dimension3, T>::GetNbDofElement(int i, const Face<Dimension3>& elt) const
  {
    int r = this->order;
    if (this->type_variable_order != this->CONSTANT_ORDER)
      r = this->order_face(i);
    
    int nb_dof =  this->number_map.GetNbDofVertex(r)*elt.GetNbVertices();
    nb_dof += this->number_map.GetNbDofEdge(r)*elt.GetNbEdges();
    
    if (elt.GetNbVertices() == 3)
      nb_dof += this->number_map.GetNbDofTriangle(r);
    else
      nb_dof += this->number_map.GetNbDofQuadrangle(r);
    
    return nb_dof;
  }

  
  //! returns maximal number of dofs on the element i
  template<class T>
  int MeshNumbering<Dimension3, T>::GetNbDofElement(int i, const Volume& elt) const
  {
    int nb_dof = 0;
    int r = this->order;
    if (this->type_variable_order != this->CONSTANT_ORDER)
      r = this->order_element(i);
    
    nb_dof += this->number_map.GetNbDofVertex(r)*elt.GetNbVertices();
    nb_dof += this->number_map.GetNbDofEdge(r)*elt.GetNbEdges();
        
    for (int j = 0; j < elt.GetNbFaces(); j++)
      {
	int nf = elt.numFace(j);
	if (this->mesh.Boundary(nf).GetNbVertices() == 3)
	  nb_dof += this->number_map.GetNbDofTriangle(r);
	else
	  nb_dof += this->number_map.GetNbDofQuadrangle(r);
      }
    
    switch (elt.GetNbVertices())
      {
      case 4:
	nb_dof += this->number_map.GetNbDofTetrahedron(r);
	break;
      case 5:
	nb_dof += this->number_map.GetNbDofPyramid(r);
	break;
      case 6:
	nb_dof += this->number_map.GetNbDofWedge(r);
	break;
      case 8:
	nb_dof += this->number_map.GetNbDofHexahedron(r);
	break;
      }
    
    return nb_dof;
  }
  

  //! has the face the same orientation in the two elements e1 and e2 ?
  template<class T>
  int MeshNumbering<Dimension3, T>
  ::GetBoundaryRotation(int f, int e1, int& e2, int pos1, int& pos2, int& rot) const
  {
    if (this->mesh.Boundary(f).GetNbElements() < 2)
      {
        e2 = -1;
        return -1;
      }
    
    e2 = this->mesh.Boundary(f).numElement(0);
    if (e1 == e2)
      e2 = this->mesh.Boundary(f).numElement(1);
    
    int way_face1 = this->mesh.Element(e1).GetOrientationFace(pos1);
    int way_face2 = way_face1;
    pos2 = this->mesh.Element(e2).GetPositionBoundary(f);
    int f2 = -1;
    if (pos2 < 0)
      {
	int f2 = this->PeriodicityBoundary(f);
        if (f2 >= 0)
          {            
            pos2 = this->mesh.Element(e2).GetPositionBoundary(f2);
            way_face2 = this->mesh.Element(e2).GetOrientationFace(pos2);
          }
      }
    else
      way_face2 = this->mesh.Element(e2).GetOrientationFace(pos2);
    
    int nv = this->mesh.Boundary(f).GetNbVertices();
    rot = this->GetRotationFace(way_face1, way_face2, nv);
    
    return f2;
  }
  
    
  //! computes the order associated for each edge, face, element
  template<class T>
  void MeshNumbering<Dimension3, T>::ComputeVariableOrder(int dg_form)
  {
    MeshNumbering_Dim<Dimension3, T>::ComputeVariableOrder(dg_form);
    if (this->type_variable_order != this->CONSTANT_ORDER)
      {	
        // order used for each face
	this->order_face.Reallocate(this->mesh.GetNbBoundary());
        for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
          {
            char r = 1;
            for (int j = 0; j < this->mesh.Boundary(i).GetNbEdges(); j++)
              r = max(r, this->order_edge(this->mesh.Boundary(i).numEdge(j)));
            
            this->order_face(i) = r;
          }
        
        if (dg_form == ElementReference_Base::CONTINUOUS)
          {
            for (int i = 0; i < this->mesh.GetNbElt(); i++)
              {
                char r = this->order_inside(i);
                for (int j = 0; j < this->mesh.Element(i).GetNbEdges(); j++)
                  r = max(r, this->order_edge(this->mesh.Element(i).numEdge(j)));
                
                this->order_element(i) = r;
              }
          }
      }

    if (dg_form == ElementReference_Base::HDG
	&& (this->type_variable_order != this->CONSTANT_ORDER))
      {
	// order for each face is equal to the order of quadrature
	for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
	  this->order_face(i) = this->order_boundary(i);
      }
  }

  
  //! checking numbering of the mesh
  template<class T>
  void MeshNumbering<Dimension3, T>::CheckMeshNumber()
  {
    if (this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    for (int i = 0; i < this->mesh.GetNbEdgesRef(); i++)
      if (this->mesh.GetEdge(i) != this->mesh.EdgeRef(i))
        {
          cout << "Referenced edges not at the beginning of edges" << endl;
          abort();
        }
    
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      if (this->mesh.BoundaryRef(i) != this->mesh.Boundary(i))
        {
          cout << "Referenced faces not at the beginning of faces" << endl;
          abort();
        }

    int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int offset_loc = 0;
        if (nb_dof_vertex > 0)
          for (int j = 0; j < this->mesh.Element(i).GetNbVertices(); j++)
            {
              int nv = this->mesh.Element(i).numVertex(j);
              if (this->OffsetDofVertexNumber(nv) != this->Element(i).GetNumberDof(offset_loc))
                {
                  cout << "OffsetDofVertex incorrect" << endl;

                  abort();
                }
              
              offset_loc += nb_dof_vertex;
            }
        
        int nb_edges_loc = this->mesh.Element(i).GetNbEdges();
        int r = this->GetOrderElement(i);
        int nb_dof_edge_nominal = this->number_map.GetNbDofEdge(r);
        for (int j = 0; j < nb_edges_loc; j++)
          {
            int ne = this->mesh.Element(i).numEdge(j);
            int re = this->GetOrderEdge(ne);
            int nb_dof_edge = this->number_map.GetNbDofEdge(re);
            if (nb_dof_edge > 0)
              {
                int offset = this->Element(i).GetNumberDof(offset_loc);
                for (int k = 1; k < nb_dof_edge; k++)
                  offset = min(offset,
                               this->Element(i).GetNumberDof(offset_loc+k));
                
                if (offset != this->OffsetDofEdgeNumber(ne))
                  {
                    cout << "OffsetDofEdge incorrect" << endl;
                    abort();
                  }
              }
            
            offset_loc += nb_dof_edge_nominal;
          }
        
        int nb_faces_loc = this->mesh.Element(i).GetNbFaces();
        int nb_dof_tri = this->number_map.GetNbDofTriangle(r);
        int nb_dof_quad = this->number_map.GetNbDofQuadrangle(r);
        for (int j = 0; j < nb_faces_loc; j++)
          {
            int nf = this->mesh.Element(i).numFace(j);
            int rf = this->GetOrderFace(nf);
            int nb_dof_face = this->number_map.GetNbDofElement(rf, this->mesh.Boundary(nf));
            if (nb_dof_face > 0)
              {
                int offset = this->Element(i).GetNumberDof(offset_loc);
                for (int k = 1; k < nb_dof_face; k++)
                  offset = min(offset,
                               this->Element(i).GetNumberDof(offset_loc+k));
                
                if (offset != this->OffsetDofFaceNumber(nf))
                  {
                    cout << "OffsetDofFace incorrect" << endl;
                    abort();
                  }
              }
            
            if (this->mesh.Boundary(nf).GetNbVertices() == 3)
              offset_loc += nb_dof_tri;
            else
              offset_loc += nb_dof_quad;
          }
      }
  }
  
  
  //! reconstructs offsets for dofs on edges, faces, elements
  template<class T>
  void MeshNumbering<Dimension3, T>::ReconstructOffsetDofs()
  {
    if (this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    if (this->number_map.FormulationDG() == ElementReference_Base::HDG)
      {
	for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  {
	    int nb_faces_loc = this->mesh.Element(i).GetNbBoundary();
	    int offset_loc = 0;
	    for (int j = 0; j < nb_faces_loc; j++)
	      {
		int nf = this->mesh.Element(i).numBoundary(j);
		int rf = this->GetOrderQuadrature(nf);
		int nb_dof_face = 0;
		if (this->mesh.Boundary(nf).GetNbVertices() == 3)
		  nb_dof_face = this->number_map.GetNbDofTriangle(rf);
		else
		  nb_dof_face = this->number_map.GetNbDofQuadrangle(rf);
		
		this->OffsetDofFaceNumber(nf) = this->Element(i).GetNumberDof(offset_loc);
		for (int k = 1; k < nb_dof_face; k++)
		  this->OffsetDofFaceNumber(nf)
		    = min(this->OffsetDofFaceNumber(nf), this->Element(i).GetNumberDof(offset_loc+k));

		offset_loc += nb_dof_face;
	      }
	  }
	
	return;
      }
    
    int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
    int nb_max_vertex = 0, nb_max_edge = 0, nb_max_face = 0;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int offset_loc = 0;
        if (nb_dof_vertex > 0)
          for (int j = 0; j < this->mesh.Element(i).GetNbVertices(); j++)
            {
              int nv = this->mesh.Element(i).numVertex(j);
              this->OffsetDofVertexNumber(nv) = this->Element(i).GetNumberDof(offset_loc);
              offset_loc += nb_dof_vertex;
              nb_max_vertex = max(nb_max_vertex, this->OffsetDofVertexNumber(nv) + nb_dof_vertex);
            }
      }
        
    nb_max_edge = nb_max_vertex; 
    this->OffsetDofEdgeNumber.Fill(nb_max_vertex);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int offset_loc = nb_dof_vertex*this->mesh.Element(i).GetNbVertices();
        int nb_edges_loc = this->mesh.Element(i).GetNbEdges();
        int r = this->GetOrderElement(i);
        int nb_dof_edge_nominal = this->number_map.GetNbDofEdge(r);
        for (int j = 0; j < nb_edges_loc; j++)
          {
            int ne = this->mesh.Element(i).numEdge(j);
            int re = this->GetOrderEdge(ne);
            int nb_dof_edge = this->number_map.GetNbDofEdge(re);
            if (nb_dof_edge > 0)
              {
                this->OffsetDofEdgeNumber(ne) = this->Element(i).GetNumberDof(offset_loc);
                for (int k = 1; k < nb_dof_edge; k++)
                  this->OffsetDofEdgeNumber(ne) = min(this->OffsetDofEdgeNumber(ne),
						      this->Element(i).GetNumberDof(offset_loc+k));
                
                nb_max_edge = max(nb_max_edge, this->OffsetDofEdgeNumber(ne) + nb_dof_edge);
              }
            
            offset_loc += nb_dof_edge_nominal;
          }
      }
       
    nb_max_face = nb_max_edge;
    this->OffsetDofFaceNumber.Fill(nb_max_edge);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int offset_loc = nb_dof_vertex*this->mesh.Element(i).GetNbVertices();
        int nb_edges_loc = this->mesh.Element(i).GetNbEdges();
        int r = this->GetOrderElement(i);
        int nb_dof_edge_nominal = this->number_map.GetNbDofEdge(r);
        offset_loc += nb_dof_edge_nominal*nb_edges_loc;
        
        int nb_faces_loc = this->mesh.Element(i).GetNbFaces();
        int nb_dof_tri = this->number_map.GetNbDofTriangle(r);
        int nb_dof_quad = this->number_map.GetNbDofQuadrangle(r);
        for (int j = 0; j < nb_faces_loc; j++)
          {
            int nf = this->mesh.Element(i).numFace(j);
            int rf = this->GetOrderFace(nf);
            int nb_dof_face = this->number_map.GetNbDofElement(rf, this->mesh.Boundary(nf));
            if (nb_dof_face > 0)
              {
                this->OffsetDofFaceNumber(nf) = this->Element(i).GetNumberDof(offset_loc);
                for (int k = 1; k < nb_dof_face; k++)
                  this->OffsetDofFaceNumber(nf) = min(this->OffsetDofFaceNumber(nf),
						      this->Element(i).GetNumberDof(offset_loc+k));
                
                nb_max_face = max(nb_max_face, this->OffsetDofFaceNumber(nf) + nb_dof_face);
              }
            
            if (this->mesh.Boundary(nf).GetNbVertices() == 3)
              offset_loc += nb_dof_tri;
            else
              offset_loc += nb_dof_quad;
          }
      }
    
    this->OffsetDofVertexNumber(this->mesh.GetNbVertices()) = nb_max_vertex;
    this->OffsetDofEdgeNumber(this->mesh.GetNbEdges()) = nb_max_edge;
    this->OffsetDofFaceNumber(this->mesh.GetNbFaces()) = nb_max_face;
  }
  
  
  //! performs numbering for faces
  template<class T>
  void MeshNumbering<Dimension3, T>::NumberFaces()
  {
    // we compute the number of dofs to add in each face
    IVect offset_face(this->mesh.GetNbFaces()+1);
    offset_face(0) = 0;
    for (int i = 0; i < this->mesh.GetNbFaces(); i++)
      {
	int r = this->order;
	if (this->type_variable_order != this->CONSTANT_ORDER)
	  r = this->order_face(i);
	
	offset_face(i+1) = offset_face(i)
	  + this->number_map.GetNbDofElement(r, this->mesh.Boundary(i));
      }
        
    // loop over elements
    IVect NegativeDof(5000);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	int r = this->order;
	if (this->type_variable_order != this->CONSTANT_ORDER)
	  r = this->order_element(i);
	
	int nb_vert_loc = this->mesh.Element(i).GetNbVertices();
	int nb_edges_loc = this->mesh.Element(i).GetNbEdges();
	int nb_faces_loc = this->mesh.Element(i).GetNbFaces();
	IVect NegativeDof_ = this->Element(i).GetNegativeDofNumber();
	int nb_negative_dof = NegativeDof_.GetM();
	for (int k = 0; k < nb_negative_dof; k++)
	  NegativeDof(k) = NegativeDof_(k);
	
	int offset_loc = this->number_map.GetNbDofVertex(r)*nb_vert_loc 
	  + this->number_map.GetNbDofEdge(r)*nb_edges_loc;
        
	for (int j = 0; j < nb_faces_loc; j++)
	  {
	    int nf = this->mesh.Element(i).numFace(j);
	    int rf = this->order;
	    if (this->type_variable_order != this->CONSTANT_ORDER)
	      rf = this->order_face(nf);
	    
	    int nb_dof_face = this->number_map.
	      GetNbDofElement(rf, this->mesh.Boundary(nf));
	    
	    int nb_dofF = this->number_map.
	      GetNbDofElement(r, this->mesh.Boundary(nf));
	    
            int rot = this->mesh.Element(i).GetOrientationFace(j);
	    
            if (nb_dof_face > 0)
              {
                if (this->number_map.DofInvariantByRotation(rf, this->mesh.Boundary(nf)))
                  {
                    // dofs invariation by rotation
                    // we can find dof numbers obtained after rotation of the face
                    for (int k = 0; k < nb_dof_face; k++)
                      {
                        int krot = this->number_map.
                          GetRotationFaceDof(rot, rf, k, this->mesh.Boundary(nf));
                        
                        int num_dof_loc = offset_loc + krot;
                        
                        this->Element(i).
                          SetNumberDof(num_dof_loc, this->nodl + offset_face(nf) + k);
                        
                        if (this->number_map.IsNegativeDof(rot, rf, k, this->mesh.Boundary(nf)))
                          {
                            NegativeDof(nb_negative_dof) = num_dof_loc;
                            nb_negative_dof++;
                          }
                      }
                  }
                else
                  {
                    // dofs not invariant by rotation
                    // local dofs are obtained by performing linear combination of global dofs
                    // we only store global dofs (without any renumbering)
                    for (int k = 0; k < nb_dof_face; k++)
                      {                    
                        int num_dof_loc = offset_loc + k;
                        
                        this->Element(i).
                          SetNumberDof(num_dof_loc, this->nodl + offset_face(nf) + k);
                        
                      }
                    
                    // and specifying an order for face
                    this->Element(i).SetOrderFace(j, rf, offset_loc);
                    
                  }
              }
	    
	    offset_loc += nb_dofF;
	  }
	
	// if negative dofs, specifying it
	if (nb_negative_dof > 0)
	  this->Element(i).SetNegativeDofNumber(nb_negative_dof, NegativeDof);
	
      }
    
    this->OffsetDofFaceNumber = offset_face;
    for (int i = 0; i <= this->mesh.GetNbFaces(); i++)
      this->OffsetDofFaceNumber(i) += this->nodl;

    this->nodl += offset_face(this->mesh.GetNbFaces());
  }


  template<class T>
  void MeshNumbering<Dimension3, T>::NumberSurfaceEdges()
  {
    bool variable_order = this->IsOrderVariable();
    int nb_elt = this->mesh.GetNbBoundaryRef();
    IVect NegativeDof(5000);
    if (!variable_order)
      {
	int nb_dof_vertex = this->number_map.GetNbDofVertex(this->order);
	int nb_dof_edge = this->number_map.GetNbDofEdge(this->order);
	this->OffsetDofEdgeNumber.Reallocate(this->mesh.GetNbEdges()+1);
	this->OffsetDofEdgeNumber.Fill();
	Mlt(nb_dof_edge, this->OffsetDofEdgeNumber);
	for (int i = 0; i <= this->mesh.GetNbEdges(); i++)
	  this->OffsetDofEdgeNumber(i) += this->nodl;
        
	for (int i = 0; i < nb_elt; i++)
	  {
	    int nb_vertices_loc = this->mesh.BoundaryRef(i).GetNbVertices();
	    int nb_edges_loc = this->mesh.BoundaryRef(i).GetNbEdges();
	    int offset_loc = nb_vertices_loc*nb_dof_vertex;
	    int nb_negative_dof = 0;
	    for (int j = 0; j < nb_edges_loc; j++)
	      {
		int ne = this->mesh.BoundaryRef(i).numEdge(j);
		int num_vertex = this->mesh.BoundaryRef(i).numVertex(j);
		
		if (this->mesh.GetEdge(ne).numVertex(0) == num_vertex)
		  {
		    // same orientation between local and global edge
		    for (int k = 0; k < nb_dof_edge; k++)
		      this->dof_element(i).SetNumberDof(offset_loc + k,
							this->nodl + ne*nb_dof_edge + k);
		  }
		else
		  {
		    // opposite orientation between local and global edge
		    for (int k = 0; k < nb_dof_edge; k++)
		      {
			int num_dof_loc
			  = offset_loc + this->number_map.GetSymmetricEdgeDof(this->order, k);
			
			this->dof_element(i).SetNumberDof(num_dof_loc, this->nodl + ne*nb_dof_edge + k);
			
			if (this->number_map.IsSkewSymmetricEdgeDof(this->order, k))
			  {
			    NegativeDof(nb_negative_dof) = num_dof_loc;
			    nb_negative_dof++;
			  }
		      }
		  }
		
		offset_loc += nb_dof_edge;
	      }
	    
	    // if negative dofs, specifying it
	    if (nb_negative_dof > 0)
	      this->dof_element(i).SetNegativeDofNumber(nb_negative_dof, NegativeDof);
	  }
	
	this->nodl += nb_dof_edge*this->mesh.GetNbEdges();
      }
    else
      {
	int nb_edges = this->mesh.GetNbEdges();
	IVect offset_edge(nb_edges+1);
	offset_edge(0) = 0;
	for (int i = 0; i < nb_edges; i++)
	  offset_edge(i+1) = offset_edge(i) + this->number_map.GetNbDofEdge(this->order_edge(i));
	
	this->OffsetDofEdgeNumber = offset_edge;
	for (int i = 0; i <= nb_edges; i++)
	  this->OffsetDofEdgeNumber(i) += this->nodl;
	
	for (int i = 0; i < nb_elt; i++)
	  {
	    int nb_vertices_loc = this->mesh.BoundaryRef(i).GetNbVertices();
	    int nb_edges_loc = this->mesh.BoundaryRef(i).GetNbEdges();
	    int r = this->order_element(i);
	    int offset_loc = nb_vertices_loc*this->number_map.GetNbDofVertex(r);
	    int nb_dofE = this->number_map.GetNbDofEdge(r);
	    
	    int nb_negative_dof = 0;
	    for (int j = 0; j < nb_edges_loc; j++)
	      {
		int ne = this->mesh.BoundaryRef(i).numEdge(j);
		int nb_dof_edge = this->number_map.GetNbDofEdge(this->order_edge(ne));
		int num_vertex = this->mesh.BoundaryRef(i).numVertex(j);
		
		if (this->mesh.GetEdge(ne).numVertex(0) == num_vertex)
		  {
		    // same orientation between local and global edge
		    for (int k = 0; k < nb_dof_edge; k++)
		      this->dof_element(i).SetNumberDof(offset_loc + k,
							this->nodl + offset_edge(ne) + k);
		  }
		else
		  {
		    // opposite orientation between local and global edge
		    for (int k = 0; k < nb_dof_edge; k++)
		      {
			int num_dof_loc = offset_loc + this->number_map.
			  GetSymmetricEdgeDof(this->order_edge(ne), k);
			
			this->dof_element(i).SetNumberDof(num_dof_loc,
							  this->nodl + offset_edge(ne) + k);
			
			if (this->number_map.IsSkewSymmetricEdgeDof(this->order_edge(ne), k))
			  {
			    NegativeDof(nb_negative_dof) = num_dof_loc;
			    nb_negative_dof++;
			  }
		      }
		  }
		
		offset_loc += nb_dofE;
	      }
	    
	    // if negative dofs, specifying it
	    if (nb_negative_dof > 0)
	      this->dof_element(i).SetNegativeDofNumber(nb_negative_dof, NegativeDof);
	    
	  }
	    
	this->nodl += offset_edge(nb_edges);
      }
  }
  
  
  //! Numbering of tetrahedral elements
  template<class T>
  void MeshNumbering<Dimension3, T>::
  ConstructTetrahedralNumbering(int order, Array3D<int>& NumNodes3D,
			       Matrix<int>& CoordinateNodes)
  {
    if (order < 0)
      {
	NumNodes3D.Clear();
	CoordinateNodes.Clear();
	return;
      }
    
    NumNodes3D.Reallocate(order+1, order+1, order+1);
    NumNodes3D.Fill(-1);
    CoordinateNodes.Reallocate(NumNodes3D.GetSize(), 4);
    // nodes on the four vertices
    NumNodes3D(0, 0, 0) = 0;     NumNodes3D(order, 0, 0) = 1;
    NumNodes3D(0, order, 0) = 2; NumNodes3D(0, 0, order) = 3;
    
    int node = 0;
    for (int i = 0; i <= order; i += order)
      for (int j = 0; j <= order; j += order)
	for (int k = 0; k <= order; k += order)
	  {
	    node = NumNodes3D(i, j, k);
	    if (node != -1)
	      {
		CoordinateNodes(node, 0) = i;
		CoordinateNodes(node, 1) = j;
		CoordinateNodes(node, 2) = k;
	      }
	  }
    
    // nodes on edges
    node = 4;
    for (int ne = 0; ne < 6; ne++)
      {
 	// we get the two extremities of the edge
	int num_vertex_1 = MeshNumbering<Dimension3, T>::FirstExtremityEdge(0, ne);
	int num_vertex_2 = MeshNumbering<Dimension3, T>::SecondExtremityEdge(0, ne);
	// we get the (i,j,k) coordinate of the first vertex
	int i1 = CoordinateNodes(num_vertex_1, 0);
	int j1 = CoordinateNodes(num_vertex_1, 1);
	int k1 = CoordinateNodes(num_vertex_1, 2);
	// the second vertex
	int i2 = CoordinateNodes(num_vertex_2, 0);
	int j2 = CoordinateNodes(num_vertex_2, 1);
	int k2 = CoordinateNodes(num_vertex_2, 2);
	
	// we are deducing the coordinates for any point on the edge
	for (int l = 1; l < order; l++)
	  {
	    int i = (l*i2 + (order-l)*i1)/order;
	    int j = (l*j2 + (order-l)*j1)/order;
	    int k = (l*k2 + (order-l)*k1)/order;
	    
	    NumNodes3D(i, j, k) = node;
	    node++;
	  }
      }

    // nodes on the faces
    // first face z = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(i, j, 0) = node++;
    
    // second face y = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(i, 0, j) = node++;
    
    // third face x = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(0, i, j) = node++;
    
    // fourth face x+y+z = 1
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(order-i-j, i, j) = node++;
    
    // nodes inside the tetrahedron
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	for (int k = 1; k < order-i-j; k++)
	  NumNodes3D(i,j,k) = node++;

    // CoordinateNodes is generated
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order-i-j; k++)
	  {
	    CoordinateNodes(NumNodes3D(i,j,k), 0) = order-i-j-k;
	    CoordinateNodes(NumNodes3D(i,j,k), 1) = i;
	    CoordinateNodes(NumNodes3D(i,j,k), 2) = j;
	    CoordinateNodes(NumNodes3D(i,j,k), 3) = k;
	  }
  }

  
  //! Numbering of pyramidal elements
  template<class T>
  void MeshNumbering<Dimension3, T>::
  ConstructPyramidalNumbering(int order, Array3D<int>& NumNodes3D,
			      Matrix<int>& CoordinateNodes)
  {
    if (order < 0)
      {
	NumNodes3D.Clear();
	CoordinateNodes.Clear();
	return;
      }

    NumNodes3D.Reallocate(order+1, order+1, order+1);
    NumNodes3D.Fill(-1);
    CoordinateNodes.Reallocate((order+1)*(order+2)*(2*order+3)/6, 3);
    // nodes on the five vertices
    NumNodes3D(0, 0, 0) = 0;     NumNodes3D(order, 0, 0) = 1;
    NumNodes3D(order, order, 0) = 2; NumNodes3D(0, order, 0) = 3;
    NumNodes3D(0, 0, order) = 4;
    
    int node = 0;
    for (int i = 0; i <= order; i += order)
      for (int j = 0; j <= order; j += order)
	for (int k = 0; k <= order; k += order)
	  {
	    node = NumNodes3D(i, j, k);
	    if (node != -1)
	      {
		CoordinateNodes(node, 0) = i;
		CoordinateNodes(node, 1) = j;
		CoordinateNodes(node, 2) = k;
	      }
	  }
    
    // nodes on edges
    node = 5;
    for (int ne = 0; ne < 8; ne++)
      {
	// we get the two extremities of the edge
	int num_vertex_1 = MeshNumbering<Dimension3, T>::FirstExtremityEdge(1, ne);
	int num_vertex_2 = MeshNumbering<Dimension3, T>::SecondExtremityEdge(1, ne);
	// we get the (i,j,k) coordinate of the first vertex
	int i1 = CoordinateNodes(num_vertex_1, 0);
	int j1 = CoordinateNodes(num_vertex_1, 1);
	int k1 = CoordinateNodes(num_vertex_1, 2);
	// the second vertex
	int i2 = CoordinateNodes(num_vertex_2, 0);
	int j2 = CoordinateNodes(num_vertex_2, 1);
	int k2 = CoordinateNodes(num_vertex_2, 2);
	
	// we are deducing the coordinates for any point on the edge
	for (int l = 1; l < order; l++)
	  {
	    int i = (l*i2 + (order-l)*i1)/order;
	    int j = (l*j2 + (order-l)*j1)/order;
	    int k = (l*k2 + (order-l)*k1)/order;
	    
	    NumNodes3D(i, j, k) = node;
	    node++;
	  }
      }
    // DISP(NumNodes3D);
    
    // nodes on the faces
    // first face z = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(i, j, 0) = node++;
    
    // second face y = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(i, 0, j) = node++;
    
    // third face x+z = 1
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(order-j, i, j) = node++;
    
    // fourth face y+z = 1
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(i, order-j, j) = node++;

    // fifth face x = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(0, i, j) = node++;
    
    // nodes inside the pyramid
    for (int k = 1; k < order; k++)
      for (int i = 1; i < order-k; i++)
	for (int j = 1; j < order-k; j++)
	  NumNodes3D(i, j, k) = node++;
    
    // DISP(NumNodes3D);
    // CoordinateNodes is generated
    for (int k = 0; k <= order; k++)
      for (int i = 0; i <= order-k; i++)
	for (int j = 0; j <= order-k; j++)
	  {
	    CoordinateNodes(NumNodes3D(i,j,k), 0) = i;
	    CoordinateNodes(NumNodes3D(i,j,k), 1) = j;
	    CoordinateNodes(NumNodes3D(i,j,k), 2) = k;
	  }
  }
  
  
  //! Numbering of tetrahedral elements
  template<class T>
  void MeshNumbering<Dimension3, T>::
  ConstructPrismaticNumbering(int order, Array3D<int>& NumNodes3D,
			      Matrix<int>& CoordinateNodes)
  {
    if (order < 0)
      {
	NumNodes3D.Clear();
	CoordinateNodes.Clear();
	return;
      }
    
    NumNodes3D.Reallocate(order+1, order+1, order+1);
    NumNodes3D.Fill(-1);
    CoordinateNodes.Reallocate(NumNodes3D.GetSize(), 4);
    // nodes on the six vertices
    NumNodes3D(0, 0, 0) = 0;     NumNodes3D(order, 0, 0) = 1;
    NumNodes3D(0, order, 0) = 2; NumNodes3D(0, 0, order) = 3;
    NumNodes3D(order, 0, order) = 4; NumNodes3D(0, order, order) = 5;
    
    int node = 0;
    for (int i = 0; i <= order; i += order)
      for (int j = 0; j <= order; j += order)
	for (int k = 0; k <= order; k += order)
	  {
	    node = NumNodes3D(i, j, k);
	    if (node != -1)
	      {
		CoordinateNodes(node, 0) = i;
		CoordinateNodes(node, 1) = j;
		CoordinateNodes(node, 2) = k;
	      }
	  }
    
    // nodes on edges
    node = 6;
    for (int ne = 0; ne < 9; ne++)
      {
	// we get the two extremities of the edge
	int num_vertex_1 = MeshNumbering<Dimension3, T>::FirstExtremityEdge(2, ne);
	int num_vertex_2 = MeshNumbering<Dimension3, T>::SecondExtremityEdge(2, ne);
	// we get the (i,j,k) coordinate of the first vertex
	int i1 = CoordinateNodes(num_vertex_1, 0);
	int j1 = CoordinateNodes(num_vertex_1, 1);
	int k1 = CoordinateNodes(num_vertex_1, 2);
	// the second vertex
	int i2 = CoordinateNodes(num_vertex_2, 0);
	int j2 = CoordinateNodes(num_vertex_2, 1);
	int k2 = CoordinateNodes(num_vertex_2, 2);
	
	// we are deducing the coordinates for any point on the edge
	for (int l = 1; l < order; l++)
	  {
	    int i = (l*i2 + (order-l)*i1)/order;
	    int j = (l*j2 + (order-l)*j1)/order;
	    int k = (l*k2 + (order-l)*k1)/order;
	    
	    NumNodes3D(i, j, k) = node;
	    node++;
	  }
      }

    // nodes on the faces
    // first face z = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(i, j, 0) = node++;
    
    // second face y = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(i, 0, j) = node++;

    // third face x+y = 1
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(order-i, i, j) = node++;
    
    // fourth face x = 0
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(0, i, j) = node++;
    
    // fifth face z = 1
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	NumNodes3D(i, j, order) = node++;
    
    // nodes inside the wedge
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order-i; j++)
	for (int k = 1; k < order; k++)
	  NumNodes3D(i, j, k) = node++;

    // CoordinateNodes is generated
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  {
	    CoordinateNodes(NumNodes3D(i,j,k), 0) = i;
	    CoordinateNodes(NumNodes3D(i,j,k), 1) = j;
	    CoordinateNodes(NumNodes3D(i,j,k), 2) = k;
	  }
  }

  
  //! Numbering of hexahedral elements
  template<class T>
  void MeshNumbering<Dimension3, T>::
  ConstructHexahedralNumbering(int order, Array3D<int>& NumNodes3D,
			       Matrix<int>& CoordinateNodes)
  {
    HexahedronQuadrature::ConstructHexahedralNumbering(order, NumNodes3D, CoordinateNodes);
  }

  
  //! computes new numbers after rotation and/or symmetry of the triangle
  template<class T>
  void MeshNumbering<Dimension3, T>::
  GetRotationTriangularFace(const Matrix<int>& NumNodes2D,
			    Matrix<int>& FacesTri_Rotation)
  {
    int r = NumNodes2D.GetM();
    FacesTri_Rotation.Reallocate(6, (r+1));
    FacesTri_Rotation.Fill(-1);
    // rotation of 120° : (x,y) -> (y, 1-x)
    //                   (i,j) -> (j, r-2-i)
    // rotation of 240°  (i,j) -> (r-2-i, r-2-j)
    // symmetry          (i,j) -> (j,i)
    // rotation of 120° and symmetry (i,j) -> (r-2-i,j)
    // rotation of 240° and symmetry (i,j) -> (r-2-j,r-2-i)
    int nb = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	{
	  // i -> lambda 2  j-> lambda3
	  int k = r - 1 - i - j;
	  FacesTri_Rotation(0, NumNodes2D(i,j)) = nb;
	  FacesTri_Rotation(1, NumNodes2D(k,i)) = nb;
	  FacesTri_Rotation(2, NumNodes2D(j,k)) = nb;
	  FacesTri_Rotation(3, NumNodes2D(j,i)) = nb;
	  FacesTri_Rotation(4, NumNodes2D(k,j)) = nb;
	  FacesTri_Rotation(5, NumNodes2D(i,k)) = nb;
	  nb++;
	}
  }
  
  
  //! computes new numbers after rotation and/or symmetry of the triangle
  template<class T>
  void MeshNumbering<Dimension3, T>::
  GetRotationTriangularFace(const VectR2& Points2D,
			    Matrix<int>& FacesTri_Rotation)
  {
    int nb_points = Points2D.GetM();
    FacesTri_Rotation.Reallocate(6, nb_points);
    for (int i = 0; i < nb_points; i++)
      {
	T x = Points2D(i)(0), y = Points2D(i)(1);
	// loop over all quadrature points to find the point after rotation
	int rot0 = -1, rot1 = -1, rot2 = -1, rot3 = -1, rot4 = -1;
	T epsilon = 1e-10;
	for (int j = 0; j < nb_points; j++)
	  {
	    T li = Points2D(j)(0), lj = Points2D(j)(1), lk = 1.0-li-lj;
	    if ((abs(lk-x) < epsilon)&&(abs(li-y) < epsilon))
	      rot0 = j;
	    
	    if ((abs(lj-x) < epsilon)&&(abs(lk-y) < epsilon))
	      rot1 = j;
	    
	    if ((abs(lj-x) < epsilon)&&(abs(li-y) < epsilon))
	      rot2 = j;
	    
	    if ((abs(lk-x) < epsilon)&&(abs(lj-y) < epsilon))
	      rot3 = j;
	    
	    if ((abs(li-x) < epsilon)&&(abs(lk-y) < epsilon))
	      rot4 = j;
	  }
	FacesTri_Rotation(0, i) = i;
	FacesTri_Rotation(1, i) = rot0;
	FacesTri_Rotation(2, i) = rot1;
        FacesTri_Rotation(3, i) = rot2;
	FacesTri_Rotation(4, i) = rot3;
	FacesTri_Rotation(5, i) = rot4;
      }
  }
  
  
  //! computes new numbers after rotation and/or symmetry of the square
  template<class T>
  void MeshNumbering<Dimension3, T>::
  GetRotationQuadrilateralFace(const Matrix<int>& NumNodes2D,
			       Matrix<int>& FacesQuad_Rotation)
  {
    int r = NumNodes2D.GetM(); 
    FacesQuad_Rotation.Reallocate(8, r*r);
    FacesQuad_Rotation.Fill(-1);
    // rotation of 90° : (x,y) -> (y,1-x)
    //                   (i,j) -> (j, r-1-i)
    // rotation of 180°  (i,j) -> (r-1-i, r-1-j)
    // rotation of 270°  (i,j) -> (r-1-j, i)
    // symmetry          (i,j) -> (j, i)
    // rotation of 90° and symmetry (i, j) -> (r-1-i,j)
    // rotation of 180° and symmetry (i, j) -> (r-1-j, r-1-i)
    // rotation of 270° and symmetry (i, j) -> (i, r-1-j)
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{
	  int node = NumNodes2D(i,j);
	  FacesQuad_Rotation(0, node) = NumNodes2D(i,j);
	  FacesQuad_Rotation(1, node) = NumNodes2D(j, r-1-i);
	  FacesQuad_Rotation(2, node) = NumNodes2D(r-1-i, r-1-j);
	  FacesQuad_Rotation(3, node) = NumNodes2D(r-1-j, i);
	  FacesQuad_Rotation(4, node) = NumNodes2D(j, i);
	  FacesQuad_Rotation(5, node) = NumNodes2D(r-1-i, j);
	  FacesQuad_Rotation(6, node) = NumNodes2D(r-1-j, r-1-i);
	  FacesQuad_Rotation(7, node) = NumNodes2D(i, r-1-j);
	}
    
  }
  
  template<class T>
  void MeshNumbering<Dimension3, T>::
  GetRotationQuadrilateralFace(const VectR2& Points2D,
                               Matrix<int>& FacesQuad_Rotation)
  {
    int nb_points = Points2D.GetM();
    FacesQuad_Rotation.Reallocate(8, nb_points);
    for (int i = 0; i < nb_points; i++)
      {
	T x = Points2D(i)(0), y = Points2D(i)(1);
	// loop over all points to find the point after rotation
	int rot1 = -1, rot2 = -1, rot3 = -1, rot4 = -1, rot5 = -1, rot6 = -1, rot7 = -1;
	T epsilon = 1e-10;
	for (int j = 0; j < nb_points; j++)
	  {
	    T xp = Points2D(j)(0), yp = Points2D(j)(1);
            // x' = 1-y   y'= x
	    if ((abs(1.0-yp - x) < epsilon)&&(abs(xp-y) < epsilon))
	      rot1 = j;
	    
            // x' = 1-x  y' = 1-y
	    if ((abs(1.0-xp - x) < epsilon)&&(abs(1.0-yp-y) < epsilon))
	      rot2 = j;
	    
            // x' = y   y' = 1-x
	    if ((abs(yp-x) < epsilon)&&(abs(1.0-xp-y) < epsilon))
	      rot3 = j;
	    
            // x' = y  y' = x
	    if ((abs(yp-x) < epsilon)&&(abs(xp-y) < epsilon))
	      rot4 = j;
	    
            // x' = 1-x  y' = y
	    if ((abs(1.0-xp-x) < epsilon)&&(abs(yp-y) < epsilon))
	      rot5 = j;

            // x' = 1-y  y' = 1-x
	    if ((abs(1.0-yp-x) < epsilon)&&(abs(1.0-xp-y) < epsilon))
	      rot6 = j;
            
            // x' = x  y' = 1-y
	    if ((abs(xp-x) < epsilon)&&(abs(1.0-yp-y) < epsilon))
	      rot7 = j;
	  }
        
	FacesQuad_Rotation(0, i) = i;
	FacesQuad_Rotation(1, i) = rot1;
	FacesQuad_Rotation(2, i) = rot2;
        FacesQuad_Rotation(3, i) = rot3;
	FacesQuad_Rotation(4, i) = rot4;
	FacesQuad_Rotation(5, i) = rot5;
        FacesQuad_Rotation(6, i) = rot6;
        FacesQuad_Rotation(7, i) = rot7;
      }
  }  
#endif
  

#ifdef MONTJOIE_WITH_THREE_DIM
  //! treatment of periodic boundaries
  template<class T>
  void MeshNumbering<Dimension3, T>::TreatPeriodicCondition()
  {
    int nb_boundaries_ref = this->mesh.GetNbBoundaryRef();
    
    // matching boundaries
    // PeriodicityBoundary(i) will be the global number of the translated boundary 
    // to the global boundary i
    this->PeriodicityBoundary.Reallocate(nb_boundaries_ref);
    this->PeriodicityBoundary.Fill(-1);
    this->type_periodicity_boundary.Reallocate(nb_boundaries_ref);
    this->type_periodicity_boundary.Fill(0);
    this->TranslationPeriodicBoundary.Reallocate(nb_boundaries_ref);
    if (this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      this->type_formulation_periodic = this->WEAK_PERIODIC;

    if (this->number_map.FormulationDG() == ElementReference_Base::HDG)
      if (this->type_formulation_periodic == this->WEAK_PERIODIC)
	this->type_formulation_periodic = this->SAME_PERIODIC_DOFS;


    if (this->mesh.GetNbPeriodicReferences() <= 0)
      return;
    
    int nodl = this->nodl;
    int nb_vertices = this->mesh.GetNbVertices();
    
    IVect PeriodicityDof(nodl);
    IVect DofUsed(nodl); DofUsed.Fill(false);
    
    int nb_dof_periodic = 0;
    // evaluating the number of periodic dofs (upper bound)
    IVect ref_cond(this->mesh.GetNbReferences()+1);
    ref_cond.Fill(0);
    for (int num_per = 0; num_per < this->mesh.GetNbPeriodicReferences(); num_per++)
      ref_cond(this->mesh.GetPeriodicReference(num_per,1)) = 1;
    
    for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
      {
        int ref = this->mesh.Boundary(i).GetReference();
        if ((ref > 0) && ref_cond(ref) == 1)
          {
            int nb_vert = this->mesh.Boundary(i).GetNbVertices();
            int r = this->GetOrder();
            for (int j = 0; j < nb_vert; j++)
              {
                //int nv = this->mesh.Boundary(i).numVertex(j);
                nb_dof_periodic += this->number_map.GetNbDofVertex(r);
                
                int ne = this->mesh.Boundary(i).numEdge(j);
                if (this->IsOrderVariable())
                  r = this->order_edge(ne);
                
                nb_dof_periodic += this->number_map.GetNbDofEdge(r);
              }
            
            r = this->GetOrderFace(i);
            if (nb_vert == 3)
              nb_dof_periodic += this->number_map.GetNbDofTriangle(r);
            else
            nb_dof_periodic += this->number_map.GetNbDofQuadrangle(r);
          }
      }
    
    // allocating arrays with that upper bound
    this->TranslationPeriodicDof.Reallocate(nb_dof_periodic);
    this->PeriodicDofs.Reallocate(nb_dof_periodic);
    this->PeriodicDofs_Original.Reallocate(nb_dof_periodic);
    this->type_periodicity_dof.Reallocate(nb_dof_periodic);

    nb_dof_periodic = 0;
    // loop over all periodic boundaries specified by the user
    for (int num_per = 0; num_per < this->mesh.GetNbPeriodicReferences(); num_per++)
      {
	// references of the two corresponding faces
	int ref1 = this->mesh.GetPeriodicReference(num_per, 0);
	int ref2 = this->mesh.GetPeriodicReference(num_per, 1);
	
	// we get all the faces of reference 1, and of reference 2
	VectBool IsFace_Ref1(nb_boundaries_ref), IsFace_Ref2(nb_boundaries_ref);
	IsFace_Ref1.Fill(false);  IsFace_Ref2.Fill(false);
	int nb_faces_ref1 = 0, nb_faces_ref2 = 0;
	for (int i = 0; i < nb_boundaries_ref; i++)
	  {
	    if (this->mesh.BoundaryRef(i).GetReference() == ref1)
	      {
		nb_faces_ref1++;
		IsFace_Ref1(i) = true;
	      }
	    
	    if (this->mesh.BoundaryRef(i).GetReference() == ref2)
	      {
		nb_faces_ref2++;
		IsFace_Ref2(i) = true;
	      }
	  }

	//  computation of centers of faces
	VectR3 barycenter_faces1(nb_faces_ref1), barycenter_faces2(nb_faces_ref2);

	// index ref 1 :  face numbers related to first reference
	// index ref 2 :  face numbers related to second reference
	IVect Index_ref1(nb_faces_ref1), Index_ref2(nb_faces_ref2);
	R3 center_ref1, center_ref2;
	nb_faces_ref1 = 0; nb_faces_ref2 = 0;
	
	// we fill TypeVertex 
	// TypeVertex(i) = 1 if vertex i belongs to ref1
	// TypeVertex(i) = 2 if vertex i belongs to ref2
	// TypeVertex(i) = 3 if vertex i belongs to both boundaries
	// we get all the vertices of reference 1 and 2
	IVect TypeVertex(nb_vertices); TypeVertex.Fill(-1);
	for (int i = 0; i < nb_boundaries_ref; i++)
	  {
	    if (this->mesh.BoundaryRef(i).GetReference() == ref1)
	      {
		Index_ref1(nb_faces_ref1) = i;
	
		for (int j = 0; j < this->mesh.BoundaryRef(i).GetNbVertices(); j++)
		  {
                    int nv = this->mesh.BoundaryRef(i).numVertex(j);
                    barycenter_faces1(nb_faces_ref1) += this->mesh.Vertex(nv);
		    // for cyclic periodicity, the vertices on the axis are both in ref1 and ref2
		    if (TypeVertex(nv) == 2)
		      TypeVertex(nv) = 3;
		    else if (TypeVertex(nv) == -1)
                      TypeVertex(nv) = 1;
		  }
		
		Mlt(T(1)/this->mesh.BoundaryRef(i).GetNbVertices(),
                    barycenter_faces1(nb_faces_ref1));
		center_ref1 += barycenter_faces1(nb_faces_ref1);
		nb_faces_ref1++;
	      }
	    
	    if (this->mesh.BoundaryRef(i).GetReference() == ref2)
	      {
		Index_ref2(nb_faces_ref2) = i;
	
		for (int j = 0; j < this->mesh.BoundaryRef(i).GetNbVertices(); j++)
		  {
                    int nv = this->mesh.BoundaryRef(i).numVertex(j);
                    barycenter_faces2(nb_faces_ref2) += this->mesh.Vertex(nv);
		    // for cyclic periodicity, the vertices on the axis are both 
                    // in ref1 and ref2 -> TypeVertex = 3
		    if (TypeVertex(nv) == 1)
		      TypeVertex(nv) = 3;
		    else if (TypeVertex(nv) == -1)
		      TypeVertex(nv) = 2;
		  }
		
		Mlt(T(1)/this->mesh.BoundaryRef(i).GetNbVertices(),
                    barycenter_faces2(nb_faces_ref2));

                center_ref2 += barycenter_faces2(nb_faces_ref2);
		nb_faces_ref2++;
	      }
	  }
	
	// center_ref1 and center_ref2 are the centers of the two boundaries ref1 and ref2
	Mlt(T(1.0/nb_faces_ref1), center_ref1);
	Mlt(T(1.0/nb_faces_ref2), center_ref2);
        
	// number of vertices on each boundary
	int nb_vertices_ref1(0), nb_vertices_ref2(0);
	for (int i = 0; i < nb_vertices; i++)
	  {
	    if (TypeVertex(i) == 1)
	      nb_vertices_ref1++;
	    else if (TypeVertex(i) == 2)
	      nb_vertices_ref2++;
	  }
	
	// vertices of reference 1 and 2, construction of indexing arrays
	VectR3 Vertices_Ref1(nb_vertices_ref1);
	VectR3 Vertices_Ref2(nb_vertices_ref2);
	
	// Index_Vert1(i) : vertex number of the vertex i belonging to ref 1
        // Index_Vert2(i) : vertex number of the vertex i belonging to ref 2
	// Vertices_Ref1(i)  : coordinates 
	// Vertices_Ref2(i) : coordinates
	IVect Index_Vert1(nb_vertices_ref1), Index_Vert2(nb_vertices_ref2);
	
	nb_vertices_ref1 = 0; nb_vertices_ref2 = 0;
	for (int i = 0; i < nb_vertices; i++)
	  {
	    // vertices belonging to boundary ref1 are retrieved
	    // vertex numbers are in Index_Vert1, coordinates in Vertices_Ref1
	    if (TypeVertex(i) == 1)
	      {
		Index_Vert1(nb_vertices_ref1) = i;
		Vertices_Ref1(nb_vertices_ref1) = this->mesh.Vertex(i);
	       
		nb_vertices_ref1++;
	      }
	    
	    // vertices belonging to boundary ref2 are retrieved
	    // vertex numbers are in Index_Vert2, coordinates in Vertices_Ref2
	    if (TypeVertex(i) == 2)
	      {
		Index_Vert2(nb_vertices_ref2) = i;
		Vertices_Ref2(nb_vertices_ref2) = this->mesh.Vertex(i);
		nb_vertices_ref2++;
	      }
	  }

	// translation vector between the two boundaries ?
	R3 vh_translat;
        int type_per = this->mesh.GetPeriodicityTypeReference(num_per);
	if (type_per == BoundaryConditionEnum::PERIODIC_THETA) 
	  {
	    // cyclic case
            vh_translat(0) = this->mesh.GetPeriodicAlpha();
	  }
	else 
	  {
	    // cartesian coordinates
	    vh_translat = center_ref2 - center_ref1;
	    
	    // all centers of faces of ref2 are translated
	    for (int i = 0; i < nb_faces_ref2; i++)
	      barycenter_faces2(i) -= vh_translat;
	    
	  }

        if ((type_per != BoundaryConditionEnum::PERIODIC_CTE)
            && !this->plane_wave_quasi_periodic)
          { 
            // storing 0 or 1 for translations
            if (abs(vh_translat(0)) < R3::threshold)
              vh_translat(0) = 0.0;
            else if (vh_translat(0) < 0)
              vh_translat(0) = -1.0;
            else if (vh_translat(0) > 0)
              vh_translat(0) = 1.0;
            
            if (abs(vh_translat(1)) < R3::threshold)
              vh_translat(1) = 0.0;
            else if (vh_translat(1) < 0)
              vh_translat(1) = -1.0;
            else if (vh_translat(1) > 0)
              vh_translat(1) = 1.0;
            
            if (abs(vh_translat(2)) < R3::threshold)
              vh_translat(2) = 0.0;
            else if (vh_translat(2) < 0)
              vh_translat(2) = -1.0;
            else if (vh_translat(2) > 0)
              vh_translat(2) = 1.0;
          }
	
	// all centers of faces are sorted, so that we can compare them
	Sort(nb_faces_ref1, barycenter_faces1, Index_ref1);
	Sort(nb_faces_ref2, barycenter_faces2, Index_ref2);
	
	if (type_per == BoundaryConditionEnum::PERIODIC_THETA) 
	  {
	    // vertices of ref2 are rotated in order to compare them
	    // with vertices of ref1
	    for (int i = 0; i < nb_vertices_ref2; i++)
	      {
		R3 temp = Vertices_Ref2(i) ; 
		// rotation de - theta
		Vertices_Ref2(i)(0)
                  = cos( -vh_translat(0) ) * temp(0) - sin(- vh_translat(0)) *  temp(1);
		Vertices_Ref2(i)(1)
                  = sin(-vh_translat(0) ) * temp(0) + cos(-vh_translat(0) ) *  temp(1);
	      }
	  }
	else
	  {
	    // all vertices of ref2 are translated
	    for (int i = 0; i < nb_vertices_ref2; i++)
	      Vertices_Ref2(i) -= vh_translat;
            
	  }
	
	// vertices are sorted to have direct matching
	Sort(nb_vertices_ref1, Vertices_Ref1, Index_Vert1);
	Sort(nb_vertices_ref2, Vertices_Ref2, Index_Vert2);
	
	IVect InverseIndex_Vert(nb_vertices);
	// InverseIndex_Vert(i) -> numero local (parmi tous les sommets de ref1 ) 
        // du sommet global i
	for (int i = 0; i < nb_vertices_ref1; i++)
	  InverseIndex_Vert(Index_Vert1(i)) = i;
	
	// InverseIndex_Vert(i) -> numero local (parmi tous les sommets de ref2 ) 
        // du sommet global i
	for (int i = 0; i < nb_vertices_ref2; i++)
	  InverseIndex_Vert(Index_Vert2(i)) = i;
	
	PeriodicityDof.Fill(-1);
	// degrees of freedom associated with vertices are treated
	int nb_dof_vertex = this->number_map.GetNbDofVertex(this->GetOrder());
	int offset_vertices = 0;
	if (nb_dof_vertex > 0)
	  {
	    if (nb_vertices_ref1 != nb_vertices_ref2)
	      {
		cout << "the mesh should be the same between periodic boundaries" << endl;
		abort();
	      }
	    
	    for (int i = 0; i < nb_vertices_ref1; i++)
	      {
		int n1 = Index_Vert1(i);
		int n2 = Index_Vert2(i);
		for (int k = 0; k < nb_dof_vertex; k++)
		  PeriodicityDof(n2*nb_dof_vertex + k) = n1*nb_dof_vertex + k;
	      }
	    
	    offset_vertices += nb_dof_vertex*nb_vertices;
	  }
	
	// degrees of freedom associated with edges are treated
        if (this->number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
          {
            Vector<bool> EdgeInRef1(this->mesh.GetNbEdges());
            Vector<bool> EdgeInRef2(this->mesh.GetNbEdges());
            IVect Index_Edge1(4*nb_faces_ref1);
            IVect Index_Edge2(4*nb_faces_ref1);
            EdgeInRef1.Fill(false);
            EdgeInRef2.Fill(false);
            int nb_edges_ref1 = 0, nb_edges_ref2 = 0;
            for (int i = 0; i < nb_faces_ref1; i++)
              {
                int nf = Index_ref1(i);
                for (int j = 0; j < this->mesh.BoundaryRef(nf).GetNbVertices(); j++)
                  {
                    int ne = this->mesh.BoundaryRef(nf).numEdge(j);
                    if (!EdgeInRef1(ne))
                      {
                        EdgeInRef1(ne) = true;
                        Index_Edge1(nb_edges_ref1) = ne;
                        nb_edges_ref1++;
                      }
                  }
              }

            for (int i = 0; i < nb_faces_ref1; i++)
              {
                int nf = Index_ref2(i);
                for (int j = 0; j < this->mesh.BoundaryRef(nf).GetNbVertices(); j++)
                  {
                    int ne = this->mesh.BoundaryRef(nf).numEdge(j);
                    if (!EdgeInRef2(ne))
                      {
                        EdgeInRef2(ne) = true;
                        Index_Edge2(nb_edges_ref2) = ne;
                        nb_edges_ref2++;
                      }
                  }
              }
            
            if (nb_edges_ref1 != nb_edges_ref2)
              {
                cout << "Periodic faces should have the same mesh " << endl;
                abort();
              }
            
            for (int i = 0; i < nb_edges_ref1; i++)
              {
                int n1 = Index_Edge1(i);
                int n2 = Index_Edge2(i);
                int r = this->GetOrderEdge(n1);
                int nb_dof_edge = this->number_map.GetNbDofEdge(r);	    
                int offset1 = this->OffsetDofEdgeNumber(n1);
                int offset2 = this->OffsetDofEdgeNumber(n2);
                
                for (int k = 0; k < nb_dof_edge; k++)
                  PeriodicityDof(offset2 + k) = offset1 + k;	    
              }
            
            offset_vertices += this->OffsetDofEdgeNumber(this->mesh.GetNbEdges());
          }
	

	if (nb_faces_ref1 != nb_faces_ref2)
	  {
	    cout << "Periodic faces should have the same mesh " << endl;
	    abort();
	  }
	
	for (int i = 0; i < nb_faces_ref1; i++)
	  {
	    int n1 = Index_ref1(i);
	    int n2 = Index_ref2(i);
            for (int k = 0; k < this->mesh.BoundaryRef(n1).GetNbVertices(); k++)
              {
                int nv1 = InverseIndex_Vert(this->mesh.BoundaryRef(n1).numVertex(k));
                int nv2 = InverseIndex_Vert(this->mesh.BoundaryRef(n2).numVertex(k));
                if (nv1 != nv2)
                  {
                    cout << "Periodic boundaries must be the same " << endl;
                    this->mesh.Write("toto.mesh");
                    abort();
                  }
              }
                
	    int r = this->GetOrderFace(n1);
	    int nb_dof_face = 0;
	    if (this->mesh.BoundaryRef(n1).GetNbVertices() == 3)
	      nb_dof_face = this->number_map.GetNbDofTriangle(r);    
	    else
	      nb_dof_face = this->number_map.GetNbDofQuadrangle(r);
	    
	    if (this->number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
              {
                int offset1 = this->OffsetDofFaceNumber(n1);
                int offset2 = this->OffsetDofFaceNumber(n2);
                
                for (int k = 0; k < nb_dof_face; k++)
                  PeriodicityDof(offset2 + k) = offset1 + k;	    
              }
            
	    int num_elem1 = this->mesh.Boundary(n1).numElement(0);
            int num_elem2 = this->mesh.Boundary(n2).numElement(0);
	    
	    this->mesh.Boundary(n1).AddElement(num_elem2);
            this->mesh.Boundary(n2).AddElement(num_elem1);
	    
            this->PeriodicityBoundary(n1) = n2;
            this->PeriodicityBoundary(n2) = n1;
            this->type_periodicity_boundary(n1) = type_per;
            this->type_periodicity_boundary(n2) = type_per;
            this->TranslationPeriodicBoundary(n2) = vh_translat;
            this->TranslationPeriodicBoundary(n1) = -vh_translat;
	  }
	
	// pour tous les ddl periodiques
	for (int i = 0; i < nodl; i++)
	  if (i != PeriodicityDof(i))
	    if ((PeriodicityDof(i) >= 0)&&(!DofUsed(i)))
	      {
		this->TranslationPeriodicDof(nb_dof_periodic) = vh_translat;
		this->PeriodicDofs(nb_dof_periodic) = i;
		this->PeriodicDofs_Original(nb_dof_periodic) = PeriodicityDof(i);
		this->type_periodicity_dof(nb_dof_periodic) = type_per;
                DofUsed(i) = true;
		nb_dof_periodic++;
	      }
      }
        
    this->PeriodicDofs.Resize(nb_dof_periodic);
    this->PeriodicDofs_Original.Resize(nb_dof_periodic);
    this->TranslationPeriodicDof.Resize(nb_dof_periodic);
    this->type_periodicity_dof.Resize(nb_dof_periodic);

    // we are searching original dofs
    IVect IndexPeriodicDof(nodl);
    IndexPeriodicDof.Fill(-1);
    for (int i = 0; i < nb_dof_periodic; i++)
      IndexPeriodicDof(this->PeriodicDofs(i)) = i;
    
    for (int i = 0; i < nb_dof_periodic; i++)
      {
        R3 v1 = this->TranslationPeriodicDof(i), v2, v3;
        int n1 = this->PeriodicDofs_Original(i);
        int n2 = IndexPeriodicDof(n1);
        int t1 = this->type_periodicity_dof(i), t2 = -1, t3 = -1;
        if (n2 >= 0)
          {
            t2 = this->type_periodicity_dof(n2);
            v2 = this->TranslationPeriodicDof(n2);
            n1 = this->PeriodicDofs_Original(n2);
            n2 = IndexPeriodicDof(n1);
            if (n2 >= 0)
              {
                t3 = this->type_periodicity_dof(n2);
                v3 = this->TranslationPeriodicDof(n2);
                n1 = this->PeriodicDofs_Original(n2);
                n2 = IndexPeriodicDof(n1);
                if (n2 >= 0)
                  {
                    abort();
                  }
              }
          }
        
        this->PeriodicDofs_Original(i) = n1;
        this->type_periodicity_dof(i) = BoundaryConditionEnum::GetCompositionPeriodicity(t1, t2, t3);
        this->TranslationPeriodicDof(i) = v1 + v2 + v3;
      }
    
    if ((this->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
        || this->type_formulation_periodic == this->WEAK_PERIODIC)
      {
        this->TranslationPeriodicDof.Clear();
        this->PeriodicDofs.Clear();
        this->PeriodicDofs_Original.Clear();
        this->type_periodicity_dof.Clear();        
      }
    else if (this->UseSameNumberForPeriodicDofs())
      {
	// we associate same dof numbers for periodic dofs
	// total number of dofs is thus decreased
	IVect NewNumber(nodl); NewNumber.Fill(-1); int nb = 0;
	for (int i = 0; i < nodl; i++)
	  if (!DofUsed(i))
	    NewNumber(i) = nb++;
        
	for (int i = 0; i < nb_dof_periodic; i++)
	  {
            int n1 = this->PeriodicDofs_Original(i);
	    NewNumber(this->PeriodicDofs(i)) = NewNumber(n1);
          }
        
        // replacing dof numbers
	for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  {
            int nb_dof = this->Element(i).GetNbDof();
            
            for (int j = 0; j < nb_dof; j++)
              {
                int num_dof = this->Element(i).GetNumberDof(j);
                this->Element(i).SetNumberDof(j, NewNumber(num_dof));
              }
            
          }
	
        this->nodl = nb;
        this->ReconstructOffsetDofs();
	nb_dof_periodic = 0;
        this->PeriodicDofs.Clear();
        this->PeriodicDofs_Original.Clear();
        this->TranslationPeriodicDof.Clear();
        this->type_periodicity_dof.Clear();
      }
    else
      {
        //this->PeriodicityBoundary.Fill(-1);
        //this->TranslationPeriodicBoundary.Clear();
      }
    
  }
#endif
  
} // namespace Montjoie

#define MONTJOIE_FILE_NUMBER_MESH_CXX
#endif
