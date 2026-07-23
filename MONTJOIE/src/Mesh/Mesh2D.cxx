#ifndef MONTJOIE_FILE_MESH2D_CXX

namespace Montjoie
{

  //! copy constructor
  template<class T>
  Mesh<Dimension2, T>::Mesh(const Mesh<Dimension2, T>& mesh)
    : Mesh_Base<Dimension2, T>(mesh), MeshBoundaries<Dimension2, T>(mesh),
      triangle_reference(mesh.triangle_reference),
      quadrangle_reference(mesh.quadrangle_reference)
  {
    UpdateReferenceElement();
  }

  
  //! Copies a mesh
  template<class T>
  Mesh<Dimension2, T>& Mesh<Dimension2, T>
  ::operator=(const Mesh<Dimension2, T>& mesh)
  {
    static_cast<Mesh_Base<Dimension2, T>& >(*this)
      = static_cast<const Mesh_Base<Dimension2, T>& >(mesh);
    
    static_cast<MeshBoundaries<Dimension2, T>& >(*this)
      = static_cast<const MeshBoundaries<Dimension2, T>& >(mesh);
    
    triangle_reference = mesh.triangle_reference;
    quadrangle_reference = mesh.quadrangle_reference;
    UpdateReferenceElement();
    
    return *this;
  }

  
  //! updates pointers stored in reference_element
  template<class T>
  void Mesh<Dimension2, T>::UpdateReferenceElement()
  {
    for (int i = 0; i < this->reference_element.GetM(); i++)
      {
        if (this->Element(i).GetNbVertices() == 3)
          this->reference_element(i) = &triangle_reference;
        else
          this->reference_element(i) = &quadrangle_reference;
      }
  }


  //! returns node numbers of the edge num_loc
  template<class T>
  IVect Mesh<Dimension2, T>::GetLocalNodalNumber(int type, int num_loc)
  {
    IVect num;
    switch (type)
      {
      case 0 :
        {
          num.Reallocate(triangle_reference.GetNbNodalBoundary(num_loc));
          for (int j = 0; j < num.GetM(); j++)
            num(j) = triangle_reference.GetNodalNumber(num_loc, j);
        }
        break;
      case 1 :
        {
          num.Reallocate(quadrangle_reference.GetNbNodalBoundary(num_loc));
          for (int j = 0; j < num.GetM(); j++)
            num(j) = quadrangle_reference.GetNodalNumber(num_loc, j);
        }
        break;
      }
    
    return num;
  }

  
  //! sets the order of approximation of the mesh
  template<class T>
  void Mesh<Dimension2, T>::SetGeometryOrder(int r, bool only_quad)
  {
    // if the order is the same, nothing to do
    if (this->lob_curve.GetOrder() == r)
      return;
    
    Globatto<Real_wp> old_lob;
    old_lob = this->lob_curve;
    
    // nodal points on Gauss-Lobatto points
    this->lob_curve.ConstructQuadrature(r, this->lob_curve.QUADRATURE_LOBATTO);
    this->lob_curve.ComputeGradPhi(1e3*epsilon_machine);
    VectReal_wp xi;
    if (r > 1)  
      {
	xi.Reallocate(r-1);
	for (int i = 1; i < r; i++)
	  xi(i-1) = this->lob_curve.Points(i);
      }
    
    this->interval_reference.ConstructFiniteElement(r);
    
    // we try to recompute PointsEdgeRef with the new order of approximation provided
    if (this->GetNbBoundaryRef() > 0)
      {
	if (r <= 1)
	  this->PointsEdgeRef.Clear();
	else
	  {
	    if (old_lob.GetOrder() <= 1)
	      {
		this->PointsEdgeRef.Reallocate(this->GetNbBoundaryRef(), r-1);
		// basic linear interpolation		
		for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		  {
		    int ref = this->BoundaryRef(i).GetReference();		    
		    if ((this->GetCurveType(ref) > 0)
                        &&(this->GetCurveType(ref) != this->CURVE_FILE))
		      {
			Vector<R_N> points_interm;
			this->GetPointsOnCurve(i, ref, xi, points_interm);
			for (int k = 0; k < xi.GetM(); k++)
			  this->PointsEdgeRef(i, k) = points_interm(k);
		      }
		    else
		      {
			R_N ptA = this->Vertex(this->BoundaryRef(i).numVertex(0));
			R_N ptB = this->Vertex(this->BoundaryRef(i).numVertex(1));
			for (int k = 0; k < xi.GetM(); k++)
			  {
			    this->PointsEdgeRef(i, k).Zero();
			    Add(1.0-xi(k), ptA, this->PointsEdgeRef(i, k));
			    Add(xi(k), ptB, this->PointsEdgeRef(i, k));
			  }
		      }
		  }
	      }
	    else
	      {
		int m = old_lob.GetOrder();
		Matrix<Real_wp> ValuePhi(m+1, r+1);
		for (int j = 0; j <= m; j++)
		  for (int k = 1; k < r; k++)
		    ValuePhi(j, k-1) = old_lob.EvaluatePhi(j, xi(k-1));
		
		Matrix<R_N> OldPointsEdge = this->PointsEdgeRef;
		this->PointsEdgeRef.Reallocate(this->GetNbBoundaryRef(), r-1);
		// interpolation with previous points
		for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		  {
		    int ref = this->BoundaryRef(i).GetReference();		    
		    if ((this->GetCurveType(ref) > 0)
                        &&(this->GetCurveType(ref) != this->CURVE_FILE))
		      {
			Vector<R_N> points_interm;
			this->GetPointsOnCurve(i, ref, xi, points_interm);
			for (int k = 0; k < xi.GetM(); k++)
			  this->PointsEdgeRef(i, k) = points_interm(k);
		      }
		    else
		      {
			VectR_N OldPoints(m+1);
			OldPoints(0) = this->Vertex(this->BoundaryRef(i).numVertex(0));
			OldPoints(m) = this->Vertex(this->BoundaryRef(i).numVertex(1));
			for (int k = 0; k < m-1; k++)
			  OldPoints(k+1) = OldPointsEdge(i, k);
						
			for (int k = 0; k < xi.GetM(); k++)
			  {
			    this->PointsEdgeRef(i, k).Zero();
			    for (int j = 0; j <= m; j++)
			      Add(ValuePhi(j, k), OldPoints(j), this->PointsEdgeRef(i, k));
			  }
		      }
		  }
	      }
	  }
      }
    
    // construction of reference elements
    MeshNumbering<Dimension2> mesh_num(*this);
    if (!only_quad)
      {
        triangle_reference.ConstructFiniteElement(r);
        triangle_reference.SetMesh(*this);
      }
    
    quadrangle_reference.ConstructFiniteElement(r);
    quadrangle_reference.SetMesh(*this);
  }
  

  /************************
   * Convenient Functions *
   ************************/

  
  //! modifies the number of referenced edges, previous edges are not kept
  template<class T>
  void Mesh<Dimension2, T>::ReallocateBoundariesRef(int i)
  {
    this->edges_ref.Reallocate(i);
    this->PointsEdgeRef.Reallocate(i, this->GetGeometryOrder()-1);
    this->value_parameter_boundary.Reallocate(i);
    this->integer_parameter_boundary.Reallocate(i);
  }


  //! modifies the number of referenced edges, previous edges are kept
  template<class T>
  void Mesh<Dimension2, T>::ResizeBoundariesRef(int i)
  {
    this->edges_ref.Resize(i);
    this->PointsEdgeRef.Resize(i, this->GetGeometryOrder()-1);

    this->value_parameter_boundary.Resize(i);
    this->integer_parameter_boundary.Resize(i);
  }


  //! fills an array containing reference elements for elements
  template<class T>
  void Mesh<Dimension2, T>
  ::GetReferenceElementVolume(Vector<const ElementGeomReference<Dimension2>* >& ref) const
  {
    int nb_tri = 0, nb_quad = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	if (this->Element(i).GetNbVertices() == 3)
	  nb_tri++;
	else
	  nb_quad++;
      }

    if ((nb_tri > 0) && (nb_quad > 0))
      ref.Reallocate(2);
    else
      ref.Reallocate(1);
    
    int nb = 0;
    if (nb_tri > 0)
      ref(nb++) = &this->triangle_reference;

    if (nb_quad > 0)
      ref(nb++) = &this->quadrangle_reference;
  }
  

  //! returns the number of triangles contained in the mesh
  template<class T>
  int Mesh<Dimension2, T>::GetNbTriangles() const
  {
    int nb_tri = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->elements(i).GetNbVertices() == 3)
	nb_tri++;
    
    return nb_tri;
  }
  

  //! returns the number of quadrilaterals contained in the mesh
  template<class T>
  int Mesh<Dimension2, T>::GetNbQuadrangles() const
  {
    int nb_quad = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->elements(i).GetNbVertices() == 4)
	nb_quad++;
    
    return nb_quad;
  }

  
    
  template<class T>
  size_t Mesh<Dimension2, T>::GetMemorySize() const
  {
    size_t taille = Mesh_Base<Dimension2, T>::GetMemorySize()
      + MeshBoundaries<Dimension2, T>::GetMemorySize() + sizeof(*this);
    
    taille += triangle_reference.GetMemorySize() + quadrangle_reference.GetMemorySize();
    return taille;
  }
  
  
  //! returns a bounding number of edges, if number of elements is known
  template<class T>
  int Mesh<Dimension2, T>::GetUpperBoundNumberOfEdges() const
  {
    int nb_quad = this->GetNbQuadrangles();
    int nb_tri = this->GetNbTriangles();
    
    return (4*nb_quad + 3*nb_tri);
  }
  
  
  //! connectivity (edges, edge numbers) is cleared
  template<class T>
  void Mesh<Dimension2, T>::ClearConnectivity()
  {
    Mesh_Base<Dimension2, T>::ClearConnectivity();
    MeshBoundaries<Dimension2, T>::ClearConnectivity();
  }
  
    
  //! all the mesh is removed
  template<class T>
  void Mesh<Dimension2, T>::Clear()
  {
    Mesh_Base<Dimension2>::Clear();
    MeshBoundaries<Dimension2, T>::Clear();

    for (int i = 0; i < this->GetNbReferences(); i++)
      if (this->GetCurveType(i) == this->CURVE_FILE)
        this->SetCurveType(i, this->NO_CURVE);
  }
  
    
  /****************
   * Input/Output *
   ****************/
  
  
  //! modification of parameters with a line of the data file
  template<class T>
  void Mesh<Dimension2, T>
  ::SetInputData(const string& description_field, const Vector<string>& parameters)
  {
    Mesh_Base<Dimension2, T>::SetInputData(description_field, parameters);
    MeshBoundaries<Dimension2, T>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("TypeMesh"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh" << endl;
	    cout << "TypeMesh needs more parameters, for instance :" << endl;
	    cout << "TypeMesh = QUAD" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	this->required_mesh_type = Mesh<Dimension2>::TRIANGULAR_MESH;
	if (!parameters(0).compare("TRI"))
	  this->required_mesh_type = Mesh<Dimension2>::TRIANGULAR_MESH;
	else if (!parameters(0).compare("QUAD"))
	  this->required_mesh_type = Mesh<Dimension2>::QUADRILATERAL_MESH;
	else if (!parameters(0).compare("RADAU"))
	  this->required_mesh_type = Mesh<Dimension2>::RADAU_MESH;
      }

  }
  
  
  //! constructs arrays useful to write high-order curved mesh
  /*!
    \param[out] PosNodes all the nodes of the mesh (even internal nodes due to high-order)
    \param[out] Nodle Numbering of all the nodes (as for continuous finite elements)
    \param[in] rmax maximal allowed order (if you want only nodes for first order mesh, put 1)
  */
  template<class T>
  int Mesh<Dimension2, T>
  ::GetNodesCurvedMesh(VectR_N& PosNodes, Vector<IVect>& Nodle,
                       int rmax, bool lobatto) const
  {
    int r = this->lob_curve.GetOrder();
    int order = min(rmax, r);
    if (order == 1)
      {
	PosNodes = this->Vertices;
	Nodle.Reallocate(this->GetNbElt());
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    int nv = this->Element(i).GetNbVertices();
	    Nodle(i).Reallocate(nv);
	    for (int j = 0; j < nv; j++)
	      Nodle(i)(j) = this->Element(i).numVertex(j);
	  }
      }
    else
      {
	Vector<VectR_N> points_div(2);
	points_div(0).Reallocate((order+1)*(order+2)/2);
	points_div(1).Reallocate((order+1)*(order+1));
	Matrix<int> NumNodes2D, coor;
        
	VectReal_wp step_x(order+1), weights_x;
	if (lobatto)
          {
            if (this->GetNbTriangles() > 0)
	      TriangleGeomReference::
		ConstructLobattoPoints(order, 0, step_x, points_div(0));
            else
              ComputeGaussLobatto(step_x, weights_x, order);
          }
        else
          {
            for (int i = 0; i <= order; i++)
              step_x(i) = Real_wp(i)/order;
            
            // regular points on triangles
            MeshNumbering<Dimension2>::ConstructTriangularNumbering(order, NumNodes2D, coor);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order-i; j++)
                points_div(0)(NumNodes2D(i, j)).Init(step_x(i), step_x(j));            
          }
	
	// tensorized points on quadrangles
	MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(order, NumNodes2D, coor);
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    points_div(1)(NumNodes2D(i, j)).Init(step_x(i), step_x(j));
        
	Vector<IVect> NodleSurf;
	this->GetMeshSubdivision(step_x, points_div, PosNodes, Nodle, NodleSurf);
      }
    
    return order;
  }
  
  
  //! reading of a mesh, several formats are available
  /*!
    \param[in] file_mesh name of the file where the mesh is stored
  */
  template<class T>
  void Mesh<Dimension2, T>::Read(const string & file_mesh) 
  { 
    // previous mesh is removed
    Clear();
    // extension of the mesh
    string extension_input = GetExtension(file_mesh);
    
    // fichier de maillage en binaire
    ifstream file_in(file_mesh.data());
    if (!file_in.is_open())
      {
	cout<<"Mesh file not found"<<endl;
	cout<<file_mesh<<endl;
	abort();
      }
    
    Read(file_in, extension_input);
    
    file_in.close();
  }
  
  
  //! reading of a mesh, several formats are available
  /*!
    \param[in] file_in stream from which the mesh is read
    \param[in] extension_input extension of the file mesh
    The extension is used to know the format of the mesh data
  */  
  template<class T>
  void Mesh<Dimension2, T>::Read(istream& file_in, const string& extension_input) 
  {    
    if (!extension_input.compare("t"))
      {
	int dimension, itmp, ref, nb_vertices_loc;
	int nb_elt, nb_vertices;
	file_in>>dimension>>nb_elt>>nb_vertices>>itmp;
	if (dimension != 2)
	  {
	    cout<<"Ce n'est pas un maillage de dimension 2"<<endl;
	    abort();
	  }
	
	// lecture des sommets
	this->ReallocateVertices(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  file_in>>this->Vertex(i)(0)>>this->Vertex(i)(1);
	
	int nb_edges_ref = 0;
	// lecture des quadrangles/triangles
	this->ReallocateElements(nb_elt); IVect num(4); num.Fill(-1);
	// on surdimensionne aretes de reference
	Vector<TinyVector<int,2> > num_edges(4*nb_elt);
	
	IVect ref_edges(4*nb_elt);
	for (int i = 0; i < nb_elt; i++)
	  {
	    file_in>>nb_vertices_loc>>ref>>itmp;
	    // sommets de l'element
	    for (int j = 0; j < nb_vertices_loc; j++)
	      file_in>>num(j);
	    
	    if (nb_vertices_loc == 3)
	      this->Element(i).InitTriangular(num(0)-1, num(1)-1, num(2)-1, ref);
	    else
	      this->Element(i).InitQuadrangular(num(0)-1, num(1)-1, num(2)-1, num(3)-1, ref);
	    
	    // autres ddl, on passe
	    for (int j = nb_vertices_loc; j < itmp; j++)
	      file_in>>itmp;
	    
	    // references des aretes
	    for (int j = 0; j < nb_vertices_loc; j++)
	      {
		file_in>>ref;
		if (ref != 0)
		  {
		    // arete de reference
		    int n0 = num(j)-1; int n1 = num((j+1)%nb_vertices_loc)-1;
		    Sort(n0, n1);
		    num_edges(nb_edges_ref)(0) = n0;
		    num_edges(nb_edges_ref)(1) = n1;
		    ref_edges(nb_edges_ref++) = ref;
		  }
	      }
	  }
	
	// on elimine doublons
	RemoveDuplicate(nb_edges_ref, num_edges, ref_edges);
	
	this->ReallocateBoundariesRef(nb_edges_ref);
	for (int i = 0; i < nb_edges_ref; i++)
	  this->BoundaryRef(i).Init(num_edges(i)(0), num_edges(i)(1), ref_edges(i));
	
      }
    else if (!extension_input.compare("vtk"))
      {
	string ligne;
	getline(file_in, ligne);
	getline(file_in, ligne);
	getline(file_in, ligne);
	bool ascii = false;
	if (ligne == "ASCII")
	  ascii = true;
	
	getline(file_in, ligne);
	if (ligne.find("UNSTRUCTURED_GRID") == string::npos)
	  {
	    cout << "Only unstructured grid accepted " << endl;
	    abort();
	  }
	
	string keyword, param, type_data;
	file_in >> keyword >> param >> type_data;
	if (keyword != "POINTS")
	  {
	    cout << "POINTS expected" << endl;
	    abort();
	  }
	
	bool double_prec = false;
	if (type_data == "double")
	  double_prec = true;
	else if (type_data == "float")
	  double_prec = false;
	else
	  {
	    cout << "Unexpected type of data : " << type_data << endl;
	    abort();
	  }
	
	int nb_vert = to_num<int>(param);
	this->ReallocateVertices(nb_vert);
	Real_wp z = 0;
	if (ascii)
	  for (int i = 0; i < this->GetNbVertices(); i++)
	    file_in >> this->Vertex(i)(0) >> this->Vertex(i)(1) >> z;
	else
	  {
	    getline(file_in, ligne);
	    VectReal_wp Pts(3*this->GetNbVertices());
	    ReadBinaryDoubleOrFloat(Pts, file_in, double_prec, false, true);
	    for (int i = 0; i < this->GetNbVertices(); i++)
	      this->Vertex(i).Init(Pts(3*i), Pts(3*i+1));	    
	  }
	
	// reading CELLS
	int nb_entity = 0, nb_all_points = 0;
	file_in >> keyword >> nb_entity >> nb_all_points;
	if (keyword != "CELLS")
	  {
	    cout << "Expected CELLS but got " << keyword << endl;
	    abort();
	  }
	
	IVect AllNum(nb_all_points);
	if (ascii)
	  for (int i = 0; i < nb_all_points; i++)
	    file_in >> AllNum(i);
	else
	  {
	    getline(file_in, ligne);
	    ReadBinaryInteger(AllNum, file_in, false, true);
	  }
	
	// reading CELL_TYPES
	file_in >> keyword >> nb_entity;
	if (keyword != "CELL_TYPES")
	  {
	    cout << "Expected CELL_TYPES but got " << keyword << endl;
	    abort();
	  }

	IVect TypeEntity(nb_entity), RefEntity(nb_entity);
	RefEntity.Fill(1);
	if (ascii)
	  for (int i = 0; i < nb_entity; i++)
	    file_in >> TypeEntity(i);
	else
	  {
	    getline(file_in, ligne);
	    ReadBinaryInteger(TypeEntity, file_in, false, true);
	  }

	file_in >> keyword >> nb_entity;
	if (file_in.good() && (keyword == "CELL_DATA"))
	  {
	    // references
	    // lines SCALARS cell_scalars int 1 and LOOKUP_TABLE default
	    getline(file_in, ligne);
	    getline(file_in, ligne);
	    getline(file_in, ligne);
	    if (ascii)
	      for (int i = 0; i < nb_entity; i++)
		file_in >> RefEntity(i);
	    else
	      ReadBinaryInteger(RefEntity, file_in, false, true);	
	  }
	
	// counting edges, triangles and quadrangles
	int nb_edge_ref = 0, nb_tri = 0, nb_quad = 0;	
	for (int i = 0; i < nb_entity; i++)
	  {
	    if (TypeEntity(i) == 3)
	      nb_edge_ref++;
	    else if (TypeEntity(i) == 5)
	      nb_tri++;
	    else if (TypeEntity(i) == 9)
	      nb_quad++;
	  }
	
	this->ReallocateBoundariesRef(nb_edge_ref);
	this->ReallocateElements(nb_tri + nb_quad);
	int nb = 0;
	nb_edge_ref = 0; int nb_elt = 0;
	for (int i = 0; i < nb_entity; i++)
	  {
	    // int nb_vert = AllNum(nb);
	    if (TypeEntity(i) == 3)
	      {
		this->BoundaryRef(nb_edge_ref).Init(AllNum(nb+1), AllNum(nb+2), RefEntity(i));
		nb += 3; nb_edge_ref++;
	      }
	    else if (TypeEntity(i) == 5)
	      {
		this->Element(nb_elt).InitTriangular(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), RefEntity(i));
		nb += 4; nb_elt++;
	      }
	    else if (TypeEntity(i) == 9)
	      {
		this->Element(nb_elt)
                  .InitQuadrangular(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), AllNum(nb+4), RefEntity(i));
		nb += 5; nb_elt++;
	      }
	  }
      }
    else if (!extension_input.compare("ele"))
      {
	abort();
	/*
	// on lit d'abord les elements
	int itmp, n1, n2, n3, ref, nb_elt;
	file_in>>nb_elt>>itmp>>itmp;
	this->ReallocateElements(nb_elt);
	for (int i = 0; i < nb_elt; i++)
	{
	file_in>>itmp>>n1>>n2>>n3>>ref;
	this->Element(i).InitTriangular(n1-1, n2-1, n3-1, ref);
	}
	
	// puis les points
	string name_file = GetBaseString(file_mesh) + string(".node");
	file_in.open(name_file.data());
	if (!file_in.is_open())
	{
	cout<<"Mesh file not found"<<endl;
	cout<<file_mesh<<endl;
	abort();
	}
        
	int nb_vertices;
	file_in>>nb_vertices>>itmp>>itmp>>itmp;
	this->ReallocateVertices(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	{
	file_in>>itmp>>this->Vertex(i)(0)>>this->Vertex(i)(1)>>ref;
	}
	
	file_in.close();
	
	// puis les aretes
	name_file = GetBaseString(file_mesh) + string(".edge");
	file_in.open(name_file.data());
	if (!file_in.is_open())
	{
	cout<<"Mesh file not found"<<endl;
	cout<<file_mesh<<endl;
	abort();
	}
	int nb_edge_all; file_in>>nb_edge_all>>itmp;
	this->ReallocateBoundariesRef(nb_edge_all);
	int nb_edges_ref = 0;
	for (int i = 0; i < nb_edge_all; i++)
	{
	file_in>>itmp>>n1>>n2>>ref;
	if (ref != 0)
	{
	this->BoundaryRef(nb_edges_ref).Init(n1-1, n2-1, ref);
	nb_edges_ref++;
	}
	}
	this->ResizeBoundariesRef(nb_edges_ref);
	
	file_in.close();
	*/
      }
    else if (!extension_input.compare("geod"))
      {
        string ligne;
        int num, nb_elt, nb_pts;
        Real_wp x, y, z;
        int n1, n2, n3, ref;
        
        // en-tete region 1
        getline(file_in, ligne);
        getline(file_in, ligne);
        getline(file_in, ligne);
        
        istringstream data_ligne(ligne);
        data_ligne >> nb_elt >> nb_pts;
        
        this->ReallocateVertices(nb_pts);
        this->ReallocateElements(nb_elt);
        
        // lecture sommets region 1
        for (int k = 0; k < nb_pts; k++)
          {      
            getline(file_in, ligne);
            for (unsigned i = 0; i < ligne.size(); i++)
              if (ligne[i] == 'D')
                ligne[i] = 'E';
            
            istringstream data_ligne_reel(ligne);
            data_ligne_reel >> num >> x >> y >> z;
            this->Vertex(k).Init(x, y);
          }
        
        // lecture elements region 1
        for (int k = 0; k < nb_elt; k++)
          {
            getline(file_in, ligne);
            istringstream data_ligne_entier(ligne);
            data_ligne_entier >> num >> n1 >> n2 >> n3 >> ref;
            
            n1--; n2--; n3--;
            this->Element(k).InitTriangular(n1, n2, n3, 1);
          }
        
        // en-tete region 2
        getline(file_in, ligne);
        getline(file_in, ligne);
        getline(file_in, ligne);
        
        istringstream data_ligne2(ligne);
        data_ligne2 >> nb_elt >> nb_pts;
        
        Mesh<Dimension2> mesh2;
        mesh2.ReallocateVertices(nb_pts);
        mesh2.ReallocateElements(nb_elt);
        
        // lecture sommets region 2
        for (int k = 0; k < nb_pts; k++)
          {      
            getline(file_in, ligne);
            for (unsigned i = 0; i < ligne.size(); i++)
              if (ligne[i] == 'D')
                ligne[i] = 'E';
            
            istringstream data_ligne_reel(ligne);
            data_ligne_reel >> num >> x >> y >> z;
            mesh2.Vertex(k).Init(x, y);
          }
        
        // lecture elements region 2
        for (int k = 0; k < nb_elt; k++)
          {
            getline(file_in, ligne);
            istringstream data_ligne_entier(ligne);
            data_ligne_entier >> num >> n1 >> n2 >> n3 >> ref;
            
            n1--; n2--; n3--;      
            mesh2.Element(k).InitTriangular(n1, n2, n3, 2);
          }
        
        
        // en-tete region 3
        getline(file_in, ligne);
        getline(file_in, ligne);
        getline(file_in, ligne);
        
        istringstream data_ligne3(ligne);
        data_ligne3 >> nb_elt >> nb_pts;
        
        Mesh<Dimension2> mesh3;
        mesh3.ReallocateVertices(nb_pts);
        mesh3.ReallocateElements(nb_elt);
        
        // lecture sommets region 3
        for (int k = 0; k < nb_pts; k++)
          {      
            getline(file_in, ligne);
            for (unsigned i = 0; i < ligne.size(); i++)
              if (ligne[i] == 'D')
                ligne[i] = 'E';
            
            istringstream data_ligne_reel(ligne);
            data_ligne_reel >> num >> x >> y >> z;
            mesh3.Vertex(k).Init(x, y);
          }
        
        // lecture elements region 3
        for (int k = 0; k < nb_elt; k++)
          {
            getline(file_in, ligne);
            istringstream data_ligne_entier(ligne);
            data_ligne_entier >> num >> n1 >> n2 >> n3 >> ref;
            
            n1--; n2--; n3--;
            mesh3.Element(k).InitTriangular(n1, n2, n3, 2);
          }
        
        this->ReorientElements();
        this->FindConnectivity();
        this->AddBoundaryEdges();
        this->RedistributeReferences();
        this->ProjectPointsOnCurves();
        
        this->AppendMesh(mesh2);
        this->AppendMesh(mesh3);
        
        this->ReorientElements();            
        this->FindConnectivity();
        this->AddBoundaryEdges();
        
        this->RedistributeReferences();
        this->Param_known.Fill(false);
        this->ProjectPointsOnCurves();
        this->Write("original.mesh");
      }
    else if (!extension_input.compare("exo"))
      {
	// ascii exodus mesh
	
	// cout<<"lecture fichier exo"<<endl;
	string ligne; Vector<string> parameters;
	int nb_elt = 0; int nb_vertices = 0; int nb_edges_ref = 0;
	
	// reading dimensions first
	string dim_chaine = "dimensions:";
	string variables_chaine = "variables:";
	string data_chaine = "data:";
	bool test_loop = true;
	while (test_loop)
	  {
	    getline(file_in, ligne);
	    if (!ligne.compare(0, dim_chaine.size(), dim_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	test_loop = true;
	IVect nb_nodes_per_block(10); nb_nodes_per_block.Fill(4);
	IVect nb_faces_block(10); nb_faces_block.Fill(0);
	IVect nb_elt_block(10); nb_elt_block.Fill(0);
	while (test_loop)
	  {
	    getline(file_in, ligne);
	    if (!ligne.compare(0,variables_chaine.size(), variables_chaine))
	      test_loop = false;
	    else
	      {
		// tokenization of the line
		if ((ligne.size() > 0) && (ligne[0] != '#') && (ligne[0] != '/'))
		  {
		    StringTokenize(ligne, parameters, string(" \t"));
		    if (!parameters(0).compare("num_nodes"))
		      nb_vertices = to_num<int>(parameters(2));
		    else if (!parameters(0).compare("num_elem"))
		      nb_elt = to_num<int>(parameters(2));
		    else if (!parameters(0).compare(0,14,"num_nod_per_el"))
		      {
			int ref = to_num<int>(parameters(0).substr(14));
			nb_nodes_per_block(ref) = to_num<int>(parameters(2));
		      }
		    else if (!parameters(0).compare("num_el_blk"))
		      {
			int nb_ref = to_num<int>(parameters(2));
			nb_elt_block.Reallocate(nb_ref+1); nb_elt_block.Fill(0);
			nb_nodes_per_block.Reallocate(nb_ref+1);
			nb_nodes_per_block.Fill(4);
		      }
		    else if (!parameters(0).compare("num_side_sets"))
		      {
			int nb_ref = to_num<int>(parameters(2));
			nb_faces_block.Reallocate(nb_ref+1);
			nb_faces_block.Fill(0);
		      }
		    else if (!parameters(0).compare(0,11,"num_side_ss"))
		      {
			int ref = to_num<int>(parameters(0).substr(11));
			nb_faces_block(ref) = to_num<int>(parameters(2));
			nb_edges_ref += nb_faces_block(ref);
		      }
		    else if (!parameters(0).compare(0,13,"num_el_in_blk"))
		      {
			int ref = to_num<int>(parameters(0).substr(13));
			nb_elt_block(ref) = to_num<int>(parameters(2));
		      }
		  }
	      }
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	this->ReallocateVertices(nb_vertices);
	this->ReallocateElements(nb_elt);
	this->ReallocateBoundariesRef(nb_edges_ref);
	
	test_loop = true;
	// we skip variables: field
	while (test_loop)
	  {
	    getline(file_in, ligne);
	    if (!ligne.compare(0,data_chaine.size(), data_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	
	IVect ElemFaceRef(nb_edges_ref); IVect LocalFaceRef(nb_edges_ref);
	IVect ReferenceFace(nb_edges_ref); IVect elem_map(nb_elt); char comma;
	nb_edges_ref = 0; int num_local = 0; nb_elt = 0;
	test_loop = true;
	
	// now treating datas
	while (test_loop)
	  {
	    getline(file_in, ligne);
	    if ((ligne.size() != 0) && (ligne[0] != '#') && (ligne[0]!='/'))
	      {
		StringTokenize(ligne, parameters, string(" ,\t;"));
		if ((parameters.GetM() >= 2)&&(!parameters(1).compare("=")))
		  {
		    string keyword = parameters(0);

		    for (int i = 2; i < parameters.GetM(); i++)
		      parameters(i-2) = parameters(i);
		    
		    parameters.Resize(parameters.GetM()-2);
		    
		    if (!keyword.compare("elem_map"))
		      {
			// reading map number
			for (int i = 0; i < parameters.GetM(); i++)
			  elem_map(i) = to_num<int>(parameters(i))-1;
			
			for (int i = parameters.GetM(); i < elem_map.GetM(); i++)
			  {
			    file_in>>elem_map(i)>>comma;
			    elem_map(i)--;
			  }
		      }
		    else if (!keyword.compare(0,7,"elem_ss"))
		      {
			// reading referenced faces of specified reference ref
			// here element numbers that belong to faces
			int ref = to_num<int>(keyword.substr(7));
			for (int i = 0; i < parameters.GetM(); i++)
			  {
			    ElemFaceRef(nb_edges_ref) = to_num<int>(parameters(i))-1;
			    ReferenceFace(nb_edges_ref) = ref;
			    nb_edges_ref++;
			  }
			
			for (int i = parameters.GetM(); i < nb_faces_block(ref); i++)
			  {
			    file_in>>ElemFaceRef(nb_edges_ref)>>comma;
			    ElemFaceRef(nb_edges_ref)--;
			    ReferenceFace(nb_edges_ref) = ref;
			    nb_edges_ref++;
			  }
		      }
		    else if (!keyword.compare(0,12,"dist_fact_ss"))
		      {
			int ref = to_num<int>(keyword.substr(12)); int itmp;
			for (int i = parameters.GetM(); i < nb_faces_block(ref); i++)
			  file_in>>itmp>>comma;
		      }
		    else if (!keyword.compare(0,7,"side_ss"))
		      {
			// reading local number in the element of the referenced faces
			int ref = to_num<int>(keyword.substr(7));
			for (int i = 0; i < parameters.GetM(); i++)
			  {
			    LocalFaceRef(num_local) = to_num<int>(parameters(i))-1;
			    num_local++;
			  }
			
			for (int i = parameters.GetM(); i < nb_faces_block(ref); i++)
			  {
			    file_in>>LocalFaceRef(num_local)>>comma;
			    LocalFaceRef(num_local)--;
			    num_local++;
			  }
			
		      }
		    else if (!keyword.compare(0,7,"connect"))
		      {
			// local to global numbering
			int ref = to_num<int>(keyword.substr(7));
			int nb_vert = nb_nodes_per_block(ref);
			int N = nb_vert*nb_elt_block(ref); 
			IVect all_number(N);
			for (int i = 0; i < parameters.GetM(); i++)
			  all_number(i) = to_num<int>(parameters(i))-1;
			 			
			for (int i = parameters.GetM(); i < N; i++)
			  {
			    file_in>>all_number(i)>>comma; all_number(i)--;
			  }
			
			IVect num(nb_vert); int nb = 0;
			for (int i = 0; i < all_number.GetM()/nb_vert; i++)
			  {
			    for (int j = 0; j < nb_vert; j++)
			      num(j) = all_number(nb++);
			    
			    
			    if (nb_vert == 9)
			      {
				cout<<"case not implemented "<<endl;
				abort();				
			      }
			    else if (nb_vert == 4)
			      {
				this->elements(nb_elt).Init(num, ref);
			      }
			    else if (nb_vert == 3)
			      {
				this->elements(nb_elt).Init(num, ref);
			      }
			    
			    nb_elt++;
			  }
		      }
		    else if (!keyword.compare("coord"))
		      {
			// coordinates of all the vertices
			for (int i = 0; i < parameters.GetM(); i++)
			  this->Vertex(i)(0) = to_num<T>(parameters(i));
			
			for (int i = parameters.GetM(); i < nb_vertices; i++)
			  file_in>>this->Vertex(i)(0)>>comma;
			
			for (int i = 0; i < nb_vertices; i++)
			  file_in>>this->Vertex(i)(1)>>comma;
			
		      }
		  }
		
	      }
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	// now reconstructing faces_ref
	// local numbers of montjoie and exodus are different	
	for (int i = 0; i < nb_edges_ref; i++)
	  {
	    int num_elem = ElemFaceRef(i);
	    int type_loc = this->Element(num_elem).GetHybridType();
	    int num_loc = LocalFaceRef(i);
	    int ref = ReferenceFace(i);
	    
	    int n0 = this->Element(num_elem).numVertex
	      (MeshNumbering<Dimension2>::FirstExtremityEdge(type_loc, num_loc));
	    int n1 = this->Element(num_elem).numVertex
	      (MeshNumbering<Dimension2>::SecondExtremityEdge(type_loc, num_loc));
	    
	    this->BoundaryRef(i).Init(n0, n1, ref);
	  }	
      }
    else if (!extension_input.compare("msh"))
      {
	// Gmsh file
	string ligne; int itmp;
	file_in>>ligne;
	StringTrim(ligne); int version = 1;
	if (!ligne.compare("$MeshFormat"))
	  {
	    version = 2; double xtmp;
	    file_in>>xtmp>>itmp>>itmp;
	    file_in>>ligne;
	  }
	
        while ((ligne != "$Nodes") && (ligne != "$NOD"))
          {
            file_in >> ligne;
            if (file_in.eof())
              {
                cout << "Problem while reading mesh" << endl;
                abort();
              }
          }
        
	int nb_all_nodes;
	file_in >> nb_all_nodes;
	
	IVect node_numbers(nb_all_nodes);
	VectR_N AllNodes(nb_all_nodes); T z;
	// we read all the nodes
	for (int i = 0; i < nb_all_nodes; i++)
	  file_in >> node_numbers(i) >> AllNodes(i)(0) >> AllNodes(i)(1) >> z;
	
	int numero_max = 0;
	for (int i = 0; i < nb_all_nodes; i++)
	  numero_max = max(node_numbers(i), numero_max);
        
        // line $EndNodes
        file_in>>ligne;
        
        while ((ligne != "$Elements") && (ligne != "$ELM"))
          {
            file_in >> ligne;
            if (file_in.eof())
              {
                cout << "Problem while reading mesh" << endl;
                abort();
              }
          }
        
	// we read all the elements and edges
	VectBool UsedVertices(numero_max + 1); UsedVertices.Fill(false);
	int nb_entites;
	file_in>>nb_entites;
	IVect type_entite(nb_entites), ref_entite(nb_entites), geom_entite(nb_entites);
	Vector<IVect> num_nodes(nb_entites); IVect nb_vertices_entite(nb_entites);
	int nb_elements = 0, nb_boundaries = 0; bool curved_mesh = false; int order_mesh = -1;
	for (int i = 0; i < nb_entites; i++)
	  {
	    int nb_nodes(0), nb_tags, nb_vert(0);
	    if (version == 2)
	      {
		file_in>>itmp>>type_entite(i)>>nb_tags>>ref_entite(i)>>geom_entite(i);
		for (int j = 2; j < nb_tags; j++)
		  file_in>>itmp;		
		
		nb_nodes = NumberMap::NbNodesGmshEntity(type_entite(i));
	      }
	    else
	      file_in>>itmp>>type_entite(i)>>ref_entite(i)>>geom_entite(i)>>nb_nodes;
	    
	    nb_vert = NumberMap::NbVerticesGmshEntity(type_entite(i));
	    nb_vertices_entite(i) = nb_vert;
	    if (nb_vert == 2)
	      nb_boundaries++;
	    else
	      nb_elements++;
	    
	    if (nb_nodes > 4)
	      {
		curved_mesh = true;
		int r = NumberMap::OrderGmshEntity(type_entite(i));
		if (order_mesh == -1)
		  order_mesh = r;
		else
		  if (order_mesh != r)
		    {
		      cout<<"case not treated"<<endl;
                      cout<<"Did you write Mesh.SecondOrderIncomplete = 0; ? " << endl;
		      abort();
		    }
		
	      }
	    
	    num_nodes(i).Reallocate(nb_nodes);
	    for (int j = 0; j < nb_nodes; j++)
	      file_in>>num_nodes(i)(j);
	    
	    // specifying used vertices
	    for (int j = 0; j < nb_vert; j++)
	      UsedVertices(num_nodes(i)(j)) = true;
	    
	  }
	
	// we retrieve actual number of vertices
	int nb_real_vertices = 0;
	IVect perm(numero_max + 1); perm.Fill(-1);
	for (int i = 0; i < nb_all_nodes; i++)
	  if (UsedVertices(node_numbers(i)))
	    perm(node_numbers(i)) = nb_real_vertices++;
	
	this->ReallocateVertices(nb_real_vertices);
	nb_real_vertices = 0;
	for (int i = 0; i < nb_all_nodes; i++)
	  if (UsedVertices(node_numbers(i)))
	    this->Vertex(nb_real_vertices++) = AllNodes(i);
		
	// on construit les elements H1
	this->ReallocateElements(nb_elements);
	this->ReallocateBoundariesRef(nb_boundaries);
	IVect NumEntiteBoundary(nb_boundaries), NumEntiteElement(nb_elements);
	nb_elements = 0; nb_boundaries = 0; IVect num(4);
	for (int i = 0; i < nb_entites; i++)
	  {
	    num.Fill(-1);
	    // permutation
	    for (int j = 0; j < nb_vertices_entite(i); j++)
	      num(j) = perm(num_nodes(i)(j));
	    
	    if (nb_vertices_entite(i) == 2)
	      {
		this->BoundaryRef(nb_boundaries).Init(num(0), num(1), ref_entite(i));
		this->BoundaryRef(nb_boundaries).SetGeometryReference(geom_entite(i));
		NumEntiteBoundary(nb_boundaries) = i;
		nb_boundaries++;
	      }
	    else if (nb_vertices_entite(i) == 3)
	      {
		this->Element(nb_elements).InitTriangular(num(0), num(1), num(2), ref_entite(i));
		NumEntiteElement(nb_elements) = i;
		nb_elements++;
	      }
	    else if (nb_vertices_entite(i) == 4)
	      {
		this->Element(nb_elements).InitQuadrangular(num(0), num(1), num(2), num(3),
							    ref_entite(i));
		NumEntiteElement(nb_elements) = i;
		nb_elements++;
	      }
	  }	
        
	if (this->lob_curve.GetOrder() > 1)
	  if (curved_mesh)
	    {
	      // construction of reciprocal array for node_numbers
	      for (int i = 0; i < AllNodes.GetM(); i++)
		perm(node_numbers(i)) = i;
	      
	      // we search an available reference
	      int ref_max = 0;
	      for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		ref_max = max(ref_max, this->BoundaryRef(i).GetReference());
	      
	      Vector<bool> RefUsed(ref_max+2); RefUsed.Fill(false);
	      for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		RefUsed(this->BoundaryRef(i).GetReference()) = true;
	      
	      int new_ref = this->GetNewReference();
              
	      RefUsed.Fill(false);
	      int r = this->lob_curve.GetOrder();
	      // computation of ValPhi for interpolation
	      VectReal_wp RegPoints(order_mesh+1);
	      RegPoints.Fill(); Mlt(Real_wp(1)/order_mesh, RegPoints);
	      Globatto<Real_wp> lob_reg; lob_reg.AffectPoints(RegPoints);
	      Matrix<Real_wp> ValPhi(order_mesh+1, r-1);
	      for (int i = 0; i <= order_mesh; i++)
		for (int j = 1; j < r; j++)
		  ValPhi(i, j-1) = lob_reg.EvaluatePhi(i, this->lob_curve.Points(j));
              
              // we get nodes for all referenced edges
              VectR_N OldPoints(order_mesh+1);
              for (int iref = 0; iref < this->GetNbBoundaryRef(); iref++)
                {
                  // this is a referenced edge
                  OldPoints(0) = this->Vertex(this->BoundaryRef(iref).numVertex(0));
                  OldPoints(order_mesh) = this->Vertex(this->BoundaryRef(iref).numVertex(1));
                  int ne = NumEntiteBoundary(iref);
		  
                  bool edge_curved = false;
                  for (int k = 1; k < order_mesh; k++)
                    {
                      OldPoints(k) = AllNodes(perm(num_nodes(ne)(k+1)));
                      R_N ptM; Add(1.0-RegPoints(k), OldPoints(0), ptM);
                      Add(RegPoints(k), OldPoints(order_mesh), ptM);
                      if (ptM != OldPoints(k))
                        edge_curved = true;
                    }
                  
                  if (edge_curved)
                    RefUsed(this->BoundaryRef(iref).GetReference()) = true;
                  
                  for (int k = 1; k < r; k++)
                    {
                      this->PointsEdgeRef(iref, k-1).Zero();
                      for (int j = 0; j <= order_mesh; j++)
                        Add(ValPhi(j, k-1), OldPoints(j), this->PointsEdgeRef(iref, k-1));
                    }
                }
              
              // sorting referenced edges
              this->SortBoundariesRef();
              
	      // we find connectivity
	      FindConnectivity();
	      
	      // now we are retrieving internal points of edges if the mesh is curved
	      nb_boundaries = this->GetNbBoundaryRef();
	      this->ResizeBoundariesRef(this->GetNbBoundary());	      
	      for (int i = this->GetNbBoundaryRef(); i < this->GetNbBoundary(); i++)
		{
		  // we retrieve the points of the edge
		  OldPoints(0) = this->Vertex(this->Boundary(i).numVertex(0));
		  OldPoints(order_mesh) = this->Vertex(this->Boundary(i).numVertex(1));
		  int num_elem = this->Boundary(i).numElement(0);
		  int num_loc = this->Element(num_elem).GetPositionBoundary(i);
		  int nv = this->Element(num_elem).GetNbVertices();
		  int ne = NumEntiteElement(num_elem);
		  if (this->Element(num_elem).GetOrientationEdge(num_loc))
		    for (int k = 1; k < order_mesh; k++)
		      OldPoints(k)
			= AllNodes(perm(num_nodes(ne)(nv + (order_mesh-1)*num_loc + k-1)));
		  else
		    for (int k = 1; k < order_mesh; k++)
		      OldPoints(k) = AllNodes(perm(num_nodes(ne)(nv + (order_mesh-1)*num_loc
								 + order_mesh - k - 1)));
		  
		  // the edge is curved ?
		  bool edge_curved = false;
		  for (int k = 1; k < order_mesh; k++)
		    {
		      R_N ptM;
		      Add(1.0-RegPoints(k), OldPoints(0), ptM);
		      Add(RegPoints(k), OldPoints(order_mesh), ptM);
		      if (ptM != OldPoints(k))
			edge_curved = true;
		    }
		  
		  if (edge_curved)
		    {
		      this->BoundaryRef(nb_boundaries) = this->Boundary(i);
		      for (int k = 1; k < r; k++)
			{
			  this->PointsEdgeRef(nb_boundaries, k-1).Zero();
			  for (int j = 0; j <= order_mesh; j++)
			    Add(ValPhi(j, k-1), OldPoints(j),
				this->PointsEdgeRef(nb_boundaries, k-1));
			  
			} 
		      this->BoundaryRef(nb_boundaries).SetReference(new_ref);
		      RefUsed(new_ref) = true;
		      nb_boundaries++;
		    }
		}
	      
	      this->ResizeBoundariesRef(nb_boundaries);
	      for (int i = 1; i < RefUsed.GetM(); i++)
		if (RefUsed(i))		  
		  this->SetCurveType(i, this->CURVE_FILE);
	    }
      }
    else if (!extension_input.compare("mh2"))
      {
	// geompack mesh file
	int nb_vertices;
	file_in>>nb_vertices; int ref;
	this->ReallocateVertices(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  file_in>>this->Vertex(i)(0)>>this->Vertex(i)(1)>>ref;
	
	int nb_extra_vert; double xtmp;
	file_in>>nb_extra_vert;
	for (int i = 0; i < nb_extra_vert; i++)
	  file_in>>ref>>ref>>xtmp;
	
	int nb_nodes, nb_elt; 
	file_in>>nb_nodes>>nb_elt; IVect num(20); num.Fill(-1);
	this->ReallocateElements(nb_elt);
	for (int i = 0; i < nb_elt; i++)
	  {
	    for (int j = 0; j < abs(nb_nodes); j++)
	      {
		file_in>>num(j);
		num(j)--;
	      }
	    
	    if (num(3) == -1)
	      this->Element(i).InitTriangular(num(0), num(1), num(2), 1);
	    else
	      this->Element(i).InitQuadrangular(num(0), num(1), num(2), num(3), 1);
	  }
	
      }
    else if (!extension_input.compare("meshb"))
      {	
	int nb_elt = 0;
	// endian and version are read
	int version, endian, dimension;
	file_in.read(reinterpret_cast<char*>(&endian), sizeof(int));
	file_in.read(reinterpret_cast<char*>(&version), sizeof(int));
	
	// boucle sur tous les mots-cles
	int mot_cle, new_position = 1; int ref; IVect num(8);
	int nb_triangles(0), nb_quadrangles(0);
	while ((!file_in.eof())&&(new_position!=0))
	  {
	    // next keyword is read
	    file_in.read(reinterpret_cast<char*>(&mot_cle), sizeof(int));
	    file_in.read(reinterpret_cast<char*>(&new_position), sizeof(int));
	    
	    if (mot_cle == 3)
	      {
		// dimension
		file_in.read(reinterpret_cast<char*>(&dimension), sizeof(int));
		if (dimension != 2)
		  {
		    cout<<"This is a 3-D mesh, not a 2-D mesh"<<endl;
		    abort();
		  }
	      }
	    else if (mot_cle == 4)
	      {
		// liste des sommets
		int nb_vertices = 0;
		file_in.read(reinterpret_cast<char*>(&nb_vertices), sizeof(int));
		this->ReallocateVertices(nb_vertices);
		TinyVector<float, 2> point;
		for (int i = 0; i < nb_vertices; i++)
		  {
		    // coordonnees
		    file_in.read(reinterpret_cast<char*>(&point(0)), 2*sizeof(float));
		    this->Vertex(i) = point;
		    // reference du point
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		  }
	      }
	    else if (mot_cle == 5)
	      {
		// liste des aretes
		int nb_edges_ref = 0;
		file_in.read(reinterpret_cast<char*>(&nb_edges_ref), sizeof(int));
		this->ReallocateBoundariesRef(nb_edges_ref);
		for (int i = 0; i < nb_edges_ref; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 2*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    this->BoundaryRef(i).Init(num(0)-1, num(1)-1, ref);
		  }
	      }
	    else if (mot_cle == 6)
	      {
		// liste des triangles
		file_in.read(reinterpret_cast<char*>(&nb_triangles), sizeof(int));
		this->ResizeElements(this->GetNbElt() + nb_triangles);
		for (int i = 0; i < nb_triangles; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 3*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    this->Element(nb_elt++).InitTriangular(num(0)-1, num(1)-1, num(2)-1, ref);
		  }
	      }
	    else if (mot_cle == 7)
	      {
		// liste des quadrilateres
		file_in.read(reinterpret_cast<char*>(&nb_quadrangles), sizeof(int));
		this->ResizeElements(this->GetNbElt() + nb_quadrangles);
		for (int i = 0; i < nb_quadrangles; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 4*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    this->Element(nb_elt++).InitQuadrangular(num(0)-1, num(1)-1, num(2)-1,
							     num(3)-1, ref);
		  }
	      }
	    
	    // on passe a la nouvelle position
	    file_in.seekg(new_position, ios::beg);
	  }
      }
    else if (!extension_input.compare("mesh"))
      {
	string keyword_section;
	string ligne; int ref, n0, n1, n2, n3;
	int nb_triangles(0), nb_quadrangles(0);
	while (!file_in.eof())
	  {
	    getline(file_in, keyword_section);
	    if ((keyword_section.size() > 0)&&(keyword_section[0] != '#'))
	      {
		if (keyword_section.find("Vertices") != string::npos)
		  {
		    // all the vertices
		    int nb_vertices = 0;
		    file_in>>nb_vertices;
		    this->ReallocateVertices(nb_vertices);
		    for (int i = 0; i < nb_vertices; i++)
		      file_in>>this->Vertex(i)(0)>>this->Vertex(i)(1)>>ref;
		  }
		else if (keyword_section.find("Edges") != string::npos)
		  {
		    int nb_edges_ref = 0;
		    file_in>>nb_edges_ref;
		    this->ReallocateBoundariesRef(nb_edges_ref);
		    for (int i = 0; i < nb_edges_ref; i++)
		      {
			file_in>>n0>>n1>>ref;
			this->BoundaryRef(i).Init(n0-1, n1-1, ref);
		      }
		  }
		else if (keyword_section.find("Triangles") != string::npos)
		  {
		    file_in>>nb_triangles; 
		    int nb_old_elt = this->GetNbElt();
		    this->ResizeElements(this->GetNbElt() + nb_triangles);
		    for (int i = nb_old_elt; i < this->GetNbElt(); i++)
		      {
			file_in>>n0>>n1>>n2>>ref;
			this->Element(i).InitTriangular(n0-1, n1-1, n2-1, ref);
		      }
		  }
		else if (keyword_section.find("Quadrilaterals") != string::npos)
		  {
		    file_in>>nb_quadrangles; 
		    int nb_old_elt = this->GetNbElt();
		    this->ResizeElements(nb_old_elt + nb_quadrangles);
		    for (int i = nb_old_elt; i < this->GetNbElt(); i++)
		      {
			file_in>>n0>>n1>>n2>>n3>>ref;
			this->Element(i).InitQuadrangular(n0-1, n1-1, n2-1, n3-1, ref);
		      }
		  }
	      }
	  }
      }
    else
      {
	cout<<"Montjoie does not know how to read this format ."<<extension_input<<endl;
	cout<<"Try to use a different mesh format "<<endl;
	abort();
      }
    
    // maximal reference of edges ?
    int ref_max = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      ref_max = max(this->BoundaryRef(i).GetReference(), ref_max);
    
    // resizing arrays if necessary
    if (ref_max > this->GetNbReferences())
      this->ResizeNbReferences(ref_max);
    
    // we change number of vertices for periodic boundaries
    this->SortPeriodicBoundaries();
    
    // we sort referenced boundaries
    this->SortBoundariesRef();
    
    this->ReorientElements();            
    this->FindConnectivity();

    // now we add to edges_ref, the edges which are on the boundary of the domain
    this->AddBoundaryEdges();
    
    this->ProjectPointsOnCurves();    
  }
  
  
  //! writing of a mesh, several formats are available
  /*!
    \param[in] file_name name of the file where the mesh is stored
    \param[in] double_prec if true the mesh is written in double precision
    \param[in] ascii if false the mesh file is written in binary (if possible)
    Some formats (like .vtk) allow the writing of the file in ascii
    format or binary format, and the same extension is used.
  */
  template<class T>
  void Mesh<Dimension2, T>::Write(const string & file_name, bool double_prec, bool ascii) const
  {
    if (this->print_level >= 2)
      cout << " we write mesh on file "<<file_name<<endl;
    
    ofstream file_out(file_name.data());
    if (!file_out.is_open())
      {
	cout << "Unable to open file " << file_name << endl;
	abort();
      }
    
    string extension_output = GetExtension(file_name);
    
    Write(file_out, extension_output, double_prec, ascii);
    
    file_out.close();
  }
  
  
  //! writing of a mesh, several formats are available
  /*!
    \param[in] file_out stream where the mesh is written
    \param[in] extension_output format of the mesh (usually extension of the file)
    \param[in] double_prec mesh written in single or double precision ?
    \param[in] ascii mesh written in ascii or binary format ?
  */  
  template<class T>
  void Mesh<Dimension2, T>::Write(ostream& file_out, const string& extension_output,
				  bool double_prec, bool ascii) const
  {
    if (double_prec)
      file_out.precision(15);
    else
      file_out.precision(7);
    
    file_out.setf(ios_base::scientific, ios_base::floatfield);
    
  
    IVect ref_vertex; this->CreateReferenceVertices(ref_vertex);
    if (!extension_output.compare("t"))
      {
	file_out<<"2  "<<this->GetNbElt()<<"  "
                <<this->GetNbVertices()<<"   "<<this->GetNbVertices()<<'\n';
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertex(i)(0)<<' '<<this->Vertex(i)(1)<<'\n';
	
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    file_out<<this->Element(i).GetNbVertices()<<"  "<<this->Element(i).GetReference()
		    <<"  "<<this->Element(i).GetNbVertices()<<"  ";
	    
	    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	      file_out<<(this->Element(i).numVertex(j)+1)<<"  ";
	     
	    for (int j = 0; j < this->Element(i).GetNbEdges(); j++)
	      {
		int num_edge = this->Element(i).numEdge(j);
		file_out<<abs(this->GetEdge(num_edge).GetReference())<<"  ";
	      }
	    
	    file_out<<endl;
	  }
	
      }
    else if (!extension_output.compare("vtk"))
      {
	file_out << "# vtk DataFile Version 2.0 " << endl;
	file_out << "Montjoie 2-D mesh" << endl;
	if (ascii)
	  file_out << "ASCII" << endl;
	else
	  file_out << "BINARY" << endl;
	
	file_out << "DATASET UNSTRUCTURED_GRID" << endl;
	if (double_prec)
	  file_out << "POINTS " << this->GetNbVertices() << " double " << endl;
	else
	  file_out << "POINTS " << this->GetNbVertices() << " float " << endl;
	
	if (ascii)
	  {
	    for (int i = 0; i < this->GetNbVertices(); i++)
	      file_out << this->Vertex(i)(0) << ' ' << this->Vertex(i)(1) << " 0 \n";
	  }
	else
	  {
	    Vector<double> Pts(3*this->GetNbVertices());
	    for (int i = 0; i < this->GetNbVertices(); i++)
	      {
		Pts(3*i) = toDouble(this->Vertex(i)(0));
		Pts(3*i+1) = toDouble(this->Vertex(i)(1));
		Pts(3*i+2) = 0;
	      }
	    
	    WriteBinaryDoubleOrFloat(Pts, file_out, double_prec, false, true);
	    file_out << '\n';
	  }
	
	int nb_tri = this->GetNbTriangles();
	int nb_quad = this->GetNbElt() - nb_tri;
	int nb_edges_ref = this->GetNbBoundaryRef();
	int size_all = nb_tri*4 + nb_quad*5 + 3*nb_edges_ref;
	file_out << "CELLS " << this->GetNbElt() + nb_edges_ref << " " << size_all << endl;
	if (ascii)
	  {
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      file_out << "2 " << this->BoundaryRef(i).numVertex(0) 
		       << " " << this->BoundaryRef(i).numVertex(1) << '\n';

	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		file_out<<this->Element(i).GetNbVertices()<<' ';
		for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
		  file_out<<(this->Element(i).numVertex(j))<<' ';
		
		file_out<<'\n';
	      }
	  }
	else
	  {
	    IVect num_vert(size_all);
	    int nb = 0;
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      {
		num_vert(nb) = 2;
		num_vert(nb+1) = this->BoundaryRef(i).numVertex(0);
		num_vert(nb+2) = this->BoundaryRef(i).numVertex(1);
		nb += 3;
	      }
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		num_vert(nb++) = this->Element(i).GetNbVertices();
		for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
		  num_vert(nb++) = this->Element(i).numVertex(j);
		
	      }
	    
	    WriteBinaryInteger(num_vert, file_out, false, true);
	    file_out << '\n';
	  }
	
	file_out << "CELL_TYPES " << this->GetNbElt() + nb_edges_ref << endl;
	if (ascii)
	  {
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      file_out << "3\n";
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (this->Element(i).GetNbVertices() == 3)
		  file_out << "5\n";
		else
		  file_out << "9\n";
	      }

	    // writing references in a section CELL_DATA
	    file_out << "CELL_DATA " << this->GetNbElt()+nb_edges_ref << '\n';
	    file_out << "SCALARS cell_scalars int 1" << '\n';
	    file_out << "LOOKUP_TABLE default\n";
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      file_out << this->BoundaryRef(i).GetReference() << '\n';
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      file_out << this->Element(i).GetReference() << '\n';	    
	  }
	else
	  {
	    IVect cell_type(this->GetNbElt() + nb_edges_ref);
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      cell_type(i) = 3;
	    
	    int nb = this->GetNbBoundaryRef();
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (this->Element(i).GetNbVertices() == 3)
		  cell_type(nb+i) = 5;
		else
		  cell_type(nb+i) = 9;
	      }
	    
	    WriteBinaryInteger(cell_type, file_out, false, true);
	    file_out << '\n';

	    // writing references in a section CELL_DATA
	    file_out << "CELL_DATA " << this->GetNbElt()+nb_edges_ref << '\n';
	    file_out << "SCALARS cell_scalars int 1" << '\n';
	    file_out << "LOOKUP_TABLE default\n";
	    nb = 0;
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      cell_type(nb++) = this->BoundaryRef(i).GetReference();
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      cell_type(nb++) = this->Element(i).GetReference();   
	    
	    WriteBinaryInteger(cell_type, file_out, false, true);
	    file_out << '\n';
	  }	
      }
    else if (!extension_output.compare("exo"))
      {
	if (this->GetNbVertices() < 2)
	  return;
	
	// we retrieve nodes of the mesh (different from vertices if the mesh is curved)
	VectR_N PosNodes; Vector<IVect> Nodle;
	int order = GetNodesCurvedMesh(PosNodes, Nodle, 2);
	int nb_nodes = PosNodes.GetM();
	
	// we count referenced edges with respect to their reference
	int ref_max_edge = 0;
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  ref_max_edge = max(ref_max_edge, this->BoundaryRef(i).GetReference());
	
	IVect nb_edges_per_ref(ref_max_edge+1); nb_edges_per_ref.Fill(0);
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  nb_edges_per_ref(this->BoundaryRef(i).GetReference())++;
	
	// we count triangles and quadrilaterals with respect to their references
	int ref_max_elt = 0;
	for (int i = 0; i < this->GetNbElt(); i++)
	  ref_max_elt = max(ref_max_elt, this->Element(i).GetReference());
	
	IVect nb_tri_per_ref(ref_max_elt+1); nb_tri_per_ref.Fill(0);
	IVect nb_quad_per_ref(ref_max_elt+1); nb_quad_per_ref.Fill(0);
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    int ref = this->Element(i).GetReference();
	    if (this->Element(i).GetNbVertices() == 3)
	      nb_tri_per_ref(ref)++;
	    else
	      nb_quad_per_ref(ref)++;
	  }
	
	// number of blocks ?
	int nb_blocks = 0;
	for (int i = 0; i < nb_tri_per_ref.GetM(); i++)
	  if (nb_tri_per_ref(i) > 0)
	    nb_blocks++;
	
	for (int i = 0; i < nb_quad_per_ref.GetM(); i++)
	  if (nb_quad_per_ref(i) > 0)
	    nb_blocks++;
	
	// writing records related to dimensions
	file_out<<"netcdf mesh {\n";
	file_out<<"dimensions:\n";
	file_out<<'\t'<<"len_string = 33 ;\n";
	file_out<<'\t'<<"len_line = 256 ; \n";
	file_out<<'\t'<<"four = 4 ;\n";
	file_out<<'\t'<<"time_step = UNLIMITED ;\n";
	file_out<<'\t'<<"num_dim = 2 ;\n";
	file_out<<'\t'<<"num_nodes = "<<nb_nodes<<" ;\n";
	file_out<<'\t'<<"num_elem = "<<this->GetNbElt()<<" ;\n";
	file_out<<'\t'<<"num_el_blk = "<<nb_blocks<<" ;\n";
	file_out<<'\t'<<"num_qa_rec = 1 ; \n";
	file_out<<'\t'<<"num_node_sets = 1 ;\n";
	file_out<<'\t'<<"num_nod_ns1 = "<<nb_nodes<<" ;\n";
	
	int nb_side_sets = 0;
	for (int ref = 0; ref <= ref_max_edge; ref++)
	  if (nb_edges_per_ref(ref) > 0)
	    nb_side_sets++;
	
	file_out<<'\t'<<"num_side_sets = "<<nb_side_sets<<" ;\n";
	// counting nodes on boundaries
	int nb = 1;
	for (int ref = 0; ref <= ref_max_edge; ref++)
	  if (nb_edges_per_ref(ref) > 0)
	    {
	      Vector<bool> VertexUsed(this->GetNbVertices()); VertexUsed.Fill(false);
	      int nb_vert_ref = 0;
	      for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		if (this->BoundaryRef(i).GetReference() == ref)
		  {
		    for (int k = 0; k < 2; k++)
		      {
			int nv = this->BoundaryRef(i).numVertex(k);
			if (!VertexUsed(nv))
			  {
			    VertexUsed(nv) = true;
			    nb_vert_ref++;
			  }
		      }
		  }
	      
	      file_out<<'\t'<<"num_side_ss"<<nb<<" = "<<nb_edges_per_ref(ref)<<" ;\n"<<endl;
	      file_out<<'\t'<<"num_df_ss"<<nb<<" = "
                      <<nb_vert_ref + (order-1)*nb_edges_per_ref(ref) <<" ;\n"<<endl;
	      nb++;
	    }
	
	// informations for each block
	nb = 1;
	for (int ref = 0; ref < nb_tri_per_ref.GetM(); ref++)
	  if (nb_tri_per_ref(ref) > 0)
	    {
	      file_out<<'\t'<<"num_el_in_blk"<<nb<<" = "<<nb_tri_per_ref(ref)<<" ;\n";
	      file_out<<'\t'<<"num_nod_per_el"<<nb<<" = ";
	      if (order == 1)
		file_out<<"3 ;\n";
	      else
		file_out<<"6 ;\n";
	      nb++;
	    }
	
	for (int ref = 0; ref < nb_quad_per_ref.GetM(); ref++)
	  if (nb_quad_per_ref(ref) > 0)
	    {
	      file_out<<'\t'<<"num_el_in_blk"<<nb<<" = "<<nb_quad_per_ref(ref)<<" ;\n";
	      file_out<<'\t'<<"num_nod_per_el"<<nb<<" = ";
	      if (order == 1)
		file_out<<"4 ;\n";
	      else
		file_out<<"9 ;\n";
	      
	      nb++;
	    }
	
	// variables (very basic, maybe not compatible with all readers
	file_out<<"variables:\n";
	file_out<<"double time_whole(time_step) ;\n";
	
	// data
	file_out<<"data:\n";
	file_out<<'\n'<<" qa_records = \n"<<"  \"CUBIT\",\n"<<"  \"10.2\",\n"
                <<"  \"06/19/2009\",\n"<<"  \"15:26:41\" ;\n";
	file_out<<'\n'<<" coor_names = \n"<<"  \"x\",\n"<<"  \"y\" ;\n";
	file_out<<'\n'<<" ns_status = 1 ;\n\n"<<" ns_prop1 = 1 ;\n\n";
	
	file_out<<" node_ns1 = 1";
	for (int i = 2; i < this->GetNbVertices(); i++)
	  {
	    if (i%20 == 0)
	      file_out<<",\n    "<<i;
	    else
	      file_out<<", "<<i;
	  }
	    
	file_out<<", "<<this->GetNbVertices()<<" ;\n\n";
	
	file_out<<" ss_status = 1";
	for (int i = 1; i < nb_side_sets; i++)
	  file_out<<", 1";
	
	file_out<<" ;\n\n";
	file_out<<" ss_prop1 = 1";
	for (int i = 2; i <= nb_side_sets; i++)
	  file_out<<", "<<i;
	
	file_out<<" ;\n\n";
	nb = 1;
	// writing all the edges on the boundary
	for (int ref = 0; ref <= ref_max_edge; ref++)
	  if (nb_edges_per_ref(ref) > 0)
	    {
	      file_out<<" elem_ss"<<nb<<" = ";
	      int cpt = 1;
	      for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		if (this->BoundaryRef(i).GetReference() == ref)
		  {
		    int num_elem = this->BoundaryRef(i).numElement(0);
		    if (cpt == 1)
		      file_out<<num_elem+1;
		    else
		      {
			if (cpt%20 == 0)
			  file_out<<",\n     "<<num_elem+1;
			else
			  file_out<<", "<<num_elem+1;
		      }
		    cpt++;
		  }
	      
	      file_out<<" ;\n\n";
	      file_out<<" side_ss"<<nb<<" = ";
	      cpt = 1;
	      for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		if (this->BoundaryRef(i).GetReference() == ref)
		  {
		    int num_elem = this->BoundaryRef(i).numElement(0);
		    int num_loc = this->Element(num_elem).GetPositionBoundary(i);
		    if (cpt == 1)
		      file_out<<num_loc+1;
		    else
		      {
			if (cpt%20 == 0)
			  file_out<<",\n     "<<num_loc+1;
			else
			  file_out<<", "<<num_loc+1;
		      }
		    cpt++;
		  }
	      
	      file_out<<" ;\n\n";
	    }
	
	file_out<<" elem_map = 1";
	for (int i = 2; i < this->GetNbElt(); i++)
	  {
	    if (i%20 == 0)
	      file_out<<",\n    "<<i;
	    else
	      file_out<<", "<<i;
	  }
	
	file_out<<", "<<this->GetNbElt()<<" ;\n\n";
	
	file_out<<" eb_status = 1";
	for (int i = 1; i < nb_blocks; i++)
	  file_out<<", 1";
	
	file_out<<" ;\n\n";
	file_out<<" eb_prop1 = 1";
	for (int i = 2; i <= nb_blocks; i++)
	  file_out<<", "<<i;
	
	file_out<<" ;\n\n";
	// writing connectivity of all the elements
	nb = 1;
	for (int ref = 0; ref < nb_tri_per_ref.GetM(); ref++)
	  if (nb_tri_per_ref(ref) > 0)
	    {
	      file_out<<" connec"<<nb<<" = \n";
	      int cpt = 1; TinyVector<int, 3> num;
	      for (int i = 0; i < this->GetNbElt(); i++)
		if ((this->Element(i).GetNbVertices()==3)
                    &&(this->Element(i).GetReference() == ref))
		  {
		    for (int k = 0; k < 3; k++)
		      num(k) = this->Element(i).numVertex(k)+1;
		    
		    if (cpt == 1)
		      {
			file_out<<num(0)<<", "<<num(1)<<", "<<num(2);
		      }
		    else
		      {
			if (cpt%6 == 0)
			  file_out<<",\n    "<<num(0);
			else
			  file_out<<", "<<num(0);
			
			file_out<<", "<<num(1)<<", "<<num(2);
		      }
		    cpt++;
		  }
	      
	      file_out<<" ;\n\n";
	      nb++;
	    }
	
	for (int ref = 0; ref < nb_quad_per_ref.GetM(); ref++)
	  if (nb_quad_per_ref(ref) > 0)
	    {
	      file_out<<" connec"<<nb<<" = \n";
	      int cpt = 1; TinyVector<int, 4> num;
	      for (int i = 0; i < this->GetNbElt(); i++)
		if ((this->Element(i).GetNbVertices()==4)
                    &&(this->Element(i).GetReference() == ref))
		  {
		    for (int k = 0; k < 4; k++)
		      num(k) = this->Element(i).numVertex(k)+1;
		    
		    if (cpt == 1)
		      {
			file_out<<num(0)<<", "<<num(1)<<", "<<num(2)<<", "<<num(3);
		      }
		    else
		      {
			if (cpt%5 == 0)
			  file_out<<",\n    "<<num(0);
			else
			  file_out<<", "<<num(0);
			
			file_out<<", "<<num(1)<<", "<<num(2)<<", "<<num(3);
		      }
		    cpt++;
		  }
	      
	      file_out<<" ;\n\n";
	      nb++;
	    }
	
	// now writing all the nodes
	file_out<<" coord = \n";
	for (int i = 0; i < PosNodes.GetM()-1; i++)
	  file_out<<PosNodes(i)(0)<<", "<<PosNodes(i)(1)<<", \n";
	
	int N = PosNodes.GetM()-1;
	file_out<<PosNodes(N)(0)<<", "<<PosNodes(N)(1)<<" ;\n";
	
	file_out<<"}\n";	      
      }
    else if (!extension_output.compare("geod"))
      {
        file_out << "##FACE##\n";
        file_out << "    0     0\n";
        file_out << std::setw(7) << this->GetNbElt() << "   " << std::setw(7) << this->GetNbVertices() << '\n';
        for (int i = 0; i < this->GetNbVertices(); i++)
          file_out << std::setw(6) << i+1 << "    " << std::fixed << std::setw(22) << std::setprecision(16) << this->Vertex(i)(0) << "    "
                   << std::setw(22) << std::setprecision(16) << this->Vertex(i)(1) << "    "
		   << std::setw(22) << std::setprecision(16) << 0.0 << "\n";
        
        for (int i = 0; i < this->GetNbElt(); i++)
          file_out << std::setw(7) << i+1 << "   " << "  " << std::setw(7) << this->Element(i).numVertex(0)+1
                   << "  " << std::setw(7) << this->Element(i).numVertex(1)+1
                   << "  " << std::setw(7) << this->Element(i).numVertex(2)+1
		   << "   " << std::setw(6) << 18 << "\n";
      }
    else if (!extension_output.compare("msh"))
      {
	VectR_N PosNodes; Vector<IVect> Nodle; 
	int order = GetNodesCurvedMesh(PosNodes, Nodle, this->GetGeometryOrder());
	
	// numbering different in Gmsh and Montjoie
	IVect num_tri, num_quad;
	NumberMap::GetGmshTriangularNumbering(order, num_tri);
	NumberMap::GetGmshQuadrilateralNumbering(order, num_quad);
	
	// sommets
	file_out<<"$MeshFormat\n"<<"2.1 0 8\n"<<"$EndMeshFormat\n"<<"$Nodes\n";
	file_out<<PosNodes.GetM()<<'\n';
	for (int i = 0; i < PosNodes.GetM(); i++)
	  {
	    file_out<<(i+1)<<' '<<PosNodes(i)(0)<<' '<<PosNodes(i)(1)<<" 0.0 \n";
	  }
	
	int nb = this->GetNbBoundaryRef() + this->GetNbElt();
	file_out<<"$EndNodes"<<'\n'<<"$Elements"<<'\n'<<nb<<'\n';
	
	nb = 1;
	// aretes de reference
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  {
	    int type_entite = NumberMap::GetGmshEntityNumber(this->BoundaryRef(i), order);
	    file_out<<nb<<" "<<type_entite<<" 3 "<<this->BoundaryRef(i).GetReference()<<" "
		    << this->BoundaryRef(i).GetGeometryReference() << " 0 ";
	    file_out<<(this->BoundaryRef(i).numVertex(0)+1)<<' '
                    << (this->BoundaryRef(i).numVertex(1)+1);

	    for (int k = 0; k < order-1; k++)
	      file_out << ' ' << this->GetNbVertices() + (order-1)*i + k + 1;
	    
	    file_out<<'\n';
	    nb++;
	  }
	
	// elements
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    file_out<<nb<<' '<<NumberMap::GetGmshEntityNumber(this->Element(i), order)<<" 3 "<<
	      this->Element(i).GetReference()<<" 0 0 ";
	    
	    int nb_vert = this->Element(i).GetNbVertices();
	    if (nb_vert == 3)
	      {
		for (int j = 0; j < nb_vert; j++)
		  file_out<<' '<<Nodle(i)(num_tri(j)) + 1;
		
		for (int j = 0; j < nb_vert; j++)
		  for (int k = 0; k < order-1; k++)
		    file_out << ' ' << Nodle(i)(num_tri(nb_vert + (order-1)*j + k)) + 1;
		
		for (int j = 3*order; j < Nodle(i).GetM(); j++)
		  file_out << ' ' << Nodle(i)(num_tri(j)) + 1;
	      }
	    else
	      {
		for (int j = 0; j < nb_vert; j++)
		  file_out<<' '<<Nodle(i)(num_quad(j)) + 1;
		
		for (int j = 0; j < nb_vert; j++)
		  for (int k = 0; k < order-1; k++)
		    file_out << ' ' << Nodle(i)(num_quad(nb_vert + (order-1)*j + k)) + 1;
		
		for (int j = 4*order; j < Nodle(i).GetM(); j++)
		  file_out << ' ' << Nodle(i)(num_quad(j)) + 1;
	      }
	    
	    file_out<<'\n'; nb++;
	  }
	
	file_out<<"$EndElements"<<endl;
      }
    else if (!extension_output.compare("mh2"))
      {
	// geompack mesh file
	file_out<<this->GetNbVertices()<<'\n';
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertex(i)(0)<<" "<<this->Vertex(i)(1)<<" "<<ref_vertex(i)<<'\n';
	
	file_out<<"0\n";
	
	int nb_nodes = -4;
	if (this->GetNbTriangles() == 0)
	  nb_nodes = 4;
	if (this->GetNbQuadrangles() == 0)
	  nb_nodes = 3;
	file_out<<nb_nodes<<" "<<this->GetNbElt()<<'\n';
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	      file_out<<this->Element(i).numVertex(j)+1<<" ";
	    
	    if ((this->Element(i).GetNbVertices() == 3)&&(nb_nodes == -4))
	      file_out<<"0 ";
	    
	    file_out<<'\n';
	  }
      }
    else if (!extension_output.compare("meshb"))
      {
	int itmp = 1;
	// ecriture endian et version
	file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	
	// ecriture dimension
	itmp = 3; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	itmp = 5*sizeof(int);
	file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	itmp = 2; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	
	// ecriture sommets
	itmp = 4; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	int offset = 8*sizeof(int) + this->GetNbVertices()*(2*sizeof(float) + sizeof(int));
	file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	itmp = this->GetNbVertices(); file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	float xtmp;
	for (int i = 0; i < this->GetNbVertices(); i++)
	  {
	    for (int j = 0; j < 2; j++)
	      {
		xtmp = toDouble(this->Vertex(i)(j));
		file_out.write(reinterpret_cast<char*>(&xtmp), sizeof(float));
	      }
	    itmp = ref_vertex(i);
	    file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	  }
	
	IVect num(8);
	// ecriture aretes
	if (this->GetNbBoundaryRef() > 0)
	  {
	    itmp = 5; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + this->GetNbBoundaryRef()*sizeof(int)*3;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    itmp = this->GetNbBoundaryRef();
	    file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      {
		num(0) = this->BoundaryRef(i).numVertex(0)+1;
		num(1) = this->BoundaryRef(i).numVertex(1)+1;
		file_out.write(reinterpret_cast<char*>(&num(0)), 2*sizeof(int));
		itmp = this->BoundaryRef(i).GetReference();
		file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	      }
	  }
	
	// ecriture triangles
	int nb_triangles = GetNbTriangles();
	if (nb_triangles > 0)
	  {
	    itmp = 6; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + nb_triangles*sizeof(int)*4;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    file_out.write(reinterpret_cast<char*>(&nb_triangles), sizeof(int));
	    for (int i = 0; i < this->GetNbElt(); i++)
	      if (this->Element(i).GetNbVertices() == 3)
		{
		  num(0) = this->Element(i).numVertex(0)+1;
		  num(1) = this->Element(i).numVertex(1)+1;
		  num(2) = this->Element(i).numVertex(2)+1;
		  file_out.write(reinterpret_cast<char*>(&num(0)), 3*sizeof(int));
		  itmp = this->Element(i).GetReference();
		  file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
		}
	  }
	
	// ecriture quadrangles
	int nb_quadrangles = GetNbTriangles();
	if (nb_quadrangles > 0)
	  {
	    itmp = 7; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + nb_quadrangles*sizeof(int)*5;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    file_out.write(reinterpret_cast<char*>(&nb_quadrangles), sizeof(int));
	    for (int i = 0; i < this->GetNbElt(); i++)
	      if (this->Element(i).GetNbVertices() == 4)
		{
		  num(0) = this->Element(i).numVertex(0)+1;
		  num(1) = this->Element(i).numVertex(1)+1;
		  num(2) = this->Element(i).numVertex(2)+1;
		  num(3) = this->Element(i).numVertex(3)+1;
		  file_out.write(reinterpret_cast<char*>(&num(0)), 4*sizeof(int));
		  itmp = this->Element(i).GetReference();
		  file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
		}
	  }
	
	// fin maillage
	itmp = 54; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	itmp = 0; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
      }
    else if (!extension_output.compare("mesh"))
      {
	file_out<<"MeshVersionFormatted 1"<<'\n';
	file_out<<"Dimension "<<'\n'<<"2"<<'\n'<<"Vertices"<<'\n'<<this->GetNbVertices()<<'\n';
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertices(i)(0)<<"  "<<this->Vertices(i)(1)<<"  "<<ref_vertex(i)<<'\n';
	
	if (this->GetNbBoundaryRef() > 0)
	  {
	    file_out<<"Edges"<<'\n'<<this->GetNbBoundaryRef()<<'\n';
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      file_out<<(this->BoundaryRef(i).numVertex(0)+1)<<"  "<<
		(this->BoundaryRef(i).numVertex(1)+1)<<"  "<<
		this->BoundaryRef(i).GetReference()<<'\n';
	  }
	
	int nb_triangles = GetNbTriangles();
	if (nb_triangles > 0)
	  {
	    file_out<<"Triangles"<<'\n'<<nb_triangles<<'\n';
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (this->Element(i).GetNbVertices()==3)
		  {
		    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
		      file_out<<(this->Element(i).numVertex(j)+1)<<"  ";
		    file_out<<this->Element(i).GetReference()<<'\n';
		  }
	      }
	  }
	
	int nb_quadrangles = GetNbQuadrangles();
	if (nb_quadrangles > 0)
	  {
	    file_out<<"Quadrilaterals"<<'\n'<<nb_quadrangles<<'\n';
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (this->Element(i).GetNbVertices()==4)
		  {
		    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
		      file_out<<(this->Element(i).numVertex(j)+1)<<"  ";
		    file_out<<abs(this->Element(i).GetReference())<<'\n';
		  }
	      }
	  }
	
	file_out<<"End"<<'\n';
      }
    else
      {
	cout<<"Montjoie does not know how to write this format ."<<extension_output<<endl;
	cout<<"Try to use a different mesh format "<<endl;
	abort();
      }	
  }

  
  /******************
   * Mesh Treatment *
   ******************/
  
  
  //! Permutation of numbers to have positive jacobians
  /*!
    the local numbers are permuted so that elements are oriented counterclockwise
    By doing that, we have always positive jacobians
  */
  template<class T>
  void Mesh<Dimension2, T>::ReorientElements()
  {    
    IVect num(4); VectR_N s;
    R_N vec_u, vec_v;
    T w; this->affine_element.Reallocate(this->GetNbElt());
    this->reference_element.Reallocate(this->GetNbElt());
    ElementGeomReference<Dimension2>* null_ptr = NULL;
    this->reference_element.Fill(null_ptr);
    this->affine_element.Fill(false);
    for(int i = 0; i < this->GetNbElt(); i++)
      {
	// we get the numbers of the face
	this->GetVerticesElement(i, s);
	int nb_vert = this->Element(i).GetNbVertices();
	for (int j = 0; j < nb_vert; j++)
	  num(j) = this->Element(i).numVertex(j);
	
	// detecting affine elements
	if (this->Element(i).IsCurved())
	  this->affine_element(i) = false;
	else
	  this->affine_element(i) = ElementGeomReference<Dimension2>::IsAffineTransformation(s);
        
	// we compute the vectorial product vec_u ^ vec_v
	// where vec_u is the vector from the first vertex to the second
	// and vec_v from the first to the last
	vec_u = s(1) - s(0);
	vec_v = s(this->Element(i).GetNbVertices()-1) - s(0);
	w = vec_u(0)*vec_v(1) - vec_u(1)*vec_v(0);
	if (w < 0)
	  {
	    // the orientation is not correct
	    // we inverse the orientation
	    Vector<bool> way(nb_vert); IVect num_edge(nb_vert);
            IVect numv(nb_vert);
	    for (int j = 0; j < nb_vert; j++)
	      {
		way(j) = !this->Element(i).GetOrientationEdge(nb_vert-1-j);
		num_edge(j) = this->Element(i).numEdge(nb_vert-1-j);
	      }
	    
            numv(0) = num(0);
	    for (int j = 1; j < nb_vert; j++)
              numv(j) = num(nb_vert-j);
            
            this->Element(i).Init(numv, this->Element(i).GetReference());
	    
	    for (int j = 0; j < nb_vert; j++)
              this->Element(i).SetEdge(j, num_edge(j));
	  }
	
	// setting the pointer to the reference element
	if (nb_vert == 3)
	  this->reference_element(i) = &triangle_reference;
	else
	  this->reference_element(i) = &quadrangle_reference;
      }

    // updating bounding box for the mesh
    this->xmin_ = 1e300; this->xmax_ = -1e300;
    this->ymin_ = 1e300; this->ymax_ = -1e300;
    for (int i = 0; i < this->GetNbVertices(); i++)
      {
	vec_u = this->Vertex(i);
	this->xmin_ = min(this->xmin_, vec_u(0));
	this->xmax_ = max(this->xmax_, vec_u(0));
	this->ymin_ = min(this->ymin_, vec_u(1));
	this->ymax_ = max(this->ymax_, vec_u(1));
      }
  }
  
  
  //! we find boundaries of all the elements (edges in 2-D)
  template<class T>
  void Mesh<Dimension2, T>::FindConnectivity()
  {
    // we delete previous connectivity
    ClearConnectivity();
        
    // we create edges from referenced edges
    int nb_edges = 0;
    this->ReallocateEdges(this->GetNbBoundaryRef());
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      this->AddSingleEdge(nb_edges, this->BoundaryRef(i));
    
    if (nb_edges < this->GetNbEdges())
      this->ResizeEdges(nb_edges);
    
    // then we add inside edges
    this->AddEdgesFromElements();
    
    // relation between referenced edges and edges
    this->ConstructCrossReferenceBoundaryRef();
  }
  
  
  //! we find vertices with the same coordinates, and remove duplicates
  template<class T>
  void Mesh<Dimension2, T>::RemoveDuplicateVertices()
  {
    // vertices with the same positions are removed
    int nb_vertices = this->GetNbVertices();
    if (nb_vertices <= 0)
      return;
    
    IVect permutation(nb_vertices); permutation.Fill();
    // new vertex numbers
    IVect new_num(nb_vertices); new_num.Fill(-1);    
    // vertices are sorted
    Sort(nb_vertices, this->Vertices, permutation);
    
    int nb = 1; R_N prec = this->Vertices(0);
    new_num(permutation(0)) = 0;
    for (int i = 1; i < nb_vertices; i++)
      {
	if (this->Vertices(i) != prec)
	  {
	    this->Vertices(nb) = this->Vertices(i);
	    new_num(permutation(i)) = nb;
	    prec = this->Vertices(nb);
	    nb++;
	  }
	else
	  {
	    // we have found a duplicate, we don't increment nb
	    // so that current vertex will be removed
	    new_num(permutation(i)) = nb-1;
	  }
      }
    
    if (nb != nb_vertices)
      this->ResizeVertices(nb);
    
    // vertex numbers of edges are modified
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      this->BoundaryRef(i).Init(new_num(this->BoundaryRef(i).numVertex(0)),
				new_num(this->BoundaryRef(i).numVertex(1)));
   
    // vertex numbers of elements are modified
    for (int i = 0; i < this->GetNbElt(); i++)
      {
        int nb_vert = this->Element(i).GetNbVertices();
        IVect num(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          num(j) = new_num(this->Element(i).numVertex(j));
        
        this->Element(i).Init(num, this->Element(i).GetReference());
      }
    
    this->SortBoundariesRef();
    
    // we sort boundaries in order to remove duplicate boundaries
    if (this->GetNbBoundaryRef() > 0)
      {
	nb = this->edges_ref.GetM(); permutation.Reallocate(nb); permutation.Fill();
	Sort(nb, this->edges_ref, permutation);
	
	// assembling boundaries 
	nb = 1; Edge<Dimension2> edge_prec = this->BoundaryRef(0);
	Matrix<R_N> OldPoints = this->PointsEdgeRef;
	for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
	  this->PointsEdgeRef(0, k) = OldPoints(permutation(0), k);

	for (int i = 1; i < this->GetNbBoundaryRef(); i++)
	  {
	    if (this->BoundaryRef(i) != edge_prec)
	      {
		this->BoundaryRef(nb) = this->BoundaryRef(i);
		for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
		  this->PointsEdgeRef(nb, k) = OldPoints(permutation(i), k);

		edge_prec = this->BoundaryRef(nb);
		nb++;
	      }
	    else
	      {
		// we have found a duplicate, we don't increment nb
		// so that current boundary will be removed
	      }
	  }
	
	this->ResizeBoundariesRef(nb);
      }
    
    // edges may have changed
    this->FindConnectivity();
    this->ComputeParameterValueOnCurves();
  }
  
  
  //! modification of the mesh in order to force "coherence"
  /*!
    The unused vertices (no elements refer to them) are removed
  */
  template<class T>
  void Mesh<Dimension2, T>::ForceCoherenceMesh(bool keep_edge)
  {
    // we eliminate vertices that are not used by elements
    int nb_vertices = this->GetNbVertices(), nb_elt = this->GetNbElt();
    VectBool VertexUsed(nb_vertices);
    VertexUsed.Fill(false);
    for (int i = 0; i < nb_elt; i++)
      for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	VertexUsed(this->Element(i).numVertex(j)) = true;

    if (keep_edge)
      for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	{
	  VertexUsed(this->BoundaryRef(i).numVertex(0)) = true;
	  VertexUsed(this->BoundaryRef(i).numVertex(1)) = true;
	}
    
    int nb_vertices_unused = 0;
    for (int i = 0; i < nb_vertices; i++)
      if (!VertexUsed(i))
	nb_vertices_unused++;
    
    if (nb_vertices_unused > 0)
      {
	IVect NewNumber(nb_vertices); NewNumber.Fill(-1);
	int nb = 0;
	for (int i = 0; i < nb_vertices; i++)
	  if (VertexUsed(i))
	    NewNumber(i) = nb++;
	
	VectR_N Vertices_tmp(nb_vertices);
	Vertices_tmp = this->Vertices;
	this->ReallocateVertices(nb);
	nb = 0;
	for (int i = 0; i < nb_vertices; i++)
	  if (VertexUsed(i))
	    this->Vertices(nb++) = Vertices_tmp(i);
	
	nb_vertices = nb;
	for (int i = 0; i < nb_elt; i++)
	  {
            int nb_vert = this->Element(i).GetNbVertices();
            IVect num(nb_vert);
            for (int j = 0; j < nb_vert; j++)
              num(j) = NewNumber(this->Element(i).numVertex(j));
            
            this->Element(i).Init(num, this->Element(i).GetReference());
          }
	
	// unactive referenced edges are deleted
	Vector<Edge<Dimension2> > old_edge;
	old_edge = this->edges_ref;
	int nb_edges_ref = 0;
	for (int i = 0; i < old_edge.GetM(); i++)
	  if ((NewNumber(old_edge(i).numVertex(0)) >= 0)&&
	      (NewNumber(old_edge(i).numVertex(1)) >= 0))
	    nb_edges_ref++;
	
	Matrix<R_N> OldPoints = this->PointsEdgeRef;
	Vector<R_N> old_value = this->value_parameter_boundary;
	Vector<TinyVector<int, 2> > old_integer = this->integer_parameter_boundary;
	this->ReallocateBoundariesRef(nb_edges_ref);
	nb_edges_ref = 0;
	for (int i = 0; i < old_edge.GetM(); i++)
	  if ((NewNumber(old_edge(i).numVertex(0)) >= 0)&&
	      (NewNumber(old_edge(i).numVertex(1)) >= 0))
	    {
	      this->BoundaryRef(nb_edges_ref) = old_edge(i);
	      this->BoundaryRef(nb_edges_ref).Init(NewNumber(old_edge(i).numVertex(0)),
						   NewNumber(old_edge(i).numVertex(1)));
	      
	      for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
		this->PointsEdgeRef(nb_edges_ref, k) = OldPoints(i, k);
	      
	      this->value_parameter_boundary(nb_edges_ref) = old_value(i);
	      this->integer_parameter_boundary(nb_edges_ref) = old_integer(i);
	      nb_edges_ref++;
	    }	
	
	// edges may have changed
	this->FindConnectivity();
      }

    if (keep_edge)
      return;
    
    // we eliminate referenced edges not used by elements
    Vector<Edge<Dimension2> > old_edge;
    old_edge = this->edges_ref;
    int nb_edges_ref = 0;
    for (int i = 0; i < old_edge.GetM(); i++)
      if (this->edges_ref(i).GetNbElements() > 0)
        nb_edges_ref++;
    
    if (nb_edges_ref != old_edge.GetM())
      {
	Matrix<R_N> OldPoints = this->PointsEdgeRef;
	Vector<R_N> old_value = this->value_parameter_boundary;
	Vector<TinyVector<int, 2> > old_integer = this->integer_parameter_boundary;
	this->ReallocateBoundariesRef(nb_edges_ref);
	nb_edges_ref = 0;
	for (int i = 0; i < old_edge.GetM(); i++)
	  if (old_edge(i).GetNbElements() > 0)
	    {
	      this->BoundaryRef(nb_edges_ref) = old_edge(i);
	      for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
		this->PointsEdgeRef(nb_edges_ref, k) = OldPoints(i, k);
	      
	      this->value_parameter_boundary(nb_edges_ref) = old_value(i);
	      this->integer_parameter_boundary(nb_edges_ref) = old_integer(i);
	      nb_edges_ref++;
	    }	
        
	// edges may have changed
	this->FindConnectivity();
      }    
  }
  
  
  //! Removes referenced boundaries with zero reference
  template<class T>
  void Mesh<Dimension2, T>::ClearZeroBoundaryRef()
  {
    int nb_edges_ref = 0;
    for (int i = 0; i < this->edges_ref.GetM(); i++)
      if (this->edges_ref(i).GetReference() != 0)
        nb_edges_ref++;
    
    if (nb_edges_ref == this->edges_ref.GetM())
      return;

    // unactive referenced edges are deleted
    Vector<Edge<Dimension2> > old_edge;
    old_edge = this->edges_ref;

    Matrix<R_N> OldPoints = this->PointsEdgeRef;
    Vector<R_N> old_value = this->value_parameter_boundary;
    Vector<TinyVector<int, 2> > old_integer = this->integer_parameter_boundary;
    this->ReallocateBoundariesRef(nb_edges_ref);
    nb_edges_ref = 0;
    for (int i = 0; i < old_edge.GetM(); i++)
      if (old_edge(i).GetReference() != 0)
        {
          this->BoundaryRef(nb_edges_ref) = old_edge(i);
          
          for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
            this->PointsEdgeRef(nb_edges_ref, k) = OldPoints(i, k);
          
          this->value_parameter_boundary(nb_edges_ref) = old_value(i);
          this->integer_parameter_boundary(nb_edges_ref) = old_integer(i);
          nb_edges_ref++;
        }	
    
    // edges may have changed
    this->FindConnectivity();
  }
  
  
  //! local refinement of the mesh near a vertex
  /*!
    \param[in] num_vertex vertex where a refinement is asked
    \param[in] nb_levels depth level of refinement
    \param[in] coef refinement ratio
    works for quadrilateral or triangular elements
  */
  template<class T>
  void Mesh<Dimension2, T>::ApplyLocalRefinement(int num_vertex,
						 int nb_levels, const T& coef)
  {
    // we get all the elements and edges around the given vertex
    IVect ListeElem, ListeEdge;
    for (int num_edge = this->head_minv(num_vertex); num_edge != -1;
	 num_edge = this->next_edge(num_edge))
      {
	ListeEdge.PushBack(num_edge);
	for (int j = 0; j < this->edges(num_edge).GetNbElements(); j++)
	  ListeElem.PushBack(this->edges(num_edge).numElement(j));
      }
    
    for (int num_edge = this->head_maxv(num_vertex); num_edge != -1;
	 num_edge = this->next_edge_arr(num_edge))
      {
	ListeEdge.PushBack(num_edge);
	for (int j = 0; j < this->edges(num_edge).GetNbElements(); j++)
	  ListeElem.PushBack(this->edges(num_edge).numElement(j));
      }
    
    // elimination of duplicate elements and edges
    int nb_elem_near = ListeElem.GetM();
    Assemble(nb_elem_near, ListeElem);
    
    int nb_edges_near = ListeEdge.GetM();
    Assemble(nb_edges_near, ListeEdge);
   
    // for each element near the vertex to refine, we permute numbering
    // so that first vertex is the vertex to refine
    for (int i = 0; i < nb_elem_near; i++)
      {
	int num_elem = ListeElem(i);
	int nb_vert = this->Element(num_elem).GetNbVertices();
	IVect num(nb_vert), num_edge(nb_vert);
	Vector<bool> way_edge(nb_vert);
	int j = -1;
	for (int k = 0; k < nb_vert; k++)
	  {
	    num(k) = this->Element(num_elem).numVertex(k);
	    num_edge(k) = this->Element(num_elem).numEdge(k);
	    way_edge(k) = this->Element(num_elem).GetOrientationEdge(k);
	    if (num(k) == num_vertex)
	      j = k;
	  }
	
	if (j > 0)
	  {
	    // permutation
	    for (int k = 0; k < nb_vert; k++)
	      {
		int kloc = (j+k)%nb_vert;
		this->Element(num_elem).SetVertex(k, num(kloc));
		this->Element(num_elem).SetEdge(k, num_edge(kloc));
              }
	  }
      }
    
    int r = this->GetGeometryOrder();
    Vector<R_N> interm_points(1+2*(r-1));
    Edge<Dimension2> new_edge;

    // counting triangles and quadrangles to refine
    int nb_tri_near = 0, nb_quad_near = 0;
    for (int i = 0; i < nb_elem_near; i++)
      {
	if (this->Element(ListeElem(i)).GetNbVertices() == 3)
	  nb_tri_near++;
	else
	  nb_quad_near++;
      }
    
    // on each edge, we will create a new vertex at the middle
    // on each quadrilateral, a new vertex is created in the center
    int nb_old_vertices = this->GetNbVertices();
    int nb_new_vertices = nb_old_vertices + nb_levels*(nb_edges_near + nb_quad_near);
    this->ResizeVertices(nb_new_vertices);

    // new elements and new edges are created
    int nb_old_edges = this->GetNbEdges();
    int nb_new_edges = nb_old_edges + (nb_edges_near + nb_tri_near*2 + nb_quad_near*3)*nb_levels;
    this->ResizeEdges(nb_new_edges);
    
    int nb_old_edges_ref = this->GetNbBoundaryRef();
    int nb_new_edges_ref = nb_old_edges_ref + (nb_edges_near 
                                               + nb_tri_near*2 + nb_quad_near*3)*nb_levels;
    this->ResizeBoundariesRef(nb_new_edges_ref);
    
    int nb_old_elt = this->GetNbElt();
    int nb_new_elt = nb_old_elt + 2*nb_elem_near*nb_levels;
    this->ResizeElements(nb_new_elt);

    // now we init nb_new... and increment then
    nb_new_elt = nb_old_elt; nb_new_vertices = nb_old_vertices;
    nb_new_edges_ref = nb_old_edges_ref; nb_new_edges = nb_old_edges;
    while (nb_levels > 0)
      {
	// loop over neighboring edges
	for (int i = 0; i < nb_edges_near; i++)
	  {
	    int num_edge = ListeEdge(i);
	    int n1 = this->GetEdge(num_edge).numVertex(0);
	    int n2 = this->GetEdge(num_edge).numVertex(1);
	    int nm = nb_new_vertices + i;
	    int ref = this->GetEdge(num_edge).GetReference();
	    T lambda = 0.5;
	    // with a specific ratio if required
	    if (n2 == num_vertex)
	      lambda = 1.0 - 1.0/coef;
	    else
	      lambda = 1.0/coef;
	    
	    Vector<T> xi(1+2*(r-1)); xi(0) = lambda;
	    for (int k = 0; k < r-1; k++)
	      {
		Real_wp xi_k = this->GetNodalPointInsideEdge(k);
		xi(k+1) = lambda*xi_k;
		xi(k+r) = xi_k*lambda + 1.0 - xi_k;
	      }
	    
	    // computing positions of internal points on this edge
	    if (num_edge < nb_old_edges_ref)
	      this->GetPointsOnCurve(num_edge, ref, xi, interm_points);
	    else
	      interm_points(0) = (1.0-lambda)*this->Vertex(n1) + lambda*this->Vertex(n2);
	    
	    // main vertex
	    this->Vertex(nm) = interm_points(0);
	    
	    // then other points
	    if (n1 == num_vertex)
	      {
		// we change the current edge to the edge close to num_vertex
		this->RemoveSingleEdge(num_edge);
		new_edge.Init(n1, nm, ref); int p = num_edge;
		this->AddSingleEdge(p, new_edge);
		
		// and we add a new edge
		new_edge.Init(n2, nm, ref);
		this->AddSingleEdge(nb_new_edges, new_edge);
		if (num_edge < nb_old_edges_ref )
		  {
		    this->BoundaryRef(num_edge).Init(n1, nm, ref);
		    this->BoundaryRef(nb_new_edges_ref).Init(n2, nm, ref);
		    for (int k = 0; k < r-1; k++)
		      {
			this->SetPointInsideEdge(num_edge, k, interm_points(k+1));
			this->SetPointInsideEdge(nb_new_edges_ref, k, interm_points(k+r));
		      }
		    
		    nb_new_edges_ref++;
		  }
	      }
	    else
	      {
		// we change the current edge to the edge close to num_vertex
		this->RemoveSingleEdge(num_edge);
		new_edge.Init(n2, nm, ref); int p = num_edge;
		this->AddSingleEdge(p, new_edge);
		
		// and we add a new edge
		new_edge.Init(n1, nm, ref);
		this->AddSingleEdge(nb_new_edges, new_edge);
		if (num_edge < nb_old_edges_ref)
		  {
		    this->BoundaryRef(num_edge).Init(n2, nm, ref);
		    this->BoundaryRef(nb_new_edges_ref).Init(n1, nm, ref);
		    for (int k = 0; k < r-1; k++)
		      {
			this->SetPointInsideEdge(num_edge, k, interm_points(k+r));
			this->SetPointInsideEdge(nb_new_edges_ref, k, interm_points(k+1));
		      }
		    nb_new_edges_ref++;
		  }		
	      }
	    
	  }
	
	// loop over elements
	nb_quad_near = 0; nb_tri_near = 0; R2 ptA; Face<Dimension2> new_elt;
	for (int i = 0; i < nb_elem_near; i++)
	  {
	    int num_elem = ListeElem(i);
	    int nb_vert = this->Element(num_elem).GetNbVertices();
	    if (nb_vert == 4)
	      {
		int n1 = this->Element(num_elem).numVertex(0);
		int n2 = this->Element(num_elem).numVertex(1);
		int n3 = this->Element(num_elem).numVertex(2);
		int n4 = this->Element(num_elem).numVertex(3);
		int ref = this->Element(num_elem).GetReference();
		
		int ne0 = this->Element(num_elem).numEdge(0);
		int ne1 = this->Element(num_elem).numEdge(3);
		int e0 = -1, e1 = -1;
		for (int j = 0; j < nb_edges_near; j++)
		  {
		    if (ne0 == ListeEdge(j))
		      e0 = nb_new_vertices + j;
		    
		    if (ne1 == ListeEdge(j))
		      e1 = nb_new_vertices + j;
		  }
		
		// use of Gordon Hall transformation
		T lambda = 1.0/coef;
		
		ptA = (1.0-lambda)*(this->Vertex(e0) + this->Vertex(e1));
		ptA -= (1.0-lambda)*(1.0-lambda)*this->Vertex(n1) - lambda*lambda*this->Vertex(n3);
		int nc = nb_new_vertices + nb_edges_near + nb_quad_near;
		this->Vertex(nc) = ptA;

		// first sub-quadrilateral
		this->Element(num_elem).InitQuadrangular(num_vertex, e0, nc, e1, ref);
		for (int j = 0; j < 4; j++)
		  this->AddEdgeOfElement(nb_new_edges, num_elem, j);
		
		// new sub-quadrilaterals
		new_elt.InitQuadrangular(e0, n2, n3, nc, ref);
		this->AddElement(nb_new_elt, nb_new_edges, new_elt);
		
		new_elt.InitQuadrangular(e1, nc, n3, n4, ref);
		this->AddElement(nb_new_elt, nb_new_edges, new_elt);

		nb_quad_near++;
	      }
	    else
	      {
		// triangle
		int n1 = this->Element(num_elem).numVertex(0);
		int n2 = this->Element(num_elem).numVertex(1);
		int n3 = this->Element(num_elem).numVertex(2);
		
		int ne0 = this->Element(num_elem).numEdge(0);
		int ne1 = this->Element(num_elem).numEdge(2);
		int e0 = -1, e1 = -1;
		for (int j = 0; j < nb_edges_near; j++)
		  {
		    if (ne0 == ListeEdge(j))
		      e0 = nb_new_vertices + j;
		    
		    if (ne1 == ListeEdge(j))
		      e1 = nb_new_vertices + j;
		  }
		
		int ref = this->Element(num_elem).GetReference();
		
		// first sub-triangle
		this->Element(num_elem).InitTriangular(n1, e0, e1, ref);
		for (int j = 0; j < 3; j++)
		  this->AddEdgeOfElement(nb_new_edges, num_elem, j);
		
		// new sub-triangles
		new_elt.InitTriangular(e0, n2, e1, ref);
		this->AddElement(nb_new_elt, nb_new_edges, new_elt);
		
		new_elt.InitTriangular(e1, n2, n3, ref);
		this->AddElement(nb_new_elt, nb_new_edges, new_elt);
		nb_tri_near++;
	      }
	  }
	
	nb_new_vertices += nb_edges_near + nb_quad_near;
	this->ComputeParameterValueOnCurves();
	nb_levels--;
      }
    
    this->ReorientElements();
    this->ResizeBoundariesRef(nb_new_edges_ref);    
    this->FindConnectivity();
    //this->ConstructCrossReferenceBoundaryRef();
    this->ProjectPointsOnCurves();
  }

  
  //! local refinement of the mesh near edges
  /*!
    \param[in] num_vertex vertices where a refinement is asked
    \param[in] nb_levels depth level of refinement
    \param[in] coef refinement ratio
    works for quadrilateral or triangular elements
  */
  template<class T>
  void Mesh<Dimension2, T>::ApplyRefinementEdge(const IVect& ref_cond, int ref_target, bool anisotrope, 
                                                int nb_levels, const T& coef)
  {
    Vector<int> VertexUsed(this->GetNbVertices());
    VertexUsed.Fill(false);
    int nb_edges_to_spare = 0, nb_vertices_to_refine = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      {
        int ref = this->BoundaryRef(i).GetReference();
        if ((ref < ref_cond.GetM()) && (ref_cond(ref) == ref_target))
          {
            if (anisotrope)
              nb_edges_to_spare++;

            int n1 = this->BoundaryRef(i).numVertex(0);
            if (!VertexUsed(n1))
              {
                VertexUsed(n1) = true;
                nb_vertices_to_refine++;
              }

            int n2 = this->BoundaryRef(i).numVertex(1);
            if (!VertexUsed(n2))
              {
                VertexUsed(n2) = true;
                nb_vertices_to_refine++;
              }
          }
      }
    
    Vector<int> edge_to_spare(nb_edges_to_spare);
    Vector<int> num_vertex(nb_vertices_to_refine);
    nb_vertices_to_refine = 0; nb_edges_to_spare = 0;
    VertexUsed.Fill(false);
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      {
        int ref = this->BoundaryRef(i).GetReference();
        if ((ref < ref_cond.GetM()) && (ref_cond(ref) == ref_target))
          {
            if (anisotrope)
              edge_to_spare(nb_edges_to_spare++) = i;

            int n1 = this->BoundaryRef(i).numVertex(0);
            if (!VertexUsed(n1))
              {
                VertexUsed(n1) = true;
                num_vertex(nb_vertices_to_refine++) = n1;
              }

            int n2 = this->BoundaryRef(i).numVertex(1);
            if (!VertexUsed(n2))
              {
                VertexUsed(n2) = true;
                num_vertex(nb_vertices_to_refine++) = n2;
              }
          }
      }
    
    // counting the number of elements/edges close to the vertices to refine
    int nb_elem_near = 0, nb_edges_near = 0;
    Vector<bool> ElemUsed(this->GetNbElt()); ElemUsed.Fill(false);
    Vector<bool> EdgeUsed(this->GetNbEdges()); EdgeUsed.Fill(false);
    for (int k = 0; k < num_vertex.GetM(); k++)
      {
        for (int num_edge = this->head_minv(num_vertex(k)); num_edge != -1;
             num_edge = this->next_edge(num_edge))
          {
            if (!EdgeUsed(num_edge))
              {
                nb_edges_near++;
                EdgeUsed(num_edge) = true;
                for (int j = 0; j < this->edges(num_edge).GetNbElements(); j++)
                  {
                    int ne = this->edges(num_edge).numElement(j);
                    if (!ElemUsed(ne))
                      {
                        nb_elem_near++;
                        ElemUsed(ne) = true;
                      }
                  }
              }
          }
        
        for (int num_edge = this->head_maxv(num_vertex(k)); num_edge != -1;
             num_edge = this->next_edge_arr(num_edge))
          {
            if (!EdgeUsed(num_edge))
              {
                nb_edges_near++;
                EdgeUsed(num_edge) = true;
                for (int j = 0; j < this->edges(num_edge).GetNbElements(); j++)
                  {
                    int ne = this->edges(num_edge).numElement(j);
                    if (!ElemUsed(ne))
                      {
                        nb_elem_near++;
                        ElemUsed(ne) = true;
                      }
                  }
              }
          }
      }
    
    // sparing required edges
    for (int i = 0; i < edge_to_spare.GetM(); i++)
      EdgeUsed(edge_to_spare(i)) = false;
    
    // we get all the elements and edges around the vertices to refine 
    IVect ListeElem(nb_elem_near), ListeEdge(nb_edges_near);
    ListeElem.Zero(); ListeEdge.Zero();
    nb_edges_near = 0; nb_elem_near = 0;
    for (int i = 0; i < this->GetNbEdges(); i++)
      if (EdgeUsed(i))
        ListeEdge(nb_edges_near++) = i;

    for (int i = 0; i < this->GetNbElt(); i++)
      if (ElemUsed(i))
        ListeElem(nb_elem_near++) = i;
    
    int r = this->GetGeometryOrder();
    Vector<R_N> interm_points(1+2*(r-1));

    // counting the number of new points and new elements
    int nb_new_vertices = 0, nb_new_elt = 0;
    Vector<bool> AddVertex(this->GetNbElt()); AddVertex.Fill(false);
    for (int i = 0; i < ListeElem.GetM(); i++)
      {
        int num_elem = ListeElem(i);
        int ne0 = this->Element(num_elem).numEdge(0);
        int ne1 = this->Element(num_elem).numEdge(1);
        int ne2 = this->Element(num_elem).numEdge(2);
        int nb_edges_refine = 0;
        if (EdgeUsed(ne0))
          nb_edges_refine++;

        if (EdgeUsed(ne1))
          nb_edges_refine++;

        if (EdgeUsed(ne2))
          nb_edges_refine++;

        if (this->Element(num_elem).GetNbVertices() == 3)
          {
            if (nb_edges_refine == 3)
              nb_new_elt += 3;
            else
              nb_new_elt++;
          }
        else
          {
            int ne3 = this->Element(num_elem).numEdge(3);
            if (EdgeUsed(ne3))
              nb_edges_refine++;

            bool edge_opp = false;
            if (EdgeUsed(ne0) && EdgeUsed(ne2))
              edge_opp = true;

            if (EdgeUsed(ne1) && EdgeUsed(ne3))
              edge_opp = true;

            switch (nb_edges_refine)
              {
              case 1 : nb_new_elt++; break;
              case 2: 
                if (edge_opp)
                  nb_new_elt++;
                else
                  { nb_new_elt += 2; nb_new_vertices++; AddVertex(num_elem) = true; }
                break;
              case 3:
              case 4:
                nb_new_elt += 3; nb_new_vertices++; AddVertex(num_elem) = true; break;
              }
          }
      }
    
    // counting edges
    nb_new_vertices += nb_edges_near;
    int nb_new_edges = 0;
    Vector<int> IndexEdge(this->GetNbBoundary());
    IndexEdge.Fill(-1);
    for (int i = 0; i < nb_edges_near; i++)
      {
        IndexEdge(ListeEdge(i)) = i;
        if (ListeEdge(i) < this->GetNbBoundaryRef())
          nb_new_edges++;
      }
    
    // resizing arrays
    int nb_old_vertices = this->GetNbVertices();
    this->ResizeVertices(nb_new_vertices + nb_old_vertices);

    int nb_old_edges = this->GetNbBoundaryRef();
    this->ResizeBoundariesRef(nb_new_edges + nb_old_edges);

    int nb_old_elt = this->GetNbElt();
    this->ResizeElements(nb_new_elt + nb_old_elt);

    // loop over edges to refine
    VectReal_wp EdgeRatio(nb_edges_near);
    nb_new_vertices = nb_old_vertices; nb_new_edges = nb_old_edges;
    for (int i = 0; i < nb_edges_near; i++)
      {
        int num_edge = ListeEdge(i);
        int n1 = this->GetEdge(num_edge).numVertex(0);
        int n2 = this->GetEdge(num_edge).numVertex(1);
        int ref = this->GetEdge(num_edge).GetReference();
        T lambda = 0.5;
        // with a specific ratio if required
        if (VertexUsed(n2) && (!VertexUsed(n1)))
          lambda = 1.0 - 1.0/coef;
        
        if (VertexUsed(n1) && (!VertexUsed(n2)))
          lambda = 1.0/coef;
        
        EdgeRatio(i) = lambda;
        Vector<T> xi(1+2*(r-1)); xi(0) = lambda;
        for (int k = 0; k < r-1; k++)
          {
            Real_wp xi_k = this->GetNodalPointInsideEdge(k);
            xi(k+1) = lambda*xi_k;
            xi(k+r) = xi_k*lambda + 1.0 - xi_k;
          }
        
        // computing positions of internal points on this edge
        if (num_edge < nb_old_edges)
          this->GetPointsOnCurve(num_edge, ref, xi, interm_points);
        else
          interm_points(0) = (1.0-lambda)*this->Vertex(n1) + lambda*this->Vertex(n2);
        
        // new vertex
        int nm = nb_new_vertices;
        this->Vertex(nb_new_vertices) = interm_points(0);
        nb_new_vertices++;
        
        // then other points
        if (num_edge < nb_old_edges)
          {
            this->BoundaryRef(num_edge).Init(n1, nm, ref);
            this->BoundaryRef(nb_new_edges).Init(n2, nm, ref);
            for (int k = 0; k < r-1; k++)
              {
                this->SetPointInsideEdge(num_edge, k, interm_points(k+1));
                this->SetPointInsideEdge(nb_new_edges, k, interm_points(k+r));
              }
            
            nb_new_edges++;
          }
      }

    // loop over elements to refine
    nb_new_elt = nb_old_elt;
    for (int i = 0; i < nb_elem_near; i++)
      {
        int num_elem = ListeElem(i);
        int nb_vert = this->Element(num_elem).GetNbVertices();
        int n1 = this->Element(num_elem).numVertex(0);
        int n2 = this->Element(num_elem).numVertex(1);
        int n3 = this->Element(num_elem).numVertex(2);

        int ne0 = this->Element(num_elem).numEdge(0);
        int ne1 = this->Element(num_elem).numEdge(1);
        int ne2 = this->Element(num_elem).numEdge(2);
        int ref = this->Element(num_elem).GetReference();

        int e0 = nb_old_vertices + IndexEdge(ne0);
        int e1 = nb_old_vertices + IndexEdge(ne1);
        int e2 = nb_old_vertices + IndexEdge(ne2);
        
        if (nb_vert == 4)
          {
            int ne3 = this->Element(num_elem).numEdge(3);
            int n4 = this->Element(num_elem).numVertex(3);
            int e3 = nb_old_vertices + IndexEdge(ne3);
            int nm = -1;
            
            // use of Gordon Hall transformation to compute the new vertex
            if (AddVertex(num_elem))
              {
                R2 ptE0, ptE1, ptE2, ptE3;
                T x(0.5), y(0.5);
                
                if (EdgeUsed(ne0))
                  {
                    ptE0 = this->Vertex(e0);
                    if (this->Element(num_elem).GetOrientationEdge(0))
                      x = EdgeRatio(IndexEdge(ne0));
                    else
                      x = 1.0-EdgeRatio(IndexEdge(ne0));
                  }
                else
                  ptE0 = (1.0-x)*this->Vertex(n1) + x*this->Vertex(n2);

                if (EdgeUsed(ne1))
                  {
                    ptE1 = this->Vertex(e1);
                    if (this->Element(num_elem).GetOrientationEdge(1))
                      y = EdgeRatio(IndexEdge(ne1));
                    else
                      y = 1.0-EdgeRatio(IndexEdge(ne1));
                  }
                else
                  ptE1 = (1.0-y)*this->Vertex(n2) + y*this->Vertex(n3);

                if (EdgeUsed(ne2))
                  {
                    ptE2 = this->Vertex(e2);
                    if (this->Element(num_elem).GetOrientationEdge(2))
                      x = 1.0-EdgeRatio(IndexEdge(ne2));
                    else
                      x = EdgeRatio(IndexEdge(ne2));
                  }
                else
                  ptE2 = (1.0-x)*this->Vertex(n4) + x*this->Vertex(n3);

                if (EdgeUsed(ne3))
                  {
                    ptE3 = this->Vertex(e3);
                    if (this->Element(num_elem).GetOrientationEdge(3))
                      y = 1.0-EdgeRatio(IndexEdge(ne3));
                    else
                      y = EdgeRatio(IndexEdge(ne3));
                  }
                else
                  ptE3 = (1.0-y)*this->Vertex(n1) + y*this->Vertex(n4);
                
                R2 ptA = (1.0-x)*ptE3 + x*ptE1 + (1.0-y)*ptE0 + y*ptE2;
                ptA -= (1.0-x)*(1.0-y)*this->Vertex(n1) + x*(1.0-y)*this->Vertex(n2) 
                  + x*y*this->Vertex(n3) + (1.0-x)*y*this->Vertex(n4);
                
		this->Vertex(nb_new_vertices) = ptA;
                nm = nb_new_vertices;
                nb_new_vertices++;
              }

            if (EdgeUsed(ne0))
              {
                // ne0 refined
                if (EdgeUsed(ne1))
                  {
                    // ne0, ne1 refined 
                    if (EdgeUsed(ne2))
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e0, nm, e3, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e0, n2, e1, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e1, n3, e2, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e2, n4, e3, nm, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(e0, n2, e1, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e1, n3, e2, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(n1, nm, e2, n4, ref);
                            this->Element(nb_new_elt++).InitTriangular(n1, e0, nm, ref);
                          }
                      }
                    else
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e0, nm, e3, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e0, n2, e1, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e3, nm, n3, n4, ref);
                            this->Element(nb_new_elt++).InitTriangular(e1, n3, nm, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e0, nm, n4, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e0, n2, e1, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e1, n3, n4, nm, ref);
                          }
                      }
                  }
                else
                  {
                    // ne0 refined , ne1 not refined 
                    if (EdgeUsed(ne2))
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e0, nm, e3, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e3, nm, e2, n4, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(n2, n3, e2, nm, ref);
                            this->Element(nb_new_elt++).InitTriangular(e0, n2, nm, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e0, e2, n4, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e0, n2, n3, e2, ref);
                          }
                      }
                    else
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e0, nm, e3, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e0, n2, n3, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e3, nm, n3, n4, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(e0, n2, n3, n4, ref);
                            this->Element(nb_new_elt++).InitTriangular(n1, e0, n4, ref);
                          }
                      }
                  }
              } // end if EdgeUsed(ne0)
            else
              {
                // cases where ne0 is not refined
                if (EdgeUsed(ne1))
                  {
                    // ne0 not refined, ne1 refined 
                    if (EdgeUsed(ne2))
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, n2, e1, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e1, n3, e2, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e3, nm, e2, n4, ref);
                            this->Element(nb_new_elt++).InitTriangular(n1, nm, e3, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(n1, n2, e1, nm, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(nm, e1, n3, e2, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(n1, nm, e2, n4, ref);
                          }
                      }
                    else
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, n2, e1, e3, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e1, n3, n4, e3, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(n1, e1, n3, n4, ref);
                            this->Element(nb_new_elt++).InitTriangular(n1, n2, e1, ref);
                          }
                      }
                  }
                else
                  {
                    // ne0 not refined, ne1 not refined
                    if (EdgeUsed(ne2))
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(n1, n2, nm, e3, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(nm, n2, n3, e2, ref);
                            this->Element(nb_new_elt++).InitQuadrangular(e3, nm, e2, n4, ref);
                          }
                        else
                          {
                            this->Element(num_elem).InitQuadrangular(n1, n2, n3, e2, ref);
                            this->Element(nb_new_elt++).InitTriangular(n1, e2, n4, ref);
                          }
                      }
                    else
                      {
                        if (EdgeUsed(ne3))
                          {
                            this->Element(num_elem).InitQuadrangular(e3, n2, n3, n4, ref);
                            this->Element(nb_new_elt++).InitTriangular(n1, n2, e3, ref);
                          }
                      }
                  }
              }
          }
        else
          {
            // triangular case
            if (EdgeUsed(ne0))
              {
                // ne0 refined
                if (EdgeUsed(ne1))
                  {
                    // ne0, ne1 refined 
                    if (EdgeUsed(ne2))
                      {
                        this->Element(num_elem).InitTriangular(n1, e0, e2, ref);
                        this->Element(nb_new_elt++).InitTriangular(e0, n2, e1, ref);
                        this->Element(nb_new_elt++).InitTriangular(e2, e1, n3, ref);
                        this->Element(nb_new_elt++).InitTriangular(e0, e1, e2, ref);
                      }
                    else
                      {
                        this->Element(num_elem).InitTriangular(e0, n1, e1, ref);
                        this->Element(nb_new_elt++).InitQuadrangular(n1, e0, e1, n3, ref);
                      }
                  }
                else
                  {
                    // ne0 refined, n1 not refined
                    if (EdgeUsed(ne2))
                      {
                        this->Element(num_elem).InitTriangular(n1, e0, e2, ref);
                        this->Element(nb_new_elt++).InitQuadrangular(e0, n2, n3, e2, ref);
                      }
                    else
                      {
                        this->Element(num_elem).InitTriangular(n1, e0, n3, ref);
                        this->Element(nb_new_elt++).InitTriangular(e0, n2, n3, ref);
                      }
                  }
              }
            else
              {
                // ne0 not refined
                if (EdgeUsed(ne1))
                  {
                    // ne0 not refined , ne1 refined 
                    if (EdgeUsed(ne2))
                      {
                        this->Element(num_elem).InitTriangular(e2, e1, n3, ref);
                        this->Element(nb_new_elt++).InitQuadrangular(n1, n2, e1, e2, ref);
                      }
                    else
                      {
                        this->Element(num_elem).InitTriangular(n1, n2, e1, ref);
                        this->Element(nb_new_elt++).InitTriangular(n1, e1, n3, ref);
                      }
                  }
                else
                  {
                    // ne0 refined, ne1 not refined
                    if (EdgeUsed(ne2))
                      {
                        this->Element(num_elem).InitTriangular(n1, n2, e2, ref);
                        this->Element(nb_new_elt++).InitTriangular(e2, n2, n3, ref);
                      }
                  }
              }
          }
      }
    
    this->ReorientElements();
    this->FindConnectivity();
    //this->ConstructCrossReferenceBoundaryRef();
    this->ProjectPointsOnCurves();

    // recursive call
    if (nb_levels > 0)
      ApplyRefinementEdge(ref_cond, ref_target, anisotrope, nb_levels-1, coef);
  }


  //! adds a vertex to the mesh and creates elements such that the mesh stays conformal
  template<class T>
  void Mesh<Dimension2, T>::AddVertexConformal(const TinyVector<T, 2>& pt)
  {
    // searching the element where the point is
    Real_wp dist_min = this->Vertex(0).Distance(pt);
    int nv = 0;
    for (int i = 1; i < this->GetNbVertices(); i++)
      {
	Real_wp dist = this->Vertex(i).Distance(pt);
	if (dist < dist_min)
	  {
	    dist_min = dist;
	    nv = i;
	  }
      }
    
    if (dist_min <= TinyVector<T, 2>::threshold)
      {
        // the point already exists
        return;
      }
    
    IVect neighbor_elt;
    this->GetNeighboringElementsAroundVertex(nv, neighbor_elt);
    int nb_elt = neighbor_elt.GetM();
    
    // then loop over potential neighbor elements
    R2 pt_loc;
    for (int i = 0; i < nb_elt; i++)
      {
	int ne = neighbor_elt(i);
	FjInverseProblem<Dimension2> inverseFj(*this, ne);
        bool pt_inside = inverseFj.Solve(pt, pt_loc);
        if (pt_inside)
          {
            int nb_vert = this->GetNbVertices();
            int nb_elt = this->GetNbElt();
            int ref = this->Element(ne).GetReference();
            
            Real_wp dist_edge = this->GetDistanceToBoundary(pt_loc, ne);
            
            this->ResizeVertices(nb_vert+1);            
            this->Vertex(nb_vert) = pt;
            
            if (dist_edge <= TinyVector<T, 2>::threshold)
              {
                int num_loc = this->ProjectPointOnBoundary(pt_loc, ne);
                int num_edge = this->Element(ne).numEdge(num_loc);
                int ref_edge = this->GetEdge(num_edge).GetReference();
                if (this->Type_curve(ref_edge) >= 1)
                  {
                    cout << "Case not handled" << endl;
                    abort();
                  }

		int nb_edges_ref = this->GetNbBoundaryRef();
                if (num_edge < nb_edges_ref)
                  {
                    int n0 = this->BoundaryRef(num_edge).numVertex(0);
                    int n1 = this->BoundaryRef(num_edge).numVertex(1);
                    this->ResizeBoundariesRef(nb_edges_ref+1);
                    this->BoundaryRef(num_edge).Init(n0, nb_vert, ref_edge);
                    this->BoundaryRef(nb_edges_ref).Init(n1, nb_vert, ref_edge);
                  }
                  
                for (int n = 0; n < this->GetEdge(num_edge).GetNbElements(); n++)
                  {
                    int num_elem = this->GetEdge(num_edge).numElement(n);
                    int ref = this->Element(num_elem).GetReference();
                    int num_loc2 = this->Element(num_elem).GetPositionBoundary(num_edge);
                    if (this->Element(num_elem).GetNbVertices() == 3)
                      {
                        int n0 = this->Element(num_elem).numVertex(0);
                        int n1 = this->Element(num_elem).numVertex(1);
                        int n2 = this->Element(num_elem).numVertex(2);
                        
                        this->ResizeElements(nb_elt+1);
                        switch (num_loc2)
                          {
                          case 0 :
                            this->Element(num_elem).InitTriangular(n0, nb_vert, n2, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n1, n2, ref);
                            break;
                          case 1 :
                            this->Element(num_elem).InitTriangular(n1, nb_vert, n0, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n2, n0, ref);
                            break;
                          case 2 :
                            this->Element(num_elem).InitTriangular(n2, nb_vert, n1, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n0, n1, ref);
                            break;
                          }
                        
                        nb_elt++;
                      }
                    else
                      {
                        int n0 = this->Element(num_elem).numVertex(0);
                        int n1 = this->Element(num_elem).numVertex(1);
                        int n2 = this->Element(num_elem).numVertex(2);
                        int n3 = this->Element(num_elem).numVertex(3);

                        this->ResizeElements(nb_elt+2);
                        
                        switch (num_loc2)
                          {
                          case 0 :
                            this->Element(num_elem).InitTriangular(n0, nb_vert, n3, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n1, n2, ref);
                            this->Element(nb_elt+1).InitTriangular(nb_vert, n2, n3, ref);
                            break;
                          case 1 :
                            this->Element(num_elem).InitTriangular(n1, nb_vert, n0, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n2, n3, ref);
                            this->Element(nb_elt+1).InitTriangular(nb_vert, n3, n0, ref);
                            break;
                          case 2 :
                            this->Element(num_elem).InitTriangular(n2, nb_vert, n1, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n3, n0, ref);
                            this->Element(nb_elt+1).InitTriangular(nb_vert, n0, n1, ref);
                            break;
                          case 3 :
                            this->Element(num_elem).InitTriangular(n3, nb_vert, n2, ref);
                            this->Element(nb_elt).InitTriangular(nb_vert, n0, n1, ref);
                            this->Element(nb_elt+1).InitTriangular(nb_vert, n1, n2, ref);
                            break;
                          }
                        
                        nb_elt += 2;
                      }
                  }
              }
            else
              {    
                int n0 = this->Element(ne).numVertex(0);
                int n1 = this->Element(ne).numVertex(1);
                int n2 = this->Element(ne).numVertex(2);
                
                if (nb_vert == 3)
                  {
                    this->ResizeElements(nb_elt+2);
                    this->Element(ne).InitTriangular(n0, n1, nb_vert, ref);
                    this->Element(nb_elt).InitTriangular(n1, n2, nb_vert, ref);
                    this->Element(nb_elt+1).InitTriangular(n2, n0, nb_vert, ref);
                  }
                else
                  {
                    this->ResizeElements(nb_elt+3);
                    int n3 = this->Element(ne).numVertex(3);
                    this->Element(ne).InitTriangular(n0, n1, nb_vert, ref);
                    this->Element(nb_elt).InitTriangular(n1, n2, nb_vert, ref);
                    this->Element(nb_elt+1).InitTriangular(n2, n3, nb_vert, ref);
                    this->Element(nb_elt+2).InitTriangular(n3, n0, nb_vert, ref);
                  }
              }
            
            break;
          }
      }
    
    this->ReorientElements();
    this->FindConnectivity();
    this->ComputeParameterValueOnCurves();    
  }
  
  
  //! The mesh is symmetrized with respected to x-axis or y-axis
  /*!
    the axis can also be given by normale and xc if num_coor = -1
  */
  template<class T>
  void Mesh<Dimension2, T>::SymmetrizeMesh(int num_coor, T xc,
                                           TinyVector<T, 2> normale, bool keep_edge_on_plane)
  {
    // normalizing normale
    T coef = Norm2(normale);
    if (num_coor == -1)
      {
        coef = T(1)/coef;
        normale *= coef;
        xc *= coef;
      }
    
    int nb_elt = this->GetNbElt();
    this->ResizeElements(2*nb_elt);
    int nb_old_edges = this->GetNbBoundaryRef();
    // counting vertices
    int nb_old_vertices = this->GetNbVertices();
    int nb_vertices = nb_old_vertices;
    // IndexVert(i) -> new vertex number
    IVect IndexVert(nb_vertices); IndexVert.Fill(-1);
    for (int i = 0; i < nb_old_vertices; i++)
      {
	// if vertex on axis of symmetry, we keep the same number
	// (no new vertex created)
        if (num_coor == -1)
          {
            if (abs(DotProd(this->Vertices(i), normale) + xc) > R_N::threshold)
              IndexVert(i) = nb_vertices++;
            else
              IndexVert(i) = i;
          }
        else
          {
            if (abs(this->Vertex(i)(num_coor)-xc) >= R_N::threshold)
              IndexVert(i) = nb_vertices++;
            else
              IndexVert(i) = i;
          }
      }
    
    // creating symmetrized vertices
    this->ResizeVertices(nb_vertices);
    for (int i = 0; i < nb_old_vertices; i++)
      {
	this->Vertices(IndexVert(i)) = this->Vertices(i);
	if (IndexVert(i) != i)
          {
            if (num_coor == -1)
              SymmetrizePointPlane(this->Vertices(i), normale, xc, this->Vertices(IndexVert(i)));
            else
              this->Vertices(IndexVert(i))(num_coor) = 2.0*xc - this->Vertices(i)(num_coor);
          }
      }
    
    // creating symmetrized elements
    for (int i = 0; i < nb_elt; i++)
      {
	int nv = this->Element(i).GetNbVertices();
	IVect num(nv);
	for (int j = 0; j < nv; j++)
	  num(j) = IndexVert(this->Element(i).numVertex(j));
	
	this->Element(nb_elt+i).Init(num, this->Element(i).GetReference());
      }
    
    Matrix<R_N> OldPointsEdge = this->PointsEdgeRef;
    Vector<Edge<Dimension2> > old_edges = this->edges_ref;
    this->ReallocateBoundariesRef(2*nb_old_edges);
    int nb_edges_ref = 0; R_N point;
    for (int i = 0; i < nb_old_edges; i++)
      {
	int n0 = old_edges(i).numVertex(0);
	int n1 = old_edges(i).numVertex(1);
        int ref = old_edges(i).GetReference();
	int n0_sym = IndexVert(n0);
	int n1_sym = IndexVert(n1);
        if ((n0_sym != n0)||(n1_sym != n1))
	  {
            this->BoundaryRef(nb_edges_ref).Init(n0, n1, ref);
            
	    for (int k = 0; k < this->GetGeometryOrder()-1; k++)
              this->PointsEdgeRef(nb_edges_ref, k) = OldPointsEdge(i, k);
            
	    nb_edges_ref++;
            
	    this->BoundaryRef(nb_edges_ref).Init(n0_sym, n1_sym, ref);
            
	    for (int k = 0; k < this->GetGeometryOrder()-1; k++)
	      {
		point = OldPointsEdge(i, k);
                if (num_coor == -1)
                  SymmetrizePointPlane(OldPointsEdge(i, k), normale, xc, point);
                else
                  point(num_coor) = 2.0*xc - point(num_coor);
                
		this->SetPointInsideEdge(nb_edges_ref, k, point);
	      }
            
	    nb_edges_ref++;
	  }
        else
          {
            if (keep_edge_on_plane)
              {
                this->BoundaryRef(nb_edges_ref).Init(n0, n1, ref);
                
                for (int k = 0; k < this->GetGeometryOrder()-1; k++)
                  this->PointsEdgeRef(nb_edges_ref, k) = OldPointsEdge(i, k);
                
                nb_edges_ref++;
              }
          }
	
      }
    
    this->ResizeBoundariesRef(nb_edges_ref);
    this->SortBoundariesRef();
    
    this->ReorientElements();
    this->FindConnectivity();
    this->ComputeParameterValueOnCurves();
  }

  
  //! mesh is periodized
  template<class T>
  void Mesh<Dimension2, T>::PeriodizeMeshTeta(const R_N& center)
  {
    T teta, teta_min, teta_max;
    GetTetaMinMax(center, teta_min, teta_max);
    
    // number of periodicity 
    teta = teta_max - teta_min;
    int nper = toInteger(round(2.0*pi_wp/teta));
    
    R_N vec_u;
    int nb_old_vert = this->GetNbVertices();
    int nb_old_elt = this->GetNbElt();
    int nb_old_edges = this->GetNbBoundaryRef();
    int nb_vertices = nb_old_vert, nb_elt = nb_old_elt, nb_edges_ref = nb_old_edges;
    this->ResizeVertices(nb_vertices*nper);
    this->ResizeElements(nb_elt*nper);
    this->ResizeBoundariesRef(nb_edges_ref*nper);
    T cos_nteta, sin_nteta; IVect num; R_N point;
    // loop over each new periodized cell
    for (int i = 1; i < nper; i++)
      {
	cos_nteta = cos(teta*i); sin_nteta = sin(teta*i);
	// new vertices
	for (int j = 0; j < nb_old_vert; j++)
	  {
	    vec_u = this->Vertex(j) - center;
	    this->Vertex(nb_vertices + j)(0) = center(0) + cos_nteta*vec_u(0)
	      - sin_nteta*vec_u(1);
	    this->Vertex(nb_vertices + j)(1) = center(1) + sin_nteta*vec_u(0) 
	      + cos_nteta*vec_u(1);
	  }
	
	// new elements
	for (int j = 0; j < nb_old_elt; j++)
	  {
	    int nv = this->Element(j).GetNbVertices(); num.Reallocate(nv);
	    for (int k = 0; k < nv; k++)
	      num(k) = this->Element(j).numVertex(k) + nb_vertices;
	    
	    this->Element(nb_elt + j).Init(num, this->Element(j).GetReference());
	  }
	
	// new edges
	for (int j = 0; j < nb_old_edges; j++)
	  {
	    int n0 = this->BoundaryRef(j).numVertex(0) + nb_vertices;
	    int n1 = this->BoundaryRef(j).numVertex(1) + nb_vertices;
	    
	    this->BoundaryRef(nb_edges_ref + j).Init(n0, n1, this->BoundaryRef(j).GetReference());
	    for (int k = 0; k < this->GetGeometryOrder()-1; k++)
	      {
		vec_u = this->PointsEdgeRef(j, k) - center;
		point(0) = center(0) + cos_nteta*vec_u(0) - sin_nteta*vec_u(1);
		point(1) = center(1) + sin_nteta*vec_u(0) + cos_nteta*vec_u(1);
		this->PointsEdgeRef(nb_edges_ref+j, k) = point;
	      }
	  }
	
	nb_vertices += nb_old_vert;
	nb_elt += nb_old_elt;
	nb_edges_ref += nb_old_edges;
      }
    
    this->ReorientElements();
    this->RemoveDuplicateVertices();
    this->ComputeParameterValueOnCurves();
  }

  
  //! the mesh is split in order to have only quadrilaterals
  /*!
    Each triangle is cut in three quadrilaterals and each quadrilateral in four
  */
  template<class T>
  void Mesh<Dimension2, T>::SplitIntoQuadrilaterals()
  {
    if (this->print_level >= 2)
      cout << " we cut each triangle in three quadrilaterals "
           << "and each quadrilateral in four"<<endl;
    
    ReorientElements(); 
    R_N point;
    
    // we place a new point on each middle of edge, and face
    int nb_old_vertices = this->GetNbVertices();
    int nb_old_edges = this->GetNbEdges();
    int nb_old_edges_ref = this->GetNbBoundaryRef();
    int nb_old_elt = this->GetNbElt();
    int nb_new_vertices = nb_old_vertices + nb_old_edges + nb_old_elt;
    this->ResizeVertices(nb_new_vertices);
    
    int r = this->GetGeometryOrder();
    VectReal_wp xi(1+2*(r-1));
    Vector<R_N> interm_points(1+2*(r-1));
    xi(0) = 0.5;
    for (int i = 0; i < r-1; i++)
      {
	xi(i+1) = 0.5*this->GetNodalPointInsideEdge(i);
	xi(i+r) = 1.0 - xi(i+1);
      }
    
    this->ResizeBoundariesRef(2*nb_old_edges_ref);
    
    // first we divide all the edges
    for (int i = 0; i < nb_old_edges; i++)
      {
	int n1 = this->GetEdge(i).numVertex(0);
	int n2 = this->GetEdge(i).numVertex(1);
	if (i < nb_old_edges_ref)
	  {
	    int ref = this->GetEdge(i).GetReference();
	    this->GetPointsOnCurve(i, ref, xi, interm_points);
	    this->Vertex(nb_old_vertices + i) = interm_points(0);
	    this->BoundaryRef(i).Init(n1, nb_old_vertices+i, ref); 
	    this->BoundaryRef(i+nb_old_edges_ref).Init(n2, nb_old_vertices+i, ref); 
	    for (int k = 0; k < r-1; k++)
	      {
		this->SetPointInsideEdge(i, k, interm_points(k+1));
		this->SetPointInsideEdge(i+nb_old_edges_ref, k, interm_points(k+r));
	      }
	  }
	else
	  {
	    point = T(0.5)*(this->Vertex(n1) + this->Vertex(n2));
	    this->Vertex(nb_old_vertices + i) = point;
	  }
      }
    
    VectR_N s(4), f_edge(4);
    int nb_new_elt = 3*this->GetNbTriangles() + 4*this->GetNbQuadrangles();
    this->ResizeElements(nb_new_elt);
    nb_new_elt = nb_old_elt;
    
    // then we split elements
    for (int i = 0; i < nb_old_elt; i++)
      {
	int nb_vert = this->Element(i).GetNbVertices();
	this->GetVerticesElement(i, s);
	IVect num(nb_vert), num_edge(nb_vert);
	for (int j = 0; j < nb_vert; j++)
	  {
	    num(j) = this->Element(i).numVertex(j);
	    num_edge(j) = this->Element(i).numEdge(j);
	    f_edge(j) = this->Vertex(nb_old_vertices + num_edge(j));
	    num_edge(j) += nb_old_vertices;
	  }
	
	int nc = nb_old_vertices + nb_old_edges + i;
	int ref = this->Element(i).GetReference();
	if (nb_vert == 3)
	  {
	    point = T(-1.0)/9*( s(0) + s(1) + s(2))
	      + T(4.0)/9*(f_edge(0) + f_edge(1) + f_edge(2));
	    
	    this->Vertex(nc) = point;
	    
	    this->Element(i).
	      InitQuadrangular(num(0), num_edge(0), nc, num_edge(2), ref);
	    
	    this->Element(nb_new_elt++)
	      .InitQuadrangular(num_edge(0), num(1), num_edge(1), nc, ref);
	    
	    this->Element(nb_new_elt++)
	      .InitQuadrangular(nc, num_edge(1), num(2), num_edge(2), ref);
	    
	  }
	else if (nb_vert == 4)
	  {
	    point = T(-0.25)*( s(0) + s(1) + s(2) + s(3))
	      + T(0.5)*(f_edge(0) + f_edge(1) + f_edge(2) + f_edge(3));
	    
	    this->Vertex(nc) = point;

	    this->Element(i).
	      InitQuadrangular(num(0), num_edge(0), nc, num_edge(3), ref);
	    
	    this->Element(nb_new_elt++)
	      .InitQuadrangular(num_edge(0), num(1), num_edge(1), nc, ref);
	    
	    this->Element(nb_new_elt++)
	      .InitQuadrangular(nc, num_edge(1), num(2), num_edge(2), ref);

	    this->Element(nb_new_elt++)
	      .InitQuadrangular(num_edge(3), nc, num_edge(2), num(3), ref);
	  }
	
      }
    
    this->ReorientElements();            
    this->FindConnectivity();
    this->ComputeParameterValueOnCurves();
  }
  
  
  //! the mesh is split in order to have only triangles
  /*!
    Each quadrilateral is split into two triangles
  */
  template<class T>
  void Mesh<Dimension2, T>::SplitIntoTriangles()
  {
    if (this->print_level >= 2)
      cout << " we cut each quadrangle in two triangles "
           << "and each triangle in four"<<endl;
    
    ReorientElements(); 
    
    int nb_new_elt = this->GetNbTriangles() + 2*this->GetNbQuadrangles();
    int nb_old_edges = this->GetNbEdges();
    this->ResizeEdges(nb_old_edges + this->GetNbQuadrangles());
    int nb_old_elt = this->GetNbElt();
    this->ResizeElements(nb_new_elt);    
    nb_new_elt = nb_old_elt;
    Face<Dimension2> e;
    
    for (int i = 0; i < nb_old_elt; i++)
      {
	int nb_vertices_loc = this->Element(i).GetNbVertices();
	
	if (nb_vertices_loc==4)
	  {
	    int n1 = this->Element(i).numVertex(0);
	    int n2 = this->Element(i).numVertex(1);
	    int n3 = this->Element(i).numVertex(2);
	    int n4 = this->Element(i).numVertex(3);
	    int ref = this->Element(i).GetReference();
	    int nb = i;
	    e.InitTriangular(n1, n2, n4, ref);
	    this->AddElement(nb, nb_old_edges, e);
	    
	    e.InitTriangular(n2, n3, n4, ref);
	    this->AddElement(nb_new_elt, nb_old_edges, e);
	  }
      }
    
    this->ReorientElements();            
    this->FindConnectivity();
    this->ConstructCrossReferenceBoundaryRef();
    this->ComputeParameterValueOnCurves();
  }
  
  
  //! returns all the nodes of the mesh and numbering
  /*!
    \param[in] outside_subdiv intern points, including extremities 
    (ie 0  0.2 0.5 0.6 0.8  1.0, for a regular subdivision in 5 intervals)
    \param[out] Nodle mesh numbering
  */
  template<class T>
  void Mesh<Dimension2, T>
  ::GetMeshSubdivision(const Vector<T>& outside_subdiv, const Vector<Vector<R_N> >& points_div,
                       Vector<R_N>& PosNodes, Vector<IVect>& Nodle, Vector<IVect>& NodleSurf) const
  {
    if (outside_subdiv.GetM() < 2)
      return;
    
    if (outside_subdiv.GetM() == 2)
      {
	PosNodes = this->Vertices;
	Nodle.Reallocate(this->GetNbElt());
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    int nv = this->Element(i).GetNbVertices();
	    Nodle(i).Reallocate(nv);
	    for (int j = 0; j < nv; j++)
	      Nodle(i)(j) = this->Element(i).numVertex(j);	    
	  }
	
	NodleSurf.Reallocate(this->GetNbBoundaryRef());
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  {
	    NodleSurf(i).Reallocate(this->BoundaryRef(i).GetNbVertices());
	    for (int j = 0; j < this->BoundaryRef(i).GetNbVertices(); j++)
	      NodleSurf(i)(j) = this->BoundaryRef(i).numVertex(j);
	  }
	
	return;
      }

    // step_subdiv is an array containing all the intermediary points on an edge
    Vector<T> step_subdiv(outside_subdiv.GetM()-2);
    for (int i = 0; i < step_subdiv.GetM(); i++)
      step_subdiv(i) = outside_subdiv(i+1);    
    
    int r = step_subdiv.GetM() + 1;
        
    // counting the number of nodes
    int nb_nodes = this->GetNbVertices() + this->GetNbEdges()*(r-1)
      + this->GetNbTriangles()*(r-1)*(r-2)/2
      + this->GetNbQuadrangles()*(r-1)*(r-1);
    
    PosNodes.Reallocate(nb_nodes);
    SetPoints<Dimension2> PointsElem;
    VectR_N s, InterpPoints;

    // constructing interpolation to points_div
    FiniteElementInterpolator interp;
    
    // constructing 2-D points
    interp.InitProjection(this->reference_element, outside_subdiv, points_div);
    
    for (int i = 0; i < this->GetNbVertices(); i++)
      PosNodes(i) = this->Vertex(i);
 
    // 1-D interpolation
    int order = this->GetGeometryOrder();
    Matrix<Real_wp> ValPhi(order+1, r-1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j < r-1; j++)
	ValPhi(i, j) = this->lob_curve.EvaluatePhi(i, step_subdiv(j));
    
    // computation of nodes for all the referenced edges
    VectR2 OldPoints(order+1);;
    NodleSurf.Reallocate(this->GetNbBoundaryRef());
    int offset = this->GetNbVertices();
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      {
	NodleSurf(i).Reallocate(r+1);
	NodleSurf(i)(0) = this->BoundaryRef(i).numVertex(0);
	NodleSurf(i)(r) = this->BoundaryRef(i).numVertex(1);

	for (int k = 1; k < r; k++)
	  NodleSurf(i)(k) = offset + i*(r-1) + k-1;
	
	int ref = this->BoundaryRef(i).GetReference();
	OldPoints(0) = this->Vertex(this->BoundaryRef(i).numVertex(0));
	OldPoints(order) = this->Vertex(this->BoundaryRef(i).numVertex(1));
	int offset = this->GetNbVertices() + i*(r-1);
	if (this->GetCurveType(ref) > 0)
	  {
	    for (int k = 1; k < order; k++)
	      OldPoints(k) = this->PointsEdgeRef(i, k-1);
	    	    
	    for (int k = 0; k < r-1; k++)
	      {
		PosNodes(offset + k).Zero();
		for (int j = 0; j <= order; j++)
		  Add(ValPhi(j, k), OldPoints(j), PosNodes(offset+k));
	      }
	  }
	else
	  {
	    for (int k = 0; k < r-1; k++)
	      {
		PosNodes(offset + k).Zero();
		Add(1.0-step_subdiv(k), OldPoints(0), PosNodes(offset+k));
		Add(step_subdiv(k), OldPoints(order), PosNodes(offset+k));
	      }  
	  }
      }
    
    // construction of Nodle and computation of nodes
    Nodle.Reallocate(this->GetNbElt());
    int nb = offset + (r-1)*this->GetNbEdges();
    for (int i = 0; i < this->GetNbElt(); i++)
      {	
	int nb_dof = 0;
	int nb_vert = this->Element(i).GetNbVertices();
	if (nb_vert == 3)
	  nb_dof = (r+1)*(r+2)/2;
	else
	  nb_dof = (r+1)*(r+1);
	
	Nodle(i).Reallocate(nb_dof);
	// vertices
	for (int j = 0; j < nb_vert; j++)
	  Nodle(i)(j) = this->Element(i).numVertex(j);
	
	// then edges
	for (int j = 0; j < nb_vert; j++)
	  {
	    int ne = this->Element(i).numEdge(j);
	    if (this->Element(i).GetOrientationEdge(j))
	      for (int k = 0; k < r-1; k++)
		Nodle(i)(nb_vert + j*(r-1) + k) = offset + ne*(r-1) + k;
	    else
	      for (int k = 0; k < r-1; k++)
		Nodle(i)(nb_vert + j*(r-1) + k) = offset + ne*(r-1) + r - 2 - k;
	  }
	
	// the interior
	for (int k = nb_vert*r; k < nb_dof; k++)
	  Nodle(i)(k) = nb++;
	
	// computing the points
	this->GetVerticesElement(i, s);
	this->FjElemNodal(s, PointsElem, *this, i);
	
	// interpolation
	interp.Project(PointsElem.GetPointNodal(), InterpPoints, this->GetTypeElement(i));
	
	for (int k = 0; k < nb_dof; k++)
	  PosNodes(Nodle(i)(k)) = InterpPoints(k);
	
      }
  }  
  
  
  //! the mesh is split into sub-quads/sub-triangles
  /*!
    \param[in] step_subdiv intern points, including extremities 
    (ie 0  0.2 0.5 0.6 0.8  1.0, for a regular subdivision in 5 intervals)
  */
  template<class T>
  void Mesh<Dimension2, T>::SubdivideMesh(const Vector<T>& step_subdiv)
  {
    Vector<IVect> Nodle, NodleSurf; Vector<VectR_N> points_div;
    Vector<Vector<T> > points_surf;
    SubdivideMesh(step_subdiv, points_surf, points_div, Nodle, NodleSurf);
  }    
  
  
  //! the mesh is split into sub-quads/sub-triangles
  /*!
    \param[in] outside_subdiv intern points, including extremities 
    (ie 0  0.2 0.5 0.6 0.8  1.0, for a regular subdivision in 5 intervals)
    \param[out] Nodle mesh numbering of the split mesh
    \param[out] points_div "split" points on unit triangle and unit square
  */
  template<class T>
  void Mesh<Dimension2, T>::SubdivideMesh(const Vector<T>& outside_subdiv,
                                          Vector<Vector<T> >& points_surf,
					  Vector<VectR_N>& points_div, Vector<IVect>& Nodle,
					  Vector<IVect>& NodleSurf)
  {
    if (outside_subdiv.GetM() < 2)
      return;
    
    if ((abs(outside_subdiv(0)) > R_N::threshold)||
	(abs(1.0-outside_subdiv(outside_subdiv.GetM()-1)) > R_N::threshold))
      {
	cout << "Extremities should be contained in the provided subdivision " <<endl;
	abort();
      }
    
    if (outside_subdiv.GetM() == 2)
      {
	NodleSurf.Reallocate(this->GetNbBoundaryRef());
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  {
	    NodleSurf(i).Reallocate(this->BoundaryRef(i).GetNbVertices());
	    for (int j = 0; j < this->BoundaryRef(i).GetNbVertices(); j++)
	      NodleSurf(i)(j) = this->BoundaryRef(i).numVertex(j);
	  }

	Nodle.Reallocate(this->GetNbElt());
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    Nodle(i).Reallocate(this->Element(i).GetNbVertices());
	    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	      Nodle(i)(j) = this->Element(i).numVertex(j);
	  }
	
	points_div.Reallocate(2);
	points_div(0).Reallocate(3); points_div(1).Reallocate(4);
	    
	points_div(0)(0)(0) = 0.0; points_div(0)(0)(1) = 0.0;
	points_div(0)(1)(0) = 1.0; points_div(0)(1)(1) = 0.0;
	points_div(0)(2)(0) = 0.0; points_div(0)(2)(1) = 1.0;
	points_div(1)(0)(0) = 0.0; points_div(1)(0)(1) = 0.0;
	points_div(1)(1)(0) = 1.0; points_div(1)(1)(1) = 0.0;
	points_div(1)(2)(0) = 1.0; points_div(1)(2)(1) = 1.0;
	points_div(1)(3)(0) = 0.0; points_div(1)(3)(1) = 1.0;
	return;
      }
    
    if (this->print_level >= 4)
      cout<< " we divide the initial mesh "<<endl;
    
    Vector<T> step_subdiv(outside_subdiv.GetM()-2);
    for (int i = 0; i < step_subdiv.GetM(); i++)
      step_subdiv(i) = outside_subdiv(i+1);
    
    // step_subdiv is an array containing all the intermediary points to add on a edge
    // exemple step_subdiv={0.5} -> we will split each quad in 4 sub-quad
    // and each triangle in 4 sub-triangles
    int nb_div = step_subdiv.GetM();
    Vector<T> Points1D(nb_div+2);
    points_div.Reallocate(2);
    Vector<Matrix<int> > NumNodes2D(2);
    Matrix<int> CoordinateNodes;
    
    int r = nb_div+1;
    
    Points1D(0) = T(0);  Points1D(r) = T(1);
    for (int i = 0; i < nb_div; i++)
      Points1D(i+1) = step_subdiv(i);
    
    // we get NumNodes2D for triangle and quadrilateral
    MeshNumbering<Dimension2>::
      ConstructTriangularNumbering(nb_div+1, NumNodes2D(0), CoordinateNodes);
    
    points_div(0).Reallocate((r+1)*(r+2)/2);
  
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(nb_div+1, NumNodes2D(1), CoordinateNodes);
    
    points_div(1).Reallocate((r+1)*(r+1));
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  points_div(1)(NumNodes2D(1)(i,j)).Init(Points1D(i),Points1D(j));
	  if (NumNodes2D(0)(i,j) != -1)
	    points_div(0)(NumNodes2D(0)(i,j)).Init(Points1D(i), Points1D(j));
	}
    
    // now we compute all the nodes and fill Nodle
    VectR_N PosNodes;
    // reconstructing edges (in the case where it is not done)
    GetMeshSubdivision(outside_subdiv, points_div, PosNodes, Nodle, NodleSurf);
    
    int nb_old_vert = this->GetNbVertices();
    this->ReallocateVertices(PosNodes.GetM());
    for (int i = 0; i < PosNodes.GetM(); i++)
      this->Vertex(i) = PosNodes(i);
    
    PosNodes.Clear();
    
    int nb_old_elt = this->GetNbElt();
    int nb_new_elt = nb_old_elt*r*r;
    Vector<Face<Dimension2> > OldElement = this->elements;
    this->ReallocateElements(nb_new_elt);
    nb_new_elt = 0;
    
    for (int i = 0; i < nb_old_elt; i++)
      {
	int ref = OldElement(i).GetReference();
	IVect num(Nodle(i).GetM());
	for (int j = 0; j < Nodle(i).GetM(); j++)
	  num(j) = Nodle(i)(j);
	
	// now we create small elements
	if (OldElement(i).GetNbVertices() == 3)
	  {
	    // a triangle is subdivided into small triangles
	    for (int j = 0; j <= nb_div; j++)
	      for (int k = 0; k <= nb_div; k++)
		{
		  if ((NumNodes2D(0)(j+1,k)!=-1)&(NumNodes2D(0)(j,k+1)!=-1))
		    {
		      this->Element(nb_new_elt++).
			InitTriangular(num(NumNodes2D(0)(j,k)),
				       num(NumNodes2D(0)(j+1,k)),
				       num(NumNodes2D(0)(j,k+1)), ref);
		      
		      if (NumNodes2D(0)(j+1,k+1)!=-1)
			this->Element(nb_new_elt++).
			  InitTriangular(num(NumNodes2D(0)(j+1,k)),
					 num(NumNodes2D(0)(j+1,k+1)),
					 num(NumNodes2D(0)(j,k+1)), ref);
		    }
		}
	  }
	else
	  for (int j = 0; j <= nb_div; j++)
	    for (int k = 0; k <= nb_div; k++)
	      {
		this->Element(nb_new_elt++).
		  InitQuadrangular(num(NumNodes2D(1)(j,k)),
				   num(NumNodes2D(1)(j+1,k)),
				   num(NumNodes2D(1)(j+1,k+1)),
				   num(NumNodes2D(1)(j,k+1)), ref);
	      }
      }
    
    int nb_old_edges_ref = this->GetNbBoundaryRef();
    int nb_new_edges_ref = nb_old_edges_ref*r;
    int rg = this->lob_curve.GetOrder(); 
    // we compute internal points on referenced edges (curved mesh)
    Matrix<R_N> NewPoints;
    if (rg > 1)
      {
	NewPoints.Reallocate(nb_new_edges_ref, rg-1);
	Vector<R_N> OldPos(rg+1);
        Matrix<Real_wp> ValPhi(rg+1, r*(rg-1)); VectReal_wp Xi(r*(rg-1));
	for (int i = 0; i < r; i++)
	  for (int j = 1; j < rg; j++)
	    {
	      int num_point = i*(rg-1) + j-1;
	      Real_wp lambda = (1.0-this->lob_curve.Points(j))*outside_subdiv(i) + 
		this->lob_curve.Points(j)*outside_subdiv(i+1);
	      if (i == r-1)
		lambda = this->lob_curve.Points(j)*outside_subdiv(i) + 
		  (1.0-this->lob_curve.Points(j))*outside_subdiv(i+1);
	      
	      Xi(num_point) = lambda;
	      for (int k = 0; k <= rg; k++)
		ValPhi(k, num_point) = this->lob_curve.EvaluatePhi(k, lambda);
	    }
	
	for (int i = 0; i < nb_old_edges_ref; i++)
	  {
	    OldPos(0) = this->Vertex(this->BoundaryRef(i).numVertex(0));
	    OldPos(rg) = this->Vertex(this->BoundaryRef(i).numVertex(1));
	    for (int k = 1; k < rg; k++)
	      OldPos(k) = this->PointsEdgeRef(i, k-1);
	    
	    int ref = this->BoundaryRef(i).GetReference();
	    if (this->GetCurveType(ref) != this->CURVE_FILE )
	      {
		// predefined curve, we compute the image of intermediary points
		Vector<R_N> points_interm;
		this->GetPointsOnCurve(i, ref, Xi, points_interm);
		for (int is = 0; is < r; is++)
		  for (int j = 1; j < rg; j++)
		    {
		      int num_point = is*(rg-1) + j-1;
		      int ne = i*r + is;
		      NewPoints(ne, j-1) = points_interm(num_point);
		    }
	      }
	    else
	      {
		// curved mesh, we perform an interpolation
		for (int is = 0; is < r; is++)
		  for (int j = 1; j < rg; j++)
		    {
		      int num_point = is*(rg-1) + j-1;
		      int ne = i*r + is;
		      NewPoints(ne, j-1).Zero();
		      for (int k = 0; k <= rg; k++) 
			Add(ValPhi(k, num_point), OldPos(k), NewPoints(ne, j-1));
		    }
	      }
	  }
      }
    
    // creation of new referenced edges
    Vector<Edge<Dimension2> > OldEdgeRef = this->edges_ref;

    this->ReallocateBoundariesRef(nb_new_edges_ref);
    IVect num(r+1); nb_new_edges_ref = 0;
    for (int i = 0; i < nb_old_edges_ref; i++)
      {
	int ref = OldEdgeRef(i).GetReference();	
	num(0) = OldEdgeRef(i).numVertex(0);
	num(r) = OldEdgeRef(i).numVertex(1);
	for (int j = 1; j < r; j++)
	  num(j) = nb_old_vert + i*(r-1) + j-1;
	
	for (int j = 0; j <= nb_div; j++)
	  {
	    if (num(j) < num(j+1))
	      this->BoundaryRef(nb_new_edges_ref++).Init(num(j), num(j+1), ref);
	    else
	      this->BoundaryRef(nb_new_edges_ref++).Init(num(j+1), num(j), ref);
	  }
      }
    
    if (rg > 1)
      this->PointsEdgeRef = NewPoints;    
    
    this->ReorientElements();
    FindConnectivity();
    this->ProjectPointsOnCurves();
     
  }

  
  //! if the mesh is with negative x, a symmetry is applied to have only x > 0
  /*!
    This method is used for axisymmetric computations
    \warning if a mesh has both negative x and positive x,
    the number of elements is set to 0
  */
  template<class T>
  void Mesh<Dimension2, T>::SymmetrizeMeshToGetPositiveX()
  {
    int nb_vert_plus = 0; int nb_vert_minus = 0;
    for (int i = 0; i < this->GetNbVertices(); i++)
      {
	if (abs(this->Vertex(i)(0)) <= R_N::threshold)
	  this->Vertex(i)(0) = T(0);
	else if (this->Vertex(i)(0) < -R_N::threshold)
	  nb_vert_minus++;
	else
	  nb_vert_plus++;
      }
    
    if ((nb_vert_minus >0)&&(nb_vert_plus>0))
      {
	cout << "Mesh crosses y-axis" << endl;
	abort();
      }
    
    if (nb_vert_minus > 0)
      {
	for (int i = 0; i < this->GetNbVertices(); i++)
	  this->Vertex(i)(0) = -this->Vertex(i)(0);
	
	for (int i = 0; i < this->PointsEdgeRef.GetM(); i++)
	  for (int j = 0; j < this->PointsEdgeRef.GetN(); j++)
	    this->PointsEdgeRef(i, j)(0) *= -1.0;
	
	this->ReorientElements();
	this->ComputeParameterValueOnCurves();
      }
    
  }
  
  
  /**********************
   * Creation of meshes *
   **********************/
  
  
  //! creation of a regular mesh of a rectangular domain [xmin xmax ymin ymax]
  /*!
    \param[in] ptA point [xmin ymin]
    \param[in] ptB point [xmax ymax]
    \param[in] nbPoints_x the number of points along the edge xmin -> xmax
    \param[in] ref_domain reference of the inside domain
    \param[in] ref_boundary  reference of south, east,
    north and west side of the rectangle
    \param[in] triangular_mesh if true, we do a mesh in triangles
  */
  template<class T>
  void Mesh<Dimension2, T>::
  CreateRegularMesh(const R_N& ptA, const R_N& ptB, const TinyVector<int, 2>& nbPoints,
		    int ref_domain, const TinyVector<int, 4>& ref_boundary,
		    int type_mesh, R2 ratio)
  {
    // previous mesh is cleared
    Clear();
        
    Real_wp xmin = ptA(0); Real_wp ymin = ptA(1);
    Real_wp xmax = ptB(0); Real_wp ymax = ptB(1);
    
    if ((xmax <= xmin)||(ymax <= ymin))
      return;

    bool geom_x = false;
    if (ratio(0) != 1.0)
      geom_x = true;
    
    bool geom_y = false;
    if (ratio(1) != 1.0)
      geom_y = true;
    
    int nbPoints_x = nbPoints(0), nbPoints_y = nbPoints(1);
    T step_x, step_y;
    if (nbPoints_x<=1)
      {
	step_x = xmax - xmin;
	nbPoints_x = 2;
      }
    else
      step_x = (xmax - xmin)/(nbPoints_x-1); 
    
    if (nbPoints_y<=1)
      {
	step_y = ymax - ymin;
	nbPoints_y = 2;
      }
    else
      step_y = (ymax - ymin)/(nbPoints_y-1); 
    
    if (geom_x)
      step_x = (xmax-xmin)*(1.0-ratio(0))/(1.0 - pow(ratio(0), nbPoints_x-1));
    
    if (geom_y)
      step_y = (ymax-ymin)*(1.0-ratio(1))/(1.0 - pow(ratio(1), nbPoints_y-1));
    
    T pos_y, pos_x;
    bool triangular_mesh = (type_mesh == this->TRIANGULAR_MESH);
    if (type_mesh == this->RADAU_MESH)
      triangular_mesh = true;
    
    int nb_elt = 0;
    if (triangular_mesh)
      nb_elt = 2*(nbPoints_y-1)*(nbPoints_x-1);
    else
      nb_elt = (nbPoints_y-1)*(nbPoints_x-1); 
    
    int nb_vertices = nbPoints_y*nbPoints_x; 
    this->ReallocateVertices(nb_vertices);
    this->ReallocateElements(nb_elt);
    this->ReallocateBoundariesRef(2*(nbPoints_x+nbPoints_y-2));
    int ne = 0;
    nb_elt = 0; int nb_edges_ref = 0;
    bool radau = (type_mesh == this->RADAU_MESH);
    // loop on y
    for (int j = 0; j < nbPoints_y; j++)
      {
	if (j == (nbPoints_y-1))
	  pos_y = ymax;
	else
	  {
	    if (geom_y)
	      {
		if (j == 0)
		  pos_y = ymin;
		else
		  pos_y = ymin + step_y*(1.0-pow(ratio(1), j))/(1.0-ratio(1));
	      }
	    else
	      pos_y = ymin + step_y*j;
	  }
	
	// loop on x
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    if (i==(nbPoints_x-1))
	      pos_x = xmax;
	    else
	      {
		if (geom_x)
                  {
                    if (i == 0)
                      pos_x = xmin;
                    else
                      pos_x = xmin + step_x*(1.0-pow(ratio(0), i))/(1.0-ratio(0));
                  }
                else
		  pos_x = xmin + step_x*i;
	      }
	    
	    this->Vertex(ne).Init(pos_x, pos_y);
	    // inside point -> we create an element
	    if ((i<nbPoints_x-1)&&(j<nbPoints_y-1)) 
	      {
		if (triangular_mesh)
		  {
		    if (radau)
		      {
                        int n0 = ne, n1 = ne+1, n2 = ne+nbPoints_x, n3 = ne+nbPoints_x;
			this->Element(nb_elt++).
			  InitQuadrangular(n0, n1, n2, n3, ref_domain);
			
                        n0 = ne+1; n1 = ne+nbPoints_x+1; n2 = ne+nbPoints_x; n3 = n2;
			this->Element(nb_elt).
			  InitQuadrangular(n0, n1, n2, n3, ref_domain);

		      }
		    else
		      {
                        int n0 = ne, n1 = ne+1, n2 = ne+nbPoints_x;
			this->Element(nb_elt++).
			  InitTriangular(n0, n1, n2, ref_domain);
			
                        n0 = ne+1; n1 = ne+nbPoints_x+1; n2 = ne+nbPoints_x;
			this->Element(nb_elt).
			  InitTriangular(n0, n1, n2, ref_domain);
		      }
		  }
		else
		  {
                    int n0 = ne, n1 = ne+1, n2 = ne+nbPoints_x+1, n3 = ne+nbPoints_x;
                    this->Element(nb_elt).
                      InitQuadrangular(n0, n1, n2, n3, ref_domain); 
                  }
		
		nb_elt++;
	      }
	    
	    if ((i==0)&&(j<nbPoints_y-1))
	      {
		// point on east side
		this->edges_ref(nb_edges_ref).Init(ne, ne+nbPoints_x, ref_boundary(3));
		nb_edges_ref++;
	      }
	    
	    if ((j==0)&&(i<nbPoints_x-1)) 
	      {
		// point on south side
		this->edges_ref(nb_edges_ref).Init(ne, ne+1, ref_boundary(0));
		nb_edges_ref++;
	      }
	    
	    if ((i==nbPoints_x-1)&&(j<nbPoints_y-1))
	      {
		// point on east side
		this->edges_ref(nb_edges_ref).Init(ne, ne+nbPoints_x, ref_boundary(1));
		nb_edges_ref++;    
	      }
	    if ((j==nbPoints_y-1)&&(i<nbPoints_x-1)) 
	      {
		// point on north side
		this->edges_ref(nb_edges_ref).Init(ne, ne+1, ref_boundary(2));
		nb_edges_ref++;    
	      }
	    ne++;
	  }
      }
    
    this->SortPeriodicBoundaries();
    this->SortBoundariesRef();
    
    this->ReorientElements();            
    this->FindConnectivity();
    this->ProjectPointsOnCurves();    
  }


  //! creation of a regular mesh of a rectangular domain [xmin xmax ymin ymax]
  //! with user defined subdivision and references
  template<class T>
  void Mesh<Dimension2, T>::
  CreateStructuredMesh(const Vector<T>& pos_layer_x, const Vector<T>& pos_layer_y,
                       Vector<int>& nb_layers_x, Vector<int>& nb_layers_y,
                       const Matrix<int>& ref_domain_layer, const TinyVector<int, 4>& ref_boundary)
  {
    // previous mesh is cleared
    Clear();

    int ref_max = ref_boundary(0);
    for (int i = 1; i < 4; i++)
      ref_max = max(ref_max, ref_boundary(i));
    
    if (ref_max > this->GetNbReferences())
      this->ResizeNbReferences(ref_max);

    int nbx = 1;
    for (int i = 0; i < nb_layers_x.GetM(); i++)
      nbx += nb_layers_x(i);
    
    int nby = 1;
    for (int i = 0; i < nb_layers_y.GetM(); i++)
      nby += nb_layers_y(i);

    int nb_vertices = nbx*nby;
    int nb_elt = 0;
    for (int i = 0; i < nb_layers_x.GetM(); i++)
      for (int j = 0; j < nb_layers_y.GetM(); j++)
        if (ref_domain_layer(i, j) > 0)
          nb_elt += nb_layers_x(i)*nb_layers_y(j);
    
    this->ReallocateVertices(nb_vertices);    
    this->ReallocateElements(nb_elt);
    this->ReallocateBoundariesRef(2*(nbx+nby-2));
    int ne = 0;
    nb_elt = 0; int nb_edges_ref = 0;

    Matrix<int> ref_domain(nbx-1, nby-1);    
    Vector<T> pos_x(nbx), pos_y(nby);
    pos_x(0) = pos_layer_x(0); int nb = 1;
    int offset_ix = 0;
    for (int i = 0; i < nb_layers_x.GetM(); i++)
      {
        T x0 = pos_layer_x(i), x1 = pos_layer_x(i+1);
        for (int j = 1; j < nb_layers_x(i); j++)
          {
            T lam = T(j) / nb_layers_x(i);
            pos_x(nb) = (T(1)-lam)*x0 + lam*x1;
            nb++;
          }
        
        pos_x(nb) = x1; nb++;
        
        int offset_iy = 0;
        for (int iy = 0; iy < nb_layers_y.GetM(); iy++)
          {
            for (int j = 0; j < nb_layers_x(i); j++)                
              for (int jy = 0; jy < nb_layers_y(iy); jy++)                  
                {
                  int ig = offset_ix + j;
                  int jg = offset_iy + jy;
                  ref_domain(ig, jg) = ref_domain_layer(i, iy);
                }
            offset_iy += nb_layers_y(iy);
          }
        
        offset_ix += nb_layers_x(i);
      }
    
    pos_y(0) = pos_layer_y(0); nb = 1;
    for (int i = 0; i < nb_layers_y.GetM(); i++)
      {
        T x0 = pos_layer_y(i), x1 = pos_layer_y(i+1);
        for (int j = 1; j < nb_layers_y(i); j++)
          {
            T lam = T(j) / nb_layers_y(i);
            pos_y(nb) = (T(1)-lam)*x0 + lam*x1;
            nb++;
          }
        
        pos_y(nb) = x1; nb++;
      }
    
    Matrix<int> num(2, 2);
    num(0, 0) = 0;
    num(1, 0) = 1;
    num(1, 1) = nbx+1;
    num(0, 1) = nbx;
    
    nb_elt = 0; nb_edges_ref = 0;
    for (int j = 0; j < nby; j++)
      for (int i = 0; i < nbx; i++)
        {            
          if ((i<(nbx-1))&&(j<(nby-1)))
            {
              if (ref_domain(i, j) > 0)
                this->elements(nb_elt++).
                  InitQuadrangular(num(0,0), num(1,0), num(1,1),
                                   num(0,1), ref_domain(i, j));
              
              // boundary x = xmin
              if (i==0)
                {
                  this->BoundaryRef(nb_edges_ref++).
                      Init(num(0,0), num(0,1), ref_boundary(3));
                }
              
              // boundary x = xmax
              if (i==(nbx-2))
                {
                  this->BoundaryRef(nb_edges_ref++).
                    Init(num(1,0), num(1,1), ref_boundary(1));
                }
              
              // boundary y = ymin
              if (j==0)
                this->BoundaryRef(nb_edges_ref++).Init(num(0,0),num(1,0), ref_boundary(0));
              
              // boundary y = ymax
              if (j==(nby-2))
                this->BoundaryRef(nb_edges_ref++)
                  .Init(num(0,1), num(1,1), ref_boundary(2));
            }
	    
          this->Vertices(num(0,0)).Init(pos_x(i), pos_y(j));
          
          for (int l = 0; l < 2; l++)
            for (int m = 0; m < 2; m++)
                num(l, m)++;
        }
    
    this->SortPeriodicBoundaries();
    this->SortBoundariesRef();
    
    this->ReorientElements();            
    this->FindConnectivity();
    this->ProjectPointsOnCurves();    
  }


  //! extrusion of the mesh along a coordinate (x, y or z)
  /*!
    \param[in] num_coord coordinate
    \param[in] nb_layes number of cells to add in the direction of the extrusion
    \param[in] pos position of the plane x = pos (for the first coordinate)
    \param[in] delta height of extrusion
    this function is used to add PML layers (extern boundary is a rectangle)
  */  
  template<class T>
  void Mesh<Dimension2, T>::ExtrudeCoordinate(int num_coord, int nb_layers,
					      const T& pos, const T& delta)
  {    
    // points on the boundary are found
    IVect Point_to_extrude(this->GetNbVertices()); Point_to_extrude.Fill(-1);
    
    int nb_new_vertices = 0;
    for (int i = 0; i < this->GetNbVertices(); i++)
      if (abs(pos-this->Vertices(i)(num_coord)) <= R_N::threshold)
	{
	  Point_to_extrude(i) = this->GetNbVertices() + nb_new_vertices;
	  nb_new_vertices++;
	}

    bool quad_layer = true;
    
    // points are sorted
    IVect IndexPoints(nb_new_vertices);
    VectR_N AllPoints(nb_new_vertices); nb_new_vertices = 0;
    for (int i = 0; i < this->GetNbVertices(); i++)
      if (Point_to_extrude(i) >= 0)
	{
	  IndexPoints(nb_new_vertices) = i;
	  AllPoints(nb_new_vertices) = this->Vertices(i);
	  nb_new_vertices++;
	}
    
    Sort(nb_new_vertices, AllPoints, IndexPoints);
    
    // intermediary steps between xmax and xmax + delta (if num_coor = 0)
    VectReal_wp step_subdiv(nb_layers);
    for (int i = 0; i < nb_layers; i++)
      step_subdiv(i) = pos + Real_wp(i+1)*delta/nb_layers;
    
    // creation of vertices
    this->ResizeVertices(this->GetNbVertices() + nb_new_vertices*nb_layers);
    R_N ptM;
    for (int k = 0; k < nb_new_vertices; k++)
      {
	ptM = this->Vertex(IndexPoints(k));
	for (int m = 0; m < nb_layers; m++)
	  {
	    ptM(num_coord) = step_subdiv(m);
	    this->Vertex(Point_to_extrude(IndexPoints(k)) + m*nb_new_vertices) = ptM;
	  }
      }

    int nb_edges_ref = this->GetNbBoundaryRef();
    VectBool edge_on_top(nb_edges_ref);
    int nb_new_edges_ref = 0;
    // loop over all edges of the boundary
    // if edge_on_top(i), edge i has to be extruded
    for (int i = 0; i < nb_edges_ref; i++)
      {
	edge_on_top(i) = true;
	for (int j = 0; j < 2; j++)
	  {
	    int n0 = this->edges_ref(i).numVertex(j);
	    if (Point_to_extrude(n0) == -1)
	      edge_on_top(i) = false;
	  }
      }
    
    // we count the number of elements for each vertex
    int nb_new_elt = 0;
    Matrix<int> Edge_Point(this->GetNbVertices(), 2); Edge_Point.Fill(-1);
    for (int i = 0; i < nb_edges_ref; i++)
      if (edge_on_top(i))
	{
	  for (int j = 0; j < 2; j++)
	    {
	      int n0 = this->BoundaryRef(i).numVertex(j);
	      if (Edge_Point(n0,0) == -1)
		Edge_Point(n0, 0) = i;
	      else
		Edge_Point(n0, 1) = i;
	    }
	   
	  if (quad_layer)
	    nb_new_elt += 1;
	  else
	    nb_new_elt += 2;
	}
    
    // for edges on the extremity, we have to insert new lateral edges
    nb_new_edges_ref = 0;
    for (int i = 0; i < this->GetNbVertices(); i++)
      if ((Edge_Point(i, 0) != -1)&&(Edge_Point(i, 1) == -1))
	nb_new_edges_ref++;
    
    int nb_elt = this->GetNbElt();
    // changing the size of arrays
    this->ResizeBoundariesRef(nb_edges_ref + nb_new_edges_ref*nb_layers);
    this->ResizeElements(nb_elt + nb_new_elt*nb_layers);

    IVect num(4), numb(4); VectR2 ptQ(4);
    nb_new_elt = nb_elt;
    // loop over each edge to extrude
    for (int i = 0; i < nb_edges_ref; i++)
      if (edge_on_top(i))
	{
	  int ref = this->edges_ref(i).GetReference();
	  int type = this->elements(this->edges_ref(i).numElement(0)).GetReference();
	  num(0) = this->edges_ref(i).numVertex(0);
	  num(1) = this->edges_ref(i).numVertex(1);
	  num(2) = Point_to_extrude(num(1));
	  num(3) = Point_to_extrude(num(0));
	  if (quad_layer)
	    {
	      // we create quadrilaterals
	      for (int m = 0; m < nb_layers; m++)
		{
                  numb(0) = num(0); numb(1) = num(1); numb(2) = num(2); numb(3) = num(3);
		  ptQ(0) = this->Vertex(num(0));
                  ptQ(1) = this->Vertex(num(1));
                  ptQ(2) = this->Vertex(num(2));
                  ptQ(3) = this->Vertex(num(3));
                  Sort(ptQ, numb);
                  this->Element(nb_new_elt).InitQuadrangular(numb(0), numb(2), numb(3), numb(1), type);
		  num(0) = num(3); num(1) = num(2);
		  num(2) += nb_new_vertices; num(3) += nb_new_vertices;
                  nb_new_elt++;
		}	      
	    }
	  else
	    {
	      // we create triangles
	      for (int m = 0; m < nb_layers; m++)
		{
		  this->Element(nb_new_elt++).InitTriangular(num(0), num(1), num(2), type);
		  this->Element(nb_new_elt++).InitTriangular(num(0), num(2), num(3), type);
		  num(0) = num(3); num(1) = num(2);
		  num(2) += nb_new_vertices; num(3) += nb_new_vertices;
		}
	    }
	  
	  // on change l'arete
	  this->BoundaryRef(i).Init(num(0), num(1), ref);
	}
    
    // for extremities, we add referenced edges
    for (int i = 0; i < this->GetNbVertices(); i++)
      if ((Edge_Point(i,0) != -1)&&(Edge_Point(i,1) == -1))
	{
	  int ref = this->BoundaryRef(Edge_Point(i,0)).GetReference();
	  // a-t-on d'autres references possibles ?
	  int other_ref;
	  this->GetPointWithOtherReference(i, ref, ptM, other_ref);
	  
	  // we check that the neighbor point ptM is in the good direction
	  R_N vec_u = this->Vertex(i) - ptM;
	  if ( ((num_coord == 0) && (abs(vec_u(1)) <= 10*R_N::threshold))
	       || ((num_coord == 1) && (abs(vec_u(0)) <= 10*R_N::threshold)))
	    {
	      // if it is the case, we take this reference
	      if (other_ref != 0)
		ref = other_ref;
	    }

	  int n1 = i; int n2 = Point_to_extrude(i);
	  for (int m = 0; m < nb_layers; m++)
	    {
	      this->BoundaryRef(nb_edges_ref + m*nb_new_edges_ref).Init(n1, n2, ref);
	      n1 = n2;
	      n2 += nb_new_vertices;
	    }
	  nb_edges_ref++;
	}    
    
    // we reconstruct relation between referenced edges and edges
    this->ReorientElements();
    this->FindConnectivity();
    this->ProjectPointsOnCurves();    
  }
  
  
  //! general construction of a mesh
  /*!
    \param[in] type_mesh (hybrid, purely quadrilateral or purely triangular)
    \param[in] vars object containing problem datas
    \param[in] parameters matching line of the data file
    this method can read or generate mesh, depending what the user wants
    the mesh is then split, depending the type. Boundaries are computed,
    and if asked, PML layers are added
  */
  template<class T>
  void Mesh<Dimension2, T>::ConstructMesh(int type_mesh, const VectString& parameters)
  {
    if (this->required_mesh_type >= 0)
      type_mesh = this->required_mesh_type;
    
    int factor_refinement_mesh = 1;
    
    // predefined mesh or file ?
    if (!parameters(0).compare("REGULAR"))
      {
	int ref_domain = 1;
	TinyVector<int, 4> ref_boundary; ref_boundary.Fill(3);

	// regular mesh of a rectangle
	TinyVector<int, 2> nbPoints;
	int nbPoints_x = to_num<int>(parameters(1));
	R_N ptMin, ptMax; ptMax.Fill(1);
	if (parameters.GetM() >=6)
	  {
	    to_num(parameters(2), ptMin(0)); to_num(parameters(3), ptMax(0));
	    to_num(parameters(4), ptMin(1)); to_num(parameters(5), ptMax(1));
	  }

	int nbPoints_y = toInteger(ceil((ptMax(1) - ptMin(1))/
					(ptMax(0) - ptMin(0))*nbPoints_x));
	
	nbPoints(0) = nbPoints_x; nbPoints(1) = nbPoints_y;
	
	if (parameters.GetM() >= 8)
	  {
	    to_num(parameters(6), ref_domain);
	    to_num(parameters(7), ref_boundary(0));
	    if (parameters.GetM() >= 11)
	      {
		to_num(parameters(8), ref_boundary(1));
		to_num(parameters(9), ref_boundary(2));
		to_num(parameters(10), ref_boundary(3));
	      }
	    else
	      ref_boundary.Fill(ref_boundary(0));
	  }

	// creation of a regular mesh (tri/tet or quad/hex depending on finite element)
        if (this->print_level >= 4)
          cout<<"Creation of a regular mesh "<<endl;
	
	CreateRegularMesh(ptMin, ptMax, nbPoints,
			  ref_domain, ref_boundary, type_mesh);
      }
    else if (!parameters(0).compare("REGULAR_ANISO"))
      {
	int ref_domain = 1;
	TinyVector<int, 4> ref_boundary; ref_boundary.Fill(3);

	// regular mesh of a rectangle
	TinyVector<int, 2> nbPoints;
	int nbPoints_x = to_num<int>(parameters(1));
	int nbPoints_y = to_num<int>(parameters(2));
	R_N ptMin, ptMax; ptMax.Fill(1);
	if (parameters.GetM() >= 7)
	  {
	    to_num(parameters(3), ptMin(0)); to_num(parameters(4), ptMax(0));
	    to_num(parameters(5), ptMin(1)); to_num(parameters(6), ptMax(1));
	  }

	nbPoints(0) = nbPoints_x; nbPoints(1) = nbPoints_y;
	
	if (parameters.GetM() >= 9)
	  {
	    to_num(parameters(7), ref_domain);
	    to_num(parameters(8), ref_boundary(0));
	    if (parameters.GetM() >= 12)
	      {
		to_num(parameters(9), ref_boundary(1));
		to_num(parameters(10), ref_boundary(2));
		to_num(parameters(11), ref_boundary(3));
	      }
	    else
	      ref_boundary.Fill(ref_boundary(0));
	  }

	// creation of a regular mesh (tri/tet or quad/hex depending on finite element)
        if (this->print_level >= 4)
          cout<<"Creation of a regular mesh "<<endl;
	
	CreateRegularMesh(ptMin, ptMax, nbPoints,
			  ref_domain, ref_boundary, type_mesh);
      }
    else
      {
	// file_mesh = name of the file where the mesh is stored
	string file_mesh = string(this->mesh_path);
	file_mesh += parameters(0);
	// otherwise the mesh is read
	Read(file_mesh);
	for (int i = 1; i < (parameters.GetM()-1); i++)
	  {
	    if (!parameters(i).compare("REFINED"))
	      to_num(parameters(i+1), factor_refinement_mesh);
	  }
      }
    
    if (factor_refinement_mesh < 1)
      factor_refinement_mesh = 1;
    
    // adding asked vertices
    for (int i = 0; i < this->ListeVertices_to_be_added.GetM(); i++)
      AddVertexConformal(this->ListeVertices_to_be_added(i));

    // local refinement around a set of curves
    if (this->ParamRefinement_Area.GetM() > 0)
      {
	const Vector<string>& param = this->ParamRefinement_Area;
	int level = to_num<int>(param(0));
	T ratio = to_num<T>(param(1));
	int nb = 2;
        if (param(nb) == "Ref")
          {
            nb++;
            bool anisotrope = false;
            if (param(nb) == "Anisotrope")
              anisotrope = true;
            
            nb++;
            IVect ref_cond(this->GetNbReferences()+1);
            ref_cond.Zero();
            for (int j = nb; j < param.GetM(); j++)
              {
                int ref = to_num<int>(param(j));
                if (ref < ref_cond.GetM())
                  ref_cond(ref) = 1;
              }
            
            this->ApplyRefinementEdge(ref_cond, 1, anisotrope, level, ratio);
          }
      }
    
    int nb_vertices_to_be_refined = this->ListeVertices_to_be_refined.GetM();
    if (this->print_level >= 4)
      cout<<"Number of vertices to be refined "<<nb_vertices_to_be_refined<<endl;
    
    // the user can ask a local refinement around a vertex
    for (int i = 0; i < nb_vertices_to_be_refined; i++)
      {
	int num_vertex = this->FindVertexNumber(this->ListeVertices_to_be_refined(i));
	if (num_vertex != -1)
	  ApplyLocalRefinement(num_vertex, this->LevelRefinement_Vertex(i),
			       this->RatioRefinement_Vertex(i));
	else if (this->print_level >= 0)
	  {
	    cout<<"we didn't find the vertex " << this->ListeVertices_to_be_refined(i)
		<<" in the mesh"<<endl;
	    abort();
	  }
      }
    
    // Mesh is split if necessary
    // for example if finite element method accepts only quadrilaterals
    // and the mesh comprises triangles and quadrilaterals
    // then each triangle is split in three and quadrilateral in four
    if ((type_mesh == this->TRIANGULAR_MESH)&&(!IsOnlyTriangular()))
      SplitIntoTriangles();
    else if ((type_mesh == this->QUADRILATERAL_MESH)&&(!IsOnlyQuadrilateral()))
      SplitIntoQuadrilaterals();
    else
      {
	// if the user wants very distorted mesh, irregular_mesh_double_cut = true
	// then the mesh is split twice
	if ((this->irregular_mesh_double_cut)
	    &&(type_mesh == this->TRIANGULAR_MESH))
	  {
	    SplitIntoQuadrilaterals(); 
	    SplitIntoTriangles(); 
	  }
	else if ((this->irregular_mesh_double_cut)&&
		 (type_mesh == this->QUADRILATERAL_MESH))
	  {
	    SplitIntoTriangles();
	    SplitIntoQuadrilaterals();
	  }
      }
    
    // if the user wants a finer mesh, it can be split in two or more intervals
    if (factor_refinement_mesh != 1)
      {
	// each edge is split in n intervals where n = factor_mesh_refinement
	// and consequently elements are split...
	Vector<T> step_subdiv(factor_refinement_mesh+1);
	for (int i = 0; i <= factor_refinement_mesh; i++)
	  step_subdiv(i) = T(i)/factor_refinement_mesh;
	
	SubdivideMesh(step_subdiv);
      }
    
    // If the user wants extra-layers (PML layers) around the computational
    // they are added to the existing mesh
    for (int i = 0; i < this->pml_areas.GetM(); i++)
      this->pml_areas(i).AddPML(i, *this);
    
    // all curves are transformed to "file" curves
    // such that parallel computation will be based on the same geometrical mesh
    for (int ref = 0; ref < this->Type_curve.GetM(); ref++)
      if (this->Type_curve(ref) > 0)
        this->Type_curve(ref) = this->CURVE_FILE;
    
    // triangles are placed first and then quadrilaterals
    this->RearrangeElements();
  }

  
  //! triangles are placed first, then quadrilaterals
  template<class T>
  void Mesh<Dimension2, T>::RearrangeElements()
  {
    IVect permut(this->GetNbElt());
    int nb = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->Element(i).GetNbVertices() == 3)
        permut(nb++) = i;
    
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->Element(i).GetNbVertices() == 4)
        permut(nb++) = i;
    
    bool permut_elt = false;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (permut(i) != i)
        permut_elt = true;
    
    if (permut_elt)
      {
        Vector<Face<Dimension2> > OldElement = this->elements;
        Vector<bool> OldAffine = this->affine_element;
        IVect inv_permut(this->GetNbElt());
        for (int i = 0; i < this->GetNbElt(); i++)
          inv_permut(permut(i)) = i;
        
        for (int i = 0; i < this->GetNbElt(); i++)
          {
            this->Element(i) = OldElement(permut(i));
            this->affine_element(i) = OldAffine(permut(i));
          }
        
        for (int i = 0; i < this->edges.GetM(); i++)
          {
            int ne0 = this->edges(i).numElement(0);
            int ne1 = this->edges(i).numElement(1);
            if (ne0 >= 0)
              ne0 = inv_permut(ne0);

            if (ne1 >= 0)
              ne1 = inv_permut(ne1);
            
            this->edges(i).SetElement(ne0, ne1);
          }

        for (int i = 0; i < this->edges_ref.GetM(); i++)
          {
            int ne0 = this->edges_ref(i).numElement(0);
            int ne1 = this->edges_ref(i).numElement(1);
            if (ne0 >= 0)
              ne0 = inv_permut(ne0);

            if (ne1 >= 0)
              ne1 = inv_permut(ne1);
            
            this->edges_ref(i).SetElement(ne0, ne1);
          }
      }

    UpdateReferenceElement();
  }
  
  
  //! creation of a sub-mesh with vertices and elements given by the user
  /*!
    \param[out] sub_mesh resulting sub-mesh
    \param[in] nb_vertices_subdomain vertex numbers to be taken
    \param[in] nb_elt_subdomain element numbers to be taken
    \param[in] VertexOnSubdomain if true, the vertex i is on the sub-mesh
    \param[in] ElementOnSubdomain if true, the element i is on the sub-mesh
   */
  template<class T>
  void Mesh<Dimension2, T>::
  CreateSubmesh(Mesh<Dimension2, T>& sub_mesh, int nb_vertices_subdomain,
		int nb_elt_subdomain, const VectBool& VertexOnSubdomain,
		const VectBool& ElementOnSubdomain) const 
  {
    // sub mesh is cleared
    sub_mesh.Clear();
    
    // updating mesh
    IVect NewNumberVertex(this->GetNbVertices()); NewNumberVertex.Fill(-1);
    // vertices
    sub_mesh.ReallocateVertices(nb_vertices_subdomain);
    nb_vertices_subdomain = 0;
    for (int n = 0; n < this->GetNbVertices(); n++)
      if (VertexOnSubdomain(n))
	{
	  NewNumberVertex(n) = nb_vertices_subdomain;
	  sub_mesh.Vertex(nb_vertices_subdomain++) = this->Vertex(n);
	}
    
    IVect NewNumberElement(this->GetNbElt()); NewNumberElement.Fill(-1);
    // elements
    sub_mesh.ReallocateElements(nb_elt_subdomain);
    nb_elt_subdomain = 0;
    for (int e = 0; e < this->GetNbElt(); e++)
      if (ElementOnSubdomain(e))
	{
          int nb_vert = this->Element(e).GetNbVertices();
          IVect num(nb_vert);
	  for (int n = 0; n < nb_vert; n++)
	    num(n) = NewNumberVertex(this->Element(e).numVertex(n));
          
          sub_mesh.Element(nb_elt_subdomain).Init(num, this->Element(e).GetReference());
	  sub_mesh.Element(nb_elt_subdomain)
            .SetPML(this->Element(e).GetNumberPML(), this->Element(e).GetTypePML());
          
	  nb_elt_subdomain++;
	}
    
    sub_mesh.FindConnectivity();
    
    // edges
    int nb_edges_subdomain = 0;
    for (int e = 0; e < this->GetNbBoundaryRef(); e++)
      {
        int n1 = NewNumberVertex(this->BoundaryRef(e).numVertex(0));
        int n2 = NewNumberVertex(this->BoundaryRef(e).numVertex(1));
        if ((n1 >= 0) && (n2 >= 0))
          {
            int ne = sub_mesh.FindEdgeWithExtremities(n1, n2);
            if (ne >= 0)
              nb_edges_subdomain++;
          }
      }
    
    bool only_quad = sub_mesh.IsOnlyQuadrilateral();
    sub_mesh.SetGeometryOrder(this->GetGeometryOrder(), only_quad);
    sub_mesh.ReallocateBoundariesRef(nb_edges_subdomain);
    nb_edges_subdomain = 0;
    for (int e = 0; e < this->GetNbBoundaryRef(); e++)
      {
        int n1 = NewNumberVertex(this->BoundaryRef(e).numVertex(0));
        int n2 = NewNumberVertex(this->BoundaryRef(e).numVertex(1));
        if ((n1 >= 0) && (n2 >= 0))
          {
            int ne = sub_mesh.FindEdgeWithExtremities(n1, n2);
            if (ne >= 0)
              {
                sub_mesh.BoundaryRef(nb_edges_subdomain).
                  Init(n1, n2, this->BoundaryRef(e).GetReference());
                
                for (int k = 0; k < this->GetGeometryOrder()-1; k++)
                  sub_mesh.PointsEdgeRef(nb_edges_subdomain, k) = this->PointsEdgeRef(e, k);
                
                nb_edges_subdomain++;
              }
          }
      }
    
    sub_mesh.Type_curve = this->Type_curve; sub_mesh.Cond_curve = this->Cond_curve;
    sub_mesh.Body_curve = this->Body_curve; sub_mesh.Param_curve = this->Param_curve;
    sub_mesh.Param_known = this->Param_known;
    sub_mesh.reference_original_neighbor = this->reference_original_neighbor;
    
    sub_mesh.ReorientElements();
    sub_mesh.FindConnectivity();
    sub_mesh.ProjectPointsOnCurves();

  }

  
  //! a mesh is appended to the current mesh
  /*!
    \param[in] mesh the mesh to add
    \param[in] elimine if true, this method eliminates vertices in double
  */
  template<class T>
  void Mesh<Dimension2, T>::AppendMesh(Mesh<Dimension2, T>& mesh, bool elimine,
                                       R2 vec_u)
  {
    if (mesh.GetNbReferences() > this->GetNbReferences())
      this->ResizeNbReferences(mesh.GetNbReferences()+1);
    
    Matrix<Real_wp> ValPhi;
    // interpolation necessary if the order of two meshes are different
    if (this->GetGeometryOrder() != mesh.GetGeometryOrder())
      {
	ValPhi.Reallocate(mesh.GetGeometryOrder()+1, this->GetGeometryOrder()-1);
	for (int i = 0; i < ValPhi.GetM(); i++)
	  for (int j = 0; j < ValPhi.GetN(); j++)
	    ValPhi(i, j) = mesh.lob_curve.EvaluatePhi(i, this->lob_curve.Points(j+1));
	
      }
    
    // rajout des sommets
    int nb_old = this->GetNbVertices(), offset = nb_old;
    this->ResizeVertices(nb_old + mesh.GetNbVertices());
    for (int i = nb_old; i < this->GetNbVertices(); i++)
      this->Vertices(i) = mesh.Vertices(i-nb_old) + vec_u;
    
    // rajouts des faces de reference
    nb_old = this->GetNbBoundaryRef();
    this->ResizeBoundariesRef(nb_old + mesh.GetNbBoundaryRef());
    for (int i = nb_old; i < this->GetNbBoundaryRef(); i++)
      {
	this->BoundaryRef(i) = mesh.BoundaryRef(i-nb_old);
	// decalage des sommets
	this->BoundaryRef(i).Init(this->BoundaryRef(i).numVertex(0)+offset,
				  this->BoundaryRef(i).numVertex(1)+offset);
	
	if (this->GetGeometryOrder() != mesh.GetGeometryOrder())
	  {
	    // interpolation
	    VectR2 OldPoints(mesh.GetGeometryOrder()+1);
	    OldPoints(0) = mesh.Vertex(mesh.BoundaryRef(i-nb_old).numVertex(0)) + vec_u;
	    OldPoints(mesh.GetGeometryOrder())
              = mesh.Vertex(mesh.BoundaryRef(i-nb_old).numVertex(1)) + vec_u;
	    for (int k = 1; k < mesh.GetGeometryOrder(); k++)
	      OldPoints(k) = mesh.PointsEdgeRef(i-nb_old, k-1) + vec_u;
	    
	    for (int k = 0; k < this->GetGeometryOrder()-1; k++)
	      {
		this->PointsEdgeRef(i, k).Zero();
		for (int j = 0; j < OldPoints.GetM(); j++)
		  Add(ValPhi(j, k), OldPoints(j), this->PointsEdgeRef(i, k));
	      }
	  }
	else
	  for (int k = 0; k < this->GetGeometryOrder()-1; k++)
	    this->PointsEdgeRef(i, k) = mesh.PointsEdgeRef(i-nb_old, k) + vec_u;
      }
    
    // rajouts des elements volumiques
    nb_old = this->GetNbElt();
    this->ResizeElements(nb_old + mesh.GetNbElt());
    for (int i = nb_old; i < this->GetNbElt(); i++)
      {
        int nb_vert = mesh.Element(i-nb_old).GetNbVertices();
	IVect num(nb_vert);
	// decalage des sommets
	for (int j = 0; j < nb_vert; j++)
	  num(j) = mesh.Element(i-nb_old).numVertex(j)+offset;
        
        this->Element(i).Init(num, mesh.Element(i-nb_old).GetReference());
      }
    
    this->ReorientElements();
    // on elimine tous les sommets doublons
    if (elimine)
      this->RemoveDuplicateVertices();
    else
      this->FindConnectivity();
    
    this->ProjectPointsOnCurves();
  }
  
  
  /********************
   * Internal methods *
   ********************/
  
  
  //! finds minimum and maximum angle of mesh
  template<class T>
  void Mesh<Dimension2, T>::GetTetaMinMax(const R_N& center,
					  T& teta_min, T& teta_max) const
  {
    teta_min = T(10); teta_max = T(-10);
    // computation of teta_min, teta_max
    T teta, r; R_N vec_u;
    for (int i = 0; i < this->GetNbVertices(); i++)
      {
	// polar coordinates (r, teta)
	vec_u = this->Vertices(i) - center;
	r = this->Vertices(i).Distance(center);
	if (r > R_N::threshold)
	  {
	    if (vec_u(0) >= r)
	      teta = 0;
	    else if (vec_u(0) <= -r)
	      teta = pi_wp;
	    else
	      {
		teta = acos(vec_u(0)/r);
		if (vec_u(1) < 0)
		  teta = -teta;
	      }
	    
	    teta_min = min(teta_min, teta);
	    teta_max = max(teta_max, teta);
	  }
      }
    
    if ((teta_max-teta_min) > pi_wp)
      {
	// inversion of teta_max and teta_min + 2pi
	teta = teta_max;
	teta_max = teta_min+2.0*pi_wp;
	teta_min = teta;
      }
  }
  
  
  //! finds the edges and vertices on the Y-axis (x = 0)
  /*!
    \param[out] list_vertices vertices on the axis x = 0
    \param[out] list_edges edges on the axis
    \param[out] Vertex_On_Axe is a vertex on the axis x = 0 ?
   */
  template<class T>
  void Mesh<Dimension2, T>::
  GetEdgesOnZaxis(IVect& list_vertices, IVect& list_edges, VectBool& Vertex_On_Axe)
  {
    int nb_vertices_axis = 0;
    list_vertices.Reallocate(this->GetNbVertices());
    Vertex_On_Axe.Reallocate(this->GetNbVertices());
    Vertex_On_Axe.Fill(false);
    for (int i = 0; i < this->GetNbVertices(); i++)
      if (abs(this->Vertex(i)(0)) <= R_N::threshold)
	{
	  list_vertices(nb_vertices_axis) = i; 
	  Vertex_On_Axe(i) = true;
	  nb_vertices_axis++;
	}
    list_vertices.Resize(nb_vertices_axis);
    
    int nb_edges_axis = 0;
    list_edges.Reallocate(2*this->GetNbVertices());
    for (int i = 0; i < this->GetNbEdges(); i++)
      {
	int n1 = this->edges(i).numVertex(0), n2 = this->edges(i).numVertex(1);
	if (Vertex_On_Axe(n1)&&Vertex_On_Axe(n2))
	  list_edges(nb_edges_axis++) = i;
      }
    list_edges.Resize(nb_edges_axis);
  }
  
  
  //! adds PML layers along a coordinate of the mesh
  /*!
    if circle_layer is false, adds elements in direction of axis
    if circle_layer is true, adds circular layers, axis is the center of the circle
   */
  template<class T>
  void ExtrudePMLLayer(Mesh<Dimension2, T>& mesh, const TinyVector<T, 2>& axis,
		       const T& scal_max, const Vector<int>& ref_cond, int ref_target,
		       const T& thickness, int nb_layers, bool circle_layer)
  {            
    //mesh.Write("initial.mesh");
    
    // now we select all the vertices on this same "boundary line"
    int N = mesh.GetNbVertices();
    IVect ListeVertices(N), IndexVert(N);
    int nb_vert = 0; IndexVert.Fill(-1); T prod_scal;
    if (ref_target == -1)
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	{
	  if (circle_layer)
	    prod_scal = mesh.Vertex(i).Distance(axis);
	  else
	    prod_scal = DotProd(mesh.Vertex(i), axis);
	  
	  if (abs(prod_scal - scal_max) <= TinyVector<T, 2>::threshold)
	    {
	      ListeVertices(nb_vert) = i;
	      IndexVert(i) = nb_vert;
	      nb_vert++;
	    }
	}
    else
      {
	for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	  if (ref_cond(mesh.BoundaryRef(i).GetReference()) == ref_target)
	    {
	      int n1 = mesh.BoundaryRef(i).numVertex(0);
	      int n2 = mesh.BoundaryRef(i).numVertex(1);
	      if (IndexVert(n1) == -1)
		{
		  IndexVert(n1) = 1;
		  nb_vert++;
		}
	      
	      if (IndexVert(n2) == -1)
		{
		  IndexVert(n2) = 1;
		  nb_vert++;
		}
	    }
	
	nb_vert = 0;
	for (int i = 0; i < mesh.GetNbVertices(); i++)
	  if (IndexVert(i) > 0)
	    {
	      ListeVertices(nb_vert) = i;
	      IndexVert(i) = nb_vert;
	      nb_vert++;
	    }
      }
    
    ListeVertices.Resize(nb_vert);

    // now we detect vertices, that have already been extruded
    // for each vertex, we store the numbers of the translated points
    Matrix<int> NumTranslat(nb_layers, nb_vert);
    NumTranslat.Fill(-1);
    int nb_new_vertices = 0;
    Vector<bool> OriginalVertex(nb_vert);
    OriginalVertex.Fill(true);
    for (int i = 0; i < nb_vert; i++)
      {
	int nv = ListeVertices(i);
	Vector<int> num_elem;
	mesh.GetNeighboringElementsAroundVertex(nv, num_elem);
	int pml_elt = -1;
	for (int k = 0; k < num_elem.GetM(); k++)
	  if (mesh.Element(num_elem(k)).IsPML())
	    {
	      pml_elt = num_elem(k);
	      break;
	    }
	
	if (pml_elt == -1)
	  {
	    // vertex that will create newly extruded vertices
	    nb_new_vertices++;
	  }
	else
	  {
	    OriginalVertex(i) = false;
	    // we take numbers of the adjacent PML
	    int n_current = nv;
	    for (int k = 0; k < nb_layers; k++)
	      {
		int nloc = -1;
		int nb_v = mesh.Element(pml_elt).GetNbVertices();
		for (int j = 0; j < nb_v; j++)
		  if (mesh.Element(pml_elt).numVertex(j) == n_current)
		    nloc = j;
		
		int num_edge = mesh.Element(pml_elt).numBoundary(nloc);
		int n_next = -1;
		if (mesh.Boundary(num_edge).GetNbElements() == 1)
		  n_next = mesh.Element(pml_elt).numVertex((nloc+1)%nb_v);
		else
		  {
		    num_edge = mesh.Element(pml_elt).numBoundary((nloc+nb_v-1)%nb_v);
		    n_next = mesh.Element(pml_elt).numVertex((nloc+nb_v-1)%nb_v);
		  }
		
		// the reference is set to 0 for this shared boundary
		mesh.Boundary(num_edge).SetReference(0);
		mesh.BoundaryRef(num_edge).SetReference(0);
		
		NumTranslat(k, i) = n_next;
		n_current = n_next;
		
		// now we search the next element
		for (int j = 0; j < nb_v; j++)
		  {
		    num_edge = mesh.Element(pml_elt).numBoundary(j);
		    if (mesh.Boundary(num_edge).GetNbElements() == 2)
		      {
			int ne = mesh.Boundary(num_edge).numElement(0);
			if (ne == pml_elt)
			  ne = mesh.Boundary(num_edge).numElement(1);
			
			bool pt_found = false;
			for (int j2 = 0; j2 < mesh.Element(ne).GetNbVertices(); j2++)
			  if (mesh.Element(ne).numVertex(j2) == n_next)
			    pt_found = true;
			
			if (pt_found)
			  {
			    pml_elt = ne;
			    break;
			  }
		      }
		  }		
	      }
	  }
      }
    
    // creating the new vertices
    mesh.ResizeVertices(N + nb_new_vertices*nb_layers);
    int offset = N;
    for (int k = 0; k < nb_layers; k++)
      {
	TinyVector<T, 2> translat = axis;
	Mlt(T(k+1)*thickness/(nb_layers*Norm2(translat)), translat);
	for (int i = 0; i < nb_vert; i++)
	  if (OriginalVertex(i))
	    {
	      int num_vertex = ListeVertices(i);
	      mesh.Vertex(offset) = mesh.Vertex(num_vertex);
	      if (circle_layer)
		{
		  translat = mesh.Vertex(num_vertex) - axis;
		  Mlt(thickness*T(k+1)/(nb_layers*Norm2(translat)), translat);
		  
		  mesh.Vertex(offset)(0) += translat(0);
		  mesh.Vertex(offset)(1) += translat(1);
		}
	      else
		{
                  mesh.Vertex(offset)(0) += translat(0);
                  mesh.Vertex(offset)(1) += translat(1);
		}
              
	      NumTranslat(k, i) = offset;
	      offset++;
	    }
      }
    
    // couting the number of elements to add
    int Nbound = mesh.GetNbBoundaryRef();
    int nb_new_elt = 0;
    Vector<int> RefUsed(mesh.GetNbReferences()+1);
    RefUsed.Fill(0);
    for (int i = 0; i < Nbound; i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
        int ref = mesh.BoundaryRef(i).GetReference();
	if ((IndexVert(n1) >= 0)&&(IndexVert(n2) >= 0))
	  {
            nb_new_elt++;
            RefUsed(ref) = 1;
          }
      }
    
    // changing parameters of the circle
    if (circle_layer)
      for (int ref = 0; ref <= mesh.GetNbReferences(); ref++)
        if (RefUsed(ref) > 0)
          {
            Vector<T> param = mesh.GetCurveParameter(ref);            
            if (mesh.GetCurveType(ref) == mesh.CURVE_CIRCLE)
              {
                int new_ref = mesh.GetNewReference();
		if (new_ref >= RefUsed.GetM())
		  {
		    RefUsed.Resize(new_ref+1);
		    RefUsed(new_ref) = 0;
		  }
		
                RefUsed(ref) = new_ref;
                mesh.SetCurveType(new_ref, mesh.CURVE_CIRCLE);
                mesh.SetCurveParameter(new_ref, param);

                param(2) += thickness;
                mesh.SetCurveParameter(ref, param);
              }
          }
    
    // creating new elements
    int Nelt = mesh.GetNbElt();
    mesh.ResizeElements(Nelt + nb_new_elt*nb_layers);
    nb_new_elt = 0;
    for (int i = 0; i < Nbound; i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	int num_elem = mesh.BoundaryRef(i).numElement(0);
	int ref_inside = mesh.Element(num_elem).GetReference();
	if ((IndexVert(n1) >= 0)&&(IndexVert(n2) >= 0))
	  {
	    int n1_ = n1, n2_ = n2;
	    int n1_b = NumTranslat(0, IndexVert(n1));	    
	    int n2_b = NumTranslat(0, IndexVert(n2));
	    
	    for (int k = 0; k < nb_layers; k++)
	      {
		mesh.Element(Nelt + nb_new_elt + k).InitQuadrangular(n1_, n2_, n2_b, n1_b, ref_inside);
		n1_ = n1_b; n2_ = n2_b;
		if (k < nb_layers-1)
		  {
		    n1_b = NumTranslat(k+1, IndexVert(n1));
		    n2_b = NumTranslat(k+1, IndexVert(n2));
		  }
	      }
            
	    nb_new_elt += nb_layers;
	  }
      }
    
    int nb_new_bound = 0;
    for (int i = 0; i < Nbound; i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	if ((IndexVert(n1) >= 0) && (IndexVert(n2) >= 0))
          {
            if (circle_layer)
              nb_new_bound++;
          }
        else if (IndexVert(n1) >= 0)
          {
	    if (OriginalVertex(IndexVert(n1)))
	      nb_new_bound += nb_layers;
	  }
        else if (IndexVert(n2) >= 0)
          {
	    if (OriginalVertex(IndexVert(n2)))
	      nb_new_bound += nb_layers;
	  }
      }

    // creating the new edges
    mesh.ResizeBoundariesRef(Nbound + nb_new_bound);
    int Nbound_old = Nbound;
    for (int i = 0; i < Nbound_old; i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	int ref = mesh.BoundaryRef(i).GetReference();
	if ((IndexVert(n1) >= 0)&&(IndexVert(n2) >= 0))
	  {
	    // edge translated
	    int n1b = NumTranslat(nb_layers-1, IndexVert(n1));
	    int n2b = NumTranslat(nb_layers-1, IndexVert(n2));
            if (circle_layer)
              {
                mesh.BoundaryRef(i).SetReference(RefUsed(ref));
                mesh.BoundaryRef(Nbound).Init(n1b, n2b, ref);
                Nbound++;
              }
            else
              mesh.BoundaryRef(i).Init(n1b, n2b, ref);
	  }
	else if (IndexVert(n1) >= 0)
	  {
	    if (OriginalVertex(IndexVert(n1)))
	      {
		// new boundary edges
		int n1b = n1;
		int n2b = NumTranslat(0, IndexVert(n1));
		for (int k = 0; k < nb_layers; k++)
		  {
		    mesh.BoundaryRef(Nbound+k).Init(n1b, n2b, ref);
		    n1b = n2b;
		    if (k < nb_layers-1)
		      n2b = NumTranslat(k+1, IndexVert(n1));
		  }
		
		Nbound += nb_layers;
	      }
	  }
	else if (IndexVert(n2) >= 0)
	  {
	    if (OriginalVertex(IndexVert(n2)))
	      {
		// new boundary edges
		int n1b = n2;
		int n2b = NumTranslat(0, IndexVert(n2));
		for (int k = 0; k < nb_layers; k++)
		  {
		    mesh.BoundaryRef(Nbound+k).Init(n1b, n2b, ref);
		    n1b = n2b;
		    if (k < nb_layers-1)
		      n2b = NumTranslat(k+1, IndexVert(n2));
		  }
		
		Nbound += nb_layers;
	      }
	  }
      }
    
    mesh.ReorientElements();
    mesh.FindConnectivity();
    mesh.ProjectPointsOnCurves();
    //mesh.Write("toto.mesh");
  }
  
  
  //! displays details about class Mesh2D
  template<class T>
  ostream& operator <<(ostream& out, const Mesh<Dimension2, T>& mesh)
  {
    cout << "Mesh with " << mesh.GetNbElt() << " elements and " << mesh.GetNbBoundaryRef()
	 << " referenced edges, " << mesh.GetNbEdges() << " edges " << endl;
    
    return out;
  }
    
}

#define MONTJOIE_FILE_MESH2D_CXX
#endif
