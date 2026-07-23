#ifndef MONTJOIE_FILE_MESH_BASE_CXX

namespace Montjoie
{

  /********************
   * GetDimensionMesh *
   ********************/
  
  
  //! returns dimension of a mesh (2 or 3)
  /*!
    \param[in] name file name where the mesh is written
    \param[in] extension extension of the file name
   */
  int GetDimensionMesh(const string& name, const string& extension)
  {
    if (!extension.compare("mit"))
      {
	return 3;
      }
    else if (!extension.compare("dat"))
      {
	return 3;
      }
    else if (!extension.compare("gts"))
      {
	return 3;
      }
    else if (!extension.compare("nas"))
      {
	return 3;
      }
    else if (!extension.compare("stl"))
      {
	return 3;
      }
    else if (!extension.compare("mh3"))
      {
	return 3;
      }
    else if (!extension.compare("mhs"))
      {
	return 3;
      }
    else if (!extension.compare("vtk"))
      {
	ifstream file_in(name.data());
	if (!file_in.is_open())
	  {
	    cout << "Mesh file not found" << endl;
	    cout << name << endl;
	    abort();
	  }
	
	string ligne;
	getline(file_in, ligne);
	getline(file_in, ligne);
	file_in.close();
	if (ligne.find("2-D") != string::npos)
	  return 2;
	 
	return 3;
      }
    else if (!extension.compare("neu"))
      {
	return 3;
      }
    else if (!extension.compare("neub"))
      {
	return 3;
      }
    else if (!extension.compare("hmascii"))
      {
	return 3;
      }
    else if (!extension.compare("exo"))
      {
	ifstream file_in(name.data());
	if (!file_in.is_open())
	  {
	    cout << "Mesh file not found" << endl;
	    cout << name << endl;
	    abort();
	  }
	
	string ligne; Vector<string> parameters;
	bool test_loop = true; int dim_N = 2;
	while (test_loop)
	  {
	    getline(file_in, ligne);
	    
	    if ((ligne.size()!=0) && (ligne[0]!='#') && (ligne[0]!='/'))
	      {
		StringTokenize(ligne, parameters, string(" \t"));
		if (!parameters(0).compare("num_dim"))
		  {
		    dim_N = to_num<int>(parameters(2));
		    test_loop = false;
		  }
	      }
	  }
	
	return dim_N;	    
      }
    else if (!extension.compare("vf"))
      {
	return 3;
      }
    else if (!extension.compare("obj"))
      {
	return 3;
      }
    else if (!extension.compare("msh"))
      {
	ifstream file_in(name.data());
	if (!file_in.is_open())
	  {
	    cout << "Mesh file not found" << endl;
	    cout << name << endl;
	    abort();
	  }
	
	string ligne; double xtmp; int itmp, nb_vertices;
	file_in >> ligne;
	if (!ligne.compare("$MeshFormat"))
	  {
	    file_in >> xtmp >> itmp >> itmp;
	    file_in >> ligne;
	    file_in >> ligne;
	  }
	
	file_in >> nb_vertices;
	// lecture des sommets
	double x, y, z, zmoy(0); bool all_vertices_on_plane = true;
	for (int i = 0; i < nb_vertices; i++)
	  {
	    file_in >> itmp >> x >> y >> z;
	    if (i == 0)
	      zmoy = z;
	    else
	      {
		if (abs(zmoy-z) > 1e-12)
		  all_vertices_on_plane = false;
	      }
	  }
        
	file_in.close();
	
	if (all_vertices_on_plane)
	  return 2;
	else
	  return 3;
      }
    else if (!extension.compare("meshb"))
      {
	ifstream file_in(name.data());
	if (!file_in.is_open())
	  {
	    cout << "Mesh file not found" << endl;
	    cout << name << endl;
	    abort();
	  }
	
	// on lit le "endian" et la version
	int version, endian, dimension;
	file_in.read(reinterpret_cast<char*>(&endian), sizeof(int));
	file_in.read(reinterpret_cast<char*>(&version), sizeof(int));
	
	// boucle sur tous les mots-cles
	int mot_cle, new_position = 1;
	while ((!file_in.eof())&&(new_position!=0))
	  {
	    // on lit le mot cle et la position du fichier pour atteindre nouveau mot cle
	    file_in.read(reinterpret_cast<char*>(&mot_cle), sizeof(int));
	    file_in.read(reinterpret_cast<char*>(&new_position), sizeof(int));
	    
	    if (mot_cle == 3)
	      {
		// dimension
		file_in.read(reinterpret_cast<char*>(&dimension), sizeof(int));
	      }
	    
	    // on passe a la nouvelle position
	    file_in.seekg(new_position, ios::beg);
	  }
	
	file_in.close();
	return dimension;
      }
    else if (!extension.compare("mesh"))
      {
	ifstream file_in(name.data()); int dimension = 2;
	if (!file_in.is_open())
	  {
	    cout << "Mesh file not found" << endl;
	    cout << name << endl;
	    abort();
	  }
	
	string ligne;
	while (!file_in.eof())
	  {
	    getline(file_in, ligne);
	    DeleteSpaceAtExtremityOfString(ligne);
	    
	    if (!ligne.compare("Dimension"))
	      {
		file_in >> dimension;
		file_in.close();
		return dimension;
	      }
	  }
	
	file_in.close();
	return dimension;
	
      }
    
    return 2;
  }
  
  
  //! returns dimension of a mesh (2 or 3)
  /*!
    \param[in] name file name where the mesh is written
  */
  int GetDimensionMesh(const string& name)
  {
    return GetDimensionMesh(name, GetExtension(name));
  }

    
  /****************
   * MESH_BASE    *
   ****************/
  

  //! default constructor
  template<class Dimension, class T>
  Mesh_Base<Dimension, T>::Mesh_Base()
  {
    print_level = -1;
    
    mesh_path = string("./");
    irregular_mesh_double_cut = false;
    
    required_mesh_type = -1;

    reference_original_neighbor = 0;
    xmin_ = 0; xmax_ = 0;
    ymin_ = 0; ymax_ = 0;
    zmin_ = 0; zmax_ = 0;
  }
    
    
  //! returns size of memory used by the object in bytes
  template<class Dimension, class T>
  size_t Mesh_Base<Dimension, T>::GetMemorySize() const
  {
    size_t taille = Vertices.GetMemorySize();
    taille += sizeof(Vector<T>)*3;
    for (int i = 0; i < edges.GetM(); i++)
      taille += edges(i).GetMemorySize();

    for (int i = 0; i < edges_ref.GetM(); i++)
      taille += edges_ref(i).GetMemorySize();

    for (int i = 0; i < elements.GetM(); i++)
      taille += elements(i).GetMemorySize();
    
    taille += reference_element.GetMemorySize();
    taille += affine_element.GetMemorySize();
    taille += ListeVertices_to_be_added.GetMemorySize() +
      ListeVertices_to_be_refined.GetMemorySize() +
      LevelRefinement_Vertex.GetMemorySize() + RatioRefinement_Vertex.GetMemorySize();
    
    taille += head_minv.GetMemorySize() + head_maxv.GetMemorySize();
    taille += next_edge.GetMemorySize() + next_edge_arr.GetMemorySize();
    taille += reference_original_neighbor.GetMemorySize();
    
    taille += Seldon::GetMemorySize(Glob_invSqrtJacobian)
      + Seldon::GetMemorySize(Glob_invSqrtJacobianNeighbor);

    taille += Seldon::GetMemorySize(Glob_GradJacobian);
    
    taille += GlobElementNumber_Subdomain.GetMemorySize();
    taille += GlobFaceNumber_Subdomain.GetMemorySize();
    taille += GlobEdgeNumber_Subdomain.GetMemorySize();
    taille += GlobVertexNumber_Subdomain.GetMemorySize();
    return taille;
  }
  
  
  //! returns vertex number with coordinates equal to pt_glob
  /*!
    \param[in] pt_glob coordinates of point to find
  */
  template<class Dimension, class T>
  int Mesh_Base<Dimension, T>::FindVertexNumber(const R_N& pt_glob) const
  {
    for (int i = 0; i < this->GetNbVertices(); i++)
      if ( pt_glob == this->Vertex(i))
	return i;
    
    return -1;
  }
  

  //! changes the number of edges, previous edges are lost
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::ReallocateEdges(int nb_edges)
  {
    // next_edge is a linked list
    // next_edge(e) is the following edge leaving from the same point than the edge e
    // if there is no other edge, next_edge(e) will be equal to -1
    next_edge.Reallocate(nb_edges); next_edge.Fill(-1); head_minv.Fill(-1);
    // next_edge_arr(e) is the following edge arriving to the same point than the edge e
    // if there is no other edge, next_edge_arr(e) will be equal to -1
    next_edge_arr.Reallocate(nb_edges); next_edge_arr.Fill(-1); head_maxv.Fill(-1);
    // The arrays of edges
    edges.Reallocate(nb_edges);
  }


  //! changes the number of edges, previous edges are kept
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::ResizeEdges(int nb_edges)
  {
    int nb_old = edges.GetM();
    
    // next_edge is a linked list
    // next_edge(e) is the following edge leaving from the same point than the edge e
    // if there is no other edge, next_edge(e) will be equal to -1
    next_edge.Resize(nb_edges);
    // next_edge_arr(e) is the following edge arriving to the same point than the edge e
    // if there is no other edge, next_edge_arr(e) will be equal to -1
    next_edge_arr.Resize(nb_edges);
    // The arrays of edges
    edges.Resize(nb_edges);
    
    for (int i = nb_old; i < nb_edges; i++)
      {
	next_edge(i) = -1;
	next_edge_arr(i) = -1;
      }
  }
  

  //! changes the size of the array containing the vertices, previous vertices are lost
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::ReallocateVertices(int i)
  {
    Vertices.Reallocate(i);
    // we allocate arrays
    // head_minv will contain for each point
    // the number of the first edge who leaves this point
    head_minv.Reallocate(i); head_minv.Fill(-1);
    // head_maxv will contain for each point
    // the number of the first edge who comes to this point
    head_maxv.Reallocate(i); head_maxv.Fill(-1);
  }

  
  //! changes the size of the array containing the vertices, previous vertices are kept
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::ResizeVertices(int nb_vert)
  {
    int nb_old = Vertices.GetM();
    Vertices.Resize(nb_vert);
    head_minv.Resize(nb_vert);
    head_maxv.Resize(nb_vert);
    for (int i = nb_old; i < nb_vert; i++)
      {
	head_minv(i) = -1;
	head_maxv(i) = -1;
      }
  }
  
  
  //! gets list of vertices of element num_elt
  /*!
    \param[in] num_elt element number
    \param[out] points vertices of the element
   */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::GetVerticesElement(int num_elt,
						   VectR_N& points) const
  {
    int nb_v = elements(num_elt).GetNbVertices();
    points.Reallocate(nb_v);
    for (int i = 0; i < nb_v; i++)
      points(i) = Vertices(elements(num_elt).numVertex(i));
  }  

  
  //! gets list of "dof" points of element num_elt
  /*!
    \param[in] num_elt element number
    \param[out] points position of dofs (for interpolatory basis)
    \param[in] FaceBasis finite element 
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::GetDofsElement(int num_elt, VectR_N& points,
					       const ElementGeomReference<Dimension>& FaceBasis) const
  {
    VectR_N s;
    SetPoints<Dimension> PointsElem;
    
    this->GetVerticesElement(num_elt, s);
    FaceBasis.FjElemDof(s, PointsElem, this->GetLeafClass(), num_elt);
    
    points.Reallocate(FaceBasis.GetNbPointsDof());
    for (int i = 0; i < FaceBasis.GetNbPointsDof(); i++)
      points(i) = PointsElem.GetPointDof(i);    
  }
    
  
  //! modifies parameters of the mesh with a line of the data file
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("MeshPath"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh_Base" << endl;
	    cout << "MeshPath needs more parameters, for instance :" << endl;
	    cout << "MeshPath = path" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// the path where the mesh files are
	mesh_path = string(parameters(0));
      }
    else if (!description_field.compare("IrregularMesh"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh_Base" << endl;
	    cout << "IrregularMesh needs more parameters, for instance :" << endl;
	    cout << "IrregularMesh = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("YES"))
	  irregular_mesh_double_cut = true;
	else
	  irregular_mesh_double_cut = false;
      }
    else if (!description_field.compare("RefinementVertex"))
      {
	if (parameters.GetM() <= Dimension::dim_N+1)
	  {
	    cout << "In SetInputData of Mesh_Base" << endl;
	    cout << "RefinementVertex needs more parameters, for instance :" << endl;
	    cout << "RefinementVertex = xA yA level ratio" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int n = this->ListeVertices_to_be_refined.GetM();
	this->ListeVertices_to_be_refined.Resize(n + 1);
	this->LevelRefinement_Vertex.Resize(n + 1);
	this->RatioRefinement_Vertex.Resize(n + 1);
	for (int k = 0; k < Dimension::dim_N; k++)
	  this->ListeVertices_to_be_refined(n)(k) = to_num<T>(parameters(k));
	
	int nb = Dimension::dim_N;
	this->LevelRefinement_Vertex(n) = to_num<int>(parameters(nb++));
	this->RatioRefinement_Vertex(n) = to_num<Real_wp>(parameters(nb++));
      }
    else if (!description_field.compare("RefinementArea"))
      {
	ParamRefinement_Area = parameters;
      }
    else if (!description_field.compare("AddVertex"))
      {
	if (parameters.GetM() < Dimension::dim_N)
	  {
	    cout << "In SetInputData of Mesh_Base" << endl;
	    cout << "AddVertex needs more parameters, for instance :" << endl;
	    cout << "AddVertex = xA yA" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        typename Dimension::R_N pt;
	for (int k = 0; k < Dimension::dim_N; k++)
	  pt(k) = to_num<T>(parameters(k));
        
	this->ListeVertices_to_be_added.PushBack(pt);
      }
    else if (!description_field.compare("ThresholdMesh"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh_Base" << endl;
	    cout << "ThresholdMesh needs more parameters, for instance :" << endl;
	    cout << "ThresholdMesh = eps" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	R_N::threshold = to_num<T>(parameters(0));
      }
    else if (!description_field.compare("PrintLevel"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh_Base" << endl;
	    cout << "PrintLevel needs more parameters, for instance :" << endl;
	    cout << "PrintLevel = level" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	print_level = to_num<int>(parameters(0));
      }
  }
  
  
  // return true if two vertices are on the same edge
  template<class Dimension, class T>
  bool Mesh_Base<Dimension, T>::PointsOnSameEdge(int i, int ni, int& ne) const
  {
    // loop on edges leaving point i
    for (int e = head_minv(i); e != -1; e = next_edge(e))
      for (int j = 0; j < edges(e).GetNbElements(); j++)
	if ((edges(e).numVertex(0) == ni)||(edges(e).numVertex(1) == ni))
	  {
	    ne = e;
	    return true;
	  }
    
    // loop on edges coming to point i
    for (int e = head_maxv(i); e != -1; e = next_edge_arr(e))
      for (int j = 0; j < edges(e).GetNbElements(); j++)
	if ((edges(e).numVertex(0) == ni)||(edges(e).numVertex(1) == ni))
	  {
	    ne = e;
	    return true;
	  }
    
    return false;
  }

  
  //! returns list of elements near a vertex of the mesh
  /*!
    \param[in] i  vertex numbers
    \param[out] num list of elements
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::GetNeighboringElementsAroundVertex(int i, IVect& num) const
  {
    num.Clear();
    std::vector<int> num_v; num_v.reserve(100);
    // loop on edges leaving point i
    for (int e = head_minv(i); e != -1; e = next_edge(e))
      for (int j = 0; j < edges(e).GetNbElements(); j++)
	num_v.push_back(edges(e).numElement(j));
    
    // loop on edges coming to point i
    for (int e = head_maxv(i); e != -1; e = next_edge_arr(e))
      for (int j = 0; j < edges(e).GetNbElements(); j++)
	num_v.push_back(edges(e).numElement(j));
    
    if (num_v.size() <= 0)
      return;
    
    num.SetData(num_v.size(), &num_v[0]);
    int nb_elt_around = num_v.size();
    // eliminating redondant numbers
    Assemble(nb_elt_around, num);
    // we fill num
    num.Nullify(); num.Reallocate(nb_elt_around);
    for (int i = 0; i < nb_elt_around; i++)
      num(i) = num_v[i];
    
  }
  
  
  //! returns list of edges near a vertex of the mesh
  /*!
    \param[in] i  vertex numbers
    \param[out] num list of edges
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::GetNeighboringEdgesAroundVertex(int i, IVect& num) const
  {
    num.Clear();
    std::vector<int> num_v; num_v.reserve(100);
    // loop on edges leaving point i
    for (int e = head_minv(i); e != -1; e = next_edge(e))
      num_v.push_back(e);
    
    // loop on edges coming to point i
    for (int e = head_maxv(i); e != -1; e = next_edge_arr(e))
      num_v.push_back(e);
    
    if (num_v.size() <= 0)
      return;
    
    num.Reallocate(num_v.size());
    for (int i = 0; i < num.GetM(); i++)
      num(i) = num_v[i];
    
  }
  

  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::SetVerticesToBeAdded(const VectR_N& points)
  {
    ListeVertices_to_be_added = points;
  }

  
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::SetVerticesToBeRefined(const VectR_N& points, int lvl, const T& alpha)
  {
    ListeVertices_to_be_refined = points;
    LevelRefinement_Vertex.Reallocate(points.GetM());
    LevelRefinement_Vertex.Fill(lvl);
    RatioRefinement_Vertex.Reallocate(points.GetM());
    RatioRefinement_Vertex.Fill(alpha);
  }

  
  //! we find a reference for each vertex according to references of edges/faces
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::CreateReferenceVertices(IVect& ref_vertex) const
  {
    const Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    ref_vertex.Reallocate(this->GetNbVertices());
    ref_vertex.Fill(0);
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      for (int j = 0; j < mesh_leaf.BoundaryRef(i).GetNbVertices(); j++)
	ref_vertex(mesh_leaf.BoundaryRef(i).numVertex(j))
	  = mesh_leaf.BoundaryRef(i).GetReference();
    
  }
  
  
  //! returns the average length of the edges of the mesh
  template<class Dimension, class T>
  T Mesh_Base<Dimension, T>::GetMeshSize() const
  {
    if (this->GetNbEdges() == 0)
      {
	cout << "No edges in the mesh " << endl;
	cout << "Did you call FindConnectivity() ? " << endl;
	abort();
      }
    
    T space_step(0);
    for (int i = 0; i < this->GetNbEdges(); i++)
      {
	int n0 = edges(i).numVertex(0), n1 = edges(i).numVertex(1);
	space_step += Vertices(n0).Distance(Vertices(n1));
      }
    
    space_step = space_step/this->GetNbEdges();
    return space_step;
  }
  
  
  //! connectivity of edges is cleared
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::ClearConnectivity()
  {
    for (int i = 0; i < this->GetNbElt(); i++)
      this->Element(i).ClearConnectivity();
    
    head_minv.Fill(-1); head_maxv.Fill(-1);
    next_edge.Clear(); next_edge_arr.Clear();
    edges.Clear(); 
    for (int i = 0; i < edges_ref.GetM(); i++)
      edges_ref(i).ClearConnectivity();
    
  }


  //! clears input data used to construct the mesh
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::ClearInput()
  {
    irregular_mesh_double_cut = false;
    required_mesh_type = 0;
    ListeVertices_to_be_added.Clear();
    ListeVertices_to_be_refined.Clear();
    LevelRefinement_Vertex.Clear();
    RatioRefinement_Vertex.Clear();
  }
  
  
  //! all the mesh is removed
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::Clear()
  {
    head_minv.Clear(); head_maxv.Clear();
    next_edge.Clear(); next_edge_arr.Clear();
    
    Vertices.Clear(); edges.Clear(); edges_ref.Clear();
    elements.Clear(); affine_element.Clear();    
  }
  
  
  //! internal method
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::AddEdgeOfElement(int& nb_edges, int i, int na)
  {
    int type_elt = elements(i).GetHybridType();
    // local numbers
    int iloc = MeshNumbering<Dimension, T>::FirstExtremityEdge(type_elt, na);
    int jloc = MeshNumbering<Dimension, T>::SecondExtremityEdge(type_elt, na);
    
    // we get global numbers of the two extremities of the local edge
    int iglob = elements(i).numVertex(iloc);
    int jglob = elements(i).numVertex(jloc);
    
    // We consider the edge iglob -> jglob 
    // we sort the numbers
    if (iglob > jglob)
      {
	int itmp = iglob;
	iglob = jglob;
	jglob = itmp;
      }
    
    bool existe = false;
    // we look at the list of edges leaving from iglob
    // if we find the same edge than iglob->jglob we add the face i
    // to this edge
    for (int e = head_minv(iglob); e != -1; e = next_edge(e))
      {
	// the edge iglob->jglob already exists
	// we update edges(e) and faces(i)
	if ( edges(e).numVertex(1) == jglob)
	  {
	    edges(e).AddElement(i);
	    elements(i).SetEdge(na, e);
	    existe = true;
	  }
      }
    
    // if the edge has not been found, we create it
    if (!existe) 
      {
	edges(nb_edges).AddElement(i);
	// we initialize the edge 
	edges(nb_edges).Init(iglob,jglob,0);
	// we add it to the element i
	elements(i).SetEdge(na, nb_edges);
		
	// we update the arrays next_edge, head_minv, etc
	next_edge(nb_edges) = head_minv(iglob);
	head_minv(iglob) = nb_edges;
	next_edge_arr(nb_edges) = head_maxv(jglob);
	head_maxv(jglob) = nb_edges;
	nb_edges++;
      }
  }
  
  
  //! internal method
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::RemoveSingleEdge(int num_edge)
  {
    int n1 = this->GetEdge(num_edge).numVertex(0);
    int n2 = this->GetEdge(num_edge).numVertex(1);
    // we remove num_edge from linked list (head_minv, next_edge)
    if (head_minv(n1) == num_edge)
      head_minv(n1) = next_edge(num_edge);
    else
      {
	int eprec = head_minv(n1), enext = next_edge(eprec);
	while ((enext != num_edge)&&(enext != -1))
	  {
	    eprec = enext;
	    enext = next_edge(eprec);
	  }
	next_edge(eprec) = next_edge(enext);
      }
    
    // we remove num_edge from linked list (head_maxv, next_edge_arr)
    if (head_maxv(n2) == num_edge)
      head_maxv(n2) = next_edge_arr(num_edge);
    else
      {
	int eprec = head_maxv(n2), enext = next_edge_arr(eprec);
	while ((enext != num_edge)&&(enext != -1))
	  {
	    eprec = enext;
	    enext = next_edge_arr(eprec);
	  }
	next_edge_arr(eprec) = next_edge_arr(enext);
      }
  }
  
  
  //! internal method
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::AddSingleEdge(int& nb_edges, const Edge<Dimension>& edge)
  {
    int n1 = edge.numVertex(0);
    int n2 = edge.numVertex(1);
    
    bool existe = false;
    // we look at the list of edges leaving from n1
    // if we find the same edge than n1->n2 we add the edge
    for (int e = head_minv(n1); e != -1; e = next_edge(e))
      {
	if ( edges(e).numVertex(1) == n2)
	  existe = true;
      }
    
    // if the edge has not been found, we create it
    if (!existe) 
      {
	// we initialize the edge 
	edges(nb_edges) = edge;
	
	// we update the arrays next_edge, head_minv, etc
	next_edge(nb_edges) = head_minv(n1);
	head_minv(n1) = nb_edges;
	next_edge_arr(nb_edges) = head_maxv(n2);
	head_maxv(n2) = nb_edges;
	nb_edges++;
      }
  }
  
  
  //! seeks all the edges in the mesh
  /*!
    This methods finds all the edges of the mesh,
    from the list of the elements. The method
    is quasi-optimal (F. Hecht's algorithm).
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::AddEdgesFromElements()
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    
    // previous edges
    int nb_old = mesh_leaf.GetNbEdges();
    
    // we compute a max value for the number of edges due to elements
    int nb_edges = mesh_leaf.GetUpperBoundNumberOfEdges();
    
    mesh_leaf.ResizeEdges(nb_old + nb_edges);
    nb_edges = nb_old;
    
    if (print_level >= 6)
      {
        cout<< " Maximum number of edges "<<edges.GetM()<<endl;
        cout<<" Number of elements "<<this->GetNbElt()<<endl;
      }
	
    // we make a loop on all elements
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	int nb_edges_loc = elements(i).GetNbEdges();
	
	for (int na = 0; na < nb_edges_loc; na++)
	  AddEdgeOfElement(nb_edges, i, na);
	
      }
    
    // we resize the arrays
    ResizeEdges(nb_edges);
    
    if (print_level >= 6)
      cout << " Edges constructed "<<endl;
    
    if (print_level >= 6)
      cout << " Number of edges in the mesh  "<<nb_edges<<endl;
  } 
     
  
  //! extends a subdomain with neighboring elements
  /*!
    \param[in] delta_overlapping depth level, if equal to one, only direct neighbors are added
    \param[inout] nb_vertices_subdomain number of vertices before and after
    \param[inout] nb_elt_subdomain number of elements before and after
    \param[inout] VertexOnSubdomain VertexOnSubdomain(i) is true 
    if the vertex i is part of the subdomain
    \param[inout] ElementOnSubdomain ElementOnSubdomain(i) is true 
    if the element i is part of the subdomain
   */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::
  AddOverlapRegion(int delta_overlapping, int& nb_vertices_subdomain, int& nb_elt_subdomain, 
		   Vector<bool>& VertexOnSubdomain, Vector<bool>& ElementOnSubdomain) const
  {
    // eventual overlapping
    for (int num_overlap = 0; num_overlap < delta_overlapping; num_overlap++)
      {
	// creation of vertices of the subdomain
	IVect ListVertex(nb_vertices_subdomain);
	int index = 0;
	for (int i = 0; i < this->GetNbVertices(); i++)
	  if (VertexOnSubdomain(i))
	    ListVertex(index++) = i;
	
	IVect num;
	// loop on all vertices of the subdomain
	for (index = 0; index < ListVertex.GetM(); index++)
	  {
	    int i = ListVertex(index);
	    // we get element numbers near this vertex
	    GetNeighboringElementsAroundVertex(i, num);
	    for (int j = 0; j < num.GetM(); j++)
	      if (!ElementOnSubdomain(num(j)))
		{
		  // we get an element outside the subdomain
		  // but with a point on the subdomain
		  // we add it to the subdomain
		  ElementOnSubdomain(num(j)) = true;
		  int nb_vertices_elt = elements(num(j)).GetNbVertices();
		  // we add also the vertices of the element
		  for (int n = 0; n < nb_vertices_elt; n++)
		    if (!VertexOnSubdomain(elements(num(j)).numVertex(n)))
		      {
			VertexOnSubdomain(elements(num(j)).numVertex(n)) = true;
			nb_vertices_subdomain++;
		      }
		  nb_elt_subdomain++;
		}
	  }
      }
  }
  

  //! the mesh is split into small meshes, thanks to a regular grid
  /*!
    \param[in] nb_subdivX number of subdivisions of the grid along x
    \param[in] nb_subdivY number of subdivisions of the grid along y
    \param[in] nb_subdivZ number of subdivisions of the grid along z
    \param[out] nb_subdomains number of subdomains created
    \param[out] NumElement_Subdomain list of elements for each subdomain
    \param[out] sub_mesh list of the submeshes
    \param[out] delta_overlapping level of overlapping
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::
  SplitIntoBoxes(int nb_subdivX, int nb_subdivY, int nb_subdivZ, int& nb_subdomains,
		 Vector<IVect>& NumElement_Subdomain, Vector<Mesh<Dimension, T> >& sub_mesh,
		 int delta_overlapping, bool change_periodic_ref)
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    VectBool ElementUsed(this->GetNbElt()); ElementUsed.Fill(false);
    VectBool ElementOnSubdomain(this->GetNbElt());
    T step_x, step_y, step_z, zmin(0), zmax(0);
    
    Real_wp xmin = this->GetXmin(); Real_wp xmax = this->GetXmax();
    Real_wp ymin = this->GetYmin(); Real_wp ymax = this->GetYmax();
    if (Dimension::dim_N == 3)
      {
	zmin = this->GetZmin();
	zmax = this->GetZmax();
      }
    
    if (nb_subdivX <= 1)
      step_x = xmax-xmin;
    else
      step_x = (xmax-xmin)/nb_subdivX;
    
    if (nb_subdivY <= 1)
      step_y = ymax-ymin;
    else
      step_y = (ymax-ymin)/nb_subdivY;
    
    if (nb_subdivZ <= 1)
      step_z = zmax-zmin;
    else
      step_z = (zmax-zmin)/nb_subdivZ;
    
    T x0, y0, z0, x1, y1, z1;
    VectR_N s(8); R_N center;
    
    VectBool VertexOnSubdomain(this->GetNbVertices());
    
    nb_subdomains = nb_subdivX*nb_subdivY*nb_subdivZ;
    
    sub_mesh.Reallocate(nb_subdomains);
    NumElement_Subdomain.Reallocate(nb_subdomains);
    int num_subdomain = 0;
    
    for (int i = 0; i < nb_subdivX; i++)
      for (int j = 0; j < nb_subdivY; j++)
	for (int k = 0; k < nb_subdivZ; k++)
	  {
	    // definition of the box where we get subdomain
	    x0 = xmin + i*step_x;
	    y0 = ymin + j*step_y;
	    z0 = zmin + k*step_z;
	    x1 = x0 + step_x;
	    y1 = y0 + step_y;
	    z1 = z0 + step_z;
	    
	    ElementOnSubdomain.Fill(false);
	    VertexOnSubdomain.Fill(false);
	    int nb_elt_subdomain = 0, nb_vertices_subdomain = 0;
	    for (int e = 0; e < this->GetNbElt(); e++)
	      {
		GetVerticesElement(e, s);
		int nb_vertices_elt = elements(e).GetNbVertices();
		center.Fill(0);
		for (int n = 0; n < nb_vertices_elt; n++)
		  center += s(n);
		
		Mlt(T(1)/nb_vertices_elt, center);
		if ((center(0) >= x0)&&(center(0) < x1))
		  if ((center(1) >= y0)&&(center(1) < y1))
		    if ( (Dimension::dim_N == 2)|| ((center(2) >= z0)&&(center(2) < z1)) )
		      {
			// DISP(e); DISP(ElementUsed(e));
			// element in the subdomain
			if (!ElementUsed(e))
			  {
			    ElementOnSubdomain(e) = true;
			    ElementUsed(e) = true;
			    for (int n = 0; n < nb_vertices_elt; n++)
			      if (!VertexOnSubdomain(elements(e).numVertex(n)))
				{
				  VertexOnSubdomain(elements(e).numVertex(n)) = true;
				  nb_vertices_subdomain++;
				}
			    nb_elt_subdomain++;
			  }
		      }
	      }
	    
	    this->AddOverlapRegion(delta_overlapping, nb_vertices_subdomain, nb_elt_subdomain, 
				   VertexOnSubdomain, ElementOnSubdomain);
	    
	    // DISP(nb_elt_subdomain);
	    if (nb_elt_subdomain > 0)
	      {
		NumElement_Subdomain(num_subdomain).Reallocate(nb_elt_subdomain);
		nb_elt_subdomain = 0;
		for (int e = 0; e < this->GetNbElt(); e++)
		  if (ElementOnSubdomain(e))
		    NumElement_Subdomain(num_subdomain)(nb_elt_subdomain++) = e;
		
		mesh_leaf.CreateSubmesh(sub_mesh(num_subdomain), nb_vertices_subdomain,
					nb_elt_subdomain,
					VertexOnSubdomain, ElementOnSubdomain);
		
		// sub_mesh(num_subdomain).Write("subdomain.mesh");
		// cout<<"we wait"<<endl; int test_input; cin>>test_input;
		
		num_subdomain++;
	      }
	  }
    
    nb_subdomains = num_subdomain;
  }
  
  
  //! partitioning of a mesh with metis functions
  /*!
    \param[in] nb_parts wished number of sub-domains
    \param[out] NumElement_Subdomain list of element numbers for each subdomain
    \param[out] sub_mesh list of the small meshes
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::
  SplitMetis(int nb_parts, const IVect& weight_elt, Vector<IVect>& NumElement_Subdomain,
	     Vector<Mesh<Dimension, T> >& sub_mesh,
	     int delta_overlapping, bool change_periodic_ref)
  {
#ifndef MONTJOIE_WITHOUT_METIS
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    if (this->GetNbElt() <= 0)
      return;
    
    Vector<metis_int_t> epart(this->GetNbElt());
    
    // constructing adjacency graph for Metis (equivalent to CSR format)
    Vector<metis_int_t> xadj(this->GetNbElt()+1); xadj(0) = 0; int nf = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	xadj(i+1) = xadj(i);
	for (int j = 0; j < mesh_leaf.Element(i).GetNbBoundary(); j++)
	  {
	    nf = mesh_leaf.Element(i).numBoundary(j);
	    if (mesh_leaf.Boundary(nf).GetNbElements() == 2)
	      xadj(i+1)++;
	  }
      }
    
    Vector<metis_int_t> adjncy(xadj(this->GetNbElt()));
    int nb = 0, num_elem2;
    for (int i = 0; i < this->GetNbElt(); i++)
      for (int j = 0; j < mesh_leaf.Element(i).GetNbBoundary(); j++)
	{
	  nf = mesh_leaf.Element(i).numBoundary(j);
	  if (mesh_leaf.Boundary(nf).GetNbElements() == 2)
	    {
	      num_elem2 = mesh_leaf.Boundary(nf).numElement(0);
	      if (num_elem2 == i)
		num_elem2 = mesh_leaf.Boundary(nf).numElement(1);
	      
	      adjncy(nb) = num_elem2;
	      nb++;
	    }
	}
    
    Vector<metis_int_t> weight_element(weight_elt.GetM());
    for (int i = 0; i < weight_elt.GetM(); i++)
      weight_element(i) = weight_elt(i);

    metis_int_t nb_part = nb_parts;
    metis_int_t nb_elt = this->GetNbElt(), ncon = 1;
    metis_int_t total_comm;
    METIS_PartGraphKway(&nb_elt, &ncon, xadj.GetData(), adjncy.GetData(),
                        weight_element.GetData(), NULL, NULL,
                        &nb_part, NULL, NULL, NULL, &total_comm, epart.GetData());
    
    weight_element.Clear();
    Vector<int> Epart(epart.GetM());
    for (int i = 0; i < epart.GetM(); i++)
      Epart(i) = epart(i);
    
    PartMesh(nb_parts, Epart, NumElement_Subdomain, sub_mesh, delta_overlapping, change_periodic_ref);
#else
    cout << "Recompile with metis " << endl;
    abort();
#endif
  }
  
  
  //! partitioning of a mesh with Scotch functions
  /*!
    \param[in] nb_parts wished number of sub-domains
    \param[out] NumElement_Subdomain list of element numbers for each subdomain
    \param[out] sub_mesh list of the small meshes
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::
  SplitScotch(int nb_parts, const IVect& weight_elt, Vector<IVect>& NumElement_Subdomain,
              Vector<Mesh<Dimension, T> >& sub_mesh, int delta_overlap,
	      bool change_periodic_ref)
  {
#ifdef MONTJOIE_WITH_SCOTCH
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    if (this->GetNbElt() <= 0)
      return;
    
    Vector<SCOTCH_Idx> epart_int(this->GetNbElt());
    Vector<SCOTCH_Idx> weight_elt_int(weight_elt.GetM());
    for (int i = 0; i < weight_elt.GetM(); i++)
      weight_elt_int(i) = weight_elt(i);
    
    // constructing adjacency graph
    SCOTCH_Idx baseval = 0;
    SCOTCH_Idx vertnbr = this->GetNbElt();    
    //IVect weight_edge(this->GetNbBoundary());
    //weight_edge.Fill(1);
    
    Vector<SCOTCH_Idx> verttab(this->GetNbElt()+1); verttab(0) = 0;
    int nf = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	verttab(i+1) = verttab(i);
	for (int j = 0; j < mesh_leaf.Element(i).GetNbBoundary(); j++)
	  {
	    nf = mesh_leaf.Element(i).numBoundary(j);
	    if (mesh_leaf.Boundary(nf).GetNbElements() == 2)
	      verttab(i+1)++;
	  }
      }
    
    int edgenbr = verttab(this->GetNbElt());
    Vector<SCOTCH_Idx> edgetab(edgenbr);
    int nb = 0, num_elem2;
    for (int i = 0; i < this->GetNbElt(); i++)
      for (int j = 0; j < mesh_leaf.Element(i).GetNbBoundary(); j++)
	{
	  nf = mesh_leaf.Element(i).numBoundary(j);
	  if (mesh_leaf.Boundary(nf).GetNbElements() == 2)
	    {
	      num_elem2 = mesh_leaf.Boundary(nf).numElement(0);
	      if (num_elem2 == i)
		num_elem2 = mesh_leaf.Boundary(nf).numElement(1);
	      
	      edgetab(nb) = num_elem2;
	      nb++;
	    }
	}
 
    SCOTCH_Graph grafptr;
    SCOTCH_Strat straptr;
    
    SCOTCH_stratInit(&straptr);
    SCOTCH_graphBuild(&grafptr, baseval, vertnbr, verttab.GetData(),
                      &verttab(1), weight_elt_int.GetData(), NULL,
                      edgenbr, edgetab.GetData(), NULL);
    
    SCOTCH_Idx nb_parts_int = nb_parts;
    SCOTCH_graphPart(&grafptr, nb_parts_int, &straptr, epart_int.GetData());
    
    SCOTCH_graphExit(&grafptr);
    SCOTCH_stratExit(&straptr);
    
    Vector<int> epart(epart_int.GetM());
    for (int i = 0; i < epart.GetM(); i++)
      epart(i) = epart_int(i);
    
    PartMesh(nb_parts, epart, NumElement_Subdomain, sub_mesh, delta_overlap, change_periodic_ref);
#else
    cout << "Recompile Montjoie with Scotch" << endl;
    abort();
#endif
  }
  
  
  //! partitioning a mesh with concentric circles (spheres in 3-D)
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>
  ::SplitConcentric(int nb_parts, const VectReal_wp& radius,
                    Vector<IVect>& NumElement_Subdomain,
                    Vector<Mesh<Dimension, T> >& sub_mesh,
		    int delta_overlap, bool change_periodic_ref)
  {
    IVect epart(this->GetNbElt());    
    typename Dimension::R_N center;
    for (int i = 0; i < epart.GetM(); i++)
      {
        center.Fill(0);
        for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
          center += this->Vertex(this->Element(i).numVertex(j));
        
        Real_wp coef = 1.0/this->Element(i).GetNbVertices();
        center *= coef;
        Real_wp r = Norm2(center);
        
        int k = 0;
        while ((k < radius.GetM()) && (r > radius(k)))
          k++;
        
        epart(i) = k;
      }
    
    PartMesh(nb_parts, epart, NumElement_Subdomain, sub_mesh, delta_overlap, change_periodic_ref);
  }


  //! partitioning a mesh with layered mesh
  /*!
    ref_layer(i) = references of the elements belonging to the i-th part
  */
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>
  ::SplitLayered(int nb_parts, const Vector<IVect>& ref_layer,
                 Vector<IVect>& NumElement_Subdomain,
                 Vector<Mesh<Dimension, T> >& sub_mesh, int delta_overlap,
		 bool change_periodic_ref)
  {
    IVect epart(this->GetNbElt());    
    typename Dimension::R_N center;
    // on prend la reciproque de ref_layer
    int ref_max = 0;
    for (int i = 0; i < ref_layer.GetM(); i++)
      ref_max = max(ref_max, ref_layer(i).GetNormInf());
    
    IVect inv_ref(ref_max+1);
    inv_ref.Fill(-1);
    for (int i = 0; i < ref_layer.GetM(); i++)
      for (int j = 0; j < ref_layer(i).GetM(); j++)
        inv_ref(ref_layer(i)(j)) = i;
    
    // filling epart
    for (int i = 0; i < epart.GetM(); i++)
      {
        if (this->Element(i).GetReference() >= inv_ref.GetM())
          {
            cout << "This reference cannot be treated " << endl;
            cout << "Too many references in the mesh ? " << endl;
            abort();
          }
        
        epart(i) = inv_ref(this->Element(i).GetReference());
        
        if (epart(i) == -1)
          {
            cout << "This element does not belong to any processor" << endl;
            abort();
          }
      }
    
    PartMesh(nb_parts, epart, NumElement_Subdomain, sub_mesh, delta_overlap, change_periodic_ref);
  }

  
  //! partitioning of a mesh by giving directly the epart array
  /*!
    \param[in] nb_parts wished number of sub-domains
    \param[in] epart epart(i) is the part number for element i
    \param[out] NumElement_Subdomain list of element numbers for each subdomain
    \param[out] sub_mesh list of the small meshes
  */  
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::
  PartMesh(int nb_parts, const IVect& epart, Vector<IVect>& NumElement_Subdomain,
           Vector<Mesh<Dimension, T> >& sub_mesh, int delta_overlapping,
	   bool change_periodic_ref)
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    int lvl = mesh_leaf.print_level;
    mesh_leaf.print_level = 0;

    // creating faces which are on the interface between subdomains
    int nb_face_on_interface = 0;
    for (int i = 0; i < mesh_leaf.GetNbBoundary(); i++)
      if (mesh_leaf.Boundary(i).GetNbElements() == 2)
        {
          int n1 = mesh_leaf.Boundary(i).numElement(0);
          int n2 = mesh_leaf.Boundary(i).numElement(1);
          int p1 = epart(n1);
          int p2 = epart(n2);
          if (p1 != p2)
            {
              if (i >= mesh_leaf.GetNbBoundaryRef())
                nb_face_on_interface++;
            }
        }
    
    int nb_old_faces = mesh_leaf.GetNbBoundaryRef();
    mesh_leaf.ResizeBoundariesRef(nb_old_faces + nb_face_on_interface);
    nb_face_on_interface = 0;
    IVect ref_to_change(mesh_leaf.GetNbReferences()+1);
    ref_to_change.Fill(0);
    int new_ref = 0, ref_m2 = 0, ref_m1 = 0;
    int cond_neighbor = BoundaryConditionEnum::LINE_NEIGHBOR;
    for (int i = 0; i < mesh_leaf.GetNbBoundary(); i++)
      if (mesh_leaf.Boundary(i).GetNbElements() == 2)
        {
          int n1 = mesh_leaf.Boundary(i).numElement(0);
          int n2 = mesh_leaf.Boundary(i).numElement(1);
          int p1 = epart(n1);
          int p2 = epart(n2);

          if (p1 != p2)
            {
              if (i >= nb_old_faces)
                {
                  if (new_ref == 0)
                    {
                      new_ref = mesh_leaf.GetNewReference();
                      reference_original_neighbor(new_ref) = 0;
                      mesh_leaf.SetBoundaryCondition(new_ref, cond_neighbor);
                    }
                  
                  int nf = nb_face_on_interface + nb_old_faces;                  
                  mesh_leaf.BoundaryRef(nf) = mesh_leaf.Boundary(i);
                  mesh_leaf.BoundaryRef(nf).SetReference(new_ref);
                  nb_face_on_interface++;
                }
              else
                {
                  int ref = mesh_leaf.BoundaryRef(i).GetReference();
                  if (ref >= 0)
                    {
                      int ref2 = ref_to_change(ref);
                      if (ref2 == 0)
                        {
                          ref2 = mesh_leaf.GetNewReference();
                          ref_to_change(ref) = ref2;
                          reference_original_neighbor(ref2) = ref;
                          mesh_leaf.Type_curve(ref2) = mesh_leaf.Type_curve(ref);
                          mesh_leaf.Body_curve(ref2) = mesh_leaf.Body_curve(ref);
                          if (mesh_leaf.Cond_curve(ref) > 0)
                            mesh_leaf.Cond_curve(ref2) = mesh_leaf.Cond_curve(ref);
                          else
                            mesh_leaf.Cond_curve(ref2) = cond_neighbor;
                          
                          mesh_leaf.Param_curve(ref2) = mesh_leaf.Param_curve(ref);
                          mesh_leaf.Param_known(ref2) = mesh_leaf.Param_known(ref);
                        }
                      
                      mesh_leaf.BoundaryRef(i).SetReference(ref2);
                      mesh_leaf.Boundary(i).SetReference(ref2);
                    }
                  else if (ref == -1)
                    {
                      if (ref_m1 == 0)
                        {
                          ref_m1 = mesh_leaf.GetNewReference();
                          reference_original_neighbor(ref_m1) = -1;
                          mesh_leaf.SetBoundaryCondition(ref_m1, cond_neighbor);
                        }
                      
                      mesh_leaf.Type_curve(ref_m1) = mesh_leaf.CURVE_FILE;
                      mesh_leaf.BoundaryRef(i).SetReference(ref_m1);
                      mesh_leaf.Boundary(i).SetReference(ref_m1);                      
                    }
                  else if (ref == -2)
                    {
                      if (ref_m2 == 0)
                        {
                          ref_m2 = mesh_leaf.GetNewReference();
                          reference_original_neighbor(ref_m2) = -2;
                          mesh_leaf.SetBoundaryCondition(ref_m2, cond_neighbor);
                        }
                      
                      mesh_leaf.Type_curve(ref_m2) = mesh_leaf.CURVE_FILE;
                      mesh_leaf.BoundaryRef(i).SetReference(ref_m2);
                      mesh_leaf.Boundary(i).SetReference(ref_m2);                      
                    }
                }
            }
        }

    // periodic references are set to neighbor
    if (change_periodic_ref)
      for (int i = 0; i < mesh_leaf.periodicity_references.GetM(); i++)
	{
	  mesh_leaf.SetBoundaryCondition(mesh_leaf.periodicity_references(i)(0), cond_neighbor);
	  mesh_leaf.SetBoundaryCondition(mesh_leaf.periodicity_references(i)(1), cond_neighbor);
	}
    
    // updating edges/faces
    mesh_leaf.FindConnectivity();
    
    // creating sub-meshes
    int nb_elt = mesh_leaf.GetNbElt();
    NumElement_Subdomain.Reallocate(nb_parts); sub_mesh.Reallocate(nb_parts);
    VectBool VertexOnSubdomain(this->GetNbVertices()), ElementOnSubdomain(nb_elt);
    for (int i = 0; i < nb_parts; i++)
      {
	VertexOnSubdomain.Fill(false); ElementOnSubdomain.Fill(false);
	// counting number of elements, vertices
	int nb_elt_subdomain = 0, nb_vert_subdomain = 0;
	for (int j = 0; j < nb_elt; j++)
	  if (epart(j) == i)
	    {
	      ElementOnSubdomain(j) = true;
	      for (int k = 0; k < elements(j).GetNbVertices(); k++)
		{
		  int num_vertex = elements(j).numVertex(k);
		  if (!VertexOnSubdomain(num_vertex))
		    {
		      VertexOnSubdomain(num_vertex) = true;
		      nb_vert_subdomain++;
		    }
		}
	      nb_elt_subdomain++;
	    }
        
	this->AddOverlapRegion(delta_overlapping, nb_vert_subdomain, nb_elt_subdomain, 
			       VertexOnSubdomain, ElementOnSubdomain);
	
 	// now filling NumElement_Subdomain
	NumElement_Subdomain(i).Reallocate(nb_elt_subdomain);
	nb_elt_subdomain = 0;
	for (int j = 0; j < nb_elt; j++)
	  if (ElementOnSubdomain(j))
	    NumElement_Subdomain(i)(nb_elt_subdomain++) = j;
	
	
	// creating the submesh
        if (lvl >= 10)
          sub_mesh(i).print_level = 6;
        else
          sub_mesh(i).print_level = 0;
        
	mesh_leaf.CreateSubmesh(sub_mesh(i), nb_vert_subdomain, nb_elt_subdomain,
				VertexOnSubdomain, ElementOnSubdomain);
        
      }

    mesh_leaf.print_level = lvl;
  }


  //! boundaries with periodic references are modified if needed (parallel computation only)
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::UpdatePeriodicReferenceSplit(const IVect& Epart, const IVect& NumLoc,
							     Vector<Mesh<Dimension, T> >& sub_mesh)
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    if (mesh_leaf.GetNbPeriodicReferences() <= 0)
      return;
    
    IVect ref_to_change(mesh_leaf.GetNbReferences()+1);
    ref_to_change.Fill(0);
    
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (mesh_leaf.Boundary(i).GetNbElements() == 2)
	{
	  int n1 = mesh_leaf.Boundary(i).numElement(0);
	  int num_loc = mesh_leaf.Element(n1).GetPositionBoundary(i);
	  int n2 = mesh_leaf.Boundary(i).numElement(1);
	  int p1 = Epart(n1);
	  int p2 = Epart(n2);
	  int ref = mesh_leaf.BoundaryRef(i).GetReference();
	  if ((ref >= 0) && (p1 != p2)
	      && (mesh_leaf.GetBoundaryCondition(ref) != BoundaryConditionEnum::LINE_NEIGHBOR))
	    {
	      int ref2 = ref_to_change(ref);
	      if (ref2 == 0)
		{
		  ref2 = mesh_leaf.GetNewReference();
		  ref_to_change(ref) = ref2;
		  reference_original_neighbor(ref2) = ref;
		  mesh_leaf.Type_curve(ref2) = mesh_leaf.Type_curve(ref);
		  mesh_leaf.Body_curve(ref2) = mesh_leaf.Body_curve(ref);
		  mesh_leaf.Cond_curve(ref2) = BoundaryConditionEnum::LINE_NEIGHBOR;
		  mesh_leaf.Param_curve(ref2) = mesh_leaf.Param_curve(ref);
		  mesh_leaf.Param_known(ref2) = mesh_leaf.Param_known(ref);
		}

	      mesh_leaf.BoundaryRef(i).SetReference(ref2);
	      mesh_leaf.Boundary(i).SetReference(ref2);

	      int nf1 = sub_mesh(p1).Element(NumLoc(n1)).numBoundary(num_loc);
	      sub_mesh(p1).Boundary(nf1).SetReference(ref2);
	      sub_mesh(p1).BoundaryRef(nf1).SetReference(ref2);
	    }
	}

    for (int i = 0; i < sub_mesh.GetM(); i++)
      {
	sub_mesh(i).Type_curve = mesh_leaf.Type_curve; sub_mesh(i).Cond_curve = mesh_leaf.Cond_curve;
	sub_mesh(i).Body_curve = mesh_leaf.Body_curve; sub_mesh(i).Param_curve = mesh_leaf.Param_curve;
	sub_mesh(i).Param_known = mesh_leaf.Param_known;
	sub_mesh(i).reference_original_neighbor = mesh_leaf.reference_original_neighbor;
      }
  }
  

  //! writes the mesh, reference of each element is set to the order of the element
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::WriteOrder(const string & file_name,
					   const MeshNumbering<Dimension, T>& mesh_num)
  {
    Mesh<Dimension, T>& mesh = this->GetLeafClass();
    Vector<int> color_elt(mesh.GetNbElt());
    for (int i = 0; i < mesh.GetNbElt(); i++)
      color_elt(i) = mesh_num.GetOrderElement(i);
    
    this->WriteColor(file_name, color_elt);
  }


  //! writes the mesh, reference of each element is set to the color given in argument
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>::WriteColor(const string & file_name,
					   const Vector<int>& color_elt)
  {
    Mesh<Dimension, T>& mesh = this->GetLeafClass();
    IVect RefElt(mesh.GetNbElt()), RefFace(mesh.GetNbBoundaryRef());
    // changing references of elements and faces with the order of approximation
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        RefElt(i) = mesh.Element(i).GetReference();
        mesh.Element(i).SetReference(color_elt(i));
      }

    /*
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        RefFace(i) = mesh.BoundaryRef(i).GetReference();
        mesh.BoundaryRef(i).SetReference(mesh_num.GetOrderQuadrature(i));
      }
    */
    
    // writing mesh with these references
    mesh.Write(file_name);
    
    // old references are put back
    for (int i = 0; i < mesh.GetNbElt(); i++)
      mesh.Element(i).SetReference(RefElt(i));
      
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      mesh.BoundaryRef(i).SetReference(RefFace(i));
    
  }
  
  
  //! computes the bounding box of element i
  template<class Dimension, class T>
  void Mesh_Base<Dimension, T>
  ::GetBoundingBox(int i, const VectR_N& s, const SetPoints<Dimension>& PointsElem,
		   TinyVector<R_N, 2>& enveloppe) const
  {
    enveloppe(0) = s(0); enveloppe(1) = s(0);
    R_N pointA;
    if (!this->IsElementAffine(i))
      for (int j = 0; j < PointsElem.GetNbPointsNodal(); j++)
	{
	  pointA = PointsElem.GetPointNodal(j);
	  UpdateMinimum(pointA, enveloppe(0));
	  UpdateMaximum(pointA, enveloppe(1));
	}
    else
      {
	for (int j = 1; j < s.GetM(); j++)
	  {
	    UpdateMinimum(s(j), enveloppe(0));
	    UpdateMaximum(s(j), enveloppe(1));
	  }
      }
  }
  

  //! modifying arrays associated with nodal points
  template<class Dimension, class T>
  void SurfacicMesh_Base<Dimension, T>::
  SetBoundaryNodal(int i, SetPoints<Dimension>& pts,
		   SetMatrices<Dimension>& mat) const
  {
    pts.ReallocatePointsQuadratureBoundary(PointsNodal.GetN());
    mat.ReallocatePointsQuadratureBoundary(PointsNodal.GetN());
    for (int j = 0; j < PointsNodal.GetN(); j++)
      {
	pts.SetPointQuadratureBoundary(j, PointsNodal(i,j));
	mat.SetPointQuadratureBoundary(j, MatricesNodal(i,j));
	mat.SetNormaleQuadratureBoundary(j, NormaleNodal(i,j));
	mat.SetDsQuadratureBoundary(j, DsjNodal(i,j));
      }
  }

  
  //! allocation of nodal arrays
  template<class Dimension, class T>
  void SurfacicMesh_Base<Dimension, T>::ReallocateNodal(int m, int n)
  {
    PointsNodal.Reallocate(m,n);
    NormaleNodal.Reallocate(m,n);
    DsjNodal.Reallocate(m,n); DsjNodal.Fill(0);
    MatricesNodal.Reallocate(m,n);
  }

  
  WrongMesh::WrongMesh(string function, string comment)  :
    Error("Degenerated mesh", function, comment)
  {
    this->CoutWhat();
    abort();
  }  
  
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_BASE_CXX
#endif


