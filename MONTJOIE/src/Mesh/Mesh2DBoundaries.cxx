#ifndef MONTJOIE_FILE_MESH2D_BOUNDARIES_CXX

namespace Montjoie
{
  
  /**********************
   * PmlRegionParameter *
   **********************/


  //! returns average length of the extern boundary
  /*!
    this method is useful to select automatically number of layers in the PML
   */
  template<class T>
  T PmlRegionParameter<Dimension2, T>
  ::GetAverageLengthOnBoundary(const Mesh<Dimension2>& mesh_leaf) const
  {
    T mean = 0.0; int nb = 0;
    if (this->ref_pml.GetM() > 0)
      {
	Vector<int> ref_cond(mesh_leaf.GetNbReferences()+1);
	ref_cond.Zero();
	for (int i = 0; i < this->ref_pml.GetM(); i++)
	  ref_cond(this->ref_pml(i)) = 1;
	
	for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
	  if (ref_cond(mesh_leaf.BoundaryRef(i).GetReference()) == 1)
	    {
	      nb++;
	      mean += Distance(mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(0)),
			       mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(1)));

	    }
      }
    else
      {
	VectBool vertex_on_boundary(mesh_leaf.GetNbVertices());
	vertex_on_boundary.Fill(false);    
	if (this->pml_add_axis_x ==  this->PML_RADIAL)
	  {
	    // center and radius of the circle
	    T Rmax;
	    TinyVector<T, 2> origin = mesh_leaf.GetCenterRadialPML(Rmax);
	    
	    // circular boundary
	    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
	      if (abs(mesh_leaf.Vertex(i).Distance(origin)-Rmax) <= TinyVector<T, 2>::threshold)
		vertex_on_boundary(i) = true;
	  }
	else
	  {
	    // rectangular boundary
	    Real_wp xmin = mesh_leaf.GetXmin(); Real_wp xmax = mesh_leaf.GetXmax();
	    Real_wp ymin = mesh_leaf.GetYmin(); Real_wp ymax = mesh_leaf.GetYmax();
	    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
	      {
		if (abs(mesh_leaf.Vertex(i)(0) - xmin) < TinyVector<T, 2>::threshold)
		  vertex_on_boundary(i) = true;
		
		if (abs(mesh_leaf.Vertex(i)(0) - xmax) < TinyVector<T, 2>::threshold)
		  vertex_on_boundary(i) = true;
		
		if (abs(mesh_leaf.Vertex(i)(1) - ymin) < TinyVector<T, 2>::threshold)
		  vertex_on_boundary(i) = true;
		
		if (abs(mesh_leaf.Vertex(i)(1) - ymax) < TinyVector<T, 2>::threshold)
		  vertex_on_boundary(i) = true;
	      }
	  }
	
	// boucle sur toutes les aretes de ref
	for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
	  {
	    // arete sur le bord exterieur ?
	    bool edge_on_boundary = true;
	    for (int j = 0; j < 2; j++)
	      if (!vertex_on_boundary(mesh_leaf.BoundaryRef(i).numVertex(j)))
		edge_on_boundary = false;
	    
	    if (edge_on_boundary)
	      {
		nb++;
		mean += Distance(mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(0)),
				 mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(1)));
	      }
	  }
      }
    
    if (nb == 0)
      mean = 0.0;
    else
      mean = mean/nb;
    
    return mean;
  }


  //! PML layers are added around the mesh
  /*!
    the method assumes that the extern boundary is rectangular (or circular)
    \param[in] nb_div number of layers
   */
  template<class T>
  void PmlRegionParameter<Dimension2, T>::AddPMLElements(int nb_div, int num_pml,
							 Mesh<Dimension2>& mesh_leaf)
  {
    int nb_elt = mesh_leaf.GetNbElt();
    
    // the mesh is extruded
    T delta = this->thicknessPML;
    
    if (this->ref_pml.GetM() > 0)
      {
	// case where the PML is extruded from a set of references
	int ref = this->ref_pml(0);
	Vector<int> ref_cond(mesh_leaf.GetNbReferences()+1);
	ref_cond.Zero();
	for (int i = 0; i < this->ref_pml.GetM(); i++)
	  ref_cond(this->ref_pml(i)) = 1;
	
	if (mesh_leaf.GetCurveType(ref) == mesh_leaf.CURVE_CIRCLE)
	  {
	    Real_wp x0 = mesh_leaf.GetCurveParameter(ref)(0);
	    Real_wp y0 = mesh_leaf.GetCurveParameter(ref)(1);
	    this->radiusPML = mesh_leaf.GetCurveParameter(ref)(2);
	    this->originPML.Init(x0, y0);
	    
	    //DISP(nb_div); DISP(this->radiusPML); DISP(this->originPML); DISP(delta);
	    ExtrudePMLLayer(mesh_leaf, this->originPML, this->radiusPML,
			    ref_cond, 1, delta, nb_div, true);
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
      }

    if (this->pml_add_axis_x ==  this->PML_RADIAL)
      {
        // center of the circle
        T Rmax;
        TinyVector<T, 2> origin = mesh_leaf.GetCenterRadialPML(Rmax);
        this->radiusPML = Rmax;
        this->originPML = origin;

	Vector<int> ref_cond;
        ExtrudePMLLayer(mesh_leaf, origin, Rmax,
			ref_cond, -1, delta, nb_div, true);
      }
    
    Real_wp xmin = mesh_leaf.GetXmin(), xmax = mesh_leaf.GetXmax();
    Real_wp ymin = mesh_leaf.GetYmin(), ymax = mesh_leaf.GetYmax();

    if ((this->pml_add_axis_x == this->PML_BOTH_SIDES)||
        (this->pml_add_axis_x == this->PML_NEGATIVE_SIDE))
      mesh_leaf.ExtrudeCoordinate(0, nb_div, mesh_leaf.GetXmin(), -delta);
    
    if ((this->pml_add_axis_x == this->PML_BOTH_SIDES)
        ||(this->pml_add_axis_x == this->PML_POSITIVE_SIDE))
      mesh_leaf.ExtrudeCoordinate(0, nb_div, mesh_leaf.GetXmax(), delta);

    if ((this->pml_add_axis_y == this->PML_BOTH_SIDES)
        ||(this->pml_add_axis_y == this->PML_NEGATIVE_SIDE))
      mesh_leaf.ExtrudeCoordinate(1, nb_div, mesh_leaf.GetYmin(), -delta);
    
    if ((this->pml_add_axis_y == this->PML_BOTH_SIDES)
        ||(this->pml_add_axis_y == this->PML_POSITIVE_SIDE))
      mesh_leaf.ExtrudeCoordinate(1, nb_div, mesh_leaf.GetYmax(), delta);
    
    for (int i = nb_elt; i < mesh_leaf.GetNbElt(); i++)
      {
        int type_pml = 0;
        if (this->pml_add_axis_x !=  this->PML_RADIAL)
          {
            int nb_vert = mesh_leaf.Element(i).GetNbVertices();
            R2 center;
            for (int j = 0; j < nb_vert; j++)
              center += mesh_leaf.Vertex(mesh_leaf.Element(i).numVertex(j));

            center *= Real_wp(1) / nb_vert;
            if (center(0) < xmin)
              type_pml += PML_Xm;

            if (center(0) > xmax)
              type_pml += PML_Xp;

            if (center(1) < ymin)
              type_pml += PML_Ym;

            if (center(1) > ymax)
              type_pml += PML_Yp;
          }
        
        mesh_leaf.Element(i).SetPML(num_pml, type_pml);    
      }

  }


  //! Returns the number of parameters
  template<class T>
  int PmlRegionParameter<Dimension2, T>::GetNbParameters() const
  {
    return 3;
  }
  
    
  //! Fills parameters associated with the current region
  template<class T>
  void PmlRegionParameter<Dimension2, T>
  ::FillParameters(VectReal_wp& all_param, int& nb_param) const
  {
    all_param(nb_param) = this->GetRadiusPML(); nb_param++;
    all_param(nb_param) = this->GetOriginRadialPML()(0); nb_param++;
    all_param(nb_param) = this->GetOriginRadialPML()(1); nb_param++;
  }

  
  //! Sets the current region with parameters contained in all_param
  template<class T>
  void PmlRegionParameter<Dimension2, T>
  ::SetRegion(const VectReal_wp& all_param, int& nb_param)
  {
    Real_wp radius_pml = all_param(nb_param); nb_param++;
    this->SetRadiusPML(radius_pml);
    Real_wp x0_pml = all_param(nb_param), y0_pml = all_param(nb_param+1); nb_param += 2;
    this->SetOriginRadialPML(R2(x0_pml, y0_pml));
  }
  

  /******************
   * MeshBoundaries *
   ******************/
  

  //! returns the size of memory used to store the object
  template<class T>
  size_t MeshBoundaries<Dimension2, T>::GetMemorySize() const
  {
    size_t taille = MeshBoundaries_Base<Dimension2, T>::GetMemorySize();
    taille += PointsEdgeRef.GetMemorySize();
    taille += value_parameter_boundary.GetMemorySize();
    taille += integer_parameter_boundary.GetMemorySize();
    taille += interval_reference.GetMemorySize();
    return taille;
  }
  
  
  //! returns the following vertex after n1, n2, and belonging to a referenced curve
  /*!
    \param[in] n1 first point
    \param[in] n2 second point
    \param[in] ref reference of edge to be found
    \param[out] e edge number found
   */
  template<class T>
  int MeshBoundaries<Dimension2, T>
  ::FindFollowingVertex(int n1, int n2,
                        int ref, const IVect& ref_cond, int& e) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    e = mesh_leaf.head_minv(n2);
    const Vector<Edge<Dimension2> >& edges = mesh_leaf.edges;
    // Loop on all edges who leaves this vertex
    while (e != -1)
      {
	// if we find an edge with reference ref 
	// which pass by another vertex than n1 or n2
	// we returns this vertex
	if (ref_cond(edges(e).GetReference()) == ref) 
	  {
	    if ((edges(e).numVertex(0)!=n1)&&(edges(e).numVertex(0)!=n2))
	      return edges(e).numVertex(0);
	    if ((edges(e).numVertex(1)!=n1)&&(edges(e).numVertex(1)!=n2))
	      return edges(e).numVertex(1);
	  }
	
	e = mesh_leaf.next_edge(e);
      }
    // e is the first edge who comes to the vertex n2
    e = mesh_leaf.head_maxv(n2);
    // Loop on all edges who comes to this vertex
    while (e != -1)
      {
	// if we find an edge with reference ref 
	// which pass by another vertex than n1 or n2
	// we returns this vertex
	if (ref_cond(edges(e).GetReference()) == ref)
	  {
	    if ((edges(e).numVertex(0) != n1)&&(edges(e).numVertex(0) != n2))
	      return edges(e).numVertex(0);
	    if ((edges(e).numVertex(1) != n1)&&(edges(e).numVertex(1) != n2))
	      return edges(e).numVertex(1);
	  }
	e = mesh_leaf.next_edge_arr(e);
      }
    
    return -1;
  }
  
  
  //! returns the following vertex after n1, n2, and belonging to a referenced curve
  /*!
    \param[in] n1 first point
    \param[in] n2 second point
    \param[in] ref reference of edge to be found
    \param[out] e edge number found
   */
  template<class T>
  int MeshBoundaries<Dimension2, T>
  ::FindFollowingVertex(int n1, int n2,
                        int ref, const IVect& ref_cond,
                        const Vector<bool>& used_edge, int& e) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    e = mesh_leaf.head_minv(n2);
    const Vector<Edge<Dimension2> >& edges = mesh_leaf.edges;
    
    // Loop on all edges who leaves this vertex
    while (e != -1)
      {
	// if we find an edge with reference ref 
	// which pass by another vertex than n1 or n2
	// we returns this vertex
	if (!used_edge(e))
	  if (ref_cond(edges(e).GetReference()) == ref) 
	    {
	      if ((edges(e).numVertex(0) != n1)&&(edges(e).numVertex(0) != n2))
		return edges(e).numVertex(0);
	      if ((edges(e).numVertex(1) != n1)&&(edges(e).numVertex(1) != n2))
		return edges(e).numVertex(1);
	    }	
	e = mesh_leaf.next_edge(e);
      }
    
    // e is the first edge who comes to the vertex n2
    e = mesh_leaf.head_maxv(n2);
    // Loop on all edges who comes to this vertex
    while (e!=-1)
      {
	// if we find an edge with reference ref 
	// which pass by another vertex than n1 or n2
	// we returns this vertex
	if (!used_edge(e))
	  if (ref_cond(edges(e).GetReference()) == ref) 
	    {
	      if ((edges(e).numVertex(0) != n1)&&(edges(e).numVertex(0) != n2))
		return edges(e).numVertex(0);
	      if ((edges(e).numVertex(1) != n1)&&(edges(e).numVertex(1) != n2))
		return edges(e).numVertex(1);
	    }
	
	e = mesh_leaf.next_edge_arr(e);
      }
    return -1;
  }
  

  //! returns true if the points M is on the referenced edge num_glob_edge_ref
  /*!
    \param[in] num_glob_edge_ref referenced edge number
    \param[in] M point to test
   */
  template<class T>
  bool MeshBoundaries<Dimension2, T>::
  IsPointOnBoundaryRef(int num_glob_edge_ref, const R_N& M) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    R_N A = mesh_leaf.Vertices(mesh_leaf.edges_ref(num_glob_edge_ref).numVertex(0));
    R_N B = mesh_leaf.Vertices(mesh_leaf.edges_ref(num_glob_edge_ref).numVertex(1));
    bool presence=false;
    T det=(M(0)-A(0))*(B(1)-A(1))-(M(1)-A(1))*(B(0)-A(0));
    //DISP(det);
    T ps = (M(0)-A(0))*(B(0)-A(0))+(M(1)-A(1))*(B(1)-A(1));
    if (abs(det) <= R_N::threshold)
      {
	T t = ps/(A.Distance(B)*A.Distance(B));
	
	if ((t>=0)&&(t<=1))
	  presence = true;
      }
    
    return presence;
  }

        
  //! retrieves all the references of edges surrounding a point
  /*!
    \param[in] num_point vertex number
    \param[out] ref list of references
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GetAllReferencesAroundPoint(int num_point, IVect& ref) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    ref.Clear();
    int e = mesh_leaf.head_minv(num_point);
    // Loop on all edges who leaves this vertex
    while (e != -1)
      {
	ref.PushBack(mesh_leaf.edges(e).GetReference()); 
	e = mesh_leaf.next_edge(e);
      }
    
    e = mesh_leaf.head_maxv(num_point);
    // Loop on all edges who comes to this vertex
    while (e != -1)
      {
	ref.PushBack(mesh_leaf.edges(e).GetReference()); 
	e = mesh_leaf.next_edge_arr(e);
      }
    
    int n = ref.GetM();
    // duplicate entries are removed
    RemoveDuplicate(n, ref);
  }
  
  
  //! gets a point on a referenced edge whose reference is different from ref
  template<class T>
  int MeshBoundaries<Dimension2, T>::GetPointWithOtherReference(int num_point, int ref,
								R_N& pt, int& other_ref) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    int e = mesh_leaf.head_minv(num_point); other_ref = 0;    
    // Loop on all edges who leaves this vertex
    while (e != -1)
      {
	other_ref = mesh_leaf.edges(e).GetReference();
	if ((other_ref > 0)&&(other_ref != ref))
	  {
	    if (mesh_leaf.edges(e).numVertex(0) == num_point)
	      {
		pt = mesh_leaf.Vertices(mesh_leaf.edges(e).numVertex(1));
		return mesh_leaf.edges(e).numVertex(1);
	      }
	    else
	      {
		pt = mesh_leaf.Vertices(mesh_leaf.edges(e).numVertex(0));
		return mesh_leaf.edges(e).numVertex(0);
	      }
	  }
	e = mesh_leaf.next_edge(e);
      }
    
    e = mesh_leaf.head_maxv(num_point);
    // Loop on all edges who comes to this vertex
    while (e != -1)
      {
	other_ref = mesh_leaf.edges(e).GetReference();
	if ((other_ref > 0)&&(other_ref != ref))
	  {
	    if (mesh_leaf.edges(e).numVertex(0) == num_point)
	      {
		pt = mesh_leaf.Vertices(mesh_leaf.edges(e).numVertex(1));
		return mesh_leaf.edges(e).numVertex(1);
	      }
	    else
	      {
		pt = mesh_leaf.Vertices(mesh_leaf.edges(e).numVertex(0));
		return mesh_leaf.edges(e).numVertex(0);
	      }
	  }
	e = mesh_leaf.next_edge_arr(e);
      }
    return -1;
  }
    
  
  /**************************************
   * Treatment of referenced boundaries *
   **************************************/
  
  
  //! removes referenced edges with a same reference
  template<class T>
  void MeshBoundaries<Dimension2, T>::RemoveReference(int ref)
  {
    Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    int nb_edges_ref = mesh_leaf.GetNbBoundaryRef();
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (mesh_leaf.edges_ref(i).GetReference() == ref)
	{
	  mesh_leaf.Boundary(i).SetReference(0);
	  nb_edges_ref--;
	}
    
    if (nb_edges_ref < mesh_leaf.GetNbBoundaryRef())
      {
	if (nb_edges_ref == 0)
	  {
	    mesh_leaf.edges_ref.Clear();
	    PointsEdgeRef.Clear();
	  }
	else
	  {
	    // storing previous informations about referenced edges
	    Vector<Edge<Dimension2> > old_edges;
	    old_edges = mesh_leaf.edges_ref;
	    Matrix<R2> OldPoints = PointsEdgeRef;
	    Vector<R_N> old_value = value_parameter_boundary;
	    Vector<TinyVector<int, 2> > old_integer = integer_parameter_boundary;
	    
	    // we change the size
	    mesh_leaf.ResizeBoundariesRef(nb_edges_ref);
	    
	    int nb = 0;
	    // we copy the old edges in the new arrays
	    for (int i = 0; i < old_edges.GetM(); i++)
	      if (old_edges(i).GetReference() != ref)
		{
		  for (int k = 0; k < PointsEdgeRef.GetN(); k++)
		    PointsEdgeRef(nb, k) = OldPoints(i, k); 
		  
		  value_parameter_boundary(nb) = old_value(i);
		  integer_parameter_boundary(nb) = old_integer(i);
		  mesh_leaf.edges_ref(nb++) = old_edges(i);
		}
	  }

	mesh_leaf.FindConnectivity();
      }
  }
  
  
  //! permutes numbering of elements so that periodic edges have the same orientation
  template<class T>
  void MeshBoundaries<Dimension2, T>::SortPeriodicBoundaries()
  {
    Mesh<Dimension2, T>& mesh = this->GetLeafClass();
    
    if (this->periodicity_references.GetM() <= 0)
      return;
    
    bool cyclic_domain = false;
    for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
      if (this->type_coord_periodicity(num_per) == BoundaryConditionEnum::PERIODIC_THETA)
        cyclic_domain = true;
    
    int nb_vertices = mesh.GetNbVertices();
    IVect inv_perm(nb_vertices);
    inv_perm.Fill();    
    if (cyclic_domain)
      {
	R2 coeffDirect1, coeffDirect2, ptC;
	int nb_vertices_ref1(0), nb_vertices_ref2(0);
	int nb_vertices_ref3 = 0;
	IVect TypeVertex(nb_vertices);
	TypeVertex.Fill(-1);
	
	for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
	  {
	    int ref1 = this->periodicity_references(num_per)(0);
	    int ref2 = this->periodicity_references(num_per)(1);
	    
	    // TypeVertex(i) = 1 if vertex i belongs to ref1
	    // TypeVertex(i) = 2 if vertex i belongs to ref2
	    // TypeVertex(i) = 3 if vertex i belongs to both boundaries
	    // we get all the vertices of reference 1 and 2
	    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	      {
		if (mesh.BoundaryRef(i).GetReference() == ref1)
		  {
		    int n1 = mesh.BoundaryRef(i).numVertex(0);
		    int n2 = mesh.BoundaryRef(i).numVertex(1);
		    if (TypeVertex(n1) == 2)
		      TypeVertex(n1) = 3;
		    else if (TypeVertex(n1) == -1)
		      TypeVertex(n1) = 1;

		    if (TypeVertex(n2) == 2)
		      TypeVertex(n2) = 3;
		    else if (TypeVertex(n2) == -1)
		      TypeVertex(n2) = 1;
		  }
		
		if (mesh.BoundaryRef(i).GetReference() == ref2)
		  {
		    int n1 = mesh.BoundaryRef(i).numVertex(0);
		    int n2 = mesh.BoundaryRef(i).numVertex(1);
		    if (TypeVertex(n1) == 1)
		      TypeVertex(n1) = 3;
		    else if (TypeVertex(n1) == -1)
		      TypeVertex(n1) = 2;
		    
		    if (TypeVertex(n2) == 1)
		      TypeVertex(n2) = 3;
		    else if (TypeVertex(n2) == -1)
		      TypeVertex(n2) = 2;
		  }
	      }
          }
        
        // number of vertices on each boundary	    	    
        R2 ptA_ref1, ptB_ref1, ptA_ref2, ptB_ref2;
        for (int i = 0; i < nb_vertices; i++)
          {
            if (TypeVertex(i) == 1)
              {
                if (nb_vertices_ref1 == 0)
                  {
                    ptA_ref1 = mesh.Vertex(i);
                    ptB_ref1 = mesh.Vertex(i);
                  }
                else
                  {
                    if (mesh.Vertex(i) < ptA_ref1)
                      ptA_ref1 = mesh.Vertex(i);
                    
                    if (mesh.Vertex(i) > ptB_ref1)
                      ptB_ref1 = mesh.Vertex(i);
                  }
                
                nb_vertices_ref1++;
              }
            else if (TypeVertex(i) == 2)
              {
                if (nb_vertices_ref2 == 0)
                  {
                    ptA_ref2 = mesh.Vertex(i);
                    ptB_ref2 = mesh.Vertex(i);
                  }
                else
                  {
                    if (mesh.Vertex(i) < ptA_ref2)
                      ptA_ref2 = mesh.Vertex(i);
                    
                    if (mesh.Vertex(i) > ptB_ref2)
                      ptB_ref2 = mesh.Vertex(i);
                  }
                
                nb_vertices_ref2++;
              }
            else if (TypeVertex(i) == 3)
              {
                ptC = mesh.Vertex(i);
                nb_vertices_ref3++;
              }
          }

        // coeffDirect1 : slope of the line corresponding to the boundary ref1
        coeffDirect1(0) =  ptB_ref1(0) - ptA_ref1(0);
        coeffDirect1(1) =  ptB_ref1(1) - ptA_ref1(1);
	
        // coeffDirect2 : slope of the line corresponding to the boundary ref2
        coeffDirect2(0) =  ptB_ref2(0) - ptA_ref2(0);
        coeffDirect2(1) =  ptB_ref2(1) - ptA_ref2(1);
        
        if (nb_vertices_ref3  == 0)
          {	
            Real_wp k1, k2;
            // we retrieve the center of circle
            // method : C(0) = sommetref1(0) + k1 coeffDirect1(0)
            //          C(1) = sommetref1(1) + k1 coeffDirect1(1)
            // 	     C(0) = sommetref2(0) + k2 coeffDirect2(0)
            //       C(1) = sommetref2(1) + k2 coeffDirect2(1)
            // --> il faut distinguer plusieurs cas pour ne pas diviser par 0
            // -- > on calcule k1 et k2 puis le centre du cercle
            if (  (abs(coeffDirect1(0)) <= 1e-10) &&  (abs(coeffDirect2(0)) > 1e-10) )
              {
                ptC(0) = ptA_ref1(0);
                k2 = (ptC(0) - ptA_ref2(0))/ coeffDirect2(0);
                ptC(1) =  ptA_ref2(1) + k2 * coeffDirect2(1);
              }
            else if( (abs(coeffDirect2(0)) <= 1e-10) && (abs(coeffDirect1(0)) >1e-10))
              {
                ptC(0) =  ptA_ref2(0);
                k1 = (ptC(0) - ptA_ref1(0))/ coeffDirect1(0);
                ptC(1) =  ptA_ref1(1) + k1 * coeffDirect1(1);
              }
            else if( (abs(coeffDirect2(0)) <= 1e-10) && (abs(coeffDirect1(0)) < 1e-10))
              {
                k1 = 0 ;
                k2 = 0 ;
                ptC(1) =  0;
                ptC(0) =  0;
              }
            else if( (abs(coeffDirect2(1)) <= 1e-10) && (abs(coeffDirect1(1)) < 1e-10))
              {
                k1 = 0 ;
                k2 = 0 ;
                ptC(1) =  0;
                ptC(0) =  0;
              }
            else if( (abs(coeffDirect2(1)) <= 1e-10) && (abs(coeffDirect1(1)) > 1e-10 )) 
              {
                ptC(1) =  ptA_ref2(1);
                k1 = (ptC(1) - ptA_ref1(1))/ coeffDirect1(1);
                ptC(0) =  ptA_ref1(0) + k1 * coeffDirect1(0);
              }
            else if( (abs(coeffDirect1(1)) <= 1e-10) && (abs(coeffDirect2(1)) > 1e-10 ) )
              {
                ptC(1) = ptA_ref1(1);
                k2 = (ptC(1) - ptA_ref2(1))/ coeffDirect2(1);
                ptC(0) =  ptA_ref2(0) + k2  *coeffDirect2(0);
		
              }
            else
              {
                k2 = ( coeffDirect1(1) * ptA_ref1(0)  - coeffDirect1(0)*ptA_ref1(1)
                       - coeffDirect1(1)* ptA_ref2(0) + coeffDirect1(0)* ptA_ref2(1))/
                  ( coeffDirect1(1)*coeffDirect2(0) - coeffDirect1(0)*coeffDirect2(1));
                
                k1 = (ptA_ref2(0) - ptA_ref1(0)  + k2* coeffDirect2(0))/coeffDirect1(0);
		
                ptC(0) = ptA_ref1(0) + k1 * coeffDirect1(0);
                ptC(1) = ptA_ref1(1) + k1 * coeffDirect1(1);
              }
          }
        
        this->periodic_center_polar = ptC;
        
	Mlt(1.0/Norm2(coeffDirect1), coeffDirect1);
	Mlt(1.0/Norm2(coeffDirect2), coeffDirect2);
	
	// on repere l angle oriente theta 
	// signe du cosinus
	Real_wp signecos = DotProd(coeffDirect1, coeffDirect2);
	
	// theta compris entre (-pi et pi)
       	Real_wp theta = asin( coeffDirect1(0) * coeffDirect2(1)
			      - coeffDirect1(1) * coeffDirect2(0)) ;
	// si cos < 0 , l angle est en fait = pi -theta
	if (signecos < 0 ) 
	  theta = pi_wp - theta ;

	this->periodic_alpha_polar = theta;
        
	// we are retrieving vertices of Ref1 and Ref2
	// (vertices of Ref2 are rotated)
	IVect permutation(nb_vertices);
	permutation.Fill();
	int num_vertex3 = -1;
	Vector<Real_wp> RadiusRef1(nb_vertices_ref1), RadiusRef2(nb_vertices_ref2);
	IVect permut1(nb_vertices_ref1), permut2(nb_vertices_ref2);
        nb_vertices_ref1 = 0;
        nb_vertices_ref2 = 0;
        nb_vertices_ref3 = 0;
	for (int i = 0; i < nb_vertices; i++)
	  {
	    if (TypeVertex(i) == 1)
	      {
		RadiusRef1(nb_vertices_ref1) = ptC.Distance(mesh.Vertex(i));
		permut1(nb_vertices_ref1) = i;
		nb_vertices_ref1++;
	      }
	    else if (TypeVertex(i) == 2)
	      {
		RadiusRef2(nb_vertices_ref2) = ptC.Distance(mesh.Vertex(i));
		permut2(nb_vertices_ref2) = i;
		nb_vertices_ref2++;
	      }
	    else if (TypeVertex(i) == 3)
	      {
		num_vertex3 = i;
		nb_vertices_ref3++;
	      }
	  }
	
	// sorting points by ascending radius
	Sort(nb_vertices_ref1, RadiusRef1, permut1);
	Sort(nb_vertices_ref2, RadiusRef2, permut2);
	
	// modification of permutation	
	int nb = 0;
	if (nb_vertices_ref3 > 0)
	  permutation(nb++) = num_vertex3;
	
	for (int i = 0; i < nb_vertices_ref1; i++)
	  permutation(nb++) = permut1(i);
	
	for (int i = 0; i < nb_vertices_ref2; i++)
	  permutation(nb++) = permut2(i);
	
	for (int i = 0; i < nb_vertices; i++)
	  if (TypeVertex(i) < 0)
	    permutation(nb++) = i;
	
        VectR2 OldVert = mesh.Vertex();
	for (int i = 0; i < nb_vertices; i++)
	  {
            inv_perm(permutation(i)) = i;
            mesh.Vertex(i) = OldVert(permutation(i));
          }
        
      }	    
    else
      {
	// if cartesian periodicity, we sort all the vertices of periodic boundaries
	Vector<bool> VertexToSort(nb_vertices);
	IVect permutation(nb_vertices);
	VertexToSort.Fill(false);
	permutation.Fill();
	
	IVect ref_cond(mesh.GetNbReferences()+1);
	ref_cond.Fill(0);
	for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
	  {
	    int ref1 = this->periodicity_references(num_per)(0);
	    int ref2 = this->periodicity_references(num_per)(1);
	    ref_cond(ref1) = 1;
	    ref_cond(ref2) = 1;
	  }
	
	for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	  {
	    int ref = mesh.BoundaryRef(i).GetReference();
	    int n1 = mesh.BoundaryRef(i).numVertex(0);
	    int n2 = mesh.BoundaryRef(i).numVertex(1);
	    if (ref_cond(ref) == 1)
	      {
		VertexToSort(n1) = true;
		VertexToSort(n2) = true;
	      }
	  }
	
	int nb = 0;
	for (int i = 0; i < nb_vertices; i++)
	  if (VertexToSort(i))
	    nb++;
	
	Vector<R_N> VertBound(nb);
	IVect permut(nb);
	nb = 0;
	for (int i = 0; i < nb_vertices; i++)
	  if (VertexToSort(i))
	    {
	      VertBound(nb) = mesh.Vertex(i);
	      permut(nb) = i;
	      nb++;
	    }
	
	Sort(nb, VertBound, permut);
	
	// injecting permut into permutation
	nb = 0;
	for (int i = 0; i < nb_vertices; i++)
	  if (VertexToSort(i))
	    {
	      permutation(i) = permut(nb);
	      mesh.Vertex(i) = VertBound(nb);
	      nb++;
	    }
	
	// inverse of permutation
	for (int i = 0; i < nb_vertices; i++)
	  inv_perm(permutation(i)) = i;
      }
    
    // changing edge numbers
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	int ref = mesh.BoundaryRef(i).GetReference();
	n1 = inv_perm(n1);
	n2 = inv_perm(n2);
	mesh.BoundaryRef(i).Init(n1, n2, ref);
      }
    
    // changing element numbers
    IVect num;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	int nb_vert = mesh.Element(i).GetNbVertices();
	int ref = mesh.Element(i).GetReference();
	num.Reallocate(nb_vert);
	for (int j = 0; j < nb_vert; j++)
	  num(j) = inv_perm(mesh.Element(i).numVertex(j));
	
	mesh.Element(i).Init(num, ref);
      }    
  }
  
  
  //! permutes extremity numbers of referenced edges so that first extremity number is lowest
  /*!
    You should ProjectPointsOnCurves after this method, otherwise value_parameter_boundary
    and integer_parameter_boundary won't be correctly initialized
  */
  template<class T>
  void MeshBoundaries<Dimension2, T>::SortBoundariesRef()
  {
    Mesh<Dimension2, T>& mesh = this->GetLeafClass();
    
    // sorting vertex numbers for referenced edges
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	int ref = mesh.BoundaryRef(i).GetReference();
	if (n1 > n2)
	  {
	    mesh.BoundaryRef(i).Init(n2, n1, ref);
	    if (PointsEdgeRef.GetN() > 0)
	      {
		int r = this->lob_curve.GetOrder();
		Vector<R_N> OldPoints(r-1);
		for (int k = 0; k < r-1; k++)
		  OldPoints(k) = PointsEdgeRef(i, k);
		
		for (int k = 0; k < r-1; k++)
		  PointsEdgeRef(i, k) = OldPoints(r-2-k); 
	      }
	  }
      }
  }
  

  //! checks that referenced edges correspond to actual edges of the mesh
  template<class T>
  void MeshBoundaries<Dimension2, T>::ConstructCrossReferenceBoundaryRef(bool check_mesh)
  { 
    Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    
    // we put in the array num_element of each edge :
    // num_element(0) : lower reference domain
    // num_element(1) : higher reference domain
    // This sort is done, in order to have always an inward normal at the interface
    // if we  compute normal with the first element
    // If reference are equal, we consider that num_element(0)
    // is the element with the higher number
    for (int i = 0; i < mesh_leaf.GetNbBoundary(); i++)
      if ((mesh_leaf.GetEdge(i).GetNbElements() > 1)&&(mesh_leaf.GetEdge(i).numElement(1) >= 0))
	{
	  if (mesh_leaf.Element(mesh_leaf.GetEdge(i).numElement(0)).GetReference()
	      > mesh_leaf.Element(mesh_leaf.GetEdge(i).numElement(1)).GetReference())
	    mesh_leaf.GetEdge(i).SetElement(mesh_leaf.GetEdge(i).numElement(1),
                                            mesh_leaf.GetEdge(i).numElement(0));
	  else if (mesh_leaf.Element(mesh_leaf.GetEdge(i).numElement(0)).GetReference()
		   == mesh_leaf.Element(mesh_leaf.GetEdge(i).numElement(1)).GetReference())
	    {
	      if (mesh_leaf.GetEdge(i).numElement(0) < mesh_leaf.GetEdge(i).numElement(1))
		mesh_leaf.GetEdge(i).SetElement(mesh_leaf.GetEdge(i).numElement(1),
                                                mesh_leaf.GetEdge(i).numElement(0));
	    }
	}
    
    // loop on all the referenced edges of the mesh
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	// extremities of the referenced edge i
	int iglob = mesh_leaf.BoundaryRef(i).numVertex(0);
	int jglob = mesh_leaf.BoundaryRef(i).numVertex(1);
        if (iglob > jglob)
          {
            cout<< " Call SortBoundariesRef in order to have sorted numbers " <<endl;
            abort();
          }
	
	// loop on all global edges, which leave vertex iglob
	bool edge_found = false;
	for (int e = mesh_leaf.head_minv(iglob); e != -1; e = mesh_leaf.next_edge(e))
	  {
	    if (mesh_leaf.edges(e).numVertex(1) == jglob)
	      {
		// we have found the edge referenced
		// we update edges_ref(i) and edges(e)
		edge_found = true;
		mesh_leaf.edges(e).Init(iglob, jglob, mesh_leaf.BoundaryRef(i).GetReference());
		if (i != e)
		  {
		    cout << "A problem occured during FindConnectivity" << endl;
		    abort();
		  }
		
		mesh_leaf.BoundaryRef(i).SetElement(mesh_leaf.GetEdge(e).numElement(0),
                                                    mesh_leaf.GetEdge(e).numElement(1));
		break;
		// edges_ref are equal to edges
		// so that the extremities are in the same way
	      }
	  }
	
	if ((!edge_found)&&(check_mesh))
	  {
	    // we didn't find the referenced edge
	    // the mesh has a problem, we stop program
	    cout<<"Referenced edges are wrong. Check your mesh "<<endl;
	    mesh_leaf.Write("test.mesh");
	    abort();
	  }
      }
  }
  
  
  //! this method adds eventual edges on the boundary
  //! that are not already stored in the array edges_ref_
  template<class T>
  void MeshBoundaries<Dimension2, T>::AddBoundaryEdges()
  {
    Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    
    // edges on the boundary have only one element (triangle or quadrangle)
    // edges inside have two elements !
    int nb_old_edges_ref = mesh_leaf.GetNbBoundaryRef();
    int nb_isolated_edges = 0, nb_unref_edges =0;
    for (int i = nb_old_edges_ref; i < mesh_leaf.GetNbBoundary(); i++)
      {
	// edge with a reference in array edges and not found in the array edges_ref
	// it is an unreferenced edge
	if (mesh_leaf.Boundary(i).GetReference() != 0)
	  nb_unref_edges++;
	
	// edge on a boundary and not in array edges_ref, it is an edge on the boundary
	else if (mesh_leaf.Boundary(i).GetNbElements() == 1)
	  nb_isolated_edges++;
      }
    
    // if we find isolated or unreferenced edges, we add them to edges_ref_
    if (nb_isolated_edges + nb_unref_edges > 0)
      {
	int nb_old = nb_old_edges_ref;
	mesh_leaf.ResizeBoundariesRef(nb_old_edges_ref + nb_isolated_edges + nb_unref_edges);
	
	// we add a new boundary condition, we consider that the edge is a straight line 
	// with no boundary condition
	int new_reference = this->GetNewReference();
	
	// we update the array edges_ref_
	for (int i = nb_old; i < mesh_leaf.GetNbBoundary(); i++)
	  {
	    if (mesh_leaf.Boundary(i).GetReference() != 0)
	      {
		mesh_leaf.BoundaryRef(nb_old_edges_ref) = mesh_leaf.Boundary(i);
		nb_old_edges_ref++;
	      }
	    else if (mesh_leaf.Boundary(i).GetNbElements() == 1)
	      {
		mesh_leaf.BoundaryRef(nb_old_edges_ref) = mesh_leaf.Boundary(i);
		mesh_leaf.BoundaryRef(nb_old_edges_ref).SetReference(new_reference);
		nb_old_edges_ref++;
	      }
	  }

	mesh_leaf.FindConnectivity();
      }
  }
  
  
  //! this methods has several aims
  /*!
    detection of boundaries which are not in the referenced edges
    construction of links between referenced edges and global edges
    splitting of boundaries in different curves, (list of successive points)
    detection of corners (an angle between two edges less than pi-alpha )
    projection of points on circles/ellipses/splines etc 
  */
  template <class T>
  void MeshBoundaries<Dimension2, T>::ProjectPointsOnCurves()
  {
    Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
        
    // we find parameters of curves
    this->FindParametersCurve();
    
    // We search points on the curves and project them exactly on the curves
    IVect VerticesToBeProjected(mesh_leaf.GetNbVertices());
    VerticesToBeProjected.Fill(-1);
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	int ref = mesh_leaf.BoundaryRef(i).GetReference();
	if ((ref >= 1)&&(this->GetCurveType(ref) > 0))
	  for (int j = 0; j < 2; j++)
	    {
	      VerticesToBeProjected(mesh_leaf.BoundaryRef(i).numVertex(j)) = ref;
	    }
      }
    
    // projection of vertices on curved edges
    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
      if (VerticesToBeProjected(i) > 0)
	{
	  ProjectToCurve(mesh_leaf.Vertex(i), VerticesToBeProjected(i));
        }
    
    this->ComputeParameterValueOnCurves();
      
    if (this->lob_curve.GetOrder() > 1)
      {
	VectReal_wp xi(this->lob_curve.GetOrder()-1);
	for (int i = 1; i < this->lob_curve.GetOrder(); i++)
	  xi(i-1) = this->lob_curve.Points(i);
	
	if (this->PointsEdgeRef.GetM() != mesh_leaf.GetNbBoundaryRef())
	  this->PointsEdgeRef.Resize(mesh_leaf.GetNbBoundaryRef(), this->lob_curve.GetOrder()-1);
	
	for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
	  {
	    int ref = mesh_leaf.BoundaryRef(i).GetReference();		    
	    if (this->GetCurveType(ref) == CURVE_FILE)
	      {
		// nothing, we consider that in that case,
		// PointsEdgeRef is filled during the reading or the mesh
		// or filled by the user
	      }
	    else if (this->GetCurveType(ref) > 0)
	      {
		// predefined curve, we project the internal points 
		Vector<R_N> points_interm;
		this->GetPointsOnCurve(i, ref, xi, points_interm);
		for (int k = 0; k < xi.GetM(); k++)
                  this->PointsEdgeRef(i, k) = points_interm(k);
                
	      }
	    else 
	      {
		// straight boundary
		R_N ptA = mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(0));
		R_N ptB = mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(1));
		for (int k = 0; k < xi.GetM(); k++)
		  {
		    this->PointsEdgeRef(i, k).Zero();
		    Add(1.0-xi(k), ptA, this->PointsEdgeRef(i, k));
		    Add(xi(k), ptB, this->PointsEdgeRef(i, k));
		  }
	      }
	  }
      }       

    // updating bounding box for the mesh
    R_N vec_u;
    for (int i = 0; i < this->PointsEdgeRef.GetM(); i++)
      for (int j = 0; j < this->PointsEdgeRef.GetN(); j++)
        {
          vec_u = this->PointsEdgeRef(i, j);
          mesh_leaf.xmin_ = min(mesh_leaf.xmin_, vec_u(0));
          mesh_leaf.xmax_ = max(mesh_leaf.xmax_, vec_u(0));
          mesh_leaf.ymin_ = min(mesh_leaf.ymin_, vec_u(1));
          mesh_leaf.ymax_ = max(mesh_leaf.ymax_, vec_u(1));
        }
  }
  
  
  //! affecting a different reference number per connected curve
  template<class T>
  void MeshBoundaries<Dimension2, T>::RedistributeReferences()
  {
    Mesh<Dimension2>& mesh_leaf = this->GetLeafClass();
    
    // first, constructing list of edges for each point
    IVect NbEdgesPerPoint(mesh_leaf.GetNbVertices());
    Vector<TinyVector<int, 4> > ListEdgesPerPoint(mesh_leaf.GetNbVertices());
    //Vector<bool> EdgeUsed(mesh_leaf.GetNbBoundaryRef());
    NbEdgesPerPoint.Fill(0);
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
        int n0 = mesh_leaf.BoundaryRef(i).numVertex(0);
        int n1 = mesh_leaf.BoundaryRef(i).numVertex(1);
        if (NbEdgesPerPoint(n0) < 4)
          {
            ListEdgesPerPoint(n0)(NbEdgesPerPoint(n0)) = i;
            NbEdgesPerPoint(n0)++;
          }
        else
          {
            cout << "impossible" << endl;
            abort();
          }

        if (NbEdgesPerPoint(n1) < 4)
          {
            ListEdgesPerPoint(n1)(NbEdgesPerPoint(n1)) = i;
            NbEdgesPerPoint(n1)++;
          }
        else
          {
            cout << "impossible" << endl;
            abort();
          }
      }
    
    int nb_remaining_face = mesh_leaf.GetNbBoundaryRef();
    VectBool FaceUsed(nb_remaining_face);  FaceUsed.Fill(false);
    IVect ListeFace(nb_remaining_face); ListeFace.Fill(-1);
    int ref = 1;
    // while there is a remaining referenced face
    while (nb_remaining_face > 0)
      {
	// we select the first remaining face
	int nf = -1;
	for (int i = 0; i < FaceUsed.GetM(); i++)
	  if (!FaceUsed(i))
	    {
	      nf = i;
	      break;
	    }
        
	//  ListeFace will contain the list of faces belonging to the same
	// connex surface as face nf
	int nb_new_face = 0;
	IVect new_face(1); new_face(0) = nf; FaceUsed(nf) = true;
	while (new_face.GetM() > 0)
	  {
	    // we add new_face to ListeFace
	    for (int i = 0; i < new_face.GetM(); i++)
	      {
		ListeFace(nb_new_face) = new_face(i);
		nb_new_face++; nb_remaining_face--;
	      }
	    
	    // loop over all the edges of new faces to retrieve surrounding faces
	    IVect num_face;
	    for (int i = 0; i < new_face.GetM(); i++)
	      {
		for (int j = 0; j < 2; j++)
		  {
		    int nv = mesh_leaf.BoundaryRef(new_face(i)).numVertex(j);
		    for (int k = 0; k < NbEdgesPerPoint(nv); k++)
		      {
			int ne = ListEdgesPerPoint(nv)(k);
                        if (!FaceUsed(ne))
                          {
                            num_face.PushBack(ne);
                            FaceUsed(ne) = true;
                          }
		      }
		  }
	      }
	    
	    new_face = num_face;
	  }
	
	// affecting new reference to those faces
	for (int i = 0; i < nb_new_face; i++)
	  mesh_leaf.BoundaryRef(ListeFace(i)).SetReference(ref);
	
	// next reference
	ref++;
      }
  }
      

  /**********************************
   * Treatment of curved boundaries *
   **********************************/
  
  
  //! checks if the parameters given by the user are correct
  template<class T>
  void MeshBoundaries<Dimension2, T>::CheckCurveParameter(int ref)
  {    
    Vector<T>& param = this->Param_curve(ref);
    switch (this->GetCurveType(ref))
      {
      case CURVE_CIRCLE:
        {
	  if (param.GetM() < 3)
	    {
	      cout << "A circle is determined by three parameters (Ox, Oy) and the radius" << endl;
	      cout << "Only " << param.GetM() << " are given " << endl;
	      cout << "Current parameters of reference " << ref << " : " << param << endl;
	      abort();
	    }
	  
	  if (param(2) <= 0)
	    {
	      cout << "The radius of the circle must be strictly positive but is equal to " << param(2) << endl;
	      abort();
	    }
	}
        break;
      case CURVE_ELLIPSE:
        {
	  if (param.GetM() < 4)
	    {
	      cout << "An ellipse is determined by three parameters (Ox, Oy) and the two radii" << endl;
	      cout << "Only " << param.GetM() << " are given " << endl;
	      cout << "Current parameters of reference " << ref << " : " << param << endl;
	      abort();
	    }
	  
	  if ((param(2) <= 0) || (param(3) <= 0))
	    {
	      cout << "The radii of the ellipse must be strictly positive but are equal to " << param(2) << " and " << param(3) << endl;
	      abort();
	    }
	}
        break;
      case CURVE_PEANUT:
        {
	  if (param.GetM() < 5)
	    {
	      cout << "A peanut is determined by five parameters (Ox, Oy) and a,b,c" << endl;
	      cout << "Only " << param.GetM() << " are given " << endl;
	      cout << "Current parameters of reference " << ref << " : " << param << endl;
	      abort();
	    }
	  
	  if ((param(2) <= 0) || (param(3) <= 0) || (param(4) <= 0))
	    {
	      cout << "The parameters a, b, c of the peanut must be strictly positive but are equal to "
		   << param(2) << ", " << param(3) << " and " << param(4) << endl;
	      
	      abort();
	    }
	}
        break;
      }
  }
  
  
  //! reads type of curve
  /*!
    \param[in] ref references of all the curves with this type
    \param[in] description keyword to describe type of curve (CIRCLE, PARABOLE...)
    This methods modifies array Type_curve
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  ReadCurveType(const IVect& ref, const string& description)
  {
    if (!description.compare("CIRCLE"))
      for (int i = 0; i < ref.GetM(); i++)
	this->Type_curve(ref(i)) = this->CURVE_CIRCLE;
    else if (!description.compare("ELLIPSE"))
      for (int i = 0; i < ref.GetM(); i++)
	this->Type_curve(ref(i)) = this->CURVE_ELLIPSE;
    else if (!description.compare("LOCAL_SPLINE"))
      for (int i = 0; i < ref.GetM(); i++)
	this->Type_curve(ref(i)) = this->CURVE_LOCAL_SPLINE;
    else if (!description.compare("SPLINE"))
      for (int i = 0; i < ref.GetM(); i++)
	this->Type_curve(ref(i)) = this->CURVE_SPLINE;
    else if (!description.compare("PEANUT"))
      for (int i = 0; i < ref.GetM(); i++)
	this->Type_curve(ref(i)) = this->CURVE_PEANUT;
    else
      for (int i = 0; i < ref.GetM(); i++)
	this->Type_curve(ref(i)) = 0;
    
  }
  
  
  //! projection of a point on the curve
  /*!
    \param[in,out] pt point to be projected
    \param[in] ref reference of the curved surface
    \param[in] n1 first vertex
    \param[in] n2 second vertex
    ptM is assumed to be located between two vertices n1 and n2
  */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  ProjectToCurve(R_N& pt, int ref) const
  {
    // const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    if ((ref != 0)&&(this->Type_curve(ref)>=1))
      {
	// curved edge
	if (this->Type_curve(ref) == this->CURVE_ELLIPSE)
	  {
	    // center of the ellipse
	    R_N center(this->Param_curve(ref)(0), this->Param_curve(ref)(1));
	    // a and b
	    T a = this->Param_curve(ref)(2), b = this->Param_curve(ref)(3);

	    // first estimate of teta
	    T teta = 0;
	    if (abs(pt(0)-center(0)) < R_N::threshold)
	      {
		teta = 0;
	      }
	    else
	      {
		teta = atan((pt(1)-center(1))*a/((pt(0)-center(0))*b));
		if (pt(0)-center(0) < 0)
		  teta += pi_wp;
	      }
	    
	    // then minimization	    
	    EllipseMinimizer fct(pt(0), pt(1), center, a, b, teta);
	    TinyVector<T, 1> xsol;
	    
	    // optimization to have the closest point to the ellipse
	    SolveLeastSquaresLvm(fct, xsol, 1e-16, -1, 30, 0, 0.1);
	    teta = xsol(0);
	    
	    // final point
	    pt(0) = center(0) + a*cos(teta);
	    pt(1) = center(1) + b*sin(teta);
	  }
	else if (this->Type_curve(ref)==this->CURVE_PEANUT)
	  {
	    // center of the cacahuete
	    R_N center(this->Param_curve(ref)(0), this->Param_curve(ref)(1));
	    // a and b
	    T a = this->Param_curve(ref)(2),
	      b = this->Param_curve(ref)(3), c = this->Param_curve(ref)(4);
	    
	    // first estimate of teta
	    T sin_teta = (pt(1)-center(1))/b;
	    sin_teta = min(Real_wp(1), sin_teta);
	    sin_teta = max(Real_wp(-1), sin_teta);
	    T teta = asin(sin_teta);
	    
	    T cos_teta = sqrt(Real_wp(1)-sin_teta*sin_teta);
	    if (pt(0) < center(0))
	      {
		cos_teta = -cos_teta;
		teta = pi_wp - teta;
	      }
	    
	    // minimization
	    PeanutMinimizer fct(pt(0), pt(1), center, a, b, c, teta);
	    TinyVector<T, 1> xsol;

	    SolveLeastSquaresLvm(fct, xsol, 1e-16, -1, 30, 0, 0.1);
	    teta = xsol(0);
	    cos_teta = cos(teta); sin_teta = sin(teta);
	    
	    // final point
	    pt(0) = center(0) + a*cos_teta*sqrt(1.0+c*sin_teta*sin_teta);
	    pt(1) = center(1) + b*sin_teta;
	  }
	else if (this->Type_curve(ref)==this->CURVE_SPLINE)
	  {
	    // no projection since spline relies on the points of the mesh
	  }
	else if (this->Type_curve(ref)==this->CURVE_LOCAL_SPLINE)
	  {
	    // no projection since spline relies on the points of the mesh
	  }
	else if (this->Type_curve(ref) == this->CURVE_CIRCLE)
	  {
	    // center of the circle
	    R_N center(this->Param_curve(ref)(0), this->Param_curve(ref)(1));
	    // radius of the circle
	    T radius = this->Param_curve(ref)(2);
	    //DISP(center); DISP(radius);
	    R_N diff; diff = pt-center;
	    Mlt(radius/Norm2(diff), diff);
	    Add(diff, center, pt);
	  }
      }
  }
  
  
  template<class T>
  void MeshBoundaries<Dimension2, T>::ComputeParameterValueOnCurves()
  {
    Mesh<Dimension2>& mesh = this->GetLeafClass();

    // marking curved elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      mesh.Element(i).UnsetCurved();
    
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int nb_elem_edge = mesh.BoundaryRef(i).GetNbElements();
	// Type_curve( edges_ref(i).ref ) >= 1 ->  we have a curved element
	if (this->Type_curve(mesh.edges_ref(i).GetReference()) >= 1)
	  for (int k = 0; k < nb_elem_edge; k++)
	    {
              int num_elem = mesh.edges_ref(i).numElement(k);
              mesh.Element(num_elem).SetCurved();
              mesh.affine_element(num_elem) = false;
            }
      }

    R_N ptA, ptB;
    value_parameter_boundary.Reallocate(mesh.GetNbBoundaryRef());
    integer_parameter_boundary.Reallocate(mesh.GetNbBoundaryRef());
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int ref = mesh.BoundaryRef(i).GetReference();
	const VectReal_wp& param = this->Param_curve(ref);
	int n1 = mesh.BoundaryRef(i).numVertex(0);
	int n2 = mesh.BoundaryRef(i).numVertex(1);
	if ((n1 >= 0)&&(n2 >= 0))
	  {
	    ptA = mesh.Vertex(n1); ptB = mesh.Vertex(n2);
	    switch (mesh.GetCurveType(ref))
	      {
	      case CURVE_CIRCLE :
	      case CURVE_ELLIPSE :
		{
		  // computing teta for first point
		  Real_wp x = ptA(0) - param(0);
		  Real_wp y = ptA(1) - param(1);
		  Real_wp r = param(2);
		  Real_wp tetaA(0);
		  if (x >= r)
		    tetaA = 0;
		  else if (x <= -r)
		    tetaA = pi_wp;
		  else
		    tetaA = acos(x/r);
		  
		  if (y < 0)
		    tetaA = -tetaA;
		  
		  // computing teta for second point	      
		  x = ptB(0) - param(0);
		  y = ptB(1) - param(1);
		  Real_wp tetaB(0);
		  if (x >= r)
		    tetaB = 0;
		  else if (x <= -r)
		    tetaB = pi_wp;
		  else
		    tetaB = acos(x/r);
		  
		  if (y < 0)
		    tetaB = -tetaB;
		  
		  // if the angles are too different, we try a translation of 2 pi
		  if (tetaA > tetaB+pi_wp)
		    tetaA -= 2*pi_wp;
		  else if (tetaA < tetaB-pi_wp)
		    tetaA += 2.0*pi_wp;

		  // storing tetaA and tetaB
		  value_parameter_boundary(i)(0) = tetaA;
		  value_parameter_boundary(i)(1) = tetaB;
		}
		break;
	      case CURVE_PEANUT :
		{
		  // computing teta for first point
		  Real_wp x = ptA(0) - param(0);
		  Real_wp y = ptA(1) - param(1);
		  Real_wp b = param(3);
		  Real_wp tetaA(0);
		  if (y >= b)
		    tetaA = 0.5*pi_wp;
		  else if (y <= -b)
		    tetaA = -0.5*pi_wp;
		  else
		    tetaA = asin(y/b);
		  
		  if (x < 0)
		    tetaA = pi_wp-tetaA;
		  
		  // computing teta for second point	      
		  x = ptB(0) - param(0);
		  y = ptB(1) - param(1);
		  Real_wp tetaB(0);
		  if (y >= b)
		    tetaB = 0.5*pi_wp;
		  else if (y <= -b)
		    tetaB = -0.5*pi_wp;
		  else
		    tetaB = asin(y/b);

		  if (x < 0)
		    tetaB = pi_wp-tetaB;
		  
		  // if the angles are too different, we try a translation of 2 pi
		  if (tetaA > tetaB+pi_wp)
		    tetaA -= 2.0*pi_wp;
		  else if (tetaA < tetaB-pi_wp)
		    tetaA += 2.0*pi_wp;

		  // storing tetaA and tetaB
		  value_parameter_boundary(i)(0) = tetaA;
		  value_parameter_boundary(i)(1) = tetaB;
		}
		break;
	      }
	  }
      }
    
    // specific treatment for splines
    IVect ListePoints; IVect IndexVert(mesh.GetNbVertices());
    IVect ref_cond(this->GetNbReferences()); ref_cond.Fill();
    for (int ref = 0; ref < this->Type_curve.GetM(); ref++)
      {
	if (this->Type_curve(ref) == this->CURVE_SPLINE)
	  {
	    IndexVert.Fill(-1);
	    FindParametersSpline(ref, this->Param_curve(ref), ListePoints);
	    // constructing index
	    for (int i = 0; i < this->Param_curve(ref).GetM()/4; i++)
	      IndexVert(ListePoints(i)) = i;
	    
	    // now affecting integer_parameter_boundary
	    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	      if (mesh.BoundaryRef(i).GetReference() == ref)
		for (int k = 0; k < 2; k++)
		  integer_parameter_boundary(i)(k) = 4*IndexVert(mesh.BoundaryRef(i).numVertex(k));
	    
	  }
	else if (this->Type_curve(ref) == this->CURVE_LOCAL_SPLINE)
	  {
	    int iedge = -1, n0, n3;
	    // for each edge, we find the point before and the point after
	    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	      if (mesh.BoundaryRef(i).GetReference() == ref)
		{
		  int n1 = mesh.BoundaryRef(i).numVertex(0);
		  int n2 = mesh.BoundaryRef(i).numVertex(1);
		  n0 = FindFollowingVertex(n2, n1, ref, ref_cond, iedge);
		  n3 = FindFollowingVertex(n1, n2, ref, ref_cond, iedge);
		  integer_parameter_boundary(i)(0) = n0;
		  integer_parameter_boundary(i)(1) = n3;
		}
	  }
      }
  }
  
  
  //! computation of internal points on a curved edge
  /*!
    \param[in] i referenced edge number
    \param[in] ref reference
    \param[in] lambda values of parameter s on interval [0,1]
    \param[out] points computed points Gamma(s) at the given values of parameter s
    Here Gamma(s) denotes the curve on which the edge is
    (Gamma(0) => first extremity of edge i, Gamma(1) => second extremity of edge i)
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>
  ::GetPointsOnCurve(int i, int ref, const VectReal_wp& lambda, Vector<R_N>& points) const
  {
    const Mesh<Dimension2>& mesh = this->GetLeafClass();
    const VectReal_wp& param = this->Param_curve(ref);
    points.Reallocate(lambda.GetM()); R_N point;
    switch (this->GetCurveType(ref))
      {
      case NO_CURVE :
	{
	  R2 ptA = mesh.Vertex(mesh.BoundaryRef(i).numVertex(0));
	  R2 ptB = mesh.Vertex(mesh.BoundaryRef(i).numVertex(1));
	  for (int k = 0; k < lambda.GetM(); k++)
	    {
	      points(k).Zero();
	      Add(1.0-lambda(k), ptA, points(k));
	      Add(lambda(k), ptB, points(k));
	    }
	}
	break;
      case CURVE_CIRCLE :
	{
	  Real_wp teta1 = this->value_parameter_boundary(i)(0);
	  Real_wp teta2 = this->value_parameter_boundary(i)(1);
	  for (int k = 0; k < lambda.GetM(); k++)
	    {
	      Real_wp teta = (1.0-lambda(k))*teta1 + lambda(k)*teta2;
	      
	      points(k)(0) = param(0) + param(2)*cos(teta);
	      points(k)(1) = param(1) + param(2)*sin(teta);
	    }
	  break;
	}
      case CURVE_ELLIPSE :
	{
	  Real_wp teta1 = this->value_parameter_boundary(i)(0);
	  Real_wp teta2 = this->value_parameter_boundary(i)(1);
	  
	  // length of arc
	  Real_wp L(0);
	  for (int k = 0; k <= this->lob_curve.GetOrder(); k++)
	    {
	      Real_wp teta = (1.0-this->lob_curve.Points(k))*teta1 + this->lob_curve.Points(k)*teta2;
	      Real_wp poids = this->lob_curve.Weights(k);
	      L += poids*sqrt(square(param(2)*sin(teta)) + square(param(3)*cos(teta)));
	    }

	  L *= (teta2-teta1);
	  
	  for (int k = 0; k < lambda.GetM(); k++)
	    {
	      // objective length
	      Real_wp Lobj = L*lambda(k);

	      // Newton algorithm to find the good teta
	      Real_wp err(1e200), err_prec(1e300);
	      Real_wp teta = (1.0-lambda(k))*teta1 + lambda(k)*teta2;
	      while (err < err_prec)
		{
		  Real_wp f = 0, df = 0;
		  for (int q = 0; q <= this->lob_curve.GetOrder(); q++)
		    {
		      Real_wp tetaq = (1.0-this->lob_curve.Points(q))*teta1 + this->lob_curve.Points(q)*teta;
		      Real_wp poids = this->lob_curve.Weights(q);
		      df = sqrt(square(param(2)*sin(tetaq)) + square(param(3)*cos(tetaq)));
		      f += poids*df;
		    }

		  f = (f*(teta-teta1) - Lobj); 
		  err_prec = err; err = abs(f);
		  teta -= f/df;
		}
	      
	      points(k)(0) = param(0) + param(2)*cos(teta);
	      points(k)(1) = param(1) + param(3)*sin(teta);
	    }
	  break;
	}
      case CURVE_PEANUT :
	{
	  Real_wp teta1 = this->value_parameter_boundary(i)(0);
	  Real_wp teta2 = this->value_parameter_boundary(i)(1);

	  // length of arc
	  Real_wp L(0);
	  for (int k = 0; k <= this->lob_curve.GetOrder(); k++)
	    {
	      Real_wp teta = (1.0-this->lob_curve.Points(k))*teta1 + this->lob_curve.Points(k)*teta2;
	      Real_wp poids = this->lob_curve.Weights(k);
	      Real_wp radical = sqrt(1.0+param(4)*square(sin(teta)));
	      Real_wp coef = -radical + param(4)*square(cos(teta))/radical;
	      L += poids*sqrt(square(param(2)*sin(teta)*coef) + square(param(3)*cos(teta)));
	    }
	  
	  L *= (teta2-teta1);
	  
	  for (int k = 0; k < lambda.GetM(); k++)
	    {
	      // objective length
	      Real_wp Lobj = L*lambda(k);

	      // Newton algorithm to find the good teta
	      Real_wp err(1e200), err_prec(1e300);
	      Real_wp teta = (1.0-lambda(k))*teta1 + lambda(k)*teta2;
	      while (err < err_prec)
		{
		  Real_wp f = 0, df = 0;
		  for (int q = 0; q <= this->lob_curve.GetOrder(); q++)
		    {
		      Real_wp tetaq = (1.0-this->lob_curve.Points(q))*teta1 + this->lob_curve.Points(q)*teta;
		      Real_wp poids = this->lob_curve.Weights(q);
		      Real_wp radical = sqrt(1.0+param(4)*square(sin(tetaq)));
		      Real_wp coef = -radical + param(4)*square(cos(tetaq))/radical;
		      df = sqrt(square(param(2)*sin(tetaq)*coef) + square(param(3)*cos(tetaq)));
		      f += poids*df;
		    }
		  
		  f = (f*(teta-teta1) - Lobj); 
		  err_prec = err; err = abs(f);
		  teta -= f/df;
		}
	      
	      Real_wp sin_teta = sin(teta);
	      points(k)(0) = param(0) + param(2)*cos(teta)*sqrt(1.0 + param(4)*sin_teta*sin_teta);
	      points(k)(1) = param(1) + param(3)*sin_teta;
	    }
	  break;
	}
      case CURVE_SPLINE :
	{
	  const VectReal_wp& val = this->Param_curve(ref);
	  int pos1 = this->integer_parameter_boundary(i)(0);
	  int pos2 = this->integer_parameter_boundary(i)(1);
	  for (int k = 0; k < lambda.GetM(); k++)
	    {
	      Real_wp x = lambda(k);
	      Real_wp fct1 = (1.0-x)*((1.0-x)*(1.0-x) - 1.0);
	      Real_wp fct2 = x*(x*x - 1.0);
	      points(k)(0) = val(pos1)*(1.0-x) + val(pos2)*x + val(pos1+2)*fct1 + val(pos2+2)*fct2;
	      points(k)(1)
                = val(pos1+1)*(1.0-x) + val(pos2+1)*x + val(pos1+3)*fct1 + val(pos2+3)*fct2;
	    }
	}
	break;
      case CURVE_LOCAL_SPLINE :
	{
	  int n0 = this->integer_parameter_boundary(i)(0);
	  int n1 = mesh.BoundaryRef(i).numVertex(0);
	  int n2 = mesh.BoundaryRef(i).numVertex(1);
	  int n3 = this->integer_parameter_boundary(i)(1);
	  if (n0 < 0)
	    {
	      for (int k = 0; k < lambda.GetM(); k++)
		{
		  Real_wp x = lambda(k);
		  point.Zero();
		  Add(0.5*(1.0-x)*(2.0-x), mesh.Vertex(n1), point);
		  Add(x*(2.0-x), mesh.Vertex(n2), point);
		  Add(-0.5*x*(1.0-x), mesh.Vertex(n3), point);
		  points(k) = point;
		}
	    }
	  else if (n3 < 0)
	    {
	      for (int k = 0; k < lambda.GetM(); k++)
		{
		  Real_wp x = lambda(k);
		  point.Zero();
		  Add(-0.5*x*(1.0-x), mesh.Vertex(n0), point);
		  Add((1.0+x)*(1.0-x), mesh.Vertex(n1), point);
		  Add(0.5*x*(1.0+x), mesh.Vertex(n2), point);
		  points(k) = point;
		}
	    }
	  else
	    {
	      for (int k = 0; k < lambda.GetM(); k++)
		{
		  Real_wp x = lambda(k);
		  point.Zero();
		  Add(-x*(1.0-x)*(2.0-x)/6, mesh.Vertex(n0), point);
		  Add((1.0+x)*(1.0-x)*(2.0-x)/2, mesh.Vertex(n1), point);
		  Add(x*(1.0+x)*(2.0-x)/2, mesh.Vertex(n2), point);
		  Add(-x*(1.0-x)*(1.0+x)/6, mesh.Vertex(n3), point);
		  points(k) = point;
		}
	    }
	}
        break;
      case CURVE_FILE :
	{
	  for (int k = 0; k < lambda.GetM(); k++)
	    {
	      point.Zero();
	      Add(this->lob_curve.EvaluatePhi(0, lambda(k)),
                  mesh.Vertex(mesh.BoundaryRef(i).numVertex(0)), point);
	      Add(this->lob_curve.EvaluatePhi(this->lob_curve.GetOrder(), lambda(k)),
		  mesh.Vertex(mesh.BoundaryRef(i).numVertex(1)), point);
	      for (int j = 1; j < this->lob_curve.GetOrder(); j++)
		Add(this->lob_curve.EvaluatePhi(j, lambda(k)), PointsEdgeRef(i, j-1), point);
	      
	      points(k) = point;
	    }	  
	}
	break;	
      default :
	abort();
	
      }
  }
  
  
  //! retrieves all the points belonging to a referenced surface, computes normale
  /*!
    \param[in] ref reference of the considered surface
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GetAllPointsOnReference(int ref, VectR_N& Points, VectR_N& Normale) const
  {
    IVect ref_cond(this->GetNbReferences()+1);
    ref_cond.Fill(0);
    ref_cond(ref) = 1;
    
    GetAllPointsOnReference(1, ref_cond, Points, Normale);
  }
  
  
  //! retrieves all the points belonging to a referenced surface, computes normale
  /*!
    \param[in] ref reference of the considered surface
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GetAllPointsOnReference(int ref, const IVect& ref_cond, VectR_N& Points, VectR_N& Normale) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    VectBool VertexOnCircle(mesh_leaf.GetNbVertices());
    VertexOnCircle.Fill(false);
    
    int nb_vertex_circle = 0;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (ref_cond(mesh_leaf.BoundaryRef(i).GetReference()) == ref)
	{
	  if (!VertexOnCircle(mesh_leaf.BoundaryRef(i).numVertex(0)))
	    {
	      VertexOnCircle(mesh_leaf.BoundaryRef(i).numVertex(0)) = true;
	      nb_vertex_circle++;
	    }
	  if (!VertexOnCircle(mesh_leaf.BoundaryRef(i).numVertex(1)))
	    {
	      VertexOnCircle(mesh_leaf.BoundaryRef(i).numVertex(1)) = true;
	      nb_vertex_circle++;
	    }
	}
    
    Points.Reallocate(nb_vertex_circle);
    IVect Index_Vertices(mesh_leaf.GetNbVertices()); Index_Vertices.Fill(-1);
    nb_vertex_circle = 0;
    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
      if (VertexOnCircle(i))
	{
	  Index_Vertices(i) = nb_vertex_circle;
	  Points(nb_vertex_circle++) = mesh_leaf.Vertex(i);
	}
    
    // normales are computed with a mean of normale on each edge
    VectR_N s(2); R_N vec_u, vec_w;
    Normale.Reallocate(nb_vertex_circle); Normale.Fill(vec_u);
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (ref_cond(mesh_leaf.BoundaryRef(i).GetReference()) == ref)
	{
	  for (int j = 0; j < 2; j++)
	    s(j) = mesh_leaf.Vertices(mesh_leaf.BoundaryRef(i).numVertex(j));
	  
	  vec_u = s(1)- s(0);
	  vec_w(0) = vec_u(1); vec_w(1) = -vec_u(0); 
	  
	  for (int j = 0; j < 2; j++)
	    {
	      int nb_vertex_ref = Index_Vertices(mesh_leaf.BoundaryRef(i).numVertex(j));
	      if (DotProd(vec_w, Normale(nb_vertex_ref)) < 0)
		Normale(nb_vertex_ref) -= vec_w;
	      else
		Normale(nb_vertex_ref) += vec_w;
              
	    }
	}
    
    for (int i = 0; i < Normale.GetM(); i++)
      Mlt(T(1.0/Norm2(Normale(i))), Normale(i));
  }
  

  //! selects some points of the referenced surface ref
  /*!
    \param[in] ref reference of the considered surface
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
    \param[in] nb_points number of points to find
    This methods selects points so that points are enough distant
  */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GetDistantPointsOnReference(int ref, VectR_N& Points, VectR_N& Normale, int nb_points) const
  {
    IVect ref_cond(this->GetNbReferences()+1);
    ref_cond.Fill(0);
    ref_cond(ref) = 1;
    
    GetDistantPointsOnReference(1, ref_cond, Points, Normale, nb_points);
  }
  
  
  //! selects some points of the referenced surface ref
  /*!
    \param[in] ref reference of the considered surface
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
    \param[in] nb_points number of points to find
    This methods selects points so that points are enough distant
  */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GetDistantPointsOnReference(int ref, const IVect& ref_cond, VectR_N& Points,
			      VectR_N& Normale, int nb_points) const
  {
    GetAllPointsOnReference(ref, ref_cond, Points, Normale);
    if (Points.GetM() < nb_points)
      return;
    
    IVect num(nb_points);
    // first point
    num(0) = 0;
    VectR_N NewPoints(nb_points);
    VectR_N NewNormale(nb_points);
    NewPoints(0) = Points(num(0));
    NewNormale(0) = Normale(num(0));
    // choice of points
    for (int current = 1; current < nb_points; current++)
      {
	T dist_max = 0.0; T dist_vert = 0.0;
	int imax = -1;
	for (int i = 0; i < Points.GetM(); i++)
	  {
	    dist_vert = Points(i).Distance(Points(num(0)));
	    for (int j = 0; j < current; j++)
	      dist_vert = min(dist_vert, Points(i).Distance(Points(num(j))));
	    
	    if (dist_vert >= dist_max)
	      {
		imax = i;
		dist_max = dist_vert;
	      }	
	  }
	
	num(current) = imax;
	NewPoints(current) = Points(imax);
	NewNormale(current) = Normale(imax);
      }
    Points = NewPoints; Normale = NewNormale;
  
  }


  //! finds parameters of curves
  /*!
    For example, if a curve is a circle arc, this method will find the center and
    radius of the circle.
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::FindParametersCurve()
  {
    Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    // estimation of parameters if not given by the user	
    for (int i = 1; i <= this->GetNbReferences(); i++)
      {
	if ((this->Type_curve(i) > 0)&&(!this->Param_known(i)))
	  {
	    if (this->Type_curve(i)== this->CURVE_CIRCLE)
	      {
		T xc, yc, radius;
		FindParametersCircle(i, xc, yc, radius);
                if (mesh_leaf.print_level >= 4)
                  {
                    cout << " Curve "<<i<<" is a circle of center "<<
                      xc<<","<<yc<<"  and radius  "<<radius<<endl;
                  }
		
		this->Param_curve(i).Reallocate(3);
		this->Param_curve(i)(0) = xc;
		this->Param_curve(i)(1) = yc;
		this->Param_curve(i)(2) = radius;
		this->Param_known(i) = true;
	      }
	    else if (this->Type_curve(i)== this->CURVE_ELLIPSE)
	      {
		T xc, yc, a, b;
		FindParametersEllipse(i, xc, yc, a, b);
		
		this->Param_curve(i).Reallocate(4);
		this->Param_curve(i)(0) = xc;
		this->Param_curve(i)(1) = yc;
		this->Param_curve(i)(2) = a;
		this->Param_curve(i)(3) = b;
		this->Param_known(i) = true;
	      }
	    else if (this->Type_curve(i)== this->CURVE_PEANUT)
	      {
		T xc, yc, a, b, c;
		FindParametersPeanut(i, xc, yc, a, b, c);
		
		this->Param_curve(i).Reallocate(5);
		this->Param_curve(i)(0) = xc;
		this->Param_curve(i)(1) = yc;
		this->Param_curve(i)(2) = a;
		this->Param_curve(i)(3) = b;
		this->Param_curve(i)(4) = c;
		this->Param_known(i) = true;
	      }
	    else if (this->Type_curve(i) == this->CURVE_SPLINE)
	      {
		// nothing to do
	      }
	    else if (this->Type_curve(i) == this->CURVE_LOCAL_SPLINE)
	      {
		// nothing to do
	      }
	    else if (this->Type_curve(i) == this->CURVE_FILE)
	      {
		// nothing to find
	      }
	    else
	      {
		cout<<"Case not treated "<<endl;
		abort();
	      }
	  }
      }
  }
  
  
  //! the curve ref is assumed to be a line ax + by +c = 0, parameters a,b and c are estimated
  /*!
    \param[in] ref reference
    \param[out] coef_plane coefficients a, b and c
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  FindParametersPlane(int ref, Vector<T>& coef_plane) const
  {
    VectR_N Points, Normale;
    int nb_points = 2;
    GetDistantPointsOnReference(ref, Points, Normale, nb_points);
    if (Points.GetM() < nb_points)
      return;
    
    coef_plane.Reallocate(3);
    coef_plane(0) = Normale(0)(0); coef_plane(1) = Normale(0)(1);
    coef_plane(2) = -DotProd(Normale(0), Points(0));
  }

  
  //! the curve ref is supposed to be a part of circle, center and radius are estimated
  /*!
    \param[in] ref reference
    \param[out] xc abscisse of the center
    \param[out] yc ordonnee of the center
    \param[out] radius radius of the circle
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  FindParametersCircle(int ref, T& xc, T& yc, T& radius) const
  {
    xc = 0.0; yc = 0.0; radius = 0.0;
    
    VectR_N Points, Normale;
    int nb_points = 3;
    GetDistantPointsOnReference(ref, Points, Normale, nb_points);
    if (Points.GetM() < nb_points)
      return;
    
    // estimation of radius and center
    R_N middle1, middle2, pt1, pt2, center;
    Seldon::Add(Points(0), Points(1), middle1); middle1 *= 0.5;
    Seldon::Add(Points(1), Points(2), middle2); middle2 *= 0.5;
    pt1(0) = middle1(0) + Points(1)(1) - Points(0)(1);
    pt1(1) = middle1(1) - Points(1)(0) + Points(0)(0);
    pt2(0) = middle2(0) + Points(1)(1) - Points(2)(1);
    pt2(1) = middle2(1) - Points(1)(0) + Points(2)(0);
    
    int nb_intersec = IntersectionDroites(middle1, pt1, middle2, pt2, center, epsilon_machine);
    if (nb_intersec != 1)
      {
        cout << "Failed to find an intersection in FindParametersCircle" << endl;
        abort();
      }
    
    xc = center(0); yc = center(1);
    radius = center.Distance(Points(0));
  }
  
  
  //! the curve ref is supposed to be a part of ellipse, center and radii are estimated
  /*!
    \param[in] ref reference
    \param[out] xc x-coordinate of the center
    \param[out] yc y-coordinate of the center
    \param[out] a x-radius of the ellipse
    \param[out] b y-radius of the ellipse
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  FindParametersEllipse(int ref, T& xc, T& yc, T& a, T& b) const
  {
    xc = 0.0; yc = 0.0; a = 0.0; b = 0;
    abort();
  }
  
  
  //! the curve ref is supposed to be a part of peanut, center and other parameters are estimated
  /*!
    \param[in] ref reference
    \param[out] xc x-coordinate of the center
    \param[out] yc y-coordinate of the center
    \param[out] a parameter a
    \param[out] b parameter b
    \param[out] c parameter c
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  FindParametersPeanut(int ref, T& xc, T& yc, T& a, T& b, T& c) const
  {
    xc = 0.0; yc = 0.0; a = 0.0; b = 0; c = 0;
    abort();
  }
  
  
  //! finds parameters for spline
  /*!
    \param[in] ref reference
    \param[out] param parameters defining the spline
   */
  template<class T>
  void MeshBoundaries<Dimension2, T>
  ::FindParametersSpline(int ref, Vector<T>& param, IVect& Point_curve)
  {
    Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    
    // first step : we find an initial point (ie an extremity if open spline)
    IVect NbEdgesVertices(mesh_leaf.GetNbVertices());
    NbEdgesVertices.Fill(0);
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (mesh_leaf.BoundaryRef(i).GetReference() == ref)
	{
	  NbEdgesVertices(mesh_leaf.BoundaryRef(i).numVertex(0))++;
	  NbEdgesVertices(mesh_leaf.BoundaryRef(i).numVertex(1))++;
	}
    
    int num_init = -1;
    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
      if (NbEdgesVertices(i) > 0)
	{
	  if (num_init == -1)
	    num_init = i;
	  else if (NbEdgesVertices(i) == 1)
	    num_init = i;
	  else if (NbEdgesVertices(i) >= 3)
	    {
	      cout << "Spline with a branch"<<endl;
	      abort();
	    }
	}
    
    // if no initial point, the spline is empty
    if (num_init == -1)
      return;
    
    // second step : loop to get all the points of the spline
    Point_curve.Reallocate(mesh_leaf.GetNbBoundaryRef()+1);
    IVect ref_cond(this->GetNbReferences()); ref_cond.Fill();
    int ipoint = num_init, iback = -1, inext = -1, iedge = -1, nb_points = 0;
    while ((nb_points == 0)||((ipoint != num_init)&&(ipoint != -1)))
      {
	Point_curve(nb_points) = ipoint;
	inext = FindFollowingVertex(iback, ipoint, ref, ref_cond, iedge);
	iback = ipoint; ipoint = inext;
	nb_points++;
      }
    
    // storing coordinates of all the points of the spline in Pts
    VectR_N Pts(nb_points+2);
    Pts(0) = mesh_leaf.Vertex(Point_curve(nb_points-1));
    Pts(nb_points+1) = mesh_leaf.Vertex(Point_curve(0));
    for (int i = 1; i <= nb_points; i++)
      Pts(i) = mesh_leaf.Vertex(Point_curve(i-1));
    
    if (nb_points == 2)
      {
	param.Reallocate(8);
	param.Fill(0);
	param(0) = Pts(1)(0); param(1) = Pts(1)(1);
	param(4) = Pts(2)(0); param(5) = Pts(2)(1);
	return;
      }

    param.Reallocate(4*nb_points); param.Fill(0);
    
    if (NbEdgesVertices(num_init) == 2)
      {
	// closed spline
	Vector<T> beta(nb_points), delta(nb_points), alpha(nb_points), coef(nb_points);
	alpha.Fill(0); beta.Fill(0); delta.Fill(4); coef.Fill(0);
	
	// we are solving the following linear system
	// |  4   1     0   ...   1  | param =     | x_0 + 2 x_1 - x_2
	// |  1   4     1   ...   0  |             | x_1 + 2 x_2 - x_3
	// |  0   1  ...     1    0  |             |     ...	    
	// |  0   ...   1    4    1  |             | x_(i-1) + 2 x_i - x_(i+1)
	// |  1    0   ...   1    4  |             |     ...
	    
	// we perform an LU factorization of this matrix
	// we store L in the form 
	
	// | delta_0  
	// | alpha_0  delta_1  
	// |                   ...              
	// |                         delta_n-2  
	// |  coef_0   coef_1  ...    alpha_n-2  delta_n-1
	
	// factorisation LU de la matrice
	coef(0) = 1.0; Real_wp pivot;
	for (int i = 0; i < nb_points-2; i++)
	  {
	    alpha(i) = 1.0/delta(i);
	    delta(i+1) -= alpha(i);
	    pivot = coef(i)*alpha(i);
	    delta(nb_points-1) -= pivot*coef(i);
	    coef(i) = pivot;
	    coef(i+1) = -pivot;
	  }
	
	coef(nb_points-2) += 1.0;
	alpha(nb_points-2) = coef(nb_points-2)/delta(nb_points-2);
	delta(nb_points-1) -= alpha(nb_points-2)*coef(nb_points-2);
	coef(nb_points-2) = alpha(nb_points-2);
	
	// we replace diagonal by its inverse
	for (int i = 0; i < nb_points; i++)
	  delta(i) = 1.0/delta(i);
	
	// loop on the components (x and y)
	for (int n = 0; n < 2; n++)
	  {	    
	    // right hand side is stored in beta
	    for (int i = 0; i < nb_points; i++)
	      beta(i) = Pts(i)(n) - 2.0*Pts(i+1)(n) + Pts(i+2)(n);
	    
	    // on resout L x = beta
	    for (int i = 1; i < nb_points-1; i++)
	      {
		beta(i) -= alpha(i-1)*beta(i-1);
		beta(nb_points-1) -= coef(i-1)*beta(i-1);
	      }
	    
	    beta(nb_points-1) -= alpha(nb_points-2)*beta(nb_points-2);
	    // on resout U x = beta
	    beta(nb_points-1) *= delta(nb_points-1);
	    beta(nb_points-2) = beta(nb_points-2)*delta(nb_points-2) 
              - coef(nb_points-2)*beta(nb_points-1);
	    
            for (int i = nb_points-3; i >= 0; i--)
	      beta(i) = beta(i)*delta(i) - alpha(i)*beta(i+1) - coef(i)*beta(nb_points-1);
	    
	    for (int i = 0; i < nb_points; i++)
	      {
		param(4*i+n) = Pts(i+1)(n);
		param(4*i+n+2) = beta(i);
	      }
	    
	  }
      }
    else
      {
	// open spline
	Vector<T> beta(nb_points-2), delta(nb_points-2), alpha(nb_points-2);
	alpha.Fill(0); beta.Fill(0); delta.Fill(4);
	
	// we are solving the following linear system
	// |  4   1     0   ...   0  | param =     | x_0 + 2 x_1 - x_2
	// |  1   4     1   ...   0  |             | x_1 + 2 x_2 - x_3
	// |  0   1  ...     1    0  |             |     ...	    
	// |  0   ...   1    4    1  |             | x_(i-1) + 2 x_i - x_(i+1)
	// |  0    0   ...   1    4  |             |     ...
	    
	// we perform an LU factorization of this matrix
	// we store L in the form 
	
	// | delta_0  
	// | alpha_0  delta_1  
	// |                   ...              
	// |                         delta_n-2  
	// |  0                      alpha_n-2  delta_n-1
	
	// factorisation LU de la matrice
	for (int i = 0; i < nb_points-3; i++)
	  {
	    alpha(i) = 1.0/delta(i);
	    delta(i+1) -= alpha(i);
	  }
	
	// we replace diagonal by its inverse
	for (int i = 0; i < nb_points-2; i++)
	  delta(i) = 1.0/delta(i);
	
	// loop on the components (x and y)
	for (int n = 0; n < 2; n++)
	  {	    
	    // right hand side is stored in beta
	    for (int i = 1; i < nb_points-1; i++)
	      beta(i-1) = Pts(i)(n) - 2.0*Pts(i+1)(n) + Pts(i+2)(n);
	    	    
	    // on resout L x = beta
	    for (int i = 1; i < nb_points-2; i++)
	      beta(i) -= alpha(i-1)*beta(i-1);
	    
	    // on resout U x = beta
	    beta(nb_points-3) *= delta(nb_points-3);
	    for (int i = nb_points-4; i >= 0; i--)
	      beta(i) = beta(i)*delta(i) - beta(i+1)*alpha(i);
	    
	    param(n) = Pts(1)(n);
	    for (int i = 1; i < nb_points-1; i++)
	      {
		param(4*i+n) = Pts(i+1)(n);
		param(4*i+n+2) = beta(i-1);
	      }
	    param(4*(nb_points-1)+n) = Pts(nb_points)(n);
	  }
      }
  }

  
  /*****************
   * Other methods *
   *****************/
  
  
  //! returns the center of radial PML
  template<class T> TinyVector<T, 2>
  MeshBoundaries<Dimension2, T>::GetCenterRadialPML(T& Rmax) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();

    TinyVector<T, 2> origin;
    Rmax = 0;
    for (int ref = 0; ref <= mesh_leaf.GetNbReferences(); ref++)
      if (mesh_leaf.GetCurveType(ref) == mesh_leaf.CURVE_CIRCLE)
        {
          const Vector<T>& param = mesh_leaf.GetCurveParameter(ref);
          if (param(2) > Rmax)
            {
              origin.Init(param(0), param(1));
              Rmax = param(2);
            }
        }
    
    if (Rmax == 0)
      {
        for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
          Rmax = max(Rmax, mesh_leaf.Vertex(i).Distance(origin));
      }
    
    return origin;
  }


  //! extraction of surfacic mesh from a volumic mesh
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GetBoundaryMesh(int ref, SurfacicMesh<Dimension2>& mesh,
                  const IVect& ref_cond, int ref_domain) const
  {
    GetBoundaryMesh(ref, mesh, mesh.ListeBoundaries, mesh.ListeVertices,
                    mesh.NumElement, mesh.NumLocalBoundary, ref_cond, ref_domain);
  }
  
  
  //! extraction of surfacic mesh from a volumic mesh
  /*!
    \param[in] ref reference of the faces to be extracted
    \param[out] mesh_surf surfacic mesh extracted
    \param[out] Index_EdgeSurf_to_EdgeRef numbers of referenced edges extracted
    \param[out] Liste_Vertices numbers of vertices extracted
    \param[in] ref_cond we are extracting references ref_cond(i) == ref
  */
  template<class T> template<class Mesh2D_>
  void MeshBoundaries<Dimension2, T>::
  GetBoundaryMesh(int ref, Mesh2D_& mesh_surf, IVect& Index_EdgeSurf_to_EdgeRef,
		  IVect& Liste_Vertices, IVect& NumElement, IVect& NumLocalBoundary,
                  const IVect& ref_cond, int ref_domain) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    int nf = 0;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (ref_cond(mesh_leaf.BoundaryRef(i).GetReference()) == ref)
        {
          bool edge_taken = true;
          if (ref_domain > 0)
            {
              bool elt_ref = false;
              int num_elem = mesh_leaf.BoundaryRef(i).numElement(0);
              if (num_elem >= 0)                  
                if (mesh_leaf.Element(num_elem).GetReference() == ref_domain)
                  elt_ref = true;
              
              num_elem = mesh_leaf.BoundaryRef(i).numElement(1);
              if (num_elem >= 0)                  
                if (mesh_leaf.Element(num_elem).GetReference() == ref_domain)
                  elt_ref = true;
              
              if (!elt_ref)
                edge_taken = false;
            }
          
          if (edge_taken)
            nf++;
        }
    
    if (nf <= 0)
      {
	mesh_surf.Clear(); Index_EdgeSurf_to_EdgeRef.Clear();
	Liste_Vertices.Clear();
	return;
      }
    
    mesh_surf.ClearCurves(); 
    mesh_surf.SetGeometryOrder(this->GetGeometryOrder());
    mesh_surf.Type_curve = mesh_leaf.Type_curve;
    mesh_surf.Param_curve = mesh_leaf.Param_curve;
    mesh_surf.Param_known = mesh_leaf.Param_known;
    mesh_surf.ReallocateBoundariesRef(nf);
    Index_EdgeSurf_to_EdgeRef.Reallocate(nf);
    Liste_Vertices.Reallocate(2*nf);
    NumElement.Reallocate(nf);
    NumLocalBoundary.Reallocate(nf);
    int ind = 0; nf = 0;
    // we retrieve all the vertex numbers and ege numbers
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	int ref_edge = mesh_leaf.BoundaryRef(i).GetReference();
	ref_edge = ref_cond(ref_edge);
	if (ref_edge == ref)
	  {
            bool edge_taken = true;
            int num_elem = -1;
            if (ref_domain > 0)
              {
                bool elt_ref = false;
                int num_elem2 = mesh_leaf.BoundaryRef(i).numElement(0);
                if (num_elem2 >= 0)
                  if (mesh_leaf.Element(num_elem2).GetReference() == ref_domain)
                    {
                      elt_ref = true;
                      num_elem = num_elem2;
                    }
                
                num_elem2 = mesh_leaf.BoundaryRef(i).numElement(1);
                if (num_elem2 >= 0)
                  if (mesh_leaf.Element(num_elem2).GetReference() == ref_domain)
                    {
                      elt_ref = true;
                      num_elem = num_elem2;
                    }
                
                if (!elt_ref)
                  edge_taken = false;
              }
            else
              num_elem = mesh_leaf.BoundaryRef(i).numElement(0);
            
            if (edge_taken)
              {
                mesh_surf.edges_ref(nf) = mesh_leaf.edges_ref(i);
                Index_EdgeSurf_to_EdgeRef(nf) = i;
                for (int j = 0; j < 2; j++)
                  {
                    Liste_Vertices(ind) = mesh_leaf.edges_ref(i).numVertex(j);
                    ind ++;
                  }
                
                if (mesh_leaf.GetNbElt() > 0)
                  {
                    NumElement(nf) = num_elem;
                    NumLocalBoundary(nf) = mesh_leaf.elements(num_elem)
                      .GetPositionBoundary(i);
                  }
                
                nf++;
              }
          }
      }
    
    // we sort list of vertices
    Assemble(ind, Liste_Vertices);
    Liste_Vertices.Resize(ind);
    
    // affecting new vertex numbers to the list of edges
    mesh_surf.ReallocateVertices(ind);
    IVect NewNumber_Vertices(mesh_leaf.GetNbVertices());
    for (int i = 0; i < ind; i++)
      {
	mesh_surf.Vertices(i) = mesh_leaf.Vertices(Liste_Vertices(i));
	NewNumber_Vertices(Liste_Vertices(i)) = i;
      }
    
    for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
      {
	int ref_edge = mesh_surf.BoundaryRef(i).GetReference();
	
	int n0 = NewNumber_Vertices(mesh_surf.BoundaryRef(i).numVertex(0));
	int n1 = NewNumber_Vertices(mesh_surf.BoundaryRef(i).numVertex(1));
	
	mesh_surf.BoundaryRef(i).Init(n0, n1, ref_edge);
      }
    
    // edges
    mesh_surf.edges = mesh_surf.edges_ref;
    
    // internal points
    for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
      {
	int ne = Index_EdgeSurf_to_EdgeRef(i);
	for (int k = 0; k < mesh_surf.PointsEdgeRef.GetN(); k++)
	  {	    
	    mesh_surf.PointsEdgeRef(i, k) = PointsEdgeRef(ne, k);
	  }
	mesh_surf.value_parameter_boundary(i) = value_parameter_boundary(ne);
	mesh_surf.integer_parameter_boundary(i) = integer_parameter_boundary(ne);
      }
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! extrudes boundary on a 2-D mesh in order to generate a surfacic mesh of triangles
  /*! 
    \param[in] ref reference target
    \param[in] mesh_3d 3-D triangular mesh
    \param[out] Index_FaceSurf_to_EdgeRef edge number from which
    3D triangle i has been extruded
    \param[out] Index_VertexSurf_to_Vertex 2-D vertex number from which
    3-D vertex i has been extruded
    \param[out] AngleVertex angle \f$ \theta \f$ for the vertex i
    \param[in] ref_cond edges of reference i such that ref_cond(i) = ref are extruded
  */
  template<class T>
  void MeshBoundaries<Dimension2, T>::
  GenerateSurfaceOfRevolution(int ref, Mesh<Dimension3>& mesh_3d,
			      IVect& Index_FaceSurf_to_EdgeRef,
			      IVect& Index_VertexSurf_to_Vertex,
			      Vector<T>& AngleVertex, const IVect& ref_cond) const
  {
    const Mesh<Dimension2, T>& mesh_leaf = this->GetLeafClass();
    // first we get list of vertices and edges of the 2-D, which will be extruded
    // we compute too the space step of the boundary mesh
    T space_step = 0;
    IVect ListVertices(2*mesh_leaf.GetNbBoundaryRef()), ListEdges(mesh_leaf.GetNbBoundaryRef());
    IVect Index_Vertex(mesh_leaf.GetNbVertices()); Index_Vertex.Fill(-1);
    int nb_edges_boundary = 0; int nb_vertices_boundary = 0;
    int n1, n2, na, nb, num_edge;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	int ref_edge = mesh_leaf.edges_ref(i).GetReference();
	
	if (ref_cond(ref_edge) == ref)
	  {
	    ListEdges(nb_edges_boundary++) = i;
	    n1 = mesh_leaf.edges_ref(i).numVertex(0);
	    n2 = mesh_leaf.edges_ref(i).numVertex(1);
	    if (Index_Vertex(n1) == -1)
	      {
		Index_Vertex(n1) = nb_vertices_boundary;
		ListVertices(nb_vertices_boundary++) = n1;
	      }
	    if (Index_Vertex(n2) == -1)
	      {
		Index_Vertex(n2) = nb_vertices_boundary;
		ListVertices(nb_vertices_boundary++) = n2;
	      }
	    
	    space_step += mesh_leaf.Vertices(n1).Distance(mesh_leaf.Vertices(n2));
	  }
      }
    
    space_step /= nb_edges_boundary;
    ListVertices.Resize(nb_vertices_boundary);
    ListEdges.Resize(nb_edges_boundary);
    
    // We get for each 2-D vertex of the boundary the two edges ,
    // which are before and after this point
    Matrix<int> NumEdge(nb_vertices_boundary,2); NumEdge.Fill(-1);
    for (int i = 0; i < nb_edges_boundary; i++)
      {
	num_edge = ListEdges(i);
	n1 = mesh_leaf.edges_ref(num_edge).numVertex(0);
	n2 = mesh_leaf.edges_ref(num_edge).numVertex(1);
	na = Index_Vertex(n1); nb = Index_Vertex(n2);
	if (NumEdge(na,0) == -1)
	  NumEdge(na,0) = i;
	else
	  NumEdge(na,1) = i;
	
	if (NumEdge(nb,0) == -1)
	  NumEdge(nb,0) = i;
	else
	  NumEdge(nb,1) = i;
      }
    
    // we suppose that we have only one curve
    // we search initial point (point, which owns to only one referenced edge)
    int initial_point = -1;
    for (int i = 0; i < nb_vertices_boundary; i++)
      if (NumEdge(i,1) == -1)
	initial_point = i;
    
    // DamierVertex is used to alternate the theta angle of the start point
    IVect DamierVertex(nb_vertices_boundary); DamierVertex.Zero();
    int ipoint = initial_point, iprec = ipoint, inew, bascule = 0;
    while (ipoint != -1)
      {
	DamierVertex(ipoint) = bascule;
	inew = -1;
	for (int k = 0; k < 2; k++)
	  {
	    if (NumEdge(ipoint,k) != -1)
	      {
		num_edge = ListEdges(NumEdge(ipoint,k));
		for (int l = 0; l < 2; l++)
		  {
		    n1 = mesh_leaf.edges_ref(num_edge).numVertex(l);
		    na = Index_Vertex(n1);
		    if ((na != iprec)&&(na != ipoint))
		      inew = na;
		  }
	      }
	  }
	iprec = ipoint;
	ipoint = inew;
	bascule = 1-bascule;
      }
    
    T xmax_ = mesh_leaf.Vertices(0)(0);
    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
      xmax_ = std::max(mesh_leaf.Vertices(i)(0), xmax_);
    
    // we extrude the mesh by rotation, axis of rotation is Oz
    // first we create vertices
    R_N point;  T theta;
    // un majorant du nombre de sommets
    int nb_vert3d = max(toInteger(ceil(xmax_*2*pi_wp/space_step)),6) + 2;
    nb_vert3d *= nb_vertices_boundary; 
    
    mesh_3d.SetGeometryOrder(this->GetGeometryOrder());
    mesh_3d.ReallocateVertices(nb_vert3d);
    AngleVertex.Reallocate(nb_vert3d);
    Index_VertexSurf_to_Vertex.Reallocate(nb_vert3d);
    nb_vert3d = 0;
    
    IVect CumulNbPoints(nb_vertices_boundary+1); IVect NbPoints(nb_vertices_boundary);
    for (int i = 0; i < nb_vertices_boundary; i++)
      {
	CumulNbPoints(i) = nb_vert3d;
	n1 = ListVertices(i);
	point = mesh_leaf.Vertices(n1);
	// number of points on the circle depends on the radius r = point(0)
	int nb_points_circle = max(toInteger(ceil(2*pi_wp*point(0)/space_step)),5) + 1;
	// we test if the point is on the axis
	// if not          -> we put at least 3 points, to avoid degenerated triangles
	if (abs(point(0)) <= R_N::threshold)
	  {
	    // if it is a case -> we put only one point (sommet d'un cone)
	    nb_points_circle = 1;
	    AngleVertex(nb_vert3d) = T(0);
	    Index_VertexSurf_to_Vertex(nb_vert3d) = n1;
	    mesh_3d.Vertex(nb_vert3d++).Init(0,0,point(1));
	  }
	else
	  {
	    nb_points_circle = std::max(nb_points_circle,3);
	    if (DamierVertex(i) == 0)
	      {
		// first vertex is so that y = 0
		for (int j=0 ; j<nb_points_circle ; j++)
		  {
		    theta = 2*pi_wp*j/nb_points_circle;
		    AngleVertex(nb_vert3d) = theta;
		    Index_VertexSurf_to_Vertex(nb_vert3d) = n1;
		    mesh_3d.Vertex(nb_vert3d++).
		      Init(point(0)*cos(theta),point(0)*sin(theta),point(1));
		  }
	      }
	    else
	      {
		// the middle of the first vertex and the last vertex is so that y = 0
		for (int j=0 ; j<nb_points_circle ; j++)
		  {
		    theta = 2*pi_wp*(0.5+j)/nb_points_circle;
		    AngleVertex(nb_vert3d) = theta;
		    Index_VertexSurf_to_Vertex(nb_vert3d) = n1;
		    mesh_3d.Vertex(nb_vert3d++).
		      Init(point(0)*cos(theta),point(0)*sin(theta),point(1));
		  }
	      }
	  }
	NbPoints(i) = nb_points_circle;
      }
    CumulNbPoints(nb_vertices_boundary) = nb_vert3d;
    mesh_3d.ResizeVertices(nb_vert3d);
    AngleVertex.Resize(nb_vert3d);
    Index_VertexSurf_to_Vertex.Resize(nb_vert3d);
    
    // now we create list of triangles
    mesh_3d.ReallocateBoundariesRef(2*nb_vert3d);
    Index_FaceSurf_to_EdgeRef.Reallocate(2*nb_vert3d);
    int nb_faces3d = 0;
    for (int i = 0; i < nb_edges_boundary; i++)
      {
	num_edge = ListEdges(i);
	int ref_edge = mesh_leaf.edges_ref(num_edge).GetReference();
	n1 = mesh_leaf.edges_ref(num_edge).numVertex(0);
	n2 = mesh_leaf.edges_ref(num_edge).numVertex(1);
	na = Index_Vertex(n1); nb = Index_Vertex(n2);
	// DISP(n1); DISP(n2); DISP(na); DISP(nb);
	
	// int nb_ptsA = NbPoints(na); int nb_ptsB = NbPoints(nb);
	int posA = CumulNbPoints(na); int posB = CumulNbPoints(nb);
	int endA = CumulNbPoints(na+1)-1; int endB = CumulNbPoints(nb+1)-1;
	int ia = posA; int ib = posB;
	while ((ia!=endA)||(ib!=endB))
	  {
	    // we search following point
	    inew = -1; bool new_point_A = false;
	    if (ia == endA)
	      {
		inew = ib+1;
		new_point_A = false;
	      }
	    else if (ib == endB)
	      {
		inew = ia+1;
		new_point_A = true; 
	      }
	    else
	      {
		if (AngleVertex(ia+1) < AngleVertex(ib+1))
		  {
		    inew = ia+1;
		    new_point_A = true; 
		  }
		else
		  {
		    inew = ib+1;
		    new_point_A = false; 
		  }
	      }
	    // we create triangle ia, ib  and inew
	    // triangle is oriented, so that normale is inward to the boundary mesh
	    Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
	    mesh_3d.BoundaryRef(nb_faces3d++).InitTriangular(ia,ib,inew,ref_edge);
	    if (new_point_A)
	      ia = inew;
	    else
	      ib = inew;
	    
	  } // end loop on circles
	
	// we don't forget the two last triangles
	if (AngleVertex(posA) > AngleVertex(posB))
	  {
	    if (ia == posA)
	      {
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(ia,ib,posB,ref_edge);
	      }
	    else if (ib==posB)
	      {
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(ia,posB,posA,ref_edge);
	      }
	    else
	      {
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(ia,ib,posB,ref_edge);
		
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(ia,posB,posA,ref_edge);
	      }
	  }
	else
	  {
	    if (ia == posA)
	      {
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(posA,ib,posB,ref_edge);
	      }
	    else if (ib==posB)
	      {
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(ia,ib,posA,ref_edge);
	      }
	    else
	      {
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(ia,ib,posA,ref_edge);
		Index_FaceSurf_to_EdgeRef(nb_faces3d) = num_edge;
		mesh_3d.BoundaryRef(nb_faces3d++).
		  InitTriangular(posA,ib,posB,ref_edge);
	      }
	  }
	
      } // end loop on boundary edges
    
    mesh_3d.ResizeBoundariesRef(nb_faces3d);
    Index_FaceSurf_to_EdgeRef.Resize(nb_faces3d);
    // DISP(Index_FaceSurf_to_EdgeRef);
    
    // now, we reorient triangles, so that normale (computed with dF/dx \times dF/dy)
    // is outward to the computational domain
    if (mesh_leaf.GetNbElt() > 0)
      {
	TinyVector<T, 3> dF_dx,dF_dy; TinyVector<T, 3> normale,normale_rotation; R_N normale_2D;
	for (int i = 0; i < nb_faces3d; i++)
	  {
	    na = mesh_3d.BoundaryRef(i).numVertex(0);
	    nb = mesh_3d.BoundaryRef(i).numVertex(1);
	    int nc = mesh_3d.BoundaryRef(i).numVertex(2);
	    dF_dx = mesh_3d.Vertex(nb) - mesh_3d.Vertex(na);
	    dF_dy = mesh_3d.Vertex(nc) - mesh_3d.Vertex(na);
	    TimesProd(dF_dx,dF_dy,normale);
	    
	    // now we compute normale on the 2-D mesh by rotation
	    n1 = Index_VertexSurf_to_Vertex(na);
	    n2 = Index_Vertex(n1); num_edge = ListEdges(NumEdge(n2,0));
	    int num_elem = mesh_leaf.edges_ref(num_edge).numElement(0);
	    int num_loc = mesh_leaf.elements(num_elem).GetPositionBoundary(num_edge);
	    int type_elt = mesh_leaf.elements(num_elem).GetHybridType();
	    
            int num_point1 = mesh_leaf.elements(num_elem).
	      numVertex(MeshNumbering<Dimension2>::FirstExtremityEdge(type_elt, num_loc));
            
	    int num_point2 = mesh_leaf.elements(num_elem).
	      numVertex(MeshNumbering<Dimension2>::SecondExtremityEdge(type_elt,num_loc));
	    
            normale_2D(0) = mesh_leaf.Vertices(num_point2)(1) - mesh_leaf.Vertices(num_point1)(1);
	    normale_2D(1) = mesh_leaf.Vertices(num_point1)(0) - mesh_leaf.Vertices(num_point2)(0);
	    theta = AngleVertex(na);
	    normale_rotation(0) = normale_2D(0)*cos(theta);
	    normale_rotation(1) = normale_2D(0)*sin(theta);
	    normale_rotation(2) = normale_2D(1);
	    
	    T prod_scal = DotProd(normale_rotation,normale);
	    if (prod_scal < 0)
	      {
		// normale are near to be opposite, we change orientation of triangle
		mesh_3d.BoundaryRef(i).
		  InitTriangular(na,nc,nb,mesh_3d.BoundaryRef(i).GetReference());
	      }
	  }
      }
    
    // last step : for curved 2-D mesh, we generate a curved 3-D mesh
    mesh_3d.FindConnectivity();
    Vector<bool> RefCurved(mesh_3d.GetNbReferences()); RefCurved.Fill(false);
    int nb_edges_3d = mesh_3d.GetNbEdges();
    Vector<bool> EdgeUsed(nb_edges_3d); EdgeUsed.Fill(false);
    mesh_3d.ReallocateEdgesRef(nb_edges_3d);
    nb_edges_3d = 0;
    for (int i = 0; i < mesh_3d.GetNbBoundaryRef(); i++)
      {
	int ns = Index_FaceSurf_to_EdgeRef(i);
	int ref = mesh_leaf.BoundaryRef(ns).GetReference();
	if (this->Type_curve(ref) > 0)
	  {
	    RefCurved(ref) = true;
	    int n1 = mesh_3d.BoundaryRef(i).numVertex(0);
	    int n2 = mesh_3d.BoundaryRef(i).numVertex(1);
	    int n3 = mesh_3d.BoundaryRef(i).numVertex(2);
	    int n1s = Index_VertexSurf_to_Vertex(n1);
	    int n2s = Index_VertexSurf_to_Vertex(n2);
	    int n3s = Index_VertexSurf_to_Vertex(n3);
	    
	    Real_wp teta1 = AngleVertex(n1);
	    Real_wp teta2 = AngleVertex(n2);
	    Real_wp teta3 = AngleVertex(n3);
	    
	    Real_wp s1(0), s2(0), s3(0);
	    if (n1s == mesh_leaf.BoundaryRef(ns).numVertex(1))
	      s1 = 1.0;
	    
	    if (n2s == mesh_leaf.BoundaryRef(ns).numVertex(1))
	      s2 = 1.0;
	    
	    if (n3s == mesh_leaf.BoundaryRef(ns).numVertex(1))
	      s3 = 1.0;
	    
	    bool pt1_on_axis = abs(mesh_leaf.Vertex(n1s)(0) <= 1e-10);
	    bool pt2_on_axis = abs(mesh_leaf.Vertex(n2s)(0) <= 1e-10);
	    bool pt3_on_axis = abs(mesh_leaf.Vertex(n3s)(0) <= 1e-10);
	    if (pt1_on_axis)
	      {
		if (teta3 > teta2 + pi_wp)
		  teta3 -= 2.0*pi_wp;
		if (teta3 < teta2 - pi_wp)
		  teta3 += 2.0*pi_wp;
	      }
	    else if (pt2_on_axis)
	      {
		if (teta3 > teta1 + pi_wp)
		  teta3 -= 2.0*pi_wp;
		if (teta3 < teta1 - pi_wp)
		  teta3 += 2.0*pi_wp;
	      }
	    else if (pt3_on_axis)
	      {
		if (teta1 > teta2 + pi_wp)
		  teta1 -= 2.0*pi_wp;
		if (teta1 < teta2 - pi_wp)
		  teta1 += 2.0*pi_wp;
	      }
	    else
	      {
		// modification of angles so that they are always close
		if (teta1 > teta2 + pi_wp)
		  teta1 -= 2.0*pi_wp;
		if (teta1 < teta2 - pi_wp)
		  teta1 += 2.0*pi_wp;
		
		if (teta3 > teta2 + pi_wp)
		  teta3 -= 2.0*pi_wp;
		if (teta3 < teta2 - pi_wp)
		  teta3 += 2.0*pi_wp;
		
	      }
	    
	    // computation of 2-D points
	    int order = this->lob_curve.GetOrder();
	    int Np = mesh_leaf.triangle_reference.GetNbPointsNodalElt();	
	    VectReal_wp xi(Np); VectR2 points2d; VectR3 points3d;
	    for (int j = 0; j < Np; j++)
	      {
		Real_wp x = mesh_leaf.triangle_reference.PointsNodalND(j)(0);
		Real_wp y = mesh_leaf.triangle_reference.PointsNodalND(j)(1);
		xi(j) = s1*(1.0-x-y) + s2*x + s3*y;
	      }
	    
	    GetPointsOnCurve(ns, ref, xi, points2d);
	    
	    // now computation of 3-D points
	    points3d.Reallocate(Np);
	    for (int j = 0; j < Np; j++)
	      {
		Real_wp x = mesh_leaf.triangle_reference.PointsNodalND(j)(0);
		Real_wp y = mesh_leaf.triangle_reference.PointsNodalND(j)(1);
		Real_wp teta = teta1*(1.0-x-y) + teta2*x + teta3*y;
		if (j > 2)
		  {
		    if (pt1_on_axis)
		      teta = (teta2*x + teta3*y)/(x+y);
		    else if (pt2_on_axis)
		      teta = (teta1*(1.0-x-y) + teta3*y)/(1.0-x);
		    else if (pt3_on_axis)
		      teta = (teta1*(1.0-x-y) + teta2*x)/(1.0-y);
		  }
		
		points3d(j).Init(points2d(j)(0)*cos(teta),
                                 points2d(j)(0)*sin(teta), points2d(j)(1));

		//DISP(j); DISP(x); DISP(y); DISP(teta);
		//DISP(teta1); DISP(teta2); DISP(teta3);
		//DISP(points2d(j)); DISP(points3d(j));
	      }
	    
	    // modification of edges
	    for (int j = 0; j < 3; j++)
	      {
		int num_edge = mesh_3d.BoundaryRef(i).numEdge(j);
		if (!EdgeUsed(num_edge))
		  {
		    EdgeUsed(num_edge) = true;
		    mesh_3d.EdgeRef(nb_edges_3d) = mesh_3d.GetEdge(num_edge);
		    mesh_3d.EdgeRef(nb_edges_3d).SetReference(ref);
		    if (mesh_3d.BoundaryRef(i).GetOrientationEdge(j))
		      for (int k = 0; k < order-1; k++)
			mesh_3d.SetPointInsideEdge(nb_edges_3d, k, points3d(3+j*(order-1)+k));
		    else
		      for (int k = 0; k < order-1; k++)
			mesh_3d.SetPointInsideEdge(nb_edges_3d, k,
                                                   points3d(3+j*(order-1)+order-2-k));
		    
		    nb_edges_3d++;
		  }
	      }
	    
	    // modification of internal points of triangles
	    for (int k = 3*order; k < Np; k++)
	      mesh_3d.SetPointInsideFace(i, k-3*order, points3d(k));
	  }
      }
    
    for (int i = 1; i < RefCurved.GetM(); i++)
      if (RefCurved(i))
	mesh_3d.SetCurveType(i, mesh_3d.CURVE_FILE);

    mesh_3d.ResizeEdgesRef(nb_edges_3d);
    mesh_3d.ConstructCrossReferenceEdgeRef();
    mesh_3d.SortBoundariesRef();
    mesh_3d.FindConnectivity();
  }  
#endif


  EllipseMinimizer::EllipseMinimizer(const Real_wp& x, const Real_wp& y, const R2& C,
				     const Real_wp& a0, const Real_wp& b0, const Real_wp& t0)
    : x0(x), y0(y), center(C), a(a0), b(b0), teta0(t0)
  {
    this->m_ = 2; this->n_ = 1;
  }


  void EllipseMinimizer::SetInitialAngle(const Real_wp& t0)
  {
    teta0 = t0;
  }

  
  void EllipseMinimizer::FindInitGuess(TinyVector<Real_wp, 1>& x)
  {
    x(0) = teta0;
  }


  void EllipseMinimizer::EvaluateF(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f)
  {
    Real_wp teta = x(0);
    f(0) = x0 - (center(0) + a*cos(teta));
    f(1) = y0 - (center(1) + b*sin(teta));
  }


  void EllipseMinimizer
  ::EvaluateJacobian(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f,
		     TinyMatrix<Real_wp, General, 2, 1>& df)
  {
    Real_wp teta = x(0);
    f(0) = x0 - (center(0) + a*cos(teta));
    f(1) = y0 - (center(1) + b*sin(teta));
    df(0, 0) = a*sin(teta);
    df(1, 0) = -b*cos(teta);
  }


  PeanutMinimizer::PeanutMinimizer(const Real_wp& x, const Real_wp& y, const R2& C,
				   const Real_wp& a0, const Real_wp& b0, const Real_wp& c0, const Real_wp& t0)
    : x0(x), y0(y), center(C), a(a0), b(b0), c(c0), teta0(t0)
  {
    this->m_ = 2; this->n_ = 1;
  }


  void PeanutMinimizer::SetInitialAngle(const Real_wp& t0)
  {
    teta0 = t0;
  }

  
  void PeanutMinimizer::FindInitGuess(TinyVector<Real_wp, 1>& x)
  {
    x(0) = teta0;
  }


  void PeanutMinimizer::EvaluateF(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f)
  {
    Real_wp teta = x(0);
    Real_wp sin_teta = sin(teta);
    f(0) = x0 - (center(0) + a*cos(teta)*sqrt(Real_wp(1) + c*sin_teta*sin_teta));
    f(1) = y0 - (center(1) + b*sin_teta);
  }


  void PeanutMinimizer
  ::EvaluateJacobian(const TinyVector<Real_wp, 1>& x, TinyVector<Real_wp, 2>& f,
		     TinyMatrix<Real_wp, General, 2, 1>& df)
  {
    Real_wp teta = x(0);
    Real_wp sin_teta = sin(teta), cos_teta = cos(teta);
    Real_wp radical = sqrt(Real_wp(1) + c*sin_teta*sin_teta);
    f(0) = x0 - (center(0) + a*cos_teta*radical);
    f(1) = y0 - (center(1) + b*sin_teta);
    df(0, 0) = a*sin_teta*(radical - cos_teta*c*cos_teta/radical);
    df(1, 0) = -b*cos(teta);
  }
    
}

#define MONTJOIE_FILE_MESH2D_BOUNDARIES_CXX
#endif
