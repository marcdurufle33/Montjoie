#ifndef MONTJOIE_FILE_MESH3D_CXX

namespace Montjoie
{

  //! copy constructor
  template<class T>
  Mesh<Dimension3, T>::Mesh(const Mesh<Dimension3, T>& mesh)
    : Mesh_Base<Dimension3, T>(mesh), MeshBoundaries<Dimension3, T>(mesh),
      faces(mesh.faces), faces_ref(mesh.faces_ref),
      next_face(mesh.next_face), head_faces_minv(mesh.head_faces_minv),
      tetrahedron_reference(mesh.tetrahedron_reference),
      pyramid_reference(mesh.pyramid_reference), wedge_reference(mesh.wedge_reference),
      hexahedron_reference(mesh.hexahedron_reference)
  {
    regular_mesh_with_slight_modification = mesh.regular_mesh_with_slight_modification;
    UpdateReferenceElement();
  }

  
  //! Copies a mesh
  template<class T>
  Mesh<Dimension3, T>& Mesh<Dimension3, T>
  ::operator=(const Mesh<Dimension3, T>& mesh)
  {
    static_cast<Mesh_Base<Dimension3, T>& >(*this)
      = static_cast<const Mesh_Base<Dimension3, T>& >(mesh);
    
    static_cast<MeshBoundaries<Dimension3, T>& >(*this)
      = static_cast<const MeshBoundaries<Dimension3, T>& >(mesh);

    faces = mesh.faces;
    faces_ref = mesh.faces_ref;
    next_face = mesh.next_face;
    head_faces_minv = mesh.head_faces_minv;
    
    tetrahedron_reference = mesh.tetrahedron_reference;
    pyramid_reference = mesh.pyramid_reference;
    wedge_reference = mesh.wedge_reference;
    hexahedron_reference = mesh.hexahedron_reference;

    regular_mesh_with_slight_modification = mesh.regular_mesh_with_slight_modification;

    UpdateReferenceElement();
    
    return *this;
  }

  
  //! updates pointers stored in reference_element
  template<class T>
  void Mesh<Dimension3, T>::UpdateReferenceElement()
  {
    for (int i = 0; i < this->reference_element.GetM(); i++)
      {
        int nb_vertices_loc = this->Element(i).GetNbVertices();
	switch (nb_vertices_loc)
	  {
	  case 4:
	    this->reference_element(i) = &tetrahedron_reference; break;
	  case 5:
	    this->reference_element(i) = &pyramid_reference; break;
	  case 6 :
	    this->reference_element(i) = &wedge_reference; break;
	  case 8 :
	    this->reference_element(i) = &hexahedron_reference; break;
	  }
      }
  }

  
  //! returns node numbers of the face num_loc
  template<class T>
  IVect Mesh<Dimension3, T>::GetLocalNodalNumber(int type, int num_loc)
  {
    IVect num;
    switch (type)
      {
      case 0 :
        {
          num.Reallocate(tetrahedron_reference.GetNbNodalBoundary(num_loc));
          for (int j = 0; j < num.GetM(); j++)
            num(j) = tetrahedron_reference.GetNodalNumber(num_loc, j);
        }
        break;
      case 1 :
        {
          num.Reallocate(pyramid_reference.GetNbNodalBoundary(num_loc));
          for (int j = 0; j < num.GetM(); j++)
            num(j) = pyramid_reference.GetNodalNumber(num_loc, j);
        }
        break;
      case 2 :
        {
          num.Reallocate(wedge_reference.GetNbNodalBoundary(num_loc));
          for (int j = 0; j < num.GetM(); j++)
            num(j) = wedge_reference.GetNodalNumber(num_loc, j);
        }
        break;
      case 3 :
        {
          num.Reallocate(hexahedron_reference.GetNbNodalBoundary(num_loc));
          for (int j = 0; j < num.GetM(); j++)
            num(j) = hexahedron_reference.GetNodalNumber(num_loc, j);
        }
        break;
      }
    
    return num;
  }
  
  
  //! returns the number of tetrahedra in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbTetrahedra() const
  {
    int nb_tetras = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->Element(i).GetNbVertices() == 4)
	nb_tetras++;
    
    return nb_tetras;
  }
  
  
  //! returns the number of pyramids in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbPyramids() const
  {
    int nb_pyramids = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->Element(i).GetNbVertices() == 5)
	nb_pyramids++;
    
    return nb_pyramids;
  }
  
  
  //! returns the number of triangular prisms in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbWedges() const
  {
    int nb_wedges = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->Element(i).GetNbVertices() == 6)
	nb_wedges++;
    
    return nb_wedges;
  }
  
  
  //! returns the number of hexahedra in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbHexahedra() const
  {
    int nb_hexas = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      if (this->Element(i).GetNbVertices() == 8)
	nb_hexas++;
    
    return nb_hexas;
  }
  
  
  //! returns the number of referenced triangles in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbTrianglesRef() const
  {
    int nb_tri = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      if (this->BoundaryRef(i).GetNbVertices() == 3)
	nb_tri++;
    
    return nb_tri;
  }


  //! returns the number of referenced quadrilaterals in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbQuadranglesRef() const
  {
    int nb_quad = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      if (this->BoundaryRef(i).GetNbVertices() == 4)
	nb_quad++;
    
    return nb_quad;
  }


  //! returns the number of triangles in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbTriangles() const
  {
    int nb_tri = 0;
    for (int i = 0; i < this->GetNbFaces(); i++)
      if (this->faces(i).GetNbVertices() == 3)
	nb_tri++;
    
    return nb_tri;
  }


  //! returns the number of quadrilaterals in the mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNbQuadrangles() const
  {
    int nb_quad = 0;
    for (int i = 0; i < this->GetNbFaces(); i++)
      if (this->faces(i).GetNbVertices() == 4)
	nb_quad++;
    
    return nb_quad;
  }


  //! we change the number of vertices, previous vertices are lost
  template<class T>
  void Mesh<Dimension3, T>::ReallocateVertices(int nb_vertices)
  {
    Mesh_Base<Dimension3, T>::ReallocateVertices(nb_vertices);
    head_faces_minv.Reallocate(nb_vertices); head_faces_minv.Fill(-1);    
  }
  
   
  //! we change the number of vertices, previous vertices are saved
  template<class T>
  void Mesh<Dimension3, T>::ResizeVertices(int nb_vertices)
  {
    int nb_old = this->GetNbVertices();
    Mesh_Base<Dimension3, T>::ResizeVertices(nb_vertices);
    head_faces_minv.Resize(nb_vertices);
    for (int i = nb_old; i < nb_vertices; i++)
      head_faces_minv(i) = -1;    
  }
  
  
  //! we change the number of faces, previous faces are lost
  template<class T>
  void Mesh<Dimension3, T>::ReallocateFaces(int nb_faces)
  {
    faces.Reallocate(nb_faces);
    next_face.Reallocate(nb_faces);
    head_faces_minv.Fill(-1);
  }
    
  
  //! we change the number of vertices, previous faces are saved
  template<class T>
  void Mesh<Dimension3, T>::ResizeFaces(int nb_faces)
  {
    int nb_old = this->GetNbFaces();
    faces.Resize(nb_faces);
    next_face.Resize(nb_faces);
    for (int i = nb_old; i < nb_faces; i++)
      next_face(i) = -1;
  }
  
  
  //! we change the number of referenced faces, previous referenced faces are lost
  template<class T>
  void Mesh<Dimension3, T>::ReallocateBoundariesRef(int nb_faces_ref)
  {
    faces_ref.Reallocate(nb_faces_ref);
    this->PointsFaceRef.Reallocate(nb_faces_ref);
    this->value_parameter_boundary.Reallocate(nb_faces_ref);
    this->integer_parameter_boundary.Reallocate(nb_faces_ref);
  }

  
  //! we change the number of referenced faces, previous referenced faces are saved
  template<class T>
  void Mesh<Dimension3, T>::ResizeBoundariesRef(int nb_faces_ref)
  {
    faces_ref.Resize(nb_faces_ref);
    this->PointsFaceRef.Resize(nb_faces_ref);
    this->value_parameter_boundary.Resize(nb_faces_ref);
    this->integer_parameter_boundary.Resize(nb_faces_ref);
  }
  
  
  //! we change the number of referenced edges, previous referenced edges are lost
  template<class T>
  void Mesh<Dimension3, T>::ReallocateEdgesRef(int nb_edges_ref)
  {
    this->edges_ref.Reallocate(nb_edges_ref);
    this->PointsEdgeRef.Reallocate(nb_edges_ref, this->GetGeometryOrder()-1);
    this->value_parameter_edge.Reallocate(nb_edges_ref);
    this->integer_parameter_edge.Reallocate(nb_edges_ref);
  }
  
  
  //! we change the number of referenced edges, previous referenced edges are kept
  template<class T>
  void Mesh<Dimension3, T>::ResizeEdgesRef(int nb_edges_ref)
  {
    this->edges_ref.Resize(nb_edges_ref);
    this->PointsEdgeRef.Resize(nb_edges_ref, this->GetGeometryOrder()-1);
    this->value_parameter_edge.Resize(nb_edges_ref);
    this->integer_parameter_edge.Resize(nb_edges_ref);
  }
  

  //! fills an array containing reference elements for elements
  template<class T>
  void Mesh<Dimension3, T>
  ::GetReferenceElementVolume(Vector<const ElementGeomReference<Dimension3>* >& ref) const
  {
    int nb_tet = 0, nb_pyr = 0, nb_prism = 0, nb_hex = 0;
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	switch (this->Element(i).GetNbVertices())
	  {
	  case 4 : nb_tet++; break;
	  case 5 : nb_pyr++; break;
	  case 6 : nb_prism++; break;
	  case 8 : nb_hex++; break;
	  }
      }

    int nb = 0;
    if (nb_tet > 0)
      nb++;

    if (nb_pyr > 0)
      nb++;

    if (nb_prism > 0)
      nb++;

    if (nb_hex > 0)
      nb++;
    
    ref.Reallocate(nb);
    nb = 0;
    if (nb_tet > 0)
      ref(nb++) = &this->tetrahedron_reference;

    if (nb_pyr > 0)
      ref(nb++) = &this->pyramid_reference;

    if (nb_prism > 0)
      ref(nb++) = &this->wedge_reference;

    if (nb_hex > 0)
      ref(nb++) = &this->hexahedron_reference;
  }

    
  //! modifies the order of approximation for the geometry
  /*!
    \param[in] r the new order of approximation
    An interpolation is performed to project points of curved faces to points such that
    they belong to polynomials of degrees up to r
    This method should be called before reading a curved mesh, otherwise the curved
    mesh will be read as a first-order mesh (i.e not curved), since the default order is 1.
    
    Mesh<Dimension3> mesh;
    
    mesh.SetGeometryOrder(4);
    mesh.Read("test.mesh");
    
    is a good way to proceed
   */
  template<class T>
  void Mesh<Dimension3, T>::SetGeometryOrder(int r, bool only_quad)
  {
    if (r <= 0)
      {
        cout << "you can't have a order less than 1" << endl;
        abort();
      }

    // if the order is the same, nothing to do
    if (this->lob_curve.GetOrder() == r)
      return;
    
    // constructing new finite elements
    Globatto<Real_wp> new_lob;
    TriangleGeomReference new_tri;
    QuadrangleGeomReference new_quad;
    new_lob.ConstructQuadrature(r, this->lob_curve.QUADRATURE_LOBATTO);
    if (!only_quad)
      new_tri.ConstructFiniteElement(r);    
    
    new_quad.ConstructFiniteElement(r);
    
    // specification of numbering scheme
    this->nodal_number_map.SetNbDofVertex(r, 1);
    this->nodal_number_map.SetNbDofEdge(r, r-1);
    this->nodal_number_map.SetNbDofTriangle(r, (r-1)*(r-2)/2);
    this->nodal_number_map.SetNbDofQuadrangle(r, (r-1)*(r-1));
    this->nodal_number_map.SetNbDofTetrahedron(r, (r-3)*(r-2)*(r-1)/6);
    this->nodal_number_map.SetNbDofPyramid(r, (r-1)*(r-2)*(2*r-3)/6);
    this->nodal_number_map.SetNbDofWedge(r, (r-1)*(r-2)*(r-1)/2);
    this->nodal_number_map.SetNbDofHexahedron(r, (r-1)*(r-1)*(r-1));
    this->nodal_number_map.SetOppositeEdgesDofSymmetry(r, r-1);        
    
    VectReal_wp xi;
    if (r > 1)  
      {
	xi.Reallocate(r-1);
	for (int i = 1; i < r; i++)
	  xi(i-1) = new_lob.Points(i);
      }

    // we try to recompute PointsEdgeRef with the new order of approximation provided
    Matrix<R_N> NewPointsEdge;
    if (this->GetNbEdgesRef() > 0)
      {
	if (r <= 1)
	  NewPointsEdge.Clear();
	else
	  {
	    if (this->lob_curve.GetOrder() <= 1)
	      {
		NewPointsEdge.Reallocate(this->GetNbEdgesRef(), r-1);
		// basic linear interpolation		
		for (int i = 0; i < this->GetNbEdgesRef(); i++)
		  {
		    int ref = this->EdgeRef(i).GetReference();		    
		    if ((this->GetCurveType(ref) > 0)
                        &&(this->GetCurveType(ref) != this->CURVE_FILE))
		      {
			Vector<R_N> points_interm;
			this->GetPointsOnCurvedEdge(i, ref, xi, points_interm);
			for (int k = 0; k < xi.GetM(); k++)
			  NewPointsEdge(i, k) = points_interm(k);
		      }
		    else
		      {
			R_N ptA = this->Vertex(this->EdgeRef(i).numVertex(0));
			R_N ptB = this->Vertex(this->EdgeRef(i).numVertex(1));
			for (int k = 0; k < xi.GetM(); k++)
			  {
			    NewPointsEdge(i, k).Zero();
			    Add(1.0-xi(k), ptA, NewPointsEdge(i, k));
			    Add(xi(k), ptB, NewPointsEdge(i, k));
			  }
		      }
		  }
	      }
	    else
	      {
		int m = this->lob_curve.GetOrder();
		Matrix<Real_wp> ValuePhi(m+1, r+1);
		for (int j = 0; j <= m; j++)
		  for (int k = 1; k < r; k++)
		    ValuePhi(j, k-1) = this->lob_curve.EvaluatePhi(j, xi(k-1));
				
		NewPointsEdge.Reallocate(this->GetNbEdgesRef(), r-1);
		// interpolation with previous points
		for (int i = 0; i < this->GetNbEdgesRef(); i++)
		  {
		    int ref = this->EdgeRef(i).GetReference();		    
		    if ((this->GetCurveType(ref) > 0)
                        &&(this->GetCurveType(ref) != this->CURVE_FILE))
		      {
			Vector<R_N> points_interm;
			this->GetPointsOnCurvedEdge(i, ref, xi, points_interm);
			for (int k = 0; k < xi.GetM(); k++)
			  NewPointsEdge(i, k) = points_interm(k);
		      }
		    else
		      {
			VectR3 OldPoints(m+1);
			OldPoints(0) = this->Vertex(this->EdgeRef(i).numVertex(0));
			OldPoints(m) = this->Vertex(this->EdgeRef(i).numVertex(1));
			for (int k = 0; k < m-1; k++)
			  OldPoints(k+1) = this->PointsEdgeRef(i, k);
						
			for (int k = 0; k < xi.GetM(); k++)
			  {
			    NewPointsEdge(i, k).Zero();
			    for (int j = 0; j <= m; j++)
			      Add(ValuePhi(j, k), OldPoints(j), NewPointsEdge(i, k));
			  }
		      }
		  }
	      }
	  }
      }
    
    VectR2 xi_tri, xi_quad;
    int Ntri = new_tri.GetNbPointsNodalElt() - 3*r;
    int Nquad = new_quad.GetNbPointsNodalElt() - 4*r;
    if (only_quad)
      Ntri = 0;
    
    if (r > 1)  
      {
	Matrix<int> FacesDof_Rotation_Tri;

	if (!only_quad)
          {
            xi_tri.Reallocate(Ntri);
            for (int i = 0; i < Ntri; i++)
              xi_tri(i) = new_tri.PointsNodalND(i+3*r);

            MeshNumbering<Dimension3, T>::
              GetRotationTriangularFace(xi_tri, FacesDof_Rotation_Tri);
		
            this->nodal_number_map.SetFacesDofRotationTri(r, FacesDof_Rotation_Tri);
          }
        
	xi_quad.Reallocate(Nquad);
	for (int i = 0; i < Nquad; i++)
	  xi_quad(i) = new_quad.PointsNodalND(4*r+i);
	
	Matrix<int> NumNodes2D_InsideQuad(r-1, r-1);
	for (int i = 1; i < r; i++)
	  for (int j = 1; j < r; j++)
	    NumNodes2D_InsideQuad(i-1, j-1) = new_quad.GetNumNodes2D(i, j) - 4*r;
	
	Matrix<int> FacesDof_Rotation_Quad;
	MeshNumbering<Dimension3, T>::
	  GetRotationQuadrilateralFace(NumNodes2D_InsideQuad, FacesDof_Rotation_Quad);
	
	this->nodal_number_map.SetFacesDofRotationQuad(r, FacesDof_Rotation_Quad);
      }
    
    // we try to recompute PointsFaceRef with the new order of approximation provided
    Vector<Vector<R_N> > NewPointsFace;
    if (this->GetNbBoundaryRef() > 0)
      {
	if (r <= 1)
	  NewPointsFace.Clear();
	else
	  {
            NewPointsFace.Reallocate(this->GetNbBoundaryRef());
	    if (this->lob_curve.GetOrder() <= 1)
	      {		
		// basic linear interpolation		
		for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		  {
		    int ref = this->BoundaryRef(i).GetReference();		    
		    bool tri = (this->BoundaryRef(i).GetNbVertices() == 3);
		    if ((this->GetCurveType(ref) > 0)
                        &&(this->GetCurveType(ref) != this->CURVE_FILE))
		      {
			Vector<R_N> points_interm;
			if (tri)
			  {
			    NewPointsFace(i).Reallocate(Ntri);
			    this->GetPointsOnCurvedFace(i, ref, xi_tri, points_interm);
			    for (int k = 0; k < xi_tri.GetM(); k++)
			      NewPointsFace(i)(k) = points_interm(k);
			  }
			else
			  {
			    NewPointsFace(i).Reallocate(Nquad);
			    this->GetPointsOnCurvedFace(i, ref, xi_quad, points_interm);
			    for (int k = 0; k < xi_quad.GetM(); k++)
			      NewPointsFace(i)(k) = points_interm(k);
			  }
		      }
		    else
		      {
			R_N ptA = this->Vertex(this->BoundaryRef(i).numVertex(0));
			R_N ptB = this->Vertex(this->BoundaryRef(i).numVertex(1));
			R_N ptC = this->Vertex(this->BoundaryRef(i).numVertex(2));
			if (tri)
			  {
			    NewPointsFace(i).Reallocate(Ntri);
			    for (int k = 0; k < xi_tri.GetM(); k++)
			      {
				NewPointsFace(i)(k).Zero();
				Add(1.0-xi_tri(k)(0)-xi_tri(k)(1), ptA, NewPointsFace(i)(k));
				Add(xi_tri(k)(0), ptB, NewPointsFace(i)(k));
				Add(xi_tri(k)(1), ptC, NewPointsFace(i)(k));
			      }
			  }
			else
			  {
			    R_N ptD = this->Vertex(this->BoundaryRef(i).numVertex(3));
			    NewPointsFace(i).Reallocate(Nquad);
			    for (int k = 0; k < xi_quad.GetM(); k++)
			      {
				NewPointsFace(i)(k).Zero();
				Add((1.0-xi_quad(k)(0))*(1.0-xi_quad(k)(1)),
                                    ptA, NewPointsFace(i)(k));
				
                                Add(xi_quad(k)(0)*(1.0-xi_quad(k)(1)), ptB, NewPointsFace(i)(k));
				Add(xi_quad(k)(0)*xi_quad(k)(1), ptC, NewPointsFace(i)(k));
				Add((1.0-xi_quad(k)(0))*xi_quad(k)(1), ptD, NewPointsFace(i)(k));
			      }
			  }
		      }
		  }
	      }
	    else
	      {
		Matrix<Real_wp> ValuePhiTri;
                if (!only_quad)
                  ValuePhiTri.Reallocate(this->triangle_reference.GetNbPointsNodalElt(), Ntri);
                
		Matrix<Real_wp> ValuePhiQuad(this->quadrangle_reference.GetNbPointsNodalElt(),
                                             Nquad);
		VectReal_wp phi;
		for (int k = 0; k < Ntri; k++)
		  {
		    this->triangle_reference.ComputeValuesPhiNodalRef(xi_tri(k), phi);
		    for (int j = 0; j < this->triangle_reference.GetNbPointsNodalElt(); j++)
		      ValuePhiTri(j, k) = phi(j);
		  }
		
		for (int k = 0; k < Nquad; k++)
		  {
		    this->quadrangle_reference.ComputeValuesPhiNodalRef(xi_quad(k), phi);
		    for (int j = 0; j < this->quadrangle_reference.GetNbPointsNodalElt(); j++)
		      ValuePhiQuad(j, k) = phi(j);
		  }
		
		// interpolation with previous points
		for (int i = 0; i < this->GetNbBoundaryRef(); i++)
		  {
		    int ref = this->BoundaryRef(i).GetReference();		    
		    bool tri = (this->BoundaryRef(i).GetNbVertices() == 3);
		    if ((this->GetCurveType(ref) > 0)
                        &&(this->GetCurveType(ref) != this->CURVE_FILE))
		      {
			Vector<R_N> points_interm;
			if (tri)
			  {
			    NewPointsFace(i).Reallocate(Ntri);
			    this->GetPointsOnCurvedFace(i, ref, xi_tri, points_interm);
			    for (int k = 0; k < xi_tri.GetM(); k++)
			      NewPointsFace(i)(k) = points_interm(k);
			  }
			else
			  {
			    NewPointsFace(i).Reallocate(Nquad);
			    this->GetPointsOnCurvedFace(i, ref, xi_quad, points_interm);
			    for (int k = 0; k < xi_quad.GetM(); k++)
			      NewPointsFace(i)(k) = points_interm(k);
			  }
		      }
		    else
		      {
			if (tri)
			  {
			    VectR3 OldPoints;
                            this->GetNodalPointsOnFace(i, OldPoints);
                            
			    NewPointsFace(i).Reallocate(Ntri);
			    for (int k = 0; k < xi_tri.GetM(); k++)
			      {
				NewPointsFace(i)(k).Zero();
				for (int j = 0; j < OldPoints.GetM(); j++)
				  Add(ValuePhiTri(j, k), OldPoints(j), NewPointsFace(i)(k));
			      }
			  }
			else
			  {
			    VectR3 OldPoints;
                            this->GetNodalPointsOnFace(i, OldPoints);
                            
			    NewPointsFace(i).Reallocate(Nquad);
			    for (int k = 0; k < xi_quad.GetM(); k++)
			      {
				NewPointsFace(i)(k).Zero();
				for (int j = 0; j < OldPoints.GetM(); j++)
				  Add(ValuePhiQuad(j, k), OldPoints(j), NewPointsFace(i)(k));
			      }
			  }
		      }
		  }
	      }
	  }
      }
        
    this->lob_curve = new_lob;
    MeshNumbering<Dimension3> mesh_num(*this);
    if (!only_quad)
      new_tri.SetMesh(*this);
    
    new_quad.SetMesh(*this);
    
    if (!only_quad)
      this->triangle_reference = new_tri;
    
    this->quadrangle_reference = new_quad;
    this->PointsEdgeRef = NewPointsEdge;
    this->PointsFaceRef = NewPointsFace;

    if (!only_quad)
      {
        tetrahedron_reference.ConstructFiniteElement(r);
        tetrahedron_reference.SetMesh(*this);
        
        wedge_reference.ConstructFiniteElement(r);
        wedge_reference.SetMesh(*this);
        
        pyramid_reference.ConstructFiniteElement(r);
        pyramid_reference.SetMesh(*this);
      }
    
    hexahedron_reference.ConstructFiniteElement(r);
    hexahedron_reference.SetMesh(*this);    
  }
  
  
  /************************
   * Convenient Functions *
   ************************/
  

  //! returns true if the mesh is hybrid
  //! (containing a mix of hexahedra, tetrahedra, prisms and pyramids)
  template<class T>
  bool Mesh<Dimension3, T>::HybridMesh() const
  {
    int nb_tri = GetNbTriangles();
    int nb_quad = GetNbQuadrangles();
    if ((nb_tri >0) && (nb_quad > 0))
      return true;
    
    return false;
  }

  
  //! returns true if the mesh is made only of tetrahedra
  template<class T>
  bool Mesh<Dimension3, T>::IsOnlyTetrahedral() const
  {
    int nb_pyramids = GetNbPyramids();
    int nb_wedges = GetNbWedges();
    int nb_hexas = GetNbHexahedra();
    if ((nb_hexas > 0)||(nb_pyramids > 0)||(nb_wedges > 0))
      return false;
    
    return true;
  }
  
  
  //! returns true if the mesh is made only of hexahedra
  template<class T>
  bool Mesh<Dimension3, T>::IsOnlyHexahedral() const
  {
    int nb_pyramids = GetNbPyramids();
    int nb_wedges = GetNbWedges();
    int nb_tetras = GetNbTetrahedra();
    if ((nb_tetras > 0)||(nb_pyramids > 0)||(nb_wedges > 0))
      return false;
    
    return true;
  }
  
  
  //! returns size of memory used by the object in bytes
  template<class T>
  size_t Mesh<Dimension3, T>::GetMemorySize() const
  {
    size_t taille = Mesh_Base<Dimension3, T>::GetMemorySize();
    taille += MeshBoundaries<Dimension3, T>::GetMemorySize() + sizeof(*this);
    
    taille += sizeof(int)*(next_face.GetM() + head_faces_minv.GetM());
    for (int i = 0; i < faces.GetM(); i++)
      taille += faces(i).GetMemorySize();

    for (int i = 0; i < faces_ref.GetM(); i++)
      taille += faces_ref(i).GetMemorySize();
    
    taille += tetrahedron_reference.GetMemorySize();
    taille += pyramid_reference.GetMemorySize();
    taille += wedge_reference.GetMemorySize();
    taille += hexahedron_reference.GetMemorySize();
    return taille;
  }
  
  
  //! returns a bounding number of edges, if number of elements is known
  template<class T>
  int Mesh<Dimension3, T>::GetUpperBoundNumberOfEdges() const
  {
    int nb_hexas = GetNbHexahedra();
    int nb_tetras = GetNbTetrahedra();
    int nb_pyramids = GetNbPyramids();
    int nb_wedges = GetNbWedges();
    return (12*nb_hexas + 6*nb_tetras + 8*nb_pyramids + 9*nb_wedges);
    // theorical asymptote is 3 for hexahedras and 1.1755 for tetrahedras
    // a bounding number is 12 for hexahedras and 6 for tetrahedras
  }
  
  
  //! returns a bounding number of faces
  //! knowning only number of hexahedra, tetrahedra and referenced faces
  template<class T>
  int Mesh<Dimension3, T>::GetUpperBoundNumberOfFaces() const
  {
    int nb_hexas = GetNbHexahedra();
    int nb_tetras = GetNbTetrahedra();
    int nb_pyramids = GetNbPyramids();
    int nb_wedges = GetNbWedges();
    if (this->GetNbElt()<40)
      return (6*nb_hexas+4*nb_tetras+5*(nb_pyramids+nb_wedges));
    else if (this->GetNbElt() < 1000)
      return (6*nb_hexas + 4*nb_tetras+5*(nb_pyramids+nb_wedges));
    
    // theorical asymptote is 3 for hexahedras and 1.1755 for tetrahedras
    // a bounding number is 6 for hexahedras and 4 for tetrahedras 
    
    return (4*nb_hexas + 3*nb_tetras + 3*nb_pyramids + 3*nb_wedges + GetNbBoundaryRef());
  }
  
  
  //! modifies parameters of the mesh with a line of the data file
  template<class T>
  void Mesh<Dimension3, T>
  ::SetInputData(const string& description_field, const Vector<string>& parameters)
  {
    Mesh_Base<Dimension3, T>::SetInputData(description_field, parameters);
    MeshBoundaries<Dimension3, T>::SetInputData(description_field, parameters);
    if (!description_field.compare("TypeMesh"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh3D" << endl;
	    cout << "TypeMesh needs more parameters, for instance :" << endl;
	    cout << "TypeMesh = HEXA" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	this->required_mesh_type = Mesh<Dimension3>::TETRAHEDRAL_MESH;
	if (!parameters(0).compare("TETRA"))
	  this->required_mesh_type = Mesh<Dimension3>::TETRAHEDRAL_MESH;
	else if (!parameters(0).compare("HEXA"))
	  this->required_mesh_type = Mesh<Dimension3>::HEXAHEDRAL_MESH;
	else if (!parameters(0).compare("PYRAMID"))
	  this->required_mesh_type = Mesh<Dimension3>::PYRAMID_MESH;
	else if (!parameters(0).compare("WEDGE"))
	  this->required_mesh_type = Mesh<Dimension3>::WEDGE_MESH;
	else if (!parameters(0).compare("HYBRID"))
	  this->required_mesh_type = Mesh<Dimension3>::HYBRID_MESH;
	else if (!parameters(0).compare("CRUSH"))
	  this->required_mesh_type = Mesh<Dimension3>::CRUSH_MESH;    
        // degenerated hexahedra into pyramids
	
      }
    else if (!description_field.compare("SlightModificationOnRegularMesh"))
      {       
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of Mesh3D" << endl;
	    cout << "SlightModificationOnRegularMesh needs more parameters, for instance :" 
                 << endl;
	    cout << "SlightModificationOnRegularMesh = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	regular_mesh_with_slight_modification = REGULAR_MESH_NOT_MODIFIED;
	if (!parameters(0).compare("YES"))
	  regular_mesh_with_slight_modification = REGULAR_MESH_SLIGHTLY_MODIFIED;
        else if (!parameters(0).compare("RANDOM"))
          regular_mesh_with_slight_modification = REGULAR_MESH_RANDOMLY_MODIFIED;
      }
  }
  

  //! constructs arrays useful to write high-order curved mesh
  template<class T>
  int Mesh<Dimension3, T>::GetNodesCurvedMesh(Vector<R_N>& PosNodes, Vector<IVect>& Nodle,
					      int rmax, bool lobatto) const
  {
    Vector<IVect> NodleSurf;
    return GetNodesCurvedMesh(PosNodes, Nodle, NodleSurf, rmax, lobatto);
  }
  

  //! constructs arrays useful to write high-order curved mesh
  /*!
    \param[out] PosNodes all the nodes of the mesh (even internal nodes due to high-order)
    \param[out] Nodle Numbering of all the nodes (as for continuous finite elements)
    \param[in] rmax maximal allowed order (if you want only nodes for first order mesh, put 1)
  */
  template<class T>
  int Mesh<Dimension3, T>::GetNodesCurvedMesh(Vector<R_N>& PosNodes, Vector<IVect>& Nodle,
					      Vector<IVect>& NodleSurf, int rmax, bool lobatto) const
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

	NodleSurf.Reallocate(this->GetNbBoundaryRef());
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  {
	    int nv = this->BoundaryRef(i).GetNbVertices();
	    NodleSurf(i).Reallocate(nv);
	    for (int j = 0; j < nv; j++)
	      NodleSurf(i)(j) = this->BoundaryRef(i).numVertex(j);
	  }
      }
    else
      {
	VectReal_wp step_x(order+1);
	
	Vector<VectR2> points_surf(2);
	points_surf(0).Reallocate((order+1)*(order+2)/2);
	points_surf(1).Reallocate((order+1)*(order+1));
	
	Matrix<int> NumNodes2D_tri, NumNodes2D_quad, coor;
	if (lobatto)
          {
	    // Lobatto points on triangle
            TriangleGeomReference::
              ConstructLobattoPoints(order, 0, step_x, points_surf(0));
          }
	else
	  {
	    for (int i = 0; i <= order; i++)
	      step_x(i) = Real_wp(i)/order;
	    
	    // regular points on triangles
	    MeshNumbering<Dimension2, T>::ConstructTriangularNumbering(order, NumNodes2D_tri, coor);
	    for (int i = 0; i <= order; i++)
	      for (int j = 0; j <= order-i; j++)
		points_surf(0)(NumNodes2D_tri(i, j)).Init(step_x(i), step_x(j));
	  }
	
	// tensorized points on quadrangles
	MeshNumbering<Dimension2, T>::ConstructQuadrilateralNumbering(order, NumNodes2D_quad, coor);
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    points_surf(1)(NumNodes2D_quad(i, j)).Init(step_x(i), step_x(j));
	
	Vector<VectR3> points_div(4);
	points_div(0).Reallocate((order+1)*(order+2)*(order+3)/6);
	points_div(1).Reallocate((order+1)*(order+2)*(2*order+3)/6);
	points_div(2).Reallocate((order+1)*(order+2)*(order+1)/2);
	points_div(3).Reallocate((order+1)*(order+1)*(order+1));
	
	Array3D<int> NumNodes3D;
	if (lobatto)
	  {
	    VectReal_wp points1d; VectR2 points2d;
	    TetrahedronGeomReference::
	      ConstructLobattoPoints(order, points1d, points2d, points_div(0));
	  }
	else
	  {
	    // regular points on tetrahedron
	    MeshNumbering<Dimension3, T>::ConstructTetrahedralNumbering(order, NumNodes3D, coor);
	    for (int i = 0; i <= order; i++)
	      for (int j = 0; j <= order-i; j++)
		for (int k = 0; k <= order-i-j; k++)
		  points_div(0)(NumNodes3D(i, j, k)).Init(step_x(i), step_x(j), step_x(k));
	  }
	
	if (lobatto)
	  {
	    VectReal_wp points1d; VectR2 points2d, points2d_q;
	    PyramidGeomReference::
	      ConstructLobattoPoints(order, points1d, points2d,
				     points2d_q, NumNodes2D_quad, points_div(1));
	  }
	else
	  {
	    // regular points on pyramid
	    MeshNumbering<Dimension3, T>::ConstructPyramidalNumbering(order, NumNodes3D, coor);
	    for (int k = 0; k <= order; k++)
	      for (int i = 0; i <= order-k; i++)
		for (int j = 0; j <= order-k; j++)
		  {
		    Real_wp z = step_x(k);
		    Real_wp x(0), y(0);
		    if (k < order)
		      {
			x = step_x(i)/step_x(order-k);
			y = step_x(j)/step_x(order-k);
		      }
		    
		    points_div(1)(NumNodes3D(i, j, k))
		      .Init((1.0-z)*(2.0*x-1.0), (1.0-z)*(2.0*y-1.0), z);
		  }
	  }
	
	if (lobatto)
	  {
	    VectReal_wp points1d; VectR2 points2d, points2d_q;
	    MeshNumbering<Dimension2, T>::ConstructTriangularNumbering(order, NumNodes2D_tri, coor);
	    MeshNumbering<Dimension3, T>::ConstructPrismaticNumbering(order, NumNodes3D, coor);
	    WedgeGeomReference::
	      ConstructLobattoPoints(order, points1d, points2d, NumNodes2D_tri,
				     points2d_q, NumNodes2D_quad, points_div(2), NumNodes3D);	     
	  }
	else
	  {
	    // regular points on prism
	    MeshNumbering<Dimension3, T>::ConstructPrismaticNumbering(order, NumNodes3D, coor);
	    for (int i = 0; i <= order; i++)
	      for (int j = 0; j <= order-i; j++)
		for (int k = 0; k <= order; k++)
		  points_div(2)(NumNodes3D(i, j, k)).Init(step_x(i), step_x(j), step_x(k));
	  }
	
	// tensorized points on hexahedron
	MeshNumbering<Dimension3, T>::ConstructHexahedralNumbering(order, NumNodes3D, coor);
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      points_div(3)(NumNodes3D(i, j, k)).Init(step_x(i), step_x(j), step_x(k));
	
	this->GetMeshSubdivision(step_x, points_div, points_surf, PosNodes, Nodle, NodleSurf);
      }
    
    return order;
  }

  
  //! reading of a mesh, several formats are available
  /*!
    \param[in] file_mesh file name where the mesh is stored
  */
  template<class T>
  void Mesh<Dimension3, T>::Read(const string & file_mesh) 
  { 
    Clear();
    string extension_input = GetExtension(file_mesh);
    
    // mesh file is opened
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
  void Mesh<Dimension3, T>::Read(istream& file_in, const string& extension_input) 
  {
    int new_ref = -1;
    
    if (!extension_input.compare("mit"))
      {
	int nb_vertices, nb_faces_ref;
	int itmp = 0;
	file_in>>nb_vertices>>nb_faces_ref>>itmp;
	this->ReallocateVertices(nb_vertices);
	// lecture des sommets
	for (int i = 0; i < nb_vertices; i++)
	  file_in>>this->Vertices(i)(0)>>this->Vertices(i)(1)>>this->Vertices(i)(2);
	
	string ligne;
	Vector<string> param;
	
	// lecture des faces de references
	getline(file_in, ligne);
	ReallocateBoundariesRef(nb_faces_ref);
	IVect num(4);
	for (int i = 0; i < nb_faces_ref; i++)
	  {
	    getline(file_in, ligne);
	    StringTokenize(ligne, param, " ");
	    int nb_vertices_loc = param.GetM();
	    for (int j = 0; j < nb_vertices_loc; j++)
	      num(j) = to_num<int>(param(j));
	    
	    if (nb_vertices_loc == 3)
	      faces_ref(i).InitTriangular(num(0), num(1), num(2), 0);
	    else
	      faces_ref(i).InitQuadrangular(num(0), num(1), num(2), num(3), 0);
	  }
	
	// affectation des faces de reference
	for (int i = 0; i < nb_faces_ref; i++)
	  file_in>>itmp;
	
	for (int i = 0; i < nb_faces_ref; i++)
	  {
	    file_in>>itmp;	                
            faces_ref(i).SetReference(itmp);
	  }
	
      }
    else if (!extension_input.compare("gts"))
      {
	
	// on tokenize chaque ligne
	Vector<string> parameters;
	
	int num_ligne = 0, type_champ = 0; int n0, n1, n2, nb; string ligne;
	IVect num(6); int nb_edges = 0, nb_vertices = 0, nb_faces_ref = 0;
	while (!file_in.eof())
	  {
	    parameters.Clear();
	    getline(file_in, ligne);
	    // si une ligne commence avec #, c'est un commentaire 
	    if ((ligne.size() > 0) && (ligne[0] !='#'))
	      {
		StringTokenize(ligne, parameters);
		
		if (type_champ == 0)
		  {
		    // premiere ligne, nombre sommets, aretes et triangles
		    nb_vertices = to_num<int>(parameters(0));
		    nb_edges = to_num<int>(parameters(1));
		    nb_faces_ref = to_num<int>(parameters(2));
		    // allocation
		    this->ReallocateVertices(nb_vertices);
		    this->ReallocateEdges(nb_edges);
		    this->ReallocateBoundariesRef(nb_faces_ref);
		    type_champ++;
		  }
		else if (type_champ == 1)
		  {
		    // sommet
		    this->Vertices(num_ligne)(0) = to_num<T>(parameters(0));
		    this->Vertices(num_ligne)(1) = to_num<T>(parameters(1));
		    this->Vertices(num_ligne)(2) = to_num<T>(parameters(2));
		    num_ligne++;
		    if (num_ligne == nb_vertices)
		      {
			type_champ++;
			num_ligne = 0;
		      }
		  }
		else if (type_champ == 2)
		  {
		    // arete
		    n0 = to_num<int>(parameters(0));
		    n1 = to_num<int>(parameters(1));
		    this->edges(num_ligne).Init(n0-1, n1-1, 0);
		    num_ligne++;
		    if (num_ligne == nb_edges)
		      {
			type_champ++;
			num_ligne = 0;
		      }
		  }
		else if (type_champ == 3)
		  {
		    // triangle
		    n0 = to_num<int>(parameters(0))-1;
		    n1 = to_num<int>(parameters(1))-1;
		    n2 = to_num<int>(parameters(2))-1;
		    num(0) = this->edges(n0).numVertex(0); num(1) = this->edges(n0).numVertex(1);
		    num(2) = this->edges(n1).numVertex(0); num(3) = this->edges(n1).numVertex(1);
		    num(4) = this->edges(n2).numVertex(0); num(5) = this->edges(n2).numVertex(1);
		    nb = 6; Assemble(nb, num);
		    faces_ref(num_ligne).InitTriangular(num(0), num(1), num(2), 1);     
		    num_ligne++;
		    if (num_ligne == nb_faces_ref)
		      type_champ++;
		  }
	      }
	  }
      }
    else if (!extension_input.compare("vf"))
      {
	int itmp = 0;
	int nb_vertices = 0, nb_faces = 0;
	file_in>>nb_vertices;
	this->ReallocateVertices(nb_vertices);
	// lecture des sommets
	for (int i = 0; i < nb_vertices; i++)
	  file_in>>this->Vertices(i)(0)>>this->Vertices(i)(1)>>this->Vertices(i)(2);
	
	file_in>>nb_faces;
	
	// lecture des faces
	faces.Reallocate(nb_faces);
	IVect num(20);
	for (int i = 0; i < nb_faces; i++)
	  {
	    int nb_vertices_loc;
	    file_in>>nb_vertices_loc;
	    for (int j = 0; j < nb_vertices_loc; j++)
	      {
		file_in>>num(j);
		num(j)--;
	      }
	    file_in>>itmp;
	    
	    if (nb_vertices_loc == 3)
	      faces(i).InitTriangular(num(0), num(1), num(2), itmp);
	    else
	      faces(i).InitQuadrangular(num(0), num(1), num(2), num(3), itmp);
	  }
	
	// extraction faces references
	int nb_faces_ref = 0;
	for (int i = 0; i < nb_faces; i++)
	  if ((faces(i).GetReference() != 0)&&(faces(i).GetReference() != 0))
	    nb_faces_ref++;
	
	ReallocateBoundariesRef(nb_faces_ref);
	nb_faces_ref = 0;
	for (int i = 0; i < nb_faces; i++)
	  if ((faces(i).GetReference() != 0)&&(faces(i).GetReference() != 0))
	    {
	      faces_ref(nb_faces_ref) = faces(i);
	      nb_faces_ref++;
	    }
	
	// lecture tetraedres/hexahedres
	int nb_elt = 0;
	file_in>>nb_elt;
	
	this->ReallocateElements(nb_elt);
	IVect IndexVertexOnFace(nb_vertices); IndexVertexOnFace.Fill(-1);
	for (int i = 0; i < nb_elt; i++)
	  {
	    int nb_faces_loc, ref;
	    file_in >> nb_faces_loc;
	    for (int j = 0; j < nb_faces_loc; j++)
	      {
		file_in>>num(j);
		num(j)--;
	      }
	    file_in>>ref>>itmp;
	    
	    if (nb_faces_loc == 4)
	      {
		IVect num_vertex(12);
		int nb = 0;
		// liste des sommets
		for (int j = 0; j < nb_faces_loc; j++)
		  for (int k = 0; k < 3; k++)
		    num_vertex(nb++) = faces(num(j)).numVertex(k);
		
		// on trie, on elimine les doublons
		nb = 12;
		Assemble(nb, num_vertex); num_vertex.Resize(nb);
		
		this->elements(i).InitTetrahedral(num_vertex, itmp);
	      }
	    else
	      {
		IVect num_vertex(8);
		// quatre premiers sommets de la premiere face
		for (int k = 0; k < 4; k++)
		  {
		    num_vertex(k) = faces(num(0)).numVertex(k);
		    IndexVertexOnFace(num_vertex(k)) = k;
		  }
		
		// sommets opposes ?
		// boucle sur les autres faces
		for (int j = 1; j < nb_faces_loc; j++)
		  // if (j != num_face_opp)
		  {
		    for (int k = 0; k < 4; k++)
		      {
			int n0 = faces(num(j)).numVertex(k);
			int n1 = faces(num(j)).numVertex((k+1)%4);
			int pos0 = IndexVertexOnFace(n0);
			int pos1 = IndexVertexOnFace(n1);
			if ((pos0 != -1)&& (pos1 == -1))
			  num_vertex(pos0+4) = n1;
			
			if ((pos0 == -1)&& (pos1 != -1))
			  num_vertex(pos1+4) = n0;
		      }
		  }
		
		for (int k = 0; k < 4; k++)
		  IndexVertexOnFace(num_vertex(k)) = -1;
		
		this->elements(i).InitHexahedral(num_vertex, ref);
	      }
	  }
	
      }
    else if (!extension_input.compare("stl"))
      {
	string ligne;
	getline(file_in, ligne);
	std::list<R3> all_vertices; R3 point;
	while (!file_in.eof())
	  {
	    getline(file_in, ligne);
	    if (ligne.find("facet ") != string::npos)
	      {
		getline(file_in, ligne);
		for (int k = 0; k < 3; k++)
		  {
		    file_in>>ligne>>point(0)>>point(1)>>point(2);
		    all_vertices.push_back(point);
		  }
	      }
	  }
	
	this->ReallocateVertices(all_vertices.size());
	list<R3>::iterator iter_pt; int nb_vertices = 0;
	for (iter_pt = all_vertices.begin(); iter_pt != all_vertices.end(); iter_pt++)
	  this->Vertices(nb_vertices++) = *iter_pt;
	
	// faces
	int nb_faces_ref = all_vertices.size()/3;
	ReallocateBoundariesRef(nb_faces_ref);
	for (int i = 0; i < nb_faces_ref; i++)
	  faces_ref(i).InitTriangular(3*i, 3*i+1, 3*i+2, 1);
	
	// elimination of duplicate vertices
	this->RemoveDuplicateVertices();
      }
    else if (!extension_input.compare("nas"))
      {
	// on lit toutes les lignes
	string ligne; Vector<string> mots;
	std::list<R3> all_vertices; T x, y, z;
	std::list<TinyVector<int, 5> > all_faces;
	std::list<TinyVector<int, 9> > all_elt;
	TinyVector<int, 9> elt;
	while (!file_in.eof())
	  {
	    getline(file_in, ligne);
	    if (!file_in.good())
	      break;
	    
	    StringTokenize(ligne, mots, " \t");
	    //cout << mots(0) << endl;
	    //int test_input; cout << "waiting" << endl; cin >> test_input;
	    if (mots(0) == "GRID")
	      {
		// a vertex
		x = to_num<T>(mots(2).substr(0, 8));
		y = to_num<T>(mots(2).substr(8, 8));
		z = to_num<T>(mots(2).substr(16, 8));
		all_vertices.push_back(R3(x, y, z));
	      }
	    else if (mots(0) == "CTRIA3")
	      {
		int ref = to_num<int>(mots(2));
		int n1 = to_num<int>(mots(3))-1;
		int n2 = to_num<int>(mots(4))-1;
		int n3 = to_num<int>(mots(5))-1;
		all_faces.push_back(TinyVector<int, 5>(n1, n2, n3, -1, ref));
	      }
	    else if (mots(0) == "CQUAD4")
	      {
		int ref = to_num<int>(mots(2));
		int n1 = to_num<int>(mots(3))-1;
		int n2 = to_num<int>(mots(4))-1;
		int n3 = to_num<int>(mots(5))-1;
		int n4 = to_num<int>(mots(6))-1;
		all_faces.push_back(TinyVector<int, 5>(n1, n2, n3, n4, ref));
	      }
	    else if (mots(0) == "CHEXA")
	      {
		elt(0) = to_num<int>(mots(2));
		elt(1) = to_num<int>(mots(3))-1;
		elt(2) = to_num<int>(mots(4))-1;
		elt(3) = to_num<int>(mots(5))-1;
		elt(4) = to_num<int>(mots(6))-1;
		elt(5) = to_num<int>(mots(7))-1;
		elt(6) = to_num<int>(mots(8))-1;
		getline(file_in, ligne);
		StringTokenize(ligne, mots, " \t");
		elt(7) = to_num<int>(mots(1))-1;
		elt(8) = to_num<int>(mots(2))-1;
		all_elt.push_back(elt);
	      }
	    else if (mots(0) == "CTETRA")
	      {
		elt(0) = to_num<int>(mots(2));
		elt(1) = to_num<int>(mots(3))-1;
		elt(2) = to_num<int>(mots(4))-1;
		elt(3) = to_num<int>(mots(5))-1;
		elt(4) = to_num<int>(mots(6))-1;
		elt(5) = -1; elt(6) = -1; elt(7) = -1; elt(8) = -1;
		all_elt.push_back(elt);
	      }
	    else if (mots(0) == "CPENTA")
	      {
		elt(0) = to_num<int>(mots(2));
		elt(1) = to_num<int>(mots(3))-1;
		elt(2) = to_num<int>(mots(4))-1;
		elt(3) = to_num<int>(mots(5))-1;
		elt(4) = to_num<int>(mots(6))-1;
		elt(5) = to_num<int>(mots(7))-1;
		elt(6) = to_num<int>(mots(8))-1;
		elt(7) = -1; elt(8) = -1;
		if (elt(3) == elt(6))
		  {
		    // pyramid
		    int n1 = elt(1), n2 = elt(2), n3 = elt(5), n4 = elt(4), n5 = elt(3);
		    elt(1) = n1; elt(2) = n2; elt(3) = n3; elt(4) = n4; elt(5) = n5;
		    elt(6) = -1;
		  }
		
		all_elt.push_back(elt);
	      }
	    
	  }

	this->ReallocateVertices(all_vertices.size());
	list<R3>::iterator iter_pt; int nb_vertices = 0;
	for (iter_pt = all_vertices.begin(); iter_pt != all_vertices.end(); iter_pt++)
	  this->Vertices(nb_vertices++) = *iter_pt;

	all_vertices.clear();
	
	// faces
	ReallocateBoundariesRef(all_faces.size());
	TinyVector<int, 5> num;
	list<TinyVector<int, 5> >::iterator iter_face; int nface = 0;
	for (iter_face = all_faces.begin(); iter_face != all_faces.end(); iter_face++)
	  {
            num = *iter_face;
            if (num(3) == -1)
              this->BoundaryRef(nface).InitTriangular(num(0), num(1), num(2), num(4));
            else
              this->BoundaryRef(nface).InitQuadrangular(num(0), num(1), num(2), num(3), num(4));
            
            nface++;
          }

	all_faces.clear();

	// elements
	this->ReallocateElements(all_elt.size());
	list<TinyVector<int, 9> >::iterator iter_elt; int nelt = 0;
	for (iter_elt = all_elt.begin(); iter_elt != all_elt.end(); iter_elt++)
	  {
            elt = *iter_elt;
	    if (elt(8) == -1)
	      {
		if (elt(5) == -1)
		  this->Element(nelt).InitTetrahedral(elt(1), elt(2), elt(3), elt(4), elt(0));
		else if (elt(6) == -1)
		  this->Element(nelt).InitPyramidal(elt(1), elt(2), elt(3), elt(4), elt(5), elt(0));
		else
		  this->Element(nelt).InitWedge(elt(1), elt(2), elt(3), elt(4),
						elt(5), elt(6), elt(0));		  
	      }
	    else
	      this->Element(nelt).InitHexahedral(elt(1), elt(2), elt(3), elt(4),
						 elt(5), elt(6), elt(7), elt(8), elt(0));

	    nelt++;
	  }		
      }
    else if (!extension_input.compare("obj"))
      {
	string ligne;
	std::list<R3> all_vertices; R3 point;
        std::list<TinyVector<int, 4> > all_faces;
        TinyVector<int, 4> num;
        Vector<string> param, list_num;
	while (!file_in.eof())
	  {
            getline(file_in, ligne);
	    if (ligne.size() > 2)
              {
                if (ligne.substr(0, 2) == "v ")
                  {
                    istringstream input(ligne.substr(2));
                    input >> point(0) >> point(1) >> point(2);
		    all_vertices.push_back(point);
		  }
                else if (ligne.substr(0, 2) == "f ")
                  {
                    StringTokenize(ligne.substr(2), param, string(" "));
                    num.Fill(-1);
                    for (int i = 0; i < param.GetM(); i++)
                      {
                        size_t pos = param(i).find("/",0);
                        if (pos != string::npos)
                          num(i) = to_num<int>(param(i).substr(0, pos))-1;   
                      }
                    
                    all_faces.push_back(num);
                  }
	      }
	  }
	
	this->ReallocateVertices(all_vertices.size());
	list<R3>::iterator iter_pt; int nb_vertices = 0;
	for (iter_pt = all_vertices.begin(); iter_pt != all_vertices.end(); iter_pt++)
	  this->Vertices(nb_vertices++) = *iter_pt;
	
	// faces
	ReallocateBoundariesRef(all_faces.size());
	list<TinyVector<int, 4> >::iterator iter_face; int nface = 0;
	for (iter_face = all_faces.begin(); iter_face != all_faces.end(); iter_face++)
	  {
            num = *iter_face;
            if (num(3) == -1)
              this->BoundaryRef(nface).InitTriangular(num(0), num(1), num(2), 1);
            else
              this->BoundaryRef(nface).InitQuadrangular(num(0), num(1), num(2), num(3), 1);
            
            nface++;
          }
        
        this->FindConnectivity();
        this->RedistributeReferences();
      }
    else if (!extension_input.compare("hmascii"))
      {
        string node_chaine = "BEGIN NODES";
	bool test_loop = true;
        string ligne;
	while (test_loop)
	  {
	    getline(file_in, ligne, '\n');
	    if (!ligne.compare(0, node_chaine.size(), node_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	test_loop = true;
        
        // we read vertices
        node_chaine = "END NODES";
        list<R3> points; list<int> num_point;
        Vector<string> param; string delim = "(,)";
        Real_wp x, y, z;
        while (test_loop)
          {
            getline(file_in, ligne, '\n');
            
            if (!ligne.compare(0, node_chaine.size(), node_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
            
            if (test_loop)
              {
                StringTokenize(ligne, param, delim);
                
                num_point.push_back(to_num<int>(param(1)));
                
                x = to_num<Real_wp>(param(2));
                y = to_num<Real_wp>(param(3));
                z = to_num<Real_wp>(param(4));
                points.push_back(R3(x, y, z));
              }
	    
          }
        
        IVect NumPoint, InverseNumPoint; VectR3 AllCoord;
        Copy(points, AllCoord);
        Copy(num_point, NumPoint);
        
        this->ReallocateVertices(AllCoord.GetM());
        for (int i = 0; i < AllCoord.GetM(); i++)
          this->Vertex(i) = AllCoord(i);

        int num_max = 0;
        for (int i = 0; i < NumPoint.GetM(); i++)
          num_max = max(num_max, NumPoint(i));
        
        InverseNumPoint.Reallocate(num_max+1);
        for (int i = 0; i < NumPoint.GetM(); i++)
          InverseNumPoint(NumPoint(i)) = i;
        
        num_point.clear(); points.clear();
        
        test_loop = true;
        node_chaine = "BEGIN COMPONENTS";
	while (test_loop)
	  {
	    getline(file_in, ligne, '\n');
            if (!ligne.compare(0, node_chaine.size(), node_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	test_loop = true;
        
        // then we read elements
        node_chaine = "END COMPONENTS";
        list<TinyVector<int, 10> > num_nodes;
        TinyVector<int, 10> num; 
        int nb_elt = 0, nb_faces = 0;
        while (test_loop)
          {
            getline(file_in, ligne, '\n');
                        
            if (!ligne.compare(0, node_chaine.size(), node_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
            
            if (test_loop)
              {
                if (ligne.find("component") != string::npos)
                  {
                    // we don't treat this line
                  }
                else
                  {                    
                    StringTokenize(ligne, param, delim);
                    
                    if (ligne.find("tria3") != string::npos)
                      {
                        num(0) = 2;
                        num(1) = InverseNumPoint(to_num<int>(param(3)));
                        num(2) = InverseNumPoint(to_num<int>(param(4)));
                        num(3) = InverseNumPoint(to_num<int>(param(5)));
                        nb_faces++;
                      }
                    else if (ligne.find("hexa8") != string::npos)
                      {
                        num(0) = 8;                        
                        num(1) = InverseNumPoint(to_num<int>(param(3)));
                        num(2) = InverseNumPoint(to_num<int>(param(4)));
                        num(3) = InverseNumPoint(to_num<int>(param(5)));
                        num(4) = InverseNumPoint(to_num<int>(param(6)));
                        num(5) = InverseNumPoint(to_num<int>(param(7)));
                        num(6) = InverseNumPoint(to_num<int>(param(8)));
                        num(7) = InverseNumPoint(to_num<int>(param(9)));
                        num(8) = InverseNumPoint(to_num<int>(param(10)));
                        nb_elt++;
                      }
                    else if (ligne.find("pyramid5") != string::npos)
                      {
                        num(0) = 5;                        
                        num(1) = InverseNumPoint(to_num<int>(param(3)));
                        num(2) = InverseNumPoint(to_num<int>(param(4)));
                        num(3) = InverseNumPoint(to_num<int>(param(5)));
                        num(4) = InverseNumPoint(to_num<int>(param(6)));
                        num(5) = InverseNumPoint(to_num<int>(param(7)));
                        nb_elt++;
                      }
                    else if (ligne.find("tetra4") != string::npos)
                      {
                        num(0) = 4;                        
                        num(1) = InverseNumPoint(to_num<int>(param(3)));
                        num(2) = InverseNumPoint(to_num<int>(param(4)));
                        num(3) = InverseNumPoint(to_num<int>(param(5)));
                        num(4) = InverseNumPoint(to_num<int>(param(6)));
                        nb_elt++;
                      }
                    
                    num_nodes.push_back(num);
                  }
              }	    
          }        
        
        this->ReallocateElements(nb_elt);
        this->ReallocateBoundariesRef(nb_faces);
        nb_elt = 0; nb_faces = 0;
        list<TinyVector<int, 10> >::iterator iter;
        for (iter = num_nodes.begin(); iter != num_nodes.end(); ++iter)
          {
            num = *iter;
            switch (num(0))
              {
              case 2 :
                {
                  this->BoundaryRef(nb_faces).InitTriangular(num(1), num(2), num(3), 1);
                  nb_faces++;
                }
                break;
              case 4 :
                {
                  this->Element(nb_elt).InitTetrahedral(num(1), num(2), num(3), num(4), 1);
                  nb_elt++;
                }
                break;
              case 5 :
                {
                  this->Element(nb_elt).InitPyramidal(num(1), num(2), num(3), num(4), num(5), 1);
                  nb_elt++;
                }
                break;
              case 8 :
                {
                  this->Element(nb_elt).InitHexahedral(num(1), num(2), num(3), num(4),
                                                       num(5), num(6), num(7), num(8), 1);
                  nb_elt++;
                }
                break;
              }
          }
        
        this->FindConnectivity();
        
        // we create referenced faces
	this->AddBoundaryFaces();
	this->CreateConnexionEdgeToFace();
	this->ConstructCrossReferenceBoundaryRef();
	
	// now we affect references by splitting boundaries in connex surfaces
	this->RedistributeReferences();
	
      }
    else if (!extension_input.compare("dat"))
      {
        int nb_points, nb_elt, n1, n2, n3;
        file_in >> nb_points >> nb_elt >> n1 >> n2 >> n3;
        // we read nodes        
        Vector<R3> points(nb_points);
        for (int i = 0; i < nb_points; i++)
          file_in >> n1 >> points(i)(0) >> points(i)(1) >> points(i)(2);
        
        string type_elt;
        Vector<TinyVector<int, 8> > num_hex(nb_elt);
        IVect num(20);
        Vector<int> IndexNode(nb_points);
        IndexNode.Fill(-1);
        int nb_vertices = 0;
        for (int i = 0; i < nb_elt; i++)
          {
            file_in >> n1 >> n2 >> type_elt;
            for (int j = 0; j < 20; j++)
              {
                file_in >> num(j);
                num(j)--;
              }
            
            for (int j = 0; j < 8; j++)
              {
                num_hex(i)(j) = num(j);
                if (IndexNode(num(j)) == -1)
                  {
                    IndexNode(num(j)) = 0;
                    nb_vertices++;
                  }
              }
          }
        
        this->ReallocateVertices(nb_vertices);
        nb_vertices = 0;
        for (int i = 0; i < nb_points; i++)
          if (IndexNode(i) == 0)
            {
              IndexNode(i) = nb_vertices;
              this->Vertex(nb_vertices) = points(i);
              nb_vertices++;
            }
        
        this->ReallocateElements(nb_elt);
        num.Reallocate(8);
        int ref = 1;
        for (int i = 0; i < nb_elt; i++)
          {
            for (int j = 0; j < 8; j++)
              num(j) = IndexNode(num_hex(i)(j));
            
            this->Element(i).Init(num, ref);
          }
        
        this->FindConnectivity();
        
        // we create referenced faces
	this->AddBoundaryFaces();
	this->CreateConnexionEdgeToFace();
	this->ConstructCrossReferenceBoundaryRef();
	
	// now we affect references by splitting boundaries in connex surfaces
	this->RedistributeReferences();
	
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
	if (ascii)
	  for (int i = 0; i < this->GetNbVertices(); i++)
	    file_in >> this->Vertex(i)(0) >> this->Vertex(i)(1) >> this->Vertex(i)(2);
	else
	  {
	    getline(file_in, ligne);
	    VectReal_wp Pts;
	    Pts.SetData(3*this->GetNbVertices(), 
			reinterpret_cast<Real_wp*>(this->Vertices.GetData()));
	    
	    ReadBinaryDoubleOrFloat(Pts, file_in, double_prec, false, true);
	    Pts.Nullify();
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
	
	// counting number of faces and elements
	int nb_faces_ref = 0, nb_elt = 0;
	for (int i = 0; i < nb_entity; i++)
	  {
	    switch (TypeEntity(i))
	      {
	      case 5 :
	      case 9 :
		nb_faces_ref++;
		break;
	      case 10 :
	      case 12 :
	      case 13 :
	      case 14 :
		nb_elt++;
		break;
	      }
	  }
	
	this->ReallocateBoundariesRef(nb_faces_ref);
	this->ReallocateElements(nb_elt);
	int nb = 0;
	nb_faces_ref = 0; nb_elt = 0;
	for (int i = 0; i < nb_entity; i++)
	  {
	    // int nb_vert = AllNum(nb);
	    switch (TypeEntity(i))
	      {
	      case 5 :
		this->BoundaryRef(nb_faces_ref).
		  InitTriangular(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), RefEntity(i));
		
		nb += 4; nb_faces_ref++;
		break;
	      case 9 :
		this->BoundaryRef(nb_faces_ref).
		  InitQuadrangular(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), AllNum(nb+4), RefEntity(i));
		
		nb += 5; nb_faces_ref++;
		break;
	      case 10 :
		this->Element(nb_elt).
		  InitTetrahedral(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), 
				  AllNum(nb+4), RefEntity(i));
		
		nb += 5; nb_elt++;
		break;
	      case 14 :
		this->Element(nb_elt).
		  InitPyramidal(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), 
				AllNum(nb+4), AllNum(nb+5), RefEntity(i));
		
		nb += 6; nb_elt++;
		break;
	      case 13 :
		this->Element(nb_elt).
		  InitWedge(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), 
			    AllNum(nb+4), AllNum(nb+5), AllNum(nb+6), RefEntity(i));
		
		nb += 7; nb_elt++;
		break;
	      case 12 :
		this->Element(nb_elt).
		  InitHexahedral(AllNum(nb+1), AllNum(nb+2), AllNum(nb+3), 
				 AllNum(nb+4), AllNum(nb+5), AllNum(nb+6),
				 AllNum(nb+7), AllNum(nb+8), RefEntity(i));
		
		nb += 9; nb_elt++;
		break;
	      }
	  }
      }
    else if (!extension_input.compare("exo"))
      {
	// ascii exodus mesh
	
	// cout<<"lecture fichier exo"<<endl;
	string ligne; Vector<string> parameters;
	int nb_elt = 0, nb_faces_ref = 0;
	
	// reading dimensions first
	string dim_chaine = "dimensions:";
	string variables_chaine = "variables:";
	string data_chaine = "data:";
	bool test_loop = true;
	while (test_loop)
	  {
	    getline(file_in, ligne, '\n');
	    if (!ligne.compare(0,dim_chaine.size(), dim_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	test_loop = true;
	IVect nb_nodes_per_block(10); nb_nodes_per_block.Fill(4);
	IVect nb_faces_block(10); nb_faces_block.Fill(0);
	IVect nb_elt_block(10); nb_elt_block.Fill(0);
	int nodl = 0;
	while (test_loop)
	  {
	    getline(file_in, ligne, '\n');
	    if (!ligne.compare(0,variables_chaine.size(), variables_chaine))
	      test_loop = false;
	    else
	      {
		// tokenization of the line
		if ((ligne.size()!=0)&&(ligne[0]!='#')&&(ligne[0]!='/'))
		  {
		    StringTokenize(ligne, parameters, string(" \t"));
		    if (!parameters(0).compare("num_nodes"))
		      nodl = to_num<int>(parameters(2));
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
			nb_nodes_per_block.Reallocate(nb_ref+1); nb_nodes_per_block.Fill(4);
		      }
		    else if (!parameters(0).compare("num_side_sets"))
		      {
			int nb_ref = to_num<int>(parameters(2));
			nb_faces_block.Reallocate(nb_ref+1); nb_faces_block.Fill(0);
		      }
		    else if (!parameters(0).compare(0,11,"num_side_ss"))
		      {
			int ref = to_num<int>(parameters(0).substr(11));
			nb_faces_block(ref) = to_num<int>(parameters(2));
			nb_faces_ref += nb_faces_block(ref);
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
	
	this->ReallocateVertices(nodl);
	this->ReallocateElements(nb_elt); ReallocateBoundariesRef(nb_faces_ref);
	
	test_loop = true;
	// we skip variables: field
	while (test_loop)
	  {
	    getline(file_in, ligne, '\n');
	    if (!ligne.compare(0, data_chaine.size(), data_chaine))
	      test_loop = false;
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	Vector<IVect> Nodle(nb_elt);
	IVect ElemFaceRef(nb_faces_ref); IVect LocalFaceRef(nb_faces_ref);
	IVect ReferenceFace(nb_faces_ref); IVect elem_map(nb_elt); char comma;
	nb_faces_ref = 0; int num_local = 0; nb_elt = 0;
	test_loop = true;
	
	// matching array between numbering of exodus 
        // and numbering of Montjoie for second-order hex
	IVect corres_hex27(27); corres_hex27.Fill();
	corres_hex27(20) = 23; corres_hex27(21) = 25; corres_hex27(22) = 21;
	corres_hex27(23) = 22; corres_hex27(24) = 26; corres_hex27(25) = 24;
	corres_hex27(26) = 20;
	
	// cout<<"Now we read datas "<<endl;
	// now treating datas
        bool file_mesh_curved = false;
	while (test_loop)
	  {
	    getline(file_in, ligne, '\n');
	    if ((ligne.size()!=0)&&(ligne[0]!='#')&&(ligne[0]!='/'))
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
			// here number of elements that belong to faces
			int ref = to_num<int>(keyword.substr(7));
			for (int i = 0; i < parameters.GetM(); i++)
			  {
			    ElemFaceRef(nb_faces_ref) = to_num<int>(parameters(i))-1;
			    ReferenceFace(nb_faces_ref) = ref;
			    nb_faces_ref++;
			  }
			
			for (int i = parameters.GetM(); i < nb_faces_block(ref); i++)
			  {
			    file_in>>ElemFaceRef(nb_faces_ref)>>comma; ElemFaceRef(nb_faces_ref)--;
			    ReferenceFace(nb_faces_ref) = ref;
			    nb_faces_ref++;
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
			    file_in>>LocalFaceRef(num_local)>>comma; LocalFaceRef(num_local)--;
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
			    			    
			    if (nb_vert == 27)
			      {
				file_mesh_curved = true;
				// eight first nodes are vertices
				this->elements(nb_elt).InitHexahedral(num, ref);
				
				// other nodes are "dofs"
				Nodle(nb_elt).Reallocate(27);
				for (int k = 0; k < 27; k++)
				  Nodle(nb_elt)(k) = num(corres_hex27(k));
			      }
			    else if (nb_vert == 8)
			      {
				this->elements(nb_elt).InitHexahedral(num, ref);
			      }
			    else if (nb_vert == 4)
			      {
				this->elements(nb_elt).InitTetrahedral(num, ref);
			      }
			    
			    nb_elt++;
			  }
		      }
		    else if (!keyword.compare("coord"))
		      {
			// coordinates of all the vertices
			for (int i = 0; i < parameters.GetM(); i++)
			  this->Vertices(i)(0) = to_num<T>(parameters(i));
			
			for (int i = parameters.GetM(); i < nodl; i++)
			  file_in>>this->Vertices(i)(0)>>comma;
			
			for (int i = 0; i < nodl; i++)
			  file_in>>this->Vertices(i)(1)>>comma;
			
			for (int i = 0; i < nodl; i++)
			  file_in>>this->Vertices(i)(2)>>comma;
		      }
		  }
		
	      }
	    
	    if (file_in.eof())
	      test_loop = false;
	  }
	
	// now reconstructing faces_ref
	// local numbers of montjoie and exodus are different
	Matrix<int> corres_face(9,6); corres_face.Fill(-1);
	corres_face(8,3) = 0; corres_face(8,0) = 1; corres_face(8,4) = 2;
	corres_face(8,5) = 3; corres_face(8,2) = 4; corres_face(8,1) = 5; 
	
	corres_face(4,3) = 0; corres_face(4,0) = 1; corres_face(4,2) = 2; corres_face(4,1) = 3;    
	for (int i = 0; i < nb_faces_ref; i++)
	  {
	    int num_elem = ElemFaceRef(i);
	    int nb_vertices = this->elements(num_elem).GetNbVertices();
	    int num_loc = corres_face(nb_vertices, LocalFaceRef(i));
	    int ref = ReferenceFace(i);
	    int nb_vertices_face = 3;
            int type_elt = this->Element(num_elem).GetHybridType();
	    if (nb_vertices == 8)
	      nb_vertices_face = 4;
	    
	    IVect num(nb_vertices_face);
	    num(0) = this->elements(num_elem)
	      .numVertex(MeshNumbering<Dimension3>::FirstExtremityFace(type_elt, num_loc));
	    
	    num(1) = this->elements(num_elem)
	      .numVertex(MeshNumbering<Dimension3>::SecondExtremityFace(type_elt, num_loc));
	    
	    num(2) = this->elements(num_elem)
	      .numVertex(MeshNumbering<Dimension3>::ThirdExtremityFace(type_elt, num_loc));
	    
	    if (nb_vertices_face == 4)
	      num(3) = this->elements(num_elem)
		.numVertex(MeshNumbering<Dimension3>::FourthExtremityFace(type_elt, num_loc));
	    
	    faces_ref(i).Init(num, ref);
	  }
        
	if (file_mesh_curved)
	  {
	    // now treating elements Hex27 (curved elements)
	    this->FindConnectivity();
	    int old_order = this->GetGeometryOrder();
	    // the mesh is curved (second order)
	    this->SetGeometryOrder(2);
	    
	    // filling arrays PointsEdge and PointsFace
	    int nb_edges = this->GetNbEdges(); int nb_faces = this->GetNbFaces();
	    VectR3 PointsEdge(nb_edges), PointsFace(nb_faces);
	    VectBool EdgeTreated(nb_edges), FaceTreated(nb_faces);
	    EdgeTreated.Fill(false); FaceTreated.Fill(false);
	    VectBool EdgeCurved(nb_edges), FaceCurved(nb_faces);
	    EdgeCurved.Fill(false); FaceCurved.Fill(false);
	    
	    file_mesh_curved = false;
	    // loop on each element
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (Nodle(i).GetM() == 27)
		  {
		    R3 barycenter;
		    for (int j = 0; j < 12; j++)
		      {
			int num_edge = this->elements(i).numEdge(j);
			if (!EdgeTreated(num_edge))
			  {
			    EdgeTreated(num_edge) = true;
			    PointsEdge(num_edge) = this->Vertices(Nodle(i)(j+8));
			    // is this edge curved ?
			    barycenter.Zero();
			    Add(Real_wp(0.5), this->Vertices(this->edges(num_edge).numVertex(0)),
                                barycenter);
			    Add(Real_wp(0.5), this->Vertices(this->edges(num_edge).numVertex(1)),
                                barycenter);
			    
			    if (barycenter != PointsEdge(num_edge))
			      {
				EdgeCurved(num_edge) = true;
				file_mesh_curved = true;
			      }
			  }
		      }
		    			
		    for (int j = 0; j < 6; j++)
		      {
			int num_face = this->elements(i).numFace(j);
			if (!FaceTreated(num_face))
			  {
			    FaceTreated(num_face) = true;
			    PointsFace(num_face) = this->Vertices(Nodle(i)(j+20));
			    // is this face curved
			    barycenter.Zero();
			    for (int k = 0; k < 4; k++)
			      Add(Real_wp(0.25), this->Vertices(faces(num_face).numVertex(k)), barycenter);
			    
			    // PointsFace(num_face) = barycenter;
                            if (barycenter != PointsFace(num_face))
			      {
				FaceCurved(num_face) = true;
				file_mesh_curved = true;
			      }
			  }
		      }
		  }
	      }
	    
	    if (file_mesh_curved)
	      {
		new_ref = this->GetNewReference();
		this->SetCurveType(new_ref, this->CURVE_FILE);
	      }
	    
	    CreateConnexionEdgeToFace();
	    this->ConstructCrossReferenceBoundaryRef();
	    
	    // adding referenced faces if necessary
	    int nb_new_faces_ref = 0;
	    for (int i = 0; i < nb_faces; i++)
	      if (FaceCurved(i))
                {
                  if (i >= nb_faces_ref)
                    {
		      nb_new_faces_ref++;
		      this->Boundary(i).SetReference(new_ref);
		    }
                  else
                    {
                      int ref = this->BoundaryRef(i).GetReference();
                      this->SetCurveType(ref, this->CURVE_FILE);
                    }
                }
	    
            EdgeTreated.Clear(); FaceTreated.Clear();
            IVect ListeFace(nb_faces_ref + nb_new_faces_ref);
            ListeFace.Fill();
	    if (nb_new_faces_ref > 0)
	      {
		// new faces found
		ResizeBoundariesRef(nb_faces_ref+nb_new_faces_ref);
		for (int i = 0; i < nb_faces; i++)
		  if ((FaceCurved(i))&&(i >= nb_faces_ref))
		    {
		      faces_ref(nb_faces_ref) = faces(i);
		      faces_ref(nb_faces_ref).SetReference(new_ref);
                      ListeFace(nb_faces_ref) = i;
		      nb_faces_ref++;
		    }
	      }
	    
	    // creating references edges
	    int nb_edges_ref = 0;
	    for (int i = 0; i < nb_edges; i++)
	      if (EdgeCurved(i))
		nb_edges_ref++;
	    
	    this->ReallocateEdgesRef(nb_edges_ref);
            IVect ListeEdge(nb_edges_ref);
	    nb_edges_ref = 0;
	    for (int i = 0; i < nb_edges; i++)
	      if (EdgeCurved(i))
		{
		  this->edges_ref(nb_edges_ref) = this->edges(i);
                  int eref = -1;
                  bool edge_curved = false;
                  for (int j = 0; j < this->GetEdge(i).GetNbFaces(); j++)
                    if (this->Boundary(this->GetEdge(i).numFace(j)).GetReference() > 0)
                      {
                        if (!edge_curved)
                          {
                            eref = this->Boundary(this->GetEdge(i).numFace(j)).GetReference();
                            if (this->GetCurveType(eref) > 0)
                              edge_curved = true;
                          }
                      }
                  
                  if (eref == -1)
                    this->edges_ref(nb_edges_ref).SetReference(new_ref);
                  else
                    this->edges_ref(nb_edges_ref).SetReference(eref);
                  
                  ListeEdge(nb_edges_ref) = i;
		  nb_edges_ref++;
		}
	    
	    // now filling the arrays PointsFaceRef and PointsEdgeRef
	    this->PointsFaceRef.Reallocate(nb_faces_ref);
	    this->PointsEdgeRef.Reallocate(nb_edges_ref, 1);
	    for (int i = 0; i < nb_edges_ref; i++)
	      this->PointsEdgeRef(i, 0) = PointsEdge(ListeEdge(i));
	    
	    for (int i = 0; i < nb_faces_ref; i++)
	      {
		this->PointsFaceRef(i).Reallocate(1);
		this->PointsFaceRef(i)(0) = PointsFace(ListeFace(i));
	      }
	    
            this->ReorientElements();
            this->FindConnectivity();
	    
	    // eliminating internal nodes
	    ForceCoherenceMesh();
	    
	    // returning to the previous order of approximation
	    this->SetGeometryOrder(old_order);
	    
	  } // end if file_mesh_curved
      }
    else if (!extension_input.compare("msh"))
      {
	// Gmsh format
	string ligne; int itmp;
	file_in>>ligne;
        
	// two different version for gmsh files
	StringTrim(ligne); int version = 1;
	if (!ligne.compare("$MeshFormat"))
	  {
	    version = 2; double xtmp;
	    file_in>>xtmp>>itmp>>itmp;
	    file_in>>ligne;
	    file_in>>ligne;
	  }
        
        // we read all the nodes
	int nb_all_nodes;
	file_in>>nb_all_nodes;
	IVect node_numbers(nb_all_nodes);
	VectR3 AllNodes(nb_all_nodes);
	
	for (int i = 0; i < nb_all_nodes; i++)
	  file_in>>node_numbers(i)>>AllNodes(i)(0)>>AllNodes(i)(1)>>AllNodes(i)(2);
	
	// maximum node number
	int numero_max = 0;
	for (int i = 0; i < nb_all_nodes; i++)
	  numero_max = max(node_numbers(i), numero_max);
	
	// we read all the elements, faces and edges
	VectBool UsedVertices(numero_max + 1); UsedVertices.Fill(false);
	file_in>>ligne; file_in>>ligne;
	int nb_entites;
	file_in>>nb_entites;
	IVect type_entite(nb_entites), ref_entite(nb_entites), geom_entite(nb_entites);
	Vector<IVect> num_nodes(nb_entites); IVect nb_vertices_entite(nb_entites);
	IVect dimension_entite(nb_entites);
	int nb_elements = 0, nb_faces = 0, nb_edges = 0;
	bool curved_mesh = false; int order_mesh = -1;
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
	    // incrementing elements, edges, faces
	    int dim = NumberMap::DimensionGmshEntity(type_entite(i));
	    dimension_entite(i) = dim;
	    if (dim == 1)
	      nb_edges++;
	    else if ( dim == 2)
	      nb_faces++;
	    else if (dim == 3)
	      nb_elements++;
	    
	    // checking if the order of approximation is the same for all elements
	    int r = NumberMap::OrderGmshEntity(type_entite(i));
	    if (r > 1)
	      {
		curved_mesh = true;
		if (order_mesh == -1)
		  order_mesh = r;
		else
		  if (order_mesh != r)
		    {
		      cout<<"case not treated"<<endl;
		      abort();
		    }
	      }
	    
	    // node numbers are read
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
	
	// we construct edges, faces and elements
	this->ReallocateElements(nb_elements);
	this->ReallocateBoundariesRef(nb_faces);
	this->ReallocateEdgesRef(nb_edges);
	IVect NumEntiteFace(nb_faces), NumEntiteEdge(nb_edges), NumEntiteElement(nb_elements);
	nb_elements = 0; nb_faces = 0; nb_edges = 0;
	IVect num(8);
	for (int i = 0; i < nb_entites; i++)
	  {
	    num.Reallocate(nb_vertices_entite(i));;
	    // permutation
	    for (int j = 0; j < nb_vertices_entite(i); j++)
	      num(j) = perm(num_nodes(i)(j));
	    
	    if (dimension_entite(i) == 1)
	      {
		// edge
		this->EdgeRef(nb_edges).Init(num(0), num(1), ref_entite(i));
		NumEntiteEdge(nb_edges) = i;
		nb_edges++;
	      }
	    else if (dimension_entite(i) == 2)
	      {
		// face
		this->BoundaryRef(nb_faces).Init(num, ref_entite(i));
		this->BoundaryRef(nb_faces).SetGeometryReference(geom_entite(i));
		NumEntiteFace(nb_faces) = i;
		nb_faces++;
	      }
	    else if (dimension_entite(i) == 3)
	      {
		// element
		this->Element(nb_elements).Init(num, ref_entite(i));
		NumEntiteElement(nb_elements) = i;
		nb_elements++;
	      }
	  }
	
	if (this->GetGeometryOrder() > 1)
	  if (curved_mesh)
	    {
	      // construction of reciprocal array for node_numbers
	      for (int i = 0; i < AllNodes.GetM(); i++)
		perm(node_numbers(i)) = i;
	      
	      // available reference
	      new_ref = this->GetNewReference();
	      Vector<bool> RefUsed(this->GetNbReferences()+1); RefUsed.Fill(false);
	      
              // gmsh element numbering
              Vector<IVect> invNum_elt(4);
              Vector<IVect> num_elt(4);
              NumberMap::GetGmshTetrahedralNumbering(order_mesh, invNum_elt(0));
              NumberMap::GetGmshPyramidalNumbering(order_mesh, invNum_elt(1));
              NumberMap::GetGmshPrismaticNumbering(order_mesh, invNum_elt(2));
              NumberMap::GetGmshHexahedralNumbering(order_mesh, invNum_elt(3));
              for (int k = 0; k < 4; k++)
                {
                  num_elt(k).Reallocate(invNum_elt(k).GetM());
                  for (int i = 0; i < num_elt(k).GetM(); i++)
                    num_elt(k)(invNum_elt(k)(i)) = i;
                }
              
	      // construction of regular points on triangle/quadrangle
	      int r = this->GetGeometryOrder();
	      int rm = order_mesh;
	      VectR2 PointsRegTri((rm+1)*(rm+2)/2);
	      VectR2 PointsRegQuad((rm+1)*(rm+1));
	      VectReal_wp PointsReg(rm+1); PointsReg.Fill(); Mlt(Real_wp(1)/rm, PointsReg);
	      Matrix<int> NumNodes2D_quad(rm+1, rm+1), NumNodes2D_tri(rm+1, rm+1);
	      NumNodes2D_quad.Fill(-1); NumNodes2D_tri.Fill(-1);
	      NumNodes2D_quad(0,0) = 0; NumNodes2D_quad(rm, 0) = 1;
	      NumNodes2D_quad(rm,rm) = 2; NumNodes2D_quad(0, rm) = 3;
	      NumNodes2D_tri(0,0) = 0; NumNodes2D_tri(rm, 0) = 1; NumNodes2D_tri(0, rm) = 2;
	      for (int i = 1; i < rm; i++)
		{
		  NumNodes2D_quad(i, 0) = 3+i;
		  NumNodes2D_quad(rm, i) = 2+rm+i;
		  NumNodes2D_quad(rm-i, rm) = 1+2*rm+i;
		  NumNodes2D_quad(0, rm-i) = 3*rm+i;
		  
		  NumNodes2D_tri(i, 0) = 2+i;
		  NumNodes2D_tri(rm-i, i) = 1+rm+i;
		  NumNodes2D_tri(0,rm-i) = 2*rm+i;
		}
	      
	      int node2 = 3*rm;
	      Vector<R2> Points2D_InsideTri((rm-1)*(rm-2)/2);
	      for (int j = 1; j < rm; j++)
		for (int i = 1; i < rm-j; i++)
		  {
		    NumNodes2D_tri(i, j) = node2;
		    Points2D_InsideTri(node2-3*rm).Init(PointsReg(i), PointsReg(j));
		    node2++;
		  }
	      
	      node2 = 4*rm;
	      Matrix<int> NumNodes2D_InsideQuad(rm-1, rm-1);
	      for (int j = 1; j < rm; j++)
		for (int i = 1; i < rm; i++)
		  {
		    NumNodes2D_quad(i, j) = node2;
		    NumNodes2D_InsideQuad(i-1, j-1) = node2-4*rm;
		    node2++;
		  }
	      
	      for (int i = 0; i <= rm; i++)
		{
		  for (int j = 0; j <= rm; j++)
		    PointsRegQuad(NumNodes2D_quad(i, j))
		      = R2(PointsReg(i), PointsReg(j));
		  
		  for (int j = 0; j <= rm-i; j++)
		    PointsRegTri(NumNodes2D_tri(i, j))
		      = R2(PointsReg(i), PointsReg(j));
		}
              
              // computation of correspondence between gmsh NumNodes2D and our NumNodes2D
              IVect SymNodesTri((rm-1)*(rm-2)/2), SymNodesQuad((rm-1)*(rm-1));
              node2 = 0;
              for (int i = 1; i < rm; i++)
                for (int j = 1; j < rm; j++)		
		  {
                    SymNodesQuad(node2) = NumNodes2D_quad(i, j) - 4*rm;                    
                    node2++;
                  }
              
              node2 = 0;
              for (int i = 1; i < rm; i++)
                for (int j = 1; j < rm-i; j++)		
		  {
                    SymNodesTri(node2) = NumNodes2D_tri(i, j) - 3*rm;                    
                    node2++;
                  }
              
              // rotation of Gmsh nodes on faces
	      Matrix<int> FaceRotTri, FaceRotQuad;
	      MeshNumbering<Dimension3, T>::
		GetRotationQuadrilateralFace(NumNodes2D_InsideQuad, FaceRotQuad);
	      
	      MeshNumbering<Dimension3, T>::
		GetRotationTriangularFace(Points2D_InsideTri, FaceRotTri);
	      
	      // construction of reference triangle/quadrangle with those regular points
	      TriangleGeomReference tri_reg;
	      QuadrangleGeomReference quad_reg;
	      tri_reg.ConstructFiniteElement(rm);
	      quad_reg.ConstructFiniteElement(rm);
	      tri_reg.SetNodalPoints(rm, PointsRegTri);
	      quad_reg.SetNodalPoints(rm, PointsRegQuad);
	      
	      // computation of ValPhi for interpolation from regular points
	      Matrix<Real_wp> ValPhiTri((rm+1)*(rm+2)/2, (r-1)*(r-2)/2),
                ValPhiQuad((rm+1)*(rm+1), (r-1)*(r-1));
	      VectReal_wp phi;
	      for (int i = 3*r; i < this->triangle_reference.GetNbPointsNodalElt(); i++)
		{
                  tri_reg.ComputeValuesPhiNodalRef(this->triangle_reference.PointsNodalND(i), phi);
		  for (int j = 0; j < ValPhiTri.GetM(); j++)
		    ValPhiTri(j, i-3*r) = phi(j);
		}
              
	      for (int i = 4*r; i < this->quadrangle_reference.GetNbPointsNodalElt(); i++)
		{
		  quad_reg.ComputeValuesPhiNodalRef(this->quadrangle_reference
                                                    .PointsNodalND(i), phi);
		  for (int j = 0; j < ValPhiQuad.GetM(); j++)
		    ValPhiQuad(j, i-4*r) = phi(j);
		}
	      
	      Globatto<Real_wp> lob_reg; lob_reg.AffectPoints(PointsReg);
	      Matrix<Real_wp> ValPhi(rm+1, r-1);
	      for (int i = 0; i <= rm; i++)
		for (int j = 1; j < r; j++)
		  ValPhi(i, j-1) = lob_reg.EvaluatePhi(i, this->lob_curve.Points(j));
	      
	      // we find connectivity
	      FindConnectivity();
	      
	      // now we are retrieving internal points of edges if the mesh is curved
	      nb_edges = this->GetNbEdgesRef();
	      int nb_edges_ref = nb_edges;
	      int nb_faces_ref = this->GetNbBoundaryRef();
	      this->ResizeEdgesRef(this->GetNbEdges());
	      VectR3 OldPoints(order_mesh+1);
	      for (int i = 0; i < this->GetNbEdges(); i++)
		{
		  if (i >= nb_edges_ref)
		    {
		      // we retrieve the points of the edge
		      OldPoints(0) = this->Vertex(this->GetEdge(i).numVertex(0));
		      OldPoints(rm) = this->Vertex(this->GetEdge(i).numVertex(1));
		      bool linear_edge = true;
		      if (this->GetEdge(i).GetNbElements() > 0)
			{
			  linear_edge = false;
			  int num_elem = this->GetEdge(i).numElement(0);
			  int num_loc = this->Element(num_elem).GetPositionEdge(i);
			  int nv = this->Element(num_elem).GetNbVertices();
                          int type = this->Element(num_elem).GetHybridType();
			  int ne = NumEntiteElement(num_elem);
			  if (this->Element(num_elem).GetOrientationEdge(num_loc))
			    for (int k = 1; k < rm; k++)
			      OldPoints(k)
                                = AllNodes(perm(num_nodes(ne)(num_elt(type)(nv + (order_mesh-1)
                                                                            *num_loc + k-1))));
			  else
			    for (int k = 1; k < order_mesh; k++)
			      OldPoints(k)
                                = AllNodes(perm(num_nodes(ne)(num_elt(type)(nv + (order_mesh-1)
                                                                            *num_loc + order_mesh
                                                                            - k - 1))));
			}
		      else if (this->GetEdge(i).GetNbFaces() > 0)
			{
			  for (int j = 0; j < this->GetEdge(i).GetNbFaces(); j++)
			    {
			      int num_face = this->GetEdge(i).numFace(j);
			      if (num_face < nb_faces)
				{
				  linear_edge = false;
				  int ne = NumEntiteFace(num_face);
				  int num_loc = this->BoundaryRef(num_face).GetPositionBoundary(i);
				  int nv = this->BoundaryRef(num_face).GetNbVertices();
				  if (this->BoundaryRef(num_face).GetOrientationEdge(num_loc))
				    for (int k = 1; k < rm; k++)
				      OldPoints(k)
                                        = AllNodes(perm(num_nodes(ne)(nv + (order_mesh-1)*num_loc
                                                                      + k-1)));
				  else
				    for (int k = 1; k < order_mesh; k++)
				      OldPoints(k)
                                        = AllNodes(perm(num_nodes(ne)(nv + (order_mesh-1)*num_loc
                                                                      + order_mesh - k - 1)));
				}
			    }
			}
		      
		      if (linear_edge)
			for (int k = 1; k < order_mesh; k++)
			  OldPoints(k) = (1.0-PointsReg(k))*OldPoints(0)
                            + PointsReg(k)*OldPoints(rm);		      
		      
		      // the edge is curved ?
		      bool edge_curved = false;
		      for (int k = 1; k < order_mesh; k++)
			{
			  R_N ptM;
			  Add(1.0-PointsReg(k), OldPoints(0), ptM);
			  Add(PointsReg(k), OldPoints(order_mesh), ptM);
			  if (ptM != OldPoints(k))
			    edge_curved = true;
			}
		      
		      if (edge_curved)
			{
			  this->EdgeRef(nb_edges) = this->GetEdge(i);
			  for (int k = 1; k < r; k++)
			    {
			      this->PointsEdgeRef(nb_edges, k-1).Zero();
			      for (int j = 0; j <= order_mesh; j++)
				Add(ValPhi(j, k-1), OldPoints(j),
                                    this->PointsEdgeRef(nb_edges, k-1));
                              
			    } 
                          
                          int eref = -1;
                          for (int j = 0; j < this->GetEdge(i).GetNbFaces(); j++)
                            if (this->Boundary(this->GetEdge(i).numFace(j)).GetReference() > 0)
                              eref = this->Boundary(this->GetEdge(i).numFace(j)).GetReference();
                          
                          if (eref == -1)
                            {
                              this->EdgeRef(nb_edges).SetReference(new_ref);
                              RefUsed(new_ref) = true;
                            }
                          else
                            {
                              this->EdgeRef(nb_edges).SetReference(eref);
                              RefUsed(eref) = true;
                            }
                          
			  nb_edges++;
			}
		    }
		  else
		    {
		      // this is a referenced edge
		      OldPoints(0) = this->Vertex(this->EdgeRef(i).numVertex(0));
		      OldPoints(order_mesh) = this->Vertex(this->EdgeRef(i).numVertex(1));
		      int ne = NumEntiteEdge(i);
                      bool edge_curved = false;
		      for (int k = 1; k < order_mesh; k++)
			{
			  OldPoints(k) = AllNodes(perm(num_nodes(ne)(k+1)));
			  R_N ptM; Add(1.0-PointsReg(k), OldPoints(0), ptM);
			  Add(PointsReg(k), OldPoints(order_mesh), ptM);
			  if (ptM != OldPoints(k))
			    edge_curved = true;
			}
		      
		      if (edge_curved)
			RefUsed(this->EdgeRef(i).GetReference()) = true;
		      
		      for (int k = 1; k < r; k++)
			{
			  this->PointsEdgeRef(i, k-1).Zero();
			  for (int j = 0; j <= order_mesh; j++)
			    Add(ValPhi(j, k-1), OldPoints(j), this->PointsEdgeRef(i, k-1));
                          
			}
		    }
		}
	      
	      this->ResizeEdgesRef(nb_edges);
	      	      
	      nb_faces = this->GetNbBoundaryRef();
	      nb_faces_ref = nb_faces;
	      this->ResizeBoundariesRef(this->GetNbBoundary());
	      for (int i = 0; i < this->GetNbBoundary(); i++)
		{
		  bool tri = (this->Boundary(i).GetNbVertices() == 3);
		  if (tri)
		    OldPoints.Reallocate((rm+1)*(rm+2)/2);
		  else
		    OldPoints.Reallocate((rm+1)*(rm+1));
		  
		  if (i >= nb_faces_ref)
		    {
		      // we retrieve the points of the edge
		      for (int k = 0; k < this->Boundary(i).GetNbVertices(); k++)
			OldPoints(k) = this->Vertex(this->Boundary(i).numVertex(k));
		      
		      int num_elem = this->Boundary(i).numElement(0);
		      int num_loc = this->Element(num_elem).GetPositionBoundary(i);
		      int nv = this->Element(num_elem).GetNbVertices();
		      int ne = NumEntiteElement(num_elem);
		      int ind = this->Boundary(i).GetNbVertices();
		      int type = this->Element(num_elem).GetHybridType();
                      for (int j = 0; j < this->Boundary(i).GetNbEdges(); j++)
			{
			  int num_edge = this->Boundary(i).numEdge(j);
			  int pos = this->Element(num_elem).GetPositionEdge(num_edge);
			  if (this->Element(num_elem).GetOrientationEdge(pos)
                              == this->Boundary(i).GetOrientationEdge(j))
			    for (int k = 1; k < rm; k++)
			      OldPoints(ind++) = AllNodes(perm(num_nodes(ne)
                                                               (num_elt(type)(nv + (rm-1)*pos
                                                                              + k-1))));
			  else
			    for (int k = 1; k < rm; k++)
			      OldPoints(ind++) = AllNodes(perm(num_nodes(ne)
                                                               (num_elt(type)(nv + (rm-1)*pos
                                                                              + rm - k - 1))));
			}
		      
		      int offset = nv + this->Element(num_elem).GetNbEdges()*(rm-1);
		      for (int j = 0; j < num_loc; j++)
			{
			  int nf = this->Element(num_elem).numFace(j);
			  if (this->Boundary(nf).GetNbVertices() == 3)
			    offset += (rm-1)*(rm-2)/2;
			  else
			    offset += (rm-1)*(rm-1);
			}
		      
		      int rot = this->Element(num_elem).GetOrientationFace(num_loc);
                      // opposite orientation :
                      rot = MeshNumbering<Dimension3, T>::GetOppositeOrientationFace(rot, tri);
                                            
		      int nb_pts_bound = ind;
		      if (tri)
			for (int k = 0; k < (rm-1)*(rm-2)/2; k++)
			  OldPoints(ind + FaceRotTri(rot, k))
                            = AllNodes(perm(num_nodes(ne)(num_elt(type)(offset
                                                                        + SymNodesTri(k)))));
		      else
			for (int k = 0; k < (rm-1)*(rm-1); k++)
			  OldPoints(ind+FaceRotQuad(rot, k))
                            = AllNodes(perm(num_nodes(ne)(num_elt(type)(offset
                                                                        + SymNodesQuad(k)))));
		                            
		      // the face is curved ?
		      bool face_curved = false;
		      if (tri)
			for (int k = nb_pts_bound; k < OldPoints.GetM(); k++)
			  {
			    R_N ptM;
			    Add(1.0-PointsRegTri(k)(0)-PointsRegTri(k)(1), OldPoints(0), ptM);
			    Add(PointsRegTri(k)(0), OldPoints(1), ptM);
			    Add(PointsRegTri(k)(1), OldPoints(2), ptM);                            
			    if (ptM != OldPoints(k))
			      {
                                //DISP(k); DISP(ptM); DISP(OldPoints(k));
                                face_curved = true;                       
                                //abort();
                              }
			  }
		      else
			for (int k = nb_pts_bound; k < OldPoints.GetM(); k++)
			  {
			    R_N ptM;
			    Add((1.0-PointsRegQuad(k)(0))*(1.0-PointsRegQuad(k)(1)),
                                OldPoints(0), ptM);
			    Add(PointsRegQuad(k)(0)*(1.0-PointsRegQuad(k)(1)), OldPoints(1), ptM);
			    Add(PointsRegQuad(k)(0)*PointsRegQuad(k)(1), OldPoints(2), ptM);
			    Add((1.0-PointsRegQuad(k)(0))*PointsRegQuad(k)(1), OldPoints(3), ptM);
			    if (ptM != OldPoints(k))
			      face_curved = true;
			  }
			
		      if (face_curved)
			{
			  this->BoundaryRef(nb_faces) = this->Boundary(i);			  
			  if (tri)
                            {
                              this->PointsFaceRef(nb_faces).Reallocate(ValPhiTri.GetN());
                              for (int k = 0; k < ValPhiTri.GetN(); k++)
                                {
                                  this->PointsFaceRef(nb_faces)(k).Zero();
                                  
                                  for (int j = 0; j < ValPhiTri.GetM(); j++)
                                    Add(ValPhiTri(j, k), OldPoints(j),
                                        this->PointsFaceRef(nb_faces)(k));
                                }
                            }
                          else
                            {
                              this->PointsFaceRef(nb_faces).Reallocate(ValPhiQuad.GetN());
                              for (int k = 0; k < ValPhiQuad.GetN(); k++)
				{
                                  this->PointsFaceRef(nb_faces)(k).Zero();
                                  for (int j = 0; j < ValPhiQuad.GetM(); j++)
                                    Add(ValPhiQuad(j, k), OldPoints(j),
                                        this->PointsFaceRef(nb_faces)(k));
                                }
			    } 
			  
                          this->BoundaryRef(nb_faces).SetReference(new_ref);
			  RefUsed(new_ref) = true;
			  nb_faces++;
			}
		    }
		  else
		    {
		      // this is a referenced face
		      for (int k = 0; k < this->BoundaryRef(i).GetNbVertices(); k++)
			OldPoints(k) = this->Vertex(this->BoundaryRef(i).numVertex(k));
		      
		      int ne = NumEntiteFace(i);
                      bool face_curved = false;
		      		      
		      if (tri)
			{
			  for (int k = 3; k < 3*rm; k++)
			    OldPoints(k) = AllNodes(perm(num_nodes(ne)(k)));
			  
                          if (rm > 2)
                            if (num_nodes(ne).GetM() < (rm+1)*(rm+2)/2)
                              {
                                cout << "this triangle has incomplete "
                                     << "order with no points inside" << endl;
                                cout << "Did you write : " <<endl;
                                cout << "Mesh.SecondOrderIncomplete = 0; "<< endl;
                                cout << "in your geo file ?" <<endl;
                                abort();
                              }
                          
			  for (int k = 0; k < (rm-1)*(rm-2)/2; k++)
			    {
			      int k2 = k+3*rm;
			      OldPoints(k2) = AllNodes(perm(num_nodes(ne)(k2)));
			      R_N ptM;
			      Add(1.0-PointsRegTri(k2)(0)-PointsRegTri(k2)(1), OldPoints(0), ptM);
			      Add(PointsRegTri(k2)(0), OldPoints(1), ptM);
			      Add(PointsRegTri(k2)(1), OldPoints(2), ptM);
			      
			      if (ptM != OldPoints(k2))
				face_curved = true;
			    }
			}
		      else
			{
			  for (int k = 4; k < 4*rm; k++)
			    OldPoints(k) = AllNodes(perm(num_nodes(ne)(k)));
			  
                          if (rm > 2)
                            if (num_nodes(ne).GetM() < (rm+1)*(rm+1))
                              {
                                cout << "this quadrangle has incomplete order "
                                     << "with no points inside" << endl;
                                cout << "Did you write : " <<endl;
                                cout << "Mesh.SecondOrderIncomplete = 0; "<< endl;
                                cout << "in your geo file ?" <<endl;
                                abort();
                              }
                          
			  for (int k = 0; k < (rm-1)*(rm-1); k++)
			    {
			      int k2 = k+4*rm;
			      OldPoints(k2) = AllNodes(perm(num_nodes(ne)(k2)));
			      R_N ptM;
			      Add((1.0-PointsRegQuad(k2)(0))*(1.0-PointsRegQuad(k2)(1)),
                                  OldPoints(0), ptM);
			      Add(PointsRegQuad(k2)(0)*(1.0-PointsRegQuad(k2)(1)),
                                  OldPoints(1), ptM);
			      Add(PointsRegQuad(k2)(0)*PointsRegQuad(k2)(1), OldPoints(2), ptM);
			      Add((1.0-PointsRegQuad(k2)(0))*PointsRegQuad(k2)(1),
                                  OldPoints(3), ptM);
			      
			      if (ptM != OldPoints(k2))
				face_curved = true;
			    }
			}
		      
		      if (face_curved)
			RefUsed(this->BoundaryRef(i).GetReference()) = true;
		      
		      if (tri)
			{
			  this->PointsFaceRef(i).Reallocate(ValPhiTri.GetN());
                          for (int k = 0; k < ValPhiTri.GetN(); k++)
			    {
			      this->PointsFaceRef(i)(k).Zero();
			      
			      for (int j = 0; j < ValPhiTri.GetM(); j++)
				Add(ValPhiTri(j, k), OldPoints(j), this->PointsFaceRef(i)(k));
			      
			    }
                        }
		      else
			{
			  this->PointsFaceRef(i).Reallocate(ValPhiQuad.GetN());
			  for (int k = 0; k < ValPhiQuad.GetN(); k++)
			    {
			      this->PointsFaceRef(i)(k).Zero();
			      
			      for (int j = 0; j < ValPhiQuad.GetM(); j++)
				Add(ValPhiQuad(j, k), OldPoints(j), this->PointsFaceRef(i)(k));
			      
			    } 
			}
			  
		    }
		}
              
              this->ResizeBoundariesRef(nb_faces);
	      
	      for (int i = 1; i < RefUsed.GetM(); i++)
		if (RefUsed(i))		  
		  if (this->GetCurveType(i) <= 0)
                    this->SetCurveType(i, this->CURVE_FILE);
	      
	    }
      }
    else if (!extension_input.compare("neub"))
      {
	// mps neutral format
	int nb_vertices;
	file_in >> nb_vertices;
	
	this->ReallocateVertices(nb_vertices);
	R3 point;
	for (int i = 0; i < nb_vertices; i++)
	  {
	    file_in >> point(0) >> point(1) >> point(2);
	    this->Vertex(i) = point;
	  }
	
	int nb_elt, n0, n1, n2, n3, n4, n5, n6, n7, ref;
	file_in >> nb_elt;
	this->ReallocateElements(nb_elt);
	for (int i = 0; i < nb_elt; i++)
	  {
	    file_in >> ref >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7;
	    this->Element(i).InitHexahedral(n0-1, n1-1, n2-1, n3-1, n4-1, n5-1, n6-1, n7-1, ref);
	  }
	
	int nb_faces;
	file_in >> nb_faces;
	this->ReallocateBoundariesRef(nb_faces);
	for (int i = 0; i < nb_faces; i++)
	  {
	    file_in >> ref >> n0 >> n1 >> n2 >> n3;
	    
	    this->BoundaryRef(i).InitQuadrangular(n0-1, n1-1, n2-1, n3-1, ref);
	  }
      }
    else if (!extension_input.compare("neu"))
      {
	int nb_groups = 0, nb_sets = 0, ncoor = 3, nvector = 3, num_vertex = 0, num_elem = 0;
        R3 point;
	IVect NumCells; VectBool IsCellFace;
	int nb_vertices = 0, nb_elt = 0, nb_faces_ref = 0;
	// Gambit neutral file format
	while (!file_in.eof())
	  {
	    string keyword_section; bool test_loop = true; int num_record = 0;
	    string ligne;
	    getline(file_in, keyword_section);
	    if (keyword_section.find("CONTROL INFO") != string::npos)
	      {
		// section with general informations
		test_loop = !file_in.eof(); num_record = 0;
		while (test_loop)
		  {
		    getline(file_in, ligne);
		    if ((ligne.size() > 0)&&(ligne[0] != '/'))
		      {
			// new record
			num_record++;
			if (num_record == 6)
			  {
			    std::istringstream reader_string(ligne);
			    reader_string>>nb_vertices>>nb_elt>>nb_groups>>nb_sets>>ncoor>>nvector;
			  }
			
			// end of section ?
			if (ligne.find("ENDOFSECTION") != string::npos)
			  test_loop = false;
		      }
		    
		    if (file_in.eof())
		      test_loop = false;
		  }
	      }
	    else if (keyword_section.find("NODAL COORDINATES") != string::npos)
	      {
		// vertex coordinates of the mesh
		if (ncoor != 3)
		  {
		    cout<<"Is it a 3-D mesh ??"<<endl;
		    cout<<"Vertices with only "<<ncoor<<" coordinates"<<endl;
		    abort();
		  }
		this->ReallocateVertices(nb_vertices);
		test_loop = !file_in.eof(); num_record = 0; num_vertex = 0; point.Zero();
		while (test_loop)
		  {
		    getline(file_in, ligne);
		    if ((ligne.size() > 0)&&(ligne[0] != '/'))
		      {
			if (ligne.find("ENDOFSECTION") != string::npos)
			  test_loop = false;
			else
			  {
			    std::istringstream reader_string(ligne);
			    reader_string>>num_vertex>>point(0)>>point(1)>>point(2);
			    this->Vertices(num_vertex-1) = point;
			  }
		      }
		    
		    if (file_in.eof())
		      test_loop = false;
		  }
	      }
	    else if (keyword_section.find("ELEMENTS/CELLS") != string::npos)
	      {
		// section containing elements connectivity
		this->ReallocateElements(nb_elt); ReallocateBoundariesRef(nb_elt);
		NumCells.Reallocate(nb_elt); IsCellFace.Reallocate(nb_elt);
		test_loop = !file_in.eof(); num_record = 0; num_elem = 0; point.Zero();
		int ntype = 0, nb_nodes = 0; IVect num;
		int nb_tetras = 0, nb_hexas = 0, nb_pyramids = 0, nb_wedges = 0;
		int nb_triangles_ref = 0, nb_quadrangles_ref = 0;
		nb_elt = 0; nb_faces_ref = 0;
		while (test_loop)
		  {
		    getline(file_in, ligne);
		    if ((ligne.size() > 0)&&(ligne[0] != '/'))
		      {
			if (ligne.find("ENDOFSECTION") != string::npos)
			  test_loop = false;
			else
			  {
			    std::istringstream reader_string(ligne);
			    reader_string>>num_elem>>ntype>>nb_nodes;
			    num.Reallocate(nb_nodes);
			    if (nb_nodes > 7)
			      {
				for (int j = 0; j < 7; j++)
				  reader_string>>num(j);
						
				// last number on the next line
				getline(file_in, ligne);
				num(7) = to_num<int>(ligne);
			      }
			    else
			      for (int j = 0; j < nb_nodes; j++)
				reader_string>>num(j);
			    
			    // nodes begin at 0 and not 1 in Montjoie
			    for (int j = 0; j < nb_nodes; j++)
			      num(j)--;
			    
			    // ntype = 1 -> edge
			    if (ntype == 1)
			      {
				// nothing to do, we don't store edges for 3-D meshes
			      }
			    else if (ntype == 2)
			      {
				// quadrilateral
				faces_ref(nb_faces_ref).Init(num, 1);
				NumCells(num_elem-1) = nb_faces_ref; IsCellFace(num_elem-1) = true;
				nb_quadrangles_ref++; nb_faces_ref++;
			      }
			    else if (ntype == 3)
			      {
				// triangle
				faces_ref(nb_faces_ref).Init(num, 1);
				NumCells(num_elem-1) = nb_faces_ref; IsCellFace(num_elem-1) = true;
				nb_triangles_ref++; nb_faces_ref++;
			      }
			    else if (ntype == 4)
			      {
				// hexahedron
				this->elements(nb_elt)
                                  .InitHexahedral(num(0), num(1), num(3), num(2),
                                                  num(4), num(5), num(7), num(6), 1);
				NumCells(num_elem-1) = nb_elt; IsCellFace(num_elem-1) = false;
				nb_hexas++; nb_elt++;
			      }
			    else if (ntype == 5)
			      {
				// prism
				this->elements(nb_elt).InitWedge(num, 1);
				NumCells(num_elem-1) = nb_elt; IsCellFace(num_elem-1) = false;
				nb_wedges++; nb_elt++;
			      }
			    else if (ntype == 6)
			      {
				// tetrahedron
				this->elements(nb_elt).InitTetrahedral(num, 1);
				NumCells(num_elem-1) = nb_elt; IsCellFace(num_elem-1) = false;
				nb_tetras++; nb_elt++;
			      }
			    else if (ntype == 7)
			      {
				// pyramid
				this->elements(nb_elt)
                                  .InitPyramidal(num(0), num(1), num(3), num(2), num(4), 1);
				NumCells(num_elem-1) = nb_elt; IsCellFace(num_elem-1) = false;
				nb_pyramids++; nb_elt++;
			      }			    
			    
			  }
		      }
		    
		    if (file_in.eof())
		      test_loop = false;
		  }
		
		this->ResizeElements(nb_elt); ResizeBoundariesRef(nb_faces_ref);
	      }
	    else if (keyword_section.find("ELEMENT GROUP") != string::npos)
	      {
		test_loop = !file_in.eof();
		while (test_loop)
		  {
		    getline(file_in, ligne);
		    		    
		    if ((ligne.size() > 0)&&(ligne[0] != '/'))
		      if (ligne.find("ENDOFSECTION") != string::npos)
			test_loop = false;
		    
		    if (file_in.eof())
		      test_loop = false;

		    if (test_loop)
		      {
			string chaine; int ngroup, nb_cells, material, nflags;
			std::istringstream reader_string(ligne);
			reader_string>>chaine>>ngroup>>chaine
                                     >>nb_cells>>chaine>>material>>chaine>>nflags;
			
			getline(file_in, ligne);
			for (int i = 0; i < max(1, nflags); i++)
			  file_in>>material;
			
			for (int i = 0; i < nb_cells; i++)
			  {
			    file_in>>num_elem;
			    num_elem--;
			    this->Element(num_elem).SetReference(ngroup);
			  }
			
			getline(file_in, ligne);
		      }
		  }		
	      }
	    else
	      {
		test_loop = !file_in.eof(); num_record = 0;
		while (test_loop)
		  {
		    getline(file_in, ligne);
		    if ((ligne.size() > 0)&&(ligne[0] != '/'))
		      if (ligne.find("ENDOFSECTION") != string::npos)
			test_loop = false;
		    
		    if (file_in.eof())
		      test_loop = false;
		  }
	      }
	    
	  }
	
	// permutation to have first tetra, then pyramids, wedges, and finally hexas
	Vector<Volume> old_elt = this->elements;
	int nb = 0;
	for (int i = 0; i < nb_elt; i++)
	  if (old_elt(i).GetNbVertices() == 4)
	    this->elements(nb++) = old_elt(i);
	
	for (int i = 0; i < nb_elt; i++)
	  if (old_elt(i).GetNbVertices() == 5)
	    this->elements(nb++) = old_elt(i);
	
	for (int i = 0; i < nb_elt; i++)
	  if (old_elt(i).GetNbVertices() == 6)
	    this->elements(nb++) = old_elt(i);

	for (int i = 0; i < nb_elt; i++)
	  if (old_elt(i).GetNbVertices() == 8)
	    this->elements(nb++) = old_elt(i);
		
	// we find edges faces
	this->FindConnectivity();
	
	// faces on the interface between domains are referenced
	int new_refb = this->GetNewReference();
	for (int i = 0; i < this->GetNbFaces(); i++)
	  if (this->faces(i).GetNbElements() > 1)
	    {
	      int n1 = this->faces(i).numElement(0);
	      int n2 = this->faces(i).numElement(1);
	      if (this->elements(n1).GetReference() != this->elements(n2).GetReference())
		this->faces(i).SetReference(new_refb);
	    }
	
	// we create referenced faces
	this->AddBoundaryFaces();
	this->CreateConnexionEdgeToFace();
	this->ConstructCrossReferenceBoundaryRef();
	
	// now we affect references by splitting boundaries in connex surfaces
	this->RedistributeReferences();
	
      }
    else if (!extension_input.compare("mh3"))
      {
	// geompack mesh file
	int nb_vertices = 0, nb_elt = 0;
	file_in>>nb_vertices; int ref;
	this->ReallocateVertices(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  file_in>>this->Vertices(i)(0)>>this->Vertices(i)(1)>>this->Vertices(i)(2)>>ref;
	
	int nb_extra_vert; double xtmp;
	file_in>>nb_extra_vert;
	for (int i = 0; i < nb_extra_vert; i++)
	  file_in>>ref>>ref>>xtmp>>xtmp;
	
	int nb_nodes; 
	file_in>>nb_nodes>>nb_elt; IVect num(20); num.Fill(-1);
	this->ReallocateElements(nb_elt);
	for (int i = 0; i < nb_elt; i++)
	  {
	    for (int j = 0; j < abs(nb_nodes); j++)
	      {
		file_in>>num(j);
		num(j)--;
	      }
	    
	    if (num(4) == -1)
	      this->elements(i).InitTetrahedral(num, 1);
	    else if (num(5) == -1)
	      this->elements(i).InitPyramidal(num, 1);
	    else if (num(6) == -1)
	      this->elements(i).InitWedge(num, 1);
	    else
	      this->elements(i).InitHexahedral(num, 1);
	  }
	
      }
    else if (!extension_input.compare("mhs"))
      {
	// geompack mesh file
	int nb_vertices = 0;
	file_in>>nb_vertices; int ref;
	this->ReallocateVertices(nb_vertices);
	for (int i = 0; i < nb_vertices; i++)
	  file_in>>this->Vertices(i)(0)>>this->Vertices(i)(1)>>this->Vertices(i)(2)>>ref;
	
	int nb_extra_vert; double xtmp;
	file_in>>nb_extra_vert;
	for (int i = 0; i < nb_extra_vert; i++)
	  file_in>>ref>>ref>>xtmp>>xtmp;
	
	int nb_nodes, nb_bounds; 
	
	file_in>>nb_nodes>>nb_bounds; IVect num(4); num.Fill(-1);
	if (nb_nodes < 0)
	  nb_nodes += 100;
	else
	  nb_nodes -= 100;

	this->ReallocateBoundariesRef(nb_bounds);
	for (int i = 0; i < nb_bounds; i++)
	  {
	    for (int j = 0; j < abs(nb_nodes); j++)
	      {
		file_in>>num(j);
		num(j)--;
	      }
	    
	    if (num(3) == -1)
	      this->BoundaryRef(i).InitTriangular(num(0), num(1), num(2), ref);
	    else
	      this->BoundaryRef(i).InitQuadrangular(num(0), num(1), num(2), num(3), ref);
	  }
	
      }
    else if (!extension_input.compare("meshb"))
      {
	// fichier de maillage en binaire
	int nb_elt = 0, nb_faces_ref = 0, nb_vertices = 0;
	// on lit le "endian" et la version
	int version, endian, dimension;
	file_in.read(reinterpret_cast<char*>(&endian), sizeof(int));
	file_in.read(reinterpret_cast<char*>(&version), sizeof(int));
	
	// boucle sur tous les mots-cles
	int mot_cle, new_position = 1; int ref; IVect num(8);
	int nb_triangles_ref = 0, nb_quadrangles_ref = 0;
	int nb_hexas = 0, nb_tetras = 0, nb_edges_ref = 0;
	while ((!file_in.eof())&&(new_position!=0))
	  {
	    // on lit le mot cle et la position du fichier pour atteindre nouveau mot cle
	    file_in.read(reinterpret_cast<char*>(&mot_cle), sizeof(int));
	    file_in.read(reinterpret_cast<char*>(&new_position), sizeof(int));
	    
	    if (mot_cle == 3)
	      {
		// dimension
		file_in.read(reinterpret_cast<char*>(&dimension), sizeof(int));
		if (dimension != 3)
		  {
		    cout<<"Ceci est un fichier de maillage 2-D et non 3-D"<<endl;
		    abort();
		  }
	      }
	    else if (mot_cle == 4)
	      {
		// liste des sommets
		file_in.read(reinterpret_cast<char*>(&nb_vertices), sizeof(int));
		this->ReallocateVertices(nb_vertices);
		TinyVector<float,3> point;
		for (int i = 0; i < nb_vertices; i++)
		  {
		    // coordonnees
		    file_in.read(reinterpret_cast<char*>(&point(0)), 3*sizeof(float));
		    this->Vertices(i) = point;
		    // reference du point
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		  }
	      }
	    else if (mot_cle == 5)
	      {
		// liste des aretes
		file_in.read(reinterpret_cast<char*>(&nb_edges_ref), sizeof(int));
		this->ReallocateEdgesRef(nb_edges_ref);
		for (int i = 0; i < nb_edges_ref; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 2*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    this->edges_ref(i).Init(num(0)-1, num(1)-1, ref);
		  }
	      }
	    else if (mot_cle == 6)
	      {
		// liste des triangles
		file_in.read(reinterpret_cast<char*>(&nb_triangles_ref), sizeof(int));
		ResizeBoundariesRef(nb_faces_ref + nb_triangles_ref);
		for (int i = 0; i < nb_triangles_ref; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 3*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    faces_ref(nb_faces_ref++).InitTriangular(num(0)-1, num(1)-1, num(2)-1, ref);
		  }
	      }
	    else if (mot_cle == 7)
	      {
		// liste des quadrilateres
		file_in.read(reinterpret_cast<char*>(&nb_quadrangles_ref), sizeof(int));
		ResizeBoundariesRef(nb_faces_ref + nb_quadrangles_ref);
		for (int i = 0; i < nb_quadrangles_ref; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 4*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    faces_ref(nb_faces_ref++)
                      .InitQuadrangular(num(0)-1, num(1)-1, num(2)-1, num(3)-1, ref);
		  }		
	      }
	    else if (mot_cle == 8)
	      {
		// liste des tetraedres
		file_in.read(reinterpret_cast<char*>(&nb_tetras), sizeof(int));
		this->ResizeElements(nb_elt + nb_tetras);
		for (int i = 0; i < nb_tetras; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 4*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    this->elements(nb_elt++)
                      .InitTetrahedral(num(0)-1, num(1)-1, num(2)-1, num(3)-1, ref);
		  }
	      }
	    else if (mot_cle == 10)
	      {
		// liste des hexaedres
		file_in.read(reinterpret_cast<char*>(&nb_hexas), sizeof(int));
		this->ResizeElements(nb_elt + nb_hexas);
		for (int i = 0; i < nb_hexas; i++)
		  {
		    file_in.read(reinterpret_cast<char*>(&num(0)), 8*sizeof(int));
		    file_in.read(reinterpret_cast<char*>(&ref), sizeof(int));
		    for (int j = 0; j < 8; j++)
		      num(j)--;
                    
                    if (num(7) == num(6))
                      {
                        if (num(5) == num(4))
                          this->Element(nb_elt).InitPyramidal(num(0), num(1), num(2),
                                                              num(3), num(4), ref);
                        else
                          this->Element(nb_elt).InitWedge(num(0), num(1), num(2),
                                                          num(4), num(5), num(6), ref);
                      }
                    else
                      this->Element(nb_elt).InitHexahedral(num, ref);
                    
                    nb_elt++;
		  }
	      }
	    
	    // on passe a la nouvelle position
	    file_in.seekg(new_position, ios::beg);
	  }
	
      }
    else if (!extension_input.compare("mesh"))
      {
	string keyword_section;
	string ligne; int ref, n0, n1, n2, n3,n4,n5,n6,n7;
	int nb_elt = 0, nb_vertices = 0, nb_edges_ref = 0, nb_faces_ref = 0;
	int nb_triangles_ref = 0, nb_quadrangles_ref = 0;
	int nb_tetras = 0, nb_hexas = 0;
	while (!file_in.eof())
	  {
	    getline(file_in, keyword_section);
	    if ((keyword_section.size() > 0)&&(keyword_section[0] != '#'))
	      {
		if (keyword_section.find("Vertices") != string::npos)
		  {
		    // all the vertices
		    file_in>>nb_vertices;
		    this->ReallocateVertices(nb_vertices);
		    for (int i = 0; i < nb_vertices; i++)
		      file_in>>this->Vertices(i)(0)>>this->Vertices(i)(1)
                             >>this->Vertices(i)(2)>>ref;
		  }
		else if (keyword_section.find("Edges") != string::npos)
		  {
		    file_in>>nb_edges_ref;
		    this->ReallocateEdgesRef(nb_edges_ref);
		    for (int i = 0; i < nb_edges_ref; i++)
		      {
			file_in>>n0>>n1>>ref;
			this->edges_ref(i).Init(n0-1, n1-1, ref);
		      }
		  }
		else if (keyword_section.find("Triangles") != string::npos)
		  {
		    file_in>>nb_triangles_ref; 
		    int nb_old_faces = nb_faces_ref; nb_faces_ref += nb_triangles_ref;
		    ResizeBoundariesRef(nb_faces_ref);
		    for (int i = nb_old_faces; i < nb_faces_ref; i++)
		      {
			file_in>>n0>>n1>>n2>>ref;
			faces_ref(i).InitTriangular(n0-1, n1-1, n2-1, ref);
		      }
		  }
		else if (keyword_section.find("Quadrilaterals") != string::npos)
		  {
		    file_in>>nb_quadrangles_ref; 
		    int nb_old_faces = nb_faces_ref; nb_faces_ref += nb_quadrangles_ref;
		    ResizeBoundariesRef(nb_faces_ref);
		    for (int i = nb_old_faces; i < nb_faces_ref; i++)
		      {
			file_in>>n0>>n1>>n2>>n3>>ref;
			faces_ref(i).InitQuadrangular(n0-1, n1-1, n2-1, n3-1, ref);
		      }
		  }
		else if (keyword_section.find("Tetrahedra") != string::npos)
		  {
		    file_in>>nb_tetras; 
		    int nb_old_elt = nb_elt; nb_elt += nb_tetras;
		    this->ResizeElements(nb_elt);
		    for (int i = nb_old_elt; i < nb_elt; i++)
		      {
			file_in>>n0>>n1>>n2>>n3>>ref;
			this->elements(i).InitTetrahedral(n0-1, n1-1, n2-1, n3-1, ref);
		      }
		  }
		else if (keyword_section.find("Hexahedra") != string::npos)
		  {
		    file_in>>nb_hexas; 
		    int nb_old_elt = nb_elt; nb_elt += nb_hexas;
		    this->ResizeElements(nb_elt);
		    for (int i = nb_old_elt; i < nb_elt; i++)
		      {
			file_in>>n0>>n1>>n2>>n3>>n4>>n5>>n6>>n7>>ref;
			// for degenerate hexahedra, detection of pyramids and wedges
			if (n7==n6)
			  {
			    if (n5 == n4)
			      this->elements(i).InitPyramidal(n0-1, n1-1, n2-1, n3-1, n4-1, ref);
			    else
			      this->elements(i).InitWedge(n0-1, n1-1, n2-1, n4-1, n5-1, n6-1, ref);
			  }
			else
			  this->elements(i).InitHexahedral(n0-1, n1-1, n2-1, n3-1,
                                                           n4-1, n5-1, n6-1, n7-1, ref);
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

    // checking if a new reference has been used for curved faces not on the boundary
    if (new_ref > 0)
      {
        // changing edges/faces with this reference with -1 reference
        for (int i = 0; i < this->GetNbEdgesRef(); i++)
          if (this->EdgeRef(i).GetReference() == new_ref)
            this->EdgeRef(i).SetReference(-1);

        for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() == new_ref)
            this->BoundaryRef(i).SetReference(-1);
      }

    // maximal reference of edges ?
    int ref_max = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      ref_max = max(this->BoundaryRef(i).GetReference(), ref_max);
    
    // resizing arrays if necessary
    if (ref_max > this->GetNbReferences())
      this->ResizeNbReferences(ref_max);
    
    // we sort referenced edges and faces
    this->SortPeriodicBoundaries();
    this->SortEdgesRef();
    this->SortBoundariesRef();
    
    this->ReorientElements(); 
    this->FindConnectivity();
    
    // now we add to faces_ref, the faces which are on the boundary of the domain
    this->AddBoundaryFaces();
    
    this->ProjectPointsOnCurves();
    
    // on compte le nombre de tetras/ hexas, triangles/quadrangles
    if (this->print_level >= 4)
      {
        cout << " Mesh containing "<<this->GetNbElt()<<" this->elements "<<endl;
        cout<<GetNbTetrahedra()<<" tetrahedra, "<<GetNbHexahedra()<<" hexahedra, "<<
          GetNbPyramids()<<" pyramids, "<<GetNbWedges()<<" wedges "<<endl;
      }    
  }
  
  
  //! Permutation of local numbers of elements
  template<class T>
  void Mesh<Dimension3, T>::ChangeLocalNumberingMesh()
  {
    // for tetrahedra, we use warburton algorithm so that they will be correctly oriented
    IVect num(4), num2(4);
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	if (this->Element(i).GetNbVertices() == 4)
	  {
	    for (int j = 0; j < 4; j++)
	      num(j) = this->Element(i).numVertex(j);
	    
	    Sort(num);
	    for (int j = 0; j < 4; j++)
	      num2(3-j) = num(j);
	    
	    this->Element(i).Init(num2, this->Element(i).GetReference());
	  }
	else
	  {
	    cout << "Algorithm only implemented for purely tetrahedral meshes " << endl;
	    abort();
	  }
      }
  }
  
  
  //! writing the mesh, several formats are available
  /*!
    \param[in] file_name name of the file where the mesh is stored
  */
  template<class T>
  void Mesh<Dimension3, T>::Write(const string & file_name, bool double_prec, bool ascii) const
  {
    if (this->print_level >= 2)
      cout << " the mesh will be written on file "<<file_name<<endl;
    
    string extension_output = GetExtension(file_name);
    
    ofstream file_out(file_name.data());
    if (!file_out.is_open())
      {
	cout << "Unable to open file " << file_name << endl;
	abort();
      }
    
    Write(file_out, extension_output, double_prec, ascii);
    
    file_out.close();
  }
  

  //! writing the mesh, several formats are available
  /*!
    \param[in] file_out stream where the mesh is written
    \param[in] extension_output format of the mesh (usually extension of the file)
    \param[in] double_prec mesh written in single or double precision ?
    \param[in] ascii mesh written in ascii or binary format ?
  */  
  template<class T>
  void Mesh<Dimension3, T>::Write(ostream& file_out, const string& extension_output,
				  bool double_prec, bool ascii) const
  {  
    if (double_prec)
      file_out.precision(15);
    else
      file_out.precision(7);
    
    file_out.setf(ios_base::scientific,ios_base::floatfield);
    
    int nb_triangles_ref = 0;
    int nb_quadrangles_ref = 0;
    int nb_edges_ref = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      if (this->BoundaryRef(i).GetReference() > 0)
        {
          if (this->BoundaryRef(i).IsTriangular())
            nb_triangles_ref++;
          else
            nb_quadrangles_ref++;
        }
    
    for (int i = 0; i < this->GetNbEdgesRef(); i++)
      if (this->EdgeRef(i).GetReference() > 0)
        nb_edges_ref++;
    
    int nb_faces_ref = nb_triangles_ref + nb_quadrangles_ref;
    int nb_tetras = GetNbTetrahedra();
    int nb_pyramids = GetNbPyramids();
    int nb_wedges = GetNbWedges();
    int nb_hexas = GetNbHexahedra();

    if (!extension_output.compare("mit"))
      {
	// nb de points, nb de faces et 3*nb_triangles+4*nb_quadrangles
	file_out<<this->GetNbVertices()<<'\n'<<nb_faces_ref<<'\n'
                <<(3*nb_triangles_ref+4*nb_quadrangles_ref)<<'\n';
	// ecriture des sommets
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertices(i)(0)<<" "
                  <<this->Vertices(i)(1)<<" "<<this->Vertices(i)(2)<<'\n';
	
	// ecriture des faces
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() > 0)
            {
              for (int j = 0; j < faces_ref(i).GetNbVertices(); j++)
                file_out<<(faces_ref(i).numVertex(j)+1)<<" ";
              
              file_out<<'\n';
            }
	
	// offset pour les faces
	int offset = 0;
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() > 0)
            {
              offset += faces_ref(i).GetNbVertices();
              file_out<<offset<<'\n';
            }
	
	// reference des faces
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  if (this->BoundaryRef(i).GetReference() > 0)
            file_out<<faces_ref(i).GetReference()<<'\n';
	
      }
    else if (!extension_output.compare("gts"))
      {
	file_out<<this->GetNbVertices()<<' '<<this->GetNbEdges()<<' '<<nb_triangles_ref<<'\n';
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertices(i)(0)<<" "
                  <<this->Vertices(i)(1)<<" "<<this->Vertices(i)(2)<<'\n';
	  
	for (int i = 0; i < this->GetNbEdges(); i++)
	  file_out<<this->edges(i).numVertex(0)+1<<" "<<this->edges(i).numVertex(1)+1<<'\n';
	
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  if (this->BoundaryRef(i).GetReference() > 0)
            if (faces_ref(i).GetNbVertices() == 3)
              file_out<<faces_ref(i).numEdge(0)+1<<" "<<faces_ref(i).numEdge(1)+1
                      <<" "<<faces_ref(i).numEdge(2)+1<<'\n';
	
      }
    else if (!extension_output.compare("vf"))
      {
	file_out<<this->GetNbVertices()<<'\n';
	
	// ecriture des sommets
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertices(i)(0)<<" "<<this->Vertices(i)(1)<<" "<<this->Vertices(i)(2)<<'\n';
	
	// ecriture des faces
	file_out<<GetNbFaces()<<'\n';
	for (int i = 0; i < GetNbFaces(); i++)
	  {
	    int nb_vertices_loc = faces(i).GetNbVertices();
	    file_out<<nb_vertices_loc<<" ";
	    for (int j = 0; j < nb_vertices_loc; j++)
	      file_out<<(faces(i).numVertex(j)+1)<<" ";
	    
	    file_out<<faces(i).GetReference()<<'\n';
	  }
	
	// ecriture des elements
	file_out<<this->GetNbElt()<<'\n';
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    int nb_faces_loc = this->Element(i).GetNbFaces();
	    file_out<<nb_faces_loc<<" ";
	    for (int j = 0; j < nb_faces_loc; j++)
	      file_out<<(this->Element(i).numFace(j)+1)<<" ";
	    
	    file_out<<this->Element(i).GetReference()<<" 1\n";
	  }
	
      }
    else if (!extension_output.compare("stl"))
      {
	file_out<<"solid test"<<endl;
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() > 0)
            {
              file_out<<"facet normal 0 0 0\n";
              file_out<<"  outer loop\n";
              for (int j = 0; j < faces_ref(i).GetNbVertices(); j++)
                {
                  int nv = faces_ref(i).numVertex(j);
                  file_out<<"    vertex "<<this->Vertex(nv)(0)<<' '
                          <<this->Vertex(nv)(1)<<' '<<this->Vertex(nv)(2)<<'\n';
                }
              
              file_out<<"  endloop\n";
              file_out<<"endfacet\n";
            }
	file_out<<"endsolid test"<<endl;
      }
    else if (!extension_output.compare("vtk"))
      {
	// vtk file
	file_out<<"# vtk DataFile Version 2.0"<<endl;
	file_out<<"Montjoie 3-D mesh"<<endl;
	if (ascii)
	  file_out<<"ASCII"<<endl;
	else
	  file_out<<"BINARY"<<endl;
	
	file_out<<"DATASET UNSTRUCTURED_GRID"<<endl;
	if (double_prec)
	  file_out<<"POINTS "<<this->GetNbVertices()<<" double "<<endl;
	else
	  file_out<<"POINTS "<<this->GetNbVertices()<<" float "<<endl;
	
	if (ascii)
	  for (int i = 0; i < this->GetNbVertices(); i++)
	    file_out<<this->Vertex(i)(0)<<' '<<this->Vertex(i)(1)<<' '<<this->Vertex(i)(2)<<'\n';
	else
	  {
	    VectReal_wp Pts;
	    Pts.SetData(3*this->GetNbVertices(),
			reinterpret_cast<Real_wp*>(this->Vertices.GetData()));
	    WriteBinaryDoubleOrFloat(Pts, file_out, double_prec, false, true);
	    Pts.Nullify();
	    file_out << '\n';
	  }
	
	int size_all = nb_hexas*9 + nb_wedges*7 + nb_pyramids*6 
          + nb_tetras*5 + nb_quadrangles_ref*5 + nb_triangles_ref*4;
	file_out<<"CELLS "<<this->GetNbElt()+nb_faces_ref<<" "<<size_all<<endl;
	if (ascii)
	  {
	    for (int i = 0; i < GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
		{
		  file_out<<faces_ref(i).GetNbVertices()<<' ';
		  for (int j = 0; j < faces_ref(i).GetNbVertices(); j++)
		    file_out<<(faces_ref(i).numVertex(j))<<' ';
		  
		  file_out<<'\n';
		}
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		file_out<<this->elements(i).GetNbVertices()<<' ';
		for (int j = 0; j < this->elements(i).GetNbVertices(); j++)
		  file_out<<(this->elements(i).numVertex(j))<<' ';
		
		file_out<<'\n';
	      }
	  }
	else
	  {
	    IVect num(size_all);
	    int nb = 0;
	    for (int i = 0; i < GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
		{
		  num(nb++) = faces_ref(i).GetNbVertices();
		  for (int j = 0; j < faces_ref(i).GetNbVertices(); j++)
		    num(nb++) = faces_ref(i).numVertex(j);
		}
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		num(nb++) = this->elements(i).GetNbVertices();
		for (int j = 0; j < this->elements(i).GetNbVertices(); j++)
		  num(nb++) = this->elements(i).numVertex(j);
	      }
	    
	    WriteBinaryInteger(num, file_out, false, true);
	    file_out << '\n';
	  }
        
	file_out<<"CELL_TYPES "<<this->GetNbElt()+nb_faces_ref<<endl;
	if (ascii)
	  {
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
		{
		  if (faces_ref(i).GetNbVertices() == 3)
		    file_out<<"5\n";
		  else
		    file_out<<"9\n";
		}
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (this->elements(i).GetNbVertices() == 4)
		  file_out<<"10\n";
		else if (this->elements(i).GetNbVertices() == 5)
		  file_out<<"14\n";
		else if (this->elements(i).GetNbVertices() == 6)
		  file_out<<"13\n";
		else
		  file_out<<"12\n";
	      }

	    // writing references in a section CELL_DATA
	    file_out << "CELL_DATA " << this->GetNbElt()+nb_faces_ref << '\n';
	    file_out << "SCALARS cell_scalars int 1" << '\n';
	    file_out << "LOOKUP_TABLE default\n";
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
		file_out << this->BoundaryRef(i).GetReference() << '\n';
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      file_out << this->Element(i).GetReference() << '\n';	    
	  }
	else
	  {
	    IVect num(this->GetNbElt()+nb_faces_ref);
	    int nb = 0;
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
		{		  
		  if (faces_ref(i).GetNbVertices() == 3)
		    num(nb++) = 5;
		  else
		    num(nb++) = 9;
		}
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		switch (this->elements(i).GetNbVertices())
		  {
		  case 4 : num(nb++) = 10; break;
		  case 5 : num(nb++) = 14; break;
		  case 6 : num(nb++) = 13; break;
		  case 8 : num(nb++) = 12; break;
		  }
	      }
	    
	    WriteBinaryInteger(num, file_out, false, true);
	    file_out << '\n';

	    // writing references in a section CELL_DATA
	    file_out << "CELL_DATA " << this->GetNbElt()+nb_faces_ref << '\n';
	    file_out << "SCALARS cell_scalars int 1" << '\n';
	    file_out << "LOOKUP_TABLE default\n";
	    nb = 0;
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
		num(nb++) = this->BoundaryRef(i).GetReference();
	    
	    for (int i = 0; i < this->GetNbElt(); i++)
	      num(nb++) = this->Element(i).GetReference();   
	    
	    WriteBinaryInteger(num, file_out, false, true);
	    file_out << '\n';
	  }	
      }
    else if (!extension_output.compare("mh3"))
      {
	// geompack mesh file
	file_out<<this->GetNbVertices()<<'\n'; int ref = 0;
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertex(i)(0)<<" "<<this->Vertex(i)(1)
                  <<" "<<this->Vertex(i)(2)<<" "<<ref<<'\n';
	
	file_out<<'\n'<<"0"<<"\n\n";
	
	file_out<<"-8 "<<this->GetNbElt()<<'\n';
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	      file_out<<this->Element(i).numVertex(j)+1<<" ";
	    
	    for (int j = this->Element(i).GetNbVertices(); j < 8; j++)
	      file_out<<"0 ";
	    
	    file_out<<'\n';
	  }
	
      }
    else if (!extension_output.compare("mhs"))
      {
	// geompack mesh file
	file_out<<this->GetNbVertices()<<'\n'; int ref = 1;
	for (int i = 0; i < this->GetNbVertices(); i++)
	  file_out<<this->Vertex(i)(0)<<" "<<this->Vertex(i)(1)
                  <<" "<<this->Vertex(i)(2)<<" "<<ref<<'\n';
	
	file_out<<'\n'<<"0\n\n";
	
	int nb_nodes = 103, node_max = 3;
	if (nb_quadrangles_ref > 0)
	  {
	    nb_nodes = -104;
	    node_max = 4;
	  }
	
	if (nb_triangles_ref == 0)
	  nb_nodes = 104;
	
	file_out<<nb_nodes<<" "<<nb_faces_ref<<'\n';
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() > 0)
            {
              for (int j = 0; j < this->BoundaryRef(i).GetNbVertices(); j++)
                file_out<<this->BoundaryRef(i).numVertex(j)+1<<" ";
              
              for (int j = this->BoundaryRef(i).GetNbVertices(); j < node_max; j++)
                file_out<<" 0 ";
              
              file_out<<'\n';
            }
	
	file_out<<'\n';
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() > 0)
            {
              file_out<<"1 ";
              for (int j = 0; j < node_max; j++)
                file_out<<" 1";
              
              file_out<<'\n';
            }
	
      }
    else if (!extension_output.compare("msh"))
      {
	Vector<R_N> PosNodes; Vector<IVect> Nodle, NodleSurf; 
	int order = GetNodesCurvedMesh(PosNodes, Nodle, NodleSurf, this->GetGeometryOrder());
	
	// numbering different in Gmsh and Montjoie
	IVect num_tri, num_quad, num_hex, num_tet, num_pyr, num_prism;
	NumberMap::GetGmshTriangularNumbering(order, num_tri);
	NumberMap::GetGmshQuadrilateralNumbering(order, num_quad);
	NumberMap::GetGmshTetrahedralNumbering(order, num_tet);
	NumberMap::GetGmshPyramidalNumbering(order, num_pyr);
	NumberMap::GetGmshPrismaticNumbering(order, num_prism);
	NumberMap::GetGmshHexahedralNumbering(order, num_hex);
	
	// sommets
	file_out<<"$MeshFormat\n"<<"2 0 8\n"<<"$EndMeshFormat\n"<<"$Nodes\n";
	file_out<<PosNodes.GetM()<<'\n';
	for (int i = 0; i < PosNodes.GetM(); i++)
	  file_out<<(i+1)<<' '<<PosNodes(i)(0)<<' '<<PosNodes(i)(1)<<' '<<PosNodes(i)(2)<<'\n';
	
	int nb = nb_edges_ref + nb_faces_ref + this->GetNbElt();
	file_out<<"$EndNodes"<<'\n'<<"$Elements"<<'\n'<<nb<<'\n';
	
	nb = 1;
	// referenced edges
	for (int i = 0; i < this->GetNbEdgesRef(); i++)
          if (this->EdgeRef(i).GetReference() > 0)
            {
              int type_entite = NumberMap::GetGmshEntityNumber(this->EdgeRef(i), order);
              file_out<<nb<<" "<<type_entite<<" 3 "<<this->EdgeRef(i).GetReference()<<" "
                      << " 0 " << " 0 ";
              
              file_out<<(this->EdgeRef(i).numVertex(0)+1)<<' '<< (this->EdgeRef(i).numVertex(1)+1);
              for (int k = 0; k < order-1; k++)
                file_out << ' ' << this->GetNbVertices() + (order-1)*i + k + 1;
              
              file_out<<'\n';
              nb++;
            }
	
	// referenced faces
	for (int i = 0; i < GetNbBoundaryRef(); i++)
          if (this->BoundaryRef(i).GetReference() > 0)
            {
              file_out<<nb<<' '
                      <<NumberMap::GetGmshEntityNumber(this->BoundaryRef(i), order)<<" 3 "<<
                this->BoundaryRef(i).GetReference()<<" " 
                      << this->BoundaryRef(i).GetGeometryReference() <<" 0 ";
              
              int nb_vert = this->BoundaryRef(i).GetNbVertices();
              if (nb_vert == 3)
                {
                  for (int j = 0; j < nb_vert; j++)
                    file_out<<' '<<NodleSurf(i)(num_tri(j)) + 1;
                  
                  for (int j = 0; j < nb_vert; j++)
                    for (int k = 0; k < order-1; k++)
                      file_out << ' ' << NodleSurf(i)(num_tri(nb_vert + (order-1)*j + k)) + 1;
                  
                  for (int j = 3*order; j < NodleSurf(i).GetM(); j++)
                    file_out << ' ' << NodleSurf(i)(num_tri(j)) + 1;
                }
              else
                {
                  for (int j = 0; j < nb_vert; j++)
                    file_out<<' '<<NodleSurf(i)(num_quad(j)) + 1;
                  
                  for (int j = 0; j < nb_vert; j++)
                    for (int k = 0; k < order-1; k++)
                      file_out << ' ' << NodleSurf(i)(num_quad(nb_vert + (order-1)*j + k)) + 1;
                  
                  for (int j = 4*order; j < NodleSurf(i).GetM(); j++)
                    file_out << ' ' << NodleSurf(i)(num_quad(j)) + 1;
                }
              
              file_out<<'\n'; nb++;
            }
	
	// elements
	int nb_inside_tri = (order-1)*(order-2)/2;
	int nb_inside_quad = (order-1)*(order-1);
	int nb_inside_tet = (order-1)*(order-2)*(order-3)/6;
	int nb_inside_pyr = (order-1)*(order-2)*(2*order-3)/6;
	int nb_inside_wedge = nb_inside_tri*(order-1);
	int nb_inside_hex = (order-1)*(order-1)*(order-1);
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    file_out<<nb<<' '<<NumberMap::GetGmshEntityNumber(this->Element(i), order)<<" 3 "<<
	      this->Element(i).GetReference()<<" 0 0 ";
	    
	    int nb_vert = this->Element(i).GetNbVertices();
	    int nb_edges = this->Element(i).GetNbEdges();
	    int nb_faces = this->Element(i).GetNbFaces();
	    if (nb_vert == 4)
	      {
		// tetrahedra
		for (int j = 0; j < nb_vert; j++)
		  file_out<<' '<<Nodle(i)(num_tet(j)) + 1;
		
		int ind = nb_vert;
		for (int j = 0; j < nb_edges; j++)
		  for (int k = 0; k < order-1; k++)
		    file_out << ' ' << Nodle(i)(num_tet(ind++)) + 1;
		
		for (int j = 0; j < nb_faces; j++)
		  for (int k = 0; k < nb_inside_tri; k++)
		    file_out << ' ' << Nodle(i)(num_tet(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_tet; k++)
		  file_out << ' ' << Nodle(i)(num_tet(ind++)) + 1;
	      }
	    else if (nb_vert == 5)
	      {
		// pyramid
		for (int j = 0; j < nb_vert; j++)
		  file_out<<' '<<Nodle(i)(num_pyr(j)) + 1;
		
		int ind = nb_vert;
		for (int j = 0; j < nb_edges; j++)
		  for (int k = 0; k < order-1; k++)
		    file_out << ' ' << Nodle(i)(num_pyr(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_quad; k++)
		  file_out << ' ' << Nodle(i)(num_pyr(ind++)) + 1;
		
		for (int j = 0; j < 4; j++)
		  for (int k = 0; k < nb_inside_tri; k++)
		    file_out << ' ' << Nodle(i)(num_pyr(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_pyr; k++)
		  file_out << ' ' << Nodle(i)(num_pyr(ind++)) + 1;
	      }
	    else if (nb_vert == 6)
	      {
		// prism
		for (int j = 0; j < nb_vert; j++)
		  file_out<<' '<<Nodle(i)(num_prism(j)) + 1;
		
		int ind = nb_vert;
		for (int j = 0; j < nb_edges; j++)
		  for (int k = 0; k < order-1; k++)
		    file_out << ' ' << Nodle(i)(num_prism(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_tri; k++)
		  file_out << ' ' << Nodle(i)(num_prism(ind++)) + 1;
 
		for (int j = 0; j < 3; j++)
		  for (int k = 0; k < nb_inside_quad; k++)
		    file_out << ' ' << Nodle(i)(num_prism(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_tri; k++)
		  file_out << ' ' << Nodle(i)(num_prism(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_wedge; k++)
		  file_out << ' ' << Nodle(i)(num_prism(ind++)) + 1;
	      }
	    else
	      {
		// hexahedron
		for (int j = 0; j < nb_vert; j++)
		  file_out<<' '<<Nodle(i)(num_hex(j)) + 1;
		
		int ind = nb_vert;
		for (int j = 0; j < nb_edges; j++)
		  for (int k = 0; k < order-1; k++)
		    file_out << ' ' << Nodle(i)(num_hex(ind++)) + 1;
		
		for (int j = 0; j < nb_faces; j++)
		  for (int k = 0; k < nb_inside_quad; k++)
		    file_out << ' ' << Nodle(i)(num_hex(ind++)) + 1;
		
		for (int k = 0; k < nb_inside_hex; k++)
		  file_out << ' ' << Nodle(i)(num_hex(ind++)) + 1;
	      }
	    
	    file_out<<'\n'; nb++;
	  }
	
	file_out<<"$EndElements"<<endl;
      }
    else if (!extension_output.compare("meshb"))
      {
	int itmp = 1;
	// ecriture endian et version
	file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	
	// ecriture dimension
	itmp = 3; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	itmp = 5*sizeof(int); file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	itmp = 3; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	
	// ecriture sommets
	itmp = 4; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	int offset = 8*sizeof(int) + this->GetNbVertices()*(3*sizeof(float) + sizeof(int));
	file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	itmp = this->GetNbVertices();
	file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	float xtmp;
	for (int i = 0; i < this->GetNbVertices(); i++)
	  {
	    for (int j = 0; j < 3; j++)
	      {
		xtmp = toDouble(this->Vertices(i)(j));
		file_out.write(reinterpret_cast<char*>(&xtmp), sizeof(float));
	      }
	    itmp = 0; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	  }
	
	IVect num(8);
	// ecriture aretes
	if (this->GetNbEdgesRef() > 0)
	  {
	    itmp = 5; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + this->GetNbEdgesRef()*sizeof(int)*3;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    itmp = this->GetNbEdgesRef(); 
	    file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    for (int i = 0; i < this->GetNbEdgesRef(); i++)
	      {
		num(0) = this->edges_ref(i).numVertex(0)+1;
                num(1) = this->edges_ref(i).numVertex(1)+1;
		file_out.write(reinterpret_cast<char*>(&num(0)), 2*sizeof(int));
		itmp = this->edges_ref(i).GetReference();
                file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	      }
	  }
	
	// ecriture triangles
	if (nb_triangles_ref > 0)
	  {
	    itmp = 6; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + nb_triangles_ref*sizeof(int)*4;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    file_out.write(reinterpret_cast<char*>(&nb_triangles_ref), sizeof(int));
	    for (int i = 0; i < GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
                if (faces_ref(i).GetNbVertices() == 3)
                  {
                    num(0) = faces_ref(i).numVertex(0)+1; num(1) = faces_ref(i).numVertex(1)+1;
                    num(2) = faces_ref(i).numVertex(2)+1;
                    file_out.write(reinterpret_cast<char*>(&num(0)), 3*sizeof(int));
                    itmp = faces_ref(i).GetReference();
                    file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
                  }
	  }
	
	// ecriture quadrangles
	if (nb_quadrangles_ref > 0)
	  {
	    itmp = 7; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + nb_quadrangles_ref*sizeof(int)*5;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    file_out.write(reinterpret_cast<char*>(&nb_quadrangles_ref), sizeof(int));
	    for (int i = 0; i < GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
                if (faces_ref(i).GetNbVertices() == 4)
                  {
                    num(0) = faces_ref(i).numVertex(0)+1; num(1) = faces_ref(i).numVertex(1)+1;
                    num(2) = faces_ref(i).numVertex(2)+1; num(3) = faces_ref(i).numVertex(3)+1;
                    file_out.write(reinterpret_cast<char*>(&num(0)), 4*sizeof(int));
                    itmp = faces_ref(i).GetReference();
                    file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
                  }
	  }
	
	// ecriture tetraedres
	if (nb_tetras > 0)
	  {
	    itmp = 8; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	    offset += 3*sizeof(int) + nb_tetras*sizeof(int)*5;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    file_out.write(reinterpret_cast<char*>(&nb_tetras), sizeof(int));
	    for (int i = 0; i < this->GetNbElt(); i++)
	      if (this->elements(i).GetNbVertices() == 4)
		{
		  for (int j = 0; j < 4; j++)
		    num(j) = this->elements(i).numVertex(j)+1;
		  
		  file_out.write(reinterpret_cast<char*>(&num(0)), 4*sizeof(int));
		  itmp = this->elements(i).GetReference();
                  file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
		}
	  }
	
	// ecriture hexaedres
	if ((nb_hexas > 0)||(nb_pyramids > 0)||(nb_wedges > 0))
	  {
	    itmp = 10; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
            int nb_elt_ = nb_hexas + nb_pyramids + nb_wedges;
	    offset += 3*sizeof(int) + nb_elt_*sizeof(int)*9;
	    file_out.write(reinterpret_cast<char*>(&offset), sizeof(int));
	    file_out.write(reinterpret_cast<char*>(&nb_elt_), sizeof(int));
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
                if (this->Element(i).GetNbVertices() > 4)
                  {
                    switch (this->Element(i).GetNbVertices())
                      {
                      case 8 :
                        for (int j = 0; j < 8; j++)
                          num(j) = this->elements(i).numVertex(j)+1;
                        break;
                      case 6 :
                        for (int j = 0; j < 3; j++)
                          num(j) = this->elements(i).numVertex(j)+1;
                        
                        for (int j = 0; j < 3; j++)
                          num(4+j) = this->elements(i).numVertex(3+j)+1;
                        
                        num(3) = num(2); num(7) = num(6);
                        break;
                      case 5 :
                        for (int j = 0; j < 5; j++)
                          num(j) = this->elements(i).numVertex(j)+1;
                        
                        num(5) = num(4); num(6) = num(4); num(7) = num(4);
                        break;
                      }
                    
                    file_out.write(reinterpret_cast<char*>(&num(0)), 8*sizeof(int));
                    itmp = this->Element(i).GetReference();
                    file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
                  }
              }
	  }
	
	// fin maillage
	itmp = 54; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
	itmp = 0; file_out.write(reinterpret_cast<char*>(&itmp), sizeof(int));
      }
    else if (!extension_output.compare("mesh"))
      {
	file_out<<"MeshVersionFormatted 1"<<'\n'<<"Dimension"<<'\n'<<"3"<<'\n';
	file_out<<"Vertices"<<'\n'<<this->GetNbVertices()<<'\n';
	if (ref_vertices.GetM() == this->GetNbVertices())
	  {
	    for (int i = 0; i < this->GetNbVertices(); i++)
	      file_out<<this->Vertices(i)(0)<<' '<<this->Vertices(i)(1)
		      <<' '<<this->Vertices(i)(2)<<" " << ref_vertices(i) << "\n";
	  }
	else
	  for (int i = 0; i < this->GetNbVertices(); i++)
	    file_out<<this->Vertices(i)(0)<<' '<<this->Vertices(i)(1)
		    <<' '<<this->Vertices(i)(2)<<" 1\n";
	
	if (nb_quadrangles_ref > 0)
	  {
	    file_out<<"Quadrilaterals"<<'\n'<<nb_quadrangles_ref<<'\n';
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
                if (faces_ref(i).GetNbVertices()==4)
                  file_out<<(faces_ref(i).numVertex(0)+1)<<' '<<(faces_ref(i).numVertex(1)+1)
                          <<' '<<(faces_ref(i).numVertex(2)+1)<<' '
                          <<(faces_ref(i).numVertex(3)+1)<<' '
                          <<(faces_ref(i).GetReference())<<'\n';
	  }
	if (nb_triangles_ref > 0)
	  {
	    file_out<<"Triangles"<<'\n'<<nb_triangles_ref<<'\n';
	    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	      if (this->BoundaryRef(i).GetReference() > 0)
                if (faces_ref(i).GetNbVertices()==3)
                  file_out<<(faces_ref(i).numVertex(0)+1)
                          <<' '<<(faces_ref(i).numVertex(1)+1)<<' '
                          <<(faces_ref(i).numVertex(2)+1)<<' '
                          <<(faces_ref(i).GetReference())<<'\n';
	    
	  }
	
	if (nb_tetras > 0)
	  {
	    file_out<<"Tetrahedra"<<endl<<nb_tetras<<endl;
	    for (int i = 0; i < this->GetNbElt(); i++)
	      if (this->elements(i).GetNbVertices()==4)
		{
		  file_out<<(this->elements(i).numVertex(0)+1)<<' '
                          <<(this->elements(i).numVertex(1)+1)<<' ';
		  file_out<<(this->elements(i).numVertex(2)+1)<<' '
                          <<(this->elements(i).numVertex(3)+1)<<' ';
		  file_out<<this->elements(i).GetReference()<<'\n';
		}
	  }

	if ((nb_hexas > 0)||(nb_wedges > 0)||(nb_pyramids > 0))
	  {
	    // file_out<<"Hexahedra"<<'\n'<<nb_hexas<<'\n';
	    file_out<<"Hexahedra"<<'\n'<<nb_hexas+nb_pyramids+nb_wedges<<'\n';
	    for (int i = 0; i < this->GetNbElt(); i++)
	      {
		if (this->elements(i).GetNbVertices()==8)
		  {
		    file_out<<(this->elements(i).numVertex(0)+1)<<' '
                            <<(this->elements(i).numVertex(1)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(2)+1)<<' '
                            <<(this->elements(i).numVertex(3)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(4)+1)<<' '
                            <<(this->elements(i).numVertex(5)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(6)+1)<<' '
                            <<(this->elements(i).numVertex(7)+1)<<' ';
		    file_out<<this->elements(i).GetReference()<<'\n';
		  }
		
		if (this->elements(i).GetNbVertices()==5)
		  {
		    file_out<<(this->elements(i).numVertex(0)+1)<<' '
                            <<(this->elements(i).numVertex(1)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(2)+1)<<' '
                            <<(this->elements(i).numVertex(3)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(4)+1)<<' '
                            <<(this->elements(i).numVertex(4)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(4)+1)<<' '
                            <<(this->elements(i).numVertex(4)+1)<<' ';
		    file_out<<this->elements(i).GetReference()<<'\n';
		  }
		
		if (this->elements(i).GetNbVertices()==6)
		  {
		    file_out<<(this->elements(i).numVertex(0)+1)<<' '
                            <<(this->elements(i).numVertex(1)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(2)+1)<<' '
                            <<(this->elements(i).numVertex(2)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(3)+1)<<' '
                            <<(this->elements(i).numVertex(4)+1)<<' ';
		    file_out<<(this->elements(i).numVertex(5)+1)<<' '
                            <<(this->elements(i).numVertex(5)+1)<<' ';
		    file_out<<this->elements(i).GetReference()<<'\n';
		  }
	      }
	  }
	
	file_out<<"End"<<endl;
      }
    else
      {
	cout<<"Montjoie does not know how to write this format ."<<extension_output<<endl;
	cout<<"Try to use a different mesh format "<<endl;
	abort();
      }    
  }
    
  
  //! destruction of connectivity arrays (edges and faces)
  template<class T>
  void Mesh<Dimension3, T>::ClearConnectivity()
  {
    Mesh_Base<Dimension3, T>::ClearConnectivity();
    MeshBoundaries<Dimension3, T>::ClearConnectivity();
    faces.Clear();
    
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      faces_ref(i).ClearConnectivity();
    
    next_face.Clear(); head_faces_minv.Fill(-1);
  }
  
  
  //! removes all the mesh
  template<class T>
  void Mesh<Dimension3, T>::Clear()
  {
    Mesh_Base<Dimension3, T>::Clear();
    MeshBoundaries<Dimension3, T>::Clear();
    
    faces.Clear(); faces_ref.Clear();
    next_face.Clear(); head_faces_minv.Clear();
    this->PointsFaceRef.Clear();
    this->PointsEdgeRef.Clear();
    
    for (int i = 0; i < this->GetNbReferences(); i++)
      if (this->GetCurveType(i) == this->CURVE_FILE)
        this->SetCurveType(i, this->NO_CURVE);
    
  }

  
  /******************
   * Mesh Treatment *
   ******************/
  
  
  //! permutes local vertices of each element so that the jacobian Ji is always positive
  /*!
    We consider that FindConnectivity will be called after a call to this method
  */
  template<class T>
  void Mesh<Dimension3, T>::ReorientElements()
  {    
    VectR3 s;
    R3 vec_u, vec_v, vec_w;
    T w; this->affine_element.Reallocate(this->GetNbElt());
    this->reference_element.Reallocate(this->GetNbElt());
    this->affine_element.Fill(false);
    int nb_tet_affine = 0, nb_hex_affine = 0, nb_pyr_affine = 0, nb_prism_affine = 0;
    int nb_tet = 0, nb_hex = 0, nb_pyr = 0, nb_prism = 0;
    for(int i = 0; i < this->GetNbElt(); i++)
      {
	// vertices of the element
	int nb_vertices_loc = this->elements(i).GetNbVertices();
	this->GetVerticesElement(i, s);
        IVect num(nb_vertices_loc);
	for (int j = 0; j < nb_vertices_loc; j++)
	  num(j) = this->elements(i).numVertex(j);
	
	// detecting affine elements
	if (this->Element(i).IsCurved())
	  this->affine_element(i) = false;
	else
	  this->affine_element(i) = ElementGeomReference<Dimension3>::IsAffineTransformation(s);
        
        if (this->affine_element(i))
          {
            switch (nb_vertices_loc)
              {
              case 4:
                nb_tet_affine++; break;
              case 5:
                nb_pyr_affine++; break;
              case 6:
                nb_prism_affine++; break;
              case 8:
                nb_hex_affine++; break;
              }
          }

        switch (nb_vertices_loc)
          {
          case 4:
            nb_tet++; break;
          case 5:
            nb_pyr++; break;
          case 6:
            nb_prism++; break;
          case 8:
            nb_hex++; break;
          }
        
	// we compute the mixed product vec_u ^ vec_v . vec_w (determinant)
	// where vec_u is the vector from the first vertex to the second
	// and vec_v from the first to the third
	// vec_w from the first to the fifth for hexahedra
	//                      to the fourth for tetrahedra
	vec_u = s(1) - s(0);
	if ((nb_vertices_loc == 4)||(nb_vertices_loc == 6))
	  vec_v = s(2) - s(0);
	if ((nb_vertices_loc == 4)||(nb_vertices_loc == 6))
	  vec_w = s(3) - s(0);

	if (nb_vertices_loc == 5)
	  vec_v = s(3) - s(0);
	if (nb_vertices_loc == 5)
	  vec_w = s(4) - s(0);

	if (nb_vertices_loc == 8)
	  vec_v = s(3) - s(0);
	if (nb_vertices_loc == 8)
	  vec_w = s(4) - s(0);
	
	w = Determinant(vec_u, vec_v, vec_w);
	if (w < 0)
	  {
            IVect numb = num;
            int ref = this->Element(i).GetReference();
            
	    // we are in a wrong orientation
	    // we change the local numbering
	    if (nb_vertices_loc == 8)
	      {
		// hexahedron, vertices 1, 3 are swapped and vertices 5, 7
                numb(1) = num(3); numb(3) = num(1);
                numb(5) = num(7); numb(7) = num(5);
		this->Element(i).Init(numb, ref);
              }
	    else if (nb_vertices_loc == 6)
	      {
		// prism
		// vertices 1, 2 are swapped, and vertices 4, 5
                numb(1) = num(2); numb(2) = num(1);
                numb(4) = num(5); numb(5) = num(4);
		this->Element(i).Init(numb, ref);
	      }
	    else if (nb_vertices_loc == 5)
	      {
		// pyramid
		// vertices 1 and 3 are swapped
                numb(1) = num(3); numb(3) = num(1);
                this->Element(i).Init(numb, ref);		
              }
	    else
	      {
		// tetrahedron
		// vertices 1, 2 are swapped
                numb(1) = num(2); numb(2) = num(1);
                this->Element(i).Init(numb, ref);                
              }
	    
	  }
	
	switch (nb_vertices_loc)
	  {
	  case 4:
	    this->reference_element(i) = &tetrahedron_reference; break;
	  case 5:
	    this->reference_element(i) = &pyramid_reference; break;
	  case 6 :
	    this->reference_element(i) = &wedge_reference; break;
	  case 8 :
	    this->reference_element(i) = &hexahedron_reference; break;
	  }
      }
    
    // updating bounding box for the mesh
    this->xmin_ = 1e300; this->xmax_ = -1e300;
    this->ymin_ = 1e300; this->ymax_ = -1e300;
    this->zmin_ = 1e300; this->zmax_ = -1e300;
    for (int i = 0; i < this->GetNbVertices(); i++)
      {
	vec_u = this->Vertex(i);
	this->xmin_ = min(this->xmin_, vec_u(0));
	this->xmax_ = max(this->xmax_, vec_u(0));
	this->ymin_ = min(this->ymin_, vec_u(1));
	this->ymax_ = max(this->ymax_, vec_u(1));
	this->zmin_ = min(this->zmin_, vec_u(2));
	this->zmax_ = max(this->zmax_, vec_u(2));
      }

    if (this->GetNbElt() > 0)
      if (this->print_level >= 4)
        {
          cout << " The mesh contains " << nb_tet_affine 
               << " affine tetrahedra among " << nb_tet << " tetrahedra "<< endl;
          cout << " The mesh contains " << nb_pyr_affine 
               << " affine pyramids among " << nb_pyr << " pyramids "<< endl;
          cout << " The mesh contains " << nb_prism_affine 
               << " affine prisms among " << nb_prism << " prisms "<< endl;
          cout << " The mesh contains " << nb_hex_affine 
               << " affine hexahedra among " << nb_hex << " hexahedra "<< endl;
        }
    // Write("toto.mesh");
  }
  
  
  //! reorientation of referenced faces in counterclockwise order (for a plane surface)
  /*!
    \param[in] normale the normale of each face is reoriented
    in order to be near from this input vector
  */
  template<class T>
  void Mesh<Dimension3, T>::ReorientFaces(int refs, const IVect& ref_cond, const R3& normale)
  {
    VectR3 s(4); IVect num(4);
    R3 u, v, w;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)      
      {
	int ref = faces_ref(i).GetReference();
        if (ref_cond(ref) == refs)
          { 
            int nb_vertices_elt = faces_ref(i).GetNbVertices();
            for (int j = 0; j < nb_vertices_elt; j++)
              {
                num(j) = faces_ref(i).numVertex(j);
                s(j) = this->Vertices(num(j));
              }
            
            u = s(1) - s(0);
            v = s(nb_vertices_elt-1) - s(0);
            
            TimesProd(u, v, w);
            if (DotProd(w, normale) < 0)
              {
                // face badly oriented
                if (nb_vertices_elt == 3)
                  faces_ref(i).InitTriangular(num(0), num(2), num(1), ref);
                else
                  faces_ref(i).InitQuadrangular(num(0), num(3), num(2), num(1), ref);
              }
          }
      }
    
    FindConnectivity();
  }
  

  //! Permutes numbers of a parallelepedic element such that the jacobian matrix will be diagonal
  /*!
    \param[inout] num on input the numbers of the vertices of the element,
    on output the new numbers
    \param[in] ptE vertices of the element
   */
  template<class T>
  void Mesh<Dimension3, T>::ReorientElementCartesian(IVect& num, VectR3& ptE)
  {
    IVect n(num);
    // on calcule DF au point (0, 0, 0)
    R3 dF_dx, dF_dy, dF_dz;
    dF_dx = ptE(1) - ptE(0);
    dF_dy = ptE(3) - ptE(0);
    dF_dz = ptE(4) - ptE(0);
    Real_wp jacob = Determinant(dF_dx, dF_dy, dF_dz);
    
    // on cherche un vecteur oriente suivant ez
    Real_wp epsilon = R3::threshold*Norm2(dF_dx);
    if ((abs(dF_dx(0)) < epsilon) && (abs(dF_dx(1)) < epsilon) )
      {
        if ( (abs(dF_dy(2)) < epsilon) && (abs(dF_dz(2)) < epsilon))
          {
            // on cherche les vecteurs orientes suivant ex, ey
            if ((abs(dF_dy(0)) < epsilon) && (abs(dF_dz(1) < epsilon)))
              {
                // dF_dx = e_z, dF_dy = e_y, dF_dz = e_x
                if (jacob < 0)
                  {
                    num(0) = n(0); num(1) = n(4); num(2) = n(7); num(3) = n(3);
                    num(4) = n(1); num(5) = n(5); num(6) = n(6); num(7) = n(2);
                  }
                else
                  {
                    num(0) = n(4); num(1) = n(0); num(2) = n(3); num(3) = n(7);
                    num(4) = n(5); num(5) = n(1); num(6) = n(2); num(7) = n(6);
                  }
              }
            else
              {
                // dF_dx = e_z, dF_dy = e_x, dF_dz = e_y
                if (jacob < 0)
                  {
                    num(0) = n(3); num(1) = n(0); num(2) = n(4); num(3) = n(7);
                    num(4) = n(2); num(5) = n(1); num(6) = n(5); num(7) = n(6);
                  }
                else
                  {
                    num(0) = n(0); num(1) = n(3); num(2) = n(7); num(3) = n(4);
                    num(4) = n(1); num(5) = n(2); num(6) = n(6); num(7) = n(5);
                  }
              }
          }
      }
    else if ((abs(dF_dy(0)) < epsilon) && (abs(dF_dy(1)) < epsilon) )
      {
        if ( (abs(dF_dx(2)) < epsilon) && (abs(dF_dz(2)) < epsilon))
          {
            // on cherche les vecteurs orientes suivant ex, ey
            if ((abs(dF_dx(0)) < epsilon) && (abs(dF_dz(1) < epsilon)))
              {
                // dF_dx = e_y, dF_dy = e_z, dF_dz = e_x
                if (jacob < 0)
                  {
                    num(0) = n(4); num(1) = n(0); num(2) = n(1); num(3) = n(5);
                    num(4) = n(7); num(5) = n(3); num(6) = n(2); num(7) = n(6);
                  }
                else
                  {
                    num(0) = n(0); num(1) = n(4); num(2) = n(5); num(3) = n(1);
                    num(4) = n(3); num(5) = n(7); num(6) = n(6); num(7) = n(2);
                  }
              }
            else
              {
                // dF_dx = e_x, dF_dy = e_z, dF_dz = e_y
                if (jacob < 0)
                  {
                    num(0) = n(0); num(1) = n(1); num(2) = n(5); num(3) = n(4);
                    num(4) = n(3); num(5) = n(2); num(6) = n(6); num(7) = n(7);
                  }
                else
                  {
                    num(0) = n(1); num(1) = n(0); num(2) = n(4); num(3) = n(5);
                    num(4) = n(2); num(5) = n(3); num(6) = n(7); num(7) = n(6);
                  }
              }
          }
      }
    else if ((abs(dF_dz(0)) < epsilon) && (abs(dF_dz(1)) < epsilon) )
      {
        if ( (abs(dF_dx(2)) < epsilon) && (abs(dF_dy(2)) < epsilon))
          {
            // on cherche les vecteurs orientes suivant ex, ey
            if ((abs(dF_dx(0)) < epsilon) && (abs(dF_dz(1) < epsilon)))
              {
                // dF_dx = e_y, dF_dy = e_x, dF_dz = e_z
                if (jacob < 0)
                  {
                    num(0) = n(0); num(1) = n(3); num(2) = n(2); num(3) = n(1);
                    num(4) = n(4); num(5) = n(7); num(6) = n(6); num(7) = n(5);
                  }
                else
                  {
                    num(0) = n(3); num(1) = n(0); num(2) = n(1); num(3) = n(2);
                    num(4) = n(7); num(5) = n(4); num(6) = n(5); num(7) = n(6);
                  }
              }
            else
              {
                // dF_dx = e_x, dF_dy = e_y, dF_dz = e_z
                if (jacob < 0)
                  {
                    num(0) = n(1); num(1) = n(0); num(2) = n(3); num(3) = n(2);
                    num(4) = n(5); num(5) = n(4); num(6) = n(7); num(7) = n(6);
                  }
              }
          }
      }    
  }
  
  
  //! we add a single face to the mesh while updating connectivity
  template<class T>
  void Mesh<Dimension3, T>::AddSingleFace(int& nb_faces, const Face<Dimension3>& face)
  {
    int k0 = 0;
    int n1 = face.numVertex(0);
    int nb_vert = face.GetNbVertices();
    IVect num(nb_vert); int ref = face.GetReference();
    IVect num_edge(nb_vert);
    for (int k = 1; k < nb_vert; k++)
      if (face.numVertex(k) < n1)
	{
	  n1 = face.numVertex(k);
	  k0 = k;
	}
    
    for (int k = 0; k < nb_vert; k++)
      {
	num(k) = face.numVertex((k0+k)%nb_vert);
	num_edge(k) = face.numEdge((k0+k)%nb_vert);
      }
    
    bool existe = false;
    // we look at the list of edges leaving from n1
    // if we find the same edge than n1->n2 we add the edge
    for (int e = head_faces_minv(n1); e != -1; e = next_face(e))
      {
	if (faces(e).GetNbVertices() == nb_vert)
	  {
	    bool all_vert_equal = true;
	    for (int k = 1; k < nb_vert; k++)
	      if (faces(e).numVertex(k) != num(k))
		all_vert_equal = false;
	    
	    if (!all_vert_equal)
	      {
		all_vert_equal = true;
		for (int k = 1; k < nb_vert; k++)
		  if (faces(e).numVertex(k) != num(nb_vert-k))
		    all_vert_equal = false;
	      }
	    
	    if (all_vert_equal)
	      existe = true;
	  }
      }
    
    // if the face has not been found, we create it
    if (!existe) 
      {
	// we initialize the face
	faces(nb_faces).Init(num, ref);
	for (int k = 0; k < nb_vert; k++)
          faces(nb_faces).SetEdge(k, num_edge(k));
        
	// we update the arrays next_edge, head_minv, etc
	next_face(nb_faces) = head_faces_minv(n1);
	head_faces_minv(n1) = nb_faces;
	nb_faces++;
      }
  }
  
  
  //! we add a face belonging to an element while updating connectivity of the mesh
  template<class T>
  void Mesh<Dimension3, T>::AddFaceOfElement(int& nb_faces, int i, int j)
  {
    int type_elt = this->Element(i).GetHybridType();
    // quadrilateral_face ?
    int nb_vertices_face = 3;
    bool quad_face = MeshNumbering<Dimension3, T>::IsLocalFaceQuadrilateral(type_elt, j);
    if (quad_face)
      nb_vertices_face = 4;
    
    // we get the four local numbers of the face
    int iloc = MeshNumbering<Dimension3, T>::FirstExtremityFace(type_elt, j);
    int jloc = MeshNumbering<Dimension3, T>::SecondExtremityFace(type_elt, j);
    int kloc = MeshNumbering<Dimension3, T>::ThirdExtremityFace(type_elt, j);
    
    // and global numbers
    int iglob = this->Element(i).numVertex(iloc);
    int jglob = this->Element(i).numVertex(jloc);
    int kglob = this->Element(i).numVertex(kloc);
    
    int iglob2 = iglob, jglob2 = jglob, kglob2 = kglob, lglob2 = -1;
    int rot = 0;
    // sorting numbers
    if (quad_face)
      {
	int lloc = MeshNumbering<Dimension3, T>::FourthExtremityFace(type_elt, j);
	int lglob = this->Element(i).numVertex(lloc);
	lglob2 = lglob;	
	if (jglob < min(min(iglob, kglob), lglob))
	  {
	    iglob2 = jglob;
	    jglob2 = kglob;
	    kglob2 = lglob;
	    lglob2 = iglob;
	    rot = 3;
	  }
	else if (kglob < min(min(iglob, jglob), lglob))
	  {
	    iglob2 = kglob;
	    jglob2 = lglob;
	    kglob2 = iglob;
	    lglob2 = jglob;
	    rot = 2;
	  } 
	else if (lglob < min(min(iglob, jglob), kglob))
	  {
	    iglob2 = lglob;
	    jglob2 = iglob;
	    kglob2 = jglob;
	    lglob2 = kglob;
	    rot = 1;
	  }
      }
    else
      {
	if (jglob < min(iglob, kglob))
	  {
	    iglob2 = jglob;
	    jglob2 = kglob;
	    kglob2 = iglob;
	    rot = 2;
	  }
	else if (kglob < min(iglob, jglob))
	  {
	    iglob2 = kglob;
	    jglob2 = iglob;
	    kglob2 = jglob;
	    rot = 1;
	  }
      }
    
    bool existe = false;
    // let us search if this face has already been created
    // we do a loop on all the faces, which leaves vertex iglob2
    for (int e = head_faces_minv(iglob2); e != -1; e = next_face(e))
      if (faces(e).GetNbVertices() == nb_vertices_face)
	{
	  // the four vertices of face e
	  int n1 = faces(e).numVertex(1);
	  int n2 = faces(e).numVertex(2);
	  bool face_equal = false;
	  if (quad_face)
	    {
	      int n3 = faces(e).numVertex(3);
	      if ((n1 == jglob2)&&(n2==kglob2)&&(n3==lglob2))
		face_equal = true;
	      
	      if ((n3 == jglob2)&&(n2==kglob2)&&(n1==lglob2))
		face_equal = true;
	    }
	  else
	    {
	      if ((n1 == jglob2)&&(n2==kglob2))
		face_equal = true;
	      
	      if ((n2 == jglob2)&&(n1==kglob2))
		face_equal = true;
	    }
	  
	  
	  if (face_equal)
	    {
	      // they are equal -> the face is already created
	      // we update the connectivity of the face and element
	      faces(e).AddElement(i);
	      this->elements(i).SetFace(j, e);
	      this->elements(i).SetOrientationFace(j, 0);
	      
	      if (quad_face)
		{
		  // we set the way of face
		  // same way except a rotation
		  if ((faces(e).numVertex(1) == iglob)&&(faces(e).numVertex(2) == jglob))
		    this->elements(i).SetOrientationFace(j,1);
		  else if ((faces(e).numVertex(2) == iglob)&&(faces(e).numVertex(3) == jglob))
		    this->elements(i).SetOrientationFace(j,2);
		  else if ((faces(e).numVertex(3) == iglob)&&(faces(e).numVertex(0) == jglob))
		    this->elements(i).SetOrientationFace(j,3);
		  
		  // opposite way and eventually a rotation
		  if ((faces(e).numVertex(0) == iglob)&&(faces(e).numVertex(3) == jglob))
		    this->elements(i).SetOrientationFace(j,4);
		  else if ((faces(e).numVertex(1) == iglob)&&(faces(e).numVertex(0) == jglob))
		    this->elements(i).SetOrientationFace(j,5);
		  else if ((faces(e).numVertex(2) == iglob)&&(faces(e).numVertex(1) == jglob))
		    this->elements(i).SetOrientationFace(j,6);
		  else if ((faces(e).numVertex(3) == iglob)&&(faces(e).numVertex(2) == jglob))
		    this->elements(i).SetOrientationFace(j,7);		  
		}
	      else
		{
		  // we set the way of face
		  // same way except a rotation
		  if ((faces(e).numVertex(1) == iglob)&&(faces(e).numVertex(2) == jglob))
		    this->elements(i).SetOrientationFace(j,1);
		  else if ((faces(e).numVertex(2) == iglob)&&(faces(e).numVertex(0) == jglob))
		    this->elements(i).SetOrientationFace(j,2);
		  
		  // opposite way and eventually a rotation
		  if ((faces(e).numVertex(0) == iglob)&&(faces(e).numVertex(2) == jglob))
		    this->elements(i).SetOrientationFace(j,3);
		  else if ((faces(e).numVertex(1) == iglob)&&(faces(e).numVertex(0) == jglob))
		    this->elements(i).SetOrientationFace(j,4);
		  else if ((faces(e).numVertex(2) == iglob)&&(faces(e).numVertex(1) == jglob))
		    this->elements(i).SetOrientationFace(j,5);
		  
		}
	      existe = true;
	      break;
	    }
	}
	    
    if (!existe)
      {
	// the face is created
	faces(nb_faces).AddElement(i);
	// initialization of the quadrangular face
	// the last number is the reference
	if (quad_face)
	  faces(nb_faces).InitQuadrangular(iglob2, jglob2, kglob2, lglob2, 0);
	else
	  faces(nb_faces).InitTriangular(iglob2, jglob2, kglob2, 0);
	
	this->elements(i).SetFace(j, nb_faces); // we add it to the element i
	this->elements(i).SetOrientationFace(j, rot);
	// the arrays next_face and head_faces_minv are updated
	next_face(nb_faces) = head_faces_minv(iglob2);
	head_faces_minv(iglob2) = nb_faces; 
	// we increment the number of faces
	nb_faces++;
      }
  }
  
  
  //! Find all the faces of the mesh, even inside
  /*!
    Hecht's algorithm (quasi-optimal)
    arrays head_faces_minv and next_face are created
  */
  template<class T>
  void Mesh<Dimension3, T>::AddFacesFromElements()
  {
    int nb_old = this->GetNbFaces();
    // a majorant for the number of faces
    int nb_faces = GetUpperBoundNumberOfFaces(); 
    
    this->ResizeFaces(nb_faces + nb_old);
    
    nb_faces = nb_old;
    // loop over elements of the mesh
    for (int i = 0; i < this->GetNbElt(); i++)
      {
	int nb_faces_loc = this->elements(i).GetNbFaces();
	// loop over all the faces of element i
	for (int j = 0; j < nb_faces_loc; j++)
	    AddFaceOfElement(nb_faces, i, j);
	
      }
    
    this->ResizeFaces(nb_faces);
    
    if (this->print_level >= 6)
      cout << " Faces constructed "<<endl;
    
    if (this->print_level >= 6)
      cout << " Number of faces in the mesh  "<<nb_faces<<endl;
  }
  
  
  //! we recompute connectivity of the mesh (relation between elements, faces and edges)
  template<class T>
  void Mesh<Dimension3, T>::FindConnectivity()
  {
    ClearConnectivity();
        
    // we put referenced edges
    int nb_edges = 0;
    this->ReallocateEdges(this->GetNbEdgesRef());
    for (int i = 0; i < this->GetNbEdgesRef(); i++)
      this->AddSingleEdge(nb_edges, this->EdgeRef(i));
    
    if (nb_edges < this->GetNbEdges())
      this->ResizeEdges(nb_edges);
    
    // then we add edges from referenced faces
    this->AddEdgesFromReferencedFaces();
    
    // and edges from elements
    this->AddEdgesFromElements();
    
    // we put referenced faces
    int nb_faces = 0;
    this->ReallocateFaces(this->GetNbBoundaryRef());
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      this->AddSingleFace(nb_faces, this->BoundaryRef(i));
    
    if (nb_faces < this->GetNbFaces())
      this->ResizeFaces(nb_faces);
    
    // then we add faces from elements
    this->AddFacesFromElements();
    
    this->CreateConnexionEdgeToFace();

    this->ConstructCrossReferenceEdgeRef();
    this->ConstructCrossReferenceBoundaryRef();
  }
  
  
  //! we add an edge belonging to a face, while updating connectivity
  template<class T>
  void Mesh<Dimension3, T>::AddEdgeOfFace(int& nb_edges, int i, int na)
  {
    int nb_edges_loc = this->BoundaryRef(i).GetNbEdges();
    
    // local numbers of the two extremities
    int iloc = na;
    int jloc = (na+1)%nb_edges_loc;
    
    // we get global numbers of the two extremities of the local edge
    int iglob = this->BoundaryRef(i).numVertex(iloc);
    int jglob = this->BoundaryRef(i).numVertex(jloc);
    
    // We consider the edge iglob -> jglob 
    // we sort the numbers
    if (iglob>jglob)
      {
	int itmp = iglob;
	iglob = jglob;
	jglob = itmp;
      }
    
    bool existe = false;
    // we look at the list of edges leaving from iglob
    // if we find the same edge than iglob->jglob we add the face i
    // to this edge
    for (int e = this->head_minv(iglob); e != -1; e = this->next_edge(e))
      {
	// the edge iglob->jglob already exists
	// we update edges(e) and faces(i)
	if ( this->edges(e).numVertex(1) == jglob)
	  {
	    this->BoundaryRef(i).SetEdge(na, e);
            this->GetEdge(e).AddFace(i);
	    existe = true;
	  }
      }
    
    // if the edge has not been found, we create it
    if (!existe) 
      {
	// we initialize the edge 
	this->edges(nb_edges).Init(iglob,jglob,0);
	this->BoundaryRef(i).SetEdge(na, nb_edges); // we add it to the face i
        this->GetEdge(nb_edges).AddFace(i);
		
	// we update the arrays next_edge, head_minv, next_edge_arr and head_maxv
	this->next_edge(nb_edges) = this->head_minv(iglob);
	this->head_minv(iglob) = nb_edges;
	this->next_edge_arr(nb_edges) = this->head_maxv(jglob);
	this->head_maxv(jglob) = nb_edges;
	nb_edges++;
      }
  }
  
  
  //! we find the edges of the referenced faces
  /*!
    The inside edges owning to volumes are ignored
    Only referenced faces are considered
  */
  template<class T>
  void Mesh<Dimension3, T>::AddEdgesFromReferencedFaces()
  {
    int nb_old = this->GetNbEdges();
    // we compute a max value for the number of edges in the mesh
    int nb_edges = 4*GetNbBoundaryRef();
    
    // changing size of arrays
    this->ResizeEdges(nb_edges + nb_old);
    
    // loop over all referenced faces
    nb_edges = nb_old;
    for(int i = 0; i < this->GetNbBoundaryRef(); i++)
      {
	int nb_edges_loc = faces_ref(i).GetNbEdges();
	
	// loop on local edges of the face i
	for (int na = 0; na < nb_edges_loc; na++) 
	  this->AddEdgeOfFace(nb_edges, i, na);
      }
    
    this->ResizeEdges(nb_edges);
    
  } 
  
  
  //! fills connectivity arrays edges <-> faces
  template<class T>
  void Mesh<Dimension3, T>::CreateConnexionEdgeToFace()
  {    
    // loop on the faces
    IVect numero_edge(12),num_vertex1(12),num_vertex2(12);
    for (int i = 0; i < this->GetNbFaces(); i++)
      {
	// we get the list of the edges around the face i
	int nb_edges_around = 0;
	int num_elem = faces(i).numElement(0);
	if (num_elem >= 0)
	  {
	    int nb_edges_elem = this->elements(num_elem).GetNbEdges();
	    for (int k = 0; k < nb_edges_elem; k++)
	      {
		int num_edge = this->elements(num_elem).numEdge(k);
		numero_edge(nb_edges_around) = num_edge;
		num_vertex1(nb_edges_around) = this->edges(num_edge).numVertex(0);
		num_vertex2(nb_edges_around) = this->edges(num_edge).numVertex(1);
		nb_edges_around++;
	      }
	    
	    // now we identify edges of the face with edges of the list numero_edge
	    int nb_vertex_face = faces(i).GetNbVertices();
	    for (int j = 0; j < nb_vertex_face; j++)
	      {
		// two extremities of the local edge
		int n0 = faces(i).numVertex(j);
		int n1 = faces(i).numVertex( (j+1)%nb_vertex_face );
		// we sort the two extremities
		if (n0 > n1)
		  {
                    int ntmp = n0;
		    n0 = n1;
		    n1 = ntmp;
		  }
		
		bool trouve = false;
		for (int k = 0; k < nb_edges_around; k++)
		  {
		    if ((num_vertex1(k) == n0)&&(num_vertex2(k) == n1))
		      {
			faces(i).SetEdge(j, numero_edge(k));
			this->edges(numero_edge(k)).AddFace(i);
			trouve = true;
		      }
		  }
		
		if (!trouve)
		  {
		    Write("test.mesh");
		    // DISP(num_elem); DISP(elements(num_elem)); DISP(i); DISP(faces(i));
		    // we didn't find the edge
		    cout<<"Faces or edges are wrong. Check your mesh "<<endl;
		    abort();
		  }
	      }
	  }
      }
  }

  
  //! symmetrize the mesh (you can choose to keep or not referenced faces 
  //! which are on the symmetry plane)
  /*!
    The symmetry is performed with respect to the plane x = xc if num_coor=0,
    y = xc if num_coor=1 or z = xc if num_coor = 2
    if num_coor is equal to -1, the symmetry plane is
    a x + b y + c z + d = where (a, b, c, d) are stored in normale and xc
    if keep_face_on_plane is true, the referenced faces belonging to the symmetry plane are kept,
    otherwise they are removed
   */
  template<class T>
  void Mesh<Dimension3, T>::SymmetrizeMesh(int num_coor, T xc,
                                           TinyVector<T, 3> normale,
                                           bool keep_face_on_plane)
  {
    // normalizing normale
    T coef = Norm2(normale);
    if (num_coor == -1)
      {
        coef = T(1)/coef;
        normale *= coef;
        xc *= coef;
      }
    
    // counting vertices
    int nb_old_vertices = this->GetNbVertices();
    // IndexVert(i) -> new vertex number
    IVect IndexVert(this->GetNbVertices()); IndexVert.Fill(-1);
    int nb_vertices = nb_old_vertices;
    for (int i = 0; i < nb_old_vertices; i++)
      {
	// if vertex on symmetry plane, we keep the same number (no new vertex created)
	if (num_coor == -1)
          {
            if (abs(DotProd(this->Vertices(i), normale) + xc) > R3::threshold)
              IndexVert(i) = nb_vertices++;
            else
              IndexVert(i) = i;
          }
        else
          {
            if (abs(this->Vertices(i)(num_coor)-xc) > R3::threshold)
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
    
    // now we count referenced edges
    int nb_old_edges_ref = this->GetNbEdgesRef();
    int nb_edges_ref = nb_old_edges_ref;
    for (int i = 0; i < nb_old_edges_ref; i++)
      {
        int n1 = this->EdgeRef(i).numVertex(0);
        int n2 = this->EdgeRef(i).numVertex(1);
        if ((IndexVert(n1) < nb_old_vertices)&&(IndexVert(n2) < nb_old_vertices))
          {
            // edge on symmetry plane, no new edge
          }
        else
          nb_edges_ref++;
      }
    
    // now we create new edges
    this->ResizeEdgesRef(nb_edges_ref);
    nb_edges_ref = nb_old_edges_ref;
    for (int i = 0; i < nb_old_edges_ref; i++)
      {
        int n1 = this->EdgeRef(i).numVertex(0);
        int n2 = this->EdgeRef(i).numVertex(1);
        if ((IndexVert(n1) < nb_old_vertices)&&(IndexVert(n2) < nb_old_vertices))
          {
            // edge on symmetry plane, no new edge
          }
        else
          {
            this->EdgeRef(nb_edges_ref).Init(IndexVert(n1), IndexVert(n2), 
                                             this->EdgeRef(i).GetReference());
            
            for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
              {
                this->PointsEdgeRef(nb_edges_ref, k) = this->PointsEdgeRef(i, k);
                if (num_coor == -1)
                  SymmetrizePointPlane(this->PointsEdgeRef(i, k), normale, xc,
                                       this->PointsEdgeRef(nb_edges_ref, k));
                else
                  this->PointsEdgeRef(nb_edges_ref, k)(num_coor)
                    = 2.0*xc - this->PointsEdgeRef(nb_edges_ref, k)(num_coor);
              }
            
            nb_edges_ref++;
          }
      }

    // counting referenced faces
    int nb_old_faces_ref = this->GetNbBoundaryRef();
    int nb_faces_ref = nb_old_faces_ref;
    IVect num(4); 
    for (int i = 0; i < nb_old_faces_ref; i++)
      {
        int nb_vert = this->BoundaryRef(i).GetNbVertices();
        bool face_on_plane = true;
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = IndexVert(this->BoundaryRef(i).numVertex(j));
            if (num(j) >= nb_old_vertices)
              face_on_plane = false;
          }
        
        if (face_on_plane)
          nb_faces_ref++;                
      }
    
    // we create new referenced faces
    Vector<Vector<R_N> > OldPointsFace = this->PointsFaceRef;
    Vector<Face<Dimension3> > old_faces = this->faces_ref;
    ReallocateBoundariesRef(2*this->GetNbBoundaryRef());
    nb_faces_ref = 0;
    for (int i = 0; i < nb_old_faces_ref; i++)
      {
        int nb_vert = old_faces(i).GetNbVertices();
        bool face_on_plane = true;
        num.Reallocate(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = IndexVert(old_faces(i).numVertex(j));
            if (num(j) >= nb_old_vertices)
              face_on_plane = false;
          }
        
        if (!face_on_plane)
          {
            // initial face
            this->BoundaryRef(nb_faces_ref) = old_faces(i);
            this->PointsFaceRef(nb_faces_ref) = OldPointsFace(i);
            nb_faces_ref++;
            
            // symmetrized face
            this->BoundaryRef(nb_faces_ref).Init(num, old_faces(i).GetReference());
            this->PointsFaceRef(nb_faces_ref) = OldPointsFace(i);
            if (num_coor == -1)
              for (int k = 0; k < OldPointsFace(i).GetM(); k++)
                SymmetrizePointPlane(OldPointsFace(i)(k), normale, xc,
                                     this->PointsFaceRef(nb_faces_ref)(k));
            else
              for (int k = 0; k < OldPointsFace(i).GetM(); k++)
                this->PointsFaceRef(nb_faces_ref)(k)(num_coor) 
                  = 2.0*xc - this->PointsFaceRef(nb_faces_ref)(k)(num_coor);
            
            nb_faces_ref++;
          }
        
        if (keep_face_on_plane)
          {
            // initial face
            this->BoundaryRef(nb_faces_ref) = old_faces(i);
            this->PointsFaceRef(nb_faces_ref) = OldPointsFace(i);
            nb_faces_ref++;
          }
      }
    
    ResizeBoundariesRef(nb_faces_ref);
        
    // creating symmetrized elements
    int nb_elt = this->GetNbElt();
    this->ResizeElements(2*nb_elt);
    for (int i = 0; i < nb_elt; i++)
      {
	int nb_vert = this->elements(i).GetNbVertices();
        num.Reallocate(nb_vert);
	for (int j = 0; j < nb_vert; j++)
	  num(j) = IndexVert(this->Element(i).numVertex(j));
	
	this->Element(nb_elt + i).Init(num, this->Element(i).GetReference());
      }
    
    // and we reorient elements and find connectivity
    this->ReorientElements();
    this->FindConnectivity();
    this->ProjectPointsOnCurves();
  }
  
  
  //! returns a subdivision of the current mesh through points and numbering
  /*!
    \param[in] outside_subdiv 1-D subdivision, including extremities 
    (ie 0  0.2 0.4 0.6 0.8  1.0, for a regular subdivision in 5 intervals)
    \param[in] points_div subdivision for each reference element (0 -> tet, 1 -> pyramid, etc)
    \param[in] point_surf subdivision for unit triangle and quadrangle
    \param[out] PosNodes vertices of the subdivided mesh
    \param[out] Nodle numbering of the subdivided mesh
    \param[out] NodleSurf surface numbering of the subdivided mesh
  */
  template<class T>
  void Mesh<Dimension3, T>::
  GetMeshSubdivision(const Vector<T>& outside_subdiv,
                     const Vector<Vector<R_N> >& points_div,
                     const Vector<VectR2>& points_surf, Vector<R_N>& PosNodes,
                     Vector<IVect>& Nodle, Vector<IVect>& NodleSurf) const
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
	    int nv = this->BoundaryRef(i).GetNbVertices();
	    NodleSurf(i).Reallocate(nv);
	    for (int j = 0; j < nv; j++)
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
    int nb_inside_tri = (r-1)*(r-2)/2, nb_tri = (r+1)*(r+2)/2;
    int nb_inside_quad = (r-1)*(r-1), nb_quad = (r+1)*(r+1);
    int nb_inside_tet = (r-1)*(r-2)*(r-3)/6, nb_tet = (r+1)*(r+2)*(r+3)/6;
    int nb_inside_pyr = (r-1)*(r-2)*(2*r-3)/6, nb_pyr = (r+1)*(r+2)*(2*r+3)/6;;
    int nb_inside_wed = (r-1)*nb_inside_tri, nb_wed = (r+1)*nb_tri;
    int nb_inside_hex = (r-1)*nb_inside_quad, nb_hex = (r+1)*nb_quad;
    int nb_nodes = this->GetNbVertices() + this->GetNbEdges()*(r-1)
      + this->GetNbTriangles()*nb_inside_tri
      + this->GetNbQuadrangles()*(r-1)*(r-1) + this->GetNbTetrahedra()*nb_inside_tet
      + this->GetNbPyramids()*nb_inside_pyr
      + this->GetNbWedges()*nb_inside_wed + this->GetNbHexahedra()*nb_inside_hex;
    
    PosNodes.Reallocate(nb_nodes);
    SetPoints<Dimension3> PointsElem;
    Vector<R_N> s, InterpPoints;
    
    // constructing interpolation to points_surf
    FiniteElementInterpolator interp_surf;
    Vector<ElementGeomReference<Dimension2>* > elt_surf(2);
    elt_surf(0) = &const_cast<TriangleGeomReference&>(this->triangle_reference);
    elt_surf(1) = &const_cast<QuadrangleGeomReference&>(this->quadrangle_reference);

    // constructing 2-D points
    interp_surf.InitProjection(elt_surf, outside_subdiv, points_surf);

    Matrix<int> FaceRotationTri, FaceRotationQuad;
    VectR2 Points2D(nb_inside_tri);
    for (int i = 0; i < nb_inside_tri; i++)
      Points2D(i) = points_surf(0)(3*r+i);
    
    MeshNumbering<Dimension3, T>::GetRotationTriangularFace(Points2D, FaceRotationTri);
    
    Points2D.Reallocate(nb_inside_quad);
    for (int i = 0; i < nb_inside_quad; i++)
      Points2D(i) = points_surf(1)(4*r+i);
    
    MeshNumbering<Dimension3, T>::GetRotationQuadrilateralFace(Points2D, FaceRotationQuad);
    
    // constructing interpolation to points_div
    FiniteElementInterpolator interp_vol;
    
    // constructing 3-D points
    interp_vol.InitProjection(this->reference_element, outside_subdiv, points_div);    
    
    // first we copy vertices
    for (int i = 0; i < this->GetNbVertices(); i++)
      PosNodes(i) = this->Vertex(i);
    
    // 1-D interpolation
    int order = this->GetGeometryOrder();
    Matrix<Real_wp> ValPhi(order+1, r-1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j < r-1; j++)
	ValPhi(i, j) = this->lob_curve.EvaluatePhi(i, step_subdiv(j));

    // computation of nodes for all the referenced edges
    VectR3 OldPoints(order+1);;
    for (int i = 0; i < this->GetNbEdgesRef(); i++)
      {
	OldPoints(0) = this->Vertex(this->EdgeRef(i).numVertex(0));
	OldPoints(order) = this->Vertex(this->EdgeRef(i).numVertex(1));
	for (int k = 1; k < order; k++)
	  OldPoints(k) = this->PointsEdgeRef(i, k-1);
	
	int offset = this->GetNbVertices() + i*(r-1);
	for (int k = 0; k < r-1; k++)
	  {
	    PosNodes(offset + k).Zero();
	    for (int j = 0; j <= order; j++)
	      Add(ValPhi(j, k), OldPoints(j), PosNodes(offset+k));	    
	  }
      }
    
    // construction of NodleSurf and computation of surfacic nodes
    NodleSurf.Reallocate(this->GetNbBoundaryRef());
    int offset = this->GetNbVertices();
    int nb = offset + (r-1)*this->GetNbEdges();
    IVect OffsetFace(this->GetNbBoundary()+1); OffsetFace(0) = nb;
    for (int i = 0; i < this->GetNbBoundary(); i++)
      {
	if (this->Boundary(i).GetNbVertices() == 3)
	  OffsetFace(i+1) = OffsetFace(i) + nb_inside_tri;
	else
	  OffsetFace(i+1) = OffsetFace(i) + nb_inside_quad;
      }
	
    VectR3 PtEdge(order-1);
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      {	
	int nb_dof = 0;
	int nb_vert = this->BoundaryRef(i).GetNbVertices();
	if (nb_vert == 3)
	  {
            nb_dof = nb_tri;
            s.Reallocate((order+1)*(order+2)/2);
          }
	else
	  {
            nb_dof = nb_quad;
            s.Reallocate((order+1)*(order+1));
          }
	
	NodleSurf(i).Reallocate(nb_dof);
	// vertices
	for (int j = 0; j < nb_vert; j++)
	  {
	    NodleSurf(i)(j) = this->BoundaryRef(i).numVertex(j);
	    s(j) = this->Vertex(NodleSurf(i)(j));
	  }
	
	// then edges
	for (int j = 0; j < nb_vert; j++)
	  {
	    int ne = this->BoundaryRef(i).numEdge(j);
	    if (ne >= this->GetNbEdgesRef())
	      {
		R3 ptA = this->Vertex(this->GetEdge(ne).numVertex(0));
		R3 ptB = this->Vertex(this->GetEdge(ne).numVertex(1));
		for (int k = 0; k < order-1; k++)
		  PtEdge(k) = (1.0-this->lob_curve.Points(k+1))*ptA 
                    + this->lob_curve.Points(k+1)*ptB;
	      }
	    else
	      for (int k = 0; k < order-1; k++)
		PtEdge(k) = this->PointsEdgeRef(ne, k);
	    
	    if (this->BoundaryRef(i).GetOrientationEdge(j))
	      {
		for (int k = 0; k < r-1; k++)
		  NodleSurf(i)(nb_vert + j*(r-1) + k) = offset + ne*(r-1) + k;
		
		for (int k = 0; k < order-1; k++)
		  s(nb_vert + j*(order-1) + k) = PtEdge(k);
	      }
	    else
	      {
		for (int k = 0; k < r-1; k++)
		  NodleSurf(i)(nb_vert + j*(r-1) + k) = offset + ne*(r-1) + r - 2 - k;
		
		for (int k = 0; k < order-1; k++)
		  s(nb_vert + j*(order-1) + k) = PtEdge(order-2-k);
	      }
	  }
	
	// the interior of face
	for (int k = nb_vert*r; k < nb_dof; k++)
	  NodleSurf(i)(k) = OffsetFace(i) + k-nb_vert*r;
	
	int nb_pts_bound = nb_vert*order;
	for (int k = 0; k < this->PointsFaceRef(i).GetM(); k++)
	  s(nb_pts_bound + k) = this->PointsFaceRef(i)(k);
	
	// interpolation
	interp_surf.Project(s, InterpPoints, nb_vert-3);
	
	for (int k = 0; k < nb_dof; k++)
	  PosNodes(NodleSurf(i)(k)) = InterpPoints(k);
        
      }
    
    // construction of Nodle and computation of nodes
    Nodle.Reallocate(this->GetNbElt());
    int offset_vol = OffsetFace(this->GetNbBoundary());
    for (int i = 0; i < this->GetNbElt(); i++)
      {	
	int nb_dof = 0;
	int nb_vert = this->Element(i).GetNbVertices();
	int nb_edges = this->Element(i).GetNbEdges();
        switch (nb_vert)
	  {
	  case 4:
	    nb_dof = nb_tet;
            break;
	  case 5 :
	    nb_dof = nb_pyr;
            break;
	  case 6 :
	    nb_dof = nb_wed;
            break;
	  case 8 :
	    nb_dof = nb_hex;
            break;
	  }
	
	Nodle(i).Reallocate(nb_dof);
	// vertices
	for (int j = 0; j < nb_vert; j++)
          Nodle(i)(j) = this->Element(i).numVertex(j);
        
	// computing the points
	this->GetVerticesElement(i, s);
	this->FjElemNodal(s, PointsElem, *this, i);
	
	// interpolation
	interp_vol.Project(PointsElem.GetPointNodal(), InterpPoints, this->GetTypeElement(i));
  	
	// then edges
	for (int j = 0; j < nb_edges; j++)
	  {
	    int ne = this->Element(i).numEdge(j);
            
	    if (this->Element(i).GetOrientationEdge(j))
	      {
                for (int k = 0; k < r-1; k++)
                  Nodle(i)(nb_vert + j*(r-1) + k) = offset + ne*(r-1) + k;
              }
	    else
	      {
                for (int k = 0; k < r-1; k++)
                  Nodle(i)(nb_vert + j*(r-1) + k) = offset + ne*(r-1) + r - 2 - k;
              }            
	  }
	
	// then faces
        int offset_elt = nb_vert + (r-1)*nb_edges;
        for (int j = 0; j < this->Element(i).GetNbFaces(); j++)
          {
	    int ne = this->Element(i).numFace(j);
            int rot = this->Element(i).GetOrientationFace(j);
            if (this->Boundary(ne).GetNbVertices() == 3)
              {
                for (int k = 0; k < nb_inside_tri; k++)
                  {
                    int krot = FaceRotationTri(rot, k);
                    Nodle(i)(offset_elt + krot) = OffsetFace(ne) + k;
                  }
                
                if (ne < this->GetNbBoundaryRef())
                  {
                    for (int k = 0; k < nb_inside_tri; k++)
                      if (Nodle(i)(offset_elt+FaceRotationTri(rot, k)) != NodleSurf(ne)(k+3*r))
                        {
                          DISP(i); DISP(ne); DISP(rot); DISP(offset_elt);
                          DISP(Nodle(i));
                          DISP(Nodle(i)(offset_elt+FaceRotationTri(rot, k)));
                          DISP(NodleSurf(ne)(k+3*r));
			  cout << "Inconsistent numbering in the mesh" << endl;
                          abort();
                        }
                  }
                
                offset_elt += nb_inside_tri;
              }
	    else
	      {
                for (int k = 0; k < nb_inside_quad; k++)
                  {
                    int krot = FaceRotationQuad(rot, k);
                    Nodle(i)(offset_elt + krot) = OffsetFace(ne) + k;
                  }
                
                offset_elt += nb_inside_quad;
              }
            
	  }
        
        // then volumes
        for (int j = offset_elt; j < nb_dof; j++)
          Nodle(i)(j) = offset_vol++;        
        
	for (int k = 0; k < nb_dof; k++)
	  PosNodes(Nodle(i)(k)) = InterpPoints(k);	
      }
    
  }  
  
    
  //! we find vertices with the same coordinates, and remove duplicates
  template<class T>
  void Mesh<Dimension3, T>::RemoveDuplicateVertices()
  {
    // vertices with the same positions are removed
    int nb_vertices = this->GetNbVertices();
    IVect permutation(nb_vertices); permutation.Fill();
    // new vertex numbers
    IVect new_num(nb_vertices); new_num.Fill(-1);    
    // vertices are sorted
    Sort(nb_vertices, this->Vertices, permutation);
    
    int nb = 1; R3 prec = this->Vertex(0);
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
	  new_num(permutation(i)) = nb-1;
      }
    
    this->ResizeVertices(nb);
    
    // updating edge numbers
    for (int i = 0; i < this->GetNbEdgesRef(); i++)
      {
        int n1 = new_num(this->EdgeRef(i).numVertex(0));
        int n2 = new_num(this->EdgeRef(i).numVertex(1));
        int ref = this->EdgeRef(i).GetReference();
        this->EdgeRef(i).Init(n1, n2, ref);
      }
    
    // updating face numbers
    IVect num;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      {
        int nb_vert = this->BoundaryRef(i).GetNbVertices();
        int ref = this->BoundaryRef(i).GetReference();
        num.Reallocate(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          num(j) = new_num(this->BoundaryRef(i).numVertex(j));
        
        this->BoundaryRef(i).Init(num, ref);
      }
        
    // updating element numbers
    for (int i = 0; i < this->GetNbElt(); i++)
      {
        int nb_vert = this->Element(i).GetNbVertices();
        int ref = this->Element(i).GetReference();
        num.Reallocate(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          num(j) = new_num(this->Element(i).numVertex(j));
        
        this->Element(i).Init(num, ref);
      }

    // due to change of numbers we need to sort edges and faces
    this->SortEdgesRef();
    this->SortBoundariesRef();
    
    // duplicate edges are removed
    if (this->GetNbEdgesRef() > 0)
      {
	nb = this->edges_ref.GetM();
        permutation.Reallocate(nb); permutation.Fill();
	Sort(nb, this->edges_ref, permutation);
	
	// assembling edges
	nb = 1; Edge<Dimension3> edge_prec = this->EdgeRef(0);
	Matrix<R_N> OldPoints = this->PointsEdgeRef;
	for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
	  this->PointsEdgeRef(0, k) = OldPoints(permutation(0), k);

	for (int i = 1; i < this->GetNbEdgesRef(); i++)
	  {
	    if (this->EdgeRef(i) != edge_prec)
	      {
		this->EdgeRef(nb) = this->EdgeRef(i);
		for (int k = 0; k < this->PointsEdgeRef.GetN(); k++)
		  this->PointsEdgeRef(nb, k) = OldPoints(permutation(i), k);

		edge_prec = this->EdgeRef(nb);
		nb++;
	      }
	    else
	      {
		// we have found a duplicate, we don't increment nb
		// so that current edge will be removed
	      }
	  }
	
	this->ResizeEdgesRef(nb);
      }

    // duplicate faces are removed
    if (this->GetNbBoundaryRef() > 0)
      {
	nb = this->faces_ref.GetM();
        permutation.Reallocate(nb); permutation.Fill();
	Sort(nb, this->faces_ref, permutation);
	
	// assembling edges
	nb = 1; Face<Dimension3> face_prec = this->BoundaryRef(0);
	Vector<Vector<R_N> > OldPoints = this->PointsFaceRef;
        this->PointsFaceRef(0) = OldPoints(permutation(0));

	for (int i = 1; i < this->GetNbBoundaryRef(); i++)
	  {
	    if (this->BoundaryRef(i) != face_prec)
	      {
		this->BoundaryRef(nb) = this->BoundaryRef(i);
                this->PointsFaceRef(nb) = OldPoints(permutation(i));
                
		face_prec = this->BoundaryRef(nb);
		nb++;
	      }
	    else
	      {
		// we have found a duplicate, we don't increment nb
		// so that current face will be removed
	      }
	  }
	
	this->ResizeBoundariesRef(nb);
      }
    
    this->FindConnectivity();
    this->ProjectPointsOnCurves();
  }
  
  
  //! the coherence of mesh is enforced (elimination of unused faces and vertices)
  template<class T>
  void Mesh<Dimension3, T>::ForceCoherenceMesh()
  {
    int nb_vertices = this->GetNbVertices();
    int nb_elt = this->GetNbElt();
    Vector<bool> VertexUsed(nb_vertices); VertexUsed.Fill(false);
    for (int i = 0; i < nb_elt; i++)
      for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	VertexUsed(this->Element(i).numVertex(j)) = true;
    
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
        
        // modifying vertices
        VectR3 OldVertices = this->Vertices;
	this->ReallocateVertices(nb); nb = 0;
	for (int i = 0; i < OldVertices.GetM(); i++)
	  if (VertexUsed(i))
	    this->Vertices(nb++) = OldVertices(i);

        // counting edges
        int nb_edges_ref = 0;
        for (int i = 0; i < this->GetNbEdgesRef(); i++)
          {
            int n1 = this->EdgeRef(i).numVertex(0);
            int n2 = this->EdgeRef(i).numVertex(0);
            if (VertexUsed(n1) && VertexUsed(n2))
              nb_edges_ref++;
          }
        
        // modifying edges
        Vector<Edge<Dimension3> > old_edge = this->edges_ref;
        Matrix<R3> OldPoints = this->PointsEdgeRef;
        this->ReallocateEdgesRef(nb_edges_ref);
        nb_edges_ref = 0;
        for (int i = 0; i < old_edge.GetM(); i++)
          {
            int n1 = NewNumber(old_edge(i).numVertex(0));
            int n2 = NewNumber(old_edge(i).numVertex(1));
            int ref = old_edge(i).GetReference();
            if ((n1 >= 0) && (n2 >= 0))
              {
                this->EdgeRef(nb_edges_ref).Init(n1, n2, ref);
                for (int k = 0; k < OldPoints.GetN(); k++)
                  this->PointsEdgeRef(nb_edges_ref, k) = OldPoints(i, k);
                
                nb_edges_ref++;
              }
          }
        
        // counting faces
        int nb_faces_ref = 0;
        for (int i = 0; i < this->GetNbBoundaryRef(); i++)
          {
            bool face_used = true;
            for (int j = 0; j < this->BoundaryRef(i).GetNbVertices(); j++)
              if (!VertexUsed(this->BoundaryRef(i).numVertex(j)))
                face_used = false;
            
            if (face_used)
              nb_faces_ref++;
          }
        
        // modifying faces
        Vector<Face<Dimension3> > old_face = this->faces_ref;
        Vector<Vector<R_N> > OldPointsF = this->PointsFaceRef;
        this->ReallocateBoundariesRef(nb_faces_ref);
        nb_faces_ref = 0; IVect num;
        for (int i = 0; i < old_face.GetM(); i++)
          {
            int nb_vert = old_face(i).GetNbVertices();
            bool face_used = true;
            num.Reallocate(nb_vert);
            for (int j = 0; j < nb_vert; j++)
              {
                num(j) = NewNumber(old_face(i).numVertex(j));
                if (num(j) < 0)
                  face_used = false;
              }
            
            if (face_used)
              {
                int ref = old_face(i).GetReference();
                this->BoundaryRef(nb_faces_ref).Init(num, ref);
                this->PointsFaceRef(nb_faces_ref) = OldPointsF(i);
                nb_faces_ref++;
              }
          }
        
        // and modification of elements
        for (int i = 0; i < nb_elt; i++)
          {
            int nb_vert = this->Element(i).GetNbVertices();
            num.Reallocate(nb_vert);
            for (int j = 0; j < nb_vert; j++)
              num(j) = NewNumber(this->Element(i).numVertex(j));
            
            this->Element(i).Init(num, this->Element(i).GetReference());
          }
      }
    
    FindConnectivity();
    this->ProjectPointsOnCurves();
  }

  
  //! Removes referenced boundaries with zero reference
  template<class T>
  void Mesh<Dimension3, T>::ClearZeroBoundaryRef()
  {
    // counting faces
    int nb_faces_ref = 0;
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      if (this->BoundaryRef(i).GetReference() != 0)
        nb_faces_ref++;
    
    if (nb_faces_ref == this->GetNbBoundaryRef())
      return;
        
    // modifying faces
    Vector<Face<Dimension3> > old_face = this->faces_ref;
    Vector<Vector<R_N> > OldPointsF = this->PointsFaceRef;
    this->ReallocateBoundariesRef(nb_faces_ref);
    nb_faces_ref = 0;
    for (int i = 0; i < old_face.GetM(); i++)
      if (old_face(i).GetReference() != 0)
        {
          this->faces_ref(nb_faces_ref) = old_face(i);
          this->PointsFaceRef(nb_faces_ref) = OldPointsF(i);
          nb_faces_ref++;
        }
    
    FindConnectivity();
    this->ProjectPointsOnCurves();
  }
  
  
  //! splits tetrahedra into hexahedra
  /*!
    Each tetrahedron is split in four hexahedra
    \warning No hybrid mesh (with wedges, pyramids and hex) is allowed
  */
  template<class T>
  void Mesh<Dimension3, T>::SplitIntoHexahedra()
  {
    if (!IsOnlyTetrahedral())
      {
	if (this->print_level >= 0)
	  cout<<"The mesh must only comprise tetrahedra"<<endl;
	
	abort();
      }
    
    if (this->print_level >= 2)
      cout<< " we split each tetrahedron into four hexahedra"<<endl;
    
    // we place a new point on each middle of edge, face and tetrahedral
    int nb_vertices = this->GetNbVertices();
    int nb_edges = this->GetNbEdges();
    int nb_edges_ref = this->GetNbEdgesRef();
    int nb_faces = this->GetNbFaces();
    int nb_faces_ref = this->GetNbBoundaryRef();
    int nb_elt = this->GetNbElt();
    int nb_new_vertices = nb_vertices + nb_edges + nb_faces  + nb_elt;
    this->ResizeVertices(nb_new_vertices);
    
    R3 ptA; VectReal_wp middle_edge(1); VectR3 points_interm(1);
    middle_edge.Fill(0.5);
    // we take middle of each edge of the mesh
    for (int i = 0; i < nb_edges; i++)
      {
        if (i >= this->GetNbEdgesRef())
          {
            ptA.Fill(0);
            Add(T(0.5), this->Vertex(this->edges(i).numVertex(0)), ptA);
            Add(T(0.5), this->Vertex(this->edges(i).numVertex(1)), ptA);
          }
        else
          {
            int ref = this->EdgeRef(i).GetReference();
            this->GetPointsOnCurvedEdge(i, ref, middle_edge, points_interm);
            ptA = points_interm(0);
          }
        
	this->Vertex(nb_vertices + i) = ptA;
      }
    
    // we take middle of each face of the mesh
    T one_third = T(1)/3;
    Vector<R2> middle_tri(1); middle_tri(0).Init(one_third, one_third);
    for (int i = 0; i < nb_faces; i++)
      {
        if (i >= nb_faces_ref)
          {
            ptA.Fill(0);
            Add(one_third, this->Vertex(this->Boundary(i).numVertex(0)), ptA);
            Add(one_third, this->Vertex(this->Boundary(i).numVertex(1)), ptA);
            Add(one_third, this->Vertex(this->Boundary(i).numVertex(2)), ptA);
          }
        else
          {
            int ref = this->BoundaryRef(i).GetReference();
            this->GetPointsOnCurvedFace(i, ref, middle_tri, points_interm);
            ptA = points_interm(0);
          }
        
	// middle of the face
	this->Vertex(nb_vertices + nb_edges + i) = ptA;
      }
    
    // we compute each center of tetrahedron
    R3 center_tetra;
    for (int i = 0; i < nb_elt; i++)
      {
	T weight_edge = T(-1)/T(8);
	T weight_face = T(27)/T(64);
	T weight_vertex = T(1)/T(64);
	
	center_tetra.Fill(0);
	
	for (int j = 0; j < 4; j++)
	  Add(weight_vertex, this->Vertex(this->elements(i).numVertex(j)), center_tetra);
	
	for (int j = 0; j < 6; j++)
	  Add(weight_edge, this->Vertex(nb_vertices+this->elements(i).numEdge(j)), center_tetra);
	
	for (int j = 0; j < 4; j++)
	  Add(weight_face, this->Vertex(nb_vertices+nb_edges+this->elements(i).numFace(j)),
              center_tetra);
	
	this->Vertex(nb_vertices+nb_edges+nb_faces+i) = center_tetra;
      }
    
    // sub-points needed for curved approximation
    int order = this->GetGeometryOrder();
    VectReal_wp xi(order-1), XiEdge; VectR2 xi_quad, XiQuad;
    int Nquad = (order-1)*(order-1);
    Matrix<Real_wp> ValPhi, ValPhiQuad;
    if (order > 1)
      {
        for (int k = 0; k < order-1; k++)
          xi(k) = this->lob_curve.Points(k+1);
        
        xi_quad.Reallocate(Nquad);
        for (int k = 0; k < Nquad; k++)
          xi_quad(k) = this->quadrangle_reference.PointsNodalND(4*order+k);
        
        // position of internal points when dividing an edge into two parts
        XiEdge.Reallocate(2*(order-1));
        for (int k = 0; k < order-1; k++)
          {
            XiEdge(k) = 0.5*xi(k);
            XiEdge(k+order-1) = 1.0 - 0.5*xi(k);
          }
        
        XiQuad.Reallocate(3*(order-1)+3*Nquad);
        for (int k = 0; k < order-1; k++)
          {
            // edge (1/2, 0) -> (1/3, 1/3)
            XiQuad(k).Init((1.0-xi(k))*0.5 + xi(k)*one_third, xi(k)*one_third);
            // edge (0, 1/2) -> (1/3, 1/3)
            XiQuad(k+order-1).Init(xi(k)*one_third, (1.0-xi(k))*0.5 + xi(k)*one_third);
            // edge (1/2, 1/2) -> (1/3, 1/3)
            XiQuad(k+2*(order-1)).Init((1.0-xi(k))*0.5 + xi(k)*one_third,
                                       (1.0-xi(k))*0.5 + xi(k)*one_third);
          }
        
        // then inside of the three quads
        int offset = 3*(order-1);
        for (int k = 0; k < Nquad; k++)
          {
            T coef1 = (1.0-xi_quad(k)(0))*(1.0-xi_quad(k)(1));
            T coef2 = xi_quad(k)(0)*(1.0-xi_quad(k)(1));
            T coef3 = xi_quad(k)(0)*xi_quad(k)(1);
            T coef4 = (1.0-xi_quad(k)(0))*xi_quad(k)(1);
            
            // first quad : (0, 0) -> (1/2, 0) -> (1/3, 1/3) -> (0, 1/2)
            XiQuad(offset+k).Init(0.5*coef2 + one_third*coef3, one_third*coef3+0.5*coef4);
            
            // second quad : (1, 0) -> (1/2, 1/2) -> (1/3, 1/3) -> (1/2, 0)
            XiQuad(offset+Nquad+k).Init(coef1 + 0.5*coef2 + one_third*coef3 + 0.5*coef4,
                                        0.5*coef2 + one_third*coef3);
            
            // third quad : (0, 1) -> (0, 1/2) -> (1/3, 1/3) -> (1/2, 1/2)
            XiQuad(offset+2*Nquad+k).Init(one_third*coef3 + 0.5*coef4,
                                          coef1 + 0.5*coef2 + one_third*coef3 + 0.5*coef4);
            
          }
        
        // computing phi_i(XiEdge) and phi^2D_i(XiQuad)
        ValPhi.Reallocate(order+1, 2*(order-1));
        for (int j = 0; j < XiEdge.GetM(); j++)
          for (int i = 0; i <= order; i++)
            ValPhi(i, j) = this->lob_curve.EvaluatePhi(i, XiEdge(j));
        
        ValPhiQuad.Reallocate((order+1)*(order+2)/2, XiQuad.GetM());
        VectReal_wp phi;
        for (int j = 0; j < XiQuad.GetM(); j++)
          {
            this->triangle_reference.ComputeValuesPhiNodalRef(XiQuad(j), phi);
            for (int i = 0; i < ValPhiQuad.GetM(); i++)
              ValPhiQuad(i, j) = phi(i);
          }
        
      }

    // creating new referenced edges
    nb_edges_ref = this->GetNbEdgesRef();
    int nb_new_edges_ref = 2*nb_edges_ref + 3*this->GetNbBoundaryRef();
    Matrix<R_N> OldPoints = this->PointsEdgeRef;
    this->ResizeEdgesRef(nb_new_edges_ref);
    nb_new_edges_ref = nb_edges_ref;
    VectR3 PtEdge(order+1);
    for (int i = 0; i < nb_edges_ref; i++)
      {
        int n1 = this->EdgeRef(i).numVertex(0);
        int n2 = this->EdgeRef(i).numVertex(1);
        int ref = this->EdgeRef(i).GetReference();
        
        // two edges
        this->EdgeRef(i).Init(n1, nb_vertices+i, ref);
        this->EdgeRef(nb_new_edges_ref).Init(n2, nb_vertices + i, ref);
        
        // for curved boundaries computation of internal points
        if (order > 1)
          {
            PtEdge(0) = this->Vertex(n1);
            PtEdge(order) = this->Vertex(n2);
            for (int k = 1; k < order; k++)
              PtEdge(k) = this->PointsEdgeRef(i, k-1);
            
            for (int k = 0; k < order - 1; k++)
              {
                ptA.Zero();
                for (int j = 0; j < ValPhi.GetM(); j++)
                  Add(ValPhi(j, k), PtEdge(j), ptA);
                
                this->PointsEdgeRef(i, k) = ptA;
              }

            for (int k = 0; k < order - 1; k++)
              {
                ptA.Zero();
                for (int j = 0; j < ValPhi.GetM(); j++)
                  Add(ValPhi(j, k+order-1), PtEdge(j), ptA);
                
                this->PointsEdgeRef(nb_new_edges_ref, k) = ptA;
              }
          } 
            
        nb_new_edges_ref++;
      }
    

    points_interm.Reallocate(ValPhiQuad.GetN()); 
    
    // three times more quadrilaterals than triangles
    nb_faces_ref = this->GetNbBoundaryRef();
    int nb_new_faces_ref = 3*nb_faces_ref;
    this->ResizeBoundariesRef(nb_new_faces_ref);
    
    nb_new_faces_ref = nb_faces_ref;
    VectR3 PtFace((order+1)*(order+2)/2);
    for (int i = 0; i < nb_faces_ref; i++)
      {
	int ref = this->BoundaryRef(i).GetReference();
        
        int n1 = this->BoundaryRef(i).numVertex(0);
        int n2 = this->BoundaryRef(i).numVertex(1);
        int n3 = this->BoundaryRef(i).numVertex(2);
        
        int ne1 = nb_vertices + this->BoundaryRef(i).numEdge(0);
        int ne2 = nb_vertices + this->BoundaryRef(i).numEdge(1);
        int ne3 = nb_vertices + this->BoundaryRef(i).numEdge(2);
        
        int nc = nb_vertices + nb_edges + i;
        
        // retrieving points of the face
        if (order > 1)
          {
            PtFace(0) = this->Vertex(n1);
            PtFace(1) = this->Vertex(n2);
            PtFace(2) = this->Vertex(n3);
            for (int j = 0; j < 3; j++)
              {
                int ne = this->BoundaryRef(i).numEdge(j);
                if (ne >= nb_edges_ref)
                  {
                    for (int k = 0; k < order-1; k++)
                      {
                        ptA = (1.0-xi(k))*PtFace(j) + xi(k)*PtFace((j+1)%3);
                        PtFace(3+k+j*(order-1)) = ptA;
                      }
                  }
                else
                  {
                    if (this->BoundaryRef(i).GetOrientationEdge(j))
                      for (int k = 0; k < order-1; k++)
                        PtFace(3+k+j*(order-1)) = OldPoints(ne, k);
                    else
                      for (int k = 0; k < order-1; k++)
                        PtFace(3+k+j*(order-1)) = OldPoints(ne, order-2-k);
                  }
              }
            
            for (int k = 0; k < this->PointsFaceRef(i).GetM(); k++)
              PtFace(3*order + k) = this->PointsFaceRef(i)(k);
            
            // computing internal points
            for (int k = 0; k < ValPhiQuad.GetN(); k++)
              {
                points_interm(k).Zero();
                for (int j = 0; j < ValPhiQuad.GetM(); j++)
                  Add(ValPhiQuad(j, k), PtFace(j), points_interm(k));
              }
          }        

        // creating intern edges
        this->EdgeRef(nb_new_edges_ref).Init(ne1, nc, ref);
        for (int k = 0; k < order-1; k++)
          this->PointsEdgeRef(nb_new_edges_ref, k) = points_interm(k);
        
        nb_new_edges_ref++;
        this->EdgeRef(nb_new_edges_ref).Init(ne3, nc, ref);
        for (int k = 0; k < order-1; k++)
          this->PointsEdgeRef(nb_new_edges_ref, k) = points_interm(k+order-1);
        
        nb_new_edges_ref++;
        this->EdgeRef(nb_new_edges_ref).Init(ne2, nc, ref);
        for (int k = 0; k < order-1; k++)
          this->PointsEdgeRef(nb_new_edges_ref, k) = points_interm(k+2*(order-1));
        
        nb_new_edges_ref++;
        
        // then intern quadrilaterals
        int offset = 3*(order-1);
        this->BoundaryRef(i).InitQuadrangular(n1, ne1, nc, ne3, ref);
        if (order > 1)
          {
            this->PointsFaceRef(i).Reallocate(Nquad);
            for (int k = 0; k < Nquad; k++)
              this->PointsFaceRef(i)(k) = points_interm(offset+k);
          }
        
        this->BoundaryRef(nb_new_faces_ref).InitQuadrangular(n2, ne2, nc, ne1, ref);
        if (order > 1)
          {
            this->PointsFaceRef(nb_new_faces_ref).Reallocate(Nquad);
            for (int k = 0; k < Nquad; k++)
              this->PointsFaceRef(nb_new_faces_ref)(k) = points_interm(offset+k+Nquad);
          }
        
        nb_new_faces_ref++;
        this->BoundaryRef(nb_new_faces_ref).InitQuadrangular(n3, ne3, nc, ne2, ref);
        if (order > 1)
          {
            this->PointsFaceRef(nb_new_faces_ref).Reallocate(Nquad);
            for (int k = 0; k < Nquad; k++)
              this->PointsFaceRef(nb_new_faces_ref)(k) = points_interm(offset+k+2*Nquad);
          }
        
        nb_new_faces_ref++;
      }

    // okay we have now four times more hexahedra 
    int nb_new_elt = 4*nb_elt;
    this->ResizeElements(nb_new_elt);
    nb_new_elt = nb_elt;
    
    IVect num_hexas(8);    
    for (int i = 0; i  < nb_elt; i++)
      {
	// First Hexahedral
	// vertex 0 -> middle edge 0-3 -> middle edge 0-2 -> middle face 0-2-3
	// middle edge 0-1 -> middle face 0-1-3 -> middle face 0-1-2 -> middle tetrahedron
	num_hexas(0) = this->elements(i).numVertex(0);
	num_hexas(1) = nb_vertices+this->elements(i).numEdge(2);
	num_hexas(3) = nb_vertices+this->elements(i).numEdge(1);
	num_hexas(2) = nb_vertices+nb_edges+this->elements(i).numFace(2);
	num_hexas(4) = nb_vertices+this->elements(i).numEdge(0);
	num_hexas(5) = nb_vertices+nb_edges+this->elements(i).numFace(1);
	num_hexas(7) = nb_vertices+nb_edges+this->elements(i).numFace(0);
	num_hexas(6) = nb_vertices+nb_edges+nb_faces+i;
	this->Element(nb_new_elt).InitHexahedral(num_hexas, this->elements(i).GetReference());
	nb_new_elt++;
	
	// Second Hexahedral
	// middle edge 0-1 -> middle face 0-1-3 -> middle face 0-1-2 -> middle tetrahedron
	// vertex 1 -> middle edge 1-3 -> middle edge 1-2 -> middle face 1-2-3
	num_hexas(0) = nb_vertices+this->elements(i).numEdge(0);
	num_hexas(1) = nb_vertices+nb_edges+this->elements(i).numFace(1);
	num_hexas(3) = nb_vertices+nb_edges+this->elements(i).numFace(0);
	num_hexas(2) = nb_vertices+nb_edges+nb_faces+i;
	num_hexas(4) = this->elements(i).numVertex(1);
	num_hexas(5) = nb_vertices+this->elements(i).numEdge(4);
	num_hexas(7) = nb_vertices+this->elements(i).numEdge(3);
	num_hexas(6) = nb_vertices+nb_edges+this->elements(i).numFace(3);
	this->Element(nb_new_elt).InitHexahedral(num_hexas, this->elements(i).GetReference());
	nb_new_elt++;
	
	// Third Hexahedral
	// middle edge 0-2 -> middle face 0-2-3 -> vertex 2 -> middle edge 2-3
	// middle face 0-1-2 -> middle tetrahedron -> middle edge 1-2 -> middle face 1-2-3
	num_hexas(0) = nb_vertices+this->elements(i).numEdge(1);
	num_hexas(1) = nb_vertices+nb_edges+this->elements(i).numFace(2);
	num_hexas(3) = this->elements(i).numVertex(2);
	num_hexas(2) = nb_vertices+this->elements(i).numEdge(5);
	num_hexas(4) = nb_vertices+nb_edges+this->elements(i).numFace(0);
	num_hexas(5) = nb_vertices+nb_edges+nb_faces+i;
	num_hexas(7) = nb_vertices+this->elements(i).numEdge(3);
	num_hexas(6) = nb_vertices+nb_edges+this->elements(i).numFace(3);
	this->Element(nb_new_elt).InitHexahedral(num_hexas,this->elements(i).GetReference());
	nb_new_elt++;
	
	// Fourth Hexahedral
	// middle edge 0-3 -> Vertex 3 -> middle face 0-2-3 -> middle edge 2-3
	// middle face 0-1-3 -> middle edge 1-3 -> middle tetrahedron -> middle face 1-2-3
	num_hexas(0) = nb_vertices+this->elements(i).numEdge(2);
	num_hexas(1) = this->elements(i).numVertex(3);
	num_hexas(3) = nb_vertices+nb_edges+this->elements(i).numFace(2);
	num_hexas(2) = nb_vertices+this->elements(i).numEdge(5);
	num_hexas(4) = nb_vertices+nb_edges+this->elements(i).numFace(1);
	num_hexas(5) = nb_vertices+this->elements(i).numEdge(4);
	num_hexas(7) = nb_vertices+nb_edges+nb_faces+i;
	num_hexas(6) = nb_vertices+nb_edges+this->elements(i).numFace(3);
	this->Element(i).InitHexahedral(num_hexas,this->elements(i).GetReference());
      }
    
    ReorientElements();
    FindConnectivity();
    this->ProjectPointsOnCurves();
  }
  
  
  //! splits an hybrid mesh into tetrahedra
  /*!
    Quadrilateral faces are split into two triangles
    Then usually, each hexahedron is split into 6 tetrahedra,
    each prism into 3 tets, each pyramid into 2 tets
    However, some pathologic cases may appear, and in that case a point is created at the center of
    hex or prism, and the hexahedron is split into 12 test, and the prism split into 8 tets
  */
  template<class T>
  void Mesh<Dimension3, T>::SplitIntoTetrahedra()
  {
    if (this->print_level >= 2)
      cout<< " we split each hexahedron in six tetrahedra,"
          <<" each pyramid in two tets, each wedge in 3 tets"<<endl;
    
    int nb_elt = this->GetNbElt();
    
    // sub-points needed for curved approximation
    int order = this->GetGeometryOrder();
    VectReal_wp xi; VectR2 xi_tri, XiTri;
    int Ntri = (order-1)*(order-2)/2;
    Matrix<Real_wp> ValPhi;
    R3 ptA;
    if (order > 1)
      {
        xi.Reallocate(order-1);
        for (int k = 0; k < order-1; k++)
          xi(k) = this->lob_curve.Points(k+1);
        
        xi_tri.Reallocate(Ntri);
        for (int k = 0; k < Ntri; k++)
          xi_tri(k) = this->triangle_reference.PointsNodalND(3*order+k);
        
        XiTri.Reallocate(order-1+2*Ntri);
        // edge (1, 0) -> (0, 1)
        for (int k = 0; k < order-1; k++)
          XiTri(k).Init((1.0-xi(k)), xi(k));
        
        // then inside of the two triangles
        for (int k = 0; k < Ntri; k++)
          {
            T coef1 = 1.0 - xi_tri(k)(0) - xi_tri(k)(1);
            T coef2 = xi_tri(k)(0);
            T coef3 = xi_tri(k)(1);
            
            // first triangle : (0, 0) -> (1, 0) -> (0, 1)
            XiTri(order-1+k).Init(coef2, coef3);
            
            // second triangle : (1, 0) -> (1, 1) -> (0, 1)
            XiTri(order-1+Ntri+k).Init(coef1+coef2, coef2+coef3);
          }
        
        // computing phi^2D_i(XiTri)
        ValPhi.Reallocate((order+1)*(order+1), XiTri.GetM());
        VectReal_wp phi;
        for (int j = 0; j < XiTri.GetM(); j++)
          {
            this->quadrangle_reference.ComputeValuesPhiNodalRef(XiTri(j), phi);
            for (int i = 0; i < ValPhi.GetM(); i++)
              ValPhi(i, j) = phi(j);
          }
      }

    // creating new referenced edges
    int nb_edges_ref = this->GetNbEdgesRef();
    int nb_new_edges_ref = nb_edges_ref + this->GetNbQuadranglesRef();
    this->ResizeEdgesRef(nb_new_edges_ref);
    nb_new_edges_ref = nb_edges_ref;

    VectR3 points_interm(ValPhi.GetN());     
    // splitting each quadrilateral face into two triangles
    int nb_faces_ref = this->GetNbBoundaryRef();
    int nb_new_faces_ref = nb_faces_ref + this->GetNbQuadranglesRef();
    this->ResizeBoundariesRef(nb_new_faces_ref);
    
    nb_new_faces_ref = nb_faces_ref;
    VectR3 PtFace((order+1)*(order+1));
    for (int i = 0; i < nb_faces_ref; i++)
      if (this->BoundaryRef(i).GetNbVertices() == 4)
        {
          int ref = this->BoundaryRef(i).GetReference();
          
          int n1 = this->BoundaryRef(i).numVertex(0);
          int n2 = this->BoundaryRef(i).numVertex(1);
          int n3 = this->BoundaryRef(i).numVertex(2);
          int n4 = this->BoundaryRef(i).numVertex(3);
          
          // retrieving points of the face
          if (order > 1)
            {
              this->GetNodalPointsOnFace(i, PtFace);
              
              // computing internal points
              for (int k = 0; k < ValPhi.GetN(); k++)
                {
                  points_interm(k).Zero();
                  for (int j = 0; j < ValPhi.GetM(); j++)
                    Add(ValPhi(j, k), PtFace(j), points_interm(k));
                }
            }        
          
          // creating intern edge
          if (n2 < n4)
            {
              this->EdgeRef(nb_new_edges_ref).Init(n2, n4, ref);
              for (int k = 0; k < order-1; k++)
                this->PointsEdgeRef(nb_new_edges_ref, k) = points_interm(k);
            }
          else
            {
              this->EdgeRef(nb_new_edges_ref).Init(n4, n2, ref);
              for (int k = 0; k < order-1; k++)
                this->PointsEdgeRef(nb_new_edges_ref, k) = points_interm(order-2-k);
            }
          
          nb_new_edges_ref++;
          
          // then triangles
          this->BoundaryRef(i).InitTriangular(n1, n2, n4, ref);
          if (order > 1)
            {
              this->PointsFaceRef(i).Reallocate(Ntri);
              for (int k = 0; k < Ntri; k++)
                this->PointsFaceRef(i)(k) = points_interm(order-1+k);
            }
        
          this->BoundaryRef(nb_new_faces_ref).InitTriangular(n2, n3, n4, ref);
          if (order > 1)
            {
              this->PointsFaceRef(nb_new_faces_ref).Reallocate(Ntri);
              for (int k = 0; k < Ntri; k++)
                this->PointsFaceRef(nb_new_faces_ref)(k) = points_interm(order-1+k+Ntri);
            }
        
          nb_new_faces_ref++;
        }
    
    // we may insert new vertices in order to split each hex in tets
    int nb_vertices = this->GetNbVertices();
    int nb_new_vertices = nb_vertices;
    this->ResizeVertices(nb_vertices + nb_elt);
    // and for elements, worst case will be 12 tet per hex, 8 per wedge
    int nb_new_elt = this->GetNbTetrahedra() + 2*this->GetNbPyramids()
      + 8*this->GetNbWedges() + 12*this->GetNbHexahedra();
    
    Vector<Volume> old_elt = this->elements;
    this->ReallocateElements(nb_new_elt);
    nb_new_elt = 0;
    for (int i = 0; i < nb_elt; i++)
      if (old_elt(i).GetNbVertices() > 4)
        {
          // only non-tetrahedral elements are concerned
          int nb_vert = old_elt(i).GetNbVertices();
          int type = old_elt(i).GetHybridType();
          int ref = old_elt(i).GetReference();
          IVect num(nb_vert); R3 center;
          for (int j = 0; j < nb_vert; j++)
            {
              num(j) = old_elt(i).numVertex(j);
              center += this->Vertex(num(j));
            }
          
          center *= T(1)/nb_vert;
          
          // for each quadrilateral face, we check if the diagonal is the same
          // for the element and the global face, if not wf(j) = false
          int nb_faces = old_elt(i).GetNbFaces();
          Vector<bool> wf(nb_faces); wf.Fill(true);
          IVect nf(nb_faces); nf.Fill(-1);
          for (int j = 0; j < nb_faces; j++)
            if (MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type, j))
              {
                nf(j) = old_elt(i).numFace(j);
                int ne = old_elt(i).numVertex(MeshNumbering<Dimension3>::
                                              FirstExtremityFace(type, j));
                if ((ne == this->Boundary(nf(j)).numVertex(0))
                    ||(ne == this->Boundary(nf(j)).numVertex(2)))
                  wf(j) = true;
                else
                  wf(j) = false;
              }
          
          // different treatment for pyramids, wedges and hexa
          switch (nb_vert)
            {
            case 5 :
              {
                // pyramid is split into two tets
                if (wf(0))
                  {
                    this->Element(nb_new_elt++)
                      .InitTetrahedral(num(0), num(1), num(3), num(4), ref);
                    this->Element(nb_new_elt++)
                      .InitTetrahedral(num(1), num(2), num(3), num(4), ref);
                  }
                else
                  {
                    this->Element(nb_new_elt++)
                      .InitTetrahedral(num(0), num(1), num(2), num(4), ref);
                    this->Element(nb_new_elt++)
                      .InitTetrahedral(num(2), num(3), num(0), num(4), ref);
                  }
              }
              break;
            case 6 :
              {
                // wedge is split into three tets (or six tets for pathologic configurations)
                if (wf(1))
                  {
                    if (wf(2))
                      {
                        if (wf(3))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(3), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(4), num(5), ref);
                          }
                        else
                          {
                            int nc = nb_new_vertices; this->Vertex(nb_new_vertices++) = center;
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(4), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(4), num(5), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(5), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(5), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), nc, ref);
                          }
                      }
                    else
                      {
                        if (wf(3))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(3), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                          }
                        else
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(5), ref);
                          }
                      }
                  }
                else
                  {
                    if (wf(2))
                      {
                        if (wf(3))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(2), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(4), num(5), ref);
                          }
                        else
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(4), num(5), ref);
                          }
                      }
                    else
                      {
                        if (wf(3))
                          {
                            int nc = nb_new_vertices; this->Vertex(nb_new_vertices++) = center;
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(3), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(5), num(4), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(5), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(0), num(3), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), nc, ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), nc, ref);
                          }
                        else
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(5), ref); 
                          }
                      }
                  }
              }
              break;
            case 8 :
              {
                // if wf(3) false, we permute number to have wf(3) always true
                if (!wf(3))
                  {
                    IVect old_num = num;
                    Vector<bool> old_wf = wf;
                    num(0) = old_num(1); num(1) = old_num(2);
                    num(2) = old_num(3); num(3) = old_num(0);
                    num(4) = old_num(5); num(5) = old_num(6);
                    num(6) = old_num(7); num(7) = old_num(4);
                    
                    wf(0) = !old_wf(1); wf(1) = old_wf(5); wf(2) = !old_wf(2);
                    wf(3) = !old_wf(3); wf(4) = old_wf(0); wf(5) = !old_wf(4);
                  }
                
                // hexahedron
                bool pathologic = false;
                if (wf(2))
                  {
                    if (wf(0) && wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(7), ref);
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                      }
                    else if (!wf(0)&&wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(7), ref);
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            pathologic = true;
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                      }
                    else if (wf(0)&&!wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            pathologic = true;
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                      }
                    else if (!wf(0) && !wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(5), num(6), num(7), ref);
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(3), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                          }
                      }
                  }
                else
                  {
                    if (wf(0) && wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(4), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            pathologic = true;
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(4), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(5), ref);
                          }
                      }
                    else if (!wf(0)&&wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(7), ref);
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            pathologic = true;
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(3), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(5), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(2), num(3), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(5), ref);
                          }
                      }
                    else if (wf(0)&&!wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            pathologic = true;
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            pathologic = true;
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            pathologic = true;
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            pathologic = true;
                          }
                      }
                    else if (!wf(0) && !wf(5))
                      {
                        if (wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(4), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(7), ref);
                          }
                        else if (!wf(1)&&wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(1), num(2), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(7), ref);
                          }
                        else if (wf(1)&&!wf(4))
                          {
                            pathologic = true;
                          }
                        else if (!wf(1)&&!wf(4))
                          {
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(5), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(5), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(4), num(5), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(3), num(6), num(7), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(1), num(2), num(6), ref);
                            this->Element(nb_new_elt++)
                              .InitTetrahedral(num(0), num(2), num(3), num(6), ref);
                          }
                      }
                  }
                
                if (pathologic)
                  {
                    int nc = nb_new_vertices; this->Vertex(nb_new_vertices++) = center;
                    if (wf(0))
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(3), num(4), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(3), num(4), num(7), nc, ref);
                      }
                    else
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(3), num(7), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(4), num(7), nc, ref);
                      }
                    
                    if (wf(1))
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(1), num(4), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(1), num(4), num(5), nc, ref);
                      }
                    else
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(1), num(5), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(4), num(5), nc, ref);
                      }

                    if (wf(2))
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(1), num(3), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(1), num(2), num(3), nc, ref);
                      }
                    else
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(1), num(2), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(0), num(2), num(3), nc, ref);
                      }

                    if (wf(3))
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(4), num(5), num(7), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(5), num(6), num(7), nc, ref);
                      }
                    else
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(4), num(5), num(6), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(4), num(6), num(7), nc, ref);
                      }

                    if (wf(4))
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(2), num(3), num(7), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(2), num(6), num(7), nc, ref);
                      }
                    else
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(2), num(3), num(6), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(3), num(6), num(7), nc, ref);
                      }

                    if (wf(5))
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(1), num(2), num(5), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(2), num(5), num(6), nc, ref);
                      }
                    else
                      {
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(1), num(2), num(6), nc, ref);
                        this->Element(nb_new_elt++)
                          .InitTetrahedral(num(1), num(5), num(6), nc, ref);
                      }
                    
                  }
                
              }
              break;              
            }
          
        }  // end loop over elements
    
    this->ResizeElements(nb_new_elt);
    this->ResizeVertices(nb_new_vertices);
    
    ReorientElements();
    FindConnectivity();
    this->ProjectPointsOnCurves();
  }
  
  
  //! the mesh is split with subdivisions specified by the user
  /*!
    \param[in] step_subdiv subdivisions (by example 0 0.5 1)
  */
  template<class T>
  void Mesh<Dimension3, T>::SubdivideMesh(const Vector<T>& step_subdiv)
  {
    Vector<IVect> Nodle; Vector<VectR3> points_div;
    Vector<IVect> NodleSurf; Vector<VectR2> points_surf;
    SubdivideMesh(step_subdiv, points_surf, points_div, Nodle, NodleSurf);
  }


  //! the mesh is split into sub-hexahedra/sub-tetrahedra
  /*!
    \param[in] outside_subdiv intern points, including extremities 
    (ie 0  0.2 0.5 0.6 0.8  1.0, for a regular subdivision in 5 intervals)
    \param[out] Nodle mesh numbering of the split mesh
    \param[out] points_div subdivision points on unit tetrahedron, cube, etc
    \param[out] points_surf subdivision points on unit triangle, quadrangle
    \param[out] NodleSurf numbering for referenced faces of the mesh
  */
  template<class T>
  void Mesh<Dimension3, T>::
  SubdivideMesh(const Vector<T>& outside_subdiv,
                Vector<VectR2>& points_surf, Vector<VectR3>& points_div,
                Vector<IVect>& Nodle, Vector<IVect>& NodleSurf)
  {    
    if (outside_subdiv.GetM() < 2)
      return;
    else if (outside_subdiv.GetM() == 2)
      {
	NodleSurf.Reallocate(this->GetNbBoundaryRef());
	for (int i = 0; i < this->GetNbBoundaryRef(); i++)
	  {
	    NodleSurf(i).Reallocate(this->BoundaryRef(i).GetNbVertices());
	    for (int j = 0; j < this->BoundaryRef(i).GetNbVertices(); j++)
	      NodleSurf(i)(j) = this->BoundaryRef(i).numVertex(j);
	  }
        
	points_surf.Reallocate(2);
	points_surf(0).Reallocate(3); points_surf(1).Reallocate(4);
	    
	points_surf(0)(0)(0) = 0.0; points_surf(0)(0)(1) = 0.0;
	points_surf(0)(1)(0) = 1.0; points_surf(0)(1)(1) = 0.0;
	points_surf(0)(2)(0) = 0.0; points_surf(0)(2)(1) = 1.0;
	points_surf(1)(0)(0) = 0.0; points_surf(1)(0)(1) = 0.0;
	points_surf(1)(1)(0) = 1.0; points_surf(1)(1)(1) = 0.0;
	points_surf(1)(2)(0) = 1.0; points_surf(1)(2)(1) = 1.0;
	points_surf(1)(3)(0) = 0.0; points_surf(1)(3)(1) = 1.0;
        
	Nodle.Reallocate(this->GetNbElt());
	for (int i = 0; i < this->GetNbElt(); i++)
	  {
	    Nodle(i).Reallocate(this->Element(i).GetNbVertices());
	    for (int j = 0; j < this->Element(i).GetNbVertices(); j++)
	      Nodle(i)(j) = this->Element(i).numVertex(j);
	  }
	
	points_div.Reallocate(4);
	// points of the unit tetrahedron
	points_div(0).Reallocate(4);
	points_div(0)(0).Init(0.0, 0.0, 0.0);
	points_div(0)(1).Init(1.0, 0.0, 0.0);
	points_div(0)(2).Init(0.0, 1.0, 0.0);
	points_div(0)(3).Init(0.0, 0.0, 1.0);
	
	// points of the unit pyramid
	points_div(1).Reallocate(5);
	points_div(1)(0).Init(-1.0, -1.0, 0.0);
	points_div(1)(1).Init(1.0, -1.0, 0.0);
	points_div(1)(2).Init(1.0, 1.0, 0.0);
	points_div(1)(3).Init(-1.0, 1.0, 1.0);
	points_div(1)(4).Init(0.0, 0.0, 1.0);
	
	// points of the unit wedge
	points_div(2).Reallocate(6);
	points_div(2)(0).Init(0.0, 0.0, 0.0);
	points_div(2)(1).Init(1.0, 0.0, 0.0);
	points_div(2)(2).Init(0.0, 1.0, 0.0);
	points_div(2)(3).Init(0.0, 0.0, 1.0);
	points_div(2)(4).Init(1.0, 0.0, 1.0);
	points_div(2)(5).Init(0.0, 1.0, 1.0);
	
	// points of the unit cube
	points_div(3).Reallocate(8);
	points_div(3)(0).Init(0.0, 0.0, 0.0);
	points_div(3)(1).Init(1.0, 0.0, 0.0);
	points_div(3)(2).Init(1.0, 1.0, 0.0);
	points_div(3)(3).Init(0.0, 1.0, 0.0);
	points_div(3)(4).Init(0.0, 0.0, 1.0);
	points_div(3)(5).Init(1.0, 0.0, 1.0);
	points_div(3)(6).Init(1.0, 1.0, 1.0);
	points_div(3)(7).Init(0.0, 1.0, 1.0);
	
	return;
      }
    
    Vector<T> step_subdiv(outside_subdiv.GetM()-2);
    for (int i = 0; i < step_subdiv.GetM(); i++)
      step_subdiv(i) = outside_subdiv(i+1);

    if (this->print_level >= 2)
      cout<< " we divide the initial mesh "<<endl;
        
    // step_subdiv is an array containing all the intermediary points to add on a edge
    // exemple step_subdiv={0.5} -> we will split each quad in 4 sub-quad
    // and each triangle in 4 sub-triangles
    int nb_div = step_subdiv.GetM();
    Vector<T> Points1D(nb_div+2);
    points_surf.Reallocate(2);
    Vector<Matrix<int> > NumNodes2D(2);
    Matrix<int> CoordinateNodes;
    
    int r = nb_div+1;
    
    Points1D(0) = T(0);  Points1D(r) = T(1);
    for (int i = 0; i < nb_div; i++)
      Points1D(i+1) = step_subdiv(i);
    
    // we get NumNodes2D for triangle and quadrilateral
    MeshNumbering<Dimension2, T>::
      ConstructTriangularNumbering(nb_div+1, NumNodes2D(0), CoordinateNodes);
    points_surf(0).Reallocate((r+1)*(r+2)/2);
  
    MeshNumbering<Dimension2, T>::
      ConstructQuadrilateralNumbering(nb_div+1, NumNodes2D(1), CoordinateNodes);
    points_surf(1).Reallocate((r+1)*(r+1));
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  points_surf(1)(NumNodes2D(1)(i,j)).Init(Points1D(i),Points1D(j));
	  if (NumNodes2D(0)(i,j) != -1)
	    points_surf(0)(NumNodes2D(0)(i,j)).Init(Points1D(i), Points1D(j));
	}

    points_div.Reallocate(4);
    points_div(0).Reallocate((r+1)*(r+2)*(r+3)/6);
    points_div(1).Reallocate((r+1)*(r+2)*(2*r+3)/6);
    points_div(2).Reallocate((r+1)*(r+2)*(r+1)/2);
    points_div(3).Reallocate((r+1)*(r+1)*(r+1));
    
    // points on tetrahedron
    Vector<Array3D<int> > NumNodes3D(4); Matrix<int> coor;
    MeshNumbering<Dimension3, T>::ConstructTetrahedralNumbering(r, NumNodes3D(0), coor);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= r-i-j; k++)
          points_div(0)(NumNodes3D(0)(i, j, k)).Init(Points1D(i), Points1D(j), Points1D(k));
    
    // points on pyramid
    MeshNumbering<Dimension3, T>::ConstructPyramidalNumbering(r, NumNodes3D(1), coor);
    for (int k = 0; k <= r; k++)
      for (int i = 0; i <= r-k; i++)
        for (int j = 0; j <= r-k; j++)
          {
            Real_wp z = Points1D(k);
            Real_wp x(0), y(0);
            if (k < r)
              {
                x = Points1D(i)/Points1D(r-k);
                y = Points1D(j)/Points1D(r-k);
              }
            
            points_div(1)(NumNodes3D(1)(i, j, k))
              .Init((1.0-z)*(2.0*x-1.0), (1.0-z)*(2.0*y-1.0), z);
          }
    
    // points on prism
    MeshNumbering<Dimension3, T>::ConstructPrismaticNumbering(r, NumNodes3D(2), coor);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= r; k++)
          points_div(2)(NumNodes3D(2)(i, j, k)).Init(Points1D(i), Points1D(j), Points1D(k));
    
    // points on hexahedron
    MeshNumbering<Dimension3, T>::ConstructHexahedralNumbering(r, NumNodes3D(3), coor);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        for (int k = 0; k <= r; k++)
          points_div(3)(NumNodes3D(3)(i, j, k)).Init(Points1D(i), Points1D(j), Points1D(k));
    
    // now we compute all the nodes and fill Nodle
    Vector<R_N> PosNodes;
    // reconstructing edges (in the case where it is not done)
    GetMeshSubdivision(outside_subdiv, points_div, points_surf, PosNodes, Nodle, NodleSurf);
    
    int nb_old_vert = this->GetNbVertices();
    // new vertices of the subdivided mesh
    this->ReallocateVertices(PosNodes.GetM());
    for (int i = 0; i < PosNodes.GetM(); i++)
      this->Vertex(i) = PosNodes(i);
    
    PosNodes.Clear();
    
    int nb_triangles = this->GetNbTrianglesRef();
    int nb_quadrangles = this->GetNbQuadranglesRef();
    int nb_pyramids = this->GetNbPyramids();
    int nb_old_elt = this->GetNbElt();
    // order of approximation for geometry
    int rg = this->lob_curve.GetOrder(); 
    
    // now we subdivide referenced edges
    int nb_old_edges_ref = this->GetNbEdgesRef();
    int nb_new_edges_ref = nb_old_edges_ref*r;
    // new edges due to subdivision of referenced faces
    if (rg > 1)
      {
        nb_new_edges_ref += 3*(r-1)*r*nb_triangles/2;
        nb_new_edges_ref += 2*(r-1)*r*nb_quadrangles;
      }
    Vector<Edge<Dimension3> > new_edge(nb_new_edges_ref);
    
    // first, we compute internal points on referenced edges (curved mesh)
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
	    OldPos(0) = this->Vertex(this->EdgeRef(i).numVertex(0));
	    OldPos(rg) = this->Vertex(this->EdgeRef(i).numVertex(1));
	    for (int k = 1; k < rg; k++)
	      OldPos(k) = this->PointsEdgeRef(i, k-1);
            
	    int ref = this->EdgeRef(i).GetReference();
	    if (this->GetCurveType(ref) != this->CURVE_FILE )
	      {
		// predefined curve, we compute the image of intermediary points
		Vector<R_N> points_interm;
		this->GetPointsOnCurvedEdge(i, ref, Xi, points_interm);
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

    IVect num(r+1); nb_new_edges_ref = 0;
    for (int i = 0; i < nb_old_edges_ref; i++)
      {
	int ref = this->EdgeRef(i).GetReference();	
	num(0) = this->EdgeRef(i).numVertex(0);
	num(r) = this->EdgeRef(i).numVertex(1);
	for (int j = 1; j < r; j++)
	  num(j) = nb_old_vert + i*(r-1) + j-1;
	
	for (int j = 0; j <= nb_div; j++)
	  {
	    if (num(j) < num(j+1))
	      new_edge(nb_new_edges_ref++).Init(num(j), num(j+1), ref);
	    else
	      new_edge(nb_new_edges_ref++).Init(num(j+1), num(j), ref);
            
	  }
      }
    
    // now we subdivide referenced faces
    int nb_old_faces_ref = this->GetNbBoundaryRef();
    int nb_new_faces_ref = nb_old_faces_ref*r*r;
    
    // first, we compute internal points on referenced edges (curved mesh)
    Vector<Vector<R_N> > NewPointsF(nb_new_faces_ref);

    if (rg > 1)
      {
	NewPointsF.Reallocate(nb_new_faces_ref);
        int Ntri = (rg-1)*(rg-2)/2;
        int Nquad = (rg-1)*(rg-1);
        int nb_edges_tri = 3*(r-1)*r/2;
        int nb_edges_quad = 2*(r-1)*r;
        nb_new_faces_ref = 0;
        
        // 1-D points
        VectReal_wp xi(rg-1); R2 pt;
        for (int k = 1; k < rg; k++)
          xi(k-1) = this->lob_curve.Points(k);
        
        VectR2 xi_tri(Ntri), xi_quad(Nquad);
        for (int k = 0; k < Ntri; k++)
          xi_tri(k) = this->triangle_reference.PointsNodalND(3*rg + k);
        
        for (int k = 0; k < Nquad; k++)
          xi_quad(k) = this->quadrangle_reference.PointsNodalND(4*rg + k);
        
        // construction of internal points for triangles
	Vector<R_N> OldPosTri((rg+1)*(rg+2)/2);
        Matrix<Real_wp> ValPhiTri((rg+1)*(rg+2)/2, (rg-1)*nb_edges_tri + Ntri*r*r);
        VectR2 XiTri(ValPhiTri.GetN());
        
        // creation of points on sub-edges of triangles
        int nb = 0;
        for (int i = 0; i < r-1; i++)
          for (int j = 0; j < r-1-i; j++)
            {
              // we consider the triangle
              // (i+1, j) -> (i+1, j+1) -> (i, j+1)
              int n1 = NumNodes2D(0)(i+1, j);
              int n2 = NumNodes2D(0)(i+1, j+1);
              int n3 = NumNodes2D(0)(i, j+1);
              
              // points on the edges
              for (int k = 0; k < rg-1; k++)
                {
                  pt.Zero();
                  Add(1.0-xi(k), points_surf(0)(n1), pt);
                  Add(xi(k), points_surf(0)(n2), pt);
                  XiTri(nb++) = pt;
                }

              for (int k = 0; k < rg-1; k++)
                {
                  pt.Zero();
                  Add(1.0-xi(k), points_surf(0)(n2), pt);
                  Add(xi(k), points_surf(0)(n3), pt);
                  XiTri(nb++) = pt;
                }

              for (int k = 0; k < rg-1; k++)
                {
                  pt.Zero();
                  Add(1.0-xi(k), points_surf(0)(n3), pt);
                  Add(xi(k), points_surf(0)(n1), pt);
                  XiTri(nb++) = pt;
                }
            }
              
        // creating of points inside sub-triangles of triangle
        for (int i = 0; i < r; i++)
          {
            for (int j = 0; j < r-1-i; j++)
              {
                // we consider the triangle
                // (i+1, j) -> (i+1, j+1) -> (i, j+1)
                int n1 = NumNodes2D(0)(i+1, j);
                int n2 = NumNodes2D(0)(i+1, j+1);
                int n3 = NumNodes2D(0)(i, j+1);
                
                for (int k = 0; k < Ntri; k++)
                  {
                    pt.Zero();
                    Add(1.0-xi_tri(k)(0)-xi_tri(k)(1), points_surf(0)(n1), pt);
                    Add(xi_tri(k)(0), points_surf(0)(n2), pt);
                    Add(xi_tri(k)(1), points_surf(0)(n3), pt);
                    XiTri(nb++) = pt;
                  }
                
                // and triangle
                // (i, j) -> (i+1, j) -> (i, j+1)
                
                n1 = NumNodes2D(0)(i, j);
                n2 = NumNodes2D(0)(i+1, j);
                n3 = NumNodes2D(0)(i, j+1);
                
                for (int k = 0; k < Ntri; k++)
                  {
                    pt.Zero();
                    Add(1.0-xi_tri(k)(0)-xi_tri(k)(1), points_surf(0)(n1), pt);
                    Add(xi_tri(k)(0), points_surf(0)(n2), pt);
                    Add(xi_tri(k)(1), points_surf(0)(n3), pt);
                    XiTri(nb++) = pt;
                  }
              }
            
            // last triangle (i, r-1-i) -> (i+1, r-1-i) -> (i, r-i)
            int n1 = NumNodes2D(0)(i, r-1-i);
            int n2 = NumNodes2D(0)(i+1, r-1-i);
            int n3 = NumNodes2D(0)(i, r-i);
            
            for (int k = 0; k < Ntri; k++)
              {
                pt.Zero();
                Add(1.0-xi_tri(k)(0)-xi_tri(k)(1), points_surf(0)(n1), pt);
                Add(xi_tri(k)(0), points_surf(0)(n2), pt);
                Add(xi_tri(k)(1), points_surf(0)(n3), pt);
                XiTri(nb++) = pt;
              }
          }
                
        // then we compute ValPhiTri
        VectReal_wp phi;
        for (int j = 0; j < ValPhiTri.GetN(); j++)
          {
            this->triangle_reference.ComputeValuesPhiNodalRef(XiTri(j), phi);
            for (int i = 0; i < ValPhiTri.GetM(); i++)
              ValPhiTri(i, j) = phi(i);
          }
                
        // same stuff for quadrangle
	Vector<R_N> OldPosQuad((rg+1)*(rg+1));
        Matrix<Real_wp> ValPhiQuad((rg+1)*(rg+1), (rg-1)*nb_edges_quad + Nquad*r*r);
        VectR2 XiQuad(ValPhiQuad.GetN());

        // creation of points on sub-edges of quadrangles
        nb = 0;
        for (int i = 0; i < r-1; i++)
          for (int j = 0; j < r; j++)
            {
              // we consider edge (i+1, j) -> (i+1, j+1)
              int n1 = NumNodes2D(1)(i+1, j);
              int n2 = NumNodes2D(1)(i+1, j+1);
              for (int k = 0; k < rg-1; k++)
                {
                  pt.Zero();
                  Add(1.0-xi(k), points_surf(1)(n1), pt);
                  Add(xi(k), points_surf(1)(n2), pt);
                  XiQuad(nb++) = pt;
                }

              // we consider edge (j, i+1) -> (j+1, i+1)
              n1 = NumNodes2D(1)(j, i+1);
              n2 = NumNodes2D(1)(j+1, i+1);
              for (int k = 0; k < rg-1; k++)
                {
                  pt.Zero();
                  Add(1.0-xi(k), points_surf(1)(n1), pt);
                  Add(xi(k), points_surf(1)(n2), pt);
                  XiQuad(nb++) = pt;
                }
            }
              
        // creating of points inside sub-quadrilaterals of quadrangle
        for (int i = 0; i < r; i++)
          {
            for (int j = 0; j < r; j++)
              {
                // we consider the quadrangle
                // (i, j) -> (i+1, j) -> (i+1, j+1) -> (i, j+1)
                int n1 = NumNodes2D(1)(i, j);
                int n2 = NumNodes2D(1)(i+1, j);
                int n3 = NumNodes2D(1)(i+1, j+1);
                int n4 = NumNodes2D(1)(i, j+1);
                
                for (int k = 0; k < Nquad; k++)
                  {
                    pt.Zero();
                    Add((1.0-xi_quad(k)(0))*(1.0-xi_quad(k)(1)), points_surf(1)(n1), pt);
                    Add(xi_quad(k)(0)*(1.0-xi_quad(k)(1)), points_surf(1)(n2), pt);
                    Add(xi_quad(k)(0)*xi_quad(k)(1), points_surf(1)(n3), pt);
                    Add((1.0-xi_quad(k)(0))*xi_quad(k)(1), points_surf(1)(n4), pt);
                    XiQuad(nb++) = pt;
                  }
              }
          }
        
        // then we compute ValPhiQuad
        for (int j = 0; j < ValPhiQuad.GetN(); j++)
          {
            this->quadrangle_reference.ComputeValuesPhiNodalRef(XiQuad(j), phi);
            for (int i = 0; i < ValPhiQuad.GetM(); i++)
              ValPhiQuad(i, j) = phi(i);
          }
        
        // then loop over all referenced faces
	for (int n = 0; n < nb_old_faces_ref; n++)
	  {
            int ref = this->BoundaryRef(n).GetReference();
            if ( this->BoundaryRef(n).IsTriangular())
              {
                this->GetNodalPointsOnFace(n, OldPosTri);
                                
                Vector<R_N> points_interm(ValPhiTri.GetN());
                if (this->GetCurveType(ref) != this->CURVE_FILE )
                  {
                    // predefined curve, we compute the image of intermediary points
                    this->GetPointsOnCurvedFace(n, ref, XiTri, points_interm);
                  }
                else
                  {
                    // curved mesh, we perform an interpolation
                    for (int j = 0; j < ValPhiTri.GetN(); j++)
                      {
                        points_interm(j).Zero();
                        for (int i = 0; i < ValPhiTri.GetM(); i++)
                          Add(ValPhiTri(i, j), OldPosTri(i), points_interm(j));
                        
                      }
                  }
                
                // we create interior edges and affect internal points
                nb = 0;
                for (int i = 0; i < r-1; i++)
                  for (int j = 0; j < r-1-i; j++)
                    {
                      // we consider the triangle
                      // (i+1, j) -> (i+1, j+1) -> (i, j+1)
                      int n1 = NumNodes2D(0)(i+1, j);
                      int n2 = NumNodes2D(0)(i+1, j+1);
                      int n3 = NumNodes2D(0)(i, j+1);
                      
                      new_edge(nb_new_edges_ref).Init(NodleSurf(n)(n1), NodleSurf(n)(n2), ref);
                      for (int k = 0; k < rg-1; k++)
                        NewPoints(nb_new_edges_ref, k) = points_interm(nb++);
                      
                      nb_new_edges_ref++;
                      new_edge(nb_new_edges_ref).Init(NodleSurf(n)(n2), NodleSurf(n)(n3), ref);
                      for (int k = 0; k < rg-1; k++)
                        NewPoints(nb_new_edges_ref, k) = points_interm(nb++);
                      
                      nb_new_edges_ref++;
                      new_edge(nb_new_edges_ref).Init(NodleSurf(n)(n3), NodleSurf(n)(n1), ref);
                      for (int k = 0; k < rg-1; k++)
                        NewPoints(nb_new_edges_ref, k) = points_interm(nb++);
                      
                      nb_new_edges_ref++;
                    }
                
                // we affect points on faces
                // creating of points inside sub-triangles of triangle
                for (int i = 0; i < r; i++)
                  {
                    for (int j = 0; j < r-1-i; j++)
                      {
                        // we consider the triangle
                        // (i+1, j) -> (i+1, j+1) -> (i, j+1)
                        //int n1 = NumNodes2D(0)(i+1, j);
                        //int n2 = NumNodes2D(0)(i+1, j+1);
                        //int n3 = NumNodes2D(0)(i, j+1);
                        
                        NewPointsF(nb_new_faces_ref).Reallocate(Ntri);
                        for (int k = 0; k < Ntri; k++)
                          NewPointsF(nb_new_faces_ref)(k) = points_interm(nb++);
                        
                        // and triangle
                        // (i, j) -> (i+1, j) -> (i, j+1)
                        
                        //n1 = NumNodes2D(0)(i, j);
                        //n2 = NumNodes2D(0)(i+1, j);
                        //n3 = NumNodes2D(0)(i, j+1);
                        
                        nb_new_faces_ref++;
                        NewPointsF(nb_new_faces_ref).Reallocate(Ntri);
                        for (int k = 0; k < Ntri; k++)
                          NewPointsF(nb_new_faces_ref)(k) = points_interm(nb++);
                        
                        nb_new_faces_ref++;
                      }

            
                    // last triangle (i, r-1-i) -> (i+1, r-1-i) -> (i, r-i)
                    NewPointsF(nb_new_faces_ref).Reallocate(Ntri);
                    for (int k = 0; k < Ntri; k++)
                      NewPointsF(nb_new_faces_ref)(k) = points_interm(nb++);
                    
                    nb_new_faces_ref++;
                  }
	      }
            else
              {
                this->GetNodalPointsOnFace(n, OldPosQuad);
                
                Vector<R_N> points_interm(Nquad);
                if (this->GetCurveType(ref) != this->CURVE_FILE )
                  {
                    // predefined curve, we compute the image of intermediary points
                    this->GetPointsOnCurvedFace(n, ref, XiQuad, points_interm);
                  }
                else
                  {
                    points_interm.Reallocate(ValPhiQuad.GetN());
                    // curved mesh, we perform an interpolation
                    for (int j = 0; j < ValPhiQuad.GetN(); j++)
                      {
                        points_interm(j).Zero();
                        for (int i = 0; i < ValPhiQuad.GetM(); i++)
                          Add(ValPhiQuad(i, j), OldPosQuad(i), points_interm(j));
                      }
                  }
                
                // creation of sub-edges of quadrangles
                nb = 0;
                for (int i = 0; i < r-1; i++)
                  for (int j = 0; j < r; j++)
                    {
                      // we consider edge (i+1, j) -> (i+1, j+1)
                      int n1 = NumNodes2D(1)(i+1, j);
                      int n2 = NumNodes2D(1)(i+1, j+1);
                      
                      new_edge(nb_new_edges_ref).Init(NodleSurf(n)(n1), NodleSurf(n)(n2), ref);
                      for (int k = 0; k < rg-1; k++)
                        NewPoints(nb_new_edges_ref, k) = points_interm(nb++);
                      
                      nb_new_edges_ref++;
                      
                      // we consider edge (j, i+1) -> (j+1, i+1)
                      n1 = NumNodes2D(1)(j, i+1);
                      n2 = NumNodes2D(1)(j+1, i+1);
                      
                      new_edge(nb_new_edges_ref).Init(NodleSurf(n)(n1), NodleSurf(n)(n2), ref);
                      for (int k = 0; k < rg-1; k++)
                        NewPoints(nb_new_edges_ref, k) = points_interm(nb++);
                      
                      nb_new_edges_ref++;
                    }
              
                // affecting points inside sub-quadrilaterals of quadrangle
                for (int i = 0; i < r; i++)
                  {
                    for (int j = 0; j < r; j++)
                      {
                        // we consider the quadrangle
                        // (i, j) -> (i+1, j) -> (i+1, j+1) -> (i, j+1)
                        
                        NewPointsF(nb_new_faces_ref).Reallocate(Nquad);
                        for (int k = 0; k < Nquad; k++)
                          NewPointsF(nb_new_faces_ref)(k) = points_interm(nb++);
                        
                        nb_new_faces_ref++;
                      }
                  }
              }
	  }
      }
    
    this->ReallocateEdgesRef(nb_new_edges_ref);    
    for (int i = 0; i < new_edge.GetM(); i++)
      this->EdgeRef(i) = new_edge(i);
    
    Vector<Face<Dimension3> > old_face = this->faces_ref;
    this->ReallocateBoundariesRef(nb_new_faces_ref);
    
    if (rg > 1)
      {
        this->PointsEdgeRef = NewPoints;    
        this->PointsFaceRef = NewPointsF;
      }
    
    // then creation of new referenced faces
    nb_new_faces_ref = 0;
    for (int n = 0; n < nb_old_faces_ref; n++)
      {
	int ref = old_face(n).GetReference();
	
        if (old_face(n).IsQuadrangular())
          {
	    // we create now sub-quadrilaterals
            for (int i = 0; i < r; i++)
              {
                for (int j = 0; j < r; j++)
                  {
                    // we consider the quadrangle
                    // (i, j) -> (i+1, j) -> (i+1, j+1) -> (i, j+1)
                    int n1 = NumNodes2D(1)(i, j);
                    int n2 = NumNodes2D(1)(i+1, j);
                    int n3 = NumNodes2D(1)(i+1, j+1);
                    int n4 = NumNodes2D(1)(i, j+1);
                    
                    this->BoundaryRef(nb_new_faces_ref)
                      .InitQuadrangular(NodleSurf(n)(n1), NodleSurf(n)(n2),
                                        NodleSurf(n)(n3), NodleSurf(n)(n4), ref);
                    
                    nb_new_faces_ref++;
                  }
              }                
	  }
	else
	  {
            // we affect points on faces
            // creating of sub-triangles of triangle
            for (int i = 0; i < r; i++)
              {
                for (int j = 0; j < r-1-i; j++)
                  {
                    // we consider the triangle
                    // (i+1, j) -> (i+1, j+1) -> (i, j+1)
                    int n1 = NumNodes2D(0)(i+1, j);
                    int n2 = NumNodes2D(0)(i+1, j+1);
                    int n3 = NumNodes2D(0)(i, j+1);
                    
                    this->BoundaryRef(nb_new_faces_ref).
                      InitTriangular(NodleSurf(n)(n1), NodleSurf(n)(n2), NodleSurf(n)(n3), ref);
                    
                    nb_new_faces_ref++;
                    
                    // and triangle
                    // (i, j) -> (i+1, j) -> (i, j+1)
                    
                    n1 = NumNodes2D(0)(i, j);
                    n2 = NumNodes2D(0)(i+1, j);
                    n3 = NumNodes2D(0)(i, j+1);
                        
                    this->BoundaryRef(nb_new_faces_ref).
                      InitTriangular(NodleSurf(n)(n1), NodleSurf(n)(n2), NodleSurf(n)(n3), ref);
                                        
                    nb_new_faces_ref++;
                  }

                
                // last triangle (i, r-1-i) -> (i+1, r-1-i) -> (i, r-i)
                int n1 = NumNodes2D(0)(i, r-1-i);
                int n2 = NumNodes2D(0)(i+1, r-1-i);
                int n3 = NumNodes2D(0)(i, r-i);
                
                this->BoundaryRef(nb_new_faces_ref).
                  InitTriangular(NodleSurf(n)(n1), NodleSurf(n)(n2), NodleSurf(n)(n3), ref);
                                        
                nb_new_faces_ref++;
              }
          }
      }
    
    // sorting the new edges and faces
    this->SortEdgesRef();
    this->SortBoundariesRef();
    
    old_face.Clear();
    Vector<Volume> old_elt = this->elements;
    
    // dividing elements
    int nb_new_elt = (nb_old_elt-nb_pyramids)*r*r*r + nb_pyramids*(4*r*r*r-r)/3;
    this->ReallocateElements(nb_new_elt);
    nb_new_elt = 0; TinyVector<int, 8> nx;
    for (int i = 0; i < nb_old_elt; i++)
      {
	int ref = old_elt(i).GetReference();
        const IVect& num = Nodle(i);
        
        switch (old_elt(i).GetHybridType())
          {
          case 0 :
            // now we create small tetrahedra
            for (int j = 0; j < r; j++)
              for (int k = 0; k < r-j; k++)
                for (int n = 0; n < r-j-k; n++)
                  {
                    nx.Fill(-1);
                    if (NumNodes3D(0)(j, k, n) != -1)
                      nx(0) = num(NumNodes3D(0)(j, k, n));
                    
                    if (NumNodes3D(0)(j, k, n+1) != -1)
                      nx(1) = num(NumNodes3D(0)(j, k, n+1));
                    
                    if (NumNodes3D(0)(j, k+1, n) != -1)
                      nx(2) = num(NumNodes3D(0)(j, k+1, n));
                    
                    if (NumNodes3D(0)(j, k+1, n+1)!=-1)
                      nx(3) = num(NumNodes3D(0)(j,k+1,n+1));
                    
                    if (NumNodes3D(0)(j+1, k, n) != -1)
                      nx(4) = num(NumNodes3D(0)(j+1,k,n));
                    
                    if (NumNodes3D(0)(j+1, k, n+1) != -1)
                      nx(5) = num(NumNodes3D(0)(j+1, k, n+1));
                    
                    if (NumNodes3D(0)(j+1, k+1, n) != -1)
                      nx(6) = num(NumNodes3D(0)(j+1, k+1, n));
                    
                    if (NumNodes3D(0)(j+1, k+1, n+1)!=-1)
                      nx(7) = num(NumNodes3D(0)(j+1, k+1, n+1));
		    
		    if ((nx(0) != -1)&&(nx(1) != -1))
		      {
			// first tetrahedron
			this->Element(nb_new_elt).InitTetrahedral(nx(0), nx(4), nx(1), nx(2), ref);
			nb_new_elt++;
			
			if (nx(3) != -1)
			  {
			    // four other tetrahedra
			    this->Element(nb_new_elt)
                              .InitTetrahedral(nx(1), nx(4), nx(2), nx(5), ref); nb_new_elt++;
			    this->Element(nb_new_elt)
                              .InitTetrahedral(nx(1), nx(2), nx(3), nx(5), ref); nb_new_elt++;
			    this->Element(nb_new_elt)
                              .InitTetrahedral(nx(4), nx(6), nx(2), nx(5), ref); nb_new_elt++;
			    this->Element(nb_new_elt)
                              .InitTetrahedral(nx(2), nx(6), nx(3), nx(5), ref); nb_new_elt++;
			    
			    if (nx(7) != -1)
			      {
				// last tetrahedron
				this->Element(nb_new_elt)
                                  .InitTetrahedral(nx(5), nx(6), nx(3), nx(7), ref); nb_new_elt++;
			      }
			  }
		      }
                    
		  }
	    break;
	  case 1 :
	    // splitting of a pyramid
	    for (int j = 0; j <= nb_div; j++)
	      for (int k = 0; k <= nb_div; k++)
		for (int l = 0; l <= nb_div; l++)
		  {
		    // small hex
		    nx.Fill(-1);
		    if (NumNodes3D(1)(j, k, l) != -1)
		      nx(0) = num(NumNodes3D(1)(j, k, l));
		    
		    if (NumNodes3D(1)(j+1, k, l) != -1)
		      nx(1) = num(NumNodes3D(1)(j+1, k, l));
		    
		    if (NumNodes3D(1)(j+1, k+1, l) != -1)
		      nx(2) = num(NumNodes3D(1)(j+1, k+1, l));
		    
		    if (NumNodes3D(1)(j, k+1, l) != -1)
		      nx(3) = num(NumNodes3D(1)(j, k+1, l));
		    
		    if (NumNodes3D(1)(j, k, l+1) != -1)
		      nx(4) = num(NumNodes3D(1)(j, k, l+1));
		    
		    if (NumNodes3D(1)(j+1, k, l+1) != -1)
		      nx(5) = num(NumNodes3D(1)(j+1, k, l+1));
		    
		    if (NumNodes3D(1)(j+1, k+1, l+1) != -1)
		      nx(6) = num(NumNodes3D(1)(j+1, k+1, l+1));
		    
		    if (NumNodes3D(1)(j, k+1, l+1) != -1)
		      nx(7) = num(NumNodes3D(1)(j, k+1, l+1));
		    
		    // we try to create each subelement if all nodes are existing
		    if ((nx(0) != -1)&&(nx(1) != -1)&&(nx(2) != -1)
			&&(nx(3) != -1)&&(nx(4) != -1))
		      this->Element(nb_new_elt++)
                        .InitPyramidal(nx(0), nx(1), nx(2), nx(3), nx(4), ref);
		    
		    if ((nx(4) != -1)&&(nx(5) != -1)&&(nx(6) != -1)
			&&(nx(7) != -1)&&(nx(2) != -1))
		      this->Element(nb_new_elt++)
                        .InitPyramidal(nx(4), nx(5), nx(6), nx(7), nx(2), ref);

		    if ((nx(1) != -1)&&(nx(2) != -1)&&(nx(4) != -1)&&(nx(5) != -1))
		      this->Element(nb_new_elt++)
                        .InitTetrahedral(nx(1), nx(2), nx(4), nx(5), ref);

		    if ((nx(3) != -1)&&(nx(2) != -1)&&(nx(4) != -1)&&(nx(7) != -1))
		      this->Element(nb_new_elt++)
                        .InitTetrahedral(nx(3), nx(2), nx(4), nx(7), ref);
		  }
	    break;
	  case 2:
	    // splitting of a prism
	    for (int j = 0; j <= nb_div; j++)
	      for (int k = 0; k <= nb_div-j; k++)
		for (int l = 0; l <= nb_div; l++)
		  {
		    // small hex
		    nx.Fill(-1);
		    if (NumNodes3D(2)(j, k, l) != -1)
		      nx(0) = num(NumNodes3D(2)(j, k, l));
		    
		    if (NumNodes3D(2)(j+1, k, l) != -1)
		      nx(1) = num(NumNodes3D(2)(j+1, k, l));
		    
		    if (NumNodes3D(2)(j+1, k+1, l) != -1)
		      nx(2) = num(NumNodes3D(2)(j+1, k+1, l));
		    
		    if (NumNodes3D(2)(j, k+1, l) != -1)
		      nx(3) = num(NumNodes3D(2)(j, k+1, l));
		    
		    if (NumNodes3D(2)(j, k, l+1) != -1)
		      nx(4) = num(NumNodes3D(2)(j, k, l+1));
		    
		    if (NumNodes3D(2)(j+1, k, l+1) != -1)
		      nx(5) = num(NumNodes3D(2)(j+1, k, l+1));
		    
		    if (NumNodes3D(2)(j+1, k+1, l+1) != -1)
		      nx(6) = num(NumNodes3D(2)(j+1, k+1, l+1));
		    
		    if (NumNodes3D(2)(j, k+1, l+1) != -1)
		      nx(7) = num(NumNodes3D(2)(j, k+1, l+1));
		    
		    // we try to create each subelement if all nodes are existing
		    if ((nx(0) != -1)&&(nx(1) != -1)&&(nx(3) != -1)
			&&(nx(4) != -1)&&(nx(5) != -1)&&(nx(7) != -1))
		      this->Element(nb_new_elt++)
                        .InitWedge(nx(0), nx(1), nx(3), nx(4), nx(5), nx(7), ref);
		    
		    if ((nx(1) != -1)&&(nx(2) != -1)&&(nx(3) != -1)
			&&(nx(5) != -1)&&(nx(6) != -1)&&(nx(7) != -1))
                      this->Element(nb_new_elt++)
                        .InitWedge(nx(1), nx(2), nx(3), nx(5), nx(6), nx(7), ref);
		  }
	    break;
	  case 3 :
	    for (int j = 0; j <= nb_div; j++)
	      for (int k = 0; k <= nb_div; k++)
		for (int l = 0; l <= nb_div; l++)
		  {
		    // now we create small hexahedra
		    nx(0) = num(NumNodes3D(3)(j,k,l));
		    nx(1) = num(NumNodes3D(3)(j+1,k,l));
		    nx(2) = num(NumNodes3D(3)(j+1,k+1,l));
		    nx(3) = num(NumNodes3D(3)(j,k+1,l));
		    nx(4) = num(NumNodes3D(3)(j,k,l+1));
		    nx(5) = num(NumNodes3D(3)(j+1,k,l+1));
		    nx(6) = num(NumNodes3D(3)(j+1,k+1,l+1));
		    nx(7) = num(NumNodes3D(3)(j,k+1,l+1));
		    this->Element(nb_new_elt).InitHexahedral(nx(0), nx(1), nx(2), nx(3),
                                                             nx(4), nx(5), nx(6), nx(7), ref);
		    nb_new_elt++;
		  }
	    break;
	  }
      }    
    
    this->ReorientElements();
    FindConnectivity();    
    this->ProjectPointsOnCurves();
  }    

  
  /*****************
   * Mesh creation *
   *****************/
  
  
  //! creation of a regular mesh of a parallelepipedic domain
  //! [xmin xmax ymin ymax zmin zmax]
  /*!
    \param[in] ptA point [xmin ymin zmin]
    \param[in] ptB point [xmax ymax zmax]
    \param[in] nbx the number of point along the edge xmin -> xmax
    \param[in] ref_domain reference of the inside domain
    \param[in] ref_boundary  reference of the six faces of the parallelepiped
    (x = xmin, y = min, z = zmin, z = zmax, y = ymax and x = xmax)
    \param[in] type_mesh hexahedral, tetrahedral or hybrid mesh
  */
  template<class T>
  void Mesh<Dimension3, T>::
  CreateRegularMesh(const R_N& ptA, const R_N& ptB, const TinyVector<int, 3>& nb_points,
		    int ref_domain, const TinyVector<int, 6>& ref_boundary, int type_mesh,
                    R3 ratio)
  {
    // we clear previous mesh
    Clear();
    
    int ref_max = ref_boundary(0);
    for (int i = 1; i < 6; i++)
      ref_max = max(ref_max, ref_boundary(i));
    
    if (ref_max > this->GetNbReferences())
      this->ResizeNbReferences(ref_max);
    
    Real_wp xmin = ptA(0);
    Real_wp ymin = ptA(1);
    Real_wp zmin = ptA(2);
    Real_wp xmax = ptB(0);
    Real_wp ymax = ptB(1);
    Real_wp zmax = ptB(2);
    
    if (this->print_level >= 2)
      cout<< " we create a regular mesh "<<endl;
    
    if (xmax <= xmin)
      return;
    
    if (ymax <= ymin)
      return;

    if (zmax <= zmin)
      return;
    
    bool geom_x = false;
    if (ratio(0) != 1.0)
      geom_x = true;
    
    bool geom_y = false;
    if (ratio(1) != 1.0)
      geom_y = true;

    bool geom_z = false;
    if (ratio(2) != 1.0)
      geom_z = true;
    
    int nbx = nb_points(0), nby = nb_points(1), nbz = nb_points(2);
    T step_x = (xmax-xmin)/(nbx-1);
    T step_y = (ymax-ymin)/(nby-1);
    T step_z = (zmax-zmin)/(nbz-1);
    if (geom_x)
      step_x = (xmax-xmin)*(1.0-ratio(0))/(1.0 - pow(ratio(0), nbx-1));
    
    if (geom_y)
      step_y = (ymax-ymin)*(1.0-ratio(1))/(1.0 - pow(ratio(1), nby-1));

    if (geom_z)
      step_z = (zmax-zmin)*(1.0-ratio(2))/(1.0 - pow(ratio(2), nbz-1));
    
    T pos_x, pos_y, pos_z;
    int nb_vertices = nbx*nby*nbz;
    int nb_elt = (nbx-1)*(nby-1)*(nbz-1);
    if ((type_mesh == this->PYRAMID_MESH)||(type_mesh == this->CRUSH_MESH))
      nb_vertices += nb_elt;
    
    this->ReallocateVertices(nb_vertices);
    
    int offset_hyb = nb_elt*2;
    if (type_mesh == this->TETRAHEDRAL_MESH)
      nb_elt *= 6;
    if (type_mesh == this->PYRAMID_MESH)
      nb_elt *= 6;
    if (type_mesh == this->WEDGE_MESH)
      nb_elt *= 2;
    if (type_mesh == this->HYBRID_MESH)
      nb_elt *= 4;
    if (type_mesh == this->CRUSH_MESH)
      nb_elt *= 6;
    
    this->ReallocateElements(nb_elt);
    Array3D<int> num(2, 2, 2);
    num(0, 0, 0) = 0;
    num(1, 0, 0) = 1;
    num(1, 1, 0) = nbx+1;
    num(0, 1, 0) = nbx;
    num(0, 0, 1) = nbx*nby;
    num(1, 0, 1) = nbx*nby + 1;
    num(1, 1, 1) = nbx*nby + nbx + 1;
    num(0, 1, 1) = nbx*nby + nbx;
    int nb_faces_ref = 2*(nbx-1)*(nby-1) + 2*(nbx-1)*(nbz-1) + 2*(nby-1)*(nbz-1);
    if (type_mesh == this->TETRAHEDRAL_MESH)
      nb_faces_ref *= 2;
    else if (type_mesh == this->HYBRID_MESH)
      nb_faces_ref = 2*(nbx-1)*(nby-1) + 4*(nbz-1)*(nbx+nby-2);
    else if (type_mesh == this->WEDGE_MESH)
      nb_faces_ref = 2*(nbz-1)*(nby-1) + 2*(nbz-1)*(nbx-1) + 4*(nbx-1)*(nby-1);

    int offset_face_xmin = 0, offset_face_xmax = 0, offset_face_ymin = 0,
      offset_face_ymax = 0, offset_face_zmin = 0, offset_face_zmax = 0;
    
    if ( (type_mesh == this->TETRAHEDRAL_MESH) || (type_mesh == this->HYBRID_MESH))
      {
        offset_face_xmax = offset_face_xmin + 2*(nby-1)*(nbz-1);
        offset_face_ymin = offset_face_xmax + 2*(nby-1)*(nbz-1);
        offset_face_ymax = offset_face_ymin + 2*(nbx-1)*(nbz-1);
        offset_face_zmin = offset_face_ymax + 2*(nbx-1)*(nbz-1);
      } 
    else
      {
        offset_face_xmax = offset_face_xmin + (nby-1)*(nbz-1);
        offset_face_ymin = offset_face_xmax + (nby-1)*(nbz-1);
        offset_face_ymax = offset_face_ymin + (nbx-1)*(nbz-1);
        offset_face_zmin = offset_face_ymax + (nbx-1)*(nbz-1);
      }
    
    if ( (type_mesh == this->TETRAHEDRAL_MESH) || (type_mesh == this->WEDGE_MESH))
      offset_face_zmax = offset_face_zmin + 2*(nbx-1)*(nby-1);
    else
      offset_face_zmax = offset_face_zmin + (nbx-1)*(nby-1);
    
    ReallocateBoundariesRef(nb_faces_ref);
    nb_elt = 0; nb_faces_ref = 0; int nb_cell = nbx*nby*nbz;
    for (int k = 0; k < nbz; k++)
      for (int j = 0; j < nby; j++)
	for (int i = 0; i < nbx; i++)
	  {
	    if (i == (nbx-1))
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
                  pos_x = xmin + i*step_x;
              }
            
	    if (j == (nby-1))
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
                  pos_y = ymin + j*step_y;
              }
            
	    if (k == (nbz-1))
	      pos_z = zmax;
	    else
	      {
                if (geom_z)
                  {
                    if (k == 0)
                      pos_z = zmin;
                    else
                      pos_z = zmin + step_z*(1.0-pow(ratio(2), k))/(1.0-ratio(2));
                  }
                else
                  pos_z = zmin + k*step_z;
              }
            
	    if ((i<(nbx-1))&&(j<(nby-1))&&(k<(nbz-1)))
	      {
		if (type_mesh == this->TETRAHEDRAL_MESH)
		  {
		    this->elements(nb_elt++).
		      InitTetrahedral(num(0,0,0), num(1,0,0),
				      num(0,1,0), num(0,0,1), ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(0,1,0), num(0,0,1),
				      num(1,0,0), num(1,0,1), ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(0,0,1),num(0,1,0),num(0,1,1),num(1,0,1),ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(1,1,0), num(1,1,1),
				      num(0,1,1), num(1,0,1), ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(1,0,0), num(1,1,0),
				      num(0,1,0), num(1,0,1), ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(0,1,0), num(1,0,1),
				      num(1,1,0), num(0,1,1), ref_domain);
		  }
		else if (type_mesh == this->HYBRID_MESH)
		  {
		    this->elements(offset_hyb++).
		      InitPyramidal(num(0,0,0), num(0,1,0), num(1,1,0),
				    num(1,0,0), num(0,1,1), ref_domain);
		    
		    this->elements(offset_hyb++).
		      InitPyramidal(num(0,0,1), num(0,1,1), num(1,1,1),
				    num(1,0,1), num(1,0,0), ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(1,1,0), num(1,0,0),
				      num(1,1,1), num(0,1,1),ref_domain);
		    
		    this->elements(nb_elt++).
		      InitTetrahedral(num(0,0,1), num(0,1,1),
				      num(0,0,0), num(1,0,0),ref_domain);
		  }
		else if (type_mesh == this->PYRAMID_MESH)
		  {
		    this->elements(nb_elt++).
		      InitPyramidal(num(0,0,0), num(0,1,0), num(0,1,1),
				    num(0,0,1), nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitPyramidal(num(1,0,0), num(1,1,0), num(1,1,1),
				    num(1,0,1), nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitPyramidal(num(0,0,0), num(1,0,0), num(1,0,1), num(0,0,1),
				    nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitPyramidal(num(0,1,0), num(1,1,0),
				    num(1,1,1), num(0,1,1), nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitPyramidal(num(0,0,0), num(1,0,0),
				    num(1,1,0), num(0,1,0), nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitPyramidal(num(0,0,1), num(1,0,1),
				    num(1,1,1), num(0,1,1), nb_cell, ref_domain);
		  }
		else if (type_mesh == this->WEDGE_MESH)
		  {
		    this->elements(nb_elt++).
		      InitWedge(num(0,0,0), num(1,0,0), num(0,1,0),
				num(0,0,1), num(1,0,1), num(0,1,1), ref_domain);

		    this->elements(nb_elt++).
		      InitWedge(num(1,0,0), num(1,1,0), num(0,1,0),
				num(1,0,1), num(1,1,1), num(0,1,1), ref_domain);
		  }
		else if (type_mesh == this->CRUSH_MESH)
		  {
		    this->elements(nb_elt++).
		      InitHexahedral(num(0,0,0), num(0,1,0), num(0,1,1),
				     num(0,0,1), nb_cell, nb_cell,
				     nb_cell, nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitHexahedral(num(1,0,0), num(1,1,0), num(1,1,1),
				     num(1,0,1), nb_cell, nb_cell,
				     nb_cell, nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitHexahedral(num(0,0,0), num(1,0,0), num(1,0,1), 
		                     num(0,0,1), nb_cell, nb_cell,
		                     nb_cell, nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitHexahedral(num(0,1,0), num(1,1,0), num(1,1,1),
				     num(0,1,1), nb_cell, nb_cell,
				     nb_cell, nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitHexahedral(num(0,0,0), num(1,0,0), num(1,1,0),
				     num(0,1,0), nb_cell, nb_cell,
				     nb_cell, nb_cell, ref_domain);
		    
		    this->elements(nb_elt++).
		      InitHexahedral(num(0,0,1), num(1,0,1), num(1,1,1),
				     num(0,1,1), nb_cell, nb_cell,
				     nb_cell, nb_cell, ref_domain);
		  }
		else
		  this->elements(nb_elt++).
		    InitHexahedral(num(0,0,0), num(1,0,0), num(1,1,0),
				   num(0,1,0), num(0,0,1), num(1,0,1),
				   num(1,1,1), num(0,1,1), ref_domain);
		
		// middle point for mesh made of pyramids only
		if ((type_mesh == this->PYRAMID_MESH)||(type_mesh == this->CRUSH_MESH))
                  this->Vertices(nb_cell++).Init(pos_x+0.5*step_x, pos_y+0.5*step_y,
                                                 pos_z+0.5*step_z);
		   
		
		// boundary x = xmin
		if (i==0)
		  {
		    if (type_mesh == this->TETRAHEDRAL_MESH)
		      {
			faces_ref(offset_face_xmin++).
			  InitTriangular(num(0,0,0), num(0,1,0),
					 num(0,0,1), ref_boundary(0));
			
			faces_ref(offset_face_xmin++).
			  InitTriangular(num(0,1,0), num(0,1,1),
					 num(0,0,1), ref_boundary(0));
		      }
		    else if (type_mesh == this->HYBRID_MESH)
		      {
		        faces_ref(nb_faces_ref++).
			  InitTriangular(num(0,0,0), num(0,0,1),
					 num(0,1,1), ref_boundary(0));
		        
			faces_ref(nb_faces_ref++).
			  InitTriangular(num(0,0,0), num(0,1,0),
					 num(0,1,1), ref_boundary(0));
		      }
		    else
		      faces_ref(nb_faces_ref++).
			InitQuadrangular(num(0,0,0), num(0,1,0),
					 num(0,1,1), num(0,0,1), ref_boundary(0));
		  }
		
		// boundary x = xmax
		if (i==(nbx-2))
		  {
		    if (type_mesh == this->TETRAHEDRAL_MESH)
		      {
			faces_ref(offset_face_xmax++).
			  InitTriangular(num(1,0,0), num(1,1,0),
					 num(1,0,1), ref_boundary(5));
			
			faces_ref(offset_face_xmax++).
			  InitTriangular(num(1,1,0), num(1,1,1),
					 num(1,0,1), ref_boundary(5));
		      }
		    else if (type_mesh == this->HYBRID_MESH)
		      {
		        faces_ref(nb_faces_ref++).
			  InitTriangular(num(1,0,0), num(1,0,1),
					 num(1,1,1), ref_boundary(5));
		        faces_ref(nb_faces_ref++).
			  InitTriangular(num(1,0,0), num(1,1,0),
					 num(1,1,1), ref_boundary(5));
		      }
		    else
		      faces_ref(nb_faces_ref++).
			InitQuadrangular(num(1,0,0), num(1,1,0),
					 num(1,1,1), num(1,0,1), ref_boundary(5));
		  }
		
		// boundary y = ymin
		if (j==0)
		  {
		    if (type_mesh == this->TETRAHEDRAL_MESH)
		      {
			faces_ref(offset_face_ymin++).
			  InitTriangular(num(0,0,0), num(1,0,0),
					 num(0,0,1), ref_boundary(1));
			
			faces_ref(offset_face_ymin++).
			  InitTriangular(num(1,0,0), num(1,0,1),
					 num(0,0,1), ref_boundary(1));
		      }
		    else if (type_mesh == this->HYBRID_MESH)
		      {
		        faces_ref(nb_faces_ref++).
			  InitTriangular(num(0,0,1), num(1,0,0),
					 num(1,0,1), ref_boundary(1));
			
		        faces_ref(nb_faces_ref++).
			  InitTriangular(num(0,0,0), num(0,0,1),
					 num(1,0,0), ref_boundary(1));
		      }
		    else
		      faces_ref(nb_faces_ref++).InitQuadrangular(num(0,0,0),num(1,0,0),num(1,0,1),
                                                                 num(0,0,1),ref_boundary(1));
		  }
		
		// boundary y = ymax
		if (j==(nby-2))
		  {
		    if (type_mesh == this->TETRAHEDRAL_MESH)
		      {
			faces_ref(offset_face_ymax++)
                          .InitTriangular(num(0,1,0),num(1,1,0),num(0,1,1),ref_boundary(4));
			faces_ref(offset_face_ymax++)
                          .InitTriangular(num(1,1,0),num(1,1,1),num(0,1,1),ref_boundary(4));
		      }
		    else if (type_mesh == this->HYBRID_MESH)
		      {
		        faces_ref(nb_faces_ref++)
                          .InitTriangular(num(0,1,0),num(0,1,1),num(1,1,0),ref_boundary(4));
		        faces_ref(nb_faces_ref++)
                          .InitTriangular(num(0,1,1),num(1,1,0),num(1,1,1),ref_boundary(4));
		      }
		    else
		      faces_ref(nb_faces_ref++)
                        .InitQuadrangular(num(0,1,0),num(1,1,0),
                                          num(1,1,1),num(0,1,1),ref_boundary(4));
		  }
		
		// boundary z = zmin
		if (k==0)
		  {
		    if ((type_mesh == this->TETRAHEDRAL_MESH)||(type_mesh == this->WEDGE_MESH))
		      {
			faces_ref(offset_face_zmin++)
                          .InitTriangular(num(0,0,0),num(1,0,0),num(0,1,0),ref_boundary(2));
			faces_ref(offset_face_zmin++)
                          .InitTriangular(num(1,0,0),num(1,1,0),num(0,1,0),ref_boundary(2));
		      }
		    else
		      faces_ref(nb_faces_ref++).InitQuadrangular(num(0,0,0),num(1,0,0),num(1,1,0),
                                                                 num(0,1,0),ref_boundary(2));
		  }
		
		// boundary z = zmax
		if (k==(nbz-2))
		  {
		    if ((type_mesh == this->TETRAHEDRAL_MESH)||(type_mesh == this->WEDGE_MESH))
		      {
			faces_ref(offset_face_zmax++)
                          .InitTriangular(num(0,0,1),num(1,0,1),num(0,1,1),ref_boundary(3));
			faces_ref(offset_face_zmax++)
                          .InitTriangular(num(1,0,1),num(1,1,1),num(0,1,1),ref_boundary(3));
		      }
		    else
		      faces_ref(nb_faces_ref++).InitQuadrangular(num(0,0,1),num(1,0,1),num(1,1,1),
                                                                 num(0,1,1),ref_boundary(3));
		  }
	      }
	    
	    if (regular_mesh_with_slight_modification == REGULAR_MESH_SLIGHTLY_MODIFIED)
	      {
                if ((i%2 == 1)&&(j%2 == 1)&&(k%2 == 1))
                  {
                    pos_x += 0.18*step_x;
                    pos_y += 0.2*step_y;
                    pos_z += 0.15*step_z;
                  }
              }
            else if (regular_mesh_with_slight_modification == REGULAR_MESH_RANDOMLY_MODIFIED)
              {
                pos_x += 0.18*(1.0 - 2.0*rand()/RAND_MAX)*step_x;
                pos_y += 0.2*(1.0 - 2.0*rand()/RAND_MAX)*step_y;
                pos_z += 0.15*(1.0 - 2.0*rand()/RAND_MAX)*step_z;
              }
	    
	    this->Vertices(num(0,0,0)).Init(pos_x, pos_y, pos_z);
	    
	    for (int l = 0; l < 2; l++)
	      for (int m = 0; m < 2; m++)
		for (int n = 0; n < 2; n++)
		  num(l, m, n)++;
	  }

    
    if (type_mesh == this->TETRAHEDRAL_MESH)
      {
        offset_face_xmin = 0;
        offset_face_xmax = nbx*nbz*(nby-1);
        offset_face_ymin = offset_face_xmax + (nbx-1)*nby*nbz;
        offset_face_ymax = offset_face_ymin + (nbz-1)*nbx*nby;
        offset_face_zmin = offset_face_ymax + (nby-1)*(nbz-1)*nbx;
        offset_face_zmax = offset_face_zmin + (nbx-1)*(nbz-1)*nby;
        
        num(0, 0, 0) = 0;
        num(1, 0, 0) = 1;
        num(1, 1, 0) = nbx+1;
        num(0, 1, 0) = nbx;
        num(0, 0, 1) = nbx*nby;
        num(1, 0, 1) = nbx*nby + 1;
        num(1, 1, 1) = nbx*nby + nbx + 1;
        num(0, 1, 1) = nbx*nby + nbx;
        
        this->ReallocateEdgesRef(offset_face_zmax + (nbx-1)*(nby-1)*nbz);
                                 
        for (int k = 0; k < nbz; k++)
          for (int j = 0; j < nby; j++)
            for (int i = 0; i < nbx; i++)
              {                
                if (j < nby-1)
                  this->EdgeRef(offset_face_xmin++).Init(num(0, 0, 0), num(0, 1, 0), 1);
                
                if (i < nbx-1)
                  this->EdgeRef(offset_face_xmax++).Init(num(0, 0, 0), num(1, 0, 0), 1);
                
                if (k < nbz-1)
                  this->EdgeRef(offset_face_ymin++).Init(num(0, 0, 0), num(0, 0, 1), 1);
                
                if ((j < nby-1) && (k < nbz-1))
                  this->EdgeRef(offset_face_ymax++).Init(num(0, 1, 0), num(0, 0, 1), 1);

                if ((i < nbx-1) && (k < nbz-1))
                  this->EdgeRef(offset_face_zmin++).Init(num(1, 0, 0), num(0, 0, 1), 1);
                
                if ((i < nbx-1) && (j < nby-1))
                  this->EdgeRef(offset_face_zmax++).Init(num(1, 0, 0), num(0, 1, 0), 1);

                for (int l = 0; l < 2; l++)
                  for (int m = 0; m < 2; m++)
                    for (int n = 0; n < 2; n++)
                      num(l, m, n)++;                
              }
      }
    

    this->SortPeriodicBoundaries();
    this->SortBoundariesRef();
    
    this->ReorientElements();
    this->FindConnectivity();
    this->ProjectPointsOnCurves();
    // this->Write("cube.mesh");
  }


  template<class T>
  void Mesh<Dimension3, T>::
  CreateStructuredMesh(const Vector<T>& pos_layer_x, const Vector<T>& pos_layer_y, const Vector<T>& pos_layer_z,
                       Vector<int>& nb_layers_x, Vector<int>& nb_layers_y, Vector<int>& nb_layers_z,
                       const Array3D<int>& ref_domain_layer, const TinyVector<int, 6>& ref_boundary)
  {
    // we clear previous mesh
    Clear();
    
    int ref_max = ref_boundary(0);
    for (int i = 1; i < 6; i++)
      ref_max = max(ref_max, ref_boundary(i));
    
    if (ref_max > this->GetNbReferences())
      this->ResizeNbReferences(ref_max);
    
    if (this->print_level >= 2)
      cout<< " we create a structured mesh "<<endl;
    
    int nbx = 1;
    for (int i = 0; i < nb_layers_x.GetM(); i++)
      nbx += nb_layers_x(i);

    int nby = 1;
    for (int i = 0; i < nb_layers_y.GetM(); i++)
      nby += nb_layers_y(i);

    int nbz = 1;
    for (int i = 0; i < nb_layers_z.GetM(); i++)
      nbz += nb_layers_z(i);
    
    int nb_vertices = nbx*nby*nbz;
    int nb_elt = 0;
    for (int i = 0; i < nb_layers_x.GetM(); i++)
      for (int j = 0; j < nb_layers_y.GetM(); j++)
        for (int k = 0; k < nb_layers_z.GetM(); k++)
          if (ref_domain_layer(i, j, k) > 0)
            nb_elt += nb_layers_x(i)*nb_layers_y(j)*nb_layers_z(k);
    
    this->ReallocateVertices(nb_vertices);    
    this->ReallocateElements(nb_elt);

    Array3D<int> ref_domain(nbx-1, nby-1, nbz-1);    
    Vector<T> pos_x(nbx), pos_y(nby), pos_z(nbz);
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
            int  offset_iz = 0;
            for (int iz = 0; iz < nb_layers_z.GetM(); iz++)
              {
                for (int j = 0; j < nb_layers_x(i); j++)                
                  for (int jy = 0; jy < nb_layers_y(iy); jy++)                  
                    for (int jz = 0; jz < nb_layers_z(iz); jz++)
                      {
                        int ig = offset_ix + j;
                        int jg = offset_iy + jy;
                        int kg = offset_iz + jz;
                        ref_domain(ig, jg, kg) = ref_domain_layer(i, iy, iz);
                      }

                offset_iz += nb_layers_z(iz);
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

    pos_z(0) = pos_layer_z(0); nb = 1;
    for (int i = 0; i < nb_layers_z.GetM(); i++)
      {
        T x0 = pos_layer_z(i), x1 = pos_layer_z(i+1);
        for (int j = 1; j < nb_layers_z(i); j++)
          {
            T lam = T(j) / nb_layers_z(i);
            pos_z(nb) = (T(1)-lam)*x0 + lam*x1;
            nb++;
          }
        
        pos_z(nb) = x1; nb++;
      }
    
    Array3D<int> num(2, 2, 2);
    num(0, 0, 0) = 0;
    num(1, 0, 0) = 1;
    num(1, 1, 0) = nbx+1;
    num(0, 1, 0) = nbx;
    num(0, 0, 1) = nbx*nby;
    num(1, 0, 1) = nbx*nby + 1;
    num(1, 1, 1) = nbx*nby + nbx + 1;
    num(0, 1, 1) = nbx*nby + nbx;
    
    int nb_faces_ref = 2*(nbx-1)*(nby-1) + 2*(nbx-1)*(nbz-1) + 2*(nby-1)*(nbz-1);
        
    ReallocateBoundariesRef(nb_faces_ref);
    nb_elt = 0; nb_faces_ref = 0;
    for (int k = 0; k < nbz; k++)
      for (int j = 0; j < nby; j++)
	for (int i = 0; i < nbx; i++)
	  {            
	    if ((i<(nbx-1))&&(j<(nby-1))&&(k<(nbz-1)))
	      {
                if (ref_domain(i, j, k) > 0)
                  this->elements(nb_elt++).
                    InitHexahedral(num(0,0,0), num(1,0,0), num(1,1,0),
                                   num(0,1,0), num(0,0,1), num(1,0,1),
                                   num(1,1,1), num(0,1,1), ref_domain(i, j, k));
                
		// boundary x = xmin
		if (i==0)
		  {
                    faces_ref(nb_faces_ref++).
                      InitQuadrangular(num(0,0,0), num(0,1,0),
                                       num(0,1,1), num(0,0,1), ref_boundary(0));
		  }
                
		// boundary x = xmax
		if (i==(nbx-2))
		  {
                    faces_ref(nb_faces_ref++).
                      InitQuadrangular(num(1,0,0), num(1,1,0),
                                       num(1,1,1), num(1,0,1), ref_boundary(5));
		  }
		
		// boundary y = ymin
		if (j==0)
		  {
                    faces_ref(nb_faces_ref++).InitQuadrangular(num(0,0,0),num(1,0,0),num(1,0,1),
                                                               num(0,0,1),ref_boundary(1));
		  }
		
		// boundary y = ymax
		if (j==(nby-2))
		  {
                    faces_ref(nb_faces_ref++)
                      .InitQuadrangular(num(0,1,0),num(1,1,0),
                                        num(1,1,1),num(0,1,1),ref_boundary(4));
		  }
		
		// boundary z = zmin
		if (k==0)
		  {
                    faces_ref(nb_faces_ref++).InitQuadrangular(num(0,0,0),num(1,0,0),num(1,1,0),
                                                               num(0,1,0),ref_boundary(2));
		  }
		
		// boundary z = zmax
		if (k==(nbz-2))
		  {
                    faces_ref(nb_faces_ref++).InitQuadrangular(num(0,0,1),num(1,0,1),num(1,1,1),
                                                               num(0,1,1),ref_boundary(3));
		  }
	      }
	    
	    this->Vertices(num(0,0,0)).Init(pos_x(i), pos_y(j), pos_z(k));
	    
	    for (int l = 0; l < 2; l++)
	      for (int m = 0; m < 2; m++)
		for (int n = 0; n < 2; n++)
		  num(l, m, n)++;
	  }
    
    this->SortPeriodicBoundaries();
    this->SortBoundariesRef();
    
    this->ReorientElements();
    this->FindConnectivity();
    this->AddBoundaryFaces();
    this->ProjectPointsOnCurves();
    this->FindConnectivity();
    // this->Write("cube.mesh");
  }
  
  
  //! creation of a regular mesh of a ball
  /*!
    \param[in] a radius on the internal sphere
    \param[in] b radius of the external sphere
    \param[in] nb_points_by_wave_length number of subdivisions in unit interval (1/h)
    \param[in] ref_inside reference of the elements inside the ball
    \param[in] extern_sphere if false external boundary is a cube
    this methods creates a mesh of a ball of radius a, and add spheric layers between
    sphere of radius a and sphere of radius b
  */
  template<class T>
  void Mesh<Dimension3, T>
  ::CreateSphericalBall(const T& a, const T& b, const T& nb_points_by_wave_length,
                        int ref_inside, bool extern_sphere, int type_mesh, int ref_outside_domain,
                        int nb_pts_sphere, int nb_pts_line, int ref_extern_sphere)
  {    
    // modification of the mesh, so that exterior boundary is a sphere of radius a
    // estimation of the radius of the sphere, where we can put the mesh into
    T radius_mesh(0);
    for (int i = 0; i < this->GetNbVertices(); i++)
      radius_mesh = max(radius_mesh, Norm2(this->Vertices(i)));
    
    T coef = a/radius_mesh;
    for (int i = 0; i < this->GetNbVertices(); i++)
      Mlt(coef, this->Vertices(i));

    // modification of the mesh so that,
    // referenced faces on the boundary have a new reference number
    int new_ref = this->GetNewReference();
    for (int i = 0; i < GetNbBoundaryRef(); i++)
      if (faces_ref(i).GetReference() == ref_inside )
	faces_ref(i).SetReference(new_ref);

    for (int i = 0; i < GetNbEdgesRef(); i++)
      if (this->edges_ref(i).GetReference() == ref_inside )
	this->edges_ref(i).SetReference(new_ref);
    
    this->SetNewReference(new_ref, this->CURVE_SPHERE, 0);
    this->Param_curve(new_ref).Reallocate(4); this->Param_curve(new_ref).Fill(0);
    this->Param_curve(new_ref)(3) = a; this->Param_known(new_ref) = true;

    // If the mesh is tetrahedral, we cut it
    if (type_mesh == HEXAHEDRAL_MESH)
      if (GetNbTetrahedra() > 0)
	SplitIntoHexahedra();

    T step_h = this->GetMeshSize();
    
    int nb_points_div = toInteger( ceil( T(step_h*nb_points_by_wave_length) ) );
    if (nb_pts_sphere > 0)
      nb_points_div = nb_pts_sphere;
    
    IVect prime_factor; FindPrimeFactorization(nb_points_div, prime_factor, 7);
    if (nb_points_div > 1)
      {
	for (int k = 0; k < prime_factor.GetM(); k++)
	  {
	    Vector<T> step_subdiv(prime_factor(k) + 1);
	    for (int i = 0; i <= prime_factor(k); i++)
	      step_subdiv(i) = T(i)/prime_factor(k);
	    
	    SubdivideMesh(step_subdiv);
	  }
      }

    this->SetCurveType(new_ref, 0);
    
    // on revient a la reference initiale
    for (int i = 0; i < GetNbBoundaryRef(); i++)
      if (faces_ref(i).GetReference() == new_ref )
	faces_ref(i).SetReference(ref_inside);
    
    for (int i = 0; i < GetNbEdgesRef(); i++)
      if (this->edges_ref(i).GetReference() == new_ref)
	this->edges_ref(i).SetReference(ref_inside);
    
    this->FindConnectivity();
    
    this->Param_curve(ref_inside).Reallocate(4); this->Param_curve(ref_inside).Fill(0);
    this->Param_curve(ref_inside)(3) = a; 
    this->Param_known(ref_inside) = true;
    this->ProjectPointsOnCurves();

    // on extrude le maillage
    int nb_intervals = toInteger(T(ceil((b-a)*nb_points_by_wave_length))) + 1;
    if (nb_pts_line > 0)
      nb_intervals = nb_pts_line;
    
    int ref_outside = ref_extern_sphere;
    this->Param_known(ref_outside) = true;
    this->Param_curve(ref_outside).Reallocate(4); this->Param_curve(ref_outside).Fill(0);
    this->Param_curve(ref_outside)(3) = b;
    
    SurfacicMesh<Dimension3> mesh;
    mesh.SetGeometryOrder(this->GetGeometryOrder());
    IVect ref_cond(this->GetNbReferences()); ref_cond.Fill();
    this->GetBoundaryMesh(ref_inside, mesh, ref_cond);
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
	T coef = b/a;
	Mlt(coef, mesh.Vertex(i));
        if (!extern_sphere)
	  {
	    Mlt(sqrt(3), mesh.Vertex(i));
	    T x = mesh.Vertex(i)(0), y = mesh.Vertex(i)(1), z = mesh.Vertex(i)(2);
            coef = 1.0;
            if (abs(x) > b + R_N::threshold)
	      coef = min(coef, b/abs(x));
            
	    if (abs(y) > b + R_N::threshold)
	      coef = min(coef, b/abs(y));
	    
            if (abs(z) > b + R_N::threshold)
	      coef = min(coef, b/abs(z));
	    
            Mlt(coef, mesh.Vertex(i));
          }
      }

    // Write("mesh_vol.mesh");
    // mesh.Write("mesh_surf.mesh");
    IVect nb_interval(1); nb_interval(0) = nb_intervals;
    IVect ref_surf(2); ref_surf(0) = ref_inside; ref_surf(1) = ref_outside;
    IVect ref_domain(1); ref_domain(0) = ref_outside_domain; Vector<T> step_sphere;
    ExtrudeOrigin(nb_interval, mesh, ref_surf, ref_domain, step_sphere);    
  }
  
  
  //! creation of the mesh between two spheres of radius a and b
  /*!
    \param[in] a radius on the internal sphere
    \param[in] b radius of the external sphere
    \param[in] nb_points_per_unit number of subdivisions in unit interval (1/h)
    \param[in] intern_sphere if false internal boundary is a cube
    \param[in] extern_sphere if false external boundary is a cube
  */
  template<class T>
  void Mesh<Dimension3, T>
  ::CreateSphericalCrown(const T& a, const T& b, const T& nb_points_per_unit,
                         bool intern_sphere, bool extern_sphere, const Vector<T>& step_subdiv, bool all_ref)
  {
    T x = a;
    int ref_inside = 1; TinyVector<int, 6> ref_boundary; ref_boundary.Fill(1);
    Mesh<Dimension3, T> mesh_sb;
    mesh_sb.SetGeometryOrder(this->GetGeometryOrder());
    
    bool curved_intern = false;
    if (this->GetCurveType(1) == this->CURVE_SPHERE)
      {
        curved_intern = true;
        this->SetCurveType(1, this->NO_CURVE);
      }
    
    Vector<bool> curved_subdiv(step_subdiv.GetM());
    curved_subdiv.Fill(false);
    for (int k = 4; k < 4+step_subdiv.GetM(); k++)
      if (this->GetCurveType(k) == this->CURVE_SPHERE)
        {
          curved_subdiv(k-4) = true;
          this->SetCurveType(k, this->NO_CURVE);
        }

    // maillage de la sphere interne -> reference 1
    TinyVector<int, 3> nb_pts_; nb_pts_.Fill(2);
    CreateRegularMesh(R_N(-x, -x, -x), R_N(x, x, x), nb_pts_,
		      ref_inside, ref_boundary, this->HEXAHEDRAL_MESH);
    
    this->ClearElements();
    
    bool curved_extern = false;
    if (this->GetCurveType(3) == this->CURVE_SPHERE)
      {
        curved_extern = true;
        this->SetCurveType(3, this->NO_CURVE);
      }
   
    // maillage de la sphere externe -> reference 3
    x = b;  ref_boundary.Fill(3);
    int nb_offset_ref = 4;
    if (all_ref)
      {
	nb_offset_ref = 8;
	for (int i = 0; i < 6; i++)
	  ref_boundary(i) = i+2;
      }
    
    mesh_sb.CreateRegularMesh(R_N(-x, -x, -x), R_N(x, x, x), nb_pts_,
			      ref_inside, ref_boundary, this->HEXAHEDRAL_MESH);
    
    mesh_sb.ClearElements();
    
    // on subdivise le maillage
    T step_mesh = T(1)/nb_points_per_unit;
    // le nombre de points par longueur d'onde est vrai pour une section moyenne
    int nb_intervals_ = toInteger(T(ceil(sqrt(4*pi_wp)*(a+b)/(6*step_mesh))));
    Vector<T> step_subdiv_(nb_intervals_+1); step_subdiv_.Fill();
    Mlt(T(1)/nb_intervals_, step_subdiv_);
    SubdivideMesh(step_subdiv_);
    mesh_sb.SubdivideMesh(step_subdiv_);
    
    // on projette sur la sphere de rayon a, et de rayon b, si demande
    for (int i = 0; i < this->GetNbVertices(); i++)
      {
	if (intern_sphere)
	  {
	    T coef = a/Norm2(this->Vertices(i));
	    Mlt(coef, this->Vertices(i));
	  }
	
	if (extern_sphere)
	  {
	    T coef = b/Norm2(mesh_sb.Vertices(i));
	    Mlt(coef, mesh_sb.Vertices(i));
	  }
      }
    
    // numbering internal spheres 
    int nb_div = step_subdiv.GetM();
    IVect ref_surf(nb_div+2), ref_domain(nb_div+1);
    ref_surf(0) = 1; ref_surf(nb_div+1) = 3; ref_domain(nb_div) = 1;
    for (int i = 0; i < nb_div; i++)
      {
	ref_surf(i+1) = nb_offset_ref+i;
	ref_domain(i) = i+2;
      }
    
    // number of cells for each layer
    IVect nb_intervals(nb_div+1);
    if (nb_div == 0)
      nb_intervals(0) = max(toInteger((b-a)/step_mesh),1);
    else
      {
	nb_intervals(0) = max(toInteger((step_subdiv(0)-a)/step_mesh),1);
	nb_intervals(nb_div) = max(toInteger((b-step_subdiv(nb_div-1))/step_mesh),1);
	for (int i = 1; i < nb_div; i++)
	  nb_intervals(i) = max(toInteger((step_subdiv(i)-step_subdiv(i-1))/step_mesh),1);
      }
    
    // changing step_subdiv to fit in interval [0,1]
    VectReal_wp step_sphere = step_subdiv;
    for (int i = 0; i < nb_div; i++)
      step_sphere(i) = (step_sphere(i)-a)/(b-a);
    
    if (abs(T(b-a)) > R_N::threshold)
      ExtrudeOrigin(nb_intervals, mesh_sb, ref_surf, ref_domain, step_sphere, all_ref);

    Vector<T> param(4); param.Fill(0);
    if (curved_intern)
      {
        param(3) = a;
        this->SetCurveType(1, this->CURVE_SPHERE);
        this->SetCurveParameter(1, param);
      }

    if (curved_extern)
      {
        param(3) = b;
        this->SetCurveType(3, this->CURVE_SPHERE);
        this->SetCurveParameter(3, param);
      }

    for (int k = 4; k < 4+step_subdiv.GetM(); k++)
      if (curved_subdiv(k-4))
        {
          param(3) = a + step_sphere(k-4)*(b-a);
          this->SetCurveType(k, this->CURVE_SPHERE);
          this->SetCurveParameter(k, param);
        }
    
    this->ProjectPointsOnCurves();
    // Write("sphere.mesh");
    //int test_input; cout<<"On attend"<<endl; cin>>test_input;
    // Print();
    
  }
  

  //! layers are added along a cartesian direction (x, y or z)
  /*!
    \param[in] num_coord coordinate where the extrusion is done
    \param[in] nb_layers number of layers to add
    \param[in] pos position of the extern boundary to extrude
    \param[in] delta thickness of all the layers
  */
  template<class T>
  void Mesh<Dimension3, T>::ExtrudeCoordinate(int num_coord, int nb_layers,
					      const T& pos, const T& delta)
  {
    // we find all points that are on the face to be extruded
    IVect Point_to_extrude(this->GetNbVertices()); Point_to_extrude.Fill(-1);
    
    int nb_new_vertices = 0;
    for (int i = 0; i < this->GetNbVertices(); i++)
      if (abs(pos-this->Vertices(i)(num_coord)) <= R_N::threshold)
	{
	  Point_to_extrude(i) = this->GetNbVertices()+nb_new_vertices;
	  nb_new_vertices++;
	}
        
    Vector<T> subdiv(nb_layers);
    for (int j = 0; j < nb_layers; j++)
      subdiv(j) = delta*Real_wp(j+1)/nb_layers;
 
    // we add all those vertices
    int nb_old_vertices = this->GetNbVertices();
    int nb_vertices = nb_old_vertices;
    this->ResizeVertices(nb_vertices + nb_new_vertices*nb_layers);
    for (int j = 0; j < nb_layers; j++)
      for (int i = 0; i < nb_old_vertices; i++)
	if (Point_to_extrude(i) != -1)
	  {
	    this->Vertices(nb_vertices) = this->Vertices(i);
	    this->Vertices(nb_vertices)(num_coord) += subdiv(j);
	    nb_vertices++;
	  }
    
    // now, we search referenced faces that will be shifted to the new hyperplane x(num_coor) = pos
    VectBool face_on_top(GetNbBoundaryRef());
    int nb_new_faces_ref = 0, nb_new_elt = 0;
    // loop over referenced faces to find faces on the extruded boundary
    for (int i = 0; i < GetNbBoundaryRef(); i++)
      {
	face_on_top(i) = true;
	for (int j = 0; j < faces_ref(i).GetNbVertices(); j++)
	  {
	    int nv = faces_ref(i).numVertex(j);
	    if (Point_to_extrude(nv) == -1)
	      face_on_top(i) = false;
	  }
      }
    
    // counting all edges which are on the boundary of the plane to extrude
    int nb_edges_ref_on_top = 0;
    for (int i = 0; i < GetNbBoundaryRef(); i++)
      if (face_on_top(i))
	{
	  nb_new_elt++;
	  
	  // loop over edges
	  for (int j = 0; j < faces_ref(i).GetNbEdges(); j++)
	    {
	      int num_edge = faces_ref(i).numEdge(j);
	      int nb_faces_on_top = 0;
              for (int k = 0; k < this->GetEdge(num_edge).GetNbFaces(); k++)
                {
                  int num_face_ref
                    = this->GetEdge(num_edge).numFace(k);
                  
                  if (num_face_ref < this->GetNbBoundaryRef())
                    if (face_on_top(num_face_ref))
                      nb_faces_on_top++;
                }
              
	      if (nb_faces_on_top == 1)
		{
                  nb_new_faces_ref++;
                  if (num_edge < GetNbEdgesRef())
                    nb_edges_ref_on_top++;
                }
	    }
	}
    
    // creating referenced edges
    int rg = this->GetGeometryOrder();
    int nb_old_faces_ref = this->GetNbBoundaryRef();
    int nb_old_edges_ref = this->GetNbEdgesRef();
    int nb_old_edges = this->GetNbEdges();
    this->ResizeEdgesRef(nb_old_edges_ref + nb_edges_ref_on_top*nb_layers);
    this->ResizeEdges(nb_old_edges + nb_edges_ref_on_top*nb_layers);
    int nb_new_edges_ref = nb_old_edges_ref;
    int nb_new_edges = nb_old_edges; R_N ptA;
    for (int i = 0; i < nb_old_faces_ref; i++)
      {
	if (face_on_top(i))
	  {
            int nb_vert = faces_ref(i).GetNbVertices();
            // loop over edges to detect if an edge is a boundary edge
	    // in that case, extruded edges are created
	    for (int j = 0; j < nb_vert; j++)
	      {
                int num_edge = faces_ref(i).numEdge(j);
                if (num_edge < nb_old_edges_ref)
                  {
                    int nb_faces_on_top = 0;
                    for (int k = 0; k < this->GetEdge(num_edge).GetNbFaces(); k++)
                      {                        
                        int num_face_ref
                          = this->GetEdge(num_edge).numFace(k);
                        
                        if (num_face_ref < nb_old_faces_ref)
                          if (face_on_top(num_face_ref))
                            nb_faces_on_top++;
                      }
                    
                    if (nb_faces_on_top == 1)
                      {
                        int n1 = this->EdgeRef(num_edge).numVertex(0);
                        int n2 = this->EdgeRef(num_edge).numVertex(1);
                        int n1b = Point_to_extrude(n1); 
                        int n2b = Point_to_extrude(n2);
                        int ref = this->EdgeRef(num_edge).GetReference();
                        
                        for (int k = 0; k < nb_layers; k++)
                          {
                            int ne = nb_new_edges + k*nb_edges_ref_on_top;
                            int ne_ref = nb_new_edges_ref + k*nb_edges_ref_on_top;
                            this->EdgeRef(ne_ref).
                              Init(n1b+k*nb_new_vertices, n2b+k*nb_new_vertices, ref);
                            
                            this->GetEdge(ne)
                              = this->EdgeRef(nb_new_edges_ref + k*nb_edges_ref_on_top);
                            
                            for (int p = 0; p < rg-1; p++)
                              {
                                ptA = this->PointsEdgeRef(num_edge, p);
                                ptA(num_coord) += subdiv(k);
                                this->PointsEdgeRef(ne_ref, p) = ptA;
                              }                                                
                          }    
                        
                        nb_new_edges++; nb_new_edges_ref++;
                      }                
                  }
              }
          }
      }
    
    
    // reallocating arrays
    ResizeBoundariesRef(GetNbBoundaryRef()+nb_new_faces_ref*nb_layers);
    int nb_elt = this->GetNbElt();
    this->ResizeElements(nb_elt+nb_new_elt*nb_layers);
    
    nb_new_edges = nb_old_edges; 
    int nb_faces_ref = nb_old_faces_ref;
    IVect num(8), numb(8); Vector<IVect> num_face(nb_layers+1);
    VectR3 ptE(8), ptEb(8); R3 vec_u, vec_v, vec_w;
    // loop over referenced faces
    for (int i = 0; i < nb_old_faces_ref; i++)
      {
	if (face_on_top(i))
	  {
	    int nb_vert = faces_ref(i).GetNbVertices();
	    for (int k = 0; k <= nb_layers; k++)
	      num_face(k).Reallocate(nb_vert);
	   
	    // number of vertices of the original face that will be extruded 
	    for (int j = 0; j < nb_vert; j++)
	      num_face(0)(j) = faces_ref(i).numVertex(j);
	    
	    // number of the extruded vertices
	    for (int k = 0; k < nb_layers; k++)
	      for (int j = 0; j < nb_vert; j++)
		num_face(k+1)(j) = Point_to_extrude(num_face(0)(j)) + nb_new_vertices*k;
	    
	    // reference of this face 
	    int ref = faces_ref(i).GetReference();
	    // reference of the volume
	    int type = this->Element(faces_ref(i).numElement(0)).GetReference();
	    
	    // loop over edges to detect if an edge is a boundary edge
	    // in that case, lateral faces are created
	    for (int j = 0; j < nb_vert; j++)
	      {
		int num_edge = faces_ref(i).numEdge(j);
                int nb_faces_on_top = 0;
                for (int k = 0; k < this->GetEdge(num_edge).GetNbFaces(); k++)
                  {
                    int num_face_ref
                      = this->GetEdge(num_edge).numFace(k);
                    
                    if (num_face_ref < nb_old_faces_ref)
                      if (face_on_top(num_face_ref))
                        nb_faces_on_top++;
                  }
                
                if (nb_faces_on_top == 1)
		  {
		    // boundary edge
		    int jp1 = (j+1)%nb_vert;
		    // looking for another reference 
		    int ref_true = ref;
		    for (int k = 0; k < this->edges(num_edge).GetNbFaces(); k++)
		      {
                        int num_face_ref
                          = this->edges(num_edge).numFace(k);
                        
                        if (num_face_ref < nb_old_faces_ref)
                          if (!face_on_top(num_face_ref))
                            {
                              if (this->BoundaryRef(num_face_ref).GetReference() > 0)
                                ref_true = faces_ref(num_face_ref).GetReference(); 
                            }
                      }
                    
		    // new quadrilaterals
		    for (int k = 0; k < nb_layers; k++)
		      {
			faces_ref(nb_faces_ref)
                          .InitQuadrangular(num_face(k)(j), num_face(k)(jp1),
                                            num_face(k+1)(jp1), num_face(k+1)(j), ref_true);
                        
                        if (num_edge < nb_old_edges_ref)
                          {
                            if (k == 0)
                              faces_ref(nb_faces_ref).SetEdge(0, num_edge);
                            else
                              faces_ref(nb_faces_ref)
                                .SetEdge(0, nb_new_edges+(k-1)*nb_edges_ref_on_top);
                            
                            faces_ref(nb_faces_ref)
                              .SetEdge(2, nb_new_edges+k*nb_edges_ref_on_top);    
                          }
                        
                        this->ComputeInternalPointsFaceWithCurvedEdges(nb_faces_ref);
                        
                        nb_faces_ref++;
		      }
                    
                    if (num_edge < nb_old_edges_ref)
                      nb_new_edges++;
		  }
	      }
		
	    // new elements are created (hexa or wedge)
	    num.Reallocate(2*nb_vert);
	    for (int k = 0; k < nb_layers; k++)
	      {
		for (int j = 0; j < nb_vert; j++)
		  {
		    num(j) = num_face(k)(j);
                    ptE(j) = this->Vertex(num(j));
		    num(j+nb_vert) = num_face(k+1)(j);
                    ptE(j+nb_vert) = this->Vertex(num(j+nb_vert));
		  }
                
		if (nb_vert == 4)
                  {
                    // element is reoriented in order to fit to the cartesian coordinates
                    Mesh<Dimension3>::ReorientElementCartesian(num, ptE);
                    
                    this->elements(nb_elt).InitHexahedral(num(0), num(1), num(2), num(3),
                                                          num(4), num(5), num(6), num(7), type);
                  }
                else
                  this->elements(nb_elt).InitWedge(num(0), num(1), num(2),
                                                   num(3), num(4), num(5), type);
                
		nb_elt++;
	      }
	    
	    // referenced face is changed in order to have top face referenced
	    faces_ref(i).Init(num_face(nb_layers), ref);
            for (int k = 0; k < this->PointsFaceRef(i).GetM(); k++)
              this->PointsFaceRef(i)(k)(num_coord) += delta;
	    
	  }
      }
    
    this->SortEdgesRef();
    this->SortBoundariesRef();
    this->ReorientElements();
    this->FindConnectivity();
    this->ProjectPointsOnCurves();
  }

  
  //! a mesh is appended to the current mesh
  /*!
    \param[in] mesh the mesh to add
    \param[in] elimine if true, this method removes duplicate vertices
  */
  template<class T>
  void Mesh<Dimension3, T>::AppendMesh(const Mesh<Dimension3, T>& mesh, bool elimine)
  {
    if (mesh.GetNbVertices() <= 0)
      return;
    
    if (mesh.GetNbReferences() > this->GetNbReferences())
      this->ResizeNbReferences(mesh.GetNbReferences());
    
    for (int ref = 1; ref < mesh.GetNbReferences(); ref++)
      {
        int type1 = mesh.GetCurveType(ref);
        int type2 = this->GetCurveType(ref);
        if (type2 == this->NO_CURVE)
          {
            this->SetCurveType(ref, type1);
            this->SetCurveParameter(ref, mesh.GetCurveParameter(ref));
          }
        else if (type2 == this->CURVE_FILE)
          {
          }
        else
          {
            if ((type1 != mesh.NO_CURVE)&&(type1 != mesh.CURVE_FILE))
              {
                if (type1 != type2)
                  {
                    cout<<"curves different "<<endl;
                    cout<<"references should be changed"<<endl;
                    abort();
                  }
              }
          }
      }

    Matrix<Real_wp> ValPhi, ValPhiTri, ValPhiQuad;
    // interpolation necessary if the order of two meshes are different
    int r = this->GetGeometryOrder();
    int rm = mesh.GetGeometryOrder();
    if (r != rm)
      {
	ValPhi.Reallocate(rm+1, r-1);
	for (int i = 0; i < ValPhi.GetM(); i++)
	  for (int j = 0; j < ValPhi.GetN(); j++)
	    ValPhi(i, j) = mesh.lob_curve.EvaluatePhi(i, this->lob_curve.Points(j+1));
	
        VectReal_wp phi;
	ValPhiTri.Reallocate((rm+1)*(rm+2)/2, (r-1)*(r-2)/2);
	for (int j = 0; j < ValPhiTri.GetN(); j++)
	  {
            mesh.triangle_reference
              .ComputeValuesPhiNodalRef(this->triangle_reference.PointsNodalND(3*r+j), phi);
            
            for (int i = 0; i < ValPhiTri.GetM(); i++)
              ValPhiTri(i, j) = phi(i);
          }

	ValPhiQuad.Reallocate((rm+1)*(rm+1), (r-1)*(r-1));
	for (int j = 0; j < ValPhiQuad.GetN(); j++)
	  {
            mesh.quadrangle_reference
              .ComputeValuesPhiNodalRef(this->quadrangle_reference.PointsNodalND(4*r+j), phi);
            
            for (int i = 0; i < ValPhiQuad.GetM(); i++)
              ValPhiQuad(i, j) = phi(i);
          }
      }

    // vertices are merged
    int nb_old = this->GetNbVertices(); int offset = nb_old;
    int nb_vertices = nb_old + mesh.GetNbVertices();
    this->ResizeVertices(nb_vertices);
    for (int i = nb_old; i < nb_vertices; i++)
      this->Vertices(i) = mesh.Vertices(i-nb_old);
    
    // edges are merged
    nb_old = this->GetNbEdgesRef();
    int nb_edges_ref = nb_old + mesh.GetNbEdgesRef();
    this->ResizeEdgesRef(nb_edges_ref);
    for (int i = nb_old; i < nb_edges_ref; i++)
      {
        this->EdgeRef(i) = mesh.EdgeRef(i-nb_old);
        this->EdgeRef(i).Init(mesh.EdgeRef(i-nb_old).numVertex(0)+offset,
                              mesh.EdgeRef(i-nb_old).numVertex(1)+offset);
        
        if (rm != r)
	  {
	    // interpolation
	    VectR3 OldPoints(rm+1);
	    OldPoints(0) = mesh.Vertex(mesh.EdgeRef(i-nb_old).numVertex(0));
	    OldPoints(rm) = mesh.Vertex(mesh.EdgeRef(i-nb_old).numVertex(1));
	    for (int k = 1; k < rm; k++)
	      OldPoints(k) = mesh.PointsEdgeRef(i-nb_old, k-1);
	    
	    for (int k = 0; k < r-1; k++)
	      {
		this->PointsEdgeRef(i, k).Zero();
		for (int j = 0; j < OldPoints.GetM(); j++)
		  Add(ValPhi(j, k), OldPoints(j), this->PointsEdgeRef(i, k));
	      }
	  }
	else
	  for (int k = 0; k < r-1; k++)
	    this->PointsEdgeRef(i, k) = mesh.PointsEdgeRef(i-nb_old, k);
        
      }
    
    // faces are merged
    int nb_faces_ref = this->GetNbBoundaryRef();
    nb_old = nb_faces_ref;
    nb_faces_ref += mesh.GetNbBoundaryRef();
    ResizeBoundariesRef(nb_faces_ref);
    for (int i = nb_old; i < nb_faces_ref; i++)
      {
	this->BoundaryRef(i) = mesh.BoundaryRef(i-nb_old);
        int ref = this->BoundaryRef(i).GetReference();
        int nb_vert = this->BoundaryRef(i).GetNbVertices();
        IVect num(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          num(j) = this->BoundaryRef(i).numVertex(j) + offset;
        
        this->BoundaryRef(i).Init(num, ref);
        
        if (rm != r)
	  {
	    // interpolation
            if (nb_vert == 3)
              {
                VectR3 OldPoints;
                mesh.GetNodalPointsOnFace(i-nb_old, OldPoints);
                
                this->PointsFaceRef(i).Reallocate(ValPhiTri.GetN());
                for (int k = 0; k < ValPhiTri.GetN(); k++)
                  {
                    this->PointsFaceRef(i)(k).Zero();
                    for (int j = 0; j < OldPoints.GetM(); j++)
                      Add(ValPhiTri(j, k), OldPoints(j), this->PointsFaceRef(i)(k));
                  }
              }
            else
              {
                VectR3 OldPoints;
                mesh.GetNodalPointsOnFace(i-nb_old, OldPoints);
                                
                this->PointsFaceRef(i).Reallocate(ValPhiQuad.GetN());
                for (int k = 0; k < ValPhiQuad.GetN(); k++)
                  {
                    this->PointsFaceRef(i)(k).Zero();
                    for (int j = 0; j < OldPoints.GetM(); j++)
                      Add(ValPhiQuad(j, k), OldPoints(j), this->PointsFaceRef(i)(k));
                  }
              }
	  }
	else
          this->PointsFaceRef(i) = mesh.PointsFaceRef(i-nb_old);
        
      }
    
    // rajouts des elements volumiques
    nb_old = this->GetNbElt();
    int nb_elt = nb_old + mesh.GetNbElt();
    this->ResizeElements(nb_elt);
    for (int i = nb_old; i < nb_elt; i++)
      {
	this->Element(i) = mesh.Element(i-nb_old);
        int ref = this->Element(i).GetReference();
        int nb_vert = this->Element(i).GetNbVertices();
        IVect num(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          num(j) = this->Element(i).numVertex(j) + offset;
        
        this->Element(i).Init(num, ref);
      }
    
    this->ReorientElements();

    // on elimine tous les sommets doublons
    if (elimine)
      this->RemoveDuplicateVertices();
    else
      {
        this->FindConnectivity();
        this->ProjectPointsOnCurves();
      }
  }
  
    
  //! adds spherical layers to the mesh
  /*!
    Several layers are added to the extern boundary of a sphere
    \param[in] mesh external boundary of the mesh
    \param[in] nb_intervals number of cells for each layer
    \param[in] ref_surf reference for each sphere
    \param[in] ref_domain reference of elements in layers
    \param[in] step_subdiv positions of internal spheres
  */
  template<class T>
  void Mesh<Dimension3, T>
  ::ExtrudeOrigin(const IVect& nb_intervals, const Mesh<Dimension3>& mesh,
                  const IVect& ref_surf, const IVect& ref_domain, const Vector<T>& step_subdiv, bool all_ref)
  {
    if (mesh.GetGeometryOrder() != this->GetGeometryOrder())
      {
        cout << " case not treated " <<endl;
        abort();
      }
    
    int r = this->GetGeometryOrder();
    // extraction of the extern boundary of the initial mesh
    SurfacicMesh<Dimension3> mesh_surfacic; 
    IVect Ref_identity(this->GetNbReferences()+1); Ref_identity.Fill();
    this->GetBoundaryMesh(ref_surf(0), mesh_surfacic, Ref_identity);
    
    // increasing the number of vertices 
    int nb_vertices = this->GetNbVertices();
    int nb_old_vertices = this->GetNbVertices();
    int nb_div = step_subdiv.GetM();
    for (int i = 0; i <= nb_div; i++)
      nb_vertices += mesh_surfacic.GetNbVertices()*nb_intervals(i);
    
    this->ResizeVertices(nb_vertices);
        
    R3 pointA, pointB, pointC, pointD, point;
    
    nb_vertices = nb_old_vertices;
    // creation of vertices for all the layers
    for (int k = 0; k <= nb_div; k++)
      for (int j = 0; j < nb_intervals(k); j++)
	for (int i = 0; i < mesh_surfacic.GetNbVertices(); i++)
	  {
	    // barycenter between internal boundary and external boundary
	    pointC = mesh_surfacic.Vertex(i); 
	    pointD = mesh.Vertex(i);
	    T li(0), lip1(0);
	    if (k == 0)
	      li = 0;
	    else
	      li = step_subdiv(k-1);
	    
	    if (k==nb_div)
	      lip1 = 1;
	    else
	      lip1 = step_subdiv(k);
	    
	    // pointA and pointB are on the interface between layers
	    pointA.Zero(); Add(li, pointD, pointA); Add(1.0-li, pointC, pointA);
	    pointB.Zero(); Add(lip1, pointD, pointB); Add(1.0-lip1, pointC, pointB);
	    
	    point = pointA; Mlt(Real_wp(nb_intervals(k)-j-1)/nb_intervals(k),point);
	    Add(Real_wp(j+1)/nb_intervals(k),pointB,point);
	    
	    this->Vertices(nb_vertices++) = point;
	  }
    
    int offset = nb_old_vertices-mesh_surfacic.GetNbVertices();
    // offset_cumul(k)(j) is the number of vertices of the mesh before adding sublayer j of layer k
    Vector<IVect> offset_cumul(nb_div+1);
    for (int k = 0; k <= nb_div; k++)
      {
	offset_cumul(k).Reallocate(nb_intervals(k)+1);
	offset_cumul(k)(0) = offset;
	
	for (int j = 0; j < nb_intervals(k); j++)
	  offset_cumul(k)(j+1) = offset_cumul(k)(j) + mesh_surfacic.GetNbVertices();
	
	offset = offset_cumul(k)(nb_intervals(k));
      }
    
    // creation of referenced edges
    int nb_old_edges_ref = this->GetNbEdgesRef();
    int nb_edges_ref = nb_old_edges_ref + mesh.GetNbEdgesRef()*(nb_div+1);
    this->ResizeEdgesRef(nb_edges_ref);
    nb_edges_ref = nb_old_edges_ref;
    for (int k = 0; k <= nb_div; k++)
      {
        T lip1(0);
        if (k==nb_div)
          lip1 = 1;
        else
          lip1 = step_subdiv(k);
        
        for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
          {
            int n0 = mesh.EdgeRef(i).numVertex(0);
            int n1 = mesh.EdgeRef(i).numVertex(1);
            int ref = mesh.EdgeRef(i).GetReference();
            if (ref > 0)
              if (mesh.GetCurveType(ref) == mesh.CURVE_FILE)
                this->SetCurveType(ref_surf(k+1), this->CURVE_FILE);
            
            offset = offset_cumul(k)(nb_intervals(k));
            this->EdgeRef(nb_edges_ref).Init(offset+n0, offset+n1, ref_surf(k+1));
            for (int k = 0; k < r-1; k++)
              {	    
                pointC = mesh_surfacic.PointsEdgeRef(i, k);
                pointD = mesh.PointsEdgeRef(i, k);
                
                // position of the interface
                pointB.Zero(); Add(lip1, pointD, pointB); Add(1.0-lip1, pointC, pointB);
                this->PointsEdgeRef(nb_edges_ref, k) = pointB;
              }
            
            nb_edges_ref++;
          }
      }
    
    // creation of referenced faces
    int nb_old_faces_ref = this->GetNbBoundaryRef();
    int nb_faces_ref = nb_old_faces_ref + mesh.GetNbBoundaryRef()*(nb_div+1);
    this->ResizeBoundariesRef(nb_faces_ref);
    nb_faces_ref = nb_old_faces_ref;
    for (int k = 0; k <= nb_div; k++)
      {
        T lip1(0);
        if (k==nb_div)
          lip1 = 1;
        else
          lip1 = step_subdiv(k);
        
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
            IVect num(nb_vert);
            offset = offset_cumul(k)(nb_intervals(k));
            for (int j = 0; j < nb_vert; j++)
              num(j) = offset + mesh.BoundaryRef(i).numVertex(j);
            
            int ref = mesh.BoundaryRef(i).GetReference();
            if (ref > 0)
              if (mesh.GetCurveType(ref) == mesh.CURVE_FILE)
                this->SetCurveType(ref_surf(k+1), this->CURVE_FILE);
            
            if ((all_ref) && (k == nb_div))
	      this->BoundaryRef(nb_faces_ref).Init(num, ref);
	    else
	      this->BoundaryRef(nb_faces_ref).Init(num, ref_surf(k+1));
	    
            int Nquad = mesh.PointsFaceRef(i).GetM();
            this->PointsFaceRef(nb_faces_ref).Reallocate(Nquad);
            for (int k = 0; k < Nquad; k++)
              {	    
                pointC = mesh_surfacic.PointsFaceRef(i)(k);
                pointD = mesh.PointsFaceRef(i)(k);
                
                // position of the interface
                pointB.Zero(); Add(lip1, pointD, pointB); Add(1.0-lip1, pointC, pointB);
                this->PointsFaceRef(nb_faces_ref)(k) = pointB;
              }
            
            nb_faces_ref++;
          }
      }
    
    
    // creation des hexahedres/prismes		       
    // increasing the number of elements
    int nb_elt = this->GetNbElt();
    int nb_old_elt = this->GetNbElt();
    for (int i = 0; i <= nb_div; i++)
      nb_elt += mesh_surfacic.GetNbBoundaryRef()*nb_intervals(i);
    
    if (nb_old_elt == 0)
      this->ReallocateElements(nb_elt);
    else
      this->ResizeElements(nb_elt);

    nb_elt = nb_old_elt;
    for (int k = 0; k <= nb_div; k++)
      {
	for (int i = 0; i < mesh_surfacic.GetNbBoundaryRef(); i++)
	  {
	    int nv = mesh_surfacic.BoundaryRef(i).GetNbVertices();
	    // vertex numbers of bottom and top
	    IVect num_bottom(nv), num_top(nv), nums(nv);
	    // num_elt -> vertex numbers of each element that will be created (hexahedron or wedge)
	    IVect num_elt(2*nv);
	    // nums -> vertex numbers of the initial face
	    offset = offset_cumul(k)(1);
	    for (int j = 0; j < nv; j++)
	      {
		nums(j) = mesh_surfacic.BoundaryRef(i).numVertex(j); 
		if (k == 0)
                  num_bottom(j) = mesh_surfacic.ListeVertices(nums(j));
                else
                  num_bottom(j) = offset_cumul(k)(0) + nums(j);
                
		num_top(j) = offset + nums(j);
	      }
	    
	    // we create elements
	    for (int j = 0; j < nb_intervals(k); j++)
	      {
		for (int m = 0; m < nv; m++)
		  {
		    num_elt(m) = num_bottom(m);
		    num_elt(m+nv) = num_top(m);
		  }
		
		this->elements(nb_elt++).Init(num_elt, ref_domain(k));
		
		// updating num_bottom and num_top
		if (j < (nb_intervals(k)-1))
		  for (int m = 0; m < nv; m++)
		    {
		      num_bottom(m) = num_top(m);
		      num_top(m) = nums(m) + offset_cumul(k)(j+2);
		    }
	      }	    
	  }
      }
    
    this->ReorientElements();
    this->FindConnectivity();
    this->ProjectPointsOnCurves();

    // we create referenced faces (if needed)
    this->AddBoundaryFaces();
    this->CreateConnexionEdgeToFace();
    this->ConstructCrossReferenceBoundaryRef();
  }

  
  //! extrusion of a 2-D mesh along a vector, with possible holes or dielectric objects
  /*!
    \param[in] mesh_2d 2-D mesh to extrude
    \param[in] nbCells_layer number of cells to be created along axis of extrusion, for each layer
    \param[in] step_layer x-coordinate for each interface
    \param[in] ref_surfacic reference of extruded surfacic shapes
    ref_surfacic(2*i) -> new references of lateral faces of layer i
    ref_surfacic(2*i+1) -> new references of top faces of layer i
    \param[in] ref_volumic reference of extruded volumic shapes
    \param[in] vec_e1 first vector of plate
    \param[in] vec_e2 second vector of plate
    \param[in] vec_e3 axis of extrusion
  */
  template<class T>
  void Mesh<Dimension3, T>::
  ExtrudeSurfaceMesh(Mesh<Dimension2, T>& mesh_2d, const IVect& nbCells_layer,
		     const Vector<T>& step_layer, const Vector<IVect>& ref_surfacic,
                     const Vector<IVect>& ref_volumic,
		     const R3& vec_e1, const R3& vec_e2, const R3& vec_e3)
  {
    // previous mesh is cleared
    Clear();
    
    // we get the maximum reference number for the surfaces
    int ref_max = mesh_2d.GetNbReferences();
    for (int i = 0; i < ref_surfacic.GetM(); i++)
      for (int j = 0; j < ref_surfacic(i).GetM(); j++)
	ref_max = max(ref_max, ref_surfacic(i)(j));
    
    if (this->GetNbReferences() < ref_max)
      this->ResizeNbReferences(ref_max);
    
    // counting layers
    int nb_layers = nbCells_layer.GetM(); R3 point;
    int nb_all_layers = 0;
    for (int i = 0; i < nb_layers; i++)
      nb_all_layers += nbCells_layer(i);
    
    // we create vertices
    int nb_vert2d = mesh_2d.GetNbVertices();
    int nb_vertices = nb_vert2d*(nb_all_layers+1);
    this->ReallocateVertices(nb_vertices);
    for (int i = 0; i < nb_vert2d; i++)
      {
	point.Zero(); 
	Add(mesh_2d.Vertex(i)(0), vec_e1, point);
	Add(mesh_2d.Vertex(i)(1), vec_e2, point);
	Add(step_layer(0), vec_e3, point);
	this->Vertices(i) = point;
      }
    
    int offset = nb_vert2d;     T coef, lambda;
    for (int num_layer = 0; num_layer < nb_layers; num_layer++)
      {
	// creation of vertices of the layer
	for (int j = 0; j < nbCells_layer(num_layer); j++)
	  {
	    // translation of new vertices
	    lambda = T(j+1)/nbCells_layer(num_layer);
	    coef = (1.0-lambda)*step_layer(num_layer) + lambda*step_layer(num_layer+1);
	    for (int i = 0; i < nb_vert2d; i++)
	      {
		point.Zero(); 
		Add(mesh_2d.Vertex(i)(0), vec_e1, point);
		Add(mesh_2d.Vertex(i)(1), vec_e2, point);
		Add(coef, vec_e3, point);
		this->Vertices(offset + i) = point;
	      }
            
            offset += nb_vert2d;
	  }
      }
    
    int nb_elt2d = mesh_2d.GetNbElt();
    int nb_edge2d = mesh_2d.GetNbBoundaryRef();
    
    int order = this->GetGeometryOrder();
    int rm = mesh_2d.GetGeometryOrder();
    Matrix<Real_wp> ValPhi; R2 pt2d; VectR2 OldPoints(rm+1);
    if (rm != order)
      {
        ValPhi.Reallocate(rm+1, order-1);
        for (int i = 0; i <= rm; i++)
          for (int j = 0; j < order-1; j++)
            ValPhi(i, j) = mesh_2d.GetInterpolate1D(i, this->GetNodalPointInsideEdge(j));
      }
    
    // creating referenced edges
    int nb_edges_ref = nb_edge2d*(nb_all_layers+1);
    this->ReallocateEdgesRef(nb_edges_ref);
    nb_edges_ref = 0;
    for (int i = 0; i < nb_edge2d; i++)
      {
        int ref = mesh_2d.BoundaryRef(i).GetReference();
        if (mesh_2d.GetCurveType(ref) > 0)
          {
            this->SetCurveType(ref, this->CURVE_FILE);
            int n0 = mesh_2d.BoundaryRef(i).numVertex(0);
            int n1 = mesh_2d.BoundaryRef(i).numVertex(1);
            this->EdgeRef(nb_edges_ref).Init(n0, n1, ref);
            
            if (rm != order)
              {
                OldPoints(0) = mesh_2d.Vertex(n0);
                OldPoints(rm) = mesh_2d.Vertex(n1);
                for (int k = 0; k < rm-1; k++)
                  OldPoints(k+1) = mesh_2d.GetPointInsideEdge(i, k);
                
                for (int k = 0; k < order-1; k++)
                  {
                    pt2d.Zero();
                    for (int j = 0; j <= rm; j++)
                      Add(ValPhi(j, k), OldPoints(j), pt2d);
                    
                    point.Zero();
                    Add(pt2d(0), vec_e1, point);
                    Add(pt2d(1), vec_e2, point);
                    this->PointsEdgeRef(nb_edges_ref, k) = point;
                  }
              }
            else
              for (int k = 0; k < order-1; k++)
                {
                  point.Zero();
                  Add(mesh_2d.GetPointInsideEdge(i, k)(0), vec_e1, point);
                  Add(mesh_2d.GetPointInsideEdge(i, k)(1), vec_e2, point);
                  this->PointsEdgeRef(nb_edges_ref, k) = point;
                }
            nb_edges_ref++;
          }
      }
    
    offset = nb_vert2d;
    for (int num_layer = 0; num_layer < nb_layers; num_layer++)
      {
        for (int j = 0; j < nbCells_layer(num_layer); j++)
          {
            lambda = T(j+1)/nbCells_layer(num_layer);
	    coef = (1.0-lambda)*step_layer(num_layer) + lambda*step_layer(num_layer+1);
            for (int i = 0; i < nb_edge2d; i++)
              {
                int ref = mesh_2d.BoundaryRef(i).GetReference();
                if (mesh_2d.GetCurveType(ref) > 0)
                  {
                    int ref_edge = ref_surfacic(2*num_layer)(ref);
                    if (ref_edge <= 0)
                      ref_edge = ref_surfacic(2*num_layer+1)(ref);
                
                    if (ref_edge > 0)
                      {
                        this->SetCurveType(ref_edge, this->CURVE_FILE);
                        int n0 = mesh_2d.BoundaryRef(i).numVertex(0);
                        int n1 = mesh_2d.BoundaryRef(i).numVertex(1);
                        this->EdgeRef(nb_edges_ref).Init(offset+n0, offset+n1, ref_edge);
                        
                        // computation of internal nodes of the edge
                        if (rm != order)
                          {
                            OldPoints(0) = mesh_2d.Vertex(n0);
                            OldPoints(rm) = mesh_2d.Vertex(n1);
                            for (int k = 0; k < rm-1; k++)
                              OldPoints(k+1) = mesh_2d.GetPointInsideEdge(i, k);
                            
                            for (int k = 0; k < order-1; k++)
                              {
                                pt2d.Zero();
                                for (int j = 0; j <= rm; j++)
                                  Add(ValPhi(j, k), OldPoints(j), pt2d);
                                
                                point.Zero();
                                Add(pt2d(0), vec_e1, point);
                                Add(pt2d(1), vec_e2, point);
                                Add(coef, vec_e3, point);
                                this->PointsEdgeRef(nb_edges_ref, k) = point;
                              }
                          }
                        else
                          for (int k = 0; k < order-1; k++)
                            {
                              point.Zero();
                              Add(mesh_2d.GetPointInsideEdge(i, k)(0), vec_e1, point);
                              Add(mesh_2d.GetPointInsideEdge(i, k)(1), vec_e2, point);
                              Add(coef, vec_e3, point);
                              this->PointsEdgeRef(nb_edges_ref, k) = point;
                            }
                        nb_edges_ref++;
                      }
                  }
              }
            
            offset += nb_vert2d;
          }
      }
    
    this->ResizeEdgesRef(nb_edges_ref);

    int nb_faces_ref = nb_elt2d*(nb_layers+1) + nb_edge2d*nb_all_layers;
    // referenced faces on basis of the extruded object
    // we keep the references of the 2-D mesh
    ReallocateBoundariesRef(nb_faces_ref); IVect num;
    for (int i = 0; i < nb_elt2d; i++)
      {
	int nb_vert = mesh_2d.Element(i).GetNbVertices();
	int ref = mesh_2d.Element(i).GetReference();
        num.Reallocate(nb_vert);
	for (int j = 0; j < nb_vert; j++)
	  num(j) = mesh_2d.Element(i).numVertex(j);
	
	faces_ref(i).Init(num, ref);
      }
    
    // loop on layers to create referenced faces
    nb_faces_ref = nb_elt2d; offset = 0;
    for (int num_layer = 0; num_layer < nb_layers; num_layer++)
      {
	// lateral faces with reference : ref_surfaces(2*num_layer)
	for (int i = 0; i < nb_edge2d; i++)
	  {
	    int ref = mesh_2d.BoundaryRef(i).GetReference();
	    int n0 = mesh_2d.BoundaryRef(i).numVertex(0);
	    int n1 = mesh_2d.BoundaryRef(i).numVertex(1);
	    if (ref_surfacic(2*num_layer)(ref) > 0)
	      for (int j = 0; j < nbCells_layer(num_layer); j++)
		{
		  num.Reallocate(4);
		  num(0) = offset + j*nb_vert2d + n0;
		  num(1) = offset + j*nb_vert2d + n1;
		  num(2) = offset + (j+1)*nb_vert2d + n1;
		  num(3) = offset + (j+1)*nb_vert2d + n0;
		  faces_ref(nb_faces_ref).Init(num, ref_surfacic(2*num_layer)(ref));
		  nb_faces_ref++;
		}
	  }
        
        offset += nb_vert2d*nbCells_layer(num_layer);
	
	// faces of the top section (reference : ref_surfacic(2*num_layer+1)
	for (int i = 0; i < nb_elt2d; i++)
	  {
	    int nb_vert = mesh_2d.Element(i).GetNbVertices();
	    int ref = mesh_2d.Element(i).GetReference(); num.Reallocate(nb_vert);
	    for (int j = 0; j < nb_vert; j++)
	      num(j) = offset + mesh_2d.Element(i).numVertex(j);
	    
	    if (ref_surfacic(2*num_layer+1)(ref) > 0)
	      {
		faces_ref(nb_faces_ref).Init(num, ref_surfacic(2*num_layer+1)(ref));
		nb_faces_ref++;
	      }
	  }
      }
    
    this->ResizeBoundariesRef(nb_faces_ref);
    
    // now creating volume elements
    int nb_elt = nb_elt2d*nb_all_layers;
    this->ReallocateElements(nb_elt);
    nb_elt = 0; offset = 0;
    for (int num_layer = 0; num_layer < nb_layers; num_layer++)
      {
        for (int i = 0; i < nb_elt2d; i++)
	  {
	    int ref = mesh_2d.Element(i).GetReference();
	    if (ref_volumic(num_layer)(ref) > 0)
	      {
		int nb_vert = mesh_2d.Element(i).GetNbVertices();
		num.Reallocate(2*nb_vert);
		for (int j = 0; j < nb_vert; j++)
		  {
		    num(j) = offset + mesh_2d.Element(i).numVertex(j);
		    num(nb_vert+j) = num(j) + nb_vert2d;
		  }
		
		for (int j = 0; j < nbCells_layer(num_layer); j++)
		  {		    
		    this->Element(nb_elt).Init(num, ref_volumic(num_layer)(ref));
		    for (int k = 0; k < 2*nb_vert; k++)
		      num(k) += nb_vert2d;
		    
		    nb_elt++;
		  }
	      }
	  }
        
        offset += nb_vert2d*nbCells_layer(num_layer);
      }
    
    this->ResizeElements(nb_elt);
    this->ReorientElements();
    // removing unused vertices
    ForceCoherenceMesh();    
    
    this->SortPeriodicBoundaries();
    this->SortBoundariesRef();
    this->FindConnectivity();
    
    
    // replacing PointsFaceRef
    for (int i = 0; i < this->GetNbBoundaryRef(); i++)
      this->ComputeInternalPointsFaceWithCurvedEdges(i);
    
    if (this->print_level >= 2)
      Write("fluide.mesh");
    
    // exit(0);
  }
  
  
  void FctScalingExtrudeMesh::EvaluateCoefficient(const Real_wp& z,
						  Real_wp& coef_x, Real_wp& coef_y)
  {
    Real_wp lambda = 2.0*(z - 0.5*(z0+z1))/(z1-z0);
    coef_x = Real_wp(1);
    coef_y = Real_wp(1);
    if (abs(lambda) <= 1)
      {
	coef_x = 1.0 + (max_ratio-1.0)*square(cos(0.5*pi_wp*lambda));
	coef_y = 1.0/coef_x;
      }
  }

  
  //! general construction of a mesh
  /*!
    \param[in] type_mesh (hybrid, purely hexahedral or purely tetrahedral)
    \param[in] vars object containing problem datas
    \param[in] parameters matching line of the data file
    this method can read or generate mesh, depending what the user wants
    the mesh is then split, depending the type. Boundaries are computed,
    and if asked, PML layers are added
    In the data file, it refers to the line :
    FileMesh = parameters
  */
  template<class T>
  void Mesh<Dimension3, T>::ConstructMesh(int type_mesh, const VectString& parameters)
  {
    if (this->required_mesh_type >= 0)
      type_mesh = this->required_mesh_type;
    
    if (parameters.GetM() <= 0)
      return;
    
    // for predefined shapes, parameters(0) = REGULAR, SPHERE, SPHERE_CUBE,
    // BALL, BALL_CUBE or EXTRUDE
    // for any mesh parameters(0) = name of the mesh file
    if (!parameters(0).compare("REGULAR"))
      {
	// regular mesh of a paralleliped
	TinyVector<int, 3> nbPoints;
	int nbx = to_num<int>(parameters(1));
	if (nbx <= 1)
	  nbx = 2;
	
	nbPoints.Fill(nbx);

	R3 ptMin, ptMax; ptMax.Fill(1);
	if (parameters.GetM() >= 8)
	  {
	    to_num(parameters(2), ptMin(0)); to_num(parameters(3), ptMax(0));
	    to_num(parameters(4), ptMin(1)); to_num(parameters(5), ptMax(1));
	    to_num(parameters(6), ptMin(2)); to_num(parameters(7), ptMax(2));
	  }
	  	
	int nby = toInteger( T( (ptMax(1)-ptMin(1))/(ptMax(0)-ptMin(0))*nbx ) );
	int nbz = toInteger( T( (ptMax(2)-ptMin(2))/(ptMax(0)-ptMin(0))*nbx ) );
	if (nby <= 1)
	  nby = 2;
	
	if (nbz <= 1)
	  nbz = 2;
	
	nbPoints(1) = nby; nbPoints(2) = nbz;
	int ref_domain(1);
	TinyVector<int, 6> ref_boundary; ref_boundary.Fill(3);
	if (parameters.GetM() >= 10)
	  {
	    to_num(parameters(8), ref_domain);
	    to_num(parameters(9), ref_boundary(0));
	    if (parameters.GetM() >= 15)
	      {
		to_num(parameters(10), ref_boundary(1));
		to_num(parameters(11), ref_boundary(2));
		to_num(parameters(12), ref_boundary(3));
		to_num(parameters(13), ref_boundary(4));
		to_num(parameters(14), ref_boundary(5));
	      }
	    else
	      ref_boundary.Fill(ref_boundary(0));
	  }
	
	CreateRegularMesh(ptMin, ptMax, nbPoints, 1, ref_boundary, type_mesh);
      }
    else if (!parameters(0).compare("REGULAR_ANISO"))
      {
	// regular mesh of a paralleliped
	TinyVector<int, 3> nbPoints;
	int nbx = to_num<int>(parameters(1));
	int nby = to_num<int>(parameters(2));
	int nbz = to_num<int>(parameters(3));	
	if (nbx <= 1)
	  nbx = 2;
	
	if (nby <= 1)
	  nby = 2;
	
	if (nbz <= 1)
	  nbz = 2;
	
	nbPoints(0) = nbx; nbPoints(1) = nby; nbPoints(2) = nbz;

	R3 ptMin, ptMax; ptMax.Fill(1);
	if (parameters.GetM() >= 10)
	  {
	    to_num(parameters(4), ptMin(0)); to_num(parameters(5), ptMax(0));
	    to_num(parameters(6), ptMin(1)); to_num(parameters(7), ptMax(1));
	    to_num(parameters(8), ptMin(2)); to_num(parameters(9), ptMax(2));
	  }
	
	int ref_domain(1);
	TinyVector<int, 6> ref_boundary; ref_boundary.Fill(3);
	if (parameters.GetM() >= 12)
	  {
	    to_num(parameters(10), ref_domain);
	    to_num(parameters(11), ref_boundary(0));
	    if (parameters.GetM() >= 17)
	      {
		to_num(parameters(12), ref_boundary(1));
		to_num(parameters(13), ref_boundary(2));
		to_num(parameters(14), ref_boundary(3));
		to_num(parameters(15), ref_boundary(4));
		to_num(parameters(16), ref_boundary(5));
	      }
	    else
	      ref_boundary.Fill(ref_boundary(0));
	  }
	
	CreateRegularMesh(ptMin, ptMax, nbPoints, 1, ref_boundary, type_mesh);
      }
    else if (!parameters(0).compare("SPHERE"))
      {
	// regular mesh of a spherical crown
	T radius_sphere_inside = to_num<T>(parameters(1).data());
	T radius_sphere_outside = to_num<T>(parameters(2).data());
	T nb_points_by_wavelength = to_num<T>(parameters(3).data());
	Vector<T> step_sphere;
	if (parameters.GetM() > 4)
	  {
	    step_sphere.Reallocate(parameters.GetM()-4);
	    for (int k = 0; k < (parameters.GetM()-4); k++)
	      step_sphere(k) = to_num<T>(parameters(k+4));
	  }
	
	// creation of the spherical domain a <= r <= b
	CreateSphericalCrown(radius_sphere_inside, radius_sphere_outside,
			     nb_points_by_wavelength,
			     true, true, step_sphere, false);
      }
    else if (!parameters(0).compare("SPHERE_CUBE"))
      {
	// regular mesh of domain between a cube and a sphere
	T radius_sphere_inside = to_num<T>(parameters(1));
	T radius_sphere_outside = to_num<T>(parameters(2));
	T nb_points_by_wavelength = to_num<T>(parameters(3));
	
	Vector<T> step_sphere;
	
	bool all_ref = false;
	if (parameters.GetM() > 4)
	  all_ref = true;
	
	// creation of the spherical domain a <= r <= b
	CreateSphericalCrown(radius_sphere_inside, radius_sphere_outside,
			     nb_points_by_wavelength,
			     true, false, step_sphere, all_ref);
      }
    else if ((!parameters(0).compare("BALL"))||(!parameters(0).compare("BALL_CUBE")))
      {
	// regular mesh of a ball containing another ball (dielectric sphere)
	// the extern boundary can be a cube instead of a sphere
	bool is_ext_sphere = true;
	if (!parameters(0).compare("BALL_CUBE"))
	  is_ext_sphere = false;
        
	T radius_sphere_inside = to_num<T>(parameters(1));
	T radius_sphere_outside = to_num<T>(parameters(2));
	T nb_points_by_wavelength = to_num<T>(parameters(3));
	int reference_outside_boundary = to_num<int>(parameters(4));
        bool curved_intern
          = (this->GetCurveType(reference_outside_boundary) == this->CURVE_SPHERE);
        bool curved_extern = (this->GetCurveType(3) == this->CURVE_SPHERE);
        if (curved_intern)
          this->SetCurveType(reference_outside_boundary, this->NO_CURVE);
        
        if (curved_extern)
          this->SetCurveType(3, this->NO_CURVE);
        
        if (!parameters(5).compare("REGULAR"))
	  {
	    if (type_mesh == HEXAHEDRAL_MESH)
	      {
		this->ReallocateVertices(7);
		this->Vertices(0).Init(T(0), T(0), T(0));
		this->Vertices(1).Init(T(1), T(0), T(0));
		this->Vertices(2).Init(T(-1), T(0), T(0));
		this->Vertices(3).Init(T(0), T(1), T(0));
		this->Vertices(4).Init(T(0), T(-1), T(0));
		this->Vertices(5).Init(T(0), T(0), T(1));
		this->Vertices(6).Init(T(0), T(0), T(-1));
		this->ReallocateElements(8); int ref_inside = 2;
		this->elements(0).InitTetrahedral(0, 1, 3, 5, ref_inside);
		this->elements(1).InitTetrahedral(0, 4, 1, 5, ref_inside);
		this->elements(2).InitTetrahedral(0, 4, 6, 1, ref_inside);
		this->elements(3).InitTetrahedral(0, 6, 3, 1, ref_inside);
		this->elements(4).InitTetrahedral(0, 3, 2, 5, ref_inside);
		this->elements(5).InitTetrahedral(0, 2, 4, 5, ref_inside);
		this->elements(6).InitTetrahedral(0, 6, 4, 2, ref_inside);
		this->elements(7).InitTetrahedral(0, 3, 6, 2, ref_inside);
		ReallocateBoundariesRef(8);
		faces_ref(0).InitTriangular(1, 3, 5, reference_outside_boundary);
		faces_ref(1).InitTriangular(1, 6, 3, reference_outside_boundary);
		faces_ref(2).InitTriangular(1, 4, 5, reference_outside_boundary);
		faces_ref(3).InitTriangular(1, 6, 4, reference_outside_boundary);
		faces_ref(4).InitTriangular(2, 3, 5, reference_outside_boundary);
		faces_ref(5).InitTriangular(2, 5, 4, reference_outside_boundary);
		faces_ref(6).InitTriangular(2, 6, 3, reference_outside_boundary);
		faces_ref(7).InitTriangular(2, 4, 6, reference_outside_boundary);
                
                this->ReorientElements();
                this->FindConnectivity();
                this->ProjectPointsOnCurves();
	      }
	    else
	      {
		int ref_inside = 2; TinyVector<int, 6> ref_boundary;
		TinyVector<int, 3> nb_pts_; nb_pts_.Fill(3);
		ref_boundary.Fill(reference_outside_boundary);
		this->CreateRegularMesh(R3(-1,-1,-1), R3(1,1,1), nb_pts_,
					ref_inside, ref_boundary, type_mesh);
	      }
            
            R_N point; Real_wp norme_init, coef_max;
            for (int i = 0; i < this->GetNbVertices(); i++)
              {
                norme_init = Norm2(this->Vertex(i));
                if (norme_init > 1e-5)
                  {
                    point = this->Vertex(i);
                    coef_max = max(abs(point(0)), abs(point(1)));
                    coef_max = max(abs(point(2)), coef_max);
                    point *= 1.0/coef_max;
                    this->Vertex(i) *= 1.0/Norm2(point);
                  }
              }
          }
	else
	  {
	    string file_mesh(this->mesh_path);
	    file_mesh += parameters(5);
	    // we read initial mesh
	    Read(file_mesh);
	  }
	
	CreateSphericalBall(radius_sphere_inside, radius_sphere_outside,
			    nb_points_by_wavelength, reference_outside_boundary,
			    is_ext_sphere, type_mesh);
        
        if (curved_intern)
          this->SetCurveType(reference_outside_boundary, this->CURVE_SPHERE);
	else
	  this->SetCurveType(reference_outside_boundary, this->NO_CURVE);

        if (curved_extern)
          this->SetCurveType(3, this->CURVE_SPHERE);
	else
	  this->SetCurveType(3, this->NO_CURVE);
	
	this->ProjectPointsOnCurves();
      }
    else if ((parameters(0) == "EXTRUDE") || (parameters(0) == "EXTRUDE_X")
             || (parameters(0) == "EXTRUDE_Y"))
      {
        IVect type_curve = this->GetCurveType();
        Vector<VectReal_wp> param_curve = this->GetCurveParameter();
        
        for (int ref = 1; ref <= this->GetNbReferences(); ref++)
          this->SetCurveType(ref, this->NO_CURVE);
        
	// extrusion of a 2-D mesh
	int nb = 1; int nb_layers = to_num<int>(parameters(nb++));
	// z-coordinate for each interface
	Vector<T> step_layer(nb_layers+1);
	for (int i = 0; i <= nb_layers; i++)
	  step_layer(i) = to_num<T>(parameters(nb++));
		
	// 2-D mesh file is read
	string file_mesh(this->mesh_path);
	file_mesh += parameters(nb++);
	Mesh<Dimension2, T> mesh_2d;
	mesh_2d.Read(file_mesh);
	if (type_mesh == HEXAHEDRAL_MESH)
	  if (mesh_2d.GetNbTriangles() > 0)
	    mesh_2d.SplitIntoQuadrilaterals();
	
	// number of subdivisions for each layer
	T hmean = 0.6*mesh_2d.GetMeshSize();
	IVect nbCells_layer(nb_layers);
	if (!parameters(nb).compare("AUTO"))
	  {
	    for (int i = 0; i < nb_layers; i++)
	      nbCells_layer(i) = toInteger(floor((step_layer(i+1)-step_layer(i))/hmean))+1;
	    
	    nb++;
	  }
	else
	  {
	    for (int i = 0; i < nb_layers; i++)
	      nbCells_layer(i) = to_num<int>(parameters(nb++));
	  }
	
        // surface references
        int nb_ref = to_num<int>(parameters(nb++));
	Vector<IVect> ref_surface(2*nb_layers);
	for (int i = 0; i < 2*nb_layers; i++)
	  {
	    ref_surface(i).Reallocate(nb_ref+1);
	    ref_surface(i).Fill(-1);
	    for (int j = 1; j <= nb_ref; j++)
	      ref_surface(i)(j) = to_num<int>(parameters(nb++));
	  }
	
	// volume references
	nb_ref = to_num<int>(parameters(nb++));
	Vector<IVect> ref_volume(nb_layers);
	for (int i = 0; i < nb_layers; i++)
	  {
	    ref_volume(i).Reallocate(nb_ref+1);
	    ref_volume(i).Fill(-1);
	    for (int j = 1; j <= nb_ref; j++)
	      ref_volume(i)(j) = to_num<int>(parameters(nb++));
	  }
		
	// ez -> axis of extrusion
	R3 vec_e1, vec_e2, vec_e3;
	vec_e1(0) = 1.0; vec_e2(1) = 1.0; vec_e3(2) = 1.0;
        if (parameters(0) == "EXTRUDE_X")
          {
            vec_e1(1) = 1.0;
            vec_e2(2) = 1.0;
            vec_e3(0) = 1.0;
          }

        if (parameters(0) == "EXTRUDE_Y")
          {
            vec_e1(0) = 1.0;
            vec_e2(2) = 1.0;
            vec_e3(1) = 1.0;
          }
        
	ExtrudeSurfaceMesh(mesh_2d, nbCells_layer, step_layer,
			   ref_surface, ref_volume, vec_e1, vec_e2, vec_e3);

        for (int ref = 1; ref < type_curve.GetM(); ref++)
          {
            this->SetCurveType(ref, type_curve(ref));
            this->SetCurveParameter(ref, param_curve(ref));
          }
        
        this->ProjectPointsOnCurves();

        if (regular_mesh_with_slight_modification == REGULAR_MESH_SLIGHTLY_MODIFIED)
          {
            FctScalingExtrudeMesh fct;
            fct.z0 = 4.0;
            fct.z1 = 6.0;
            fct.max_ratio = 1.4;
            ScaleVariableMesh(*this, fct);
          }
      }
    else
      {
	string file_mesh(this->mesh_path);
	file_mesh += parameters(0);
	Read(file_mesh);
      }
    
    // we change local numbering for tetrahedra (warburton algorithm)
    // ChangeLocalNumberingMesh();
      
    // refinement to apply to the initial mesh
    int factor_refinement_mesh = 1;
    
    for (int i = 1; i< parameters.GetM(); i++)
      {
	if (!parameters(i).compare("REFINED"))
	  to_num(parameters(i+1),factor_refinement_mesh);
      }
    
    if (factor_refinement_mesh < 1)
      factor_refinement_mesh = 1;

    // the user can ask a local refinement around a vertex
    int nb_vertices_to_be_refined = this->ListeVertices_to_be_refined.GetM();
    for (int i = 0; i < nb_vertices_to_be_refined; i++)
      {
	int num_vertex = this->FindVertexNumber(this->ListeVertices_to_be_refined(i));
	IVect num(1); num(0) = num_vertex; IVect num_edge;
	if (num_vertex != -1)
	  ApplyLocalRefinement(num, num_edge, this->LevelRefinement_Vertex(i),
			       this->RatioRefinement_Vertex(i));
	else
	  {
	    if (this->print_level >= 0)
	      cout<<"we didn't find the vertex " << this->ListeVertices_to_be_refined(i)
		  <<" in the mesh"<<endl;
	    
	    abort();
	  }
      }

    if (this->ParamRefinement_Area.GetM() > 0)
      {
	const Vector<string>& param = this->ParamRefinement_Area;
	int level_refinement = to_num<int>(param(0));
	T ratio_refinement = to_num<T>(param(1));
	int nb = 2; IVect num_vertex, num_edge;
	if (param(nb) == "Points")
	  {
	    nb++;
	    int n = to_num<int>(param(nb++));
	    num_vertex.Reallocate(n);
	    for (int i = 0; i < n; i++)
	      {
		R3 pt(to_num<T>(param(nb)), to_num<T>(param(nb+1)), to_num<T>(param(nb+2))); nb += 3;
		num_vertex(i) = this->FindVertexNumber(pt);
		if (num_vertex(i) == -1)
		  {
		    if (this->print_level >= 0)
		      cout<<"we didn't find the vertex " << pt <<" in the mesh"<<endl;
		  }		
	      }
	  }

	if (param(nb) == "Edges")
	  {
	    nb++;
	    int n = to_num<int>(param(nb++));
	    
	    for (int i = 0; i < n; i++)
	      {
		R3 pt0(to_num<T>(param(nb)), to_num<T>(param(nb+1)), to_num<T>(param(nb+2))); nb += 3;
		R3 pt1(to_num<T>(param(nb)), to_num<T>(param(nb+1)), to_num<T>(param(nb+2))); nb += 3;
		R3 vec_u = pt1 - pt0;
		for (int j = 0; j < this->GetNbEdges(); j++)
		  {
		    R3 ptA = this->Vertex(this->GetEdge(j).numVertex(0));
		    R3 ptB = this->Vertex(this->GetEdge(j).numVertex(1));
		    R3 ptC = T(0.5)*(ptA + ptB);
		    R3 vec_v = pt0 - ptC, vec_w = pt1 - ptC, vec_prod;
		    TimesProd(vec_v, vec_w, vec_prod);
		    if (NormInf(vec_prod) <= R3::threshold)
		      if (DotProd(vec_v, vec_w) < 0)
			{
			  vec_v = ptB - ptA;
			  TimesProd(vec_u, vec_v, vec_prod);
			  if (NormInf(vec_prod) <= R3::threshold)
			    num_edge.PushBack(j);
			}
		  }
	      }
	  }

	int new_ref = 0;
	if (param(nb) == "Ref")
	  {	    
	    nb++;
	    new_ref = to_num<int>(param(nb++));
	  }
	
	RemoveDuplicate(num_edge);

	ApplyLocalRefinement(num_vertex, num_edge, level_refinement, ratio_refinement, new_ref);
      }

    // SPLITTING OF MESH IF NEEDED
    // the mesh read is cut if finite element method accepts only triangles/tetrahedras 
    // or quadrilaterals/hexahedras
    if ((type_mesh == this->TETRAHEDRAL_MESH)&&(!IsOnlyTetrahedral()))
      SplitIntoTetrahedra();
    else if ((type_mesh == this->HEXAHEDRAL_MESH)&&(!IsOnlyHexahedral()))
      SplitIntoHexahedra();
    else
      {
	// DISP(irregular_mesh_double_cut); DISP(type_mesh);
	if ((this->irregular_mesh_double_cut)&&(type_mesh == this->TETRAHEDRAL_MESH))
	  {
	    SplitIntoHexahedra(); 
	    SplitIntoTetrahedra(); 
	  }
	else if ((this->irregular_mesh_double_cut)&&(type_mesh == this->HEXAHEDRAL_MESH))
	  {
	    SplitIntoTetrahedra();
	    SplitIntoHexahedra();
	  }
      }

    if (factor_refinement_mesh != 1)
      {
	// regular subdivision of the mesh
	// each edge is split into the same number of intervals
	Vector<T> step_subdiv(factor_refinement_mesh+1);
	for (int i = 0; i <= factor_refinement_mesh; i++)
	  step_subdiv(i) = T(i)/factor_refinement_mesh;
	
	SubdivideMesh(step_subdiv);
      }

    // PML layers are added if required
    for (int i = 0; i < this->pml_areas.GetM(); i++)
      this->pml_areas(i).AddPML(i, *this);
    
    RearrangeElements();
  }
  

  template<class T>
  void Mesh<Dimension3, T>::RearrangeElements()
  {
    UpdateReferenceElement();
  }

  
  //! creation of a new mesh with given vertices and elements
  /*!
    \param[out] sub_mesh created mesh
    \param[in] nb_vertices_subdomain number of vertices to be taken
    \param[in] nb_elt_subdomain number of elements to be taken
    \param[in] VertexOnSubdomain if true, the vertex i is on the sub-mesh
    \param[in] ElementOnSubdomain if true, the element i is on the sub-mesh
   */
  template<class T>
  void Mesh<Dimension3, T>::CreateSubmesh(Mesh<Dimension3>& sub_mesh, int nb_vertices_subdomain,
                                          int nb_elt_subdomain, const VectBool& VertexOnSubdomain,
                                          const VectBool& ElementOnSubdomain)
  {
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
	  sub_mesh.Vertices(nb_vertices_subdomain++) = this->Vertices(n);
	}
    
    IVect NewNumberElement(this->GetNbElt()); NewNumberElement.Fill(-1);
    // elements
    sub_mesh.ReallocateElements(nb_elt_subdomain);
    nb_elt_subdomain = 0;
    for (int e = 0; e < this->GetNbElt(); e++)
      if (ElementOnSubdomain(e))
	{
	  sub_mesh.Element(nb_elt_subdomain) = this->elements(e);
	  sub_mesh.Element(nb_elt_subdomain).ClearConnectivity();
          int nb_vert = this->Element(e).GetNbVertices();
          IVect num(nb_vert); 
	  for (int n = 0; n < nb_vert; n++)
            num(n) = NewNumberVertex(this->Element(e).numVertex(n));
          
          sub_mesh.Element(nb_elt_subdomain).Init(num, this->Element(e).GetReference());
	  
	  nb_elt_subdomain++;
	}
    
    // referenced edges
    int nb_edges_subdomain = 0;
    for (int e = 0; e < this->GetNbEdgesRef(); e++)
      {
        bool edge_on_subdomain = false;
	for (int n = 0; n < this->EdgeRef(e).GetNbElements(); n++)
	  if (ElementOnSubdomain(this->EdgeRef(e).numElement(n)))
	    edge_on_subdomain = true;
	
	if (edge_on_subdomain)
          nb_edges_subdomain++;
      }
    
    bool only_quad = sub_mesh.IsOnlyHexahedral();
    sub_mesh.SetGeometryOrder(this->GetGeometryOrder(), only_quad);
    sub_mesh.ReallocateEdgesRef(nb_edges_subdomain);
    nb_edges_subdomain = 0;
    for (int e = 0; e < this->GetNbEdgesRef(); e++)
      {
        bool edge_on_subdomain = false;
	for (int n = 0; n < this->EdgeRef(e).GetNbElements(); n++)
	  if (ElementOnSubdomain(this->EdgeRef(e).numElement(n)))
	    edge_on_subdomain = true;
	
        if (edge_on_subdomain)
          {
            int n1 = NewNumberVertex(this->EdgeRef(e).numVertex(0));
            int n2 = NewNumberVertex(this->EdgeRef(e).numVertex(1));
            
            sub_mesh.EdgeRef(nb_edges_subdomain).
              Init(n1, n2, this->EdgeRef(e).GetReference());
            
            for (int k = 0; k < this->GetGeometryOrder()-1; k++)
              sub_mesh.PointsEdgeRef(nb_edges_subdomain, k) = this->PointsEdgeRef(e, k);
            
            nb_edges_subdomain++;
          }
      }
    
    sub_mesh.Type_curve = this->Type_curve; sub_mesh.Cond_curve = this->Cond_curve;
    sub_mesh.Body_curve = this->Body_curve; sub_mesh.Param_curve = this->Param_curve;
    sub_mesh.Param_known = this->Param_known;
    sub_mesh.reference_original_neighbor = this->reference_original_neighbor;
    
    // referenced faces
    int nb_faces_subdomain = 0;
    VectBool FaceOnSubdomain(this->GetNbBoundaryRef()); FaceOnSubdomain.Fill(false);
    for (int e = 0; e < this->GetNbBoundaryRef(); e++)
      {
	bool face_on_subdomain = false;
	for (int n = 0; n < faces_ref(e).GetNbElements(); n++)
	  if (ElementOnSubdomain(faces_ref(e).numElement(n)))
	    face_on_subdomain = true;
	
	if (face_on_subdomain)
	  {
	    FaceOnSubdomain(e) = true;
	    nb_faces_subdomain++;
	  }
      }
    
    sub_mesh.ReallocateBoundariesRef(nb_faces_subdomain);
    nb_faces_subdomain = 0;
    IVect num;
    for (int e = 0; e < this->GetNbBoundaryRef(); e++)
      if (FaceOnSubdomain(e))
	{
	  num.Reallocate(faces_ref(e).GetNbVertices());
	  for (int n = 0; n < faces_ref(e).GetNbVertices(); n++)
	    num(n) = NewNumberVertex(faces_ref(e).numVertex(n));
	  
	  sub_mesh.faces_ref(nb_faces_subdomain).Init(num, faces_ref(e).GetReference());
          sub_mesh.PointsFaceRef(nb_faces_subdomain) = this->PointsFaceRef(e);
	  nb_faces_subdomain++;
	}
    
    sub_mesh.ReorientElements();
    sub_mesh.FindConnectivity();
    sub_mesh.ProjectPointsOnCurves();
  }
    
  
  //! returns minimal and maximal value of theta (in cylindrical coordinates)
  template<class T>
  void Mesh<Dimension3, T>::
  GetTetaMinMax(const R_N& center, T& teta_min, T& teta_max) const
  {
    teta_min = T(10); teta_max = T(-10);
    // computation of teta_min, teta_max
    T teta, r; R_N vec_u;
    for (int i = 0; i < this->GetNbVertices(); i++)
      {
	// polar coordinates (r, teta)
	vec_u = this->Vertices(i) - center;
	r = sqrt(vec_u(0)*vec_u(0) + vec_u(1)*vec_u(1));
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
    
  
  //! periodize a section to obtain the global mesh by cyclicity
  template<class T>
  void Mesh<Dimension3, T>::PeriodizeMeshTeta(const R_N& center)
  {
    T teta, teta_min, teta_max;
    GetTetaMinMax(center, teta_min, teta_max);
    
    // number of periodicity 
    teta = teta_max - teta_min;
    int nper = toInteger(round(2.0*pi_wp/teta));
    
    R_N vec_u;
    int nb_old_vert = this->GetNbVertices();
    int nb_old_elt = this->GetNbElt();
    int nb_old_faces = this->GetNbBoundaryRef();
    int nb_old_edges = this->GetNbEdgesRef();
    int nb_vertices = nb_old_vert, nb_elt = nb_old_elt, nb_edges_ref = nb_old_edges,
      nb_faces_ref = nb_old_faces;
    this->ResizeVertices(nb_vertices*nper);
    this->ResizeElements(nb_elt*nper);
    this->ResizeEdgesRef(nb_edges_ref*nper);
    this->ResizeBoundariesRef(nb_faces_ref*nper);
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
            
            this->Vertex(nb_vertices + j)(2) = this->Vertex(j)(2); 
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
	    int n0 = this->EdgeRef(j).numVertex(0) + nb_vertices;
	    int n1 = this->EdgeRef(j).numVertex(1) + nb_vertices;
	    
	    this->EdgeRef(nb_edges_ref + j).Init(n0, n1, this->EdgeRef(j).GetReference());
	    for (int k = 0; k < this->GetGeometryOrder()-1; k++)
	      {
		vec_u = this->PointsEdgeRef(j, k) - center;
		point(0) = center(0) + cos_nteta*vec_u(0) - sin_nteta*vec_u(1);
		point(1) = center(1) + sin_nteta*vec_u(0) + cos_nteta*vec_u(1);
                point(2) = this->PointsEdgeRef(j, k)(2);
		this->PointsEdgeRef(nb_edges_ref+j, k) = point;
	      }
	  }

	// new faces
	for (int j = 0; j < nb_old_faces; j++)
	  {
            int nb_vert = this->BoundaryRef(j).GetNbVertices();
            IVect numf(nb_vert);
            for (int k = 0; k < nb_vert; k++)
              numf(k) = this->BoundaryRef(j).numVertex(k) + nb_vertices;
	    	    
	    this->BoundaryRef(nb_faces_ref + j).Init(numf, this->BoundaryRef(j).GetReference());
	    for (int k = 0; k < this->PointsFaceRef(j).GetM(); k++)
	      {
		vec_u = this->PointsFaceRef(j)(k) - center;
		point(0) = center(0) + cos_nteta*vec_u(0) - sin_nteta*vec_u(1);
		point(1) = center(1) + sin_nteta*vec_u(0) + cos_nteta*vec_u(1);
                point(2) = this->PointsFaceRef(j)(k)(2);
		this->PointsFaceRef(nb_faces_ref+j)(k) = point;
	      }
	  }
	
	nb_vertices += nb_old_vert;
	nb_elt += nb_old_elt;
	nb_faces_ref += nb_old_faces;
	nb_edges_ref += nb_old_edges;
      }
    
    this->ReorientElements();
    this->RemoveDuplicateVertices();
    this->ProjectPointsOnCurves();
  }
  
  
  //! local refinement of the mesh near vertices and edges
  /*!
    \param[in] num_vertex vertices where a refinement is asked
    \param[in] num_edge edges where a refinement is asked
    \param[in] nb_levels depth level of refinement
    \param[in] coef refinement ratio
    works for some hexahedral meshes
  */
  template<class T>
  void Mesh<Dimension3, T>::ApplyLocalRefinement(const Vector<int>& num_vertex,
						 const Vector<int>& num_edge_raff,
						 int nb_levels, const T& coef,
						 int new_ref)
  {
    Vector<char> TypeElement(this->GetNbElt());
    Vector<int> TypeVertex(this->GetNbVertices());
    Vector<int> TypeFace(this->GetNbFaces());
    Vector<int> TypeEdge(this->GetNbEdges());
    TypeVertex.Zero();
    TypeElement.Fill(0);
    TypeEdge.Fill(0);
    TypeFace.Fill(0);

    for (int i = 0; i < num_vertex.GetM(); i++)
      TypeVertex(num_vertex(i)) = 2;
    
    // we mark the edges to be refined as 1
    Vector<int> num_vertex_raff = num_vertex;
    for (int i = 0; i < num_edge_raff.GetM(); i++)
      {
	TypeEdge(num_edge_raff(i)) = 1;
	int n0 = this->GetEdge(num_edge_raff(i)).numVertex(0);
	int n1 = this->GetEdge(num_edge_raff(i)).numVertex(1);
	if (TypeVertex(n0) == 0)
	  {
	    TypeVertex(n0) = 1;
	    num_vertex_raff.PushBack(n0);
	  }

	if (TypeVertex(n1) == 0)
	  {
	    TypeVertex(n1) = 1;
	    num_vertex_raff.PushBack(n1);
	  }
      }
    
    // we get edges around the vertices and edges to be refined
    // edges that will be split are marked 2 + index
    Vector<int> num_edge, num_edge_tmp; int off_tmp = 0;
    for (int i = 0; i < num_vertex_raff.GetM(); i++)
      {
	this->GetNeighboringEdgesAroundVertex(num_vertex_raff(i), num_edge_tmp);
	for (int j = 0; j < num_edge_tmp.GetM(); j++)
	  if (TypeEdge(num_edge_tmp(j)) == 0)
	    {
	      TypeEdge(num_edge_tmp(j)) = 2 + off_tmp;
	      num_edge.PushBack(num_edge_tmp(j));
	      off_tmp++;
	    }
      }
    
    // faces in contact with edges to be refined
    IVect num_face; int nb_face2 = 0;
    for (int i = 0; i < num_edge.GetM(); i++)
      {
	int ne = num_edge(i);
	for (int j = 0; j < this->GetEdge(ne).GetNbFaces(); j++)
	  {
	    int nf = this->GetEdge(ne).numFace(j);
	    if (TypeFace(nf) == 0)
	      {
		TypeFace(nf) = 1;
		int num1 = -1, num2 = -1;
		for (int k = 0; k < this->Boundary(nf).GetNbEdges(); k++)
		  {
		    int ne0 = this->Boundary(nf).numEdge(k);
		    if (TypeEdge(ne0) >= 2)
		      {
			if (num1 == -1)
			  num1 = k;
			else
			  num2 = k;
		      }
		  }
		
		if (abs(num1-num2) != 2)
		  {
		    TypeFace(nf) = 2 + nb_face2;
		    nb_face2++;
		  }
		
		num_face.PushBack(nf);
	      }
	  }
      }

    // elements in contact with edges
    IVect num_elt; int nb_elt2 = 0;
    for (int i = 0; i < num_edge.GetM(); i++)
      {
	int ne = num_edge(i);
	for (int j = 0; j < this->GetEdge(ne).GetNbElements(); j++)
	  {
	    int nelt = this->GetEdge(ne).numElement(j);
	    if (TypeElement(nelt) == 0)
	      {
		TypeElement(nelt) = 1;
		int n2 = 0;
		for (int k = 0; k < this->Element(nelt).GetNbBoundary(); k++)
		  if (TypeFace(this->Element(nelt).numBoundary(k)) >= 2)
		    n2++;

		if (n2 == 3)
		  {
		    TypeElement(nelt) = 2;
		    nb_elt2++;
		  }
		
		num_elt.PushBack(nelt);
	      }
	  }
      }
    
    // on compte les nouveaux sommets
    int off_vert = this->GetNbVertices(), off_edge = num_edge.GetM(), off_face = num_edge.GetM() + nb_face2;
    off_edge += off_vert; off_face += off_vert;
    int nb_new_vertices = num_edge.GetM() + nb_face2 + nb_elt2;
    
    this->ResizeVertices(off_vert + nb_new_vertices);
    
    // creating points on the edges
    Real_wp L = 1.0/coef;
    for (int i = 0; i < num_edge.GetM(); i++)
      {
	int n1 = this->GetEdge(num_edge(i)).numVertex(0), n2 = this->GetEdge(num_edge(i)).numVertex(1);
	if (TypeVertex(n2) > 0)
	  {
	    n1 = n2;
	    n2 = this->GetEdge(num_edge(i)).numVertex(0);
	  }
	
	R3 new_pt = (1.0-L)*this->Vertex(n1) + L*this->Vertex(n2);
	this->Vertex(off_vert + i) = new_pt;
      }

    int nb_new_faces = 0;
    for (int i = 0; i < num_face.GetM(); i++)
      if (num_face(i) < this->GetNbBoundaryRef())
	{
	  if (TypeFace(num_face(i)) == 1)
	    nb_new_faces++;
	  else
	    nb_new_faces += 2;
	}
    
    int nb_old_face = this->GetNbBoundaryRef();
    this->ResizeBoundariesRef(nb_old_face + nb_new_faces);
    
    // creating points on the faces
    TinyVector<int, 4> nv; nb_new_faces = nb_old_face;
    off_tmp = off_edge;
    for (int i = 0; i < num_face.GetM(); i++)
      {
	int nf = num_face(i);
	int nb_edges_refined = 0;
	nv.Fill(-1);
	for (int j = 0; j < this->Boundary(nf).GetNbEdges(); j++)
	  {
	    int ne = this->Boundary(nf).numEdge(j);
	    int ne_loc = TypeEdge(ne);
	    if (ne_loc >= 2)
	      {
		nv(j) = off_vert + ne_loc-2;
		nb_edges_refined++;
	      }
	  }
	
	if ((nb_edges_refined != 2) || (this->Boundary(nf).GetNbEdges() != 4))
	  {
	    cout << "Case not handled" << endl;
	    abort();
	  }

	R3 pt0, pt1, pt2, pt3;
	int n0 = this->Boundary(nf).numVertex(0);
	int n1 = this->Boundary(nf).numVertex(1);
	int n2 = this->Boundary(nf).numVertex(2);
	int n3 = this->Boundary(nf).numVertex(3);
	pt0 = this->Vertex(n0);
	pt1 = this->Vertex(n1);
	pt2 = this->Vertex(n2);
	pt3 = this->Vertex(n3);

	int ref = this->Boundary(nf).GetReference();
	if (TypeFace(nf) >= 2)
	  {
	    // face with two neighboring refined edges
	    // => three sub-faces and a new point
	    Real_wp x = L, y = L;	
	    if ((nv(0) >= 0) && (nv(3) >= 0))
	      {
		// vertex 0
		if (nf < nb_old_face)
		  {
		    this->BoundaryRef(nf).InitQuadrangular(n0, nv(0), off_tmp, nv(3), ref);
		    this->BoundaryRef(nb_new_faces).InitQuadrangular(nv(0), n1, n2, off_tmp, ref);
		    this->BoundaryRef(nb_new_faces+1).InitQuadrangular(nv(3), off_tmp, n2, n3, ref);
		    nb_new_faces += 2;
		  }
	      }
	    else if ((nv(0) >= 0) && (nv(1) >= 0))
	      {
		// vertex 1
		x = 1.0-L;
		if (nf < nb_old_face)
		  {
		    this->BoundaryRef(nf).InitQuadrangular(nv(0), n1, nv(1), off_tmp, ref);
		    this->BoundaryRef(nb_new_faces).InitQuadrangular(n0, nv(0), off_tmp, n3, ref);
		    this->BoundaryRef(nb_new_faces+1).InitQuadrangular(off_tmp, nv(1),  n2, n3, ref);
		    nb_new_faces += 2;
		  }
	      }
	    else if ((nv(1) >= 0) && (nv(2) >= 0))
	      {
		// vertex 2
		x = 1.0-L; y = 1.0-L;
		if (nf < nb_old_face)
		  {
		    this->BoundaryRef(nf).InitQuadrangular(off_tmp, nv(1), n2, nv(2), ref);
		    this->BoundaryRef(nb_new_faces).InitQuadrangular(n0, n1, nv(1), off_tmp, ref);
		    this->BoundaryRef(nb_new_faces+1).InitQuadrangular(n0, off_tmp, nv(2), n3, ref);
		    nb_new_faces += 2;
		  }
	      }
	    else if ((nv(2) >= 0) && (nv(3) >= 0))
	      {
		// vertex 3
		y = 1.0-L;
		if (nf < nb_old_face)
		  {
		    this->BoundaryRef(nf).InitQuadrangular(nv(3), off_tmp, nv(2), n3, ref);
		    this->BoundaryRef(nb_new_faces).InitQuadrangular(n0, n1, off_tmp, nv(3), ref);
		    this->BoundaryRef(nb_new_faces+1).InitQuadrangular(off_tmp, n1, n2, nv(2), ref);
		    nb_new_faces += 2;
		  }
	      }
	    else
	      {
		cout << "Case not handled" << endl;
		abort();
	      }
	    
	    this->Vertex(off_tmp) = (1.0-x)*(1.0-y)*pt0 + x*(1.0-y)*pt1 + x*y*pt2 + (1.0-x)*y*pt3;
	    off_tmp++;
	  }
	else if (nf <= this->GetNbBoundaryRef())
	  {
	    // face with two opposite refined edges
	    // => two sub-faces
	    if ((nv(0) >= 0) && (nv(2) >= 0))
	      {
		this->BoundaryRef(nf).InitQuadrangular(n0, nv(0), nv(2), n3, ref);
		this->BoundaryRef(nb_new_faces).InitQuadrangular(nv(0), n1, n2, nv(2), ref);
		nb_new_faces++;
	      }
	    else if ((nv(1) >= 0) && (nv(3) >= 0))
	      {
		this->BoundaryRef(nf).InitQuadrangular(n0, n1, nv(1), nv(3), ref);
		this->BoundaryRef(nb_new_faces).InitQuadrangular(nv(3), nv(1), n2, n3, ref);
		nb_new_faces++;
	      }
	    else
	      {
		cout << "Case not handled" << endl;
		abort();
	      }
	  }
      }
    
    // creating points on the elements
    int nb_new_elt = num_elt.GetM()*2 + nb_elt2;
    int offset_elt = this->GetNbElt();
    this->ResizeElements(offset_elt + nb_new_elt);
    TinyVector<int, 12> nvb; TinyVector<int, 6> nfb;
    for (int i = 0; i < num_elt.GetM(); i++)
      {
	int n = num_elt(i);
	int nb_edges_refined = 0;
	nvb.Fill(-1);
	for (int j = 0; j < this->Element(n).GetNbEdges(); j++)
	  {
	    int ne = this->Element(n).numEdge(j);
	    int ne_loc = TypeEdge(ne);
	    
	    if (ne_loc >= 2)
	      {
		nvb(j) = off_vert + ne_loc-2;
		nb_edges_refined++;
	      }
	  }

	nfb.Fill(-1);
	int nb_faces_refined = 0, nb_faces_one = 0;
	for (int j = 0; j < this->Element(n).GetNbFaces(); j++)
	  {
	    int nf = this->Element(n).numFace(j);
	    int nf_loc = TypeFace(nf);
	    if (nf_loc >= 2)
	      {
		nfb(j) = off_edge + nf_loc-2;
		nb_faces_refined++;
	      }
	    else if (nf_loc == 1)
	      nb_faces_one++;
	  }
	
	R3 pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7;
	int n0 = this->Element(n).numVertex(0);
	int n1 = this->Element(n).numVertex(1);
	int n2 = this->Element(n).numVertex(2);
	int n3 = this->Element(n).numVertex(3);
	int n4 = this->Element(n).numVertex(4);
	int n5 = this->Element(n).numVertex(5);
	int n6 = this->Element(n).numVertex(6);
	int n7 = this->Element(n).numVertex(7);
	pt0 = this->Vertex(n0);
	pt1 = this->Vertex(n1);
	pt2 = this->Vertex(n2);
	pt3 = this->Vertex(n3);
	pt4 = this->Vertex(n4);
	pt5 = this->Vertex(n5);
	pt6 = this->Vertex(n6);
	pt7 = this->Vertex(n7);
        
	int ref = this->Element(n).GetReference();
	int ref_b = ref;
	if ((nb_levels <= 1) && (new_ref != 0))
	  ref_b = new_ref;

	if ((nb_edges_refined == 4) && (nb_faces_refined == 2) && (nb_faces_one == 2))
	  {
	    if ((nvb(3) >= 0) && (nvb(4) >= 0) && (nvb(1) >= 0) && (nvb(5) >= 0))
	      {
		this->Element(n).InitHexahedral(n0, n1, nvb(1), nvb(3), nvb(4), nvb(5), nfb(5), nfb(0), ref_b);
		this->Element(offset_elt).InitHexahedral(nvb(4), nvb(5), nfb(5), nfb(0), n4, n5, n6, n7, ref);
		this->Element(offset_elt+1).InitHexahedral(nvb(3), nvb(1), n2, n3, nfb(0), nfb(5), n6, n7, ref);
	      }
	    else if ((nvb(3) >= 0) && (nvb(7) >= 0) && (nvb(1) >= 0) && (nvb(6) >= 0))
	      {
		this->Element(n).InitHexahedral(nvb(3), nvb(1), n2, n3, nfb(0), nfb(5), nvb(6), nvb(7), ref_b);
		this->Element(offset_elt).InitHexahedral(n0, n1, nvb(1), nvb(3), n4, n5, nfb(5), nfb(0), ref);
		this->Element(offset_elt+1).InitHexahedral(nfb(0), nfb(5), nvb(6), nvb(7), n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(11) >= 0) && (nvb(7) >= 0) && (nvb(9) >= 0) && (nvb(6) >= 0))
	      {
		this->Element(n).InitHexahedral(nfb(0), nfb(5), nvb(6), nvb(7), nvb(11), nvb(9), n6, n7, ref_b);
		this->Element(offset_elt).InitHexahedral(n0, n1, nfb(5), nfb(0), n4, n5, nvb(9), nvb(11), ref);
		this->Element(offset_elt+1).InitHexahedral(n0, n1, n2, n3, nfb(0), nfb(5), nvb(6), nvb(7), ref);
	      }
	    else if ((nvb(11) >= 0) && (nvb(4) >= 0) && (nvb(9) >= 0) && (nvb(5) >= 0))
	      {
		this->Element(n).InitHexahedral(nvb(4), nvb(5), nfb(5), nfb(0), n4, n5, nvb(9), nvb(11), ref_b);
		this->Element(offset_elt).InitHexahedral(nfb(0), nfb(5), n2, n3, nvb(11), nvb(9), n6, n7, ref);
		this->Element(offset_elt+1).InitHexahedral(n0, n1, n2, n3, nvb(4), nvb(5), nfb(5), nfb(0), ref);
	      }
	    else if ((nvb(0) >= 0) && (nvb(4) >= 0) && (nvb(2) >= 0) && (nvb(7) >= 0))
	      {
		this->Element(n).InitHexahedral(n0, nvb(0), nvb(2), n3, nvb(4), nfb(1), nfb(4), nvb(7), ref_b);
	    	this->Element(offset_elt).InitHexahedral(nvb(0), n1, n2, nvb(2), nfb(1), n5, n6, nfb(4), ref);
	    	this->Element(offset_elt+1).InitHexahedral(nvb(4), nfb(1), nfb(4), nvb(7), n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(0) >= 0) && (nvb(5) >= 0) && (nvb(2) >= 0) && (nvb(6) >= 0))
	      {
		this->Element(n).InitHexahedral(nvb(0), n1, n2, nvb(2), nfb(1), nvb(5), nvb(6), nfb(4), ref_b);
	    	this->Element(offset_elt).InitHexahedral(n0, nvb(0), nvb(2), n3, n4, nfb(1), nfb(4), n7, ref);
	    	this->Element(offset_elt+1).InitHexahedral(nfb(1), nvb(5), nvb(6), nfb(4), n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(8) >= 0) && (nvb(5) >= 0) && (nvb(10) >= 0) && (nvb(6) >= 0))
	      {
		this->Element(n).InitHexahedral(nfb(1), nvb(5), nvb(6), nfb(4), nvb(8), n5, n6, nvb(10), ref_b);
	    	this->Element(offset_elt).InitHexahedral(n0, n1, n2, n3, nfb(1), nvb(5), nvb(6), nfb(4), ref);
	    	this->Element(offset_elt+1).InitHexahedral(n0, nfb(1), nfb(4), n3, n4, nvb(8), nvb(10), n7, ref);
	      }
	    else if ((nvb(8) >= 0) && (nvb(4) >= 0) && (nvb(10) >= 0) && (nvb(7) >= 0))
	      {
		this->Element(n).InitHexahedral(nvb(4), nfb(1), nfb(4), nvb(7), n4, nvb(8), nvb(10), n7, ref_b);
	    	this->Element(offset_elt).InitHexahedral(n0, n1, n2, n3, nvb(4), nfb(1), nfb(4), nvb(7), ref);
	    	this->Element(offset_elt+1).InitHexahedral(nfb(1), n1, n2, nfb(4), nvb(8), n5, n6, nvb(10), ref);
	      }
	    else if ((nvb(0) >= 0) && (nvb(3) >= 0) && (nvb(8) >= 0) && (nvb(11) >= 0))
	      {
		this->Element(n).InitHexahedral(n0, nvb(0), nfb(2), nvb(3), n4, nvb(8), nfb(3), nvb(11), ref_b);
	    	this->Element(offset_elt).InitHexahedral(nvb(0), n1, n2, nfb(2), nvb(8), n5, n6, nfb(3), ref);
	    	this->Element(offset_elt+1).InitHexahedral(nvb(3), nfb(2), n2, n3, nvb(11), nfb(3), n6, n7, ref);
	      }
	    else if ((nvb(0) >= 0) && (nvb(1) >= 0) && (nvb(8) >= 0) && (nvb(9) >= 0))
	      {
		this->Element(n).InitHexahedral(nvb(0), n1, nvb(1), nfb(2), nvb(8), n5, nvb(9), nfb(3), ref_b);
	    	this->Element(offset_elt).InitHexahedral(n0, nvb(0), nfb(2), n3, n4, nvb(8), nfb(3), n7, ref);
	    	this->Element(offset_elt+1).InitHexahedral(nfb(2), nvb(1), n2, n3, nfb(3), nvb(9), n6, n7, ref);
	      }
	    else if ((nvb(2) >= 0) && (nvb(1) >= 0) && (nvb(10) >= 0) && (nvb(9) >= 0))
	      {
		this->Element(n).InitHexahedral(nfb(2), nvb(1), n2, nvb(2), nfb(3), nvb(9), n6, nvb(10), ref_b);
	    	this->Element(offset_elt).InitHexahedral(n0, n1, nvb(1), nfb(2), n4, n5, nvb(9), nfb(3), ref);
	    	this->Element(offset_elt+1).InitHexahedral(n0, nfb(2), nvb(2), n3, n4, nfb(3), nvb(10), n7, ref);
	      }
	    else if ((nvb(2) >= 0) && (nvb(3) >= 0) && (nvb(10) >= 0) && (nvb(11) >= 0))
	      {
		this->Element(n).InitHexahedral(nvb(3), nfb(2), nvb(2), n3, nvb(11), nfb(3), nvb(10), n7, ref_b);
	    	this->Element(offset_elt).InitHexahedral(n0, n1, nfb(2), nvb(3), n4, n5, nfb(3), nvb(11), ref);
	    	this->Element(offset_elt+1).InitHexahedral(nfb(2), n1, n2, nvb(2), nfb(3), n5, n6, nvb(10), ref);
	      }
	    else
	      {
		cout << "Case not handled" << endl;
		DISP(nvb); DISP(nfb);
		abort();
	      }
	    
	    offset_elt += 2;
	  }
	else
	  {
	    if ((nb_edges_refined != 3) || (nb_faces_refined != 3) ||
		(this->Element(n).GetNbVertices() != 8) || (nb_faces_one != 0))
	      {
                DISP(nb_edges_refined); DISP(nb_faces_refined); DISP(nb_faces_one);
		cout << "Case not handled" << endl;
		abort();
	      }
	    
	    Real_wp x = L, y = L, z = L;
	    if ((nvb(0) >= 0) && (nvb(3) >= 0) && (nvb(4) >= 0))
	      {
		// vertex 0
		this->Element(n).InitHexahedral(n0, nvb(0), nfb(2), nvb(3),
						nvb(4), nfb(1), off_face + i, nfb(0), ref_b);
		
		this->Element(offset_elt).InitHexahedral(nvb(0), n1, n2, nfb(2),
							 nfb(1), n5, n6, off_face+i, ref);
		
		this->Element(offset_elt+1).InitHexahedral(nvb(3), nfb(2), n2, n3,
							   nfb(0), off_face+i, n6, n7, ref);
		
		this->Element(offset_elt+2).InitHexahedral(nvb(4), nfb(1), off_face+i, nfb(0),
							   n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(0) >= 0) && (nvb(1) >= 0) && (nvb(5) >= 0))
	      {
		// vertex 1
		x = 1.0 - L;
		this->Element(n).InitHexahedral(nvb(0), n1, nvb(1), nfb(2),
						nfb(1), nvb(5), nfb(5), off_face+i, ref_b);
		
		this->Element(offset_elt).InitHexahedral(n0, nvb(0), nfb(2), n3,
							 n4, nfb(1), off_face+i, n7, ref);
		
		this->Element(offset_elt+1).InitHexahedral(nfb(2), nvb(1), n2, n3,
							   off_face+i, nfb(5), n6, n7, ref);
	    
		this->Element(offset_elt+2).InitHexahedral(nfb(1), nvb(5), nfb(5), off_face+i,
							   n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(1) >= 0) && (nvb(2) >= 0) && (nvb(6) >= 0))
	      {
		// vertex 2
		x = 1.0 - L; y = 1.0 - L;
		this->Element(n).InitHexahedral(nfb(2), nvb(1), n2, nvb(2),
						off_face+i, nfb(5), nvb(6), nfb(4), ref_b);
		
		this->Element(offset_elt).InitHexahedral(n0, n1, nvb(1), nfb(2),
							 n4, n5, nfb(5), off_face+i, ref);
	    
		this->Element(offset_elt+1).InitHexahedral(n0, nfb(2), nvb(2), n3,
							   n4, off_face+i, nfb(4), n7, ref);
	    
		this->Element(offset_elt+2).InitHexahedral(off_face+i, nfb(5), nvb(6), nfb(4),
							   n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(2) >= 0) && (nvb(3) >= 0) && (nvb(7) >= 0))
	      {
		// vertex 3
		y = 1.0 - L;
		this->Element(n).InitHexahedral(nvb(3), nfb(2), nvb(2), n3,
						nfb(0), off_face+i, nfb(4), nvb(7), ref_b);
		
		this->Element(offset_elt).InitHexahedral(n0, n1, nfb(2), nvb(3),
							 n4, n5, off_face+i, nfb(0), ref);
		
		this->Element(offset_elt+1).InitHexahedral(nfb(2), n1, n2, nvb(2),
							   off_face+i, n5, n6, nfb(4), ref);
		
		this->Element(offset_elt+2).InitHexahedral(nfb(0), off_face+i, nfb(4), nvb(7),
							   n4, n5, n6, n7, ref);
	      }
	    else if ((nvb(4) >= 0) && (nvb(8) >= 0) && (nvb(11) >= 0))
	      {
		// vertex 4
		z = 1.0 - L;
		this->Element(n).InitHexahedral(nvb(4), nfb(1), off_face+i, nfb(0),
						n4, nvb(8), nfb(3), nvb(11), ref_b);
		
		this->Element(offset_elt).InitHexahedral(n0, n1, n2, n3,
							     nvb(4), nfb(1), off_face+i, nfb(0), ref);
		this->Element(offset_elt+1).InitHexahedral(nfb(1), n1, n2, off_face+i,
							   nvb(8), n5, n6, nfb(3), ref);
		
		this->Element(offset_elt+2).InitHexahedral(nfb(0), off_face+i, n2, n3,
							   nvb(11), nfb(3), n6, n7, ref);
	      }
	    else if ((nvb(5) >= 0) && (nvb(8) >= 0) && (nvb(9) >= 0))
	      {
		// vertex 5
		x = 1.0 - L; z = 1.0 - L;
		this->Element(n).InitHexahedral(nfb(1), nvb(5), nfb(5), off_face+i,
						nvb(8), n5, nvb(9), nfb(3), ref_b);
	    
		this->Element(offset_elt).InitHexahedral(n0, n1, n2, n3,
							 nfb(1), nvb(5), nfb(5), off_face+i, ref);
		this->Element(offset_elt+1).InitHexahedral(n0, nfb(1), off_face+i, n3,
							   n4, nvb(8), nfb(3), n7, ref);
		
		this->Element(offset_elt+2).InitHexahedral(off_face+i, nfb(5), n2, n3,
							   nfb(3), nvb(9), n6, n7, ref);
	      }
	    else if ((nvb(6) >= 0) && (nvb(9) >= 0) && (nvb(10) >= 0))
	      {
		// vertex 6
		x = 1.0 - L; y = 1.0 - L; z = 1.0 - L;
		this->Element(n).InitHexahedral(off_face+i, nfb(5), nvb(6), nfb(4),
						nfb(3), nvb(9), n6, nvb(10), ref_b);
		
		this->Element(offset_elt).InitHexahedral(n0, n1, n2, n3,
							 off_face+i, nfb(5), nvb(6), nfb(4), ref);
	    
		this->Element(offset_elt+1).InitHexahedral(n0, n1, nfb(5), off_face+i,
							   n4, n5, nvb(9), nfb(3), ref);
	    
		this->Element(offset_elt+2).InitHexahedral(n0, off_face+i, nfb(4), n3,
							   n4, nfb(3), nvb(10), n7, ref);
	      }
	    else if ((nvb(7) >= 0) && (nvb(10) >= 0) && (nvb(11) >= 0))
	      {
		// vertex 7
		y = 1.0 - L; z = 1.0 - L;
		this->Element(n).InitHexahedral(nfb(0), off_face+i, nfb(4), nvb(7),
						nvb(11), nfb(3), nvb(10), n7, ref_b);
		
		this->Element(offset_elt).InitHexahedral(n0, n1, n2, n3,
							 nfb(0), off_face+i, nfb(4), nvb(7), ref);
		
		this->Element(offset_elt+1).InitHexahedral(n0, n1, off_face+i, nfb(0),
							   n4, n5, nfb(3), nvb(11), ref);
		
		this->Element(offset_elt+2).InitHexahedral(off_face+i, n1, n2, nfb(4),
							   nfb(3), n5, n6, nvb(10), ref);
	      }
	    else
	      {
		cout << "Case not handled" << endl;
		abort();
	      }

	    offset_elt += 3;
	    this->Vertex(off_tmp) = (1.0-x)*(1.0-y)*(1.0-z)*pt0 + x*(1.0-y)*(1.0-z)*pt1
	      + x*y*(1.0-z)*pt2 + (1.0-x)*y*(1.0-z)*pt3 + (1.0-x)*(1.0-y)*z*pt4
	      + x*(1.0-y)*z*pt5 + x*y*z*pt6 + (1.0-x)*y*z*pt7;

	    off_tmp++;
	  }
      }

    // on stocke les aretes
    Vector<TinyVector<int, 2> > connec_edge(num_edge_raff.GetM());
    for (int i = 0; i < num_edge_raff.GetM(); i++)
      {
	int ne = num_edge_raff(i);
	connec_edge(i)(0) = this->GetEdge(ne).numVertex(0);
	connec_edge(i)(1) = this->GetEdge(ne).numVertex(1);
      }
    
    this->ReorientElements();
    this->FindConnectivity();

    if (nb_levels > 1)
      {
	Vector<int> num_edge_new(num_edge_raff.GetM());
	for (int i = 0; i < num_edge_raff.GetM(); i++)
	  num_edge_new(i) = this->FindEdgeWithExtremities(connec_edge(i)(0), connec_edge(i)(1));
	
	ApplyLocalRefinement(num_vertex, num_edge_new, nb_levels-1, coef, new_ref);
	return;
      }
    
    this->ProjectPointsOnCurves();

    //this->Write("refine.mesh");
  }

  
  //! displays details about class Mesh<Dimension3>
  template<class T>
  ostream& operator <<(ostream& out, const Mesh<Dimension3, T>& mesh)
  {
    out<<"3-D Mesh with "<<mesh.GetNbElt()<<" elements : "<<mesh.GetNbHexahedra()<<
      " hexahedras,  "<<mesh.GetNbPyramids()<<" pyramids, " << mesh.GetNbWedges() <<
       " wedges and "<<mesh.GetNbTetrahedra()<<" tetrahedras"<<endl;
    
    out<<"There are "<<mesh.GetNbFaces()<<" faces : "<<mesh.GetNbQuadrangles()<<
      " quadrangles and "<<mesh.GetNbTriangles()<<" triangles"<<endl;
    
    out<<"There are "<<mesh.GetNbEdges()<<" edges"<<endl;
    
    out<<"There are "<<mesh.GetNbVertices()<<" vertices"<<endl;

    out<<"There are "<<mesh.GetNbBoundaryRef()<<" referenced faces : "<<
      mesh.GetNbQuadranglesRef()<<" referenced quadrangles and "
       <<mesh.GetNbTrianglesRef()<<" referenced triangles"<<endl;
    
    return out;
  }

} // end namespace

#define MONTJOIE_FILE_MESH3D_CXX
#endif




