#ifndef MONTJOIE_FILE_MESH3D_BOUNDARIES_CXX

namespace Montjoie
{

  /**********************
   * PmlRegionParameter *
   **********************/


  //! returns average length of the extern boundary
  template<class T>
  T PmlRegionParameter<Dimension3, T>
  ::GetAverageLengthOnBoundary(const Mesh<Dimension3>& mesh_leaf) const
  {
    Real_wp xmin = mesh_leaf.GetXmin(); Real_wp xmax = mesh_leaf.GetXmax();
    Real_wp ymin = mesh_leaf.GetYmin(); Real_wp ymax = mesh_leaf.GetYmax();
    Real_wp zmin = mesh_leaf.GetZmin(); Real_wp zmax = mesh_leaf.GetZmax();
    
    // computation of the number of layers if necessary
    VectBool VertexOnBoundary(mesh_leaf.GetNbVertices());
    VertexOnBoundary.Fill(false);
    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
      if ((abs(mesh_leaf.Vertex(i)(0) - xmax) <= TinyVector<T, 3>::threshold)
	  ||(abs(mesh_leaf.Vertex(i)(0) - xmin) <= TinyVector<T, 3>::threshold)
	  ||(abs(mesh_leaf.Vertex(i)(1) - ymax) <= TinyVector<T, 3>::threshold)
	  ||(abs(mesh_leaf.Vertex(i)(1) - ymin) <= TinyVector<T, 3>::threshold)
	  ||(abs(mesh_leaf.Vertex(i)(2) - zmax) <= TinyVector<T, 3>::threshold)
	  ||(abs(mesh_leaf.Vertex(i)(2) - zmin) <= TinyVector<T, 3>::threshold))
	VertexOnBoundary(i) = true;
    
    // on fait la moyenne de la longueur des aretes du bord
    int nb_edges_boundary = 0; T step_h = 0.0;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      for (int j = 0; j < mesh_leaf.BoundaryRef(i).GetNbVertices(); j++)
	{
	  int n0 = mesh_leaf.BoundaryRef(i).numVertex(j);
	  int n1 = mesh_leaf.BoundaryRef(i).numVertex((j+1)%mesh_leaf.BoundaryRef(i).GetNbVertices());
	  if (VertexOnBoundary(n0)&&VertexOnBoundary(n1))
	    {
	      nb_edges_boundary++;
	      step_h += mesh_leaf.Vertex(n0).Distance(mesh_leaf.Vertex(n1));
	    }
	}
    
    if (nb_edges_boundary != 0)
      step_h /= nb_edges_boundary;
    
    return step_h;
  }
  

  //! Add PML layers to the mesh
  /*!
    \param[in] nb_div number of layers to add
    \warning the extern boundary must be a rectangular box
   */
  template<class T>
  void PmlRegionParameter<Dimension3, T>::AddPMLElements(int nb_div, int num_pml,
							 Mesh<Dimension3>& mesh_leaf)
  {
    int nb_elt = mesh_leaf.GetNbElt();
    if ((abs(this->thicknessPML) < R3::threshold)||(nb_div <= 0))
      return;

    Real_wp xmin = mesh_leaf.GetXmin(), xmax = mesh_leaf.GetXmax();
    Real_wp ymin = mesh_leaf.GetYmin(), ymax = mesh_leaf.GetYmax();
    Real_wp zmin = mesh_leaf.GetZmin(), zmax = mesh_leaf.GetZmax();
    
    // on extrude suivant chacune des coordonnees
    T delta = this->thicknessPML;
    if ((this->pml_add_axis_x == this->PML_BOTH_SIDES)
        ||(this->pml_add_axis_x == this->PML_NEGATIVE_SIDE))
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
    
    if ((this->pml_add_axis_z == this->PML_BOTH_SIDES)
        ||(this->pml_add_axis_z == this->PML_NEGATIVE_SIDE))
      mesh_leaf.ExtrudeCoordinate(2, nb_div, mesh_leaf.GetZmin(), -delta);
    
    if ((this->pml_add_axis_z == this->PML_BOTH_SIDES)
        ||(this->pml_add_axis_z == this->PML_POSITIVE_SIDE))
      mesh_leaf.ExtrudeCoordinate(2, nb_div, mesh_leaf.GetZmax(), delta);
    
    for (int i = nb_elt; i < mesh_leaf.GetNbElt(); i++)
      {
        int type_pml = 0;
        int nb_vert = mesh_leaf.Element(i).GetNbVertices();
        R3 center;
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
        
        if (center(2) < zmin)
          type_pml += PML_Zm;
        
        if (center(2) > zmax)
          type_pml += PML_Zp;

        mesh_leaf.Element(i).SetPML(num_pml, type_pml);    
      }
  }
  
  
  /******************
   * MeshBoundaries *
   ******************/
  

  template<class T>
  bool MeshBoundaries<Dimension3, T>::
  IsBoundaryCurved(int ne) const
  {
    if (ne >= this->GetLeafClass().GetNbBoundaryRef())
      return false;

    int ref = this->GetLeafClass().faces_ref(ne).GetReference();
    if (ref < 0)
      return true;
      
    if (this->Type_curve(ref) >= 1)
      return true;
    
    return false;
  }

  
  //! returns size of memory used by the object in bytes
  template<class T>
  size_t MeshBoundaries<Dimension3, T>::GetMemorySize() const
  {
    size_t taille = MeshBoundaries_Base<Dimension3, T>::GetMemorySize();
    taille += PointsEdgeRef.GetMemorySize();
    taille += Seldon::GetMemorySize(PointsFaceRef) + Seldon::GetMemorySize(value_parameter_boundary);
    taille += integer_parameter_boundary.GetMemorySize();
    taille += value_parameter_edge.GetMemorySize();
    taille += integer_parameter_edge.GetMemorySize();
    taille += triangle_reference.GetMemorySize() + quadrangle_reference.GetMemorySize();
    return taille;
  }
  
  
  //! modifies parameters of the mesh with a line of the data file
  template<class T>
  void MeshBoundaries<Dimension3, T>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    MeshBoundaries_Base<Dimension3, T>::SetInputData(description_field, parameters);
  }
  
  
  //! removes referenced edges with a null reference
  template<class T>
  void MeshBoundaries<Dimension3, T>::RemoveReference(int ref)
  {
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();
    int nb_edges_ref = mesh_leaf.edges_ref.GetM();
    for (int i = 0; i < mesh_leaf.edges_ref.GetM(); i++)
      if (mesh_leaf.edges_ref(i).GetReference() == ref)
	nb_edges_ref--;

    bool mesh_modified = false;
    if (nb_edges_ref < mesh_leaf.edges_ref.GetM())
      {
	if (nb_edges_ref == 0)
	  mesh_leaf.ReallocateEdgesRef(0);
	else
	  {
	    Vector<Edge<Dimension3> > old_edges; old_edges = mesh_leaf.edges_ref;
	    Matrix<R3> OldPoints = mesh_leaf.PointsEdgeRef;
	    mesh_leaf.ResizeEdgesRef(nb_edges_ref);
	    
	    int nb = 0;
	    for (int i = 0; i < old_edges.GetM(); i++)
	      if (old_edges(i).GetReference() != ref)
		{
		  for (int k = 0; k < OldPoints.GetN(); k++)
		    PointsEdgeRef(nb, k) = OldPoints(i, k);
		
		  mesh_leaf.edges_ref(nb++) = old_edges(i);
		}
	  }

	mesh_modified = true;
      }
    
    int nb_faces_ref = mesh_leaf.faces_ref.GetM();
    for (int i = 0; i < mesh_leaf.faces_ref.GetM(); i++)
      if (mesh_leaf.faces_ref(i).GetReference() == ref)
	nb_faces_ref--;
    
    if (nb_faces_ref < mesh_leaf.faces_ref.GetM())
      {
	if (nb_faces_ref == 0)
	  mesh_leaf.ReallocateBoundariesRef(0);
	else
	  {
	    Vector<Face<Dimension3> > old_faces; old_faces = mesh_leaf.faces_ref;
	    Vector<VectR3> OldPoints = PointsFaceRef;
	    Vector<VectReal_wp> old_value = value_parameter_boundary;
	    Vector<TinyVector<int, 4> > old_integer = integer_parameter_boundary;

	    mesh_leaf.ResizeBoundariesRef(nb_faces_ref);
	    int nb = 0;
	    for (int i = 0; i < old_faces.GetM(); i++)
	      if (old_faces(i).GetReference() != ref)
		{
		  PointsFaceRef(nb) = OldPoints(i); 
		  
		  value_parameter_boundary(nb) = old_value(i);
		  integer_parameter_boundary(nb) = old_integer(i); 
		  mesh_leaf.faces_ref(nb++) = old_faces(i);
		}
	  }
	mesh_modified = true;
      }

    if (mesh_modified)
      mesh_leaf.FindConnectivity();
  }
  
  
  /**************************************
   * Treatment of referenced boundaries *
   **************************************/
  

  //! permutes numbering of elements so that periodic facs have the same orientation  
  template<class T>
  void MeshBoundaries<Dimension3, T>::SortPeriodicBoundaries()
  {
    Mesh<Dimension3, T>& mesh = this->GetLeafClass();
    
    if (this->periodicity_references.GetM() <= 0)
      return;

    
    bool cyclic_domain = false;
    for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
      if (this->type_coord_periodicity(num_per) == BoundaryConditionEnum::PERIODIC_THETA)
        cyclic_domain = true;
    
    int nb_vertices = mesh.GetNbVertices();
    IVect permutation(nb_vertices), inv_perm(nb_vertices);
    inv_perm.Fill(-1);
    
    if (cyclic_domain)
      {
        // we assume that axis is Oz, and center is (0, 0, 0)
        R3 axis, center;
        center.Fill(0); axis.Init(0, 0, 1);
        
        // TypeVertex(i) will be 1 for reference1, 2 for reference 2
        // and 3 for common vertices (on the axis)
	int nb_vertices_ref1(0), nb_vertices_ref2(0);
	int nb_vertices_ref3 = 0;        
	IVect TypeVertex(nb_vertices);
	TypeVertex.Fill(-1);
	
	for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
          if (this->type_coord_periodicity(num_per) == BoundaryConditionEnum::PERIODIC_THETA) 
            {
              int ref1 = this->periodicity_references(num_per)(0);
              int ref2 = this->periodicity_references(num_per)(1);
              
              // TypeVertex(i) = 1 if vertex i belongs to ref1
              // TypeVertex(i) = 2 if vertex i belongs to ref2
              // TypeVertex(i) = 3 if vertex i belongs to both boundaries
              for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
                {
                  int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
                  if (mesh.BoundaryRef(i).GetReference() == ref1)
                    {
                      for (int j = 0; j < nb_vert; j++)
                        {
                          int nv = mesh.BoundaryRef(i).numVertex(j);
                          if (TypeVertex(nv) == 2)
                            TypeVertex(nv) = 3;
                          else if (TypeVertex(nv) == -1)
                            TypeVertex(nv) = 1;
                        }
                    }
                  
                  if (mesh.BoundaryRef(i).GetReference() == ref2)
                    {
                      for (int j = 0; j < nb_vert; j++)
                        {
                          int nv = mesh.BoundaryRef(i).numVertex(j);
                          if (TypeVertex(nv) == 1)
                            TypeVertex(nv) = 3;
                          else if (TypeVertex(nv) == -1)
                            TypeVertex(nv) = 2;
                        }
                    }
                }
            }
        
        // number of vertices on each boundary 
        R3 vec_u; Real_wp radius, rmax1 = 0, rmax2 = 0;
        R3 pt1, pt2;
        for (int i = 0; i < nb_vertices; i++)
          {
            if (TypeVertex(i) == 1)
              {
                radius = abs(mesh.Vertex(i)(0)) + abs(mesh.Vertex(i)(1));
                if (radius > rmax1)
                  {
                    rmax1 = radius;
                    pt1 = mesh.Vertex(i);
                  }
                
                nb_vertices_ref1++;
              }
            else if (TypeVertex(i) == 2)
              {
                radius = abs(mesh.Vertex(i)(0)) + abs(mesh.Vertex(i)(1));
                if (radius > rmax2)
                  {
                    rmax2 = radius;
                    pt2 = mesh.Vertex(i);
                  }
                
                nb_vertices_ref2++;
              }
            else if (TypeVertex(i) == 3)
              {
                TimesProd(mesh.Vertex(i), axis, vec_u);
                if (Norm2(vec_u) > R3::threshold)
                  {
                    cout << " axis for cyclic domains is Oz " << endl;
                    abort();
                  }
                
                nb_vertices_ref3++;
              }
          }
        
        this->periodic_center_polar = center;
        
        // cancelling z-component of points so that they are on the same plane
        pt1(2) = 0; pt2(2) = 0;
        Mlt(1.0/Norm2(pt1), pt1);
        Mlt(1.0/Norm2(pt2), pt2);
        
        // on repere l angle oriente theta 
	// signe du cosinus
	Real_wp signecos = DotProd(pt1, pt2);
	
	// theta compris entre (-pi et pi)
       	Real_wp theta = asin(pt1(0)*pt2(1) - pt1(1)*pt2(0));
	// si cos < 0 , l angle est en fait = pi -theta
	if (signecos < 0 ) 
	  theta = pi_wp - theta ;
	
	this->periodic_alpha_polar = theta;
        
	// we are retrieving vertices of Ref1 and Ref2
	// (vertices of Ref2 are rotated)	
	permutation.Fill();
        
        // polar coordinates of points (we omit theta)
	Vector<R2> RadiusRef1(nb_vertices_ref1), RadiusRef2(nb_vertices_ref2);
	IVect permut1(nb_vertices_ref1), permut2(nb_vertices_ref2);
        IVect permut3(nb_vertices_ref3);
        VectReal_wp PosVertAxe(nb_vertices_ref3);
        nb_vertices_ref1 = 0;
        nb_vertices_ref2 = 0;
        nb_vertices_ref3 = 0;
	for (int i = 0; i < nb_vertices; i++)
	  {
	    if (TypeVertex(i) == 1)
	      {
		RadiusRef1(nb_vertices_ref1)(0)
                  = sqrt(square(mesh.Vertex(i)(0)) + square(mesh.Vertex(i)(1)));
                
                RadiusRef1(nb_vertices_ref1)(1) = mesh.Vertex(i)(2);
                                
		permut1(nb_vertices_ref1) = i;
		nb_vertices_ref1++;
	      }
	    else if (TypeVertex(i) == 2)
	      {
		RadiusRef2(nb_vertices_ref2)(0)
                  = sqrt(square(mesh.Vertex(i)(0)) + square(mesh.Vertex(i)(1)));
                
                RadiusRef2(nb_vertices_ref2)(1) = mesh.Vertex(i)(2);
                
		permut2(nb_vertices_ref2) = i;
		nb_vertices_ref2++;
	      }
	    else if (TypeVertex(i) == 3)
	      {
                PosVertAxe(nb_vertices_ref3) = mesh.Vertex(i)(2);
                permut3(nb_vertices_ref3) = i;
		nb_vertices_ref3++;
	      }
	  }
	
	// sorting points with respect to radius (and z for equal radius)
        //DISP(RadiusRef1); DISP(RadiusRef2);
	Sort(nb_vertices_ref1, RadiusRef1, permut1);
	Sort(nb_vertices_ref2, RadiusRef2, permut2);
	Sort(nb_vertices_ref3, PosVertAxe, permut3);
        
	// modification of permutation	
	int nb = 0;
        for (int i = 0; i < nb_vertices_ref3; i++)
	  permutation(nb++) = permut3(i);
	
	for (int i = 0; i < nb_vertices_ref1; i++)
	  permutation(nb++) = permut1(i);
	
	for (int i = 0; i < nb_vertices_ref2; i++)
	  permutation(nb++) = permut2(i);
	
	for (int i = 0; i < nb_vertices; i++)
	  {
            if (TypeVertex(i) < 0)
              permutation(nb++) = i;
          }
        
        VectR3 OldVert = mesh.Vertex();
	for (int i = 0; i < nb_vertices; i++)
          mesh.Vertex(i) = OldVert(permutation(i));
        
        for (int i = 0; i < nb_vertices; i++)
          inv_perm(permutation(i)) = i;
        
        // changing edge numbers
        for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
          {
            int n1 = mesh.EdgeRef(i).numVertex(0);
            int n2 = mesh.EdgeRef(i).numVertex(1);
            int ref = mesh.EdgeRef(i).GetReference();
            n1 = inv_perm(n1);
            n2 = inv_perm(n2);
            mesh.EdgeRef(i).Init(n1, n2, ref);
          }
        
        // changing face numbers
        IVect num, num2;
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
            int ref = mesh.BoundaryRef(i).GetReference();
            num.Reallocate(nb_vert);
            for (int j = 0; j < nb_vert; j++)
              num(j) = inv_perm(mesh.BoundaryRef(i).numVertex(j));
            
            mesh.BoundaryRef(i).Init(num, ref);
          }
        
        // changing element numbers
        for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            int nb_vert = mesh.Element(i).GetNbVertices();
            int ref = mesh.Element(i).GetReference();
            num.Reallocate(nb_vert);
            for (int j = 0; j < nb_vert; j++)
              num(j) = inv_perm(mesh.Element(i).numVertex(j));
            
            mesh.Element(i).Init(num, ref);
          }
        
        inv_perm.Fill(-1);
        int Nv = nb_vertices_ref1 + nb_vertices_ref2 + nb_vertices_ref3;
        for (int i = 0; i < Nv; i++)
          inv_perm(i) = 1;
        
      }
    
    // if cartesian periodicity, we sort all the vertices of periodic boundaries
    bool cartesian_periodic = false;
    for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
      if (this->type_coord_periodicity(num_per) != BoundaryConditionEnum::PERIODIC_THETA) 
        cartesian_periodic = true;
    
    if (cartesian_periodic)
      {
        Vector<bool> VertexToSort(nb_vertices);
        VertexToSort.Fill(false);
        permutation.Fill();
	
        IVect ref_cond(mesh.GetNbReferences()+1);
        ref_cond.Fill(0);
        for (int num_per = 0; num_per < this->periodicity_references.GetM(); num_per++)
          if (this->type_coord_periodicity(num_per) != BoundaryConditionEnum::PERIODIC_THETA) 
            {
              int ref1 = this->periodicity_references(num_per)(0);
              int ref2 = this->periodicity_references(num_per)(1);
              ref_cond(ref1) = 1;
              ref_cond(ref2) = 1;
            }
        
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            int ref = mesh.BoundaryRef(i).GetReference();
            if (ref >= 0)
              if (ref_cond(ref) == 1)
                for (int k = 0; k < mesh.BoundaryRef(i).GetNbVertices(); k++)
                  {
                    int n1 = mesh.BoundaryRef(i).numVertex(k);
                    if (inv_perm(n1) < 0)
                      VertexToSort(n1) = true;
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
	
        // changing edge numbers
        for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
          {
            int n1 = mesh.EdgeRef(i).numVertex(0);
            int n2 = mesh.EdgeRef(i).numVertex(1);
            int ref = mesh.EdgeRef(i).GetReference();
            n1 = inv_perm(n1);
            n2 = inv_perm(n2);
            mesh.EdgeRef(i).Init(n1, n2, ref);
          }
        
        // changing face numbers
        IVect num, num2;
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          {
            int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
            int ref = mesh.BoundaryRef(i).GetReference();
            num.Reallocate(nb_vert);
            for (int j = 0; j < nb_vert; j++)
              num(j) = inv_perm(mesh.BoundaryRef(i).numVertex(j));
        
            mesh.BoundaryRef(i).Init(num, ref);
          }
        
        // changing element numbers
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
  }
    
  
  //! permutes extremity numbers of referenced edges so that first extremity number is lowest
  template<class T>
  void MeshBoundaries<Dimension3, T>::SortEdgesRef()
  {
    Mesh<Dimension3, T>& mesh = this->GetLeafClass();
    
    // sorting vertex numbers for referenced edges
    for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
      {
	int n1 = mesh.EdgeRef(i).numVertex(0);
	int n2 = mesh.EdgeRef(i).numVertex(1);
	int ref = mesh.EdgeRef(i).GetReference();
	if (n1 > n2)
	  {
	    mesh.EdgeRef(i).Init(n2, n1, ref);
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
  
  
  //! returns nodal numbers of triangles and quadrangles after rotation
  template<class T>
  void MeshBoundaries<Dimension3, T>::
  GetGeometryFaceRotation(Matrix<int>& RotationTri, Matrix<int>& RotationQuad) const
  {
    MeshNumbering<Dimension3, T>::
      GetRotationQuadrilateralFace(this->quadrangle_reference.PointsNodalND(), RotationQuad);
    
    MeshNumbering<Dimension3, T> ::
      GetRotationTriangularFace(this->triangle_reference.PointsNodalND(), RotationTri);
  }
  
  
  //! permutes vertex numbers of referenced faces so that first vertex number is the lowest
  /*!
    you should call ProjectPointsOnCurves in the case where permutation has been made on some faces
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::SortBoundariesRef(bool keep_orientation)
  {
    Mesh<Dimension3, T>& mesh = this->GetLeafClass();
    
    // sorting vertex numbers for referenced faces
    IVect num;
    int r = this->GetGeometryOrder();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	// retrieving vertex numbers
	int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
	num.Reallocate(nb_vert);
	for (int j = 0; j < nb_vert; j++)
	  num(j) = mesh.BoundaryRef(i).numVertex(j);
	
	int ref = mesh.BoundaryRef(i).GetReference();
	// we find the lowest vertex
	int rot = 0; int nmin = num(0);        
	for (int j = 1; j < nb_vert; j++)
	  if (num(j) < nmin)
	    {
	      rot = j;
	      nmin = num(j);
	    }
        
        if (num((rot+1)%nb_vert) > num((rot+nb_vert-1)%nb_vert))
          rot += nb_vert;
	
	// then we change the face
	if (rot > 0)
	  {
	    IVect new_num(nb_vert);
	    if (rot >= nb_vert)
              for (int j = 0; j < nb_vert; j++)
                new_num(j) = num((rot+nb_vert-j)%nb_vert);
            else
              for (int j = 0; j < nb_vert; j++)
                new_num(j) = num((rot+j)%nb_vert);
	    
	    mesh.BoundaryRef(i).Init(new_num, ref);
	    if (PointsFaceRef(i).GetM() > 0)
	      {
		Vector<R_N> OldPoints = PointsFaceRef(i);
		
		for (int k = 0; k < PointsFaceRef(i).GetM(); k++)
		  PointsFaceRef(i)(nodal_number_map
                                   .GetRotationFaceDof(rot, r, k, mesh.BoundaryRef(i)))
                    = OldPoints(k); 
	      }
	  }
      }
  }
  
  
  //! constructs bijection between faces and faces_ref
  template<class T>
  void MeshBoundaries<Dimension3, T>::ConstructCrossReferenceBoundaryRef()
  {
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();

    // we put in the array num_element of each face :
    // num_element(0) : lower reference domain
    // num_element(1) : higher reference domain
    // this sort is done, in order to have always an inward normal at the interface,
    // if we compute normal with the first element
    // if reference are equal, we consider that num_element(0) is the element
    // with the higher number
    for (int i = 0; i < mesh_leaf.GetNbFaces(); i++)
      if ((mesh_leaf.faces(i).GetNbElements() > 1)
	  &&(mesh_leaf.faces(i).numElement(1) >= 0))
	{
	  if (mesh_leaf.elements(mesh_leaf.faces(i).numElement(0)).GetReference()
	       > mesh_leaf.elements(mesh_leaf.faces(i).numElement(1)).GetReference())
	    {
	      mesh_leaf.faces(i).SetElement(mesh_leaf.faces(i).numElement(1),
					      mesh_leaf.faces(i).numElement(0));
	    }
	  else if (mesh_leaf.elements(mesh_leaf.faces(i).numElement(0)).GetReference()
		   == mesh_leaf.elements(mesh_leaf.faces(i).numElement(1)).GetReference())
	    {
	      if (mesh_leaf.faces(i).numElement(0) < mesh_leaf.faces(i).numElement(1))
		mesh_leaf.faces(i).SetElement(mesh_leaf.faces(i).numElement(1),
						mesh_leaf.faces(i).numElement(0));
	    }	  
	}
    

    // we make now the correspondance between faces and faces_ref
    int end_list=-1;
    
    int iglob, jglob, kglob, lglob,
      iglob2, jglob2, kglob2, lglob2, n0, n1, n2, n3;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	// we get points of the face, and sort them
	iglob = mesh_leaf.faces_ref(i).numVertex(0);
	jglob = mesh_leaf.faces_ref(i).numVertex(1);
	kglob = mesh_leaf.faces_ref(i).numVertex(2);
	int ref = mesh_leaf.faces_ref(i).GetReference();
	iglob2 = iglob;
	jglob2 = jglob;
	kglob2 = kglob; 	
	lglob = -1; lglob2 = -1;
	bool quad = (mesh_leaf.faces_ref(i).GetNbVertices() == 4);
	if (quad)
	  {
	    lglob = mesh_leaf.faces_ref(i).numVertex(3);
	    lglob2 = lglob;
	    Sort(iglob2, jglob2, kglob2, lglob2);
	  }
	else
	  Sort(iglob2, jglob2, kglob2);
	
	bool face_found = false;
	for (int e = mesh_leaf.head_faces_minv(iglob2); e != end_list;
	     e = mesh_leaf.next_face(e))
	  if (mesh_leaf.faces(e).GetNbVertices() == mesh_leaf.faces_ref(i).GetNbVertices())
	    {
	      n0 = mesh_leaf.faces(e).numVertex(0);
	      n1 = mesh_leaf.faces(e).numVertex(1);
	      n2 = mesh_leaf.faces(e).numVertex(2); n3 = -1;
	      if (quad)
		{
		  n3 = mesh_leaf.faces(e).numVertex(3);
		  Sort(n0,n1,n2,n3);
		}
	      else
		Sort(n0, n1, n2);
	      
	      // cout<<e<<" face of vertices "<<n0<<"  "<<n1<<"  "<<n2<<endl;
	      // we have found the face referenced
	      // we update faces_ref(i) and faces(e)
	      if (( n1 == jglob2)&&( n2 == kglob2)&&( n3 == lglob2))
		{
		  face_found = true;
		  mesh_leaf.faces(e).SetReference(ref);
		  mesh_leaf.faces(e).SetGeometryReference(mesh_leaf.faces_ref(i)
                                                          .GetGeometryReference());
                  mesh_leaf.faces_ref(i).SetElement(mesh_leaf.faces(e).numElement(0),
                                                    mesh_leaf.faces(e).numElement(1));
		  if (i != e)
		    {
		      cout << "A problem occured during FindConnectivity" << endl;
		      abort();
		    }
		}
	    }
	
	if (!face_found)
	  {
	    mesh_leaf.Write("test.mesh");
	    cout << "The referenced face " << i << " has not been found as a face of the mesh" << endl;	    
            abort();
	  }
      }
  }
  
  
  //! constructs bijection between edges and edges_ref
  template<class T>
  void MeshBoundaries<Dimension3, T>::ConstructCrossReferenceEdgeRef()
  { 
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();

    // loop on all the edges of the mesh
    for (int i = 0; i < mesh_leaf.GetNbEdgesRef(); i++)
      {
	// extremities of the referenced edge i
	int iglob = mesh_leaf.EdgeRef(i).numVertex(0);
	int jglob = mesh_leaf.EdgeRef(i).numVertex(1);
	
	// loop on all global edges, which leave vertex iglob
	bool edge_found = false;
	for (int e = mesh_leaf.head_minv(iglob); e != -1; e = mesh_leaf.next_edge(e))
	  {
	    if (mesh_leaf.edges(e).numVertex(1) == jglob)
	      {
		edge_found = true;
		// we have found the edge referenced
		// we update edges_ref(i) and edges(e)
		mesh_leaf.edges(e).Init(iglob, jglob, mesh_leaf.edges_ref(i).GetReference());
		mesh_leaf.edges_ref(i) = mesh_leaf.edges(e);
		
		// edges_ref are equal to edges, so that the extremities are in the same way
		if (e != i)
		  {
		    cout << "A problem occured during FindConnectivity" << endl;
		    abort();
		  }
	      }
	  }
	
	if (!edge_found)
	  {
	    mesh_leaf.Write("test.mesh");
	    // we didn't find the referenced edge
	    // the mesh has a problem, we stop program
	    //cout<<"Referenced edges are wrong. Check your mesh "<<endl;
	    //throw WrongMesh("MeshBoundaries<Dimension3, T>::ConstructCrossReferenceEdgeRef()",
            //	    string("Referenced edges are wrong. Check your mesh"));
	  }
      }
  }
  
  
  //! adds eventual faces on the boundary, which are not already stored in the array faces_ref_
  /*! 
    \param[in,out] nb_faces_ref_ number of referenced faces, modified by the method
    \param[in,out] faces_ref_ list of referenced faces, modified by the method
  */
  template<class T>
  void MeshBoundaries<Dimension3, T>::AddBoundaryFaces()
  {
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();

    // faces on the boundary have only one element (tetrahedron or hexahedron)
    // faces inside have two elements !
    int nb_old_faces_ref = mesh_leaf.GetNbBoundaryRef();
    int nb_isolated_faces = 0, nb_unref_faces =0;
    for (int i = nb_old_faces_ref; i < mesh_leaf.GetNbFaces(); i++)
      {
	// face with a reference in array faces and not found in the array faces_ref
	// it is an unreferenced face
	if (mesh_leaf.faces(i).GetReference() > 0)
	  nb_unref_faces++;          
	else if (mesh_leaf.faces(i).GetNbElements() == 1)
	  {
	    // face on a boundary and not in array faces_ref, it is an isolated face
	    nb_isolated_faces++;              
	  }
      }
    
    int nb_faces_ref = nb_old_faces_ref + nb_isolated_faces + nb_unref_faces;
    // if we find isolated or unreferenced edges, we add them to faces_ref_
    if (nb_faces_ref > nb_old_faces_ref)
      {
	int nb_old = nb_old_faces_ref;
	mesh_leaf.ResizeBoundariesRef(nb_faces_ref);
        
        // we get an available reference number
        int new_reference = this->GetNewReference();
	
	nb_faces_ref = nb_old_faces_ref;
	// we update the array faces_ref_
	for (int i = nb_old; i < mesh_leaf.GetNbFaces(); i++)
	  {
	    if (mesh_leaf.faces(i).GetReference() > 0)
	      {
		mesh_leaf.BoundaryRef(nb_faces_ref) = mesh_leaf.faces(i);
		mesh_leaf.BoundaryRef(nb_faces_ref)
		  .SetReference(mesh_leaf.faces(i).GetReference());
		nb_faces_ref++;
	      }
	    else if (mesh_leaf.faces(i).GetNbElements() == 1)
	      {
		mesh_leaf.faces(i).SetReference(new_reference);		  
		mesh_leaf.BoundaryRef(nb_faces_ref) = mesh_leaf.faces(i);
		nb_faces_ref++; 
	      }
	  }
      }
    
    // an available reference number
    int new_reference = this->GetNewReference();
    int nb_affected_faces = 0;
    
    // last check : boundary face with -1 reference -> we change this reference
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	if (mesh_leaf.Boundary(i).GetNbElements() == 1)
	  if (mesh_leaf.Boundary(i).GetReference() == -1)
	    {
	      mesh_leaf.Boundary(i).SetReference(new_reference);
	      mesh_leaf.BoundaryRef(i).SetReference(new_reference);
	      nb_affected_faces++;
	    }
      }
    
    // for these faces, curvature has been given in the file
    if (nb_affected_faces > 0)
      this->SetCurveType(new_reference, CURVE_FILE);
    
    if ((nb_faces_ref > nb_old_faces_ref) || (nb_affected_faces > 0))
      mesh_leaf.FindConnectivity();
  }
  

  //! projects vertices on curves
  template<class T>
  void MeshBoundaries<Dimension3, T>::ProjectVerticesOnCurves()
  {
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();   
    
    // We search points on the curves and project them exactly on the curves
    IVect RefVertices(mesh_leaf.GetNbVertices());
    IVect FaceVertices(mesh_leaf.GetNbVertices());
    IVect RefOldVertices(mesh_leaf.GetNbVertices());
    VectR3 OldVertex(mesh_leaf.Vertex());
    RefVertices.Fill(-1); RefOldVertices.Fill(-1); FaceVertices.Fill(-1);
    VectR3 PtToProject; Vector<int> RefSphereToProject, RefToProject;
    Vector<int> FaceSphereToProject, FaceToProject, NumToProject;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
	int ref = mesh_leaf.BoundaryRef(i).GetReference();
	if ((ref >= 1)&&(this->GetCurveType(ref) > 0)&&(this->GetCurveType(ref) != CURVE_FILE))
	  for (int j = 0; j < mesh_leaf.BoundaryRef(i).GetNbVertices(); j++)
	    {
              int nv = mesh_leaf.BoundaryRef(i).numVertex(j);
              int ref_before = RefVertices(nv);
              int ref_prev = RefOldVertices(nv);
              if (ref_before == -1)
                {
                  RefVertices(nv) = ref;
                  FaceVertices(nv) = i;
                  ProjectToCurve(mesh_leaf.Vertex(nv), ref);
                }
              else if ((ref != ref_before) && (ref_prev == -1))
                {
                  if (this->AreCurveEqual(ref, ref_before))
                    continue;
                  
                  RefOldVertices(nv) = ref;
                  
                  PtToProject.PushBack(OldVertex(nv));
                  NumToProject.PushBack(nv);
                  if (this->GetCurveType(ref) == CURVE_SPHERE)
                    {
                      RefToProject.PushBack(ref_before);
                      FaceToProject.PushBack(FaceVertices(nv));
                      RefSphereToProject.PushBack(ref);
                      FaceSphereToProject.PushBack(i);
                    }
                  else if (this->GetCurveType(ref_before) == CURVE_SPHERE)
                    {
                      mesh_leaf.Vertex(nv) = OldVertex(nv);
                      ProjectToCurve(mesh_leaf.Vertex(nv), ref);
                      
                      RefToProject.PushBack(ref);
                      FaceToProject.PushBack(i);
                      RefSphereToProject.PushBack(ref_before);
                      FaceSphereToProject.PushBack(FaceVertices(nv));
                    }
                  else
                    {
                      //cout << "Case not implemented" << endl;
                      //abort();
                    }
                }
              else if ((ref != ref_before) && (ref != ref_prev))
                {
                  //cout << "Case not implemented" << endl;
                  //abort();
                }
            }
      }

    // we compute parameters associated with vertices
    this->ComputeParameterValueOnCurves();
    
    // specific case of vertices on intersection of two curved surfaces
    for (int p = 0; p < NumToProject.GetM(); p++)
      {
        int nv = NumToProject(p);
        int nf = FaceToProject(p);
        int ref1 = RefToProject(p);
        int ref2 = RefSphereToProject(p);
                
        R3 center(this->GetCurveParameter(ref2)(0), this->GetCurveParameter(ref2)(1),
                  this->GetCurveParameter(ref2)(2));
        
        Real_wp radius = this->GetCurveParameter(ref2)(3);
        
        // vertices that belong to sphere are marked 0
        RefVertices.Fill(-1);
        for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
          if (mesh_leaf.BoundaryRef(i).GetReference() == ref2)
            for (int j = 0; j < mesh_leaf.BoundaryRef(i).GetNbVertices(); j++)
              RefVertices(mesh_leaf.BoundaryRef(i).numVertex(j)) = 0;

        // we find a neighbor vertex not in the sphere
        int pos = -1;
        int nb_vert = mesh_leaf.BoundaryRef(nf).GetNbVertices();
        for (int j = 0; j < nb_vert; j++)
          if (mesh_leaf.BoundaryRef(nf).numVertex(j) == nv)
            pos = j;

        int pos2 = (pos+1)%nb_vert;
        if (RefVertices(mesh_leaf.BoundaryRef(nf).numVertex(pos2)) == 0)
          {
            pos2 = (pos+nb_vert-1)%nb_vert;
            if (RefVertices(mesh_leaf.BoundaryRef(nf).numVertex(pos2)) == 0)
              {
                cout << "Impossible case : boundary on the sphere" << endl;
                abort();
              }
          }

        int n1 = nv, n2 = mesh_leaf.BoundaryRef(nf).numVertex(pos2);
        Sort(n1, n2);
        int num_edge = -1;
        for (int j = 0; j < nb_vert; j++)
          {
            int ne = mesh_leaf.BoundaryRef(nf).numEdge(j);
            if ((mesh_leaf.GetEdge(ne).numVertex(0) == n1)
                && (mesh_leaf.GetEdge(ne).numVertex(1) == n2))
              num_edge = ne;
          }
          
        IntersectionCurvedEdgeFace_System<T> sys(mesh_leaf, num_edge,
                                                 ref1, center, radius);
        
        // Newton method to find the intersection between the edge and the sphere
        T x(0); R3 point_current = mesh_leaf.Vertex(nv);
        if (n1 != nv)
          x = T(1);
        
        Real_wp f = sys.EvaluateFunction(x);
        Real_wp err = abs(f), err_prec = 2*err;
        if (err < 1e-12)
          err = err_prec;
        
        int nb_iter = 0;
        while ((err < err_prec) && (nb_iter < 50))
          {
            Real_wp df = sys.EvaluateDerivative(x);
            point_current = sys.GetPointOnSurface(x);
            if (df != 0)
              x = x - f / df;
            
            f = sys.EvaluateFunction(x);
            err_prec = err;
            err = abs(f);
            nb_iter++;
          }
        
        mesh_leaf.Vertex(nv) = point_current;
      }
  }
  
  
  //! detection of boundaries, projection of points on curved surfaces ...
  /*!
    \param[in] FaceBasis finite element (used for curved boundaries)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::ProjectPointsOnCurves()
  {  
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();   
    
    // we find parameters of curves
    this->FindParametersCurve();
    
    // we project vertices on curves
    this->ProjectVerticesOnCurves();
    
    if (this->lob_curve.GetOrder() > 1)
      {
	// computing PointsEdgeRef for predefined curves
	VectReal_wp xi(this->lob_curve.GetOrder()-1);
	for (int i = 1; i < this->lob_curve.GetOrder(); i++)
	  xi(i-1) = this->lob_curve.Points(i);
	
	if (this->PointsEdgeRef.GetM() != mesh_leaf.GetNbEdgesRef())
	  this->PointsEdgeRef.Resize(mesh_leaf.GetNbEdgesRef(), this->lob_curve.GetOrder()-1);
        
	for (int i = 0; i < mesh_leaf.GetNbEdgesRef(); i++)
	  {
            bool intersection_edge = false;
            int ref1 = -1, ref2 = -1;
            int num_loc1 = -1, num_loc2 = -1;
            for (int j = 0; j < mesh_leaf.EdgeRef(i).GetNbFaces(); j++)
              {
                int num_face = mesh_leaf.EdgeRef(i).numFace(j);
                int ref_face = mesh_leaf.Boundary(num_face).GetReference();
                if (this->GetCurveType(ref_face) > 0)
                  {
                    if (ref1 == -1)
                      {
                        ref1 = ref_face;
                        num_loc1 = j;
                      }
                    else 
                      {
                        if (ref_face != ref1)
                          {
                            ref2 = ref_face;
                            num_loc2 = j;
                            intersection_edge = true;
                            break;
                          }
                      }
                  }
              }
            
            if (intersection_edge)
              {
                if (this->GetCurveType(ref1) == CURVE_SPHERE)
                  {
                    int reft = ref1, num_loct = num_loc1;
                    ref1 = ref2; num_loc1 = num_loc2;
                    ref2 = reft; num_loc2 = num_loct;
                  }                
              }
            
	    int ref = mesh_leaf.EdgeRef(i).GetReference();
            if (this->GetCurveType(ref) == CURVE_FILE)
	      {
		// nothing, we consider that in that case,
		// PointsEdgeRef is filled during the reading or the mesh
		// or filled by the user
	      }
            else if (intersection_edge)
              {
                // treating intersection between predefined surfaces
                if (this->GetCurveType(ref2) != CURVE_SPHERE)
                  {
                    continue;
                    //cout << "Case not implemented" << endl;
                    //abort();
                  }
                
                R3 center(this->GetCurveParameter(ref2)(0), this->GetCurveParameter(ref2)(1),
                          this->GetCurveParameter(ref2)(2));
                
                Real_wp radius = this->GetCurveParameter(ref2)(3);

                int num_face1 = mesh_leaf.EdgeRef(i).numFace(num_loc1);
                int num_edge_loc1 = mesh_leaf.Boundary(num_face1).GetPositionBoundary(i);                
                for (int k = 0; k < xi.GetM(); k++)
                  {
                    Real_wp pos_edge = xi(k);
                    if (!mesh_leaf.Boundary(num_face1).GetOrientationEdge(num_edge_loc1))
                      pos_edge = T(1) - pos_edge;
                    
                    // Newton algorithm to find intersection
                    IntersectionCurvedEdgeFace_System<T> sys(mesh_leaf, num_face1, num_edge_loc1, pos_edge,
                                                             ref1, center, radius);
                    
                    Real_wp x(0); R3 point_current = sys.GetPointOnSurface(x);
                    Real_wp f = sys.EvaluateFunction(x);
                    Real_wp err = abs(f), err_prec = 2*err;
                    if (err < 1e-12)
                      err_prec = err;
                    
                    int nb_iter = 0;
                    while ((err < err_prec) && (nb_iter < 50))
                      {
                        Real_wp df = sys.EvaluateDerivative(x);
                        point_current = sys.GetPointOnSurface(x);
                        if (df != 0)
                          x = x - f / df;
                        
                        f = sys.EvaluateFunction(x);
                        err_prec = err;
                        err = abs(f);
                        nb_iter++;
                      }
                    
                    if (err > 1e-6)
                      {
                        cout << "Intersection not found" << endl;
                        abort();
                      }
                    
                    //DISP(k); DISP(point_current);
                    this->PointsEdgeRef(i, k) = point_current;
                  }
              }
	    else if (this->GetCurveType(ref) > 0)
	      {
		// predefined curve, we project the internal points 
		Vector<R_N> points_interm;
		//R_N ptA = mesh_leaf.Vertex(mesh_leaf.EdgeRef(i).numVertex(0));
		//R_N ptB = mesh_leaf.Vertex(mesh_leaf.EdgeRef(i).numVertex(1));
                this->GetPointsOnCurvedEdge(i, ref, xi, points_interm);
		for (int k = 0; k < xi.GetM(); k++)
		  this->PointsEdgeRef(i, k) = points_interm(k);
	      }
	    else 
	      {
		// straight boundary
		R_N ptA = mesh_leaf.Vertex(mesh_leaf.EdgeRef(i).numVertex(0));
		R_N ptB = mesh_leaf.Vertex(mesh_leaf.EdgeRef(i).numVertex(1));
		for (int k = 0; k < xi.GetM(); k++)
		  {
		    this->PointsEdgeRef(i, k).Zero();
		    Add(1.0-xi(k), ptA, this->PointsEdgeRef(i, k));
		    Add(xi(k), ptB, this->PointsEdgeRef(i, k));
		  }
	      }
	  }
	
	// computing PointsFaceRef for predefined curves
	int r = this->lob_curve.GetOrder();
	int Ntri = this->triangle_reference.GetNbPointsNodalElt() - 3*r;
	int Nquad = this->quadrangle_reference.GetNbPointsNodalElt() - 4*r;
	VectR2 xi_tri, xi_quad(Nquad);        
	if (Ntri >= 0)
          {
            xi_tri.Reallocate(Ntri);
            for (int i = 0; i < Ntri; i++)
              xi_tri(i) = this->triangle_reference.PointsNodalND(i+3*r);
          }
	
	for (int i = 0; i < Nquad; i++)
	  xi_quad(i) = this->quadrangle_reference.PointsNodalND(i+4*r);
	
	if (this->PointsFaceRef.GetM() != mesh_leaf.GetNbBoundaryRef())
	  this->PointsFaceRef.Resize(mesh_leaf.GetNbBoundaryRef());
	
	for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
	  {
	    bool tri = (mesh_leaf.BoundaryRef(i).GetNbVertices() == 3);
	    if (tri)
	      {
		if (this->PointsFaceRef(i).GetM() != Ntri)
		  this->PointsFaceRef(i).Resize(Ntri);
	      }
	    else
	      {
		if (this->PointsFaceRef(i).GetM() != Nquad)
		  this->PointsFaceRef(i).Resize(Nquad);
	      }

	    int ref = mesh_leaf.BoundaryRef(i).GetReference();		    
	    if (this->GetCurveType(ref) == CURVE_FILE)
	      {
		// nothing, we consider that in that case,
		// PointsFaceRef is filled during the reading or the mesh
		// or filled by the user
	      }
	    else if (this->GetCurveType(ref) > 0)
	      {
		// predefined curve, we project the internal points 
                /*this->ComputeInternalPointsFaceWithCurvedEdges(i);
                for (int k = 0; k < xi_tri.GetM(); k++)
                  ProjectToCurve(this->PointsFaceRef(i)(k), mesh_leaf.BoundaryRef(i).GetReference());
                */
                
		Vector<R_N> points_interm;
		if (tri)
		  {
		    this->GetPointsOnCurvedFace(i, ref, xi_tri, points_interm);
		    for (int k = 0; k < xi_tri.GetM(); k++)
		      this->PointsFaceRef(i)(k) = points_interm(k);
		  }
		else
		  {
		    this->GetPointsOnCurvedFace(i, ref, xi_quad, points_interm);
		    for (int k = 0; k < xi_quad.GetM(); k++)
		      this->PointsFaceRef(i)(k) = points_interm(k);
		  }
                
                //DISP(tri);
                //cout << "waiting" << endl; int test_input; cin >> test_input;
	      }
	    else 
	      {
                this->ComputeInternalPointsFaceWithCurvedEdges(i);
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
          mesh_leaf.zmin_ = min(mesh_leaf.zmin_, vec_u(2));
          mesh_leaf.zmax_ = max(mesh_leaf.zmax_, vec_u(2));
        }

    for (int i = 0; i < this->PointsFaceRef.GetM(); i++)
      for (int j = 0; j < this->PointsFaceRef(i).GetM(); j++)
        {
          vec_u = this->PointsFaceRef(i)(j);
          mesh_leaf.xmin_ = min(mesh_leaf.xmin_, vec_u(0));
          mesh_leaf.xmax_ = max(mesh_leaf.xmax_, vec_u(0));
          mesh_leaf.ymin_ = min(mesh_leaf.ymin_, vec_u(1));
          mesh_leaf.ymax_ = max(mesh_leaf.ymax_, vec_u(1));
          mesh_leaf.zmin_ = min(mesh_leaf.zmin_, vec_u(2));
          mesh_leaf.zmax_ = max(mesh_leaf.zmax_, vec_u(2));
        }
  }
  
  
  //! affecting a different reference number per connected surface
  template<class T>
  void MeshBoundaries<Dimension3, T>::RedistributeReferences()
  {
    Mesh<Dimension3>& mesh_leaf = this->GetLeafClass();
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
		for (int j = 0; j < mesh_leaf.faces_ref(new_face(i)).GetNbEdges(); j++)
		  {
		    int ne = mesh_leaf.faces_ref(new_face(i)).numEdge(j);
		    for (int k = 0; k < mesh_leaf.edges(ne).GetNbFaces(); k++)
		      {
			nf = mesh_leaf.edges(ne).numFace(k);
			if (nf < mesh_leaf.GetNbBoundaryRef())
			  if (!FaceUsed(nf))
			    {
			      num_face.PushBack(nf);
			      FaceUsed(nf) = true;
			    }
		      }
		  }
	      }
	    
	    new_face = num_face;
	  }
	
	// affecting new reference to those faces
	for (int i = 0; i < nb_new_face; i++)
	  mesh_leaf.faces_ref(ListeFace(i)).SetReference(ref);
	
	// next reference
	ref++;
      }
  }
    
  
  /**********************************
   * Treatment of curved boundaries *
   **********************************/
  
  
  //! checks if the parameters given by the user are correct
  template<class T>
  void MeshBoundaries<Dimension3, T>::CheckCurveParameter(int ref)
  {    
    Vector<T>& param = this->Param_curve(ref);
    switch (this->GetCurveType(ref))
      {
      case CURVE_SPHERE :
        {
	  if (param.GetM() < 4)
	    {
	      cout << "The sphere is determined by four parameters (Ox, Oy, Oz) and the radius" << endl;
	      cout << "Only " << param.GetM() << " are given " << endl;
	      cout << "Current parameters of reference " << ref << " : " << param << endl;
	      abort();
	    }
	  
	  if (param(3) <= 0)
	    {
	      cout << "The radius of the sphere must be strictly positive but is equal to " << param(3) << endl;
	      abort();
	    }
	}
        break;
      case CURVE_CYLINDER :
        {
	  if (param.GetM() < 7)
	    {
	      cout << "A cylinder is determined by seven parameters (Ox, Oy, Oz), (u_x, u_y, u_z) and the radius " << endl;
	      cout << "Only " << param.GetM() << " are given " << endl;
	      cout << "Current parameters of reference " << ref << " : " << param << endl;
	      abort();
	    }

	  R3 vec_u(param(3), param(4), param(5));
	  Mlt(1.0/Norm2(vec_u), vec_u);
	  param(3) = vec_u(0); param(4) = vec_u(1); param(5) = vec_u(2);
	  
	  if (param(6) <= 0)
	    {
	      cout << "The radius of the cylinder must be strictly positive but is equal to " << param(6) << endl;
	      abort();
	    }
        }
        break;
      case CURVE_CONIC :
        {
	  if (param.GetM() < 7)
	    {
	      cout << "A cone is determined by seven parameters (Ox, Oy, Oz), (u_x, u_y, u_z) and the angle " << endl;
	      cout << "Only " << param.GetM() << " are given " << endl;
	      cout << "Current parameters of reference " << ref << " : " << param << endl;
	      abort();
	    }

	  R3 vec_u(param(3), param(4), param(5));
	  Mlt(1.0/Norm2(vec_u), vec_u);
	  param(3) = vec_u(0); param(4) = vec_u(1); param(5) = vec_u(2);

	  if (param(6) <= 0)
	    {
	      cout << "The angle of the cone must be strictly positive but is equal to " << param(6) << endl;
	      abort();
	    }

	  // conversion in radians
	  param(6) = param(6)*pi_wp/180.0;
        }
        break;
      }
  }


  //! fills parameters s for a point of a predefined surface
  /*!
    \param[in] ref reference of the surface
    \param[in] pt point for which we search parameters 
    \param[in] param caracteristics of the surface (e.g. center and radius for a sphere)
    \param[out] s parameters of the surface such that F(s) = pt
    For example, if the surface is a sphere, s will be (teta, phi)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::
  GetParameterValueOnPoint(int ref, const R_N& pt, const VectReal_wp& param, R_N& s) const
  {
    switch (this->GetCurveType(ref))
      {
      case CURVE_SPHERE :
        {
          // x = r sin(teta) cos(phi)
          // y = r sin(teta) sin(phi)
          // z = r cos(teta)
          Real_wp r, teta, phi;
          r = param(3);
          Real_wp x = pt(0) - param(0);
          Real_wp y = pt(1) - param(1);
          Real_wp z = pt(2) - param(2);
          if (z >= r)
            teta = 0.0;
          else if (z <= -r)
            teta = pi_wp;
          else
            teta = acos(z/r);
          
          if (abs(x) <= epsilon_machine)
            {
              if (y < 0)
                phi = -0.5*pi_wp;
              else
                phi = 0.5*pi_wp;
            }
          else
            {
              phi = atan(y/x);
              if (x < 0)
                phi += pi_wp;
              else
                {
                  if (y < 0)
                    phi += 2.0*pi_wp;
                }
            }
          
          s(0) = teta;
          s(1) = phi;
        }
        break;
      case CURVE_CYLINDER :
        {
          R3 axis, u1, u2, PM;          
          axis(0) = param(3); axis(1) = param(4); axis(2) = param(5);
          GetVectorPlane(axis, u1, u2);
          Real_wp alpha = (pt(0) - param(0))*axis(0)
            + (pt(1) - param(1))*axis(1) + (pt(2) - param(2))*axis(2);
          
          PM = pt; Add(-alpha, axis, PM);
          Real_wp x = DotProd(u1, PM);
          Real_wp y = DotProd(u2, PM);
          Real_wp r = param(6);
          Real_wp teta;
          if (x >= r)
            teta = 0;
          else if ( x <= -r)
            teta = pi_wp;
          else
            teta = acos(x/r);
          
          if (y < 0)
            teta = -teta;
          
          s(0) = alpha;
          s(1) = teta;
        }
        break;
      case CURVE_CONIC :
        {
          R3 axis, u1, u2, PM;          
          axis(0) = param(3); axis(1) = param(4); axis(2) = param(5);
          GetVectorPlane(axis, u1, u2);
          Real_wp alpha = (pt(0) - param(0))*axis(0)
            + (pt(1) - param(1))*axis(1) + (pt(2) - param(2))*axis(2);
          
          PM = pt; Add(-alpha, axis, PM);
          Real_wp x = DotProd(u1, PM);
          Real_wp r = Norm2(PM);
          Real_wp teta;
          if (x >= r)
            teta = 0;
          else if ( x <= r)
            teta = pi_wp;
          else
            teta = acos(x/r);
          
          s(0) = alpha;
          s(1) = teta;
        }
        break;
      }
  }
  

  //! computation of parameters for curved surfaces
  template<class T>
  void MeshBoundaries<Dimension3, T>::ComputeParameterValueOnCurves()
  {
    Mesh<Dimension3>& mesh = this->GetLeafClass();

    // unmarking elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      mesh.Element(i).UnsetCurved();      
    
    // creation of referenced edges for predefined curves
    int nb_new_edges = 0; int nb_edges_ref = mesh.GetNbEdgesRef();
    Vector<bool> EdgeToAdd(mesh.GetNbEdges()); EdgeToAdd.Fill(false);
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int ref = mesh.BoundaryRef(i).GetReference();
	int type = this->GetCurveType(ref);
        if ((type != CURVE_FILE) && (type > 0))
	  {
	    // predefined curve
	    for (int j = 0; j < mesh.BoundaryRef(i).GetNbEdges(); j++)
	      {
		int num_edge = mesh.BoundaryRef(i).numEdge(j);
		if ((num_edge >= nb_edges_ref)&&(!EdgeToAdd(num_edge)))
		  {
		    // new edge
		    EdgeToAdd(num_edge) = true;
		    mesh.GetEdge(num_edge).SetReference(ref);
		    nb_new_edges++;
		  }
	      }
	  }
      }
    
    if (nb_new_edges > 0)
      {
	// we add the new referenced edges
	int nb_old = mesh.GetNbEdgesRef();
	mesh.ResizeEdgesRef(nb_new_edges + nb_old);
	for (int i = 0; i < mesh.GetNbEdges(); i++)
	  if (EdgeToAdd(i))
	    {
	      mesh.EdgeRef(nb_old) = mesh.GetEdge(i);
              nb_old++;
	    }
      }
    
    // trying to find new referenced faces
    int nb_new_faces = 0;
    Vector<bool> FaceToAdd(mesh.GetNbFaces()); FaceToAdd.Fill(false);
    for (int i = mesh.GetNbBoundaryRef(); i < mesh.GetNbBoundary(); i++)
      {
	bool curved_face = false;
	for (int j = 0; j < mesh.Boundary(i).GetNbEdges(); j++)
	  {
	    int num_edge = mesh.Boundary(i).numEdge(j);
	    int ref = mesh.GetEdge(num_edge).GetReference();
	    int type = this->GetCurveType(ref);
	    if (type > 0)
	      curved_face = true;
	  }
	
	if (curved_face)
	  {
	    FaceToAdd(i) = true;
	    nb_new_faces++;
	  }
      }
    
    if (nb_new_faces > 0)
      {
	// we add the new referenced faces
	int nb_old = mesh.GetNbBoundaryRef();
	mesh.ResizeBoundariesRef(nb_new_faces + nb_old);
	for (int i = 0; i < mesh.GetNbBoundary(); i++)
	  if (FaceToAdd(i))
	    {
	      mesh.BoundaryRef(nb_old) = mesh.Boundary(i);
              mesh.BoundaryRef(nb_old).SetReference(-2);
	      nb_old++;
	    }
      }
    
    // then marking elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	bool curved_element = false;
	for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	  {
	    int num_edge = mesh.Element(i).numEdge(j);
	    int ref = mesh.GetEdge(num_edge).GetReference();
	    int type = this->GetCurveType(ref);
	    if (type > 0)
	      curved_element = true;
	  }
	
	if (curved_element)
	  {
            mesh.Element(i).SetCurved();
            mesh.affine_element(i) = false;
          }
      }

    if ((nb_new_edges > 0) || (nb_new_faces > 0))
      mesh.FindConnectivity();

    // computing parameters (e.g. teta, phi for a sphere) for edges
    value_parameter_edge.Reallocate(mesh.GetNbEdgesRef());
    integer_parameter_edge.Reallocate(mesh.GetNbEdgesRef());
    for (int i = 0; i < mesh.GetNbEdgesRef(); i++)
      {
	int ref = mesh.EdgeRef(i).GetReference();
        if (ref > 0)
          {
            const VectReal_wp& param_ = this->Param_curve(ref);
            bool edge_valid = true;
            int n1 = mesh.EdgeRef(i).numVertex(0);
            int n2 = mesh.EdgeRef(i).numVertex(1);
            R_N Param1, Param2;
            if (n1 < 0 )
              edge_valid = false;
            else
              GetParameterValueOnPoint(ref, mesh.Vertex(n1), param_, Param1);
            
            if (n2 < 0 )
              edge_valid = false;
            else
              GetParameterValueOnPoint(ref, mesh.Vertex(n2), param_, Param2);
            
            if (edge_valid)
              {
                switch (mesh.GetCurveType(ref))
                  {
                  case CURVE_SPHERE :
                  case CURVE_CYLINDER :
                  case CURVE_CONIC :
                    {
                      value_parameter_edge(i)(0) = Param1(0);
                      value_parameter_edge(i)(1) = Param1(1);
                      value_parameter_edge(i)(2) = Param2(0);
                      
                      // we force angles phi/teta to be close
                      if (Param2(1) < Param1(1) - pi_wp)
                        Param2(1) += 2.0*pi_wp;
                      else if (Param2(1) > Param1(1) + pi_wp)
                        Param2(1) -= 2.0*pi_wp;
                      
                      value_parameter_edge(i)(3) = Param2(1);
                    }
                    break;
                  }
              }
          }
      }
    
    value_parameter_boundary.Reallocate(mesh.GetNbBoundaryRef());
    integer_parameter_boundary.Reallocate(mesh.GetNbBoundaryRef());
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int ref = mesh.BoundaryRef(i).GetReference();
        if (ref > 0)
          {
            const VectReal_wp& param_ = this->Param_curve(ref);
            int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
            IVect num(nb_vert); Vector<R_N> Pts(nb_vert), Param(nb_vert);
            bool face_valid = true;
            for (int j = 0; j < nb_vert; j++)
              {
                num(j) = mesh.BoundaryRef(i).numVertex(j);
                if (num(j) < 0)
                  face_valid = false;
                else
                  {
                    Pts(j) = mesh.Vertex(num(j));                
                    GetParameterValueOnPoint(ref, Pts(j), param_, Param(j));
                  }
              }
            
            if (face_valid)
              {
                switch (mesh.GetCurveType(ref))
                  {
                  case CURVE_SPHERE :
                  case CURVE_CYLINDER :
                  case CURVE_CONIC :
                    {
                      value_parameter_boundary(i).Reallocate(2*nb_vert);
                      for (int j = 0; j < nb_vert; j++)
                        {
                          if (j > 0)
                            {
                              // we force angles phi to be close
                              if (Param(j)(1) < Param(0)(1) - pi_wp)
                                Param(j)(1) += 2.0*pi_wp;
                              else if (Param(j)(1) > Param(0)(1) + pi_wp)
                                Param(j)(1) -= 2.0*pi_wp;
                            }
                          
                          value_parameter_boundary(i)(2*j) = Param(j)(0);
                          value_parameter_boundary(i)(2*j+1) = Param(j)(1);
                        }
                    }
                    break;
                  }
              }
          }
      }
  }
  
  
  //! returns the nodal points (real positions) of a face
  /*!
    \param[in] i referenced face number
    \param[out] PtFace nodal points of the face
    Here PtFace contains all the nodal points (vertices and points on edges are included)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::GetNodalPointsOnFace(int i, Vector<R_N>& PtFace) const
  {
    const Mesh<Dimension3, T>& mesh = this->GetLeafClass();
    int r = this->lob_curve.GetOrder();
    const VectReal_wp& xi = this->lob_curve.Points();
    if (mesh.BoundaryRef(i).IsTriangular())
      {
        PtFace.Reallocate(this->triangle_reference.GetNbPointsNodalElt());
        for (int j = 0; j < 3; j++)
          {
            int nv = mesh.BoundaryRef(i).numVertex(j);
            PtFace(j) = mesh.Vertex(nv);
          }
        
        for (int j = 0; j < 3; j++)
          {
            int ne = mesh.BoundaryRef(i).numEdge(j);
            if (ne >= mesh.GetNbEdgesRef())
              for (int k = 0; k < r-1; k++)
                PtFace(3+k+j*(r-1)) = (1.0-xi(k+1))*PtFace(j) + xi(k+1)*PtFace((j+1)%3);
            else
              {
                if (mesh.BoundaryRef(i).GetOrientationEdge(j))
                  for (int k = 0; k < r-1; k++)
                    PtFace(3+k+j*(r-1)) = this->PointsEdgeRef(ne, k);
                else
                  for (int k = 0; k < r-1; k++)
                    PtFace(3+k+j*(r-1)) = this->PointsEdgeRef(ne, r-2-k);
              }
          }
        
        for (int k = 0; k < this->PointsFaceRef(i).GetM(); k++)
          PtFace(3*r + k) = this->PointsFaceRef(i)(k);
      }
    else
      {
        PtFace.Reallocate(this->quadrangle_reference.GetNbPointsNodalElt());
        for (int j = 0; j < 4; j++)
          {
            int nv = mesh.BoundaryRef(i).numVertex(j);
            PtFace(j) = mesh.Vertex(nv);
          }
        
        for (int j = 0; j < 4; j++)
          {
            int ne = mesh.BoundaryRef(i).numEdge(j);
            if (ne >= mesh.GetNbEdgesRef())
              for (int k = 0; k < r-1; k++)
                PtFace(4+k+j*(r-1)) = (1.0-xi(k+1))*PtFace(j) + xi(k+1)*PtFace((j+1)%4);
            else
              {
                if (mesh.BoundaryRef(i).GetOrientationEdge(j))
                  for (int k = 0; k < r-1; k++)
                    PtFace(4+k+j*(r-1)) = this->PointsEdgeRef(ne, k);
                else
                  for (int k = 0; k < r-1; k++)
                    PtFace(4+k+j*(r-1)) = this->PointsEdgeRef(ne, r-2-k);
              }
          }
        
        for (int k = 0; k < this->PointsFaceRef(i).GetM(); k++)
          PtFace(4*r + k) = this->PointsFaceRef(i)(k);
      }
  }
  
  
  //! Modifies nodal points of edge i such that the edge is a straight edge
  template<class T>
  void MeshBoundaries<Dimension3, T>::SetStraightEdge(int i)
  {
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();
    R3 ptA = mesh_leaf.Vertex(mesh_leaf.EdgeRef(i).numVertex(0));
    R3 ptB = mesh_leaf.Vertex(mesh_leaf.EdgeRef(i).numVertex(1));
    T lambda; R3 ptC;
    for (int k = 0; k < this->lob_curve.GetOrder()-1; k++)
      {
        lambda = this->lob_curve.Points(k+1);
        ptC.Zero();
        Add(1.0-lambda, ptA, ptC);
        Add(lambda, ptB, ptC);
        this->PointsEdgeRef(i, k) = ptC;
      }
  }
  

  //! Modifies nodal points of face i such that the curvature comes only from edges
  /*!
    Here we consider that the expression of the curved face can be computed only
    with curved edges (e.g. with Gordon-Hall transformation for quadrilateral faces)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::SetStraightFace(int i)
  {
    ComputeInternalPointsFaceWithCurvedEdges(i);
  }
  
  
  //! Modifies nodal points of face i such that the curvature comes only from edges
  /*!
    Here we consider that the expression of the curved face can be computed only
    with curved edges (e.g. with Gordon-Hall transformation for quadrilateral faces)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::ComputeInternalPointsFaceWithCurvedEdges(int i)
  {
    Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();
    // here we consider that the face is curved because of surrounding edges
    int nb_vert = mesh_leaf.BoundaryRef(i).GetNbVertices();
    VectR3 s(nb_vert); SetPoints<Dimension3> PointsElem;
    for (int j = 0; j < nb_vert; j++)
      s(j) = mesh_leaf.Vertex(mesh_leaf.BoundaryRef(i).numVertex(j));
    
    int r = this->GetGeometryOrder();
    int Ntri = (r-1)*(r-2)/2; 
    int Nquad = (r-1)*(r-1);
    if (nb_vert == 3)
      {
        this->triangle_reference
          .FjElemNodalCurve(s, PointsElem, mesh_leaf, i, mesh_leaf.BoundaryRef(i));
        
        this->PointsFaceRef(i).Reallocate(Ntri);
        for (int k = 0; k < Ntri; k++)
          this->PointsFaceRef(i)(k) = PointsElem.GetPointNodal(3*r+k);
      }
    else
      {
        this->quadrangle_reference
          .FjElemNodalCurve(s, PointsElem, mesh_leaf, i, mesh_leaf.BoundaryRef(i));
        
        this->PointsFaceRef(i).Reallocate(Nquad);
        for (int k = 0; k < Nquad; k++)
          this->PointsFaceRef(i)(k) = PointsElem.GetPointNodal(4*r+k);
      }
  }
  
  
  //! computes the position of points of a curved edge
  /*!
    \param[in] i referenced edge number
    \param[in] ref reference of the surface containing the edge
    \param[in] xi local parameters on the edge (xi \in [0,1])
    \param[out] points points = F(xi)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::GetPointsOnCurvedEdge(int i, int ref, const VectReal_wp& xi, VectR3& points) const
  {
    const Mesh<Dimension3>& mesh = this->GetLeafClass();
    points.Reallocate(xi.GetM()); R_N point;
    switch (this->GetCurveType(ref))
      {
      case CURVE_SPHERE :
      case CURVE_CYLINDER :
      case CURVE_CONIC :
        {
          int n1 = mesh.EdgeRef(i).numVertex(0);
          int n2 = mesh.EdgeRef(i).numVertex(1);
          for (int k = 0; k < xi.GetM(); k++)
            {
              point.Zero();
              Add(1.0-xi(k), mesh.Vertex(n1), point);
              Add(xi(k), mesh.Vertex(n2), point);
              ProjectToCurve(point, ref);
              points(k) = point;
            }          
        }
        break;
      default :
        {
          // use of PointsEdgeRef
          int n1 = mesh.EdgeRef(i).numVertex(0);
          int n2 = mesh.EdgeRef(i).numVertex(1);
          int r = this->lob_curve.GetOrder();
          for (int k = 0; k < xi.GetM(); k++)
	    {
	      point.Zero();
	      Add(this->lob_curve.EvaluatePhi(0, xi(k)), mesh.Vertex(n1), point);
	      Add(this->lob_curve.EvaluatePhi(r, xi(k)), mesh.Vertex(n2), point);
	      for (int j = 1; j < r; j++)
		Add(this->lob_curve.EvaluatePhi(j, xi(k)), PointsEdgeRef(i, j-1), point);
	      
	      points(k) = point;
	    }	  
        }
      }
  }
  
  
  //! computes the position of points of a curved face
  /*!
    \param[in] i referenced face number
    \param[in] ref reference of the surface containing the edge
    \param[in] xi local parameters on the face (xi \in [0,1]^2 for quadrilaterals)
    \param[out] points points = F(xi)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::GetPointsOnCurvedFace(int i, int ref, const VectR2& xi, VectR3& points,
                          bool straight_face) const
  {
    const Mesh<Dimension3>& mesh = this->GetLeafClass();
    points.Reallocate(xi.GetM()); R_N point;
    if (mesh.BoundaryRef(i).IsTriangular())
      {
        switch (this->GetCurveType(ref))
          {
          case CURVE_SPHERE :
          case CURVE_CYLINDER :
          case CURVE_CONIC :
            {
              int n1 = mesh.BoundaryRef(i).numVertex(0);
              int n2 = mesh.BoundaryRef(i).numVertex(1);
              int n3 = mesh.BoundaryRef(i).numVertex(2);
              if (straight_face)
                {
                  for (int k = 0; k < xi.GetM(); k++)
                    {
                      point.Zero();
                      Add(1.0-xi(k)(0)-xi(k)(1), mesh.Vertex(n1), point);
                      Add(xi(k)(0), mesh.Vertex(n2), point);
                      Add(xi(k)(1), mesh.Vertex(n3), point);
                      ProjectToCurve(point, ref);
                      points(k) = point;
                    }
                  
                  return;
                }
              
              // Nielson transformation
              // contribution from vertices
              for (int k = 0; k < xi.GetM(); k++)
                {
                  points(k).Zero();
                  Add(-(1.0-xi(k)(0)-xi(k)(1))/3, mesh.Vertex(n1), points(k));
                  Add(-xi(k)(0)/3, mesh.Vertex(n2), points(k));
                  Add(-xi(k)(1)/3, mesh.Vertex(n3), points(k));
                }
              
              int r = this->lob_curve.GetOrder();
              VectR3 PtsEdge(r+1); VectReal_wp phi_loc(r+1);
              R3 ptM;

              // contribution from first edge
              int num_edge = mesh.BoundaryRef(i).numEdge(0);
              if (num_edge < mesh.GetNbEdgesRef())
                {
                  PtsEdge(0) = mesh.Vertex(n1);
                  PtsEdge(r) = mesh.Vertex(n2);
                  if (mesh.BoundaryRef(i).GetOrientationEdge(0))
                    for (int k = 1; k < r; k++)
                      PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, k-1);
                  else
                    for (int k = 1; k < r; k++)
                      PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, r-k-1);
                }
              
              for (int k = 0; k < xi.GetM(); k++)
                {
                  // part (2-x-2y)/3 e0(x)
                  Real_wp x = xi(k)(0), y = xi(k)(1);
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(x, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (1.0-x)*mesh.Vertex(n1) + x*mesh.Vertex(n2);
                  
                  Add((2.0-x-2*y)/3, ptM, points(k));
                  
                  // part (1+x-y)/3 e0(x+y)
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(x+y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (1.0-(x+y))*mesh.Vertex(n1) + (x+y)*mesh.Vertex(n2);
                  
                  Add((1.0 + x-y)/3, ptM, points(k));
                  
                  // part -x/3 e0(1-y)                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(1.0-y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = y*mesh.Vertex(n1) + (1.0-y)*mesh.Vertex(n2);

                  Add(-x/3, ptM, points(k));
                  
                  // part -(1-x-y)/3 e0(y)                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (1.0-y)*mesh.Vertex(n1) + y*mesh.Vertex(n2);
                  
                  Add(-(1.0-x-y)/3, ptM, points(k));
                }              

              // contribution from second edge
              num_edge = mesh.BoundaryRef(i).numEdge(1);
              if (num_edge < mesh.GetNbEdgesRef())
                {
                  PtsEdge(0) = mesh.Vertex(n2);
                  PtsEdge(r) = mesh.Vertex(n3);
                  if (mesh.BoundaryRef(i).GetOrientationEdge(1))
                    for (int k = 1; k < r; k++)
                      PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, k-1);
                  else
                    for (int k = 1; k < r; k++)
                      PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, r-k-1);
                }
              
              for (int k = 0; k < xi.GetM(); k++)
                {
                  // part (x+2y)/3 e1(1-x)
                  Real_wp x = xi(k)(0), y = xi(k)(1);
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(1.0-x, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = x*mesh.Vertex(n2) + (1.0-x)*mesh.Vertex(n3);
                  
                  Add((x+2*y)/3, ptM, points(k));
                  
                  // part (2x+y)/3 e1(y)
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (1.0-y)*mesh.Vertex(n2) + y*mesh.Vertex(n3);
                  
                  Add((2.0*x+y)/3, ptM, points(k));
                  
                  // part -y/3 e1(x+y)                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(x+y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (1.0-x-y)*mesh.Vertex(n2) + (x+y)*mesh.Vertex(n3);
                  
                  Add(-y/3, ptM, points(k));
                  
                  // part -x/3 e1(1-x-y)                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(1.0-x-y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (x+y)*mesh.Vertex(n2) + (1.0-x-y)*mesh.Vertex(n3);
                  
                  Add(-x/3, ptM, points(k));
                }              

              // contribution from third edge
              num_edge = mesh.BoundaryRef(i).numEdge(2);
              if (num_edge < mesh.GetNbEdgesRef())
                {
                  PtsEdge(0) = mesh.Vertex(n3);
                  PtsEdge(r) = mesh.Vertex(n1);
                  if (mesh.BoundaryRef(i).GetOrientationEdge(2))
                    for (int k = 1; k < r; k++)
                      PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, k-1);
                  else
                    for (int k = 1; k < r; k++)
                      PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, r-k-1);
                }
              
              for (int k = 0; k < xi.GetM(); k++)
                {
                  // part (2-2x-y)/3 e2(1-y)
                  Real_wp x = xi(k)(0), y = xi(k)(1);
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(1.0-y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = y*mesh.Vertex(n3) + (1.0-y)*mesh.Vertex(n1);
                  
                  Add((2.0-2.0*x-y)/3, ptM, points(k));
                  
                  // part (1+y-x)/3 e2(1-x-y)
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(1.0-x-y, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (x+y)*mesh.Vertex(n3) + (1.0-x-y)*mesh.Vertex(n1);
                  
                  Add((1.0+y-x)/3, ptM, points(k));
                  
                  // part -(1-x-y)/3 e2(1-x)                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(1.0-x, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = x*mesh.Vertex(n3) + (1.0-x)*mesh.Vertex(n1);
                  
                  Add(-(1.0-x-y)/3, ptM, points(k));
                  
                  // part -y/3 e2(x)                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      this->lob_curve.ComputeValuesPhiRef(x, phi_loc);
                      
                      ptM.Zero();
                      for (int j = 0; j <= r; j++)
                        ptM += phi_loc(j)*PtsEdge(j);
                    }
                  else
                    ptM = (1.0-x)*mesh.Vertex(n3) + x*mesh.Vertex(n1);
                  
                  Add(-y/3, ptM, points(k));
                }              

              // projection to surface
              for (int k = 0; k < xi.GetM(); k++)
                ProjectToCurve(points(k), ref);
            }
            break;
          default :
            {
              // use of PointsFaceRef
              VectR3 PtFace;
              this->GetNodalPointsOnFace(i, PtFace);
              VectReal_wp phi;
              for (int k = 0; k < xi.GetM(); k++)
                {
                  point.Zero();
                  this->triangle_reference.ComputeValuesPhiNodalRef(xi(k), phi);
                  for (int j = 0; j < PtFace.GetM(); j++)
                    Add(phi(j), PtFace(j), point);
                  
                  points(k) = point;
                }	  
            }
          }
      }
    else
      {
        switch (this->GetCurveType(ref))
          {
          case CURVE_SPHERE :
          case CURVE_CYLINDER :
          case CURVE_CONIC :
            {
              int n1 = mesh.BoundaryRef(i).numVertex(0);
              int n2 = mesh.BoundaryRef(i).numVertex(1);
              int n3 = mesh.BoundaryRef(i).numVertex(2);
              int n4 = mesh.BoundaryRef(i).numVertex(3);

              VectR3 s(4);
              s(0) = mesh.Vertex(n1);
              s(1) = mesh.Vertex(n2);
              s(2) = mesh.Vertex(n3);
              s(3) = mesh.Vertex(n4);
              
              if (straight_face)
                {
                  for (int k = 0; k < xi.GetM(); k++)
                    {
                      point.Zero();
                      Add((1.0-xi(k)(0))*(1.0-xi(k)(1)), mesh.Vertex(n1), point);
                      Add(xi(k)(0)*(1.0-xi(k)(1)), mesh.Vertex(n2), point);
                      Add(xi(k)(0)*xi(k)(1), mesh.Vertex(n3), point);
                      Add((1.0-xi(k)(0))*xi(k)(1), mesh.Vertex(n4), point);
                      ProjectToCurve(point, ref);
                      points(k) = point;
                    }
                  
                  return;
                }
              
              // we compute e0(x), e1(y), e2(1-x) and e3(1-y)
              int r = this->lob_curve.GetOrder();
              VectR3 PtsEdge(r+1); VectReal_wp phi_loc(r+1);
              TinyVector<VectR3, 4> proj_edge;              
              VectReal_wp xi_loc(xi.GetM());
              for (int ne_loc = 0; ne_loc < 4; ne_loc++)
                {
                  proj_edge(ne_loc).Reallocate(xi.GetM());
                  int num_edge = mesh.BoundaryRef(i).numEdge(ne_loc);
                  switch (ne_loc)
                    {
                    case 0 :
                      {
                        for (int k = 0; k < xi.GetM(); k++)
                          xi_loc(k) = xi(k)(0);
                      }
                      break;
                    case 1 :
                      {
                        for (int k = 0; k < xi.GetM(); k++)
                          xi_loc(k) = xi(k)(1);
                      }
                      break;
                    case 2 :
                      {
                        for (int k = 0; k < xi.GetM(); k++)
                          xi_loc(k) = 1.0 - xi(k)(0);
                      }
                      break;
                    case 3 :
                      {
                        for (int k = 0; k < xi.GetM(); k++)
                          xi_loc(k) = 1.0-xi(k)(1);                        
                      }
                      break;
                    }
                  
                  if (num_edge < mesh.GetNbEdgesRef())
                    {
                      PtsEdge(0) = s(ne_loc);
                      PtsEdge(r) = s((ne_loc+1)%4);
                      if (mesh.BoundaryRef(i).GetOrientationEdge(ne_loc))
                        for (int k = 1; k < r; k++)
                          PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, k-1);
                      else
                        for (int k = 1; k < r; k++)
                          PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, r-k-1);                      
                      
                      for (int k = 0; k < xi.GetM(); k++)
                        {
                          this->lob_curve.ComputeValuesPhiRef(xi_loc(k), phi_loc);
                          
                          proj_edge(ne_loc)(k).Zero();
                          for (int j = 0; j <= r; j++)
                            proj_edge(ne_loc)(k) += phi_loc(j)*PtsEdge(j);
                        }
                    }
                  else
                    for (int k = 0; k < xi.GetM(); k++)
                      proj_edge(ne_loc)(k) = (1.0-xi_loc(k))*s(ne_loc) + xi_loc(k)*s((ne_loc+1)%4);
                
                }
              
              for (int k = 0; k < xi.GetM(); k++)
                {
                  point.Zero();
                  // Gordon-Hall transformation
                  // Fi = (1-y) e0(x) + y e2(1-x) + x e1(y) + (1-x) e3(1-y)
                  //      - (1-x)(1-y) S0 - x(1-y) S1 - xy S2 - (1-x)y S3
                  point = (1.0 - xi(k)(1))*proj_edge(0)(k) + xi(k)(1)*proj_edge(2)(k)
                    + (1.0 - xi(k)(0))*proj_edge(3)(k) + xi(k)(0)*proj_edge(1)(k);

                  Add(-(1.0-xi(k)(0))*(1.0-xi(k)(1)), mesh.Vertex(n1), point);
                  Add(-xi(k)(0)*(1.0-xi(k)(1)), mesh.Vertex(n2), point);
                  Add(-xi(k)(0)*xi(k)(1), mesh.Vertex(n3), point);
                  Add(-(1.0-xi(k)(0))*xi(k)(1), mesh.Vertex(n4), point);

                  ProjectToCurve(point, ref);
                  points(k) = point;
                }
            }
            break;
          default :
            {
              // use of PointsFaceRef
              VectR3 PtFace;
              this->GetNodalPointsOnFace(i, PtFace);
              
              VectReal_wp phi;
              for (int k = 0; k < xi.GetM(); k++)
                {
                  point.Zero();
                  this->quadrangle_reference.ComputeValuesPhiNodalRef(xi(k), phi);
                  for (int j = 0; j < PtFace.GetM(); j++)
                    Add(phi(j), PtFace(j), point);
                  
                  points(k) = point;
                }	  
            }
          }
      }
  }
  
  
  //! reads type of curve
  /*!
    \param[in] ref references of all the curves with this type
    \param[in] description keyword to describe type of curve (CIRCLE, PARABOLE...)
    This methods modifies array Type_curve
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::ReadCurveType(const IVect& ref, const string& description)
  {
    if (!description.compare("SPHERE"))
      for (int i = 0; i < ref.GetM(); i++)
	this->SetCurveType(ref(i), CURVE_SPHERE);
    else if (!description.compare("CYLINDER"))
      for (int i = 0; i < ref.GetM(); i++)
	this->SetCurveType(ref(i), CURVE_CYLINDER);
    else if (!description.compare("CONE"))
      for (int i = 0; i < ref.GetM(); i++)
	this->SetCurveType(ref(i), CURVE_CONIC);
    else
      for (int i = 0; i < ref.GetM(); i++)
	this->SetCurveType(ref(i), NO_CURVE);
  }
  
  
  //! projection of ptA on a sphere
  /*!
    \param[in,out] ptA point to be projected
    \param[in] center center of the sphere
    \param[in] radius radius of the sphere
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::ProjectToSphere(R3& ptA, const R3& center, const T& radius)
  {
    Add(T(-1), center, ptA);
    T coef = radius/Norm2(ptA);
    Mlt(coef, ptA); Add(T(1), center, ptA);
  }
  
  
  //! projection of a point on a conic cylinder
  /*!
    \param[in,out] ptM point to be projected
    \param[in] ptC center of the cone
    \param[in] vdirector axis of the cone
    \param[in] theta angle of the cone
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::ProjectToConicCylinder(R3& ptM, const R3& ptC, const R3& vdirector, const T& theta)
  {
    // ptP is the projection of M to the axis of cylinder
    R3 CM, ptP; CM = ptM - ptC;
    T alpha = DotProd(CM, vdirector);
    ptP = ptC; Add(alpha, vdirector, ptP);
    
    // evaluating radius
    T rp(0);
    CM = ptC - ptP;
    rp = tan(theta)*Norm2(CM);
    
    // now projecting ptM to the cylinder/conic cylinder
    Add(T(-1), ptP, ptM); 
    if (Norm2(ptM) > R3::threshold)
      {
	T coef = rp/Norm2(ptM);
        Mlt(coef, ptM);
	Add(T(1), ptP, ptM);
      }
    else
      ptM = ptC;
  }
  
  
  //! projection of a point on a cylinder
  /*!
    \param[in,out] ptM point to be projected
    \param[in] ptC point on the axis
    \param[in] vdirector axis of the cylinder
    \param[in] radius radius of the cylinder
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::ProjectToCylinder(R3& ptM, const R3& ptC, const R3& vdirector, const T& radius)
  {
    // ptP is the projection of M to the axis of cylinder
    R3 CM, ptP; CM = ptM - ptC;
    T alpha = DotProd(CM, vdirector);
    ptP = ptC; Add(alpha, vdirector, ptP);
    
    // evaluating radius
    T rp(0);
    // cylinder
    rp = radius;
    
    // now projecting ptM to the cylinder/conic cylinder
    Add(T(-1), ptP, ptM); T coef = rp/Norm2(ptM);
    Mlt(coef, ptM); Add(T(1), ptP, ptM);
  }
  
  
  //! projection of a point on the curve
  /*!
    \param[in,out] ptM point to be projected
    \param[in] ref reference of the curved surface
  */
  template<class T>
  void MeshBoundaries<Dimension3, T>::ProjectToCurve(R3& ptM, int ref) const
  {
    if ((ref > 0)&&(this->GetCurveType(ref) >=1 ))
      {
	// curved face
	switch (this->GetCurveType(ref))
	  {
          case CURVE_CYLINDER :
	    {
              // a point A on the axis
              R3 pointA(this->Param_curve(ref)(0),
                        this->Param_curve(ref)(1),this->Param_curve(ref)(2));
              // unitary vector along axis
              R3 vec_u(this->Param_curve(ref)(3),
                       this->Param_curve(ref)(4),this->Param_curve(ref)(5));
              // radius of the cylinder
              T radius = this->Param_curve(ref)(6);
              
              ProjectToCylinder(ptM, pointA, vec_u, radius);
            }
            break;
          case CURVE_CONIC :
            {
              // center of the cone
              R3 pointA(this->Param_curve(ref)(0),
                        this->Param_curve(ref)(1),this->Param_curve(ref)(2));
              // unitary vector along axis
              R3 vec_u(this->Param_curve(ref)(3),
                       this->Param_curve(ref)(4),this->Param_curve(ref)(5));
              // angle
              T theta = pi_wp*this->Param_curve(ref)(6)/180;
              
              ProjectToConicCylinder(ptM, pointA, vec_u, theta);
            }
            break;
          case CURVE_SPHERE :
            {
              // center of the sphere
              R3 center(this->Param_curve(ref)(0),
                        this->Param_curve(ref)(1),this->Param_curve(ref)(2));
              // radius of the sphere
              T radius = this->Param_curve(ref)(3);
	    
              ProjectToSphere(ptM, center, radius);
            }
            break;
          }
      }
  }
  
  
  //! retrieves all the points belonging to a referenced surface, computes normale
  /*!
    \param[in] ref reference of the considered surface
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::GetAllPointsOnReference(int ref, VectR3& Points, VectR3& Normale) const
  {
    IVect ref_cond(this->GetLeafClass().GetNbReferences()+1);
    ref_cond.Fill();
    GetAllPointsOnReference(ref, ref_cond, Points, Normale);
  }
  
  
  //! retrieves all the points belonging to a referenced surface, computes normale
  /*!
    \param[in] ref surfaces selected will be such that ref_cond(ref_surface) = ref
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::GetAllPointsOnReference(int ref, const IVect& ref_cond,
                            VectR3& Points, VectR3& Normale) const
  {
    const Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();
    IVect VertexOnRef(mesh_leaf.GetNbVertices()), NbFacesRef(mesh_leaf.GetNbVertices());
    VertexOnRef.Fill(-1);
    NbFacesRef.Zero(); IVect Index_Vertices(mesh_leaf.GetNbVertices()); Index_Vertices.Fill(-1);
    
    int nb_vertex_ref = 0;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
        int ref_face = mesh_leaf.faces_ref(i).GetReference();
        if (ref_face >= 0)
          if (ref_cond(ref_face) == ref)
            {
              for (int j = 0; j < mesh_leaf.faces_ref(i).GetNbVertices(); j++)
                {
                  if (VertexOnRef(mesh_leaf.faces_ref(i).numVertex(j)) == -1)
                    {
                      VertexOnRef(mesh_leaf.faces_ref(i).numVertex(j)) = i;
                      nb_vertex_ref++;
                    }
                  NbFacesRef(mesh_leaf.faces_ref(i).numVertex(j))++;
                }
            }
      }
    
    if (nb_vertex_ref == 0)
      {
	Points.Clear();
	return;
      }

    Points.Reallocate(nb_vertex_ref);
    Normale.Reallocate(nb_vertex_ref);
    nb_vertex_ref = 0;
    VectR3 s(4); R3 vec_u, vec_v, vec_w;
    for (int i = 0; i < mesh_leaf.GetNbVertices(); i++)
      if (VertexOnRef(i) != -1)
	{
	  Normale(nb_vertex_ref).Fill(0);
	  Index_Vertices(i) = nb_vertex_ref;
	  Points(nb_vertex_ref++) = mesh_leaf.Vertices(i);
	}
    
    // normales are computed with a mean of normale on each face
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
        int ref_face = mesh_leaf.faces_ref(i).GetReference();
        if (ref_face >= 0)
          if (ref_cond(ref_face) == ref)
            {
              int nv = mesh_leaf.faces_ref(i).GetNbVertices();
              for (int j = 0; j < nv; j++)
                s(j) = mesh_leaf.Vertices(mesh_leaf.faces_ref(i).numVertex(j));
              
              vec_u = s(1) - s(0);
              vec_v = s(nv-1) - s(0);
              TimesProd(vec_u, vec_v, vec_w);
              Mlt(T(1)/Norm2(vec_w), vec_w);
              
              for (int j = 0; j < mesh_leaf.faces_ref(i).GetNbVertices(); j++)
                {
                  nb_vertex_ref = Index_Vertices(mesh_leaf.faces_ref(i).numVertex(j));
                  T coef = T(1)/T(NbFacesRef(mesh_leaf.faces_ref(i).numVertex(j)));
                  if (DotProd(vec_w, Normale(nb_vertex_ref)) < 0)
                    Add(-coef, vec_w, Normale(nb_vertex_ref));
                  else
                    Add(coef, vec_w, Normale(nb_vertex_ref));
                }
            }
      }
  } 
  
  //! selects some points of the referenced surface ref
  /*!
    \param[in] ref reference of the considered surface
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
    \param[in] nb_points number of points to select
    This methods selects points so that points are enough distant and not coplanars
  */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::GetDistantPointsOnReference(int ref, VectR3& Points, VectR3& Normale, int nb_points) const
  {
    IVect ref_cond(this->GetLeafClass().GetNbReferences()+1);
    ref_cond.Fill();
    GetDistantPointsOnReference(ref, ref_cond, Points, Normale, nb_points);
  }
  

  //! selects some points of the referenced surface ref
  /*!
    \param[in] ref selected surfaces will be such that ref_cond(ref_surface) == ref
    \param[out] Points list of points of the considered surface
    \param[out] Normale normale on each point
    \param[in] nb_points number of points to select
    This methods selects points so that points are enough distant and not coplanars
  */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::GetDistantPointsOnReference(int ref, const IVect& ref_cond,
                                VectR3& Points, VectR3& Normale, int nb_points) const
  {
    GetAllPointsOnReference(ref, ref_cond, Points, Normale);
    if ((nb_points <= 0)||(Points.GetM() < nb_points))
      {
	Points.Clear();
	return;
      }
    
    IVect num(nb_points);
    // first point
    num(0) = 0;
    
    VectR3 NewPoints(nb_points); T c(0); R3 u, v, normale;
    VectR3 NewNormale(nb_points);
    NewPoints(0) = Points(num(0));
    NewNormale(0) = Normale(num(0));
    // loop for choosing points
    for (int current = 1; current < nb_points; current++)
      {
	T dist_max = 0.0; T dist_vert = 0.0;
	int imax = -1;
	for (int i = 0; i < Points.GetM(); i++)
	  {
	    dist_vert = Points(i).Distance(Points(num(0)));
	    for (int j = 0; j < current; j++)
	      dist_vert = min(dist_vert, Points(i).Distance(Points(num(j))));
	    	    
	    if (current >= 3)
	      {
		// we want any point that is not coplanar with three others
		for (int j = 0; j < current; j++)
		  for (int k = (j+1); k < current; k++)
		    for (int l = (k+1); l < current; l++)
		      {
			u = NewPoints(j) - NewPoints(k);
			v = NewPoints(j) - NewPoints(l);
			TimesProd(u, v, normale);
			c = DotProd(normale, NewPoints(j));
			dist_vert = min(dist_vert, T(abs(T(DotProd(normale, Points(i))-c))));
		      }
	      }
	    
	    if (dist_vert >= dist_max)
	      {
		imax = i;
		dist_max = dist_vert;
	      }	
	  }
	// we found a good point
	num(current) = imax;
	NewPoints(current) = Points(imax);
	NewNormale(current) = Normale(imax);
      }
	    
    Points.Clear(); Points.Copy(NewPoints);    
    Normale.Clear(); Normale.Copy(NewNormale);    
  }
  

  //! finds parameters of curved surfaces
  /*!
    For example, if a curved surface is a sphere, this method will find the center and
    radius of the sphere.
   */  
  template<class T>
  void MeshBoundaries<Dimension3, T>::FindParametersCurve()
  {
    Mesh<Dimension3>& mesh_leaf = this->GetLeafClass();
    // estimation of parameters if not given by the user
    for (int i = 1; i <= this->GetNbReferences(); i++)
      {
	if ((this->GetCurveType(i) > 0)&&(!this->Param_known(i)))
	  {
	    if (this->GetCurveType(i) == this->CURVE_FILE)
	      {
		// nothing to do
	      }
	    else if (this->GetCurveType(i) == this->CURVE_SPHERE)
	      {
		R3 center; T radius;
		FindParametersSphere(i, center, radius);
		//DISP(i); DISP(center); DISP(radius);
                if (mesh_leaf.print_level >= 2)
                  cout<<"Curve "<<i<<" is a sphere of center "
                      <<center<<"  and radius  "<<radius<<endl;
		
		this->Param_curve(i).Reallocate(4);
		this->Param_curve(i)(0) = center(0);
		this->Param_curve(i)(1) = center(1);
		this->Param_curve(i)(2) = center(2);
		this->Param_curve(i)(3) = radius;
		this->Param_known(i) = true;
		
	      }
	    else if (this->GetCurveType(i) == CURVE_CYLINDER)
	      {
		R3 ptC, vdirector; T radius;
		FindParametersCylinder(i, ptC, vdirector, radius);
		
                if (mesh_leaf.print_level >= 2)
                  cout<<"Curve "<<i<<" is a cylinder of axis "<<vdirector
                      <<" with the point "<<ptC<<" and radius  "<<radius<<endl;
		
		this->Param_curve(i).Reallocate(7);
		this->Param_curve(i)(0) = ptC(0);
		this->Param_curve(i)(1) = ptC(1);
		this->Param_curve(i)(2) = ptC(2);
		this->Param_curve(i)(3) = vdirector(0);
		this->Param_curve(i)(4) = vdirector(1);
		this->Param_curve(i)(5) = vdirector(2);
		this->Param_curve(i)(6) = radius;
		this->Param_known(i) = true;
		
	      }
	    else if (this->GetCurveType(i) == CURVE_CONIC)
	      {
		R3 ptC, vdirector; T theta(0);
		FindParametersConic(i, ptC, vdirector, theta);
		
                if (mesh_leaf.print_level >= 2)
                  cout<<"Curve "<<i<<" is a cone of axis "<<vdirector
                      <<" and center "<<ptC<<" and angle  "<<theta*180/pi_wp<<endl;
		
		this->Param_curve(i).Reallocate(7);
		this->Param_curve(i)(0) = ptC(0);
		this->Param_curve(i)(1) = ptC(1);
		this->Param_curve(i)(2) = ptC(2);
		this->Param_curve(i)(3) = vdirector(0);
		this->Param_curve(i)(4) = vdirector(1);
		this->Param_curve(i)(5) = vdirector(2);
		this->Param_curve(i)(6) = theta*180/pi_wp;
		this->Param_known(i) = true;
	      }
	    else
	      {
		cout<<"Case not treated "<<endl;
		abort();
	      }
	  }
      }
  }
  
  
  //! the surface ref is supposed to be a plane ax + by + cz + d = 0, 
  //! parameters a,b,c and d are estimated
  /*!
    \param[in] ref reference of the surface
    \param[out] coef_plane coefficients a, b, c
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::FindParametersPlane(int ref, Vector<T>& coef_plane) const
  {
    VectR3 Points, Normale;
    int nb_points = 3;
    GetDistantPointsOnReference(ref, Points, Normale, nb_points);
    if (Points.GetM() < nb_points)
      return;
    
    coef_plane.Reallocate(4);
    coef_plane(0) = Normale(0)(0); coef_plane(1) = Normale(0)(1); coef_plane(2) = Normale(0)(2);
    coef_plane(3) = -DotProd(Normale(0), Points(0));
  }
  
  
  //! estimation of the center and radius of the "potential" sphere
  /*!
    \param[in] ref reference of the surface
    \param[out] center estimation of the center of the sphere
    \param[out] radius estimation of the radius of the sphere
    \return 0 if successful
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::FindParametersSphere(int ref, R3& center, T& radius) const
  {
    center.Fill(0); radius = 0.0;
    VectR3 Points, Normale; R3 normale, middle, rhs; T c;
    Matrix3_3 mat, inv_mat;
    int nb_points = 4;
    GetDistantPointsOnReference(ref, Points, Normale, nb_points);
    if (Points.GetM() < nb_points)
      {
	cout<<"Not enough points to find parameters of the sphere "<< endl;
	cout<<"Reference of the surface " << ref;
	return;
      }
    
    // equation of perpendicular bisector plane of the edge AB
    // n_x * x + n_y * y + n_z * z  = c
    // where (n_x, n_y, n_z) is the normale to the plane (oriented along AB)
    // c = n \cdot middle(A,B) 
    
    // A = Points(0) B = Points(1)
    normale = Points(0) - Points(1);
    Add(Points(0), Points(1), middle); Mlt(0.5, middle);
    c = DotProd(normale, middle);
    mat(0,0) = normale(0); mat(0,1) = normale(1); mat(0,2) = normale(2); rhs(0) = c;
    
    // A = Points(0) B = Points(2)
    normale = Points(0) - Points(2);
    Add(Points(0), Points(2), middle); Mlt(0.5, middle);
    c = DotProd(normale, middle);
    mat(1,0) = normale(0); mat(1,1) = normale(1); mat(1,2) = normale(2); rhs(1) = c;
    
    // A = Points(0) B = Points(3)
    normale = Points(0) - Points(3);
    Add(Points(0), Points(3), middle); Mlt(0.5, middle);
    c = DotProd(normale, middle);
    mat(2,0) = normale(0); mat(2,1) = normale(1); mat(2,2) = normale(2); rhs(2) = c;
    
    // center of the sphere is the intersection of the three planes
    T det = Det(mat);
    if (abs(det) < 1e-15)
      {
	cout<<"Singular matrix for detection of parameters of the sphere "<<endl;
	cout<<"Four points are perhaps coplanar "<<endl;
	cout<<"Points used to determine sphere "<<Points<<endl;
	cout<<"Linear system to be inversed "<<mat<<endl;
	cout<<"Right hand side "<<rhs<<endl;
	throw WrongMesh("MeshBoundaries<Dimension3, T>::FindParametersSphere",
			"Singular matrix");
      }
    
    GetInverse(mat, inv_mat);
    Mlt(inv_mat, rhs, center);
    
    // radius is trivial
    radius = center.Distance(Points(0));
  }
  
  //! estimation of axis and angle of conic cylinder
  /*!
    \param[in] ref reference of the surface
    \param[out] ptC point on the axis
    \param[out] vdirector axis
    \param[out] theta angle of the cone
    \warning this methods fails often, be careful
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::FindParametersConic(int ref, R3& ptC, R3& vdirector, T& theta) const
  {
    ConicParam<T> fct;
    Vector<T> parameters;
    FindParametersGeneric(ref, parameters, fct);
    ptC(0) = parameters(0); ptC(1) = parameters(1); ptC(2) = parameters(2);
    vdirector(0) = parameters(3); vdirector(1) = parameters(4); vdirector(2) = parameters(5);
    theta = atan(parameters(6));
  }
  
  
  //! estimation of axis and angle of cylinder
  /*!
    \param[in] ref reference of the surface
    \param[out] ptC point on the axis
    \param[out] vdirector axis
    \param[out] radius radius of cylinder
    \warning this methods fails often, be careful
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>
  ::FindParametersCylinder(int ref, R3& ptC, R3& vdirector, T& radius) const
  {
    CylinderParam<T> fct;
    Vector<T> parameters;
    FindParametersGeneric(ref, parameters, fct);
    ptC(0) = parameters(0); ptC(1) = parameters(1); ptC(2) = parameters(2);
    vdirector(0) = parameters(3); vdirector(1) = parameters(4); vdirector(2) = parameters(5);
    radius = parameters(6);
  }
  
    
  //! generic method to find the parameters of a curved surface
  /*!
    \param[in] ref reference of the surface
    \param[out] xsol found parameters of the surface
    \param[in] fct given description of the curved surface
   */
  template<class T> 
  void MeshBoundaries<Dimension3, T>
  ::FindParametersGeneric(int ref, Vector<T>& xsol, SurfaceParamMinimizer<T>& fct) const
  {
    GetDistantPointsOnReference(ref, fct.Points, fct.Normales, fct.nb_points);
    if (fct.Points.GetM() < fct.nb_points)
      {
	cout<<"Not enough points to find parameters of the curved surface "<< endl;
	cout<<"Reference of the surface " << ref;
	return;
      }
    
    Real_wp residual;

    // we use Levenberg-Marquardt algorithm
    residual = Montjoie::SolveLeastSquaresLvm(fct, xsol, 1e2*epsilon_machine);
    
    if (residual > 1e-10)
      {
        cout << "Failed to find correct parameters" << endl;
        cout << "Parameters found = " << xsol << endl;
        cout << "Residual = " << residual << endl;
      }
  }
  
  
  /*****************
   * Other methods *
   *****************/
  
  
  //! clearings informations linking faces and edges
  template<class T>
  void MeshBoundaries<Dimension3, T>::ClearConnectivity()
  {
    MeshBoundaries_Base<Dimension3, T>::ClearConnectivity();
  }
  
  
  //! extraction of surfacic mesh from a volumic mesh
  /*!
    \param[in] ref reference of the faces to be extracted
    \param[out] mesh_surf surfacic mesh extracted
    \param[in] ref_cond if true, ref is related to a boundary condition
    (instead of a mesh reference)
   */
  template<class T>
  void MeshBoundaries<Dimension3, T>::GetBoundaryMesh(int ref, SurfacicMesh<Dimension3>& mesh,
                                                      const IVect& ref_cond, int ref_domain) const
  {
    GetBoundaryMesh(ref, mesh, mesh.ListeBoundaries, mesh.ListeVertices,
                    mesh.NumElement, mesh.NumLocalBoundary, ref_cond, ref_domain);
  }
  
  
  //! extraction of surfacic mesh from a volumic mesh
  /*!
    \param[in] ref reference of the faces to be extracted
    \param[out] mesh_surf surfacic mesh extracted
    \param[out] Index_FaceSurf_to_FaceRef numbers of referenced faces extracted
    \param[out] Liste_Vertices numbers of vertices extracted
    \param[out] NumElement element number for each face
    \param[out] NumLocalBoundary local position on the element
    \param[in] ref_cond extracts all the refence such that ref_cond(i) = ref
    \param[in] ref_domain if > 0, we extract only surface mesh in contact with elements whose
                          volume reference is ref_domain
   */
  template<class T> template<class Mesh3D_>
  void MeshBoundaries<Dimension3, T>::
  GetBoundaryMesh(int ref, Mesh3D_& mesh_surf, IVect& Index_FaceSurf_to_FaceRef,
                  IVect& Liste_Vertices, IVect& NumElement, IVect& NumLocalBoundary,
                  const IVect& ref_cond, int ref_domain) const
  {
    const Mesh<Dimension3, T>& mesh_leaf = this->GetLeafClass();    

    mesh_surf.Clear();
    mesh_surf.SetGeometryOrder(mesh_leaf.GetGeometryOrder());
    mesh_surf.Type_curve = this->Type_curve; mesh_surf.Cond_curve = this->Cond_curve;
    mesh_surf.Body_curve = this->Body_curve; mesh_surf.Param_curve = this->Param_curve;
    mesh_surf.Param_known = this->Param_known;

    // couting how many faces in the boundary mesh
    int nf = 0;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      {
        int ref_b = mesh_leaf.BoundaryRef(i).GetReference();
        if (ref_b > 0)
          if (ref_cond(ref_b) == ref)
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
      }
    
    // allocating arrays for the surfacic mesh
    mesh_surf.ReallocateBoundariesRef(nf);
    Index_FaceSurf_to_FaceRef.Reallocate(nf);
    Liste_Vertices.Reallocate(4*nf);
    NumElement.Reallocate(nf);
    NumLocalBoundary.Reallocate(nf);
    Vector<bool> EdgeUsed(mesh_leaf.GetNbEdgesRef());
    EdgeUsed.Fill(false);
    int ind = 0; nf = 0;
    // loop over all faces of the initial 3-D mesh
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      if (mesh_leaf.faces_ref(i).GetReference() > 0)
        {
          int ref_face = mesh_leaf.faces_ref(i).GetReference();
          ref_face = ref_cond(ref_face);
          if (ref_face == ref)
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
                  // the face is a face with the good boundary condition (array ref_cond)
                  mesh_surf.faces_ref(nf) = mesh_leaf.faces_ref(i);
                  Index_FaceSurf_to_FaceRef(nf) = i;
                  // adding the vertices to ListeVertices
                  for (int j = 0; j < mesh_leaf.faces_ref(i).GetNbVertices(); j++)
                    {
                      Liste_Vertices(ind) = mesh_leaf.faces_ref(i).numVertex(j);
                      int ne = mesh_leaf.BoundaryRef(i).numEdge(j);
                      if (ne < mesh_leaf.GetNbEdgesRef())
                        EdgeUsed(ne) = true;
                      
                      ind ++;
                    }
                  
                  // adding the face to NumElement
                  if (mesh_leaf.GetNbElt() > 0)
                    {
		      NumElement(nf) = num_elem;
		      NumLocalBoundary(nf) = mesh_leaf.elements(num_elem).
			GetPositionBoundary(i);
		    }
		  
                  nf++;
                }
            }
        }
    
    // we sort list of vertices
    Assemble(ind, Liste_Vertices);
    Liste_Vertices.Resize(ind);
    mesh_surf.ReallocateVertices(ind);
    IVect NewNumber_Vertices(mesh_leaf.GetNbVertices());
    for (int i = 0; i < ind; i++)
      {
	mesh_surf.Vertex(i) = mesh_leaf.Vertices(Liste_Vertices(i));
	NewNumber_Vertices(Liste_Vertices(i)) = i;
      }
    
    // now the new numbering of vertices modifies the connectivity of faces
    IVect num(4);
    for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
      {
        int ref_face = mesh_surf.faces_ref(i).GetReference();
        
        for (int j = 0; j < mesh_surf.faces_ref(i).GetNbVertices(); j++)
          num(j) = NewNumber_Vertices(mesh_surf.faces_ref(i).numVertex(j));
        
        if (mesh_surf.faces_ref(i).IsTriangular())
          mesh_surf.faces_ref(i).InitTriangular(num(0),num(1),num(2),ref_face);
        else
          mesh_surf.faces_ref(i).InitQuadrangular(num(0),num(1),num(2),num(3),ref_face);
        
        mesh_surf.PointsFaceRef(i) = mesh_leaf.PointsFaceRef(Index_FaceSurf_to_FaceRef(i));
      }
    
    // extracting referenced edges as well
    int nb_edges_ref = 0;
    for (int i = 0; i < EdgeUsed.GetM(); i++)
      if (EdgeUsed(i))
        nb_edges_ref++;
    
    mesh_surf.ReallocateEdgesRef(nb_edges_ref);
    nb_edges_ref = 0;
    for (int i = 0; i < EdgeUsed.GetM(); i++)
      if (EdgeUsed(i))
        {
          int n1 = mesh_leaf.EdgeRef(i).numVertex(0);
          int n2 = mesh_leaf.EdgeRef(i).numVertex(1);
          int ref = mesh_leaf.EdgeRef(i).GetReference();
          mesh_surf.EdgeRef(nb_edges_ref)
            .Init(NewNumber_Vertices(n1), NewNumber_Vertices(n2), ref);
          
          for (int k = 0; k < mesh_leaf.PointsEdgeRef.GetN(); k++)
            mesh_surf.PointsEdgeRef(nb_edges_ref, k) = mesh_leaf.PointsEdgeRef(i, k);
          
          nb_edges_ref++;
        }
    
    // last step (useful when a modal source is computed)
    // we find reference edges (i.e. isolated edges)
    // on the initial 3-D, these edges are not isolated, so we use the reference of adjacent face
    mesh_surf.FindConnectivity();
    int nb_old = nb_edges_ref;
    for (int i = mesh_surf.GetNbEdgesRef(); i < mesh_surf.GetNbEdges(); i++)
      if (mesh_surf.edges(i).GetNbFaces() == 1)
	nb_edges_ref++;
    
    if (nb_edges_ref > nb_old)
      {
        mesh_surf.ResizeEdgesRef(nb_edges_ref);
        nb_edges_ref = nb_old;
        for (int i = nb_old; i < mesh_surf.GetNbEdges(); i++)
	  if (mesh_surf.edges(i).GetNbFaces() == 1)
	    {
	      int n0 = mesh_surf.edges(i).numVertex(0);
	      int n1 = mesh_surf.edges(i).numVertex(1);

	      // which reference ?
	      nf = mesh_surf.edges(i).numFace(0);
	      int ref = mesh_surf.faces_ref(nf).GetReference();

	      // edge number and face number on the initial 3-D mesh
	      int num_loc = mesh_surf.faces_ref(nf).GetPositionBoundary(i);
	      int nf2 = Index_FaceSurf_to_FaceRef(nf);
	      int ne2 = mesh_leaf.faces_ref(nf2).numEdge(num_loc);
	      for (int m = 0; m < mesh_leaf.edges(ne2).GetNbFaces(); m++)
		{
		  int nfb = mesh_leaf.edges(ne2).numFace(m);
		  if ((nfb != nf2) && (nfb <= mesh_leaf.GetNbBoundaryRef())
		      &&(mesh_leaf.faces(nfb).GetReference() > 0))
		    ref = mesh_leaf.faces(nfb).GetReference();
		}
	      
	      mesh_surf.edges_ref(nb_edges_ref).Init(n0, n1, ref);
	      nb_edges_ref++;
	    }
	
	mesh_surf.FindConnectivity();
      }

    mesh_surf.ProjectPointsOnCurves();
  }
  

  /*************************************
   * IntersectionCurvedEdgeFace_System *
   *************************************/

  
  //! constructor
  /*!
    mesh0 : working mesh
    nf : face number where we search the intersection
    num_loc : local edge number on the face where the intersection will occur
    pos : local position of the point on the edge num_loc
    ref : reference of the face
    pt0 : center of the sphere
    radius : radius of the sphere
   */
  template<class T>
  IntersectionCurvedEdgeFace_System<T>
  ::IntersectionCurvedEdgeFace_System(const Mesh<Dimension3>& mesh0, int nf, int num_loc, const Real_wp& pos,
                                      int ref, const R3& pt0, const Real_wp& r)
    : mesh(mesh0), num_face1(nf), num_edge_loc1(num_loc), pos_edge(pos), ref1(ref), center(pt0), radius(r)
  {
    num_edge1 = -1;
    pt_on_edge = false;
  }


  template<class T>
  IntersectionCurvedEdgeFace_System<T>
  ::IntersectionCurvedEdgeFace_System(const Mesh<Dimension3>& mesh0, int ne,
                                      int ref, const R3& pt0, const Real_wp& r)
    : mesh(mesh0), num_edge1(ne), ref1(ref), center(pt0), radius(r)
  {
    num_face1 = -1; num_edge_loc1 = -1;
    pt_on_edge = true;
  }

  
  //! returns the point on the face with local abscissa x
  template<class T>
  R3 IntersectionCurvedEdgeFace_System<T>::GetPointOnSurface(const Real_wp& x) const
  {
    VectR3 points(1);
    if (pt_on_edge)
      {
        VectReal_wp xi_edge(1); xi_edge(0) = x;
        mesh.GetPointsOnCurvedEdge(num_edge1, ref1, xi_edge, points);
      }
    else
      {
        VectR2 xi_face(1); 
        bool tri = (mesh.Boundary(num_face1).GetNbVertices() == 3);
        if (tri)
          {
            // 0 : on balaie la ligne verticale x_loc = pos_edge => pt = (pos_edge, 0)
            //    quand x = 0, on est bien sur le point de l'arete 0 : (pos_edge, 0)
            // 1 : on balaie la ligne horizontale y_loc = pos_edge => pt = (1-x-pos_edge, pos_edge)
            //    quand x = 0, on est bien sur le point de l'arete 1 : (1-pos_edge, pos_edge)
            // 2 : on balaie la ligne diagonale x_loc + y_loc = 1- pos_edge => pt = (x, 1-pos_edge-x)
            //   quand x = 0, on est bien sur le point de l'arete 2 : (0, 1-pos_edge)
            switch (num_edge_loc1)
              {
              case 0: xi_face(0).Init(pos_edge, x); break;
              case 1: xi_face(0).Init(T(1) - x - pos_edge, pos_edge); break;
              case 2: xi_face(0).Init(x, Real_wp(1)-pos_edge-x); break;
              } 
          }
        else
          {
            // 0 : on balaie la ligne verticale x_loc = pos_edge => pt = (pos_edge, x)
            //     quand x = 0, on est bien sur le point de l'arete 0 : (pos_edge, 0)
            // 1 : on balaie la ligne horizontale y_loc = pos_edge => pt = (1-x, pos_edge)
            //     quand x = 0, on est bien sur le point de l'arete 1 (1, pos_edge)
            // 2 : on balaie la ligne verticale x_loc = 1-pos_edge => pt = (1-pos_edge, 1-x)
            //     quand x = 0, on est bien sur le point de l'arete 2 : (1-pos_edge, 1)
            // 3 : on balaie la ligne horizontale y_loc = 1-pos_edge => pt = (x, 1-pos_edge)
            //     quand x = 0, on est bien sur le point de l'arete 3 : (0, 1-pos_edge)
            switch (num_edge_loc1)
              {
              case 0: xi_face(0).Init(pos_edge, x); break; 
              case 1: xi_face(0).Init(T(1) - x, pos_edge); break;
              case 2: xi_face(0).Init(Real_wp(1)-pos_edge, T(1) - x); break;
              case 3: xi_face(0).Init(x, Real_wp(1)-pos_edge); break;
              }
          }
        
        mesh.GetPointsOnCurvedFace(num_face1, ref1, xi_face, points, true);
      }

    return points(0);
  }
  

  //! evaluates the distance between a point of the curved edge and the sphere
  template<class T>
  Real_wp IntersectionCurvedEdgeFace_System<T>::EvaluateFunction(const Real_wp& x) const
  { 
    R3 point = GetPointOnSurface(x);
    Real_wp f = Norm2(point-center) - radius;
    return f;
  }

  
  //! evalutes the derivative of EvaluateFunction
  template<class T>
  Real_wp IntersectionCurvedEdgeFace_System<T>::EvaluateDerivative(const Real_wp& x) const
  {
    // numerical derivative
    Real_wp h = 1e-6;
    Real_wp fp = EvaluateFunction(x+h);
    Real_wp fm = EvaluateFunction(x-h);
    return (fp - fm) / (2.0*h);
  }
  
    
  /*****************
   * CylinderParam *
   *****************/
  
  
  //! default constructor
  template<class T>
  CylinderParam<T>::CylinderParam()
  {
    this->n = 7;
    this->nb_points = 12;
    
    this->m_ = 1 + this->nb_points;
    this->n_ = this->n;
  }
  
  
  //! evaluation of the global "distance" of points to the cylinder
  /*!
    \param[in] x parameters of the cylinder (radius, axis ...)
    \param[in] param set of points
    \param[out] feval "distance" of points to the cylinder
    If all the points belong to the cylinder, the distance will be equal to 0
   */
  template<class T>
  void CylinderParam<T>::EvaluateFunction(const Vector<T>& x, T& feval)
  {
    R3 ptC, u, pt_current, CM;
    // point of the cylinder
    ptC.Init(x(0), x(1), x(2));
    // axis of the cylinder
    u.Init(x(3), x(4), x(5));
    T radius = x(6), r2, alpha, norme_u2, sq_norme_u2, evalg;
    norme_u2 = DotProd(u,u); r2 = radius*radius;
    sq_norme_u2 = norme_u2*norme_u2;
    feval = 0;
    // loop on points of the cylinder
    for (int i = 0; i < this->nb_points; i++)
      {
	pt_current = this->Points(i);
	
	CM = pt_current - ptC;
	alpha = DotProd(CM,u);
	evalg = DotProd(CM,CM)*sq_norme_u2 - alpha*alpha*norme_u2 - r2*sq_norme_u2;
	feval += evalg*evalg;
      }
    
    feval += (1.0-norme_u2)*(1.0-norme_u2);
  }
  
  
  //! gradient of the global "distance" of points to the cylinder
  /*!
    \param[in] x parameters of the cylinder (radius, axis ...)
    \param[in] param set of points
    \param[out] feval "distance" of points to the cylinder
    \param[out] gradient  gradient of "distance" of points to the cylinder
    If all the points belong to the cylinder, the distance will be equal to 0
    The derivation is done with respect to the parameters of the cylinder
    (radius, axis, point of the axis)
   */
  template<class T>
  void CylinderParam<T>
  ::EvaluateFunctionGradient(const Vector<T>& x, T& feval, Vector<T>& gradient)
  {
    R3 ptC, u, pt_current, CM;
    // point of the cylinder
    ptC.Init(x(0), x(1), x(2));
    // axis of the cylinder
    u.Init(x(3), x(4), x(5));
    T radius = x(6), r2, alpha, norme_u2, sq_norme_u2, evalg;
    norme_u2 = DotProd(u,u); r2 = radius*radius;
    sq_norme_u2 = norme_u2*norme_u2;
    feval = 0; gradient.Fill(0);
    // loop on points of the cylinder
    for (int i = 0; i < this->nb_points; i++)
      {
	pt_current = this->Points(i);
	
	CM = pt_current - ptC;
	alpha = DotProd(CM,u);
	evalg = DotProd(CM,CM)*sq_norme_u2 - alpha*alpha*norme_u2 - r2*sq_norme_u2;
	feval += evalg*evalg;
	
	for (int k = 0; k < 3; k++)
	  {
	    // derivative along C
	    gradient(k) += 2*evalg*(-2.0*CM(k)*sq_norme_u2 + 2.0*alpha*norme_u2*u(k));
	    // derivative along u
	    gradient(k+3) += 2*evalg*(4.0*DotProd(CM,CM)*norme_u2*u(k) - 2.0*alpha*alpha*u(k)
                                      - 2.0*alpha*CM(k)*norme_u2 - 4.0*r2*norme_u2*u(k));
	  }
	
	gradient(6) -= 4.0*evalg*radius*sq_norme_u2;
      }
    
    feval += (1.0-norme_u2)*(1.0-norme_u2);
    for (int k = 0; k < 3; k++)
      gradient(k+3) -= 4.0*(1.0-norme_u2)*u(k);
  }
  
  
  //! we find an initial guess for radius, axis, and point on the axis of a cylinder
  /*!
    \param[out] param found guess
   */
  template<class T>
  void CylinderParam<T>::FindInitGuess(Vector<T>& param)
  {
    param.Reallocate(this->n);
    param.Fill(T(0));
    // direction of axis is estimated by mean of cross products of normals
    R3 ptC, vdirector, vec_u, vec_v, vec_w;
    T sin_teta; vdirector.Fill(0); int nb_samples = 0;
    int imax = -1, jmax = -1; T angle_max(0), sum_coef(0);
    for (int i = 0; i < this->nb_points; i++)
      for (int j = (i+1); j < this->nb_points; j++)
	{
	  vec_u = this->Normales(i);
	  vec_v = this->Normales(j);
	  TimesProd(vec_u, vec_v, vec_w);
	  sin_teta = Norm2(vec_w);
	  T coef = abs(sin_teta);
	  if (coef > 1e-5)
	    {
	      Mlt(T(1)/sin_teta, vec_w);
	      if (DotProd(vdirector, vec_w) < 0)
		Add(T(-coef), vec_w, vdirector);
	      else
		Add(T(coef), vec_w, vdirector);
	      
	      sum_coef += coef;
	      nb_samples++;
	      if (sin_teta > angle_max)
		{
		  angle_max = sin_teta;
		  imax = i; jmax = j;
		}
	    }
	}
    
    Mlt(T(1)/sum_coef, vdirector);
        
    // equations of two planes 
    TimesProd(vdirector, this->Normales(imax), vec_u);
    TimesProd(vdirector, this->Normales(jmax), vec_v);

    // pointC is intersection of three planes
    Matrix3_3 mat, inv_mat; R3 rhs;
    mat(0,0) = vec_u(0); mat(0,1) = vec_u(1); mat(0,2) = vec_u(2);
    mat(1,0) = vec_v(0); mat(1,1) = vec_v(1); mat(1,2) = vec_v(2);
    mat(2,0) = vdirector(0); mat(2,1) = vdirector(1); mat(2,2) = vdirector(2);
    rhs(0) = DotProd(vec_u, this->Points(imax));
    rhs(1) = DotProd(vec_v, this->Points(jmax));
    rhs(2) = T(0);
    T det = Det(mat);
    if (abs(det) < 1e-15)
      {
	cout<<"Singular matrix for detection of parameters of the cylinder "<<endl;
	cout<<"Four points are perhaps coplanar "<<endl;
	cout<<"Points used to determine cylinder "<<this->Points<<endl;
	cout<<"Linear system to be inversed "<<mat<<endl;
	cout<<"Right hand side "<<rhs<<endl;
	throw WrongMesh("CylinderParam::FindInitGuess()", "Singular matrix");
      }
    
    GetInverse(mat, inv_mat);
    Mlt(inv_mat, rhs, ptC);
    
    param(0) = ptC(0); param(1) = ptC(1); param(2) = ptC(2);
    // we put a non-null radius, axis and vortho for each point
    param(3) = vdirector(0); param(4) = vdirector(1); param(5) = vdirector(2);
    // estimation of radius
    vec_u = this->Points(imax) - ptC;
    TimesProd(vec_u, vdirector, vec_v); 
    T radius = Norm2(vec_v);
    param(6) = radius;
    // exit(0);
  }


  template<class T>
  void CylinderParam<T>::EvaluateF(const Vector<T>& x, Vector<T>& feval)
  {
    R3 ptC, u, pt_current, CM;
    // point of the cylinder
    ptC.Init(x(0), x(1), x(2));

    // axis of the cylinder
    u.Init(x(3), x(4), x(5));
    T norme_u2 = DotProd(u, u), alpha;

    T radius = x(6);
    
    // loop on points of the cylinder
    int nb = 0;
    for (int i = 0; i < this->nb_points; i++)
      {
	pt_current = this->Points(i);	
	CM = pt_current - ptC;
        alpha = DotProd(CM, u); 
        feval(nb++) = DotProd(CM, CM) - alpha*alpha - radius*radius;
      }

    // unitary axis
    feval(nb++) = norme_u2 - 1.0;
  }
  
  
  template<class T>
  void CylinderParam<T>::EvaluateJacobian(const Vector<T>& x, Vector<T>& feval,
                                          Matrix<T, General, ColMajor>& fjac)
  {
    EvaluateF(x, feval);

    Vector<T> feval_m(this->m_), feval_p(this->m_);
    Vector<T> xtmp(x);
    T h = 1e-6;
    for (int i = 0; i < this->n_; i++)
      {
        xtmp(i) += h;
        EvaluateF(xtmp, feval_p); 
        
        xtmp(i) -= 2.0*h;
        EvaluateF(xtmp, feval_m); 
        
        for (int j = 0; j < this->n; j++)
          fjac(j, i) = (feval_p(j) - feval_m(j)) / (2.0*h);
        
        xtmp(i) = x(i);
      }
  }
  
  
  /**************
   * ConicParam *
   **************/
  
  
  //! default constructor
  template<class T>
  ConicParam<T>::ConicParam()
  {
    this->n = 7; 
    this->nb_points = 12;

    this->m_ = 1 + this->nb_points;
    this->n_ = this->n;
  }
  
  
  //! evaluation of the global "distance" of points to the conic cylinder
  /*!
    \param[in] xsol parameters of the cylinder (point on the axis, axis, angle)
    \param[in] param set of points
    \param[out] feval "distance" of points to the conic cylinder
    If all the points belong to the conic cylinder, the distance will be equal to 0
   */
  template<class T>
  void ConicParam<T>::EvaluateFunction(const Vector<T>& xsol, T& feval)
  {
    R3 ptC, u, CM, pt_current;
    ptC.Init(xsol(0), xsol(1), xsol(2));
    u.Init(xsol(3), xsol(4), xsol(5));
    T tan_theta = xsol(6), alpha, norme_u2, sq_norme_u2, evalg;
    norme_u2 = DotProd(u,u);
    sq_norme_u2 = norme_u2*norme_u2;;
    //theta = atan(tan_theta);
    T tan2_theta = tan_theta*tan_theta;
    feval = 0;
    // loop on points of the conic cylinder
    for (int i = 0; i < this->nb_points; i++)
      {
	pt_current = this->Points(i);
	
	CM = pt_current - ptC;
	alpha = DotProd(CM,u);
	evalg = DotProd(CM,CM)*sq_norme_u2 - alpha*alpha*(norme_u2+tan2_theta);
	feval += evalg*evalg;
      }
    
    feval += (1.0-norme_u2)*(1.0-norme_u2);
  }
  
  //! gradient of the global "distance" of points to the conic cylinder
  /*!
    \param[in] xsol parameters of the conic cylinder (point on the axis, axis, angle)
    \param[in] param set of points
    \param[out] feval "distance" of points to the conic cylinder
    \param[out] gradient  gradient of "distance" of points to the conic cylinder
    If all the points belong to the conic cylinder, the distance will be equal to 0
    The derivation is done with respect to the parameters of the conic cylinder
    (point on the axis, axis, angle)
   */
  template<class T>
  void ConicParam<T>::EvaluateFunctionGradient(const Vector<T>& xsol,
                                               T& feval, Vector<T>& gradient)
  {
    R3 ptC, u, CM, pt_current;
    ptC.Init(xsol(0), xsol(1), xsol(2));
    u.Init(xsol(3), xsol(4), xsol(5));
    T tan_theta = xsol(6), alpha, norme_u2, sq_norme_u2, evalg;
    norme_u2 = DotProd(u,u);
    sq_norme_u2 = norme_u2*norme_u2;;
    //theta = atan(tan_theta);
    T tan2_theta = tan_theta*tan_theta;
    feval = 0; gradient.Fill(0);
    // loop on points of the conic cylinder
    for (int i = 0; i < this->nb_points; i++)
      {
	pt_current = this->Points(i);
	CM = pt_current - ptC;
	
	alpha = DotProd(CM,u);
	evalg = DotProd(CM,CM)*sq_norme_u2 - alpha*alpha*(norme_u2+tan2_theta);
	feval += evalg*evalg;
	
	for (int k = 0; k < 3; k++)
	  {
	    // derivative along C
	    gradient(k) += 2*evalg*(-2.0*CM(k)*sq_norme_u2 + 2.0*alpha*(norme_u2+tan2_theta)*u(k));
	    // derivative along u
	    gradient(k+3) += 2*evalg*(4.0*DotProd(CM,CM)*norme_u2*u(k) - 2.0*alpha*alpha*u(k)
                                      - 2.0*alpha*CM(k)*(norme_u2+tan2_theta));
	  }
	
	gradient(6) -= 4.0*evalg*tan_theta*alpha*alpha;
      }
    
    feval += (1.0-norme_u2)*(1.0-norme_u2);
    for (int k = 0; k < 3; k++)
      gradient(k+3) -= 4.0*(1.0-norme_u2)*u(k);
    
    
    // numeric gradient
    /* Vector<T> grad_bis(7), xp; T h = 1e-5, fp, fm;
    for (int i = 0; i < 7; i++)
      {
	xp = xsol; xp(i) += h;
	EvaluateFunction(xp, param, fp);
	xp = xsol; xp(i) -= h;
	EvaluateFunction(xp, param, fm);
	grad_bis(i) = (fp-fm)/(2.0*h);
      }
    
      DISP(xsol); DISP(feval); DISP(gradient); DISP(grad_bis);*/
  }
  
  //! we find an initial guess for parameters of the conic cylinder
  //! (point on the axis, axis, angle) 
  /*!
    \param[out] param found guess 
   */
  template<class T>
  void ConicParam<T>::FindInitGuess(Vector<T>& param)
  {
    param.Reallocate(this->n);
    param.Fill(T(0));
    Matrix3_3 mat, inv_mat; R3 rhs;
    T det, det_max (0), sum_det(0);
    R3 ptC, vdirector, vec_u, vec_v, vec_w;
    int nb_samples = 0; int imax = -1, jmax = -1, kmax = -1;
    // estimation of center C
    ptC.Fill(0);
    for (int i = 0; i < this->nb_points; i++)
      for (int j = (i+1); j < this->nb_points; j++)
	for (int k = (j+1); k < this->nb_points; k++)
	  {
	    mat(0,0) = this->Normales(i)(0); mat(0,1) = this->Normales(i)(1); mat(0,2) = this->Normales(i)(2);
	    mat(1,0) = this->Normales(j)(0); mat(1,1) = this->Normales(j)(1); mat(1,2) = this->Normales(j)(2);
	    mat(2,0) = this->Normales(k)(0); mat(2,1) = this->Normales(k)(1); mat(2,2) = this->Normales(k)(2);
	    det = Det(mat);
	    if (abs(det) > 1e-4)
	      {
		rhs(0) = DotProd(this->Normales(i), this->Points(i));
		rhs(1) = DotProd(this->Normales(j), this->Points(j));
		rhs(2) = DotProd(this->Normales(k), this->Points(k));
		GetInverse(mat, inv_mat);
		Mlt(inv_mat, rhs, vec_u);
		Add(det*det, vec_u, ptC);
		sum_det += det*det;
		nb_samples++;
		if (abs(det) > det_max)
		  {
		    det_max = abs(det);
		    imax = i; jmax = j; kmax = k;
		  }
	      }
	  }
    if (nb_samples == 0)
      {
	cout<<"Unable to detect center of the cone"<<endl;
	cout<<"Four points are perhaps coplanar "<<endl;
	cout<<"Points used to determine cone "<<this->Points<<endl;
	return;
      }
    
    Mlt(T(1)/T(sum_det), ptC);
    SubConicParam<T> sub;
    vec_u = this->Points(imax) - ptC;
    vec_v = this->Points(jmax) - ptC;
    vec_w = this->Points(kmax) - ptC;
    Mlt(T(1)/Norm2(vec_u), vec_u);
    Mlt(T(1)/Norm2(vec_v), vec_v);
    Mlt(T(1)/Norm2(vec_w), vec_w);
    if (DotProd(vec_u, vec_v) < 0)
      Mlt(T(-1), vec_v);
    if (DotProd(vec_u, vec_w) < 0)
      Mlt(T(-1), vec_w);
    sub.u1 = vec_u; sub.u2 = vec_v; sub.u3 = vec_w;
    
    // direction of axis is estimated
    Vector<T> vdir_param(3);
    T theta = sub.Solve(vdir_param);
    vdirector(0) = vdir_param(0);
    vdirector(1) = vdir_param(1);
    vdirector(2) = vdir_param(2);
    
    param(0) = ptC(0); param(1) = ptC(1); param(2) = ptC(2);
    param(3) = vdirector(0); param(4) = vdirector(1); param(5) = vdirector(2);
    param(6) = tan(theta);
  }


  template<class T>
  void ConicParam<T>::EvaluateF(const Vector<T>& xsol, Vector<T>& feval)
  {
    R3 ptC, u, CM, pt_current;
    ptC.Init(xsol(0), xsol(1), xsol(2));
    u.Init(xsol(3), xsol(4), xsol(5));
    T tan_theta = xsol(6), alpha;
    T norme_u2 = DotProd(u, u);
    T tan2_theta = tan_theta*tan_theta;
    // loop on points of the conic cylinder
    for (int i = 0; i < this->nb_points; i++)
      {
	pt_current = this->Points(i);	
	CM = pt_current - ptC;
	alpha = DotProd(CM, u);
	feval(i) = DotProd(CM, CM) - alpha*alpha*(1.0 + tan2_theta);
      }
    
    feval(this->nb_points) = norme_u2 - 1.0;
  }
  
  
  template<class T>
  void ConicParam<T>::EvaluateJacobian(const Vector<T>& x, Vector<T>& feval,
                                       Matrix<T, General, ColMajor>& fjac)
  {
    EvaluateF(x, feval);

    Vector<T> feval_m(this->m_), feval_p(this->m_);
    Vector<T> xtmp(x);
    T h = 1e-6;
    for (int i = 0; i < this->n_; i++)
      {
        xtmp(i) += h;
        EvaluateF(xtmp, feval_p); 
        
        xtmp(i) -= 2.0*h;
        EvaluateF(xtmp, feval_m); 
        
        for (int j = 0; j < this->n; j++)
          fjac(j, i) = (feval_p(j) - feval_m(j)) / (2.0*h);
        
        xtmp(i) = x(i);
      }
  }
  
    
  //! evaluation of function
  template<class T>
  void SubConicParam<T>::EvaluateFunction(const Vector<T>& param, Vector<T>& feval)
  {
    R3 v,v_u1, v_u2, v_u3;
    T x1, x2, x3;
    v(0) = param(0); v(1) = param(1); v(2) = param(2);
    TimesProd(v, u1, v_u1); x1 = Norm2(v_u1)*Norm2(v_u1);
    TimesProd(v, u2, v_u2); x2 = Norm2(v_u2)*Norm2(v_u2);
    TimesProd(v, u3, v_u3); x3 = Norm2(v_u3)*Norm2(v_u3);
    
    // two equations
    feval(0) = x1-x2;
    feval(1) = x1-x3;
    // equation || v || = 1
    feval(2) = v(0)*v(0) + v(1)*v(1) + v(2)*v(2) - T(1);
  }
  
  
  //! not used
  template<class T>
  void SubConicParam<T>::FindInitGuess(Vector<T>& param) const
  {
    param.Reallocate(3);
    for (int i = 0; i < 3; i++)
      param(i) = (u1(i)+u2(i)+u3(i))/3;
  }
  
  
  //! not used
  template<class T>
  void SubConicParam<T>::EvaluateJacobian(const Vector<T>& param, Matrix<T>& fjac)
  {
    R3 v,v_u1, v_u2, v_u3, u1_v_u1, u2_v_u2, u3_v_u3;
    v(0) = param(0); v(1) = param(1); v(2) = param(2);
    TimesProd(v, u1, v_u1);
    TimesProd(v, u2, v_u2);
    TimesProd(v, u3, v_u3);
    TimesProd(u1, v_u1, u1_v_u1);
    TimesProd(u2, v_u2, u2_v_u2);
    TimesProd(u3, v_u3, u3_v_u3);
    
    // two equations
    for (int i = 0; i < 3; i++)
      {
	fjac(0,i) = T(2)*(u1_v_u1(i) - u2_v_u2(i));
	fjac(1,i) = T(2)*(u1_v_u1(i) - u3_v_u3(i));
	fjac(2,i) = 2*v(i);
      }
  }
  
  
  //! not used
  template<class T>
  T SubConicParam<T>::Solve(Vector<T>& param)
  {
    Vector<T> RControl(10); RControl.Fill(T(0));
    IVect Control(10);
    int N = 3;
    Vector<T> fvec(N), scale(N); param.Reallocate(N);
    Matrix<T> fjac(N,N);
    
    // we set control parameters
    RControl(0) = 1e-15;
    RControl(1) = T(1);
    
    Control(0) = 10;
    Control(3) = 1;
    Control(4) = 0;
    
    FindInitGuess(param);
    
    fvec.Fill(T(0)); scale.Fill(T(0)); fjac.Fill(T(0));
    Montjoie::SolveMinpack(*this, param, fvec, fjac, scale, Control, RControl);
    R3 v,w;
    v(0) = param(0); v(1) = param(1); v(2) = param(2);
    TimesProd(v, u1, w); 
    T theta = asin(Norm2(w));
    
    return theta;
  }

  
}

#define MONTJOIE_FILE_MESH3D_BOUNDARIES_CXX
#endif
