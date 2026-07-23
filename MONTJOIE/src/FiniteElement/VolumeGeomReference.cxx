#ifndef MONTJOIE_FILE_VOLUME_GEOM_REFERENCE_CXX

namespace Montjoie
{
  //! returns size of memory used by the object
  size_t ElementGeomReference<Dimension3>::GetMemorySize() const
  {
    size_t taille = ElementGeomReference_Base<Dimension3>::GetMemorySize();
    taille += sizeof(Real_wp)*(weights3d.GetM()+3*points3d.GetM()+3*points_dof3d.GetM()
                               +3*this->points_nodal_nd.GetM()+2*points_nodal2d_tri.GetM()
                               +2*points_nodal2d_quad.GetM()+2*points2d_tri.GetM()
                               +2*points2d_quad.GetM()+weights2d_tri.GetM()+weights2d_quad.GetM());
    
    taille += sizeof(Real_wp)*(2*points_dof2d_tri.GetM()+2*points_dof2d_quad.GetM());
    
    taille += sizeof(int)*(FacesNodal.GetDataSize()+NumNodes3D.GetDataSize()
                           +CoordinateNodes.GetDataSize()+NumNodes2D_tri.GetDataSize()
                           +CoordinateNodes2D_tri.GetDataSize()+NumNodes2D_quad.GetDataSize()
                           +CoordinateNodes2D_quad.GetDataSize());
    
    taille += sizeof(Real_wp)*CoefCurveTri.GetDataSize();
    taille += sizeof(int)*FacesNodal_Rotation_Tri.GetDataSize();
    taille += sizeof(int)*FacesNodal_Rotation_Quad.GetDataSize();
    taille += sizeof(Real_wp)*RhTriGeomTrans.GetDataSize();
    taille += sizeof(int)*(NumQuad2D.GetDataSize());
    taille += tangente_loc_x.GetMemorySize() + tangente_loc_y.GetMemorySize();
    taille += sizeof(*this) - sizeof(ElementGeomReference_Base<Dimension3>);
    return taille;
  }


  //! returns quadrature points associated with 2-D integrals only
  VectR3 ElementGeomReference<Dimension3>
  ::PointsQuadInsideND(int nb_points_quadrature_inside) const
  {
    VectR3 pts_quad(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      pts_quad(i) = this->points3d(i);
    
    return pts_quad;
  }

  
  /*******************************************************
   * Computation of geometric quantities on the boundary *
   *******************************************************/
  
  
  //! computes enclosing polyedron
  void ElementGeomReference<Dimension3>
  ::GetBoundingBox(const VectR3& s,
		   const Real_wp coef, VectR3& box, TinyVector<R3, 2>& enveloppe)
  {
    int nb_vertices = s.GetM();
    // enveloppe(0) is the left bottom corner
    // of the smallest straight parallepipede bounding the polyhedra defined by s
    // enveloppe(1) is the right top corner of the 
    // smallest straight parallepipede bounding the polyhedra defined by s
    enveloppe(0) = s(0);
    enveloppe(1) = s(0); 
    for (int i = 0; i < nb_vertices; i++)
      {
        box(i) = s(i);
        enveloppe(0).Init(std::min(s(i)(0), enveloppe(0)(0)),
                          std::min(s(i)(1), enveloppe(0)(1)),
                          std::min(s(i)(2), enveloppe(0)(2)));
        
        enveloppe(1).Init(std::max(s(i)(0), enveloppe(1)(0)),
                          std::max(s(i)(1), enveloppe(1)(1)),
                          std::max(s(i)(2), enveloppe(1)(2)));
      }

    if (coef != Real_wp(0))
      {
        Real_wp dx = enveloppe(1)(0) - enveloppe(0)(0);
        Real_wp dy = enveloppe(1)(1) - enveloppe(0)(1);
        Real_wp dz = enveloppe(1)(2) - enveloppe(0)(2);
        enveloppe(0)(0) -= coef*dx; enveloppe(1)(0) += coef*dx;
        enveloppe(0)(1) -= coef*dy; enveloppe(1)(1) += coef*dy;
        enveloppe(0)(2) -= coef*dz; enveloppe(1)(2) += coef*dz;
      }
  }
  

  void ElementGeomReference<Dimension3>
  ::InitTriangularCurvature(const TriangleGeomReference& tri)
  {
    // matrix used for computation of curvature 
    RhTriGeomTrans.Reallocate(2*this->Points2D_tri().GetM(),
			      this->PointsNodal2D_tri().GetM());
    RhTriGeomTrans.Fill(0);
    VectR2 grad_phi;
    for (int i = 0; i < this->Points2D_tri().GetM(); i++)
      {
        tri.ComputeGradientPhiNodalRef(this->Points2D_tri()(i), grad_phi);
        for (int j = 0; j < this->PointsNodal2D_tri().GetM(); j++)
          {
            RhTriGeomTrans(2*i, j) = grad_phi(j)(0);
            RhTriGeomTrans(2*i+1, j) = grad_phi(j)(1);
          }
      }
  }

  
  void ElementGeomReference<Dimension3>::SetMesh(Mesh<Dimension3>& mesh)
  {
    ElementGeomReference_Base<Dimension3>::SetMesh(mesh);
    
    int r = mesh.GetGeometryOrder();
    // on triangles
    if ((this->order_geom > 2) && (this->points_nodal2d_tri.GetM() > 0))
      {
	TriangleGeomReference tri;
	tri.ConstructFiniteElement(r);
	// constructing nodal points from nodal points provided by the mesh
	VectR2 nodal_points(tri.GetNbPointsNodalElt());
	nodal_points(0) = R2(0,0);
	nodal_points(1) = R2(1,0);
	nodal_points(2) = R2(0,1);
	// points on edges
	for (int i = 1; i < r; i++)
	  {
	    Real_wp x = mesh.GetNodalPointInsideEdge(i-1);
	    nodal_points(i+2) = R2(x, 0);
	    nodal_points(i+r+1) = R2(1-x, x);
	    nodal_points(i+2*r) = R2(0, 1-x);
	  }
        
	// points inside triangle
	for (int i = 3*r; i < tri.GetNbPointsNodalElt(); i++)
	  nodal_points(i) = mesh.GetNodalPointInsideTriangle(i-3*r);
	
	tri.SetNodalPoints(r, nodal_points);
	
	CoefCurveTri.Reallocate((this->order_geom-1)*(this->order_geom-2)/2, (r+1)*(r+2)/2);
	VectReal_wp phi;
	for (int i = 0; i < CoefCurveTri.GetM(); i++)
	  {
	    tri.ComputeValuesPhiNodalRef(points_nodal2d_tri(3*this->order_geom+i), phi);
	    for (int j = 0; j < tri.GetNbPointsNodalElt(); j++)
	      CoefCurveTri(i, j) = phi(j);
	  }
      }
  }

  
  //! transformation Fi for nodal points 
  /*!
    \param[in] s vertices of the element
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
  */
  void ElementGeomReference<Dimension3>::
  FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad) const
  {
    // here we treat only points on edges and faces (interior depends on the element)
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    // we consider that the order of approximation of curves in the mesh
    // can be different from order of finite element
    int r = mesh.GetGeometryOrder();
    
    // vertices
    int nv = mesh.Element(nquad).GetNbVertices();
    int type = mesh.Element(nquad).GetHybridType();
    for (int k = 0; k < nv; k++)
      res.SetPointNodal(k, s(k));
    
    R3 res_n; Real_wp x, y;
    
    if ((this->CoefCurve1D.GetM() != this->order_geom-1)||(this->CoefCurve1D.GetN() != r+1))
      {
	cout<<"Did you call SetMesh after ConstructFiniteElement ?"<<endl;
	abort();
      }
    
    // edges
    int nb_edges = mesh.Element(nquad).GetNbEdges();
    for (int j = 0; j < nb_edges; j++)
      {
	// we get the edge number and the reference
	int num_edge = mesh.Element(nquad).numEdge(j);
	// DISP(num_edge_ref); DISP(PointsNodal1D);
	if (num_edge >= mesh.GetNbEdgesRef())
	  {
	    // case of a straight edge
	    // two extremities of the edge
	    R3 ptA = s(MeshNumbering<Dimension3>::FirstExtremityEdge(type, j));
	    R3 ptB = s(MeshNumbering<Dimension3>::SecondExtremityEdge(type, j));	    
	    for (int k = 1; k < this->order_geom; k++)
	      {
		res_n.Zero();
		// linear interpolation
		Add(Real_wp(1)-this->points_nodal1d(k), ptA, res_n); 
		Add(this->points_nodal1d(k), ptB, res_n);
		res.SetPointNodal(j*(this->order_geom-1)+k+nv-1, res_n); 
                //DISP(j); DISP(k); DISP(j*(order-1)+k+nv-1); DISP(res_n);
	      }
	  }
	else
	  {
	    // first step : we get the points on the curved edge
	    VectR3 PtsEdge(r+1);
	    PtsEdge(0) = s(MeshNumbering<Dimension3>::FirstExtremityEdge(type, j));
	    PtsEdge(r) = s(MeshNumbering<Dimension3>::SecondExtremityEdge(type, j));
	    if (mesh.Element(nquad).GetOrientationEdge(j))
	      for (int k = 1; k < r; k++)
		PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, k-1);
	    else
	      for (int k = 1; k < r; k++)
		PtsEdge(k) = mesh.GetPointInsideEdge(num_edge, r-k-1);
	    
	    // second step : interpolation to nodal points of the finite element
	    for (int k = 1; k < this->order_geom; k++)
	      {
		res_n.Zero();
		for (int p = 0; p <= r; p++)
		  Add(this->CoefCurve1D(k-1, p), PtsEdge(p), res_n);
		
		res.SetPointNodal(j*(this->order_geom-1)+k+nv-1, res_n);
	      }
	  }
      }
    
    int nb = nv + nb_edges*(this->order_geom-1);
    int Ntri = (this->order_geom-1)*(this->order_geom-2)/2;
    // faces
    for (int j = 0; j < mesh.Element(nquad).GetNbFaces(); j++)
      {
	// we get the face number and the reference
	int num_face = mesh.Element(nquad).numFace(j);
	if (mesh.Boundary(num_face).GetNbVertices() == 3)
	  {
	    // triangle
	    if (num_face >= mesh.GetNbBoundaryRef())
	      {
		// case of a straight face
		R3 ptA = s(MeshNumbering<Dimension3>::FirstExtremityFace(type, j));
		R3 ptB = s(MeshNumbering<Dimension3>::SecondExtremityFace(type, j));	    
		R3 ptC = s(MeshNumbering<Dimension3>::ThirdExtremityFace(type, j));
		for (int k = 0; k < Ntri; k++)
		  {
		    res_n.Zero();
		    // linear interpolation
		    x = points_nodal2d_tri(3*this->order_geom+k)(0);
		    y = points_nodal2d_tri(3*this->order_geom+k)(1);
		    Add(Real_wp(1)-x-y, ptA, res_n); 
		    Add(x, ptB, res_n); Add(y, ptC, res_n);
		    // DISP(k); DISP(nb+k); DISP(res_n);
		    res.SetPointNodal(nb + k, res_n); 
		  }
	      }
	    else
	      {
		// case of a curved face
		// first step : we get the points on the curved face
		Vector<R3> PtsFace((r+1)*(r+2)/2);
		// vertices
		PtsFace(0) = s(MeshNumbering<Dimension3>::FirstExtremityFace(type, j));
		PtsFace(1) = s(MeshNumbering<Dimension3>::SecondExtremityFace(type, j));
		PtsFace(2) = s(MeshNumbering<Dimension3>::ThirdExtremityFace(type, j));
		// points on edges
		for (int num_loc = 0; num_loc < 3; num_loc++)
		  {
		    int num_edge = mesh.BoundaryRef(num_face).numEdge(num_loc);
		    if (num_edge >= mesh.GetNbEdgesRef())
		      {
			for (int k = 1; k < r; k++)
			  {
			    res_n.Zero();
			    Add(Real_wp(1)-mesh.GetNodalPointInsideEdge(k-1),
                                PtsFace(num_loc), res_n);
			    Add(mesh.GetNodalPointInsideEdge(k-1), PtsFace((num_loc+1)%3), res_n);
			    PtsFace((r-1)*num_loc + 2 + k) = res_n;
			  }
		      }
		    else
		      {
			if (mesh.EdgeRef(num_edge).numVertex(0)
                            == mesh.BoundaryRef(num_face).numVertex(num_loc))
			  for (int k = 1; k < r; k++)
			    PtsFace((r-1)*num_loc + 2 + k)
                              = mesh.GetPointInsideEdge(num_edge, k-1);
			else
			  for (int k = 1; k < r; k++)
			    PtsFace((r-1)*num_loc + 2 + k)
                              = mesh.GetPointInsideEdge(num_edge, r-k-1);
		      }
		  }
		
		// points inside face
		for (int k = 0; k < (r-1)*(r-2)/2; k++)
		  PtsFace(3*r+k) = mesh.GetPointInsideFace(num_face, k);
                
		// second step : interpolation
		int rot = mesh.Element(nquad).GetOrientationFace(j);
		for (int k = 0; k < (this->order_geom-1)*(this->order_geom-2)/2; k++)
		  {
		    res_n.Zero();
		    for (int p = 0; p < (r+1)*(r+2)/2; p++)
		      Add(CoefCurveTri(k, p), PtsFace(p), res_n);
                    
		    res.SetPointNodal(nb + FacesNodal_Rotation_Tri(rot, k), res_n);
		  }
	      }
	    nb += (this->order_geom-1)*(this->order_geom-2)/2;
	  }
	else
	  {
	    Real_wp one(1);
	    // quadrangle
	    if (num_face >= mesh.GetNbBoundaryRef())
	      {
		// case of a straight face
		R3 ptA = s(MeshNumbering<Dimension3>::FirstExtremityFace(type, j));
		R3 ptB = s(MeshNumbering<Dimension3>::SecondExtremityFace(type, j));	    
		R3 ptC = s(MeshNumbering<Dimension3>::ThirdExtremityFace(type, j));
		R3 ptD = s(MeshNumbering<Dimension3>::FourthExtremityFace(type, j)); 
		for (int k1 = 1; k1 < this->order_geom; k1++)
		  for (int k2 = 1; k2 < this->order_geom; k2++)
		    {
		      res_n.Zero();
		      // bilinear interpolation
		      x = this->points_nodal1d(k1);
		      y = this->points_nodal1d(k2);
		      Add((one-x)*(one-y), ptA, res_n); 
		      Add(x*(one-y), ptB, res_n);
		      Add(x*y, ptC, res_n);
		      Add((one-x)*y, ptD, res_n);
		      res.SetPointNodal(nb + (k1-1)*(this->order_geom-1) + k2-1, res_n); 
		    }
	      }
	    else
	      {
		// case of a curved face
		// first step : we get the points on the curved face
		Vector<R3> PtsFace((r+1)*(r+1)), PtsInterm((r+1)*(this->order_geom+1));
		Matrix<int> num_nodes, coord_nodes;
		// vertices
		for (int num_loc = 0; num_loc < 4; num_loc++)
		  PtsFace(num_loc) = mesh.Vertex(mesh.Boundary(num_face).numVertex(num_loc));
		
		MeshNumbering<Dimension2>::
                  ConstructQuadrilateralNumbering(r, num_nodes, coord_nodes);
		
		// points on edges
		for (int num_loc = 0; num_loc < 4; num_loc++)
		  {		      
		    int num_edge = mesh.BoundaryRef(num_face).numEdge(num_loc);
		    if (num_edge >= mesh.GetNbEdgesRef())
		      {
			for (int k = 1; k < r; k++)
			  {
			    res_n.Zero();
			    Add(one-mesh.GetNodalPointInsideEdge(k-1), PtsFace(num_loc), res_n);
			    Add(mesh.GetNodalPointInsideEdge(k-1), PtsFace((num_loc+1)%4), res_n);
			    PtsFace((r-1)*num_loc + 3 + k) = res_n;
			  }
		      }
		    else
		      {
			if (mesh.EdgeRef(num_edge).numVertex(0)
                            == mesh.BoundaryRef(num_face).numVertex(num_loc))
			  for (int k = 1; k < r; k++)
			    PtsFace((r-1)*num_loc + 3 + k)
                              = mesh.GetPointInsideEdge(num_edge, k-1);
			else
			  for (int k = 1; k < r; k++)
			    PtsFace((r-1)*num_loc + 3 + k)
                              = mesh.GetPointInsideEdge(num_edge, r-k-1);
		      }
		  }
		// points inside face
		for (int k = 0; k < (r-1)*(r-1); k++)
		  PtsFace(4*r+k) = mesh.GetPointInsideFace(num_face, k);
		
		// second step : interpolation
		int rot = mesh.Element(nquad).GetOrientationFace(j);
		for (int k1 = 1; k1 < this->order_geom; k1++)
		  for (int p2 = 0; p2 <= r; p2++)
		    {
		      PtsInterm(k1*(r+1) + p2).Zero();
		      for (int p1 = 0; p1 <= r; p1++)
			Add(this->CoefCurve1D(k1-1, p1),
                            PtsFace(num_nodes(p1, p2)), PtsInterm(k1*(r+1)+p2));
		    }
		
		for (int k1 = 1; k1 < this->order_geom; k1++)
		  for (int k2 = 1; k2 < this->order_geom; k2++)
		    {
		      res_n.Zero();
		      for (int p2 = 0; p2 <= r; p2++)
			Add(this->CoefCurve1D(k2-1, p2), PtsInterm(k1*(r+1)+p2), res_n);
		      
		      int k = (k1-1)*(this->order_geom-1) + k2-1;
		      //DISP(k); DISP(nb+FacesNodal_Rotation_Quad(rot,k)); DISP(res_n);
		      res.SetPointNodal(nb+FacesNodal_Rotation_Quad(rot, k), res_n);
		    }
	      }
	    nb += (this->order_geom-1)*(this->order_geom-1);
	  }
      }
    // DISP(res);
  }
  

  //! computation of curvature tensor on quadrature points of the face
  void ElementGeomReference<Dimension3>::
  GetCurvatureOnSurface(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			SetMatrices<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int num_elem,
			int num_face, int num_loc) const
  {
    if (mesh.Element(num_elem).IsCurved())
      {
        // the curvature tensor is computed only for curved elements
        // C = grad n where n is the normale
        int Ntri = this->PointsNodal2D_tri().GetM();
        int Nquad = this->GetNbQuadBoundary(num_loc);
        
        // nabla_nx, nabla_ny, nabla_nz : gradient for the reference face (\hat{nabla})
        VectReal_wp nabla_nx(2*Nquad), nabla_ny(2*Nquad), nabla_nz(2*Nquad);
        Matrix3_3 mat_dfj, dfjm1; R3 normale_fj;

        // curvature tensor C
        Matrix<Real_wp> tensorC(3, 3);

        // we compute the gradient of nx, ny, nz with respect
        // to \hat{x} \hat{y} where \hat{x}, \hat{y}
        // are coordinates on the unit square or unit triangle (depending on the face)
        // the result is stored in nabla_nx, nabla_ny, nabla_nz
        if (this->IsLocalFaceQuadrilateral(num_loc))
          {
            R3 dn_dx, dn_dy;
	    
	    // 1-D basis functions based on quadrature points
	    // this object is used to compute the gradient of n on
	    // the quadrilateral face from values of n on quadrature points
	    Globatto<Real_wp> lob; lob.AffectPoints(this->points1d);
	    lob.ComputeGradPhi(1e3*epsilon_machine);
                    
            // quadrilateral base
            // we are exploiting tensorization to compute
            // efficiently the gradient with respect to the unit square
            for (int i = 0; i < this->points1d.GetM(); i++)
              for (int j = 0; j < this->points1d.GetM(); j++)
                {
                  int npoint = this->NumQuad2D(i, j);
                  dn_dx.Zero(); dn_dy.Zero();
                  for (int k = 0; k < points1d.GetM(); k++)
                    {
                      // derivative with respect to x
                      int n = this->NumQuad2D(k, j);
                      Add(lob.GradPhi(k, i), res.GetNormaleQuadratureBoundary(n), dn_dx);

                      // derivative with respect to y
                      n = this->NumQuad2D(i, k);
                      Add(lob.GradPhi(k, j), res.GetNormaleQuadratureBoundary(n), dn_dy);
                    }
                  
                  // we store the values on nabla_nx, nabla_ny, nabla_nz
                  nabla_nx(2*npoint) = dn_dx(0);
                  nabla_nx(2*npoint+1) = dn_dy(0);

                  nabla_ny(2*npoint) = dn_dx(1);
                  nabla_ny(2*npoint+1) = dn_dy(1);

                  nabla_nz(2*npoint) = dn_dx(2);
                  nabla_nz(2*npoint+1) = dn_dy(2);
                }
          }
        else
          {
            // triangular face
            // we retrieve values of nx, ny, nz on nodal points of the face
            VectReal_wp nx(Ntri), ny(Ntri), nz(Ntri);
            for (int i = 0; i < Ntri; i++)
              {
                mat_dfj = res.GetPointNodal(this->FacesNodal(i, num_loc));
                Real_wp deter = Det(mat_dfj);
                GetInverse(mat_dfj, dfjm1);
                MltTrans(dfjm1, this->NormaleLoc(num_loc), normale_fj);
                Mlt(deter, normale_fj);
                Real_wp dsj = Norm2(normale_fj);
                Mlt(Real_wp(Real_wp(1)/dsj), normale_fj);
                nx(i) = normale_fj(0);
                ny(i) = normale_fj(1);
                nz(i) = normale_fj(2);
              }
            
            // then the gradient is computed
            Mlt(RhTriGeomTrans, nx, nabla_nx);
            Mlt(RhTriGeomTrans, ny, nabla_ny);
            Mlt(RhTriGeomTrans, nz, nabla_nz);
          }
        

        R3 d_nx, d_ny, d_nz;
        R3 grad_nx, grad_ny, grad_nz;
        for (int i = 0; i < Nquad; i++)
          {
            // we compute the gradient of nx ny and nz on the 3-D face
            // from gradient on the unit triangle or unit square
            // the result is stored in d_nx, d_ny and d_nz            
	    this->GetGradient3D_FromGradient2D(num_loc, i, nabla_nx, nabla_ny, nabla_nz,
					       d_nx, d_ny, d_nz);

            // we compute the gradient of n
            // with respect to coordinates x, y, z of the real element
            // by using jacobian matrix DF (associated with the
            // transformation F between reference element and real element K)
            GetInverse(res.GetPointQuadratureBoundary(i), dfjm1);
            MltTrans(dfjm1, d_nx, grad_nx);
            MltTrans(dfjm1, d_ny, grad_ny);
            MltTrans(dfjm1, d_nz, grad_nz);
            
            // last step, the gradient of n is shaped into a matrix
            // tensorC, and we compute eigenvalues of this matrix
            // we keep only the two non-null eigenvalues (called K1 and K2)
            // h = 0.5(k1+k2) is the mean curvature, g = k1 k2 the gaussian curvature
            tensorC(0, 0) = grad_nx(0);
            tensorC(0, 1) = grad_nx(1);
            tensorC(0, 2) = grad_nx(2);
            tensorC(1, 0) = grad_ny(0);
            tensorC(1, 1) = grad_ny(1);
            tensorC(1, 2) = grad_ny(2);
            tensorC(2, 0) = grad_nz(0);
            tensorC(2, 1) = grad_nz(1);
            tensorC(2, 2) = grad_nz(2);              
            
            VectReal_wp lambda(3), lambda_imag(3);
            GetEigenvalues(tensorC, lambda, lambda_imag);
            
            if (abs(lambda(2)) <= min(abs(lambda(0)), abs(lambda(1))))
              {
                res.SetK1QuadratureBoundary(i, lambda(0));
                res.SetK2QuadratureBoundary(i, lambda(1));
              }
            
            if (abs(lambda(1)) < min(abs(lambda(0)), abs(lambda(2))))
              {
                res.SetK1QuadratureBoundary(i, lambda(0));
                res.SetK2QuadratureBoundary(i, lambda(2));
              }
            
            if (abs(lambda(0)) < min(abs(lambda(1)), abs(lambda(2))))
              {
                res.SetK1QuadratureBoundary(i, lambda(1));
                res.SetK2QuadratureBoundary(i, lambda(2));
              }
          }
      }
    else
      {
        // straight element : null curvature
        for (int k = 0; k < PTReel.GetNbPointsQuadratureBoundary(); k++)
          {
            res.SetK1QuadratureBoundary(k, 0.0);
            res.SetK2QuadratureBoundary(k, 0.0);
          }
      }
  }


  /*****************
   * Other methods *
   *****************/
  
  
  //! returns true if the transformation Fi is affine
  /*!
    The transformation Fi is affine for tetrahedra,
    and for other elements if quadrilateral faces are parallelogramms
   */
  bool ElementGeomReference<Dimension3>::IsAffineTransformation(const VectR3& s)
  {
    R3 vec_u, vec_v, vec_w;
    switch (s.GetM())
      {
      case 4 :
	// straight tetrahedron -> Fi affine
	return true;
      case 5 :
	// pyramid, we check the quadrilateralbase
	vec_u = s(1) - s(0);
	vec_v = s(2) - s(3);
	if (vec_u == vec_v)
	  return true;
	
	return false;
      case 6 :
	// wedge, we check the three quadrilateral faces
	vec_u = s(3) - s(0);
	vec_v = s(4) - s(1);
	vec_w = s(5) - s(2);
	if ((vec_v == vec_u)&&(vec_w == vec_v))
	  return true;
	
	return false;
      case 8 :
	// hexahedron, first we check face z = 0
	vec_u = s(1) - s(0);
	vec_v = s(2) - s(3);
	if (vec_u != vec_v)
	  return false;
	
	// then lateral faces
	vec_u = s(4) - s(0);
	vec_v = s(5) - s(1);
	if (vec_u != vec_v)
	  return false;
	
	vec_v = s(6) - s(2);
	if (vec_u != vec_v)
	  return false;
	
	vec_v = s(7) - s(3);
	if (vec_u != vec_v)
	  return false;
	
	return true;
      }
    
    return false;
  }
  
  
  //! displays information about class VolumeReference<type>
  ostream& operator <<(ostream& out, const ElementGeomReference<Dimension3>& e)
  {
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_VOLUME_GEOM_REFERENCE_CXX
#endif
