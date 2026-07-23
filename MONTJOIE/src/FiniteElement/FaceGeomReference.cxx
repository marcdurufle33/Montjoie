#ifndef MONTJOIE_FILE_FACE_GEOM_REFERENCE_CXX

namespace Montjoie
{
  //! returns size of memory used by the object
  size_t ElementGeomReference<Dimension2>::GetMemorySize() const
  {
    size_t taille = ElementGeomReference_Base<Dimension2>::GetMemorySize();
    taille += sizeof(*this) - sizeof(ElementGeomReference_Base<Dimension2>);
    taille += sizeof(Real_wp)*(weights2d.GetM()+2*points2d.GetM()+2*points_dof2d.GetM()
                               +2*this->points_nodal_nd.GetM());
    taille += sizeof(int)*(EdgesNode.GetDataSize()+NumNodes2D.GetDataSize()
                           +CoordinateNodes.GetDataSize());
    taille += lob_geom.GetMemorySize() + edge_geom.GetMemorySize();
    return taille;
  }
  

  //! returns quadrature points associated with 2-D integrals only
  VectR2 ElementGeomReference<Dimension2>
  ::PointsQuadInsideND(int nb_points_quadrature_inside) const
  {
    VectR2 pts_quad(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      pts_quad(i) = this->points2d(i);
    
    return pts_quad;
  }
  
  
  /*******************************************************
   * computation of geometric quantities on the boundary *
   *******************************************************/
  
      
  //! computes the curvature on the local edge of an element of a mesh
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel projection of "reference" points on the element
    \param res jacobian matrices on "reference" points, and curvatures
    \param[in] mesh considered mesh
    \param[in] num_elem number of the element
    \param[in] num_edge global number of the edge
    \param[in] num_loc_edge local number of the edge
  */
  void ElementGeomReference<Dimension2>::
  GetCurvatureOnSurface(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			int num_elem, int num_edge, int num_loc_edge) const
  {
    if (mesh.Element(num_elem).IsCurved())
      {
        if (axisym_geometry)
          GetCurvatureOnAxisymmetricSurface(s, PTReel, res, mesh,
                                            num_elem, num_edge, num_loc_edge);
        else
          GetCurvatureOnCurvedSurface(s, PTReel, res, mesh, num_elem, num_edge, num_loc_edge);
      }
    else
      {
        if (axisym_geometry)
          GetCurvatureOnAxisymmetricLinearSurface(s, PTReel, res);
        else
          GetCurvatureOnLinearSurface(s, PTReel, res, mesh, num_elem, num_edge, num_loc_edge);
      }
  }
  
  
  //! computes the curvature on the local edge of an element of a mesh
  //! in case of LINEAR elements
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel projection of "reference" points on the element
    \param res jacobian matrices on "reference" points, and curvatures
    \param[in] mesh considered mesh
    \param[in] num_elem number of the element
    \param[in] num_edge global number of the edge
    \param[in] num_loc_edge local number of the edge
  */
  void ElementGeomReference<Dimension2>::
  GetCurvatureOnLinearSurface(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			      SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			      int num_elem, int num_edge, int num_loc_edge) const
  {
    int hyb = mesh.Element(num_elem).GetHybridType();
    int vertex2 = mesh.Element(num_elem).
      numVertex(MeshNumbering<Dimension2>::FirstExtremityEdge(hyb, num_loc_edge));
    
    int vertex3 = mesh.Element(num_elem).
      numVertex(MeshNumbering<Dimension2>::SecondExtremityEdge(hyb, num_loc_edge));
    
    int ref = mesh.Boundary(num_edge).GetReference(), iedge;
    IVect ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill();
    int vertex1 = mesh.FindFollowingVertex(vertex3, vertex2, ref, ref_cond, iedge);
    int vertex4 = mesh.FindFollowingVertex(vertex2, vertex3, ref, ref_cond, iedge);
    Real_wp h_courbure,x,dsj; R2 normale_fj,dn_ds;
    // DISP(vertex1); DISP(vertex2); DISP(vertex3); DISP(vertex4); DISP(ref);
    if ((vertex1 < 0)&&(vertex4 < 0))
      {
	// we don't find neighbourging vertices, we assume that we have an isolated edge
	// normale is constant and the courbure is null
	h_courbure = Real_wp(0);
	for (int i = 0; i < this->points1d.GetM(); i++)
	  res.SetK1QuadratureBoundary(i, h_courbure);
        
	return;
      }
    
    if (vertex1 < 0)
      {
	// vertex2 is the first vertex of a curve
	// then we get phi_2 = 1/2(1-x)(2-x) phi_3 = x(2-x) phi_4 = -1/2 (1-x) x
	// dphi_2 = 1/2 (-3+2x) dphi_3 = 2-2x dphi_4 = 1/2(2x -1)
	// d2phi_2 = 1 d2phi_3 = -2 d2phi_4 = 1
	
	R2 pt2,pt3,pt4,n2,n3,n4,pt; pt2 = mesh.Vertex(vertex2);
	pt3 = mesh.Vertex(vertex3); pt4 = mesh.Vertex(vertex4);
	n2.Init(pt2(1),-pt2(0)); n3.Init(pt3(1),-pt3(0)); n4.Init(pt4(1),-pt4(0));
	
	for (int i = 0; i < this->points1d.GetM(); i++)
	  {
	    pt.Fill(0); normale_fj.Fill(0); h_courbure = 0; dn_ds.Fill(0);
	    x = this->points1d(i);
	    Add(0.5*(1.0-x)*(2.0-x),pt2,pt);
	    Add(x*(2.0-x),pt3,pt);
	    Add(-0.5*x*(1.0-x),pt4,pt);
	    
	    Add(-1.5+x,n2,normale_fj);
	    Add(2.0-2.0*x,n3,normale_fj);
	    Add(x-0.5,n4,normale_fj);
	    
            dn_ds = n2 - Real_wp(2)*n3 + n4;
            
	    dsj = Norm2(normale_fj);
	    Mlt(1.0/dsj,normale_fj); Mlt(1.0/dsj,dn_ds);
	    h_courbure = dn_ds(0)*normale_fj(1) - dn_ds(1)*normale_fj(0);
	    h_courbure *= -1.0/dsj;
	    //PTReel.SetPointQuadratureBoundary(i,pt);
	    //res.SetNormaleQuadratureBoundary(i,normale_fj);
	    //res.SetDsQuadratureBoundary(i,dsj);
	    // h_courbure = 0.0; 
	    res.SetK1QuadratureBoundary(i,h_courbure);
	  }
	return;
      }
    if (vertex4<0)
      {
	// vertex4 is the first vertex of a curve
	// then we get phi_1 = -1/2x(1-x) phi_2 = (1+x)(1-x) phi_3 = 1/2 (1+x) x
	// dphi_1 = 1/2 (2x-1) dphi_2 = -2x dphi_3 = 1/2(2x+1)
	R2 pt1,pt2,pt3,n1,n2,n3,pt; pt1 = mesh.Vertex(vertex1);
	pt2 = mesh.Vertex(vertex2); pt3 = mesh.Vertex(vertex3);
	n1.Init(pt1(1),-pt1(0)); n2.Init(pt2(1),-pt2(0)); n3.Init(pt3(1),-pt3(0));
	
	for (int i = 0; i < this->points1d.GetM(); i++)
	  {
	    pt.Fill(0); normale_fj.Fill(0); h_courbure = 0.0; dn_ds.Fill(0);
	    x = this->points1d(i);
	    Add(-0.5*x*(1.0-x),pt1,pt);
	    Add((1.0+x)*(1.0-x),pt2,pt);
	    Add(0.5*x*(1.0+x),pt3,pt);
	    
	    Add(x-1.0,n1,normale_fj);
	    Add(-2.0*x,n2,normale_fj);
	    Add(0.5+x,n3,normale_fj);
	    
            dn_ds = n1 - Real_wp(2)*n2 + n3;
	    
	    dsj = Norm2(normale_fj);
	    Mlt(1.0/dsj, normale_fj); Mlt(1.0/dsj, dn_ds);

	    h_courbure = dn_ds(0)*normale_fj(1) - dn_ds(1)*normale_fj(0);
	    h_courbure *= -1.0/dsj;
	    //  PTReel.SetPointQuadratureBoundary(i,pt);
	    //res.SetNormaleQuadratureBoundary(i,normale_fj);
	    //res.SetDsQuadratureBoundary(i,dsj);
	    // h_courbure = 0.0; 
	    res.SetK1QuadratureBoundary(i,h_courbure);
	  }
	return;
      }
    
    // we get phi_1 = -x(1.0-x)(2.0-x)/6 phi_2 = (1+x)*(1.0-x)*(2.0-x)/2
    // phi_3 =(1+x) x(2-x)/2 phi_4 = -(1+x)(1.0-x)x/6
    // dphi_1 = 1/6 (-3x^2+6x-2) dphi_2 = 1/2 (-1-4x+3x^2)
    // dphi_3 = 1/2 (2+2x-3x^2) dphi_4 = 1/6 (3x^2-1)
    // d2phi_1 = 1-x dphi_2 = 3*x-2 dphi_3 = 1-3*x dphi_4 = x
    R2 n1, n2, n3, n4, pt1, pt2, pt3, pt4, pt;
    pt1 = mesh.Vertex(vertex1);
    pt2 = mesh.Vertex(vertex2);
    pt3 = mesh.Vertex(vertex3);
    pt4 = mesh.Vertex(vertex4);
    n1.Init(pt1(1),-pt1(0));
    n2.Init(pt2(1),-pt2(0));
    n3.Init(pt3(1),-pt3(0));
    n4.Init(pt4(1),-pt4(0));
    
    // DISP(n1); DISP(n2); DISP(n3); DISP(n4);
    for (int i = 0; i < this->points1d.GetM(); i++)
      {
	pt.Fill(0); normale_fj.Fill(0); dn_ds.Fill(0);
	h_courbure = 0.0;
	x = this->points1d(i);
	Add(-x*(1.0-x)*(2.0-x)/6,pt1,pt);
	Add((1.0+x)*(1.0-x)*(2.0-x)/2,pt2,pt);
	Add(x*(1.0+x)*(2.0-x)/2,pt3,pt);
	Add(-x*(1.0-x)*(1.0+x)/6,pt4,pt);
	
	Add((-2.0+6.0*x-3.0*x*x)/6,n1,normale_fj);
	Add((-1.0-4.0*x+3.0*x*x)/2,n2,normale_fj);
	Add((2.0+2.0*x-3.0*x*x)/2,n3,normale_fj);
	Add((-1.0+3.0*x*x)/6,n4,normale_fj);
	
	Add(1.0-x,n1,dn_ds); Add(3.0*x-2.0,n2,dn_ds);
	Add(1.0-3.0*x,n3,dn_ds);
	Add(x,n4,dn_ds);
	
	dsj = Norm2(normale_fj);
	Mlt(1.0/dsj,normale_fj); Mlt(1.0/dsj,dn_ds);
	h_courbure = dn_ds(0)*normale_fj(1) - dn_ds(1)*normale_fj(0);
	h_courbure *= -1.0/dsj;
	//PTReel.SetPointQuadratureBoundary(i,pt);
	//res.SetNormaleQuadratureBoundary(i,normale_fj);
	//res.SetDsQuadratureBoundary(i,dsj);
	// h_courbure = 0.0; 
	res.SetK1QuadratureBoundary(i,h_courbure);
	// DISP(i); DISP(x); DISP(normale_fj); DISP(dn_ds); DISP(dsj); DISP(h_courbure);
      }
  }
  
  
  //! computes the curvature on the local edge of an element of a mesh 
  //! in case of CURVED elements
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel projection of "reference" points on the element
    \param res jacobian matrices on "reference" points, and curvatures
    \param[in] mesh considered mesh
    \param[in] num_elem number of the element
    \param[in] num_edge global number of the edge
    \param[in] num_loc_edge local number of the edge
  */
  void ElementGeomReference<Dimension2>::
  GetCurvatureOnCurvedSurface(const VectR2& s, const SetPoints<Dimension2>& PTReel,
			      SetMatrices<Dimension2>& res, const Mesh<Dimension2>& mesh,
			      int num_elem, int num_edge, int num_loc_edge) const
  {
    Globatto<Real_wp> lob; lob.AffectPoints(this->points1d);
    lob.ComputeGradPhi(1e3*epsilon_machine);
    R2 dn_ds, normale_fj; Real_wp h_courbure = 0;
    for (int i = 0; i < this->points1d.GetM(); i++)
      {
	dn_ds.Fill(0);
	for (int j = 0; j < this->points1d.GetM(); j++)
	  Add(lob.GradPhi(j, i), res.GetNormaleQuadratureBoundary(j), dn_ds);
	
	normale_fj = res.GetNormaleQuadratureBoundary(i);
	h_courbure = dn_ds(0)*normale_fj(1) - dn_ds(1)*normale_fj(0);
	h_courbure *= -Real_wp(1)/res.GetDsQuadratureBoundary(i);
	res.SetK1QuadratureBoundary(i, h_courbure);
      }
  }
  
  
  //! computes curvature on nodal points
  void ElementGeomReference<Dimension2>::
  GetCurvatureNodal(const Globatto<Real_wp>& lob,
		    const VectR2& normale, const VectReal_wp& ds,
		    VectReal_wp& curvature)
  {
    curvature.Reallocate(lob.GetOrder()+1);
    R2 dn_ds, normale_fj; Real_wp h_courbure(0);
    for (int i = 0; i <= lob.GetOrder(); i++)
      {
	dn_ds.Fill(0);
	for (int j = 0; j <= lob.GetOrder(); j++)
	  Add(lob.GradPhi(j, i), normale(j), dn_ds);
	
	normale_fj = normale(i);
	h_courbure = dn_ds(0)*normale_fj(1) - dn_ds(1)*normale_fj(0);
	h_courbure *= -Real_wp(1)/ds(i);
	curvature(i) = h_courbure;
      }
  }
  
  
  //! computes the two curvatures on the local edge
  //! of an element of a mesh in axisymmetric case
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel projection of "reference" points on the element
    \param res jacobian matrices on "reference" points, and curvatures
    \param[in] mesh considered mesh
    \param[in] num_elem number of the element
    \param[in] num_edge global number of the edge
    \param[in] num_loc_edge local number of the edge
    case of a curved element
  */
  void ElementGeomReference<Dimension2>::
  GetCurvatureOnAxisymmetricSurface(const VectR2& s,
				    const SetPoints<Dimension2>& PTReel,
				    SetMatrices<Dimension2>& res,
				    const Mesh<Dimension2>& mesh,
				    int num_elem, int num_edge, int num_loc_edge) const
  {
    Globatto<Real_wp> gauss_;
    int r = this->points1d.GetM()-1;
    gauss_.ConstructQuadrature(r, gauss_.QUADRATURE_LOBATTO);
    gauss_.AffectPoints(this->points1d);
    gauss_.ComputeGradPhi(1e3*epsilon_machine);
    // DISP(lob.valGrad);
    R2 dn_ds, normale_fj; Real_wp ri, k1_curve = 0, k2_curve = 0;
    for (int i = 0; i < this->points1d.GetM(); i++)
      {
	dn_ds.Fill(0);
	for (int j = 0; j < this->points1d.GetM(); j++)
	  Add(gauss_.GradPhi(j,i), res.GetNormaleQuadratureBoundary(j), dn_ds);
	
	normale_fj = res.GetNormaleQuadratureBoundary(i);
	k1_curve = dn_ds(0)*normale_fj(1) - dn_ds(1)*normale_fj(0);
	k1_curve *= Real_wp(1)/res.GetDsQuadratureBoundary(i); // DISP(k1_curve);
	ri = PTReel.GetPointQuadratureBoundary(i)(0); // DISP(ri);
	if (abs(ri) <= 100.0*epsilon_machine*mesh.GetXmax())
	  {
	    if (i==0)
	      k2_curve = -dn_ds(0)/res.GetDsQuadratureBoundary(i);
	    else
	      k2_curve = dn_ds(0)/res.GetDsQuadratureBoundary(i);
	  }
	else
	  k2_curve = -normale_fj(0)/ri;
	
	res.SetK1QuadratureBoundary(i, -k1_curve);
	res.SetK2QuadratureBoundary(i, -k2_curve);  // DISP(k2_curve);
      }
  }
  
  
  //! computes the two curvatures on the local edge
  //! of an element of a mesh in axisymmetric case
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel projection of "reference" points on the element
    \param res jacobian matrices on "reference" points, and curvatures
    case of a straight element (one curvature is null)
  */
  void ElementGeomReference<Dimension2>::
  GetCurvatureOnAxisymmetricLinearSurface(const VectR2& s,
					  const SetPoints<Dimension2>& PTReel,
					  SetMatrices<Dimension2>& res) const
  {
    Real_wp k2_curve, k1_curve, ri;
    for (int i = 0; i < this->points1d.GetM(); i++)
      {
	ri = PTReel.GetPointQuadratureBoundary(i)(0);
	k2_curve = -res.GetNormaleQuadratureBoundary(i)(0)/ri;
	k1_curve = 0.0;
	res.SetK1QuadratureBoundary(i, -k1_curve);
	res.SetK2QuadratureBoundary(i, -k2_curve);
      }
  }
  
  
  /* methods used to compute F_i^{-1} */
  
  
  //! computes a bounding box of the polygon defined by the vertices s
  /*!
    \param[in] s list of vertices of the polygon
    \param[in] coef minimal distance between the bounding box and the polygon
    \param[out] box bounding polygon, that contains s
    \param[out] enveloppe bounding rectangular box
   */
  void ElementGeomReference<Dimension2>::GetBoundingBox(const VectR2& s,
                                                        const Real_wp& coef, VectR2& box,
                                                        TinyVector<R2, 2>& enveloppe)
  {
    int nb_vertices = s.GetM();
    // enveloppe(0) is the left bottom corner of the smallest rectangle
    // bounding the polygon defined by s
    
    // enveloppe(1) is the right top corner of the smallest rectangle
    // bounding the polygon defined by s
    enveloppe(0) = s(0);
    enveloppe(1) = s(0);
    if (coef == 0)
      {
	for (int i = 0; i < nb_vertices; i++)
	  {
	    box(i) = s(i);
	    enveloppe(0).Init(min(box(i)(0), enveloppe(0)(0)),
			      min(box(i)(1), enveloppe(0)(1)));
	    
	    enveloppe(1).Init(max(box(i)(0), enveloppe(1)(0)),
			      max(box(i)(1), enveloppe(1)(1)));
	  }
      }
    else
      {
	int num_next;
	R2 pointB0,pointB1;
	R2 pointC0(s(nb_vertices-1)(0) + (s(0)(1)-s(nb_vertices-1)(1))*coef,
		   s(nb_vertices-1)(1) - (s(0)(0)-s(nb_vertices-1)(0))*coef);
	
	R2 pointC1(s(0)(0) + (s(0)(1)-s(nb_vertices-1)(1))*coef,
		   s(0)(1) - (s(0)(0)-s(nb_vertices-1)(0))*coef);
	
        Real_wp delta = max(abs(s(0)(0)), abs(s(0)(1)));
        for (int i = 1; i < nb_vertices; i++)
          delta = max(delta, max(abs(s(i)(0)), abs(s(i)(1))));

	for (int i = 0; i < nb_vertices; i++)
	  {
	    num_next = i+1;
	    if (num_next == nb_vertices)
	      num_next = 0;
	    
	    pointB0.Init(s(i)(0) + (s(num_next)(1)-s(i)(1))*coef,
			 s(i)(1) - (s(num_next)(0)-s(i)(0))*coef);
	    pointB1.Init(s(num_next)(0) + (s(num_next)(1)-s(i)(1))*coef,
			 s(num_next)(1) - (s(num_next)(0)-s(i)(0))*coef);
	    
	    IntersectionDroites(pointB0, pointB1, pointC0, pointC1,
                                box(i), 10.0*epsilon_machine*delta);
	    
            enveloppe(0)(0) = min(box(i)(0), enveloppe(0)(0));
	    enveloppe(0)(1) = min(box(i)(1), enveloppe(0)(1));
	    enveloppe(1)(0) = max(box(i)(0), enveloppe(1)(0));
	    enveloppe(1)(1) = max(box(i)(1), enveloppe(1)(1));
	    pointC0 = pointB0;
	    pointC1 = pointB1;
	  }    
      }
  }
  
  
  //! Intersection between the edge pointA-pointB and a polygon defined by box 
  //! (quadrangle or triangle)
  /*!
    \param[in] nb_vertices number of vertices of the polygon box
    \param[in] box list of vertices of the polygon
    \param[in] pointA first extremity of the edge
    \param[in] pointB second extremity of the edge
    \param[out] res intersection points of the polygon with the edge
    \return True if an intersection has been found
   */
  bool ElementGeomReference<Dimension2>::ComputeIntersection(int nb_vertices, const VectR2& box,
							     const R2& pointA, const R2& pointB,
							     VectR2& res)
  {
    R2 temp;
    VectR2 point_intersec(4);
    int ind=0,result,num_next;
    for (int i = 0; i < nb_vertices; i++)
      {
	if (i == (nb_vertices-1))
	  num_next = 0;
	else
	  num_next = i+1;
        
	result = IntersectionEdges(box(i), box(num_next), pointA,
				   pointB, temp, epsilon_machine);

	if (result == 2)
	  {
	    R2 point1(box(i));
	    R2 point2(box(num_next)); 
	    Sort(point1, point2);
	    if ((pointB >= box(i))&&(pointA <= box(num_next)))
	      {
		res(0) = max(pointA, point1);
		res(1) = min(pointB, point2);
		return true;
	      }
	    return false;
	  }   
	else if (result == 1)
	  {
	    point_intersec(ind) = temp;
	    ind++;
	  }	  
      }
    
    if (ind == 1)
      {
	res(0) = point_intersec(0);
	res(1) = point_intersec(0);
	return true;
      }
    else if (ind > 1)
      {
	Assemble(ind,point_intersec);
	res(0) = point_intersec(0);
	res(1) = point_intersec(1);
	return true;
      }
    
    return false;
  }
  
    
  /*****************
   * Other methods *
   *****************/
  
  
  //! returns true if the element is affine
  /*!
    \param[in] s vertices of the element
    triangles are always affine, parallelograms are affine
   */
  bool ElementGeomReference<Dimension2>::IsAffineTransformation(const VectR2& s)
  {
    R2 vec_u, vec_v;
    switch (s.GetM())
      {
      case 3 :
	// straight triangle -> Fi affine
	return true;
      case 4 :
	// quadrilateral
	vec_u = s(1) - s(0);
	vec_v = s(2) - s(3);
	if (vec_u == vec_v)
	  return true;
	
	return false;
      }
    
    return false;
  }
  
  
  //! displays informations about the class FaceReference<type>
  ostream& operator <<(ostream& out, const ElementGeomReference<Dimension2>& e)
  {
    out<<"Order of geometry "<<e.order_geom<<endl;
    return out;
  }
  
}
  
#define MONTJOIE_FILE_FACE_GEOM_REFERENCE_CXX
#endif
