#ifndef MONTJOIE_FILE_HEXAHEDRON_GEOM_REFERENCE_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronGeomReference::HexahedronGeomReference() : ElementGeomReference<Dimension3>()
  {
    this->normale.Reallocate(6);
    this->normale(0).Init(-1,0,0);
    this->normale(1).Init(0,-1,0);
    this->normale(2).Init(0,0,-1);
    this->normale(3).Init(0,0,1);
    this->normale(4).Init(0,1,0);
    this->normale(5).Init(1,0,0);

    this->tangente_loc_x.Reallocate(6);
    this->tangente_loc_y.Reallocate(6);
    this->tangente_loc_x(0).Init(0,1,0);
    this->tangente_loc_y(0).Init(0,0,1);
    this->tangente_loc_x(1).Init(1,0,0);
    this->tangente_loc_y(1).Init(0,0,1);
    this->tangente_loc_x(2).Init(1,0,0);
    this->tangente_loc_y(2).Init(0,1,0);
    this->tangente_loc_x(3).Init(1,0,0);
    this->tangente_loc_y(3).Init(0,1,0);
    this->tangente_loc_x(4).Init(1,0,0);
    this->tangente_loc_y(4).Init(0,0,1);
    this->tangente_loc_x(5).Init(0,1,0);
    this->tangente_loc_y(5).Init(0,0,1);
    
    this->hybrid_type_elt = 3;
    this->nb_vertices_elt = 8;
    this->nb_edges_elt = 12;
    this->nb_boundaries_elt = 6;
    
    this->is_local_face_quad.Reallocate(6);
    this->is_local_face_quad.Fill(true);
  }
  

  //! returns size of memory used by the object
  size_t HexahedronGeomReference::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension3>::GetMemorySize();
    taille += function_geom_quad.GetMemorySize();
    taille += sizeof(Real_wp)*coefFi.GetDataSize();
    taille += sizeof(int)*(10+3*coefFi_curve.GetM() + 3*numberFi_curve.GetM());
    for (int i = 0; i < coefFi_curve.GetM(); i++)
      taille += sizeof(Real_wp)*coefFi_curve(i).GetM();

    for (int i = 0; i < numberFi_curve.GetM(); i++)
      taille += sizeof(int)*numberFi_curve(i).GetM();
    
    return taille;
  }
  
  
  //! construction of the finite element
  void HexahedronGeomReference::ConstructFiniteElement(int rgeom)
  {
    this->order_geom = rgeom;

    ConstructNodalShapeFunctions(rgeom);
    
    ComputeCurvedTransformation();
    ComputeCoefficientTransformation();
  }
    
  
  //! construction of shape functions
  /*!
    By calling this function only (and not ConstructFiniteElement)
    shape functions are constructed and you can call ComputeValuesPhiNodalRef
   */
  void HexahedronGeomReference::ConstructNodalShapeFunctions(int r)
  {
    // clearing previous points
    this->points_nodal1d.Clear(); this->points_nodal2d_quad.Clear();
    this->points_nodal2d_tri.Clear(); this->points_nodal_nd.Clear();
    
    function_geom_quad.ConstructFiniteElement(r);
    
    // NumNodes3D(i,j,k) contains the number of the nodal points
    // whose coordinates are (\xi_i,\xi_j,\xi_k)
    MeshNumbering<Dimension3>::
      ConstructHexahedralNumbering(r, this->NumNodes3D, this->CoordinateNodes);
    
    // NumNodes2D_quad(i,j) contains the number of quadrature points on face
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(r, this->NumNodes2D_quad,
				      this->CoordinateNodes2D_quad);

    // 1-D nodal points are Gauss-Lobatto points
    lob_geom.ConstructQuadrature(r, lob_geom.QUADRATURE_LOBATTO);
    lob_geom.ComputeGradPhi(1e3*epsilon_machine);
    this->points_nodal1d = lob_geom.Points();
    
    // 2-D nodal points
    this->points_nodal2d_quad.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        this->points_nodal2d_quad(this->NumNodes2D_quad(i, j))
          .Init(this->points_nodal1d(i), this->points_nodal1d(j));
    
    // 3-D nodal points
    this->points_nodal_nd.Reallocate((r+1)*(r+1)*(r+1));    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	for (int k = 0; k <= r; k++)
          this->points_nodal_nd(this->NumNodes3D(i, j, k))
            .Init(this->points_nodal1d(i), this->points_nodal1d(j), this->points_nodal1d(k));
    
    // filling FacesNodal
    this->FacesNodal.Reallocate((r+1)*(r+1), 6);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	this->FacesNodal(this->NumNodes2D_quad(i, j), 0) = this->NumNodes3D(0, i, j);
    
    // quadrature points on second face
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	this->FacesNodal(this->NumNodes2D_quad(i, j), 1) = this->NumNodes3D(i, 0, j);
    
    // quadrature points on third face
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	this->FacesNodal(this->NumNodes2D_quad(i, j), 2) = this->NumNodes3D(i, j, 0);
    
    // quadrature points on fourth face
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	this->FacesNodal(this->NumNodes2D_quad(i, j), 3) = this->NumNodes3D(i, j, r);
 
    // quadrature points on fifth face
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	this->FacesNodal(this->NumNodes2D_quad(i, j), 4) = this->NumNodes3D(i, r, j);
    
    // quadrature points on sixth face
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	this->FacesNodal(this->NumNodes2D_quad(i, j), 5) = this->NumNodes3D(r, i, j);
    
    if (r > 1)
      {
        Matrix<int> NumNodes2D_InsideQuad(r-1, r-1);
        for (int i = 1;  i < r; i++)
          for (int j = 1;  j < r; j++)
            NumNodes2D_InsideQuad(i-1, j-1) = (i-1)*(r-1) + j-1;
        
        MeshNumbering<Dimension3>::
          GetRotationQuadrilateralFace(NumNodes2D_InsideQuad, this->FacesNodal_Rotation_Quad);
      }
  }
  

  //! computation of coefficients in order to obtain faster evaluations of Fi and DFi  
  void HexahedronGeomReference::ComputeCurvedTransformation()
  {
    // for curvilinear case
    numberFi_curve.Reallocate(this->points_nodal_nd.GetM());
    coefFi_curve.Reallocate(this->points_nodal_nd.GetM());
    
    int r = this->order_geom;
    Real_wp x,y,z, one(1);
    // use curved transformation to get nodal points inside the hexahedron
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)      
	for (int k = 1; k < r; k++)      
	  {
	    int node = this->NumNodes3D(i, j, k);
	    x = this->points_nodal_nd(node)(0);
	    y = this->points_nodal_nd(node)(1);
	    z = this->points_nodal_nd(node)(1);
	  
	    numberFi_curve(node).Reallocate(26);
	    coefFi_curve(node).Reallocate(26);
	    // vertices
	    numberFi_curve(node)(0) = this->NumNodes3D(0, 0, 0);
	    coefFi_curve(node)(0) = (one-x)*(one-y)*(one-z);
	    numberFi_curve(node)(1) = this->NumNodes3D(r, 0, 0);
	    coefFi_curve(node)(1) = x*(one-y)*(one-z);
	    numberFi_curve(node)(2) = this->NumNodes3D(r, r, 0);
	    coefFi_curve(node)(2) = x*y*(one-z);
	    numberFi_curve(node)(3) = this->NumNodes3D(0, r, 0);
	    coefFi_curve(node)(3) = (one-x)*y*(one-z);
	    numberFi_curve(node)(4) = this->NumNodes3D(0, 0, r);
	    coefFi_curve(node)(4) = (one-x)*(one-y)*z;
	    numberFi_curve(node)(5) = this->NumNodes3D(r, 0, r);
	    coefFi_curve(node)(5) = x*(one-y)*z;
	    numberFi_curve(node)(6) = this->NumNodes3D(r, r, r);
	    coefFi_curve(node)(6) = x*y*z;
	    numberFi_curve(node)(7) = this->NumNodes3D(0, r, r);
	    coefFi_curve(node)(7) = (one-x)*y*z;
	    
	    // edges
	    numberFi_curve(node)(8) = this->NumNodes3D(i, 0, 0);
	    coefFi_curve(node)(8) = -(one-y)*(one-z);
	    numberFi_curve(node)(9) = this->NumNodes3D(r, j, 0);
	    coefFi_curve(node)(9) = -x*(one-z);
	    numberFi_curve(node)(10) = this->NumNodes3D(i, r, 0);
	    coefFi_curve(node)(10) = -y*(one-z);
	    numberFi_curve(node)(11) = this->NumNodes3D(0, j, 0);
	    coefFi_curve(node)(11) = -(one-x)*(one-z);
	    numberFi_curve(node)(12) = this->NumNodes3D(0, 0, k);
	    coefFi_curve(node)(12) = -(one-x)*(one-y);
	    numberFi_curve(node)(13) = this->NumNodes3D(r, 0, k);
	    coefFi_curve(node)(13) = -x*(one-y);
	    numberFi_curve(node)(14) = this->NumNodes3D(r, r, k);
	    coefFi_curve(node)(14) = -x*y;
	    numberFi_curve(node)(15) = this->NumNodes3D(0, r, k);
	    coefFi_curve(node)(15) = -(one-x)*y;
	    numberFi_curve(node)(16) = this->NumNodes3D(i, 0, r);
	    coefFi_curve(node)(16) = -(one-y)*z;
	    numberFi_curve(node)(17) = this->NumNodes3D(r, j, r);
	    coefFi_curve(node)(17) = -x*z;
	    numberFi_curve(node)(18) = this->NumNodes3D(i, r, r);
	    coefFi_curve(node)(18) = -y*z;
	    numberFi_curve(node)(19) = this->NumNodes3D(0, j, r);
	    coefFi_curve(node)(19) = -(one-x)*z;
	    
	    // faces
	    numberFi_curve(node)(20) = this->NumNodes3D(0, j, k);
	    coefFi_curve(node)(20) = one-x;
	    numberFi_curve(node)(21) = this->NumNodes3D(i, 0, k);
	    coefFi_curve(node)(21) = one-y;
	    numberFi_curve(node)(22) = this->NumNodes3D(i, j, 0);
	    coefFi_curve(node)(22) = one-z;
	    numberFi_curve(node)(23) = this->NumNodes3D(i, j, r);
	    coefFi_curve(node)(23) = z;
	    numberFi_curve(node)(24) = this->NumNodes3D(i, r, k);
	    coefFi_curve(node)(24) = y;
	    numberFi_curve(node)(25) = this->NumNodes3D(r, j, k);
	    coefFi_curve(node)(25) = x;
	  }
  }

  
  //! computation of coefficients in order to obtain faster evaluations of Fi and DFi  
  void HexahedronGeomReference::ComputeCoefficientTransformation()
  {
    // for linear case
    int nb_points_all = this->points_nodal_nd.GetM() + this->points_dof3d.GetM()
      + this->points3d.GetM();
    
    coefFi.Reallocate(nb_points_all,8);
    Real_wp x,y,z, one(1);
    // F_i = (1-x)*(1-y)*(1-z) A_0 + (1-x)*(1-y)*z A_1 + (1-x)*y*(1-z) A_2
    //   + (1-x)*y*z A_3 + x*(1-y)*(1-z) A_4 + x*(1-y)*z A_5 + x*y*(1-z) A_6 + x*y*z A_7
    // where A_i are the eight vertices of the hexahedron
    VectR3 all_points = this->points_nodal_nd;
    all_points.PushBack(this->points3d);
    all_points.PushBack(this->points_dof3d);
    for (int i = 0; i < nb_points_all; i++)
      {
	x = all_points(i)(0);
	y = all_points(i)(1);
	z = all_points(i)(2);
		
	coefFi(i,0) = (one-x)*(one-y)*(1-z);
	coefFi(i,1) =   x  *(one-y)*(1-z);
	coefFi(i,2) =   x  *   y *(1-z);
	coefFi(i,3) = (one-x)*   y *(1-z);
	coefFi(i,4) = (one-x)*(one-y)*   z ;
	coefFi(i,5) =   x  *(one-y)*   z ;
	coefFi(i,6) =   x  *   y *   z ;
	coefFi(i,7) = (one-x)*   y *   z ;
      }
    
    // for curvilinear case
    numberFi_curve.Resize(nb_points_all);
    coefFi_curve.Resize(nb_points_all);
    
    VectReal_wp value;
    for (int i = this->points_nodal_nd.GetM(); i < nb_points_all; i++)
      {
	ComputeValuesPhiNodalRef(all_points(i), value);
	// first test to count number of non-zeros
	int nnz = 0;
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  if (abs(value(j)) > 1e5*epsilon_machine)
	    nnz++;
	
	coefFi_curve(i).Reallocate(nnz);
	numberFi_curve(i).Reallocate(nnz);
	nnz = 0;
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  if (abs(value(j)) > 1e5*epsilon_machine)
	    {
	      coefFi_curve(i)(nnz) = value(j);
	      numberFi_curve(i)(nnz) = j;
	      nnz++;
	    }
      }
  }
  

  void HexahedronGeomReference
  ::GetGradient3D_FromGradient2D(int num_loc, int i, const VectReal_wp& nabla_nx,
				 const VectReal_wp& nabla_ny, const VectReal_wp& nabla_nz,
				 R3& d_nx, R3& d_ny, R3& d_nz) const
  {
    switch (num_loc)
      {
      case 0:
      case 5:
	d_nx(0) = 0; d_nx(1) = nabla_nx(2*i); d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = 0; d_ny(1) = nabla_ny(2*i); d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = 0; d_nz(1) = nabla_nz(2*i); d_nz(2) = nabla_nz(2*i+1);
	break;
      case 1:
      case 4:
	d_nx(0) = nabla_nx(2*i); d_nx(1) = 0; d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = nabla_ny(2*i); d_ny(1) = 0; d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = nabla_nz(2*i); d_nz(1) = 0; d_nz(2) = nabla_nz(2*i+1);
	break;
      case 2:
      case 3:
	d_nx(0) = nabla_nx(2*i); d_nx(1) = nabla_nx(2*i+1); d_nx(2) = 0;
	d_ny(0) = nabla_ny(2*i); d_ny(1) = nabla_ny(2*i+1); d_ny(2) = 0;
	d_nz(0) = nabla_nz(2*i); d_nz(1) = nabla_nz(2*i+1); d_nz(2) = 0;
      }
  }


  //! returns a tangential vector on a face of the hexahedron from a 2-D vector on the square
  R3 HexahedronGeomReference::GetTangentialVector(int num_loc, const R2& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
      case 5:
	return R3(Real_wp(0), vec_u(0), vec_u(1));
      case 1:
      case 4:
	return R3(vec_u(0), Real_wp(0), vec_u(1));
      case 2:
      case 3:
	return R3(vec_u(0), vec_u(1), Real_wp(0));
      }

    return R3();
  }


  //! transpose operation of GetTangentialVector
  R2 HexahedronGeomReference::TransposeTangentialVector(int num_loc, const R3& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
      case 5:
	return R2(vec_u(1), vec_u(2));
      case 1:
      case 4:
	return R2(vec_u(0), vec_u(2));
      case 2:
      case 3:
	return R2(vec_u(0), vec_u(1));
      }

    return R2();
  }


  //! transpose operation of GetTangentialVector
  R2_Complex_wp HexahedronGeomReference::TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
      case 5:
	return R2_Complex_wp(vec_u(1), vec_u(2));
      case 1:
      case 4:
	return R2_Complex_wp(vec_u(0), vec_u(2));
      case 2:
      case 3:
	return R2_Complex_wp(vec_u(0), vec_u(1));
      }

    return R2_Complex_wp();
  }

    
  //! not implemented
  void HexahedronGeomReference
  ::ComputeCoefJacobian(const VectR3& s, VectReal_wp& CoefJacobian) const
  {
    abort();
  }
  
  
  /****************
   * Fj transform *
   ****************/
  

  //! transformation Fi in the case of straight hexahedron  
  void HexahedronGeomReference::FjLinear(const VectR3& s,
					 const R3& point, R3& res) const
  {
    Real_wp one(1);
    Real_wp x = point(0), y = point(1), z = point(2);
    // Fi = (1-x)*{ (1-y)*[ (1-z) S0 + z S4 ] + y * [ (1-z) S3 + z S7 ] } + 
    //          x*{ (1-y)*[ (1-z) S1 + z S5 ] + y * [ (1-z) S2 + z S6 ] }
    for (int j = 0; j < 3; j++)
      res(j) = (one-x) * ( (one-y) * ( (one-z)*s(0)(j) + z*s(4)(j)  ) + 
	               	     y     * ( (one-z)*s(3)(j) + z*s(7)(j)) ) + 
	         x    *  ( (one-y) * ( (one-z)*s(1)(j) + z*s(5)(j)  ) + 
			     y     * ( (one-z)*s(2)(j) + z*s(6)(j)  ) );
  }
  
  
  //! transformation DFi in the case of straight hexahedron
  void HexahedronGeomReference::DFjLinear(const VectR3& s, const R3& point,
					  Matrix3_3& res) const
  {
    Real_wp x = point(0), y = point(1), z = point(2), one(1);
    for (int j = 0; j < 3; j++)
      {
	// dFi/dx = (1-y)(1-z) (S1-S0) + (1-y)z (S5-S4) + y(1-z) (S2-S3) + yz (S6-S7)
	// dFi/dy = (1-x)(1-z) (S3-S0) + (1-x)z (S7-S4) + x(1-z) (S2-S1) + xz (S6-S5)
	// dFi/dz = (1-x)(1-y) (S4-S0) + (1-x)y (S7-S3) + x(1-y) (S5-S1) + xy (S6-S2)
	res(j,0) = (one-y)*(one-z)*(s(1)(j)-s(0)(j)) + (one-y)*z*(s(5)(j)-s(4)(j)) + 
	  y*(one-z)*(s(2)(j)-s(3)(j)) + y*z*(s(6)(j)-s(7)(j));
	
	res(j,1) = (one-x)*(one-z)*(s(3)(j)-s(0)(j)) + (one-x)*z*(s(7)(j)-s(4)(j)) + 
	  x*(one-z)*(s(2)(j)-s(1)(j)) + x*z*(s(6)(j)-s(5)(j));
	
	res(j,2) = (1-x)*(1-y)*(s(4)(j)-s(0)(j)) + (one-x)*y*(s(7)(j)-s(3)(j))+
	  x*(one-y)*(s(5)(j)-s(1)(j)) + x*y*(s(6)(j)-s(2)(j));
      }
  }
  
  
  //! returns the minimal length of the element
  Real_wp HexahedronGeomReference
  ::GetMinimalSize(const VectR3& s) const
  {
    Real_wp h = Distance(s(0), s(1));
    h = min(h, Distance(s(0), s(3)));
    h = min(h, Distance(s(0), s(4)));
    h = min(h, Distance(s(1), s(2)));
    h = min(h, Distance(s(1), s(5)));
    h = min(h, Distance(s(2), s(3)));
    h = min(h, Distance(s(2), s(6)));
    h = min(h, Distance(s(3), s(7)));
    h = min(h, Distance(s(4), s(5)));
    h = min(h, Distance(s(4), s(7)));
    h = min(h, Distance(s(5), s(6)));
    h = min(h, Distance(s(6), s(7)));
    return h;
  }
  

  //! returns true if the point is outside the unit cube
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
   */
  bool HexahedronGeomReference::
  OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const
  {
    Real_wp Rmax(0);
    for (int i = 0; i < 8; i++)
      Rmax = max(max(max(Rmax, abs(s(i)(0))), abs(s(i)(1))), abs(s(i)(2)));
    
    Real_wp dx = s(0).Distance(s(1)), dy = s(0).Distance(s(3)), dz = s(0).Distance(s(4));
    if ((Rmax == Real_wp(0)) || (dx == Real_wp(0)) || (dy == Real_wp(0)) || (dz == Real_wp(0)))
      return true;

    if (dx*Xn(0) < -epsilon*Rmax)
      return true;

    if (dy*Xn(1) < -epsilon*Rmax)
      return true;

    if (dz*Xn(2) < -epsilon*Rmax)
      return true;

    if (dx*(Xn(0) - Real_wp(1)) > epsilon*Rmax)
      return true;

    if (dy*(Xn(1) - Real_wp(1)) > epsilon*Rmax)
      return true;

    if (dz*(Xn(2) - Real_wp(1)) > epsilon*Rmax)
      return true;
    
    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
   */
  Real_wp HexahedronGeomReference::GetDistanceToBoundary(const R3& pointloc) const
  {
    Real_wp distance = pointloc(0);
    distance = min(distance,Real_wp(Real_wp(1) - pointloc(0)));
    
    distance = min(distance,pointloc(1));
    distance = min(distance,Real_wp(Real_wp(1) - pointloc(1)));
    distance = min(distance,pointloc(2));
    distance = min(distance,Real_wp(Real_wp(1) - pointloc(2)));
    
    return distance;
  }


  Real_wp HexahedronGeomReference::GetDistanceToBoundary(const R3& pointloc, int n) const
  {
    Real_wp distance(0);
    switch (n)
      {
      case 0:
	{
	  distance = abs(pointloc(0));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	  distance = max(distance, -pointloc(2));
	  distance = max(distance, pointloc(2)-Real_wp(1));
	}
	break;
      case 1 :
	{
	  distance = abs(pointloc(1));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	  distance = max(distance, -pointloc(2));
	  distance = max(distance, pointloc(2)-Real_wp(1));
	}
	break;
      case 2 :
	{
	  distance = abs(pointloc(2));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	}
	break;
      case 3 :
	{
	  distance = abs(pointloc(2)-Real_wp(1));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	}
	break;
      case 4 :
	{
	  distance = abs(pointloc(1)-Real_wp(1));
	  distance = max(distance, -pointloc(0));
	  distance = max(distance, pointloc(0)-Real_wp(1));
	  distance = max(distance, -pointloc(2));
	  distance = max(distance, pointloc(2)-Real_wp(1));
	}
	break;
      case 5 :
	{
	  distance = abs(pointloc(0)-Real_wp(1));
	  distance = max(distance, -pointloc(1));
	  distance = max(distance, pointloc(1)-Real_wp(1));
	  distance = max(distance, -pointloc(2));
	  distance = max(distance, pointloc(2)-Real_wp(1));
	}
	break;
      }

    return distance;
  }
  
  
  //! project the local point on the boundary, if outside the unit cube
  int HexahedronGeomReference::ProjectPointOnBoundary(R3& pointloc) const
  {
    pointloc(0) = max( Real_wp(0), pointloc(0) );
    pointloc(0) = min( Real_wp(1), pointloc(0) ); 
    pointloc(1) = max( Real_wp(0), pointloc(1) );
    pointloc(1) = min( Real_wp(1), pointloc(1) ); 
    pointloc(2) = max( Real_wp(0), pointloc(2) );
    pointloc(2) = min( Real_wp(1), pointloc(2) ); 
    
    return -1;
  }
  
  
  //! returns 3-D coordinates of a point on a face
  /*!
    \param[in] num_loc local face number in the quad
    \param[in] point_loc coordinate on the face
    \param[out] res local coordinates of the point in the unit square
  */
  void HexahedronGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const
  {
    switch(num_loc)
      {
      case 0:
	res.Init(0, point_loc(0), point_loc(1));break;
      case 1:
	res.Init(point_loc(0), 0, point_loc(1));break;
      case 2:
	res.Init(point_loc(0), point_loc(1), 0);break;
      case 3:
	res.Init(point_loc(0), point_loc(1), Real_wp(1));break;
      case 4:
	res.Init(point_loc(0), Real_wp(1), point_loc(1));break;
      case 5:
	res.Init(Real_wp(1), point_loc(0), point_loc(1));break;
      }
  }


  //! computes 2-D coordinates of a point on a face
  /*!
    \param[in] num_loc local face number in the quad
    \param[out] point_loc coordinate on the face
    \param[in] pt local coordinates of the point in the unit cube
  */
  void HexahedronGeomReference::
  GetLocalCoordOnBoundary(int num_loc, const R3& pt, R2& point_loc) const
  {
    switch(num_loc)
      {
      case 0:
	point_loc.Init(pt(1), pt(2)); break;
      case 1:
	point_loc.Init(pt(0), pt(2)); break;
      case 2:
	point_loc.Init(pt(0), pt(1)); break;
      case 3:
	point_loc.Init(pt(0), pt(1)); break;
      case 4:
	point_loc.Init(pt(0), pt(2)); break;
      case 5:
	point_loc.Init(pt(1), pt(2)); break;
      }
  }
  

  //! transformation Fi in the case of curved hexahedron
  void HexahedronGeomReference::
  FjCurve(const SetPoints<Dimension3>& PTReel, const R3& pointloc, R3& res) const
  {
    VectReal_wp phi;
    ComputeValuesPhiNodalRef(pointloc, phi);
    
    res.Zero();
    for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
      Add(phi(k), PTReel.GetPointNodal(k), res);
    
  }
  
  
  //! transformation DFi in the case of curved hexahedron
  void HexahedronGeomReference::
  DFjCurve(const SetPoints<Dimension3>& PTReel,
	   const R3& pointloc, Matrix3_3& res) const
  {
    int r = this->order_geom;
    VectR3 value_phi(r+1), value_dphi(r+1);
    // we compute Gauss-Lobatto functions at point x, y and z
    for (int k1 = 0; k1 <= r; k1++)
      {
	for (int j = 0; j < 3; j++)
	  value_phi(k1)(j) = lob_geom.EvaluatePhi(k1, pointloc(j));
      }
    
    // we compute derivative of gauss-lobatto functions at point x, y and z
    for (int k1 = 0; k1 <= r; k1++)
      {
	value_dphi(k1).Fill(0);
	for (int j = 0; j <= r; j++)
	  Add(lob_geom.GradPhi(k1,j), value_phi(j), value_dphi(k1));
      }
    
    TinyVector<R3, 3> dfi;
    for (int k1 = 0; k1 <= r; k1++)
      for (int k2 = 0; k2 <= r; k2++)
	for (int k3 = 0; k3 <= r; k3++)
	  {
	    Add(value_dphi(k1)(0)*value_phi(k2)(1)*value_phi(k3)(2),
		PTReel.GetPointNodal(this->NumNodes3D(k1,k2,k3)), dfi(0));
	    
	    Add(value_phi(k1)(0)*value_dphi(k2)(1)*value_phi(k3)(2),
		PTReel.GetPointNodal(this->NumNodes3D(k1,k2,k3)), dfi(1));
	    
	    Add(value_phi(k1)(0)*value_phi(k2)(1)*value_dphi(k3)(2),
		PTReel.GetPointNodal(this->NumNodes3D(k1,k2,k3)), dfi(2));
	  }
    
    res = Matrix3_3(dfi);
    
  }
  
  
  //! computes res = Fi(point) for all reference points
  /*!
    \param[in] s list of vertices of the hexahedron
    \param[out] res "reference points" after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  void HexahedronGeomReference::
  FjElem(const VectR3& s, SetPoints<Dimension3>& res,
	 const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureCurve(s, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      {
	FjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = Fi(point) for all nodal points
  void HexahedronGeomReference::
  FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
	      const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
    else
      FjElemNodalLinear(s, res);
  }

  
  //! computes res = Fi(point) for all quadrature points
  void HexahedronGeomReference::
  FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s, res, mesh, nquad);
      }
    else
      FjElemQuadratureLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all dof points
  void HexahedronGeomReference::
  FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
	    const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      FjElemDofLinear(s, res);
  }
  
  
  //! computes res = DFi(point) for all reference points
  /*!
    \param[in] s list of vertices of the hexahedron
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  void HexahedronGeomReference::
  DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	  SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      {
	DFjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = DFi(point) for all nodal points
  void HexahedronGeomReference::
  DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	       SetMatrices<Dimension3>& res,
	       const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
    else
      DFjElemNodalLinear(s, res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points
  void HexahedronGeomReference::
  DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
		    int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemQuadratureLinear(s, res);
  }
  
  
  //! computes res = DFi(point) for all dof points
  void HexahedronGeomReference::
  DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
	     SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemDofLinear(s, res);
  }
  
  
  //! nodal points in the case of straight hexahedron
  void HexahedronGeomReference::
  FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;     res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	res_n.Fill(0);
	for (int j = 0; j < 8; j++)
	  Add(coefFi(i,j), s(j), res_n);
	
	res.SetPointNodal(i, res_n);
      }
  }
  
  
  //! quadrature points in the case of straight hexahedron
  void HexahedronGeomReference::
  FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int j = 0; j < 8; j++)
	  Add(coefFi(this->points_nodal_nd.GetM()+i,j), s(j), res_n);
	
	res.SetPointQuadrature(i, res_n);
      }
  }
  
  
  //! dof points in the case of straight hexahedron
  void HexahedronGeomReference::
  FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n; int offset = this->points_nodal_nd.GetM()+this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int j = 0; j < 8; j++)
	  Add(coefFi(offset+i,j), s(j), res_n);
	    
	res.SetPointDof(i, res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of straight hexahedron
  void HexahedronGeomReference::
  DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	DFjLinear(s, this->points_nodal_nd(i), mat);
	res.SetPointNodal(i, mat);
      }
  }
  
  
  //! transformation DFi for quadrature points in the case of straight hexahedron
  void HexahedronGeomReference::
  DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {  
    Matrix3_3 mat;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	DFjLinear(s, this->points3d(i), mat);
	res.SetPointQuadrature(i, mat);
      }
  }
  
  
  //! transformation DFi for dof points in the case of straight hexahedron
  void HexahedronGeomReference::
  DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	DFjLinear(s, this->points_dof3d(i), mat);
	res.SetPointDof(i, mat);
      }
  }
  
  
  //! transformation Fi for nodal points in the case of curved hexahedron
  /*!
    \param[in] s vertices of the element
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
  */
  void HexahedronGeomReference::
  FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
		   const Mesh<Dimension3>& mesh, int nquad, const Volume& elt) const
  {
    ElementGeomReference<Dimension3>::FjElemNodalCurve(s, res, mesh, nquad);
    
    R3 res_n;
    // for other nodal points, we use coefFi_curve
    for (int node = 6*this->order_geom*this->order_geom+2;
         node < this->points_nodal_nd.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < numberFi_curve(node).GetM(); k++)
	  Add(coefFi_curve(node)(k), res.GetPointNodal(numberFi_curve(node)(k)), res_n);
	
	res.SetPointNodal(node, res_n);
      }
  }
  
  
  //! transformation Fi for quadrature points in the case of curved hexahedron
  void HexahedronGeomReference::
  FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
			const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < coefFi_curve(this->points_nodal_nd.GetM()+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(this->points_nodal_nd.GetM()+i)(node);
	    Add(coef, res.
		GetPointNodal(numberFi_curve(this->points_nodal_nd.GetM()+i)(node)), res_n);
	  }
	res.SetPointQuadrature(i, res_n);
      }
  }
  
  
  //! transformation Fi for dof points in the case of curved hexahedron
  void HexahedronGeomReference::
  FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
		 const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n; int offset = this->points_nodal_nd.GetM() + this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < numberFi_curve(offset+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(offset+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(offset+i)(node)), res_n);
	  }
	res.SetPointDof(i, res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of curved hexahedron
  void HexahedronGeomReference::
  DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		    SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
		    int nquad) const
  {
    Matrix3_3 res_n;
    TinyVector<R3, 3> dfi;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    int r = this->order_geom;
    for (int j1 = 0 ; j1 <= r; j1++)
      for (int j2 = 0; j2 <= r; j2++)
	for (int j3 = 0; j3 <= r; j3++)
	  {
	    int node = this->NumNodes3D(j1,j2,j3);
	    dfi(0).Fill(0); dfi(1).Fill(0); dfi(2).Fill(0);
	    for (int k1 = 0; k1 <= r; k1++)
	      {		    
		int node2 = this->NumNodes3D(k1,j2,j3);
		Add(lob_geom.GradPhi(k1, j1), PTReel.GetPointNodal(node2), dfi(0));
	      }
            
	    for (int k2 = 0; k2 <= r; k2++)
	      {		    
		int node2 = this->NumNodes3D(j1,k2,j3);
		Add(lob_geom.GradPhi(k2, j2), PTReel.GetPointNodal(node2), dfi(1));
	      }
            
	    for (int k3 = 0; k3 <= r; k3++)
	      {		    
		int node2 = this->NumNodes3D(j1,j2,k3);
		Add(lob_geom.GradPhi(k3, j3), PTReel.GetPointNodal(node2), dfi(2));
	      }
            
	    res_n = Matrix3_3(dfi);
	    res.SetPointNodal(node, res_n);
	  }
  }
  
  
  //! transformation DFi for quadrature points in the case of curved hexahedron
  void HexahedronGeomReference::
  DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
			 SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
			 int nquad) const
  {
    Matrix3_3 res_n; int offset = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < numberFi_curve(offset+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(offset+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(offset+i)(node)), res_n);
	  }
	res.SetPointQuadrature(i, res_n);
      }
  }
  
  
  //! transformation DFi for dof points in the case of curved hexahedron
  void HexahedronGeomReference::
  DFjElemDofCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
		  SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
		  int nquad) const
  {
    Matrix3_3 res_n;
    int offset = this->points_nodal_nd.GetM()+this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < numberFi_curve(offset+i).GetM(); node++)
	  {
	    Real_wp coef = coefFi_curve(offset+i)(node);
	    Add(coef, res.GetPointNodal(numberFi_curve(offset+i)(node)), res_n);
	  }
	res.SetPointDof(i, res_n);
      }
  }

    
  /*******************
   * Other functions *
   *******************/
  
  
  //! evaluation of 1-D nodal functions on a point
  /*!
    \param[in] x point where functions should be evaluated
    \param[out] phi values of 1-D nodal functions
   */
  void HexahedronGeomReference::ComputeValuesNodalPhi1D(const Real_wp& x, VectReal_wp& phi) const
  {
    phi.Reallocate(this->order_geom+1);
    for (int i = 0; i <= this->order_geom; i++)
      phi(i) = lob_geom.EvaluatePhi(i, x);
  }
  
  
  //! evaluation of nodal shape functions on a point
  /*!
    \param[in] point_loc point where functions should be evaluated
    \param[out] res values of nodal shape functions
   */
  void HexahedronGeomReference
  ::ComputeValuesPhiNodalRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(this->points_nodal_nd.GetM());
    int r = this->order_geom;
    // computation of 1-D basis functions
    VectReal_wp phix(r+1), phiy(r+1), phiz(r+1);
    for (int i = 0; i <= r; i++)
      {
	phix(i) = lob_geom.EvaluatePhi(i, point_loc(0));
	phiy(i) = lob_geom.EvaluatePhi(i, point_loc(1));
	phiz(i) = lob_geom.EvaluatePhi(i, point_loc(2));
      }
    
    // now tensorization
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	for (int k = 0; k <= r; k++)
	  res(this->NumNodes3D(i, j, k)) = phix(i)*phiy(j)*phiz(k);
  }
  
  
  //! evaluating gradient of nodal shape functions on a point
  /*!
    \param[in] point_loc point where functions should be evaluated
    \param[out] res gradient of nodal shape functions
   */
  void HexahedronGeomReference::
  ComputeGradientPhiNodalRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(this->points_nodal_nd.GetM());
    int r = this->order_geom;
    // computation of 1-D basis functions
    VectReal_wp phix(r+1), phiy(r+1), phiz(r+1);
    VectReal_wp dphix(r+1), dphiy(r+1), dphiz(r+1);
    for (int i = 0; i <= r; i++)
      {
	phix(i) = lob_geom.EvaluatePhi(i, point_loc(0));
	dphix(i) = lob_geom.EvaluatePhiGrad(i, point_loc(0));
	phiy(i) = lob_geom.EvaluatePhi(i, point_loc(1));
	dphiy(i) = lob_geom.EvaluatePhiGrad(i, point_loc(1));
	phiz(i) = lob_geom.EvaluatePhi(i, point_loc(2));
	dphiz(i) = lob_geom.EvaluatePhiGrad(i, point_loc(2));
      }
    
    // now tensorization
    R3 grad;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	for (int k = 0; k <= r; k++)
	  {
	    grad(0) = dphix(i)*phiy(j)*phiz(k);
	    grad(1) = phix(i)*dphiy(j)*phiz(k);
	    grad(2) = phix(i)*phiy(j)*dphiz(k);
	    res(this->NumNodes3D(i, j, k)) = grad;
	  }
    
  }
  
  
  //! evaluating gradient of u on nodal points
  /*!
    \param[in] Uloc values of u on nodal points
    \param[out] gradU gradient of u on nodal points
   */
  template<class T0>
  void HexahedronGeomReference::
  ComputeNodalGradientRefT(const Vector<T0>& Uloc, Vector<TinyVector<T0, 3> >& gradU) const
  {
    int r = this->order_geom;
    gradU.Reallocate(this->points_nodal_nd.GetM());
    for (int i1 = 0; i1 <= r; i1++)
      for (int i2 = 0; i2 <= r; i2++)
	for (int i3 = 0; i3 <= r; i3++)
	  {
	    int i = this->NumNodes3D(i1,i2,i3);
		
	    gradU(i).Fill(0);
	    for (int j = 0; j <= r; j++)
	      {
		gradU(i)(0) += Uloc(this->NumNodes3D(j,i2,i3))*lob_geom.GradPhi(j,i1);		  
		gradU(i)(1) += Uloc(this->NumNodes3D(i1,j,i3))*lob_geom.GradPhi(j,i2);
		gradU(i)(2) += Uloc(this->NumNodes3D(i1,i2,j))*lob_geom.GradPhi(j,i3);
	      }
	  }
  }
  
    
  void HexahedronGeomReference
  ::ComputeNodalGradientRef(const Vector<Real_wp>& Uloc, VectR3& gradU) const
  {
    ComputeNodalGradientRefT(Uloc, gradU);
  }
  
  
  void HexahedronGeomReference
  ::ComputeNodalGradientRef(const Vector<Complex_wp>& Uloc, Vector<R3_Complex_wp>& gradU) const
  {
    ComputeNodalGradientRefT(Uloc, gradU);
  }


  //! displays informations class HexahedronGeomReference
  ostream& operator <<(ostream& out, const HexahedronGeomReference& e)
  {
    out << static_cast<const ElementGeomReference<Dimension3>&>(e);
    out<<"coefficients used for Fi"<<endl<<e.coefFi<<endl;
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_GEOM_REFERENCE_CXX
#endif
