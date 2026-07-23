#ifndef MONTJOIE_FILE_WEDGE_GEOM_REFERENCE_CXX

namespace Montjoie
{
  
  /****************************
   * Initialization functions *
   ****************************/
  

  //! default constructor
  WedgeGeomReference::WedgeGeomReference() : ElementGeomReference<Dimension3>()
  {
    // normales of the five faces
    this->normale.Reallocate(5);
    this->normale(0).Init(0,0,-1.0); 
    this->normale(1).Init(0,-1,0);
    this->normale(2).Init(1,1,0);
    this->normale(3).Init(-1,0,0); 
    this->normale(4).Init(0,0,1.0);

    this->tangente_loc_x.Reallocate(5);
    this->tangente_loc_y.Reallocate(5);
    this->tangente_loc_x(0).Init(1,0,0);
    this->tangente_loc_y(0).Init(0,1,0);
    this->tangente_loc_x(1).Init(1,0,0);
    this->tangente_loc_y(1).Init(0,0,1);
    this->tangente_loc_x(2).Init(-1,1,0);
    this->tangente_loc_y(2).Init(0,0,1);
    this->tangente_loc_x(3).Init(0,1,0);
    this->tangente_loc_y(3).Init(0,0,1);
    this->tangente_loc_x(4).Init(1,0,0);
    this->tangente_loc_y(4).Init(0,1,0);
    
    this->hybrid_type_elt = 2;
    this->nb_vertices_elt = 6;
    this->nb_edges_elt = 9;
    this->nb_boundaries_elt = 5;

    this->is_local_face_quad.Reallocate(5);
    this->is_local_face_quad.Fill(true);
    this->is_local_face_quad(0) = false;
    this->is_local_face_quad(4) = false;
  }


  size_t WedgeGeomReference::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension3>::GetMemorySize();
    taille += function_geom_tri.GetMemorySize();
    taille += function_quad.GetMemorySize();
    taille += EdgesNodal.GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    taille += coefFi_curve.GetMemorySize();
    taille += coefFi.GetMemorySize() + coefDFi_dx.GetMemorySize() + coefDFi_dy.GetMemorySize()
      + coefDFi_dz.GetMemorySize();
    taille += NumNodesTri.GetMemorySize();
    taille += CoordinateNodes.GetMemorySize();
    return taille;
  }

  
  //! constructing finite element
  void WedgeGeomReference::
  ConstructFiniteElement(int rgeom)
  {
    this->order_geom = rgeom;
    int r = rgeom;
    
    // nodal functions
    ConstructNodalShapeFunctions(rgeom);

    ConstructOrthogonalFunctions(r);
  }


  //! constructs Legendre polynomials
  void WedgeGeomReference::ConstructOrthogonalFunctions(int r)
  {
    // Legendre functions
    GetJacobiPolynomial(LegendrePolynom, r, Real_wp(0), Real_wp(0));

    CoefLegendre.Reallocate(r+1); CoefLegendre.Fill(0);
    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, r);
    for (int q = 0; q <= r; q++)
      {
	VectReal_wp Pn;
	EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*xi(q) - 1.0, Pn);
	for (int i = 0; i <= r; i++)
	  CoefLegendre(i) += omega(q)*Pn(i)*Pn(i);
      }
    
    // we keep 1/sqrt(omega)
    for (int i = 0; i <= r; i++)
      CoefLegendre(i) = 1.0/sqrt(CoefLegendre(i));    
  }
  
  
  //! construction of nodal points, which are used for interpolation  
  void WedgeGeomReference::ConstructNodalShapeFunctions(int r)
  {
    // clearing previous points
    this->points_nodal1d.Clear(); this->points_nodal2d_quad.Clear();
    this->points_nodal2d_tri.Clear(); this->points_nodal_nd.Clear();
   
    MeshNumbering<Dimension3>::
      ConstructPrismaticNumbering(r, this->NumNodes3D, this->CoordinateNodes);
    
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(r, this->NumNodes2D_quad, this->CoordinateNodes2D_quad);
    MeshNumbering<Dimension2>::
      ConstructTriangularNumbering(r, this->NumNodes2D_tri, this->CoordinateNodes2D_tri);

    this->EdgesNodal.Reallocate(r+1, 9);
    for (int i = 0; i <= r; i++)
      {
	this->EdgesNodal(i, 0) = this->NumNodes3D(i, 0, 0);
	this->EdgesNodal(i, 1) = this->NumNodes3D(r-i, i, 0);
	this->EdgesNodal(i, 2) = this->NumNodes3D(0, i, 0);
	this->EdgesNodal(i, 3) = this->NumNodes3D(0, 0, i);
	this->EdgesNodal(i, 4) = this->NumNodes3D(r, 0, i);
	this->EdgesNodal(i, 5) = this->NumNodes3D(0, r, i);
	this->EdgesNodal(i, 6) = this->NumNodes3D(i, 0, r);
	this->EdgesNodal(i, 7) = this->NumNodes3D(r-i, i, r);
	this->EdgesNodal(i, 8) = this->NumNodes3D(0, i, r);
      }
    
    this->FacesNodal.Reallocate((r+1)*(r+1), 5);
    this->FacesNodal.Fill(-1);
    // quadrangular faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  int node = this->NumNodes2D_quad(i,j);
	  this->FacesNodal(node, 1) = this->NumNodes3D(i, 0, j);
	  this->FacesNodal(node, 2) = this->NumNodes3D(r-i, i, j);
	  this->FacesNodal(node, 3) = this->NumNodes3D(0, i, j);
	}
    
    // triangular faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node = this->NumNodes2D_tri(i,j);
	  this->FacesNodal(node,0) = this->NumNodes3D(i, j, 0);
	  this->FacesNodal(node,4) = this->NumNodes3D(i, j, r);
	}
        
    lob_geom.ConstructQuadrature(r, lob_geom.QUADRATURE_LOBATTO);
 
    ConstructLobattoPoints(r, this->points_nodal1d, this->points_nodal2d_tri,
                           this->NumNodes2D_tri,
                           this->points_nodal2d_quad, this->NumNodes2D_quad,
                           this->points_nodal_nd, this->NumNodes3D);
    
    function_geom_tri.ConstructFiniteElement(r);
    
    this->NumNodesTri.Reallocate(this->points_nodal2d_tri.GetM(), r+1);
    this->CoordinateNodes.Reallocate(this->points_nodal_nd.GetM(), 2);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= r; k++)
          {
            int node = this->NumNodes3D(i, j, k);
            this->CoordinateNodes(node, 0) = this->NumNodes2D_tri(i, j);
            this->CoordinateNodes(node, 1) = k;
            this->NumNodesTri(this->NumNodes2D_tri(i, j), k) = node;
          }
        
    if (r > 2)
      {
        VectR2 points_nodal2d_inside_tri((r-1)*(r-2)/2);
        for (int i = 0; i < points_nodal2d_inside_tri.GetM(); i++)
          points_nodal2d_inside_tri(i) = this->points_nodal2d_tri(i+3*r);
        
        MeshNumbering<Dimension3>::
          GetRotationTriangularFace(points_nodal2d_inside_tri, this->FacesNodal_Rotation_Tri);
      }
    
    if (r > 1)
      {
        Matrix<int> NumNodes2D_InsideQuad(r-1, r-1);
        for (int i = 1;  i < r; i++)
          for (int j = 1;  j < r; j++)
            NumNodes2D_InsideQuad(i-1, j-1) = (i-1)*(r-1) + j-1;
        
	MeshNumbering<Dimension3>::
	  GetRotationQuadrilateralFace(NumNodes2D_InsideQuad, this->FacesNodal_Rotation_Quad);
      }
    
    VectR3 grad_phi;
    this->GradientPhi_Nodal.Reallocate(this->points_nodal_nd.GetM(), this->points_nodal_nd.GetM());
    for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
      {
        ComputeGradientPhiNodalRef(this->points_nodal_nd(j), grad_phi);
        
        for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
          this->GradientPhi_Nodal(k, j) = grad_phi(k);
      }
    
    ComputeCurvedTransformation();
    ComputeCoefficientTransformation();
  }
  
  
  //! nodal points are set on Pascal's Wedge, regularly spaced points  
  void WedgeGeomReference::
  ConstructRegularPoints(int r, VectReal_wp& points1d_, VectR2& points2d_tri_,
                         const Matrix<int>& NumNodes2D_tri_,
                         VectR2& points2d_quad_, const Matrix<int>& NumNodes2D_quad_,
                         VectR3& points3d_, const Array3D<int>& NumNodes3D_)
  {
    abort();
  }
  
  
  //! constructing points on the wedge   
  void WedgeGeomReference::
  ConstructLobattoPoints(int r, VectReal_wp& points1d_, VectR2& points2d_tri_,
                         const Matrix<int>& NumNodes2D_tri_,
                         VectR2& points2d_quad_, const Matrix<int>& NumNodes2D_quad_,
                         VectR3& points3d_, const Array3D<int>& NumNodes3D_)
  {    
    // we construct lobatto points
    TriangleGeomReference::ConstructLobattoPoints(r, 0, points1d_, points2d_tri_);
    
    points2d_quad_.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	points2d_quad_(NumNodes2D_quad_(i, j)) = R2(points1d_(i), points1d_(j));
    
    points3d_.Reallocate(points2d_tri_.GetM()*(r+1));
    for (int k = 0; k <= r; k++)
      for (int i = 0; i <= r; i++)
	for (int j = 0; j <= r-i; j++)
        {
	  int node = NumNodes3D_(i, j, k);
	  int node2D = NumNodes2D_tri_(i, j);
          points3d_(node)(0) = points2d_tri_(node2D)(0);
	  points3d_(node)(1) = points2d_tri_(node2D)(1);
	  points3d_(node)(2) = points1d_(k);
	}    
  }

  
  //! computation of coefFicients to compute quickly Fi and DFi  
  void WedgeGeomReference::ComputeCurvedTransformation()
  {
    coefFi_curve.Reallocate(this->points_nodal_nd.GetM(), this->points_nodal_nd.GetM());
    coefFi_curve.Fill(0); VectReal_wp phi;
    // use of transformation Fi defined in Solin book for inside points
    TinyVector<Real_wp, 3> lambda; Real_wp x, y, z, zeta, coef, coefs;
    for (int i = 4*this->order_geom*this->order_geom+2; i < this->points_nodal_nd.GetM(); i++)
      {
	lambda(1) = this->points_nodal_nd(i)(0);
	lambda(2) = this->points_nodal_nd(i)(1);
	x = lambda(1); y = lambda(2);
	z = this->points_nodal_nd(i)(2);
	lambda(0) = 1.0 - lambda(1) - lambda(2);
	
	// part due to vertices
	coefFi_curve(i, 0) = lambda(0)*(1-z);
	coefFi_curve(i, 1) = lambda(1)*(1-z);
	coefFi_curve(i, 2) = lambda(2)*(1-z);
	coefFi_curve(i, 3) = lambda(0)*z;
	coefFi_curve(i, 4) = lambda(1)*z;
	coefFi_curve(i, 5) = lambda(2)*z;
        
	// part due to edges on triangles
	for (int j = 0; j < 9; j++)
	  if ((j < 3)||(j>=6))
	    {
	      int n1 = MeshNumbering<Dimension3>::FirstExtremityEdge(2, j);
	      int n2 = MeshNumbering<Dimension3>::SecondExtremityEdge(2, j);
	      zeta = 0.5*(lambda(n2%3)-lambda(n1%3)) + 0.5; 
	      coef = lambda(n1%3)*lambda(n2%3)/(zeta*(1.0-zeta));
	      if (j < 3)
		coef *= 1.0-z;
	      else
		coef *= z;
	      
	      // contribution of points of the edge
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, this->EdgesNodal(k, j) ) += lob_geom.EvaluatePhi(k, zeta)*coef;
	      
	      // contribution of the extremities of the edge
	      coefFi_curve(i, this->EdgesNodal(0,j) ) -= (1.0-zeta)*coef;
	      coefFi_curve(i, this->EdgesNodal(this->order_geom,j) ) -= zeta*coef;
	    }
	
	// part due to edges on quadrilaterals
	for (int j = 3; j < 6; j++)
	  {
	    coef = lambda(j-3);
	    // contribution of points of the edge
	    for (int k = 0; k <= this->order_geom; k++)
	      coefFi_curve(i, this->EdgesNodal(k, j) ) += lob_geom.EvaluatePhi(k, z)*coef;
	    
	    // contribution of the extremities of the edge
	    coefFi_curve(i, this->EdgesNodal(0,j) ) -= (1.0-z)*coef;
	    coefFi_curve(i, this->EdgesNodal(this->order_geom,j) ) -= z*coef;
	  }
	
	
	// part due to first triangular face
	coef = (1.0-z);
	function_geom_tri.ComputeValuesPhiNodalRef(R2(x, y), phi);
	for (int k = 0; k < this->points_nodal2d_tri.GetM(); k++)
	  coefFi_curve(i, this->FacesNodal(k, 0)) += coef*phi(k);
	
	// part due to the edges of the face
	for (int j = 0; j < 3; j++)
	  {
	    zeta = 0.5*(lambda((j+1)%3)-lambda(j)) + 0.5; 
	    coefs = lambda((j+1)%3)*lambda(j)/(zeta*(1.0-zeta));
	    if (j == 2)
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, this->EdgesNodal(this->order_geom-k, j))
                  -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	    else
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, this->EdgesNodal(k, j)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	    // contribution of the extremities
	    coefFi_curve(i, j) += (1.0-zeta)*coef*coefs;
	    coefFi_curve(i, (j+1)%3) += zeta*coef*coefs;
	  }
	
	coefFi_curve(i, 0) -= coef*lambda(0);
	coefFi_curve(i, 1) -= coef*lambda(1);
	coefFi_curve(i, 2) -= coef*lambda(2);
		
	// part due to second triangular face
	coef = z;
	function_geom_tri.ComputeValuesPhiNodalRef(R2(x, y), phi);
	for (int k = 0; k < this->points_nodal2d_tri.GetM(); k++)
	  coefFi_curve(i, this->FacesNodal(k, 4)) += coef*phi(k);
	
	// part due to the edges of the face
	for (int j = 0; j < 3; j++)
	  {
	    zeta = 0.5*(lambda((j+1)%3)-lambda(j)) + 0.5; 
	    coefs = lambda((j+1)%3)*lambda(j)/(zeta*(1.0-zeta));
	    if (j == 2)
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, this->EdgesNodal(this->order_geom-k, j+6))
                  -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	    else
	      for (int k = 0; k <= this->order_geom; k++)
		coefFi_curve(i, this->EdgesNodal(k, j+6)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	    // contribution of the extremities
	    coefFi_curve(i, j+3) += (1.0-zeta)*coef*coefs;
	    coefFi_curve(i, (j+1)%3+3) += zeta*coef*coefs;
	  }
	
	coefFi_curve(i, 3) -= coef*lambda(0);
	coefFi_curve(i, 4) -= coef*lambda(1);
	coefFi_curve(i, 5) -= coef*lambda(2);
	
	// part due to first quadrilateral face
	zeta = x + 0.5*y;
	coef = lambda(0)*lambda(1)/(zeta*(1.0-zeta));
	for (int j = 0; j <= this->order_geom; j++)
	  for (int k = 0; k <= this->order_geom; k++)
	    coefFi_curve(i, this->FacesNodal(this->NumNodes2D_quad(j, k), 1))
              += coef*lob_geom.EvaluatePhi(j, zeta)*lob_geom.EvaluatePhi(k, z);
	
	for (int k = 0; k <= this->order_geom; k++)
	  {
	    coefFi_curve(i, this->EdgesNodal(k, 0)) -= coef*(1.0-z)*lob_geom.EvaluatePhi(k, zeta);
	    coefFi_curve(i, this->EdgesNodal(k, 4)) -= coef*zeta*lob_geom.EvaluatePhi(k, z);
	    coefFi_curve(i, this->EdgesNodal(k, 6)) -= coef*z*lob_geom.EvaluatePhi(k, zeta);
	    coefFi_curve(i, this->EdgesNodal(k, 3)) -= coef*(1.0-zeta)*lob_geom.EvaluatePhi(k, z);
	  }
	
	coefFi_curve(i, 0) += coef*(1.0-zeta)*(1.0-z);
	coefFi_curve(i, 1) += coef*zeta*(1.0-z);
	coefFi_curve(i, 4) += coef*zeta*z;
	coefFi_curve(i, 3) += coef*(1.0-zeta)*z;

	// part due to second quadrilateral face
	zeta = 0.5*(y-x) + 0.5;
	coef = lambda(1)*lambda(2)/(zeta*(1.0-zeta));
	for (int j = 0; j <= this->order_geom; j++)
	  for (int k = 0; k <= this->order_geom; k++)
	    coefFi_curve(i, this->FacesNodal(this->NumNodes2D_quad(j, k), 2))
              += coef*lob_geom.EvaluatePhi(j, zeta)*lob_geom.EvaluatePhi(k, z);
	
	for (int k = 0; k <= this->order_geom; k++)
	  {
	    coefFi_curve(i, this->EdgesNodal(k, 1)) -= coef*(1.0-z)*lob_geom.EvaluatePhi(k, zeta);
	    coefFi_curve(i, this->EdgesNodal(k, 5)) -= coef*zeta*lob_geom.EvaluatePhi(k, z);
	    coefFi_curve(i, this->EdgesNodal(k, 7)) -= coef*z*lob_geom.EvaluatePhi(k, zeta);
	    coefFi_curve(i, this->EdgesNodal(k, 4)) -= coef*(1.0-zeta)*lob_geom.EvaluatePhi(k, z);
	  }
	
	coefFi_curve(i, 1) += coef*(1.0-zeta)*(1.0-z);
	coefFi_curve(i, 2) += coef*zeta*(1.0-z);
	coefFi_curve(i, 5) += coef*zeta*z;
	coefFi_curve(i, 4) += coef*(1.0-zeta)*z;
	
	// part due to third quadrilateral face
	zeta = y + 0.5*x;
	coef = lambda(0)*lambda(2)/(zeta*(1.0-zeta));
	for (int j = 0; j <= this->order_geom; j++)
	  for (int k = 0; k <= this->order_geom; k++)
	    coefFi_curve(i, this->FacesNodal(this->NumNodes2D_quad(j, k), 3))
              += coef*lob_geom.EvaluatePhi(j, zeta)*lob_geom.EvaluatePhi(k, z);
	
	for (int k = 0; k <= this->order_geom; k++)
	  {
	    coefFi_curve(i, this->EdgesNodal(k, 2)) -= coef*(1.0-z)*lob_geom.EvaluatePhi(k, zeta);
	    coefFi_curve(i, this->EdgesNodal(k, 5)) -= coef*zeta*lob_geom.EvaluatePhi(k, z);
	    coefFi_curve(i, this->EdgesNodal(k, 8)) -= coef*z*lob_geom.EvaluatePhi(k, zeta);
	    coefFi_curve(i, this->EdgesNodal(k, 3)) -= coef*(1.0-zeta)*lob_geom.EvaluatePhi(k, z);
	  }
	
	coefFi_curve(i, 0) += coef*(1.0-zeta)*(1.0-z);
	coefFi_curve(i, 2) += coef*zeta*(1.0-z);
	coefFi_curve(i, 5) += coef*zeta*z;
	coefFi_curve(i, 3) += coef*(1.0-zeta)*z;
	
      }

    const TriangleGeomReference& tri = this->GetTriangularSurfaceFiniteElement();
    this->InitTriangularCurvature(tri);    
  }
  

  //! computation of coefFicients to compute quickly Fi and DFi    
  void WedgeGeomReference::ComputeCoefficientTransformation()
  {
    
    VectR3 all_points = this->points_nodal_nd;
    all_points.PushBack(this->points3d); all_points.PushBack(this->points_dof3d);   
    int nb_points_all = all_points.GetM();
    
    // straight Wedge
    coefFi.Reallocate(nb_points_all,6);
    coefDFi_dx.Reallocate(nb_points_all,6);
    coefDFi_dy.Reallocate(nb_points_all,6);
    coefDFi_dz.Reallocate(nb_points_all,6);
    VectReal_wp phi; VectR3 dphi;
    
    for (int i = 0; i < nb_points_all; i++) 
      {
        phi.Reallocate(6); ComputeValuesPhiFirstOrder(all_points(i), phi);   
        dphi.Reallocate(6); ComputeGradientPhiFirstOrder(all_points(i), dphi);  
        for (int j = 0; j < 6; j++)
	  {
	    coefFi(i,j) = phi(j);
	    coefDFi_dx(i,j) = dphi(j)(0);
	    coefDFi_dy(i,j) = dphi(j)(1);
	    coefDFi_dz(i,j) = dphi(j)(2);
	  }
      }
    
    coefFi_curve.Resize(nb_points_all, this->points_nodal_nd.GetM());
    // for other points, we use lagrangian interpolation
    for (int i = this->points_nodal_nd.GetM(); i < nb_points_all; i++)
      {
	ComputeValuesPhiNodalRef(all_points(i), phi);
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  coefFi_curve(i, j) = phi(j);
      }
    
  }
    

  void WedgeGeomReference
  ::GetGradient3D_FromGradient2D(int num_loc, int i, const VectReal_wp& nabla_nx,
				 const VectReal_wp& nabla_ny, const VectReal_wp& nabla_nz,
				 R3& d_nx, R3& d_ny, R3& d_nz) const
  {
    switch (num_loc)
      {
      case 0:
      case 4:
	d_nx(0) = nabla_nx(2*i); d_nx(1) = nabla_nx(2*i+1); d_nx(2) = 0;
	d_ny(0) = nabla_ny(2*i); d_ny(1) = nabla_ny(2*i+1); d_ny(2) = 0;
	d_nz(0) = nabla_nz(2*i); d_nz(1) = nabla_nz(2*i+1); d_nz(2) = 0;
	break;
      case 1:
	d_nx(0) = nabla_nx(2*i); d_nx(1) = 0; d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = nabla_ny(2*i); d_ny(1) = 0; d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = nabla_nz(2*i); d_nz(1) = 0; d_nz(2) = nabla_nz(2*i+1);
	break;
      case 2 :
	// DF = [-1 0; 1 0; 0 1] and DF*-1 = [-1/2 0; 1/2 0; 0 1]
	d_nx(0) = -0.5*nabla_nx(2*i); d_nx(1) = -d_nx(0); d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = -0.5*nabla_ny(2*i); d_ny(1) = -d_ny(0); d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = -0.5*nabla_nz(2*i); d_nz(1) = -d_nz(0); d_nz(2) = nabla_nz(2*i+1);
	break;
      case 3 :
	d_nx(0) = 0; d_nx(1) = nabla_nx(2*i); d_nx(2) = nabla_nx(2*i+1);
	d_ny(0) = 0; d_ny(1) = nabla_ny(2*i); d_ny(2) = nabla_ny(2*i+1);
	d_nz(0) = 0; d_nz(1) = nabla_nz(2*i); d_nz(2) = nabla_nz(2*i+1);
	break;
      }
  }


  //! returns a tangential vector on a face of the hexahedron from a 2-D vector on the square
  R3 WedgeGeomReference::GetTangentialVector(int num_loc, const R2& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
      case 4:
	return R3(vec_u(0), vec_u(1), Real_wp(0));
      case 1:
	return R3(vec_u(0), Real_wp(0), vec_u(1));
      case 2 :
	// DF = [-1 0; 1 0; 0 1] and DF*-1 = [-1/2 0; 1/2 0; 0 1]
	return R3(-0.5*vec_u(0), 0.5*vec_u(0), vec_u(1));
      case 3 :
	return R3(Real_wp(0), vec_u(0), vec_u(1));
      }

    return R3();
  }


  //! transpose operation of GetTangentialVector
  R2 WedgeGeomReference::TransposeTangentialVector(int num_loc, const R3& vec_u) const
  {
  switch (num_loc)
      {
      case 0:
      case 4:
	return R2(vec_u(0), vec_u(1));
      case 1:
	return R2(vec_u(0), vec_u(2));
      case 2 :
        return R2(-0.5*(vec_u(0) - vec_u(1)), vec_u(2));
      case 3 :
	return R2(vec_u(1), vec_u(2));
      }

    return R2();
  }


  //! transpose operation of GetTangentialVector
  R2_Complex_wp WedgeGeomReference::TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const
  {
  switch (num_loc)
      {
      case 0:
      case 4:
	return R2_Complex_wp(vec_u(0), vec_u(1));
      case 1:
	return R2_Complex_wp(vec_u(0), vec_u(2));
      case 2 :
        return R2_Complex_wp(-0.5*(vec_u(0) - vec_u(1)), vec_u(2));
      case 3 :
	return R2_Complex_wp(vec_u(1), vec_u(2));
      }

    return R2_Complex_wp();
  }
  
  
  /****************
   * Fj transform *
   ****************/
  

  //! transformation Fi in the case of straight Wedge
  void WedgeGeomReference::FjLinear(const VectR3& s,const R3& point,R3& res) const
  {
    // Fi = phi1*S0 + phi2*S1 + phi3*S2 + phi4*S3 + phi5*S4 + phi6*S5
    res.Zero();
    VectReal_wp phi; ComputeValuesPhiFirstOrder(point,phi);

    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 3; j++)
        res(j) += phi(i)*s(i)(j);
  }
  
  
  //! transformation DFi in the case of straight Wedge
  void WedgeGeomReference::DFjLinear(const VectR3& s,const R3& point,Matrix3_3& res) const
  {
    res.Zero();
    VectR3 dphi; ComputeGradientPhiFirstOrder(point,dphi);

    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 3; j++)
	{
	  res(j,0) += dphi(i)(0)*s(i)(j);
	  res(j,1) += dphi(i)(1)*s(i)(j);
	  res(j,2) += dphi(i)(2)*s(i)(j);
	}
    
  }

  
  //! returns the minimal length of the element
  Real_wp WedgeGeomReference
  ::GetMinimalSize(const VectR3& s) const
  {
    Real_wp h = Distance(s(0), s(1));
    h = min(h, Distance(s(0), s(2)));
    h = min(h, Distance(s(0), s(3)));
    h = min(h, Distance(s(1), s(2)));
    h = min(h, Distance(s(1), s(4)));
    h = min(h, Distance(s(2), s(5)));
    h = min(h, Distance(s(3), s(4)));
    h = min(h, Distance(s(3), s(5)));
    h = min(h, Distance(s(4), s(5)));
    return h;
  }
  

  //! returns true if the point is outside the unit Wedge
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
   */  
  bool WedgeGeomReference
  ::OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const
  {
    Real_wp Rmax = max(abs(s(0)(0)), abs(s(0)(1))); Rmax = max(Rmax, abs(s(0)(2)));
    Rmax = max(Rmax, abs(s(1)(0))); Rmax = max(Rmax, abs(s(1)(1))); Rmax = max(Rmax, abs(s(1)(2)));
    Rmax = max(Rmax, abs(s(2)(0))); Rmax = max(Rmax, abs(s(2)(1))); Rmax = max(Rmax, abs(s(2)(2)));
    Rmax = max(Rmax, abs(s(3)(0))); Rmax = max(Rmax, abs(s(3)(1))); Rmax = max(Rmax, abs(s(3)(2)));
    Rmax = max(Rmax, abs(s(4)(0))); Rmax = max(Rmax, abs(s(4)(1))); Rmax = max(Rmax, abs(s(4)(2)));
    Rmax = max(Rmax, abs(s(5)(0))); Rmax = max(Rmax, abs(s(5)(1))); Rmax = max(Rmax, abs(s(5)(2)));
    
    Real_wp dx = s(0).Distance(s(1)), dy = s(0).Distance(s(2)), dz = s(0).Distance(s(3));
    if ((Rmax == Real_wp(0)) || (dx == Real_wp(0)) || (dy == Real_wp(0)) || (dz == Real_wp(0)))
      return true;
    
    Real_wp one(1);
    if (dx*Xn(0) < -epsilon*Rmax)
      return true;

    if (dy*Xn(1) < -epsilon*Rmax)
      return true;

    if (dz*Xn(2) < -epsilon*Rmax)
      return true;

    if (dz*(Xn(2) - one) > epsilon*Rmax)
      return true;

    if (min(dx, dy)*(Xn(0) + Xn(1) - one) > epsilon*Rmax)
      return true;
    
    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
   */
  Real_wp WedgeGeomReference::GetDistanceToBoundary(const R3& pointloc) const
  {
    Real_wp distance = pointloc(2); 
    distance = min(distance, pointloc(0));
    distance = min(distance, pointloc(1));  
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)));
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(0)-pointloc(1)));    
    return distance;
  }
  
  
  Real_wp WedgeGeomReference::GetDistanceToBoundary(const R3& pointloc, int n) const
  {
    cout << "Not implemented" << endl;
    abort();
    return Real_wp(0);
  }


  //! project the local point on the boundary, if outside the unit Wedge
  int WedgeGeomReference::ProjectPointOnBoundary(R3& pointloc) const
  {
    pointloc(2) = max( Real_wp(0), pointloc(2) );
    Real_wp zero(0),one(1);
    if (pointloc(0)<zero)
      pointloc(0) = Real_wp(pointloc(2)-one);
    else if (pointloc(1)<zero)
      pointloc(0) = Real_wp(one-pointloc(2));
    if (one-pointloc(1)<zero)
      pointloc(1) = Real_wp(pointloc(2)-one);
    else if (Real_wp(1)-pointloc(0)-pointloc(1)<zero)
      pointloc(1) = Real_wp(one-pointloc(2));
    
    return -1;
  }
  

  //! evaluates lowest order shape functions at a point in the element  
  void WedgeGeomReference
  ::ComputeValuesPhiFirstOrder(const R3& point, VectReal_wp & phi) const
  {
    Real_wp x = point(0), y = point(1), z = point(2);
    phi.Reallocate(6);
    Real_wp one(1);
    
    phi(0) = (one-x-y)*(one-z);
    phi(1) = x*(one-z);
    phi(2) = y*(one-z);
    phi(3) = (one-x-y)*z;
    phi(4) = x*z;
    phi(5) = y*z;

  }
  

  //! evaluates gradient of lowest order shape functions at a point in the element      
  void WedgeGeomReference::ComputeGradientPhiFirstOrder(const R3& point, VectR3& dphi) const
  {
    Real_wp x = point(0), y = point(1), z = point(2);
    dphi.Reallocate(6);
    Real_wp zero(0),one(1);

    dphi(0)(0) = z-one;   dphi(1)(0) = one-z; dphi(2)(0) = zero;  dphi(3)(0) = -z;
    dphi(4)(0) = z;    dphi(5)(0) = zero;
    dphi(0)(1) = z-one;   dphi(1)(1) = zero;  dphi(2)(1) = one-z; dphi(3)(1) = -z;
    dphi(4)(1) = zero; dphi(5)(1) = z;
    dphi(0)(2) = x+y-one; dphi(1)(2) = -x;    dphi(2)(2) = -y;    dphi(3)(2) = one-x-y;
    dphi(4)(2) = x;    dphi(5)(2) = y;
  
  }

  
  //! returns 3-D coordinates of a point on a face
  /*!
    \param[in] num_loc local face number in the quad
    \param[in] point_loc coordinate on the face
    \param[out] res local coordinates of the point in the unit square
  */
  void WedgeGeomReference
  ::GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const
  {
    switch(num_loc)
      {
      case 0:
	res.Init(point_loc(0), point_loc(1), 0);break;
      case 1:
	res.Init(point_loc(0),Real_wp(0),point_loc(1));break;
      case 2:
	res.Init(Real_wp(1)-point_loc(0),point_loc(0),point_loc(1));break;
      case 3:
	res.Init(Real_wp(0),point_loc(0),point_loc(1));break;
      case 4:
	res.Init(point_loc(0),point_loc(1),Real_wp(1));break;
      }
  }


  //! computes 2-D coordinates of a point on a wedge
  /*!
    \param[in] num_loc local face number in the quad
    \param[out] point_loc coordinate on the face
    \param[in] pt local coordinates of the point in the unit prism
  */
  void WedgeGeomReference
  ::GetLocalCoordOnBoundary(int num_loc, const R3& pt, R2& point_loc) const
  {
    switch(num_loc)
      {
      case 0:
	point_loc.Init(pt(0), pt(1)); break;
      case 1:
	point_loc.Init(pt(0), pt(2)); break;
      case 2:
	point_loc.Init(pt(1), pt(2)); break;
      case 3:
	point_loc.Init(pt(0), pt(2)); break;
      case 4:
	point_loc.Init(pt(0), pt(1)); break;
      }
  }

  
  //! transformation Fi in the case of curved Wedge
  void WedgeGeomReference
  ::FjCurve(const SetPoints<Dimension3>& PTReel,const R3& pointloc,R3& res) const
  {
    res.Fill(0); VectReal_wp phi;
    ComputeValuesPhiNodalRef(pointloc, phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      Add(phi(node), PTReel.GetPointNodal(node), res);
  }
  
  
  //! transformation DFi in the case of curved Wedge
  void WedgeGeomReference::
  DFjCurve(const SetPoints<Dimension3>& PTReel,const R3& pointloc,Matrix3_3& res) const
  {
    R3 grad; VectR3 grad_phi;
    res.Fill(0); ComputeGradientPhiNodalRef(pointloc, grad_phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      {
	grad = grad_phi(node);
	res(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	res(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
	res(2,0) += PTReel.GetPointNodal(node)(2)*grad(0);
	res(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	res(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
	res(2,1) += PTReel.GetPointNodal(node)(2)*grad(1);
	res(0,2) += PTReel.GetPointNodal(node)(0)*grad(2);
	res(1,2) += PTReel.GetPointNodal(node)(1)*grad(2);
	res(2,2) += PTReel.GetPointNodal(node)(2)*grad(2);
      }
  }


  //! computing decomposition of jacobian in polynomials  
  void WedgeGeomReference
  ::ComputeCoefJacobian(const VectR3& s, VectReal_wp& CoefJacobian) const
  {
    CoefJacobian.Reallocate(7);
    R3 A1 = -s(0) + s(1);
    R3 A2 = -s(0) + s(2);
    R3 A3 = -s(0) + s(3);
    R3 C1 = s(0) - s(1) - s(3) + s(4);
    R3 C2 = s(0) - s(2) - s(3) + s(5);
    
    CoefJacobian(0) = Determinant(A1,A2,A3);
    CoefJacobian(1) = Determinant(A1,A2,C1);
    CoefJacobian(2) = Determinant(A1,A2,C2);
    CoefJacobian(3) = Determinant(C1,A2,A3) + Determinant(A1,C2,A3);
    CoefJacobian(4) = Determinant(A1,C2,C1);
    CoefJacobian(5) = Determinant(C1,A2,C2);
    CoefJacobian(6) = Determinant(C1,C2,A3);
  }
  
  
  /**********************
   * FjElem and DFjElem *
   **********************/
  
  
  //! computes res = Fi(point) for all reference points
  /*!
    \param[in] s list of vertices of the element
    \param[out] res "reference points" after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  void WedgeGeomReference::FjElem(const VectR3& s, SetPoints<Dimension3>& res,
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
  void WedgeGeomReference
  ::FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
                const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
    else
      FjElemNodalLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all quadrature points  
  void WedgeGeomReference
  ::FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
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
  void WedgeGeomReference
  ::FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
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
    \param[in] s list of vertices of the element
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */  
  void WedgeGeomReference
  ::DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
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
  void WedgeGeomReference
  ::DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                 SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
    else
      DFjElemNodalLinear(s, res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points  
  void WedgeGeomReference
  ::DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                      SetMatrices<Dimension3>& res, 
                      const Mesh<Dimension3>& mesh, int nquad) const
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
  void WedgeGeomReference
  ::DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
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
  
  
  //! nodal points in the case of straight Wedge  
  void WedgeGeomReference::FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int n = 0; n < this->points_nodal_nd.GetM(); n++)
      {
        res_n.Fill(0);
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 3; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);

	res.SetPointNodal(n,res_n);
      }
  }
  
  
  //! quadrature points in the case of straight Wedge  
  void WedgeGeomReference
  ::FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {  
    R3 res_n;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int n = this->points_nodal_nd.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM(); n++)
      {
        res_n.Fill(0);
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 3; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointQuadrature(n-this->points_nodal_nd.GetM(),res_n);
      }
  }
  
  
  //! dof points in the case of straight Wedge
  void WedgeGeomReference::FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;     res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int n = this->points_nodal_nd.GetM()+this->points3d.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM()+this->points_dof3d.GetM(); n++)
      {
        res_n.Fill(0);
	for (int i = 0; i < 6;i++)
          for (int j = 0; j < 3;j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointDof(n-this->points_nodal_nd.GetM()-this->points3d.GetM(),res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of straight Wedge  
  void WedgeGeomReference::DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int n = 0; n < this->points_nodal_nd.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 6; i++)
         for (int j = 0; j < 3; j++)
	   {
	     mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	     mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	     mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	   }
        res.SetPointNodal(n,mat);	  
      }
  }
  
  
  //! transformation DFi for quadrature points in the case of straight Wedge  
  void WedgeGeomReference
  ::DFjElemQuadratureLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int n = this->points_nodal_nd.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 6; i++)
	  for (int j = 0; j < 3; j++)
	    {
	      mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	      mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	      mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	    }
        res.SetPointQuadrature(n-this->points_nodal_nd.GetM(),mat);	  
      }
  }
  
  
  //! transformation DFi for dof points in the case of straight Wedge  
  void WedgeGeomReference::DFjElemDofLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int n = this->points_nodal_nd.GetM()+this->points3d.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM()+this->points_dof3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 6; i++)
	  for (int j = 0; j < 3; j++)
	    {
	      mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	      mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	      mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	    }
        res.SetPointDof(n-this->points_nodal_nd.GetM()-this->points3d.GetM(),mat);	  
      }
  }
  
    
  //! transformation Fi for nodal points in the case of curved Wedge
  /*!
    \param[in] s four vertices of the element
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
   */
  void WedgeGeomReference::
  FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
                   const Mesh<Dimension3>& mesh, int nquad, const Volume& elt) const
  {
    ElementGeomReference<Dimension3>::FjElemNodalCurve(s, res, mesh, nquad);
    
    // DISP(s); DISP(nquad); DISP(res);
    R3 res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    // for other nodal points, we use coefFi_curve
    for (int node = 4*this->order_geom*this->order_geom+2;
         node < this->points_nodal_nd.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < 4*this->order_geom*this->order_geom+2; k++)
	  Add(coefFi_curve(node, k), res.GetPointNodal(k), res_n);
	
	res.SetPointNodal(node, res_n);
      }
  }
  
  
  //! transformation Fi for quadrature points in the case of curved Wedge
  void WedgeGeomReference::
  FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
                        const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n; int offset = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }
  
  //! transformation Fi for dof points in the case of curved Wedge
  void WedgeGeomReference::
  FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
                 const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n; int offset = this->points_nodal_nd.GetM()+this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node) ,res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  //! transformation DFi for nodal points in the case of curved Wedge
  void WedgeGeomReference::
  DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                    SetMatrices<Dimension3>& res,
                    const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 tmp;
    R3 grad; grad.Fill(0);
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int i = 0; i < this->points_nodal_nd.GetM(); i++)
      {
	tmp.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  {
	    grad = this->GradientPhi_Nodal(node, i);
	    tmp(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	    tmp(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
            tmp(2,0) += PTReel.GetPointNodal(node)(2)*grad(0);
	    
            tmp(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	    tmp(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
            tmp(2,1) += PTReel.GetPointNodal(node)(2)*grad(1);
            
            tmp(0,2) += PTReel.GetPointNodal(node)(0)*grad(2);
	    tmp(1,2) += PTReel.GetPointNodal(node)(1)*grad(2);
            tmp(2,2) += PTReel.GetPointNodal(node)(2)*grad(2);
	  }
        
	res.SetPointNodal(i,tmp);
      } 
  }
  
  
  //! transformation DFi for quadrature points in the case of curved Wedge  
  void WedgeGeomReference
  ::DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                           SetMatrices<Dimension3>& res,
                           const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n; int offset = this->points_nodal_nd.GetM();
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }

  }
  
  
  //! transformation DFi for dof points in the case of curved Wedge
  void WedgeGeomReference::
  DFjElemDofCurve(const VectR3& s,const SetPoints<Dimension3>& PTReel,
                  SetMatrices<Dimension3>& res,const Mesh<Dimension3>& mesh,int nquad) const
  {
    Matrix3_3 res_n; int offset = this->points_nodal_nd.GetM() + this->points3d.GetM();
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  /*******************
   * Other functions *
   *******************/
  

  //! Evaluating orthogonal functions at point    
  void WedgeGeomReference
  ::ComputeValuesPhiOrthoRef(int r, const Array3D<int>& NumOrtho,
                             const VectReal_wp& InvWeightFct,
                             const R3& point, VectReal_wp & phi) const
  {
    phi.Reallocate(InvWeightFct.GetM());
    VectReal_wp Pn, Pn_tri; R2 pt_tri(point(0), point(1));
    EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*point(2) - 1.0, Pn);
    function_geom_tri.ComputeValuesPhiOrthoRef(r, function_geom_tri.NumOrtho2D,
                                               function_geom_tri.InvWeightPolynomial,
                                               pt_tri, Pn_tri);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r; k++)
	  {
	    int node = NumOrtho(i, j, k);
	    int node2d = function_geom_tri.NumOrtho2D(i, j);
	    phi(node) = Pn_tri(node2d)*Pn(k)*CoefLegendre(k);
	  }
    
  }
  

  //! evaluating gradient of orthogonal functions at point    
  void WedgeGeomReference
  ::ComputeGradientPhiOrthoRef(int r, const Array3D<int>& NumOrtho,
                               const VectReal_wp& InvWeight,
                               const R3& point, VectR3 & dphi) const
  {
    dphi.Reallocate(InvWeight.GetM());
    VectReal_wp Pn, dPn, Pn_tri; VectR2 dPn_tri; R2 pt_tri(point(0), point(1));
    EvaluateJacobiPolynomial(LegendrePolynom, r, 2.0*point(2) - 1.0, Pn, dPn);
    function_geom_tri.ComputeValuesPhiOrthoRef(r, function_geom_tri.NumOrtho2D,
                                               function_geom_tri.InvWeightPolynomial,
                                               pt_tri, Pn_tri);
    
    function_geom_tri.ComputeGradientPhiOrthoRef(r, function_geom_tri.NumOrtho2D,
                                               function_geom_tri.InvWeightPolynomial,
                                                 pt_tri, dPn_tri);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r; k++)
	  {
	    int node = NumOrtho(i, j, k);
	    int node2d = function_geom_tri.NumOrtho2D(i, j);
	    Real_wp vz = Pn(k)*CoefLegendre(k);
	    Real_wp dvz = 2.0*dPn(k)*CoefLegendre(k);
	    dphi(node).Init(dPn_tri(node2d)(0)*vz, dPn_tri(node2d)(1)*vz, Pn_tri(node2d)*dvz);
	  }
  }
  

  void WedgeGeomReference::ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const
  {
    cout << "not used " << endl;
    abort();
  }
  
  
  //! Evaluating shape functions on a point of the element
  /*!
    \param[in] point local point where functions are evaluated
    \param[out] phi values of basis functions on point
  */    
  void WedgeGeomReference::ComputeValuesPhiNodalRef(const R3& point, VectReal_wp & phi) const
  {
    phi.Reallocate(this->points_nodal_nd.GetM());
    int r = this->order_geom;
    VectReal_wp Pz(r+1), phi2d;
    for (int i = 0; i <= r; i++)
      Pz(i) = lob_geom.EvaluatePhi(i, point(2));
    
    R2 pt_tri(point(0), point(1));
    function_geom_tri.ComputeValuesPhiNodalRef(pt_tri, phi2d);   
    
    // final result
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node2d = this->NumNodes2D_tri(i, j);
	  for (int k = 0; k <= r; k++)
	    phi(this->NumNodes3D(i, j, k)) = phi2d(node2d)*Pz(k);
	}
  }
  
    
  //! Evaluating gradient of shape functions on a point of the element
  /*!
    \param[in] point local point where functions are evaluated
    \param[out] dphi gradient of shape functions on point
  */  
  void WedgeGeomReference::ComputeGradientPhiNodalRef(const R3& point, VectR3 & dphi) const
  {
    dphi.Reallocate(this->points_nodal_nd.GetM());
    int r = this->order_geom;
    VectReal_wp Pz(r+1), dPz(r+1), phi2d; VectR2 grad_phi2d;
    for (int i = 0; i <= r; i++)
      {
	Pz(i) = lob_geom.EvaluatePhi(i, point(2));
	dPz(i) = lob_geom.EvaluatePhiGrad(i, point(2));
      }
    
    R2 pt_tri(point(0), point(1));
    function_geom_tri.ComputeValuesPhiNodalRef(pt_tri, phi2d);   
    function_geom_tri.ComputeGradientPhiNodalRef(pt_tri, grad_phi2d);   
    
    // final result
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node2d = this->NumNodes2D_tri(i, j);
	  for (int k = 0; k <= r; k++)
	    dphi(this->NumNodes3D(i, j, k)).Init(grad_phi2d(node2d)(0)*Pz(k),
                                                 grad_phi2d(node2d)(1)*Pz(k),
                                                 phi2d(node2d)*dPz(k));
	}    
  }
  
  
  //! displays informations class WedgeGeomReference
  ostream& operator <<(ostream& out, const WedgeGeomReference& e)
  {
    out << static_cast<const ElementGeomReference<Dimension3>& >(e);
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_GEOM_REFERENCE_CXX
#endif
