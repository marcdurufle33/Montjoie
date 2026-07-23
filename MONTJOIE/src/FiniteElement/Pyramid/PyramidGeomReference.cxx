#ifndef MONTJOIE_FILE_PYRAMID_GEOM_REFERENCE_CXX

namespace Montjoie
{

  //! default constructor
  PyramidGeomReference::PyramidGeomReference() : ElementGeomReference<Dimension3>()
  {
    type_interpolation_nodal = LOBATTO_BASIS;
        
    // normales of the five faces
    this->normale.Reallocate(5);
    this->normale(0).Init(0,0,-4);
    this->normale(1).Init(0,-2,2);
    this->normale(2).Init(2,0,2);
    this->normale(3).Init(0,2,2);
    this->normale(4).Init(-2,0,2);

    this->tangente_loc_x.Reallocate(5);
    this->tangente_loc_y.Reallocate(5);
    this->tangente_loc_x(0).Init(2,0,0);
    this->tangente_loc_y(0).Init(0,2,0);
    this->tangente_loc_x(1).Init(2,0,0);
    this->tangente_loc_y(1).Init(1,1,1);
    this->tangente_loc_x(2).Init(0,2,0);
    this->tangente_loc_y(2).Init(-1,1,1);
    this->tangente_loc_x(3).Init(2,0,0);
    this->tangente_loc_y(3).Init(1,-1,1);
    this->tangente_loc_x(4).Init(0,2,0);
    this->tangente_loc_y(4).Init(1,1,1);
    
    this->hybrid_type_elt = 1;
    this->nb_vertices_elt = 5;
    this->nb_edges_elt = 8;
    this->nb_boundaries_elt = 5;

    this->is_local_face_quad.Reallocate(5);
    this->is_local_face_quad.Fill(false);
    this->is_local_face_quad(0) = true;
  }
  

  size_t PyramidGeomReference::GetMemorySize() const
  {
    size_t taille = ElementGeomReference<Dimension3>::GetMemorySize();
    taille += function_tri.GetMemorySize();
    taille += function_quad.GetMemorySize();
    taille += EdgesNodal.GetMemorySize();
    taille += coefFi_curve.GetMemorySize();
    taille += coefFi.GetMemorySize() + coefDFi_dx.GetMemorySize() + coefDFi_dy.GetMemorySize()
      + coefDFi_dz.GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += Seldon::GetMemorySize(EvenJacobiPolynom);
    taille += InvWeightPolynomial.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    taille += CoefEvenJacobi.GetMemorySize();
    taille += InverseVDM.GetMemorySize();
    taille += NumOrtho3D.GetMemorySize();
    taille += GradientPhiNodal_Quad.GetMemorySize();
    taille += GradientPhiNodal_Dof.GetMemorySize();
    return taille;
  }

    
  //! constructing finite element
  void PyramidGeomReference::
  ConstructFiniteElement(int rgeom)
  {
    this->order_geom = rgeom;
    
    // 2-D basis functions
    if (type_interpolation_nodal == REGULAR_BASIS)
      function_tri.type_interpolation_nodal = function_tri.REGULAR_BASIS;
    else
      function_tri.type_interpolation_nodal = function_tri.LOBATTO_BASIS;
    
    function_tri.ConstructFiniteElement(rgeom);
    
    // orthogonal functions
    ComputeOrthogonalFunctions(rgeom);
    
    // shape functions (used for Fi)
    ConstructNodalShapeFunctions(rgeom);
  }
  
  
  //! construction of shape functions, which are used for definition of Fi for curved elements  
  void PyramidGeomReference::ConstructNodalShapeFunctions(int r)
  {
    // clearing previous points
    this->points_nodal1d.Clear(); this->points_nodal2d_quad.Clear();
    this->points_nodal2d_tri.Clear();
    
    // number of nodal points
    MeshNumbering<Dimension3>::ConstructPyramidalNumbering(r, this->NumNodes3D,
							   this->CoordinateNodes);

    this->points_nodal_nd.Reallocate((r+1)*(r+2)*(2*r+3)/6);
        
    this->EdgesNodal.Reallocate(r+1, 8);
    for (int i = 0; i <= r; i++)
      {
	this->EdgesNodal(i, 0) = this->NumNodes3D(i, 0, 0);
	this->EdgesNodal(i, 1) = this->NumNodes3D(r, i, 0);
	this->EdgesNodal(i, 2) = this->NumNodes3D(i, r, 0);
	this->EdgesNodal(i, 3) = this->NumNodes3D(0, i, 0);
	this->EdgesNodal(i, 4) = this->NumNodes3D(0, 0, i);
	this->EdgesNodal(i, 5) = this->NumNodes3D(r-i, 0, i);
	this->EdgesNodal(i, 6) = this->NumNodes3D(r-i, r-i, i);
	this->EdgesNodal(i, 7) = this->NumNodes3D(0, r-i, i);
      }

    Matrix<int> coor;
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, this->NumNodes2D_quad, coor);
    
    this->FacesNodal.Reallocate((r+1)*(r+1), 5);
    this->FacesNodal.Fill(-1);
    // quadrangular face
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  int node = this->NumNodes2D_quad(i,j);
	  this->FacesNodal(node, 0) = this->NumNodes3D(i, j, 0);
	}
    
    MeshNumbering<Dimension2>::ConstructTriangularNumbering(r, this->NumNodes2D_tri, coor);
    
    // triangular faces
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node = this->NumNodes2D_tri(i,j);
	  this->FacesNodal(node,1) = this->NumNodes3D(i, 0, j);
	  this->FacesNodal(node,2) = this->NumNodes3D(r-j, i, j);
	  this->FacesNodal(node,3) = this->NumNodes3D(i, r-j, j);
	  this->FacesNodal(node,4) = this->NumNodes3D(0, i, j);
	}

    lob_geom.ConstructQuadrature(r, lob_geom.QUADRATURE_LOBATTO);
    
    ConstructLobattoPoints(r, this->points_nodal1d, this->points_nodal2d_tri,
                           this->points_nodal2d_quad, this->NumNodes2D_quad,
                           this->points_nodal_nd);
    
    InvWeightPolynomial.Reallocate(this->points_nodal_nd.GetM());
    InvWeightPolynomial.Fill(0);
    NumOrtho3D.Reallocate(r+1, r+1, r+1); NumOrtho3D.Fill(-1);
    int nb = 0;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        for (int k = 0; k <= r-max(i,j); k++)
          {
            NumOrtho3D(i, j, k) = nb;
            InvWeightPolynomial(NumOrtho3D(i, j, k))
              = CoefLegendre(i)*CoefLegendre(j)*CoefEvenJacobi(max(i,j), k);
            nb++;
          }
        
    ComputeLagrangianFunctions(r);
    
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
  
  
  //! regularly spaced nodal points
  void PyramidGeomReference::
  ConstructRegularPoints(int r, VectReal_wp& points_nodal1d_, VectR2& points_nodal2d_tri_,
                         const Matrix<int>& NumNodes2D_tri_,
                         VectR2& points_nodal2d_quad_, const Matrix<int>& NumNodes2D_quad_,
                         VectR3& points_nodal3d_, const Array3D<int>& NumNodes3D_)
  {    
    points_nodal1d_.Reallocate(r+1);
    for (int i = 0; i <= r; i++)
      points_nodal1d_(i) = Real_wp(i)/r;
    
    points_nodal2d_quad_.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  points_nodal2d_quad_(NumNodes2D_quad_(i, j)).
	    Init(Real_wp(i)/r, Real_wp(j)/r);
	}

    points_nodal2d_tri_.Reallocate((r+1)*(r+2)/2);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  points_nodal2d_tri_(NumNodes2D_tri_(i, j)).
	    Init(Real_wp(i)/r, Real_wp(j)/r);
	}
    
    points_nodal3d_.Reallocate((r+1)*(r+2)*(2*r+3)/6);
    for (int k = 0; k <= r; k++)
      {
	Real_wp z = Real_wp(k)/r;
	if (k == 0)
	  z = 0;
	else if (k == r)
	  z = 1;
	
	for (int i = 0; i <= r-k; i++)
	  for (int j = 0; j <= r-k; j++)
	    {
	      Real_wp x = 0, y = 0;
	      if (k!= r)
		{
		  x = (1-z)*(Real_wp(2)*Real_wp(i)/(r-k) - Real_wp(1));
		  y = (1-z)*(Real_wp(2)*Real_wp(j)/(r-k) - Real_wp(1));
		  points_nodal3d_(NumNodes3D_(i, j, k)).Init(x, y, z);
		}
	      else
		points_nodal3d_(NumNodes3D_(i, j, k)).Init(0, 0, 1);
	    }
      }
  }
  

  //! construcing non-regular interpolation such that restriction 
  //! on quadrilateral base is Gauss-Lobatto points
  void PyramidGeomReference::
  ConstructLobattoPoints(int r, VectReal_wp& points_nodal1d_, VectR2& points_nodal2d_tri_,
                         VectR2& points_nodal2d_quad_, const Matrix<int>& NumNodes2D_quad_,
                         VectR3& points_nodal3d_)
  {
    // we construct 3d nodal points
    points_nodal3d_.Reallocate((r+1)*(r+2)*(2*r+3)/6);
    
    // vertices
    int nb = 0;
    points_nodal3d_(nb) = R3(-1,-1,0);
    GenerateSymPts(points_nodal3d_(nb),
                   points_nodal3d_(nb+1), points_nodal3d_(nb+2), points_nodal3d_(nb+3));
    nb += 4;
    // apex
    points_nodal3d_(nb) = R3(0,0,1);
    nb++;
    
    // we construct lobatto 1d nodal points and triangular nodal points
    TriangleGeomReference::ConstructLobattoPoints(r, 0, points_nodal1d_,
						  points_nodal2d_tri_);
    
    points_nodal2d_quad_.Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	points_nodal2d_quad_(NumNodes2D_quad_(i, j)) = R2(points_nodal1d_(i), points_nodal1d_(j));
    
    if (r > 1)
      {
        // edges on the quadrangle
	for (int i = 1;  i < r; i++)
          {
            points_nodal3d_(nb) = R3(2.0*points_nodal1d_(i)-1, -1, 0);
            points_nodal3d_(nb+r-1) = R3(1.0, 2.0*points_nodal1d_(i)-1, 0);
            points_nodal3d_(nb+2*(r-1)) = R3(2.0*points_nodal1d_(i)-1, 1, 0);
            points_nodal3d_(nb+3*(r-1)) = R3(-1.0, 2.0*points_nodal1d_(i)-1, 0);
	    nb++;
	  }
	nb += 3*(r-1);
	
	// edges on the triangle
	for (int i = 1;  i < r; i++)
          {	  
            points_nodal3d_(nb) = R3(points_nodal1d_(i)-1,
                                     points_nodal1d_(i)-1, points_nodal1d_(i));
	    GenerateSymPts(points_nodal3d_(nb), points_nodal3d_(nb+r-1),
                           points_nodal3d_(nb+2*(r-1)), points_nodal3d_(nb+3*(r-1)));
	    nb++;
	  }
	nb += 3*(r-1);
	
        // quadrangular face
	for (int i = 1;  i < r; i++)
	  for (int j = 1;  j < r; j++)
            points_nodal3d_(nb++) = R3(2*points_nodal1d_(i)-1, 2*points_nodal1d_(j)-1, 0);
      }
    
    if (r > 2)
      {
        // triangular face
        VectR2 points_nodal2d_inside_tri((r-1)*(r-2)/2);
        for (int i = 0; i < points_nodal2d_inside_tri.GetM(); i++)
          points_nodal2d_inside_tri(i) = points_nodal2d_tri_(i+3*r);
	
        int nb_points_nodal_inside_tri = points_nodal2d_inside_tri.GetM();
        for (int i = 0;  i < points_nodal2d_inside_tri.GetM(); i++)
          {
	    points_nodal3d_(nb)
              = R3(2*points_nodal2d_inside_tri(i)(0)+points_nodal2d_inside_tri(i)(1)-1,
                   points_nodal2d_inside_tri(i)(1)-1,points_nodal2d_inside_tri(i)(1));

	    points_nodal3d_(nb+nb_points_nodal_inside_tri)
              = R3(Real_wp(1)-points_nodal2d_inside_tri(i)(1),
                   2*points_nodal2d_inside_tri(i)(0)+points_nodal2d_inside_tri(i)(1)-1,
                   points_nodal2d_inside_tri(i)(1));
            
	    points_nodal3d_(nb+2*nb_points_nodal_inside_tri)
              = R3(2*points_nodal2d_inside_tri(i)(0)+points_nodal2d_inside_tri(i)(1)-1,
                   Real_wp(1) - points_nodal2d_inside_tri(i)(1), points_nodal2d_inside_tri(i)(1));
	    
	    points_nodal3d_(nb+3*nb_points_nodal_inside_tri)
              = R3(points_nodal2d_inside_tri(i)(1)-Real_wp(1),
                   2*points_nodal2d_inside_tri(i)(0)+points_nodal2d_inside_tri(i)(1)-1,
                   points_nodal2d_inside_tri(i)(1));
	    nb++;
	  }
        
	nb += 3*points_nodal2d_inside_tri.GetM();
        
        // interior points
        // we construct regular 1d nodal points
	Globatto<Real_wp> lob_int;
	VectReal_wp points_nodal1d_reg, points_nodal1d_int;
        points_nodal1d_reg.Reallocate(r);
        points_nodal1d_reg.Fill();
        Mlt(Real_wp(1)/(r-1), points_nodal1d_reg);
	
	points_nodal3d_(nb) = R3(0, 0, points_nodal1d_reg(r-2));
	nb++;
	
	if (r > 3)
          {
	    for (int k = 2;  k < r-1; k++)
              {
                // we construct gauss-lobatto 1d nodal points
                lob_int.ConstructQuadrature(k+1, lob_int.QUADRATURE_LOBATTO);
		points_nodal1d_int.Reallocate(k);
		for (int i = 0; i < k; i++)
		  points_nodal1d_int(i) = lob_int.Points(i+1);
		
                Real_wp coeff = points_nodal1d_reg(r-k-1);
                
                for (int i = 0;  i < k; i++)
                  for (int j = 0; j < k; j++)
		    points_nodal3d_(nb++)
                      = R3((2*points_nodal1d_int(i)-1)*(1-coeff),
                           (2*points_nodal1d_int(j)-1)*(1-coeff),coeff);
		
              }
	  }
      }
    // points_nodal3d.WriteText("pts.dat");
  }
  
  
  //! constructing orthogonal functions
  void PyramidGeomReference::ComputeOrthogonalFunctions(int r)
  {
    // generating Jacobi polynomials P_m^{2i,0}
    GetJacobiPolynomial(LegendrePolynom, r+1, Real_wp(0), Real_wp(0));
    EvenJacobiPolynom.Reallocate(r+2);
    for (int i = 0; i <= r+1; i++)
      GetJacobiPolynomial(EvenJacobiPolynom(i), r+1, Real_wp(2*i+2), Real_wp(0));
    
    // we compute weight of polynomials in order to obtain orthonormal expansion
    CoefLegendre.Reallocate(r+2); CoefLegendre.Fill(0);
    CoefEvenJacobi.Reallocate(r+2, r+2); CoefEvenJacobi.Fill(0);
    VectReal_wp xi, omega;
    ComputeGaussLegendre(xi, omega, r+1);
    for (int q = 0; q <= r+1; q++)
      {
	VectReal_wp Pn;
	EvaluateJacobiPolynomial(LegendrePolynom, r+1, 2.0*xi(q) - 1.0, Pn);
	for (int i = 0; i <= r+1; i++)
	  CoefLegendre(i) += 2.0*omega(q)*Pn(i)*Pn(i);
	
	for (int i = 0; i <= r+1; i++)
	  {
	    EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r+1, 2.0*xi(q) - 1.0, Pn);
	    for (int j = 0; j <= r+1; j++)
	      CoefEvenJacobi(i,j) += omega(q)*Pn(j)*Pn(j)*pow(1.0-xi(q), 2*i+2); 
	  }
      }

    // storing 1/sqrt(poids) for orthonormalization
    for (int i = 0; i <= r; i++)
      {
	CoefLegendre(i) = 1.0/sqrt(CoefLegendre(i));
	for (int j = 0; j <= r-i; j++)
	  CoefEvenJacobi(i, j) = 1.0/sqrt(CoefEvenJacobi(i, j));
      }

  }
  
  
  //! computation of basis functions
  void PyramidGeomReference::ComputeLagrangianFunctions(int r)
  {
    // vandermonde computation to force delta_ij on orthogonal polynomials
    InverseVDM.Reallocate(this->points_nodal_nd.GetM(), this->points_nodal_nd.GetM());
    InverseVDM.Fill(0);
        
    // computing and inverting VDM
    VectReal_wp phi;
    for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
      {
	ComputeValuesPhiOrthoRef(r, NumOrtho3D, InvWeightPolynomial, this->points_nodal_nd(k), phi);
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  InverseVDM(j, k) = phi(j);
      }

    GetInverse(InverseVDM);
  }

  
  //! computation of coefFicients to compute quickly Fi and DFi
  void PyramidGeomReference::ComputeCurvedTransformation()
  {    
    coefFi_curve.Reallocate(this->points_nodal_nd.GetM(), this->points_nodal_nd.GetM());
    coefFi_curve.Fill(0);
    Real_wp x, y, z, xs, ys, zeta, coef, coefs;
    Real_wp one(1), half(0.5), quarter(0.25);
    int r = this->order_geom; VectReal_wp phi;
    for (int i = 3*r*r+2; i < this->points_nodal_nd.GetM(); i++)
      {
	x = this->points_nodal_nd(i)(0);
	y = this->points_nodal_nd(i)(1);
	z = this->points_nodal_nd(i)(2);
	
	// linear part
	ComputeValuesPhiFirstOrder(this->points_nodal_nd(i), phi);
	for (int k = 0; k < 5; k++)
	  coefFi_curve(i, k) = phi(k);
	
	// first edge
	zeta = half + half*x/(one - z);
	coef = half*(one - y - z);
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 0)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) -= (one-zeta)*coef;
	coefFi_curve(i, 1) -= zeta*coef;
	
	// second edge
	zeta = half + half*y/(one - z);
	coef = half*(one + x - z);
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 1)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) -= (one-zeta)*coef;
	coefFi_curve(i, 2) -= zeta*coef;

	// third edge
	zeta = half + half*x/(one-z);
	coef = half*(one + y - z);
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 2)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) -= (one-zeta)*coef;
	coefFi_curve(i, 2) -= zeta*coef;
	
	// fourth edge
	zeta = half + half*y/(one-z);
	coef = half*(one - x - z);
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 3)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) -= (one-zeta)*coef;
	coefFi_curve(i, 3) -= zeta*coef;
	
	// fifth edge
	zeta = z;
	coef = quarter*(one-x-z)*(one-y-z)/((one-z)*(one-z));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 4)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// sixth edge
	zeta = z;
	coef = quarter*(one+x-z)*(one-y-z)/((one-z)*(one-z));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 5)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// seventh edge
	zeta = z;
	coef = quarter*(one+x-z)*(one+y-z)/((one-z)*(one-z));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 6)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 2) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// eigth edge
	zeta = z;
	coef = quarter*(one-x-z)*(one+y-z)/((one-z)*(one-z));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 7)) += coef*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// quadrangular base
	xs = half + half*x/(one-z);
	ys = half + half*y/(one-z);
	coef = one-z;
	for (int j = 0; j <= r; j++)
	  for (int k = 0; k <= r; k++)
	    coefFi_curve(i, this->NumNodes3D(j, k, 0))
              += coef*lob_geom.EvaluatePhi(j, xs)*lob_geom.EvaluatePhi(k, ys);

	for (int k = 0; k <= r; k++)
	  {
	    coefFi_curve(i, this->EdgesNodal(k, 0)) -= coef*(one-ys)*lob_geom.EvaluatePhi(k, xs);
	    coefFi_curve(i, this->EdgesNodal(k, 1)) -= coef*xs*lob_geom.EvaluatePhi(k, ys);
	    coefFi_curve(i, this->EdgesNodal(k, 2)) -= coef*ys*lob_geom.EvaluatePhi(k, xs);
	    coefFi_curve(i, this->EdgesNodal(k, 3)) -= coef*(one-xs)*lob_geom.EvaluatePhi(k, ys);
	  }
	
	coefFi_curve(i, 0) += coef*(one-xs)*(one-ys);
	coefFi_curve(i, 1) += coef*xs*(one-ys);
	coefFi_curve(i, 2) += coef*xs*ys;
	coefFi_curve(i, 3) += coef*(one-xs)*ys;

	// first triangular face
	xs = half*(one+x-z);
	ys = z;
	coef = half*(one-y-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < this->points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, this->FacesNodal(j, 1)) += coef*phi(j);
	
	coefFi_curve(i, 0) -= (one-xs-ys)*coef;
	coefFi_curve(i, 1) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 0)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 1) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 5)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 4)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	// second triangular face
	xs = half*(one+y-z);
	ys = z;
	coef = half*(one+x-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < this->points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, this->FacesNodal(j, 2)) += coef*phi(j);
	
	coefFi_curve(i, 1) -= (one-xs-ys)*coef;
	coefFi_curve(i, 2) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 1)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 2) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 6)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 2) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 5)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	// third triangular face
	xs = half*(one+x-z);
	ys = z;
	coef = half*(one+y-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < this->points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, this->FacesNodal(j, 3)) += coef*phi(j);
	
	coefFi_curve(i, 3) -= (one-xs-ys)*coef;
	coefFi_curve(i, 2) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 2)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 2) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 7)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 6)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 2) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	// fourth triangular face
	xs = half*(one+y-z);
	ys = z;
	coef = half*(one-x-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < this->points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, this->FacesNodal(j, 4)) += coef*phi(j);
	
	coefFi_curve(i, 0) -= (one-xs-ys)*coef;
	coefFi_curve(i, 3) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 3)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 3) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 7)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= r; k++)
	  coefFi_curve(i, this->EdgesNodal(k, 4)) -= coef*coefs*lob_geom.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
      }

    const TriangleGeomReference& tri = this->GetTriangularSurfaceFiniteElement();
    this->InitTriangularCurvature(tri);
  }
  
  
  //! computing coefficients for a fast evaluation of Fi for straight elements
  void PyramidGeomReference::ComputeCoefficientTransformation()
  {    

    VectR3 all_points = this->points_nodal_nd;
    all_points.PushBack(this->points3d);
    all_points.PushBack(this->points_dof3d);
    int nb_points_all = all_points.GetM();
    
    // straight pyramid
    coefFi.Reallocate(nb_points_all, 5);
    coefDFi_dx.Reallocate(nb_points_all, 5);
    coefDFi_dy.Reallocate(nb_points_all, 5);
    coefDFi_dz.Reallocate(nb_points_all, 5);
    VectReal_wp phi(5);
    VectR3 dphi(5);

    for (int i = 0; i < nb_points_all; i++) 
      {
        ComputeValuesPhiFirstOrder(all_points(i), phi);   
        ComputeGradientPhiFirstOrder(all_points(i), dphi);  
        for (int j = 0; j < 5; j++)
	  {
	    coefFi(i,j) = phi(j);
	    coefDFi_dx(i,j) = dphi(j)(0);
	    coefDFi_dy(i,j) = dphi(j)(1);
	    coefDFi_dz(i,j) = dphi(j)(2);
	  }
      }

    // for other points, we use lagrangian interpolation
    coefFi_curve.Resize(nb_points_all, this->points_nodal_nd.GetM());
    for (int i = this->points_nodal_nd.GetM(); i < nb_points_all; i++)
      {
	ComputeValuesPhiNodalRef(all_points(i), phi);
	for (int j = 0; j < this->points_nodal_nd.GetM(); j++)
	  coefFi_curve(i, j) = phi(j);
      }

    GradientPhiNodal_Quad.Reallocate(this->points_nodal_nd.GetM(), this->points3d.GetM());
    GradientPhiNodal_Dof.Reallocate(this->points_nodal_nd.GetM(), this->points_dof3d.GetM());
    for (int j = 0; j < this->points3d.GetM(); j++)
      {
        ComputeGradientPhiNodalRef(this->points3d(j), dphi);
        
        for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
          GradientPhiNodal_Quad(k, j) = dphi(k);
      }

    for (int j = 0; j < this->points_dof3d.GetM(); j++)
      {
        ComputeGradientPhiNodalRef(this->points_dof3d(j), dphi);
        
        for (int k = 0; k < this->points_nodal_nd.GetM(); k++)
          GradientPhiNodal_Dof(k, j) = dphi(k);
      }
  }


  void PyramidGeomReference
  ::GetGradient3D_FromGradient2D(int num_loc, int i, const VectReal_wp& nabla_nx,
				 const VectReal_wp& nabla_ny, const VectReal_wp& nabla_nz,
				 R3& d_nx, R3& d_ny, R3& d_nz) const
  {
    switch (num_loc)
      {
      case 0:
	d_nx(0) = 0.5*nabla_nx(2*i); d_nx(1) = 0.5*nabla_nx(2*i+1); d_nx(2) = 0;
	d_ny(0) = 0.5*nabla_ny(2*i); d_ny(1) = 0.5*nabla_ny(2*i+1); d_ny(2) = 0;
	d_nz(0) = 0.5*nabla_nz(2*i); d_nz(1) = 0.5*nabla_nz(2*i+1); d_nz(2) = 0;
	break;
      case 1:
	// DF = [2 1; 0 1; 0 1] and DF*-1 = [1/2 0; -1/4 1/2; -1/4 1/2]
	d_nx(0) = 0.5*nabla_nx(2*i);
	d_ny(0) = 0.5*nabla_ny(2*i);
	d_nz(0) = 0.5*nabla_nz(2*i);
        
	d_nx(1) = -0.25*nabla_nx(2*i) + 0.5*nabla_nx(2*i+1);
	d_ny(1) = -0.25*nabla_ny(2*i) + 0.5*nabla_ny(2*i+1);
	d_nz(1) = -0.25*nabla_nz(2*i) + 0.5*nabla_nz(2*i+1);
        
	d_nx(2) = d_nx(1);
	d_ny(2) = d_ny(1);
	d_nz(2) = d_nz(1);
	break;
      case 2 :
	// DF = [0 -1; 2 1; 0 1] and DF*-1 = [1/4 -1/2; 1/2 0; -1/4 1/2]
	d_nx(0) = 0.25*nabla_nx(2*i) - 0.5*nabla_nx(2*i+1);
	d_ny(0) = 0.25*nabla_ny(2*i) - 0.5*nabla_ny(2*i+1);
	d_nz(0) = 0.25*nabla_nz(2*i) - 0.5*nabla_nz(2*i+1);
        
	d_nx(1) = 0.5*nabla_nx(2*i);
	d_ny(1) = 0.5*nabla_ny(2*i);
	d_nz(1) = 0.5*nabla_nz(2*i);
        
	d_nx(2) = -d_nx(0);
	d_ny(2) = -d_ny(0);
	d_nz(2) = -d_nz(0);
	break;
      case 3 :
	// DF = [2 1; 0 -1; 0 1] and DF*-1 = [1/2 0; 1/4 -1/2; -1/4 1/2]
	d_nx(0) = 0.5*nabla_nx(2*i);
	d_ny(0) = 0.5*nabla_ny(2*i);
	d_nz(0) = 0.5*nabla_nz(2*i);
        
	d_nx(1) = 0.25*nabla_nx(2*i) - 0.5*nabla_nx(2*i+1);
	d_ny(1) = 0.25*nabla_ny(2*i) - 0.5*nabla_ny(2*i+1);
	d_nz(1) = 0.25*nabla_nz(2*i) - 0.5*nabla_nz(2*i+1);
        
	d_nx(2) = -d_nx(1);
	d_ny(2) = -d_ny(1);
	d_nz(2) = -d_nz(1);
	break;
      case 4 :
	// DF = [0 1; 2 1; 0 1] and DF*-1 = [-1/4 1/2; 1/2 0; -1/4 1/2]
	d_nx(0) = -0.25*nabla_nx(2*i) + 0.5*nabla_nx(2*i+1);
	d_ny(0) = -0.25*nabla_ny(2*i) + 0.5*nabla_ny(2*i+1);
	d_nz(0) = -0.25*nabla_nz(2*i) + 0.5*nabla_nz(2*i+1);
        
	d_nx(1) = 0.5*nabla_nx(2*i);
	d_ny(1) = 0.5*nabla_ny(2*i);
	d_nz(1) = 0.5*nabla_nz(2*i);
        
	d_nx(2) = d_nx(0);
	d_ny(2) = d_ny(0);
	d_nz(2) = d_nz(0);
      }
  }
    

  //! returns a tangential vector on a face of the pyramid from a 2-D vector on the square/triangle
  R3 PyramidGeomReference::GetTangentialVector(int num_loc, const R2& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
	return R3(0.5*vec_u(0), 0.5*vec_u(1), Real_wp(0));
      case 1:
	// DF = [2 1; 0 1; 0 1] and DF*-1 = [1/2 0; -1/4 1/2; -1/4 1/2]
	{
	  Real_wp coef = -0.25*vec_u(0)+0.5*vec_u(1);
	  return R3(0.5*vec_u(0), coef, coef);
	}
      case 2 :
	// DF = [0 -1; 2 1; 0 1] and DF*-1 = [1/4 -1/2; 1/2 0; -1/4 1/2]
	{
	  Real_wp coef = 0.25*vec_u(0) - 0.5*vec_u(1);
	  return R3(coef, 0.5*vec_u(0), -coef);
	}
	break;
      case 3 :
	// DF = [2 1; 0 -1; 0 1] and DF*-1 = [1/2 0; 1/4 -1/2; -1/4 1/2]
	{
	  Real_wp coef = 0.25*vec_u(0) - 0.5*vec_u(1);
	  return R3(0.5*vec_u(0), coef, -coef);
	}
	break;
      case 4 :
	// DF = [0 1; 2 1; 0 1] and DF*-1 = [-1/4 1/2; 1/2 0; -1/4 1/2]
	{
	  Real_wp coef = -0.25*vec_u(0) + 0.5*vec_u(1);
	  return R3(coef, 0.5*vec_u(0), coef);
	}
      }

    return R3();
  }
  

  //! transpose operation of GetTangentialVector
  R2 PyramidGeomReference::TransposeTangentialVector(int num_loc, const R3& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
	return R2(0.5*vec_u(0), 0.5*vec_u(1));
      case 1:
        return R2(0.5*vec_u(0) -0.25*(vec_u(1)+vec_u(2)),
                  0.5*(vec_u(1) + vec_u(2)));
      case 2 :
        return R2(0.25*(vec_u(0) - vec_u(2)) + 0.5*vec_u(1),
                  -0.5*(vec_u(0) - vec_u(2)));
      case 3 :
        return R2(0.5*vec_u(0) - 0.25*(vec_u(1)-vec_u(2)),
                  -0.5*(vec_u(1)-vec_u(2)));
      case 4 :
        return R2(-0.25*(vec_u(0) + vec_u(2)) + 0.5*vec_u(1),
                  0.5*(vec_u(0) + vec_u(2)));
      }
    
    return R2();
  }


  //! transpose operation of GetTangentialVector
  R2_Complex_wp PyramidGeomReference::TransposeTangentialVector(int num_loc, const R3_Complex_wp& vec_u) const
  {
    switch (num_loc)
      {
      case 0:
	return R2_Complex_wp(0.5*vec_u(0), 0.5*vec_u(1));
      case 1:
        return R2_Complex_wp(0.5*vec_u(0) -0.25*(vec_u(1)+vec_u(2)),
                  0.5*(vec_u(1) + vec_u(2)));
      case 2 :
        return R2_Complex_wp(0.25*(vec_u(0) - vec_u(2)) + 0.5*vec_u(1),
                  -0.5*(vec_u(0) - vec_u(2)));
      case 3 :
        return R2_Complex_wp(0.5*vec_u(0) - 0.25*(vec_u(1)-vec_u(2)),
                  -0.5*(vec_u(1)-vec_u(2)));
      case 4 :
        return R2_Complex_wp(-0.25*(vec_u(0) + vec_u(2)) + 0.5*vec_u(1),
                  0.5*(vec_u(0) + vec_u(2)));
      }
    
    return R2_Complex_wp();
  }
  
  
  /****************
   * Fj transform *
   ****************/


  //! transformation Fi in the case of straight pyramid   
  void PyramidGeomReference::FjLinear(const VectR3& s, const R3& point, R3& res) const
  {
    // Fi = phi1*S0 + phi2*S1 + phi3*S2 + phi4*S3 + z*S4
    res.Zero();
    VectReal_wp phi; ComputeValuesPhiFirstOrder(point, phi);

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 3; j++)
        res(j) += phi(i)*s(i)(j);	
  }


  //! transformation DFi in the case of straight pyramid   
  void PyramidGeomReference::DFjLinear(const VectR3& s,const R3& point,Matrix3_3& res) const
  {
    res.Zero();
    VectR3 dphi; ComputeGradientPhiFirstOrder(point,dphi);

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 3; j++)
	{
	  res(j,0) += dphi(i)(0)*s(i)(j);
	  res(j,1) += dphi(i)(1)*s(i)(j);
	  res(j,2) += dphi(i)(2)*s(i)(j);
	}
  }


  //! returns the minimal length of the element
  Real_wp PyramidGeomReference
  ::GetMinimalSize(const VectR3& s) const
  {
    Real_wp h = Distance(s(0), s(1));
    h = min(h, Distance(s(0), s(3)));
    h = min(h, Distance(s(0), s(4)));
    h = min(h, Distance(s(1), s(2)));
    h = min(h, Distance(s(1), s(4)));
    h = min(h, Distance(s(2), s(3)));
    h = min(h, Distance(s(2), s(4)));
    h = min(h, Distance(s(3), s(4)));
    return h;
  }
  

  //! returns true if the point is outside the unit pyramid
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
  */  
 bool PyramidGeomReference
 ::OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const
  {
    Real_wp Rmax = max(abs(s(0)(0)), abs(s(0)(1))); Rmax = max(Rmax, abs(s(0)(2)));
    Rmax = max(Rmax, abs(s(1)(0))); Rmax = max(Rmax, abs(s(1)(1))); Rmax = max(Rmax, abs(s(1)(2)));
    Rmax = max(Rmax, abs(s(2)(0))); Rmax = max(Rmax, abs(s(2)(1))); Rmax = max(Rmax, abs(s(2)(2)));
    Rmax = max(Rmax, abs(s(3)(0))); Rmax = max(Rmax, abs(s(3)(1))); Rmax = max(Rmax, abs(s(3)(2)));
    Rmax = max(Rmax, abs(s(4)(0))); Rmax = max(Rmax, abs(s(4)(1))); Rmax = max(Rmax, abs(s(4)(2)));
    
    Real_wp dx = s(0).Distance(s(1)), dy = s(0).Distance(s(3)), dz = s(0).Distance(s(4));
    if ((Rmax == Real_wp(0)) || (dx == Real_wp(0)) || (dy == Real_wp(0)) || (dz == Real_wp(0)))
      return true;
    
    Real_wp one(1);
    if (dz*Xn(2) < -epsilon*Rmax)
      return true;

    if (min(dx, dz)*(Xn(2) - Xn(0) - one) > epsilon*Rmax)
      return true;

    if (min(dy, dz)*(Xn(2) - Xn(1) - one) > epsilon*Rmax)
      return true;

    if (min(dx, dz)*(Xn(2) + Xn(0) - one) > epsilon*Rmax)
      return true;

    if (min(dy, dz)*(Xn(2) + Xn(1) - one) > epsilon*Rmax)
      return true;
    
    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
   */  
  Real_wp PyramidGeomReference::GetDistanceToBoundary(const R3& pointloc) const
  {
    Real_wp distance = pointloc(2); 
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)+pointloc(0)));
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)-pointloc(0))); 
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)+pointloc(1)));
    distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)-pointloc(1)));
    return distance;
  }
  
  
  Real_wp PyramidGeomReference::GetDistanceToBoundary(const R3& pointloc, int n) const
  {
    cout << "Not implemented" << endl;
    abort();
    return Real_wp(0);
  }


  //! project the local point on the boundary, if outside the unit pyramid  
  int PyramidGeomReference::ProjectPointOnBoundary(R3& pointloc) const
  {
    pointloc(2) = max( Real_wp(0), pointloc(2) );
    Real_wp one(1);
    if (pointloc(2)-pointloc(0)<one)
      pointloc(0) = Real_wp(pointloc(2)-one);
    else if (pointloc(2)+pointloc(0)<one)
      pointloc(0) = Real_wp(one-pointloc(2));
    if (pointloc(2)-pointloc(1)<one)
      pointloc(1) = Real_wp(pointloc(2)-one);
    else if (pointloc(2)+pointloc(1)<one)
      pointloc(1) = Real_wp(one-pointloc(2));
    
    return -1;
  }
  
  
  //! evaluates lowest order shape functions at a point in the element
  void PyramidGeomReference
  ::ComputeValuesPhiFirstOrder(const R3& point, VectReal_wp & phi) const
  {
    Real_wp x = point(0), y = point(1), z = point(2);
    phi.Reallocate(5);
    Real_wp one(1);
    if (abs(z-Real_wp(1)) > epsilon_machine)
      {
        phi(0) = 0.25*(one-x-z)*(one-y-z)/(one-z);
        phi(1) = 0.25*(one+x-z)*(one-y-z)/(one-z);
        phi(2) = 0.25*(one+x-z)*(one+y-z)/(one-z);
        phi(3) = 0.25*(one-x-z)*(one+y-z)/(one-z);
        phi(4) = z;
      }
    else
      {
        phi.Fill(Real_wp(0));
	phi(4) = one;
      }
  }
  

  //! evaluates gradient of lowest order shape functions at a point in the element    
  void PyramidGeomReference
  ::ComputeGradientPhiFirstOrder(const R3& point, VectR3& dphi) const
  {
    Real_wp x = point(0), y = point(1), z = point(2);
    dphi.Reallocate(5);
    Real_wp one(1);
    if (abs(z-one) > epsilon_machine)
      {
        dphi(0)(0) = -0.25*(one-z-y)/(one-z); dphi(0)(1) = -0.25*(one-z-x)/(one-z);
        dphi(0)(2) = -0.25*(one-2*z+z*z-x*y)/((one-z)*(one-z));
        dphi(1)(0) = 0.25*(one-z-y)/(one-z); dphi(1)(1) = -0.25*(one-z+x)/(one-z);
        dphi(1)(2) = -0.25*(one-2*z+z*z+x*y)/((one-z)*(one-z));
        dphi(2)(0) = 0.25*(one-z+y)/(one-z); dphi(2)(1) = 0.25*(one-z+x)/(one-z);
        dphi(2)(2) = -0.25*(one-2*z+z*z-x*y)/((one-z)*(one-z));
        dphi(3)(0) = -0.25*(one-z+y)/(one-z); dphi(3)(1) = 0.25*(one-z-x)/(one-z);
        dphi(3)(2) = -0.25*(one-2*z+z*z+x*y)/((one-z)*(one-z));
        dphi(4)(0) = Real_wp(0); dphi(4)(1) = Real_wp(0); dphi(4)(2) = Real_wp(1);
      }
    else
      {
	// we take x = y = 0
	dphi(0)( 0) = -0.25; dphi(0)( 1) = -0.25; dphi(0)( 2) = -0.25;
	dphi(1)( 0) = 0.25; dphi(1)( 1) = -0.25; dphi(1)( 2) = -0.25;
	dphi(2)( 0) = 0.25; dphi(2)( 1) = 0.25; dphi(2)( 2) = -0.25;
	dphi(3)( 0) = -0.25; dphi(3)( 1) = 0.25; dphi(3)( 2) = -0.25;
	dphi(4)( 0) = -0.0; dphi(4)( 1) = -0.0; dphi(4)( 2) = 1.0;
      }
  }
  
  
  //! returns 3-D coordinates of a point on a face
  /*!
    \param[in] num_loc local face number in the quad
    \param[in] point_loc coordinate on the face
    \param[out] res local coordinates of the point in the unit square
  */
  void PyramidGeomReference
  ::GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const
  {
    switch(num_loc)
      {
      case 0:
	res.Init(2*point_loc(0)-1.0, 2*point_loc(1)-1.0, 0);break;
      case 1:
	res.Init(2*point_loc(0)+point_loc(1)-Real_wp(1),point_loc(1)-Real_wp(1),point_loc(1));
        break;
      case 2:
	res.Init(Real_wp(1)-point_loc(1),2*point_loc(0)+point_loc(1)-Real_wp(1),point_loc(1));
        break;
      case 3:
	res.Init(2*point_loc(0)+point_loc(1)-Real_wp(1),Real_wp(1)-point_loc(1),point_loc(1));
        break;
      case 4:
	res.Init(point_loc(1)-Real_wp(1),2*point_loc(0)+point_loc(1)-Real_wp(1),point_loc(1));
        break;
      }
  }


  //! computes 2-D coordinates of a point on the symmetric pyramid
  /*!
    \param[in] num_loc local face number in the quad
    \param[out] point_loc coordinate on the face
    \param[in] pt local coordinates of the point in the unit square
  */
  void PyramidGeomReference
  ::GetLocalCoordOnBoundary(int num_loc, const R3& pt, R2& point_loc) const
  {
    switch(num_loc)
      {
      case 0:
	point_loc.Init(0.5*pt(0)+0.5, 0.5*pt(1)+0.5); break;
      case 1:
	point_loc.Init(0.5*(pt(0)-pt(1)), pt(2)); break;
      case 2:
	point_loc.Init(0.5*(pt(1)+pt(0)), pt(2)); break;
      case 3:
	point_loc.Init(0.5*(pt(0)+pt(1)), pt(2)); break;
      case 4:
	point_loc.Init(0.5*(pt(1)-pt(0)), pt(2)); break;
      }
  }
  
  
  //! transformation Fi in the case of curved pyramid
  void PyramidGeomReference::
  FjCurve(const SetPoints<Dimension3>& PTReel,const R3& pointloc,R3& res) const
  {
    res.Fill(0); VectReal_wp phi;
    ComputeValuesPhiNodalRef(pointloc, phi);
    for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
      Add(phi(node), PTReel.GetPointNodal(node), res);
  }
  
  
  //! transformation DFi in the case of curved pyramid
  void PyramidGeomReference::
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
  void PyramidGeomReference
  ::ComputeCoefJacobian(const VectR3& s, VectReal_wp& CoefJacobian) const
  {
    CoefJacobian.Reallocate(4);
    R3 A1 = -s(0)+s(1)+s(2)-s(3); A1 *= 0.25;
    R3 A2 = -s(0)-s(1)+s(2)+s(3); A2 *= 0.25;
    R3 A3 = Real_wp(4)*s(4)-s(0)-s(1)-s(2)-s(3); A3 *= 0.25;
    R3 C = s(0)-s(1)+s(2)-s(3); C *= 0.25;
    
    CoefJacobian(0) = Determinant(A1,A2,A3);
    CoefJacobian(1) = Determinant(A1,C,A3);
    CoefJacobian(2) = Determinant(C,A2,A3);
    CoefJacobian(3) = Determinant(A1,A2,C);
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
  void PyramidGeomReference::FjElem(const VectR3& s, SetPoints<Dimension3>& res,
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
  void PyramidGeomReference::FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
					 const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad, mesh.Element(nquad));
    else
      FjElemNodalLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all quadrature points  
  void PyramidGeomReference
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
  void PyramidGeomReference::FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
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
      FjElemDofLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all reference points
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */  
  void PyramidGeomReference
  ::DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
            SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
            int nquad) const
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
  void PyramidGeomReference
  ::DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                 SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s,PTReel,res,mesh,nquad);
    else
      DFjElemNodalLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points  
  void PyramidGeomReference
  ::DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                      SetMatrices<Dimension3>& res,
                      const Mesh<Dimension3>& mesh,int nquad) const
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
      DFjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all dof points  
  void PyramidGeomReference
  ::DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
               SetMatrices<Dimension3>& res,
               const Mesh<Dimension3>& mesh, int nquad) const
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
  
  
  //! nodal points in the case of straight pyramid  
  void PyramidGeomReference
  ::FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int n = 0; n < this->points_nodal_nd.GetM(); n++)
      {
        res_n.Fill(0);
        for (int i = 0; i < 5; i++)
          for (int j = 0; j < 3; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointNodal(n,res_n);
      }
  }
  
  
  //! quadrature points in the case of straight pyramid  
  void PyramidGeomReference
  ::FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {  
    R3 res_n;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int n = this->points_nodal_nd.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM(); n++)
      {
        res_n.Fill(0);
        for (int i = 0; i < 5; i++)
          for (int j = 0; j< 3 ; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointQuadrature(n-this->points_nodal_nd.GetM(),res_n);
      }
  }
  
  
  //! dof points in the case of straight pyramid  
  void PyramidGeomReference
  ::FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    //DISP(points_dof3d.GetM());
    for (int n = this->points_nodal_nd.GetM()+this->points3d.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM()+this->points_dof3d.GetM(); n++)
      {
        res_n.Fill(0);
	for (int i = 0; i < 5; i++)
          for (int j = 0; j < 3; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointDof(n-this->points_nodal_nd.GetM()-this->points3d.GetM(),res_n);
      }

  }
  
  
  //! transformation DFi for nodal points in the case of straight pyramid  
  void PyramidGeomReference
  ::DFjElemNodalLinear(const VectR3& s,SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    for (int n = 0; n < this->points_nodal_nd.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 5; i++)
         for (int j = 0; j < 3; j++)
          {
	    mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	    mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	    mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	  }
        res.SetPointNodal(n,mat);	  
      }
  }
  
  
  //! transformation DFi for quadrature points in the case of straight pyramid
  void PyramidGeomReference
  ::DFjElemQuadratureLinear(const VectR3& s,SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int n = this->points_nodal_nd.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 5; i++)
         for (int j = 0; j < 3; j++)
          {
	    mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	    mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	    mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	  }
        res.SetPointQuadrature(n-this->points_nodal_nd.GetM(),mat);	  
      }
  }
  
  
  //! transformation DFi for dof points in the case of straight pyramid  
  void PyramidGeomReference
  ::DFjElemDofLinear(const VectR3& s,SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int n = this->points_nodal_nd.GetM()+this->points3d.GetM();
         n < this->points_nodal_nd.GetM()+this->points3d.GetM()+this->points_dof3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 5; i++)
         for (int j = 0; j < 3; j++)
          {
	    mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	    mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	    mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	  }
        res.SetPointDof(n-this->points_nodal_nd.GetM()-this->points3d.GetM(),mat);	  
      }
  }
  
    
  //! transformation Fi for nodal points in the case of curved pyramid
  /*!
    \param[in] s four vertices of the element
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
   */  
  void PyramidGeomReference
  ::FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
                     const Mesh<Dimension3>& mesh, int nquad, const Volume& elt) const
  {
    ElementGeomReference<Dimension3>::FjElemNodalCurve(s, res, mesh, nquad);
    
    R3 res_n;
    res.ReallocatePointsNodal(this->points_nodal_nd.GetM());
    // for other nodal points, we use coefFi_curve
    int nb_inside = 3*this->order_geom*this->order_geom+2;
    for (int node = nb_inside; node < this->points_nodal_nd.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < nb_inside; k++)
	  Add(coefFi_curve(node, k), res.GetPointNodal(k), res_n);
	
	res.SetPointNodal(node, res_n);
      }
  }
  
  
  //! transformation Fi for quadrature points in the case of curved pyramid  
  void PyramidGeomReference
  ::FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
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
  
  
  //! transformation Fi for dof points in the case of curved pyramid  
  void PyramidGeomReference
  ::FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
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
  
  
  //! transformation DFi for nodal points in the case of curved pyramid  
  void PyramidGeomReference
  ::DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
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
  
  
  //! transformation DFi for quadrature points in the case of curved pyramid  
  void PyramidGeomReference
  ::DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                           SetMatrices<Dimension3>& res,
                           const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n;
    res.ReallocatePointsQuadrature(this->points3d.GetM());
    for (int i = 0; i < this->points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Rank1Update(1.0, PTReel.GetPointNodal(node), GradientPhiNodal_Quad(node, i), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }
  
  
  //! transformation DFi for dof points in the case of curved pyramid  
  void PyramidGeomReference
  ::DFjElemDofCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                    SetMatrices<Dimension3>& res,
                    const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n;
    res.ReallocatePointsDof(this->points_dof3d.GetM());
    for (int i = 0; i < this->points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < this->points_nodal_nd.GetM(); node++)
	  Rank1Update(1.0, PTReel.GetPointNodal(node), GradientPhiNodal_Dof(node, i), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  /*******************
   * Other functions *
   *******************/
  
  
  //! Evaluating orthogonal functions at pointloc
  void PyramidGeomReference
  ::ComputeValuesPhiOrthoRef(int r, const Array3D<int>& NumOrtho,
                             const VectReal_wp& InvWeightFct,
                             const R3& pointloc, VectReal_wp& phi) const
  {
    // orthogonal polynomials :
    // L_i(x/(1-z)) L_j(y/(1-z)) P_k^{2 max(i,j) + 2}(2z-1) (1-z)^(max(i,j))
    // = L_i(a) L_j(b) P_k^(2i+2j+2)(c) ((1-c)/2)^(max(i,j))
    // where (a, b, c) is on the cube [-1,1]^3 and (x,y,z) on the unit tetrahedron
    
    phi.Reallocate(InvWeightFct.GetM());
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0);
    if (abs(1.0 - z) > epsilon_machine)
      {
        a = x/(1.0-z);
        b = y/(1.0-z);
      }
    
    /* for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  int m = max(i, j);
	  Real_wp vloc = 0;
          if (m == 0)
            vloc = 1.0;
          else if (m == 1)
            {
              if (i == 0)
                vloc = y;
              else if (j == 0)
                vloc = x;
              else
                vloc = a*y;
            }
          else
	    vloc = pow(a, i)*pow(b, j)*pow(1.0-z, m);
	  
	  for (int k = 0; k <= r-max(i, j); k++)
	    phi(NumOrtho(i, j, k)) = vloc*pow(z, k);
	}
    
    return;
    */
    VectReal_wp Px, Py, pow_oneMinusZ(r+1);
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, r, b, Py);
    Vector<VectReal_wp> Pz(r+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, 2.0*z-1, Pz(i));
    
    for (int i = 0; i < r; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          int m = max(i, j);
          if (m == 0)
            val = 1.0;
          else if (m == 1)
            {
              if (i == 0)
                val = y;
              else if (j == 0)
                val = x;
              else
                val = a*y;
            }
          else
            val = Px(i)*Py(j)*pow_oneMinusZ(m);
          
          for (int k = 0; k <= r-max(i, j); k++)
            {
              int node = NumOrtho(i, j, k);
              phi(node) = val*Pz(m)(k);
            }
        }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeightFct(i);
  }
  
  

  //! evaluating gradient of orthogonal functions at pointloc
  void PyramidGeomReference
  ::ComputeGradientPhiOrthoRef(int r, const Array3D<int>& NumOrtho,
                               const VectReal_wp& InvWeightFct,
                               const R3& pointloc, VectR3& grad_phi) const
  {
    // orthogonal polynomials :
    // L_i(x/(1-z)) L_j(y/(1-z)) P_k^{2 max(i,j) + 2}(2z-1) (1-z)^(max(i,j))
    // = L_i(a) L_j(b) P_k^(2i+2j+2)(c) ((1-c)/2)^(max(i,j))
    // where (a, b, c) is on the cube [-1,1]^3 and (x,y,z) on the unit tetrahedron
    
    // derivative with respect to a, b, c is equal to :
    // L'_i(a) L_j(b) P_k^(2i+2j+2)(c) ((1-c)/2)^(max(i,j))
    // L_i(a) L'_j(b) P_k^(2i+2j+2)(c) ((1-c)/2)^(max(i,j))
    // L_i(a) L_j(b) ((1-c)/2)^(max(i,j)-1) ( (1-c)/2 P'_k^(2i+2j+2)(c)
    //   - max(i,j)/2 P_k^(2i+2j+2)(c) )
    
    grad_phi.Reallocate(InvWeightFct.GetM());
    
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0), da_dx(0), da_dz(0), db_dy(0), db_dz(0), one(1);
    if (abs(1-z) >  1e3*epsilon_machine)
      {
	a = x/(one-z);
	b = y/(one-z);
	da_dx = one/(one-z); da_dz = x/square(one-z);
	db_dy = one/(one-z); db_dz = y/square(one-z);
      }
    
    VectReal_wp Px, Py, dPx, dPy, pow_oneMinusZ(r+1);
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, r, b, Py, dPy);
    Vector<VectReal_wp> Pz(r+1), dPz(r+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= r; i++)      
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, 2.0*z-1, Pz(i), dPz(i));
    
    for (int i = 0; i < r; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val, dval_dx, dval_dz, dval_dy, valz;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          int m = max(i, j);
          if (m == 0)
            {
              val = 1.0;
              dval_dx = 0.0; dval_dy = 0; dval_dz = 0;
            }
          else if (m == 1)
            {
              if (i == 0)
                {
                  val = y;
                  dval_dx = 0.0; dval_dy = 1.0; dval_dz = 0;
                }
              else if (j == 0)
                {
                  val = x;
                  dval_dx = 1.0; dval_dy = 0; dval_dz = 0;
                }
              else
                {
                  // pathologic case
                  // gradient is equal to (y/(1-z), x/(1-z), xy/(1-z)^2)
                  // which takes different values on the faces of the pyramid
                  val = a*y;
                  dval_dx = b; dval_dy = a; dval_dz = a*b;
                }
            }
          else
            {
              val = Px(i)*Py(j)*pow_oneMinusZ(m);
              dval_dx = dPx(i)*Py(j)*da_dx*pow_oneMinusZ(m);
              dval_dy = dPy(j)*Px(i)*db_dy*pow_oneMinusZ(m);
              dval_dz = dPx(i)*Py(j)*da_dz*pow_oneMinusZ(m)
                + Px(i)*dPy(j)*db_dz*pow_oneMinusZ(m) - m*Px(i)*Py(j)*pow_oneMinusZ(m-1); 
            }
          
          for (int k = 0; k <= r-max(i, j); k++)
            {
              int node = NumOrtho(i, j, k);
              valz = Pz(m)(k);
              grad_phi(node).Init(dval_dx*valz, dval_dy*valz, dval_dz*valz + 2.0*val*dPz(m)(k));
            }
        }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < grad_phi.GetM(); i++)
      grad_phi(i) *= InvWeightFct(i);
    
  }
  
  
  void PyramidGeomReference::ComputeValuesNodalPhi1D(const Real_wp&, VectReal_wp&) const
  {
    cout << "not used " << endl;
    abort();
  }

  
  //! Evaluating shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */
  void PyramidGeomReference
  ::ComputeValuesPhiNodalRef(const R3 &pointloc, VectReal_wp& phi) const
  {
    VectReal_wp psi; phi.Reallocate(this->points_nodal_nd.GetM()); phi.Fill(0);
    ComputeValuesPhiOrthoRef(this->order_geom, NumOrtho3D,
                             InvWeightPolynomial, pointloc, psi);
    Mlt(InverseVDM, psi, phi);
  }
  
  
  //! Evaluating gradient of shape functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of shape functions on pointloc
  */  
  void PyramidGeomReference
  ::ComputeGradientPhiNodalRef(const R3& pointloc, VectR3& grad_phi) const
  {
    VectReal_wp psi, phi; VectR3 grad_psi; 
    psi.Reallocate(this->points_nodal_nd.GetM()); psi.Fill(0);
    phi.Reallocate(this->points_nodal_nd.GetM()); phi.Fill(0);
    grad_phi.Reallocate(this->points_nodal_nd.GetM());
    ComputeGradientPhiOrthoRef(this->order_geom, NumOrtho3D,
                               InvWeightPolynomial, pointloc, grad_psi);
    for (int i = 0; i < psi.GetM(); i++)
      psi(i) = grad_psi(i)(0);
	
    Mlt(InverseVDM, psi, phi);
    
    for (int i = 0; i < psi.GetM(); i++)
      {
	grad_phi(i)(0) = phi(i);
	psi(i) = grad_psi(i)(1);
      }
    
    Mlt(InverseVDM, psi, phi);
    for (int i = 0; i < psi.GetM(); i++)
      {
	grad_phi(i)(1) = phi(i);
	psi(i) = grad_psi(i)(2);
      }
    
    Mlt(InverseVDM, psi, phi);
    for (int i = 0; i < psi.GetM(); i++)
      grad_phi(i)(2) = phi(i);
  }
  
  
  //! displays informations class PyramidGeomReference  
  ostream& operator <<(ostream& out,const PyramidGeomReference& e)
  {
    out << static_cast<const ElementGeomReference<Dimension3>& >(e);
    out<<"Type interpolation "<<e.type_interpolation_nodal<<endl;
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_GEOM_REFERENCE_CXX
#endif
