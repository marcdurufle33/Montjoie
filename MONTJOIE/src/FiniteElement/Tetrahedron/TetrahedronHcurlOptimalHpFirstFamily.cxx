#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronHcurlOptimalHpFirstFamily::TetrahedronHcurlOptimalHpFirstFamily()
    : TetrahedronReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! how to number mesh
  void TetrahedronHcurlOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofTetrahedron(order, nb_dof_loc-nb_dof_boundaries);
    nmap.SetNbDofHexahedron(order, 0);
    
    nmap.SetEqualEdgesDofSymmetry(order, order+1);
    
    // FindHcurlSignEdge(*this, nmap);
    
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_tri_surf);
  }


  size_t TetrahedronHcurlOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<2>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += NumDofsX.GetMemorySize() + NumDofsY.GetMemorySize() + NumDofsZ.GetMemorySize();
    return taille;
  }
  
  
  //! construction of finite element
  void TetrahedronHcurlOptimalHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHcurlElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(false);
  }
  
  
  //! construction of basis functions
  void TetrahedronHcurlOptimalHpFirstFamily::ConstructFunctions()
  {
    int nb_dof_inside_tri = order*(order-1);
    int nb_dof_inside_elt = order*(order-1)*(order-2)/2;
	  
    nb_dof_boundaries = 4*nb_dof_inside_tri + 6*order;
    nb_dof_loc = nb_dof_inside_elt + nb_dof_boundaries;
    nb_dof_tri = 3*order + nb_dof_inside_tri;
    nb_dof_quad = 0;
	  
    /*points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points_dof3d(i) = points3d(i);
    */
    
    VectReal_wp points_lob;
    VectR2 points_tri, points_quad; VectR3 points_inside;
    Matrix<int> NumTri, NumQuad; Array3D<int> NumNodes;
    Fb_geom.ConstructLobattoPoints(order+1, points_lob, points_tri, points_inside);

    VectReal_wp points_dof1d; VectR3 points_dof3d;
    VectReal_wp weights_lob;
    ComputeGaussLegendre(points_dof1d, weights_lob, order-1);
    
    /*points_dof3d.Reallocate(nb_dof_loc);
    tangente_dof.Reallocate(nb_dof_loc);
    // dofs on edges
    for (int i = 0; i < order; i++)
      {
	points_dof3d(i).Init(points_dof1d(i), 0, 0);
	tangente_dof(i).Init(1, 0, 0);
	
	points_dof3d(order+i).Init(0, points_dof1d(i), 0);
	tangente_dof(order+i).Init(0, 1, 0);
            
	points_dof3d(2*order+i).Init(0, 0, points_dof1d(i));
	tangente_dof(2*order+i).Init(0, 0, 1);
	
	points_dof3d(3*order+i).Init(1.0-points_dof1d(i), points_dof1d(i), 0);
	tangente_dof(3*order+i).Init(-1, 1, 0);
        
	points_dof3d(4*order+i).Init(1.0-points_dof1d(i), 0, points_dof1d(i));
	tangente_dof(4*order+i).Init(-1, 0, 1);
        
	points_dof3d(5*order+i).Init(0, 1.0-points_dof1d(i), points_dof1d(i));
	tangente_dof(5*order+i).Init(0, -1, 1);
      }
    
    // dofs on faces
    int offset = 6*order; int Ntri = (order-1)*order;
    for (int i = 0; i < (order-1)*order/2; i++)
      {
	Real_wp x = points_tri(3*(order+1) + i)(0);
	Real_wp y = points_tri(3*(order+1) + i)(1);
	points_dof3d(offset).Init(x, y, 0);
	points_dof3d(offset + 1).Init(x, y, 0);
        
	tangente_dof(offset).Init(1, 0, 0);
	tangente_dof(offset + 1).Init(0, 1, 0);
	
	points_dof3d(Ntri+offset).Init(x, 0, y);
	points_dof3d(Ntri+offset + 1).Init(x, 0, y);
        
	tangente_dof(Ntri+offset).Init(1, 0, 0);
	tangente_dof(Ntri+offset + 1).Init(0, 0, 1);
        
	points_dof3d(2*Ntri+offset).Init(0, x, y);
	points_dof3d(2*Ntri+offset + 1).Init(0, x, y);
        
	tangente_dof(2*Ntri+offset).Init(0, 1, 0);
	tangente_dof(2*Ntri+offset + 1).Init(0, 0, 1);
	
	points_dof3d(3*Ntri+offset).Init(1.0-x-y, x, y);
	points_dof3d(3*Ntri+offset + 1).Init(1.0-x-y, x, y);
        
	tangente_dof(3*Ntri+offset).Init(-1, 1, 0);
	tangente_dof(3*Ntri+offset + 1).Init(-1, 0, 1);
        
	offset += 2;
      }
    
    // dofs inside
    offset += 3*Ntri;
    for (int i = 0; i < (order-2)*(order-1)*order/6; i++)
      {
	R3 pt = points_inside(2*(order+1)*(order+1)+2+i);
	points_dof3d(offset) = pt;
	points_dof3d(offset+1) = pt;
	points_dof3d(offset+2) = pt;
        
	tangente_dof(offset).Init(1, 0, 0);
	tangente_dof(offset+1).Init(0, 1, 0);
	tangente_dof(offset+2).Init(0, 0, 1);
	offset += 3;
      }
    */
    
    points_dof1d = this->Points1D();
    points_dof3d = this->PointsND();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;

    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDof2D_tri(this->Points2D_tri());
    this->SetPointsDofND(points_dof3d);
    
    TriangleHcurlOptimalHpFirstFamily* Fb_tri = new TriangleHcurlOptimalHpFirstFamily();
    element_tri_surf = Fb_tri;
    Fb_tri->ConstructFiniteElement(order);
    
    if (order > 2)
      {
	NumDofsX.Reallocate(order-2, order-2, order-2);
	NumDofsY.Reallocate(order-2, order-2, order-2);
	NumDofsZ.Reallocate(order-2, order-2, order-2);
	NumDofsX.Fill(-1); NumDofsY.Fill(-1); NumDofsZ.Fill(-1);
	int offset = 6*order + 4*(order-1)*order;
	for (int diag = 0; diag <= order-3; diag++)
	  for (int i = 0; i <= diag; i++)
	    for (int j = 0; j <= diag-i; j++)
	      for (int k = 0; k <= diag-i-j; k++)
		if (i+j+k == diag)
		  {
		    NumDofsX(i, j, k) = offset++;
		    NumDofsY(i, j, k) = offset++;
		    NumDofsZ(i, j, k) = offset++;
		  }
      }
    
  }
  

  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void TetrahedronHcurlOptimalHpFirstFamily
  ::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    //TetrahedronReference<2>::ComputeProjectionDofRef(feval, contrib);
    
    const VectReal_wp& weights3d = this->WeightsND();    
    Vector1 feval_weight = feval;
    for (int i = 0; i < feval_weight.GetM()/3; i++)
      {
	feval_weight(3*i) *= weights3d(i);
	feval_weight(3*i+1) *= weights3d(i);
	feval_weight(3*i+2) *= weights3d(i);
      }
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }

  
  void TetrahedronHcurlOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	contrib(2*this->order+i) = -contrib(2*this->order+i);
  }


  void TetrahedronHcurlOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	contrib(2*this->order+i) = -contrib(2*this->order+i);
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point_loc
  */
  void TetrahedronHcurlOptimalHpFirstFamily
  ::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;
    
    VectReal_wp P1, P2, P3, P4, P5, P6, Px, Py, Pz;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x+y+z-1, P1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y-x, P2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x+2*y+z-1, P3);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x+y+2*z-1, P4);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, z-x, P5);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, z-y, P6);
    
    // edges
    // Edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = (1-y-z)*P1(i);
        phi(node)(1) = x*P1(i);
        phi(node)(2) = x*P1(i);
        node++;
      }
    
    // Edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = y*P3(i);
        phi(node)(1) = (1-x-z)*P3(i);
        phi(node)(2) = y*P3(i);
        node++;
      }
    
    // Edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = z*P4(i);
        phi(node)(1) = z*P4(i);
        phi(node)(2) = (1-x-y)*P4(i);
        node++;
      }
    
    // Edge 4
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = -y*P2(i);
        phi(node)(1) = x*P2(i);
        phi(node)(2) = 0;
        node++;
      }
    
    // Edge 5
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = -z*P5(i);
        phi(node)(1) = 0;
        phi(node)(2) = x*P5(i);
        node++;
      }
    
    // Edge 6
    for (int i = 0; i <= order-1; i++)
      {	
        phi(node)(0) = 0;
        phi(node)(1) = -z*P6(i);
        phi(node)(2) = y*P6(i);
        node++;
      }
    
    // face triangulaire 1
    int offset = 3*order;
    const TriangleHcurlOptimalHpFirstFamily& Fb_tri
      = static_cast<const TriangleHcurlOptimalHpFirstFamily& >(*element_tri_surf);
    
    const Matrix<int>& NumDofsX_tri = Fb_tri.GetNumDofsX();
    const Matrix<int>& NumDofsY_tri = Fb_tri.GetNumDofsY();
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {		
	  node = offset + NumDofsX_tri(i, j);
          phi(node)(0) = x*y*P1(i)*P3(j);
          phi(node)(1) = x*(1-x-z)*P1(i)*P3(j);
          phi(node)(2) = x*y*P1(i)*P3(j);
          
	  node = offset + NumDofsY_tri(i, j);
          phi(node)(0) = y*(1-y-z)*P1(i)*P3(j);
          phi(node)(1) = y*x*P1(i)*P3(j);
          phi(node)(2) = y*x*P1(i)*P3(j);
        }	
    
    // face triangulaire 2
    offset += (order-1)*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
	  node = offset + NumDofsX_tri(i, j);
          phi(node)(0) = x*z*P1(i)*P4(j);
          phi(node)(1) = x*z*P1(i)*P4(j);
          phi(node)(2) = x*(1-x-y)*P1(i)*P4(j);
	  
	  node = offset + NumDofsY_tri(i, j);
          phi(node)(0) = z*(1-y-z)*P1(i)*P4(j);
          phi(node)(1) = z*x*P1(i)*P4(j);
          phi(node)(2) = z*x*P1(i)*P4(j);
	}
    
    // face triangulaire 3
    offset += (order-1)*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
	  node = offset + NumDofsX_tri(i, j);
          phi(node)(0) = y*z*P3(i)*P4(j);
          phi(node)(1) = y*z*P3(i)*P4(j);
          phi(node)(2) = y*(1-x-y)*P3(i)*P4(j);
          
	  node = offset + NumDofsY_tri(i, j);
          phi(node)(0) = z*y*P3(i)*P4(j);
          phi(node)(1) = z*(1-x-z)*P3(i)*P4(j);
          phi(node)(2) = z*y*P3(i)*P4(j);
        }
    
    // face triangulaire 4
    offset += (order-1)*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
	  node = offset + NumDofsX_tri(i, j);
          phi(node)(0) = -y*z*P2(i)*P5(j);
          phi(node)(1) = 0;
          phi(node)(2) = y*x*P2(i)*P5(j);
          
          node = offset + NumDofsY_tri(i, j);
          phi(node)(0) = -y*z*P2(i)*P5(j);
          phi(node)(1) = z*x*P2(i)*P5(j);
          phi(node)(2) = 0;
        }
    
    // interior
    for (int i = 0; i <= order-3; i++)
      for (int j = 0; j <= order-3-i; j++)
        for (int k = 0; k <= order-3-i-j; k++)
          {
	    node = NumDofsX(i, j, k);
            phi(node)(0) = x*y*z*P1(i)*P3(j)*P4(k);
            phi(node)(1) = x*y*z*P1(i)*P3(j)*P4(k);
            phi(node)(2) = x*y*(1-x-y)*P1(i)*P3(j)*P4(k);
            
	    node = NumDofsY(i, j, k);
            phi(node)(0) = y*z*(1-y-z)*P1(i)*P3(j)*P4(k);
            phi(node)(1) = y*z*x*P1(i)*P3(j)*P4(k);
            phi(node)(2) = y*z*x*P1(i)*P3(j)*P4(k);
            
	    node = NumDofsZ(i, j, k);
            phi(node)(0) = x*z*y*P1(i)*P3(j)*P4(k);
            phi(node)(1) = x*z*(1-x-z)*P1(i)*P3(j)*P4(k);
            phi(node)(2) = x*z*y*P1(i)*P3(j)*P4(k);
          }	
    
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void TetrahedronHcurlOptimalHpFirstFamily
  ::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    FillZero(res);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;
    
    VectReal_wp P1, P2, P3, P4, P5, P6;
    VectReal_wp dP1, dP2, dP3, dP4, dP5, dP6;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x+y+z-1, P1, dP1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y-x, P2, dP2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x+2*y+z-1, P3, dP3);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x+y+2*z-1, P4, dP4);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, z-x, P5, dP5);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, z-y, P6, dP6);

    const TriangleHcurlOptimalHpFirstFamily& Fb_tri
      = static_cast<const TriangleHcurlOptimalHpFirstFamily& >(*element_tri_surf);
    
    const Matrix<int>& NumDofsX_tri = Fb_tri.GetNumDofsX();
    const Matrix<int>& NumDofsY_tri = Fb_tri.GetNumDofsY();
    
    // edges
    // Edge 1
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = 0;
        res(node)(1) = dP1(i)-z*dP1(i)-y*dP1(i)-2*P1(i)-2*x*dP1(i);
        res(node)(2) = -dP1(i)+2*x*dP1(i)+y*dP1(i)+2*P1(i)+z*dP1(i);
        node++;
      }
    
    // Edge 2
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = -dP3(i)+2*y*dP3(i)+x*dP3(i)+2*P3(i)+z*dP3(i);
        res(node)(1) = 0;
        res(node)(2) =  dP3(i)-z*dP3(i)-x*dP3(i)-2*P3(i)-2*y*dP3(i);
        node++;
      }
    
    // Edge 3
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = dP4(i)-y*dP4(i)-x*dP4(i)-2*P4(i)-2*z*dP4(i);
        res(node)(1) = -dP4(i)+2*z*dP4(i)+x*dP4(i)+2*P4(i)+y*dP4(i);
        res(node)(2) = 0;
        node++;
      }
    
    // Edge 4
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = 0;
        res(node)(1) = 0;
        res(node)(2) = -x*dP2(i)+2*P2(i)+y*dP2(i);
        node++;
      }
    
    // Edge 5
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = 0;
        res(node)(1) = -z*dP5(i)-2*P5(i)+x*dP5(i);
        res(node)(2) = 0;
        node++;
      }
    
    // Edge 6
    for (int i = 0; i <= order-1; i++)
      {	
        res(node)(0) = -y*dP6(i)+2*P6(i)+z*dP6(i);
        res(node)(1) = 0;
        res(node)(2) = 0;
        node++;
      }
    
    // face triangulaire 1
    int offset = 3*order;
    R3 grad_PiPj, grad_lambda, vec_u, curl; Real_wp lambda;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {		
          lambda = x*P1(i)*P3(j);
          grad_PiPj(0) = 2.0*dP1(i)*P3(j) + P1(i)*dP3(j);
          grad_PiPj(1) = dP1(i)*P3(j) + 2.0*P1(i)*dP3(j);
          grad_PiPj(2) = dP1(i)*P3(j) + P1(i)*dP3(j);
          
          grad_lambda = x*grad_PiPj;
          grad_lambda(0) += P1(i)*P3(j);
          
          vec_u.Init(y, 1.0-x-z, y);
          
          TimesProd(grad_lambda, vec_u, curl);
                    
	  node = offset + NumDofsX_tri(i, j);
          res(node)(0) = 2.0*lambda + curl(0); 
          res(node)(1) = curl(1);
          res(node)(2) = -2.0*lambda + curl(2);
          
          lambda = y*P1(i)*P3(j);
          vec_u.Init(1.0-y-z, x, x);
          
          grad_lambda = y*grad_PiPj;
          grad_lambda(1) += P1(i)*P3(j);
          
          TimesProd(grad_lambda, vec_u, curl);
                    
	  node = offset + NumDofsY_tri(i, j);
          res(node)(0) = curl(0); 
          res(node)(1) = -2.0*lambda + curl(1);
          res(node)(2) = 2.0*lambda + curl(2);
        }	
    
    // face triangulaire 2
    offset += (order-1)*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
          lambda = x*P1(i)*P4(j);
          grad_PiPj(0) = 2.0*dP1(i)*P4(j) + P1(i)*dP4(j);
          grad_PiPj(1) = dP1(i)*P4(j) + P1(i)*dP4(j);
          grad_PiPj(2) = dP1(i)*P4(j) + 2.0*P1(i)*dP4(j);
          
          grad_lambda = x*grad_PiPj;
          grad_lambda(0) += P1(i)*P4(j);
          
          vec_u.Init(z, z, 1.0-x-y);
          
          TimesProd(grad_lambda, vec_u, curl);
          
	  node = offset + NumDofsX_tri(i, j);
          res(node)(0) = -2.0*lambda + curl(0); 
          res(node)(1) = 2.0*lambda + curl(1);
          res(node)(2) = curl(2);
          node++;
          
          grad_lambda = z*grad_PiPj;
          grad_lambda(2) += P1(i)*P4(j);
          
          lambda = z*P1(i)*P4(j);
          vec_u.Init(1.0-y-z, x, x);
          
          TimesProd(grad_lambda, vec_u, curl);
                    
	  node = offset + NumDofsY_tri(i, j);
          res(node)(0) = curl(0); 
          res(node)(1) = -2.0*lambda + curl(1);
          res(node)(2) = 2.0*lambda + curl(2);
        }
    
    // face triangulaire 3
    offset += (order-1)*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
          lambda = y*P3(i)*P4(j);
          grad_PiPj(0) = dP3(i)*P4(j) + P3(i)*dP4(j);
          grad_PiPj(1) = 2.0*dP3(i)*P4(j) + P3(i)*dP4(j);
          grad_PiPj(2) = dP3(i)*P4(j) + 2.0*P3(i)*dP4(j);
          
          grad_lambda = y*grad_PiPj;
          grad_lambda(1) += P3(i)*P4(j);
          
          vec_u.Init(z, z, 1.0-x-y);
          
          TimesProd(grad_lambda, vec_u, curl);
          
	  node = offset + NumDofsX_tri(i, j);
          res(node)(0) = -2.0*lambda + curl(0); 
          res(node)(1) = 2.0*lambda + curl(1);
          res(node)(2) = curl(2);
          
          grad_lambda = z*grad_PiPj;
          grad_lambda(2) += P3(i)*P4(j);
          
          lambda = z*P3(i)*P4(j);
          vec_u.Init(y, 1.0-x-z, y);
          
          TimesProd(grad_lambda, vec_u, curl);
          
	  node = offset + NumDofsY_tri(i, j);
          res(node)(0) = 2.0*lambda + curl(0); 
          res(node)(1) = curl(1);
          res(node)(2) = -2.0*lambda + curl(2);
        }
    
    // face triangulaire 4
    offset += (order-1)*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
          lambda = y*P2(i)*P5(j);
          grad_PiPj(0) = -dP2(i)*P5(j) - P2(i)*dP5(j);
          grad_PiPj(1) = dP2(i)*P5(j);
          grad_PiPj(2) = P2(i)*dP5(j);
          
          grad_lambda = y*grad_PiPj;
          grad_lambda(1) += P2(i)*P5(j);
          
          vec_u.Init(-z, 0, x);
          
          TimesProd(grad_lambda, vec_u, curl);
          
	  node = offset + NumDofsX_tri(i, j);
          res(node)(0) = curl(0); 
          res(node)(1) = -2.0*lambda + curl(1);
          res(node)(2) = curl(2);          
          
          lambda = z*P2(i)*P5(j);
          
          grad_lambda = z*grad_PiPj;
          grad_lambda(2) += P2(i)*P5(j);
          
          vec_u.Init(-y, x, 0);
          
          TimesProd(grad_lambda, vec_u, curl);
          
	  node = offset + NumDofsY_tri(i, j);
          res(node)(0) = curl(0); 
          res(node)(1) = curl(1);
          res(node)(2) = 2.0*lambda + curl(2);        
          node++;
        }
    
    // interior
    Real_wp PiPj;
    for (int i = 0; i <= order-3; i++)
      for (int j = 0; j <= order-3-i; j++)
        for (int k = 0; k <= order-3-i-j; k++)
          {
            PiPj = P1(i)*P3(j)*P4(k);
            
            grad_PiPj(0) = 2.0*dP1(i)*P3(j)*P4(k) + P1(i)*dP3(j)*P4(k) + P1(i)*P3(j)*dP4(k);
            grad_PiPj(1) = dP1(i)*P3(j)*P4(k) + 2.0*P1(i)*dP3(j)*P4(k) + P1(i)*P3(j)*dP4(k);
            grad_PiPj(2) = dP1(i)*P3(j)*P4(k) + P1(i)*dP3(j)*P4(k) + 2.0*P1(i)*P3(j)*dP4(k);
            
            lambda = x*y*PiPj;
            vec_u.Init(z, z, 1.0-x-y);
            
            grad_lambda = x*y*grad_PiPj;
            grad_lambda(0) += y*PiPj;
            grad_lambda(1) += x*PiPj;
            
            TimesProd(grad_lambda, vec_u, curl);
	    
	    node = NumDofsX(i, j, k);
            res(node)(0) = -2.0*lambda + curl(0); 
            res(node)(1) = 2.0*lambda + curl(1);
            res(node)(2) = curl(2);
                      
            lambda = y*z*PiPj;
            vec_u.Init(1.0-y-z, x, x);
            
            grad_lambda = y*z*grad_PiPj;
            grad_lambda(1) += z*PiPj;
            grad_lambda(2) += y*PiPj;
            
            TimesProd(grad_lambda, vec_u, curl);
          
	    node = NumDofsY(i, j, k);
            res(node)(0) = curl(0); 
            res(node)(1) = -2.0*lambda + curl(1);
            res(node)(2) = 2.0*lambda + curl(2);            
            
            lambda = x*z*PiPj;
            vec_u.Init(y, 1.0-x-z, y);
            
            grad_lambda = x*z*grad_PiPj;
            grad_lambda(0) += z*PiPj;
            grad_lambda(2) += x*PiPj;
            
            TimesProd(grad_lambda, vec_u, curl);
          
	    node = NumDofsZ(i, j, k);
            res(node)(0) = 2.0*lambda + curl(0); 
            res(node)(1) = curl(1);
            res(node)(2) = -2.0*lambda + curl(2);              
          }	
	  
  }


  //! displays details of class TetrahedronHcurlOptimalHpFirstFamily
  ostream& operator <<(ostream& out, const TetrahedronHcurlOptimalHpFirstFamily& e)
  {
    out<<static_cast<const TetrahedronReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
