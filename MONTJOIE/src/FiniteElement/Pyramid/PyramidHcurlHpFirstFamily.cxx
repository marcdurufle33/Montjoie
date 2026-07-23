#ifndef MONTJOIE_FILE_PYRAMID_HCURL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidHcurlHpFirstFamily::PyramidHcurlHpFirstFamily() : PyramidReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    type_basis = OPTIMAL;
    //type_basis = NIGAM_PHILLIPS_1;
  }
  

  //! how to number mesh
  void PyramidHcurlHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 2*order*(order-1));
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);
    nmap.SetNbDofHexahedron(order, 0);

    nmap.SetEqualEdgesDofSymmetry(order, order+1);
    
    FindHcurlSignEdge(nmap);
    
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_quad_surf);
    
  }


  size_t PyramidHcurlHpFirstFamily::GetMemorySize() const
  {
    size_t taille = PyramidReference<2>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize() + JacobiPolynom11.GetMemorySize();
    taille += Seldon::GetMemorySize(JacobiPolynomMij) + OperatorDofTri.GetMemorySize();
    taille += NumDofsX.GetMemorySize() + NumDofsY.GetMemorySize() + NumDofsZ.GetMemorySize();
    return taille;
  }
  
    
  //! constructing finite element
  void PyramidHcurlHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    // type_quad = PyramidQuadrature::QUADRATURE_JACOBI1;    
    PyramidReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
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
  void PyramidHcurlHpFirstFamily::ConstructFunctions()
  {
    //nb_dof_loc = order*(7 + 7*order + 2*order*order)/2;
    nb_dof_loc = order*(order+1)*(2*order+1)/2 + order + (order-1)*order + 2*order*(order+1);
    
    nb_dof_boundaries = 8*order + 4*order*(order-1) + 2*order*(order-1);
    nb_dof_quad = 4*order + 2*order*(order-1);
    nb_dof_tri = 3*order + order*(order-1);

    VectReal_wp points_dof1d; VectR3 points_dof3d;
    //points_dof3d.Reallocate(nb_points_quadrature_inside);
    //for (int i = 0; i < nb_points_quadrature_inside; i++)
    //points_dof3d(i) = points3d(i);
    
    VectR3 points_inside;
    VectR2 points_tri, points_quad;
    VectReal_wp points_lob, weights_lob;
    Matrix<int> NumQuad(order+2, order+2);
    NumQuad.Fill();
    ComputeGaussLegendre(points_dof1d, weights_lob, order-1);
    ComputeGaussLobatto(points_lob, weights_lob, order+1);
    Fb_geom.ConstructLobattoPoints(order+1, points_lob, points_tri,
				   points_quad, NumQuad, points_inside);
    
    VectR2 points_tri_bis;
    VectR3 points_inside_xy;
    NumQuad.Reallocate(order+1, order+1); NumQuad.Fill();
    ComputeGaussLobatto(points_lob, weights_lob, order);
    Fb_geom.ConstructLobattoPoints(order, points_lob, points_tri_bis,
				   points_quad, NumQuad, points_inside_xy);
    
    /* points_dof3d.Reallocate(nb_dof_loc);
    tangente_dof.Reallocate(nb_dof_loc);
    // dofs on edges
    for (int i = 0; i < order; i++)
      {
        Real_wp x = 2.0*points_dof1d(i) - 1.0;
        Real_wp xc = points_dof1d(i);
        points_dof3d(i).Init(x, -1.0, 0);
        tangente_dof(i).Init(2, 0, 0);
        
        points_dof3d(order+i).Init(1.0, x, 0);
        tangente_dof(order+i).Init(0, 2, 0);

        points_dof3d(2*order+i).Init(x, 1.0, 0);
        tangente_dof(2*order+i).Init(2, 0, 0);
        
        points_dof3d(3*order+i).Init(-1.0, x, 0);
        tangente_dof(3*order+i).Init(0, 2, 0);
        
        points_dof3d(4*order+i).Init(xc - 1.0, xc - 1.0, xc);
        tangente_dof(4*order+i).Init(1, 1, 1);
        
        points_dof3d(5*order+i).Init(1.0-xc, xc - 1.0, xc);
        tangente_dof(5*order+i).Init(-1, 1, 1);
        
        points_dof3d(6*order+i).Init(1.0-xc, 1.0-xc, xc);
        tangente_dof(6*order+i).Init(-1, -1, 1);
        
        points_dof3d(7*order+i).Init(xc-1.0, 1.0-xc, xc);
        tangente_dof(7*order+i).Init(1, -1, 1);
      }
    
    // dofs on quadrilateral base
    int offset = 8*order;
    for (int i = 0; i < order; i++)
      for (int j = 1; j < order; j++)
        {
          Real_wp x = 2.0*points_dof1d(i) - 1.0;
          Real_wp y = 2.0*points_lob(j) - 1.0;
          
          points_dof3d(offset).Init(x, y, 0);
          tangente_dof(offset).Init(2, 0, 0);
          offset++;
        }

    for (int i = 1; i < order; i++)
      for (int j = 0; j < order; j++)
        {
          Real_wp x = 2.0*points_lob(i) - 1.0;
          Real_wp y = 2.0*points_dof1d(j) - 1.0;
          
          points_dof3d(offset).Init(x, y, 0);
          tangente_dof(offset).Init(0, 2, 0);
          offset++;
        }
    
    // dofs on triangular faces
    int Ntri = (order-1)*order;
    for (int i = 0; i < (order-1)*order/2; i++)
      {
        Real_wp x = points_tri(3*(order+1) + i)(0);
        Real_wp y = points_tri(3*(order+1) + i)(1);
        points_dof3d(offset).Init(2.0*x + y -1.0, y-1.0, y);
        points_dof3d(offset + 1).Init(2.0*x + y -1.0, y-1.0, y);
        
        tangente_dof(offset).Init(2, 0, 0);
        tangente_dof(offset + 1).Init(1, 1, 1);
        
        points_dof3d(Ntri+offset).Init(1.0-y, 2.0*x + y -1.0, y);
        points_dof3d(Ntri+offset + 1).Init(1.0-y, 2.0*x + y -1.0, y);
        
        tangente_dof(Ntri+offset).Init(0, 2, 0);
        tangente_dof(Ntri+offset + 1).Init(-1, 1, 1);
        
        points_dof3d(2*Ntri+offset).Init(2.0*x+y-1.0, 1.0-y, y);
        points_dof3d(2*Ntri+offset + 1).Init(2.0*x+y-1.0, 1.0-y, y);
        
        tangente_dof(2*Ntri+offset).Init(2, 0, 0);
        tangente_dof(2*Ntri+offset + 1).Init(1, -1, 1);
        
        points_dof3d(3*Ntri+offset).Init(y-1.0, 2.0*x + y -1.0, y);
        points_dof3d(3*Ntri+offset + 1).Init(y-1.0, 2.0*x + y -1.0, y);
        
        tangente_dof(3*Ntri+offset).Init(0, 2, 0);
        tangente_dof(3*Ntri+offset + 1).Init(1, 1, 1);
        
        offset += 2;
      }
    
    // dofs inside
    offset += 3*Ntri;
    for (int i = 0; i < (order-1)*order*(2*order-1)/6; i++)
      {
	R3 pt = points_inside(3*(order+1)*(order+1)+2+i);
	points_dof3d(offset) = pt;
	tangente_dof(offset).Init(1, 0, 0);
	points_dof3d(offset+1) = pt;
	tangente_dof(offset+1).Init(0, 1, 0);
	points_dof3d(offset+2) = pt;
	tangente_dof(offset+2).Init(0, 0, 1);
	offset += 3;
      }
    */
    
    points_dof3d = this->PointsND();
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDof2D_tri(this->Points2D_tri());
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(points_dof3d);
    
    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;

    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    GetJacobiPolynomial(JacobiPolynom11, order, Real_wp(1), Real_wp(1));
    JacobiPolynomMij.Reallocate(order);
    for (int mij = 0; mij < order; mij++)
      GetJacobiPolynomial(JacobiPolynomMij(mij), order-2, Real_wp(2*mij+2), Real_wp(0));
    
    if (type_basis == NIGAM_PHILLIPS_1)
      {
        // nb_dof_loc = order*(7 + 7*order + 2*order*order)/2;
        nb_dof_loc = 3*order*order*order + 5*order;
      }
    
    TriangleHcurlOptimalHpFirstFamily* Fb_tri = new TriangleHcurlOptimalHpFirstFamily();
    QuadrangleHcurlHpFirstFamily* Fb_quad = new QuadrangleHcurlHpFirstFamily();
    Fb_tri->ConstructFiniteElement(order);
    Fb_quad->ConstructFiniteElement(order);
    
    element_tri_surf = Fb_tri;
    element_quad_surf = Fb_quad;

    NumDofsX.Reallocate(order-1, order-1, order-1);
    NumDofsY.Reallocate(order-1, order-1, order-1);
    NumDofsZ.Reallocate(order-1, order-1, order-1);
    int offset = nb_dof_boundaries;
    for (int diag = 0; diag <= order-2; diag++)
      for (int i = 0; i <= order-2; i++)
	for (int j = 0; j <= order-2; j++)
	  for (int k = 0; k <= order-2-max(i, j); k++)
	    if (max(i, j)+k == diag)
	      {
		NumDofsX(i, j, k) = offset++;
		NumDofsY(i, j, k) = offset++;
		NumDofsZ(i, j, k) = offset++;
	      }
  }
  
  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void PyramidHcurlHpFirstFamily
  ::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    //PyramidReference<2>::ComputeProjectionDofRef(feval, contrib);

    const VectReal_wp& weights3d = this->WeightsND();    
    Vector1 feval_weight = feval;
    for (int i = 0; i < feval_weight.GetM(); i+= 3)
      {
        int j = i/3;
        feval_weight(i) *= weights3d(j);
        feval_weight(i+1) *= weights3d(j);
        feval_weight(i+2) *= weights3d(j);
      }
    
    contrib.Reallocate(nb_dof_loc);
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  

  void PyramidHcurlHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    bool quad = (num_loc == 0);
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	{
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
 	  if (quad)
	    contrib(3*this->order+i) = -contrib(3*this->order+i);
	}
  }


  void PyramidHcurlHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    bool quad = (num_loc == 0);
    for (int i = 0; i < this->order; i++)
      if (i%2 == 0)
	{
	  contrib(2*this->order+i) = -contrib(2*this->order+i);
	  if (quad)
	    contrib(3*this->order+i) = -contrib(3*this->order+i);
	}
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point_loc
  */
  void PyramidHcurlHpFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    Real_wp omz, invOmz(0);
    omz = 1.0 - z;

    // parameters for the pyramid
    Real_wp beta1 = 0.5*(1.0-x-z);
    Real_wp beta2 = 0.5*(1.0-y-z);
    Real_wp beta3 = 0.5*(1.0+x-z);
    Real_wp beta4 = 0.5*(1.0+y-z);
        
    // coordinates on the cube
    Real_wp a(0), b(0), c;
    Real_wp lambda1(0), lambda2(0), lambda3(0), lambda4(0);
    if (abs(omz) > epsilon_machine)
      {
        invOmz = 1.0/(1.0-z);
        a = x*invOmz;
        b = y*invOmz;

        lambda1 = beta1*beta2*invOmz;
        lambda2 = beta2*beta3*invOmz;
        lambda3 = beta3*beta4*invOmz;
        lambda4 = beta4*beta1*invOmz;
      }      
    
    c = 2.0*z - 1.0;
    
    Real_wp lambda5 = z;
    Real_wp gamma1 = z + 0.5*(x+y);
    Real_wp gamma2 = z + 0.5*(-x+y);
    Real_wp gamma3 = z + 0.5*(-x-y);
    Real_wp gamma4 = z + 0.5*(x-y);
    
    VectReal_wp powOneMinusZ(order);
    powOneMinusZ(0) = 1.0;
    for (int i = 0; i < order-1; i++)
      powOneMinusZ(i+1) = powOneMinusZ(i)*omz;
    
    // values of Legendre polynomials for edges
    VectReal_wp Px, Py, P1, P2, P3, P4, Jx, Jy, Px2, Py2, dPx, dPy;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x, Px2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y, Py2);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, b, Py);

    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    Vector<VectReal_wp> Pz(order);
    for (int i = 0; i <= order-2; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-2, c, Pz(i));
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma1, P1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma2, P2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma3, P3);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma4, P4);
    
    // elementary functions for horizontal edges
    R3 v1, v2, v3, v4;
    v1.Init(0.5*(lambda1+lambda2), 0, 0.5*(-lambda1+lambda2));
    v2.Init(0, 0.5*(lambda2+lambda3), 0.5*(-lambda2+lambda3));
    v3.Init(0.5*(lambda3+lambda4), 0, -0.5*(-lambda3+lambda4));
    v4.Init(0, 0.5*(lambda4+lambda1), -0.5*(-lambda4+lambda1));
    
    // hierarchical functions on horizontal edges
    int node = 0;
    for (int i = 0; i < order; i++)
      {
        phi(node) = v1;
        Mlt(Px2(i), phi(node));
        node++;
      }
    
    for (int i = 0; i < order; i++)
      {
        phi(node) = v2;
        Mlt(Py2(i), phi(node));
        node++;
      }
    
    for (int i = 0; i < order; i++)
      {
        phi(node) = v3;
        Mlt(Px2(i), phi(node));
        node++;
      }

    for (int i = 0; i < order; i++)
      {
        phi(node) = v4;
        Mlt(Py2(i), phi(node));
        node++;
      }
    
    // elementary functions for vertical edges
    R3 grad_L1(-0.5*beta2, -0.5*beta1, -0.5*(beta2+beta1));
    grad_L1 *= invOmz; grad_L1(2) += beta1*beta2*invOmz*invOmz;
    
    R3 grad_L2(0.5*beta2, -0.5*beta3, -0.5*(beta3+beta2));
    grad_L2 *= invOmz; grad_L2(2) += beta2*beta3*invOmz*invOmz;
    
    R3 grad_L3(0.5*beta4, 0.5*beta3, -0.5*(beta4+beta3));
    grad_L3 *= invOmz; grad_L3(2) += beta3*beta4*invOmz*invOmz;
    
    R3 grad_L4(-0.5*beta4, 0.5*beta1, -0.5*(beta1+beta4));
    grad_L4 *= invOmz; grad_L4(2) += beta4*beta1*invOmz*invOmz;

    R3 v5, v6, v7, v8;
    v5 = grad_L1; Mlt(-z, v5); v5(2) += lambda1;
    v6 = grad_L2; Mlt(-z, v6); v6(2) += lambda2;
    v7 = grad_L3; Mlt(-z, v7); v7(2) += lambda3;
    v8 = grad_L4; Mlt(-z, v8); v8(2) += lambda4;
    
    // hierarchical functions on vertical edges
    for (int i = 0; i < order; i++)
      {
        phi(node) = v5;
        Mlt(P1(i), phi(node));
        node++;
      }
    
    for (int i = 0; i < order; i++)
      {
        phi(node) = v6;
        Mlt(P2(i), phi(node));
        node++;
      }
    
    for (int i = 0; i < order; i++)
      {
        phi(node) = v7;
        Mlt(P3(i), phi(node));
        node++;
      }

    for (int i = 0; i < order; i++)
      {
        phi(node) = v8;
        Mlt(P4(i), phi(node));
        node++;
      }

    const TriangleHcurlOptimalHpFirstFamily& Fb_tri
      = static_cast<const TriangleHcurlOptimalHpFirstFamily& >(*element_tri_surf);
    
    const Matrix<int>& NumDofsX_tri = Fb_tri.GetNumDofsX();
    const Matrix<int>& NumDofsY_tri = Fb_tri.GetNumDofsY();
        
    const QuadrangleHcurlHpFirstFamily& Fb_quad
      = static_cast<const QuadrangleHcurlHpFirstFamily& >(*element_quad_surf);
    
    const Matrix<int>& NumDofsX_quad = Fb_quad.GetNumDofsX();
    const Matrix<int>& NumDofsY_quad = Fb_quad.GetNumDofsY();

    // evaluating functions on quadrilateral base
    EvaluateJacobiPolynomial(JacobiPolynom11, order-1, a, Jx);
    EvaluateJacobiPolynomial(JacobiPolynom11, order-1, b, Jy);
    int offset = 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j < order-1; j++)
        {	
	  node = offset + NumDofsX_quad(i, j);
          phi(node) = v1;
          phi(node) *= beta4*Px(i)*Jy(j)*powOneMinusZ(max(i,j))*invOmz;
          
	  node = offset + NumDofsY_quad(j, i);
          phi(node) = v4;
          phi(node) *= beta3*Jx(j)*Py(i)*powOneMinusZ(max(i,j))*invOmz;
        }
    
    // functions on triangular faces
    offset += Fb_quad.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
          phi(node) = v5;
          phi(node) *= beta3*Px2(i)*P1(j);
	  
	  node = offset + NumDofsY_tri(i, j);
          phi(node) = v1;
          phi(node) *= lambda5*Px2(i)*P1(j);
	}

    offset += Fb_tri.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
          phi(node) = v6;
          phi(node) *= beta4*Py2(i)*P2(j);
          
	  node = offset + NumDofsY_tri(i, j);
          phi(node) = v2;
          phi(node) *= lambda5*Py2(i)*P2(j);
	}
    
    offset += Fb_tri.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
          phi(node) = v8;
          phi(node) *= beta3*Px2(i)*P4(j);
          
	  node = offset + NumDofsY_tri(i, j);
          phi(node) = v3;
          phi(node) *= lambda5*Px2(i)*P4(j);
	}
    
    offset += Fb_tri.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
          phi(node) = v5;
          phi(node) *= beta4*Py2(i)*P1(j);
          
	  node = offset + NumDofsY_tri(i, j);
          phi(node) = v4;
          phi(node) *= lambda5*Py2(i)*P1(j);
	}
    
    // interior functions
    if (type_basis == OPTIMAL)
      {
        for (int i = 0; i <= order-2; i++)
          for (int j = 0; j <= order-2; j++)
            {
              int m = max(i, j);
              for (int k = 0; k <= order-2-m; k++)
                {
                  Real_wp coef = Px(i)*Py(j)*powOneMinusZ(m)*Pz(m)(k)*invOmz;
                  
		  node = NumDofsX(i, j, k);
                  phi(node) = v1;
                  phi(node) *= coef*beta4*lambda5;              
                  
		  node = NumDofsY(i, j, k);
                  phi(node) = v4;
                  phi(node) *= coef*beta3*lambda5;
                  
		  node = NumDofsZ(i, j, k);
                  phi(node) = v5;
                  phi(node) *= coef*beta4*beta3;
		}
            }        
      }
    else if (type_basis == NIGAM_PHILLIPS_1)
      {
        Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
        Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
        Real_wp zi = z/(1.0-z);
        
        // DF^-1
        Matrix3_3 dfjm1; R3 vec_u;
        dfjm1(0, 0) = 0.5*invOmz;
        dfjm1(1, 1) = 0.5*invOmz;
        dfjm1(2, 2) = 1.0*invOmz*invOmz;
        dfjm1(0, 2) = 0.5*x*invOmz*invOmz;
        dfjm1(1, 2) = 0.5*y*invOmz*invOmz;
        
		
        for (int i = 0; i <= order-1; i++)
          for (int j = 0; j <= order-1; j++)
            for (int k = 0; k <= order-2; k++)
              {
                
                Real_wp coef = Px(i)*Py(j)*powOneMinusZ(order-1-k)*Pz(order-2-k)(k);
                Real_wp vx = 0.25*(1.0+a)*(1.0-a)*z;
                Real_wp vy = 0.25*(1.0+b)*(1.0-b)*z;
                Real_wp vz = (1.0-a)*(1.0+a)*(1.0-b)*(1.0+b)/16;
                
                if (j < order-1)
                  {
                    //vec_u.Init(yi*(1-yi)*zi/pow(1+zi,order+1)
                    // *pow(xi,i)*pow(yi,j)*pow(zi,k), 0, 0);
                    //MltTrans(dfjm1, vec_u, phi(node));
                    phi(node).Init(vy*coef, 0, a*vy*coef);
                    // phi(node).Init(vy*coef, 0, 0);
                    node++;
                  }
                
                if (i < order-1)
                  {
                    //vec_u.Init(0,xi*(1-xi)*zi/pow(1+zi,order+1)*pow(xi,i)*pow(yi,j)*pow(zi,k),0);
                    //MltTrans(dfjm1, vec_u, phi(node));
                    phi(node).Init(0, vx*coef, b*vx*coef);
                    node++;
                  }
                
                if ((i < order-1) && (j < order-1))
                  {
                    //vec_u.Init(0,0,xi*(1-xi)*yi*(1-yi)/pow(1+zi,order+1)
                    // *pow(xi,i)*pow(yi,j)*pow(zi,k));
                    //MltTrans(dfjm1, vec_u, phi(node));
                    phi(node).Init(0, 0, vz*coef);
                    node++;
                  }
              }
        
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*xi - 1.0, Px, dPx);
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*yi - 1.0, Py, dPy);

        for (int i = 0; i <= order-2; i++)
          for (int j = 0; j <= order-2; j++)
            {
              Real_wp coef = pow(zi/(1.0+zi), order-1)/square(1.0+zi);
              Real_wp valR = xi*(1.0-xi)*yi*(1.0-yi)*Px(i)*Py(j);
              Real_wp Rx = yi*(1.0-yi)*Py(j)*((1.0-2.0*xi)*Px(i) + 2.0*(1.0-xi)*xi*dPx(i));
              Real_wp Ry = xi*(1.0-xi)*Px(i)*((1.0-2.0*yi)*Py(j) + 2.0*(1.0-yi)*yi*dPy(j));
              
              vec_u.Init(Rx*zi, Ry*zi, -valR);                
              MltTrans(dfjm1, vec_u, phi(node));
              phi(node) *= coef;
              node++;
            }
      }
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void PyramidHcurlHpFirstFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
	  
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
	  
    Real_wp omz, invOmz(0);
    omz = 1-z;
	  
    // parameters for the pyramid
    Real_wp beta1 = 0.5*(1.0-x-z);
    Real_wp beta2 = 0.5*(1.0-y-z);
    Real_wp beta3 = 0.5*(1.0+x-z);
    Real_wp beta4 = 0.5*(1.0+y-z);
        
    // coordinates on the cube
    Real_wp a(0), b(0), da_dx(0), da_dz(0), db_dy(0), db_dz(0), c;
    Real_wp lambda1(0), lambda2(0), lambda3(0), lambda4(0);
    if (abs(omz) > epsilon_machine)
      {
        invOmz = 1.0/(1.0-z);
        a = x*invOmz;
        b = y*invOmz;
        da_dx = invOmz;
        db_dy = invOmz;
        da_dz = a*invOmz;
        db_dz = b*invOmz;
        
        lambda1 = beta1*beta2*invOmz;
        lambda2 = beta2*beta3*invOmz;
        lambda3 = beta3*beta4*invOmz;
        lambda4 = beta4*beta1*invOmz;
      }      
    
    c = 2.0*z - 1.0;
    
    Real_wp lambda5 = z;
    Real_wp gamma1 = z + 0.5*(x+y);
    Real_wp gamma2 = z + 0.5*(-x+y);
    Real_wp gamma3 = z + 0.5*(-x-y);
    Real_wp gamma4 = z + 0.5*(x-y);
    
    R3 dgamma1(0.5, 0.5, 1.0);
    R3 dgamma2(-0.5, 0.5, 1.0);
    R3 dgamma3(-0.5, -0.5, 1.0);
    R3 dgamma4(0.5, -0.5, 1.0);
    
    VectReal_wp powOneMinusZ(order), dpowOneMinusZ(order);
    powOneMinusZ(0) = 1.0;
    dpowOneMinusZ(0) = 0.0;
    for (int i = 0; i < order-1; i++)
      {
	powOneMinusZ(i+1) = powOneMinusZ(i)*omz;
	dpowOneMinusZ(i+1) = -Real_wp(i+1)*powOneMinusZ(i);
      }
    
    // values of Legendre polynomials for edges
    VectReal_wp Px, Py, P1, P2, P3, P4, Jx, Jy, Px2, Py2;
    VectReal_wp dPx, dPy, dP1, dP2, dP3, dP4, dJx, dJy, dPx2, dPy2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x, Px2, dPx2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y, Py2, dPy2);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, b, Py, dPy);

    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    Vector<VectReal_wp> Pz(order), dPz(order);
    for (int i = 0; i <= order-2; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), order-2, c, Pz(i), dPz(i));
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma1, P1, dP1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma2, P2, dP2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma3, P3, dP3);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, gamma4, P4, dP4);
    
    // elementary functions for horizontal edges
    R3 v1, v2, v3, v4;
    v1.Init(0.5*(lambda1+lambda2), 0, 0.5*(-lambda1+lambda2));
    v2.Init(0, 0.5*(lambda2+lambda3), 0.5*(-lambda2+lambda3));
    v3.Init(0.5*(lambda3+lambda4), 0, -0.5*(-lambda3+lambda4));
    v4.Init(0, 0.5*(lambda4+lambda1), -0.5*(-lambda4+lambda1));
        
    // elementary functions for vertical edges
    R3 grad_L1(-0.5*beta2, -0.5*beta1, -0.5*(beta2+beta1));
    grad_L1 *= invOmz; grad_L1(2) += beta1*beta2*invOmz*invOmz;
    
    R3 grad_L2(0.5*beta2, -0.5*beta3, -0.5*(beta3+beta2));
    grad_L2 *= invOmz; grad_L2(2) += beta2*beta3*invOmz*invOmz;
    
    R3 grad_L3(0.5*beta4, 0.5*beta3, -0.5*(beta4+beta3));
    grad_L3 *= invOmz; grad_L3(2) += beta3*beta4*invOmz*invOmz;
    
    R3 grad_L4(-0.5*beta4, 0.5*beta1, -0.5*(beta1+beta4));
    grad_L4 *= invOmz; grad_L4(2) += beta4*beta1*invOmz*invOmz;
    
    R3 v5, v6, v7, v8;
    v5 = grad_L1; Mlt(-z, v5); v5(2) += lambda1;
    v6 = grad_L2; Mlt(-z, v6); v6(2) += lambda2;
    v7 = grad_L3; Mlt(-z, v7); v7(2) += lambda3;
    v8 = grad_L4; Mlt(-z, v8); v8(2) += lambda4;
    
    // curl of elementary functions
    R3 curl_v1, curl_v2, curl_v3, curl_v4;
    R3 curl_v5, curl_v6, curl_v7, curl_v8;
    curl_v1.Init(0.5*(-grad_L1(1)+grad_L2(1)), 0.5*(grad_L1(2)+grad_L2(2)+grad_L1(0)-grad_L2(0)),
                 -0.5*(grad_L1(1)+grad_L2(1)));
    
    curl_v2.Init(0.5*(-grad_L2(1)+grad_L3(1)-grad_L2(2)-grad_L3(2)), 0.5*(grad_L2(0)-grad_L3(0)),
                 0.5*(grad_L2(0)+grad_L3(0)));
    
    curl_v3.Init(0.5*(grad_L3(1)-grad_L4(1)), 0.5*(grad_L3(2)+grad_L4(2)-grad_L3(0)+grad_L4(0)),
                 -0.5*(grad_L3(1)+grad_L4(1)));
    
    curl_v4.Init(0.5*(grad_L4(1)-grad_L1(1)-grad_L4(2)-grad_L1(2)), 0.5*(-grad_L4(0)+grad_L1(0)),
                 0.5*(grad_L4(0)+grad_L1(0)));
    
    curl_v5(0) = 2.0*grad_L1(1); curl_v5(1) = -2.0*grad_L1(0);
    curl_v6(0) = 2.0*grad_L2(1); curl_v6(1) = -2.0*grad_L2(0);
    curl_v7(0) = 2.0*grad_L3(1); curl_v7(1) = -2.0*grad_L3(0);
    curl_v8(0) = 2.0*grad_L4(1); curl_v8(1) = -2.0*grad_L4(0);
    
    // hierarchical functions on horizontal edges
    int node = 0;
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v1;
        Mlt(Px2(i), res(node));
        res(node)(1) -= dPx2(i)*v1(2);
        res(node)(2) += dPx2(i)*v1(1);
        node++;
      }
    
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v2;
        Mlt(Py2(i), res(node));
        res(node)(0) += dPy2(i)*v2(2);
        res(node)(2) -= dPy2(i)*v2(0);
        node++;
      }
    
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v3;
        Mlt(Px2(i), res(node));
        res(node)(1) -= dPx2(i)*v3(2);
        res(node)(2) += dPx2(i)*v3(1);
        node++;
      }

    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v4;
        Mlt(Py2(i), res(node));
        res(node)(0) += dPy2(i)*v4(2);
        res(node)(2) -= dPy2(i)*v4(0);
        node++;
      }

    // hierarchical functions on vertical edges
    Real_wp val; R3 grad_val;
    TimesProd(dgamma1, v5, grad_val); 
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v5;
        Mlt(P1(i), res(node));
	Add(dP1(i), grad_val, res(node));
        node++;
      }
    
    TimesProd(dgamma2, v6, grad_val); 
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v6;
        Mlt(P2(i), res(node));
	Add(dP2(i), grad_val, res(node));
        node++;
      }
    
    TimesProd(dgamma3, v7, grad_val); 
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v7;
        Mlt(P3(i), res(node));
	Add(dP3(i), grad_val, res(node));
        node++;
      }

    TimesProd(dgamma4, v8, grad_val); 
    for (int i = 0; i < order; i++)
      {
        res(node) = curl_v8;
        Mlt(P4(i), res(node));
	Add(dP4(i), grad_val, res(node));
        node++;
      }

    const TriangleHcurlOptimalHpFirstFamily& Fb_tri
      = static_cast<const TriangleHcurlOptimalHpFirstFamily& >(*element_tri_surf);
    
    const Matrix<int>& NumDofsX_tri = Fb_tri.GetNumDofsX();
    const Matrix<int>& NumDofsY_tri = Fb_tri.GetNumDofsY();

    const QuadrangleHcurlHpFirstFamily& Fb_quad
      = static_cast<const QuadrangleHcurlHpFirstFamily& >(*element_quad_surf);
    
    const Matrix<int>& NumDofsX_quad = Fb_quad.GetNumDofsX();
    const Matrix<int>& NumDofsY_quad = Fb_quad.GetNumDofsY();
    
    // evaluating functions on quadrilateral base
    EvaluateJacobiPolynomial(JacobiPolynom11, order-1, a, Jx, dJx);
    EvaluateJacobiPolynomial(JacobiPolynom11, order-1, b, Jy, dJy);
    Real_wp vz = beta4*invOmz;
    R3 dvz(0, 0.5*invOmz, invOmz*(-0.5+beta4*invOmz));

    Real_wp vz2 = beta3*invOmz;
    R3 dvz2(0.5*invOmz, 0, invOmz*(-0.5+beta3*invOmz));
    int offset = 4*order;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j < order-1; j++)
        {	
	  int m = max(i,j);
	  node = offset + NumDofsX_quad(i, j);
	  val = vz*Px(i)*Jy(j)*powOneMinusZ(m);	  
	  grad_val.Init(da_dx*dPx(i)*vz*Jy(j)*powOneMinusZ(m),
			Px(i)*powOneMinusZ(m)*(dvz(1)*Jy(j) + vz*dJy(j)*db_dy),
			vz*powOneMinusZ(m)*(dPx(i)*da_dz*Jy(j) + Px(i)*dJy(j)*db_dz)
			+ Px(i)*Jy(j)*(vz*dpowOneMinusZ(m) + dvz(2)*powOneMinusZ(m)));
	  
	  TimesProd(grad_val, v1, res(node));
	  Add(val, curl_v1, res(node));
          
	  node = offset + NumDofsY_quad(j, i);
	  val = vz2*Jx(j)*Py(i)*powOneMinusZ(m);	  
	  grad_val.Init(Py(i)*powOneMinusZ(m)*(da_dx*dJx(j)*vz2 + Jx(j)*dvz2(0)),
			Jx(j)*powOneMinusZ(m)*vz2*dPy(i)*db_dy,
			vz2*powOneMinusZ(m)*(dJx(j)*da_dz*Py(i) + Jx(j)*dPy(i)*db_dz)
			+ Jx(j)*Py(i)*(vz2*dpowOneMinusZ(m) + dvz2(2)*powOneMinusZ(m)));
	  
	  TimesProd(grad_val, v4, res(node));
	  Add(val, curl_v4, res(node));
        }
    
    // functions on triangular faces
    offset += Fb_quad.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
	  val = beta3*Px2(i)*P1(j);
	  grad_val = dgamma1; Mlt(beta3*Px2(i)*dP1(j), grad_val);
	  grad_val(0) += P1(j)*(0.5*Px2(i) + beta3*dPx2(i));
	  grad_val(2) -= 0.5*Px2(i)*P1(j);
          
	  TimesProd(grad_val, v5, res(node));
	  Add(val, curl_v5, res(node));
	  
	  node = offset + NumDofsY_tri(i, j);
          val = lambda5*Px2(i)*P1(j);
	  grad_val = dgamma1; Mlt(lambda5*Px2(i)*dP1(j), grad_val);
	  grad_val(0) += P1(j)*lambda5*dPx2(i);
	  grad_val(2) += Px2(i)*P1(j);
	  
	  TimesProd(grad_val, v1, res(node));
	  Add(val, curl_v1, res(node));
        }
    
    offset += Fb_tri.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
	  val = beta4*Py2(i)*P2(j);
	  grad_val = dgamma2; Mlt(beta4*Py2(i)*dP2(j), grad_val);
	  grad_val(1) += P2(j)*(0.5*Py2(i) + beta4*dPy2(i));
	  grad_val(2) -= 0.5*Py2(i)*P2(j);
	  
	  TimesProd(grad_val, v6, res(node));
	  Add(val, curl_v6, res(node));
	  
	  node = offset + NumDofsY_tri(i, j);
          val = lambda5*Py2(i)*P2(j);
	  grad_val = dgamma2; Mlt(lambda5*Py2(i)*dP2(j), grad_val);
	  grad_val(1) += P2(j)*lambda5*dPy2(i);
	  grad_val(2) += Py2(i)*P2(j);
          
	  TimesProd(grad_val, v2, res(node));
	  Add(val, curl_v2, res(node));
	}
    
    offset += Fb_tri.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
          val = beta3*Px2(i)*P4(j);
	  grad_val = dgamma4; Mlt(beta3*Px2(i)*dP4(j), grad_val);
	  grad_val(0) += P4(j)*(0.5*Px2(i) + beta3*dPx2(i));
	  grad_val(2) -= 0.5*Px2(i)*P4(j);
	  
	  TimesProd(grad_val, v8, res(node));
	  Add(val, curl_v8, res(node));
          
	  node = offset + NumDofsY_tri(i, j);
          val = lambda5*Px2(i)*P4(j);
	  grad_val = dgamma4; Mlt(lambda5*Px2(i)*dP4(j), grad_val);
	  grad_val(0) += P4(j)*lambda5*dPx2(i);
	  grad_val(2) += Px2(i)*P4(j);
	  
	  TimesProd(grad_val, v3, res(node));
	  Add(val, curl_v3, res(node));
	}
    
    offset += Fb_tri.GetNbDof() - 3*order;
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {
	  node = offset + NumDofsX_tri(i, j);
	  val = beta4*Py2(i)*P1(j);
	  grad_val = dgamma1; Mlt(beta4*Py2(i)*dP1(j), grad_val);
	  grad_val(1) += P1(j)*(0.5*Py2(i) + beta4*dPy2(i));
	  grad_val(2) -= 0.5*Py2(i)*P1(j);
	  
	  TimesProd(grad_val, v5, res(node));
	  Add(val, curl_v5, res(node));
          
	  node = offset + NumDofsY_tri(i, j);
          val = lambda5*Py2(i)*P1(j);
	  grad_val = dgamma1; Mlt(lambda5*Py2(i)*dP1(j), grad_val);
	  grad_val(1) += P1(j)*lambda5*dPy2(i);
	  grad_val(2) += Py2(i)*P1(j);
          
	  TimesProd(grad_val, v4, res(node));
	  Add(val, curl_v4, res(node));
        }
    
    // interior functions
    if (type_basis == OPTIMAL)
      {
        for (int i = 0; i <= order-2; i++)
          for (int j = 0; j <= order-2; j++)
            {
              int m = max(i, j);
              vz = Px(i)*Py(j)*powOneMinusZ(m)*invOmz;
              dvz.Init(powOneMinusZ(m)*invOmz*dPx(i)*Py(j)*da_dx,
                       powOneMinusZ(m)*invOmz*Px(i)*dPy(j)*db_dy,
                       powOneMinusZ(m)*invOmz*(dPx(i)*Py(j)*da_dz + Px(i)*dPy(j)*db_dz)
                       + Px(i)*Py(j)*(dpowOneMinusZ(m)*invOmz + powOneMinusZ(m)*invOmz*invOmz));
              
              for (int k = 0; k <= order-2-m; k++)
                {
                  vz2 = vz*Pz(m)(k);
                  dvz2 = dvz; dvz2 *= Pz(m)(k);
                  dvz2(2) += vz*2.0*dPz(m)(k);
                  
                  //if (j+k < order - 2)
                  //{
		  node = NumDofsX(i, j, k);
                  val = vz2*beta4*lambda5;
                  grad_val = dvz2; grad_val *= beta4*lambda5;
                  grad_val(1) += 0.5*vz2*lambda5;
                  grad_val(2) += vz2*(-0.5*lambda5 + beta4);
                  
                  TimesProd(grad_val, v1, res(node));
                  Add(val, curl_v1, res(node));
                  node++;
                  //}
              
                  //if (i+k < order - 2)
                  //{
		  node = NumDofsY(i, j, k);
                  val = vz2*beta3*lambda5;
                  grad_val = dvz2; grad_val *= beta3*lambda5;
                  grad_val(0) += 0.5*vz2*lambda5;
                  grad_val(2) += vz2*(-0.5*lambda5 + beta3);
                  
                  TimesProd(grad_val, v4, res(node));
                  Add(val, curl_v4, res(node));
                  node++;
                  //}
                  
		  node = NumDofsZ(i, j, k);
                  val = vz2*beta3*beta4;
                  grad_val = dvz2; grad_val *= beta3*beta4;
                  grad_val(0) += 0.5*vz2*beta4;
                  grad_val(1) += 0.5*vz2*beta3;
                  grad_val(2) += -0.5*vz2*(beta3+beta4);
                  
                  TimesProd(grad_val, v5, res(node));
                  Add(val, curl_v5, res(node));
                  node++;
                }
            }
        
      }
    else if (type_basis == NIGAM_PHILLIPS_1)
      {
        Real_wp xi = 0.5*(x/(1.0-z) + 1.0);
        Real_wp yi = 0.5*(y/(1.0-z) + 1.0);
        Real_wp zi = z/(1.0-z);
        
        // DF^-1
        Matrix3_3 dfjm1, invJi_dfj; R3 vec_u;
        dfjm1(0, 0) = 0.5*invOmz;
        dfjm1(1, 1) = 0.5*invOmz;
        dfjm1(2, 2) = 1.0*invOmz*invOmz;
        dfjm1(0, 2) = 0.5*x*invOmz*invOmz;
        dfjm1(1, 2) = 0.5*y*invOmz*invOmz;
        
        invJi_dfj = dfjm1; GetInverse(invJi_dfj);
        Real_wp invJi = 1.0/Det(invJi_dfj);
        Mlt(invJi, invJi_dfj);

        Real_wp vx = 0.25*(1.0+a)*(1.0-a)*z;
        Real_wp vy = 0.25*(1.0+b)*(1.0-b)*z;
        Real_wp vz = (1.0-a)*(1.0+a)*(1.0-b)*(1.0+b)/16;
        R3 grad_vx(da_dx*(-0.5*a*z), 0, da_dz*(-0.5*a*z) + 0.25*(1.0+a)*(1.0-a));
        R3 grad_vy(0, db_dy*(-0.5*b*z), db_dz*(-0.5*b*z) + 0.25*(1.0+b)*(1.0-b));
        R3 grad_vz(-da_dx*a*(1.0-b)*(1.0+b)/8, -db_dy*b*(1.0-a)*(1.0+a)/8,
                   -da_dz*a*(1.0-b)*(1.0+b)/8 -db_dz*b*(1.0-a)*(1.0+a)/8);
        
        for (int i = 0; i <= order-1; i++)
          for (int j = 0; j <= order-1; j++)
            for (int k = 0; k <= order-2; k++)
              {
                Real_wp valz = powOneMinusZ(order-1-k)*Pz(order-2-k)(k);
                Real_wp dvalz = -Real_wp(order-1-k)*powOneMinusZ(order-2-k)*Pz(order-2-k)(k)
                  + 2.0*powOneMinusZ(order-1-k)*dPz(order-2-k)(k);
                
                Real_wp coef = Px(i)*Py(j)*valz;
                R3 grad_coef(da_dx*dPx(i)*Py(j)*valz,
                             db_dy*Px(i)*dPy(j)*valz,
                             valz*(da_dz*dPx(i)*Py(j)+db_dz*Px(i)*dPy(j)) + dvalz*Px(i)*Py(j)); 
                                
                if (j < order-1)
                  {
                    //res(node).Init(0, grad_vy(2)*coef + grad_coef(2)*vy,
                    //-grad_vy(1)*coef - grad_coef(1)*vy);
                    res(node).Init(a*(grad_vy(1)*coef + grad_coef(1)*vy),
                                   -a*(grad_vy(0)*coef + grad_coef(0)*vy)
                                   - da_dx*vy*coef + grad_vy(2)*coef + grad_coef(2)*vy,
                                   -grad_vy(1)*coef - grad_coef(1)*vy);
                    node++;
                  }
                
                if (i < order-1)
                  {
                    res(node).Init(b*(grad_vx(1)*coef + grad_coef(1)*vx) + db_dy*vx*coef
                                   - grad_vx(2)*coef - grad_coef(2)*vx,
                                   -b*(grad_vx(0)*coef + grad_coef(0)*vx),
                                   grad_vx(0)*coef + grad_coef(0)*vx);
                    node++;
                  }
                
                if ((i < order-1) && (j < order-1))
                  {
                    res(node).Init(grad_vz(1)*coef + grad_coef(1)*vz,
                                   -grad_vz(0)*coef - grad_coef(0)*vz, 0);
                    node++;
                  }
              }
        
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*xi - 1.0, Px, dPx);
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*yi - 1.0, Py, dPy);

        for (int i = 0; i <= order-2; i++)
          for (int j = 0; j <= order-2; j++)
            {
              Real_wp coef = pow(zi/(1.0+zi), order-1)/square(1.0+zi);
              Real_wp der_coef = pow(zi/(1.0+zi), order-2)/pow(1.0+zi, 4)
                *(Real_wp(order-1) - 2.0*zi);
              
              Real_wp Rx = yi*(1.0-yi)*Py(j)*((1.0-2.0*xi)*Px(i) + 2.0*(1.0-xi)*xi*dPx(i));
              Real_wp Ry = xi*(1.0-xi)*Px(i)*((1.0-2.0*yi)*Py(j) + 2.0*(1.0-yi)*yi*dPy(j));
              vec_u.Init(-2.0*Ry*coef - Ry*zi*der_coef, 2.0*Rx*coef + Rx*zi*der_coef, 0);
              Mlt(invJi_dfj, vec_u, res(node));              
              node++;
            }
      }
  }


  //! displays details of class PyramidHcurlHpFirstFamily
  ostream& operator <<(ostream& out, const PyramidHcurlHpFirstFamily& e)
  {
    out<<static_cast<const PyramidReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HCURL_HP_FIRST_FAMILY_CXX
#endif
