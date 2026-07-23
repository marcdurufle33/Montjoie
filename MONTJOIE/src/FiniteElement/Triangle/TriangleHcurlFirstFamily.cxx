#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TriangleHcurlFirstFamily::TriangleHcurlFirstFamily() : TriangleReference<2>()
  {
    // dof, nodal and quadrature points are all different
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    // nodal points used to express basis functions
    // NODAL_REGULAR -> regular points (Graglia expressions)
    // NODAL_GAUSS -> gauss points along edges (for conformity with quadrilaterals)
    type_nodal_basis = NODAL_GAUSS;
    //type_nodal_basis = NODAL_REGULAR;
  }
  
  
  //! how to number the mesh
  void TriangleHcurlFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TriangleReference<2>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 0);
    nmap.SetNbDofTriangle(order, nb_dof_loc - 3*order);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t TriangleHcurlFirstFamily::GetMemorySize() const
  {
    size_t taille = TriangleReference<2>::GetMemorySize();
    taille += NumDofs_S1.GetMemorySize();
    taille += NumDofs_S2.GetMemorySize();
    taille += NumDofs_S3.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += InverseVDM_Edge.GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    return taille;
  }
  
  
  //! construction of finite element
  void TriangleHcurlFirstFamily::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			 int rsurf, int type_surf)
  {
    TriangleReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    ConstructFunctions();    
    
    this->Fb_geom.ComputeCoefficientTransformation();
    ConstructHcurlElementaryMatrix();    

    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r-1, r, r, type_quad, EdgeGauss::GAUSS);
    edge->SetDofPoints(this->PointsDof1D());
    this->element_surface = edge;

  }
  

  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  void TriangleHcurlFirstFamily
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
  }
  
  
  //! construction of basis functions
  void TriangleHcurlFirstFamily::ConstructFunctions()
  {
    // dofs on the three edges and inside dofs (two by inside nodes)
    nb_dof_loc = order*3 + (order-1)*order;
    nb_dof_boundaries = 3*order;
    // this->nb_dof_H = order*(order+1)/2; uncomment when ComputeValuesPhiHRef is completed
    
    EdgesDof.Reallocate(order, 3);
    NumDofs_S1.Reallocate(order+2, order+2);
    NumDofs_S1.Fill(-1);
    NumDofs_S2.Reallocate(order+2, order+2);
    NumDofs_S2.Fill(-1);
    NumDofs_S3.Reallocate(order+2, order+2);
    NumDofs_S3.Fill(-1);
    CoordinateDofs.Reallocate(nb_dof_loc, 4);
    CoordinateDofs.Fill(-1);
    
    // First Edge
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 0) = i;
        CoordinateDofs(i, 0) = 3;
	CoordinateDofs(i, 1) = order-i;
	CoordinateDofs(i, 2) = i+1;
	CoordinateDofs(i, 3) = 0;
	NumDofs_S3(CoordinateDofs(i, 2), CoordinateDofs(i, 3)) = i;
      }
    
    // Second Edge
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 1) = i + order;
        CoordinateDofs(i+order, 0) = 1;
	CoordinateDofs(i+order, 1) = 0;
	CoordinateDofs(i+order, 2) = order-i;
	CoordinateDofs(i+order, 3) = i+1;
	NumDofs_S1(CoordinateDofs(i+order, 2),
		   CoordinateDofs(i+order, 3)) = i + order;
      }
    
    // Third Edge
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 2) = i + 2*order;	      
        CoordinateDofs(i+2*order, 0) = 2;
	CoordinateDofs(i+2*order, 1) = i+1;
	CoordinateDofs(i+2*order, 2) = 0;
	CoordinateDofs(i+2*order, 3) = order-i;
	NumDofs_S2(CoordinateDofs(i+2*order, 2),
		   CoordinateDofs(i+2*order, 3)) = i+2*order;
      }
    
    int numero = 0;
    numero = 3*order;
    nb_dof_boundaries = numero;
    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);    

    for (int i = 1; i < order; i++)
      for (int j = 1; j < (order+1-i); j++)
	{
	  NumDofs_S3(i,j) = numero++;
	  NumDofs_S2(i,j) = numero++;
	  CoordinateDofs(numero-1,0) = 2;
	  CoordinateDofs(numero-2,0) = 3;
	  CoordinateDofs(numero-1,1) = order+1-i-j;
	  CoordinateDofs(numero-2,1) = order+1-i-j;
	  CoordinateDofs(numero-1,2) = i;
	  CoordinateDofs(numero-2,2) = i;
	  CoordinateDofs(numero-1,3) = j;
	  CoordinateDofs(numero-2,3) = j;
        }
    
    VectReal_wp points_lob;
    VectReal_wp points_dof1d; VectR2 points_dof2d;
    points_dof1d.Reallocate(order);    
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();    
    if (type_nodal_basis != NODAL_REGULAR)
      {
        // location of dofs on the edges
        VectReal_wp poids;
        ComputeGaussLegendre(points_dof1d, poids, order-1);
        if (type_nodal_basis == NODAL_LOBATTO)
          {
            ComputeGaussLobatto(points_lob, poids, order+1);
            points_dof1d.Reallocate(order);
            for (int i = 0; i < order; i++)
              points_dof1d(i) = points_lob(i+1);
          }
        
        // weight for orthonormalization
        InvWeightBasisRr.Reallocate(nb_dof_loc);
        int num = 0;
        for (int i = 0; i < order; i++)
          {
            for (int j = 0; j <= order-1-i; j++)
              {
                Real_wp coef = CoefLegendre(i)*CoefOddJacobi(i, j);
                InvWeightBasisRr(num) = coef;
                InvWeightBasisRr(num+1) = coef;
                num += 2;
              }
            
            InvWeightBasisRr(num) = CoefLegendre(i)*CoefOddJacobi(i, order-1-i);
            num++;
          }
        
        InvWeightBasisRr.Fill(1.0);
        VectReal_wp Scale(nb_dof_loc);
        VectR2 phi(nb_dof_loc); Scale.Fill(0);
        for (int k = 0; k < nb_points_quadrature_inside; k++)
          {
            ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
            for (int i = 0; i < nb_dof_loc; i++)
              Scale(i) += DotProd(phi(i), phi(i))*this->WeightsND(k);
          }
        
        for (int i = 0; i < nb_dof_loc; i++)
          InvWeightBasisRr(i) = 1.0/sqrt(Scale(i));
        
        // points where dofs are defined, and tangent on these points
        VectR2 points_inside;
	tangente_dof.Reallocate(nb_dof_loc);
        Fb_geom.ConstructLobattoPoints(order+1, 0, points_lob, points_inside);
        points_dof2d.Reallocate(nb_dof_loc);
        for (int i = 0; i < order; i++)
          {
            points_dof2d(i).Init(points_dof1d(i), 0);
            tangente_dof(i).Init(1, 0);
            
            points_dof2d(order+i).Init(1.0-points_dof1d(i), points_dof1d(i));
            tangente_dof(order+i).Init(-1, 1);

            points_dof2d(2*order+i).Init(0, 1.0-points_dof1d(i));
            tangente_dof(2*order+i).Init(0, -1);
          }
        
        for (int i = 0; i < (order-1)*order/2; i++)
          {
            points_dof2d(3*order + 2*i) = points_inside(3*(order+1) + i);
            points_dof2d(3*order + 2*i+1) = points_inside(3*(order+1) + i);

            tangente_dof(3*order + 2*i).Init(1, 0);
            tangente_dof(3*order + 2*i+1).Init(0, 1);
          }
        
        // computation of VDM = psi_j(xi_k) \cdot t_k
        InverseVDM_Edge.Reallocate(nb_dof_loc, nb_dof_loc);
        for (int k = 0; k < nb_dof_loc; k++)
          {
            ComputeValuesPhiOrthoRef(points_dof2d(k), phi);
            for (int j = 0; j < nb_dof_loc; j++)
              InverseVDM_Edge(j, k) = DotProd(phi(j), tangente_dof(k));
          }
        
        //InverseVDM_Edge.WriteText("VdmTri"+to_str(order)+".dat");
        GetInverse(InverseVDM_Edge);
        
        // then checking phi_i(xi_k) \cdot t_k = delta_{i, k}
        for (int k = 0; k < nb_dof_loc; k++)
          {
            ComputeValuesPhiRef(points_dof2d(k), phi);            
            for (int j = 0; j < nb_dof_loc; j++)
              {
                Real_wp val = DotProd(phi(j), tangente_dof(k));
                Real_wp val_exact = 0;
                if (j == k)
                  val_exact = 1.0;
                
                if (abs(val - val_exact) > 1e4*epsilon_machine)
                  {
                    DISP(val); DISP(val_exact);
                    abort();
                  }
              }
          }
      }
    else
      {
        for (int i = 0; i < order; i++)
          points_dof1d(i) = Real_wp(i+1)/(order+1);

        points_dof2d.Reallocate(nb_dof_loc);
        tangente_dof.Reallocate(nb_dof_loc);
	VectR2 phi;
        for (int i = 0; i < order; i++)
          {
            points_dof2d(i).Init(points_dof1d(i), 0);
            tangente_dof(i).Init(1, 0);
            
            points_dof2d(order+i).Init(1.0-points_dof1d(i), points_dof1d(i));
            tangente_dof(order+i).Init(-1, 1);

            points_dof2d(2*order+i).Init(0, 1.0-points_dof1d(i));
            tangente_dof(2*order+i).Init(0, -1);
          }
        
        int i = 0;
        for (int j = 0; j < order-1; j++)
          for (int k = 0; k < order-1-j; k++)
            {
              points_dof2d(3*order + 2*i).Init(points_dof1d(j), points_dof1d(k));
              points_dof2d(3*order + 2*i+1).Init(points_dof1d(j), points_dof1d(k));
              
              tangente_dof(3*order + 2*i).Init(1, 0);
              tangente_dof(3*order + 2*i+1).Init(0, 1);
              i++;
            }
        
        // then checking phi_i(xi_k) \cdot t_k = delta_{i, k}
        for (int k = 0; k < nb_dof_loc; k++)
          {
            ComputeValuesPhiRef(points_dof2d(k), phi);
            for (int j = 0; j < nb_dof_loc; j++)
              {
                Real_wp val = DotProd(phi(j), tangente_dof(k));
                Real_wp val_exact = 0;
                if (j == k)
                  val_exact = 1.0;
                
                if (abs(val - val_exact) > 1e4*epsilon_machine)
                  {
                    DISP(val); DISP(val_exact);
                    abort();
                  }
              }
          }
      }

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);
    //ConstructLobattoPoints(order, points_lob, points_dof2d, numb);

    this->nb_points_dof_inside = points_dof2d.GetM();
    this->num_dof_points_surf.Reallocate(3);
    for (int n = 0; n < 3; n++)
      {
	this->num_dof_points_surf(n).Reallocate(order);
	for (int i = 0; i < order; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }    
  }
  
  
  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void TriangleHcurlFirstFamily::
  ComputeValuesPhiOrthoRef(const R2& pointloc, VectR2& phi) const
  {
    Real_wp xtilde(0), ytilde(1);    
    if (abs(pointloc(1) - 1.0) > epsilon_machine)
      {
	xtilde = 2.0*pointloc(0)/(1.0-pointloc(1)) - 1.0;
	ytilde = 2.0*pointloc(1) - 1.0;
      }
    
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    
    int num = 0;
    VectReal_wp Li, Pj; phi.Reallocate(nb_dof_loc);
    EvaluateJacobiPolynomial(LegendrePolynom, order, xtilde, Li);
    Real_wp pow_oneMinusytilde = 1.0, oneMinusy = 0.5*(1.0-ytilde), vali;
    for (int i = 0; i < order; i++)
      {
        if (i == 0)
          vali = 1.0;
        else if (i == 1)
          vali = 2.0*pointloc(0) + pointloc(1) -1.0;
        else
          vali = pow_oneMinusytilde*Li(i);
        
	EvaluateJacobiPolynomial(OddJacobiPolynom(i), order-i, ytilde, Pj);
        // part due to P_{r-1}^2
	for (int j = 0; j <= order-1-i; j++)
	  {
            phi(num).Init(Pj(j)*vali, 0);
            phi(num+1).Init(0, Pj(j)*vali);
            num += 2;
          }
        
        // part due to S_r
        phi(num++).Init(vali*Pj(order-1-i)*pointloc(1), -vali*Pj(order-1-i)*pointloc(0));
        
	pow_oneMinusytilde *= oneMinusy;
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) *= InvWeightBasisRr(i);
  }
  
  
  //! Evaluating divergence of nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] curl_phi curl of nearly orthogonal functions on pointloc
  */
  void TriangleHcurlFirstFamily::
  ComputeCurlPhiOrthoRef(const R2& pointloc, VectReal_wp& curl_phi) const
  {
    Real_wp xtilde(0), ytilde(1), dxtilde_dx(0), dxtilde_dy(0);
    if (abs(pointloc(1) - 1.0) > epsilon_machine)
      {
	xtilde = 2.0*pointloc(0)/(1.0-pointloc(1)) - 1.0;
        dxtilde_dx = 2.0/(1.0-pointloc(1));
        dxtilde_dy = dxtilde_dx*pointloc(0)/(1.0-pointloc(1));
        
        ytilde = 2.0*pointloc(1) - 1.0;
      }

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    
    int num = 0;
    VectReal_wp Li, Pj, dLi, dPj;
    curl_phi.Reallocate(nb_dof_loc);
    EvaluateJacobiPolynomial(LegendrePolynom, order, xtilde, Li, dLi);
    Real_wp powm1_oneMinusytilde = 0.0, pow_oneMinusytilde = 1.0,
      oneMinusy = 0.5*(1.0-ytilde), vali;
    
    R2 gradi;
    for (int i = 0; i < order; i++)
      {
        if (i == 0)
          {
            vali = 1.0;
            gradi.Fill(0);
          }
        else if (i == 1)
          {
            vali = 2.0*pointloc(0) + pointloc(1) -1.0;
            gradi.Init(2.0, 1.0);
          }
        else
          {
            vali = pow_oneMinusytilde*Li(i);
            gradi(0) = pow_oneMinusytilde*dLi(i)*dxtilde_dx;
            gradi(1) = -i*powm1_oneMinusytilde*Li(i) + pow_oneMinusytilde*dLi(i)*dxtilde_dy;
          }
        
	EvaluateJacobiPolynomial(OddJacobiPolynom(i), order-i, ytilde, Pj, dPj);
        // part due to P_{r-1}^2
	for (int j = 0; j <= order-1-i; j++)
	  {
            curl_phi(num) = -(2.0*vali*dPj(j) + gradi(1)*Pj(j));
            curl_phi(num+1) = Pj(j)*gradi(0);
            num += 2;
          }
        
        // part due to S_r
        curl_phi(num) = -vali*Pj(order-1-i) - pointloc(0)*Pj(order-1-i)*gradi(0);
        curl_phi(num) -= vali*Pj(order-1-i) + pointloc(1)
          *(2.0*vali*dPj(order-1-i) + gradi(1)*Pj(order-1-i));
        num++;
                
        powm1_oneMinusytilde = pow_oneMinusytilde;
	pow_oneMinusytilde *= oneMinusy;
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      curl_phi(i) *= InvWeightBasisRr(i);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TriangleHcurlFirstFamily::
  ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc);
    if (type_nodal_basis != NODAL_REGULAR)
      {
        VectR2 psi;
        ComputeValuesPhiOrthoRef(point_loc, psi);
        FillZero(res);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            Add(InverseVDM_Edge(i, j), psi(j), res(i));
        
        return;
      }
    
    Real_wp l2 = point_loc(0), l3 = point_loc(1), l1 = 1.0-l2-l3;
    for (int num_dof = 0; num_dof < nb_dof_loc; num_dof++)
      {
	int i = CoordinateDofs(num_dof, 1);
	int j = CoordinateDofs(num_dof, 2);
	int k = CoordinateDofs(num_dof, 3);
	int num_edge = CoordinateDofs(num_dof,0);
	int nbddls = order+1;
	
	if (num_edge==3)
	  {
	    // first edge, opposite to the vertex 3
	    res(num_dof).Init(1.0-l3,l2);
	    for (int l = 0; l < k; l++)
	      res(num_dof) *= (l3*nbddls-l)/(l+1);
	    
	    for (int l = 1; l < j; l++)
	      res(num_dof) *= (l2*nbddls-l)/l;
	    
	    for (int l = 1; l < i; l++)
	      res(num_dof) *= (l1*nbddls-l)/l;
	  }
	else if (num_edge==2)
	  {
	    // third edge, opposite to the vertex 2
	    res(num_dof).Init(-l3,l2-1.0);
	    
	    for (int l = 0; l < j; l++)
	      res(num_dof) *= (l2*nbddls-l)/(l+1);
	    
	    for (int l = 1; l < k; l++)
	      res(num_dof) *= (l3*nbddls-l)/l;
	    
	    for (int l = 1; l < i; l++)
	      res(num_dof) *= (l1*nbddls-l)/l;
	  }
	else if (num_edge==1)
	  {
	    // second edge, opposite to the vertex 1
	    res(num_dof).Init(-l3,l2);
	    for (int l = 0; l < i; l++)
	      res(num_dof) *= (l1*nbddls-l)/(l+1);
	    
	    for (int l = 1; l < j; l++)
	      res(num_dof) *= (l2*nbddls-l)/l;
	    
	    for (int l = 1; l < k; l++)
	      res(num_dof) *= (l3*nbddls-l)/l;
	  }
      }
    
    // DISP(res);
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void TriangleHcurlFirstFamily::ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    if (type_nodal_basis != NODAL_REGULAR)
      {
        VectReal_wp curl_psi;
        ComputeCurlPhiOrthoRef(point_loc, curl_psi);
        FillZero(res);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            res(i) += InverseVDM_Edge(i, j)*curl_psi(j);
        
        return;
      }

    Real_wp l2 = point_loc(0), l3 = point_loc(1), l1 = 1.0-l2-l3;
    for (int num_dof = 0; num_dof < nb_dof_loc; num_dof++)
      {
	int i = CoordinateDofs(num_dof,1);
	int j = CoordinateDofs(num_dof,2);
	int k = CoordinateDofs(num_dof,3);
	R2 grad;
	Real_wp curl(0), val1, val2, val3, val, val_tmp;
	Real_wp l1_mul_n,l2_mul_n,l3_mul_n;
	
	int nbddls = order+1;
	int num_edge = CoordinateDofs(num_dof,0);
	l1_mul_n = l1*nbddls;
	l2_mul_n = l2*nbddls;
	l3_mul_n = l3*nbddls;
    
	// DISP(l1); DISP(l2); DISP(l3);
	// DISP(num_dof); DISP(i); DISP(j); DISP(k); DISP(num_edge);
	
	if (num_edge==3)
	  {
	    // first edge, opposite to the vertex 3
	    val1 = 1.0; val2 = 1.0; val3 = 1.0;
	    for (int l = 0; l < k; l++)
	      val3 *= (l3_mul_n-l)/(l+1);
	    
	    for (int l = 1; l < j; l++)
	      val2 *= (l2_mul_n-l)/l;
	    
	    for (int l = 1; l < i; l++)
	      val1 *= (l1_mul_n-l)/l;
	    
	    // now we do val1 <- val2*val3 ; val2 <- val1*val3 and val3 <- val1*val2  
	    val = val1; val_tmp = val2;
	    val1 = val2*val3;
	    val2 = val*val3;
	    val3 = val*val_tmp;
	    val *= val1;
	    
	    // DISP(val1); DISP(val2); DISP(val3); DISP(val);
	    curl = val*2;
	    val = 0.0;
	    for (int l = 1; l < i; l++)
	      {
		val_tmp = l1_mul_n*val1;
		for (int m = 1; m < i; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l1_mul_n-m)/m;
		    else
		      val_tmp /= m;
		  }
		val += val_tmp;
	      }
	    // DISP(val);
	    curl += val;
	    
	    val = 0.0;
	    for (int l = 1; l < j; l++)
	      {
		val_tmp = l2_mul_n*val2;
		for (int m = 1; m < j; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l2_mul_n-m)/m;
		    else
		      val_tmp /= m;
		  }
		val += val_tmp;
	      }
	    // DISP(val);
	    curl += val; val= 0.0;
	    
	    for (int l = 0; l < k; l++)
	      {
		val_tmp = (l3-1.0)*nbddls*val3;
		for (int m = 0; m < k; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l3_mul_n-m)/(m+1);
		    else
		      val_tmp /= (m+1);
		  }
		val += val_tmp;
	      }
	    // DISP(val);
	    curl += val;
	    
	  }
	else if (num_edge==2)
	  {
	    // third edge, opposite to the vertex 2
	    val1 = 1.0; val2 = 1.0; val3 = 1.0;
	    for (int l = 0; l < j; l++)
	      val2 *= (l2_mul_n-l)/(l+1);
	    
	    for (int l = 1; l < k; l++)
	      val3 *= (l3_mul_n-l)/l;
	    
	    for (int l = 1; l < i; l++)
	      val1 *= (l1_mul_n-l)/l;
	    
	    // now we do val1 <- val2*val3 ; val2 <- val1*val3 and val3 <- val1*val2  
	    val = val1; val_tmp = val2;
	    val1 = val2*val3;
	    val2 = val*val3;
	    val3 = val*val_tmp;
	    val *= val1;
	    
	    curl = val*2;
	    val = 0.0;
	    for (int l = 1; l < i; l++)
	      {
		val_tmp = l1_mul_n*val1;
		for (int m = 1; m < i; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l1_mul_n-m)/m;
		    else
		      val_tmp /= m;
		  }
		val += val_tmp;
	      }
	    curl += val;
	    
	    val = 0.0;
	    for (int l = 0; l < j; l++)
	      {
		val_tmp = (l2-1.0)*nbddls*val2;
		for (int m = 0; m < j; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l2_mul_n-m)/(m+1);
		    else
		      val_tmp /= (m+1);
		  }
		val += val_tmp;
	      }
	    curl += val; val= 0.0;
	    
	    for (int l = 1; l < k; l++)
	      {
		val_tmp = l3_mul_n*val3;
		for (int m = 1; m < k; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l3_mul_n-m)/m;
		    else
		      val_tmp /= m;
		  }
		val += val_tmp;
	      }
	    curl += val;
	  }
	else if (num_edge==1)
	  {
	    // second edge, opposite to the vertex 1
	    val1 = 1.0; val2 = 1.0; val3 = 1.0;
	    for (int l = 0; l < i; l++)
	      val1 *= (l1_mul_n-l)/(l+1);
	    
	    for (int l = 1; l < j; l++)
	      val2 *= (l2_mul_n-l)/l;
	    
	    for (int l = 1; l < k; l++)
	      val3 *= (l3_mul_n-l)/l;
	    
	    // now we do val1 <- val2*val3 ; val2 <- val1*val3 and val3 <- val1*val2  
	    val = val1; val_tmp = val2;
	    val1 = val2*val3;
	    val2 = val*val3;
	    val3 = val*val_tmp;
	    val *= val1;
	    
	    curl = val*2;
	    val = 0.0;
	    for (int l = 0; l < i; l++)
	      {
		val_tmp = (l1-1.0)*nbddls*val1;
		for (int m = 0; m < i; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l1_mul_n-m)/(m+1);
		    else
		      val_tmp /= (m+1);
		  }
		val += val_tmp;
	      }
	    curl += val;
	    
	    val = 0.0;
	    for (int l = 1; l < j; l++)
	      {
		val_tmp = l2_mul_n*val2;
		for (int m = 1; m < j; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l2_mul_n-m)/m;
		    else
		      val_tmp /= m;
		  }
		val += val_tmp;
	      }
	    curl += val; val= 0.0;
	    
	    for (int l = 1; l < k; l++)
	      {
		val_tmp = l3_mul_n*val3;
		for (int m = 1; m < k; m++)
		  {
		    if (m!=l)
		      val_tmp *= (l3_mul_n-m)/m;
		    else
		      val_tmp /= m;
		  }
		val += val_tmp;
	      }
	    curl += val;
	  }
	
	res(num_dof) = curl;
      }

  }
  
  
  //! displays informations about class TriangleHcurlFirstFamily
  ostream& operator <<(ostream& out, const TriangleHcurlFirstFamily& e)
  {
    out << static_cast<const TriangleReference<2>&>(e);
    out << "Triangle of Nedelec's first family " << endl;
    out << "Number of local degrees of freedom " << e.GetNbDof() << endl;
    return out;
  }


  //! displays informations about class TriangleHcurlOptimalFirstFamily
  ostream& operator <<(ostream& out, const TriangleHcurlOptimalFirstFamily& e)
  {
    out<<static_cast<const TriangleReference<2>&>(e);
    out<<"Triangle of Nedelec's first family "<<endl;
    out<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    return out;
  }

}
  
#define MONTJOIE_FILE_TRIANGLE_HCURL_FIRST_FAMILY_CXX
#endif
