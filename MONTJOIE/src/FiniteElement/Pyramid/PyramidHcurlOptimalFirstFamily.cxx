#ifndef MONTJOIE_FILE_PYRAMID_HCURL_OPTIMAL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidHcurlOptimalFirstFamily::PyramidHcurlOptimalFirstFamily() : PyramidReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    type_basis = LOBATTO_POINTS;
    //type_basis = REGULAR_POINTS;

    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  
    
  //! how to number mesh
  void PyramidHcurlOptimalFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 2*order*order);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);

    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
    
    // sign changes on edges
    //FindHcurlSignEdge(*this, nmap);
    
    // rotation of dofs on faces
    FindHcurlLinearCombinationRotation(nmap, *element_tri_surf, *element_quad_surf);
    
  }


  size_t PyramidHcurlOptimalFirstFamily::GetMemorySize() const
  {
    size_t taille = PyramidReference<2>::GetMemorySize();
    taille += InverseVDM_Edge.GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    taille += GL_GX.GetMemorySize();
    taille += NumQuad2D_gauss.GetMemorySize();
    return taille;
  }


  //! construction nearly-orthogonal functions psi_i used to generate the correct space
  void PyramidHcurlOptimalFirstFamily
  ::ConstructOrthogonalBasis(int r, int rgeom, int rquad, int type_quad)
  {  
    if (rquad <= r)
      rquad = r+1;
    
    PyramidReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad, r, r+1,
                                                TriangleQuadrature::QUADRATURE_GAUSS,
                                                Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    nb_dof_loc = order*(order+3)*(2*order+3)/2;
 
    nb_dof_boundaries = 8*order + 4*order*(order-1) + 2*order*order;
    nb_dof_quad = 4*order + 2*order*order;
    nb_dof_tri = 3*order + order*(order-1);
    
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_dof_loc);
    
    // weight for orthonormalization
    InvWeightBasisRr.Reallocate(nb_dof_loc);
    InvWeightBasisRr.Fill(1); int num = 0;
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();

    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        for (int k = 0; k < r-max(i, j); k++)
          {
            Real_wp coef = CoefLegendre(i)*CoefLegendre(j)*CoefEvenJacobi(max(i, j), k);
            InvWeightBasisRr(num) = coef;
            InvWeightBasisRr(num+1) = coef;
            InvWeightBasisRr(num+2) = coef;
            num += 3;
          }
    
    for (int p = 0; p < r; p++)
      InvWeightBasisRr(num++) = CoefLegendre(p)*CoefLegendre(p);

    for (int n = 0; n < r-1; n++)
      for (int m = 0; m <= n; m++)
        {
          Real_wp coef = CoefLegendre(m)*CoefLegendre(n+2);
          InvWeightBasisRr(num) = coef;
          InvWeightBasisRr(num+1) = coef;
          num += 2;
        }

    for (int p = 0; p < r; p++)
      for (int q = 0; q <= r+1; q++)
        {
          Real_wp coef = CoefLegendre(p)*CoefLegendre(q);
          InvWeightBasisRr(num) = coef;
          InvWeightBasisRr(num+1) = coef;
          num += 2;
        }

    this->SetPointsDofND(points_dof3d);    
    // calcul de la matrice de rigidite
    /* Matrix<Real_wp, Symmetric, RowSymPacked> Kh(nb_dof_loc, nb_dof_loc);
    Kh.Fill(0); VectR3 curl_phi;
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeCurlPhiOrthoRef(points3d(k), curl_phi);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = i; j < nb_dof_loc; j++)
            Kh(i, j) += weights3d(k)*DotProd(curl_phi(i), curl_phi(j));
      }
    
      Kh.WriteText(string("Kh")+to_str(order)+".dat"); */
    
  }
    
  //! construction of finite element
  void PyramidHcurlOptimalFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    ConstructOrthogonalBasis(r, rgeom, rquad, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    this->nb_points_dof_inside = this->nb_dof_loc;

    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHcurlElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace(true);

    TriangleHcurlOptimalFirstFamily* Fb_tri = new TriangleHcurlOptimalFirstFamily();
    QuadrangleHcurlOptimalFirstFamily* Fb_quad = new QuadrangleHcurlOptimalFirstFamily(false);
    Fb_tri->ConstructFiniteElement(order);
    Fb_quad->ConstructFiniteElement(order);
    element_tri_surf = Fb_tri;
    element_quad_surf = Fb_quad;

    VectR2 points_dof2d = Fb_tri->PointsDofND();
    this->SetPointsDof2D_tri(points_dof2d);    

    VectR2 points_dof2d_quad = Fb_quad->PointsDofND();
    this->SetPointsDof2D_quad(points_dof2d_quad);
    
    this->num_dof_points_surf.Reallocate(5);
    for (int n = 0; n < 5; n++)
      {
	this->num_dof_points_surf(n).Reallocate(this->GetNbDofBoundary(n));
	for (int i = 0; i < this->GetNbDofBoundary(n); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    
  }
  
  
  //! construction of basis functions
  void PyramidHcurlOptimalFirstFamily::ConstructFunctions()
  {    
    // points where dofs are defined, and tangent on these points
    VectR3 phi, points_inside;
    tangente_dof.Reallocate(nb_dof_loc);
    VectR2 points_tri, points_quad;
    VectReal_wp points_lob, weights_lob;
    Matrix<int> NumQuad(order+2, order+2);
    NumQuad.Fill();
    ComputeGaussLobatto(points_lob, weights_lob, order+1);
    VectReal_wp points_dof1d;
    points_dof1d.Reallocate(order);
    for (int i = 0; i < order; i++)
      points_dof1d(i) = points_lob(i+1);
    
    Fb_geom.ConstructLobattoPoints(order+1, points_lob, points_tri,
				   points_quad, NumQuad, points_inside);
    
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_dof_loc);
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
      for (int j = 1; j <= order; j++)
        {
          Real_wp x = 2.0*points_dof1d(i) - 1.0;
          Real_wp y = 2.0*points_lob(j) - 1.0;
          
          points_dof3d(offset).Init(x, y, 0);
          tangente_dof(offset).Init(2, 0, 0);
          offset++;
        }

    for (int i = 1; i <= order; i++)
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
        points_dof3d(offset+1) = pt;
        points_dof3d(offset+2) = pt;
            
        tangente_dof(offset).Init(1, 0, 0);
        tangente_dof(offset+1).Init(0, 1, 0);
        tangente_dof(offset+2).Init(0, 0, 1);
        offset += 3;
      }

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof3d);

    InvWeightBasisRr.Fill(1.0);
    VectReal_wp Scale(nb_dof_loc); Scale.Fill(0);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
	for (int i = 0; i < nb_dof_loc; i++)
          Scale(i) += DotProd(phi(i), phi(i))*this->WeightsND(k);
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      InvWeightBasisRr(i) = 1.0/sqrt(Scale(i));
    

    Matrix<Real_wp> Mh(nb_dof_loc, nb_dof_loc);
    Mh.Fill(0);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            Mh(i, j) += this->WeightsND(k)*DotProd(phi(i), phi(j));
      }
    
    //Mh.WriteText("MassPyrOpt"+to_str(order)+".dat");

    // DISP(points_dof3d); DISP(tangente_dof);        
    // computation of VDM = psi_j(xi_k) \cdot t_k
    InverseVDM_Edge.Reallocate(nb_dof_loc, nb_dof_loc);
    for (int k = 0; k < nb_dof_loc; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsDofND(k), phi);
        for (int j = 0; j < nb_dof_loc; j++)
          InverseVDM_Edge(j, k) = DotProd(phi(j), tangente_dof(k));
      }
    
    //InverseVDM_Edge.WriteText("VdmPyrOpt"+to_str(order)+".dat");
    GetInverse(InverseVDM_Edge);
    
    Matrix<int> coor;
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(order, NumQuad2D_gauss, coor);
    Globatto<Real_wp> lob_xy;
    lob_xy.ConstructQuadrature(order+1, lob_xy.QUADRATURE_LOBATTO);
    VectReal_wp points_gauss, weights_gauss;
    ComputeGaussLegendre(points_gauss, weights_gauss, order);
    GL_GX.Reallocate(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
	GL_GX(i, j) = lob_xy.EvaluatePhi(i, points_gauss(j));
    
    // then checking phi_i(xi_k) \cdot t_k = delta_{i, k}
    /* for (int k = 0; k < nb_dof_loc; k++)
      {
        ComputeValuesPhiRef(points_dof3d(k), phi);            
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
    */
  }
  

  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
    (\omega, \xi_k) are here Gauss points and not Gauss-Lobatto points
  */  
  template<class Vector1, class Vector2>
  void PyramidHcurlOptimalFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    if (num_loc == 0)
      {
	const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
	Vector1 fx(3*(order+1)*(order+2)), fy(3*(order+2)*(order+2));
	FillZero(fx); FillZero(fy);
	// we express feval on Gauss-Lobatto points
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int node = 3*(i*(order+1) + j);
	      for (int m = 0; m <= order; m++)
		{
		  int p = 3*NumQuad2D_gauss(m, j);
		  fx(node) += GL_GX(i, m)*feval(p);
		  fx(node+1) += GL_GX(i, m)*feval(p+1);
		  fx(node+2) += GL_GX(i, m)*feval(p+2);
		}
	    }
	
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      int node = 3*NumQuad2D(i, j);
	      for (int m = 0; m <= order; m++)
		{
		  int p = 3*(i*(order+1) + m);
		  fy(node) += GL_GX(j, m)*fx(p);
		  fy(node+1) += GL_GX(j, m)*fx(p+1);
		  fy(node+2) += GL_GX(j, m)*fx(p+2);
		}
	    }
	
	// then we use integration with Gauss-Lobatto points
	ComputeIntegralSurfaceRef(fy, res, num_loc);
      }
    else
      ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void PyramidHcurlOptimalFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectR3 psi;
    ComputeValuesPhiOrthoRef(point_loc, psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
        Add(InverseVDM_Edge(i, j), psi(j), res(i));
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] curl curl of basis functions on point_loc
  */
  void PyramidHcurlOptimalFirstFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectR3 curl_psi;
    ComputeCurlPhiOrthoRef(point_loc, curl_psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
        res(i) += InverseVDM_Edge(i, j)*curl_psi(j);
  }


  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void PyramidHcurlOptimalFirstFamily::
  ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    // part due to C_{r-1}^3
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0), c = 2.0*z - 1.0;
    if (abs(1.0 - z) > epsilon_machine)
      {
        a = x/(1.0-z);
        b = y/(1.0-z);
      }

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    VectReal_wp Px, Py, pow_oneMinusZ(r+1);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, b, Py);
    Vector<VectReal_wp> Pz(r+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i));
    
    for (int i = 0; i < r; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val, vloc; int num = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
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
          
          for (int k = 0; k < r-max(i, j); k++)
            {
              vloc = val*Pz(m)(k);
              phi(num).Init(vloc, 0, 0);
              phi(num+1).Init(0, vloc, 0);
              phi(num+2).Init(0, 0, vloc);
              num += 3;
            }
        }
    
    // part due to x^p y^p (1-z)^p (y, x, xy)
    for (int p = 0; p < r; p++)
      {
        val = Px(p)*Py(p)*pow_oneMinusZ(p);
        phi(num).Init(b*val, a*val, a*b*val);
        num++;
      }

    // part due to x^m y^(n+1) (1-z)^(n+1) (y, 0, x y) and x^(n+1) y^m (1-z)^(n+1) (0, x, x y)
    for (int n = 0; n < r-1; n++)
      for (int m = 0; m <= n; m++)
      {
        //val = Px(m)*Py(n+1)*pow_oneMinusZ(n+1);
        //phi(num).Init(b*val, 0, a*b*val);
        
        //val = Px(n+1)*Py(m)*pow_oneMinusZ(n+1);
        //phi(num+1).Init(0, a*val, a*b*val);

        val = Px(m)*Py(n+2)*pow_oneMinusZ(n+1);
        phi(num).Init(val, 0, a*val);
        
        val = Px(n+2)*Py(m)*pow_oneMinusZ(n+1);
        phi(num+1).Init(0, val, b*val);
        num += 2;
      }
    
    // part due to x^p y^q (1-z)^r (1, 0, x) and  x^q y^p (1-z)^r (0, 1, y)
    for (int p = 0; p < r; p++)
      for (int q = 0; q <= r+1; q++)
        {
          val = Px(p)*Py(q)*pow_oneMinusZ(r);
          phi(num).Init(val, 0, a*val);
          
          val = Px(q)*Py(p)*pow_oneMinusZ(r);
          phi(num+1).Init(0, val, b*val);
          num += 2;
        }    
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeightBasisRr(i);
    
  }

  
  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] curl_phi curl of nearly orthogonal functions on pointloc
  */
  void PyramidHcurlOptimalFirstFamily::
  ComputeCurlPhiOrthoRef(const R3& pointloc, VectR3& curl_phi) const
  {
    curl_phi.Reallocate(nb_dof_loc);
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0), da_dx(0), da_dz(0), db_dy(0), db_dz(0), one(1);
    if (abs(1-z) >  1e3*epsilon_machine)
      {
	a = x/(one-z);
	b = y/(one-z);
	da_dx = one/(one-z); da_dz = x/square(one-z);
	db_dy = one/(one-z); db_dz = y/square(one-z);
      }

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    Real_wp c = 2.0*z - 1.0;
    int r = order;
    VectReal_wp Px, Py, dPx, dPy, pow_oneMinusZ(r+1);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, b, Py, dPy);
    Vector<VectReal_wp> Pz(r+1), dPz(r+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= r; i++)      
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i), dPz(i));
    
    for (int i = 0; i < r; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val, dval_dx, dval_dz, dval_dy, valz;
    R3 grad_phi; int num = 0;
    // part due to C_{r-1}^3
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
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
          
          for (int k = 0; k < r-max(i, j); k++)
            {
              valz = Pz(m)(k);
              grad_phi.Init(dval_dx*valz, dval_dy*valz, dval_dz*valz + 2.0*val*dPz(m)(k));
              curl_phi(num).Init(0, grad_phi(2), -grad_phi(1));
              curl_phi(num+1).Init(-grad_phi(2), 0, grad_phi(0));
              curl_phi(num+2).Init(grad_phi(1), -grad_phi(0), 0);
              num += 3;
            }
        }

    // part due to x^p y^p (1-z)^p (y, x, xy)
    R3 vec_u(b, a, a*b), vec_v, rot_u, rot_v;
    for (int p = 0; p < r; p++)
      {
        val = Px(p)*Py(p)*pow_oneMinusZ(p);
        grad_phi(0) = dPx(p)*Py(p)*da_dx*pow_oneMinusZ(p);
        grad_phi(1) = dPy(p)*Px(p)*db_dy*pow_oneMinusZ(p);
        if (p == 0)
          grad_phi(2) = dPx(p)*Py(p)*da_dz*pow_oneMinusZ(p)
            + Px(p)*dPy(p)*db_dz*pow_oneMinusZ(p);
        else
          grad_phi(2) = dPx(p)*Py(p)*da_dz*pow_oneMinusZ(p)
            + Px(p)*dPy(p)*db_dz*pow_oneMinusZ(p) - p*Px(p)*Py(p)*pow_oneMinusZ(p-1); 
        
        //DISP(a); DISP(b); DISP(vec_u); DISP(grad_phi);
        TimesProd(grad_phi, vec_u, curl_phi(num));
        //DISP(num); DISP(curl_phi(num));
        num++;
      }

    // part due to x^m y^(n+1) (1-z)^(n+1) (y, 0, x y) and x^(n+1) y^m (1-z)^(n+1) (0, x, x y)
    //vec_u.Init(b, 0, a*b);
    //rot_u.Init(da_dz, 0, -da_dx);    
    //vec_v.Init(0, a, a*b);
    //rot_v.Init(0, -db_dz, db_dy);  
    vec_u.Init(1.0, 0, a);
    rot_u.Init(0, -da_dx, 0);
    vec_v.Init(0, 1.0, b);
    rot_v.Init(db_dy, 0, 0);
    for (int n = 0; n < r-1; n++)
      for (int m = 0; m <= n; m++)
      {
        /*val = Px(m)*Py(n+1)*pow_oneMinusZ(n+1);
        grad_phi(0) = dPx(m)*Py(n+1)*da_dx*pow_oneMinusZ(n+1);
        grad_phi(1) = dPy(n+1)*Px(m)*db_dy*pow_oneMinusZ(n+1);
        grad_phi(2) = dPx(m)*Py(n+1)*da_dz*pow_oneMinusZ(n+1)
          + Px(m)*dPy(n+1)*db_dz*pow_oneMinusZ(n+1) - (n+1)*Px(m)*Py(n+1)*pow_oneMinusZ(n); 
        
        TimesProd(grad_phi, vec_u, curl_phi(num));
        Add(val, rot_u, curl_phi(num));
        
        val = Px(n+1)*Py(m)*pow_oneMinusZ(n+1);
        grad_phi(0) = dPx(n+1)*Py(m)*da_dx*pow_oneMinusZ(n+1);
        grad_phi(1) = dPy(m)*Px(n+1)*db_dy*pow_oneMinusZ(n+1);
        grad_phi(2) = dPx(n+1)*Py(m)*da_dz*pow_oneMinusZ(n+1)
          + Px(n+1)*dPy(m)*db_dz*pow_oneMinusZ(n+1) - (n+1)*Px(n+1)*Py(m)*pow_oneMinusZ(n); 
        
        TimesProd(grad_phi, vec_v, curl_phi(num+1));
        Add(val, rot_v, curl_phi(num+1)); */
        
        val = Px(m)*Py(n+2)*pow_oneMinusZ(n+1);
        grad_phi(0) = dPx(m)*Py(n+2)*da_dx*pow_oneMinusZ(n+1);
        grad_phi(1) = dPy(n+2)*Px(m)*db_dy*pow_oneMinusZ(n+1);
        grad_phi(2) = dPx(m)*Py(n+2)*da_dz*pow_oneMinusZ(n+1)
          + Px(m)*dPy(n+2)*db_dz*pow_oneMinusZ(n+1) - (n+1)*Px(m)*Py(n+2)*pow_oneMinusZ(n); 
        
        TimesProd(grad_phi, vec_u, curl_phi(num));
        Add(val, rot_u, curl_phi(num));
        
        val = Px(n+2)*Py(m)*pow_oneMinusZ(n+1);
        grad_phi(0) = dPx(n+2)*Py(m)*da_dx*pow_oneMinusZ(n+1);
        grad_phi(1) = dPy(m)*Px(n+2)*db_dy*pow_oneMinusZ(n+1);
        grad_phi(2) = dPx(n+2)*Py(m)*da_dz*pow_oneMinusZ(n+1)
          + Px(n+2)*dPy(m)*db_dz*pow_oneMinusZ(n+1) - (n+1)*Px(n+2)*Py(m)*pow_oneMinusZ(n); 
        
        TimesProd(grad_phi, vec_v, curl_phi(num+1));
        Add(val, rot_v, curl_phi(num+1));
        
        num += 2;
      }
    
    // part due to x^p y^q (1-z)^r (1, 0, x) and  x^q y^p (1-z)^r (0, 1, y)
    vec_u.Init(1.0, 0, a);
    rot_u.Init(0, -da_dx, 0);
    vec_v.Init(0, 1.0, b);
    rot_v.Init(db_dy, 0, 0);
    for (int p = 0; p < r; p++)
      for (int q = 0; q <= r+1; q++)
        {
          val = Px(p)*Py(q)*pow_oneMinusZ(r);
          grad_phi(0) = dPx(p)*Py(q)*da_dx*pow_oneMinusZ(r);
          grad_phi(1) = dPy(q)*Px(p)*db_dy*pow_oneMinusZ(r);
          grad_phi(2) = dPx(p)*Py(q)*da_dz*pow_oneMinusZ(r)
            + Px(p)*dPy(q)*db_dz*pow_oneMinusZ(r) - r*Px(p)*Py(q)*pow_oneMinusZ(r-1); 
          
          TimesProd(grad_phi, vec_u, curl_phi(num));
          Add(val, rot_u, curl_phi(num));

          val = Px(q)*Py(p)*pow_oneMinusZ(r);
          grad_phi(0) = dPx(q)*Py(p)*da_dx*pow_oneMinusZ(r);
          grad_phi(1) = dPy(p)*Px(q)*db_dy*pow_oneMinusZ(r);
          grad_phi(2) = dPx(q)*Py(p)*da_dz*pow_oneMinusZ(r)
            + Px(q)*dPy(p)*db_dz*pow_oneMinusZ(r) - r*Px(q)*Py(p)*pow_oneMinusZ(r-1); 
          
          TimesProd(grad_phi, vec_v, curl_phi(num+1));
          Add(val, rot_v, curl_phi(num+1));
          
          num += 2;
        }    
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < nb_dof_loc; i++)
      curl_phi(i) *= InvWeightBasisRr(i);
    
  }
  

  //! displays details of class PyramidHcurlFirstFamily
  ostream& operator <<(ostream& out, const PyramidHcurlOptimalFirstFamily& e)
  {
    out<<static_cast<const PyramidReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HCURL_OPTIMAL_FIRST_FAMILY_CXX
#endif
