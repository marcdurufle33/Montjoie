#ifndef MONTJOIE_FILE_TRIANGLE_DG_ORTHO_CXX

namespace Montjoie
{
  //! default constructor
  TriangleDgOrtho::TriangleDgOrtho() : TriangleReference<1>()
  {
    this->discontinuous_element = true;
    use_quadrature_for_rh = true;
    use_quadrature_for_sh = true;
    type_quadrature = TriangleQuadrature::QUADRATURE_TENSOR;
  }


  size_t TriangleDgOrtho::GetMemorySize() const
  {
    size_t taille = TriangleReference<1>::GetMemorySize();
    taille += NumFct2D.GetMemorySize() + InvWeightFct.GetMemorySize();
    taille += ChX.GetMemorySize();
    taille += ChY.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += points1d_y.GetMemorySize() + ValGaussX0.GetMemorySize()
      + ValGaussX1.GetMemorySize() + ValGaussY0.GetMemorySize();
    taille += DerivGaussX0.GetMemorySize() + DerivGaussX1.GetMemorySize()
      + DerivGaussY0.GetMemorySize();
    taille += Gy_G.GetMemorySize() + dGy_G.GetMemorySize();
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDy.GetMemorySize();
    taille += gauss_x.GetMemorySize() + gauss_y.GetMemorySize();
    taille += Seldon::GetMemorySize(ProjOperatorOrderY);
    return taille;
  }

  
  //! constructing finite element
  void TriangleDgOrtho::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					       int rsurf, int type_surf)
  {    
    if (type_quad == -1)
      type_quad = TriangleQuadrature::QUADRATURE_TENSOR;

    if (rquad == 0)
      rquad = r;

    if (r != rquad)
      {
        cout << "A quadrature order different from order not implemented for TriangleDgOrtho : " << r << " not equal to " << rquad << endl;
        abort();
      }
    
    TriangleReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    type_quadrature = type_quad;
    
    order = r;
    nb_dof_loc = (r+1)*(r+2)/2;
    VectR2 points_dof2d;
    points_dof2d = this->PointsND();
    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    this->elt_geom.dof_equal_quadrature = true;
    
    this->SetPointsDof1D(this->Points1D());
    this->SetPointsDofND(points_dof2d);
    
    InvWeightFct.Reallocate((r+1)*(r+2)/2);
    NumFct2D.Reallocate(r+1, r+1); NumFct2D.Fill(-1);
    int nb = 0;
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();	
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        {
          NumFct2D(i, j) = nb;
          InvWeightFct(nb) = CoefLegendre(i)*CoefOddJacobi(i, j);
          nb++;
        }
    
    this->ConstructMassMatrix();
    this->ConstructStiffnessMatrix();
    this->Fb_geom.ComputeCoefficientTransformation();    
    
    this->ConstructElementaryMatrix(*this);
  }

  
  void TriangleDgOrtho::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    EdgeHierarchicReference* edge = new EdgeHierarchicReference();
    edge->ConstructFiniteElement(r, r, rquad, type_quad, EdgeHierarchic::LEGENDRE);
    this->element_surface = edge;    
  }
  
  
  //! constructing mass matrix
  void TriangleDgOrtho::ConstructMassMatrix()
  {
    VectReal_wp weights1d, Pn;
    ComputeGaussJacobi(points1d_y, weights1d, order, Real_wp(1), Real_wp(0));
    
    Matrix<Real_wp, General, ArrayRowSparse> ch;
    ch.Reallocate(nb_dof_loc, (order+1)*(order+1));
    Real_wp vloc(0);
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    for (int i = 0; i <= order; i++)
      for (int k = 0; k <= order; k++)
	{
	  EvaluateJacobiPolynomial(OddJacobiPolynom(i), order-i, 2.0*points1d_y(k) - 1.0, Pn); 
	  for (int j = 0; j <= order-i; j++)
	    {
	      int num_dof = NumFct2D(i, j);
	      vloc = Pn(j)*CoefOddJacobi(i, j)*pow(1.0-points1d_y(k), i);
	      ch.AddInteraction(num_dof, i*(order+1) + k, vloc);
	    }	    
	}
    
    Copy(ch, ChY);
    
    ch.Clear();
    ch.Reallocate((order+1)*(order+1), (order+1)*(order+1));
    const VectReal_wp& points1d = this->Points1D();
    for (int i = 0; i <= order; i++)
      for (int k = 0; k <= order; k++)
	{
	  EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*points1d(k) - 1.0, Pn); 
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = i*(order+1) + j;
	      vloc = Pn(i)*CoefLegendre(i);
	      ch.AddInteraction(num_dof, k*(order+1) + j, vloc);
	    }	    
	}
    
    Copy(ch, ChX);
  }
  
  
  //! constructing stiffness matrix
  void TriangleDgOrtho::ConstructStiffnessMatrix()
  {
    const VectReal_wp& points1d = this->Points1D();
    gauss_x.AffectPoints(points1d);
    gauss_y.AffectPoints(points1d_y);
    gauss_x.ComputeGradPhi(1e3*epsilon_machine);
    gauss_y.ComputeGradPhi(1e3*epsilon_machine);
    
    Matrix<Real_wp, General, ArrayRowSparse> rh;
    rh.Reallocate((order+1)*(order+1), 2*(order+1)*(order+1));
    Real_wp vloc(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    vloc = gauss_x.GradPhi(i, k);
	    rh.AddInteraction(i*(order+1)+j, 2*(k*(order+1)+j), vloc);
	    
	    vloc = gauss_y.GradPhi(j, k);
	    rh.AddInteraction(i*(order+1)+j, 2*(i*(order+1)+k)+1, vloc);
	  }
    
    Copy(rh, rh_loc);
    
    ValGaussX0.Reallocate(order+1);
    ValGaussX1.Reallocate(order+1);
    ValGaussY0.Reallocate(order+1);
    DerivGaussX0.Reallocate(order+1);
    DerivGaussX1.Reallocate(order+1);
    DerivGaussY0.Reallocate(order+1);
    int N = (order+1)*(order+1);
    DerivDxtildeDx.Reallocate(N+3*(order+1));
    DerivDxtildeDy.Reallocate(N+3*(order+1));
    Gy_G.Reallocate(order+1, order+1);
    dGy_G.Reallocate(order+1, order+1);
    const VectR2& points2d = this->PointsND();
    for (int i = 0; i <= order; i++)
      {
	ValGaussX0(i) = gauss_x.EvaluatePhi(i, Real_wp(0));
	ValGaussX1(i) = gauss_x.EvaluatePhi(i, Real_wp(1));
	ValGaussY0(i) = gauss_y.EvaluatePhi(i, Real_wp(0));

	DerivGaussX0(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(0));
	DerivGaussX1(i) = gauss_x.EvaluatePhiGrad(i, Real_wp(1));
	DerivGaussY0(i) = gauss_y.EvaluatePhiGrad(i, Real_wp(0));
	for (int j = 0; j <= order; j++)
	  {
	    Gy_G(i, j) = gauss_y.EvaluatePhi(i, points1d(j));
            dGy_G(i, j) = gauss_y.EvaluatePhiGrad(i, points1d(j));
	    int node = i*(order+1) + j;
	    Real_wp x = points2d(node)(0), y = points2d(node)(1);
	    DerivDxtildeDx(node) = 1.0/(1.0-y);
	    DerivDxtildeDy(node) = x/square(1.0-y);
	  }
      }
    
    for (int i = 0; i <= order; i++)
      {
        Real_wp x = points1d(i), y = 0;
        int node = N + i;
        DerivDxtildeDx(node) = 1.0/(1.0-y);
        DerivDxtildeDy(node) = x/square(1.0-y);

        x = 1.0-points1d(i); y = points1d(i);
        node = N + (order+1) + i;
        DerivDxtildeDx(node) = 1.0/(1.0-y);
        DerivDxtildeDy(node) = x/square(1.0-y);

        x = 0; y = 1.0 - points1d(i);
        node = N + 2*(order+1) + i;
        DerivDxtildeDx(node) = 1.0/(1.0-y);
        DerivDxtildeDy(node) = x/square(1.0-y);
      }
  }
  

  //! Integration against gradient of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \varphi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(this->nb_points_quadrature_inside);
    ApplyRhQuadratureGen(Uh, Uquad);
    ApplyChGen(Uquad, Vh);
  }
  
  
  //! computation of gradient of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j grad phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = grad phi_i(xi_j)
  */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(this->nb_points_quadrature_inside);
    ApplyChTransposeGen(Uh, Uquad); 
    ApplyRhQuadratureTransposeGen(Uquad, Vh);
  }
  
  
  //! Integration against gradient of basis functions associated with quadrature points
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \psi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
  */  
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 dUh(2*nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dUh(2*i) = DerivDxtildeDx(i)*Uh(2*i) + DerivDxtildeDy(i)*Uh(2*i+1);
	dUh(2*i+1) = Uh(2*i+1);
      }
    
    Mlt(rh_loc, dUh, Vh);
  }
    
  
  //! Computation of gradient on quadrature points from values on quadrature points
  /*!
    \param[in] Uh values of u on quadrature points
    \param[out] Vh gradient of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Rh* Uh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
  */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);

    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Vh(2*i+1) += Vh(2*i)*DerivDxtildeDy(i);
	Vh(2*i) *= DerivDxtildeDx(i);
      }
  }

  
  //! Integration against basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Ch Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
  */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (type_quadrature == TriangleQuadrature::QUADRATURE_TENSOR)
      {
	Vector1 Ux(this->nb_points_quadrature_inside);
	Mlt(ChX, Uh, Ux);
	Mlt(ChY, Ux, Vh);
      }
    else
      ElementReference<Dimension2, 1>::ApplyChGen(Uh, Vh);
  }
  
  
  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (type_quadrature == TriangleQuadrature::QUADRATURE_TENSOR)
      {
	Vector1 Ux(this->nb_points_quadrature_inside);
	Mlt(SeldonTrans, ChY, Uh, Ux);
	Mlt(SeldonTrans, ChX, Ux, Vh);
      }
    else
      ElementReference<Dimension2, 1>::ApplyChTransposeGen(Uh, Vh);
  }
  
  
  //! computation of gradient of u on quadrature points of an edge
  /*!
    \param[in] num_loc edge number
    \param[in] Uh values of u on quadrature points
    \param[out] Vh values of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
    and \varphi_i are functions associated with quadrature points
  */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho
  ::ApplyShQuadratureTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0)||(r == order))
      {
	switch (num_loc)
	  {
	  case 0 :
	    {
	      int nb = 0;
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  Vh(i) += ValGaussY0(j)*Uh(nb++);
	    }
	    break;
	  case 1 :
	    {
	      typename Vector2::value_type vloc;
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int i = 0; i <= order; i++)
		    vloc += ValGaussX1(i)*Uh(i*(order+1) + j);
		  
		  for (int k = 0; k <= order; k++)
		    Vh(k) += Gy_G(j, k)*vloc;
		}
	    }
	    break;
	  case 2 :
	    {
	      typename Vector2::value_type vloc;
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int i = 0; i <= order; i++)
		    vloc += ValGaussX0(i)*Uh(i*(order+1) + j);
		  
		  for (int k = 0; k <= order; k++)
		    Vh(order-k) += Gy_G(j, k)*vloc;
		}
	    }
	    break;
	  }
      }
    else
      {
	switch (num_loc)
	  {
	  case 0 :
	    {
	      int nb = 0;
	      const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	      typename Vector2::value_type vloc;
	      for (int i = 0; i <= order; i++)
		{
		  vloc = 0;
		  for (int j = 0; j <= order; j++)
		    vloc += ValGaussY0(j)*Uh(nb++);
		  
		  for (int k = 0; k <= r; k++)
		    Vh(k) += Ah(i, k)*vloc;
		}
	    }
	    break;
	  case 1 :
	    {
	      typename Vector2::value_type vloc;
	      const Matrix<Real_wp>& Ah = ProjOperatorOrderY(r);
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int i = 0; i <= order; i++)
		    vloc += ValGaussX1(i)*Uh(i*(order+1) + j);
		  
		  for (int k = 0; k <= r; k++)
		    Vh(k) += Ah(j, k)*vloc;
		}
	    }
	    break;
	  case 2 :
	    {
	      typename Vector2::value_type vloc;
	      const Matrix<Real_wp>& Ah = ProjOperatorOrderY(r);
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int i = 0; i <= order; i++)
		    vloc += ValGaussX0(i)*Uh(i*(order+1) + j);
		  
		  for (int k = 0; k <= r; k++)
		    Vh(r-k) += Ah(j, k)*vloc;
		}
	    }
	    break;
	  }
      }
  }
    
  
  //! integration against quadrature functions on an edge
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc edge number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    \psi_i are functions associated with quadrature points
    (\omega_k, \xi_k) is a quadrature rule on the edge
   */
  template<class T0, class Vector1, class Vector2>
  void TriangleDgOrtho
  ::ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order))
      {
	switch (num_loc)
	  {
	  case 0 :
	    {
	      int nb = 0;
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  Vh(nb++) += alpha*ValGaussY0(j)*Uh(i);
	    }
	    break;
	  case 1 :
	    {
	      typename Vector2::value_type vloc;
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int k = 0; k <= order; k++)
		    vloc += Gy_G(j, k)*Uh(k);
		  
		  vloc *= alpha;
		  for (int i = 0; i <= order; i++)
		    Vh(i*(order+1) + j) += ValGaussX1(i)*vloc;
		  		  
		}
	    }
	    break;
	  case 2 :
	    {
	      typename Vector2::value_type vloc;
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int k = 0; k <= order; k++)
		    vloc += Gy_G(j, k)*Uh(order-k);
		  
		  vloc *= alpha;
		  for (int i = 0; i <= order; i++)
		    Vh(i*(order+1) + j) += ValGaussX0(i)*vloc;		  
		}
	    }
	    break;
	  }
      }
    else
      {
	switch (num_loc)
	  {
	  case 0 :
	    {
	      int nb = 0;
	      typename Vector2::value_type vloc;
	      const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	      for (int i = 0; i <= order; i++)
		{
		  vloc = 0;
		  for (int k = 0; k <= r; k++)
		    vloc += Ah(i, k)*Uh(k);
		  
		  vloc *= alpha;
		  for (int j = 0; j <= order; j++)
		    Vh(nb++) += ValGaussY0(j)*vloc;
		}
	    }
	    break;
	  case 1 :
	    {
	      typename Vector2::value_type vloc;
	      const Matrix<Real_wp>& Ah = ProjOperatorOrderY(r);
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int k = 0; k <= r; k++)
		    vloc += Ah(j, k)*Uh(k);
		  
		  vloc *= alpha;
		  for (int i = 0; i <= order; i++)
		    Vh(i*(order+1) + j) += ValGaussX1(i)*vloc;
		  		  
		}
	    }
	    break;
	  case 2 :
	    {
	      typename Vector2::value_type vloc;
	      const Matrix<Real_wp>& Ah = ProjOperatorOrderY(r);
	      for (int j = 0; j <= order; j++)
		{
		  vloc = 0;
		  for (int k = 0; k <= r; k++)
		    vloc += Ah(j, k)*Uh(r-k);
		  
		  vloc *= alpha;
		  for (int i = 0; i <= order; i++)
		    Vh(i*(order+1) + j) += ValGaussX0(i)*vloc;
		}
	    }
	    break;
	  }
      }
  }

  
  //! computation of gradient of u on quadrature points of an edge
  /*!
    \param[in] num_loc edge number
    \param[in] Uh values of u on quadrature points
    \param[out] Vh gradient of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    V = Sh* U
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
    and \varphi_i functions associated with quadrature points
  */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& Uh,
							   Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    Vector2 Vquad(2*(order+1)); Vquad.Fill(0);
    
    typename Vector2::value_type vx, vy;
    switch (num_loc)
      {
      case 0 :
        {
          int nb = 0;
          for (int i = 0; i <= order; i++)
            {
              vx = 0; vy = 0;
              for (int j = 0; j <= order; j++)
                {
                  vx += ValGaussY0(j)*Uh(nb);
                  vy += DerivGaussY0(j)*Uh(nb);
                  nb++;
                }
              
              Vquad(2*i+1) = vy;
              for (int j = 0; j <= order; j++)
                Vquad(2*j) += gauss_x.GradPhi(i, j)*vx;
            }
        }
        break;
      case 1 :
        {
          int node;
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              for (int i = 0; i <= order; i++)
                {
                  node = i*(order+1) + j;
                  vx += DerivGaussX1(i)*Uh(node);
                  vy += ValGaussX1(i)*Uh(node);
                }
              
              for (int k = 0; k <= order; k++)
                {
                  Vquad(2*k) += Gy_G(j, k)*vx;
                  Vquad(2*k+1) += dGy_G(j, k)*vy;
                }
            }
        }
        break;
      case 2 :
        {
          int node;
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              for (int i = 0; i <= order; i++)
                {
                  node = i*(order+1) + j;
                  vx += DerivGaussX0(i)*Uh(node);
                  vy += ValGaussX0(i)*Uh(node);
                }
              
              for (int k = 0; k <= order; k++)
                {
                  Vquad(2*(order-k)) += Gy_G(j, k)*vx;
                  Vquad(2*(order-k)+1) += dGy_G(j, k)*vy;
                }
            }
        }
        break;
      }
    
    int offset = nb_points_quadrature_inside + (order+1)*num_loc;
    for (int i = 0; i <= order; i++)
      {
        Vquad(2*i+1) += Vquad(2*i)*DerivDxtildeDy(offset+i);
        Vquad(2*i) *= DerivDxtildeDx(offset+i);
      }
    
    if ( (r == 0) || (r == order))
      Copy(Vquad, Vh);
    else
      {
        const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
        for (int i = 0; i <= order; i++)
          for (int k = 0; k <= r; k++)
            {
              Vh(2*k) += Ah(i, k)*Vquad(2*i);
              Vh(2*k+1) += Ah(i, k)*Vquad(2*i+1);
            }
      }
  }
    
  
  //! integration against gradient of quadrature functions on an edge
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc edge number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    \psi_i are functions associated with quadrature points
    (\omega_k, \xi_k) is a quadrature rule on the edge
   */
  template<class T0, class Vector1, class Vector2>
  void TriangleDgOrtho::ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh,
						  Vector2& Vh, int r) const
  {    
    Vector2 Vquad(2*(order+1));
    
    if ( (r == 0) || (r == order))
      Copy(Uh, Vquad);
    else
      {
        const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
        Vquad.Fill(0);
        for (int i = 0; i <= order; i++)
          for (int k = 0; k <= r; k++)
            {
              Vquad(2*i) += Ah(i, k)*Uh(2*k);
              Vquad(2*i+1) += Ah(i, k)*Uh(2*k+1);
            }
      }

    int offset = nb_points_quadrature_inside + (order+1)*num_loc;
    for (int i = 0; i <= order; i++)
      Vquad(2*i) = Vquad(2*i)*DerivDxtildeDx(offset+i) + Vquad(2*i+1)*DerivDxtildeDy(offset+i);
    
    Mlt(alpha, Vquad);
    typename Vector2::value_type vx, vy;
    switch (num_loc)
      {
      case 0 :
        {
          int nb = 0;
          for (int i = 0; i <= order; i++)
            {
              vx = 0; vy = 0;
              vy = Vquad(2*i+1);
              for (int j = 0; j <= order; j++)
                vx += gauss_x.GradPhi(i, j)*Vquad(2*j);
              
              for (int j = 0; j <= order; j++)
                {
                  Vh(nb) += ValGaussY0(j)*vx;
                  Vh(nb) += DerivGaussY0(j)*vy;
                  nb++;
                }
            }
        }
        break;
      case 1 :
        {
          int node;
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              for (int k = 0; k <= order; k++)
                {
                  vx += Gy_G(j, k)*Vquad(2*k);
                  vy += dGy_G(j, k)*Vquad(2*k+1);
                }
              
              for (int i = 0; i <= order; i++)
                {
                  node = i*(order+1) + j;
                  Vh(node) += DerivGaussX1(i)*vx;
                  Vh(node) += ValGaussX1(i)*vy;
                }
            }
        }
        break;
      case 2 :
        {
          int node;
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              for (int k = 0; k <= order; k++)
                {
                  vx += Gy_G(j, k)*Vquad(2*(order-k)) ;
                  vy += dGy_G(j, k)*Vquad(2*(order-k)+1);
                }
              
              for (int i = 0; i <= order; i++)
                {
                  node = i*(order+1) + j;
                  Vh(node) += DerivGaussX0(i)*vx;
                  Vh(node) += ValGaussX0(i)*vy;
                }
            }
        }
        break;
      }    
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] contrib dof components
   */
  template<class Vector1, class Vector2>
  void TriangleDgOrtho::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights2d = this->WeightsND();
    for (int i = 0; i < this->nb_points_dof_inside; i++)
      feval_weight(i) *= weights2d(i);
    
    ApplyCh(feval_weight, contrib);
  }


  //! computation of projection between finite element of different orders    
  void TriangleDgOrtho::
  ComputeInterpolationProjectorOrder(const IVect& order_elt,
				     const ElementReference<Dimension1, 1>& Fb,
				     const Vector<VectReal_wp>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi(this->order+1);
    ProjOperatorOrder.Reallocate(rmax+1);
    ProjOperatorOrderY.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorOrder(r).Reallocate(this->order+1, Pts(i).GetM());
	ProjOperatorOrderY(r).Reallocate(this->order+1, Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int k = 0; k <= this->order; k++)
	      {
		ProjOperatorOrder(r)(k, j) = gauss_x.EvaluatePhi(k, Pts(i)(j));
		ProjOperatorOrderY(r)(k, j) = gauss_y.EvaluatePhi(k, Pts(i)(j));
	      }
	  }
      }
  }
    
}

#define MONTJOIE_FILE_TRIANGLE_DG_ORTHO_CXX
#endif
