#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_LEGENDRE_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronDgLegendre::HexahedronDgLegendre() : HexahedronReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    this->discontinuous_element = true;
    
    mass_lumping = false;
    diagonal_mass = false;
    
#ifdef MONTJOIE_FAST_HEXA 
    use_quadrature_free_sh = true;
#endif
    
  }


  size_t HexahedronDgLegendre::GetMemorySize() const
  {
    size_t taille = HexahedronReference<1>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    taille += InvWeightPolynomial.GetMemorySize();
    taille += NumOrtho3D.GetMemorySize();
    taille += rh_const.GetMemorySize() + rh_constX.GetMemorySize() +
      rh_constY.GetMemorySize() + rh_constZ.GetMemorySize();
    taille += sh1_loc.GetMemorySize() + sh2_loc.GetMemorySize() + sh3_loc.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize();
    taille += ValLi.GetMemorySize() + ValDLi.GetMemorySize() +
      Li_Geom.GetMemorySize() + dLi_Geom.GetMemorySize();
    taille += CoefDLi.GetMemorySize() + stiff1D.GetMemorySize();
    taille += Seldon::GetMemorySize(NumProjOperator);

    taille += Seldon::GetMemorySize(ProjOperatorDerivQuadOrder);
    return taille;
  }

  
  //! construction of finite element
  void HexahedronDgLegendre::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						    int rsurf_tri, int rsurf_quad,
						    int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    HexahedronReference<1>::
      ConstructFiniteElement(r, rgeom, rquad, Globatto<Real_wp>::QUADRATURE_GAUSS);

    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    
    ch1_loc.Clear(); ch2_loc.Clear(); ch3_loc.Clear();
    
    // we get Legendre polynomials
    GetJacobiPolynomial(LegendrePolynom, r, Real_wp(0), Real_wp(0));
    
    // then coefficient for orthonormalization
    CoefLegendre.Reallocate(r+1);
    CoefLegendre.Fill(0);
    
    VectReal_wp xi, omega, Pn;
    ComputeGaussLegendre(xi, omega, order);
    for (int q = 0; q <= order; q++)
      {
	EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi(q) - 1.0, Pn);
	for (int i = 0; i <= order; i++)
	  CoefLegendre(i) += omega(q)*Pn(i)*Pn(i);
      }

    for (int i = 0; i <= order; i++)
      CoefLegendre(i) = 1.0/sqrt(CoefLegendre(i));

    nb_dof_loc = (order+1)*(order+2)*(order+3)/6;
    nb_dof_boundaries = 0;
    
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points_dof3d(i) = this->PointsND(i);
    
    this->SetPointsDofND(points_dof3d);
    
    this->nb_points_dof_inside = nb_points_quadrature_inside;

    InvWeightPolynomial.Reallocate(nb_dof_loc);
    int node = 0;
    NumOrtho3D.Reallocate(order+1, order+1, order+1);
    NumOrtho3D.Fill(-1);
    CoordinateDofs.Reallocate(nb_dof_loc, 3);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++)
          {
            NumOrtho3D(i, j, k) = node;
	    CoordinateDofs(node, 0) = i;
	    CoordinateDofs(node, 1) = j;
	    CoordinateDofs(node, 2) = k;
            InvWeightPolynomial(node) = CoefLegendre(i)*CoefLegendre(j)*CoefLegendre(k);
            node++;
          }
    
    // decomposition of dLi/dx on functions Li
    CoefDLi.Reallocate(order+1, order+1);
    VectReal_wp dPn; CoefDLi.Fill(0);
    for (int q = 0; q <= order; q++)
      {
	EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi(q) - 1.0, Pn, dPn);
	for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            CoefDLi(i, j) += 2.0*omega(q)*dPn(i)*Pn(j)*CoefLegendre(i)*CoefLegendre(j);
      }    
    
    stiff1D.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          stiff1D(i, j) = 0; 
          for (int k = 0; k <= order; k++)
            stiff1D(i, j) += CoefDLi(i, k)*CoefDLi(j, k);
        }
    
    // computation of stiffness matrix for affine elements
    Matrix<Real_wp, General, ArrayRowSparse> Rh_const(nb_dof_loc, 3*nb_dof_loc);
    Matrix<Real_wp, General, ArrayRowSparse> Rh_constX(nb_dof_loc, 3*nb_dof_loc);
    Matrix<Real_wp, General, ArrayRowSparse> Rh_constY(nb_dof_loc, 3*nb_dof_loc);
    Matrix<Real_wp, General, ArrayRowSparse> Rh_constZ(nb_dof_loc, 3*nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++)
	  {
	    node = NumOrtho3D(i, j, k);
	    for (int m = 0; m <= order; m++)
	      {
		if (NumOrtho3D(m, j, k) >= 0)
		  if (abs(CoefDLi(m, i)) > 10*epsilon_machine)
                    {
                      Rh_const.AddInteraction(NumOrtho3D(m, j, k), 3*node, CoefDLi(m, i));
                      Rh_constX.AddInteraction(NumOrtho3D(m, j, k), 3*node, CoefDLi(m, i));
                    }
                   		
		if (NumOrtho3D(i, m, k) >= 0)
		  if (abs(CoefDLi(m, j)) > 10*epsilon_machine)
                    {
                      Rh_const.AddInteraction(NumOrtho3D(i, m, k), 3*node+1, CoefDLi(m, j));
                      Rh_constY.AddInteraction(NumOrtho3D(i, m, k), 3*node+1, CoefDLi(m, j));
                    }
                   		
		if (NumOrtho3D(i, j, m) >= 0)
		  if (abs(CoefDLi(m, k)) > 10*epsilon_machine)
                    {
                      Rh_const.AddInteraction(NumOrtho3D(i, j, m), 3*node+2, CoefDLi(m, k));
                      Rh_constZ.AddInteraction(NumOrtho3D(i, j, m), 3*node+2, CoefDLi(m, k));
                    }
                
              }
	  }
    
    Copy(Rh_const, rh_const);
    Copy(Rh_constX, rh_constX);
    Copy(Rh_constY, rh_constY);
    Copy(Rh_constZ, rh_constZ);
    
    // computing Legendre polynomials at -1 and 1
    VectReal_wp Pn0, Pn1;
    EvaluateJacobiPolynomial(LegendrePolynom, order, Real_wp(-1), Pn0);
    EvaluateJacobiPolynomial(LegendrePolynom, order, Real_wp(1), Pn1);
    for (int i = 0; i <= order; i++)
      {
	Pn0(i) *= CoefLegendre(i);
	Pn1(i) *= CoefLegendre(i);
      }
    
    // computation of matrices for flux terms
    Matrix<Real_wp, General, ArrayRowSparse> Sh1, Sh2, Sh3;
    Matrix<int> num2D(order+1, order+1); num2D.Fill(-1);
    node = 0; int Ntri = (order+1)*(order+2)/2;
    Sh1.Reallocate(nb_dof_loc, 6*Ntri);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        {
	  num2D(i, j) = node;
	  for (int k = 0; k <= order-i-j; k++)
	    {
	      Sh1.AddInteraction(NumOrtho3D(k, i, j), node, Pn0(k));
	      Sh1.AddInteraction(NumOrtho3D(i, k, j), Ntri + node, Pn0(k));
	      Sh1.AddInteraction(NumOrtho3D(i, j, k), 2*Ntri + node, Pn0(k));
	      Sh1.AddInteraction(NumOrtho3D(i, j, k), 3*Ntri + node, Pn1(k));
	      Sh1.AddInteraction(NumOrtho3D(i, k, j), 4*Ntri + node, Pn1(k));
	      Sh1.AddInteraction(NumOrtho3D(k, i, j), 5*Ntri + node, Pn1(k));
	    }
	  
	  node++;
	}
    
    ValLi.Reallocate(order+1, order+1);
    ValDLi.Reallocate(order+1, order+1);
    for (int k = 0; k <= order; k++)
      {
	EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi(k)-1.0, Pn, dPn);
	for (int i = 0; i <= order; i++)
	  {
            ValLi(i, k) = Pn(i)*CoefLegendre(i);
            ValDLi(i, k) = 2.0*dPn(i)*CoefLegendre(i);
          }
      }
    
    int Nquad = (order+1)*(order+1);
    Sh2.Reallocate(6*Ntri, 6*Nquad);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        {
	  node = num2D(i, j);
	  for (int k = 0; k <= order; k++)
	    {
	      int n = NumQuad2D(i, k);
	      Sh2.AddInteraction(node, n, ValLi(j, k));
	      Sh2.AddInteraction(Ntri+node, Nquad+n, ValLi(j, k));
	      Sh2.AddInteraction(2*Ntri+node, 2*Nquad+n, ValLi(j, k));
	      Sh2.AddInteraction(3*Ntri+node, 3*Nquad+n, ValLi(j, k));
	      Sh2.AddInteraction(4*Ntri+node, 4*Nquad+n, ValLi(j, k));
	      Sh2.AddInteraction(5*Ntri+node, 5*Nquad+n, ValLi(j, k));
	    }
	}

    Sh3.Reallocate(6*Nquad, 6*Nquad);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
	  node = NumQuad2D(i, j);
	  for (int k = 0; k <= order; k++)
	    {
	      int n = NumQuad2D(k, j);
	      Sh3.AddInteraction(node, n, ValLi(i, k));
	      Sh3.AddInteraction(Nquad+node, Nquad+n, ValLi(i, k));
	      Sh3.AddInteraction(2*Nquad+node, 2*Nquad+n, ValLi(i, k));
	      Sh3.AddInteraction(3*Nquad+node, 3*Nquad+n, ValLi(i, k));
	      Sh3.AddInteraction(4*Nquad+node, 4*Nquad+n, ValLi(i, k));
	      Sh3.AddInteraction(5*Nquad+node, 5*Nquad+n, ValLi(i, k));
	    }
	}
    
    Copy(Sh1, sh1_loc);
    Copy(Sh2, sh2_loc);
    Copy(Sh3, sh3_loc);
    
    Array3D<int> NumPrism3D(order+1, order+1, order+1);
    NumPrism3D.Fill(-1); node = 0;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  NumPrism3D(i, j, k) = node++;


    int nb_points_prism = (order+1)*(order+1)*(order+2)/2;
    Matrix<Real_wp, General, ArrayRowSparse> ch;
    ch.Reallocate(nb_dof_loc, nb_points_prism);
    
    Matrix<Real_wp> Pv(order+1, order+1); VectReal_wp Pm;
    Pv.Fill(0);
    const VectReal_wp& points1d = this->Points1D();
    for (int k = 0; k <= order; k++)
      {
	Real_wp xi_k = 2.0*points1d(k) - 1.0;
	EvaluateJacobiPolynomial(LegendrePolynom, order, xi_k, Pm);
	for (int j = 0; j < Pm.GetM(); j++)
	  Pv(j, k) = Pm(j)*CoefLegendre(j);
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	for (int k = 0; k <= order; k++)
	  {
            node = NumPrism3D(i, j, k);
	    for (int m = 0; m <= order-i-j; m++)
	      ch.AddInteraction(NumOrtho3D(i, j, m), node, Pv(m, k));
          }
    
    Copy(ch, ch1_loc);
    
    ch.Clear(); ch.Reallocate(nb_points_prism, nb_points_quadrature_inside);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
            node = NumQuad3D(i, j, k);
	    for (int m = 0; m <= order-i; m++)
	      ch.AddInteraction(NumPrism3D(i, m, k), node, Pv(m, j));
          }
    
    Copy(ch, ch2_loc);
    
    ch.Clear(); ch.Reallocate(nb_points_quadrature_inside, nb_points_quadrature_inside);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
            node = NumQuad3D(i, j, k);
	    for (int m = 0; m <= order; m++)
	      ch.AddInteraction(NumQuad3D(m, j, k), node, Pv(m, i));
          }
    
    Copy(ch, ch3_loc);
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    int order_geom = this->GetGeometryOrder();
    Li_Geom.Reallocate(order+1, order_geom+1);
    dLi_Geom.Reallocate(order+1, order_geom+1);
    for (int j = 0; j <= order_geom; j++)
      {
        EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*this->PointsNodal1D(j)-1.0, Pn, dPn);
        for (int i = 0; i <= order; i++)
          {
            Li_Geom(i, j) = Pn(i)*CoefLegendre(i);
            dLi_Geom(i, j) = 2.0*dPn(i)*CoefLegendre(i);
          }
      }

    QuadrangleDgGauss* function_basis_quad = new QuadrangleDgGauss();
    function_basis_quad->ConstructFiniteElement(r, rgeom, rquad, type_quad);
    element_quad_surf = function_basis_quad;
  }
  
  
  //! computation of U on nodal points
  /*!
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
   */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::
  ComputeNodalValuesGen(const Vector1& Uloc, Vector2& Uloc_node) const
  {
    int r = order, rg = this->GetGeometryOrder();
    Vector1 Ux( (r+1)*(r+1)*(rg+1));
    Vector1 Uy( (r+1)*(rg+1)*(rg+1));
    Ux.Fill(0); Uy.Fill(0); Uloc_node.Fill(0);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= r-i-j; k++)
          {
            int p = NumOrtho3D(i, j, k);
            for (int m = 0; m <= rg; m++)
              Ux( (rg+1)*(i*(r+1) + j) + m) += Li_Geom(k, m)*Uloc(p);
          }

    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= rg; k++)
          {
            int p = (rg+1)*(i*(r+1) + j) + k;
            for (int m = 0; m <= rg; m++)
              Uy( (rg+1)*(i*(rg+1) + m) + k) += Li_Geom(j, m)*Ux(p);
          }

    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= rg; j++)
        for (int k = 0; k <= rg; k++)
          {
            int p = (rg+1)*(i*(rg+1) + j) + k;
            for (int m = 0; m <= rg; m++)
              Uloc_node(this->GetNumNodes3D(m, j, k)) += Li_Geom(i, m)*Uy(p);
          }
  }
  
    
  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
   */  
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  //! computation of u on nodal points of a face
  /*!
    \param[in] feval dof components of u
    \param[out] res values of u on nodal points of a face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::
  ComputeValueBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    int r = order, rg = this->GetGeometryOrder();
    Vector<typename Vector2::value_type> Ux((r+1)*(r+1)), Uy((rg+1)*(r+1));
    Ux.Fill(0); Uy.Fill(0); res.Fill(0);
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0 :
      case 5:
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*feval(NumOrtho3D(k, i, j));
        
        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*feval(NumOrtho3D(i, k, j));
        
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*feval(NumOrtho3D(i, j, k));
        
        break;
      }
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
        for (int k = 0; k <= rg; k++)
          Uy(i*(rg+1) + k) += Li_Geom(j, k)*Ux(NumQuad2D(i, j));
    
    const Matrix<int>& NumNodes2D_quad = this->Fb_geom.GetNumNodes2D_quad();
    for (int i = 0; i <= rg; i++)
      for (int j = 0; j <= rg; j++)
        for (int k = 0; k <= r; k++)
          res(NumNodes2D_quad(i, j)) += Li_Geom(k, i)*Uy(k*(rg+1) + j);
  }

  
  //! computation of grad u on nodal points of a face
  /*!
    \param[in] feval dof components of u
    \param[out] res gradient of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */  
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::
  ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    int r = order, rg = this->GetGeometryOrder();
    Vector<typename Vector2::value_type> Ux((r+1)*(r+1)), dUx((r+1)*(r+1));
    Vector<typename Vector2::value_type> Uy(3*(rg+1)*(r+1));
    Ux.Fill(0); dUx.Fill(0); Uy.Fill(0); res.Fill(0);
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    const Matrix<int>& NumNodes2D_quad = this->Fb_geom.GetNumNodes2D_quad();    
    switch (num_loc)
      {
      case 0 :
      case 5:
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*feval(NumOrtho3D(k, i, j));
                dUx(NumQuad2D(i, j)) += dLi_Geom(k, ext)*feval(NumOrtho3D(k, i, j));
              }
        
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= rg; k++)
              {
                int p = 3*( i*(rg+1) + k);
                int q = NumQuad2D(i, j);
                Uy(p) += Li_Geom(j, k)*dUx(q);
                Uy(p+1) += Li_Geom(j, k)*Ux(q);
                Uy(p+2) += dLi_Geom(j, k)*Ux(q);
              }
        
        for (int i = 0; i <= rg; i++)
          for (int j = 0; j <= rg; j++)
            for (int k = 0; k <= r; k++)
              {
                int p = 3*NumNodes2D_quad(i, j);
                int q = 3*( k*(rg+1) + j);
                res(p) += Li_Geom(k, i)*Uy(q);
                res(p+1) += dLi_Geom(k, i)*Uy(q+1);
                res(p+2) += Li_Geom(k, i)*Uy(q+2);
              }

        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*feval(NumOrtho3D(i, k, j));
                dUx(NumQuad2D(i, j)) += dLi_Geom(k, ext)*feval(NumOrtho3D(i, k, j));
              }
        
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= rg; k++)
              {
                int p = 3*( i*(rg+1) + k);
                int q = NumQuad2D(i, j);
                Uy(p) += Li_Geom(j, k)*Ux(q);
                Uy(p+1) += Li_Geom(j, k)*dUx(q);
                Uy(p+2) += dLi_Geom(j, k)*Ux(q);
              }
        
        for (int i = 0; i <= rg; i++)
          for (int j = 0; j <= rg; j++)
            for (int k = 0; k <= r; k++)
              {
                int p = 3*NumNodes2D_quad(i, j);
                int q = 3*( k*(rg+1) + j);
                res(p) += dLi_Geom(k, i)*Uy(q);
                res(p+1) += Li_Geom(k, i)*Uy(q+1);
                res(p+2) += Li_Geom(k, i)*Uy(q+2);
              }
        
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*feval(NumOrtho3D(i, j, k));
                dUx(NumQuad2D(i, j)) += dLi_Geom(k, ext)*feval(NumOrtho3D(i, j, k));
              }
        
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= rg; k++)
              {
                int p = 3*( i*(rg+1) + k);
                int q = NumQuad2D(i, j);
                Uy(p) += Li_Geom(j, k)*Ux(q);
                Uy(p+1) += dLi_Geom(j, k)*Ux(q);
                Uy(p+2) += Li_Geom(j, k)*dUx(q);
              }
        
        for (int i = 0; i <= rg; i++)
          for (int j = 0; j <= rg; j++)
            for (int k = 0; k <= r; k++)
              {
                int p = 3*NumNodes2D_quad(i, j);
                int q = 3*( k*(rg+1) + j);
                res(p) += dLi_Geom(k, i)*Uy(q);
                res(p+1) += Li_Geom(k, i)*Uy(q+1);
                res(p+2) += Li_Geom(k, i)*Uy(q+2);
              }
        break;
      }
    
    

  }
  
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      feval_weight(i) *= weights3d(i);
    
    ApplyChGen(feval_weight, contrib);
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
  void HexahedronDgLegendre::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector<typename Vector2::value_type> Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
    Mlt(ch3_loc, Uh, Ux);
    Mlt(ch2_loc, Ux, Uy);
    Mlt(ch1_loc, Uy, Vh);
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
  void HexahedronDgLegendre::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector<typename Vector2::value_type> Ux(ch1_loc.GetN()), Uy(ch2_loc.GetN());
    Mlt(SeldonTrans, ch1_loc, Uh, Ux);
    Mlt(SeldonTrans, ch2_loc, Ux, Uy);
    Mlt(SeldonTrans, ch3_loc, Uy, Vh);
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
  void HexahedronDgLegendre::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Vquad(nb_points_quadrature_inside);
    typename Vector2::value_type vloc;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int p = NumQuad3D(i, j, k);
            vloc = 0;
            for (int m = 0; m <= order; m++)
              {
                vloc += lob_quad.GradPhi(i, m)*Uh(3*NumQuad3D(m, j, k));
                vloc += lob_quad.GradPhi(j, m)*Uh(3*NumQuad3D(i, m, k)+1);
                vloc += lob_quad.GradPhi(k, m)*Uh(3*NumQuad3D(i, j, m)+2);
              }
            
            Vquad(p) = vloc;
          }
    
    ApplyCh(Vquad, Vh);
  }
  

  //! Integration against derivatives of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vx Vx_i = \int_K f d/dx(\varphi_i)  dx
    \param[out] Vy Vy_i = \int_K f d/dy(\varphi_i)  dx
    \param[out] Vz Vz_i = \int_K f d/dz(\varphi_i)  dx
    This operation is equivalent to a matrix vector product
    Vx = Rh^1 Uh, Vy = Rh^2 Uh, Vz = Rh^3 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
    where (Rh)^3_{i,j} = d/dz(\varphi_i(\xi_j) )
   */  
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre
  ::ApplyRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const
  {
    Vector2 Vx_quad(nb_points_quadrature_inside);
    Vector2 Vy_quad(nb_points_quadrature_inside);
    Vector2 Vz_quad(nb_points_quadrature_inside);
    typename Vector2::value_type vloc_x, vloc_y, vloc_z;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int p = NumQuad3D(i, j, k);
            vloc_x = 0; vloc_y = 0; vloc_z = 0;
            for (int m = 0; m <= order; m++)
              {
                vloc_x += lob_quad.GradPhi(i, m)*Uh(3*NumQuad3D(m, j, k));
                vloc_y += lob_quad.GradPhi(j, m)*Uh(3*NumQuad3D(i, m, k)+1);
                vloc_z += lob_quad.GradPhi(k, m)*Uh(3*NumQuad3D(i, j, m)+2);
              }
            
            Vx_quad(p) = vloc_x;
            Vy_quad(p) = vloc_y;
            Vz_quad(p) = vloc_z;
          }
    
    ApplyCh(Vx_quad, Vx);
    ApplyCh(Vy_quad, Vy);
    ApplyCh(Vz_quad, Vz);
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
  void HexahedronDgLegendre::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Uquad(nb_points_quadrature_inside);
    ApplyChTranspose(Uh, Uquad);
    typename Vector2::value_type vloc;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int p = 3*NumQuad3D(i, j, k);
            vloc = 0;
            for (int m = 0; m <= order; m++)
              vloc += lob_quad.GradPhi(m, i)*Uquad(NumQuad3D(m, j, k));
            
            Vh(p) = vloc;
            vloc = 0;
            for (int m = 0; m <= order; m++)
              vloc += lob_quad.GradPhi(m, j)*Uquad(NumQuad3D(i, m, k));
            
            Vh(p+1) = vloc;
            vloc = 0;
            for (int m = 0; m <= order; m++)
              vloc += lob_quad.GradPhi(m, k)*Uquad(NumQuad3D(i, j, m));
            
            Vh(p+2) = vloc;
          }    
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Uh = Rh Vh
   */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::ApplyConstantRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    Mlt(rh_const, Vh, Uh);
  }


  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    (Rh^3)_{i, j} = \int_K \varphi_j d/dz(\varphi_i) dx
    then Ux = Rh^1 Uh, Uy = Rh^2 Uh, Uz = Rh^3 Uh
   */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre
  ::ApplyConstantRhSplitGen(const Vector1& Vh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
    Mlt(rh_constX, Vh, Ux);
    Mlt(rh_constY, Vh, Uy);
    Mlt(rh_constZ, Vh, Uz);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_const, Uh, Vh);
  }

  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the face
    \param[in] r2 order of quadrature rule of the face
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r2) const
  {
    int r = order;
    Vector<typename Vector2::value_type> Ux((r+1)*(r+1));
    Ux.Fill(0); Vh.Fill(0);
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0 :
      case 5:
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*Uh(NumOrtho3D(k, i, j));
        
        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*Uh(NumOrtho3D(i, k, j));
        
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*Uh(NumOrtho3D(i, j, k));
        
        break;
      }
    
    if ((r2 == 0) || (r2 == order))
      {
	typename Vector1::value_type vloc;
	for (int i = 0; i <= r; i++)
	  for (int k = 0; k <= r; k++)
	    {
	      vloc = 0;
	      for (int j = 0; j <= r-i; j++)
		vloc += ValLi(j, k)*Ux(NumQuad2D(i, j));
	      
	      for (int j = 0; j <= r; j++)
		Vh(NumQuad2D(j, k)) += ValLi(i, j)*vloc;
	    }
      }
    else
      {
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r2);
	const Matrix<int>& Num = NumProjOperator(r2);
	for (int i = 0; i <= r; i++)
	  for (int k = 0; k <= r2; k++)
	    {
	      vloc = 0;
	      for (int j = 0; j <= r-i; j++)
		vloc += Ah(j, k)*Ux(NumQuad2D(i, j));
	      
	      for (int j = 0; j <= r2; j++)
		Vh(Num(j, k)) += Ah(i, j)*vloc;
	    }
      }
  }
  
  
  //! integration against basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r2 order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronDgLegendre
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r2) const
  {
    int r = order;
    Vector<typename Vector2::value_type> Ux((r+1)*(r+1));
    Ux.Fill(0);
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((r2 == 0) || (r2 == order))
      {
	typename Vector1::value_type vloc;
	for (int j = 0; j <= r; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      vloc = 0;
	      for (int i = 0; i <= r; i++)
		vloc += ValLi(k, i)*Uh(NumQuad2D(i, j));
	      
	      for (int i = 0; i <= r-k; i++)
		Ux(NumQuad2D(k, i)) += ValLi(i, j)*vloc;
	    }
      }
    else
      {
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r2);
	const Matrix<int>& Num = NumProjOperator(r2);
	for (int j = 0; j <= r2; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      vloc = 0;
	      for (int i = 0; i <= r2; i++)
		vloc += Ah(k, i)*Uh(Num(i, j));
	      
	      for (int i = 0; i <= r-k; i++)
		Ux(NumQuad2D(k, i)) += Ah(i, j)*vloc;
	    }
      }
    
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    switch (num_loc)
      {
      case 0 :
      case 5:
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Vh(NumOrtho3D(k, i, j)) += alpha*Li_Geom(k, ext)*Ux(NumQuad2D(i, j));
        
        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Vh(NumOrtho3D(i, k, j)) += alpha*Li_Geom(k, ext)*Ux(NumQuad2D(i, j));
        
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              Vh(NumOrtho3D(i, j, k)) += alpha*Li_Geom(k, ext)*Ux(NumQuad2D(i, j));
        
        break;
      }
  }

  
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] n boundary number
    \param[in] res components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r2 order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronDgLegendre
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r2) const
  {
    int r = order;
    Vector<typename Vector2::value_type> Ux((r+1)*(r+1)), dUx((r+1)*(r+1));
    Ux.Fill(0); Vh.Fill(0); dUx.Fill(0);
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    int p0 = 0, p1 = 1, p2 = 2;
    switch (num_loc)
      {
      case 0 :
      case 5:
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*Uh(NumOrtho3D(k, i, j));
                dUx(NumQuad2D(i, j)) += dLi_Geom(k, ext)*Uh(NumOrtho3D(k, i, j));
              }
        
        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*Uh(NumOrtho3D(i, k, j));
                dUx(NumQuad2D(i, j)) += dLi_Geom(k, ext)*Uh(NumOrtho3D(i, k, j));
              }
        
        p0 = 1; p1 = 0; p2 = 2;
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Ux(NumQuad2D(i, j)) += Li_Geom(k, ext)*Uh(NumOrtho3D(i, j, k));
                dUx(NumQuad2D(i, j)) += dLi_Geom(k, ext)*Uh(NumOrtho3D(i, j, k));
              }
        
        p0 = 2; p1 = 0; p2 = 1;
        break;
      }
    
    if ((r2 == 0) || (r2 == order))
      {
	typename Vector1::value_type vx, vy, vz;
	for (int i = 0; i <= r; i++)
	  for (int k = 0; k <= r; k++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int j = 0; j <= r-i; j++)
		{
                  vx += ValLi(j, k)*dUx(NumQuad2D(i, j));
                  vy += ValLi(j, k)*Ux(NumQuad2D(i, j));
                  vz += ValDLi(j, k)*Ux(NumQuad2D(i, j));
                }
	      
	      for (int j = 0; j <= r; j++)
		{
                  Vh(3*NumQuad2D(j, k)+p0) += ValLi(i, j)*vx;
                  Vh(3*NumQuad2D(j, k)+p1) += ValDLi(i, j)*vy;
                  Vh(3*NumQuad2D(j, k)+p2) += ValLi(i, j)*vz;
                }
	    }
      }
    else
      {
	typename Vector1::value_type vx, vy, vz;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r2);
	const Matrix<Real_wp>& dAh = ProjOperatorDerivQuadOrder(r2);
	const Matrix<int>& Num = NumProjOperator(r2);
	for (int i = 0; i <= r; i++)
	  for (int k = 0; k <= r2; k++)
	    {
              vx = 0; vy = 0; vz = 0;
	      for (int j = 0; j <= r-i; j++)
		{
                  vx += Ah(j, k)*dUx(NumQuad2D(i, j));
                  vy += Ah(j, k)*Ux(NumQuad2D(i, j));
                  vz += dAh(j, k)*Ux(NumQuad2D(i, j));
                }
	      
	      for (int j = 0; j <= r2; j++)
		{
                  Vh(3*Num(j, k)+p0) += Ah(i, j)*vx;
                  Vh(3*Num(j, k)+p1) += dAh(i, j)*vy;
                  Vh(3*Num(j, k)+p2) += Ah(i, j)*vz;
                }
	    }
      }
  }
  
  
  //! integration against gradient of basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r2 order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronDgLegendre
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r2) const
  {
    int r = order;
    Vector<typename Vector2::value_type> Ux((r+1)*(r+1)), dUx((r+1)*(r+1));
    Ux.Fill(0); dUx.Fill(0);
    
    int p0 = 0, p1 = 1, p2 = 2;
    switch (num_loc)
      {
      case 1 :
      case 4:
        p0 = 1; p1 = 0; p2 = 2;
        break;
      case 2 :
      case 3:
        p0 = 2; p1 = 0; p2 = 1;
        break;
      }

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();        
    if ((r2 == 0) || (r2 == order))
      {
        typename Vector1::value_type vx, vy, vz;
	for (int i = 0; i <= r; i++)
	  for (int k = 0; k <= r; k++)
	    {
	      vx = 0; vy = 0; vz = 0;	      
	      for (int j = 0; j <= r; j++)
		{
                  vx += ValLi(i, j)*Uh(3*NumQuad2D(j, k)+p0);
                  vy += ValDLi(i, j)*Uh(3*NumQuad2D(j, k)+p1);
                  vz += ValLi(i, j)*Uh(3*NumQuad2D(j, k)+p2);
                }

              for (int j = 0; j <= r-i; j++)
		{
                  dUx(NumQuad2D(i, j)) += ValLi(j, k)*vx;
                  Ux(NumQuad2D(i, j)) += ValLi(j, k)*vy;
                  Ux(NumQuad2D(i, j)) += ValDLi(j, k)*vz;
                }
	      
	    }
      }
    else
      {
	typename Vector1::value_type vx, vy, vz;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r2);
	const Matrix<Real_wp>& dAh = ProjOperatorDerivQuadOrder(r2);
	const Matrix<int>& Num = NumProjOperator(r2);
	for (int j = 0; j <= r2; j++)
	  for (int k = 0; k <= r; k++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int i = 0; i <= r2; i++)
		{
                  vx += Ah(k, i)*Uh(3*Num(i, j)+p0);
                  vy += dAh(k, i)*Uh(3*Num(i, j)+p1);
                  vz += Ah(k, i)*Uh(3*Num(i, j)+p2);
                }
	      
	      for (int i = 0; i <= r-k; i++)
		{
                  dUx(NumQuad2D(k, i)) += Ah(i, j)*vx;
                  Ux(NumQuad2D(k, i)) += Ah(i, j)*vy;
                  Ux(NumQuad2D(k, i)) += dAh(i, j)*vz;
                }
	    }
      }
    
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    Mlt(alpha, Ux); Mlt(alpha, dUx);
    switch (num_loc)
      {
      case 0 :
      case 5:
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Vh(NumOrtho3D(k, i, j)) += Li_Geom(k, ext)*Ux(NumQuad2D(i, j));
                Vh(NumOrtho3D(k, i, j)) += dLi_Geom(k, ext)*dUx(NumQuad2D(i, j));
              }
        
        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Vh(NumOrtho3D(i, k, j)) += Li_Geom(k, ext)*Ux(NumQuad2D(i, j));
                Vh(NumOrtho3D(i, k, j)) += dLi_Geom(k, ext)*dUx(NumQuad2D(i, j));
              }
        
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r-i-j; k++)
              {
                Vh(NumOrtho3D(i, j, k)) += Li_Geom(k, ext)*Ux(NumQuad2D(i, j));
                Vh(NumOrtho3D(i, j, k)) += dLi_Geom(k, ext)*dUx(NumQuad2D(i, j));
              }
        
        break;
      }
  }

  
  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T>
  void HexahedronDgLegendre
  ::AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      A.AddInteraction(m+i, n+i, mass);
  }
  
  
  //! we add constant elementary matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass, C, D, E coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + mass M + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p + \sum_p  C(p, q) S^{p,q}
    where M is the mass matrix
    M_{i, j} \int_K \phi_j \phi_i  dx
    where S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
    where R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void HexahedronDgLegendre
  ::AddConstantElemMatrixGen(int m, int n, const T& mass, const TinyMatrix<T, Prop, 3, 3>& C,
			     const TinyVector<T, 3>& D, const TinyVector<T, 3>& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {        
    Real_wp* data = rh_const.GetData();
    long* ptr = rh_const.GetPtr();
    int* ind = rh_const.GetInd();
    
    for (int i = 0; i < nb_dof_loc; i++)
      {
	if (!null_term(0))
	  A.AddInteraction(m+i, n+i, mass);
	
	if (!null_term(2) || !null_term(3))
	  for (int k = ptr[i]; k < ptr[i+1]; k++)
	    {
	      int j = ind[k]/3;
	      int ncomp = ind[k]%3;
 	      if (!null_term(2))
		A.AddInteraction(m+i, n+j, data[k]*D(ncomp));
	      
	      if (!null_term(3))
		A.AddInteraction(m+j, n+i, data[k]*E(ncomp));
	    }
      }

    if (!null_term(1))
      for (int i1 = 0; i1 <= order; i1++)
	for (int i2 = 0; i2 <= order-i1; i2++)
	  for (int i3 = 0; i3 <= order-i1-i2; i3++)
	    {
	      int irow = NumOrtho3D(i1, i2, i3);
	      for (int k = 0; k <= order-i2-i3; k++)
		{
		  int jcol = NumOrtho3D(k, i2, i3);
		  A.AddInteraction(m+irow, n+jcol, C(0, 0)*stiff1D(i1, k));
		}
	      
	      for (int k = 0; k <= order-i1-i3; k++)
		{
		  int jcol = NumOrtho3D(i1, k, i3);
		  A.AddInteraction(m+irow, n+jcol, C(1, 1)*stiff1D(i2, k));
		}
	      
	      for (int k = 0; k <= order-i2-i1; k++)
		{
		  int jcol = NumOrtho3D(i1, i2, k);
		  A.AddInteraction(m+irow, n+jcol, C(2, 2)*stiff1D(i3, k));
		}
	      
	      for (int j = 0; j <= order; j++)
		for (int k = 0; k <= order-j-i3; k++)
		  {
		    int jcol = NumOrtho3D(j, k, i3);
		    A.AddInteraction(m+irow, n+jcol, C(0, 1)*CoefDLi(i1, j)*CoefDLi(k, i2)
				     + C(1, 0)*CoefDLi(i2, k)*CoefDLi(j, i1));
		  }
	      
	      for (int j = 0; j <= order; j++)
		for (int k = 0; k <= order-j-i2; k++)
		  {
		    int jcol = NumOrtho3D(j, i2, k);
		    A.AddInteraction(m+irow, n+jcol, C(0, 2)*CoefDLi(i1, j)*CoefDLi(k, i3)
				     + C(2, 0)*CoefDLi(i3, k)*CoefDLi(j, i1));
		  }
	      
	      for (int j = 0; j <= order; j++)
		for (int k = 0; k <= order-j-i1; k++)
		  {
		    int jcol = NumOrtho3D(i1, j, k);
		    A.AddInteraction(m+irow, n+jcol, C(1, 2)*CoefDLi(i2, j)*CoefDLi(k, i3)
				     + C(2, 1)*CoefDLi(i3, k)*CoefDLi(j, i2));
		  }
	    }
  }
  
  
  //! adds elementary matrix for a variable jacobian J_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A, C, D, E coefficients multiplied by omega_i J_i for each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + S + R + Rt
    where M_ij = \int A varphi_j varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
    where R_ij = \int E varphi_i  \nabla varphi_j
    where Rt_ij = \int D varphi_j  \nabla varphi_i
   */
  template<class T, class Prop>
  void HexahedronDgLegendre
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& A,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			     const Vector<TinyVector<T, 3> >& D,
			     const Vector<TinyVector<T, 3> >& E,
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension3, 1>::
      AddVariableElemMatrixOpt(off_row, off_col, A, C, D, E, null_term, mat);
  }
    
  
  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
   */
  template<class T>
  void HexahedronDgLegendre
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension3, 1>::AddVariableMassMatrixOpt(off_row, off_col, A, mat);
  }


  //! computation of projection between finite element of different orders  
  void HexahedronDgLegendre::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>&,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp Pn, dPn; VectR2 Points; IVect perm;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    ProjOperatorDerivQuadOrder.Reallocate(rmax+1);
    NumProjOperator.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate(order+1, r+1);
        ProjOperatorDerivQuadOrder(r).Reallocate(order+1, r+1);
	NumProjOperator(r).Reallocate(r+1, r+1);
	Points = Pts(i); perm.Reallocate(Points.GetM()); perm.Fill();
	Sort(Points, perm);
	int nb = 0;
	for (int j1 = 0; j1 <= r; j1++)
	  for (int j2 = 0; j2 <= r; j2++)
	    NumProjOperator(r)(j1, j2) = perm(nb++);
	
	for (int i2 = 0; i2 <= r; i2++)
	  {
	    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*Points(i2)(1) - 1.0, Pn, dPn);
	    for (int i1 = 0; i1 <= order; i1++)
	      {
                ProjOperatorQuadOrder(r)(i1, i2) = Pn(i1)*CoefLegendre(i1);
                ProjOperatorDerivQuadOrder(r)(i1, i2) = 2.0*dPn(i1)*CoefLegendre(i1);
              }
	  }	
      }
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point
    */
  void HexahedronDgLegendre::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    VectReal_wp Px, Py, Pz;
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*point_loc(0)-1.0, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*point_loc(1)-1.0, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*point_loc(2)-1.0, Pz);
    //DISP(Px); DISP(Py); DISP(Pz);
    phi.Reallocate(nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++)
          {
            int node = NumOrtho3D(i, j, k);
            phi(node) = Px(i)*Py(j)*Pz(k)*InvWeightPolynomial(node);
          }
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] dphi gradient of basis functions on point
  */
  void HexahedronDgLegendre::ComputeGradientPhiRef(const R3& point_loc, VectR3& dphi) const
  {
    VectReal_wp Px, Py, Pz, dPx, dPy, dPz;
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*point_loc(0)-1.0, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*point_loc(1)-1.0, Py, dPy);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*point_loc(2)-1.0, Pz, dPz);
    dphi.Reallocate(nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
        for (int k = 0; k <= order-i-j; k++)
          {
            int node = NumOrtho3D(i, j, k);
            dphi(node).Init(dPx(i)*Py(j)*Pz(k), Px(i)*dPy(j)*Pz(k), Px(i)*Py(j)*dPz(k));
            dphi(node) *= 2.0*InvWeightPolynomial(node);
          }
  }
  

  //! retrieves values of a single basis function on all quadrature points
  void HexahedronDgLegendre
  ::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    int i1 = CoordinateDofs(i, 0);
    int i2 = CoordinateDofs(i, 1);
    int i3 = CoordinateDofs(i, 2);
    for (int k1 = 0; k1 <= order; k1++)
      for (int k2 = 0; k2 <= order; k2++)
	for (int k3 = 0; k3 <= order; k3++)
	  phi(NumQuad3D(k1, k2, k3)) = ValLi(i1, k1)*ValLi(i2, k2)*ValLi(i3, k3);
  }

  
  //! retrieves values and gradients of a single basis function on all quadrature points
  void HexahedronDgLegendre
  ::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR3& grad_phi) const
  {
    int i1 = CoordinateDofs(i, 0);
    int i2 = CoordinateDofs(i, 1);
    int i3 = CoordinateDofs(i, 2);
    for (int k1 = 0; k1 <= order; k1++)
      for (int k2 = 0; k2 <= order; k2++)
	for (int k3 = 0; k3 <= order; k3++)
	  {
	    int n = NumQuad3D(k1, k2, k3);
	    phi(n) = ValLi(i1, k1)*ValLi(i2, k2)*ValLi(i3, k3);
	    grad_phi(n)(0) = ValDLi(i1, k1)*ValLi(i2, k2)*ValLi(i3, k3);
	    grad_phi(n)(1) = ValLi(i1, k1)*ValDLi(i2, k2)*ValLi(i3, k3);
	    grad_phi(n)(2) = ValLi(i1, k1)*ValLi(i2, k2)*ValDLi(i3, k3);
	  }
  }
  

  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */  
  void HexahedronDgLegendre::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad3D(k, 0);
        int k2 = CoordinateQuad3D(k, 1);
        int k3 = CoordinateQuad3D(k, 2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order-i; j++)
            for (int k = 0; k <= order-i-j; k++)
              phi(NumOrtho3D(i, j, k)) = ValLi(i, k1)*ValLi(j, k2)*ValLi(k, k3);
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/((order+1)*(order+1));
        int ks = (k-nb_points_quadrature_inside)%((order+1)*(order+1));
        int k1 = CoordinateQuad2D(ks, 0);
        int k2 = CoordinateQuad2D(ks, 1);
        int ext = 0;
        if (num_loc > 2)
          ext = this->GetGeometryOrder();
        
        switch (num_loc)
          {
          case 0 :
          case 5 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order-i; j++)
                  for (int k = 0; k <= order-i-j; k++)
                    phi(NumOrtho3D(k, i, j)) = ValLi(i, k1)*ValLi(j, k2)*Li_Geom(k, ext);
            }
            break;
          case 1:
          case 4 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order-i; j++)
                  for (int k = 0; k <= order-i-j; k++)
                    phi(NumOrtho3D(i, k, j)) = ValLi(i, k1)*ValLi(j, k2)*Li_Geom(k, ext);
            }
            break;
          case 2:
          case 3 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order-i; j++)
                  for (int k = 0; k <= order-i-j; k++)
                    phi(NumOrtho3D(i, j, k)) = ValLi(i, k1)*ValLi(j, k2)*Li_Geom(k, ext);
            }
            break;
          }
      }
  }
  
  
  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
   */
  void HexahedronDgLegendre::GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad3D(k, 0);
        int k2 = CoordinateQuad3D(k, 1);
        int k3 = CoordinateQuad3D(k, 2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order-i; j++)
            for (int k = 0; k <= order-i-j; k++)
              grad_phi(NumOrtho3D(i, j, k)).Init(ValDLi(i, k1)*ValLi(j, k2)*ValLi(k, k3),
                                                 ValLi(i, k1)*ValDLi(j, k2)*ValLi(k, k3),
                                                 ValLi(i, k1)*ValLi(j, k2)*ValDLi(k, k3));
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/((order+1)*(order+1));
        int ks = (k-nb_points_quadrature_inside)%((order+1)*(order+1));
        int k1 = CoordinateQuad2D(ks, 0);
        int k2 = CoordinateQuad2D(ks, 1);
        int ext = 0;
        if (num_loc > 2)
          ext = this->GetGeometryOrder();
        
        switch (num_loc)
          {
          case 0 :
          case 5 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order-i; j++)
                  for (int k = 0; k <= order-i-j; k++)
                    grad_phi(NumOrtho3D(k, i, j)).Init(ValLi(i, k1)*ValLi(j, k2)*dLi_Geom(k, ext),
                                                       ValDLi(i, k1)*ValLi(j, k2)*Li_Geom(k, ext),
                                                       ValLi(i, k1)*ValDLi(j, k2)*Li_Geom(k, ext));
            }
            break;
          case 1:
          case 4 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order-i; j++)
                  for (int k = 0; k <= order-i-j; k++)
                    grad_phi(NumOrtho3D(i, k, j)).Init(ValDLi(i, k1)*ValLi(j, k2)*Li_Geom(k, ext),
                                                       ValLi(i, k1)*ValLi(j, k2)*dLi_Geom(k, ext),
                                                       ValLi(i, k1)*ValDLi(j, k2)*Li_Geom(k, ext));
            }
            break;
          case 2:
          case 3 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order-i; j++)
                  for (int k = 0; k <= order-i-j; k++)
                    grad_phi(NumOrtho3D(i, j, k)).Init(ValDLi(i, k1)*ValLi(j, k2)*Li_Geom(k, ext),
                                                       ValLi(i, k1)*ValDLi(j, k2)*Li_Geom(k, ext),
                                                       ValLi(i, k1)*ValLi(j, k2)*dLi_Geom(k, ext));
            }
            break;
          }
      }
  }

} // namespace Montjoie
  
#define MONTJOIE_FILE_HEXAHEDRON_DG_LEGENDRE_CXX
#endif
