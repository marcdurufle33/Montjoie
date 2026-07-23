#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_GAUSS_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronDgGauss::HexahedronDgGauss() : HexahedronGauss()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    mass_lumping = true;
    diagonal_mass = true;
    this->discontinuous_element = true;
    use_quadrature_for_rh = false;
    use_quadrature_for_sh = false;
  }
  

  size_t HexahedronDgGauss::GetMemorySize() const
  {
    size_t taille = HexahedronGauss::GetMemorySize();
    taille += G_Geom.GetMemorySize();
    taille += dG_Geom.GetMemorySize();
    return taille;
  }

  
  //! construction of finite element class
  void HexahedronDgGauss::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf_tri, int rsurf_quad,
						 int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    HexahedronGauss::ConstructFiniteElement(r, rgeom, r, type_quad);
    
    this->elt_geom.dof_equal_nodal = false;

    lob_basis = lob_quad;
    
    this->SetPointsDof1D(this->Points1D());
    this->SetPointsDof2D_quad(this->Points2D_quad());
    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_points_quadrature_inside);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points_dof3d(i) = this->PointsND(i);
    
    this->SetPointsDofND(points_dof3d);
    NumDofs3D = NumQuad3D;
    CoordinateDofs = CoordinateQuad3D;
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    int order_geom = this->GetGeometryOrder();
    G_Geom.Reallocate(order+1, order_geom+1);
    dG_Geom.Reallocate(order+1, order_geom+1);
    const VectReal_wp& points_nodal1d = this->PointsNodal1D();
    for (int i = 0; i <= order; i++)
      {
        for (int j = 0; j <= order_geom; j++)
	  {
            G_Geom(i, j) = lob_basis.EvaluatePhi(i, points_nodal1d(j));
            dG_Geom(i, j) = lob_basis.EvaluatePhiGrad(i, points_nodal1d(j));
          }
      }

    this->ConstructStiffnessMatrix();
    
    Matrix<Real_wp, General, ArrayRowSparse> sh0, sh1, sh2, sh3, sh4, sh5;
    int N = (order+1)*(order+1)*(order+1);
    int Nquad = (order+1)*(order+1);
    sh0.Reallocate(N, Nquad); sh1.Reallocate(N, Nquad); sh2.Reallocate(N, Nquad);
    sh3.Reallocate(N, Nquad); sh4.Reallocate(N, Nquad); sh5.Reallocate(N, Nquad);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int nv = NumDofs3D(i,j,k);
            int ns = NumDofs2D(j, k);
            sh0.AddInteraction(nv, ns, G_Geom(i, 0));
            sh5.AddInteraction(nv, ns, G_Geom(i, order_geom));

            ns = NumDofs2D(i, k);
            sh1.AddInteraction(nv, ns, G_Geom(j, 0));
            sh4.AddInteraction(nv, ns, G_Geom(j, order_geom));
            
            ns = NumDofs2D(i, j);
            sh2.AddInteraction(nv, ns, G_Geom(k, 0));
            sh3.AddInteraction(nv, ns, G_Geom(k, order_geom));
          }
    
    Copy(sh0, sparse_const_sh(0)); Copy(sh1, sparse_const_sh(1));
    Copy(sh2, sparse_const_sh(2)); Copy(sh3, sparse_const_sh(3));
    Copy(sh4, sparse_const_sh(4)); Copy(sh5, sparse_const_sh(5));

    delete function_basis_quad;
    function_basis_quad = new QuadrangleDgGauss();
    function_basis_quad->ConstructFiniteElement(r, rgeom, rquad, type_quad);
    element_quad_surf = function_basis_quad;
  }
  
  
  //! computation of U on nodal points
  /*!
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
   */
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::
  ComputeNodalValuesGen(const Vector1& Uloc, Vector2& Uloc_node) const
  {
    int order_geom = this->GetGeometryOrder();
    Vector2 U1((order_geom+1)*(order+1)*(order+1)), U2((order_geom+1)*(order_geom+1)*(order+1));
    FillZero(U1); FillZero(U2);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order_geom; k++)
	  {
	    int node = (order_geom+1)*(i*(order+1) + j) + k;
	    for (int m = 0; m <= order; m++)
	      U1(node) += G_Geom(m, k)*Uloc(NumDofs3D(i, j, m));
	  }

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k <= order_geom; k++)
	  {
            int node = (order_geom+1)*(i*(order_geom+1) + j) + k;
	    for (int m = 0; m <= order; m++)
	      U2(node) += G_Geom(m, j)*U1((order_geom+1)*(i*(order+1) + m) + k);
	  }
    
    Uloc_node.Reallocate(this->GetNbPointsNodalElt());
    FillZero(Uloc_node);
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k <= order_geom; k++)
	  {
	    int node = this->GetNumNodes3D(i, j, k);
	    for (int m = 0; m <= order; m++)
              {
                int p = (order_geom+1)*(m*(order_geom+1) + j) + k;
                Uloc_node(node) += G_Geom(m, i)*U2(p);
              }
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
  void HexahedronDgGauss::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  //! computation of u on nodal points of a face
  /*!
    \param[in] feval dof components of u
    \param[out] res values of u on quadrature points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::
  ComputeValueBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    int ext = 0;
    int order_geom = this->GetGeometryOrder();
    if (num_loc >= 3)
      ext = order_geom;
    
    int Nquad = (order+1)*(order+1);
    Vector2 Udof(Nquad);
    Udof.Fill(0);
    switch (num_loc)
      {
      case 0 :
      case 5:
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      {
		int num_dof = NumDofs3D(i,j,k);
		Udof(j*(order+1)+k) += G_Geom(i,ext)*feval(num_dof);
	      }
        break;
      case 1 :
      case 4 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      {
		int num_dof = NumDofs3D(i,j,k);
		Udof(i*(order+1)+k) += G_Geom(j,ext)*feval(num_dof);
	      }
        break;
      case 2 :
      case 3 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      {
		int num_dof = NumDofs3D(i,j,k);
		Udof(i*(order+1)+j) += G_Geom(k,ext)*feval(num_dof);
	      }
        break;
      }
    
    res.Fill(0.0);
    Vector2 Ux((order_geom+1)*(order+1));
    Ux.Fill(0);
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order; j++)
        {
          int p = i*(order+1) + j;
          for (int k = 0; k <= order; k++)
            Ux(p) += G_Geom(k, i)*Udof(k*(order+1)+j);
        }
    
    const Matrix<int>& NumNodes2D_quad = this->Fb_geom.GetNumNodes2D_quad();
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
        {
          int p = NumNodes2D_quad(i, j);
          for (int k = 0; k <= order; k++)
            res(p) += G_Geom(k, j)*Ux(i*(order+1)+k);
        }
  }

  
  //! computation of grad u on nodal points of a face
  /*!
    \param[in] feval dof components of u
    \param[out] res gradient of u on quadrature points of the face
    \param[in] num_loc local position of the face in the element
  */  
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::
  ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    int ext = 0;
    int order_geom = this->GetGeometryOrder();
    if (num_loc >= 3)
      ext = order_geom;

    const Matrix<int>& NumNodes2D_quad = this->Fb_geom.GetNumNodes2D_quad();    
    int Nquad = (order+1)*(order+1);
    Vector1 Udof(Nquad), dUdof(Nquad);
    Vector1 Ux(3*(order_geom+1)*(order+1));
    Udof.Fill(0); dUdof.Fill(0);
    switch (num_loc)
      {
      case 0 :
      case 5:
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      {
		int num_dof = NumDofs3D(i,j,k);
		Udof(j*(order+1)+k) += G_Geom(i,ext)*feval(num_dof);
                dUdof(j*(order+1)+k) += dG_Geom(i,ext)*feval(num_dof);
	      }
        
        res.Fill(0.0);        
        Ux.Fill(0);
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            {
              int p = 3*(i*(order+1) + j);
              for (int k = 0; k <= order; k++)
                {
                  Ux(p) += G_Geom(k, i)*dUdof(k*(order+1)+j);
                  Ux(p+1) += dG_Geom(k, i)*Udof(k*(order+1)+j);
                  Ux(p+2) += G_Geom(k, i)*Udof(k*(order+1)+j);
                }
            }
    
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            {
              int p = 3*NumNodes2D_quad(i, j);
              for (int k = 0; k <= order; k++)
                {
                  int q = 3*(i*(order+1)+k);
                  res(p) += G_Geom(k, j)*Ux(q);
                  res(p+1) += G_Geom(k, j)*Ux(q+1);
                  res(p+2) += dG_Geom(k, j)*Ux(q+2);
                }
            }
        
        break;
      case 1 :
      case 4 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      {
		int num_dof = NumDofs3D(i,j,k);
		Udof(i*(order+1)+k) += G_Geom(j,ext)*feval(num_dof);
                dUdof(i*(order+1)+k) += dG_Geom(j,ext)*feval(num_dof);
	      }
        
        res.Fill(0.0);        
        Ux.Fill(0);
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            {
              int p = 3*(i*(order+1) + j);
              for (int k = 0; k <= order; k++)
                {
                  Ux(p) += dG_Geom(k, i)*Udof(k*(order+1)+j);
                  Ux(p+1) += G_Geom(k, i)*dUdof(k*(order+1)+j);
                  Ux(p+2) += G_Geom(k, i)*Udof(k*(order+1)+j);
                }
            }
    
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            {
              int p = 3*NumNodes2D_quad(i, j);
              for (int k = 0; k <= order; k++)
                {
                  int q = 3*(i*(order+1)+k);
                  res(p) += G_Geom(k, j)*Ux(q);
                  res(p+1) += G_Geom(k, j)*Ux(q+1);
                  res(p+2) += dG_Geom(k, j)*Ux(q+2);
                }
            }
        break;
      case 2 :
      case 3 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    for (int k = 0; k <= order; k++)
	      {
		int num_dof = NumDofs3D(i,j,k);
		Udof(i*(order+1)+j) += G_Geom(k,ext)*feval(num_dof);
                dUdof(i*(order+1)+j) += dG_Geom(k,ext)*feval(num_dof);
	      }

        res.Fill(0.0);        
        Ux.Fill(0);
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            {
              int p = 3*(i*(order+1) + j);
              for (int k = 0; k <= order; k++)
                {
                  Ux(p) += dG_Geom(k, i)*Udof(k*(order+1)+j);
                  Ux(p+1) += G_Geom(k, i)*Udof(k*(order+1)+j);
                  Ux(p+2) += G_Geom(k, i)*dUdof(k*(order+1)+j);
                }
            }
    
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            {
              int p = 3*NumNodes2D_quad(i, j);
              for (int k = 0; k <= order; k++)
                {
                  int q = 3*(i*(order+1)+k);
                  res(p) += G_Geom(k, j)*Ux(q);
                  res(p+1) += dG_Geom(k, j)*Ux(q+1);
                  res(p+2) += G_Geom(k, j)*Ux(q+2);
                }
            }
        break;
      }
    

  }
  
  
  //! x is overwritten by M^-1 x where M is the mass matrix
  template<class Vector1>
  void HexahedronDgGauss::SolveMassMatrixGen(Vector1& x) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      x(i) *= invWeights3d(i);
  }
  
  
  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix
  template<class Vector1>
  void HexahedronDgGauss::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      x(i) *= invSqrtWeights3d(i);
  }
    
  
  //! x is overwritten by M x where M is the mass matrix
  template<class Vector1>
  void HexahedronDgGauss::MltMassMatrixGen(Vector1& x) const
  {
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      x(i) *= weights3d(i);
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
  void HexahedronDgGauss::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    // mass lumping, Ch = Id
    Copy(Uh, Vh);
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
  void HexahedronDgGauss::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    // mass lumping, Ch = Id
    Copy(Uh, Vh);
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
  void HexahedronDgGauss::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(rh_loc, Uh, Vh);
  }
  
  
  //! Integration against derivatives of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Ux Ux_i = \int_K f d/dx(\varphi_i)  dx
    \param[out] Uy Uy_i = \int_K f d/dy(\varphi_i)  dx
    \param[out] Uz Uz_i = \int_K f d/dz(\varphi_i)  dx
    This operation is equivalent to a matrix vector product
    Ux = Rh^1 Uh, Uy = Rh^2 Uh, Uz = Rh^3 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
    where (Rh)^3_{i,j} = d/dz(\varphi_i(\xi_j) )
   */  
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::
  ApplyRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
    Mlt(rh_locX, Uh, Ux);
    Mlt(rh_locY, Uh, Uy);
    Mlt(rh_locZ, Uh, Uz);
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
  void HexahedronDgGauss::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Uh = Rh Vh
   */
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::ApplyConstantRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    Mlt(rh_weight, Vh, Uh);
  }
  
  
  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    (Rh^3)_{i, j} = \int_K \varphi_j d/dz(\varphi_i) dx
    then Ux = Rh^1 Uh, Uy = Rh^2 Uh, Uz = Rh^3 Uh
   */
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::
  ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
    Mlt(rh_weightX, Uh, Ux);
    Mlt(rh_weightY, Uh, Uy);
    Mlt(rh_weightZ, Uh, Uz);
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_weight, Uh, Vh);
  }

  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the face
    \param[in] r order of quadrature rule of the face
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronDgGauss
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    int ext = 0;
    if (num_loc >= 3)
      ext = this->GetGeometryOrder();
    
    if ((r == 0) || (r == order))
      {
        Mlt(SeldonTrans, sparse_const_sh(num_loc), Uh, Vh);
        return;
        
        Vh.Fill(0);
	switch (num_loc)
	  {
	  case 0 :
	  case 5 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      Vh(NumDofs2D(j,k)) += G_Geom(i, ext)*Uh(num_dof);
		    }
	    }
	    break;
	  case 1 :
	  case 4 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      Vh(NumDofs2D(i,k)) += G_Geom(j, ext)*Uh(num_dof);
		    }
	    }
	    break;
	  case 2 :
	  case 3 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      Vh(NumDofs2D(i,j)) += G_Geom(k, ext)*Uh(num_dof);
		    }
	    }
	    break;
	  }
      }
    else
      {
	Vector2 Vquad((order+1)*(order+1));
	Vh.Fill(0); Vquad.Fill(0);
	switch (num_loc)
	  {
	  case 0 :
	  case 5 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      Vquad(NumDofs2D(j,k)) += G_Geom(i, ext)*Uh(num_dof);
		    }
	    }
	    break;
	  case 1 :
	  case 4 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      Vquad(NumDofs2D(i,k)) += G_Geom(j, ext)*Uh(num_dof);
		    }
	    }
	    break;
	  case 2 :
	  case 3 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      Vquad(NumDofs2D(i,j)) += G_Geom(k, ext)*Uh(num_dof);
		    }
	    }
	    break;
	  }
	
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vloc = 0;
	      for (int j = 0; j <= order; j++)
		vloc += Ah(j, m)*Vquad(NumDofs2D(i, j));
	      
	      for (int n = 0; n <= r; n++)
		Vh(Num(n, m)) += Ah(i, n)*vloc;
	    }
      }
  }
  
   
  //! integration against basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] feval vector containing \omega_k f(\xi_k) 
    \param[out] res res_i = res_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronDgGauss::ApplyShGen(const T0& alpha, int num_loc, const Vector1& feval,
				     Vector2& res, int r) const
  {
    int ext = 0;
    if (num_loc >= 3)
      ext = this->GetGeometryOrder();
    
    if ((r == 0) || (r == order))
      {
        MltAdd(alpha, sparse_const_sh(num_loc), feval, T0(1), res);
        return;
        
	switch (num_loc)
	  {
	  case 0 :
	  case 5 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      res(num_dof) += alpha*G_Geom(i, ext)*feval(NumDofs2D(j,k));
		    }
	    }
	    break;
	  case 1 :
	  case 4 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      res(num_dof) += alpha*G_Geom(j, ext)*feval(NumDofs2D(i,k));
		    }
	    }
	    break;
	  case 2 :
	  case 3 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      res(num_dof) += alpha*G_Geom(k, ext)*feval(NumDofs2D(i,j));
		    }
	    }
	    break;
	  }
      }
    else
      {
	Vector1 Vquad((order+1)*(order+1));
	Vquad.Fill(0);
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vloc = 0;
	      for (int n = 0; n <= r; n++)
		vloc += Ah(i, n)*feval(Num(n, m));
	      
	      for (int j = 0; j <= order; j++)
		Vquad(NumDofs2D(i, j)) += Ah(j, m)*vloc;
	    }
	
	switch (num_loc)
	  {
	  case 0 :
	  case 5 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      res(num_dof) += alpha*G_Geom(i, ext)*Vquad(NumDofs2D(j,k));
		    }
	    }
	    break;
	  case 1 :
	  case 4 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      res(num_dof) += alpha*G_Geom(j, ext)*Vquad(NumDofs2D(i,k));
		    }
	    }
	    break;
	  case 2 :
	  case 3 :
	    {
	      for (int i = 0; i <= order; i++)
		for (int j = 0; j <= order; j++)
		  for (int k = 0; k <= order; k++)
		    {
		      int num_dof = NumDofs3D(i,j,k);
		      res(num_dof) += alpha*G_Geom(k, ext)*Vquad(NumDofs2D(i,j));
		    }
	    }
	    break;
	  }		
      }
  }
  
  
  //! computation of gradient of u on quadrature points of a face
  /*!
    \param[in] num_loc face number
    \param[in] res components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the face
    \param[in] r order of quadrature rule of the face
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>
  void HexahedronDgGauss::ApplyNablaShTransposeGen(int num_loc, const Vector1& res,
						   Vector2& Vh, int r) const
  {
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    int Nquad = (order+1)*(order+1);
    Vector2 Udof(Nquad), feval(3*Nquad);
    Udof.Fill(0);
    feval.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0 :
      case 5 :
        {
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(i, j, k);
                    Udof(p) += G_Geom(i, ext)*res(node);
                    feval(3*p) += dG_Geom(i, ext)*res(node);
                  }
              }
          
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    feval(3*NumQuad2D(k, j) + 1) += lob_basis.GradPhi(i, k)*Udof(num);
                    feval(3*NumQuad2D(i, k) + 2) += lob_basis.GradPhi(j, k)*Udof(num);
                  }
              }
        }
        break;
      case 1 :
      case 4 :
        {
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, i, k);
                    feval(3*p+1) += dG_Geom(i, ext)*res(node);
                    Udof(p) += G_Geom(i, ext)*res(node);
                  }
              }
          
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    feval(3*NumQuad2D(k, j)) += lob_basis.GradPhi(i, k)*Udof(num);
                    feval(3*NumQuad2D(i, k) + 2) += lob_basis.GradPhi(j, k)*Udof(num);
                  }
              }
        }
        break;
      case 2 :
      case 3 :
        {
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, k, i);
                    feval(3*p+2) += dG_Geom(i, ext)*res(node);
                    Udof(p) += G_Geom(i, ext)*res(node);
                  }
              }
          
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    feval(3*NumQuad2D(k, j)) += lob_basis.GradPhi(i, k)*Udof(num);
                    feval(3*NumQuad2D(i, k) + 1) += lob_basis.GradPhi(j, k)*Udof(num);
                  }
              }
        }
        break;
      }

    if ((r != 0) && (r != this->order)) 
      {
        typename Vector1::value_type vx, vy, vz;
	Vh.Fill(0);
        const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
        for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int j = 0; j <= order; j++)
		{
                  vx += Ah(j, m)*feval(3*NumDofs2D(i, j));
                  vy += Ah(j, m)*feval(3*NumDofs2D(i, j)+1);
                  vz += Ah(j, m)*feval(3*NumDofs2D(i, j)+2);
                }
	      
	      for (int n = 0; n <= r; n++)
		{
                  Vh(3*Num(n, m)) += Ah(i, n)*vx;
                  Vh(3*Num(n, m)+1) += Ah(i, n)*vy;
                  Vh(3*Num(n, m)+2) += Ah(i, n)*vz;
                }
	    }
      }    
    else
      Copy(feval, Vh);
  }
  
  
  //! integration against gradient of basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] res res_i = res_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is the quadrature rule of face num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronDgGauss::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
					  Vector2& res, int r) const
  {
    int Nquad = (order+1)*(order+1);
    Vector2 feval(3*Nquad);

    if ((r != 0) && (r != this->order)) 
      {
        typename Vector1::value_type vx, vy, vz;
        feval.Fill(0);
        const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int n = 0; n <= r; n++)
		{
                  vx += Ah(i, n)*Uh(3*Num(n, m));
                  vy += Ah(i, n)*Uh(3*Num(n, m)+1);
                  vz += Ah(i, n)*Uh(3*Num(n, m)+2);
                }
	      
	      for (int j = 0; j <= order; j++)
		{
                  feval(3*NumDofs2D(i, j)) += Ah(j, m)*vx;
                  feval(3*NumDofs2D(i, j)+1) += Ah(j, m)*vy;
                  feval(3*NumDofs2D(i, j)+2) += Ah(j, m)*vz;
                }
	    }
      }    
    else
      Copy(Uh, feval);
    
    Mlt(alpha, feval);
    
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();
    
    Vector2 Udof(Nquad);
    Udof.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0 :
      case 5 :
        {
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    Udof(num) += lob_basis.GradPhi(i, k)*feval(3*NumQuad2D(k, j) + 1);
                    Udof(num) += lob_basis.GradPhi(j, k)*feval(3*NumQuad2D(i, k) + 2);
                  }
              }
          
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(i, j, k);
                    res(node) += dG_Geom(i, ext)*feval(3*p) + G_Geom(i, ext)*Udof(p);
                  }
              }

        }
        break;
      case 1 :
      case 4 :
        {
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    Udof(num) += lob_basis.GradPhi(i, k)*feval(3*NumQuad2D(k, j));
                    Udof(num) += lob_basis.GradPhi(j, k)*feval(3*NumQuad2D(i, k) + 2);
                  }
              }
          
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, i, k);
                    res(node) += dG_Geom(i, ext)*feval(3*p+1) + G_Geom(i, ext)*Udof(p);
                  }
              }
        }
        break;
      case 2 :
      case 3 :
        {
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    Udof(num) += lob_basis.GradPhi(i, k)*feval(3*NumQuad2D(k, j));
                    Udof(num) += lob_basis.GradPhi(j, k)*feval(3*NumQuad2D(i, k) + 1);
                  }
              }
          
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, k, i);
                    res(node) += dG_Geom(i, ext)*feval(3*p+2) + G_Geom(i, ext)*Udof(p);
                  }
              }
        }
        break;
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
  void HexahedronDgGauss
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
  void HexahedronDgGauss::AddVariableMassMatrixGen(int off_row, int off_col,
						   const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      mat.AddInteraction(i+off_row, i+off_col, A(i));
  }


  //! retrieves values of a single basis function on all quadrature points
  void HexahedronDgGauss::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    phi.Zero();
    phi(i) = Real_wp(1);
  }

  
  //! retrieves values and gradients of a single basis function on all quadrature points
  void HexahedronDgGauss::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR3& grad_phi) const
  {
    phi.Zero();
    grad_phi.Zero();
    phi(i) = Real_wp(1);
    int i1 = CoordinateDofs(i, 0);
    int i2 = CoordinateDofs(i, 1);
    int i3 = CoordinateDofs(i, 2);
    for (int k = 0; k <= order; k++)
      {
	grad_phi(NumQuad3D(k, i2, i3))(0) = lob_basis.GradPhi(i1, k);
	grad_phi(NumQuad3D(i1, k, i3))(1) = lob_basis.GradPhi(i2, k);
	grad_phi(NumQuad3D(i1, i2, k))(2) = lob_basis.GradPhi(i3, k);
      }
  }
  

  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void HexahedronDgGauss::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    if (k < nb_points_quadrature_inside)
      {
        // interior quadrature point, we are using mass lumping property
        phi(k) = 1.0;
      }
    else
      {
        // quadrature point on the boundary
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
                phi(NumDofs3D(i, k1, k2)) = G_Geom(i, ext);
            }
            break;
          case 1:
          case 4 :
            {
              for (int i = 0; i <= order; i++)
                phi(NumDofs3D(k1, i, k2)) = G_Geom(i, ext);
            }
            break;
          case 2:
          case 3 :
            {
              for (int i = 0; i <= order; i++)
                phi(NumDofs3D(k1, k2, i)) = G_Geom(i, ext);
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
  void HexahedronDgGauss::GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    FillZero(grad_phi);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad3D(k, 0);
        int k2 = CoordinateQuad3D(k, 1);
        int k3 = CoordinateQuad3D(k, 2);
        for (int i = 0; i <= order; i++)
          {
            grad_phi(NumDofs3D(i, k2, k3))(0) = lob_basis.GradPhi(i, k1);
            grad_phi(NumDofs3D(k1, i, k3))(1) = lob_basis.GradPhi(i, k2);
            grad_phi(NumDofs3D(k1, k2, i))(2) = lob_basis.GradPhi(i, k3);
          }
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
          case 0:
          case 5:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumDofs3D(i, k1, k2))(0) = dG_Geom(i, ext);
                  for (int j = 0; j <= order; j++)
                    {
                      grad_phi(NumDofs3D(i, j, k2))(1) = G_Geom(i, ext)*lob_basis.GradPhi(j, k1);
                      grad_phi(NumDofs3D(i, k1, j))(2) = G_Geom(i, ext)*lob_basis.GradPhi(j, k2);
                    }
                }
              
            }
            break;
          case 1:
          case 4:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumDofs3D(k1, i, k2))(1) = dG_Geom(i, ext);
                  for (int j = 0; j <= order; j++)
                    {
                      grad_phi(NumDofs3D(j, i, k2))(0) = G_Geom(i, ext)*lob_basis.GradPhi(j, k1);
                      grad_phi(NumDofs3D(k1, i, j))(2) = G_Geom(i, ext)*lob_basis.GradPhi(j, k2);
                    }
                }
            }
            break;
          case 2:
          case 3:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumDofs3D(k1, k2, i))(2) = dG_Geom(i, ext);
                  for (int j = 0; j <= order; j++)
                    {
                      grad_phi(NumDofs3D(j, k2, i))(0) = G_Geom(i, ext)*lob_basis.GradPhi(j, k1);
                      grad_phi(NumDofs3D(k1, j, i))(1) = G_Geom(i, ext)*lob_basis.GradPhi(j, k2);
                    }
                }
            }
            break;
          }
      }
  }

} // namespace Montjoie
  
#define MONTJOIE_FILE_HEXAHEDRON_DG_GAUSS_CXX
#endif
