#ifndef MONTJOIE_FILE_QUADRANGLE_DG_GAUSS_CXX

namespace Montjoie
{
  
  //! default constructor
  QuadrangleDgGauss::QuadrangleDgGauss() : QuadrangleGauss()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    this->discontinuous_element = true;
  }
  
  
  size_t QuadrangleDgGauss::GetMemorySize() const
  {
    size_t taille = QuadrangleGauss::GetMemorySize();
    taille += ValGauss_Extremity.GetMemorySize();
    taille += dG_GL.GetMemorySize();
    taille += G_Geom.GetMemorySize();
    taille += dG_Geom.GetMemorySize();
    return taille;
  }

  
  //! constructing finite element
  void QuadrangleDgGauss::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf, int type_surf)
  {
    QuadrangleGauss::ConstructFiniteElement(r, rgeom, rquad, Globatto<Real_wp>::QUADRATURE_GAUSS);

    // we choose gauss points for degrees of freedom
    lob_basis.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    lob_basis.ComputeGradPhi();
    
    int order_geom = this->GetGeometryOrder();
    G_GL.Reallocate(order+1, order_geom+1);
    dG_GL.Reallocate(order+1, order_geom+1);

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        G_GL(i,j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        dG_GL(i,j) = lob_basis.EvaluatePhiGrad(i, this->PointsNodal1D(j));
    
    // 1-D stiffness matrix
    stiffness_matrix1D.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  stiffness_matrix1D(i,j) = 0;
	  for (int k = 0; k <= order; k++)
            stiffness_matrix1D(i, j) += lob_basis.Weights(k)*lob_basis.GradPhi(i,k)*lob_basis.GradPhi(j,k);
	}
    
    G_Geom.Reallocate(order + 1, order_geom + 1);
    dG_Geom.Reallocate(order + 1, order_geom + 1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        {
          G_Geom(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
          dG_Geom(i, j) = lob_basis.EvaluatePhiGrad(i, this->PointsNodal1D(j));
        }
    
    // Gauss points for dofs
    this->SetPointsDof1D(lob_basis.Points());
    
    VectR2 points_dof2d;
    points_dof2d = this->PointsND();
    this->nb_points_dof_inside = (order+1)*(order+1);
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    
    if (order == order_quad)
      {
        NumFct2D = NumQuad2D;
        CoordinateDofs = CoordinateQuad;
        weights2d = this->WeightsND();
        this->num_dof_points_surf = this->num_quad_points_surf;
        this->elt_geom.dof_equal_quadrature = true;
        diagonal_mass = true;
        mass_lumping = true;
        this->use_quadrature_for_rh = false;
      }
    else
      {
        diagonal_mass = false;
        mass_lumping = false;
        this->use_quadrature_for_rh = true;

        weights2d.Reallocate(nb_dof_loc);
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            {
              weights2d(NumFct2D(i, j)) = lob_basis.Weights(i)*lob_basis.Weights(j);
              points_dof2d(this->NumFct2D(i,j))
                = R2(this->PointsDof1D(i), this->PointsDof1D(j));
            }
        
        invWeights2d.Reallocate(nb_dof_loc);
        invSqrtWeights2d.Reallocate(nb_dof_loc);
        for (int k = 0; k < nb_dof_loc; k++)
          {
            invWeights2d(k) = 1.0/weights2d(k);
            invSqrtWeights2d(k) = 1.0/sqrt(weights2d(k));
          }
        
        // num_dof_points_surf already filled in QuadrangleGauss
        this->elt_geom.dof_equal_quadrature = false;

        this->ConstructMassMatrix();
      }
    
    this->SetPointsDofND(points_dof2d);
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ValGauss_Extremity.Reallocate(order+1, 2);
    for (int i = 0; i <= order; i++)
      {
	ValGauss_Extremity(i, 0) = lob_basis.EvaluatePhi(i, Real_wp(0));
	ValGauss_Extremity(i, 1) = lob_basis.EvaluatePhi(i, Real_wp(1));
      }
    
    ConstructElementaryMatrix(*this);
    ConstructStiffnessMatrix();
  }


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  void QuadrangleDgGauss
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r, r, rquad, type_quad,
				 EdgeGauss::GAUSS);
    
    this->element_surface = edge;
  }
  
  
  //! constructing stiffness matrix
  void QuadrangleDgGauss::ConstructStiffnessMatrix()
  {
    // matrix R_h and transpose
    int Nquad = nb_points_quadrature_inside;
    Matrix<Real_wp, General, ArrayRowSparse> Rh(Nquad, 2*Nquad);

    int numero, numpsi;
    // loop on all quadrature points of the square (i,j)
    for (int i = 0; i <= order_quad; i++)
      {
	for (int j = 0; j <= order_quad; j++)
	  {
	    // considered dof
	    int num_loc = NumQuad2D(i,j);
	    
	    // derivatives with respect to y
	    for (int l = 0; l <= order_quad; l++)
	      {
		numero = NumQuad2D(i,l);
		numpsi = 2*numero+1; // vectorial dof number (unknown v)
		
		// Rh -> no weight
		if (abs(lob_quad.GradPhi(j, l)) >= 1e4*epsilon_machine)
                  Rh.AddInteraction(num_loc, numpsi, lob_quad.GradPhi(j,l));
	      }
	    
	    // derivative with respect to x
	    for (int k = 0; k <= order_quad; k++)
	      {
		numero = NumQuad2D(k,j);
		numpsi = 2*numero; // vectorial dof number (unknown v)
		
		// Rh -> no weight
		if (abs(lob_quad.GradPhi(i,k)) >= 1e4*epsilon_machine)
                  Rh.AddInteraction(num_loc, numpsi, lob_quad.GradPhi(i,k));
                
	      }
	  }
      }
    
    // conversion to Column Sparse Row format
    // to have more efficient computations
    Seldon::Copy(Rh, rh_loc);
  }
  

  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleDgGauss::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    int order_geom = this->GetGeometryOrder();
    Vector2 w((order+1)*(order_geom+1));
    w.Fill(0);
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          w((order+1)*i+j) += G_GL(k, i)*Un(NumFct2D(k, j));
    
    const Matrix<int>& NumNodes2D = this->GetNumNodes2D();
    Unode.Reallocate((order_geom+1)*(order_geom+1));
    Unode.Fill(0);
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
        for (int k = 0; k <= order; k++)
          Unode(NumNodes2D(i, j)) += G_GL(k, j)*w((order+1)*i+k);
    
  }
  
  
  //! Integration against basis functions on an edge
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc edge number
  */  
  template<class Vector1,class Vector2>
  void QuadrangleDgGauss::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  //! x is overwritten by M x where M is the mass matrix
  template<class Vector2>
  void QuadrangleDgGauss::MltMassMatrixGen(Vector2& x) const
  {
    for (int i = 0; i < this->nb_dof_loc; i++)
      x(i) *= weights2d(i);
  }
  
  
  //! x is overwritten by M^-1 x where M is the mass matrix
  template<class Vector2>
  void QuadrangleDgGauss::SolveMassMatrixGen(Vector2& x) const
  {
    for (int i = 0; i < invWeights2d.GetM(); i++)
      x(i) *= invWeights2d(i);
  }
  
  
  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix
  template<class Vector1>
  void QuadrangleDgGauss::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    for (int i = 0; i < invSqrtWeights2d.GetM(); i++)
      x(i) *= invSqrtWeights2d(i);
  }
  
  
  //! computation of u on nodal points of an edge
  /*!
    \param[in] feval dof components of u
    \param[out] Unode values of u on quadrature points of the edge
    \param[in] num_loc local position of the edge in the element
  */
  template<class Vector1, class Vector2>
  void QuadrangleDgGauss
  ::ComputeValueBoundaryGen(const Vector1& feval, Vector2 & Unode, int num_loc) const
  {
    Vector1 Udof(order+1); Udof.Fill(0);
    switch (num_loc)
      {
      case 0 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
	      Udof(i) += ValGauss_Extremity(j,0)*feval(num_dof);
	    }
        break;
      case 1 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
	      Udof(j) += ValGauss_Extremity(i,1)*feval(num_dof);
	    }
        break;
      case 2 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
	      Udof(order-i) += ValGauss_Extremity(j,1)*feval(num_dof);
	    }
        break;
      case 3 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
	      Udof(order-j) += ValGauss_Extremity(i,0)*feval(num_dof);
	    }
        break;
      }
    
    int order_geom = this->GetGeometryOrder();
    Unode.Reallocate(order_geom+1);
    for (int j = 0; j <= order_geom; j++)
      {
        Unode(j) = 0;
        for (int k = 0; k <= order; k++)
          Unode(j) += G_Geom(k, j)*Udof(k);
      }
    
  }
  

  //! computation of grad u on nodal points of an edge
  /*!
    \param[in] feval dof components of u
    \param[out] Unode gradient of u on quadrature points of the edge
    \param[in] num_loc local position of the edge in the element
  */  
  template<class Vector1, class Vector2>
  void QuadrangleDgGauss
  ::ComputeGradientBoundaryGen(const Vector1& feval, Vector2 & Unode, int num_loc) const
  {
    int order_geom = this->GetGeometryOrder();
    Vector1 Udof(2*(order+1)); Udof.Fill(0);
    Unode.Reallocate(2*(order_geom+1));
    switch (num_loc)
      {
      case 0 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
	      Udof(2*i) += ValGauss_Extremity(j, 0)*feval(num_dof);
              Udof(2*i+1) += dG_Geom(j, 0)*feval(num_dof);
	    }
        
        for (int j = 0; j <= order_geom; j++)
          {
            Unode(2*j) = 0; Unode(2*j+1) = 0;
            for (int k = 0; k <= order; k++)
              {
                Unode(2*j) += dG_Geom(k, j)*Udof(2*k);
                Unode(2*j+1) += G_Geom(k, j)*Udof(2*k+1);
              }
          }
        
        break;
      case 1 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
              Udof(2*j) += dG_Geom(i, order_geom)*feval(num_dof);
	      Udof(2*j+1) += ValGauss_Extremity(i, 1)*feval(num_dof);
	    }
        
        for (int j = 0; j <= order_geom; j++)
          {
            Unode(2*j) = 0; Unode(2*j+1) = 0;
            for (int k = 0; k <= order; k++)
              {
                Unode(2*j) += G_Geom(k, j)*Udof(2*k);
                Unode(2*j+1) += dG_Geom(k, j)*Udof(2*k+1);
              }
          }

        break;
      case 2 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(order-i,j);
	      Udof(2*i) += ValGauss_Extremity(j,1)*feval(num_dof);
              Udof(2*i+1) += dG_Geom(j, order_geom)*feval(num_dof);
	    }

        for (int j = 0; j <= order_geom; j++)
          {
            Unode(2*j) = 0; Unode(2*j+1) = 0;
            for (int k = 0; k <= order; k++)
              {
                Unode(2*j) -= dG_Geom(k, j)*Udof(2*k);
                Unode(2*j+1) += G_Geom(k, j)*Udof(2*k+1);
              }
          }
        break;
      case 3 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      int num_dof = NumFct2D(i,j);
	      Udof(2*(order-j)) += dG_Geom(i, 0)*feval(num_dof);
              Udof(2*(order-j)+1) += ValGauss_Extremity(i,0)*feval(num_dof);
	    }

        for (int j = 0; j <= order_geom; j++)
          {
            Unode(2*j) = 0; Unode(2*j+1) = 0;
            for (int k = 0; k <= order; k++)
              {
                Unode(2*j) += G_Geom(k, j)*Udof(2*k);
                Unode(2*j+1) -= dG_Geom(k, j)*Udof(2*k+1);
              }
          }
        break;
      }
        
  }


  //! computation of u on quadrature points of an edge
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>  
  void QuadrangleDgGauss
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {    
    if ((r != 0) && (r != this->order_quad)) 
      {
	Vector2 Vquad(order+1);
	ElementReference<Dimension2, 1>::ApplyShTranspose(num_loc, Uh, Vquad);
	Vh.Reallocate(r+1); Vh.Fill(0);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= order; j++)
	    Vh(i) += ProjOperatorOrder(r)(j, i)*Vquad(j);
      }    
    else
      ElementReference<Dimension2, 1>::ApplyShTranspose(num_loc, Uh, Vh);
  }  

  
  //! integration against basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uquad vector containing \omega_k f(\xi_k) 
    \param[out] Vh res_i = res_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the edge
  */
  template<class T0, class Vector1, class Vector2>
  void QuadrangleDgGauss
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uquad, Vector2& Vh, int r) const
  {
    Vector1 Uh;
    if ((r != 0) && (r != this->order_quad))
      {
	Uh.Reallocate(order+1); Uh.Fill(0);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= order; j++)
	    Uh(j) += ProjOperatorOrder(r)(j, i)*Uquad(i);
      }
    else
      {
        if (this->order != this->order_quad)
          return ElementReference<Dimension2, 1>::ApplyShGen(alpha, num_loc, Uquad, Vh, r);
        
        Uh = Uquad;
      }
    
    typename Vector1::value_type vloc;
    switch (num_loc)
      {
      case 0 :	  
	for (int i = 0; i <= order; i++)
	  {
	    vloc = alpha*Uh(i);
	    for (int j = 0; j <= order; j++)
	      Vh(NumFct2D(i,j)) += ValGauss_Extremity(j, 0)*vloc;
	  }
	break;
      case 1 :
	for (int j = 0; j <= order; j++)
	  {
	    vloc = alpha*Uh(j);
	    for (int i = 0; i <= order; i++)
	      Vh(NumFct2D(i,j)) += ValGauss_Extremity(i, 1)*vloc;
	  }
        break;
      case 2 :	  
	for (int i = 0; i <= order; i++)
	  {
	    vloc = alpha*Uh(order-i);
	    for (int j = 0; j <= order; j++)
	      Vh(NumFct2D(i,j)) += ValGauss_Extremity(j, 1)*vloc;
	  }
	break;
      case 3 :
	for (int j = 0; j <= order; j++)
	  {
	    vloc = alpha*Uh(order-j);
	    for (int i = 0; i <= order; i++)
	      Vh(NumFct2D(i,j)) += ValGauss_Extremity(i, 0)*vloc;
	  }
      }
  }
  
  
  //! computation of gradient of u on quadrature points of an edge
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>  
  void QuadrangleDgGauss
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {    
    if (this->order != this->order_quad)
      return ElementReference<Dimension2, 1>::ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);

    Vector2 Vquad(2*(order+1));
    Vquad.Fill(0);
    
    int order_geom = this->GetGeometryOrder();
    typename Vector2::value_type vx, vy;
    switch (num_loc)
      {
      case 0 :
        {
          for (int i = 0; i <= order; i++)
            {
              vx = 0; vy = 0;
              for (int j = 0; j <= order; j++)
                {
                  vx += G_Geom(j, 0)*Uh(NumFct2D(i, j));
                  vy += dG_Geom(j, 0)*Uh(NumFct2D(i, j));
                }
              
              Vquad(2*i+1) = vy;
              for (int j = 0; j <= order; j++)
                Vquad(2*j) += lob_basis.GradPhi(i, j)*vx;
            }
        }
        break;
      case 1 :
        {
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              for (int i = 0; i <= order; i++)
                {
                  vx += dG_Geom(i, order_geom)*Uh(NumFct2D(i, j));
                  vy += G_Geom(i, order_geom)*Uh(NumFct2D(i, j));
                }
              
              Vquad(2*j) = vx;
              for (int i = 0; i <= order; i++)
                Vquad(2*i+1) += lob_basis.GradPhi(j, i)*vy;
            }
        }
        break;
      case 2 :
        {
          for (int i = 0; i <= order; i++)
            {
              vx = 0; vy = 0;
              for (int j = 0; j <= order; j++)
                {
                  vx += G_Geom(j, order_geom)*Uh(NumFct2D(i, j));
                  vy += dG_Geom(j, order_geom)*Uh(NumFct2D(i, j));
                }
              
              Vquad(2*(order-i)+1) = vy;
              for (int j = 0; j <= order; j++)
                Vquad(2*(order-j)) += lob_basis.GradPhi(i, j)*vx;
            }
        }
        break;
      case 3 :
        {
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              for (int i = 0; i <= order; i++)
                {
                  vx += dG_Geom(i, 0)*Uh(NumFct2D(i, j));
                  vy += G_Geom(i, 0)*Uh(NumFct2D(i, j));
                }
              
              Vquad(2*(order-j)) = vx;
              for (int i = 0; i <= order; i++)
                Vquad(2*(order-i)+1) += lob_basis.GradPhi(j, i)*vy;
            }
        }
        break;
      }    
    
    if ((r != 0) && (r != this->order)) 
      {
	Vh.Fill(0);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= order; j++)
	    {
              Vh(2*i) += ProjOperatorOrder(r)(j, i)*Vquad(2*j);
              Vh(2*i+1) += ProjOperatorOrder(r)(j, i)*Vquad(2*j+1);
            }
      }    
    else
      Copy(Vquad, Vh);
  }  

  
  //! integration against gradient of basis functions on an edge
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc edge number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is the quadrature rule of edge num_loc
  */
  template<class T0, class Vector1, class Vector2>
  void QuadrangleDgGauss
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (this->order != this->order_quad)
      return ElementReference<Dimension2, 1>::ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
    
    Vector2 Vquad(2*(order+1));
    Vquad.Fill(0);

    if ((r != 0) && (r != this->order_quad)) 
      {
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= order; j++)
	    {
              Vquad(2*j) += ProjOperatorOrder(r)(j, i)*Uh(2*i);
              Vquad(2*j+1) += ProjOperatorOrder(r)(j, i)*Uh(2*i+1);
            }
      }    
    else
      Copy(Uh, Vquad);
    
    Mlt(alpha, Vquad);

    int order_geom = this->GetGeometryOrder();
    typename Vector2::value_type vx, vy;
    switch (num_loc)
      {
      case 0 :
        {
          for (int i = 0; i <= order; i++)
            {
              vx = 0; vy = 0;
              vy = Vquad(2*i+1);
              for (int j = 0; j <= order; j++)
                vx += lob_basis.GradPhi(i, j)*Vquad(2*j);
              
              for (int j = 0; j <= order; j++)
                {
                  Vh(NumFct2D(i, j)) += G_Geom(j, 0)*vx;
                  Vh(NumFct2D(i, j)) += dG_Geom(j, 0)*vy;
                }
            }
        }
        break;
      case 1 :
        {
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              vx = Vquad(2*j);
              for (int i = 0; i <= order; i++)
                vy += lob_basis.GradPhi(j, i)*Vquad(2*i+1);
              
              for (int i = 0; i <= order; i++)
                {
                  Vh(NumFct2D(i, j)) += dG_Geom(i, order_geom)*vx;
                  Vh(NumFct2D(i, j)) += G_Geom(i, order_geom)*vy;
                }
            }
        }
        break;
      case 2 :
        {
          for (int i = 0; i <= order; i++)
            {
              vx = 0; vy = 0;
              vy = Vquad(2*(order-i)+1);
              for (int j = 0; j <= order; j++)
                vx += lob_basis.GradPhi(i, j)*Vquad(2*(order-j));
              
              for (int j = 0; j <= order; j++)
                {
                  Vh(NumFct2D(i, j)) += G_Geom(j, order_geom)*vx;
                  Vh(NumFct2D(i, j)) += dG_Geom(j, order_geom)*vy;
                }
            }
        }
        break;
      case 3 :
        {
          for (int j = 0; j <= order; j++)
            {
              vx = 0; vy = 0;
              vx = Vquad(2*(order-j));
              for (int i = 0; i <= order; i++)
                vy += lob_basis.GradPhi(j, i)*Vquad(2*(order-i)+1);
              
              for (int i = 0; i <= order; i++)
                {
                  Vh(NumFct2D(i, j)) += dG_Geom(i, 0)*vx;
                  Vh(NumFct2D(i, j)) += G_Geom(i, 0)*vy;
                }
            }
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
  void QuadrangleDgGauss::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    T zero; SetComplexZero(zero);
    if (mass == zero)
      return;
    
    for (int i = 0; i < nb_dof_loc; i++)
      A.AddInteraction(m+i, n+i, mass*weights2d(i));
  }
  
  
  //! we add constant elementary matrix to the matrix A
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + mass M + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p + \sum_p  C(p, q) S^{p,q}
    where M is the mass matrix equal to
    M_{i, j} = \int_K \phi_j \phi_i dx
    R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
    S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void QuadrangleDgGauss::
  AddConstantElemMatrixGen(int m, int n, const T& mass,
			   const TinyMatrix<T, Prop, 2, 2>& C,
			   const TinyVector<T, 2>& D, 
			   const TinyVector<T, 2>& E,
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
        {
          int i = NumFct2D(i1, i2);
	  if (!null_term(0))
	    A.AddInteraction(m+i, n+i, mass*weights2d(i));
	  
          for (int k = 0; k <= order; k++)
            {
	      if (!null_term(2))
		{
		  A.AddInteraction(i+m, n+NumFct2D(k, i2),
				   D(0)*lob_basis.GradPhi(i1, k)
				   *lob_basis.Weights(i2)*lob_basis.Weights(k));
		  
		  A.AddInteraction(i+m, n+NumFct2D(i1, k),
				   D(1)*lob_basis.GradPhi(i2, k)
				   *lob_basis.Weights(i1)*lob_basis.Weights(k));
		}
	      
              if (!null_term(3))
                {
		  A.AddInteraction(m+NumFct2D(k, i2), n+i,
				   E(0)*lob_basis.GradPhi(i1, k)
				   *lob_basis.Weights(i2)*lob_basis.Weights(k));
		  
		  A.AddInteraction(m+NumFct2D(i1, k), n+i,
				   E(1)*lob_basis.GradPhi(i2, k)
				   *lob_basis.Weights(i1)*lob_basis.Weights(k));
		}
            }
	  
	  if (null_term(1))
	    continue;
	  
          for (int j1 = 0; j1 <= order; j1++)
	    A.AddInteraction(m+i, n+NumFct2D(j1, i2), C(0,0)*stiffness_matrix1D(i1, j1)*lob_basis.Weights(i2));
          
          for (int j2 = 0; j2 <= order; j2++)
	    A.AddInteraction(m+i, n+NumFct2D(i1, j2),  C(1,1)*stiffness_matrix1D(i2, j2)*lob_basis.Weights(i1));
              
          for (int j1 = 0; j1 <= order; j1++)
            for (int j2 = 0; j2 <= order; j2++)
              {
                int j = NumFct2D(j1, j2);
		A.AddInteraction(m+i, n+j, C(0,1)*lob_basis.GradPhi(i1, j1)*lob_basis.GradPhi(j2, i2)
				 *lob_basis.Weights(j1)*lob_basis.Weights(i2)
				 + C(1,0)*lob_basis.GradPhi(j1, i1)*lob_basis.GradPhi(i2, j2)
				 *lob_basis.Weights(i1)*lob_basis.Weights(j2));
              }
        }
  }
  

  //! adds variable elementary matrix to mat
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + R + Rt + S
    M_ij = \int mass \varphi_j \varphi_i
    R_ij = \int D \cdot \nabla \varphi_i \varphi_j
    Rt_ij = \int E \cdot \nabla \varphi_j \varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
   */
  template<class T, class Prop>
  void QuadrangleDgGauss
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
			     const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			     const Vector<TinyVector<T, 2> >& D,
			     const Vector<TinyVector<T, 2> >& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension2, 1>::
      AddVariableElemMatrixOpt(off_row, off_col, mass, C, D, E, null_term, mat);
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
  void QuadrangleDgGauss::AddVariableMassMatrixGen(int off_row, int off_col,
						   const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    if (order != order_quad)
      return ElementReference<Dimension2, 1>::
        AddVariableMassMatrixOpt(off_row, off_col, A, mat);
    
    for (int i = 0; i < nb_dof_loc; i++)
      mat.AddInteraction(i+off_row, i+off_col, A(i));
  }


  //! retrieves values of a single basis functions on all quadrature points
  void QuadrangleDgGauss::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    if (order != order_quad)
      return ElementReference<Dimension2, 1>::GetValueSinglePhiQuadrature(i, phi);

    phi.Zero();
    phi(i) = Real_wp(1);
  }


  //! retrieves values of a single basis functions on all quadrature points
  void QuadrangleDgGauss
  ::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR2& grad_phi) const
  {
    if (order != order_quad)
      return ElementReference<Dimension2, 1>::GetGradientSinglePhiQuadrature(i, phi, grad_phi);
    
    phi.Zero();
    grad_phi.Zero();
    phi(i) = Real_wp(1);
    int i1 = CoordinateDofs(i, 0);
    int i2 = CoordinateDofs(i, 1);
    for (int k = 0; k <= order; k++)
      {
	grad_phi(NumQuad2D(k, i2))(0) = lob_basis.GradPhi(i1, k);
	grad_phi(NumQuad2D(i1, k))(1) = lob_basis.GradPhi(i2, k);
      }
  }


  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void QuadrangleDgGauss::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    if (order != order_quad)
      return ElementReference<Dimension2, 1>::GetValuePhiOnQuadraturePoint(k, phi);

    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    if (k < nb_points_quadrature_inside)
      {
        phi(k) = 1.0;
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/(order+1);
        int k2 = (k-nb_points_quadrature_inside)%(order+1);
        switch (num_loc)
          {
          case 0:
            {
              for (int j = 0; j <= order; j++)
                phi(NumFct2D(k2, j)) = ValGauss_Extremity(j, 0);
            }
            break;
          case 1:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(i, k2)) = ValGauss_Extremity(i, 1);
            }
            break;
          case 2:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(order-k2, i)) = ValGauss_Extremity(i, 1);
            }
            break;
          case 3:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(i, order-k2)) = ValGauss_Extremity(i, 0);
            }
            break;
          }
      }
  }
    
}
  
#define MONTJOIE_FILE_QUADRANGLE_DG_GAUSS_CXX
#endif
