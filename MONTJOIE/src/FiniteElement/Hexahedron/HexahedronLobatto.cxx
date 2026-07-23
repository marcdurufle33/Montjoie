#ifndef MONTJOIE_FILE_HEXAHEDRON_LOBATTO_CXX

namespace Montjoie
{

  //! default constructor
  HexahedronLobatto::HexahedronLobatto() : HexahedronGauss()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = true;
    
    mass_lumping = true;
    diagonal_mass = true;
    use_quadrature_for_sh = false;
    use_quadrature_for_rh = false;
    optimized_elem_matrix = true;

    type_integration_edge = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  
  
  //! constructing finite element
  void HexahedronLobatto::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf_tri, int rsurf_quad,
						 int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    HexahedronGauss::ConstructFiniteElement(r, rgeom, r, Globatto<Real_wp>::QUADRATURE_LOBATTO);

    if (this->GetGeometryOrder() == rquad)
      this->Fb_geom.quadrature_equal_nodal = true;
    else
      this->Fb_geom.quadrature_equal_nodal = false;
        
    points3d_all.Reallocate(6*(order+1)*(order+1) + (order+1)*(order+1)*(order+1));
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points3d_all(i) = this->PointsND(i);
    
    int nb = nb_points_quadrature_inside;
    for (int num_loc = 0; num_loc < 6; num_loc++)
      {
        for (int k = 0; k < (order+1)*(order+1); k++)
          points3d_all(nb+k) = this->PointsND(this->num_quad_points_surf(num_loc)(k));
        
        nb += (order+1)*(order+1);
      }

#ifdef MONTJOIE_WITH_META_PROGRAMMING
    switch (r)
      {
      case 1:
	elt_r1.ConstructFiniteElement(*this);
	break;
      case 2:
	elt_r2.ConstructFiniteElement(*this);
	break;
      case 3:
	elt_r3.ConstructFiniteElement(*this);
	break;
      case 4:
	elt_r4.ConstructFiniteElement(*this);
	break;
      case 5:
	elt_r5.ConstructFiniteElement(*this);
	break;
      }
#endif

  }
  
  
  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
   */  
  template<class Vector1,class Vector2>
  void HexahedronLobatto::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    res.Fill(0);
    HexahedronGauss::ApplyShGen(1.0, num_loc, feval, res);
  }
  
  
  //! computation of grad u on nodal points of a face
  /*!
    \param[in] feval dof components of u
    \param[out] res gradient of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */  
  template<class Vector1, class Vector2>
  void HexahedronLobatto::
  ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    res.Reallocate(3*(order+1)*(order+1));
    int node = 0; int ext = 0;
    if (num_loc >= 3)
      ext = order;
    
    if ((num_loc == 0)||(num_loc == 5))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    node = NumDofs2D(i,j);
	    res(3*node) = 0; res(3*node+1) = 0; res(3*node+2) = 0;
	    for (int k = 0; k <= order; k++)
              {
                res(3*node) += lob_basis.GradPhi(k,ext)*feval(NumDofs3D(k,i,j));
                res(3*node+1) += lob_basis.GradPhi(k,i)*feval(NumDofs3D(ext,k,j));
                res(3*node+2) += lob_basis.GradPhi(k,j)*feval(NumDofs3D(ext,i,k));
              }
	  }
    else if ((num_loc == 1)||(num_loc == 4))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    node = NumDofs2D(i,j);
            res(3*node) = 0; res(3*node+1) = 0; res(3*node+2) = 0;
	    for (int k = 0; k <= order; k++)
              {
                res(3*node) += lob_basis.GradPhi(k,i)*feval(NumDofs3D(k,ext,j));
                res(3*node+1) += lob_basis.GradPhi(k,ext)*feval(NumDofs3D(i,k,j));
                res(3*node+2) += lob_basis.GradPhi(k,j)*feval(NumDofs3D(i,ext,k));
              }
	  }
    else if ((num_loc == 2)||(num_loc == 3))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    node = NumDofs2D(i,j);
            res(3*node) = 0; res(3*node+1) = 0; res(3*node+2) = 0;
	    for (int k = 0; k <= order; k++)
              {
                res(3*node) += lob_basis.GradPhi(k,i)*feval(NumDofs3D(k,j,ext));
                res(3*node+1) += lob_basis.GradPhi(k,j)*feval(NumDofs3D(i,k,ext));
                res(3*node+2) += lob_basis.GradPhi(k,ext)*feval(NumDofs3D(i,j,k));
              }
	  }
    
  }

  
  //! x is overwritten by M^-1 x where M is the mass matrix
  template<class Vector1>
  void HexahedronLobatto::SolveMassMatrixGen(Vector1& Uh) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      Uh(i) *= invWeights3d(i);
  }
   
  
  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix
  template<class Vector1>
  void HexahedronLobatto::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& Uh) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      Uh(i) *= invSqrtWeights3d(i);
  }
  

  //! Vh is overwritten by M Vh where M is the mass matrix  
  template<class Vector2>
  void HexahedronLobatto::MltMassMatrixGen(Vector2& Vh) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      Vh(i) *= this->WeightsND(i);
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
  void HexahedronLobatto::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    // mass lumping => Ch = I
    Vh = Uh;
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
  void HexahedronLobatto::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    // mass lumping => Ch = I
    Vh = Uh;
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
  void HexahedronLobatto::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    switch (order)
      {
      case 1 :
        elt_r1.ApplyRh(Uh, Vh);
        return;
      case 2 :
        elt_r2.ApplyRh(Uh, Vh);
        return;
      case 3 :
        elt_r3.ApplyRh(Uh, Vh);
        return;
      case 4 :
        elt_r4.ApplyRh(Uh, Vh);
        return;
      case 5 :
        elt_r5.ApplyRh(Uh, Vh);
        return;
      }
#endif
    
    Mlt(rh_loc, Uh, Vh);
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
  void HexahedronLobatto::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    switch (order)
      {
      case 1 :
        elt_r1.ApplyRhTranspose(Uh, Vh);
        return;
      case 2 :
        elt_r2.ApplyRhTranspose(Uh, Vh);
        return;
      case 3 :
        elt_r3.ApplyRhTranspose(Uh, Vh);
        return;
      case 4 :
        elt_r4.ApplyRhTranspose(Uh, Vh);
        return;
      case 5 :
        elt_r5.ApplyRhTranspose(Uh, Vh);
        return;
      } 
#endif

    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }


  //! transpose of ApplyRhBoundaryTransposeGen
  template<class Vector1, class Vector2>
  void HexahedronLobatto::ApplyRhBoundaryGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(rh_boundary, Uh, Vh);
  }


  //! computes v = grad u only for boundary dofs
  template<class Vector1, class Vector2>
  void HexahedronLobatto::ApplyRhBoundaryTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_boundary, Uh, Vh);
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
  void HexahedronLobatto::
  ApplyRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    switch (order)
      {
      case 1 :
        elt_r1.ApplyRhSplit(Uh, Ux, Uy, Uz);
        return;
      case 2 :
        elt_r2.ApplyRhSplit(Uh, Ux, Uy, Uz);
        return;
      case 3 :
        elt_r3.ApplyRhSplit(Uh, Ux, Uy, Uz);
        return;
      case 4 :
        elt_r4.ApplyRhSplit(Uh, Ux, Uy, Uz);
        return;
      case 5 :
        elt_r5.ApplyRhSplit(Uh, Ux, Uy, Uz);
        return;
      }
#endif

    Mlt(rh_locX, Uh, Ux);
    Mlt(rh_locY, Uh, Uy);
    Mlt(rh_locZ, Uh, Uz);
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
  void HexahedronLobatto
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    int Nquad = (this->order_quad+1)*(this->order_quad+1);
    if ((r == 0) || (r == order))
      {
	for (int i = 0; i < Nquad; i++)
	  Vh(i) = Uh(FacesDof(i, num_loc));    
      }
    else
      {
	Vector2 Vlob(Nquad);
	for (int i = 0; i < Nquad; i++)
	  Vlob(i) = Uh(FacesDof(i, num_loc));
	
	Vh.Fill(0);
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vloc = 0;
	      for (int j = 0; j <= order; j++)
		vloc += Ah(j, m)*Vlob(NumDofs2D(i, j));
	      
	      for (int n = 0; n <= r; n++)
		Vh(Num(n, m)) += Ah(i, n)*vloc;
	    }
      }
  }
  
  
  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronLobatto::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    int Nquad = (order+1)*(order+1);
    if ((r == 0) || (r == order))
      {
	int num_dof;
	for (int i = 0; i < Nquad; i++)
	  {
	    num_dof = FacesDof(i, num_loc);
	    Vh(num_dof) += alpha*Uh(i);
	  }
      }
    else
      {
	Vector1 Vlob(Nquad);
	Vlob.Fill(0);
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vloc = 0;
	      for (int n = 0; n <= r; n++)
		vloc += Ah(i, n)*Uh(Num(n, m));
	      
	      for (int j = 0; j <= order; j++)
		Vlob(NumDofs2D(i, j)) += Ah(j, m)*vloc;
	    }
	
	int num_dof;
	for (int i = 0; i < Nquad; i++)
	  {
	    num_dof = FacesDof(i, num_loc);
	    Vh(num_dof) += alpha*Vlob(i);
	  }
      }
  }
  
  
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] res components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronLobatto
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& res, Vector2& Vh, int r) const
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
                Udof(p) = res(NumDofs3D(ext, j, k));
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(i, j, k);
                    feval(3*p) += dGL_Geom(i, ext)*res(node);
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
                Udof(p) = res(NumDofs3D(j, ext, k));
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, i, k);
                    feval(3*p+1) += dGL_Geom(i, ext)*res(node);                    
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
                Udof(p) += res(NumDofs3D(j, k, ext));
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, k, i);
                    feval(3*p+2) += dGL_Geom(i, ext)*res(node);                    
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

    if ((r == 0) || (r == order))
      {
        Copy(feval, Vh);
      }
    else
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
  }
  
  
  //! integration against gradient of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] res res_i = res_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronLobatto::
  ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& res, int r) const
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
                res(NumDofs3D(ext, j, k)) += Udof(p);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(i, j, k);
                    res(node) += dGL_Geom(i, ext)*feval(3*p);
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
                res(NumDofs3D(j, ext, k)) += Udof(p);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, i, k);
                    res(node) += dGL_Geom(i, ext)*feval(3*p+1);
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
                res(NumDofs3D(j, k, ext)) += Udof(p);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, k, i);
                    res(node) += dGL_Geom(i, ext)*feval(3*p+2);
                  }
              }
        }
        break;
      }
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then y = Rh x
   */
  template<class Vector1, class Vector2>
  void HexahedronLobatto::ApplyConstantRhGen(const Vector1& x, Vector2& y) const
  {
    Mlt(rh_weight, x, y);
  }
  
    
  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    (Rh^3)_{i, j} = \int_K \varphi_j d/dz(\varphi_i) dx
    then Ux = Rh^1 Uh, Uy = Rh^2 Uh, Uz = Rh^3 Uh
   */
  template<class Vector1, class Vector2>
  void HexahedronLobatto::
  ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
    Mlt(rh_weightX, Uh, Ux);
    Mlt(rh_weightY, Uh, Uy);
    Mlt(rh_weightZ, Uh, Uz);
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then y = Rh* x
   */  
  template<class Vector1, class Vector2>
  void HexahedronLobatto::ApplyConstantRhTransposeGen(const Vector1& x, Vector2& y) const
  {
    Mlt(SeldonTrans, rh_weight, x, y);
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
  void HexahedronLobatto::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    T zero; SetComplexZero(zero);
    const VectReal_wp& weights3d = this->WeightsND();
    if (mass != zero)
      for (int i = 0; i < nb_dof_loc; i++)
	A.AddInteraction(m+i, n+i, mass*weights3d(i));    
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
  void HexahedronLobatto
  ::AddConstantElemMatrixGen(int m, int n, const T& mass, const TinyMatrix<T, Prop, 3, 3>& C,
			     const TinyVector<T, 3>& D, const TinyVector<T, 3>& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;
    
    T vloc;
    const VectReal_wp& weights3d = this->WeightsND();
    const VectReal_wp& weights1d = this->Weights1D();    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
        for (int i3 = 0; i3 <= order; i3++)
          {
	    int i = NumDofs3D(i1, i2, i3);
	    val.Zero();
	    
	    // mass term
	    if (!null_term(0))
	      {
		vloc = mass*weights3d(i);
		val(i) = vloc;
	      }

	    if (!null_term(2))
	      for (int k = 0; k <= order; k++)
		{
		  // gradient terms
		  int j = NumDofs3D(k, i2, i3);
		  vloc = D(0)*lob_basis.GradPhi(i1, k)*weights3d(j);
		  
		  val(j) += vloc;		
		  
		  j = NumDofs3D(i1, k, i3);
		  vloc = D(1)*lob_basis.GradPhi(i2, k)*weights3d(j);
		  
		  val(j) += vloc;
		  
		  j = NumDofs3D(i1, i2, k);
		  vloc = D(2)*lob_basis.GradPhi(i3, k)*weights3d(j);
		  
		  val(j) += vloc;
		}

	    if (!null_term(3))
	      for (int k = 0; k <= order; k++)
		{
		  // gradient terms
		  int j = NumDofs3D(k, i2, i3);
		  vloc = E(0)*lob_basis.GradPhi(k, i1)*weights3d(i);
		  
		  val(j) += vloc;		
		  
		  j = NumDofs3D(i1, k, i3);
		  vloc = E(1)*lob_basis.GradPhi(k, i2)*weights3d(i);
		  
		  val(j) += vloc;
		  
		  j = NumDofs3D(i1, i2, k);
		  vloc = E(2)*lob_basis.GradPhi(k, i3)*weights3d(i);
		  
		  val(j) += vloc;
		}
	    	    
	    if (!null_term(1))
	      {
		for (int j1 = 0; j1 <= order; j1++)
		  val(NumDofs3D(j1, i2, i3))
		    += C(0,0)*stiffness_matrix1D(i1, j1)*weights1d(i2)*weights1d(i3);
		
		for (int j2 = 0; j2 <= order; j2++)
		  val(NumDofs3D(i1, j2, i3))
		    += C(1,1)*stiffness_matrix1D(i2, j2)*weights1d(i1)*weights1d(i3);
		
		for (int j3 = 0; j3 <= order; j3++)
		  val(NumDofs3D(i1, i2, j3))
		    += C(2,2)*stiffness_matrix1D(i3, j3)*weights1d(i1)*weights1d(i2);
		
		for (int j1 = 0; j1 <= order; j1++)
		  for (int j2 = 0; j2 <= order; j2++)
		    {
		      int j = NumDofs3D(j1, j2, i3);
                      val(j) += C(0,1)*lob_basis.GradPhi(i1, j1)
                        *lob_basis.GradPhi(j2, i2)*weights1d(j1)*weights1d(i2)*weights1d(i3);
                      
                      val(j) += C(1,0)*lob_basis.GradPhi(j1, i1)
                        *lob_basis.GradPhi(i2, j2)*weights1d(i1)*weights1d(j2)*weights1d(i3);
		      
		      j = NumDofs3D(j1, i2, j2);
                      val(j) += C(0,2)*lob_basis.GradPhi(i1, j1)
                        *lob_basis.GradPhi(j2, i3)*weights1d(j1)*weights1d(i3)*weights1d(i2);
                      
                      val(j) += C(2,0)*lob_basis.GradPhi(j1, i1)
                        *lob_basis.GradPhi(i3, j2)*weights1d(i1)*weights1d(j2)*weights1d(i2);
                  
		      j = NumDofs3D(i1, j1, j2);
		      val(j) += C(1,2)*lob_basis.GradPhi(i2, j1)
			*lob_basis.GradPhi(j2, i3)*weights1d(j1)*weights1d(i3)*weights1d(i1);
		  
		      val(j) += C(2,1)*lob_basis.GradPhi(j1, i2)
			*lob_basis.GradPhi(i3, j2)*weights1d(i2)*weights1d(j2)*weights1d(i1);
                    }
	      }
	    
	    A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
	  }
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
  void HexahedronLobatto::
  AddVariableMassMatrixGen(int off_row, int off_col,
			   const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      mat.AddInteraction(i+off_row, i+off_col, A(i));
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
  void HexahedronLobatto
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& A,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			     const Vector<TinyVector<T, 3> >& D,
			     const Vector<TinyVector<T, 3> >& E,
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc, zero;
    SetComplexZero(zero);
    int num_loc, num_loc_2;
    const Matrix<Real_wp>& dGL_GL = lob_basis.GradPhi();
    
    // Loop on the local dof
    // we use a tensorial numbering
    for (int l = 0; l <= order; l++)
      for (int m = 0; m <= order; m++)
	for (int n = 0; n <= order; n++)
	  {
	    // row
	    num_loc_2 = NumDofs3D(l, m, n);
	    val.Zero();
	    
	    if (!null_term(0))
	      val(num_loc_2) = A(num_loc_2);
	    
	    // Loop on the local dof
	    // we use a tensorial numbering
	    for (int i = 0; i <= order; i++)
	      for (int j = 0; j <= order; j++)
		for (int k = 0; k <= order; k++)
		  {
		    // column
		    num_loc = NumDofs3D(i, j, k);
		    
		    // Now we compute the interaction between num_loc and num_loc2
		    // in the stiffness matrix
		    vloc = zero;
		    if (i==l)
		      {
			if (k == n)
			  {
			    // du/dy dv/dy
			    if (!null_term(1))
			      for (int q = 0; q <= order; q++)
				vloc += C(NumQuad3D(l, q, n))(1, 1)*dGL_GL(j, q)*dGL_GL(m, q);
			    
			    // du/dy v
			    if (!null_term(3))
			      vloc += E(NumQuad3D(l, m, n))(1)*dGL_GL(j, m);

			    // u dv/dy
			    if (!null_term(2))
			      vloc += D(NumQuad3D(i, j, k))(1)*dGL_GL(m, j);
			  }

			if (j==m)
			  {
			    // du/dz dv/dz
			    if (!null_term(1))
			      for (int q = 0; q <= order; q++)
				vloc += C(NumQuad3D(l, m, q))(2, 2)*dGL_GL(k, q)*dGL_GL(n, q);
			    
			    // du/dz v
			    if (!null_term(3))
			      vloc += E(NumQuad3D(l, m, n))(2)*dGL_GL(k, n);

			    // u dv/dz
			    if (!null_term(2))
			      vloc += D(NumQuad3D(i, j, k))(2)*dGL_GL(n, k);
			  }
			
			if (!null_term(1))
			  {
			    // du/dy dv/dz
			    vloc += C(NumQuad3D(l, m, k))(2, 1)*dGL_GL(j, m)*dGL_GL(n, k);
			    
			    // du/dz dv/dy
			    vloc += C(NumQuad3D(l, j, n))(1, 2)*dGL_GL(m, j)*dGL_GL(k, n);
			  }			
		      } // end case i==l
		    
		    if (j==m)
		      {
			if (k==n)
			  {
			    // du/dx dv/dx
			    if (!null_term(1))
			      for (int q = 0; q <= order; q++)
				vloc += C(NumQuad3D(q, m, n))(0, 0)*dGL_GL(i, q)*dGL_GL(l, q);
			    
			    // du/dx v
			    if (!null_term(3))
			      vloc += E(NumQuad3D(l, m, n))(0)*dGL_GL(i, l);

			    // u dv/dx
			    if (!null_term(2))
			      vloc += D(NumQuad3D(i, j, k))(0)*dGL_GL(l, i);
			  }
			
			if (!null_term(1))
			  {
			    // du/dx dv/dz
			    vloc += C(NumQuad3D(l, m, k))(2, 0)*dGL_GL(i, l)*dGL_GL(n, k);
			    
			    // du/dz dv/dx
			    vloc += C(NumQuad3D(i, m, n))(0, 2)*dGL_GL(l, i)*dGL_GL(k, n);
			  }
			
		      } // end case j==m 
		    
		    if (!null_term(1))
		      if (k==n)
			{
			  // du/dx dv/dy
			  vloc += C(NumQuad3D(l, j, n))(1, 0)*dGL_GL(i, l)*dGL_GL(m, j);
			  
			  // du/dy dv/dx
			  vloc += C(NumQuad3D(i, m, n))(0, 1)*dGL_GL(l, i)*dGL_GL(j, m);
			}
		    
		    val(num_loc) += vloc;
		  } // end dof2
	    
	    mat.AddInteractionRow(off_row+num_loc_2, this->nb_dof_loc, num, val);
	  } // end dof
  }
      
    
  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */  
  void HexahedronLobatto::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    phi(k) = 1.0;
  }
  
  
  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
   */
  void HexahedronLobatto::GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    FillZero(grad_phi);
    
    int k1 = CoordinateQuad3D(k, 0);
    int k2 = CoordinateQuad3D(k, 1);
    int k3 = CoordinateQuad3D(k, 2);
    for (int i = 0; i <= order; i++)
      grad_phi(NumDofs3D(i, k2, k3))(0) = lob_basis.GradPhi(i, k1);
    
    for (int i = 0; i <= order; i++)
      grad_phi(NumDofs3D(k1, i, k3))(1) = lob_basis.GradPhi(i, k2);
    
    for (int i = 0; i <= order; i++)
      grad_phi(NumDofs3D(k1, k2, i))(2) = lob_basis.GradPhi(i, k3);
  }


  //! retrieves values of a single basis function on all quadrature points
  void HexahedronLobatto::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    phi.Zero();
    phi(i) = Real_wp(1);
  }


  //! retrieves values and gradients of a single basis function on all quadrature points
  void HexahedronLobatto
  ::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR3& grad_phi) const
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


  //! displays details about class HexahedronLobatto
  ostream& operator <<(ostream& out, const HexahedronLobatto& e)
  {
    out << static_cast<const HexahedronReference<1>&>(e);
    return out;
  }
  
}
  
#define MONTJOIE_FILE_HEXAHEDRON_LOBATTO_CXX
#endif
