#ifndef MONTJOIE_FILE_QUADRANGLE_LOBATTO_CXX

namespace Montjoie
{
  
#ifdef MONTJOIE_WITH_META_PROGRAMMING
  template<int r>
  void OptQuadrangleLobatto<r>::ConstructFiniteElement(const QuadrangleLobatto& Fb)
  {
    for (int i = 0; i < Fb.GetNbDof(); i++)
      Weights_ND(i) = Fb.WeightsND(i);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	dphi(i, j) = Fb.lob_basis.GradPhi(i,j);
  }
#endif
  
  //! default constructor
  QuadrangleLobatto::QuadrangleLobatto() : QuadrangleGauss()
  {
    optimized_elem_matrix = true;
    type_integration_edge = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    this->use_quadrature_for_rh = false;
  }
  

  size_t QuadrangleLobatto::GetMemorySize() const
  {
    size_t taille = QuadrangleGauss::GetMemorySize();
    taille += rh_weight.GetMemorySize();
    taille += rh_weightX.GetMemorySize();
    taille += rh_weightY.GetMemorySize();
    taille += points2d_all.GetMemorySize();
    return taille;
  }

  
  //! constructing finite element
  void QuadrangleLobatto::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						 int rsurf, int type_surf)
  {    
    if (rquad == 0)
      rquad = r;

    QuadrangleGauss::ConstructFiniteElement(r, rgeom, rquad, Globatto<Real_wp>::QUADRATURE_LOBATTO);

    if (this->GetGeometryOrder() == rquad)
      this->Fb_geom.quadrature_equal_nodal = true;
    else
      this->Fb_geom.quadrature_equal_nodal = false;
    
    if (this->GetGeometryOrder() == r)
      this->Fb_geom.dof_equal_nodal = true;
    else
      this->Fb_geom.dof_equal_nodal = false;
    
    if (r == rquad)
      {
        this->Fb_geom.dof_equal_quadrature = true;
        mass_lumping = true;
        diagonal_mass = true;
      }
    else
      this->Fb_geom.dof_equal_quadrature = false;
    
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    
    ConstructMassMatrix();
    this->Fb_geom.ComputeCoefficientTransformation();    
    
    if (r == rquad)
      ConstructStiffnessMatrix();

    points2d_all.Reallocate(4*(order_quad+1) + (order_quad+1)*(order_quad+1));
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      points2d_all(i) = this->PointsND(i);
    
    int nb = nb_points_quadrature_inside;
    for (int num_loc = 0; num_loc < 4; num_loc++)
      {
        for (int k = 0; k < order_quad+1; k++)
          points2d_all(nb+k) = this->PointsND(this->num_quad_points_surf(num_loc)(k));
	
        nb += order+1;
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
      case 6:
	elt_r6.ConstructFiniteElement(*this);
	break;
      case 7:
	elt_r7.ConstructFiniteElement(*this);
	break;
      case 8:
	elt_r8.ConstructFiniteElement(*this);
	break;
      }
#endif

  }
  

  void QuadrangleLobatto::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    EdgeLobattoReference* edge = new EdgeLobattoReference();
    edge->ConstructFiniteElement(r, r, rquad);
    this->element_surface = edge;
  }
  
  
  //! constructing mass matrix
  void QuadrangleLobatto::ConstructMassMatrix()
  {
  }    
  
  
  //! computation of stiffness matrix R_h
  /*!
    this matrix is used to compute gradients
   */
  void QuadrangleLobatto::ConstructStiffnessMatrix()
  {    
    // local matrices
    // matrix R_h and transpose
    Matrix<Real_wp, General, ArrayRowSparse> Rh(nb_dof_loc, 2*nb_dof_loc);
    Matrix<Real_wp, General, ArrayRowSparse> RhPoids(nb_dof_loc, 2*nb_dof_loc);
    Matrix<Real_wp, General, ArrayRowSparse> RhPoidsX(nb_dof_loc, 2*nb_dof_loc),
      RhPoidsY(nb_dof_loc, 2*nb_dof_loc);

    Matrix<Real_wp, General, ArrayRowSparse> RhB(nb_dof_loc, 2*nb_dof_boundaries);
    
    int numero, numpsi;
    // loop on all dofs of the square (i,j)
    for (int i = 0; i <= order; i++)
      {
	for (int j = 0; j <= order; j++)
	  {
	    // considered dof
	    int num_loc = NumQuad2D(i,j);
	    
	    // derivatives with respect to y
	    for (int l = 0; l <= order; l++)
	      {
		numero = NumQuad2D(i,l);
		numpsi = 2*numero+1; // vectorial dof number (unknown v)
		Real_wp poids = this->WeightsND(numero);
                
		// Rh -> no weight
		if (abs(lob_basis.GradPhi(j,l)) >= 1e3*epsilon_machine)
		  {
		    Rh.AddInteraction(num_loc, numpsi, lob_basis.GradPhi(j,l));
                    RhPoids.AddInteraction(num_loc, numpsi, poids*lob_basis.GradPhi(j,l));
                    RhPoidsY.AddInteraction(num_loc, numpsi, poids*lob_basis.GradPhi(j,l));
		    if (numpsi < 2*this->nb_dof_boundaries)
		      RhB.AddInteraction(num_loc, numpsi, lob_basis.GradPhi(j,l));
		  }
	      }
	    
	    // derivative with respect to x
	    for (int k = 0; k <= order; k++)
	      {
		numero = NumQuad2D(k,j);
		numpsi = 2*numero; // vectorial dof number (unknown v)
		Real_wp poids = this->WeightsND(numero);
                
		// Rh -> no weight
		if (abs(lob_basis.GradPhi(i,k)) >= 1e3*epsilon_machine)
		  {
		    Rh.AddInteraction(num_loc, numpsi, lob_basis.GradPhi(i,k));
                    RhPoids.AddInteraction(num_loc, numpsi, poids*lob_basis.GradPhi(i,k));
                    RhPoidsX.AddInteraction(num_loc, numpsi, poids*lob_basis.GradPhi(i,k));
		    if (numpsi < 2*this->nb_dof_boundaries)
		      RhB.AddInteraction(num_loc, numpsi, lob_basis.GradPhi(i,k));
		  }
	      }
	  }
      }
    
    // conversion to Column Sparse Row format
    // to have more efficient computations
    Seldon::Copy(Rh, rh_loc); Seldon::Copy(RhPoids, rh_weight);
    Seldon::Copy(RhPoidsX, rh_weightX); Seldon::Copy(RhPoidsY, rh_weightY);
    Seldon::Copy(RhB, rh_boundary);      
  }
  

  //! Uh is overwritten by M^-1 Uh where M is the mass matrix
  template<class Vector1>
  void QuadrangleLobatto::SolveMassMatrixGen(Vector1& Uh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::SolveMassMatrix(Uh);

    for (int i = 0; i < nb_dof_loc; i++)
      Uh(i) *= invWeights2d(i);
  }
   
  
  //! Uh is overwritten by L^-1 Uh or L^-T Uh where M = L L^t is the mass matrix  
  template<class Vector1>
  void QuadrangleLobatto::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& Uh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::SolveCholesky(TransA, Uh);

    for (int i = 0; i < nb_dof_loc; i++)
      Uh(i) *= invSqrtWeights2d(i);
  }
  
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector2>
  void QuadrangleLobatto::MltMassMatrixGen(Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::MltMassMatrix(Vh);

    const VectReal_wp& weights2d = this->WeightsND();
    for (int i = 0; i < nb_dof_loc; i++)
      Vh(i) *= weights2d(i);
  }
  
  
  //! Integration against basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh res_i = \int_K f \varphi_i dx
   */
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyCh(Uh, Vh);

    // Ch = Id (quadrature points = interpolation points)
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
  void QuadrangleLobatto::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyChTranspose(Uh, Vh);

    // Ch = Id (quadrature points = interpolation points)
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
  void QuadrangleLobatto::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyRh(Uh, Vh);

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
  void QuadrangleLobatto::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyRhTranspose(Uh, Vh);

    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  } 
  

  //! transpose of ApplyRhBoundaryTransposeGen
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::ApplyRhBoundaryGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyRhBoundary(Uh, Vh);

    Mlt(rh_boundary, Uh, Vh);
  }


  //! computes v = grad u only for boundary dofs
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::ApplyRhBoundaryTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyRhBoundaryTranspose(Uh, Vh);

    Mlt(SeldonTrans, rh_boundary, Uh, Vh);
  } 

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyConstantRh(Uh, Vh);

    Mlt(rh_weight, Uh, Vh);
  }
  
  
  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    then Vx = Rh^1 Uh, Vy = Rh^2 Uh
   */
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyConstantRhSplit(Uh, Vx, Vy);

    Mlt(rh_weightX, Uh, Vx);
    Mlt(rh_weightY, Uh, Vy);
  }
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyConstantRhTransposeGen(Uh, Vh);

    Mlt(SeldonTrans, rh_weight, Uh, Vh);
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
  void QuadrangleLobatto
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyShTranspose(num_loc, Uh, Vh, r);

    if ((r == 0) || (r == order))
      {
	for (int i = 0; i <= order; i++)
	  Vh(i) = Uh(EdgesDof(i, num_loc));
      }
    else
      {
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	Vh.Fill(0);
	for (int i = 0; i <= order; i++)
	  {
	    int num_dof = EdgesDof(i, num_loc);
	    vloc = Uh(num_dof);
	    for (int j = 0; j <= r; j++)
	      Vh(j) += Ah(i, j)*vloc;
	  }
      }
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
  void QuadrangleLobatto
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplySh(alpha, num_loc, Uh, Vh, r);

    int num_dof;
    if ((r == 0) || (r == order))
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    Vh(num_dof) += alpha*Uh(i);
	  }
      }
    else
      {
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    vloc = 0;
	    for (int j = 0; j <= r; j++)
	      vloc += Ah(i, j)*Uh(j);
	    
	    Vh(num_dof) += alpha*vloc;
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
  void QuadrangleLobatto
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyNablaShTranspose(num_loc, Uh, Vh, r);

    Vector2 Uquad(2*(order+1));
    Uquad.Fill(0);
    switch (num_loc)
      {
      case 0 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*i) += Uh(NumFct2D(k, 0)) * lob_basis.GradPhi(k, i);
                Uquad(2*i+1) += Uh(NumFct2D(i, k)) * lob_basis.GradPhi(k, 0);
              }            
          }
        break;
      case 1 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*i) += Uh(NumFct2D(k, i)) * lob_basis.GradPhi(k, order);
                Uquad(2*i+1) += Uh(NumFct2D(order, k)) * lob_basis.GradPhi(k, i);
              }
          }
        break;
      case 2 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*(order-i)) += Uh(NumFct2D(k, order)) * lob_basis.GradPhi(k, i);
                Uquad(2*(order-i)+1) += Uh(NumFct2D(i, k)) * lob_basis.GradPhi(k, order);
              }
          }
        break;
      case 3 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*(order-i)) += Uh(NumFct2D(k, i)) * lob_basis.GradPhi(k, 0);
                Uquad(2*(order-i)+1) += Uh(NumFct2D(0, k)) * lob_basis.GradPhi(k, i);
              }
          }
        break;
      }
    
    if ((r == 0) || (r == order))
      {
        Copy(Uquad, Vh);
      }
    else
      {
        const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	Vh.Fill(0);
	for (int i = 0; i <= order; i++)
	  {
            for (int j = 0; j <= r; j++)
	      {
                Vh(2*j) += Ah(i, j)*Uquad(2*i);
                Vh(2*j+1) += Ah(i, j)*Uquad(2*i+1);
              }
	  }
      }
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
  void QuadrangleLobatto
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (order_quad != order)
      return QuadrangleGauss::ApplyNablaSh(alpha, num_loc, Uh, Vh, r);

    Vector2 Uquad(2*(order+1));
    if ((r == 0) || (r == order))
      {
        Copy(Uh, Uquad);
      }
    else
      {
        const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	Uquad.Fill(0);
	for (int i = 0; i <= order; i++)
	  {
            for (int j = 0; j <= r; j++)
	      {
                Uquad(2*i) += Ah(i, j)*Uh(2*j);
                Uquad(2*i+1) += Ah(i, j)*Uh(2*j+1);
              }
	  }
      }
    
    Mlt(alpha, Uquad);
    typedef typename Vector2::value_type Complexe;
    switch (num_loc)
      {
      case 0 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*i);
            Complexe fy = Uquad(2*i+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, 0)) += fx*lob_basis.GradPhi(k, i);
                Vh(NumFct2D(i, k)) += fy*lob_basis.GradPhi(k, 0);
              }
          }
        break;
      case 1 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*i);
            Complexe fy = Uquad(2*i+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, i)) += fx*lob_basis.GradPhi(k, order);
                Vh(NumFct2D(order, k)) += fy*lob_basis.GradPhi(k, i);
              }
          }
        break;
      case 2 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*(order-i));
            Complexe fy = Uquad(2*(order-i)+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, order)) += fx*lob_basis.GradPhi(k, i);
                Vh(NumFct2D(i, k)) += fy*lob_basis.GradPhi(k, order);
              }
          }
        break;
      case 3 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*(order-i));
            Complexe fy = Uquad(2*(order-i)+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, i)) += fx*lob_basis.GradPhi(k, 0);
                Vh(NumFct2D(0, k)) += fy*lob_basis.GradPhi(k, i);
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
  void QuadrangleLobatto::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    if (order_quad != order)
      return QuadrangleGauss::AddConstantMassMatrix(m, n, mass, A);

    T zero; SetComplexZero(zero);
    if (mass == zero)
      return;
    
    for (int i = 0; i < nb_dof_loc; i++)
      A.AddInteraction(m+i, n+i, mass*this->WeightsND(i));
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
  void QuadrangleLobatto::
  AddConstantElemMatrixGen(int m, int n, const T& mass,
			   const TinyMatrix<T, Prop, 2, 2>& C,
			   const TinyVector<T, 2>& D, 
			   const TinyVector<T, 2>& E,
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    if (order_quad != order)
      return QuadrangleGauss::AddConstantElemMatrix(m, n, mass, C, D, E, null_term, A);

    const VectReal_wp& weights1d = this->Weights1D();
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
        {
          int i = NumFct2D(i1, i2);

	  // mass term
	  if (!null_term(0))
	    A.AddInteraction(m+i, n+i, mass*this->WeightsND(i));
	  
	  // stiffness term
	  if (!null_term(1))
	    {
	      for (int j1 = 0; j1 <= order; j1++)
		A.AddInteraction(m+i, n+NumFct2D(j1, i2), C(0,0)*stiffness_matrix1D(i1, j1)*weights1d(i2));
          
	      for (int j2 = 0; j2 <= order; j2++)
		A.AddInteraction(m+i, n+NumFct2D(i1, j2), C(1,1)*stiffness_matrix1D(i2, j2)*weights1d(i1));
              
	      for (int j1 = 0; j1 <= order; j1++)
		for (int j2 = 0; j2 <= order; j2++)
		  {
		    int j = NumFct2D(j1, j2);
		    A.AddInteraction(m+i, n+j, C(0,1)*lob_basis.GradPhi(i1, j1)*lob_basis.GradPhi(j2, i2)
				     *weights1d(j1)*weights1d(i2)
				     + C(1,0)*lob_basis.GradPhi(j1, i1)*lob_basis.GradPhi(i2, j2)
				     *weights1d(i1)*weights1d(j2));
		  }
	    }
	  
	  // gradient terms
          for (int k = 0; k <= order; k++)
            {
	      if (!null_term(2))
                {
		  A.AddInteraction(i+m, n+NumFct2D(k, i2), D(0)*lob_basis.GradPhi(i1, k)
				   *weights1d(i2)*weights1d(k));
		  
		  A.AddInteraction(i+m, n+NumFct2D(i1, k), D(1)*lob_basis.GradPhi(i2, k)
				   *weights1d(i1)*weights1d(k));
		}
	      
	      if (!null_term(3))
		{
		  A.AddInteraction(m+NumFct2D(k, i2), n+i, E(0)*lob_basis.GradPhi(i1, k)
				   *weights1d(i2)*weights1d(k));
              
		  A.AddInteraction(m+NumFct2D(i1, k), n+i, E(1)*lob_basis.GradPhi(i2, k)
				   *weights1d(i1)*weights1d(k));
		}
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
  void QuadrangleLobatto
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
			     const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			     const Vector<TinyVector<T, 2> >& D,
			     const Vector<TinyVector<T, 2> >& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    if (order_quad != order)
      return QuadrangleGauss::AddVariableElemMatrix(off_row, off_col, mass, C, D, E, null_term, mat);

    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc, zero;
    bool symmetric = mat.IsSymmetric();
    SetComplexZero(zero);
    const Matrix<Real_wp>& dGL_GL = lob_basis.GradPhi();

    // Loop on the local dof
    // we use a tensorial numbering
    for (int k = 0; k <= order; k++)
      for (int l = 0; l <= order; l++)
        {
	  val.Zero();
	  int row = NumFct2D(k, l) + off_row;
          
	  if (!null_term(0))
	    val(row-off_row) = mass(row-off_row);
	  
          // Second loop on the local dof
          // we use a tensorial numbering
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j <= order; j++)	      
              {
		int col = NumFct2D(i, j) + off_col;
                if (!symmetric||(row <= col))
                  {
                    // Now we compute the interaction between row and col
                    // in the stiffness matrix
                    vloc = zero; 
                    
		    if (!null_term(1))
		      {
			if (j == l)
			  // du/dx dv/dx
			  for (int n = 0; n <= order; n++)
			    vloc += C(NumQuad2D(n, j))(0, 0)*dGL_GL(i, n)*dGL_GL(k, n);
			
			if (i == k)
			  // du/dy dv/dy
			  for (int n = 0; n <= order; n++)
			    vloc += C(NumQuad2D(i, n))(1, 1)*dGL_GL(j, n)*dGL_GL(l, n);
                    
			// du/dx dv/dy
			vloc += C(NumQuad2D(k, j))(1, 0)*dGL_GL(i, k)*dGL_GL(l, j);
			
			// du/dy dv/dx
			vloc += C(NumQuad2D(i, l))(0, 1)*dGL_GL(k, i)*dGL_GL(j, l);
		      }
		    
		    if (!null_term(2))
		      {
			if (j == l)
			  // u dv/dx
			  vloc += D(NumQuad2D(i, j))(0)*dGL_GL(k, i);
			
			if (i == k)
			  // u dv/dy
			  vloc += D(NumQuad2D(i, j))(1)*dGL_GL(l, j);
		      }

		    if (!null_term(3))
		      {
			if (j == l)
			  // du/dx v
			  vloc += E(NumQuad2D(k, l))(0)*dGL_GL(i, k);
			
			if (i == k)
			  // du/dy v
			  vloc += E(NumQuad2D(k, l))(1)*dGL_GL(j,l);
		      }
		    
		    val(col-off_col) += vloc;
                  } // end if row <= col
              } // end loop on columns

	  mat.AddInteractionRow(row, this->nb_dof_loc, num, val);
	  
        } // end loop on rows
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
  void QuadrangleLobatto::
  AddVariableMassMatrixGen(int off_row, int off_col,
			   const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    if (order_quad != order)
      return QuadrangleGauss::AddVariableMassMatrix(off_row, off_col, A, mat);

    for (int i = 0; i < nb_dof_loc; i++)
      mat.AddInteraction(i+off_row, i+off_col, A(i));
  }

    
  //! computation of gradient on quadrature points of the edge of the square
  /*!
    \param[in] pts reference points
    \param[in] mat jacobian matrices
    \param[in] feval components of u (by dofs)
    \param[out] res value of gradient of u on quadrature points of the edge
    \param[in] vars given problem
    \param[in] FaceBasis given finite element 
    \param[in] n element number
    \param[in] num_loc local number of the edge
   */
  template<class Vector1, class Vector2>
  void QuadrangleLobatto::
  ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    res.Reallocate(2*(order+1));
    if (num_loc == 0)
      for (int i = 0; i <= order; i++)
	{
	  res(2*i) = 0; res(2*i+1) = 0;
	  for (int j = 0; j <= order; j++)
            {
              res(2*i) += lob_basis.GradPhi(j,i)*feval(NumFct2D(j,0));
              res(2*i+1) += lob_basis.GradPhi(j,0)*feval(NumFct2D(i,j));
            }
	}
    else if (num_loc == 1)
      for (int i = 0; i <= order; i++)
	{
          res(2*i) = 0; res(2*i+1) = 0;
	  for (int j = 0; j <= order; j++)
            {
              res(2*i) += lob_basis.GradPhi(j,order)*feval(NumFct2D(j,i));
              res(2*i+1) += lob_basis.GradPhi(j,i)*feval(NumFct2D(order,j));
            }
	}
    else if (num_loc == 2)
      for (int i = 0; i <= order; i++)
	{
          res(2*(order-i)) = 0; res(2*(order-i)+1) = 0;
	  for (int j = 0; j <= order; j++)
            {
              res(2*(order-i)) += lob_basis.GradPhi(j,i)*feval(NumFct2D(j,order));
              res(2*(order-i)+1) += lob_basis.GradPhi(j,order)*feval(NumFct2D(i,j));
            }
	}
    else if (num_loc == 3)
      for (int i = 0; i <= order; i++)
	{
          res(2*(order-i)) = 0; res(2*(order-i)+1) = 0;
	  for (int j = 0; j <= order; j++)
            {
              res(2*(order-i)) += lob_basis.GradPhi(j,0)*feval(NumFct2D(j,i));
              res(2*(order-i)+1) += lob_basis.GradPhi(j,i)*feval(NumFct2D(0,j));
            }
	}
  }
    
  
  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void QuadrangleLobatto::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    if (order_quad != order)
      return QuadrangleGauss::GetValuePhiOnQuadraturePoint(k, phi);

    phi.Reallocate(nb_dof_loc);
    phi.Fill(0);
    phi(k) = 1.0;
  }
  
  
  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
  */
  void QuadrangleLobatto::GetGradientPhiOnQuadraturePoint(int k, VectR2& grad_phi) const
  {
    if (order_quad != order)
      return QuadrangleGauss::GetGradientPhiOnQuadraturePoint(k, grad_phi);

    grad_phi.Reallocate(nb_dof_loc);
    FillZero(grad_phi);
    int k1 = CoordinateDofs(k, 0);
    int k2 = CoordinateDofs(k, 1);
    for (int i = 0; i <= order; i++)
      {
        grad_phi(NumFct2D(i, k2))(0) = lob_basis.GradPhi(i, k1);
        grad_phi(NumFct2D(k1, i))(1) = lob_basis.GradPhi(i, k2);
      }
  }

  
  //! retrieves values of a single basis function on all quadrature points
  void QuadrangleLobatto::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    if (order_quad != order)
      return QuadrangleGauss::GetValueSinglePhiQuadrature(i, phi);

    phi.Zero();
    phi(i) = Real_wp(1);
  }


  //! retrieves values and gradients of a single basis function on all quadrature points
  void QuadrangleLobatto
  ::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR2& grad_phi) const
  {
    if (order_quad != order)
      return QuadrangleGauss::GetGradientSinglePhiQuadrature(i, phi, grad_phi);

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
  
  
  //! selecting dofs near current dof
  /*!
    \param[in] pos position of the current dof
    \param[in] DofUsed array to know which dofs are already used
    \param[out] ListeDof dofs near the current dof and not used
    \param[in] nb_dof number of dofs we want to pick
   */
  void QuadrangleLobatto::PickNearDofs(int pos, const VectBool& DofUsed,
				       IVect& ListeDof, int nb_dof) const
  {
    if (nb_dof <= 0)
      return;
    
    const Matrix<int>& CoordinateNodes = this->GetCoordinateNodes2D();
    int i0 = CoordinateNodes(pos,0); int j0 = CoordinateNodes(pos,1);
    ListeDof.Reallocate(nb_dof); ListeDof.Fill(-1);
    int node, k = 1, nb = 0;
    if (!DofUsed(pos))
      ListeDof(nb++) = pos;
    
    while (nb < nb_dof)
      {
	// loop on concentric square at distance k
	if ((j0-k) >= 0)
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(m,j0-k);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	if (((i0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(i0+k,m);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	if (((j0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(m,j0+k);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	if (((i0-k) >= 0)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	    if (nb < nb_dof)
	      {
		node = NumFct2D(i0-k,m);
		if (!DofUsed(node))
		  ListeDof(nb++) = node;
	      }
	
	k++;
	if (k > order)
	  nb = nb_dof;
	
      }
	
  }

  
  //! displays informations about class QuadrangleLobatto
  ostream& operator <<(ostream& out, const QuadrangleLobatto& e)
  {
    out<<static_cast<const QuadrangleReference<1>&>(e);
    out<<"Number of degrees of freedom on a face "<<e.nb_dof_loc<<endl;
    return out;
  }
  
}
  
#define MONTJOIE_FILE_QUADRANGLE_LOBATTO_CXX
#endif
