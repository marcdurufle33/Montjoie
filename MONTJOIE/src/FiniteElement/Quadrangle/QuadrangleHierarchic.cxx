#ifndef MONTJOIE_FILE_QUADRANGLE_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  QuadrangleHierarchic::QuadrangleHierarchic() : QuadrangleReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
  }
  

  //! method useful to number mesh
  void QuadrangleHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<1>::ConstructNumberMap(nmap, dg);
    
    if (dg == ElementReference_Base::HDG)
      {
	nmap.SetFormulationDG(ElementReference_Base::HDG);
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, this->order+1);
	nmap.SetNbDofTriangle(this->order, 0);

	nmap.SetEqualEdgesDofSymmetry(this->order, this->order+1);
	nmap.SetEvenEdgesDofToSkewSymmetric(this->order);

	nmap.SetEdgesDofSymmetry(this->order, 0, this->order, false);
	nmap.SetEdgesDofSymmetry(this->order, this->order, 0, false);
	
	return;
      }

    nmap.SetFormulationDG(ElementReference_Base::CONTINUOUS);
    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
    nmap.SetNbDofQuadrangle(order, (order-1)*(order-1));
    nmap.SetNbDofTriangle(order, (order-1)*(order-2)/2);
    
    nmap.SetEqualEdgesDofSymmetry(order, order-1);
    nmap.SetOddEdgesDofToSkewSymmetric(order);
  }


  size_t QuadrangleHierarchic::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<1>::GetMemorySize();
    taille += NumDofsQuad.GetMemorySize();
    taille += jacobi_11_pol.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += ch1_loc.GetMemorySize();
    taille += ch2_loc.GetMemorySize();
    taille += const_rh.GetMemorySize();
    taille += MhLoc.GetMemorySize();
    taille += ShLoc.GetMemorySize();
    taille += InvWeightFct.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    return taille;
  }
  
  
  //! constructing finite element
  void QuadrangleHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						    int rsurf, int type_surf)
  {
    QuadrangleReference<1>::
      ConstructFiniteElement(r, rgeom, rquad, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();    
    
    ConstructElementaryMatrix(*this);

    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    delete this->element_surface;
    EdgeHierarchicReference* edge = new EdgeHierarchicReference();
    edge->ConstructFiniteElement(r, r, rquad, type_surf);
    this->element_surface = edge;    
  }
  
  
  //! constructions of basis functions
  void QuadrangleHierarchic::ConstructFunctions()
  {
    nb_dof_loc = (order+1)*(order+1);
    nb_dof_boundaries = 4*order;

    EdgesDof.Reallocate(order+1, 4);
    // changing EdgesDof
    EdgesDof(0, 0) = 0; EdgesDof(order, 0) = 1;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 0) = 3+i;
    
    EdgesDof(0, 1) = 1; EdgesDof(order, 1) = 2;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 1) = 2+i+order;

    EdgesDof(0, 2) = 2; EdgesDof(order, 2) = 3;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 2) = 1+i+2*order;
    
    EdgesDof(0, 3) = 3; EdgesDof(order, 3) = 0;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 3) = i+3*order;

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    VectR2 points_dof2d; VectReal_wp points_dof1d;
    points_dof2d = this->PointsND();
    this->nb_points_dof_inside = nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    this->elt_geom.dof_equal_quadrature = true;
    
    points_dof1d = this->Points1D();

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    // computing P_m^{1,1}
    GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
    CoefLegendre.Reallocate(order-1); CoefLegendre.Fill(0);
    VectReal_wp Pn;
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*this->Points1D(i) - 1.0, Pn);
        for (int j = 1; j < order; j++)
          CoefLegendre(j-1) += this->Weights1D(i)*square(this->Points1D(i)*(1.0-this->Points1D(i))*Pn(j-1));
      }
    
    for (int j = 1; j < order; j++)
      CoefLegendre(j-1) = 1.0/sqrt(CoefLegendre(j-1));
    
    InvWeightFct.Reallocate(nb_dof_loc); InvWeightFct.Fill(1);
    
    // hierarchical numbering for internal nodes
    NumDofsQuad.Reallocate(order+1, order+1);
    NumDofsQuad(0, 0) = 0;
    NumDofsQuad(order, 0) = 1;
    NumDofsQuad(order, order) = 2;
    NumDofsQuad(0, order) = 3;
    for (int i = 1; i < order; i++)
      {
        NumDofsQuad(i, 0) = 3+i;
        NumDofsQuad(order, i) = 2+order+i;
        NumDofsQuad(order-i, order) = 1+2*order+i;
        NumDofsQuad(0, order-i) = 3*order + i;
      }
    
    int node = 4*order;
    for (int diag = 1; diag < order; diag++)
      {
        for (int i = 1; i < diag; i++)
          NumDofsQuad(diag, i) = node++;
        
        for (int i = 1; i <= diag; i++)
          NumDofsQuad(i, diag) = node++;
      }
  }
  
  
  //! constructing mass matrix
  void QuadrangleHierarchic::ConstructMassMatrix()
  {    
    ShLoc.Reallocate(order+1, order_quad+1);
    VectReal_wp Pn;
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*this->Points1D(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLegendre(j);
        
        ShLoc(0, i) = 1.0 - this->Points1D(i);
        ShLoc(order, i) = this->Points1D(i);
        for (int j = 1; j < order; j++)
          ShLoc(j, i) = (1.0-this->Points1D(i))*this->Points1D(i)*Pn(j-1);        
      }
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2;
    Ch1.Reallocate(nb_dof_loc, (order+1)*(order_quad+1));
    Ch2.Reallocate((order+1)*(order_quad+1), nb_points_quadrature_inside);
    for (int i = 0; i <= order_quad; i++)
      {
        /* First interpolation along y */
        /* we compute u(x, \xi_j), where xi_j are quadrature points in y-coordinate */
        
        // contribution of vertices
        Ch1.AddInteraction(0, i, 1.0-this->Points1D(i));
        Ch1.AddInteraction(1, (order_quad+1)*order+i, 1.0-this->Points1D(i));
        Ch1.AddInteraction(3, i, this->Points1D(i));
        Ch1.AddInteraction(2, (order_quad+1)*order+i, this->Points1D(i));
        
        // contribution of edges        
        for (int j = 1; j < order; j++)
          {
            Ch1.AddInteraction(3+j, (order_quad+1)*j + i, 1.0-this->Points1D(i));
            Ch1.AddInteraction(2+order+j, (order_quad+1)*order + i, ShLoc(j, i));
            if (j%2 == 0)
              {
                Ch1.AddInteraction(1+2*order+j, (order_quad+1)*j + i, -this->Points1D(i));
                Ch1.AddInteraction(3*order+j, i, -ShLoc(j, i));
              }
            else
              {
                Ch1.AddInteraction(1+2*order+j, (order_quad+1)*j + i, this->Points1D(i));
                Ch1.AddInteraction(3*order+j, i, ShLoc(j, i));
              }
            
          }
        
        // contribution of the interior
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            {
              int node = NumDofsQuad(j, k);
              Ch1.AddInteraction(node, (order_quad+1)*j + i, ShLoc(k, i));
            }

        /* Second interpolation along x */
        /* we compute u(\xi_i, \xi_j), where xi_i are quadrature points in x-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            Ch2.AddInteraction((order_quad+1)*j + k, NumQuad2D(i, k), ShLoc(j, i));
        
      }
    
    Copy(Ch1, ch1_loc);
    Copy(Ch2, ch2_loc);
    
  }    
  
  
  //! computation of stiffness matrix R_h
  /*!
    this matrix is used to compute gradients
   */
  void QuadrangleHierarchic::ConstructStiffnessMatrix()
  {    
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    Matrix<Real_wp, General, ArrayRowSparse> Rh;
    Rh.Reallocate(nb_points_quadrature_inside, 2*nb_points_quadrature_inside);
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
        {
          int node = NumQuad2D(i, j);
          for (int k = 0; k <= order_quad; k++)
            {
              Rh.AddInteraction(node, 2*NumQuad2D(k, j), lob_quad.GradPhi(i, k));
              Rh.AddInteraction(node, 2*NumQuad2D(i, k)+1, lob_quad.GradPhi(j, k));
            }
        }
    
    Copy(Rh, rh_loc);
    ConvertToSparse(const_grad_matrix, const_rh, 1e3*epsilon_machine);
    ConvertToSparse(mass_matrix, MhLoc, 1e3*epsilon_machine);
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] res dof components
   */
  template<class Vector1>
  void QuadrangleHierarchic::ComputeProjectionDofGen(const Vector1& feval, Vector1& contrib) const
  {
    Vector1 feval_weight = feval;
    for (int i = 0; i < feval_weight.GetM(); i++)
      feval_weight(i) *= this->WeightsND(i);
    
    ApplyChGen(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }

    
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector1>
  void QuadrangleHierarchic::MltMassMatrixGen(Vector1& Vh) const
  {
    Vector1 Uh(Vh);
    Mlt(MhLoc, Uh, Vh);
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
  void QuadrangleHierarchic::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch2_loc.GetM());
    Mlt(ch2_loc, Uh, Ux);
    Mlt(ch1_loc, Ux, Vh);
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
  void QuadrangleHierarchic::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch2_loc.GetM());
    Mlt(SeldonTrans, ch1_loc, Uh, Ux);
    Mlt(SeldonTrans, ch2_loc, Ux, Vh);
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
  void QuadrangleHierarchic::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
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
  void QuadrangleHierarchic::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
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
  void QuadrangleHierarchic::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(rh_loc, Uh, Vh);
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
  void QuadrangleHierarchic::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void QuadrangleHierarchic::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_rh, Uh, Vh);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void QuadrangleHierarchic::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, const_rh, Uh, Vh);
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
  void QuadrangleHierarchic
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    int num_dof;
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= order_quad; j++)
	      Vh(j) += ShLoc(i, j)*Uh(num_dof);
	  }
      }
    else
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= r; j++)
	      Vh(j) += ProjOperatorOrder(r)(i, j)*Uh(num_dof);
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
  void QuadrangleHierarchic::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    int num_dof;
    typename Vector2::value_type vloc;
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    vloc = 0;
	    for (int j = 0; j <= order_quad; j++)
	      vloc += ShLoc(i, j)*Uh(j);

	    Vh(num_dof) += alpha*vloc;
	  }
      }
    else
      {
	for (int i = 0; i <= order; i++)
	  {
	    num_dof = EdgesDof(i, num_loc);
	    vloc = 0;
	    for (int j = 0; j <= r; j++)
	      vloc += ProjOperatorOrder(r)(i, j)*Uh(j);
	    
	    Vh(num_dof) += alpha*vloc;
	  }
      }
  }
  
  
  //! computation of prolongation operator inside an element
  /*!
    \param[in,out] proj prolongation operator
    \param[in,out] LocalProlongation prolongation operator
    \param[in] FaceCoarse coarse finite element
    \param[in] FaceFine fine finite element
   */
  void QuadrangleHierarchic
  ::ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
			     const ElementReference_Dim<Dimension2>& FaceCoarse,
			     const ElementReference_Dim<Dimension2>& FaceFine) const    
  {
    int rc = FaceCoarse.GetOrder();
    int rf = FaceFine.GetOrder();
    IVect num_coarse(FaceCoarse.GetNbDof()), num_fine(FaceCoarse.GetNbDof());
    num_coarse.Fill();
    num_fine.Fill(-1);
    // vertices
    num_fine(0) = 0; num_fine(1) = 1;
    num_fine(2) = 2; num_fine(3) = 3;
    
    // edges
    int node = 4;
    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 4+i;

    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 4+(rf-1)+i;

    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 4+2*(rf-1)+i;

    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 4+3*(rf-1)+i;
    
    // interior
    const QuadrangleHierarchic& Fb_c = static_cast<const QuadrangleHierarchic&>(FaceCoarse);
    const QuadrangleHierarchic& Fb_f = static_cast<const QuadrangleHierarchic&>(FaceFine);
    for (int i = 0; i < rc-1; i++)
      for (int j = 0; j < rc-1; j++)
        num_fine(Fb_c.NumDofsQuad(i, j)) = Fb_f.NumDofsQuad(i, j);
    
    proj.SetIdentity(num_coarse, num_fine, FaceFine.GetNbDof());
  }
  

  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point_loc
  */
  void QuadrangleHierarchic::ComputeValuesPhiRef(const R2& point_loc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    
    phi(0) = (1.0-point_loc(0))*(1.0-point_loc(1));
    phi(1) = point_loc(0)*(1.0-point_loc(1));
    phi(2) = point_loc(0)*point_loc(1);
    phi(3) = (1.0-point_loc(0))*point_loc(1);
    
    // edge y = 0
    VectReal_wp Pn, Pn2; int num_dof = 4;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(0)-1, Pn);
    for (int j = 0; j < order-1; j++)
      Pn(j) *= CoefLegendre(j);
    
    for (int i = 1; i < order; i++)
      phi(num_dof++) = point_loc(0)*(1.0-point_loc(0))*(1.0-point_loc(1))*Pn(i-1);
    
    // edge x = 1
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(1)-1, Pn);
    for (int j = 0; j < order-1; j++)
      Pn(j) *= CoefLegendre(j);
        
    for (int i = 1; i < order; i++)
      phi(num_dof++) = point_loc(1)*(1.0-point_loc(1))*point_loc(0)*Pn(i-1);
    
    // edge y = 1
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 1.0-2*point_loc(0), Pn);
    for (int j = 0; j < order-1; j++)
      Pn(j) *= CoefLegendre(j);
    
    for (int i = 1; i < order; i++)
      phi(num_dof++) = point_loc(0)*(1.0-point_loc(0))*point_loc(1)*Pn(i-1);
    
    // edge x = 0
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 1.0-2*point_loc(1), Pn);
    for (int j = 0; j < order-1; j++)
      Pn(j) *= CoefLegendre(j);
    
    for (int i = 1; i < order; i++)
      phi(num_dof++) = point_loc(1)*(1.0-point_loc(1))*(1.0-point_loc(0))*Pn(i-1);
    
    // interior of quadrangle
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(0)-1, Pn);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(1)-1, Pn2);
    for (int j = 0; j < order-1; j++)
      {
        Pn(j) *= CoefLegendre(j);
        Pn2(j) *= CoefLegendre(j);
      }
    
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  phi(NumDofsQuad(i, j)) = point_loc(0)*(1.0-point_loc(0))*Pn(i-1)
            *point_loc(1)*(1.0-point_loc(1))*Pn2(j-1);
	}
    
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) *= InvWeightFct(i);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res gradient of basis functions on point_loc
  */
  void QuadrangleHierarchic::ComputeGradientPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc);
    res(0).Init(point_loc(1)-1, point_loc(0)-1);
    res(1).Init(1-point_loc(1), -point_loc(0));
    res(2).Init(point_loc(1), point_loc(0));
    res(3).Init(-point_loc(1), 1-point_loc(0));
    
    // edge y = 0
    VectReal_wp Pn, dPn, Pn2, dPn2; int num_dof = 4;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(0)-1, Pn, dPn);
    for (int j = 0; j < order-1; j++)
      {
        Pn(j) *= CoefLegendre(j);
        dPn(j) *= CoefLegendre(j);
      }
    
    for (int i = 1; i < order; i++)
      {
	res(num_dof)(0) = (1.0-2.0*point_loc(0))*(1.0-point_loc(1))*Pn(i-1) 
	  + 2*point_loc(0)*(1-point_loc(0))*(1.0-point_loc(1))*dPn(i-1);
	
	res(num_dof)(1) = -point_loc(0)*(1-point_loc(0))*Pn(i-1);
	num_dof++;
      }
    
    // edge x = 1
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(1)-1, Pn, dPn);
    for (int j = 0; j < order-1; j++)
      {
        Pn(j) *= CoefLegendre(j);
        dPn(j) *= CoefLegendre(j);
      }
    for (int i = 1; i < order; i++)
      {
	res(num_dof)(0) = point_loc(1)*(1.0-point_loc(1))*Pn(i-1);
	res(num_dof)(1) = (1.0-2*point_loc(1))*point_loc(0)*Pn(i-1)
	  + 2*point_loc(1)*(1.0-point_loc(1))*point_loc(0)*dPn(i-1);
	num_dof++;
      }
    
    // edge y = 1
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 1.0-2*point_loc(0), Pn, dPn);
    for (int j = 0; j < order-1; j++)
      {
        Pn(j) *= CoefLegendre(j);
        dPn(j) *= CoefLegendre(j);
      }
    
    for (int i = 1; i < order; i++)
      {
	res(num_dof)(0) = (1.0-2.0*point_loc(0))*point_loc(1)*Pn(i-1)
	  - 2*point_loc(0)*(1.0-point_loc(0))*point_loc(1)*dPn(i-1);
	res(num_dof)(1) = point_loc(0)*(1.0-point_loc(0))*Pn(i-1);
	num_dof++;
      }
    
    // edge x = 0
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 1.0-2*point_loc(1), Pn, dPn);
    for (int j = 0; j < order-1; j++)
      {
        Pn(j) *= CoefLegendre(j);
        dPn(j) *= CoefLegendre(j);
      }
    
    for (int i = 1; i < order; i++)
      {
	res(num_dof)(0) = -point_loc(1)*(1.0-point_loc(1))*Pn(i-1);
	res(num_dof)(1) = (1.0-2.0*point_loc(1))*(1.0-point_loc(0))*Pn(i-1)
	  - 2.0*point_loc(1)*(1.0-point_loc(1))*(1.0-point_loc(0))*dPn(i-1);
	num_dof++;
      }

    Real_wp vx(1), vy(1), dvx, dvy;
    
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(0)-1, Pn, dPn);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*point_loc(1)-1, Pn2, dPn2);
    for (int j = 0; j < order-1; j++)
      {
        Pn(j) *= CoefLegendre(j);
        dPn(j) *= CoefLegendre(j);
        Pn2(j) *= CoefLegendre(j);
        dPn2(j) *= CoefLegendre(j);
      }
    
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  vx = point_loc(0)*(1.0-point_loc(0))*Pn(i-1);
	  dvx = (1.0-2*point_loc(0))*Pn(i-1) + 2*point_loc(0)*(1.0-point_loc(0))*dPn(i-1); 

	  vy = point_loc(1)*(1.0-point_loc(1))*Pn2(j-1);
	  dvy = (1.0-2.0*point_loc(1))*Pn2(j-1) + 2*point_loc(1)*(1.0-point_loc(1))*dPn2(j-1);
	  
          num_dof = NumDofsQuad(i, j);
	  res(num_dof)(0) = dvx*vy; res(num_dof)(1) = vx*dvy;
	}

    for (int i = 0; i < nb_dof_loc; i++)
      res(i) *= InvWeightFct(i);
  }
  
  
  //! displays informations about class QuadrangleHierarchic
  ostream& operator <<(ostream& out, const QuadrangleHierarchic& e)
  {
    out<<static_cast<const QuadrangleReference<1>&>(e);
    out<<"Number of degrees of freedom on a face "<<e.nb_dof_loc<<endl;
    return out;
  }
  
}
  
#define MONTJOIE_FILE_QUADRANGLE_HIERARCHIC_CXX
#endif
