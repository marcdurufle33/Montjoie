#ifndef MONTJOIE_FILE_HEXAHEDRON_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHierarchic::HexahedronHierarchic() : HexahedronReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
  }
  

  size_t HexahedronHierarchic::GetMemorySize() const
  {
    size_t taille = HexahedronReference<1>::GetMemorySize();
    taille += jacobi_11_pol.GetMemorySize();
    taille += CoefLegendre.GetMemorySize();
    taille += ShLoc.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize();
    taille += sh1_loc.GetMemorySize() + sh2_loc.GetMemorySize()
      + rh_loc.GetMemorySize() + const_rh.GetMemorySize();
    taille += MhLoc.GetMemorySize();
    taille += NumDofs3D.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize() + NumDofs2D.GetMemorySize() 
      + CoordinateDofs_quad.GetMemorySize();
    taille += Seldon::GetMemorySize(NumProjOperator);
    return taille;
  }

  //! how to number mesh
  void HexahedronHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg != ElementReference_Base::CONTINUOUS)
      return HexahedronReference<1>::ConstructNumberMap(nmap, dg);

    nmap.SetFormulationDG(dg);
    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
    nmap.SetNbDofQuadrangle(order, (order-1)*(order-1));
        
    nmap.SetNbDofHexahedron(order, (order-1)*(order-1)*(order-1));
    	
    nmap.SetEqualEdgesDofSymmetry(order, order-1);
    nmap.SetOddEdgesDofToSkewSymmetric(order);
    
    // rotation of dofs on faces
    Matrix<int> FacesDof_Rotation(8, (order-1)*(order-1));
    Matrix<bool> SignDof_Rotation(8, (order-1)*(order-1));
    SignDof_Rotation.Fill(false);
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  int node = NumDofs2D(i, j) - 4*order;
	  int node2 = NumDofs2D(j, i) - 4*order;
	  FacesDof_Rotation(0, node) = node;
	  FacesDof_Rotation(1, node) = node2;
	  SignDof_Rotation(1, node) = (i%2 == 0);
	  FacesDof_Rotation(2, node) = node;
	  SignDof_Rotation(2, node) = ((i%2 == 0)^(j%2 == 0));
	  FacesDof_Rotation(3, node) = node2;
	  SignDof_Rotation(3, node) = (j%2 == 0);
	  FacesDof_Rotation(4, node) = node2;
	  FacesDof_Rotation(5, node) = node;
	  SignDof_Rotation(5, node) = (i%2 == 0);
	  FacesDof_Rotation(6, node) = node2;
	  SignDof_Rotation(6, node) = ((i%2 == 0)^(j%2 == 0));
	  FacesDof_Rotation(7, node) = node;
	  SignDof_Rotation(7, node) = (j%2 == 0);
	}

    nmap.SetFacesDofRotationQuad(order, FacesDof_Rotation);
    nmap.SetSignDofRotationQuad(order, SignDof_Rotation);

  }
  
  
  //! constructing finite element
  void HexahedronHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						    int rsurf_tri, int rsurf_quad,
						    int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    HexahedronReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    this->ConstructElementaryMatrix(*this);

    ConstructMassMatrix();
    ConstructStiffnessMatrix();
  }
  
  
  //! construction of basis functions
  void HexahedronHierarchic::ConstructFunctions()
  {
    int r = order;    
    MeshNumbering<Dimension3>::ConstructHexahedralNumbering(r, NumDofs3D, CoordinateDofs);

    function_basis_quad = new QuadrangleHierarchic();   
    function_basis_quad->ConstructFiniteElement(r);
    element_quad_surf = function_basis_quad;
	
    NumDofs2D = function_basis_quad->GetNumDofs2D();
    
    // changing dofs on face
    int offset = 8 + 12*(r-1) - 4*r;
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  NumDofs3D(0, i, j) = offset + NumDofs2D(i, j);
	  NumDofs3D(i, 0, j) = offset + (r-1)*(r-1) + NumDofs2D(i, j);
	  NumDofs3D(i, j, 0) = offset + 2*(r-1)*(r-1) + NumDofs2D(i, j);
	  NumDofs3D(i, j, r) = offset + 3*(r-1)*(r-1) + NumDofs2D(i, j);
	  NumDofs3D(i, r, j) = offset + 4*(r-1)*(r-1) + NumDofs2D(i, j);
	  NumDofs3D(r, i, j) = offset + 5*(r-1)*(r-1) + NumDofs2D(i, j);
	}
    
    // interior
    offset += 4*r + 6*(r-1)*(r-1);
    for (int diag = 1; diag < order; diag++)
      for (int i = 1; i <= diag; i++)
	for (int j = 1; j <= diag; j++)
	  for (int k = 1; k <= diag; k++)
	    if ((i==diag) || (j==diag) || (k == diag))
	      NumDofs3D(i, j, k) = offset++;
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    CoordinateDofs(NumDofs3D(i, j, k), 0) = i;
	    CoordinateDofs(NumDofs3D(i, j, k), 1) = j;
	    CoordinateDofs(NumDofs3D(i, j, k), 2) = k;
	  }
    
    Globatto<Real_wp> lob_basis;
    lob_basis.ConstructQuadrature(r, lob_basis.QUADRATURE_LOBATTO);
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);

    nb_dof_quad = (r+1)*(r+1);
    nb_dof_loc = (r+1)*(r+1)*(r+1);
    nb_dof_boundaries = 6*r*r + 2;
    
    VectReal_wp points_dof1d; VectR2 points_dof2d_quad; VectR3 points_dof3d;
    points_dof1d = this->Points1D();
    points_dof2d_quad = this->Points2D_quad();
    points_dof3d = this->PointsND();

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDof2D_quad(points_dof2d_quad);
    this->SetPointsDofND(points_dof3d);

    this->elt_geom.dof_equal_nodal = false;

    this->nb_points_dof_inside = nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;    
    
    FacesDof.Reallocate((r+1)*(r+1), 6);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D(i, j), 0) = NumDofs3D(0, i, j);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D(i, j), 1) = NumDofs3D(i, 0, j);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D(i, j), 2) = NumDofs3D(i, j, 0);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D(i, j), 3) = NumDofs3D(i, j, r);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D(i, j), 4) = NumDofs3D(i, r, j);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	FacesDof(NumDofs2D(i, j), 5) = NumDofs3D(r, i, j);

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);

    // computing P_m^{1,1}
    GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
    
    CoefLegendre.Reallocate(order-1); CoefLegendre.Fill(0);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*points1d(i) - 1.0, Pn);
        for (int j = 1; j < order; j++)
          CoefLegendre(j-1) += weights1d(i)*square(points1d(i)*(1.0-points1d(i))*Pn(j-1));
      }
    
    for (int j = 1; j < order; j++)
      CoefLegendre(j-1) = 1.0/sqrt(CoefLegendre(j-1));
    
  }
  
  
  //! constructing mass matrix
  void HexahedronHierarchic::ConstructMassMatrix()
  {
    ShLoc.Reallocate(order+1, order_quad+1);
    VectReal_wp Pn;
    const VectReal_wp& points1d = this->Points1D();
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*points1d(i) - 1.0, Pn);
        for (int j = 0; j < order-1; j++)
          Pn(j) *= CoefLegendre(j);
        
        ShLoc(0, i) = 1.0 - points1d(i);
        ShLoc(order, i) = points1d(i);
        for (int j = 1; j < order; j++)
          ShLoc(j, i) = (1.0-points1d(i))*points1d(i)*Pn(j-1);        
      }
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch1.Reallocate(nb_dof_loc, (order+1)*(order_quad+1)*(order+1));
    Ch2.Reallocate((order+1)*(order+1)*(order_quad+1), (order+1)*(order_quad+1)*(order_quad+1));
    Ch3.Reallocate((order+1)*(order_quad+1)*(order_quad+1), nb_points_quadrature_inside);
    
    Matrix<Real_wp, General, ArrayRowSparse> Sh1, Sh2;
    Sh1.Reallocate((order+1)*(order+1), (order_quad+1)*(order+1));
    Sh2.Reallocate((order_quad+1)*(order+1), (order_quad+1)*(order_quad+1));
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i <= order_quad; i++)
      {
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            Sh1.AddInteraction(NumDofs2D(j, k), (order_quad+1)*j + i, ShLoc(k, i));
        
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            Sh2.AddInteraction((order_quad+1)*j + k, NumQuad2D(i, k), ShLoc(j, i));
        
        /* First interpolation along z */
        /* we compute u(x, y, \xi_j), where xi_j are quadrature points for z-coordinate */
        
        // contribution of vertices
        /*
        Ch1.AddInteraction(0, i, 1.0-points1d(i));
        Ch1.AddInteraction(1, (order_quad+1)*(order+1)*order + i, 1.0-points1d(i));
        Ch1.AddInteraction(2, (order_quad+1)*(order+2)*order + i, 1.0-points1d(i));
        Ch1.AddInteraction(3, (order_quad+1)*order + i, 1.0-points1d(i));

        Ch1.AddInteraction(4, i, points1d(i));
        Ch1.AddInteraction(5, (order_quad+1)*(order+1)*order + i, points1d(i));
        Ch1.AddInteraction(6, (order_quad+1)*(order+2)*order + i, points1d(i));
        Ch1.AddInteraction(7, (order_quad+1)*order + i, points1d(i));
        
        // contribution of edges        
        for (int j = 1; j < order; j++)
          {
            // edges on face z = 0
            Ch1.AddInteraction(3 + j, (order_quad+1)*(order+1)*j + i, 1.0-points1d(i));
            Ch1.AddInteraction(3 + (order-1) + j,
            (order_quad+1)*((order+1)*order + j) + i, 1.0-points1d(i));
            Ch1.AddInteraction(3 + 2*(order-1) + j,
            (order_quad+1)*(order + (order+1)*j) + i, 1.0-points1d(i));
            Ch1.AddInteraction(3 + 3*(order-1) + j, (order_quad+1)*j + i, 1.0-points1d(i));
            
            // vertical edges
            Ch1.AddInteraction(3 + 4*(order-1) + j, i, ShLoc(j, i));
            Ch1.AddInteraction(3 + 5*(order-1) + j, (order_quad+1)*(order+1)*order + i,
            ShLoc(j, i));
            Ch1.AddInteraction(3 + 6*(order-1) + j, (order_quad+1)*(order+2)*order + i,
            ShLoc(j, i));
            Ch1.AddInteraction(3 + 7*(order-1) + j, (order_quad+1)*order + i, ShLoc(j, i));
            
            // edges on face z = 1
            Ch1.AddInteraction(3 + 8*(order-1) + j, (order_quad+1)*(order+1)*j + i, points1d(i));
            Ch1.AddInteraction(3 + 9*(order-1) + j, (order_quad+1)*((order+1)*order + j) + i,
            points1d(i));
            Ch1.AddInteraction(3 + 10*(order-1) + j, (order_quad+1)*(order + (order+1)*j) + i,
            points1d(i));
            Ch1.AddInteraction(3 + 11*(order-1) + j, (order_quad+1)*j + i, points1d(i));
          }
        
        // contribution of faces
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            {
              int node = 8 + 6*(order-1) + (j-1)*(order-1) + k-1;
              Ch1.AddInteraction(node, (order_quad+1)*j + i, ShLoc(k, i));
              
              node += (order-1)*(order-1);
              Ch1.AddInteraction(node, (order_quad+1)*(order+1)*j + i, ShLoc(k, i));

              node += (order-1)*(order-1);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + k) + i, 1.0-points1d(i));

              node += (order-1)*(order-1);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + k) + i, points1d(i));

              node += (order-1)*(order-1);
              Ch1.AddInteraction(node, (order_quad+1)*(order + (order+1)*j) + i, ShLoc(k, i));

              node += (order-1)*(order-1);
              Ch1.AddInteraction(node, (order_quad+1)*((order+1)*order + j) + i, ShLoc(k, i));
            }
        
        
        // contribution of the interior
        for (int j = 1; j < order; j++)
          for (int k = 1; k < order; k++)
            for (int m = 1; m < order; m++)
              {
                int node = 4*order*order + 2 + (order-1)*((order-1)*(j-1) + (k-1)) + m-1;
                Ch1.AddInteraction(node, (order_quad+1)*((order+1)*j + k) + i, ShLoc(m, i));
              }
        */
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            for (int m = 0; m <= order; m++)
              Ch1.AddInteraction(NumDofs3D(j, k, m),
                                 (order_quad+1)*((order+1)*j+k)+i, ShLoc(m, i));
        
        /* Second interpolation along y */
        /* we compute u(x, \xi_i, \xi_j), where xi_i are quadrature points in y-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            for (int m = 0; m <= order_quad; m++)
              Ch2.AddInteraction((order_quad+1)*((order+1)*j + k) + m,
                                 (order_quad+1)*((order_quad+1)*j+i)+m, ShLoc(k, i));
        
        /* Third interpolation along x */
        /* we compute u(\xi_i, \xi_j, \xi_k), where xi_i are quadrature points in x-coordinate */
        
        // contribution of 1-D basis functions on quadrature points
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order_quad; k++)
            for (int m = 0; m <= order_quad; m++)
              Ch3.AddInteraction((order_quad+1)*((order_quad+1)*j + k) + m,
                                 NumQuad3D(i, k, m), ShLoc(j, i));
        
      }
    
    Copy(Ch1, ch1_loc);
    Copy(Ch2, ch2_loc);
    Copy(Ch3, ch3_loc);

    Copy(Sh1, sh1_loc);
    Copy(Sh2, sh2_loc);
  }
  
  
  //! constructing stiffness matrix
  void HexahedronHierarchic::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    Matrix<Real_wp, General, ArrayRowSparse> Rh;
    Rh.Reallocate(nb_points_quadrature_inside, 3*nb_points_quadrature_inside);
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
        for (int k = 0; k <= order_quad; k++)
          {
            int node = NumQuad3D(i, j, k);
            for (int m = 0; m <= order_quad; m++)
              {
                Rh.AddInteraction(node, 3*NumQuad3D(m, j, k), lob_quad.GradPhi(i, m));
                Rh.AddInteraction(node, 3*NumQuad3D(i, m, k)+1, lob_quad.GradPhi(j, m));
                Rh.AddInteraction(node, 3*NumQuad3D(i, j, m)+2, lob_quad.GradPhi(k, m));
            }
        }
    
    Copy(Rh, rh_loc);
    ConvertToSparse(const_grad_matrix, const_rh, 1e3*epsilon_machine);
    ConvertToSparse(mass_matrix, MhLoc, 1e3*epsilon_machine);
  }
  

  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector1>
  void HexahedronHierarchic::MltMassMatrixGen(Vector1& Vh) const
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
  void HexahedronHierarchic::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
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
  void HexahedronHierarchic::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Ux(ch3_loc.GetM()), Uy(ch2_loc.GetM());
    Mlt(SeldonTrans, ch1_loc, Uh, Uy);
    Mlt(SeldonTrans, ch2_loc, Uy, Ux);
    Mlt(SeldonTrans, ch3_loc, Ux, Vh);
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
  void HexahedronHierarchic::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHierarchic::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHierarchic::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHierarchic::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }

  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void HexahedronHierarchic::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_rh, Uh, Vh);
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void HexahedronHierarchic::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, const_rh, Uh, Vh);
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
  void HexahedronHierarchic
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
        Vector1 Uquad(this->nb_dof_quad), Ux(sh1_loc.GetM());
        for (int i = 0; i < this->nb_dof_quad; i++)
          Uquad(i) = Uh(FacesDof(i, num_loc));
        
        Mlt(SeldonTrans, sh1_loc, Uquad, Ux);
        Mlt(SeldonTrans, sh2_loc, Ux, Vh);
      }
    else
      {
	int Nquad = this->Points2D_quad().GetM();
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
  void HexahedronHierarchic::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0) || (r == order_quad))
      {
        Vector1 Uquad(this->nb_dof_quad), Ux(sh1_loc.GetM());
        Mlt(sh2_loc, Uh, Ux);
        Mlt(sh1_loc, Ux, Uquad);

        for (int i = 0; i < this->nb_dof_quad; i++)
          Vh(FacesDof(i, num_loc)) += alpha*Uquad(i);
      }
    else
      {
	Vector1 Vlob((order+1)*(order+1));
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
	int Nquad = this->Points2D_quad().GetM();
	for (int i = 0; i < Nquad; i++)
	  {
	    num_dof = FacesDof(i, num_loc);
	    Vh(num_dof) += alpha*Vlob(i);
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
  void HexahedronHierarchic::
  ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
			   const ElementReference_Dim<Dimension3>& FaceCoarse,
			   const ElementReference_Dim<Dimension3>& FaceFine) const
  {
    int rc = FaceCoarse.GetOrder();
    int rf = FaceFine.GetOrder();
    IVect num_coarse(FaceCoarse.GetNbDof()), num_fine(FaceCoarse.GetNbDof());
    num_coarse.Fill();
    num_fine.Fill(-1);
    
    // vertices
    for (int i = 0; i < FaceCoarse.GetNbVertices(); i++)
      num_fine(i) = i;
    
    int nc = FaceCoarse.GetNbVertices();
    int nf = nc;
    
    // edges
    for (int e = 0; e < FaceCoarse.GetNbEdges(); e++)
      for (int i = 0; i < rc-1; i++)
        num_fine(nc + (rc-1)*e + i) = nf + (rf-1)*e + i;
    
    nc += (rc-1)*FaceCoarse.GetNbEdges();
    nf += (rf-1)*FaceCoarse.GetNbEdges();
    
    // faces
    int type = FaceCoarse.GetHybridType();
    for (int f = 0; f < FaceCoarse.GetNbBoundaries(); f++)
      {
        if (MeshNumbering<Dimension3>::IsLocalFaceQuadrilateral(type, f))
          {
            for (int i = 0; i < rc-1; i++)
              for (int j = 0; j < rc-1; j++)
                num_fine(nc + (rc-1)*i + j) = nf + (rf-1)*i + j;
            
            nc += (rc-1)*(rc-1);
            nf += (rf-1)*(rf-1);
          }
        else
          {
            for (int i = 1; i < rf; i++)
              for (int j = 1; j < rf-i; j++)
                {
                  if ((i+j) < rc)
                    num_fine(nc++) = nf;
                  
                  nf++;
                }
          }
      }
    
    // interior
    for (int i = 0; i < rc-1; i++)
      for (int j = 0; j < rc-1; j++)
        for (int k = 0; k < rc-1; k++)
          num_fine(nc++) = nf + (rf-1)*((rf-1)*i + j) + k;
    
    proj.SetIdentity(num_coarse, num_fine, FaceFine.GetNbDof());
  }
  
    
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point
  */
  void HexahedronHierarchic::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    // dofs on vertices
    phi(0) = (1.0-x)*(1.0-y)*(1.0-z);
    phi(1) = x*(1.0-y)*(1.0-z);
    phi(2) = x*y*(1.0-z);
    phi(3) = (1.0-x)*y*(1.0-z);
    phi(4) = (1.0-x)*(1.0-y)*z;
    phi(5) = x*(1.0-y)*z;
    phi(6) = x*y*z;
    phi(7) = (1.0-x)*y*z;
    
    if (order <= 1)
      return;
    
    // dof on edges
    int node = 8;
    VectReal_wp Pn1, Pn2, Pn3;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*x-1, Pn1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*y-1, Pn2);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*z-1, Pn3);
    for (int i = 0; i < order-1; i++)
      {
        Pn1(i) *= CoefLegendre(i);
        Pn2(i) *= CoefLegendre(i);
        Pn3(i) *= CoefLegendre(i);
      }
    
    for (int i = 1; i < order; i++)
      phi(node++) = x*(1.0-x)*(1-y)*(1-z)*Pn1(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = x*y*(1.0-y)*(1-z)*Pn2(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = x*(1.0-x)*y*(1-z)*Pn1(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = (1.0-x)*y*(1.0-y)*(1-z)*Pn2(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = (1.0-x)*(1.0-y)*z*(1-z)*Pn3(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = x*(1.0-y)*z*(1-z)*Pn3(i-1);

    for (int i = 1; i < order; i++)
      phi(node++) = x*y*z*(1-z)*Pn3(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = (1.0-x)*y*z*(1-z)*Pn3(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = x*(1.0-x)*(1-y)*z*Pn1(i-1);

    for (int i = 1; i < order; i++)
      phi(node++) = x*y*(1.0-y)*z*Pn2(i-1);

    for (int i = 1; i < order; i++)
      phi(node++) = x*(1.0-x)*y*z*Pn1(i-1);
    
    for (int i = 1; i < order; i++)
      phi(node++) = (1.0-x)*y*(1.0-y)*z*Pn2(i-1);
    
    //  dofs on faces
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	phi(NumDofs3D(0, i, j)) = (1.0-x)*y*(1.0-y)*Pn2(i-1)*z*(1.0-z)*Pn3(j-1);
	
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	phi(NumDofs3D(i, 0, j)) = x*(1.0-x)*Pn1(i-1)*z*(1.0-z)*Pn3(j-1)*(1.0-y);
    
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	phi(NumDofs3D(i, j, 0)) = x*(1.0-x)*Pn1(i-1)*y*(1.0-y)*Pn2(j-1)*(1.0-z);
    
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	phi(NumDofs3D(i, j, order)) = x*(1.0-x)*Pn1(i-1)*y*(1.0-y)*Pn2(j-1)*z;

    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	phi(NumDofs3D(i, order, j)) = x*(1.0-x)*Pn1(i-1)*z*(1.0-z)*Pn3(j-1)*y;
    
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	phi(NumDofs3D(order, i, j)) = y*(1.0-y)*Pn2(i-1)*z*(1.0-z)*Pn3(j-1)*x;
    
    // dof inside the hexahedron
    Real_wp bubble = x*(1.0-x)*y*(1.0-y)*z*(1.0-z);
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	for (int k = 1; k < order; k++)
	  phi(NumDofs3D(i, j, k)) = bubble*Pn1(i-1)*Pn2(j-1)*Pn3(k-1);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res gradient of basis functions on point
  */
  void HexahedronHierarchic::ComputeGradientPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    int r = this->order;
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    // dofs on vertices
    res(0).Init(-(1.0-y)*(1.0-z), -(1.0-x)*(1.0-z), -(1.0-x)*(1.0-y));
    res(1).Init((1.0-y)*(1.0-z), -x*(1.0-z), -x*(1.0-y));
    res(2).Init(y*(1.0-z), x*(1.0-z), -x*y);
    res(3).Init(-y*(1.0-z), (1.0-x)*(1.0-z), -(1.0-x)*y);
    res(4).Init(-(1.0-y)*z, -(1.0-x)*z, (1.0-x)*(1.0-y));
    res(5).Init((1.0-y)*z, -x*z, x*(1.0-y));
    res(6).Init(y*z, x*z, x*y);
    res(7).Init(-y*z, (1.0-x)*z, (1.0-x)*y);
    
    if (order <= 1)
      return;
    
    // dofs on edges
    VectReal_wp Pn1, dPn1, Pn2, dPn2, Pn3, dPn3;
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*x-1, Pn1, dPn1);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*y-1, Pn2, dPn2);
    EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*z-1, Pn3, dPn3);
    for (int i = 0; i < order-1; i++)
      {
        Pn1(i) *= CoefLegendre(i); dPn1(i) *= CoefLegendre(i);
        Pn2(i) *= CoefLegendre(i); dPn2(i) *= CoefLegendre(i); 
        Pn3(i) *= CoefLegendre(i); dPn3(i) *= CoefLegendre(i);
      }
    
    int node = 8;
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = (1-y)*(1-z)*( (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1));
	res(node)(1) = -x*(1.0-x)*(1-z)*Pn1(i-1);
	res(node)(2) = -x*(1.0-x)*(1-y)*Pn1(i-1);
	node++;
      }

    for (int i = 1; i < r; i++)
      {
	res(node)(0) = y*(1.0-y)*(1-z)*Pn2(i-1);
	res(node)(1) = x*(1.0-z)*( (1.0-2*y)*Pn2(i-1) + 2*y*(1.0-y)*dPn2(i-1));
	res(node)(2) = -x*y*(1.0-y)*Pn2(i-1);
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = y*(1-z)*( (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1) );
	res(node)(1) = x*(1.0-x)*(1.0-z)*Pn1(i-1);
	res(node)(2) = -x*(1.0-x)*y*Pn1(i-1);
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = -y*(1.0-y)*(1-z)*Pn2(i-1);
	res(node)(1) = (1.0-x)*(1-z)*( (1.0-2*y)*Pn2(i-1) + 2*y*(1.0-y)*dPn2(i-1));
	res(node)(2) = -(1.0-x)*y*(1.0-y)*Pn2(i-1);
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = -(1.0-y)*z*(1-z)*Pn3(i-1);
	res(node)(1) = -(1.0-x)*z*(1-z)*Pn3(i-1);
	res(node)(2) = (1.0-x)*(1.0-y)*( (1.0-2*z)*Pn3(i-1) + 2*z*(1.0-z)*dPn3(i-1));
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = (1.0-y)*z*(1-z)*Pn3(i-1);
	res(node)(1) = -x*z*(1-z)*Pn3(i-1);
	res(node)(2) = x*(1.0-y)*( (1.0-2*z)*Pn3(i-1) + 2*z*(1.0-z)*dPn3(i-1));
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = y*z*(1-z)*Pn3(i-1);
	res(node)(1) = x*z*(1-z)*Pn3(i-1);
	res(node)(2) = x*y*( (1.0-2*z)*Pn3(i-1) + 2*z*(1-z)*dPn3(i-1));
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = -y*z*(1-z)*Pn3(i-1);
	res(node)(1) = (1.0-x)*z*(1-z)*Pn3(i-1);
	res(node)(2) = (1.0-x)*y*( (1.0-2*z)*Pn3(i-1) + 2*z*(1-z)*dPn3(i-1));
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = (1-y)*z*( (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1));
	res(node)(1) = -x*(1.0-x)*z*Pn1(i-1);
	res(node)(2) = x*(1.0-x)*(1-y)*Pn1(i-1);
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = y*(1.0-y)*z*Pn2(i-1);
	res(node)(1) = x*z*( (1.0-2*y)*Pn2(i-1) + 2*y*(1.0-y)*dPn2(i-1));
	res(node)(2) = x*y*(1.0-y)*Pn2(i-1);
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = y*z*( (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1));
	res(node)(1) = x*(1.0-x)*z*Pn1(i-1);
	res(node)(2) = x*(1.0-x)*y*Pn1(i-1);
	node++;
      }
    
    for (int i = 1; i < r; i++)
      {
	res(node)(0) = -y*(1.0-y)*z*Pn2(i-1);
	res(node)(1) = (1.0-x)*z*( (1.0-2*y)*Pn2(i-1) + 2*y*(1.0-y)*dPn2(i-1));
	res(node)(2) = (1.0-x)*y*(1.0-y)*Pn2(i-1);
	node++;
      }
    
    //  dofs on faces
    Real_wp v1, dv1, v2, dv2, v3, dv3;
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  v1 = y*(1.0-y)*Pn2(i-1);
	  dv1 = (1.0-2*y)*Pn2(i-1) + 2*y*(1.0-y)*dPn2(i-1);

	  v2 = z*(1.0-z)*Pn3(j-1);
	  dv2 = (1.0-2*z)*Pn3(j-1) + 2*z*(1.0-z)*dPn3(j-1);
	
	  node = NumDofs3D(0, i, j);
	  res(node).Init(-v1*v2, (1.0-x)*dv1*v2, (1.0-x)*v1*dv2);
	}

    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  v1 = x*(1.0-x)*Pn1(i-1);
	  dv1 = (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1);
	  
	  v2 = z*(1.0-z)*Pn3(j-1);
	  dv2 = (1.0-2*z)*Pn3(j-1) + 2*z*(1.0-z)*dPn3(j-1);
	  
	  node = NumDofs3D(i, 0, j);
	  res(node).Init((1.0-y)*dv1*v2, -v1*v2, (1.0-y)*v1*dv2);
	}
    
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  v1 = x*(1.0-x)*Pn1(i-1);
	  dv1 = (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1);
	  
	  v2 = y*(1.0-y)*Pn2(j-1);
	  dv2 = (1.0-2*y)*Pn2(j-1) + 2*y*(1.0-y)*dPn2(j-1);
	  
	  node = NumDofs3D(i, j, 0);
	  res(node).Init((1.0-z)*dv1*v2, (1.0-z)*v1*dv2, -v1*v2);
	}

    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  v1 = x*(1.0-x)*Pn1(i-1);
	  dv1 = (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1);
	  
	  v2 = y*(1.0-y)*Pn2(j-1);
	  dv2 = (1.0-2*y)*Pn2(j-1) + 2*y*(1.0-y)*dPn2(j-1);
	  
	  node = NumDofs3D(i, j, order);
	  res(node).Init(z*dv1*v2, z*v1*dv2, v1*v2);
	}

    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  v1 = x*(1.0-x)*Pn1(i-1);
	  dv1 = (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1);

	  v2 = z*(1.0-z)*Pn3(j-1);
	  dv2 = (1.0-2*z)*Pn3(j-1) + 2*z*(1.0-z)*dPn3(j-1);
	  
	  node = NumDofs3D(i, order, j);
	  res(node).Init(y*dv1*v2, v1*v2, y*v1*dv2);
	}

    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	{
	  v1 = y*(1.0-y)*Pn2(i-1);
	  dv1 = (1.0-2*y)*Pn2(i-1) + 2*y*(1.0-y)*dPn2(i-1);
	  
	  v2 = z*(1.0-z)*Pn3(j-1);
	  dv2 = (1.0-2*z)*Pn3(j-1) + 2*z*(1.0-z)*dPn3(j-1);
	  
	  node = NumDofs3D(order, i, j);
	  res(node).Init(v1*v2, x*dv1*v2, x*v1*dv2);	
	}
    
    // dofs inside the hexahedron
    for (int i = 1; i < r; i++)
      for (int j = 1; j < r; j++)
	for (int k = 1; k < r; k++)
	  {
	    v1 = x*(1.0-x)*Pn1(i-1);
	    dv1 = (1.0-2*x)*Pn1(i-1) + 2*x*(1.0-x)*dPn1(i-1);

	    v2 = y*(1.0-y)*Pn2(j-1);
	    dv2 = (1.0-2*y)*Pn2(j-1) + 2*y*(1.0-y)*dPn2(j-1);

	    v3 = z*(1.0-z)*Pn3(k-1);
	    dv3 = (1.0-2*z)*Pn3(k-1) + 2*z*(1.0-z)*dPn3(k-1);
	    
	    node = NumDofs3D(i, j, k);
	    res(node).Init(dv1*v2*v3, v1*dv2*v3, v1*v2*dv3);
	  }
  }
  
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void HexahedronHierarchic
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    VectComplex_wp feval_weight(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      feval_weight(i) = feval(i)*this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  
  
  //! projection from values on "dof points" to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void HexahedronHierarchic
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    VectReal_wp feval_weight(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      feval_weight(i) = feval(i)*this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }


  void HexahedronHierarchic::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    VectReal_wp contrib_orig(contrib);
    for (int i = 1; i < this->order; i++)
      {
	if (i%2 == this->order%2)
	  {
	    contrib(2*this->order+i+1) = -contrib_orig(3*this->order-i+1);
	    contrib(3*this->order+i) = -contrib_orig(4*this->order-i);
	  }
	else
	  {
	    contrib(2*this->order+i+1) = contrib_orig(3*this->order-i+1);
	    contrib(3*this->order+i) = contrib_orig(4*this->order-i);
	  }	
      }
  }


  void HexahedronHierarchic::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    VectComplex_wp contrib_orig(contrib);
    for (int i = 1; i < this->order; i++)
      {
	if (i%2 == this->order%2)
	  {
	    contrib(2*this->order+i+1) = -contrib_orig(3*this->order-i+1);
	    contrib(3*this->order+i) = -contrib_orig(4*this->order-i);
	  }
	else
	  {
	    contrib(2*this->order+i+1) = contrib_orig(3*this->order-i+1);
	    contrib(3*this->order+i) = contrib_orig(4*this->order-i);
	  }	
      }
  }

  
  //! computation of projection between finite element of different orders  
  void HexahedronHierarchic::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi; VectR2 Points; IVect perm;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    NumProjOperator.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate(order+1, r+1);
	NumProjOperator(r).Reallocate(r+1, r+1);
	Points = Pts(i); perm.Reallocate(Points.GetM()); perm.Fill();
	Sort(Points, perm);
	int nb = 0;
	for (int j1 = 0; j1 <= r; j1++)
	  for (int j2 = 0; j2 <= r; j2++)
	    NumProjOperator(r)(j1, j2) = perm(nb++);
	
        VectReal_wp Pn;
	for (int i2 = 0; i2 <= r; i2++)
	  {
            ProjOperatorQuadOrder(r)(0, i2) = 1.0 - Points(i2)(1);
            ProjOperatorQuadOrder(r)(order, i2) = Points(i2)(1);
            
            EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2.0*Points(i2)(1) - 1.0, Pn);
            for (int i1 = 1; i1 < order; i1++)
              ProjOperatorQuadOrder(r)(i1, i2) = Pn(i1-1)*CoefLegendre(i1-1);
          }
      }
  }


  //! displays details about class HexahedronHierarchic
  ostream& operator <<(ostream& out, const HexahedronHierarchic& e)
  {
    out<<static_cast<const HexahedronReference<1>&>(e);
    return out;
  }
  
}
  
#define MONTJOIE_FILE_HEXAHEDRON_HIERARCHIC_CXX
#endif
