#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHcurlOptimalFirstFamily::HexahedronHcurlOptimalFirstFamily() : HexahedronReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  

  //! returns size of memory used by the object
  size_t HexahedronHcurlOptimalFirstFamily::GetMemorySize() const
  {
    size_t taille = HexahedronReference<2>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize() + NumDofs_Y.GetMemorySize()
      +NumDofs_Z.GetMemorySize();
    
    taille += GLint_GL.GetMemorySize() + GL_GX.GetMemorySize() + Geom_GL.GetMemorySize()
      + dGL_GL.GetMemorySize()+stiff1d.GetMemorySize();
    
    taille += ValGaussExt.GetMemorySize() + lob_int.GetMemorySize();
    taille += NumQuad2D_gauss.GetMemorySize();
    taille += rh_loc.GetMemorySize()+ch_loc.GetMemorySize();
    taille += ch3_node.GetMemorySize()+ch1_node.GetMemorySize()+ch2_node.GetMemorySize();
    taille += WeightsMassGL.GetMemorySize() + invWeightsMassGL.GetMemorySize()
      + invSqrtWeightsMassGL.GetMemorySize();
    taille += const_grad_matrix.GetMemorySize();
    return taille;
  }
  
    
  //! in order to number a mesh
  void HexahedronHcurlOptimalFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return HexahedronReference<2>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    // 2 r (r-1) dofs on each face
    int nb_dof_inside_quad = 2*order*order;
    nmap.SetNbDofQuadrangle(order, nb_dof_inside_quad);
    
    // 3 r (r-1) (r-1) dofs inside the hexahedron
    nmap.SetNbDofHexahedron(order, nb_dof_loc - nb_dof_boundaries);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
        
    //FindHcurlLinearCombinationRotation(*this, nmap, Fb_quad, Fb_quad);
    
    element_quad_surf->FindHcurlRotationQuad(nmap);
  }
  
  
  //! construction of finite element
  void HexahedronHcurlOptimalFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (rquad <= r)
      rquad = r+1;

    HexahedronReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad); 
    
    lob_basis.ConstructQuadrature(order+1, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    ConstructFunctions();
    
    ConstructStiffnessMatrix();
    ConstructMassMatrix();
    
    FindDofsOnFace(true);
    
    this->Fb_geom.ComputeCoefficientTransformation();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
  }

  
  //! construction of basis functions
  void HexahedronHcurlOptimalFirstFamily::ConstructFunctions()
  {
    nb_dof_tri = 0;
    nb_dof_quad = 2*order*(order+2);
    nb_dof_loc = 3*order*(order+1)*(order+1);
    
    int nb_dof_inside_quad = order*order*2;
    nb_dof_boundaries = 6*nb_dof_inside_quad + 12*order;
    
    QuadrangleHcurlOptimalFirstFamily* Fb_quad_ptr = new QuadrangleHcurlOptimalFirstFamily();
    QuadrangleHcurlOptimalFirstFamily& Fb_quad = *Fb_quad_ptr;
    element_quad_surf = Fb_quad_ptr;

    Fb_quad.ConstructFiniteElement(order);
    
    FacesDof.Reallocate(nb_dof_quad, 6);
    NumDofs_X.Reallocate(order, order+2, order+2);
    NumDofs_Y.Reallocate(order+2, order, order+2);
    NumDofs_Z.Reallocate(order+2, order+2, order);
    NumDofs_X.Fill(-1); NumDofs_Y.Fill(-1); NumDofs_Z.Fill(-1);
    
    nb_dof_loc = 0;
    const Matrix<int>& CoordinateNodes = this->Fb_geom.GetCoordinateNodes3D();
    int order_geom = this->GetGeometryOrder();
    // nodes on edges
    for (int ne = 0; ne < 12 ; ne++)
      {
	// we get the two extremities of the edge
	int num_vertex_1 = MeshNumbering<Dimension3>::FirstExtremityEdge(3, ne);
	int num_vertex_2 = MeshNumbering<Dimension3>::SecondExtremityEdge(3, ne);
	// we get the (i,j,k) coordinate of the first vertex
	int i = CoordinateNodes(num_vertex_1,0);
	int j = CoordinateNodes(num_vertex_1,1);
	int k = CoordinateNodes(num_vertex_1,2);
	
	int j2 = CoordinateNodes(num_vertex_2,1);
	int k2 = CoordinateNodes(num_vertex_2,2);
        if (i == order_geom)
          i = order+1;
        
        if (j == order_geom)
          j = order+1;
        
        if (k == order_geom)
          k = order+1;

        if (j2 == order_geom)
          j2 = order+1;
        
        if (k2 == order_geom)
          k2 = order+1;
        
	if (k2 != k)
	  for (int l = 0 ; l < order ; l++)
	    {
	      // we have an edge along Oz direction
	      NumDofs_Z(i,j,l) = nb_dof_loc;
	      nb_dof_loc++;
	    }
	else if (j2 != j)
	  for (int l = 0 ; l < order ; l++)
	    {
	      // we have an edge along Oy direction
	      NumDofs_Y(i,l,k) = nb_dof_loc;
	      nb_dof_loc++;
	    }
	else
	  for (int l = 0; l < order; l++)
	    {
	      // we have an edge along Ox direction
	      NumDofs_X(l,j,k) = nb_dof_loc;
	      nb_dof_loc++;
	    }
      }
    
    nb_dof_loc = 3*order*(order+2)*(order+2);
    
    // dofs on first face
    int offset = 8*order;
    for (int j = 1 ; j <= order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_Y(0, i, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(0, j, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on second face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j <= order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, 0, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(j, 0, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on third face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j <= order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, j, 0) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Y(j, i, 0) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on fourth face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j <= order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, j, order+1) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Y(j, i, order+1) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on fifth face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j <= order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, order+1, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(j, order+1, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on sixth face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j <= order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_Y(order+1, i, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(order+1, j, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}    
    
    // dofs inside the hexahedron
    offset += Fb_quad.GetNbDof();
    for (int i = 0 ; i < order ; i++)
      for (int j = 1 ; j <= order ; j++)
	for (int k = 1 ; k <= order ; k++)
	  {
	    NumDofs_X(i, j, k) = offset++;
	    NumDofs_Y(j, i, k) = offset++;
	    NumDofs_Z(j, k, i) = offset++;
	  }
    
    VectReal_wp points_dof1d; VectR3 points_dof3d;
    points_dof3d = this->PointsND();
    
    points_dof1d.Reallocate(order);
    for (int i = 0; i < order; i++)
      points_dof1d(i) = lob_basis.Points(i+1);
    
    lob_int.AffectPoints(points_dof1d);    

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(points_dof3d);
  }
  
  
  //! computation of some arrays, G_GL, GL_G and dGL_G
  //! needed for the computation of stiffness matrix
  void HexahedronHcurlOptimalFirstFamily::ConstructStiffnessMatrix()
  {
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    dGL_GL = lob_basis.GradPhi();
    
    stiff1d.Reallocate(order+2, order+2);
    stiff1d.Fill(0);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order+1; j++)
        for (int k = 0; k <= order+1; k++)
          stiff1d(i, j) += this->Weights1D(k)*dGL_GL(i, k)*dGL_GL(j, k);
    
    Matrix<Real_wp, General, ArrayRowSparse>
      Rh(3*nb_points_quadrature_inside, 3*nb_points_quadrature_inside);
    
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  for (int m = 0; m <= order+1; m++)
	    {
	      int col = 3*NumQuad3D(i,j,k);
	      int row = 3*NumQuad3D(i,j,m)+1;
	      Rh.AddInteraction(row, col, -dGL_GL(m,k));
	      	      
	      row = 3*NumQuad3D(i,m,k)+2;
	      Rh.AddInteraction(row, col, dGL_GL(m,j));
	      	      
	      col = 3*NumQuad3D(i,j,k)+1;
	      row = 3*NumQuad3D(i,j,m);
	      Rh.AddInteraction(row, col, dGL_GL(m,k));
	      	      
	      row = 3*NumQuad3D(m,j,k)+2;
	      Rh.AddInteraction(row, col, -dGL_GL(m,i));
	      	      
	      col = 3*NumQuad3D(i,j,k)+2;
	      row = 3*NumQuad3D(i,m,k);
	      Rh.AddInteraction(row, col, -dGL_GL(m,j));
	      	      
	      row = 3*NumQuad3D(m,j,k)+1;
	      Rh.AddInteraction(row, col, dGL_GL(m,i));
            }
    
    Seldon::Copy(Rh, rh_loc);
    
  }
  
  
  //! mass matrix
  void HexahedronHcurlOptimalFirstFamily::ConstructMassMatrix()
  {
    ValGaussExt.Reallocate(order);
    for (int i = 0; i < order; i++)
      ValGaussExt(i).Init(lob_int.EvaluatePhi(i, 0.0), lob_int.EvaluatePhi(i, 1.0));
    
    int N = nb_points_quadrature_inside;
    // constructing operator for interpolation from dof points to Gauss-Lobatto points
    Matrix<Real_wp, General, ArrayRowSparse> Ch(nb_dof_loc, 3*N);
    
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  for (int m = 0; m < order; m++)
	    {
	      int node = NumQuad3D(i, j, k);
	      int num_dof = NumDofs_X(m, j, k);
              
              if (i == 0)
                Ch.AddInteraction(num_dof, 3*node, ValGaussExt(m)(0));
              else if (i == order+1)
                Ch.AddInteraction(num_dof, 3*node, ValGaussExt(m)(1));
              else if (i == m+1)
                Ch.AddInteraction(num_dof, 3*node, 1.0);
              
	      num_dof = NumDofs_Y(i, m, k);
	      if (j == 0)
                Ch.AddInteraction(num_dof, 3*node+1, ValGaussExt(m)(0));
              else if (j == order+1)
                Ch.AddInteraction(num_dof, 3*node+1, ValGaussExt(m)(1));
              else if (j == m+1)
                Ch.AddInteraction(num_dof, 3*node+1, 1.0);
	      	      
	      num_dof = NumDofs_Z(i, j, m);
	      if (k == 0)
                Ch.AddInteraction(num_dof, 3*node+2, ValGaussExt(m)(0));
              else if (k == order+1)
                Ch.AddInteraction(num_dof, 3*node+2, ValGaussExt(m)(1));
              else if (k == m+1)
                Ch.AddInteraction(num_dof, 3*node+2, 1.0);
              
            }

    WeightsMassGL.Reallocate(order+2, order+2);
    invWeightsMassGL.Reallocate(order+2, order+2);
    invSqrtWeightsMassGL.Reallocate(order+2, order+2);
    for (int j = 0; j <= order+1; j++)
      for (int k = 0; k <= order+1; k++)
        {
          WeightsMassGL(j, k) = lob_basis.Weights(j)*lob_basis.Weights(k);
          invWeightsMassGL(j, k) = 1.0/(lob_basis.Weights(j)*lob_basis.Weights(k));
          invSqrtWeightsMassGL(j, k) = 1.0/sqrt(WeightsMassGL(j, k));
        }
    
    Seldon::Copy(Ch, ch_loc);    
    GLint_GL.Reallocate(order, order+2);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        GLint_GL(i, j) = lob_int.EvaluatePhi(i, lob_basis.Points(j));
    
    mass_matrix.Reallocate(order, order);
    mass_matrix.Fill(0);
    for (int i = 0; i < order; i++)
      for (int j = i; j < order; j++)
        for (int k = 0; k <= order+1; k++)
          mass_matrix(i, j) += lob_basis.Weights(k)*GLint_GL(i, k)*GLint_GL(j, k);
    
    mass_matrix_chol = mass_matrix;
    GetCholesky(mass_matrix_chol);
    
    const_grad_matrix.Reallocate(order, order+2);
    const_grad_matrix.Fill(0);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        for (int k = 0; k <= order+1; k++)
          const_grad_matrix(i, j) += lob_basis.Weights(k)*GLint_GL(i, k)*dGL_GL(j, k);
        
    Matrix<int> coor;
    MeshNumbering<Dimension2>::
      ConstructQuadrilateralNumbering(order, NumQuad2D_gauss, coor);    
    
    VectReal_wp points_gauss, omega;
    ComputeGaussLegendre(points_gauss, omega, order);
    GL_GX.Reallocate(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
        GL_GX(i, j) = lob_basis.EvaluatePhi(i, points_gauss(j));
    
    int order_geom = this->GetGeometryOrder();
    Matrix<Real_wp> G_Geom(order, order+1), GL_Geom(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order_geom; j++)
        GL_Geom(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order_geom; j++)
        G_Geom(i, j) = lob_int.EvaluatePhi(i, this->PointsNodal1D(j));
    
    const Globatto<Real_wp>& lob_geom = this->GetNodalShapeFunctions1D();
    Geom_GL.Reallocate(order_geom+1, order+2);
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order+1; j++)
        Geom_GL(i, j) = lob_geom.EvaluatePhi(i, lob_basis.Points(j));
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch1.Reallocate(3*order*(order_geom+1)*(order_geom+1), 3*this->GetNbPointsNodalElt());
    Ch2.Reallocate(3*order*(order+2)*(order_geom+1), 3*order*(order_geom+1)*(order_geom+1));
    Ch3.Reallocate(3*order*(order+2)*(order+2), 3*order*(order+2)*(order_geom+1));
    const Array3D<int>& NumNodes3D = this->Fb_geom.GetNumNodes3D();
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k <= order_geom; k++)
	  for (int m = 0; m < order; m++)
	    {
	      int node = 3*NumNodes3D(i, j, k);
	      int irow = 3*((order_geom+1)*( m*(order_geom+1) + j) + k);
	      Ch1.AddInteraction(irow, node, G_Geom(m, i));
	      	      
              irow = 3*((order_geom+1)*( i*order + m) + k);
	      Ch1.AddInteraction(irow+1, node+1, G_Geom(m, j));
	      	      
	      irow = 3*(order*( i*(order_geom+1) + j) + m);
	      Ch1.AddInteraction(irow+2, node+2, G_Geom(m, k));
            }
    
    
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k < order; k++)
	  for (int m = 0; m <= order+1; m++)
	    {
              int icol = (order_geom+1)*(k*(order_geom+1)+i) + j;
              int irow = (order_geom+1)*(k*(order+2) + m) + j;
	      Ch2.AddInteraction(3*irow, 3*icol, GL_Geom(m, i));
	      	      
              icol = (order_geom+1)*(i*order + k) + j;
              irow = (order_geom+1)*(m*order + k) + j;
	      Ch2.AddInteraction(3*irow+1, 3*icol+1, GL_Geom(m, i));
	      	      
	      icol = order*(i*(order_geom+1) + j) + k;
              irow = order*(m*(order_geom+1) + j) + k;
	      Ch2.AddInteraction(3*irow+2, 3*icol+2, GL_Geom(m, i));
            }

    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k < order; k++)
	  for (int m = 0; m <= order+1; m++)
	    {    
              int icol = (order_geom+1)*(k*(order+2) + i) + j;
              int irow = NumDofs_X(k, i, m);
	      Ch3.AddInteraction(irow, 3*icol, GL_Geom(m, j));
	      	      
	      icol = (order_geom+1)*(i*order + k) + j;
              irow = NumDofs_Y(i, k, m);
	      Ch3.AddInteraction(irow, 3*icol+1, GL_Geom(m, j));
	      	      
	      icol = order*(i*(order_geom+1) + j) + k;
              irow = NumDofs_Z(i, m, k);
	      Ch3.AddInteraction(irow, 3*icol+2, GL_Geom(m, j));
            }
    
    Seldon::Copy(Ch1, ch1_node);
    Seldon::Copy(Ch2, ch2_node);
    Seldon::Copy(Ch3, ch3_node);
  }

  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    contrib.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  {
	    contrib(NumDofs_X(i, j, k)) = feval(3*NumQuad3D(i+1, j, k));
	    contrib(NumDofs_Y(j, i, k)) = feval(3*NumQuad3D(j, i+1, k)+1);
	    contrib(NumDofs_Z(j, k, i)) = feval(3*NumQuad3D(j, k, i+1)+2);
	  }
  }
  
  
  //! projection from values on quadrature points to dof components
  /*!
    \param[in] Equad values on quadrature points
    \param[out] Edof dof components
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const
  {
    Edof.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  {
	    // basic interpolation
	    int num_dof = NumDofs_X(i, j, k);
	    Edof(num_dof) = Equad(3*NumQuad3D(i+1, j, k));
	    
	    num_dof = NumDofs_Y(j, i, k);
	    Edof(num_dof) = Equad(3*NumQuad3D(j, i+1, k)+1);
	    
	    num_dof = NumDofs_Z(k, j, i);
	    Edof(num_dof) = Equad(3*NumQuad3D(k, j, i+1)+2);
	  }
  }
  
  
  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] feval vector containing \omega_k f(\xi_k) 
    \param[out] res res_i = res_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& feval, Vector2& res, int r) const
  {
    if ((r != 0) && (r != order))
      {
        abort();
      }
    
    //typedef typename Vector2::value_type T2;
    int ext = 0, num_dof;
    if (num_loc > 2)
      ext = order+1;

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	// Face x = 0 or Face x = 1  (ext = 0 or 1)
	
        // NumDofs_X(i, j, k)
        // integration of \int phi_i^G(ext) phi_j^{GL}(y) phi_k^{GL}(z)  f_x
        //                = phi_i^G(ext) f(ext, \xi_j, \xi_k)
 
	// NumDofs_Y(j, i, k)
	// integration of \int \phi_j^GL(ext) \phi_i^G(y)  \phi_k^{GL}(z) f_y
	// = \sum_m G_GL(i,m) \delta_{ext,j} f_y(ext, \xi_m^{GL}, \xi_k^{GL})
	
	// NumDofs_Z(k, j, i)
	// integration of \int \phi_k^{GL}(ext)  \phi_j^{GL}(y) \phi_i^G(z) f_z
	// = \sum_n G_GL(i,n) \delta{k, ext} f_z (ext, \xi_j^{GL}, \xi_n^{GL})
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                res(NumDofs_X(i, j, k)) += alpha*GLint_GL(i, ext)*feval(3*NumQuad2D(j, k));
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Y(j, i, k);
                    for (int m = 0; m <= order+1; m++)
                      res(num_dof) += alpha*GLint_GL(i, m)*feval(3*NumQuad2D(m, k)+1);
                  }

                if (k == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order+1; n++)
                      res(num_dof) += alpha*GLint_GL(i, n)*feval(3*NumQuad2D(j, n)+2);
                  }
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
        // face y = 0 or y = 1 (ext = 0 or 1)
        
        // NumDofs_X(i, j, k)
        // integration of \int phi_i^G(x) phi_j^{GL}(ext) phi_k^{GL}(z)  f_x
        //                = \sum_m G_GL(i, m) \delta_{j,ext} f_x(\xi_m, ext, \xi_k)
 
	// NumDofs_Y(j, i, k)
	// integration of \int \phi_j^GL(x) \phi_i^G(ext)  \phi_k^{GL}(z) f_y
	// = G_GL(i, ext) f_y(\xi_j^{GL}, ext, \xi_k^{GL})
	
	// NumDofs_Z(k, j, i)
	// integration of \int \phi_k^{GL}(x)  \phi_j^{GL}(ext) \phi_i^G(z) f_z
	// = \sum_n G_GL(i,n) \delta{j, ext} f_z (\xi_k^{GL}, ext, \xi_n^{GL})
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                if (j == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order+1; m++)
                      res(num_dof) += alpha*GLint_GL(i, m)*feval(3*NumQuad2D(m, k));
                  }
                
                num_dof = NumDofs_Y(j, i, k);
                res(num_dof) += alpha*GLint_GL(i, ext)*feval(3*NumQuad2D(j, k)+1);
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order+1; n++)
                      res(num_dof) += alpha*GLint_GL(i, n)*feval(3*NumQuad2D(k, n)+2);
                  }
              }
      }
    else
      {
        // face z = 0 or z = 1 (ext = 0 or 1)
        
        // NumDofs_X(i, j, k)
        // integration of \int phi_i^G(x) phi_j^{GL}(y) phi_k^{GL}(ext)  f_x
        //                = \sum_m G_GL(i, m) \delta_{k,ext} f_x(\xi_m, \xi_j, 0)
 
	// NumDofs_Y(j, i, k)
	// integration of \int \phi_j^GL(x) \phi_i^G(ext)  \phi_k^{GL}(ext) f_y
	// = \sum_n G_GL(i, n) \delta_{k, ext} f_y(\xi_j^{GL}, \xi_n, ext)
	
	// NumDofs_Z(k, j, i)
	// integration of \int \phi_k^{GL}(x)  \phi_j^{GL}(y) \phi_i^G(ext) f_z
	// = G_GL(i, ext) f_z (\xi_k^{GL}, \xi_j^{GL}, ext)
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                if (k == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order+1; m++)
                      res(num_dof) += alpha*GLint_GL(i, m)*feval(3*NumQuad2D(m, j));
                    
                    num_dof = NumDofs_Y(j, i, k);
                    for (int n = 0; n <= order+1; n++)
                      res(num_dof) += alpha*GLint_GL(i, n)*feval(3*NumQuad2D(j, n)+1);
                  }
                
                num_dof = NumDofs_Z(k, j, i);
                res(num_dof) += alpha*GLint_GL(i, ext)*feval(3*NumQuad2D(k, j)+2);
              }
      }
    
  }
  
  
  //! computation of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] res components of u on degrees of freedom
    \param[out] feval values of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    feval = Sh* res
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ApplyShTransposeGen(int num_loc, const Vector1& res, Vector2& feval, int r) const
  {
    feval.Fill(0);
    
    if ((r != 0) && (r != order))
      {
        abort();
      }
    
    //typedef typename Vector2::value_type T2;
    int ext = 0, num_dof;
    if (num_loc > 2)
      ext = order+1;

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                feval(3*NumQuad2D(j, k)) += GLint_GL(i, ext)*res(NumDofs_X(i, j, k));
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Y(j, i, k);
                    for (int m = 0; m <= order+1; m++)
                      feval(3*NumQuad2D(m, k)+1) += GLint_GL(i, m)*res(num_dof);
                  }

                if (k == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order+1; n++)
                      feval(3*NumQuad2D(j, n)+2) += GLint_GL(i, n)*res(num_dof);
                  }
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                if (j == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order+1; m++)
                      feval(3*NumQuad2D(m, k)) += GLint_GL(i, m)*res(num_dof);
                  }
                
                num_dof = NumDofs_Y(j, i, k);
                feval(3*NumQuad2D(j, k)+1) += GLint_GL(i, ext)*res(num_dof);
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order+1; n++)
                      feval(3*NumQuad2D(k, n)+2) += GLint_GL(i, n)*res(num_dof);
                  }
              }
      }
    else
      {
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                if (k == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order+1; m++)
                      feval(3*NumQuad2D(m, j)) += GLint_GL(i, m)*res(num_dof);
                    
                    num_dof = NumDofs_Y(j, i, k);
                    for (int n = 0; n <= order+1; n++)
                      feval(3*NumQuad2D(j, n)+1) += GLint_GL(i, n)*res(num_dof);
                  }
                
                num_dof = NumDofs_Z(k, j, i);
                feval(3*NumQuad2D(k, j)+2) += GLint_GL(i, ext)*res(num_dof);
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
    (\omega, \xi_k) are here Gauss points and not Gauss-Lobatto points
  */  
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    Vector1 fx(3*(order+1)*(order+2)), fy(3*(order+2)*(order+2));
    FillZero(fx); FillZero(fy);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
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
  
  
  //! integration against curl of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] feval vector containing \omega_k f(\xi_k) 
    \param[out] res res_i = res_i + alpha \int_{\partial K} f curl(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
  */
  template<class T0, class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& feval, Vector2& res, int r) const
  {
    if ((r != 0) && (r != order))
      {
        abort();
      }
    
    typedef typename Vector2::value_type T2;
    T2 vloc;
    int ext = 0, num_dof;
    if (num_loc > 2)
      ext = order+1;
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	// Face x = 0 or x = 1	
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                // contribution for curl phi_i,j,k^x
                vloc = T2(0);
                for (int n = 0; n <= order+1; n++)
                  vloc += dGL_GL(k, n)*feval(3*NumQuad2D(j, n)+1);
		                
                for (int m = 0; m <= order+1; m++)
                  vloc -= dGL_GL(j, m)*feval(3*NumQuad2D(m, k)+2);
                
                res(NumDofs_X(i, j, k) ) += alpha*GLint_GL(i, ext)*vloc;
               
                // contribution for curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k); 
                if (j == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        vloc -= GLint_GL(i, m)*dGL_GL(k, n)*feval(3*NumQuad2D(m, n));
                    
                    res(num_dof) += alpha*vloc;
                  }
                    
                vloc = T2(0);
                for (int m = 0; m <= order+1; m++)
                  vloc += GLint_GL(i, m)*feval(3*NumQuad2D(m, k)+2);
		                
                res(num_dof) += alpha*dGL_GL(j, ext)*vloc;
                
                // contribution for curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (k == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        vloc += GLint_GL(i, n)*dGL_GL(j, m)*feval(3*NumQuad2D(m, n));
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int n = 0; n <= order+1; n++)
                  vloc += GLint_GL(i, n)*feval(3*NumQuad2D(j, n)+1);
                
                res(num_dof) -= alpha*dGL_GL(k, ext)*vloc;
                
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
        // Face y = 0 or y = 1
         for (int i = 0; i < order; i++)
          for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (j == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        vloc += GLint_GL(i, m)*dGL_GL(k, n)*feval(3*NumQuad2D(m, n)+1);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int m = 0; m <= order+1; m++)
                  vloc += GLint_GL(i, m)*feval(3*NumQuad2D(m, k)+2);
                
                res(num_dof) -= alpha*dGL_GL(j, ext)*vloc;
                              
                // contribution for curl phi_j,i,k^y 
                vloc = T2(0) ;
                for (int n = 0; n <= order+1; n++)
                  vloc -= dGL_GL(k, n)*feval(3*NumQuad2D(j, n));
                
                for (int m = 0; m <= order+1; m++)
                  vloc += dGL_GL(j, m)*feval(3*NumQuad2D(m, k)+2);
                
                res(NumDofs_Y(j, i, k)) += alpha*GLint_GL(i, ext)*vloc;
                
                //contribution for  curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (j == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        vloc -= GLint_GL(i, n)*dGL_GL(k, m)*feval(3*NumQuad2D(m, n)+1);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int n = 0; n <= order+1; n++)
                  vloc += GLint_GL(i, n)*feval(3*NumQuad2D(k, n));
                
                res(num_dof) += alpha*dGL_GL(j, ext)*vloc;
              }
      }
    else if ((num_loc == 2)||(num_loc == 3))
      {
        // Face z = 0 or z = 1
         
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {

                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (k == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        vloc -= GLint_GL(i, m)*dGL_GL(j, n)*feval(3*NumQuad2D(m, n)+2);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int m = 0; m <= order+1; m++)
                  vloc += GLint_GL(i, m)*feval(3*NumQuad2D(m, j)+1);
                
                res(num_dof) += alpha*dGL_GL(k, ext)*vloc;
                
                // contribution for  curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k);
                if (k == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        vloc += GLint_GL(i, n)*dGL_GL(j, m)*feval(3*NumQuad2D(m, n)+2);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int n = 0; n <= order+1; n++)
                  vloc += GLint_GL(i, n)*feval(3*NumQuad2D(j, n));
                
                res(num_dof) -= alpha*dGL_GL(k, ext)*vloc;
                
                // contribution for curl phi_k,j,i^z 
                vloc = T2(0) ;
                for (int n = 0; n <= order+1; n++)
                  vloc += dGL_GL(j, n)*feval(3*NumQuad2D(k, n));
                
                for (int m = 0; m <= order+1; m++)
                  vloc -= dGL_GL(k, m)*feval(3*NumQuad2D(m, j)+1);
                
                res(NumDofs_Z(k, j, i)) += alpha*GLint_GL(i, ext)*vloc;
                
              }
      }
    
  }
  
  
  //! computation of curl of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] res components of u on degrees of freedom
    \param[out] feval curl of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    feval = Sh* res
    where (Sh)_{i,j} = curl \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ApplyNablaShTransposeGen(int num_loc, const Vector1& res, Vector2& feval, int r) const
  {
    if ((r != 0) && (r != order))
      {
        abort();
      }
    
    typedef typename Vector2::value_type T2;
    T2 vloc;
    int ext = 0, num_dof;
    if (num_loc > 2)
      ext = order+1;
    
    feval.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	// Face x = 0 or x = 1	
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                // contribution for curl phi_i,j,k^x
                vloc = GLint_GL(i, ext)*res(NumDofs_X(i, j, k) );

                for (int n = 0; n <= order+1; n++)
                  feval(3*NumQuad2D(j, n)+1) += dGL_GL(k, n)*vloc;
		                
                for (int m = 0; m <= order+1; m++)
                  feval(3*NumQuad2D(m, k)+2) -= dGL_GL(j, m)*vloc;
                               
                // contribution for curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k); 
                if (j == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        feval(3*NumQuad2D(m, n)) -= GLint_GL(i, m)*dGL_GL(k, n)*vloc;
                  }
                    
                vloc = dGL_GL(j, ext)*res(num_dof);
                for (int m = 0; m <= order+1; m++)
                  feval(3*NumQuad2D(m, k)+2) += GLint_GL(i, m)*vloc;
                
                // contribution for curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (k == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        feval(3*NumQuad2D(m, n)) += GLint_GL(i, n)*dGL_GL(j, m)*vloc;
                    
                  }
                
                vloc = -dGL_GL(k, ext)*res(num_dof);
                for (int n = 0; n <= order+1; n++)
                  feval(3*NumQuad2D(j, n)+1) += GLint_GL(i, n)*vloc;
                
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
        // Face y = 0 or y = 1
         for (int i = 0; i < order; i++)
          for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {
                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (j == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        feval(3*NumQuad2D(m, n)+1) += GLint_GL(i, m)*dGL_GL(k, n)*vloc;
                  }
                
                vloc = -dGL_GL(j, ext)*res(num_dof);
                for (int m = 0; m <= order+1; m++)
                  feval(3*NumQuad2D(m, k)+2) += GLint_GL(i, m)*vloc;
                              
                // contribution for curl phi_j,i,k^y 
                vloc = GLint_GL(i, ext)*res(NumDofs_Y(j, i, k));
                for (int n = 0; n <= order+1; n++)
                  feval(3*NumQuad2D(j, n)) -= dGL_GL(k, n)*vloc;
                
                for (int m = 0; m <= order+1; m++)
                  feval(3*NumQuad2D(m, k)+2) += dGL_GL(j, m)*vloc;
                                
                // contribution for  curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (j == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        feval(3*NumQuad2D(m, n)+1) -= GLint_GL(i, n)*dGL_GL(k, m)*vloc;
                  }
                
                vloc = dGL_GL(j, ext)*res(num_dof);
                for (int n = 0; n <= order+1; n++)
                  feval(3*NumQuad2D(k, n)) += GLint_GL(i, n)*vloc;
              }
      }
    else if ((num_loc == 2)||(num_loc == 3))
      {
        // Face z = 0 or z = 1
         
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order+1; j++)
            for (int k = 0; k <= order+1; k++)
              {

                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (k == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        feval(3*NumQuad2D(m, n)+2) -= GLint_GL(i, m)*dGL_GL(j, n)*vloc;
                  }
                
                vloc = dGL_GL(k, ext)*res(num_dof);
                for (int m = 0; m <= order+1; m++)
                  feval(3*NumQuad2D(m, j)+1) += GLint_GL(i, m)*vloc;
                
                // contribution for  curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k);
                if (k == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order+1; m++)
                      for (int n = 0; n <= order+1; n++)
                        feval(3*NumQuad2D(m, n)+2) += GLint_GL(i, n)*dGL_GL(j, m)*vloc;
                    
                  }
                
                vloc = -dGL_GL(k, ext)*res(num_dof);
                for (int n = 0; n <= order+1; n++)
                  feval(3*NumQuad2D(j, n)) += GLint_GL(i, n)*vloc;
                
                // contribution for curl phi_k,j,i^z 
                vloc = GLint_GL(i, ext)*res(NumDofs_Z(k, j, i));
                for (int n = 0; n <= order+1; n++)
                  feval(3*NumQuad2D(k, n)) += dGL_GL(j, n)*vloc;
                
                for (int m = 0; m <= order+1; m++)
                  feval(3*NumQuad2D(m, j)+1) -= dGL_GL(k, m)*vloc;
                
              }
      }
    
  }

  
  //! computation of U on nodal points (on reference element)
  /*!
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ComputeNodalValuesGen(const Vector1& Uloc, Vector2& Uloc_node) const
  {
    Vector2 Uy(ch3_node.GetN()), Ux(ch2_node.GetN());
    Mlt(SeldonTrans, ch3_node, Uloc, Uy);
    Mlt(SeldonTrans, ch2_node, Uy, Ux);
    Mlt(SeldonTrans, ch1_node, Ux, Uloc_node);
  }
  
    
  //! computation of u on nodal points of a face
  /*!
    \param[in] u_loc dof components of u
    \param[out] u_boundary values of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ComputeValueBoundaryGen(const Vector1& u_loc, Vector2& u_boundary, int num_loc) const
  {
    typedef typename Vector1::value_type Complexe;
    Complexe Ex, Ey, Ez;
    u_boundary.Reallocate(3*(order+2)*(order+2));
    int ext = 0;
    if (num_loc > 2)
      ext = order+1;
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc==5))
      {
	// face x = 0 or x = 1
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      // interpolation
	      Ey = 0; Ez = 0; Ex = 0;
	      for (int k = 0; k < order; k++)
		{
		  Ex += GLint_GL(k,ext)*u_loc(NumDofs_X(k,i,j));
		  Ey += GLint_GL(k, i)*u_loc(NumDofs_Y(ext,k,j));
		  Ez += GLint_GL(k, j)*u_loc(NumDofs_Z(ext,i,k));
		}
	      
	      // storing the value
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
    else if ((num_loc == 1)||(num_loc==4))
      {
	// face y = 0 or y = 1
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      // dofs along y and z
	      Ey = 0; Ez = 0; Ex = 0;
	      for (int k = 0; k < order; k++)
		{
		  Ex += GLint_GL(k, i)*u_loc(NumDofs_X(k,ext,j));
		  Ey += GLint_GL(k, ext)*u_loc(NumDofs_Y(i, k, j));
		  Ez += GLint_GL(k, j)*u_loc(NumDofs_Z(i,ext,k));
		}
	      
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
    else if ((num_loc == 2)||(num_loc==3))
      {
	// face z = 0 or z = 1
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      // dofs along x, y and z
	      Ey = 0; Ez = 0; Ex = 0;
	      for (int k = 0; k < order; k++)
		{
		  Ex += GLint_GL(k,i)*u_loc(NumDofs_X(k,j,ext));
		  Ey += GLint_GL(k,j)*u_loc(NumDofs_Y(i,k,ext));
		  Ez += GLint_GL(k,ext)*u_loc(NumDofs_Z(i,j,k));
		}
	      
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
    
  }
  
  
  //! computation of curl(u) on nodal points of a face
  /*!
    \param[in] u_loc dof components of u
    \param[out] u_boundary curl of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlOptimalFirstFamily::
  ComputeCurlBoundaryGen(const Vector1& u_loc,
			 Vector2& u_boundary, int num_loc) const
  {
    typedef typename Vector1::value_type Complexe;
    Complexe Ex, Ey, Ez;
    u_boundary.Reallocate(3*(order+2)*(order+2));
    int ext = 0;
    if (num_loc > 2)
      ext = order+1;

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
    if ((num_loc == 0)||(num_loc==5))
      {
	// face x = 0 or x = 1
	Matrix<Complexe> Ex_GL(order+2,order+2), Ey_GL(order+2,order+2),
	  Ez_GL(order+2,order+2);
	Matrix<Complexe> dEy_G(order,order+2), dEz_G(order,order+2);
	Ex_GL.Fill(0); Ey_GL.Fill(0); Ez_GL.Fill(0); dEy_G.Fill(0); dEz_G.Fill(0);
	// first step, interpolation at GL points
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      for (int k = 0; k < order; k++)
		{
		  Ex_GL(i,j) += GLint_GL(k,ext)*u_loc(NumDofs_X(k,i,j));
		  Ey_GL(i,j) += GLint_GL(k,i)*u_loc(NumDofs_Y(ext,k,j));
		  Ez_GL(i,j) += GLint_GL(k,j)*u_loc(NumDofs_Z(ext,i, k));
		  dEy_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Y(j,k,i));
		  dEz_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Z(j,i,k));
		}
	    }
	
	// second step : derivative
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      Ey = 0; Ex = 0; Ez = 0;
	      for (int k = 0; k <= order+1; k++)
		{
		  Ey += dGL_GL(k,j)*Ex_GL(i,k);
		  Ez -= dGL_GL(k,i)*Ex_GL(k,j);
		  Ex -= dGL_GL(k,j)*Ey_GL(i,k);
		  Ex += dGL_GL(k,i)*Ez_GL(k,j);
		}
	      
	      for (int k = 0; k < order; k++)
		{
		  Ez += GLint_GL(k,i)*dEy_G(k,j);
		  Ey -= GLint_GL(k,j)*dEz_G(k,i);
		}
	      // storing the value
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
    else if ((num_loc == 1)||(num_loc==4))
      {
	// face y = 0 or y = 1
	Matrix<Complexe> Ex_GL(order+2,order+2), Ey_GL(order+2,order+2),
	  Ez_GL(order+2,order+2);
	Matrix<Complexe> dEx_G(order,order+2), dEz_G(order,order+2);
	Ex_GL.Fill(0); Ey_GL.Fill(0); Ez_GL.Fill(0); dEx_G.Fill(0); dEz_G.Fill(0);
	// first step, interpolation at GL points
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      for (int k = 0; k < order; k++)
		{
		  Ex_GL(i,j) += GLint_GL(k,i)*u_loc(NumDofs_X(k,ext,j));
		  Ey_GL(i,j) += GLint_GL(k,ext)*u_loc(NumDofs_Y(i,k,j));
		  Ez_GL(i,j) += GLint_GL(k,j)*u_loc(NumDofs_Z(i,ext,k));
		  dEx_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_X(k,j,i));
		  dEz_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Z(i,j,k));
		}
	    }
	
	// second step : derivative
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      Ey = 0; Ex = 0; Ez = 0;
	      for (int k = 0; k <= order+1; k++)
		{
		  Ey += dGL_GL(k,j)*Ex_GL(i,k);
		  Ex -= dGL_GL(k,j)*Ey_GL(i,k);
		  Ez += dGL_GL(k,i)*Ey_GL(k,j);
		  Ey -= dGL_GL(k,i)*Ez_GL(k,j);
		}
	      
	      for (int k = 0; k < order; k++)
		{
		  Ez -= GLint_GL(k,i)*dEx_G(k,j);
		  Ex += GLint_GL(k,j)*dEz_G(k,i);
		}
	      
	      // storing the value
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
    else
      {
	// face z = 0 or z = 1
	Matrix<Complexe> Ex_GL(order+2,order+2), Ey_GL(order+2,order+2),
	  Ez_GL(order+2,order+2);
	Matrix<Complexe> dEx_G(order,order+2), dEy_G(order,order+2);
	Ex_GL.Fill(0); Ey_GL.Fill(0); Ez_GL.Fill(0); dEx_G.Fill(0); dEy_G.Fill(0);
	// first step, interpolation at GL points
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      for (int k = 0; k < order; k++)
		{
		  Ex_GL(i,j) += GLint_GL(k,i)*u_loc(NumDofs_X(k,j,ext));
		  Ey_GL(i,j) += GLint_GL(k,j)*u_loc(NumDofs_Y(i,k,ext));
		  Ez_GL(i,j) += GLint_GL(k,ext)*u_loc(NumDofs_Z(i,j,k));
		  dEx_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_X(k,i,j));
		  dEy_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Y(i,k,j));
		}
	    }
	
	// second step : derivative
	for (int i = 0; i <= order+1; i++)
	  for (int j = 0; j <= order+1; j++)
	    {
	      Ey = 0; Ex = 0; Ez = 0;
	      for (int k = 0; k <= order+1; k++)
		{
		  Ez -= dGL_GL(k,j)*Ex_GL(i,k);
		  Ez += dGL_GL(k,i)*Ey_GL(k,j);
		  Ex += dGL_GL(k,j)*Ez_GL(i,k);
		  Ey -= dGL_GL(k,i)*Ez_GL(k,j);
		}
	      
	      for (int k = 0; k < order; k++)
		{
		  Ey += GLint_GL(k,i)*dEx_G(k,j);
		  Ex -= GLint_GL(k,j)*dEy_G(k,i);
		}
	      
	      // storing the value
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
  }
  
  
  //! x is overwritten by M^-1 x where M is the mass matrix  
  template<class Vector1>
  void HexahedronHcurlOptimalFirstFamily::SolveMassMatrixGen(Vector1& x) const
  {
    Vector1 xloc(order);
    for (int j = 0; j <= order+1; j++)
      for (int k = 0; k <= order+1; k++)
        {
          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_X(i, j, k));
          
          Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, xloc);
          Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, xloc);
          
          Real_wp alpha = invWeightsMassGL(j, k);
          for (int i = 0; i < order; i++)
            x(NumDofs_X(i, j, k)) = alpha*xloc(i);
          
          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_Y(j, i, k));
          
          Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, xloc);
          Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, xloc);
          
          for (int i = 0; i < order; i++)
            x(NumDofs_Y(j, i, k)) = alpha*xloc(i);
          
          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_Z(j, k, i));
          
          Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, xloc);
          Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, xloc);
          
          for (int i = 0; i < order; i++)
            x(NumDofs_Z(j, k, i)) = alpha*xloc(i);
          
        }          
  }

  
  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix  
  template<class Vector1>
  void HexahedronHcurlOptimalFirstFamily::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    Vector1 xloc(order);
    for (int j = 0; j <= order+1; j++)
      for (int k = 0; k <= order+1; k++)
        {
          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_X(i, j, k));
          
          Seldon::SolveCholesky(TransA, mass_matrix_chol, xloc);
          
          Real_wp alpha = invSqrtWeightsMassGL(j, k);
          for (int i = 0; i < order; i++)
            x(NumDofs_X(i, j, k)) = alpha*xloc(i);

          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_Y(j, i, k));
          
          Seldon::SolveCholesky(TransA, mass_matrix_chol, xloc);
          
          for (int i = 0; i < order; i++)
            x(NumDofs_Y(j, i, k)) = alpha*xloc(i);

          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_Z(j, k, i));
          
          Seldon::SolveCholesky(TransA, mass_matrix_chol, xloc);
          
          for (int i = 0; i < order; i++)
            x(NumDofs_Z(j, k, i)) = alpha*xloc(i);
        }          

  }
    
  
  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector2>
  void HexahedronHcurlOptimalFirstFamily::MltMassMatrixGen(Vector2& x) const
  {
    Vector2 xloc(order), yloc(order);
    yloc.Fill(0);
    for (int j = 0; j <= order+1; j++)
      for (int k = 0; k <= order+1; k++)
        {
          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_X(i, j, k));
          
          Mlt(mass_matrix, xloc, yloc);
          
          Real_wp alpha = WeightsMassGL(j, k);
          for (int i = 0; i < order; i++)
            x(NumDofs_X(i, j, k)) = alpha*yloc(i);

          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_Y(j, i, k));
          
          Mlt(mass_matrix, xloc, yloc);
          
          for (int i = 0; i < order; i++)
            x(NumDofs_Y(j, i, k)) = alpha*yloc(i);

          for (int i = 0; i < order; i++)
            xloc(i) = x(NumDofs_Z(j, k, i));
          
          Mlt(mass_matrix, xloc, yloc);
          
          for (int i = 0; i < order; i++)
            x(NumDofs_Z(j, k, i)) = alpha*yloc(i);
        }
  }

  
  //! Integration against curl of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \nabla \times \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = curl(\varphi_i)(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    Vector1 Uquad(3*this->nb_points_quadrature_inside);
    Seldon::Mlt(rh_loc, Vh, Uquad);
    Mlt(ch_loc, Uquad, Uh);
  }
  
  
  //! computation of curl of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh curl of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j curl phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = curl phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(3*this->nb_points_quadrature_inside);
    Mlt(SeldonTrans, ch_loc, Uh, Uquad);
    Seldon::Mlt(SeldonTrans, rh_loc, Uquad, Vh);
  }
  

  //! Integration against curl of basis functions associated with quadrature points
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f curl(\psi_i) dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = \nabla \times \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::ApplyRhQuadratureGen(const Vector1& Vh, Vector2& Uh) const
  {
    Seldon::Mlt(rh_loc, Vh, Uh);
  }
  
  
  //! Computation of curl on quadrature points from values on quadrature points
  /*!
    \param[in] Uh values of u on quadrature points
    \param[out] Vh curl of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Rh* Uh
    where (Rh)_{i,j} = \nabla \times \psi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily
  ::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Seldon::Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }
  
  
  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronHcurlOptimalFirstFamily::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    Seldon::Mlt(ch_loc, Vh, Uh);
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
  void HexahedronHcurlOptimalFirstFamily::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Seldon::Mlt(SeldonTrans, ch_loc, Uh, Vh);
  }
  
  
  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C tensor
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + M
    where M is the mass matrix, M_ij = \int C \varphi_j \varphi_i dx
   */
  template<class T, class Prop>
  void HexahedronHcurlOptimalFirstFamily::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    T vloc;
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        for (int k = 0; k <= order+1; k++)
          {
            int ix = NumDofs_X(i, j, k);
            int iy = NumDofs_Y(j, i, k);
            int iz = NumDofs_Z(j, k, i);
            
            // diagonal part
            for (int ip = 0; ip < order; ip++)
              {
                int jx = NumDofs_X(ip, j, k);
		A.AddInteraction(m+ix, n+jx, C(0, 0)*WeightsMassGL(j, k)*mass_matrix(i, ip));
                
                int jy = NumDofs_Y(j, ip, k);
		A.AddInteraction(m+iy, n+jy, C(1, 1)*WeightsMassGL(j, k)*mass_matrix(i, ip));
                
                int jz = NumDofs_Z(j, k, ip);
		A.AddInteraction(m+iz, n+jz, C(2, 2)*WeightsMassGL(j, k)*mass_matrix(i, ip));
              }
	    
            // extra-diagonal part
            for (int ip = 0; ip < order; ip++)
              for (int jp = 0; jp <= order+1; jp++)
                {
                  vloc = weights1d(jp)*weights1d(j)*GLint_GL(i, jp)*GLint_GL(ip, j)*weights1d(k);
                  ix = NumDofs_X(i, j, k);
                  int jy = NumDofs_Y(jp, ip, k);
		  A.AddInteraction(m+ix, n+jy, C(0, 1)*vloc);
		  A.AddInteraction(m+jy, n+ix, C(0, 1)*vloc);
                  
                  ix = NumDofs_X(i, k, j);
                  int jz = NumDofs_Z(jp, k, ip);
		  A.AddInteraction(m+ix, n+jz, C(0, 2)*vloc);  
		  A.AddInteraction(m+jz, n+ix, C(0, 2)*vloc);
		  
                  iy = NumDofs_Y(k, i, j);
                  jz = NumDofs_Z(k, jp, ip);		  
		  A.AddInteraction(m+iy, n+jz, C(1, 2)*vloc);
		  A.AddInteraction(m+jz, n+iy, C(1, 2)*vloc);                  
		}
          }    
  }
  
  
  //! we add constant stiffness matrix \int_K  C curl(phi_j) curl(phi_i)
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C tensor
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + \sum_p  C(p, q) S^{p,q}
    where S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K curl(\phi_j)_q curl(\phi_i)_p  dx
   */  
  template<class T, class Prop>
  void HexahedronHcurlOptimalFirstFamily::
  AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    T vloc;
    int ix, iy, iz, jx, jy, jz;
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        for (int k = 0; k <= order+1; k++)
          {
            for (int ip = 0; ip < order; ip++)
              for (int kp = 0; kp <= order+1; kp++)
                {
                  // diagonal part
                  iz = NumDofs_Z(j, k, i);
                  jz = NumDofs_Z(j, kp, ip);
		  A.AddInteraction(m+iz, n+jz, C(0, 0)*mass_matrix(i, ip)*weights1d(j)*stiff1d(k, kp));
		  
                  iy = NumDofs_Y(j, i, k);
                  jy = NumDofs_Y(j, ip, kp);
		  A.AddInteraction(m+iy, n+jy, C(0, 0)*mass_matrix(i, ip)*weights1d(j)*stiff1d(k, kp));

                  ix = NumDofs_X(i, j, k);
                  jx = NumDofs_X(ip, j, kp);
		  A.AddInteraction(m+ix, n+jx, C(1, 1)*mass_matrix(i, ip)*weights1d(j)*stiff1d(k, kp));
                
                  iz = NumDofs_Z(k, j, i);
                  jz = NumDofs_Z(kp, j, ip);
		  A.AddInteraction(m+iz, n+jz, C(1, 1)*mass_matrix(i, ip)*weights1d(j)*stiff1d(k, kp));
                
                  ix = NumDofs_X(i, k, j);
                  jx = NumDofs_X(ip, kp, j);
		  A.AddInteraction(m+ix, n+jx, C(2, 2)*mass_matrix(i, ip)*weights1d(j)*stiff1d(k, kp));
                
                  iy = NumDofs_Y(k, i, j);
                  jy = NumDofs_Y(kp, ip, j);
		  A.AddInteraction(m+iy, n+jy, C(2, 2)*mass_matrix(i, ip)*weights1d(j)*stiff1d(k, kp));
                
                  iz = NumDofs_Z(j, k, i);
                  jy = NumDofs_Y(j, ip, kp);
                  vloc = C(0, 0)*const_grad_matrix(i, kp)*const_grad_matrix(ip, k)*weights1d(j);
		  A.AddInteraction(m+iz, n+jy, -vloc);
		  A.AddInteraction(m+jy, n+iz, -vloc);
                  
                  ix = NumDofs_X(i, j, k);
                  jz = NumDofs_Z(kp, j, ip);
                  vloc = C(1, 1)*const_grad_matrix(i, kp)*const_grad_matrix(ip, k)*weights1d(j);
		  A.AddInteraction(m+ix, n+jz, -vloc);
		  A.AddInteraction(m+jz, n+ix, -vloc);
                  
                  ix = NumDofs_X(i, k, j);
                  jy = NumDofs_Y(kp, ip, j);
                  vloc = C(2, 2)*const_grad_matrix(i, kp)*const_grad_matrix(ip, k)*weights1d(j);
		  A.AddInteraction(m+ix, n+jy, -vloc);
		  A.AddInteraction(m+jy, n+ix, -vloc);
                }
            
            // extra-diagonal part
            ix = NumDofs_X(i, j, k);
            iy = NumDofs_Y(j, i, k);
            iz = NumDofs_Z(j, k, i);
            for (int ip = 0; ip < order; ip++)
              for (int jp = 0; jp <= order+1; jp++)
                for (int kp = 0; kp <= order+1; kp++)
                  {
                    jx = NumDofs_X(ip, jp, kp);
                    jy = NumDofs_Y(jp, ip, kp);
                    jz = NumDofs_Z(jp, kp, ip);
                    
                    // c10
                    vloc = C(0, 1)*GLint_GL(ip, j)*weights1d(j)*dGL_GL(k, jp)
                      *weights1d(jp)*const_grad_matrix(i, kp);
                    
		    A.AddInteraction(m+iz, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iz, vloc);
		    
                    vloc = -C(0, 1)*mass_matrix(i, ip)*dGL_GL(jp, j)*weights1d(j)
                      *dGL_GL(k, kp)*weights1d(kp);
                    
		    A.AddInteraction(m+iz, n+jz, vloc);
		    A.AddInteraction(m+jz, n+iz, vloc);
                    
                    vloc = -C(0, 1)*GLint_GL(ip, j)*weights1d(j)*GLint_GL(i, jp)
                      *weights1d(jp)*stiff1d(k, kp);
                    
		    A.AddInteraction(m+iy, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iy, vloc);
                    
                    vloc = C(0, 1)*dGL_GL(jp, j)*weights1d(j)*GLint_GL(i, kp)
                      *weights1d(kp)*const_grad_matrix(ip, k);
                    
		    A.AddInteraction(m+iy, n+jz, vloc);
		    A.AddInteraction(m+jz, n+iy, vloc);
                    
                    // c20
                    vloc = C(0, 2)*dGL_GL(jp, j)*weights1d(j)*const_grad_matrix(ip, k)
                      *GLint_GL(i, kp)*weights1d(kp);
                    
		    A.AddInteraction(m+iz, n+jy, vloc);
		    A.AddInteraction(m+jy, n+iz, vloc);
                    
                    vloc = -C(0, 2)*mass_matrix(i, ip)*dGL_GL(jp, j)*weights1d(j)
                      *dGL_GL(k, kp)*weights1d(kp);
                    
		    A.AddInteraction(m+iy, n+jy, vloc);
		    A.AddInteraction(m+jy, n+iy, vloc);
                    
                    vloc = -C(0, 2)*GLint_GL(ip, j)*weights1d(j)*GLint_GL(i, kp)
                      *weights1d(kp)*stiff1d(jp, k);
                    
		    A.AddInteraction(m+iz, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iz, vloc);
                    
                    vloc = C(0, 2)*dGL_GL(k, kp)*weights1d(kp)*GLint_GL(ip, j)
                      *weights1d(j)*const_grad_matrix(i, jp);
                    
		    A.AddInteraction(m+iy, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iy, vloc);
                    
                    // c21
                    vloc = C(1, 2)*dGL_GL(k, kp)*weights1d(kp)*const_grad_matrix(i, jp)
                      *GLint_GL(ip, j)*weights1d(j);
                    
		    A.AddInteraction(m+ix, n+jy, vloc);
		    A.AddInteraction(m+jy, n+ix, vloc);
                    
                    vloc = -C(1, 2)*mass_matrix(i, ip)*dGL_GL(jp, j)*weights1d(j)
                      *dGL_GL(k, kp)*weights1d(kp);
                    
		    A.AddInteraction(m+ix, n+jx, vloc);
		    A.AddInteraction(m+jx, n+ix, vloc);
                    
                    vloc = -C(1, 2)*GLint_GL(ip, k)*weights1d(k)*GLint_GL(i, kp)
                      *weights1d(kp)*stiff1d(jp, j);
                    
		    A.AddInteraction(m+iz, n+jy, vloc);
		    A.AddInteraction(m+jy, n+iz, vloc);
                    
                    vloc = C(1, 2)*dGL_GL(jp, k)*weights1d(k)*GLint_GL(i, kp)
                      *weights1d(kp)*const_grad_matrix(ip, j);
                    
		    A.AddInteraction(m+iz, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iz, vloc);  
                  }
          }    
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point
  */
  void HexahedronHcurlOptimalFirstFamily
  ::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    FillZero(res);
    VectReal_wp phiGx(order), phiGy(order), phiGz(order);
    VectReal_wp phiGLx(order+2), phiGLy(order+2), phiGLz(order+2);
    for (int i = 0; i < order; i++)
      {
	phiGx(i) = lob_int.EvaluatePhi(i, point_loc(0));
	phiGy(i) = lob_int.EvaluatePhi(i, point_loc(1));
	phiGz(i) = lob_int.EvaluatePhi(i, point_loc(2));
      }

    for (int i = 0; i <= order+1; i++)
      {
	phiGLx(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiGLy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiGLz(i) = lob_basis.EvaluatePhi(i, point_loc(2));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  {
	    int nx = NumDofs_X(i, j, k);
	    res(nx)(0) = phiGx(i)*phiGLy(j)*phiGLz(k);
	    
	    int ny = NumDofs_Y(j, i, k);
	    res(ny)(1) = phiGLx(j)*phiGy(i)*phiGLz(k);
	    
	    int nz = NumDofs_Z(k, j, i);
	    res(nz)(2) = phiGLx(k)*phiGLy(j)*phiGz(i);
	  }
  }

  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] curl curl of basis functions on point
  */
  void HexahedronHcurlOptimalFirstFamily
  ::ComputeCurlPhiRef(const R3& point_loc, VectR3& curl) const
  {
    curl.Reallocate(nb_dof_loc);
    FillZero(curl);
    VectReal_wp phiGx(order), phiGy(order), phiGz(order),
      phiGLx(order+2), phiGLy(order+2), phiGLz(order+2);
    VectReal_wp dphiGLx(order+2), dphiGLy(order+2), dphiGLz(order+2);
    for (int i = 0; i < order; i++)
      {
	phiGx(i) = lob_int.EvaluatePhi(i, point_loc(0));
	phiGy(i) = lob_int.EvaluatePhi(i, point_loc(1));
	phiGz(i) = lob_int.EvaluatePhi(i, point_loc(2));
      }

    for (int i = 0; i <= order+1; i++)
      {
	phiGLx(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiGLy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiGLz(i) = lob_basis.EvaluatePhi(i, point_loc(2));

	dphiGLx(i) = lob_basis.EvaluatePhiGrad(i, point_loc(0));
	dphiGLy(i) = lob_basis.EvaluatePhiGrad(i, point_loc(1));
	dphiGLz(i) = lob_basis.EvaluatePhiGrad(i, point_loc(2));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k <= order+1; k++)
	  {
	    int nx = NumDofs_X(i, j, k);
	    curl(nx)(1) = dphiGLz(k)*phiGx(i)*phiGLy(j);
	    curl(nx)(2) = -dphiGLy(j)*phiGx(i)*phiGLz(k);
	    
	    int ny = NumDofs_Y(j, i, k);
	    curl(ny)(0) = -dphiGLz(k)*phiGy(i)*phiGLx(j);
	    curl(ny)(2) = dphiGLx(j)*phiGy(i)*phiGLz(k);

	    int nz = NumDofs_Z(k, j, i);
	    curl(nz)(0) = dphiGLy(j)*phiGz(i)*phiGLx(k);
	    curl(nz)(1) = -dphiGLx(k)*phiGz(i)*phiGLy(j);
	  }
    
  }

  
  //! computation of val_phi \f$ = \varphi_{node} (\xi_{num_point})  \f$
  /*!
    \param[in] num_point quadrature point number
    \param[out] phi values of basis functions
    (xi_i) are Gauss-Lobatto quadrature points
  */  
  void HexahedronHcurlOptimalFirstFamily
  ::GetValuePhiOnQuadraturePoint(int num_point, VectR3& phi) const
  {
    // for example, x-basis functions :
    // then \phi_{node}(\xi_{num_point}) = \phi_{i1}^G(\xi_{k1}^{GL}) \phi_{i2}^{GL}(\xi_{k2}^{GL})
    //                                       \phi_{i3}^{GL} (\xi_{k3}^{GL}) e_x
    //                                    = G_GL(i1,k1) \delta_{i2,k2} \delta_{i3,k3} e_x
    
    phi.Reallocate(nb_dof_loc); FillZero(phi);
    int k1 = CoordinateQuad3D(num_point, 0);
    int k2 = CoordinateQuad3D(num_point, 1);
    int k3 = CoordinateQuad3D(num_point, 2);
    
    for (int i = 0; i < order; i++)
      {
        phi(NumDofs_X(i, k2, k3))(0) = GLint_GL(i, k1);
        phi(NumDofs_Y(k1, i, k3))(1) = GLint_GL(i, k2);
        phi(NumDofs_Z(k1, k2, i))(2) = GLint_GL(i, k3);
      }
    
  }
  
  
  //! computation of curl_phi = \f$ \nabla \times \phi_{node} (\xi_{num\_point}) \f$ 
  /*!
    \param[in] num_point quadrature point number
    \param[out] curl_phi curl of basis functions 
  */ 
  void HexahedronHcurlOptimalFirstFamily
  ::GetCurlPhiOnQuadraturePoint(int num_point, VectR3& curl_phi) const
  {
    curl_phi.Reallocate(nb_dof_loc); FillZero(curl_phi);
    int k1 = CoordinateQuad3D(num_point, 0);
    int k2 = CoordinateQuad3D(num_point, 1);
    int k3 = CoordinateQuad3D(num_point, 2);
        
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
        {
          // for x-basis functions
          // curl(\phi_{i, j, k}) (\xi_{num_point}) =
          // |  0 
          // |  \phi_i^G (\xi_k1^{GL}) \phi_j^{GL} (\xi_k2^{GL}) \phi_k^{'GL} (\xi_k3^{GL})
          // |  -\phi_i^G (\xi_k1^{GL}) \phi_j^{'GL} (\xi_k2^{GL}) \phi_k^{GL} (\xi_k3^{GL})
          //
          // = |  0
          //   |  G_GL(i,k1) \delta_{j,k2} dGL_GL(k,k3)
          //   |  -G_GL(i,k1) \delta_{k,k3} dGL_GL(j,k2)
          
          curl_phi(NumDofs_X(i, k2, j))(1) = dGL_GL(j, k3)*GLint_GL(i, k1);
          curl_phi(NumDofs_X(i, j, k3))(2) = -dGL_GL(j, k2)*GLint_GL(i, k1);
          
          // for y-basis functions
          // curl(\phi_{j,i,k}) (\xi_{num_point}) =
          // |  -\phi_j^{GL} (\xi_k1^{GL}) \phi_i^G (\xi_k2^{GL})  \phi_k^{'GL} (\xi_k3^{GL}) 
          // |  0
          // |  \phi_j^{'GL} (\xi_k1^{GL}) \phi_i^G (\xi_k2^{GL})  \phi_k^{GL} (\xi_k3^{GL})
          //
          // = |  -G_GL(i,k2) \delta_{j,k1} dGL_GL(k,k3)
          //   |  0
          //   |  G_GL(i,k2) \delta_{k,k3} dGL_GL(j,k1)
    
          curl_phi(NumDofs_Y(k1, i, j))(0) = -dGL_GL(j, k3)*GLint_GL(i, k2);
          curl_phi(NumDofs_Y(j, i, k3))(2) = dGL_GL(j, k1)*GLint_GL(i, k2);

          // for z-basis functions
          // curl(\phi_{j,k,i}) (\xi_{num_point}) =
          // |  \phi_j^{GL} (\xi_k1^{GL}) \phi_k^{'GL} (\xi_k2^{GL})  \phi_i^G (\xi_k3^{GL}) 
          // |  -\phi_j^{'GL} (\xi_k1^{GL}) \phi_k^{GL} (\xi_k2^{GL}) \phi_i^G (\xi_k3^{GL})
          // |  0
          //
          // = |  G_GL(i,k3) \delta_{j,k1} dGL_GL(k,k2)
          //   |  -G_GL(i,k3) \delta_{k,k2} dGL_GL(j,k1)
          //   |  0
          
          curl_phi(NumDofs_Z(k1, j, i))(0) = dGL_GL(j, k2)*GLint_GL(i, k3);
          curl_phi(NumDofs_Z(j, k2, i))(1) = -dGL_GL(j, k1)*GLint_GL(i, k3);
        }
    
  }
    

  //! displays details of class HexahedronHcurlOptimalFirstFamily
  ostream& operator <<(ostream& out, const HexahedronHcurlOptimalFirstFamily& e)
  {
    out<<static_cast<const HexahedronReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_OPTIMAL_FIRST_FAMILY_CXX
#endif
