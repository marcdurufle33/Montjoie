#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHcurlFirstFamily::HexahedronHcurlFirstFamily() : HexahedronReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
    optimized_elem_matrix = true;
    this->mass_lumping_ortho = true;
 
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;    
  }
  

  //! returns size of memory used by the object
  size_t HexahedronHcurlFirstFamily::GetMemorySize() const
  {
    size_t taille = HexahedronReference<2>::GetMemorySize();
    taille += gauss.GetMemorySize();
    taille += NumDofs_X.GetMemorySize() + NumDofs_Y.GetMemorySize()
                           +NumDofs_Z.GetMemorySize()+CoordinateDofs.GetMemorySize();
    
    taille += G_GL.GetMemorySize()+GL_G.GetMemorySize()+dGL_G.GetMemorySize()
      +dGL_GL.GetMemorySize()+stiff1d.GetMemorySize();
    
    taille += rh_loc.GetMemorySize()+ch_loc.GetMemorySize()+ch1_loc.GetMemorySize();
    taille += ch2_loc.GetMemorySize()+ch3_loc.GetMemorySize();
    taille += rh_exact.GetMemorySize()+ch1_node.GetMemorySize()+ch2_node.GetMemorySize();
    taille += ch3_node.GetMemorySize() + gauss_exact.GetMemorySize();
    taille += G_GX.GetMemorySize()+GL_GX.GetMemorySize()+dGX_GX.GetMemorySize()
      +weights3d_Exact.GetMemorySize()+this->weights_dof.GetMemorySize()+invWeightsMassG.GetMemorySize()
      +invSqrtWeightsMassG.GetMemorySize();
    return taille;
  }
  
  
  //! in order to number a mesh
  void HexahedronHcurlFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return HexahedronReference<2>::ConstructNumberMap(nmap, dg);
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofQuadrangle(this->order, 2*this->order*(this->order+1));
	nmap.SetNbDofHexahedron(this->order, 0);
	
	// rotation of dofs on faces
	element_quad_surf->FindHcurlRotationQuad(nmap, 0);
	
	return;
      }
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    // 2 r (r-1) dofs on each face
    int nb_dof_quad_inside = 2*order*(order+1) - 4*order;
    nmap.SetNbDofQuadrangle(order, nb_dof_quad_inside);
    
    // 3 r (r-1) (r-1) dofs inside the hexahedron
    nmap.SetNbDofHexahedron(order, nb_dof_loc - 12*order
			   - 6*nb_dof_quad_inside);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
        
    // FindHcurlLinearCombinationRotation(*this, nmap, Fb_quad, Fb_quad);
    
    element_quad_surf->FindHcurlRotationQuad(nmap);
  }

  
  //! construction of finite element
  void HexahedronHcurlFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                         int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    HexahedronReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad); 
    
    lob_basis.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    gauss.ConstructQuadrature(order-1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    gauss_exact.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    ConstructFunctions();
    
    ConstructStiffnessMatrix();
    ConstructMassMatrix();
    
    FindDofsOnFace(true);
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    this->nb_points_dof_inside = this->nb_dof_loc;
    this->num_dof_points_surf.Reallocate(6);
    for (int n = 0; n < 6; n++)
      {
	this->num_dof_points_surf(n).Reallocate(nb_dof_quad);
	for (int i = 0; i < nb_dof_quad; i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    
  }

  
  //! construction of basis functions
  void HexahedronHcurlFirstFamily::ConstructFunctions()
  {
    nb_dof_tri = 0;
    nb_dof_quad = 2*order*(order+1);
    nb_dof_loc = 3*order*(order+1)*(order+1);
    
    int nb_dof_inside_quad = order*(order-1)*2;
    nb_dof_boundaries = 6*nb_dof_inside_quad + 12*order;

    QuadrangleHcurlFirstFamily* Fb_quad_ptr = new QuadrangleHcurlFirstFamily();
    QuadrangleHcurlFirstFamily& Fb_quad = *Fb_quad_ptr;
    element_quad_surf = Fb_quad_ptr;
    
    Fb_quad.ConstructFiniteElement(order);
    
    NumDofs_X.Reallocate(order, order+1, order+1);
    NumDofs_Y.Reallocate(order+1, order, order+1);
    NumDofs_Z.Reallocate(order+1, order+1, order);
    
    nb_dof_loc = 0;
    // nodes on edges
    const Matrix<int>& CoordinateNodes = this->Fb_geom.GetCoordinateNodes3D();
    int order_geom = this->GetGeometryOrder();
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
          i = order;
        
        if (j == order_geom)
          j = order;
        
        if (k == order_geom)
          k = order;

        if (j2 == order_geom)
          j2 = order;
        
        if (k2 == order_geom)
          k2 = order;
        
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
    
    nb_dof_loc = 3*order*(order+1)*(order+1);
    
    // dofs on first face
    int offset = 8*order;
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_Y(0, i, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(0, j, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on second face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, 0, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(j, 0, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on third face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, j, 0) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Y(j, i, 0) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on fourth face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, j, order) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Y(j, i, order) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on fifth face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_X(i, order, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(j, order, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}
    
    // dofs on sixth face
    offset += Fb_quad.GetNbDof() - 4*order;
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i < order ; i++)
	{
	  NumDofs_Y(order, i, j) = offset + Fb_quad.GetXdofNumber(i, j);
	  NumDofs_Z(order, j, i) = offset + Fb_quad.GetYdofNumber(j, i);
	}    
    
    // dofs inside the hexahedron
    offset += Fb_quad.GetNbDof();
    for (int i = 0 ; i < order ; i++)
      for (int j = 1 ; j < order ; j++)
	for (int k = 1 ; k < order ; k++)
	  {
	    NumDofs_X(i, j, k) = offset++;
	    NumDofs_Y(j, i, k) = offset++;
	    NumDofs_Z(j, k, i) = offset++;
	  }
    
    VectR2 points_dof2d = Fb_quad.PointsDofND();
    this->SetPointsDof2D_quad(points_dof2d);

    VectR3 points_dof3d;
    points_dof3d.Reallocate(nb_dof_loc);
    CoordinateDofs.Reallocate(nb_dof_loc,4);
    
    // DISP(NumDofs_X); DISP(NumDofs_Y); DISP(NumDofs_Z);
    for (int i = 0 ; i < order ; i++)
      for (int j = 0 ; j <= order ; j++)
	for (int k = 0 ; k <= order ; k++)
	  {
	    points_dof3d(NumDofs_X(i,j,k))
              .Init(gauss.Points(i),lob_basis.Points(j),lob_basis.Points(k));
	    
            points_dof3d(NumDofs_Y(j,i,k))
              .Init(lob_basis.Points(j),gauss.Points(i),lob_basis.Points(k));
	    
            points_dof3d(NumDofs_Z(k,j,i))
              .Init(lob_basis.Points(k),lob_basis.Points(j),gauss.Points(i));
	    
	    CoordinateDofs(NumDofs_X(i,j,k),0) = i;
	    CoordinateDofs(NumDofs_X(i,j,k),1) = j;
	    CoordinateDofs(NumDofs_X(i,j,k),2) = k;
	    CoordinateDofs(NumDofs_X(i,j,k),3) = 0;
	    CoordinateDofs(NumDofs_Y(j,i,k),0) = j;
	    CoordinateDofs(NumDofs_Y(j,i,k),1) = i;
	    CoordinateDofs(NumDofs_Y(j,i,k),2) = k;
	    CoordinateDofs(NumDofs_Y(j,i,k),3) = 1;
	    CoordinateDofs(NumDofs_Z(k,j,i),0) = k;
	    CoordinateDofs(NumDofs_Z(k,j,i),1) = j;
	    CoordinateDofs(NumDofs_Z(k,j,i),2) = i;
	    CoordinateDofs(NumDofs_Z(k,j,i),3) = 2;
	  }
    
    this->SetPointsDofND(points_dof3d);

    
  }
  
  
  //! computation of some arrays, G_GL, GL_G and dGL_G
  //! needed for the computation of stiffness matrix
  void HexahedronHcurlFirstFamily::ConstructStiffnessMatrix()
  {
    G_GL.Reallocate(order,order+1); GL_G.Reallocate(order+1,order);
    dGL_G.Reallocate(order+1,order);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  G_GL(i, j) = gauss.EvaluatePhi(i, lob_basis.Points(j));
	  GL_G(j, i) = lob_basis.EvaluatePhi(j, gauss.Points(i));
	  dGL_G(j, i) = lob_basis.EvaluatePhiGrad(j, gauss.Points(i));
	}
    
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    dGL_GL = lob_basis.GradPhi();
    
    stiff1d.Reallocate(order+1, order+1);
    stiff1d.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          stiff1d(i, j) += this->Weights1D(k)*dGL_GL(i, k)*dGL_GL(j, k);

    Matrix<Real_wp, General, ArrayRowSparse>
      Rh(3*nb_points_quadrature_inside, 3*nb_points_quadrature_inside);
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m <= order; m++)
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
    
    G_GX.Reallocate(order, order+1); GL_GX.Reallocate(order+1, order+1);
    dGX_GX.Reallocate(order+1, order+1);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	G_GX(i,j) = gauss.EvaluatePhi(i, gauss_exact.Points(j));
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  GL_GX(i,j) = lob_basis.EvaluatePhi(i, gauss_exact.Points(j));
	  dGX_GX(i,j) = gauss_exact.EvaluatePhiGrad(i, gauss_exact.Points(j));
	}
    
    int N = nb_points_quadrature_inside;
    weights3d_Exact.Reallocate(N);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  weights3d_Exact(NumQuad3D(i,j,k)) = gauss_exact.Weights(i)
	    *gauss_exact.Weights(j)*gauss_exact.Weights(k);
    
    // for exact integration
    Matrix<Real_wp, General, ArrayRowSparse> RhG(3*N, 3*N);
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m <= order; m++)
	    {
	      int col = 3*NumQuad3D(i,j,k);
	      int row = 3*NumQuad3D(i,j,m)+1;
	      RhG.AddInteraction(row, col, -dGX_GX(m,k));
	      	      
	      row = 3*NumQuad3D(i,m,k)+2;
	      RhG.AddInteraction(row, col, dGX_GX(m,j));
	      	      
	      col = 3*NumQuad3D(i,j,k)+1;
	      row = 3*NumQuad3D(i,j,m);
	      RhG.AddInteraction(row, col, dGX_GX(m,k));
	      	      
	      row = 3*NumQuad3D(m,j,k)+2;
	      RhG.AddInteraction(row, col, -dGX_GX(m,i));
	      	      
	      col = 3*NumQuad3D(i,j,k)+2;
	      row = 3*NumQuad3D(i,m,k);
	      RhG.AddInteraction(row, col, -dGX_GX(m,j));
	      	      
	      row = 3*NumQuad3D(m,j,k)+1;
	      RhG.AddInteraction(row, col, dGX_GX(m,i));
            }
    
    Seldon::Copy(RhG, rh_exact);
  }
  
  
  //! mass matrix
  void HexahedronHcurlFirstFamily::ConstructMassMatrix()
  {
    int N = nb_points_quadrature_inside;
    // constructing operator for interpolation from dof points to Gauss-Lobatto points
    Matrix<Real_wp, General, ArrayRowSparse> Ch(nb_dof_loc, 3*N);
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m < order; m++)
	    {
	      int node = NumQuad3D(i,j,k);
	      int num_dof = NumDofs_X(m,j,k);
	      Ch.AddInteraction(num_dof, 3*node, G_GL(m,i));
	      	      
	      num_dof = NumDofs_Y(i,m,k);
	      Ch.AddInteraction(num_dof, 3*node+1, G_GL(m,j));
	      	      
	      num_dof = NumDofs_Z(i,j,m);
	      Ch.AddInteraction(num_dof, 3*node+2, G_GL(m,k));
            }

    Seldon::Copy(Ch, ch_loc);
    
    // and operator if exact integration is required
    Matrix<Real_wp, General, ArrayRowSparse> Ch1(nb_dof_loc, 3*N);
    Matrix<Real_wp, General, ArrayRowSparse> Ch2(3*N, 3*N);
    Matrix<Real_wp, General, ArrayRowSparse> Ch3(3*N, 3*N);
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m < order; m++)
	    {
	      int node = NumQuad3D(i,j,k);
	      int num_dof = NumDofs_X(m,j,k);
	      Ch1.AddInteraction(num_dof, 3*node, G_GX(m,i));
	      	      
	      num_dof = NumDofs_Y(i,m,k);
	      Ch1.AddInteraction(num_dof, 3*node+1, G_GX(m,j));
	      	      
	      num_dof = NumDofs_Z(i,j,m);
	      Ch1.AddInteraction(num_dof, 3*node+2, G_GX(m,k));
            }
    
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m <= order; m++)
	    {
	      int node1 = NumQuad3D(i,j,k);
	      int node2 = NumQuad3D(i,m,k);
	      Ch2.AddInteraction(3*node2, 3*node1, GL_GX(m,j));
	      	      
	      node2 = NumQuad3D(m,j,k);
	      Ch2.AddInteraction(3*node2+1, 3*node1+1, GL_GX(m,i));
	      	      
	      node2 = NumQuad3D(m,j,k);
	      Ch2.AddInteraction(3*node2+2, 3*node1+2, GL_GX(m,i));
	      	      
	      node2 = NumQuad3D(i,j,m);
	      Ch3.AddInteraction(3*node2, 3*node1, GL_GX(m,k));
	      	      
	      node2 = NumQuad3D(i,j,m);
	      Ch3.AddInteraction(3*node2+1, 3*node1+1, GL_GX(m,k));
	      	      
	      node2 = NumQuad3D(i,m,k);
	      Ch3.AddInteraction(3*node2+2, 3*node1+2, GL_GX(m,j));
            }
    
    Seldon::Copy(Ch1, ch1_loc);
    Seldon::Copy(Ch2, ch2_loc);
    Seldon::Copy(Ch3, ch3_loc);
    
    this->weights_dof.Reallocate(nb_dof_loc);
    invWeightsMassG.Reallocate(nb_dof_loc);
    invSqrtWeightsMassG.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int num_dof = NumDofs_X(i, j, k);
            this->weights_dof(num_dof) = gauss.Weights(i)*lob_basis.Weights(j)*lob_basis.Weights(k);
            invWeightsMassG(num_dof) = 1.0/this->weights_dof(num_dof);
            invSqrtWeightsMassG(num_dof) = 1.0/sqrt(this->weights_dof(num_dof));
            
            num_dof = NumDofs_Y(j, i, k);
            this->weights_dof(num_dof) = gauss.Weights(i)*lob_basis.Weights(j)*lob_basis.Weights(k);
            invWeightsMassG(num_dof) = 1.0/this->weights_dof(num_dof);
            invSqrtWeightsMassG(num_dof) = 1.0/sqrt(this->weights_dof(num_dof));
            
            num_dof = NumDofs_Z(k, j, i);
            this->weights_dof(num_dof) = gauss.Weights(i)*lob_basis.Weights(j)*lob_basis.Weights(k);
            invWeightsMassG(num_dof) = 1.0/this->weights_dof(num_dof);
            invSqrtWeightsMassG(num_dof) = 1.0/sqrt(this->weights_dof(num_dof));
          }
    
    int order_geom = this->GetGeometryOrder();
    const VectReal_wp& points_nodal1d = this->PointsNodal1D();
    Matrix<Real_wp> G_Geom(order, order+1), GL_Geom(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        GL_Geom(i, j) = lob_basis.EvaluatePhi(i, points_nodal1d(j));
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order_geom; j++)
        G_Geom(i, j) = gauss.EvaluatePhi(i, points_nodal1d(j));
    
    Ch1.Clear(); Ch2.Clear(); Ch3.Clear();
    Ch1.Reallocate(3*order*(order_geom+1)*(order_geom+1), 3*this->GetNbPointsNodalElt());
    Ch2.Reallocate(3*order*(order+1)*(order_geom+1), 3*order*(order_geom+1)*(order_geom+1));
    Ch3.Reallocate(3*order*(order+1)*(order+1), 3*order*(order+1)*(order_geom+1));
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
	  for (int m = 0; m <= order; m++)
	    {
              int icol = (order_geom+1)*(k*(order_geom+1)+i) + j;
              int irow = (order_geom+1)*(k*(order+1) + m) + j;
	      Ch2.AddInteraction(3*irow, 3*icol, GL_Geom(m, i));
	      	      
              icol = (order_geom+1)*(i*order + k) + j;
              irow = (order_geom+1)*(m*order + k) + j;
	      Ch2.AddInteraction(3*irow+1, 3*icol+1, GL_Geom(m, i));
	      	      
	      icol = order*(i*(order_geom+1) + j) + k;
              irow = order*(m*(order_geom+1) + j) + k;
	      Ch2.AddInteraction(3*irow+2, 3*icol+2, GL_Geom(m, i));
            }

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k < order; k++)
	  for (int m = 0; m <= order; m++)
	    {    
              int icol = (order_geom+1)*(k*(order+1) + i) + j;
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
  void HexahedronHcurlFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    contrib.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      {
	contrib(i) = feval(3*i + CoordinateDofs(i,3));
      }
  }
  
  
  //! projection from values on quadrature points to dof components
  /*!
    \param[in] Equad values on quadrature points
    \param[out] Edof dof components
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlFirstFamily::
  ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const
  {
    Edof.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    // basic interpolation
	    int num_dof = NumDofs_X(i,j,k);
	    Edof(num_dof) = 0;
            for (int l = 0; l <= order; l++)
	      Edof(num_dof) += GL_G(l, i)*Equad(3*NumQuad3D(l,j,k));
	    
	    num_dof = NumDofs_Y(j,i,k);
	    Edof(num_dof) = 0;
	    for (int l = 0; l <= order; l++)
	      Edof(num_dof) += GL_G(l, i)*Equad(3*NumQuad3D(j,l,k)+1);
	    
	    num_dof = NumDofs_Z(k,j,i);
	    Edof(num_dof) = 0;
	    for (int l = 0; l<= order; l++)
	      Edof(num_dof) += GL_G(l, i)*Equad(3*NumQuad3D(k,j,l)+2);
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
  void HexahedronHcurlFirstFamily::
  ApplyShGen(const T0& alpha, int num_loc, const Vector1& feval, Vector2& res, int r) const
  {
    if ((r != 0) && (r != order))
      {
        abort();
      }
    
    //typedef typename Vector2::value_type T2;
    int ext = 0, num_dof;
    if (num_loc > 2)
      ext = order;

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
	  for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                res(NumDofs_X(i, j, k)) += alpha*G_GL(i, ext)*feval(3*NumQuad2D(j, k));
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Y(j, i, k);
                    for (int m = 0; m <= order; m++)
                      res(num_dof) += alpha*G_GL(i, m)*feval(3*NumQuad2D(m, k)+1);
                  }

                if (k == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order; n++)
                      res(num_dof) += alpha*G_GL(i, n)*feval(3*NumQuad2D(j, n)+2);
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
	  for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                if (j == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order; m++)
                      res(num_dof) += alpha*G_GL(i, m)*feval(3*NumQuad2D(m, k));
                  }
                
                num_dof = NumDofs_Y(j, i, k);
                res(num_dof) += alpha*G_GL(i, ext)*feval(3*NumQuad2D(j, k)+1);
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order; n++)
                      res(num_dof) += alpha*G_GL(i, n)*feval(3*NumQuad2D(k, n)+2);
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
	  for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                if (k == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order; m++)
                      res(num_dof) += alpha*G_GL(i, m)*feval(3*NumQuad2D(m, j));
                    
                    num_dof = NumDofs_Y(j, i, k);
                    for (int n = 0; n <= order; n++)
                      res(num_dof) += alpha*G_GL(i, n)*feval(3*NumQuad2D(j, n)+1);
                  }
                
                num_dof = NumDofs_Z(k, j, i);
                res(num_dof) += alpha*G_GL(i, ext)*feval(3*NumQuad2D(k, j)+2);
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
  void HexahedronHcurlFirstFamily::
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
      ext = order;

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                feval(3*NumQuad2D(j, k)) += G_GL(i, ext)*res(NumDofs_X(i, j, k));
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Y(j, i, k);
                    for (int m = 0; m <= order; m++)
                      feval(3*NumQuad2D(m, k)+1) += G_GL(i, m)*res(num_dof);
                  }

                if (k == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order; n++)
                      feval(3*NumQuad2D(j, n)+2) += G_GL(i, n)*res(num_dof);
                  }
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                if (j == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order; m++)
                      feval(3*NumQuad2D(m, k)) += G_GL(i, m)*res(num_dof);
                  }
                
                num_dof = NumDofs_Y(j, i, k);
                feval(3*NumQuad2D(j, k)+1) += G_GL(i, ext)*res(num_dof);
                
                if (j == ext)
                  {
                    num_dof = NumDofs_Z(k, j, i);
                    for (int n = 0; n <= order; n++)
                      feval(3*NumQuad2D(k, n)+2) += G_GL(i, n)*res(num_dof);
                  }
              }
      }
    else
      {
	for (int i = 0; i < order; i++)
	  for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                if (k == ext)
                  {
                    num_dof = NumDofs_X(i, j, k);
                    for (int m = 0; m <= order; m++)
                      feval(3*NumQuad2D(m, j)) += G_GL(i, m)*res(num_dof);
                    
                    num_dof = NumDofs_Y(j, i, k);
                    for (int n = 0; n <= order; n++)
                      feval(3*NumQuad2D(j, n)+1) += G_GL(i, n)*res(num_dof);
                  }
                
                num_dof = NumDofs_Z(k, j, i);
                feval(3*NumQuad2D(k, j)+2) += G_GL(i, ext)*res(num_dof);
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
  void HexahedronHcurlFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    int Nquad = this->Points2D_quad().GetM();
    Vector1 fx(3*Nquad), fy(3*Nquad);
    FillZero(fx); FillZero(fy);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    // we express feval on Gauss-Lobatto points
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          int node = 3*NumQuad2D(i, j);
          for (int m = 0; m <= order; m++)
            {
              int p = 3*NumQuad2D(m, j);
              fx(node) += GL_GX(i, m)*feval(p);
              fx(node+1) += GL_GX(i, m)*feval(p+1);
              fx(node+2) += GL_GX(i, m)*feval(p+2);
            }
        }

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          int node = 3*NumQuad2D(i, j);
          for (int m = 0; m <= order; m++)
            {
              int p = 3*NumQuad2D(i, m);
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
  void HexahedronHcurlFirstFamily::
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
      ext = order;
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	// Face x = 0 or x = 1	
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                // contribution for curl phi_i,j,k^x
                vloc = T2(0);
                for (int n = 0; n <= order; n++)
                  vloc += dGL_GL(k, n)*feval(3*NumQuad2D(j, n)+1);
		                
                for (int m = 0; m <= order; m++)
                  vloc -= dGL_GL(j, m)*feval(3*NumQuad2D(m, k)+2);
                
                res(NumDofs_X(i, j, k) ) += alpha*G_GL(i, ext)*vloc;
               
                // contribution for curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k); 
                if (j == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        vloc -= G_GL(i, m)*dGL_GL(k, n)*feval(3*NumQuad2D(m, n));
                    
                    res(num_dof) += alpha*vloc;
                  }
                    
                vloc = T2(0);
                for (int m = 0; m <= order; m++)
                  vloc += G_GL(i, m)*feval(3*NumQuad2D(m, k)+2);
		                
                res(num_dof) += alpha*dGL_GL(j, ext)*vloc;
                
                // contribution for curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (k == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        vloc += G_GL(i, n)*dGL_GL(j, m)*feval(3*NumQuad2D(m, n));
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int n = 0; n <= order; n++)
                  vloc += G_GL(i, n)*feval(3*NumQuad2D(j, n)+1);
                
                res(num_dof) -= alpha*dGL_GL(k, ext)*vloc;
                
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
        // Face y = 0 or y = 1
         for (int i = 0; i < order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (j == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        vloc += G_GL(i, m)*dGL_GL(k, n)*feval(3*NumQuad2D(m, n)+1);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int m = 0; m <= order; m++)
                  vloc += G_GL(i, m)*feval(3*NumQuad2D(m, k)+2);
                
                res(num_dof) -= alpha*dGL_GL(j, ext)*vloc;
                              
                // contribution for curl phi_j,i,k^y 
                vloc = T2(0) ;
                for (int n = 0; n <= order; n++)
                  vloc -= dGL_GL(k, n)*feval(3*NumQuad2D(j, n));
                
                for (int m = 0; m <= order; m++)
                  vloc += dGL_GL(j, m)*feval(3*NumQuad2D(m, k)+2);
                
                res(NumDofs_Y(j, i, k)) += alpha*G_GL(i, ext)*vloc;
                
                //contribution for  curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (j == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        vloc -= G_GL(i, n)*dGL_GL(k, m)*feval(3*NumQuad2D(m, n)+1);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int n = 0; n <= order; n++)
                  vloc += G_GL(i, n)*feval(3*NumQuad2D(k, n));
                
                res(num_dof) += alpha*dGL_GL(j, ext)*vloc;
              }
      }
    else if ((num_loc == 2)||(num_loc == 3))
      {
        // Face z = 0 or z = 1
         
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {

                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (k == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        vloc -= G_GL(i, m)*dGL_GL(j, n)*feval(3*NumQuad2D(m, n)+2);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int m = 0; m <= order; m++)
                  vloc += G_GL(i, m)*feval(3*NumQuad2D(m, j)+1);
                
                res(num_dof) += alpha*dGL_GL(k, ext)*vloc;
                
                // contribution for  curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k);
                if (k == ext)
                  {
                    vloc = T2(0);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        vloc += G_GL(i, n)*dGL_GL(j, m)*feval(3*NumQuad2D(m, n)+2);
                    
                    res(num_dof) += alpha*vloc;
                  }
                
                vloc = T2(0);
                for (int n = 0; n <= order; n++)
                  vloc += G_GL(i, n)*feval(3*NumQuad2D(j, n));
                
                res(num_dof) -= alpha*dGL_GL(k, ext)*vloc;
                
                // contribution for curl phi_k,j,i^z 
                vloc = T2(0) ;
                for (int n = 0; n <= order; n++)
                  vloc += dGL_GL(j, n)*feval(3*NumQuad2D(k, n));
                
                for (int m = 0; m <= order; m++)
                  vloc -= dGL_GL(k, m)*feval(3*NumQuad2D(m, j)+1);
                
                res(NumDofs_Z(k, j, i)) += alpha*G_GL(i, ext)*vloc;
                
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
  void HexahedronHcurlFirstFamily::
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
      ext = order;
    
    feval.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc == 5))
      {
	// Face x = 0 or x = 1	
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                // contribution for curl phi_i,j,k^x
                vloc = G_GL(i, ext)*res(NumDofs_X(i, j, k) );

                for (int n = 0; n <= order; n++)
                  feval(3*NumQuad2D(j, n)+1) += dGL_GL(k, n)*vloc;
		                
                for (int m = 0; m <= order; m++)
                  feval(3*NumQuad2D(m, k)+2) -= dGL_GL(j, m)*vloc;
                               
                // contribution for curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k); 
                if (j == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        feval(3*NumQuad2D(m, n)) -= G_GL(i, m)*dGL_GL(k, n)*vloc;
                  }
                    
                vloc = dGL_GL(j, ext)*res(num_dof);
                for (int m = 0; m <= order; m++)
                  feval(3*NumQuad2D(m, k)+2) += G_GL(i, m)*vloc;
                
                // contribution for curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (k == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        feval(3*NumQuad2D(m, n)) += G_GL(i, n)*dGL_GL(j, m)*vloc;
                    
                  }
                
                vloc = -dGL_GL(k, ext)*res(num_dof);
                for (int n = 0; n <= order; n++)
                  feval(3*NumQuad2D(j, n)+1) += G_GL(i, n)*vloc;
                
              }
      }
    else if ((num_loc == 1)||(num_loc == 4))
      {
        // Face y = 0 or y = 1
         for (int i = 0; i < order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (j == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        feval(3*NumQuad2D(m, n)+1) += G_GL(i, m)*dGL_GL(k, n)*vloc;
                  }
                
                vloc = -dGL_GL(j, ext)*res(num_dof);
                for (int m = 0; m <= order; m++)
                  feval(3*NumQuad2D(m, k)+2) += G_GL(i, m)*vloc;
                              
                // contribution for curl phi_j,i,k^y 
                vloc = G_GL(i, ext)*res(NumDofs_Y(j, i, k));
                for (int n = 0; n <= order; n++)
                  feval(3*NumQuad2D(j, n)) -= dGL_GL(k, n)*vloc;
                
                for (int m = 0; m <= order; m++)
                  feval(3*NumQuad2D(m, k)+2) += dGL_GL(j, m)*vloc;
                                
                // contribution for  curl phi_k,j,i^z
                num_dof = NumDofs_Z(k, j, i);
                if (j == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        feval(3*NumQuad2D(m, n)+1) -= G_GL(i, n)*dGL_GL(k, m)*vloc;
                  }
                
                vloc = dGL_GL(j, ext)*res(num_dof);
                for (int n = 0; n <= order; n++)
                  feval(3*NumQuad2D(k, n)) += G_GL(i, n)*vloc;
              }
      }
    else if ((num_loc == 2)||(num_loc == 3))
      {
        // Face z = 0 or z = 1
         
        for (int i = 0; i < order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {

                // contribution for curl phi_i,j,k^x
                num_dof = NumDofs_X(i, j, k);
                if (k == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        feval(3*NumQuad2D(m, n)+2) -= G_GL(i, m)*dGL_GL(j, n)*vloc;
                  }
                
                vloc = dGL_GL(k, ext)*res(num_dof);
                for (int m = 0; m <= order; m++)
                  feval(3*NumQuad2D(m, j)+1) += G_GL(i, m)*vloc;
                
                // contribution for  curl phi_j,i,k^y
                num_dof = NumDofs_Y(j, i, k);
                if (k == ext)
                  {
                    vloc = res(num_dof);
                    for (int m = 0; m <= order; m++)
                      for (int n = 0; n <= order; n++)
                        feval(3*NumQuad2D(m, n)+2) += G_GL(i, n)*dGL_GL(j, m)*vloc;
                    
                  }
                
                vloc = -dGL_GL(k, ext)*res(num_dof);
                for (int n = 0; n <= order; n++)
                  feval(3*NumQuad2D(j, n)) += G_GL(i, n)*vloc;
                
                // contribution for curl phi_k,j,i^z 
                vloc = G_GL(i, ext)*res(NumDofs_Z(k, j, i));
                for (int n = 0; n <= order; n++)
                  feval(3*NumQuad2D(k, n))+= dGL_GL(j, n)*vloc;
                
                for (int m = 0; m <= order; m++)
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
  void HexahedronHcurlFirstFamily::
  ComputeNodalValuesGen(const Vector1& Uloc, Vector2& Uloc_node) const
  {
    if (order == this->GetGeometryOrder())
      {
        Uloc_node.Fill(0);
	const Array3D<int>& NumNodes3D = this->Fb_geom.GetNumNodes3D();
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int node = NumNodes3D(i, j, k);
                for (int m = 0; m < order; m++)
                  {
                    Uloc_node(3*node) += G_GL(m,i)*Uloc(NumDofs_X(m,j,k));
                    Uloc_node(3*node+1) += G_GL(m,j)*Uloc(NumDofs_Y(i,m,k));
                    Uloc_node(3*node+2) += G_GL(m,k)*Uloc(NumDofs_Z(i,j,m));
                  }
              }
      }
    else
      {
        Vector2 Uy(ch3_node.GetN()), Ux(ch2_node.GetN());
        Mlt(SeldonTrans, ch3_node, Uloc, Uy);
        Mlt(SeldonTrans, ch2_node, Uy, Ux);
        Mlt(SeldonTrans, ch1_node, Ux, Uloc_node);
      }
  }
  
    
  //! computation of u on nodal points of a face
  /*!
    \param[in] u_loc dof components of u
    \param[out] u_boundary values of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlFirstFamily::
  ComputeValueBoundaryGen(const Vector1& u_loc, Vector2& u_boundary, int num_loc) const
  {
    typedef typename Vector1::value_type Complexe;
    int Nquad = this->PointsNodal2D_quad().GetM();
    Complexe Ex, Ey, Ez;
    u_boundary.Reallocate(3*Nquad);
    int ext = 0;
    if (num_loc > 2)
      ext = order;
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((num_loc == 0)||(num_loc==5))
      {
	// face x = 0 or x = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      // interpolation
	      Ey = 0; Ez = 0; Ex = 0;
	      for (int k = 0; k < order; k++)
		{
		  Ex += G_GL(k,ext)*u_loc(NumDofs_X(k,i,j));
		  Ey += G_GL(k,i)*u_loc(NumDofs_Y(ext,k,j));
		  Ez += G_GL(k,j)*u_loc(NumDofs_Z(ext,i,k));
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
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      // dofs along y and z
	      Ey = 0; Ez = 0; Ex = 0;
	      for (int k = 0; k < order; k++)
		{
		  Ex += G_GL(k, i)*u_loc(NumDofs_X(k,ext,j));
		  Ey += G_GL(k, ext)*u_loc(NumDofs_Y(i, k, j));
		  Ez += G_GL(k, j)*u_loc(NumDofs_Z(i,ext,k));
		}
	      
	      u_boundary(3*NumQuad2D(i,j)) = Ex;
	      u_boundary(3*NumQuad2D(i,j)+1) = Ey;
	      u_boundary(3*NumQuad2D(i,j)+2) = Ez;
	    }
      }
    else if ((num_loc == 2)||(num_loc==3))
      {
	// face z = 0 or z = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      // dofs along x, y and z
	      Ey = 0; Ez = 0; Ex = 0;
	      for (int k = 0; k < order; k++)
		{
		  Ex += G_GL(k,i)*u_loc(NumDofs_X(k,j,ext));
		  Ey += G_GL(k,j)*u_loc(NumDofs_Y(i,k,ext));
		  Ez += G_GL(k,ext)*u_loc(NumDofs_Z(i,j,k));
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
  void HexahedronHcurlFirstFamily::
  ComputeCurlBoundaryGen(const Vector1& u_loc,
			 Vector2& u_boundary, int num_loc) const
  {
    typedef typename Vector1::value_type Complexe;
    Complexe Ex, Ey, Ez;
    int Nquad = this->PointsNodal2D_quad().GetM();
    u_boundary.Reallocate(3*Nquad);
    int ext = 0;
    if (num_loc > 2)
      ext = order;

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
    if ((num_loc == 0)||(num_loc==5))
      {
	// face x = 0 or x = 1
	Matrix<Complexe> Ex_GL(order+1,order+1), Ey_GL(order+1,order+1),
	  Ez_GL(order+1,order+1);
	Matrix<Complexe> dEy_G(order,order+1), dEz_G(order,order+1);
	Ex_GL.Fill(0); Ey_GL.Fill(0); Ez_GL.Fill(0); dEy_G.Fill(0); dEz_G.Fill(0);
	// first step, interpolation at GL points
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      for (int k = 0; k < order; k++)
		{
		  Ex_GL(i,j) += G_GL(k,ext)*u_loc(NumDofs_X(k,i,j));
		  Ey_GL(i,j) += G_GL(k,i)*u_loc(NumDofs_Y(ext,k,j));
		  Ez_GL(i,j) += G_GL(k,j)*u_loc(NumDofs_Z(ext,i, k));
		  dEy_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Y(j,k,i));
		  dEz_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Z(j,i,k));
		}
	    }
	
	// second step : derivative
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      Ey = 0; Ex = 0; Ez = 0;
	      for (int k = 0; k <= order; k++)
		{
		  Ey += dGL_GL(k,j)*Ex_GL(i,k);
		  Ez -= dGL_GL(k,i)*Ex_GL(k,j);
		  Ex -= dGL_GL(k,j)*Ey_GL(i,k);
		  Ex += dGL_GL(k,i)*Ez_GL(k,j);
		}
	      
	      for (int k = 0; k < order; k++)
		{
		  Ez += G_GL(k,i)*dEy_G(k,j);
		  Ey -= G_GL(k,j)*dEz_G(k,i);
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
	Matrix<Complexe> Ex_GL(order+1,order+1), Ey_GL(order+1,order+1),
	  Ez_GL(order+1,order+1);
	Matrix<Complexe> dEx_G(order,order+1), dEz_G(order,order+1);
	Ex_GL.Fill(0); Ey_GL.Fill(0); Ez_GL.Fill(0); dEx_G.Fill(0); dEz_G.Fill(0);
	// first step, interpolation at GL points
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      for (int k = 0; k < order; k++)
		{
		  Ex_GL(i,j) += G_GL(k,i)*u_loc(NumDofs_X(k,ext,j));
		  Ey_GL(i,j) += G_GL(k,ext)*u_loc(NumDofs_Y(i,k,j));
		  Ez_GL(i,j) += G_GL(k,j)*u_loc(NumDofs_Z(i,ext,k));
		  dEx_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_X(k,j,i));
		  dEz_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Z(i,j,k));
		}
	    }
	
	// second step : derivative
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      Ey = 0; Ex = 0; Ez = 0;
	      for (int k = 0; k <= order; k++)
		{
		  Ey += dGL_GL(k,j)*Ex_GL(i,k);
		  Ex -= dGL_GL(k,j)*Ey_GL(i,k);
		  Ez += dGL_GL(k,i)*Ey_GL(k,j);
		  Ey -= dGL_GL(k,i)*Ez_GL(k,j);
		}
	      
	      for (int k = 0; k < order; k++)
		{
		  Ez -= G_GL(k,i)*dEx_G(k,j);
		  Ex += G_GL(k,j)*dEz_G(k,i);
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
	Matrix<Complexe> Ex_GL(order+1,order+1), Ey_GL(order+1,order+1),
	  Ez_GL(order+1,order+1);
	Matrix<Complexe> dEx_G(order,order+1), dEy_G(order,order+1);
	Ex_GL.Fill(0); Ey_GL.Fill(0); Ez_GL.Fill(0); dEx_G.Fill(0); dEy_G.Fill(0);
	// first step, interpolation at GL points
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      for (int k = 0; k < order; k++)
		{
		  Ex_GL(i,j) += G_GL(k,i)*u_loc(NumDofs_X(k,j,ext));
		  Ey_GL(i,j) += G_GL(k,j)*u_loc(NumDofs_Y(i,k,ext));
		  Ez_GL(i,j) += G_GL(k,ext)*u_loc(NumDofs_Z(i,j,k));
		  dEx_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_X(k,i,j));
		  dEy_G(k,i) += dGL_GL(j,ext)*u_loc(NumDofs_Y(i,k,j));
		}
	    }
	
	// second step : derivative
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      Ey = 0; Ex = 0; Ez = 0;
	      for (int k = 0; k <= order; k++)
		{
		  Ez -= dGL_GL(k,j)*Ex_GL(i,k);
		  Ez += dGL_GL(k,i)*Ey_GL(k,j);
		  Ex += dGL_GL(k,j)*Ez_GL(i,k);
		  Ey -= dGL_GL(k,i)*Ez_GL(k,j);
		}
	      
	      for (int k = 0; k < order; k++)
		{
		  Ey += G_GL(k,i)*dEx_G(k,j);
		  Ex -= G_GL(k,j)*dEy_G(k,i);
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
  void HexahedronHcurlFirstFamily::SolveMassMatrixGen(Vector1& x) const
  {
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= invWeightsMassG(i);
  }
  

  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix  
  template<class Vector1>
  void HexahedronHcurlFirstFamily::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= invSqrtWeightsMassG(i);
  }
    

  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector2>
  void HexahedronHcurlFirstFamily::MltMassMatrixGen(Vector2& x) const
  {
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= this->weights_dof(i);
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
  void HexahedronHcurlFirstFamily::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
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
  void HexahedronHcurlFirstFamily::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHcurlFirstFamily::ApplyRhQuadratureGen(const Vector1& Vh, Vector2& Uh) const
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
  void HexahedronHcurlFirstFamily::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHcurlFirstFamily::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
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
  void HexahedronHcurlFirstFamily::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHcurlFirstFamily::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    T vloc;
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int ix = NumDofs_X(i, j, k);
            int iy = NumDofs_Y(j, i, k);
            int iz = NumDofs_Z(j, k, i);
            
            // diagonal part
	    A.AddInteraction(m+ix, n+ix, C(0, 0)*this->weights_dof(ix));
	    A.AddInteraction(m+iy, n+iy, C(1, 1)*this->weights_dof(iy));
	    A.AddInteraction(m+iz, n+iz, C(2, 2)*this->weights_dof(iz));
	    
            // extra-diagonal part
            for (int ip = 0; ip < order; ip++)
              for (int jp = 0; jp <= order; jp++)
                {
                  vloc = weights1d(jp)*weights1d(j)*G_GL(i, jp)*G_GL(ip, j)*weights1d(k);
                  ix = NumDofs_X(i, j, k);
                  int jy = NumDofs_Y(jp, ip, k);
		  A.AddInteraction(m+ix, n+jy, C(0, 1)*vloc);
		  A.AddInteraction(m+jy, n+ix, C(1, 0)*vloc);
                  
                  ix = NumDofs_X(i, k, j);
                  int jz = NumDofs_Z(jp, k, ip);
		  A.AddInteraction(m+ix, n+jz, C(0, 2)*vloc);
		  A.AddInteraction(m+jz, n+ix, C(2, 0)*vloc);

                  iy = NumDofs_Y(k, i, j);
                  jz = NumDofs_Z(k, jp, ip);
		  A.AddInteraction(m+iy, n+jz, C(1, 2)*vloc);
		  A.AddInteraction(m+jz, n+iy, C(2, 1)*vloc);
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
  void HexahedronHcurlFirstFamily::
  AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    T vloc;
    const VectReal_wp& weights1d = this->Weights1D();    
    int ix, iy, iz, jx, jy, jz;
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            for (int kp = 0; kp <= order; kp++)
              {
                // diagonal part
                iz = NumDofs_Z(j, k, i);
                jz = NumDofs_Z(j, kp, i);
		A.AddInteraction(m+iz, n+jz, C(0, 0)*gauss.Weights(i)*weights1d(j)*stiff1d(k, kp));
                
                iy = NumDofs_Y(j, i, k);
                jy = NumDofs_Y(j, i, kp);
		A.AddInteraction(m+iy, n+jy, C(0, 0)*gauss.Weights(i)*weights1d(j)*stiff1d(k, kp));
		
                ix = NumDofs_X(i, j, k);
                jx = NumDofs_X(i, j, kp);
		A.AddInteraction(m+ix, n+jx, C(1, 1)*gauss.Weights(i)*weights1d(j)*stiff1d(k, kp));
                
                iz = NumDofs_Z(k, j, i);
                jz = NumDofs_Z(kp, j, i);
		A.AddInteraction(m+iz, n+jz, C(1, 1)*gauss.Weights(i)*weights1d(j)*stiff1d(k, kp));
                
                ix = NumDofs_X(i, k, j);
                jx = NumDofs_X(i, kp, j);
		A.AddInteraction(m+ix, n+jx, C(2, 2)*gauss.Weights(i)*weights1d(j)*stiff1d(k, kp));
                
                iy = NumDofs_Y(k, i, j);
                jy = NumDofs_Y(kp, i, j);
		A.AddInteraction(m+iy, n+jy, C(2, 2)*gauss.Weights(i)*weights1d(j)*stiff1d(k, kp));
		
                for (int ip = 0; ip < order; ip++)
                  {
                    iz = NumDofs_Z(j, k, i);
                    jy = NumDofs_Y(j, ip, kp);
                    vloc = C(0, 0)*gauss.Weights(i)*gauss.Weights(ip)
                      *weights1d(j)*dGL_G(k, ip)*dGL_G(kp, i);
                    
		    A.AddInteraction(m+iz, n+jy, -vloc);
		    A.AddInteraction(m+jy, n+iz, -vloc);
                    
                    ix = NumDofs_X(i, j, k);
                    jz = NumDofs_Z(kp, j, ip);
                    vloc = C(1, 1)*gauss.Weights(i)*gauss.Weights(ip)
                      *weights1d(j)*dGL_G(k, ip)*dGL_G(kp, i);
                    
		    A.AddInteraction(m+ix, n+jz, -vloc);
		    A.AddInteraction(m+jz, n+ix, -vloc);
		    
                    ix = NumDofs_X(i, k, j);
                    jy = NumDofs_Y(kp, ip, j);
                    vloc = C(2, 2)*gauss.Weights(i)*gauss.Weights(ip)
                      *weights1d(j)*dGL_G(k, ip)*dGL_G(kp, i);
                    
		    A.AddInteraction(m+ix, n+jy, -vloc);
		    A.AddInteraction(m+jy, n+ix, -vloc);
                  }
              }
            
            // extra-diagonal part
            ix = NumDofs_X(i, j, k);
            iy = NumDofs_Y(j, i, k);
            iz = NumDofs_Z(j, k, i);
            for (int ip = 0; ip < order; ip++)
              for (int jp = 0; jp <= order; jp++)
                for (int kp = 0; kp <= order; kp++)
                  {
                    jx = NumDofs_X(ip, jp, kp);
                    jy = NumDofs_Y(jp, ip, kp);
                    jz = NumDofs_Z(jp, kp, ip);
                    
                    // c10
                    vloc = C(0, 1)*G_GL(ip, j)*weights1d(j)*dGL_GL(k, jp)
                      *weights1d(jp)*dGL_G(kp, i)*gauss.Weights(i);
                    
		    A.AddInteraction(m+iz, n+jx, vloc);                    
		    A.AddInteraction(m+jx, n+iz, vloc);
                                        
                    if (i == ip)
                      {
                        vloc = -C(0, 1)*gauss.Weights(i)*dGL_GL(jp, j)
                          *weights1d(j)*dGL_GL(k, kp)*weights1d(kp);
                        
			A.AddInteraction(m+iz, n+jz, vloc);
			A.AddInteraction(m+jz, n+iz, vloc);
                      }
                    
                    vloc = -C(0, 1)*G_GL(ip, j)*weights1d(j)*G_GL(i, jp)
                      *weights1d(jp)*stiff1d(k, kp);
                    
		    A.AddInteraction(m+iy, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iy, vloc);
                    
                    vloc = C(0, 1)*dGL_GL(jp, j)*weights1d(j)*G_GL(i, kp)
                      *weights1d(kp)*dGL_G(k, ip)*gauss.Weights(ip);
                    
		    A.AddInteraction(m+iy, n+jz, vloc);
		    A.AddInteraction(m+jz, n+iy, vloc);
                    
                    // c20
                    vloc = C(0, 2)*dGL_GL(jp, j)*weights1d(j)*dGL_G(k, ip)
                      *gauss.Weights(ip)*G_GL(i, kp)*weights1d(kp);
                    
		    A.AddInteraction(m+iz, n+jy, vloc);
		    A.AddInteraction(m+jy, n+iz, vloc);
                    
                    if (i == ip)
                      {
                        vloc = -C(0, 2)*gauss.Weights(i)*dGL_GL(jp, j)
                          *weights1d(j)*dGL_GL(k, kp)*weights1d(kp);
                        
			A.AddInteraction(m+iy, n+jy, vloc);
			A.AddInteraction(m+jy, n+iy, vloc);
                      }
                    
                    vloc = -C(0, 2)*G_GL(ip, j)*weights1d(j)*G_GL(i, kp)
                      *weights1d(kp)*stiff1d(jp, k);
                    
		    A.AddInteraction(m+iz, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iz, vloc);
                    
                    vloc = C(0, 2)*dGL_GL(k, kp)*weights1d(kp)*G_GL(ip, j)
                      *weights1d(j)*dGL_G(jp, i)*gauss.Weights(i);
                    
		    A.AddInteraction(m+iy, n+jx, vloc);
		    A.AddInteraction(m+jx, n+iy, vloc);
                    
                    // c21
                    vloc = C(1, 2)*dGL_GL(k, kp)*weights1d(kp)*dGL_G(jp, i)
                      *gauss.Weights(i)*G_GL(ip, j)*weights1d(j);
                    
		    A.AddInteraction(m+ix, n+jy, vloc);
		    A.AddInteraction(m+jy, n+ix, vloc);
                    
                    if (i == ip)
                      {
                        vloc = -C(1, 2)*gauss.Weights(i)*dGL_GL(jp, j)
                          *weights1d(j)*dGL_GL(k, kp)*weights1d(kp);
                        
			A.AddInteraction(m+ix, n+jx, vloc);
			A.AddInteraction(m+jx, n+ix, vloc);
                      }
                    
                    vloc = -C(1, 2)*G_GL(ip, k)*weights1d(k)*G_GL(i, kp)
                      *weights1d(kp)*stiff1d(jp, j);
                    
		    A.AddInteraction(m+iz, n+jy, vloc);
		    A.AddInteraction(m+jy, n+iz, vloc);
                    
                    vloc = C(1, 2)*dGL_GL(jp, k)*weights1d(k)*G_GL(i, kp)
                      *weights1d(kp)*dGL_G(j, ip)*gauss.Weights(ip);
                    
		    A.AddInteraction(m+iz, n+jx, vloc); 
		    A.AddInteraction(m+jx, n+iz, vloc);  
                  }
          }
    
  }
  
  
  //! adds \int A curl(varphi_j) curl(varphi_i)
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A tensor A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A curl(varphi_j) curl(varphi_i)
   */
  template<class T, class Prop>
  void HexahedronHcurlFirstFamily
  ::AddVariableStiffnessMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& A,
				  VirtualMatrix<T>& mat) const
  {
    bool non_symmetric = !mat.IsSymmetric();
    int num_dof2, num_dof;    
    T vloc, zero; SetComplexZero(zero);
    
    const Array3D<int>& NumNodes3D = this->Fb_geom.GetNumNodes3D();
    for (int j1 = 0; j1 < order; j1++)
      for (int j2 = 0; j2 <= order; j2++)
        for (int j3 = 0; j3 <= order; j3++)
          {
            ////////////////////
            // CURL CURL PART //
            ////////////////////
            
            // interaction rot(phi_{j1,j2,j3}^1) rot(phi_{k1,k2,k3}^1)
            num_dof = NumDofs_X(j1, j2, j3);
            for (int k1 = 0;  k1 < order; k1++)
              for (int k2 = 0; k2 <= order ; k2++)
                for (int k3 = 0; k3 <= order ; k3++)
                  {
                    num_dof2 = NumDofs_X(k1, k2, k3);
                    vloc = zero;
                    if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                      {
                        if (j2 == k2)
                          for (int m = 0;  m <= order; m++)
                            for (int n = 0; n <= order; n++)
                              vloc += A(NumNodes3D(m,j2,n))(1,1)
                                *G_GL(j1,m)*dGL_GL(j3,n)*G_GL(k1,m)*dGL_GL(k3,n);
                        
                        for (int m = 0;  m <= order; m++)
                          {
                            vloc -= A(NumNodes3D(m,k2,j3))(2,1)
                              *G_GL(j1,m)*dGL_GL(j2,k2)*G_GL(k1,m)*dGL_GL(k3,j3);
                            vloc -= A(NumNodes3D(m,j2,k3))(2,1)
                              *G_GL(j1,m)*dGL_GL(j3,k3)*G_GL(k1,m)*dGL_GL(k2,j2);
                          }
                        
                        if (j3 == k3)
                          for (int m = 0;  m <= order; m++)
                            for (int n = 0; n <= order; n++)
                              vloc += A(NumNodes3D(m,n,j3))(2,2)
                                *G_GL(j1,m)*dGL_GL(j2,n)*G_GL(k1,m)*dGL_GL(k2, n);
                        
                        mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                      }
                  }
            
            // interaction rot(phi_{j1,j2,j3}^1) rot(phi_{k1,k2,k3}^2)
            for (int k1 = 0; k1 <= order; k1++)
              for (int k2 = 0; k2 < order; k2++)
                for (int k3 = 0; k3 <= order; k3++)
                  {
                    num_dof2 = NumDofs_Y(k1,k2,k3);
                    vloc = zero;
                    for (int m = 0; m <= order; m++)
                      {
                        vloc -= A(NumNodes3D(k1,j2,m))(1,0)
                          *G_GL(j1,k1)*dGL_GL(j3,m)*G_GL(k2,j2)*dGL_GL(k3,m);
                        vloc += A(NumNodes3D(k1,m,j3))(2,0)
                          *G_GL(j1,k1)*dGL_GL(j2,m)*G_GL(k2,m)*dGL_GL(k3,j3);
                        vloc += A(NumNodes3D(m,j2,k3))(2,1)
                          *G_GL(j1,m)*dGL_GL(j3,k3)*G_GL(k2,j2)*dGL_GL(k1,m);
                      }
                    
                    if (j3 == k3)
                      for (int m = 0;  m <= order; m++)
                        for (int n = 0; n <= order; n++)
                          vloc -= A(NumNodes3D(m,n,j3))(2,2)
                            *G_GL(j1,m)*dGL_GL(j2,n)*G_GL(k2,n)*dGL_GL(k1,m);
                    
		    mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
		    mat.AddInteraction(off_row+num_dof2, off_col+num_dof, vloc);
                  }
            
            // interaction rot(phi_{j1,j2,j3}^1) rot(phi_{k1,k2,k3}^3)
            for (int k1 = 0; k1 <= order; k1++)
              for (int k2 = 0; k2 <= order; k2++)
                for (int k3 = 0; k3 < order; k3++)
                  {
                    num_dof2 = NumDofs_Z(k1,k2,k3);
                    vloc = zero;
                    for (int m = 0; m <= order; m++)
                      {
                        vloc += A(NumNodes3D(k1,j2,m))(1,0)
                          *G_GL(j1,k1)*dGL_GL(j3,m)*G_GL(k3,m)*dGL_GL(k2,j2);
                        vloc -= A(NumNodes3D(k1,m,j3))(2,0)
                          *G_GL(j1,k1)*dGL_GL(j2,m)*G_GL(k3,j3)*dGL_GL(k2,m);
                        vloc += A(NumNodes3D(m,k2,j3))(2,1)
                          *G_GL(j1,m)*dGL_GL(j2,k2)*G_GL(k3,j3)*dGL_GL(k1,m);
                      }
                    
                    if (j2 == k2)
                      for (int m = 0;  m <= order; m++)
                        for (int n = 0; n <= order; n++)
                          vloc -= A(NumNodes3D(m,j2,n))(1,1)
                            *G_GL(j1,m)*dGL_GL(j3,n)*G_GL(k3,n)*dGL_GL(k1,m);
                    
		    mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
		    mat.AddInteraction(off_row+num_dof2, off_col+num_dof, vloc);
                  }
          }
    
    for (int j1 = 0; j1 <= order; j1++)
      for (int j2 = 0; j2 < order; j2++)
        for (int j3 = 0; j3 <= order; j3++)
          {
            // interaction rot(phi_{j1,j2,j3}^2) rot(phi_{k1,k2,k3}^2)
            num_dof = NumDofs_Y(j1,j2,j3);
            for (int k1 = 0; k1 <= order; k1++)
              for (int k2 = 0; k2 < order; k2++)
                for (int k3 = 0; k3 <= order; k3++)
                  {
                    num_dof2 = NumDofs_Y(k1,k2,k3);
                    vloc = zero;
                    
                    if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                      {
                        if (j1 == k1)
                          for (int m = 0;  m <= order; m++)
                            for (int n = 0; n <= order; n++)
                              vloc += A(NumNodes3D(j1,m,n))(0,0)
                                *G_GL(j2,m)*dGL_GL(j3,n)*G_GL(k2,m)*dGL_GL(k3,n);
                        
                        for (int m = 0;  m <= order; m++)
                          {
                            vloc -= A(NumNodes3D(k1,m,j3))(2,0)
                              *G_GL(j2,m)*dGL_GL(j1,k1)*G_GL(k2,m)*dGL_GL(k3,j3);
                            vloc -= A(NumNodes3D(j1,m,k3))(2,0)
                              *G_GL(j2,m)*dGL_GL(j3,k3)*G_GL(k2,m)*dGL_GL(k1,j1);
                          }
                        
                        if (j3 == k3)
                          for (int m = 0;  m <= order; m++)
                            for (int n = 0; n <= order; n++)
                              vloc += A(NumNodes3D(m,n,j3))(2,2)
                                *G_GL(j2,n)*dGL_GL(j1,m)*G_GL(k2,n)*dGL_GL(k1,m);
                        
                        mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                      }
                  }
            
            // interaction rot(phi_{j1,j2,j3}^2) rot(phi_{k1,k2,k3}^3)
            for (int k1 = 0; k1 <= order; k1++)
              for (int k2 = 0; k2 <= order; k2++)
                for (int k3 = 0; k3 < order; k3++)
                  {
                    num_dof2 = NumDofs_Z(k1,k2,k3);
                    vloc = zero;
                    for (int m = 0; m <= order; m++)
                      {
                        vloc += A(NumNodes3D(k1,m,j3))(2,0)
                          *G_GL(j2,m)*dGL_GL(j1,k1)*G_GL(k3,j3)*dGL_GL(k2,m);
                        vloc += A(NumNodes3D(j1,k2,m))(1,0)
                          *G_GL(j2,k2)*dGL_GL(j3,m)*G_GL(k3,m)*dGL_GL(k1,j1);
                        vloc -= A(NumNodes3D(m,k2,j3))(2,1)
                          *G_GL(j2,k2)*dGL_GL(j1,m)*G_GL(k3,j3)*dGL_GL(k1,m);
                      }
                    
                    if (j1 == k1)
                      for (int m = 0;  m <= order; m++)
                        for (int n = 0; n <= order; n++)
                          vloc -= A(NumNodes3D(j1,m,n))(0,0)
                            *G_GL(j2,m)*dGL_GL(j3,n)*G_GL(k3,n)*dGL_GL(k2,m);
                                        
		    mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
		    mat.AddInteraction(off_row+num_dof2, off_col+num_dof, vloc);
                  }
          }
    
    for (int j1 = 0; j1 <= order; j1++)
      for (int j2 = 0; j2 <= order; j2++)
        for (int j3 = 0; j3 < order; j3++)
          {
            // interaction rot(phi_{j1,j2,j3}^3) rot(phi_{k1,k2,k3}^3)
            num_dof = NumDofs_Z(j1,j2,j3);
            for (int k1 = 0; k1 <= order; k1++)
              for (int k2 = 0; k2 <= order; k2++)
                for (int k3 = 0; k3 < order; k3++)
                  {
                    num_dof2 = NumDofs_Z(k1,k2,k3);
                    vloc = zero;
                    
                    if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                      {
                        if (j1 == k1)
                          for (int m = 0;  m <= order; m++)
                            for (int n = 0; n <= order; n++)
                              vloc += A(NumNodes3D(j1,m,n))(0,0)
                                *G_GL(j3,n)*dGL_GL(j2,m)*G_GL(k3,n)*dGL_GL(k2,m);
                        
                        for (int m = 0;  m <= order; m++)
                          {
                            vloc -= A(NumNodes3D(k1,j2,m))(1,0)
                              *G_GL(j3,m)*dGL_GL(j1,k1)*G_GL(k3,m)*dGL_GL(k2,j2);
                            vloc -= A(NumNodes3D(j1,k2,m))(1,0)
                              *G_GL(j3,m)*dGL_GL(j2,k2)*G_GL(k3,m)*dGL_GL(k1,j1);
                          }
                        
                        if (j2 == k2)
                          for (int m = 0;  m <= order; m++)
                            for (int n = 0; n <= order; n++)
                              vloc += A(NumNodes3D(m,j2,n))(1,1)
                                *G_GL(j3,n)*dGL_GL(j1,m)*G_GL(k3,n)*dGL_GL(k1,m);
                        
                        mat.AddInteraction(off_row+num_dof,off_col+num_dof2, vloc);
                      }
                  }
          }
  }
  
  
  //! adds \int B varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] B tensor B on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int B varphi_j  varphi_i
   */
  template<class T, class Prop>
  void HexahedronHcurlFirstFamily
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& B,
			     VirtualMatrix<T>& mat) const
  {    
    int num_dof2, num_dof;    
    bool non_symmetric = !mat.IsSymmetric();
    T vloc, zero; SetComplexZero(zero);

    const Array3D<int>& NumNodes3D = this->Fb_geom.GetNumNodes3D();    
    for (int l = 0; l < order; l++)
      for (int m = 0; m <= order; m++)
        for (int n = 0; n <= order; n++)
          { 
            //////////////////////
            // MASS MATRIX PART //
            //////////////////////
            
            // interaction phi_{i,j,k}^1 phi_{l,m,n}^1
            num_dof2 = NumDofs_X(l,m,n);
            for (int i = 0; i < order; i++)
              {
                num_dof = NumDofs_X(i,m,n); vloc = zero;
                if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                  {
                    for (int p = 0; p <= order; p++)
                      vloc += B(NumNodes3D(p,m,n))(0,0)*G_GL(i,p)*G_GL(l,p);
                    
                    mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                  }
              }
            
            // interaction phi_{i,j,k}^1 phi_{m,l,n}^2
            num_dof2 = NumDofs_Y(m,l,n);
            for (int i = 0; i < order; i++)
              for (int j = 0; j <= order; j++)
                {
                  num_dof = NumDofs_X(i,j,n); vloc = zero;
                  if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                    {
                      vloc += B(NumNodes3D(m,j,n))(1,0)*G_GL(i,m)*G_GL(l,j);
                      
                      mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                    }
                }
            
            // interaction phi_{i,j,k}^1 phi_{n,m,l}^3
            num_dof2 = NumDofs_Z(n,m,l);
            for (int i = 0; i < order; i++)
              for (int k = 0; k <= order; k++)
                {
                  num_dof = NumDofs_X(i,m,k); vloc = zero;
                  if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                    {
                      vloc += B(NumNodes3D(n,m,k))(2,0)*G_GL(i,n)*G_GL(l,k);
                      
                      mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                    }
                }
            
            // interaction phi_{j,i,k}^2 phi_{l,m,n}^1
            num_dof2 = NumDofs_X(l,m,n);
            for (int i = 0; i < order; i++)
              for (int j = 0; j <= order; j++)
                {
                  num_dof = NumDofs_Y(j,i,n); vloc = zero;
                  if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                    {
                      vloc += B(NumNodes3D(j,m,n))(1,0)*G_GL(i,m)*G_GL(l,j);
                      
                      mat.AddInteraction(off_row+num_dof,off_col+num_dof2, vloc);
                    }
                }
            
            // interaction phi_{j,i,k}^2 phi_{m,l,n}^2
            num_dof2 = NumDofs_Y(m,l,n);
            for (int i = 0; i < order; i++)
              {
                num_dof = NumDofs_Y(m,i,n); vloc = zero;
                if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                  {
                    for (int p = 0; p <= order; p++)
                      vloc += B(NumNodes3D(m,p,n))(1,1)*G_GL(i,p)*G_GL(l,p);
                    
                    mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                  }
              }
            
            // interaction phi_{j,i,k}^2 phi_{n,m,l}^3
            num_dof2 = NumDofs_Z(n,m,l);
            for (int i = 0; i < order; i++)
              for (int k = 0; k <= order; k++)
                {
                  num_dof = NumDofs_Y(n,i,k); vloc = zero;
                  if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                    {
                      vloc += B(NumNodes3D(n,m,k))(2,1)*G_GL(i,m)*G_GL(l,k);
                      
                      mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                    }
                }
            
            // interaction phi_{k,j,i}^3 phi_{l,m,n}^1
            num_dof2 = NumDofs_X(l,m,n);
            for (int i = 0; i < order; i++)
              for (int k = 0; k <= order; k++)
                {
                  num_dof = NumDofs_Z(k,m,i); vloc = zero;
                  if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                    {
                      vloc += B(NumNodes3D(k,m,n))(2,0)*G_GL(i,n)*G_GL(l,k);
                      
                      mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                    }
                }
            
            // interaction phi_{k,j,i}^3 phi_{m,l,n}^2
            num_dof2 = NumDofs_Y(m,l,n);
            for (int i = 0; i < order; i++)
              for (int j = 0; j <= order; j++)
                {
                  num_dof = NumDofs_Z(m,j,i); vloc = zero;
                  if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                    {
                      vloc += B(NumNodes3D(m,j,n))(2,1)*G_GL(i,n)*G_GL(l,j);
                      
                      mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                    }
                }
            
            // interaction phi_{k,j,i}^3 phi_{n,m,l}^3
            num_dof2 = NumDofs_Z(n,m,l);
            for (int i = 0; i < order; i++)
              {
                num_dof = NumDofs_Z(n,m,i); vloc = zero;
                if (non_symmetric || (off_row+num_dof <= off_col+num_dof2))
                  {
                    for (int p = 0; p <= order; p++)
                      vloc += B(NumNodes3D(n,m,p))(2,2)*G_GL(i,p)*G_GL(l,p);
                    
                    mat.AddInteraction(off_row+num_dof, off_col+num_dof2, vloc);
                  }
              }
          } // end loop on l,m and n
  }
  

  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point
  */
  void HexahedronHcurlFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(R3(0,0,0));
    VectReal_wp phiGx(order), phiGy(order), phiGz(order),
      phiGLx(order+1), phiGLy(order+1), phiGLz(order+1);
    for (int i = 0; i < order; i++)
      {
	phiGx(i) = gauss.EvaluatePhi(i, point_loc(0));
	phiGy(i) = gauss.EvaluatePhi(i, point_loc(1));
	phiGz(i) = gauss.EvaluatePhi(i, point_loc(2));
      }

    for (int i = 0; i <= order; i++)
      {
	phiGLx(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiGLy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiGLz(i) = lob_basis.EvaluatePhi(i, point_loc(2));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
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
  void HexahedronHcurlFirstFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& curl) const
  {
    curl.Reallocate(nb_dof_loc); curl.Fill(R3());
    VectReal_wp phiGx(order), phiGy(order), phiGz(order),
      phiGLx(order+1), phiGLy(order+1), phiGLz(order+1);
    VectReal_wp dphiGLx(order+1), dphiGLy(order+1), dphiGLz(order+1);
    for (int i = 0; i < order; i++)
      {
	phiGx(i) = gauss.EvaluatePhi(i, point_loc(0));
	phiGy(i) = gauss.EvaluatePhi(i, point_loc(1));
	phiGz(i) = gauss.EvaluatePhi(i, point_loc(2));
      }

    for (int i = 0; i <= order; i++)
      {
	phiGLx(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiGLy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiGLz(i) = lob_basis.EvaluatePhi(i, point_loc(2));

	dphiGLx(i) = lob_basis.EvaluatePhiGrad(i, point_loc(0));
	dphiGLy(i) = lob_basis.EvaluatePhiGrad(i, point_loc(1));
	dphiGLz(i) = lob_basis.EvaluatePhiGrad(i, point_loc(2));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
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
  void HexahedronHcurlFirstFamily::GetValuePhiOnQuadraturePoint(int num_point, VectR3& phi) const
  {
    // for example, x-basis functions :
    // then \phi_{node}(\xi_{num_point})= \phi_{i1}^G(\xi_{k1}^{GL}) \phi_{i2}^{GL} (\xi_{k2}^{GL})
    //                                       \phi_{i3}^{GL} (\xi_{k3}^{GL}) e_x
    //                                    = G_GL(i1,k1) \delta_{i2,k2} \delta_{i3,k3} e_x
    
    phi.Reallocate(nb_dof_loc); FillZero(phi);
    int k1 = CoordinateQuad3D(num_point, 0);
    int k2 = CoordinateQuad3D(num_point, 1);
    int k3 = CoordinateQuad3D(num_point, 2);
    
    for (int i = 0; i < order; i++)
      {
        phi(NumDofs_X(i, k2, k3))(0) = G_GL(i, k1);
        phi(NumDofs_Y(k1, i, k3))(1) = G_GL(i, k2);
        phi(NumDofs_Z(k1, k2, i))(2) = G_GL(i, k3);
      }
    
  }
  
  
  //! computation of curl_phi = \f$ \nabla \times \phi_{node} (\xi_{num\_point}) \f$ 
  /*!
    \param[in] num_point quadrature point number
    \param[out] curl_phi curl of basis functions 
  */ 
  void HexahedronHcurlFirstFamily
  ::GetCurlPhiOnQuadraturePoint(int num_point, VectR3& curl_phi) const
  {
    curl_phi.Reallocate(nb_dof_loc); FillZero(curl_phi);
    int k1 = CoordinateQuad3D(num_point, 0);
    int k2 = CoordinateQuad3D(num_point, 1);
    int k3 = CoordinateQuad3D(num_point, 2);
        
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
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
          
          curl_phi(NumDofs_X(i, k2, j))(1) = dGL_GL(j, k3)*G_GL(i, k1);
          curl_phi(NumDofs_X(i, j, k3))(2) = -dGL_GL(j, k2)*G_GL(i, k1);
          
          // for y-basis functions
          // curl(\phi_{j,i,k}) (\xi_{num_point}) =
          // |  -\phi_j^{GL} (\xi_k1^{GL}) \phi_i^G (\xi_k2^{GL})  \phi_k^{'GL} (\xi_k3^{GL}) 
          // |  0
          // |  \phi_j^{'GL} (\xi_k1^{GL}) \phi_i^G (\xi_k2^{GL})  \phi_k^{GL} (\xi_k3^{GL})
          //
          // = |  -G_GL(i,k2) \delta_{j,k1} dGL_GL(k,k3)
          //   |  0
          //   |  G_GL(i,k2) \delta_{k,k3} dGL_GL(j,k1)
    
          curl_phi(NumDofs_Y(k1, i, j))(0) = -dGL_GL(j, k3)*G_GL(i, k2);
          curl_phi(NumDofs_Y(j, i, k3))(2) = dGL_GL(j, k1)*G_GL(i, k2);

          // for z-basis functions
          // curl(\phi_{j,k,i}) (\xi_{num_point}) =
          // |  \phi_j^{GL} (\xi_k1^{GL}) \phi_k^{'GL} (\xi_k2^{GL})  \phi_i^G (\xi_k3^{GL}) 
          // |  -\phi_j^{'GL} (\xi_k1^{GL}) \phi_k^{GL} (\xi_k2^{GL}) \phi_i^G (\xi_k3^{GL})
          // |  0
          //
          // = |  G_GL(i,k3) \delta_{j,k1} dGL_GL(k,k2)
          //   |  -G_GL(i,k3) \delta_{k,k2} dGL_GL(j,k1)
          //   |  0
          
          curl_phi(NumDofs_Z(k1, j, i))(0) = dGL_GL(j, k2)*G_GL(i, k3);
          curl_phi(NumDofs_Z(j, k2, i))(1) = -dGL_GL(j, k1)*G_GL(i, k3);
        }
    
  }
    

  //! computation of prolongation operator using low-order cells inside an element
  /*!
    \param[in,out] LocalProlongation prolongation operator
    \param[in] FaceCoarse coarse finite element
    prolongation operator is computed by using subdivided cube and low-order approximation
   */
  void HexahedronHcurlFirstFamily::
  ComputeLocalProlongationLowOrder(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
				   const ElementReference_Dim<Dimension3>& FaceCoarse) const
  {
    int nb_dof_loc = this->nb_dof_loc;
    int nb_dof_coarse = FaceCoarse.GetNbDof();
    LocalProlongation.Reallocate(nb_dof_loc, nb_dof_coarse);
    LocalProlongation.Fill(0);
    Real_wp eps = 1e3*epsilon_machine; R3 phi_coarse;
    HexahedronHcurlFirstFamily FaceOrder1;
    FaceOrder1.ConstructQuadrature(1, 0);
    FaceOrder1.ConstructFunctions();
    
    VectReal_wp Points(this->PointsNodal1D());
    VectReal_wp PoidsCoarse(FaceCoarse.GetNbDof());
    VectReal_wp PoidsFine(this->nb_dof_loc);
    R3 ext1, ext2, pt_coarse, phi; IVect num(4);
    
    for (int i1 = 0; i1 < this->order; i1++)
      for (int i2 = 0; i2 < this->order; i2++)
	for (int i3 = 0; i3 < this->order; i3++)
	  {
	    ext1(0) = Points(i1)-eps; ext1(1) = Points(i2)-eps;
	    ext1(2) = Points(i3)-eps;
	    ext2(0) = Points(i1+1)+eps; ext2(1) = Points(i2+1)+eps;
	    ext2(2) = Points(i3+1)+eps;
	    
	    for (int p_c = 0; p_c < 12; p_c++)
	      {
		int num_dof_fine = -1;
		for (int l = 0; l < 4; l++)
		  num(l) = FaceOrder1.CoordinateDofs(p_c,l);
		
		if (num(3) == 0)
		  {
		    num_dof_fine = this->NumDofs_X(i1, i2+num(1), i3+num(2));
		    PoidsFine(num_dof_fine) = 1.0/(ext2(0)-ext1(0));
		  }
		else if (num(3) == 1)
		  {
		    num_dof_fine = this->NumDofs_Y(i1+num(0), i2, i3+num(2));
		    PoidsFine(num_dof_fine) = 1.0/(ext2(1)-ext1(1));
		  }
		else
		  {
		    num_dof_fine = this->NumDofs_Z(i1+num(0), i2+num(1), i3);
		    PoidsFine(num_dof_fine) = 1.0/(ext2(2)-ext1(2));
		  }
		
	      }
	  }
    
    Points = FaceCoarse.PointsNodal1D();
    VectR3 val_phi;
    const HexahedronHcurlFirstFamily& Fc = static_cast<const HexahedronHcurlFirstFamily&>(FaceCoarse);
    const Array3D<int>& FcNumDofs_X = Fc.GetNumDofsX();
    const Array3D<int>& FcNumDofs_Y = Fc.GetNumDofsY();
    const Array3D<int>& FcNumDofs_Z = Fc.GetNumDofsZ();
    const Matrix<int>& FcCoordinateDofs = Fc.GetCoordinateDofs();
    for (int i1 = 0; i1 < FaceCoarse.GetOrder(); i1++)
      for (int i2 = 0; i2 < FaceCoarse.GetOrder(); i2++)
	for (int i3 = 0; i3 < FaceCoarse.GetOrder(); i3++)
	  {
	    ext1(0) = Points(i1)-eps; ext1(1) = Points(i2)-eps;
	    ext1(2) = Points(i3)-eps;
	    ext2(0) = Points(i1+1)+eps; ext2(1) = Points(i2+1)+eps;
	    ext2(2) = Points(i3+1)+eps;
	    // value = ext1(0)*ext2(0);
	    
	    for (int p_c = 0; p_c < 12; p_c++)
	      {
		int num_dof_coarse = -1;
		for (int l = 0; l < 4; l++)
		  num(l) = FaceOrder1.CoordinateDofs(p_c,l);
		
		if (num(3) == 0)
		  {
		    num_dof_coarse = FcNumDofs_X(i1, i2+num(1), i3+num(2));
		    PoidsCoarse(num_dof_coarse) = 1.0/(ext2(0)-ext1(0));
		  }
		else if (num(3) == 1)
		  {
		    num_dof_coarse = FcNumDofs_Y(i1+num(0), i2, i3+num(2));
		    PoidsCoarse(num_dof_coarse) = 1.0/(ext2(1)-ext1(1));
		  }
		else
		  {
		    num_dof_coarse = FcNumDofs_Z(i1+num(0), i2+num(1), i3);
		    PoidsCoarse(num_dof_coarse) = 1.0/(ext2(2)-ext1(2));
		  }
		
		for (int n = 0; n < nb_dof_loc; n++)
		  {
		    phi = this->PointsDofND(n);
		    if ((phi(0) >= ext1(0))&&(phi(0) <= ext2(0)))
		      if ((phi(1) >= ext1(1))&&(phi(1) <= ext2(1)))
			if ((phi(2) >= ext1(2))&&(phi(2) <= ext2(2)))
			  if (FcCoordinateDofs(num_dof_coarse, 3)
			      == this->CoordinateDofs(n,3))
			    {
			      for (int l = 0; l < 3; l++)
				pt_coarse(l) = (phi(l)-ext1(l))/(ext2(l)-ext1(l));
			      
			      FaceOrder1.ComputeValuesPhiRef(pt_coarse, val_phi);
			      LocalProlongation(n, num_dof_coarse) =
				PoidsCoarse(num_dof_coarse)/PoidsFine(n)*val_phi(p_c)(num(3));
			    }
		  }
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
  void HexahedronHcurlFirstFamily::
  ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
			   const ElementReference_Dim<Dimension3>& FaceCoarse,
			   const ElementReference_Dim<Dimension3>& FaceFine) const 
  {
    const HexahedronHcurlFirstFamily& Fc = static_cast<const HexahedronHcurlFirstFamily&>(FaceCoarse);
    R3 val; VectR3 val_phi;
    LocalProlongation.Reallocate(FaceFine.GetNbDof(), FaceCoarse.GetNbDof());
    for (int j = 0; j < FaceFine.GetNbDof(); j++)
      {
	Fc.ComputeValuesPhiRef(this->PointsDofND(j), val_phi);
	for (int k = 0; k < FaceCoarse.GetNbDof(); k++)
	  LocalProlongation(j, k) = val_phi(k)(this->CoordinateDofs(j,3));
      }
    
    const Array3D<int>& FcNumDofs_X = Fc.GetNumDofsX();
    const Array3D<int>& FcNumDofs_Y = Fc.GetNumDofsY();
    const Array3D<int>& FcNumDofs_Z = Fc.GetNumDofsZ();
    
    int rc = FaceCoarse.GetOrder(), rf = this->order;
    Matrix<Real_wp> G_GF(rc,rf), GL_GF(rc+1,rf+1);
    for (int i = 0; i < rc; i++)
      for (int j = 0; j < rf; j++)
	G_GF(i,j) = Fc.gauss.EvaluatePhi(i, gauss.Points(j));
    
    for (int i = 0; i <= rc; i++)
      for (int j = 0; j <= rf; j++)
	GL_GF(i,j) = Fc.lob_basis.EvaluatePhi(i, lob_basis.Points(j));
    
    Matrix<Real_wp, General, ArrayRowSparse> Ctmp;
    Ctmp.Reallocate(3*rf*(rc+1)*(rc+1), 3*rc*(rc+1)*(rc+1));
    for (int i = 0; i < rc; i++)
      for (int j = 0; j <= rc; j++)
	for (int k = 0; k <= rc; k++)
	  for (int m = 0; m < rf; m++)
	    {
	      int col = FcNumDofs_X(i,j,k);
	      int row = m*(rc+1)*(rc+1) + j*(rc+1) + k;
	      Ctmp.AddInteraction(row, col, G_GF(i,m));
	      
	      col = FcNumDofs_Y(j,i,k);
	      row = rf*(rc+1)*(rc+1) + m*(rc+1)*(rc+1) + j*(rc+1) + k;
	      Ctmp.AddInteraction(row, col, G_GF(i,m));
	      
	      col = FcNumDofs_Z(k,j,i);
	      row = 2*rf*(rc+1)*(rc+1) + m*(rc+1)*(rc+1) + j*(rc+1) + k;
	      Ctmp.AddInteraction(row, col, G_GF(i,m));
	    }
    
    TensorizedProjector<Dimension3>& true_proj = dynamic_cast<TensorizedProjector<Dimension3>& >(proj);
    true_proj.SetChOperator(0, Ctmp);
    
    Ctmp.Clear();
    Ctmp.Reallocate(3*rf*(rf+1)*(rc+1), 3*rf*(rc+1)*(rc+1));
    for (int i = 0; i < rf; i++)
      for (int j = 0; j <= rc; j++)
	for (int k = 0; k <= rc; k++)
	  for (int m = 0; m <= rf; m++)
	    {
	      int col = i*(rc+1)*(rc+1) + j*(rc+1) + k;
	      int row = i*(rf+1)*(rc+1) + m*(rc+1) + k;
	      Ctmp.AddInteraction(row, col, GL_GF(j,m));
	      
	      col = rf*(rc+1)*(rc+1) + i*(rc+1)*(rc+1) + j*(rc+1) + k;
	      row = rf*(rf+1)*(rc+1) + i*(rf+1)*(rc+1) + m*(rc+1) + k;
	      Ctmp.AddInteraction(row, col, GL_GF(j,m));
	      
	      col = 2*rf*(rc+1)*(rc+1) + i*(rc+1)*(rc+1) + j*(rc+1) + k;
	      row = 2*rf*(rf+1)*(rc+1) + i*(rf+1)*(rc+1) + m*(rc+1) + k;
	      Ctmp.AddInteraction(row, col, GL_GF(j,m));
	    }
    
    true_proj.SetChOperator(1, Ctmp);
    
    Ctmp.Clear();
    Ctmp.Reallocate(3*rf*(rf+1)*(rf+1), 3*rf*(rf+1)*(rc+1));
    for (int i = 0; i < rf; i++)
      for (int j = 0; j <= rf; j++)
	for (int k = 0; k <= rc; k++)
	  for (int m = 0; m <= rf; m++)
	    {
	      int col = i*(rf+1)*(rc+1) + j*(rc+1) + k;
	      int row = this->NumDofs_X(i,j,m);
	      Ctmp.AddInteraction(row, col, GL_GF(k,m));
	      
	      col = rf*(rf+1)*(rc+1) + i*(rf+1)*(rc+1) + j*(rc+1) + k;
	      row = this->NumDofs_Y(j,i,m);
	      Ctmp.AddInteraction(row, col, GL_GF(k,m));
	      
	      col = 2*rf*(rf+1)*(rc+1) + i*(rf+1)*(rc+1) + j*(rc+1) + k;
	      row = this->NumDofs_Z(m,j,i);
	      Ctmp.AddInteraction(row, col, GL_GF(k,m));
	    }
    
    true_proj.SetChOperator(2, Ctmp);
    true_proj.AllocateTemporaryVectors();
  }
  
  
  //! selects dofs near the current dof
  /*!
    \param[in] pos current dof
    \param[in] DofUsed some dofs are already in use and can't be selected
    \param[out] ListeDof selected dofs
    \param[in] nb_dof number of dofs to select
   */
  void HexahedronHcurlFirstFamily::PickNearDofs(int pos, const VectBool& DofUsed,
					     IVect& ListeDof, int nb_dof) const
  {
    if (nb_dof <= 0)
      return;
    
    int i0 = CoordinateDofs(pos,0); int j0 = CoordinateDofs(pos,1),
				      k0 = CoordinateDofs(pos,2);
    ListeDof.Reallocate(nb_dof); ListeDof.Fill(-1);
    int k = 1, nb = 0, node;
    if (i0 < order)
      if (!DofUsed(NumDofs_X(i0,j0,k0)))
	ListeDof(nb++) = NumDofs_X(i0,j0,k0);
    
    if (j0 < order)
      if (!DofUsed(NumDofs_Y(i0,j0,k0)))
	if (nb < nb_dof)
	  ListeDof(nb++) = NumDofs_Y(i0,j0,k0);
    
    if (k0 < order)
      if (!DofUsed(NumDofs_Z(i0,j0,k0)))
	if (nb < nb_dof)
	  ListeDof(nb++) = NumDofs_Z(i0,j0,k0);
	
    // DISP(DofUsed); DISP(i0); DISP(j0); DISP(k0); DISP(NumQuad3D);
    while (nb < nb_dof)
      {
	// loop on concentric cube at distance k
	if ((j0-k) >= 0)
	  {
	    for (int m = max(0,i0-k); m <= min(order-1,i0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_X(m,j0-k,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Y(m,j0-k,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order-1,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Z(m,j0-k,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	  }
	
	if (((j0+k) <= order)&&(nb < nb_dof))
	  {
	    for (int m = max(0,i0-k); m <= min(order-1,i0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_X(m,j0+k,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    if ((j0+k) < order)
	      for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
		for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		  if (nb < nb_dof)
		    {
		      node = NumDofs_Y(m,j0+k,n);
		      if (!DofUsed(node))
			ListeDof(nb++) = node;
		    }
	    
	    for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order-1,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Z(m,j0+k,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	  }
	
	if (((i0-k) >= 0)&&(nb < nb_dof))
	  {
	    for (int m = max(0,j0-k+1); m < min(order+1,j0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_X(i0-k,m,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Y(i0-k,m,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
		    
	    for (int m = max(0,j0-k+1); m < min(order+1,j0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order-1,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Z(i0-k,m,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	  }
	
	if (((i0+k) <= order)&&(nb < nb_dof))
	  {
	    if ((i0+k) < order)
	      for (int m = max(0,j0-k+1); m < min(order+1,j0+k); m++)
		for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		  if (nb < nb_dof)
		    {
		      node = NumDofs_X(i0+k,m,n);
		      if (!DofUsed(node))
			ListeDof(nb++) = node;
		    }
	    
	    for (int m = max(0,j0-k+1); m < min(order,j0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Y(i0+k,m,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    for (int m = max(0,j0-k+1); m < min(order+1,j0+k); m++)
	      for (int n = max(0,k0-k); n <= min(order-1,k0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Z(i0+k,m,n);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	  }
	
	if (((k0-k) >= 0)&&(nb < nb_dof))
	  {
	    for (int m = max(0,i0-k+1); m < min(order,i0+k); m++)
	      for (int n = max(0,j0-k+1); n < min(order+1,j0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_X(m,n,k0-k);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    for (int m = max(0,i0-k+1); m < min(order+1,i0+k); m++)
	      for (int n = max(0,j0-k+1); n < min(order,j0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Y(m,n,k0-k);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
		    
	    for (int m = max(0,i0-k+1); m < min(order+1,i0+k); m++)
	      for (int n = max(0,j0-k+1); n < min(order+1,j0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Z(m,n,k0-k);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	  }
	
	if (((k0+k) <= order)&&(nb < nb_dof))
	  {
	    for (int m = max(0,i0-k+1); m < min(order,i0+k); m++)
	      for (int n = max(0,j0-k+1); n < min(order+1,j0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_X(m,n,k0+k);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    for (int m = max(0,i0-k+1); m < min(order+1,i0+k); m++)
	      for (int n = max(0,j0-k+1); n < min(order,j0+k); n++)
		if (nb < nb_dof)
		  {
		    node = NumDofs_Y(m,n,k0+k);
		    if (!DofUsed(node))
		      ListeDof(nb++) = node;
		  }
	    
	    if ((k0+k) < order)
	      for (int m = max(0,i0-k+1); m < min(order+1,i0+k); m++)
		for (int n = max(0,j0-k+1); n < min(order+1,j0+k); n++)
		  if (nb < nb_dof)
		    {
		      node = NumDofs_Z(m,n,k0+k);
		      if (!DofUsed(node))
			ListeDof(nb++) = node;
		    }
	  }
	
	k++;
	if (k > (order+1))
	  {
	    nb = nb_dof;
	    DISP(ListeDof); DISP(DofUsed);
	    for (int j = 0; j < nb_dof_loc; j++)
	      if (!DofUsed(j))
		{
		  bool non_used_dof = true;
		  for (int k = 0; k < nb_dof; k++)
		    if (ListeDof(k) == j)
		      non_used_dof = false;
		  
		  if (non_used_dof)
		    {
		      DISP(i0); DISP(j0); DISP(k0);
		      DISP(j); DISP(CoordinateDofs(j,0));
		      DISP(CoordinateDofs(j,1)); DISP(CoordinateDofs(j,2));
		      DISP(CoordinateDofs(j,3));
		    }
		}
	    
	    exit(0);
	  }
	
      }
	
  }
  

  //! displays details of class HexahedronHcurlFirstFamily
  ostream& operator <<(ostream& out, const HexahedronHcurlFirstFamily& e)
  {
    out<<static_cast<const HexahedronReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_FIRST_FAMILY_CXX
#endif
