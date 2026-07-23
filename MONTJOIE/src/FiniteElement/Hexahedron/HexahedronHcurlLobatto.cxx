#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_LOBATTO_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHcurlLobatto::HexahedronHcurlLobatto() : HexahedronReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    this->mass_lumping = true;
    this->mass_lumping_ortho = true;
    
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  

  //! returns size of memory used by the object
  size_t HexahedronHcurlLobatto::GetMemorySize() const
  {
    size_t taille = HexahedronReference<2>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize() + NumDofs_Y.GetMemorySize() + NumDofs_Z.GetMemorySize();
    taille += Node_to_DofX.GetMemorySize() + Node_to_DofY.GetMemorySize() + Node_to_DofZ.GetMemorySize();
    taille += ch1_node.GetMemorySize() + ch2_node.GetMemorySize() + ch3_node.GetMemorySize()
      + rh_permut.GetMemorySize() + CoordinateDofs.GetMemorySize() + GL_GX.GetMemorySize()
      + rh_loc.GetMemorySize() + stiff1d.GetMemorySize();
    
    taille += this->weights_dof.GetMemorySize() + invWeightsMassG.GetMemorySize()
      + invSqrtWeightsMassG.GetMemorySize();
    
    return taille;
  }

  
  //! in order to number a mesh
  void HexahedronHcurlLobatto::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    // dg formulation
    nmap.SetFormulationDG(dg);
    
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return HexahedronReference<2>::ConstructNumberMap(nmap, dg);
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofQuadrangle(this->order, 2*(this->order+1)*(this->order+1));
	nmap.SetNbDofHexahedron(this->order, 0);
	
	// rotation of dofs on faces
	element_quad_surf->FindHcurlRotationQuad(nmap, 0);
	return;
      }
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order+1);
    nmap.SetNbDofQuadrangle(order, 2*(order+1)*(order-1));
    nmap.SetNbDofHexahedron(order, nb_dof_loc - 12*(order+1)
			    - 12*(order+1)*(order-1));
    
    nmap.SetOppositeEdgesDofSymmetry(order, order+1);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
    
    // rotation of dofs on faces
    Matrix<bool> OppositeSigneDof; //!< dof signs after rotation
    Matrix<int> FacesDof_Rotation;
    
    // for dofs inside the face, we compute FacesDof_Rotation
    int nb_dof_inside_quad = 2*(order-1)*(order+1);
    FacesDof_Rotation.Reallocate(8, nb_dof_inside_quad);
    OppositeSigneDof.Reallocate(8, nb_dof_inside_quad);
    OppositeSigneDof.Fill(false);
    for (int i = 0 ; i <= order ; i++)
      for (int j = 0 ; j < (order-1) ; j++)
	{
	  FacesDof_Rotation(0, j*(order+1)+i) = j*(order+1)+i;
	  FacesDof_Rotation(0, (order-1)*(order+1)+j*(order+1)+i)
            = (order-1)*(order+1)+j*(order+1)+i;
	  
	  // rotation of 90° of the face and direct way (x,y) -> (y,1-x)
	  //                                            (i,j) -> (j,r-2-i)
	  FacesDof_Rotation(1,j*(order+1)+i) = (order-1)*(order+1) 
	    + (order+1)*j + order-i;
	  FacesDof_Rotation(1,(order-1)*(order+1)+j*(order+1)+i) =
	    (order+1)*(order-2-j) + i;
	  OppositeSigneDof(1,j*(order+1)+i) = true;
	  OppositeSigneDof(1,(order-1)*(order+1)+j*(order+1)+i) = false;
	  
	  // rotation of 180° of the face and direct way (x,y) -> (1-x,1-y)
	  //                                             (i,j) -> (r-2-i,r-2-j)
	  FacesDof_Rotation(2,j*(order+1)+i) = (order+1)*(order-2-j) + order-i;
	  FacesDof_Rotation(2,(order-1)*(order+1)+j*(order+1)+i) =
	    (order-1)*(order+1) + (order+1)*(order-2-j) + order-i;
	  OppositeSigneDof(2,j*(order+1)+i) = true;
	  OppositeSigneDof(2,(order-1)*(order+1)+j*(order+1)+i) = true;
	  
	  // rotation of 270° of the face and direct way (x,y) -> (1-y,x)
	  //                                             (i,j) -> (r-2-j,i)
	  FacesDof_Rotation(3,j*(order+1)+i) = (order-1)*(order+1)
	    + (order+1)*(order-2-j) + i;
	  FacesDof_Rotation(3,(order-1)*(order+1)+j*(order+1)+i) =
	    (order+1)*j + order-i;
	  OppositeSigneDof(3, j*(order+1)+i) = false;
	  OppositeSigneDof(3, (order-1)*(order+1)+j*(order+1)+i) = true;
	  
	  // rotation of 0° of the face and opposite way (x,y) -> (y,x)
	  //                                             (i,j) -> (j,i)
	  FacesDof_Rotation(4, j*(order+1)+i) = (order-1)*(order+1) + (order+1)*j + i;
	  FacesDof_Rotation(4, (order-1)*(order+1)+j*(order+1)+i) = (order+1)*j + i;
	  OppositeSigneDof(4, j*(order+1)+i) = false;
	  OppositeSigneDof(4, (order-1)*(order+1)+j*(order+1)+i) = false;
	  
	  // rotation of 90° of the face and opposite way (x,y) -> (1-x,y)
	  //                                              (i,j) -> (r-2-i,j)
	  FacesDof_Rotation(5,j*(order+1)+i) = (order+1)*j + order-i;
	  FacesDof_Rotation(5,(order-1)*(order+1)+j*(order+1)+i) = 
	    (order-1)*(order+1)+(order+1)*(order-2-j) + i;
	  OppositeSigneDof(5, j*(order+1)+i) = true;
	  OppositeSigneDof(5, (order-1)*(order+1)+j*(order+1)+i) = false;
	  
	  // rotation of 180° of the face and opposite way (x,y) -> (1-y,1-x)
	  //                                               (i,j) -> (r-2-j,r-2-i)
	  FacesDof_Rotation(6,j*(order+1)+i) = (order-1)*(order+1)
	    + (order+1)*(order-2-j) + order-i;
	  FacesDof_Rotation(6,(order-1)*(order+1)+j*(order+1)+i) =
	    (order+1)*(order-2-j) + order-i;
	  OppositeSigneDof(6, j*(order+1)+i) = true;
	  OppositeSigneDof(6, (order-1)*(order+1)+j*(order+1)+i) = true;
	  
	  // rotation of 270° of the face and opposite way (x,y) -> (x,1-y)
	  //                                               (i,j) -> (i,r-2-j) 
	  FacesDof_Rotation(7,j*(order+1)+i) = (order+1)*(order-2-j) + i;
	  FacesDof_Rotation(7,(order-1)*(order+1)+j*(order+1)+i) = (order-1)*(order+1)
	    + (order+1)*j + order-i;
	  OppositeSigneDof(7, j*(order+1)+i) = false;
	  OppositeSigneDof(7, (order-1)*(order+1)+j*(order+1)+i) = true;
	  
	}
    
    nmap.SetFacesDofRotationQuad(order, FacesDof_Rotation);
    nmap.SetSignDofRotationQuad(order, OppositeSigneDof);
  }
  

  //! construction of finite element  
  void HexahedronHcurlLobatto::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						      int rsurf_tri, int rsurf_quad,
						      int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    HexahedronReference<2>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();

#ifdef MONTJOIE_HEXAHEDRON_HCURL_LOBATTO_GAUSS    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(r-1, gauss.QUADRATURE_GAUSS);
    points3d.Resize(nb_points_quadrature_inside + r*r*r);
    weights3d.Resize(nb_points_quadrature_inside + r*r*r);
    int nb = nb_points_quadrature_inside;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	for (int k = 0; k < r; k++)
	  {
	    points3d(nb).Init(gauss.Points(i), gauss.Points(j), gauss.Points(k));
	    weights3d(nb) = gauss.Weights(i)*gauss.Weights(j)*gauss.Weights(k);
	    nb++;
	  }
    
    nb_points_quadrature_inside = nb;
#endif
    
    FindDofsOnFace(true);
    
    this->Fb_geom.ComputeCoefficientTransformation();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    ConstructStiffnessMatrix();
    ConstructMassMatrix();

    QuadrangleHcurlLobatto* Fb_quad = new QuadrangleHcurlLobatto();
    Fb_quad->ConstructFiniteElement(r);
    element_quad_surf = Fb_quad;
  }

  
  //! construction of basis functions
  void HexahedronHcurlLobatto::ConstructFunctions()
  {
    nb_dof_quad = 2*(order+1)*(order+1);
    nb_dof_tri = 0;
    nb_dof_loc = 3*this->GetNbPointsNodalElt();
    int nb_dof_inside_quad = 2*(order+1)*(order-1);
    nb_dof_boundaries = 6*nb_dof_inside_quad + 12*(order+1);
    
    lob_basis.ConstructQuadrature(order, lob_basis.QUADRATURE_LOBATTO);
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    
    NumDofs_X.Reallocate(order+1, order+1, order+1);
    NumDofs_Y.Reallocate(order+1, order+1, order+1);
    NumDofs_Z.Reallocate(order+1, order+1, order+1);
    
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
	//int i2 = CoordinateNodes(num_vertex_1,0);
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
	  for (int l = 0 ; l <= order ; l++)
	    {
	      // we have an edge along Oz direction
	      NumDofs_Z(i,j,l) = nb_dof_loc;
	      nb_dof_loc++;
	    }
	else if (j2 != j)
	  for (int l = 0 ; l <= order ; l++)
	    {
	      // we have an edge along Oy direction
	      NumDofs_Y(i,l,k) = nb_dof_loc;
	      nb_dof_loc++;
	    }
	else
	  for (int l = 0; l <= order; l++)
	    {
	      // we have an edge along Ox direction
	      NumDofs_X(l,j,k) = nb_dof_loc;
	      nb_dof_loc++;
	    }
      }
    
    // dofs on first face
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i <= order ; i++)
	{
	  NumDofs_Y(0,i,j) = nb_dof_loc;
	  NumDofs_Z(0,j,i) = (order+1)*(order-1) + nb_dof_loc++;
	}
    nb_dof_loc += (order+1)*(order-1);
    // dofs on second face
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i <= order ; i++)
	{
	  NumDofs_X(i,0,j) = nb_dof_loc;
	  NumDofs_Z(j,0,i) = (order+1)*(order-1) + nb_dof_loc++;
	}
    nb_dof_loc += (order+1)*(order-1);
    // dofs on third face
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i <= order ; i++)
	{
	  NumDofs_X(i,j,0) = nb_dof_loc;
	  NumDofs_Y(j,i,0) = (order+1)*(order-1) + nb_dof_loc++;
	}
    nb_dof_loc += (order+1)*(order-1);
    // dofs on fourth face
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i <= order ; i++)
	{
	  NumDofs_X(i,j,order) = nb_dof_loc;
	  NumDofs_Y(j,i,order) = (order+1)*(order-1) + nb_dof_loc++;
	}
    nb_dof_loc += (order+1)*(order-1);
    // dofs on fifth face
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i <= order ; i++)
	{
	  NumDofs_X(i,order,j) = nb_dof_loc;
	  NumDofs_Z(j,order,i) = (order+1)*(order-1) + nb_dof_loc++;
	}
    nb_dof_loc += (order+1)*(order-1);
    // dofs on sixth face
    for (int j = 1 ; j < order ; j++)
      for (int i = 0; i <= order ; i++)
	{
	  NumDofs_Y(order,i,j) = nb_dof_loc;
	  NumDofs_Z(order,j,i) = (order+1)*(order-1)+nb_dof_loc++;
	}    
    nb_dof_loc += (order+1)*(order-1);
    // dofs inside the hexahedron
    for (int i = 0 ; i <= order ; i++)
      for (int j = 1 ; j < order ; j++)
	for (int k = 1 ; k < order ; k++)
	  {
	    NumDofs_X(i,j,k) = nb_dof_loc++;
	    NumDofs_Y(j,i,k) = nb_dof_loc++;
	    NumDofs_Z(k,j,i) = nb_dof_loc++;
	  }
    
    CoordinateDofs.Reallocate(nb_dof_loc, 4);
    
    for (int i = 0 ; i <= order ; i++)
      for (int j = 0 ; j <= order ; j++)
	for (int k = 0 ; k <= order ; k++)
	  {
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

    Dof_to_Node.Reallocate(nb_dof_loc);
    Node_to_DofX.Reallocate(nb_points_quadrature_inside);
    Node_to_DofY.Reallocate(nb_points_quadrature_inside);
    Node_to_DofZ.Reallocate(nb_points_quadrature_inside);
    for (int i = 0 ; i <= order ; i++)
      for (int j = 0 ; j <= order ; j++)
	for (int k = 0 ; k <= order ; k++)
	  {
	    Node_to_DofX(NumQuad3D(i,j,k)) = NumDofs_X(i,j,k);
	    Node_to_DofY(NumQuad3D(i,j,k)) = NumDofs_Y(i,j,k);
	    Node_to_DofZ(NumQuad3D(i,j,k)) = NumDofs_Z(i,j,k);
	    Dof_to_Node(NumDofs_X(i, j, k)) = NumQuad3D(i, j, k);
	    Dof_to_Node(NumDofs_Y(i, j, k)) = NumQuad3D(i, j, k);
	    Dof_to_Node(NumDofs_Z(i, j, k)) = NumQuad3D(i, j, k);
	  }
    
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(this->PointsND());

  }
  
    
  //! computation of stiffness matrix
  void HexahedronHcurlLobatto::ConstructStiffnessMatrix()
  {
    Matrix<Real_wp, General, ArrayRowSparse> Rh_curl(nb_dof_loc, nb_dof_loc);
    rh_permut.Reallocate(nb_dof_loc, nb_dof_loc);
    int num_dof, num_dof2, node; Real_wp vloc;
    const VectReal_wp& weights3d = this->WeightsND();
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  {
	    num_dof = NumDofs_X(i1, i2, i3);
	    for (int n = 0; n <= order; n++)
	      {
		node = NumQuad3D(i1, i2, n);
		num_dof2 = 3*node+1; vloc = lob_basis.GradPhi(i3,n);
		Rh_curl.AddInteraction(num_dof, num_dof2, vloc);		
		rh_permut.AddInteraction(num_dof, NumDofs_Y(i1, i2, n), vloc*weights3d(node));
				
		node = NumQuad3D(i1,n,i3);
		num_dof2 = 3*node+2; vloc = -lob_basis.GradPhi(i2,n);
		Rh_curl.AddInteraction(num_dof, num_dof2, vloc);
		rh_permut.AddInteraction(num_dof, NumDofs_Z(i1, n, i3), vloc*weights3d(node));
              }
	    
	    num_dof = NumDofs_Y(i1,i2,i3);
	    for (int n = 0; n <= order; n++)
	      {
		node = NumQuad3D(i1,i2,n);
		num_dof2 = 3*node; vloc = -lob_basis.GradPhi(i3,n);
		Rh_curl.AddInteraction(num_dof, num_dof2, vloc);
		rh_permut.AddInteraction(num_dof, NumDofs_X(i1, i2, n), vloc*weights3d(node));
				
		node = NumQuad3D(n,i2,i3);
		num_dof2 = 3*node+2; vloc = lob_basis.GradPhi(i1,n);
		Rh_curl.AddInteraction(num_dof, num_dof2, vloc);
		rh_permut.AddInteraction(num_dof, NumDofs_Z(n, i2, i3), vloc*weights3d(node));
              }
	    
	    num_dof = NumDofs_Z(i1,i2,i3);
	    for (int n = 0; n <= order; n++)
	      {
		node = NumQuad3D(i1,n,i3);
		num_dof2 = 3*node; vloc = lob_basis.GradPhi(i2,n);
		Rh_curl.AddInteraction(num_dof, num_dof2, vloc);
		rh_permut.AddInteraction(num_dof, NumDofs_X(i1, n, i3), vloc*weights3d(node));
				
		node = NumQuad3D(n,i2,i3);
		num_dof2 = 3*node+1; vloc = -lob_basis.GradPhi(i1,n);
		Rh_curl.AddInteraction(num_dof, num_dof2, vloc);
		rh_permut.AddInteraction(num_dof, NumDofs_Y(n, i2, i3), vloc*weights3d(node));
              }
	  }

    Seldon::Copy(Rh_curl, rh_loc);
    /*
    for (int num_loc = 0; num_loc < 6; num_loc++)
      {
	int nx = 0, ny = 1;
	if ((num_loc == 0) || (num_loc == 5))
	  { nx = 1; ny = 2; }
	else if ((num_loc == 1) || (num_loc == 4))
	  { nx = 0; ny = 2; }	

	for (int k = 0; k < (order+1)*(order+1); k++)
	  {
	    int npoint = this->GetQuadNumber(num_loc, k);	
	    int kx, ky, kz;
	    this->GetDofNumber_FromPointNode(npoint, kx, ky, kz);
	    
	    int jx = kx, jy = ky;
	    if ((num_loc == 0) || (num_loc == 5))
	      { jx = ky; jy = kz; }
	    else if ((num_loc == 1) || (num_loc == 4))
	      { jx = kx; jy = kz; }

	    Real_wp coef = this->WeightsQuadratureBoundary(k, num_loc);
	    if ((num_loc == 0) || (num_loc == 2) || (num_loc == 4))
	      coef = -coef;
	    
	    Rh.AddInteraction(jy, 3*npoint + nx, coef);
	    Rh.AddInteraction(jx, 3*npoint + ny, -coef);
	  }
      }

    Vector<int> permut(this->nb_dof_loc);
    for (int k = 0; k < this->nb_points_quadrature_inside; k++)
      {
	int kx, ky, kz;
	this->GetDofNumber_FromPointNode(k, kx, ky, kz);
	permut(kx) = 3*k;
	permut(ky) = 3*k+1;
	permut(kz) = 3*k+2;
      }

    permut.WriteText("permut.dat");
    */
    
    //Seldon::Copy(Rh, rh_weight);
    //Rh.WriteText("RhB.dat");

    stiff1d.Reallocate(order+1, order+1);
    stiff1d.Fill(0);    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          stiff1d(i, j) += this->Weights1D(k)*dGL_GL(i, k)*dGL_GL(j, k);
  }
  
  
  //! constructing mass matrix
  void HexahedronHcurlLobatto::ConstructMassMatrix()
  {
    Globatto<Real_wp> gauss_exact;
    gauss_exact.ConstructQuadrature(order, gauss_exact.QUADRATURE_GAUSS);
    GL_GX.Reallocate(order+1, order+1);
        
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        GL_GX(i,j) = lob_basis.EvaluatePhi(i, gauss_exact.Points(j));
    
    int order_geom = this->GetGeometryOrder();
    const VectReal_wp& weights3d = this->WeightsND();
    this->weights_dof.Reallocate(nb_dof_loc);
    invWeightsMassG.Reallocate(nb_dof_loc);
    invSqrtWeightsMassG.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        this->weights_dof(Node_to_DofX(i)) = weights3d(i);
        this->weights_dof(Node_to_DofY(i)) = weights3d(i);
        this->weights_dof(Node_to_DofZ(i)) = weights3d(i);
      }

    for (int num_dof = 0; num_dof < nb_dof_loc; num_dof++)
      {
        invWeightsMassG(num_dof) = 1.0/this->weights_dof(num_dof);
        invSqrtWeightsMassG(num_dof) = 1.0/sqrt(this->weights_dof(num_dof));
      }
    
    Matrix<Real_wp> GL_Geom(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        GL_Geom(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
    
    
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch1.Clear(); Ch2.Clear(); Ch3.Clear();
    Ch1.Reallocate(3*(order+1)*(order_geom+1)*(order_geom+1), 3*this->GetNbPointsNodalElt());
    Ch2.Reallocate(3*(order+1)*(order+1)*(order_geom+1),
                   3*(order+1)*(order_geom+1)*(order_geom+1));
    Ch3.Reallocate(3*(order+1)*(order+1)*(order+1), 3*(order+1)*(order+1)*(order_geom+1));
    const Array3D<int>& NumNodes3D = this->Fb_geom.GetNumNodes3D();
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k <= order_geom; k++)
	  for (int m = 0; m <= order; m++)
	    {
	      int node = 3*NumNodes3D(i, j, k);
	      int irow = 3*((order_geom+1)*( m*(order_geom+1) + j) + k);
	      Ch1.AddInteraction(irow, node, GL_Geom(m, i));
              Ch1.AddInteraction(irow+1, node+1, GL_Geom(m, i));
              Ch1.AddInteraction(irow+2, node+2, GL_Geom(m, i));
            }
        
    for (int i = 0; i <= order_geom; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m <= order; m++)
	    {
              int icol = 3*((order_geom+1)*(k*(order_geom+1)+i) + j);
              int irow = 3*((order_geom+1)*(k*(order+1) + m) + j);
	      Ch2.AddInteraction(irow, icol, GL_Geom(m, i));
              Ch2.AddInteraction(irow+1, icol+1, GL_Geom(m, i));
              Ch2.AddInteraction(irow+2, icol+2, GL_Geom(m, i));
            }	      	      
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
	for (int k = 0; k <= order; k++)
	  for (int m = 0; m <= order; m++)
	    {    
              int icol = (order_geom+1)*(k*(order+1) + i) + j;
	      Ch3.AddInteraction(NumDofs_X(k, i, m), 3*icol, GL_Geom(m, j));
              Ch3.AddInteraction(NumDofs_Y(k, i, m), 3*icol+1, GL_Geom(m, j));
              Ch3.AddInteraction(NumDofs_Z(k, i, m), 3*icol+2, GL_Geom(m, j));
            }
    
    Seldon::Copy(Ch1, ch1_node);
    Seldon::Copy(Ch2, ch2_node);
    Seldon::Copy(Ch3, ch3_node);

  }
  
  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are placed
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlLobatto::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    contrib.Reallocate(nb_dof_loc);
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  {
	    int node = NumQuad3D(i1,i2,i3);
	    contrib(NumDofs_X(i1,i2,i3)) = feval(3*node);
	    contrib(NumDofs_Y(i1,i2,i3)) = feval(3*node+1);
	    contrib(NumDofs_Z(i1,i2,i3)) = feval(3*node+2);
	  }
  }
  

  void HexahedronHcurlLobatto
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {
    for (int i = 2*(this->order+1); i < 4*(this->order+1); i++)
      contrib(i) = -contrib(i);
  }


  void HexahedronHcurlLobatto
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    for (int i = 2*(this->order+1); i < 4*(this->order+1); i++)
      contrib(i) = -contrib(i);
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
  void HexahedronHcurlLobatto::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    Vector1 fx(3*this->Points2D_quad().GetM());
    fx.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    // we express feval on Gauss-Lobatto points
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          int node = NumQuad2D(i, j);
          for (int m = 0; m <= order; m++)
            {
              fx(3*node) += GL_GX(i, m)*feval(3*NumQuad2D(m, j));
              fx(3*node+1) += GL_GX(i, m)*feval(3*NumQuad2D(m, j)+1);
              fx(3*node+2) += GL_GX(i, m)*feval(3*NumQuad2D(m, j)+2);
            }
        }

    Vector1 fy(3*this->Points2D_quad().GetM());
    fy.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          int node = NumQuad2D(i, j);
          for (int m = 0; m <= order; m++)
            {
              fy(3*node) += GL_GX(j, m)*fx(3*NumQuad2D(i, m));
              fy(3*node+1) += GL_GX(j, m)*fx(3*NumQuad2D(i, m)+1);
              fy(3*node+2) += GL_GX(j, m)*fx(3*NumQuad2D(i, m)+2);
            }
        }
    
    // then we use integration with Gauss-Lobatto points
    ComputeIntegralSurfaceRef(fy, res, num_loc);
  }
  

  //! computation of U on nodal points (on reference element)
  /*!
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
  */
  template<class Vector1, class Vector2>
  void HexahedronHcurlLobatto::
  ComputeNodalValuesGen(const Vector1& Uloc, Vector2& Uloc_node) const
  {    
    if (order == this->GetGeometryOrder())
      for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
        {
          Uloc_node(3*i) = Uloc(Node_to_DofX(i));
          Uloc_node(3*i+1) = Uloc(Node_to_DofY(i));
          Uloc_node(3*i+2) = Uloc(Node_to_DofZ(i));
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
    \param[out] u_boundary values of u on quadrature points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlLobatto::
  ComputeValueBoundaryGen(const Vector1& u_loc, Vector2& u_boundary, int num_loc) const
  {
    u_boundary.Reallocate(3*this->PointsNodal2D_quad().GetM());
    int ext = 0;
    if (num_loc > 2)
      ext = order;
    
    const Matrix<int>& NumNodes2D_quad = this->Fb_geom.GetNumNodes2D_quad();
    if ((num_loc == 0)||(num_loc==5))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    u_boundary(3*NumNodes2D_quad(i,j)) = u_loc(NumDofs_X(ext,i,j));
	    u_boundary(3*NumNodes2D_quad(i,j)+1) = u_loc(NumDofs_Y(ext,i,j));
	    u_boundary(3*NumNodes2D_quad(i,j)+2) = u_loc(NumDofs_Z(ext,i,j));
	  }
    else if ((num_loc == 1)||(num_loc==4))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    u_boundary(3*NumNodes2D_quad(i,j)) = u_loc(NumDofs_X(i,ext,j));
	    u_boundary(3*NumNodes2D_quad(i,j)+1) = u_loc(NumDofs_Y(i,ext,j));
	    u_boundary(3*NumNodes2D_quad(i,j)+2) = u_loc(NumDofs_Z(i,ext,j));
	  }
    else
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    u_boundary(3*NumNodes2D_quad(i,j)) = u_loc(NumDofs_X(i,j,ext));
	    u_boundary(3*NumNodes2D_quad(i,j)+1) = u_loc(NumDofs_Y(i,j,ext));
	    u_boundary(3*NumNodes2D_quad(i,j)+2) = u_loc(NumDofs_Z(i,j,ext));
	  }
  }
  
  
  //! computation of curl(u) on nodal points of a face
  /*!
    \param[in] u_loc dof components of u
    \param[out] u_boundary curl of u on quadrature points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1,class Vector2>
  void HexahedronHcurlLobatto::
  ComputeCurlBoundaryGen(const Vector1& u_loc, Vector2& u_boundary, int num_loc) const
  {
    u_boundary.Reallocate(3*this->PointsNodal2D_quad().GetM());
    typedef typename Vector1::value_type Complexe;
    TinyVector<Complexe,3> val;
    int ext = 0;
    if (num_loc > 2)
      ext = order;
    
    const Matrix<int>& NumNodes2D_quad = this->Fb_geom.GetNumNodes2D_quad();    
    const Matrix<Real_wp>& dGL = lob_basis.GradPhi();
    if ((num_loc == 0)||(num_loc==5))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    val.Zero();
	    for (int k = 0; k <= order; k++)
	      {
		val(1) += dGL(k,j)*u_loc(NumDofs_X(ext,i,k));
		val(2) -= dGL(k,i)*u_loc(NumDofs_X(ext,k,j));
		val(0) -= dGL(k,j)*u_loc(NumDofs_Y(ext,i,k));
		val(2) += dGL(k,ext)*u_loc(NumDofs_Y(k,i,j));
		val(0) += dGL(k,i)*u_loc(NumDofs_Z(ext,k,j));
		val(1) -= dGL(k,ext)*u_loc(NumDofs_Z(k,i,j));
	      }
	    u_boundary(3*NumNodes2D_quad(i,j)) = val(0);
            u_boundary(3*NumNodes2D_quad(i,j)+1) = val(1);
            u_boundary(3*NumNodes2D_quad(i,j)+2) = val(2);
	  }
    else if ((num_loc == 1)||(num_loc==4))
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    val.Zero();
	    for (int k = 0; k <= order; k++)
	      {
		val(1) += dGL(k,j)*u_loc(NumDofs_X(i,ext,k));
		val(2) -= dGL(k,ext)*u_loc(NumDofs_X(i,k,j));
		val(0) -= dGL(k,j)*u_loc(NumDofs_Y(i,ext,k));
		val(2) += dGL(k,i)*u_loc(NumDofs_Y(k,ext,j));
		val(0) += dGL(k,ext)*u_loc(NumDofs_Z(i,k,j));
		val(1) -= dGL(k,i)*u_loc(NumDofs_Z(k,ext,j));
	      }
            u_boundary(3*NumNodes2D_quad(i,j)) = val(0);
            u_boundary(3*NumNodes2D_quad(i,j)+1) = val(1);
            u_boundary(3*NumNodes2D_quad(i,j)+2) = val(2);
	  }
    else
      for (int i = 0; i <= order; i++)
	for (int j = 0; j <= order; j++)
	  {
	    val.Zero();
	    for (int k = 0; k <= order; k++)
	      {
		val(1) += dGL(k,ext)*u_loc(NumDofs_X(i,j,k));
		val(2) -= dGL(k,j)*u_loc(NumDofs_X(i,k,ext));
		val(0) -= dGL(k,ext)*u_loc(NumDofs_Y(i,j,k));
		val(2) += dGL(k,i)*u_loc(NumDofs_Y(k,j,ext));
		val(0) += dGL(k,j)*u_loc(NumDofs_Z(i,k,ext));
		val(1) -= dGL(k,i)*u_loc(NumDofs_Z(k,j,ext));
	      }
            u_boundary(3*NumNodes2D_quad(i,j)) = val(0);
            u_boundary(3*NumNodes2D_quad(i,j)+1) = val(1);
            u_boundary(3*NumNodes2D_quad(i,j)+2) = val(2);
          }
  }
  

  //! x is overwritten by M^-1 x where M is the mass matrix    
  template<class Vector1>
  void HexahedronHcurlLobatto::SolveMassMatrixGen(Vector1& x) const
  {
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= invWeightsMassG(i);
  }
  

  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix    
  template<class Vector1>
  void HexahedronHcurlLobatto::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= invSqrtWeightsMassG(i);
  }
    
  
  //! x is overwritten by M x where M is the mass matrix  
  template<class Vector2>
  void HexahedronHcurlLobatto::MltMassMatrixGen(Vector2& x) const
  {
    for (int i = 0; i < x.GetM(); i++)
      x(i) *= this->weights_dof(i);
  }


  //! Integration against curl of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Vh_i = \int_K f \nabla \times \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = \nabla \times \varphi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronHcurlLobatto::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    Seldon::Mlt(rh_loc, Vh, Uh);
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
  void HexahedronHcurlLobatto::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  void HexahedronHcurlLobatto::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        Uh(Node_to_DofX(i)) = Vh(3*i);
        Uh(Node_to_DofY(i)) = Vh(3*i+1);
        Uh(Node_to_DofZ(i)) = Vh(3*i+2);
      }
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
  void HexahedronHcurlLobatto::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        Vh(3*i) = Uh(Node_to_DofX(i));
        Vh(3*i+1) = Uh(Node_to_DofY(i));
        Vh(3*i+2) = Uh(Node_to_DofZ(i));
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
  void HexahedronHcurlLobatto
  ::ApplyShTransposeGen(int num_loc, const Vector1& res, Vector2& feval, int r) const
  {
    int ext = 0;
    if (num_loc > 2)
      ext = order;
    
    int num_dof;
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0:
      case 5 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(ext, i, j);
	      feval(3*NumQuad2D(i,j)) = res(num_dof); 
	      num_dof = NumDofs_Y(ext, i, j);
	      feval(3*NumQuad2D(i,j)+1) = res(num_dof);
	      num_dof = NumDofs_Z(ext, i, j);
	      feval(3*NumQuad2D(i,j)+2) = res(num_dof);
	    }
        break;
      case 1 :
      case 4 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(i, ext, j);
	      feval(3*NumQuad2D(i,j)) = res(num_dof);
	      num_dof = NumDofs_Y(i, ext, j);
	      feval(3*NumQuad2D(i,j)+1) = res(num_dof);
	      num_dof = NumDofs_Z(i, ext, j);
	      feval(3*NumQuad2D(i,j)+2) = res(num_dof);
	    }
        break;
      case 2 :
      case 3 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(i, j, ext);
	      feval(3*NumQuad2D(i,j)) = res(num_dof);
              num_dof = NumDofs_Y(i, j, ext);
	      feval(3*NumQuad2D(i,j)+1) = res(num_dof);
	      num_dof = NumDofs_Z(i, j, ext);
	      feval(3*NumQuad2D(i,j)+2) = res(num_dof);
	    }
        break;
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
  void HexahedronHcurlLobatto
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& feval, Vector2& res, int r) const
  {
    int ext = 0;
    if (num_loc > 2)
      ext = order;
    
    int num_dof;
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0:
      case 5 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(ext, i, j);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j));
	      num_dof = NumDofs_Y(ext, i, j);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j)+1);
	      num_dof = NumDofs_Z(ext, i, j);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j)+2);
	    }
        break;
      case 1 :
      case 4 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(i, ext, j);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j));
	      num_dof = NumDofs_Y(i, ext, j);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j)+1);
	      num_dof = NumDofs_Z(i, ext, j);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j)+2);
	    }
        break;
      case 2 :
      case 3 :
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(i, j, ext);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j));
              num_dof = NumDofs_Y(i, j, ext);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j)+1);
	      num_dof = NumDofs_Z(i, j, ext);
	      res(num_dof) += alpha*feval(3*NumQuad2D(i,j)+2);
	    }
        break;
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
  void HexahedronHcurlLobatto
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& res, Vector2& feval, int r) const
  {
    int num_dof;
    int ext = 0;
    typename Vector2::value_type vloc;
    if (num_loc > 2)
      ext = order;
    
    const Matrix<Real_wp>& dGL = lob_basis.GradPhi();
    feval.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)      
      {
      case 0 :
      case 5 :
        // face x = 0  or x = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(ext, i, j);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                {
                  feval(3*NumQuad2D(i, m)+1) += vloc*dGL(j, m);
                  feval(3*NumQuad2D(m, j)+2) -= vloc*dGL(i, m);
                }
              
	      num_dof = NumDofs_Y(ext, i, j);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                feval(3*NumQuad2D(i, m)) -= vloc*dGL(j, m);              
              
              vloc = 0;
              for (int k = 0; k <= order; k++)
                vloc += res(NumDofs_Y(k, i, j))*dGL(k, ext);
              
              feval(3*NumQuad2D(i, j)+2) += vloc;
              
	      num_dof = NumDofs_Z(ext, i, j);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                feval(3*NumQuad2D(m, j)) += vloc*dGL(i, m);
                 
              vloc = 0;
              for (int k = 0; k <= order; k++)
                vloc -= res(NumDofs_Z(k, i, j))*dGL(k, ext);

              feval(3*NumQuad2D(i, j)+1) += vloc;
	    }
        break;
      case 1 :
      case 4 :
        // face y = 0  or y = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_Y(i, ext, j);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                {
                  feval(3*NumQuad2D(m, j)+2) += vloc*dGL(i, m);
                  feval(3*NumQuad2D(i, m)) -= vloc*dGL(j, m);
                }

	      num_dof = NumDofs_Z(i, ext, j);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                feval(3*NumQuad2D(m, j)+1) -= vloc*dGL(i, m);
              
              vloc = 0;
              for (int k = 0; k <= order; k++)
                 vloc += res(NumDofs_Z(i, k, j))*dGL(k, ext);
              
              feval(3*NumQuad2D(i, j)) += vloc;
              
	      num_dof = NumDofs_X(i, ext, j);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                feval(3*NumQuad2D(i, m)+1) += vloc*dGL(j, m);
              
              vloc = 0;
              for (int k = 0; k <= order; k++)
                vloc -= res(NumDofs_X(i, k, j))*dGL(k, ext);
              
              feval(3*NumQuad2D(i, j)+2) += vloc;
	    }
        break;
      case 2 :
      case 3 :
        // face z = 0  or z = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_Z(i, j, ext);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                {
                  feval(3*NumQuad2D(i, m)) += vloc*dGL(j, m);
                  feval(3*NumQuad2D(m, j)+1) -= vloc*dGL(i, m);
                }
              
	      num_dof = NumDofs_X(i, j, ext);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                feval(3*NumQuad2D(i, m)+2) -= vloc*dGL(j, m);
              
              vloc = 0;
              for (int k = 0; k <= order; k++)
                vloc += res(NumDofs_X(i, j, k))*dGL(k, ext);
              
              feval(3*NumQuad2D(i, j)+1) += vloc;
              
	      num_dof = NumDofs_Y(i, j, ext);
              vloc = res(num_dof);
              for (int m = 0; m <= order; m++)
                feval(3*NumQuad2D(m, j)+2) += vloc*dGL(i, m);
              
              vloc = 0;
              for (int k = 0; k <= order; k++)
                vloc -= res(NumDofs_Y(i, j, k))*dGL(k, ext);
              
              feval(3*NumQuad2D(i, j)) += vloc;
	    }
        break;
      }
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
  void HexahedronHcurlLobatto
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& feval, Vector2& res, int r) const
  {
    int num_dof;
    int ext = 0;
    typename Vector2::value_type vloc;
    if (num_loc > 2)
      ext = order;
    
    const Matrix<Real_wp>& dGL = lob_basis.GradPhi();
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)      
      {
      case 0 :
      case 5 :
        // face x = 0  or x = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_X(ext, i, j);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                {
                  vloc += feval(3*NumQuad2D(i, m)+1)*dGL(j, m);
                  vloc -= feval(3*NumQuad2D(m, j)+2)*dGL(i, m);
                }
              res(num_dof) += alpha*vloc;
              
	      num_dof = NumDofs_Y(ext, i, j);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                vloc -= feval(3*NumQuad2D(i, m))*dGL(j, m);
              
              res(num_dof) += alpha*vloc;
              
              vloc = alpha*feval(3*NumQuad2D(i, j)+2);
              for (int k = 0; k <= order; k++)
                res(NumDofs_Y(k, i, j)) += vloc*dGL(k, ext);
              
	      num_dof = NumDofs_Z(ext, i, j);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                vloc += feval(3*NumQuad2D(m, j))*dGL(i, m);
              
              res(num_dof) += alpha*vloc;
              
              vloc = alpha*feval(3*NumQuad2D(i, j)+1);
              for (int k = 0; k <= order; k++)
                res(NumDofs_Z(k, i, j)) -= vloc*dGL(k, ext);
	    }
        break;
      case 1 :
      case 4 :
        // face y = 0  or y = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_Y(i, ext, j);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                {
                  vloc += feval(3*NumQuad2D(m, j)+2)*dGL(i, m);
                  vloc -= feval(3*NumQuad2D(i, m))*dGL(j, m);
                }
              res(num_dof) += alpha*vloc;
              
	      num_dof = NumDofs_Z(i, ext, j);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                vloc -= feval(3*NumQuad2D(m, j)+1)*dGL(i, m);
              
              res(num_dof) += alpha*vloc;
              
              vloc = alpha*feval(3*NumQuad2D(i, j));
              for (int k = 0; k <= order; k++)
                res(NumDofs_Z(i, k, j)) += vloc*dGL(k, ext);
              
	      num_dof = NumDofs_X(i, ext, j);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                vloc += feval(3*NumQuad2D(i, m)+1)*dGL(j, m);
              
              res(num_dof) += alpha*vloc;
              
              vloc = alpha*feval(3*NumQuad2D(i, j)+2);
              for (int k = 0; k <= order; k++)
                res(NumDofs_X(i, k, j)) -= vloc*dGL(k, ext);
	    }
        break;
      case 2 :
      case 3 :
        // face z = 0  or z = 1
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      num_dof = NumDofs_Z(i, j, ext);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                {
                  vloc += feval(3*NumQuad2D(i, m))*dGL(j, m);
                  vloc -= feval(3*NumQuad2D(m, j)+1)*dGL(i, m);
                }
              res(num_dof) += alpha*vloc;
              
	      num_dof = NumDofs_X(i, j, ext);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                vloc -= feval(3*NumQuad2D(i, m)+2)*dGL(j, m);
              
              res(num_dof) += alpha*vloc;
              
              vloc = alpha*feval(3*NumQuad2D(i, j)+1);
              for (int k = 0; k <= order; k++)
                res(NumDofs_X(i, j, k)) += vloc*dGL(k, ext);
              
	      num_dof = NumDofs_Y(i, j, ext);
              vloc = 0;
              for (int m = 0; m <= order; m++)
                vloc += feval(3*NumQuad2D(m, j)+2)*dGL(i, m);
              
              res(num_dof) += alpha*vloc;
              
              vloc = alpha*feval(3*NumQuad2D(i, j));
              for (int k = 0; k <= order; k++)
                res(NumDofs_Y(i, j, k)) -= vloc*dGL(k, ext);
	    }
        break;
      }
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
  void HexahedronHcurlLobatto::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int ix = NumDofs_X(i, j, k);
            int iy = NumDofs_Y(i, j, k);
            int iz = NumDofs_Z(i, j, k);
            
            // mass lumping
	    A.AddInteraction(m+ix, n+ix, C(0, 0)*this->weights_dof(ix));
	    A.AddInteraction(m+iy, n+iy, C(1, 1)*this->weights_dof(iy));
	    A.AddInteraction(m+iz, n+iz, C(2, 2)*this->weights_dof(iz));
                
	    A.AddInteraction(m+iy, n+ix, C(1, 0)*this->weights_dof(ix));
	    A.AddInteraction(m+iz, n+iy, C(2, 1)*this->weights_dof(iy));
	    A.AddInteraction(m+iz, n+ix, C(2, 0)*this->weights_dof(iz));
                
	    A.AddInteraction(m+ix, n+iy, C(0, 1)*this->weights_dof(ix));
	    A.AddInteraction(m+iy, n+iz, C(1, 2)*this->weights_dof(iy));
	    A.AddInteraction(m+ix, n+iz, C(0, 2)*this->weights_dof(iz));
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
  void HexahedronHcurlLobatto::
  AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const
  {
    T vloc;
    int ix, iy, iz, jx, jy, jz;
    const VectReal_wp& weights1d = this->Weights1D();
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            for (int kp = 0; kp <= order; kp++)
              {
                // diagonal part
                iz = NumDofs_Z(j, k, i);
                jz = NumDofs_Z(j, kp, i);
		A.AddInteraction(m+iz, n+jz, C(0, 0)*weights1d(i)*weights1d(j)*stiff1d(k, kp));
		
                iy = NumDofs_Y(j, i, k);
                jy = NumDofs_Y(j, i, kp);
		A.AddInteraction(m+iy, n+jy, C(0, 0)*weights1d(i)*weights1d(j)*stiff1d(k, kp));

                ix = NumDofs_X(i, j, k);
                jx = NumDofs_X(i, j, kp);
		A.AddInteraction(m+ix, n+jx, C(1, 1)*weights1d(i)*weights1d(j)*stiff1d(k, kp));
                
                iz = NumDofs_Z(k, j, i);
                jz = NumDofs_Z(kp, j, i);
		A.AddInteraction(m+iz, n+jz, C(1, 1)*weights1d(i)*weights1d(j)*stiff1d(k, kp));
                
                ix = NumDofs_X(i, k, j);
                jx = NumDofs_X(i, kp, j);
		A.AddInteraction(m+ix, n+jx, C(2, 2)*weights1d(i)*weights1d(j)*stiff1d(k, kp));
                
                iy = NumDofs_Y(k, i, j);
                jy = NumDofs_Y(kp, i, j);
		A.AddInteraction(m+iy, n+jy, C(2, 2)*weights1d(i)*weights1d(j)*stiff1d(k, kp));
                
                for (int ip = 0; ip <= order; ip++)
                  {
                    iz = NumDofs_Z(j, k, i);
                    jy = NumDofs_Y(j, ip, kp);
                    vloc = C(0, 0)*weights1d(i)*weights1d(ip)*weights1d(j)
                      *dGL_GL(k, ip)*dGL_GL(kp, i);
                    
		    A.AddInteraction(m+iz, n+jy, -vloc);
		    A.AddInteraction(m+jy, n+iz, -vloc);
                    
                    ix = NumDofs_X(i, j, k);
                    jz = NumDofs_Z(kp, j, ip);
                    vloc = C(1, 1)*weights1d(i)*weights1d(ip)*weights1d(j)
                      *dGL_GL(k, ip)*dGL_GL(kp, i);
                    
		    A.AddInteraction(m+ix, n+jz, -vloc);
		    A.AddInteraction(m+jz, n+ix, -vloc);
		    
                    ix = NumDofs_X(i, k, j);
                    jy = NumDofs_Y(kp, ip, j);
                    vloc = C(2, 2)*weights1d(i)*weights1d(ip)*weights1d(j)
                      *dGL_GL(k, ip)*dGL_GL(kp, i);
                    
		    A.AddInteraction(m+ix, n+jy, -vloc);
		    A.AddInteraction(m+jy, n+ix, -vloc);
                  }
              }
            
            // extra-diagonal part
            ix = NumDofs_X(i, j, k);
            iy = NumDofs_Y(j, i, k);
            iz = NumDofs_Z(j, k, i);
            for (int ip = 0; ip <= order; ip++)
              for (int jp = 0; jp <= order; jp++)
                for (int kp = 0; kp <= order; kp++)
                  {
                    jx = NumDofs_X(ip, jp, kp);
                    jy = NumDofs_Y(jp, ip, kp);
                    jz = NumDofs_Z(jp, kp, ip);
                    
                    // c10
                    if (ip == j)
                      {
                        vloc = C(0, 1)*weights1d(j)*dGL_GL(k, jp)*weights1d(jp)
                          *dGL_GL(kp, i)*weights1d(i);
                        
			A.AddInteraction(m+iz, n+jx, vloc);
			A.AddInteraction(m+jx, n+iz, vloc);
                      }
                    
                    if (i == ip)
                      {
                        vloc = -C(0, 1)*weights1d(i)*dGL_GL(jp, j)*weights1d(j)
                          *dGL_GL(k, kp)*weights1d(kp);
                        
			A.AddInteraction(m+iz, n+jz, vloc);                        
			A.AddInteraction(m+jz, n+iz, vloc);
                      }
                    
                    if ((ip == j) && (i == jp))
                      {
                        vloc = -C(0, 1)*weights1d(j)*weights1d(jp)*stiff1d(k, kp);
			A.AddInteraction(m+iy, n+jx, vloc);
			A.AddInteraction(m+jx, n+iy, vloc);
                      }
                    
                    if (i == kp)
                      {
                        vloc = C(0, 1)*dGL_GL(jp, j)*weights1d(j)*weights1d(kp)
                          *dGL_GL(k, ip)*weights1d(ip);
                        
			A.AddInteraction(m+iy, n+jz, vloc);
			A.AddInteraction(m+jz, n+iy, vloc);
                      }
                    
                    // c20
                    if (i == kp)
                      {
                        vloc = C(0, 2)*dGL_GL(jp, j)*weights1d(j)*dGL_GL(k, ip)
                          *weights1d(ip)*weights1d(kp);
                        
			A.AddInteraction(m+iz, n+jy, vloc);
			A.AddInteraction(m+jy, n+iz, vloc);
                      }
                    
                    if (i == ip)
                      {
                        vloc = -C(0, 2)*weights1d(i)*dGL_GL(jp, j)*weights1d(j)
                          *dGL_GL(k, kp)*weights1d(kp);
                        
			A.AddInteraction(m+iy, n+jy, vloc);
			A.AddInteraction(m+jy, n+iy, vloc);
                      }
                    
                    if ((i == kp)&&(ip == j))
                      {
                        vloc = -C(0, 2)*weights1d(j)*weights1d(kp)*stiff1d(jp, k);
			A.AddInteraction(m+iz, n+jx, vloc);
			A.AddInteraction(m+jx, n+iz, vloc);
                      }
                    
                    if (ip == j)
                      {
                        vloc = C(0, 2)*dGL_GL(k, kp)*weights1d(kp)*weights1d(j)
                          *dGL_GL(jp, i)*weights1d(i);
                        
			A.AddInteraction(m+iy, n+jx, vloc);
			A.AddInteraction(m+jx, n+iy, vloc);
                      }
                    
                    // c21
                    if (ip == j)
                      {
                        vloc = C(1, 2)*dGL_GL(k, kp)*weights1d(kp)*dGL_GL(jp, i)
                          *weights1d(i)*weights1d(j);
                        
			A.AddInteraction(m+ix, n+jy, vloc);
			A.AddInteraction(m+jy, n+ix, vloc);
                      }
                    
                    if (i == ip)
                      {
                        vloc = -C(1, 2)*weights1d(i)*dGL_GL(jp, j)*weights1d(j)
                          *dGL_GL(k, kp)*weights1d(kp);
                        
			A.AddInteraction(m+ix, n+jx, vloc);
			A.AddInteraction(m+jx, n+ix, vloc);
                      }
                    
                    if ((i == kp) && (k == ip))
                      {
                        vloc = -C(1, 2)*weights1d(k)*weights1d(kp)*stiff1d(jp, j);

			A.AddInteraction(m+iz, n+jy, vloc);
			A.AddInteraction(m+jy, n+iz, vloc);
                      }
                    
                    if (i == kp)
                      {
                        vloc = C(1, 2)*dGL_GL(jp, k)*weights1d(k)*weights1d(kp)
                          *dGL_GL(j, ip)*weights1d(ip);
                        
			A.AddInteraction(m+iz, n+jx, vloc);
			A.AddInteraction(m+jx, n+iz, vloc);
                      }
                  }
          }
    
  }


  //! Adds mass matrix for variable case
  template<class T, class Prop>
  void HexahedronHcurlLobatto
  ::AddVariableMassMatrixGen(int m, int n,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			     VirtualMatrix<T>& A) const
  {
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          {
            int ix = NumDofs_X(i, j, k);
            int iy = NumDofs_Y(i, j, k);
            int iz = NumDofs_Z(i, j, k);
            int p = NumQuad3D(i, j, k);
	    
            // mass lumping
	    A.AddInteraction(m+ix, n+ix, C(p)(0, 0));
	    A.AddInteraction(m+iy, n+iy, C(p)(1, 1));
	    A.AddInteraction(m+iz, n+iz, C(p)(2, 2));
	    
	    A.AddInteraction(m+iy, n+ix, C(p)(1, 0));
	    A.AddInteraction(m+iz, n+iy, C(p)(2, 1));
	    A.AddInteraction(m+iz, n+ix, C(p)(2, 0));
                
	    A.AddInteraction(m+ix, n+iy, C(p)(0, 1));
	    A.AddInteraction(m+iy, n+iz, C(p)(1, 2));
	    A.AddInteraction(m+ix, n+iz, C(p)(0, 2));
	  }
  }
  

  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point
  */
  void HexahedronHcurlLobatto::
  ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(R3());
    VectReal_wp phix(order+1), phiy(order+1), phiz(order+1);
    for (int i = 0; i <= order; i++)
      {
	phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiz(i) = lob_basis.EvaluatePhi(i, point_loc(2));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    Real_wp val = phix(i)*phiy(j)*phiz(k);
	    res(NumDofs_X(i,j,k))(0) = val;
	    res(NumDofs_Y(i,j,k))(1) = val;
	    res(NumDofs_Z(i,j,k))(2) = val;
	  }
  }

  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point
  */
  void HexahedronHcurlLobatto::
  ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(R3());
    VectReal_wp phix(order+1), phiy(order+1), phiz(order+1);
    VectReal_wp dphix(order+1), dphiy(order+1), dphiz(order+1);
    for (int i = 0; i <= order; i++)
      {
	phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	phiz(i) = lob_basis.EvaluatePhi(i, point_loc(2));

	dphix(i) = lob_basis.EvaluatePhiGrad(i, point_loc(0));
	dphiy(i) = lob_basis.EvaluatePhiGrad(i, point_loc(1));
	dphiz(i) = lob_basis.EvaluatePhiGrad(i, point_loc(2));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int nx = NumDofs_X(i, j, k);
	    res(nx).Init(0, dphiz(k)*phix(i)*phiy(j), -dphiy(j)*phix(i)*phiz(k));
	    
	    int ny = NumDofs_Y(i, j, k);
	    res(ny).Init(-dphiz(k)*phiy(j)*phix(i), 0, dphix(i)*phiy(j)*phiz(k));
	    
	    int nz = NumDofs_Z(i, j, k);
	    res(nz).Init(dphiy(j)*phiz(k)*phix(i), -dphix(i)*phiz(k)*phiy(j), 0);
	  }
  }
  

  //! computation of val_phi \f$ = \varphi_{node} (\xi_{k})  \f$
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions
    (xi_i) are Gauss-Lobatto quadrature points
  */  
  void HexahedronHcurlLobatto::GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const
  {
    int k1 = CoordinateQuad3D(k, 0);
    int k2 = CoordinateQuad3D(k, 1);
    int k3 = CoordinateQuad3D(k, 2);
    phi.Reallocate(nb_dof_loc); FillZero(phi);
    
    phi(NumDofs_X(k1, k2, k3))(0) = 1.0;
    phi(NumDofs_Y(k1, k2, k3))(1) = 1.0;
    phi(NumDofs_Z(k1, k2, k3))(2) = 1.0;
  }
  
  
  //! computation of curl_phi = \f$ \nabla \times \phi_{node} (\xi_{k}) \f$ 
  /*!
    \param[in] k quadrature point number
    \param[out] curl_phi curl of basis functions 
  */ 
  void HexahedronHcurlLobatto::GetCurlPhiOnQuadraturePoint(int k, VectR3& curl_phi) const
  {
    int k1 = CoordinateQuad3D(k, 0);
    int k2 = CoordinateQuad3D(k, 1);
    int k3 = CoordinateQuad3D(k, 2);
    curl_phi.Reallocate(nb_dof_loc); FillZero(curl_phi);
    
    for (int i = 0; i <= order; i++)
      {
        curl_phi(NumDofs_X(k1, k2, i))(1) = lob_basis.GradPhi(i, k3);
        curl_phi(NumDofs_X(k1, i, k3))(2) = -lob_basis.GradPhi(i, k2);
        
        curl_phi(NumDofs_Y(k1, k2, i))(0) = -lob_basis.GradPhi(i, k3);
        curl_phi(NumDofs_Y(i, k2, k3))(2) = lob_basis.GradPhi(i, k1);
        
        curl_phi(NumDofs_Z(k1, i, k3))(0) = lob_basis.GradPhi(i, k2);
        curl_phi(NumDofs_Z(i, k2, k3))(1) = -lob_basis.GradPhi(i, k1);
      }
  }
  
    
  //! H and E are discretized on Gauss-Lobatto points
  template<class Vector1, class Vector2>
  void HexahedronHcurlLobatto::
  ComputeNodalValuesH(const Vector1& Un, Vector2 & Unode)
  {
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
	Unode(i)(0) = Un(3*i);
	Unode(i)(1) = Un(3*i+1);
	Unode(i)(2) = Un(3*i+2);
      }
  }
  

  //! displays details about class HexahedronHcurlLobatto
  ostream& operator <<(ostream& out, const HexahedronHcurlLobatto& e)
  {
    out<<static_cast<const HexahedronReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_LOBATTO_CXX
#endif
