#ifndef MONTJOIE_FILE_TETRAHEDRON_SECOND_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronHcurlSecondFamily::TetrahedronHcurlSecondFamily() : TetrahedronReference<2>()
  {
    quadrature_equal_nodal = false;
    dof_equal_nodal = true;
    dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void TetrahedronHcurlSecondFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<2>::ConstructNumberMap(nmap, dg);
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order+1);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 0);
    // triangular face
    nmap.SetNbDofTriangle(order, 3*(order-1)*(order-2)/2 + 3*(order-1));
    
    // dofs inside the tetrahedron
    nmap.SetNbDofTetrahedron(order, nb_dof_loc-nb_dof_boundaries);
    nmap.SetNbDofHexahedron(order, 0);

    nmap.SetOppositeEdgesDofSymmetry(order, order+1);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
    
    // rotation of dofs on faces
    nmap.SetFacesDofRotationTri(order, FacesDof_Rotation_Tri);
    
    nmap.SetNbPointsQuadBoundariesTetrahedron(order, nb_points_quadrature_boundaries);
  }

  
  //! constructing finite element
  void TetrahedronHcurlSecondFamily::ConstructFiniteElement(int order_, int rsurf_tri, int rsurf_quad,
							    int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronReference<2>::ConstructFiniteElement(order_);
    ComputeCoefficientTransformation();
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
  }
  
  
  //! construction of basis functions
  void TetrahedronHcurlSecondFamily::ConstructFunctions()
  {
    int nb_dof_inside_tri = 3*(order-1)*(order-2)/2 + 3*(order-1);
    int nb_dof_elt_inside = (order-1)*(order-2)*(order-3)/2;
    nb_dof_tri = (order+1)*3 + nb_dof_inside_tri;
    nb_dof_boundaries = 4*nb_dof_inside_tri + 6*(order+1);
    nb_dof_loc = nb_dof_elt_inside + nb_dof_boundaries;
    nb_dof_quad = 0;
    
    FacesDof.Reallocate(nb_dof_tri, 4);
    NodalDof.Reallocate(nb_dof_loc);
    DirectionDof.Reallocate(nb_dof_loc);
    
    // dofs on the first edge
    nb_dof_loc = 0;
    for (int i = 0; i <= order; i++)
      {
	FacesDof(i,0) = nb_dof_loc;
	FacesDof(i,1) = nb_dof_loc;
	if (i == order)
	  DirectionDof(nb_dof_loc).Init(1.0,1.0,1.0);
	else if (i == 0)
	  DirectionDof(nb_dof_loc).Init(1.0,0.0,0.0);
	else
	  DirectionDof(nb_dof_loc).Init(1.0,0.5,0.5);
	NodalDof(nb_dof_loc) = NumNodes3D(i,0,0);
	nb_dof_loc++;
      }
    
    // second edge
    for (int i = 0; i <= order; i++)
      {
	FacesDof(2*(order+1)+i,0) = nb_dof_loc;
	FacesDof(i,2) = nb_dof_loc;
	if (i == order)
	  DirectionDof(nb_dof_loc).Init(1.0,1.0,1.0);
	else if (i == 0)
	  DirectionDof(nb_dof_loc).Init(0.0,1.0,0.0);
	else
	  DirectionDof(nb_dof_loc).Init(0.5,1.0,0.5);
	NodalDof(nb_dof_loc) = NumNodes3D(0,i,0);
	nb_dof_loc++;
      }
    
    // third edge
    for (int i = 0; i <= order; i++)
      {
	FacesDof(2*(order+1)+i,1) = nb_dof_loc;
	FacesDof(2*(order+1)+i,2) = nb_dof_loc;
	if (i == order)
	  DirectionDof(nb_dof_loc).Init(1.0,1.0,1.0);
	else if (i == 0)
	  DirectionDof(nb_dof_loc).Init(0.0,0.0,1.0);
	else
	  DirectionDof(nb_dof_loc).Init(0.5,0.5,1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,0,i);
	nb_dof_loc++;
      }
    
    // fourth edge
    for (int i = 0; i <= order; i++)
      {
	FacesDof(order+1+i,0) = nb_dof_loc;
	FacesDof(i,3) = nb_dof_loc;
	if (i == 0)
	  DirectionDof(nb_dof_loc).Init(0.0,1.0,0.0);
	else if (i == order)
	  DirectionDof(nb_dof_loc).Init(-1.0,0.0,0.0);
	else
	  DirectionDof(nb_dof_loc).Init(-0.5,0.5,0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(order-i,i,0);
	nb_dof_loc++;
      }
    
    // fifth edge
    for (int i = 0; i <= order; i++)
      {
	FacesDof(order+1+i,1) = nb_dof_loc;
	FacesDof(2*(order+1)+i,3) = nb_dof_loc;
	if (i == 0)
	  DirectionDof(nb_dof_loc).Init(0.0,0.0,1.0);
	else if (i == order)
	  DirectionDof(nb_dof_loc).Init(-1.0,0.0,0.0);
	else
	  DirectionDof(nb_dof_loc).Init(-0.5,0.0,0.5);
	NodalDof(nb_dof_loc) = NumNodes3D(order-i,0,i);
	nb_dof_loc++;
      }
    
    // sixth edge
    for (int i = 0; i <= order; i++)
      {
	FacesDof(order+1+i,2) = nb_dof_loc;
	FacesDof(order+1+i,3) = nb_dof_loc;
	if (i == 0)
	  DirectionDof(nb_dof_loc).Init(0.0,0.0,1.0);
	else if (i == order)
	  DirectionDof(nb_dof_loc).Init(0.0,-1.0,0.0);
	else
	  DirectionDof(nb_dof_loc).Init(0.0, -0.5, 0.5);
	NodalDof(nb_dof_loc) = NumNodes3D(0,order-i,i);
	nb_dof_loc++;
      }
    
    // edge-based face dofs
    int nb_dof_edges_face = 3*(order+1);
    // first face
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(0.0, 1.0, 0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(i,0,0);
	FacesDof(nb_dof_edges_face+i-1,0) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(-1.0, -1.0, -1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(order-i,i,0);
	FacesDof(nb_dof_edges_face+(order-1)+i-1,0) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(1.0, 0.0, 0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,i,0);
	FacesDof(nb_dof_edges_face+2*(order-1)+i-1,0) = nb_dof_loc;
	nb_dof_loc++;
      }
    // second face
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(0.0, 0.0, 1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(i,0,0);
	FacesDof(nb_dof_edges_face+i-1,1) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(-1.0, -1.0, -1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(order-i,0,i);
	FacesDof(nb_dof_edges_face+(order-1)+i-1,1) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(1.0, 0.0, 0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,0,i);
	FacesDof(nb_dof_edges_face+2*(order-1)+i-1,1) = nb_dof_loc;
	nb_dof_loc++;
      }
    // third face
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(0.0, 0.0, 1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,i,0);
	FacesDof(nb_dof_edges_face+i-1,2) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(-1.0, -1.0, -1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,order-i,i);
	FacesDof(nb_dof_edges_face+(order-1)+i-1,2) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(0.0, 1.0, 0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,0,i);
	FacesDof(nb_dof_edges_face+2*(order-1)+i-1,2) = nb_dof_loc;
	nb_dof_loc++;
      }
    // fourth face
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(0.0, 0.0, 1.0);
	NodalDof(nb_dof_loc) = NumNodes3D(order-i,i,0);
	FacesDof(nb_dof_edges_face+i-1,3) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(1.0, 0.0, 0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(0,order-i,i);
	FacesDof(nb_dof_edges_face+(order-1)+i-1,3) = nb_dof_loc;
	nb_dof_loc++;
      }
    for (int i = 1; i < order; i++)
      {
	DirectionDof(nb_dof_loc).Init(0.0, 1.0, 0.0);
	NodalDof(nb_dof_loc) = NumNodes3D(order-i,0,i);
	FacesDof(nb_dof_edges_face+2*(order-i)+i-1,3) = nb_dof_loc;
	nb_dof_loc++;
      }    

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    // first face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	if ((i+j) < order)
	  {
	  }
    
    // dofs inside the tetrahedron
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	for (int k = 1; k < order; k++)
	  if ((i+j+k) < order)
	    {
	      
	    }
    
    FacesDof_Rotation_Tri.Reallocate(5, nb_dof_inside_tri);
    SignFaceDof.Reallocate(6,nb_dof_inside_tri); SignFaceDof.Fill(true);
    for (int i = 1; i < order; i++)
      {
	FacesDof_Rotation_Tri(0, i-1) = 3*(order-1)-i;
	FacesDof_Rotation_Tri(1, i-1) = (order-1)+i-1;
	FacesDof_Rotation_Tri(2, i-1) = 2*(order-1)+i-1;
	FacesDof_Rotation_Tri(3, i-1) = (order-1)-i;
	FacesDof_Rotation_Tri(4, i-1) = 2*(order-1)-i;
      }
    
    for (int i = 1; i < order; i++)
      {
	int j = (order-1)+i-1;
	FacesDof_Rotation_Tri(0, j) = i-1;
	FacesDof_Rotation_Tri(1, j) = 3*(order-1)-i;
	FacesDof_Rotation_Tri(2, j) = 2*(order-1)-i;
	FacesDof_Rotation_Tri(3, j) = 2*(order-1)+i-1;
	FacesDof_Rotation_Tri(4, j) = (order-1)-i;
      }
      
    for (int i = 1; i < order; i++)
      {
	int j = 2*(order-1)+i-1;
	FacesDof_Rotation_Tri(0, j) = 2*(order-1)-i;
	FacesDof_Rotation_Tri(1, j) = (order-1)-i;
	FacesDof_Rotation_Tri(2, j) = i-1;
	FacesDof_Rotation_Tri(3, j) = (order-1)+i-1;
	FacesDof_Rotation_Tri(4, j) = 3*(order-1)-i;
      }
    
    ListeDof_Node.Reallocate(points_nodal3d.GetM(), 3);
    IVect nb_dof_node(points_nodal3d.GetM()); nb_dof_node.Fill(0);
    for (int i = 0; i < nb_dof_loc; i++)
      {
	int node = NodalDof(i);
	ListeDof_Node(node, nb_dof_node(node)) = i;
	nb_dof_node(node)++;
      }
    
    AhDof.Reallocate(points_nodal3d.GetM());
    Matrix3_3 A_tmp;
    for (int i = 0; i < points_nodal3d.GetM(); i++)
      {
	for (int j = 0; j < 3; j++)
	  for (int k = 0; k < 3; k++)
	    A_tmp(j,k) = DirectionDof(ListeDof_Node(i,k))(j);
	
	GetInverse(A_tmp, AhDof(i));
      }
    
    Value_Phi.Reallocate(nb_dof_loc, nb_points_quadrature_inside);
    VectR3 val_phi;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
	ComputeValuesPhiRef(points3d(j), val_phi);
	for (int i = 0; i < nb_dof_loc; i++)
	  Value_Phi(i, j) = val_phi(i);
      }
    
    Curl_Phi.Reallocate(nb_dof_loc, nb_points_quadrature_inside);
    VectR3 curl_phi;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
	ComputeCurlPhiRef(points3d(j), curl_phi);
	for (int i = 0; i < nb_dof_loc; i++)
	  Curl_Phi(i,j) = curl_phi(i);
      }
    
  }
  
  
  //! computation of elementary stiffness matrix on the reference element
  void TetrahedronHcurlSecondFamily::ConstructStiffnessMatrix()
  {
    stiffness_matrix00.Reallocate(nb_dof_loc, nb_dof_loc);
    stiffness_matrix11.Reallocate(nb_dof_loc, nb_dof_loc);
    stiffness_matrix22.Reallocate(nb_dof_loc, nb_dof_loc);
    stiffness_matrix01.Reallocate(nb_dof_loc, nb_dof_loc);
    stiffness_matrix02.Reallocate(nb_dof_loc, nb_dof_loc);
    stiffness_matrix12.Reallocate(nb_dof_loc, nb_dof_loc);
    Real_wp vloc1, vloc2, vloc3, vloc4, vloc5, vloc6, zero(0);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  vloc1 = zero; vloc2 = zero; vloc3 = zero;
	  vloc4 = zero; vloc5 = zero; vloc6 = zero;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    {
	      vloc1 += Curl_Phi(i,k)(0)*Curl_Phi(j,k)(0)*weights3d(k);
	      vloc2 += Curl_Phi(i,k)(1)*Curl_Phi(j,k)(1)*weights3d(k);
	      vloc3 += Curl_Phi(i,k)(2)*Curl_Phi(j,k)(2)*weights3d(k);
	      
	      vloc4 += Curl_Phi(i,k)(0)*Curl_Phi(j,k)(1)*weights3d(k);
	      vloc4 += Curl_Phi(i,k)(1)*Curl_Phi(j,k)(0)*weights3d(k);
	      vloc5 += Curl_Phi(i,k)(0)*Curl_Phi(j,k)(2)*weights3d(k);
	      vloc5 += Curl_Phi(i,k)(2)*Curl_Phi(j,k)(0)*weights3d(k);
	      vloc6 += Curl_Phi(i,k)(1)*Curl_Phi(j,k)(2)*weights3d(k);
	      vloc6 += Curl_Phi(i,k)(2)*Curl_Phi(j,k)(1)*weights3d(k);
	    }
	  stiffness_matrix00(i,j) = vloc1;
	  stiffness_matrix11(i,j) = vloc2;
	  stiffness_matrix22(i,j) = vloc3;
	  stiffness_matrix01(i,j) = vloc4;
	  stiffness_matrix02(i,j) = vloc5;
	  stiffness_matrix12(i,j) = vloc6;
	}
  }
  
  
  //! computation of elementary mass matrix on the reference element
  void TetrahedronHcurlSecondFamily::ConstructMassMatrix()
  {
    Real_wp zero(0);
    mass_matrix00.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix11.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix22.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix01.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix02.Reallocate(nb_dof_loc, nb_dof_loc);
    mass_matrix12.Reallocate(nb_dof_loc, nb_dof_loc);
    Real_wp vloc1, vloc2, vloc3, vloc4, vloc5, vloc6;
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  vloc1 = zero; vloc2 = zero; vloc3 = zero;
	  vloc4 = zero; vloc5 = zero; vloc6 = zero;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    {
	      vloc1 += Value_Phi(i,k)(0)*Value_Phi(j,k)(0)*weights3d(k);
	      vloc2 += Value_Phi(i,k)(1)*Value_Phi(j,k)(1)*weights3d(k);
	      vloc3 += Value_Phi(i,k)(2)*Value_Phi(j,k)(2)*weights3d(k);
	      
	      vloc4 += Value_Phi(i,k)(0)*Value_Phi(j,k)(1)*weights3d(k);
	      vloc4 += Value_Phi(i,k)(1)*Value_Phi(j,k)(0)*weights3d(k);
	      vloc5 += Value_Phi(i,k)(0)*Value_Phi(j,k)(2)*weights3d(k);
	      vloc5 += Value_Phi(i,k)(2)*Value_Phi(j,k)(0)*weights3d(k);
	      vloc6 += Value_Phi(i,k)(1)*Value_Phi(j,k)(2)*weights3d(k);
	      vloc6 += Value_Phi(i,k)(2)*Value_Phi(j,k)(1)*weights3d(k);
	    }
	  mass_matrix00(i,j) = vloc1;
	  mass_matrix11(i,j) = vloc2;
	  mass_matrix22(i,j) = vloc3;
	  mass_matrix01(i,j) = vloc4;
	  mass_matrix02(i,j) = vloc5;
	  mass_matrix12(i,j) = vloc6;
	}
  }
  
  
  //! returns stiffness interaction \f$ \int coef curl(\varphi_i) curl(\varphi_j) \f$
  Real_wp TetrahedronHcurlSecondFamily::
  ComputeStiffnessInteraction(int i, int j, const Matrix3_3sym& coef) const
  {
    Real_wp vloc;
    // stiffness part
    vloc = stiffness_matrix00(i,j) * coef(0,0) + stiffness_matrix11(i,j) * coef(1,1) + 
      stiffness_matrix01(i,j) * coef(0,1) + stiffness_matrix22(i,j) * coef(2,2) + 
      stiffness_matrix02(i,j) * coef(0,2) + stiffness_matrix12(i,j) * coef(1,2);
    
    return vloc;
  }
  
  
  //! returns mass interaction \f$ \int coef \varphi_i \varphi_j \f$
  Real_wp TetrahedronHcurlSecondFamily::
  ComputeMassInteraction(int i, int j, const Matrix3_3sym& coef) const
  {
    Real_wp vloc;
    // mass part
    vloc = mass_matrix00(i,j) * coef(0,0) + mass_matrix11(i,j) * coef(1,1) + 
      mass_matrix01(i,j) * coef(0,1) + mass_matrix22(i,j) * coef(2,2) + 
      mass_matrix02(i,j) * coef(0,2) + mass_matrix12(i,j) * coef(1,2);
    
    return vloc;
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TetrahedronHcurlSecondFamily::
  ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc); VectReal_wp phi;
    ComputeValuesPhiNodalRef(point_loc, phi);
    for (int node = 0; node < nb_dof_loc; node++)
      {
	res(node) = DirectionDof(node);
	Mlt(phi(node), res(node));
      }
  }

  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void TetrahedronHcurlSecondFamily::ComputeCurlPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc); VectR3 grad_phi;
    ComputeGradientPhiNodalRef(point_loc, grad_phi);
    for (int node = 0; node < nb_dof_loc; node++)
      TimesProd(grad_phi(node), DirectionDof(node), res(node));
  }
  
  
    //! returns value of \f$ \varphi_{node} (\xi_{num\_point}) \f$ 
  const R3& TetrahedronHcurlSecondFamily::GetValuePhi(int node, int num_point) const
  {
    return Value_Phi(node, num_point);
  }
  
   
  //! returns n \times phi
  R3 TetrahedronHcurlSecondFamily::
  GetValuePhiOnBoundary(int num_dof, int k, int num_loc) const
  {
    R3 phi = Value_Phi(num_dof, this->num_quad_points_surf(num_loc)(k));
    switch (num_loc)
      {
      case 0 :
	return R3(phi(1), -phi(0), 0);
      case 1 :
	return R3(-phi(2), 0, phi(0));
      case 2 :
	return R3(0, phi(2), -phi(1));
      case 3 :
	return R3(phi(2)-phi(1), phi(0)-phi(2), phi(1)-phi(0));
      }
    return R3();
  }
  
    
  //! computes \f$ curl(\varphi_{node}) (\xi_{num\_point}) \f$
  const R3& TetrahedronHcurlSecondFamily::GetCurlPhi(int node, int num_point) const
  {
    return Curl_Phi(node, num_point);
  }
  
   
  //! returns n \times  curl(phi)
  R3 TetrahedronHcurlSecondFamily::
  GetCurlPhiOnBoundary(int num_dof, int k, int num_loc) const
  {
    R3 phi = Curl_Phi(num_dof, this->num_quad_points_surf(num_loc)(k));
    switch (num_loc)
      {
      case 0 :
	return R3(phi(1), -phi(0), 0);
      case 1 :
	return R3(-phi(2), 0, phi(0));
      case 2 :
	return R3(0, phi(2), -phi(1));
      case 3 :
	return R3(phi(2)-phi(1), phi(0)-phi(2), phi(1)-phi(0));
      }
    
    return phi;
  }
  
  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are placed
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeProjectionDofRef(const Vector1& feval, Vector2& contrib) const 
  {
    typedef typename Vector2::value_type Complexe;
    TinyVector<Complexe,3> vec_u, projection_u;
    for (int i = 0; i < points_nodal3d.GetM(); i++)
      {
	for (int k = 0; k < 3; k++)
	  vec_u(k) = feval(i)(k);
	
	Mlt(AhDof(i), vec_u, projection_u);
	for (int j = 0; j < 3; j++)
	  contrib(ListeDof_Node(i,j)) = projection_u(j);
      }
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \varphi_i \f$ for all i
  /*!
    \param[in] feval values of f on quadrature points
    \param[out] res result vector
   */
  template<class Vector1, class Vector2>
  void TetrahedronHcurlSecondFamily::ComputeIntegralRef(const Vector1& feval,
							  Vector2& res) const
  {
    res.Fill(0);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int k = 0; k < nb_points_quadrature_inside; k++)
	res(i) += DotProd(feval(k), Value_Phi(i,k));
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f curl(\varphi_i) \f$ for all i
  /*!
    \param[in] feval values of f on quadrature points
    \param[out] res result vector
   */
  template<class Vector1, class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeIntegralCurlRef(const Vector1& feval, Vector2& res) const
  {
    res.Fill(0);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int k = 0; k < nb_points_quadrature_inside; k++)
	res(i) += DotProd(feval(k), Curl_Phi(i,k));
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \varphi_i \times n \f$ for all i
  /*!
    \param[in] feval values of f on quadrature points of the face
    \param[out] res result vector (dof components)
    \param[in] num_loc local face number
   */
  template<class Vector1, class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeIntegralSurfaceRef(const Vector1& feval, Vector2& res, int num_loc) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    for (int i = 0; i < FacesDof.GetM(); i++)
      {
	int num_dof = FacesDof(i,num_loc);
	for (int k = 0; k < points2d_tri.GetM(); k++)
	  res(num_dof) += DotProd(feval(k), GetValuePhiOnBoundary(num_dof, k, num_loc));
      }
  }
  
  
  template<class Vector1, class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeGaussIntegralSurfaceRef(const Vector1& feval, Vector2& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  template<class Vector1, class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeIntegralSurfaceCurlRef(const Vector1& feval, Vector2& res, int num_loc) const
  {
    abort();
  }
  
  
  //! computation of U on nodal points (on reference element)
  /*!
    \param[in] u_loc components of U on dofs
    \param[out] u_boundary values of U on nodal points
    \param[in] num_loc local face number
  */
  template<class Vector1,class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeValueBoundaryRef(const Vector1& u_loc, Vector2& u_boundary, int num_loc) const
  {
    abort();
  }
  
  
  //! not implemented
  template<class Vector1,class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeCurlBoundaryRef(const Vector1& u_loc, Vector2& u_boundary,
			 int num_loc) const
  {
    abort();
  }
  
  
  //! not implemented
  template<class Vector1,class Vector2>
  void TetrahedronHcurlSecondFamily::ComputeNodalValuesRef(const Vector1& u_loc,
						      Vector2& u_boundary) const
  {
    abort();
  }
  
  
  //! not implemented
  template<class Vector1,class Vector2>
  void TetrahedronHcurlSecondFamily::
  ComputeNodalCurlRef(const Vector1& u_loc, Vector2& u_boundary) const
  {
    abort();
  }

  
  //! displays details of class TetrahedronHcurlSecondFamily
  ostream& operator <<(ostream& out, const TetrahedronHcurlSecondFamily& e)
  {
    out<<static_cast<const TetrahedronReference<2>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_SECOND_FAMILY_CXX
#endif
