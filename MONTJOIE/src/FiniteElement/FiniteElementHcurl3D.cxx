#ifndef MONTJOIE_FILE_FINITE_ELEMENT_HCURL3D_CXX

namespace Montjoie
{
  
  //! compute \f$ \varphi_i(\mbox{pointloc}) \f$
  /*!
    \param[in] point_loc local coordinates of the point on the unit element
    \param[out] res values of basis functions at local point
    \param[in] dfjm1 inverse of the jacobian matrix
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
  */
  void FiniteElementHcurl<Dimension3>::
  ComputeValuesPhi(const R3& point_loc, Vector<R3>& res, const Matrix3_3& dfjm1,
                   const MeshNumbering<Dimension3>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
    
    res.Reallocate(FaceBasis.GetNbDof());
    FaceBasis.ComputeValuesPhiRef(point_loc, res);
    R3 phi;
    if (FaceBasis.UsePiolaTransform())
      for (int i = 0; i < FaceBasis.GetNbDof(); i++)
	{
	  phi = res(i);
	  MltTrans(dfjm1, phi, res(i));
	}
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \nabla \times \varphi_i (\mbox{pointloc}) \f$
  /*!
    \param[in] point_loc local coordinates of the point on the unit element
    \param[out] res curl of basis functions at local point
    \param[in] dfjm1 inverse of the jacobian matrix
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
  */
  void FiniteElementHcurl<Dimension3>::
  ComputeValuesGradientPhi(const R3& point_loc, VectR3& res, const Matrix3_3& dfjm1,
			   const MeshNumbering<Dimension3>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    res.Reallocate(nb_dof_loc); R3 curl, vec;
    if (FaceBasis.UsePiolaTransform())
      {
	Matrix3_3 dfj; GetInverse(dfjm1, dfj); Real_wp invJacob = Det(dfjm1);
	FaceBasis.ComputeCurlPhiRef(point_loc, res);
	for (int node = 0; node < nb_dof_loc; node++)
	  {
	    curl = res(node);
	    Mlt(dfj, curl, res(node));
	    Mlt(invJacob, res(node));
	  }
      }
    else
      {
	cout << "not implemented" << endl;
	abort();
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! Computes curl U from the gradient of U
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  GetCurlFromGradient(const Vector<Vector<T> >& dU, Vector<Vector<T> >& gradU)
  {
    int nb_unknowns = dU.GetM() / 9;
    gradU.Reallocate(3*nb_unknowns);
    for (int m = 0; m < nb_unknowns; m++)
      {
	int nb_nodes = dU(0).GetM();
	gradU(3*m).Reallocate(nb_nodes);
	gradU(3*m+1).Reallocate(nb_nodes);
	gradU(3*m+2).Reallocate(nb_nodes);
        for (int j = 0; j < nb_nodes; j++)
          {
            // curl U = (dUz/dy - dUy/dz, dUx/dz - dUz/dx, dUy/dx - dUx/dy)
            gradU(3*m)(j) = dU(9*m+7)(j) - dU(9*m+5)(j);
            gradU(3*m+1)(j) = dU(9*m+2)(j) - dU(9*m+6)(j);
            gradU(3*m+2)(j) = dU(9*m+3)(j) - dU(9*m+1)(j);
          }
      }
  }
  
  
  //! computation of basis functions on quadrature points on the edge
  /*!
    \param[out] ValuePhi_Boundary values
    \param[in] MatricesElem jacobian matrices on the element
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
    \param[in] num_loc local face number
   */
  // computes ValuePhi_Boundary(i,k) = \phi_i(\xi_k),
  // where \xi_k is a quadrature point on the face num_loc
  void FiniteElementHcurl<Dimension3>::
  ComputeValuesPhiBoundary(Matrix<R3>& ValuePhi_Boundary,
                           const SetMatrices<Dimension3>& MatricesElem,
			   const MeshNumbering<Dimension3>& mesh_num, int nquad, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    ValuePhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    FillZero(ValuePhi_Boundary);
    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }
    
    Vector<R3> val_phi_chap(nb_dof_loc);
    R3 val_phi; Matrix3_3 dfjm1;
    for (int k = 0; k < nb_quad; k++)
      {        
        GetInverse(MatricesElem.GetPointQuadratureBoundary(k), dfjm1);
        
        int num_point = FaceBasis.GetQuadNumber(num_loc, k);
        FaceBasis.GetValuePhiOnQuadraturePoint(num_point, val_phi_chap);
        for (int num_dof = 0; num_dof < nb_dof_loc; num_dof++)
          {
            MltTrans(dfjm1, val_phi_chap(num_dof), val_phi);            
	    ValuePhi_Boundary(num_dof, k) = val_phi; 
	  }
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, ValuePhi_Boundary, nquad);
  }
  
  
  //! computation of curl of basis functions on quadrature points on the edge/face
  /*!
    \param[out] GradientPhi_Boundary values of derivative
    \param[in] MatricesElem jacobian matrices on the element
    \param[in] num_loc local edge/face number
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
    \param[in] num_loc local face number
   */
  void FiniteElementHcurl<Dimension3>::
  ComputeValuesGradientPhiBoundary(Matrix<R3>& CurlPhi_Boundary,
				   const SetMatrices<Dimension3>& MatricesElem, 
				   const MeshNumbering<Dimension3>& mesh_num,
				   int nquad, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    CurlPhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    FillZero(CurlPhi_Boundary);
    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }
    
    VectR3 curl_phi_chap(nb_dof_loc);
    R3 curl_phi; Real_wp jacob;

    for (int k = 0; k < nb_quad; k++)
      {
        jacob = 1.0/Det(MatricesElem.GetPointQuadratureBoundary(k));
        int num_point = FaceBasis.GetQuadNumber(num_loc, k);
        FaceBasis.GetCurlPhiOnQuadraturePoint(num_point, curl_phi_chap);
        
        for (int j = 0; j < nb_dof_loc; j++)
          {
	    Mlt(MatricesElem.GetPointQuadratureBoundary(k), curl_phi_chap(j), curl_phi);
            curl_phi *= jacob;
            
	    CurlPhi_Boundary(j, k) = curl_phi; 
	  }
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, CurlPhi_Boundary, nquad);
  }
  
  
  //! computation of U on nodal points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
    \param[in] mesh considered mesh
    \param[in] num_elem element number in the mesh
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeNodalValues(const SetMatrices<Dimension3>& MatricesElem,
		     const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
		     const Mesh<Dimension3>& mesh, int num_elem) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int nb_points_nodal = FaceBasis.GetNbPointsNodalElt();
    int nb_unknowns = Uloc.GetM();
    Uloc_node.Reallocate(3*nb_unknowns);
    for (int m = 0; m < 3*nb_unknowns; m++)
      Uloc_node(m).Reallocate(nb_points_nodal);
    
    Matrix3_3 dfjm1; TinyVector<T, 3> vec_u, vec_v;
    Vector<T> Unode(3*nb_points_nodal);
    
    // loop on all the unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// computation on reference element
	FaceBasis.ComputeNodalValuesRef(Uloc(m), Unode);
	
	for (int i = 0; i < nb_points_nodal; i++)
	  {
	    Uloc_node(3*m)(i) = Unode(3*i);
	    Uloc_node(3*m+1)(i) = Unode(3*i+1);
	    Uloc_node(3*m+2)(i) = Unode(3*i+2);
	  }
      }
    
    // application of transformation DF_i^{*-1}
    if (FaceBasis.UsePiolaTransform())
      for (int i = 0; i < nb_points_nodal; i++)
	{
	  GetInverse(MatricesElem.GetPointNodal(i), dfjm1);
	  
	  for (int m = 0; m < nb_unknowns; m++)
	    {
	      vec_u(0) = Uloc_node(3*m)(i);
	      vec_u(1) = Uloc_node(3*m+1)(i);
	      vec_u(2) = Uloc_node(3*m+2)(i);
	      
	      MltTrans(dfjm1, vec_u, vec_v);
	      
	      Uloc_node(3*m)(i) = vec_v(0);
	      Uloc_node(3*m+1)(i) = vec_v(1);
	      Uloc_node(3*m+2)(i) = vec_v(2);            
	    }
	}
  }
  

  //! computation of U on quadrature points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on quadrature points
    \param[out] dUloc_node curl of U on quadrature points
    \param[in] compute_u if true, value of u is computed
    \param[in] compute_curl if true, curl of u is computed
    \param[in] mesh considered mesh
    \param[in] iquad element number in the mesh
  */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeQuadratureValues(const SetMatrices<Dimension3>& MatricesElem,
			  const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
			  Vector<Vector<T> >& dUloc_node, bool compute_u, bool compute_curl,
			  const Mesh<Dimension3>& mesh, int iquad) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int nb_points_quad = FaceBasis.GetNbPointsQuadratureInside();
    int nb_unknowns = Uloc.GetM();
    if (compute_u)
      {
	Uloc_node.Reallocate(3*nb_unknowns);
	for (int m = 0; m < 3*nb_unknowns; m++)
	  Uloc_node(m).Reallocate(nb_points_quad);
      }
    
    if (compute_curl)
      {
	dUloc_node.Reallocate(3*nb_unknowns);
	for (int m = 0; m < 3*nb_unknowns; m++)
	  dUloc_node(m).Reallocate(nb_points_quad);
      }

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }
    
    Vector<T> Uquad(3*nb_points_quad), Ucurl(3*nb_points_quad);
    TinyVector<T, 3> val, vloc; Matrix3_3 dfjm1, invJacob_DF;
    Real_wp invJacob;
    
    // loop on all unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// projection on quadrature points
	FaceBasis.ApplyChTranspose(Uloc(m), Uquad);
	if (compute_curl)
	  {
	    if (FaceBasis.UseQuadraturePointsForRh())
              FaceBasis.ApplyRhQuadratureTranspose(Uquad, Ucurl);
	    else
	      FaceBasis.ApplyRhTranspose(Uloc(m), Ucurl);
	  }
	
	// storing values on Uloc_node
	if (compute_u)	 
	  for (int i = 0; i < nb_points_quad; i++)
	    {
	      Uloc_node(3*m)(i) = Uquad(3*i);
	      Uloc_node(3*m+1)(i) = Uquad(3*i+1);
	      Uloc_node(3*m+2)(i) = Uquad(3*i+2);
	    }
        
	if (compute_curl)
	  for (int i = 0; i < nb_points_quad; i++)
	    {
	      dUloc_node(3*m)(i) = Ucurl(3*i);
	      dUloc_node(3*m+1)(i) = Ucurl(3*i+1);
	      dUloc_node(3*m+2)(i) = Ucurl(3*i+2);
	    }
      }
    
    // applying transformation DF_i^{*-1}
    for (int i = 0; i < nb_points_quad; i++)
      {
	GetInverse(MatricesElem.GetPointQuadrature(i), dfjm1);
	invJacob = Det(dfjm1);
	invJacob_DF = MatricesElem.GetPointQuadrature(i);
	Mlt(invJacob, invJacob_DF);
        
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    if (compute_u)
	      {
		val(0) = Uloc_node(3*m)(i);
		val(1) = Uloc_node(3*m+1)(i);
		val(2) = Uloc_node(3*m+2)(i);
		
		MltTrans(dfjm1, val, vloc);
		
		Uloc_node(3*m)(i) = vloc(0);
		Uloc_node(3*m+1)(i) = vloc(1);
		Uloc_node(3*m+2)(i) = vloc(2);
	      }
	    
	    if (compute_curl)
	      {
		val(0) = dUloc_node(3*m)(i);
		val(1) = dUloc_node(3*m+1)(i);
		val(2) = dUloc_node(3*m+2)(i);
		
		Mlt(invJacob_DF, val, vloc);
		
		dUloc_node(3*m)(i) = vloc(0);
		dUloc_node(3*m+1)(i) = vloc(1);
		dUloc_node(3*m+2)(i) = vloc(2);
	      }
	  }
      }
	    
  }

  
  //! project function f on degrees of freedom
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on dof points
    \param[out] res components on degrees of freedom
    \param[in] mesh_num mesh numbering
    \param[in] n element number in the mesh
  */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeProjectionDof(const SetMatrices<Dimension3>& MatricesElem,
		       Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
		       const MeshNumbering<Dimension3>& mesh_num, int n) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    // we apply DF_i^t to the evaluation of function f
    int nb_dof_loc = FaceBasis.GetNbDof(), nb_points_dof = FaceBasis.GetNbPointsDofInside();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    TinyVector<T, 3> E1, E2;
    for (int m = 0; m < nb_unknowns; m++)
      {
	for (int k = 0; k < nb_points_dof; k++)
	  {
            E1(0) = feval(m)(3*k);
            E1(1) = feval(m)(3*k+1);
            E1(2) = feval(m)(3*k+2);
	    
	    MltTrans(MatricesElem.GetPointDof(k), E1, E2);
	    
            feval(m)(3*k) = E2(0);
            feval(m)(3*k+1) = E2(1);
            feval(m)(3*k+2) = E2(2);
	  }
	
	// then we project on reference element
	res(m).Reallocate(nb_dof_loc);
	FaceBasis.ComputeProjectionDofRef(feval(m), res(m));        
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.GetGlobalUnknownVector(mesh_num, res, n);
  }


  //! project function f on degrees of freedom
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on dof points
    \param[out] res components on degrees of freedom
    \param[in] mesh_num mesh numbering
    \param[in] n element number in the mesh
    \param[in] num_loc local face number
  */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeProjectionSurfaceDof(const SetMatrices<Dimension3>& MatricesElem,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			      const MeshNumbering<Dimension3>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
    
    // we apply DF_i^t to the evaluation of function f
    int nb_dof = FaceBasis.GetNbDofBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    Vector<Vector<T> > u_loc(1);
    u_loc(0).Reallocate(FaceBasis.GetNbDof());
    u_loc(0).Zero();
    res.Reallocate(nb_unknowns);

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    if (FaceBasis.DiscontinuousElement())
      {
	cout << "This method is not compatible with discontinuous element" << endl;
	abort();
      }

    TinyVector<T, 3> E1, E2;
    for (int m = 0; m < nb_unknowns; m++)
      {
	for (int k = 0; k < FaceBasis.GetNbPointsDofSurface(num_loc); k++)
	  {
            E1(0) = feval(m)(3*k);
            E1(1) = feval(m)(3*k+1);
            E1(2) = feval(m)(3*k+2);
	    
	    MltTrans(MatricesElem.GetPointDofBoundary(k), E1, E2);
	    
            feval(m)(3*k) = E2(0);
            feval(m)(3*k+1) = E2(1);
            feval(m)(3*k+2) = E2(2);
	  }
	
	// then we project on reference element
	res(m).Reallocate(nb_dof);
	FaceBasis.ComputeProjectionSurfaceDofRef(feval(m), res(m), num_loc);        

	// using a vector with all degrees of freedom
	// in order to take into account signs
	for (int j = 0; j < nb_dof; j++)
	  u_loc(0)(FaceBasis.GetLocalNumber(num_loc, j)) = res(m)(j);

	mesh_num.number_map.GetGlobalUnknownVector(mesh_num, u_loc, n);
	
	for (int j = 0; j < nb_dof; j++)
	  res(m)(j) = u_loc(0)(FaceBasis.GetLocalNumber(num_loc, j));
      }
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \varphi_i \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral agains each basis function
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number in the mesh
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeIntegral(const SetMatrices<Dimension3>& mat,
		  Vector<Vector<T> >& feval, Vector<Vector<T> > & res,
		  const MeshNumbering<Dimension3>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    int Nquad = FaceBasis.GetNbPointsQuadratureInside();
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    TinyVector<T, 3> E1, E2; 
    if (!FaceBasis.UsePiolaTransform())
      {
	for (int k = 0; k < Nquad; k++)
	  {
	    Real_wp coef = FaceBasis.WeightsND(k)*Det(mat.GetPointQuadrature(k));
	    
	    for (int m = 0; m < nb_unknowns; m++)
	      {
		CopyVector(feval(m), k, E1);
		E1 *= coef;
		CopyVector(E1, k, feval(m));
	      }
	  }
      }
    else
      {
	Matrix3_3 dfjm1;
	
	// we apply J_i DF_i^{-t} \omega to the functions
	for (int k = 0; k < Nquad; k++)
	  {
	    GetInverse(mat.GetPointQuadrature(k), dfjm1);
	    Real_wp coef = FaceBasis.WeightsND(k)*Det(mat.GetPointQuadrature(k));
	    
	    for (int m = 0; m < nb_unknowns; m++)
	      {
		E1(0) = feval(m)(3*k);
		E1(1) = feval(m)(3*k+1);
		E1(2) = feval(m)(3*k+2);
		
		Mlt(dfjm1, E1, E2);	    
		Mlt(coef, E2);
		
		feval(m)(3*k) = E2(0);
		feval(m)(3*k+1) = E2(1);
		feval(m)(3*k+2) = E2(2);
	      }
	  }
      }
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
	// then we integrate on the reference element
	FaceBasis.ComputeIntegralRef(feval(m), res(m));
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \cdot \nabla \times \varphi \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral against curl of each basis function
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number in the mesh
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeIntegralGradient(const SetMatrices<Dimension3>& mat,
			  Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
			  const MeshNumbering<Dimension3>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int Nquad = FaceBasis.GetNbPointsQuadratureInside();
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    TinyVector<T, 3> E1, E2;
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we apply DF_i^t \omega to the functions
	for (int i = 0; i < Nquad; i++)
	  {
            E1(0) = feval(m)(3*i);
            E1(1) = feval(m)(3*i+1);
            E1(2) = feval(m)(3*i+2);
	    
	    MltTrans(mat.GetPointQuadrature(i), E1, E2);
	    
	    Real_wp coef = FaceBasis.WeightsND(i);
	    feval(m)(3*i) = coef*E2(0);
            feval(m)(3*i+1) = coef*E2(1);
            feval(m)(3*i+2) = coef*E2(2);
	  }
	
	res(m).Reallocate(nb_dof_loc);
	// then we integrate on reference element
	FaceBasis.ComputeIntegralCurlRef(feval(m), res(m));
      }

    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i \times n ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains each basis function
    \param[in] mesh_num mesh numbering
    \param[in] n element number in the mesh
    \param[in] num_loc local edge number
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeIntegralSurface(const SetMatrices<Dimension3>& mat,
			 Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
			 const MeshNumbering<Dimension3>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_points_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    // we apply DF_i^-1 on evaluation of f
    TinyVector<T, 3> E1, E2;
    Matrix3_3 dfjm1;
    for (int i = 0; i < nb_points_quad; i++)
      {
	GetInverse(mat.GetPointQuadratureBoundary(i), dfjm1);
	Real_wp coef = FaceBasis.WeightsQuadratureBoundary(i, num_loc);
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    E1(0) = feval(m)(3*i);
            E1(1) = feval(m)(3*i+1);
            E1(2) = feval(m)(3*i+2);
	    
	    Mlt(dfjm1, E1, E2);
	    Mlt(mat.GetDsQuadratureBoundary(i)*coef, E2);
            
	    feval(m)(3*i) = E2(0); 
            feval(m)(3*i+1) = E2(1); 
            feval(m)(3*i+2) = E2(2); 
	  }
      }
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
        FaceBasis.ComputeIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the face
    \param[out] res integral agains each basis function
    \param[in] mesh_num mesh numbering
    \param[in] n element number in the mesh
    \param[in] num_loc local edge number
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeGaussIntegralSurface(const Vector<Matrix3_3>& dfjm1,
                              const VectReal_wp& weights, const VectReal_wp& ds,
                              Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
                              const MeshNumbering<Dimension3>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
        
    int nb_points_quad = weights.GetM();
    int nb_unknowns = feval.GetM();

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    // we apply DF_i^-1 on evaluation of f
    TinyVector<T, 3> E1, E2;
    res.Reallocate(nb_unknowns);
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(FaceBasis.GetNbDof());
	for (int i = 0; i < nb_points_quad; i++)
	  {
            E1(0) = feval(m)(3*i)*weights(i);
            E1(1) = feval(m)(3*i+1)*weights(i);
            E1(2) = feval(m)(3*i+2)*weights(i);
	    
	    Mlt(dfjm1(i), E1, E2);
	    
	    feval(m)(3*i) = E2(0);
            feval(m)(3*i+1) = E2(1);
            feval(m)(3*i+2) = E2(2);
	  }
	
        FaceBasis.ComputeGaussIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }
  

  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains each basis function
    \param[in] num_loc local face number
    \param[in] mesh_num mesh numbering
    \param[in] n element number
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeIntegralSurfaceHDG(const SetMatrices<Dimension3>& mat,
			    Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
			    const MeshNumbering<Dimension3>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& Fb_vol
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    const ElementReferenceTrace& Fb_surf = Fb_vol.GetSurfaceFiniteElement(num_loc);
    
    // initializations
    int nb_quad = Fb_surf.GetNbPointsQuadratureInside();
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);

    int nb_dof = Fb_surf.GetNbDof();
    
    // loop on each unknown
    Matrix3_3 dfjm1; TinyVector<T, 3> vec_u, vec_v; TinyVector<T, 2> vec_vs;
    Vector<T> feval_surf(nb_quad*2);
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we multiply the function f by the weight of integration (ds \omega_i)
        // and DF_i^{-1}
	for (int j = 0; j < nb_quad; j++)
          {
            CopyVector(feval(m), j, vec_u);
            
            GetInverse(mat.GetPointQuadratureBoundary(j), dfjm1);
            Mlt(dfjm1, vec_u, vec_v);
            Mlt(Fb_surf.WeightsND(j)*mat.GetDsQuadratureBoundary(j), vec_v);
            
            vec_vs = Fb_vol.TransposeTangentialVector(num_loc, vec_v);
            CopyVector(vec_vs, j, feval_surf);
          }
	
	res(m).Reallocate(nb_dof);
        
	// integration on the reference boundary
	Fb_surf.ApplyCh(feval_surf, res(m));
      }
  }

  
  //! computation of \int_\Gamma f \curl phi
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeIntegralSurfaceGradient(const SetMatrices<Dimension3>& mat,
				 Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
				 const MeshNumbering<Dimension3>& mesh_num,
                                 int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    int nb_points_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    
    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    // we apply 1/Ji DF_i^t on evaluation of f
    TinyVector<T, 3> E1, E2; Real_wp jacob;
    res.Reallocate(nb_unknowns);
    for (int i = 0; i < nb_points_quad; i++)
      {
	jacob = Det(mat.GetPointQuadratureBoundary(i));
	Real_wp coef = FaceBasis.WeightsQuadratureBoundary(i, num_loc)/jacob;
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
            E1(0) = feval(m)(3*i);
            E1(1) = feval(m)(3*i+1);
            E1(2) = feval(m)(3*i+2);
	    
            MltTrans(mat.GetPointQuadratureBoundary(i), E1, E2);
            Mlt(mat.GetDsQuadratureBoundary(i)*coef, E2);
	    
            feval(m)(3*i) = E2(0); 
            feval(m)(3*i+1) = E2(1); 
            feval(m)(3*i+2) = E2(2); 
	  } 
      }
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(FaceBasis.GetNbDof());
        FaceBasis.ComputeIntegralSurfaceCurlRef(feval(m), res(m), num_loc);
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }
  
  
  //! Computation of the values of u on nodal points of the boundary
  //! from values of u on degrees of freedom
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval nodal values of u on the element
    \param[out] res nodal values of u on the boundary
    \param[in] mesh considered mesh
    \param[in] n element number
    \param[in] num_loc local boundary number
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeValueBoundary(const SetMatrices<Dimension3>& mat,
		       const Vector<Vector<T> >& feval, Vector<Vector<TinyVector<T, 3> > >& res,
		       const Mesh<Dimension3>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    int nb_points_nodal_face = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    Vector<T> u_boundary;
    TinyVector<T, 3> E1, E2;
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_points_nodal_face);
	u_boundary.SetData(3*res(m).GetM(), reinterpret_cast<T*>(res(m).GetData()));
	// computation on the reference element
	FaceBasis.ComputeValueBoundaryRef(feval(m), u_boundary, num_loc);
	
	// application of DF_i^*-1
	for (int j = 0; j < nb_points_nodal_face; j++)
	  {
            E1(0) = u_boundary(3*j);
            E1(1) = u_boundary(3*j+1);
            E1(2) = u_boundary(3*j+2);
	    
	    MltTrans(mat.GetInversePointNodalBoundary(j), E1, E2);
	    
            u_boundary(3*j) = E2(0);
            u_boundary(3*j+1) = E2(1);
            u_boundary(3*j+2) = E2(2);
	  }
	
	u_boundary.Nullify();
      }
  }
  
  
  //! Computation of the curl of u on nodal points of the boundary
  //! from values of u on degrees of freedom
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval nodal values of u on the element
    \param[out] res nodal values of curl u on the boundary
    \param[in] mesh considered mesh
    \param[in] num_loc local boundary number
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeGradientBoundary(const SetMatrices<Dimension3>& mat,
			  const Vector<Vector<T> >& feval,
			  Vector<Vector<TinyVector<T, 3> > >& res,
			  const Mesh<Dimension3>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);

    int nb_points_nodal_face = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    if (!FaceBasis.UsePiolaTransform())
      {
	cout << "not implemented" << endl;
	abort();
      }

    Vector<T> grad_boundary;    
    TinyVector<T, 3> E1, E2; Real_wp jacob; 
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_points_nodal_face);
	grad_boundary.SetData(3*res(m).GetM(), reinterpret_cast<T*>(res(m).GetData()));
	
	// computation on the reference element
	FaceBasis.ComputeCurlBoundaryRef(feval(m), grad_boundary, num_loc);
	
	// application of 1/Ji DF_i
	for (int j = 0; j < nb_points_nodal_face; j++)
	  {
	    jacob = 1.0/Det(mat.GetPointNodalBoundary(j));
            E1(0) = grad_boundary(3*j);
            E1(1) = grad_boundary(3*j+1);
            E1(2) = grad_boundary(3*j+2);
            
            Mlt(mat.GetPointNodalBoundary(j), E1, E2);
	    
            grad_boundary(3*j) = jacob*E2(0);
            grad_boundary(3*j+1) = jacob*E2(1);
            grad_boundary(3*j+2) = jacob*E2(2);
	  }
	
	grad_boundary.Nullify();
      }
  }
  
  
  //! Extraction of the nodal values of on the boundary, from nodal values on all the element
  /*!
    \param[in] feval nodal values of u on the element
    \param[out] res nodal values of u on the boundary
    \param[in] num_loc local boundary number
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ComputeValueNodalBoundary(const Vector<T>& feval, Vector<T>& res, int num_loc) const
  {
    const ElementReference<Dimension3, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension3, 2>& >(*this);
    
    int nb_dof = FaceBasis.GetNbNodalBoundary(num_loc);
    res.Reallocate(nb_dof);
    for (int j = 0; j < nb_dof; j++)
      {
	int num_dof = FaceBasis.GetNodalNumber(num_loc, j);
	res(j) = feval(num_dof);
      }
  }
  
  
  //! Rotation of values of edge elements because of a rotation
  /*!
    \param[in] n mode number
    \param[in] fft_interface fft object
    \param[inout] val_u rotation of u of an angle n alpha
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ApplyRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                      TinyVector<T, 3>& val_u)
  {
    Real_wp cos_nalpha(1), sin_nalpha(0);
    fft_interface.GetCosSinAlpha(n, cos_nalpha, sin_nalpha);
    
    TinyVector<T, 3> val_v(val_u);
    val_u(0) = cos_nalpha * val_v(0) - sin_nalpha * val_v(1);
    val_u(1) = sin_nalpha * val_v(0) + cos_nalpha * val_v(1);
  }
  

  //! Rotation of values of edge elements because of a rotation
  /*!
    \param[in] n mode number
    \param[in] fft_interface fft object
    \param[inout] val_u rotation of u of an angle -n alpha
   */
  template<class T>
  void FiniteElementHcurl<Dimension3>::
  ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                             TinyVector<T, 3>& val_u)
  {
    Real_wp cos_nalpha(1), sin_nalpha(0);
    fft_interface.GetCosSinAlpha(n, cos_nalpha, sin_nalpha);
    
    TinyVector<T, 3> val_v(val_u);
    val_u(0) = cos_nalpha * val_v(0) + sin_nalpha * val_v(1);
    val_u(1) = -sin_nalpha * val_v(0) + cos_nalpha * val_v(1);
  }
  

  template<class T>
  void FiniteElementHcurl<Dimension3>
  ::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
			       Vector<T>& val_u)
  {
    int N = val_u.GetM() / 3;
    TinyVector<T, 3> vec_u;
    for (int k = 0; k < N; k++)
      {
	CopyVector(val_u, k, vec_u);
	ApplyInverseRotationCyclic(n, fft_interface, vec_u);
	CopyVector(vec_u, k, val_u);
      }
  }

}

#define MONTJOIE_FILE_FINITE_ELEMENT_HCURL3D_CXX
#endif
