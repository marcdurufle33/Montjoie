#ifndef MONTJOIE_FILE_FINITE_ELEMENT_HCURL2D_CXX

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
  void FiniteElementHcurl<Dimension2>::
  ComputeValuesPhi(const R2& point_loc, Vector<R2>& res, const Matrix2_2& dfjm1,
		   const MeshNumbering<Dimension2>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    R2 phi;
    FaceBasis.ComputeValuesPhiRef(point_loc, res);
    for (int i = 0; i < FaceBasis.GetNbDof(); i++)
      {
	phi = res(i);
	// applying transformation DF_i^{*-1}
	MltTrans(dfjm1, phi, res(i));	
      }
    
    // we take in account signs of basis functions 
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
  void FiniteElementHcurl<Dimension2>::
  ComputeValuesGradientPhi(const R2& point_loc, Vector<TinyVector<Real_wp, 1> >& res,
                           const Matrix2_2& dfjm1,
			   const MeshNumbering<Dimension2>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    Real_wp curl;
    res.Reallocate(nb_dof_loc);
    VectReal_wp curl_phi;
    FaceBasis.ComputeCurlPhiRef(point_loc, curl_phi);
    // computation of the curl of each basis function at point_loc
    for (int node = 0; node < nb_dof_loc; node++)
      {
	curl = Det(dfjm1)*curl_phi(node);
	res(node)(0) = curl;
      }
    
    // we take in account signs of basis functions 
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! Computes curl U from the gradient of U
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  GetCurlFromGradient(const Vector<Vector<T> >& dU, Vector<Vector<T> >& gradU)
  {
    int nb_unknowns = dU.GetM() / 4;
    gradU.Reallocate(nb_unknowns);
    for (int m = 0; m < nb_unknowns; m++)
      {
	int nb_nodes = dU(0).GetM();
	gradU(m).Reallocate(nb_nodes);
        for (int j = 0; j < nb_nodes; j++)
          {
            // curl U = dUy/dx - dUx/dy
            gradU(m)(j) = dU(4*m+2)(j) - dU(4*m+1)(j);
          }
      }
  }
  
  
  //! computation of tangential components of functions on quadrature points on the edge
  /*!
    \param[out] ValuePhi_Boundary values
    \param[in] MatricesElem jacobian matrices on the element
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
    \param[in] num_loc local edge number
   */
  // computes ValuePhi_Boundary(i,k) = \phi_i(\xi_k)
  // where \xi_k is a quadrature point on the edge num_loc
  void FiniteElementHcurl<Dimension2>::
  ComputeValuesPhiBoundary(Matrix<R2>& ValuePhi_Boundary,
			   const SetMatrices<Dimension2>& MatricesElem,
			   const MeshNumbering<Dimension2>& mesh_num,
			   int nquad, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_dof_loc = FaceBasis.GetNbDof();
    ValuePhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    FillZero(ValuePhi_Boundary);
    
    Vector<R2> val_phi_chap(nb_dof_loc);
    R2 val_phi; Matrix2_2 dfjm1;
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
    
    // we take into account signs of basis functions 
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, ValuePhi_Boundary, nquad);
  }
  
  
  //! computation of curl of basis functions on quadrature points on the edge/face
  /*!
    \param[out] GradientPhi_Boundary values of curl
    \param[in] MatricesElem jacobian matrices on the element
    \param[in] num_loc local edge/face number
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
   */
  void FiniteElementHcurl<Dimension2>::
  ComputeValuesGradientPhiBoundary(Matrix<TinyVector<Real_wp, 1> >& CurlPhi_Boundary,
				   const SetMatrices<Dimension2>& MatricesElem,
				   const MeshNumbering<Dimension2>& mesh_num,
				   int nquad, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    CurlPhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    FillZero(CurlPhi_Boundary);
    
    VectReal_wp curl_phi_chap(nb_dof_loc);
    Real_wp curl_phi; Real_wp jacob;

    for (int k = 0; k < nb_quad; k++)
      {
        jacob = 1.0/Det(MatricesElem.GetPointQuadratureBoundary(k));
        int num_point = FaceBasis.GetQuadNumber(num_loc, k);
        FaceBasis.GetCurlPhiOnQuadraturePoint(num_point, curl_phi_chap);
        
        for (int j = 0; j < nb_dof_loc; j++)
          {
            curl_phi = curl_phi_chap(j)*jacob;            
	    CurlPhi_Boundary(j, k) = curl_phi; 
	  }
      }
    
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, CurlPhi_Boundary, nquad);
  }
  
  
  //! Extraction of the nodal values of on the boundary, from nodal values on all the element
  /*!
    \param[in] feval nodal values of u on the element
    \param[out] res nodal values of u on the boundary
    \param[in] num_loc local boundary number
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeValueNodalBoundary(const Vector<T>& feval, Vector<T>& res, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_dof = FaceBasis.GetNbNodalBoundary(num_loc);
    res.Reallocate(nb_dof);
    for (int j = 0; j < nb_dof; j++)
      {
	int num_dof = FaceBasis.GetNodalNumber(num_loc, j);
	res(j) = feval(num_dof);
      } 
  }
  
  
  //! computation of U on nodal points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
    \param[in] mesh_num mesh numbering
    \param[in] iquad element number in the mesh
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeNodalValues(const SetMatrices<Dimension2>& MatricesElem,
		     const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
		     const Mesh<Dimension2>& mesh, int iquad) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_points_nodal = FaceBasis.GetNbPointsNodalElt();
    int nb_unknowns = Uloc.GetM();
    Uloc_node.Reallocate(2*nb_unknowns);
    for (int m = 0; m < 2*nb_unknowns; m++)
      Uloc_node(m).Reallocate(nb_points_nodal);
    
    Vector<T> Unode(2*nb_points_nodal);
    TinyVector<T, 2> val, vloc; Matrix2_2 dfjm1;
    
    // loop on all unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// projection on nodal points
	FaceBasis.ComputeNodalValuesRef(Uloc(m), Unode);
	
	// storing values on Uloc_node
	for (int i = 0; i < nb_points_nodal; i++)
	  {
	    Uloc_node(2*m)(i) = Unode(2*i);
	    Uloc_node(2*m+1)(i) = Unode(2*i+1);
	  }
      }
    
    // applying transformation DF_i^{*-1}
    for (int i = 0; i < nb_points_nodal; i++)
      {
	GetInverse(MatricesElem.GetPointNodal(i), dfjm1);
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    val(0) = Uloc_node(2*m)(i);
	    val(1) = Uloc_node(2*m+1)(i);
	    
	    MltTrans(dfjm1, val, vloc);
	    
	    Uloc_node(2*m)(i) = vloc(0);
	    Uloc_node(2*m+1)(i) = vloc(1);
	  }
      }	    
  }
  

  //! computation of U on quadrature points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on quadrature points
    \param[out] dUloc_node curls of U on quadrature points
    \param[in] compute_u if true, values of U are computed
    \param[in] compute_curl if true, curls of U are computed
    \param[in] mesh considered mesh
    \param[in] iquad element number in the mesh
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeQuadratureValues(const SetMatrices<Dimension2>& MatricesElem,
			  const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
			  Vector<Vector<T> >& dUloc_node, bool compute_u, bool compute_curl,
			  const Mesh<Dimension2>& mesh, int iquad) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_points_quad = FaceBasis.GetNbPointsQuadratureInside();
    int nb_unknowns = Uloc.GetM();
    if (compute_u)
      Uloc_node.Reallocate(2*nb_unknowns);
    
    if (compute_curl)
      dUloc_node.Reallocate(nb_unknowns);
      
    Vector<T> Uquad(2*nb_points_quad);
    TinyVector<T, 2> val, vloc; Matrix2_2 dfjm1; Real_wp invJacob;
    
    // loop on all unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// projection on quadrature points
	FaceBasis.ApplyChTranspose(Uloc(m), Uquad);
	if (compute_curl)
	  {
	    dUloc_node(m).Reallocate(nb_points_quad);
	    if (FaceBasis.UseQuadraturePointsForRh())
	      FaceBasis.ApplyRhQuadratureTranspose(Uquad, dUloc_node(m));
	    else
	      FaceBasis.ApplyRhTranspose(Uloc(m), dUloc_node(m));
	  }
	
	// storing values on Uloc_node
	if (compute_u)
	  {
	    Uloc_node(2*m).Reallocate(nb_points_quad);
	    Uloc_node(2*m+1).Reallocate(nb_points_quad);
	    for (int i = 0; i < nb_points_quad; i++)
	      {
		Uloc_node(2*m)(i) = Uquad(2*i);
		Uloc_node(2*m+1)(i) = Uquad(2*i+1);
	      }
	  }
      }
    
    // applying transformation DF_i^{*-1}
    for (int i = 0; i < nb_points_quad; i++)
      {
	GetInverse(MatricesElem.GetPointQuadrature(i), dfjm1);
	invJacob = Det(dfjm1);
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    if (compute_u)
	      {
		val(0) = Uloc_node(2*m)(i);
		val(1) = Uloc_node(2*m+1)(i);
		
		MltTrans(dfjm1, val, vloc);
		
		Uloc_node(2*m)(i) = vloc(0);
		Uloc_node(2*m+1)(i) = vloc(1);
	      }
	    
	    if (compute_curl)
	      dUloc_node(m)(i) *= invJacob;
	  }
      }    
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
  void FiniteElementHcurl<Dimension2>::
  ComputeValueBoundary(const SetMatrices<Dimension2>& mat,
		       const Vector<Vector<T> >& feval, Vector<Vector<TinyVector<T, 2> > >& res,
		       const Mesh<Dimension2>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);
        
    int nb_points_nodal_face = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    Vector<T> u_boundary;
    TinyVector<T, 2> E1, E2;
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_points_nodal_face);
	u_boundary.SetData(2*res(m).GetM(), reinterpret_cast<T*>(res(m).GetData()));
	// computation on the reference element
	FaceBasis.ComputeValueBoundaryRef(feval(m), u_boundary, num_loc);
	
	// application of DF_i^*-1
	for (int j = 0; j < nb_points_nodal_face; j++)
	  {
	    E1(0) = u_boundary(2*j);
            E1(1) = u_boundary(2*j+1);
            
	    MltTrans(mat.GetInversePointNodalBoundary(j), E1, E2);

            u_boundary(2*j) = E2(0);
            u_boundary(2*j+1) = E2(1);
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
    \param[in] n element number
    \param[in] num_loc local boundary number
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeGradientBoundary(const SetMatrices<Dimension2>& mat,
			  const Vector<Vector<T> >& feval, Vector<Vector<TinyVector<T, 1> > >& res,
			  const Mesh<Dimension2>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);
        
    int nb_points_nodal_face = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    Vector<T> grad_boundary;
    Real_wp jacob; 
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_points_nodal_face);
	grad_boundary.SetData(res(m).GetM(), reinterpret_cast<T*>(res(m).GetData()));
	
	// computation on the reference element
	FaceBasis.ComputeCurlBoundaryRef(feval(m), grad_boundary, num_loc);
	
	// application of 1/Ji DF_i
	for (int j = 0; j < nb_points_nodal_face; j++)
	  {
	    jacob = 1.0/Det(mat.GetPointNodalBoundary(j));
            grad_boundary(j) *= jacob;
	  }
	
	grad_boundary.Nullify();
      }
  }
  
  
  //! project function f on degrees of freedom
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on dof points
    \param[out] res components on degrees of freedom
    \param[in] mesh_num mesh numbering
    \param[in] n element number
  */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeProjectionDof(const SetMatrices<Dimension2>& MatricesElem,
		       Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
		       const MeshNumbering<Dimension2>& mesh_num, int n) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof(), nb_points_dof = FaceBasis.GetNbPointsDofInside();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    TinyVector<T, 2> vec_u, vec_v;
    for (int m = 0; m < nb_unknowns; m++)
      {
        for (int i = 0; i < nb_points_dof; i++)
	  {
	    vec_u(0) = feval(m)(2*i);
	    vec_u(1) = feval(m)(2*i+1);
	    MltTrans(MatricesElem.GetPointDof(i), vec_u, vec_v);
            feval(m)(2*i) = vec_v(0);
            feval(m)(2*i+1) = vec_v(1);
	  }
        
	res(m).Reallocate(nb_dof_loc);
	FaceBasis.ComputeProjectionDofRef(feval(m), res(m));
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }


  //! project function f on degrees of freedom
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on dof points
    \param[out] res components on degrees of freedom
    \param[in] mesh_num mesh numbering
    \param[in] n element number
    \param[in] num_loc local edge number
  */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeProjectionSurfaceDof(const SetMatrices<Dimension2>& MatricesElem,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			      const MeshNumbering<Dimension2>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);
        
    int nb_dof = FaceBasis.GetNbDofBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    Vector<Vector<T> > u_loc(1);
    u_loc(0).Reallocate(FaceBasis.GetNbDof());
    u_loc(0).Zero();
    res.Reallocate(nb_unknowns);
    
    if (FaceBasis.DiscontinuousElement())
      {
	cout << "This method is not compatible with discontinuous element" << endl;
	abort();
      }

    TinyVector<T, 2> vec_u, vec_v;
    for (int m = 0; m < nb_unknowns; m++)
      {
        for (int i = 0; i < FaceBasis.GetNbPointsDofSurface(num_loc); i++)
	  {
	    vec_u(0) = feval(m)(2*i);
	    vec_u(1) = feval(m)(2*i+1);
	    MltTrans(MatricesElem.GetPointDofBoundary(i), vec_u, vec_v);
            feval(m)(2*i) = vec_v(0);
            feval(m)(2*i+1) = vec_v(1);
	  }
        
	res(m).Reallocate(nb_dof);
	FaceBasis.ComputeProjectionSurfaceDofRef(feval(m), res(m), num_loc);

	// using a vector with all degrees of freedom
	// in order to take into account signs
	for (int j = 0; j < nb_dof; j++)
	  u_loc(0)(FaceBasis.GetLocalNumber(num_loc, j)) = res(m)(j);

	mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, u_loc, n);
	
	for (int j = 0; j < nb_dof; j++)
	  res(m)(j) = u_loc(0)(FaceBasis.GetLocalNumber(num_loc, j));
      }
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \varphi_i \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral against each basis function
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeIntegral(const SetMatrices<Dimension2>& mat,
		  Vector<Vector<T> >& feval, Vector<Vector<T> > & res,
		  const MeshNumbering<Dimension2>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_points_quadrature = FaceBasis.GetNbPointsQuadratureInside();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    Matrix2_2 dfjm1;
    TinyVector<T, 2> fvec, vec_u;
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	GetInverse(mat.GetPointQuadrature(i), dfjm1);
	Mlt(FaceBasis.WeightsND(i)*Det(mat.GetPointQuadrature(i)), dfjm1);
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    fvec(0) = feval(m)(2*i);
	    fvec(1) = feval(m)(2*i+1);
	    Mlt(dfjm1, fvec, vec_u);
	    feval(m)(2*i) = vec_u(0);
	    feval(m)(2*i+1) = vec_u(1);
	  }	
      }

    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
	
	FaceBasis.ComputeIntegralRef(feval(m), res(m));
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \nabla \times \varphi \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral agains gradient of each basis function
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeIntegralGradient(const SetMatrices<Dimension2>& mat,
			  Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			  const MeshNumbering<Dimension2>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);
        
    int nb_points_quadrature = FaceBasis.GetNbPointsQuadratureInside();
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
	for (int i = 0; i < nb_points_quadrature; i++)
	  feval(m)(i) *= FaceBasis.WeightsND(i);
	
	FaceBasis.ComputeIntegralCurlRef(feval(m), res(m));
      }

    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains each basis function
    \param[in] mesh_num mesh numbering
    \param[in] n element number
    \param[in] num_loc local edge number
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeIntegralSurface(const SetMatrices<Dimension2>& mat,
			 Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			 const MeshNumbering<Dimension2>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    
    Matrix2_2 dfjm1;
    TinyVector<T, 2> vec_u, vec_v;
    
    for (int j = 0; j < nb_quad; j++)
      {
	Real_wp poids = FaceBasis.Weights1D(j)*mat.GetDsQuadratureBoundary(j);
	GetInverse(mat.GetPointQuadratureBoundary(j), dfjm1);
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    vec_u(0) = feval(m)(2*j);
	    vec_u(1) = feval(m)(2*j+1);
	    
	    Mlt(dfjm1, vec_u, vec_v);
	    
	    feval(m)(2*j) = vec_v(0)*poids;
	    feval(m)(2*j+1) = vec_v(1)*poids;
	  }
      }
    
    res.Reallocate(nb_unknowns); 
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc); 
	FaceBasis.ComputeIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
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
  void FiniteElementHcurl<Dimension2>::
  ComputeGaussIntegralSurface(const Vector<Matrix2_2>& dfjm1,
                              const VectReal_wp& weights, const VectReal_wp& ds,
                              Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
                              const MeshNumbering<Dimension2>& mesh_num,
			      int n, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_points_quad = weights.GetM();
    int nb_unknowns = feval.GetM();
    
    // we apply DF_i^-1 on evaluation of f
    TinyVector<T, 2> E1, E2;
    res.Reallocate(nb_unknowns);
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(FaceBasis.GetNbDof());
	for (int i = 0; i < nb_points_quad; i++)
	  {
            E1(0) = feval(m)(2*i)*weights(i);
            E1(1) = feval(m)(2*i+1)*weights(i);
	    
	    Mlt(dfjm1(i), E1, E2);
	    
	    feval(m)(2*i) = E2(0);
            feval(m)(2*i+1) = E2(1);
          }
        
        FaceBasis.ComputeGaussIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }

  
  //! computation of \f$ \int_\Sigma f \phi_i' ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral against surfacic derivative of each basis function
    \param[in] mesh_num mesh numbering
    \param[in] n element number
    \param[in] num_loc local edge number
  */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ComputeIntegralSurfaceGradient(const SetMatrices<Dimension2>& mat,
				 Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
				 const MeshNumbering<Dimension2>& mesh_num,
				 int n, int num_loc) const
  {
    const ElementReference<Dimension2, 2>& FaceBasis
      = static_cast<const ElementReference<Dimension2, 2>& >(*this);

    int nb_points_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    
    // we apply 1/Ji on evaluation of f
    Real_wp jacob;
    for (int i = 0; i < nb_points_quad; i++)
      {
	jacob = Det(mat.GetPointQuadratureBoundary(i));
        
	Real_wp coef = FaceBasis.WeightsQuadratureBoundary(i, num_loc)
	  *mat.GetDsQuadratureBoundary(i)/jacob;
	
	for (int m = 0; m < nb_unknowns; m++)
	  feval(m)(i) *= coef;
      } 
    
    res.Reallocate(nb_unknowns);
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(FaceBasis.GetNbDof());
        FaceBasis.ComputeIntegralSurfaceCurlRef(feval(m), res(m), num_loc);
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);    
  }


  //! Rotation of values of edge elements because of a rotation
  /*!
    \param[in] n mode number
    \param[in] fft_interface fft object
    \param[inout] val_u rotation of u of an angle n alpha
   */
  template<class T>
  void FiniteElementHcurl<Dimension2>::
  ApplyRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                      TinyVector<T, 2>& val_u)
  {
    Real_wp cos_nalpha(1), sin_nalpha(0);
    fft_interface.GetCosSinAlpha(n, cos_nalpha, sin_nalpha);
    
    TinyVector<T, 2> val_v(val_u);
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
  void FiniteElementHcurl<Dimension2>::
  ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
                             TinyVector<T, 2>& val_u)
  {
    Real_wp cos_nalpha(1), sin_nalpha(0);
    fft_interface.GetCosSinAlpha(n, cos_nalpha, sin_nalpha);
    
    TinyVector<T, 2> val_v(val_u);
    val_u(0) = cos_nalpha * val_v(0) + sin_nalpha * val_v(1);
    val_u(1) = -sin_nalpha * val_v(0) + cos_nalpha * val_v(1);
  }

  
  template<class T>
  void FiniteElementHcurl<Dimension2>
  ::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
			       Vector<T>& val_u)
  {
    int N = val_u.GetM() / 2;
    TinyVector<T, 2> vec_u;
    for (int k = 0; k < N; k++)
      {
	CopyVector(val_u, k, vec_u);
	ApplyInverseRotationCyclic(n, fft_interface, vec_u);
	CopyVector(vec_u, k, val_u);
      }
  }
  
}

#define MONTJOIE_FILE_FINITE_ELEMENT_HCURL2D_CXX
#endif

