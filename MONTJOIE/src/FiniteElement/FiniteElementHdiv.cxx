#ifndef MONTJOIE_FILE_FINITE_ELEMENT_HDIV_CXX

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
  template<class Dimension>
  void FiniteElementHdiv<Dimension>::
  ComputeValuesPhi(const R_N& point_loc,
		   Vector<R_N>& res, const MatrixN_N& dfjm1, 
		   const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    // we compute the values of all basis functions at point_loc
    res.Reallocate(FaceBasis.GetNbDof());
    FaceBasis.ComputeValuesPhiRef(point_loc, res);
    
    // multiplication by 1/Ji DFi (Piola transformation)
    MatrixN_N dfj; R_N vec_u;
    GetInverse(dfjm1, dfj);
    Real_wp jacob = Det(dfj);
    Mlt(1.0/jacob, dfj);
    for (int i = 0; i < res.GetM(); i++)
      {
	Mlt(dfj, res(i), vec_u);
	res(i) = vec_u;
      }
    
    // we take into account signs of basis functions 
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \nabla \varphi_i (\mbox{pointloc}) \f$
  /*!
    \param[in] point_loc local coordinates of the point on the unit element
    \param[out] res gradient of basis functions at local point
    \param[in] dfjm1 inverse of the jacobian matrix
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
  */
  template<class Dimension>
  void FiniteElementHdiv<Dimension>::
  ComputeValuesGradientPhi(const R_N& point_loc,
			   Vector<TinyVector<Real_wp, 1> >& res, const MatrixN_N& dfjm1,
                           const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);
        
    // we compute the gradient of all basis functions at point_loc
    int nb_dof_loc = FaceBasis.GetNbDof();
    res.Reallocate(nb_dof_loc);
    VectReal_wp div_phi(nb_dof_loc);
    FaceBasis.ComputeDivPhiRef(point_loc, div_phi);
    
    // applying relation div phi = 1/Ji \hat{div} \hat{phi}
    Real_wp invJacob = Det(dfjm1);
    for (int node = 0; node < nb_dof_loc; node++)
      res(node)(0) = div_phi(node)*invJacob;
    
    // we take into account signs of basis functions 
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! Computes div U from the gradient of U
  template<> template<class T>
  void FiniteElementHdiv<Dimension2>::
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
            // div U = dUx/dx + dUy/dy
            gradU(m)(j) = dU(4*m)(j) + dU(4*m+3)(j);
          }
      }    
  }


  //! Computes div U from the gradient of U
  template<> template<class T>
  void FiniteElementHdiv<Dimension3>::
  GetCurlFromGradient(const Vector<Vector<T> >& dU, Vector<Vector<T> >& gradU)
  {
    int nb_unknowns = dU.GetM() / 9;
    gradU.Reallocate(nb_unknowns);
    for (int m = 0; m < nb_unknowns; m++)
      {
	int nb_nodes = dU(0).GetM();
	gradU(m).Reallocate(nb_nodes);
        for (int j = 0; j < nb_nodes; j++)
          {
            // div U = dUx/dx + dUy/dy + dUz/dz
            gradU(m)(j) = dU(9*m)(j) + dU(9*m+4)(j) + dU(9*m+8)(j);
          }
      }    
  }
  
  
  //! computation of basis functions on quadrature points on the edge/face
  /*!
    \param[out] ValuePhi_Boundary values
    \param[in] MatricesElem jacobian matrices on the element
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
    \param[in] num_loc local edge number
  */
  // computes ValuePhi_Boundary(i,k) = \phi_i(\xi_k),
  // where \xi_k is a quadrature point on the edge num_loc
  template<class Dimension>
  void FiniteElementHdiv<Dimension>::
  ComputeValuesPhiBoundary(Matrix<R_N>& ValuePhi_Boundary,
			   const SetMatrices<Dimension>& MatricesElem,
			   const MeshNumbering<Dimension>& mesh_num,
			   int nquad, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    ValuePhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    Vector<R_N> phi; MatrixN_N dfj; Real_wp invJacob;
    for (int k = 0; k < nb_quad; k++)
      {
        int num_point = FaceBasis.GetQuadNumber(num_loc, k);
	dfj = MatricesElem.GetPointQuadratureBoundary(k);
	invJacob = 1.0/Det(dfj); Mlt(invJacob, dfj);
        FaceBasis.GetValuePhiOnQuadraturePoint(num_point, phi);	
        for (int i = 0; i < nb_dof_loc; i++)
	  Mlt(dfj, phi(i), ValuePhi_Boundary(i, k));
      }
    
    // we take into account signs of basis functions 
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, ValuePhi_Boundary, nquad);
  }
    
  
  //! computation of derivative of basis functions
  //! on quadrature points on the edge/face
  /*!
    \param[out] GradPhi_Boundary values of derivative
    \param[in] MatricesElem jacobian matrices on the element
    \param[in] num_loc local edge/face number
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
  */
  template<class Dimension>
  void FiniteElementHdiv<Dimension>::
  ComputeValuesGradientPhiBoundary(Matrix<TinyVector<Real_wp, 1> >& GradientPhi_Boundary,
				   const SetMatrices<Dimension>& MatricesElem,
				   const MeshNumbering<Dimension>& mesh_num,
				   int nquad, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    GradientPhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    Vector<Real_wp> div_phi;
    for (int j = 0; j < nb_quad; j++)
      {
	Real_wp invJacob = 1.0/Det(MatricesElem.GetPointQuadratureBoundary(j));
        int num_point = FaceBasis.GetQuadNumber(num_loc, j);
        FaceBasis.GetDivPhiOnQuadraturePoint(num_point, div_phi);
	for (int i = 0; i < nb_dof_loc; i++)
	  GradientPhi_Boundary(i, j)(0) = div_phi(i)*invJacob;
      }
    
    // we take into account signs of basis functions 
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, GradientPhi_Boundary, nquad);
  }
  
  
  //! computation of U on quadrature points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on quadrature points
    \param[out] dUloc_node divergence of U on quadrature points
    \param[in] compute_u if true, value of u is computed
    \param[in] compute_div if true, divergence of u is computed
    \param[in] mesh considered mesh
    \param[in] iquad element number in the mesh
   */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeQuadratureValues(const SetMatrices<Dimension>& MatricesElem,
			  const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
			  Vector<Vector<T> >& dUloc_node, bool compute_u, bool compute_div,
			  const Mesh<Dimension>& mesh, int iquad) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_points_quad = FaceBasis.GetNbPointsQuadratureInside();
    int nb_unknowns = Uloc.GetM();
    if (compute_u)
      Uloc_node.Reallocate(Dimension::dim_N*nb_unknowns);
    
    if (compute_div)
      dUloc_node.Reallocate(nb_unknowns);

    for (int m = 0; m < Uloc_node.GetM(); m++)
      Uloc_node(m).Reallocate(nb_points_quad);
      
    Vector<T> Uquad(Dimension::dim_N*nb_points_quad);
    TinyVector<T, Dimension::dim_N> val, vloc; Real_wp invJacob;
    typename Dimension::MatrixN_N dfj;
    
    // loop on all unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// projection on quadrature points
	FaceBasis.ApplyChTranspose(Uloc(m), Uquad);
	if (compute_div)
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
	    for (int i = 0; i < nb_points_quad; i++)
	      {
		CopyVector(Uquad, i, val);
		CopyVector(val, i, Dimension::dim_N*m, Uloc_node);
	      }
	  }
      }
    
    // applying transformation DF_i^{*-1}
    for (int i = 0; i < nb_points_quad; i++)
      {
	dfj = MatricesElem.GetPointQuadrature(i);
	invJacob = 1.0 / Det(dfj);
	Mlt(invJacob, dfj);
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    if (compute_u)
	      {
		ExtractVector(Uloc_node, i, Dimension::dim_N*m, val);
		
		MltTrans(dfj, val, vloc);
		
		CopyVector(vloc, i, Dimension::dim_N*m, Uloc_node);
	      }
	    
	    if (compute_div)
	      dUloc_node(m)(i) *= invJacob;
	  }
      }
  }
  
  
  //! computation of U on nodal points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
    \param[in] mesh considered mesh
    \param[in] i element number in the mesh
  */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeNodalValues(const SetMatrices<Dimension>& MatricesElem,
		     const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
		     const Mesh<Dimension>& mesh, int nquad) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_points_nodal = FaceBasis.GetNbPointsNodalElt();
    int nb_unknowns = Uloc.GetM();
    Uloc_node.Reallocate(Dimension::dim_N*nb_unknowns);
    for (int m = 0; m < Dimension::dim_N*nb_unknowns; m++)
      Uloc_node(m).Reallocate(nb_points_nodal);
    
    Vector<T> Unode(Dimension::dim_N*nb_points_nodal);
    TinyVector<T, Dimension::dim_N> val, vloc; MatrixN_N dfjm1;

    // loop over unknowns
    Real_wp invJacob;
    for (int m = 0; m < nb_unknowns; m++)
      {
	// projection on nodal points
	FaceBasis.ComputeNodalValuesRef(Uloc(m), Unode);
	
	// storing values on Uloc_node
	for (int i = 0; i < nb_points_nodal; i++)
	  {
	    CopyVector(Unode, i, val);
	    CopyVector(val, i, Dimension::dim_N*m, Uloc_node);
	  }
      }

    // applying Piola transform
    for (int i = 0; i < nb_points_nodal; i++)
      {	    
	invJacob = 1.0/Det(MatricesElem.GetPointNodal(i));
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    ExtractVector(Uloc_node, i, Dimension::dim_N*m, val);
	    Mlt(MatricesElem.GetPointNodal(i), val, vloc);
	    Mlt(invJacob, vloc);

	    CopyVector(vloc, i, Dimension::dim_N*m, Uloc_node);	    
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
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeProjectionDof(const SetMatrices<Dimension>& MatricesElem,
		       Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
		       const MeshNumbering<Dimension>& mesh_num, int n) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_points_dof = FaceBasis.GetNbPointsDofInside();
    int nb_dof = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    // taking into account Piola transform
    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    MatrixN_N dfjm1; Real_wp jacob;
    for (int i = 0; i < nb_points_dof; i++)
      {
	GetInverse(MatricesElem.GetPointDof(i), dfjm1);
	jacob = Det(MatricesElem.GetPointDof(i));
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    CopyVector(feval(m), i, vec_u);
	    
	    Mlt(dfjm1, vec_u, vec_v);
	    Mlt(jacob, vec_v);
	    
	    CopyVector(vec_v, i, feval(m));
	  }
      }
    
    // then projection on reference element
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof);
	FaceBasis.ComputeProjectionDofRef(feval(m), res(m));
      }
    
    mesh_num.number_map.GetGlobalUnknownVector(mesh_num, res, n);
  }


  //! project function f on degrees of freedom
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on dof points
    \param[out] res components on degrees of freedom
    \param[in] mesh_num mesh numbering
    \param[in] n element number in the mesh
    \param[in] num_loc local boundary number
  */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeProjectionSurfaceDof(const SetMatrices<Dimension>& MatricesElem,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			      const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_dof = FaceBasis.GetNbDofBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    Vector<Vector<T> > u_loc(1);
    u_loc(0).Reallocate(FaceBasis.GetNbDof());
    u_loc(0).Zero();
    
    if (FaceBasis.DiscontinuousElement())
      {
	cout << "This method is not compatible with discontinuous element" << endl;
	abort();
      }

    // taking into account Piola transform
    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    MatrixN_N dfjm1; Real_wp jacob;
    for (int i = 0; i < FaceBasis.GetNbPointsDofSurface(num_loc); i++)
      {
	GetInverse(MatricesElem.GetPointDofBoundary(i), dfjm1);
	jacob = Det(MatricesElem.GetPointDofBoundary(i));
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    CopyVector(feval(m), i, vec_u);
	    
	    Mlt(dfjm1, vec_u, vec_v);
	    Mlt(jacob, vec_v);
	    
	    CopyVector(vec_v, i, feval(m));
	  }
      }
    
    // then projection on reference element
    for (int m = 0; m < nb_unknowns; m++)
      {
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
    \param[out] res integral against each basis function
    \param[in] mesh_num mesh numbering
    \param[in] nquad eleemnt number
   */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeIntegral(const SetMatrices<Dimension>& mat,
		  Vector<Vector<T> >& feval, Vector<Vector<T> > & res,
		  const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_points_quadrature = FaceBasis.GetNbPointsQuadratureInside();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    // initializations
    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    
    // loop on each unknown
    for (int m = 0; m < nb_unknowns; m++)
      {
	// in feval we store the evaluation on the quadrature
	// multiplied by the weight of integration (\omega_i J_i)
	for (int i = 0; i < nb_points_quadrature; i++)
	  {
	    CopyVector(feval(m), i, vec_u);

	    MltTrans(mat.GetPointQuadrature(i), vec_u, vec_v);
	    Mlt(FaceBasis.WeightsND(i), vec_v);
	    
	    CopyVector(vec_v, i, feval(m));
          }
	
	// integral over the reference element \hat{K}
	res(m).Reallocate(nb_dof_loc);
	FaceBasis.ComputeIntegralRef(feval(m), res(m));
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \nabla \varphi \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral against gradient of each basis function
    \param[in] mesh_num mesh numbering
    \param[in] nquad element number
  */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeIntegralGradient(const SetMatrices<Dimension>& mat,
			  Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			  const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);
        
    int nb_points_quadrature = FaceBasis.GetNbPointsQuadratureInside();
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);

    // now, we do integration over the reference element
    for (int m = 0; m < nb_unknowns; m++)
      {
	for (int i = 0; i < nb_points_quadrature; i++)
	  feval(m)(i) *= FaceBasis.WeightsND(i);
	
	res(m).Reallocate(nb_dof_loc);
	FaceBasis.ComputeIntegralDivRef(feval(m), res(m));
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains each basis function
    \param[in] num_loc local edge number
    \param[in] mesh_num mesh numbering
    \param[in] n element number
   */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeIntegralSurface(const SetMatrices<Dimension>& mat,
			 Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			 const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    // initializations
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    Real_wp invJacob;
    
    for (int j = 0; j < nb_quad; j++)
      {
	invJacob = 1.0/Det(mat.GetPointQuadratureBoundary(j));
	invJacob *= FaceBasis.WeightsQuadratureBoundary(j, num_loc)*mat.GetDsQuadratureBoundary(j);
	
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    CopyVector(feval(m), j, vec_u);
	    
	    MltTrans(mat.GetPointQuadratureBoundary(j), vec_u, vec_v);
	    Mlt(invJacob, vec_v);
	    
	    CopyVector(vec_v, j, feval(m));
	  }
      }

    // integration on the reference boundary    
    for (int m = 0; m < nb_unknowns; m++)    
      {
	res(m).Reallocate(FaceBasis.GetNbDof());
	FaceBasis.ComputeIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
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
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeGaussIntegralSurface(const Vector<typename Dimension::MatrixN_N>& dfjm1,
                              const VectReal_wp& weights, const VectReal_wp& ds,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
                              const MeshNumbering<Dimension>& mesh_num,
			      int n, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    
    // initializations
    int nb_quad = weights.GetM();
    res.Reallocate(nb_unknowns);
    TinyVector<T, Dimension::dim_N> vec_u, vec_v;
    Real_wp invJacob; MatrixN_N dfj;
	
    for (int j = 0; j < nb_quad; j++)
      {
	GetInverse(dfjm1(j), dfj);
	invJacob = Det(dfjm1(j));
	for (int m = 0; m < nb_unknowns; m++)
	  {	    
	    CopyVector(feval(m), j, vec_u);
	    
	    MltTrans(dfj, vec_u, vec_v);	    
	    Mlt(invJacob*weights(j), vec_v);
	    
	    CopyVector(vec_v, j, feval(m));
	  }
      }

    for (int m = 0; m < nb_unknowns; m++)	
      {
	res(m).Reallocate(nb_dof_loc);
	// integration on the reference boundary
	FaceBasis.ComputeGaussIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }  
  
  
  //! computation of \f$ \int_\Sigma f \phi_i' ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains surfacic derivative of each basis function
    \param[in] num_loc local edge number
    \param[in] mesh_num mesh numbering
    \param[in] n element number
  */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeIntegralSurfaceGradient(const SetMatrices<Dimension>& mat,
				 Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
                                 const MeshNumbering<Dimension>& mesh_num,
				 int n, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    for (int j = 0; j < nb_quad; j++)
      {
	Real_wp invJacob = 1.0/Det(mat.GetPointQuadratureBoundary(j));
	Real_wp coef = FaceBasis.WeightsQuadratureBoundary(j, num_loc)
	  *mat.GetDsQuadratureBoundary(j)*invJacob;
	
	for (int m = 0; m < nb_unknowns; m++)
	  feval(m)(j) *= coef;
      }

    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
	// integration over reference boundary
	FaceBasis.ComputeIntegralSurfaceDivRef(feval(m), res(m), num_loc);
      }

    mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }

    
  //! Extraction of the nodal values of u on the boundary
  //! from nodal values on all the element
  /*!
    \param[in] feval nodal values of u on the element
    \param[out] res nodal values of u on the boundary
    \param[in] num_loc local boundary number
   */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeValueNodalBoundary(const Vector<T>& feval, Vector<T>& res, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_dof = FaceBasis.GetNbNodalBoundary(num_loc);
    res.Reallocate(nb_dof);
    for (int j = 0; j < nb_dof; j++)
      {
	int num_dof = FaceBasis.GetNodalNumber(num_loc, j);
	res(j) = feval(num_dof);
      }
  }
  
  
  //! Computation of the values of u on nodal points of the boundary
  //! from values of u on degrees of freedom
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval components of u on the element
    \param[out] res nodal values of u on the boundary
    \param[in] mesh considered mesh
    \param[in] n element number
    \param[in] num_loc local boundary number
   */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeValueBoundary(const SetMatrices<Dimension>& mat,
		       const Vector<Vector<T> >& feval,
                       Vector<Vector<TinyVector<T, Dimension::dim_N> > >& res,
		       const Mesh<Dimension>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);

    int nb_points_nodal_face = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    Vector<T> u_boundary;
    TinyVector<T, Dimension::dim_N> E1, E2;
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_points_nodal_face);
	u_boundary.SetData(Dimension::dim_N*res(m).GetM(),
			   reinterpret_cast<T*>(res(m).GetData()));
	
	// computation on the reference element
	FaceBasis.ComputeValueBoundaryRef(feval(m), u_boundary, num_loc);
	
	// application of 1/Ji DF_i
	for (int j = 0; j < nb_points_nodal_face; j++)
	  {
	    CopyVector(u_boundary, j, E1);
	    
	    Real_wp invJacob = 1.0/Det(mat.GetPointNodalBoundary(j));
	    Mlt(mat.GetPointNodalBoundary(j), E1, E2);
	    Mlt(invJacob, E2);
	    
	    CopyVector(E2, j, u_boundary);
	  }
	
	u_boundary.Nullify();
      }    
  }
  
  
  //! Computation of the values of gradient u on nodal points of the boundary
  //! from values of u on degrees of freedom
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval components of u on the element
    \param[out] res nodal values of gradient u on the boundary
    \param[in] mesh considered mesh
    \param[in] n element number
    \param[in] num_loc local boundary number
   */
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::
  ComputeGradientBoundary(const SetMatrices<Dimension>& mat,
			  const Vector<Vector<T> >& feval,
			  Vector<Vector<TinyVector<T, 1> > >& res,
			  const Mesh<Dimension>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension, 3>& FaceBasis
      = static_cast<const ElementReference<Dimension, 3>& >(*this);
        
    int nb_points_nodal_face = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    res.Reallocate(nb_unknowns);
    
    Vector<T> grad_boundary;
    Real_wp invJacob; 
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_points_nodal_face);
	grad_boundary.SetData(res(m).GetM(), reinterpret_cast<T*>(res(m).GetData()));
	
	// computation on the reference element
	FaceBasis.ComputeDivBoundaryRef(feval(m), grad_boundary, num_loc);
	
	// application of 1/Ji
	for (int j = 0; j < nb_points_nodal_face; j++)
	  {
	    invJacob = 1.0/Det(mat.GetPointNodalBoundary(j));	    
	    grad_boundary(j) *= invJacob;
	  }
	
	grad_boundary.Nullify();
      }
  }  
  
  
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
							 TinyVector<T, Dimension::dim_N>& val_u)
  {
    Real_wp cos_nalpha(1), sin_nalpha(0);
    fft_interface.GetCosSinAlpha(n, cos_nalpha, sin_nalpha);

    TinyVector<T, Dimension::dim_N> val_v(val_u);
    val_u(0) = cos_nalpha * val_v(0) - sin_nalpha * val_v(1);
    val_u(1) = sin_nalpha * val_v(0) + cos_nalpha * val_v(1);
  }
    
  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
								TinyVector<T, Dimension::dim_N>& val_u)
  {
    Real_wp cos_nalpha(1), sin_nalpha(0);
    fft_interface.GetCosSinAlpha(n, cos_nalpha, sin_nalpha);
    
    TinyVector<T, Dimension::dim_N> val_v(val_u);
    val_u(0) = cos_nalpha * val_v(0) + sin_nalpha * val_v(1);
    val_u(1) = -sin_nalpha * val_v(0) + cos_nalpha * val_v(1);

  }

  template<class Dimension> template<class T>
  void FiniteElementHdiv<Dimension>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
								Vector<T>& val_u)
  {
    int N = val_u.GetM() / Dimension::dim_N;
    TinyVector<T, Dimension::dim_N> vec_u;
    for (int k = 0; k < N; k++)
      {
	CopyVector(val_u, k, vec_u);
	ApplyInverseRotationCyclic(n, fft_interface, vec_u);
	CopyVector(vec_u, k, val_u);
      }
  }    
  
}

#define MONTJOIE_FILE_FINITE_ELEMENT_HDIV_CXX
#endif
