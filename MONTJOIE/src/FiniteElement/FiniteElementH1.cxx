#ifndef MONTJOIE_FILE_FINITE_ELEMENT_H1_CXX

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
  void FiniteElementH1<Dimension>::
  ComputeValuesPhi(const R_N& point_loc, Vector<TinyVector<Real_wp, 1> >& res,
                   const MatrixN_N& dfjm1, const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
    
    // we compute the values of all basis functions at point_loc
    res.Reallocate(FaceBasis.GetNbDof());
    VectReal_wp phi;
    FaceBasis.ComputeValuesPhiRef(point_loc, phi);

    if ((!mesh.IsElementAffine(nquad)) && (ElementReference_Base::use_warburton_trick))
      {
        Real_wp invSqrtJacob = 0;
        if (FaceBasis.UseQuadraturePointsForRh())
          {
            // evaluation of 1/sqrt(J) from values on quadrature points
            VectReal_wp phi_quad;
            FaceBasis.ComputeValuesPhiQuadratureRef(point_loc, phi_quad);
            
            int N = FaceBasis.GetNbPointsNodalElt();
            for (int j = 0; j < phi_quad.GetM(); j++)
              invSqrtJacob += phi_quad(j)*mesh.Glob_invSqrtJacobian(nquad)(N+j);
          }
        else
          {
            // evaluation of 1/sqrt(J) from values on nodal points
            VectReal_wp phi_nodal;
            FaceBasis.ComputeValuesPhiNodalRef(point_loc, phi_nodal);
	    
            for (int j = 0; j < phi_nodal.GetM(); j++)
              invSqrtJacob += phi_nodal(j)*mesh.Glob_invSqrtJacobian(nquad)(j);
          }
        
        phi *= invSqrtJacob;
      }

    for (int i = 0; i < res.GetM(); i++)
      res(i)(0) = phi(i);
    
    // we take into account signs of basis functions 
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
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
  void FiniteElementH1<Dimension>::
  ComputeValuesGradientPhi(const R_N& point_loc, Vector<R_N>& res, const MatrixN_N& dfjm1,
                           const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    // we compute the gradient of all basis functions at point_loc
    // need to apply transformation DF_i^{*-1} because
    // grad(u) = DF_i^{*-1} \hat{\nabla} u
    int nb_dof_loc = FaceBasis.GetNbDof();
    res.Reallocate(nb_dof_loc);
    FaceBasis.ComputeGradientPhiRef(point_loc, res); R_N grad;
    for (int node = 0; node < nb_dof_loc; node++)
      {
	MltTrans(dfjm1, res(node), grad);
	res(node) = grad;
      }

    if (ElementReference_Base::use_warburton_trick)
      abort();
    
    // we take into account signs of basis functions 
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! For compatibility with H^1 elements
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  GetCurlFromGradient(const Vector<Vector<T> >& dU, Vector<Vector<T> >& gradU)
  {
    gradU = dU;
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
  void FiniteElementH1<Dimension>::
  ComputeValuesPhiBoundary(Matrix<TinyVector<Real_wp, 1> >& ValuePhi_Boundary,
			   const SetMatrices<Dimension>& MatricesElem,
			   const MeshNumbering<Dimension>& mesh_num,
			   int nquad, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    ValuePhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    VectReal_wp phi;
    for (int k = 0; k < nb_quad; k++)
      {
        int num_point = FaceBasis.GetQuadNumber(num_loc, k);
        FaceBasis.GetValuePhiOnQuadraturePoint(num_point, phi);
        for (int i = 0; i < nb_dof_loc; i++)
          ValuePhi_Boundary(i, k)(0) = phi(i);
      }
    
    if ((!mesh.IsElementAffine(nquad)) && (ElementReference_Base::use_warburton_trick))
      {
        abort();
      }
    
    // we take into account signs of basis functions 
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
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
  void FiniteElementH1<Dimension>::
  ComputeValuesGradientPhiBoundary(Matrix<R_N>& GradientPhi_Boundary,
				   const SetMatrices<Dimension>& MatricesElem,
				   const MeshNumbering<Dimension>& mesh_num,
				   int nquad, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    GradientPhi_Boundary.Reallocate(nb_dof_loc, nb_quad);
    typename Dimension::R_N grad; typename Dimension::MatrixN_N dfjm1;
    Vector<R_N> grad_phi;
    for (int j = 0; j < nb_quad; j++)
      {
	GetInverse(MatricesElem.GetPointQuadratureBoundary(j), dfjm1);
        int num_point = FaceBasis.GetQuadNumber(num_loc, j);
        FaceBasis.GetGradientPhiOnQuadraturePoint(num_point, grad_phi);
	for (int i = 0; i < nb_dof_loc; i++)
	  {
	    MltTrans(dfjm1, grad_phi(i), grad);
	    GradientPhi_Boundary(i, j) = grad;
	  }
      }

    if (ElementReference_Base::use_warburton_trick)
      abort();
    
    // we take into account signs of basis functions 
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, GradientPhi_Boundary, nquad);
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
  void FiniteElementH1<Dimension>::
  ComputeNodalValues(const SetMatrices<Dimension>& MatricesElem,
		     const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
		     const Mesh<Dimension>& mesh, int nquad) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_points_nodal = FaceBasis.GetNbPointsNodalElt();
    
    // loop over unknowns
    int nb_unknowns = Uloc.GetM();
    if (Uloc_node.GetM() < nb_unknowns)
      Uloc_node.Reallocate(nb_unknowns);
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	Uloc_node(m).Reallocate(nb_points_nodal);
	
	// projection on nodal points
	FaceBasis.ComputeNodalValuesRef(Uloc(m), Uloc_node(m));
	
        if ((!mesh.IsElementAffine(nquad)) && (ElementReference_Base::use_warburton_trick))
          {
            for (int j = 0; j < nb_points_nodal; j++)
              Uloc_node(m)(j) *= mesh.Glob_invSqrtJacobian(nquad)(j);
          }
      }
  }
  
  
  //! computation of U on quadrature points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on quadrature points
    \param[out] dUloc_node gradients of U on quadrature points
    \param[in] compute_u if true, values of U are computed in Uloc_node
    \param[in] compute_grad if true, gradients of U are computed in dUloc_node
    \param[in] mesh considered mesh
    \param[in] iquad element number in the mesh
   */
  template<class Dimension> template<class T>  
  void FiniteElementH1<Dimension>::
  ComputeQuadratureValues(const SetMatrices<Dimension>& MatricesElem,
			  const Vector<Vector<T> >& Uloc, Vector<Vector<T> >& Uloc_node,
			  Vector<Vector<T> >& dUloc_node, bool compute_u, bool compute_grad,
			  const Mesh<Dimension>& mesh, int iquad) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_points_quad = FaceBasis.GetNbPointsQuadratureInside();
    TinyVector<T, Dimension::dim_N> grad_chap, grad_real;
    typename Dimension::MatrixN_N dfjm1;
    
    // loop on all unknowns
    int nb_unknowns = Uloc.GetM();
    if ((compute_u) && (Uloc_node.GetM() < nb_unknowns))
      Uloc_node.Reallocate(nb_unknowns);

    if ((compute_grad) && (dUloc_node.GetM() < nb_unknowns))
      dUloc_node.Reallocate(nb_unknowns);

    for (int m = 0; m < nb_unknowns; m++)
      {   
	if (compute_u)
	  Uloc_node(m).Reallocate(nb_points_quad);
    	
	// projection on quadrature points
	FaceBasis.ApplyChTranspose(Uloc(m), Uloc_node(m));
	if (compute_grad)
	  {
	    dUloc_node(m).Reallocate(Dimension::dim_N*nb_points_quad);
	    if (FaceBasis.UseQuadraturePointsForRh())
	      FaceBasis.ApplyRhQuadratureTranspose(Uloc_node(m), dUloc_node(m));
	    else
	      FaceBasis.ApplyRhTranspose(Uloc(m), dUloc_node(m));
	  }	
      }
    
    // applying transformation DF_i^{*-1}
    if (compute_grad)
      for (int i = 0; i < nb_points_quad; i++)
	{
	  GetInverse(MatricesElem.GetPointQuadrature(i), dfjm1);
	  
	  for (int m = 0; m < nb_unknowns; m++)
	    {
	      CopyVector(dUloc_node(m), i, grad_chap);
	      MltTrans(dfjm1, grad_chap, grad_real);
	      CopyVector(grad_real, i, dUloc_node(m));
	    }
	}
  }

  
  //! project function f on degrees of freedom
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on dof points
    \param[out] res components on degrees of freedom
    \param[in] n element number in the mesh
    \param[in] mesh_num mesh numbering
  */
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeProjectionDof(const SetMatrices<Dimension>& MatricesElem,
		       Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
		       const MeshNumbering<Dimension>& mesh_num, int n) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
        
    int nb_points_dof = FaceBasis.GetNbPointsDofInside();
    int nb_dof = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof);
        
        if ((!mesh.IsElementAffine(n)) && (ElementReference_Base::use_warburton_trick))
          {
            if (FaceBasis.DofEqualNodal())
              for (int i = 0 ; i < nb_points_dof; i++)
                feval(m)(i) /= mesh.Glob_invSqrtJacobian(n)(i);
            else
              abort();
          }
	
	FaceBasis.ComputeProjectionDofRef(feval(m), res(m));	
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.GetGlobalUnknownVector(mesh_num, res, n);
  }
  
  
  //! project function f on degrees of freedom (only surface dofs)
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] feval values of the function f on surface dof points
    \param[out] res components on degrees of freedom of the surface num_loc
    \param[in] mesh_num mesh numbering
    \param[in] n element number in the mesh
    \param[in] num_loc local face number
  */
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeProjectionSurfaceDof(const SetMatrices<Dimension>& MatricesElem,
			      Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
			      const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);
        
    int nb_dof = FaceBasis.GetNbDofBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    Vector<Vector<T> > u_loc(1);
    u_loc(0).Reallocate(FaceBasis.GetNbDof());
    u_loc(0).Zero();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);

    if (FaceBasis.DiscontinuousElement())
      {
	cout << "This method is not compatible with discontinuous element" << endl;
	abort();
      }

    for (int m = 0; m < nb_unknowns; m++)
      {       
	res(m).Reallocate(nb_dof);        
	FaceBasis.ComputeProjectionSurfaceDofRef(feval(m), res(m), num_loc);	
	
	// using a vector with all degrees of freedom
	// in order to take into account signs or linear combinations
	for (int j = 0; j < nb_dof; j++)
	  u_loc(0)(FaceBasis.GetLocalNumber(num_loc, j)) = res(m)(j);

	if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
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
  void FiniteElementH1<Dimension>::
  ComputeIntegral(const SetMatrices<Dimension>& mat,
		  Vector<Vector<T> > & feval, Vector<Vector<T> > & res,
		  const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    const Mesh<Dimension>& mesh = mesh_num.GetMesh();
        
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_points_quadrature = FaceBasis.GetNbPointsQuadratureInside();

    // initializations
    Real_wp jacob;
    bool divide_by_sqrt = ((ElementReference_Base::use_warburton_trick) && (!mesh.IsElementAffine(nquad)));
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    // in feval we store the evaluation on the quadrature
    // multiplied by the weight of integration (\omega_i J_i)
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	jacob = Det(mat.GetPointQuadrature(i));
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    if (divide_by_sqrt)
              feval(m)(i) *= FaceBasis.WeightsND(i)*sqrt(jacob);
            else
              feval(m)(i) *= FaceBasis.WeightsND(i)*jacob;
          }
      }
    
    for (int m = 0; m < nb_unknowns; m++)
      {	
	res(m).Reallocate(nb_dof_loc);
	// integral over the reference element \hat{K}
	FaceBasis.ComputeIntegralRef(feval(m), res(m));
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
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
  void FiniteElementH1<Dimension>::
  ComputeIntegralGradient(const SetMatrices<Dimension>& mat,
			  Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
			  const MeshNumbering<Dimension>& mesh_num, int nquad) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);
        
    int nb_points_quadrature = FaceBasis.GetNbPointsQuadratureInside();
    MatrixN_N dfjm1;
    TinyVector<T, Dimension::dim_N> fvec, vloc;
    
    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);

    // first step, we apply DF_i^{-1} transformation on f because
    // \int_K f \cdot \nabla \phi
    // \, = \, \int_{\hat{K}} J_i DF_i^{-1} f \cdot \hat{\nabla} \hat{\phi}
    // we multiply also by the weight of integration (\omega_i J_i)
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	GetInverse(mat.GetPointQuadrature(i), dfjm1);
	Real_wp coef = FaceBasis.WeightsND(i)*Det(mat.GetPointQuadrature(i));
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    CopyVector(feval(m), i, fvec);
	    Mlt(dfjm1, fvec, vloc);
	    vloc *= coef;
	    CopyVector(vloc, i, feval(m)); 
	  }
      }
    
    // now, we do integration over the reference element
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
	
	FaceBasis.ComputeIntegralGradientRef(feval(m), res(m));
      }
    
    if (ElementReference_Base::use_warburton_trick)
      abort();
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, nquad);
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral against each basis function
    \param[in] num_loc local edge number
    \param[in] mesh_num mesh numbering
    \param[in] FaceBasis given finite element
    \param[in] n element number
   */
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeIntegralSurface(const SetMatrices<Dimension>& mat,
			 Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
			 const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    
    // initializations
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    // loop on each unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we multiply the function f by the weight of integration (ds \omega_i)
	for (int j = 0; j < nb_quad; j++)
	  feval(m)(j) *= FaceBasis.WeightsQuadratureBoundary(j, num_loc)
	    *mat.GetDsQuadratureBoundary(j);
	
        if (ElementReference_Base::use_warburton_trick)
          abort();
        
	res(m).Reallocate(nb_dof_loc);
	// integration on the reference boundary
	FaceBasis.ComputeIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }  
  

  //! computes the integral of f against basis functions (points and weights are provided)
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeIntegralSurface(const Vector<typename Dimension::MatrixN_N>& dfjm1,
                         const VectReal_wp& weights, const VectReal_wp& ds,
                         Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
                         const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    
    // initializations
    int nb_quad = weights.GetM();
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    // loop on each unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we multiply the function f by the weight of integration (ds \omega_i)
	for (int j = 0; j < nb_quad; j++)
	  feval(m)(j) *= weights(j);

	res(m).Reallocate(nb_dof_loc);
	// integration on the reference boundary
	FaceBasis.ComputeIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    if (ElementReference_Base::use_warburton_trick)
      abort();
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
  }  


  //! computes the integral of f against basis functions (accurate integration)
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeGaussIntegralSurface(const Vector<typename Dimension::MatrixN_N>& dfjm1,
                              const VectReal_wp& weights, const VectReal_wp& ds,
                              Vector<Vector<T> >& feval, Vector<Vector<T> >& res,
                              const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    
    // initializations
    int nb_quad = weights.GetM();
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    // loop on each unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we multiply the function f by the weight of integration (ds \omega_i)
	for (int j = 0; j < nb_quad; j++)
	  feval(m)(j) *= weights(j);
	
	res(m).Reallocate(nb_dof_loc);
	// integration on the reference boundary
	FaceBasis.ComputeGaussIntegralSurfaceRef(feval(m), res(m), num_loc);
      }
    
    if (ElementReference_Base::use_warburton_trick)
      abort();
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
      mesh_num.number_map.ModifyLocalUnknownVector(mesh_num, res, n);
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
  void FiniteElementH1<Dimension>::
  ComputeIntegralSurfaceHDG(const SetMatrices<Dimension>& mat,
			    Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
                            const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& Fb_vol
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    const ElementReferenceTrace& Fb_surf = Fb_vol.GetSurfaceFiniteElement(num_loc);
    
    // initializations
    int nb_quad = Fb_surf.GetNbPointsQuadratureInside();
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);

    int nb_dof = Fb_surf.GetNbDof();
    
    // loop on each unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we multiply the function f by the weight of integration (ds \omega_i)
	for (int j = 0; j < nb_quad; j++)
	  feval(m)(j) *= Fb_surf.WeightsND(j)*mat.GetDsQuadratureBoundary(j);
	
        if (ElementReference_Base::use_warburton_trick)
          abort();
        
	res(m).Reallocate(nb_dof);
	// integration on the reference boundary
	Fb_surf.ApplyCh(feval(m), res(m));
      }
  }

  
  //! computation of \f$ \int_\Sigma f \phi_i' ds \f$
  /*!
    \param[in] mat jacobian matrices
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral against surfacic derivative of each basis function
    \param[in] num_loc local edge number
    \param[in] mesh_num mesh numbering
    \param[in] n element number
  */
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeIntegralSurfaceGradient(const SetMatrices<Dimension>& mat,
				 Vector<Vector<T> > & feval, Vector<Vector<T> >& res,
                                 const MeshNumbering<Dimension>& mesh_num, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_dof_loc = FaceBasis.GetNbDof();
    int nb_quad = FaceBasis.GetNbQuadBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    TinyVector<T, Dimension::dim_N> E1, E2;
    TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> dfjm1;
    
    // multiplying by ds \omega_i 
    for (int j = 0; j < nb_quad; j++)
      {
	Real_wp coef = FaceBasis.WeightsQuadratureBoundary(j, num_loc)
	  *mat.GetDsQuadratureBoundary(j);
	
	// we apply DF_i^-1 on evaluation of f
	GetInverse(mat.GetPointQuadratureBoundary(j), dfjm1);
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    CopyVector(feval(m), j, E1);
	    Mlt(dfjm1, E1, E2);
	    Mlt(coef, E2);
	    CopyVector(E2, j, feval(m));
	  }
      }
    
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_dof_loc);
    
	// integration over reference boundary
	FaceBasis.ComputeIntegralSurfaceGradientRef(feval(m), res(m), num_loc);
      }

    if (ElementReference_Base::use_warburton_trick)
      abort(); 
    
    if (mesh_num.number_map.FormulationDG() == ElementReference_Base::CONTINUOUS)
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
  void FiniteElementH1<Dimension>::
  ComputeValueNodalBoundary(const Vector<T>& feval, Vector<T>& res,
			    int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

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
    \param[in] pts reference points on the element
    \param[in] mat jacobian matrices
    \param[in] feval components of u on the element
    \param[out] res nodal values of u on the boundary
    \param[in] mesh considered mesh
    \param[in] n element number
    \param[in] num_loc local boundary number
   */
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::
  ComputeValueBoundary(const SetMatrices<Dimension>& mat,
		       const Vector<Vector<T> >& feval,
		       Vector<Vector<TinyVector<T, 1> > >& res,
                       const Mesh<Dimension>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_nodes = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);

    Vector<T> vec;
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_nodes);
	vec.SetData(nb_nodes, reinterpret_cast<T*>(res(m).GetData()));
        FaceBasis.ComputeValueBoundaryRef(feval(m), vec, num_loc);
	vec.Nullify();
      }
    
    if ((!mesh.IsElementAffine(n)) && (ElementReference_Base::use_warburton_trick))
      {
        for (int k = 0; k < nb_nodes; k++)
          {
            int num = FaceBasis.GetNodalNumber(num_loc, k);
            Real_wp invSqrtJacob = mesh.Glob_invSqrtJacobian(n)(num);
	    for (int m = 0; m < nb_unknowns; m++)
	      res(m)(k) *= invSqrtJacob;
          }
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
  void FiniteElementH1<Dimension>::
  ComputeGradientBoundary(const SetMatrices<Dimension>& mat,
			  const Vector<Vector<T> >& feval,
			  Vector<Vector<TinyVector<T, Dimension::dim_N> > >& res,
			  const Mesh<Dimension>& mesh, int n, int num_loc) const
  {
    const ElementReference<Dimension, 1>& FaceBasis
      = static_cast<const ElementReference<Dimension, 1>& >(*this);

    int nb_nodes = FaceBasis.GetNbNodalBoundary(num_loc);
    int nb_unknowns = feval.GetM();
    if (res.GetM() < nb_unknowns)
      res.Reallocate(nb_unknowns);
    
    Vector<T> vec;
    TinyVector<T, Dimension::dim_N> grad_chap;
    for (int m = 0; m < nb_unknowns; m++)
      {
	res(m).Reallocate(nb_nodes);
	vec.SetData(nb_nodes*Dimension::dim_N, reinterpret_cast<T*>(res(m).GetData()));
        FaceBasis.ComputeGradientBoundaryRef(feval(m), vec, num_loc);
	vec.Nullify();
    
	for (int j = 0; j < nb_nodes; j++)
	  {
	    grad_chap = res(m)(j);
	    MltTrans(mat.GetInversePointNodalBoundary(j), grad_chap, res(m)(j));
          }
      }
    
    if (ElementReference_Base::use_warburton_trick)
      abort();
  }  
  
  
  //! Local prolongation from an element to another (different orders)
  /*!
    \param[out] LocalProlongation prolongation operator
    \param[in] FaceCoarse "coarse" finite element to prolongate
    \param[in] FaceFine "fine" finite element
   */
  template<class Dimension>
  void FiniteElementH1<Dimension>::
  ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                           const ElementReference_Dim<Dimension>& FaceCoarse0,
			   const ElementReference_Dim<Dimension>& FaceFine0) const
  {
    const ElementReference<Dimension, 1>& FaceCoarse = dynamic_cast<const ElementReference<Dimension, 1>& >(FaceCoarse0);
    const ElementReference<Dimension, 1>& FaceFine = dynamic_cast<const ElementReference<Dimension, 1>& >(FaceFine0);
    if (FaceFine.GetOrder() <= FaceFine.GetMaximalOrderRestriction())
      {
        // computing prolongation operation as a matrix
        LocalProlongation.Reallocate(FaceFine.GetNbDof(), FaceCoarse.GetNbDof());
        VectReal_wp phi;
        for (int j = 0; j < FaceFine.GetNbDof(); j++)
          {
            FaceCoarse.ComputeValuesPhiRef(FaceFine.PointsDofND(j), phi);
            for (int k = 0; k < FaceCoarse.GetNbDof(); k++)
              LocalProlongation(j, k) = phi(k);
          }

        return;
      }
    
    // for quadrilaterals/hexahedra, we can have a fast operator (product of sparse matrices)
    proj.Init(FaceCoarse.GetGeometricElement(), FaceFine.PointsDof1D(), FaceFine.PointsDofND());
  }

  
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
						       TinyVector<T, 1>& val_u)
  {
  }
  
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
							      TinyVector<T, 1>& val_u)
  {
  }
    
  template<class Dimension> template<class T>
  void FiniteElementH1<Dimension>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>& fft_interface,
							      Vector<T>& val_u)
  {
  }
  
}

#define MONTJOIE_FILE_FINITE_ELEMENT_H1_CXX
#endif
