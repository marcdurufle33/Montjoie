#ifndef MONTJOIE_FILE_ELEMENT_REFERENCE_BASE_CXX

namespace Montjoie
{

  bool ElementReference_Base::use_warburton_trick(false);
    

  /*************************
   * ElementReference_Base *
   *************************/
  
  //! default constructor
  ElementReference_Base::ElementReference_Base(const VectReal_wp& poids)
    : weights_quad_inside(poids)
  {
    print_level = 0;
    order = -1;
    order_quad = -1;
    
    nb_points_quadrature_inside = 0;
    nb_points_quadrature_boundaries = 0;
    nb_dof_boundaries = 0;
    nb_dof_loc = 0;    
    nb_points_dof_inside = 0;
    
    mass_lumping = false;
    mass_lumping_ortho = false;
    diagonal_mass = false;
    discontinuous_element = false;
    use_piola_transform = true;
    
    optimized_mass_matrix = false;
    optimized_elem_matrix = false;
    sparse_mass_matrix = false;
    linear_sparse_mass_matrix = false;
    use_quadrature_for_sh = false;
    use_quadrature_for_rh = false;
    use_quadrature_free_sh = false;
    
    type_integration_edge = Globatto<Real_wp>::QUADRATURE_GAUSS;
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    type_integration_tri = TriangleQuadrature::QUADRATURE_GAUSS;
  }
  

  //! returns size of memory used by the object
  size_t ElementReference_Base::GetMemorySize() const
  {
    size_t size = sizeof(*this);
    size += Seldon::GetMemorySize(num_dof_points_surf);
    size += Seldon::GetMemorySize(num_quad_points_surf);
    size += sizeof(Real_wp)*this->MatProjectionDof.GetDataSize();
    size += sizeof(Real_wp)*this->tauProjectionDof.GetM();
    size += offset_faceSh.GetMemorySize() + power_two_face.GetMemorySize() + PosDofOnFace.GetMemorySize();
    size += sizeof(Real_wp)*mass_matrix_chol.GetDataSize();
    size += sizeof(Real_wp)*mass_matrix.GetDataSize();
    for (int i = 0; i < 6; i++)
      {
        size += sizeof(Real_wp)*PoidsFlux(i).GetM();
        size += sizeof(Real_wp)*const_nabla_sh(i).GetDataSize();
        size += sizeof(Real_wp)*sparse_const_sh(i).GetDataSize();
      }
    
    return size;
  }
  
  
  //! overwrites x by M^-1 x with M mass matrix on the reference element
  void ElementReference_Base::SolveMassMatrix(VectReal_wp& x) const
  {
    Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, x);
    Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, x);
  }


  //! overwrites x by M^-1 x with M mass matrix on the reference element
  void ElementReference_Base::SolveMassMatrix(VectComplex_wp& x) const
  {
    VectReal_wp xreal(x.GetM()), ximag(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
        xreal(i) = real(x(i));
        ximag(i) = imag(x(i));
      }
    
    Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, xreal);
    Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, xreal);

    Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, ximag);
    Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, ximag);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = Complex_wp(xreal(i), ximag(i));
  }
  
  
  //! overwrites x by L^-1 x or L^-t x with L Cholesky factor of mass matrix
  void ElementReference_Base
  ::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    Seldon::SolveCholesky(TransA, mass_matrix_chol, x);
  }


  //! overwrites x by L^-1 x or L^-t x with L Cholesky factor of mass matrix
  void ElementReference_Base
  ::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    Seldon::SolveCholesky(TransA, mass_matrix_chol, x);
  }

  
  //! overwrites x by M x with M mass matrix on the reference element
  void ElementReference_Base::MltMassMatrix(VectReal_wp& x) const
  {
    Seldon::MltCholesky(SeldonTrans, mass_matrix_chol, x);
    Seldon::MltCholesky(SeldonNoTrans, mass_matrix_chol, x);
  }
  
  
  //! overwrites x by M x with M mass matrix on the reference element
  void ElementReference_Base::MltMassMatrix(VectComplex_wp& x) const
  {
    VectReal_wp xreal(x.GetM()), ximag(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
        xreal(i) = real(x(i));
        ximag(i) = imag(x(i));
      }
    
    Seldon::MltCholesky(SeldonTrans, mass_matrix_chol, xreal);
    Seldon::MltCholesky(SeldonNoTrans, mass_matrix_chol, xreal);

    Seldon::MltCholesky(SeldonTrans, mass_matrix_chol, ximag);
    Seldon::MltCholesky(SeldonNoTrans, mass_matrix_chol, ximag);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = Complex_wp(xreal(i), ximag(i));
  }
  
  
  //! computation of mass matrix knowing decomposition of Jacobian for non-affine elements
  /*!
    \param[out] A computed mass matrix
    \param[in] coef decomposition of jacobian |DF|
   */
  void ElementReference_Base
  ::ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
                      const VectReal_wp & coef) const
  {
    // method overloaded if LinearSparseMassMatrix is true
    cout << "ComputeMassMatrix of base class, did you overload in the leaf class ?" << endl;
    abort();
  }
  
  
  //! computation of mass matrix knowing weighted jacobian on each quadrature point
  /*!
    \param[out] A computed mass matrix
    \param[in] coef values of jacobian |DF| on quadrature points
   */
  void ElementReference_Base
  ::IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
                        const VectReal_wp & coef) const
  {
    // method overloaded if SparseMassMatrix is true
    cout << "IntegrateMassMatrix of base class, did you overload in the leaf class ?" << endl;
    abort();
  }
  

  //! selection of near dofs
  /*!
    \param[in] pos position of the current dof
    \param[in] DofUsed array to know which dofs are already used
    \param[out] ListeDof dofs near the current dof and not used
    \param[in] nb_dof number of dofs we want to pick
   */
  void ElementReference_Base::PickNearDofs(int pos, const VectBool& DofUsed,
					   IVect& ListeDof, int nb_dof) const
  {
    ListeDof.Reallocate(nb_dof); ListeDof.Fill(-1); int nb = 0;
    for (int i = 0; i < DofUsed.GetM(); i++)
      if (!DofUsed(i))
	ListeDof(nb++) = i;
  }


  //! Integration against basis functions
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_K f \varphi_i dx
   */
  template<class Vector1,class Vector2>
  void ElementReference_Base::
  ComputeIntegralRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyCh(feval, res);
  }
   

  //! Integration against basis functions on a boundary
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc boundary number
   */  
  template<class Vector1,class Vector2>
  void ElementReference_Base::
  ComputeIntegralSurfaceRef(const Vector1 & feval,
                            Vector2& res, int num_loc) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    this->ApplySh(one, num_loc, feval, res);
  }
  

  //! Integration against basis functions on boundaries
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc boundary number
    The quadrature rules (\omega_k, \xi_k) are here classical Gauss rules
    and not the quadrature rules of the considered class
   */        
  template<class Vector1,class Vector2>
  void ElementReference_Base::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                 Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  //! Integration against gradient of basis functions associated with quadrature points
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \nabla \psi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void ElementReference_Base
  ::ApplyRhQuadratureGen(const Vector1& Vh, Vector2& Uh) const
  {
    abort();
  }
   

  //! Computation of gradient on quadrature points from values on quadrature points
  /*!
    \param[in] Vh values of u on quadrature points
    \param[out] Uh gradient of u on quadrature points
    This operation is equivalent to a matrix vector product
    Uh = Rh* Vh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference_Base
  ::ApplyRhQuadratureTransposeGen(const Vector1& Vh, Vector2& Uh) const
  {
    abort();
  }

  
  //! computation of u on quadrature points of a boundary
  /*!
    \param[in] n boundary number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyShTransposeGen(int n, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Mlt(SeldonTrans, sparse_const_sh(n), Uh, Vh);
  }
  

  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule of the boundary num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyShGen(const T0& alpha, int n, const Vector1& Uh, Vector2& Vh, int r) const
  {
    MltAdd(alpha, sparse_const_sh(n), Uh, T0(1), Vh);
  }
  
  
  //! computation of u on quadrature points of a boundary
  /*!
    \param[in] n boundary number
    \param[in] Uh values of u on quadrature points
    \param[out] Vh values of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \psi_i(\xi_j)  where \xi_j are quadrature points
    and \psi functions associated with quadrature points
   */  
  template<class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyShQuadratureTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    cout << "You can not use ShQuadrature for that finite element " << endl;
    abort();
  }
  
  
  //! integration against functions associated with quadrature points
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule of the boundary num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    cout << "You can not use ShQuadrature for that finite element " << endl;
    abort();
  }
  
  
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] n boundary number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyNablaShTransposeGen(int n, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Mlt(SeldonTrans, const_nabla_sh(n), Uh, Vh);
  }
  
  
  //! computation of grad u on quadrature points of a boundary
  /*!
    \param[in] n boundary number
    \param[in] Uh values of u on quadrature points
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \nabla \psi_i(\xi_j)  where \xi_j are quadrature points
    and \psi functions associated with quadrature points
   */  
  template<class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyNablaShQuadratureTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    cout << "You can not use NablaShQuadrature for that finite element " << endl;
    abort();
  }
  
  
  //! integration against gradient of functions associated with quadrature points
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule of the boundary num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyNablaShQuadratureGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    cout << "You can not use NablaShQuadrature for that finite element " << endl;
    abort();
  }
  
  
  //! integration against gradient of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule of the boundary num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void ElementReference_Base::
  ApplyNablaShGen(const T0& alpha, int n, const Vector1& Uh, Vector2& Vh, int r) const
  {
    MltAdd(alpha, const_nabla_sh(n), Uh, T0(1), Vh);
  }

}  

#define MONTJOIE_FILE_ELEMENT_REFERENCE_BASE_CXX
#endif
