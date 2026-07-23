#ifndef MONTJOIE_FILE_QUADRANGLE_DG_GAUSS_INLINE_CXX

namespace Montjoie
{

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
  inline void QuadrangleDgGauss::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order != order_quad)
      return QuadrangleGauss::ApplyChGen(Uh, Vh);

    Vh = Uh;
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
  inline void QuadrangleDgGauss::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order != order_quad)
      return QuadrangleGauss::ApplyChTransposeGen(Uh, Vh);

    Vh = Uh;
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
  inline void QuadrangleDgGauss::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order != order_quad)
      return QuadrangleGauss::ApplyRhGen(Uh, Vh);

    Mlt(rh_loc, Uh, Vh);
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
  inline void QuadrangleDgGauss::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (order != order_quad)
      return QuadrangleGauss::ApplyRhTransposeGen(Uh, Vh);

    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  } 
  
  
  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi gradient of basis functions on quadrature point k
  */
  inline void QuadrangleDgGauss::GetGradientPhiOnQuadraturePoint(int k, VectR2& phi) const
  {
    ElementReferenceType<Dimension2, 1>::GetGradientPhiOnQuadraturePoint(k, phi);
  }


  inline void QuadrangleDgGauss::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }

  inline void QuadrangleDgGauss::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
    
  inline void QuadrangleDgGauss
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void QuadrangleDgGauss
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void QuadrangleDgGauss::SolveMassMatrix(VectReal_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void QuadrangleDgGauss::SolveMassMatrix(VectComplex_wp& x) const
  {
    SolveMassMatrixGen(x);
  }
  
  inline void QuadrangleDgGauss::SolveCholesky(const SeldonTranspose& TransA, VectReal_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
  
  inline void QuadrangleDgGauss::SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp& x) const
  {
    SolveCholeskyGen(TransA, x);
  }
    
  inline void QuadrangleDgGauss::MltMassMatrix(VectReal_wp& x) const
  {
    MltMassMatrixGen(x);
  }
  
  inline void QuadrangleDgGauss::MltMassMatrix(VectComplex_wp& x) const
  {
    MltMassMatrixGen(x);
  }

  inline void QuadrangleDgGauss::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleDgGauss::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleDgGauss::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
  
  inline void QuadrangleDgGauss::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  inline void QuadrangleDgGauss::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleDgGauss::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }

  inline void QuadrangleDgGauss::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleDgGauss::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }

  inline void QuadrangleDgGauss
  ::ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void QuadrangleDgGauss
  ::ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeValueBoundaryGen(Un, Unode, num_loc);
  }
    
  inline void QuadrangleDgGauss
  ::ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }
  
  inline void QuadrangleDgGauss
  ::ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const
  {
    ComputeGradientBoundaryGen(Un, Unode, num_loc);
  }

  inline void QuadrangleDgGauss
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleDgGauss
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleDgGauss
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleDgGauss
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleDgGauss
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleDgGauss
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  inline void QuadrangleDgGauss
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleDgGauss
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
        
  inline void QuadrangleDgGauss
  ::AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
  
  inline void QuadrangleDgGauss
  ::AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantMassMatrixGen(m, n, mass, A);
  }
    
  inline void QuadrangleDgGauss
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			  const TinyMatrix<Complex_wp, General, 2, 2>& C,
			  const R2_Complex_wp& D, const R2_Complex_wp& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Complex_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }

  inline void QuadrangleDgGauss
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			  const TinyMatrix<Real_wp, General, 2, 2>& C,
			  const R2& D, const R2& E,
			  const TinyVector<bool, 4>& null_term,
			  VirtualMatrix<Real_wp>& A) const
  {
    this->AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
  }
  
  inline void QuadrangleDgGauss
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			  const Vector<R2>& D, const Vector<R2>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
  inline void QuadrangleDgGauss
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			  const Vector<R2_Complex_wp>& D, const Vector<R2_Complex_wp>& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    this->AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, A);
  }
  
    
  inline void QuadrangleDgGauss
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }
  
  inline void QuadrangleDgGauss
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }  
  
}

#define MONTJOIE_FILE_QUADRANGLE_DG_GAUSS_INLINE_CXX
#endif
