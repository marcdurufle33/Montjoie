#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{
  
  inline Real_wp QuadrangleHcurlFirstFamily::WeightsGauss2D(int k) const
  {
    return WeightsGauss(k);
  }
  
  inline int QuadrangleHcurlFirstFamily::GetXdofNumber(int i, int j) const
  {
    return NumDofs_X(i, j);
  }
  
  inline int QuadrangleHcurlFirstFamily::GetYdofNumber(int i, int j) const
  {
    return NumDofs_Y(i, j);
  }
  
  inline const Matrix<Real_wp>& QuadrangleHcurlFirstFamily::GetGgl() const
  {
    return G_GL;
  }
  
  inline const Matrix<Real_wp>& QuadrangleHcurlFirstFamily::GetdGlg() const
  {
    return dGL_G;
  }
  
  inline int QuadrangleHcurlFirstFamily::GetNumNodes2DGauss(int i, int j) const
  {
    return NumNodes2D_gauss(i, j);
  }


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
  inline void QuadrangleHcurlFirstFamily::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(ch_loc, Uh, Vh);
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
  inline void QuadrangleHcurlFirstFamily::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, ch_loc, Uh, Vh);
  }
  

  //! Integration against curl of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \nabla \times \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = \nabla \times \varphi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  inline void QuadrangleHcurlFirstFamily::ApplyRhGauss(const Vector1& Vh, Vector2& Uh) const
  {
    Mlt(rh_loc, Vh, Uh);
  }
  
  
  //! computation of curl u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh curl of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j curl phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = curl phi_i(xi_j)
  */
  template<class Vector1, class Vector2>
  inline void QuadrangleHcurlFirstFamily
  ::ApplyRhGaussTranspose(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }
  
    
  inline void QuadrangleHcurlFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    this->ComputeProjectionDofGen(feval, contrib);
  }


  inline void QuadrangleHcurlFirstFamily
  ::ProjectQuadratureToDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    this->ProjectQuadratureToDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::ProjectQuadratureToDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    this->ProjectQuadratureToDofGen(feval, contrib);
  }
    
  inline void QuadrangleHcurlFirstFamily
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleHcurlFirstFamily
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleHcurlFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    this->ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void QuadrangleHcurlFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    this->ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void QuadrangleHcurlFirstFamily
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    this->ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    this->ComputeNodalValuesGen(Un, Unode);
  }

  inline void QuadrangleHcurlFirstFamily::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlFirstFamily::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    this->ApplyChGen(Vh, Uh);
  }

  inline void QuadrangleHcurlFirstFamily
  ::ApplyChTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyChTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::ApplyChTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyChTransposeGen(Uh, Vh);
  }

  inline void QuadrangleHcurlFirstFamily
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			  VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::AddVariableMassMatrix(int off_row, int off_col,
			  const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			  VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableMassMatrixGen(off_row, off_col, C, mat);
  }
  
  inline void QuadrangleHcurlFirstFamily
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const
  {
    this->AddVariableStiffnessMatrixGen(off_row, off_col, A, mat);
  }

  inline void QuadrangleHcurlFirstFamily
  ::AddVariableStiffnessMatrix(int off_row, int off_col,
			       const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const
  {
    this->AddVariableStiffnessMatrixGen(off_row, off_col, A, mat);
  }
  
}

#define MONTJOIE_FILE_QUADRANGLE_HCURL_FIRST_FAMILY_INLINE_CXX
#endif
