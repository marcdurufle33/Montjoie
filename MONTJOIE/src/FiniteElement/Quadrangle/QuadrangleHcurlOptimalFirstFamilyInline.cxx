#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_FIRST_FAMILY_INLINE_CXX

namespace Montjoie
{

  inline int QuadrangleHcurlOptimalFirstFamily::GetXdofNumber(int i, int j) const
  {
    return NumDofs_X(i, j);
  }
  
  inline int QuadrangleHcurlOptimalFirstFamily::GetYdofNumber(int i, int j) const
  {
    return NumDofs_Y(i, j);
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
  inline void QuadrangleHcurlOptimalFirstFamily::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
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
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, ch_loc, Uh, Vh);
  }


  inline void QuadrangleHcurlOptimalFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }

  inline void QuadrangleHcurlOptimalFirstFamily
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, Uh, Vh, num_loc, r);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, Uh, Vh, num_loc, r);
  }

  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, Uh, Vh, num_loc, r);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, Uh, Vh, num_loc, r);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplyChTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyChTransposeGen(Uh, Vh);
  }

  inline void QuadrangleHcurlOptimalFirstFamily
  ::ApplyChTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyChTransposeGen(Uh, Vh);
  }
    
}

#define MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_FIRST_FAMILY_INLINE_CXX
#endif
