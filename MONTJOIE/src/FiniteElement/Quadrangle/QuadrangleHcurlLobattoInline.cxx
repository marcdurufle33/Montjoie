#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_LOBATTO_INLINE_CXX

namespace Montjoie
{

  inline int QuadrangleHcurlLobatto::GetXdofNumber(int i) const
  {
    return Node_to_DofX(i);
  }
  
  inline int QuadrangleHcurlLobatto::GetYdofNumber(int i) const
  {
    return Node_to_DofY(i);
  }
  
  inline int QuadrangleHcurlLobatto::GetXdofNumber(int i, int j) const
  {
    return NumDofs_X(i, j);
  }
  
  inline int QuadrangleHcurlLobatto::GetYdofNumber(int i, int j) const
  {
    return NumDofs_Y(i, j);
  }
  
  inline const Real_wp& QuadrangleHcurlLobatto::dGL_GL(int i, int j) const
  {
    return lob_basis.GradPhi(i, j);
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
  inline void QuadrangleHcurlLobatto::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
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
  inline void QuadrangleHcurlLobatto::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }


  //! projection from quadrature points (Gauss-Lobatto points) to dof components
  /*!
    \param[in] Equad values on quadrature points
    \param[out] Edof dof components
  */
  template<class Vector1, class Vector2>
  inline void QuadrangleHcurlLobatto::
  ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const
  {
    ComputeProjectionDofRef(Equad, Edof);
  }
  
  
  //! H and E are discretized on Gauss-Lobatto points
  template<class Vector1, class Vector2>
  inline void QuadrangleHcurlLobatto::
  ComputeNodalValuesH(const Vector1& Un, Vector2 & Unode) const
  {
    Unode = Un;
  }  

  inline void QuadrangleHcurlLobatto::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyChTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyChTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyChTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyChTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyRhGen(Vh, Uh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
  
  inline void QuadrangleHcurlLobatto::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    ApplyRhTransposeGen(Uh, Vh);
  }
    
  inline void QuadrangleHcurlLobatto
  ::ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }
  
  inline void QuadrangleHcurlLobatto
  ::ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const
  {
    ComputeNodalValuesGen(Un, Unode);
  }

  inline void QuadrangleHcurlLobatto
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  inline void QuadrangleHcurlLobatto
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

  inline void QuadrangleHcurlLobatto
  ::ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }
  
  inline void QuadrangleHcurlLobatto
  ::ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const
  {
    ProjectQuadratureToDofGen(Equad, Edof);
  }
  
  inline void QuadrangleHcurlLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
				   VectReal_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
    
  inline void QuadrangleHcurlLobatto
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
				   VectComplex_wp& res, int num_loc) const
  {
    ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void QuadrangleHcurlLobatto
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
	    VectReal_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleHcurlLobatto
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
	    VectComplex_wp& Vh, int r) const
  {
    ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }

  inline void QuadrangleHcurlLobatto
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  inline void QuadrangleHcurlLobatto
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
      
}

#define MONTJOIE_FILE_QUADRANGLE_HCURL_LOBATTO_INLINE_CXX
#endif
