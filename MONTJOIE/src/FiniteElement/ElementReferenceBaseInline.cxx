#ifndef MONTJOIE_FILE_ELEMENT_REFERENCE_BASE_INLINE_CXX

namespace Montjoie
{

  /*************************
   * ElementReference_Base *
   *************************/

  
  //! destructor
  inline ElementReference_Base::~ElementReference_Base()
  {
  }
  
  
  //! returns true if the mass is diagonal (or block-diagonal with 2x2 or 3x3 blocks)
  inline bool ElementReference_Base::LumpedMassMatrix() const
  {
    return mass_lumping;
  }
  

  //! returns true if the mass is diagonal for orthogonal element
  inline bool ElementReference_Base::MassLumpingOrthogonalElement() const
  {
    return mass_lumping_ortho;
  }
  
  
  //! returns true if the mass matrix is diagonal for any geometry
  inline bool ElementReference_Base::DiagonalMassMatrix() const
  {
    return diagonal_mass;
  }
  
  
  //! returns true if the element is discontinuous
  inline bool ElementReference_Base::DiscontinuousElement() const
  {
    return discontinuous_element;
  }
  

  //! returns true if Piola transform is used
  inline bool ElementReference_Base::UsePiolaTransform() const
  {
    return use_piola_transform;
  }


  //! returns true if Piola transform is used
  inline void ElementReference_Base::SetPiolaTransform(bool flag)
  {
    use_piola_transform = flag;
  }


  //! true if there is an optimized function for computing the mass matrix
  inline bool ElementReference_Base::OptimizedComputationMassMatrix() const
  {
    return optimized_mass_matrix;
  }
  
  
  //! true if there is an optimized function for computing elementary matrix
  inline bool ElementReference_Base::OptimizedComputationElementaryMatrix() const
  {
    return optimized_elem_matrix;
  }
  
  
  //! true if the mass matrix is sparse (even for a curved element)
  inline bool ElementReference_Base::SparseMassMatrix() const
  {
    return sparse_mass_matrix;
  }
  
  
  //! true if the mass matrix is sparse only for a non-curved element (but possibly non-affine)
  inline bool ElementReference_Base::LinearSparseMassMatrix() const
  {
    return linear_sparse_mass_matrix;
  }
  
  
  //! true if multiplication by operator Sh is using values on quadrature points
  inline bool ElementReference_Base::UseQuadraturePointsForSh() const
  {
    return use_quadrature_for_sh;
  }
  
  
  //! true if multiplication by operator Rh is using values on quadrature points
  inline bool ElementReference_Base::UseQuadraturePointsForRh() const
  {
    return use_quadrature_for_rh;
  }
  
  
  //! true if operator Sh can be directly applied without quadrature
  inline bool ElementReference_Base::UseQuadratureFreeSh() const
  {
    return use_quadrature_free_sh;
  }

  
  //! return the maximal order used to store the restriction matrix
  inline int ElementReference_Base::GetMaximalOrderRestriction() const
  {
    return 1;
  }
  
  
  //! returns order of approximation
  inline int ElementReference_Base::GetOrder() const
  {
    return order;
  }

  
  //! returns order used for geometry
  inline int ElementReference_Base::GetQuadratureOrder() const
  {
    return order_quad;
  }
  
  
  //! returns number of degrees of freedom
  inline int ElementReference_Base::GetNbDof() const
  {
    return nb_dof_loc;
  }

  
  //! returns number of quadrature points for volume integrals
  inline int ElementReference_Base::GetNbPointsQuadratureInside() const
  {
    return nb_points_quadrature_inside;
  }
  
  
  //! returns number of quadrature for surface integrals
  inline int ElementReference_Base::GetNbPointsQuadBoundaries() const
  {
    return nb_points_quadrature_boundaries;
  }

  
  //! returns number of points used for fluxes terms in DG
  inline int ElementReference_Base::GetNbPointsUsedForSh() const
  {
    return nb_points_quadrature_boundaries;
  }
  
  
  //! returns number of degrees of freedom located on the boundary
  inline int ElementReference_Base::GetNbDofBoundaries() const
  {
    return nb_dof_boundaries;
  }


  //! returns the number of dof points used to project a volume function
  inline int ElementReference_Base::GetNbPointsDofInside() const
  {
    return nb_points_dof_inside;
  }


  //! returns the number of dof points used to project a surface function
  inline int ElementReference_Base::GetNbPointsDofSurface(int num_loc) const
  {
    return num_dof_points_surf(num_loc).GetM();
  }

  
  //! returns the number of the dof point k of the edge/face num_loc
  inline int ElementReference_Base::GetPointDofNumber(int num_loc, int k) const
  {
    return num_dof_points_surf(num_loc)(k);
  }


  //! returns the number of the quadrature point k of the edge/face num_loc
  /*!
    \param[in] k local quadrature point number on the edge/face
    \param[in] num_loc edge/face number
   */
  inline int ElementReference_Base::GetQuadNumber(int num_loc, int k) const
  {
    return num_quad_points_surf(num_loc)(k);
  }


  inline const Vector<IVect>& ElementReference_Base::GetQuadNumbersBoundary() const
  {
    return num_quad_points_surf;
  }
  
  
  //! returns offset to acces quadrature points located on the boundary num_loc
  inline int ElementReference_Base::GetOffsetSh(int num_loc) const
  {
    return offset_faceSh(num_loc);
  }
    

  //! sets the numbers of quadrature points on the boundaries
  inline void ElementReference_Base::SetQuadNumbersBoundary(const Vector<IVect>& num_quad)
  {
    num_quad_points_surf = num_quad;
  }


  //! sets the numbers of dof points on the boundaries
  inline void ElementReference_Base::SetDofNumbersBoundary(const Vector<IVect>& num_dof)
  {
    num_dof_points_surf = num_dof;
  }
  
  
  //! returns the quadrature rule used on edges
  inline int ElementReference_Base::GetTypeIntegrationEdge() const
  {
    return type_integration_edge;
  }
    
  
  //! returns the quadrature rule used on triangulars
  inline int ElementReference_Base::GetTypeIntegrationTriangle() const
  {
    return type_integration_tri; 
  }
  
  
  //! returns the quadrature rule used on quadrilaterals
  inline int ElementReference_Base::GetTypeIntegrationQuadrangle() const
  {
    return type_integration_quad;
  }
    
  
  //! returns true if the dof belongs to the surface num_loc
  inline bool ElementReference_Base::IsTangentialDof(int j, int num_loc) const
  {
    if (j >= PosDofOnFace.GetM())
      return false;
    
    if (PosDofOnFace(j) & power_two_face(num_loc))
      return true;
    
    return false;
  }


  inline const Real_wp& ElementReference_Base::WeightsDofND(int i) const
  {
    return weights_dof(i);
  }

  
  inline const VectReal_wp& ElementReference_Base::WeightsDofND() const
  {
    return weights_dof;
  }
  

  //! returns weights omega_k for volume integrals
  inline const VectReal_wp& ElementReference_Base::WeightsND() const
  {
    return weights_quad_inside;
  }

    
  //! returns 1/2 omega_k where omega_k are integration weights on the face num_loc
  inline const VectReal_wp& ElementReference_Base
  ::GetFluxWeight(int num_loc) const
  {
    return PoidsFlux(num_loc);
  }
  

  /*******************
   * Virtual methods *
   *******************/
  
  
  //  inline void ElementReference_Base::CancelHighOrderDofs(IVect&, int r)
  // {
  //}
  
  
  inline void ElementReference_Base
  ::ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const
  {
    this->ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  
  inline void ElementReference_Base
  ::ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const
  {
    this->ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }

  inline void ElementReference_Base
  ::ApplyRh(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    abort();
  }
  
  
  inline void ElementReference_Base
  ::ApplyRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    abort();
  }

  
  inline void ElementReference_Base
  ::ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    abort();
  }
  
  
  inline void ElementReference_Base
  ::ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    abort();
  }

  
  inline void ElementReference_Base
  ::ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyRhQuadratureGen(Uh, Vh);
  }
  
  
  inline void ElementReference_Base
  ::ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyRhQuadratureGen(Uh, Vh);
  }
    
  
  inline void ElementReference_Base
  ::ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    this->ApplyRhQuadratureTransposeGen(Uh, Vh);
  }
  
  
  inline void ElementReference_Base
  ::ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    this->ApplyRhQuadratureTransposeGen(Uh, Vh);
  }


  inline void ElementReference_Base
  ::ApplyRhBoundary(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    abort();
  }
  
  
  inline void ElementReference_Base
  ::ApplyRhBoundary(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    abort();
  }
    
  
  inline void ElementReference_Base
  ::ApplyRhBoundaryTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const
  {
    abort();
  }
  
  
  inline void ElementReference_Base
  ::ApplyRhBoundaryTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const
  {
    abort();
  }

  
  inline void ElementReference_Base
  ::ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyShTransposeGen(num_loc, Uh, Vh, r);
  }
    
  
  inline void ElementReference_Base
  ::ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyShGen(alpha, num_loc, Uh, Vh, r);
  }
    
  
  inline void ElementReference_Base
  ::ApplyShQuadratureTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
    
  
  inline void ElementReference_Base
  ::ApplyShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r) const
  {
    this->ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

  
  inline void ElementReference_Base
  ::ApplyShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r) const
  {
    this->ApplyShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  
   inline void ElementReference_Base
  ::ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaShQuadratureTranspose(int num_loc, const VectReal_wp& Uh,
				    VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaShQuadratureTranspose(int num_loc, const VectComplex_wp& Uh,
				    VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShQuadratureTransposeGen(num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaShQuadrature(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			   VectReal_wp& Vh, int r) const
  {
    this->ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }
  
  
  inline void ElementReference_Base
  ::ApplyNablaShQuadrature(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			   VectComplex_wp& Vh, int r) const
  {
    this->ApplyNablaShQuadratureGen(alpha, num_loc, Uh, Vh, r);
  }

} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_REFERENCE_BASE_INLINE_CXX
#endif
