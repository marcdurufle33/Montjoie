#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_JACOBI_HXX

namespace Montjoie
{
  
  //! nodal finite element using Gauss-Jacobi points and DG method
  /*!
    test class, do not use it
   */
  class HexahedronDgJacobi : public HexahedronDgGauss
  {
    
  protected :
    Globatto<Real_wp> base1D, face1D, z1D, lobatto1D;

    VectReal_wp Points1D_B, Points1D_F, Points1D_Z, Points1D_L;
    //!< 1-D integration points (integration on edge)
    VectReal_wp Weights1D_B, Weights1D_F, Weights1D_Z; //!< 1-D integration weights

    Matrix<Real_wp> ValGauss_Extremity,ValJacobi_Extremity;
    Matrix<Real_wp> J_GL; //!< \f$ \varphi^{J}_i(\xi^{GL}_j) \f$
    
  public :
    HexahedronDgJacobi();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const;
    
    Real_wp GetValuePhi(int num_dof, int num_point) const;
    Real_wp GetValuePhiOnBoundary(int k, int j, int) const;
    
    const R3& GetGradientPhi(int num_dof, int num_point) const;
    const R3& GetGradientPhiOnBoundary(int node, int num_point, int num_loc) const;
    
    template<class TypeElement, class Vector1, class Vector2>
    void ComputeValueBoundary(const SetPoints<Dimension3>& pts,
			      const SetMatrices<Dimension3>& mat,
			      const Vector1& feval, Vector2& res,
			      TypeElement& FaceBasis, const Mesh<Dimension3>& mesh,
                              int n, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalValues(const SetMatrices<Dimension>& MatricesElem,
			    const Vector1& Uloc, Vector2& Uloc_node,
			    const Mesh<Dimension3>& mesh, int i) const;
    
    template<class Vector1, class Vector2>
    void ComputeValueBoundaryRef(const Vector1& feval, Vector2& res, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceRef(const Vector1 & feval,
				   Vector2& res, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceGradientRef(const Vector1 & feval,
					   Vector2& res, int num_loc) const;
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_DG_JACOBI_HXX
#endif
