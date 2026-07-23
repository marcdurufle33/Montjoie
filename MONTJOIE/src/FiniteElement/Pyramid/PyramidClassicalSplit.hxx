#ifndef MONTJOIE_FILE_PYRAMID_CLASSICAL_SPLIT_HXX

namespace Montjoie
{
  
  //! nodal finite element on classic pyramid
  /*!
    test class with nodal basis functions when the pyramid is viewed as two tetrahedra
    Do not use this finite element
   */
  class PyramidClassicalSplit : public PyramidReferenceSplit, public FiniteElementH1<Dimension3>
  {
    friend class FiniteElementH1<Dimension3>;
    
  public :
    typedef DenseProjector<Dimension3> ProjectionOperator; //!< volumetric projector
    typedef DenseProjector<Dimension2> ProjectionSurfacicOperator; //!< surfacic projector
    
    static const int ORDER_MAX_RESTRICTION_STD = 1; //!< matrix is stored
    static const int ORDER_MAX_MLT_STD = 1; //!< matrix is stored
    
  protected :    
    //! stiffness matrixij = \f$ \int \frac{\partial u}{\partial x_i}
    //! \frac{\partial v}{\partial x_j} \f$
    Matrix<Real_wp> stiffness_matrix00, stiffness_matrix11, stiffness_matrix22,
      stiffness_matrix01, stiffness_matrix02, stiffness_matrix12,
      stiffness_matrix10, stiffness_matrix20, stiffness_matrix21;
    //! gradient_matrixi = \f$ \int \frac{\partial u}{\partial x_i} v \f$
    Matrix<Real_wp> gradient_matrix0, gradient_matrix1, gradient_matrix2;
    
  public :
    PyramidClassicalSplit();
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();

    template<class T, class Matrix1, class Vector1>
    void ComputeElementaryInteraction(int i, int j, const T& mass, const Matrix1& C,
				      const Vector1& D, const Vector1& E, T& vloc) const;
    
    void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const;
    
    Real_wp GetValuePhi(int num_dof, int num_point) const;
    Real_wp GetValuePhiOnBoundary(int k, int j, int) const;
    
    const R3& GetGradientPhi(int num_dof, int num_point) const;
    const R3& GetGradientPhiOnBoundary(int node, int num_point, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralRef(const Vector1 & feval, Vector2 & res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralGradientRef(const Vector1 & feval, Vector2 & res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceRef(const Vector1 & feval, Vector2 & res, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceGradientRef(const Vector1 & feval,
                                           Vector2 & res, int num_loc) const;

    friend ostream& operator <<(ostream& out,const PyramidClassicalSplit& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_PYRAMID_CLASSICAL_SPLIT_HXX
#endif
