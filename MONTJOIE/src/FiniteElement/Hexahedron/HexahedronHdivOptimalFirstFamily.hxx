#ifndef MONTJOIE_FILE_HEXAHEDRON_HDIV_OPTIMAL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for optimal Nedelec's first family on hexhaedra
  /*!
    Nodal basis functions are used basd on Gauss and Gauss-Lobatto points
    phi_i^GL(x) phi_j^G(y) phi_k^G(z)  e_x
    phi_j^G(x) phi_i^GL(y) phi_k^G(z)  e_y
    phi_j^G(x) phi_k^G(y) phi_i^GL(z)  e_z
    with i <= r+2, j <= r, k <= r
    Gauss-Lobatto points are used for integration
   */
  class HexahedronHdivOptimalFirstFamily : public HexahedronReference<3>
  {
  protected :
    
    Globatto<Real_wp> gauss; //!< Gauss formula
    //! tensorial <-> scalar matching array for dofs orientated along x
    Seldon::Array3D<int> NumDofs_X;
    //! tensorial <-> scalar matching array for dofs orientated along y 
    Seldon::Array3D<int> NumDofs_Y;
    //! tensorial <-> scalar matching array for dofs orientated along z
    Seldon::Array3D<int> NumDofs_Z;
    Matrix<int> CoordinateDofs; //!< tensorial <-> scalar matching array for dofs
    
    // G = Gauss, GL = Gauss-Lobatto
    Matrix<Real_wp> G_GL; 
    Matrix<Real_wp> GL_G; 
    Matrix<Real_wp> dGL_G;
    Matrix<Real_wp> dGL_GL, stiff1d;

    Matrix<Real_wp, General, RowSparse> rh_loc; //!< stiffness matrix
    Matrix<Real_wp, General, RowSparse> ch_loc; //!< mass matrix
    
    //!< stiffness and mass matrices
    Matrix<Real_wp, General, RowSparse> ch1_loc; //!< mass matrix
    Matrix<Real_wp, General, RowSparse> ch2_loc; //!< mass matrix
    Matrix<Real_wp, General, RowSparse> ch3_loc; //!< mass matrix
    Matrix<Real_wp, General, RowSparse> rh_exact; //!< stiffness matrix
    Matrix<Real_wp, General, RowSparse> ch1_node, ch2_node, ch3_node;
    
    //! interpolation from r Gauss points to r+1 Gauss points
    Matrix<Real_wp> G_GX; 
    //! interpolation from r+1 Gauss-Lobatto points to r+1 Gauss points
    Matrix<Real_wp> GL_GX;
    //! derivative from r+1 Gauss points to r+1 Gauss points
    Matrix<Real_wp> dGX_GX; 
    //! integration weights for r+1 Gauss points
    VectReal_wp weights3d_Exact; 
    //! r+1 Gauss points
    Globatto<Real_wp> gauss_exact; 
    int nb_points_quadrature_gauss; //!< number of points for H
    VectReal_wp WeightsMassG, invWeightsMassG, invSqrtWeightsMassG;
    //! coefficients used to evaluate Legendre polynomials
    Matrix<Real_wp> LegendrePolynom;
    bool super_optimal;
    
  public :
    
    HexahedronHdivOptimalFirstFamily();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    // construction of orthogonal functions
    void ConstructOrthogonalBasis(int r, bool super_opt = false);
    void ComputeValuesPhiOrthoRef(const R3& x, VectR3& phi);
    
  protected:
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();        
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& res) const;

    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const;

  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiRef(const R3& ptloc, VectReal_wp& dphi) const;
        
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
    
  };
  
  ostream& operator <<(ostream& out, const HexahedronHdivOptimalFirstFamily& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_HDIV_OPTIMAL_FIRST_FAMILY_HXX
#endif
