#ifndef MONTJOIE_FILE_HEXAHEDRON_HDIV_OPTIMAL_HP_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for optimal Nedelec's first family on hexahedra
  /*!
    Hierarchical basis functions are used 
    and generate the space Q_{r+2, r, r} x Q_{r, r+2, r} x Q_{r, r, r+2}
    They are based on Jacobi polynomials
    Quadrature rules are Gauss points
   */
  class HexahedronHdivOptimalHpFirstFamily : public HexahedronReference<3>
  {
  protected :
    Array3D<int> NumDofsX, NumDofsY, NumDofsZ;
    
    Matrix<Real_wp> LegendrePolynom, JacobiPolynom;    
    Matrix<Real_wp> ValLeg, ValBubbleJac, DerBubbleJac;
    Matrix<Real_wp, General, RowSparse> ch1_node, ch2_node, ch3_node,
      ch1_loc, ch2_loc, ch3_loc, rh_loc;
    
    TinyVector<Matrix<Real_wp, General, RowSparse>, 6> sh_loc, sh_locX, sh_locY,
      sh_curl_loc, sh_curl_locX, sh_curl_locY, sh_nodeX, sh_nodeY, sh_curl_nodeX, sh_curl_nodeY;
    
    TinyMatrix<Matrix<Real_wp, General, ArrayRowSparse>, General, 3, 3>
    sparse_mass_matrix, sparse_stiff_matrix;
    
    Matrix<Real_wp, Symmetric, RowSymSparse> sparse_matrix_chol;
    Matrix<Real_wp> ValuePhi2D;
    
  public :

    int type_basis;
    enum {OPTIMAL, SUPER_OPTIMAL};

    HexahedronHdivOptimalHpFirstFamily();

    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;

  protected:
    void ConstructFunctions();    
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
  public:
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
        
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiRef(const R3& ptloc, VectReal_wp& dphi) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
    
  };
  
  ostream& operator <<(ostream& out, const HexahedronHdivOptimalHpFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_HDIV_OPTIMAL_HP_FIRST_FAMILY_HXX
#endif

