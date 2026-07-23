#ifndef MONTJOIE_FILE_TETRAHEDRON_HDIV_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on tetrahedra
  /*!
    An interpolatory basis is constructed with Vandermonde matrix
    VDM_ij = psi_i(xi_j) \cdot t_j
    where xi_j are location of dofs and t_j their direction
    psi_i are nearly orthogonal functions of
    D_r = P_{r-1}^3 + [x, y, z] Ptilde_{r-1}
    where Ptilde_r is the polynomial space of degree exactly equal to r
   */
  class TetrahedronHdivFirstFamily : public TetrahedronReference<3>
  {
  protected :
    
    Matrix<int> CoordinateDofs; //!< tensorial <-> scalar matching array for dofs
    
    //! regular or non-regular nodal points ?
    int type_basis;
    enum {NODAL_REGULAR, NODAL_GAUSS, NODAL_LOBATTO};
    Matrix<Real_wp> InverseVDM_Edge;
    VectReal_wp InvWeightBasisRr;
    
  public :
    TetrahedronHdivFirstFamily();

    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected:
    void ConstructFunctions();
    
    void ConstructOrthogonalBasis(int r, int rgeom = 0, int rquad = 0, int type_quad = -1);
    
    void ComputeValuesPhiOrthoRef(const R3& point_loc, VectR3& res) const;
    void ComputeDivPhiOrthoRef(const R3& pointloc, VectReal_wp& div_phi) const;

  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeDivPhiRef(const R3& ptloc, VectReal_wp& dphi) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;        
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;

  };
  
  ostream& operator <<(ostream& out, const TetrahedronHdivFirstFamily& e);

  //! class for facet finite element of Nedelec's first family
  /*!
    This class is the same as TetrahedronHdivFirstFamily except that it
    is linked with optimal hexahedra, prisms and pyramids
   */
  class TetrahedronHdivOptimalFirstFamily : public TetrahedronHdivFirstFamily
  {
  public:
    TetrahedronHdivOptimalFirstFamily();
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_HDIV_FIRST_FAMILY_HXX
#endif

