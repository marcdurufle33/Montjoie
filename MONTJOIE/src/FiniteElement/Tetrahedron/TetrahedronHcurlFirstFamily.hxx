#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on tetrahedra
  /*!
    An interpolatory basis is constructed with vandermonde matrix
    VDM_ij = psi_i(xi_j) \cdot t_j
    where xi_j are location of dofs and t_j direction of dofs
    psi_i are a set of nearly orthogonal functions generating R_r
    the classical Nedelec's first family space
   */
  class TetrahedronHcurlFirstFamily : public TetrahedronReference<2>
  {
  protected :
    
    Matrix<int> CoordinateDofs; //!< tensorial <-> scalar matching array for dofs
    
    //! basis functions
    Vector<MultivariatePolynomial<Real_wp> > Poly_Phi_X, Poly_Phi_Y, Poly_Phi_Z;
    //! curl of basis functions
    Vector<MultivariatePolynomial<Real_wp> > Poly_CurlPhi_X, Poly_CurlPhi_Y, Poly_CurlPhi_Z;
    //! regular or non-regular nodal points ?
    int type_basis;
    enum {NODAL_REGULAR, NODAL_GAUSS, NODAL_LOBATTO};
    Matrix<Real_wp> InverseVDM_Edge;
    VectReal_wp InvWeightBasisRr;
    
  public :
    TetrahedronHcurlFirstFamily();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    void ConstructOrthogonalBasis(int r, int rgeom = 0, int rquad = 0, int type_quad = -1);

  protected:
    void ConstructFunctions();
    
  public:
    void ComputeValuesPhiOrthoRef(const R3& point_loc, VectR3& res) const;
    void ComputeCurlPhiOrthoRef(const R3& pointloc, VectR3& curl_phi) const;

    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
        
  };
  
  ostream& operator <<(ostream& out, const TetrahedronHcurlFirstFamily& e);

  //! same class as TetrahedronHcurlFirstFamily except it is linked 
  //! with optimal hexahedra, prisms and pyramids
  class TetrahedronHcurlOptimalFirstFamily : public TetrahedronHcurlFirstFamily
  {
  public:
    TetrahedronHcurlOptimalFirstFamily();
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_FIRST_FAMILY_HXX
#endif

