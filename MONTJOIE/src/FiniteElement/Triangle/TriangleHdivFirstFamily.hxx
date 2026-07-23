#ifndef MONTJOIE_FILE_TRIANGLE_HDIV_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! class for facet finite element of Nedelec's first family
  /*!
    An interpolatory basis is constructed with Vandermonde matrix
    VDM_ij = psi_i(xi_j) \cdot t_j
    where xi_j are location of dofs and t_j their direction
    psi_i are nearly orthogonal functions of
    D_r = P_{r-1}^3 + [x, y] Ptilde_{r-1}
    where Ptilde_r is the polynomial space of degree exactly equal to r
   */
  class TriangleHdivFirstFamily : public TriangleReference<3>
  {    
  public:
    //! dof points on edges can be chosen as regular points or gauss points
    int type_nodal_basis;
    enum {NODAL_REGULAR, NODAL_GAUSS, NODAL_LOBATTO};
    
  protected :
    //! inverse of Vandermonde matrix for expression of nodal functions
    Matrix<Real_wp> InverseVDM_Edge;
    //! coefficient for orthonormalisation of orthogonal basis of Dr
    VectReal_wp InvWeightBasisRr;  
    
  public :
    TriangleHdivFirstFamily();

    // how to number mesh
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
  protected:
    void ConstructFunctions();

  public:
    void ComputeValuesPhiOrthoRef(const R2& point_loc, VectR2& res) const;
    void ComputeDivPhiOrthoRef(const R2& pointloc, VectReal_wp& curl_phi) const;
    
    void ComputeValuesPhiRef(const R2& pointloc, VectR2& ) const;
    void ComputeDivPhiRef(const R2& pointloc, VectReal_wp& ) const;

  };
  
  ostream& operator <<(ostream& out, const TriangleHdivFirstFamily& e);

  //! class for facet finite element of Nedelec's first family
  /*!
    This class is the same as TriangleHdivFirstFamily except that it
    is linked with optimal quadrilaterals
   */
  class TriangleHdivOptimalFirstFamily : public TriangleHdivFirstFamily
  {    
  public :
    TriangleHdivOptimalFirstFamily();
    
  };

  ostream& operator <<(ostream& out, const TriangleHdivOptimalFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_HDIV_FIRST_FAMILY_HXX
#endif

