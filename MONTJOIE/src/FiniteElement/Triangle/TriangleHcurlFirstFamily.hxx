#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! class for edge finite element of Nedelec's first family
  /*!
    An interpolatory basis is constructed with vandermonde matrix
    VDM_ij = psi_i(xi_j) \cdot t_j
    where xi_j are location of dofs and t_j direction of dofs
    psi_i are a set of nearly orthogonal functions generating R_r
    the classical Nedelec's first family space
   */
  class TriangleHcurlFirstFamily : public TriangleReference<2>
  {    
  public :
    //! dof points on edges can be chosen as regular points or gauss points
    int type_nodal_basis;
    enum {NODAL_REGULAR, NODAL_GAUSS, NODAL_LOBATTO};
    
  protected :
    //! dof numbers, tensorial <-> scalar
    /*! NumDofs_S1(i,j) returns the degree of freedom number
      with barycentric coordinate (i/(r+1),j/(r+1),(r+1-i-j)/(r+1))
      and oriented along Omega_1
      where \f$ \Omega_1 = \lambda_1 \nabla \lambda_2 - \lambda_2 \nabla \lambda_1 \f$
      analog for S2 and S3
      CoordinatesDofs(i,0) returns 1, 2 or 3, depending on
      // which vector is used to orient degree of freedom
      CoordinatesDof(i,1) (i,2) and (i,3) returns the barycentric coordinates
      multiplied by (r+1)
    */
    Matrix<int> NumDofs_S1;
    Matrix<int> NumDofs_S2; //!< dof numbers orientated along Omega2
    Matrix<int> NumDofs_S3; //!< dof numbers orientated along Omega3
    Matrix<int> CoordinateDofs; //!< scalar -> tensorial
    //! inverse of Vandermonde matrix for expression of nodal functions
    Matrix<Real_wp> InverseVDM_Edge;
    //! coefficient for orthonormalisation of orthogonal basis of Rr
    VectReal_wp InvWeightBasisRr;  
    
  public :
    TriangleHcurlFirstFamily();

    // how to number mesh
    void ConstructNumberMap(NumberMap& nmap, int dg) const;

    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
    void ConstructFunctions();

    void ComputeValuesPhiOrthoRef(const R2& point_loc, VectR2& res) const;
    void ComputeCurlPhiOrthoRef(const R2& pointloc, VectReal_wp& curl_phi) const;
    
    void ComputeValuesPhiRef(const R2& pointloc, VectR2& ) const;
    void ComputeCurlPhiRef(const R2& pointloc, VectReal_wp& ) const;

    void ComputeValuesPhiHRef(const R2& point_loc, VectReal_wp& res) const;
    
  };
  
  ostream& operator <<(ostream& out, const TriangleHcurlFirstFamily& e);

  class TriangleHcurlOptimalFirstFamily : public TriangleHcurlFirstFamily
  {
  public:
    TriangleHcurlOptimalFirstFamily();
    
  };

  ostream& operator <<(ostream& out, const TriangleHcurlOptimalFirstFamily& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_HCURL_FIRST_FAMILY_HXX
#endif

