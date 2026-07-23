#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_GAUSS_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! edge finite element of Nedelec's first family on quadrilateral elements
  /*!
    Gauss and Gauss-Lobatto points are used for dof definitions
    (r+1)^2 Gauss points are used for integration
   */
  class QuadrangleHcurlGaussFirstFamily : public QuadrangleReference<2>
  {    
  protected :
    //! tensorial to scalar numbering for dofs orientated along X
    Matrix<int> NumDofs_X;
    //! tensorial to scalar numbering for dofs orientated along Y
    Matrix<int> NumDofs_Y;
    //! On each dof, normale
    VectR2 normale_dof;
    Matrix<int> CoordinateDofs;
    Globatto<Real_wp> lob_ortho;
    
  public :
    QuadrangleHcurlGaussFirstFamily();
    
    // how to number mesh
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
  protected:
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

  public:    
    void ComputeValuesPhiRef(const R2& pointloc, VectR2& ) const;
    void ComputeValuesPhiHRef(const R2& pointloc, VectReal_wp& ) const;
    void ComputeCurlPhiRef(const R2& pointloc, VectReal_wp& ) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
  };
  
  ostream& operator <<(ostream& out, const QuadrangleHcurlGaussFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HCURL_GAUSS_FIRST_FAMILY_HXX
#endif

