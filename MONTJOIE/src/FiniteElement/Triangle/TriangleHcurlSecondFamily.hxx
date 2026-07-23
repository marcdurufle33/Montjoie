#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_SECOND_FAMILY_HXX

namespace Montjoie
{
  //! class for edge finite element of Nedelec's second family
  /*!
    An interpolatory basis is used
    same interpolation points than for TriangleClassical
  */
  class TriangleHcurlSecondFamily : public TriangleReference<2>
  {
  protected :
    //! nodal point associated to each dof
    IVect NodalDof;
    //! direction of each basis function
    VectR2 DirectionDof; VectR2 normale_dof; //!< normale for each dof
    //! dof list for each node
    Matrix<int> ListeDof_Node, NumDofs2D;
    Matrix<int> CoordinateDofs;
    Vector<Matrix2_2> AhDof; //!< 2x2 matrix to solve on each node, for projection
    Matrix<Real_wp> InverseBasisVDM;
    VectReal_wp InvWeightFct;
    Matrix<int> NumFct2D;
    
  public :
    TriangleHcurlSecondFamily();
    
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
  protected:
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

  public:
    void ComputeValuesPhiRef(const R2& pointloc, VectR2& ) const;
    void ComputeCurlPhiRef(const R2& pointloc, VectReal_wp& ) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
  };
  
  ostream& operator <<(ostream& out,const TriangleHcurlSecondFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_HCURL_SECOND_FAMILY_HXX
#endif
