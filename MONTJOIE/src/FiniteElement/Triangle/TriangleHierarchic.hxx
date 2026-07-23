#ifndef MONTJOIE_FILE_TRIANGLE_HIERARCHIC_HXX

namespace Montjoie
{
  
  //! hp finite element for triangles
  /*!
    Hierarchical basis functions are constructed with Jacobi polynomials
    and generate the polynomial space P_r
   */
  class TriangleHierarchic : public TriangleReference<1>
  {
  protected :
    //! matrix to project on basis functions
    Matrix<Real_wp> jacobi_11_pol;
    Vector<Matrix<Real_wp> > jacobi_2ip1_pol;
    int type_basis;
    VectReal_wp CoefLeg11;
    Matrix<Real_wp> CoefJacobi;
    Matrix<int> NumDofsTri;
    
  public :
    enum { TENSOR_BASIS, INVARIANT_BASIS};
    
    TriangleHierarchic();

    void SetBasisType(int type);
    const Matrix<int>& GetNumDofs2D() const;
    
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    
  protected :
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    template<class Vector1>
    void ComputeProjectionDofGen(const Vector1& feval, Vector1& contrib) const;

  public :    
    void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                                  const ElementReference_Dim<Dimension2>& FaceCoarse,
                                  const ElementReference_Dim<Dimension2>& FaceFine) const;

    void ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& ) const;
    void ComputeGradientPhiRef(const R2& pointloc, VectR2& ) const;
    
    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
    friend ostream& operator <<(ostream& out, const TriangleHierarchic& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_HIERARCHIC_HXX
#endif
