#ifndef MONTJOIE_FILE_QUADRANGLE_HIERARCHIC_HXX

namespace Montjoie
{

  //! hp finite element on quadrilaterals
  /*!
    Hierarchical functions are constructed with Jacobi polynomials
    and generate the usual space Q_r
   */
  class QuadrangleHierarchic : public QuadrangleReference<1>
  {    
  protected :
    Matrix<int> NumDofsQuad;

    //! coefficients for Jacobi polynomials P_m^{1,1}
    Matrix<Real_wp> jacobi_11_pol;
    //! sparse stiffness matrix
    Matrix<Real_wp, General, RowSparse> rh_loc, ch1_loc, ch2_loc, const_rh;
    Matrix<Real_wp, Symmetric, RowSymSparse> MhLoc;
    Matrix<Real_wp> ShLoc;
    VectReal_wp InvWeightFct;    
    VectReal_wp CoefLegendre;

  public :
    QuadrangleHierarchic();
    
    const Matrix<int>& GetNumDofs2D() const;

    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    size_t GetMemorySize() const;

    // construction of finite element
    void ConstructFiniteElement(int r, int rquad = 0, int rgeom = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
    
  protected :
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();

    template<class Vector1>
    void ComputeProjectionDofGen(const Vector1& feval, Vector1& contrib) const;
        
    template<class Vector1>
    void MltMassMatrixGen(Vector1& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
        
  public :
    void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                                  const ElementReference_Dim<Dimension2>& FaceCoarse,
                                  const ElementReference_Dim<Dimension2>& FaceFine) const;

    void ComputeValuesPhiRef(const R2& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R2& x, VectR2& phi) const;

    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void MltMassMatrix(VectReal_wp&) const;
    void MltMassMatrix(VectComplex_wp&) const;

    void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhQuadrature(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyRhQuadrature(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyConstantRh(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRh(const VectComplex_wp&, VectComplex_wp&) const;

    void ApplyConstantRhTranspose(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRhTranspose(const VectComplex_wp&, VectComplex_wp&) const;

    void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r = 0) const;
        
    friend ostream& operator <<(ostream& out, const QuadrangleHierarchic& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HIERARCHIC_HXX
#endif

