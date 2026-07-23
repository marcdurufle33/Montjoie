#ifndef MONTJOIE_FILE_WEDGE_HIERARCHIC_HXX

namespace Montjoie
{
  //! hierarchical finite element for triangular prisms
  /*!
    Hierarchical basis functions are constructed with Jacobi polynomials
    and generate the polynomial space P_r(x, y) P_r(z)
   */
  class WedgeHierarchic : public WedgeReference<1>
  {
  protected :
    Matrix<Real_wp> jacobi_11_pol;
    Vector<Matrix<Real_wp> > jacobi_2ip1_pol;
    VectReal_wp CoefLeg11;
    Matrix<Real_wp> CoefJacobi;

    Array3D<int> NumDofs3D;
    Matrix<int> CoordinateDofs3D, NumDofs2D_tri, NumDofs2D_quad;
    
    Globatto<Real_wp> lob_quad, lob_y;
    VectReal_wp points1d_y, weights1d_y;
    
    Matrix<Real_wp> ShLoc, ShTri;
    Matrix<Real_wp, General, RowSparse> sh1_loc, sh2_loc;
    Matrix<Real_wp, General, RowSparse> ch1_loc, ch2_loc, ch3_loc, const_rh, rh_loc;
    Matrix<Real_wp, Symmetric, RowSymSparse> MhLoc;
    Vector<Matrix<int> > NumProjOperator;
    
    VectReal_wp DerivDxtildeDx, DerivDxtildeDy;
    TriangleHierarchic* function_basis_tri;
    QuadrangleHierarchic* function_basis_quad;
    
  public :
    WedgeHierarchic();
    
    void ConstructNumberMap(NumberMap& nmap, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);

  protected :    
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    template<class Vector1>
    void MltMassMatrixGen(Vector1& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
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
                                  const ElementReference_Dim<Dimension3>& FaceCoarse,
                                  const ElementReference_Dim<Dimension3>& FaceFine) const;

    void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const;
    
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
    
    void ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                                   const ElementReference<Dimension2, 1>& Fb,
                                                   const Vector<VectR2>& Pts);
    
    void MltMassMatrix(VectReal_wp&) const;
    void MltMassMatrix(VectComplex_wp&) const;

    void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyConstantRh(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyConstantRh(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyConstantRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyConstantRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		 VectReal_wp& Vh, int r = 0) const;
    
    void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r = 0) const;
        
    friend ostream& operator <<(ostream& out, const WedgeHierarchic& e);
    
  };

 } // end namespace

#define MONTJOIE_FILE_WEDGE_HIERARCHIC_HXX
#endif
