#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_LEGENDRE_HXX

namespace Montjoie
{
  
  //! nodal finite element using Legendre polynomials and DG method
  /*!
    Basis functions are product of Legendre polynomials :
    L_m(2x-1)  L_n(2y-1)  L_k(2z-1)   with  m+n+k <= r
    where L_i are the classical Legendre polynomials
    quadrature points are Gauss points
   */
  class HexahedronDgLegendre : public HexahedronReference<1>
  {
  protected :
    Matrix<Real_wp> LegendrePolynom; VectReal_wp CoefLegendre;
    VectReal_wp InvWeightPolynomial;
    Array3D<int> NumOrtho3D; Matrix<int> CoordinateDofs;
    Matrix<Real_wp, General, RowSparse> rh_const, rh_constX, rh_constY, rh_constZ;
    Matrix<Real_wp, General, RowSparse> sh1_loc, sh2_loc, sh3_loc;
    Matrix<Real_wp, General, RowSparse> ch1_loc, ch2_loc, ch3_loc;
    Matrix<Real_wp> ValLi, ValDLi, Li_Geom, dLi_Geom;
    Matrix<Real_wp> CoefDLi, stiff1D;
    Vector<Matrix<int> > NumProjOperator;

    Vector<Matrix<Real_wp> > ProjOperatorDerivQuadOrder;
    
  public :
    HexahedronDgLegendre();
   
#ifdef MONTJOIE_FAST_HEXA 
    int GetNbPointsUsedForSh() const;
    int GetNbPointsUsedForFaceSh(int num_loc) const;
    int GetOffsetSh(int num_loc) const;
#endif

    size_t GetMemorySize() const;
    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);

  protected :

    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;
    
    template<class Vector1,class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2& res, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeValueBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;
    
    template<class Vector1>
    void SolveMassMatrixGen(Vector1& x) const;
    
    template<class Vector1>
    void MltMassMatrixGen(Vector1& y) const;

    template<class Vector1>
    void SolveCholeskyGen(const SeldonTranspose& TransA, Vector1&) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&, Vector2&) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhGen(const Vector1&, Vector2& ) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&, Vector2&) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhTransposeGen(const Vector1&, Vector2& ) const;
    
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;

    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
			 Vector2& Vh, int r = 0) const;
    
    template<class T>
    void AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const;
    
    template<class T, class Prop>
    void AddConstantElemMatrixGen(int m, int n, const T& mass, const TinyMatrix<T, Prop, 3, 3>& C,
				  const TinyVector<T, 3>& D, const TinyVector<T, 3>& E, 
				  const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const;

    template<class T>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<T>& A, VirtualMatrix<T>& mat) const;

    template<class T, class Prop>
    void AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& A,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
				  const Vector<TinyVector<T, 3> >& D,
				  const Vector<TinyVector<T, 3> >& E,
				  const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const;
    
  public :
    void ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                                   const ElementReference<Dimension2, 1>&,
                                                   const Vector<VectR2>& Pts);
    
    void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const;
    void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const;
    
    void GetValueSinglePhiQuadrature(int, VectReal_wp&) const;
    void GetGradientSinglePhiQuadrature(int, VectReal_wp&, VectR3&) const;
    
    void GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const;
    void GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const;

    void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;
    
    void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const;
    void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const;

    void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;
    
    void ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    void ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;

    void SolveMassMatrix(VectReal_wp&) const;
    void SolveMassMatrix(VectComplex_wp& x) const;

    void SolveCholesky(const SeldonTranspose& TransA, VectReal_wp&) const;
    void SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp&) const;
    
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

    void ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp&) const;
    void ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
		      VectComplex_wp& Vy, VectComplex_wp& Vz) const;

    void ApplyConstantRh(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRh(const VectComplex_wp&, VectComplex_wp&) const;
    
    void ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp&) const;
    void ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
			      VectComplex_wp& Vy, VectComplex_wp&) const;
    
    void ApplyConstantRhTranspose(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRhTranspose(const VectComplex_wp&, VectComplex_wp&) const;

    void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		 VectComplex_wp& Vh, int r = 0) const;

    void ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
		      VectReal_wp& Vh, int r = 0) const;

    void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
		      VectComplex_wp& Vh, int r = 0) const;

    void AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const;
    void AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			       const TinyMatrix<Real_wp, General, 3, 3>& C,
			       const R3& D, const R3& E,
			       const TinyVector<bool, 4>& null_term,
			       VirtualMatrix<Real_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			       const TinyMatrix<Complex_wp, General, 3, 3>& C,
			       const R3_Complex_wp& D, const R3_Complex_wp& E,
			       const TinyVector<bool, 4>& null_term,
			       VirtualMatrix<Complex_wp>& A) const;
    
    void AddVariableMassMatrix(int off_row, int off_col,
			       const Vector<Real_wp>& A,
			       VirtualMatrix<Real_wp>& mat) const;
    
    void AddVariableMassMatrix(int off_row, int off_col,
			       const Vector<Complex_wp>& A,
			       VirtualMatrix<Complex_wp>& mat) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			       const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
			       const Vector<R3>& D, const Vector<R3>& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			       const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
			       const Vector<R3_Complex_wp>& D, const Vector<R3_Complex_wp>& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_DG_LEGENDRE_HXX
#endif
