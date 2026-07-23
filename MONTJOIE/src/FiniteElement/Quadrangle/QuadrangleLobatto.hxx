#ifndef MONTJOIE_FILE_QUADRANGLE_LOBATTO_HXX

namespace Montjoie
{

  //! Gauss-Lobatto quadrilateral finite element
  /*!
    Gauss-Lobatto points are used both for interpolation and quadrature
   */
  class QuadrangleLobatto : public QuadrangleGauss
  {
  protected :
    Matrix<Real_wp, General, RowSparse> rh_weight, rh_weightX, rh_weightY, rh_boundary;
    VectR2 points2d_all;
    
  public :
    QuadrangleLobatto();
    size_t GetMemorySize() const;

    // construction of finite element    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
  protected :
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    template<class Vector1>
    void SolveMassMatrixGen(Vector1&) const;
    
    template<class Vector1>
    void SolveCholeskyGen(const SeldonTranspose& TransA, Vector1&) const;
    
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
    void ApplyRhBoundaryGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhBoundaryTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

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
    void AddConstantElemMatrixGen(int m, int n, const T& mass, const TinyMatrix<T, Prop, 2, 2>& C, 
				  const TinyVector<T, 2>& D, const TinyVector<T, 2>& E, 
				  const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const;

    template<class T>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<T>& A, VirtualMatrix<T>& mat) const;

    template<class T, class Prop>
    void AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& A,
				  const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
				  const Vector<TinyVector<T, 2> >& D,
				  const Vector<TinyVector<T, 2> >& E,
				  const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const;
    
    template<class Vector1, class Vector2>
    void ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const;

  public :
    const VectR2& PointsQuadND() const;
    
    void GetValueSinglePhiQuadrature(int, VectReal_wp&) const;
    void GetGradientSinglePhiQuadrature(int, VectReal_wp&, VectR2&) const;

    void GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const;
    void GetGradientPhiOnQuadraturePoint(int k, VectR2& phi) const;

    void PickNearDofs(int pos, const VectBool& DofUsed, IVect& ListeDof, int nb_dof) const;    

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
    
    void ApplyRhBoundary(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyRhBoundary(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhBoundaryTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhBoundaryTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyConstantRh(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRh(const VectComplex_wp&, VectComplex_wp&) const;

    void ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const;
    void ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const;

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

    void AddConstantMassMatrix(int m, int n, const Real_wp& mass,
			       VirtualMatrix<Real_wp>& A) const;
    
    void AddConstantMassMatrix(int m, int n, const Complex_wp& mass,
			       VirtualMatrix<Complex_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Real_wp& mass,
			       const TinyMatrix<Real_wp, General, 2, 2>& C,
			       const R2& D, const R2& E,
			       const TinyVector<bool, 4>& null_term,
			       VirtualMatrix<Real_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
			       const TinyMatrix<Complex_wp, General, 2, 2>& C,
			       const R2_Complex_wp& D, const R2_Complex_wp& E,
			       const TinyVector<bool, 4>& null_term,
			       VirtualMatrix<Complex_wp>& A) const;
    
    void AddVariableMassMatrix(int off_row, int off_col,
			       const Vector<Real_wp>& A,
			       VirtualMatrix<Real_wp>& mat) const;
    
    void AddVariableMassMatrix(int off_row, int off_col,
			       const Vector<Complex_wp>& A,
			       VirtualMatrix<Complex_wp>& mat) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			       const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
			       const Vector<R2>& D, const Vector<R2>& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			       const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
			       const Vector<R2_Complex_wp>& D, const Vector<R2_Complex_wp>& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
    
    void ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    void ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;
    
    friend ostream& operator <<(ostream& out, const QuadrangleLobatto& e);
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_LOBATTO_HXX
#endif

