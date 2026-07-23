#ifndef MONTJOIE_FILE_QUADRANGLE_GAUSS_HXX

namespace Montjoie
{
  
  //! class for nodal finite element with Gauss-Lobatto points
  //! but with integration on Gauss points
  /*!
    This class is particularly useful for the class QuadrangleDgGauss
    otherwise it has been used directly to solve Helmholtz equation
    in order to evaluate the advantage to use more accurate integration formulas
    And as a result, it is clearly better to use QuadrangleLobatto and not QuadrangleGauss
    because the improvement is quite neglectible while the overcost in time is large
   */
  class QuadrangleGauss : public QuadrangleReference<1>
  {
  protected :
    
    // arrays to interpolate from Gauss (G) to Gauss-Lobatto points (GL)
    Matrix<Real_wp> GL_G, invGL_G, GL_Gquad; //!< interpolation from Gauss-Lobatto points to Gauss points
    Matrix<Real_wp> G_GL; //!< interpolation from Gauss points to Gauss-Lobatto points
    Matrix<Real_wp> dGL_G; //!< derivate from Gauss-Lobatto points to Gauss points
    Matrix<Real_wp> GL_Geom, dGL_Geom;
    Matrix<int> NumFct2D, CoordinateDofs;
    VectReal_wp invWeights2d, invSqrtWeights2d;
    Matrix<Real_wp> stiffness_matrix1D, mass_matrix1D, gradient_matrix1D;    
    TensorizedProjector<Dimension2> proj_dof_to_nodal;
    //! sparse stiffness matrix
    Matrix<Real_wp, General, RowSparse> rh_loc, ch1_node, ch2_node, rh_locX, rh_locY;
    Matrix<Real_wp, General, RowSparse> ch1_loc, ch2_loc;
    
  public :
    static int type_quadrature_default;
    
    QuadrangleGauss();

    size_t GetMemorySize() const;
    void SetRadauInterpolationPoints();
    int GetNumDofs(int i, int j) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);
  protected :
    void ConstructFunctions();
    void ConstructMassMatrix();
    
    template<class Vector1,class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2& res, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;

    template<class Vector1, class Vector2>
    void ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;
    
    template<class Vector1>
    void SolveMassMatrixGen(Vector1&) const;
    
    template<class Vector1>
    void SolveCholeskyGen(const SeldonTranspose& TransA, Vector1&) const;
    
    template<class Vector2>
    void MltMassMatrixGen(Vector2& Vh) const;    

    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&) const;

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
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
    
  public :
    void ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& ) const;
    void ComputeGradientPhiRef(const R2& pointloc, VectR2& ) const;

    void GetValueSinglePhiQuadrature(int, VectReal_wp&) const;
    void GetGradientSinglePhiQuadrature(int, VectReal_wp&, VectR2&) const;
    
    void GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const;
    void GetGradientPhiOnQuadraturePoint(int k, VectR2& phi) const;

    void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const;
    void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const;

    void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;

    void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;
    
    void ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    void ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    void SolveCholesky(const SeldonTranspose& TransA, VectReal_wp&) const;
    void SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp&) const;
    
    void MltMassMatrix(VectReal_wp&) const;
    void MltMassMatrix(VectComplex_wp&) const;

    void SolveMassMatrix(VectReal_wp&) const;
    void SolveMassMatrix(VectComplex_wp&) const;

    void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    void ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const;
    void ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const;

    void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;
    
    void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh, VectComplex_wp& Vh, int r = 0) const;

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
    
  };
  
}

#define MONTJOIE_FILE_QUADRANGLE_GAUSS_HXX
#endif
