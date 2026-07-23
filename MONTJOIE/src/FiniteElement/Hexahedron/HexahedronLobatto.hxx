#ifndef MONTJOIE_FILE_HEXAHEDRON_LOBATTO_HXX

namespace Montjoie
{
  
  //! nodal finite element with Gauss-Lobatto points
  /*!
    Basis functions are Lagrange functions associated with Gauss-Lobatto points
    phi_i^GL(x) phi_j^GL(y) phi_k^GL(z)
    Gauss-Lobatto points are also used as quadrature points
   */
  class HexahedronLobatto : public HexahedronGauss
  {
  protected:
    VectR3 points3d_all;
    
  public :
    HexahedronLobatto();
    
    const Vector<R3>& PointsQuadND() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected :
    template<class Vector1,class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2 & res, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const;

    template<class Vector1>
    void SolveMassMatrixGen(Vector1&) const;

    template<class Vector1>
    void SolveCholeskyGen(const SeldonTranspose& TransA, Vector1&) const;
    
    template<class Vector2>
    void MltMassMatrixGen(Vector2&) const;

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
    void ApplyRhSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&, Vector2&) const;

    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;

    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
			 Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhGen(const Vector1&, Vector2& ) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhTransposeGen(const Vector1&, Vector2& ) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&, Vector2&) const;
    
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
    void GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const;    
    void GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const;

    void GetValueSinglePhiQuadrature(int, VectReal_wp&) const;
    void GetGradientSinglePhiQuadrature(int, VectReal_wp&, VectR3&) const;

    void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval, VectReal_wp& res, int num_loc) const;
    void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval, VectComplex_wp& res, int num_loc) const;

    void ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    void ComputeGradientBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

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

    void ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp&) const;
    void ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
		      VectComplex_wp& Vy, VectComplex_wp& Vz) const;

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

    void ApplyConstantRh(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRh(const VectComplex_wp&, VectComplex_wp&) const;
    
    void ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp&) const;
    void ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
			      VectComplex_wp& Vy, VectComplex_wp&) const;
    
    void ApplyConstantRhTranspose(const VectReal_wp&, VectReal_wp&) const;
    void ApplyConstantRhTranspose(const VectComplex_wp&, VectComplex_wp&) const;

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
  
  ostream& operator <<(ostream& out, const HexahedronLobatto& e);

} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_LOBATTO_HXX
#endif
