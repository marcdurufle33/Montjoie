#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on hexhaedra
  /*!
    Nodal basis functions are used basd on Gauss and Gauss-Lobatto points
    phi_i^G(x) phi_j^GL(y) phi_k^GL(z)  e_x
    phi_j^GL(x) phi_i^G(y) phi_k^GL(z)  e_y
    phi_j^GL(x) phi_k^GL(y) phi_i^G(z)  e_z
    with i <= r-1, j <= r, k <= r
    Gauss-Lobatto points are used for integration so that Ch is very sparse
   */
  class HexahedronHcurlFirstFamily : public HexahedronReference<2>
  {
  protected :
    
    Globatto<Real_wp> gauss; //!< Gauss formula
    //! tensorial <-> scalar matching array for dofs orientated along x
    Seldon::Array3D<int> NumDofs_X;
    //! tensorial <-> scalar matching array for dofs orientated along y 
    Seldon::Array3D<int> NumDofs_Y;
    //! tensorial <-> scalar matching array for dofs orientated along z
    Seldon::Array3D<int> NumDofs_Z;
    Matrix<int> CoordinateDofs; //!< tensorial <-> scalar matching array for dofs
    
    // G = Gauss, GL = Gauss-Lobatto
    Matrix<Real_wp> G_GL; 
    Matrix<Real_wp> GL_G; 
    Matrix<Real_wp> dGL_G;
    Matrix<Real_wp> dGL_GL, stiff1d;

    Matrix<Real_wp, General, RowSparse> rh_loc; //!< stiffness matrix
    Matrix<Real_wp, General, RowSparse> ch_loc; //!< mass matrix
    
    //!< stiffness and mass matrices
    Matrix<Real_wp, General, RowSparse> ch1_loc; //!< mass matrix
    Matrix<Real_wp, General, RowSparse> ch2_loc; //!< mass matrix
    Matrix<Real_wp, General, RowSparse> ch3_loc; //!< mass matrix
    Matrix<Real_wp, General, RowSparse> rh_exact; //!< stiffness matrix
    Matrix<Real_wp, General, RowSparse> ch1_node, ch2_node, ch3_node;
    
    //! interpolation from r Gauss points to r+1 Gauss points
    Matrix<Real_wp> G_GX; 
    //! interpolation from r+1 Gauss-Lobatto points to r+1 Gauss points
    Matrix<Real_wp> GL_GX;
    //! derivative from r+1 Gauss points to r+1 Gauss points
    Matrix<Real_wp> dGX_GX; 
    //! integration weights for r+1 Gauss points
    VectReal_wp weights3d_Exact; 
    //! r+1 Gauss points
    Globatto<Real_wp> gauss_exact; 
    int nb_points_quadrature_gauss; //!< number of points for H
    VectReal_wp invWeightsMassG, invSqrtWeightsMassG;
    
  public :
    
    HexahedronHcurlFirstFamily();

    const Matrix<int>& GetCoordinateDofs() const;
    const Array3D<int>& GetNumDofsX() const;
    const Array3D<int>& GetNumDofsY() const;
    const Array3D<int>& GetNumDofsZ() const;
    const Matrix<Real_wp>& GetGradGL_GL() const;
    const Matrix<Real_wp>& GetGradGL_G() const;
    const Matrix<Real_wp>& GetValGauss_GL() const;
        
    size_t GetMemorySize() const;

    void ConstructNumberMap(NumberMap& map, int dg) const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
  protected:
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();        

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& res) const;
    
    template<class Vector1, class Vector2>
    void ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const;

    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;

    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2& res, int num_loc) const;

    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
                      Vector2& Vh, int r = 0) const;
                
    template<class Vector1,class Vector2>
    void ComputeNodalValuesGen(const Vector1& u_loc, Vector2& u_boundary) const;
    
    template<class Vector1,class Vector2>
    void ComputeValueBoundaryGen(const Vector1& u_loc,
				  Vector2& u_boundary, int num_loc) const;
    
    template<class Vector1,class Vector2>
    void ComputeCurlBoundaryGen(const Vector1& u_loc,
				Vector2& u_boundary, int num_loc) const;
    
    template<class Vector1>
    void SolveMassMatrixGen(Vector1&) const;
    
    template<class Vector1>
    void SolveCholeskyGen(const SeldonTranspose& TransA, Vector1&) const;
    
    template<class Vector2>
    void MltMassMatrixGen(Vector2&) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class T, class Prop>
    void AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const;
    
    template<class T, class Prop>
    void AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const;

    template<class T, class Prop>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& B,
				  VirtualMatrix<T>& mat) const;
    
    template<class T, class Prop>
    void AddVariableStiffnessMatrixGen(int off_row, int off_col,
				       const Vector<TinyMatrix<T, Prop, 3, 3> >& B,
				       VirtualMatrix<T>& mat) const;
    
  public:
    int GetCoordinateDof(int i) const;
    int GetPointNode_FromDofNumber(int i) const;
    
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
    
    void GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const;    
    void GetCurlPhiOnQuadraturePoint(int k, VectR3& grad_phi) const;

    void ComputeLocalProlongationLowOrder(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
					  const ElementReference_Dim<Dimension3>& FaceCoarse) const;
    
    void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
				  const ElementReference_Dim<Dimension3>& FaceCoarse,
				  const ElementReference_Dim<Dimension3>& FaceFine) const;
    
    void PickNearDofs(int pos, const VectBool& DofUsed, IVect& ListeDof, int nb_dof) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const;
    virtual void ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const;

    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    virtual void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    virtual void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;

    virtual void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    virtual void ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    virtual void SolveMassMatrix(VectReal_wp&) const;
    virtual void SolveMassMatrix(VectComplex_wp& x) const;
    
    virtual void SolveCholesky(const SeldonTranspose& TransA, VectReal_wp&) const;
    virtual void SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp&) const;

    virtual void MltMassMatrix(VectReal_wp&) const;
    virtual void MltMassMatrix(VectComplex_wp&) const;

    virtual void ApplyRhQuadrature(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhQuadrature(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ApplyRhQuadratureTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhQuadratureTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			 VectReal_wp& Vh, int r = 0) const;

    virtual void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			 VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    virtual void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
				  VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			      VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			      VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
				       VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
				       VectComplex_wp& Vh, int r = 0) const;

    virtual void AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& mass,
					    VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantStiffnessMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& mass,
					    VirtualMatrix<Complex_wp>& A) const;

    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
				       VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
				       VirtualMatrix<Complex_wp>& A) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				      const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
				       VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
				       VirtualMatrix<Complex_wp>& mat) const;

    virtual void AddVariableStiffnessMatrix(int off_row, int off_col,
					    const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
					    VirtualMatrix<Real_wp>& mat) const;
    
    virtual void AddVariableStiffnessMatrix(int off_row, int off_col,
					    const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
					    VirtualMatrix<Complex_wp>& mat) const;
    
  };
  
  ostream& operator <<(ostream& out, const HexahedronHcurlFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_FIRST_FAMILY_HXX
#endif
