#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_LOBATTO_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's second family on hexahedron
  /*!
    Nodal basis functions are used basd on Gauss-Lobatto points
    phi_i^GL(x) phi_j^GL(y) phi_k^GL(z)  e_x
    phi_j^GL(x) phi_i^GL(y) phi_k^GL(z)  e_y
    phi_j^GL(x) phi_k^GL(y) phi_i^GL(z)  e_z
    with i <= r, j <= r, k <= r
    Gauss-Lobatto points are also used for integration, inducing mass lumping
   */
  class HexahedronHcurlLobatto : public HexahedronReference<2>
  {    
  protected:
    IVect Node_to_DofX; //!< for each node, dof along x
    IVect Node_to_DofY; //!< for each node, dof along y
    IVect Node_to_DofZ; //!< for each node, dof along z
    IVect Dof_to_Node;
    Seldon::Array3D<int> NumDofs_X; //!< dofs on hexahedron oriented along X
    Seldon::Array3D<int> NumDofs_Y; //!< dofs on hexahedron oriented along Y
    Seldon::Array3D<int> NumDofs_Z; //!< dofs on hexahedron oriented along Z
    
    Matrix<Real_wp, General, RowSparse> rh_loc; //!< stiffness matrix
    Matrix<Real_wp, General, ArrayRowSparse> rh_permut;
    Matrix<int> CoordinateDofs; //!< tensorial <-> scalar matching array
    Matrix<Real_wp> GL_GX, stiff1d;
    Matrix<Real_wp, General, RowSparse> ch1_node, ch2_node, ch3_node;
    VectReal_wp invWeightsMassG, invSqrtWeightsMassG;
    
  public :
    HexahedronHcurlLobatto();

    const Real_wp& dGL_GL(int i, int j) const;
    const Matrix<Real_wp>& GradPhi() const;
    inline const Matrix<Real_wp, General, ArrayRowSparse>& GetPermutedRh() const { return rh_permut; }
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
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
    void ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                        Vector2& res, int num_loc) const;
    
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
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
        
    template<class Vector1, class Vector2>
    void ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class Vector1, class Vector2>
    void ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r = 0) const;
    
    template<class T0, class Vector1, class Vector2>
    void ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh,
			 Vector2& Vh, int r = 0) const;
    
    template<class T, class Prop>
    void AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const;
    
    template<class T, class Prop>
    void AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const;

    template<class T, class Prop>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
				  VirtualMatrix<T>& mat) const;
    
  public:
    int GetCoordinateDof(int i) const;
    int GetPointNode_FromDofNumber(int i) const;
    void GetDofNumber_FromPointNode(int i, int& ix, int& iy, int& iz) const;
    
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;

    void GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const;    
    void GetCurlPhiOnQuadraturePoint(int k, VectR3& phi) const;

    template<class Vector1, class Vector2>
    void ComputeNodalValuesH(const Vector1& Un, Vector2 & Unode);

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;

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

  };
  
  ostream& operator <<(ostream& out, const HexahedronHcurlLobatto& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_LOBATTO_HXX
#endif

