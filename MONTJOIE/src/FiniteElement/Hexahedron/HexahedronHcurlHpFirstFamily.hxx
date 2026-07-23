#ifndef MONTJOIE_FILE_HEXAHEDRON_HCURL_HP_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! finite element class for Nedelec's first family on hexahedra
  /*!
    Hierarchical basis functions are used (see Bergot thesis)
    and generate the space Q_{r-1, r, r} x Q_{r, r-1, r} x Q_{r, r, r-1}
    They are based on Jacobi polynomials
    Quadrature rules are Gauss points
   */
  class HexahedronHcurlHpFirstFamily : public HexahedronReference<2>
  {
    
  protected :
    Array3D<int> NumDofsX, NumDofsY, NumDofsZ;
    
    Matrix<Real_wp> LegendrePolynom, JacobiPolynom;    
    Matrix<Real_wp> ValLeg, ValBubbleJac, DerBubbleJac;
    Matrix<Real_wp, General, RowSparse> ch1_node, ch2_node, ch3_node,
      ch1_loc, ch2_loc, ch3_loc, rh_loc;
    
    TinyVector<Matrix<Real_wp, General, RowSparse>, 6> sh_loc, sh_locX, sh_locY,
      sh_curl_loc, sh_curl_locX, sh_curl_locY, sh_nodeX, sh_nodeY, sh_curl_nodeX, sh_curl_nodeY;
    
    TinyMatrix<Matrix<Real_wp, General, ArrayRowSparse>, General, 3, 3>
    sparse_mass_matrix, sparse_stiff_matrix;
    
    Matrix<Real_wp, Symmetric, RowSymSparse> sparse_matrix_chol;
    
  public :
    HexahedronHcurlHpFirstFamily();
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf_tri = 0, int rsurf_quad = 0,
                                int type_surf_tri = -1, int type_surf_quad = -1, int gauss_z = -1);
    
    //void CancelHighOrderDofs(IVect&, const IVect& re, const IVect& rf, int ri);
    
  protected:
    void ConstructFunctions();
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    void ComputeProjectionMatrix(Matrix<Real_wp, General, ArrayRowSparse>& Ch1,
                                 Matrix<Real_wp, General, ArrayRowSparse>& Ch2,
                                 Matrix<Real_wp, General, ArrayRowSparse>& Ch3,
                                 const Array3D<int>& Num, const VectReal_wp& points,
                                 Matrix<Real_wp>& ValLegendre, Matrix<Real_wp>& ValBubble);
    
    void ComputeSurfaceProjection(TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shX,
                                  TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shY,
                                  const VectReal_wp& points, const Matrix<Real_wp>& ValLegendre,
                                  const Matrix<Real_wp>& ValBubble, const Matrix<int>& Num);
    
    void ComputeSurfaceCurlProjection(TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shX,
				      TinyVector<Matrix<Real_wp, General, RowSparse>, 6>& shY,
				      const VectReal_wp& points,
				      const Matrix<Real_wp>& DerBubble,
				      const Matrix<Real_wp>& ValLegendre,
				      const Matrix<Real_wp>& ValBubble, const Matrix<int>& Num);
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;
    
    template<class Vector1, class Vector2>
    void ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const;
            
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

  public:
    void ComputeValuesPhiRef(const R3& ptloc, VectR3& phi) const;
    void ComputeCurlPhiRef(const R3& ptloc, VectR3& dphi) const;
    
    void GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const;    
    void GetCurlPhiOnQuadraturePoint(int k, VectR3& grad_phi) const;

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;

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
        
  };
  
  ostream& operator <<(ostream& out, const HexahedronHcurlHpFirstFamily& e);
  
} // namespace Montjoie

#define MONTJOIE_FILE_HEXAHEDRON_HCURL_HP_FIRST_FAMILY_HXX
#endif

