#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! edge finite element of Nedelec's first family on quadrilateral elements
  /*!
    Gauss and Gauss-Lobatto points are used for dof definitions
    Gauss-Lobatto points are used quadrature formula for integration of mass matrix
    Gauss points are used for integration of stiffness matrix
   */
  class QuadrangleHcurlFirstFamily : public QuadrangleReference<2>
  {       
  protected :
    //! tensorial to scalar numbering for dofs orientated along X
    Matrix<int> NumDofs_X;
    //! tensorial to scalar numbering for dofs orientated along Y
    Matrix<int> NumDofs_Y;
    //! scalar to tensorial matching arry
    Matrix<int> CoordinateDofs;
    //! tensorial to scalar numbering for Gauss quadrature points 
    Matrix<int> NumNodes2D_gauss;
    //! On each dof, normale
    VectR2 normale_dof;
    //! \f$ \hat{\varphi}_i^{GL} ( \hat{\xi}_j^G) \f$
    Matrix<Real_wp> dGL_G, dGL_GL;
    //! \f$ \hat{\varphi}_i^{G} ( \hat{\xi}_j^GL) \f$
    Matrix<Real_wp> G_GL, GL_Gexact;
    //! \f$ \hat{\varphi}_i^{GL} ( \hat{\xi}_j^G) \f$
    Matrix<Real_wp> GL_G;
    //! \f$ \hat{\varphi}_i^{G}(0) \mbox{ and } \hat{\varphi}_i^{G}(0) \f$
    Matrix<Real_wp> ValGaussExtremity;
    //! sparse gradient operator
    Matrix<Real_wp, General, RowSparse> rh_loc;
    //! sparse interpolation operator
    Matrix<Real_wp, General, RowSparse> ch_loc;
    VectReal_wp WeightsGauss;
    
  public :
    QuadrangleHcurlFirstFamily();
    
    Real_wp WeightsGauss2D(int k) const;
    int GetXdofNumber(int i, int j) const;
    int GetYdofNumber(int i, int j) const;
    const Matrix<Real_wp>& GetGgl() const;
    const Matrix<Real_wp>& GetdGlg() const;
    int GetNumNodes2DGauss(int i, int j) const;
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                                int rsurf = 0, int type_surf = -1);

    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
  protected:
    void ConstructMassMatrix();
    void ConstructStiffnessMatrix();
    
    void ConstructFunctions();

    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;
    
    template<class Vector1, class Vector2>
    void ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const;

    template<class T>
    void ApplyShGen(const T& alpha, int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const;
    
    template<class T>
    void ApplyNablaShGen(const T& alpha, int num_loc, const Vector<T>& Uh, Vector<T>& Vh, int r) const;
    
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval,
                                        Vector2& res, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalCurlGen(const Vector1& Un, Vector2 & Unode) const;

    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class T, class Prop>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 2, 2> >& B,
				  VirtualMatrix<T>& mat) const;
    
    template<class T>
    void AddVariableStiffnessMatrixGen(int off_row, int off_col, const Vector<T>& A,
				       VirtualMatrix<T>& mat) const;
    
  public:
    void ComputeValuesPhiRef(const R2& x, VectR2& phi) const;
    void ComputeCurlPhiRef(const R2& x, VectReal_wp& phi) const;
    
    void ComputeValuesPhiH(const R2& point_loc, const R2& point_glob,
			   Vector<TinyVector<Real_wp, 1> >& res, const Matrix2_2& dfjm1,
			   const Mesh<Dimension2>& mesh, int nquad) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalValuesH(const Vector1& Un, Vector2 & Unode) const;

    template<class Vector1, class Vector2>
    void ApplyRhGauss(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGaussTranspose(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1>
    void InterpolateHQuadrature(const SetPoints<Dimension2>& PointsElem,
                                const SetMatrices<Dimension2>& MatricesElem,
                                const Vector1& Hloc, Vector1& Hloc_node,
                                const ElementReference<Dimension2, 2>& FaceBasis,
                                const Mesh<Dimension2>& mesh, int nquad) const;

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ProjectQuadratureToDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ProjectQuadratureToDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
    virtual void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			 VectReal_wp& Vh, int r = 0) const;

    virtual void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			 VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			      VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			      VectComplex_wp& Vh, int r = 0) const;

    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    virtual void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    virtual void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;

    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				      const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
				      VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
				       VirtualMatrix<Complex_wp>& mat) const;
    
    virtual void AddVariableStiffnessMatrix(int off_row, int off_col,
					    const Vector<Real_wp>& A, VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableStiffnessMatrix(int off_row, int off_col,
					    const Vector<Complex_wp>& A, VirtualMatrix<Complex_wp>& mat) const;
    
  };
  
  ostream& operator <<(ostream& out, const QuadrangleHcurlFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HCURL_FIRST_FAMILY_HXX
#endif
