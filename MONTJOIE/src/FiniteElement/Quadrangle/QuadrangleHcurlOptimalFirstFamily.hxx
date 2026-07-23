#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_FIRST_FAMILY_HXX

namespace Montjoie
{
  
  //! edge finite element of optimal Nedelec's first family on quadrilateral elements
  /*!
    Gauss-Lobatto points are used for dof definitions and quadrature
    Basis functions are equal to :
    phi_i^GLint(x) phi_j^GL(y)  e_x   i <= r-1, j <= r+1 
    phi_j^GL(x) phi_i^GLint(y)  e_y   i <= r-1, j <= r+1 
    where phi_i^GL are Lagrange functions associated with Gauss-Lobatto points
    and phi_i^GLint are associated with interior Gauss-Lobatto points
    (the two extremities are removed)
  */
  class QuadrangleHcurlOptimalFirstFamily : public QuadrangleReference<2>
  {    
  protected:
    //! tensorial to scalar numbering for dofs orientated along X
    Matrix<int> NumDofs_X;
    //! tensorial to scalar numbering for dofs orientated along Y
    Matrix<int> NumDofs_Y;

    //! scalar to tensorial matching arry
    Matrix<int> CoordinateDofs;
    //! On each dof, normale
    VectR2 normale_dof;
    //! \f$ \hat{\varphi}_i^{GL} ( \hat{\xi}_j^G) \f$
    Matrix<Real_wp> dGL_GL;
    //! \f$ \hat{\varphi}_i^{G} ( \hat{\xi}_j^GL) \f$
    Matrix<Real_wp> GL_Gexact;
    //! \f$ \hat{\varphi}_i^{G}(0) \mbox{ and } \hat{\varphi}_i^{G}(0) \f$
    Matrix<Real_wp> ValGaussExtremity;
    //! sparse gradient operator
    Matrix<Real_wp, General, RowSparse> rh_loc;
    //! sparse interpolation operator
    Matrix<Real_wp, General, RowSparse> ch_loc;
    bool project_with_quadrature;

  public :
    QuadrangleHcurlOptimalFirstFamily(bool proj = true);

    int GetXdofNumber(int i, int j) const;
    int GetYdofNumber(int i, int j) const;
    
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
    
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval,
                                        Vector2& res, int num_loc) const;
    
    template<class T>
    void ApplyShGen(const T& alpha, const Vector<T>& feval, Vector<T>& res, int num_loc, int r) const;
    
    template<class T>
    void ApplyNablaShGen(const T& alpha, const Vector<T>& feval, Vector<T>& res, int n, int r) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;
    
    template<class Vector1, class Vector2>
    void ComputeNodalCurlGen(const Vector1& Un, Vector2 & Unode) const;

    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

  public:
    void ComputeValuesPhiRef(const R2& x, VectR2& phi) const;
    void ComputeCurlPhiRef(const R2& x, VectReal_wp& phi) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    virtual void ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const;
    virtual void ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const;

    virtual void ComputeGaussIntegralSurfaceRef(const VectReal_wp & feval,
						VectReal_wp& res, int num_loc) const;
    
    virtual void ComputeGaussIntegralSurfaceRef(const VectComplex_wp & feval,
						VectComplex_wp& res, int num_loc) const;

    virtual void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			 VectReal_wp& Vh, int r = 0) const;

    virtual void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			 VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			      VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			      VectComplex_wp& Vh, int r = 0) const;

    virtual void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    virtual void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;

    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
  };
  
  ostream& operator <<(ostream& out, const QuadrangleHcurlOptimalFirstFamily& e);

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_FIRST_FAMILY_HXX
#endif
