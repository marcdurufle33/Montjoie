#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_LOBATTO_HXX

namespace Montjoie
{
  
  //! edge finite element of Nedelec's second family on quadrilaterals
  /*!
    Use of Gauss-Lobatto points, to have a mass lumping 
    and to have more efficient computation of matrix 
    Solution of Maxwell's equations with this finite element
    exhibit spurious modes, the first family is more attractive
    since it is spurious-free.
   */
  class QuadrangleHcurlLobatto : public QuadrangleReference<2>
  {    
  protected :
    //! not used
    int nb_points_quadrature_gauss;
    //! tensorial to scalar numbering for dofs oriented along X
    Matrix<int> NumDofs_X;
    //! tensorial to scalar numbering for dofs oriented along Y
    Matrix<int> NumDofs_Y;
    //! scalar to tensorial matching array
    Matrix<int> CoordinateDofs;
    //! sparse gradient operator
    Matrix<Real_wp, General, RowSparse> rh_loc;
    //! for each nodes, dof numbers for x-dof
    IVect Node_to_DofX;
    //! for each nodes, dof numbers for y-dof
    IVect Node_to_DofY;
    //! sign between x-dof and y-dof on each node
    IVect Sign_CrossMassMatrix;
    //! Gauss-Lobatto weights
    VectReal_wp WeightGauss;
    Matrix<Real_wp> GL_Gexact;
    
  public :
    QuadrangleHcurlLobatto();
    
    int GetXdofNumber(int i) const;
    int GetYdofNumber(int i) const;
    int GetXdofNumber(int i, int j) const;
    int GetYdofNumber(int i, int j) const;
    const Real_wp& dGL_GL(int i, int j) const;
    
    void ConstructNumberMap(NumberMap& map, int dg) const;
    size_t GetMemorySize() const;
    
    // construction of finite element
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int rsurf = 0, int type_surf = -1);

    void ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad);
    
  protected:
    void ConstructFunctions();

    void ConstructStiffnessMatrix();
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;
    
    template<class Vector1, class Vector2>
    void ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const;

    template<class Vector1, class Vector2>
    void ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const;
    
    template<class Vector1, class Vector2>
    void ComputeGaussIntegralSurfaceGen(const Vector1& feval,
                                        Vector2& res, int num_loc) const;
    
    template<class T>
    void ApplyShGen(const T& alpha, int num_loc, const Vector<T>& feval, Vector<T>& res, int r) const;
    
    template<class T>
    void ApplyNablaShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const;

  public:
    void ComputeValuesPhiRef(const R2& x, VectR2& phi) const;
    void ComputeCurlPhiRef(const R2& x, VectReal_wp& phi) const;    

    template<class Vector1, class Vector2>
    void ComputeNodalValuesH(const Vector1& Un, Vector2 & Unode) const;

    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    virtual void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

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
    
  };

} // namespace Montjoie

#define MONTJOIE_FILE_QUADRANGLE_HCURL_LOBATTO_HXX
#endif

