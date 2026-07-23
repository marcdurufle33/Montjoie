#ifndef MONTJOIE_FILE_VOLUME_REFERENCE_HXX

#include "VolumeGeomReference.hxx"

namespace Montjoie
{
  
  //! base class for 3-D finite element
  template<int type>
  class VolumeReference : public ElementReference_Dim<Dimension3>,
			  public ElementReferenceType<Dimension3, type>
  {
  public :
    typedef Dimension3 Dimension; //!< dimension 3
    
  protected :
    
    int nb_dof_tri, nb_dof_quad;
    // local number of degrees of freedom on each face
    //! FacesDof(k, i) -> dof k on face i
    Matrix<int> FacesDof, FacesNode;
        
    //! projection operator for interpolation on different orders (for DG)
    Vector<Matrix<Real_wp> > ProjOperatorTriOrder, ProjOperatorQuadOrder;
    
  public :
    
    VolumeReference(ElementGeomReference<Dimension3>&);
 
    static void SetInputData(const string&, const Vector<string>&);

        
    /************************
     * Convenient functions *
     ************************/
    
    
    size_t GetMemorySize() const;        
    int GetLocalNumber(int j, int k) const;    
    int GetNbDofBoundary(int) const;
    
    virtual void ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const;
    virtual void ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const;
    
    virtual void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    virtual void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;


    /*****************
     * Other methods *
     *****************/
        
    void SetVariableOrder(const Mesh<Dimension3>& mesh,
			  const MeshNumbering<Dimension3>& mesh_num);

    virtual void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
					int rsurf_tri = 0, int rsurf_quad = 0, int type_surf_tri = -1,
					int type_surf_quad = -1, int gauss_z = -1) = 0;

    virtual void ComputeTriangularInterpolationProjectorOrder(const IVect& order_elt,
							      const Vector<VectR2>& Pts);
    
    virtual void ComputeQuadrangularInterpolationProjectorOrder(const IVect& order_elt,
								const Vector<VectR2>& Pts);

  protected:
    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;
    
    template<int t>
    friend ostream& operator <<(ostream& out, const VolumeReference<t>& e);
    
  };

  template<int t>
  ostream& operator <<(ostream& out, const VolumeReference<t>& e);
  
#ifdef MONTJOIE_WITH_NODAL_H1
  //! 3-D H^1 elements
  template<>
  class ElementReference<Dimension3, 1> : public VolumeReference<1>, public FiniteElementH1<Dimension3>
  {
  protected :

    ElementReference<Dimension2, 1>* element_tri_surf, *element_quad_surf;

  public:
    ElementReference(ElementGeomReference<Dimension3>&);
    virtual ~ElementReference();
    
    const ElementReference<Dimension2, 1>& GetTriangularSurfaceFiniteElement() const;
    const ElementReference<Dimension2, 1>& GetQuadrangularSurfaceFiniteElement() const;
    const ElementReference<Dimension2, 1>& GetSurfaceFiniteElement(int n) const;
    
    size_t GetMemorySize() const;
    
    virtual void ComputeTriangularInterpolationProjectorOrder(const IVect& order_elt,
							      const Vector<VectR2>& Pts);
    
    virtual void ComputeQuadrangularInterpolationProjectorOrder(const IVect& order_elt,
								const Vector<VectR2>& Pts);

  protected:
    virtual void ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                                         const ElementReference<Dimension2, 1>& Fb,
                                                         const Vector<VectR2>& Pts);
    
    virtual void ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                                           const ElementReference<Dimension2, 1>& Fb,
                                                           const Vector<VectR2>& Pts);
    
    template<class Vector1, class Vector2>
    void ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;

    void FindH1LinearCombinationRotation(NumberMap& nmap,
                                         const ElementReference<Dimension2, 1>&,
					 const ElementReference<Dimension2, 1>&,
                                         bool display_message = true) const;
    
    void FindH1SignEdge(NumberMap& nmap) const;

    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhQuadratureSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&, Vector2&) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Vh, Vector2& Uh) const;

    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhGen(const Vector1&, Vector2& ) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhTransposeGen(const Vector1&, Vector2& ) const;

    template<class Vector1, class Vector2>
    void ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const;
    
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
    
    template<class T>
    void AddVariableMassMatrixOpt(int off_row, int off_col,
                                  const Vector<T>& A, VirtualMatrix<T>& mat) const;

    template<class T, class Prop>
    void AddVariableElemMatrixOpt(int off_row, int off_col, const Vector<T>& A,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
				  const Vector<TinyVector<T, 3> >& D,
				  const Vector<TinyVector<T, 3> >& E,
				  const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const;

  public :    
    template<class Vector1,class Vector2>
    void ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const;

    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceGradientRef(const Vector1 & feval, Vector2& res, int n) const;
    
    virtual void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;
    
    virtual void ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeGradientBoundaryRef(const VectComplex_wp& Un,
					    VectComplex_wp& Unode, int num_loc) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    virtual void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    virtual void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
        
    virtual void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;
    
    virtual void ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy, VectReal_wp&) const;
    virtual void ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
			      VectComplex_wp& Vy, VectComplex_wp&) const;

    virtual void ApplyRhQuadratureSplit(const VectReal_wp& Uh, VectReal_wp& Vh,
					VectReal_wp&, VectReal_wp&) const;

    virtual void ApplyRhQuadratureSplit(const VectComplex_wp& Uh, VectComplex_wp& Vh,
					VectComplex_wp&, VectComplex_wp&) const;

    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyConstantRh(const VectReal_wp&, VectReal_wp&) const;
    virtual void ApplyConstantRh(const VectComplex_wp&, VectComplex_wp&) const;
    
    virtual void ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx,
				      VectReal_wp& Vy, VectReal_wp&) const;

    virtual void ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx,
				      VectComplex_wp& Vy, VectComplex_wp&) const;
    
    virtual void ApplyConstantRhTranspose(const VectReal_wp&, VectReal_wp&) const;
    virtual void ApplyConstantRhTranspose(const VectComplex_wp&, VectComplex_wp&) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const;
    virtual void AddConstantMassMatrix(int m, int n, const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const;
    
    virtual void AddConstantElemMatrix(int m, int n, const Real_wp& mass,
				       const TinyMatrix<Real_wp, General, 3, 3>& C,
				       const R3& D, const R3& E,
				       const TinyVector<bool, 4>& null_term,
				       VirtualMatrix<Real_wp>& A) const;

    virtual void AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
				       const TinyMatrix<Complex_wp, General, 3, 3>& C,
				       const R3_Complex_wp& D, const R3_Complex_wp& E,
				       const TinyVector<bool, 4>& null_term,
				       VirtualMatrix<Complex_wp>& A) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<Real_wp>& A,
				       VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<Complex_wp>& A,
				       VirtualMatrix<Complex_wp>& mat) const;
    
    virtual void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
				       const Vector<TinyMatrix<Real_wp, General, 3, 3> >& C,
				       const Vector<R3>& D, const Vector<R3>& E,
				       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;

    virtual void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
				       const Vector<TinyMatrix<Complex_wp, General, 3, 3> >& C,
				       const Vector<R3_Complex_wp>& D, const Vector<R3_Complex_wp>& E,
				       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
    
    void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                                  const ElementReference_Dim<Dimension3>& FaceCoarse,
                                  const ElementReference_Dim<Dimension3>& FaceFine) const;
  };
#endif


#ifdef MONTJOIE_WITH_NODAL_HCURL
  class VolumeHcurlReference;
  
  //! 3-D H(curl) elements
  template<>
  class ElementReference<Dimension3, 2> : public VolumeReference<2>, public FiniteElementHcurl<Dimension3>
  {
    friend class VolumeHcurlReference;
    
  protected:
    ElementReference<Dimension2, 2>* element_tri_surf, *element_quad_surf;
    ElementReference<Dimension3, 1>* element_scal_vol;
    
  public :
    ElementReference(ElementGeomReference<Dimension3>&);
    virtual ~ElementReference();

    const ElementReference<Dimension2, 2>& GetTriangularSurfaceFiniteElement() const;
    const ElementReference<Dimension2, 2>& GetQuadrangularSurfaceFiniteElement() const;
    const ElementReference<Dimension2, 2>& GetSurfaceFiniteElement(int n) const;
    const ElementReference<Dimension3, 1>& GetScalarElement() const;
    
    size_t GetMemorySize() const;
    
  protected:    
    void ConstructHcurlElementaryMatrix();
    
    void FindDofsOnFace(bool nodal_element);

    void FindHcurlLinearCombinationRotation(NumberMap& nmap,
                                            const ElementReference<Dimension2, 2>&,
					    const ElementReference<Dimension2, 2>&,
                                            bool take_all_dofs = false, bool display_message = true) const;

    void FindHcurlSignEdge(NumberMap& nmap) const;
    
    template<class Vector1, class Vector2>
    void ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeCurlBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Vh, Vector2& Uh) const;

    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class T, class Prop>
    void AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const;
    
    template<class T, class Prop>
    void AddConstantStiffnessMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& C, VirtualMatrix<T>& A) const;

    template<class T, class Prop>
    void AddVariableStiffnessMatrixGen(int off_row, int off_col,
				       const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
				       VirtualMatrix<T>& mat) const;
    
    template<class T, class Prop>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
				  VirtualMatrix<T>& mat) const;
    
  public:    
    template<class Vector1,class Vector2>
    void ComputeIntegralCurlRef(const Vector1 & feval, Vector2& res) const;

    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceCurlRef(const Vector1 & feval, Vector2& res, int n) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    virtual void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    virtual void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;

    virtual void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    virtual void ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeCurlBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    virtual void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

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


  //! base class linking H^1 elements with H(curl) elements
  class VolumeHcurlReference
  {
  protected:
    ElementReference<Dimension3, 1>& elt_H1;
    ElementReference<Dimension3, 2>& elt_Hcurl;

  public:
    VolumeHcurlReference(ElementReference<Dimension3, 1>& elt,
			 ElementReference<Dimension3, 2>& elt_hc);

    void ConstructFiniteElement();

    void ComputeValuesPhiRef(const R3& x, VectR3& phi) const;

    void GetValuePhiOnQuadraturePoint(int, VectR3& phi) const;
    
    template<class T>
    void ApplyChGen(const Vector<T>& U, Vector<T>& V) const;
    
  };
    
#endif


#ifdef MONTJOIE_WITH_NODAL_HDIV
  //! 3-D H(div) elements
  template<>
  class ElementReference<Dimension3, 3> : public VolumeReference<3>, public FiniteElementHdiv<Dimension3>
  {
  protected:
    ElementReference<Dimension2, 1>* element_tri_surf, *element_quad_surf;
    
  public :
    ElementReference(ElementGeomReference<Dimension3>&);
    virtual ~ElementReference();
    
    const ElementReference<Dimension2, 1>& GetTriangularSurfaceFiniteElement() const;
    const ElementReference<Dimension2, 1>& GetQuadrangularSurfaceFiniteElement() const;
    const ElementReference<Dimension2, 1>& GetSurfaceFiniteElement(int n) const;
    
    size_t GetMemorySize() const;
    
  protected:
    void ConstructHdivElementaryMatrix();
    void FindDofsOnFace();

    void FindHdivLinearCombinationRotation(NumberMap& nmap,
					   const VectR2& PointsTri, const VectR2& PointsQuad,
					   const Matrix<Real_wp>& ValuePhiTri,
					   const Matrix<Real_wp>& ValuePhiQuad,
					   bool display_message = true) const;
    
    template<class Vector1, class Vector2>
    void ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;

    template<class Vector1, class Vector2>
    void ComputeDivBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChGen(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const;
    
    template<class T>
    void ApplyShGen(const T& alpha, int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const;

    template<class T>
    void ApplyNablaShGen(const T& alpha, int num_loc, const Vector<T>& Uh, Vector<T>& Vh, int r) const;

    template<class T, class Prop>
    void AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 3, 3>& mass, VirtualMatrix<T>& A) const;

    template<class T>
    void AddConstantStiffnessMatrixGen(int m, int n, const T& C, VirtualMatrix<T>& A) const;
    
  public:
    template<class Vector1,class Vector2>
    void ComputeIntegralDivRef(const Vector1 & feval, Vector2& res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceDivRef(const Vector1 & feval, Vector2& res, int n) const;
    
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    virtual void ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const;
    virtual void ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const;
        
    virtual void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;

    virtual void ComputeDivBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeDivBoundaryRef(const VectComplex_wp& Un,
				       VectComplex_wp& Unode, int num_loc) const;

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

    virtual void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			      VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			      VectComplex_wp& Vh, int r = 0) const;

    virtual void AddConstantStiffnessMatrix(int m, int n, const Real_wp& mass,
					    VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantStiffnessMatrix(int m, int n, const Complex_wp& mass,
					    VirtualMatrix<Complex_wp>& A) const;

    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 3, 3>& C,
				       VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 3, 3>& C,
				       VirtualMatrix<Complex_wp>& A) const;

  };
#endif
  
} // end namespace

#define MONTJOIE_FILE_VOLUME_REFERENCE_HXX
#endif
