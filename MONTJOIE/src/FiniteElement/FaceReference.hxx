#ifndef MONTJOIE_FILE_FACE_REFERENCE_HXX

#include "FaceGeomReference.hxx"

namespace Montjoie
{
  // declaration of finite element classes in 2-D (triangular or quadrangular)
  // only nodal finite elements (in a H1 or DG formulation)
  
  //! Base class for 2-D finite element
  /*!
    this class contains common structures to triangles and quadrangles
    quadrature points, nodal points, computation of curvatures
    definition of methods in file Face_Reference.cxx
  */
  template<int type>
  class FaceReference : public ElementReference_Dim<Dimension2>,
			public ElementReferenceType<Dimension2, type>
  {
  protected :    
    //! numbers of the dof points on the edges of the element
    Matrix<int> EdgesDof, EdgesNodal;

    ElementReference<Dimension1, 1>* element_surface;

  public :        
    //! projection operator for interpolation on different orders (for DG)
    Vector<Matrix<Real_wp> > ProjOperatorOrder;
    
    FaceReference(ElementGeomReference<Dimension2>&);
    virtual ~FaceReference();
    
    static void SetInputData(const string&, const Vector<string>&);

    
    /**********************
     * Convenient methods *
     **********************/
    
    size_t GetMemorySize() const;    
    int GetLocalNumber(int j, int k) const;
    int GetNbDofBoundary(int n) const;

    const ElementReference<Dimension1, 1>& GetSurfaceFiniteElement(int n) const;
    
    virtual void ProjectQuadratureToDofRef(const VectReal_wp& Equad, VectReal_wp& Edof) const;
    virtual void ProjectQuadratureToDofRef(const VectComplex_wp& Equad, VectComplex_wp& Edof) const;

    virtual void ComputeNodalValuesRef(const VectReal_wp& Un, VectReal_wp& Unode) const;
    virtual void ComputeNodalValuesRef(const VectComplex_wp& Un, VectComplex_wp& Unode) const;

    
    /*****************
     * Other methods *
     *****************/

    virtual void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
					int rsurf = 0, int type_surf = -1) = 0;
                    
    void SetVariableOrder(const Mesh<Dimension2>& mesh,
			  const MeshNumbering<Dimension2>& mesh_num);
    
    void ComputeVariableInterpolation(int nb_order,
				      Vector<bool>& order_present);

  protected:
    virtual void ComputeInterpolationProjectorOrder(const IVect& order_elt,
						    const ElementReference<Dimension1, 1>& Fb,
						    const Vector<VectReal_wp>& Pts);
    
    template<class Vector1, class Vector2>
    void ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const;

    template<int t>
    friend ostream& operator <<(ostream& out, const FaceReference<t>& e);
    
  };

  template<int t>
  ostream& operator <<(ostream& out, const FaceReference<t>& e);

  
#ifdef MONTJOIE_WITH_NODAL_H1
  //! 2-D H^1 elements
  template<>
  class ElementReference<Dimension2, 1> : public FaceReference<1>, public FiniteElementH1<Dimension2>
  {
  protected :
    template<class Vector1, class Vector2>
    void ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;
    
    template<class Vector1, class Vector2>
    void ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const;

    template<class Vector1, class Vector2>
    void ApplyRhGen(const Vector1& Vh, Vector2& Uh) const;
    
    template<class Vector1, class Vector2>
    void ApplyRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy) const;

    template<class Vector1, class Vector2>
    void ApplyRhQuadratureSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&) const;

    template<class Vector1, class Vector2>
    void ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const;

    template<class T1>
    void ApplyChGen(const Vector<T1>& Vh, Vector<T1>& Uh) const;

    template<class T1>
    void ApplyChTransposeGen(const Vector<T1>& Uh, Vector<T1>& Vh) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhGen(const Vector1&, Vector2& ) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy) const;
    
    template<class Vector1, class Vector2>
    void ApplyConstantRhTransposeGen(const Vector1&, Vector2& ) const;
    
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
    
    template<class T>
    void AddVariableMassMatrixOpt(int off_row, int off_col,
                                  const Vector<T>& A, VirtualMatrix<T>& mat) const;

    template<class T, class Prop>
    void AddVariableElemMatrixOpt(int off_row, int off_col, const Vector<T>& A,
				  const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
				  const Vector<TinyVector<T, 2> >& D,
				  const Vector<TinyVector<T, 2> >& E,
				  const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const;
    
  public :
    ElementReference(ElementGeomReference<Dimension2>&);

    template<class Vector1,class Vector2>
    void ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const;
    
    template<class Vector1,class Vector2>
    void ComputeIntegralSurfaceGradientRef(const Vector1 & feval, Vector2& res, int n) const;
    
  protected:
    virtual void ComputeInterpolationProjectorOrder(const IVect& order_elt,
						    const ElementReference<Dimension1, 1>& Fb,
						    const Vector<VectReal_wp>& Pts);

  public:
    virtual void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;
    
    virtual void ComputeGradientBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeGradientBoundaryRef(const VectComplex_wp& Un,
					    VectComplex_wp& Unode, int num_loc) const;
        
    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    virtual void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;
    
    virtual void ApplyRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const;
    virtual void ApplyRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const;

    virtual void ApplyRhQuadratureSplit(const VectReal_wp& Uh, VectReal_wp& Vh, VectReal_wp&) const;
    virtual void ApplyRhQuadratureSplit(const VectComplex_wp& Uh, VectComplex_wp& Vh, VectComplex_wp&) const;

    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;
    
    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyConstantRh(const VectReal_wp&, VectReal_wp&) const;
    virtual void ApplyConstantRh(const VectComplex_wp&, VectComplex_wp&) const;
    
    virtual void ApplyConstantRhSplit(const VectReal_wp& Uh, VectReal_wp& Vx, VectReal_wp& Vy) const;
    virtual void ApplyConstantRhSplit(const VectComplex_wp& Uh, VectComplex_wp& Vx, VectComplex_wp& Vy) const;
    
    virtual void ApplyConstantRhTranspose(const VectReal_wp&, VectReal_wp&) const;
    virtual void ApplyConstantRhTranspose(const VectComplex_wp&, VectComplex_wp&) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const Real_wp& mass,
				       VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const Complex_wp& mass,
				       VirtualMatrix<Complex_wp>& A) const;
    
    virtual void AddConstantElemMatrix(int m, int n, const Real_wp& mass,
				       const TinyMatrix<Real_wp, General, 2, 2>& C,
				       const R2& D, const R2& E,
				       const TinyVector<bool, 4>& null_term,
				       VirtualMatrix<Real_wp>& A) const;

    virtual void AddConstantElemMatrix(int m, int n, const Complex_wp& mass,
				       const TinyMatrix<Complex_wp, General, 2, 2>& C,
				       const R2_Complex_wp& D, const R2_Complex_wp& E,
				       const TinyVector<bool, 4>& null_term,
				       VirtualMatrix<Complex_wp>& A) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<Real_wp>& A,
				       VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<Complex_wp>& A,
 				       VirtualMatrix<Complex_wp>& mat) const;
    
    virtual void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
				       const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
				       const Vector<R2>& D, const Vector<R2>& E,
				       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;

    virtual void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
				       const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
				       const Vector<R2_Complex_wp>& D, const Vector<R2_Complex_wp>& E,
				       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;

    void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
                                  const ElementReference_Dim<Dimension2>& FaceCoarse,
                                  const ElementReference_Dim<Dimension2>& FaceFine) const;
    
    static void FindH1RotationTri(int rf, const VectR2& Points2D, const VectReal_wp& Weights2D,
                                  const Matrix<Real_wp>& ValuePhi, NumberMap& nmap);    
    
  };
#endif

  
#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! 2-D H(curl) elements
  template<>
  class ElementReference<Dimension2, 2> : public FaceReference<2>, public FiniteElementHcurl<Dimension2>
  {
  protected:
    int nb_dof_H;
    Matrix<Real_wp> Ch_Hloc;
    
  public :
    ElementReference(ElementGeomReference<Dimension2>&);

    int GetNbDofH() const;
    void FindHcurlRotationQuad(NumberMap& nmap, int offset = -1) const;
    void FindHcurlRotationTri(NumberMap& nmap, int offset = -1) const;
    
  protected:
    void ConstructHcurlElementaryMatrix();
    void FindDofsOnEdge();

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

    template<class T>
    void ApplyShGen(const T& alpha, int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const;

    template<class T>
    void ApplyShTransposeGen(int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const;

    template<class T>
    void ApplyNablaShGen(const T& alpha, int num_loc, const Vector<T>& Uh, Vector<T>& Vh, int r) const;
    
    template<class T>
    void ApplyNablaShTransposeGen(int num_loc, const Vector<T>& Uh, Vector<T>& Vh, int r) const;
    
    template<class T, class Prop>
    void AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 2, 2>& mass, VirtualMatrix<T>& A) const;
    
    template<class T>
    void AddConstantStiffnessMatrixGen(int m, int n, const T& C, VirtualMatrix<T>& A) const;
    
    template<class T>
    void AddVariableStiffnessMatrixGen(int off_row, int off_col,
				       const Vector<T>& C, VirtualMatrix<T>& mat) const;

    template<class T, class Prop>
    void AddVariableMassMatrixGen(int off_row, int off_col,
				  const Vector<TinyMatrix<T, Prop, 2, 2> >& A,
				  VirtualMatrix<T>& mat) const;
    
  public:

    inline void GetValueSinglePhiQuadrature(int, VectReal_wp&) const {}
    
    template<class Vector1, class Vector2>
    void ComputeIntegralCurlRef(const Vector1 & feval, Vector2& res) const;
    
    template<class Vector1, class Vector2>
    void ComputeIntegralSurfaceCurlRef(const Vector1 & feval, Vector2& res, int n) const;

    virtual void ComputeValuesPhiHRef(const R2& point_loc, VectReal_wp& res) const;

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;  
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    virtual void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    virtual void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    virtual void ComputeValueBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeValueBoundaryRef(const VectComplex_wp& Un, VectComplex_wp& Unode, int num_loc) const;
    
    virtual void ComputeCurlBoundaryRef(const VectReal_wp& Un, VectReal_wp& Unode, int num_loc) const;
    virtual void ComputeCurlBoundaryRef(const VectComplex_wp& Un,
					VectComplex_wp& Unode, int num_loc) const;
    
    virtual void ApplyRh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyRh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyRhTranspose(const VectReal_wp& Uh, VectReal_wp& Vh) const;
    virtual void ApplyRhTranspose(const VectComplex_wp& Uh, VectComplex_wp& Vh) const;

    virtual void ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChH(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChH(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;

    virtual void ApplyChTransposeH(const VectReal_wp& Vh, VectReal_wp& Uh) const;
    virtual void ApplyChTransposeH(const VectComplex_wp& Vh, VectComplex_wp& Uh) const;
    
    virtual void ApplyShTranspose(int num_loc, const VectReal_wp& Uh, VectReal_wp& Vh, int r = 0) const;
    virtual void ApplyShTranspose(int num_loc, const VectComplex_wp& Uh,
				  VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplySh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			 VectReal_wp& Vh, int r = 0) const;

    virtual void ApplySh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			 VectComplex_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShTranspose(int num_loc, const VectReal_wp& Uh,
				       VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaShTranspose(int num_loc, const VectComplex_wp& Uh,
				       VectComplex_wp& Vh, int r = 0) const;
    
    virtual void ApplyNablaSh(const Real_wp& alpha, int num_loc, const VectReal_wp& Uh,
			      VectReal_wp& Vh, int r = 0) const;

    virtual void ApplyNablaSh(const Complex_wp& alpha, int num_loc, const VectComplex_wp& Uh,
			      VectComplex_wp& Vh, int r = 0) const;

    virtual void AddConstantStiffnessMatrix(int m, int n, const Real_wp& mass, VirtualMatrix<Real_wp>& A) const;
    virtual void AddConstantStiffnessMatrix(int m, int n,
					    const Complex_wp& mass, VirtualMatrix<Complex_wp>& A) const;

    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 2, 2>& C,
				       VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 2, 2>& C,
				       VirtualMatrix<Complex_wp>& A) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
                                       const Vector<TinyMatrix<Real_wp, General, 2, 2> >& C,
                                       VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableMassMatrix(int off_row, int off_col,
				       const Vector<TinyMatrix<Complex_wp, General, 2, 2> >& C,
				       VirtualMatrix<Complex_wp>& mat) const;
    
    virtual void AddVariableStiffnessMatrix(int off_row, int off_col,
					    const Vector<Real_wp>& A,
					    VirtualMatrix<Real_wp>& mat) const;

    virtual void AddVariableStiffnessMatrix(int off_row, int off_col,
					    const Vector<Complex_wp>& A,
					    VirtualMatrix<Complex_wp>& mat) const;

  };
#endif


#ifdef MONTJOIE_WITH_NODAL_HDIV
  //! 2-D H(div) elements
  template<>
  class ElementReference<Dimension2, 3> : public FaceReference<3>, public FiniteElementHdiv<Dimension2>
  {
  public :
    ElementReference(ElementGeomReference<Dimension2>&);

    static void FindHdivRotationQuad(int rf, const Matrix<int>& NumQuad2D, 
				     const Matrix<Real_wp>&, NumberMap& nmap);

    static void FindHdivRotationTri(int rf, const VectR2& Points2D, 
				    const VectReal_wp& Weights2D,
				    const Matrix<Real_wp>&, NumberMap& nmap);
    
  protected:    
    void ConstructHdivElementaryMatrix();
    void FindDofsOnEdge();

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
    void AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 2, 2>& mass, VirtualMatrix<T>& A) const;
    
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

    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Real_wp, General, 2, 2>& C,
				       VirtualMatrix<Real_wp>& A) const;
    
    virtual void AddConstantMassMatrix(int m, int n, const TinyMatrix<Complex_wp, General, 2, 2>& C,
				       VirtualMatrix<Complex_wp>& A) const;
    
  };
#endif
  
} // end namespace


#define MONTJOIE_FILE_FACE_REFERENCE_HXX
#endif
