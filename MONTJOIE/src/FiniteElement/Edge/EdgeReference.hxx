#ifndef MONTJOIE_FILE_EDGE_REFERENCE_HXX

namespace Montjoie
{
  
  //! Base class for finite element in 1-D
  template<>
  class ElementGeomReference<Dimension1>
  {
    friend class ElementReference<Dimension1, 1>;
    friend class ElementReference_Dim<Dimension1>;
  protected :
    //! informations to print
    int print_level;
    //! order of approximation (r)
    int order;
    //! order of quadrature
    int order_quad;
    //! order of approximation for geometry
    int order_geom;
    //! number of quadrature points
    int nb_points;
    //! number of dofs per edge
    int nb_dof_loc;
    //! dof points
    VectReal_wp points_dof;
    //! true for finite elements achieving mass lumping (diagonal mass matrix)
    bool mass_lumping;

    //! stiffness_matrix(i,j) = \f$ \int \hat{\varphi}_i'(x) \hat{\varphi}_j'(x) dx  \f$
    //! on the unit edge [0,1]
    Matrix<Real_wp, Symmetric, RowSymPacked> stiffness_matrix, mass_matrix;
    
    //! gradient_matrix(i,j) = \f$ \int \hat{\varphi}_j'(x) \hat{\varphi}_i(x) dx  \f$
    //! on the unit edge [0,1]
    Matrix<Real_wp> gradient_matrix;

  public :
    
    VectReal_wp Weights; //!< integration weights
    VectReal_wp Points; //!< integration points
    
    ElementGeomReference<Dimension1>();
    virtual ~ElementGeomReference<Dimension1>();

    bool LumpedMassMatrix() const;
    int GetHybridType() const;
    
    size_t GetMemorySize() const;
    int GetOrder() const;
    int GetGeometryOrder() const;
    const VectReal_wp& PointsDof() const;

  protected:
    // construction of quadrature formulas
    void ConstructQuadrature(int r, int type_quadrature =
			     Globatto<Real_wp>::QUADRATURE_GAUSS);

    virtual void SetDofPoints(const VectReal_wp&) = 0;
    
  public:
    int GetNbPointsQuadratureInside() const;
    int GetNbDof() const;
    
    int GetNbPointsNodalElt() const;

    FiniteElementProjector* GetNewNodalInterpolation() const;

    virtual void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
					int type_func = -1) = 0;
    
    virtual void GetPolynomialFunctions(Vector<UnivariatePolynomial<Real_wp> >& Phi,
					Vector<UnivariatePolynomial<Real_wp> >& dPhi) const;

    virtual void ComputeValuesPhiRef(const Real_wp& pointloc, VectReal_wp& res) const = 0;
    virtual void ComputeGradientPhiRef(const Real_wp& pointloc, VectReal_wp& res) const = 0;

    virtual Real_wp GetValuePhi1D(int i, const Real_wp& pointloc) const = 0;
    virtual Real_wp GetGradientPhi1D(int i, const Real_wp& pointloc) const = 0;
    virtual void GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const = 0;
    
    void ComputeValuesPhiNodalRef(const Real_wp& x, VectReal_wp& res) const;
    
    // x = \f$ F_j(\hat{x}) \f$
    // where \f$ \hat{x} \f$ is defined on the unit edge [0,1]
    void FjElem(const VectReal_wp& s, VectReal_wp& points) const; 
    void FjElemDof(const VectReal_wp& s, VectReal_wp& points) const; 

    const Matrix<Real_wp, Symmetric, RowSymPacked>& GetStiffnessMatrix() const;
    const Matrix<Real_wp, Symmetric, RowSymPacked>& GetMassMatrix() const;
    const Matrix<Real_wp>& GetGradientMatrix() const;
    
    virtual void ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const = 0;
    virtual void ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const = 0;

    virtual void ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const = 0;
    virtual void ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const = 0;

    virtual void ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const = 0;
    virtual void ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const = 0;

    // compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i'(x) dx \f$
    // feval(j) is the evalution of f on the quadrature point j
    virtual void ComputeIntegralGradientRef(const VectReal_wp& feval, VectReal_wp& contrib) const = 0;
    virtual void ComputeIntegralGradientRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const = 0;

    virtual void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const = 0;
    virtual void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const = 0;

    virtual void AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C, const Real_wp& D,
				       const Real_wp& E, const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const = 0;

    virtual void AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C, const Complex_wp& D,
				       const Complex_wp& E, const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const = 0;

    virtual void AddVariableElemMatrix(int m, int n, const VectReal_wp& mass, const VectReal_wp& C, const VectReal_wp& D,
				       const VectReal_wp& E, const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const = 0;

    virtual void AddVariableElemMatrix(int m, int n, const VectComplex_wp& mass, const VectComplex_wp& C, const VectComplex_wp& D,
				       const VectComplex_wp& E, const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const = 0;
    
    // to display object ElementGeomReference<Dimension1>
    friend ostream& operator<<(ostream& out, const ElementGeomReference<Dimension1>& e);

  };


  //! Gauss-Lobatto finite element in 1-D with Gauss quadrature
  class EdgeGauss : public ElementGeomReference<Dimension1>
  {
    
  protected :
    //! value of basis functions on quadrature points
    Matrix<Real_wp> Value_Phi;
    //! value of gradient of basis functions on quadrature points
    Matrix<Real_wp> Gradient_Phi;
    //! object defining basis functions
    Globatto<Real_wp> lob_geom, lob_basis;
    
  public :
    enum{ GAUSS, LOBATTO, LOBATTO_INT, RADAU};
    
    EdgeGauss();

    size_t GetMemorySize() const;
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int type_func = -1);
    
  protected :
    void SetDofPoints(const VectReal_wp&);
    
    // computation of stiffness matrix
    void ConstructFunctions(int type_func = -1);
    void ConstructStiffnessMatrix();

    template<class Vector1>
    void ComputeIntegralGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ApplyChGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ApplyChTransposeGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ComputeIntegralGradientGen(const Vector1& feval, Vector1& contrib) const;
    
    template<class T0, class Matrix1>
    void AddConstantElemMatrixGen(int m, int n, const T0& mass, const T0& C,
				  const T0& D, const T0& E,
				  const TinyVector<bool, 4>& null_term, Matrix1& A) const;
        
    template<class T, class Matrix1>
    void AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
				  const Vector<T>& C, const Vector<T>& D, const Vector<T>& E,
				  const TinyVector<bool, 4>& null_term, Matrix1& mat) const;
    
  public :
    // how to number the 1D-mesh
    void ConstructNumberMap(NumberMap& map, int) const;
    
    // compute res(i) = \f$ \hat{\varphi}_i(\mbox{pointloc}) \f$
    void ComputeValuesPhiRef(const Real_wp& pointloc, VectReal_wp& res) const;
    void ComputeGradientPhiRef(const Real_wp& pointloc, VectReal_wp& res) const;
    void ComputeValuesPhiNodalRef(const Real_wp& pointloc, VectReal_wp& res) const;
    void GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const;
    
    Real_wp GetValuePhi1D(int i, const Real_wp& pointloc) const;
    Real_wp GetGradientPhi1D(int i, const Real_wp& pointloc) const;
    void GetPolynomialFunctions(Vector<UnivariatePolynomial<Real_wp> >& Phi,
				Vector<UnivariatePolynomial<Real_wp> >& dPhi) const;
    
    // compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i(x) dx \f$
    // feval(j) is the evalution of f on the quadrature point j
    void ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    // compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i'(x) dx \f$
    // feval(j) is the evalution of f on the quadrature point j
    void ComputeIntegralGradientRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeIntegralGradientRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C,
			       const Real_wp& D, const Real_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C,
			       const Complex_wp& D, const Complex_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			       const VectReal_wp& C, const VectReal_wp& D, const VectReal_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const;

    void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			       const VectComplex_wp& C, const VectComplex_wp& D, const VectComplex_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const;
    
    friend ostream& operator<<(ostream& out, const EdgeGauss& e);
  };
  
  
  //! Gauss-Lobatto finite element in 1-D with Lobatto quadrature
  class EdgeLobatto : public EdgeGauss
  {
  private:
    template<class Vector1>
    void ComputeIntegralGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ApplyChGen(const Vector1& feval, Vector1& contrib) const;

  public :
    EdgeLobatto();
    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int type_func = -1);

    void GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const;
    
    // compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i(x) dx \f$
    // feval(j) is the evalution of f on the quadrature point j
    void ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

  };


  //! Gauss-Lobatto finite element in 1-D with Radau quadrature
  class EdgeRadau : public EdgeLobatto
  {
  public:
    EdgeRadau();

    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int type_func = -1);
    
  };

  
  //! Hierarchical basis functions over interval [0,1]
  class EdgeHierarchic : public ElementGeomReference<Dimension1>
  {
    
  protected :
    //! value of basis functions on quadrature points
    Matrix<Real_wp> Value_Phi;
    //! value of gradient of basis functions on quadrature points
    Matrix<Real_wp> Gradient_Phi;
    Matrix<Real_wp> jacobi_11_pol;
    VectReal_wp CoefLeg11;
    Matrix<Real_wp> legendre_pol;
    LegendrePolynomial<Real_wp> Leg_pol;
    int type_function;
    
  public :
    enum {LEGENDRE, JACOBI_11, LEGENDRE_COMBINED, SHEN}; 

    EdgeHierarchic();

    size_t GetMemorySize() const;
      
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int type_func = -1);

  protected :
    void SetDofPoints(const VectReal_wp&);
    
    void ConstructFunctions(int type_func);

    template<class Vector1>
    void ComputeIntegralGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ApplyChGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ApplyChTransposeGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void ComputeIntegralGradientGen(const Vector1& feval, Vector1& contrib) const;

    template<class Vector1>
    void SolveMassMatrix(Vector1& contrib) const;

  public :
    // how to number the 1D-mesh
    void ConstructNumberMap(NumberMap& map, int) const;

    Real_wp GetValuePhi1D(int i, const Real_wp& pointloc) const;
    Real_wp GetGradientPhi1D(int i, const Real_wp& pointloc) const;
    
    // compute res(i) = \f$ \hat{\varphi}_i(\mbox{pointloc}) \f$
    void ComputeValuesPhiRef(const Real_wp& pointloc, VectReal_wp& res) const;
    void ComputeValuesPhiNodalRef(const Real_wp& pointloc, VectReal_wp& res) const;
    void ComputeGradientPhiRef(const Real_wp& pointloc, VectReal_wp& res) const;
    void GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const;
    
    // compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i(x) dx \f$
    // feval(j) is the evalution of f on the quadrature point j
    void ComputeIntegralRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeIntegralRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    // compute contrib(i) = \f$ \int_{[0,1]} f(x) \hat{\varphi}_i'(x) dx \f$
    // feval(j) is the evalution of f on the quadrature point j
    void ComputeIntegralGradientRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeIntegralGradientRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;
    
    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C,
			       const Real_wp& D, const Real_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C,
			       const Complex_wp& D, const Complex_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			       const VectReal_wp& C, const VectReal_wp& D, const VectReal_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const;

    void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			       const VectComplex_wp& C, const VectComplex_wp& D, const VectComplex_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const;
    
    friend ostream& operator<<(ostream& out, const EdgeHierarchic& e);
  };

  
  template<>
  class ElementReference_Dim<Dimension1> : public ElementReference_Base
  {
  protected :
    ElementGeomReference<Dimension1>& elt_geom;
    
  public :
    ElementReference_Dim(ElementGeomReference<Dimension1>& elt);
    virtual ~ElementReference_Dim();
    
    inline int GetHybridType() const { return 0; }
    int GetGeometryOrder() const;
    size_t GetMemorySize() const;

    inline const ElementGeomReference<Dimension1>& GetGeometricElement() const { return elt_geom; }
    
    void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
				int type_func = -1);

    void SetDofPoints(const VectReal_wp&);
    const VectReal_wp& PointsDof() const;
    
    void ComputeValuesPhiRef(const Real_wp& pointloc, VectReal_wp& res) const;
    Real_wp GetValuePhi1D(int i, const Real_wp& pointloc) const;
    Real_wp GetGradientPhi1D(int i, const Real_wp& pointloc) const;
    void ComputeGradientPhiRef(const Real_wp& pointloc, VectReal_wp& res) const;

    void GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const;
    
    const VectReal_wp& Points() const;
    const VectReal_wp& Weights() const;

    const Real_wp& Points(int i) const;
    const Real_wp& Weights(int i) const;
    const Real_wp& WeightsND(int i) const;
    const VectReal_wp& WeightsND() const;

    inline int GetNbPointsNodalElt() const { return elt_geom.PointsDof().GetM(); }
    inline const Real_wp& PointsNodalND(int i) const { return elt_geom.PointsDof()(i); }
    
    const Matrix<Real_wp, Symmetric, RowSymPacked>& GetStiffnessMatrix() const;
    const Matrix<Real_wp>& GetGradientMatrix() const;

    void GetPolynomialFunctions(Vector<UnivariatePolynomial<Real_wp> >& Phi,
				Vector<UnivariatePolynomial<Real_wp> >& dPhi) const;

    void FjElem(const VectReal_wp& s, VectReal_wp& points) const;

    void ComputeIntegralRef(const VectReal_wp& u, VectReal_wp& v) const;
    void ComputeIntegralGradientRef(const VectReal_wp& u, VectReal_wp& v) const;

    void ComputeIntegralRef(const VectComplex_wp& u, VectComplex_wp& v) const;
    void ComputeIntegralGradientRef(const VectComplex_wp& u, VectComplex_wp& v) const;

    void ConstructNumberMap(NumberMap& nmap, int) const;
    
    void ApplyCh(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyCh(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ApplyChTranspose(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ApplyChTranspose(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void SolveMassMatrix(VectReal_wp&) const;
    void SolveMassMatrix(VectComplex_wp& x) const;

    void SolveCholesky(const SeldonTranspose& TransA, VectReal_wp&) const;
    void SolveCholesky(const SeldonTranspose& TransA, VectComplex_wp&) const;

    void MltMassMatrix(VectReal_wp&) const;
    void MltMassMatrix(VectComplex_wp&) const;

    void ComputeMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
			   const VectReal_wp & coef) const;

    void IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
			     const VectReal_wp & coef) const;

    FiniteElementProjector* GetNewNodalInterpolation() const;
    void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
				  const ElementReference_Dim<Dimension1>& FaceCoarse,
				  const ElementReference_Dim<Dimension1>& FaceFine) const;
    
    void ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const;
    void ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const;

    void ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const;
    void ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const;

    void AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C,
			       const Real_wp& D, const Real_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const;
    
    void AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C,
			       const Complex_wp& D, const Complex_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const;
    
    void AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			       const VectReal_wp& C, const VectReal_wp& D, const VectReal_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const;

    void AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			       const VectComplex_wp& C, const VectComplex_wp& D, const VectComplex_wp& E,
			       const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const;
    
  };

  template<>
  class ElementReference<Dimension1, 1> : public ElementReference_Dim<Dimension1>
  {
  public:
    inline ElementReference(ElementGeomReference<Dimension1>& elt) :
      ElementReference_Dim<Dimension1>(elt) {}
    
    typedef ElementReference<Dimension1, 1> ElementReferenceTrace;
  };

  
  class EdgeLobattoReference : public ElementReference<Dimension1, 1>
  {
    EdgeLobatto edge;
    
  public :
    EdgeLobattoReference();
    
  };

  class EdgeRadauReference : public ElementReference<Dimension1, 1>
  {
    EdgeRadau edge;
    
  public :
    EdgeRadauReference();
    
  };

  
  class EdgeGaussReference : public ElementReference<Dimension1, 1>
  {
    EdgeGauss edge;
    
  public :
    EdgeGaussReference();

  };


  class EdgeHierarchicReference : public ElementReference<Dimension1, 1>
  {
    EdgeHierarchic edge;

  public :
    EdgeHierarchicReference();
    
  };
  
}

#define MONTJOIE_FILE_EDGE_REFERENCE_HXX
#endif


