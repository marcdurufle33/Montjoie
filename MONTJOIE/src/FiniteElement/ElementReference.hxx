#ifndef MONTJOIE_FILE_ELEMENT_REFERENCE_HXX

#include "ElementGeomReference.hxx"

namespace Montjoie
{

  template<class Dimension>
  class ElementReference_Dim;

  template<class Dimension, int type>
  class ElementReference { public: enum{nb_components_u=1, nb_components_grad=Dimension::dim_N};
    typedef ElementReference<Dimension, type> ElementReferenceTrace; };
  
  template<class Dim, int type>
  class ElementReferenceType
  { };

#ifdef MONTJOIE_WITH_NODAL_H1  
  //! scalar finite element
  template<class Dim>
  class ElementReferenceType<Dim, 1>
  {
  protected :
    typedef typename Dim::R_N R_N;
    typedef typename Dim::VectR_N VectR_N;
    
    //! Value_Phi(i,j) = \f$ \hat{\varphi}_i(\hat{\xi}_j) \f$
    Matrix<Real_wp> Value_Phi;
    //! Gradient_Phi(i,j) = \f$ \nabla \hat{\varphi}_i(\hat{\xi}_j) \f$
    Matrix<R_N> Gradient_Phi;

    //! values of basis functions on nodal points
    Matrix<Real_wp> ValuePhi_Nodal;
    //! gradient of basis functions on nodal points
    Matrix<R_N> GradPhi_Nodal;
    
    //! \int_K \grad phi_i phi_j
    Matrix<Real_wp> const_grad_matrix;

    //! \int_K \grad phi_i grad phi_j
    Matrix<Real_wp, Symmetric, RowSymPacked> const_stiff_matrix;

    template<class T1>
    void ProjectQuadratureToDofGen(const ElementReference<Dim, 1>& Fb,
				   const Vector<T1>& Equad, Vector<T1>& Edof) const;
    
  public :
    size_t GetMemorySize() const;
    
    virtual void GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const;
    virtual void GetGradientPhiOnQuadraturePoint(int k, VectR_N& phi) const;

    virtual void GetValueSinglePhiQuadrature(int, VectReal_wp&) const;
    virtual void GetGradientSinglePhiQuadrature(int, VectReal_wp&, VectR_N&) const;
    
    virtual void ComputeValuesPhiRef(const R_N&, VectReal_wp&) const = 0;
    virtual void ComputeGradientPhiRef(const R_N&, VectR_N&) const = 0;

    inline const Matrix<Real_wp>& GetValuePhi() { return Value_Phi; }

  protected :
    // computation of stiffness matrix, Value_Phi, Gradient_Phi
    void ConstructElementaryMatrix(ElementReference<Dim, 1>& Fb);
    
  };
#endif


#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! H(curl) finite element in 2-D
  template<>
  class ElementReferenceType<Dimension2, 2>
  {
  protected :
    //! \int \varphi_i^m \varphi_j^n for edge elements
    TinyMatrix<Matrix<Real_wp>, General, 2, 2> const_mass_matrix;
    //! \int curl(\varphi_i)^m curl(\varphi_j)^n for edge elements
    Matrix<Real_wp> const_curl_matrix;

    //! values of basis functions on quadrature points
    Matrix<R2> Value_PhiVec;
    //! values of basis functions on nodal points
    Matrix<Real_wp> ValuePhi_Nodal;
    //! curl of basis functions on nodal points (for edge elements)
    Matrix<Real_wp> CurlPhi_Nodal;
    //! curl of basis functions on quadrature points (2-D edge elements)
    Matrix<Real_wp> Curl_Phi2D;
    Vector<R2> tangente_dof;
    
    template<class T1>
    void ProjectQuadratureToDofGen(const ElementReference<Dimension2, 2>& Fb,
				   const Vector<T1>& Equad, Vector<T1>& Edof) const;
    
  public :
    size_t GetMemorySize() const;
    
    const Vector<R2>& GetOrientationDofs() const;
    
    virtual void GetValuePhiOnQuadraturePoint(int k, VectR2& phi) const;
    virtual void GetCurlPhiOnQuadraturePoint(int k, VectReal_wp& phi) const;

    virtual void ComputeValuesPhiRef(const R2&, VectR2&) const = 0;
    virtual void ComputeCurlPhiRef(const R2&, VectReal_wp&) const = 0;
    
  };


#ifdef MONTJOIE_WITH_THREE_DIM
  //! H(curl) finite element in 3-D
  template<>
  class ElementReferenceType<Dimension3, 2>
  {
  protected : 
    //! \int \varphi_i^m \varphi_j^n for edge elements
    TinyMatrix<Matrix<Real_wp>, General, 3, 3> const_mass_matrix;
    //! \int curl(\varphi_i)^m curl(\varphi_j)^n for edge elements
    TinyMatrix<Matrix<Real_wp>, General, 3, 3> const_curl_matrix;

    //! values of basis functions on nodal points
    Matrix<Real_wp> ValuePhi_Nodal;
    //! curl of basis functions on nodal points (for edge elements)
    Matrix<Real_wp> CurlPhi_Nodal;
    
    //! values of basis functions on quadrature points (for edge/facet elements)
    Matrix<R3> Value_PhiVec;
    //! curl of basis functions on quadrature points (3-D edge elements)
    Matrix<R3> Curl_Phi;

    Vector<R3> tangente_dof;

    template<class T1>
    void ProjectQuadratureToDofGen(const ElementReference<Dimension3, 2>& Fb,
				   const Vector<T1>& Equad, Vector<T1>& Edof) const;

  public :
    size_t GetMemorySize() const;

    virtual int GetCoordinateDof(int i) const;
    virtual int GetPointNode_FromDofNumber(int i) const;
    
    virtual void GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const;
    virtual void GetCurlPhiOnQuadraturePoint(int k, VectR3& phi) const;

    virtual void ComputeValuesPhiRef(const R3&, VectR3&) const = 0;
    virtual void ComputeCurlPhiRef(const R3&, VectR3&) const = 0;    
  };
#endif

#endif
  

#ifdef MONTJOIE_WITH_NODAL_HDIV
  //! H(div) finite element
  template<class Dim>
  class ElementReferenceType<Dim, 3>
  {
  protected : 
    typedef typename Dim::R_N R_N;
    typedef typename Dim::VectR_N VectR_N;
    typedef typename Dim::DimensionBoundary DimensionB;

    //! values of basis functions on nodal points
    Matrix<Real_wp> ValuePhi_Nodal;
    //! divergence of basis functions on nodal points (for facet elements)
    Matrix<Real_wp> DivPhi_Nodal;

    //! \int \varphi_i^m \varphi_j^n for facet elements
    TinyMatrix<Matrix<Real_wp>, General, Dim::dim_N, Dim::dim_N> const_mass_matrix;

    //! \int_K div phi_i div phi_j for facet elements
    Matrix<Real_wp, Symmetric, RowSymPacked> const_div_matrix;

    //! values of basis functions on quadrature points (for edge/facet elements)
    Matrix<R_N> Value_PhiVec;
 
    //! divergence of basis functions on quadrature points (facet elements)
    Matrix<Real_wp> Div_Phi;
    
    Vector<R_N> tangente_dof;

    ElementReference<DimensionB, 1>* element_surface, *element_tri_surf, *element_quad_surf;    
    
  protected :
    template<class T1>
    void ProjectQuadratureToDofGen(const ElementReference<Dim, 3>& Fb,
				   const Vector<T1>& Equad, Vector<T1>& Edof) const;
        
  public :
    size_t GetMemorySize() const;
    
    virtual void GetValuePhiOnQuadraturePoint(int k, VectR_N& phi) const;
    virtual void GetDivPhiOnQuadraturePoint(int k, VectReal_wp& phi) const;

    virtual void ComputeValuesPhiRef(const R_N&, VectR_N&) const = 0;
    virtual void ComputeDivPhiRef(const R_N&, VectReal_wp&) const = 0;
    
  };
#endif


  //! Generic class for finite element (with H1, Hcurl or H(div) transformation)
  template<class Dim>
  class ElementReference_Dim : public ElementReference_Base, public ElementGeomReferenceContainer<Dim>
  {
    template<class DimB, int t> friend class ElementReferenceType;

    typedef Dim Dimension;
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::VectR_N VectR_N;
    typedef typename Dimension::DimensionBoundary::VectR_N VectR_Nm1;
    
  public:
    ElementReference_Dim(ElementGeomReference<Dimension>&);    

    // pure virtual methods
    virtual void SetVariableOrder(const Mesh<Dimension>&, const MeshNumbering<Dimension>&) = 0;

    const VectReal_wp& WeightsND() const;
    const Real_wp& WeightsND(int j) const;
    
    // other methods
    VectR_N PointsQuadInsideND() const;
    
    void ComputeValuesPhiQuadratureRef(const R_N& pt, VectReal_wp& phi) const;
    
    void FjSurfaceElem(const VectR_N& s, SetPoints<Dimension>& res,
		       const Mesh<Dimension>&, int nquad, int ne) const;
    
    void DFjSurfaceElem(const VectR_N& s, const SetPoints<Dimension>& res,
			SetMatrices<Dimension>&,
			const Mesh<Dimension>&, int nquad, int ne) const;
    
    void FjSurfaceElemDof(const VectR_N& s, SetPoints<Dimension>& res,
			  const Mesh<Dimension>&, int nquad, int ne) const;
    
    void DFjSurfaceElemDof(const VectR_N& s, const SetPoints<Dimension>& res,
			   SetMatrices<Dimension>&,
			   const Mesh<Dimension>&, int nquad, int ne) const;

    virtual void ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
					  const ElementReference_Dim<Dimension>& FaceCoarse,
					  const ElementReference_Dim<Dimension>& FaceFine) const;

    virtual void ComputeLocalProlongationLowOrder(FiniteElementProjector& proj,
						  Matrix<Real_wp>& LocalProlongation,
						  const ElementReference_Dim<Dimension>& FaceCoarse) const;

    virtual void ComputeProjectionPointsSurf(int num_loc, const VectR_Nm1& ptsA,
					     const VectR_Nm1& ptsB, Matrix<Real_wp>& proj) const;
    
    virtual int GetLocalNumber(int j, int k) const = 0;
    virtual int GetNbDofBoundary(int n) const = 0;

  };


  void FillPositionDofBoundaries(const Matrix<int>&, Vector<int>& power_two, Vector<int>& pos);
  
} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_REFERENCE_HXX
#endif
