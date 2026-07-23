#ifndef MONTJOIE_FILE_BOUNDARY_INTEGRAL_EQUATION_HXX

namespace Montjoie
{

  //! base class for integral equations
  class BoundaryIntegralEquation_Base : public InputDataProblem_Base
  {
  protected :
    Real_wp omega, omega2, frequency;
    Vector<VectString> mesh_data;
    int default_order;
    int nb_unknowns_u, nodl;
    int print_level;
    int order_over_quadrature;
    int type_quadrature_element;
    int order_integration_singular;
    
  public:
    BoundaryIntegralEquation_Base();
    
    void SetInputData(const string& description_field, const VectString& parameters);

    void ConstructAll(const string& input_file, const string& name_elt);
    void ConstructAll(const Vector<string>& lines_data_file, const string& name_elt);

    virtual void InitIndices(int N) = 0;
    virtual void SetIndices(int, const Vector<string>&) = 0;
    virtual void SetPhysicalIndex(const string&, int ref, const Vector<string>&) = 0;

    void RunAll(const string& input_file, const string& name_elt);

    virtual void ComputeMeshAndFiniteElement(const string& name_elt) = 0;
    virtual void ComputeGeometryQuantities() = 0;
    virtual int GetNbBoundary() const = 0;

    void ComputeIntegralMatrix(Matrix<Complex_wp>& A) const;

    virtual void ComputeRightHandSide(Vector<Complex_wp>& b) const = 0;
    virtual void WriteSolution(const Vector<Complex_wp>& b) const = 0;
    virtual void ComputeElementaryMatrix(int, int, Vector<int>& row, Vector<int>& col,
					 Matrix<Complex_wp>& mat) const = 0;
    
  };


  //! class storing jacobian, normales on nodal/quadrature points
  template<class Dimension>
  class GeometryQuantityBEM
  {
    typedef typename Dimension::DimensionBoundary DimensionB;
    typedef typename Dimension::R_N R_N;
    
  protected:
    Vector<R_N> PointsNodal, NormaleNodal, PointsQuad, NormaleQuad;
    Vector<TinyMatrix<Real_wp, General, 3, 2> > MatDfNodal, MatDfQuad;
    VectReal_wp JacobNodal, JacobQuad;
    Vector<int> OffsetNodal, OffsetQuad;

  public:
    void ComputeGeometry(const Mesh<Dimension>& mesh,
			 const Vector<const ElementReference_Dim<DimensionB>* >& fem);

    inline const TinyMatrix<Real_wp, General, 3, 2>&
    GetDfNodal(int i, int j) const { return MatDfNodal(OffsetNodal(i) + j); }

    inline Real_wp
    GetJacobNodal(int i, int j) const { return JacobNodal(OffsetNodal(i) + j); }
    
    void FillJacobianMatQuad(int i, Vector<R_N>& PointsQuadI,
			     Vector<TinyMatrix<Real_wp, General, 3, 2> >& MatDfI) const;

    void FillNormaleQuad(int i, Vector<R_N>& PointsQuadI, Vector<R_N>& NormaleQuadI) const;

    void FillNormaleDfQuad(int i, Vector<R_N>& PointsQuadI, Vector<R_N>& NormaleQuadI,
			   Vector<TinyMatrix<Real_wp, General, 3, 2> >& MatDfI) const;
    
  };  
    
  //! class for integral equations for a given dimension
  template<class Dimension>
  class BoundaryIntegralEquation_Dim : public BoundaryIntegralEquation_Base
  {    
  public:
    Mesh<Dimension> mesh;
    MeshNumbering<Dimension> mesh_num;
    GeometryQuantityBEM<Dimension> var_geom;
    R3 wave_vector;
    
  public:
    BoundaryIntegralEquation_Dim();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void ConstructMesh(const Vector<string>& param, Mesh<Dimension>& mesh);
    void CheckOrientationNormales(Mesh<Dimension>& mesh, const R3& center);
    void ComputeMeshAndFiniteElement(const string& name_elt);
    
    virtual void ConstructFiniteElement(const string& name_elt) = 0;

    inline int GetNbBoundary() const { return mesh.GetNbBoundaryRef(); }
  };

  //! class for integral equations for a given dimension and given finite element
  template<class Dimension, int type>
  class BoundaryIntegralEquation_Fem : public BoundaryIntegralEquation_Dim<Dimension>
  {
    typedef typename Dimension::DimensionBoundary DimensionB;
    
  protected:
    ElementReference<DimensionB, type>* tri_element, *quad_element;
    Vector<const ElementReference<DimensionB, type>* > reference_element;
    
  public:
    BoundaryIntegralEquation_Fem();
    ~BoundaryIntegralEquation_Fem();
    
    void ConstructFiniteElement(const string& name_elt);
    void ComputeGeometryQuantities();

    inline const ElementReference<DimensionB, type>& GetReferenceElement(int i) const { return *reference_element(i); }
    
  };


  class BoundaryIntegralEquation_Maxwell3D : public BoundaryIntegralEquation_Fem<Dimension3, 3>
  {
  protected:
    Vector<Complex_wp> ref_epsilon, ref_mu, ref_sigma;
    R3 polarization;
    
  public:
    BoundaryIntegralEquation_Maxwell3D();
    
    virtual void InitIndices(int N);
    virtual void SetIndices(int, const Vector<string>&);
    virtual void SetPhysicalIndex(const string&, int ref, const Vector<string>&);

    void SetInputData(const string& description_field, const VectString& parameters);
    
    void ComputeElementaryMatrix(int, int, Vector<int>& row, Vector<int>& col,
				 Matrix<Complex_wp>& mat) const;

    void ComputeRightHandSide(Vector<Complex_wp>& b) const;
    void WriteSolution(const Vector<Complex_wp>& b) const;

  };
  
}

#define MONTJOIE_FILE_BOUNDARY_INTEGRAL_EQUATION_HXX
#endif
