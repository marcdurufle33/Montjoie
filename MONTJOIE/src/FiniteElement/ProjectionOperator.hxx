#ifndef MONTJOIE_FILE_PROJECTION_OPERATOR_HXX

namespace Montjoie
{
  
  class FiniteElementProjector
  {
  protected :
    enum {SPARSE, IDENTITY};
    int type_projector;
    IVect num_row; int nb_row;
    int nb_nodes_u, nb_nodes_v;

  public :
    FiniteElementProjector();
    virtual ~FiniteElementProjector();
    
    void SetIdentity(const IVect& icol, const IVect& irow, int n);
    
    template<class T1, class T2, int nb_unknowns>
    void Project(const Vector<TinyVector<T1, nb_unknowns> >& Un,
                 Vector<TinyVector<T2, nb_unknowns> >& Vn) const;
    
    template<class T1, class T2, int p, int q>
    void Project(const Vector<TinyMatrix<T1, General, p, q> >& Un,
                 Vector<TinyMatrix<T2, General, p, q> >& Vn) const;
    
    template<class T1, class T2, int nb_unknowns>
    void TransposeProject(const Vector<TinyVector<T1,nb_unknowns> >& Un,
                          Vector<TinyVector<T2, nb_unknowns> >& Vn) const;
    
    template<class T1, class T2>
    void Project(const Vector<T1>& Un, Vector<T2>& Vn) const;

    template<class T1, class T2>
    void TransposeProject(const Vector<T1>& Un, Vector<T2>& Vn) const;
    
    virtual size_t GetMemorySize() const = 0;
    
    virtual void Init(const ElementGeomReference<Dimension1>& FaceBasis,
		      const VectReal_wp& points_div) = 0;
    
    virtual void Init(const ElementGeomReference<Dimension1>&, const VectReal_wp&,
		      const VectReal_wp& points_div) = 0;

    virtual void Init(const ElementGeomReference<Dimension2>& FaceBasis,
		      const VectR2& points_div) = 0;
    
    virtual void Init(const ElementGeomReference<Dimension2>&, const VectReal_wp&,
		      const VectR2& points_div) = 0;

    virtual void Init(const ElementGeomReference<Dimension3>& FaceBasis,
		      const VectR3& points_div) = 0;
    
    virtual void Init(const ElementGeomReference<Dimension3>&, const VectReal_wp&,
		      const VectR3& points_div) = 0;

    virtual void ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const = 0;
    virtual void ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const = 0;

    virtual void TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const = 0;
    virtual void TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const = 0;
    
  };
  
  
  //! base class to interpolate on tensorized domains (quad/hexs)
  /*!
    It is used to compute \f$ u(\zeta_{i_1}, \zeta_{i_2}) \, = \, 
    \sum_{j_1,j_2} \varphi_{j_1}(\zeta_{i_1}) \varphi_{j_2}(\zeta_{i_2}) u_{j_1,j_2} \f$
    double sums (in 2-D) and triple sums (in 3-D) can be decomposed in two simple sums,
    and three simple sums.
    A simple sum is equivalent to a matrix vector product, with sparse matrix.
   */
  template<class Dimension>
  class TensorizedProjector_Base : public FiniteElementProjector
  {
  protected :
    //! interpolation operator on each direction
    TinyVector<Matrix<Real_wp, General, RowSparse>, Dimension::dim_N> ch_loc;
    
  public :
    TensorizedProjector_Base();
    
    size_t GetMemorySize() const;
    
    template<class MatrixSparse>
    void SetChOperator(int m, const MatrixSparse& A);
    
  };
  
  
  //! base class to interpolate on any domains (tri/tet)
  /*!
    It is used to compute \f$ u(\zeta_{i}) \, = \, \sum_{j} \varphi_{j}(\zeta_{i}) u_j \f$
    There is no efficient structures in points, like in quads.
    This interpolation is a matrix vector product with a full matrix
   */
  template<class Dimension>
  class DenseProjector : public FiniteElementProjector
  {
  protected :
    //! operator projection is a full matrix
    Matrix<Real_wp> ch_loc;
    
  public :
    DenseProjector();
    
    size_t GetMemorySize() const;

  protected:
    void InitProjector(const ElementGeomReference<Dimension>& FaceBasis,
                       const typename Dimension::VectR_N& points_div);
    
    void InitProjector(const ElementGeomReference<Dimension>& FaceBasis, const VectReal_wp& step_subdiv,
                       const typename Dimension::VectR_N& points_div);

    template<class Dim>
    void InitProjector(const ElementGeomReference<Dim>&, const typename Dim::VectR_N& Pts);

    template<class Dim>
    void InitProjector(const ElementGeomReference<Dim>&, const VectReal_wp& subdiv,
                       const typename Dim::VectR_N& Pts);

  public:
    void Init(const ElementGeomReference<Dimension1>& FaceBasis,
              const VectReal_wp& points_div);
    
    void Init(const ElementGeomReference<Dimension1>&, const VectReal_wp&,
              const VectReal_wp& points_div);
    
    void Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectR2& points_div);
    
    void Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectReal_wp& step_subdiv,
              const VectR2& points_div);
    
    void Init(const ElementGeomReference<Dimension3>& FaceBasis, const VectR3& points_div);
    
    void Init(const ElementGeomReference<Dimension3>& FaceBasis,
              const VectReal_wp& step_subdiv, const VectR3& points_div);

    void ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const;
    void ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const;

    void TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const;
    void TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const;
    
  };
  

  //! class to interpolate on quadrilaterals
  /*!
    It is used to compute \f$ u(\zeta_{i_1}, \zeta_{i_2}) \, = \,
    \sum_{j_1,j_2} \varphi_{j_1}(\zeta_{i_1}) \varphi_{j_2}(\zeta_{i_2}) u_{j_1,j_2} \f$
    This double sum is decomposed in :
    \f$ v_{i_1,j_2} = \sum_{j_1} \, \varphi_{j_1}(\zeta_{i_1}) u_{j_1,j_2} \f$
    \f$ w_{i_1,i_2} = \sum_{j_2} \, \varphi_{j_2}(\zeta_{i_2}) v_{i_1,j_2} \f$
    Theses operations are matrix-vector products :
    \f$ V \, = \, C_1 \, U \quad W \, = \, C_2 \, V \f$
    the matrices C1 and C2 are stored in ch_loc
   */
  template<>
  class TensorizedProjector<Dimension2> : public TensorizedProjector_Base<Dimension2>
  {
  protected :
    //! tensorial <-> scalar numerotation for target points
    Matrix<int> NumNodes;
    Vector<Real_wp> Ux_real;
    Vector<Complex_wp> Ux_cplx;
    
  public :
    TensorizedProjector();

    size_t GetMemorySize() const;

    void Init(const ElementGeomReference<Dimension1>& FaceBasis,
              const VectReal_wp& points_div);
    
    void Init(const ElementGeomReference<Dimension1>&, const VectReal_wp&,
              const VectReal_wp& points_div);
    
    void Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectR2& points_div);
    
    void Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectReal_wp& step_subdiv,
              const VectR2& points_div);

    void Init(const ElementGeomReference<Dimension3>& FaceBasis,
		      const VectR3& points_div);
    
    void Init(const ElementGeomReference<Dimension3>&, const VectReal_wp&,
              const VectR3& points_div);

  protected:
    void InitProjector(const ElementGeomReference<Dimension2>& FaceBasis, const VectReal_wp& step_x,
                       const VectReal_wp& step_y, const VectR2& points_div);
    
    template<class TypeFunction>
    void InitProjector(const TypeFunction& lob, const Matrix<int>& NumNodes2D,
                       const VectReal_wp& step_subdiv, const VectR2& points_div);
    
    template<class TypeFunction>
    void InitProjector(const TypeFunction& lob, const Matrix<int>& NumNodes2D,
                       const VectReal_wp& step_x, const VectReal_wp& step_y, const VectR2& points_div);

  public:
    void ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const;
    void ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const;

    void TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const;
    void TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const;
    
  };
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! class to interpolate on hexahedrals
  /*!
    It is used to compute \f$ u(\zeta_{i_1}, \zeta_{i_2}, \zeta_{i,3}) \, = \,
    \sum_{j_1,j_2,j_3} \varphi_{j_1}(\zeta_{i_1}) \varphi_{j_2}(\zeta_{i_2})
    \varphi_{j_3}(\zeta_{i_3}) u_{j_1,j_2,j_3} \f$
    This triple sum is decomposed in :
    \f$ v_{i_1,j_2,j_3} = \sum_{j_1} \, \varphi_{j_1}(\zeta_{i_1}) u_{j_1,j_2,j_3} \f$
    \f$ w_{i_1,i_2,j_3} = \sum_{j_2} \, \varphi_{j_2}(\zeta_{i_2}) v_{i_1,j_2,j_3} \f$
    \f$ z_{i_1,i_2,i_3} = \sum_{j_3} \, \varphi_{j_3}(\zeta_{i_3}) w_{i_1,i_2,j_3} \f$
    Theses operations are matrix-vector products :
    \f$ V \, = \, C_1 \, U \quad W \, = \, C_2 \, V \quad Z \, = \, C_3 \, W \f$
    the matrices C1, C2 and C3 are stored in ch_loc
   */
  template<>
  class TensorizedProjector<Dimension3> : public TensorizedProjector_Base<Dimension3>
  {
  protected :
    //! tensorial <-> scalar numerotation for target points
    Array3D<int> NumNodes;
    Vector<Real_wp> Ux_real, Uy_real;
    Vector<Complex_wp> Ux_cplx, Uy_cplx;
    
  public :

    TensorizedProjector();

    size_t GetMemorySize() const;

    void Init(const ElementGeomReference<Dimension1>& FaceBasis,
              const VectReal_wp& points_div);
    
    void Init(const ElementGeomReference<Dimension1>&, const VectReal_wp&,
              const VectReal_wp& points_div);
    
    void Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectR2& points_div);
    
    void Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectReal_wp& step_subdiv,
              const VectR2& points_div);
    
    void Init(const ElementGeomReference<Dimension3>& FaceBasis, const VectR3& points_div);
    
    void Init(const ElementGeomReference<Dimension3>& FaceBasis,
              const VectReal_wp& step_subdiv, const VectR3& points_div);

  protected:
    void InitProjector(const ElementGeomReference<Dimension3>& FaceBasis, const VectReal_wp& step_x,
                       const VectReal_wp& step_y, const VectReal_wp& step_z, const VectR3& points_div);

    template<class TypeFunction>
    void InitProjector(const TypeFunction& lob, const Array3D<int>& NumNodes3D,
                       const VectReal_wp& step_subdiv, const VectR3& points_div);

    template<class TypeFunction>
    void InitProjector(const TypeFunction& lob, const Array3D<int>& NumNodes3D, const VectReal_wp& step_x,
                       const VectReal_wp& step_y, const VectReal_wp& step_z, const VectR3& points_div);

  public:
    void AllocateTemporaryVectors();
    
    void ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const;
    void ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const;

    void TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const;
    void TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const;

  };
#endif
  
  //! class for interpolation on a set of finite elements
  /*!
    For a given order of approximation, this interpolator
    admits different types of elements (tet, hexa, prism, pyramid in 3-D)
   */
  class FiniteElementInterpolator
  {
  protected:
    Vector<FiniteElementProjector*> stored_proj;
    
  public :
    ~FiniteElementInterpolator();
    
    size_t GetMemorySize() const;
    
    template<class FeElement, class Vector1>
    void InitProjection(const Vector<FeElement*>&, const Vector1& points_div);

    template<class FeElement, class Vector1, class Vector2>
    void InitProjection(const Vector<FeElement*>&,
			const Vector1& step_subdiv, const Vector2& points_div);
    
    template<class T1, class T2, int nb_unknowns>
    void Project(const Vector<TinyVector<T1, nb_unknowns> >& U,
                 Vector<TinyVector<T2, nb_unknowns> >& V, int type_elt) const;
    
    template<class T1, class T2, int nb_unknowns>
    void TransposeProject(const Vector<TinyVector<T1, nb_unknowns> >& U,
                          Vector<TinyVector<T2, nb_unknowns> >& V, int type_elt) const;

    template<class T1, class T2, int p, int q>
    void Project(const Vector<TinyMatrix<T1, General, p, q> >& U,
                 Vector<TinyMatrix<T2, General, p, q> >& V, int type_elt) const;

    template<class Vector1, class Vector2>
    void ProjectScalar(const Vector1& U, Vector2& V, int type_elt) const;
    
    template<class Vector1, class Vector2>
    void TransposeProjectScalar(const Vector1& U, Vector2& V, int type_elt) const;

    template<class T1, class T2>
    void Project(const Vector<T1>& U, Vector<T2>& V, int type_elt) const;
    
    template<class T1, class T2>
    void TransposeProject(const Vector<T1>& U, Vector<T2>& V, int type_elt) const;
    
    template<class Dimension>
    void ComputeLocalProlongationLowOrder(const Vector<const ElementReference_Dim<Dimension>* >& elt_fine,
					  int rf,
					  const Vector<const ElementReference_Dim<Dimension>* >& elt_coarse,
					  int rc,
					  TinyVector<Matrix<Real_wp>, 4>& ProlongationMatrix);

    template<class Dimension>
    void ComputeLocalProlongation(const Vector<const ElementReference_Dim<Dimension>* >& elt_fine,
				  int rf,
				  const Vector<const ElementReference_Dim<Dimension>* >& elt_coarse,
				  int rc,
				  TinyVector<Matrix<Real_wp>, 4>& ProlongationMatrix);
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_PROJECTION_OPERATOR_HXX
#endif

  
