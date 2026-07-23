#ifndef MONTJOIE_FILE_PHYSICAL_PROPERTY_HXX

#include "Algebra/TinySymmetricTensor.cxx"
#include "OneDimPhysicalIndex.hxx"

namespace Montjoie
{
  
  //! cut-off function on a parallepiped box
  /*!
    cut-off function considered :
    tronc(x) = exp(-coef/(x-xmin)^2)*exp(-coef/(x-xmax)^2)
    so that tronc(x) = 0 for x = xmin and x = xmax
    and tronc(x) = 1 when x is far from extremities
  */
  template<class Dimension, class T>
  class RectangleCutOff
  {
  protected :
    typedef TinyVector<T, Dimension::dim_N> R_N;
    TinyVector<R_N, 2> envelope; //!< bounding box
    R_N coef_gaussian; //!< parameters for the cut-off function
    int type_cut_off; //!< type of cut-off function to use
    //! available types of cut-off
    enum {NO_CUTOFF, GAUSSIAN_CUTOFF};
    
  public :
    RectangleCutOff();
    
    size_t GetMemorySize() const;
    void Init(const R_N& xmin_, const R_N & xmax_, const R_N& coef_);
    
    template<class T0>
    void ApplyCutOff(const R_N& point, T0& val) const;
    
    template<class T0>
    void ApplyCutOff(const R_N& point, T0& val,
                     TinyVector<T0, Dimension::dim_N>& grad) const;
    
  };
  
  
  //! sinusoidal field
  /*!
    Physical indexes (like epsilon, mu, ...) may vary sinusoidally
    inside a parallelepipedic domain
    The field will be equal to
    mu0 + amplitude*tronc(x)*tronc(y)*tronc(z)*sin(kx*x)*sin(ky*y)*sin(kz*z)
    tronc is a cut-off function, (kx,ky,kz) the wave vector
  */
  template<class Dimension, class T>
  class PhysicalSinusoidalMedia : public RectangleCutOff<Dimension, Real_wp>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    R_N kwave; //!< wave vector
    T amplitude; //!< coefficient
    T mu0; //!< coefficient
    bool inverse; //!< true if inverse has to be computed
    
  public :
    friend class PhysicalVaryingMedia<Dimension, T>; //!< friend class
    
    PhysicalSinusoidalMedia();
    
    // Inline methods
    bool IsComplex() const;
    void GetInverse();
    void Init(const R_N& xmin_, const R_N & xmax_, const R_N& coef_,
	      const T& offset, const T& ampli, const R_N& kwave_);
    
    void Mlt(const T& coef);
    
    // Other methods
    size_t GetMemorySize() const;
    void ComputeIndexAtPoint(const R_N& point, T& val,
                             TinyVector<T, Dimension::dim_N>& grad) const;
    
    void ComputeIndexAtPoint(const R_N& point, T& val) const;
    
    void ComputeIndex(const Mesh<Dimension>& mesh, const IVect&,
                      int ref, Vector<Vector<R_N> >& Points,
                      bool, bool, Vector<Vector<T> >& rho,
                      Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
		      Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& hess_rho);
    
  };
  
  
  //! field interpolated from values defined on a regular grid
  /*!
    Physical indexes (like epsilon, mu, ...) may be defined on a regular grid
    The field will be known by discrete values U_{i,j}
    and be equal to mu0 + amplitude*U_{i,j}
  */
  template<class Dimension, class T>
  class PhysicalRegularMedia : public RectangleCutOff<Dimension, Real_wp>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    TinyVector<int, Dimension::dim_N> nbx; //!< number of points of the regular grid
    T amplitude; //!< coefficient
    T mu0; //!< coefficient
    Vector<T> nu; //!< discrete values on regular grid
    bool inverse; //!< true if inverse has to be computed
    bool complex_index; //!< complex index ?
    int periodic; //!< quasi-periodic media ?
    VectReal_wp xstep, ystep, zstep; //!< not used
    R_N step_delta_x; //!< space step along each coordinate
    
  public :
    enum { USUAL, PERIODICITY, QUASI_PERIODICITY};
    
    friend class PhysicalVaryingMedia<Dimension, T>; //!< friend class
    
    PhysicalRegularMedia();    
    
    // Inline methods
    void SetNoPeriodic();
    void SetPeriodic();
    void SetQuasiPeriodic();
    
    bool IsComplex() const;
    void GetInverse();
    void Mlt(const T& coef);

  protected:        
    Real_wp GetIndexBounds(int i, int& imin,
			   int& imax, int n, const Real_wp& dx) const;
    
  public:
    // Other methods
    size_t GetMemorySize() const;
    
    void Init(const R_N& xmin_, const R_N & xmax_, const R_N& coef_, const T& offset,
              const T& ampli, int dim_N, const string& data_file,
	      bool complex_number, bool double_prec);
    
    void Clear();
    
    void ReadFile(const string& file_name, int dim_N, bool complex_number,
		  bool double_prec, const Real_wp& cte);

    void ReadFile(const string& file_name, int dim_N, bool complex_number,
		  bool double_prec, const Complex_wp& cte);

    void ComputeIndexAtPoint(const R_N& point, T& coef) const;
    
    void ComputeIndexAtPoint(const R_N& point, T& coef, 
                             TinyVector<T, Dimension::dim_N>& grad_coef) const;
    
    void ComputeIndex(const Mesh<Dimension>& mesh, const IVect&, int ref,
                      Vector<Vector<R_N> >& Points, bool, bool, Vector<Vector<T> >& rho,
                      Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
		      Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& hess_rho);
    
  };
  
  
  //! field interpolated from values defined on a mesh
  /*!
    Physical indexes (like epsilon, mu, ...) may be defined on a regular grid
    The field will be known by discrete values U_i
    and be equal to mu0 + amplitude*U_i
  */
  template<class Dimension, class T>
  class PhysicalMeshMedia
  {
  protected :    
    //! components of the indice on the nodal points of the mesh
    Vector<Vector<T> > nu;
    //! gradient of nu on nodal points of the mesh
    Vector<Vector<TinyVector<T, Dimension::dim_N> > > grad_nu;
    //! hessian of nu on nodal points of the mesh
    Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > > hessian_nu;
    //! index has to be inverted ?
    bool inverse;
    bool complex_index; //!< complex index ?
    //! file name of the mesh
    string mesh_file_name;
    //! order of the mesh
    int mesh_order;
    
  public :    
    typedef typename Dimension::R_N R_N;
    friend class PhysicalVaryingMedia<Dimension, T>; //!< friend class
    friend class GenericPhysicalIndice<Dimension, T>; //!< friend class
    
    PhysicalMeshMedia();    
    
    // Inline methods
    bool IsComplex() const;
    void GetInverse();
    void Mlt(const T& coef);
    
    const string& GetMeshName() const;
    int GetOrder() const;
    
    int GetNbElt() const;
    int GetNbNodes(int i) const;
    
    const T& GetValue(int i, int j) const;
    const TinyVector<T, Dimension::dim_N>& GetGradient(int i, int j) const;
    TinyVector<T, Dimension::dim_N> GetHessian(int i, int j, int r) const;
    
    // Other methods
    size_t GetMemorySize() const;
    void Clear();

    void ReadFile(const string & name_mesh, const string& file_name, int order);
    
    void Init(const string& file_mesh, const T& offset, const T& amplitude,
	      int dim_N, const string& data_file, int order);
    
    void RemoveGlobalElements(const IVect& NumLoc);
     
    void CheckMesh(const Mesh<Dimension>& mesh);
    void ComputeGradient(const Mesh<Dimension>&, bool hessian = false);

    void ComputeIndexAtPoint(int i, const VectReal_wp& phi, T& coef) const;
    
  };
  
  
  //! radial field
  template<class Dimension, class T>
  class RadialVaryingMedia
  {
  protected :
    //! 1-D index
    VariableParameter1D_Base<T>* index_radial;
    //! index has to be inverted ?
    bool inverse;
    //! complex index ?
    bool complex;
    
  public :
    typedef typename Dimension::R_N R_N;
    
    RadialVaryingMedia();
    ~RadialVaryingMedia();

    RadialVaryingMedia(const RadialVaryingMedia<Dimension, T>&);
    RadialVaryingMedia<Dimension, T>& operator=(const RadialVaryingMedia<Dimension, T>&);
    
    // Other methods
    size_t GetMemorySize() const;
    void InitSpline(const string& data_file, const T& offset,
                    const T& amplitude, T& cte_infty);
    
    bool IsComplex() const;
    void GetInverse();
    void Mlt(const T& coef);
    void Clear();
    
    void ComputeIndexAtPoint(const R_N& point, T& coef) const;
    
    void ComputeIndexAtPoint(const R_N& point, T& coef, 
                             TinyVector<T, Dimension::dim_N>& grad_coef,
                             TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>& hess_coef,
                             bool hessian) const;
    
    void ComputeIndex(const Mesh<Dimension>& mesh, const IVect&, int ref,
                      Vector<Vector<R_N> >& Points, bool, bool, Vector<Vector<T> >& rho,
                      Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
		      Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& hess_rho);
    
  };
  
  
  //! varying field 
  /*!
    The field can be constant, sinusoidal, defined on a grid or a mesh
    This class is used for varying physical indexes (like epsilon, mu)
  */
  template<class Dimension, class T>
  class PhysicalVaryingMedia
  {
    friend class GenericPhysicalIndice<Dimension, T>; //!< friend class
    
  protected :
    //! type of variation of the field
    int type;
    //! coefficient on quadrature points
    Vector<Vector<T> > eval_coef;
    //! gradient of coefficient (if needed)
    Vector<Vector<TinyVector<T, Dimension::dim_N> > > grad_coef;
    //! hessian of coefficient (if needed)
    Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > > hessian_coef;
    T cte_coef; //!< constant coefficient if constant field
    T offset_coef, amplitude_coef; //!< for user-defined fields
    TinyVector<T, Dimension::dim_N> zero_grad; //!< null vector
    TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> zero_hess;
    bool same_mesh; //!< true if the mesh is the same as the one used to discretize the unknown
        
  public :   
    typedef typename Dimension::R_N R_N;
    
    PhysicalRegularMedia<Dimension, T> index_regular; //!< field defined on a regular grid
    PhysicalSinusoidalMedia<Dimension, T> index_sinus; //!< sinusoidal field
    PhysicalMeshMedia<Dimension, T> index_mesh; //!< field defined on a mesh
    RadialVaryingMedia<Dimension, T> index_radial; //!< field versus the radius r
    
    // available types of variation
    enum {CONSTANT, SINUSOIDE, REGULAR_GRID, MESH, RADIAL, USER};
    
    PhysicalVaryingMedia();    
    PhysicalVaryingMedia(const T& cte);
    PhysicalVaryingMedia(int cte);
    
    // Inline methods
    bool IsVarying() const;
    bool IsZero() const;
    bool IsSameMesh() const;
    
    void SetMediaType(int);
    int GetMediaType() const;
        
    T GetConstant() const;
    const T& GetOffset() const;
    const T& GetAmplitude() const;

    void InitConstant(const T& constante);
    
    void InitSinus(const R_N& xmin, const R_N& xmax, const R_N& coef,
		   const T& offset, const T& amplitude, const R_N& kwave);
   
    void InitRandom(const R_N& xmin, const R_N& xmax, const R_N& coef,
		    const T& offset, const T& amplitude, int dim_N,
		    const string& data_file, bool, bool);
    
    void InitMesh(const string& mesh_file, const T& offset, const T& amplitude,
		  int dim_N, const string& data_file, int order);
    
    void InitRadial(const string& data_file, const T& offset,
                    const T& amplitude, T& cte);
    
    void InitUser(const T& offset, const T& amplitude, const T& cte);
    
    void SetSameMesh();
    
    Vector<Vector<T> >& GetValue();
    Vector<Vector<TinyVector<T, Dimension::dim_N> > >& GetGradient();
    Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& GetHessian();
    
    void SetCoefficient(int i, int j, const T& coef);
    void SetCoefGradient(int i, int j, const TinyVector<T, Dimension::dim_N>& grad);
    void SetCoefHessian(int i, int j, int r, const TinyVector<T, Dimension::dim_N>& row);

    const T& GetCoefficient(int i, int j) const;
    const TinyVector<T, Dimension::dim_N>& GetCoefGradient(int i, int j) const;
    const TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>&
    GetCoefHessian(int i, int j) const;
    
    // Other methods
    size_t GetMemorySize() const;
    void Reallocate(int i, bool, bool compute_hess = false);
    void ReallocateOnElement(int i, int j, bool, bool compute_hess = false);
    
    bool IsComplex() const;
    
    void GetInverse();
    void Mlt(const T& coef);
    void Zero();

    void Clear();
    void ClearHessianGradient();
        
    T GetIndexAtPoint(int i, const VectReal_wp& phi, const R_N& pt_glob) const;
    
  };
  
  
  //! base class to read a varying field in the data file
  template<class Dimension, class T>
  class GenericPhysicalIndice
  {
  public :    
    typedef typename Dimension::R_N R_N;
    
    static void SetInputVaryingMedia(int& nb,
                                     PhysicalVaryingMedia<Dimension, T>& index,
                                     T& cte, const VectString& parameters);
  };
  

  //! class for a varying scalar field 
  /*!
    The field comprises only one component,
    it can be constant or varying
  */
  template<class Dimension, class T>
  class ScalarPhysicalIndice : public GenericPhysicalIndice<Dimension, T>
  {
  protected :
    int type_media; //!< constant or varying field ?
    enum {CONSTANT, VARYING};
    PhysicalVaryingMedia<Dimension, T> fct_rho; //!< varying field
    T cte_rho; //!< constant value
    
  public : 
    static Real_wp threshold; //!< threshold
    
    ScalarPhysicalIndice();
    
    // Inline methods
    operator T() const;      
    
    bool IsZero() const;
    bool IsComplex() const;
    bool IsVarying() const;

    void SetMediaType(int);
    int GetMediaType() const;
    inline PhysicalVaryingMedia<Dimension, T>& GetVaryingIndex() {  return fct_rho; }

    void SetIdentity();
    const T& GetConstant() const;
    void SetConstant(const T& coef);
    void Zero();
    void GetInverse();
    void Mlt(const T& coef);

    T GetCoefficient(const VarPhysicalProblem& var, int i, int j) const;
    
    void GetCoefGradient(const VarPhysicalProblem& var, int i, int j,
                         T&, TinyVector<T, Dimension::dim_N>& ) const;
    
    void GetCoefHessian(const VarPhysicalProblem& var, int i, int j,
			T&, TinyVector<T, Dimension::dim_N>&,
			TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>&);

    void SetInputData(int& nb, const VectString& parameters,
		      const string& keyword);
        
    // Other methods
    size_t GetMemorySize() const;
    int GetNbVaryingMedia() const;
    void GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho, IVect &);
    
    void Clear();			
    
    //! friend method
    template<class Dim, class T0>
    friend ostream& operator <<(ostream& out, const ScalarPhysicalIndice<Dim, T0>& var);
    
  };
  
  template<class Dimension, class T>
  Real_wp ScalarPhysicalIndice<Dimension, T>::threshold(1e-30);
  
  
  //! class for a varying vectorial field 
  /*!
    The field comprises only two or three components, depending on the dimension
    it can be constant or varying
  */
  template<class Dimension, int m, class T>
  class VectorPhysicalIndice : public GenericPhysicalIndice<Dimension, T>
  {
  protected :
    int type_media; //!< varying or constant field ?
    enum {CONSTANT, VARYING};
    TinyVector<PhysicalVaryingMedia<Dimension, T>, m> fct_rho; //!< varying field
    TinyVector<T, m> cte_rho; //!< constant vector
    static Real_wp threshold; //!< threshold
    
  public : 
    VectorPhysicalIndice();
    
    // Inline methods
    const T& operator()(int i) const;
    const TinyVector<T, m>& GetConstant() const;
      
    int GetMediaType() const;
    int GetMediaType(int) const;
    
    void SetIdentity();
    void SetConstant(const TinyVector<T, m>& coef);
    void Zero();
    void Mlt(const T& coef);

    inline PhysicalVaryingMedia<Dimension, T>& GetVaryingIndex(int p) {  return fct_rho(p); }
    TinyVector<T, m> GetCoefficient(const VarPhysicalProblem& var, int i, int j) const;
    
    // Other methods
    size_t GetMemorySize() const;
    bool IsComplex() const;
    bool IsVarying() const;
    bool IsZero() const;
    
    void SetMediaType(int);
    void SetMediaType(int, int);

    int GetNbVaryingMedia() const;
    void GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho, IVect &);

    void Clear();

    void GetCoefGradient(const VarPhysicalProblem& var, int i, int j, TinyVector<T, m>&,
			 TinyMatrix<T, General, m, Dimension::dim_N>& ) const;
    
    void GetCoefHessian(const VarPhysicalProblem& var, int i, int j,
			TinyVector<T, m>&, TinyMatrix<T, General, m, Dimension::dim_N>&,
			TinyVector<TinyMatrix<T, Symmetric,
			Dimension::dim_N, Dimension::dim_N>, m>&);
    
    void SetInputData(int& nb,
		      const VectString& parameters, const string& keyword);
    
    //! friend method
    template<class Dimension0, int m0, class T0>
    friend ostream& operator <<(ostream& out,
				const VectorPhysicalIndice<Dimension0, m0, T0>& var);    
  };
  
  template<class Dimension, int m, class T>
  Real_wp VectorPhysicalIndice<Dimension, m, T>::threshold(1e-30);
  
  
  //! class for a varying tensorial field 
  /*!
    The field comprises three or six components
    d x d symmetric matrix, where d is the dimension 
    it can be constant or varying
  */
  template<class Dimension, int m, class T>
  class TensorPhysicalIndice  : public GenericPhysicalIndice<Dimension, T>
  {
  protected :
    int type_media; //!< constant or varying field ?
    int type_anisotropy; //!< isotrope or anisotrope field ?

    //! varying field
    TinyMatrix<PhysicalVaryingMedia<Dimension, T>, Symmetric, m, m> fct_mu;
    //! constant tensor
    TinyMatrix<T, Symmetric, m, m> cte_mu;
    static Real_wp threshold; //!< tensor
    
  public : 
    enum {CONSTANT, VARYING};
    enum {ISOTROPE, ORTHOTROPE, ANISOTROPE};

    TensorPhysicalIndice();
    
    // Inline methods
    const TinyMatrix<T, Symmetric, m, m>& GetConstant() const;
    
    const T& operator()(int i, int j) const;
    
    int GetAnisotropy() const;
    void SetAnisotropy(int type);
    
    int GetMediaType() const;
    int GetMediaType(int, int) const;

    void SetIdentity();
    int GetSize() const;

    inline PhysicalVaryingMedia<Dimension, T>& GetVaryingIndex(int p, int q) {  return fct_mu(p, q); }
    
    // Other methods
    size_t GetMemorySize() const;
    void SetMediaType(int);
    void SetMediaType(int, int, int);

    template<int p, class T0, class T1>
    void MltMatrix(const VarPhysicalProblem& var, int i, int j,
                   const TinyVector<T0, p>& u, TinyVector<T1, p>& v) const;

    template<class T0, class T1>
    void MltMatrix(const VarPhysicalProblem& var, int i, int j,
                   const TinyVector<T0, 2>& u, TinyVector<T1, 2>& v) const;

    template<class T0, class T1>
    void MltMatrix(const VarPhysicalProblem& var, int i, int j,
                   const TinyVector<T0, 3>& u, TinyVector<T1, 3>& v) const;

    int GetNbVaryingMedia() const;
    void GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho,
                         IVect &);

    void SetDiagonal(const T& coef);
    void SetConstant(const TinyMatrix<T, Symmetric, m, m>& coef);
    
    void Clear();
    void Zero();
    void GetInverse();
    void Mlt(const T& coef);
    
    bool IsComplex() const;
    bool IsVarying() const;
    bool IsZero() const;    
    
    void MltMatrixPoint(int i1, const VectReal_wp& phi,
			const typename Dimension::R_N& pt_glob,
			const TinyVector<T, m>& u, TinyVector<T, m>& v) const;

    TinyMatrix<T, Symmetric, m, m>
    GetMatrixPoint(int i1, const VectReal_wp& phi, 
		   const typename Dimension::R_N& pt_glob) const;
    
    void SetInputData(int& nb,
		      const VectString& parameters, const string& keyword);
    
    
    TinyMatrix<T, Symmetric, m, m>
    GetCoefficient(const VarPhysicalProblem& var, int i, int j) const;
    
    T GetCoefficient(const VarPhysicalProblem& var, int i, int j, int p, int q) const;
    
    void GetCoefGradient(const VarPhysicalProblem& var, int i, int j,
                         TinyMatrix<T, Symmetric, m, m>&,
			 TinyArray3D<T, m, m, Dimension::dim_N>& ) const;
    
    //! friend method
    template<class Dimension0, int m0, class T0>
    friend ostream& operator <<(ostream& out,
				const TensorPhysicalIndice<Dimension0, m0, T0>& var);

  };
  
  template<class Dimension, int m, class T>
  Real_wp TensorPhysicalIndice<Dimension, m, T>::threshold(1e-30);
  
  
  //! class for "elastic-tensor" field
  /*!
    Tensor C for elastodynamic problems
  */
  template<class Dimension, int m, class T>
  class ElasticPhysicalIndice : public GenericPhysicalIndice<Dimension, T>
  {
  protected :
    enum {nb_components = m*(m+1)/2};
    int type_media; //!< constant or varying field ?
    int type_anisotropy; //!< isotrope or anisotrope field ?

    //! full tensor C in anisotropic case
    TinySymmetricTensor<T, m> cte_C;
    //! for varying media
    TinySymmetricTensor<PhysicalVaryingMedia<Dimension, T>, m> fct_C;
    TinyVector<int, m> indic;

    template<class Vector1>
    void MltMatrixVar(int i1, int j, const TinyVector<int, 2>&,
                      const TinyVector<Vector1, 2>& dU, TinyVector<Vector1, 2>& dV) const;

    template<class Vector1>
    void MltMatrixVar(int i1, int j, const TinyVector<int, 3>&,
                      const TinyVector<Vector1, 3>& dU, TinyVector<Vector1, 3>& dV) const;

    template<class Vector1>
    void MltVectorVar(int i1, int j, const TinyVector<int, 2>&,
                      const Vector1& dU, Vector1& dV) const;

    template<class Vector1>
    void MltVectorVar(int i1, int j, const TinyVector<int, 3>&,
                      const Vector1& dU, Vector1& dV) const;
    
  public :
    typedef typename Dimension::R_N R_N;
    enum {CONSTANT, VARYING};
    enum {ISOTROPE, ORTHOTROPE, ANISOTROPE};

    ElasticPhysicalIndice();

    // Inline methods
    void SetMediaType(int);
    int GetMediaType() const;
    int GetAnisotropyType() const;
    
    const T& operator()(int i, int j, int k, int l) const;

    void Zero();    
    void SetIdentity();
    const TinySymmetricTensor<T, m>& GetConstant() const;

    template<class Vector1>
    void MltMatrix(const VarPhysicalProblem& var, int i, int j,
		   const TinyVector<Vector1, m>& dU, TinyVector<Vector1, m>& dV) const;

    template<class Vector1>
    void MltVector(const VarPhysicalProblem& var, int i, int j,
		   const Vector1& dU, Vector1& dV) const;

    bool IsVarying() const;

    inline PhysicalVaryingMedia<Dimension, T>& GetVaryingIndex(int i, int j, int k, int l)
    {  return fct_C(i, j, k, l); }
    
    // Other methods
    size_t GetMemorySize() const;
    static void InitStaticData();
    
    int GetNbVaryingMedia() const;
    void GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho, IVect &);

    template<class T0>
    void GetCoefficient(const VarPhysicalProblem& var, int i, int j,
			TinyMatrix<TinyMatrix<T0, General, m, m>, General, m, m>& ) const;
    
    void Mlt(const T& coef);
    void SetConstant(const TinySymmetricTensor<T, m>& coef);
        
    void GetInverse();

    template<class Vector1>
    void MltVectorPoint(int i1, const VectReal_wp& phi,
                        const R2& pt_glob, const Vector1& dU, Vector1& dV) const;

    template<class Vector1>
    void MltVectorPoint(int i1, const VectReal_wp& phi,
                        const R3& pt_glob, const Vector1& dU, Vector1& dV) const;
    
    bool IsComplex() const;
    void Clear();
    
    void SetInputData(int& nb, const VectString& parameters,
		      const string& keyword);
    
    //! friend method
    template<class Dimension0, int m0, class T0>
    friend ostream& operator <<(ostream& out,
				const ElasticPhysicalIndice<Dimension0, m0, T0>& var);
    
  };
  
}

#define MONTJOIE_FILE_PHYSICAL_PROPERTY_HXX
#endif

