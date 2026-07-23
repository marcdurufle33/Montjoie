#ifndef MONTJOIE_FILE_DEFINE_SOURCE_ELLIPTIC_HXX

namespace Montjoie
{
  
  //! base class for source computation
  class VarSourceProblem_Base
  {
  public :
    // different types of sources
    enum {SRC_NULL, SRC_VOLUME, SRC_TOTAL_FIELD, SRC_DIFFRACTED_FIELD,
	  SRC_USER, SRC_DIRAC};
    
    // incident fields
    enum {INCIDENT_PLANE_WAVE, INCIDENT_PLANE_WAVE_CPLX, INCIDENT_GAUSSIAN_BEAM,
          INCIDENT_HANKEL, INCIDENT_LAYERED_PLANE_WAVE, INCIDENT_NONE};

  private :
    VarComputationProblem& var_computation;
    VarBoundaryCondition_Base& var_boundary;
    VarProblem_Base& var_problem;
    
    void InitDefaultValues();
    
  protected :
    //! type of the source (diffracted wave, point source)
    Vector<int> type_source;
    Vector<int> incident_source; //!< type of incident wave (e.g. plane waves, gaussian beams)

    //! threshold used to neglect higher modes
    Real_wp threshold_rhs;
    //! list of parameters for the spatial sources
    Vector<Vector<VectString> > source_space_param;
    
  public :
    template<class TypeEquation>
    VarSourceProblem_Base(EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    bool OnlyOneSource() const;
    int GetNbRhs() const;
    int GetSourceType(int n) const;
    Real_wp GetThresholdSource() const;
    int GetIncidentFieldType(int n) const;
    int GetNbParameterSource(int n) const;
    const VectString& GetParameterSource(int n, int num) const;
    
    void GetMemoryUsed(map<string, size_t>& var) const;

    // Other methods    
    void SetInputData(const string& description_field,
 		      const VectString& parameters);    

    void SetParameterSource(const Vector<VectString>& param);
    static int GetIncidentFieldType(const string& name);
    
    void UpdateSourceParam(const VectString& parameters, Vector<int>& type, Vector<int>& incident);

    template<class Dimension>
    void InitGaussianParameter(GaussianSource<Dimension>& f, const VectString& param, int& nb) const;

    template<class Dimension>
    void InitRandomGaussianParameter(GaussianSource<Dimension>& f, const VectString& param, int nb) const;
    
    void InitFftComputation(FftInterface<Complex_wp>& fft_interface) const;
    
    void CopyInputData(const VarSourceProblem_Base& var);

    // almost pure virtual methods
    virtual void ModifyVolumetricSource(int i, int j, const R2&,
					const VirtualSourceField<Real_wp, Dimension2>&,
					Vector<Real_wp>&) const;

    virtual void ModifyVolumetricSource(int i, int j, const R2&,
					const VirtualSourceField<Complex_wp, Dimension2>&,
					Vector<Complex_wp>&) const;

    virtual void ModifyVolumetricSource(int i, int j, const R3&,
					const VirtualSourceField<Real_wp, Dimension3>&,
					Vector<Real_wp>&) const;

    virtual void ModifyVolumetricSource(int i, int j, const R3&,
					const VirtualSourceField<Complex_wp, Dimension3>&,
					Vector<Complex_wp>&) const;
    
    virtual void AddIncidentWave(const Real_wp&, Vector<VectReal_wp>&) = 0;
    virtual void AddIncidentWave(const Complex_wp&, Vector<VectComplex_wp>&) = 0;

    virtual void InitIncidentField() = 0;
    
    virtual void ComputeRightHandSide(VectReal_wp& b_source, bool assemble = false) = 0;    
    virtual void ComputeRightHandSide(VectComplex_wp& b_source, bool assemble = false) = 0;        

    virtual void GetCoefAB_Infinity( Real_wp& a_infty, Real_wp& b_infty) {}
  };
  

  //! class handling source with finite element classes
  template<class Dimension>
  class VarSourceProblem_Dim
  {
    typedef typename Dimension::R_N R_N;
    
  private :
    VarProblem<Dimension>& var_problem;
    VarSourceProblem_Base& var_source;
    
  public :    
    template<class TypeEquation>
    VarSourceProblem_Dim(EllipticProblem<TypeEquation>& var);

    void ReadIncidentWaveParam(const VectString& parameters, R_N& kwave, TinyVector<Complex_wp, Dimension::dim_N>&, R_N& origin, Real_wp& omega, Real_wp& w) const;

    virtual IncidentWaveField<Real_wp, Dimension>* GetNewIncidentField(int n, const Vector<VectString>& param,
								       const Real_wp& val) const;
    
    virtual IncidentWaveField<Complex_wp, Dimension>* GetNewIncidentField(int n, const Vector<VectString>& param,
									  const Complex_wp& val) const;

    virtual IncidentWaveField<Real_wp, Dimension>* GetIncidentField(int n, const Real_wp& val) const = 0;
    virtual IncidentWaveField<Complex_wp, Dimension>* GetIncidentField(int n, const Complex_wp& val) const = 0;

    virtual IncidentWaveProjector<Real_wp, Dimension>* GetIncidentWaveProjector(int n, const Real_wp&) const = 0;    
    virtual IncidentWaveProjector<Complex_wp, Dimension>* GetIncidentWaveProjector(int n, const Complex_wp&) const = 0;

    virtual void InitIncidentField() = 0;
    
    virtual void AddVolumeProjection(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				     Vector<VirtualProjectorFEM<Real_wp, Dimension>* >& f) const = 0;

    virtual void AddVolumeProjection(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				     Vector<VirtualProjectorFEM<Complex_wp, Dimension>* >& f) const = 0;
    
    virtual void AddSurfaceSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				  Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const = 0;

    virtual void AddSurfaceSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				  Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const = 0;

    virtual void SetSurfaceProjection(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				      Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const = 0;
    
    virtual void SetSurfaceProjection(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				      Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const = 0;

    virtual void AddVolumeSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				 Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const = 0;

    virtual void AddVolumeSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				 Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const = 0;
    
  };


  //! source class without dependence on the type of finite element
  template<class Complexe, class Dimension>
  class VarSourceProblem_Cplx : public VarSourceProblem_Dim<Dimension>
  {
  protected:
    Complexe coefficient_volumic_source;

    Vector<IncidentWaveProjector<Complexe, Dimension>* > incident_wave_proj;
    Vector<IncidentWaveField<Complexe, Dimension>* > incident_wave_field;
    
  private:
    VarProblem<Dimension>& var_problem;
    VarSourceProblem_Base& var_source;

  public:
    template<class TypeEquation>
    VarSourceProblem_Cplx(EllipticProblem<TypeEquation>&);

    ~VarSourceProblem_Cplx();
    
  protected:
    void ClearIncidentArray();
    
  public:
    // Inline methods
    Complexe GetCoefficientVolumeSource() const;

  protected:
    template<class T>
    IncidentWaveField<T, Dimension>* GetIncidentFieldGen(int n, const T& val) const;
    
    template<class T>
    IncidentWaveProjector<T, Dimension>* GetIncidentWaveProjectorGen(int n, const T&) const;

  public:
    virtual IncidentWaveField<Real_wp, Dimension>* GetIncidentField(int n, const Real_wp& val) const;
    virtual IncidentWaveProjector<Real_wp, Dimension>* GetIncidentWaveProjector(int n, const Real_wp&) const;

    virtual IncidentWaveField<Complex_wp, Dimension>* GetIncidentField(int n, const Complex_wp& val) const;
    virtual IncidentWaveProjector<Complex_wp, Dimension>* GetIncidentWaveProjector(int n, const Complex_wp&) const;

  protected:
    // other methods
    IncidentWaveField<Complexe, Dimension>* GetIncidentFieldGen(int n, const Complexe& val) const;
    IncidentWaveProjector<Complexe, Dimension>* GetIncidentWaveProjectorGen(int n, const Complexe&) const;

  public:
    virtual IncidentWaveProjector<Complexe, Dimension>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complexe, Dimension>& u_inc) const = 0;
        
    virtual void ConstructVolumeSourceFunctions(const Vector<VectString>&,
                                                VolumetricSource_Base<Complexe, Dimension>& var) const;

    virtual VirtualSourceField<Complexe, Dimension>* 
    GetNewVolumeSourceFunction(const IVect&, const VectString&, int&, Vector<Complexe>&,
                               VolumetricSource_Base<Complexe, Dimension>&) const;
    
    virtual void ConstructSurfaceSourceFunctions(const Vector<VectString>&,
						 VolumetricSource_Base<Complexe, Dimension>& var) const;

    virtual VirtualSourceField<Complexe, Dimension>* 
    GetNewSurfaceSourceFunction(const IVect&, const VectString&, int&, Vector<Complexe>&,
				VolumetricSource_Base<Complexe, Dimension>&) const;

    void FillVariableSource(const IVect& ref, typename Dimension::VectR_N& Points,
			    Matrix<Complexe>& val, Vector<Vector<Vector<Complexe> > >& eval) const;
        
    virtual ModalSourceBoundary_Dim<Complexe, Dimension>* GetNewModalSourceEquation() const = 0;
    virtual VirtualSourceFEM<Complexe, Dimension>* GetNewSourceEquationObject(int n) = 0;
    
  };
  
  
  //! class handling source with finite element classes
  template<class Dimension>
  class VarSourceProblem_Fem : public VarSourceProblem_Base
  {
  protected :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
  private :
    VarProblem<Dimension>& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    VarSourceProblem_Dim<Dimension>& var_source;
    
  public :
    template<class TypeEquation>
    VarSourceProblem_Fem(EllipticProblem<TypeEquation>& var);

    template<class T>
    void AddVolumetricSource(const T& alpha, Vector<Vector<T> > & b_source,
			     Vector<VirtualSourceFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void SetDirichletSource(Vector<Vector<T> >& b_source,
			    Vector<VirtualSourceFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void AddVolumetricProjection(const T& alpha, Vector<Vector<T> >& b_src,
				 Vector<VirtualProjectorFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void SetSurfacicProjection(const T& alpha, Vector<Vector<T> >& b_src,
			       Vector<VirtualSourceFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void AddSurfacicSource(const T& alpha, Vector<Vector<T> >& b_source,
			   Vector<VirtualSourceFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void AddDiracSourceGen(const T& alpha, Vector<Vector<T> > & b_source,
			   Vector<VirtualSourceFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void AddDiracSourceAxisym(const T& alpha, const Vector<bool>& vec_unknown, int m,
			      Vector<Vector<T> > & b_source, const VectR3& pt,
			      Vector<VirtualSourceFEM<T, Dimension>* >& f) const;
    
    template<class T>
    void AddDiracSourceGen(const T& alpha, const VectR_N& points_source,
			   const Vector<Vector<T> >& polar,
			   const IVect& point_number, Vector<Vector<T> > & b_source,
			   Vector<VirtualSourceFEM<T, Dimension>* >& f) const;

  public:
    template<class T>
    void ComputeGenericSource(Vector<Vector<T> > & b_source,
			      Vector<VirtualSourceFEM<T, Dimension>* >& f, bool assemble) const;

    // Inline methods
    virtual void AddDiracSource(const Real_wp&, Vector<VectReal_wp>&,
				Vector<VirtualSourceFEM<Real_wp, Dimension>* >&) const;
    
    virtual void AddDiracSource(const Complex_wp&, Vector<VectComplex_wp>&,
				Vector<VirtualSourceFEM<Complex_wp, Dimension>* >&) const;
    
  };
  

  template<class Complexe, class Dimension>
  class VarSourceProblem : public VarSourceProblem_Fem<Dimension>,
			   public VarSourceProblem_Cplx<Complexe, Dimension>
  {
  protected:
    typedef Vector<Complexe> VectComplexe;
        
  private:
    VarProblem<Dimension>& var_problem;
    
  public:
    template<class TypeEquation>
    VarSourceProblem(EllipticProblem<TypeEquation>&);
    
    // inline methods
    template<class T>
    void AddIncidentWaveGen(const T&, Vector<Vector<T> >&);

    virtual void AddIncidentWave(const Real_wp&, Vector<VectReal_wp>&);
    virtual void AddIncidentWave(const Complex_wp&, Vector<VectComplex_wp>&);
    
    // other methods
    virtual void InitIncidentField();
    
    void ComputeRightHandSideGen(VectComplexe & b_source, bool assemble = false);    
    virtual void ComputeRightHandSide(Vector<VectComplexe> & b_source, bool assemble = false);

    void ComputeRightHandSide(Matrix<Complexe, General, ColMajor> & b_source,
			      bool assemble = false);    
    
    void AddIncidentWaveGen(const Complexe&, Vector<VectComplexe>&);    

    void CopyInputData(const VarSourceProblem_Base& var);
    
    // inline methods
    template<class T>
    void ComputeRightHandSideGen(Vector<T>& b, bool assemble = false);
    
    virtual void ComputeRightHandSide(VectReal_wp& b_source, bool assemble = false);    
    virtual void ComputeRightHandSide(VectComplex_wp& b_source, bool assemble = false);        
    
    virtual void AddVolumeProjection(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				     Vector<VirtualProjectorFEM<Real_wp, Dimension>* >& f) const;

    virtual void AddVolumeProjection(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				     Vector<VirtualProjectorFEM<Complex_wp, Dimension>* >& f) const;

    virtual void AddSurfaceSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				  Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const;

    virtual void AddSurfaceSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				  Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const;

    virtual void SetSurfaceProjection(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				      Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const;

    virtual void SetSurfaceProjection(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				      Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const;

    virtual void AddVolumeSource(const Real_wp& alpha, Vector<VectReal_wp>& b_src,
				 Vector<VirtualSourceFEM<Real_wp, Dimension>* >& f) const;

    virtual void AddVolumeSource(const Complex_wp& alpha, Vector<VectComplex_wp>& b_src,
				 Vector<VirtualSourceFEM<Complex_wp, Dimension>* >& f) const;

  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_DEFINE_SOURCE_ELLIPTIC_HXX
#endif
