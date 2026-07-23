#ifndef MONTJOIE_FILE_SOURCE_SPATIALE_HXX

namespace Montjoie
{

  //! Abstract class for a source
  template<class T, class Dimension>
  class VirtualSourceField
  {
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    virtual ~VirtualSourceField();

    // Pure virtual methods that need to be overloaded
    virtual void EvaluateFunction(const R_N& x, Vector<T>& f) const = 0;
    virtual void EvaluateGradient(const R_N& x, Vector<T>& df) const;
    
  };
  
  
  //! gaussian source around a point
  /*!
    the considered function is equal to
    \f$ f(r) = beta \mbox{exp}(-alpha r^2) \f$
    where \f$ alpha = -\frac{\log(1e-6)}{r_0^2} \f$
    $r_0$ is the radius of the gaussian.
    r is the distance to the source point
  */
  template<class Dimension>
  class GaussianSource
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    enum{dim_N = Dimension::dim_N};
    
    Real_wp RadiusSource; //!< radius of gaussian
    Real_wp RadiusSource2; //!< square of the radius
    Real_wp alpha, beta;
    Real_wp RadiusSourceCutOff;
    R_N origin;
    
  public :
    GaussianSource();
    
    GaussianSource(const R_N&, const Real_wp&);
        
    Real_wp GetRadius() const;
    Real_wp GetCutOffRadius() const;
    const R_N& GetCenter() const;
    
    void Init(const R_N&, const Real_wp&, const Real_wp&);
    
    Real_wp GetAmplitude(const R_N& x) const;
    void GetGradAmplitude(const R_N& x, Real_wp&, R_N&) const;

    void GetHessianAmplitude(const R_N& x, Real_wp& aj, R_N& grad_aj,
                             TinyMatrix<Real_wp, Symmetric, dim_N, dim_N>& hess_aj) const;
  };
  
  
  //! Field with a gaussian source
  template<class T, class Dimension>
  class GaussianSourceField : public GaussianSource<Dimension>,
			      public VirtualSourceField<T, Dimension>
  {
  protected:
    Vector<T> polar;
    
  public :    
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    void SetPolarization(const Vector<T>& P);
    void EvaluateFunction(const R_N& x, Vector<T>& f) const;    
    void EvaluateGradient(const R_N& x, Vector<T>& df) const;
    
  };


  //! Uniform field
  template<class T, class Dimension>
  class UniformSourceField : public VirtualSourceField<T, Dimension>
  {
  protected:
    Vector<T> polar;
    
  public :    
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    void SetPolarization(const Vector<T>& P);
    void EvaluateFunction(const R_N& x, Vector<T>& f) const;    
    void EvaluateGradient(const R_N& x, Vector<T>& df) const;
    
  };
  

  //! Abstract class for an incident field
  template<class T, class Dimension>
  class IncidentWaveField
  {
  protected:
    //! pulse in time 
    VirtualTimeSource<T>* time_pulse;
    
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3

    IncidentWaveField();
    virtual ~IncidentWaveField();

    VirtualTimeSource<T>* GetTimeSource();
    void SetTimeSource(VirtualTimeSource<T>* fsrc);
    
    virtual void InitElement(int i, const typename Dimension::VectR_N& s);

    // Pure virtual methods that need to be overloaded
    virtual void EvaluateFunction(const R_N& x, T& f) const = 0;
    virtual void EvaluateFunctionGradient(const R_N& x, T& f,
					  TinyVector<T, Dimension::dim_N>& df) const = 0;
    
    virtual void UpdateCoefAB(T& a, T& b) const;
    
  };
  
  
  //! empty class, specialized
  template<class T, class Dimension>
  class PlaneWaveIncidentField
  {
  };
  
  
  //! not implemented for real numbers
  template<class Dimension>
  class PlaneWaveIncidentField<Real_wp, Dimension>
    : public IncidentWaveField<Real_wp, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    PlaneWaveIncidentField(const R_N&, const R_N&);
    
    void EvaluateFunction(const R_N& x, Real_wp& f) const;
    void EvaluateFunctionGradient(const R_N& x, Real_wp& f,
				  TinyVector<Real_wp, Dimension::dim_N>& df) const;
    
  };
  
  
  //! plane wave incident field
  /*!
    the function f is equal to \f$ f = \mbox{exp}( i k \cdot (x-x_0) ) \f$
    where k is the wave vector, and x_0 the origine of the phase     
  */
  template<class Dimension>
  class PlaneWaveIncidentField<Complex_wp, Dimension>
    : public IncidentWaveField<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    R_N origin; //!< origin of the plane wave
    R_N k_wave; //!< wave vector
    
  public :
    PlaneWaveIncidentField(const R_N&, const R_N&);

    void Init(const R_N&, const R_N&);
    
    void EvaluateFunction(const R_N& pt, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R_N& pt, Complex_wp& f,
				  TinyVector<Complex_wp, Dimension::dim_N>& df) const;
    
  };


  //! empty class, specialized
  template<class T, class Dimension>
  class PlaneWaveIncidentFieldComplex
  {
  };


  //! not implemented for real numbers
  template<class Dimension>
  class PlaneWaveIncidentFieldComplex<Real_wp, Dimension>
    : public IncidentWaveField<Real_wp, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    PlaneWaveIncidentFieldComplex(const R_N&, const TinyVector<Complex_wp, Dimension::dim_N>&);
    
    void EvaluateFunction(const R_N& x, Real_wp& f) const;
    void EvaluateFunctionGradient(const R_N& x, Real_wp& f,
				  TinyVector<Real_wp, Dimension::dim_N>& df) const;
    
  };
  
  
  //! plane wave incident field
  /*!
    the function f is equal to \f$ f = \mbox{exp}( i k \cdot (x-x_0) ) \f$
    where k is the wave vector, and x_0 the origine of the phase     
  */
  template<class Dimension>
  class PlaneWaveIncidentFieldComplex<Complex_wp, Dimension>
    : public IncidentWaveField<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    R_N origin; //!< origin of the plane wave
    TinyVector<Complex_wp, Dimension::dim_N> k_wave; //!< wave vector
    
  public :
    PlaneWaveIncidentFieldComplex(const R_N&, const TinyVector<Complex_wp, Dimension::dim_N>&);
    
    void Init(const R_N&, const TinyVector<Complex_wp, Dimension::dim_N>&);
    
    void EvaluateFunction(const R_N& pt, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R_N& pt, Complex_wp& f,
				  TinyVector<Complex_wp, Dimension::dim_N>& df) const;
    
  };


  //! empty class, specialized
  template<class T, class Dimension>
  class HankelIncidentField
  {
  };


  //! not implemented for real numbers
  template<class Dimension>
  class HankelIncidentField<Real_wp, Dimension>
    : public IncidentWaveField<Real_wp, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3

    HankelIncidentField(const R_N&, const R_N&);
    
    void EvaluateFunction(const R_N& x, Real_wp& f) const;
    void EvaluateFunctionGradient(const R_N& x, Real_wp& f,
				  TinyVector<Real_wp, Dimension::dim_N>& df) const;
    
  };
  
  
  //! Hankel incident field (spherical wave in 3-D)
  template<class Dimension>
  class HankelIncidentField<Complex_wp, Dimension>
    : public IncidentWaveField<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    R_N origin;
    R_N k_wave; //!< wave vector
    Real_wp norm_kwave; //!< wave number
    
  public :
    HankelIncidentField(const R_N&, const R_N&);    
    
    void Init(const R_N&, const R_N&);    

    void EvaluateFunction(const R_N& pt, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R_N& pt, Complex_wp&,
				  TinyVector<Complex_wp, Dimension::dim_N>& f) const;
    
  };
  

  //! empty class, overloaded
  template<class T, class Dimension>
  class GaussianBeamIncidentField
  {
  };


  //! not implemented for real numbers
  template<class Dimension>
  class GaussianBeamIncidentField<Real_wp, Dimension>
    : public IncidentWaveField<Real_wp, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N; //!< R2 or R3

    GaussianBeamIncidentField(const R_N&, const R_N&, const Real_wp&);
    
    void EvaluateFunction(const R_N& x, Real_wp& f) const;
    void EvaluateFunctionGradient(const R_N& x, Real_wp& f,
				  TinyVector<Real_wp, Dimension::dim_N>& df) const;
    
  };

  
#ifdef MONTJOIE_WITH_TWO_DIM
  //! gaussiam beam
  /*!
    gaussian beam, with direction of propagation (k_wave)
    and width of the beam (waist)
  */
  template<>
  class GaussianBeamIncidentField<Complex_wp, Dimension2>
    : public IncidentWaveField<Complex_wp, Dimension2>
  {
  protected :
    R2 origin;
    R2 k_wave; //!< wave vector
    Real_wp omega; //!< pulsation
    Real_wp waist; //!< waist of the gaussian beam
    
  public :
    GaussianBeamIncidentField(const R2& pt, const R2& u, const Real_wp&);
    
    void Init(const R2& pt, const R2& u, const Real_wp&);
    
    void EvaluateFunction(const R2& x, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R2& x, Complex_wp& f,
				  TinyVector<Complex_wp, 2>& df) const;
    
  };
#endif
#ifdef MONTJOIE_WITH_THREE_DIM
  //! gaussian beam
  /*! gaussian beam, with direction of propagation (k_wave)
    and width of the beam (waist_x, waist_y)
  */
  template<>
  class GaussianBeamIncidentField<Complex_wp, Dimension3>
    : public IncidentWaveField<Complex_wp, Dimension3>
  {
  protected :
    R3 origin;
    R3 k_wave; //!< wave vector
    R3 k_unit; //!< unit vector with same orientation than k_wave
    Real_wp omega; //!< pulsation
    Real_wp waist; //!< waist of the gaussian beam
    
  public :
    GaussianBeamIncidentField(const R3& pt, const R3& u, const Real_wp&);
    
    void Init(const R3& pt, const R3& u, const Real_wp&);
    
    void EvaluateFunction(const R3& x, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R3& x, Complex_wp& f,
				  TinyVector<Complex_wp, 3>& df) const;
    
  };
#endif

  
  //! empty class, specialized
  template<class T, class Dimension>
  class LayeredPlaneWaveIncidentField
  {
  };


  //! plane wave incident field in layered medium
  /*!
  */
  template<class Dimension>
  class LayeredPlaneWaveIncidentField<Complex_wp, Dimension>
    : public IncidentWaveField<Complex_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    R_N origin; //!< origin of the plane wave
    VectComplex_wp A, B;
    Real_wp kx;
    VectComplex_wp ky;
    VectComplex_wp a, b; VectReal_wp d;
    Real_wp omega;
    int num_layer;

  public :
    LayeredPlaneWaveIncidentField(const R_N& pt0, const R_N& kwave, const Real_wp& omega,
                                  const Real_wp& a_infty, const Real_wp& b_infty,
                                  const VectString& parameters);
    
    void Init(const R_N& pt0, const R_N& kwave, const Real_wp& omega,
              const Real_wp& a_infty, const Real_wp& b_infty,
	      const VectString& parameters);
    
    void InitElement(int i, const typename Dimension::VectR_N& s);
    
    void EvaluateFunction(const R_N& pt, Complex_wp& f) const;
    void EvaluateFunctionGradient(const R_N& pt, Complex_wp& f,
				  TinyVector<Complex_wp, Dimension::dim_N>& df) const;
    
    void UpdateCoefAB(Complex_wp& a, Complex_wp& b) const;
    
  };

  
  //! Field given as u_inc polarization where u_inc is as incident wave
  template<class T, class Dimension>
  class IncidentWaveSourceField : public VirtualSourceField<T, Dimension>
  {
  protected:
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    
    IncidentWaveField<T, Dimension>* incident_wave;
    Vector<T> polar;
    
  public :    
    IncidentWaveSourceField(IncidentWaveField<T, Dimension>*);
    ~IncidentWaveSourceField();
    
    void SetPolarization(const Vector<T>& P);
    void EvaluateFunction(const R_N& x, Vector<T>& f) const;    
    void EvaluateGradient(const R_N& x, Vector<T>& df) const;
    
  };


  template<class T>
  class VirtualProjectorFEM_Base
  {
  protected:
    //! object for ffts
    FftInterface<Complex_wp> fft_interface;
    
    void ApplyFFT_Point(int nx, int ny, int nz, Vector<Vector<Vector<Vector<T> > > >& feval);
    void ApplyFFT_Global(Vector<Vector<Vector<Vector<T> > > >& feval);
    
  };
    
  
  //! base class to project a function in finite element space
  template<class T, class Dimension>
  class VirtualProjectorFEM : public VirtualProjectorFEM_Base<T>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    
    //! number of unknowns
    int nb_unknowns;
    
  private :
    const VarBoundaryCondition_Base& var_boundary;
    const DistributedProblem<Dimension>& var_problem;
    
  public :    
    template<class TypeEquation>
    VirtualProjectorFEM(const EllipticProblem<TypeEquation>&);
    
    virtual ~VirtualProjectorFEM();

    int GetNbUnknowns() const;
    virtual size_t GetMemorySize() const;
    
    virtual void InitElement(int i, const VectR_N& s) = 0;

    virtual void EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f) = 0;

    void ApplyFFT(Vector<Vector<Vector<Vector<T> > > >& feval);    
    
    void ModifyPoints(int nsrc, VectR_N& s, SetPoints<Dimension>& Pts,
		      SetMatrices<Dimension>& Mat);

    void ModifyEvaluationProjection(Vector<Vector<Vector<Vector<T> > > >& feval,
                                    Vector<bool>& is_f_vec);
    
  };
  

  //! Class used to project plane wave into finite element class
  template<class T, class Dimension>
  class IncidentWaveProjector : public VirtualProjectorFEM<T, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    
    IncidentWaveField<T, Dimension>& incident_wave;
    Vector<T> polarization;
    Real_wp t; //!< current time
    
  public :
    template<class TypeEquation>
    IncidentWaveProjector(const EllipticProblem<TypeEquation>&,
			  IncidentWaveField<T, Dimension>&);

    virtual ~IncidentWaveProjector();

    const IncidentWaveField<T, Dimension>& GetIncidentWaveField() const;

    virtual void Init(const Real_wp& t0);
    
    void SetPolarization(const Vector<T>& p);
    void SetTime(const Real_wp& tn, int n);
    
    void InitElement(int i, const typename Dimension::VectR_N& s);
    void EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f);
    
    virtual void EvaluateFunctionGradient(int i, int j, const R_N& x, Vector<T>& f, Vector<T>& df);
    
  };

  
  //! base class for sources in class EllipticProblem
  /*!
    In order to define a new type of source for elliptic problems,
    the user can derive this class, and overload some methods.
    By default these methods return a null source
  */
  template<class T, class Dimension>
  class VirtualSourceFEM : public VirtualProjectorFEM<T, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points

    //! referenced edge number (in array edges_ref or faces_ref)
    int num_bound_ref_; 
    int num_boundary_; //!< edge number (in array faces or edges)
    int num_loc_; //!< local number of the boundary
    int num_elem_; //!< element number
    int ref_boundary_; //!< reference of the edge/face

    Vector<T> polarization, polarization_grad;
    bool grad_dirac;
    
  private:
    const DistributedProblem<Dimension>& var_problem;
    const VarBoundaryCondition_Base& var_boundary;
    
  public :
    template<class TypeEquation>
    VirtualSourceFEM(const EllipticProblem<TypeEquation>&);

    virtual ~VirtualSourceFEM();
    
    virtual const Vector<T>& GetPolarization() const;
    virtual const Vector<T>& GetPolarizationGrad() const;
    virtual const R_N& GetOrigin() const;
    
    void SetPolarization(const Vector<T>&);
    void SetPolarizationGrad(const Vector<T>&);
    bool IsGradientDirac() const;
    
    // virtual methods that can be overloaded
    virtual Real_wp Init(const Real_wp& t, const Real_wp& dt, int print_level,
			 int nb_deriv, bool scalar_eq);
    
    virtual void InitElement(int i, const VectR_N& s);
    virtual void EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f);
    
    virtual void EvaluateVolumetricSource(int i, int j, const R_N& x, Vector<T>& f);
    virtual bool IsNonNullVolumetricSource(const VectR_N& s);
    
    virtual void EvaluateGradientSource(int i, int j, const R_N& x, Vector<T>& f);
    virtual bool IsNonNullGradientSource(const VectR_N& s);

    virtual void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    virtual bool IsNonNullSurfacicSource(int cond_ref);
    virtual void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
					const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);

    virtual bool IsNonNullSurfacicSourceGradient(int cond_ref);
    virtual void EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
						const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);
    
    virtual bool IsNonNullDirichletSource(int cond_ref);
    virtual bool PresenceDirichlet() const;
    
    virtual bool IsDiracSource() const;
    virtual T GetCoefficientVolume() const;
    
    void ModifyEvaluationVolume(Vector<bool>& int_phi, Vector<bool>& int_grad,
				Vector<Vector<Vector<Vector<T> > > >& feval,
				Vector<Vector<Vector<Vector<T> > > >& feval_grad,
                                Vector<bool>&, Vector<bool>&);

    void ModifyEvaluationSurface(bool& int_phi, bool& int_grad,
				 Vector<Vector<Vector<Vector<T> > > >& feval,
                                 Vector<Vector<Vector<Vector<T> > > >& feval_diff,
                                 Vector<bool>&, Vector<bool>&);
    
  };  
  
  
  //! General source for EllipticProblem object
  /*!
    The right hand side is formed with terms
    \int_\Omega f \varphi dx + \int_\Omega g \nabla \varphi dx + \int_\Gamma f_surf \varphi ds
    An inhomogeneous Dirichlet condition can also be set u = h
   */
  template<class T, class Dimension>
  class VolumetricSource_Base : public VirtualSourceFEM<T, Dimension>
  {
  private :
    typedef typename Dimension::R_N R_N; //!< R2 or R2
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    
    //! object representing f (for each reference)
    Vector<VirtualSourceField<T, Dimension>* > f_vol;

    //! object representing g (for each reference)
    Vector<VirtualSourceField<T, Dimension>* > g_vol;

    //! surface source (for each reference) (used both for f_surf or h)
    Vector<VirtualSourceField<T, Dimension>* > f_surf;

    //! surface source (for each reference) for gradient terms
    Vector<VirtualSourceField<T, Dimension>* > g_surf;
    
    //! reference to the object VarSourceProblem
    const VarSourceProblem_Base& var_source;
    const DistributedProblem<Dimension>& var_problem;
    
    Vector<bool> variable_fvol, variable_gvol;
    Vector<bool> variable_fsurf;
    
    void InitDefaultValues();

  public :
    //! evaluation of f on quadrature points for \f$ \int_\Omega f \varphi \f$
    Vector<Vector<Vector<T> > > evalS;
    //! evaluation of g on quadrature points for \f$ \int_\Omega g \nabla \varphi \f$ 
    Vector<Vector<Vector<T> > > evalG;
    //! evaluation of f_surf on quadrature points for \f$ \int_\Gamma f_{surf} \varphi \f$
    Vector<Vector<Vector<T> > > evalSurf;
    
  public :
    T coef_vol; //!< coefficient
    
    template<class TypeEquation>
    VolumetricSource_Base(const EllipticProblem<TypeEquation>& var,
			  const Vector<VectString>& param);
    
    ~VolumetricSource_Base();

    size_t GetMemorySize() const;
    
    void SetVolumeSource(const IVect& ref, VirtualSourceField<T, Dimension>* f);
    void SetVolumeSourceGrad(const IVect& ref, VirtualSourceField<T, Dimension>* f);
    void SetVolumeSourceFunction(VirtualSourceField<T, Dimension>& f);
    void NullifyVolumeSourceFunction();

    void SetModalSource(const IVect& ref);
    void SetSurfaceSource(const IVect& ref, VirtualSourceField<T, Dimension>* f);
    void SetSurfaceSourceGrad(const IVect& ref, VirtualSourceField<T, Dimension>* f);
    
    void SetVariableSource(const IVect& ref);
    void SetVariableSurfaceSource(const IVect& ref);
    void SetVariableGradientSource(const IVect& ref);
    
    bool IsNonNullVolumetricSource(const VectR_N& s);
    void EvaluateVolumetricSource(int i, int j, const R_N& x, Vector<T>& f);

    bool IsNonNullGradientSource(const VectR_N& s);    
    void EvaluateGradientSource(int i, int j, const R_N& x, Vector<T>& f);

    bool IsNonNullDirichletSource(int cond_ref);        
    void EvaluateFunction(int i, int j, const R_N& x, Vector<T>& f);

    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);

    bool IsNonNullSurfacicSourceGradient(int cond_ref);
    void EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
                                        const SetMatrices<Dimension>& MatricesElem, Vector<T>& f);
    
    T GetCoefficientVolume() const;
    
  };


  template<class TypeEquation>
  class VolumetricSource : public VolumetricSource_Base<typename TypeEquation::Complexe,
							typename TypeEquation::Dimension>
  {
  public :
    VolumetricSource(const EllipticProblem<TypeEquation>&,
		     const Vector<VectString>& param);
    
  };
  
  
  //! Dirac source
  template<class TypeEquation>
  class DiracSource : public VirtualSourceFEM<typename TypeEquation::Complexe,
					      typename TypeEquation::Dimension>
  {
  protected:
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;

    typename Dimension::R_N origin;
    Complexe coef_vol;
    
  public :
    DiracSource(const EllipticProblem<TypeEquation>& var,
		const Vector<VectString>& param);
    
    bool IsDiracSource() const;
    
    const typename Dimension::R_N& GetOrigin() const;

    Complexe GetCoefficientVolume() const;
    
  };
  
}

#define MONTJOIE_FILE_SOURCE_SPATIALE_HXX
#endif
