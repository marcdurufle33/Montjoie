#ifndef MONTJOIE_FILE_VAR_LAPLACE_HXX

namespace Montjoie
{
  
  //! class to specify the resolution of Laplace equation
  /*!
    The Laplace equation solved is
    \f$  \rho \, u  + \sigma u - \mbox{div}( \mu \nabla u) \, = \, f \f$
  */
  template<class Dimension>
  class LaplaceEquation : public HelmholtzEquation_Base<Real_wp, Dimension>
  {
  };
  
  
  
  //! base class to solve real Helmholtz equation
  template<class TypeEquation>
  class VarLaplace_Base : public VarHelmholtz_Eq<TypeEquation>
  {
  public :
    typedef typename TypeEquation::Dimension Dimension; //!< dimension
    typedef typename TypeEquation::Complexe Complexe;
    
    // Inline methods
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
                           IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
			   Vector<bool>& compute_hess);
        
  };


  //! class to compute the solution of wave equation outside the computational domain 
  //! with retarded potentials
  template<class Dimension>
  class VarComputationRCS_Laplace : public VarComputationRCS_Base<Dimension>
  {
  protected:
    //! u(t) will contribute to u^{far}(t') for t' belonging to the interval [t+tmin, t+tmax]
    Real_wp tmin, tmax;
    
    //! inverse of the wave speed at infinity
    /*!
      The media outside the computational domain is assumed to be homogeneous
      with a wave speed c_inf
      invCinfinity = 1/c_inf
    */
    Real_wp invCinfinity;
    
    //! if true values u(\tau) du/dn(\tau) are stored 
    /*!
      There are two ways to complete the computation of u^{far} :
      1- store u(\tau), du/dn(\tau) on the surface of integration
       for tau in [t - tmax, t-tmin]
       and make the computation of u^{far}(t+tmin), once all needed
       values of u(tau) and du/dn(tau) are available
       
       2- do not store u(tau), du/dn(tau), but store u^{far}(t')
        for t' in [t+tmin, t+tmax]. At each time step, the values
        u and du/dn are used to increment u^{far}(t').
        
        The computational cost is the same, but the second alternative
        is cheaper in memory if the number of points where u^{far} is computed
        is relatively small.
     */
    bool store_values_on_surface;
    
    //! values u(tau) on the surface of integration
    Vector<VectReal_wp> delayed_trace_Pn;
    //! values du/dt(tau) on the surface of integration
    Vector<VectReal_wp> delayed_trace_PnPoint;
    //! values du/dn(tau) on the surface of integration
    Vector<VectReal_wp> delayed_trace_dPnDn;
    
    //! u, du/dt, du_dn on the surface of integration for the current time t
    VectReal_wp Pn_current, PnPoint_current, dPnDn_current;
    
    //! stored values of u^{far}(t') for t' in [t+tmin, t+tmax]
    Vector<VectReal_wp> valU_outside;
    
    //! initial time of computation
    Real_wp t0;
    //! 1/(4 pi)
    Real_wp coef_u;
    //! initial time of computation
    Real_wp t0_start;
    //! u(t) is stored in delayed_trace_Pn(n0_current)
    int n0_current;
    //! number of time steps since the initial time
    int num_t0;
    
    //! times associated with delayed_trace_Pn
    VectReal_wp delayed_trace_time;
    
    //! if true, the L^2 norm \int p^{far}(t)^2 dt is computed
    bool compute_L2_norm;
    //! L^2 norm \int p^{far}(t)^2 dt
    VectReal_wp normL2;
    //! file where the L^2 norm is written
    string file_name_normL2;
    //! the interval of integration for L^2 norm is [t_begin_calculL2, t]
    Real_wp t_begin_calculL2;
    
    //! buffer used to write p^{far}
    WriteOnTheGoWithBuffer<Real_wp> output_buffer;
    //! size of the buffer
    int size_buffer;
    //! if true, the first writing of p^{far} on the file has not been done
    bool first_output_ext;

    typedef typename Dimension::R_N R_N;
    
  private:
    void InitDefaultValues();
    
    VarProblem<Dimension>& var_problem;
    VarOutputProblem<Dimension>& var_output;

  public:
    //! default constructor
    template<class TypeEquation>
    VarComputationRCS_Laplace(EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    void ComputeRCS(const VectReal_wp& );
    void ComputeRCS(const VectComplex_wp& );

    const string& GetFileNameNormL2() const;
    void SetFileNameNormL2(const string& nom);
    
    // Other methods
    void InitComputationRCS(bool assemble = false);    

    void LoadDatas(ParameterOutputReprise& );
    void SaveDatas(ParameterOutputReprise& );

    void WriteOutput(int nt, const Real_wp& t, const Real_wp& dt,
                     const VectReal_wp& Pn, const VectReal_wp& PnPoint);

    void SetL2NormOutput(bool compute_L2_norm_, const Real_wp& t_begin_calculL2_,
			 const string& file_name_normL2_);
    
    void ComputeTimeIntegralRepresentation(int nt, const Real_wp& t, const Real_wp& dt,
					   const MeshInterpolationFEM<Dimension>& mesh_, 
					   const R_N& pointX, Real_wp& scal_u) const;

    void UpdateTimeIntegralRepresentation(int nt, const Real_wp& t, const Real_wp& dt,
                                          const MeshInterpolationFEM<Dimension>& mesh_,
					  const R_N& pointX, VectReal_wp& scal_u) const;
    
  };

  
  //! class to compute the solution of wave equation outside the computational domain 
  //! with retarded potentials
  template<class Dimension>
  class VarComputationRCS<LaplaceEquation<Dimension> >
    : public VarComputationRCS_Laplace<Dimension>
  {
  public :
    template<class TypeEquation>
    VarComputationRCS(EllipticProblem<TypeEquation>& var);

  };
  

#ifdef MONTJOIE_WITH_TRANSMISSION
  //! class for transmission conditions for Laplace equation
  template<>
  class VarTransmission<LaplaceEquation<Dimension2> > 
    : public VarTransmission_Helm<Real_wp, Dimension2>
  {
  public :
    template<class TypeEquation>
    VarTransmission(EllipticProblem<TypeEquation>& var);

  };
#endif


  //! base class to compute incident field in time-domain and Laplace equation
  template<class Dimension>
  class IncidentWaveProjector_Laplace : public IncidentWaveProjector<Real_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    R_N k_wave; //!< wave vector
    Real_wp c0; //!< wave velocity
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    
    VirtualTimeSource<Real_wp>* wave_pulse;
    
  private:
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    IncidentWaveProjector_Laplace(const EllipticProblem<TypeEquation>& vars_helm,
				  IncidentWaveField<Real_wp, Dimension>&);
    
    void Init(const Real_wp& t0);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f);    
    void EvaluateFunctionGradient(int i, int j, const R_N& x, VectReal_wp& f, VectReal_wp& df);
    
  };
  
  
  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
   */
  template<class Dimension>
  class DiffractedWaveSource_Laplace : public VirtualSourceFEM<Real_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    R_N k_wave; //!< wave vector
    Real_wp c0; //!< wave velocity
    bool element_insidePML; //!< element of a PML layer ?
    Real_wp t; //!< current time
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    int num_deriv;
    bool scalar_eq;
    bool dirichlet_cond;
    
  private:
    IncidentWaveField<Real_wp, Dimension>& incident_wave;
    VirtualTimeSource<Real_wp>* wave_pulse;
    
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Real_wp, Dimension>& var_laplace;

  public :
    template<class TypeEquation>
    DiffractedWaveSource_Laplace(const EllipticProblem<TypeEquation>& vars_helm,
                                 IncidentWaveField<Real_wp, Dimension>&);

    Real_wp Init(const Real_wp&, const Real_wp&, int p, int n, bool);

    bool PresenceDirichlet() const;
    void EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f);    
    
    void InitElement(int num_elem, const VectR_N& s);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    bool IsNonNullGradientSource(const VectR_N& s);
    
    void EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f);    
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  //! source for the scattering of an incident wave
  template<class Dimension>
  class DiffractedWaveSource<LaplaceEquation<Dimension> >
    : public DiffractedWaveSource_Laplace<Dimension>
  {
  public :
    //! constructor with given problem
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Real_wp, Dimension>& fsrc_);

  };

  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  template<class Dimension>
  class TotalWaveSource_Laplace : public VirtualSourceFEM<Real_wp, Dimension>
  {
  protected:
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    R_N k_wave; //!< wave vector
    Real_wp c0; //!< wave velocity
    bool element_insidePML; //!< element of a PML layer ?
    Real_wp t; //!< current time
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    int num_deriv;
    bool scalar_eq;
    
  private:
    IncidentWaveField<Real_wp, Dimension>& incident_wave;
    VirtualTimeSource<Real_wp>* wave_pulse;
    
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Real_wp, Dimension>& var_laplace;

  public :
    template<class TypeEquation>
    TotalWaveSource_Laplace(const EllipticProblem<TypeEquation>& vars_helm,
			    IncidentWaveField<Real_wp, Dimension>&);
    
    Real_wp Init(const Real_wp&, const Real_wp&, int p, int n, bool);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f);
    
    void InitElement(int num_elem, const VectR_N& s);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    bool IsNonNullGradientSource(const VectR_N& s);
    
    void EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  
  //! source for the scattering of an incident wave
  template<class Dimension>
  class TotalWaveSource<LaplaceEquation<Dimension> >
    : public TotalWaveSource_Laplace<Dimension>
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Real_wp, Dimension>& fsrc);
  };

  
  //! class used to solve Poisson equation with continuous elements
  template<class Dimension>
  class EllipticProblem<LaplaceEquation<Dimension> >
    : public VarLaplace_Base<LaplaceEquation<Dimension> >
  {
  public :
    IncidentWaveProjector<Real_wp, Dimension>*
    GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Real_wp, Dimension>& u_inc) const;
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

  };
  

  template<class Dimension>
  class FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Real_wp, LaplaceEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<LaplaceEquation<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  template<class Dimension>
  class FemMatrixFreeClass<Complex_wp, LaplaceEquation<Dimension> >
    : public FemMatrixFreeClass_Eq<Complex_wp, LaplaceEquation<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<LaplaceEquation<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  template<class T, class Dimension>
  class CondensationBlockSolver<T, LaplaceEquation<Dimension> >
    : public CondensationBlockSolver_Helm<T, LaplaceEquation<Dimension> >
  {
  public:
    CondensationBlockSolver(EllipticProblem<LaplaceEquation<Dimension> >& var);
    
  };
  
  
  //! class to specify the resolution of real Helmholtz equation with dg method
  /*!
    The Helmholtz equation solved is
    \f$ \rho \, u - \mbox{div}( v) \, = \, f \f$
    \f$ \frac{1}{\mu}  v + \mbox{grad}(v) \, = \, g \f$
   */
  template<class Dimension>
  class LaplaceEquationDG : public HelmholtzEquationDG_Base<Real_wp, Dimension>
  {
  public :
    template<class TypeEquation, class T, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& var,
			      int i, int j, const GlobalGenericMatrix<T>& nat_mat, int ref, MatMass& mass);
    
    template<class TypeEquation, class T, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int k,
				const GlobalGenericMatrix<T>& nat_mat, int ref, Vector1& Un, Vector1& Vn);
    
    template<class T, class Matrix1, class TypeEquation>
    static void GetNabc(Matrix1& Nabc, typename Dimension::R_N& normale,
			int ref, int iquad, int k, const GlobalGenericMatrix<T>& nat_mat, int ref_d,
			const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb);
    
    template<class Matrix1, class TypeEquation, class T>
    static void GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad,
			   int k, int nf, const GlobalGenericMatrix<T>& nat_mat, int ref, int ref2,
			   const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);
    
    template<class Vector1, class T, class TypeEquation>
    static void MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector1& Un,
			   int i, int n, int nf, const GlobalGenericMatrix<T>& nat_mat,
			   int ref, int ref2, const EllipticProblem<TypeEquation>& vars,
			   const ElementReference<Dimension, 1>& Fb);
  
    template<class Vector1, class T, class TypeEquation>
    static void MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
			int num_elem1, int num_point, const GlobalGenericMatrix<T>& nat_mat, int ref_d,
			const EllipticProblem<TypeEquation>& vars,
			const ElementReference<Dimension, 1>& Fb);
    
  };


  //! base class to compute incident field in time-domain and Laplace equation
  template<class Dimension>
  class IncidentWaveProjector_LaplaceDG : public IncidentWaveProjector_Laplace<Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;

  public :
    template<class TypeEquation>
    IncidentWaveProjector_LaplaceDG(const EllipticProblem<TypeEquation>& var,
				    IncidentWaveField<Real_wp, Dimension>&);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f);    
    void EvaluateFunctionGradient(int i, int j, const R_N& x, VectReal_wp& f, VectReal_wp& df);
    
  };


  //! source for the scattering of an incident wave
  /*!
    the diffracted field is computed
   */
  template<class Dimension>
  class DiffractedWaveSource_LaplaceDG : public VirtualSourceFEM<Real_wp, Dimension>
  {
  protected:
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    R_N k_wave; //!< wave vector
    Real_wp c0; //!< wave velocity
    bool element_insidePML; //!< element of a PML layer ?
    Real_wp t; //!< current time
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    int num_deriv;
    int reference_condition;
    bool scalar_eq;

  private:
    IncidentWaveField<Real_wp, Dimension>& incident_wave;
    VirtualTimeSource<Real_wp>* wave_pulse;
    
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Real_wp, Dimension>& var_laplace;
    
  public :
    template<class TypeEquation>
    DiffractedWaveSource_LaplaceDG(const EllipticProblem<TypeEquation>& vars_helm,
				   IncidentWaveField<Real_wp, Dimension>&);
    
    Real_wp Init(const Real_wp&, const Real_wp&, int p, int n, bool);
    
    void InitElement(int num_elem, const VectR_N& s);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    bool IsNonNullGradientSource(const VectR_N& s);
    
    void EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                                const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };
  
  //! source for the scattering of an incident wave
  template<class Dimension>
  class DiffractedWaveSource<LaplaceEquationDG<Dimension> >
    : public DiffractedWaveSource_LaplaceDG<Dimension>
  {
  public :
    template<class TypeEquation>
    DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Real_wp, Dimension>& fsrc_);
  };

  
  //! source for the scattering of an incident wave
  /*!
    the total field is computed
   */
  template<class Dimension>
  class TotalWaveSource_LaplaceDG : public VirtualSourceFEM<Real_wp, Dimension>
  {
  protected :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    
    R_N k_wave; //!< wave vector
    Real_wp c0; //!< wave velocity
    bool element_insidePML; //!< element of a PML layer ?
    Real_wp t; //!< current time
    Real_wp offset_arg; //!< f will be evaluated for k.x - c t - offset_arg
    int num_deriv;
    bool scalar_eq;
    
  private:
    IncidentWaveField<Real_wp, Dimension>& incident_wave;
    VirtualTimeSource<Real_wp>* wave_pulse;
    
    const VarSourceProblem_Base& var_source;
    const VarBoundaryCondition_Base& var_boundary;
    const VarProblem<Dimension>& var_problem;
    const VarHelmholtz_Base<Real_wp, Dimension>& var_laplace;
    
  public :
    template<class TypeEquation>
    TotalWaveSource_LaplaceDG(const EllipticProblem<TypeEquation>& vars_helm,
			      IncidentWaveField<Real_wp, Dimension>&);
    
    Real_wp Init(const Real_wp&, const Real_wp&, int p, int n, bool);
    
    void EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f);
    
    void InitElement(int num_elem, const VectR_N& s);
    bool IsNonNullVolumetricSource(const VectR_N& s);
    
    void EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    bool IsNonNullGradientSource(const VectR_N& s);
    
    void EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f);
    
    void InitSurface(int i, int num_face, int num_elem, int num_loc);
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				 const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f);
    
    bool IsNonNullSurfacicSource(int cond_ref);
    
  };

  
  //! source for the scattering of an incident wave
  template<class Dimension>
  class TotalWaveSource<LaplaceEquationDG<Dimension> >
    : public TotalWaveSource_LaplaceDG<Dimension>
  {
  public :
    template<class TypeEquation>
    TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Real_wp, Dimension>& fsrc_);    
  };


  template<class Dimension>
  class VolumetricSource<LaplaceEquationDG<Dimension> >
    : public VolumetricSource_HelmDG<Real_wp, Dimension>
  {
  public :    
    VolumetricSource(const EllipticProblem<LaplaceEquationDG<Dimension> >&,
		     const Vector<VectString>& param);
    
  };

  
  //! class for first-order absorbing boundary condition
  template<class Complexe, class Dimension>
  class ImpedanceABC_Laplace
    : public ImpedanceFunction_Base<Complexe, Dimension>
  {
  private:
    const VarHelmholtz_Base<Real_wp, Dimension>& var_helm;
    const VarProblem<Dimension>& var_problem;

  public :
    template<class TypeEquation>
    ImpedanceABC_Laplace(const EllipticProblem<TypeEquation>& var);
    
    //! impedance
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
                              const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat);  
  };


  //! class for first-order absorbing boundary condition for Helmholtz equation
  template<class Complexe, class Dimension>
  class ImpedanceABC<Complexe, LaplaceEquation<Dimension> > 
    : public ImpedanceABC_Laplace<Complexe, Dimension>
  {    
  public :
    ImpedanceABC(const EllipticProblem<LaplaceEquation<Dimension> >& var);
    
  };


  //! class for first-order absorbing boundary condition for Helmholtz equation
  template<class Complexe, class Dimension>
  class ImpedanceABC<Complexe, LaplaceEquationDG<Dimension> > 
    : public ImpedanceABC_Laplace<Complexe, Dimension>
  {    
  public :
    ImpedanceABC(const EllipticProblem<LaplaceEquationDG<Dimension> >& var);
    
  };

  
  template<class T, class Dimension>
  class CondensationBlockSolver<T, LaplaceEquationDG<Dimension> >
    : public CondensationBlockSolver_HelmDG<T, Real_wp, Dimension>
  {
  public:
    CondensationBlockSolver(EllipticProblem<LaplaceEquationDG<Dimension> >& var);
    
  };


  //! class used to solve Poisson equation with DG method
  template<class Dimension>
  class EllipticProblem<LaplaceEquationDG<Dimension> >
    : public VarLaplace_Base<LaplaceEquationDG<Dimension> >
  {    
  public:
    IncidentWaveProjector<Real_wp, Dimension>*
    GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Real_wp, Dimension>& u_inc) const;

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);
    
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
    void AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			       const GlobalGenericMatrix<Real_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

    void AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			       const GlobalGenericMatrix<Complex_wp>& nat_mat,
			       int offset_row = 0, int offset_col = 0);

  };    


  template<class Dimension>
  class FemMatrixFreeClass<Real_wp, LaplaceEquationDG<Dimension> >
    : public FemMatrixFreeClass_Eq<Real_wp, LaplaceEquationDG<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<LaplaceEquationDG<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };


  template<class Dimension>
  class FemMatrixFreeClass<Complex_wp, LaplaceEquationDG<Dimension> >
    : public FemMatrixFreeClass_Eq<Complex_wp, LaplaceEquationDG<Dimension> >
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<LaplaceEquationDG<Dimension> >& var_);
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };

}

#define MONTJOIE_FILE_VAR_LAPLACE_HXX
#endif

