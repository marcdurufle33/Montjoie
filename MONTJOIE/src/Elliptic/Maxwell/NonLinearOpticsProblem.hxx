#ifndef MONTJOIE_FILE_NON_LINEAR_OPTICS_PROBLEM_HXX

namespace Montjoie
{
  
  //! base class for Raman effect
  class RamanEffect
  {
  protected :
    //! characteristic times
    Real_wp tau1, tau2, Tmax;
    //! part in Raman effect
    Real_wp alpha;
    Real_wp coef, invTau1, invTau2;
    
    // for manual convolution
    //! stored values of |E|^2 at different times
    Vector<VectReal_wp> moduleE;
    Real_wp deltat, length_interval;
    int num_current; Real_wp t_current;
    int level_romberg, level_global;
    VectReal_wp hEval;
    
  public :
    
    RamanEffect();
    
    void Init(const Real_wp&, const Real_wp&,
	      const Real_wp&, const Real_wp&);

    Real_wp GetAlpha() const;
    
    Real_wp EvaluateH(const Real_wp& t) const;
    
    void GetCoefPDE(Real_wp&, Real_wp&, Real_wp&) const;
    
    // for manual convolution
    void Construct(int N, const Real_wp& dt);
    
    template<class T>
    void StoreNewField(const Real_wp& t, const Vector<T>& E);

    int GetNumberOfIterates() const;
    VectReal_wp& GetIterate(int k);
    
    void InitCoefficient(const Real_wp& t);
    void UpdateCoefficient(const Real_wp& t, int i, Real_wp& gamma, Real_wp&, Real_wp&);
    
  };

  
  //! base class for non-linear optics
  class NonLinearOpticsProblem
  {
  public :
    //! Newton solver used to solve the non-linear system
    NewtonSolver<Complex_wp> newton_solver;
    
    //! non-linear coefficient (Kerr effect)
    Real_wp gammaNL;

    //! object for Raman effect
    RamanEffect var_raman;
    VectComplex_wp raman_fct_hat;
    
    //! Source for Newton solver
    /*!
      The non-linear system is written in the form f(E^{n+1}) + source_newton = 0
    */
    Vector<Complex_wp> source_newton;
    
    //! material used
    LorentzMaterial index;
    
    //! adimensionalization parameters
    Real_wp z0_adim, E0_adim, t0_adim;
    
    //! parameters of time impulsion
    Real_wp tau_impulse, Tmax, amplitude_impulse;
    
    //! number of points in time, and padding
    int nb_points_time, padding_time;
    
    //! interval for interpolation grid
    Real_wp zmin_display, zmax_display;

    //! theta parameter for time scheme
    Real_wp theta, theta_polar;
    //! time step
    Real_wp dt;
    Real_wp gamma_implicit;
    //for DIRK method
    VectReal_wp vect_dt_gamma_implicit;
    // for pade scheme (to handle complex roots)
    VectComplex_wp vect_dt_gamma_implicit_cplx;
    //! initial and final time
    Real_wp initial_time, final_time;
    
    //! parameters for time-scheme
    Vector<string> param_time_scheme;

    //! interpolation grid for seismogramms
    GridInterpolation<Dimension1> grid_sismo;
    //! files where seismogramms are written
    string file_output_sismo;
    //! nominal number of points on each seismogramm
    int nb_points_sismo;
    //! time step for seismogramms
    Real_wp dt_sismo;
    //! threshold for seismogramms
    Real_wp threshold_sismo;
    //! initial time for each seismogramm
    VectReal_wp time_begin_sismo;
    //! last sismo point written
    IVect last_sismo_point_number;
    IVect previous_proc_sismo;
    
    //! checking if the mesh should be translated ?
    int nb_iterations_check_mesh;
    //! threshold to know if the mesh should be translated
    Real_wp threshold_translate;
    //! velocity of signal and beginning of translation
    Real_wp velocity_translate, t0_begin_translate;
    Real_wp dz_translation; int nb_iterations_translation;

    //! times for which E will be computed
    VectReal_wp time_interval;
    
    //! central pulsation 
    Real_wp omega_L;
    
    //! pulsations for which E will be computed
    VectReal_wp omega;
    
    //! k(omega) for each pulsation omega
    VectComplex_wp kwave;

    //! FFT in time
    FftInterface<Complex_wp> fft_envelope;
    FftRealInterface fft_time;
    bool harmonic_resolution;
    
    //! type of polarization used
    enum {LINEAR, CIRCULAR, ELLIPTIC};
    int type_polarization;

    //! polarization and phase of initial condition
    TinyVector<Real_wp, 2> polarization_init, phase_init;
    
    //! verbosity level
    int print_level;
    
    //! type of spectrum to compute
    int type_spectrum;
    
    //! different choices of spectrum
    /*!
      CENTERED : from 0 until omega_max
      SHIFTED_SINGLE : from omega_L - omega_max until omega_L + omega_max)
      ODD_FREQUENCIES : enveloppes around odd frequencies omega_L, 3 omega_L, 5 omega_L, etc
    */
    enum {CENTERED, SHIFTED_SINGLE, ODD_FREQUENCIES};
    
    //! number of odd frequencies
    int number_odd_frequencies;
    
    //! if true the Newton's algorithm is inialized with the solution of 
    //! non-linear Schrodinger equation
    bool initialisation_with_schrodinger;
    
    //! coefficient in order to compute the true Fourier transform
    Real_wp coef_fft;
    
    //! pulsation for which Fourier transform should be computed
    VectReal_wp omega_fourier;    
    
  public :
    NonLinearOpticsProblem();
    
    Real_wp GetAdimTime() const;
    void SetFrequency(const Real_wp& f);
    
    void SetInputData(const string& description_field, const VectString& parameters);

    void SetHarmonicResolution();
    
    void ConstructSpectrum(const string&);
    
    Complex_wp GetPulseHarmonic(const Real_wp& omega);
    Complex_wp GetPulseTime(const Real_wp& t);
    void GetPulseTime(const Real_wp& t, Complex_wp&);
    void GetPulseTime(const Real_wp& t, Real_wp&);
    
    void GetElectricFieldFromDisplacement(const Real_wp& t, const VectComplex_wp& D, VectComplex_wp& E);
    void GetDeriveElectricFieldFromDisplacement(const Real_wp& t, const VectComplex_wp& D,
						const VectComplex_wp& dD_dt,
						const VectComplex_wp& E, VectComplex_wp& dE_dt);

    void GetElectricFieldFromDisplacement(const Real_wp& t, const VectReal_wp& D, VectReal_wp& E);
    void GetDeriveElectricFieldFromDisplacement(const Real_wp& t, const VectReal_wp& D,
						const VectReal_wp& dD_dt,
						const VectReal_wp& E, VectReal_wp& dE_dt);
    
  };

  
  //! parameters for a fiber
  class OpticalFiber
  {
  protected :
    Real_wp a, n1, n2;
    int type_solver;
    enum {NEWTON_SOLVER, DICHOTOMY_SOLVER};
    
  public :
    OpticalFiber();
    
    void SetCircularFiber(const Real_wp& b, const Real_wp& n1_, const Real_wp& n2_);
    
    Real_wp EvaluateFunction(int m, const Real_wp& k0, const Real_wp& beta,
                             const Real_wp& Q, const Real_wp& constte);
      
    Real_wp GetKomega(int m, const Real_wp& omega, const Real_wp& c0);
    Real_wp GetAeff(int m, const Real_wp& beta, const Real_wp& omega, const Real_wp& c0);
    Real_wp GetNomega(int m, const Real_wp& beta, const Real_wp& omega, const Real_wp& c0);
    
  };
}

#define MONTJOIE_FILE_NON_LINEAR_OPTICS_PROBLEM_HXX
#endif
