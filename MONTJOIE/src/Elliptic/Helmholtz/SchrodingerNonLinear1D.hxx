#ifndef MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_1D_HXX

namespace Montjoie
{
  
  //! main class for solving 1-D non-linear schrodinger equation (with kerr effect)
  /*!
    We consider 1-D non-linear schrodinger equation :
    UPPE: 2 i k(omega) dE/dz + gammaNL/(epsilon0 c^2) omega^2 F( |Et|^2 Et ) exp(-i k(omega) z) = 0
    FOP: dE/dz - i k(omega_L)*(omega/omega_L-1) E 
          - i omega k_L/(2 omega_L)* gammaNL/(epsilon0) F( |Et|^2 Et ) exp(-i k_L z) = 0 
    //changer
    NEE: dE/dz - i (k(omega)-k_L-1/vg*(omega-omega_L))E 
    - i gammaNL/(2 epsilon0 c^2) omega_L^2/k_L (1+(omega-omega_L)/omega_L)
    F( |Et|^2 Et ) exp(-i k_L z) = 0 
    NEE constant(NEEC):dE/dz - i (k(omega)-k_L-1/vg*(omega-omega_L))E
    - i gammaNL/(2 epsilon0 c^2) omega_L^2/k_L  F( |Et|^2 Et ) exp(-i k_L z) = 0 
    NLS: dE/dz + i sigma (omega-omega_L)^2 E - i gammaNL/(2 epsilon0 c^2)
    omega_L^2/k_L  F( |Et|^2 Et ) exp(-i k_L z) = 0
    MIRO:2 i k_L dE/dz - 2 k_L sigma omega^2 E
    + gammaNL/(epsilon0 c^2) omega^2 F( |Et|^2 Et ) exp(-i k_L z) = 0
    
    with Et = F^{-1}(E exp(i k(omega) z)) for the UPPE model
    or Et= F^{-1}(E exp(i k(omega_L) z)) for the others
    (Et is the electrical field in time domain)
    whereas E is the enveloppe of electrical field in frequential domain
    sigma is the group dispersion and vg the group velocity
    
    k(omega = n(omega) omega / c
    n(omega) is given by Sellmeier's relations :
    n(omega)^2 = eps_inf + i c/ omega alpha_damping + \sum_k alpha_k
    / ( 1 - i omega/(omega_k^2 T_k) - (omega/omega_k)^2 )    
    
    initial condition is equal to E(z=0, t) = amplitude exp(-i omega_L t) exp(-1/2 (t/tau)^2) + c.c
    Fourier transform is equal to E(z=0, omega)
    = amplitude tau exp(-1/2 tau^2 (omega - omega_L)^2 )
  */
  class KerrProblem : public VirtualOdeSystem<Complex_wp>, public InputDataProblem_Base
  {
  public :
    //! object used to compute FFTs in time
    FftInterface<Complex_wp> fft;
    FftRealInterface fft_time;
    //! number of points in fft and padding
    int nb_points_time, padding_time;
    //! force the number of points to be close to a power of 2, 3, 5 and/or 7 ?
    int type_forcing_number_fft_points;
    enum {POWER_NONE, POWER_TWO, POWER_TWO_THREE,
          POWER_TWO_THREE_FIVE, POWER_TWO_THREE_FIVE_SEVEN};
    
    //! non-linear coefficient
    Real_wp gammaNL;

    //! object for Raman effect
    RamanEffect var_raman;
    VectComplex_wp raman_fct_hat;
    Real_wp coef_convolution;
    
    //! pulsations for which E will be computed
    VectReal_wp omega;
    //! k(omega) for each pulsation omega
    VectComplex_wp kwave;
    //! central pulsation 
    Real_wp omega_L;
    //! k(omega_L)
    Complex_wp k_L;
    Complex_wp group_velocity;
    
    //! parameters of impulsion in time
    Real_wp tau_impulse, Tmax, amplitude_impulse;
    //! times for which E will be computed
    VectReal_wp time_interval;
    //! space step used
    Real_wp dz;
    //! solution is computed from z=0 until z=zmax
    Real_wp zmax;
    //! Lorentz material used (alpha, omega, sigma, eps_inf)
    LorentzMaterial index;
    //! for a monomode fiber
    OpticalFiber fiber;
    bool fiber_computation;
    
    //! adimensionalization parameters
    Real_wp z0_adim, E0_adim, t0_adim;
    
    //! solution can be written from z0 until z1 at regular intervals
    ParamOutputClass output_grid_z;
    //! directory where output files are written
    string DOSSIER;
    //! coefficient in order to compute the true Fourier transform
    Real_wp coef_fft, coef_fft_time;
    
    //! type of model used
    enum {SCHRODINGER_CONSTANT, UPPE, FOP, NEE, NEEC, NLS, MIRO};
    int type_model;
    
    //! type of polarization used
    enum {LINEAR, CIRCULAR, ELLIPTIC};
    int type_polarization;
    
    //! use a local enveloppe ?
    bool local_enveloppe;
    Real_wp z0_current, threshold_phase_imag;
    VectReal_wp arg_phaseZ_stored;
    Vector<VectComplex_wp> phaseStored;
    Vector<VectReal_wp> phaseRealStored, invPhaseRealStored;

    //! polarization components
    TinyVector<Complex_wp, 2> polarization;
    
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

    //! time scheme to use
    All_TimeScheme<Complex_wp> time_scheme;
    
    //! verbosity level
    int print_level;
    
    // variables used in EvaluateFunction
    VectReal_wp omegaAux, n2omega, Aeff_aux, Aeff_om;
    VectComplex_wp kwaveAux; VectComplex_wp scaleAux;
    VectComplex_wp coefAlphaEpsilon;
    //! if true, the envelope is displayed instead of the real solution
    bool affiche_enveloppe;
    
    // variables for time refinement
    FftInterface<Complex_wp> fft_refinement;
    FftRealInterface fft_refinement_time; Real_wp coef_fft_refinement;
    int factor_refinement_time;
    Real_wp threshold_time_extraction;
    int size_refinement_output;
    FftInterface<Complex_wp> fft_envelope_time; Real_wp coef_fft_envelope_time;
    
    // variables to check viability of the solution
    int nb_iterations_between_two_tests;
    bool check_overflow_scheme;
    Real_wp threshold_negative_epsilon, threshold_last_frequency;
    Real_wp epsilon_min_treated, epsilon_max_treated;
    Vector<bool> kwave_to_consider;
    Real_wp coef_stability_damping;
    
    // default constructor
    KerrProblem();
  
    int GetNbPointsTime() const;
    
    void SetInputData(const string& keyword, const Vector<string>& parameters);
    
    void SetFrequency(const Real_wp& f);

    Complex_wp GetPhaseKz(const Complex_wp& k, const Real_wp& z);
    void GetPhaseImagKz(const Complex_wp& k, const Real_wp& z,
			Real_wp& phase_real, Real_wp& inv_phase_real);
    
    void GetLocalPhase(const Vector<Complex_wp>& kw, const Real_wp& z,
		       Vector<Complex_wp>& phase, VectReal_wp& phase_real,
		       VectReal_wp& inv_phase_real);
    
    void ConstructAll();
    
    void Init(MaxwellProblem1D&);

    template<class T>
    void ApplyRamanEffect(const VectReal_wp& moduleU2,
			  Vector<T>& utime, Real_wp c=1);
    
    template<class T>
    void ApplyRamanEffect(const VectReal_wp& moduleU2,
			  Vector<T>& utime, Vector<T>&, Real_wp c=1);
    
    void EvaluateFunction(const Real_wp& zn, const VectComplex_wp& u,
			  VectComplex_wp& scheme, bool invert_mass = true, bool add_source = true);

    Real_wp GetEnergy(VectComplex_wp& En) const;
    
    void GiveIterate(int n, const Real_wp& t, VectComplex_wp& En);
    void ModifyIterate(int n, const Real_wp& t, VectComplex_wp& En);
    
    // fonctions for other time schemes
    void SetDirichletCondition(const Real_wp& t, int nb_deriv,
                               VectComplex_wp& x, Real_wp alpha = 1);
  
    void SourceOnlyTime(const Real_wp& t0, int nb_deriv, Complex_wp& pulse);
    
    Real_wp GetFrequency() const;
    
    Real_wp GetFinalTimeSource() const;

    void SolveSplitOperator(const Real_wp& t, const Real_wp& dt,
                            const VectComplex_wp& Un, VectComplex_wp& Un_next, int num);
    
    // initial condition
    Complex_wp GetHarmonicPulse(const Real_wp&);
    
    void GetInitialCondition(Vector<Complex_wp>& E0);
    
    template<class T>
    void ExtractNonZeroComponents(Real_wp& t0, Real_wp& t1, Vector<T>& utime,
                                  const Real_wp& threshold);

    bool NonLinearPhase() const;
    
    void WriteOutputFile(const string& file_name, const string& file_name_time,
			 const Real_wp& zn, Vector<Complex_wp>& En);

    void WriteOutputFileScalar(const string& file_name, const string& file_name_time,
			       const Real_wp& zn, Vector<Complex_wp>& En);
    
  };

}

#define MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_1D_HXX
#endif
