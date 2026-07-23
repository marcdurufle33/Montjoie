#ifndef MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_3D_HXX

namespace Montjoie
{
  
  //! solves Schrodinger equation with FFTs
  /*!
    We consider 3-D non-linear schrodinger equation :
    UPPE : 2 i k_z(omega) dE/dz
          + gammaNL / (epsilon0 c^2) omega^2 F( |Et|^2 E) / exp(i k_z(omega) z)
    
    with Et = F^{-1}(E exp(i k(omega) z)) for the UPPE model
    (Et is the electrical field in time domain)
    whereas E is the envelope of the eletrical in spectral domain
    
    initial condition is equal to E(z=0, t) = amplitude exp(-i omega_L t) exp(-1/2 (t/tau)^2) + c.c.
   */
  class SchrodingerProblem3D : public InputDataProblem_Base, public VirtualOdeSystem<Complex_wp>
  {
  public :
    //! object to compute FFTs in 2-D (in space)
    FftInterface<Complex_wp> fft2d;
    //! object used to compute FFTs in time
    FftInterface<Complex_wp> fft;
    FftRealInterface fft_time;
    //! number of points in fft and padding
    int nb_points_time, padding_time;
    //! force the number of points to be close to a power of 2, 3, 5 and/or 7 ?
    int type_forcing_number_fft_points;
    enum {POWER_NONE, POWER_TWO, POWER_TWO_THREE, POWER_TWO_THREE_FIVE,
          POWER_TWO_THREE_FIVE_SEVEN};

    //! non-linear coefficient
    Real_wp gammaNL;

    //! pulsations omega for which E will be computed
    VectReal_wp omega;
    //! k(omega) for each pulsation omega
    VectComplex_wp kwave;

    //! central pulsation
    Real_wp omega_L;
    //! k(omega_L)
    Complex_wp k_L;

    //! Lorentz material considered
    LorentzMaterial index;
  
    //! parameters of impulsion in time
    Real_wp tau_impulse, Tmax, amplitude_impulse;
    //! times for which E will be computed
    VectReal_wp time_interval;
    //! space step used
    Real_wp dz;
    //! solution is computed from z=0 until z=zmax
    Real_wp zmax;

    //! number of points on each direction
    int nb_points_x, nb_points_y;
    //! space step for each coordinate
    Real_wp dx, dy;
    //! subdivisions on each direction
    VectReal_wp step_x, step_y;
      
    //! computational domain : [xmin, xmax] x [ymin, ymax]
    Real_wp xmin, xmax, ymin, ymax;
    //! waist of the gaussian beam
    Real_wp waist;

    //! adimensionalization parameters
    Real_wp z0_adim, E0_adim, t0_adim;

    //! solution can be written from z0 until z1 at regular intervals
    ParamOutputClass output_grid_z;
    //! directory where output files are written
    string DOSSIER;
    //! interpolation grid on each direction
    GridInterpolation<Dimension1> grid_x, grid_y;
    //! coefficient in order to compute the true Fourier transform
    Real_wp coef_fft;

    //! type of model used
    enum {UPPE_PARAXIAL, UPPE};
    int type_model;

    //! spatial frequencies in spectral domain (kx, ky)
    VectReal_wp nu_x, nu_y;
    //! -Laplacian coefficient (kx^2 + ky^2)
    VectReal_wp coef_Gamma2;
    //! Kz coefficients : sqrt(k(omega)^2 - kx^2 - ky^2)
    VectComplex_wp coef_Kz, coef_InvKz;
    //! coefficient used for split-scheme
    VectComplex_wp coef_Laplacian;
  
    //! type of polarization used
    enum {LINEAR, CIRCULAR, ELLIPTIC};
    int type_polarization;
    
    //! polarization components
    TinyVector<Complex_wp, 2> polarization;
    
    //! use a local enveloppe ?
    bool local_enveloppe;
    Real_wp z0_current, threshold_phase_imag;
    VectReal_wp arg_phaseZ_stored;
    Vector<VectComplex_wp> phaseStored;
    Vector<VectReal_wp> phaseRealStored, invPhaseRealStored;

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
    
    //! if true, the envelope is written in snapshot files
    bool affiche_enveloppe;
    
    // variables for time refinement
    FftInterface<Complex_wp> fft_refinement;
    FftRealInterface fft_refinement_time; Real_wp coef_fft_refinement;
    int factor_refinement_time;
    Real_wp threshold_time_extraction;
    int size_refinement_output;
    FftInterface<Complex_wp> fft_envelope_time; Real_wp coef_fft_envelope_time;

    // variables to check viability of the solution
    bool check_overflow_scheme;
    Real_wp threshold_negative_epsilon, threshold_last_frequency;
    Real_wp epsilon_min_treated, epsilon_max_treated;
    Vector<bool> kwave_to_consider;
    Real_wp coef_stability_damping;
    
    //! how Ez is solved ?
    int type_approximation_Ez;
    enum {EZ_NULL, EZ_DIVERGENCE, EZ_UNKNOWN};

  public:    
    SchrodingerProblem3D();

    Real_wp GetDz() const;
    Real_wp GetZmax() const;

    int GetNbPointsGrid() const;
    int GetNbPointsTime() const;
    
    void SetInputData(const string& keyword, const Vector<string>& parameters);
  
    void SetFrequency(const Real_wp& f);
    
    void SetRegularGrid(const Real_wp& x0, const Real_wp& xN, const Real_wp& y0,
			const Real_wp& yN, int nbx, int nby);
    
    Complex_wp GetPhaseKz(const Complex_wp& k, const Real_wp& z);
    void GetPhaseImagKz(const Complex_wp& k, const Real_wp& z,
			Real_wp& phase_real, Real_wp& inv_phase_real);
    
    void GetLocalPhase(const Vector<Complex_wp>& kw, const Real_wp& z,
		       Vector<Complex_wp>& phase, VectReal_wp& phase_real,
		       VectReal_wp& inv_phase_real);
    
    void ConstructAll();
    
    void InitGrid(const Real_wp& x0, const Real_wp& xN, int Nd, int nb_points,
                  GridInterpolation<Dimension1>& grid);
    
    void InitGrid(const Real_wp& x0, const Real_wp& xN,
                  const Real_wp& y0, const Real_wp& yN, int Nx, int Ny);
    
    void GetInterpolateUn(VectComplex_wp& Un, const VectComplex_wp& xsol) const;
    void GetInterpolateFreqUn(Matrix<Complex_wp>& Un, const VectComplex_wp& xsol) const;
    
    void EvaluateFunction(const Real_wp& z, const Vector<Complex_wp>& En,
                          Vector<Complex_wp>& ProdEn,
                          bool invert_mass = true, bool source = true);

    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectComplex_wp& X,
				    VectComplex_wp& Y, bool invert_mass = true, bool source = true);
    
    void SetDirichletCondition(const Real_wp& t, int nb_deriv,
                               VectComplex_wp& x, Real_wp alpha = 1);
    
    void GiveIterate(int n, const Real_wp& t, VectComplex_wp& En);
    void ModifyIterate(int n, const Real_wp& t, VectComplex_wp& En);
    
    Real_wp GetFrequency() const;
    Real_wp GetFinalTimeSource() const;
    
    // for split scheme
    void SolveSplitOperator(const Real_wp& z, const Real_wp& deltaz,
			    const VectComplex_wp& En, VectComplex_wp& En_next, int op_num);
 
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

#define MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_3D_HXX
#endif
