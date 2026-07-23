#ifndef MONTJOIE_FILE_PHYSICAL_CONSTANT_HXX

namespace Montjoie
{    
  
  class PhysicalConstant
  {
  public :
    // several physical constants (used for the adimensionalization)
    static const Real_wp speed_light; //!< speed of light

    //! 1/c^2 where c is the speed of light
    static const Real_wp invSquare_speed_light; 
    static const Real_wp elementary_charge; //!< elementary charge
    static const Real_wp mass_electron; //!< electron mass
    static const Real_wp mass_proton; //!< proton mass
    static const Real_wp mass_neutron; //!< neutron mass

    //! vacuum dielectric permittivity
    static const Real_wp epsilon0_permittivity;

    //! vacuum magnetic permeability
    static const Real_wp mu0_permeability;
    static const Real_wp impedance0;

    //! do we use adimensionalized variables or physical variables ?
    /*! adim_one corresponds to the case where we consider eps0 = mu0 = 1
                and the data file is supposed to have too 
		adimensionalized inputs
       adim_yes corresponds to the case where we have eps0 = mu0 = 1
                but for the user, these values are real values
                and consequently the data file is supposed
		to have physical inputs
                and the outputs are modified to be physical too
                this case is useful in order to have
		values near of one for the computations
                and the IO are physical so that the physicians are not lost
      adim_no corresponds to the case where eps0 and mu0 are real values
                and data file and outputs are physical too
    */
    static int adimensionalization;

    //! maximal number of physical indices
    static int nb_max_indices;
    
    //! types of adimensionalizations
    enum {ADIM_NO, ADIM_YES, ADIM_ONE};

  };

  void getElement_Equation(string& file_name, string & type_element,
			   string& type_equation);
  
  class LorentzMaterial
  {
  public :
    //! velocity c0
    Real_wp c0;
    //! parameters omega_k, alpha_k (for each linear polarization)
    VectReal_wp omega_polarization, sigma_polarization, alpha_polarization;
    //! relative permittivity for infinite frequency
    Real_wp eps_inf;
    //! damping coefficient
    Real_wp alpha_damping;
    //! if true negative values of n2 are forbiddedn
    bool forbid_negative_index;
    
    LorentzMaterial();

    Complex_wp GetSqrtComplex(const Complex_wp& n2);
    
    Complex_wp GetWaveNumber(const Real_wp & w);
    Complex_wp GetGroupVelocity(const Real_wp & w);
    Complex_wp GetDispersion(const Real_wp & w);
    
    void SetSigma(const Real_wp& sigma, const Real_wp& z0_adim);
    void SetRefractiveIndex(const string& material, const Vector<string>& parameters);
  };


  //! class handling Drude's model
  /*!
    this class implements the following model of permittivity :
    eps(omega) = eps_inf ( 1 - \sum_k (omega_p^2 - i omega sigma) / (omega^2 + i gamma omega - omega_0^2))
    the coefficients eps_inf sigma, eps_inf omega_p^2, omega_0^2, gamma are stored as arrays
    eps_inf is a scalar
   */
  class DrudeMaterial
  {
  public:
    VectReal_wp eps_omega_p2, omega_02, gamma, eps_sigma;
    Real_wp eps_inf;
    int type_mode;
    enum {NONE, TE, TM};
    
    DrudeMaterial();

    void SetInputData(int nb, const Vector<string>& parameters, const string& material);    
    void SetRefractiveIndex(const string& material, const Vector<string>& parameters);
    void FillPadeExpansionEV(VectComplex_wp& Om, VectComplex_wp& Sig,
                             Real_wp gam = Real_wp(0), Real_wp sig = Real_wp(0));
    
    VectComplex_wp GetPoleOmega() const;
    VectComplex_wp GetZeroOmega(Real_wp cte = 0) const;
    
    void Adimensionalize(const Real_wp& L0);
    Complex_wp GetEpsilon(const Real_wp& omega) const;
    Complex_wp GetDeltaEpsilon(const Complex_wp& omega) const;
    
    bool IsEnabled() const;
    inline bool IsModeTE() const { if (type_mode == TE) return true; return false; }
    inline bool IsModeTM() const { if (type_mode == TM) return true; return false; }
    
    void ModifyCoefficientHelmholtz(const Real_wp&, Complex_wp& rho, Complex_wp& sigma, TinyMatrix<Complex_wp, Symmetric, 2, 2>& mu) const;
    void ModifyCoefficientHelmholtz(const Real_wp&, Complex_wp& rho, Complex_wp& sigma, TinyMatrix<Complex_wp, Symmetric, 3, 3>& mu) const;
    void ModifyCoefficientHelmholtz(const Real_wp&, Real_wp& rho, Real_wp& sigma, TinyMatrix<Real_wp, Symmetric, 2, 2>& mu) const;
    void ModifyCoefficientHelmholtz(const Real_wp&, Real_wp& rho, Real_wp& sigma, TinyMatrix<Real_wp, Symmetric, 3, 3>& mu) const;

    void ModifyCoefficientMaxwell(const Real_wp&, TinyMatrix<Complex_wp, Symmetric, 3, 3>& rho, TinyMatrix<Complex_wp, Symmetric, 3, 3>& sigma,
                                  TinyMatrix<Complex_wp, Symmetric, 3, 3>& mu) const;

    void ModifyCoefficientMaxwell(const Real_wp&, TinyMatrix<Real_wp, Symmetric, 3, 3>& rho, TinyMatrix<Real_wp, Symmetric, 3, 3>& sigma,
                                  TinyMatrix<Real_wp, Symmetric, 3, 3>& mu) const;
    
  };
  
}

#define MONTJOIE_FILE_PHYSICAL_CONSTANT_HXX
#endif

