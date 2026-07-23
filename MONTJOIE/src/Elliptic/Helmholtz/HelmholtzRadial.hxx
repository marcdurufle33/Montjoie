#ifndef MONTJOIE_FILE_HELMHOLTZ_RADIAL_HXX

namespace Montjoie
{
  class HelmholtzEquationRadial : public HelmholtzEquation1D
  {
  };

  
  //! class for atmosphere condition
  class AtmosphereImpedanceCondition
  {
  protected:
    int order;
    Real_wp alpha;
    bool exact_condition;
    bool abc_condition;
    bool ha_condition;
    Real_wp threshold;

  public:
    AtmosphereImpedanceCondition();
    
    int ReadParameters(int& nb, const VectString& parameters, int side);
    void UpdateImpedance(const Real_wp& rt, const Complex_wp& kinf, Complex_wp& impedance);
    bool LaplacianPresent() const;
    bool ExactCondition() const;
    
    Complex_wp GetLaplacianImpedance(const Real_wp& factor_laplace, const Real_wp& rt, 
				     const Complex_wp& kinf);
    
    void ComputeExactImpedance(int Lmax, const Real_wp& rt, const Complex_wp& kinf,
			       VectComplex_wp& impedance);
    
  };
  
  
  //! base class for solving Helmholtz equation in spherical coordinates
  template<class Complexe>
  class VarHelmholtz_Radial : public VarHelmholtz_1D<Complexe>
  {
  public :
    typedef Vector<Complexe> VectComplexe; //!< "complex" vector
    
    //! maximal mode number (l)
    int Lmax;
    //! threshold used to stop the computation
    Real_wp threshold_mode;
    //! if true the number of involved modes is automatically computed
    bool number_mode_to_be_computed;
    Vector<ParamOutputClass> output_grid_param; //!< output parameters
    //! predefined display grids
    Vector<GridInterpolationFull<Dimension1> > var_grid;
    //! diagonal matrix for \int mu \varphi_j \varphi_i
    Vector<Complexe> mass_matrix;
    //! matrix for -\omega^2 \int \rho r^2 \varphi_j \varphi_i
    //! + \int \mu r^2 d/dr \varphi_j d/dr \varphi_i
    Matrix<Complex_wp, General, BandedCol> stiffness_matrix;
    //! object used to compute spherical harmonics
    AssociatedLegendrePolynomial<Real_wp> fct_leg;
    //! type of right hand side
    int type_source;
    enum{SRC_DIFFRACTED_FIELD, SRC_TOTAL_FIELD, SRC_DIRAC, SRC_GAUSSIAN, SRC_GAUSSIAN_R};
    //! direction of plane wave
    R3 wave_vector;
    //! origin of Dirac
    R3 origin_dirac; Real_wp radius_gaussian_src, radius_gaussian_cut_off;
    //! quadrature rules used to perform an integration in theta-coordinate
    Globatto<Real_wp> gauss_teta;
    //! impedance for transparent condition
    Vector<Complexe> coef_impedance_transparent;
    //! object for atmosphere condition
    AtmosphereImpedanceCondition atmo_condition;
    
    //! scaling parameters
    bool diagonal_scaling; VectReal_wp diagonal_scale_left, diagonal_scale_right;

    bool display_du_dr;
    
    VarHelmholtz_Radial();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    int ReadBoundaryCondition(int& nb, const VectString& parameters, int side);
    void TreatAtmosphericCondition();
    void AddLaplacianAtmosphericCondition(const Real_wp& factor_laplace, VirtualMatrix<Complex_wp>& mat);
    
    void ConstructAll(const string&, const string&, MontjoieTimer&);
    void ComputeVarGrid();
    
    void PerformFactorizationStep(Matrix<Complex_wp, General, BandedCol>& mat_lu,
                                  IVect& pivot, int l);
    
    void ComputeSolution(Matrix<Complex_wp, General, BandedCol>& mat_lu,
                         IVect& pivot, Vector<Complex_wp>& x);
    
    void ComputeIncidentField(const Real_wp& kinf, const Real_wp& kr, Vector<VectComplex_wp>& u_inc,
                              Vector<VectComplex_wp>& du_inc_dr, bool compute_grad);
    
    void ComputeRightHandSide(Vector<Vector<Vector<Complexe> > >& rhs, MontjoieTimer&);
    void ComputeRhsDiffractedField(Vector<Vector<Vector<Complexe> > >& rhs);
    void ComputeRhsTotalField(Vector<Vector<Vector<Complexe> > >& rhs);
    
    void WriteDatas(Vector<Vector<Vector<Complexe> > >& sol);
    
    // Inline methods
    int GetLmax() const;
    
    int GetIndexOutputFiles(int type);

    Complex_wp GetYnm(int l, int m2, const Vector<VectReal_wp>& P,
                      const Complex_wp& phase);
    
  };
  
  
  template<>
  class EllipticProblem<HelmholtzEquationRadial> : public VarHelmholtz_Radial<Complex_wp>
  {
  };
  
}

#define MONTJOIE_FILE_HELMHOLTZ_RADIAL_HXX
#endif
