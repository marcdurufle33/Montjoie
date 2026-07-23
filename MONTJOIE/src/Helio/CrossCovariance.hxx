#ifndef MONTJOIE_FILE_CROSS_COVARIANCE_HXX

namespace Montjoie
{

  class MySource : public VirtualSourceField<Complex_wp, Dimension3>
  {
  private :
    int surfRef;
    R3 pt_source;
    int m; Real_wp omega;
    int Lmax;
    Vector<Real_wp> diracAmp; 
    Vector<VectReal_wp> Psrc; // values of P_l(cos(thetaSrc))

  public :
    mutable AssociatedLegendrePolynomial<Real_wp> assocLegendre;

    MySource(const R3& pt, const Real_wp& om, int m0, const Vector<Real_wp> diracAmp);

    void EvaluateFunction(const R3& x, Vector<Complex_wp>& f) const; 

  };

  class RandomSource : public VirtualSourceField<Complex_wp, Dimension3>
  {
  private :
    int m; Real_wp omega;
    int Lmax;
    Vector<VectReal_wp> Psrc; // values of P_l(cos(thetaSrc))
    string source_power_type;
    Real_wp Ts;  // source correlation time
    Real_wp source_power_freq;
    Vector<Real_wp> depth_sources;
    Vector<Real_wp> delta_depth_sources;
    Vector<Real_wp> source_power_depth;

  public :
    mutable AssociatedLegendrePolynomial<Real_wp> assocLegendre;
    Vector<Vector<Complex_wp> > realisations; // random realisations as a function of depth_sources and l

    RandomSource(const Real_wp& om, int m0, int lmax, string sourcePowerType, Real_wp ts, Vector<Real_wp> depthSources, Vector<Real_wp> deltaDepthSources, Vector<Real_wp> sourcePowerDepth);

    Real_wp GetSourcePower(const Real_wp& omega);

    void EvaluateFunction(const R3& x, Vector<Complex_wp>& f) const; 

  };

  


  class CrossCovarianceParameter : public MultiFrequencyParameter
  {
  public :
    VectR3 pt_source;
    Vector<string> srcLines; 
    int    nSrc, lMax;
    // for the filtered source
    bool filteredSrc;
    string typeSource;
    Vector<Vector<Real_wp> > diracAmps; // amplitude of each l coefficients for each source
    int ref_surface_observation; // reference of the  observation height

    EllipticProblem<HelmholtzEquationAxi>* var;
    MultiFrequencyProblem* var_multi;
    bool display_solution;
    Real_wp omega_stored;
    string damping_type;
    Real_wp gamma0, omega0, beta; // damping parameters
    Real_wp gammaOld; // saved value of the damping at previous frequency

    // source power parameter
    string source_power_type;
    Real_wp Ts;  // source correlation time
    Vector<Real_wp> depth_sources;
    Vector<Real_wp> delta_depth_sources;
    Vector<Real_wp> source_power_depth;

    Real_wp RSUN;
    VectComplex_wp x_solFull; // to store for solution for all m if required

    
#ifdef SELDON_WITH_MPI
    MPI_Comm comm_freq;
#endif
    
  public:
    CrossCovarianceParameter();

    void SetPulsation(int m, int mg, const Real_wp& omega,
                      VarComputationProblem& var_pb, MultiFrequencyProblem& pb);

    Real_wp GetDamping(const Real_wp& omega);

    void SetLmax(int lMax);
  
    int GetRealM(int mg);

    void SetParameters(int k, int kg, int m, int mg,
                       VarComputationProblem& var_pb, MultiFrequencyProblem&);

    void InitSource(Vector<string> lines_data_file);

    void InitDamping(VectString dampingString);

    void CreateFilter(VectString filterString);

    void InitSourcePower(VectString sourcePowerString);
    
    //void SetInputData(const string& keyword, const Vector<string>& params);

    void SetSourceFrequency(int k, int kg, int m, int mg, const VarComputationProblem& var_pb,
                            Matrix<Complex_wp, General, ColMajor>& sol);
  
    void GiveSolutionFrequency(int k, int kg, int m, int mg, const VarComputationProblem& var_exp,
                               Matrix<Complex_wp, General, ColMajor>& sol_exp);
  
  };
  
  class CrossCovarianceProblem : public MultiFrequencyProblem
  {
  public:
    Real_wp RSUN; // solar radius
    VectR3& pt_source;
    CrossCovarianceParameter& all_param;
    
    bool compute_covariance; // forward model computation
    bool compute_kernel; // kernel computation
    Real_wp rmin, rmax; // if the sources are located only between rmin and rmax
    Vector<string> typeOfKernels; // type of kernels to be computed
    string directoryOutput; // where to store the result
   
    // Kernel
    Matrix<Complex_wp> mat_kernel; //  matrix (ID_pt*len(typeOfKernels)+ID_type, ID_pair)
    VectReal_wp coef_freq;
    int num_grid;

    // parameters for refinement
    bool refinement_vertex_source;
    Real_wp ratio_refinement;
    int lvl_refinement;
    
    // matrice de covariance stockee sous forme creuse
    // row_mat, col_mat : numeros de lignes et de colonnes
    // mat_covariance(i, :) : valeurs non nulles pour la frequence i
    Vector<int> row_mat, col_mat;
    Matrix<Complex_wp> mat_covariance;
    
    int nb_modes;
    MontjoieTimer chrono;
    string filename;
    int ref_surface_observation; // reference utilisee pour l'integrale de surface

        
    VectReal_wp mat_mass;
    
  public :
    CrossCovarianceProblem(CrossCovarianceParameter& param);
    
    void SetInputData(const string& keyword, const Vector<string>& params);

    void SetProblem(DistributedProblem<Dimension2>& var);
    
    void InitCovariance(VarOutputProblem_Dim<Dimension2>& var_output);
    
    template<class TypeEquation> 
    void getPhysicalCoefficients(EllipticProblem<TypeEquation>& var, int ID_elt, int ID_pt, int ref, Complex_wp& rho, Complex_wp& c, Complex_wp& gamma);

    template<class TypeEquation>
    void ComputeMassMatrix(EllipticProblem<TypeEquation>& var);
    
    template<class TypeEquation>
    void computeXS(EllipticProblem<TypeEquation>& var, int k, int mg, int ns,
                   Matrix<Complex_wp, General, ColMajor>& all_sol);

    template<class TypeEquation>
    void computeKernel(EllipticProblem<TypeEquation>& var, int k, int mg, int ns,
                       Matrix<Complex_wp, General, ColMajor>& all_sol);
    
    void WriteCovariance();
    void WriteKernel();
    
  };
  
}

#define MONTJOIE_FILE_CROSS_COVARIANCE_HXX
#endif
