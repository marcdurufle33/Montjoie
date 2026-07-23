#ifndef MONTJOIE_FILE_VAR_AXISYM_PROBLEM_HXX

namespace Montjoie
{

#ifdef MONTJOIE_WITH_TWO_DIM
  class VarAxisymProblem;
  
  class VarComputationRCS_Axi : public VarComputationRCS_Base<Dimension2>
  {
  private:
    DistributedProblem<Dimension2>& var_problem;
    VarAxisymProblem& var_axisym;

  public:
    VectR3 incident_direction;

  public :    
    template<class TypeEquation>
    VarComputationRCS_Axi(EllipticProblem<TypeEquation>& var);
    
    void InitIncidentDirection();
    void GetWaveVectorMonostatic(int k, const Real_wp& omega, R3& kwave, R3& polar);
  };
  
  
  //! base class for axisymmetric problems
  class VarAxisymProblem
  {
  protected:
    IVect Dof_On_Axe;

    //! 3-D wave vector
    R3 kwave3D;
    R3 origine_phase3D; //!< phase origin in definition of plane wave
    bool number_mode_to_be_computed; //!< the number of modes has to be computed
    Real_wp threshold_mode; //!< threshold used to stop the computation

    //! number of modes involved can depend on the incidence
    IVect maximum_mode_for_each_incidence;
    
    //! evaluation of bessel functions enables faster computation, but induces additional storage
    Vector<VectReal_wp> Jn_precomputed;
    
    //! wave number of the envelope
    Complex_wp kwave_envelope;    
    
    Vector<bool> Vertex_On_Axe, Element_On_Axe;
    bool force_diagonal_mass_invertible;
    
  private:
    DistributedProblem<Dimension2>& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    VarSourceProblem_Base& var_source;
    VarGeneralizedImpedance_Base& var_gibc_base;
    
    void InitDefaultValues();

  public:
    template<class TypeEquation>  
    VarAxisymProblem(EllipticProblem<TypeEquation>& var);
    
    // Inline methods
    void SetWaveVector(const R3& k);
    void UpdateWaveVector(const Real_wp& omega);
    void SetPhaseOrigin(const R3& pt); 
    const R3& GetWaveVector() const;
    const R3& GetPhaseOrigin() const;
    const Complex_wp& GetWaveNumberInfinity() const;
    virtual int GetNbModesSource() const;
    
    void GetFourierMode(const Real_wp& teta, Complex_wp& val) const;
    void GetFourierMode(const Real_wp& teta, Real_wp& val) const;

    bool IsVertexOnAxis(int n) const;
    bool IsElementNearAxis(int iquad) const;
    bool NumberOfModesToBeComputed() const;
    Real_wp GetModeThreshold() const;

    Real_wp GetBessel_Value(int n, int num_point) const;
    
    // other methods
    void GetMemoryUsed(map<string, size_t>& var) const;
    void SetInputData(const string& description_field, const VectString& parameters);

    void CheckSectionMeshAxi();
    
    void ComputeDofOnAxe(const VarProblem<Dimension2>& var_pb);

    static void Get_KwavePerp_Kz_Phase(const Complex_wp& rho_tilde, const Complex_wp& mu, int m,
				       const R3& kwave, const Real_wp& omega,
				       Real_wp& k_perp, Real_wp& kz, bool& incidence_axial,
				       Complex_wp& phase);

    // methods needed to compute the radar cross section (monostatic/bistatic)
    void ComputeListMode(VarComputationRCS_Axi& rcs_param);
    int ComputeNbModes_Generic(const Real_wp& max_number);
    
    void InitBesselArray(VarComputationRCS_Axi& rcs_param);
    void InitRcs(VarComputationRCS_Axi& rcs_param);
    int GetNbRightHandSide(VarComputationRCS_Axi& rcs_param) const;    
    
  };
#endif
  
}

#define MONTJOIE_FILE_VAR_AXISYM_PROBLEM_HXX
#endif
