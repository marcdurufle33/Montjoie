#ifndef MONTJOIE_FILE_HARMONIC_MAXWELL_2D_HXX

namespace Montjoie
{

  //! class to specify the resolution of 2-D time-harmonic Maxwell equations
  /*!
    The time-harmonic maxwell equation solved is
    \f$ -\omega^2 \epsilon \, E + \mbox{curl}( \frac{1}{\mu} \mbox{curl} E ) \, = \, f \f$
  */
  class HarmonicMaxwellEquation_2D : public GenericEquation<Complex_wp>
  {
  public :
    typedef Dimension2 Dimension;
    
    enum {nb_unknowns = 1, nb_unknowns_scal = 1, nb_unknowns_hdg=1,
	  nb_components_en = 1, nb_components_hn = 1, nb_unknowns_vec = 1,
    type_element = 2};
    
    static inline bool ComputeDFjm1() { return false; }
    static inline bool SymmetricGlobalMatrix() { return true; }
    static inline bool SymmetricElementaryMatrix() { return false; }

    template<class TypeEquation>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
				  int num_elem, const ElementReference_Dim<Dimension2>& Fb);
    
    template<class TypeEquation>
    static void ComputeMassMatrixQuad(EllipticProblem<TypeEquation>& var,
				      int num_elem, const QuadrangleHcurlFirstFamily& Fb);
    
    template<class TypeEquation, class T0>
    static void GetPenalizationCoef(const EllipticProblem<TypeEquation>& var,
				    int i, const GlobalGenericMatrix<T0>& nat_mat, int ref, T0& coef);
    
  };
  
  
  //! class to solve 2-D time-harmonic Maxwell equations
  /*!
    The time-harmonic maxwell equation solved is
    \f$ -\omega^2 \epsilon \, E + \mbox{curl}( \frac{1}{\mu} \mbox{curl} E ) \, = \, f \f$
    Edge finite elements are used
  */
  template<class Complexe>
  class HarmonicMaxwell_2D
  {
  public :
    typedef Vector<Complexe> VectComplexe; //!< vector of floats or complex floats
    
    //! magnetic permeability
    Vector<ScalarPhysicalIndice<Dimension2, Complexe> > ref_mu;
    Vector<ScalarPhysicalIndice<Dimension2, Complexe> > ref_invMu;
    //!< inverse of magnetic permeability
    //! dieletric permittivity
    Vector<TensorPhysicalIndice<Dimension2, 2, Complexe> > ref_epsilon;
    Vector<TensorPhysicalIndice<Dimension2, 2, Complexe> > ref_invEpsilon;
    //!< inverse of dieletric permittivity
    Vector<TensorPhysicalIndice<Dimension2, 2, Complexe> > ref_sigma;
    //! magnetic permeability in empty media (used for definition of the plane wave)
    Complexe mu0;
    Complexe invMu0; //!< inverse of magnetic permeability in empty media
    //! dieletric permittivity in empty media
    Complexe epsilon0;
    //!< inverse of dieletric permittivity in empty media
    //! geometric coefficients, J_i DF_i^{-1} DF_i^{*-1}
    Vector<Vector<TinyMatrix<Complexe, Symmetric, 2, 2> > > Glob_matMass_Bh;
    //! geometric coefficients 1/J_i
    Vector<Vector<Complexe> > Glob_matMass_Dh;
    //! impedances for first order absorbing condition for each physical media
    Vector<Complexe> coefficient_impedance_absorbing;
    
  private:
    VarProblem<Dimension2>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension2>& var_boundary;

  public :    
    template<class TypeEquation>
    HarmonicMaxwell_2D(EllipticProblem<TypeEquation>&);

    void CopyInputData(const HarmonicMaxwell_2D<Complexe>& var_base);
    
    template<class T0, class Prop>
    void GetMassMatrixCoef(int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int m, int n,
			   TinyMatrix<T0, Prop, 2, 2>& mass) const;
    
    template<class T0>
    void GetStiffMatrixCoef(int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			    int m, int n, T0& stiff) const;
    
    //! allocation of indices arrays
    void InitIndices(int n);
    int GetNbPhysicalIndices() const;
    //! filling of indices arrays depending the data file
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complexe>* >& rho_complex,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void FinalizeComputationVaryingIndices();
    bool IsVaryingMedia(int i) const;
    
    void ModifyPhysicalCoefPML(TinyMatrix<Complexe, Symmetric, 2, 2>& eps,
			       Complexe& mu, TinyMatrix<Complexe, Symmetric, 2, 2>& sigma, int, int) const;
    
    void AllocateMassMatrices();
    
  };


  //! class for first-order absorbing boundary condition
  template<>
  class ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_2D>
    : public ImpedanceFunction_Base<Complex_wp, Dimension2>
  {
  private:
    const VarProblem<Dimension2>& var_problem;
    const HarmonicMaxwell_2D<Complex_wp>& var_maxwell;

  public :
    template<class TypeEquation>
    ImpedanceABC(const EllipticProblem<TypeEquation>& var);
    
    //! impedance
    void EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                              const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                              const SetPoints<Dimension2>& Pts, const SetMatrices<Dimension2>& Mat);  
    
  };
  

  //! base class for solving 2-D Maxwell's equations
  template<class TypeEquation>
  class HarmonicMaxwell_2D_Eq : public HarmonicMaxwell_2D<typename TypeEquation::Complexe>,
				public VarHarmonic<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;

  public:
    HarmonicMaxwell_2D_Eq();

    void InitIndices(int n);
    int GetNbPhysicalIndices() const;

    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    Real_wp GetVelocityOfMedia(int ref) const;
    Real_wp GetVelocityOfInfinity() const;
    
    void FinalizeComputationVaryingIndices();
    bool IsVaryingMedia(int i) const;
    bool IsVaryingMedia(int m, int i) const;
    
    void AllocateMassMatrices();
    
  };


  //! base class for solving time-harmonic 2-D Maxwell's equations
  template<class TypeEquation>
  class HarmonicMaxwell_2D_Cplx : public HarmonicMaxwell_2D_Eq<TypeEquation>
  {
    typedef typename TypeEquation::Complexe Complexe;
    
  public:
    void GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
			   Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
			   IVect& num_ref, IVect& num_index, IVect& num_component,
			   Vector<bool>& compute_grad, Vector<bool>& compute_hess);
    
  };
  

  //! class to solve 2-D time-harmonic Maxwell equations
  template<>
  class EllipticProblem<HarmonicMaxwellEquation_2D>
    : public HarmonicMaxwell_2D_Cplx<HarmonicMaxwellEquation_2D>
  {
  public:
    IncidentWaveProjector<Complex_wp, Dimension2>*
    GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) const;

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);
    
  };

  template<>
  class FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_2D>
    : public FemMatrixFreeClass_Eq<Complex_wp, HarmonicMaxwellEquation_2D>
  {
  public:
    
    FemMatrixFreeClass(const EllipticProblem<HarmonicMaxwellEquation_2D>& var_);

    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Real_wp>& X, Vector<Real_wp>& Y) const;
    
    virtual void MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
			    const SeldonTranspose&, int lvl, 
			    const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const;
    
  };
  
} // end namespace

#define MONTJOIE_FILE_HARMONIC_MAXWELL_2D_HXX
#endif
