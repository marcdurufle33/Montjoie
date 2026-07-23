#ifndef MONTJOIE_FILE_VLASOV_MAXWELL_HXX

namespace Montjoie
{
  
  //! class to specifiy a beam of particles or space charge limited emission
  template<class Dimension>
  class Vlasov_BeamCurrent
  {
    public :
    typedef typename Dimension::R_N R_N;
    
    //! extremities of the beam
    R_N ptA, ptB, ptC, ptD;
    //! center of the beam
    R_N center;
    //! velocity of the particles of the beam
    R_N v0;
    //! current of the beam
    Real_wp I0;
    //! charge of each particle
    Real_wp charge;
    //! ratio q/m for each particle
    Real_wp ratio_q0_over_m0;
    //! coefficient used in the case of gaussian beam (exp^(-alpha r^2) )
    Real_wp alpha;
    //! breakdown field for space charge limited emission
    Real_wp breakdown_field;
    //! reference of the boundary where the creation of particles is allowed
    int ref;
    //! number of particles to create for each meter
    Real_wp nb_part_per_unit;
    //! weights
    VectReal_wp WeightBeam;
    //! type of beam (uniform or gaussian)
    int type_beam;
    // possible choices
    enum {UNIFORM_BEAM, GAUSSIAN_BEAM}; //!< gaussian beam
    //! number of iterations between each creation of particles 
    int delta_iter;
    //! particles are created uniformly in space ? or randomly ?
    bool uniform_spatial;
    //! for space charge limited emission, gamv_normal_min < p . n < gamv_normal_max
    Real_wp gamv_normal_min, gamv_normal_max;
    //! for space charge limited emission, gamv_tangential_min <= p . t <= gamv_tangential_max
    Real_wp gamv_tangential_min, gamv_tangential_max;
    //! distance of the particles from the boundary ,  hmin <= h <= hmax
    Real_wp height_min, height_max;
    
  public:
    Vlasov_BeamCurrent();
    
    //! reading of a line of the data file BeamCurrent = UNIFORM ....
    void SetInputData(const string& description_field, const VectString& parameters);
    //! reading of a line of the data file BreakdownCurrent = UNIFORM ....
    void SetInputBreakdown(const string& description_field, const VectString& parameters);
    
    //! select the appropriate mass and charge depending the input parameters
    static void SetInputSpecies(Real_wp& charge, Real_wp& mass, const VectString& param, int& nb);
    
    //! return 1 for an uniform beam, or a gaussian for a gaussian beam
    Real_wp GetFunction(const Real_wp& lambda);
    
    //! same for 3-D case
    Real_wp GetFunction(const Real_wp& lambda, const Real_wp& mu);
    
  };
  
  //! class to specify source for the current J
  template<class Dimension>
  class CurrentSource_Vlasov
  {
  prublic :
    //!< stationary problem
    typedef typename Dimension::R_N R_N; //!< point in R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    
    //! evaluation of J
    VectR_N evalJ;
    //! evaluation of rho
    Vector<Real_wp> evalRho;
    //! coefficients \f$ \hat{S}(|x_k - \xi_m|) \f$
    VectReal_wp CoefficientShape, coef_q0;
    //! coefficients used for the computation of the shape function
    Real_wp alpha, beta;
    //! variable R (i.e. influence radius of the particles)
    Real_wp cut_off_Sshape;
    //! 1/R 
    Real_wp invRadius;
    //! 1/R^2
    Real_wp invR2;

  private:
    VlasovMaxwell<Dimension>& var_time;
    Mesh<Dimension>& mesh;
    Vector<VectReal_wp>& Glob_jacobian;
    Real_wp xmin, xmax, ymin, ymax, zmin, zmax;
    
    void InitDefaultValues();
    
  public:
    template<class TypeEquation>
    CurrentSource_Vlasov(HyperbolicProblem<TypeEquation>& var);
    
    void UpdateCurrent(const Real_wp& poids_J, const R_N& vk, const R_N& xk,
		       int num_elem, const R_N& pt_loc);
    
    Real_wp GetShapeFunction(const Real_wp& dist);
    
  };
  
  
  //! overloaded class
  template<class Dimension>
  class VlasovMaxwell
  {
  };

  
  //! abstract class for Vlasov-Maxwell equations
  template<class Dimension>
  class VlasovMaxwell_Base
  {
  public :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    
    //////////////////////////
    // PARTICLES PARAMETERS //
    
    //! increment when adding particles
    int increment_particles;
    //! number of particles for each kind
    IVect nb_particles;
    //! maximum number of particles for each kind
    IVect nb_max_particles;
    //! weight for each particle
    Vector<VectReal_wp> weight_particle;
    //! charge for each kind
    VectReal_wp charge_species;
    //! ratio q/m for each kind
    VectReal_wp ratio_qm_species;
    
    //! localization of particles inside the mesh
    GridInterpolation<Dimension> loc_particle;
    //! initial parameters for the particles
    Vector<VectString> initial_particles_param;
    
    //! name of the file where the tracking of particles have to be written
    string name_file_position;
    //! if true, particles are considered relativist
    bool relativist_particle;
    
    //! eventual beams
    Vector<Vlasov_BeamCurrent<Dimension> > Liste_beam_current;
    //! boundaries with space charge limited emission
    Vector<Vlasov_BeamCurrent<Dimension> > breakdown_currents;
    //! algorithm used to detect the intervals where E > E_{breakdown}
    int type_space_charge_limited_emission;
    //! available algorithms
    enum{FIND_INTERVAL, SUBDIVISION_INTERVAL, FULL_INTERVAL};
    
    //! computation of energy ?
    bool computation_energy;
    //! file where the energy history is stored
    string file_output_energy;
    
    // PARTICLES PARAMETERS //
    //////////////////////////
    
    ///////////////////////////
    // RESOLUTION PARAMETERS //
    
    //! list of quadrature points
    Vector<VectR_N> PointsQuadrature, PointsInterpolation;
    VectReal_wp WeightsQuadrature;
    //! order of integration for J
    int order_integration_current; int nb_points_quadrature_current;
    //! order of interpolation
    int order_interpolation_current; int nb_points_interpolation_current;
    //! jacobian at quadrature points
    Vector<VectReal_wp> Glob_jacobian;
    Vector<Vector<MatrixN_N> > Glob_DFjm1, Glob_DFjm1_Boundary;
    //! normales on each edge/face
    Vector<Vector<R_N> > Glob_normale; Vector<Vector<R_N> > Glob_point_boundary;
    VectReal_wp inv_mass_matrix;
        
    //! type of correction in order to enforce charge conservation law
    int type_correction_divergence;
    //! parameters \f$ \chi \mbox{ and } \sigma\f$ for the hyperbolic correction
    Real_wp xsi_hyperbolic_correction, sigma_hyperbolic;
    
    // available corrections to enforce charge conservation law
    enum{CORRECTION_NONE, CORRECTION_BORIS, CORRECTION_HYPERBOLIC, CORRECTION_INITIAL};
    
    //! object used for Lagrange multiplier for the divergence correction
    EllipticProblem<LaplaceEquation<Dimension> > var_laplace;
    All_LinearSolver laplace_solver;
    
    //! translation vectors to deal with the periodicity
    /*!
      By example for a full 2-D periodicity, there will
      eight vectors (a,0) (b,0), (-a,0), (-b,0)
      (a,b) (a,-b) (-a,b) and (-a,-b)
     */
    VectR_N periodic_vector;
    
    //! radius of influence for each particle
    Real_wp radius_influence_particle;
    //! parameter alpha for each particle
    Real_wp alpha_influence_particle;
    //! maximum allowed velocity
    Real_wp maximum_velocity;
    //! true if the current J is first interpolated and then integrated
    bool regular_interpolation_current; int nb_subdiv_interpolation; 
    //! operators useful to project to nodal, quadrature or interpolation points
    FiniteElementProjector
    proj_interpolate_to_quadrature, proj_nodal_to_quadrature;
    
    //! localization of neighbouring quadrature points
    Vector<IVect> CloudElem_QuadraturePoint;
    
    //! parameters for writing of rho and div(varepsilon E) - rho
    ParamOutputClass output_rho_param;
    //! parameters for writing the positions and speeds of particles
    ParamOutputClass output_particle_param;
    VarRandomGenerator var_random;
        
    // RESOLUTION PARAMETERS //
    ///////////////////////////

  private:
    DistributedProblem<Dimension>& var_harmonic;
    VarInstationary_Base& var_time;
    VlasovMaxwell<Dimension>& var_vlasov;

  public:
    VlasovMaxwell_Base(HyperbolicProblem<TypeElement, TypeEquation>& var);
    
    void SetInputData(const string& description_field, const VectString& parameters);
    void InitTimeIterations();
          
    // counting the particles
    int GetMaximumNumberParticles() const;
    int GetNumberParticles() const;
    int GetInitialNumberParticles() const;
    void CreateNewSpecies(const Real_wp& q0, const Real_wp& m0, VectReal_wp& Y);
    
    // how to get v from p, and vice versa
    void GetMomentum_FromVelocity(const R_N& speed, R_N& momentum) const;
    void GetVelocity_FromMomentum(const R_N& momentum, R_N& speed) const;
    void GetVelocity_FromMomentum(const R_N& momentum, R_N& speed, Real_wp&, Real_wp&) const;
  };
  
} // end namespace

#define MONTJOIE_FILE_VLASOV_MAXWELL_HXX
#endif
