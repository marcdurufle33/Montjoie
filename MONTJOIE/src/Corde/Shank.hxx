#ifndef MONTJOIE_FILE_SHANK_HXX

namespace Montjoie
{
  
  //! Class for the shank
  template<int nb_base, int nb_quad, class TypeEquation>
  class Shank : public SpaceStringScheme<nb_base, nb_quad, TypeEquation>,
		public InputDataProblem_Base
  {
  public :
    enum {dimension = TypeEquation::dimension};
    
    //! equation solved for the motion of the shank
    TypeEquation var_eq; 
    //! directory where the outputs are written
    string DOSSIER;
    
    Real_wp epsilon; //!< tolerance pour les differences finies

    //! initial value for theta
    Real_wp initial_teta;
    //! center of the impact between the shank and the jack
    Real_wp point_impact_jack;
    //! parameters of the jack (thicknesss and slope used to construct the delta function)
    Real_wp thickness_jack, slope_jack;
    //! maximal value of theta for which the jack is pushing
    Real_wp theta_end_jack, distance_min_jack;
    //! if true the jack is pushing until a given theta
    bool push_jack_until_theta;
    //! vecteur de repartition de Fjack (creux)
    Vector<Real_wp, VectSparse> repartition_sparse_jack;        
    //! position des ddls du maillage
    VectReal_wp PositionDofs;
    // Force du jack lue dans le fichier
    Vector<Real_wp> ForceJackExt, TangentialForceJackExt, TorqueJackExt;
    // Force du jack apres interpolation (en temps)
    Vector<Real_wp> NewForceJack, NewTangentialForceJack, NewTorqueJack;
    //! Fjack est un creneau entre t_begin_jack et t_end_jack ou DeltatforceJack si on utilise un fichier
    Real_wp t_begin_jack, t_end_jack, amplitude_jack, DeltatForceJack, DeltatTangentialForceJack, DeltatTorqueJack;
    //! function in time for jack
    int type_function_jack, type_function_jack_tangential, type_function_jack_torque;
    enum {JACK_RECTANGLE, JACK_COS_SQUARE, JACK_X_SQRT_COS, FROM_FILE};
    // true si la force exercée par le jack a une composante tengentielle
    bool tangential_component, torque_component;
    //! type of condition for stopping jack pressure
    int type_condition_jack;
    enum {JACK_THETA, JACK_DISTANCE};
    
    Real_wp E; //!< module d'Young
    Real_wp rho; //!< densite lineique du manche
    Real_wp A; //!< section (m^2)
    Real_wp I; //!< module d'inertie
    Real_wp G; //!< module de cisaillement
    Real_wp k_prime; //!< param de Timoshenko
    Real_wp T0; //!< parametre nul
    Real_wp L; //!< longueur du manche
    Real_wp gravity; //!< gravite (constante g)
    R3 amo, amoB2; //!< amortissement sur w et phi
    Real_wp amo_theta; //!< amortissement en teta
    //! pour t > time_attrape, on met un amortissement plus grand
    Real_wp time_attrape, amo_theta_after_attrape;
    
    //! temps initial, final, temps courant, et Delta t entre deux affichages
    Real_wp t_begin, t_end, t_courant, t_display;
    Real_wp Deltat; //!< pas de temps
    Real_wp Deltat2; //!< carre du pas de temps
    Real_wp invDeltat; //!< inverse du pas de temps
    Real_wp invDeltat2; //!< 1 / (Delta t)^2
    Real_wp invDeuxDeltat; //!< 1 / (2 Delta t)
    
    Real_wp alpha; //!< theta parameter used in the theta-scheme
    MegaNewtonSolver Newton_solver; //!< non-linear solver used
    bool neglect_w_square; //!< si vrai, les termes en w^2 sont negliges (petites deplacements)
    
    //! parametres de sortie pour calculer l'interpolee du deplacement
    Real_wp t_begin_interp , t_end_interp , deltat_interp;
    //! parametres de sortie d'un sismogramme sur le manche
    Real_wp t_begin_sismo, t_end_sismo, deltat_sismo;

    GridInterpolation<Dimension1> grid_interp_time;// pour eventuelle force exterieure interpolee
    Globatto<Real_wp> base_time; // fonctions de base ( quadrature lobatto )
    Mesh<Dimension1> mesh_time; // maillage, points etc...
    
    //! on sort t, w(L), phi(L) et theta dans le sismo
    WriteOnTheGoWithBuffer<Real_wp> output_shank;
    // on sort F_jack dans un sismo special
    WriteOnTheGoWithTinyBuffer<Real_wp, 2> output_Fjack;
    string file_output_sismo, file_output_jack, file_output_force;
    
    //! u(offset_shank) renvoie a la premiere inconnue du manche
    int offset_shank;
    //! numero utilise pour reperer l'inconnue teta
    int teta_dof_number;
    //! numero utilise pour reperer l'inconnue lambda utilisee pour forcer \int w = 0
    int lambda_dof_number;
    //! deplacement du marteau ?
    Real_wp val_xsi;
    
    //! valeur de \int x phi_i pour chaque element du maillage
    Vector<TinyVector<Real_wp, nb_base> > ValueIx_Mesh;
    //! matrice de base elementaire \int phi_i phi_j
    TinyMatrix<Real_wp, Symmetric, nb_base, nb_base> MassMatrixBase;
    //! poids des fonctions de base \int phi_i
    TinyVector<Real_wp, nb_base> ValueI1_Base;
    
    //! partie lineaire de la matrice differentielle
    Matrix<Real_wp, General, ArrayRowSparse> linear_matrix;
    
    //! parametres de la tete de marteau
    Real_wp H_head, mass_head, Hprime_head, d0, y_wall;
    //! distance entre la corde et le centre de la liaison pivot du manche
    Real_wp distance_corde_pivot;
    //! d : distance du marteau a la corde, Fhead : force du marteau
    Real_wp last_d_stored, last_Fhead_stored, last_Fjack_stored, last_tangential_Fjack_stored, last_torque_jack_stored;
    //! marteau
    Hammer& hammer;
    //! si vrai, interaction du manche avec la corde
    bool is_interacting_string;
    //! si vrai, on utilise une formulation avec multiplicateur de Lagrange
    bool formulation_with_multiplier;
    //! numero du ddl pour la force Fcoupl (couplage tete marteau <-> manche),
    //! et le deplacement du marteau
    int Fcoupl_dof_number, xi_dof_number;
    int size_sismo_buffer;
    
    //! buffer de sortie pour d et la force du marteau
    WriteOnTheGoWithTinyBuffer<Real_wp, 3> output_Force_d;
    //! energy of the shank when it is removed from global system string+shank
    Real_wp free_kinetic_energy, free_potential_energy, free_internal_energy, free_shank_energy;
    //! first contact happened ?
    bool has_first_contact_happened;
    
    int print_level; //!< verbosity level
    
  public:
    Shank(Hammer& h);

    int GetM() const;

    void ComputeRightHandSide(VectReal_wp & scheme, Matrix<Real_wp> & val2,
                              bool newton_rhs);

    void ConstructAll(int num, bool double_prec = false, bool binary = false);
    void ConstructLinearMatrix();
    void FlushBuffer();

    void WriteInterp(int n, const Matrix<Real_wp>& val);
    void WriteSnapshots(int nb_iter, Real_wp t, const Matrix<Real_wp>& val);

    void SetTimeStep(const Real_wp& dt);
    bool AllowedKeyword(const string& keyword);

    void SetInputData(const string & keyword, const Vector<string>& parameters);

    void Init(int nb_iterations, Matrix<Real_wp>& val);

    void Advance(int nt, const Real_wp& t, Matrix<Real_wp>& val);

    void Display(ostream& file_out, bool write_on_file);

    Real_wp GetThetaShank(const Matrix<Real_wp>& val, int n);
    void SetThetaShank(Vector<Real_wp>& scheme, const Real_wp& x);

    Real_wp deriv_sinc(const Real_wp& x);
    Real_wp delta_cos(const Real_wp& a, const Real_wp& b);
    Real_wp deriv_delta_cos(const Real_wp& a, const Real_wp& b);
    Real_wp delta_sin(const Real_wp& a, const Real_wp& b);
    Real_wp deriv_delta_sin(const Real_wp& a, const Real_wp& b);

    void ComputeTM(int ne, int instant, const Matrix<Real_wp> & val,
                   TinyMatrix<Real_wp, General, dimension, nb_quad> & TM,
                   TinyMatrix<Real_wp, General, dimension, nb_quad> & dTM);

    void ComputeWbase(int ne, int instant, const Matrix<Real_wp> & val,
                      TinyVector<Real_wp, nb_base> & ubase,
                      TinyVector<Real_wp, nb_quad> & TM);

    void ComputeScheme(VectReal_wp & scheme, Matrix<Real_wp> & val);
    void ComputeScheme(VectReal_wp & scheme, VectReal_wp& rhs, Matrix<Real_wp> & val);

    void AddScheme(VectReal_wp & scheme, Matrix<Real_wp> & val,
                   int hammer_dof, const Real_wp& interaction_hammer, 
                   const string& file_name = string(""));

    void ForceConstraint(Matrix<Real_wp>& x, int hammer_dof = -1);

    Real_wp ComputeDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val);

    void AddDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val,
                 int hammer_dof, int force_hammer_dof);

    Real_wp Fct_Delta(Real_wp x);

    template<class GenericPb>
    void ComputeInterpolatedJackForce(GenericPb & var);

    template<class GenericPb>
    void ComputeInterpolatedTangentialJackForce(GenericPb & var);

    template<class GenericPb>
    void ComputeInterpolatedTorqueJack(GenericPb & var);
  
    Real_wp GetEnergy(const Matrix<Real_wp>& val, Real_wp& kinetic_energy,
                      Real_wp& potential_energy, Real_wp& internal_energy);
    
  };  
  
}

#define MONTJOIE_FILE_SHANK_HXX
#endif
