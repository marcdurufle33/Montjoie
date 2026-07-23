#ifndef MONTJOIE_FILE_MULTI_STRING_HXX

namespace Montjoie
{

  class MultiString_Base : public InputDataProblem_Base
  {
  private :
    int dimension_;
    
  public :
    //! taille totale du systeme non-lineaire a resoudre
    int size;
    //! nombre de multiplicateurs de lagrange x nombre de cordes
    int size_LM;
    //! offsets utilises pour acceder aux differentes inconnues du systeme
    //! u(offset_hammer) : deplacement du marteau
    int offset_hammer, offset_hammer_aux, offset_force_hammer, offset_shank;
    int n_LaMu;
    
    //! nombre de multiplicateurs de Lagrange
    int nb_lagrange_mult;

    //! u(offset_lagrange_mult(n)): valeur du multiplicateur de Lagrange n 
    //! Le premier multiplicateur de Lagrange (n = 0) est la Force du chevalet dans la direction nu
    //! Le second multiplicateur de Lagrange (n = 1) est Force du chevalet dans la direction tau
    Vector<int> offset_lagrange_mult;

    //! nombre de cordes totales et nombre de cordes frappees par le marteau
    //! certaines cordes peuvent ne pas etre frappees par le marteau mais vibrer par sympathie
    int nb_corde, nb_corde_struck;

    //! offsets pour acceder aux inconnues des cordes
    Vector<int> nb_dof_before;

    Real_wp Deltat; //!< time step
    Real_wp Deltat2; //!< square of time step
    Real_wp invDeltat; //!< inverse of time step
    Real_wp invDeltat2; //!< square of the inverse of time step
    Real_wp invDeuxDeltat; //!< 1/(2 Delta t)
    
  public:
    MultiString_Base(int d);
    
    virtual bool IsAttachedString(int g) = 0;
    //~ virtual Real_wp GetNuString(int g, int k) = 0;
    //~ virtual Real_wp GetNuOrthogonalString(int g, int k) = 0;
    //~ virtual Real_wp GetNuOrthogonal2String(int g, int k) = 0;
    virtual Real_wp GetTransmissionChevaletString(int g, int k1, int k2) = 0;
    virtual Real_wp GetBridgeHeight(int i) = 0;
    virtual int GetNbStringDofs(int g) = 0;
    virtual int GetStringType() = 0;
    int GetNbStrings() const;
    int GetNbUnknownsPerString() const;
    
    int GetStringDofNumber(int i, int k, int j);
    int GetLMDofNumber(int i);
    int GetLaMuDofNumber(int n_LaMu, int i);
    int GetLMGDofNumber(int i);
    int GetLagrangeM_DofNumber(int i, int p);
    int GetHammerDofNumber();
    int GetHammerAuxDofNumber(int i);
	
    
    Real_wp GetValue(int i, int k, int j, int n, Matrix<Real_wp> & val2);
    Real_wp GetLM(int i, int n, Matrix<Real_wp> & val2);
    Real_wp GetLaMu(int n_LaMu, int i, int n, Matrix<Real_wp> & val2);
    Real_wp GetLMG(int i, int n, Matrix<Real_wp> & val2);
    Real_wp GetLagrangeM(int i, int p, int n, Matrix<Real_wp>& val2);
    Real_wp GetHammer(int n, Matrix<Real_wp> & val2);
    Real_wp GetHammerAux(int i, int n, Matrix<Real_wp> & val2); 

    
    void SetValue(int i, int k, int j, int n, Real_wp valeur, Matrix<Real_wp>& val2);
    void SetLM(int i, int n, Real_wp valeur, Matrix<Real_wp>& val2);
    void SetLaMu(int n_LaMu, int i, int n, Real_wp valeur, Matrix<Real_wp>& val2);
    void SetHammer(int n, Real_wp valeur, Matrix<Real_wp>& val2);
    void SetHammerAux(int i, int n, Real_wp valeur, Matrix<Real_wp>& val2);

    
    void SetSchemeString(int i, int k, int j, Real_wp valeur, VectReal_wp& scheme);
    void SetSchemeLM(int i, Real_wp valeur, VectReal_wp& scheme);
    void SetSchemeLaMu(int n_LaMu, int i, Real_wp valeur, VectReal_wp& scheme);
    void SetSchemeLMG(int i, Real_wp valeur, VectReal_wp& scheme);
    void SetSchemeHammer(Real_wp valeur, VectReal_wp& scheme);
    void SetSchemeHammerAux(int i, Real_wp valeur, VectReal_wp& scheme);
    
    void AddSchemeString(int i, int k, int j, Real_wp valeur, VectReal_wp& scheme);
    void AddSchemeLM(int i, Real_wp valeur, VectReal_wp& scheme);
    void AddSchemeLaMu(int n_LaMu, int i, Real_wp valeur, VectReal_wp& scheme);
    void AddSchemeHammer(Real_wp valeur, VectReal_wp& scheme);
    void AddSchemeHammerAux(int i, Real_wp valeur, VectReal_wp& scheme);
    
    void SetDiffStringString(int i1, int k1, int j1, Real_wp valeur,
                             int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffLMString(int i, Real_wp valeur, int i2,
                         int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffLMGString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffLaMuString(int n_LaMu, int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
        
    void SetDiffStringLM(int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffStringLMG(int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffStringLaMu(int n_LaMu, int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);     
    
    void SetDiffStringHammer(int i, int k, int j, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix);
    void SetDiffStringHammerAux(int i, int k, int j, Real_wp valeur, int l, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffHammerString(Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix);
    void SetDiffHammerAuxString(int l, Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffLMLM(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void SetDiffLMGLMG(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix); 
    
    void SetDiffLaMuLaMu(int n_LaMu, int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix);
        
    void SetDiffHammerHammer(Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix);
    void SetDiffHammerHammerAux(Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);
    void SetDiffHammerAuxHammer(int i, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix);
    void SetDiffHammerAuxHammerAux(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffStringString(int i1, int k1, int j1, Real_wp valeur,
                             int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffLMString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffLMGString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffLaMuString(int n_LaMu, int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffStringLM(int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffStringLaMu(int n_LaMu, int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffStringHammer(int i, int k, int j, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix);
    void AddDiffStringHammerAux(int i, int k, int j, Real_wp valeur, int l, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffHammerString(Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix);
    void AddDiffHammerAuxString(int l, Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix);
    
    void AddDiffHammerHammer(Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix);
    void AddDiffHammerHammerAux(Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix);
    void AddDiffHammerAuxHammer(int i, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix);
    void AddDiffHammerAuxHammerAux(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix);
    
  };
  

  //! class handling several strings struck by a same hammer
  template<int nb_base, int nb_quad, class TypeEquation>
  class MultiString : public MultiString_Base
  {
  public :
    //! number of components for the displacement of each string
    enum{dimension = TypeEquation::dimension};
    
    //! bridge coupled with the strings
    Bridge bridge;
    //! hammer striking the strings
    Hammer hammer;
    //! shank of the hammer
    Shank<nb_base, nb_quad, WaveEquationStiffString> shank;
    //! parameters of all the strings
    ParameterMultistring param_corde_data;
    //! list of strings struck by the same hammer
    Vector<PianoString<nb_base, nb_quad, TypeEquation> > piano_strings;
    
    //! initial time, final time, current time, and frequency of displays
    Real_wp t_begin, t_end, t_courant, t_display;
    //! directory where the outputs are written
    string DOSSIER;
    //! name of the data file
    string name_data_file;
    //! initial time, final time, Delta t between each snapshot
    Real_wp t_begin_interp , t_end_interp , deltat_interp;
    //! initial time, final time, Delta t between each sismogramm
    Real_wp t_begin_sismo, t_end_sismo, deltat_sismo;
    //! initial time, final time, Delta t between each computation of the energy
    Real_wp t_begin_ener, t_end_ener, deltat_ener;
    //! initial time, final time, Delta t between each output of the hammer crush
    Real_wp t_begin_crush, t_end_crush, deltat_crush;
    //! number of iterations in time
    int nb_max_iter;

    //! donnees pour reprendre le calcul apres une interruption
    bool load_reprise, save_reprise;
    //! donnees pour reprendre le calcul apres une interruption
    int load_iter_reprise, save_iter_reprise, bascule_reprise;
    //! chemin et fichier utilises pour les reprises
    string path_reprise, file_storage_param;
    
    Real_wp k_chirp;
    
    //! Somme des forces pour les differentes cordes
    Real_wp somme_interaction_hammer;
    
    //! buffer de sortie pour le marteau et la force du marteau
    WriteOnTheGoWithTinyBuffer<Real_wp,2> output_hammer, output_force;
    //! buffer de sortie pour l'energie
    WriteOnTheGoWithBuffer<Real_wp> output_ener;
    //! buffer de sortie pour les multiplicateurs de Lagrange
    WriteOnTheGoWithBuffer<Real_wp> output_LM;
    //! taille des buffers de sortie
    int size_buffer, size_buffer_ener;
    
    //! numero de l'iteration en temps
    int nb_iter_courant;
    
    //! nombre de combinaisons pour le schema temporel du terme non-lineaire
    int nb_combinaisons;
    //! PERMUT_ONE : une seule permutation utilisee pour le terme non-lineaire
    //! PERMUT_TWO_WAY : une permutation et son inverse
    //!         sont utilisees pour obtenir un schema symetrique (donc d'ordre 2)
    //! PERMUT_ALL : toutes les permutations sont utilisees et moyennes (schema d'ordre 2 aussi)
    enum {PERMUT_ONE, PERMUT_TWO_WAY, PERMUT_ALL};
    //! type de permutation utilisee (PERMUT_ONE, PERMUT_TWO_WAY ou PERMUT_ALL)
    int permutations;
    Matrix<int> Sigma; //!< matrice des combinaisons
    //! coefficient utilise pour chaque permutation
    Vector<Real_wp> Vect_zeta;
    
    //! schemas temporels disponibles
    enum{EXPLICIT, THETA_NL_HALF, THETA_NL_QUARTER};
    //! type de schema utilise (EXPLICIT, THETA_NL_HALF, etc) 
    int solver;
    //! parametres theta pour le schema en temps
    VectReal_wp theta;
    
    //! vaut 1 si on met une donnee initiale, 0 si on veut la corde au repos.
    bool existence_donnee_ini;
    //! type de conditions initiales disponibles
    enum{NONE, SINUS, DOUBLE_SINUS, GAUSSIAN, SINUSMULTI, NONZERO, PLUCK, PLUCKANGLE};
    int type_initial_condition; //!< condition initiale utilisee (NONE, SINUS ou GAUSSIAN)
    Real_wp amplitude_ini; //!< amplitude de la condition initiale
    int nb_mode_ini; //!< numero du mode pour la condition initiale
    Real_wp pluckangle_ini;
    //! composante du deplacement de la corde sur lequel on a la condition initiale
    int dimension_ini;
    Real_wp pluck_height, pluck_position;
    Vector<Real_wp> vect_amplitude_ini;
    Vector<int> vect_nb_mode_ini;
    Vector<int> vect_dimension_ini;
    
    int number_of_duplex_string;
    Vector<int> vect_duplex_string_ID;
    
    //! coefficients pour le couplage avec la table d'harmonie
    //     Real_wp scal_a, scal_b;
    //! coefficient de couplage ?
    //Real_wp coef_couplage;
 
    //! solution pour trois temps u^{n-1} = val(:, 0), u^{n} = val(:, 1), u^{n+1} = val(:, 2)
    Matrix<Real_wp> val; 
    
    //! solveur utilise pour la resolution du systeme non-lineaire
    MegaNewtonSolver Newton_solver;
    int print_level; //!< verbosity level
    //! measured times for the beginning and end of the force
    Real_wp beginning_time_force, last_time_force; 
    Real_wp final_time_shank; //!< last time for which the shank is resolved
    bool drop_hammer; //!< if true, the hammer is absent
    bool drop_hammer_when_possible; //!< if true, the hammer will be removed as soon as possible
    bool use_naive_non_linear_term; //!< if true, a naive scheme is used for non-linear terms
    bool separate_energy_uv; //!< if true the energy for u and v are separated
    bool shank_present_from_t0; //!< if true the shank is at the beginning of the simulation
    bool output_sismo_binary, output_sismo_double; // type of output for seismograms
    string file_output_LM, file_output_hammer, file_output_force, file_output_energy;

    bool BridgeHasHeight;
    bool useQuaternion;

    MultiString();
    
    void SetInputData(const string& keyword, const Vector<string>& param);
    
    void ConstructAll(int num = -1);
	
    void Display(ostream& file_out, bool write_on_file);
    
    // ***********
    //    Acces   *
    // ***********


    bool IsAttachedString(int g);
    //~ Real_wp GetNuString(int g, int k2);
    //~ Real_wp GetNuOrthogonalString(int g, int k2);
    //~ Real_wp GetNuOrthogonal2String(int g, int k2);
    Real_wp GetTransmissionChevaletString(int g, int k1, int k2);
    Real_wp GetBridgeHeight(int i);
    int GetStringType();
    int GetNbStringDofs(int g);
    
    void FlushBuffer();

    void Init(bool init_val = true);
    
    void SaveDatas(int nt);
    void RunAll();

    void Advance(Real_wp temps, int nt,
                    bool save_snapshot = true, bool compute_time = true);
    
    Real_wp WriteEnergy(bool comes_from_elsewhere = false);
    Real_wp ComputeEnergy(VectReal_wp &);
    Real_wp ComputeCouplingEnergy(int i);
    
    void ComputeTransmissionChevaletQuaternion(int i_, Real_wp angleStringPolarisation, Real_wp angle_chevalet, Real_wp angle_lateral = 0);  
     
    void WriteSnapshots(int nb_iter, Real_wp t);
    
    /// *************
    // * Schema   *
    // ************

    void ForceConstraint(Matrix<Real_wp>& x);    
    void ComputeScheme(VectReal_wp & scheme, VectReal_wp& rhs, Matrix<Real_wp> & val2);
    void ComputeRightHandSide(VectReal_wp & scheme, Matrix<Real_wp> & val2,
                              bool newton_rhs); 

    void ComputeDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val2); 

    void UpdateDiff(VirtualMatrix<Real_wp>& DiffMatrix, Matrix<Real_wp> & val2); 
    
    void ConstructLinearMatrix();
    
    Real_wp ComputeHammerInteraction(PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                                     Matrix<Real_wp> & val2);
    Vector<Real_wp> ComputeHammerInteractionDerivative(PianoString<nb_base, nb_quad,
                                               TypeEquation> & piano_string,
                                               Matrix<Real_wp> & val2);
    
    void ComputeNablaApproxU(TinyMatrix<Real_wp, General, dimension, nb_quad> & Qplus,
			     TinyMatrix<Real_wp, General, dimension, nb_quad> &Qminus,
			     TinyMatrix<Real_wp, General, dimension, nb_quad> &nablaU,
			     PianoString<nb_base, nb_quad, TypeEquation> & piano_string);
    
    void ComputeJacobApproxU(TinyMatrix<Real_wp, General, dimension, nb_quad> &Qplus,
			     TinyMatrix<Real_wp, General, dimension, nb_quad> &Qminus,
			     TinyArray3D<Real_wp, dimension, dimension, nb_quad> & JacobU,
			     PianoString<nb_base, nb_quad, TypeEquation> & piano_string);
    
    void FillQsigma(TinyMatrix<Real_wp, General, dimension, nb_quad> &Qplus,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & Qminus,
		    TinyVector<int, dimension-1>& sigma,
                    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma, int k);
    
    void Fillsigma(TinyVector<int, dimension-1>& sigma, Real_wp & zeta, int s, int k);
    
  };
  
  template<class GenericPb>
  void GenerateCombinaisons(GenericPb & var);
    
  void GenereCombinaisons(Matrix<int>& Sigma);
  
  Real_wp fact(int n);
  
  void CalculeZeta(Matrix<int>& sigma, Vector<Real_wp>& Vect_zeta);

}

#define MONTJOIE_FILE_MULTI_STRING_HXX
#endif

