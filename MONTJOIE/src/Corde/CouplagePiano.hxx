#ifndef MONTJOIE_FILE_COUPLAGE_PIANO_HXX

namespace Montjoie
{
  
  //! classe de couplage entre la table d'harmonie, les cordes et l'air
  class CouplagePiano : public InputDataProblem_Base
  {
  public:
    //! inverse de la matrice de masse fluide
    /*!
      MhCouplage(i) = 1 / MhMinus(i) + 1 / MhPlus(i)
      ou MhPlus(i) = 1/(rho_f c^2) \int \varphi_i \varphi_i pour ddls de Gamma+
      MhMinus(i) = 1/(rho_f c^2) \int \varphi_i \varphi_i pour ddls de Gamma-
     */
    VectReal_wp MhCouplage;
    //! inverse de la matrice de masse fluide pour les ddls \Gamma+, et les ddls \Gamma-
    VectReal_wp invMhSurfMinus, invMhSurfPlus;
    //! vecteur Ch \Lambda^{n-1} pour les ddls de Gamma- et Gamma+
    VectReal_wp ChLambdaNm1_Minus, ChLambdaNm1_Plus;

#ifdef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
    //! objet pour la propagation du son dans l'air
    FluideExterieur var_fluide;
#endif
    
    //! chronometre utilise pour les temps de calcul
    MontjoieTimer chrono;
    //! si vrai, on affiche le temps max pour un seul proc
    bool display_max_time_single_proc;
    //! on affiche les temps de calcul tous les xxx iterations
    int nb_iter_for_each_computation_time;
    //! differents chronometres
    int CHRONO_NOTE, CHRONO_PLATE, CHRONO_COUPLAGE, CHRONO_ENERGY, CHRONO_FEM;
    //! objet pour la table d'harmonie
    Soundboard soundboard;
    //! si vrai, les modes propres de la table d'harmonie sont ecrits 
    bool display_eigenvectors_plate;
    //! energie (de l'air, de la plaque, corde, etc)
    TinyVector<Real_wp, 4> previous_energy;
    //! fonction Chi de distribution de l'attache du chevalet sur la table 
    SourceChi source_coupling_string;
    //! prefixe pour l'ecriture des fichiers de reprise
    string prefix_save;
    //! si vrai, les fichiers de reprise sont ecrits dans un repertoire local puis copies
    bool reprise_on_temporary_path;
    //! parametres de reprise
    ParameterOutputReprise output_reprise;
    
    //! cordes des pianos
#ifdef PIANO_STIFF_STRING
    MultiString<5, 5, WaveEquationStiffString> note;
#elif defined(PIANO_LINEAR_STIFF_STRING)
    MultiString<5, 5, WaveEquationStiffLINString> note;
#elif defined(PIANO_NONLINEAR_STIFF_STRING)
    MultiString<5, 5, WaveEquationStiffNLString> note;
#elif defined(PIANO_LINEAR_STRING)
    MultiString<5, 5, WaveEquationString> note;
#elif defined(PIANO_LINEAR_STRING_2T)
    MultiString<5, 5, WaveEquationStringLIN2T> note; 
#elif defined(PIANO_NONLINEAR_STRING_2T)
    MultiString<5, 5, WaveEquationStringNL2T> note;   
#elif defined(PIANO_NONLINEAR_STIFF_STRING_2T)
    MultiString<5, 5, WaveEquationStringStiffNL2T> note; 
#elif defined(PIANO_LINEAR_STIFF_STRING_2T)
    MultiString<5, 5, WaveEquationStringStiffLIN2T> note; 
#else
    MultiString<5, 5, WaveEquationNLString> note;
#endif
    
    //! buffer pour les sorties de l'energie
    WriteOnTheGoWithTinyBuffer<Real_wp, 4> output_ener_all_but_note;
    string file_output_energy;

    CouplagePiano(int argc, char**argv);
    
    void SetInputData(const string & keyword, const VectString & param);
    void Display();

    void ConstructAll(const string & data_file);
    void SaveDatas(int nt);

    void ConstructAllImpedance(const string & data_file);
    void ConstructAllStringSoundboard(const string & data_file);
    Real_wp WriteEnergy();

    void DisplayTimeChronos(int nt);
    void Advance(int nt, Real_wp temps, bool save_snapshot = true);
    void AdvancePlateFluid(int nt, Real_wp temps, bool save_snapshot = true);
    
    void AdvanceStringPlate(int nt, Real_wp temps, bool save_snapshot = true);

  };
  
}

#define MONTJOIE_FILE_COUPLAGE_PIANO_HXX
#endif
