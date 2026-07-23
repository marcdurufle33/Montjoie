#ifndef MONTJOIE_FILE_HAMMER_HXX

namespace Montjoie
{
  
  //! Class for the hammer
  class Hammer
  {
  public :
    //! if true the hammer is interacting with the string
    /*!
      and a differential equation is resolved for the hammer
     otherwise the hammer is like a right hand side (the force is given)
    */
    bool is_interacting;
    bool strike_at_an_angle;

    // bool pour l'utilisation du schéma avec variable auxiliaire
    bool with_auxiliary_variable;
    
    // COUPLING : interaction between the hammer and the string
    // through a law : M^H d^2 xsi/dt^2 = - F^H
    // where F^H = K^H phi(<u> - xsi) + R^H d/dt( phi(<u> - xsi) )
    // where phi(d) = Positive part(xsi_barre - d)^p
    // M^H : mass of the hammer
    // K^H : stiffness of the hammer
    // R^H : damping of the hammer
    // p : exponent of the hammer
    // <u> : displacement of the string
    // xsi : displacement of the hammer
    
    // EXPERIMENT : F^H is given by a file .dat
    //             usually values stored in the file come
    //             from experimental measures

    // RICKER, SINUS : not implemented
    
    // NONE : hammer is not present
    enum{COUPLING, EXPERIMENT, RICKER, SINUS, CHIRP, NONE};
    // type of interaction between the hammer and the string
    // (COUPLING, EXPERIMENT, etc)
    int type_interaction;
    
    // mass of the hammer M^H
    Real_wp masse_marteau;
    
    // center of the impact between the hammer and the string
    Real_wp point_impact_marteau;
    
    // initial speed of the hammer
    Real_wp vitesse_marteau;
    Real_wp exposant; // exposant de la loi
    Real_wp epaisseur_marteau; // epaisseur du marteau
    Real_wp pente_marteau; // pour la fonction de troncature du marteau
    Real_wp delta; // taille du marteau au repos
    Real_wp eloignement_marteau; // position initiale du marteau (xsi initial)
    Real_wp strike_angle; //strike_angle defines the angle where the hammer meets the 1st transverse string.

    
    Real_wp epsilon;// tolerance pour les differences finies

    enum{ORIGINAL, ROBOTRAN};
    int contact_model_with_string;

    bool regularised_contact_model;
    
    Vector<Real_wp> Ki; //!< les coefficients du marteau
    Vector<Real_wp> Ri; //!< les amortissements du marteau
    Vector<Real_wp> Ai, Bi, Ci, Di; //!< les coefficients du marteau dans le cas robotran
    Vector<Real_wp> ForceMarteauExt; // Force du marteau lue dans le fichier
    Vector<Real_wp> ForceMarteau; // Force du marteau apres interpolation (en temps)
    Real_wp DeltatForce; // \Delta t utilise dans le fichier ou la force est ecrite
    int nb_base_time; // ordre utilise pour interpoler la force du marteau
    
    GridInterpolation<Dimension1> grid_interp_time;// pour eventuelle force exterieure interpolee
    Globatto<Real_wp> base_time; // fonctions de base ( quadrature lobatto )
    Mesh<Dimension1> mesh_time; // maillage, points etc...
    Real_wp energy_free_hammer;
    
    Real_wp sin_amplitude;
    Real_wp sin_freq;
    Real_wp chirp_amplitude;
    Real_wp chirp_fstart;
    Real_wp chirp_fend;
    Real_wp chirp_phase;
    
  public:
    Hammer();

    void SetInputData(const string & description_field, const Vector<string> &parameters);

    template<int p, int q, class TypeEq>
    void SetInputData(const string & description_field, const Vector<string> &parameters,
                      PianoString<p, q, TypeEq>& string_p, ParameterMultistring& param_corde);

    void Display(ostream& file_out, bool write_on_file);

    template<class GenericPb, int p, int q, class TypeEq>
    void SetInitialConditions(GenericPb & var, Vector<PianoString<p, q, TypeEq> > & piano_strings);

    template<class GenericPb>
    Real_wp ComputeKineticEnergy(GenericPb & var);

    template<class GenericPb>
    Real_wp ComputePotentialEnergy(GenericPb & var, int i, Real_wp ux0p, Real_wp ux0m);

    template<class GenericPb>
    void ComputeInterpolatedForce(GenericPb & var);
    
    Real_wp Fct_Delta(Real_wp x);
    Real_wp Phi(Real_wp u);
    Real_wp PhiPrime(Real_wp u);
    Real_wp Psi(Real_wp u);

    Real_wp DiffFiniePsi(Real_wp u, Real_wp v);
    Real_wp DiffFiniePsiOPT(Real_wp u, Real_wp v);
    Real_wp DeriveeDiffFiniePsi(Real_wp u, Real_wp v );
    Real_wp DeriveeDiffFiniePsiOPT(Real_wp u, Real_wp v);

    // Guillaume //
    Real_wp Phi1Aux(Real_wp u, int i);
    Real_wp Phi2Aux(Real_wp u, int i);
    Real_wp GAux(Real_wp u, int i);
    Real_wp DeriveeGAux(Real_wp u, int i);
    Real_wp Echelon(Real_wp x);
  };
}

#define MONTJOIE_FILE_HAMMER_HXX
#endif
