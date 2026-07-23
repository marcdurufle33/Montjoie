#ifndef MONTJOIE_FILE_PIANO_STRING_HXX

namespace Montjoie
{
  
  //! classe gerant une seule corde de piano
  template<int nb_base, int nb_quad, class TypeEquation>
  class PianoString : public SpaceStringScheme<nb_base, nb_quad, TypeEquation> 
  {
  public :
    
    enum{dimension = TypeEquation::dimension};
    enum{NONE, SINUS, DOUBLE_SINUS, GAUSSIAN, SINUSMULTI,NONZERO, PLUCK, PLUCKANGLE};
    // si is_struck est vrai, la corde est frappee par un marteau
    // si is_attached est vrai, la corde est attachee
    bool is_struck, is_attached;
    //! numero de corde
    int numero_string;
    Real_wp L; // taille de la corde
    Real_wp E; // module d'Young
    Real_wp rho; // densite lineique de la corde
    Real_wp A; // section (m^2)
    Real_wp I; // module d'inertie
    Real_wp G; // module de cisaillement
    Real_wp k_prime; // param de Timoshenko
    TinyVector<Real_wp, 3> amo; // les coef d'amo fluide (var_eq.dimension est au plus 3)
    TinyVector<Real_wp, 3> amoB2; // les coef d'amo fluide dependant de la frequence (en omega^2)
    Real_wp T0; // tension de la corde
    Real_wp angle_chevalet, angleStringPolarisation, angle_chevalet_lateral; // angle du chevalet
    Real_wp rho_detuning; //density detuning, used to specify two different eigenfrequencies for two polarisations
    Real_wp bridge_height;
    Vector<Real_wp> vect_bridge_height;
    //~ TinyVector<Real_wp, dimension> NU, NUORTH, NUORTH2; // le vecteur de projection du chevalet sur la corde
    TinyMatrix<Real_wp, General, dimension, dimension> transmission_chevalet;
    //! matrices that will store linear terms of the scheme applied to U^n and U^(n-1)
    TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
               General, dimension, dimension> mat_Unm1, mat_Un, mat_Unp1;
    
    bool is_duplex;
    
    // objet specifiant l'equation de la corde
    TypeEquation var_eq;
    
    Real_wp ecrasement, force_marteau; // on les stocke pour les ecrire
    
    PianoString();
    
    void InitDefaultValues();
    void SetInputData(const string & description_field, const Vector<string> &parameters,
                      ParameterMultistring& param_corde_data);
    
    PianoString(const Vector<string> &parameters, Hammer& hammer,
                ParameterMultistring& param_corde_data);
    
    void FillPianoString(const Vector<string> &parameters, Hammer& hammer,
                         ParameterMultistring& param_corde_data);
    
    void Display(ostream&, bool);
    
    template<class GenericPb>
    void SetInitialConditions(int type_initial_condition, Real_wp amplitude,
                              int nb_mode, int dimension_ini, Real_wp Deltat, GenericPb & var);
    
    template<class GenericPb>
    void WriteStringSnapshot(int nb, Real_wp t, GenericPb & var, int i);
    
    template<class GenericPb>
    Real_wp ComputeScheme(VectReal_wp & scheme,Matrix<Real_wp> & val, GenericPb & var);

    template<class GenericPb>
    void ComputeRightHandSide(VectReal_wp & scheme, Matrix<Real_wp> & val,
                              GenericPb & var, bool newton_rhs);
    
    template<class GenericPb>
    void GetScalarProductAtBridge(Matrix<Real_wp>& val, GenericPb& var,
                                  Real_wp& somme0, Real_wp& somme2);
    
    template<class GenericPb>
    Real_wp ComputeDiff(VirtualMatrix<Real_wp>& DiffMatrix,
                        Matrix<Real_wp>& val, int num_xsi, GenericPb& var);

    template<class GenericPb>
    void UpdateDiff(VirtualMatrix<Real_wp>& DiffMatrix,
                        Matrix<Real_wp>& val, GenericPb& var);
   
    template<class GenericPb>
    void ComputeLinearMatrix(GenericPb& var);
    
    template<class GenericPb>
    Real_wp ComputeEnergy(GenericPb & var, Real_wp& ener_u, Real_wp& ener_v);
    
  };
  
}

#define MONTJOIE_FILE_PIANO_STRING_HXX
#endif

