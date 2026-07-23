#ifndef MONTJOIE_FILE_OSCILLATOR_HXX

namespace Montjoie
{

  class Oscillator_Base
  {
  public :
    //! pas de temps
    Real_wp Deltat;
    //! nombre de valeurs propres
    int nev;
    //! Fichiers ou sont stockes les coefficients de masse, rigidite et amortissement
    string mass_file_name, lambda_file_name, amo_file_name, vect_second_membre_file_name;
    bool vect_second_membre_is_given;

    //    Real_wp scal_b; //!< contient sum (MatriceResolExactSecondMembre)
    //!< ou sum(lambdaJ*VectSecondMembreModal) si ConstruitLambdaSchur est apelle. 
    //Real_wp scal_a; //!< coefficient de couplage analogue a scal_b

    //! vecteurs ou sont stockes les coefficients de masse, rigidite et amortissement
    VectReal_wp vect_lambda, vect_amortissement, vect_masse;

    VectReal_wp  VectSecondMembreModal; //!< tr(P) (J//0)
    bool write_oscillator_param_bool;    
    VectReal_wp Lambda; //!< contient Lambda^{n-1/2}
    VectReal_wp LambdaPoint; //!< contient \partial_t Lambda^{n-1/2}

    VectReal_wp Lambda0, Lambda0Point; //!< contient R_{\Delta t} transpose(P_h) (J_h,0)

    VectReal_wp MatriceResolExacte0, MatriceResolExacte0Point;//!< matrice contenant S^0_{\Delta t}
    VectReal_wp MatriceResolExacte1, MatriceResolExacte1Point;//!< matrice contenant S^1_{\Delta t}
    VectReal_wp MatriceResolExacteSecondMembre, MatriceResolExacteSecondMembrePoint;
    //!< matrice contenant R_{\Delta t}

    VectReal_wp RDeltaDemi ; //!< contient R_{∆t}^{1/2}

    //! buffer de sortie des Lambda^(n-1/2)
    WriteOnTheGoWithBuffer<Real_wp> output_lambda;
    //! parametres de sortie des Lambda^(n-1/2)
    Real_wp t_begin_output_lambda, t_end_output_lambda, deltat_output_lambda;
    //! taille du buffer de sortie des Lambda^(n-1/2)
    int size_buffer_lambda, N_output_Lambda;

  public:
    Oscillator_Base();
    
    void Display(ostream& file_out, bool write_on_file);
    
    void SetInputData(const string & keyword, const VectString & param);
    
    void ReadImpedance(const string& DOSSIER);
    
    Real_wp ConstruitMatricesResolExacte();
    void ConstruitLambda0etPoint();
    
    Real_wp CalculeScal_a();
    
    Real_wp GetKineticEnergy();
    Real_wp GetPotentialEnergy();

    void ActualiseLambda(Real_wp coef);
    
    void WriteOutputLambda(int nt, Real_wp temps);
  };

}

#define MONTJOIE_FILE_OSCILLATOR_HXX
#endif

