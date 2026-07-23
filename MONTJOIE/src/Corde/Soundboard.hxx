#ifndef MONTJOIE_FILE_SOUNDBOARD_HXX

#include "Oscillator.hxx"

namespace Montjoie
{
  
  //! source used for the coupling between the soundboard and the string
  class SourceChi : public VirtualSourceField<Real_wp, Dimension2>
  {
  public :
    R2 x0;
    Real_wp r0, alpha, ampl, beta;
    
    SourceChi();
    
    void EvaluateFunction(const R2& x, Vector<Real_wp>& f) const;
        
    void SetInputData(const string& description_field, const VectString& parameters,
                      ParameterMultistring& param_corde);
    
  };
  
  
  //! main class for handling the soundboard (computation and evolution of eigen-modes)
  class Soundboard : public Oscillator_Base, public InputDataProblem_Base
  {
  public : 

    typedef DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> MatrixDirect;
    
    //! nombre de ddls elements finis sur la plaque (pour une inconnue)
    int nFEM;
    string init_file; //!< contient le nom du fichier ini pour le soundboard.
    //! second membre (source mecanique sur la plaque)
    /*!
      RhsMecanic est la decomposition de la source sur les ddls
      RhsMecanicModal est la decomposition de la source sur les modes propres
     */
    VectReal_wp RhsMecanic, RhsMecanicModal;
    
    //!< et \dot{R}_{\Delta t} transpose(P_h) (J_h,0)
    VectReal_wp LambdaVolLambda, LambdaJ; //!< pour MAJ \Lambda^{n+1/2} avec juste coef_couplage
    VectReal_wp VectSecondMembre; //!< le vecteur (J//0)
    //! coefficients de la loi d'amortissement des modes de la plaque
    /*!
      L'amortissement pour chaque mode est calcule suivant la loi
      amo = a f^2 + b f + c
      a = VectParamAmo(0)
      b = VectParamAmo(1)
      c = VectParamAmo(2)
     */
    TinyVector<Real_wp, 3> VectParamAmo;
    //! chronometre
    MontjoieTimer chrono;
    
    //! table d'harmonie
    EllipticProblem<ReissnerMindlinEquation> VarSoundboard;
    //! ValLambdaSismoPoint(i, k) = mode_k(xi_i)
    /*!
      mode_k est le k-eme mode propre
      xi_i est le i-eme point ou on demande un sismogramme
      Cette matrice permet de calculer rapidement les valeurs du deplacement
      de la table sur un nombre faible de points.
     */
    Matrix<Real_wp> ValLambdaSismoPoint;
    
    //! les vecteurs propres sont charges, sauvegardes ?
    bool save_eigenvectors_asked, load_eigenvectors_asked;
    //! repertoire ou sont stockes les vecteurs propres
    string path_eigenvectors;
    
    //! nombre de valeurs propres calcules par run
    /*!
      Il est souvent plus efficace de calculer les valeurs propres
      en differents runs (par exemple 200 valeurs propres)
      L'ensemble des valeurs propres est ensuite reconstitue
     */
    int nb_eigenvalues_per_run;
    //! valeurs propres et parties imaginaires (nulles)
    VectReal_wp	eigen_values, lambda_imag;
    //! vecteurs propres
    Matrix<Real_wp, General, ColMajor> eigen_vectors;
    Matrix<Real_wp, General, RowMajor> Pu; //!< les vect propres que pour la partie sur U
    
#ifdef SELDON_WITH_SCALAPACK
    //! matrice contenant la factorisation LDL^T du complement de schur
    //DistributedMatrix<Real_wp, General, ColMajor> LUSchurP, Xloc_schur;
    Matrix<Real_wp, Symmetric, RowSymPacked> LUSchurP;
    //! taille du bloc utilise dans Scalapack
    int size_block_schur;
#else
    Matrix<Real_wp, Symmetric, RowSymPacked> LUSchurP;
    //!< matrice contenant la facto LDL^T de la matrice du complement de Schur
#endif
        
    //! buffer de sortie des sismogrammes mecaniques
    WriteOnTheGoWithBuffer<Real_wp> output_mecanic_sismo;
    //! taille du buffer de sortie des sismos mecaniques
    int size_buffer_mecanic_sismo;
    
    Soundboard();
    
    void SetInputData(const string & keyword, const VectString & param);
        
    void ActualiseLambda(Real_wp coef, Real_wp pulse = 0);
    
    void RunSoundboardDiagonalization(string DOSSIER, SourceChi&,
                                      Vector<string>&, bool remove_file = true);
    void PerformDiagonalization(EllipticProblem<ReissnerMindlinEquation>& var_soundboard,
                                const string& prefix, const string& prefix_save,
                                const Vector<string>&, bool, bool, SourceChi&,
                                bool display_eigen = false);
    
    void AssembleModalVector(VectReal_wp& VectSecondMembreModal);
    
    void ConstruitVectSecondMembreModal();
    Real_wp ConstruitLambdaSchur();
    Real_wp ConstructAll();
    void Initialise(Real_wp rhodelta);
    
    void ActualiseLambdaSchur(Real_wp coef, VectReal_wp &, VectReal_wp&);
    void ActualiseLambda(Real_wp coef,VectReal_wp &SecondMembreP);
    
    template<class VibroPb>
    void CalculeLambda(const Real_wp& coef, const Real_wp& pulse,
                       const VectReal_wp& Gvol_minus,
		       const VectReal_wp& Gvol_plus, const VectReal_wp& Mh,
		       VibroPb& var_vibro, VectReal_wp &LambdaNew);
    
    template<class VibroPb>
    Real_wp ComputeScal_a(const VectReal_wp & Gvol_minus, const VectReal_wp&,
                          VectReal_wp&, VectReal_wp&, VibroPb& var_vibro);
    
    template<class VibroVar>
    void CalculeLU(VectReal_wp& Mh, const VibroVar& var);
    
    void SolveSchur(VectReal_wp& X);

    void InitOutput(EllipticProblem<ReissnerMindlinEquation> & var_soundboard);
    
    void WriteOutputFEM(string toto, EllipticProblem<ReissnerMindlinEquation>& var_soundboard);
    
    void WriteOutputFEM(int nt, Real_wp temps, EllipticProblem<ReissnerMindlinEquation> & var_soundboard);
    
    void WriteOutputVect(VectReal_wp& Vect, string toto,
			 EllipticProblem<ReissnerMindlinEquation> & var_soundboard);
    
  };
  
}

#define MONTJOIE_FILE_SOUNDBOARD_HXX
#endif
