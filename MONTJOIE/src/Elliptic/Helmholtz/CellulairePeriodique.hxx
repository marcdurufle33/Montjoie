#ifndef MONTJOIE_FILE_CELLULAIRE_PERIODIQUE_HXX

namespace Montjoie
{
  template<class TypeEquation>
  class VarHelmholtz_CellulairePeriodique : public VarHelmholtz_Cplx<TypeEquation>
  {

  public :
    typedef typename TypeEquation::Dimension Dimension;
       
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    // typedef typename Dimension::ArrayVectR_N ArrayVectR_N;
    //typedef typename Dimension::MatrixFullR_N MatrixFullR_N;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym;
    typedef typename TypeEquation::Complexe Complexe;
    typedef Vector<Complexe> VectComplexe;
    typedef Vector<VectComplexe> ArrayVectComplexe;
    typedef Matrix<Complexe, General, ArrayRowSparse> MatrixComplexe;
    
    // indices physiques (Maxwell)  (epsilon = rho for 2-D, and Helmholtz equation)
    Vector<ScalarPhysicalIndice<Dimension, Complexe> > ref_epsilon_max, ref_mu_max, ref_sigma_max;
    
    // epsilon, mu, sigma at the infinite
    Complexe eps_max_infty, mu_max_infty, sigma_max_infty;
    Complexe mu_infty;
    
    // nouveaux champs
    // troncature pour l operateur DtN
    int troncature ;
    
    // R0 rayon 
    Real_wp R0 ; 
    Real_wp ParametrePenalisation ;
    // Cas symetrique ou non 
    bool SymmetricCase ;
    bool CasDirichlet ;
    Complexe Sautldirichlet, Moyenneldirichlet ;
    
    // Pour les sources 
    int type_source_cellule ;
    enum {SRC_MU, SRC_ABSORBING, SRC_MU_X, SRC_MU_Y} ;
    
    // pour savoir si on a un comportement en z ou |z| (vrai si c est en |z|)
    int source_valeur_absolue ;
    // Periodes suivant X et Y
    Real_wp PeriodeX, PeriodeY ;
    
    // fichier ou ecrire les constantes
    string file_name_constant;
    
  public:
    // constructeurs et destructeurs 
    VarHelmholtz_CellulairePeriodique();
    
    // prise en compte des nouveaux champs R0 et troncature
    void SetInputData(const string& description_field, const VectString& parameters);
    
    inline string GetConstantFileName() const { return file_name_constant; }

    // surchage des fonctions init indices et setIndices
    // pour pouvoir crer ref_rho2 ( dans lequel sont stockees le vrai valeurs de rho)
    // alors que  ref_rho sera identiquement nul (elliptique = helmholtz + rho =0) 
    void InitIndices(int n);
    void SetIndices(int i, const VectString& parameters);
    void SetPhysicalIndex(const string& name_media, int i, const VectString& parameters);
    string GetPhysicalIndexName(int m) const;
    
    // calcul de B1 et B2
    // B1 = (omega^2 + R0 / \mu_inf) (int_{omega 1/2}  rho_inf -rho)
    // B2 = (1/mu_inf) (int_{omega 1/2} (mu_inf/R0) -mu(dU/dZ +1/R0))
    void CalculB1B2(const VectComplexe & U0, Complexe& B1, Complexe& B2);
    
    // Pour le 3D 
    // D11 = 1/T ( int_{omega 1/2} (mu dU/DX + mu) ) -1
    // D22 = 1/T ( int_{omega 1/2} (mu dU/DY + mu) ) -1
    // A3D = 1/T ( int_{omega 1/2} DU/DZ  ) -1
    void CalculD11D22A3D(const VectComplexe & U0, Complexe& D11, Complexe& D22, Complexe& A3D);

    // calcul de A 
    // A = 2 A_1^1/R0 = 2/R_0 (int_{gamma+} U(R+,theta) )- 2 R+/R_0 
    void CalculA(const VectComplexe & U0, Complexe& A, bool saut = false);
    
    //Calcul de D11
    // D11 = 1/T ( int_{omega 1/2} (mu dU/DX + mu) ) -1
    //void CalculD11(const VectComplexe & U0);

     //Calcul de D22
    // D22 = 1/T ( int_{omega 1/2} (mu dU/DY + mu) ) -1
    //void CalculD22(const VectComplexe & U0);

    //Calcul de A3D
    // A3D = 1/T ( int_{omega 1/2} DU/DZ  ) -1
    //void CalculA3D(const VectComplexe & U0);
    
    // operateur DtN specifique ( resolution de problemes elliptiques 
    // avec conditions de croissance quand R tend vers l infini 
    void AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>& mat_sp,
				   const GlobalGenericMatrix<Complex_wp>& nat_mat,
				   int offset_row = 0, int offset_col = 0);

    void AddBoundaryConditionTerms(VirtualMatrix<Real_wp>& mat_sp,
				   const GlobalGenericMatrix<Real_wp>& nat_mat,
				   int offset_row = 0, int offset_col = 0);
    
    // construction du second membre
    void ComputeRightHandSide(VectReal_wp & b_source, bool assemble = true); 
    void ComputeRightHandSide(VectComplex_wp & b_source, bool assemble = true); 
    void ComputeRightHandSide(Vector<VectComplexe> & b_source, bool assemble = true);  

    //  void ComputeRightHandSide(VectComplexe & b_source, bool assemble = true) const;
    // void ComputeRightHandSide(VectComplexe & b_source) const;   
    
  };
  
  template<class Dimension>
  class HelmholtzEquation_CellulairePeriodique: public HelmholtzEquation<Dimension>
  {
  public : 
    static inline bool ComputeDFjm1() { return true; }
    
  };
  
  template<class Dimension>
  class EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension> >
    : public  VarHelmholtz_CellulairePeriodique<HelmholtzEquation_CellulairePeriodique<Dimension> >
  {
  public:
    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
				 CondensationBlockSolver_Base<Real_wp>&,
				 const GlobalGenericMatrix<Real_wp>& nat_mat);

    void ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
				 CondensationBlockSolver_Base<Complex_wp>&,
				 const GlobalGenericMatrix<Complex_wp>& nat_mat);

  };


  
  template<class Dimension>
  class MuSource : public VirtualSourceFEM<Complex_wp, Dimension>
  {
  public :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    bool Source_x;
    
  private:    
    const EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension> >& vars;

  public:
    template<class TypeEquation>
    MuSource(const EllipticProblem<TypeEquation>& var);
    
    bool IsNonNullGradientSource(const VectR_N& s);
    
    void EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f);
    
  };
  
  template<class Dimension>
  class SourceAbsorbante : public VirtualSourceFEM<Complex_wp, Dimension>
  {
    public :
    typedef typename Dimension::R_N R_N;
    typedef typename Dimension::VectR_N VectR_N;
    bool valeurAbsolue ; 
    int ReferenceLineAbsorbing ;

  private:    
    const EllipticProblem<HelmholtzEquation_CellulairePeriodique<Dimension> >& vars;

  public:   
    template<class TypeEquation>
    SourceAbsorbante(const EllipticProblem<TypeEquation>& var);

    bool IsNonNullSurfacicSource ( int ref );
    
    void EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
				const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f);
    
  };
  
}
#define MONTJOIE_FILE_CELLULAIRE_PERIODIQUE_HXX
#endif
