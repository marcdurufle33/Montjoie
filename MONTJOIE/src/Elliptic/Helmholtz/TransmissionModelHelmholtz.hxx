#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_HELMHOLTZ_HXX

namespace Montjoie
{

  //! base class for transmission conditions and Helmholtz equation (see Delourme's thesis)
  template<class Complexe, class Dimension>
  class VarTransmission_Helm : public VarTransmission_Base<Dimension>
  {
  protected :
    typedef Vector<typename Dimension::R_N> VectR_N;
    
    //! coefficients A, B1, B2 
    Complexe AalphaDelta, B1alphaDelta, B2alphaDelta ;
    bool ModeleCentre; //!< if true a centered model is used
    bool CasDirichlet; //!< if true the model will approximate obstacles with Dirichlet condition
    Complexe A, B, C, mu0, rho0;
    bool ModeleTest;
    
  private:
    VarProblem<Dimension>& var_problem;
    VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary;
    VarHelmholtzIndex_Base<Complexe, Dimension>& var_helm;
    
    void InitDefaultValues();
    
  public :
    template<class TypeEquation>
    VarTransmission_Helm(EllipticProblem<TypeEquation>&) ;
    
    void SetInputData(const string&, const VectString& param);
    
    void InitTransmission();
    
    void AddTransmissionTerms(const Complex_wp& alpha,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      VirtualMatrix<Complex_wp>& mat_sp, int, int);

    void AddTransmissionTerms(const Real_wp& alpha,
			      const GlobalGenericMatrix<Real_wp>& nat_mat,
			      VirtualMatrix<Real_wp>& mat_sp, int, int);
    
  };


  //! class for transmission conditions and Helmholtz equation (see Delourme's thesis)
  template<>
  class VarTransmission<HelmholtzEquation<Dimension2> > 
    : public VarTransmission_Helm<Complex_wp, Dimension2>
  {
  public :
    VarTransmission(EllipticProblem<HelmholtzEquation<Dimension2> >& var);
    
  };
  
} // end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_HELMHOLTZ_HXX
#endif

  
