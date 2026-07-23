#ifndef MONTJOIE_FILE_TIME_ACOUSTIC_HXX

namespace Montjoie
{ 

  //! base class for acoustic equation
  template<class Dimension>
  class AcousticEquation
  {
  public :
    typedef LaplaceEquation<Dimension> TypeEquationStationary; //!< stationary equation
    
  };

  
  //! base class for acoustic equation and discontinuous galerkin
  template<class Dimension>
  class AcousticEquationDG
  {
  public :
    typedef LaplaceEquationDG<Dimension> TypeEquationStationary; //!< stationary equation
    
  };
  
  
  //! class to solve acoustic equation
  /*!
    acoustic equation can be written as
    \f$ \rho \frac{\partial^2 u}{\partial t^2} + \sigma \frac{\partial u}{\partial t} 
    - div( \mu grad(u)) = 0 \f$
    + boundary conditions + PML layers
  */
  template<class Dimension>
  class TimeAcoustic_Base
  {
  public :
    // type declarations
    typedef typename Dimension::R_N R_N; //!< R2 or R3
    typedef typename Dimension::VectR_N VectR_N; //!< vector of points
    typedef typename Dimension::MatrixN_N MatrixN_N; //!< 2x2 or 3x3 matrix
    typedef typename Dimension::MatrixN_Nsym MatrixN_Nsym; //!< 2x2 or 3x3 symmetric matrix
    
    TinyVector<VectReal_wp, Dimension::dim_N> extrapolV, Vneighbor;
    TinyVector<VectReal_wp, 1> extrapolU, Uneighbor;
    Vector<TinyVector<VectReal_wp, Dimension::dim_N> > evalUn_quad;
    Vector<TinyVector<VectReal_wp, Dimension::dim_N> > evalVn_quad;
    VectReal_wp evalLambda_n;
    
#ifdef MONTJOIE_WITH_TIME_REVERSAL
    // object dealing with time-reversal experiments
    TimeAcousticReversal<Dimension> time_reversal;
#endif
    
  private:
    VarInstationary_Fem<Dimension>& var_time;
    VarProblem<Dimension>& var_problem;
    VarHelmholtz_Base<Real_wp, Dimension>& var_helm;
    VarOutputProblem_Dim<Dimension>& var_output;
    VarBoundaryCondition_Base& var_boundary;
    VarSourceProblem_Cplx<Real_wp, Dimension>& var_source;
    
  public:
    template<class TypeEquation>
    TimeAcoustic_Base(HyperbolicProblem<TypeEquation>& var);
    
    bool IsDampedMedia(int ref);
    void InitTimeIterations();
    
    void SetInputData(const string& description_field, const VectString& parameters);
    void GetMemoryUsed(map<string, size_t>& var) const;

    // Functions used to do interface with general schemes like Runge-Kutta
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
				    const VectReal_wp& X, VectReal_wp& Y,
				    bool invert = true, bool source = true);
    
    void EvaluateFunctionS(const Real_wp& tn,const VectReal_wp& Y, const VectReal_wp& Yp,
                          VectReal_wp& ProdY, bool invert = true, bool source = true);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
                                    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert = true, bool source = false);
    
    void ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Vh,
			       const Real_wp& beta, VectReal_wp& Yh, bool extrapol = true);
    
    void ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
				  const Real_wp& beta, VectReal_wp& Yh, bool extrapol = true);
    
    void GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void GiveVectorialIterate(int nb_iter, const Real_wp& t, VectReal_wp Vh);

    virtual void GetExtrapolationAcoustic(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
					  bool extrapol_u, bool extrapol_v) = 0;

    virtual void MltAddStiffness(const Real_wp& alpha, int level, int level2,
				 const VectReal_wp& Uh, const VectReal_wp& Vh,
				 const VectReal_wp& Lambda, const Real_wp& beta,
				 VectReal_wp& ProdUh, VectReal_wp& ProdVh);
    
    virtual void GetExtrapolationLambda(const VectReal_wp& Uh, const VectReal_wp& Vh, int level, VectReal_wp& Lambda);

    void AddTimeSourceHDG(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& Y, VectReal_wp& Lambda);
    
  };

  
  //! mass matrix for discontinuous elements and acoustic equation
  template<class Dimension>
  class DiscontinuousMassMatrixAcoustic : public DiscontinuousUnsteadyMassMatrix<Dimension>
  {
  private:
    VarProblem<Dimension>& var_problem;
    VarInstationary_Dim<Dimension>& var_time;
    VarBoundaryCondition_Base& var_boundary;
    VarHelmholtz_Base<Real_wp, Dimension>& var_helm;
    
  public:
    template<class TypeEquation>
    DiscontinuousMassMatrixAcoustic(HyperbolicProblem<TypeEquation>& var);
    
    void Init(bool compute_time = true);    
  };


  //! mass matrix for continuous elements and acoustic equation
  template<class Dimension>
  class ContinuousMassMatrixAcoustic : public ContinuousUnsteadyMassMatrix<Dimension>
  {
  private:
    VarProblem<Dimension>& var_problem;
    VarInstationary_Dim<Dimension>& var_time;
    VarBoundaryCondition_Base& var_boundary;
    
  public:
    template<class TypeEquation>
    ContinuousMassMatrixAcoustic(HyperbolicProblem<TypeEquation>& var);

    void Init(bool compute_time = true);    
  };
  
  
  //! class for solving acoustics equation
  template<class TypeEquation>
  class TimeAcoustic_Eq : public VarInstationary<TypeEquation>,
			  public TimeAcoustic_Base<typename TypeEquation::TypeEquationStationary::Dimension>
  {
    typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  public:
    TimeAcoustic_Eq();

    void GetMemoryUsed(map<string, size_t>& var) const;
    
    void InitTimeIterations();
    
    void AddScalarSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y);
    void AddVectorialSourceAtTime(const Real_wp& alpha, const Real_wp& tn, int n, VectReal_wp& Y);
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    // Functions used to do interface with general schemes like Runge-Kutta
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
				    const VectReal_wp& X, VectReal_wp& Y,
				    bool invert = true, bool source = true);
    
    void EvaluateFunctionS(const Real_wp& tn,const VectReal_wp& Y, const VectReal_wp& Yp,
                          VectReal_wp& ProdY, bool invert = true, bool source = true);

    void EvaluateFunction(const Real_wp& tn,const VectReal_wp& Y,
                          VectReal_wp& ProdY, bool invert = true, bool source = true);
    
    void EvaluateFunction(const Real_wp& tn, const Real_wp& alpha, int level,
                          const VectReal_wp& X, VectReal_wp& Y,
                          bool invert_mass = true, bool source = false);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
                                    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert = true, bool source = false);
    
    void ApplyOperatorRhScalar(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Vh,
			       const Real_wp& beta, VectReal_wp& Yh, bool extrapol = true);
    
    void ApplyOperatorRhVectorial(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
				  const Real_wp& beta, VectReal_wp& Yh, bool extrapol = true);
    
    void GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y);
    void GiveVectorialIterate(int nb_iter, const Real_wp& t, VectReal_wp& Vh);
    
    void ComputeExtrapolationVectorial(const VectReal_wp& Uh, const VectReal_wp& Vh);
    void ComputeExtrapolationScalar(const VectReal_wp& Uh, const VectReal_wp& Vh);
    
    bool SplitSystem() const;
    
    VirtualMassMatrix* GetNewMassMatrix();
    
    void GetExtrapolationAcoustic(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
				  bool extrapol_u, bool extrapol_v);
    
  };
  
  
  //! class to solve acoustic equation with H1 elements
  template<class Dimension>
  class HyperbolicProblem<AcousticEquation<Dimension> >
    : public TimeAcoustic_Eq<AcousticEquation<Dimension> >
  {
  public :
    void AddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
                                   const VectReal_wp& Uh, VectReal_wp& Vh,
                                   const Real_wp& coef_mu, const Real_wp& coef_sigma);

    void SolveElementStiffnessTensorGen(int i, int& offset, const Real_wp& alpha,
                                        const VectReal_wp& Uh, VectReal_wp& Vh,
                                        const Real_wp& coef_sigma, bool add);

    void SolveElementStiffnessTensor(int i, int& offset, VectReal_wp& Vh,
                                     const Real_wp& coef_sigma);
    
    void SolveAddElementStiffnessTensor(int i, int& offset, const Real_wp& alpha,
                                        const VectReal_wp& Uh, VectReal_wp& Vh,
                                        const Real_wp& coef_sigma);
    
    void MltAddStiffnessScalar(const Real_wp&, int lvl, const VectReal_wp&,
			       const Real_wp&, VectReal_wp&);
    
    void MltAddStiffnessVectorial(const Real_wp&, int lvl, const VectReal_wp&,
				  const Real_wp&, VectReal_wp&);

  };  
  
  
  //! class to solve acoustic equation with LDG formulation
  template<class Dimension>
  class HyperbolicProblem<AcousticEquationDG<Dimension> >
    : public TimeAcoustic_Eq<AcousticEquationDG<Dimension> >
  {
  public :
    void MltAddStiffnessScalar(const Real_wp&, int lvl, const VectReal_wp&,
			       const Real_wp&, VectReal_wp&);
    
    void MltAddStiffnessVectorial(const Real_wp&, int lvl, const VectReal_wp&,
				  const Real_wp&, VectReal_wp&);
    
    void MltAddStiffness(const Real_wp& alpha, int level, int lvl2, const VectReal_wp& Uh, const VectReal_wp& Vh,
			 const VectReal_wp& Lambda, const Real_wp& beta, VectReal_wp& ProdUh, VectReal_wp& ProdVh);
    
    void GetExtrapolationLambda(const VectReal_wp& Uh, const VectReal_wp& Vh, int level, VectReal_wp& Lambda);
    
    void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp&, Real_wp&);

    Real_wp GetCoefficientTauHDG(int ref) const;

    void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
				int nb_deriv, VectReal_wp& b_src);
    
    void GetExtrapolationLambdaHDG(int level, const VectReal_wp& Y, VectReal_wp& ProdL);
    void MltAddStiffnessHDG(const Real_wp& alpha, int level, int lvl2, const VectReal_wp& Y, const VectReal_wp& L,
			    const Real_wp& beta, VectReal_wp& ProdY);

    void MltAddLambdaHDG(const Real_wp& alpha, int level, const VectReal_wp& L, VectReal_wp& ProdY);
    void MltAddLambdaHDG(const Complex_wp& alpha, int level, const VectComplex_wp& L, VectComplex_wp& ProdY);
    
  };  
    
}

#define MONTJOIE_FILE_TIME_ACOUSTIC_HXX
#endif

