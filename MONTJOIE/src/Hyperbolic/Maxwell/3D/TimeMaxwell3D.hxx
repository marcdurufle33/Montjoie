#ifndef MONTJOIE_FILE_TIME_MAXWELL3D_HXX

namespace Montjoie
{
  //! class to specify 3-D time-domain Maxwell equations and discontinuous galerkin
  class TimeMaxwellEquation_3D_DG
  {
  public :
    typedef StaticMaxwellEquation_3D_DG TypeEquationStationary; //!< stationary equation
    
  };
  
  //! class to specify 3-D time-domain Maxwell equations for edge finite elements
  class TimeMaxwellEquation_3D
  {
  public :
    typedef StaticMaxwellEquation_3D TypeEquationStationary; //!< stationary equation
    
  };


  //! class to specify 3-D time-domain Maxwell equations for HDG formulation
  class TimeMaxwellEquationHdg_3D
  {
  public :
    typedef StaticMaxwellEquationHdg_3D TypeEquationStationary; //!< stationary equation
    
  };

  //! base class to solve time-domain Maxwell equations in 3-D
  template<class TypeEquation>
  class TimeMaxwell_3D: public VarInstationary<TypeEquation>
  {
  public:
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary; 
    //!< stationary equation
    
    TinyVector<VectReal_wp, 3> extrapolH, Hneighbor;
    TinyVector<VectReal_wp, 3> extrapolE, Eneighbor;
    Vector<TinyVector<VectReal_wp, 3> > evalEn_quad, evalHn_quad;
    
    void InitTimeIterations();
    void RunTimeIterations();

    virtual void GetExtrapolationMaxwell3D(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
					   bool extrapol_u, bool extrapol_v);
    
  };  

    
  //! class to solve time-domain Maxwell equations with edge finite elements in 3-D
  template<>
  class HyperbolicProblem<TimeMaxwellEquation_3D> 
    : public TimeMaxwell_3D<TimeMaxwellEquation_3D>
  {
  public :
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X, VectReal_wp& Y,
				    bool invert_mass = true, bool source = true);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
				    const VectReal_wp& X, VectReal_wp& Y,
				    bool invert_mass = true, bool source = false);
    
    void MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			       const Real_wp& beta, VectReal_wp& C);
    
    void MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
				  const Real_wp& beta, VectReal_wp& C);

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
    
  };

  //! class to solve time-domain Maxwell equations with DG method
  template<>
  class HyperbolicProblem<TimeMaxwellEquation_3D_DG> 
    : public TimeMaxwell_3D<TimeMaxwellEquation_3D_DG>
  {
  public :
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X, VectReal_wp& Y,
				    bool invert_mass = true, bool source = true);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
				    const VectReal_wp& X, VectReal_wp& Y,
				    bool invert_mass = true, bool source = false);
    
    void ComputeExtrapolationVectorial(const VectReal_wp& En, const VectReal_wp& Hn);
    void ComputeExtrapolationScalar(const VectReal_wp& En, const VectReal_wp& Hn);
    
    void GetExtrapolationMaxwell3D(const VectReal_wp& U, const VectReal_wp& V, int lvl, 
				   bool extrapol_u, bool extrapol_v);
    
    void MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			       const Real_wp& beta, VectReal_wp& C);
    
    void MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
				  const Real_wp& beta, VectReal_wp& C);
    
    void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma);
    
  };


  //! class to solve time-domain Maxwell equations with HDG method
  template<>
  class HyperbolicProblem<TimeMaxwellEquationHdg_3D> 
    : public TimeMaxwell_3D<TimeMaxwellEquationHdg_3D>
  {
  protected:
    VectReal_wp evalLambda_n;
    
  public :
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X, VectReal_wp& Y,
				    bool invert_mass = true, bool source = true);
    
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
				    const VectReal_wp& X, VectReal_wp& Y,
				    bool invert_mass = true, bool source = false);

    void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp& rho, Real_wp& sigma);
    
    Real_wp GetCoefficientTauHDG(int ref) const;

    void AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
				int nb_deriv, VectReal_wp& b_src);
    
    void MltAddStiffness(const Real_wp& alpha, int level, int level2,
			 const VectReal_wp& Uh, const VectReal_wp& Vh,
			 const VectReal_wp& Lambda, const Real_wp& beta,
			 VectReal_wp& ProdUh, VectReal_wp& ProdVh);
    
    void GetExtrapolationLambda(const VectReal_wp& Uh, const VectReal_wp& Vh, int level, VectReal_wp& Lambda);

    void AddTimeSourceHDG(const Real_wp& alpha, const Real_wp& t, int n, VectReal_wp& Y, VectReal_wp& Lambda);

    VirtualMassMatrix* GetNewMassMatrix();
    
  };


  //! mass matrix for HDG and Maxwell's equations
  class DiscontinuousMassMatrixHdgMaxwell3D : public DiscontinuousUnsteadyMassMatrix<Dimension3>
  {
  private:
    VarProblem<Dimension3>& var_problem;
    VarInstationary_Dim<Dimension3>& var_time;
    VarBoundaryCondition_Base& var_boundary;
    HarmonicMaxwell_3D<Real_wp>& var_maxwell;

    int type_matrix;
    
  public:
    template<class TypeEquation>
    DiscontinuousMassMatrixHdgMaxwell3D(HyperbolicProblem<TypeEquation>& var);
    
    void Init(bool compute_time = true);

    void ComputeOperatorCh();
    void ComputeOperatorBhGeom();
    
    void SolveOperatorCh(VectReal_wp& Y);
    
  };

}

#define MONTJOIE_FILE_TIME_MAXWELL3D_HXX
#endif
