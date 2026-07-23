#ifndef MONTJOIE_FILE_TIME_ELASTIC_HXX

namespace Montjoie
{ 
  //! base class for elastic equation in time-domain
  template<class Dimension>
  class TimeElasticEquation
  {
  public :
    typedef ElasticEquation<Dimension> TypeEquationStationary; //!< stationary equation
    
  };
  
  
  //! class to solve elaticity equations in time-domain
  template<class TypeEquation> 
  class TimeElastic_Base : public VarInstationary<TypeEquation>
  {
  public :
    
    // type declarations
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
    typedef typename TypeEquationStationary::Dimension Dimension; //!< dimension
    
    TimeElastic_Base();

    bool SplitSystem() const;
    bool IsDampedMedia(int ref);

    void InitTimeIterations();
    
    // Functions used to do interface with general schemes like Runge-Kutta
    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
				    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = true);

    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                                    const Real_wp& alpha, int level,
				    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = false);
    
    VirtualMassMatrix* GetNewMassMatrix();
    
  };
  
 
  //! mass matrix for continuous elements and acoustic equation
  template<class Dimension>
  class ContinuousMassMatrixElastic : public ContinuousUnsteadyMassMatrix<Dimension>
  {
  private:
    VarProblem<Dimension>& var_problem;
    VarBoundaryCondition<Real_wp, Dimension>& var_boundary;
    VarInstationary_Dim<Dimension>& var_time;
    
  public:
    template<class TypeEquation>
    ContinuousMassMatrixElastic(HyperbolicProblem<TypeEquation>& var);
    
    void Init(bool compute_time = true);    
  };
  
   
  //! class to solve elaticity equations in time-domain with H1 formulation
  template<class Dimension>
  class HyperbolicProblem<TimeElasticEquation<Dimension> >
    : public TimeElastic_Base<TimeElasticEquation<Dimension> >
  {
  public:
    void MltAddStiffnessScalar(const Real_wp& alpha, int level, const VectReal_wp& B,
			       const Real_wp& beta, VectReal_wp& C);
    
    void MltAddStiffnessVectorial(const Real_wp& alpha, int level, const VectReal_wp& B,
				  const Real_wp& beta, VectReal_wp& C);
    
  };
  
}

#define MONTJOIE_FILE_TIME_ELASTIC_HXX
#endif

