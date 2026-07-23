#ifndef MONTJOIE_FILE_TIME_MAXWELL_AXI_HXX

namespace Montjoie
{
  
  //! specifications for time-domain Maxwell equations in axisymmetric domain
  class TimeMaxwellEquation_Axi_DG
  {
  public :
    typedef StaticMaxwellEquation_Axi_DG TypeEquationStationary; //!< stationary equation

  };
  
  //! specifications for time-domain Maxwell equations in axisymmetric domain
  class TimeMaxwellEquation_HcurlAxi
  {
  public :
    typedef StaticMaxwellEquation_HcurlAxi TypeEquationStationary; //!< stationary equation
    
  };
  
  
  //! class to solve time-domain Maxwell equations in axisymmetric domain
  template<class TypeEquation>
  class TimeMaxwell_Axi : public VarInstationary<TypeEquation>
  {
  public:    
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
    //!< stationary equation
    
    //! multiplication by -1 for dofs E_theta
    template<class Matrix1>
    void SetMinus_DofTeta(const VectBool& IsDofTeta, Matrix1& A) const;
    
    //! launching all the simulation
    void RunAll(const string& input_file, const string&, const string&);
        
  };

  //! class to solve time-domain Maxwell equations with DG method
  template<>
  class HyperbolicProblem<TimeMaxwellEquation_Axi_DG> 
    : public TimeMaxwell_Axi<TimeMaxwellEquation_Axi_DG>
  {
    
  };

  //! class to solve time-domain Maxwell equations with DG method
  template<>
  class HyperbolicProblem<TimeMaxwellEquation_HcurlAxi> 
    : public TimeMaxwell_Axi<TimeMaxwellEquation_HcurlAxi>
  {
  };

}

#define MONTJOIE_FILE_TIME_MAXWELL_AXI_HXX
#endif
