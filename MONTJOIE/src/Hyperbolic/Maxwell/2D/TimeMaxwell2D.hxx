#ifndef MONTJOIE_FILE_TIME_MAXWELL_2D_HXX

namespace Montjoie
{

  //! class to specify 2-D time Maxwell equation for edge finite element
  class TimeMaxwellEquation_2D
  {
  public :
    typedef StaticMaxwellEquation_2D TypeEquationStationary; //!< stationary equation
    
  };
  
  //! class to specify 2-D maxwell equation and discontinuous galerkin
  class TimeMaxwellEquation_2D_DG
  {
  public :
    typedef StaticMaxwellEquation_2D_DG TypeEquationStationary; //!< stationary equation
    
  };
  
  
  //! class to solve time-domain Maxwell equations in 2-D
  template<class TypeEquation>
  class TimeMaxwell_2D: public VarInstationary<TypeEquation>
  {
  public:
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
    //!< stationary equation
    

    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
				    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = true);

    void EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp&, int level,
				    const VectReal_wp& X, VectReal_wp& Y,
                                    bool invert_mass = true, bool source = true);
    
  };


  //! class to solve time-domain Maxwell equations with edge finite elements in 3-D
  template<>
  class HyperbolicProblem<TimeMaxwellEquation_2D> 
    : public TimeMaxwell_2D<TimeMaxwellEquation_2D>
  {
  };

}

#define MONTJOIE_FILE_TIME_MAXWELL_2D_HXX
#endif
