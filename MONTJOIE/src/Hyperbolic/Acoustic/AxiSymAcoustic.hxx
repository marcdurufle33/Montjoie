#ifndef MONTJOIE_FILE_AXISYM_ACOUSTIC_HXX

namespace Montjoie
{
  
  //! Axisymmetric acoustic equation
  class AcousticEquationAxi
  {
  public :
    typedef LaplaceEquationAxi TypeEquationStationary; //!< stationary equation
    
  };


  //! class to solve time-domain acoustic equations in axisymmetric domain
  template<class TypeEquation>
  class TimeAcoustic_Axi : public VarInstationary<TypeEquation>
  {
  public:    
    void SetInputData(const string& keyword, const Vector<string>& param);

  };
  
  //! class to solve acoustic equations with finite element method
  template<>
  class HyperbolicProblem<AcousticEquationAxi> 
    : public TimeAcoustic_Axi<AcousticEquationAxi>
  {
  };

}

#define MONTJOIE_FILE_AXISYM_ACOUSTIC_HXX
#endif
