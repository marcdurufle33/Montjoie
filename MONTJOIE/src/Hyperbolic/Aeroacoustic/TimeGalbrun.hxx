#ifndef MONTJOIE_FILE_TIME_GALBRUN_HXX

namespace Montjoie
{
  
  //! base class for time aero-acoustic equation
  template<class Dimension>
  class TimeGalbrunEquation : public AcousticEquationDG<Dimension>
  {
  public :
    typedef GalbrunStationaryEquationDG<Dimension> TypeEquationStationary; //!< stationary equation
    
  };
  
  
  //! class to solve Galbrun equation in time-domain
  template<class TypeEquation>
  class TimeGalbrun_Base : public VarInstationary<TypeEquation>
  {    
  };

    
  //! class used to solve aeroacoustic equation in time-domain with DG method
  template<class Dimension>
  class HyperbolicProblem<TimeGalbrunEquation<Dimension> >
    : public TimeGalbrun_Base<TimeGalbrunEquation<Dimension> >
  {
  public:
    void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp&, Real_wp&);
    
  };
  
}

#define MONTJOIE_FILE_TIME_GALBRUN_HXX
#endif
