#ifndef MONTJOIE_FILE_TIME_AERO_ACOUSTIC_HXX

namespace Montjoie
{
  
  //! base class for time aero-acoustic equation
  template<class Dimension>
  class TimeLinearizedEulerEquation
  {
  public :
    typedef StationaryLinearizedEulerEquation<Dimension> TypeEquationStationary;
    //!< stationary equation
    
  };
  
  //! base class for time aero-acoustic equation
  template<class Dimension>
  class TimeAeroAcousticEquation
  {
  public :
    typedef AeroStationaryEquation<Dimension> TypeEquationStationary; //!< stationary equation
    
  };
  
  //! class to solve aeroacoustic equation in time-domain
  template<class TypeEquation>
  class TimeAeroAcoustic_Base : public VarInstationary<TypeEquation>
  {
  };

  //! class used to solve aeroacoustic equation in time-domain and uniform flow with DG method
  template<class Dimension>
  class HyperbolicProblem<TimeLinearizedEulerEquation<Dimension> >
    : public TimeAeroAcoustic_Base<TimeLinearizedEulerEquation<Dimension> >
  {
  public:
    void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp&, Real_wp&);
  };

  //! class used to solve aeroacoustic equation in time-domain with DG method
  template<class Dimension>
  class HyperbolicProblem<TimeAeroAcousticEquation<Dimension> >
    : public TimeAeroAcoustic_Base<TimeAeroAcousticEquation<Dimension> >
  {
  public:
    void GetMassDampingCoefficient(int i, int j, int num, int ref, Real_wp&, Real_wp&);
  };
  
}

#define MONTJOIE_FILE_TIME_AERO_ACOUSTIC_HXX
#endif
