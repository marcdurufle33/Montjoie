#ifndef MONTJOIE_FILE_PHYSICAL_CONSTANTS_HXX

namespace Montjoie
{

  //! Class containing physical constants
  /*!
    Physical constants for acoustical computation
  */
  class Physics
  {
  public:
    Real_wp c; // celerity
    Real_wp rho; // density
    Real_wp gamma; // 
    Real_wp mu; // 
    Real_wp kappa; // 
    Real_wp Cp; // 
    Real_wp temp; // temperature

    Physics(const Real_wp temperature)
    {
      // formulas from Chaigne and Kergomard 2016 
      Real_wp T0 = 273.15;// absolute temperature
      Real_wp Temp = temperature + T0;
      temp = temperature;
      c = 331.45 * sqrt(Temp/T0);
      rho = 1.2929*T0/Temp;
      gamma = 1.402;
      mu = 1.708e-5 * (1.0 + 0.0029 * temperature);
      Cp = 240;
      kappa = 5.77 * 1e-3 * (1 + 0.0033 * temperature);
    };

    void Display()
    {
      cout << "Physical Parameters :" <<endl;
      DISP(temp);
      DISP(c);
      DISP(rho);
      DISP(gamma);
      DISP(mu);
      DISP(kappa);
      DISP(Cp);
    };

  };
}

#define MONTJOIE_FILE_PHYSICAL_CONSTANTS_HXX
#endif