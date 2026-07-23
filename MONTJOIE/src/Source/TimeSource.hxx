#ifndef MONTJOIE_FILE_TIME_SOURCE_HXX

namespace Montjoie
{
  
  //! class to evaluate Ricker
  class TimeRickerSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    Real_wp f0; //!< frequency
    
  public :
    TimeRickerSource();    
    TimeRickerSource(const Real_wp& freq);
    
    void SetFrequency(const Real_wp& freq);
    
    Real_wp Evaluate(const Real_wp&);
    Real_wp EvaluateDerivative(const Real_wp&);
    Real_wp EvaluateSecondDerivative(const Real_wp&);
    
  };


  //! class to evalute "derivative" of Ricker
  class DerivativeTimeRickerSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    TimeRickerSource f;
    
  public :
    DerivativeTimeRickerSource();    
    DerivativeTimeRickerSource(const Real_wp& freq);
    
    void SetFrequency(const Real_wp& freq);
    
    Real_wp Evaluate(const Real_wp&);
    Real_wp EvaluateDerivative(const Real_wp&);
    
  };


  //! class to evaluate Ricker
  class TimeModifiedRickerSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    Real_wp tc, f0;
    
  public :
    TimeModifiedRickerSource();
    TimeModifiedRickerSource(const Real_wp& freq, const Real_wp& t0);
    
    void Init(const Real_wp& freq, const Real_wp& t0);

    Real_wp Evaluate(const Real_wp&);
    Real_wp EvaluateDerivative(const Real_wp&);
    
  };
  

  //! class to evaluate Gaussian
  class TimeGaussianSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    Real_wp f0; //!< frequency
    Real_wp width_fct, t_begin, period;
    
  public:
    TimeGaussianSource();
    TimeGaussianSource(const Real_wp& freq, const Real_wp& t0, const Real_wp& width);   
    
    void Init(const Real_wp& freq, const Real_wp& t0, const Real_wp& width);
    
    Real_wp Evaluate(const Real_wp&);
    Real_wp EvaluateDerivative(const Real_wp&);
    
  };
  
  
  //! class to evaluate modulated cardinal sinus
  class TimeModulatedCardinalSinusSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    Real_wp T0, f0, Bw, Corr;
    
  public :
    TimeModulatedCardinalSinusSource();

    TimeModulatedCardinalSinusSource(const Real_wp& center_time, const Real_wp& frequency, const Real_wp& bandwidth,
				     const Real_wp& correlation);
    
    void Init(const Real_wp& center_time, const Real_wp& frequency, const Real_wp& bandwidth,
	      const Real_wp& correlation);
  
    Real_wp sinc(const Real_wp& t0);
    
    Real_wp Evaluate(const Real_wp&);
    
  };

  
  //! class to evaluate modulated cardinal cosinus
  class TimeModulatedCardinalCosinusSource : public TimeModulatedCardinalSinusSource
  {
  public :
    TimeModulatedCardinalCosinusSource();
    
    TimeModulatedCardinalCosinusSource(const Real_wp& center_time, const Real_wp& frequency,
				       const Real_wp& bandwidth, const Real_wp& correlation);
    
    Real_wp Evaluate(const Real_wp&);
    
  };

  //! class to evaluate modulated cardinal sinus
  class TimeRandomModulatedCardinalSource : public VirtualTimeSource<Real_wp>
  {
  protected :
    VectReal_wp values;
    Globatto<Real_wp> lob; //!< interpolation data
    Real_wp frequency; //!< frequency
    Real_wp bandwidth; //!< Bandwidth
    Real_wp center_time; //!< Center time of the function
    Real_wp correlation; //!< Correlation term
    Real_wp initial_time, final_time, dt;

  public:
    TimeRandomModulatedCardinalSource(const Real_wp& freqy, const Real_wp& bw,
				      const Real_wp& corr, const Real_wp& t0,
				      const Real_wp& t_init, const Real_wp& t_final,
				      const Real_wp deltat, bool take_cos);

    void Init(const Real_wp& freqy, const Real_wp& bw, const Real_wp& corr, const Real_wp& t0,
	      const Real_wp& t_init, const Real_wp& t_final, const Real_wp deltat, bool take_cos);
    
    Real_wp Evaluate(const Real_wp&);
    
  };


  //! class to evaluate sinus
  class TimeHarmonicSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    Real_wp f0, tend;
    
  public :
    TimeHarmonicSource();
    TimeHarmonicSource(const Real_wp& freq, const Real_wp& t_final);

    void Init(const Real_wp& freq, const Real_wp& t_final);
    
    Real_wp Evaluate(const Real_wp&);
    Real_wp EvaluateDerivative(const Real_wp&);

  };
  
  //! class to evaluate sinus modulated by a gaussian
  class TimeSinusGaussianSource : public VirtualTimeSource<Real_wp>
  {
  protected:
    Real_wp f0, toff, b;
    
  public :
    TimeSinusGaussianSource();
    TimeSinusGaussianSource(const Real_wp& freq, const Real_wp& t_offset, const Real_wp& alpha);

    void Init(const Real_wp& freq, const Real_wp& t_offset, const Real_wp& alpha);
    
    Real_wp Evaluate(const Real_wp&);
    Real_wp EvaluateDerivative(const Real_wp&);
    
  };
  
  //! source in time, whose values are read in a file
  class TimeFileSource : public VirtualTimeSource<Real_wp>
  {
  protected :
    VectReal_wp values;
    Globatto<Real_wp> lob; //!< interpolation data
    Real_wp frequency; //!< frequency
    Real_wp initial_time, final_time, dt; 
   
  public:
    TimeFileSource();

    TimeFileSource(const Real_wp& freq, const string& file_source_time,
		   const Real_wp& tinit, const Real_wp& t_final);
    
    void Init(const Real_wp& freq, const string& file_source_time,
	      const Real_wp& tinit, const Real_wp& t_final);
    
    Real_wp Evaluate(const Real_wp&);
    
  };


  //! class for sources in time 
  class TimeSourceHyperbolic
  {
  protected :
    //! order of interpolation
    int order_interpolation;
    //! interpolating points
    Globatto<Real_wp> lob;
    // interval of time
    Real_wp t0; //!< initial time
    Real_wp t1; //!< final time
    Real_wp epsilon_time; //!< threshold to distinguish two different times
    //! number of derivatives to be calculated
    int nb_deriv;
    //! stored values of source, and derivatives
    Matrix<Real_wp> ValueSource;
    
    VirtualTimeSource<Real_wp>* fsrc; //!< object describing temporal function
    Real_wp frequency; //!< frequency related to the source
    
  public :
    TimeSourceHyperbolic();
   
    ~TimeSourceHyperbolic();
    void Clear();
    
    // initialization of some variables, depending the time-problem var
    void Init(VirtualTimeSource<Real_wp>* f,
	      const Real_wp& freq, const Real_wp& eps, const Real_wp& t_init);

  protected:
    // computation of f and derivatives on interpolation points
    void GenerateValues();

  public:
    // computation of nth-derivative of f at any point
    void EvaluateDerivative(const Real_wp& t, int nb_deriv, Real_wp& pulse);
    
  };
    
} // namespace Montjoie

#define MONTJOIE_FILE_TIME_SOURCE_HXX
#endif
