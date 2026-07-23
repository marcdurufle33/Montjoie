#ifndef MONTJOIE_FILE_TIME_SOURCE_CXX

namespace Montjoie
{

  /********************
   * TimeRickerSource *
   ********************/
  
  
  //! default constructor
  TimeRickerSource::TimeRickerSource()
  {
    f0 = 0.0; 
  }
  
  
  //! constructor with the frequency
  TimeRickerSource::TimeRickerSource(const Real_wp& freq)
  {
    f0 = freq;
  }
  
  
  //! modifies the frequency
  void TimeRickerSource::SetFrequency(const Real_wp& freq)
  {
    f0 = freq;
  }


  //! Returns primitive of a Ricker at time t0 with central frequency f0
  Real_wp TimeRickerSource::Evaluate(const Real_wp& t0)
  {
    Real_wp T0 = 1.0/f0;
    if ((t0 <= -0.9*T0) || (t0 >= 2.9*T0))
      return Real_wp(0);
    
    Real_wp pulse = 0;

    Real_wp tmp = f0*t0-1.0;
    Real_wp pi2 = pi_wp*pi_wp;
    Real_wp f0tm1_2 = tmp*tmp;
    Real_wp gaussian_term = exp(-f0tm1_2*pi2); 
      
    // primitive of a Ricker
    // -(f0 t-1) exp(-(f0 t -1)^2 pi^2 ) / f0
    pulse = -tmp*gaussian_term / f0;
    return pulse;
  }
  
  
  //! Derivative of the source
  Real_wp TimeRickerSource::EvaluateDerivative(const Real_wp& t0)
  {
    Real_wp T0 = 1.0/f0;
    if ((t0 <= -0.9*T0) || (t0 >= 2.9*T0))
      return Real_wp(0);

    // actual expression of the Ricker
    Real_wp der_pulse = 0;
    Real_wp tmp = f0*t0-1.0;
    Real_wp f0tm1_2 = square(tmp*pi_wp);
    der_pulse = (2.0*f0tm1_2-1.0)*exp(-f0tm1_2);
    return der_pulse;
  }


  //! second-derivative of the source
  Real_wp TimeRickerSource::EvaluateSecondDerivative(const Real_wp& t0)
  {
    Real_wp T0 = 1.0/f0;
    if ((t0 <= -0.9*T0) || (t0 >= 2.9*T0))
      return Real_wp(0);

    // derivative of the Ricker
    Real_wp der_pulse = 0;
    Real_wp tmp = f0*t0-1.0;
    Real_wp f0tm1_2 = square(tmp*pi_wp);
    der_pulse = 2*pi_wp*pi_wp*f0*tmp*(3.0 - 2.0*f0tm1_2)*exp(-f0tm1_2);
    return der_pulse;
  }


  /******************************
   * DerivativeTimeRickerSource *
   ******************************/
  

  //! default constructor
  DerivativeTimeRickerSource::DerivativeTimeRickerSource()
    : f()
  {
  }


  //! constructor with a given frequency
  DerivativeTimeRickerSource::DerivativeTimeRickerSource(const Real_wp& freq)
    : f(freq)
  {
  }
  

  //! sets the frequency of the Ricket
  void DerivativeTimeRickerSource::SetFrequency(const Real_wp& freq)
  {
    f.SetFrequency(freq);
  }


  //! evaluates the Ricker
  Real_wp DerivativeTimeRickerSource::Evaluate(const Real_wp& t)
  {
    return f.EvaluateDerivative(t);
  }


  //! evaluate the derivative of the Ricker
  Real_wp DerivativeTimeRickerSource::EvaluateDerivative(const Real_wp& t)
  {
    return f.EvaluateSecondDerivative(t);
  }
  

  /****************************
   * TimeModifiedRickerSource *
   ****************************/

  
  //! default constructor
  TimeModifiedRickerSource::TimeModifiedRickerSource()
  {
    tc = 0.0;
    f0 = 0; 
  }
  
  
  //! constructor with the frequency and center time
  TimeModifiedRickerSource::TimeModifiedRickerSource(const Real_wp& freq, const Real_wp& t0)
  {
    f0 = freq;
    tc = t0;
  }
  
  
  //! Sets the frequency and center time of the source
  void TimeModifiedRickerSource::Init(const Real_wp& freq, const Real_wp& t0)
  {
    f0 = freq; tc = t0; 
  }


  //! Returns the value of a Ricker at time t0 with central Fourier frequency f0 and center time tc
  Real_wp TimeModifiedRickerSource::Evaluate(const Real_wp& t0)
  {
    Real_wp pulse = 0;

    Real_wp tmp = f0*(t0-tc);
    Real_wp f0tm1_2 = 2*square(tmp*pi_wp);
    Real_wp gaussian_term = exp(-f0tm1_2); 
      
    // primitive of a Rickert
    // -(t-tc) exp(- 2 f0^2 (t -tc)^2 pi^2 )
    pulse = -(t0-tc)*gaussian_term;
        
    return pulse;
  }


  //! Returns the derivative of the Ricker
  Real_wp TimeModifiedRickerSource::EvaluateDerivative(const Real_wp& t0)
  {
    Real_wp der_pulse = 0;

    Real_wp tmp = f0*(t0-tc);
    Real_wp f0tm1_2 = 2*square(tmp*pi_wp);
    Real_wp gaussian_term = exp(-f0tm1_2); 
      
    der_pulse = (2*f0tm1_2-1)*gaussian_term;
    
    return der_pulse;    
  }


  /**********************
   * TimeGaussianSource *
   **********************/
  

  //! default constructor
  TimeGaussianSource::TimeGaussianSource()
  {
    f0 = 0.0; period = 0.0;
    width_fct = 0.0; t_begin = 0.0;
  }
    

  //! constructor with the frequency, beginning time and width of the plateau
  TimeGaussianSource::TimeGaussianSource(const Real_wp& freq, const Real_wp& t0, const Real_wp& width)
  { 
    f0 = freq; 
    period = 1.0/freq; 
    width_fct = width; 
    t_begin = t0; 
  }
  
  
  //! Inits the frequency, beginning time and width of the plateau
  void TimeGaussianSource::Init(const Real_wp& freq, const Real_wp& t0, const Real_wp& width)
  { 
    f0 = freq; 
    period = 1.0/freq; 
    width_fct = width; 
    t_begin = t0; 
  }
    

  //! Returns the value of a gaussian at time t0 with frequency f0
  Real_wp TimeGaussianSource::Evaluate(const Real_wp& t)
  {
    // we add a plateau of value 1 and of width width_fct
    Real_wp pulse = 0, t0(t);
    if (t < t_begin+period)
      t0 = t-t_begin;
    else if (t >= t_begin+period+width_fct)
      t0 = t-t_begin-width_fct;
    else
      return Real_wp(1);
    
    Real_wp tmp = f0*t0-1.0;
    Real_wp f0tm1_2 = square(pi_wp*tmp);
    pulse = exp(-f0tm1_2); 
    
    return pulse;
  }

  
  //! Returns the derivative of the gaussian
  Real_wp TimeGaussianSource::EvaluateDerivative(const Real_wp& t)
  {
    Real_wp der_pulse = 0, t0(t);
    if (t < t_begin+period)
      t0 = t-t_begin;
    else if (t >= t_begin+period+width_fct)
      t0 = t-t_begin-width_fct;
    else
      return Real_wp(0);
    
    Real_wp tmp = pi_wp*(f0*t0-1.0);
    Real_wp f0tm1_2 = tmp*tmp;
    der_pulse = -2.0*pi_wp*tmp*f0*exp(-f0tm1_2); 
    
    return der_pulse;
  }
  

  /************************************
   * TimeModulatedCardinalSinusSource *
   ************************************/
  
  
  //! default constructor
  TimeModulatedCardinalSinusSource::TimeModulatedCardinalSinusSource()
  { 
    T0 = 0.0; f0 = 0.0;
    Bw = 0.0; Corr = 0.0; 
  }
  
  
  //! constructor with the center time, frequency, bandwith and correlation
  TimeModulatedCardinalSinusSource
  ::TimeModulatedCardinalSinusSource(const Real_wp& center_time, const Real_wp& frequency, const Real_wp& bandwidth,
				     const Real_wp& correlation)
  { 
    T0 = center_time; f0 = frequency; 
    Bw = bandwidth; Corr = correlation; 
  }
  
  
  //! Sets the center time, frequency, bandwith and correlation
  void TimeModulatedCardinalSinusSource
  ::Init(const Real_wp& center_time, const Real_wp& frequency, const Real_wp& bandwidth,
	 const Real_wp& correlation)
  { 
    T0 = center_time; f0 = frequency;
    Bw = bandwidth; Corr = correlation; 
  }


  //! Evaluation of cardinal sinus
  Real_wp TimeModulatedCardinalSinusSource::sinc(const Real_wp& t0)
  {
    if (t0 == 0.0) return 1.0;
    return sin(t0) / t0;
  }
  

  //! Evaluation of modified cardinal sinus using a gaussian function
  Real_wp TimeModulatedCardinalSinusSource::Evaluate(const Real_wp& t0)
  {
    Real_wp pulse = sinc(Bw*(t0-T0));
    Real_wp tmp = sin(2. * pi_wp * f0 * (t0-T0));
    Real_wp tmp2 = exp(-(t0-T0)*(t0-T0)/(2.*Corr*Corr));
    pulse = pulse * tmp * tmp2;
    
    return pulse;
  }


  //! default constructor
  TimeModulatedCardinalCosinusSource::TimeModulatedCardinalCosinusSource()
  {
  }
  
  
  //! constructor with the center time, frequency, bandwith and correlation  
  TimeModulatedCardinalCosinusSource
  ::TimeModulatedCardinalCosinusSource(const Real_wp& center_time, const Real_wp& frequency,
				       const Real_wp& bandwidth, const Real_wp& correlation)
    : TimeModulatedCardinalSinusSource(center_time, frequency, bandwidth, correlation) 
  {
  }
  
  
  //! Evaluation of modified cardinal cosinus using a gaussian function
  Real_wp TimeModulatedCardinalCosinusSource::Evaluate(const Real_wp& t0)
  {
    Real_wp pulse = sinc(Bw*(t0-T0));
    Real_wp tmp = cos(2. * pi_wp * f0 * (t0-T0));
    Real_wp tmp2 = exp(-(t0-T0)*(t0-T0)/(2.*Corr*Corr));
    pulse = pulse * tmp * tmp2;
    
    return pulse;
  }
  

  /*************************************
   * TimeRandomModulatedCardinalSource *
   *************************************/
  
  
  //! constructor with parameters
  TimeRandomModulatedCardinalSource
  ::TimeRandomModulatedCardinalSource(const Real_wp& freq, const Real_wp& bw,
				      const Real_wp& corr, const Real_wp& t0,
				      const Real_wp& t_init, const Real_wp& t_final,
				      const Real_wp deltat, bool take_cos)
  {
    Init(freq, bw, corr, t0, t_init, t_final, deltat, take_cos);
  }
  

  //! kind of constructor
  void TimeRandomModulatedCardinalSource
  ::Init(const Real_wp& freq, const Real_wp& bw, const Real_wp& corr, const Real_wp& t0,
	 const Real_wp& t_init, const Real_wp& t_final, const Real_wp deltat, bool take_cos)
  {
    frequency = freq;
    bandwidth = bw;
    correlation = corr;
    center_time = t0;

    // time steps
    initial_time = t_init;
    final_time = t_final;
    dt = deltat;
    
    VirtualTimeSource<Real_wp>* fsrc;
    if (take_cos)
      fsrc = new TimeModulatedCardinalCosinusSource(t0, freq, bw, corr);
    else
      fsrc = new TimeModulatedCardinalSinusSource(t0, freq, bw, corr);

#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#endif
        
    int Ntime = toInteger((final_time - initial_time) / dt)+1;
    // Be sure that Ntime can be divided by 2
    if ((Ntime % 2) !=0) Ntime++;
    // Here we're sure that Ntime is even

    values.Resize(Ntime);
    for (int counter = 0; counter < Ntime; counter++)
      values(counter) = fsrc->Evaluate(Real_wp(counter)*dt);
    
    delete fsrc;
    
    // Allocation of a complex vector whose size is size of the real vector
    Vector<Complex_wp> cvalues;
    cvalues.Resize(Ntime);

    // Define and initialize the Fourier transform
    FftInterface<Complex_wp> Fourier;
    Fourier.Init(Ntime);

    for(int counter=0; counter < Ntime; counter++)
      {
        cvalues(counter) = values(counter);
      }

    // Use Forward Fourier Transform (overwriting in cvalues)
    Fourier.ApplyForward(cvalues);

    // Force recall of seed and arbitrary generation of 3 random numbers
    srand(Seed);
    for(int k=0; k<3; k++)
      {
        rand();
      }

    for(int counter=1; counter < Ntime/2; counter++)
      {
        Real_wp random_amplitude = Real_wp(rand())/Real_wp(RAND_MAX);
        Real_wp random_phase = 2 * pi_wp * Real_wp(rand())/Real_wp(RAND_MAX);
#ifdef SELDON_WITH_MPI
        if (rank_proc == 0)
          {
            for (int processor_to_send = 1;
                 processor_to_send < nb_proc;
                 processor_to_send++)
              {
                MPI_Send(&random_amplitude, 1,
                         MPI_DOUBLE, processor_to_send,
                         667, MPI_COMM_WORLD);
              }
          }
        else
          {
            MPI_Status status;
            MPI_Recv(&random_amplitude, 1, MPI_DOUBLE,
                     0, 667, MPI_COMM_WORLD, &status);
          }
        if (rank_proc == 0)
          {
            for (int processor_to_send = 1;
                 processor_to_send < nb_proc;
                 processor_to_send++)
              {
                MPI_Send(&random_phase, 1,
                         MPI_DOUBLE, processor_to_send,
                         668, MPI_COMM_WORLD);
              }
          }
        else
          {
            MPI_Status status;
            MPI_Recv(&random_phase, 1, MPI_DOUBLE,
                     0, 668, MPI_COMM_WORLD, &status);
          }
#endif
        cvalues(counter) *= random_amplitude * exp(Iwp*random_phase);
        cvalues(Ntime-counter) *= random_amplitude * exp(-Iwp*random_phase);
      }

    Real_wp random_amplitude_0 = Real_wp(rand())/Real_wp(RAND_MAX);
#ifdef SELDON_WITH_MPI
    if (rank_proc == 0)
      {
        for (int processor_to_send = 1;
             processor_to_send < nb_proc;
             processor_to_send++)
          {
            MPI_Send(&random_amplitude_0, 1,
                     MPI_DOUBLE, processor_to_send,
                     667, MPI_COMM_WORLD);
          }
      }
    else
      {
        MPI_Status status;
        MPI_Recv(&random_amplitude_0, 1, MPI_DOUBLE,
                 0, 667, MPI_COMM_WORLD, &status);
      }
#endif

    cvalues(0)*= 2.0*random_amplitude_0-1.0;

    Real_wp random_amplitude_half = Real_wp(rand())/Real_wp(RAND_MAX);
#ifdef SELDON_WITH_MPI
    if (rank_proc == 0)
      {
        for (int processor_to_send = 1;
             processor_to_send < nb_proc;
             processor_to_send++)
          {
            MPI_Send(&random_amplitude_half, 1,
                     MPI_DOUBLE, processor_to_send,
                     667, MPI_COMM_WORLD);
          }
      }
    else
      {
        MPI_Status status;
        MPI_Recv(&random_amplitude_half, 1, MPI_DOUBLE,
                 0, 667, MPI_COMM_WORLD, &status);
      }
#endif

    cvalues(Ntime/2)*= 2.0*random_amplitude_half-1.0;

    // Use Inverse Fourier Transform (overwriting in cvalues)
    Fourier.ApplyInverse(cvalues);

    // Get back the real value
    for(int counter=0; counter < Ntime; counter++)
      {
        values(counter) = real(cvalues(counter));
      }

    // Modify the value vector to have the mean-value equal to 0 and the
    // variance equal to 1.
    Real_wp average_values=0., variance_values=0.;
    for(int counter=0; counter < Ntime; counter++)
      {
        average_values += values(counter);
      }
    average_values /= Real_wp(Ntime);
    for(int counter=0; counter < Ntime; counter++)
      {
        values(counter) -= average_values;
      }

    for(int counter=0; counter < Ntime; counter++)
      {
        variance_values += square(values(counter));
      }
    variance_values /= Real_wp(Ntime);
    for(int counter=0; counter < Ntime; counter++)
      {
        values(counter) /= sqrt(variance_values);
      }

    // cubic interpolation
    Real_wp one(1);
    VectReal_wp Points(4); Points.Fill(); Mlt(one/3, Points);
    lob.AffectPoints(Points);
    
  }

  
  //! evaluation of the modulated cardinal sinus
  Real_wp TimeRandomModulatedCardinalSource::Evaluate(const Real_wp& t0)
  {
    if ((t0 < initial_time)||(t0 > final_time))
      return 0.0;
    
    // finding interval where t is
    int n0 = toInteger(floor((t0-initial_time)/(3.0*dt))); n0 = max(0, n0);
    int i = 3*n0, j = i + 3;
    if (j >= values.GetM())
      {
	j = values.GetM() - 1; i = j-3;
      }
    
    Real_wp lambda = (t0 - dt*i-initial_time)/(3.0*dt);
    // cubic interpolation
    Real_wp f = 0;
    for (int k = 0; k < 4; k++)
      f += lob.EvaluatePhi(k, lambda)*values(i+k);
    
    return f;
  }  


  /**********************
   * TimeHarmonicSource *
   **********************/


  //! default constructor
  TimeHarmonicSource::TimeHarmonicSource()
  {
    f0 = 0.0; 
    tend = 0.0; 
  }


  //! constructor with the frequency and final time
  TimeHarmonicSource::TimeHarmonicSource(const Real_wp& freq, const Real_wp& t_final) 
  {
    f0 = freq;
    tend = t_final;
  }
  
  
  //! Sets the frequency and final time
  void TimeHarmonicSource::Init(const Real_wp& freq, const Real_wp& t_final) 
  {
    f0 = freq; 
    tend = t_final; 
  }
  
  
  //! sinus with a cut-off function so that it is infinitely smooth
  Real_wp TimeHarmonicSource::Evaluate(const Real_wp& t0)
  {
    Real_wp pulse = 0;
    if ((t0 <= 0)||(t0 >= tend))
      return pulse;
    
    // Harmonic source
    // sin(2 \pi f0 t) * exp(-1/(pi t^2) ) * exp(-1 / (pi (t-T)^2 ) )
    Real_wp omega = 2.0*pi_wp*f0; Real_wp phase = omega*t0;
    Real_wp f0t_2 = f0*f0*t0*t0;
    Real_wp f0_tmT_2 = f0*f0*(t0-tend)*(t0-tend);
    Real_wp exp_minus = exp(-1.0/(pi_wp*f0t_2)) , exp_plus = exp(-1.0/(pi_wp*f0_tmT_2)); 
    pulse = exp_minus*sin(phase)*exp_plus;
    return pulse;
  }

  
  //! derivative of HarmonicTime
  Real_wp TimeHarmonicSource::EvaluateDerivative(const Real_wp& t0)
  {
    Real_wp der_pulse = 0;
    if ((t0 <= 0)||(t0 >= tend))
      return der_pulse;
    
    Real_wp omega = 2.0*pi_wp*f0; Real_wp phase = omega*t0;
    Real_wp f0t_2 = square(f0*t0);
    Real_wp f0_tmT_2 = square(f0*(t0-tend));
    Real_wp exp_minus = exp(-1.0/(pi_wp*f0t_2)) , exp_plus = exp(-1.0/(pi_wp*f0_tmT_2)); 
    der_pulse = 2.0*exp_minus*exp_plus
      *( pi_wp*f0*cos(phase) 
         + sin(phase)*(1.0/(pi_wp*f0t_2*t0) + 1.0/(pi_wp*f0_tmT_2*(t0-tend))));
    
    return der_pulse;
  }

  
  /***************************
   * TimeSinusGaussianSource *
   ***************************/
  
  
  //! default constructor
  TimeSinusGaussianSource::TimeSinusGaussianSource()
  {
    f0 = 0.0;
    toff = 0.0;
    b = 0.0;
  }
  
  
  //! constructor with the frequency, center time and gaussian parameter
  TimeSinusGaussianSource
  ::TimeSinusGaussianSource(const Real_wp& freq, const Real_wp& t_offset, const Real_wp& alpha) 
  {
    f0 = freq;
    toff = t_offset;
    b = alpha; 
  }
  
  
  //! Sets the frequency, center time and gaussian parameter
  void TimeSinusGaussianSource::Init(const Real_wp& freq, const Real_wp& t_offset, const Real_wp& alpha) 
  { 
    f0 = freq; 
    toff = t_offset; 
    b = alpha; 
  }


  //! evaluates a sinus modulated by a gaussian
  Real_wp TimeSinusGaussianSource::Evaluate(const Real_wp& t0)
  {
    Real_wp pulse = 0;
    Real_wp t0_toff = t0-toff, t0_toff_2 = t0_toff*t0_toff;

    // Harmonic source
    // sin(2 \pi f0 t) * exp(-b*(t-toff)**2)
    Real_wp omega = 2.0*pi_wp*f0; Real_wp phase = omega*t0;
    Real_wp gaussian = exp(-b*t0_toff_2);
    pulse = sin(phase)*gaussian;
    
    return pulse;
  }

  
  //! derivative of the modulated sinus
  Real_wp TimeSinusGaussianSource::EvaluateDerivative(const Real_wp& t0)
  {
    Real_wp der_pulse = 0;
    Real_wp t0_toff = t0-toff, t0_toff_2 = t0_toff*t0_toff;

    Real_wp omega = 2.0*pi_wp*f0; Real_wp phase = omega*t0;
    Real_wp gaussian = exp(-b*t0_toff_2);
    der_pulse = gaussian*(-2*b*t0_toff*sin(phase) + omega*cos(phase));
    
    return der_pulse;
  }

  
  /******************
   * TimeFileSource *
   ******************/

  
  //! default constructor
  TimeFileSource::TimeFileSource()
  {
    frequency = 0.0;
    initial_time = 0.0;
    final_time = 0.0;
    dt = 0.0;
  }
  

  //! constructor with the file containing values
  TimeFileSource::TimeFileSource(const Real_wp& freq, const string& file_source_time,
				 const Real_wp& tinit, const Real_wp& t_final)
  {
    frequency = 0.0;
    initial_time = 0.0;
    final_time = 0.0;
    dt = 0.0;

    Init(freq, file_source_time, tinit, t_final);
  }


  //! initialization of the computation of a time source with values contained in a file
  void TimeFileSource::Init(const Real_wp& freq, const string& file_source_time,
			    const Real_wp& tinit, const Real_wp& t_final)
  {
    frequency = freq;
    values.ReadText(file_source_time);
    // cubic interpolation
    Real_wp one(1);
    VectReal_wp Points(4); Points.Fill(); Mlt(one/3, Points);
    lob.AffectPoints(Points);
    
    // time steps
    initial_time = tinit;
    final_time = t_final;
    dt = (final_time-initial_time)/(values.GetM()-1);
    
    // DISP(initial_time); DISP(final_time); DISP(dt); DISP(nb_values); DISP(values(50));
  }
  
  
  //! evaluation of the source at time t
  /*!
    \param[in] t fixed time where the source has to be computed
    \returns value of the source
   */
  Real_wp TimeFileSource::Evaluate(const Real_wp& t)
  {
    if ((t < initial_time)||(t>final_time))
      return 0.0;
    
    // finding interval where t is
    int n0 = toInteger(floor((t-initial_time)/(3.0*dt))); n0 = max(0, n0);
    int i = 3*n0, j = i + 3;
    if (j >= values.GetM())
      {
	j = values.GetM() - 1; i = j-3;
      }
    
    Real_wp lambda = (t - dt*i-initial_time)/(3.0*dt);
    // cubic interpolation
    Real_wp f = 0;
    for (int k = 0; k < 4; k++)
      f += lob.EvaluatePhi(k, lambda)*values(i+k);
    
    return f;
  }  
  
  
  /************************
   * TimeSourceHyperbolic *
   ************************/
  
  
  //! default constructor
  TimeSourceHyperbolic::TimeSourceHyperbolic()
  {
    fsrc = NULL;
    
    // use of an interpolation method to evaluate the derivatives
    order_interpolation = 24;
    lob.ConstructQuadrature(order_interpolation, lob.QUADRATURE_LOBATTO);
    lob.ComputeGradPhi(1e3*epsilon_machine);
    t0 = 0; t1 = -1;
    frequency = 1.0;
    epsilon_time = 1e-10;
    nb_deriv = 10;
    ValueSource.Reallocate(nb_deriv, order_interpolation+1);
    ValueSource.Fill(0.0);
  }
  
  
  //! destructor
  TimeSourceHyperbolic::~TimeSourceHyperbolic()
  {
    Clear();
  }


  //! releases memory used by the object
  void TimeSourceHyperbolic::Clear()
  {
    if (fsrc != NULL)
      {
	delete fsrc;
	fsrc = NULL;
      }
  }
  
  
  //! initialization before computation of source and derivatives
  void TimeSourceHyperbolic::Init(VirtualTimeSource<Real_wp>* f,
				  const Real_wp& freq, const Real_wp& eps, const Real_wp& t_init) 
  {
    fsrc = f;
    frequency = freq;
    epsilon_time = eps;
    t0 = t_init; t1 = t0 + 0.5/frequency;
    ValueSource.Reallocate(nb_deriv, order_interpolation+1);
    GenerateValues();
  }
  
  
  //! computation of source and derivatives on interpolation points (Gauss-Lobatto points)
  void TimeSourceHyperbolic::GenerateValues()
  {
    if (fsrc == NULL)
      {
	ValueSource.Fill(0);
	return;
      }
    
    // evaluation of the function on interpolation points
    for (int i = 0; i <= order_interpolation; i++)
      {
	Real_wp t = (1.0-lob.Points(i))*t0 + lob.Points(i)*t1;
	ValueSource(0, i) = fsrc->Evaluate(t);
      }
    
    // generation of derivatives
    for (int i = 1; i < nb_deriv; i++)
      {
	for (int j = 0; j <= order_interpolation; j++)
	  {
	    ValueSource(i, j) = 0.0;
	    for (int k = 0; k <= order_interpolation; k++)
	      ValueSource(i, j) += lob.GradPhi(k,j)*ValueSource(i-1, k);
	    
	    ValueSource(i,j) *= 1.0/(t1-t0);
	  }
      }
    
    // DISP(t0); DISP(t1); DISP(ValueSource);
  }
  
  
  //! evaluation of derivative of source by interpolation
  void TimeSourceHyperbolic
  ::EvaluateDerivative(const Real_wp& t, int nb_deriv, Real_wp& pulse)
  {
    if ((t > t1+epsilon_time) && ( t1 > t0) )
      {
	// new interval
	while (t > t1+epsilon_time)
	  {
	    Real_wp tnext = t1 + (t1-t0);
	    t0 = t1; t1 = tnext;
	  }
	
	GenerateValues();
      }
   
    // computation of the source by interpolation
    pulse = 0.0;
    for (int i = 0; i <= order_interpolation; i++)
      {
	Real_wp lambda = (t-t0)/(t1-t0);
	Real_wp coef = lob.EvaluatePhi(i, lambda);
	pulse += coef*ValueSource(nb_deriv,i);
      }
    
    // DISP(t); DISP(nb_deriv); DISP(pulse);
  }  
  
}

#define MONTJOIE_FILE_TIME_SOURCE_CXX
#endif
