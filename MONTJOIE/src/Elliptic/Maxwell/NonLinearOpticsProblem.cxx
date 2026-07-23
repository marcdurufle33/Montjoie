#ifndef MONTJOIE_FILE_NON_LINEAR_OPTICS_PROBLEM_CXX

namespace Montjoie
{
  //! default constructor
  RamanEffect::RamanEffect()
  {
    tau1 = 0.0; tau2 = 0.0; Tmax = 0.0;
    invTau1 = 0.0; invTau2 = 0.0;
    alpha = 0.0; coef = 0.0;
    level_romberg = 10;
    deltat = 0; length_interval = 0;
    num_current = -1; t_current = 0;
    level_global = 0;    
  }
  
  
  //! returns parameter alpha
  Real_wp RamanEffect::GetAlpha() const
  {
    return alpha;
  }
  
  
  //! sets tau1, tau2, Tmax and alpha
  void RamanEffect::Init(const Real_wp& t1, const Real_wp& t2,
			 const Real_wp& tmax, const Real_wp& a)
  {
    tau1 = t1; invTau1 = 1.0/t1;
    tau2 = t2; invTau2 = 1.0/t2;
    Tmax = tmax;
    alpha = a;
    coef = (t1*t1 + t2*t2) / (t1*t2*t2);
  }
  
  
  //! returns function h_R(t)
  Real_wp RamanEffect::EvaluateH(const Real_wp& t) const
  {
    if (t < 0)
      return Real_wp(0);
    
    Real_wp res = coef*exp(-t*invTau2)*sin(t*invTau1);
    
    // smooth function
    // Real_wp res = coef*exp(-0.5*square((t-8.0*tau2)/tau2));
    return res;
  }
  
  
  //! fills coefficients of the associated PDE
  /*!
    P_R satisfies the following PDE :
    d^2 P_R/dt^2 + gamma_r dP_R/dt + beta_r P_R = alpha_R |E|^2
    and P_R is present in the expression of the displacement :
    D = \eps_inf E + gamma (1-alpha) E |E|^2 + gamma alpha E P_R
   */
  void RamanEffect::GetCoefPDE(Real_wp& alpha_r, Real_wp& beta_r, Real_wp& gamma_r) const
  {
    alpha_r = (tau1*tau1 + tau2*tau2) * (invTau1*invTau1*invTau2*invTau2);
    beta_r = invTau1*invTau1 + invTau2*invTau2;
    gamma_r = 2.0*invTau2;
  }
  
  
  //! allocates arrays to store |E|^2 at different times
  void RamanEffect::Construct(int N, const Real_wp& dt)
  {
    if (alpha == Real_wp(0))
      return;
    
    // counting the number of intervals (power of 2)
    deltat = dt;
    length_interval = dt;
    level_global = 0;
    int nb_intervals = 1;
    while (length_interval < Tmax)
      {
	level_global++;
	length_interval *= 2;
	nb_intervals *= 2;
      }
    
    //DISP(level_global); DISP(length_interval);
    
    // allocating the arrays
    //DISP(nb_intervals);
    moduleE.Reallocate(nb_intervals+1);
    for (int i = 0; i <= nb_intervals; i++)
      {
	moduleE(i).Reallocate(N);
	moduleE(i).Fill(0);
      }
    
    num_current = -1;
  }

  
  //! stores |E|^2 and increments current time
  template<class T>
  void RamanEffect::StoreNewField(const Real_wp& t, const Vector<T>& E)
  {
    if (alpha == Real_wp(0))
      return;

    num_current = (num_current+1)%moduleE.GetM();
    t_current = t;
    for (int i = 0; i < moduleE(num_current).GetM(); i++)
      moduleE(num_current)(i) = absSquare(E(i));
  }
  
  
  //! returns the number of stored iterates
  int RamanEffect::GetNumberOfIterates() const
  {
    return moduleE.GetM();
  }
  
  
  //! returns iterate k
  VectReal_wp& RamanEffect::GetIterate(int k)
  {
    return moduleE(k);
  }
  

  //! computes fonction h_R(t) for different times
  void RamanEffect::InitCoefficient(const Real_wp& t)
  {
    if (alpha == Real_wp(0))
      return;

    int N = moduleE.GetM();
    hEval.Reallocate(N);
    for (int i = 0; i < N; i++)
      hEval(i) = EvaluateH(t-t_current+i*deltat);    
  }
  
  
  //! updates coefficients epsilon and gamma due to Raman contribution
  /*!
    gamma is replaced by gamma (1-alpha)
    and eps by eps + gamma alpha \int_{-\infty}^t R(t-s) |E|^2(s) ds
   */
  void RamanEffect::UpdateCoefficient(const Real_wp& t, int i, Real_wp& gamma,
				      Real_wp& eps, Real_wp& invEps)
  {
    if (alpha == Real_wp(0))
      return;

    int N = moduleE.GetM();
    
    // trapeze rule to obtain \int_{-\infty}^{t^n} R(t-s) |E(s)|^2 ds
    VectReal_wp EvalTrapeze(level_global+1);
    EvalTrapeze.Fill(0);
    
    Real_wp feval0 = hEval(0)*moduleE(num_current)(i);
    Real_wp fevalN = hEval(N-1)*moduleE((num_current+1)%N)(i);
    Real_wp length = length_interval;
    EvalTrapeze(0) = 0.5*(feval0+fevalN)*length;
    
    int pas = (N-1);
    for (int l = 1; l <= level_global; l++)
      {
	for (int k = pas/2; k < N-1; k += pas)
	  EvalTrapeze(l) += hEval(k)*moduleE((num_current-k+N)%N)(i);
	
	length /= 2; pas /= 2;
	EvalTrapeze(l) *= length;
	EvalTrapeze(l) += 0.5*EvalTrapeze(l-1);
      }
    
    // Romberg method to refine the result
    Real_wp powFour = 4.0;
    for (int m = 1; m <= min(level_romberg, level_global); m++)
      {
	Real_wp denom = 1.0/(powFour - 1);
	for (int n = level_global; n >= level_global-level_romberg+m; n--)
	  EvalTrapeze(n) = denom*(powFour*EvalTrapeze(n) - EvalTrapeze(n-1));
	
	powFour *= 4.0;
      }
    
    // final result for \int_{-\infty}^{t^n} R(t-s) |E(s)|^2 ds
    Real_wp res = EvalTrapeze(level_global);
    
    if (t < t_current-1e-8*deltat)
      {
	cout << "Time t must be greater than current time" << endl;
	DISP(t); DISP(t_current);
	abort();
      }
    
    // adding \int_{t^n}^t R(t-s) |E(s)|^2 ds 
    // trapezoid rule => 0.5 (t-t^n) R(t-t^n) |E(t^n)|^2
    res += 0.5*(t-t_current)*feval0;    
    
    eps += gamma*alpha*res;
    invEps = 1.0/eps;
    gamma *= 1.0-alpha;
  }
  
  
  //! default problem
  NonLinearOpticsProblem::NonLinearOpticsProblem()
  {
    gammaNL = 0.0;
    
    z0_adim = 1.0;
    E0_adim = 1.0;
    t0_adim = 1.0;
    
    tau_impulse = 1.0; Tmax = 5.0; amplitude_impulse = 1.0;
    nb_points_time = 1; padding_time = 0;

    theta = 0.25;
    theta_polar = 0.25;
    gamma_implicit = 0.0;

    vect_dt_gamma_implicit.Reallocate(1);
    vect_dt_gamma_implicit.Fill(0.0);

    vect_dt_gamma_implicit_cplx.Reallocate(1);
    vect_dt_gamma_implicit_cplx.Fill(Complex_wp(0,0));

    dt = 0.0;
    initial_time = 0.0; final_time = 0.0;
    
    zmin_display = 0.0;
    zmax_display = 0.0;    

    file_output_sismo = "sismo.dat";
    nb_points_sismo = 0;
    dt_sismo = 1.0;
    threshold_sismo = 100.0*epsilon_machine;

    nb_iterations_check_mesh = 0;
    threshold_translate = 1e-6;    
    velocity_translate = -1.0;
    t0_begin_translate = 0.0;
    dz_translation = 0.0; nb_iterations_translation = 0;
    
    omega_L = 0.0;
    
    coef_fft = 1.0;
    type_spectrum = CENTERED;
    number_odd_frequencies = 1;
    
    // circular case
    type_polarization = CIRCULAR;
    phase_init.Init(0.0, 0.5*pi_wp);
    polarization_init.Init(1.0, 1.0);

    initialisation_with_schrodinger = true;    
    harmonic_resolution = false;

  }
  

  Real_wp NonLinearOpticsProblem::GetAdimTime() const
  {
    return t0_adim;
  }

  
  void NonLinearOpticsProblem::SetFrequency(const Real_wp& f)
  {
    omega_L = 2.0*pi_wp*f;
  }
  
  
  //! modifies parameters of the problem with a line of the data file
  void NonLinearOpticsProblem
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "WavelengthAdim")
      {
        // WavelengthAdim = z0
        if (parameters.GetM() <= 0)
          {
	    cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "WaveLengthAdim needs one parameter, for instance :" << endl;
	    cout << "WaveLengthAdim = z0_adim" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        // the user can provide z0 : the characteristic length, which is usually the wavelength
        // then the characteristic time is equal to t0 = z0/c
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          {
            z0_adim = to_num<Real_wp>(parameters(0));
            t0_adim = z0_adim/PhysicalConstant::speed_light;
          }
      }
    else if (description_field == "NonLinearSusceptibility")
      {
        // NonLinearSusceptibility = gammaNL
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "NonLinearSusceptibility needs one parameter, for instance :" << endl;
	    cout << "NonLinearSusceptibility = sigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        gammaNL = to_num<Real_wp>(parameters(0));
        
        // if gammaNL is non-null, we set E0 = sqrt(epsilon0/gammaNL) 
        // in order to set the coefficient of non-linear term to 1
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          {
            if (gammaNL == 0)
              {
                E0_adim = 1.0;
              }
            else
              {
                E0_adim = sqrt(abs(PhysicalConstant::epsilon0_permittivity/gammaNL));
                if (gammaNL < 0)
                  gammaNL = -1.0;
                else
                  gammaNL = 1.0;
              }
          }
        else if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
          gammaNL /= PhysicalConstant::epsilon0_permittivity;
	
	if (parameters.GetM() > 1)
	  {
	    if (parameters(1) == "Raman")
	      {
		if (parameters.GetM() <= 4)
		  {
		    cout << "In SetInputData of NonLinearOpticsProblem" << endl;
		    cout << "NonLinearSusceptibility needs more parameters, for instance :" << endl;
		    cout << "NonLinearSusceptibility = gamma Raman tau1 tau2 alpha" << endl;
		    cout << "Current parameters are : " << endl << parameters << endl;
		    abort();
		  }
		
		Real_wp tau1 = to_num<Real_wp>(parameters(2));
		Real_wp tau2 = to_num<Real_wp>(parameters(3));
		Real_wp alpha = to_num<Real_wp>(parameters(4));
		tau1 /= t0_adim;
		tau2 /= t0_adim;
		var_raman.Init(tau1, tau2, 25.0*tau2, alpha);
	      }
	  }
      }
    else if (description_field == "DissipativeCoefficient")
      {
        // DissipativeCoefficient = sigma (conductivity)
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "DissipativeCoefficient needs one parameter, for instance :" << endl;
	    cout << "DissipativeCoefficient = sigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }
        
        Real_wp sigma = to_num<Real_wp>(parameters(0));
	index.SetSigma(sigma, z0_adim);
      }
    else if (description_field == "RefractiveIndex")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "RefractiveIndex needs at least one parameter, for instance :" << endl;
	    cout << "RefractiveIndex = Material" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        // in file PhysicalConstant.cxx
        index.SetRefractiveIndex(parameters(0), parameters);
	
        // adimensionalization step
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          {
            Mlt(t0_adim, index.omega_polarization);
            Mlt(1.0/t0_adim, index.sigma_polarization);
          }
      }
    else if (description_field == "TimeStep")
      {
        // TimeStep = dt
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "TimeStep needs one parameter, for instance :" << endl;
	    cout << "TimeStep = dt" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        dt = to_num<Real_wp>(parameters(0));
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          dt /= t0_adim;
      }
    else if (description_field == "TimeInterval")
      {
        // TimeInterval = initial_time final_time
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "TimeInterval needs two parameters, for instance :" << endl;
	    cout << "TimeInterval = t0 tmax" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        initial_time = to_num<Real_wp>(parameters(0));
        final_time = to_num<Real_wp>(parameters(1));
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          {
            initial_time /= t0_adim;
            final_time /= t0_adim;
          }
      }
    else if (description_field == "TimeScheme")
      {
        // TimeScheme = theta theta_polar
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "TimeScheme needs two parameters, for instance :" << endl;
	    cout << "TimeScheme = theta theta_polar" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        theta = to_num<Real_wp>(parameters(0));
        theta_polar = to_num<Real_wp>(parameters(1));
      }
    else if (description_field == "TimeImpulsion")
      {
        // TimeImpulsion = N N_padding Gaussian tau Tmax amplitude
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "TimeImpulsion needs at least three parameters, for instance :" << endl;
	    cout << "TimeImpulsion = N Npadding type_impulsion" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
	
        nb_points_time = to_num<int>(parameters(0));
        padding_time = to_num<int>(parameters(1));
        if ((parameters(2) == "Gaussian") || (parameters(2) == "GaussianHalf") || (parameters(2) == "GaussianMiro"))
          {
            if (parameters.GetM() <= 5)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "TimeImpulsion needs at least six parameters, for instance :" << endl;
                cout << "TimeImpulsion = N Npadding Gaussian tau Tmax amplitude" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
	    if (parameters(2) == "GaussianHalf")
              {
                Real_wp fwhm_impulse = to_num<Real_wp>(parameters(3));
                //Relation between fwhm and tau
                tau_impulse =fwhm_impulse /(Real_wp(2) * sqrt(Real_wp(2)*log(Real_wp(2)))); 
              }
	    else if (parameters(2) == "GaussianMiro")
              {
                Real_wp fwhm_impulse = to_num<Real_wp>(parameters(3));
                //Relation between fwhm and tau
                tau_impulse = fwhm_impulse /(Real_wp(2) * sqrt(log(Real_wp(2)))); 
              }
            else
              tau_impulse = to_num<Real_wp>(parameters(3));
            
            Tmax = to_num<Real_wp>(parameters(4));
            if (parameters.GetM() > 5)
              amplitude_impulse = to_num<Real_wp>(parameters(5));
            
            if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
              {
                tau_impulse /= t0_adim;
                Tmax /= t0_adim;
                amplitude_impulse /= E0_adim;
              }

	    //DISP(tau_impulse); DISP(Tmax);
          }
        else if (parameters(2) == "Monochromatic")
          {
            nb_points_time = 1;
          }
        else
          {
	    cout << "Impulsion not known" << endl;
            abort();
          }
      }
    else if (description_field == "FourierTransform")
      {
        // FourierTransform = type_spectrum Tmax nb_points
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
            cout << "FourierTransform needs at least three parameters, for instance :" << endl;
            cout << "FourierTransform = type_spectrum Tmax nb_points" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        Real_wp fmin, fmax, dt;
        int nb_omega(2);
        if (parameters(0) == "TIME_CENTER")
          {
	    type_spectrum = CENTERED;
            Real_wp Tm = to_num<Real_wp>(parameters(1));
            nb_omega = to_num<int>(parameters(2));
            dt = 2.0*Tm / (t0_adim*nb_omega);
            fmin = omega_L/(2.0*pi_wp) - 0.5/dt;
            fmax = fmin + 1.0/dt;
          }
        else if (parameters(0) == "ODD_FREQUENCIES")
          {
            if (parameters.GetM() <= 3)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "FourierTransform needs four parameters, for instance :" << endl;
                cout << "FourierTransform = ODD_FREQUENCIES Tmax nb_points p" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }

            type_spectrum = ODD_FREQUENCIES;
            Real_wp Tm = to_num<Real_wp>(parameters(1));
            nb_omega = to_num<int>(parameters(2));
            number_odd_frequencies = to_num<int>(parameters(3));
            dt = 2.0*Tm / (t0_adim*nb_omega);
            fmin = omega_L/(2.0*pi_wp) - 0.5/dt;
            int N = 2*(number_odd_frequencies-1)+1;
            fmax = N*omega_L/(2.0*pi_wp) + 0.5/dt;
          }
        else
          {
	    type_spectrum = SHIFTED_SINGLE;
            Real_wp Tm = to_num<Real_wp>(parameters(1));
            nb_omega = to_num<int>(parameters(2));
            dt = 2.0*Tm / (t0_adim*nb_omega);
            fmin = omega_L/(2.0*pi_wp) - 0.5/dt;
            fmax = fmin + 1.0/dt;
            if (fmin < 0.0)
              fmin = Real_wp(0);
          }
        
        Linspace(fmin, fmax, nb_omega, omega_fourier);
        Mlt(2.0*pi_wp, omega_fourier);
      }
    else if (description_field == "Wavelength")
      {
        // wavelength in the vacuum
	// Wavelength = lambda
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "Wavelength needs one parameter, for instance :" << endl;
	    cout << "Wavelength = lambda" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

	Real_wp lambda = to_num<Real_wp>(parameters(0));
	Real_wp f;
	if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_ONE)
	  f = index.c0/lambda;
	else
	  f = PhysicalConstant::speed_light/lambda*t0_adim;
	
	SetFrequency(f);
      }
    else if (description_field == "Frequency")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "Frequency needs one parameter, for instance :" << endl;
	    cout << "Frequency = f" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        Real_wp f = to_num<Real_wp>(parameters(0));
        SetFrequency(f*t0_adim);
      }
    else if (description_field == "SismoPoint")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "SismoPoint needs at least one parameter, for instance :" << endl;
	    cout << "SismoPoint = x" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        if (parameters(0) == "LINE")
          {
            // SismoPoint = LINE z0 z1 dz
            if (parameters.GetM() <= 3)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "SismoPoint needs four parameters, for instance :" << endl;
                cout << "SismoPoint = LINE z0 z1 dz" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }

            Real_wp z0 = to_num<Real_wp>(parameters(1))/z0_adim;
            Real_wp z1 = to_num<Real_wp>(parameters(2))/z0_adim;
            Real_wp dz = to_num<Real_wp>(parameters(3))/z0_adim;
            
            int N = toInteger(round((z1 - z0)/dz)) + 1;
            VectReal_wp points;
            if (N <= 1)
              {
                points.Reallocate(2);
                points(0) = z0;
                points(1) = z1;
              }
            else
              {
                points.Reallocate(N);
                for (int i = 0; i < N; i++)
                  points(i) = z0 + i*dz;
              }
            
            for (int i = 0; i < points.GetM(); i++)
              grid_sismo.AddPoint(points(i));
          }
        else
          {
            // SismoPoint = z
            grid_sismo.AddPoint(to_num<Real_wp>(parameters(0))/z0_adim);
          }
      }
    else if (description_field == "FileOutputSismoPoint")
      {
        // FileOutputSismoPoint = file_name dt
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "FileOutputSismoPoint needs at least two parameters, for instance :" << endl;
	    cout << "FileOutputSismoPoint = file_name dt" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        file_output_sismo = parameters(0);
        dt_sismo = to_num<Real_wp>(parameters(1))/t0_adim;
        if (parameters.GetM() > 2)
          threshold_sismo = to_num<Real_wp>(parameters(2));
	
	if (parameters.GetM() > 3)
	  nb_points_sismo = to_num<int>(parameters(3));
      }
    else if (description_field == "TranslateMesh")
      {
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
	    cout << "Translate needs at least two parameters, for instance :" << endl;
	    cout << "TranslateMesh = 100 epsilon" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        nb_iterations_check_mesh = to_num<int>(parameters(0));
        velocity_translate = -1.0;
        if (parameters(1) == "VELOCITY")
          {
            if (parameters.GetM() <= 4)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "TranslateMesh needs five parameters, for instance :" << endl;
                cout << "TranslateMesh = 100 VELOCITY v0 t0 epsilon" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }

            velocity_translate = to_num<Real_wp>(parameters(2));
            t0_begin_translate = to_num<Real_wp>(parameters(3));
            threshold_translate = to_num<Real_wp>(parameters(4));
            velocity_translate *= t0_adim/z0_adim;
            t0_begin_translate /= t0_adim;
          }
        else if (parameters(1) == "AUTO")
          {
           if (parameters.GetM() <= 3)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "TranslateMesh needs four parameters, for instance :" << endl;
                cout << "TranslateMesh = 100 AUTO t0 epsilon" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }

            velocity_translate = 0.0;
            t0_begin_translate = to_num<Real_wp>(parameters(2));
            threshold_translate = to_num<Real_wp>(parameters(3));
            t0_begin_translate /= t0_adim;
          }
        else
          threshold_translate = to_num<Real_wp>(parameters(1));
      }
    else if (description_field == "Polarization")
      {	
        // Polarization = LINEAR
        // Polarization = CIRCULAR
        // Polarization = ELLIPTIC Ex Ey phix phiy
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
            cout << "Polarization needs at least one parameter, for instance :" << endl;
            cout << "Polarization = type_polarization" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        // type of polarization
    	if (parameters(0) == "Linear")
    	  {
            type_polarization = LINEAR;
            phase_init.Fill(0.0);
            polarization_init.Init(1.0, 0.0);
          }
    	else if (parameters(0)== "Circular")
    	  {
            type_polarization = CIRCULAR;
            phase_init.Init(0.0, 0.5*pi_wp);
            polarization_init.Init(1.0, 1.0);
          }
    	else
    	  {
            if (parameters.GetM() <= 4)
              {
                cout << "In SetInputData of NonLinearOpticsProblem" << endl;
                cout << "Polarization needs four parameters, for instance :" << endl;
                cout << "Polarization = Elliptic Ex Ey phi_x phi_y" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
            type_polarization = ELLIPTIC;
            polarization_init(0) = to_num<Real_wp>(parameters(1));
            polarization_init(1) = to_num<Real_wp>(parameters(2));
            phase_init(0) = pi_wp/180.0*to_num<Real_wp>(parameters(3));
            phase_init(1) = pi_wp/180.0*to_num<Real_wp>(parameters(4));
          }
      }
    else if (description_field == "InitGuessNewton")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
            cout << "InitGuessNewton needs one parameter, for instance :" << endl;
            cout << "InitGuessNewton = Schrodinger" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
            
	if (parameters(0) == "Schrodinger")
	  initialisation_with_schrodinger = true;
	else
	  initialisation_with_schrodinger = false;
      }
    else if (description_field == "PrintLevel")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of NonLinearOpticsProblem" << endl;
            cout << "PrintLevel needs one parameter, for instance :" << endl;
            cout << "PrintLevel = level" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        print_level = to_num<int>(parameters(0));
      }
    else if (description_field == "OrderTimeScheme")
      param_time_scheme = parameters;
  }
  
  
  //! sets resolution in time-harmonic
  void NonLinearOpticsProblem::SetHarmonicResolution()
  {
    harmonic_resolution = true;
  }
   
  
  //! constructs spectrum to be computed
  void NonLinearOpticsProblem::ConstructSpectrum(const string& DOSSIER)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    // impulsion is defined on interval [-Tmax, Tmax]
    omega.Reallocate(nb_points_time);
    kwave.Reallocate(nb_points_time);
    time_interval.Reallocate(nb_points_time);
    Real_wp dt = 2.0*Tmax/Real_wp(nb_points_time);
    time_interval.Fill(0);
    for (int i = 0; i < nb_points_time/2; i++)
      {
        time_interval(i) = Real_wp(i)*dt;
        time_interval(i+nb_points_time/2) = Real_wp(i-nb_points_time/2)*dt;
      }

    // calcul de l'indice pour une gamme de longueurs d'ondes
    if (rank_proc == 0)
      if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
        {
          int Nl = 1000;
          VectReal_wp Lambda;
          VectReal_wp indice_optique(Nl);
          Linspace(Real_wp(0.4e-6), Real_wp(3.5e-6), Nl, Lambda);
          for (int i = 0; i < Lambda.GetM(); i++)
            {
              Real_wp pulsation = (2.0*pi_wp*index.c0/Lambda(i))*t0_adim;
              indice_optique(i) = index.c0*real(index.GetWaveNumber(pulsation))/pulsation;
            }
          
          Lambda.WriteText(DOSSIER + "wavelength.dat");
          indice_optique.WriteText(DOSSIER + "indice.dat");
        }
    
    // frequencies are in the interval [-1/2 dt, 1/2 dt]
    VectReal_wp Omega;
    Linspace(-pi_wp/dt, pi_wp/dt, nb_points_time+1, Omega);    
    
    if (type_spectrum == SHIFTED_SINGLE)
      {
        omega.Reallocate(nb_points_time);
        omega.Fill(omega_L); kwave.Fill(0.0);
	for (int i = 0; i < nb_points_time/2; i++)
          {
            omega(i) = omega_L + Omega(nb_points_time/2+i);
            kwave(i) = index.GetWaveNumber(omega(i));
            
            Real_wp omega_ = omega_L + Omega(i);
            if ((omega_L > 0.0) && (omega_ < 0))
              {
                DISP(omega_L); DISP(dt); DISP(Omega.GetNormInf()); 
                cout << "The spectrum is assumed to be around the value " << omega_L << endl;
                cout << "Therefore the pulsation omega should be always positive" << endl;
                cout << "But it is equal to " << omega_ << endl;
                abort();
              }
            
            omega(i+nb_points_time/2) = omega_;
            kwave(i+nb_points_time/2) = index.GetWaveNumber(omega_);
	  }
        
        // writing pulsations and wave numbers on the disk
        if (rank_proc == 0)
          {
            VectReal_wp omega_display(nb_points_time);
            VectComplex_wp kwave_display(nb_points_time);
            for (int i = 0; i < nb_points_time/2; i++)
              {
                omega_display(nb_points_time/2+i) = omega(i);
                omega_display(i) = omega(nb_points_time/2+i);
                kwave_display(nb_points_time/2+i) = kwave(i);
                kwave_display(i) = kwave(nb_points_time/2+i);
              }
            
            Mlt(1.0/t0_adim, omega_display);
            Mlt(1.0/z0_adim, kwave_display);
            
            omega_display.Write(DOSSIER + "omega.dat");
            kwave_display.Write(DOSSIER + "kwave.dat");    
          }

        // initializing FFT
        fft_envelope.Init(nb_points_time + padding_time);
        coef_fft = 2.0*Tmax/sqrt(2.0*pi_wp);
      }
    else if (type_spectrum == ODD_FREQUENCIES)
      {
        int N = number_odd_frequencies;
        omega.Reallocate(nb_points_time*N);
        kwave.Reallocate(nb_points_time*N);
        omega.Fill(omega_L); kwave.Fill(0.0);
        // pulsations to be written
        VectReal_wp omega_display(nb_points_time*N);
        VectComplex_wp kwave_display(nb_points_time*N);
        int offset = 0;
        for (int k = 0; k < N; k++)
          {
            for (int i = 0; i < nb_points_time/2; i++)
              {
                omega(offset+i) = (2*k+1)*omega_L + Omega(nb_points_time/2+i);
                kwave(offset+i) = index.GetWaveNumber(omega(offset+i));
                
                Real_wp omega_ = (2*k+1)*omega_L + Omega(i);
                if ((omega_L > 0.0) && (omega_ < 0))
                  {
                    DISP(omega_L); DISP(dt); DISP(Omega.GetNormInf()); 
                    cout << "The spectrum is assumed to be around the value " << omega_L << endl;
                    cout << "Therefore the pulsation omega should be always positive" << endl;
                    cout << "But it is equal to " << omega_ << endl;
                    abort();
                  }
                
                omega(offset+i+nb_points_time/2) = omega_;
                kwave(offset+i+nb_points_time/2) = index.GetWaveNumber(omega_);
                
                omega_display(offset+nb_points_time/2+i) = omega(offset+i);
                omega_display(offset+i) = omega(offset+nb_points_time/2+i);
                kwave_display(offset+nb_points_time/2+i) = kwave(offset+i);
                kwave_display(offset+i) = kwave(offset+nb_points_time/2+i);
              }
            offset += nb_points_time;
          }
        
        Mlt(1.0/t0_adim, omega_display);
        Mlt(1.0/z0_adim, kwave_display);
        
        if (rank_proc == 0)
          {
            omega_display.Write(DOSSIER + "omega.dat");
            kwave_display.Write(DOSSIER + "kwave.dat");    
          }
        
        // initializing FFT
        fft_envelope.Init(nb_points_time + padding_time);
        coef_fft = 2.0*Tmax/sqrt(2.0*pi_wp);
      }
    else
      {
        omega.Reallocate(nb_points_time/2+1);
	kwave.Reallocate(nb_points_time/2+1);
        for (int i = 1; i <= nb_points_time/2; i++)
          {
            omega(i) = Omega(nb_points_time/2+i);
            kwave(i) = index.GetWaveNumber(omega(i));
	  }
	
	omega(0) = 0; kwave(0) = 0;
	
        // writing pulsations and wave numbers on the disk
        if (rank_proc == 0)
          {
            VectReal_wp omega_display(omega);
            VectComplex_wp kwave_display(kwave);
            
            Mlt(1.0/t0_adim, omega_display);
            Mlt(1.0/z0_adim, kwave_display);
            
            omega_display.Write(DOSSIER + "omega.dat");
            kwave_display.Write(DOSSIER + "kwave.dat");    
          }

        // initializing FFT
        fft_time.Init(nb_points_time + padding_time);
        coef_fft = 4.0*Tmax/(nb_points_time+padding_time)/sqrt(2.0*pi_wp); 
	
	// Fourier transform of Raman function
	if (var_raman.GetAlpha() != Real_wp(0))
	  {
	    VectReal_wp raman_time;
	    raman_time.Reallocate(nb_points_time);
	    raman_fct_hat.Reallocate(nb_points_time/2+1);
	    
	    for (int i = 0; i < nb_points_time; i++)
	      raman_time(i) = var_raman.EvaluateH(time_interval(i));
	    
	    raman_time.Write("raman.dat");
	    fft_time.ApplyForward(raman_time, raman_fct_hat);
	    Mlt(coef_fft, raman_fct_hat);
	    
	    for (int i = 0; i < raman_fct_hat.GetM(); i++)
	      raman_fct_hat(i) = conj(raman_fct_hat(i));
	  }
      }
  }
  
  
  //! returns non-homogeneous Dirichlet condition f(t^n)
  Complex_wp NonLinearOpticsProblem::GetPulseTime(const Real_wp& t)
  {
    Complex_wp pulse;
    GetPulseTime(t, pulse);
    return pulse;
  }
  
  
  //! returns non-homogeneous Dirichlet condition f(t^n)
  void NonLinearOpticsProblem::GetPulseTime(const Real_wp& t, Complex_wp& pulse_c)
  {
    pulse_c =  Complex_wp(0.0, 0.0);
    if (t > 2.0*Tmax)
      return;

    Real_wp pulse = exp(-0.5*square((t-Tmax)/tau_impulse));
    Real_wp pulse_x = pulse*polarization_init(0)*cos(omega_L*t + phase_init(0));
    Real_wp pulse_y = pulse*polarization_init(1)*cos(omega_L*t + phase_init(1));
    pulse_x *= amplitude_impulse;    pulse_y *= amplitude_impulse;    
    pulse_c = Complex_wp(pulse_x, pulse_y);
    //DISP(t); DISP(Tmax); DISP(tau_impulse);  DISP(pulse); DISP(omega_L); DISP(phase_init(0));
    //DISP(polarization_init); DISP(pulse_x); DISP(pulse_y);
  }


  //! returns non-homogeneous Dirichlet condition f(t^n) (linear polarization)
  void NonLinearOpticsProblem::GetPulseTime(const Real_wp& t, Real_wp& pulse_c)
  {
    pulse_c = 0;
    if (t > 2.0*Tmax)
      return;
    
    Real_wp pulse = exp(-0.5*square((t-Tmax)/tau_impulse));
    pulse_c = pulse*polarization_init(0)*cos(omega_L*t + phase_init(0));
    pulse_c *= amplitude_impulse; DISP(amplitude_impulse);
  }
  
  
  //! returns Fourier transform of non-homogeneous Dirichlet condition
  Complex_wp NonLinearOpticsProblem::GetPulseHarmonic(const Real_wp& omega)
  {
    Complex_wp coef = amplitude_impulse*tau_impulse;
    Complex_wp pulse = coef*exp(-0.5*square(tau_impulse*(omega-omega_L)));
    return pulse;
  }
  

  //! computes the electric field E from displacement D
  void NonLinearOpticsProblem
  ::GetElectricFieldFromDisplacement(const Real_wp& t, const VectComplex_wp& D,
                                     VectComplex_wp& E)
  {
    // performance tests completed in file newton_vec.cc
    // we keep the most efficient choice
    Real_wp err, Ex, Ey;
    TinyMatrix<Real_wp, Symmetric, 2, 2> mat_DF;
    TinyVector<Real_wp, 2> vec_F, vec_V;
    Real_wp beta_x, beta_y, E_module, inv_det;    
    Real_wp invEps = 1.0/index.eps_inf, eps = index.eps_inf;
    Real_wp gamma = gammaNL;
    int N = E.GetM();
    if (gammaNL == 0)
      {
        #pragma omp parallel for
	for (int i = 0; i < N; i++)
	  {
	    beta_x = real(D(i)); beta_y = imag(D(i));
	    
	    // subtracting polarization part
	    // beta = D - \sum_k P_k
	    for (int k = 0; k < index.omega_polarization.GetM(); k++)
	      {
		beta_x -= real(D(N*(k+1) + i));
		beta_y -= imag(D(N*(k+1) + i));
	      }
        
	    Ex = beta_x*invEps;
	    Ey = beta_y*invEps;
	    E(i) = Complex_wp(Ex, Ey);	
	  }
	
	return ;
      }

    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::FLUX);
    
    int offset_raman = N*(index.omega_polarization.GetM()+1);
    
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
      {
        // we want to solve the following equation on point i :
        // eps_inf E + \sum_k P_k + \gamma_{NL} |E|^2 E = D 
        beta_x = real(D(i)); beta_y = imag(D(i));
        
        // subtracting polarization part
        // beta = D - \sum_k P_k
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          {
            beta_x -= real(D(N*(k+1) + i));
            beta_y -= imag(D(N*(k+1) + i));
          }
        
	// updating eps and gamma due to Raman effect
	gamma = gammaNL; eps = index.eps_inf;
	if (var_raman.GetAlpha() != Real_wp(0))
	  {
	    eps += gammaNL*var_raman.GetAlpha()*real(D(offset_raman+i));
	    gamma *= 1.0-var_raman.GetAlpha();
	    invEps = 1.0/eps;
	  }
		
        // applying newton's Method to solve
        // eps_inf E + \gamma_{NL} |E|^2 E = beta
        Ex = beta_x*invEps;
        Ey = beta_y*invEps;
	
        E_module = eps + gamma*(Ex*Ex + Ey*Ey);
        vec_F(0) = Ex*E_module - beta_x;
        vec_F(1) = Ey*E_module - beta_y;
	
        mat_DF(0, 0) = eps + gamma*(Ey*Ey + 3.0*Ex*Ex);
        mat_DF(0, 1) = 2.0*gamma*Ex*Ey;
        mat_DF(1, 1) = eps + gamma*(Ex*Ex + 3.0*Ey*Ey);
        inv_det = 1.0/Det(mat_DF);
	
        vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
        vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
	
        Ex -= vec_V(0);
        Ey -= vec_V(1);
	
        E_module = eps + gamma*(Ex*Ex + Ey*Ey);
        vec_F(0) = Ex*E_module - beta_x;
        vec_F(1) = Ey*E_module - beta_y;
	
        mat_DF(0, 0) = eps + gamma*(Ey*Ey + 3.0*Ex*Ex);
        mat_DF(0, 1) = 2.0*gamma*Ex*Ey;
        mat_DF(1, 1) = eps + gamma*(Ex*Ex + 3.0*Ey*Ey);
        inv_det = 1.0/Det(mat_DF);
	
        vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
        vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
	
        Ex -= vec_V(0);
        Ey -= vec_V(1);
	
        // this last test is usually not needed, since Newton's method
        // should have converged, but for intense fields, as a protection...
        E_module = eps + gamma*(Ex*Ex + Ey*Ey);
        vec_F(0) = Ex*E_module - beta_x;
        vec_F(1) = Ey*E_module - beta_y;
        err = Norm1(vec_F);
        if (err > 1e-14)
          {
            mat_DF(0, 0) = eps + gamma*(Ey*Ey + 3.0*Ex*Ex);
            mat_DF(0, 1) = 2.0*gamma*Ex*Ey;
            mat_DF(1, 1) = eps + gamma*(Ex*Ex + 3.0*Ey*Ey);
            inv_det = 1.0/Det(mat_DF);
            
            vec_V(0) = inv_det*(mat_DF(1, 1)*vec_F(0) - mat_DF(0, 1)*vec_F(1));
            vec_V(1) = inv_det*(mat_DF(0, 0)*vec_F(1) - mat_DF(0, 1)*vec_F(0));
            
            Ex -= vec_V(0);
            Ey -= vec_V(1);	      
          }
        
        E(i) = Complex_wp(Ex, Ey);
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::FLUX);
  }
    
  
  //! gets dE/dt from dD/dt and E
  void NonLinearOpticsProblem
  ::GetDeriveElectricFieldFromDisplacement(const Real_wp& t, const VectComplex_wp& D,
                                           const VectComplex_wp& dD_dt,
					   const VectComplex_wp& E, VectComplex_wp& dE_dt)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::FLUX);
    
    TinyVector<Real_wp, 2> beta, vec_V;
    TinyMatrix<Real_wp, Symmetric, 2, 2> mat_DF;
    Real_wp Ex, Ey, inv_det;
    Real_wp eps = index.eps_inf;
    Real_wp gamma = gammaNL;
    int N = E.GetM();
    int offset_raman = N*(index.omega_polarization.GetM()+1);
    
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
      {
	// the derivative in time of the equation
        // eps_inf E + \sum_k P_k + \gamma_{NL} |E|^2 E = D 
	// is equal to :
	// eps_inf dE/dt + \sum _k dP_k/dt + \gamma_{NL} Jacob(|E|^2 E) dE/dt
        beta(0) = real(dD_dt(i)); beta(1) = imag(dD_dt(i));
        
        // subtracting polarization part
        // beta = dD/dt - \sum_k dP_k/dt
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          {
            beta(0) -= real(dD_dt(N*(k+1) + i));
            beta(1) -= imag(dD_dt(N*(k+1) + i));
          }
	
	// modification of gamma and epsilon due to Raman effect
	if (var_raman.GetAlpha() != Real_wp(0))
	  {
	    eps += gammaNL*var_raman.GetAlpha()*real(D(offset_raman+i));
	    gamma *= 1.0-var_raman.GetAlpha();
	    
	    Real_wp coef_raman = real(dD_dt(offset_raman+i))*gammaNL*var_raman.GetAlpha();
	    beta(0) -= coef_raman*real(E(i));
	    beta(1) -= coef_raman*imag(E(i));
	  }
	
	Ex = real(E(i)); Ey = imag(E(i));
	
	// evaluating Jacobian
	mat_DF(0, 0) = eps + gamma*(Ey*Ey + 3.0*Ex*Ex);
        mat_DF(0, 1) = 2.0*gamma*Ex*Ey;
        mat_DF(1, 1) = eps + gamma*(Ex*Ex + 3.0*Ey*Ey);
	
	// then solving the system to obtain dE/dt
	inv_det = 1.0/Det(mat_DF);
	
        vec_V(0) = inv_det*(mat_DF(1, 1)*beta(0) - mat_DF(0, 1)*beta(1));
        vec_V(1) = inv_det*(mat_DF(0, 0)*beta(1) - mat_DF(0, 1)*beta(0));
	
	// storing values of dE/dt
	dE_dt(i) = Complex_wp(vec_V(0), vec_V(1));
      }

    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::FLUX);
  }


  //! computes the electric field E from displacement D (linear polarization)
  void NonLinearOpticsProblem
  ::GetElectricFieldFromDisplacement(const Real_wp& t, const VectReal_wp& D, VectReal_wp& E)
  {
    // performance tests completed in file newton.cc
    // we keep the most efficient choice
    Real_wp err, Ex, DF, F;
    Real_wp beta, E_module;    
    Real_wp invEps = 1.0/index.eps_inf, eps = index.eps_inf;
    Real_wp gamma = gammaNL;
    int N = E.GetM();
    if (gammaNL == 0)
      {
        #pragma omp parallel for
	for (int i = 0; i < N; i++)
	  {
	    beta = D(i);
	    
	    // subtracting polarization part
	    // beta = D - \sum_k P_k
	    for (int k = 0; k < index.omega_polarization.GetM(); k++)
              beta -= D(N*(k+1) + i);
        
	    E(i) = beta*invEps;
	  }
	
	return ;
      }

    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::FLUX);

    int offset_raman = N*(index.omega_polarization.GetM()+1);
    
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
      {
        // we want to solve the following equation on point i :
        // eps_inf E + \sum_k P_k + \gamma_{NL} |E|^2 E = D 
        beta = D(i);
        
        // subtracting polarization part
        // beta = D - \sum_k P_k
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          beta -= D(N*(k+1) + i);
        
	// updating eps and gamma due to Raman effect
	gamma = gammaNL; eps = index.eps_inf;
	if (var_raman.GetAlpha() != Real_wp(0))
	  {
	    eps += gammaNL*var_raman.GetAlpha()*D(offset_raman+i);
	    gamma *= 1.0-var_raman.GetAlpha();
	  }
	
	invEps = 1.0/eps;
	
        // applying newton's Method to solve
        // eps_inf E + \gamma_{NL} |E|^2 E = beta
        Ex = beta*invEps;
        
        E_module = index.eps_inf + gammaNL*Ex*Ex;
        F = Ex*E_module - beta;        	
        DF = index.eps_inf + 3.0*gammaNL*Ex*Ex;
        Ex -= F/DF;

        E_module = index.eps_inf + gammaNL*Ex*Ex;
        F = Ex*E_module - beta;        	
        DF = index.eps_inf + 3.0*gammaNL*Ex*Ex;
        Ex -= F/DF;

        // this last test is usually not needed, since Newton's method
        // should have converged, but for intense fields, as a protection...
        E_module = index.eps_inf + gammaNL*Ex*Ex;
        F = Ex*E_module - beta;        	
        err = abs(F);
        if (err > 1e-14)
          {
            DF = index.eps_inf + 3.0*gammaNL*Ex*Ex;
            Ex -= F/DF;
          }
        
        E(i) = Ex;
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::FLUX);
  }
    
  
  //! gets dE/dt from dD/dt and E
  void NonLinearOpticsProblem
  ::GetDeriveElectricFieldFromDisplacement(const Real_wp& t, const VectReal_wp& D,
                                           const VectReal_wp& dD_dt,
					   const VectReal_wp& E, VectReal_wp& dE_dt)
  {
    if (print_level >= 2)
      glob_chrono.Start(VirtualTimer::FLUX);
    
    Real_wp eps = index.eps_inf;
    Real_wp gamma = gammaNL;
    
    Real_wp beta, Ex, DF;
    int N = E.GetM();
    int offset_raman = N*(index.omega_polarization.GetM()+1);
 
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
      {
	// the derivative in time of the equation
        // eps_inf E + \sum_k P_k + \gamma_{NL} |E|^2 E = D 
	// is equal to :
	// eps_inf dE/dt + \sum _k dP_k/dt + \gamma_{NL} Jacob(|E|^2 E) dE/dt
        beta = dD_dt(i);
        
        // subtracting polarization part
        // beta = dD/dt - \sum_k dP_k/dt
        for (int k = 0; k < index.omega_polarization.GetM(); k++)
          beta -= dD_dt(N*(k+1) + i);

	// modification of gamma and epsilon due to Raman effect
	if (var_raman.GetAlpha() != Real_wp(0))
	  {
	    eps += gammaNL*var_raman.GetAlpha()*D(offset_raman+i);
	    gamma *= 1.0-var_raman.GetAlpha();
	    
	    Real_wp coef_raman = dD_dt(offset_raman+i)*gammaNL*var_raman.GetAlpha();
	    beta -= coef_raman*E(i);
	  }
	
	Ex = E(i);
	
	// evaluating Jacobian
        DF = eps + 3.0*gamma*Ex*Ex;
        
	// then solving the system to obtain dE/dt
	dE_dt(i) = beta / DF;
      }
    
    if (print_level >= 2)
      glob_chrono.Stop(VirtualTimer::FLUX);
  }


  /****************
   * OpticalFiber *
   ****************/


  OpticalFiber::OpticalFiber()
  {
    a = 1.0; n1 = 1.5; n2 = 1.4;
    type_solver = DICHOTOMY_SOLVER;
    //type_solver = NEWTON_SOLVER;
  }
  
  
  void OpticalFiber::SetCircularFiber(const Real_wp& b, const Real_wp& n1_, const Real_wp& n2_)
  {
    a = b;
    n1 = n1_;
    n2 = n2_;
  }
  
  
  Real_wp OpticalFiber::EvaluateFunction(int m, const Real_wp& k0, const Real_wp& beta,
					 const Real_wp& Q, const Real_wp& constte)
  {
    Real_wp u, w, s, b1, b1p, b2, b2p, f1, f2, f;
    u = sqrt(n1*n1 - beta*beta);
    w = sqrt(beta*beta - n2*n2);
    s = k0*a;
    
    b1 = jn(m, s*u);
    b1p = (jn(m-1, s*u) - jn(m+1, s*u))/2;
    
    b2 = GetKn(m, s*w);
    b2p = -(GetKn(m-1, s*w) + GetKn(m+1, s*w))/2;
    
    f1 = b1p/((s*u)*b1);
    f2 = b2p/((s*w)*b2);
    
    f = (f1+f2)*(f1+Q*f2) - constte*(beta*beta)/pow(u*w, 4);
    return f;
  }
  
  
  Real_wp OpticalFiber::GetKomega(int m, const Real_wp& omega, const Real_wp& c0)
  {
    Real_wp Q = square(n2/n1);
    Real_wp k0 = omega/c0;
    
    Real_wp V = k0*a*sqrt(n1*n1-n2*n2); 
    if ((V > 2.405) || (std::isnan(V)))
      {
	cout << "Not a monomode fiber" << endl;
	abort();
      }
    
    Real_wp constte = square(m*n1*(1.0-Q)/square(k0*a));
    Real_wp f1, f2, f1prime, f2prime, f, fprime;
    
    Real_wp v1 = 1e-5*(n1-n2);
    if (type_solver == DICHOTOMY_SOLVER)
      {
	// first, we find a point where f > 0
	Real_wp val[7] = {1e-12, 1e-9, 1e-6, 1e-3, 1e-2, 0.1, 0.25};
	v1 = 0;
	for (int k = 0; k < 7; k++)
	  {
	    f1 = EvaluateFunction(m, k0, n2+val[k]*(n1-n2), Q, constte);
	    if (f1 > 0)
	      v1 = val[k]*(n1-n2);
	  }
	
	// second we find a point where f < 0
	Real_wp v2 = 0;
	for (int k = 0; k < 7; k++)
	  {
	    f2 = EvaluateFunction(m, k0, n1+val[k]*(n2-n1), Q, constte);
	    if (f2 < 0)
	      v2 = val[k]*(n2-n1);
	  }
	
	if (v1 == 0)
	  {
	    // in that case, we assume that the solution
	    // is too close to n2 to see a positive value of f
	    return n2;
	  }
	
	if (v2 == 0)
	  {
	    // the solution is too close to n1 to see a negative value of f
	    return n1;
	  }
	
	// Dichotomy algorithm to find the solution
	Real_wp a_sol, b_sol, c_sol(1);
	Real_wp threshold = 10.0*epsilon_machine;
	a_sol = n2 + v1;
	b_sol = n1 + v2;
	while ((b_sol - a_sol) > threshold)
	  {
	    c_sol = (a_sol+b_sol)/2;
	    f = EvaluateFunction(m, k0, c_sol, Q, constte);
	    if (f < 0)
	      b_sol = c_sol;
	    else
	      a_sol = c_sol;
	  }
	
	return c_sol;
      }
    
    // Newton's algorithm
    Real_wp beta = n2+v1;
    Real_wp f_prec = 1e100;
    int k = 0;
    
    Real_wp u, w, s, b1, b1p, b1s, b2, b2p, b2s;
    f = f_prec / 2;
    while ((k <= 100) && (abs(f) < abs(f_prec)))
      {
	u = sqrt(n1*n1 - beta*beta);
	w = sqrt(beta*beta - n2*n2);
	s = k0*a;
        
	b1 = jn(m, s*u);
	b1p = (jn(m-1, s*u) - jn(m+1, s*u))/2;
	b1s = -b1p/(s*u) -(1.0-square(m/(s*u)))*b1;
        
	b2 = GetKn(m, s*w);
	b2p = -(GetKn(m-1, s*w) + GetKn(m+1, s*w))/2;
	b2s = -b2p/(s*w) + (1.0+square(m/(s*w)))*b2;
	
	f1 = b1p/((s*u)*b1);
	f2 = b2p/((s*w)*b2);
        
	f1prime = -beta*b1s/(b1*u*u)+ f1*beta/(u*u) + beta*s*s*f1*f1;
	f2prime = beta*b2s/(b2*w*w) - f2*beta/(w*w) - beta*s*s*f2*f2;
	
	f_prec = f;
	f = (f1+f2)*(f1+Q*f2) - constte*(beta*beta)/pow(u*w, 4);
	fprime = (f1prime+f2prime)*(f1+Q*f2)+(f1+f2)*(f1prime+Q*f2prime)-2*constte*beta/pow(u*w, 4)
	  + 4*constte*pow(beta, 3)*(u*u - w*w)/pow(u*w, 6);
	
	beta = beta - f/fprime;          
        
	k = k+1;
      }
    
    if ((abs(f) > 1e-6) || std::isnan(f))
      {
	cout << "No convergence for omega = " << omega << endl;
	DISP(n1); DISP(n2); DISP(a); DISP(f);
	abort();
      }
    
    return beta;
  } 
  
  
  Real_wp OpticalFiber
  ::GetAeff(int m, const Real_wp& beta, const Real_wp& omega, const Real_wp& c0)
  {
    Real_wp k0 = omega/c0;
    Real_wp Kappa = k0*sqrt(n1*n1-beta*beta);
    Real_wp Gamma = k0*sqrt(beta*beta - n2*n2);
    
    Real_wp C = GetJn(0, Kappa*a) / GetKn(0, Gamma*a);
    
    Real_wp integrale1 = 0;
    Real_wp integrale2 = 0;
    Real_wp integrale3 = 0;
    
    Real_wp numer = 2.0*pi_wp/square(Kappa)
      *(0.5*a*Kappa*(square(GetJn(0, Kappa*a)) + square(GetJn(1, Kappa*a))));
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(12);
    
    // integral in the intern disc
    Real_wp x, poids;
    for (int k = 0; k <= gauss.GetOrder(); k++)
      {
	x = Kappa*a*gauss.Points(k);
	poids = Kappa*a*gauss.Weights(k);
	integrale2 += poids*x*pow(GetJn(0, x), 4);
      }
    
    // outside integrals
    int num = 1;
    for (int n = 0; n < 8; n++)
      {
	Real_wp xmin = Gamma*a*num, xmax = 2.0*xmin;
	Real_wp length_interval = xmax-xmin;
	for (int k = 0; k <= gauss.GetOrder(); k++)
	  {
	    x = xmin + gauss.Points(k)*(xmax-xmin);
	    poids = gauss.Weights(k)*length_interval;
            
	    integrale1 += poids*x*square(GetKn(0, x));
	    integrale3 += poids*x*pow(GetKn(0, x), 4);
	  }
	
	num *= 2;
      }
    
    numer += 2.0*pi_wp/square(Gamma)*square(C)*integrale1;
    
    Real_wp denom = 2.0*pi_wp/square(Kappa)*integrale2;
    denom += 2.0*pi_wp/square(Gamma)*pow(C, 4)*integrale3;
    
    Real_wp Aeff = square(numer) / denom;
    return Aeff;
  }
  
  
  Real_wp OpticalFiber
  ::GetNomega(int m, const Real_wp& beta, const Real_wp& omega, const Real_wp& c0)
  {
    Real_wp k0 = omega/c0;
    Real_wp Kappa = k0*sqrt(n1*n1-beta*beta);
    Real_wp Gamma = k0*sqrt(beta*beta - n2*n2);
    
    Real_wp Btilde = Gamma*(GetJn(0, Kappa*a) - GetJn(2, Kappa*a))*GetKn(1, Gamma*a)*n1*n1; 
    Btilde = Btilde - Kappa*GetJn(1, Kappa*a)*(GetKn(0, Gamma*a) + GetKn(2, Gamma*a))*n1*n1; 
    Btilde = 0.5*a*Kappa*Gamma*Btilde;
    // Btilde est le quotient entre le champ B et le champ E;
    //quo = omega*(n1*n1 -n2*n2)*beta*j1(Kappa*a)*GetKn(1, Gamma*a) ;
    
    Real_wp coeff1 = square(Btilde) + n1*n1 ;
    coeff1 = coeff1*square(Kappa);
    coeff1 = coeff1*k0;
    coeff1 = coeff1/pow(Kappa, 4);
    
    Real_wp coeff2 = square(Btilde) + n1*n1;
    coeff2 = coeff2*beta*k0;
    
    Real_wp coeff3 = Kappa*square(k0);
    coeff3 = coeff3*n1*n1;
    coeff3 = coeff3*Btilde;
    coeff3 = coeff3 + Kappa*Btilde*square(beta);
    coeff3 = -2*coeff3/pow(Kappa, 4);
    
    //  Les 4 integrales intervenant dans le calcul de la puissance traversant le coeur
    Real_wp integrale1 = 1.5*(square(GetJn(0, Kappa*a)) + square(GetJn(1, Kappa*a)))*a*a;
    integrale1 = integrale1 + 0.5*(square(GetJn(2, Kappa*a)) - GetJn(1, Kappa*a)*GetJn(3, Kappa*a))*a*a;
    integrale1 = integrale1 - 4*(1-square(GetJn(0, Kappa*a)))/square(Kappa);
    integrale1 = 0.25*integrale1;
    
    Real_wp integrale2 = 0.5*(1-square(GetJn(0, Kappa*a)) - square(GetJn(1, Kappa*a)));
    Real_wp integrale3 = 0.5*(square(GetJn(1, Kappa*a)) -1)/Kappa;
    
    Real_wp NmOmega = coeff1*integrale1 + coeff2*integrale2 + coeff3*integrale3;
    NmOmega = 2.0*pi_wp*NmOmega;
    
    // On calcule maintenant la puissance traversant la gaine      
    Real_wp Ctilde = jn(m, Kappa*a)/GetKn(m, Gamma*a);
    Real_wp Dtilde = Btilde*jn(m, Kappa*a)/GetKn(m, Gamma*a);
    
    Real_wp coeff1bis = square(Ctilde*n2)*k0*beta*square(Gamma);
    coeff1bis = coeff1bis + k0*beta*square(Dtilde*Gamma);
    coeff1bis = coeff1bis/pow(Gamma, 4);
    
    Real_wp coeff2bis = beta*k0*Dtilde*Dtilde;
    coeff2bis = coeff2bis + beta*k0*square(Ctilde*n2);
    coeff2bis = coeff2bis*m*m/pow(Gamma, 4);
    
    Real_wp coeff3bis = -2*beta*beta*Gamma*Ctilde*Dtilde;
    coeff3bis = coeff3bis - Gamma*k0*Ctilde*Dtilde*n2*n2;
    coeff3bis = coeff3bis*m/pow(Gamma, 4);
    
    Real_wp integrale1bis(0), integrale2bis(0), integrale3bis(0);
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(10);
    
    int num = 1;
    Real_wp x, feval, poids;
    for (int n = 0; n < 7; n++)
      {
	Real_wp xmin = Gamma*a*num, xmax = 2.0*xmin;
	Real_wp length_interval = xmax-xmin;
	for (int k = 0; k <= gauss.GetOrder(); k++)
	  {
	    x = xmin + gauss.Points(k)*(xmax-xmin);
	    poids = gauss.Weights(k)*length_interval;
	    feval = x*square(GetKn(m+1, x) + GetKn(m-1, x));
	    integrale1bis += poids*feval;
            
	    feval = square(GetKn(m, x))/x;
	    integrale2bis += poids*feval;              
	  }
	
	num *= 2;
      }
    
    integrale1bis = integrale1bis/(Gamma*Gamma);
    integrale3bis = -0.5*square(GetKn(m, Gamma*a))/(2.0*Gamma);
    
    NmOmega = NmOmega + 2*pi_wp*(coeff1bis*integrale1bis +coeff2bis*integrale2bis + coeff3bis*integrale3bis);
    
    return NmOmega;
  }
  
}

#define MONTJOIE_FILE_NON_LINEAR_OPTICS_PROBLEM_CXX
#endif
