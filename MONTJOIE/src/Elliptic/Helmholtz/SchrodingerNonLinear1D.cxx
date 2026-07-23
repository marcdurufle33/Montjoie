#ifndef MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_1D_CXX

namespace Montjoie
{
  //! default constructor
  KerrProblem::KerrProblem()
  {
    nb_points_time = 1;
    padding_time = 0;
    type_forcing_number_fft_points = POWER_NONE;    
    
    gammaNL = 0.0;
    
    omega_L = 0.0; k_L = 0.0;
    group_velocity = 1.0;
    
    tau_impulse = 1.0;
    amplitude_impulse = 1.0;
    Tmax = 5.0;    
    
    dz = 1.0;
    zmax = 1.0;
    
    z0_adim = 1.0;
    E0_adim = 1.0;
    t0_adim = 1.0;
    coef_fft = 1.0;
    coef_fft_time = 1.0;
    coef_convolution = 1.0;
    
    // default directory is the current directory
    DOSSIER = "./";
    
    type_model = SCHRODINGER_CONSTANT;
    type_polarization = CIRCULAR;
    polarization.Init(Complex_wp(1, 0), Complex_wp(0, 0));
    
    type_spectrum = SHIFTED_SINGLE;
    number_odd_frequencies = 1;
    
    nb_iterations_between_two_tests = 10;
    check_overflow_scheme = true;
    factor_refinement_time = 1;
    size_refinement_output = 0;
    threshold_time_extraction = 0;
    threshold_last_frequency = 1e-4;
    threshold_negative_epsilon = 1e-6;
    affiche_enveloppe = false;
    coef_fft_envelope_time = 1.0;
    
    local_enveloppe = false;
    z0_current = 0.0; threshold_phase_imag = 1e-100;
    
    epsilon_min_treated = -1e300; epsilon_max_treated = 1e300;
    coef_stability_damping = 0.0;

    print_level = 0;    
    fiber_computation = false;
  }
  
  
  //! returns the number of points in time
  int KerrProblem::GetNbPointsTime() const
  {
    return nb_points_time;
  }
  
  
  //! modifying parameters with a line of the datafile
  void KerrProblem::SetInputData(const string& keyword, const Vector<string>& parameters)
  {
    if (keyword == "WavelengthAdim")
      {
	// WavelengthAdim = z0
        if (parameters.GetM() <= 0)
          {
	    cout << "In SetInputData of KerrProblem" << endl;
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
    else if (keyword == "DissipativeCoefficient")
      {
        // DissipativeCoefficient = sigma (conductivity)
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "DissipativeCoefficient needs one parameter, for instance :" << endl;
	    cout << "DissipativeCoefficient = sigma" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }
        
        Real_wp sigma = to_num<Real_wp>(parameters(0));
	index.SetSigma(sigma, z0_adim);
      }
    else if (keyword == "NonLinearSusceptibility")
      {
	// NonLinearSusceptibility = gammaNL
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "NonLinearSusceptibility needs one parameter, for instance :" << endl;
	    cout << "NonLinearSusceptibility = gammaNL" << endl;
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
    else if (keyword == "RefractiveIndex")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
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
    else if (keyword == "Fiber")
      {
        if (parameters.GetM() <= 3)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "Fiber needs more parameters" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        fiber_computation = true;
        if (parameters(0) == "Circle")
          {
            Real_wp a = to_num<Real_wp>(parameters(1));
            Real_wp n1 = to_num<Real_wp>(parameters(2));
            Real_wp n2 = to_num<Real_wp>(parameters(3));
            a /= z0_adim;
            fiber.SetCircularFiber(a, n1, n2);
          }
      }
    else if (keyword == "DisplayZ")
      {
	// DisplayZ = z0 z1 dz_display
        output_grid_z.SetInputData(parameters);
	output_grid_z.AdimTime(z0_adim);
      }
    else if (keyword == "FileOutputDisplayZ")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "FileOutputDisplayZ needs one parameter, for instance :" << endl;
	    cout << "FileOutputDisplayZ = file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

	// FileOutputDisplayZ = name_file
        output_grid_z.SetTotalFieldFile(parameters(0));
      }
    else if (keyword == "PropagationStep")
      {
	// PropagationStep = dz
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "PropagationStep needs one parameter, for instance :" << endl;
	    cout << "PropagationStep = dz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        dz = to_num<Real_wp>(parameters(0))/z0_adim;
      }
    else if (keyword == "PropagationInterval")
      {
	// PropagationInterval = zmax
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "PropagationInterval needs one parameter, for instance :" << endl;
	    cout << "PropagationInterval = zmax" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        zmax = to_num<Real_wp>(parameters(0))/z0_adim;
      }
    else if (keyword == "Frequency")
      {
	// Frequency = f_L
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "Frequency needs one parameter, for instance :" << endl;
	    cout << "Frequency = f" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        Real_wp f = to_num<Real_wp>(parameters(0));
        SetFrequency(f*t0_adim);
      }
    else if (keyword == "Wavelength")
      {
	// wavelength in the vacuum
	// Wavelength = lambda
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
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
    else if (keyword == "SchrodingerModel")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
	    cout << "Wavelength needs one parameter, for instance :" << endl;
	    cout << "Wavelength = file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

	if (parameters(0) == "Constant")
	  type_model = SCHRODINGER_CONSTANT;
	else if (parameters(0) == "UPPE")
	  type_model = UPPE;
	else if (parameters(0) == "FOP")
	  type_model = FOP;
	else if (parameters(0) == "NEE")
	  type_model = NEE;
	else if (parameters(0) == "NEEConstant")
	  type_model = NEEC;
	else if (parameters(0) == "NLS")
	  type_model = NLS;
	else
	  type_model = MIRO;
      }
    else if (keyword == "TimeImpulsion")
      {
	// TimeImpulsion = N Npadding Gaussian tau Tmax amplitude
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of KerrProblem" << endl;
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
                cout << "In SetInputData of KerrProblem" << endl;
                cout << "TimeImpulsion needs six parameters, for instance :" << endl;
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
            amplitude_impulse = to_num<Real_wp>(parameters(5));
	    
	    DISP(tau_impulse);
	    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
	      {
		tau_impulse /= t0_adim;
		Tmax /= t0_adim;
		amplitude_impulse /= E0_adim;
	      }
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
    else if (keyword == "Polarization")
      {	
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "Polarization needs at least one parameter, for instance :" << endl;
            cout << "Polarization = type_polarization" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        // Polarization = Linear
        // Polarization = Circular
        // Polarization = Elliptic Ex Ey phix phiy
        // type of polarization
	R2 phase_init, polarization_init;
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
                cout << "In SetInputData of KerrProblem" << endl;
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
	
	polarization(0) = polarization_init(0)*exp(-Iwp*phase_init(0));
	polarization(1) = polarization_init(1)*exp(-Iwp*phase_init(1));
      }
    else if (keyword == "FourierTransform")
      {
        // FourierTransform = TIME_CENTER
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "FourierTransform needs at least one parameter, for instance :" << endl;
            cout << "FourierTransform = TIME_CENTER" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        if (parameters(0) == "TIME_CENTER")
          type_spectrum = CENTERED;
        else if (parameters(0) == "ODD_FREQUENCIES")
          {
            if (parameters.GetM() <= 1)
              {
                cout << "In SetInputData of KerrProblem" << endl;
                cout << "FourierTransform needs two parameter, for instance :" << endl;
                cout << "FourierTransform = ODD_FREQUENCIES p" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }

            type_spectrum = ODD_FREQUENCIES;
            number_odd_frequencies = to_num<int>(parameters(1));
          }
        else
          type_spectrum = SHIFTED_SINGLE;
      }
    else if (keyword == "ForbidNegativeIndex")
      index.forbid_negative_index = true;
    else if (keyword == "ForceNiceNumberPoints")
      {
        // ForceNiceNumberPoints = 2357
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "ForceNiceNumberPoints needs one parameter, for instance :" << endl;
            cout << "ForceNiceNumberPoints = 235" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        if (parameters(0) == "2")
          type_forcing_number_fft_points = POWER_TWO;
        else if (parameters(0) == "23")
          type_forcing_number_fft_points = POWER_TWO_THREE;
        else if (parameters(0) == "235")
          type_forcing_number_fft_points = POWER_TWO_THREE_FIVE;
        else if (parameters(0) == "2357")
          type_forcing_number_fft_points = POWER_TWO_THREE_FIVE_SEVEN;
        else
          type_forcing_number_fft_points = POWER_NONE;
      }
    else if (keyword == "OrderTimeScheme")
      time_scheme.SetInputData(keyword, parameters);
    else if (keyword == "TimeRefinement")
      {
        // TimeRefinement = r threshold size_min
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "TimeRefinement needs at least one parameter, for instance :" << endl;
            cout << "TimeRefinement = factor" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        factor_refinement_time = to_num<int>(parameters(0));
        if (parameters.GetM() >= 2)
          threshold_time_extraction = to_num<Real_wp>(parameters(1));
        if (parameters.GetM() >= 3)
          size_refinement_output = to_num<int>(parameters(2));
      }
    else if (keyword == "AfficheEnveloppe")
      {
        // AfficheEnveloppe = YES
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "AfficheEnveloppe needs one parameter, for instance :" << endl;
            cout << "AfficheEnveloppe = YES" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        if (parameters(0) == "YES")
          affiche_enveloppe = true;
        else
          affiche_enveloppe = false;
      }
    else if (keyword == "IntervalEpsilon")
      {
        // IntervalEpsilon = eps_min eps_max
        if (parameters.GetM() <= 1)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "IntervalEpsilon needs two parameters, for instance :" << endl;
            cout << "IntervalEpsilon = eps_min eps_max" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        epsilon_min_treated = to_num<Real_wp>(parameters(0));
        epsilon_max_treated = to_num<Real_wp>(parameters(1));
      }
    else if (keyword == "CoefficientDampingCFL")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "CoefficientDampingCFL needs one parameter, for instance :" << endl;
            cout << "CoefficientDampingCFL = CFL" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        if (parameters(0) == "LOCAL")
	  {
	    local_enveloppe = true;
	    if (parameters.GetM() > 1)
	      threshold_phase_imag = to_num<Real_wp>(parameters(1));
	    else
	      threshold_phase_imag = 1e-100;
	  }
	else
	  coef_stability_damping = to_num<Real_wp>(parameters(0));	    
      }
    else if (keyword == "PrintLevel")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of KerrProblem" << endl;
            cout << "PrintLevel needs one parameter, for instance :" << endl;
            cout << "PrintLevel = 2" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        print_level = to_num<int>(parameters(0));
        if (parameters.GetM() > 1)
          nb_iterations_between_two_tests = to_num<int>(parameters(1));
      }
    else if (keyword == "DirectoryOutput")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of KerrProblem" << endl;
	    cout << "DirectoryOutput needs one parameter, for instance :" << endl;
	    cout << "DirectoryOutput = chemin" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	if (parameters(0) == "AUTO")
          DOSSIER.clear();
        else
          DOSSIER = parameters(0);
      }
  }
  
  
  //! setting the frequency
  void KerrProblem::SetFrequency(const Real_wp& f)
  {
    omega_L = 2.0*pi_wp*f;
  }

  //! constructs arrays needed for time iterations
  void KerrProblem::ConstructAll()
  {
    if (!fiber_computation)
      k_L = index.GetWaveNumber(omega_L);
    else
      {
        k_L = fiber.GetKomega(1, omega_L, index.c0);
	k_L *= omega_L/index.c0;
      }
    
    group_velocity = index.GetGroupVelocity(omega_L);
    
    // for faster fft, multiples of 2, 3, 5 and 7 can be forced
    if (type_forcing_number_fft_points == POWER_TWO)
      {
	nb_points_time = FindClosestPow2(nb_points_time);
	padding_time = FindClosestPow2(padding_time);
	factor_refinement_time = FindClosestPow2(factor_refinement_time);
	size_refinement_output = FindClosestPow2(size_refinement_output);
      }
    else if (type_forcing_number_fft_points == POWER_TWO_THREE)
      {
	nb_points_time = FindClosestPow23(nb_points_time);
	padding_time = FindClosestPow23(padding_time);
	factor_refinement_time = FindClosestPow23(factor_refinement_time);
	size_refinement_output = FindClosestPow23(size_refinement_output);
      }
    else if (type_forcing_number_fft_points == POWER_TWO_THREE_FIVE)
      {
	nb_points_time = FindClosestPow235(nb_points_time);
	padding_time = FindClosestPow235(padding_time);
	factor_refinement_time = FindClosestPow235(factor_refinement_time);
	size_refinement_output = FindClosestPow235(size_refinement_output);
      }
    else if (type_forcing_number_fft_points == POWER_TWO_THREE_FIVE_SEVEN)
      {
	nb_points_time = FindClosestPow2357(nb_points_time);
	padding_time = FindClosestPow2357(padding_time);
	factor_refinement_time = FindClosestPow2357(factor_refinement_time);
	size_refinement_output = FindClosestPow2357(size_refinement_output);
      }
    
    if (print_level >= 0)
      {
	cout << "Effective number of points = " << nb_points_time << endl;
	cout << "Padding (unused) = " << padding_time << endl;
	cout << "Refinement factor = " << factor_refinement_time << endl;
	cout << "Number of extracted points = " << size_refinement_output << endl;
      }

    // impulsion is defined on interval [-Tmax, Tmax]
    omega.Reallocate(nb_points_time);
    kwave.Reallocate(nb_points_time);
    if (fiber_computation)
      Aeff_om.Reallocate(nb_points_time);
    
    time_interval.Reallocate(nb_points_time);
    Real_wp dt = 2.0*Tmax/Real_wp(nb_points_time);
    time_interval.Fill(0);
    for (int i = 0; i < nb_points_time/2; i++)
      {
        time_interval(i) = Real_wp(i)*dt;
        time_interval(i+nb_points_time/2) = Real_wp(i-nb_points_time/2)*dt;
      }
    
    // calcul de l'indice pour une gamme de longueurs d'ondes
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
	omega.Fill(omega_L); kwave.Fill(k_L);
	for (int i = 0; i < nb_points_time/2; i++)
	  {
	    omega(i) = omega_L + Omega(nb_points_time/2+i);
	    if (!fiber_computation)
              kwave(i) = index.GetWaveNumber(omega(i));
            else
              {
		Real_wp beta = fiber.GetKomega(1, omega(i), index.c0);
                kwave(i) = beta;
                kwave(i) *= omega(i)/index.c0;
		
		Real_wp Aeff = fiber.GetAeff(1, beta, omega(i), index.c0);
		Aeff_om(i) = Aeff*square(z0_adim);
              }
	    
	    Real_wp omega_ = omega_L + Omega(i);
	    if (omega_ < 0)
	      {
		DISP(omega_L); DISP(dt); DISP(Omega.GetNormInf()); 
		cout << "The spectrum is assumed to be around the value " << omega_L << endl;
		cout << "Therefore the pulsation omega should be always positive" << endl;
		cout << "But it is equal to " << omega_ << endl;
		abort();
	      }
	    
	    omega(i+nb_points_time/2) = omega_;
	    if (!fiber_computation)
              kwave(i+nb_points_time/2) = index.GetWaveNumber(omega_);
            else
              {
		Real_wp beta = fiber.GetKomega(1, omega_, index.c0);	    
                kwave(i+nb_points_time/2) = beta;
                kwave(i+nb_points_time/2) *= omega_/index.c0;

		Real_wp Aeff = fiber.GetAeff(1, beta, omega_, index.c0);
		Aeff_om(i+nb_points_time/2) = Aeff*square(z0_adim);
              }
	  }
	
	// writing pulsations and wave numbers on the disk
	VectReal_wp omega_display(nb_points_time);
	VectReal_wp Aeff_display(nb_points_time);
	VectComplex_wp kwave_display(nb_points_time);
	for (int i = 0; i < nb_points_time/2; i++)
	  {
	    omega_display(nb_points_time/2+i) = omega(i);
	    omega_display(i) = omega(nb_points_time/2+i);
	    
	    kwave_display(nb_points_time/2+i) = kwave(i);
	    kwave_display(i) = kwave(nb_points_time/2+i);
	    
	    if (fiber_computation)
	      {
		Aeff_display(nb_points_time/2+i) = Aeff_om(i);
		Aeff_display(i) = Aeff_om(nb_points_time/2+i);
	      }
	  }
	
	Mlt(1.0/t0_adim, omega_display);
	Mlt(1.0/z0_adim, kwave_display);
	
	omega_display.Write(DOSSIER + "omega.dat");
	kwave_display.Write(DOSSIER + "kwave.dat");
	if (fiber_computation)
	  Aeff_display.Write(DOSSIER + "Aeff.dat");
	
	// initializing FFT
	fft.Init(nb_points_time + padding_time);
	coef_fft = 2.0*Tmax/sqrt(2.0*pi_wp);

        fft_time.Init(nb_points_time + padding_time);
        coef_fft_time = 4.0*Tmax/(nb_points_time+padding_time)/sqrt(2.0*pi_wp);

        // refinement in time ?
        if (factor_refinement_time > 1)
          {
            int Nt = nb_points_time*factor_refinement_time;
            fft_refinement.Init(Nt);
            coef_fft_refinement = 2.0*Tmax/sqrt(2.0*pi_wp);
            time_interval.Reallocate(Nt);
            dt = 2.0*Tmax/Real_wp(Nt);
            time_interval.Fill(0);
            for (int i = 0; i < Nt/2; i++)
              {
                time_interval(i) = Real_wp(i)*dt;
                time_interval(i+Nt/2) = Real_wp(i-Nt/2)*dt;
              }
          }

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
        
        omega_display.Write(DOSSIER + "omega.dat");
        kwave_display.Write(DOSSIER + "kwave.dat");    
        
        // initializing FFT
        fft.Init(nb_points_time + padding_time);
        coef_fft = 2.0*Tmax/sqrt(2.0*pi_wp);
      }
    else
      {
	omega.Reallocate(nb_points_time/2+1);
	kwave.Reallocate(nb_points_time/2+1);
	Aeff_om.Reallocate(nb_points_time/2+1);
        for (int i = 1; i <= nb_points_time/2; i++)
          {
            omega(i) = Omega(nb_points_time/2+i);
            if (!fiber_computation)
	      kwave(i) = index.GetWaveNumber(omega(i));
	    else
	      {
		Real_wp beta = fiber.GetKomega(1, omega(i), index.c0);
                kwave(i) = beta;
                kwave(i) *= omega(i)/index.c0;
		
		Real_wp Aeff = fiber.GetAeff(1, beta, omega(i), index.c0);
		Aeff_om(i) = Aeff*square(z0_adim);
	      }
	  }
	
	omega(0) = 0; kwave(0) = 0;
	
        // writing pulsations and wave numbers on the disk
        VectReal_wp omega_display(omega);
        VectComplex_wp kwave_display(kwave);
        
        Mlt(1.0/t0_adim, omega_display);
        Mlt(1.0/z0_adim, kwave_display);
	
	VectReal_wp Aeff_display;
	if (fiber_computation)
	  Aeff_display = Aeff_om;
        
        omega_display.Write(DOSSIER + "omega.dat");
        kwave_display.Write(DOSSIER + "kwave.dat");    
	if (fiber_computation)
	  Aeff_display.Write(DOSSIER + "Aeff.dat");

        // initializing FFT
        fft_time.Init(nb_points_time + padding_time);
        coef_fft = 4.0*Tmax/(nb_points_time+padding_time)/sqrt(2.0*pi_wp);
	coef_fft_time = coef_fft;
	
        // refinement in time ?
        int Nt = nb_points_time;
        if (factor_refinement_time > 1)
          {
            Nt = nb_points_time*factor_refinement_time;
            fft_refinement_time.Init(Nt);
            coef_fft_refinement = 4.0*Tmax/Nt/sqrt(2.0*pi_wp);
            time_interval.Reallocate(Nt);
            dt = 2.0*Tmax/Real_wp(Nt);
            time_interval.Fill(0);
            for (int i = 0; i < Nt/2; i++)
              {
                time_interval(i) = Real_wp(i)*dt;
                time_interval(i+Nt/2) = Real_wp(i-Nt/2)*dt;
              }
          }
        
        // fft to obtain the enveloppe 
        if (affiche_enveloppe)
          {
            fft_envelope_time.Init(Nt);
            coef_fft_envelope_time = 2.0*Tmax/Nt/sqrt(2.0*pi_wp);
          }
      }
    
    // Fourier transform of Raman function
    if (var_raman.GetAlpha() != Real_wp(0))
      {
	VectReal_wp raman_time;
	raman_time.Reallocate(nb_points_time);
	raman_fct_hat.Reallocate(nb_points_time/2+1);
	
	for (int i = 0; i < nb_points_time; i++)
	  raman_time(i) = var_raman.EvaluateH(time_interval(i));
	
	time_interval.Write("time.dat");
	raman_time.Write("raman.dat");
	fft_time.ApplyForward(raman_time, raman_fct_hat);
	Mlt(coef_fft_time, raman_fct_hat);
	
	for (int i = 0; i < raman_fct_hat.GetM(); i++)
	  raman_fct_hat(i) = conj(raman_fct_hat(i));
	
	coef_convolution = sqrt(pi_wp/2);
	Mlt(coef_convolution, raman_fct_hat);

	raman_fct_hat.Write("raman_hat.dat");
      }
    
    // appel de EvaluateFunction pour initialiser les tableaux restants
    VectComplex_wp u;
    GetInitialCondition(u);
    VectComplex_wp scheme(u);
    EvaluateFunction(0.0, u, scheme); 
    
    // on regarde si epsilon est dans le bon intervalle
    kwave_to_consider.Reallocate(kwaveAux.GetM());
    kwave_to_consider.Fill(true);
    n2omega.Reallocate(kwaveAux.GetM());
    n2omega.Fill(0);
    for (int i = 0; i < kwaveAux.GetM(); i++)
      if (omegaAux(i) != Real_wp(0))
        {
          Complex_wp eps = square(kwaveAux(i)*index.c0/omegaAux(i));
          if ((real(eps) < epsilon_min_treated) || (real(eps) > epsilon_max_treated))
            kwave_to_consider(i) = false;
          else
            kwave_to_consider(i) = true;
	  
          n2omega(i) = real(eps);
	  if (local_enveloppe)
	    {
	      Real_wp phase_real, inv_phase_real;
	      GetPhaseImagKz(kwaveAux(i), dz, phase_real, inv_phase_real);
	      if (phase_real == Real_wp(0))
		kwave_to_consider(i) = false;
	    }
        }    
    
    cout << "The group velocity is equal to " << group_velocity*(z0_adim/t0_adim) << endl;
  }
  
  
  //! input parameters are copied from 1-D maxwell's problem
  void KerrProblem::Init(MaxwellProblem1D& var)
  {
    nb_points_time = var.nb_points_time;
    padding_time = var.padding_time;
    gammaNL = var.gammaNL;
    omega = var.omega;
    kwave = var.kwave;
    omega_L = var.omega_L;
    k_L = var.index.GetWaveNumber(var.omega_L);
    index.alpha_damping = var.index.alpha_damping;
    tau_impulse = var.tau_impulse;
    Tmax = var.Tmax;
    amplitude_impulse = var.amplitude_impulse;
    time_interval = var.time_interval;
    index.c0 = var.index.c0;
    dz = 0;
    zmax = 0;
    index.omega_polarization = var.index.omega_polarization;
    index.sigma_polarization = var.index.sigma_polarization;
    index.alpha_polarization = var.index.alpha_polarization;
    index.eps_inf = var.index.eps_inf;
    z0_adim = var.z0_adim;
    E0_adim = var.E0_adim;
    t0_adim = var.t0_adim;
    coef_fft = var.coef_fft;
    type_model = UPPE;
    type_polarization = var.type_polarization;
    polarization(0) = var.polarization_init(0)*exp(-Iwp*var.phase_init(0));
    polarization(1) = var.polarization_init(1)*exp(-Iwp*var.phase_init(1));
    type_spectrum = var.type_spectrum;
    number_odd_frequencies = var.number_odd_frequencies;
    time_scheme.SetInputData(string("OrderTimeScheme"), var.param_time_scheme);
    print_level = 0;
    
    fft.Init(nb_points_time+padding_time);
    coef_fft = 2.0*Tmax/sqrt(2.0*pi_wp);
  }
  
  
  //! returns exp(i real(k) z) 
  Complex_wp KerrProblem::GetPhaseKz(const Complex_wp& k, const Real_wp& z)
  {
    return Complex_wp(cos(real(k)*z), sin(real(k)*z));
  }
  
  
  //! computes exp(- imag(k) z) and its inverse
  void KerrProblem::GetPhaseImagKz(const Complex_wp& k, const Real_wp& z,
				   Real_wp& phase_real, Real_wp& inv_phase_real)
  {
    //phase_real = 1.0; inv_phase_real = 1.0;
    //return;
    
    phase_real = exp(-imag(k)*z);
    if (phase_real <= threshold_phase_imag)
      {
	phase_real = 0.0;
	inv_phase_real = 0.0;
      }
    else
      inv_phase_real = 1.0/phase_real;
  }
  
  
  //! computes phase for a given argument that may already be computed
  void KerrProblem::GetLocalPhase(const Vector<Complex_wp>& kw, const Real_wp& z,
				  Vector<Complex_wp>& phase, VectReal_wp& phase_real,
				  VectReal_wp& inv_phase_real)
  {
    int n = -1;
    for (int p = 0; p < arg_phaseZ_stored.GetM(); p++)
      if (abs(arg_phaseZ_stored(p) - z) < 1e-8*dz)
	n = p;
    
    if (n == -1)
      {
	// new argument
        #pragma omp parallel for
	for (int k = 0; k < kw.GetM(); k++)	      
	  {
	    phase(k) = GetPhaseKz(kw(k), z);
	    GetPhaseImagKz(kw(k), z, phase_real(k), inv_phase_real(k));
	  }
	
	arg_phaseZ_stored.PushBack(z);
	phaseStored.PushBack(phase);
	phaseRealStored.PushBack(phase_real);
	invPhaseRealStored.PushBack(inv_phase_real);
      }
    else
      {
        #pragma omp parallel for
	for (int k = 0; k < kw.GetM(); k++)	      
	  {
	    //phase(k) = GetPhaseKz(kw(k), z);
	    //GetPhaseImagKz(kw(k), z, phaseImag(k), inv_phase_real(k));
	    
	    phase(k) = phaseStored(n)(k);
	    phase_real(k) = phaseRealStored(n)(k);
	    inv_phase_real(k) = invPhaseRealStored(n)(k);
	  }
      }
  }
  

  //! Applies Raman effect to |E|^2
  template<class T>
  void KerrProblem::ApplyRamanEffect(const VectReal_wp& moduleU2,
				     Vector<T>& utime, Real_wp coef_corr)
  {
    Vector<T> utime_y;
    ApplyRamanEffect(moduleU2, utime, utime_y, coef_corr);
  }
  
  
  //! Applies Raman effect to |E|^2
  template<class T>
  void KerrProblem::ApplyRamanEffect(const VectReal_wp& moduleU2,
				     Vector<T>& utime, Vector<T>& utime_y,
				     Real_wp coef_corr)
  {
    if (this->var_raman.GetAlpha() == Real_wp(0))
      {
	// instantaneous Kerr effect
	// forming |E|^2 E
	#pragma omp parallel for
	for (int j = 0; j < utime.GetM(); j++)
	  {
	    utime(j) *= moduleU2(j);
	    if (utime_y.GetM() > 0)
	      utime_y(j) *= moduleU2(j);
	  }
      }
    else
      {
	VectComplex_wp moduleU2_chap(nb_points_time/2+1);
	VectReal_wp moduleU2_conv(moduleU2.GetM());
	// R = (1-\alpha) \delta + \alpha h_R
	// treating convolution with h_R
	fft_time.ApplyForward(moduleU2, moduleU2_chap);
	
	#pragma omp parallel for
	for (int i = 0; i < moduleU2_chap.GetM(); i++)
	  moduleU2_chap(i) = moduleU2_chap(i)*conj(raman_fct_hat(i));
	
	fft_time.ApplyInverse(moduleU2_chap, moduleU2_conv);
	
	// then forming (1-alpha) |E|^2 E + alpha (h_R * |E|^2) E
	Real_wp coef_conv = this->var_raman.GetAlpha(), coef_kerr = 1.0-coef_conv;
	coef_conv *= coef_corr;
	Real_wp coef_g;
	#pragma omp parallel for
	for (int j = 0; j < nb_points_time; j++)
	  {
	    coef_g = coef_kerr*moduleU2(j) + coef_conv*moduleU2_conv(j);
	    utime(j) *= coef_g;
	    if (utime_y.GetM() > 0)
	      utime_y(j) *= coef_g;
	  }
      }
  }
  

  //! evaluation of scheme = g(u) where the evolution system is du/dz = g(u)
  void KerrProblem::EvaluateFunction(const Real_wp& zn, const VectComplex_wp& u,
				     VectComplex_wp& scheme, bool invert_mass, bool source)
  {
    if (print_level >= 4)
      glob_chrono.Start(VirtualTimer::JACOBIAN);
    
    // on calcule u en temps (avec du zero-padding)
    int N = u.GetM();
    Complex_wp one, zero;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    //on les utilise pour éviter les répétitions de code 
    // (à cause de la polarisation) dans la partie du choix du modele
    bool init = false;
    if (omegaAux.GetM() == 0)
      {
        init = true;
        omegaAux.Reallocate(N); omegaAux.Fill(0);
        kwaveAux.Reallocate(N); kwaveAux.Fill(zero);        
	Aeff_aux.Reallocate(N); Aeff_aux.Fill(1);
      }

    // on calcule la polarization non lineaire (P_NL)
    if (type_polarization == ELLIPTIC)
      {
	Vector<Complex_wp> phase(N/2);
	Vector<Real_wp> phaseImag(N/2), invPhaseImag(N/2);

        if (init)
          for (int i = 0; i < N/2; i++)
            {
              omegaAux(i) = omega(i);
              omegaAux(N/2+i) = omega(i);
              kwaveAux(i) = kwave(i);
              kwaveAux(N/2+i) = kwave(i);
	      if (fiber_computation)
		{
		  Aeff_aux(i) = Aeff_om(i);
		  Aeff_aux(N/2+i) = Aeff_om(i);
		}
            }
        
        if (print_level >= 4)
          glob_chrono.Start(VirtualTimer::PML);

	phaseImag.Fill(1.0); invPhaseImag.Fill(1.0);
	
	if (NonLinearPhase())
	  {
	    if (local_enveloppe)
	      GetLocalPhase(kwave, zn-z0_current, phase, phaseImag, invPhaseImag);
	    else
	      {
                #pragma omp parallel for
                for (int k = 0; k < N/2 ; k++)	      
                  phase(k) = GetPhaseKz(kwave(k), zn);
              }
	  }
	else
	  {
	    if (local_enveloppe)
	      {
		cout << "Not implemented for this model " << endl;
		abort();
	      }
	    
	    // if the phase is affine, there is no need
	    // to consider it in the non-linear term
	    phase.Fill(one);
	  }

        if (print_level >= 4)
          glob_chrono.Stop(VirtualTimer::PML);
	
	if (type_spectrum == CENTERED)
	  {
            if (print_level >= 4)
              glob_chrono.Start(VirtualTimer::EXTRAPOL);

	    VectComplex_wp ux(omega.GetM()), uy(omega.GetM());
	    VectReal_wp utime_x(nb_points_time), utime_y(nb_points_time);
            
            #pragma omp parallel for
	    for (int k = 0; k < omega.GetM(); k++)
	      {
		ux(k) = conj(u(k)*phase(k))*phaseImag(k);
		uy(k) = conj(u(k+N/2)*phase(k))*phaseImag(k);
	      }
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::EXTRAPOL);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    // returning back in time
	    fft_time.ApplyInverse(ux, utime_x);
	    fft_time.ApplyInverse(uy, utime_y);

            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::PROD);
              }

	    Mlt(1.0/coef_fft, utime_x);
	    Mlt(1.0/coef_fft, utime_y);
	    
	    // forming |E|^2 E
	    VectReal_wp moduleU2(nb_points_time);
            #pragma omp parallel for
	    for (int i = 0; i < nb_points_time; i++)
	      moduleU2(i) = square(utime_x(i)) + square(utime_y(i));
	    
	    this->ApplyRamanEffect(moduleU2, utime_x, utime_y);
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::PROD);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    // going fo Fourier domain
	    fft_time.ApplyForward(utime_x, ux);
	    fft_time.ApplyForward(utime_y, uy);
            
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::EXTRAPOL);
              }

	    Mlt(coef_fft, ux); Mlt(coef_fft, uy);
         
            #pragma omp parallel for
	    for (int k = 0; k < omega.GetM(); k++)
	      {
		// x-component
		scheme(k) = conj(ux(k)*phase(k))*invPhaseImag(k);
		scheme(k+N/2) = conj(uy(k)*phase(k))*invPhaseImag(k);
	      }

            if (print_level >= 4)
              glob_chrono.Stop(VirtualTimer::EXTRAPOL);
	  }
        else if (type_spectrum == ODD_FREQUENCIES)
          {
            if (print_level >= 4)
              glob_chrono.Start(VirtualTimer::FLUX);
            
            // not multithreaded
            int Nf = number_odd_frequencies;
            Vector<VectComplex_wp> VecEx(Nf), VecEy(Nf);
            VectComplex_wp PolarEx(Nf), PolarEy(Nf);
            for (int n = 0; n < Nf; n++)
              {
                VecEx(n).Reallocate(nb_points_time);
                VecEy(n).Reallocate(nb_points_time);
              }
            
            // extracting Ex, Ey and returning back in time
            int offset = 0;
            for (int n = 0; n < Nf; n++)
              {
                for (int k = 0; k < nb_points_time; k++)
                  {
                    VecEx(n)(k) = u(offset + k)*phase(offset+k)*phaseImag(offset+k);
                    VecEy(n)(k) = u(offset + k + N/2)*phase(offset+k)*phaseImag(offset+k);
                  }
                
                fft.ApplyForward(VecEx(n));
                fft.ApplyForward(VecEy(n));
                
                Mlt(1.0/coef_fft, VecEx(n));
                Mlt(1.0/coef_fft, VecEy(n));
                
                offset += nb_points_time;
              }
                
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::PROD);
              }

            // forming |E|^2 E with odd frequencies
            Real_wp coef_gamma = 1.0;
            Complex_wp Ek_dot_El, Ex_k, Ey_k, Ex_l, Ey_l;
            int Nmax = 2*Nf-1;
            for (int j = 0; j < nb_points_time; j++)
              {
                PolarEx.Fill(zero); PolarEy.Fill(zero);
                for (int k2 = -Nf; k2 < Nf; k2++)
                  for (int l2 = k2; l2 < Nf; l2++)
                    {
                      int k = 2*k2+1; int l = 2*l2+1;
                      if (k2 == l2)
                        coef_gamma = 0.25;
                      else
                        coef_gamma = 0.5;
                      
                      // forming coef_gamma Ek . El
                      if (k < 0)
                        {
                          Ex_k = conj(VecEx(-k2-1)(j));
                          Ey_k = conj(VecEy(-k2-1)(j));
                        }
                      else
                        {
                          Ex_k = VecEx(k2)(j);
                          Ey_k = VecEy(k2)(j);
                        }
                      
                      if (l < 0)
                        {
                          Ex_l = conj(VecEx(-l2-1)(j));
                          Ey_l = conj(VecEy(-l2-1)(j));
                        }
                      else
                        {
                          Ex_l = VecEx(l2)(j);
                          Ey_l = VecEy(l2)(j);
                        }
                      
                      Ek_dot_El = coef_gamma*(Ex_k*Ex_l + Ey_k*Ey_l);
                      
                      // then adding terms E_k \cdot E_l E_{n-k-l}
                      for (int n2 = 0; n2 < Nf; n2++)
                        {
                          int n = 2*n2+1;
                          int p = n-k-l;
                          
                          if (p < 0)
                            {
                              if (p >= -Nmax)
                                {
                                  int p2 = -(p+1)/2;
                                  PolarEx(n2) += Ek_dot_El*conj(VecEx(p2)(j));
                                  PolarEy(n2) += Ek_dot_El*conj(VecEy(p2)(j));
                                }
                            }
                          else
                            {
                              if (p <= Nmax)
                                {
                                  int p2 = (p-1)/2;
                                  PolarEx(n2) += Ek_dot_El*VecEx(p2)(j);
                                  PolarEy(n2) += Ek_dot_El*VecEy(p2)(j);
                                }
                            }
                        }                      
                    }
                
                // putting non-linear polarization in VecEx, VecEy
                for (int n = 0; n < Nf; n++)
                  {
                    VecEx(n)(j) = PolarEx(n);
                    VecEy(n)(j) = PolarEy(n);
                  }
              }

            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::PROD);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }
            
            // going to Fourier domain
            offset = 0;
            for (int n = 0; n < Nf; n++)
              {
                fft.ApplyInverse(VecEx(n));
                fft.ApplyInverse(VecEy(n));
		
                Mlt(coef_fft, VecEx(n)); Mlt(coef_fft, VecEy(n));
		
                for (int k = 0; k < nb_points_time; k++)
                  {
                    int p = offset+k;
                    scheme(p) = VecEx(n)(k)*conj(phase(p))*invPhaseImag(p);
                    scheme(p+N/2) = VecEy(n)(k)*conj(phase(p))*invPhaseImag(p);
                  }
		
                offset += nb_points_time;
              }
	    
            if (print_level >= 4)
              glob_chrono.Stop(VirtualTimer::FLUX);
          }
	else
	  {
            if (print_level >= 4)
              glob_chrono.Start(VirtualTimer::EXTRAPOL);
            
	    VectComplex_wp utime_x(nb_points_time + padding_time);
	    VectComplex_wp utime_y(nb_points_time + padding_time);
	    utime_x.Fill(zero);
	    utime_y.Fill(zero);
            
            #pragma omp parallel for
	    for (int k = 0; k < N/4; k++)
	      {
		// x-component
		utime_x(k) = u(k)*phase(k)*phaseImag(k);
		utime_x(k + N/4 + padding_time) = u(k+N/4)*phase(k+N/4)*phaseImag(k+N/4);
		// y-component
		utime_y(k) = u(k+N/2)*phase(k)*phaseImag(k);
		utime_y(k + N/4 + padding_time) = u(k+3*N/4)*phase(k+N/4)*phaseImag(k+N/4);
	      }
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::EXTRAPOL);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    // returning back in time
	    fft.ApplyForward(utime_x);    	    
	    fft.ApplyForward(utime_y);    

            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::PROD);
              }
            
            Mlt(1.0/coef_fft, utime_x);
	    Mlt(1.0/coef_fft, utime_y);
            	    
	    // forming 2 |E|^2 E + E.E E*
	    Complex_wp coef1, coef2;
	    if (this->var_raman.GetAlpha() == Real_wp(0))
	      {
		#pragma omp parallel for
		for (int i = 0; i < nb_points_time; i++)
		  {
		    coef1 = absSquare(utime_x(i)) + absSquare(utime_y(i));
		    coef2 = utime_x(i)*utime_x(i) + utime_y(i)*utime_y(i);
		    utime_x(i) = Real_wp(2) * coef1 * utime_x(i) + coef2 * conj(utime_x(i));
		    utime_y(i) = Real_wp(2) * coef1 * utime_y(i) + coef2 * conj(utime_y(i));
		  }
	      }
	    else
	      {
	        VectReal_wp moduleU2(nb_points_time);
		#pragma omp parallel for
	        for (int i = 0; i < nb_points_time; i++)
		  moduleU2(i) = absSquare(utime_x(i)) + absSquare(utime_y(i));
		
		VectComplex_wp moduleU2_chap(nb_points_time/2+1);
		VectReal_wp moduleU2_conv(moduleU2.GetM());
		// R = (1-\alpha) \delta + \alpha h_R
		// treating convolution with h_R
		fft_time.ApplyForward(moduleU2, moduleU2_chap);
		
                #pragma omp parallel for
		for (int i = 0; i < moduleU2_chap.GetM(); i++)
		  moduleU2_chap(i) = moduleU2_chap(i)*conj(raman_fct_hat(i));
		
		fft_time.ApplyInverse(moduleU2_chap, moduleU2_conv);
		
		Real_wp coef_conv = this->var_raman.GetAlpha(), coef_kerr = 1.0-coef_conv;
		
		Complex_wp coef_ux, coef_uxbar;
		#pragma omp parallel for		
		for (int i = 0; i < nb_points_time; i++)
		  {
		    coef2 = utime_x(i)*utime_x(i) + utime_y(i)*utime_y(i);
		    coef_uxbar = coef2*coef_kerr;
		    coef_ux = 2.0*(moduleU2(i)*coef_kerr + moduleU2_conv(i)*coef_conv);
		    utime_x(i) = utime_x(i)*coef_ux + conj(utime_x(i))*coef_uxbar;
		    utime_y(i) = utime_y(i)*coef_ux + conj(utime_y(i))*coef_uxbar;		}
	      }
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::PROD);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    fft.ApplyInverse(utime_x);	    
	    fft.ApplyInverse(utime_y);

            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::EXTRAPOL);
              }

            Mlt(coef_fft, utime_x);    
	    Mlt(coef_fft, utime_y);
            
            #pragma omp parallel for
	    for (int k = 0; k < N/4; k++)
	      {
		// x-component
		scheme(k) = utime_x(k)*conj(phase(k))*invPhaseImag(k);
		scheme(k+N/4) = utime_x(k+N/4+padding_time)*conj(phase(k+N/4))*invPhaseImag(k+N/4);
		// y-component
		scheme(k+N/2) = utime_y(k)*conj(phase(k))*invPhaseImag(k);
		scheme(k+3*N/4) = utime_y(k+N/4+padding_time)*conj(phase(k+N/4))*invPhaseImag(k+N/4);
	      }

            if (print_level >= 4)
              glob_chrono.Stop(VirtualTimer::EXTRAPOL);
	  }
      }
    else
      {	
	Vector<Complex_wp> phase(N);
	Vector<Real_wp> phaseImag(N), invPhaseImag(N);
	if (init)
          {
            omegaAux = omega;
            kwaveAux = kwave;
	    if (fiber_computation)
	      Aeff_aux = Aeff_om;
          }
	
        if (print_level >= 4)
          glob_chrono.Start(VirtualTimer::PML);
	
	phaseImag.Fill(1.0); invPhaseImag.Fill(1.0);
	
        if (NonLinearPhase())
	  {
	    if (local_enveloppe)
	      GetLocalPhase(kwave, zn-z0_current, phase, phaseImag, invPhaseImag);
	    else
	      {
                #pragma omp parallel for
                for (int k = 0; k < N; k++)
                  phase(k) = GetPhaseKz(kwave(k), zn);
              }
	  }
	else
	  {
	    if (local_enveloppe)
	      {
		cout << "Not implemented for this model " << endl;
		abort();
	      }

	    // if the phase is affine, there is no need
	    // to consider it in the non-linear term
	    phase.Fill(one);
	  }
        
        if (print_level >= 4)
          glob_chrono.Stop(VirtualTimer::PML);
	  
	if (type_spectrum == CENTERED)
	  {
            if (print_level >= 4)
              glob_chrono.Start(VirtualTimer::EXTRAPOL);

	    VectComplex_wp ux(omega.GetM());
	    VectReal_wp utime(nb_points_time), moduleU2(nb_points_time);	  
            #pragma omp parallel for
	    for (int k = 0; k < omega.GetM(); k++)
	      ux(k) = conj(u(k)*phase(k))*phaseImag(k);
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::EXTRAPOL);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    fft_time.ApplyInverse(ux, utime);

            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::PROD);
              }
            
	    Mlt(1.0/coef_fft, utime);
            
	    if (type_polarization == CIRCULAR)
	      {
		VectComplex_wp uy(ux);
		VectReal_wp utime_y(nb_points_time);
		Mlt(Iwp, uy);
		
                if (print_level >= 4)
                  {
                    glob_chrono.Stop(VirtualTimer::PROD);	    
                    glob_chrono.Start(VirtualTimer::FLUX);
                  }

		fft_time.ApplyInverse(uy, utime_y);

                if (print_level >= 4)
                  {
                    glob_chrono.Stop(VirtualTimer::FLUX);	    
                    glob_chrono.Start(VirtualTimer::PROD);
                  }
                
		Mlt(1.0/coef_fft, utime_y);
		
		// forming |E|^2
                #pragma omp parallel for
		for (int j = 0; j < nb_points_time; j++)
		  moduleU2(j) = (square(utime(j)) + square(utime_y(j)));
	      }
	    else
	      {
		// forming |E|^2
                #pragma omp parallel for
		for (int j = 0; j < nb_points_time; j++)
		  moduleU2(j) = square(utime(j));
	      }
	    
	    this->ApplyRamanEffect(moduleU2, utime);
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::PROD);
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    fft_time.ApplyForward(utime, ux);

            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::EXTRAPOL);
              }
	    
	    Mlt(coef_fft, ux);
	    
            #pragma omp parallel for
	    for (int k = 0; k < omega.GetM(); k++)
	      scheme(k) = conj(ux(k)*phase(k))*invPhaseImag(k);	    

            if (print_level >= 4)
              glob_chrono.Stop(VirtualTimer::EXTRAPOL);
	    //DISP(scheme.GetNormInf());
	  }
        else if (type_spectrum == ODD_FREQUENCIES)
          {
            // not multithreaded
            if (print_level >= 4)
              glob_chrono.Start(VirtualTimer::FLUX);

            int Nf = number_odd_frequencies;
            Vector<VectComplex_wp> VecEx(Nf);
            VectComplex_wp PolarEx(Nf);
            for (int n = 0; n < Nf; n++)
              VecEx(n).Reallocate(nb_points_time);
            
            // extracting Ex and returning back in time
            int offset = 0;
            for (int n = 0; n < Nf; n++)
              {
                for (int k = 0; k < nb_points_time; k++)
                  VecEx(n)(k) = u(offset + k)*phase(offset+k)*phaseImag(offset+k);
                                  
                fft.ApplyForward(VecEx(n));
                Mlt(1.0/coef_fft, VecEx(n));
                                
                offset += nb_points_time;
              }
                
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::PROD);
              }

            // forming |E|^2 E with odd frequencies
            Real_wp coef_gamma = 1.0;
            Complex_wp Ek_dot_El, Ex_k, Ex_l;
            int Nmax = 2*Nf-1;
            for (int j = 0; j < nb_points_time; j++)
              {
                PolarEx.Fill(zero);
                for (int k2 = -Nf; k2 < Nf; k2++)
                  for (int l2 = k2; l2 < Nf; l2++)
                    {
                      int k = 2*k2+1; int l = 2*l2+1;
                      if (k2 == l2)
                        coef_gamma = 0.25;
                      else
                        coef_gamma = 0.5;
                      
                      // forming coef_gamma Ek . El
                      if (k < 0)
                        Ex_k = conj(VecEx(-k2-1)(j));
                      else
                        Ex_k = VecEx(k2)(j);
                                              
                      if (l < 0)
                        Ex_l = conj(VecEx(-l2-1)(j));
                      else
                        Ex_l = VecEx(l2)(j);
                                              
                      Ek_dot_El = coef_gamma*Ex_k*Ex_l;
                      
                      // then adding terms E_k \cdot E_l E_{n-k-l}
                      for (int n2 = 0; n2 < Nf; n2++)
                        {
                          int n = 2*n2+1;
                          int p = n-k-l;
                          
                          if (p < 0)
                            {
                              if (p >= -Nmax)
                                {
                                  int p2 = -(p+1)/2;
                                  PolarEx(n2) += Ek_dot_El*conj(VecEx(p2)(j));
                                }
                            }
                          else
                            {
                              if (p <= Nmax)
                                {
                                  int p2 = (p-1)/2;
                                  PolarEx(n2) += Ek_dot_El*VecEx(p2)(j);
                                }
                            }
                        }              
                    }
                
                // putting non-linear polarization in VecEx, VecEy
                for (int n = 0; n < Nf; n++)
                  VecEx(n)(j) = PolarEx(n);
              }
            
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::PROD);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

            // going to Fourier domain
            offset = 0;
            for (int n = 0; n < Nf; n++)
              {
                fft.ApplyInverse(VecEx(n));
                Mlt(coef_fft, VecEx(n));
		
                for (int k = 0; k < nb_points_time; k++)
                  {
                    int p = offset+k;
                    scheme(p) = VecEx(n)(k)*conj(phase(p))*invPhaseImag(p);
                  }

                offset += nb_points_time;
              }

            if (print_level >= 4)
              glob_chrono.Stop(VirtualTimer::FLUX);            
          }
	else
	  {
            if (print_level >= 4)
              glob_chrono.Start(VirtualTimer::EXTRAPOL);
	    
            VectComplex_wp utime(nb_points_time + padding_time);
	    utime.Fill(zero);
            
	    #pragma omp parallel for
	    for (int k = 0; k < N/2; k++)
	      {
		utime(k) = u(k)*phase(k)*phaseImag(k);
		utime(k + N/2 + padding_time) = u(k+N/2)*phase(k+N/2)*phaseImag(k+N/2);
	      }
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::EXTRAPOL);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }
            
	    fft.ApplyForward(utime);    
            
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::PROD);
              }
            
	    Mlt(1.0/coef_fft, utime);
            
	    // on constitue |E|^2 E
	    VectReal_wp moduleU2(nb_points_time);
            #pragma omp parallel for
	    for (int i = 0; i < nb_points_time+padding_time; i++)
	      moduleU2(i) = absSquare(utime(i));
	    
	    // the term E.E is neglected in Raman effect because
	    // the delay function attenuates highly variable fields => coefficient 2/3
	    Real_wp coef_corrective(1);
	    if (type_polarization == LINEAR)
	      coef_corrective = Real_wp(2)/3;
	    
	    this->ApplyRamanEffect(moduleU2, utime, coef_corrective);
	    
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::PROD);	    
                glob_chrono.Start(VirtualTimer::FLUX);
              }

	    fft.ApplyInverse(utime);    
            
            if (print_level >= 4)
              {
                glob_chrono.Stop(VirtualTimer::FLUX);	    
                glob_chrono.Start(VirtualTimer::EXTRAPOL);
              }
            
	    Mlt(coef_fft, utime);
	    
            #pragma omp parallel for
	    for (int k = 0; k < N/2; k++)
	      {
		scheme(k) = utime(k)*conj(phase(k))*invPhaseImag(k);
		scheme(k+N/2) = utime(k+N/2+padding_time)*conj(phase(k+N/2))*invPhaseImag(k+N/2);
	      }
            
            if (print_level >= 4)
              glob_chrono.Stop(VirtualTimer::EXTRAPOL);
	    //DISP(scheme.GetNormInf());
	  }
      }
    
    // coefficient of the non linear part depending on the polarization
    Real_wp coefPol(1);
    if (type_spectrum == SHIFTED_SINGLE)
      {
	if (type_polarization == LINEAR)
	  coefPol = Real_wp(3)/4;
	else if(type_polarization == CIRCULAR)
	  coefPol = 1.0 ;
	else
	  coefPol = Real_wp(1)/4;
      }
    
    Complex_wp coefNL = -coefPol * gammaNL / (Real_wp(2)*Iwp*index.c0*index.c0);
    Complex_wp coefAlpha = -index.alpha_damping / Real_wp(2);
    //Complex_wp alpha = coefNL*square(omega_L)/k_L*square(amplitude_impulse);
    //DISP(alpha);

    //Complex_wp rho = coefNL*square(omega_L)/k_L*square(amplitude_impulse);
    //if (zn*z0_adim <= 1e-6)
    //DISP(abs(z0_adim/rho));
    
    if (print_level >= 4)
      glob_chrono.Start(VirtualTimer::STIFFNESS);
    
    if (type_model == SCHRODINGER_CONSTANT)
      {
        if (init)
          {
            coefAlphaEpsilon.Reallocate(N);
	    for (int k = 0; k < N; k++)
	      {
		coefAlphaEpsilon(k) = -imag(kwaveAux(k));
		if (coef_stability_damping > 0)
		  {
		    if (abs(coefAlphaEpsilon(k)) > coef_stability_damping/dz)
		      coefAlphaEpsilon(k) = -coef_stability_damping/dz;
		  }
	      }
	    
	    if (local_enveloppe)
	      coefAlphaEpsilon.Fill(0.0);
          }
        
	coefNL *= square(omega_L)/k_L;
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  scheme(k) = coefAlphaEpsilon(k)*u(k) + scheme(k)*coefNL;
      }
    else if (type_model == UPPE)
      {
        if (init)
          {
            scaleAux.Reallocate(N);
            coefAlphaEpsilon.Reallocate(N);
	    for (int k = 0; k < N; k++)
	      {
		scaleAux(k) = coefNL*square(omegaAux(k))/kwaveAux(k);
		if (fiber_computation)
		  scaleAux(k) /= Aeff_aux(k);
		
		coefAlphaEpsilon(k) = -imag(kwaveAux(k));
		if (coef_stability_damping > 0)
		  {
		    if (abs(coefAlphaEpsilon(k)) > coef_stability_damping/dz)
		      coefAlphaEpsilon(k) = -coef_stability_damping/dz;
		  }
	      }
	    
	    if (local_enveloppe)
	      coefAlphaEpsilon.Fill(0.0);
          }
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  {
	    if (kwaveAux(k) == zero)
	      scheme(k) = 0;
            else
              scheme(k) = scheme(k)*scaleAux(k) + coefAlphaEpsilon(k)*u(k);
	  }
      }
    else if (type_model == FOP)
      {
	coefNL = -gammaNL/(2.0*Iwp)*coefPol;
	Complex_wp coefE = 1.0/Iwp;
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  scheme(k) = (coefAlpha + coefE*(1-omegaAux(k)/omega_L)*k_L)*u(k)
            + coefNL* omegaAux(k)* k_L / omega_L* scheme(k);
	
      }
     else if (type_model == NEE )
      {
	VectComplex_wp D(N);
	
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  D(k) = (kwaveAux(k) - k_L - 1.0/group_velocity * (omegaAux(k)-omega_L))/Iwp;
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  scheme(k)= (coefAlpha - D(k))*u(k) + coefNL*square(omega_L)
            *(1.0+(omegaAux(k)-omega_L)/omega_L)* scheme(k)/ k_L;
      }
    else if (type_model == NEEC )
      {
	VectComplex_wp D(N);
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  D(k) = (kwaveAux(k) - k_L - 1.0/group_velocity * (omegaAux(k)-omega_L))/Iwp;
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  scheme(k) = (coefAlpha - D(k))*u(k) + coefNL*square(omega_L)* scheme(k)/ k_L;
      }
    else if (type_model == NLS)
      {
	Complex_wp sigma;
	sigma = index.GetDispersion(omega_L);

	VectComplex_wp D(N);
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  D(k) = sigma*square(omegaAux(k)-omega_L)/Iwp;
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  scheme(k) = (coefAlpha + D(k))*u(k) + coefNL*square(omega_L)* scheme(k)/ k_L;
      }
    else
      {
	//modele Miro
	Complex_wp sigma, coefE;
	sigma = index.GetDispersion(omega_L);
	coefE = sigma/Iwp;
        
        #pragma omp parallel for
	for (int k = 0; k < N; k++)
	  scheme(k) = (coefAlpha + coefE*square(omegaAux(k)))*u(k)
            + coefNL*square(omega_L)* scheme(k)/ k_L;
      }


    SetDirichletCondition(zn, 1, scheme);

    if (print_level >= 4)
      {
        glob_chrono.Stop(VirtualTimer::STIFFNESS);    
        glob_chrono.Stop(VirtualTimer::JACOBIAN);
      }

    //int test_input; cout << "Waiting ..." << endl; cin >> test_input;
  }
  
  
  //! sets u = 0 for forbidden frequencies
  void KerrProblem::SetDirichletCondition(const Real_wp& z, int nb_deriv,
                                          VectComplex_wp& x, Real_wp alpha)
  {
    Complex_wp zero; SetComplexZero(zero);
    #pragma omp parallel for
    for (int i = 0; i < kwave_to_consider.GetM(); i++)
      if (!kwave_to_consider(i))
        x(i) = zero;
  }
  

  void KerrProblem::SourceOnlyTime(const Real_wp& t0, int nb_deriv, Complex_wp& pulse)
  {
    pulse = 0; 
  }
    
  
  Real_wp KerrProblem::GetFrequency() const
  {
    return omega_L/(2.0*pi_wp);
  }
  
  
  Real_wp KerrProblem::GetFinalTimeSource() const
  {
    return Real_wp(0);
  }

  
  //! solves the num-th operator for a splitting scheme
  void KerrProblem::SolveSplitOperator(const Real_wp& z, const Real_wp& delta_z,
                                       const VectComplex_wp& En, VectComplex_wp& En_next, int num)
  {
    // calcul de En_next en fonction de En avec un pas de temps delta_z

    if(num == 0){

      if (print_level >= 4)
	glob_chrono.Start(VirtualTimer::JACOBIAN);
    
      // on calcule u en temps (avec du zero-padding)
      int N = En.GetM();
      Complex_wp zero;
      SetComplexZero(zero);
          
      //on les utilise pour éviter les répétitions de code dans la partie du choix du modele
      //VectReal_wp omegaAux;
      //VectComplex_wp kwaveAux;
      bool init = false;
      if (omegaAux.GetM() == 0)
	{
	  init = true;
	  omegaAux.Reallocate(N); omegaAux.Fill(0);
	  kwaveAux.Reallocate(N); kwaveAux.Fill(zero);        
	}

      if(type_polarization == ELLIPTIC)
	{
	 if (init)
	   for(int i=0;i<N/2;i++)
	    {
	      omegaAux(i) = omega(i);
	      omegaAux(N/2+i) = omega(i);
	      kwaveAux(i) = kwave(i);
	      kwaveAux(N/2+i) = kwave(i);
	    }
	}
      else
	{
	  if (init)
	    {
	      omegaAux = omega;
	      kwaveAux = kwave;
	    }
	}
	
      if (print_level >= 4)
	glob_chrono.Start(VirtualTimer::PML);
		
      if(type_model == NEEC)
	{
	  
	  Complex_wp D;
	    
	  for(int k=0;k<N;k++)
	    {
	      D = kwaveAux(k)- k_L- 1.0/group_velocity *(omegaAux(k)-omega_L);
	      En_next(k) = En(k)*GetPhaseKz(D,delta_z);
	    }
	}
      else //type_model==NLS
	{
	  Complex_wp sigma,D;
	  sigma = index.GetDispersion(omega_L);
	      
	  for(int k=0;k<N;k++)
	    {
	      D = - sigma * square(omegaAux(k)-omega_L);
	      En_next(k) = En(k)*GetPhaseKz(D,delta_z);
	    }
	}
    }
    else //if num==1
      { 
	int N = En.GetM();
	Complex_wp zero;
	SetComplexZero(zero);

	if(type_polarization == ELLIPTIC)
	  {
	    if(type_spectrum == SHIFTED_SINGLE)
	      {
		if (print_level >= 4)
		  glob_chrono.Start(VirtualTimer::FLUX);

		VectComplex_wp utime_x(nb_points_time + padding_time);
		VectComplex_wp utime_y(nb_points_time + padding_time);
		utime_x.Fill(zero);
		utime_y.Fill(zero);

		for (int k = 0; k < N/4; k++)
		  {
		    // x-component
		    utime_x(k) = En(k);
		    utime_x(k + N/4 + padding_time) = En(k+N/4);
		    // y-component
		    utime_y(k) = En(k+N/2);
		    utime_y(k + N/4 + padding_time) = En(k+3*N/4);
		  }
	    
		// returning back in time
		fft.ApplyForward(utime_x);    
		Mlt(1.0/coef_fft, utime_x);
		fft.ApplyForward(utime_y);    
		Mlt(1.0/coef_fft, utime_y);

		if (print_level >= 4)
		  {
		    glob_chrono.Stop(VirtualTimer::FLUX);	    
		    glob_chrono.Start(VirtualTimer::PROD);
		  }
            	    
		// forming the temporal solution
		Complex_wp coefNL,coefE1,coefE2,coef1,coef2,coef3,coef4,l1,l2,utime_x0,utime_y0;
		coefNL = gammaNL/(Real_wp(2)*index.c0*index.c0);
		
		for(int i = 0; i< nb_points_time+padding_time; i++)
		  {
		    coef1 = absSquare(utime_x(i)) + absSquare(utime_y(i));
		    coef2 = imag(utime_x(i)*conj(utime_y(i)));
		    l1 = 3.0 * coef1 + 2.0 * coef2;
		    l2 = 3.0 * coef1 - 2.0 * coef2;
		    coefE1 = 1.0/4.0 * coefNL * square(omega_L) / k_L * l1;
		    coefE2 = 1.0/4.0 * coefNL * square(omega_L) / k_L * l2;
		    coef3 = 0.5 * (GetPhaseKz(coefE1,delta_z) + GetPhaseKz(coefE2, delta_z));
		    coef4 = 0.5 * (GetPhaseKz(coefE1,delta_z) - GetPhaseKz(coefE2, delta_z));
		    utime_x0 = utime_x(i);
		    utime_y0 = utime_y(i);
		    utime_x(i) = coef3 * utime_x0 - Iwp * coef4 * utime_y0;
		    utime_y(i) = Iwp * coef4 * utime_x0 + coef3 * utime_y0;
		  }
	    
		if (print_level >= 4)
		  {
		    glob_chrono.Stop(VirtualTimer::PROD);	    
		    glob_chrono.Start(VirtualTimer::FLUX);
		  }

		fft.ApplyInverse(utime_x);
		Mlt(coef_fft, utime_x);
		fft.ApplyInverse(utime_y);    
		Mlt(coef_fft, utime_y);
         
		for (int k = 0; k < N/4; k++)
		  {
		    // x-component
		    En_next(k) = utime_x(k);
		    En_next(k+N/4) = utime_x(k+N/4+padding_time);
		    // y-component
		    En_next(k+N/2) = utime_y(k);
		    En_next(k+3*N/4) = utime_y(k+N/4+padding_time);
		  }

		if (print_level >= 4)
		  glob_chrono.Stop(VirtualTimer::FLUX);

		
	      }
	  }
	else //circular or linear polarisation
	  {
	    if(type_spectrum == SHIFTED_SINGLE)
	      {
		if (print_level >= 4)
		  glob_chrono.Start(VirtualTimer::FLUX);
	    
		VectComplex_wp utime(nb_points_time + padding_time);
		utime.Fill(zero);
	    
		for (int k = 0; k < N/2; k++)
		  {
		    utime(k) = En(k);
		    utime(k + N/2 + padding_time) = En(k+N/2);
		  }
	    
		fft.ApplyForward(utime);    
		Mlt(1.0/coef_fft, utime);

		if (print_level >= 4)
		  {
		    glob_chrono.Stop(VirtualTimer::FLUX);	    
		    glob_chrono.Start(VirtualTimer::PROD);
		  }
            
		// on constitue la solution en temporel
		//DISP(utime.GetNormInf());
		Complex_wp coefNL , coefE, coefPol;
		if(type_polarization == LINEAR)
		  coefPol= Real_wp(3)/4;
		else //type_polarization == CIRCULAR
		  coefPol=1.0;
		
		coefNL = gammaNL/(Real_wp(2)*index.c0*index.c0);
		
		for (int i = 0; i < nb_points_time+padding_time; i++)
		  {
		    coefE = coefNL*coefPol*square(omega_L)*absSquare(utime(i))/k_L ;
		    utime(i) = utime(i)* GetPhaseKz(coefE,delta_z);
		  }
		
		if (print_level >= 4)
		  {
		    glob_chrono.Stop(VirtualTimer::PROD);	    
		    glob_chrono.Start(VirtualTimer::FLUX);
		  }

		fft.ApplyInverse(utime);    
		Mlt(coef_fft, utime);
	    
		//solution en fréquentiel
		for (int k = 0; k < N/2; k++)
		  {
		    En_next(k) = utime(k);
		    En_next(k+N/2) = utime(k+N/2+padding_time);
		  }

		if (print_level >= 4)
		  glob_chrono.Stop(VirtualTimer::FLUX);
		//DISP(scheme.GetNormInf());
	      }
	    //else type_spectrum
	  }
      }
  }
  
  
  //! returns the harmonic pulse
  Complex_wp KerrProblem::GetHarmonicPulse(const Real_wp& w)
  {
    // basic gaussian
    // E(z=0, omega) = amplitude tau exp( -1/2 tau^2 (omega - omega_L)^2)*polarization
    return exp(-0.5*square(tau_impulse*(w-omega_L)));
  }
  
  
  //! Sets E0 to the initial condition E(z=0, omega)
  void KerrProblem::GetInitialCondition(Vector<Complex_wp>& E0)
  {
    if (type_polarization == ELLIPTIC)
      {
	E0.Reallocate(2*omega.GetM());
	Complex_wp coef = amplitude_impulse*tau_impulse, pulse;
	for (int i = 0; i < omega.GetM(); i++)
	  {
	    pulse = GetHarmonicPulse(omega(i));
	    E0(i) = coef*pulse*polarization(0);
	    E0(omega.GetM()+i) = coef*pulse*polarization(1);
	  }
      }
    else
      {
	E0.Reallocate(omega.GetM());

	Complex_wp coef = amplitude_impulse*tau_impulse, pulse;
	for (int i = 0; i < omega.GetM(); i++)
	  {
	    pulse = GetHarmonicPulse(omega(i));
	    E0(i) = coef*pulse*polarization(0);
	  }
      }
  }
  
  
  //! returns energy associated with En
  Real_wp KerrProblem::GetEnergy(VectComplex_wp& En) const
  {
    Real_wp norme = 0;
    for (int i = 0; i < En.GetM(); i++)
      norme += absSquare(En(i))*n2omega(i);
    
    return norme;
  }
  
  
  //! iterate E^n is given to this object  
  void KerrProblem::GiveIterate(int n, const Real_wp& t, VectComplex_wp& En)
  {
    if (n == 0)
      return;
    
    if (print_level >= 4)
      glob_chrono.Start(VirtualTimer::MASS);
    
    if (n%nb_iterations_between_two_tests != 0)
      return;
    
    // testing overflow    
    if (check_overflow_scheme)
      for (int i = 0; i < En.GetM(); i++)
        if ((abs(En(i)) > 1e100) || std::isnan(abs(En(i))))
          {
            DISP(i); DISP(En(i));
            cout << "Instable scheme ?" << endl;
            abort();
          }
    
    // testing that the solution does not touch too much
    // the domain where epsilon is negative
    Real_wp max_En = 0, norm_En = En.GetNormInf();
    for (int i = 0; i < En.GetM(); i++)
      if (real(square(kwaveAux(i))) < Real_wp(0))
        max_En = max(max_En, abs(En(i)));
    
    if (print_level >= 2)
      if (n%100 == 0)
        cout << "Part of the solution in negative epsilon = " << max_En/norm_En << endl;
    
    if (max_En/norm_En > threshold_negative_epsilon)
      {
        //DISP(max_En/amplitude_impulse);
        //cout << "The solution expands to a forbidden spectrum" << endl;
        //abort();
      }

    // testing the last frequency
    int N = kwave.GetM();
    int offset_last = 0;
    if (type_spectrum == SHIFTED_SINGLE)
      offset_last = -N/2;
    
    max_En = 0;
    for (int k = 1; k <= En.GetM()/N; k++)
      max_En = max(max_En, abs(En(k*N+offset_last-1)));
    
    if (print_level >= 2)
      if (n%100 == 0)
        cout << "amplitude for the last frequency = " << max_En/norm_En << endl;
    
    if (max_En/norm_En > threshold_last_frequency)
      {
        // DISP(max_En/amplitude_impulse);
        //cout << "The solution expands to a forbidden spectrum" << endl;
        //abort();
      }    

    if (print_level >= 4)
      glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! modifies iterates
  void KerrProblem::ModifyIterate(int n, const Real_wp& z, VectComplex_wp& En)
  {
    if (print_level >= 4)
      glob_chrono.Start(VirtualTimer::PML);

    if (NonLinearPhase() && local_enveloppe)
      {
	z0_current = z;
	if (type_polarization == ELLIPTIC)
	  {
	    abort();
	  }
	else
	  {
	    int N = En.GetM();
	    Vector<Complex_wp> phase(N);
	    VectReal_wp phase_real(N), invPhase_real(N);
	    
	    GetLocalPhase(kwave, dz, phase, phase_real, invPhase_real);
	    
	    for (int k = 0; k < En.GetM(); k++)
	      En(k) *= phase(k)*phase_real(k);
	  }
      }

    if (print_level >= 4)
      glob_chrono.Stop(VirtualTimer::PML);
  }
  
  
  //! extracts values that are above a given threshold
  template<class T> void KerrProblem::
  ExtractNonZeroComponents(Real_wp& t0, Real_wp& t1, Vector<T>& utime, const Real_wp& threshold)
  {
    int Nt = utime.GetM();
    int i0(0), i1(Nt-1);
    Real_wp t_begin(t0), t_end(t1);
    for (int i = 0; i < Nt; i++)
      if (abs(utime(i)) > threshold)
        {
          i0 = i;
          break;
        }
    
    for (int i = Nt-1; i >= 0; i--)
      if (abs(utime(i)) > threshold)
        {
          i1 = i;
          break;
        }
    
    // if all the values are below the given threshold, the function is left
    // without changing any arrays
    if (i0 > i1)
      return;
    
    // same thing if no value is dropped
    if ((i0 == 0) && (i1 == Nt-1))
      return;
    
    int N = 0;
    if (size_refinement_output == 0)
      {
	// we search the closest power of two
	N = FindClosestPow2(i1-i0+1);
      }
    else
      {
	// otherwise a multiple of size_refinement_output is searched
	int N = size_refinement_output;
	if (N < (i1-i0+1))
	  {
	    while (N < (i1-i0+1))
	      N *= 2;
	  }		
      }

    // if it is too large, no extraction
    if (N >= Nt)	  
      return;
    
    // otherwise we extend i0, i1
    ExtendIntegerInterval(i0, i1, Nt, N);
    
    Vector<T> u_sub(i1-i0+1);
    for (int i = i0; i <= i1; i++)
      u_sub(i-i0) = utime(i);
    
    Real_wp dt = (t1-t0)/(Nt-1);
    t_begin = t0 + i0*dt; t_end = t0 + i1*dt;
    
    t0 = t_begin; t1 = t_end;
    SwapPointer(u_sub, utime);
  }
  
  
  //! writes En in real space and Fourier space
  void KerrProblem::
  WriteOutputFile(const string& file_name, const string& file_name_time,
		  const Real_wp& zn, Vector<Complex_wp>& En)
  {
    if (print_level >= 4)
      glob_chrono.Start(VirtualTimer::OUTPUT);
    
    if (type_polarization == ELLIPTIC)
      {
	int N = omega.GetM();
	VectComplex_wp Ex(N), Ey(N);
	for (int i = 0; i < N; i++)
	  {
	    Ex(i) = En(i);
	    Ey(i) = En(i+N);
	  }
	
	string name_file = GetBaseString(file_name);
	string name_x = name_file + "_X.dat";
	string name_time_x = name_file + "_timeX.dat";
	WriteOutputFileScalar(name_x, name_time_x, zn, Ex);
	
	string name_y = name_file + "_Y.dat";
	string name_time_y = name_file + "_timeY.dat";
	WriteOutputFileScalar(name_y, name_time_y, zn, Ey);
      }
    else
      WriteOutputFileScalar(file_name, file_name_time, zn, En);

    if (print_level >= 4)
      glob_chrono.Stop(VirtualTimer::OUTPUT);
  }
  
  
  //! returns true if k(omega) is a non-affine function of omega
  bool KerrProblem::NonLinearPhase() const
  {
    if ((type_model == UPPE) || (type_model == SCHRODINGER_CONSTANT))
      return true;
    
    return false;
  }
  
  
  //! writes En in real space and Fourier space
  void KerrProblem::
  WriteOutputFileScalar(const string& file_name, const string& file_name_time,
			const Real_wp& zn, Vector<Complex_wp>& En)
  {    
    Complex_wp zero; SetComplexZero(zero);
    Vector<Complex_wp> u(En);
    Mlt(E0_adim, u);
    int N = nb_points_time;
    Complex_wp phase;	
    Complex_wp vg = group_velocity;
    Real_wp omega0(0), omega1(0);
    cout << "Decalage en temps : " << zn/real(vg)*t0_adim << endl;

    GridInterpolationFull<Dimension1> grid;
    int type = GridInterpolationFull<Dimension1>::LINE;
    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
    bool ascii = false;
    Real_wp threshold = threshold_time_extraction*amplitude_impulse*E0_adim;
	
    if (type_spectrum == CENTERED)
      {
	omega0 = omega(0);
	omega1 = omega(omega.GetM()-1);
	
        int Nt = nb_points_time;
        if (factor_refinement_time > 1)
          {
            Nt *= factor_refinement_time;
            u.Resize(Nt/2+1);
            // appending 0 => refinement in time
            for (int i = En.GetM(); i < u.GetM(); i++)
              u(i) = zero;
          }
        
	// returning back in time
        Complex_wp phase_cte = GetPhaseKz(k_L, zn);
        if (NonLinearPhase())
          {
            for (int j = 0; j < En.GetM(); j++)
              {                
		if (local_enveloppe)
		  phase = GetPhaseKz(-omega(j)/vg, zn);
		else
		  phase = GetPhaseKz(kwave(j)-omega(j)/vg, zn);
		
                u(j) = conj(u(j)*phase);
              }
          }
        else
          {
            phase = phase_cte;
            for (int j = 0; j < En.GetM(); j++)
              u(j) = conj(u(j)*phase);
	  }
	
        Real_wp t0 = -Tmax + zn/real(vg), t1 = zn/real(vg) + time_interval(Nt/2-1);
	if (affiche_enveloppe)
          {
            VectComplex_wp utime(Nt);
            utime.Fill(zero);
            utime(0) = u(Nt/2);
            for (int i = 0; i < Nt/2; i++)
              utime(i) = u(i);
            
            fft_envelope_time.ApplyInverse(utime);
            Mlt(1.0/coef_fft_envelope_time, utime);
            
            // iffshift on utime 
            Complex_wp vloc;
            for (int i = 0; i < Nt/2; i++)
              {
                vloc = utime(i);
                utime(i) = utime(Nt/2+i);
                utime(Nt/2+i) = vloc;
              }
            
            if (threshold_time_extraction > Real_wp(0))
              ExtractNonZeroComponents(t0, t1, utime, threshold);
            
            // using loadND format instead
            grid.Init(type, t0*t0_adim, t1*t0_adim, utime.GetM());
            WriteMatlab(utime, grid, DOSSIER + file_name_time, double_prec, ascii);
          }
        else
          {
            VectReal_wp utime(Nt);
            
            if (factor_refinement_time > 1)
              {
                fft_refinement_time.ApplyInverse(u, utime);	
                Mlt(1.0/coef_fft_refinement, utime);
              }
            else
              {
                fft_time.ApplyInverse(u, utime);
                Mlt(1.0/coef_fft, utime);
              }
            
            // iffshift on utime 
            Real_wp vloc;
            for (int i = 0; i < Nt/2; i++)
              {
                vloc = utime(i);
                utime(i) = utime(Nt/2+i);
                utime(Nt/2+i) = vloc;
              }
            
            if (threshold_time_extraction > Real_wp(0))
              ExtractNonZeroComponents(t0, t1, utime, threshold);
            
            // using loadND format instead
            grid.Init(type, t0*t0_adim, t1*t0_adim, utime.GetM());
            WriteMatlab(utime, grid, DOSSIER + file_name_time, double_prec, ascii);
          }
        
        // writing u in Fourier space
        if (factor_refinement_time > 1)
          u.Resize(En.GetM());
        
        // if the envelope is required, we multiply by the constant phase exp(i k_L z)
	if (affiche_enveloppe)
          for (int j = 0; j < En.GetM(); j++)
            u(j) = conj(u(j)*phase_cte);
        else
          for (int j = 0; j < En.GetM(); j++)
            u(j) = conj(u(j));

	
	Mlt(t0_adim, u);    

	grid.Init(type, omega0/t0_adim, omega1/t0_adim, u.GetM());
	WriteMatlab(u, grid, DOSSIER + file_name, double_prec, ascii);	
      }
    else if (type_spectrum == ODD_FREQUENCIES)
      {
        Complex_wp phase_cte = GetPhaseKz(k_L, zn);
        if (NonLinearPhase())
          {
            for (int k = 0; k < omega.GetM(); k++)
              {
		if (local_enveloppe)
		  phase = GetPhaseKz(-omega(k)/vg, zn);
		else
		  phase = GetPhaseKz(kwave(k)-omega(k)/vg, zn);
		
                if (!affiche_enveloppe)
                  u(k) *= phase;
                else
                  u(k) *= phase*conj(phase_cte);
              }            
          }
        else
          {
            if (!affiche_enveloppe)
              Mlt(phase_cte, u);
          }

	// writing u in Fourier space
	Mlt(t0_adim, u);    
	
        VectComplex_wp uchap(nb_points_time);
        int offset = 0;
	for (int n = 0; n < number_odd_frequencies; n++)
          {
            string nom = GetBaseString(file_name) + "_Freq" + to_str(n) + ".dat";
            for (int k = 0; k < nb_points_time/2; k++)
              {
                uchap(k) = u(nb_points_time/2+offset+k);
                uchap(nb_points_time/2+k) = u(offset+k);
              }
            
            omega0 = omega(offset+nb_points_time/2); omega1 = omega(offset+nb_points_time/2-1);
            grid.Init(type, omega0/t0_adim, omega1/t0_adim, uchap.GetM());
            WriteMatlab(uchap, grid, DOSSIER + nom, double_prec, ascii);
            offset += nb_points_time;
          }
        
        // output in time not implemented
      }
    else
      {	
        int Nt = nb_points_time;
        if (factor_refinement_time > 1)
          Nt *= factor_refinement_time;
        
	// then computing u in real space
	VectComplex_wp utime(Nt);
	utime.Fill(zero);
	Complex_wp phase_cte = GetPhaseKz(k_L, zn);	
	if (NonLinearPhase())
	  {
	    for (int k = 0; k < N/2; k++)
	      {
		if (local_enveloppe)
		  phase = GetPhaseKz(-omega(k)/vg, zn);
		else
		  phase = GetPhaseKz(kwave(k)-omega(k)/vg, zn);
		
		utime(k) = u(k)*phase;
                u(k) = utime(k);

		if (local_enveloppe)
		  phase = GetPhaseKz(-omega(k+N/2)/vg, zn);
		else
		  phase = GetPhaseKz(kwave(k+N/2)-omega(k+N/2)/vg, zn);
		
		utime(k+Nt-N/2) = u(k+N/2)*phase;
                u(k+N/2) = utime(k+Nt-N/2);
	      }
	  }
	else
	  {
            phase = phase_cte;
	    for (int k = 0; k < N/2; k++)
	      {
		u(k) *= phase;
                utime(k) = u(k);

		u(k+N/2) *= phase;
                utime(k+Nt-N/2) = u(k+N/2);
	      }
	  }
	
	if (factor_refinement_time > 1)
          {
            fft_refinement.ApplyForward(utime);    
            Mlt(1.0/coef_fft_refinement, utime);
          }
        else
          {
            fft.ApplyForward(utime);    
            Mlt(1.0/coef_fft, utime);
          }
	
	// computing utime = ifftshift(utime)
	Complex_wp vloc;
	for (int k = 0; k < Nt/2; k++)
	  {
	    vloc = utime(k);
	    utime(k) = utime(Nt/2+k);
	    utime(Nt/2+k) = vloc;
	  }
	
	Real_wp t0 = -Tmax + zn/real(vg), t1 = zn/real(vg) + time_interval(Nt/2-1);
        if (threshold_time_extraction > Real_wp(0))
          ExtractNonZeroComponents(t0, t1, utime, threshold);

	// using loadND format instead
	grid.Init(type, t0*t0_adim, t1*t0_adim, utime.GetM());
	WriteMatlab(utime, grid, DOSSIER + file_name_time, double_prec, ascii);
	
	// writing u in Fourier space
	Mlt(t0_adim, u);    
	
        // if the enveloppe is required we multiply by exp(-i k_L z)
        Mlt(conj(phase_cte), u);
        
	// computing u = ifftshift(u)
	for (int k = 0; k < N/2; k++)
	  {
	    vloc = u(k);
	    u(k) = u(N/2+k);
	    u(N/2+k) = vloc;
	  }
	
	omega0 = omega(nb_points_time/2); omega1 = omega(nb_points_time/2-1);
	grid.Init(type, omega0/t0_adim, omega1/t0_adim, u.GetM());
	WriteMatlab(u, grid, DOSSIER + file_name, double_prec, ascii);
      }
  }
  
}

#define MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_1D_CXX
#endif

