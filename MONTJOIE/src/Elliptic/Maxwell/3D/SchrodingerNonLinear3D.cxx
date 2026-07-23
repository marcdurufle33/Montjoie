#ifndef MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_3D_CXX

namespace Montjoie
{
  
  //! default constructor
  SchrodingerProblem3D::SchrodingerProblem3D()
  {
    nb_points_time = 1;
    padding_time = 0;
    type_forcing_number_fft_points = POWER_NONE;    
    
    gammaNL = 0.0;
    
    omega_L = 0.0; k_L = 0;
    
    tau_impulse = 1.0;
    amplitude_impulse = 1.0;
    Tmax = 1.0;
    
    dz = 0.0;
    zmax = 1.0;
    
    dx = 0.0; dy = 0.0;
    nb_points_x = 0;
    nb_points_y = 0;
        
    xmin = 0.0; xmax = 0.0;
    ymin = 0.0; ymax = 0.0;    
    waist = 0.0;
    
    z0_adim = 1.0; E0_adim = 1.0; t0_adim = 1.0;
    
    // default directory is the current directory
    DOSSIER = "./";
    
    coef_fft = 1.0;
    type_model = UPPE;
    
    type_polarization = CIRCULAR;
    polarization.Init(Complex_wp(1, 0), Complex_wp(0, 0));
    
    local_enveloppe = false;
    z0_current = 0.0; threshold_phase_imag = 1e-100;
    
    type_spectrum = SHIFTED_SINGLE;
    number_odd_frequencies = 1;
    
    affiche_enveloppe = false;
    
    factor_refinement_time = 1;
    threshold_time_extraction = 0;
    size_refinement_output = 1;
    
    coef_fft_envelope_time = 1.0;
    check_overflow_scheme = true;
    threshold_last_frequency = 1e-4;
    threshold_negative_epsilon = 1e-6;        
    
    epsilon_min_treated = -1e100; epsilon_max_treated = 1e100;
    coef_stability_damping = 0.0;

    type_approximation_Ez = EZ_NULL;
    
    print_level = 0;
  }
  
  
  //! returns the space step along z
  Real_wp SchrodingerProblem3D::GetDz() const
  {
    return dz;
  }
  
  
  //! returns the final value of z 
  Real_wp SchrodingerProblem3D::GetZmax() const
  {
    return zmax;
  }
  
  
  //! returns the number of points in the 2-D interpolation grid
  int SchrodingerProblem3D::GetNbPointsGrid() const
  {
    return grid_x.GetNbPointsGrid()*grid_y.GetNbPointsGrid();
  }
  
  
  //! returns the number of points in time
  int SchrodingerProblem3D::GetNbPointsTime() const
  {
    return nb_points_time;
  }
  
  
  //! modifying parameters with a line of the datafile
  void SchrodingerProblem3D::SetInputData(const string& keyword, const Vector<string>& parameters)
  {
    if (keyword == "WavelengthAdim")
      {
	// WavelengthAdim = z0
        if (parameters.GetM() <= 0)
          {
	    cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
      }
    else if (keyword == "RefractiveIndex")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
    else if (keyword == "ComputationalGrid")
      {
        // ComputationalGrid = xmin xmax ymin ymax nb_points_x nb_points_y
        if (parameters.GetM() <= 5)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
            cout << "ComputationalGrid needs six parameters, for instance :" << endl;
            cout << "ComputationalGrid = x0 xN y0 yN Nx Ny" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        Real_wp x0 = to_num<Real_wp>(parameters(0))/z0_adim;
        Real_wp xN = to_num<Real_wp>(parameters(1))/z0_adim;
        Real_wp y0 = to_num<Real_wp>(parameters(2))/z0_adim;
        Real_wp yN = to_num<Real_wp>(parameters(3))/z0_adim;
        int Nx = to_num<int>(parameters(4));
        int Ny = to_num<int>(parameters(5));
        SetRegularGrid(x0, xN, y0, yN, Nx, Ny);
      }
    else if (keyword == "DisplayGrid")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "Provide more parameters" << endl;
            abort();
          }
        
        if (parameters(0) == "AUTO")
          {            
            InitGrid(xmin, xmax, ymin, ymax, nb_points_x+1, nb_points_y+1);
          }
	else if (parameters(0) == "AUTO_X")
	  {
	    InitGrid(xmin, xmax, (ymin+ymax)/2, (ymin+ymax)/2, nb_points_x+1, 1);
	  }
        else
          {
            if (parameters.GetM() <= 5)
              {
                cout << "In SetInputData of SchrodingerProblem3D" << endl;
                cout << "DisplayGrid needs six parameters, for instance :" << endl;
                cout << "DisplayGrid = x0 xN y0 yN Nx Ny" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
            Real_wp x0 = to_num<Real_wp>(parameters(0))/z0_adim;
            Real_wp xN = to_num<Real_wp>(parameters(1))/z0_adim;
            Real_wp y0 = to_num<Real_wp>(parameters(2))/z0_adim;
            Real_wp yN = to_num<Real_wp>(parameters(3))/z0_adim;
            int Nx = to_num<int>(parameters(4));
            int Ny = to_num<int>(parameters(5));
            InitGrid(x0, xN, y0, yN, Nx, Ny);
          }
      }
    else if (keyword == "DisplayZ")
      {
	// DisplayZ = z0 z1 dz_display
        output_grid_z.SetInputData(keyword, parameters);
	output_grid_z.AdimTime(z0_adim);
      }
    else if (keyword == "FileOutputDisplayZ")
      {
	// FileOutputDisplayZ = name_file
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
	    cout << "FileOutputDisplayZ needs one parameter, for instance :" << endl;
	    cout << "FileOutputDisplayZ = file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
          }

        output_grid_z.SetTotalFieldFile(parameters(0));
      }
    else if (keyword == "PropagationStep")
      {
	// PropagationStep = dz
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
	    cout << "Wavelength needs one parameter, for instance :" << endl;
	    cout << "Wavelength = file_name" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

	if (parameters(0) == "UPPE")
	  type_model = UPPE;
	else if (parameters(0) == "UPPE_PARAXIAL")
	  type_model = UPPE_PARAXIAL;
      }
    else if (keyword == "TimeImpulsion")
      {
	// TimeImpulsion = N Npadding Gaussian tau Tmax amplitude
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
                cout << "In SetInputData of SchrodingerProblem3D" << endl;
                cout << "TimeImpulsion needs six parameters, for instance :" << endl;
                cout << "TimeImpulsion = N Npadding Gaussian tau Tmax amplitude" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }
            
	    if (parameters(2) == "GaussianHalf")
              {
                Real_wp fwhm_impulse = to_num<Real_wp>(parameters(3));
                //Relation between fwhm and tau
                tau_impulse = fwhm_impulse /(Real_wp(2) * sqrt(Real_wp(2)*log(Real_wp(2)))); 
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
            if (parameters.GetM() <= 3)
              {
                cout << "In SetInputData of SchrodingerProblem3D" << endl;
                cout << "TimeImpulsion needs four parameters, for instance :" << endl;
                cout << "TimeImpulsion = N Npadding Monochromatic amplitude" << endl;
                cout << "Current parameters are : " << endl << parameters << endl;
                abort();
              }

            nb_points_time = 1;
            amplitude_impulse = to_num<Real_wp>(parameters(3));

            if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
              amplitude_impulse /= E0_adim;            
          }
        else
          {
            abort();
          }
      }
    else if (keyword == "Waist")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
            cout << "Waist needs one parameter, for instance :" << endl;
            cout << "Waist = w" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        waist = to_num<Real_wp>(parameters(0))/z0_adim;
      }
    else if (keyword == "Polarization")
      {	
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
                cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
                cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
            cout << "In SetInputData of SchrodingerProblem3D" << endl;
            cout << "PrintLevel needs one parameter, for instance :" << endl;
            cout << "PrintLevel = 2" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        print_level = to_num<int>(parameters(0));
      }
    else if (keyword == "DirectoryOutput")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of SchrodingerProblem3D" << endl;
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
  void SchrodingerProblem3D::SetFrequency(const Real_wp& f)
  {
    omega_L = 2.0*pi_wp*f;
  }
  
  //! setting the computational domain
  void SchrodingerProblem3D
  ::SetRegularGrid(const Real_wp& x0, const Real_wp& xN, const Real_wp& y0,
                   const Real_wp& yN, int nbx, int nby)
  {
    xmin = x0;
    xmax = xN;
    ymin = y0;
    ymax = yN;
    nb_points_x = nbx;
    nb_points_y = nby;

    // condition periodique en x et y, on ne discretise
    // pas le dernier point 
    dx = (xN - x0)/nb_points_x;
    dy = (yN - y0)/nb_points_x;
    step_x.Reallocate(nb_points_x);
    step_y.Reallocate(nb_points_y);
    for (int i = 0; i < nb_points_x; i++)
      step_x(i) = x0 + i*dx;

    for (int i = 0; i < nb_points_y; i++)
      step_y(i) = y0 + i*dy;
  }
  
  
  //! constructs needed arrays to complete the simulation
  void SchrodingerProblem3D::ConstructAll()
  {
    Complex_wp zero; SetComplexZero(zero);
    k_L = index.GetWaveNumber(omega_L);
    
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
	    kwave(i) = index.GetWaveNumber(omega(i));
	    
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
	    kwave(i+nb_points_time/2) = index.GetWaveNumber(omega_);
	  }
	
	// writing pulsations and wave numbers on the disk
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
    
	// initializing FFT
	fft.Init(nb_points_time + padding_time);
	coef_fft = 2.0*Tmax/sqrt(2.0*pi_wp);

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
        for (int i = 1; i <= nb_points_time/2; i++)
          {
            omega(i) = Omega(nb_points_time/2+i);
            kwave(i) = index.GetWaveNumber(omega(i));
	  }
	
	omega(0) = 0; kwave(0) = 0;
	
        // writing pulsations and wave numbers on the disk
        VectReal_wp omega_display(omega);
        VectComplex_wp kwave_display(kwave);
        
        Mlt(1.0/t0_adim, omega_display);
        Mlt(1.0/z0_adim, kwave_display);
        
        omega_display.Write(DOSSIER + "omega.dat");
        kwave_display.Write(DOSSIER + "kwave.dat");    

        // initializing FFT
        fft_time.Init(nb_points_time + padding_time);
        coef_fft = 4.0*Tmax/(nb_points_time+padding_time)/sqrt(2.0*pi_wp);

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
    
    fft2d.Init(nb_points_x, nb_points_y);
    
    //EnHat.Reallocate(nb_points_x*nb_points_y*nb_points_time);
    //EnHat.Fill(0);
    
    // "frequences" ou on evalue la fft de E
    nu_x.Reallocate(nb_points_x);
    nu_y.Reallocate(nb_points_y);
    for (int i = 0; i < nb_points_x/2; i++)
      {
        nu_x(i) = 2.0*pi_wp*Real_wp(i)/(xmax-xmin);
        nu_x(nb_points_x-1-i) = -2.0*pi_wp*Real_wp(i+1)/(xmax-xmin);
      }

    for (int i = 0; i < nb_points_y/2; i++)
      {
        nu_y(i) = 2.0*pi_wp*Real_wp(i)/(ymax-ymin);
        nu_y(nb_points_y-1-i) = -2.0*pi_wp*Real_wp(i+1)/(ymax-ymin);
      }
    
    // coef_Gamma2 = kx^2 + ky^2
    // kz = \sqrt(k^2 - kx^2 -ky^2)
    int N = nb_points_x*nb_points_y;
    coef_Gamma2.Reallocate(nb_points_x*nb_points_y);
    coef_Kz.Reallocate(nb_points_x*nb_points_y*nb_points_time);
    coef_InvKz.Reallocate(nb_points_x*nb_points_y*nb_points_time);
    Real_wp kx2_ky2;
    for (int i = 0; i < nb_points_x; i++)
      for (int j = 0; j < nb_points_y; j++)
        {
	  kx2_ky2 = nu_x(i)*nu_x(i) + nu_y(j)*nu_y(j);
	  coef_Gamma2(i*nb_points_y + j) = kx2_ky2;
	  for (int k = 0; k < nb_points_time; k++)
            {
              if (kx2_ky2 < real(kwave(k)*kwave(k)))
                coef_Kz(N*k + i*nb_points_y + j) = sqrt(kwave(k)*kwave(k) - kx2_ky2);            
              else
                coef_Kz(N*k + i*nb_points_y + j) = zero;
            }
	}
    
    for (int i = 0; i < coef_Kz.GetM(); i++)
      {
        if (coef_Kz(i) != zero)
          coef_InvKz(i) = 1.0/coef_Kz(i);
        else
          coef_InvKz(i) = zero;
      }
    
    // on regarde si epsilon est dans le bon intervalle
    kwave_to_consider.Reallocate(kwave.GetM());
    kwave_to_consider.Fill(true);
    for (int i = 0; i < kwave.GetM(); i++)
      if (omega(i) != Real_wp(0))
        {
          Complex_wp eps = square(kwave(i)*index.c0/omega(i));
          if ((real(eps) < epsilon_min_treated) || (real(eps) > epsilon_max_treated))
            kwave_to_consider(i) = false;
          else
            kwave_to_consider(i) = true;
        }
    
    // pre-computing coefficient for splitting schemes
    if (time_scheme.GetTimeSchemeType() == TimeSchemeEnum::SPLIT_EXACT)
      {
        if (type_model == UPPE_PARAXIAL)
          {
            coef_Laplacian.Reallocate(coef_Gamma2.GetM());
            Complex_wp beta = dz/(2.0*Iwp*k_L);
            for (int i = 0; i < coef_Gamma2.GetM(); i++)
              coef_Laplacian(i) = exp(0.5*beta*coef_Gamma2(i));
          }
        else if (type_model == UPPE)
          {
            for (int i = 0; i < coef_Gamma2.GetM(); i++)
              coef_Laplacian(i) = exp(0.5*Iwp*(coef_Kz(i) - real(k_L))*dz);
          }
      }
  }


  //! returns exp(i real(k) z) 
  Complex_wp SchrodingerProblem3D::GetPhaseKz(const Complex_wp& k, const Real_wp& z)
  {
    return Complex_wp(cos(real(k)*z), sin(real(k)*z));
  }
  
  
  //! computes exp(- imag(k) z) and its inverse
  void SchrodingerProblem3D::GetPhaseImagKz(const Complex_wp& k, const Real_wp& z,
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
  void SchrodingerProblem3D::GetLocalPhase(const Vector<Complex_wp>& kw, const Real_wp& z,
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
  
  
  //! setting 1-D interpolation grid
  void SchrodingerProblem3D::InitGrid(const Real_wp& x0, const Real_wp& xN, int Nd,
                                      int nb_points, GridInterpolation<Dimension1>& grid)
  {
    Mesh<Dimension1> mesh;
    mesh.CreateRegularMesh(xmin, xmax, nb_points+1, 1);
    
    grid.Init(x0, xN, Nd);
    grid.LocalizePoints(mesh);
  }
  
  
  //! setting 2-D interpolation grid
  void SchrodingerProblem3D::InitGrid(const Real_wp& x0, const Real_wp& xN, const Real_wp& y0,
                                      const Real_wp& yN, int Nx, int Ny)
  {
    InitGrid(x0, xN, Nx, nb_points_x, grid_x);
    InitGrid(y0, yN, Ny, nb_points_y, grid_y);
  }
  
  
  //! Computes the interpolation of xsol on the interpolation grid
  /*!
    \param[in] xsol vector to interpolate
    \param[out] Un interpolation of xsol on the interpolation grid
  */
  void SchrodingerProblem3D::GetInterpolateUn(VectComplex_wp& Un, const VectComplex_wp& xsol) const
  {
    // linear interpolation
    /* if ( (nb_points_x != grid_x.GetNbPointsGrid()-1)
       || (nb_points_y != grid_y.GetNbPointsGrid()-1))
       {
       cout << "Forget about linear interpolation" << endl;
       abort();
       }
    
       Un.Reallocate((nb_points_x+1)*(nb_points_y+1));
       for (int i = 0; i < nb_points_x; i++)
       for (int j = 0; j < nb_points_y; j++)
       Un(i*(nb_points_y+1) + j) = xsol(i*nb_points_y+j);
    
       for (int i = 0; i < nb_points_x; i++)
       Un(i*(nb_points_y+1) + nb_points_y) = xsol(i*nb_points_y);
    
       for (int j = 0; j < nb_points_y; j++)
       Un(nb_points_x*(nb_points_y+1) + j) = xsol(j);
    
       return; */
    
    Un.Reallocate(grid_x.GetNbPointsGrid()*grid_y.GetNbPointsGrid());
    for (int i = 0; i < grid_x.GetNbPointsGrid(); i++)
      for (int j = 0; j < grid_y.GetNbPointsGrid(); j++)
        {
          int num = i*grid_y.GetNbPointsGrid() + j;
          int ne_x = grid_x.GetElementNumber(i);
          int ne_y = grid_y.GetElementNumber(j);
          if ((ne_x >= 0) && (ne_y >= 0))
            {
              Real_wp x = grid_x.GetLocalCoordinate(i);
              Real_wp y = grid_y.GetLocalCoordinate(j);
              DISP(x); DISP(y);
              int ne_plus_x = (ne_x + 1)%nb_points_x;
              int ne_plus_y = (ne_y + 1)%nb_points_y;              
              Un(num) = (1.0-x)*(1.0-y)*xsol(ne_x*nb_points_y + ne_y)
                + x*(1.0-y)*xsol(ne_plus_x*nb_points_y + ne_y)
                + x*y*xsol(ne_plus_x*nb_points_y + ne_plus_y)
                + (1.0-x)*y*xsol(ne_x*nb_points_y + ne_plus_y);
            }
	  else
	    Un(num) = 0.0;
	}
  }
  

  //! Computes the interpolation of xsol on the interpolation grid for all the frequencies
  /*!
    \param[in] xsol vector to interpolate
    \param[out] Un interpolation of xsol on the interpolation grid
  */  
  void SchrodingerProblem3D::GetInterpolateFreqUn(Matrix<Complex_wp>& Un,
                                                  const VectComplex_wp& xsol) const  
  {
    int Ng = grid_x.GetNbPointsGrid()*grid_y.GetNbPointsGrid();
    Un.Reallocate(nb_points_time, Ng);
    
    VectComplex_wp Un_tmp(Ng), xsol_tmp;
    int N = nb_points_x*nb_points_y;
    for (int k = 0; k < nb_points_time; k++)
      {
        xsol_tmp.SetData(N, const_cast<Complex_wp*>(&xsol(k*N)));
        
        GetInterpolateUn(Un_tmp, xsol_tmp);
        
        for (int i = 0; i < Ng; i++)
          Un(k, i) = Un_tmp(i);
        
        xsol_tmp.Nullify();
      }
  }
  
  
  //! evaluates ProdEn = g(En) where g is the function such that du/dt = g(u)
  void SchrodingerProblem3D::EvaluateFunction(const Real_wp& zn, const Vector<Complex_wp>& En,
                                              Vector<Complex_wp>& ProdEn,
					      bool invert_mass, bool source)
  {
    if (type_approximation_Ez != EZ_NULL)
      {
        cout << "not implemented " << endl;
        abort();
      }
    
    if (type_polarization == ELLIPTIC)
      {
        cout << "not implemented " << endl;
        abort();
      }
    
    if (type_spectrum != SHIFTED_SINGLE)
      {
        cout << "not implemented " << endl;
        abort();
      }
    
    int N = nb_points_x*nb_points_y;
    if (nb_points_time == 1)
      {
        Vector<Complex_wp> phase;
        // on passe dans le domaine spatial reel x, y 
        if (type_model == UPPE)
          {
            phase.Reallocate(N);
            for (int i = 0; i < N; i++)
              {
                phase(i) = GetPhaseKz(coef_Kz(i), zn);
                ProdEn(i) = En(i)*phase(i);
              }
          }
        else
          Copy(En, ProdEn);
        
        fft2d.ApplyInverse(ProdEn);
        
        // on multiplie par -omega^2/c^2 gammaNL |E|^2 E
        Real_wp alpha = -gammaNL*square(omega_L/index.c0);
        for (int i = 0; i < En.GetM(); i++)
          ProdEn(i) = alpha*absSquare(ProdEn(i))*ProdEn(i);
        
        // on passe dans le domaine spectral
        fft2d.ApplyForward(ProdEn);
        
        if (type_model == UPPE)
          {
            // on multiplie par l'inverse de la phase
            for (int i = 0; i < N; i++)
              ProdEn(i) *= conj(phase(i));

            // on divise par 2 i k_z(omega)
            Complex_wp beta = 1.0/(2.0*Iwp);
            for (int i = 0; i < ProdEn.GetM(); i++)
              ProdEn(i) *= beta*coef_InvKz(i);
            
            for (int i = 0; i < N; i++)
              ProdEn(i) -= imag(coef_Kz(i))*En(i);
          }
        else
          {
            // on rajoute le laplacien, diagonal en domaine de Fourier
            for (int i = 0; i < N; i++)
              ProdEn(i) += coef_Gamma2(i)*En(i);
            
            // on divise par 2 i k(omega)
            Complex_wp beta = 1.0/(2.0*Iwp*k_L);
            Mlt(beta, ProdEn);
          }
        
        return;
      }
    
    if (local_enveloppe)
      {
        cout << "not implemented" << endl;
        abort();
      }

    // partie non-lineaire -gammaNL omega^2/c^2 F( |F^{-1}( E exp(i k(omega) z) )|^2
    //     F^{-1}( E exp(i k(omega) z) ) ) exp(-i k(omega) z)
    Vector<Complex_wp> phase;
    Vector<Real_wp> phaseImag, invPhaseImag;

    // on multiplie par la phase exp(i kz(omega) z)    
    if (type_model == UPPE)
      {
        phase.Reallocate(N*nb_points_time);
        //phaseImag.Reallocate(N*nb_points_time);
        //invPhaseImag.Reallocate(N*nb_points_time);
        for (int k = 0; k < nb_points_time; k++)
          for (int i = 0; i < N; i++)
            phase(k*N+i) = GetPhaseKz(coef_Kz(k*N+i), zn);
        
        for (int j = 0; j < En.GetM(); j++)
          ProdEn(j) = En(j)*phase(j);
      }
    else
      {
        phase.Reallocate(nb_points_time);
        for (int k = 0; k < nb_points_time; k++)
          phase(k) = GetPhaseKz(kwave(k), zn);
        
        for (int k = 0; k < nb_points_time; k++)
          for (int i = 0; i < N; i++)
            ProdEn(k*N+i) = En(k*N+i)*phase(k);
      }

    // on passe dans le domaine reel (en espace)
    Complex_wp* ProdEn_ptr = ProdEn.GetData();
    VectComplex_wp ProdEn_tmp;
    for (int k = 0; k < nb_points_time; k++)
      {
        ProdEn_tmp.SetData(N, ProdEn_ptr);            
        fft2d.ApplyInverse(ProdEn_tmp);            
        ProdEn_ptr += N; ProdEn_tmp.Nullify();
      }
    
    // on passe en temporel pour exprimer le terme non-lineaire
    VectComplex_wp EnTmp(nb_points_time + padding_time);
    for (int i = 0; i < N; i++)
      {
        EnTmp.Fill(0);
        for (int k = 0; k < nb_points_time/2; k++)
          {
            EnTmp(k) = ProdEn(k*N + i);
            int k2 = nb_points_time/2 + k;
            EnTmp(k2+padding_time) = ProdEn(k2*N + i);
          }
        
        // on passe en temporel
        fft.ApplyForward(EnTmp);
        Mlt(1.0/coef_fft, EnTmp);
            
        // on calcule |E|^2 E
        for (int k = 0; k < EnTmp.GetM(); k++)
          EnTmp(k) = absSquare(EnTmp(k))*EnTmp(k);
        
        // on repasse en Fourier
        fft.ApplyInverse(EnTmp);
        Mlt(coef_fft, EnTmp);
        
        // et on met a jour ProdEn
        for (int k = 0; k < nb_points_time/2; k++)
          {
            ProdEn(k*N + i) = EnTmp(k);
            int k2 = nb_points_time/2 + k;
            ProdEn(k2*N + i) = EnTmp(k2+padding_time);
          }
      }
    
    // on passe dans le domaine spectral (en espace)
    ProdEn_ptr = ProdEn.GetData();
    for (int k = 0; k < nb_points_time; k++)
      {
        ProdEn_tmp.SetData(N, ProdEn_ptr);            
        fft2d.ApplyForward(ProdEn_tmp);            
        ProdEn_ptr += N; ProdEn_tmp.Nullify();
      }
    
    // on multiplie par l'inverse de la phase (le conjugue ici)
    if (type_model == UPPE)
      {
        for (int j = 0; j < En.GetM(); j++)
          ProdEn(j) *= conj(phase(j));
      }
    else
      {
        for (int k = 0; k < nb_points_time; k++)
          for (int i = 0; i < N; i++)
            ProdEn(k*N+i) *= conj(phase(k));
      }
    
    // maintenant on calcule la fonction f(E) tel que dE/dz = f(E)
    Real_wp coefPol(1);
    if (type_polarization == LINEAR)
      coefPol = Real_wp(3)/4;
    else if(type_polarization == CIRCULAR)
      coefPol = 1.0 ;
    else
      coefPol = Real_wp(1)/4;
    
    VectComplex_wp coefNL(nb_points_time);
    for (int k = 0; k < nb_points_time; k++)
      {
        coefNL(k) = -coefPol*gammaNL*square(omega(k)/index.c0) /(Real_wp(2)*Iwp);
        if (type_model == UPPE_PARAXIAL)
          coefNL(k) /= kwave(k);
      }
    
    if (type_model == UPPE)
      {
        for (int k = 0; k < nb_points_time; k++)
          for (int i = 0; i < N; i++)
            {
              int p = k*N+i;
              ProdEn(p) = -imag(coef_Kz(p))*En(p) + coef_InvKz(p)*coefNL(k)*ProdEn(p);
            }
      }
    else
      {
        for (int k = 0; k < nb_points_time; k++)
          {
            Complex_wp beta = 1.0/kwave(k);
            for (int i = 0; i < N; i++)
              {
                int p = k*N+i;
                ProdEn(p) = coef_Gamma2(i)*beta*En(p) + coefNL(k)*ProdEn(p);
              }
          }
      }
    
    SetDirichletCondition(zn, 1, ProdEn);
  }
  

  //! sets u = 0 for forbidden frequencies
  void SchrodingerProblem3D::SetDirichletCondition(const Real_wp& z, int nb_deriv,
                                                   VectComplex_wp& x, Real_wp alpha)
  {
    Complex_wp zero; SetComplexZero(zero);
    int N = nb_points_x*nb_points_y;
    for (int i = 0; i < kwave_to_consider.GetM(); i++)
      {
        if (!kwave_to_consider(i))
          {
            for (int j = 0; j < N; j++)
              x(i*N + j) = zero;
          }
        else
          {
            for (int j = 0; j < N; j++)
              if (coef_Kz(i*N+j) == zero)
                x(i*N+j) = zero;
          }
      }
  }
  

  //! applying a resolution with a split operator
  void SchrodingerProblem3D::SolveSplitOperator(const Real_wp& z, const Real_wp& deltaz,
						const VectComplex_wp& En, VectComplex_wp& En_next,
                                                int op_num)
  {
    if (nb_points_time > 1)
      {
        cout << "Not implemented" << endl;
        abort();
      }
    
    if (op_num == 0)
      {
        // laplacian part
        Copy(En, En_next);
        for (int i = 0; i < En.GetM(); i++)
          En_next(i) *= coef_Laplacian(i);
      }
    else
      {
        Complex_wp beta = deltaz/(2.0*Iwp*k_L);
        Complex_wp alpha = -beta*gammaNL*square(omega_L/index.c0);
        
        Copy(En, En_next);
        fft2d.ApplyInverse(En_next);
        
        // non-linear part
        for (int i = 0; i < En_next.GetM(); i++)
          {
            Complex_wp coef_NL = exp(-alpha*absSquare(En(i)));
            En_next(i) *= coef_NL;
          }

        fft2d.ApplyForward(En_next);
      }
  }
  

  //! returns the harmonic pulse
  Complex_wp SchrodingerProblem3D::GetHarmonicPulse(const Real_wp& w)
  {
    // basic gaussian E(z=0, omega) = amplitude tau exp( -1/2 tau^2 (omega - omega_L)^2)
    //                               *polarization
    return exp(-0.5*square(tau_impulse*(w-omega_L)));
  }
  
    
  //! returns the initial condition
  void SchrodingerProblem3D::GetInitialCondition(Vector<Complex_wp>& E0)
  {
    Complex_wp zero, one; SetComplexZero(zero); SetComplexOne(one);
    
    // en temps, une gaussienne exp(-1/2 (t/tau)^2) -> exp( -tau^2 omega^2)
    VectComplex_wp pulse_time(omega.GetM());
    if (nb_points_time > 1)
      for (int i = 0; i < omega.GetM(); i++)
        pulse_time(i) = GetHarmonicPulse(omega(i));
    else
      pulse_time.Fill(one);
  
    int nb_comp = 1;
    if (type_polarization == ELLIPTIC)
      nb_comp++;
    
    if (type_approximation_Ez == EZ_UNKNOWN)
      nb_comp++;
    
    // en espace une gaussienne exp(-x^2/w^2)
    Real_wp alpha = 1.0/square(waist);
    int N = nb_points_x*nb_points_y;
    VectComplex_wp pulse_space(N);
    for (int i = 0; i < nb_points_x; i++)
      for (int j = 0; j < nb_points_y; j++)
        pulse_space(i*nb_points_y + j) = exp(-alpha*(step_x(i)*step_x(i) + step_y(j)*step_y(j)));
    
    // initial condition in spectral domain (in space)
    fft2d.ApplyForward(pulse_space);
    Complex_wp coef = amplitude_impulse*tau_impulse;
    if (nb_points_time == 1)
      coef = amplitude_impulse;
    
    E0.Reallocate(nb_comp*N*omega.GetM());
    E0.Fill(zero);    
    for (int i = 0; i < N; i++)
      {
        if (type_polarization == ELLIPTIC)
          {
            for (int k = 0; k < omega.GetM(); k++)
              {
                E0(k*N + i) = coef*pulse_space(i)*pulse_time(k)*polarization(0);
                E0((k+omega.GetM())*N + i) = coef*pulse_space(i)*pulse_time(k)*polarization(1);
              }
          }
        else
          {
            for (int k = 0; k < nb_points_time; k++)
              E0(k*N + i) = coef*pulse_space(i)*pulse_time(k)*polarization(0);
          }
      }
  }


  //! iterate E^n is given to this object  
  void SchrodingerProblem3D::GiveIterate(int n, const Real_wp& t, VectComplex_wp& En)
  {
    if (n == 0)
      return;
    
    int N = nb_points_x*nb_points_y;
    // testing overflow    
    if (check_overflow_scheme)
      for (int i = 0; i < En.GetM(); i++)
        if ((abs(En(i)) > 1e100) || isnan(abs(En(i))))
          {
            DISP(i); DISP(En(i));
            cout << "Instable scheme ?" << endl;
            abort();
          }
    
    // testing that the solution does not touch too much
    // the domain where epsilon is negative
    Real_wp max_En = 0, norm_En = En.GetNormInf();
    for (int i = 0; i < kwave.GetM(); i++)
      if (real(square(kwave(i))) < Real_wp(0))
        {
          for (int j = 0; j < N; j++)
            max_En = max(max_En, abs(En(i*N + j)));
        }
    
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
    int offset_last = kwave.GetM()-1;
    if (type_spectrum == SHIFTED_SINGLE)
      offset_last = kwave.GetM()/2-1;

    max_En = 0;
    //for (int k = 1; k <= En.GetM()/(kwave.GetM()*N); k++)
    if (nb_points_time > 1)
      {
        for (int i = 0; i < N; i++)
          max_En = max(max_En, abs(En(offset_last*N + i)));
        
        if (print_level >= 2)
          if (n%100 == 0)
            cout << "amplitude for the last frequency = " << max_En/norm_En << endl;
      }
    
    if (max_En/norm_En > threshold_last_frequency)
      {
        // DISP(max_En/amplitude_impulse);
        //cout << "The solution expands to a forbidden spectrum" << endl;
        //abort();
      }    
  }

  
  //! modifies iterates
  void SchrodingerProblem3D::ModifyIterate(int n, const Real_wp& z, VectComplex_wp& En)
  {
  }
  

  //! computes Y = f(X)
  void SchrodingerProblem3D
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectComplex_wp& X,
                               VectComplex_wp& Y, bool invert_mass, bool source)
  {
    EvaluateFunction(tn, X, Y);
  }
  
  
  //! returns the central frequency
  Real_wp SchrodingerProblem3D::GetFrequency() const
  {
    return omega_L/(2.0*pi_wp);
  }
    
  
  //! returns the final time for source
  Real_wp SchrodingerProblem3D::GetFinalTimeSource() const
  {
    return Real_wp(0);
  }
  

  //! extracts values that are above a given threshold
  template<class T> void SchrodingerProblem3D::
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
  void SchrodingerProblem3D::
  WriteOutputFile(const string& file_name, const string& file_name_time,
		  const Real_wp& zn, Vector<Complex_wp>& En)
  {
    if (print_level >= 4)
      glob_chrono.Start(VirtualTimer::OUTPUT);
    
    if (type_polarization == ELLIPTIC)
      {
	int N = omega.GetM()*nb_points_x*nb_points_y;
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
  bool SchrodingerProblem3D::NonLinearPhase() const
  {
    if ((type_model == UPPE) || (type_model == UPPE_PARAXIAL))
      return true;
    
    return false;
  }
  
  
  //! writes En in real space and Fourier space
  void SchrodingerProblem3D::
  WriteOutputFileScalar(const string& file_name, const string& file_name_time,
			const Real_wp& zn, Vector<Complex_wp>& En)
  {    
    Complex_wp zero; SetComplexZero(zero);
    Vector<Complex_wp> u(En);
    Mlt(E0_adim, u);
    int N = nb_points_time;
    int Ns = nb_points_x*nb_points_y;    
    Complex_wp vg, vloc;
    vg = index.GetGroupVelocity(omega_L);
    Real_wp omega0(0), omega1(0);

    // u is computed in real space
    Complex_wp phase;	
    for (int k = 0; k < N; k++)
      {
        if (type_model == UPPE)
          {
            for (int i = 0; i < Ns; i++)
              {
                phase = GetPhaseKz(coef_Kz(k*Ns + i), zn);
                u(Ns*k + i) *= phase;
              }
          }
        else
          {
            phase = GetPhaseKz(kwave(k), zn);
            for (int i = 0; i < Ns; i++)
              u(Ns*k + i) *= phase;
          }       
      }
    
    if (nb_points_time == 1)
      {
        fft2d.ApplyInverse(u);
      }
    else
      {
        Complex_wp* u_ptr = u.GetData();
        VectComplex_wp u_tmp;
        for (int k = 0; k < nb_points_time; k++)
          {
            u_tmp.SetData(Ns, u_ptr);            
            fft2d.ApplyInverse(u_tmp);            
            u_ptr += Ns; u_tmp.Nullify();
          }        
      }

    GridInterpolationFull<Dimension3> grid;
    //int type = GridInterpolationFull<Dimension3>::VOLUME;
    R3 ptMin_grid, ptMax_grid;
    int nx = nb_points_x, ny = nb_points_y, nz = 0;
    ptMin_grid(0) = step_x(0)*z0_adim; ptMin_grid(1) = step_y(0)*z0_adim;
    ptMax_grid(0) = step_x(nb_points_x-1)*z0_adim; ptMax_grid(1) = step_y(nb_points_y-1)*z0_adim;
    int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
    bool ascii = false;
    // Real_wp threshold = threshold_time_extraction*amplitude_impulse*E0_adim;
	
    int Nt = nb_points_time;
    if (factor_refinement_time > 1)
      Nt *= factor_refinement_time;
    
    // then computing u in real time
    Vector<VectComplex_wp> utime(Ns);
    for (int i = 0; i < Ns; i++)
      {
        utime(i).Reallocate(Nt);
        utime(i).Fill(zero);
      }
    
    for (int k = 0; k < N/2; k++)
      {
        phase = GetPhaseKz(-omega(k)/vg, zn);
        for (int i = 0; i < Ns; i++)
          {
            utime(i)(k) = u(Ns*k + i)*phase;
            u(Ns*k + i) = utime(i)(k);
          }
        
        phase = GetPhaseKz(-omega(k+N/2)/vg, zn);
        for (int i = 0; i < Ns; i++)
          {
            utime(i)(k+Nt-N/2) = u(Ns*(k+N/2) + i)*phase;
            u(Ns*(k+N/2) + i) = utime(i)(k+Nt-N/2);
          }
      }
    
    if (nb_points_time > 1)
      {
        if (factor_refinement_time > 1)
          {
            for (int i = 0; i < Ns; i++)
              {
                fft_refinement.ApplyForward(utime(i));    
                Mlt(1.0/coef_fft_refinement, utime(i));
              }
          }
        else
          {
            for (int i = 0; i < Ns; i++)
              {
                fft.ApplyForward(utime(i));    
                Mlt(1.0/coef_fft, utime(i));
              }
          }
        
        // computing utime = ifftshift(utime)
        VectComplex_wp utime_vec(Ns*Nt);
        for (int i = 0; i < Ns; i++)
          for (int k = 0; k < Nt/2; k++)
            {
              vloc = utime(i)(k);
              utime_vec(k*Ns + i) = utime(i)(Nt/2+k);
              utime_vec((Nt/2+k)*Ns + i) = vloc;
            }
        
        utime.Clear();
        Real_wp t0 = -Tmax + zn/real(vg), t1 = zn/real(vg) + time_interval(Nt/2-1);
    
        // using loadND format    
        nz = Nt; ptMin_grid(2) = t0*t0_adim; ptMax_grid(2) = t1*t0_adim;
        grid.SetVolumeOutput(ptMin_grid, ptMax_grid, nx, ny, nz);
        WriteMatlab(utime_vec, grid, file_name_time, double_prec, ascii);
      }
	
    // writing u in Fourier space
    Mlt(t0_adim, u);    
    
    // computing u = ifftshift(u)
    for (int i = 0; i < Ns; i++)
      for (int k = 0; k < N/2; k++)
        {
          vloc = u(Ns*k + i);
          u(Ns*k + i) = u((N/2+k)*Ns + i);
          u(Ns*(N/2+k) + i) = vloc;
        }
    
    if (nb_points_time > 1)
      {
        omega0 = omega(nb_points_time/2);
        omega1 = omega(nb_points_time/2-1);
      }
    else
      {
        omega0 = omega_L; omega1 = omega_L;
      }
    
    nz = N; ptMin_grid(2) = omega0/t0_adim; ptMax_grid(2) = omega1/t0_adim;
    grid.SetVolumeOutput(ptMin_grid, ptMax_grid, nx, ny, nz);
    WriteMatlab(u, grid, file_name, double_prec, ascii);
  }

}

#define MONTJOIE_FILE_SCHRODINGER_NON_LINEAR_3D_CXX
#endif
