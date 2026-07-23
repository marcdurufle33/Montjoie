#ifndef MONTJOIE_FILE_PHYSICAL_CONSTANT_CXX

namespace Montjoie
{

#ifndef SELDON_WITH_COMPILED_LIBRARY

#ifdef MONTJOIE_WITH_MPFR
  // physical constants
  const Real_wp PhysicalConstant::speed_light(299792458.0, MONTJOIE_DEFAULT_PRECISION);
  const Real_wp PhysicalConstant::invSquare_speed_light = Real_wp(1, MONTJOIE_DEFAULT_PRECISION)
    /(PhysicalConstant::speed_light*PhysicalConstant::speed_light);
  const Real_wp PhysicalConstant::elementary_charge(1.602176487e-19, MONTJOIE_DEFAULT_PRECISION);

  const Real_wp PhysicalConstant::mass_electron(9.10938215e-31, MONTJOIE_DEFAULT_PRECISION);
  const Real_wp PhysicalConstant::mass_proton(1.672621637e-27, MONTJOIE_DEFAULT_PRECISION);
  const Real_wp PhysicalConstant::mass_neutron(1.674927211e-27, MONTJOIE_DEFAULT_PRECISION);
  
  const Real_wp PhysicalConstant::epsilon0_permittivity(8.854187817620389e-12,
                                                        MONTJOIE_DEFAULT_PRECISION);
  
  const Real_wp PhysicalConstant::mu0_permeability(1.256637061435917e-6,
                                                   MONTJOIE_DEFAULT_PRECISION);
  
#else
  // physical constants
  const Real_wp PhysicalConstant::speed_light = 299792458.0;
  const Real_wp PhysicalConstant::invSquare_speed_light 
  = 1.0/(PhysicalConstant::speed_light*PhysicalConstant::speed_light);
  
  const Real_wp PhysicalConstant::elementary_charge = 1.602176487e-19;
  
  const Real_wp PhysicalConstant::mass_electron = 9.10938215e-31;
  const Real_wp PhysicalConstant::mass_proton = 1.672621637e-27;
  const Real_wp PhysicalConstant::mass_neutron = 1.674927211e-27;
  
  const Real_wp PhysicalConstant::epsilon0_permittivity =  8.854187817620389e-12;
  const Real_wp PhysicalConstant::mu0_permeability = 1.256637061435917e-6;
  
#endif

  const Real_wp PhysicalConstant::impedance0
  = sqrt(PhysicalConstant::mu0_permeability/PhysicalConstant::epsilon0_permittivity);
  
  int PhysicalConstant::adimensionalization(PhysicalConstant::ADIM_ONE);

  int PhysicalConstant::nb_max_indices(20);
#endif
  

  //! gets type of element and type of equation in the input file
  /*!
    \param[in] file_name name of the input file
    \param[out] type_element type of element
    \param[out] type_equation type of equation
  */
  void getElement_Equation(string& file_name, string & type_element, string& type_equation)
  {
    ifstream file_in(file_name.data());
    string description_field, keyword;
    if (!file_in.is_open())
      {
	cout<<"Data file not found"<<endl;
	cout<<file_name<<" doesn't exist !"<<endl;
	abort();
      }
    
    Vector<string> param; string ligne;
    while (!file_in.eof())
      {
	getlineM(file_in, ligne);
	
	// si une ligne commence avec #, c'est un commentaire 
	if ((ligne.size() > 0) && (ligne[0] != '#'))
	  {
	    StringTokenize(ligne, param, string(" \t"));
	    if (param.GetM() >= 3)
	      {
		description_field = param(0); keyword = param(2);
		if (!description_field.compare("TypeElement"))
		  type_element = keyword;
		
		if (!description_field.compare("TypeEquation"))
		  type_equation = keyword;
		
                if (description_field == "Timer")
                  {
                    if (keyword == "Accurate")
                      MontjoieTimer::default_timer = MontjoieTimer::ACCURATE_TIMER;
                    else if (keyword == "Real")
                      MontjoieTimer::default_timer = MontjoieTimer::REAL_TIMER;
                    else
                      MontjoieTimer::default_timer = MontjoieTimer::BASIC_TIMER;
                  }

		if (description_field == "NumberPhysicalMedia")
		  {
		    PhysicalConstant::nb_max_indices = to_num<int>(keyword);
		  }
		
		if (!description_field.compare("Adimensionalization"))
		  {
		    if (!keyword.compare("YES"))
		      PhysicalConstant::adimensionalization = PhysicalConstant::ADIM_YES;
		    else if (!keyword.compare("NO"))
		      PhysicalConstant::adimensionalization = PhysicalConstant::ADIM_NO;
		    else
		      PhysicalConstant::adimensionalization = PhysicalConstant::ADIM_ONE;
		  }
		
		if (is_seed_defined == false)
		  {
#ifdef SELDON_WITH_MPI
                    int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
                    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
                    int rank_proc(0);
#endif
                    
		    // We determine if there exists a seed in the initialization file
		    if (rank_proc == 0)
		      cout << "Looking for seed" << endl;
		    
		    if (!description_field.compare("Seed"))
		      {
			if (keyword == "Time")
			  Seed = (unsigned int)(time(0));
			
			Seed = to_num<int>(keyword);
			if (rank_proc == 0)
			  cout << "Seed found and equal to " << Seed << endl;
		      }
		    else
		      {
			Seed = 0;
		      }
		    
#ifdef SELDON_WITH_MPI
		    if (rank_proc == 0)
		      {
			cout << rank_proc << " Sending seed " << Seed
			     << " to other processors." << endl;
			for (int processor_to_send = 1;
			     processor_to_send < nb_proc;
			     processor_to_send++)
			  {
			    MPI_Ssend(&Seed, 1, MPI_INTEGER, processor_to_send, 667, MPI_COMM_WORLD);
			  }
		      }
		    else
		      {
			MPI_Status status;
			MPI_Recv(&Seed, 1, MPI_INTEGER, 0, 667, MPI_COMM_WORLD, &status);
		      }
#endif

		    // initializing srand with the provided Seed in the ini file
		    srand(Seed);
		    
		    is_seed_defined = true;
		  }
	      }
 	  }
      }
    
    file_in.close();
  }
  
  
  LorentzMaterial::LorentzMaterial()
  {
    // if adimensionalization is performed, c is equal to 1
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
      c0 = PhysicalConstant::speed_light;
    else
      c0 = 1.0;
    
    eps_inf = 1.0;
    alpha_damping = 0.0;
    forbid_negative_index = false;
  }
  
  //! returns the complex square root of n2
  Complex_wp LorentzMaterial::GetSqrtComplex(const Complex_wp& n2)
  {
    // positive imaginary part
    Complex_wp sqrt_n2 = sqrt(n2);
    if (imag(sqrt_n2) < 0)
      sqrt_n2 = -sqrt_n2;
    
    return sqrt_n2;
  }
  
    
  //! returns the wave number k(w) from pulsation w
  Complex_wp LorentzMaterial::GetWaveNumber(const Real_wp & w)
  {
    if (sigma_polarization.GetM() > 0)
      {
        Complex_wp n2 = eps_inf + Iwp/w*c0*alpha_damping;
        
        for (int k = 0; k < omega_polarization.GetM(); k++)
          n2 += alpha_polarization(k)
            / (Real_wp(1) - Iwp*w*sigma_polarization(k) - square(w/omega_polarization(k)));
        
        return GetSqrtComplex(n2)*w/c0;
      }
    
    // wave number k = n omega/c
    Real_wp n2 = eps_inf;
    // computing square of optical index
    for (int k = 0; k < omega_polarization.GetM(); k++)
      n2 += alpha_polarization(k) / (Real_wp(1) - square(w/omega_polarization(k)));
    
    if (forbid_negative_index)
      {
	if (n2 < 0)
	  {
	    cout << "The chosen pulsation is too high" << endl;
	    cout << "It is equal to " << w
                 << " but the square of the optical index is equal to " << n2 << endl;
	    cout << "Maybe choose an impulsion with a narrower spectrum" << endl;
	    abort();
	  }
      }

    if (alpha_damping != Real_wp(0))
      {
        Complex_wp n_2 = n2 + Iwp/w*c0*alpha_damping;
        return GetSqrtComplex(n_2)*w/c0;
      }
    
    return GetSqrtComplex(Complex_wp(n2, 0))*w/c0;
  }

  
  //! returns the group velocity from omega
  Complex_wp LorentzMaterial::GetGroupVelocity(const Real_wp & w)
  {
    // computed considering that the damping is small such that the group
    // velocity is slightly affected by its value
    // group velocity vg = 1/ k'(omega_L)
    // k'(omega) = (n'(omega) omega + n(omega)) / c0
    Real_wp n2 = eps_inf, sum1 = 0;
    for (int k = 0; k < omega_polarization.GetM(); k++)
      {
	n2 += alpha_polarization(k) / (Real_wp(1) - square(w/omega_polarization(k)));
	sum1 += alpha_polarization(k) 
          /square( (Real_wp(1) - square(w/omega_polarization(k)))) * w/square(omega_polarization(k));
      }
    
    Complex_wp sqrt_n2 = GetSqrtComplex(n2);
    Complex_wp ndot = sum1 / sqrt_n2;

    return c0/(ndot*w + sqrt_n2);
  }

  
  //! returns the dispersion of velocity group (sigma)
  Complex_wp LorentzMaterial::GetDispersion(const Real_wp & w)
  {
    // dispersion sigma = -k''(omega_L)/2
    // xi(omega) = \eps_inf + \sum alpha_k / (1 - (w/w_k)^2)
    // sum1 = 1/2 xi(omega)' = \sum omega / (omega_k)^2 /  [ (1 - (w/w_k)^2) ]^2
    // sum2 = xi(omega)'' = \sum 2 / omega_k^2
    //             [ (1 - (w/w_k)^2) ]^2 + 8 (omega/omega_k)^2 / omega_k [(1 - (w/w_k)^2) ]^3
    Real_wp n2 = eps_inf, sum1(0), sum2(0), denom, denom2, denom3, wk, wk2;
    for (int k = 0; k < omega_polarization.GetM(); k++)
      {
	wk = omega_polarization(k); wk2 = wk*wk;
	denom = (Real_wp(1) - square(w/wk) );
	denom2 = denom*denom; denom3 = denom2*denom;
	n2 += alpha_polarization(k) / denom;
	sum1 += alpha_polarization(k) * w / (wk2*denom2);
	sum2 += alpha_polarization(k) * ( 2.0 / (wk2*denom2) + 8.0*square(w/wk)/(wk2*denom3));
      }
    
    Complex_wp sqrt_n2 = GetSqrtComplex(n2);
    Complex_wp ndot = sum1 / sqrt_n2;
    Complex_wp ndotdot = 0.5*sum2 / sqrt_n2 - sum1*sum1/(n2*sqrt_n2);

    return -(ndotdot*w + Real_wp(2.0)*ndot) / (2.0*c0);
  }
  
  
  //! sets conductivity of the material
  void LorentzMaterial::SetSigma(const Real_wp& sigma, const Real_wp& z0_adim)
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
      alpha_damping = sigma*PhysicalConstant::mu0_permeability;
    else if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      alpha_damping = sigma*PhysicalConstant::mu0_permeability
        *PhysicalConstant::speed_light*z0_adim;
    else
      alpha_damping = sigma;
    
    alpha_damping *= c0;
  }
  

  //! fills coefficients eps_inf, omega_k, sigma_k, alpha_k depending on the material
  /*!
    \param[in] material if equal to a float number, the coefficients alpha_k, omega_k, sigma_k
      are assumed to be given in parameters
    \param[in] parameters additional parameters describing the material
    epsilon = eps_inf + \sum_k \alpha_k / (1 - i omega / (omega_k^2 T_k) - (omega/omega_k)^2)
   */
  void LorentzMaterial::SetRefractiveIndex(const string& material,
                                           const Vector<string>& parameters)
  {
    sigma_polarization.Clear();
    if (material == "Quartz")
      {
	// the index is equal to 1.28604141 + 1.07044083 lambda^2 / (lambda^2 - 1.00585997e-2)
        //     + 1.10202242 lambda^2 / (lambda^2 - 100)
	// each term is equal to alpha lambda^2 / (lambda^2 - beta) 
        //   = alpha / (1 - beta/(2 pi c)^2 omega^2)
	// then omega_k = 2 pi c / sqrt(beta)
	omega_polarization.Reallocate(2);
	alpha_polarization.Reallocate(2);
	eps_inf = 1.28604141;
	alpha_polarization(0) = 1.07044083;
	omega_polarization(0) = 2.0*pi_wp*PhysicalConstant::speed_light
          / (1e-6*sqrt(Real_wp(1.00585997e-2)));
	
        alpha_polarization(1) = 1.10202242;
	omega_polarization(1) = 2.0*pi_wp*PhysicalConstant::speed_light
          / (1e-6*sqrt(Real_wp(100.0)));
      }
    else if (material == "Silica")
      {
	omega_polarization.Reallocate(3);
	alpha_polarization.Reallocate(3);
	eps_inf = 1.0;
	alpha_polarization(0) = 0.6961663;
	omega_polarization(0) = 2.0*pi_wp*PhysicalConstant::speed_light / (1e-6*0.0684043);
	alpha_polarization(1) = 0.4079426;
	omega_polarization(1) = 2.0*pi_wp*PhysicalConstant::speed_light / (1e-6*0.1162414);
	alpha_polarization(2) = 0.8974794;
	omega_polarization(2) = 2.0*pi_wp*PhysicalConstant::speed_light / (1e-6*9.896161);
      }
    else if (material == "DampedMaterial")
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In GetRefractiveIndex, at least three parameters are needed" << endl;
            cout << "RefractiveIndex = DampedMaterial eps 0" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

       	eps_inf = to_num<Real_wp>(parameters(1));
	int N = to_num<int>(parameters(2));
	omega_polarization.Reallocate(N);
	alpha_polarization.Reallocate(N);
        sigma_polarization.Reallocate(N);
	for (int k = 0; k < N; k++)
	  {
	    alpha_polarization(k) = to_num<Real_wp>(parameters(3+3*k));
	    Real_wp lambda_res = to_num<Real_wp>(parameters(4+3*k));
	    omega_polarization(k) = 2.0*pi_wp*PhysicalConstant::speed_light / (1e-6*lambda_res);
            Real_wp Tk = to_num<Real_wp>(parameters(5+3*k));
            sigma_polarization(k) = 1.0/(square(omega_polarization(k))*Tk);
	  }
      }
    else
      {
	eps_inf = to_num<Real_wp>(parameters(0));
	int N = to_num<int>(parameters(1));
	omega_polarization.Reallocate(N);
	alpha_polarization.Reallocate(N);
	for (int k = 0; k < N; k++)
	  {
	    alpha_polarization(k) = to_num<Real_wp>(parameters(2+2*k));
	    Real_wp lambda_res = to_num<Real_wp>(parameters(3+2*k));
	    omega_polarization(k) = 2.0*pi_wp*PhysicalConstant::speed_light / (1e-6*lambda_res);
	  }
      }
  }


  /*****************
   * DrudeMaterial *
   *****************/


  //! default constructor
  DrudeMaterial::DrudeMaterial()
  {
    eps_inf = Real_wp(1);
    type_mode = NONE;
  }
  

  //! Sets indices from parameters
  /*!
    \param[in] nb starting index for parameters
    \param[in] parameters parameters defining the material
    \param[in] mode TE or TM mode
   */
  void DrudeMaterial::SetInputData(int nb, const Vector<string>& parameters,
				   const string& mode)
  {
    // the mode is affected
    if (mode == "TE")
      type_mode = TE;
    else
      type_mode = TM;

    // type of material is retrieved
    string material(parameters(nb));
    nb++;

    // then other parameters are stored in param
    Vector<string> param(parameters.GetM()-nb);
    for (int i = nb; i < parameters.GetM(); i++)
      param(i-nb) = parameters(i);

    // the index is initialized with the given material and parameters
    SetRefractiveIndex(material, param);
  }


  //! Sets the material with given parameters
  void DrudeMaterial::SetRefractiveIndex(const string& material,
					 const Vector<string>& parameters)
  {
    VectReal_wp omega_p, omega_0, sigma;
    // for known materials, indices are filled automatically
    if (material == "Silver")
      {
	eps_inf = 1.0;
	omega_p.Reallocate(1);
	omega_0.Reallocate(1);
	gamma.Reallocate(1);
        sigma.Reallocate(1); sigma.Zero();
	omega_p(0) = 2.0*pi_wp*PhysicalConstant::speed_light / (0.138e-6);
	omega_0(0) = 0.0;
	gamma(0) = 0.0023*omega_p(0);

        if (parameters.GetM() > 0)
          {
            omega_p.Clear(); omega_0.Clear(); sigma.Clear();
            int L = to_num<int>(parameters(0));
            VectComplex_wp Om(L), Sig(L);
            Real_wp gam(0), sig(0);
            if (L == 0)
              {
                eps_inf = 3.906;
                gam = 8.66887e-02;
                sig = 8.8346e+02;                
              }
            else if (L == 1)
              {
                eps_inf = 2.65057e+00;
                gam = 8.66887e-02;
                sig = 8.83462e+02;
                Om(0) = Complex_wp(4.5236e+00, -5.1847e-01);
                Sig(0) = Complex_wp(4.7539e-01, 2.1011e+00);
              }
            else if (L == 2)
              {
                eps_inf = 1.2135;
                gam = 8.66887e-2;
                sig = 8.8346e+02;
                Om(0) = Complex_wp(4.249, -2.1153);
                Sig(0) = Complex_wp(3.5366, 4.5884);
                Om(1) = Complex_wp(4.0097, -3.0102e-01);
                Sig(1) = Complex_wp(7.4248e-01, -1.9950e-01);
              }

            FillPadeExpansionEV(Om, Sig, gam, sig);
          }
      }
    else if (material == "Gold")
      {
	eps_inf = 1.0;
	omega_p.Reallocate(1);
	omega_0.Reallocate(1);
	gamma.Reallocate(1);
        sigma.Reallocate(1); sigma.Zero();
	omega_p(0) = 2.0*pi_wp*PhysicalConstant::speed_light / (0.1496e-6);
	omega_0(0) = 0.0;
	gamma(0) = 0.0112*omega_p(0);
        
        if (parameters.GetM() > 0)
          {
            omega_p.Clear(); omega_0.Clear(); sigma.Clear();
            int L = to_num<int>(parameters(0));
            VectComplex_wp Om(L), Sig(L);
            Real_wp gam(0), sig(0);
            if (L == 0)
              {
                gam = 92.8e-3;
                sig = 744;
                eps_inf = 1.0;
              }
            else if (L == 1)
              {
                eps_inf = 3.9199;
                gam = 0.0893;
                sig = 875.79;
                Om(0) = Complex_wp(2.7326, -0.69021);
                Sig(0) = Complex_wp(3.0701, 2.9306);                
              }
            else if (L == 2)
              {
                eps_inf = 2.6585;
                gam = 0.07247;
                sig = 1056.9;
                
                Om(0) = Complex_wp(2.5509, -0.27427);
                Sig(0) = Complex_wp(0.57604, 0.18443);
                Om(1) = Complex_wp(2.8685, -1.2195);
                Sig(1) = Complex_wp(4.1891, 4.2426);
              }
            
            FillPadeExpansionEV(Om, Sig, gam, sig);
          }
      }
    else if (material == "Silicon")
      {
        if (parameters.GetM() > 0)
          {
            eps_inf = 1.12648273;
            int L = to_num<int>(parameters(0));
            VectComplex_wp Om(1), Sig(1);
            if (L == 1)
              {
                Om(0) = Complex_wp(3.95095,-0.190893);
                Sig(0) = Complex_wp(2.17595,20.77585);
              }
            else
              {
                cout << "Only 1 accepted " << endl;
                abort();
              }
            
            FillPadeExpansionEV(Om, Sig);
          }
        else
          {
            eps_inf = 0.81568;
            VectComplex_wp Om(5), Sig(5);
            Om(0) = Complex_wp(3.3736, -0.11402);
            Sig(0) = Complex_wp(1.6934, 2.084);
            Om(1) = Complex_wp(3.6519, -0.52378);
            Sig(1) = Complex_wp(5.2573, 8.0106);
            Om(2) = Complex_wp(4.287, -0.21116);
            Sig(2) = Complex_wp(-1.7164, 5.9939);
            Om(3) = Complex_wp(5.3188, -0.18434);
            Sig(3) = Complex_wp(-0.00528, 0.32911);
            Om(4) = Complex_wp(5.5064, -1.7892);
            Sig(4) = Complex_wp(-3.8438, 6.9298);
            FillPadeExpansionEV(Om, Sig);
          }
      }
    else if (material == "GaAs")
      {
        eps_inf = -0.54651;
        VectComplex_wp Om(5), Sig(5);
        Om(0) = Complex_wp(1.4377, -0.05948);
        Sig(0) = Complex_wp(0.01981, 0.01122);
        Om(1) = Complex_wp(2.7229, -1.2972);
        Sig(1) = Complex_wp(7.8336, 8.3274);
        Om(2) = Complex_wp(2.8922, -0.23992);
        Sig(2) = Complex_wp(2.706, 1.616);
        Om(3) = Complex_wp(4.5222, -0.42072);
        Sig(3) = Complex_wp(2.1137, 4.6445);
        Om(4) = Complex_wp(4.9278, -0.19972);
        Sig(4) = Complex_wp(-1.243, 1.4424);
        FillPadeExpansionEV(Om, Sig);
      }
    else if (material == "Germanium")
      {
        if (parameters.GetM() > 0)
          {
            eps_inf = 2.37031;
            VectComplex_wp Om(2), Sig(2);
            Om(0) = Complex_wp(2.22291,-1.02029); Om(1) = Complex_wp(2.07938,-0.19684);
            Sig(0) = Complex_wp(10.767, 12.8471); Sig(1) = Complex_wp(2.41218,0.64685);
            FillPadeExpansionEV(Om, Sig);
          }
        else
          {
            eps_inf = 0.79842;
            VectComplex_wp Om(5), Sig(5);
            Om(0) = Complex_wp(1.168, -0.33778);
            Sig(0) = Complex_wp(0.47159, 0.01002);
            Om(1) = Complex_wp(2.174, -0.28077);
            Sig(1) = Complex_wp(3.2926, 4.1239);
            Om(2) = Complex_wp(3.781, -1.1461);
            Sig(2) = Complex_wp(0.86584, 18.898);
            Om(3) = Complex_wp(4.3232, -0.20006);
            Sig(3) = Complex_wp(-1.7377, 2.5278);
            Om(4) = Complex_wp(5.6442, -0.41214);
            Sig(4) = Complex_wp(0.10451, 1.0292);
            FillPadeExpansionEV(Om, Sig);
          }
      }
    else if (material == "ManualPadeEV")
      {
        eps_inf = to_num<Real_wp>(parameters(0));
        int N = to_num<int>(parameters(1)); int nb = 2;
        VectComplex_wp Om(N), Sig(N);
        for (int i = 0; i < N; i++)
          {
            Om(i) = to_num<Complex_wp>(parameters(nb));
            Sig(i) = to_num<Complex_wp>(parameters(nb+1));
            nb += 2;
          }
        
        FillPadeExpansionEV(Om, Sig);
      }
    else if (material == "ManualRootEV")
      {
        eps_inf = to_num<Real_wp>(parameters(0));
        int N = to_num<int>(parameters(1)); int nb = 2;
        VectComplex_wp numer(N), denom(N);
        Real_wp xmin(1e300), xmax(-1e300);
        for (int i = 0; i < N; i++)
          {
            if (parameters.GetM() < nb+4)
              {
                cout << "Not enough parameters " << endl;
                DISP(i); DISP(nb); DISP(parameters);
                abort();
              }
            
            Real_wp x1 = to_num<Real_wp>(parameters(nb));
            Real_wp y1 = to_num<Real_wp>(parameters(nb+1));
            Real_wp x2 = to_num<Real_wp>(parameters(nb+2));
            Real_wp y2 = to_num<Real_wp>(parameters(nb+3));
            numer(i) = Complex_wp(x1, y1);
            denom(i) = Complex_wp(x2, y2);
            xmin = min(xmin, x2); xmax = max(xmax, x2);
            xmin = min(xmin, -x2); xmax = max(xmax, -x2);
            nb += 4;
          }

        // on recupere les coefs de Pade en faisant une integration
        // pour eviter le mauvais conditionnement du VDM
        VectReal_wp points, poids;
        ComputeGaussLegendre(points, poids, 2*N);

        Matrix<Complex_wp, Symmetric, RowSymPacked> Mass(2*N+1, 2*N+1);
        VectComplex_wp F(2*N+1);
        Mass.Zero(); F.Zero();
        for (int k = 0; k < points.GetM(); k++)
          {
            Real_wp x = xmin + (xmax-xmin)*points(k);
            Real_wp w = (xmax-xmin)*poids(k);

            // on evalue tous les polynomes au point x
            VectComplex_wp phi(2*N+1);
            phi.Fill(Complex_wp(1, 0));
            for (int i = 0; i < N; i++)
              {
                Complex_wp t1 = x-denom(i), t2 = x+conjugate(denom(i));
                Complex_wp facteur = t1*t2;
                phi(0) *= facteur;
                for (int j = 0; j < N; j++)
                  {
                    if (j == i)
                      {
                        phi(2*j+1) *= t2;
                        phi(2*j+2) *= t1;
                      }
                    else
                      {
                        phi(2*j+1) *= facteur;
                        phi(2*j+2) *= facteur;
                      }
                  }
              }

            // on evalue le numerateur
            Complex_wp psi(1, 0);
            for (int i = 0; i < N; i++)
              psi *= (x-numer(i))*(x+conjugate(numer(i)));

            // on incremente la masse et F
            for (int i = 0; i <= 2*N; i++)
              {
                for (int j = i; j <= 2*N; j++)
                  Mass(i, j) += w*phi(i)*phi(j);

                F(i) += w*phi(i)*psi;
              }            
          }

        Vector<int> pivot;
        GetLU(Mass, pivot);
        SolveLU(Mass, pivot, F);
        
        VectComplex_wp Sig(N);
        for (int i = 0; i < N; i++)
          Sig(i) = -Iwp*eps_inf*F(2*i+1);

        DISP(Sig); DISP(denom);
        // les coefficients de Pade sont dans F
        FillPadeExpansionEV(denom, Sig);        
      }
    else if (material == "Manual")
      {
        // manual parameters assuming that sigma = 0
	eps_inf = to_num<Real_wp>(parameters(0));
	int N = to_num<int>(parameters(1));
	omega_p.Reallocate(N);
	omega_0.Reallocate(N);
	gamma.Reallocate(N);
        sigma.Reallocate(N); sigma.Zero();
	int num = 2;
	for (int k = 0; k < N; k++)
	  {
	    // constants are given directly
	    omega_p(k) = to_num<Real_wp>(parameters(num));
            omega_0(k) = to_num<Real_wp>(parameters(num+1));
	    gamma(k) = to_num<Real_wp>(parameters(num+2));
	    
	    num += 3;
	  }
      }
    else if (material == "ManualSigma")
      {
        // manual parameters with non-null sigma
	eps_inf = to_num<Real_wp>(parameters(0));
	int N = to_num<int>(parameters(1));
	omega_p.Reallocate(N);
	omega_0.Reallocate(N);
	gamma.Reallocate(N);
        sigma.Reallocate(N);
	int num = 2;
	for (int k = 0; k < N; k++)
	  {
	    // constants are given directly
	    omega_p(k) = to_num<Real_wp>(parameters(num));
            omega_0(k) = to_num<Real_wp>(parameters(num+1));
	    gamma(k) = to_num<Real_wp>(parameters(num+2));
            sigma(k) = to_num<Real_wp>(parameters(num+3));
	    
	    num += 4;
	  }
      }
    else if (material == "ManualWaveLength")
      {
        // manual parameters, omega_p, omega_0, gamma and sigma are given in wavelengthes
        // we compute 2*pi*c / lambda to obtain the pulsation
	eps_inf = to_num<Real_wp>(parameters(0));
	int N = to_num<int>(parameters(1));
	omega_p.Reallocate(N);
	omega_0.Reallocate(N);
	gamma.Reallocate(N);
        sigma.Reallocate(N);
	int num = 2;
	for (int k = 0; k < N; k++)
	  {
	    // constants are given in wavelengths
            // omega = 2 pi c / lambda
	    omega_p(k) = 2.0*pi_wp*PhysicalConstant::speed_light
	      / to_num<Real_wp>(parameters(num));
	    
	    omega_0(k) = 2.0*pi_wp*PhysicalConstant::speed_light
	      / to_num<Real_wp>(parameters(num+1));
	    
	    gamma(k) = 2.0*pi_wp*PhysicalConstant::speed_light
	      / to_num<Real_wp>(parameters(num+2));

            sigma(k) = 2.0*pi_wp*PhysicalConstant::speed_light
	      / to_num<Real_wp>(parameters(num+3));
	    
	    num += 4;
	  }
      }
    else
      {
        cout << "Unknown material = " << material << endl;
        abort();
      }

    if (omega_p.GetM() > 0)
      {
        // on remplit eps_omega_p2, omega_02 et eps_sigma
        int N = omega_p.GetM();
        eps_omega_p2.Reallocate(N);
        omega_02.Reallocate(N);
        eps_sigma.Reallocate(N);
        for (int i = 0; i < N; i++)
          {
            eps_omega_p2(i) = eps_inf*square(omega_p(i));
            omega_02(i) = square(omega_0(i));
            eps_sigma(i) = eps_inf*sigma(i);
          }
      }

    /*VectReal_wp omEv; VectComplex_wp epsEv;
    int N = 1000;
    Linspace(Real_wp(0.5), Real_wp(7), N, omEv);

    Real_wp e = 1.602176634e-19;
    Real_wp hbar = 6.62607015e-34/(2.0*pi_wp);
    Real_wp coef = e/hbar;
    epsEv.Reallocate(N);
    for (int i = 0; i < N; i++)
      epsEv(i) = GetEpsilon(omEv(i)*coef);
    
    omEv.WriteText("omega.dat");
    epsEv.Write("eps.dat");*/
  }


  //! Fills coefficients omega_02, eps_omega_p2, eps_sigma and gamma from complex coefficients Omega and Sigma
  /*!
    we consider the following Pade expansion
    eps = eps_inf + \sum_k i \sigma_k / (omega - Omega_k) + i conj(sigma_k) / (omega + conj(Omega_k))    
   */
  void DrudeMaterial::FillPadeExpansionEV(VectComplex_wp& Om, VectComplex_wp& Sig,
                                          Real_wp gam, Real_wp sig)
  {
    int L = Om.GetM();
    eps_omega_p2.Reallocate(L);
    omega_02.Reallocate(L);
    gamma.Reallocate(L);
    eps_sigma.Reallocate(L);

    // Om and Sig are given in electron-volts, conversion in SI
    Real_wp e = 1.602176634e-19;
    Real_wp hbar = 6.62607015e-34/(2.0*pi_wp);
    Real_wp coef = e/hbar;
    Mlt(coef, Om);
    Mlt(coef, Sig);
    for (int k = 0; k < L; k++)
      {
        eps_omega_p2(k) = 2.0*imagpart(Sig(k)*conjugate(Om(k)));
        omega_02(k) = absSquare(Om(k));
        gamma(k) = -2.0*imagpart(Om(k));
        eps_sigma(k) = 2.0*realpart(Sig(k));
      }

    if (sig != Real_wp(0))
      {
        gam *= coef;
        sig *= coef;
        eps_omega_p2.PushBack(gam*sig);
        omega_02.PushBack(Real_wp(0));
        gamma.PushBack(gam);
        eps_sigma.PushBack(Real_wp(0));
      }    
  }
  

  //! returns zeros of permittivity or solves epsilon = cte
  VectComplex_wp DrudeMaterial::GetZeroOmega(Real_wp cte) const
  {
    UnivariatePolynomial<Complex_wp> Denom, Pol, Numer;
    Vector<UnivariatePolynomial<Complex_wp> > facteur(gamma.GetM());
    for (int i = 0; i < gamma.GetM(); i++)
      {
        Pol.SetOrder(2);
        Pol(0) = -omega_02(i);
        Pol(1) = Iwp*gamma(i);
        Pol(2) = Real_wp(1);
        if (i == 0)
          Denom = Pol;
        else
          Denom *= Pol;

        facteur(i) = Pol;
      }

    Numer = eps_inf*Denom - cte*Denom;
    for (int i = 0; i < gamma.GetM(); i++)
      {
        Pol.SetOrder(1);
        Pol(0) = eps_omega_p2(i);
        Pol(1) = -Iwp*eps_sigma(i);
        for (int j = 0; j < gamma.GetM(); j++)
          if (j != i)
            Pol *= facteur(j);

        Numer -= Pol;
      }

    VectComplex_wp R;
    SolvePolynomialEquation(Numer, R);
    return R;
  }

  
  //! returns poles of permittivity
  VectComplex_wp DrudeMaterial::GetPoleOmega() const
  {
    VectComplex_wp pole;
    for (int i = 0; i < gamma.GetM(); i++)
      {
        if (omega_02(i) == Real_wp(0))
          {
            pole.PushBack(Complex_wp(0, 0));
            pole.PushBack(Complex_wp(0, -gamma(i)));
          }
        else
          {
            Real_wp delta = -square(gamma(i)) + 4.0*omega_02(i);
            Complex_wp z1(0), z2(0);
            if (delta < 0)
              {
                Real_wp a = sqrt(-delta);
                z1 = Complex_wp(0, -0.5*(gamma(i) + a));
                z2 = Complex_wp(0, -0.5*(gamma(i) - a));                
              }
            else
              {
                Real_wp a = sqrt(delta);
                z1 = Complex_wp(-0.5*a, -0.5*gamma(i));
                z2 = Complex_wp(0.5*a, -0.5*gamma(i));
              }
            
            pole.PushBack(z1); pole.PushBack(z2);
          }
      }
    
    return pole;
  }


  //! performs adimensionalizations (indices are divided by speed light
  void DrudeMaterial::Adimensionalize(const Real_wp& L0)
  {
    Real_wp coef(1);
    if (PhysicalConstant::adimensionalization != PhysicalConstant::ADIM_NO)
      {
	coef = L0/PhysicalConstant::speed_light;
	Mlt(coef*coef, eps_omega_p2);
	Mlt(coef*coef, omega_02);
	Mlt(coef, gamma);
        Mlt(coef, eps_sigma);
      }

    DISP(eps_omega_p2); DISP(eps_inf); DISP(omega_02); DISP(gamma); DISP(eps_sigma);
  }
  

  //! returns eps(omega)
  Complex_wp DrudeMaterial::GetEpsilon(const Real_wp& omega) const
  {
    Complex_wp eps(eps_inf, 0);
    for (int i = 0; i < gamma.GetM(); i++)
      eps -= (eps_omega_p2(i) - Iwp*eps_sigma(i)*omega) / (omega*omega - omega_02(i) + Iwp*gamma(i)*omega);
    
    return eps;
  }


  //! returns eps(omega) - eps_inf
  Complex_wp DrudeMaterial::GetDeltaEpsilon(const Complex_wp& omega) const
  {
    Complex_wp deps(0, 0);
    for (int i = 0; i < gamma.GetM(); i++)
      deps -= (eps_omega_p2(i) - Iwp*eps_sigma(i)*omega) / (omega*omega - omega_02(i) + Iwp*gamma(i)*omega);
    
    return deps;
  }
  

  //! returns true if a TE or TM mode is set
  bool DrudeMaterial::IsEnabled() const
  {
    if (type_mode == NONE)
      return false;
    
    return true;
  }
  

  //! modifies rho, mu of Helmholtz equation 
  void DrudeMaterial::ModifyCoefficientHelmholtz(const Real_wp& omega, Complex_wp& rho,
						 Complex_wp& sigma, TinyMatrix<Complex_wp, Symmetric, 2, 2>& mu) const
  {
    if (type_mode == NONE)
      return;

    SetComplexZero(sigma);
    if (type_mode == TM)
      {
	rho = GetEpsilon(omega);
	mu.SetIdentity();
      }
    else
      {
	SetComplexOne(rho);
	Complex_wp eps = GetEpsilon(omega);
	mu.SetDiagonal(Real_wp(1)/eps);
      }
  }

  
  void DrudeMaterial::ModifyCoefficientHelmholtz(const Real_wp& omega, Complex_wp& rho,
						 Complex_wp& sigma, TinyMatrix<Complex_wp, Symmetric, 3, 3>& mu) const
  {
    if (type_mode == NONE)
      return;

    SetComplexZero(sigma);
    if (type_mode == TM)
      {
	rho = GetEpsilon(omega);
	mu.SetIdentity();
      }
    else
      {
	SetComplexOne(rho);
	Complex_wp eps = GetEpsilon(omega);
	mu.SetDiagonal(Real_wp(1)/eps);
      }
  }

  
  void DrudeMaterial::ModifyCoefficientHelmholtz(const Real_wp& omega, Real_wp& rho,
						 Real_wp& sigma, TinyMatrix<Real_wp, Symmetric, 2, 2>& mu) const
  {
    return;
  }
  

  void DrudeMaterial::ModifyCoefficientHelmholtz(const Real_wp& omega, Real_wp& rho,
						 Real_wp& sigma, TinyMatrix<Real_wp, Symmetric, 3, 3>& mu) const
  {
    return;
  }

  
  void DrudeMaterial::ModifyCoefficientMaxwell(const Real_wp& omega, TinyMatrix<Complex_wp, Symmetric, 3, 3>& epsilon,
                                               TinyMatrix<Complex_wp, Symmetric, 3, 3>& sigma, TinyMatrix<Complex_wp, Symmetric, 3, 3>& mu) const
  {
    if (type_mode == NONE)
      return;

    SetComplexZero(sigma);
    Complex_wp eps = GetEpsilon(omega);
    if (type_mode == TM)
      {
        epsilon.SetDiagonal(eps);
        mu.SetIdentity();
      }
    else
      {
        epsilon.SetIdentity();
        mu.SetDiagonal(eps);
      }
  }


  void DrudeMaterial::ModifyCoefficientMaxwell(const Real_wp& omega, TinyMatrix<Real_wp, Symmetric, 3, 3>& epsilon,
                                               TinyMatrix<Real_wp, Symmetric, 3, 3>& sigma, TinyMatrix<Real_wp, Symmetric, 3, 3>& mu) const
  {
    return;
  }
  
}

#define MONTJOIE_FILE_PHYSICAL_CONSTANT_CXX
#endif
