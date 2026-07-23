#ifndef MONTJOIE_FILE_GENERIC_WATER_WAVE_EQUATION_CXX

namespace Montjoie
{
  
  //! default constructor
  GenericWaterWaveEquation::GenericWaterWaveEquation()
  {
    // diverses variables
    h0 = 1.0;
    type_model = GENTLE;    
    InitEpsilon(0.1, 0.1, 0.1, 0.1);
    SetInitialCondition(GAUSSIAN, SINUS, -10.0, 0.1);
  }
    
  
  //! copying parameters of a similar class
  template<class TypeEq>
  void GenericWaterWaveEquation::Copy(const TypeEq& var)
  {
    a = var.a;
    h0 = var.h0;
    epsilon = var.epsilon;
    mu = var.mu;
    lambda = var.lambda;
    beta = var.beta;
    alpha = var.alpha;
    slope_bottom = var.slope_bottom;
    yintercept = var.yintercept;
    c1 = var.c1;
    nu = var.nu;
    center_init = var.center_init;
    freq_init = var.freq_init;
    type_bottom = var.type_bottom;
    freq_bottom = var.freq_bottom;
    type_initial = var.type_initial;
    type_model = var.type_model;
  }

  
  //! setting the value of epsilon, mu, alpha, beta
  void GenericWaterWaveEquation::InitEpsilon(double eps, double mu_, double alpha_, double beta_)
  {
    alpha = alpha_;
    beta = beta_;    
    epsilon = eps;
    mu = mu_;
    
    a = h0*epsilon;
    lambda = h0/sqrt(mu);    
    nu = mu/epsilon;
    c1 = 0.5;
    // c = 1 + c1 epsilon
    // c1 = 2.0*nu/3.0*square(2.0*pi_wp/lambda);
    
  }
  
  
  //! setting the initial condition for the wave
  void GenericWaterWaveEquation
  ::SetInitialCondition(int type_init_, int type_bottom_, double x0, double freq)
  {
    type_bottom = type_bottom_;
    type_initial = type_init_;
    center_init = x0;
    freq_init = freq;
  }
  
  
  //! computation of the bottom b(x) and velocity c = sqrt(1 - b)
  void GenericWaterWaveEquation
  ::ComputeBottom(double xmin, double xmax, const VectReal_wp& x_,
                  VectReal_wp& b_, VectReal_wp& c_)
  {
    // frequency of bottom
    freq_bottom = alpha;
    int N_ = x_.GetM();
    // DISP(alpha); DISP(p);
    
    // calcul de b et c
    c_.Reallocate(N_); b_.Reallocate(N_);
    //DISP(type_bottom); DISP(PERIODIC_SINUS);
 
    if (type_bottom == PERIODIC_SINUS)
      {
	int p = toInteger(round( freq_bottom*(xmax - xmin)));
        freq_bottom = p/(xmax - xmin);
        //DISP(p); DISP(alpha); DISP(freq_bottom); DISP(beta);
	for (int i = 0; i < N_; i++)
	  {
	    b_(i) = sin(2*pi_wp*freq_bottom*x_(i));
	    c_(i) = sqrt(1 - beta*b_(i));
	  }
      }
    else if (type_bottom == SINUS)
      {
	//int p = toInteger(round( freq_bottom*(xmax - xmin)));
        //freq_bottom = p/(xmax - xmin); DISP(p);
        //DISP(alpha); DISP(freq_bottom); DISP(beta);
	for (int i = 0; i < N_; i++)
	  {
	    b_(i) = sin(2*pi_wp*freq_bottom*x_(i));
	    c_(i) = sqrt(1 - beta*b_(i));
	  }
      }
    else if (type_bottom == HALF_SINUS)
      {
	for (int i = 0; i < N_; i++)
	  {
	    if (x_(i) <= epsilon_machine)
	      b_(i) = 0.0;
	    else
	      b_(i) = pow(sin(2*pi_wp*freq_bottom*x_(i)), 4.0);
	    
	    c_(i) = sqrt(1 - beta*b_(i));
	  }
      }
    else if (type_bottom == LINEAR)
      {
	slope_bottom = 1.0/(xmax - xmin);
	yintercept = 0.5 - slope_bottom*x_(0);
	for (int i = 0; i < N_; i++)
	  {
	    b_(i) = slope_bottom*x_(i) + yintercept;
	    c_(i) = sqrt(1 - beta*b_(i));
	  }
      }
    else
      {
        // default case -> flat bottom
	b_.Fill(0); c_.Fill(1);
      }
    
  }
  
  
  //! computation of initial condition u0(x)
  void GenericWaterWaveEquation::ComputeInitialCondition(double xmin, double xmax,
                                                         const VectReal_wp& x_, VectReal_wp& z0)
  {
    int N_ = x_.GetM();
    z0.Reallocate(N_);
    switch (type_initial)
      {
      case GAUSSIAN :
	for (int i = 0; i < N_; i++)
	  z0(i) = exp(-square(pi_wp*freq_init*(x_(i) - center_init)));
	break;
      case INIT_SINUS :
	for (int i = 0; i < N_; i++)
	  z0(i) = sin(2.0*pi_wp*freq_init*(x_(i) - center_init));
	break;
      case SOLITON :
	for (int i = 0; i < N_; i++)
	  z0(i) = 2.0*c1/square(cosh(sqrt(3.0*c1/(2*nu))*(x_(i) - center_init)));
	break;
      case SOLITON_GN :
        {
          double alpha_s = freq_init;
          double k = sqrt( (3.0*alpha_s*epsilon)/(4.0*mu*(1.0 + alpha_s*epsilon) ) );
          //double c = sqrt(1+alpha_s*epsilon);
          for (int i = 0; i < N_; i++)
            {
              double x = x_(i);
              z0(i) = alpha_s/square(cosh(k*(x - center_init)));
            }
        }
        break;
      case CRENEAU :
	for (int i = 0; i < N_; i++)
	  {
	    if (abs(x_(i)) < 1.0)
	      z0(i) = 0.75;
	    else
	      z0(i) = -0.5;
	  }
        break;
      case BREAKING :
	for (int i = 0; i < N_; i++)
	  z0(i) = pow(abs(x_(i)-center_init), 1.5)
            *exp(-square(pi_wp*freq_init*(x_(i) - center_init)));
	
	break;
      }
  }
  

  //! computation of initial condition u0(x), derivative and second derivative
  void GenericWaterWaveEquation
  ::ComputeInitialCondition(double xmin, double xmax, const VectReal_wp& x_,
                            VectReal_wp& z0, VectReal_wp& z0_x, VectReal_wp& z0_xx)
  {
    int N_ = x_.GetM();
    z0.Reallocate(N_); z0.Fill(0);
    z0_x.Reallocate(N_); z0_x.Fill(0);
    z0_xx.Reallocate(N_); z0_xx.Fill(0);
    Real_wp omega = 2.0*freq_init*pi_wp;
    switch (type_initial)
      {
      case GAUSSIAN :
	{
	  Real_wp arg(0), darg;
	  darg = pi_wp*freq_init;
	  for (int i = 0; i < N_; i++)
	    {
	      arg = darg*(x_(i) - center_init);
	      z0(i) = exp(-arg*arg);
	      z0_x(i) = -2.0*darg*arg*z0(i);
	      z0_xx(i) = (4.0*square(arg*darg) - 2.0*darg*darg)*z0(i);
	    }
	} 
	break;
      case INIT_SINUS :
	{
	  for (int i = 0; i < N_; i++)
	    {
	      z0(i) = sin(omega*(x_(i) - center_init));
	      z0_x(i) = omega*cos(omega*(x_(i) - center_init));
	      z0_xx(i) = -omega*omega*z0(i);
	    }
	  break;
	}
      case SOLITON :
	{
	  Real_wp gamma = sqrt(3.0*c1/(2*nu));	    
	  for (int i = 0; i < N_; i++)
	    {
	      Real_wp sinh_x = sinh(gamma*(x_(i) - center_init));
	      Real_wp cosh_x = cosh(gamma*(x_(i) - center_init));
	      z0(i) = 2.0*c1/(cosh_x*cosh_x);
	      z0_x(i) = -4.0*c1*gamma*sinh_x/pow(cosh_x, 3.0);
	      z0_xx(i) = 4.0*c1*gamma*gamma*(2.0*sinh_x*sinh_x - 1.0)/pow(cosh_x, 4.0);
	    }
	}
	break;
      case SOLITON_GN :
        {
          double alpha_s = freq_init;
          double k = sqrt( (3.0*alpha_s*epsilon)/(4.0*mu*(1.0 + alpha_s*epsilon) ) );
          //double c = sqrt(1+alpha_s*epsilon);
          for (int i = 0; i < N_; i++)
            {
              double x = x_(i);
              double cosh_x = cosh(k*(x - center_init));
              double sinh_x = sinh(k*(x - center_init));
              z0(i) = alpha_s/square(cosh_x);
              z0_x(i) = -2.0*alpha_s*k*sinh_x/pow(cosh_x, 3.0);
              z0_xx(i) = alpha_s*k*k*( -2.0/pow(cosh_x, 2.0) + 6.0*sinh_x*sinh_x/pow(cosh_x, 4.0));
            }
        }
        break;
      case CRENEAU :
	for (int i = 0; i < N_; i++)
	  {
	    if (abs(x_(i)) < 1.0)
	      z0(i) = 0.75;
	    else
	      z0(i) = -0.5;
	  }
	break;
      }
  }
  
  
  //! computation of Gauss-Lobatto basis functions and values on quadrature points
  template<class GenericPb>
  void GenericWaterWaveEquation::InitOrder(int order, const GenericPb& var)
  {
    lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);
    // computation of ValPhi and GradPhi for this order
    int N = var.gauss.GetNbPointsQuad();
    ValPhi.Reallocate(order+1, N);
    GradPhi.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      {
	for (int j = 0; j <= order; j++)
	  GradPhi(i, j) = lob.EvaluatePhiGrad(i, lob.Points(j));
	
	for (int j = 0; j < N; j++)
	  ValPhi(i, j) = lob.EvaluatePhi(i, var.gauss.Points(j));
      }
  }
    
  
  //! computation of gradient u from values of u with Gauss-Lobatto functions
  template<class GenericPb>
  void GenericWaterWaveEquation
  ::GetGradient(const GenericPb& var, const VectReal_wp& Un, VectReal_wp& dUn)
  {
    int nb = 0;
    dUn.Reallocate(var.mesh.GetNbElt()*lob.GetNbPointsQuad());
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	Real_wp val = 0;
	Real_wp invH = 1.0/(var.mesh.Vertex(i+1) - var.mesh.Vertex(i));
	for (int j = 0; j <= lob.GetOrder(); j++)
	  {
	    val = 0;
	    for (int k = 0; k <= lob.GetOrder(); k++)
	      val += GradPhi(k, j)*Un(nb+k);
	    
	    dUn(nb + j) = val*invH;
	  }
	
	nb += lob.GetNbPointsQuad();
      }
  }
  

  //! computation of u on quadrature points from values of u with Gauss-Lobatto functions  
  template<class GenericPb>
  void GenericWaterWaveEquation
  ::GetUquadrature(const GenericPb& var, const VectReal_wp& Un, VectReal_wp& UnQuad)
  {
    int nb = 0;
    UnQuad.Reallocate(var.mesh.GetNbVertices() + var.mesh.GetNbElt()*var.gauss.GetNbPointsQuad());
    
    UnQuad(0) = Un(0);
    int offset = var.mesh.GetNbVertices();
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
	Real_wp val = 0;
	for (int j = 0; j < var.gauss.GetNbPointsQuad(); j++)
	  {
	    val = 0;
	    for (int k = 0; k <= lob.GetOrder(); k++)
	      val += ValPhi(k, j)*Un(nb+k);
	    
	    UnQuad(offset + j) = val;
	  }
	
	UnQuad(i+1) = Un(nb+lob.GetOrder());
	offset += var.gauss.GetNbPointsQuad();
	nb += lob.GetNbPointsQuad();
      }
  }
  
  
  //! default constructor
  InputVariablesWaterWaves::InputVariablesWaterWaves()
  {
    frequency = 0.1;
    center_source = -10.0;
    type_initial_condition = GenericWaterWaveEquation::SOLITON;
    type_bottom = GenericWaterWaveEquation::SINUS;
    type_model = GenericWaterWaveEquation::ORIGINAL;
    type_formulation = FDTD;
    
    x0 = -20.0; xN = 20.0;
    xg_0 = -20.0; xg_N = 20.0;
      
    nb_elt = 200; order = 10;  
    nb_points_fdtd = 1000;
    
    dt = 1.0/(50.0);
    Tfinal = 20;
    
    number_grid_points = 2000;
    delta_snapshot = 10; 
    
    log_eps_min = -1.0;
    log_eps_max = -1.0;
    nb_points_log_eps = 25;
    log_log_curve = false;
    
    time_order = 2;
    type_time_scheme = TimeSchemeEnum::GAUSS_RUNGE_KUTTA;
    
    mu = 0.1;
    
    cte_epsilon = 1.0;
    power_epsilon =  1.0;
    
    cte_alpha = 0.5;
    power_alpha = 1.0;
    
    cte_beta = 0.5;
    power_beta = 0.5;
    
    file_output_history = string("Un.dat");
  }
  
  
  //! setting input data from a line of the data file
  void InputVariablesWaterWaves::
  SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (!keyword.compare("Frequency"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Frequency needs more parameters, for instance :" << endl;
	    cout << "Frequency = f0" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	frequency = to_num<double>(param(0));
      }
    else if (!keyword.compare("SourceCenter"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "SourceCenter needs more parameters, for instance :" << endl;
	    cout << "SourceCenter = x0" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
	
	center_source = to_num<double>(param(0));          
      }
    else if (keyword == "Formulation")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Formulation needs more parameters, for instance :" << endl;
	    cout << "Formulation = Spectral" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	type_formulation = FDTD;
	
	if (param(0) == "Spectral")
	  type_formulation = SPECTRAL;
	else if (param(0) == "Fdtd")
	  type_formulation = FDTD;
	else if (param(0) == "Ldg")
	  type_formulation = LDG;
	else if (param(0) == "Fem")
	  type_formulation = FEM;
      }
    else if (keyword == "FileOutputGrid")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "FileOutputGrid needs more parameters, for instance :" << endl;
	    cout << "FileOutputGrid = file_output" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	file_output_history = param(0);
      }
    else if (!keyword.compare("InitialCondition"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "InitialCondition needs more parameters, for instance :" << endl;
	    cout << "InitialCondition = Soliton" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	if (!param(0).compare("Soliton"))
	  type_initial_condition = GenericWaterWaveEquation::SOLITON;
	if (!param(0).compare("SolitonGN"))
	  type_initial_condition = GenericWaterWaveEquation::SOLITON_GN;
	else if (!param(0).compare("Gaussian"))
	  type_initial_condition = GenericWaterWaveEquation::GAUSSIAN;
	else if (!param(0).compare("Creneau"))
	  type_initial_condition = GenericWaterWaveEquation::CRENEAU;
	else if (!param(0).compare("Sinus"))
	  type_initial_condition = GenericWaterWaveEquation::INIT_SINUS;
	else if (!param(0).compare("Breaking"))
	  type_initial_condition = GenericWaterWaveEquation::BREAKING;
      }
    else if (!keyword.compare("Model"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Model needs more parameters, for instance :" << endl;
	    cout << "Model = Original" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	if (!param(0).compare("Original"))
	  type_model = GenericWaterWaveEquation::ORIGINAL;
	else if (!param(0).compare("Strong"))
	  type_model = GenericWaterWaveEquation::STRONG;
	else if (!param(0).compare("Gentle"))
	  type_model = GenericWaterWaveEquation::GENTLE;
      }
    else if (!keyword.compare("Bottom"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Bottom needs more parameters, for instance :" << endl;
	    cout << "Bottom = Flat" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
	
	if (!param(0).compare("Flat"))
	  type_bottom = GenericWaterWaveEquation::FLAT;
	else if (!param(0).compare("Sinus"))
	  type_bottom = GenericWaterWaveEquation::SINUS;
	else if (!param(0).compare("Linear"))
	  type_bottom = GenericWaterWaveEquation::LINEAR;
	else if (!param(0).compare("HalfSinus"))
	  type_bottom = GenericWaterWaveEquation::HALF_SINUS;
	else if (!param(0).compare("PeriodicSinus"))
	  type_bottom = GenericWaterWaveEquation::PERIODIC_SINUS;
      }
    else if (!keyword.compare("ComputationalDomain"))
      {
	if (param.GetM() <= 2)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "ComputationalDomain needs more parameters, for instance :" << endl;
	    cout << "ComputationalDomain = x0 xN nb_elt" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	x0 = to_num<double>(param(0));
	xN = to_num<double>(param(1));
	nb_elt = to_num<int>(param(2));
	if (param.GetM() > 3)
	  nb_points_fdtd = to_num<int>(param(3));
      }
    else if (!keyword.compare("FindOrder"))
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "FindOrder needs more parameters, for instance :" << endl;
	    cout << "FindOrder = n1 n2" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
	
	n1_order = to_num<int>(param(0));
	n2_order = to_num<int>(param(1));
      }
    else if (!keyword.compare("DisplayGrid"))
      {
	if (param.GetM() <= 2)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "DisplayGrid needs more parameters, for instance :" << endl;
	    cout << "DisplayGrid = x0 xN N" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	xg_0 = to_num<double>(param(0));
	xg_N = to_num<double>(param(1));
	if (param(2) == "AUTO")
	  number_grid_points = (order+1)*nb_elt;
	else
	  number_grid_points = to_num<int>(param(2));
	
	delta_snapshot = to_num<int>(param(3));
      }
    else if (!keyword.compare("Order"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Order needs more parameters, for instance :" << endl;
	    cout << "Order = r" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	order = to_num<int>(param(0));
      }
    else if (!keyword.compare("TimeOrder"))
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "TimeOrder needs more parameters, for instance :" << endl;
	    cout << "TimeOrder = r type_scheme" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	time_order = to_num<int>(param(0));
	if (param.GetM() > 1)
	  {
	    if (param(1) == "LeapFrog")
	      type_time_scheme = TimeSchemeEnum::LEAP_FROG;
	    else if (param(1) == "LeapFrogExact")
	      type_time_scheme = TimeSchemeEnum::LEAP_FROG_EXACT;
	    else if (param(1) == "AdamsBashforth")
	      type_time_scheme = TimeSchemeEnum::ADAMS_BASHFORTH;
	    else if (param(1) == "Talezer")
	      type_time_scheme = TimeSchemeEnum::TALEZER;
	    else if (param(1) == "RungeKuttaLow")
	      type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA_LOW_STORAGE;
	    else if (param(1) == "RungeKutta")
	      type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA;
	    else if (param(1) == "RungeKuttaExact")
	      type_time_scheme = TimeSchemeEnum::RUNGE_KUTTA_EXACT;
	    else if (param(1) == "GaussRungeKutta")
	      type_time_scheme = TimeSchemeEnum::GAUSS_RUNGE_KUTTA;
	    else if (param(1) == "CrankNicolson")
	      type_time_scheme = TimeSchemeEnum::CRANK_NICOLSON;
	    else if (param(1) == "Predictor")
	      type_time_scheme = TimeSchemeEnum::PREDICTOR;
	  }
      }
    else if (!keyword.compare("TimeStep"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "TimeStep needs more parameters, for instance :" << endl;
	    cout << "TimeStep = dt" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	dt = to_num<double>(param(0));
      }
    else if (!keyword.compare("FinalTime"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "FinalTime needs more parameters, for instance :" << endl;
	    cout << "FinalTime = tf" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	Tfinal = to_num<double>(param(0));
      }
    else if (!keyword.compare("ComputeLogLogConvergence"))
      {
	if (param.GetM() <= 2)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "ComputeLogLogConvergence needs more parameters, for instance :" << endl;
	    cout << "ComputeLogLogConvergence = log_min log_max N" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }
	
	log_log_curve = true;
	log_eps_min = to_num<double>(param(0));
	log_eps_max = to_num<double>(param(1));
	nb_points_log_eps = to_num<int>(param(2));
      }
    else if (!keyword.compare("Mu"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Mu needs more parameters, for instance :" << endl;
	    cout << "Mu = mu" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	mu = to_num<double>(param(0));          
      }
    else if (!keyword.compare("Epsilon"))
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Epsilon needs more parameters, for instance :" << endl;
	    cout << "Epsilon = cte exponent" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	cte_epsilon = to_num<double>(param(0));
	power_epsilon = to_num<double>(param(1));
      }
    else if (!keyword.compare("Alpha"))
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Alpha needs more parameters, for instance :" << endl;
	    cout << "Alpha = cte exponent" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	cte_alpha = to_num<double>(param(0));
	power_alpha = to_num<double>(param(1));
      }
    else if (!keyword.compare("Beta"))
      {
	if (param.GetM() <= 1)
	  {
	    cout << "In SetInputData of InputVariablesWaterWaves" << endl;
	    cout << "Beta needs more parameters, for instance :" << endl;
	    cout << "Beta = cte exponent" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	cte_beta = to_num<double>(param(0));
	power_beta = to_num<double>(param(1));
      }
  }
  
}

#define MONTJOIE_FILE_GENERIC_WATER_WAVE_EQUATION_CXX
#endif
