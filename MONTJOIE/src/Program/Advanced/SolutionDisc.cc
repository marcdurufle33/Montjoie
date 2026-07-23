#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

//! class for solving Helmholtz equation in circular layers
/*!
  On each circular layer, the solution solves Helmholtz equation
  - \omega^2 \rho u - div(\mu grad u) = 0
  
  the parameters rho and mu are constant in each layer
  Transmission conditions between layers are given as :
  u continuous across the interface
  mu du/dn continuous across the interface
  
  The first layer can be a dielectric disc (i.e. a disc with
  constant rho and mu) or dielectric circular ring, and a Dirichlet,
  Neumann or impedance condition condition is set on the intern boundary r = r0
  
 */
class AnalyticalHelmholtzDiscSolver : public InputDataProblem_Base
{
public :
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;
  
  //! parameters for output on a predefined grid
  GridInterpolationFull<Dimension2> output_grid;
  //! type of output
  int type_output;
  //! available types of output
  enum {OUTPUT_GRID_ASCII, OUTPUT_GRID_BINARY, OUTPUT_MESH_ASCII, OUTPUT_MESH_BINARY};
  
  //! if true, we compute the total field, otherwise we compute the diffracted field
  bool total_field;
  
  //! parameters rho and mu in Helmholtz equations  
  Vector<Complex_wp> rho, mu;
  VectReal_wp sigma;
  
  //! if true an impedance condition beta u + du/dn = 0 is set, otherwise u = 0 is set
  bool impedance_cond;
  Complex_wp beta_impedance;
  
  //! if true, we replace sommerfeld condition by a first order condition placed on an external circle
  bool first_order_abc, curved_abc, grazing_abc;
  int order_abc;
  Real_wp theta_cla_coeff, gamma_cla_coeff, zeta_cla_coeff;
  // parameters for Pade approximation of sqrt
  Real_wp alpha_pade, epsilon_pade; int nb_coef_pade;

  
  //! which data must be computed ?
  int type_data;
  //! available types of data
  enum {TRACE_U, TRACE_DU_DN, FAR_FIELD};
  
  //! number of modes (number of Bessel functions) automatically computed ?
  int type_input_modes;
  enum {MODES_AUTO, MODES_PYTHON, MODES_INPUT, MODES_MANUAL};
  //! number of Bessel functions 
  int n_modes;
  
  //! radius of all the circular layers
  Vector<Real_wp> ri_circle;
  Real_wp radius_output;
  //! radius with first-order absorbing condition
  Real_wp radius_circle_outside;
  //! pulsation
  Real_wp omega, frequency;
  //! wave number k in each media
  VectComplex_wp kwave;
  
  //! name of the mesh file
  string file_mesh;
  //! names of output file
  string name_file;
  
  //! file where the radar cross section is stored
  string file_RCS;
  int nbPointsRCS; //!< number of angles for RCS
  Real_wp first_angle_RCS, last_angle_RCS; //!< first and last angle

  // parameters to recover the solution in time domain
  int nb_points_time; Real_wp Tmax, c0, offset_arg, freq_pulse;
  int print_level;
  
  // if true, the harmonic fields are written 
  bool write_harmonic_sol;
  
public:
  AnalyticalHelmholtzDiscSolver();
  
  void SetInputData(const string& description_field, const VectString& parameters);

  void ComputeSolution(const Vector<VectComplex_wp>& alpha,
                       const Vector<VectComplex_wp>& beta,
                       const VectR2& Points, VectComplex_wp& Usol);
  
  void WriteSolution(const Vector<VectComplex_wp>& alpha,
                     const Vector<VectComplex_wp>& beta);
  
  void WriteTimeSolution();
  
  void ComputeDecompositionCoef(Vector<VectComplex_wp>& alpha,
                                Vector<VectComplex_wp>& beta);
  
};


//! default constructor
AnalyticalHelmholtzDiscSolver::AnalyticalHelmholtzDiscSolver()
{
  type_output = OUTPUT_GRID_BINARY;
  impedance_cond = false;
  beta_impedance = 0.0;
  total_field = false;
  type_data = TRACE_U;
  radius_circle_outside = 2.0;
  omega = 2.0*pi_wp; frequency = 1.0;
  
  Complex_wp one; SetComplexOne(one);
  rho.Reallocate(1); rho.Fill(one);
  mu.Reallocate(1); mu.Fill(one);
  
  file_mesh = string("disc.mesh");
  name_file = string("U.dat");
  
  file_RCS = string("Rcs.dat");
  first_angle_RCS = 0.0; last_angle_RCS = 2.0*pi_wp; nbPointsRCS = 361;
  
  first_order_abc = false; radius_circle_outside = 2.0;
  type_input_modes = MODES_AUTO; n_modes = 10;
  curved_abc = false; order_abc = 1; grazing_abc = false;
  theta_cla_coeff = zeta_cla_coeff = 0.0; gamma_cla_coeff= 1.0/4.0;
  nb_coef_pade = 0; alpha_pade = 0; epsilon_pade = 0;

  nb_points_time = 0; freq_pulse = 1.0;
  Tmax = 0.0; c0 = 1.0; offset_arg = 0.0;
  print_level = 1;
  
  write_harmonic_sol = false; radius_output = -1.0;  
}


//! reads a line of the data file
void AnalyticalHelmholtzDiscSolver::SetInputData(const string& description_field, const VectString& parameters)
{
  if (!description_field.compare("TypeOutput"))
    {
      if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "TypeOutput needs more parameters, for instance :" << endl;
          cout << "TypeOutput = Matlab BINARY" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
      if (!parameters(0).compare("Matlab"))
	{
	  type_output = OUTPUT_GRID_BINARY;
	  if (!parameters(1).compare("ASCII"))
	    type_output = OUTPUT_GRID_ASCII;
	}
      else if (!parameters(0).compare("Medit"))
	{
	  type_output = OUTPUT_MESH_BINARY;
	  if (!parameters(1).compare("ASCII"))
	    type_output = OUTPUT_MESH_ASCII;
	}
      
      if (parameters.GetM() >= 3)
        {
          if (parameters(2) == "YES")
            write_harmonic_sol = true;
          else
            write_harmonic_sol = false;
        }
    }
  else if (description_field.substr(0, 5) == "Sismo")
    {
      // parameters of the grid
      output_grid.SetInputData(description_field, parameters);
    }
  else if (!description_field.compare("TotalField"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "TotalField needs more parameters, for instance :" << endl;
          cout << "TotalField = YES" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      if (!parameters(0).compare("YES"))
	total_field = true;
      else
	total_field = false;
    }
  else if (!description_field.compare("TypeCondition"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "TypeCondition needs more parameters, for instance :" << endl;
          cout << "TypeCondition = DIRICHLET" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      if (!parameters(0).compare("DIRICHLET"))
	{
          impedance_cond = false;
          beta_impedance = Complex_wp(0, 0);
        }
      else if (!parameters(0).compare("NEUMANN"))
	{
          impedance_cond = true;
          beta_impedance = Complex_wp(0, 0);
        }
      else if (parameters(0) == "IMPEDANCE")
        {
          if (parameters.GetM() <= 1)
            {
              cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
              cout << "TypeCondition needs more parameters, for instance :" << endl;
              cout << "TypeCondition = IMPEDANCE beta" << endl;
              cout << "Current parameters are : " << endl << parameters << endl;
              abort();
            }
          
          impedance_cond = true;
          beta_impedance = to_num<Complex_wp>(parameters(1));
        }
    }
  else if (!description_field.compare("Rho"))
    {
      rho.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	to_num(parameters(i), rho(i));
    }
  else if (!description_field.compare("Mu"))
    {
      mu.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	to_num(parameters(i), mu(i));
    }
  else if (!description_field.compare("Sigma"))
    {
      sigma.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	to_num(parameters(i), sigma(i));
    }
  else if (!description_field.compare("FirstOrder_ABC"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "FirstOrder_ABC needs more parameters, for instance :" << endl;
          cout << "FirstOrder_ABC = YES" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      if (!parameters(0).compare("YES"))
	first_order_abc = true;
      else if (!parameters(0).compare("CURVED"))
        {
          first_order_abc = true;
          curved_abc = true;
        }
      else if (parameters(0) == "PADE")
        {
          if (parameters.GetM() <= 3)
            {
              cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
              cout << "FirstOrder_ABC needs more parameters, for instance :" << endl;
              cout << "FirstOrder_ABC = PADE alpha n epsilon" << endl;
              cout << "Current parameters are : " << endl << parameters << endl;
              abort();
            }
          
          alpha_pade = to_num<Real_wp>(parameters(1))*pi_wp/180;
          nb_coef_pade = to_num<int>(parameters(2));
          epsilon_pade = to_num<Real_wp>(parameters(3));
          first_order_abc = true;
          return;
        }
      else
	first_order_abc = false;
      
      if (parameters.GetM() > 1)
        {
          order_abc = to_num<int>(parameters(1));
          DISP(order_abc);
        }
      
      if (parameters.GetM() > 2)
        {
          if (parameters(2) == "GRAZING")
            grazing_abc = true;
          else
            grazing_abc = false;
        }
    }
  else if(!description_field.compare("Parameters_ABC"))
    {
      if (parameters.GetM() > 2)
        {
          gamma_cla_coeff = to_num<Real_wp>(parameters(0));
          theta_cla_coeff = to_num<Real_wp>(parameters(1));
          zeta_cla_coeff  = to_num<Real_wp>(parameters(2)); //DISP(zeta_cla_coeff);
        }
      else
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "Parameters_ABC needs three parameters, for instance :" << endl;
          cout << "Parameters_ABC = 0 0 0" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
    }
  else if (!description_field.compare("TypeData"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "TypeData needs more parameters, for instance :" << endl;
          cout << "TypeData = NEAR_FIELD" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
      if (!parameters(0).compare("NEAR_FIELD"))
	type_data = TRACE_U;
      else if (!parameters(0).compare("SURFACIC_FIELD"))
	type_data = TRACE_DU_DN;
      else if (!parameters(0).compare("FAR_FIELD"))
	type_data = FAR_FIELD;
    }
  else if (!description_field.compare("NumberModes"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "NumberModes needs more parameters, for instance :" << endl;
          cout << "NumberModes = AUTO" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      if (!parameters(0).compare("AUTO"))
	{
	  type_input_modes = MODES_AUTO;
	}
      else if (!parameters(0).compare("INPUT"))
	{
	  type_input_modes = MODES_INPUT;
	}
      else if (!parameters(0).compare("PYTHON"))
	{
	  type_input_modes = MODES_PYTHON;
	}
      else
	{
	  type_input_modes = MODES_MANUAL;
	  to_num(parameters(0), n_modes);
	}
    }
   else if (!description_field.compare("RadiusCircle"))
    {
      if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "RadiusCircle needs more parameters, for instance :" << endl;
          cout << "RadiusCircle = r0 r1" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      ri_circle.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	to_num(parameters(i), ri_circle(i));
      
      radius_circle_outside = ri_circle(parameters.GetM()-1);
    }
   else if (!description_field.compare("RadiusOutput"))
   {
       if (parameters.GetM() <= 0)
       {
           cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
           cout << "RadiusOutput needs more parameters, for instance :" << endl;
           cout << "RadiusOutput = r0 " << endl;
           cout << "Current parameters are : " << endl << parameters << endl;
           abort();
       }
       radius_output = to_num<Real_wp>(parameters(0));
   }
  else if (!description_field.compare("Frequency"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "Frequency needs more parameters, for instance :" << endl;
          cout << "Frequency = freq" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      to_num(parameters(0), omega);
      frequency = omega; omega *= 2.0*pi_wp;
    }
  else if (!description_field.compare("MeshFile"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "MeshFile needs more parameters, for instance :" << endl;
          cout << "MeshFile = test.mesh" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      file_mesh = parameters(0);
    }
  else if (!description_field.compare("OutputFile"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "OutputFile needs more parameters, for instance :" << endl;
          cout << "OutputFile = file_name" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      name_file = parameters(0);
    }
  else if (!description_field.compare("AngleRCS"))
    {
      if (parameters.GetM() <= 2)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "AngleRCS needs more parameters, for instance :" << endl;
          cout << "AngleRCS = teta0 teta1 N" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      to_num(parameters(0), first_angle_RCS); first_angle_RCS = first_angle_RCS*pi_wp/180.0-pi_wp;
      to_num(parameters(1), last_angle_RCS); last_angle_RCS = last_angle_RCS*pi_wp/180.0-pi_wp;
      to_num(parameters(2), nbPointsRCS);
    }
  else if (!description_field.compare("FileRCS"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "FileRCS needs more parameters, for instance :" << endl;
          cout << "FileRCS = file_name" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      file_RCS  = string(parameters(0));
    }
  else if (!description_field.compare("TimeImpulsion"))
    {
      if (parameters.GetM() <= 3)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "TimeImpulsion needs more parameters, for instance :" << endl;
          cout << "TimeImpulsion = N Gaussian freq Tmax" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      nb_points_time = to_num<int>(parameters(0));
      if (parameters(1) == "Gaussian")
        {
          freq_pulse = to_num<Real_wp>(parameters(2));
          Tmax = to_num<Real_wp>(parameters(3));
        }
    }
  else if (!description_field.compare("WaveVelocity"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "WaveVelocity needs more parameters, for instance :" << endl;
          cout << "WaveVelocity = file_name" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      c0 = to_num<Real_wp>(parameters(0));
    }
  else if (!description_field.compare("OffsetTimePulse"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHelmholtzDiscSolver" << endl;
          cout << "OffsetTimePulse needs more parameters, for instance :" << endl;
          cout << "OffsetTimePulse = file_name" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      offset_arg = to_num<Real_wp>(parameters(0));
    }
}


//! computes the analytical solution given by coefficients alpha and beta on the points
void AnalyticalHelmholtzDiscSolver::ComputeSolution(const Vector<VectComplex_wp>& alpha,
                                                    const Vector<VectComplex_wp>& beta,
                                                    const VectR2& Points, VectComplex_wp& Usol)
{
  int N_max = n_modes;
  Real_wp begin_order = 0.0;

  Vector<Complex_wp> Jn(N_max), Hn(N_max), Hn2(N_max), Hn2_prime(N_max);
  Vector<Complex_wp> Jn_prime(N_max), Hn_prime(N_max);
  
  Vector<Complex_wp> coefIn(N_max);
  Complex_wp i_pow_n(1.0, 0.0);
  for (int i = 0; i < N_max; i++)
    {
      if (i == 0)
        coefIn(i) = Complex_wp(1, 0);
      else
        coefIn(i) = 2.0*i_pow_n;
      
      i_pow_n *= Iwp;
    }
  
  Usol.Reallocate(Points.GetM());
  Complex_wp u_loc, kr, expIteta, expIn_teta;
  Real_wp x, y, teta, r, cos_teta, sin_teta;
  for (int i = 0; i < Points.GetM(); i++)
    {
      // cartesian coordinates
      x = Points(i)(0); 
      y = Points(i)(1); 
      
      // conversion to polar coordinates
      CartesianToPolar(x, y, r, teta);
      cos_teta = cos(teta); sin_teta = sin(teta);
      expIteta = Complex_wp(cos_teta, sin_teta);
      
      // we find the layer where the point (x, y) is
      int num = ri_circle.GetM()-1;
      for (int n = 0; n < ri_circle.GetM(); n++)
        {
          Real_wp coef_debordement = 1.0-50*epsilon_machine;
          if (n == ri_circle.GetM() - 1)
            coef_debordement = 1.0+50*epsilon_machine;
          
          if (r < coef_debordement*ri_circle(n))
            {
              num = n-1;
              break;
            }
        }
      
      SetComplexZero(u_loc);
      if ((num >= 0) && (num < ri_circle.GetM()-1))
        {
          kr = kwave(num)*r;

	  // Bessel functions and derivatives are computed
        ComputeDeriveBesselAndHankel(begin_order, N_max, kr, Jn, Hn, Jn_prime, Hn_prime);
          
          // Hankel functions of the second kind
          for (int n = 0; n < N_max; n++)
            Hn2(n) = -Hn(n) + 2.0*Jn(n);
          
          if (type_data == TRACE_DU_DN)
            for (int n = 0; n < N_max; n++)
              {
                Hn2_prime(n) = kwave(num)*(-Hn_prime(n) + 2.0*Jn_prime(n));
                Jn_prime(n) *= kwave(num);
                Hn_prime(n) *= kwave(num);
              }
          
          // computing the total field
          SetComplexOne(expIn_teta);
          if (type_data == TRACE_U)
            {
              if (num == ri_circle.GetM()-2)
                for (int n = 0; n < N_max; n++)
                  {
                    u_loc += (Jn(n) + alpha(n)(num)*Hn(n)+beta(n)(num)*Hn2(n))*coefIn(n)*real(expIn_teta);
                    expIn_teta *= expIteta;
                  }
              else
                for (int n = 0; n < N_max; n++)
                  {
                    u_loc += (alpha(n)(num)*Hn(n)+beta(n)(num)*Hn2(n))*coefIn(n)*real(expIn_teta);
                    expIn_teta *= expIteta;
                  }                
            }
          else if (type_data == TRACE_DU_DN)
            {
              if (num == ri_circle.GetM()-2)
                for (int n = 0; n < N_max; n++)
                  {
                    u_loc += (alpha(n)(num) * Hn_prime(n) + beta(n)(num)*Hn2_prime(n))*coefIn(n)*real(expIn_teta);
                    expIn_teta *= expIteta;
                  }
              else
                for (int n = 0; n < N_max; n++)
                  {
                    u_loc += (Jn_prime(n) + alpha(n)(num) * Hn_prime(n) + beta(n)(num)*Hn2_prime(n))*coefIn(n)*real(expIn_teta);
                    expIn_teta *= expIteta;
                  }
            }          
          
          // then subtracting the plane wave for the diffracted field
          if (!total_field)
            u_loc -= exp(Iwp*kwave(ri_circle.GetM()-2)*x);
        }
      
      Usol(i) = u_loc;
    }
}


//! computes the solution given by coefficients alpha and beta, and writes it on the asked output
void AnalyticalHelmholtzDiscSolver
::WriteSolution(const Vector<VectComplex_wp>& alpha,
                const Vector<VectComplex_wp>& beta)
{
  VectComplex_wp Usol;
  VectR2 Points; VectReal_wp TetaPoints;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_GRID_BINARY))
    output_grid.GenerateGridPoints(Points, TetaPoints);
  else
    {
      Mesh<Dimension2> mesh;
      mesh.Read(file_mesh);
      Points = mesh.Vertex();
    }
  
  ComputeSolution(alpha, beta, Points, Usol);
  
  int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
  bool ascii = false;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_MESH_ASCII))
    ascii = true;
  
  Dimension2 dim;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_GRID_BINARY))
    WriteMatlab(Usol, output_grid, name_file, double_prec, ascii);
  else
    WriteMedit(Usol, name_file, dim, double_prec);
}


//! computes and writes the solution in time domain
void AnalyticalHelmholtzDiscSolver::WriteTimeSolution()
{
  // times t
  VectReal_wp time_interval(nb_points_time);
  Real_wp dt = 2.0*Tmax/Real_wp(nb_points_time);
  time_interval.Fill(0);
  for (int i = 0; i < nb_points_time/2; i++)
    {
      time_interval(i) = -Real_wp(i)*dt;
      time_interval(i+nb_points_time/2) = -Real_wp(i-nb_points_time/2)*dt;
    }
  
  time_interval.WriteText("time_interval.dat");
  
  // frequencies omega
  VectReal_wp Omega(nb_points_time/2+1);
  Real_wp domega = 2.0*pi_wp/(dt*nb_points_time);
  for (int i = 0; i <= nb_points_time/2; i++)
    Omega(i) = i*domega;
  
  Omega.WriteText("pulsation.dat");
  
  // evaluation of pulse in time
  VectReal_wp pulse_time(nb_points_time);
  freq_pulse /= c0;
  // offset_arg must be involved ?
  for (int i = 0; i < nb_points_time; i++)
    pulse_time(i) = exp(-square(pi_wp*(freq_pulse*(c0*time_interval(i) + offset_arg) - 1.0)));
  
  pulse_time.WriteText("pulse_time.dat");
  
  // computation of Fourier transform
  FftRealInterface fft_time;
  fft_time.Init(nb_points_time);
  
  VectComplex_wp pulse_freq(nb_points_time/2+1);
  fft_time.ApplyForward(pulse_time, pulse_freq);
  
  pulse_freq.Write("pulse_harmonic.dat");
  
  // preparing outputs
  int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
  bool ascii = false;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_MESH_ASCII))
    ascii = true;
  
  Dimension2 dim;
  
  // points where the solution must be computed
  VectR2 Points; VectReal_wp TetaPoints;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_GRID_BINARY))
    output_grid.GenerateGridPoints(Points, TetaPoints);
  else
    {
      Mesh<Dimension2> mesh;
      mesh.Read(file_mesh);
      Points = mesh.Vertex();
    }

  // computation of the solution for each frequency in Omega (except the first one)
  Vector<VectComplex_wp> alpha, beta; VectComplex_wp Uhat;
  Vector<VectComplex_wp> Usol(Points.GetM());
  for (int i = 0; i < Usol.GetM(); i++)
    {
      Usol(i).Reallocate(Omega.GetM());
      Usol(i).Fill(0);
    }
  
  Real_wp threshold = 10.0*epsilon_machine;
  Real_wp Umax = pulse_freq.GetNormInf(); DISP(Umax);
  for (int i = 1; i < Omega.GetM(); i++)
    {
      omega = Omega(i)/c0;
      frequency = omega / (2.0*pi_wp);
            
      print_level = 0;
      if (abs(pulse_freq(i)) >= threshold*Umax)
        {
          if (i%20 == 0)
            cout << "Solving frequency " << i << " : " << omega << endl;
          
          ComputeDecompositionCoef(alpha, beta);      
          ComputeSolution(alpha, beta, Points, Uhat);
          
          if (write_harmonic_sol)
            {
              string nom = GetBaseString(name_file) + NumberToString(i) + "_H.dat";
              WriteMatlab(Uhat, output_grid, nom, double_prec, ascii);
            }
          
          //Complex_wp coef = exp(-Iwp*omega*offset_arg)/c0*pulse_freq(i);
          Complex_wp coef = pulse_freq(i);
          for (int j = 0; j < Points.GetM(); j++)
            Usol(j)(i) = conj(coef*Uhat(j));
        }
    }
  
  Uhat.Clear();
  VectReal_wp Utime(nb_points_time);
  
  // computing inverse of Fourier transform to get the solution in time
  Vector<VectReal_wp> Usol_time(nb_points_time);
  for (int i = 0; i < nb_points_time; i++)
    Usol_time(i).Reallocate(Points.GetM());
  
  for (int i = 0; i < Points.GetM(); i++)
    {
      fft_time.ApplyInverse(Usol(i), Utime);
      for (int j = 0; j < nb_points_time; j++)
        Usol_time(j)(i) = Utime(j);
    }
  
  // then writing the solution
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_GRID_BINARY))
    {
      for (int i = 0; i < nb_points_time/2; i++)
        {
          string nom = GetBaseString(name_file) + NumberToString(i) + ".dat";
          WriteMatlab(Usol_time(i), output_grid, nom, double_prec, ascii);
        }
    }
  else
    {
      for (int i = 0; i < nb_points_time/2; i++)
        {
          string nom = GetBaseString(name_file) + to_str(i) + ".dat";
          WriteMedit(Usol_time(i), nom, dim, double_prec);
        }
    }
}


bool StartsWith(const string& s, const string& chaine)
{
  if (s.substr(0, chaine.size()) == chaine)
    return true;
  
  return false;
}

void WriteMatrixData(ostream& out, VectReal_wp& r)
{
  out << "matrix([";
  
  for (int i = 0; i < r.GetM(); i++)
    {
      out << "mpf('" << r(i) << "')";
      if (i < r.GetM()-1)
	out << ", ";
    }
  
  out << "])";
}

void WriteMatrixData(ostream& out, VectComplex_wp& r)
{
  out << "matrix([";
  
  for (int i = 0; i < r.GetM(); i++)
    {
      out << "mpf('" << real(r(i)) << "')+1j*mpf('" << imag(r(i)) << "')";
      if (i < r.GetM()-1)
	out << ", ";
    }
  
  out << "])";
}


//! computes coefficients alpha and beta of the solution for a plane wave
void AnalyticalHelmholtzDiscSolver
::ComputeDecompositionCoef(Vector<VectComplex_wp>& alpha,
                           Vector<VectComplex_wp>& beta)
{
  if ((rho.GetM() != mu.GetM()) || (rho.GetM() != ri_circle.GetM()-1))
    {
      cout << "There must be N-1 physical indices rho and mu where N is the number"
           << " of radii" << endl;
      
      cout << "Current radii are " << ri_circle << endl;
      cout << "Current rho are " << rho << endl;
      cout << "Current mu are " << mu << endl;
      abort();
    }
  
  // number of different layers
  int N = ri_circle.GetM()-1;
  
  // we compute the wave number for each layer
  kwave.Reallocate(N);
  for (int i = 0; i < rho.GetM(); i++)
    {
      Complex_wp rho_tilde = rho(i);
      if (i < sigma.GetM())
        rho_tilde += Iwp*sigma(i)/omega;
      
      kwave(i) = sqrt(rho_tilde/mu(i))*omega;
    }
  
  // maximum value of kr for incident wave
  Real_wp ka = abs(kwave(N-1))*ri_circle(N);
    if(radius_output>0.0)
     ka = abs(kwave(N-1))*radius_output;
        
  
  if (type_input_modes == MODES_PYTHON)
    {
      ofstream file_out("temp.py");
      ifstream file_in("MATLAB/helmholtz_circle.py");
      file_out.precision(15);
      
      if (!file_in.is_open())
	{
	  cout << "Unable to open file MATLAB/helmholtz_circle.py" << endl;
	  abort();
	}
      
      string line;
      while (getline(file_in, line))
	{
	  if (StartsWith(line, "rho = "))
	    {
	      file_out << "rho = ";
	      WriteMatrixData(file_out, rho);
	    }
	  else if (StartsWith(line, "mu = "))
	    {
	      file_out << "mu = ";
	      WriteMatrixData(file_out, mu);
	    }
	  else if (StartsWith(line, "sigma = "))
	    {
	      if (sigma.GetM() <= 0)
		{
		  sigma.Reallocate(rho.GetM());
		  sigma.Fill(0);
		}

	      file_out << "sigma = ";
	      WriteMatrixData(file_out, sigma);
	    }
	  else if (StartsWith(line, "radius = "))
	    {
	      file_out << "radius = ";
	      WriteMatrixData(file_out, ri_circle);
	    }
	  else if (StartsWith(line, "freq = "))
	    {
	      file_out << "freq = mpf('";
	      file_out << frequency << "')";	      
	    }
	  else if (StartsWith(line, "first_order_abc = "))
	    {
	      if (first_order_abc)
		file_out << "first_order_abc = True";
	      else
		file_out << "first_order_abc = False";
	    }
	  else if (StartsWith(line, "impedance_condition = "))
	    {
	      if (impedance_cond)
		file_out << "impedance_condition = True";
	      else
		file_out << "impedance_condition = False";
	    }
	  else if (StartsWith(line, "beta_impedance = "))
	    {
	      file_out << "beta_impedance = mpf('";
	      file_out << real(beta_impedance) << "') + 1j*mpf('" << imag(beta_impedance) << "')";	      
	    }	  
	  else
	    file_out << line;
	  
	  file_out << '\n';
	}
      
      file_out.close();
      
      // on execute le script python
      system("python temp.py");
    }

  if ((type_input_modes == MODES_INPUT) || (type_input_modes == MODES_PYTHON))
    {
      Matrix<Complex_wp> A, B;
      A.ReadText("alpha.dat"); B.ReadText("beta.dat");
      
      int N_max = A.GetM();
      n_modes = N_max;
      alpha.Reallocate(N_max); beta.Reallocate(N_max);
      for (int n = 0; n < N_max; n++)
	{
	  alpha(n).Reallocate(A.GetN());
	  beta(n).Reallocate(A.GetN());
	  for (int i = 0; i < A.GetN(); i++)
	    {
	      alpha(n)(i) = A(n, i);
	      beta(n)(i) = B(n, i);
	    }
	}
      
      //DISP(alpha); DISP(beta);
      return;
    }

  // number of Bessel functions to take into account
  int N_max = 0;
  if (type_input_modes == MODES_AUTO)
    {
      N_max = ComputeOrder(ka, epsilon_machine);
      if (N_max <= 0)
        N_max = 2;
      
      n_modes = N_max;
    }
  else
    N_max = n_modes;
  
  cout << "Number of modes = " << N_max << endl;
  Real_wp begin_order = 0.0;
  Complex_wp zero, kr, one; SetComplexZero(zero);
  SetComplexOne(one);
  
  // computes the spherical bessel functions for each layer (and at the two radii of the layer)
  Vector<VectComplex_wp> Jn(2*N), Hn(2*N), Hn2(2*N), Jn_prime(2*N), Hn_prime(2*N), Hn2_prime(2*N);
  //Vector<VectComplex_wp> Yn(2*N), Yn_prime(2*N);
  for (int i = 0; i < N; i++)
    {
      kr = kwave(i)*ri_circle(i);
      
      // Spherical Bessel functions and derivatives are computed
      ComputeDeriveBesselAndHankel(begin_order, N_max, kr,
                                  Jn(2*i), Hn(2*i), Jn_prime(2*i), Hn_prime(2*i));
      
         // Hankel functions of the second kind
      Hn2(2*i).Reallocate(N_max);
      Hn2_prime(2*i).Reallocate(N_max);
      for (int n = 0; n < N_max; n++)
        {
          Hn2(2*i)(n) = -Hn(2*i)(n) + 2.0*Jn(2*i)(n);
          Hn2_prime(2*i)(n) = kwave(i)*(-Hn_prime(2*i)(n) + 2.0*Jn_prime(2*i)(n));
          Jn_prime(2*i)(n) *= kwave(i);
          Hn_prime(2*i)(n) *= kwave(i);
        }

      kr = kwave(i)*ri_circle(i+1);
      
      // Spherical Bessel functions and derivatives are computed
      ComputeDeriveBesselAndHankel(begin_order, N_max, kr,
                                  Jn(2*i+1), Hn(2*i+1), Jn_prime(2*i+1), Hn_prime(2*i+1));
      
      // Spherical Hankel functions of the second kind
      Hn2(2*i+1).Reallocate(N_max);
      Hn2_prime(2*i+1).Reallocate(N_max);
      for (int n = 0; n < N_max; n++)
        {
          Hn2(2*i+1)(n) = -Hn(2*i+1)(n) + 2.0*Jn(2*i+1)(n);
          Hn2_prime(2*i+1)(n) = kwave(i)*(-Hn_prime(2*i+1)(n) + 2.0*Jn_prime(2*i+1)(n));
          Jn_prime(2*i+1)(n) *= kwave(i);
          Hn_prime(2*i+1)(n) *= kwave(i);
        }
    }

  // expressing the linear system
  //TinyBandMatrix<Complex_wp, 2> lin_sys;
  Matrix<Complex_wp> lin_sys;       IVect pivot;
  lin_sys.Reallocate(2*N, 2*N);
  Vector<Complex_wp> rhs(2*N);
  
  // each mode is decoupled
  alpha.Reallocate(N_max); beta.Reallocate(N_max);
  for (int n = 0; n < N_max; n++)
    {
      lin_sys.Fill(zero);
      
      int num = 0;
      // first equation of the linear system set by the first radius
      // first radius -> if equal to 0, we have a dielectric internal disc
      if (ri_circle(0) == Real_wp(0))
        {
          // in this case, we have alpha_n - beta_n = 0
          lin_sys(num, num) = one; lin_sys(num, num+1) = -one;
          rhs(num) = zero;
        }
      else
        {
          if (impedance_cond)
            {
              // Robin condition -du/dn + beta u = 0
              lin_sys.Get(num, num) = -beta_impedance*Hn(0)(n) + Hn_prime(0)(n);
              lin_sys.Get(num, num+1) = -beta_impedance*Hn2(0)(n) + Hn2_prime(0)(n);
              if (N == 1)
                rhs(num) = -(-beta_impedance*Jn(0)(n) + Jn_prime(0)(n));
              else
                rhs(num) = zero;
            }
          else
            {
              // Dirichlet condition u = 0
              lin_sys.Get(num, num) = Hn(0)(n);
              lin_sys.Get(num, num+1) = Hn2(0)(n);
              if (N == 1)
                rhs(num) = -Jn(0)(n);
              else
                rhs(num) = zero;
            }
        }
      
      num++;
      
      // loop over transmission conditions (two equations by interface)
      for (int i = 1; i < ri_circle.GetM()-1; i++)
        {
          // continuity of the solution
          lin_sys.Get(num, 2*i-2) = Hn(2*i-1)(n);
          lin_sys.Get(num, 2*i-1) = Hn2(2*i-1)(n);
          lin_sys.Get(num, 2*i) = -Hn(2*i)(n);
          lin_sys.Get(num, 2*i+1) = -Hn2(2*i)(n);
          if (i == N-1)
            rhs(num) = Jn(2*i)(n);
          else
            rhs(num) = zero;
          
          num++;
          
          // continuity of mu du_dn
          lin_sys.Get(num, 2*i-2) = Hn_prime(2*i-1)(n)*mu(i-1);
          lin_sys.Get(num, 2*i-1) = Hn2_prime(2*i-1)(n)*mu(i-1);
          lin_sys.Get(num, 2*i) = -Hn_prime(2*i)(n)*mu(i);
          lin_sys.Get(num, 2*i+1) = -Hn2_prime(2*i)(n)*mu(i);
          if (i == N-1)
            rhs(num) = Jn_prime(2*i)(n)*mu(i);
          else
            rhs(num) = zero;
          
          num++;
        }
      
      // last condition due to Sommerfeld condition (at finite distance or not)
      if (first_order_abc)
        {
          Complex_wp beta_abc(0);
          
          Real_wp Kappa = 1.0/ri_circle(N);
          Real_wp n2 = square(n/ri_circle(N)); Complex_wp k = kwave(N-1);
          Complex_wp ik = -Iwp*k, alpha_cla(1, 0), beta_cla(ik), delta_cla(0, 0), eta_cla(0, 0);
          
          Real_wp theta_cla = theta_cla_coeff*Kappa;
          Real_wp gamma_cla = gamma_cla_coeff*Kappa;
          Real_wp zeta_cla  =  zeta_cla_coeff*Kappa;

     
          if (order_abc == 1)//Km1 delta
            {
              if (curved_abc)
                beta_cla +=  gamma_cla + Kappa/4.0;
              
              alpha_cla += (gamma_cla-Kappa/4.0)/(ik);
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if (order_abc == 2) // Km1 omega
            {
              beta_cla += gamma_cla + Kappa/4.0;
              
              alpha_cla += (gamma_cla-Kappa/4.0)/(ik);
              
              delta_cla = 1.0/(2.0*ik);
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if (order_abc == 3)
            {// Km2 delta 2
              
              alpha_cla =-(1+(theta_cla+gamma_cla-Kappa/4)/(ik)+(pow(Kappa,2)/8+theta_cla*gamma_cla-(theta_cla+gamma_cla)*Kappa/4)/(ik*ik)-Kappa/4*(Kappa/4+zeta_cla)*(gamma_cla-Kappa/2)/(ik*ik*ik));
              
              beta_cla = -(ik+gamma_cla+theta_cla+Kappa/4+((gamma_cla+theta_cla)*Kappa/4+gamma_cla*theta_cla-pow(Kappa,2)/8)/(ik)-((Kappa/2-gamma_cla)*(zeta_cla-Kappa/4)*Kappa/4)/(ik*ik));
              
              delta_cla =-(ik-Kappa/2)/(2*(ik*ik))+(-Kappa/4*(9*Kappa/2-3*gamma_cla)+theta_cla*(gamma_cla+3*Kappa/4))/(2*(ik*ik*ik))+(-3*pow(Kappa,2)/8*(gamma_cla-Kappa)+Kappa*zeta_cla*(gamma_cla-5*Kappa/4))/(2*(ik*ik*ik*ik));
              // a verifier le eta
              eta_cla =(theta_cla+gamma_cla-3*Kappa/4)/(2*(ik*ik*ik))+(-Kappa/2*(gamma_cla-5*Kappa/4)+theta_cla*(gamma_cla-Kappa/2))/((ik*ik*ik*ik))-Kappa/4*(Kappa/4*(7*gamma_cla-13*Kappa/2)+zeta_cla*(5*gamma_cla-11*Kappa/2))/(2*(ik*ik*ik*ik*ik));
              
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if (order_abc == 4)// Km2 omega2
            {
              alpha_cla += (gamma_cla+theta_cla-Kappa/4.0)/ik +
                (Kappa*Kappa/8.0 - (theta_cla+gamma_cla)*Kappa/4.0+(theta_cla*gamma_cla))/(ik*ik);
              
              beta_cla += ( gamma_cla+theta_cla+Kappa/4.0 ) +
                ( ( theta_cla+gamma_cla)*Kappa/4.0 + gamma_cla*theta_cla - Kappa*Kappa/8.0 )/(ik)-
                (Kappa/2.0-gamma_cla)*(zeta_cla-Kappa/4.0)*Kappa/(4.0*ik*ik);
              
              delta_cla = (1.0-Kappa/(2.0*ik))/(2.0*ik);
              
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if(order_abc==5) // Km1 delta2 (censee etre incomplete)
            {
              alpha_cla =1+(gamma_cla-Kappa/4.0)/(ik);
              beta_cla = ik + gamma_cla+Kappa/4.0;
              delta_cla = 1/(2.0*ik) - (Kappa)/(4.0*ik*ik);
              eta_cla = -(gamma_cla-3.0*Kappa/4.0)/(2.0*ik*ik*ik);
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if(order_abc==6)// Km2 delta
            {
              // Km2 delta 0 et 1 (identiques) : tout pourri
              alpha_cla += (gamma_cla+theta_cla-Kappa/4.0)/ik +
                (Kappa*Kappa/8.0 - (theta_cla+gamma_cla)*Kappa/4.0+(theta_cla*gamma_cla))/(ik*ik)-
                ((gamma_cla-Kappa/2.0)*(zeta_cla+Kappa/4.0)*Kappa/4.0)/(ik*ik*ik);
              
              beta_cla += ( gamma_cla+theta_cla+Kappa/4.0 ) +
                ( ( theta_cla+gamma_cla)*Kappa/4.0 + gamma_cla*theta_cla - Kappa*Kappa/8.0 )/(ik)-
                (Kappa/2.0-gamma_cla)*(zeta_cla-Kappa/4.0)*Kappa/(4.0*ik*ik);
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if (order_abc == 7)// Km2 omega1
            {
              alpha_cla += (gamma_cla+theta_cla-Kappa/4.0)/ik ;
              
              beta_cla += ( gamma_cla+theta_cla+Kappa/4.0 ) +
                ( ( theta_cla+gamma_cla)*Kappa/4.0 + gamma_cla*theta_cla - Kappa*Kappa/8.0 )/(ik);
              
              delta_cla = (1.0)/(2.0*ik);
              
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          else if (order_abc == 8) // Km1 omega2
            {
              beta_cla += gamma_cla + Kappa/4.0;
              alpha_cla += (gamma_cla-Kappa/4.0)/(ik);
              delta_cla = (1.0-Kappa/(2.0*ik))/(2.0*ik);
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
            }
          
          if (grazing_abc)
            {
              Complex_wp coeff = pow(k*k,1/3.0)*exp(-Iwp*pi_wp/3.0)*pow(6.0*Kappa,1.0/3.0)*tgamma(2.0/3.0)/tgamma(1.0/3.0);
              if (order_abc == 0)
                {
                  beta_cla = coeff;
                }
              else
                {
                  
                  Complex_wp alpha_cla_orig = alpha_cla;
                  Complex_wp beta_cla_orig = beta_cla;
                  Complex_wp delta_cla_orig = delta_cla;
                  
                  alpha_cla = beta_cla_orig+alpha_cla_orig*(coeff-Kappa);
                  beta_cla = beta_cla_orig*coeff+alpha_cla_orig*ik*ik;
                  delta_cla = alpha_cla_orig+coeff*delta_cla_orig;
                  eta_cla = delta_cla_orig;
                  
                  
                }
              
              beta_abc = (beta_cla + n2*delta_cla) / (alpha_cla + n2*eta_cla);
              
            }
          
          if (nb_coef_pade > 0)
            {
              Complex_wp C0;
              VectComplex_wp CoefAl, CoefBl;
              ComputePadeCoefficientsSqrt(alpha_pade, nb_coef_pade,
					  C0, CoefAl, CoefBl);
              
              beta_abc = C0;
              if (epsilon_pade < 0)
                epsilon_pade = 0.4*pow(abs(k), 1.0/3)*pow(Kappa, 2.0/3);
              
              Complex_wp keps = k + Iwp*epsilon_pade;
              Complex_wp coeff = 1.0/square(keps);
              for (int l = 0; l < nb_coef_pade; l++)
                beta_abc -= CoefAl(l)*n2*coeff/(1.0 - CoefBl(l)*n2*coeff);
                            
              beta_abc *= ik; 
              
              // exact square root
              // beta_abc = ik*sqrt(1.0-n2*coeff);
            }
          
          lin_sys.Get(num, 2*N-2) = Hn_prime(2*N-1)(n) + beta_abc*Hn(2*N-1)(n);
          lin_sys.Get(num, 2*N-1) = Hn2_prime(2*N-1)(n) + beta_abc*Hn2(2*N-1)(n);
          rhs(num) = zero;
        }
      else
        {
          lin_sys.Get(num, 2*N-1) = one;
          rhs(num) = zero;
        }
      
      num++;
      
      if (print_level >= 1)
        cout << "Norm rhs at mode " << n << " = " << Norm2(rhs) << endl;
      
      if (N == 1)
        {
          Complex_wp delta = -lin_sys(1, 0)/lin_sys(1, 1);
          Complex_wp sol0 = rhs(0)/(lin_sys(0, 0) + delta*lin_sys(0, 1)); 
          Complex_wp sol1 = delta*sol0;
          Complex_wp toto(1e30, 1e-12), tutu(2e43, 1e-4);
          if (delta == Complex_wp(1))
            {
              cout << "Cancelling this mode because of unaccuracy" << endl;
              rhs.Fill(0);
            }
          else
            {
              rhs(0) = sol0;
              rhs(1) = sol1;
            }
        }
      else
        {
          //lin_sys.Factorize();
          //lin_sys.Solve(rhs);
          
          GetLU(lin_sys, pivot);
          SolveLU(lin_sys, pivot, rhs);
        }
      
      // alpha and beta are extracted from rhs
      alpha(n).Reallocate(N); beta(n).Reallocate(N);
      for (int i = 0; i < N; i++)
        {
          alpha(n)(i) = rhs(2*i);
          beta(n)(i) = rhs(2*i+1);
        }
    }
  
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout << "Entrez le nom d'un fichier de donnees!" << endl;
      cout << "Usage : ./solution_helm.x fichier_donnees.ini " << endl;
      return -1;
    }
  
  AnalyticalHelmholtzDiscSolver vars;  
  
  string name_input_file(argv[1]);
  ReadInputFile(name_input_file, vars);
  
  if (vars.nb_points_time > 0)
    vars.WriteTimeSolution();
  else
    {
      Vector<VectComplex_wp> alpha, beta;
      vars.ComputeDecompositionCoef(alpha, beta);
      
      vars.WriteSolution(alpha, beta);
    }
  
  return 0;
}
