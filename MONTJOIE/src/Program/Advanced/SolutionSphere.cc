#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

//! class for solving time-harmonic Maxwell's equation in spherical layers
/*!
  On each spherical layer, the solution solves Maxwell's equations
  - \omega^2 \varepsilon u - curl(1 / \mu curl u) = 0
  
  the parameters epsilon and mu are constant in each layer
  Transmission conditions between layers are given as :
  E \times n continuous across the interface
  1/ mu curl E \times n continuous across the interface
  
  The first layer can be a dielectric sphere (i.e. a ball of radius r1 with
  constant epsilon and mu) or dielectric spherical crown, and a Dirichlet,
  Neumann or impedance condition condition is set on the intern boundary r = r0
  
 */
class AnalyticalHarmonicMaxwellSphereSolver : public InputDataProblem_Base
{
public :
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;
  
  //! parameters for output on a predefined grid
  GridInterpolationFull<Dimension3> output_grid;
  //! type of output
  int type_output;
  //! available types of output
  enum {OUTPUT_GRID_ASCII, OUTPUT_GRID_BINARY, OUTPUT_MESH_ASCII, OUTPUT_MESH_BINARY,
        OUTPUT_AXI_ASCII, OUTPUT_AXI_BINARY};
  
  //! if true, we compute the total field, otherwise we compute the diffracted field
  bool total_field;
  
  //! if true we compute the electric field E, otherwise the magnetic field H 
  bool electric_field;
  bool compute_true_magnetic_field; //!< if false we compute -i/omega curl E instead of H

  VectComplex_wp epsilon; //!< dielectric permittivity for each layer
  VectComplex_wp mu; //!< magnetic permeability
  VectReal_wp sigma; //!< conductivity
  VectComplex_wp Z0; //!< sqrt(epsilon/mu) for each layer
  VectComplex_wp kwave; //!< wave number for each layer
  
  //! if true, impedance condition   beta E \times n + n \times H \times n = 0;
  //! if false, dirichlet E \times n = 0
  //! if the first radius is set to 0, no condition is imposed
  bool impedance_cond;  
  Complex_wp beta_impedance; //!< impedance coefficient for impedance boundary condition
  
  //! if true, we replace exact sommerfeld condition by a silver-muller condition placed in an outside sphere
  bool first_order_abc;
  
  //! which data must be computed
  int type_data;
  //! available types of data
  enum { NEAR_FIELD, SURFACIC_FIELD, FAR_FIELD};

  //! number of modes (number of Bessel functions) automatically computed ?
  int type_input_modes;
  enum {MODES_AUTO, MODES_PYTHON, MODES_INPUT, MODES_MANUAL};
  int n_modes; //!< number of modes
  int nb_modes_axi; //! number of modes in axisymmetric output
  // parameters for axisymmetric output
  Real_wp start_angle, finish_angle;
  int nb_angles_axi;

  Vector<Real_wp> ri_sphere; //!< radius for all layers
  Real_wp radius_sphere_outside; //!< radius with silver-muller condition
  //! pulsation
  Real_wp omega, frequency;
  
  //! name of the mesh file
  string file_mesh;
  //! names of output files
  Vector<string> name_file;
  
  string file_RCS; //!< file where RCS is stored
  int nbPointsRCS; //!< number of angles for RCS
  Real_wp first_angle_RCS, last_angle_RCS; //!< first and last angle
  
public:
  AnalyticalHarmonicMaxwellSphereSolver();
  
  // method called for each field of the data file
  void SetInputData(const string& description_field, const VectString& parameters);
  
  // computation of some variables depending on the physical material : Z0, kwave
  void InitMaterial();

  void ComputeSolution(const Vector<VectComplex_wp>& alpha, const Vector<VectComplex_wp>& beta,
                       const Vector<VectComplex_wp>& gamma, const Vector<VectComplex_wp>& delta,
                       const VectR3& Points, Vector<VectComplex_wp>& Usol);
  
  void WriteSolution(const Vector<VectComplex_wp>& alpha, const Vector<VectComplex_wp>& beta,
                     const Vector<VectComplex_wp>& gamma, const Vector<VectComplex_wp>& delta);
  
  void ComputeDecompositionCoef(Vector<VectComplex_wp>& alpha, Vector<VectComplex_wp>& beta,
                                Vector<VectComplex_wp>& gamma, Vector<VectComplex_wp>& delta);

  void ComputeFarField(const Vector<VectComplex_wp>& alpha, const Vector<VectComplex_wp>& beta,
                       const Vector<VectComplex_wp>& gamma, const Vector<VectComplex_wp>& delta,
                       const VectR3& vec_teta, Vector<TinyVector<Complex_wp, 3> >& far_field);
};


//! default constructor
AnalyticalHarmonicMaxwellSphereSolver::AnalyticalHarmonicMaxwellSphereSolver()
{
  type_output = OUTPUT_GRID_BINARY;
  
  total_field = false;
  electric_field = true;
  compute_true_magnetic_field = false;
  
  epsilon.Reallocate(1); epsilon.Fill(1.0);
  mu.Reallocate(1); mu.Fill(1.0);
  Z0.Reallocate(1); Z0.Fill(1.0);
  
  impedance_cond = false; 
  beta_impedance = 0;
  
  first_order_abc = false;
  
  type_data = NEAR_FIELD;
  
  type_input_modes = MODES_AUTO;
  n_modes = 20;
  
  nb_modes_axi = 0;
  start_angle = -0.5*pi_wp;
  finish_angle = 0.5*pi_wp;

  radius_sphere_outside = 2.0; 
  omega = 2.0*pi_wp; frequency = 1.0;
      
  file_mesh = string("sphere.mesh");
  
  file_RCS = string("Rcs.dat");
  nbPointsRCS = 361;
  first_angle_RCS = 0.0; last_angle_RCS = 2.0*pi_wp;
}


//! reads a line of the data file
void AnalyticalHarmonicMaxwellSphereSolver
::SetInputData(const string& description_field, const VectString& parameters)
{
  if (!description_field.compare("TypeOutput"))
    {
      if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
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
      else if (!parameters(0).compare("ModalAxi"))
	{
          if (parameters.GetM() <= 5)
            {
              cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
              cout << "TypeOutput needs more parameters, for instance :" << endl;
              cout << "TypeOutput = ModalAxi BINARY n teta0 teta_f nb_angles" << endl;
              cout << "Current parameters are : " << endl << parameters << endl;
              abort();
            }

	  type_output = OUTPUT_AXI_BINARY;
	  if (!parameters(1).compare("ASCII"))
	    type_output = OUTPUT_AXI_ASCII;
          
	  to_num(parameters(2), nb_modes_axi);
	  to_num(parameters(3), start_angle); start_angle *= pi_wp/180.0;
	  to_num(parameters(4), finish_angle); finish_angle *= pi_wp/180.0;
	  to_num(parameters(5), nb_angles_axi);
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
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
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
  else if (!description_field.compare("ElectricField"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "ElectricField needs more parameters, for instance :" << endl;
          cout << "ElectricField = YES" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      if (!parameters(0).compare("YES"))
	electric_field = true;
      else if (!parameters(0).compare("TRUE_H"))
        {
          electric_field = false;
          compute_true_magnetic_field = true;
        }
      else
	{
          electric_field = false;
          compute_true_magnetic_field = false;
        }
    }
  else if (!description_field.compare("TypeCondition"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "TypeCondition needs more parameters, for instance :" << endl;
          cout << "TypeCondition = DIRICHLET" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      if (!parameters(0).compare("DIRICHLET"))
	{
          impedance_cond = false;
          beta_impedance = 0.0;
        }
      else if (!parameters(0).compare("NEUMANN"))
        {
          impedance_cond = true;
          beta_impedance = 0.0;
        }
      else if (!parameters(0).compare("IMPEDANCE"))
	{
          if (parameters.GetM() <= 1)
            {
              cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
              cout << "TypeCondition needs more parameters, for instance :" << endl;
              cout << "TypeCondition = IMPEDANCE beta" << endl;
              cout << "Current parameters are : " << endl << parameters << endl;
              abort();
            }
          
          impedance_cond = true;
          beta_impedance = to_num<Complex_wp>(parameters(1));
	}
    }
  else if (!description_field.compare("Epsilon"))
    {
      epsilon.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	to_num(parameters(i), epsilon(i));
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
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "FirstOrder_ABC needs more parameters, for instance :" << endl;
          cout << "FirstOrder_ABC = YES" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
      if (!parameters(0).compare("YES"))
	first_order_abc = true;
      else
	first_order_abc = false;
    }
  else if (!description_field.compare("TypeData"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "TypeData needs more parameters, for instance :" << endl;
          cout << "TypeData = NEAR_FIELD" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }
      
      if (!parameters(0).compare("NEAR_FIELD"))
	type_data = NEAR_FIELD;
      else if (!parameters(0).compare("SURFACIC_FIELD"))
	type_data = SURFACIC_FIELD;
      else if (!parameters(0).compare("FAR_FIELD"))
	type_data = FAR_FIELD;
    }
  else if (!description_field.compare("NumberModes"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
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
   else if (!description_field.compare("RadiusSphere"))
    {
      if (parameters.GetM() <= 1)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "RadiusSphere needs more parameters, for instance :" << endl;
          cout << "RadiusSphere = r0 r1" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      ri_sphere.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	to_num(parameters(i), ri_sphere(i));
      
      radius_sphere_outside = ri_sphere(parameters.GetM()-1);
    } 
  else if (!description_field.compare("Frequency"))
    {
      if (parameters.GetM() <= 0)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
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
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "MeshFile needs more parameters, for instance :" << endl;
          cout << "MeshFile = test.mesh" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      file_mesh = string(parameters(0));
    }
  else if (!description_field.compare("OutputFile"))
    {
      if (parameters.GetM() <= 2)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "OutputFile needs more parameters, for instance :" << endl;
          cout << "OutputFile = FileEx FileEy FileEz" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      name_file.Reallocate(parameters.GetM());
      for (int i = 0; i < parameters.GetM(); i++)
	name_file(i) = string(parameters(i)); 
    }
  else if (!description_field.compare("AngleRCS"))
    {
      if (parameters.GetM() <= 2)
        {
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
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
          cout << "In SetInputData of AnalyticalHarmonicMaxwellSphereSolver" << endl;
          cout << "FileRCS needs more parameters, for instance :" << endl;
          cout << "FileRCS = file_name" << endl;
          cout << "Current parameters are : " << endl << parameters << endl;
          abort();
        }

      file_RCS  = string(parameters(0));
    }
}

// basic initialization
void AnalyticalHarmonicMaxwellSphereSolver::InitMaterial()
{
  if ((type_output == OUTPUT_AXI_BINARY) || (type_output == OUTPUT_AXI_ASCII))
    if (type_data != SURFACIC_FIELD)
      {
	cout<<"For modal solutions (axisymmetric decomposition), you have "<<endl;
	cout<<"to choose output of surfacic fields "<<endl;
	abort();
      }
  
  if (mu.GetM() != epsilon.GetM())
    {
      cout<<"Epsilon and mu must have the same size"<<endl;
      abort();
    }
  
  if (mu.GetM() != ri_sphere.GetM()-1)
    {
      cout<<"You have to specify " << mu.GetM()+1 << " radii " << endl;
      cout<<"You have specified only " << ri_sphere.GetM() << " radii" << endl;
      abort();
    }

  int nb_layers = mu.GetM();
  Z0.Reallocate(nb_layers);
  kwave.Reallocate(nb_layers);
  for (int i = 0;  i < nb_layers; i++)
    {
      Complex_wp eps = epsilon(i);
      if (i < sigma.GetM())
        eps += Iwp*sigma(i)/omega;
      
      Z0(i) = sqrt(mu(i)/eps);
      kwave(i) = omega*sqrt(mu(i)*eps);
    }    
}

//! computes the analytical solution given by coefficients alpha, beta, gamma and delta
//! on given points
void AnalyticalHarmonicMaxwellSphereSolver
::ComputeSolution(const Vector<VectComplex_wp>& alpha, const Vector<VectComplex_wp>& beta,
                  const Vector<VectComplex_wp>& gamma, const Vector<VectComplex_wp>& delta,
                  const VectR3& Points, Vector<VectComplex_wp>& Usol)
{  
  Complex_wp zero, one;
  SetComplexZero(zero); SetComplexOne(one);
  
  int N_max = n_modes;
  LegendrePolynomial<Real_wp> Pn(N_max);

  VectComplex_wp Jn(N_max), Hn(N_max), Hn2(N_max), Xi2_n(N_max), Xi_n(N_max), Psi_n(N_max);
  VectComplex_wp Jn_prime(N_max), Hn_prime(N_max), Xi_prime(N_max), Xi2_prime(N_max), Psi_prime(N_max);
  VectReal_wp Pn0(N_max), Pn1(N_max), dPn1_dteta(N_max);
  
  // coordinates (cartersian, spherical, sinus and cosinus of angles)
  Real_wp x, y, z, r, teta, phi, sin_teta, cos_teta, cos_phi, sin_phi, invSinTeta;
  
  //int nbPoints_phi = toInteger(ceil(2.0*pi_wp*radius_sphere_outside*omega));
  //Globatto<Real_wp> gauss_phi;
  //if ((type_output == OUTPUT_AXI_BINARY) || (type_output == OUTPUT_AXI_ASCII))
  //gauss_phi.ConstructQuadrature(nbPoints_phi-1, gauss_phi.QUADRATURE_GAUSS);
  
  // components of the electrical field
  Complex_wp E_r, E_teta, E_phi, kr, E_x, E_y, E_z;
  Real_wp begin_order = 0.5; //Real_wp r_cyl, teta_cyl; int num_phi;
  Complex_wp value_Jtm1(0, 0), value_Jbm1(0, 0), value_Jtp1(0, 0), value_Jbp1(0, 0);
  
  // coefficients precomputed for more efficiency and lisibility
  VectComplex_wp coef_Er(N_max), coef_Ez(N_max);
  Complex_wp mI_exp_n = -Iwp;
  for (int n = 1; n < N_max; n++)
    {
      Real_wp coef_n = Real_wp(2*n+1)/Real_wp(n*(n+1));
      coef_Er(n) = mI_exp_n*Real_wp(2*n+1);
      coef_Ez(n) = mI_exp_n*coef_n;
      mI_exp_n *= -Iwp;
    }

  // allocating arrays
  Usol.Reallocate(3);
  Usol(0).Reallocate(Points.GetM()); Usol(1).Reallocate(Points.GetM());
  Usol(2).Reallocate(Points.GetM());
  
  // loop on all vertices where the solution is computed
  for (int i = 0; i < Points.GetM(); i++)
    {
      // coordinate of the points
      x = Points(i)(0); 
      y = Points(i)(1); 
      z = Points(i)(2);

      // conversion to spherical coordinates
      CartesianToSpherical(x, y, z, r, teta, phi, cos_teta, sin_teta);
      cos_phi = cos(phi); sin_phi = sin(phi); invSinTeta = Real_wp(1)/sin_teta;
      
      // we find the layer where the point (x, y, z) is
      int num = ri_sphere.GetM()-1;
      for (int n = 0; n < ri_sphere.GetM(); n++)
        {
          Real_wp coef_debordement = 1.0-50*epsilon_machine;
          if (n == ri_sphere.GetM() - 1)
            coef_debordement = 1.0+50*epsilon_machine;
          
          if (r < coef_debordement*ri_sphere(n))
            {
              num = n-1;
              break;
            }
        }
      
      SetComplexZero(E_r); SetComplexZero(E_teta); SetComplexZero(E_phi);
      SetComplexZero(E_x); SetComplexZero(E_y); SetComplexZero(E_z);
      if ((num >= 0) && (num < ri_sphere.GetM()-1))
        {
	  kr = kwave(num) * r;
	  
          // Spherical Bessel functions and derivatives are computed
  	  ComputeDeriveRiccatiBessel(begin_order, N_max, kr, Jn, Hn, Jn_prime, Hn_prime,
                                     Psi_n, Xi_n, Psi_prime, Xi_prime);
	  
	  for (int n = 1; n < N_max; n++)
	    {
	      Hn2(n) = -Hn(n) + 2.0*Jn(n);
	      Xi2_n(n) = -Xi_n(n) + 2.0*Psi_n(n);
	      Xi2_prime(n) = -Xi_prime(n) + 2.0*Psi_prime(n);
	    }
          
          Pn.EvaluateDerivative(N_max, cos_teta, sin_teta, Pn0, Pn1, dPn1_dteta);
          
          Complex_wp coef_grad, coef_curl;
          
	  // electric_field
	  if (electric_field)
	    {
              Complex_wp cosPhi_kr = cos_phi / kr;
              Complex_wp sinPhi_kr = sin_phi / kr;

              if (num == ri_sphere.GetM() - 2)
                for (int n = 1; n < N_max; n++)
                  {
                    E_r += cosPhi_kr * coef_Er(n) * Iwp * (alpha(n)(num)*Hn(n) + gamma(n)(num)*Hn2(n) + Jn(n)) * Pn1(n);
                    
                    coef_grad = Iwp*coef_Ez(n)*(alpha(n)(num)*Xi_prime(n) + gamma(n)(num)*Xi2_prime(n) + Psi_prime(n));
                    coef_curl = coef_Ez(n)*(beta(n)(num)*Xi_n(n) + delta(n)(num)*Xi2_n(n) + Psi_n(n));
                    
                    E_teta += cosPhi_kr * (coef_grad * dPn1_dteta(n) + coef_curl * Pn1(n) * invSinTeta);
                    E_phi += -sinPhi_kr * (coef_grad * Pn1(n) * invSinTeta + coef_curl * dPn1_dteta(n));
                  }
              else
                for (int n = 1; n < N_max; n++)
                  {
                    E_r += cosPhi_kr * coef_Er(n) * Iwp * (alpha(n)(num)*Hn(n) + gamma(n)(num)*Hn2(n)) * Pn1(n);
                    
                    coef_grad = Iwp*coef_Ez(n)*(alpha(n)(num)*Xi_prime(n) + gamma(n)(num)*Xi2_prime(n));
                    coef_curl = coef_Ez(n)*(beta(n)(num)*Xi_n(n) + delta(n)(num)*Xi2_n(n));
                    
                    E_teta += cosPhi_kr * (coef_grad * dPn1_dteta(n) + coef_curl * Pn1(n) * invSinTeta);
                    E_phi += -sinPhi_kr * (coef_grad * Pn1(n) * invSinTeta + coef_curl * dPn1_dteta(n));
                  }
            }
          else
            {
              Complex_wp cosPhi_kr = cos_phi / (Z0(num)*kr);
              Complex_wp sinPhi_kr = sin_phi / (Z0(num)*kr);
              if (!compute_true_magnetic_field)
                {
                  // multiplication by mu in order to compute -i/omega curl E instead of H
                  cosPhi_kr *= mu(num);
                  sinPhi_kr *= mu(num);
                }
              
              if (num == ri_sphere.GetM() - 2)
                for (int n = 1; n < N_max; n++)
                  {
                    E_r += -sinPhi_kr * Iwp * coef_Er(n) * (beta(n)(num)*Hn(n) + delta(n)(num)*Hn2(n) + Jn(n)) * Pn1(n);
                    
                    coef_grad = Iwp*coef_Ez(n)*(beta(n)(num)*Xi_prime(n) + delta(n)(num)*Xi2_prime(n) + Psi_prime(n));
                    coef_curl = coef_Ez(n)*(alpha(n)(num)*Xi_n(n) + gamma(n)(num)*Xi2_n(n) + Psi_n(n));
                    
                    E_teta += sinPhi_kr * (-coef_grad * dPn1_dteta(n) - coef_curl * Pn1(n) * invSinTeta);
                    E_phi += cosPhi_kr * (-coef_grad * Pn1(n) * invSinTeta - coef_curl * dPn1_dteta(n));
                  }
              else
                for (int n = 1; n < N_max; n++)
                  {                    
                    E_r += -sinPhi_kr * Iwp * coef_Er(n) * (beta(n)(num)*Hn(n) + delta(n)(num)*Hn2(n)) * Pn1(n);
                    
                    coef_grad = Iwp*coef_Ez(n)*(beta(n)(num)*Xi_prime(n) + delta(n)(num)*Xi2_prime(n));
                    coef_curl = coef_Ez(n)*(alpha(n)(num)*Xi_n(n) + gamma(n)(num)*Xi2_n(n));
                    
                    E_teta += sinPhi_kr * (-coef_grad * dPn1_dteta(n) - coef_curl * Pn1(n) * invSinTeta);
                    E_phi += cosPhi_kr * (-coef_grad * Pn1(n) * invSinTeta - coef_curl * dPn1_dteta(n));
                  }
            }
          
          E_x = E_r*sin_teta*cos_phi + E_teta*cos_teta*cos_phi - E_phi*sin_phi;
          E_y = E_r*sin_teta*sin_phi + E_teta*cos_teta*sin_phi + E_phi*cos_phi;
          E_z = E_r*cos_teta - E_teta*sin_teta;
          
          // then subtracting the plane wave for the diffracted field
          if (!total_field)
            {
              if (electric_field)
                E_x -= exp(-Iwp*kwave(ri_sphere.GetM()-2)*z);
              else
                E_y += one/Z0(ri_sphere.GetM()-2)*exp(-Iwp*kwave(ri_sphere.GetM()-2)*z);
            }
          
          if (type_data == SURFACIC_FIELD)
            {
              // we compute E \times n or H \times n
              // n is equal to -\hat{r}
              E_x = -E_phi*cos_teta*cos_phi - E_teta*sin_phi;
              E_y = -E_phi*cos_teta*sin_phi + E_teta*cos_phi;
              E_z = E_phi*sin_teta;
            }
          
          
          if ((type_output == OUTPUT_AXI_BINARY) || (type_output == OUTPUT_AXI_ASCII))
            {
              cout << "not implemented " << endl;
              abort();
            }      
        }
      
      Usol(0)(i) = E_x; Usol(1)(i) = E_y; Usol(2)(i) = E_z;
    }

}


void AnalyticalHarmonicMaxwellSphereSolver
::ComputeFarField(const Vector<VectComplex_wp>& alpha, const Vector<VectComplex_wp>& beta,
                  const Vector<VectComplex_wp>& gamma, const Vector<VectComplex_wp>& delta,
                  const VectR3& vec_teta, Vector<TinyVector<Complex_wp, 3> >& far_field)
{
  far_field.Reallocate(vec_teta.GetM());
  
  int N_max = n_modes;
  int N = epsilon.GetM();
  LegendrePolynomial<Real_wp> Pn(N_max);
  
  Complex_wp zero; SetComplexZero(zero);
  Complex_wp E_r, E_teta, E_phi, E_x, E_y, E_z;
  VectReal_wp Pn0(N_max), Pn1(N_max), dPn1_dteta(N_max);
  Real_wp r, teta, phi, cos_teta, sin_teta, cos_phi, sin_phi, invSinTeta;
  for (int i = 0; i < vec_teta.GetM(); i++)
    {
      // conversion to spherical coordinates
      CartesianToSpherical(vec_teta(i)(0), vec_teta(i)(1), vec_teta(i)(2),
                           r, teta, phi, cos_teta, sin_teta);
      
      cos_phi = cos(phi); sin_phi = sin(phi); invSinTeta = Real_wp(1)/sin_teta;
      
      Pn.EvaluateDerivative(N_max, cos_teta, sin_teta, Pn0, Pn1, dPn1_dteta);
      
      // equivalent when r tends to the infinity (E_r is dropped ?)
      E_r = zero; E_teta = zero; E_phi = 0;
      Complex_wp mI_exp_n = -Iwp;
      Complex_wp Xi_prime_equiv = -Iwp;
      Complex_wp Xi_n_equiv = mI_exp_n*mI_exp_n;
      for (int n = 1; n < N_max; n++)
	{
	  Real_wp coef_n = Real_wp(2*n+1)/Real_wp(n*(n+1));
	  E_teta += mI_exp_n * coef_n *(-Iwp * beta(n)(N-1) * Xi_prime_equiv * dPn1_dteta(n) -
                                        alpha(n)(N-1)*Xi_n_equiv * Pn1(n) * invSinTeta ) * sin_phi;
	  
          E_phi += mI_exp_n*coef_n*(-Iwp * beta(n)(N-1)*Xi_prime_equiv * Pn1(n) * invSinTeta -
                                    alpha(n)(N-1) * Xi_n_equiv * dPn1_dteta(n) ) * cos_phi;
	  
          mI_exp_n *= -Iwp;
	  Xi_prime_equiv *= -Iwp; Xi_n_equiv *= -Iwp; 
	}
      
      E_x = E_r*sin_teta*cos_phi + E_teta*cos_teta*cos_phi - E_phi*sin_phi;
      E_y = E_r*sin_teta*sin_phi + E_teta*cos_teta*sin_phi + E_phi*cos_phi;
      E_z = E_r*cos_teta - E_teta*sin_teta;
      
      far_field(i).Init(E_x, E_y, E_z);
    }
  
}

//! computes the solution given by coefficients alpha and beta, and writes it on the asked output
void AnalyticalHarmonicMaxwellSphereSolver
::WriteSolution(const Vector<VectComplex_wp>& alpha, const Vector<VectComplex_wp>& beta,
                const Vector<VectComplex_wp>& gamma, const Vector<VectComplex_wp>& delta)
{
  if (type_data == FAR_FIELD)
    {
      VectR3 vec_teta(nbPointsRCS);
      Vector<TinyVector<Complex_wp, 3> > far_field;
      Real_wp teta_rcs, step_angle;
      step_angle = (last_angle_RCS - first_angle_RCS) / (nbPointsRCS-1);      
      for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
        {
          teta_rcs = first_angle_RCS + step_angle*num_angle;
          vec_teta(num_angle).Init(cos(teta_rcs), Real_wp(0), sin(teta_rcs));
        }
      
      ComputeFarField(alpha, beta, gamma, delta, vec_teta, far_field);

      ofstream file_out(file_RCS.data());
      Real_wp module_E, value_rcs;
      for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
        {
          teta_rcs = first_angle_RCS + step_angle*num_angle;
          
          module_E = Norm2(far_field(num_angle));
          value_rcs = 10*log10(omega*omega/(4*pi_wp)*module_E*module_E);
          file_out << 180/pi_wp*teta_rcs << " " << value_rcs << '\n';
        }
      
      file_out.close();

      return;
    }
  
  Vector<VectComplex_wp> Usol;
  VectR3 Points;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_GRID_BINARY))
    output_grid.GenerateGridPoints(Points);
  else
    {
      Mesh<Dimension3> mesh;
      mesh.Read(file_mesh);
      Points = mesh.Vertex();
    }
  
  ComputeSolution(alpha, beta, gamma, delta, Points, Usol);
  
  int double_prec = OutputTypeEnum::DOUBLE_PRECISION;
  bool ascii = false;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_MESH_ASCII))
    ascii = true;
  
  Dimension3 dim;
  if ((type_output == OUTPUT_GRID_ASCII) || (type_output == OUTPUT_GRID_BINARY))
    {
      for (int k = 0; k < Usol.GetM(); k++)
        WriteMatlab(Usol(k), output_grid, name_file(k), double_prec, ascii);
    }
  else
    {
      for (int k = 0; k < Usol.GetM(); k++)
        WriteMedit(Usol(k), name_file(k), dim, double_prec);
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
void AnalyticalHarmonicMaxwellSphereSolver
::ComputeDecompositionCoef(Vector<VectComplex_wp>& alpha, Vector<VectComplex_wp>& beta,
                           Vector<VectComplex_wp>& gamma, Vector<VectComplex_wp>& delta)
{
  InitMaterial();
  
  if (type_input_modes == MODES_PYTHON)
    {
      ofstream file_out("temp.py");
      ifstream file_in("MATLAB/maxwell_sphere.py");
      file_out.precision(15);
      
      if (!file_in.is_open())
	{
	  cout << "Unable to open file MATLAB/maxwell_sphere.py" << endl;
	  abort();
	}
      
      string line;
      while (getline(file_in, line))
	{
	  if (StartsWith(line, "epsilon = "))
	    {
	      file_out << "epsilon = ";
	      WriteMatrixData(file_out, epsilon);
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
		  sigma.Reallocate(epsilon.GetM());
		  sigma.Fill(0);
		}

	      file_out << "sigma = ";
	      WriteMatrixData(file_out, sigma);
	    }
	  else if (StartsWith(line, "radius = "))
	    {
	      file_out << "radius = ";
	      WriteMatrixData(file_out, ri_sphere);
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
      Matrix<Complex_wp> A, B, C, D;
      A.ReadText("alpha.dat"); B.ReadText("beta.dat");
      C.ReadText("gamma.dat"); D.ReadText("delta.dat");
      
      int N_max = A.GetM();
      n_modes = N_max;
      alpha.Reallocate(N_max); beta.Reallocate(N_max);
      gamma.Reallocate(N_max); delta.Reallocate(N_max);
      for (int n = 1; n < N_max; n++)
	{
	  alpha(n).Reallocate(A.GetN());
	  beta(n).Reallocate(A.GetN());
	  gamma(n).Reallocate(A.GetN());
	  delta(n).Reallocate(A.GetN());
	  for (int i = 0; i < A.GetN(); i++)
	    {
	      alpha(n)(i) = A(n, i);
	      beta(n)(i) = B(n, i);
	      gamma(n)(i) = C(n, i);
	      delta(n)(i) = D(n, i);
	    }
	}
      
      //DISP(alpha); DISP(gamma); DISP(beta); DISP(delta);
      return;
    }
  
  // number of layers
  int N = epsilon.GetM();

  // maximum value of kr for incident wave
  Real_wp ka = abs(kwave(N-1))*ri_sphere(N);
  
  // number of Bessel functions to take into account
  int N_max = 0;
  if (type_input_modes == MODES_AUTO)
    {
      N_max = ComputeOrder(ka, epsilon_machine);
      n_modes = N_max;
    }
  else
    N_max = n_modes;
  
  Real_wp begin_order = 0.5;
  Complex_wp zero, kr, one; SetComplexZero(zero);
  SetComplexOne(one);
  
  // computes the spherical bessel functions for each layer (and at the two radii of the layer)
  Vector<VectComplex_wp> Jn(2*N), Hn(2*N), Hn2(2*N), Jn_prime(2*N), Hn_prime(2*N), Hn2_prime(2*N);
  Vector<VectComplex_wp> Psi_n(2*N), Xi_n(2*N), Xi2_n(2*N), Psi_prime(2*N), Xi_prime(2*N), Xi2_prime(2*N);
  for (int i = 0; i < N; i++)
    {
      kr = kwave(i)*ri_sphere(i);
      
      // Spherical Bessel functions and derivatives are computed
      ComputeDeriveRiccatiBessel(begin_order, N_max, kr, Jn(2*i), Hn(2*i), Jn_prime(2*i), Hn_prime(2*i),
                                 Psi_n(2*i), Xi_n(2*i), Psi_prime(2*i), Xi_prime(2*i));
      
      // Spherical Hankel functions of the second kind
      Hn2(2*i).Reallocate(N_max);
      Hn2_prime(2*i).Reallocate(N_max);
      Xi2_n(2*i).Reallocate(N_max);
      Xi2_prime(2*i).Reallocate(N_max);
      for (int n = 0; n < N_max; n++)
        {
          Hn2(2*i)(n) = -Hn(2*i)(n) + 2.0*Jn(2*i)(n);
          Hn2_prime(2*i)(n) = -Hn_prime(2*i)(n) + 2.0*Jn_prime(2*i)(n);
          Xi2_n(2*i)(n) = -Xi_n(2*i)(n) + 2.0*Psi_n(2*i)(n);
          Xi2_prime(2*i)(n) = -Xi_prime(2*i)(n) + 2.0*Psi_prime(2*i)(n);
        }

      kr = kwave(i)*ri_sphere(i+1);
      
      // Spherical Bessel functions and derivatives are computed
      ComputeDeriveRiccatiBessel(begin_order, N_max, kr,
                                 Jn(2*i+1), Hn(2*i+1), Jn_prime(2*i+1), Hn_prime(2*i+1),
                                 Psi_n(2*i+1), Xi_n(2*i+1), Psi_prime(2*i+1), Xi_prime(2*i+1));
      
      // Spherical Hankel functions of the second kind
      Hn2(2*i+1).Reallocate(N_max);
      Hn2_prime(2*i+1).Reallocate(N_max);
      Xi2_n(2*i+1).Reallocate(N_max);
      Xi2_prime(2*i+1).Reallocate(N_max);
      for (int n = 0; n < N_max; n++)
        {
          Hn2(2*i+1)(n) = -Hn(2*i+1)(n) + 2.0*Jn(2*i+1)(n);
          Hn2_prime(2*i+1)(n) = -Hn_prime(2*i+1)(n) + 2.0*Jn_prime(2*i+1)(n);
          Xi2_n(2*i+1)(n) = -Xi_n(2*i+1)(n) + 2.0*Psi_n(2*i+1)(n);
          Xi2_prime(2*i+1)(n) = -Xi_prime(2*i+1)(n) + 2.0*Psi_prime(2*i+1)(n);
        }
    }
  
  Matrix<Complex_wp> sys_alpha(2*N, 2*N), sys_beta(2*N, 2*N);
  Vector<Complex_wp> rhs_alpha(2*N), rhs_beta(2*N);
  IVect pivot_alpha(2*N), pivot_beta(2*N);

  // each mode is decoupled
  alpha.Reallocate(N_max); beta.Reallocate(N_max);
  gamma.Reallocate(N_max); delta.Reallocate(N_max);
  Complex_wp coef; DISP(N_max);
  for (int n = 1; n < N_max; n++)
    {
      sys_alpha.Fill(zero); rhs_alpha.Fill(zero);
      sys_beta.Fill(zero); rhs_beta.Fill(zero);
      
      // first equation derived from the condition on r = r0
      if (ri_sphere(0) <= epsilon_machine)
        {
          // dielectric case, alpha_n = gamma_n, beta_n = delta_n
          // since the solution involves only Bessel function of the first kind psi
          sys_alpha.Get(0, 0) = one; sys_alpha.Get(0, 1) = -one;
          sys_beta.Get(0, 0) = one; sys_beta.Get(0, 1) = -one;
        }
      else if (impedance_cond)
        {
          // impedance condition n/mu \times curl E - beta E_tangential = 0
          coef = omega / Z0(0);
          sys_alpha.Get(0, 0) = beta_impedance*Xi_prime(0)(n) + coef*Xi_n(0)(n);
          sys_alpha.Get(0, 1) = beta_impedance*Xi2_prime(0)(n) + coef*Xi2_n(0)(n);
          if (N == 1)
            rhs_alpha(0) = -beta_impedance*Psi_prime(0)(n) - coef*Psi_n(0)(n);
          
          sys_beta.Get(0, 0) = -coef*Xi_prime(0)(n) + beta_impedance*Xi_n(0)(n);
          sys_beta.Get(0, 1) = -coef*Xi2_prime(0)(n) + beta_impedance*Xi2_n(0)(n);
          if (N == 1)
            rhs_beta(0) = coef*Psi_prime(0)(n) - beta_impedance*Psi_n(0)(n);          
        }
      else
        {
          // Dirichlet condition, E_tangential = 0
          sys_alpha.Get(0, 0) = Iwp*Xi_prime(0)(n);
          sys_alpha.Get(0, 1) = Iwp*Xi2_prime(0)(n);
          if (N == 1)
            rhs_alpha(0) = -Iwp*Psi_prime(0)(n);
          
          sys_beta.Get(0, 0) = Xi_n(0)(n);
          sys_beta.Get(0, 1) = Xi2_n(0)(n);
          if (N == 1)
            rhs_beta(0) = -Psi_n(0)(n);
        }
      
      int num = 1;
      
      // loop over transmission conditions (two equations by interface)
      for (int i = 1; i < ri_sphere.GetM()-1; i++)
        {
          // continuity of the tangential trace of eletric field
          sys_alpha.Get(num, 2*i-2) = Xi_prime(2*i-1)(n) / kwave(i-1);
          sys_alpha.Get(num, 2*i-1) = Xi2_prime(2*i-1)(n) / kwave(i-1);
          sys_alpha.Get(num, 2*i) = -Xi_prime(2*i)(n) / kwave(i);
          sys_alpha.Get(num, 2*i+1) = -Xi2_prime(2*i)(n) / kwave(i);
          if (i == N-1)
            rhs_alpha(num) = Psi_prime(2*i)(n) / kwave(i);

          sys_beta.Get(num, 2*i-2) = Xi_n(2*i-1)(n) / kwave(i-1);
          sys_beta.Get(num, 2*i-1) = Xi2_n(2*i-1)(n) / kwave(i-1);
          sys_beta.Get(num, 2*i) = -Xi_n(2*i)(n) / kwave(i);
          sys_beta.Get(num, 2*i+1) = -Xi2_n(2*i)(n) / kwave(i);
          if (i == N-1)
            rhs_beta(num) = Psi_n(2*i)(n) / kwave(i);
          
          num++;
          
          // continuity of the tangential trace of magnetic field
          sys_alpha.Get(num, 2*i-2) = Xi_n(2*i-1)(n) / (kwave(i-1)*Z0(i-1));
          sys_alpha.Get(num, 2*i-1) = Xi2_n(2*i-1)(n) / (kwave(i-1)*Z0(i-1));
          sys_alpha.Get(num, 2*i) = -Xi_n(2*i)(n) / (kwave(i)*Z0(i));
          sys_alpha.Get(num, 2*i+1) = -Xi2_n(2*i)(n) / (kwave(i)*Z0(i));
          if (i == N-1)
            rhs_alpha(num) = Psi_n(2*i)(n) / (kwave(i)*Z0(i));

          sys_beta.Get(num, 2*i-2) = Xi_prime(2*i-1)(n) / (kwave(i-1)*Z0(i-1));
          sys_beta.Get(num, 2*i-1) = Xi2_prime(2*i-1)(n) / (kwave(i-1)*Z0(i-1));
          sys_beta.Get(num, 2*i) = -Xi_prime(2*i)(n) / (kwave(i)*Z0(i));
          sys_beta.Get(num, 2*i+1) = -Xi2_prime(2*i)(n) / (kwave(i)*Z0(i));
          if (i == N-1)
            rhs_beta(num) = Psi_prime(2*i)(n) / (kwave(i)*Z0(i));
          
          num++;
        }

      // last condition due to Sommerfeld condition (at finite distance or not)
      if (first_order_abc)
        {
          sys_alpha.Get(num, 2*N-2) = Xi_prime(2*N-1)(n) - Iwp*Xi_n(2*N-1)(n);
          sys_alpha.Get(num, 2*N-1) = Xi2_prime(2*N-1)(n) - Iwp*Xi2_n(2*N-1)(n);
          
          sys_beta.Get(num, 2*N-2) = Xi_prime(2*N-1)(n) - Iwp*Xi_n(2*N-1)(n);
          sys_beta.Get(num, 2*N-1) = Xi2_prime(2*N-1)(n) - Iwp*Xi2_n(2*N-1)(n);
        }
      else
        {
          // exact Sommerfeld condition => no component in Xi_n^{(2)}
          sys_alpha.Get(num, 2*N-1) = one;
          sys_beta.Get(num, 2*N-1) = one;
        }
      
      num++;
            
      cout << "Norm rhs at mode " << n << " = " << Norm2(rhs_alpha) << " , " << Norm2(rhs_beta) << endl;
      //DISP(rhs_alpha); DISP(rhs_beta); DISP(sys_alpha); DISP(sys_beta);

      GetLU(sys_alpha, pivot_alpha);
      SolveLU(sys_alpha, pivot_alpha, rhs_alpha);

      GetLU(sys_beta, pivot_beta);
      SolveLU(sys_beta, pivot_beta, rhs_beta);
      
      // alpha, beta, gamma and delta are extracted from rhs
      alpha(n).Reallocate(N); beta(n).Reallocate(N);
      gamma(n).Reallocate(N); delta(n).Reallocate(N);
      //DISP(rhs_alpha); DISP(rhs_beta);
      for (int i = 0; i < N; i++)
        {
          alpha(n)(i) = rhs_alpha(2*i);
          gamma(n)(i) = rhs_alpha(2*i+1);
          beta(n)(i) = rhs_beta(2*i);
          delta(n)(i) = rhs_beta(2*i+1);
        }
    }
  
  //DISP(alpha); DISP(gamma); DISP(beta); DISP(delta);
}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc < 2)
    {
      cout<<"Entrez le nom d'un fichier de donnees!"<<endl;
      abort();
    }
  
  AnalyticalHarmonicMaxwellSphereSolver vars;
  
  string name_input_file(argv[1]);
  ReadInputFile(name_input_file, vars);
  
  Vector<VectComplex_wp> alpha, beta, gamma, delta;
  vars.ComputeDecompositionCoef(alpha, beta, gamma, delta);

  vars.WriteSolution(alpha, beta, gamma, delta);
  
  return FinalizeMontjoie();
}
