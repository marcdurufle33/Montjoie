#define MONTJOIE_WITH_ONE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Harmonic/MontjoieHarmonic.hxx"
#include "Instationary/MontjoieTime.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Maxwell/NonLinearOpticsProblem.hxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D.hxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D_Inline.cxx"
#include "Elliptic/Helmholtz/SchrodingerNonLinear1D.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#include "Elliptic/Maxwell/NonLinearOpticsProblem.cxx"
#include "Elliptic/Helmholtz/NonLinearMaxwell1D.cxx"
#include "Elliptic/Helmholtz/SchrodingerNonLinear1D.cxx"
#endif

using namespace Montjoie;

#include <xmmintrin.h>

#define CSR_FLUSH_TO_ZERO         (1 << 15)

// probleme harmonique
void RunHarmonicProblem(MaxwellProblem1D& var, const string& input_file)
{
  typedef Montjoie::Complex_wp Complex_wp;

  // setting resolution in time-harmonic domain
  var.SetHarmonicResolution();

  // constructing needed arrays to compute the solution
  var.ConstructAll(input_file);
  
  // right hand side
  Vector<Complex_wp> source_rhs, x_sol;
  var.ComputeRightHandSide(source_rhs);
  
  // computing the solution
  var.SolveNewton(source_rhs, x_sol);
  
  // writing the solution
  var.WriteDatas(x_sol);
}

// probleme temporel
void RunTimeProblem(MaxwellProblem1D& var, const string& input_file)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;

  // reading parameters of the time scheme
  All_TimeScheme<Complex_wp> scheme;
  ReadInputFile(input_file, scheme);
  bool explicit_scheme = true;
  if (scheme.GetTimeSchemeType() == TimeSchemeEnum::THETA_SCHEME)
    explicit_scheme = false;
  
  // all what is needed to complete time iterations is constructed
  var.ConstructAll(input_file);
  
  // initialisation of time scheme
  int N = var.SetInitialCondition(explicit_scheme);
  Complex_wp zero; SetComplexZero(zero);
  if (explicit_scheme)
    {
      VectComplex_wp E0(N), E0_prime;
      E0.Fill(zero);
      //var.FillInitialCondition(E0);
      
      //E0.FillRand(); Mlt(1e-200, E0);
      if (!scheme.FirstOrderScheme())
        {
          E0_prime.Reallocate(N);
          E0_prime.Fill(zero);
	  scheme.SetInitialCondition(var.GetInitialTime(), var.GetTimeStep(), E0, E0_prime, var);
          //E0_prime.FillRand(); Mlt(1e-200, E0_prime);
        }
      else
	scheme.SetInitialCondition(var.GetInitialTime(), var.GetTimeStep(), E0, var);
    }
  
  // time iterations
  int nb_iterations_time = toInteger(ceil((var.GetFinalTime()-var.GetInitialTime())/var.GetTimeStep()));
  Real_wp t; // DISP(PhysicalConstant::adimensionalization);
  MontjoieTimer chrono; chrono.Start(VirtualTimer::ALL);
  for (int nt = 0; nt <= nb_iterations_time; nt++)
    {
      // current time t^n
      t = var.GetInitialTime() + nt*var.GetTimeStep();
      
      // outputs are done if needed
      var.WriteSnapshot(nt, t, var.GetIterate());
      //DISP(Norm2(var.GetIterate()));
      //DISP(Norm2(var.GetIterate()));
      // commented lines are used for "debug" (see check_maxwell1D.py)
      /* if (abs(nt-100000) <= 1)
	{
	  DISP(var.E0_adim);
	  int p = nt-99999;
	  VectComplex_wp& Udof = scheme.GetIterate();
	  VectComplex_wp Uinterp(var.var_laplace.var_section.GetNbPointsGrid()), Usub(var.nb_points_z);
	  int Nu = Udof.GetM()/var.nb_points_z;
	  int Nvec = var.omega_polarization.GetM()+1;
	  int nb_pts = var.var_laplace.var_section.GetNbPointsGrid();
	  VectComplex_wp U(nb_pts*Nvec); U.Fill(zero);
	  for (int i = 0; i < Nvec; i++)
	    {
	      if (i == 0)
		Usub = var.GetIterate();
	      else
		{
		  if (i < Nu)
		    for (int j = 0; j < var.nb_points_z; j++)
		      Usub(j) = Udof(i*var.nb_points_z + j);
		  else
		    for (int j = 0; j < var.nb_points_z; j++)
		      Usub(j) = var.P_current((i-1)*var.nb_points_z + j);
		}
	      
	      var.var_laplace.ComputeInterpolationU(Usub, var.var_laplace.var_section, Uinterp);
	      
	      for (int j = 0; j < nb_pts; j++)
		U(i*nb_pts + j) = Uinterp(j);
	    }
	  
	  U.Write("Un" + to_str(p) + ".dat");
	  
	  if (var.var_laplace.mixed_formulation)
	    {
	      int offset = Nvec*var.nb_points_z;
	      for (int i = 0; i < Nvec; i++)
		{
		  if (i == 0)
		    {
		      Usub.Reallocate(var.H_current.GetM());
		      for (int j = 0; j < var.H_current.GetM(); j++)
			Usub(j) = Udof(offset + j);
		    }
		  else
		    {
		      Usub.Reallocate(var.nb_points_z);
		      for (int j = 0; j < var.nb_points_z; j++)
			Usub(j) = Udof(offset + j);
		    }
		  
		  if (i == 0)
		    var.var_laplace.ComputeInterpolationU(Usub, var.var_laplace.var_section, Uinterp, true);
		  else
		    var.var_laplace.ComputeInterpolationU(Usub, var.var_laplace.var_section, Uinterp);
		  
		  for (int j = 0; j < nb_pts; j++)
		    U(i*nb_pts + j) = Uinterp(j);
		  
		  if (i == 0)
		    offset += var.H_current.GetM();
		  else
		    offset += var.nb_points_z;
		}
	      
	      U.Write("Hn" + to_str(p) + ".dat");	      
	    }
	}
      */
      
      // time scheme is advanced
      var.AdvanceScheme(nt, t, scheme);
      
      //if (nt%10000 == 0)
      //cout << "real time elapsed " << chrono.GetSeconds(VirtualTimer::ALL) << endl;
    }
  
  var.CloseBuffers();
}


int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Donnez un fichier de donnees" << endl;
      abort();
    }
  
  InitMontjoie(argc, argv);
  //cout.setf(ios::scientific);
  
  // in order to avoid denormal values
  unsigned csr = __builtin_ia32_stmxcsr();
  csr |= CSR_FLUSH_TO_ZERO;
  __builtin_ia32_ldmxcsr(csr);
  
  //_mm_setcsr( _mm_getcsr() | 0x8040 );
  
  string input_file(argv[1]);

  // on recupere le parametre d'adimensionalisation
  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  // on lance la simulation
  MaxwellProblem1D var;  
  
  if (type_equation == "HELMHOLTZ")
    {
      // simulation harmonique
      RunHarmonicProblem(var, input_file);
    }
  else
    {
      // simulation temporelle
      RunTimeProblem(var, input_file);
    }
  
  return FinalizeMontjoie();
}
