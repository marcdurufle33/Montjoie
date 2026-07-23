
#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

#include "Elliptic/Thermo/ViscoThermalEquation.cxx"
#include "Harmonic/VarHarmonic.cxx"
#include "Harmonic/GenericEquation.cxx"
#include "Computation/ElementaryMatrixH1.cxx"

using namespace Montjoie;

class MultiFrequencyParamClass : public MultiFrequencyParameter
{
  void SetSourceFrequency(int k, int kg, int m, int mg,
                          const VarComputationProblem& var,
                          Matrix<Complex_wp, General, ColMajor>& sol)
  {
    VarSourceProblem_Base& var_source
      = const_cast<VarSourceProblem_Base&>(dynamic_cast<const VarSourceProblem_Base&>(var));
    
    Vector<Complex_wp> x_sol;
    var_source.ComputeRightHandSide(x_sol);
    sol.Resize(x_sol.GetM(), 1);
    SetCol(x_sol, 0, sol);
  }
    
  
  void GiveSolutionFrequency(int k, int kg, int m, int mg,
                             const VarComputationProblem& var,
                             Matrix<Complex_wp, General, ColMajor>& sol_exp)
  {
    const VarOutputProblem_Base& var_output = dynamic_cast<const VarOutputProblem_Base&>(var);
    Vector<Complex_wp> x_sol(sol_exp.GetM());
    GetCol(sol_exp, 0, x_sol);
    var_output.WriteOutputFile(x_sol, "SolOm" + to_str(kg) + "_G");
  }
  
};


// compute multi-frequency solution
template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& name_elt,
            const string& name_equation, const string& input_file)
{  
  // the sources are read
  Vector<string> lines_data_file;
  
  MultiFrequencyProblem multi_freq;

#ifdef SELDON_WITH_MPI
  ReadLinesFile(input_file, lines_data_file, MPI_COMM_WORLD);
#else
  ReadLinesFile(input_file, lines_data_file);
#endif

  //############################################
  // Initialization du probleme multi-frequences
  
  MultiFrequencyParamClass param_freq;
  ReadInputFile(lines_data_file, multi_freq);
  
  // on distribue les frequences sur les differents procs
  VectReal_wp coef;
  multi_freq.DistributeFrequency(coef);
  
#ifdef SELDON_WITH_MPI
  var.comm_group_mode = multi_freq.GetLocalCommunicator();
  int size_loc;
  MPI_Comm_size(var.comm_group_mode, &size_loc); DISP(size_loc);
#endif
  
  All_LinearSolver* glob_solver;

  // the input file is treated, the mesh constructed, solver constructed  
  var.ConstructAll(lines_data_file, input_file, name_elt, name_equation, glob_solver);
  
  // on resout toutes les frequences
  multi_freq.SolveAllFreq(var, glob_solver, param_freq);
  
}

int main(int argc, char **argv)
{  
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {      
      string inputFileName;
      inputFileName = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string typeElement,typeEquation;
      getElement_Equation(inputFileName, typeElement, typeEquation);
      
      int dim_N = 3;
      if ((typeElement.find("TRIANGLE") == 0) || (typeElement.find("QUADRANGLE") == 0))
	dim_N = 2;
      
      // Harmonic problem
      if (typeEquation == "HELMHOLTZ")
        {
          if (dim_N == 2)
            {
              EllipticProblem<HelmholtzEquation<Dimension2> > Vars;
              RunAll(Vars, typeElement, typeEquation, inputFileName);      
            }
          else
            {
              EllipticProblem<HelmholtzEquation<Dimension3> > Vars;
              RunAll(Vars, typeElement, typeEquation, inputFileName);      
            }
        }
      else if (typeEquation == "HELMHOLTZ_AXI")
        {
          EllipticProblem<HelmholtzEquationAxi> Vars;
          RunAll(Vars, typeElement, typeEquation, inputFileName);      
        }
      else if (typeEquation == "MAXWELL")
        {
          EllipticProblem<HarmonicMaxwellEquation_3D> Vars;
          RunAll(Vars, typeElement, typeEquation, inputFileName);      
        }
      else if (typeEquation == "VISCO_THERMAL")
        {
          if (dim_N == 2)
            {
              EllipticProblem<ViscoThermalEquation<Dimension2, Complex_wp> > Vars;
              RunAll(Vars, typeElement, typeEquation, inputFileName);      
            }
          else
            {
              EllipticProblem<ViscoThermalEquation<Dimension3, Complex_wp> > Vars;
              RunAll(Vars, typeElement, typeEquation, inputFileName);
            }
        }
      else
        {
          cout << "Unknown equation " << typeEquation << endl;
          abort();
        }

    }
  else
    {
      cout << "This code needs a data file in argument" << endl;
      cout << "multi_freq.x name" << endl;
      cout << "is a good syntax" << endl;
    }
  
  return FinalizeMontjoie();
}
