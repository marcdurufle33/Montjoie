#define MONTJOIE_WITH_TWO_DIM
//#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

#include "Elliptic/Helmholtz/AxiSymHelmholtz.hxx"
#include "Elliptic/Helmholtz/AxiSymHelmholtzInline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/AxiSymHelmholtz.cxx"
#endif

#include "Helio/CrossCovariance.cxx"

using namespace Montjoie;

// compute cross-covariance
template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& name_elt,
            const string& name_equation, const string& name_file)
{  
  // parametres associe au pb multi-freq
  CrossCovarianceParameter param_cross;

  // pb multi-frequences
  CrossCovarianceProblem var_cross(param_cross);  
  
  
  // Timers initialization, input file is read
  MontjoieTimer& chrono = var_cross.chrono;

  Vector<string> lines_data_file;
  
#ifdef SELDON_WITH_MPI
  ReadLinesFile(name_file, lines_data_file, var_cross.GetFrequencyCommunicator());
  chrono.SetCommunicator(var_cross.GetFrequencyCommunicator());
#else
  ReadLinesFile(name_file, lines_data_file);
#endif

  // the sources are read
  //ReadLinesFile(listSourcesFile, param_cross.srcLines); 
  //param_cross.InitSource(lines_data_file);

  //############################################
  // Initialization du probleme multi-frequences
  
  ReadInputFile(lines_data_file, var_cross);
  
  // on distribue les frequences sur les differents procs
  var_cross.DistributeFrequency(var_cross.coef_freq);
  
#ifdef SELDON_WITH_MPI
  var.comm_group_mode = var_cross.GetLocalCommunicator();
#endif
  
  chrono.SetMessage("GlobalComputation"  ,"perform the overall simulation");
  chrono.SetMessage("ConstructAll"       ,"initialize the computation");
  chrono.SetMessage("ComputationRhs"     ,"build the right hand side");
  chrono.SetMessage("FactorizationMatrix","compute and factorize finite element matrix");
  chrono.SetMessage("ComputationSolution","solve linear systems");
  chrono.SetMessage("WriteSolution"      ,"interpolate the solution on output points and write");
  
  // ================================================================================
 
  chrono.Start("GlobalComputation");
  chrono.Start("ConstructAll");
  
  // the input file is treated, the mesh constructed, solver constructed  
  var.ConstructAll(name_file, name_elt, lines_data_file);
  
  All_LinearSolver glob_solver(var);
  ReadInputFile(lines_data_file, glob_solver);

  var.ComputeMassMatrix(true);
  
  chrono.Stop("ConstructAll");
  
  // on initialise la matrice de covariance ou kernel
  var_cross.InitCovariance(var);
  
  // on force l'ecriture de la solution pour chaque mode
  var.ForceStorageModes();
  
  // on resout toutes les frequences et tous les modes
  var_cross.SolveAllFreq(var, &glob_solver, param_cross);
  
  // on sort la matrice de covariance
  if (var_cross.compute_covariance)
    var_cross.WriteCovariance();

  if (var_cross.compute_kernel)
    var_cross.WriteKernel();

  chrono.Stop("GlobalComputation");
  
  chrono.DisplayAll();
  
}

// ======================================================================================================

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
      
      // Harmonic problem
      EllipticProblem<HelmholtzEquationAxi> Vars;
      RunAll(Vars, typeElement, typeEquation, inputFileName);
      
      
      if (Vars.GetPrintLevel()>0)
        {
          cout << "End of the program" << endl; 
          cout << "===========================================================================" << endl;
        }      
    }
  else
    {
      cout << "This code needs a data file and a file of list of sources in argument" << endl;
      cout << "helmholtz_axi.x configurationFileName srcListFileName" << endl;
      cout << "is a good syntax" << endl;
    }
  
  return FinalizeMontjoie();
}
