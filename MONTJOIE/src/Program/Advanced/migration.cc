#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_HP_H1
#define MONTJOIE_WITH_NODAL_DG

#define SELDON_WITH_PRECONDITIONING

#define MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
#define MONTJOIE_WITH_THIN_SLOT_MODEL
#define MONTJOIE_WITH_TRANSMISSION

//#define MONTJOIE_WITH_WIRES

#include "Inverse/MontjoieInverse.hxx"

using namespace Montjoie;


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var_simu,
	    const string& input_file, const string& type_element,
	    const string& type_equation)
{
  // input file is read
  Vector<string> lines_data_file;
#ifdef SELDON_WITH_MPI
  ReadLinesFile(input_file, lines_data_file, var_simu.comm_group_mode);
#else
  ReadLinesFile(input_file, lines_data_file);
#endif
  
  EllipticProblem<TypeEquation> var_exp;
  VarMigration<TypeEquation> var_migration;
  
  ReadInputFile(lines_data_file, var_migration);

  // meshes, finite elements and geometry quantities are constructed
  All_LinearSolver* solver_exp, *solver_simu;
  var_exp.ConstructAll(var_migration.input_file_exp, type_element, type_equation, solver_exp);
  var_simu.ConstructAll(var_migration.input_file_simu, type_element, type_equation, solver_simu);
  
  // quadrature points on the surface where du/dn is imposed
  var_migration.ComputeQuadrature(var_simu, var_exp);

  // experimental data is generated
  Matrix<Complex_wp, General, ColMajor> data_exp;
  var_migration.ComputeExperiment(var_exp, solver_exp, data_exp);

  // simulated data is generated (with the same right hand side)
  Matrix<Complex_wp, General, ColMajor> data_simu, sol_simu;
  var_migration.ComputeSimulation(var_simu, solver_simu, data_simu, sol_simu);
  
  // the residual is back-propagated
  Matrix<Complex_wp, General, ColMajor>& residu = data_simu;
  Add(-Complex_wp(1, 0), data_exp, residu);
  Matrix<Complex_wp, General, ColMajor> sol_residu;
  var_migration.ComputeAdjoint(var_simu, solver_simu, residu, sol_residu);

  // factorisation is cleared for var_simu
  solver_simu->ClearFactorization();

  // image is constructed
  VectComplex_wp image(var_simu.GetNbDof());
  image.Zero();
  for (int k = 0; k < sol_residu.GetN(); k++)
    {
      for (int j = 0; j < var_simu.GetNbDof(); j++)
	image(j) += sol_residu(j, k)*sol_simu(j, k);
    }
  
  var_simu.WriteOutputFile(image, 1);
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      
      if (!type_equation.compare("HELMHOLTZ"))
        {
          EllipticProblem<HelmholtzEquation<Dimension2> > Vars;
          if (Vars.GetRankProcMode() == 0)
            cout<<"Helmholtz Solver with continuous formulation"<<endl; 
	  
	  RunAll(Vars, file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
      else
        {
          EllipticProblem<HarmonicElasticEquation<Dimension2> > Vars;	      
          if (Vars.GetRankProcMode() == 0)
            cout<<"Time-harmonic elastic solver with continuous formulation "<<endl; 
          
	  RunAll(Vars, file_name_data, type_element, type_equation);
          if (Vars.GetRankProcMode() == 0)
            cout<<" we destroy the variables "<<endl;
        }
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"migration.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
