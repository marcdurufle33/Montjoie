#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Inverse/MontjoieInverse.hxx"

using namespace Montjoie;

// exact functions for Helmholtz
Real_wp exact_fct0(const R3& x)
{
  //return 1.3;
  //return 1.3 + 0.3*cos(x(0))*sin(x(1));
  
  return 1.0/sqrt(1.0 + 0.5*exp(-0.5*(x(0)*x(0) + x(1)*x(1) + x(2)*x(2))));
  //return 2700.0;
}

Real_wp exact_fct1(const R3& x)
{
  //return 1.3;
  //return 1.3 + 0.3*cos(x(0))*sin(x(1));
  
  return 1.0 - 0.2*exp(-((x(0)-1.0)*(x(0)-1.0) + x(1)*x(1) + x(2)*x(2)));
  //return 2700.0;
}

// exact functions for elastodynamic
VectReal_wp val_exact_param_cte(42);

Real_wp exact_cte0(const R3& x) { return val_exact_param_cte(0); }
Real_wp exact_cte1(const R3& x) { return val_exact_param_cte(1); }
Real_wp exact_cte2(const R3& x) { return val_exact_param_cte(2); }
Real_wp exact_cte3(const R3& x) { return val_exact_param_cte(3); }
Real_wp exact_cte4(const R3& x) { return val_exact_param_cte(4); }
Real_wp exact_cte5(const R3& x) { return val_exact_param_cte(5); }
Real_wp exact_cte6(const R3& x) { return val_exact_param_cte(6); }
Real_wp exact_cte7(const R3& x) { return val_exact_param_cte(7); }
Real_wp exact_cte8(const R3& x) { return val_exact_param_cte(8); }
Real_wp exact_cte9(const R3& x) { return val_exact_param_cte(9); }
Real_wp exact_cte10(const R3& x) { return val_exact_param_cte(10); }
Real_wp exact_cte11(const R3& x) { return val_exact_param_cte(11); }
Real_wp exact_cte12(const R3& x) { return val_exact_param_cte(12); }
Real_wp exact_cte13(const R3& x) { return val_exact_param_cte(13); }
Real_wp exact_cte14(const R3& x) { return val_exact_param_cte(14); }
Real_wp exact_cte15(const R3& x) { return val_exact_param_cte(15); }
Real_wp exact_cte16(const R3& x) { return val_exact_param_cte(16); }
Real_wp exact_cte17(const R3& x) { return val_exact_param_cte(17); }
Real_wp exact_cte18(const R3& x) { return val_exact_param_cte(18); }
Real_wp exact_cte19(const R3& x) { return val_exact_param_cte(19); }
Real_wp exact_cte20(const R3& x) { return val_exact_param_cte(20); }
Real_wp exact_cte21(const R3& x) { return val_exact_param_cte(21); }
Real_wp exact_cte22(const R3& x) { return val_exact_param_cte(22); }
Real_wp exact_cte23(const R3& x) { return val_exact_param_cte(23); }
Real_wp exact_cte24(const R3& x) { return val_exact_param_cte(24); }
Real_wp exact_cte25(const R3& x) { return val_exact_param_cte(25); }
Real_wp exact_cte26(const R3& x) { return val_exact_param_cte(26); }
Real_wp exact_cte27(const R3& x) { return val_exact_param_cte(27); }
Real_wp exact_cte28(const R3& x) { return val_exact_param_cte(28); }
Real_wp exact_cte29(const R3& x) { return val_exact_param_cte(29); }
Real_wp exact_cte30(const R3& x) { return val_exact_param_cte(30); }
Real_wp exact_cte31(const R3& x) { return val_exact_param_cte(31); }
Real_wp exact_cte32(const R3& x) { return val_exact_param_cte(32); }
Real_wp exact_cte33(const R3& x) { return val_exact_param_cte(33); }
Real_wp exact_cte34(const R3& x) { return val_exact_param_cte(34); }
Real_wp exact_cte35(const R3& x) { return val_exact_param_cte(35); }
Real_wp exact_cte36(const R3& x) { return val_exact_param_cte(36); }
Real_wp exact_cte37(const R3& x) { return val_exact_param_cte(37); }
Real_wp exact_cte38(const R3& x) { return val_exact_param_cte(38); }
Real_wp exact_cte39(const R3& x) { return val_exact_param_cte(39); }
Real_wp exact_cte40(const R3& x) { return val_exact_param_cte(40); }
Real_wp exact_cte41(const R3& x) { return val_exact_param_cte(41); }

class InputDataParameterC : public InputDataProblem_Base
{
public:
  void SetInputData(const string& keyword, const VectString& param)
  {
    if (keyword == "ExactValue")
      for (int k = 0; k < min(val_exact_param_cte.GetM(), param.GetM()); k++)
        val_exact_param_cte(k) = to_num<Real_wp>(param(k));
  }
  
};

template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var_laplace,
	    const string& input_file, const string& type_element, 
	    const string& type_equation)
{
  // defining exact parameters
  BoundaryInverseProblem<TypeEquation> var;

  // data file is read
  Vector<string> lines_data_file;
#ifdef SELDON_WITH_MPI
  ReadLinesFile(input_file, lines_data_file, MPI_COMM_WORLD);
#else
  ReadLinesFile(input_file, lines_data_file);
#endif
  
  InputDataParameterC input_c;
  ReadInputFile(lines_data_file, input_c);
  
  // forward solver is initialized
  var.Construct(lines_data_file, type_element, type_equation);

  // exact measurements are computed and stored in ydelta
  VectReal_wp ydelta;
  var.EvaluateDifferently(ydelta);

  // testing derivatives
  /* VectReal_wp c0;
  var.FindInitGuess(c0);
  c0.FillRand(); c0 *= 1e-9;
  for (int i = c0.GetM()/2; i < c0.GetM(); i++)
    c0(i) *= 0.01;
  
  DISP(c0);
  //c0(0) = 1.804; c0(1) = 0.846; c0(2) = 1.68e-3; c0(3) = 1.71e-3;
  
  Real_wp h = 1e-24;
  var.CheckDerivatives(c0, h);
  exit(0); */
  
  // initial guess
  VectReal_wp x_sol;
  var.FindInitGuess(x_sol);
  
  itreg::RegMethodP<Real_wp>* param;
  itreg::RegularizationMethod<Real_wp>* iterative_solver;
  
  // regularization method is constructed
  itreg::RegularizationMethod<Real_wp>::GetNewSolver(var.GetRegularizationParameter(), var,
						     param, iterative_solver);
    
  param->print_level = 1;
  
  // we call the iterative process to find the solution
  // of the considered inverse problem
  Real_wp delta(0);
  iterative_solver->Solve(ydelta, delta, x_sol);
    
  // solutiona are written on the disk
  var.WriteDatas(x_sol);

  delete param;
  delete iterative_solver;
}

int main(int argc, char **argv) 
{

  InitMontjoie(argc, argv);

  HelmholtzEquation_Base<Complex_wp, Dimension2>::store_dfjm1 = true;
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
            
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      cout<<"Helmholtz Solver with Lobatto quadrilaterals "<<endl; 
      
      if (!type_equation.compare("HELMHOLTZ"))
	{          
	  EllipticProblem<HelmholtzEquation<Dimension3> > vars;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[0] = &exact_fct0;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[1] = &exact_fct1;
	  RunAll(vars, file_name_data, type_element, type_equation);
	}
      else
	{
	  EllipticProblem<HarmonicElasticEquation<Dimension3> > vars;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[0] = &exact_cte0;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[1] = &exact_cte1;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[2] = &exact_cte2;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[3] = &exact_cte3;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[4] = &exact_cte4;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[5] = &exact_cte5;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[6] = &exact_cte6;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[7] = &exact_cte7;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[8] = &exact_cte8;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[9] = &exact_cte9;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[10] = &exact_cte10;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[11] = &exact_cte11;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[12] = &exact_cte12;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[13] = &exact_cte13;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[14] = &exact_cte14;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[15] = &exact_cte15;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[16] = &exact_cte16;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[17] = &exact_cte17;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[18] = &exact_cte18;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[19] = &exact_cte19;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[20] = &exact_cte20;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[21] = &exact_cte21;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[22] = &exact_cte22;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[23] = &exact_cte23;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[24] = &exact_cte24;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[25] = &exact_cte25;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[26] = &exact_cte26;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[27] = &exact_cte27;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[28] = &exact_cte28;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[29] = &exact_cte29;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[30] = &exact_cte30;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[31] = &exact_cte31;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[32] = &exact_cte32;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[33] = &exact_cte33;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[34] = &exact_cte34;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[35] = &exact_cte35;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[36] = &exact_cte36;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[37] = &exact_cte37;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[38] = &exact_cte38;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[39] = &exact_cte39;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[40] = &exact_cte40;
          BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[41] = &exact_cte41;

	  RunAll(vars, file_name_data, type_element, type_equation);
	}	
      
      cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  cout<<"End of the program"<<endl; 
  
  return FinalizeMontjoie();
  
}
