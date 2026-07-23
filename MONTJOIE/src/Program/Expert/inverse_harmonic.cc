#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Inverse/MontjoieInverse.hxx"

using namespace Montjoie;

// exact functions for Helmholtz
Real_wp exact_fct0(const R2& x)
{
  //return 1.3;
  //return 1.3 + 0.3*cos(x(0))*sin(x(1));
  
  //return x(1)/5;
  return 1.0/sqrt(1.0 + 0.5*exp(-0.5*(x(0)*x(0) + x(1)*x(1))));
  //return 2700.0;
}

Real_wp exact_fct1(const R2& x)
{
  //return 1.3;
  //return 1.3 + 0.3*cos(x(0))*sin(x(1));
  
  //return -x(0)/5;
  return 1.0 - 0.2*exp(-((x(0)-1.0)*(x(0)-1.0) + x(1)*x(1)));
  //return 2700.0;
}

Real_wp exact_fct2(const R2& x)
{
  return 2.0;
}

Real_wp exact_fct3(const R2& x)
{
  return 2.0;
}

Real_wp exact_fct4(const R2& x)
{
  return 2.0;
}

Real_wp exact_fct5(const R2& x)
{
  return 2.0;
}

// 3-D Helmholtz
Real_wp exact_fct0_3D(const R3& x)
{
  return 1.0;
}

Real_wp exact_fct1_3D(const R3& x)
{
  return 1.0;
}

Real_wp exact_fct2_3D(const R3& x)
{
  return 1.0;
}

Real_wp exact_fct3_3D(const R3& x)
{
  return 1.0;
}

Real_wp exact_fct4_3D(const R3& x)
{
  return 1.0;
}

Real_wp exact_fct5_3D(const R3& x)
{
  return 1.0;
}


// exact functions for elastodynamic
VectReal_wp val_exact_param_cte(16);

Real_wp exact_cte0(const R2& x) { return val_exact_param_cte(0); }
Real_wp exact_cte1(const R2& x) { return val_exact_param_cte(1); }
Real_wp exact_cte2(const R2& x) { return val_exact_param_cte(2); }
Real_wp exact_cte3(const R2& x) { return val_exact_param_cte(3); }
Real_wp exact_cte4(const R2& x) { return val_exact_param_cte(4); }
Real_wp exact_cte5(const R2& x) { return val_exact_param_cte(5); }
Real_wp exact_cte6(const R2& x) { return val_exact_param_cte(6); }
Real_wp exact_cte7(const R2& x) { return val_exact_param_cte(7); }
Real_wp exact_cte8(const R2& x) { return val_exact_param_cte(8); }
Real_wp exact_cte9(const R2& x) { return val_exact_param_cte(9); }
Real_wp exact_cte10(const R2& x) { return val_exact_param_cte(10); }
Real_wp exact_cte11(const R2& x) { return val_exact_param_cte(11); }
Real_wp exact_cte12(const R2& x) { return val_exact_param_cte(12); }
Real_wp exact_cte13(const R2& x) { return val_exact_param_cte(13); }

// 3-D case
Real_wp exact_cte0_3D(const R3& x) { return val_exact_param_cte(0); }
Real_wp exact_cte1_3D(const R3& x) { return val_exact_param_cte(1); }
Real_wp exact_cte2_3D(const R3& x) { return val_exact_param_cte(2); }
Real_wp exact_cte3_3D(const R3& x) { return val_exact_param_cte(3); }
Real_wp exact_cte4_3D(const R3& x) { return val_exact_param_cte(4); }
Real_wp exact_cte5_3D(const R3& x) { return val_exact_param_cte(5); }
Real_wp exact_cte6_3D(const R3& x) { return val_exact_param_cte(6); }
Real_wp exact_cte7_3D(const R3& x) { return val_exact_param_cte(7); }
Real_wp exact_cte8_3D(const R3& x) { return val_exact_param_cte(8); }
Real_wp exact_cte9_3D(const R3& x) { return val_exact_param_cte(9); }
Real_wp exact_cte10_3D(const R3& x) { return val_exact_param_cte(10); }
Real_wp exact_cte11_3D(const R3& x) { return val_exact_param_cte(11); }
Real_wp exact_cte12_3D(const R3& x) { return val_exact_param_cte(12); }
Real_wp exact_cte13_3D(const R3& x) { return val_exact_param_cte(13); }

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
  
  // additional exact values are provided
  var.ProvideExactValues(val_exact_param_cte);
  
  // exact measurements are computed and stored in ydelta
  VectReal_wp ydelta;
  var.EvaluateDifferently(ydelta);
  
  // testing derivatives 
  /* VectReal_wp c0;
  var.FindInitGuess(c0);
  //c0.FillRand(); c0 *= 1e-10;
  //for (int i = c0.GetM()/2; i < c0.GetM(); i++)
    //c0(i) *= 0.01;
  
  DISP(c0);
  //c0(0) = 1.804; c0(1) = 0.846; c0(2) = 1.68e-3; c0(3) = 1.71e-3;
  
  Real_wp h = 1e-16;
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
  int rank_proc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  if (rank_proc > 0)
    param->print_level = 0;
  
  iterative_solver->SetCommunicator(MPI_COMM_WORLD);
  
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

      int dim_N = 3;
      if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
	dim_N = 2;

      if (!type_equation.compare("HELMHOLTZ"))
	{          
	  if (dim_N == 2)
	    {
	      EllipticProblem<HelmholtzEquation<Dimension2> > vars;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[0] = &exact_fct0;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[1] = &exact_fct1;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[2] = &exact_fct2;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[3] = &exact_fct3;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[4] = &exact_fct4;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[5] = &exact_fct5;
	      RunAll(vars, file_name_data, type_element, type_equation);
	    }
	  else
	    {
	      EllipticProblem<HelmholtzEquation<Dimension3> > vars;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[0] = &exact_fct0_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[1] = &exact_fct1_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[2] = &exact_fct2_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[3] = &exact_fct3_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[4] = &exact_fct4_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[5] = &exact_fct5_3D;
	      RunAll(vars, file_name_data, type_element, type_equation);
	    }
	}
      else if (!type_equation.compare("HELMHOLTZ_AXI"))
	{
	  EllipticProblem<HelmholtzEquationAxi> vars;
          BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[0] = &exact_fct0;
          BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[1] = &exact_fct1;
          BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[2] = &exact_fct2;
          BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[3] = &exact_fct3;
          BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[4] = &exact_fct4;
          BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[5] = &exact_fct5;
	  RunAll(vars, file_name_data, type_element, type_equation);
	}
      else
	{
	  if (dim_N == 2)
	    {	      
	      EllipticProblem<HarmonicElasticEquation<Dimension2> > vars;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[0] = &exact_cte0;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[1] = &exact_cte1;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[2] = &exact_cte2;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[3] = &exact_cte3;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[4] = &exact_cte4;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[5] = &exact_cte5;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[6] = &exact_cte6;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[7] = &exact_cte7;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[8] = &exact_cte8;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[9] = &exact_cte9;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[10] = &exact_cte10;
	      BoundaryInverseProblem_Base<Dimension2>::exact_fct_to_find[11] = &exact_cte11;
	      RunAll(vars, file_name_data, type_element, type_equation);
	    }
	  else
	    {	      
	      EllipticProblem<HarmonicElasticEquation<Dimension3> > vars;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[0] = &exact_cte0_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[1] = &exact_cte1_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[2] = &exact_cte2_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[3] = &exact_cte3_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[4] = &exact_cte4_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[5] = &exact_cte5_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[6] = &exact_cte6_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[7] = &exact_cte7_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[8] = &exact_cte8_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[9] = &exact_cte9_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[10] = &exact_cte10_3D;
	      BoundaryInverseProblem_Base<Dimension3>::exact_fct_to_find[11] = &exact_cte11_3D;
	      RunAll(vars, file_name_data, type_element, type_equation);
	    }	    
	}
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
  
}
