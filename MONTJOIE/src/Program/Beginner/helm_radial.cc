#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Harmonic/MontjoieHarmonic.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/HelmholtzRadial.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#include "Elliptic/Helmholtz/HelmholtzRadial.cxx"
#endif

using namespace Montjoie;

template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file, const string& name_elt)
{
  typedef Montjoie::Real_wp Real_wp;
  typedef Montjoie::Complex_wp Complex_wp;

  MontjoieTimer var_chrono;

  var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var_chrono.Start("GlobalComputation");
  
  var_chrono.SetMessage("ConstructAll", "initialize the computation");
  var_chrono.Start("ConstructAll");
  
  // initialisation of computations
  var.ConstructAll(input_file, name_elt, var_chrono);
  
  var_chrono.Stop("ConstructAll");
  
  // computation of the right hand sides
  var_chrono.SetMessage("ComputationRhs", "compute the right hand side");
  var_chrono.Start("ComputationRhs");

  Vector<Vector<Vector<Complex_wp> > > rhs;
  var.ComputeRightHandSide(rhs, var_chrono);
  
  var_chrono.Stop("ComputationRhs");

  var_chrono.SetMessage("FactorisationMatrix", "compute and factorize finite element matrix");
  var_chrono.SetMessage("ComputationSolution", "solves linear systems");
  
  // loop on the modes
  bool test_loop = true;
  int l = 0;
  Matrix<Complex_wp, General, BandedCol> mat_lu;
  Vector<int> pivot;
  Real_wp norme_max(0);
  int nb_modes_all = 0;
  while (test_loop)
    {
      var_chrono.Start("FactorisationMatrix");
      
      // computation and factorisation of finite element matrix
      var.PerformFactorizationStep(mat_lu, pivot, l);
      
      var_chrono.Stop("FactorisationMatrix");
      
      var_chrono.Start("ComputationSolution");
      
      // loop over m to solve each mode
      Real_wp x_inf = 0;
      if (var.GetPrintLevel() >= 8)
        rhs(0)(0).Write("RHS.dat");
      
      for (int m = 0; m < rhs(l).GetM(); m++)
        if (rhs(l)(m).GetM() > 0)
          {
            var.ComputeSolution(mat_lu, pivot, rhs(l)(m));
            x_inf = max(x_inf, rhs(l)(m).GetNormInf());
	    //DISP(l); DISP(m); DISP(x_inf); rhs(l)(m).Write("sol.dat");
            nb_modes_all++;
          }
      
      var_chrono.Stop("ComputationSolution");
      
      // updating the norm of x
      cout << "Infinite norm of x for mode " << l << " = " << x_inf << endl;
      
      norme_max = max(norme_max, x_inf);
      
      // testing if the computation must be continued
      if (var.number_mode_to_be_computed)
        {
          if (x_inf/norme_max <= var.threshold_mode)
            test_loop = false;
        }
      
      if (l == var.GetLmax())
        test_loop = false;
      
      l++;
    }
  
  cout << "Number of total spherical harmonics = " << nb_modes_all << endl;
  // removing modes greater than l
  int L = l;
  for (int l = L; l <= var.GetLmax(); l++)
    rhs(l).Clear();
  
  // writing the solution
  var_chrono.SetMessage("WriteSolution", "interpolate the solution on grid points");
  var_chrono.SetMessage("InterpolationGrid", "compute the interpolation grid");
  var_chrono.Start("InterpolationGrid");
 
  var.ComputeVarGrid();

  var_chrono.Stop("InterpolationGrid");
  var_chrono.Start("WriteSolution");

  var.WriteDatas(rhs);

  var_chrono.Stop("WriteSolution");
  var_chrono.Stop("GlobalComputation");
  
  var_chrono.DisplayAll();
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

      cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      cout<<"Helmholtz Solver with Lobatto quadrilaterals "<<endl; 
      
      EllipticProblem<HelmholtzEquationRadial> vars;
      RunAll(vars, file_name_data, type_element);
      
      cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
  
}
