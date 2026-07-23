#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Harmonic/MontjoieHarmonic.hxx"
#include "Regularization/Itreg.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#endif

using namespace Montjoie;

class InverseProblem1D : public InputDataProblem_Base, public itreg::ForwardOperator<Real_wp>
{
protected:
  //! exact definition of the default
  VectReal_wp parameter_exact;

  //! object solving Helmholtz equation
  EllipticProblem<HelmholtzEquation1D> var_helm;

  //! matrices storing LU factors
  Vector<Matrix<Complex_wp, General, BandedCol> > mat_lu;
  Vector<Vector<int> > pivot;

  //! stored solutions
  Vector<VectComplex_wp> sol_u;

  VectReal_wp parameter_cn;
  
  //! explored frequencies
  VectReal_wp omega;
  
  //! parameters associated with the regularization method
  Vector<string> regularization_parameter;
  string type_element;
  Vector<string> initial_mesh_parameter;
  Real_wp rho_media1, rho_media2;
  
public:
  InverseProblem1D()
  {
    parameter_exact.Reallocate(3);
    parameter_exact(0) = 7.0;
    parameter_exact(1) = 0.04;
    parameter_exact(2) = 3.0;
    
    rho_media1 = 1.0;
    rho_media2 = 1.0;
  }


  //! returns parameters associated with the iterative regularization method
  const Vector<string>& GetRegularizationParameter() const
  {
    return regularization_parameter;
  }
  

  //! reads a line of the Montjoie data file
  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "RegularizationMethod")
      regularization_parameter = param;    
    else if (keyword == "RangeFrequency")
      {
	Real_wp omega_min = 2*pi_wp*to_num<Real_wp>(param(0));
	Real_wp omega_max = 2*pi_wp*to_num<Real_wp>(param(1));
	int nb_omega = to_num<int>(param(2));
	Linspace(omega_min, omega_max, nb_omega, omega);
      }
    else if (keyword == "LayeredBackground")
      {
	Real_wp x0 = to_num<Real_wp>(initial_mesh_parameter(1));
	Real_wp xN = to_num<Real_wp>(initial_mesh_parameter(2));
	Real_wp dx = to_num<Real_wp>(initial_mesh_parameter(4));
	
	Real_wp pos_init = to_num<Real_wp>(param(0));
	Real_wp size1 = to_num<Real_wp>(param(1));
	Real_wp size2 = to_num<Real_wp>(param(2));
	int nb_cells = to_num<Real_wp>(param(3));
	rho_media1 = to_num<Real_wp>(param(4));
	rho_media2 = to_num<Real_wp>(param(5));
	initial_mesh_parameter.Reallocate(2*nb_cells+6);
	initial_mesh_parameter(0) = "LAYERED";
	initial_mesh_parameter(1) = to_str(x0);
	Real_wp pos = pos_init;
	int nb = 2;
	for (int i = 0; i <= 2*nb_cells; i++)
	  {
	    initial_mesh_parameter(nb++) = to_str(pos);
	    if (i%2 == 0)
	      pos += size1;
	    else
	      pos += size2;
	  }

	initial_mesh_parameter(nb++) = to_str(xN);
	initial_mesh_parameter(nb++) = "AUTO";
	initial_mesh_parameter(nb++) = to_str(dx);
      }

  }
  

  //! returns true if the exact solution has been stored
  bool ExactSolutionKnown() const
  {
    return true;
  }

  
  //! returns the distance to the exact solution
  Real_wp DistanceToExactSolution(const VectReal_wp& x)
  {
    return Norm2(x - parameter_exact) / Norm2(x);
  }
  
  
  //! initialization of montjoie object var_helm
  void Construct(const string& parameter_file_montjoie, const string& name_element)
  {
    // the input file is read
    ReadInputFile(parameter_file_montjoie, var_helm);
    type_element = name_element;
    initial_mesh_parameter = var_helm.GetMeshData(0);

    ReadInputFile(parameter_file_montjoie, *this);
  }

  
  //! initial guess
  void FindInitGuess(VectReal_wp& c)
  {
    c.Reallocate(3);
    c(0) = 7.0;
    c(1) = 0.2;
    c(2) = 1.0;
  }

  
  //! computes the solution for a given set of parameters cn
  void Init(const VectReal_wp& cn)
  { 
    parameter_cn = cn;
    //Real_wp pos_target = cn(0);
    Real_wp pos_target = 7.0;
    Real_wp size_target = abs(cn(1));
    Real_wp rho = abs(cn(2));
    Vector<string> data = initial_mesh_parameter;
    Real_wp ref_defaut = -1;
    int nb_ref = 0;
    if (size_target > 0)
      {
	if (data(0) == "LAYERED")
	  {
	    int position_word = -1;
	    // we look for the first parameter, with begins by a letter
	    // This parameter is the description of the space step strategy
	    for (int i = 1; i < data.GetM(); i++)
	      if (isalpha(data(i)[0]))
		{
		  position_word = i;
		  break;
		}
	    
	    if (position_word > 2)
	      {
		int nb = position_word+1;            
		if (data(position_word) != "AUTO")
		  {
		    cout << "Not implemented " << endl;
		    abort();
		  }
		
		string dx = data(nb);
		VectReal_wp pos(position_word-1);
		for (int i = 1; i < position_word; i++)
		  pos(i-1) = to_num<Real_wp>(data(i));
		
		nb_ref = pos.GetM()+1;
		int imin = 0;
		while ((imin < pos.GetM()) && (pos(imin) < pos_target-R2::threshold))
		  imin++;
		
		int imax = pos.GetM()-1;
		while ((imax >= 0) && (pos(imax) > pos_target+size_target+R2::threshold))
		  imax--;
		
		data.Reallocate(5 + imin + (pos.GetM()-1-imax));
		data(0) = "LAYERED";
		for (int i = 0; i < imin; i++)
		  data(i+1) = to_str(pos(i));
		
		ref_defaut = imin+1;
		data(imin+1) = to_str(pos_target);
		data(imin+2) = to_str(pos_target+size_target);
		nb = imin+3;
		for (int i = imax+1; i < pos.GetM(); i++)
		  data(nb++) = to_str(pos(i));
		
		data(nb++) = "AUTO";
		data(nb++) = dx;
	      }
	  }
	else
	  {
	    cout << "Not implemented for non-layered mesh" << endl;
	    abort();
	  }
      }

    var_helm.GetMeshData(0) = data;
    
    for (int k = 0; k <= nb_ref; k++)
      {
	if (k == ref_defaut)
	  var_helm.ref_rho(k).SetConstant(rho);
	else
	  {
	    if ((k == 0) || (k == nb_ref))
	      var_helm.ref_rho(k).SetConstant(1.0);
	    else if (k%2 == 1)
	      var_helm.ref_rho(k).SetConstant(rho_media1);
	    else
	      var_helm.ref_rho(k).SetConstant(rho_media2);
	  }
      }
    
    // construction of mesh and finite element
    var_helm.ComputeMeshAndFiniteElement(type_element);
    
    // physical coefficients
    var_helm.ComputePhysicalCoefficients();

    GlobalGenericMatrix<Complex_wp> nat_mat;
    int N = var_helm.GetNbDof();
    int r = var_helm.mesh.GetOrder();    

    mat_lu.Reallocate(omega.GetM());
    pivot.Reallocate(omega.GetM());
    sol_u.Reallocate(omega.GetM());
    for (int k = 0; k < omega.GetM(); k++)
      {
	var_helm.SetOmega(omega(k));
	var_helm.ComputeMassMatrix();
	
	mat_lu(k).Reallocate(N, N, r, r);
	mat_lu(k).Zero();
	
	// Finite element matrix is computed and factorized
	var_helm.AddMatrixFEM(mat_lu(k), nat_mat);
	var_helm.AddBoundaryTerms(mat_lu(k), nat_mat);
	Seldon::GetLU(mat_lu(k), pivot(k));

	var_helm.ComputeRightHandSide(sol_u(k));
	Seldon::SolveLU(mat_lu(k), pivot(k), sol_u(k));
      }
  }


  //! computes the observable vector
  void Evaluate(VectReal_wp& y)
  {
    y.Reallocate(2*omega.GetM());
    for (int k = 0; k < omega.GetM(); k++)
      {
	y(2*k) = realpart(sol_u(k)(0));
	y(2*k+1) = imagpart(sol_u(k)(0));
      }
  }


  //! Applies jacobian matrix of F to vector h, res = DF(cn) h
  void Derivative(const VectReal_wp& h, VectReal_wp& res)
  {
    Real_wp eps = 1e-5;
    VectReal_wp cp(parameter_cn), cm(parameter_cn), ym, yp;
    cp += eps*h;
    cm -= eps*h;
        
    Init(cm);
    Evaluate(ym);
    
    Init(cp);
    Evaluate(yp);

    res = (yp - ym) / (2.0*eps);
    
    parameter_cn = cm + eps*h;
  }


  //! Applies adjoint of jacobian matrix of F to vector G, res = DF*(cn) h
  void Adjoint(const VectReal_wp& G, VectReal_wp& res)
  {
    abort();
  }
  

  //! Computes the exact measures y, parameter c will
  //! solve the non-linear least-squares problem F(c) = y
  void EvaluateDifferently(VectReal_wp& y)
  {
    Init(parameter_exact);
    Evaluate(y);

    for (int k = 0; k < omega.GetM(); k++)
      {
	var_helm.GetOutputParameters().SetTotalFieldFile("SolOmega"+to_str(k)+".dat");
	var_helm.WriteDatas(sol_u(k));
      }
  }
  
};


void RunAll(const string& input_file, const string& type_element)
{
  // forward solver is initialized
  InverseProblem1D var;
  var.Construct(input_file, type_element);

  // testing derivatives
  /* VectReal_wp c0;
  var.FindInitGuess(c0);
  c0.FillRand(); c0 *= 1e-9;
  
  Real_wp h = 1e-6;
  var.CheckDerivatives(c0, h); */

  // exact measurements are computed and stored in ydelta
  VectReal_wp ydelta;
  var.EvaluateDifferently(ydelta);
  
  // initial guess
  VectReal_wp x_sol;
  var.FindInitGuess(x_sol);
  
  itreg::RegMethodP<Real_wp>* param;
  itreg::RegularizationMethod<Real_wp>* iterative_solver;
  
  // regularization method is constructed
  itreg::RegularizationMethod<Real_wp>::GetNewSolver(var.GetRegularizationParameter(), var,
						     param, iterative_solver);
    
  param->print_level = 1;
  
  glob_chrono.Start(VirtualTimer::ALL);

  // we call the iterative process to find the solution
  // of the considered inverse problem
  Real_wp delta(0);
  iterative_solver->Solve(ydelta, delta, x_sol);
  
  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Temps de resolution = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << "s" << endl;
  
  // solutiona are written on the disk
  DISP(x_sol);
  //var.WriteDatas(x_sol);

  delete param;
  delete iterative_solver;
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
      
      RunAll(file_name_data, type_element);
      
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
