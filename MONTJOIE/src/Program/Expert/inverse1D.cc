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


Real_wp pol_sin10x(Real_wp x)
{
  Real_wp aux = sin(pi_wp*x);
  Real_wp res = Real_wp(1) + aux*aux*aux*sin(10*x);
  return res;
}

Real_wp const1(Real_wp x)
{
  return Real_wp(1);
}

class InverseProblem1D : public InputDataProblem_Base, public itreg::ForwardOperator<Real_wp>
{
protected:
  //! pointer to the source fonction
  Real_wp (*f_src)(Real_wp);
  
  //! pointer to the exact c
  Real_wp (*c_exact)(Real_wp);
  
  //! object for solving the Laplacian
  EllipticProblem<LaplaceEquation1D> var_laplace;
  
  //! matrix storing LU factors
  Matrix<Real_wp, General, BandedCol> mat_lu;
  // and pivot;
  Vector<int> pivot;
  
  //! stored solution
  VectReal_wp sol_u;
  //! Dirichlet conditions
  Real_wp ua, ub;

  //! exact values of C stored
  VectReal_wp C_store_exact;
  
  //! object discretizing c
  VariableParameter1D_Base<Real_wp>* c_discrete;
  
  //! points where c are evaluated
  VectReal_wp grid_evalC;

  //! Collocation matrix : CollMatrix(i, j) = phi_j(xi_i)
  Matrix<Real_wp, General, ArrayRowSparse> CollMatrix;
  
  // mass matrix
  VectReal_wp mat_mass;
  
  // quadrature points where c needs to be evaluated
  VectReal_wp points_quad;

  //! parameters associated with the regularization method
  Vector<string> regularization_parameter;
  Vector<string> approx_c_parameter;
  
public:
  //! Default constructor
  InverseProblem1D()
  {
    ua = Real_wp(1);
    ub = Real_wp(2);
    f_src = &const1;
    c_exact = &pol_sin10x;
  }

  //! returns parameters associated with the iterative regularization method
  const Vector<string>& GetRegularizationParameter() const
  {
    return regularization_parameter;
  }
  
  //! reads a line of the Montjoie data file
  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "ApproximationC")
      approx_c_parameter = param;
    else if (keyword == "RegularizationMethod")
      regularization_parameter = param;
      
  }

  //! returns true if the exact solution has been stored
  bool ExactSolutionKnown() const
  {
    return true;
  }

  
  //! returns the distance to the exact solution
  Real_wp DistanceToExactSolution(const VectReal_wp& x)
  {
    return Norm2(x - C_store_exact) / Norm2(x);
  }


  //! writes the parameter
  void WriteDatas(const VectReal_wp& cn)
  {
    // defining an output grid
    int nb_pts_output = 1000;
    VectReal_wp x_output;
    Real_wp a = var_laplace.GetXmin(), b = var_laplace.GetXmax();    
    Linspace(a, b, nb_pts_output, x_output);

    ofstream file_out("CnApprox.dat");
    c_discrete->Init(cn); file_out.precision(15);
    for (int i = 0; i < nb_pts_output; i++)
      {
	Real_wp cn_ex = (*c_exact)(x_output(i));
	Real_wp cn_app = c_discrete->Evaluate(x_output(i));
	file_out << x_output(i) << " " << cn_app << " " << cn_ex << '\n';
      }

    file_out.close();
  }
  
  //! initialization of montjoie object var_laplace
  void Construct(const string& parameter_file_montjoie, const string& name_element)
  {
    // the input file is read
    ReadInputFile(parameter_file_montjoie, var_laplace);
    ReadInputFile(parameter_file_montjoie, *this);
    
    // construction of mesh and finite element
    var_laplace.ComputeMeshAndFiniteElement(name_element);
    
    // mass_matrix is diagonal 
    var_laplace.ComputePhysicalCoefficients();
    var_laplace.ComputeMassMatrix();
    
    var_laplace.GetMassMatrix(mat_mass);
    
    // the discretization for c is constructed
    Real_wp a = var_laplace.GetXmin(), b = var_laplace.GetXmax();    
    c_discrete = VariableParameter1D_Base<Real_wp>::GetNewApproximation(a, b, approx_c_parameter);
    
    // the exact value of c (projection in the discrete space) is computed
    VectReal_wp evalC_exact;
    c_discrete->ProjectFunction(*c_exact, evalC_exact);
    C_store_exact = evalC_exact;
    
    // we retrieve quadrature points
    // since c should be evaluated on these points
    // in order to obtain the finite element matrix
    var_laplace.ComputeQuadraturePoints(points_quad);
    
    // computation of CollMatrix that is needed
    // to compute the interpolation of c on the dof points
    int nb_pts_quad_inside = var_laplace.GetNbDof() - 2;
    int nb_basis_fct = c_discrete->GetNbDof();
    CollMatrix.Reallocate(nb_pts_quad_inside, nb_basis_fct);    
    VectReal_wp vec_phi(nb_basis_fct);    
    // here CollMatrix(i, j) = phi_j(points_dof(i))
    for (int i = 0; i < nb_pts_quad_inside; i++)
      {
	c_discrete->ComputeValuesPhiRef(var_laplace.GetCoordinateDof(i+1), vec_phi);
	
	for (int j = 0; j < nb_basis_fct; j++)
	  if (abs(vec_phi(j)) > 1e-12)
	    CollMatrix.Get(i, j) = vec_phi(j);
      }
  }


  //! initial guess
  void FindInitGuess(VectReal_wp& cn_init)
  {
    // initial guess, cn = 1
    int nGrid = c_discrete->GetNbDof();
    cn_init.Reallocate(nGrid);
    cn_init.Fill(Real_wp(1));
  }


  //! computes the solution for a given set of parameters cn
  void Init(const VectReal_wp& cn)
  {
    // c is evaluated on quadrature points
    c_discrete->Init(cn);
    VectReal_wp eval_rho(points_quad.GetM());
    Vector<Real_wp> phi;
    for (int i = 0; i < points_quad.GetM(); i++)
      eval_rho(i) = c_discrete->Evaluate(points_quad(i));
    
    var_laplace.ref_rho(1).SetUserIndex(eval_rho, Real_wp(1));
    
    // Computing geometry quantities
    var_laplace.ComputePhysicalCoefficients();
    var_laplace.ComputeMassMatrix();
    
    // Computing mass and stiffness matrices
    GlobalGenericMatrix<Real_wp> nat_mat;
    int N = var_laplace.GetNbDof();
    int r = var_laplace.mesh.GetOrder();
    mat_lu.Reallocate(N, N, r, r);
    mat_lu.Zero();

    // Finite element matrix is computed and factorized
    var_laplace.AddMatrixFEM(mat_lu, nat_mat);
    var_laplace.AddBoundaryTerms(mat_lu, nat_mat);
    Seldon::GetLU( mat_lu, pivot);
    
    // Filling u with the right hand side
    sol_u.Reallocate(N);
    sol_u.Zero();
    for (int i = 1; i < N-1; i++)
      sol_u(i) = mat_mass(i) * (*f_src)(var_laplace.GetCoordinateDof(i));
    
    // hetereogeneous Dirichlet condition
    sol_u(0) = ua;
    sol_u(N-1) = ub;

    // u is now the solution
    Seldon::SolveLU(mat_lu, pivot, sol_u);
  }
  

  //! computes the observable vector
  void Evaluate(VectReal_wp& y)
  {
    // subtracting part due to inhomogeneous Dirichlet
    y.Reallocate(var_laplace.GetNbDof()-2);
    for (int i = 1; i < var_laplace.GetNbDof()-1; i++)
      y(i-1) = sol_u(i) - (ua + (ub-ua)*var_laplace.GetCoordinateDof(i));
  }
  
  
  //! Applies jacobian matrix of F to vector h, res = DF(cn) h
  void Derivative(const VectReal_wp& h, VectReal_wp& res)
  {
    int N = var_laplace.GetNbDof();
    
    VectReal_wp rhs(N-2);
    
    Mlt(CollMatrix, h, rhs);
    VectReal_wp rhs_(N);
    rhs_(0) = 0.0;
    rhs_(N-1) = 0.0;
    for (int i = 1; i < N-1; i++)
      rhs_(i) = - mat_mass(i)*rhs(i-1)*sol_u(i);
    
    Seldon::SolveLU(mat_lu, pivot, rhs_);
    res.Reallocate(N-2);
    for( int i = 1; i < N-1; i++)
      res(i-1) = rhs_(i);

  }

  
  //! Applies adjoint of jacobian matrix of F to vector G, res = DF*(cn) h
  void Adjoint(const VectReal_wp& G, VectReal_wp& res)
  {
    int N = var_laplace.GetNbDof();
    
    VectReal_wp aux(N-2);
    VectReal_wp rhs_(N);
    rhs_(0) = 0.0;
    rhs_(N-1) = 0.0;
    for (int i = 1; i < N-1; i++)
      rhs_(i) = G(i-1);
    
    Seldon::SolveLU(mat_lu, pivot, rhs_);
    
    for( int i = 1; i < N-1; i++)
      aux(i-1) = -mat_mass(i)*rhs_(i)*sol_u(i);
    
    res.Reallocate(c_discrete->GetNbDof());
    Mlt(SeldonTrans, CollMatrix, aux, res);
  }


  //! Computes the exact measures y, parameter c will
  //! solve the non-linear least-squares problem F(c) = y
  void EvaluateDifferently(VectReal_wp& y)
  {
    // y will store the solution evaluated with exact c
    
    // evaluating the exact c on quadrature points
    VectReal_wp eval_rho(points_quad.GetM());
    for (int i = 0; i < points_quad.GetM(); i++)
      eval_rho(i) = (*c_exact)(points_quad(i));
    
    var_laplace.ref_rho(1).SetUserIndex(eval_rho, Real_wp(1));
    
    // computing geometry quantities
    var_laplace.ComputeMassMatrix();
    
    // Computing finite element matrix and factorizing it
    Montjoie::GlobalGenericMatrix<Real_wp> nat_mat;
    int N = var_laplace.GetNbDof();
    int r = var_laplace.mesh.GetOrder();
    mat_lu.Reallocate(N, N, r, r);
    
    var_laplace.AddMatrixFEM(mat_lu, nat_mat);
    var_laplace.AddBoundaryTerms(mat_lu, nat_mat);
    Seldon::GetLU(mat_lu, pivot);
    
    // computing right hand side
    VectReal_wp rhs(N);
    for (int i = 1; i < N-1; i++)
      rhs(i) = mat_mass(i) * (*f_src)(var_laplace.GetCoordinateDof(i));
    
    rhs(0) = ua;
    rhs(N-1) = ub;
    
    // solving the linear system
    Seldon::SolveLU(mat_lu, pivot, rhs);
    
    y.Reallocate(N-2);
    for (int i = 1; i < var_laplace.GetNbDof()-1; i++)
      y(i-1) = rhs(i);
    
    // subtracting part due to inhomogeneous Dirichlet
    for (int i = 1; i < var_laplace.GetNbDof()-1; i++)
      y(i-1) -= (ua + (ub-ua)*var_laplace.GetCoordinateDof(i));
  }
  
};


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var_laplace,
	    const string& input_file, const string& type_element)
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
  var.WriteDatas(x_sol);

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
      
      EllipticProblem<HelmholtzEquation1D> vars;
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
