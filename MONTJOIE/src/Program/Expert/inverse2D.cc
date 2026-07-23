#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
#include "Regularization/Itreg.hxx"

#include "Elliptic/MultiDimPhysicalIndex.cxx"

using namespace Montjoie;

// exact c
Real_wp pol_sin(const R2& x)
{
  //return 1.3;
  return 1.3 + 0.3*cos(x(0))*sin(x(1));
}

template<class Dimension>
class InverseProblem : public InputDataProblem_Base, public itreg::ForwardOperator<Real_wp>
{
protected:
  typedef typename Dimension::R_N R_N;
  
  //! pointer to the exact c
  Real_wp (*c_exact)(const R_N&);
  
  //! object for solving the Laplacian
  EllipticProblem<LaplaceEquation<Dimension> > var_laplace;
  
  //! object solving the linear system
  All_LinearSolver linear_solver;
  
  //! stored solution
  VectReal_wp sol_u;

  //! stored source
  VectReal_wp source_rhs;
  
  //! exact values of C stored
  VectReal_wp C_store_exact;
  
  //! object discretizing c
  VariableParameterND_Base<Real_wp, Dimension>* c_discrete;
  
  // quadrature points where c needs to be evaluated
  Vector<Vector<R_N> > points_quad;

  Matrix<Real_wp, General, ArrayRowSparse> CollMatrix;
  
  //! parameters associated with the regularization method
  Vector<string> regularization_parameter;
  Vector<string> approx_c_parameter;
  
public:
  //! Default constructor
  InverseProblem() : linear_solver(var_laplace)
  {
    c_exact = &pol_sin;
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
    c_discrete->Init(cn);
    
    // c is written on the output grids defined in var_laplace
    for (int n = 0; n < var_laplace.var_grid.GetM(); n++)
      {
	GridInterpolationFull<Dimension>& var_gr = var_laplace.var_grid(n);
	GridInterpolation<Dimension>& var_interp = var_laplace.all_points_display;
	const IVect& list_points = var_gr.GetPointNumber();
	const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();
	
	int nnz = list_points.GetM();
	int Nglob = var_gr.GetNbGlobalPoints();
	VectReal_wp interp_index(Nglob), exact_index(Nglob);
	interp_index.Zero(); exact_index.Zero();
	for (int i1 = 0; i1 < nnz; i1++)
	  {
	    int i = list_points(i1);
	    if (i < 0)
	      continue;
	    
	    R_N pt_glob = var_interp.GetGlobalCoordinate(i);
	    int num_elem = var_interp.GetElementNumber(i);
	    R_N pt_loc = var_interp.GetLocalCoordinate(i);
	    Real_wp cn_ex = (*c_exact)(pt_glob);
	    Real_wp cn_app = c_discrete->Evaluate(pt_glob, num_elem, pt_loc);
	    
	    int iglob = list_points_proc(0)(i);
	    interp_index(iglob) = cn_app;
	    exact_index(iglob) = cn_ex;
	  }

	WriteMatlab(exact_index, var_gr, string("Cexact_G") + to_str(n) + ".dat", 
		    var_laplace.OutputWrittenInDoublePrecision());
	
	WriteMatlab(interp_index, var_gr, string("Capp_G") + to_str(n) + ".dat", 
		    var_laplace.OutputWrittenInDoublePrecision());
      }
  }
  
  
  //! initialization of montjoie object var_laplace
  void Construct(const string& parameter_file_montjoie, const string& name_element)
  {
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(parameter_file_montjoie, lines_data_file, MPI_COMM_WORLD);
#else
    ReadLinesFile(parameter_file_montjoie, lines_data_file);
#endif

    var_laplace.InitIndices(PhysicalConstant::nb_max_indices);
    
    // the input file is read
    ReadInputFile(lines_data_file, linear_solver);
    ReadInputFile(parameter_file_montjoie, var_laplace);
    ReadInputFile(parameter_file_montjoie, *this);
    
    // construction of mesh and finite element
    bool split_mesh = true;
    if (var_laplace.GetNbProcPerMode() == 1)
      split_mesh = false;
    
    var_laplace.ComputeMeshAndFiniteElement(name_element, split_mesh);
    var_laplace.PerformOtherInitializations();
    
    // selection of solver
    Dimension dim; Symmetric property;
    linear_solver.SelectOptimalLinearSolver(var_laplace.mesh_num.GetOrder(),
					    var_laplace.GetNbDof(), dim, property);

    if (var_laplace.GetNbProcPerMode() == 1)
      var_laplace.mesh.Write("test.mesh");
    
    // computation of Ji, DFi on the mesh
    var_laplace.ComputeMassMatrix(true, false);
    points_quad = var_laplace.Glob_PointsQuadrature;
    var_laplace.Glob_PointsQuadrature.Clear();
    
    var_laplace.ComputeQuasiPeriodicPhase();
    
    // the source is computed
    var_laplace.ComputeRightHandSide(source_rhs);
    
    // the discretization for c is constructed
    c_discrete = VariableParameterND_Base<Real_wp, Dimension>::
      GetNewApproximation(var_laplace, approx_c_parameter);
    
    // the exact value of c (projection in the discrete space) is computed
    VectReal_wp evalC_exact;
    c_discrete->ProjectFunction(*c_exact, evalC_exact);
    C_store_exact = evalC_exact;
    
    int nb_pts_quad = 0;
    for (int i = 0; i < points_quad.GetM(); i++)
      nb_pts_quad += points_quad(i).GetM();
    
    // computation of CollMatrix that is needed
    // to compute the interpolation of c on the quadrature points
    int nb_basis_fct = c_discrete->GetNbDof();
    CollMatrix.Reallocate(nb_pts_quad, nb_basis_fct);    
    Vector<Real_wp, VectSparse> vec_phi;
    // here CollMatrix(i, j) = phi_j(points_quad(i))
    int irow = 0;
    for (int i = 0; i < points_quad.GetM(); i++)
      for (int j = 0; j < points_quad(i).GetM(); j++)
	{
	  const ElementReference<Dimension, 1>& Fb = var_laplace.GetReferenceElement(i);
	  c_discrete->ComputeValuesPhiRef(i, j, points_quad(i)(j),
					  Fb.PointsND(j), vec_phi);
	  
	  for (int j = 0; j < vec_phi.GetM(); j++)
	    CollMatrix.Get(irow, vec_phi.Index(j)) = vec_phi.Value(j);

	  irow++;
	}
    
    var_laplace.print_level = -1;
    linear_solver.SetPrintLevel(-1);
  }


  //! initial guess
  void FindInitGuess(VectReal_wp& cn_init)
  {
    // initial guess, cn = 1
    int nGrid = c_discrete->GetNbDof();
    cn_init.Reallocate(nGrid);
    cn_init.Fill(Real_wp(1));
  }


  //! c is given to object var_laplace
  void SetVariableIndex(const VectReal_wp& eval_c)
  {
    PhysicalVaryingMedia<Dimension, Real_wp>& rho = var_laplace.ref_rho(1).GetVaryingIndex();
    var_laplace.ref_rho(1).SetMediaType(rho.USER);
    rho.Reallocate(var_laplace.mesh.GetNbElt(), false);
    int num = 0;
    for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
      {
	rho.ReallocateOnElement(i, points_quad(i).GetM(), false);
	for (int j = 0; j < points_quad(i).GetM(); j++)
	  rho.SetCoefficient(i, j, eval_c(num++));
      }
    
    // Computing geometry quantities
    var_laplace.ComputeMassMatrix(false);
  }

  
  //! computes the solution for a given set of parameters cn
  void Init(const VectReal_wp& cn)
  {
    // c is evaluated on quadrature points
    c_discrete->Init(cn);
    VectReal_wp eval_c(CollMatrix.GetM());
    Mlt(CollMatrix, cn, eval_c);
    
    SetVariableIndex(eval_c);
    
    // Finite element matrix is computed and factorized
    GlobalGenericMatrix<Real_wp> nat_mat;
    linear_solver.PerformFactorizationStep(nat_mat);
    
    // the solution is computed
    sol_u.Reallocate(var_laplace.GetNbDof());
    sol_u = source_rhs;
    linear_solver.ComputeSolution(sol_u, nat_mat);
  }
  

  //! computes the observable vector
  void Evaluate(VectReal_wp& y)
  {
    y = sol_u;
  }

  
  void MltMass(const Real_wp& alpha, const VectReal_wp& c_quad, VectReal_wp& y)
  {
    TinyVector<VectReal_wp, 1> Uh_dof, Uh_quad;
    y.Zero();
    int offset = 0;
    for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension, 1>& Fb = var_laplace.GetReferenceElement(i);
	var_laplace.GetLocalUnknownVector(sol_u, i, Uh_dof);

	Uh_quad(0).Reallocate(Fb.GetNbPointsQuadratureInside());
	Fb.ApplyChTranspose(Uh_dof(0), Uh_quad(0));
	
	bool affine = var_laplace.mesh.IsElementAffine(i);
	if (affine)
	  for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	    Uh_quad(0)(k) *= var_laplace.Glob_jacobian(i)(0)*Fb.WeightsND(k)*c_quad(offset+k);
	else
	  for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	    Uh_quad(0)(k) *= var_laplace.Glob_jacobian(i)(k)*c_quad(offset+k);
	
	Fb.ApplyCh(Uh_quad(0), Uh_dof(0));
	var_laplace.AddLocalUnknownVector(alpha, Uh_dof, i, y);

	offset += points_quad(i).GetM();
      }
  }


  void MltMassTranspose(const Real_wp& alpha, const VectReal_wp& y, VectReal_wp& ProdC)
  {
    TinyVector<VectReal_wp, 1> Uh_dof, Uh_quad, Vh_dof, Vh_quad;
    ProdC.Zero();
    int offset = 0;
    for (int i = 0; i < var_laplace.mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension, 1>& Fb = var_laplace.GetReferenceElement(i);
	var_laplace.GetLocalUnknownVector(sol_u, i, Uh_dof);

	var_laplace.GetLocalUnknownVector(y, i, Vh_dof);

	Uh_quad(0).Reallocate(Fb.GetNbPointsQuadratureInside());
	Fb.ApplyChTranspose(Uh_dof(0), Uh_quad(0));

	Vh_quad(0).Reallocate(Fb.GetNbPointsQuadratureInside());
	Fb.ApplyChTranspose(Vh_dof(0), Vh_quad(0));
	
	bool affine = var_laplace.mesh.IsElementAffine(i);
	if (affine)
	  for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	    ProdC(offset+k) = alpha*var_laplace.Glob_jacobian(i)(0)*Fb.WeightsND(k)*Uh_quad(0)(k)*Vh_quad(0)(k);
	else
	  for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
	    ProdC(offset+k) = alpha*var_laplace.Glob_jacobian(i)(k)*Uh_quad(0)(k)*Vh_quad(0)(k);
	
	offset += points_quad(i).GetM();
      }
  }
  
  
  //! Applies jacobian matrix of F to vector h, res = DF(cn) h
  void Derivative(const VectReal_wp& h, VectReal_wp& res)
  {
    int N = var_laplace.GetNbDof();

    // evaluation of h on quadrature points
    VectReal_wp rhs(CollMatrix.GetM());
    Mlt(CollMatrix, h, rhs);

    // multiplication by mass matrix
    VectReal_wp rhs_(N);
    MltMass(-Real_wp(1), rhs, rhs_);
    
    // solving linear system
    GlobalGenericMatrix<Real_wp> nat_mat;
    res.Reallocate(N); res = rhs_;
    linear_solver.ComputeSolution(res, nat_mat);
  }

  
  //! Applies adjoint of jacobian matrix of F to vector G, res = DF*(cn) h
  void Adjoint(const VectReal_wp& G, VectReal_wp& res)
  {
    int N = var_laplace.GetNbDof();
    
    // solving linear system
    VectReal_wp rhs_(N), G_(G);
    GlobalGenericMatrix<Real_wp> nat_mat;
    rhs_ = G_;
    linear_solver.ComputeSolution(rhs_, nat_mat);
    
    // multiplication by mass matrix
    VectReal_wp aux(CollMatrix.GetM());
    MltMassTranspose(-Real_wp(1), rhs_, aux);
        
    // transpose of CollMatrix
    res.Reallocate(CollMatrix.GetN());    
    Mlt(SeldonTrans, CollMatrix, aux, res);
  }


  //! Computes the exact measures y, parameter c will
  //! solve the non-linear least-squares problem F(c) = y
  void EvaluateDifferently(VectReal_wp& y)
  {
    // y will store the solution evaluated with exact c
    
    // evaluating the exact c on quadrature points
    VectReal_wp eval_c(CollMatrix.GetM());
    int num = 0;
    for (int i = 0; i < points_quad.GetM(); i++)
      for (int j = 0; j < points_quad(i).GetM(); j++)
	eval_c(num++) = (*c_exact)(points_quad(i)(j));
    
    SetVariableIndex(eval_c);    
    
    // Computing finite element matrix and factorizing it
    Montjoie::GlobalGenericMatrix<Real_wp> nat_mat;
    linear_solver.PerformFactorizationStep(nat_mat);
    
    // the solution is computed
    y.Reallocate(var_laplace.GetNbDof());
    y = source_rhs;
    linear_solver.ComputeSolution(y, nat_mat);    
  }
  
};


template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var_laplace,
	    const string& input_file, const string& type_element)
{
  // forward solver is initialized
  InverseProblem<typename TypeEquation::Dimension> var;
  var.Construct(input_file, type_element);

  // testing derivatives
  /*VectReal_wp c0;
  var.FindInitGuess(c0);
  c0.FillRand(); c0 *= 1e-9;
  
  Real_wp h = 1e-16;
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
      
      EllipticProblem<LaplaceEquation<Dimension2> > vars;
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
