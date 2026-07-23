#ifndef MONTJOIE_FILE_POINTS_REFERENCE_CXX

namespace Montjoie
{  

#ifdef MONTJOIE_WITH_MPFR
  template<class Dimension>
  Real_wp FjInverseProblem<Dimension>::threshold_newton(Real_wp(4, MONTJOIE_DEFAULT_PRECISION)*epsilon_machine);

  template<class Dimension>
  Real_wp FjInverseProblem<Dimension>::coef_safety_solver(10, MONTJOIE_DEFAULT_PRECISION);
#else
  template<class Dimension>
  Real_wp FjInverseProblem<Dimension>::threshold_newton(4.0*epsilon_machine);

  template<class Dimension>
  Real_wp FjInverseProblem<Dimension>::coef_safety_solver(10.0);
#endif
  
  template<class Dimension>
  int FjInverseProblem<Dimension>::nb_max_iterations(50);
  
  template<class Dimension>
  int FjInverseProblem<Dimension>::non_linear_solver(0);

  template<class Dimension>
  int64_t FjInverseProblem<Dimension>::nb_iter_all(0);
  
  /************************
   * SetPoints<Dimension> *
   ************************/
  
  
  //! points are rotated of an angle teta
  template<class Dimension>
  void SetPoints<Dimension>::RotatePoints(const Real_wp& teta)
  {
    Real_wp cos_teta = cos(teta);
    Real_wp sin_teta = sin(teta);
    R_N pt1, pt2;
    
    for (int k = 0; k < PointsQuadrature.GetM(); k++)
      {
        pt1 = PointsQuadrature(k); pt2 = pt1;
        pt2(0) = cos_teta*pt1(0) - sin_teta*pt1(1);
        pt2(1) = sin_teta*pt1(0) + cos_teta*pt1(1);
        PointsQuadrature(k) = pt2;
      }
    
    for (int k = 0; k < PointsNodal.GetM(); k++)
      {
        pt1 = PointsNodal(k); pt2 = pt1;
        pt2(0) = cos_teta*pt1(0) - sin_teta*pt1(1);
        pt2(1) = sin_teta*pt1(0) + cos_teta*pt1(1);
        PointsNodal(k) = pt2;
      }

    for (int k = 0; k < PointsDof.GetM(); k++)
      {
        pt1 = PointsDof(k); pt2 = pt1;
        pt2(0) = cos_teta*pt1(0) - sin_teta*pt1(1);
        pt2(1) = sin_teta*pt1(0) + cos_teta*pt1(1);
        PointsDof(k) = pt2;
      }
              
    for (int k = 0; k < PointsBoundary.GetM(); k++)
      {
        pt1 = PointsBoundary(k); pt2 = pt1;
        pt2(0) = cos_teta*pt1(0) - sin_teta*pt1(1);
        pt2(1) = sin_teta*pt1(0) + cos_teta*pt1(1);
        PointsBoundary(k) = pt2;
      }

    for (int k = 0; k < PointsDofBoundary.GetM(); k++)
      {
        pt1 = PointsDofBoundary(k); pt2 = pt1;
        pt2(0) = cos_teta*pt1(0) - sin_teta*pt1(1);
        pt2(1) = sin_teta*pt1(0) + cos_teta*pt1(1);
        PointsDofBoundary(k) = pt2;
      }
  }
  
  
  //! points are translated from a constant vector
  template<class Dimension>
  void SetPoints<Dimension>::TranslatePoints(const R_N& vec_u)
  {
    for (int k = 0; k < PointsQuadrature.GetM(); k++)
      PointsQuadrature(k) += vec_u;
    
    for (int k = 0; k < PointsNodal.GetM(); k++)
      PointsNodal(k) += vec_u;

    for (int k = 0; k < PointsDof.GetM(); k++)
      PointsDof(k) += vec_u;
              
    for (int k = 0; k < PointsBoundary.GetM(); k++)
      PointsBoundary(k) += vec_u;

    for (int k = 0; k < PointsDofBoundary.GetM(); k++)
      PointsDofBoundary(k) += vec_u;
  }
  
  
  //! displays informations about class SetPoints
  template<class Dimension>
  ostream& operator <<(ostream& out, const SetPoints<Dimension>& e)
  {
    out<<"Nodal points "<<e.PointsNodal<<endl;
    out<<"Points on boundary "<<e.PointsBoundary<<endl;
    return out;
  }
  
  
  /**************************
   * SetMatrices<Dimension> * 
   **************************/
  
  
  //! allocate the array containing jacobian matrices on quadrature points of the boundary
  template<class Dimension>
  void SetMatrices<Dimension>::ReallocatePointsQuadratureBoundary(int N)
  {
    MatricesBoundary.Reallocate(N);
    NormaleQuadrature.Reallocate(N);
    DsQuadrature.Reallocate(N);
    K1_curve.Reallocate(N);
    K2_curve.Reallocate(N);
    K1_curve.Fill(0);
    K2_curve.Fill(0);
  }

  
  //! compute jacobian matrices DFi on nodal points of the boundary
  template<class Dimension>
  void SetMatrices<Dimension>::ComputeNodalBoundary(int num_loc, const ElementReference_Dim<Dimension>& Fb)
  {
    int nb_pts = Fb.GetNbNodalBoundary(num_loc);
    MatricesNodalBoundary.Reallocate(nb_pts);
    invMatricesNodalBoundary.Reallocate(nb_pts);    
    for (int i = 0; i < nb_pts; i++)
      {
        int node = Fb.GetNodalNumber(num_loc, i);
        MatricesNodalBoundary(i) = MatricesNodal(node);
        GetInverse(MatricesNodalBoundary(i), invMatricesNodalBoundary(i));
      }
  }
  
  
  //! normale are rotated of an angle teta
  template<class Dimension>
  void SetMatrices<Dimension>::RotateNormale(const Real_wp& teta)
  {
    Real_wp cos_teta = cos(teta);
    Real_wp sin_teta = sin(teta);
    R_N pt1, pt2;
    for (int k = 0; k < NormaleQuadrature.GetM(); k++)
      {
        pt1 = NormaleQuadrature(k); pt2 = pt1;
        pt2(0) = cos_teta*pt1(0) - sin_teta*pt1(1);
        pt2(1) = sin_teta*pt1(0) + cos_teta*pt1(1);
        NormaleQuadrature(k) = pt2;
      }
  }
  
  
  //! displays informations about class SetMatrices
  template<class Dimension>
  ostream& operator <<(ostream& out, const SetMatrices<Dimension>& e)
  {
    out<<"Jacobian matrices on quadrature points "<<e.MatricesQuadrature<<endl;
    out<<"Jacobian matrices on nodal points "<<e.MatricesNodal<<endl;
    out<<"Jacobian matrices on boundary "<<e.MatricesBoundary<<endl;
    out<<"Normales on boundary "<<e.NormaleQuadrature<<endl;
    out<<"Surfacic element Ds "<<e.DsQuadrature<<endl;
    return out;
  }
  
  
  /********************
   * FjInverseProblem *
   ********************/
  
  
  //! constructor with parameters
  template<class Dimension>  
  FjInverseProblem<Dimension>::
  FjInverseProblem(const Mesh<Dimension>& mesh_, int nquad)
    : Fb(mesh_.GetReferenceElement(nquad)), mesh(mesh_), num_elem(nquad)
  {    
    this->m_ = Dimension::dim_N;
    this->n_ = Dimension::dim_N;
    
    mesh.GetVerticesElement(nquad, Vertices);
    affine = mesh.IsElementAffine(nquad);
    curved = mesh.Element(nquad).IsCurved();
    if (affine)
      {	
	R_N pt_loc;
	Fb.FjLinear(Vertices, pt_loc, F0);
	Fb.DFjLinear(Vertices, pt_loc, invDF);
	GetInverse(invDF);
      }
    else
      {
	Fb.FjElemNodal(Vertices, PointsElem, mesh, nquad);
	if (curved)
	  Fb.DFjElemNodal(Vertices, PointsElem, MatricesElem, mesh, nquad);
      }
    
    Rmax = Real_wp(0);
    for (int j = 0; j < Vertices.GetM(); j++)
      {
        Rmax = max(abs(Vertices(j)(0)), Rmax);
        Rmax = max(abs(Vertices(j)(1)), Rmax);
        if (Dimension::dim_N == 3)
          Rmax = max(abs(Vertices(j)(2)), Rmax);
      }
    
    if (non_linear_solver == NEWTON_SOLVER)
      {
	R_N res;
	newton_solver.ForceReevaluationJacobian(true);
	// newton_solver.EnableBacktracking();
	newton_solver.SetMaxNumberOfIterations(nb_max_iterations);
	newton_solver.Init(*this, res);    
      }
    else if (non_linear_solver == LVM_SOLVER)
      {
	// norm_DF is needed for Levenberg-Marquardt to obtain a good stopping
	// criterion for g = DF^T residual
	R_N pt_loc;
	Fb.DFjLinear(Vertices, pt_loc, invDF);
	norm_DF = MaxAbs(invDF);
      }
  }
  
  
  //! modification of the class with a line of the input file
  template<class Dimension>
  void FjInverseProblem<Dimension>::SetInputData(const string& description_field,
						 const VectString& parameters)
  {
    if (!description_field.compare("NonLinearSolver"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of FjInverseProblem" << endl;
	    cout << "NonLinearSolver needs more parameters, for instance :" << endl;
	    cout << "NonLinearSolver = MINPACK AUTO nb_iterations" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
	if (!parameters(0).compare("MINPACK"))
	  non_linear_solver = MINPACK_SOLVER;
	else if (!parameters(0).compare("LVM"))
	  non_linear_solver = LVM_SOLVER;
	else
	  non_linear_solver = NEWTON_SOLVER;
	
	if (parameters(1) == "AUTO")
          threshold_newton = 5.0*epsilon_machine;
        else
          to_num(parameters(1), threshold_newton);
        
	to_num(parameters(2), nb_max_iterations);
        
        if (parameters.GetM() > 3)
          to_num(parameters(3), coef_safety_solver);
      }
  }
  
  
  //! we solve non-linear system F_i(x) = pt
  template<class Dimension>
  bool FjInverseProblem<Dimension>::Solve(const R_N& point, R_N& res, int lvl)
  {
    if (affine)
      {
	R_N evalF(F0);
	evalF -= point;
	Mlt(invDF, evalF, res);
	res = -res;
	
	return !Fb.OutsideReferenceElement(Vertices, res, 3*coef_safety_solver*threshold_newton);
      }
    
    point_global = point;
    FindInitGuess(res);
    
    int nb_iter = 0;
    Real_wp epsilon = threshold_newton * Rmax;
    if (non_linear_solver == MINPACK_SOLVER)
      {
	VectReal_wp RControl(10); RControl.Fill(Real_wp(0));
	IVect Control(10);
	R_N fvec, scale;
	typename Dimension::MatrixN_N fjac;
	        
	// we set control parameters
	RControl(0) = epsilon;
	RControl(1) = Real_wp(1);
        	
	Control(0) = nb_max_iterations;
	Control(3) = 1;
	Control(4) = 0;       
        
	// we call function SolveMinpack defined in file NonLinear_Equations.cxx
	Montjoie::SolveMinpack(*this, res, fvec, fjac, scale, Control, RControl);
      }
    else if (non_linear_solver == NEWTON_SOLVER)
      {
	newton_solver.SetPrintLevel(lvl);
	newton_solver.SetStoppingCriterion(epsilon);
	this->ComputeAndFactoriseDiff(res, res);
 	newton_solver.Newton(*this, res, nb_iter, Real_wp(1));
      }
    else if (non_linear_solver == LVM_SOLVER)
      {
	SolveLeastSquaresLvm(*this, res, epsilon*norm_DF, threshold_newton,
			     nb_max_iterations, lvl);
      }
    
    // we compute the residual reached
    R_N X;
    Fb.Fj(Vertices, PointsElem, res, X, mesh, num_elem);
    X -= point_global;
    
    Real_wp residu = Norm2(X);
    if ((nb_iter <= 2) && (non_linear_solver == NEWTON_SOLVER)
        && (residu > 1e3*coef_safety_solver*epsilon))
      {
	// if Newton solver failed quickly, we select a better initial guess
	Real_wp dist_min = 1e300, dist;
	for (int i = 0; i < PointsElem.GetNbPointsNodal(); i++)
	  {
	    dist = AbsSquare(PointsElem.GetPointNodal(i) - point_global);
	    if (dist < dist_min)
	      {
		res = Fb.PointsNodalND(i);
		dist_min = dist;
	      }
	  }
	
	this->ComputeAndFactoriseDiff(res, res);
	newton_solver.Newton(*this, res, nb_iter, Real_wp(1));

	Fb.Fj(Vertices, PointsElem, res, X, mesh, num_elem);
	X -= point_global;
	residu = Norm2(X);
      }
    
    // if residual is small enough and point inside the element
    if ((residu <= coef_safety_solver*epsilon)
        &&(!Fb.OutsideReferenceElement(Vertices, res, 3*coef_safety_solver*threshold_newton)))
      return true;
    
    // point outside the element
    return false;
  }


  //! evaluation of function to solve, i.e. F_i(x) - pt
  template<class Dimension>
  void FjInverseProblem<Dimension>::EvaluateFunction(const R_N& x, R_N& fvec)
  {
    x_eval = x;
    if (curved)
      Fb.FjCurvePhi(PointsElem, x, fvec, phi);
    else
      Fb.Fj(Vertices, PointsElem, x, fvec, mesh, num_elem);
    fvec -= point_global;
  }
  
  
  //! jacobian matrix DF_i
  template<class Dimension>
  void FjInverseProblem<Dimension>::EvaluateJacobian(const R_N& x, MatrixN_N& fjac)
  {
    nb_iter_all++;
    if (curved)
      {
	if (x != x_eval)
	  {
	    x_eval = x;
	    Fb.ComputeValuesPhiNodalRef(x, phi);
	  }

	Fb.DFjCurvePhi(MatricesElem, x, fjac, phi);
      }
    else
      Fb.DFj(Vertices, PointsElem, x, fjac, mesh, num_elem);
  }
  
} // end namespace

#define MONTJOIE_FILE_POINTS_REFERENCE_CXX
#endif
