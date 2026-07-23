#ifndef MONTJOIE_FILE_MODAL_SOURCE_LAPLACE_CXX

namespace Montjoie
{

  template<class Dimension>
  void ModalSourceLaplace_Dim<Dimension>
  ::GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectReal_wp& coef_eigenval,
				    const Real_wp& impedance, int boundary_condition)
  {
    // no modification for real numbers
    coef_eigenval.Reallocate(eigenval.GetM());
    coef_eigenval.Fill(1.0);
  }


  template<class Dimension>
  void ModalSourceLaplace_Dim<Dimension>
  ::GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectComplex_wp& coef_eigenval,
				    const Complex_wp& impedance, int boundary_condition)
  {
    coef_eigenval.Reallocate(eigenval.GetM());
    coef_eigenval.Fill(1.0);    
    for (int k = 0; k < eigenval.GetM(); k++)
      {
	// computing beta = sqrt(omega^2 - lambda)
	Complex_wp beta = sqrt(Complex_wp(var_problem.GetSquareOmega()-eigenval(k), 0));
	if (boundary_condition == BoundaryConditionEnum::LINE_NEUMANN)
	  coef_eigenval(k) = -Iwp*beta*impedance;
	else if (boundary_condition == BoundaryConditionEnum::LINE_ABSORBING)
	  coef_eigenval(k) = -Iwp*(beta+var_problem.GetOmega())*impedance;
      }
  }

  
#ifdef MONTJOIE_WITH_TWO_DIM  
  template<>
  void ModalSourceLaplace_Dim<Dimension2>
  ::EvaluateModeSection(const Mesh<Dimension2>& result_mesh, int proc_mode,
			const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
			const Vector<int>& RotElement, const Vector<VectReal_wp>& LocalQuadPoints,
			Vector<Vector<VectReal_wp> >& eval_mode, VectReal_wp& eigenval)
  {
    R2 ptA = result_mesh.Vertex(0);
    R2 ptB = result_mesh.Vertex(1);
    Real_wp norm_AB;
    norm_AB = Norm2(ptA-ptB);
    Real_wp kwave = pi_wp/norm_AB;
    
    const IVect& number_mode_combine = modal_source.GetModeNumberToCombine();
    eval_mode.Reallocate(number_mode_combine.GetM());
    for (int k = 0; k < number_mode_combine.GetM(); k++)
      eval_mode(k).Reallocate(IndexBoundary.GetM());

    bool mode_cos = true;
    switch (modal_source.GetLateralBoundaryCondition())
      {
      case BoundaryConditionEnum::LINE_DIRICHLET:
	mode_cos = false;
	break;
      case BoundaryConditionEnum::LINE_NEUMANN:
	mode_cos = true;
	break;
      default:
	{
	  cout << "General case not handled" << endl;
	  abort();
	}
      }

    eigenval.Reallocate(number_mode_combine.GetM());
    for (int k = 0; k < number_mode_combine.GetM(); k++)
      {
	int number_mode = number_mode_combine(k);
	if (mode_cos)
	  eigenval(k) = square(number_mode*kwave);
	else
	  eigenval(k) = square((number_mode+1)*kwave);
      }
    
    for (int i = 0; i < IndexBoundary.GetM(); i++)
      {
	int num_elem = IndexElement(i);
	int num_loc = IndexBoundary(i);
	R2 ptC = var_problem.mesh.Vertex(var_problem.mesh.Element(num_elem).numVertex(num_loc));
	int nv = (num_loc+1)%var_problem.mesh.Element(num_elem).GetNbVertices();
	R2 ptD = var_problem.mesh.Vertex(var_problem.mesh.Element(num_elem).numVertex(nv));
	R2 point; Real_wp norm_CD;
	int nb_points = LocalQuadPoints(i).GetM();
	for (int k = 0; k < number_mode_combine.GetM(); k++)
	  eval_mode(k)(i).Reallocate(nb_points);
	
	for (int j = 0; j < nb_points; j++)
	  {
	    Real_wp lambda = LocalQuadPoints(i)(j);
	    point.Zero();
	    Add(1.0-lambda, ptC, point); Add(lambda, ptD, point);
	    norm_CD = point.Distance(ptA);
	    
	    for (int k = 0; k < number_mode_combine.GetM(); k++)
	      {
		int number_mode = number_mode_combine(k);
		if (mode_cos)
		  eval_mode(k)(i)(j) = cos(norm_CD*kwave*number_mode);
		else
		  eval_mode(k)(i)(j) = sin(norm_CD*kwave*(number_mode+1));
	      }
	  }	    	    
      }
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM  
  template<>
  void ModalSourceLaplace_Dim<Dimension3>
  ::EvaluateModeSection(const Mesh<Dimension3>& result_mesh, int proc_mode,
			const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
			const IVect& RotElement, const Vector<VectR2>& LocalQuadPoints,
			Vector<Vector<VectReal_wp> >& eval_mode, VectReal_wp& eigenval)
  {
    // object used to compute eigenvalues
    EllipticProblem<LaplaceEquation<Dimension2> > var_eig;    
    All_LinearSolver glob_solver(var_eig);
    IVect RotationFaceElement = RotElement;
    
    // resolution parameters
    glob_solver.SetDirectSolver();
    
    // physical properties and boundary conditions
    var_eig.InitIndices(var_problem.GetNbPhysicalIndices()-1);
    var_eig.mesh.ResizeNbReferences(var_problem.mesh.GetNbReferences());
    var_eig.mesh.SetBoundaryCondition(var_problem.mesh);
    this->CopyIndices(var_eig);
    
    // mesh and finite element are constructed
    R3 vec_u, vec_v;
    modal_source.ConstructProblem(result_mesh, var_problem, var_eig, RotationFaceElement,
				  string("TRIANGLE_LOBATTO"), vec_u, vec_v);

    // eigenvalues and eigenvectors are computed
    EigenProblemMontjoie<Real_wp> eigen_solver(var_eig, glob_solver);
    VectReal_wp eigen_values;
    Matrix<Real_wp, General, ColMajor> eigen_vectors;    

    modal_source.ComputeEigenvalues(eigen_solver, var_eig, eigen_values, eigen_vectors);     
    
    // projectors are computed
    Vector<FiniteElementInterpolator> proj;
    TinyVector<Vector<Matrix<int> >, 2> RotationPoints;
    
    modal_source.ComputeProjector(var_problem, var_eig, IndexBoundary, IndexElement,
				  LocalQuadPoints, proj, RotationPoints);

    // we store the required modes
    const MeshNumbering<Dimension2>& mesh_num_eig = var_eig.GetMeshNumbering(0);
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    VectReal_wp mode_dof(mesh_num_eig.GetNbDof());
    VectReal_wp Xmode(eigen_vectors.GetM());
    const IVect& number_mode_combined = modal_source.GetModeNumberToCombine();
    eigenval.Reallocate(number_mode_combined.GetM());
    eval_mode.Reallocate(number_mode_combined.GetM());
    for (int k = 0; k < number_mode_combined.GetM(); k++)
      {
	int number_mode = number_mode_combined(k);
	eigenval(k) = eigen_values(number_mode);
	for (int i = 0; i < eigen_vectors.GetM(); i++)
	  Xmode(i) = eigen_vectors(i, number_mode); 

        eigen_solver.ExpandVector(Xmode, mode_dof, true);
        
	// uncomment the following lines to visualize the mode
	/*
	Vector<string> param(2);
	param(0) = string("DiffracModeLaplace");
	param(1) = string("ModeLaplace");
	var_eig.SetInputData(string("FileOutputPlane"), param);
	
	param.Reallocate(3);
	param(0) = "AUTO";
	param(1) = "201";
	param(2) = "201";
	var_eig.SetInputData(string("SismoPlane"), param);
	
	var_eig.InitVarGrid();
	var_eig.ComputeVarGrid();
	
	var_eig.WriteDatas(mode_dof);*/
	
	// interpolating the mode on quadrature points of each required face
	int nb_faces = IndexBoundary.GetM();
	eval_mode(k).Reallocate(nb_faces);
	VectReal_wp u_loc, u_quad;
	for (int i = 0; i < nb_faces; i++)
	  {
	    // we retrieve value of u on the element
	    int nb_dof_loc = var_eig.GetNbLocalDof(i);
	    u_loc.Reallocate(nb_dof_loc);
	    for (int j = 0; j < nb_dof_loc; j++)
	      u_loc(j) = mode_dof(mesh_num_eig.Element(i).GetNumberDof(j));
	    
	    // projection on quadrature points
	    int num_elem = IndexElement(i);
	    int num_loc = IndexBoundary(i);
	    int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	    int r = mesh_num.GetOrderQuadrature(num_face);
	    int type_elt = var_problem.mesh.Boundary(num_face).GetHybridType();	    
	    int nb_pts_quad = LocalQuadPoints(i).GetM();
	    u_quad.Reallocate(nb_pts_quad);
	    proj(r).ProjectScalar(u_loc, u_quad, type_elt);
	    
	    // eval_mode is filled
	    eval_mode(k)(i).Reallocate(nb_pts_quad);
	    for (int j = 0; j < nb_pts_quad; j++)
	      {
		const Matrix<int>& RotationNum = RotationPoints(type_elt)(r);		
		int jrot = RotationNum(RotationFaceElement(i), j);
		eval_mode(k)(i)(j) = u_quad(jrot);
	      }
	  }

      }

  }
#endif

  //! computation of the modes of the considered section
  template<class Complexe, class Dimension>
  void ModalSourceLaplace<Complexe, Dimension>
  ::EvaluateMode(const Mesh<Dimension>& mesh, int proc_mode, const Vector<int>& IndexBoundary, 
		 const Vector<int>& IndexElement, const Vector<int>& RotElement,
		 const Vector<typename DimensionB::VectR_N>& LocalQuadPoints,
		 Vector<Vector<Vector<Complexe> > >& eval)
  {
    // evaluates the required modes (considered as real-valued)
    Vector<Vector<VectReal_wp> > eval_mode; VectReal_wp eigenval;
    ModalSourceLaplace_Dim<Dimension>::EvaluateModeSection(mesh, proc_mode, IndexBoundary, IndexElement,
							   RotElement, LocalQuadPoints, eval_mode, eigenval);
    
    // we get coefficients due to boundary condition
    Vector<Complexe> coef_eigenval;
    int ref = this->GetVolumeReference();
    this->GetCoefficientBoundaryCondition(eigenval, coef_eigenval,
					  var_laplace.coefficient_impedance_absorbing(ref),
					  this->boundary_condition_mode);
    
    // multiplication by coefficients provided in the data file
    for (int k = 0; k < eval_mode.GetM(); k++)
      coef_eigenval(k) *= this->coef_mode_combined(k);

    // filling eval by combining modes stored in eval_mode
    int m = 0;
    for (int i = 0; i < LocalQuadPoints.GetM(); i++)
      {
	int num_elem = IndexElement(i);
	int num_loc = IndexBoundary(i);
	int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	int nb_points = LocalQuadPoints(i).GetM();
	eval(m)(num_face).Reallocate(nb_points);
	for (int j = 0; j < nb_points; j++)
	  {
	    SetComplexZero(eval(m)(num_face)(j));
	    for (int k = 0; k < eval_mode.GetM(); k++)
	      eval(m)(num_face)(j) += coef_eigenval(k)*eval_mode(k)(i)(j);
	  }
      }
  }


  template<>
  void ModalSourceLaplace<Complex_wp, Dimension2>
  ::CopyIndices(VarHelmholtz_Base<Real_wp, Dimension1>& var_eig)
  {
  }

  template<>
  void ModalSourceLaplace<Real_wp, Dimension2>
  ::CopyIndices(VarHelmholtz_Base<Real_wp, Dimension1>& var_eig)
  {
  }


  template<class Complexe, class Dimension>
  void ModalSourceLaplace<Complexe, Dimension>
  ::CopyIndices(VarHelmholtz_Base<Real_wp, DimensionB>& var_eig)
  {
    Real_wp zero(0), rho;
    TinyMatrix<Real_wp, Symmetric, Dimension::dim_N-1, Dimension::dim_N-1> mu;
    for (int ref = 0; ref < var_problem.GetNbPhysicalIndices(); ref++)
      {
        rho = realpart(var_laplace.ref_rho(ref).GetConstant());
	mu.SetIdentity();
	mu *= realpart(var_laplace.ref_mu(ref).GetConstant()(0, 0));
	var_eig.ref_rho(ref).SetConstant(rho);
	var_eig.ref_sigma(ref).SetConstant(zero);
	var_eig.ref_mu(ref).SetConstant(mu);
      }
  }
  
}

#define MONTJOIE_FILE_MODAL_SOURCE_LAPLACE_CXX
#endif
