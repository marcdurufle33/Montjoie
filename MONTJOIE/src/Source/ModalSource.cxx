#ifndef MONTJOIE_FILE_MODAL_SOURCE_CXX

namespace Montjoie
{
  
  /****************************
   * ModalSourceBoundary_Base *
   ****************************/

  //! sets values of attriutes to default values
  void ModalSourceBoundary_Base::InitDefaultValues(int nb_ref_mesh)
  {
    is_ref_on_mode.Reallocate(nb_ref_mesh);
    is_ref_on_mode.Fill(0);
    
    nb_modes_to_compute = 10;

    shift_eigenval_mode = 0.1;
    lateral_condition = 0;
    boundary_condition_mode = 0;
    reference_volume = 1;

    nb_modes_to_combine = 1;
  }
  
  
#ifdef MONTJOIE_WITH_TWO_DIM  
  //! retrieves the mesh of the section and quadrature points
  void ModalSourceBoundary_Base
  ::GetMeshAndQuadraturePoints(const DistributedProblem<Dimension2>& var,
			       Mesh<Dimension2>& result_mesh, int proc_mode,
			       Vector<int>& IndexBoundary, Vector<int>& IndexElement,
			       Vector<int>& RotElement, Vector<VectReal_wp>& LocalQuadPoints)
  {
    SurfacicMesh<Dimension2> mesh;
    // we get surfacic mesh associated with references of the mode
    var.mesh.GetBoundaryMesh(1, mesh, this->is_ref_on_mode);
    
    R2 ptA, ptB; 
    Real_wp xmin(1e300), xmax(-1e300), ymin(1e300), ymax(-1e300);
    int ref_vol = -1;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int num_elem = var.mesh.BoundaryRef(mesh.ListeBoundaries(i)).numElement(0);
	ref_vol = var.mesh.Element(num_elem).GetReference();
	for (int k = 0; k < 2; k++)
	  {
	    int n0 = mesh.BoundaryRef(i).numVertex(k);
	    xmin = min(xmin, mesh.Vertex(n0)(0));
	    xmax = max(xmax, mesh.Vertex(n0)(0));
	    ymin = min(ymin, mesh.Vertex(n0)(1));
	    ymax = max(ymax, mesh.Vertex(n0)(1));
	  }
      }
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&xmin, &ptA(0), 1, MPI_DOUBLE, MPI_MIN, var.comm_group_mode);
    MPI_Allreduce(&ymin, &ptA(1), 1, MPI_DOUBLE, MPI_MIN, var.comm_group_mode);
    MPI_Allreduce(&xmax, &ptB(0), 1, MPI_DOUBLE, MPI_MAX, var.comm_group_mode);
    MPI_Allreduce(&ymax, &ptB(1), 1, MPI_DOUBLE, MPI_MAX, var.comm_group_mode);
    MPI_Allreduce(&ref_vol, &reference_volume, 1, MPI_INTEGER, MPI_MAX, var.comm_group_mode);
#else
    ptA.Init(xmin, ymin);
    ptB.Init(xmax, ymax);
    reference_volume = ref_vol;
#endif
    
    result_mesh.ReallocateVertices(2);
    result_mesh.Vertex(0) = ptA;
    result_mesh.Vertex(1) = ptB;
    
    IndexBoundary.Reallocate(mesh.GetNbBoundaryRef());
    IndexElement.Reallocate(mesh.GetNbBoundaryRef());
    LocalQuadPoints.Reallocate(mesh.GetNbBoundaryRef());
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int num_elem = var.mesh.BoundaryRef(mesh.ListeBoundaries(i)).numElement(0);
	int num_face = mesh.ListeBoundaries(i);
	int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_face);	
	int ref = var.mesh.BoundaryRef(mesh.ListeBoundaries(i)).GetReference();
	if ((!var.FirstOrderFormulationDG()) &&
	    (var.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_DIRICHLET) )
	  {
	    // For heteregeneous Dirichlet, the mode needs to be evaluated on dof points
	    // of the surface, which may differ from quadrature points
	    LocalQuadPoints(i) = var.PointsDofBoundary(num_elem, num_loc);
	  }
	else
	  LocalQuadPoints(i) = var.PointsQuadratureBoundary(num_elem, num_loc);
	
	IndexBoundary(i) = num_loc;
	IndexElement(i) = num_elem;
      }
  }
#endif


#ifdef MONTJOIE_WITH_THREE_DIM  
  //! retrieves the mesh of the section and quadrature points
  void ModalSourceBoundary_Base
  ::GetMeshAndQuadraturePoints(const DistributedProblem<Dimension3>& var,
			       Mesh<Dimension3>& result_mesh, int proc_mode,
			       Vector<int>& IndexBoundary, Vector<int>& IndexElement,
			       Vector<int>& RotElement, Vector<VectR2>& LocalQuadPoints)
  {
    SurfacicMesh<Dimension3> mesh;
    // we get surfacic mesh associated with reference ref_mode
    var.mesh.GetBoundaryMesh(1, mesh, is_ref_on_mode);
    
    if (var.GetNbProcPerMode() > 1)
      {
	cout << "Not yet implemented in parallel" << endl;
	abort();
      }
    
    result_mesh = static_cast<Mesh<Dimension3>& >(mesh);

    RotElement.Reallocate(mesh.GetNbBoundaryRef());
    IndexBoundary.Reallocate(mesh.GetNbBoundaryRef());
    IndexElement.Reallocate(mesh.GetNbBoundaryRef());
    LocalQuadPoints.Reallocate(mesh.GetNbBoundaryRef());
    reference_volume = -1;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {	
	int num_elem = var.mesh.BoundaryRef(mesh.ListeBoundaries(i)).numElement(0);
	int num_face = mesh.ListeBoundaries(i);
	int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_face);
	int rot = var.mesh.Element(num_elem).GetOrientationBoundary(num_loc);
	int ref = var.mesh.BoundaryRef(mesh.ListeBoundaries(i)).GetReference();
	reference_volume = var.mesh.Element(num_elem).GetReference();
	rot = MeshNumbering<Dimension3>
          ::GetOppositeOrientationFace(rot, var.mesh.Boundary(num_face));
	
	// changing the reference of the boundary
	result_mesh.BoundaryRef(i).SetReference(reference_volume);

	IndexBoundary(i) = num_loc;
	IndexElement(i) = num_elem;
	if ((!var.FirstOrderFormulationDG()) &&
	    (var.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_DIRICHLET) )
	  {
	    // For heteregeneous Dirichlet, the mode needs to be evaluated on dof points
	    // of the surface, which may differ from quadrature points
	    LocalQuadPoints(i) = var.PointsDofBoundary(num_elem, num_loc);
          }
	else
	  LocalQuadPoints(i) = var.PointsQuadratureBoundary(num_elem, num_loc);
	
	RotElement(i) = rot;
      }
  }

  void ModalSourceBoundary_Base
  ::ConstructProblem(const Mesh<Dimension3>& result_mesh,
		     const DistributedProblem<Dimension3>& var_problem,
		     DistributedProblem<Dimension2>& var_eig, IVect& RotationFaceElement,
		     const string& name_element, R3& vec_u, R3& vec_v) const
  {
    // various initializations
    int r = var_problem.GetMeshNumbering(0).GetOrder();
    var_eig.GetMeshNumbering(0).SetOrder(r);
    var_eig.mesh.SetGeometryOrder(r);
    var_eig.SetTypeEquation("none");
    
    // retrieving unit vectors of the plane (vec_u, vec_v -> orthogonal unit vectors,
    // ptC -> origin)
    R3 ptC, ptM;
    ptC = result_mesh.Vertex(0);
    ptM(0) = coef_plane(0); ptM(1) = coef_plane(1); ptM(2) = coef_plane(2);
    GetVectorPlane(ptM, vec_u, vec_v);
    if ((vec_u(0)*vec_v(1)-vec_u(1)*vec_v(0)) < 0)
      Mlt(-1.0, vec_v);
    
    // 3-D vertices are transformed into 2-D vertices
    // by using the center ptC and tangent vectors vec_u, vec_v of the plane
    var_eig.mesh.ReallocateVertices(result_mesh.GetNbVertices());
    for (int i = 0; i < result_mesh.GetNbVertices(); i++)
      {
	ptM = result_mesh.Vertex(i) - ptC;
	var_eig.mesh.Vertex(i)(0) = DotProd(ptM, vec_u);
	var_eig.mesh.Vertex(i)(1) = DotProd(ptM, vec_v);
      }
    
    // elements
    var_eig.mesh.ReallocateElements(result_mesh.GetNbBoundaryRef());
    IVect num; R2 vec_u2, vec_v2;
    VectBool SymmetrizedElement(result_mesh.GetNbBoundaryRef());
    SymmetrizedElement.Fill(false);
    for (int i = 0; i < result_mesh.GetNbBoundaryRef(); i++)
      {
	int nv = result_mesh.BoundaryRef(i).GetNbVertices();
	num.Reallocate(nv);
	for (int j = 0; j < nv; j++)
	  num(j) = result_mesh.BoundaryRef(i).numVertex(j);
	
	// we check if the new element is oriented counterclockwise
	// if not, SymmetrizedElement(i) will be true
	// (=> the face will be the symmetric of the face belonging to the volumic mesh)
	vec_u2 = var_eig.mesh.Vertex(num(1)) - var_eig.mesh.Vertex(num(0));
	vec_v2 = var_eig.mesh.Vertex(num(nv-1)) - var_eig.mesh.Vertex(num(0));
	if ( (vec_u2(0)*vec_v2(1)-vec_u2(1)*vec_v2(0)) < 0)
	  {
	    SymmetrizedElement(i) = true;
	    int itmp = num(nv-1); num(nv-1) = num(1); num(1) = itmp;
	  }
	
	// we take the references in result_mesh
	int ref = result_mesh.BoundaryRef(i).GetReference();
	var_eig.mesh.Element(i).Init(num, ref);
      }
    
    // referenced edges 
    var_eig.mesh.ReallocateBoundariesRef(result_mesh.GetNbEdgesRef());
    for (int i = 0; i < result_mesh.GetNbEdgesRef(); i++)
      {
	int n0 = result_mesh.EdgeRef(i).numVertex(0);
	int n1 = result_mesh.EdgeRef(i).numVertex(1);
	int ref = result_mesh.EdgeRef(i).GetReference();
	var_eig.mesh.BoundaryRef(i).Init(n0, n1, ref);
	
	if (result_mesh.GetCurveType(ref) != 0)
	  {
	    var_eig.mesh.SetCurveType(ref, var_eig.mesh.CURVE_FILE);
	    for (int k = 0; k < r-1; k++)
	      {
		ptM = result_mesh.GetPointInsideEdge(i, k) - ptC;
		vec_u2(0) = DotProd(ptM, vec_u);
		vec_u2(1) = DotProd(ptM, vec_v);
		var_eig.mesh.SetPointInsideEdge(i, k, vec_u2);
	      }
	  }
      }

    var_eig.mesh.SortBoundariesRef();    
    var_eig.mesh.ReorientElements();
    var_eig.mesh.FindConnectivity();
    var_eig.mesh.ProjectPointsOnCurves();
    
    var_eig.mesh.Write("surf.mesh");
    //VectReal_wp step_subdiv(5); step_subdiv.Fill(); Mlt(1.0/4, step_subdiv);
    //var_eig.mesh.SubdivideMesh(step_subdiv);
    //var_eig.mesh.Write("surf_subdiv.mesh");

    // numbering mesh and constructing finite element
    var_eig.ConstructFiniteElement(name_element);
    var_eig.GetMeshNumbering(0).NumberMesh();
    var_eig.ComputeNumberOfDofs();
    var_eig.GetBoundaryConditionProblem().TreatDirichletCondition();
    var_eig.FindElementsInsidePML();
    
    // changing rotation if necessary
    for (int i = 0; i < result_mesh.GetNbBoundaryRef(); i++)
      if (SymmetrizedElement(i))
	{
	  int rot = RotationFaceElement(i);
	  if (var_eig.mesh.Element(i).GetNbVertices() == 4)
	    RotationFaceElement(i) = (rot+4)%8;
	  else
	    RotationFaceElement(i) = (rot+3)%6;
	}    
  }


  void ModalSourceBoundary_Base
  ::ComputeEigenvalues(EigenProblemMontjoie<Real_wp>& eigen_solver,
		       DistributedProblem<Dimension2>& var_eig, VectReal_wp& eigen_values,
		       Matrix<Real_wp, General, ColMajor>& eigen_vectors,
		       int mode_eig, Real_wp coef_shift, bool init_arnoldi) const
  {
    int nb_modes = this->GetNbModesToCompute();
    Real_wp shift = this->GetEigenvalueShift()*coef_shift;
    if (mode_eig == -1)
      mode_eig = eigen_solver.SHIFTED_MODE;
        
    // number of modes to compute
    eigen_solver.SetNbAskedEigenvalues(nb_modes);
    if (init_arnoldi)
      eigen_solver.SetNbArnoldiVectors(5+2*nb_modes);
    
    // shifted mode is the fastest computational mode
    eigen_solver.SetComputationalMode(mode_eig);
    
    // the shift is a small value in order to obtain the smallest eigenvalues
    eigen_solver.SetTypeSpectrum(eigen_solver.CENTERED_EIGENVALUES, shift);
    
    FemMatrixFreeClass_Base<Real_wp> *Kh, *Mh;
    
    Mh = var_eig.GetComputationProblem().GetNewIterativeMatrix(Real_wp(0));
    Kh = var_eig.GetComputationProblem().GetNewIterativeMatrix(Real_wp(0));
    
    var_eig.ComputeMassMatrix(); 
    eigen_solver.InitMatrix(*Kh, *Mh);

    // computation of eigenvalues and eigenvectors
    VectReal_wp eigen_imag;
    GetEigenvaluesEigenvectors(eigen_solver, eigen_values, eigen_imag, eigen_vectors);
    
    Mlt(1.0/coef_shift, eigen_values);
    DISP(eigen_values); 
    
    eigen_solver.Clear();
    delete Kh; delete Mh;
  }


  void ModalSourceBoundary_Base
  ::ComputeProjector(const DistributedProblem<Dimension3>& var_problem,
		     const DistributedProblem<Dimension2>& var_eig,
		     const IVect& IndexBoundary, const IVect& IndexElement,
		     const Vector<VectR2>& LocalQuadPoints,
		     Vector<FiniteElementInterpolator>& proj,
		     TinyVector<Vector<Matrix<int> >, 2>& RotationPoints) const
  {
    int order_max = 0;
    for (int i = 0; i < IndexBoundary.GetM(); i++)
      {
	int num_elem = IndexElement(i);
	int num_loc = IndexBoundary(i);
	int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	order_max = max(order_max, var_problem.GetMeshNumbering(0).GetOrderQuadrature(num_face));
      }
    
    proj.Reallocate(order_max+1);
    RotationPoints(0).Reallocate(order_max+1);
    RotationPoints(1).Reallocate(order_max+1);
    TinyVector<VectR2, 4> points_div;
    TinyVector<bool, 2> order_used;
    for (int rm = 1; rm <= order_max; rm++)
      {
	order_used(0) = false; order_used(1) = false;
	points_div(0).Clear(); points_div(1).Clear();
	for (int i = 0; i < IndexBoundary.GetM(); i++)
	  {
	    int num_elem = IndexElement(i);
	    int num_loc = IndexBoundary(i);
	    int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	    int r = var_problem.GetMeshNumbering(0).GetOrderQuadrature(num_face);
	    int type_elt = var_problem.mesh.Boundary(num_face).GetHybridType();	    
	    if ((r == rm) && (!order_used(type_elt)))
	      {
		order_used(type_elt) = true;
		if (type_elt == 0)
		  MeshNumbering<Dimension3>::
		    GetRotationTriangularFace(LocalQuadPoints(i), RotationPoints(0)(r));
		else
		  MeshNumbering<Dimension3>::
		    GetRotationQuadrilateralFace(LocalQuadPoints(i), RotationPoints(1)(r));
		
		points_div(type_elt) = LocalQuadPoints(i);
	      }
	  }
	
	if (order_used(0) || order_used(1))
	  proj(rm).InitProjection(var_eig.mesh.GetReferenceElement(), points_div);
      }
  }
#endif
    
  
  /****************************
   * ModalSourceBoundary_Dim *
   ****************************/
  

  //! sets attributes with parameters of the data file
  template<class T, class Dimension>
  void ModalSourceBoundary_Dim<T, Dimension>::
  ReadSourceParameters(const IVect& boundary_condition, const IVect& ref,
		       const VectString& source_param)
  {    
    // is_ref_on_mode(ref) will be 1 if ref belongs to the mode
    this->is_ref_on_mode.Fill(0);
    for (int k = 0; k < ref.GetM(); k++)
      this->is_ref_on_mode(ref(k)) = 1;
    
    int nb = 0;
    if (source_param.GetM() < 6)
      {
	cout << "Provide more parameters for the definition of the mode" << endl;
	cout << "Current parameters are " << endl << source_param << endl;
	abort();
      }
    
    keyword_mode = source_param(nb++);
    
    this->shift_eigenval_mode = to_num<Real_wp>(source_param(nb++));
    this->nb_modes_to_compute = to_num<int>(source_param(nb++));
    this->nb_modes_to_combine = to_num<int>(source_param(nb++));
    
    this->number_mode_combined.Reallocate(nb_modes_to_combine);
    coef_mode_combined.Reallocate(nb_modes_to_combine);
    for (int k = 0; k < this->nb_modes_to_combine; k++)
      {
	if (source_param.GetM() < nb+2)
	  {
	    cout << "Provide more parameters for the definition of the mode" << endl;
	    cout << "Current parameters are " << endl << source_param << endl;
	    abort();
	  }

	this->number_mode_combined(k) = to_num<int>(source_param(nb++));
	coef_mode_combined(k) = to_num<T>(source_param(nb++));
      }    

    if (source_param.GetM() < nb+1)
      {
	cout << "Provide more parameters for the definition of the mode" << endl;
	cout << "Current parameters are " << endl << source_param << endl;
	abort();
      }
    
    // boundary condition associated with lateral faces
    lateral_condition = BoundaryConditionEnum::GetInteger(source_param(nb++));
  }

  
  //! kind of constructor for modal source
  /*!
    \param[in] num_mode number of the mode among all the modes that need to be evaluated
    \param[in] ref references associated with the section associated with the mode
    \param[in] param parameters 
   */
  template<class T, class Dimension>
  void ModalSourceBoundary_Dim<T, Dimension>
  ::ComputeMode(int num_mode, const IVect& ref,
		const VectString& param, Vector<Vector<Vector<T> > >& eval)
  {
    if (ref.GetM() <= 0)
      return;
    
    // parameters of the mode are treated
    int ref_mode = ref(0);
    this->ReadSourceParameters(var_problem.mesh.GetBoundaryCondition(), ref, param);    
    boundary_condition_mode = var_problem.mesh.GetBoundaryCondition(ref_mode);
    
    // finding the coefficients of the hyperplane where the section is
    int dim_N = Dimension::dim_N;
    this->coef_plane.Reallocate(dim_N+1); this->coef_plane.Fill(Real_wp(0));
    var_problem.mesh.FindParametersPlane(ref_mode, this->coef_plane);
    
    // computes a mesh of the section
    Mesh<Dimension> section_mesh;
    int proc_mode = 0;
    IVect IndexBoundary, IndexElement, RotElement;
    Vector<typename DimensionB::VectR_N> LocalQuadPoints;
    this->GetMeshAndQuadraturePoints(var_problem, section_mesh, proc_mode,
				     IndexBoundary, IndexElement, RotElement, LocalQuadPoints);
    
    // evaluates the required mode on the section, eval is modified
    this->EvaluateMode(section_mesh, proc_mode, IndexBoundary, IndexElement,
		       RotElement, LocalQuadPoints, eval);
  }

  
  //! computes mode and fills the array eval
  /*!
    \param[in] mesh mesh of the section of the mode
    \param[in] proc_mode processor that will store the global mesh
    \param[in] IndexBoundary IndexBoundary(i) is the face number in the array eval
    \param[in] LocalQuadPoints LocalQuadPoints(i) are the quadrature points on the face i
    \param[inout] eval eval(:)(IndexBoundary(i))(:) should contain the evaluation of the mode
                       on quadrature points of face j. 
   */
  template<class T, class Dimension>
  void ModalSourceBoundary_Dim<T, Dimension>
  ::EvaluateMode(const Mesh<Dimension>& mesh, int proc_mode,
		 const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
		 const Vector<int>& RotElement, const Vector<typename DimensionB::VectR_N>& LocalQuadPoints,
		 Vector<Vector<Vector<T> > >& eval)
  {
    cout << "Not implemented for any equation" << endl;
    abort();
  }

} // end namespace

#define MONTJOIE_FILE_MODAL_SOURCE_CXX
#endif
