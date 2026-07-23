#ifndef MONTJOIE_FILE_COMPUTATION_MODE_MAXWELL_CXX

#include "ComputationModeMaxwell.hxx"

namespace Montjoie
{

  /***********************
   * ModeMaxwellEquation *
   ***********************/
  
  
  //! computation of geometrical transformations that need to be computed
  //! for mass and stiffness matrix
  template<class T> template<class TypeEquation>
  void ModeMaxwellEquation<T>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                    int num_elem, const ElementReference_Dim<Dimension2>& Fb) 
  {
    try 
      {
	const QuadrangleHcurlFirstFamily& Fb_quad
	  = dynamic_cast<const QuadrangleHcurlFirstFamily& >(Fb);
	
	ComputeMassMatrixQuad(var, num_elem, Fb_quad);
	return;
      }
    catch (const std::bad_cast&)
      {
	// other element, this current method is used
      }
    
    Matrix2_2 dfj, dfjm1;
    Matrix2_2sym A_tmp; 
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;
      
    var.Glob_matMass_Bh(num_elem).Reallocate(N);
    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    var.Glob_matMassPoids(num_elem).Reallocate(N);
    
    for (int j = 0; j < N; j++)
      {
	// we compute and store -\omega^2 weight J_i DF_i^{-1} DF_i^{*-1}
	// physical media is supposed isotrope,
        // so that we only need to multiply by scalar values epsilon mu
        Real_wp jacob, poids = Fb.WeightsND(j);
        if (N == 1)
          poids = 1.0;
        
        if (variable)
          {
            if (affine)
              {
                jacob = var.Glob_jacobian(num_elem)(0);
                dfjm1 = var.Glob_DFjm1(num_elem)(0);                
                Mlt(1.0/jacob, dfjm1);
              }
            else
              {
                jacob = var.Glob_jacobian(num_elem)(j)/Fb.WeightsND(j);
                dfjm1 = var.Glob_DFjm1(num_elem)(j);
                Mlt(1.0/jacob, dfjm1);
              }
          }
        else
          {
            jacob = var.Glob_jacobian(num_elem)(0);
            dfjm1 = var.Glob_DFjm1(num_elem)(0);
            Mlt(1.0/jacob, dfjm1);
          }
        
        MltTrans(dfjm1, dfjm1, A_tmp);
        
	Real_wp coef = -var.GetSquareOmega()*jacob*poids;
	Mlt(coef, A_tmp);
	var.Glob_matMass_Bh(num_elem)(j) = A_tmp;
        
        var.Glob_matMass_Dh(num_elem)(j) = poids/jacob;
        var.Glob_matMassPoids(num_elem)(j) = poids*var.GetOmega();        
      }    
  }
  
  
  //! computation of geometrical transformations that need to be computed 
  //! for mass and stiffness matrix
  template<class T> template<class TypeEquation>
  void ModeMaxwellEquation<T>::
  ComputeMassMatrixQuad(EllipticProblem<TypeEquation>& var,
			int num_elem, const QuadrangleHcurlFirstFamily& Fb) 
  {
    Matrix2_2 dfj, dfjm1;
    Matrix2_2sym A_tmp; 
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;
      
    var.Glob_matMass_Bh(num_elem).Reallocate(N);
    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    var.Glob_matMassPoids(num_elem).Reallocate(N);
    
    for (int j = 0; j < N; j++)
      {
	// we compute and store -\omega^2 weight J_i DF_i^{-1} DF_i^{*-1}
	// physical media is supposed isotrope,
        // so that we only need to multiply by scalar values epsilon mu
        Real_wp jacob, poids = Fb.WeightsND(j);
        if (N == 1)
          poids = 1.0;
        
        if (variable)
          {
            if (affine)
              {
                jacob = var.Glob_jacobian(num_elem)(0);
                dfjm1 = var.Glob_DFjm1(num_elem)(0);                
                Mlt(1.0/jacob, dfjm1);
              }
            else
              {
                jacob = var.Glob_jacobian(num_elem)(j)/Fb.WeightsND(j);
                dfjm1 = var.Glob_DFjm1(num_elem)(j);
                Mlt(1.0/jacob, dfjm1);
              }
          }
        else
          {
            jacob = var.Glob_jacobian(num_elem)(0);
            dfjm1 = var.Glob_DFjm1(num_elem)(0);
            Mlt(1.0/jacob, dfjm1);
          }
        
        MltTrans(dfjm1, dfjm1, A_tmp);
        
	Real_wp coef = -var.GetSquareOmega()*jacob*poids;
	Mlt(coef, A_tmp);
	var.Glob_matMass_Bh(num_elem)(j) = A_tmp;
        
        var.Glob_matMassPoids(num_elem)(j) = poids*var.GetOmega();        
      } 
    
    VectR2 s;
    var.mesh.GetVerticesElement(num_elem, s);
    
    SetPoints<Dimension2> Points;
    Fb.FjElemQuadrature(s, Points, var.mesh, num_elem);
    
    SetMatrices<Dimension2> Mat;
    Fb.DFjElemQuadrature(s, Points, Mat, var.mesh, num_elem);
    
    int Ni = Fb.GetNbPointsQuadratureInside();
    N = min(N, Fb.GetOrder()*Fb.GetOrder());
    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    Real_wp Astiff;
    for (int j = 0; j < N; j++)
      {
	Real_wp jacob = Det(Mat.GetPointQuadrature(Ni+j));
	Real_wp poids_mass = Fb.WeightsGauss2D(j);
	if (N == 1)
	  poids_mass = 1.0;
	
        Astiff = poids_mass/jacob;

	var.Glob_matMass_Dh(num_elem)(j) = Astiff;
      }
  }

  
  /******************
   * VarModeMaxwell *
   ******************/
  
  
  //! return precomputed \f$ -\omega^2 \mbox{ weight } J_i DF_i^{-1} DF_i^{*-1} \f$
  template<class TypeEquation> template<class T0, class Prop>
  void VarModeMaxwell<TypeEquation>
  ::GetMassMatrixCoef(int num_elem, int j, const GlobalGenericMatrix<T0>& nat_mat,
		      int m, int n, TinyMatrix<T0, Prop, 2, 2>& mass) const
  {
    int ref = this->mesh.Element(num_elem).GetReference();
    if (m == 0)
      {
        if (n == 0)
          {
            mass = this->Glob_matMass_Bh(num_elem)(j);
            T0 coef;
            to_complex(this->ref_epsilon(ref).GetConstant()
                       + Iwp * this->ref_sigma(ref).GetConstant()/ this->GetOmega(), coef);
            Mlt(coef*nat_mat.GetCoefStiffness(), mass);
          }
        else
          {
            mass(0, 0) = 0;
            mass(1, 1) = 0;
            mass(0, 1) = this->Glob_matMassPoids(num_elem)(j)*nat_mat.GetCoefMass();
            mass(1, 0) = -this->Glob_matMassPoids(num_elem)(j)*nat_mat.GetCoefMass();
          }
      }
    else
      {
        if (n == 1)
          {
            mass = this->Glob_matMass_Bh(num_elem)(j);
            Mlt(this->ref_mu(ref).GetConstant()*nat_mat.GetCoefStiffness(), mass);
          }
        else
          {
            mass(0, 0) = 0;
            mass(1, 1) = 0;
            mass(1, 0) = -this->Glob_matMassPoids(num_elem)(j)*nat_mat.GetCoefMass();
            mass(0, 1) = this->Glob_matMassPoids(num_elem)(j)*nat_mat.GetCoefMass();
          }
      }
  }
    
  
  //! return coefficient weight/Ji
  template<class TypeEquation> template<class T0>
  void VarModeMaxwell<TypeEquation>::
  GetStiffMatrixCoef(int num_elem, int j,
		     const GlobalGenericMatrix<T0>& nat_mat, int m, int n, T0& stiff) const
  {
    int ref = this->mesh.Element(num_elem).GetReference();
    stiff = 0;
    if (j >= this->Glob_matMass_Dh(num_elem).GetM())
      return;
    
    if (m == n)
      {
        if (m == 0)
          stiff = this->Glob_matMass_Dh(num_elem)(j)/this->ref_mu(ref).GetConstant();
        else
          {
            T0 coef;
            to_complex(this->ref_epsilon(ref).GetConstant()
                       + Iwp * this->ref_sigma(ref).GetConstant()/ this->GetOmega(), coef);

            stiff = this->Glob_matMass_Dh(num_elem)(j)/coef;
          }
        
        stiff *= nat_mat.GetCoefStiffness();
      }
  }
  

  //! initialisation of physical indices
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>::InitIndices(int n)
  {
    ref_epsilon.Reallocate(n+1);
    ref_mu.Reallocate(n+1);
    ref_sigma.Reallocate(n+1);
    for (int i = 0; i <= n; i++)
      {
	ref_epsilon(i).SetIdentity();
	ref_mu(i).SetIdentity();
        ref_sigma(i).Zero();
      }
  }
  
  
  //! returns the number of physical media
  template<class TypeEquation>
  int VarModeMaxwell<TypeEquation>::GetNbPhysicalIndices() const
  {
    return ref_epsilon.GetM();
  }
  
  
  //! filling of indices arrays depending the data file
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_epsilon.GetM())
      {
	cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
	cout << "Current reference i : " << i << " < " << ref_epsilon.GetM() << endl;
	abort();
      }
    
    ref_epsilon(i).SetInputData(nb, parameters, parameters(0));
    ref_mu(i).SetInputData(nb, parameters, parameters(0));
    ref_sigma(i).SetInputData(nb, parameters, parameters(0));
  }
  
  
  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
  */
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_epsilon.GetM())
      {
	cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
	cout << "Current reference i : " << i << " < " << ref_epsilon.GetM() << endl;
	abort();
      }
    
    if (name_media == "epsilon")
      ref_epsilon(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "mu")
      ref_mu(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "sigma")
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    else
      {
	cout << "Unknown media : " << name_media << endl;
	abort();
      }
  }
  
  
  //! returns the name associated with the physical index num
  template<class TypeEquation>
  string VarModeMaxwell<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("epsilon");
      case 1: return string("mu");
      case 2: return string("sigma");
      }

    return string();
  }  
    
    
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>
  ::UpdateShiftAdimensionalization(Real_wp& shift, Real_wp& shift_imag)
  {
    // we recompose beta from shift
    Real_wp beta = shift;
    if (!IsComplexNumber(Complexe(0)))
      beta = sqrt(abs(shift))/(2.0*pi_wp);

    DISP(beta);
    // then we multiply by omega
    shift = beta*this->GetOmega();
  }


  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>
  ::UpdateShiftAdimensionalization(Complex_wp& shift, Complex_wp& shift_imag)
  {
    // we recompose beta from shift
    Complex_wp beta = shift;
    if (!IsComplexNumber(Complexe(0)))
      beta = sqrt(abs(shift))/(2.0*pi_wp);

    DISP(beta);
    // then we multiply by omega
    shift = beta*this->GetOmega();
  }
  
  
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>
  ::UpdateEigenvaluesAdimensionalization(VectReal_wp& eigen_values, VectReal_wp& lambda_imag,
					 Matrix<Real_wp, General, ColMajor>& eigen_vectors)
  {
    // eigenvalues are divided by omega in order to obtain the effective index
    for (int i = 0; i < eigen_values.GetM(); i++)
      {
	eigen_values(i) /= this->GetOmega();
	lambda_imag(i) /= this->GetOmega();
      }
  }


  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>
  ::UpdateEigenvaluesAdimensionalization(VectComplex_wp& eigen_values, VectComplex_wp& lambda_imag,
					 Matrix<Complex_wp, General, ColMajor>& eigen_vectors)
  {
    // eigenvalues are divided by omega in order to obtain the effective index
    for (int i = 0; i < eigen_values.GetM(); i++)
      {
	eigen_values(i) /= this->GetOmega();
	lambda_imag(i) /= this->GetOmega();
      }
  }
  
  
  //! not implemented
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
		      IVect& num_ref, IVect& num_index, IVect& num_component,
		      Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
  }
  

  template<class TypeEquation>
  bool VarModeMaxwell<TypeEquation>::IsVaryingMedia(int i) const
  {
    return false;
  }


  template<class TypeEquation>
  bool VarModeMaxwell<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return false;
  }

  
  //! allocation of needed arrays for the computation of elementary matrices
  template<class TypeEquation>
  void VarModeMaxwell<TypeEquation>::AllocateMassMatrices()
  {
    this->Glob_matMass_Bh.Reallocate(this->mesh.GetNbElt());
    this->Glob_matMass_Dh.Reallocate(this->mesh.GetNbElt());
    this->Glob_matMassPoids.Reallocate(this->mesh.GetNbElt());
  }


  /*******************
   * EllipticProblem *
   *******************/
  
  
  void EllipticProblem<ModeMaxwellEquation<Real_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementHcurl(i));
  }


  void EllipticProblem<ModeMaxwellEquation<Real_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementHcurl(i));
  }

  void EllipticProblem<ModeMaxwellEquation<Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible for complex numbers (ModeMaxwellEquation)" << endl;
    abort();
  }


  void EllipticProblem<ModeMaxwellEquation<Complex_wp> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementHcurl(i));
  }


  /*****************************
   * ModalSourceMaxwell3D_Base *
   *****************************/


#ifdef MONTJOIE_WITH_THREE_DIM
  void ModalSourceMaxwell3D_Base
  ::GetCoefficientBoundaryCondition(const VectReal_wp& eigenval, VectReal_wp& coef_eigenval,
				    const Real_wp& impedance, int boundary_condition)
  {
    // no modification for real numbers
    coef_eigenval.Reallocate(eigenval.GetM());
    coef_eigenval.Fill(1.0);
  }


  void ModalSourceMaxwell3D_Base
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
  

  void ModalSourceMaxwell3D_Base
  ::EvaluateModePoisson(const Mesh<Dimension3>& result_mesh, int proc_mode,
			const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
			const IVect& RotElement, const Vector<VectR2>& LocalQuadPoints,
			Vector<VectR3>& eval_mode)
  {
    // object used to solve Poisson equation
    EllipticProblem<LaplaceEquation<Dimension2> > var_poisson;    
    All_LinearSolver glob_solver(var_poisson);
    IVect RotationFaceElement = RotElement;

    // resolution parameters
    glob_solver.SetDirectSolver();
    
    // physical properties and boundary conditions
    var_poisson.InitIndices(var_problem.GetNbPhysicalIndices()-1);
    var_poisson.mesh.ResizeNbReferences(var_problem.mesh.GetNbReferences());
    var_poisson.mesh.SetBoundaryCondition(var_problem.mesh);
    this->CopyIndices(var_poisson);
    
    // mesh and finite element are constructed
    R3 vec_u, vec_v;
    modal_source.ConstructProblem(result_mesh, var_problem, var_poisson, RotationFaceElement,
				  string("TRIANGLE_LOBATTO"), vec_u, vec_v);

    // projectors are computed
    Vector<FiniteElementInterpolator> proj;
    TinyVector<Vector<Matrix<int> >, 2> RotationPoints;
    
    modal_source.ComputeProjector(var_problem, var_poisson, IndexBoundary, IndexElement,
				  LocalQuadPoints, proj, RotationPoints);
    
    // factorisation of finite element matrix
    GlobalGenericMatrix<Real_wp> nat_mat;
    var_poisson.ComputeMassMatrix();
    glob_solver.PerformFactorizationStep(nat_mat);

    MeshNumbering<Dimension2>& mesh_num_poisson = var_poisson.GetMeshNumbering(0);
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    
    // inhomogeneous Dirichlet condition
    VectReal_wp rhs(var_poisson.GetNbDof());
    rhs.Zero();
    for (int i = 0; i < var_poisson.mesh.GetNbBoundaryRef(); i++)
      {
	int ref = var_poisson.mesh.BoundaryRef(i).GetReference();
	int p = -1;
	for (int j = 0; j < ref_mode_poisson.GetM(); j++)
	  if (ref_mode_poisson(j) == ref)
	    {
	      p = j;
	      break;
	    }
	
	if (p >= 0)
	  {
            int num_face = i;
	    int num_elem = var_poisson.mesh.BoundaryRef(i).numElement(0);
	    int num_loc = var_poisson.mesh.Element(num_elem).GetPositionBoundary(num_face);
	    Real_wp coef = coef_mode_poisson(p);

	    const ElementReference_Dim<Dimension2>& Fb = var_poisson.GetReferenceElement(num_elem);
	    int nb_dof = Fb.GetNbDofBoundary(num_loc);
	    for (int j = 0; j < nb_dof; j++)
	      { 
		int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
		int num_dof = mesh_num_poisson.Element(num_elem).GetNumberDof(num_dof_loc);
		rhs(num_dof) = coef;
	      }
	  }
      }
    
    // Poisson problem is solved
    VectReal_wp mode_dof(var_poisson.GetNbDof());
    mode_dof = rhs;
    glob_solver.ComputeSolution(mode_dof, nat_mat);

    // uncomment the following lines to visualize the mode
    
    /*Vector<string> param(2);
    param(0) = string("DiffracModeLaplace");
    param(1) = string("ModeLaplace");
    var_poisson.SetInputData(string("FileOutputPlane"), param);
    
    param.Reallocate(3);
    param(0) = "AUTO";
    param(1) = "201";
    param(2) = "201";
    var_poisson.SetInputData(string("SismoPlane"), param);
    
    var_poisson.InitVarGrid();
    var_poisson.ComputeVarGrid();
    
    var_poisson.WriteDatas(mode_dof); */

    // the mode is interpolated on quadrature points
    int nb_faces = IndexBoundary.GetM();
    eval_mode.Reallocate(nb_faces);
    VectReal_wp u_loc, u_quad; VectR2 s, grad_Uloc, grad_Uquad;
    SetPoints<Dimension2> PointsElem; SetMatrices<Dimension2> MatricesElem;
    for (int i = 0; i < nb_faces; i++)
      {
	// we retrieve value of u on the element
	int nb_dof_loc = var_poisson.GetNbLocalDof(i);
	u_loc.Reallocate(nb_dof_loc);
	for (int j = 0; j < nb_dof_loc; j++)
	  u_loc(j) = mode_dof(mesh_num_poisson.Element(i).GetNumberDof(j));

	// gradient of u is computed
	const ElementGeomReference<Dimension2>& Fb = var_poisson.GetShapeElement(i);
	var_poisson.mesh.GetVerticesElement(i, s);
	Fb.FjElemNodal(s, PointsElem, var_poisson.mesh, i);
	Fb.DFjElemNodal(s, PointsElem, MatricesElem, var_poisson.mesh, i);
	Fb.ComputeNodalGradient(MatricesElem, u_loc, grad_Uloc);
	
	// projection on quadrature points
	int num_elem = IndexElement(i);
	int num_loc = IndexBoundary(i);
	int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	int r = mesh_num.GetOrderQuadrature(num_face);
	int type_elt = var_problem.mesh.Boundary(num_face).GetHybridType();	    
	int nb_pts_quad = LocalQuadPoints(i).GetM();
	grad_Uquad.Reallocate(nb_pts_quad);
	proj(r).Project(grad_Uloc, grad_Uquad, type_elt);
	
	// eval_mode is filled
	eval_mode(i).Reallocate(nb_pts_quad);
	for (int j = 0; j < nb_pts_quad; j++)
	  {
	    const Matrix<int>& RotationNum = RotationPoints(type_elt)(r);		
	    int jrot = RotationNum(RotationFaceElement(i), j);
	    Add(grad_Uquad(jrot)(0), vec_u, eval_mode(i)(j));
	    Add(grad_Uquad(jrot)(1), vec_v, eval_mode(i)(j));
	  }
      }
  }


  void ModalSourceMaxwell3D_Base
  ::EvaluateModeMaxwell(const Mesh<Dimension3>& result_mesh, int proc_mode,
			const Vector<int>& IndexBoundary, const Vector<int>& IndexElement,
			const IVect& RotElement, const Vector<VectR2>& LocalQuadPoints,
			Vector<Vector<VectR3> >& eval_mode, VectReal_wp& eigenval)
  {
    // object used to compute eigenvalues
    EllipticProblem<ModeMaxwellEquation<Real_wp> > var_eig;    
    All_LinearSolver glob_solver(var_eig);
    IVect RotationFaceElement = RotElement;
    
    // resolution parameters
    glob_solver.SetDirectSolver();
    
    // physical properties and boundary conditions
    var_eig.InitIndices(var_problem.GetNbPhysicalIndices()-1);
    var_eig.SetOmega(var_problem.GetOmega());
    var_eig.mesh.ResizeNbReferences(var_problem.mesh.GetNbReferences());
    var_eig.mesh.SetBoundaryCondition(var_problem.mesh);
    var_eig.SetStorageFiniteElementMatrix(var_eig.MATRIX_STORED);
    this->CopyIndices(var_eig);
    
    // mesh and finite element are constructed
    R3 vec_u, vec_v;
    modal_source.ConstructProblem(result_mesh, var_problem, var_eig, RotationFaceElement,
				  string("TRIANGLE_FIRST_FAMILY"), vec_u, vec_v);

    // eigenvalues and eigenvectors are computed
    EigenProblemMontjoie<Real_wp> eigen_solver(var_eig, glob_solver);
    VectReal_wp eigen_values;
    Matrix<Real_wp, General, ColMajor> eigen_vectors;    

    modal_source.ComputeEigenvalues(eigen_solver, var_eig, eigen_values, eigen_vectors,
				    eigen_solver.INVERT_MODE, var_eig.GetOmega(), false);
    
    // projectors are computed
    Vector<FiniteElementInterpolator> proj;
    TinyVector<Vector<Matrix<int> >, 2> RotationPoints;
    
    modal_source.ComputeProjector(var_problem, var_eig, IndexBoundary, IndexElement,
				  LocalQuadPoints, proj, RotationPoints);

    // we store the required modes
    VectReal_wp mode_dof(var_eig.GetNbDof());
    VectReal_wp Xmode(eigen_vectors.GetM());
    const IVect& number_mode_combined = modal_source.GetModeNumberToCombine();
    eigenval.Reallocate(number_mode_combined.GetM());
    eval_mode.Reallocate(number_mode_combined.GetM());
    MeshNumbering<Dimension2>& mesh_num_eig = var_eig.GetMeshNumbering(0);
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    for (int k = 0; k < number_mode_combined.GetM(); k++)
      {
	int number_mode = number_mode_combined(k);
	eigenval(k) = eigen_values(number_mode);
	for (int i = 0; i < eigen_vectors.GetM(); i++)
	  Xmode(i) = eigen_vectors(i, number_mode); 

        eigen_solver.ExpandVector(Xmode, mode_dof, true);
        
	// uncomment the following lines to visualize the mode
	
	/*Vector<string> param(2);
	param(0) = string("DiffracModeMaxwell");
	param(1) = string("ModeMaxwell");
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
	Vector<VectReal_wp> u_loc(2), uloc_node(2);
	VectR2 E_loc, E_quad, s;
	SetPoints<Dimension2> PointsElem; SetMatrices<Dimension2> MatricesElem;
	int Nvol = mesh_num_eig.GetNbDof();
	for (int i = 0; i < nb_faces; i++)
	  {
	    // we retrieve value of u on the element
	    int nb_dof_loc = var_eig.GetNbLocalDof(i);
	    u_loc(0).Reallocate(nb_dof_loc); u_loc(1).Reallocate(nb_dof_loc);
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		int num_dof = mesh_num_eig.Element(i).GetNumberDof(j);
		u_loc(0)(j) = mode_dof(num_dof);
		u_loc(1)(j) = mode_dof(Nvol + num_dof);
	      }
	    
	    mesh_num_eig.number_map.
	      ModifyLocalComponentVector(mesh_num_eig, u_loc(0), i);

	    mesh_num_eig.number_map.
	      ModifyLocalComponentVector(mesh_num_eig, u_loc(1), i);
	    
	    // values on nodal points are computed
	    const ElementReference<Dimension2, 2>& Fb = var_eig.GetReferenceElementHcurl(i);
	    var_eig.mesh.GetVerticesElement(i, s);
	    Fb.FjElemNodal(s, PointsElem, var_eig.mesh, i);
	    Fb.DFjElemNodal(s, PointsElem, MatricesElem, var_eig.mesh, i);
	    
	    Fb.ComputeNodalValues(MatricesElem, u_loc, uloc_node,
				  var_eig.mesh, i);
	    
	    int nb_nodes = uloc_node(0).GetM();
	    E_loc.Reallocate(nb_nodes);
	    for (int j = 0; j < nb_nodes; j++)
	      E_loc(j).Init(uloc_node(0)(j), uloc_node(1)(j));
	    
	    // projection on quadrature points
	    int num_elem = IndexElement(i);
	    int num_loc = IndexBoundary(i);
	    int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	    int r = mesh_num.GetOrderQuadrature(num_face);
	    int type_elt = var_problem.mesh.Boundary(num_face).GetHybridType();	    
	    int nb_pts_quad = LocalQuadPoints(i).GetM();
	    E_quad.Reallocate(nb_pts_quad);
	    proj(r).Project(E_loc, E_quad, type_elt);
	    
	    // eval_mode is filled
	    eval_mode(k)(i).Reallocate(nb_pts_quad);
	    for (int j = 0; j < nb_pts_quad; j++)
	      {
		const Matrix<int>& RotationNum = RotationPoints(type_elt)(r);
		int jrot = RotationNum(RotationFaceElement(i), j);
		Add(E_quad(jrot)(0), vec_u, eval_mode(k)(i)(j));
		Add(E_quad(jrot)(1), vec_v, eval_mode(k)(i)(j));
	      }
	  }
      }
  }


  /************************
   * ModalSourceMaxwell3D *
   ************************/


  template<class T>
  void ModalSourceMaxwell3D<T>
  ::ReadSourceParameters(const IVect& boundary_condition, const IVect& ref,
			 const VectString& source_param)
  {
    // is_ref_on_mode(ref) will be 1 if ref belongs to the mode
    this->is_ref_on_mode.Fill(0);
    for (int k = 0; k < ref.GetM(); k++)
      this->is_ref_on_mode(ref(k)) = 1;

    // keyword to know the type of mode (POISSON, MODE_MAXWELL)
    int nb = 0;
    if (source_param.GetM() < 6)
      {
	cout << "Provide more parameters for the definition of the mode" << endl;
	cout << "Current parameters are " << endl << source_param << endl;
	abort();
      }

    this->keyword_mode = source_param(nb++);
    if (!this->keyword_mode.compare("POISSON"))
      {
	this->mode_to_compute = MODE_POISSON;
	int nb_pot = to_num<int>(source_param(nb++));
	this->ref_mode_poisson.Reallocate(nb_pot);
	this->coef_mode_poisson.Reallocate(nb_pot);

	if (source_param.GetM() < nb+2*nb_pot+1)
	  {
	    cout << "Provide more parameters for the definition of the mode" << endl;
	    cout << "Current parameters are " << endl << source_param << endl;
	    abort();
	  }
	
	for (int k = 0; k < nb_pot; k++)
	  {
	    this->ref_mode_poisson(k) = to_num<int>(source_param(nb++));
	    this->coef_mode_poisson(k) = to_num<Real_wp>(source_param(nb++));
	  }
	
	this->voltage = to_num<T>(source_param(nb++));
      }
    else if (!this->keyword_mode.compare("MAXWELL"))
      {
	this->mode_to_compute = MODE_MAXWELL;
	this->shift_eigenval_mode = to_num<Real_wp>(source_param(nb++));		    
	this->nb_modes_to_compute = to_num<int>(source_param(nb++));		    
	this->nb_modes_to_combine = to_num<int>(source_param(nb++));
	if (source_param.GetM() < nb+2*this->nb_modes_to_combine)
	  {
	    cout << "Provide more parameters for the definition of the mode" << endl;
	    cout << "Current parameters are " << endl << source_param << endl;
	    abort();
	  }
	
	this->number_mode_combined.Reallocate(this->nb_modes_to_combine);
	this->coef_mode_combined.Reallocate(this->nb_modes_to_combine);
	for (int k = 0; k < this->nb_modes_to_combine; k++)
	  {
	    this->number_mode_combined(k) = to_num<int>(source_param(nb++));
	    this->coef_mode_combined(k) = to_num<T>(source_param(nb++));
	  }
      }
  }  
  
  
  //! evaluates the mode on quadrature points of the mesh
  template<class T>
  void ModalSourceMaxwell3D<T>
  ::EvaluateMode(const Mesh<Dimension3>& mesh, int proc_mode, const Vector<int>& IndexBoundary, 
		 const Vector<int>& IndexElement, const Vector<int>& RotElement,
		 const Vector<VectR2>& LocalQuadPoints,
		 Vector<Vector<Vector<T> > >& eval)
  {
    // evaluates the required modes (considered as real-valued)
    Vector<Vector<VectR3> > eval_mode; VectReal_wp eigenval;
    if (this->mode_to_compute == MODE_POISSON)
      {
	eval_mode.Reallocate(1); eigenval.Reallocate(1);
	eigenval(0) = Real_wp(0);
	this->EvaluateModePoisson(mesh, proc_mode, IndexBoundary, IndexElement,
				  RotElement, LocalQuadPoints, eval_mode(0));
      }
    else
      this->EvaluateModeMaxwell(mesh, proc_mode, IndexBoundary, IndexElement,
				RotElement, LocalQuadPoints, eval_mode, eigenval);
    
    // we get coefficients due to boundary condition
    Vector<T> coef_eigenval;
    int ref = this->GetVolumeReference();
    this->GetCoefficientBoundaryCondition(eigenval, coef_eigenval,
					  var_maxwell.coefficient_impedance_absorbing(ref),
					  this->boundary_condition_mode);
    
    // multiplication by coefficients provided in the data file
    if (this->mode_to_compute == MODE_POISSON)
      coef_eigenval(0) *= voltage;
    else
      for (int k = 0; k < eval_mode.GetM(); k++)
	coef_eigenval(k) *= this->coef_mode_combined(k);
    
    // filling eval by combining modes stored in eval_mode
    for (int i = 0; i < LocalQuadPoints.GetM(); i++)
      {
	int num_elem = IndexElement(i);
	int num_loc = IndexBoundary(i);
	int num_face = var_problem.mesh.Element(num_elem).numBoundary(num_loc);
	int num_face_ref = num_face;
	int nb_points = LocalQuadPoints(i).GetM();
	eval(0)(num_face_ref).Reallocate(nb_points);
	eval(1)(num_face_ref).Reallocate(nb_points);
	eval(2)(num_face_ref).Reallocate(nb_points);
	for (int j = 0; j < nb_points; j++)
	  {
	    SetComplexZero(eval(0)(num_face_ref)(j));
	    SetComplexZero(eval(1)(num_face_ref)(j));
	    SetComplexZero(eval(2)(num_face_ref)(j));
	    for (int k = 0; k < eval_mode.GetM(); k++)
	      {
		eval(0)(num_face_ref)(j) += coef_eigenval(k)*eval_mode(k)(i)(j)(0);
		eval(1)(num_face_ref)(j) += coef_eigenval(k)*eval_mode(k)(i)(j)(1);
		eval(2)(num_face_ref)(j) += coef_eigenval(k)*eval_mode(k)(i)(j)(2);
	      }
	  }
      }
  }
  
  
  //! Sets indices of a Laplace problem with indices of maxwell's equations
  template<class T>
  void ModalSourceMaxwell3D<T>
  ::CopyIndices(VarHelmholtz_Base<Real_wp, Dimension2>& var_eig)
  {
    Real_wp zero(0);
    TinyMatrix<Real_wp, Symmetric, 2, 2> mu;
    for (int ref = 0; ref < var_problem.GetNbPhysicalIndices(); ref++)
      {
	mu.SetIdentity();
	mu *= realpart(var_maxwell.ref_epsilon(ref).GetConstant()(0, 0));
	var_eig.ref_rho(ref).SetConstant(zero);
	var_eig.ref_sigma(ref).SetConstant(zero);
	var_eig.ref_mu(ref).SetConstant(mu);
      }
  }


  //! Sets indices of a mode problem with indices of maxwell's equations
  template<class T>
  void ModalSourceMaxwell3D<T>
  ::CopyIndices(VarModeMaxwell<ModeMaxwellEquation<Real_wp> >& var_eig)
  {
    Real_wp mu, epsilon;
    for (int ref = 0; ref < var_problem.GetNbPhysicalIndices(); ref++)
      {
	epsilon = realpart(var_maxwell.ref_epsilon(ref).GetConstant()(0, 0));
	mu = realpart(var_maxwell.ref_mu(ref).GetConstant()(0, 0));
	var_eig.ref_epsilon(ref).SetConstant(epsilon);
	var_eig.ref_mu(ref).SetConstant(mu);
      }
    
    IVect comp(1);
    comp(0) = 0;    
    for (int ref = 1; ref <= var_eig.mesh.GetNbReferences(); ref++)
      if (var_eig.mesh.GetBoundaryCondition(ref) > 0)
	{
          var_eig.mesh.SetBoundaryCondition(ref, BoundaryConditionEnum::LINE_SUPPORTED);
          var_eig.SetSupportedComponents(ref, comp);
        }
  }
#endif
  
};

#define MONTJOIE_FILE_COMPUTATION_MODE_MAXWELL_CXX
#endif
