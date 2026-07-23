#ifndef MONTJOIE_FILE_STATIC_MAXWELL_2D_CXX

namespace Montjoie
{
  
  
  /****************************
   * StaticMaxwellEquation_2D *
   ****************************/
  
  
  //! computation of mass coefficients
  template<class TypeEquation>
  void StaticMaxwellEquation_2D
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
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
    
    int ref = var.mesh.Element(num_elem).GetReference();
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;
    
    var.Glob_matMass_Bh(num_elem).Reallocate(N);
    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    
    Real_wp mu(0), invMu(0);
    Matrix2_2sym epsilon;
    
    for (int j = 0; j < N; j++)
      {
	epsilon = var.ref_epsilon(ref).GetCoefficient(var, num_elem, j);
	mu = var.ref_mu(ref).GetCoefficient(var, num_elem, j);            
	
        invMu = 1.0/mu;	    
	
	Real_wp poids = Fb.WeightsND(j);
        if (N == 1)
          poids = 1.0;
        
	Real_wp jacob; Matrix2_2 dfj, dfjm1; Matrix2_2 A_tmp;
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
	
	GetInverse(dfjm1, dfj);
	
	Matrix2_2sym Bmass; Real_wp Astiff;
		
	MltTrans(epsilon, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Bmass); Mlt(jacob*poids, Bmass);
	
	if (var.FirstOrderFormulation())
	  Astiff = poids*jacob*mu;
	else
	  Astiff = poids/jacob*invMu;
	
	var.Glob_matMass_Bh(num_elem)(j) = Bmass;
	var.Glob_matMass_Dh(num_elem)(j) = Astiff;
      }
  }

#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! computation of mass coefficients
  template<class TypeEquation>
  void StaticMaxwellEquation_2D
  ::ComputeMassMatrixQuad(EllipticProblem<TypeEquation>& var,
			  int num_elem, const QuadrangleHcurlFirstFamily& Fb)
  {
    int ref = var.mesh.Element(num_elem).GetReference();
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;
    
    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    var.Glob_matMass_Bh(num_elem).Reallocate(N);
    
    Real_wp mu(0), invMu(0);
    Matrix2_2sym epsilon;
    
    for (int j = 0; j < N; j++)
      {
	epsilon = var.ref_epsilon(ref).GetCoefficient(var, num_elem, j);
	mu = var.ref_mu(ref).GetCoefficient(var, num_elem, j);            
	
        invMu = 1.0/mu;
	Real_wp poids = Fb.WeightsND(j);
        if (N == 1)
          poids = 1.0;
        
	Matrix2_2 dfj, dfjm1; Matrix2_2 A_tmp;
	Real_wp jacob;
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
	
	GetInverse(dfjm1, dfj);
	
	Matrix2_2sym Bmass;
        
	MltTrans(epsilon, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Bmass); Mlt(jacob*poids, Bmass);
	
	var.Glob_matMass_Bh(num_elem)(j) = Bmass;
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
    for (int j = 0; j < N; j++)
      {
	invMu = var.ref_invMu(ref).GetCoefficient(var, num_elem, Ni+j);
        
	Real_wp jacob = Det(Mat.GetPointQuadrature(Ni+j));
	Real_wp poids_mass = Fb.WeightsGauss2D(j);
	if (N == 1)
	  poids_mass = 1.0;
	
	Real_wp Astiff;
	if (var.FirstOrderFormulation())
	  Astiff = poids_mass*jacob*mu;
	else
	  Astiff = poids_mass/jacob*invMu;
	
	var.Glob_matMass_Dh(num_elem)(j) = Astiff;
      }
  }
#endif

  
  /*******************************
   * StaticMaxwellEquation_2D_DG *
   *******************************/

  
  //! returns the matrix M, in the integral \f$ \int_K M u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] mass matrix M
  */    
  template<class TypeEquation, class T0, class MatMass>
  void StaticMaxwellEquation_2D_DG::
  GetTensorMass(const EllipticProblem<TypeEquation>& var,
		int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    FillZero(mass);
    mass(0, 0) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 0);
    mass(0, 1) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 1);
    mass(1, 0) = mass(0, 1);
    mass(1, 1) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 1, 1);
    mass(2, 2) = -var.ref_mu(ref).GetCoefficient(var, i, j);
    mass *= nat_mat.GetCoefMass();
  }
  
  
  //! Applies matrix M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un vector to be multiplied by M
    \param[out] Vn result vector Vn = M Un
  */    
  template<class TypeEquation, class T0, class Vector1>
  void StaticMaxwellEquation_2D_DG::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    V(0) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 0)*U(0);
    V(0) += var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 1)*U(1);
    V(1) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 1)*U(0);
    V(1) += var.ref_epsilon(ref).GetCoefficient(var, i, j, 1, 1)*U(1);
    V(2) = -var.ref_mu(ref).GetCoefficient(var, i, j)*U(2);
    V *= nat_mat.GetCoefMass();
  }
  

  //! which derivatives to evaluate during matrix-vector product ?  
  template<class TypeEquation, class T0, class Vector1>
  void StaticMaxwellEquation_2D_DG::
  GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
		      const GlobalGenericMatrix<T0>& nat_mat,
		      Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    unknown_to_derive.Fill(false);
    fct_test_to_derive.Fill(false);
    unknown_to_derive(2) = true;
    fct_test_to_derive(2) = true;
  }
  

  //! fills tensors D and E appearing in the variational formulation
  /*!
    \param[in] vars considered problem
    \param[in] num_elem element where D and E must be computed
    \param[in] jloc quadrature point where D and E must be computed
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the element
    \param[out] Ephi_grad tensor E
    \param[out] Dphi_grad tensor D    
    The tensors D and E are appearing in the terms
    \int_K D \nabla u v + E u \nabla v dx
    of the variational formulation
   */  
  template<class TypeEquation, class T0, class MatStiff>
  void StaticMaxwellEquation_2D_DG::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    Ephi_grad(1,2)(0) = -nat_mat.GetCoefStiffness();
    Ephi_grad(0,2)(1) = nat_mat.GetCoefStiffness();

    Dgrad_phi(2,1)(0) = -nat_mat.GetCoefStiffness();
    Dgrad_phi(2,0)(1) = nat_mat.GetCoefStiffness();

  }
  

  //! Applying the tensor D to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Vn gradient of the unknown vector V
    \param[out] Un result to D grad(v)
  */  
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void StaticMaxwellEquation_2D_DG::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
		       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		       int ref, Vector1& Vn, Vector2& Un)
  {
    Un.Fill(0);
    Un(0) = Vn(2)(1);
    Un(1) = -Vn(2)(0);
    Un *= nat_mat.GetCoefStiffness();
  }
  
  
  //! Applying the tensor E to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un unknown vector U
    \param[out] Vn result E u
  */
  template<class TypeEquation, class T0, class Vector1, class Vector2>  
  void StaticMaxwellEquation_2D_DG::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
		       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		       int ref, Vector1& Un, Vector2& Vn)
  {
    FillZero(Vn);
    Vn(2)(0) = -Un(1)*nat_mat.GetCoefStiffness();
    Vn(2)(1) = Un(0)*nat_mat.GetCoefStiffness();
  }
  
  
  //! computation of matrix N associated to the boundary condition
  /*!
    \param[out] Nabc matrix N
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Matrix1, class R_N, class GenericPb, class T0>
  void StaticMaxwellEquation_2D_DG::
  GetNabc(Matrix1& Nabc, R_N& normale, int ref, int iquad, int k, 
	  const GlobalGenericMatrix<T0>& nat_mat, int ref2,
	  const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Nabc.Fill(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// condition Dirichlet
	Nabc(2,0) = normale(1);
	Nabc(2,1) = -normale(0);
	Nabc(0,2) = normale(1);
	Nabc(1,2) = -normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	// condition Neumann
	Nabc(2,0) = -normale(1);
	Nabc(2,1) = normale(0);
	Nabc(0,2) = -normale(1);
	Nabc(1,2) = normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	// int ref_domain = vars.mesh.Element(iquad).GetReference();
	Nabc(2,2) = 1.0;
	Nabc(0,0) = normale(1)*normale(1);
	Nabc(1,1) = normale(0)*normale(0);
	Nabc(0,1) = -normale(0)*normale(1);
	Nabc(1,0) = -normale(0)*normale(1);	  
      }
  }
  
  
  //! Computation of "penalization" matrices C
  /*!
    \param[out] Nabc penalization matrix C
    \param[in] normale outward normale
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
  */
  template<class Matrix1, class R_N, class GenericPb, class T0>
  void StaticMaxwellEquation_2D_DG::
  GetPenalDG(Matrix1& Nabc, R_N& normale, int iquad, int k, int nf,
	     const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	     const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    Nabc.Fill(0);
    Nabc(2,2) = vars.delta_penalization;
    Nabc(0,0) = vars.alpha_penalization;
    Nabc(1,1) = vars.alpha_penalization;
    // Nabc(0,0) = vars->alpha_penalization*normale(1)*normale(1);
    // Nabc(1,1) = vars->alpha_penalization*normale(0)*normale(0);
    // Nabc(0,1) = -vars->alpha_penalization*normale(0)*normale(1);
    // Nabc(1,0) = -vars->alpha_penalization*normale(0)*normale(1);
  }
  

  //! Multiplication by penalization matrices
  /*!
    \param[in] normale outward normale
    \param[in] Vn vector to multiply
    \param[out] Un result vector Un = C*Vn
    \param[in] iquad element number
    \param[in] k local quadrature point number    
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Vector1, class Vector2, class GenericPb, class T0>
  void StaticMaxwellEquation_2D_DG::
  MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un, int i, int k,
	     int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	     const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    Un(2) = vars.delta_penalization*Vn(2);
    Un(1) = vars.alpha_penalization*Vn(1);
    Un(0) = vars.alpha_penalization*Vn(0);
    // Complexe E_times_n = vars->alpha_penalization*(Un(0)*normale(1)-Un(1)*normale(0));
    // Un(0) = E_times_n*normale(1);
    // Un(1) = -E_times_n*normale(0);
  }
  
  
  //! mutliplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] num_elem1 element number
    \param[in] num_point local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Vector1, class TypeEquation, class T0>
  void StaticMaxwellEquation_2D_DG::
  MltNabc(const R2& normale, int ref, const Vector1& Vn, Vector1& Un,
	  int num_elem1, int npoint, const GlobalGenericMatrix<T0>& nat_mat, int ref2,
	  const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    typedef typename TypeEquation::Complexe Complexe;
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Un.Fill(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	Un(2) = Vn(0)*normale(1)-Vn(1)*normale(0);
	Un(0) = Vn(2)*normale(1);
	Un(1) = -Vn(2)*normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	Un(2) = -Vn(0)*normale(1)+Vn(1)*normale(0);
	Un(0) = -Vn(2)*normale(1);
	  Un(1) = Vn(2)*normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	Un(2) = Vn(2);
	Complexe E_times_n = Un(0)*normale(1)-Un(1)*normale(0);
	Un(0) = E_times_n*normale(1);
	Un(1) = -E_times_n*normale(0);
      }
    // DISP(Un); DISP(Vn);
  }    
  
  
  /****************
   * ImpedanceABC *
   ****************/
  

  //! impedance
  template<class T>
  void ImpedanceABC<T, StaticMaxwellEquation_2D>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                         const GlobalGenericMatrix<T>& nat_mat, int ref,
                         const SetPoints<Dimension2>& Pts, const SetMatrices<Dimension2>& Mat)
  {
    this->read_param_condition = false;
    this->coef_grad = 0;
    if (var_maxwell.coefficient_impedance_absorbing.GetM() != 0)
      {
        this->coef_phi = var_maxwell.coefficient_impedance_absorbing(ref)*nat_mat.GetCoefDamping();
      }
    else
      {
        abort();
      }
  }
  
  
  /********************
   * StaticMaxwell_2D *
   ********************/
  
  
  template<class TypeEquation>
  void StaticMaxwell_2D<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
		      IVect& num_ref, IVect& num_index, IVect& num_component,
		      Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    HarmonicMaxwell_2D<Complexe>::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
						    compute_grad, compute_hess);
  }


  void EllipticProblem<StaticMaxwellEquation_2D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementHcurl(i));
  }
  
  
  void EllipticProblem<StaticMaxwellEquation_2D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementHcurl(i));
  }


  void EllipticProblem<StaticMaxwellEquation_2D_DG>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementH1(i));
  }
  
  
  void EllipticProblem<StaticMaxwellEquation_2D_DG>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementH1(i));
  }
  
  
  void EllipticProblem<StaticMaxwellEquation_2D_DG>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }

  
  void EllipticProblem<StaticMaxwellEquation_2D_DG>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
  
}

#define MONTJOIE_FILE_STATIC_MAXWELL_2D_CXX
#endif

