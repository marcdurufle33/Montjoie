#ifndef MONTJOIE_FILE_STATIC_MAXWELL_3D_CXX

namespace Montjoie
{
  
  /*******************************
   * StaticMaxwellEquation_3D_DG *
   *******************************/


  bool StaticMaxwellEquation_3D_DG::store_dfjm1(true);
  
  bool StaticMaxwellEquation_3D_DG::SymmetricGlobalMatrix()
  {
    return true;
  }
  
  
  bool StaticMaxwellEquation_3D_DG::SymmetricElementaryMatrix()
  {
    return false;
  }
  
  
  //! computation of mass coefficients
  template<class TypeEquation>
  void StaticMaxwellEquation_3D_DG::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var, int i,
		    const ElementReference_Dim<Dimension3>&)
  {
  }
  

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
  void StaticMaxwellEquation_3D_DG
  ::GetTensorMass(const EllipticProblem<TypeEquation>& var,
		  int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		  int ref, MatMass& mass)
  {
    FillZero(mass);
    mass(0, 0) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 0)*nat_mat.GetCoefMass();
    mass(1, 1) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 1, 1)*nat_mat.GetCoefMass();
    mass(2, 2) = var.ref_epsilon(ref).GetCoefficient(var, i, j, 2, 2)*nat_mat.GetCoefMass();
    
    mass(3, 3) = -var.ref_mu(ref).GetCoefficient(var, i, j, 0, 0)*nat_mat.GetCoefMass();
    mass(4, 4) = -var.ref_mu(ref).GetCoefficient(var, i, j, 1, 1)*nat_mat.GetCoefMass();
    mass(5, 5) = -var.ref_mu(ref).GetCoefficient(var, i, j, 2, 2)*nat_mat.GetCoefMass();
    
  }
  

  //! Applies matrix M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U vector to be multiplied by M
    \param[out] V result vector V = M U
  */  
  template<class TypeEquation, class T0, class Vector1>
  void StaticMaxwellEquation_3D_DG
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		    const GlobalGenericMatrix<T0>& nat_mat, int ref,
		    Vector1& U, Vector1& V)
  {
    V = U;
    V(0) *= var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 0)*nat_mat.GetCoefMass();
    V(1) *= var.ref_epsilon(ref).GetCoefficient(var, i, j, 1, 1)*nat_mat.GetCoefMass();
    V(2) *= var.ref_epsilon(ref).GetCoefficient(var, i, j, 2, 2)*nat_mat.GetCoefMass();
    
    V(3) *= -var.ref_mu(ref).GetCoefficient(var, i, j, 0, 0)*nat_mat.GetCoefMass();
    V(4) *= -var.ref_mu(ref).GetCoefficient(var, i, j, 1, 1)*nat_mat.GetCoefMass();
    V(5) *= -var.ref_mu(ref).GetCoefficient(var, i, j, 2, 2)*nat_mat.GetCoefMass();
  }
  

  //! which derivatives to evaluate during matrix-vector product ?  
  template<class TypeEquation, class T0, class Vector1>
  void StaticMaxwellEquation_3D_DG
  ::GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
			const GlobalGenericMatrix<T0>& nat_mat,
			Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    // on a besoin de ne deriver que la premiere inconnue, le champ electrique E
    unknown_to_derive(0) = true;
    unknown_to_derive(1) = true;
    unknown_to_derive(2) = true;
    
    unknown_to_derive(3) = false;
    unknown_to_derive(4) = false;
    unknown_to_derive(5) = false;
    

    fct_test_to_derive(0) = true;
    fct_test_to_derive(1) = true;
    fct_test_to_derive(2) = true;
        
    fct_test_to_derive(3) = false;
    fct_test_to_derive(4) = false;
    fct_test_to_derive(5) = false;
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
  void StaticMaxwellEquation_3D_DG
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc,
		     const GlobalGenericMatrix<T0>& nat_mat, int ref,
		     MatStiff& Ephi_grad, MatStiff& Dgrad_phi)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    T0 c = nat_mat.GetCoefStiffness();
    Dgrad_phi(4,2)(0) = -c; Dgrad_phi(5,1)(0) = c;
    Dgrad_phi(3,2)(1) = c; Dgrad_phi(5,0)(1) = -c;
    Dgrad_phi(3,1)(2) = -c; Dgrad_phi(4,0)(2) = c;    

    Ephi_grad(2,4)(0) = -c; Ephi_grad(1,5)(0) = c;
    Ephi_grad(2,3)(1) = c; Ephi_grad(0,5)(1) = -c;
    Ephi_grad(1,3)(2) = -c; Ephi_grad(0,4)(2) = c;

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
  template<class TypeEquation,
	   class T0, class Vector1, class Vector2>
  void StaticMaxwellEquation_3D_DG
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& Vn, Vector2& Un)
  {
    Un.Fill(0);
    Un(3) = Vn(2)(1) - Vn(1)(2);
    Un(4) = -Vn(2)(0) + Vn(0)(2);
    Un(5) = Vn(1)(0) - Vn(0)(1);
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
  template<class TypeEquation,
	   class T0, class Vector1, class Vector2>
  void StaticMaxwellEquation_3D_DG
  ::ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& Un, Vector2& Vn)  
  {
    Vn.Fill(0);
    Vn(2)(1) = -Un(3); Vn(1)(2) = Un(3);
    Vn(2)(0) = Un(4); Vn(0)(2) = -Un(4);
    Vn(1)(0) = -Un(5); Vn(0)(1) = Un(5);
    Vn *= -nat_mat.GetCoefStiffness();
  }
  
  
  //! computation of matrix N associated to the boundary condition
  /*!
    \param[out] Nabc matrix N
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref_d reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Matrix1, class GenericPb, class T0>
  void StaticMaxwellEquation_3D_DG
  ::GetNabc(Matrix1& Nabc, R3& normale, int ref, int iquad, int k,
	    const GlobalGenericMatrix<T0>& nat_mat, int ref_d, const GenericPb& vars,
	    const ElementReference<Dimension3, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Nabc.Fill(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	Nabc(0,4) = -normale(2); Nabc(0,5) = normale(1);
	Nabc(1,5) = -normale(0); Nabc(1,3) = normale(2);
	Nabc(2,3) = -normale(1); Nabc(2,4) = normale(0);
	Nabc(3,1) = normale(2); Nabc(3,2) = -normale(1);
	Nabc(4,2) = normale(0); Nabc(4,0) = -normale(2);
	Nabc(5,0) = normale(1); Nabc(5,1) = -normale(0);        
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	Nabc(0,4) = normale(2); Nabc(0,5) = -normale(1);
	Nabc(1,5) = normale(0); Nabc(1,3) = -normale(2);
	Nabc(2,3) = normale(1); Nabc(2,4) = -normale(0);
	Nabc(3,1) = -normale(2); Nabc(3,2) = normale(1);
	Nabc(4,2) = -normale(0); Nabc(4,0) = normale(2);
	Nabc(5,0) = -normale(1); Nabc(5,1) = normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        Real_wp c = vars.coefficient_impedance_absorbing(ref_d);
	Nabc(3,3) = -c*(1.0 - normale(0)*normale(0));
        Nabc(3,4) = c*normale(0)*normale(1);
        Nabc(3,5) = c*normale(0)*normale(2);
	
        Nabc(4,4) = -c*(1.0 - normale(1)*normale(1));
        Nabc(4,3) = Nabc(3,4);Nabc(4,5) = c*normale(2)*normale(1);
	
        Nabc(5,5) = -c*(1.0 - normale(2)*normale(2));
        Nabc(5,3) = Nabc(3,5); Nabc(5,4) = Nabc(4,5);
        
        //c = -1.0/c;
        c = 1.0/c;
        Nabc(0,0) = -c*(1.0 - normale(0)*normale(0));
        Nabc(0,2) = c*normale(0)*normale(2); Nabc(0,1) = c*normale(0)*normale(1);
	
        Nabc(1,1) = -c*(1.0 - normale(0)*normale(0));
        Nabc(1,0) = Nabc(0,1); Nabc(1,2) = c*normale(2)*normale(1);
	
        Nabc(2,2) = -c*(1.0 - normale(0)*normale(0));
        Nabc(2,0) = Nabc(0,2); Nabc(2,1) = Nabc(1,2);		
      }
    
    Nabc *= nat_mat.GetCoefStiffness();
  }
  
  
  //! multiplication by matrix involved in penalization terms
  template<class Vector1, class Vector2, class GenericPb, class T0>
  void StaticMaxwellEquation_3D_DG
  ::MltPenalDG(const R3& normale, const Vector1& Vn, Vector2& Un,
	       int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
	       int ref, int ref2, const GenericPb& vars,
	       const ElementReference<Dimension3, 1>& Fb)
  {
    // condition absorbante sur les sauts
    Real_wp E_dot_n = normale(0)*Vn(0) + normale(1)*Vn(1) + normale(2)*Vn(2);
    int ref_d = vars.mesh.Element(i).GetReference();
    Real_wp c = vars.alpha_penalization
      /vars.coefficient_impedance_absorbing(ref_d)*nat_mat.GetCoefStiffness();
    
    Un(0) = -c*(Vn(0) - normale(0)*E_dot_n);
    Un(1) = -c*(Vn(1) - normale(1)*E_dot_n);
    Un(2) = -c*(Vn(2) - normale(2)*E_dot_n);
    
    E_dot_n = normale(0)*Vn(3) + normale(1)*Vn(4) + normale(2)*Vn(5);
    c = -vars.delta_penalization*vars.coefficient_impedance_absorbing(ref_d)
      *nat_mat.GetCoefStiffness();
    
    Un(3) = -c*(Vn(3) - normale(0)*E_dot_n);
    Un(4) = -c*(Vn(4) - normale(1)*E_dot_n);
    Un(5) = -c*(Vn(5) - normale(2)*E_dot_n);
  }
    
  
  //! matrix involved in penalization terms
  template<class Matrix1, class GenericPb, class T0>
  void StaticMaxwellEquation_3D_DG
  ::GetPenalDG(Matrix1& Nabc, R3& normale, int iquad, int k, int nf,
	       const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	       const GenericPb& vars, const ElementReference<Dimension3, 1>& Fb)
  {
    Nabc.Fill(0);
    
    int ref_d = vars.mesh.Element(iquad).GetReference();
    T0 c = -vars.delta_penalization*vars.coefficient_impedance_absorbing(ref_d)
      *nat_mat.GetCoefStiffness();    
    
    Nabc(3,3) = -c*(1.0 - normale(0)*normale(0));
    Nabc(3,4) = c*normale(0)*normale(1);
    Nabc(3,5) = c*normale(0)*normale(2);
    
    Nabc(4,4) = -c*(1.0 - normale(1)*normale(1));
    Nabc(4,3) = Nabc(3,4);Nabc(4,5) = c*normale(2)*normale(1);
    
    Nabc(5,5) = -c*(1.0 - normale(2)*normale(2));
    Nabc(5,3) = Nabc(3,5); Nabc(5,4) = Nabc(4,5);
    
    c = vars.alpha_penalization/vars.coefficient_impedance_absorbing(ref_d)
      *nat_mat.GetCoefStiffness();
    
    Nabc(0,0) = -c*(1.0 - normale(0)*normale(0));
    Nabc(0,2) = c*normale(0)*normale(2); Nabc(0,1) = c*normale(0)*normale(1);
    
    Nabc(1,1) = -c*(1.0 - normale(1)*normale(1));
    Nabc(1,0) = Nabc(0,1); Nabc(1,2) = c*normale(2)*normale(1);
    
    Nabc(2,2) = -c*(1.0 - normale(2)*normale(2));
    Nabc(2,0) = Nabc(0,2); Nabc(2,1) = Nabc(1,2);		
    
  }
  
  
  //! multiplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] num_elem1 element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] refd reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Vector1, class TypeEquation, class T0>
  void StaticMaxwellEquation_3D_DG
  ::MltNabc(R3& normale, int ref, const Vector1& Vn, Vector1& Un,
	    int num_elem1, int k, const GlobalGenericMatrix<T0>& nat_mat,
	    int ref_d, const EllipticProblem<TypeEquation>& vars,
	    const ElementReference<Dimension3, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Un.Fill(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	Un(0) = -Vn(4)*normale(2) + Vn(5)*normale(1);
	Un(1) = -Vn(5)*normale(0) + Vn(3)*normale(2);
	Un(2) = -Vn(3)*normale(1) + Vn(4)*normale(0);
	Un(3) = Vn(1)*normale(2) - Vn(2)*normale(1);
	Un(4) = Vn(2)*normale(0) - Vn(0)*normale(2);
	Un(5) = Vn(0)*normale(1) - Vn(1)*normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	Un(0) = Vn(4)*normale(2) - Vn(5)*normale(1);
	Un(1) = Vn(5)*normale(0) - Vn(3)*normale(2);
	Un(2) = Vn(3)*normale(1) - Vn(4)*normale(0);
	Un(3) = -Vn(1)*normale(2) + Vn(2)*normale(1);
	Un(4) = -Vn(2)*normale(0) + Vn(0)*normale(2);
	Un(5) = -Vn(0)*normale(1) + Vn(1)*normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        Real_wp E_dot_n = normale(0)*Vn(0) + normale(1)*Vn(1) + normale(2)*Vn(2);
        Real_wp c = -1.0/vars.coefficient_impedance_absorbing(ref_d);
        Un(0) = -c*(Vn(0) - normale(0)*E_dot_n);
        Un(1) = -c*(Vn(1) - normale(1)*E_dot_n);
        Un(2) = -c*(Vn(2) - normale(2)*E_dot_n);
        
        E_dot_n = normale(0)*Vn(3) + normale(1)*Vn(4) + normale(2)*Vn(5);
        c = vars.coefficient_impedance_absorbing(ref_d);
        Un(3) = -c*(Vn(3) - normale(0)*E_dot_n);
        Un(4) = -c*(Vn(4) - normale(1)*E_dot_n);
        Un(5) = -c*(Vn(5) - normale(2)*E_dot_n);
      }
    Un *= nat_mat.GetCoefStiffness();
    // DISP(Un); DISP(Vn);
  }    

  /*****************************************
   * IncidentWaveProjector_StaticMaxwell3D *
   *****************************************/


  //! initialisation of source with time t0
  void IncidentWaveProjector_StaticMaxwell3D::Init(const Real_wp& t0)
  {
    wave_pulse = this->incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        cout << "Parameters needed for source" << endl;
        abort();
      }
    
    this->t = t0;
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();    
    k_wave = var_problem.GetWaveVector();
    var_problem.GetPolarization(polar);
    Mlt(1.0/Norm2(k_wave), k_wave);
    TimesProd(k_wave, polar, polarH);
    
    if (param(1) == "AUTO")
      {    
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());        
      }
    else
      offset_arg = to_num<Real_wp>(param(1));
  }
    
  
  //! Evaluates incident field
  void IncidentWaveProjector_StaticMaxwell3D::EvaluateFunction(int i, int j, const R3& x, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;
    
    Real_wp arg = this->t + (offset_arg - DotProd(k_wave, x))/c0;
    Real_wp fs = wave_pulse->Evaluate(arg);
    f(0) = fs*polar(0);
    f(1) = fs*polar(1);
    f(2) = fs*polar(2);
  }

  
  //! curl of incident field
  void IncidentWaveProjector_StaticMaxwell3D::EvaluateFunctionGradient(int i, int j, const R3& x, VectReal_wp& f, VectReal_wp& df)
  {
    if (wave_pulse == NULL)
      return;
    
    Real_wp arg = this->t + (offset_arg - DotProd(k_wave, x))/c0;
    Real_wp fs = wave_pulse->Evaluate(arg);
    Real_wp dfs = -wave_pulse->EvaluateDerivative(arg)/c0;

    f(0) = fs*polar(0);
    f(1) = fs*polar(1);
    f(2) = fs*polar(2);
    
    df(0) = dfs*polarH(0);
    df(1) = dfs*polarH(1);
    df(2) = dfs*polarH(2);
  }
  

  /****************************************
   * DiffractedWaveSource_StaticMaxwell3D *
   ****************************************/
  

  //! initialisation of source with time t0  
  Real_wp DiffractedWaveSource_StaticMaxwell3D
  ::Init(const Real_wp& t0, const Real_wp& dt, int print_level, int n, bool scal)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    Real_wp tlimit(0);
    wave_pulse = incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        abort();
      }
    
    t = t0;
    scalar_eq = scal;
    
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();    
    k_wave = var_problem.GetWaveVector();
    var_problem.GetPolarization(polar);
    Mlt(1.0/Norm2(k_wave), k_wave);
    TimesProd(k_wave, polar, polarH);
    
    if (param(1) == "AUTO")
      {    
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());        

        Real_wp off_max = GetMaximumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
					    var_problem.GetYmin(), var_problem.GetYmax(),
					    var_problem.GetZmin(), var_problem.GetZmax());        
	
	tlimit = off_max - offset_arg;
      }
    else
      {
	offset_arg = to_num<Real_wp>(param(1));
	tlimit = to_num<Real_wp>(param(2));
      }
    
    if ((t0 <= dt+20*epsilon_machine) && (print_level >= 2) && (rank_proc == 0))
      cout << "Offset in expression of plane wave " << offset_arg << endl;
    
    num_deriv = n;
    
    if (n > 1)
      {
        cout << "Higher-order derivatives not implemented" << endl;
        abort();
      }
    
    return tlimit;
  }

   
  //! returns true for Dirichlet condition
  bool DiffractedWaveSource_StaticMaxwell3D::PresenceDirichlet() const
  {
    return dirichlet_cond;
  }
   

  //! For Dirichlet condition
  void DiffractedWaveSource_StaticMaxwell3D::EvaluateFunction(int i, int j, const R3& x, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;

    Real_wp arg = t + (offset_arg - DotProd(k_wave, x))/c0;
    Real_wp fs(0);
    if (num_deriv == 0)
      fs = wave_pulse->Evaluate(arg);
    else
      fs = wave_pulse->EvaluateDerivative(arg);
    
    f(0) = -fs*polar(0);    
    f(1) = -fs*polar(1);    
    f(2) = -fs*polar(2);
  }
    
  
  //! initialisation of element i for volume integrals
  void DiffractedWaveSource_StaticMaxwell3D::InitElement(int num_elem, const VectR3& s)
  {
  }
   
  
  //! returns true if there is a volume integral
  bool DiffractedWaveSource_StaticMaxwell3D::IsNonNullVolumetricSource(const VectR3& s)
  {
    return false;
  }
  
  
  //! fills f for volume integral \int f \varphi
  void DiffractedWaveSource_StaticMaxwell3D::EvaluateVolumetricSource(int i, int j, const R3& x, VectReal_wp& f)
  {
  }
    
  
  //! returns true if there a volume integral \int f \nabla \varphi
  bool DiffractedWaveSource_StaticMaxwell3D::IsNonNullGradientSource(const VectR3& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \nabla \varphi
  void DiffractedWaveSource_StaticMaxwell3D::EvaluateGradientSource(int i, int j, const R3& x, VectReal_wp& f)
  {
  }

  
  //! initialisation of surface i
  void DiffractedWaveSource_StaticMaxwell3D::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Real_wp, Dimension3>::InitSurface(i, num_face, num_elem, num_loc);
  }
    
  
  //! fills f for surface integral \int f \varphi
  void DiffractedWaveSource_StaticMaxwell3D::EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
                                                                    const SetMatrices<Dimension3>& MatricesElem, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;

    if ((!scalar_eq) || var_problem.InsidePML(this->num_elem_))
      {
        f.Fill(0);
        return;
      }
    
    // coef = -n \times polarH / mu
    R3 coef;
    TimesProd(polarH, MatricesElem.GetNormaleQuadratureBoundary(k), coef);
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    coef *= var_maxwell.ref_invMu(ref).GetConstant()(0, 0);
    
    // on evalue t - (k.x + offset_arg)/c
    Real_wp arg = t + (offset_arg - DotProd(k_wave, PointsElem.GetPointQuadratureBoundary(k)))/c0;
    Real_wp fs(0);

    if (num_deriv == 0)
      {
        // primitive en temps
        fs = wave_pulse->Evaluate(arg)/c0;
      }
    else
      {
        fs = wave_pulse->EvaluateDerivative(arg)/c0; 
      }
    
    f(0) = coef(0)*fs;
    f(1) = coef(1)*fs;
    f(2) = coef(2)*fs;
  }
  
  
  bool DiffractedWaveSource_StaticMaxwell3D::IsNonNullSurfacicSource(int ref)
  {
    if (var_problem.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEUMANN)
      return true;
    
    return false;
  }
  
   
  /***********************************
   * TotalWaveSource_StaticMaxwell3D *
   ***********************************/

  
  Real_wp TotalWaveSource_StaticMaxwell3D::Init(const Real_wp& t0, const Real_wp& dt, int print_level,
                                                int n, bool scal)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    Real_wp tlimit(0);
    wave_pulse = incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        abort();
      }
    
    t = t0;
    scalar_eq = scal;
    
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();    
    k_wave = var_problem.GetWaveVector();
    var_problem.GetPolarization(polar);
    Mlt(1.0/Norm2(k_wave), k_wave);
    TimesProd(k_wave, polar, polarH);
    
    if (param(1) == "AUTO")
      {    
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());        

        Real_wp off_max = GetMaximumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
					    var_problem.GetYmin(), var_problem.GetYmax(),
					    var_problem.GetZmin(), var_problem.GetZmax());        
	
	tlimit = off_max - offset_arg;
      }
    else
      {
	offset_arg = to_num<Real_wp>(param(1));
	tlimit = to_num<Real_wp>(param(2));
      }
    
    if ((t0 <= dt+20*epsilon_machine) && (print_level >= 2) && (rank_proc == 0))
      cout << "Offset in expression of plane wave " << offset_arg << endl;
    
    num_deriv = n;
    
    if (n > 1)
      {
        cout << "Higher-order derivatives not implemented" << endl;
        abort();
      }
    
    return tlimit;
  }
    
   
  //! For Dirichlet condition
  void TotalWaveSource_StaticMaxwell3D::EvaluateFunction(int i, int j, const R3& x, VectReal_wp& f)
  {
  }
    
  
  //! initialisation of element i for volume integrals
  void TotalWaveSource_StaticMaxwell3D::InitElement(int num_elem, const VectR3& s)
  {
  }
   

  //! returns true if there is a volume integral
  bool TotalWaveSource_StaticMaxwell3D::IsNonNullVolumetricSource(const VectR3& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \varphi  
  void TotalWaveSource_StaticMaxwell3D::EvaluateVolumetricSource(int i, int j, const R3& x, VectReal_wp& f)
  {
  }
    

  //! returns true if there a volume integral \int f \nabla \varphi  
  bool TotalWaveSource_StaticMaxwell3D::IsNonNullGradientSource(const VectR3& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \nabla \varphi  
  void TotalWaveSource_StaticMaxwell3D::EvaluateGradientSource(int i, int j, const R3& x, VectReal_wp& f)
  {
  }
    
   
  //! initialisation of surface i
  void TotalWaveSource_StaticMaxwell3D::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Real_wp, Dimension3>::
      InitSurface(i, num_face, num_elem, num_loc);
  }

  
  //! fills f for surface integral \int f \varphi  
  void TotalWaveSource_StaticMaxwell3D
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension3>& PointsElem,
                           const SetMatrices<Dimension3>& MatricesElem, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;
    
    if (!scalar_eq)
      {
        f.Zero();
        return;
      }
    
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();

    R3 coef, En, normale, diffE;
    normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    TimesProd(normale, polar, En);
    diffE = En - polarH;
    TimesProd(diffE, normale, coef);
    
    coef *= var_maxwell.ref_invMu(ref).GetConstant()(0, 0) / c0;
    
    Real_wp arg = t + (offset_arg - DotProd(k_wave, PointsElem.GetPointQuadratureBoundary(k)))/c0;
    Real_wp fs(0);
    if (num_deriv == 0)
      fs = wave_pulse->Evaluate(arg);
    else
      fs = wave_pulse->EvaluateDerivative(arg);

    f(0) = coef(0)*fs;
    f(1) = coef(1)*fs;
    f(2) = coef(2)*fs;
  }
    

  //! returns true if there is a surface integral \int f \varphi
  bool TotalWaveSource_StaticMaxwell3D::IsNonNullSurfacicSource(int ref)
  {
    if (var_problem.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    return false;
  }
  
 
  /********************
   * StaticMaxwell_3D *
   ********************/

  
  template<class TypeEquation>
  void StaticMaxwell_3D<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension3, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension3, Real_wp>* >& rho_real, IVect& num_ref,
		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    HarmonicMaxwell_3D<typename TypeEquation::Complexe>::
      GetVaryingIndices(rho_real, num_ref, num_index,
			num_component, compute_grad, compute_hess);
  }


  /*******************
   * EllipticProblem *
   *******************/


  IncidentWaveProjector<Real_wp, Dimension3>*
  EllipticProblem<StaticMaxwellEquation_3D>
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Real_wp, Dimension3>& u_inc) const
  {
    IncidentWaveProjector_StaticMaxwell3D* proj;
    proj = new IncidentWaveProjector_StaticMaxwell3D(*this, u_inc);
    proj->Init(0.0);
    return proj;
  }
  

  void EllipticProblem<StaticMaxwellEquation_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    HarmonicMaxwell_3D<Real_wp>::
      ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, this->GetReferenceElementHcurl(i));
  }
    
  void EllipticProblem<StaticMaxwellEquation_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    HarmonicMaxwell_3D<Real_wp>::
      ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, this->GetReferenceElementHcurl(i));
  }

  void EllipticProblem<StaticMaxwellEquation_3D>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    HarmonicMaxwell_3D<Real_wp>::
      AddElementaryFluxesSipg(mat_sp, nat_mat, offset_row, offset_col);
  }

  
  void EllipticProblem<StaticMaxwellEquation_3D>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    HarmonicMaxwell_3D<Real_wp>::
      AddElementaryFluxesSipg(mat_sp, nat_mat, offset_row, offset_col);
  }

  void EllipticProblem<StaticMaxwellEquation_3D>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& curl_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const
  {
    HarmonicMaxwell_3D_Base::ComputeEnHnQuadratureHcurl(u_quadrature, curl_quadrature,
						   num_elem, pts, normale, compute_H, En_quad, Hn_quad);
  }


  void EllipticProblem<StaticMaxwellEquation_3D>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& curl_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const
  {
  }
  
  
  void EllipticProblem<StaticMaxwellEquation_3D>
  ::ComputeEnHnNodal(Vector<VectReal_wp >& u_nodal, Vector<VectReal_wp >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<VectReal_wp >& En_nodal, Vector<VectReal_wp >& Hn_nodal) const
  {
    HarmonicMaxwell_3D_Base::ComputeEnHnNodal(u_nodal, grad_nodal, num_elem, pts, normale,
					      En_nodal, Hn_nodal);
  }


  void EllipticProblem<StaticMaxwellEquation_3D>
  ::ComputeEnHnNodal(Vector<VectComplex_wp >& u_nodal, Vector<VectComplex_wp >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<VectComplex_wp >& En_nodal, Vector<VectComplex_wp >& Hn_nodal) const
  {
  }


  void EllipticProblem<StaticMaxwellEquation_3D_DG>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }
  
  
  void EllipticProblem<StaticMaxwellEquation_3D_DG>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }
    
  
  void EllipticProblem<StaticMaxwellEquation_3D_DG>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
                          const GlobalGenericMatrix<Real_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
  
  
  void EllipticProblem<StaticMaxwellEquation_3D_DG>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }


  EllipticProblem<StaticMaxwellEquationHdg_3D>::EllipticProblem()
  {
    this->dg_formulation = ElementReference_Base::HDG;
  }
  
  void EllipticProblem<StaticMaxwellEquationHdg_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>& cond_solver,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHdg(i, num_dof, mat_elem,
				     cond_solver, nat_mat, this->GetReferenceElementHcurl(i));
  }
  
  
  void EllipticProblem<StaticMaxwellEquationHdg_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>& cond_solver,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHdg(i, num_dof, mat_elem,
				     cond_solver, nat_mat, this->GetReferenceElementHcurl(i));
  }
  
  
  /**********************
   * FemMatrixFreeClass *
   **********************/
  

  void FemMatrixFreeClass<Real_wp, StaticMaxwellEquation_3D>
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddHcurl3D(Real_wp(1), nat_mat, trans, lvl, *this, 
		  X, Real_wp(1), Y, false);
  }
  

  void FemMatrixFreeClass<Real_wp, StaticMaxwellEquation_3D>
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddHcurl3D(Complex_wp(1), nat_mat, trans, lvl, *this, 
		  X, Complex_wp(1), Y, false);
  }


  void FemMatrixFreeClass<Complex_wp, StaticMaxwellEquation_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Types not compatible" << endl;
    abort();
  }
  

  void FemMatrixFreeClass<Complex_wp, StaticMaxwellEquation_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddHcurl3D(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		  X, Complex_wp(1, 0), Y, false);
  }


  void FemMatrixFreeClass<Real_wp, StaticMaxwellEquationHdg_3D>
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddHcurlHdg3D(Real_wp(1), nat_mat, trans, lvl, *this, 
		     X, Real_wp(1), Y, false);
  }
  

  void FemMatrixFreeClass<Real_wp, StaticMaxwellEquationHdg_3D>
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    cout << "Types not compatible" << endl;
    abort();
  }


  void FemMatrixFreeClass<Complex_wp, StaticMaxwellEquationHdg_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Types not compatible" << endl;
    abort();
  }
  

  void FemMatrixFreeClass<Complex_wp, StaticMaxwellEquationHdg_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddHcurlHdg3D(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		     X, Complex_wp(1, 0), Y, false);
  }

}

#define MONTJOIE_FILE_STATIC_MAXWELL_3D_CXX
#endif

