#ifndef MONTJOIE_FILE_SYMMETRIC_HARMONIC_GALBRUN_CXX

namespace Montjoie
{
  
  /***********************
   * GalbrunEquationSipg *
   ***********************/
  

  template<class T, class Dimension>
  void GalbrunEquationSipg_Base<T, Dimension>
  ::SetIndexToCompute(VarGalbrunIndex_Base<Dimension>& var)
  {
    var.compute_grad_flow = false;
    var.compute_grad_rho = false;
    var.compute_grad_sigma = false;
    var.compute_grad_c0 = false;
    var.compute_hessian_flow = false;
    var.compute_hessian_p0 = false;
  }


  //! for compatbility purpose
  template<class T, class Dimension> template<class TypeEquation>
  void GalbrunEquationSipg_Base<T, Dimension>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                    int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }
  
  
  //! which derivatives to evaluate during matrix-vector product ?
  template<class T, class Dimension> template<class GenericPb, class T0, class Vector1>
  void GalbrunEquationSipg_Base<T, Dimension>
  ::GetNeededDerivative(const GenericPb& vars, const GlobalGenericMatrix<T0>& nat_mat,
			Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    fct_test_to_derive.Fill(true);
    unknown_to_derive.Fill(true);
  }
  
  
  //! Application of stiffness operator C
  /*!
    \param[in] var considered problem
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dU vector containing \nabla u
    \param[in] dV result C \nabla u
    Tensor C is involved in the term \int C \nabla u \nabla v
   */
  template<class T, class Dimension> 
  template<class TypeEquation, class T0, class Vector1>
  void GalbrunEquationSipg_Base<T, Dimension>::
  ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
		       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		       int ref, Vector1& dU, Vector1& dV)
  {
    // 2-D case
    TinyVector<Real_wp, Dimension::dim_N> v0 = var.eval_flow(i)(j);
    
    FillZero(dV);
    
    T0 Mu = DotProd(dU(0), v0)*nat_mat.GetCoefStiffness();
    dV(0)(0) -= Mu*v0(0);
    dV(0)(1) -= Mu*v0(1);

    Mu = DotProd(dU(1), v0)*nat_mat.GetCoefStiffness();
    dV(1)(0) -= Mu*v0(0);
    dV(1)(1) -= Mu*v0(1);
    
    // part rho c^2 div(u) div(phi)
    Mu = (dU(0)(0) + dU(1)(1))*var.eval_rhoC2(i)(j)*nat_mat.GetCoefStiffness();
    dV(0)(0) += Mu;
    dV(1)(1) += Mu;
  }
  
  
  //! Computation of stiffness tensor C
  /*!
    \param[in] vars considered problem
    \param[in] num_elem element number
    \param[in] jloc quadrature point number
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Cgrad_grad tensor C
    Tensor C is involved in the term \int C \nabla u \nabla v
   */
  template<class T, class Dimension> 
  template<class TypeEquation, class T0, class MatStiff>
  void GalbrunEquationSipg_Base<T, Dimension>::
  GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
		    int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
		    int ref, MatStiff& Cgrad_grad)
  {
    FillZero(Cgrad_grad);

    TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(num_elem)(jloc);
    TinyMatrix<Real_wp, Symmetric, Dimension::dim_N, Dimension::dim_N> Mm;    
    GetNormalProjector(v0, Mm);
    Mlt(vars.eval_rho0(num_elem)(jloc), Mm);
    
    // part -rho m \cdot grad u  m \cdot grad phi
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        Cgrad_grad(i, i) = Mm;
        Mlt(-nat_mat.GetCoefStiffness(), Cgrad_grad(i, i));
      }
    
    // part rho c^2 div(u) div(phi)
    T0 coef = nat_mat.GetCoefStiffness();
    coef *= vars.eval_rhoC2(num_elem)(jloc);
    for (int i = 0; i < Dimension::dim_N; i++)
      for (int j = 0; j < Dimension::dim_N; j++)
        Cgrad_grad(i, j)(i, j) += coef;
    
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
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class MatStiff>
  void GalbrunEquationSipg_Base<T, Dimension>::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    // part m \cdot grad u phi
    TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(num_elem)(jloc);
    Mlt(vars.eval_rho0(num_elem)(jloc), v0);    
        
    T m_iomega; vars.GetMiomega(m_iomega);
    T0 coef;
    coef = 2.0*(nat_mat.GetCoefDamping()*m_iomega +
                vars.eval_sigma(num_elem)(jloc)*nat_mat.GetCoefStiffness());
    
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        for (int p = 0; p < Dimension::dim_N; p++)
          {
            Ephi_grad(i, i)(p) += coef*v0(p);
            Ephi_grad(i, p)(p) += nat_mat.GetCoefStiffness()*vars.grad_p0(num_elem)(jloc)(i);
            Ephi_grad(i, p)(i) -= nat_mat.GetCoefStiffness()*vars.grad_p0(num_elem)(jloc)(p);
          }
      }        
  }
  

  //! Applying the tensor D to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] iquad element number
    \param[in] k quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Vn gradient of the unknown vector V
    \param[out] Un result to D grad(v)
  */  
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void GalbrunEquationSipg_Base<T, Dimension>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
		       int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat,
		       int ref, Vector1& Vn, Vector2& Un)
  {
    Un.Fill(0);
    TinyVector<Real_wp, Dimension::dim_N> v0 = var.eval_flow(iquad)(k);
    v0 *= var.eval_rho0(iquad)(k);

    T m_iomega; var.GetMiomega(m_iomega);
    T0 coef;
    coef = 2.0*(nat_mat.GetCoefDamping()*m_iomega +
                var.eval_sigma(iquad)(k)*nat_mat.GetCoefStiffness());    
    
    for (int i = 0; i < Dimension::dim_N; i++)
      {
	// flow part
	for (int p = 0; p < Dimension::dim_N; p++)
	  Un(p) += coef*v0(i)*Vn(p)(i);
      }
  }
  

  //! Applies matrix M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un vector to be multiplied by M
    \param[out] Vn result vector V = M U
  */  
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class Vector1>
  void GalbrunEquationSipg_Base<T, Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& Un, Vector1& Vn)
  {
    Vn.Fill(0);
  }
  

  //! returns the matrix M, in the integral \f$ \int_K M u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Cj matrix M
  */  
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class MatMass>
  void GalbrunEquationSipg_Base<T, Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& Cj)
  {
    Cj.Fill(0);    
  }
  

  //! computation of matrix N associated to the boundary condition
  /*!
    \param[out] Nabc matrix N
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] iquad element number
    \param[in] npoint local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */  
  template<class T, class Dimension>
  template<class Matrix1, class GenericPb, class T0>
  void GalbrunEquationSipg_Base<T, Dimension>::
  GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
	  int ref, int iquad, int npoint,
	  const GlobalGenericMatrix<T0>& nat_mat, int ref2, 
	  const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
  }
  
  
  //! mutliplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] num_elem1 element number
    \param[in] npoint local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class T, class Dimension>
  template<class Vector1, class TypeEquation, class T0>
  void GalbrunEquationSipg_Base<T, Dimension>::
  MltNabc(typename Dimension::R_N& normale, int ref,
	  const Vector1& Vn, Vector1& Un, int num_elem1,
	  int npoint, const GlobalGenericMatrix<T0>& nat_mat, int ref2, 
	  const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& )
  {
    Un.Fill(0);
  }
  
  
  //! Computation of "penalization" matrices C
  /*!
    \param[out] Nabc penalization matrix C
    \param[in] normale outward normale
    \param[in] iquad element number
    \param[in] npoint local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class T, class Dimension>
  template<class Matrix1, class GenericPb, class T0>
  void GalbrunEquationSipg_Base<T, Dimension>::
  GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int npoint,
	     int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	     const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    //Nabc.SetDiagonal(vars.alpha_penalization*nat_mat.GetCoefStiffness());
    GetNormalProjector(normale, Nabc);
    Mlt(vars.alpha_penalization*vars.eval_rhoC2(iquad)(npoint)*nat_mat.GetCoefStiffness(), Nabc);
  }
  

  //! Multiplication by penalization matrices
  /*!
    \param[in] normale outward normale
    \param[in] Vn vector to multiply
    \param[out] Un result vector Un = C*Vn
    \param[in] iquad element number
    \param[in] npoint local quadrature point number    
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */  
  template<class T, class Dimension>
  template<class Vector1, class Vector2, class GenericPb, class T0>
  void GalbrunEquationSipg_Base<T, Dimension>::
  MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
	     int iquad, int npoint, int nf, const GlobalGenericMatrix<T0>& nat_mat,
	     int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    //Un = Vn;
    //Un *= vars.alpha_penalization*nat_mat.GetCoefStiffness();
    typename Vector1::value_type v_dot_n = DotProd(Vn, normale);
    Un = normale;
    Un *= vars.alpha_penalization*vars.eval_rhoC2(iquad)(npoint)*nat_mat.GetCoefStiffness()*v_dot_n;
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
  template<class Dimension>
  template<class TypeEquation, class T0, class MatMass>
  void GalbrunStationaryEquationSipg<Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    GalbrunEquationSipg_Base<Real_wp, Dimension>::GetTensorMass(vars, i, j, nat_mat, ref, mass);
    
    T0 coef = nat_mat.GetCoefMass();
    Real_wp sigma = vars.eval_sigma(i)(j);
    coef += sigma*(2.0*nat_mat.GetCoefDamping() + sigma*nat_mat.GetCoefStiffness());
    coef *= vars.eval_rho0(i)(j);
    for (int p = 0; p < mass.GetM(); p++)
      mass(p, p) += coef;
    
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
  template<class Dimension>
  template<class TypeEquation, class T0, class Vector1>
  void GalbrunStationaryEquationSipg<Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    GalbrunEquationSipg_Base<Real_wp, Dimension>::ApplyTensorMass(var, i, j, nat_mat, ref, U, V);
    
    T0 coef = nat_mat.GetCoefMass();
    Real_wp sigma = var.eval_sigma(i)(j);
    coef += sigma*(2.0*nat_mat.GetCoefDamping() + sigma*nat_mat.GetCoefStiffness());
    coef *= var.eval_rho0(i)(j);
    Add(coef, U, V);
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
  template<class Dimension>
  template<class TypeEquation, class T0, class MatMass>
  void HarmonicGalbrunEquationSipg<Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& var,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    GalbrunEquationSipg_Base<Complex_wp, Dimension>::GetTensorMass(var, i, j, nat_mat, ref, mass);
    
    Real_wp sigma = var.eval_sigma(i)(j);
    Complex_wp coef = -var.GetSquareOmega()*nat_mat.GetCoefMass() 
      + sigma*(-2.0*Iwp*var.GetOmega()*nat_mat.GetCoefDamping() + sigma*nat_mat.GetCoefStiffness());
    
    coef *= var.eval_rho0(i)(j);
    for (int p = 0; p < mass.GetM(); p++)
      mass(p, p) += coef;
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
  template<class Dimension>
  template<class TypeEquation, class T0, class Vector1>
  void HarmonicGalbrunEquationSipg<Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    GalbrunEquationSipg_Base<Complex_wp, Dimension>::ApplyTensorMass(var, i, j, nat_mat,
                                                                     ref, U, V);
    
    Real_wp sigma = var.eval_sigma(i)(j);
    Complex_wp coef = -var.GetSquareOmega()*nat_mat.GetCoefMass() - Iwp*var.GetOmega()* 
      + sigma*(-2.0*Iwp*var.GetOmega()*nat_mat.GetCoefDamping() + sigma*nat_mat.GetCoefStiffness());
    
    coef *= var.eval_rho0(i)(j);
    Add(coef, U, V);
  }
  

  /****************
   * ImpedanceABC *
   ****************/
  

  //! evaluation of impedance coefficient
  template<class T, class Dimension>
  void ImpedanceABC<T, HarmonicGalbrunEquationSipg<Dimension> >
  ::EvaluateImpedancePhi_H1(int i, int iquad, int num_edge, int num_loc, int k,
                            const GlobalGenericMatrix<T>& nat_mat, int ref,
                            const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    if (k == 0)
      this->stored_coef_phi(0).Reallocate(Mat.GetNbPointsQuadratureBoundary());
    
    Real_wp mu = var_galbrun.ref_rho0(ref).GetConstant()*var_galbrun.ref_c0(ref).GetConstant();
    T m_iomega; var_problem.GetMiomega(m_iomega);
    this->stored_coef_phi(0)(k) = m_iomega*nat_mat.GetCoefDamping()*mu;
  }
  
  
  //! application of impedance of first-order absorbing boundary condition 
  template<class T, class Dimension>
  inline void ImpedanceABC<T, HarmonicGalbrunEquationSipg<Dimension> >::
  ApplyImpedancePhi_H1(int n, int j, int offset, const TinyVector<Real_wp, 1>& phi,
                       const typename Dimension::R_N& grad_phi, Vector<T>& f_phi)
  {
    T coef = this->vec_normale(j)(n)*phi(0)*this->stored_coef_phi(0)(j);
    TinyVector<T, Dimension::dim_N> vec_u;
    vec_u = this->vec_normale(j);
    vec_u *= coef;
    CopyVector(vec_u, offset, f_phi);
  }

  
  /************************
   * GalbrunEquationH1_2D *
   ************************/
  
  
  template<class T>
  void GalbrunEquationH1_2D<T>
  ::SetIndexToCompute(VarGalbrunIndex_Base<Dimension2>& var)
  {
    var.compute_grad_flow = true;
    var.compute_grad_rho = false;
    var.compute_grad_sigma = false;
    var.compute_grad_c0 = false;
    var.compute_hessian_flow = true;
    var.compute_hessian_p0 = true;
  }

  
  //! for compatbility purpose
  template<class T> template<class TypeEquation>
  void GalbrunEquationH1_2D<T>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                    int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }
  
  
  //! which derivatives to evaluate during matrix-vector product ?
  template<class T> template<class GenericPb, class T0, class Vector1>
  void GalbrunEquationH1_2D<T>::
  GetNeededDerivative(const GenericPb& vars, const GlobalGenericMatrix<T0>& nat_mat, 
                      Vector1& unknown_to_derive, Vector1& fct_test_to_derive) 
  {
    fct_test_to_derive.Fill(true);
    unknown_to_derive.Fill(true);
  }
  
  
  //! Application of stiffness operator C
  /*!
    \param[in] var considered problem
    \param[in] i element number
    \param[in] j quadrature point number
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dU vector containing \nabla u
    \param[in] dV result C \nabla u
    Tensor C is involved in the term \int C \nabla u \nabla v
   */
  template<class T> 
  template<class TypeEquation, class T0, class Vector1>
  void GalbrunEquationH1_2D<T>::
  ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var, int i, int j,
		       const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& dU, Vector1& dV)
  {
    abort();
  }
  
  
  //! Computation of stiffness tensor C
  /*!
    \param[in] vars considered problem
    \param[in] num_elem element number
    \param[in] jloc quadrature point number
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Cgrad_grad tensor C
    Tensor C is involved in the term \int C \nabla u \nabla v
   */
  template<class T> 
  template<class TypeEquation, class T0, class MatStiff>
  void GalbrunEquationH1_2D<T>::
  GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                    int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                    int ref, MatStiff& Cgrad_grad)
  {
    FillZero(Cgrad_grad);

    TinyVector<Real_wp, 2> v0 = vars.eval_flow(num_elem)(jloc);
    TinyMatrix<Real_wp, Symmetric, 2, 2> Mm;    
    GetNormalProjector(v0, Mm);
    Mlt(vars.eval_rho0(num_elem)(jloc), Mm);
    
    // part -rho m \cdot grad u  m \cdot grad phi
    Cgrad_grad(0, 0) = Mm;
    Mlt(-nat_mat.GetCoefStiffness(), Cgrad_grad(0, 0));

    Cgrad_grad(1, 1) = Mm;
    Mlt(-nat_mat.GetCoefStiffness(), Cgrad_grad(1, 1));
    
    if (nb_unknowns != 4)
      Cgrad_grad(2, 2) = Cgrad_grad(0, 0);
        
    // part rho c^2 div(u) div(phi)
    T0 coef = nat_mat.GetCoefStiffness();
    coef *= vars.eval_rhoC2(num_elem)(jloc);
    
    Cgrad_grad(0, 0)(0, 0) += coef;
    Cgrad_grad(0, 0)(1, 1) += coef;
    
    Cgrad_grad(1, 1)(0, 0) += coef;
    Cgrad_grad(1, 1)(1, 1) += coef;
    
    // part rho c^2 rot(u) rot(phi)
    Cgrad_grad(0, 1)(0, 1) = coef;
    Cgrad_grad(0, 1)(1, 0) = -coef;
    
    Cgrad_grad(1, 0)(1, 0) = coef;
    Cgrad_grad(1, 0)(0, 1) = -coef;    

    // part due to B and C terms (non-uniform flow)
    int N = nb_unknowns-1;
    TinyMatrix<Real_wp, General, 2, 2> dv0 = vars.grad_flow(num_elem)(jloc);
    
    // part 2 rho M \cdot \nabla \tilde(psi) \sum_j \nabla M_j \times du/dx_j )
    T0 coef_r
      = 2*nat_mat.GetCoefStiffness()*vars.eval_rho0(num_elem)(jloc);
    
    Cgrad_grad(N, 0)(0, 0) = coef_r*dv0(0, 1)*v0(0);
    Cgrad_grad(N, 0)(1, 0) = coef_r*dv0(0, 1)*v0(1);
    Cgrad_grad(N, 0)(0, 1) = coef_r*dv0(1, 1)*v0(0);
    Cgrad_grad(N, 0)(1, 1) = coef_r*dv0(1, 1)*v0(1);
    Cgrad_grad(N, 1)(0, 0) = -coef_r*dv0(0, 0)*v0(0);
    Cgrad_grad(N, 1)(1, 0) = -coef_r*dv0(0, 0)*v0(1);
    Cgrad_grad(N, 1)(0, 1) = -coef_r*dv0(1, 0)*v0(0);
    Cgrad_grad(N, 1)(1, 1) = -coef_r*dv0(1, 0)*v0(1);	
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
  template<class T>
  template<class TypeEquation, class T0, class MatStiff>
  void GalbrunEquationH1_2D<T>::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    // part m \cdot grad u phi
    TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(num_elem)(jloc);
    Real_wp rho = vars.eval_rho0(num_elem)(jloc);
        
    T m_iomega; vars.GetMiomega(m_iomega);
    T0 coef, coef_m;
    coef = 2.0*(nat_mat.GetCoefDamping()*m_iomega +
                vars.eval_sigma(num_elem)(jloc)*nat_mat.GetCoefStiffness());
    
    coef_m = coef;
    // part 2 rho (-i omega + sigma) M \cdot \nabla u 
    Ephi_grad(0, 0)(0) = coef*rho*v0(0);
    Ephi_grad(0, 0)(1) = coef*rho*v0(1);

    Ephi_grad(1, 1)(0) = coef*rho*v0(0);
    Ephi_grad(1, 1)(1) = coef*rho*v0(1);
    
    if (nb_unknowns == 4)
      {
	// term rho_0 Dpsi/dt and rho_0 Dv/dt
	Ephi_grad(2, 2) = rho*v0;
	Ephi_grad(2, 2) *= 0.5*nat_mat.GetCoefStiffness();
	Ephi_grad(3, 3) = Ephi_grad(2, 2);
	
	Dgrad_phi(2, 2) = -Ephi_grad(2, 2);
	Dgrad_phi(3, 3) = Dgrad_phi(2, 2);
      }
    else
      Ephi_grad(2, 2) = Ephi_grad(0, 0);

    // part div u grad_p0 - (grad u)^T grad_p0
    coef = nat_mat.GetCoefStiffness()*vars.grad_p0(num_elem)(jloc)(0);
    Ephi_grad(0, 1)(1) += coef;
    Ephi_grad(1, 0)(1) -= coef;
    
    coef = nat_mat.GetCoefStiffness()*vars.grad_p0(num_elem)(jloc)(1);
    Ephi_grad(1, 0)(0) += coef;
    Ephi_grad(0, 1)(0) -= coef;
    
    // term -rho c^2 psi rot(u_test)
    coef = nat_mat.GetCoefStiffness();
    coef *= vars.eval_rhoC2(num_elem)(jloc);
    
    Dgrad_phi(0, 2)(1) = coef;    
    Dgrad_phi(1, 2)(0) = -coef;
    
    // part due to B and C terms (non-uniform flow)
    int N = nb_unknowns-1;
    TinyMatrix<Real_wp, General, 2, 2> dv0 = vars.grad_flow(num_elem)(jloc);
    TinyVector<TinyMatrix<Real_wp, Symmetric, 2, 2>, 2>
      hess_v0 = vars.hessian_flow(num_elem)(jloc);
    
    // part 2 rho (-i omega + sigma) \sum_j \nabla M_j \times du/dx_j
    Ephi_grad(N, 0)(0) -= rho*coef_m*dv0(0, 1);
    Ephi_grad(N, 0)(1) -= rho*coef_m*dv0(1, 1);
    Ephi_grad(N, 1)(0) += rho*coef_m*dv0(0, 0);
    Ephi_grad(N, 1)(1) += rho*coef_m*dv0(1, 0);
    
    // part rho \sum_{j, k} dM_k/dx_j \nabla M_j \times du/dx_k
    coef = nat_mat.GetCoefStiffness();
    Ephi_grad(N, 0)(0) -= rho*coef*dv0(0, 1)*(dv0(0, 0) + dv0(1, 1));
    Ephi_grad(N, 0)(1) -= rho*coef*(dv0(0, 1)*dv0(1, 0) + dv0(1, 1)*dv0(1, 1));
    Ephi_grad(N, 1)(0) += rho*coef*(dv0(0, 1)*dv0(1, 0) + dv0(0, 0)*dv0(0, 0));
    Ephi_grad(N, 1)(1) += rho*coef*dv0(1, 0)*(dv0(0, 0) + dv0(1, 1));
    
    // part -rho \sum_{j,k} M_j \nabla( dM_k/dxj) \times du/dx_k
    Ephi_grad(N, 0)(0) += rho*coef*(v0(0)*hess_v0(0)(0, 1) + v0(1)*hess_v0(0)(1, 1));
    Ephi_grad(N, 0)(1) += rho*coef*(v0(0)*hess_v0(1)(0, 1) + v0(1)*hess_v0(1)(1, 1));
    Ephi_grad(N, 1)(0) -= rho*coef*(v0(0)*hess_v0(0)(0, 0) + v0(1)*hess_v0(0)(0, 1));
    Ephi_grad(N, 1)(1) -= rho*coef*(v0(0)*hess_v0(1)(0, 0) + v0(1)*hess_v0(1)(0, 1));
    
    // part -\sum_j \nabla( dp0_dxj) \times \nabla u_j
    TinyMatrix<Real_wp, Symmetric, 2, 2> hess_p0 = vars.hessian_p0(num_elem)(jloc);
    Ephi_grad(N, 0)(0) += coef*hess_p0(0, 1);
    Ephi_grad(N, 0)(1) -= coef*hess_p0(0, 0);
    Ephi_grad(N, 1)(0) += coef*hess_p0(1, 1);
    Ephi_grad(N, 1)(1) -= coef*hess_p0(0, 1);
    
    // part 1/(rho c^2) \sum_j dp0/dx_j grad(p0) \times grad(u_j)
    coef /= vars.eval_rhoC2(num_elem)(jloc);
    R2 dp0 = vars.grad_p0(num_elem)(jloc);
    Ephi_grad(N, 0)(0) -= coef*dp0(0)*dp0(1);
    Ephi_grad(N, 0)(1) += coef*dp0(0)*dp0(0);
    Ephi_grad(N, 1)(0) -= coef*dp0(1)*dp0(1);
    Ephi_grad(N, 1)(1) += coef*dp0(0)*dp0(1);    
  }
  

  //! Applying the tensor D to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] iquad element number
    \param[in] k quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Vn gradient of the unknown vector V
    \param[out] Un result to D grad(v)
  */  
  template<class T>
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void GalbrunEquationH1_2D<T>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                       int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat,
                       int ref, Vector1& Vn, Vector2& Un)
  {
    abort();
  }
  

  //! Applies matrix M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un vector to be multiplied by M
    \param[out] Vn result vector V = M U
  */  
  template<class T>
  template<class TypeEquation, class T0, class Vector1>
  void GalbrunEquationH1_2D<T>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& Un, Vector1& Vn)
  {
    Vn.Fill(0);
    abort();
  }
  

  //! returns the matrix M, in the integral \f$ \int_K M u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Cj matrix M
  */  
  template<class T>
  template<class TypeEquation, class T0, class MatMass>
  void GalbrunEquationH1_2D<T>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& Cj)
  {
    Cj.Fill(0);    
    Real_wp rho = vars.eval_rho0(i)(j);
    Real_wp sigma = vars.eval_sigma(i)(j);
    T m_omega2; vars.GetMomega2(m_omega2);
    T m_iomega; vars.GetMiomega(m_iomega);
    Cj(0, 0) = rho*(m_omega2*nat_mat.GetCoefMass()
                    + sigma*(2.0*m_iomega*nat_mat.GetCoefDamping() + sigma*nat_mat.GetCoefStiffness()));
    
    Cj(1, 1) = Cj(0, 0);
    
    if (nb_unknowns == 4)
      {
	Cj(2, 2) = rho*(m_iomega*nat_mat.GetCoefDamping() + sigma*nat_mat.GetCoefStiffness());
	Cj(3, 3) = Cj(2, 2);
	Cj(2, 3) = -rho*nat_mat.GetCoefStiffness();
      }
    else
      Cj(2, 2) = Cj(0, 0);
    
  }
  

  //! computation of matrix N associated to the boundary condition
  /*!
    \param[out] Nabc matrix N
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] iquad element number
    \param[in] npoint local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */  
  template<class T>
  template<class Matrix1, class GenericPb, class T0>
  void GalbrunEquationH1_2D<T>::
  GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
	  int ref, int iquad, int npoint,
          const GlobalGenericMatrix<T0>& nat_mat, int ref2, const GenericPb& vars,
	  const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
  }
  
  
  //! mutliplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] num_elem1 element number
    \param[in] npoint local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class T>
  template<class Vector1, class TypeEquation, class T0>
  void GalbrunEquationH1_2D<T>::
  MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
	  int num_elem1, int npoint, const GlobalGenericMatrix<T0>& nat_mat, int ref2,
	  const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Un.Fill(0);
    abort();
  }
  
  
  //! Computation of "penalization" matrices C
  /*!
    \param[out] Nabc penalization matrix C
    \param[in] normale outward normale
    \param[in] iquad element number
    \param[in] npoint local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class T>
  template<class Matrix1, class GenericPb, class T0>
  void GalbrunEquationH1_2D<T>::
  GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int npoint, int nf,
             const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2, const GenericPb& vars,
	     const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
    Nabc(0, 0) = vars.delta_penalization*nat_mat.GetCoefStiffness()*vars.eval_rhoC2(iquad)(npoint);
    Nabc(1, 1) = Nabc(0, 0);
    
    if (nb_unknowns == 4)
      {
	Real_wp v0_dot_n = DotProd(normale, vars.eval_flow(iquad)(npoint));
	Nabc(2, 2) = abs(v0_dot_n)*nat_mat.GetCoefStiffness()*vars.alpha_penalization;
	Nabc(2, 2) *= vars.eval_rho0(iquad)(npoint)/vars.Glob_CoefPenalDG(nf);
	Nabc(3, 3) = Nabc(2, 2);
      }
  }
  

  //! Multiplication by penalization matrices
  /*!
    \param[in] normale outward normale
    \param[in] Vn vector to multiply
    \param[out] Un result vector Un = C*Vn
    \param[in] iquad element number
    \param[in] npoint local quadrature point number    
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */  
  template<class T>
  template<class Vector1, class Vector2, class GenericPb, class T0>
  void GalbrunEquationH1_2D<T>::
  MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
	     int iquad, int npoint, int nf, const GlobalGenericMatrix<T0>& nat_mat,
             int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Un.Fill(0);    
    abort();
  }


  //! volumetric source
  template<class T> 
  bool VolumetricSource_GalbrunH1<T>::IsNonNullVolumetricSource(const VectR2& s)
  {
    return true;
  }

  
  //! evaluation of volumetric source
  template<class T>
  void VolumetricSource_GalbrunH1<T>::EvaluateVolumetricSource(int i, int j, const R2& x, Vector<T>& f)
  {
    Real_wp aj; R2 grad_aj; TinyMatrix<Real_wp, Symmetric, 2, 2> hess_aj;
    this->fsrc.GetHessianAmplitude(x, aj, grad_aj, hess_aj);
    
    TinyVector<Real_wp, 2> polar;
    var_problem.GetPolarization(polar);
    
    Real_wp rot_f = -polar(0)*grad_aj(1) + polar(1)*grad_aj(0);
      
    if (var_galbrun_base.apply_convective_derivate_source)
      {
	T m_iomega; var_problem.GetMiomega(m_iomega);
	// part rot( ((-i omega + sigma) + M \cdot \nabla) f)
	T source_psi = (m_iomega + var_galbrun.eval_sigma(i)(j))*rot_f;
	source_psi += var_galbrun.grad_flow(i)(j)(0, 0)*polar(1)*grad_aj(0);
	source_psi += var_galbrun.grad_flow(i)(j)(1, 0)*polar(1)*grad_aj(1);
	source_psi -= var_galbrun.grad_flow(i)(j)(0, 1)*polar(0)*grad_aj(0);
	source_psi -= var_galbrun.grad_flow(i)(j)(1, 1)*polar(0)*grad_aj(1);
        
	source_psi += var_galbrun.eval_flow(i)(j)(0)*(hess_aj(0, 0)*polar(1) - hess_aj(0, 1)*polar(0));
	source_psi += var_galbrun.eval_flow(i)(j)(1)*(hess_aj(0, 1)*polar(1) - hess_aj(1, 1)*polar(0));
        
	T fx = (m_iomega + var_galbrun.eval_sigma(i)(j))*polar(0)*aj
	  + var_galbrun.eval_flow(i)(j)(0)*polar(0)*grad_aj(0)
	  + var_galbrun.eval_flow(i)(j)(1)*polar(0)*grad_aj(1);
	
	T fy = (m_iomega + var_galbrun.eval_sigma(i)(j))*polar(1)*aj
	  + var_galbrun.eval_flow(i)(j)(0)*polar(1)*grad_aj(0)
	  + var_galbrun.eval_flow(i)(j)(1)*polar(1)*grad_aj(1);
	
	// part 1/(rho c^2) f \times \nabla p0
	source_psi += 1.0/(var_galbrun.eval_rhoC2(i)(j)) 
	  * (fx*var_galbrun.grad_p0(i)(j)(1) - fy*var_galbrun.grad_p0(i)(j)(0));
	
	f(0) = fx;
	f(1) = fy;
	
	if (var_problem.nb_unknowns == 4)
	  {
	    f(2) = 0;
	    f(3) = source_psi;
	  }
	else
	  f(2) = source_psi;
	
      }
    else
      {
	// part 1/(rho c^2) f \times \nabla p0
	rot_f += 1.0/(var_galbrun.eval_rhoC2(i)(j)) 
	  * (polar(0)*var_galbrun.grad_p0(i)(j)(1) - polar(1)*var_galbrun.grad_p0(i)(j)(0)) * aj;
	
	f(0) = polar(0)*aj;
	f(1) = polar(1)*aj;
	if (var_problem.nb_unknowns == 4)
	  {
	    f(2) = 0;
	    f(3) = rot_f;
	  }
	else
	  f(2) = rot_f;
      }
  }


  /*******************
   * EllipticProblem *
   *******************/


  template<class Dimension>
  void EllipticProblem<HarmonicGalbrunEquationSipg<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicGalbrunEquationSipg<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
    

  template<class Dimension>
  void EllipticProblem<HarmonicGalbrunEquationSipg<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicGalbrunEquationSipg<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));
  }


  void EllipticProblem<HarmonicGalbrunEquationH1>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicGalbrunEquationH1>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    if (this->dg_formulation == ElementReference_Base::DISCONTINUOUS)
      Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
    

  void EllipticProblem<HarmonicGalbrunEquationH1>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void EllipticProblem<HarmonicGalbrunEquationH1>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));

    //A.Write("mat_elem.dat");
    //int test_input; cout << "waiting" << endl; cin >> test_input;
  }
  
}

#define MONTJOIE_FILE_SYMMETRIC_HARMONIC_GALBRUN_CXX
#endif
