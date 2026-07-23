#ifndef MONTJOIE_FILE_AERO_ACOUSTIC_CXX

namespace Montjoie
{
  /***************************
   * LinearizedEulerEquation *
   ***************************/
  

  template<class T, class Dimension>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::SetIndexToCompute(AeroAcoustic_Base<Dimension>& var)
  {
    var.compute_grad_flow = true;
    var.compute_div_flow = true;
    var.compute_grad_rho = true;
    var.compute_grad_sigma = false;
    var.compute_grad_c0 = false;
    var.compute_hessian_flow = false;
    var.compute_hessian_p0 = false;
    var.compute_gamma = true;
    var.store_grad_rho0_c0 = false;
    var.compute_grad_gamma = true;
  }
  
  
  //! for compatbility purpose
  template<class T, class Dimension> template<class TypeEquation>
  void LinearizedEulerEquation_Base<T, Dimension>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                    int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }
  
  
  //! which derivatives to evaluate during matrix-vector product ?
  template<class T, class Dimension> template<class GenericPb, class T0, class Vector1>
  void LinearizedEulerEquation_Base<T, Dimension>::
  GetNeededDerivative(const GenericPb& vars, const GlobalGenericMatrix<T0>& nat_mat, 
                      Vector1& unknown_to_derive, Vector1& fct_test_to_derive) 
  {
    // only function-test related to p is derived
    fct_test_to_derive.Fill(false);
    fct_test_to_derive(0) = true;
    
    unknown_to_derive.Fill(true);
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
    \int_K E \nabla u v + D u \nabla v dx
    of the variational formulation
   */
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class MatStiff>
  void LinearizedEulerEquation_Base<T, Dimension>::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
        
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        // term div(c_0^2 u) in equation of p
        Dgrad_phi(0, i+2)(i) = -square(vars.eval_c0(num_elem)(jloc));
        
        // term div u in equation of rho
        Ephi_grad(1, i+2)(i) = 1.0;
        
        // term grad(p) in equation of u
        Ephi_grad(i+2, 0)(i) = 1.0;
        
        // flow M \cdot \nabla in all equations
        for (int k = 0; k < nb_unknowns; k++)
          Ephi_grad(k, k)(i) = vars.eval_flow(num_elem)(jloc)(i);
      }
    
    Dgrad_phi *= nat_mat.GetCoefStiffness();
    Ephi_grad *= nat_mat.GetCoefStiffness();
  }
  
  
  //! Applying the tensor E to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dVn gradient of the unknown vector V
    \param[out] Un result to E grad(v)
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                       int ref, const TinyVector<TinyVector<T0, 2>, 4>& dVn,
                       TinyVector<T0, 4>& Un)
  {
    // 2-D case
    T1 s = nat_mat.GetCoefStiffness();
    typename Dimension::R_N M = var.eval_flow(i)(j);
    
    // flow in p
    Un(0) = s*(M(0)*dVn(0)(0) + M(1)*dVn(0)(1));
    
    // term div u in equation of rho (and flow)
    Un(1) = s*(dVn(2)(0) + dVn(3)(1) + M(0)*dVn(1)(0) + M(1)*dVn(1)(1));
    
    // term grad(p) in equation of u (and flow)
    Un(2) = s*(dVn(0)(0) + M(0)*dVn(2)(0) + M(1)*dVn(2)(1));
    Un(3) = s*(dVn(0)(1) + M(0)*dVn(3)(0) + M(1)*dVn(3)(1));
  }


  //! Applying the tensor E to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dVn gradient of the unknown vector V
    \param[out] Un result to E grad(v)
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                       int ref, const TinyVector<TinyVector<T0, 3>, 5>& dVn,
                       TinyVector<T0, 5>& Un)
  {
    // 3-D case
    T1 s = nat_mat.GetCoefStiffness();
    typename Dimension::R_N M = var.eval_flow(i)(j);
    
    // flow in p
    Un(0) = s*(M(0)*dVn(0)(0) + M(1)*dVn(0)(1) + M(2)*dVn(0)(2));
    
    // term div u in equation of rho (and flow)
    Un(1) = s*(dVn(2)(0) + dVn(3)(1) + dVn(4)(2) + M(0)*dVn(1)(0) + M(1)*dVn(1)(1) + M(2)*dVn(1)(2));
    
    // term grad(p) in equation of u (and flow)
    Un(2) = s*(dVn(0)(0) + M(0)*dVn(2)(0) + M(1)*dVn(2)(1) + M(2)*dVn(2)(2));
    Un(3) = s*(dVn(0)(1) + M(0)*dVn(3)(0) + M(1)*dVn(3)(1) + M(2)*dVn(3)(2));
    Un(4) = s*(dVn(0)(2) + M(0)*dVn(4)(0) + M(1)*dVn(4)(1) + M(2)*dVn(4)(2));
  }
  
  
  //! Applying the tensor D to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result D u
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                       int ref, const TinyVector<T0, 4>& Un,
                       TinyVector<TinyVector<T0, 2>, 4>& Vn)
  {
    // 2-D case
    Vn(1).Zero(); Vn(2).Zero(); Vn(3).Zero();
    
    // term div (c_0^2 u)
    Real_wp c0 = square(var.eval_c0(i)(j));
    T1 s = c0*nat_mat.GetCoefStiffness();
    Vn(0).Init(-s*Un(2), -s*Un(3));
  }
  

  //! Applying the tensor D to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result D u
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T1>& nat_mat,
                       int ref, const TinyVector<T0, 5>& Un,
                       TinyVector<TinyVector<T0, 3>, 5>& Vn)
  {
    // 3-D case
    Vn(1).Zero(); Vn(2).Zero(); Vn(3).Zero(); Vn(4).Zero();
    
    // term div (c_0^2 u)
    Real_wp c0 = square(var.eval_c0(i)(j));
    T1 s = c0*nat_mat.GetCoefStiffness();
    Vn(0).Init(-s*Un(2), -s*Un(3), -s*Un(4));
  }

  
  //! computes the mass tensor A
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class MatMass>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                  int iquad, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    // term -i omega + sigma on all equations
    T m_iomega;
    vars.GetMiomega(m_iomega);
    T0 coef;
    coef = m_iomega*nat_mat.GetCoefMass() + vars.eval_sigma(iquad)(jloc)*nat_mat.GetCoefDamping();
    mass.SetDiagonal(coef);
    
    // term rho div M in equation of rho
    mass(1, 1) += vars.div_flow(iquad)(jloc)*nat_mat.GetCoefStiffness();
    
    // term gamma (div M) p in equation of p
    mass(0, 0) += vars.eval_gamma(iquad)(jloc)*vars.div_flow(iquad)(jloc)*nat_mat.GetCoefStiffness();
    
    // vec_u = (nabla M) M
    typename Dimension::R_N vec_u;
    Mlt(vars.grad_flow(iquad)(jloc), vars.eval_flow(iquad)(jloc), vec_u);
    
    coef = (vars.eval_gamma(iquad)(jloc)-1.0)/vars.eval_rho0(iquad)(jloc)*nat_mat.GetCoefStiffness();    
    T0 coefg = square(vars.eval_c0(iquad)(jloc)) / vars.eval_gamma(iquad)(jloc) * nat_mat.GetCoefStiffness();
    
    // case where the relation \rho_0 div M + M \cdot \nabla rho_0 = 0 is used
    // Real_wp divM = vars.div_flow(iquad)(jloc);
    
    // case where this relation is not used
    Real_wp divM = -DotProd(vars.grad_rho0(iquad)(jloc), vars.eval_flow(iquad)(jloc)) / vars.eval_rho0(iquad)(jloc);
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        // term u div M in equation of u
        mass(2+i, 2+i) += divM*nat_mat.GetCoefStiffness();
        
        // term (\nabla M) u in equation of u
        for (int j = 0; j < Dimension:: dim_N; j++)
          mass(2+i, 2+j) += vars.grad_flow(iquad)(jloc)(i, j)*nat_mat.GetCoefStiffness();
        
        // term (nabla M) M rho in equation of u
        mass(2+i, 1) += vec_u(i)*nat_mat.GetCoefStiffness();
        //AJOUT NATHAN
        //mass(2+i,1) += vars.grad_p0(iquad)(jloc)(i)*nat_mat.GetCoefStiffness()/vars.eval_rho0(iquad)(jloc);

        // term -(gamma-1) / rho0 u \cdot \nabla p_0 in equation of p
        mass(0, 2+i) = -coef*vars.grad_p0(iquad)(jloc)(i);
        
        // term - u \ cdot grad(gamma) p0/rho0
        mass(0, 2+i) -= coefg*vars.grad_gamma(iquad)(jloc)(i);
      }
  }
  
  
  //! applies the mass tensor A
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                    const GlobalGenericMatrix<T1>& nat_mat, int ref, const TinyVector<T0, 4>& U,
                    TinyVector<T0, 4>& V)
  {
    // 2-D case
    T m_iomega;
    var.GetMiomega(m_iomega);
    T1 coef, s = nat_mat.GetCoefStiffness();
    coef = m_iomega*nat_mat.GetCoefMass() + var.eval_sigma(i)(j)*nat_mat.GetCoefDamping();
    
    // vec_u = (nabla M) M
    typename Dimension::R_N vec_u;
    Mlt(var.grad_flow(i)(j), var.eval_flow(i)(j), vec_u);
    
    Real_wp gamma = var.eval_gamma(i)(j);
    Real_wp divM = var.div_flow(i)(j);
    Real_wp rho = var.eval_rho0(i)(j);
    V(0) = coef*U(0) + s*(gamma*divM*U(0) - ((gamma-1.0)*(var.grad_p0(i)(j)(0)*U(2)+var.grad_p0(i)(j)(1)*U(3)) + (var.grad_gamma(i)(j)(0)*U(2) + var.grad_gamma(i)(j)(1)*U(3)) * square(var.eval_c0(i)(j)) / var.eval_gamma(i)(j)) / rho );
    V(1) = (coef + divM*s)*U(1);
    V(2) = coef*U(2) + s*(divM*U(2) + var.grad_flow(i)(j)(0, 0)*U(2) + var.grad_flow(i)(j)(0, 1)*U(3) + vec_u(0)*U(1));
    V(3) = coef*U(3) + s*(divM*U(3) + var.grad_flow(i)(j)(1, 0)*U(2) + var.grad_flow(i)(j)(1, 1)*U(3) + vec_u(1)*U(1));
    
    //AJOUT NATHAN
    //V(2) = coef*U(2) + s*(divM*U(2) + var.grad_flow(i)(j)(0, 0)*U(2) + var.grad_flow(i)(j)(0, 1)*U(3) + (var.grad_p0(i)(j)(0))*U(1))/rho;
    //V(3) = coef*U(3) + s*(divM*U(3) + var.grad_flow(i)(j)(1, 0)*U(2) + var.grad_flow(i)(j)(1, 1)*U(3) + (var.grad_p0(i)(j)(1))*U(1));
  }
  

  //! applies the mass tensor A
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                    const GlobalGenericMatrix<T1>& nat_mat, int ref, const TinyVector<T0, 5>& U,
                    TinyVector<T0, 5>& V)
  {
    // 3-D case
    T m_iomega;
    var.GetMiomega(m_iomega);
    T1 coef, s = nat_mat.GetCoefStiffness();
    coef = m_iomega*nat_mat.GetCoefMass() + var.eval_sigma(i)(j)*nat_mat.GetCoefDamping();
    
    // vec_u = (nabla M) M
    typename Dimension::R_N vec_u;
    Mlt(var.grad_flow(i)(j), var.eval_flow(i)(j), vec_u);
    
    Real_wp gamma = var.eval_gamma(i)(j);
    Real_wp divM = var.div_flow(i)(j);
    Real_wp rho = var.eval_rho0(i)(j);
    V(0) = coef*U(0) + s*gamma*divM*U(0)
      -s*(gamma-1)*(var.grad_p0(i)(j)(0)*U(2)
                    +var.grad_p0(i)(j)(1)*U(3)+var.grad_p0(i)(j)(2)*U(4))/rho;
    
    V(1) = (coef + divM*s)*U(1);
    V(2) = coef*U(2) + s*(divM*U(2)+var.grad_flow(i)(j)(0, 0)*U(2) + var.grad_flow(i)(j)(0, 1)*U(3)
                          + var.grad_flow(i)(j)(0, 2)*U(4) + vec_u(0)*U(1));
    V(3) = coef*U(3) + s*(divM*U(3)+var.grad_flow(i)(j)(1, 0)*U(2) + var.grad_flow(i)(j)(1, 1)*U(3)
                          + var.grad_flow(i)(j)(1, 2)*U(4) + vec_u(1)*U(1));
    V(4) = coef*U(4) + s*(divM*U(4)+var.grad_flow(i)(j)(2, 0)*U(2) + var.grad_flow(i)(j)(2, 1)*U(3)
                          + var.grad_flow(i)(j)(2, 2)*U(4) + vec_u(2)*U(1));
    //AJOUT NATHAN
    //V(2) = coef*U(2) + s*(divM*U(2)+var.grad_flow(i)(j)(0, 0)*U(2) + var.grad_flow(i)(j)(0, 1)*U(3)
                          //+ var.grad_flow(i)(j)(0, 2)*U(4) + var.grad_p0(i)(j)(0)*U(1))/rho;
    //V(3) = coef*U(3) + s*(divM*U(3)+var.grad_flow(i)(j)(1, 0)*U(2) + var.grad_flow(i)(j)(1, 1)*U(3)
                          //+ var.grad_flow(i)(j)(1, 2)*U(4) + var.grad_p0(i)(j)(1)*U(1))/rho;
    //V(4) = coef*U(4) + s*(divM*U(4)+var.grad_flow(i)(j)(2, 0)*U(2) + var.grad_flow(i)(j)(2, 1)*U(3)
                          //+ var.grad_flow(i)(j)(2, 2)*U(4) + var.grad_p0(i)(j)(2)*U(1))/rho;
  }

  
  //! fills matrix |D| from matrix D
  template<class T, class Dimension> template<class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  GetAbsoluteD(TinyMatrix<T0, General, 4, 4>& Dtest, const R2& normale, const Real_wp& c0)
  {
    // 2-D case
    Real_wp alpha = realpart(Dtest(0, 0));
    Real_wp nx = normale(0), ny = normale(1);
    Real_wp a = abs(alpha), d = 0.5*(abs(alpha+c0) - abs(alpha-c0)), s = 0.5*(abs(alpha+c0) + abs(alpha-c0));
    
    if (alpha >= c0)
      {
        // input matrix D contains only positive eigenvalues
        // Dtest = Dtest
      }
    else if (alpha <= -c0)
      Dtest *= -1.0;
    else
      {
        Dtest(0, 0) = s; Dtest(0, 1) = 0; Dtest(0, 2) = nx*c0*d; Dtest(0, 3) = ny*c0*d;
        Dtest(1, 0) = (s-a)/(c0*c0); Dtest(1, 1) = a; Dtest(1, 2) = nx*d/c0; Dtest(1, 3) = ny*d/c0;
        Dtest(2, 0) = nx*d/c0; Dtest(2, 1) = 0; Dtest(2, 2) = a+nx*nx*(s-a); Dtest(2, 3) = nx*ny*(s-a);
        Dtest(3, 0) = ny*d/c0; Dtest(3, 1) = 0; Dtest(3, 2) = nx*ny*(s-a); Dtest(3, 3) = a+ny*ny*(s-a);
      }    
  }


  //! fills matrix |D| from matrix D
  template<class T, class Dimension> template<class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  GetAbsoluteD(TinyMatrix<T0, General, 5, 5>& Dtest, const R3& normale, const Real_wp& c0)
  {
    // 3-D case
    Real_wp alpha = realpart(Dtest(0, 0));
    Real_wp nx = normale(0), ny = normale(1), nz = normale(2);
    Real_wp a = abs(alpha), d = 0.5*(abs(alpha+c0) - abs(alpha-c0)), s = 0.5*(abs(alpha+c0) + abs(alpha-c0));
    
    if (alpha >= c0)
      {
        // input matrix D contains only positive eigenvalues
        // Dtest = Dtest
      }
    else if (alpha <= -c0)
      Dtest *= -1.0;
    else
      {
        Dtest(0, 0) = s; Dtest(0, 1) = 0; Dtest(0, 2) = nx*c0*d; Dtest(0, 3) = ny*c0*d; Dtest(0, 4) = nz*c0*d;
        Dtest(1, 0) = (s-a)/(c0*c0); Dtest(1, 1) = a; Dtest(1, 2) = nx*d/c0; Dtest(1, 3) = ny*d/c0; Dtest(1, 4) = nz*d/c0;
        Dtest(2, 0) = nx*d/c0; Dtest(2, 1) = 0; Dtest(2, 2) = a+nx*nx*(s-a); Dtest(2, 3) = nx*ny*(s-a); Dtest(2, 4) = nx*nz*(s-a);
        Dtest(3, 0) = ny*d/c0; Dtest(3, 1) = 0; Dtest(3, 2) = nx*ny*(s-a); Dtest(3, 3) = a+ny*ny*(s-a); Dtest(3, 4) = ny*nz*(s-a);
        Dtest(4, 0) = nz*d/c0; Dtest(4, 1) = 0; Dtest(4, 2) = nx*nz*(s-a); Dtest(4, 3) = ny*nz*(s-a); Dtest(4, 4) = a+nz*nz*(s-a);
      }    
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
  template<class T, class Dimension>
  template<class Matrix1, class GenericPb, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
	  int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat,
	  int ref_d, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    T0 s = nat_mat.GetCoefStiffness();
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        Nabc.Fill(0);
        Real_wp c0 = square(vars.eval_c0(iquad)(k));
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            Nabc(0, 2+i) = -c0*normale(i)*s;
            Nabc(1, 2+i) = -normale(i)*s;
            Nabc(2+i, 0) = normale(i)*s;
          }
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
        Nabc.Fill(0);
        Real_wp c0 = square(vars.eval_c0(iquad)(k));
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            Nabc(0, 2+i) = c0*normale(i)*s;
            Nabc(1, 2+i) = normale(i)*s;
            Nabc(2+i, 0) = -normale(i)*s;
          }
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        GetAbsoluteD(Nabc, normale, vars.eval_c0(iquad)(k));
        
        Nabc *= nat_mat.GetCoefStiffness();
      }
  }
  
  
  //! multiplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] refd reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class T, class Dimension>
  template<class Vector1, class TypeEquation, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>::
  MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
          int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int refd,
	  const EllipticProblem<TypeEquation>& vars,
	  const ElementReference<Dimension, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Un.Fill(0);
    T0 s = nat_mat.GetCoefStiffness();
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        Real_wp c0 = square(vars.eval_c0(iquad)(k));
        typename Vector1::value_type v_dot_n(0);
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            v_dot_n += normale(i)*Vn(2+i);
            Un(2+i) = normale(i)*s*Vn(0);
          }
        
        v_dot_n *= s;
        Un(0) = -c0*v_dot_n;
        Un(1) = -v_dot_n;
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
        Real_wp c0 = square(vars.eval_c0(iquad)(k));
        typename Vector1::value_type v_dot_n(0);
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            v_dot_n += normale(i)*Vn(2+i);
            Un(2+i) = -normale(i)*s*Vn(0);
          }
        
        v_dot_n *= s;
        Un(0) = c0*v_dot_n;
        Un(1) = v_dot_n;
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        TinyMatrix<Real_wp, General, nb_unknowns, nb_unknowns> Dtest;
        Real_wp alpha = DotProd(vars.eval_flow(iquad)(k), normale);
        Real_wp c0 = vars.eval_c0(iquad)(k);
        Dtest.SetDiagonal(alpha);
        
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            Dtest(0, 2+i) = c0*c0*normale(i);
            Dtest(1, 2+i) = normale(i);
            Dtest(2+i, 0) = normale(i);
          }
        
        GetAbsoluteD(Dtest, normale, c0);

        Mlt(Dtest, Vn, Un);
        
        Un *= nat_mat.GetCoefStiffness();
      }
  }
  
  
  //! computes penalty tensor
  template<class T, class Dimension>
  template<class Matrix1, class GenericPb, class T0>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
               int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
               const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Real_wp c0 = vars.eval_c0(iquad)(k);
    if (vars.upwind_fluxes)
      {
        GetAbsoluteD(Nabc, normale, c0);
        Nabc *= -nat_mat.GetCoefStiffness();
        return;
      }
    
    Nabc.Fill(0);
    Nabc(0, 0) = c0*vars.alpha_penalization;
    Nabc(1, 0) = vars.alpha_penalization/c0;
    // we put n . n^* to have same terms than in acoustic equations
    // in acoustic equations, it is necessary because
    // in presence of hetereogeneities (mu discontinuous), only v.n is constant across the elements
    for (int i = 0; i < Dimension::dim_N; i++)
      for (int j = 0; j < Dimension::dim_N; j++)
	Nabc(i+2,j+2) = c0*vars.delta_penalization*normale(i)*normale(j);
    
    Mlt(nat_mat.GetCoefStiffness(), Nabc);
  }
  
  
  //! applies penalty tensor
  template<class T, class Dimension>
  template<class T0, class GenericPb, class T1>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::MltPenalDG(const typename Dimension::R_N& normale,
               const TinyVector<T0, 4>& Vn, TinyVector<T0, 4>& Un,
               int iquad, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
               int ref, int ref2, const GenericPb& vars,
	       const ElementReference<Dimension, 1>& Fb)
  {
    if (vars.upwind_fluxes)
      {
        TinyMatrix<Real_wp, General, 4, 4> Dtest;
        Real_wp alpha = DotProd(vars.eval_flow(iquad)(k), normale);
        Real_wp c0 = vars.eval_c0(iquad)(k);
        Dtest.SetDiagonal(alpha);
        
        Dtest(0, 2) = c0*c0*normale(0); Dtest(0, 3) = c0*c0*normale(1);
        Dtest(1, 2) = normale(0); Dtest(1, 3) = normale(1);
        Dtest(2, 0) = normale(0); Dtest(3, 0) = normale(1);
        
        GetAbsoluteD(Dtest, normale, c0);
        
        Mlt(Dtest, Vn, Un);
        Un *= -nat_mat.GetCoefStiffness();
        return;
      }
    
    Real_wp c0 = vars.eval_c0(iquad)(k);
    // 2-D case
    Un(0) = Vn(0)*vars.alpha_penalization*c0*nat_mat.GetCoefStiffness();
    Un(1) = Vn(1)*vars.alpha_penalization*nat_mat.GetCoefStiffness()/c0;
    
    T0 vloc = normale(0)*Vn(2) + normale(1)*Vn(3);
    vloc *= c0*vars.delta_penalization*nat_mat.GetCoefStiffness();
    Un(2) = vloc*normale(0);
    Un(3) = vloc*normale(1);
    
  }


  //! applies penalty tensor
  template<class T, class Dimension>
  template<class T0, class GenericPb, class T1>
  void LinearizedEulerEquation_Base<T, Dimension>
  ::MltPenalDG(const typename Dimension::R_N& normale,
               const TinyVector<T0, 5>& Vn, TinyVector<T0, 5>& Un,
               int iquad, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
               int ref, int ref2, const GenericPb& vars,
	       const ElementReference<Dimension, 1>& Fb)
  {
    if (vars.upwind_fluxes)
      {
        TinyMatrix<Real_wp, General, 2+Dimension::dim_N, 2+Dimension::dim_N> Dtest;
        Real_wp alpha = DotProd(vars.eval_flow(iquad)(k), normale);
        Real_wp c0 = vars.eval_c0(iquad)(k);
        Dtest.SetDiagonal(alpha);
        
        Dtest(0, 2) = c0*c0*normale(0); Dtest(0, 3) = c0*c0*normale(1); Dtest(0, 4) = c0*c0*normale(2);
        Dtest(1, 2) = normale(0); Dtest(1, 3) = normale(1); Dtest(1, 4) = normale(2);
        Dtest(2, 0) = normale(0); Dtest(3, 0) = normale(1); Dtest(4, 0) = normale(2);
        
        GetAbsoluteD(Dtest, normale, c0);
        
        Mlt(Dtest, Vn, Un);
        Un *= -nat_mat.GetCoefStiffness();
        return;
      }

    // 3-D case
    Un(0) = Vn(0)*vars.alpha_penalization*nat_mat.GetCoefStiffness();
    Un(1) = Vn(1)*vars.alpha_penalization*nat_mat.GetCoefStiffness();
    
    T0 vloc = normale(0)*Vn(2) + normale(1)*Vn(3) + normale(2)*Vn(4);
    vloc *= vars.delta_penalization*nat_mat.GetCoefStiffness();
    Un(2) = vloc*normale(0);
    Un(3) = vloc*normale(1);
    Un(4) = vloc*normale(2);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, StationaryLinearizedEulerEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddVectorH1(Real_wp(1), nat_mat, trans, lvl, *this, X, Real_wp(1), Y, false);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, StationaryLinearizedEulerEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
  

  /**************************
   * AeroStationaryEquation *
   **************************/


  template<class T, class Dimension>
  void AeroStationaryEquation_Base<T, Dimension>
  ::SetIndexToCompute(AeroAcoustic_Base<Dimension>& var)
  {
    var.compute_grad_flow = false;    
    if (var.type_model == var.SIMPLIFIED_LEE)
      var.compute_grad_flow = true;
   
    var.compute_div_flow = true;
    var.compute_grad_rho = false;
    var.compute_grad_sigma = false;
    var.compute_grad_c0 = true;
    var.compute_hessian_flow = false;
    var.compute_hessian_p0 = false;
    var.compute_gamma = false;
    var.compute_grad_gamma = false;
    if (var.type_model == var.CONSERVATIVE)
      {
        var.store_grad_rho0_c0 = true;
        var.compute_grad_rho = true;
      }
    else
      var.store_grad_rho0_c0 = false;
  }

  
  //! for compatibility purpose only
  template<class T, class Dimension> template<class TypeEquation>
  void AeroStationaryEquation_Base<T, Dimension>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                    int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }

  
  //! which derivatives to evaluate during matrix-vector product ?
  template<class T, class Dimension> template<class GenericPb, class T0, class Vector1>
  void AeroStationaryEquation_Base<T, Dimension>::
  GetNeededDerivative(const GenericPb& vars, const GlobalGenericMatrix<T0>& nat_mat, 
                      Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {    
    unknown_to_derive.Fill(true);
    if (vars.UseExactIntegrationElement())
      fct_test_to_derive.Fill(false);    
    else
      fct_test_to_derive.Fill(true);
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
    \int_K E \nabla u v + D u \nabla v dx
    of the variational formulation
   */
  template<class T, class Dimension>
  template<class TypeEquation, class T0, class MatStiff>
  void AeroStationaryEquation_Base<T, Dimension>::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);

    Real_wp c0 = vars.eval_c0(num_elem)(jloc);
    TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(num_elem)(jloc);
    
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          c0 *= vars.eval_rho0(num_elem)(jloc);
          v0 *= vars.eval_rho0(num_elem)(jloc);
          if (!vars.UseExactIntegrationElement())
            {
              for (int i = 0; i < Dimension::dim_N; i++)
                {
                  // acoustic part
                  Dgrad_phi(0, i+1)(i) = -c0*nat_mat.GetCoefStiffness();
                  Ephi_grad(i+1, 0)(i) = c0*nat_mat.GetCoefStiffness();
                  
                  // flow part
                  for (int k = 0; k < nb_unknowns; k++)
                    {
                      Dgrad_phi(k, k)(i) = -0.5*v0(i)*nat_mat.GetCoefStiffness();
                      Ephi_grad(k, k)(i) = 0.5*v0(i)*nat_mat.GetCoefStiffness();
                    }
                  
                }
            }              
          else
            {
              for (int i = 0; i < Dimension::dim_N; i++)
                {	    
                  // acoustic part
                  Ephi_grad(i+1, 0)(i) = c0*nat_mat.GetCoefStiffness();
                  Ephi_grad(0, i+1)(i) = c0*nat_mat.GetCoefStiffness();
                  
                  // flow part
                  for (int p = 0; p < (Dimension::dim_N+1); p++)
                    Ephi_grad(p, p)(i) = v0(i)*nat_mat.GetCoefStiffness();	    
                }
            }
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
        {
          c0 *= c0;
          if (!vars.UseExactIntegrationElement())
            {
              for (int i = 0; i < Dimension::dim_N; i++)
                {	    
                  // acoustic part
                  Ephi_grad(i+1, 0)(i) = nat_mat.GetCoefStiffness();
                  Dgrad_phi(0, i+1)(i) = -c0*nat_mat.GetCoefStiffness();
                  
                  // flow part
                  for (int p = 0; p < (Dimension::dim_N+1); p++)
                    {
                      Ephi_grad(p, p)(i) = 0.5*v0(i)*nat_mat.GetCoefStiffness();	    
                      Dgrad_phi(p, p)(i) = -0.5*v0(i)*nat_mat.GetCoefStiffness();	    
                    }
                } 
            }
          else
            {
              for (int i = 0; i < Dimension::dim_N; i++)
                {	    
                  // acoustic part
                  Ephi_grad(i+1, 0)(i) = nat_mat.GetCoefStiffness();
                  Ephi_grad(0, i+1)(i) = c0*nat_mat.GetCoefStiffness();
                  
                  // flow part
                  for (int p = 0; p < (Dimension::dim_N+1); p++)
                    Ephi_grad(p, p)(i) = v0(i)*nat_mat.GetCoefStiffness();	    
                } 
            }
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          // replacing c0 with (rho_0 c_0)^2
          c0 *= vars.eval_rho0(num_elem)(jloc);
          c0 *= c0;
          // and v0 by rho_0 v_0
          v0 *= vars.eval_rho0(num_elem)(jloc);

          TinyVector<Real_wp, Dimension::dim_N> dp0 = vars.grad_p0(num_elem)(jloc);
          
          T m_iomega; vars.GetMiomega(m_iomega);
          Real_wp sigma = vars.eval_sigma(num_elem)(jloc);
          T coef_dp0 = nat_mat.GetCoefStiffness() / (m_iomega + sigma);
          
          if (!vars.UseExactIntegrationElement())
            {
              for (int i = 0; i < Dimension::dim_N; i++)
                {	    
                  // acoustic part
                  //Ephi_grad(i+1, 0)(i) = nat_mat.GetCoefStiffness();
                  Dgrad_phi(i+1, 0)(i) = -nat_mat.GetCoefStiffness();
                  Ephi_grad(0, i+1)(i) = c0*nat_mat.GetCoefStiffness();                  
                                    
                  // flow part
                  for (int p = 0; p < (Dimension::dim_N+1); p++)
                    {
                      Ephi_grad(p, p)(i) = 0.5*v0(i)*nat_mat.GetCoefStiffness();	    
                      Dgrad_phi(p, p)(i) = -0.5*v0(i)*nat_mat.GetCoefStiffness();	    
                    }
                  
                  // part with grad p0
                  for (int p = 0; p < Dimension::dim_N; p++)
                    {
                      Ephi_grad(1+i, 1+p)(p) += dp0(i)*coef_dp0;
                      Ephi_grad(1+i, 1+p)(i) -= dp0(p)*coef_dp0;
                    }
                } 
            }
          else
            {
              for (int i = 0; i < Dimension::dim_N; i++)
                {	    
                  // acoustic part
                  Ephi_grad(i+1, 0)(i) = nat_mat.GetCoefStiffness();
                  Ephi_grad(0, i+1)(i) = c0*nat_mat.GetCoefStiffness();
                  
                  // flow part
                  for (int p = 0; p < (Dimension::dim_N+1); p++)
                    Ephi_grad(p, p)(i) = v0(i)*nat_mat.GetCoefStiffness();	    
                  
                  // part with grad p0
                  for (int p = 0; p < Dimension::dim_N; p++)
                    {
                      Ephi_grad(1+i, 1+p)(p) += dp0(i)*coef_dp0;
                      Ephi_grad(1+i, 1+p)(i) -= dp0(p)*coef_dp0;
                    }
                } 
            }
        }
      }
  }
  
  
  //! Applying the tensor E to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Vn gradient of the unknown vector V
    \param[out] Un result to E grad(v)
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& vars,
                       int iquad, int k, const GlobalGenericMatrix<T1>& nat_mat,
                       int ref, const TinyVector<TinyVector<T0, 2>, 3>& dVn,
                       TinyVector<T0, 3>& Un)
  {
    Un.Fill(0);
    // 2-D case
    T1 s = nat_mat.GetCoefStiffness();
    TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(iquad)(k);
    Real_wp c0 = vars.eval_c0(iquad)(k);
    
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          c0 *= vars.eval_rho0(iquad)(k);
          v0 *= vars.eval_rho0(iquad)(k);
          if (!vars.UseExactIntegrationElement())
            {
              // flow part
              Un(0) = 0.5*s*(v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1));
              Un(1) = 0.5*s*(v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1));
              Un(2) = 0.5*s*(v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1));
              
              // acoustic part
              Un(1) += s*c0*dVn(0)(0);
              Un(2) += s*c0*dVn(0)(1);
            }              
          else
            {
              Un(0) = s*(c0*(dVn(1)(0)+dVn(2)(1)) + v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1));
              Un(1) = s*(c0*dVn(0)(0) + v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1));
              Un(2) = s*(c0*dVn(0)(1) + v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1));
            }
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
        {
          c0 *= c0;
          if (!vars.UseExactIntegrationElement())
            {
              // flow part
              Un(0) = 0.5*s*(v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1));
              Un(1) = 0.5*s*(v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1));
              Un(2) = 0.5*s*(v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1));
              
              // acoustic part
              Un(1) += s*dVn(0)(0);
              Un(2) += s*dVn(0)(1);
            }              
          else
            {
              Un(0) = s*(c0*(dVn(1)(0)+dVn(2)(1)) + v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1));
              Un(1) = s*(dVn(0)(0) + v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1));
              Un(2) = s*(dVn(0)(1) + v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1));
            }
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          cout << "not implemented " << endl;
          abort();
        }
        break;
      }
  }
  

  //! Applying the tensor E to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Vn gradient of the unknown vector V
    \param[out] Un result to E grad(v)
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& vars,
                       int iquad, int k, const GlobalGenericMatrix<T1>& nat_mat,
                       int ref, const TinyVector<TinyVector<T0, 3>, 4>& dVn,
                       TinyVector<T0, 4>& Un)
  {
    Un.Fill(0);
    // 3-D case
    T1 s = nat_mat.GetCoefStiffness();
    TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(iquad)(k);
    Real_wp c0 = vars.eval_c0(iquad)(k);
    
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          c0 *= vars.eval_rho0(iquad)(k);
          v0 *= vars.eval_rho0(iquad)(k);
          if (!vars.UseExactIntegrationElement())
            {
              // flow part
              Un(0) = 0.5*s*(v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1) + v0(2)*dVn(0)(2));
              Un(1) = 0.5*s*(v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1) + v0(2)*dVn(1)(2));
              Un(2) = 0.5*s*(v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1) + v0(2)*dVn(2)(2));
              
              // acoustic part
              Un(1) += s*c0*dVn(0)(0);
              Un(2) += s*c0*dVn(0)(1);
              Un(3) += s*c0*dVn(0)(2);
            }              
          else
            {
              Un(0) = s*(c0*(dVn(1)(0)+dVn(2)(1)+dVn(3)(2))
                         + v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1)) + v0(2)*dVn(0)(2);
              Un(1) = s*(c0*dVn(0)(0) + v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1) + v0(2)*dVn(1)(2));
              Un(2) = s*(c0*dVn(0)(1) + v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1) + v0(2)*dVn(2)(2));
              Un(3) = s*(c0*dVn(0)(2) + v0(0)*dVn(3)(0) + v0(1)*dVn(3)(1) + v0(2)*dVn(3)(2));
            }
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
        {
          c0 *= c0;
          if (!vars.UseExactIntegrationElement())
            {
              // flow part
              Un(0) = 0.5*s*(v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1) + v0(2)*dVn(0)(2));
              Un(1) = 0.5*s*(v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1) + v0(2)*dVn(1)(2));
              Un(2) = 0.5*s*(v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1) + v0(2)*dVn(2)(2));
              Un(3) = 0.5*s*(v0(0)*dVn(3)(0) + v0(1)*dVn(3)(1) + v0(2)*dVn(3)(2));
              
              // acoustic part
              Un(1) += s*dVn(0)(0);
              Un(2) += s*dVn(0)(1);
              Un(3) += s*dVn(0)(2);
            }              
          else
            {
              Un(0) = s*(c0*(dVn(1)(0)+dVn(2)(1)+dVn(3)(2))
                         + v0(0)*dVn(0)(0) + v0(1)*dVn(0)(1) + v0(2)*dVn(0)(2));
              Un(1) = s*(dVn(0)(0) + v0(0)*dVn(1)(0) + v0(1)*dVn(1)(1) + v0(2)*dVn(1)(2));
              Un(2) = s*(dVn(0)(1) + v0(0)*dVn(2)(0) + v0(1)*dVn(2)(1) + v0(2)*dVn(2)(2));
              Un(3) = s*(dVn(0)(2) + v0(0)*dVn(3)(0) + v0(1)*dVn(3)(1) + v0(2)*dVn(3)(2));
            }
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          cout << "not implemented " << endl;
          abort();
        }
        break;
      }
  }

  
  //! Applying the tensor D to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result D u
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& vars,
                       int iquad, int k, const GlobalGenericMatrix<T1>& nat_mat, int ref,
                       const TinyVector<T0, 3>& Un, TinyVector<TinyVector<T0, 2>, 3>& Vn)
  {
    // 2-D case
    FillZero(Vn);
    
    if (!vars.UseExactIntegrationElement())
      {
        T1 s = nat_mat.GetCoefStiffness();
        TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(iquad)(k);
        Real_wp c0 = vars.eval_c0(iquad)(k);
        switch (vars.type_model)
          {
          case EllipticProblem<TypeEquation>::CONSERVATIVE :
            {
              c0 *= vars.eval_rho0(iquad)(k);
              v0 *= vars.eval_rho0(iquad)(k);
              Vn(0)(0) = -s*(0.5*v0(0)*Un(0) + c0*Un(1));
              Vn(0)(1) = -s*(0.5*v0(1)*Un(0) + c0*Un(2));
              Vn(1)(0) = -0.5*s*v0(0)*Un(1);
              Vn(1)(1) = -0.5*s*v0(1)*Un(1);
              Vn(2)(0) = -0.5*s*v0(0)*Un(2);
              Vn(2)(1) = -0.5*s*v0(1)*Un(2);
            }
            break;
          case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
          case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
            c0 *= c0;
            Vn(0)(0) = -s*(0.5*v0(0)*Un(0) + c0*Un(1));
            Vn(0)(1) = -s*(0.5*v0(1)*Un(0) + c0*Un(2));
            Vn(1)(0) = -0.5*s*v0(0)*Un(1);
            Vn(1)(1) = -0.5*s*v0(1)*Un(1);
            Vn(2)(0) = -0.5*s*v0(0)*Un(2);
            Vn(2)(1) = -0.5*s*v0(1)*Un(2);
            break;            
          case EllipticProblem<TypeEquation>::GALBRUN :
            {
              cout << "not implemented " << endl;
              abort();
            }
            break;
          } 
      }
  }
  

  //! Applying the tensor D to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result D u
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& vars,
                       int iquad, int k, const GlobalGenericMatrix<T1>& nat_mat, int ref,
                       const TinyVector<T0, 4>& Un, TinyVector<TinyVector<T0, 3>, 4>& Vn)
  {
    // 3-D case
    FillZero(Vn);
    
    if (!vars.UseExactIntegrationElement())
      {
        T1 s = nat_mat.GetCoefStiffness();
        TinyVector<Real_wp, Dimension::dim_N> v0 = vars.eval_flow(iquad)(k);
        Real_wp c0 = vars.eval_c0(iquad)(k);
        switch (vars.type_model)
          {
          case EllipticProblem<TypeEquation>::CONSERVATIVE :
            {
              c0 *= vars.eval_rho0(iquad)(k);
              v0 *= vars.eval_rho0(iquad)(k);
              Vn(0)(0) = -s*(0.5*v0(0)*Un(0) + c0*Un(1));
              Vn(0)(1) = -s*(0.5*v0(1)*Un(0) + c0*Un(2));
              Vn(0)(2) = -s*(0.5*v0(2)*Un(0) + c0*Un(3));
              
              Vn(1)(0) = -0.5*s*v0(0)*Un(1);
              Vn(1)(1) = -0.5*s*v0(1)*Un(1);
              Vn(1)(2) = -0.5*s*v0(2)*Un(1);
              
              Vn(2)(0) = -0.5*s*v0(0)*Un(2);
              Vn(2)(1) = -0.5*s*v0(1)*Un(2);
              Vn(2)(2) = -0.5*s*v0(2)*Un(2);

              Vn(3)(0) = -0.5*s*v0(0)*Un(3);
              Vn(3)(1) = -0.5*s*v0(1)*Un(3);
              Vn(3)(2) = -0.5*s*v0(2)*Un(3);
            }
            break;
          case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
          case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
            c0 *= c0;
            Vn(0)(0) = -s*(0.5*v0(0)*Un(0) + c0*Un(1));
            Vn(0)(1) = -s*(0.5*v0(1)*Un(0) + c0*Un(2));
            Vn(0)(2) = -s*(0.5*v0(2)*Un(0) + c0*Un(3));
            
            Vn(1)(0) = -0.5*s*v0(0)*Un(1);
            Vn(1)(1) = -0.5*s*v0(1)*Un(1);
            Vn(1)(2) = -0.5*s*v0(2)*Un(1);
            
            Vn(2)(0) = -0.5*s*v0(0)*Un(2);
            Vn(2)(1) = -0.5*s*v0(1)*Un(2);
            Vn(2)(2) = -0.5*s*v0(2)*Un(2);

            Vn(3)(0) = -0.5*s*v0(0)*Un(3);
            Vn(3)(1) = -0.5*s*v0(1)*Un(3);
            Vn(3)(2) = -0.5*s*v0(2)*Un(3);
            break;
          case EllipticProblem<TypeEquation>::GALBRUN :
            {
              cout << "not implemented " << endl;
              abort();
            }
            break;
          }
      }
  }

  
  //! Applies tensor A to a vector
  /*!
    \param[in] var given problem
    \param[in] iquad number of the element where M needs to be evaluated
    \param[in] k number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un vector to be multiplied by A
    \param[out] Vn result vector Vn = A Un
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& vars, int iquad, int k,
                  const GlobalGenericMatrix<T1>& nat_mat, int ref, const TinyVector<T0, 3>& Un, TinyVector<T0, 3>& Vn)
  {
    // 2-D case
    
    // coef =  -i omega + sigma
    T m_iomega;
    vars.GetMiomega(m_iomega);
    T1 coef;
    coef = m_iomega*nat_mat.GetCoefMass() + vars.eval_sigma(iquad)(k)*nat_mat.GetCoefDamping();
    
    T1 s = nat_mat.GetCoefStiffness();
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
        {
          Vn(0) = coef*Un(0);
          Vn(1) = coef*Un(1) + s*(vars.grad_flow(iquad)(k)(0, 0)*Un(1) + vars.grad_flow(iquad)(k)(0, 1)*Un(2));
          Vn(2) = coef*Un(2) + s*(vars.grad_flow(iquad)(k)(1, 0)*Un(1) + vars.grad_flow(iquad)(k)(1, 1)*Un(2));
          
          if (vars.UseExactIntegrationElement())
            Vn(0) += s*(vars.grad_c0(iquad)(k)(0)*Un(1)+vars.grad_c0(iquad)(k)(1)*Un(2));
          else
            {
              s *= 0.5*vars.div_flow(iquad)(k);
              Vn(0) -= s*Un(0);
              Vn(1) -= s*Un(1);
              Vn(2) -= s*Un(2);
            }
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
        {
          Vn(0) = coef*Un(0);
          Vn(1) = coef*Un(1);
          Vn(2) = coef*Un(2);
          if (vars.UseExactIntegrationElement())
            Vn(0) += s*(vars.grad_c0(iquad)(k)(0)*Un(1)+vars.grad_c0(iquad)(k)(1)*Un(2));
          else
            {
              s *= 0.5*vars.div_flow(iquad)(k);
              Vn(0) -= s*Un(0);
              Vn(1) -= s*Un(1);
              Vn(2) -= s*Un(2);
            }
        }
        break;
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          coef *= vars.eval_rho0(iquad)(k);
          Vn(0) = coef*Un(0);
          Vn(1) = coef*Un(1);
          Vn(2) = coef*Un(2);
          if (vars.UseExactIntegrationElement())
            Vn(0) += s*(vars.grad_c0(iquad)(k)(0)*Un(1)+vars.grad_c0(iquad)(k)(1)*Un(2));
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          cout << "not implemented " << endl;
          abort();
        }
        break;
      }
  }
  

  //! Applies tensor A to a vector
  /*!
    \param[in] var given problem
    \param[in] iquad number of the element where M needs to be evaluated
    \param[in] k number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un vector to be multiplied by A
    \param[out] Vn result vector Vn = A Un
  */
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& vars, int iquad, int k,
                  const GlobalGenericMatrix<T1>& nat_mat, int ref, const TinyVector<T0, 4>& Un, TinyVector<T0, 4>& Vn)
  {
    // 3-D case
    
    // coef =  -i omega + sigma
    T m_iomega;
    vars.GetMiomega(m_iomega);
    T1 coef;
    coef = m_iomega*nat_mat.GetCoefMass() + vars.eval_sigma(iquad)(k)*nat_mat.GetCoefDamping();
    
    T1 s = nat_mat.GetCoefStiffness();
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
        {
          Vn(0) = coef*Un(0);
          Vn(1) = coef*Un(1) + s*(vars.grad_flow(iquad)(k)(0, 0)*Un(1)
                                  + vars.grad_flow(iquad)(k)(0, 1)*Un(2)
                                  + vars.grad_flow(iquad)(k)(0, 2)*Un(3));
          
          Vn(2) = coef*Un(2) + s*(vars.grad_flow(iquad)(k)(1, 0)*Un(1)
                                  + vars.grad_flow(iquad)(k)(1, 1)*Un(2)
                                  + vars.grad_flow(iquad)(k)(1, 2)*Un(3));

          Vn(3) = coef*Un(3) + s*(vars.grad_flow(iquad)(k)(2, 0)*Un(1)
                                  + vars.grad_flow(iquad)(k)(2, 1)*Un(2)
                                  + vars.grad_flow(iquad)(k)(2, 2)*Un(3));
          
          if (vars.UseExactIntegrationElement())
            Vn(0) += s*(vars.grad_c0(iquad)(k)(0)*Un(1)
                        + vars.grad_c0(iquad)(k)(1)*Un(2)
                        + vars.grad_c0(iquad)(k)(2)*Un(3));
          else
            {
              s *= 0.5*vars.div_flow(iquad)(k);
              Vn(0) -= s*Un(0);
              Vn(1) -= s*Un(1);
              Vn(2) -= s*Un(2);
              Vn(3) -= s*Un(3);
            }
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
        {
          Vn(0) = coef*Un(0);
          Vn(1) = coef*Un(1);
          Vn(2) = coef*Un(2);
          Vn(3) = coef*Un(3);
          if (vars.UseExactIntegrationElement())
            Vn(0) += s*(vars.grad_c0(iquad)(k)(0)*Un(1)
                        + vars.grad_c0(iquad)(k)(1)*Un(2)
                        + vars.grad_c0(iquad)(k)(2)*Un(3));
          else
            {
              s *= 0.5*vars.div_flow(iquad)(k);
              Vn(0) -= s*Un(0);
              Vn(1) -= s*Un(1);
              Vn(2) -= s*Un(2);
              Vn(3) -= s*Un(3);
            }
        }
        break;
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          coef *= vars.eval_rho0(iquad)(k);
          Vn(0) = coef*Un(0);
          Vn(1) = coef*Un(1);
          Vn(2) = coef*Un(2);
          Vn(3) = coef*Un(3);
          if (vars.UseExactIntegrationElement())
            Vn(0) += s*(vars.grad_c0(iquad)(k)(0)*Un(1)
                        + vars.grad_c0(iquad)(k)(1)*Un(2)
                        + vars.grad_c0(iquad)(k)(2)*Un(3));
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          cout << "not implemented " << endl;
          abort();
        }
        break;
      }
  }


  //! returns the tensor A, in the integral \f$ \int_K A u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] iquad number of the element where A needs to be evaluated
    \param[in] k number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] mass matrix A
  */  
  template<class T, class Dimension>
  template<class TypeEquation, class T1, class MatMass>
  void AeroStationaryEquation_Base<T, Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int iquad, int k, const GlobalGenericMatrix<T1>& nat_mat, int ref, MatMass& mass)
  {
    mass.Fill(0);
    // term -i omega + sigma
    T m_iomega;
    vars.GetMiomega(m_iomega);
    T1 coef;
    coef = m_iomega*nat_mat.GetCoefMass() + vars.eval_sigma(iquad)(k)*nat_mat.GetCoefDamping();

    T1 s = nat_mat.GetCoefStiffness();
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE:
        {
          mass.SetDiagonal(coef);
          for (int i = 0; i < Dimension::dim_N; i++)
            {
              // term grad(M) u
              for (int j = 0; j < Dimension::dim_N; j++)
                mass(1+i, 1+j) += vars.grad_flow(iquad)(k)(i, j)*s;
              
              // term grad(c0^2) \cdot u
              if (vars.UseExactIntegrationElement())
                mass(0, 1+i) = vars.grad_c0(iquad)(k)(i)*s;
              else
                {
                  // subtracting (div M)/2
                  mass(1+i, 1+i) -= 0.5*vars.div_flow(iquad)(k)*s;
                }
            }
          
          if (!vars.UseExactIntegrationElement())
            mass(0, 0) -= 0.5*vars.div_flow(iquad)(k)*s;
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
        {
          mass.SetDiagonal(coef);
          for (int i = 0; i < Dimension::dim_N; i++)
            {
              // term grad(c0^2) \cdot u
              if (vars.UseExactIntegrationElement())
                mass(0, 1+i) = vars.grad_c0(iquad)(k)(i)*s;
              else
                {
                  // subtracting (div M)/2
                  mass(1+i, 1+i) -= 0.5*vars.div_flow(iquad)(k)*s;              
                }
            }
          
          if (!vars.UseExactIntegrationElement())
            mass(0, 0) -= 0.5*vars.div_flow(iquad)(k)*s;
          
        }
        break;
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          coef *= vars.eval_rho0(iquad)(k);
          mass.SetDiagonal(coef);
          if (vars.UseExactIntegrationElement())
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                // term grad(c0 rho0) \cdot v
                mass(0, i+1) = vars.grad_c0(iquad)(k)(i)*s;
              }
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          coef *= vars.eval_rho0(iquad)(k);
          mass.SetDiagonal(coef);          
        }
        break;
      }    
  }

  
  template<class T, class Dimension> template<class T0>
  void AeroStationaryEquation_Base<T, Dimension>
  ::GetAbsoluteD(TinyMatrix<T0, General, 3, 3>& Dtest, const R2& normale,
                 const Real_wp& c0, bool conservative, bool galbrun, const T0& gamma)
  {
    // 2-D case
    Real_wp alpha = realpart(Dtest(0, 0));
    Real_wp nx = normale(0), ny = normale(1);
    Real_wp a = abs(alpha), d = 0.5*(abs(alpha+c0) - abs(alpha-c0)), s = 0.5*(abs(alpha+c0) + abs(alpha-c0));
    
    if (alpha >= c0)
      {
        // input matrix D contains only positive eigenvalues
        // Dtest = Dtest
      }
    else if (alpha <= -c0)
      Dtest *= -1.0;
    else
      {
        if (conservative)
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nx; Dtest(0, 2) = d*ny;
            Dtest(1, 0) = d*nx; Dtest(1, 1) = a + (s-a)*nx*nx; Dtest(1, 2) = nx*ny*(s-a);
            Dtest(2, 0) = d*ny; Dtest(2, 1) = nx*ny*(s-a); Dtest(2, 2) = a + (s-a)*ny*ny;        
          }
        else if (galbrun)
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nx*c0; Dtest(0, 2) = d*ny*c0;
            Dtest(1, 0) = d*nx/c0; Dtest(1, 1) = a + (s-a)*nx*nx; Dtest(1, 2) = nx*ny*(s-a);
            Dtest(2, 0) = d*ny/c0; Dtest(2, 1) = nx*ny*(s-a); Dtest(2, 2) = a + (s-a)*ny*ny;        
            
            // terms due to gamma / (-i omega + sigma)
            if (abs(alpha) > epsilon_machine)
              {
                Real_wp sa = sign(alpha), sAc = 0.5*(sign(alpha+c0)+sign(alpha-c0));
                Dtest(1, 1) += nx*ny*gamma*(sAc - sa); Dtest(1, 2) += gamma*(sa*nx*nx + sAc*ny*ny);
                Dtest(2, 1) -= gamma*(sa*ny*ny + sAc*nx*nx); Dtest(2, 2) += nx*ny*gamma*(sa-sAc);
              }
          }
        else
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nx*c0; Dtest(0, 2) = d*ny*c0;
            Dtest(1, 0) = d*nx/c0; Dtest(1, 1) = a + (s-a)*nx*nx; Dtest(1, 2) = nx*ny*(s-a);
            Dtest(2, 0) = d*ny/c0; Dtest(2, 1) = nx*ny*(s-a); Dtest(2, 2) = a + (s-a)*ny*ny;
          }
      }
  }


  template<class T, class Dimension> template<class T0>
  void AeroStationaryEquation_Base<T, Dimension>
  ::GetAbsoluteD(TinyMatrix<T0, General, 4, 4>& Dtest, const R3& normale,
                 const Real_wp& c0, bool conservative, bool galbrun, const T0& gamma)
  {
    // 3-D case
    Real_wp alpha = realpart(Dtest(0, 0));
    Real_wp nx = normale(0), ny = normale(1), nz = normale(2);
    Real_wp a = abs(alpha), d = 0.5*(abs(alpha+c0) - abs(alpha-c0)), s = 0.5*(abs(alpha+c0) + abs(alpha-c0));
    
    if (alpha >= c0)
      {
        // input matrix D contains only positive eigenvalues
        // Dtest = Dtest
      }
    else if (alpha <= -c0)
      Dtest *= -1.0;
    else
      {
        if (conservative)
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nx; Dtest(0, 2) = d*ny; Dtest(0, 3) = d*nz;
            Dtest(1, 0) = d*nx; Dtest(1, 1) = a + (s-a)*nx*nx; Dtest(1, 2) = nx*ny*(s-a); Dtest(1, 3) = nx*nz*(s-a);
            Dtest(2, 0) = d*ny; Dtest(2, 1) = nx*ny*(s-a); Dtest(2, 2) = a + (s-a)*ny*ny; Dtest(2, 3) = ny*nz*(s-a);
            Dtest(3, 0) = d*nz; Dtest(3, 1) = nx*nz*(s-a); Dtest(3, 2) = ny*nz*(s-a); Dtest(3, 3) = a + (s-a)*nz*nz;            
          }
        else if (galbrun)
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nx*c0; Dtest(0, 2) = d*ny*c0; Dtest(0, 3) = d*nz*c0;
            Dtest(1, 0) = d*nx/c0; Dtest(1, 1) = a + (s-a)*nx*nx; Dtest(1, 2) = nx*ny*(s-a); Dtest(1, 3) = nx*nz*(s-a);
            Dtest(2, 0) = d*ny/c0; Dtest(2, 1) = nx*ny*(s-a); Dtest(2, 2) = a + (s-a)*ny*ny; Dtest(2, 3) = (s-a)*ny*nz;
            Dtest(3, 0) = d*nz/c0; Dtest(3, 1) = nx*nz*(s-a); Dtest(3, 2) = (s-a)*ny*nz; Dtest(3, 3) = a + (s-a)*nz*nz;
            
            // terms due to \grad p0 not implemented...
            if (abs(gamma) > epsilon_machine)
              {
                cout << "Case not implemented" << endl;
                abort();
              }
          }
        else
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nx*c0; Dtest(0, 2) = d*ny*c0; Dtest(0, 3) = d*nz*c0;
            Dtest(1, 0) = d*nx/c0; Dtest(1, 1) = a + (s-a)*nx*nx; Dtest(1, 2) = nx*ny*(s-a); Dtest(1, 3) = nx*nz*(s-a);
            Dtest(2, 0) = d*ny/c0; Dtest(2, 1) = nx*ny*(s-a); Dtest(2, 2) = a + (s-a)*ny*ny; Dtest(2, 3) = (s-a)*ny*nz;
            Dtest(3, 0) = d*nz/c0; Dtest(3, 1) = nx*nz*(s-a); Dtest(3, 2) = (s-a)*ny*nz; Dtest(3, 3) = a + (s-a)*nz*nz;
          }
      }
  }
  
  
  template<class T, class Dimension> template<class Matrix1, class GenericPb, class T0>
  void AeroStationaryEquation_Base<T, Dimension>
  ::GetAbsoluteD(Matrix1& Nabc, const typename Dimension::R_N& normale, int iquad, int k,
                 const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars)
  {
    T0 gamma(0);
    bool conser = false, galbrun = false;
    Real_wp c0 = vars.eval_c0(iquad)(k);
    if (vars.type_model == vars.CONSERVATIVE)
      {
        c0 *= vars.eval_rho0(iquad)(k);
        conser = true;
      }
    
    if (vars.type_model == vars.GALBRUN)
      {
        c0 *= vars.eval_rho0(iquad)(k);
        T m_iomega; vars.GetMiomega(m_iomega);
        Real_wp sigma = vars.eval_sigma(iquad)(k);
        
        if (Dimension::dim_N == 2)
          gamma = vars.grad_p0(iquad)(k)(0)*normale(1)
            -vars.grad_p0(iquad)(k)(1)*normale(0);
        
        gamma /= (m_iomega + sigma);
        galbrun = true;
      }
        
    GetAbsoluteD(Nabc, normale, c0, conser, galbrun, gamma);
  }


  template<class T, class Dimension> template<class Matrix1, class GenericPb, class T0>
  void AeroStationaryEquation_Base<T, Dimension>
  ::GetAbsoluteMatrixD(Matrix1& D, const typename Dimension::R_N& normale, int iquad, int k,
                       const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars)
  {
    T0 gamma(0);
    bool conser = false, galbrun = false;
    Real_wp c0 = vars.eval_c0(iquad)(k);
    if (vars.type_model == vars.CONSERVATIVE)
      {
        c0 *= vars.eval_rho0(iquad)(k);
        conser = true;
      }
    
    if (vars.type_model == vars.GALBRUN)
      {
        c0 *= vars.eval_rho0(iquad)(k);
        T m_iomega; vars.GetMiomega(m_iomega);
        Real_wp sigma = vars.eval_sigma(iquad)(k);
        
        if (Dimension::dim_N == 2)
          gamma = vars.grad_p0(iquad)(k)(0)*normale(1)
            -vars.grad_p0(iquad)(k)(1)*normale(0);
        
        gamma /= (m_iomega + sigma);
        galbrun = true;
      }
    
    Real_wp alpha = DotProd(vars.eval_flow(iquad)(k), normale);
    if (conser || galbrun)
      alpha *= vars.eval_rho0(iquad)(k);        
    
    D.SetDiagonal(alpha);        
    if (conser)
      {
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            D(0, i+1) = c0*normale(i);
            D(i+1, 0) = c0*normale(i);
          }
      }
    else
      {
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            D(0, i+1) = c0*c0*normale(i);
            D(i+1, 0) = normale(i);
          }
        
        if ((galbrun) && (Dimension::dim_N == 2))
          {
            D(1, 2) = gamma;
            D(2, 1) = -gamma;
          }              
      }

    GetAbsoluteD(D, normale, c0, conser, galbrun, gamma);    
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
    \param[in] vars considered problem
    \param[in] Fb finite element associated with the element
   */
  template<class T, class Dimension>
  template<class Matrix1, class GenericPb, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
	  int ref, int iquad, int npoint, const GlobalGenericMatrix<T0>& nat_mat,
          int ref2, const GenericPb& vars,
	  const ElementReference<Dimension, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    
    T0 s = nat_mat.GetCoefStiffness();
    Real_wp c0 = vars.eval_c0(iquad)(npoint);
    
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        Nabc.Fill(0);
        if (vars.type_model == vars.CONSERVATIVE)
          {
            c0 *= vars.eval_rho0(iquad)(npoint);
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Nabc(0, i+1) = -c0*normale(i)*s;
                Nabc(i+1, 0) = c0*normale(i)*s;	    
              }
          }
        else if (vars.type_model == vars.GALBRUN)
          {
            c0 *= vars.eval_rho0(iquad)(npoint);
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Nabc(0, i+1) = -c0*normale(i)*s;
                Nabc(i+1, 0) = normale(i)*s;	    
              }
          }
        else
          {
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Nabc(0, i+1) = -c0*normale(i)*s;
                Nabc(i+1, 0) = normale(i)*s;	    
              }
          }
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
        Nabc.Fill(0);
        if (vars.type_model == vars.CONSERVATIVE)
          {
            c0 *= vars.eval_rho0(iquad)(npoint);
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Nabc(0, i+1) = c0*normale(i)*s;
                Nabc(i+1, 0) = -c0*normale(i)*s;	    
              }
          }
        else if (vars.type_model == vars.GALBRUN)
          {
            c0 *= vars.eval_rho0(iquad)(npoint);
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Nabc(0, i+1) = c0*normale(i)*s;
                Nabc(i+1, 0) = -normale(i)*s;	    
              }
          }
        else
          {
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Nabc(0, i+1) = c0*normale(i)*s;
                Nabc(i+1, 0) = -normale(i)*s;	    
              }
          }
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        GetAbsoluteD(Nabc, normale, iquad, npoint, nat_mat, vars);
        
        Nabc *= nat_mat.GetCoefStiffness();
      }
  }
  
  
  //! multiplication by matrix N associated to the boundary condition
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
  template<class T, class Dimension>
  template<class Vector1, class TypeEquation, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
	  int num_elem1, int npoint, const GlobalGenericMatrix<T0>& nat_mat, int ref2,
	  const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Un.Fill(0);
    
    Real_wp c0 = vars.eval_c0(num_elem1)(npoint);
    
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        if (vars.type_model == vars.CONSERVATIVE)
          {
            c0 *= vars.eval_rho0(num_elem1)(npoint);
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Un(0) += -c0*Vn(i+1)*normale(i);
                Un(i+1) = c0*Vn(0)*normale(i);
              }
          }
        else if (vars.type_model == vars.GALBRUN)
          {
            c0 *= vars.eval_rho0(num_elem1)(npoint);
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Un(0) += -c0*normale(i)*Vn(i+1);
                Un(i+1) = normale(i)*Vn(0);
              }
          }
        else
          {
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Un(0) += -c0*Vn(i+1)*normale(i);
                Un(i+1) = Vn(0)*normale(i);
              }
          }
	
	Mlt(nat_mat.GetCoefStiffness(), Un);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	if (vars.type_model == vars.CONSERVATIVE)
          {
            c0 *= vars.eval_rho0(num_elem1)(npoint);
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Un(0) += c0*Vn(i+1)*normale(i);
                Un(i+1) = -c0*Vn(0)*normale(i);
              }
          }
        else if (vars.type_model == vars.GALBRUN)
          {
            c0 *= vars.eval_rho0(num_elem1)(npoint);
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Un(0) += c0*normale(i)*Vn(i+1);
                Un(i+1) = -normale(i)*Vn(0);
              }
          }
        else
          {
            c0 *= c0;
            for (int i = 0; i < Dimension::dim_N; i++)
              {
                Un(0) += c0*Vn(i+1)*normale(i);
                Un(i+1) = -Vn(0)*normale(i);
              }
          }
	
	Mlt(nat_mat.GetCoefStiffness(), Un);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        TinyMatrix<T0, General, nb_unknowns, nb_unknowns> D;
        GetAbsoluteMatrixD(D, normale, num_elem1, npoint, nat_mat, vars);
        
        Mlt(D, Vn, Un);
        
        Un *= nat_mat.GetCoefStiffness();
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
  template<class T, class Dimension>
  template<class Matrix1, class GenericPb, class T0>
  void AeroStationaryEquation_Base<T, Dimension>::
  GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k, int nf,
	     const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2, const GenericPb& vars,
	     const ElementReference<Dimension, 1>& Fb)
  {
    if (vars.upwind_fluxes)
      {
        GetAbsoluteD(Nabc, normale, iquad, k, nat_mat, vars);
        Nabc *= -nat_mat.GetCoefStiffness();
        return;
      }
    
    Nabc.Fill(0);
    Nabc(0, 0) = vars.alpha_penalization;
    // we put n . n^* to have same terms than in acoustic equations
    for (int i = 0; i < Dimension::dim_N; i++)
      for (int j = 0; j < Dimension::dim_N; j++)
	Nabc(i+1, j+1) = vars.delta_penalization*normale(i)*normale(j);
    
    if (vars.type_model == vars.CONSERVATIVE)
      Mlt(nat_mat.GetCoefStiffness()*vars.eval_c0(iquad)(k)*vars.eval_rho0(iquad)(k), Nabc);
    else if (vars.type_model == vars.GALBRUN)
      Mlt(nat_mat.GetCoefStiffness()*vars.eval_c0(iquad)(k)*vars.eval_rho0(iquad)(k), Nabc);
    else
      Mlt(nat_mat.GetCoefStiffness()*vars.eval_c0(iquad)(k), Nabc);
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
  template<class T, class Dimension>
  template<class T0, class GenericPb, class T1>
  void AeroStationaryEquation_Base<T, Dimension>::
  MltPenalDG(const typename Dimension::R_N& normale,
             const TinyVector<T0, 3>& Vn, TinyVector<T0, 3>& Un,
             int iquad, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
             int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    if (vars.upwind_fluxes)
      {
        TinyMatrix<T1, General, nb_unknowns, nb_unknowns> D;
        GetAbsoluteMatrixD(D, normale, iquad, k, nat_mat, vars);
        
        Mlt(D, Vn, Un);
        
        Un *= -nat_mat.GetCoefStiffness();
        return;
      }

    // 2-D case
    Un(0) = Vn(0)*vars.alpha_penalization*nat_mat.GetCoefStiffness();
    
    T0 vloc = normale(0)*Vn(1) + normale(1)*Vn(2);   
    vloc *= vars.delta_penalization*nat_mat.GetCoefStiffness();
    
    Un(1) = vloc*normale(0);
    Un(2) = vloc*normale(1);

    if (vars.type_model == vars.CONSERVATIVE)
      Mlt(vars.eval_c0(iquad)(k)*vars.eval_rho0(iquad)(k), Un);
    else if (vars.type_model == vars.GALBRUN)
      Mlt(vars.eval_c0(iquad)(k)*vars.eval_rho0(iquad)(k), Un);
    else
      Mlt(vars.eval_c0(iquad)(k), Un);    
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
  template<class T, class Dimension>
  template<class T0, class GenericPb, class T1>
  void AeroStationaryEquation_Base<T, Dimension>::
  MltPenalDG(const typename Dimension::R_N& normale,
             const TinyVector<T0, 4>& Vn, TinyVector<T0, 4>& Un,
             int iquad, int k, int nf, const GlobalGenericMatrix<T1>& nat_mat,
             int ref, int ref2, const GenericPb& vars,
	     const ElementReference<Dimension, 1>& Fb)
  {
    if (vars.upwind_fluxes)
      {
        TinyMatrix<T1, General, nb_unknowns, nb_unknowns> D;
        GetAbsoluteMatrixD(D, normale, iquad, k, nat_mat, vars);
        
        Mlt(D, Vn, Un);
        
        Un *= -nat_mat.GetCoefStiffness();
        return;
      }
    
    // 3-D case
    Un(0) = Vn(0)*vars.alpha_penalization*nat_mat.GetCoefStiffness();
    
    T0 vloc = normale(0)*Vn(1) + normale(1)*Vn(2) + normale(2)*Vn(3);   
    vloc *= vars.delta_penalization*nat_mat.GetCoefStiffness();
    
    Un(1) = vloc*normale(0);
    Un(2) = vloc*normale(1);
    Un(3) = vloc*normale(2);

    if (vars.type_model == vars.CONSERVATIVE)
      Mlt(vars.eval_c0(iquad)(k)*vars.eval_rho0(iquad)(k), Un);
    else if (vars.type_model == vars.GALBRUN)
      Mlt(vars.eval_c0(iquad)(k)*vars.eval_rho0(iquad)(k), Un);
    else
      Mlt(vars.eval_c0(iquad)(k), Un);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, AeroStationaryEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddVectorH1(Real_wp(1), nat_mat, trans, lvl, *this, X, Real_wp(1), Y, false);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, AeroStationaryEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  /*********************
   * AeroAcoustic_Base *
   *********************/
  
  
  //! parameters of the data file, specific to aero-acoustic equation
  /*!
    \param[in] description_field keyword of the considered line of the data file
    \param[in] parameters list of values associated
    \param[in] nb_param number of values
   */
  template<class Dimension>
  void AeroAcoustic_Base<Dimension>::
  SetInputData(const string& description_field, const VectString& parameters)
  {    
    VarGalbrunIndex_Base<Dimension>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("EnergyConservingAeroacousticModel"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of AeroAcoustic_Base" << endl;
	    cout << "EnergyConservingAeroacousticModel "
                 << "needs at least one parameter, for instance :" << endl;
	    cout << "EnergyConservingAeroacousticModel = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        if (var_problem.nb_unknowns == 2+Dimension::dim_N)
          {
            type_model = LEE_MODEL;
            return;
          }
        
        if (!parameters(0).compare("YES"))
          type_model = CONSERVATIVE;
        else
          {
            if (parameters(0) == "BogeyBaillyJuve")
              type_model = BOGEY_BAILLY_JUVE;
            else if (parameters(0) == "Galbrun")
              type_model = GALBRUN;
            else
              type_model = SIMPLIFIED_LEE;
          }
        //std::cout << "LEE MODEL : " << type_model << std::endl;//AJOUT NATHAN
      }
    
  }


  template<class Dimension> template<class T0>
  void AeroAcoustic_Base<Dimension>
  ::ModifyVolumetricSource(int i, int j, const R_N& x,
			   const VirtualSourceField<T0, Dimension>& fsrc,
                           Vector<T0>& f) const
  {
    //if ((this->type_model != GALBRUN) && (this->type_model != CONSERVATIVE))
    //f *= 1.0/this->eval_rho0(i)(j);
  }
  

  /*******************
   * AeroAcoustic_Eq *
   *******************/


  template<class TypeEquation>
  AeroAcoustic_Eq<TypeEquation>::AeroAcoustic_Eq()
    : VarHarmonic<TypeEquation>(),
      AeroAcoustic_Base<typename TypeEquation::Dimension>(this->GetLeafClass())
  {
  }
  

  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {    
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    AeroAcoustic_Base<Dimension>::SetInputData(description_field, parameters);
  }

    
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>::InitIndices(int n)
  {
    AeroAcoustic_Base<Dimension>::InitIndices(n);
  }


  template<class TypeEquation>
  int AeroAcoustic_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return AeroAcoustic_Base<Dimension>::GetNbPhysicalIndices();
  }


  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    AeroAcoustic_Base<Dimension>::SetIndices(i, parameters);
  }


  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    AeroAcoustic_Base<Dimension>::SetPhysicalIndex(name_media, i, parameters);
  }

  
  template<class TypeEquation>
  string AeroAcoustic_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return AeroAcoustic_Base<Dimension>::GetPhysicalIndexName(m);
  }


  template<class TypeEquation>
  bool AeroAcoustic_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return AeroAcoustic_Base<Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool AeroAcoustic_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return AeroAcoustic_Base<Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  Real_wp AeroAcoustic_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return AeroAcoustic_Base<Dimension>::GetVelocityOfMedia(ref);
  }


  template<class TypeEquation>
  Real_wp AeroAcoustic_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return AeroAcoustic_Base<Dimension>::GetVelocityOfInfinity();
  }

  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real,
		      IVect& num_ref, IVect& num_index, IVect& num_component,
		      Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    return AeroAcoustic_Base<Dimension>
      ::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
			  compute_grad, compute_hess);
  }

  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    TypeEquation::SetIndexToCompute(*this);
    
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();      
    AeroAcoustic_Base<Dimension>::ComputePhysicalCoefficients();
  }
  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>::PerformOtherInitializations()
  {
    SetComplexOne(this->coefficient_volumic_source);
  }
  

  template<class TypeEquation> template<class T, class Dim>
  void AeroAcoustic_Eq<TypeEquation>
  ::ModifyVolumetricSourceGen(int i, int j, const typename Dim::R_N&,
                              const VirtualSourceField<T, Dim>&,
                              Vector<T>&) const
  {
    cout << "Not implemented for this type" << endl;
    abort();
  }
  
  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::ModifyVolumetricSourceGen(int i, int j, const typename Dimension::R_N& x,
                              const VirtualSourceField<Complexe, Dimension>& fsrc,
                              Vector<Complexe>& f) const
  {
    AeroAcoustic_Base<Dimension>::ModifyVolumetricSource(i, j, x, fsrc, f);
  }

  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R2& pt,
                           const VirtualSourceField<Real_wp, Dimension2>& src,
                           Vector<Real_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }
    
  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R2& pt,
                           const VirtualSourceField<Complex_wp, Dimension2>& src,
                           Vector<Complex_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }

  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R3& pt,
                           const VirtualSourceField<Real_wp, Dimension3>& src,
                           Vector<Real_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }
    
  
  template<class TypeEquation>
  void AeroAcoustic_Eq<TypeEquation>
  ::ModifyVolumetricSource(int i, int j, const R3& pt,
                           const VirtualSourceField<Complex_wp, Dimension3>& src,
                           Vector<Complex_wp>& f) const
  {
    ModifyVolumetricSourceGen(i, j, pt, src, f);
  }
  

  /*******************
   * EllipticProblem *
   *******************/


  template<class Dimension>
  void EllipticProblem<HarmonicLinearizedEulerEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicLinearizedEulerEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
    

  template<class Dimension>
  void EllipticProblem<HarmonicLinearizedEulerEquation<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicLinearizedEulerEquation<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));
  }


  template<class Dimension>
  void EllipticProblem<HarmonicAeroEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicAeroEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
    

  template<class Dimension>
  void EllipticProblem<HarmonicAeroEquation<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicAeroEquation<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));
  }


  template<class Dimension>
  void EllipticProblem<AeroStationaryEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }


  template<class Dimension>
  void EllipticProblem<AeroStationaryEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not implemented" << endl;
    abort();
  }
    

  template<class Dimension>
  void EllipticProblem<AeroStationaryEquation<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));
  }

  
  template<class Dimension>
  void EllipticProblem<AeroStationaryEquation<Dimension> >
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }

}

#define MONTJOIE_FILE_AERO_ACOUSTIC_CXX
#endif
