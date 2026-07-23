#ifndef MONTJOIE_FILE_LINEARIZED_EULER_EQUATION_AXISYM_CXX

namespace Montjoie
{
  
  // for compatbility purpose
  template<class T> template<class TypeEquation>
  void AeroacousticAxiEquation_Base<T>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		      int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }
  
  
  
  template<class T> template<class GenericPb, class T0, class Vector1>
  void AeroacousticAxiEquation_Base<T>
  ::GetNeededDerivative(const GenericPb& vars,
			const GlobalGenericMatrix<T0>& nat_mat,
			Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    fct_test_to_derive.Fill(true);
    unknown_to_derive.Fill(true);
  }
  
  
  template<class T> template<class TypeEquation, class T0, class Vector1>
  void AeroacousticAxiEquation_Base<T>
  ::ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& dU, Vector1& dV)
  {
    abort();
  }
  
    
  template<class T> template<class TypeEquation, class T0, class MatStiff>
  void AeroacousticAxiEquation_Base<T>
  ::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
		      int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
		      int ref, MatStiff& Cgrad_grad)
  {
    FillZero(Cgrad_grad);
  }
  
  
  template<class T> template<class TypeEquation, class T0, class MatStiff>
  void AeroacousticAxiEquation_Base<T>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi); FillZero(Ephi_grad);
    
    T0 coef = nat_mat.GetCoefStiffness();
    T radius = vars.Glob_rtilde(num_elem)(jloc);
    
    if (vars.type_model == vars.CONSERVATIVE)
      {
        coef *= vars.eval_rho(num_elem)(jloc);
        
        // term + 1/r d/dr( r rho_0 c_0 u_r) + d/dz( rho_0 c_0 u_z)  in equation of p
        Dgrad_phi(0, 1)(0) = -radius*vars.eval_c0(num_elem)(jloc)*coef;
        Dgrad_phi(0, 3)(1) = Dgrad_phi(0, 1)(0);
        
        // term rho_0 c_0 \grad p in equation of u
        Ephi_grad(1, 0)(0) = -Dgrad_phi(0, 1)(0);
        Ephi_grad(3, 0)(1) = Ephi_grad(1, 0)(0);
      }
    else if (vars.type_model == vars.GALBRUN)
      {
        Real_wp rho0 = vars.eval_rho(num_elem)(jloc);
        Real_wp rho0_c02 = square(rho0*vars.eval_c0(num_elem)(jloc));
        
        // term  \rho_0^2 c_0^2 div_{rz} u_{rz} in equation of p
        Ephi_grad(0, 1)(0) = radius*rho0_c02*coef;
        Ephi_grad(0, 3)(1) = Ephi_grad(0, 1)(0);
        
        // term \grad p in equation of u
        //Ephi_grad(1, 0)(0) = radius*coef;
        //Ephi_grad(3, 0)(1) = Ephi_grad(1, 0)(0);
        Dgrad_phi(1, 0)(0) = -radius*coef;
        Dgrad_phi(3, 0)(1) = Dgrad_phi(1, 0)(0);
        
        T m_iomega; vars.GetMiomega(m_iomega);
        Real_wp sigma = vars.eval_sigma(num_elem)(jloc);
        T coef_dp0 = nat_mat.GetCoefStiffness() / (m_iomega + sigma);
        R2 dp0 = vars.grad_p0(num_elem)(jloc);
        
        // term du_z/dz dp0/dr / (-i omega + sigma) in equation of u_r
        Ephi_grad(1, 3)(1) = radius*coef_dp0*dp0(0);

        // term -du_z/dr dp0/dz / (-i omega + sigma) in equation of u_r
        Ephi_grad(1, 3)(0) = -radius*coef_dp0*dp0(1);        

        // term du_r/dr dp0/dz / (-i omega + sigma) in equation of u_z
        Ephi_grad(3, 1)(0) = radius*coef_dp0*dp0(1);

        // term -du_r/dz dp0/dr / (-i omega + sigma) in equation of u_z
        Ephi_grad(3, 1)(1) = -radius*coef_dp0*dp0(0);        
        
        // multiplication of coef by rho such that rho M \cdot \grad is computed after
        coef *= rho0;
      }
    else
      {
        Real_wp c02 = square(vars.eval_c0(num_elem)(jloc));
        
        // term 1/r d/dr( r c_0^2 u) + d/dz( c_0^2 u_z) in equation of p
        Dgrad_phi(0, 1)(0) = -radius*c02*coef;
        Dgrad_phi(0, 3)(1) = Dgrad_phi(0, 1)(0);
        
        // term \grad p in equation of u
        Ephi_grad(1, 0)(0) = radius*coef;
        Ephi_grad(3, 0)(1) = radius*coef;    
      }
    
    // flow term :  m_r dq/dr
    Ephi_grad(0, 0)(0) = radius*vars.eval_flow(num_elem)(jloc)(0)*coef;
    Ephi_grad(1, 1)(0) = Ephi_grad(0, 0)(0);
    Ephi_grad(2, 2)(0) = Ephi_grad(0, 0)(0);
    Ephi_grad(3, 3)(0) = Ephi_grad(0, 0)(0);
    
    // flow term m_z dq/dz
    Ephi_grad(0, 0)(1) = radius*vars.eval_flow(num_elem)(jloc)(2)*coef;
    Ephi_grad(1, 1)(1) = Ephi_grad(0, 0)(1);
    Ephi_grad(2, 2)(1) = Ephi_grad(0, 0)(1);
    Ephi_grad(3, 3)(1) = Ephi_grad(0, 0)(1);
  }
  
  
  template<class T> template<class TypeEquation,
			     class T0, class Vector1, class Vector2>
  void AeroacousticAxiEquation_Base<T>
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& dU, Vector2& V)
  {
    abort();
  }
  
  
  template<class T> template<class TypeEquation, class T0, class MatMass>
  void AeroacousticAxiEquation_Base<T>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    mass.Fill(0);
    
    Complex_wp im = Iwp*Real_wp(vars.GetCurrentModeNumber());
    T radius = vars.Glob_rtilde(num_elem)(jloc);
    Real_wp c02 = square(vars.eval_c0(num_elem)(jloc));
    Real_wp sigma = vars.eval_sigma(num_elem)(jloc);
    T0 coef = nat_mat.GetCoefStiffness(), s = coef;
    
    // term (-i omega + sigma) in all equations
    T m_iomega; vars.GetMiomega(m_iomega);
    mass(0, 0) = radius*(m_iomega*nat_mat.GetCoefMass() + sigma*nat_mat.GetCoefDamping());
    mass(1, 1) = mass(0, 0);
    mass(2, 2) = mass(0, 0);
    mass(3, 3) = mass(0, 0);

    // term -i m m_theta p/r in equation of p (and similar terms for other equations)
    coef *= Complex_wp(0, -vars.GetCurrentModeNumber())*vars.eval_flow(num_elem)(jloc)(1);    
    mass(0, 0) += coef;
    mass(1, 1) += coef;
    mass(2, 2) += coef;
    mass(3, 3) += coef;
    
    switch (vars.type_model)
      {
      case EllipticProblem<TypeEquation>::CONSERVATIVE :
        {
          // term -i m rho_0 c_0/r u_\theta in equation of p
          mass(0, 2) = -im*vars.eval_c0(num_elem)(jloc)*s;
          
          // term - i m rho_0 c_0 p/r in equation of u
          mass(2, 0) = mass(0, 2);
          
          mass *= vars.eval_rho(num_elem)(jloc);
        }
        break;
      case EllipticProblem<TypeEquation>::GALBRUN :
        {
          Real_wp rho0 = vars.eval_rho(num_elem)(jloc);
          Real_wp rho0_c02 = rho0*square(vars.eval_c0(num_elem)(jloc));
          T coef_dp0 = nat_mat.GetCoefStiffness() / (m_iomega + sigma);
          R2 dp0 = vars.grad_p0(num_elem)(jloc);
          
          // term -i m rho_0^2 c_0^2/r u_\theta in equation of p
          mass(0, 2) = -im*rho0_c02*s;
          
          // term \rho_0^2 c_0^2/r u_r in equation of p
          mass(0, 1) = rho0_c02*s;
                    
          mass *= rho0;
          
          // term coming from dp/dr = d(r p)/dr - p
          mass(1, 0) = -s;
          
          // term - i m  p/r in equation of u_\theta
          mass(2, 0) = -im*s;
          
          // term (u_r - im u_\theta)/(r (-i omega + sigma)) dp0/dr in equation of u_r
          mass(1, 1) += coef_dp0*dp0(0);
          mass(1, 2) -= im*coef_dp0*dp0(0);
          
          // term (u_r - im u_\theta)/(r (-i omega + sigma)) dp0/dz in equation of u_z
          mass(3, 1) += coef_dp0*dp0(1);
          mass(3, 2) -= im*coef_dp0*dp0(1);          
          
          // term im / (r (-i omega + sigma)) (u_r dp0/dr + u_z dp0/dz) in equation of u_\theta
          mass(2, 1) += im*coef_dp0*dp0(0);
          mass(2, 3) += im*coef_dp0*dp0(1);
        }
        break;
      case EllipticProblem<TypeEquation>::BOGEY_BAILLY_JUVE :
        {
          // term - i m c_0^2/r u_\theta in equation of p
          // put - abs(m) if the flow is null (such that the matrix is the same for -m and +m)
          mass(0, 2) = -im*c02*s;
          
          // term - i m/r p in equation of u_theta
          // put + abs(m) if the flow is null (such that the matrix is the same for -m and +m)
          mass(2, 0) = -im*s;
        }
        break;
      case EllipticProblem<TypeEquation>::SIMPLIFIED_LEE :
        {
          // term - i m c_0^2/r u_\theta in equation of p
          // put - abs(m) if the flow is null (such that the matrix is the same for -m and +m)
          mass(0, 2) = -im*c02*s;
          
          // term - i m/r p in equation of u_theta
          // put + abs(m) if the flow is null (such that the matrix is the same for -m and +m)
          mass(2, 0) = -im*s;
          
          // term grad(M) u in equation of u
          coef = s*radius;
          mass(1, 1) += coef*vars.grad_flow(num_elem)(jloc)(0, 0);
          mass(1, 3) += coef*vars.grad_flow(num_elem)(jloc)(0, 1);

          mass(2, 1) += coef*vars.grad_flow(num_elem)(jloc)(1, 0);
          mass(2, 3) += coef*vars.grad_flow(num_elem)(jloc)(1, 1);

          mass(3, 1) += coef*vars.grad_flow(num_elem)(jloc)(2, 0);
	  mass(3, 3) += coef*vars.grad_flow(num_elem)(jloc)(2, 1);
        }
      }
  }
  
  
  template<class T> template<class TypeEquation, class T0, class Vector1>
  void AeroacousticAxiEquation_Base<T>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		    const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    abort();
  }
  

  template<class T> template<class T0>
  void AeroacousticAxiEquation_Base<T>
  ::GetAbsoluteD(TinyMatrix<T0, General, 4, 4>& Dtest, const R2& normale, const Real_wp& alpha,
                 const Real_wp& c0, const T0& rtilde, bool conservative, bool galbrun, const T0& gamma)
  {
    // axisymmetric case
    Real_wp nr = normale(0), nz = normale(1);
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
            Dtest(0, 0) = s; Dtest(0, 1) = d*nr; Dtest(0, 2) = 0; Dtest(0, 3) = d*nz;
            Dtest(1, 0) = d*nr; Dtest(1, 1) = a + (s-a)*nr*nr; Dtest(1, 2) = 0; Dtest(1, 3) = nr*nz*(s-a);
            Dtest(2, 0) = 0; Dtest(2, 1) = 0; Dtest(2, 2) = a; Dtest(2, 3) = 0;
            Dtest(3, 0) = d*nz; Dtest(3, 1) = nr*nz*(s-a); Dtest(3, 2) = 0; Dtest(3, 3) = a + (s-a)*nz*nz;
          }
        else if (galbrun)
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nr*c0; Dtest(0, 2) = 0; Dtest(0, 3) = d*nz*c0;
            Dtest(1, 0) = d*nr/c0; Dtest(1, 1) = a + (s-a)*nr*nr; Dtest(1, 2) = 0; Dtest(1, 3) = nr*nz*(s-a);
            Dtest(2, 0) = 0; Dtest(2, 1) = 0; Dtest(2, 2) = a; Dtest(2, 3) = 0;
            Dtest(3, 0) = d*nz/c0; Dtest(3, 1) = nr*nz*(s-a); Dtest(3, 2) = 0; Dtest(3, 3) = a + (s-a)*nz*nz;
            
            // terms due to gamma / (-i omega + sigma)
            if (abs(alpha) > epsilon_machine)
              {
                Real_wp sa = sign(alpha), sAc = 0.5*(sign(alpha+c0)+sign(alpha-c0));
                Dtest(1, 1) += nr*nz*gamma*(sAc - sa); Dtest(1, 3) += gamma*(sa*nr*nr + sAc*nz*nz);
                Dtest(3, 1) -= gamma*(sa*nz*nz + sAc*nr*nr); Dtest(3, 3) += nr*nz*gamma*(sa-sAc);
              }
          }
        else
          {
            Dtest(0, 0) = s; Dtest(0, 1) = d*nr*c0; Dtest(0, 2) = 0; Dtest(0, 3) = d*nz*c0;
            Dtest(1, 0) = d*nr/c0; Dtest(1, 1) = a + (s-a)*nr*nr; Dtest(1, 2) = 0; Dtest(1, 3) = nr*nz*(s-a);
            Dtest(2, 0) = 0; Dtest(2, 1) = 0; Dtest(2, 2) = a; Dtest(2, 3) = 0;
            Dtest(3, 0) = d*nz/c0; Dtest(3, 1) = nr*nz*(s-a); Dtest(3, 2) = 0; Dtest(3, 3) = a + (s-a)*nz*nz;
          }
        
        Dtest *= rtilde;
      }
  }

  
  template<class T> template<class Matrix1, class GenericPb, class T0>
  void AeroacousticAxiEquation_Base<T>
  ::GetAbsoluteD(Matrix1& Nabc, const R2& normale, int iquad, int k,
                 const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars)
  {
    T0 gamma(0);
    bool conser = false, galbrun = false;
    Real_wp c0 = vars.eval_c0(iquad)(k);
    if (vars.type_model == vars.CONSERVATIVE)
      {
        c0 *= vars.eval_rho(iquad)(k);
        conser = true;
      }
    
    if (vars.type_model == vars.GALBRUN)
      {
        c0 *= vars.eval_rho(iquad)(k);
        T m_iomega; vars.GetMiomega(m_iomega);
        Real_wp sigma = vars.eval_sigma(iquad)(k);
        
        gamma = vars.grad_p0(iquad)(k)(0)*normale(1)
          -vars.grad_p0(iquad)(k)(1)*normale(0);
        
        gamma /= (m_iomega + sigma);
        galbrun = true;
      }

    T r = vars.Glob_rtilde(iquad)(k);
    Real_wp nr = normale(0), nz = normale(1);
    Real_wp alpha = vars.eval_flow(iquad)(k)(0)*nr + vars.eval_flow(iquad)(k)(2)*nz;
    if (galbrun || conser)
      alpha *= vars.eval_rho(iquad)(k);
    
    GetAbsoluteD(Nabc, normale, alpha, c0, r, conser, galbrun, gamma);
  }

  
  template<class T> template<class Matrix1, class GenericPb, class T0>
  void AeroacousticAxiEquation_Base<T>
  ::GetAbsoluteMatrixD(Matrix1& D, const R2& normale, int iquad, int k,
                       const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars)
  {
    T0 gamma(0);
    bool conser = false, galbrun = false;
    Real_wp c0 = vars.eval_c0(iquad)(k);
    if (vars.type_model == vars.CONSERVATIVE)
      {
        c0 *= vars.eval_rho(iquad)(k);
        conser = true;
      }
    
    if (vars.type_model == vars.GALBRUN)
      {
        c0 *= vars.eval_rho(iquad)(k);
        T m_iomega; vars.GetMiomega(m_iomega);
        Real_wp sigma = vars.eval_sigma(iquad)(k);
        
        gamma = vars.grad_p0(iquad)(k)(0)*normale(1)
          -vars.grad_p0(iquad)(k)(1)*normale(0);
        
        gamma /= (m_iomega + sigma);
        galbrun = true;
      }
    
    Real_wp c02 = c0*c0;
    
    T r = vars.Glob_rtilde(iquad)(k);
    Real_wp nr = normale(0), nz = normale(1);
    Real_wp alpha = vars.eval_flow(iquad)(k)(0)*nr + vars.eval_flow(iquad)(k)(2)*nz;
    if (galbrun || conser)
      alpha *= vars.eval_rho(iquad)(k);
    
    if (vars.type_model == vars.CONSERVATIVE)
      {
        D(0, 0) = alpha; D(0, 1) = c0*nr; D(0, 2) = 0; D(0, 3) = c0*nz;
        D(1, 0) = c0*nr; D(1, 1) = alpha; D(1, 2) = 0; D(1, 3) = 0;
        D(2, 0) = 0; D(2, 1) = 0; D(2, 2) = alpha; D(2, 3) = 0;
        D(3, 0) = c0*nz; D(3, 1) = 0; D(3, 2) = 0; D(3, 3) = alpha;
      }
    else
      {
        D(0, 0) = alpha; D(0, 1) = c02*nr; D(0, 2) = 0; D(0, 3) = c02*nz;
        D(1, 0) = nr; D(1, 1) = alpha; D(1, 2) = 0; D(1, 3) = 0;
        D(2, 0) = 0; D(2, 1) = 0; D(2, 2) = alpha; D(2, 3) = 0;
        D(3, 0) = nz; D(3, 1) = 0; D(3, 2) = 0; D(3, 3) = alpha;
      }        
    
    Mlt(r, D);
    
    GetAbsoluteD(D, normale, alpha, c0, r, conser, galbrun, gamma);        
  }
  
  
  template<class T> template<class Matrix1, class GenericPb, class T0>
  void AeroacousticAxiEquation_Base<T>
  ::GetNabc(Matrix1& Nabc, const R2& normale,
	    int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
	    const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
    T0 s = nat_mat.GetCoefStiffness();
    Real_wp c0 = vars.eval_c0(iquad)(k);
    T radius = vars.Glob_rtilde(iquad)(k);
    switch (vars.mesh.GetBoundaryCondition(ref))
      {
      case BoundaryConditionEnum::LINE_DIRICHLET :
	{
          if (vars.type_model == vars.CONSERVATIVE)
            {
              c0 *= vars.eval_rho(iquad)(k);
              Nabc(0, 1) = -radius*c0*normale(0)*s;
              Nabc(0, 3) = -radius*c0*normale(1)*s;
              
              Nabc(1, 0) = radius*c0*normale(0)*s;
              Nabc(3, 0) = radius*c0*normale(1)*s;
            }
          else if (vars.type_model == vars.GALBRUN)
            {
              c0 *= vars.eval_rho(iquad)(k);
              Nabc(0, 1) = -radius*c0*c0*normale(0)*s;
              Nabc(0, 3) = -radius*c0*c0*normale(1)*s;
              
              Nabc(1, 0) = radius*normale(0)*s;
              Nabc(3, 0) = radius*normale(1)*s;
            }
          else
            {
              Nabc(0, 1) = -radius*c0*c0*normale(0)*s;
              Nabc(0, 3) = -radius*c0*c0*normale(1)*s;
              
              Nabc(1, 0) = radius*normale(0)*s;
              Nabc(3, 0) = radius*normale(1)*s;
            }
	}
	break;
      case BoundaryConditionEnum::LINE_NEUMANN :
	{
          if (vars.type_model == vars.CONSERVATIVE)
            {
              c0 *= vars.eval_rho(iquad)(k);
              Nabc(0, 1) = radius*c0*normale(0)*s;
              Nabc(0, 3) = radius*c0*normale(1)*s;
              
              Nabc(1, 0) = -radius*c0*normale(0)*s;
              Nabc(3, 0) = -radius*c0*normale(1)*s;
            }
          else if (vars.type_model == vars.GALBRUN)
            {
              c0 *= vars.eval_rho(iquad)(k);
              Nabc(0, 1) = radius*c0*c0*normale(0)*s;
              Nabc(0, 3) = radius*c0*c0*normale(1)*s;
              
              Nabc(1, 0) = -radius*normale(0)*s;
              Nabc(3, 0) = -radius*normale(1)*s;
            }
          else
            {
              Nabc(0, 1) = radius*c0*c0*normale(0)*s;
              Nabc(0, 3) = radius*c0*c0*normale(1)*s;
              
              Nabc(1, 0) = -radius*normale(0)*s;
              Nabc(3, 0) = -radius*normale(1)*s;	  
            }
	}
	break;
      case BoundaryConditionEnum::LINE_ABSORBING :
	{
          GetAbsoluteD(Nabc, normale, iquad, k, nat_mat, vars);
          
          Nabc *= nat_mat.GetCoefStiffness();
	}
	break;
      }
  }
  
  
  template<class T> template<class Vector1, class TypeEquation, class T0>
  void AeroacousticAxiEquation_Base<T>
  ::MltNabc(const R2& normale, int ref,
	    const Vector1& Vn, Vector1& Un, int num_elem1,
	    int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
	    const EllipticProblem<TypeEquation>& vars,
	    const ElementReference<Dimension, 1>& )
  {
    abort();
  }
  
  
  template<class T> template<class Matrix1, class TypeEquation, class T0>
  void AeroacousticAxiEquation_Base<T>
  ::GetPenalDG(Matrix1& Nabc, const R2& normale, int iquad, int k,
	       int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	       const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb)
  {
    if (vars.upwind_fluxes)
      {
        GetAbsoluteD(Nabc, normale, iquad, k, nat_mat, vars);
        Nabc *= -nat_mat.GetCoefStiffness();
        //T radius = vars.Glob_radius(iquad)(k);
        //T r = vars.Glob_rtilde(iquad)(k);
        //Nabc *= radius/r;
        return;
      }

    Nabc.Fill(0);

    T0 s = nat_mat.GetCoefStiffness();
    s *= vars.alpha_penalization;
    if ((vars.type_model == vars.CONSERVATIVE) || (vars.type_model == vars.GALBRUN))
      s *= vars.eval_rho(iquad)(k);

    Real_wp c0 = vars.eval_c0(iquad)(k);
    T radius = vars.Glob_radius(iquad)(k);

    Nabc(0, 0) = radius*c0*s;
    Nabc(1, 1) = radius*c0*s*normale(0)*normale(0);
    Nabc(1, 3) = radius*c0*s*normale(0)*normale(1);
    Nabc(3, 1) = radius*c0*s*normale(0)*normale(1);
    Nabc(3, 3) = radius*c0*s*normale(1)*normale(1);    
  }
  
  
  template<class T> template<class Vector1, class Vector2, class GenericPb,
			     class T0>
  void AeroacousticAxiEquation_Base<T>
  ::MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
	       int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
	       int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    abort();
  }


  // for compatbility purpose
  template<class T> template<class TypeEquation>
  void LinearizedEulerEquation_Axi<T>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		      int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }
  
  
  
  template<class T> template<class GenericPb, class T0, class Vector1>
  void LinearizedEulerEquation_Axi<T>
  ::GetNeededDerivative(const GenericPb& vars,
			const GlobalGenericMatrix<T0>& nat_mat,
			Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    unknown_to_derive.Fill(true);
    fct_test_to_derive.Fill(true);
  }
  
  
  template<class T> template<class TypeEquation, class T0, class Vector1>
  void LinearizedEulerEquation_Axi<T>
  ::ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& dU, Vector1& dV)
  {
    abort();
  }
  
    
  template<class T> template<class TypeEquation, class T0, class MatStiff>
  void LinearizedEulerEquation_Axi<T>
  ::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
		      int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
		      int ref, MatStiff& Cgrad_grad)
  {
    FillZero(Cgrad_grad);
  }
  
  
  template<class T> template<class TypeEquation, class T0, class MatStiff>
  void LinearizedEulerEquation_Axi<T>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi); FillZero(Ephi_grad);
    
    T0 coef = nat_mat.GetCoefStiffness();
    T radius = vars.Glob_rtilde(num_elem)(jloc);
    Real_wp c02 = square(vars.eval_c0(num_elem)(jloc));
    
    // term 1/r d/dr( r c_0^2 u_r) + d/dz( c_0^2 u_z)  in equation of p
    Dgrad_phi(0, 1)(0) = -radius*c02*coef;
    Dgrad_phi(0, 3)(1) = Dgrad_phi(0, 1)(0);
    
    // term \grad p in equation of u
    Ephi_grad(1, 0)(0) = radius*coef;
    Ephi_grad(3, 0)(1) = radius*coef;    
    
    // flow term :  m_r dq/dr
    Ephi_grad(0, 0)(0) = radius*vars.eval_flow(num_elem)(jloc)(0)*coef;
    Ephi_grad(1, 1)(0) = Ephi_grad(0, 0)(0);
    Ephi_grad(2, 2)(0) = Ephi_grad(0, 0)(0);
    Ephi_grad(3, 3)(0) = Ephi_grad(0, 0)(0);
    Ephi_grad(4, 4)(0) = Ephi_grad(0, 0)(0);
    
    // flow term :  m_z dq/dz
    Ephi_grad(0, 0)(1) = radius*vars.eval_flow(num_elem)(jloc)(2)*coef;
    Ephi_grad(1, 1)(1) = Ephi_grad(0, 0)(1);
    Ephi_grad(2, 2)(1) = Ephi_grad(0, 0)(1);
    Ephi_grad(3, 3)(1) = Ephi_grad(0, 0)(1);
    Ephi_grad(4, 4)(1) = Ephi_grad(0, 0)(1);
    
    // term 1/r d/dr( r u_r) + d/dz( u_z) in equation of rho
    Dgrad_phi(4, 1)(0) = -radius*coef;
    Dgrad_phi(4, 3)(1) = Dgrad_phi(4, 1)(0);    
  }
  
  
  template<class T> template<class TypeEquation,
			     class T0, class Vector1, class Vector2>
  void LinearizedEulerEquation_Axi<T>
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& dU, Vector2& V)
  {
    abort();
  }
  
  
  template<class T> template<class TypeEquation, class T0, class MatMass>
  void LinearizedEulerEquation_Axi<T>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    mass.Zero();
    
    Complex_wp im = Iwp*Real_wp(vars.GetCurrentModeNumber());
    T radius = vars.Glob_rtilde(num_elem)(jloc);
    Real_wp c02 = square(vars.eval_c0(num_elem)(jloc));
    Real_wp sigma = vars.eval_sigma(num_elem)(jloc);
    T0 coef = nat_mat.GetCoefStiffness(), s = coef;
    
    // term (-i omega + sigma) in all equations
    T m_iomega; vars.GetMiomega(m_iomega);
    mass(0, 0) = radius*(m_iomega*nat_mat.GetCoefMass() + sigma*nat_mat.GetCoefDamping());
    mass(1, 1) = mass(0, 0);
    mass(2, 2) = mass(0, 0);
    mass(3, 3) = mass(0, 0);
    mass(4, 4) = mass(0, 0);

    // term -i m m_theta/r p in equation of p (and similar terms for other equations)
    coef *= Complex_wp(0, -vars.GetCurrentModeNumber())*vars.eval_flow(num_elem)(jloc)(1);    
    mass(0, 0) += coef;
    mass(1, 1) += coef;
    mass(2, 2) += coef;
    mass(3, 3) += coef;
    mass(4, 4) += coef;
    
    // term - i m c_0^2/r u_\theta in equation of p
    // put - abs(m) if the flow is null (such that the matrix is the same for -m and +m)
    mass(0, 2) = -im*c02*s;
    
    // term - i m/r p in equation of u_theta
    // put + abs(m) if the flow is null (such that the matrix is the same for -m and +m)
    mass(2, 0) = -im*s;

    // term \rho div M in equation of \rho
    coef = s*radius;
    mass(4, 4) += vars.div_flow(num_elem)(jloc)*coef;

    // u div M in equation of u
    // case were the relation \rho_0 div M + M \cdot \nabla rho_0 = 0 is true
    //Real_wp divM = vars.div_flow(num_elem)(jloc);

    // case where it is not true
    R2 Mrz(vars.eval_flow(num_elem)(jloc)(0), vars.eval_flow(num_elem)(jloc)(2));
    Real_wp divM = -DotProd(vars.grad_rho(num_elem)(jloc), Mrz) / vars.eval_rho(num_elem)(jloc);
    mass(1, 1) += divM*coef;
    mass(2, 2) += divM*coef;
    mass(3, 3) += divM*coef;
 
    // term -i m/r u_theta in equation of rho
    mass(4, 2) = -im*s;
    
    // term gamma (div M) p in equation of p
    mass(0, 0) += vars.eval_gamma(num_elem)(jloc)*vars.div_flow(num_elem)(jloc)*coef;
    
    // term -(gamma-1) / rho0 u \cdot \grad p0
    T0 coef_g(coef);
    coef_g *= (vars.eval_gamma(num_elem)(jloc)-1.0)/vars.eval_rho(num_elem)(jloc);
    mass(0, 1) = -coef_g*vars.grad_p0(num_elem)(jloc)(0);
    mass(0, 3) = -coef_g*vars.grad_p0(num_elem)(jloc)(1);
    
    // term - u \ cdot grad(gamma) p0/rho0 in equation of p
    T0 coefg = square(vars.eval_c0(num_elem)(jloc)) / vars.eval_gamma(num_elem)(jloc) * coef;
    mass(0, 1) -= coefg*vars.grad_gamma(num_elem)(jloc)(0);
    mass(0, 3) -= coefg*vars.grad_gamma(num_elem)(jloc)(1);

    // vec_u = (nabla M) M
    R3 vec_u;
    Mlt(vars.grad_flow(num_elem)(jloc), Mrz, vec_u);
    
    // term (nabla M) M rho in equation of u
    //mass(1, 4) = vec_u(0)*coef;
    //mass(2, 4) = vec_u(1)*coef;
    //mass(3, 4) = vec_u(2)*coef;
    mass(2, 4) = 0;
    mass(1, 4) = -coef*vars.grad_p0(num_elem)(jloc)(0)/vars.eval_rho(num_elem)(jloc);
    mass(3, 4) = -coef*vars.grad_p0(num_elem)(jloc)(1)/vars.eval_rho(num_elem)(jloc);

    // term grad(M) u in equation of u
    mass(1, 1) += coef*vars.grad_flow(num_elem)(jloc)(0, 0);
    mass(1, 3) += coef*vars.grad_flow(num_elem)(jloc)(0, 1);
    
    mass(2, 1) += coef*vars.grad_flow(num_elem)(jloc)(1, 0);
    mass(2, 3) += coef*vars.grad_flow(num_elem)(jloc)(1, 1);
    
    mass(3, 1) += coef*vars.grad_flow(num_elem)(jloc)(2, 0);
    mass(3, 3) += coef*vars.grad_flow(num_elem)(jloc)(2, 1);

    // AJOUT NATHAN
    // term coming from material derivatives in polar coordinates
    mass(1, 2) -= 2.0*s*vars.eval_flow(num_elem)(jloc)(1);
    mass(2, 1) += s*vars.eval_flow(num_elem)(jloc)(1);
    mass(2, 2) += s*vars.eval_flow(num_elem)(jloc)(0);


    //AJOUT NATHAN
    //gravity term in equations of u
    if(vars.compute_gravity)
    {
      //mass(1, 4) -= coef*vars.eval_gravity(num_elem)(jloc)(0);
      //mass(2, 4) -= coef*vars.eval_gravity(num_elem)(jloc)(1);
      //mass(3, 4) -= coef*vars.eval_gravity(num_elem)(jloc)(2);
    }
  }
  
  
  template<class T> template<class TypeEquation, class T0, class Vector1>
  void LinearizedEulerEquation_Axi<T>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		    const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    abort();
  }
  
  
  template<class T> template<class Matrix1, class GenericPb, class T0>
  void LinearizedEulerEquation_Axi<T>
  ::GetNabc(Matrix1& Nabc, const R2& normale,
	    int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
	    const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
    T0 s = nat_mat.GetCoefStiffness();
    Real_wp c0 = vars.eval_c0(iquad)(k);
    T radius = vars.Glob_rtilde(iquad)(k);
    switch (vars.mesh.GetBoundaryCondition(ref))
      {
      case BoundaryConditionEnum::LINE_DIRICHLET :
	{
	  Nabc(0, 1) = -radius*c0*c0*normale(0)*s;
	  Nabc(0, 3) = -radius*c0*c0*normale(1)*s;
	  
	  Nabc(4, 1) = -radius*normale(0)*s;
	  Nabc(4, 3) = -radius*normale(1)*s;
          
	  Nabc(1, 0) = radius*normale(0)*s;
	  Nabc(3, 0) = radius*normale(1)*s;
	}
	break;
      case BoundaryConditionEnum::LINE_NEUMANN :
	{
	  Nabc(0, 1) = radius*c0*c0*normale(0)*s;
	  Nabc(0, 3) = radius*c0*c0*normale(1)*s;
          
	  Nabc(4, 1) = radius*normale(0)*s;
	  Nabc(4, 3) = radius*normale(1)*s;
	  
	  Nabc(1, 0) = -radius*normale(0)*s;
	  Nabc(3, 0) = -radius*normale(1)*s;	  
	}
	break;
      case BoundaryConditionEnum::LINE_ABSORBING :
	{
	  s = nat_mat.GetCoefStiffness();
          
	  Nabc(0, 0) = radius*c0*s;
	  Nabc(4, 0) = radius*s/c0;
	  Nabc(1, 1) = radius*c0*s*normale(0)*normale(0);
	  Nabc(1, 3) = radius*c0*s*normale(0)*normale(1);
	  Nabc(3, 1) = radius*c0*s*normale(0)*normale(1);
	  Nabc(3, 3) = radius*c0*s*normale(1)*normale(1);
	}
	break;
      }
  }
  
  
  template<class T> template<class Vector1, class TypeEquation, class T0>
  void LinearizedEulerEquation_Axi<T>
  ::MltNabc(const R2& normale, int ref,
	    const Vector1& Vn, Vector1& Un, int num_elem1,
	    int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
	    const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& )
  {
    abort();
  }
  
  
  template<class T> template<class Matrix1, class TypeEquation, class T0>
  void LinearizedEulerEquation_Axi<T>
  ::GetPenalDG(Matrix1& Nabc, const R2& normale, int iquad, int k,
	       int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	       const EllipticProblem<TypeEquation>& vars,
	       const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);

    T0 s = nat_mat.GetCoefStiffness();
    s *= vars.alpha_penalization;
    Real_wp c0 = vars.eval_c0(iquad)(k);
    T radius = vars.Glob_radius(iquad)(k);

    Nabc(0, 0) = radius*c0*s;
    Nabc(4, 0) = radius*s/c0;
    Nabc(1, 1) = radius*c0*s*normale(0)*normale(0);
    Nabc(1, 3) = radius*c0*s*normale(0)*normale(1);
    Nabc(3, 1) = radius*c0*s*normale(0)*normale(1);
    Nabc(3, 3) = radius*c0*s*normale(1)*normale(1);    
  }
  
  
  template<class T> template<class Vector1, class Vector2, class GenericPb,
			     class T0>
  void LinearizedEulerEquation_Axi<T>
  ::MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
	       int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
	       int ref, int ref2, const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    abort();
  }


  /***********************
   * VarAeroacoustic_Axi *
   ***********************/
  

  //! parameters of the data file, specific to aero-acoustic equation
  /*!
    \param[in] description_field keyword of the considered line of the data file
    \param[in] parameters list of values associated
    \param[in] nb_param number of values
   */
  template<class Complexe>
  void VarAeroacoustic_Axi<Complexe>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarGalbrun_Axi<Complexe>::SetInputData(description_field, parameters);
    if (!description_field.compare("EnergyConservingAeroacousticModel"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarAeroacoustic_Axi" << endl;
	    cout << "EnergyConservingAeroacousticModel "
                 << "needs at least one parameter, for instance :" << endl;
	    cout << "EnergyConservingAeroacousticModel = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        if (var_problem.nb_unknowns == 5)
          {
            this->type_model = this->LEE_MODEL;
            return;
          }
        
        if (!parameters(0).compare("YES"))
          this->type_model = this->CONSERVATIVE;
        else
          {
            if (parameters(0) == "BogeyBaillyJuve")
              this->type_model = this->BOGEY_BAILLY_JUVE;
            else if (parameters(0) == "Galbrun")
              this->type_model = this->GALBRUN;
            else
              this->type_model = this->SIMPLIFIED_LEE;
          }
      }    
  }


  /*************************
   * VarAeroacousticAxi_Eq *
   *************************/
  

  template<class TypeEquation>  
  VarAeroacousticAxi_Eq<TypeEquation>::VarAeroacousticAxi_Eq()
    : VarHarmonic<TypeEquation>(), VarAeroacoustic_Axi<Complexe>(this->GetLeafClass())
  {
  }
  

  template<class TypeEquation>  
  const R3& VarAeroacousticAxi_Eq<TypeEquation>::GetPhaseOrigin() const
  {
    return VarAeroacoustic_Axi<Complexe>::GetPhaseOrigin();
  }

  
  //! returns true if numerical integration must be used
  template<class TypeEquation>  
  bool VarAeroacousticAxi_Eq<TypeEquation>::UseNumericalIntegration(int i) const
  {
    return true;
  }


  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::AddDiracSource(const Real_wp& alpha, Vector<Vector<Real_wp> > & b_source,
		   Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  //! adds Dirac for aeroacoustics equation
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::AddDiracSource(const Complex_wp& alpha, Vector<Vector<Complex_wp> > & b_source,
		   Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const
  {
    VectR3 point_source3D(1);
    point_source3D(0) = this->origine_phase3D;
    Vector<bool> vec_unknown(TypeEquation::nb_unknowns);
    vec_unknown.Fill(false);
    vec_unknown(1) = true; vec_unknown(2) = true;
    this->AddDiracSourceAxisym(alpha, vec_unknown, this->number_mode, b_source,
			       point_source3D, f);
  }


  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u,
			int i, const GridInterpolation<Dimension2>& var_interp,
			int iquad, bool compute_grad) const
  {
  }

  
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u,
			int i, const GridInterpolation<Dimension2>& var_interp,
			int iquad, bool compute_grad) const
  {
    Real_wp teta = var_interp.GetTheta(i);
    Complex_wp coef = exp(-Iwp*Complex_wp(this->number_mode)*teta);
    Real_wp cos_teta = cos(teta), sin_teta = sin(teta);

    for (int m = 0; m < val_u.GetM(); m++)
      val_u(m) *= coef;

    if (compute_grad)
      for (int m = 0; m < grad_u.GetM(); m++)
	grad_u(m) *= coef;
    
    Complex_wp Er = val_u(1), Eteta = val_u(2), Ez = val_u(3);
        
    // expression of E in cartesian coordinates
    val_u(1) = Er*cos_teta - Eteta*sin_teta;
    val_u(2) = Er*sin_teta + Eteta*cos_teta;
    val_u(3) = Ez;
  }


  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                             int, bool, bool) const
  {
  }
    
  
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                        int, bool, bool) const
  {
  }
  

  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    VarAeroacoustic_Axi<Complexe>::SetInputData(description_field, parameters);
  }


  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>::InitIndices(int n)
  {
    VarAeroacoustic_Axi<Complexe>::InitIndices(n);
  }
  

  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    VarAeroacoustic_Axi<Complexe>::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
						     compute_grad, compute_hess);
  }
  
  
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();
    VarAeroacoustic_Axi<Complexe>::ComputePhysicalCoefficients();
  }
  

  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>::AllocateMassMatrices()
  {
    VarGalbrun_Axi<Complexe>::AllocateMassMatrices();
  }
    
  
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
                           SetPoints<Dimension2>& PointsElem,
                           SetMatrices<Dimension2>& MatricesElem,
                           IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb)
  {
    VarHarmonic<TypeEquation>::ComputeLocalMassMatrix(i, N, linear_sparse, PointsElem,
                                                      MatricesElem, OrderFace, Fb);
    
    VarGalbrun_Axi<Complexe>::ComputeLocalMassMatrix(i);
  }
  
  
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    VarAeroacoustic_Axi<Complexe>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarAeroacoustic_Axi<Complexe>::SetPhysicalIndex(name_media, i, parameters);
  }
  

  template<class TypeEquation>
  string VarAeroacousticAxi_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarAeroacoustic_Axi<Complexe>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  int VarAeroacousticAxi_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarAeroacoustic_Axi<Complexe>::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  bool VarAeroacousticAxi_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarAeroacoustic_Axi<Complexe>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool VarAeroacousticAxi_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarAeroacoustic_Axi<Complexe>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  Real_wp VarAeroacousticAxi_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return VarAeroacoustic_Axi<Complexe>::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  Real_wp VarAeroacousticAxi_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return VarAeroacoustic_Axi<Complexe>::GetVelocityOfInfinity();
  }

  //! retrieve and treat referenced edges
  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>::PerformOtherInitializations()
  {        
    this->ComputeDofOnAxe(*this);
  }
  

  template<class TypeEquation>
  void VarAeroacousticAxi_Eq<TypeEquation>::CheckInputMesh()
  {
    VarAeroacoustic_Axi<Complexe>::CheckSectionMeshAxi();
  }


  //! volumetric source
  bool VolumetricSource_AxiAero
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    return true;
  }

  
  //! Evaluation of volumetric source f (term \int f \varphi)
  void VolumetricSource_AxiAero
  ::EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& fvec)
  {
    R3 pt3D;
    Real_wp f;
    
    pt3D.Init(x(0), 0, x(1));
    f = fsrc.GetAmplitude(pt3D);
    
    Complex_wp feval(f, 0);
    feval *= coef_vol*x(0);
    if (var_boundary.GetCurrentModeNumber() == 0)
      {
	fvec(0) = feval*polar(0);
        fvec(3) = feval*polar(3);    
      }
    else
      {
        Complex_wp fr(0, 0), ftheta(0, 0);
        if (var_boundary.GetCurrentModeNumber() == -1)
          {
            fr = 0.5*feval*(polar(1) - Iwp*polar(2));
            ftheta = 0.5*feval*(Iwp*polar(1) + polar(2));            
          }
        else if (var_boundary.GetCurrentModeNumber() == 1)
          {
            fr = 0.5*feval*(polar(1) + Iwp*polar(2));
            ftheta = 0.5*feval*(-Iwp*polar(1) + polar(2));
          }
	
	fvec(1) = fr;
	fvec(2) = ftheta;	
      }
  }


  void EllipticProblem<HarmonicLinearizedEulerEquationAxi>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void EllipticProblem<HarmonicLinearizedEulerEquationAxi>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
    

  void EllipticProblem<HarmonicLinearizedEulerEquationAxi>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }
  
  void EllipticProblem<HarmonicLinearizedEulerEquationAxi>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));
  }


  void EllipticProblem<HarmonicAeroacousticAxiEquation>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicAeroacousticAxiEquation>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
    

  void EllipticProblem<HarmonicAeroacousticAxiEquation>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void EllipticProblem<HarmonicAeroacousticAxiEquation>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
				      this->GetReferenceElementH1(iquad));
  }
  
};

#define MONTJOIE_FILE_LINEARIZED_EULER_EQUATION_AXISYM_CXX
#endif
