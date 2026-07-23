#ifndef MONTJOIE_FILE_VAR_HELMHOLTZ_CXX

#include "VarHelmholtzInline.cxx"

namespace Montjoie
{

  //! default constructor
  template<class Complexe, class Dimension>
  VarHelmholtzIndex_Base<Complexe, Dimension>::VarHelmholtzIndex_Base()
  {
    SetComplexOne(rho0); SetComplexOne(mu0);
    SetComplexZero(sigma0);
    nb_dof_drude_vec = 0; nb_dof_drude_vec_all = 0;
    use_symm_drude = false; linearize_drude = false;
  }
  
  
  /*********************
   * HelmholtzEquation *
   *********************/
  
  
  //! computation of coefficients necessary to the computation of elementary matrices
  template<class T, class Dimension>
  void HelmholtzEquation_Base<T, Dimension>
  ::ComputeMassMatrixGen(VarProblem<Dimension>& var_problem,
			 VarBoundaryCondition_Dim<T, Dimension>& var_boundary,
			 VarHelmholtz_Base<T, Dimension>& var_helm,
			 int num_elem, const ElementReference<Dimension, 1>& Fb)
  {
    TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> mu, mat_stiff, invMu;
    T rho(1), sigma(0), momega2, m_iomega, beta;
    TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> mu_dfjm1;
    typename Dimension::MatrixN_N dfjm1;    
    TinyVector<T, Dimension::dim_N> v0, vec_u;
    int ref = var_problem.mesh.Element(num_elem).GetReference();
    bool variable = var_problem.UseNumericalIntegration(num_elem);
    bool affine = var_problem.mesh.IsElementAffine(num_elem);
    
    // if constant DFi and constant physical constants, we compute
    // geometric quantities only for a quadrature point
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;

    int i1 = num_elem - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
    Real_wp jacob;    
    // m_iomega = -i omega, momega2 = - omega^2 for complex numbers
    // m_iomega = 1, momega2 = 1 for real numbers
    var_problem.GetMiomega(m_iomega);
    var_problem.GetMomega2(momega2);
    
    if (!var_problem.FirstOrderFormulation())
      {
        // second-order formulation
        // solved equation : -omega^2 rho u - i omega sigma u + div(mu grad u) = f
        
        // we store Glob_matMass_Dh = -omega^2 rho Ji(xi_k) \omega_k
        // where (xi_k, omega_k) is the quadrature formula used in the element
        // Glob_matMass_DhSigma = -i omega sigma Ji(xi_k) \omega_k  
        // Glob_matMass_Bh = Ji(xi_k) DF_i^{-1} \mu DF_i^{*-1} omega_k
        // all these expressions are computed without omega_k 
        // for constant coefficients (uniform rho, mu and affine element)
        var_helm.Glob_matMass_Dh(num_elem).Reallocate(N);
        var_helm.Glob_matMass_DhSigma(num_elem).Reallocate(N);
        var_helm.Glob_matMass_Bh(num_elem).Reallocate(N);
	if (var_helm.include_flow_term)
	  var_helm.Glob_matDamp_Mh(num_elem).Reallocate(N);
        
        // loop over quadrature points
        for (int j = 0; j < N; j++)
          {
            // evaluation of rho, sigma, mu
            rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, num_elem, j);
            sigma = var_helm.ref_sigma(ref).GetCoefficient(var_problem, num_elem, j);
            mu = var_helm.ref_mu(ref).GetCoefficient(var_problem, num_elem, j); 
	    if (var_helm.include_flow_term)
              {
                v0 = var_helm.ref_v0(ref).GetCoefficient(var_problem, num_elem, j);
                if (var_helm.include_double_gradM_flow)
                  {
                    beta = var_helm.ref_beta(ref).GetCoefficient(var_problem, num_elem, j);
                    // we subtract v0 v0^T to mu
                    Rank1Update(-beta, v0, mu);
                  }
              }
            
	    if (!var_helm.linearize_drude)
              var_helm.ref_drude(ref).ModifyCoefficientHelmholtz(var_problem.GetOmega(), rho, sigma, mu);
            
            if (var_problem.InsidePML(num_elem))
              var_helm.ModifyPMLCoefficient(rho, sigma, mu, invMu, v0, num_elem, i1, j);
            
            if (variable)
              {
                if (affine)
                  {
                    jacob = var_problem.Glob_jacobian(num_elem)(0)*Fb.WeightsND(j);
                    dfjm1 = var_problem.Glob_DFjm1(num_elem)(0);                
                    Mlt(1.0/var_problem.Glob_jacobian(num_elem)(0), dfjm1);
                  }
                else
                  {
                    jacob = var_problem.Glob_jacobian(num_elem)(j);
                    dfjm1 = var_problem.Glob_DFjm1(num_elem)(j);
                    Mlt(Fb.WeightsND(j)/jacob, dfjm1);
                  }
              }
            else
              {
                jacob = var_problem.Glob_jacobian(num_elem)(0);
                dfjm1 = var_problem.Glob_DFjm1(num_elem)(0);
                Mlt(1.0/jacob, dfjm1);
              }       
            
            MltTrans(mu, dfjm1, mu_dfjm1); 
            Mlt(dfjm1, mu_dfjm1, mat_stiff); 
            Mlt(jacob, mat_stiff);
            
            var_helm.Glob_matMass_Dh(num_elem)(j) = momega2*rho*jacob;
            var_helm.Glob_matMass_DhSigma(num_elem)(j) = m_iomega*sigma*jacob;
            var_helm.Glob_matMass_Bh(num_elem)(j) = mat_stiff;
	    
	    if (var_helm.include_flow_term)
	      {
		Mlt(dfjm1, v0, vec_u);
		vec_u *= jacob*m_iomega;
		var_helm.Glob_matDamp_Mh(num_elem)(j) = vec_u;
	      }
          }
      }
    else
      {
        // first order formulation
        // solved equations : -i omega rho u     + sigma u - div v = 0
        //                    -i omega mu^{-1} v - grad u = 0
        
        // in PML, there is no physical damping sigma,
        // but artificial dampings sigma_x, sigma_y, sigma_z
        // equations in 2-D PML :
        // -i omega rho u + (sigma_x + sigma_y)/2 u + (sigma_x - sigma_y)/2 u* - div v = f
        // -i omega rho u* + (sigma_x + sigma_y)/2 u* + (sigma_x - sigma_y)/2 u
        //                                  - (dvx/dx - dvy/dy) = 0
        // -i omega mu^{-1} v + sigma v - grad(u) = 0
        
        // in 2-D, we store :
        // Glob_matMass_Dh = -i omega rho Ji(xi_k) omega_k 
        // Glob_matMass_DhSigma = sigma Ji(xi_k) omega_k or (sigma_x + sigma_y)/2 
        //                           Ji(xi_k) omega_k in PML
        // Glob_matMass_DhSigmaDiff = (sigma_x - sigma_y)/2 Ji(xi_k) omega_k in PML
        
        // in 3-D, we store in PML :
        // Glob_matMass_DhSigma = (sigma_x + sigma_y + sigma_z)/3 Ji(xi_k) omega_k in PML
        // Glob_matMass_DhSigmaDiff = (sigma_x - sigma_y)/2 Ji(xi_k) omega_k in PML
        // Glob_matMass_DhSigmaDiff2 = (sigma_x - sigma_z)/2 Ji(xi_k) omega_k in PML
        // Glob_matMass_DhSigmaDiff3 = (sigma_y - sigma_z)/2 Ji(xi_k) omega_k in PML
        
        // These expression are computed without omega_k for constant coefficients
        var_helm.Glob_matMass_Dh(num_elem).Reallocate(N);
        var_helm.Glob_matMass_DhSigma(num_elem).Reallocate(N);
        if (var_problem.InsidePML(num_elem))
          {
            var_helm.Glob_matMass_DhSigmaDiff(i1).Reallocate(N);
            if (Dimension::dim_N == 3)
              {
                var_helm.Glob_matMass_DhSigmaDiff2(i1).Reallocate(N);
                var_helm.Glob_matMass_DhSigmaDiff3(i1).Reallocate(N);
              }
          }
        
        for (int k = 0; k < N; k++)
          {
            if (variable)
              {
                if (affine)
                  jacob = var_problem.Glob_jacobian(num_elem)(0)*Fb.WeightsND(k);
                else
                  jacob = var_problem.Glob_jacobian(num_elem)(k);
              }
            else
              jacob = var_problem.Glob_jacobian(num_elem)(0);
            
            rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, num_elem, k);
            sigma = var_helm.ref_sigma(ref).GetCoefficient(var_problem, num_elem, k);
            var_helm.Glob_matMass_Dh(num_elem)(k) = m_iomega*jacob*rho;
            var_helm.Glob_matMass_DhSigma(num_elem)(k) = jacob*sigma;
          }

        if (var_problem.InsidePML(num_elem))
          {
            T sigma_x, sigma_y, sigma_z, rho_tilde;
            for (int k = 0; k < N; k++)
              {
                rho = var_helm.ref_rho(ref).GetCoefficient(var_problem, num_elem, k);
                sigma = var_helm.ref_sigma(ref).GetCoefficient(var_problem, num_elem, k);
                rho_tilde = rho + sigma/m_iomega;
                if (variable)
                  {
                    if (affine)
                      jacob = var_problem.Glob_jacobian(num_elem)(0)*Fb.WeightsND(k);
                    else
                      jacob = var_problem.Glob_jacobian(num_elem)(k);
                  }
                else
                  jacob = var_problem.Glob_jacobian(num_elem)(0);
                
                sigma_x = var_boundary.GetTauPML(i1, k, 0);
                sigma_y = var_boundary.GetTauPML(i1, k, 1);
                if (Dimension::dim_N == 2)
                  {
                    var_helm.Glob_matMass_DhSigma(num_elem)(k)
                      += Real_wp(0.5)*rho_tilde*(sigma_x + sigma_y)*jacob;
                    
                    var_helm.Glob_matMass_DhSigmaDiff(i1)(k)
		      = Real_wp(0.5)*rho_tilde*(sigma_x - sigma_y)*jacob;		    
                  }
                else
                  {
                    sigma_z = var_boundary.GetTauPML(i1, k, 2);
                    var_helm.Glob_matMass_DhSigma(num_elem)(k)
                      += rho_tilde*(sigma_x + sigma_y + sigma_z)*jacob/Real_wp(3);
		    
                    var_helm.Glob_matMass_DhSigmaDiff(i1)(k) = rho_tilde*(sigma_x - sigma_y)*jacob/Real_wp(3);
                    var_helm.Glob_matMass_DhSigmaDiff2(i1)(k) = rho_tilde*(sigma_x - sigma_z)*jacob/Real_wp(3);
                    var_helm.Glob_matMass_DhSigmaDiff3(i1)(k) = rho_tilde*(sigma_y - sigma_z)*jacob/Real_wp(3);
                  }
              }
          }
      }
  }
  

  /*********************
   * VarHelmholtz_Base *
   *********************/
  
  
  //! computation of elementary matrix for Helmholtz equation
  /*!
    \param[in] iquad element number
    \param[out] num_dof degrees of freedom numbers
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] vars considered problem
    \param[in] Fb finite element associated with element iquad
   */
  template<class Complexe, class Dimension> template<class T>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::ComputeElementaryMatrixHelm(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				const GlobalGenericMatrix<T>& nat_mat,
				const VarProblem<Dimension>& vars,
				const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
				const ElementReference<Dimension, 1>& Fb)
  {
    int nb_dof_elt = Fb.GetNbDof();
    int nb_dof_all = nb_dof_elt;
    
    //const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
        
    // number of integration points
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    int ref_domain = vars.mesh.Element(iquad).GetReference();
    
    // dof numbers
    num_dof.Reallocate(nb_dof_all);
    IVect Nodle = vars.GetDofNumberOnElement(iquad);
    for (int j = 0; j < nb_dof_elt; j++)
      num_dof(j) = Nodle(j);
    
    int offset_uloc = nb_dof_elt;
    int offset_uloc2 = nb_dof_elt;
    int offset_vloc = nb_dof_elt;
    int nb_unknowns = 1;
    int offset_Ploc = 0, offset_Qloc = 0;
    int nb_dof_vec = Dimension::dim_N*nb_points_quad;    
    
    if (vars.FirstOrderFormulation())
      {
	if (this->ref_drude(ref_domain).IsEnabled())
          {
            if (this->ref_drude(ref_domain).IsModeTM())
	      {
		if (abs(this->ref_drude(ref_domain).eps_inf - this->ref_rho(ref_domain).GetConstant()) > epsilon_machine)
		  {
		    cout << "Rho must be equal to eps_infini" << endl;
                    cout << "Rho = " << this->ref_rho(ref_domain).GetConstant() << ", eps_infini = " << this->ref_drude(ref_domain).eps_inf << endl;
		    abort();
		  }
	      }
	    else
	      {
		if (abs(this->ref_drude(ref_domain).eps_inf - this->ref_invMu(ref_domain).GetConstant()(0, 0)) > epsilon_machine)
		  {
		    cout << "Mu must be equal to 1.0/eps_infini" << endl;
                    cout << "Mu^-1 = " << this->ref_invMu(ref_domain).GetConstant() << ", eps_infini = " << this->ref_drude(ref_domain).eps_inf << endl;
                    abort();
		  }
	      }
	  }
	
	if (vars.GetLeafStaticCondensation())
	  {
	  }
	else
	  {
            int offset_u2 = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
            int offset_v = mesh_num.GetNbDof() + (Dimension::dim_N-1)*mesh_num.GetNbDofPML();
            if (this->ref_drude(ref_domain).IsEnabled())
              {
                int nPole = this->ref_drude(ref_domain).gamma.GetM();
                if (vars.InsidePML(iquad))
                  {
                    cout << "Case not implemented" << endl;
                    abort();
                  }
		
                int offset_Q = vars.GetNbDof() - this->nb_dof_drude_vec;
                int offset_P = offset_Q - this->nb_dof_drude_vec;
                offset_P += this->OffsetDofDrudeV(iquad);
                offset_Q += this->OffsetDofDrudeV(iquad);
                                
		if (this->ref_drude(ref_domain).IsModeTE())
		  {
		    nb_dof_all += 2*nb_dof_vec*nPole;
		    offset_vloc += 2*nb_dof_vec*nPole;

		    offset_Ploc = nb_dof_elt;
		    offset_Qloc = offset_Ploc + nb_dof_vec*nPole;
		    
		    num_dof.Resize(nb_dof_all);
		    for (int j = 0; j < nPole*nb_dof_vec; j++)
		      {
			num_dof(offset_Ploc + j) = offset_P + j;
			num_dof(offset_Qloc + j) = offset_Q + j;
		      }
		  }
		else
		  {
		    nb_dof_all += 2*nb_points_quad*nPole;
		    offset_vloc += 2*nb_points_quad*nPole;
		    
		    offset_Ploc = nb_dof_elt;
		    offset_Qloc = nb_dof_elt + nb_points_quad*nPole;
		    num_dof.Resize(nb_dof_all);
		    for (int j = 0; j < nPole*nb_points_quad; j++)
		      {
			num_dof(offset_Ploc + j) = offset_P + j;
			num_dof(offset_Qloc + j) = offset_Q + j;
                      }
		  }
              }
	    
	    nb_dof_all += Dimension::dim_N*nb_points_quad;
	    if (vars.InsidePML(iquad))
	      {
		nb_dof_all += nb_dof_elt*(Dimension::dim_N-1);
		nb_unknowns = Dimension::dim_N;
		offset_uloc2 += nb_dof_elt;
		offset_vloc += nb_dof_elt*(Dimension::dim_N-1);
		num_dof.Resize(nb_dof_all);
		for (int j = 0; j < nb_dof_elt; j++)
		  num_dof(offset_uloc + j) = mesh_num.GetNbDof() + mesh_num.GetDofPML(Nodle(j));
		
		if (Dimension::dim_N == 3)
		  {
		    for (int j = 0; j < nb_dof_elt; j++)
		      num_dof(offset_uloc2 + j) = offset_u2 + mesh_num.GetDofPML(Nodle(j));
		  }
	      }
	    else
	      num_dof.Resize(nb_dof_all);
            
	    offset_v += vars.GetOffsetDofV(iquad); 
	    for (int j = 0; j < nb_points_quad; j++)
	      for (int k = 0; k < Dimension::dim_N; k++)
		num_dof(offset_vloc + j*Dimension::dim_N + k) = offset_v + Dimension::dim_N*j + k;
	  }
      }
    else
      {
        if (this->linearize_drude && this->ref_drude(ref_domain).IsEnabled())
          if (this->ref_drude(ref_domain).IsModeTM() && !var_problem.GetLeafStaticCondensation())
            {
              int nPole = this->ref_drude(ref_domain).gamma.GetM();
              nb_dof_all += nPole*nb_points_quad;
              num_dof.Resize(nb_dof_all);
              int offset_P = vars.GetNbDof() - this->nb_dof_drude_vec;
              offset_P += this->OffsetDofDrudeV(iquad);
              offset_Ploc = nb_dof_elt;
              for (int j = 0; j < nb_points_quad*nPole; j++)
                num_dof(offset_Ploc + j) = offset_P + j;
            }
      }

    mat_interac.Clear();
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();

    bool variable = vars.UseNumericalIntegration(iquad);
    bool affine = vars.mesh.IsElementAffine(iquad);
    bool stiff = false;
    if (nat_mat.GetCoefStiffness() != Real_wp(0))
      stiff = true;

    T one; SetComplexOne(one);
    Complexe m_iomega;
    vars.GetMiomega(m_iomega);

    if (vars.FirstOrderFormulation() && vars.InsidePML(iquad))
      {
	// no damping with pml for first order formulation (not implemented)
	bool presence_damping = false;
	if (!this->ref_sigma(ref_domain).IsZero())
	  presence_damping = true;
	
	if (presence_damping)
	  {
	    cout << "Damping with PML and first order formulation not implemented" << endl;
	    abort();
	  }
      }
    
    if ((vars.FirstOrderFormulation()) && (vars.GetLeafStaticCondensation()))
      {
        // loop over dofs 
        Real_wp jacobian;
                
	if (this->include_flow_term)
	  {
	    cout << "not implemented " << endl;
	    abort();
	  }

	if (vars.InsidePML(iquad))
	  {
 	    // only Neumann condition is accepted
	    bool presence_damping = false;
	    for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	      {
		int num_face = vars.mesh.Element(iquad).numBoundary(num_loc);
		int ref_line = vars.mesh.Boundary(num_face).GetReference();
		if (vars.mesh.GetBoundaryCondition(ref_line) == BoundaryConditionEnum::LINE_ABSORBING)
		  presence_damping = true;
	      }

	    if (presence_damping)
	      {
		cout << "Only Neumann or Dirichlet condition is accepted with static condensation and PML" << endl;
		abort();
	      }
	  }
	
    	// avec condensation
	TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> dfjm1;
	TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> mu_dfjm1;

        T coef = square(nat_mat.GetCoefStiffness()) / (m_iomega*nat_mat.GetCoefMass());
        T coef_drude, coef_drude_vec; SetComplexOne(coef_drude); SetComplexOne(coef_drude_vec);
        if (this->ref_drude(ref_domain).IsEnabled())
          {
            Real_wp eps_inf = this->ref_drude(ref_domain).eps_inf;
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            for (int k = 0; k < nPole; k++)
              {
                Real_wp gamma = this->ref_drude(ref_domain).gamma(k);
                Real_wp coef_omega2 = this->ref_drude(ref_domain).omega_02(k);
                Real_wp coef_omegap = this->ref_drude(ref_domain).eps_omega_p2(k) / eps_inf;
                Real_wp coef_sigma = this->ref_drude(ref_domain).eps_sigma(k) / eps_inf;
                T s = nat_mat.GetCoefStiffness();
                T m = nat_mat.GetCoefMass()*m_iomega;
                coef_drude += (coef_omegap*s*s + coef_sigma*m*s) / (coef_omega2*s*s + m*m + gamma*s*m);
              }
          }
        
        if (this->ref_drude(ref_domain).IsEnabled() && this->ref_drude(ref_domain).IsModeTE())
          {
            coef_drude_vec = Real_wp(1) / coef_drude;
            SetComplexOne(coef_drude);
          }
        
	if (variable)
	  {
	    Vector<T> Amass(nb_points_quad);
            Vector<TinyMatrix<T, General,
                              Dimension::dim_N, Dimension::dim_N> > Cgrad(nb_points_quad);
            
            Vector<TinyVector<T, Dimension::dim_N> > Evec, Dvec;
	    
	    if (vars.InsidePML(iquad))
	      {
		TinyVector<bool, 4> null_term(false, false, true, true);
		if (!Fb.LumpedMassMatrix())
		  {
		    cout << "Case not implemented" << endl;
		    abort();
		  }

		int i1 = iquad - vars.mesh.GetNbElt() + var_boundary.GetNbEltPML();

		TinyVector<Complexe, Dimension::dim_N> tau;
		T invCoef = one / (m_iomega * nat_mat.GetCoefMass());
		if (Dimension::dim_N == 3)
		  invCoef *= invCoef;
		
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    if (affine)
		      {
			jacobian = var_problem.Glob_jacobian(iquad)(0);
			dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
		      }
		    else
		      {
			jacobian = var_problem.Glob_jacobian(iquad)(k) / Fb.WeightsND(k);
			dfjm1 = var_problem.Glob_DFjm1(iquad)(k);
		      }

		    Mlt(Real_wp(1) / jacobian, dfjm1);

		    tau = var_boundary.GetTauPML(i1, k);
		    T factor_x = m_iomega*nat_mat.GetCoefMass() + tau(0)*nat_mat.GetCoefDamping();
		    T factor_y = m_iomega*nat_mat.GetCoefMass() + tau(1)*nat_mat.GetCoefDamping();
		    T factor_z; SetComplexOne(factor_z);
		    if (Dimension::dim_N == 3)
		      factor_z = m_iomega*nat_mat.GetCoefMass() + tau(2)*nat_mat.GetCoefDamping();
		    
		    Complexe rho = this->ref_rho(ref_domain).GetCoefficient(vars, iquad, k);
		    Complexe sigma = this->ref_sigma(ref_domain).GetCoefficient(vars, iquad, k);
		    Complexe rho_tilde = rho + sigma/m_iomega;

		    jacobian *= Fb.WeightsND(k);
		    Amass(k) = factor_x*factor_y*factor_z*rho_tilde*jacobian*invCoef;
		    
		    MltTrans(this->ref_mu(ref_domain).GetCoefficient(vars, iquad, k), dfjm1, mu_dfjm1);
		    mu_dfjm1(0, 0) *= invCoef*factor_y*factor_z/factor_x;
		    mu_dfjm1(1, 1) *= invCoef*factor_x*factor_z/factor_y;
		    if (Dimension::dim_N == 3)
		      mu_dfjm1(2, 2) *= invCoef*factor_x*factor_y/factor_z;
		    
		    Mlt(dfjm1, mu_dfjm1, Cgrad(k));

		    Cgrad(k) *= jacobian*square(nat_mat.GetCoefStiffness());
		  }
		
		Fb.AddVariableElemMatrix(0, 0, Amass, Cgrad, Dvec, Evec, null_term, mat_interac);
		
		//mat_interac.Write("mat_test.dat"); DISP(Fb.GetNbDof());
		//cout << "waiting" << endl; int test_input; cin >> test_input;
	      }
	    else
	      {
		// computing elementary matrix with methods of the finite element class
		bool no_grad = true, no_damp = true;
                
		for (int i = 0; i < nb_points_quad; i++)
		  {
		    Amass(i) = (this->Glob_matMass_Dh(iquad)(i)*nat_mat.GetCoefMass()*coef_drude
				+ this->Glob_matMass_DhSigma(iquad)(i)*nat_mat.GetCoefDamping());
		    
		    if (affine)
		      {
			jacobian = var_problem.Glob_jacobian(iquad)(0);
			dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
		      }
		    else
		      {
			jacobian = var_problem.Glob_jacobian(iquad)(i) / Fb.WeightsND(i);
			dfjm1 = var_problem.Glob_DFjm1(iquad)(i);
		      }
		    
		    Mlt(Real_wp(1) / jacobian, dfjm1);
		    jacobian *= Fb.WeightsND(i);
		    MltTrans(this->ref_mu(ref_domain).GetCoefficient(vars, iquad, i), dfjm1, mu_dfjm1);
		    Mlt(dfjm1, mu_dfjm1, Cgrad(i));
		    
		    Cgrad(i) *= jacobian*coef*coef_drude_vec;
		    if (!Cgrad(i).IsZero())
		      no_grad = false;
		  }
		
		TinyVector<bool, 4> null_term(false, no_grad, no_damp, no_damp);
		Fb.AddVariableElemMatrix(0, 0, Amass, Cgrad, Dvec, Evec, null_term, mat_interac);
	      }
	  }
	else
	  {
	    TinyVector<bool, 4> null_term;
	    null_term.Fill(false);
	    TinyVector<T, Dimension::dim_N> Evec, Dvec;
	    TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> Cgrad;
            
	    jacobian = var_problem.Glob_jacobian(iquad)(0);
	    dfjm1 = var_problem.Glob_DFjm1(iquad)(0);                
	    Mlt(1.0/var_problem.Glob_jacobian(iquad)(0), dfjm1);
	    
	    MltTrans(this->ref_mu(ref_domain).GetConstant(), dfjm1, mu_dfjm1);
	    Mlt(dfjm1, mu_dfjm1, Cgrad);
		
	    Cgrad *= jacobian*square(nat_mat.GetCoefStiffness()) / (m_iomega*nat_mat.GetCoefMass()) * coef_drude_vec;
	    T Amass = (this->Glob_matMass_Dh(iquad)(0)*nat_mat.GetCoefMass()*coef_drude
		       + this->Glob_matMass_DhSigma(iquad)(0)*nat_mat.GetCoefDamping());
	    
	    null_term(2) = true;
	    null_term(3) = true;
	    
	    Fb.AddConstantElemMatrix(0, 0, Amass, Cgrad, Dvec, Evec, null_term, mat_interac);		
	  }	
      }
    else if (vars.FirstOrderFormulation())
      {
        // loop over dofs 
        T poids;
        VectReal_wp val_phi(nb_points_quad);
        
	if (this->include_flow_term)
	  {
	    cout << "not implemented " << endl;
	    abort();
	  }
	
        Vector<T> feval(nb_points_quad), contrib(nb_dof_elt);
        Vector<TinyVector<Real_wp, Dimension::dim_N> > grad_phi(nb_points_quad); 
        TinyVector<Real_wp, Dimension::dim_N> vec_u;
        val_phi.Zero(); grad_phi.Zero();
 	feval.Zero(); contrib.Zero();
        
        if (Fb.LumpedMassMatrix() && !stiff && !this->ref_drude(ref_domain).IsEnabled())
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                if (variable)
                  poids = this->Glob_matMass_Dh(iquad)(k)*nat_mat.GetCoefMass() +
                    this->Glob_matMass_DhSigma(iquad)(k)*nat_mat.GetCoefDamping();
                else
                  poids = Fb.WeightsND(k)*(this->Glob_matMass_Dh(iquad)(0)*nat_mat.GetCoefMass() +
                                           this->Glob_matMass_DhSigma(iquad)(0)
                                           *nat_mat.GetCoefDamping());
                
                mat_interac.AddInteraction(k, k, poids);
                if (vars.InsidePML(iquad))
                  {
                    mat_interac.AddInteraction(offset_uloc+k, offset_uloc+k, poids);
                    if (Dimension::dim_N == 3)
                      mat_interac.AddInteraction(offset_uloc2+k, offset_uloc2+k, poids);
                  }
              }

            HelmholtzElementaryMatrixClass<Dimension>
              ::ComputeMatrixVdiag(vars, var_boundary, *this, iquad, nb_points_quad, 
                                   affine, variable, offset_uloc, offset_uloc2,
                                   offset_vloc, nat_mat, mat_interac);

            return;
          }

        for (int i = 0; i < nb_dof_elt; i++)
          {
            // computation of phi_i(xi_k) and grad phi_i(xi_k) (in reference element)
	    if (stiff)
	      Fb.GetGradientSinglePhiQuadrature(i, val_phi, grad_phi);
	    else
	      Fb.GetValueSinglePhiQuadrature(i, val_phi);

	    // computation of gradient in real element
            if (stiff)
              {
                if (affine)
                  {
                    typename Dimension::MatrixN_N dfjm1 = vars.Glob_DFjm1(iquad)(0);
                    for (int k = 0; k < nb_points_quad; k++)
                      {
			vec_u = grad_phi(k);
                        MltTrans(dfjm1, vec_u, grad_phi(k));
                      }
                  }
                else
                  for (int k = 0; k < nb_points_quad; k++)
                    {
		      vec_u = grad_phi(k);
                      MltTrans(vars.Glob_DFjm1(iquad)(k), vec_u, grad_phi(k));
                    }
              }
	    
            // partie avec v
	    HelmholtzElementaryMatrixClass<Dimension>::
	      ComputeMatrixV(vars, var_boundary, *this, iquad, i, nb_points_quad, affine,
			     offset_uloc, offset_uloc2, offset_vloc, nat_mat,
			     stiff, variable, val_phi, grad_phi, feval, contrib,
			     mat_interac);

            // partie avec u
            for (int k = 0; k < nb_points_quad; k++)
              {
                if (variable)
                  poids = this->Glob_matMass_Dh(iquad)(k)*nat_mat.GetCoefMass() +
                    this->Glob_matMass_DhSigma(iquad)(k)*nat_mat.GetCoefDamping();
                else
                  poids = Fb.WeightsND(k)*(this->Glob_matMass_Dh(iquad)(0)*nat_mat.GetCoefMass() +
                                           this->Glob_matMass_DhSigma(iquad)(0)
                                           *nat_mat.GetCoefDamping());
                
                feval(k) = poids*val_phi(k);                
              }
                        
            contrib.Fill(0);
            Fb.ApplyCh(feval, contrib);
            for (int j = 0; j < nb_dof_elt; j++)
              if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                mat_interac.AddInteraction(i, j, contrib(j));
            
            if (vars.InsidePML(iquad))
              {
                for (int j = 0; j < nb_dof_elt; j++)
                  if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                    mat_interac.AddInteraction(offset_uloc+i, offset_uloc+j, contrib(j));
                
                if (Dimension::dim_N == 3)
                  for (int j = 0; j < nb_dof_elt; j++)
                    if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                      mat_interac.AddInteraction(offset_uloc2+i, offset_uloc2+j, contrib(j));
              }
          }
        
        if (this->ref_drude(ref_domain).IsEnabled())
          {
            if (!Fb.DiagonalMassMatrix() && (this->ref_drude(ref_domain).IsModeTM()))
              {
                cout << "Drude not implemented in this case" << endl;
                abort();
              }

            bool sym = vars.GetSymmetrizationUse();
            Real_wp coef_sym(1);
            if (sym)
              coef_sym = Real_wp(-1);
            
            // material with Drude's model
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            for (int kp = 0; kp < nPole; kp++)
              {
                Real_wp gamma = this->ref_drude(ref_domain).gamma(kp);
                Real_wp coef_omega2 = this->ref_drude(ref_domain).omega_02(kp);
                Real_wp coef_eps_inf = this->ref_drude(ref_domain).eps_omega_p2(kp);
                Real_wp coef_sig = this->ref_drude(ref_domain).eps_sigma(kp);
                
                Real_wp coef_Q = coef_sym / coef_eps_inf;
                Real_wp coef_P = coef_omega2 / coef_eps_inf;                
                if (!this->use_symm_drude)
                  {
                    // pas de symmetrisation pour traiter omega_0 = 0
                    coef_Q = coef_sym;
                    coef_P = Real_wp(1);
                  }
                
                if (this->ref_drude(ref_domain).IsModeTE())
                  {
                    coef_P *= coef_sym;
                    coef_Q *= coef_sym;
                  }
                
                for (int k = 0; k < nb_points_quad; k++)
                  {
                    Real_wp poids;
                    if (affine)
                      poids = Fb.WeightsND(k) * var_problem.Glob_jacobian(iquad)(0);
                    else
                      poids = var_problem.Glob_jacobian(iquad)(k);
                    
                    if (this->ref_drude(ref_domain).IsModeTM())
                      {
                        int num_dofP = offset_Ploc + k*nPole + kp;
                        int num_dofQ = offset_Qloc + k*nPole + kp; int num_dofE = k;
                        // line dP/dt - Q = 0
                        // this line is multiplied by omega_0^2 / (eps_inf omega_p^2)
                        // to obtain symmetry
                        mat_interac.SetEntry(num_dofP, num_dofP, coef_P*m_iomega*nat_mat.GetCoefMass()*poids);
                        mat_interac.SetEntry(num_dofP, num_dofQ, -coef_P*nat_mat.GetCoefStiffness()*poids);

                        // line dQ/dt + gamma Q + omega_0^2 p - eps_inf omega_p^2 u - eps_inf sigma du/dt = 0
                        // this line is divided by -eps_inf omega_p^2 to obtain symmetry
                        mat_interac.SetEntry(num_dofQ, num_dofQ, coef_Q * (m_iomega*nat_mat.GetCoefMass()
                                                                           + gamma*nat_mat.GetCoefStiffness())*poids);
		    
                        mat_interac.SetEntry(num_dofQ, num_dofP, coef_omega2*coef_Q*nat_mat.GetCoefStiffness()*poids);
                        mat_interac.SetEntry(num_dofQ, num_dofE, -coef_Q*poids*(coef_eps_inf*nat_mat.GetCoefStiffness() 
                                                                                + coef_sig*m_iomega*nat_mat.GetCoefMass()));
                        
                        // term +Q in equation for u
                        mat_interac.SetEntry(num_dofE, num_dofQ, nat_mat.GetCoefStiffness()*poids);
                      }
                    else
                      {
                        int num_dofPx = offset_Ploc + 2*k*nPole + kp, num_dofPy = offset_Ploc + (2*k+1)*nPole + kp;
                        int num_dofQx = offset_Qloc + 2*k*nPole + kp, num_dofQy = offset_Qloc + (2*k+1)*nPole + kp;
                        int num_dofHx = offset_vloc + 2*k, num_dofHy = offset_vloc + (2*k+1);
                        
                        // line dP/dt - Q = 0
                        // this line is multiplied by omega_0^2 / (eps_inf omega_p^2)
                        // to obtain symmetry
                        mat_interac.SetEntry(num_dofPx, num_dofPx, m_iomega*coef_P*nat_mat.GetCoefMass()*poids);
                        mat_interac.SetEntry(num_dofPy, num_dofPy, m_iomega*coef_P*nat_mat.GetCoefMass()*poids);
                        
                        mat_interac.SetEntry(num_dofPx, num_dofQx, -coef_P*nat_mat.GetCoefStiffness()*poids);
                        mat_interac.SetEntry(num_dofPy, num_dofQy, -coef_P*nat_mat.GetCoefStiffness()*poids);
                        
                        // line dQ/dt + gamma Q + omega_0^2 p - eps_inf omega_p^2 v - eps_inf sigma dv/dt = 0
                        // this line is divided by -eps_inf omega_p^2 to obtain symmetry
                        mat_interac.SetEntry(num_dofQx, num_dofQx, coef_Q*(m_iomega*nat_mat.GetCoefMass()
                                                                           + gamma*nat_mat.GetCoefStiffness())*poids);
                        
                        mat_interac.SetEntry(num_dofQy, num_dofQy, coef_Q*(m_iomega*nat_mat.GetCoefMass()
                                                                           + gamma*nat_mat.GetCoefStiffness())*poids);
                        
                        mat_interac.SetEntry(num_dofQx, num_dofPx, coef_omega2*coef_Q*nat_mat.GetCoefStiffness()*poids);
                        mat_interac.SetEntry(num_dofQy, num_dofPy, coef_omega2*coef_Q*nat_mat.GetCoefStiffness()*poids);

                        T val_mat = -coef_Q*poids*(coef_eps_inf*nat_mat.GetCoefStiffness()
                                                   + coef_sig*m_iomega*nat_mat.GetCoefMass());
                        
                        mat_interac.SetEntry(num_dofQx, num_dofHx, val_mat);                        
                        mat_interac.SetEntry(num_dofQy, num_dofHy, val_mat);
                        
                        // term +Q in the equation of v
                        mat_interac.SetEntry(num_dofHx, num_dofQx, coef_sym*nat_mat.GetCoefStiffness()*poids);
                        mat_interac.SetEntry(num_dofHy, num_dofQy, coef_sym*nat_mat.GetCoefStiffness()*poids);
                      }
                  }
              }
	  }
        
	if (vars.InsidePML(iquad))
	  {
	    //mat_interac.Write("mat_ref.dat"); DISP(Fb.GetNbDof());
	    //cout << "waiting" << endl; int test_input; cin >> test_input;
	  }
      }
    else
      {
        T coef_drude; SetComplexOne(coef_drude);
        if (vars.GetLeafStaticCondensation() && this->linearize_drude
            && this->ref_drude(ref_domain).IsEnabled() && this->ref_drude(ref_domain).IsModeTM())
          {
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            Real_wp eps_inf = this->ref_drude(ref_domain).eps_inf;
            T s = nat_mat.GetCoefStiffness();
            T m = nat_mat.GetCoefMass()*m_iomega*m_iomega;
            T d = nat_mat.GetCoefDamping()*m_iomega;
            
            for (int kp = 0; kp < nPole; kp++)
              {
                Real_wp gamma = this->ref_drude(ref_domain).gamma(kp);
                Real_wp coef_omega2 = this->ref_drude(ref_domain).omega_02(kp);
                Real_wp coef_omegap = this->ref_drude(ref_domain).eps_omega_p2(kp) / eps_inf;
                Real_wp coef_sig = this->ref_drude(ref_domain).eps_sigma(kp) / eps_inf;
                
                coef_drude += (coef_omegap*s + coef_sig*d) / (coef_omega2*s + m + gamma*d);
              }
          }

        if (variable)
          {
            Vector<T> Amass(nb_points_quad);
            Vector<TinyMatrix<T, General,
                              Dimension::dim_N, Dimension::dim_N> > Cgrad(nb_points_quad);
            
            Vector<TinyVector<T, Dimension::dim_N> > Evec, Dvec;
            if (this->include_flow_term)
              {
                Dvec.Reallocate(nb_points_quad);
                Evec.Reallocate(nb_points_quad);
              }
            
            // computing elementary matrix with methods of the finite element class
            bool no_grad = true, no_damp = true;
            for (int i = 0; i < nb_points_quad; i++)
              {
                Amass(i) = (this->Glob_matMass_Dh(iquad)(i)*nat_mat.GetCoefMass()
                            + this->Glob_matMass_DhSigma(iquad)(i)*nat_mat.GetCoefDamping());
                
                Amass(i) *= coef_drude;
                Cgrad(i) = this->Glob_matMass_Bh(iquad)(i);
                Cgrad(i) *= nat_mat.GetCoefStiffness();
                if (!Cgrad(i).IsZero())
                  no_grad = false;
                
                if (this->include_flow_term)
                  {
                    Evec(i) = this->Glob_matDamp_Mh(iquad)(i);
                    Evec(i) *= nat_mat.GetCoefDamping();
                    Dvec(i) = -Evec(i);
                    if (!Evec(i).IsZero())
                      no_damp = false;
                  }
              }
            
            TinyVector<bool, 4> null_term(false, no_grad, no_damp, no_damp);
            Fb.AddVariableElemMatrix(0, 0, Amass, Cgrad, Dvec, Evec, null_term, mat_interac);
          }
        else
          {
	    TinyVector<bool, 4> null_term;
	    null_term.Fill(false);
	    TinyVector<T, Dimension::dim_N> Evec, Dvec;
            TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> Cgrad;
	    Cgrad = this->Glob_matMass_Bh(iquad)(0);
            
            Cgrad *= nat_mat.GetCoefStiffness();
            T Amass = (this->Glob_matMass_Dh(iquad)(0)*nat_mat.GetCoefMass()
		       + this->Glob_matMass_DhSigma(iquad)(0)*nat_mat.GetCoefDamping());
            
            Amass *= coef_drude;            
	    if (this->include_flow_term)
	      {
		Evec = this->Glob_matDamp_Mh(iquad)(0)*nat_mat.GetCoefDamping();
		Dvec = -Evec;
	      }
	    else
	      {
		null_term(2) = true;
		null_term(3) = true;
	      }

	    if (Cgrad.IsZero())
	      null_term(1) = true;
	    
	    Fb.AddConstantElemMatrix(0, 0, Amass, Cgrad, Dvec, Evec, null_term, mat_interac);	    
          }
        
        if (this->linearize_drude && this->ref_drude(ref_domain).IsEnabled()
            && this->ref_drude(ref_domain).IsModeTM() && !var_problem.GetLeafStaticCondensation())
          {
            //DISP(iquad); DISP(num_dof); DISP(offset_Ploc);
            Complexe m_omega2 = m_iomega*m_iomega;
            Vector<Real_wp> poids(nb_points_quad);
            Vector<T> Amass(nb_points_quad);
            for (int i = 0; i < nb_points_quad; i++)
              {
                if (affine)
                  poids(i) = Fb.WeightsND(i) * var_problem.Glob_jacobian(iquad)(0);
                else
                  poids(i) = var_problem.Glob_jacobian(iquad)(i);                
              }

            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            Complexe rho0 = this->ref_rho(ref_domain).GetConstant();
            for (int k = 0; k < nPole; k++)
              {
                for (int i = 0; i < nb_points_quad; i++)
                  {
                    int dofP = offset_Ploc + i*nPole + k;
                    if (vars.InsidePML(iquad))
                      Amass(i) = this->Glob_matMass_Dh(iquad)(i)*nat_mat.GetCoefMass() / rho0;
                    else
                      Amass(i) = poids(i)*nat_mat.GetCoefMass()*m_omega2;
                    
                    // we add - omega^2 \sum P_k for equation in u
                    mat_interac.SetEntry(i, dofP, Amass(i));
                    
                    // we add equations for P_k
                    // -omega^2 P_k - i omega gamma_k P_k + omega_0^2 P_k = c_k u - i omega sigma_k u
                    Amass(i) = poids(i)*nat_mat.GetCoefMass()*m_omega2;
                    Amass(i) += poids(i)*nat_mat.GetCoefDamping()*m_iomega*this->ref_drude(ref_domain).gamma(k);
                    Amass(i) += poids(i)*nat_mat.GetCoefStiffness()*this->ref_drude(ref_domain).omega_02(k);
                    mat_interac.SetEntry(dofP, dofP, Amass(i));
                    
                    Amass(i) = -poids(i)*nat_mat.GetCoefStiffness()*this->ref_drude(ref_domain).eps_omega_p2(k);
                    Amass(i) -= poids(i)*nat_mat.GetCoefDamping()*m_iomega*this->ref_drude(ref_domain).eps_sigma(k);
                    mat_interac.SetEntry(dofP, i, Amass(i));
                  }
              }
          }
      }

    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad, nb_unknowns);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad, nb_unknowns);
  }
      
  
  //! sets values of attributes to default values
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::InitDefaultValues()
  {
    this->rho0 = 1.0; this->mu0 = 1.0; this->sigma0 = 0.0;
    include_flow_term = false;
    include_double_gradM_flow = false;
    this->use_symm_drude = false;
  }
  
  
  //! reading of a line of the data file
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("AddFlowTerm"))
      {
	if (parameters(0) == "YES")
	  {
	    include_flow_term = true;
	    var_problem.SetSymmetricElementaryMatrix(false);
	  }
	else if (parameters(0) == "GRAD")
	  {
	    include_flow_term = true;
	    include_double_gradM_flow = true;
	    var_problem.SetSymmetricElementaryMatrix(false);
	  }
	else
	  include_flow_term = false;
      }
    else if (description_field == "SymmetrizationDrude")
      {
        if (parameters(0) == "YES")
          this->use_symm_drude = true;
        else
          this->use_symm_drude = false;
      }
    else if (description_field == "LinearizationDrude")
      {
        if (parameters(0) == "YES")
          this->linearize_drude = true;
        else
          this->linearize_drude = false;
      }
    else if (!description_field.compare("AddSlot"))
      {
#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL
	//cout<<"adding slot"<<endl;
	if (parameters.GetM() <= 5)
	  {
	    cout << "In SetInputData of VarHelmholtz_Base" << endl;
	    cout << "AddSlot needs more parameters, for instance :" << endl;
	    cout << "AddSlot = xA yA xB yB N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	SlotModelParameters<Dimension2> slot;
	slot.ptA(0) = to_num<Real_wp>(parameters(0));
	slot.ptA(1) = to_num<Real_wp>(parameters(1));
	slot.ptB(0) = to_num<Real_wp>(parameters(2));
	slot.ptB(1) = to_num<Real_wp>(parameters(3));
	slot.epsilon = to_num<Real_wp>(parameters(4));
	slot.length = slot.ptA.Distance(slot.ptB);
	// straight slot
	// vec_u in the direction of the slot
	R2 vec_u = slot.ptA - slot.ptB;
	// finding the orthogonal vector
	R2 vec_v(vec_u(1), -vec_u(0));
	Mlt(1.0/Norm2(vec_v), vec_v);
	// left section of the slot
	slot.ptC(0)(0) = slot.ptA(0) + 0.5*slot.epsilon*vec_v(0);
	slot.ptC(0)(1) = slot.ptA(1) + 0.5*slot.epsilon*vec_v(1);
	slot.ptD(0)(0) = slot.ptA(0) - 0.5*slot.epsilon*vec_v(0);
	slot.ptD(0)(1) = slot.ptA(1) - 0.5*slot.epsilon*vec_v(1);
	// right section of the slot
	slot.ptC(1)(0) = slot.ptB(0) + 0.5*slot.epsilon*vec_v(0);
	slot.ptC(1)(1) = slot.ptB(1) + 0.5*slot.epsilon*vec_v(1);
	slot.ptD(1)(0) = slot.ptB(0) - 0.5*slot.epsilon*vec_v(0);
	slot.ptD(1)(1) = slot.ptB(1) - 0.5*slot.epsilon*vec_v(1);
	
	// slot.order = to_num<int>(parameters(6));
	slot.nb_points = to_num<int>(parameters(5));
	
	DISP(slot);
	list_slots.PushBack(slot);
      }
    else if (!description_field.compare("ModelSlot"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHelmholtz_Base" << endl;
	    cout << "ModelSlot needs more parameters, for instance :" << endl;
	    cout << "ModelSlot = MESH_1D" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if (!parameters(0).compare("DTN"))
	  SlotModelParameters<Dimension2>::type_model
            = SlotModelParameters<Dimension2>::MODEL_DTN;
	else if (!parameters(0).compare("MESH_1D"))
	  SlotModelParameters<Dimension2>::type_model
            = SlotModelParameters<Dimension2>::MODEL_MESH1D;
#endif
      }
    
  }


  //! copies input parameters of another similar problem
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::CopyInputData(const VarHelmholtz_Base<Complexe, Dimension>& var)
  {
    this->rho0 = var.rho0;
    this->mu0 = var.mu0;
    this->ref_drude = var.ref_drude;
    
    ref_rho = var.ref_rho;
    ref_invRho = var.ref_invRho;
    ref_mu = var.ref_mu;
    ref_invMu = var.ref_invMu;
    ref_sigma = var.ref_sigma;
    ref_v0 = var.ref_v0;
    ref_beta = var.ref_beta;    

    include_flow_term = var.include_flow_term;
    include_double_gradM_flow = var.include_double_gradM_flow;
  }
  
  
  //! initialization of physical indexes rho and mu
  /*!
    \param[in] n number of physical domains
    The indexes are initialized to rho = 1, mu = 1, sigma = 0
   */
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::InitIndices(int n)
  {
    if (n <= 1)
      {
        cout<<"Number of physical domains should be greater than 2"<<endl;
	abort();
      }
    
    Real_wp rho0, mu0;
    if (PhysicalConstant::adimensionalization != PhysicalConstant::ADIM_NO)
      {
	// non-physical units, rho_0 and mu_0 set to 1
	rho0 = 1.0;
	mu0 = 1.0;
      }
    else
      {
	// use of physical units
	rho0 = PhysicalConstant::epsilon0_permittivity;
	mu0 = 1.0/PhysicalConstant::mu0_permeability;
      }

    ref_rho.Reallocate(n+1);
    ref_invRho.Reallocate(n+1);
    ref_mu.Reallocate(n+1);
    ref_invMu.Reallocate(n+1);
    ref_sigma.Reallocate(n+1);
    ref_v0.Reallocate(n+1);
    ref_beta.Reallocate(n+1);
    this->ref_drude.Reallocate(n+1);
    
    for (int i = 0; i <= n; i++)
      {
	ref_rho(i).SetIdentity(); 
	ref_rho(i).Mlt(rho0);
	ref_invRho(i).SetIdentity();
	ref_invRho(i).Mlt(1.0/rho0);
	ref_mu(i).SetIdentity();
	ref_mu(i).Mlt(mu0);
        ref_invMu(i).SetIdentity();
	ref_invMu(i).Mlt(1.0/mu0);
	ref_sigma(i).Zero();
	ref_v0(i).Zero();
	ref_beta(i).SetIdentity(); 
      }
  }
  
  
  //! returns the number of physical media
  template<class Complexe, class Dimension>
  int VarHelmholtz_Base<Complexe, Dimension>::GetNbPhysicalIndices() const
  {
    return ref_rho.GetM();
  }

  
  //! modification of physical indexes according to parameters of data file
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

    ref_rho(i).SetInputData(nb, parameters, parameters(0));
    ref_mu(i).SetInputData(nb, parameters, parameters(0));
    ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    if (include_flow_term)
      ref_v0(i).SetInputData(nb, parameters, parameters(0));
    
    if (include_double_gradM_flow)
      ref_beta(i).SetInputData(nb, parameters, parameters(0));
  }
  
  
  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

    if (name_media == "rho")
      ref_rho(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "mu")
      ref_mu(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "sigma")
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "M")
      {
        if (include_flow_term)
          ref_v0(i).SetInputData(nb, parameters, parameters(0));
      }
    else if (name_media == "beta")
      {
        if (include_double_gradM_flow)
          ref_beta(i).SetInputData(nb, parameters, parameters(0));
      }
    else if (name_media == "Drude")
      {
	this->ref_drude(i).SetInputData(nb, parameters, parameters(0));
	this->ref_drude(i).Adimensionalize(var_problem.GetWaveLengthAdim());

        TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;
        mu.SetIdentity();
        Complexe rho(1), sigma(0);
        
        if (this->ref_drude(i).IsModeTM())
          rho = this->ref_drude(i).eps_inf;
        else
          mu.SetDiagonal(1.0/this->ref_drude(i).eps_inf);
        
        this->ref_rho(i).SetConstant(rho);
        this->ref_mu(i).SetConstant(mu);
        this->ref_sigma(i).SetConstant(sigma);
      }
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }

  
  //! returns the name associated with the physical index num
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::SetPhysicalIndexAtInfinity(const Vector<bool>&)
  {
    int ref = var_problem.GetReferenceInfinity();
    if (ref >= 0)
      {
	this->rho0 = ref_rho(ref).GetConstant();
	this->mu0 = ref_mu(ref).GetConstant()(0, 0);
	R_N kwave = var_problem.GetWaveVector();
	kwave *= realpart(sqrt(this->rho0/this->mu0));
	var_problem.SetWaveVector(kwave);
      }
  }

  
  //! returns the name associated with the physical index num
  template<class Complexe, class Dimension>
  string VarHelmholtz_Base<Complexe, Dimension>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("rho");
      case 1: return string("mu");
      case 2: return string("sigma");
      case 3: return string("M");
      case 4: return string("beta");
      }
    
    return string();
  }  
  
  
  //! fills varying indices present in the current problem
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complexe>* >& rho_complex, IVect& num_ref,
		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    int nb = 0;
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        nb += ref_rho(i).GetNbVaryingMedia();
        nb += ref_mu(i).GetNbVaryingMedia();
        nb += ref_sigma(i).GetNbVaryingMedia();
	if (include_flow_term)
	  nb += ref_v0(i).GetNbVaryingMedia();

	if (include_double_gradM_flow)
	  nb += ref_beta(i).GetNbVaryingMedia();
      }
    
    rho_complex.Reallocate(nb);
    num_ref.Reallocate(nb);
    num_index.Reallocate(nb);
    num_component.Reallocate(nb);
    compute_grad.Reallocate(nb);
    compute_hess.Reallocate(nb);
    compute_grad.Fill(false);
    compute_hess.Fill(false);
    nb = 0;
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        int nb0 = nb;
        ref_rho(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        ref_mu(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 1;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        ref_sigma(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 2;
            num_ref(j) = i;
          }

        if (include_flow_term)
	  {
	    nb0 = nb;
	    ref_v0(i).GetVaryingMedia(nb, rho_complex, num_component);
	    for (int j = nb0; j < nb; j++)
	      {
		num_index(j) = 3;
		num_ref(j) = i;
	      }
            
            if (include_double_gradM_flow)
              {
                nb0 = nb;
                ref_beta(i).GetVaryingMedia(nb, rho_complex, num_component);
                for (int j = nb0; j < nb; j++)
                  {
                    num_index(j) = 4;
                    num_ref(j) = i;
                  }
              }
	  }
      }
  }
  
  
  //! computation of invRho and invMu from values of Rho and Mu
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::FinalizeComputationVaryingIndices()
  {
    // computation of inverses    
    for (int ref = 1; ref < ref_rho.GetM(); ref++)
      {
        ref_invRho(ref) = ref_rho(ref);
        ref_invRho(ref).GetInverse();
        ref_invMu(ref) = ref_mu(ref);
        ref_invMu(ref).GetInverse();        
      }
  }
  
  
  //! returns true if physical indexes are varying inside element i
  template<class Complexe, class Dimension>
  bool VarHelmholtz_Base<Complexe, Dimension>::IsVaryingMedia(int i) const
  {
    if (ref_rho(i).IsVarying()||ref_mu(i).IsVarying()||ref_sigma(i).IsVarying())
      return true;
    
    if (include_flow_term)
      if (ref_v0(i).IsVarying())
	return true;

    if (include_double_gradM_flow)
      if (ref_beta(i).IsVarying())
	return true;
    
    return false;
  }
  

  //! returns rho + i sigma/omega on a quadrature point
  template<>
  Complex_wp VarHelmholtz_Base<Complex_wp, Dimension2>
  ::GetRhoTilde(int ref, int num_elem, int npoint) const
  {
    Complex_wp rho_tilde = this->ref_rho(ref).GetCoefficient(var_problem, num_elem, npoint);
    rho_tilde += Iwp*this->ref_sigma(ref).GetCoefficient(var_problem, num_elem, npoint)
      /var_problem.GetOmega();
    return rho_tilde;
  }

  
  //! returns rho + i sigma/omega on a quadrature point
  template<>
  Complex_wp VarHelmholtz_Base<Complex_wp, Dimension3>
  ::GetRhoTilde(int ref, int num_elem, int npoint) const
  {
    Complex_wp rho_tilde = this->ref_rho(ref).GetCoefficient(var_problem, num_elem, npoint);
    rho_tilde += Iwp*this->ref_sigma(ref).GetCoefficient(var_problem, num_elem, npoint)
      /var_problem.GetOmega();
    return rho_tilde;
  }


  //! returns rho on a quadrature point
  template<>
  Real_wp VarHelmholtz_Base<Real_wp, Dimension2>
  ::GetRhoTilde(int ref, int num_elem, int npoint) const
  {
    return this->ref_rho(ref).GetCoefficient(var_problem, num_elem, npoint);
  }


  //! returns rho on a quadrature point
  template<>
  Real_wp VarHelmholtz_Base<Real_wp, Dimension3>
  ::GetRhoTilde(int ref, int num_elem, int npoint) const
  {
    return this->ref_rho(ref).GetCoefficient(var_problem, num_elem, npoint);
  }

  
  //! returns mu n \cdot n on a quadrature point 
  template<class Complexe, class Dimension>
  Complexe VarHelmholtz_Base<Complexe, Dimension>
  ::GetMuNormale(int ref, int num_elem, int npoint, const R_N& normale) const
  {
    TinyVector<Complexe, Dimension::dim_N> mu_n;
    this->ref_mu(ref).MltMatrix(var_problem, num_elem, npoint, normale, mu_n);
    return DotProd(normale, mu_n);
  }
  
  
  //! return true if additional unknown is required for absorbing boundary condition
  template<class Complexe, class Dimension>
  bool VarHelmholtz_Base<Complexe, Dimension>::UseAdditionalUnknownForABC() const
  {
    if (var_boundary.grazing_abc)
      {        
        if (var_boundary.GetOrderAbsorbingCondition() == 2)
          return true;
        
        if (var_boundary.GetOrderAbsorbingCondition() >= 4)
          return true;
      }
    
    if (var_boundary.GetOrderAbsorbingCondition() == 6)
      return true;

    //if (var_boundary.GetOrderAbsorbingCondition() == 23)
    //return true;
    
    return false;
  }
  

  //! returns true if the finite element matrix is symmetric
  template<class Complexe, class Dimension>
  bool VarHelmholtz_Base<Complexe, Dimension>::IsSymmetricProblem(bool eigen) const
  {
    if (include_flow_term)
      return false;
    
    if (var_problem.FirstOrderFormulation())
      {
	if (eigen)
	  return false;
        
        if (!this->use_symm_drude)
          if ((this->GetNbVectorialDofDrudeAll() > 0))
            return false;
        
	if (var_problem.GetLeafStaticCondensation())
	  return true;
	
	if ((var_boundary.GetNbGlobalEltPML() > 0) || (!var_problem.GetSymmetrizationUse()))
	  return false;
      }
    else
      {
        if (var_problem.GetLeafStaticCondensation())
	  return true;
        
        if (this->linearize_drude)
          if (this->GetNbVectorialDofDrudeAll() > 0)
            return false;
      }
    
    return true;
  }


  //! returns true if the mass matrix is symmetric
  template<class Complexe, class Dimension>
  bool VarHelmholtz_Base<Complexe, Dimension>::IsSymmetricMassMatrix() const
  {
    if (var_problem.FirstOrderFormulation())
      {
        if (this->GetNbVectorialDofDrudeAll() > 0)
          return false;
      }
    
    return true;
  }


  template<>
  void VarHelmholtz_Base<Real_wp, Dimension2>::
  ModifyPMLCoefficient(Real_wp& rho, Real_wp& sigma, TinyMatrix<Real_wp, Symmetric, 2, 2>& mu,
                       TinyMatrix<Real_wp, Symmetric, 2, 2>& invMu, TinyVector<Real_wp, 2>& v0,
                       int num_elem, int i1, int i) const
  {}

  
  //! modifies rho and mu due to PML layers
  template<>
  void VarHelmholtz_Base<Complex_wp, Dimension2>::
  ModifyPMLCoefficient(Complex_wp& rho, Complex_wp& sigma, TinyMatrix<Complex_wp, Symmetric, 2, 2>& mu,
                       TinyMatrix<Complex_wp, Symmetric, 2, 2>& invMu, TinyVector<Complex_wp, 2>& v0,
                       int num_elem, int i1, int i) const
  {    
    int num_pml = var_problem.mesh.Element(num_elem).GetNumberPML();
    if (var_problem.mesh.GetPmlArea(num_pml).GetRadiusPML() > 0)
      {
        // circular PML
        Complex_wp d = 1.0/var_boundary.GetTauPML(i1, i, 0);
        Complex_wp dbar = var_boundary.GetPrimitiveTauPML(i1, i, 1);
        
        Complex_wp exp_iteta = var_boundary.GetTauPML(i1, i, 1);
        Real_wp cos_teta = real(exp_iteta), sin_teta = imag(exp_iteta);
        
        Complex_wp dbar_d = dbar/d, d_dbar = d/dbar;

        // isotropic index only
        Complex_wp mu_r = mu(0, 0);
        
        rho *= d*dbar;
        sigma *= d*dbar;
        
        mu(0, 0) = mu_r*(dbar_d*square(cos_teta) + d_dbar*square(sin_teta));
        mu(1, 1) = mu_r*(dbar_d*square(sin_teta) + d_dbar*square(cos_teta));
        mu(0, 1) = mu_r*cos_teta*sin_teta*(dbar_d - d_dbar);
        
        GetInverse(mu, invMu);
        
        return;
      }

    Complex_wp dx = 1.0/var_boundary.GetTauPML(i1, i, 0);
    Complex_wp dy = 1.0/var_boundary.GetTauPML(i1, i, 1);
    rho *= dx*dy;
    sigma *= dx*dy;

    Complex_wp ratio = dx/dy;
    invMu(0, 0) *= ratio;    
    invMu(1, 1) /= ratio;
    
    mu(0, 0) /= ratio;
    mu(1, 1) *= ratio;

    v0(0) *= dy;
    v0(1) *= dx;
  }
  

  template<>
  void VarHelmholtz_Base<Real_wp, Dimension3>::
  ModifyPMLCoefficient(Real_wp& rho, Real_wp& sigma, TinyMatrix<Real_wp, Symmetric, 3, 3>& mu,
                       TinyMatrix<Real_wp, Symmetric, 3, 3>& invMu, TinyVector<Real_wp, 3>& v0,
                       int num_elem, int i1, int i) const
  {}

  
  //! modifies rho and mu due to PML layers
  template<>
  void VarHelmholtz_Base<Complex_wp, Dimension3>::
  ModifyPMLCoefficient(Complex_wp& rho, Complex_wp& sigma, TinyMatrix<Complex_wp, Symmetric, 3, 3>& mu,
                       TinyMatrix<Complex_wp, Symmetric, 3, 3>& invMu, TinyVector<Complex_wp, 3>& v0,
                       int num_elem, int i1, int i) const
  {
    Complex_wp dx = 1.0/var_boundary.GetTauPML(i1, i, 0);
    Complex_wp dy = 1.0/var_boundary.GetTauPML(i1, i, 1);
    Complex_wp dz = 1.0/var_boundary.GetTauPML(i1, i, 2);
    rho *= dx*dy*dz;
    sigma *= dx*dy*dz;
    
    invMu(0, 0) *= dx/(dy*dz);
    invMu(1, 1) *= dy/(dx*dz);
    invMu(2, 2) *= dz/(dx*dy);

    mu(0, 0) *= (dy*dz)/dx;
    mu(1, 1) *= (dx*dz)/dy;
    mu(2, 2) *= (dx*dy)/dz;

    v0(0) *= dy*dz;
    v0(1) *= dx*dz;
    v0(2) *= dx*dy;
  }

  
  //! allocation of arrays needed to compute finite element matrix
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::AllocateMassMatrices()
  {    
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    this->Glob_matMass_Dh.Reallocate(var_problem.mesh.GetNbElt());
    this->Glob_matMass_DhSigma.Reallocate(var_problem.mesh.GetNbElt());
    if (!var_problem.FirstOrderFormulation())
      this->Glob_matMass_Bh.Reallocate(var_problem.mesh.GetNbElt());
    else
      {
        this->Glob_matMass_DhSigmaDiff.Reallocate(var_boundary.GetNbEltPML());
        if (Dimension::dim_N == 3)
          {
            this->Glob_matMass_DhSigmaDiff2.Reallocate(var_boundary.GetNbEltPML());
            this->Glob_matMass_DhSigmaDiff3.Reallocate(var_boundary.GetNbEltPML());
          }
      }
    
    if (include_flow_term)
      this->Glob_matDamp_Mh.Reallocate(var_problem.mesh.GetNbElt());
    
    this->coefficient_mu_penalty.Reallocate(this->GetNbPhysicalIndices());
    this->coefficient_impedance_absorbing.Reallocate(this->GetNbPhysicalIndices());
    this->coefficient_impedance_absorbing.Fill(1.0);
    for (int i = 0; i < this->coefficient_impedance_absorbing.GetM(); i++)
      {
	Complexe rho = this->ref_rho(i), mu = this->ref_mu(i)(0,0);
        this->coefficient_impedance_absorbing(i) = sqrt(rho*mu);
	this->coefficient_mu_penalty(i) = MaxAbs(this->ref_mu(i).GetConstant());
      }
    
    Mlt(var_boundary.GetImpedanceCoefficientABC(), this->coefficient_impedance_absorbing);

    if ((var_problem.print_level >= 6) && (nb_proc == 1))
      {
        Vector<bool> RefUsed(this->GetNbPhysicalIndices()+1);
        RefUsed.Fill(false);
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          RefUsed(var_problem.mesh.Element(i).GetReference()) = true;
        
        cout << "Frequency = " << var_problem.GetOmega()/(2.0*pi_wp) << endl;
        cout << "Pulsation = " << var_problem.GetOmega() << endl;
        cout << "Wave vector = " << var_problem.GetWaveVector() << endl;
        for (int ref = 1; ref < RefUsed.GetM(); ref++)
          if (RefUsed(ref))
            {
              TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;
              mu = this->ref_mu(ref).GetConstant();
              Complexe rho(1), sigma(0);
              rho = this->ref_rho(ref).GetConstant();
              sigma = this->ref_sigma(ref).GetConstant();
              this->ref_drude(ref).ModifyCoefficientHelmholtz(var_problem.GetOmega(), rho, sigma, mu);
              cout << "Physical Media " << ref << endl;
              cout << "Value of mu = " << mu << endl;
              cout << "Value of rho = " << rho << endl;
              cout << "Value of sigma = " << sigma << endl;
              cout << "Value of rho + I sigma/omega = " << rho + Iwp*sigma/var_problem.GetOmega() << endl;
            }
      }    
  }
  
  
  //! Compute the number of dofs for Drude's model
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::ComputeDrudeDofs()
  {
    if (var_problem.FormulationDG() != ElementReference_Base::CONTINUOUS)
      return;
    
    this->nb_dof_drude_vec = 0;
    if (!var_problem.FirstOrderFormulation() && !this->linearize_drude)
      return;
    
    this->OffsetDofDrudeV.Reallocate(var_problem.mesh.GetNbElt()+1);
    this->OffsetDofDrudeV(0) = 0;
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	int nb_dof_v = 0;
        int ref = var_problem.mesh.Element(i).GetReference();
        if (this->ref_drude(ref).IsEnabled())
          {
            int nPole = this->ref_drude(ref).gamma.GetM();
            nb_dof_v = var_problem.GetNbPointsQuadratureInside(i);
	    if (this->ref_drude(ref).IsModeTE())
              nb_dof_v *= Dimension::dim_N;

            nb_dof_v *= nPole;
            this->nb_dof_drude_vec += nb_dof_v;
	  }
        
	this->OffsetDofDrudeV(i+1) = this->OffsetDofDrudeV(i) + nb_dof_v;
      }

    this->nb_dof_drude_vec_all = this->nb_dof_drude_vec;
  }

  
  //! Adds dofs due to pml for parallel computations
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::PutOtherGlobalDofs()
  {
    if (var_problem.FirstOrderFormulation())
      {
        MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
        int offset_loc = mesh_num.GetNbDof();
        int offset_glob = var_problem.GetNbGlobalMeshDof();
        int offset_loc2 = offset_loc + mesh_num.GetNbDofPML();
        int offset_glob2 = offset_glob + var_problem.GetNbGlobalDofPML();
        int d = Dimension::dim_N;
        
        if (mesh_num.GetNbDofPML() > 0)
          {    
            mesh_num.GlobDofNumber_Subdomain.Resize(var_problem.offset_dof_unknown(1));
            for (int i = 0; i < mesh_num.GetNbDof(); i++)
              {
                int npml = mesh_num.GetDofPML(i);
                if (npml >= 0)
                  {
                    mesh_num.GlobDofNumber_Subdomain(offset_loc + npml)
                      = offset_glob + mesh_num.GlobDofPML_Subdomain(npml);
                    
                    if (Dimension::dim_N == 3)
                      mesh_num.GlobDofNumber_Subdomain(offset_loc2 + npml)
                        = offset_glob2 + mesh_num.GlobDofPML_Subdomain(npml);
                  }
              }
            
	    Vector<IVect>& MatchingDof = var_problem.GetOriginalMatchingDofNeighbor();
            for (int num = 0; num < MatchingDof.GetM(); num++)
              {
                int nb_dof = MatchingDof(num).GetM();
                int nb_dof_pml = 0;
                for (int j = 0; j < nb_dof; j++)
                  if (mesh_num.GetDofPML(MatchingDof(num)(j)) >= 0)
                    nb_dof_pml++;
                
                if (nb_dof_pml > 0)
                  {
                    int nb_old = nb_dof;
                    MatchingDof(num).Resize(nb_old + nb_dof_pml*(d-1));
                    for (int j = 0; j < nb_old; j++)
                      {
                        int npml
                          = mesh_num.GetDofPML(MatchingDof(num)(j));
                        
                        if (npml >= 0)
                          {
                            MatchingDof(num)(nb_dof++) = offset_loc + npml;
                            if (d == 3)
                              MatchingDof(num)(nb_dof++) = offset_loc2 + npml;
                          }
                      }
                  }
              }
          }        
      }

#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&this->nb_dof_drude_vec, &this->nb_dof_drude_vec_all, 1, MPI_INT, MPI_SUM, var_problem.comm_group_mode);
#else
    this->nb_dof_drude_vec_all = this->nb_dof_drude_vec;
#endif
  }


  //! retrieves dof intervals whose signs change because of symmetrization
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    i0 = 0; i1 = 0; j0 = 0; j1 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	i0 = var_problem.offset_dof_unknown(1);
	i1 = i0 + var_problem.GetOffsetDofV(var_problem.mesh.GetNbElt());

        if (this->GetNbVectorialDofDrudeAll() > 0)
          {
            j0 = var_problem.GetNbDof() - 2*this->GetNbVectorialDofDrude();
            j1 = var_problem.GetNbDof() - this->GetNbVectorialDofDrude();
          }
      }
    else if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
	i0 = mesh_num.GetNbDof();
	i1 = i0 + mesh_num.GetNbDof()*Dimension::dim_N;
      }
    else
      {
	i0 = 0;
	i1 = mesh_num.GetNbDof();
	int Nvol = var_problem.GetOffsetDofV(var_problem.mesh.GetNbElt());
	j0 = i1 + Nvol;
	j1 = j0 + Nvol*Dimension::dim_N;
      }
  }
  

  //! retrieves coefficient used to switch between Laplace and Helmholtz equation for HDG
  template<class Complexe, class Dimension> template<class T>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::GetPenalizationTauIPP(const VarProblem<Dimension>& vars, int num_face,
                          T& coef_tau) const
  {
    if (vars.FirstOrderFormulation())
      {
	if (vars.upwind_fluxes)
	  coef_tau = -vars.alpha_penalization*var_problem.Glob_CoefPenalDG(num_face);
      }
    else
      {
	// symmetric formulation
	if (sizeof(Real_wp) != sizeof(Complexe))
	  {
	    if (vars.upwind_fluxes)
	      to_complex(Iwp*vars.alpha_penalization*
                         var_problem.Glob_CoefPenalDG(num_face), coef_tau);
	    else
	      to_complex(Iwp*vars.alpha_penalization, coef_tau);
	  }
	else
	  {
	    if (vars.upwind_fluxes)
	      coef_tau = -vars.alpha_penalization*var_problem.Glob_CoefPenalDG(num_face);
	  }
      }
  }


  //! retrieves coefficient used to switch between Laplace and Helmholtz equation for HDG
  template<class Complexe, class Dimension> template<class T>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::GetCoefficientIPP(const VarProblem<Dimension>& vars,
		      Complexe& m_iomega, Real_wp& coef_ipp, T& coef_u, T& coef_v,
                      T& coef_u_damped) const
  {
    SetComplexOne(coef_u); SetComplexOne(coef_v); coef_v = -coef_v;
    SetComplexOne(coef_u_damped);
    coef_ipp = -1.0;
    vars.GetMiomega(m_iomega);
    if (vars.FirstOrderFormulation())
      {
	// non-symmetric formulation (with positive mass matrix)
	if (sizeof(Real_wp) == sizeof(Complexe))
	  coef_v = 1.0;
	else
	  {
	    coef_u = m_iomega;
            coef_u_damped = 1.0;
	    coef_v = m_iomega;
	  }
      }
    else
      {
	// symmetric formulation
	coef_ipp = 1.0;
	if (sizeof(Real_wp) != sizeof(Complexe))
	  {
	    coef_u = -vars.GetOmega();
	    coef_v = -vars.GetOmega();
            to_complex(-Iwp, coef_u_damped);
	  }
      }
  }
  

  //! modification of coefficient due to boundary conditions
  template<class Complexe, class Dimension> template<class T>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::ModifyCoefficientBC_HDG(const VarProblem<Dimension>& vars,
			    const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
			    const Mesh<Dimension>& mesh, int ref_boundary, int k,
			    const T& s, const Real_wp& coef_ipp, const T& coef_v,
			    int ref_domain, int iquad, int num_loc, T& coef_abc) const
  {
    if (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_ABSORBING)
      coef_abc *= Real_wp(2);
    else if (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
        Complexe coef_imped; SetComplexOne(coef_imped);
	var_boundary.MltParamCondition(ref_boundary, 0, coef_imped);
	var_boundary.MltMuIntegrationByParts(ref_domain, iquad, num_loc, k, coef_imped);
	coef_abc += -coef_imped/coef_v*s;
      }
  }

  
  template<class Complexe, class Dimension>
  void VarHelmholtz_Base<Complexe, Dimension>::ComputeTauCoefficient()
  {
    if (var_problem.FormulationDG() != ElementReference_Base::HDG)
      return;

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    var_problem.Glob_CoefPenalDG.Zero();
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
        int ref = var_problem.mesh.Element(i).GetReference();
        if (this->ref_rho(ref).IsVarying() || this->ref_mu(ref).IsVarying())
          {
            int offset = var_problem.GetNbPointsQuadratureInside(i);
            for (int num_loc = 0; num_loc < var_problem.mesh.Element(i).GetNbBoundary(); num_loc++)
              {
                int nf = var_problem.mesh.Element(i).numBoundary(num_loc);
                int nb_pts_face = mesh_num.GetNbPointsQuadratureBoundary(nf);
                Real_wp rho, mu;
                for (int k = 0; k < nb_pts_face; k++)
                  {
                    rho = abs(ref_rho(ref).GetCoefficient(var_problem, i, offset+k));
                    mu = abs(ref_mu(ref).GetCoefficient(var_problem, i, offset+k)(0, 0));                    
                    var_problem.Glob_CoefPenalDG(nf) = max(var_problem.Glob_CoefPenalDG(nf), sqrt(rho*mu));
                  }

                offset += nb_pts_face;
              }
          }
        else
          for (int num_loc = 0; num_loc < var_problem.mesh.Element(i).GetNbBoundary(); num_loc++)
            {
              int nf = var_problem.mesh.Element(i).numBoundary(num_loc);
              var_problem.Glob_CoefPenalDG(nf) = max(var_problem.Glob_CoefPenalDG(nf),
                                                     abs(this->coefficient_impedance_absorbing(ref)));
            }
      }
  }


  //! computes the impedance coefficient
  template<class Complexe, class Dimension>
  Complexe VarHelmholtz_Base<Complexe, Dimension>::ComputeImpedanceCoefficient(const IVect& ref, const VectComplexe& U0)
  {
    IVect ref_cond(var_problem.mesh.GetNbReferences()+1);
    ref_cond.Fill(0);
    for (int i = 0; i < ref.GetM(); i++)
      ref_cond(ref(i)) = 1;
    
    int order = var_problem.GetMeshNumbering(0).GetOrder();
    Mesh<Dimension> mesh_subdiv;
    MeshInterpolationFEM<Dimension> interp_gamma;
    interp_gamma.SetGaussQuadrature(order);    
    interp_gamma.InitProjectionSurface(var_problem.mesh);
    interp_gamma.ComputeSurfaceMesh(ref_cond, var_problem.mesh, mesh_subdiv,
				    var_problem);
    
#ifdef SELDON_WITH_MPI
    interp_gamma.GatherQuadraturePoints(var_problem.comm_group_mode, false);
#else
    interp_gamma.GatherQuadraturePoints(false);
#endif

    VectComplexe trace_En;
    VectComplexe trace_Hn;
    var_problem.ComputeEnHnOnBoundary(interp_gamma, U0, trace_En, trace_Hn);

    int offset = 0;
    Complexe int_P, int_dP_dn; Real_wp int_one(0);
    SetComplexZero(int_P); SetComplexZero(int_dP_dn);
    for (int i = 0; i < interp_gamma.GetNbBoundary(); i++)
      {
        int Nquad = interp_gamma.GetNbPointsQuadrature(i);
        for (int k = 0; k < Nquad; k++)
          {
            Real_wp poids = interp_gamma.WeightsQuadrature(i, k);
            int_P += trace_En(offset + k)*poids;
            int_one += poids;
            int_dP_dn += trace_Hn(offset + k)*poids;
          }

        offset += Nquad;
      }

#ifdef SELDON_WITH_MPI
    if (var_problem.GetNbProcPerMode() > 1)
      {
        Complexe int_P0(int_P), int_dP0(int_dP_dn); Vector<int64_t> xtmp;
        Real_wp int_one0(int_one);
        MpiReduce(var_problem.comm_group_mode, &int_P0, xtmp, &int_P, 1, MPI_SUM, 0);
        MpiReduce(var_problem.comm_group_mode, &int_one0, xtmp, &int_one, 1, MPI_SUM, 0);
        MpiReduce(var_problem.comm_group_mode, &int_dP0, xtmp, &int_dP_dn, 1, MPI_SUM, 0);
        
        if (var_problem.GetRankProcMode() == 0)
          {
            DISP(int_P); DISP(int_one); DISP(int_dP_dn);    
          }
        
      }
#endif
   
    //DISP(int_one); DISP(int_P); DISP(int_dP_dn);
    //DISP(trace_En.GetM()); DISP(offset); DISP(trace_Hn.GetM());
    Complexe imped = int_P / (int_one*int_dP_dn);
    //DISP(imped);
    //exit(0);
    
    return imped;
  }
  

  //! computation of elementary matrix for Helmholtz equation (LDG formulation)
  /*!
    \param[in] iquad element number
    \param[out] num_dof degrees of freedom numbers
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] vars considered problem
    \param[in] Fb finite element associated with element iquad
   */  
  template<class Complexe, class Dimension> template<class T>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::ComputeElementaryMatrixHelmDG(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				  CondensationBlockSolver_Base<T>& solver_c,
				  const GlobalGenericMatrix<T>& nat_mat,
				  const VarProblem<Dimension>& vars,
				  const VarBoundaryCondition_Dim<Complexe, Dimension>& var_boundary,
				  const ElementReference<Dimension, 1>& Fb)
  {
    int nb_dof_elt = Fb.GetNbDof();
    int nb_dof_all = (Dimension::dim_N+1)*nb_dof_elt;
    
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
        
    // number of integration points
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    int ref_domain = vars.mesh.Element(iquad).GetReference();
    bool treat_inside_stiff = solver_c.TreatInsideStiffness();
    
    if (this->include_flow_term)
      {
        cout << "not implemented for dg" << endl;
        abort();
      }

    bool optim_condensation = false;
    int ic = solver_c.GetCondensedElementNumber();
    if (vars.GetLeafStaticCondensation())
      {
	CondensationBlockSolver_HelmDG<T, Complexe, Dimension>& solver
	  = dynamic_cast<CondensationBlockSolver_HelmDG<T, Complexe, Dimension>& >(solver_c);
	
	if ((Fb.LumpedMassMatrix()) && (!Fb.DiscontinuousElement()))
	  optim_condensation = true;

	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  {
	    int num_face = vars.mesh.Element(iquad).numBoundary(num_loc);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    if (rf != Fb.GetOrder())
	      optim_condensation = false;
	  }

	if(treat_inside_stiff)
	  solver.SetOptimizedCondensation(ic, optim_condensation);
      }

    // dof numbers
    num_dof.Reallocate(nb_dof_all);
    IVect Nodle = vars.GetDofNumberOnElement(iquad);
    int Nvol = mesh_num.GetNbDof();
    for (int j = 0; j < nb_dof_elt; j++)
      num_dof(j) = Nodle(j);

    int offset_v = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
    int offset_u0 = 0;
    int offset_uloc = nb_dof_elt;
    int offset_uloc2 = nb_dof_elt;
    int offset_vloc = nb_dof_elt;
    int offset_Lambda = 0;
    
    if (vars.FirstOrderFormulation())
      {
        if (vars.InsidePML(iquad))
          {
	    if (vars.FormulationDG() == ElementReference_Base::HDG)
	      {
		cout << "PML not implemented for HDG" << endl;
		abort();
	      }
	    
            nb_dof_all += nb_dof_elt*(Dimension::dim_N-1);
            offset_uloc2 += nb_dof_elt;
            offset_vloc += nb_dof_elt*(Dimension::dim_N-1);
            num_dof.Resize(nb_dof_all);
            for (int j = 0; j < nb_dof_elt; j++)
              num_dof(offset_uloc + j) = Nvol + mesh_num.GetDofPML(Nodle(j));
            
            if (Dimension::dim_N == 3)
              {
                for (int j = 0; j < nb_dof_elt; j++)
                  num_dof(offset_uloc2 + j) = offset_v + mesh_num.GetDofPML(Nodle(j));
              }
          }
	
        if (Dimension::dim_N == 3)
          offset_v += mesh_num.GetNbDofPML();
      }

    if (vars.FormulationDG() == ElementReference_Base::HDG)
      {
	// first Lambda, then u and v
	int nb_dof_lambda = mesh_num.Element(iquad).GetNbDof();
	
	if (optim_condensation)
	  {
	    nb_dof_all = nb_dof_lambda;
	    num_dof.Reallocate(nb_dof_lambda);
	    for (int j = 0; j < mesh_num.Element(iquad).GetNbDof(); j++)
	      num_dof(j) = mesh_num.Element(iquad).GetNumberDof(j);	    
	  }
	else
	  {
	    offset_u0 = nb_dof_lambda;
	    offset_vloc = nb_dof_elt + nb_dof_lambda;
	    
	    nb_dof_all = (Dimension::dim_N+1)*nb_dof_elt + nb_dof_lambda;
	    num_dof.Reallocate(nb_dof_all);
	    for (int j = 0; j < mesh_num.Element(iquad).GetNbDof(); j++)
	      num_dof(offset_Lambda+j) = mesh_num.Element(iquad).GetNumberDof(j);
	    
	    for (int j = 0; j < nb_dof_elt; j++)
	      for (int k = 0; k <= Dimension::dim_N; k++)
		num_dof(offset_u0 + j + k*nb_dof_elt) = vars.GetOffsetDofUnknown(k+1) + vars.GetOffsetDofV(iquad) + j;
	  }
      }
    else
      for (int j = 0; j < nb_dof_elt; j++)
	for (int k = 0; k < Dimension::dim_N; k++)
	  num_dof(offset_vloc + j + k*nb_dof_elt) = offset_v + Nodle(j) + k*Nvol;

    mat_interac.Clear();
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();
    
    bool variable = vars.UseNumericalIntegration(iquad);
    bool affine = vars.mesh.IsElementAffine(iquad);
    bool stiff = false;
    if (nat_mat.GetCoefStiffness() != Real_wp(0))
      stiff = true;
    
    // computation of rho omega Ji, omega Ji DFi^-1, mu^-1 omega Ji
    Complexe rho, sigma, rho_tilde;
    TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> invMu, mu;
    TinyVector<Complexe, Dimension::dim_N> v0;

    // penalization and coefficients (depending on signs and formulation used)
    Complexe m_iomega; Real_wp coef_ipp;
    T coef_u, coef_v, coef_u_damped;    
    this->GetCoefficientIPP(vars, m_iomega, coef_ipp, coef_u, coef_v,
                            coef_u_damped);
    
    if (vars.GetSymmetrizationUse() && vars.FirstOrderFormulation())
      {
	coef_ipp = 1.0;
	coef_v = -coef_v;
      }
    
    Real_wp jacob; Complexe tau;
    TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> ji_dfjm1;
    int i1 = iquad - vars.mesh.GetNbElt() + var_boundary.GetNbEltPML();

    if (optim_condensation)
      {
	TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> Btmp;

	CondensationBlockSolver_HelmDG<T, Complexe, Dimension>& solver
	  = dynamic_cast<CondensationBlockSolver_HelmDG<T, Complexe, Dimension>& >(solver_c);
	
	// we compute Dh and inverse of Bh
	Vector<T> Dh(nb_points_quad);
	Vector<TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> > invBhVol(nb_points_quad);
	if (treat_inside_stiff)
	  for (int i = 0; i < nb_points_quad; i++)
	    {
	      rho = this->ref_rho(ref_domain).GetCoefficient(vars, iquad, i);
	      sigma = this->ref_sigma(ref_domain).GetCoefficient(vars, iquad, i);
	      mu = this->ref_mu(ref_domain).GetCoefficient(vars, iquad, i);
	      // indices are modified inside PML
	      if ((vars.InsidePML(iquad)) && (!vars.FirstOrderFormulation()))
		this->ModifyPMLCoefficient(rho, sigma, mu, invMu, v0, iquad, i1, i);
	      
	      Real_wp poids = Fb.WeightsND(i);
	      Real_wp jacob_weighted = vars.Glob_jacobian(iquad)(0);
	      if (affine)
		{
		  jacob_weighted *= poids;
		  ji_dfjm1 = vars.Glob_DFjm1(iquad)(0);
		  Mlt(Real_wp(1) / vars.Glob_jacobian(iquad)(0), ji_dfjm1);
		}
	      else
		{
		  jacob_weighted = vars.Glob_jacobian(iquad)(i);
		  ji_dfjm1 = vars.Glob_DFjm1(iquad)(i);
		  Mlt(poids / jacob_weighted, ji_dfjm1);
		}
	      
	      Dh(i) = jacob_weighted*(coef_u*rho*nat_mat.GetCoefMass()
				      + coef_u_damped*sigma*nat_mat.GetCoefDamping());
	      
	      MltTrans(mu, ji_dfjm1, Btmp);
	      Mlt(ji_dfjm1, Btmp, invBhVol(i));
	      
	      invBhVol(i) *= jacob_weighted*square(nat_mat.GetCoefStiffness()) / (coef_v*nat_mat.GetCoefMass());
	    }
	
	int offset = 0;
	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  offset += Fb.GetNbQuadBoundary(num_loc);

	int nb_dof_boundary = Fb.GetNbDofBoundaries();
	Vector<T> C(nb_dof_boundary), Ca;
	
	Vector<T>* Cl_ptr;
	Vector<T> Cl_void;
	if(treat_inside_stiff)
	  Cl_ptr= &solver.GetCl(ic);
	else
	  Cl_ptr = &Cl_void;
	
	Vector<T>& Cl = *Cl_ptr;

	Vector<int>* NumU_ptr;
	Vector<int> NumU_void;
	if(treat_inside_stiff)
	  NumU_ptr= &solver.GetNumU(ic);
	else
	  NumU_ptr = &NumU_void;
	
	Vector<int>& NumU = *NumU_ptr;


	Vector<TinyVector<T, Dimension::dim_N> >* Sn_ptr;
	Vector<TinyVector<T, Dimension::dim_N> > Sn_void;
	if(treat_inside_stiff)
	  Sn_ptr= &solver.GetSn(ic);
	else
	  Sn_ptr = &Sn_void;
	
	Vector<TinyVector<T, Dimension::dim_N> >& Sn = *Sn_ptr;


	Vector<TinyVector<T, Dimension::dim_N> >* SnD_ptr;
	Vector<TinyVector<T, Dimension::dim_N> > SnD_void;
	if(treat_inside_stiff)
	  SnD_ptr= &solver.GetSnD(ic);
	else
	  SnD_ptr = &SnD_void;
	
	Vector<TinyVector<T, Dimension::dim_N> >& SnD = *SnD_ptr;

	Matrix<T> schur(offset, offset);
	schur.Zero();
	
	if (treat_inside_stiff)
	  {
	    Ca.Reallocate(offset);
	    Cl.Reallocate(offset);
	    Sn.Reallocate(nb_dof_boundary);
	    SnD.Reallocate(offset);
	    NumU.Reallocate(offset);
	    C.Zero(); Sn.Zero(); SnD.Zero(); Ca.Zero(); Cl.Zero();
	  }
	
	// loop over boundaries of the element
	Complexe phase, phase_conj;
	TinyVector<Real_wp, Dimension::dim_N> vec_ur;
	offset = 0;
	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  {
	    int num_face, ref_boundary, rf, rot; bool new_face;
	    const Matrix<int>& FacesQuadRotation =
	      var_boundary.GetGeometryPhaseData(iquad, num_loc,
						num_face, ref_boundary, rf, new_face, rot,
						phase, phase_conj);
            
	    if (rf != Fb.GetOrder())
	      {
		cout << "Variable order not implemented with optimized static condensation" << endl;
		abort();
	      }

	    Real_wp dsj(1); Complexe coef_tau;
            this->GetPenalizationTauIPP(vars, num_face, coef_tau);
	    for (int k = 0; k < Fb.GetNbQuadBoundary(num_loc); k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  dsj = vars.Glob_dsj(num_face)(k);
		else
		  dsj = vars.Glob_dsj(num_face)(krot);
		
		int kvol = Fb.GetQuadNumber(num_loc, k);
		if (treat_inside_stiff)
		  {
		    vec_ur = Fb.NormaleLoc(num_loc);		
		    vec_ur *= Fb.WeightsQuadratureBoundary(k, num_loc) / Fb.WeightsND(kvol);
		    Sn(kvol) += vec_ur;
		    
		    if ((ref_boundary == 0) ||
			(mesh.GetBoundaryCondition(ref_boundary) != BoundaryConditionEnum::LINE_DIRICHLET))
		      SnD(offset+k) += vec_ur;
		  }

		T coef_abc(coef_ipp*coef_tau);
		this->ModifyCoefficientBC_HDG(vars, var_boundary, mesh, ref_boundary,
					      k, nat_mat.GetCoefStiffness(), coef_ipp, coef_v,
					      ref_domain, iquad, num_loc, coef_abc);
		
		T vloc = dsj*Fb.WeightsQuadratureBoundary(k, num_loc)*nat_mat.GetCoefStiffness();
		if (treat_inside_stiff)
		  {
		    Ca(offset + k) = coef_abc*vloc;
		    Cl(offset + k) = coef_tau*vloc;
		    C(kvol) += coef_tau*vloc;
		    NumU(offset + k) = kvol;
		  }
		else
		  schur(offset+k, offset+k) += coef_abc*vloc;
	      }
	    
	    offset += Fb.GetNbQuadBoundary(num_loc);
	  }

	if (!treat_inside_stiff)
	  {
            mat_interac.Clear();
            mat_interac.Reallocate(schur.GetM(), schur.GetM());
            Vector<T> row_schur(schur.GetM());
	    for (int i = 0; i < schur.GetM(); i++)
              {
                schur.GetDenseRow(i, row_schur);
                mat_interac.SetDenseRow(i, row_schur);
              }
            
	    return;
	  }
	
	// then we compute matrix Btilde = D_h + (R* - Sn) B_h^{-1} (R - Sn*) + C
	Vector<TinyVector<T, Dimension::dim_N> > coefD(nb_points_quad);
	TinyVector<T, Dimension::dim_N> vec_u, vec_v;
	for (int i = 0; i < nb_dof_boundary; i++)
	  {
	    Mlt(invBhVol(i), Sn(i), vec_u);
	    Dh(i) += -coef_ipp*DotProd(vec_u, Sn(i)) + C(i);
	    coefD(i) = coef_ipp*vec_u;
	  }

	Matrix<T, Symmetric, RowSymPacked>& Btilde = solver.GetBtildeMatrix(ic);
	Btilde.Reallocate(nb_dof_elt, nb_dof_elt);
	Btilde.Zero();

	if (coef_ipp == Real_wp(1))
	  for (int i = 0; i < nb_points_quad; i++)
	    invBhVol(i) = -invBhVol(i);
	
	TinyVector<bool, 4> null_term(false, false, false, false);
	Fb.AddVariableElemMatrix(0, 0, Dh, invBhVol, coefD, coefD, null_term, Btilde);
	
	coefD.Clear(); Dh.Clear();
	if (coef_ipp == Real_wp(1))
	  for (int i = 0; i < nb_points_quad; i++)
	    invBhVol(i) = -invBhVol(i);
	
	Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > & invBh
	  = solver.GetInverseBh(ic);
	invBh.Reallocate(nb_dof_boundary);
	for (int i = 0; i < nb_dof_boundary; i++)
	  invBh(i) = invBhVol(i);

	invBhVol.Clear();

	// we compute Mint = -Cl^T + Sn^D B_h^{-1} (R - Sn^T)
	Matrix<T, General, ArrayRowSparse> Mint(offset, nb_dof_elt);
	Mint.Zero();

	VectReal_wp phi(nb_dof_elt); typename Dimension::VectR_N grad_phi(nb_dof_elt);
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    Fb.GetGradientSinglePhiQuadrature(i, phi, grad_phi);
	    for (int j = 0; j < offset; j++)
	      if ((!grad_phi(NumU(j)).IsZero()) || (NumU(j) == i))
		{
		  vec_v = grad_phi(NumU(j));
		  if (NumU(j) == i)
		    vec_v -= Sn(i);
		  
		  Mlt(invBh(NumU(j)), vec_v, vec_u);
		  Mint.AddInteraction(j, i, DotProd(vec_u, SnD(j)));
		}
	  }

	for (int j = 0; j < offset; j++)
	  Mint.AddInteraction(j, NumU(j), coef_ipp*Cl(j));
	
	Matrix<T> schur_tmp(nb_dof_elt, offset);

	GetInverse(Btilde); Transpose(Mint);

	Matrix<T> Bt(nb_dof_elt, nb_dof_elt);
	for (int i = 0; i < nb_dof_elt; i++)
	  for (int j = 0; j < nb_dof_elt; j++)
	    Bt(i, j) = Btilde(i, j);
	
	Mlt(Bt, Mint, schur_tmp);
	Transpose(Mint);
	Mlt(Mint, schur_tmp, schur);
	if (coef_ipp == Real_wp(1))
	  Mlt(Real_wp(-1), schur);
	
	for (int i = 0; i < offset; i++)
	  {
	    Mlt(invBh(NumU(i)), SnD(i), vec_u);
	    schur(i, i) += Ca(i);
	    for (int j = 0; j < offset; j++)
	      if (NumU(j) == NumU(i))
		schur(j, i) -= DotProd(SnD(j), vec_u);
	  }
	
	schur_tmp.Clear();

	// the elementary matrix is the Schur complement for static condensation
        mat_interac.Clear();
        mat_interac.Reallocate(schur.GetM(), schur.GetM());
        Vector<T> row_schur(schur.GetM());
        for (int i = 0; i < schur.GetM(); i++)
          {
            schur.GetDenseRow(i, row_schur);
            mat_interac.SetDenseRow(i, row_schur);
          }
        
	//mat_interac.Write("mat_schur.dat"); int test_input;
	//DISP(iquad); DISP(mat_interac); cout << "waiting" << endl; cin >> test_input;
	
	return;
      }
    
    if (variable)
      {
	Vector<T> MassRho(nb_points_quad);
        TinyVector<Vector<T>, Dimension::dim_N> MassRhoTau;
	TinyMatrix<Vector<T>, Symmetric, Dimension::dim_N, Dimension::dim_N> MassMu;
	TinyVector<Vector<TinyVector<T, Dimension::dim_N> >, Dimension::dim_N> MassDF;
	for (int j = 0; j < Dimension::dim_N; j++)
	  {
	    MassDF(j).Reallocate(nb_points_quad);
            MassRhoTau(j).Reallocate(nb_points_quad);
            for (int k = 0; k < Dimension::dim_N; k++)
	      MassMu(j, k).Reallocate(nb_points_quad);
	  }
	
	for (int i = 0; i < nb_points_quad; i++)
	  {
	    rho = this->ref_rho(ref_domain).GetCoefficient(vars, iquad, i);
	    sigma = this->ref_sigma(ref_domain).GetCoefficient(vars, iquad, i);
	    invMu = this->ref_invMu(ref_domain).GetCoefficient(vars, iquad, i);
            // indices are modified inside PML
            if ((vars.InsidePML(iquad)) && (!vars.FirstOrderFormulation()))
              this->ModifyPMLCoefficient(rho, sigma, mu, invMu, v0, iquad, i1, i);
            
	    Real_wp poids = Fb.WeightsND(i);
	    Real_wp jacob_weighted = vars.Glob_jacobian(iquad)(0);
	    if (affine)
	      {
		jacob_weighted *= poids;
		ji_dfjm1 = vars.Glob_DFjm1(iquad)(0);
	      }
	    else
	      {
		jacob_weighted = vars.Glob_jacobian(iquad)(i);
		ji_dfjm1 = vars.Glob_DFjm1(iquad)(i);
	      }
            
	    MassRho(i) = jacob_weighted*(coef_u*rho*nat_mat.GetCoefMass()
					 + coef_u_damped*sigma*nat_mat.GetCoefDamping());
            
	    for (int j = 0; j < Dimension::dim_N; j++)	  
	      for (int k = 0; k < Dimension::dim_N; k++)
		{
		  MassDF(j)(i)(k) = coef_ipp*ji_dfjm1(k, j)*poids*nat_mat.GetCoefStiffness(); 		  
		  MassMu(j, k)(i) = invMu(j, k)*jacob_weighted*coef_v*nat_mat.GetCoefMass();
		}
            
            if (vars.InsidePML(iquad))
              {
                rho_tilde = rho + sigma/m_iomega;
                if (vars.FirstOrderFormulation())
                  for (int k = 0; k < Dimension::dim_N; k++)
                    {
                      tau = var_boundary.GetTauPML(i1, i, k);
                      MassRhoTau(k)(i) = MassRho(i);
                      MassRhoTau(k)(i) += rho_tilde*jacob_weighted*tau*nat_mat.GetCoefDamping();
                      MassMu(k, k)(i) += invMu(k,k)*jacob_weighted*tau*nat_mat.GetCoefDamping();
                    }
              }
            
	  }

        // partie masse
        if (vars.InsidePML(iquad) && vars.FirstOrderFormulation())
          {
            Fb.AddVariableMassMatrix(offset_u0, offset_u0, MassRhoTau(0), mat_interac);
            Fb.AddVariableMassMatrix(offset_uloc, offset_uloc, MassRhoTau(1), mat_interac);
            if (Dimension::dim_N == 3)
              Fb.AddVariableMassMatrix(offset_uloc2, offset_uloc2, MassRhoTau(2), mat_interac);
          }
	else
          Fb.AddVariableMassMatrix(offset_u0, offset_u0, MassRho, mat_interac);
        
	for (int j = 0; j < Dimension::dim_N; j++)	  
	  for (int k = 0; k < Dimension::dim_N; k++)
	    Fb.AddVariableMassMatrix(offset_vloc + j*nb_dof_elt, offset_vloc+k*nb_dof_elt,
				     MassMu(j, k), mat_interac);
        
        // partie rigidite
	if (stiff && treat_inside_stiff)
	  {
	    Vector<T> ZeroMass; Vector<TinyVector<T, Dimension::dim_N> > ZeroVec;
	    Vector<TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> > ZeroStiff;
	    TinyVector<bool, 4> null_term(true, true, true, false);
	    for (int k = 0; k < Dimension::dim_N; k++)
	      Fb.AddVariableElemMatrix(offset_vloc+k*nb_dof_elt, offset_u0,
				       ZeroMass, ZeroStiff, ZeroVec, MassDF(k),
				       null_term, mat_interac);
	    
            if (vars.InsidePML(iquad) && vars.FirstOrderFormulation())
              {
                for (int k = 0; k < Dimension::dim_N; k++)
                  Fb.AddVariableElemMatrix(offset_vloc+k*nb_dof_elt, offset_uloc,
					   ZeroMass, ZeroStiff, ZeroVec, MassDF(k),
					   null_term, mat_interac);
                
                if (Dimension::dim_N == 3)
                  for (int k = 0; k < Dimension::dim_N; k++)
                    Fb.AddVariableElemMatrix(offset_vloc+k*nb_dof_elt, offset_uloc2,
					     ZeroMass, ZeroStiff, ZeroVec, MassDF(k),
					     null_term, mat_interac);                
              }
            
	    for (int k = 0; k < Dimension::dim_N; k++)
	      for (int m = 0; m < MassDF(k).GetM(); m++)
		MassDF(k)(m) *= coef_ipp;
	    
	    null_term(2) = false; null_term(3) = true;
            if (vars.InsidePML(iquad) && vars.FirstOrderFormulation())
              {
                for (int k = 0; k < Dimension::dim_N; k++)
                  Fb.AddVariableElemMatrix(offset_u0+k*nb_dof_elt, offset_vloc+k*nb_dof_elt,
					   ZeroMass, ZeroStiff, MassDF(k), ZeroVec,
					   null_term, mat_interac);
              }
            else
              for (int k = 0; k < Dimension::dim_N; k++)
                Fb.AddVariableElemMatrix(offset_u0, offset_vloc+k*nb_dof_elt,
					 ZeroMass, ZeroStiff, MassDF(k), ZeroVec,
					 null_term, mat_interac);
	  }	
      }
    else
      {
	T MassRho, zero; SetComplexZero(zero);
	TinyVector<T, Dimension::dim_N> zero_vec;
	TinyMatrix<T, General, Dimension::dim_N, Dimension::dim_N> zero_mat;
	TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> MassMu;
	TinyVector<TinyVector<T, Dimension::dim_N>, Dimension::dim_N> MassDF, OppMassDF;
	rho = this->ref_rho(ref_domain).GetConstant();
	sigma = this->ref_sigma(ref_domain).GetConstant();
	invMu = this->ref_invMu(ref_domain).GetConstant();
        
	jacob = vars.Glob_jacobian(iquad)(0);	
	ji_dfjm1 = vars.Glob_DFjm1(iquad)(0);
	
	MassRho = jacob*(coef_u*rho*nat_mat.GetCoefMass() 
			 + sigma*coef_u_damped*nat_mat.GetCoefDamping());
	
	for (int j = 0; j < Dimension::dim_N; j++)	  
	  for (int k = 0; k < Dimension::dim_N; k++)
	    {
	      MassDF(j)(k) = coef_ipp*ji_dfjm1(k, j)*nat_mat.GetCoefStiffness(); 
	      OppMassDF(j)(k) = coef_ipp*MassDF(j)(k);
	      
	      MassMu(j, k) = invMu(j, k)*jacob*coef_v*nat_mat.GetCoefMass();
	    }	
	
	Fb.AddConstantMassMatrix(offset_u0, offset_u0, MassRho, mat_interac);
	if (stiff && treat_inside_stiff)
	  {
	    for (int k = 0; k < Dimension::dim_N; k++)
              Fb.AddConstantElemMatrix(offset_vloc+k*nb_dof_elt, offset_u0,
				       zero, zero_mat, zero_vec, MassDF(k),
				       TinyVector<bool, 4>(true, true, true, false), mat_interac);

	    for (int k = 0; k < Dimension::dim_N; k++)
	      Fb.AddConstantElemMatrix(offset_u0, offset_vloc+k*nb_dof_elt,
				       zero, zero_mat, OppMassDF(k), zero_vec,
				       TinyVector<bool, 4>(true, true, false, true), mat_interac);
	  }
	
	for (int j = 0; j < Dimension::dim_N; j++)	  
	  for (int k = 0; k < Dimension::dim_N; k++)
	    Fb.AddConstantMassMatrix(offset_vloc + j*nb_dof_elt, offset_vloc+k*nb_dof_elt,
				     MassMu(j, k), mat_interac);
      }


    // part for HDG formulation
    if ((stiff) && (vars.FormulationDG() == ElementReference_Base::HDG))
      {
	typedef typename Dimension::DimensionBoundary DimensionB;
	R_N normale; Real_wp dsj;
	VectReal_wp poids; TinyVector<VectReal_wp, Dimension::dim_N> poidsNormale;
	T s = nat_mat.GetCoefStiffness(); Complexe phase, phase_conj; 
	
	// loop over boundaries of the element K
	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  {
	    int num_face, ref_boundary, rf, rot; bool new_face;
	    const Matrix<int>& FacesQuadRotation =
	      var_boundary.GetGeometryPhaseData(iquad, num_loc,
						num_face, ref_boundary, rf, new_face, rot,
						phase, phase_conj);

	    const ElementReference<DimensionB, 1>& Fb_s = vars.GetSurfaceFiniteElementH1(num_face);
	    
	    int Nquad = Fb_s.GetNbPointsQuadratureInside();
	    poids.Reallocate(Nquad);
	    for (int m = 0; m < Dimension::dim_N; m++)
	      poidsNormale(m).Reallocate(Nquad);
	    
            Complexe coef_tau;
            this->GetPenalizationTauIPP(vars, num_face, coef_tau);
	    for (int k = 0; k < Fb_s.GetNbPointsQuadratureInside(); k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  {
		    normale = vars.Glob_normale(num_face)(k);
		    dsj = vars.Glob_dsj(num_face)(k);
		  }
		else
		  {
		    normale = vars.Glob_normale(num_face)(krot);
		    dsj = vars.Glob_dsj(num_face)(krot);
		    normale = -normale;
		  }
		
		poids(k) = dsj*Fb_s.WeightsND(k);
		for (int m = 0; m < Dimension::dim_N; m++)
		  poidsNormale(m)(k) = poids(k)*normale(m);
	      }

	    VectReal_wp phi(Nquad), Ones(Fb.GetNbDof()), feval(Nquad), OnesS(Fb_s.GetNbDof());
	    // loops on dofs for volume unknowns u and v
	    if (treat_inside_stiff)
	      for (int j = 0; j < Fb.GetNbDof(); j++)
		if (Fb.DiscontinuousElement() || Fb.IsTangentialDof(j, num_loc))
		  {
		    Ones.Zero(); Ones(j) = 1.0;
		    Fb.ApplyShTranspose(num_loc, Ones, phi, rf);
		    
		    feval = poids*phi;
		    Ones.Zero();
		    Fb.ApplySh(Real_wp(1.0), num_loc, feval, Ones, rf);
		    
		    // part  coef_tau \int u phi
		    for (int k = 0; k < Fb.GetNbDof(); k++)
		      if (abs(Ones(k)) > vars.GetThresholdMatrix())
			mat_interac.AddInteraction(offset_u0+j, offset_u0+k, coef_tau*s*Ones(k));
		    
		    for (int m = 0; m < Dimension::dim_N; m++)
		      {
			feval = poidsNormale(m)*phi;
			Ones.Zero();
			Fb.ApplySh(Real_wp(1.0), num_loc, feval, Ones, rf);
			
			// part -v.n phi 
			// and part -c_ipp u psi.n
			for (int k = 0; k < Fb.GetNbDof(); k++)
			  if (abs(Ones(k)) > vars.GetThresholdMatrix())
			    {
			      mat_interac.AddInteraction(offset_u0+j, offset_vloc + k + m*nb_dof_elt, -s*Ones(k));
			      mat_interac.AddInteraction(offset_vloc + k + m*nb_dof_elt, offset_u0+j, -coef_ipp*s*Ones(k));
			    }
		      }
		  }
	    
	    // loop for surface unknown lambda
	    for (int j = 0; j < Fb_s.GetNbDof(); j++)
	      {
		Fb_s.GetValueSinglePhiQuadrature(j, phi);
				
		// part -coef_tau \int lambda phi
		// and -c_ipp coef_tau \int u q
		feval = poids*phi;
		Ones.Zero();
		Fb.ApplySh(Real_wp(1), num_loc, feval, Ones, rf);
		if (treat_inside_stiff)
		  for (int k = 0; k < Fb.GetNbDof(); k++)
		    if (abs(Ones(k)) > vars.GetThresholdMatrix())
		      {
			mat_interac.AddInteraction(offset_u0 + k, offset_Lambda+j, -coef_tau*s*phase*Ones(k));
			mat_interac.AddInteraction(offset_Lambda+j, offset_u0 + k, -coef_ipp*coef_tau*s*phase_conj*Ones(k));
		      }
		
		// part c_ipp coef_tau \int Lambda q
		feval = poids*phi;
		Fb_s.ApplyCh(feval, OnesS);

		T coef_abc(coef_ipp*coef_tau);
		this->ModifyCoefficientBC_HDG(vars, var_boundary, mesh, ref_boundary,
					      0, nat_mat.GetCoefStiffness(), coef_ipp, coef_v,
					      ref_domain, iquad, num_loc, coef_abc);
		
		for (int k = 0; k < Fb_s.GetNbDof(); k++)
		  mat_interac.AddInteraction(offset_Lambda+j, offset_Lambda+k, coef_abc*s*OnesS(k));
		
		if ((ref_boundary != 0) &&
		    (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_DIRICHLET))
		  continue;
		
		// following contributions are not present for Dirichlet dofs
		if (treat_inside_stiff)
		  for (int m = 0; m < Dimension::dim_N; m++)
		    {
		      feval = poidsNormale(m)*phi;
		      Ones.Zero();
		      Fb.ApplySh(Real_wp(1.0), num_loc, feval, Ones, rf);
		      // part c_ipp psi.n Lambda
		      // and part c_ipp v.n q		  
		      for (int k = 0; k < Fb.GetNbDof(); k++)
			if (abs(Ones(k)) > vars.GetThresholdMatrix())
			  {
			    mat_interac.AddInteraction(offset_vloc + k + m*nb_dof_elt, offset_Lambda+j, coef_ipp*s*phase*Ones(k));
			    mat_interac.AddInteraction(offset_Lambda+j, offset_vloc + k + m*nb_dof_elt, coef_ipp*s*phase_conj*Ones(k));
			  }
		    }		
	      }
	    
	    offset_Lambda += Fb_s.GetNbDof();            
	  }

	mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad, 1);
	mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad, 1);
      }
    
    //mat_interac.Write("mat_elem.dat"); int test_input;
    //DISP(iquad); DISP(mat_interac); cout << "waiting" << endl; cin >> test_input;
  }
  

  //! computation of elementary fluxes for Helmholtz equation and LDG formulation
  template<class Complexe, class Dimension> template<class T, class TypeEquation>
  void VarHelmholtz_Base<Complexe, Dimension>
  ::AddElementaryFluxesHelmDG(VirtualMatrix<T>& mat_sp, const GlobalGenericMatrix<T>& nat_mat,
			      const EllipticProblem<TypeEquation>& vars,
			      int offset_row, int offset_col)
  {
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::R_N R_N;
    
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);

    R_N normale;
    Real_wp dsj; MatrixN_N dfjm1;
    bool first_order = vars.FirstOrderFormulation();
    //bool affine = vars.mesh.IsElementAffine(iquad);
    
    int Nvol = mesh_num.GetNbDof();
    int Nvol_pml = mesh_num.GetNbDofPML();
    TinyVector<int, Dimension::dim_N+1> offset_dof;

    bool sym = vars.GetSymmetrizationUse();
    int nb_neighbor = 0;
    for (int iquad = 0; iquad < vars.mesh.GetNbElt(); iquad++)
      for (int num_pos1_face = 0; num_pos1_face < vars.mesh.Element(iquad).GetNbBoundary();
           num_pos1_face++)
	{
	  int num_face = vars.mesh.Element(iquad).numBoundary(num_pos1_face);
	  
          if (vars.FaceHasToBeConsideredForBoundaryIntegral(num_face))
            {
              int num_elem2 = vars.mesh.Boundary(num_face).numElement(0);
              int ref = vars.mesh.Boundary(num_face).GetReference();
              bool new_face = vars.IsNewFace(iquad)(num_pos1_face);
              if ((num_elem2 == iquad)&&(vars.mesh.Boundary(num_face).GetNbElements()==2))
                num_elem2 = vars.mesh.Boundary(num_face).numElement(1);

	      bool neighbor_face = (vars.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR);
	      
	      const ElementReference<Dimension, 1>& Fb = vars.GetReferenceElementH1(iquad);
	      const ElementReference<Dimension, 1>* Fb2_ptr = &vars.GetReferenceElementH1(num_elem2);
	      
	      int rf = mesh_num.GetOrderQuadrature(num_face);
	      int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
	      VectReal_wp PoidsFlux = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
	      
	      int ref_d = mesh.Element(iquad).GetReference();    
	      IVect Nodle = vars.GetDofNumberOnElement(iquad);
	      IVect Nodle2 = vars.GetDofNumberOnElement(num_elem2);
	      
	      offset_dof(0) = 0;
	      offset_dof(1) = Nvol;
	      if (vars.FirstOrderFormulation())
		offset_dof(1) += (Dimension::dim_N-1)*Nvol_pml;
	      
	      offset_dof(2) = offset_dof(1) + Nvol;
	      if (Dimension::dim_N == 3)
		offset_dof(3) = offset_dof(2) + Nvol;
	      
	      T cone(1), phase(1);
	      int rot1 = mesh.Element(iquad).GetOrientationBoundary(num_pos1_face), rot2 = 0;
	      
	      int cond = vars.mesh.GetBoundaryCondition(ref);
	      bool face_on_gamma = false;
	      if (cond != BoundaryConditionEnum::LINE_INSIDE)
		{
		  if (cond != BoundaryConditionEnum::LINE_NEIGHBOR)
		    face_on_gamma = true;
		}
	      
	      int num_pos2_face = -1;
	      int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0;
	      bool element2_pml = false;
	      if (num_elem2 != iquad)
		{
		  num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(num_face);
		  if (num_pos2_face < 0)
		    {
		      int nf2 = mesh_num.GetPeriodicBoundary(num_face);
		      vars.GetPeriodicPhase(num_face, phase);
		      num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(nf2);
		    }
		  
		  rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
		  rot = mesh_num.GetRotationFace(rot1, rot2, nv);
		  if (vars.InsidePML(num_elem2))
		    element2_pml = true;
		}
	      
	      int krot, num_point = Fb.GetNbPointsQuadratureInside();
	      for (int k = 0; k < num_pos1_face; k++)
		num_point += mesh_num.GetNbPointsQuadratureBoundary(mesh.Element(iquad).numBoundary(k));
	      
	      IVect NumRotQuad(nb_points_face);
	      int ref2 = mesh.Element(num_elem2).GetReference();
#ifdef SELDON_WITH_MPI
	      TinyVector<int, Dimension::dim_N+1> offset_dof_glob;
	      IVect Nodle2_PML;
	      int proc2 = -1;
	      int Nvol_all = vars.GetNbGlobalMeshDof();
	      int Nvol_all_pml = vars.GetNbGlobalDofPML();
	      if (neighbor_face)
		{            
		  int type_per = mesh_num.GetPeriodicityTypeForBoundary(num_face);
		  if (type_per >= 0)
		    vars.GetPeriodicPhase(num_face, phase);
		  
		  proc2 = vars.GetProcessorNeighboringFace(num_face);
		  num_pos2_face = vars.GetLocalPositionNeighboringFace(num_face);
                  int pos_loc_neighbor = mesh_num.GetLocalEdgeNumberNeighborElement(num_face);
                  Fb2_ptr = dynamic_cast<const ElementReference<Dimension, 1>* >(&vars.GetNeighborReferenceElement(pos_loc_neighbor));
		  rot = vars.GetRotationNeighboringFace(num_face);
		  Nodle2 = vars.GetNodleNeighboringFace(num_face);
		  Nodle2_PML = vars.GetNodlePmlNeighboringFace(num_face);
		  ref2 = vars.GetRefDomainNeighboringFace(num_face);
		  if (Nodle2_PML.GetM() > 0)
		    element2_pml = true;
		  
		  NumRotQuad.Reallocate(nb_points_face);
		  
		  offset_dof_glob(0) = 0;
		  offset_dof_glob(1) = Nvol_all;
		  if (vars.FirstOrderFormulation())
		    offset_dof_glob(1) += (Dimension::dim_N-1)*Nvol_all_pml;
		  
		  offset_dof_glob(2) = offset_dof_glob(1) + Nvol_all;
		  if (Dimension::dim_N == 3)
		    offset_dof_glob(3) = offset_dof_glob(2) + Nvol_all;        
		}
#endif

              const ElementReference<Dimension, 1>& Fb2 = *Fb2_ptr;
              
	      // loop over quadrature points of the face
	      TinyVector<T, Dimension::dim_N> tau_normale;
	      const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
		GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
	      
	      TinyMatrix<T, General, Dimension::dim_N+1, Dimension::dim_N+1> Nabc;
	      Vector<TinyMatrix<T, General, Dimension::dim_N+1, Dimension::dim_N+1> >
		MassIntra(nb_points_face), MassExtra(nb_points_face);
	      
	      Complexe m_iomega;
	      vars.GetMiomega(m_iomega);
	      
	      for (int k = 0; k < nb_points_face; k++)
		{
		  krot = FacesQuadRotation(rot, k);
		  
		  if (new_face)
		    {
		      normale = vars.Glob_normale(num_face)(k);
		      dsj = vars.Glob_dsj(num_face)(k);
		    }
		  else
		    {
		      normale = vars.Glob_normale(num_face)(krot);
		      dsj = vars.Glob_dsj(num_face)(krot);
		      Mlt(Real_wp(-1), normale);
		    }
		  
		  tau_normale = normale;
		  
		  bool add_penal = true;
		  if (face_on_gamma)
		    add_penal = false;
		  
		  if (vars.FirstOrderFormulation())
		    {
		      if (vars.InsidePML(iquad) || element2_pml)
			add_penal = false;
		    }
		  
		  for (int i = 0; i < Dimension::dim_N; i++)
		    {
		      if (first_order && (!sym))
			{
			  MassIntra(k)(0, i+1) = -tau_normale(i)*nat_mat.GetCoefStiffness();
			  MassExtra(k)(0, i+1) = -tau_normale(i)*nat_mat.GetCoefStiffness();
			  MassIntra(k)(i+1, 0) = tau_normale(i)*nat_mat.GetCoefStiffness();
			  MassExtra(k)(i+1, 0) = -tau_normale(i)*nat_mat.GetCoefStiffness();	    
			}
		      else
			{
			  MassIntra(k)(0, i+1) = -tau_normale(i)*nat_mat.GetCoefStiffness();
			  MassExtra(k)(0, i+1) = -tau_normale(i)*nat_mat.GetCoefStiffness();
			  MassIntra(k)(i+1, 0) = -tau_normale(i)*nat_mat.GetCoefStiffness();
			  MassExtra(k)(i+1, 0) = tau_normale(i)*nat_mat.GetCoefStiffness();	    
			}
		    }
		  
		  if (add_penal)
		    {
		      TypeEquation::GetPenalDG(Nabc, normale, iquad, num_point,
					       num_face, nat_mat, ref_d, ref2, vars, Fb);
		      
		      MassIntra(k) -= Nabc;
		      MassExtra(k) += Nabc;
		    }
		  
		  if (face_on_gamma)
		    {
		      TypeEquation::GetNabc(Nabc, normale, ref, iquad,
					    num_point, nat_mat, ref_d, vars, Fb);
		      
		      MassIntra(k) += Nabc;
		    }
		  
		  MassIntra(k) *= dsj*PoidsFlux(k);
		  MassExtra(k) *= dsj*PoidsFlux(k);
		  
		  NumRotQuad(k) = krot;
		  num_point++;
		}
	      
	      Vector<T> contrib(Fb.GetNbDof()), contrib2(Fb2.GetNbDof());
	      contrib.Fill(0); contrib2.Fill(0);
	      VectReal_wp val_phi(nb_points_face), val_phi_quad(Fb.GetNbPointsQuadratureInside());
	      Vector<T> feval(nb_points_face), feval_quad(Fb.GetNbPointsQuadratureInside());
	      Vector<T> feval_quad2(Fb2.GetNbPointsQuadratureInside());
	      VectReal_wp Ones(Fb.GetNbDof());
	      feval_quad.Fill(0); feval_quad2.Fill(0);
	      feval.Fill(0); val_phi.Fill(0); val_phi_quad.Fill(0);
	      
	      // loop over dofs
	      for (int i = 0; i < Fb.GetNbDof(); i++)
		{
		  Ones.Fill(0); Ones(i) = 1.0;
		  if (Fb.UseQuadraturePointsForSh())
		    {
		      Fb.ApplyChTranspose(Ones, val_phi_quad);
		      Fb.ApplyShQuadratureTranspose(num_pos1_face, val_phi_quad, val_phi, rf);
		    }
		  else
		    Fb.ApplyShTranspose(num_pos1_face, Ones, val_phi, rf);
		  
		  if (Norm2(val_phi) > 10.0*epsilon_machine)
		    {
		      for (int p = 0; p <= Dimension::dim_N; p++)
			for (int q = 0; q <= Dimension::dim_N; q++)
			  {
			    for (int k = 0; k < nb_points_face; k++)
			      feval(k) = val_phi(k)*MassIntra(k)(p, q);
			    
			    if (Fb.UseQuadraturePointsForSh())
			      {
				feval_quad.Fill(0);
				Fb.ApplyShQuadrature(cone, num_pos1_face, feval, feval_quad, rf);
				Fb.ApplyCh(feval_quad, contrib);
			      }
			    else
			      {
				contrib.Fill(0);
				Fb.ApplySh(cone, num_pos1_face, feval, contrib, rf);
			      }
			    
			    if ((vars.InsidePML(iquad)) && (vars.FirstOrderFormulation()))
			      {
				if ((p == 0) && (q > 0))
				  {
				    for (int j = 0; j < Fb.GetNbDof(); j++)
				      if (abs(contrib(j)) > vars.GetThresholdMatrix())
					{
					  int dofp = offset_row + Nodle(i);
					  if (q > 1)
					    dofp = Nvol + mesh_num.GetDofPML(dofp) + (q-2)*Nvol_pml;
					  
					  int dofq = offset_col + Nodle(j) + offset_dof(q);
					  mat_sp.AddInteraction(dofp, dofq, contrib(j));
					}
				  }
				else if ((p>0) && (q ==0))
				  {
				    for (int j = 0; j < Fb.GetNbDof(); j++)
				      if (abs(contrib(j)) > vars.GetThresholdMatrix())
					{
					  int dofp = offset_row + Nodle(i) + offset_dof(p);
					  int dofq = offset_col + Nodle(j);
					  int dofq_pml = offset_col + Nvol + mesh_num.GetDofPML(dofq);
					  mat_sp.AddInteraction(dofp, dofq, contrib(j));
					  mat_sp.AddInteraction(dofp, dofq_pml, contrib(j));
					  if (Dimension::dim_N == 3)
					    mat_sp.AddInteraction(dofp, dofq_pml + Nvol_pml, contrib(j));
					}
				  }
				else
				  {
				    for (int j = 0; j < Fb.GetNbDof(); j++)
				      if (abs(contrib(j)) > vars.GetThresholdMatrix())
					{
					  int dofp = offset_row + Nodle(i) + offset_dof(p);
					  int dofq = offset_col + Nodle(j) + offset_dof(q);
					  mat_sp.AddInteraction(dofp, dofq, contrib(j));
					}
				  }
			      }
			    else
			      {
				for (int j = 0; j < Fb.GetNbDof(); j++)
				  if (abs(contrib(j)) > vars.GetThresholdMatrix())
				    {
				      int dofp = offset_row + Nodle(i) + offset_dof(p);
				      int dofq = offset_col + Nodle(j) + offset_dof(q);
				      mat_sp.AddInteraction(dofp, dofq, contrib(j));
				    }
			      }
			    
			    // external contribution
			    if (num_pos2_face >= 0)
			      {
				for (int k = 0; k < nb_points_face; k++)
				  feval(NumRotQuad(k)) = val_phi(k)*MassExtra(k)(p, q);
				
				if (Fb2.UseQuadraturePointsForSh())
				  {
				    feval_quad2.Fill(0);
				    Fb2.ApplyShQuadrature(cone, num_pos2_face, feval, feval_quad2, rf);
				    Fb2.ApplyCh(feval_quad2, contrib2);
				  }
				else
				  {
				    contrib2.Fill(0);
				    Fb2.ApplySh(cone, num_pos2_face, feval, contrib2, rf);
				  }
				
				if (neighbor_face)
				  {                              
#ifdef SELDON_WITH_MPI                          
				    // interactions with another processor
				    if ((p == 0) && (q > 0))
				      {
					if ((vars.InsidePML(iquad)) && (vars.FirstOrderFormulation()))
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = Nodle(i);
						  if (q > 1)
						    dofp = Nvol + mesh_num.GetDofPML(dofp)
						      + (q-2)*Nvol_pml;
						  
						  int dofq = Nodle2(j) + offset_dof_glob(q);
						  mat_sp.AddDistantInteraction(dofp, dofq,
									       proc2, phase*contrib2(j));
						}
					  }
					else
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = Nodle(i) + offset_dof(p);
						  int dofq = Nodle2(j) + offset_dof_glob(q);
						  mat_sp.AddDistantInteraction(dofp, dofq, proc2,
									       phase*contrib2(j));
						}
					  }
				      }
				    else if ((p>0) && (q ==0))
				      {
					if (element2_pml && vars.FirstOrderFormulation())
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = Nodle(i) + offset_dof(p);
						  int dofq = Nodle2(j);
						  int dofq_pml = Nvol_all + Nodle2_PML(j);
						  mat_sp.AddDistantInteraction(dofp, dofq,
									       proc2, phase*contrib2(j));
						  
						  mat_sp.AddDistantInteraction(dofp, dofq_pml,
									       proc2, phase*contrib2(j));
						  if (Dimension::dim_N == 3)
						    mat_sp.AddDistantInteraction(dofp,
										 dofq_pml + Nvol_all_pml,
										 proc2, phase*contrib2(j));
						}
					  }
					else
					  {      
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = Nodle(i) + offset_dof(p);
						  int dofq = Nodle2(j) + offset_dof_glob(q);
						  mat_sp.AddDistantInteraction(dofp, dofq, proc2,
									       phase*contrib2(j));
						}                            
					  }
				      }
				    else
				      {
					// penalisation terms present only when no pml
					for (int j = 0; j < Fb2.GetNbDof(); j++)
					  if (abs(contrib2(j)) > vars.GetThresholdMatrix())
					    {
					      int dofp = Nodle(i) + offset_dof(p);
					      int dofq = Nodle2(j) + offset_dof_glob(q);
					      mat_sp.AddDistantInteraction(dofp, dofq, proc2,
									   phase*contrib2(j));
					    }                            
					
				      }
#endif
				  }
				else
				  {
				    if ((p == 0) && (q > 0))
				      {
					if (vars.InsidePML(iquad) && vars.FirstOrderFormulation())
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = offset_row + Nodle(i);
						  if (q > 1)
						    dofp = Nvol + mesh_num.GetDofPML(dofp)
						      + (q-2)*Nvol_pml;
						  
						  int dofq = offset_col + Nodle2(j) + offset_dof(q);
						  mat_sp.AddInteraction(dofp, dofq, phase*contrib2(j));
						}
					  }
					else
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = offset_row + Nodle(i) + offset_dof(p);
						  int dofq = offset_col + Nodle2(j) + offset_dof(q);
						  //DISP(dofp); DISP(dofq); DISP(phase); DISP(contrib2(j));
						  mat_sp.AddInteraction(dofp, dofq, phase*contrib2(j));
						}
					  }
				      }
				    else if ((p>0) && (q ==0))
				      {
					if (vars.InsidePML(num_elem2) && vars.FirstOrderFormulation())
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = offset_row + Nodle(i) + offset_dof(p);
						  int dofq = offset_col + Nodle2(j);
						  int dofq_pml = offset_col + Nvol + mesh_num.GetDofPML(dofq);
						  mat_sp.AddInteraction(dofp, dofq, phase*contrib2(j));
						  mat_sp.AddInteraction(dofp, dofq_pml, phase*contrib2(j));
						  if (Dimension::dim_N == 3)
						    mat_sp.AddInteraction(dofp, dofq_pml + Nvol_pml,
									  phase*contrib2(j));
						}
					  }
					else
					  {
					    for (int j = 0; j < Fb2.GetNbDof(); j++)
					      if (abs(contrib2(j)) > vars.GetThresholdMatrix())
						{
						  int dofp = offset_row + Nodle(i) + offset_dof(p);
						  int dofq = offset_col + Nodle2(j) + offset_dof(q);
						  //DISP(dofp); DISP(dofq); DISP(phase); DISP(contrib2(j));
						  mat_sp.AddInteraction(dofp, dofq, phase*contrib2(j));
						}
					  }
				      }
				    else
				      {
					for (int j = 0; j < Fb2.GetNbDof(); j++)
					  if (abs(contrib2(j)) > vars.GetThresholdMatrix())
					    {
					      int dofp = offset_row + Nodle(i) + offset_dof(p);
					      int dofq = offset_col + Nodle2(j) + offset_dof(q);
					      // DISP(dofp); DISP(dofq); DISP(phase); DISP(contrib2(j));
					      mat_sp.AddInteraction(dofp, dofq, phase*contrib2(j));
					    }
				      }
				  }
			      }
			  }
		    }
		} 
	      
	      if (neighbor_face)
		nb_neighbor++;
	    }
	}
  }
  

  template<> template<class Complexe, class T>
  void HelmholtzElementaryMatrixClass<Dimension2>
  ::ComputeMatrixV(const VarProblem<Dimension2>& vars,
		   const VarBoundaryCondition_Dim<Complexe, Dimension2>& var_boundary,
		   const VarHelmholtz_Base<Complexe, Dimension2>& var_helm,
		   int iquad, int i, int nb_points_quad, bool affine,
		   int offset_uloc, int offset_uloc2, int offset_vloc,
		   const GlobalGenericMatrix<T>& nat_mat, bool stiff, bool variable,
		   const VectReal_wp& val_phi, const Vector<TinyVector<Real_wp, 2> >& grad_phi,
		   Vector<T>& feval, Vector<T>& contrib, VirtualMatrix<T>& mat_interac)
  {
    Real_wp jacobian; T poids;
    TinyMatrix<T, Symmetric, 2, 2> invMu;

    Complexe m_iomega; int i1 = -1;
    vars.GetMiomega(m_iomega);

    bool sym = vars.GetSymmetrizationUse();
    TinyVector<T, 2> tau;
    const ElementReference<Dimension2, 1>& Fb = vars.GetReferenceElementH1(iquad);
    int ref_domain = vars.mesh.Element(iquad).GetReference();
    int nb_dof_elt = Fb.GetNbDof();
    
    for (int k = 0; k < nb_points_quad; k++)
      {
	if (affine)
	  jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	else
	  jacobian = vars.Glob_jacobian(iquad)(k);
	
	invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	invMu *= m_iomega*jacobian*nat_mat.GetCoefMass();
	if (stiff)
	  if ((abs(grad_phi(k)(0)) > vars.GetThresholdMatrix()) || (abs(grad_phi(k)(1)) > vars.GetThresholdMatrix()))    
	    {
	      poids = Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
	      mat_interac.AddInteraction(i, offset_vloc + 2*k, poids*grad_phi(k)(0));
	      mat_interac.AddInteraction(i, offset_vloc + 2*k+1, poids*grad_phi(k)(1));
	      
	      if (sym)
		poids *= -1.0;
	      
	      mat_interac.AddInteraction(offset_vloc + 2*k, i, -poids*grad_phi(k)(0));
	      mat_interac.AddInteraction(offset_vloc + 2*k+1, i, -poids*grad_phi(k)(1));
	    }
	
	if (i == 0)
	  {
	    if (sym)
	      invMu *= -1.0;
	    
	    mat_interac.AddInteraction(offset_vloc + 2*k, offset_vloc + 2*k, invMu(0, 0));
	    mat_interac.AddInteraction(offset_vloc + 2*k, offset_vloc + 2*k+1, invMu(1, 0));
	    mat_interac.AddInteraction(offset_vloc + 2*k+1, offset_vloc + 2*k, invMu(0, 1));
	    mat_interac.AddInteraction(offset_vloc + 2*k+1, offset_vloc + 2*k+1, invMu(1, 1));
	  }
      }
    
    if (vars.InsidePML(iquad))
      {
	i1 = iquad - vars.mesh.GetNbElt() + var_boundary.GetNbEltPML();
	for (int k = 0; k < nb_points_quad; k++)	  
	  {
	    invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	    if (affine)
	      jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	    else
	      jacobian = vars.Glob_jacobian(iquad)(k);
	    
	    if (stiff)
	      if ((abs(grad_phi(k)(0)) > vars.GetThresholdMatrix())
		  || (abs(grad_phi(k)(1)) > vars.GetThresholdMatrix()))
		{
		  poids = Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
		  mat_interac.AddInteraction(offset_uloc + i, offset_vloc + 2*k, poids*grad_phi(k)(0));
		  mat_interac.AddInteraction(offset_uloc + i, offset_vloc + 2*k+1, -poids*grad_phi(k)(1));
		}
	    
	    if (i == 0)
	      {
		tau = var_boundary.GetTauPML(i1, k);
		tau *= jacobian*nat_mat.GetCoefDamping();
		mat_interac.AddInteraction(offset_vloc + 2*k, offset_vloc + 2*k, tau(0)*invMu(0, 0));
		mat_interac.AddInteraction(offset_vloc + 2*k+1, offset_vloc + 2*k+1, tau(1)*invMu(1, 1));
	      }
	    
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(0)*Fb.WeightsND(k);
	    
	    feval(k) = poids*val_phi(k);
	  }
	
	Fb.ApplyCh(feval, contrib);
	for (int j = 0; j < nb_dof_elt; j++)
	  if (abs(contrib(j)) > vars.GetThresholdMatrix())
	    {
	      mat_interac.AddInteraction(i, offset_uloc+j, contrib(j));
	      mat_interac.AddInteraction(offset_uloc+i, j, contrib(j));
	    }
      }
  }    


  template<> template<class Complexe, class T>
  void HelmholtzElementaryMatrixClass<Dimension2>
  ::ComputeMatrixVdiag(const VarProblem<Dimension2>& vars,
                       const VarBoundaryCondition_Dim<Complexe, Dimension2>& var_boundary,
                       const VarHelmholtz_Base<Complexe, Dimension2>& var_helm,
                       int iquad, int nb_points_quad, bool affine, bool variable,
                       int offset_uloc, int offset_uloc2,
                       int offset_vloc, const GlobalGenericMatrix<T>& nat_mat,
                       VirtualMatrix<T>& mat_interac)
  {
    Real_wp jacobian; T poids;
    TinyMatrix<T, Symmetric, 2, 2> invMu;

    Complexe m_iomega; int i1 = -1;
    vars.GetMiomega(m_iomega);

    bool sym = vars.GetSymmetrizationUse();
    TinyVector<T, 2> tau;
    const ElementReference<Dimension2, 1>& Fb = vars.GetReferenceElementH1(iquad);
    int ref_domain = vars.mesh.Element(iquad).GetReference();
    
    for (int k = 0; k < nb_points_quad; k++)
      {
	if (affine)
	  jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	else
	  jacobian = vars.Glob_jacobian(iquad)(k);
	
	invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	invMu *= m_iomega*jacobian*nat_mat.GetCoefMass();
        if (sym)
          invMu *= -1.0;
        
        mat_interac.AddInteraction(offset_vloc + 2*k, offset_vloc + 2*k, invMu(0, 0));
        mat_interac.AddInteraction(offset_vloc + 2*k, offset_vloc + 2*k+1, invMu(1, 0));
        mat_interac.AddInteraction(offset_vloc + 2*k+1, offset_vloc + 2*k, invMu(0, 1));
        mat_interac.AddInteraction(offset_vloc + 2*k+1, offset_vloc + 2*k+1, invMu(1, 1));
      }
    
    if (vars.InsidePML(iquad))
      {
	i1 = iquad - vars.mesh.GetNbElt() + var_boundary.GetNbEltPML();
	for (int k = 0; k < nb_points_quad; k++)	  
	  {
	    invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	    if (affine)
	      jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	    else
	      jacobian = vars.Glob_jacobian(iquad)(k);
	    
            tau = var_boundary.GetTauPML(i1, k);
            tau *= jacobian*nat_mat.GetCoefDamping();
            mat_interac.AddInteraction(offset_vloc + 2*k, offset_vloc + 2*k, tau(0)*invMu(0, 0));
            mat_interac.AddInteraction(offset_vloc + 2*k+1, offset_vloc + 2*k+1, tau(1)*invMu(1, 1));
	    
            poids = nat_mat.GetCoefDamping();
            if (variable)
              poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(k);
            else
              poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(0)*Fb.WeightsND(k);
            
            mat_interac.AddInteraction(k, offset_uloc+k, poids);
            mat_interac.AddInteraction(offset_uloc+k, k, poids);
          }
      }
  }    


#ifdef MONTJOIE_WITH_THREE_DIM
  template<> template<class Complexe, class T>
  void HelmholtzElementaryMatrixClass<Dimension3>
  ::ComputeMatrixV(const VarProblem<Dimension3>& vars,
		   const VarBoundaryCondition_Dim<Complexe, Dimension3>& var_boundary,
		   const VarHelmholtz_Base<Complexe, Dimension3>& var_helm,
		   int iquad, int i, int nb_points_quad, bool affine,
		   int offset_uloc, int offset_uloc2, int offset_vloc,
		   const GlobalGenericMatrix<T>& nat_mat, bool stiff, bool variable,
		   const VectReal_wp& val_phi, const Vector<TinyVector<Real_wp, 3> >& grad_phi,
		   Vector<T>& feval, Vector<T>& contrib, VirtualMatrix<T>& mat_interac)
  {
    Real_wp jacobian; T poids;
    TinyMatrix<T, Symmetric, 3, 3> invMu;

    Complexe m_iomega; int i1 = -1;
    vars.GetMiomega(m_iomega);

    bool sym = vars.GetSymmetrizationUse();
    TinyVector<T, 3> tau;
    const ElementReference<Dimension3, 1>& Fb = vars.GetReferenceElementH1(iquad);
    int ref_domain = vars.mesh.Element(iquad).GetReference();
    int nb_dof_elt = Fb.GetNbDof();
    
    for (int k = 0; k < nb_points_quad; k++)
      {
	if (affine)
	  jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	else
	  jacobian = vars.Glob_jacobian(iquad)(k);
	
	invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	invMu *= m_iomega*jacobian*nat_mat.GetCoefMass();
	if (stiff)
	  if ((abs(grad_phi(k)(0)) > vars.GetThresholdMatrix())
	      || (abs(grad_phi(k)(1)) > vars.GetThresholdMatrix())
	      || (abs(grad_phi(k)(2)) > vars.GetThresholdMatrix()))	    
	    {
	      poids = Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
	      mat_interac.AddInteraction(i, offset_vloc + 3*k, poids*grad_phi(k)(0));
	      mat_interac.AddInteraction(i, offset_vloc + 3*k+1, poids*grad_phi(k)(1));
	      mat_interac.AddInteraction(i, offset_vloc + 3*k+2, poids*grad_phi(k)(2));
	      
	      if (sym)
		poids *= -1.0;
	      
	      mat_interac.AddInteraction(offset_vloc + 3*k, i, -poids*grad_phi(k)(0));
	      mat_interac.AddInteraction(offset_vloc + 3*k+1, i, -poids*grad_phi(k)(1));
	      mat_interac.AddInteraction(offset_vloc + 3*k+2, i, -poids*grad_phi(k)(2));
	    }
	
	if (i == 0)
	  {	    
	    if (sym)
	      invMu *= -1.0;
	    
	    mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k, invMu(0, 0));
	    mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k+1, invMu(1, 0));
	    mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k+2, invMu(2, 0));
	    mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k, invMu(0, 1));
	    mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k+1, invMu(1, 1));
	    mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k+2, invMu(2, 1));
	    mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k, invMu(0, 2));
	    mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k+1, invMu(1, 2));
	    mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k+2, invMu(2, 2));
	  }
      }
    
    if (vars.InsidePML(iquad))
      {
	i1 = iquad - vars.mesh.GetNbElt() + var_boundary.GetNbEltPML();
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	    if (affine)
	      jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	    else
	      jacobian = vars.Glob_jacobian(iquad)(k);
	    
	    if (stiff)
	      if ((abs(grad_phi(k)(0)) > vars.GetThresholdMatrix())
		  || (abs(grad_phi(k)(1)) > vars.GetThresholdMatrix())
		  || (abs(grad_phi(k)(2)) > vars.GetThresholdMatrix()))	    		
		{
		  poids = Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
		  mat_interac.AddInteraction(offset_uloc + i, offset_vloc + 3*k, poids*grad_phi(k)(0));
		  mat_interac.AddInteraction(offset_uloc + i, offset_vloc + 3*k+1, -poids*grad_phi(k)(1));
		  
		  mat_interac.AddInteraction(offset_uloc2 + i, offset_vloc + 3*k, poids*grad_phi(k)(0));
		  mat_interac.AddInteraction(offset_uloc2 + i, offset_vloc + 3*k+2, -poids*grad_phi(k)(2));
		}
	    
	    if (i == 0)
	      {		  
		tau = var_boundary.GetTauPML(i1, k);
		tau *= jacobian*nat_mat.GetCoefDamping();
		mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k, tau(0)*invMu(0, 0));
		mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k+1, tau(1)*invMu(1, 1));
		mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k+2, tau(2)*invMu(2, 2));
	      }
	    
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(0)*Fb.WeightsND(k);
	    
	    feval(k) = poids*val_phi(k);
	  }
	
	// part with sigma_x - sigma_y
	Fb.ApplyCh(feval, contrib);
	for (int j = 0; j < nb_dof_elt; j++)
	  if (abs(contrib(j)) > vars.GetThresholdMatrix())
	    {
	      mat_interac.AddInteraction(i, offset_uloc+j,  contrib(j));
	      mat_interac.AddInteraction(offset_uloc+i, j,  contrib(j));
	      mat_interac.AddInteraction(offset_uloc+i, offset_uloc2+j, contrib(j));
	    }

	// part with sigma_x - sigma_z
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff2(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff2(i1)(0)*Fb.WeightsND(k);
	    
	    feval(k) = poids*val_phi(k);
	  }
	
	Fb.ApplyCh(feval, contrib);
	for (int j = 0; j < nb_dof_elt; j++)
	  if (abs(contrib(j)) > vars.GetThresholdMatrix())
	    {
	      mat_interac.AddInteraction(i, offset_uloc2+j, contrib(j));
	      mat_interac.AddInteraction(offset_uloc2+i, j, contrib(j));
	      mat_interac.AddInteraction(offset_uloc2+i, offset_uloc+j, contrib(j));
	    }

	// part with sigma_y - sigma_z
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff3(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff3(i1)(0)*Fb.WeightsND(k);
	    
	    feval(k) = poids*val_phi(k);
	  }

	Fb.ApplyCh(feval, contrib);
	for (int j = 0; j < nb_dof_elt; j++)
	  if (abs(contrib(j)) > vars.GetThresholdMatrix())
	    {
	      mat_interac.AddInteraction(i, offset_uloc+j, -contrib(j));
	      mat_interac.AddInteraction(i, offset_uloc2+j, contrib(j));
	      mat_interac.AddInteraction(offset_uloc+i, offset_uloc+j, contrib(j));
	      mat_interac.AddInteraction(offset_uloc2+i, offset_uloc2+j, -contrib(j));
	    }                    
      }
  }

  
  template<> template<class Complexe, class T>
  void HelmholtzElementaryMatrixClass<Dimension3>
  ::ComputeMatrixVdiag(const VarProblem<Dimension3>& vars,
                       const VarBoundaryCondition_Dim<Complexe, Dimension3>& var_boundary,
                       const VarHelmholtz_Base<Complexe, Dimension3>& var_helm,
                       int iquad, int nb_points_quad, bool affine, bool variable,
                       int offset_uloc, int offset_uloc2,
                       int offset_vloc, const GlobalGenericMatrix<T>& nat_mat,
                       VirtualMatrix<T>& mat_interac)
  {
    Real_wp jacobian; T poids;
    TinyMatrix<T, Symmetric, 3, 3> invMu;

    Complexe m_iomega; int i1 = -1;
    vars.GetMiomega(m_iomega);

    bool sym = vars.GetSymmetrizationUse();
    TinyVector<T, 3> tau;
    const ElementReference<Dimension3, 1>& Fb = vars.GetReferenceElementH1(iquad);
    int ref_domain = vars.mesh.Element(iquad).GetReference();
    
    for (int k = 0; k < nb_points_quad; k++)
      {
	if (affine)
	  jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	else
	  jacobian = vars.Glob_jacobian(iquad)(k);
	
	invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	invMu *= m_iomega*jacobian*nat_mat.GetCoefMass();
        if (sym)
          invMu *= -1.0;
        
        mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k, invMu(0, 0));
        mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k+1, invMu(1, 0));
        mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k+2, invMu(2, 0));
        mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k, invMu(0, 1));
        mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k+1, invMu(1, 1));
        mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k+2, invMu(2, 1));
        mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k, invMu(0, 2));
        mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k+1, invMu(1, 2));
        mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k+2, invMu(2, 2));
      }
    
    if (vars.InsidePML(iquad))
      {
	i1 = iquad - vars.mesh.GetNbElt() + var_boundary.GetNbEltPML();
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    invMu = var_helm.ref_invMu(ref_domain).GetCoefficient(vars, iquad, k);
	    if (affine)
	      jacobian = vars.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	    else
	      jacobian = vars.Glob_jacobian(iquad)(k);
	    
            tau = var_boundary.GetTauPML(i1, k);
            tau *= jacobian*nat_mat.GetCoefDamping();
            mat_interac.AddInteraction(offset_vloc + 3*k, offset_vloc + 3*k, tau(0)*invMu(0, 0));
            mat_interac.AddInteraction(offset_vloc + 3*k+1, offset_vloc + 3*k+1, tau(1)*invMu(1, 1));
            mat_interac.AddInteraction(offset_vloc + 3*k+2, offset_vloc + 3*k+2, tau(2)*invMu(2, 2));
            
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff(i1)(0)*Fb.WeightsND(k);
	    
            // part with sigma_x - sigma_y
            mat_interac.AddInteraction(k, offset_uloc+k,  poids);
            mat_interac.AddInteraction(offset_uloc+k, k,  poids);
            mat_interac.AddInteraction(offset_uloc+k, offset_uloc2+k, poids);
          }

	// part with sigma_x - sigma_z
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff2(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff2(i1)(0)*Fb.WeightsND(k);
            
            mat_interac.AddInteraction(k, offset_uloc2+k, poids);
            mat_interac.AddInteraction(offset_uloc2+k, k, poids);
            mat_interac.AddInteraction(offset_uloc2+k, offset_uloc+k, poids);
          }
        
	// part with sigma_y - sigma_z
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    poids = nat_mat.GetCoefDamping();
	    if (variable)
	      poids *= var_helm.Glob_matMass_DhSigmaDiff3(i1)(k);
	    else
	      poids *= var_helm.Glob_matMass_DhSigmaDiff3(i1)(0)*Fb.WeightsND(k);
	    
            mat_interac.AddInteraction(k, offset_uloc+k, -poids);
            mat_interac.AddInteraction(k, offset_uloc2+k, poids);
            mat_interac.AddInteraction(offset_uloc+k, offset_uloc+k, poids);
            mat_interac.AddInteraction(offset_uloc2+k, offset_uloc2+k, -poids);
          }                    
      }
  }
#endif
  
  
  /*******************
   * VarHelmholtz_Eq *
   *******************/
  

  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>
  ::SetTypeEquation(const string& type_equation)
  {
    if ((type_equation == "HELMHOLTZ_SIPG") || (type_equation == "LAPLACE_SIPG")
        || (type_equation == "ACOUSTIC_SIPG"))
      {
        this->dg_formulation = ElementReference_Base::DISCONTINUOUS;
        HelmholtzEquation_Base<Complexe, Dimension>::store_dfjm1 = true;
        this->compute_dfjm1 = true;
      }
    else if ((type_equation == "HELMHOLTZ_HDG") || (type_equation == "LAPLACE_HDG")
	     || (type_equation == "ACOUSTIC_HDG"))
      {
	this->dg_formulation = ElementReference_Base::HDG;
      }
    
    if (this->dg_formulation != ElementReference_Base::CONTINUOUS)
      {
	this->alpha_penalization = -Real_wp(1);
	this->delta_penalization = -Real_wp(1);
      }

    VarHarmonic<TypeEquation>::SetTypeEquation(type_equation);
  }
  

  //! computation of the number of degrees of freedom
  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>::ComputeNumberOfDofs()
  {
    VarHarmonic<TypeEquation>::ComputeNumberOfDofs();

    if (this->FormulationDG() != ElementReference_Base::HDG)
      {
        this->ComputeDrudeDofs();
        
	if (this->FirstOrderFormulation())
          {
            int Nvol = this->mesh_num.GetNbDof();
            this->nodl += this->mesh_num.GetNbDofPML()*(Dimension::dim_N-1);    
            int nodl_scalar = Nvol + this->mesh_num.GetNbDofPML()*(Dimension::dim_N-1);
            this->nodl += 2*this->nb_dof_drude_vec;
            for (int m = 0; m <= this->nb_unknowns_scal; m++)
              this->offset_dof_unknown(m) = m*nodl_scalar;
            
            for (int m = this->nb_unknowns_scal+1; m <= this->nb_unknowns; m++)
              this->offset_dof_unknown(m) = this->offset_dof_unknown(m-1) + Nvol;
          }
        else if (this->linearize_drude)
          this->nodl += this->nb_dof_drude_vec;
      }
  }


  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation() && nat_mat.GetCoefStiffness() == Real_wp(0))
      {
        Vector<bool> diag_elt;
        int type_mat = this->GetMassMatrixType(diag_elt);
        if (type_mat == FemMassMatrix::DIAGONAL)
          return true;
      }

    return VarHarmonic<TypeEquation>::IsSparseElementaryMatrix(nat_mat);
  }


  template<class TypeEquation>  
  bool VarHelmholtz_Eq<TypeEquation>::IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation() && nat_mat.GetCoefStiffness() == Complex_wp(0, 0))
      {
        Vector<bool> diag_elt;
        int type_mat = this->GetMassMatrixType(diag_elt);
        if (type_mat == FemMassMatrix::DIAGONAL)
          return true;
      }

    return VarHarmonic<TypeEquation>::IsSparseElementaryMatrix(nat_mat);
  }

  
  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation())
      return false;
    
    if ((this->include_flow_term) && (nat_mat.GetCoefDamping() != Real_wp(0)))
      return false;

    return VarHarmonic<TypeEquation>::IsSymmetricElementaryMatrix(nat_mat);
  }
  

  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation())
      return false;
    
    if ((this->include_flow_term) && (nat_mat.GetCoefDamping() != Complex_wp(0, 0)))
      return false;

    return VarHarmonic<TypeEquation>::IsSymmetricElementaryMatrix(nat_mat);
  }
  
  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsSymmetricProblem(bool eigen) const
  {
    bool sym = VarHelmholtz_Base<Complexe, Dimension>::IsSymmetricProblem(eigen);
    if (!eigen)
      if (this->FirstOrderFormulation())
	const_cast<VarHelmholtz_Eq<TypeEquation>& >(*this).SetSymmetricElementaryMatrix(sym);
    
    if (!sym)
      return false;
    
    if (this->UseAdditionalUnknownForABC() && this->var_gibc.UseUnsymmetricImplementation())
      return false;
    
    return VarHarmonic<TypeEquation>::IsSymmetricProblem(eigen);
  }


  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {    
    Vector<bool> diag_elt;
    int type_mat = this->GetMassMatrixType(diag_elt);
    if (type_mat == FemMassMatrix::DIAGONAL)
      {
        if (this->FirstOrderFormulation() && (nat_mat.GetCoefDamping() != Real_wp(0)))
          return false;

        if (nat_mat.GetCoefStiffness() == Real_wp(0))
          return true;
      }

    return false;
  }
  

  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Complex_wp>& nat_mat) const
  {    
    Vector<bool> diag_elt;
    int type_mat = this->GetMassMatrixType(diag_elt);
    if (type_mat == FemMassMatrix::DIAGONAL)
      {
        if (this->FirstOrderFormulation() && (nat_mat.GetCoefDamping() != Complex_wp(0, 0)))
          return false;
                
        if (nat_mat.GetCoefStiffness() == Complex_wp(0, 0))
          return true;
      }

    return false;

  }

  
  template<class TypeEquation>
  bool VarHelmholtz_Eq<TypeEquation>::IsSymmetricMassMatrix() const
  {
    return VarHelmholtz_Base<Complexe, Dimension>::IsSymmetricMassMatrix();
  }

  
  //! returns the name associated with the physical index num
  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>::SetPhysicalIndexAtInfinity(const Vector<bool>& Ref)
  {
    VarHelmholtz_Base<Complexe, Dimension>::SetPhysicalIndexAtInfinity(Ref);
  }


  template<class TypeEquation>
  int VarHelmholtz_Eq<TypeEquation>::GetMassMatrixType(Vector<bool>& diag_elt) const
  {
    if (this->GetNbVectorialDofDrudeAll() > 0)
      {
        if (this->FirstOrderFormulation() || this->linearize_drude)
          return FemMassMatrix::MATRIX_FREE;
      }
    
    return VarHarmonic<TypeEquation>::GetMassMatrixType(diag_elt);
  }
  
  //! treating Dirichlet condition (retrieving dof numbers)  
  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>::TreatDirichletCondition()
  {    
    VarHarmonic<TypeEquation>::TreatDirichletCondition();      
    if (this->FirstOrderFormulation())
      {
	int d = Dimension::dim_N;
	int nb_dof = 0;
	for (int i = 0; i < this->Dirichlet_dof.GetM(); i++)
	  {
	    int j = this->Dirichlet_dof(i);
	    if (this->mesh_num.GetDofPML(j) >= 0)
	      nb_dof += d-1;
	  }
	
	if (nb_dof > 0)
	  {
	    int nb_old = this->Dirichlet_dof.GetM();
	    int Nvol = this->mesh_num.GetNbDof();
	    int Nvol_pml = this->mesh_num.GetNbDofPML();
	    this->Dirichlet_dof.Resize(nb_old + nb_dof);
            this->nb_dof_dirichlet = nb_old + nb_dof;
	    nb_dof = nb_old;
	    for (int i = 0; i < nb_old; i++)
	      {
		int j = this->Dirichlet_dof(i);
		if (this->mesh_num.GetDofPML(j) >= 0)
		  {
		    int k = this->mesh_num.GetDofPML(j);
		    this->Dirichlet_dof(nb_dof) = k + Nvol;
		    if (d == 3)
		      this->Dirichlet_dof(nb_dof+1) = k + Nvol + Nvol_pml;
		    
		    nb_dof += d-1;
		  }
	      }
	    
	    for (int i = nb_old; i < this->Dirichlet_dof.GetM(); i++)
	      this->is_dof_dirichlet(this->Dirichlet_dof(i)) = true;
	  }
      }
  }
  
  
  //! we add dofs coming from thin slots
  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>::PerformOtherInitializations()
  {
    VarHarmonic<TypeEquation>::PerformOtherInitializations();

    int n = this->nodl;
    if (this->FirstOrderFormulation()
        || (this->mesh_num.GetFormulationForPeriodicCondition() == this->mesh_num.WEAK_PERIODIC))
      this->compute_dfjm1 = true;

#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL     
    for (int i = 0; i < this->list_slots.GetM(); i++)
      if (this->list_slots(i).type_model == SlotModelParameters<Dimension2>::MODEL_MESH1D)
	{
	  this->list_slots(i).order = this->mesh_num.GetOrder();
	  this->list_slots(i).var1D.SetPulsation(this->GetOmega());
	  this->list_slots(i).ComputeVar1D();
	  this->list_slots(i).offset_nodl = this->nodl;
	  n += this->list_slots(i).var1D.GetNbDof() - 2;
	}
#endif
    
    if (n != this->nodl)
      this->ResizeNbDof(n);
  }


  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>::ComputeTauCoefficient()
  {
    VarHelmholtz_Base<Complexe, Dimension>::ComputeTauCoefficient();
  }


  template<class TypeEquation>
  int VarHelmholtz_Eq<TypeEquation>
  ::GetBoundaryConditionId(const IVect& ref, int pos, const VectString& parameters, bool& periodic)
  {
    int cond = ImpedanceABC_HelmholtzBase::
      GetBoundaryConditionId(ref, pos, parameters, periodic,
                             this->order_ABC, this->gamma_cla_coef, 
                             this->take_into_account_curvature_for_abc);
    
    if (cond >= 1)
      return cond;
    
    return VarHarmonic<TypeEquation>::GetBoundaryConditionId(ref, pos, parameters, periodic);
  }
  

  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const
  {
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(nb_points_quad);

    if (!this->FirstOrderFormulationDG())
      {
        R_N grad_u;
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    En_quad(j) = u_quadrature(0)(j);
	    if (compute_H)
	      {
                ExtractVector(grad_quadrature, j, 0, grad_u);
                Hn_quad(j) = DotProd(grad_u, normale(j));
              }
	  }
      }
    else
      {
	TinyVector<Real_wp, Dimension::dim_N> val_v;
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    En_quad(j) = u_quadrature(0)(j);
	    if (compute_H)
	      {
		ExtractVector(u_quadrature, j, 1, val_v);
		Hn_quad(j) = this->GetOmega()*DotProd(val_v, normale(j));
	      }
	  }
      }
  }


  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const
  {
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(nb_points_quad);

    TinyVector<Complex_wp, Dimension::dim_N> val_v;
    if (!this->FirstOrderFormulationDG())
      {
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    En_quad(j) = u_quadrature(0)(j);
	    if (compute_H)
	      {
                ExtractVector(grad_quadrature, j, 0, val_v);
                Hn_quad(j) = DotProd(val_v, normale(j));
              }
	  }
      }
    else
      {
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    En_quad(j) = u_quadrature(0)(j);
	    if (compute_H)
	      {
		ExtractVector(u_quadrature, j, 1, val_v);
		Hn_quad(j) = this->GetOmega()*DotProd(val_v, normale(j));
	      }
	  }
      }
  }

  
  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>
  ::ComputeEnHnNodal(Vector<VectReal_wp>& u_nodal, Vector<VectReal_wp>& grad_nodal,
		     int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
		     Vector<VectReal_wp>& En_nodal, Vector<VectReal_wp>& Hn_nodal) const
  {
    int nb_points_nodal = u_nodal(0).GetM();
    En_nodal.Reallocate(1); Hn_nodal.Reallocate(1);
    En_nodal(0).Reallocate(nb_points_nodal);
    Hn_nodal(0).Reallocate(nb_points_nodal);
    
    TinyVector<Real_wp, Dimension::dim_N> val_v;
    if (!this->FirstOrderFormulationDG())
      {
	for (int j = 0; j < nb_points_nodal; j++)
	  {
	    En_nodal(0)(j) = u_nodal(0)(j);
	    ExtractVector(grad_nodal, j, 0, val_v);
	    Hn_nodal(0)(j) = DotProd(val_v, normale(j));
	  }
      }
    else
      {
	for (int j = 0; j < nb_points_nodal; j++)
	  {
	    En_nodal(0)(j) = u_nodal(0)(j);
	    ExtractVector(u_nodal, j, 1, val_v);
	    Hn_nodal(0)(j) = this->GetOmega()*DotProd(val_v, normale(j));
	  }
      }    
  }
    
  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>
  ::ComputeEnHnNodal(Vector<VectComplex_wp>& u_nodal, Vector<VectComplex_wp>& grad_nodal,
		     int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
		     Vector<VectComplex_wp>& En_nodal, Vector<VectComplex_wp>& Hn_nodal) const
  {
    int nb_points_nodal = u_nodal(0).GetM();
    En_nodal.Reallocate(1); Hn_nodal.Reallocate(1);
    En_nodal(0).Reallocate(nb_points_nodal);
    Hn_nodal(0).Reallocate(nb_points_nodal);
    
    TinyVector<Complex_wp, Dimension::dim_N> val_v;
    if (this->FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	for (int j = 0; j < nb_points_nodal; j++)
	  {
	    En_nodal(0)(j) = u_nodal(0)(j);
	    ExtractVector(grad_nodal, j, 0, val_v);
	    Hn_nodal(0)(j) = DotProd(val_v, normale(j));
	  }
      }
    else
      {
	for (int j = 0; j < nb_points_nodal; j++)
	  {
	    En_nodal(0)(j) = u_nodal(0)(j);
	    ExtractVector(u_nodal, j, 1, val_v);
	    Hn_nodal(0)(j) = this->GetOmega()*DotProd(val_v, normale(j));
	  }
      }    
  }


  template<class TypeEquation>
  void VarHelmholtz_Eq<TypeEquation>
  ::FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const
  {
    VarHelmholtz_Base<Complexe, Dimension>::FindIntervalDofSignSymmetry(i0, i1, j0, j1);
  }
  
  
  /********************************
   * CondensationBlockSolver_Helm *
   ********************************/
  

  template<class T, class TypeEquation>
  CondensationBlockSolver_Helm<T, TypeEquation>
  ::CondensationBlockSolver_Helm(EllipticProblem<TypeEquation>& var)
    : CondensationBlockSolver_Fem<T>(var), var_problem(var)
  {
    
  }

  template<>
  void CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension2> >
  ::Init(const GlobalGenericMatrix<Real_wp>& nat_mat, bool compress_matrix)
  {
    abort();
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  void CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension3> >
  ::Init(const GlobalGenericMatrix<Real_wp>& nat_mat, bool compress_matrix)
  {
    abort();
  }
#endif

  
  template<class T, class TypeEquation>
  void CondensationBlockSolver_Helm<T, TypeEquation>
  ::Init(const GlobalGenericMatrix<T>& nat_mat, bool compress_matrix)
  {
    CondensationBlockSolver_Fem<T>::Init(nat_mat, compress_matrix);
    if (!var_problem.FirstOrderFormulation())
      return;
    
    if (var_problem.GetNbGlobalEltPML() <= 0)
      return;

    IVect& MatchingNumber_Pb = var_problem.GetProcMatchingNeighbor();
    Vector<IVect>& MatchingDofOrig_Pb = var_problem.GetOriginalMatchingDofNeighbor();

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int ns = MatchingNumber_Pb.GetM();
    IVect nb_dof(ns); nb_dof.Zero();
    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    ns = 0;
    for (int i = 0; i < MatchingNumber_Pb.GetM(); i++)
      {
	for (int k = 0; k < MatchingDofOrig_Pb(i).GetM(); k++)
	  if ((MatchingDofOrig_Pb(i)(k) >= Nvol)
	      && (MatchingDofOrig_Pb(i)(k) < Nvol+Npml))
	    nb_dof(i)++;

	if (nb_dof(i) > 0)
	  ns++;
      }

    this->SharingProcPML.Reallocate(ns);
    this->SharingRowsPML.Reallocate(ns);
    ns = 0;
    for (int i = 0; i < MatchingNumber_Pb.GetM(); i++)
      if (nb_dof(i) > 0)
	{
	  this->SharingRowsPML(ns).Reallocate(nb_dof(i));
	  nb_dof(i) = 0;
	  for (int k = 0; k < MatchingDofOrig_Pb(i).GetM(); k++)
	    if ((MatchingDofOrig_Pb(i)(k) >= Nvol)
		&& (MatchingDofOrig_Pb(i)(k) < Nvol+Npml))
	      {
		this->SharingRowsPML(ns)(nb_dof(i)) = MatchingDofOrig_Pb(i)(k) - Nvol;
		nb_dof(i)++;
	      }
	  
	  this->SharingProcPML(ns) = MatchingNumber_Pb(i);
	  ns++;
	}
    
    // retrieving tau on pml dofs
    num_u.Reallocate(mesh_num.GetNbDofPML()); num_u.Fill(-1);
    tau.Reallocate(mesh_num.GetNbDofPML());
    int offset = var_problem.mesh.GetNbElt() - var_problem.GetNbEltPML();
    if (var_problem.FirstOrderFormulation())
      {
        for (int i = offset; i < var_problem.mesh.GetNbElt(); i++)
          for (int k = 0; k < mesh_num.Element(i).GetNbDof(); k++)
            {
              int num_dof = mesh_num.Element(i).GetNumberDof(k);
              int npml = mesh_num.GetDofPML(num_dof);
              Complexe tau_x = var_problem.GetTauPML(i-offset, k)(0);
              Complexe tau_y = var_problem.GetTauPML(i-offset, k)(1);
              tau(npml).Init(realpart(tau_x), realpart(tau_y));
	      if (Dimension::dim_N == 3)
		tau(npml)(2) = realpart(var_problem.GetTauPML(i-offset, k)(2));
            }
        
        for (int i = 0; i < mesh_num.GetNbDof(); i++)
          {
            int npml = mesh_num.GetDofPML(i);
            if (npml >= 0)
              num_u(npml) = i;
          }
      }
    
    // computing mass_pml
    mass_pml.Reallocate(mesh_num.GetNbDofPML());
    mass_pml.Zero();
    
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(i);
	int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
	int ref_domain = var_problem.mesh.Element(i).GetReference();
	
	if (var_problem.InsidePML(i))
	  {
	    Vector<T> mass_elt;
	    if (var_problem.InsidePML(i))
	      mass_elt.Reallocate(nb_pts_quad);

	    bool affine = var_problem.mesh.IsElementAffine(i);
	    if (affine)
	      {
		for (int k = 0; k < nb_pts_quad; k++)
		  mass_elt(k) = Fb.WeightsND(k)*var_problem.Glob_jacobian(i)(0)*
		    var_problem.ref_rho(ref_domain).GetCoefficient(var_problem, i, k);
	      }
	    else
	      {
		for (int k = 0; k < nb_pts_quad; k++)
		  mass_elt(k) = var_problem.Glob_jacobian(i)(k)*
		    var_problem.ref_rho(ref_domain).GetCoefficient(var_problem, i, k);		
	      }

	    if (var_problem.InsidePML(i))
	      {
		for (int k = 0; k < nb_pts_quad; k++)
		  {
		    int num_dof = mesh_num.Element(i).GetNumberDof(k);
		    if (num_dof >= 0)
		      {
			int npml = mesh_num.GetDofPML(num_dof);
			mass_pml(npml) += mass_elt(k);
		      }
		  }
	      }
	  }
      }

    // assembling between processors
    this->AssemblePML(mass_pml);

#ifdef SELDON_WITH_MPI
    if (var_problem.GetNbProcPerMode() > 1)
      {
        VectReal_wp tau_scal(tau.GetM());
        for (int k = 0; k < Dimension::dim_N; k++)
          {
            for (int i = 0; i < tau.GetM(); i++)
              tau_scal(i) = tau(i)(k);
    
            AssembleVector(tau_scal, MPI_MAX, this->SharingProcPML, this->SharingRowsPML,
                           var_problem.comm_group_mode, mesh_num.GetNbDofPML(),
                           1, 21);
            
            for (int i = 0; i < tau.GetM(); i++)
              tau(i)(k) = tau_scal(i);
          }
      }
#endif

  }
  
  
  template<class T, class TypeEquation>
  void CondensationBlockSolver_Helm<T, TypeEquation>::AssemblePML(Vector<T>& x) const
  {
    if (var_problem.GetNbGlobalEltPML() <= 0)
      return;

#ifdef SELDON_WITH_MPI
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    AssembleVector(x, MPI_SUM, this->SharingProcPML, this->SharingRowsPML,
		   var_problem.comm_group_mode, mesh_num.GetNbDofPML(), 1, 21);
#endif
  }
  

  template<>
  void CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension2> >
  ::RecomposeSolution(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  void CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension3> >
  ::RecomposeSolution(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }
#endif



  template<class T, class TypeEquation>
  void CondensationBlockSolver_Helm<T, TypeEquation>
  ::RecomposeSolution(const SeldonTranspose& trans, Vector<T>& x_sol, const GlobalGenericMatrix<T>& nat_mat) const
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    IVect num_ddl;
    Matrix<T> mat_elt_unsym;
    
    int nb_elt = var_problem.GetNbElt();
    Matrix<T, Symmetric, RowSymPacked> a22_sym;
    Matrix<T> a12, a21, a22_unsym;

    Vector<T> rhs_pml, rhs_pml_tilde;
    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    bool sym = this->UseSymmetrization();
    if (var_problem.FirstOrderFormulation())
      {
 	rhs_pml.Reallocate(mesh_num.GetNbDofPML());
	if (Dimension::dim_N == 3)
	  rhs_pml_tilde.Reallocate(mesh_num.GetNbDofPML());
      }

    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    
    rhs_pml.Zero(); rhs_pml_tilde.Zero();
    var_problem.SetLeafStaticCondensation(true);
    for (int i = 0; i < nb_elt; i++)
      {	
	if (!this->store_block)
	  var_problem.ComputeElementaryMatrix(i, num_ddl, mat_elt_unsym,
					      const_cast<CondensationBlockSolver_Helm<T, TypeEquation>& >(*this),
					      nat_mat);

	// solution on internal nodes is recovered
	if (this->symmetric_elem_matrix)
	  {
	    if (this->store_block)
	      this->RecomposeSolGen(SeldonNoTrans, x_sol, i, this->block_sym(i).num_ddl, mat_elt_unsym,
				     this->block_sym(i).a12, this->block_sym(i).a21, this->block_sym(i).inv_a22);
	    else
	      this->RecomposeSolGen(SeldonNoTrans, x_sol, i, num_ddl, mat_elt_unsym, a12, a21, a22_sym);
	  }
	else
	  {
	    if (this->store_block)
	      this->RecomposeSolGen(trans, x_sol, i, this->block_unsym(i).num_ddl, mat_elt_unsym,
				    this->block_unsym(i).a12, this->block_unsym(i).a21, this->block_unsym(i).inv_a22);
	    else
	      this->RecomposeSolGen(trans, x_sol, i, num_ddl, mat_elt_unsym, a12, a21, a22_unsym);
	  }
	
	// solution on vectorial dofs is recovered
        if (var_problem.FirstOrderFormulation())
	  {
	    const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(i);
	    int nb_dof_loc = Fb.GetNbDof();
	    TinyVector<Vector<T>, 1> x;
	    x(0).Reallocate(nb_dof_loc); x(0).Zero();
	    var_problem.GetLocalUnknownVector(x_sol, i, x);
	    
	    int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
	    Vector<T> y(Dimension::dim_N*nb_pts_quad), ytilde;
	    if (Dimension::dim_N == 3)
	      if (var_problem.InsidePML(i))
		ytilde.Reallocate(3*nb_pts_quad);
            
	    Fb.ApplyRhTranspose(x(0), y);
	    
	    TinyVector<T, Dimension::dim_N> vec_u, vec_v, vec_w;
	    bool affine = var_problem.mesh.IsElementAffine(i); 
	    int ref_domain = var_problem.mesh.Element(i).GetReference();
	    int offset_v = Nvol + (Dimension::dim_N-1)*Npml;
            int offset_Q = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
            int offset_P = offset_Q - var_problem.GetNbVectorialDofDrude();
            offset_P += var_problem.OffsetDofDrudeV(i);
            offset_Q += var_problem.OffsetDofDrudeV(i);
            
            bool drude_te = false;
	    VectReal_wp gamma, omega0_2; T m(1), s(1), one, coef_drude;
            Vector<T> invDenom; int nPole = 0;
	    SetComplexOne(one); SetComplexOne(coef_drude);
	    if (var_problem.ref_drude(ref_domain).IsEnabled())
	      {
                Real_wp eps_inf = var_problem.ref_drude(ref_domain).eps_inf;
                nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
		gamma = var_problem.ref_drude(ref_domain).gamma;
		omega0_2 = var_problem.ref_drude(ref_domain).omega_02;
		m = nat_mat.GetCoefMass()*m_iomega; s = nat_mat.GetCoefStiffness();
                invDenom.Reallocate(nPole);
		for (int kp = 0; kp < nPole; kp++)
                  invDenom(kp) = Real_wp(1) / (m*(m+gamma(kp)*s) + omega0_2(kp)*s*s);
                
		if (var_problem.ref_drude(ref_domain).IsModeTE())
		  {
		    drude_te = true;
                    
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        Real_wp coef_omegap = var_problem.ref_drude(ref_domain).eps_omega_p2(kp) / eps_inf;
                        Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp) / eps_inf;
                        coef_drude += (s*s*coef_omegap + s*m*coef_sig) * invDenom(kp);
                      }

                    coef_drude = one/coef_drude;
		  }
	      }
            
	    int offset = offset_v + var_problem.GetOffsetDofV(i);
	    T coef = nat_mat.GetCoefStiffness()/(m_iomega*nat_mat.GetCoefMass());
	    T coef_mass = Real_wp(1) / (m_iomega*nat_mat.GetCoefMass());
	    T coef1, coef2, coef3;
	    int i1 = i - var_problem.mesh.GetNbElt() + var_problem.GetNbEltPML();	    
	    if (affine)
	      {
		coef1 = coef/var_problem.Glob_jacobian(i)(0);
		coef_mass /= var_problem.Glob_jacobian(i)(0);
		for (int k = 0; k < nb_pts_quad; k++)
		  {
		    CopyVector(y, k, vec_u);
		    MltTrans(var_problem.Glob_DFjm1(i)(0), vec_u, vec_v);
		    
		    ExtractVector(x_sol, offset + Dimension::dim_N*k, vec_w);
		    
		    if (var_problem.InsidePML(i))
		      {
			Real_wp tau_1 = realpart(var_problem.GetTauPML(i1, k)(0));
			Real_wp tau_2 = realpart(var_problem.GetTauPML(i1, k)(1));
			coef1 = Real_wp(1) / (var_problem.Glob_jacobian(i)(0)*(m_iomega*nat_mat.GetCoefMass() + tau_1*nat_mat.GetCoefDamping()));
			coef2 = Real_wp(1) / (var_problem.Glob_jacobian(i)(0)*(m_iomega*nat_mat.GetCoefMass() + tau_2*nat_mat.GetCoefDamping()));
			
			vec_u(0) = coef1*(nat_mat.GetCoefStiffness()*vec_v(0) + vec_w(0) / Fb.WeightsND(k));
			vec_u(1) = coef2*(nat_mat.GetCoefStiffness()*vec_v(1) + vec_w(1) / Fb.WeightsND(k));
			if (Dimension::dim_N == 3)
			  {
			    Real_wp tau_3 = realpart(var_problem.GetTauPML(i1, k)(2));
			    coef3 = Real_wp(1) / (var_problem.Glob_jacobian(i)(0)*(m_iomega*nat_mat.GetCoefMass() + tau_3*nat_mat.GetCoefDamping()));
			    vec_u(2) = coef3*(nat_mat.GetCoefStiffness()*vec_v(2) + vec_w(2) / Fb.WeightsND(k));
			  }
		      }
		    else
		      {
			coef2 = coef_mass / Fb.WeightsND(k);
			if (sym)
			  coef2 = -coef2;
			
			vec_u = coef1*vec_v + coef2*vec_w;
		      }
		    
		    Mlt(var_problem.ref_mu(ref_domain).GetCoefficient(var_problem, i, k), vec_u, vec_v);
		    if (drude_te)
		      vec_v  *= coef_drude;
		    
		    if (var_problem.InsidePML(i))
		      {
			if (Dimension::dim_N == 2)
			  {
			    vec_u(0) = vec_v(0);
			    vec_u(1) = -vec_v(1);
			  }
			else
			  {
			    vec_u(0) = vec_v(0);
			    SetComplexZero(vec_u(1));
			    vec_u(2) = -vec_v(2);

			    Mlt(var_problem.Glob_DFjm1(i)(0), vec_u, vec_w);			
			    vec_w *= Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
			    CopyVector(vec_w, k, ytilde);

			    vec_u(0) = vec_v(0);			    
			    vec_u(1) = -vec_v(1);
			    SetComplexZero(vec_u(2));
			  }
			
			Mlt(var_problem.Glob_DFjm1(i)(0), vec_u, vec_w);
			vec_w *= Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
			CopyVector(vec_w, k, y);
		      }
		    
		    ExtractVector(vec_v, offset + Dimension::dim_N*k, x_sol);
		  }
	      }
	    else
	      {
		for (int k = 0; k < nb_pts_quad; k++)
		  {
		    CopyVector(y, k, vec_u);
		    MltTrans(var_problem.Glob_DFjm1(i)(k), vec_u, vec_v);
		    
		    ExtractVector(x_sol, offset + Dimension::dim_N*k, vec_w);

		    if (var_problem.InsidePML(i))
		      {
			Real_wp tau_1 = realpart(var_problem.GetTauPML(i1, k)(0));
			Real_wp tau_2 = realpart(var_problem.GetTauPML(i1, k)(1));
			coef1 = Real_wp(1) / (var_problem.Glob_jacobian(i)(k)*(m_iomega*nat_mat.GetCoefMass() + tau_1*nat_mat.GetCoefDamping()));
			coef2 = Real_wp(1) / (var_problem.Glob_jacobian(i)(k)*(m_iomega*nat_mat.GetCoefMass() + tau_2*nat_mat.GetCoefDamping()));
			
			vec_u(0) = coef1*(nat_mat.GetCoefStiffness()*Fb.WeightsND(k)*vec_v(0) + vec_w(0));
			vec_u(1) = coef2*(nat_mat.GetCoefStiffness()*Fb.WeightsND(k)*vec_v(1) + vec_w(1));
			if (Dimension::dim_N == 3)
			  {
			    Real_wp tau_3 = realpart(var_problem.GetTauPML(i1, k)(2));
			    coef3 = Real_wp(1) / (var_problem.Glob_jacobian(i)(0)*(m_iomega*nat_mat.GetCoefMass() + tau_3*nat_mat.GetCoefDamping()));
			    vec_u(2) = coef3*(nat_mat.GetCoefStiffness()*Fb.WeightsND(k)*vec_v(2) + vec_w(2));
			  }			
		      }
		    else
		      {
			coef1 = coef*Fb.WeightsND(k) / var_problem.Glob_jacobian(i)(k);
			coef2 = coef_mass / var_problem.Glob_jacobian(i)(k);
			if (sym)
			  coef2 = -coef2;
			
			vec_u = coef1*vec_v + coef2*vec_w;
		      }
		    
		    Mlt(var_problem.ref_mu(ref_domain).GetCoefficient(var_problem, i, k), vec_u, vec_v);
		    if (drude_te)
		      vec_v  *= coef_drude;
		    
		    if (var_problem.InsidePML(i))
		      {
			if (Dimension::dim_N == 2)
			  {
			    vec_u(0) = vec_v(0);
			    vec_u(1) = -vec_v(1);
			  }
			else
			  {
			    vec_u(0) = vec_v(0);
			    SetComplexZero(vec_u(1));
			    vec_u(2) = -vec_v(2);

			    Mlt(var_problem.Glob_DFjm1(i)(k), vec_u, vec_w);			
			    vec_w *= Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
			    CopyVector(vec_w, k, ytilde);

			    vec_u(0) = vec_v(0);			    
			    vec_u(1) = -vec_v(1);
			    SetComplexZero(vec_u(2));
			  }

			Mlt(var_problem.Glob_DFjm1(i)(k), vec_u, vec_w);			
			vec_w *= Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
			CopyVector(vec_w, k, y);
		      }
		    
		    ExtractVector(vec_v, offset + Dimension::dim_N*k, x_sol);
		  }
	      }
	    
	    if (var_problem.InsidePML(i))
	      {
		Fb.ApplyRh(y, x(0));
		for (int k = 0; k < nb_dof_loc; k++)
		  {
		    int num_dof = mesh_num.Element(i).GetNumberDof(k);
		    if (num_dof >= 0)
		      {
			int npml = mesh_num.GetDofPML(num_dof);
			rhs_pml(npml) -= x(0)(k);
		      }
		  }
		
		if (Dimension::dim_N == 3)
		  {
		    Fb.ApplyRh(ytilde, x(0));
		    for (int k = 0; k < nb_dof_loc; k++)
		      {
			int num_dof = mesh_num.Element(i).GetNumberDof(k);
			if (num_dof >= 0)
			  {
			    int npml = mesh_num.GetDofPML(num_dof);
			    rhs_pml_tilde(npml) -= x(0)(k);
			  }
		      }
		  }
	      }

	    // we recover p and q for Drude's model
	    if (var_problem.ref_drude(ref_domain).IsEnabled())
	      {
		Vector<T> coef_E(nPole);
                for (int kp = 0; kp < nPole; kp++)
                  coef_E(kp) = (var_problem.ref_drude(ref_domain).eps_omega_p2(kp)*s
                                + var_problem.ref_drude(ref_domain).eps_sigma(kp)*m)*m*invDenom(kp);
                
		T coef_p = s/m;
		
		if (var_problem.ref_drude(ref_domain).IsModeTM())
		  {
		    for (int j = 0; j < nb_dof_loc; j++)
		      {
			int num_dof = mesh_num.Element(i).GetNumberDof(j);
                        for (int kp = 0; kp < nPole; kp++)
                          {
                            int num_dof_p = offset_P + j*nPole + kp;
                            int num_dof_q = offset_Q + j*nPole + kp;
                            x_sol(num_dof_q) += coef_E(kp)*x_sol(num_dof);
                            x_sol(num_dof_p) += coef_p*x_sol(num_dof_q);
                          }
                      }
		  }
		else
		  {
		    for (int kp = 0; kp < nPole; kp++)
                      for (int j = 0; j < nb_dof_loc; j++)
                        {
                          int num_dof_p = offset_P + 2*j*nPole + kp, num_dof_q = offset_Q + 2*j*nPole + kp;			
                          x_sol(num_dof_q) += coef_E(kp)*x_sol(offset+2*j);
                          x_sol(num_dof_p) += coef_p*x_sol(num_dof_q);

                          num_dof_p += nPole; num_dof_q += nPole;
                          x_sol(num_dof_q) += coef_E(kp)*x_sol(offset+2*j+1);
                          x_sol(num_dof_p) += coef_p*x_sol(num_dof_q);
                        }
		  }
	      }
	  }
        else
          {
            int ref_domain = var_problem.mesh.Element(i).GetReference();
            if (var_problem.ref_drude(ref_domain).IsEnabled() && var_problem.ref_drude(ref_domain).IsModeTM() && var_problem.linearize_drude)
              {
                int offset_P = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                offset_P += var_problem.OffsetDofDrudeV(i);
                
                T sig = nat_mat.GetCoefDamping();
                T m = nat_mat.GetCoefMass()*m_iomega;
                T s = nat_mat.GetCoefStiffness();
                
                const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(i);
                int nb_dof_loc = Fb.GetNbDof();
                bool affine = var_problem.mesh.IsElementAffine(i); 
                
                int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                T coef_E, coefDiag;
                Complexe rho = var_problem.ref_rho(ref_domain).GetConstant();
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                    Real_wp omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                    Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                    Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp); 
                    if (trans.NoTrans())
                      coef_E = coef_eps_inf*s + m_iomega*sig*coef_sig;
                    else
                      coef_E = -m_iomega*m;
                    
                    coefDiag = m_iomega*m + m_iomega*sig*gamma + s*omega0_2;
                    coef_E /= coefDiag;
                    
                    for (int j = 0; j < nb_dof_loc; j++)
                      {
                        int num_dof_p = offset_P + j*nPole + kp;
                        int num_dof = mesh_num.Element(i).GetNumberDof(j);
                        T coef_Eloc = coef_E;
                        if (var_problem.InsidePML(i) && !trans.NoTrans())
                          {
                            Real_wp jacob = var_problem.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
                            coef_Eloc *= var_problem.Glob_matMass_Dh(i)(j) / (m_iomega*m_iomega*jacob*rho);
                          }                      
                        
                        x_sol(num_dof_p) += coef_Eloc*x_sol(num_dof);
                      }
                  }
              }
          }
      }
    
    var_problem.SetLeafStaticCondensation(false);

    if (!var_problem.FirstOrderFormulation())
      return;

    // assembling between processors
    this->AssemblePML(rhs_pml);
    
    // we add source in rhs_pml
    for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
      rhs_pml(i) += x_sol(Nvol + i);

    if (Dimension::dim_N == 3)
      {
	this->AssemblePML(rhs_pml_tilde);
	for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
	  rhs_pml_tilde(i) += x_sol(Nvol + Npml + i);
      }

    // forming u^* = 1/(dt + (s1 + s2)/2 ) ( (s2 - s1)/2 u + rhs_pml / mass_pml ) in 2-D
    if (Dimension::dim_N == 2)
      for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
	{
	  x_sol(Nvol+i) = ( Real_wp(0.5)*(tau(i)(1) -tau(i)(0))*nat_mat.GetCoefDamping()*x_sol(num_u(i))
			    + rhs_pml(i) / mass_pml(i) )
	    / ( m_iomega*nat_mat.GetCoefMass() + Real_wp(0.5)*(tau(i)(1) + tau(i)(0))*nat_mat.GetCoefDamping());
	}
    else
      {
	TinyVector<T, 2> vec_u, vec_v;
	TinyMatrix<T, General, 2, 2> mat;
	for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
	  {
	    Real_wp sx = realpart(tau(i)(0));
	    Real_wp sy = realpart(tau(i)(1));
	    Real_wp sz = realpart(tau(i)(2));
	    T alpha = nat_mat.GetCoefDamping();	  
	    T u = x_sol(num_u(i));
	    vec_u(0) = u*alpha*(sy-sx)/3 + rhs_pml(i) / mass_pml(i);
	    vec_u(1) = u*alpha*(sz-sx)/3 + rhs_pml_tilde(i) / mass_pml(i);
	    
	    T diag = m_iomega*nat_mat.GetCoefMass();
	    mat(0, 0) = diag + alpha*(sx+2*sy)/3;
	    mat(0, 1) = alpha*(sx-sy)/3;
	    mat(1, 0) = alpha*(sx-sz)/3;
	    mat(1, 1) = diag + alpha*(sx+2*sz)/3;
	    GetInverse(mat);
	    Mlt(mat, vec_u, vec_v);
	    
	    x_sol(Nvol+i) = vec_v(0);
	    x_sol(Nvol+Npml + i) = vec_v(1);
	  }
      }
  }
  

  template<>
  void CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension2> >
  ::ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  void CondensationBlockSolver_Helm<Real_wp, HelmholtzEquation<Dimension3> >
  ::ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }
#endif


  template<class T, class TypeEquation>
  void CondensationBlockSolver_Helm<T, TypeEquation>
  ::ModifyRhsStaticCondensation(const SeldonTranspose& trans, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    //b_source.Write("b_before.dat");
    Vector<T> b_dir(var_problem.GetNbDirichletDof());
    for (int i = 0; i < b_dir.GetM(); i++)
      b_dir(i) = b_source(var_problem.GetDirichletDofNumber(i));
    
    IVect num_ddl;
    Matrix<T> mat_elt_unsym;
    
    int nb_elt = var_problem.GetNbElt();
    Matrix<T, Symmetric, RowSymPacked> a22_sym;
    Matrix<T> a12, a21, a22_unsym;
    
    // forming (dt + sigma_2) (f+f*)/2 + (dt + sigma_1) (f-f*/2)   (for PMLs)
    // (dt+sigma_y)(dt+sigma_z) (f+f*+fd)/3 + (dt+sigma_x)(dt+sigma_z)(f+fd-2 f*)/3 + (dt+sigma_x)(dt+sigma_y)(f+f*-2 fd) in 3-D
    int Nvol = mesh_num.GetNbDof();
    int Npml = mesh_num.GetNbDofPML();
    bool sym = this->UseSymmetrization(); Real_wp two(2), one(1);
    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    if (var_problem.FirstOrderFormulation())
      for (int i = 0; i < mesh_num.GetNbDof(); i++)
	if (mesh_num.GetDofPML(i) >= 0)
	  {
	    int ipml = mesh_num.GetDofPML(i);
	    T fstar = b_source(Nvol + ipml);
	    T f = b_source(i);	    
	    T s1 = tau(ipml)(0)*nat_mat.GetCoefDamping() / (nat_mat.GetCoefMass()*m_iomega);
	    T s2 = tau(ipml)(1)*nat_mat.GetCoefDamping() / (nat_mat.GetCoefMass()*m_iomega);
	    if (Dimension::dim_N == 2)
	      b_source(i) = (Real_wp(1) + Real_wp(0.5)*(s1+s2)) * f
		+ Real_wp(0.5)*(s2 - s1) * fstar;
	    else
	      {
		T s3 = tau(ipml)(2)*nat_mat.GetCoefDamping() / (nat_mat.GetCoefMass()*m_iomega);
		T fd = b_source(Nvol + Npml + ipml);
		b_source(i) = (one + s2)*(one + s3)*(f + fstar + fd) + (one + s1)*(one + s3)*(f + fd - two*fstar)
		  + (one + s1)*(one + s2)*(f + fstar - two*fd);

		b_source(i) /= Real_wp(3);
	      }
	  }

    var_problem.SetLeafStaticCondensation(true);
    for (int i = 0; i < nb_elt; i++)
      {	
	if (!this->store_block)
	  var_problem.
	    ComputeElementaryMatrix(i, num_ddl, mat_elt_unsym,
				    const_cast<CondensationBlockSolver_Helm<T, TypeEquation>& >(*this),
				    nat_mat);
	
	// part due to elimination of v
        bool affine = var_problem.mesh.IsElementAffine(i);
	if (var_problem.FirstOrderFormulation())
	  {
	    const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(i);
	    int nb_dof_loc = Fb.GetNbDof();
	    int ref_domain = var_problem.mesh.Element(i).GetReference();

	    int offset_v = Nvol + (Dimension::dim_N-1)*Npml;
            int offset_Q = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
            int offset_P = offset_Q - var_problem.GetNbVectorialDofDrude();
            offset_P += var_problem.OffsetDofDrudeV(i);
            offset_Q += var_problem.OffsetDofDrudeV(i);
            
	    int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
	    Vector<T> y(Dimension::dim_N*nb_pts_quad);
            
	    TinyVector<T, Dimension::dim_N> vec_u, vec_v;	    
	    int offset = offset_v + var_problem.GetOffsetDofV(i);
	    T coef = nat_mat.GetCoefStiffness()/(m_iomega*nat_mat.GetCoefMass());
	    if ((Dimension::dim_N == 3) && (var_problem.InsidePML(i)))
	      coef /= nat_mat.GetCoefMass();
	    
	    if (sym)
	      coef = -coef;
	    
	    int i1 = i - var_problem.mesh.GetNbElt() + var_problem.GetNbEltPML();
	    T one; SetComplexOne(one); T coef_drude(one);
	    VectReal_wp gamma, omega0_2; T m(0), s(0);
            Vector<T> invDenom, coef_p, coef_q;
	    bool drude_te = false; VectReal_wp coef_eps_inf; int nPole = 0;
	    if (var_problem.ref_drude(ref_domain).IsEnabled())
	      {
                Real_wp eps_inf = var_problem.ref_drude(ref_domain).eps_inf;
		gamma = var_problem.ref_drude(ref_domain).gamma;
                nPole = gamma.GetM();
		omega0_2 = var_problem.ref_drude(ref_domain).omega_02;
                coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2;
                
		m = nat_mat.GetCoefMass()*m_iomega; s = nat_mat.GetCoefStiffness();
                invDenom.Reallocate(nPole);
		for (int kp = 0; kp < nPole; kp++)
                  invDenom(kp) = one / (m*(m+gamma(kp)*s) + omega0_2(kp)*s*s);
                
		if (var_problem.use_symm_drude)
                  {
                    coef_p = -s*s*coef_eps_inf * invDenom;
                    coef_q = m*s * invDenom * coef_eps_inf;
                  }
                else
                  {
                    coef_q = m*s * invDenom;
                    coef_p = -s*s*omega0_2 * invDenom;
                  }
                
                if (sym)
                  coef_q = -coef_q;
                
                if (var_problem.ref_drude(ref_domain).IsModeTE())
                  {
                    drude_te = true;
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        Real_wp coef_omegap = var_problem.ref_drude(ref_domain).eps_omega_p2(kp) / eps_inf;
                        Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp) / eps_inf;
                        coef_drude += (s*s*coef_omegap + s*m*coef_sig) * invDenom(kp);
                      }
                    
                    coef_drude = one/coef_drude;
                  }
	      }
	    
	    if (affine)
	      {
		coef /= var_problem.Glob_jacobian(i)(0);
		for (int k = 0; k < nb_pts_quad; k++)
		  {
		    ExtractVector(b_source, offset + Dimension::dim_N*k, vec_u);
		    
		    if (var_problem.InsidePML(i))
		      {			
			Real_wp tau_1 = realpart(var_problem.GetTauPML(i1, k)(0));
			Real_wp tau_2 = realpart(var_problem.GetTauPML(i1, k)(1));
			T coef_damp = nat_mat.GetCoefDamping() / m_iomega;
			if (Dimension::dim_N == 2)
			  {
			    T ratio = (nat_mat.GetCoefMass() + tau_2*coef_damp) /
			      (nat_mat.GetCoefMass() + tau_1*coef_damp);
			    
			    vec_u(0) *= ratio;
			    vec_u(1) *= one/ratio;
			  }
			else
			  {
			    Real_wp tau_3 = realpart(var_problem.GetTauPML(i1, k)(2));
			    T sx = nat_mat.GetCoefMass() + tau_1*coef_damp;
			    T sy = nat_mat.GetCoefMass() + tau_2*coef_damp;
			    T sz = nat_mat.GetCoefMass() + tau_3*coef_damp;
			    vec_u(0) *= sy*sz / sx;
			    vec_u(1) *= sx*sz / sy;
			    vec_u(2) *= sx*sy / sz;
			  }
		      }

		    if (drude_te)
                      for (int kp = 0; kp < nPole; kp++)
                        {
                          int num_dof_P = offset_P + 2*k*nPole + kp;
                          int num_dof_Q = offset_Q + 2*k*nPole + kp;
                          Real_wp poids = var_problem.Glob_jacobian(i)(0)*Fb.WeightsND(k);
                          T fp = b_source(num_dof_P);
                          if (sym)
                            fp = -fp;

                          T vloc = coef_p(kp)*b_source(num_dof_P) + coef_q(kp)*b_source(num_dof_Q);
                        
                          vec_u(0) -= vloc;
                          if (var_problem.use_symm_drude)
                            b_source(num_dof_P) = fp*coef_eps_inf(kp) / (m*omega0_2(kp)*poids);
                          else
                            b_source(num_dof_P) = fp / (m*poids);
                          
                          if (sym)
                            b_source(num_dof_Q) = -vloc / (s*poids);
                          else
                            b_source(num_dof_Q) = vloc / (s*poids);
                          
                          num_dof_P += nPole; num_dof_Q += nPole;
                          fp = b_source(num_dof_P);
                          if (sym)
                            fp = -fp;

                          vloc = coef_p(kp)*b_source(num_dof_P) + coef_q(kp)*b_source(num_dof_Q);
                          vec_u(1) -= vloc;
                          if (var_problem.use_symm_drude)
                            b_source(num_dof_P) = fp*coef_eps_inf(kp) / (m*omega0_2(kp)*poids);
                          else
                            b_source(num_dof_P) = fp / (m*poids);
                          
                          if (sym)
                            b_source(num_dof_Q) = -vloc / (s*poids);
                          else
                            b_source(num_dof_Q) = vloc / (s*poids);
                        			
                          ExtractVector(vec_u, offset + Dimension::dim_N*k, b_source);
                        }
		    
		    Mlt(var_problem.ref_mu(ref_domain).GetCoefficient(var_problem, i, k), vec_u, vec_v);
		    Mlt(var_problem.Glob_DFjm1(i)(0), vec_v, vec_u);
		    
		    vec_u *= coef*coef_drude;
		    CopyVector(vec_u, k, y);

		  }
	      }
	    else
	      {
		for (int k = 0; k < nb_pts_quad; k++)
		  {
		    ExtractVector(b_source, offset + Dimension::dim_N*k, vec_u);
		    
		    if (var_problem.InsidePML(i))
		      {
			Real_wp tau_1 = realpart(var_problem.GetTauPML(i1, k)(0));
			Real_wp tau_2 = realpart(var_problem.GetTauPML(i1, k)(1));
			T coef_damp = nat_mat.GetCoefDamping() / m_iomega;
			if (Dimension::dim_N == 2)
			  {
			    T ratio = (nat_mat.GetCoefMass() + tau_2*coef_damp) /
			      (nat_mat.GetCoefMass() + tau_1*coef_damp);
			    
			    vec_u(0) *= ratio;
			    vec_u(1) *= one/ratio;
			  }
			else
			  {
			    Real_wp tau_3 = realpart(var_problem.GetTauPML(i1, k)(2));
			    T sx = nat_mat.GetCoefMass() + tau_1*coef_damp;
			    T sy = nat_mat.GetCoefMass() + tau_2*coef_damp;
			    T sz = nat_mat.GetCoefMass() + tau_3*coef_damp;
			    vec_u(0) *= sy*sz / sx;
			    vec_u(1) *= sx*sz / sy;
			    vec_u(2) *= sx*sy / sz;
			  }
		      }

		    if (drude_te)
                      for (int kp = 0; kp < nPole; kp++)
                        {
                          int num_dof_P = offset_P + 2*k*nPole + kp;
                          int num_dof_Q = offset_Q + 2*k*nPole + kp;
                          Real_wp poids = var_problem.Glob_jacobian(i)(k);
                          T fp = b_source(num_dof_P);
                          if (sym)
                            fp = -fp;
                          
                          T vloc = coef_p(kp)*b_source(num_dof_P) + coef_q(kp)*b_source(num_dof_Q);						
                          vec_u(0) -= vloc;
                          if (var_problem.use_symm_drude)
                            b_source(num_dof_P) = fp*coef_eps_inf(kp) / (m*omega0_2(kp)*poids);
                          else
                            b_source(num_dof_P) = fp / (m*poids);

                          if (sym)
                            b_source(num_dof_Q) = -vloc / (s*poids);
                          else
                            b_source(num_dof_Q) = vloc / (s*poids);

                          num_dof_P += nPole; num_dof_Q += nPole;
                          fp = b_source(num_dof_P);
                          if (sym)
                            fp = -fp;
                          
                          vloc = coef_p(kp)*b_source(num_dof_P) + coef_q(kp)*b_source(num_dof_Q);
                          vec_u(1) -= vloc;
                          if (var_problem.use_symm_drude)
                            b_source(num_dof_P) = fp*coef_eps_inf(kp) / (m*omega0_2(kp)*poids);
                          else
                            b_source(num_dof_P) = fp / (m*poids);

                          if (sym)
                            b_source(num_dof_Q) = -vloc / (s*poids);
                          else
                            b_source(num_dof_Q) = vloc / (s*poids);
                          
                          ExtractVector(vec_u, offset + Dimension::dim_N*k, b_source);
                        }
                    
		    Mlt(var_problem.ref_mu(ref_domain).GetCoefficient(var_problem, i, k), vec_u, vec_v);
		    Mlt(var_problem.Glob_DFjm1(i)(k), vec_v, vec_u);
		    
		    vec_u *= coef*coef_drude*Fb.WeightsND(k) / var_problem.Glob_jacobian(i)(k);				
		    CopyVector(vec_u, k, y);
		  }
	      }
	    
	    
	    TinyVector<Vector<T>, 1> x;
	    x(0).Reallocate(nb_dof_loc); x(0).Zero();
	    
	    Fb.ApplyRh(y, x(0));
	    
	    // part due to elimination of p, q (Drude's model)
	    if (var_problem.ref_drude(ref_domain).IsEnabled())
	      if (var_problem.ref_drude(ref_domain).IsModeTM())
		{
		  for (int j = 0; j < nb_dof_loc; j++)
		    {
                      Real_wp jacobian = var_problem.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
                      for (int kp = 0; kp < nPole; kp++)
                        {
                          int num_dof_p = offset_P + j*nPole + kp;
                          int num_dof_q = offset_Q + j*nPole + kp;
                          T fp = b_source(num_dof_p);
                          T vloc = coef_p(kp)*b_source(num_dof_p) + coef_q(kp)*b_source(num_dof_q);
                          x(0)(j) += vloc;
                          b_source(num_dof_q) = vloc / (s*jacobian);
                          if (var_problem.use_symm_drude)
                            b_source(num_dof_p) = fp*coef_eps_inf(kp)/(m*omega0_2(kp)*jacobian);
                          else
                            b_source(num_dof_p) = fp / (m*jacobian);
                        }
		    }
		}
	    
	    var_problem.AddLocalUnknownVector(-one, x, i, b_source);
	  }
        else
          {
            int ref_domain = var_problem.mesh.Element(i).GetReference();
            if (var_problem.ref_drude(ref_domain).IsEnabled() && var_problem.ref_drude(ref_domain).IsModeTM() && var_problem.linearize_drude)
              {
                int offset_P = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                offset_P += var_problem.OffsetDofDrudeV(i);
                
                T sig = nat_mat.GetCoefDamping();
                T m = nat_mat.GetCoefMass()*m_iomega;
                T s = nat_mat.GetCoefStiffness();
                
                const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(i);
                int nb_dof_loc = Fb.GetNbDof();
                
                int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                T coef_E, coefDiag;
                Complexe rho = var_problem.ref_rho(ref_domain).GetConstant();
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                    Real_wp omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                    Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                    Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp); 
                    if (trans.NoTrans())
                      coef_E = -m_iomega*m;
                    else
                      coef_E = coef_eps_inf*s + m_iomega*sig*coef_sig;

                    coefDiag = m_iomega*m + m_iomega*sig*gamma + s*omega0_2;
                    coef_E /= coefDiag;
                    
                    for (int j = 0; j < nb_dof_loc; j++)
                      {
                        Real_wp poids = var_problem.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
                        T coef_Eloc = coef_E;
                        if (var_problem.InsidePML(i) && trans.NoTrans())
                          coef_Eloc *= var_problem.Glob_matMass_Dh(i)(j) / (m_iomega*m_iomega*poids*rho);
                        
                        int num_dof_p = offset_P + j*nPole + kp;
                        int num_dof = mesh_num.Element(i).GetNumberDof(j);
                        T f_P = b_source(num_dof_p);
                        b_source(num_dof) += coef_Eloc*f_P;
                        b_source(num_dof_p) = f_P / (coefDiag*poids);
                      }
                  }                
              }
          }
        
	// part due to elimination of internal nodes of u
	if (this->symmetric_elem_matrix)
	  {
	    if (this->store_block)
	      this->ModifyRhsGen(SeldonNoTrans, b_source, i, this->block_sym(i).num_ddl, mat_elt_unsym,
				 this->block_sym(i).a12, this->block_sym(i).a21, this->block_sym(i).inv_a22);
	    else
	      this->ModifyRhsGen(SeldonNoTrans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_sym);
	  }
	else
	  {
	    if (this->store_block)
	      this->ModifyRhsGen(trans, b_source, i, this->block_unsym(i).num_ddl, mat_elt_unsym,
				 this->block_unsym(i).a12, this->block_unsym(i).a21, this->block_unsym(i).inv_a22);
	    else
	      this->ModifyRhsGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_unsym);
	  }
      }

    for (int i = 0; i < b_dir.GetM(); i++)
      b_source(var_problem.GetDirichletDofNumber(i)) = b_dir(i);

    var_problem.SetLeafStaticCondensation(false);
    //b_source.Write("b_after.dat");
  }


  //! returns the memory used by the object in bytes
  template<class T, class TypeEquation>
  size_t CondensationBlockSolver_Helm<T, TypeEquation>::GetMemorySize() const
  {
    size_t taille = CondensationBlockSolver_Fem<T>::GetMemorySize();
    taille += num_u.GetMemorySize() + tau.GetMemorySize() + mass_pml.GetMemorySize();
    taille += SharingProcPML.GetMemorySize() + SharingRowsPML.GetMemorySize();
    return taille;
  }

  
  /**********************************
   * CondensationBlockSolver_HelmDG *
   **********************************/


  //! returns the memory used by the object in bytes
  template<class T, class Dimension>
  size_t CondensedBlockSolver_HelmDG<T, Dimension>::GetMemorySize() const
  {
    size_t taille = invBtilde.GetMemorySize() + Sn.GetMemorySize() + SnD.GetMemorySize()
      + NumU.GetMemorySize() + invBh.GetMemorySize() + Cl.GetMemorySize();
    return taille;
  }
  
  
  //! sets an optimized condensation (or not) for element i
  template<class T, class Complexe, class Dimension>
  void CondensationBlockSolver_HelmDG<T, Complexe, Dimension>::SetOptimizedCondensation(int i, bool optim)
  {
    if (optim_elt.GetM() <= 0)
      optim_elt.Reallocate(this->GetNbCondensedElt());
    
    optim_elt(i).flag = optim;
  }


  template<>
  void CondensationBlockSolver_HelmDG<Real_wp, Complex_wp, Dimension2>
  ::ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  void CondensationBlockSolver_HelmDG<Real_wp, Complex_wp, Dimension3>
  ::ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }
#endif

  template<class T, class Complexe, class Dimension>
  void CondensationBlockSolver_HelmDG<T, Complexe, Dimension>
  ::ModifyRhsStaticCondensation(const SeldonTranspose& trans, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nl = mesh_num.GetNbDof();
    Matrix<T> mat_elt_unsym;
    int Nvol = var_problem.GetNbMainUnknownDof();

    Vector<T> Fu, Fl;
    Vector<TinyVector<T, Dimension::dim_N> > Fv;
    Vector<T> Fv_, Ft;

    Real_wp coef_ipp; T coef_u, coef_v, coef_u_damp; Complexe m_iomega;
    var_helm.GetCoefficientIPP(var_problem, m_iomega,
                               coef_ipp, coef_u, coef_v, coef_u_damp);

    if (this->UseSymmetrization() && var_problem.FirstOrderFormulation())
      {
	coef_ipp = 1.0;
	coef_v = -coef_v;
      }

    // we need the inverse of coef_v
    coef_v = Real_wp(1) / coef_v;
      
    this->level_volume.SetLevel();
    int nb_elt_lvl = this->level_volume.GetNbElt();
    for (int i0 = 0; i0 < nb_elt_lvl; i0++)
      {
	int iquad = this->level_volume.GetElementNumber(i0);
	if (!optim_elt(i0).flag)
	  {
	    if (!this->store_block)
	      {
		cout << "Not implemented" << endl;
		abort();
	      }

	    if (this->symmetric_elem_matrix)
	      CondensationBlockSolver_Fem<T>::
		ModifyRhsGen(trans, b_source, iquad, this->block_sym(i0).num_ddl, mat_elt_unsym,
			     this->block_sym(i0).a12, this->block_sym(i0).a21, this->block_sym(i0).inv_a22);
	    else
	      CondensationBlockSolver_Fem<T>::
		ModifyRhsGen(trans, b_source, iquad, this->block_unsym(i0).num_ddl, mat_elt_unsym,
			     this->block_unsym(i0).a12, this->block_unsym(i0).a21, this->block_unsym(i0).inv_a22);
	    
	    continue;
	  }

	const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(iquad);
	bool affine = var_problem.mesh.IsElementAffine(iquad);
	
	int nb_dof_elt = Fb.GetNbDof();
	int offset = mesh_num.Element(iquad).GetNbDof();
	int nb_dof_boundary = Fb.GetNbDofBoundaries();
	int ref_domain = var_problem.mesh.Element(iquad).GetReference();
	int offset_u = Nl + var_problem.GetOffsetDofV(iquad);
	    
	// etapes de modification du second membre
	Fu.Reallocate(nb_dof_elt), Fl.Reallocate(offset);
	Fv.Reallocate(nb_dof_elt);
	Fv_.Reallocate(Dimension::dim_N*nb_dof_elt);
	TinyVector<T, Dimension::dim_N> vec_u, vec_v;
	
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    Fu(i) = b_source(offset_u + i);
	    for (int j = 0; j < Dimension::dim_N; j++)
	      Fv(i)(j) = b_source(offset_u + (j+1)*Nvol + i);

	    // on multiplie par DF_i^{-1} mu Fv
	    var_helm.ref_mu(ref_domain).MltMatrix(var_problem, iquad, i, Fv(i), vec_v);
	    T invPoids = coef_v / (Fb.WeightsND(i)*nat_mat.GetCoefMass());
	    if (affine)
	      {
		vec_v *= 1.0/var_problem.Glob_jacobian(iquad)(0);
		Mlt(var_problem.Glob_DFjm1(iquad)(0), vec_v, vec_u);
	      }
	    else
	      {
		vec_v *= Fb.WeightsND(i)/var_problem.Glob_jacobian(iquad)(i);
		Mlt(var_problem.Glob_DFjm1(iquad)(i), vec_v, vec_u);		
	      }
	    
	    for (int j = 0; j < Dimension::dim_N; j++)
	      b_source(offset_u + (j+1)*Nvol + i) = vec_v(j)*invPoids;
	    
	    Fv(i) = nat_mat.GetCoefStiffness() / nat_mat.GetCoefMass() * coef_v * vec_u;
	    CopyVector(Fv(i), i, Fv_);
	  }
	
	for (int i = 0; i < offset; i++)
	  Fl(i) = b_source(mesh_num.Element(iquad).GetNumberDof(i));
	
	Ft.Reallocate(nb_dof_elt);
	Fb.ApplyRh(Fv_, Ft);
	
	const Vector<TinyVector<T, Dimension::dim_N> >& Sn = this->GetSn(i0);
	const Vector<TinyVector<T, Dimension::dim_N> >& SnD = this->GetSnD(i0);
	const Vector<int>& NumU = this->GetNumU(i0);

	Ft = Fu - Ft;
	for (int i = 0; i < nb_dof_boundary; i++)
	  Ft(i) += DotProd(Sn(i), Fv(i));
	
	Mlt(this->GetBtildeMatrix(i0), Ft, Fu);
	
	for (int i = 0; i < nb_dof_elt; i++)
	  b_source(offset_u + i) = Ft(i);
	
	Fv_.Resize(nb_dof_boundary*Dimension::dim_N);
	Fb.ApplyRhBoundaryTranspose(Fu, Fv_);

	const Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > & invBh
	  = this->GetInverseBh(i0);
	
	for (int i = 0; i < nb_dof_boundary; i++)
	  {
	    CopyVector(Fv_, i, vec_u);
	    vec_u -= Sn(i)*Fu(i);
	    Mlt(invBh(i), vec_u, vec_v);
	    if (coef_ipp == Real_wp(1))
	      vec_v = -vec_v;
	    
	    vec_v += Fv(i);
	    
	    CopyVector(vec_v, i, Fv_);
	  }

	const Vector<T>& Cl = this->GetCl(i0);
	for (int i = 0; i < offset; i++)
	  {
	    CopyVector(Fv_, NumU(i), vec_u);
	    Fl(i) += coef_ipp*(-DotProd(SnD(i), vec_u) + Cl(i)*Fu(NumU(i)));
	  }
	
	for (int i = 0; i < offset; i++)
	  b_source(mesh_num.Element(iquad).GetNumberDof(i)) = Fl(i);
      }
  }


  template<>
  void CondensationBlockSolver_HelmDG<Real_wp, Complex_wp, Dimension2>
  ::RecomposeSolution(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  template<>
  void CondensationBlockSolver_HelmDG<Real_wp, Complex_wp, Dimension3>
  ::RecomposeSolution(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }
#endif


  template<class T, class Complexe, class Dimension>
  void CondensationBlockSolver_HelmDG<T, Complexe, Dimension>
  ::RecomposeSolution(const SeldonTranspose& trans, Vector<T>& x_sol, const GlobalGenericMatrix<T>& nat_mat) const
  {
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nl = mesh_num.GetNbDof();
    Matrix<T> mat_elt_unsym;
    int Nvol = var_problem.GetNbMainUnknownDof();

    Vector<TinyVector<T, Dimension::dim_N> > SnD_Lambda;
    Vector<T> Lambda, Fu, Fv_;
    Vector<T> BhSnL, Uh;
    
    Complexe m_iomega; Real_wp coef_ipp; T coef_u, coef_v, coef_u_damp;
    var_helm.GetCoefficientIPP(var_problem, m_iomega,
			       coef_ipp, coef_u, coef_v, coef_u_damp);
    
    if (this->UseSymmetrization() && var_problem.FirstOrderFormulation())
      {
	coef_ipp = 1.0;
	coef_v = -coef_v;
      }

    coef_v = nat_mat.GetCoefStiffness() / (nat_mat.GetCoefMass() * coef_v);
    
    this->level_volume.SetLevel();
    int nb_elt = this->level_volume.GetNbElt();
    for (int i0 = 0; i0 < nb_elt; i0++)
      {
	int iquad = this->level_volume.GetElementNumber(i0);
	if (!optim_elt(i0).flag)
	  {
	    if (!this->store_block)
	      {
		cout << "Not implemented" << endl;
		abort();
	      }
	    
	    if (this->symmetric_elem_matrix)
	      CondensationBlockSolver_Fem<T>::
		RecomposeSolGen(trans, x_sol, iquad, this->block_sym(i0).num_ddl, mat_elt_unsym,
				this->block_sym(i0).a12, this->block_sym(i0).a21, this->block_sym(i0).inv_a22);
	    else
	      CondensationBlockSolver_Fem<T>::
		RecomposeSolGen(trans, x_sol, iquad, this->block_unsym(i0).num_ddl, mat_elt_unsym,
				this->block_unsym(i0).a12, this->block_unsym(i0).a21, this->block_unsym(i0).inv_a22);
	    
	    continue;
	  }

	const ElementReference<Dimension, 1>& Fb = var_problem.GetReferenceElementH1(iquad);
	bool affine = var_problem.mesh.IsElementAffine(iquad);
	
	int nb_dof_elt = Fb.GetNbDof();
	int offset = mesh_num.Element(iquad).GetNbDof();
	int nb_dof_boundary = Fb.GetNbDofBoundaries();
	int ref_domain = var_problem.mesh.Element(iquad).GetReference();
	int offset_u = Nl + var_problem.GetOffsetDofV(iquad);
	
	Lambda.Reallocate(offset);
	for (int i = 0; i < offset; i++)
	  Lambda(i) = x_sol(mesh_num.Element(iquad).GetNumberDof(i));

	const Vector<int>& NumU = this->GetNumU(i0);
	const Vector<T>& Cl = this->GetCl(i0);
	const Vector<TinyVector<T, Dimension::dim_N> >& Sn = this->GetSn(i0);
	const Vector<TinyVector<T, Dimension::dim_N> >& SnD = this->GetSnD(i0);

    	// etapes de reconstruction de U et V
	SnD_Lambda.Reallocate(nb_dof_boundary);
	SnD_Lambda.Zero();
	for (int i = 0; i < offset; i++)
	  SnD_Lambda(NumU(i)) += Lambda(i)*SnD(i);

	Fu.Reallocate(nb_dof_elt);
	BhSnL.Reallocate(Dimension::dim_N*nb_dof_boundary);
	Uh.Reallocate(nb_dof_elt);
	BhSnL.Zero();

	for (int i = 0; i < nb_dof_elt; i++)
	  Fu(i) = x_sol(offset_u + i);

	const Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > & invBh
	  = this->GetInverseBh(i0);

	TinyVector<T, Dimension::dim_N> vec_u, vec_v;
	for (int i = 0; i < nb_dof_boundary; i++)
	  {
	    Mlt(invBh(i), SnD_Lambda(i), vec_u);
	    Fu(i) -= coef_ipp*DotProd(Sn(i), vec_u);
	    CopyVector(vec_u, i, BhSnL);
	  }
	
	Fb.ApplyRhBoundary(BhSnL, Uh);
	if (coef_ipp == Real_wp(1))
	  Fu += Uh;
	else
	  Fu -= Uh;

	for (int i = 0; i < offset; i++)
	  Fu(NumU(i)) += Cl(i)*Lambda(i);

	Mlt(this->GetBtildeMatrix(i0), Fu, Uh);
	
	Fv_.Reallocate(nb_dof_elt*Dimension::dim_N);
	Fb.ApplyRhTranspose(Uh, Fv_);
	
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    x_sol(offset_u + i) = Uh(i);
	    CopyVector(Fv_, i, vec_u);
	    
	    if (i < nb_dof_boundary)
	      {
		vec_u += -Sn(i)*Uh(i) + SnD_Lambda(i);
	      }
	    
	    if (affine)
	      {
		MltTrans(var_problem.Glob_DFjm1(iquad)(0), vec_u, vec_v);
		vec_v *= 1.0/var_problem.Glob_jacobian(iquad)(0);
	      }
	    else
	      {
		MltTrans(var_problem.Glob_DFjm1(iquad)(i), vec_u, vec_v);
		vec_v *= Fb.WeightsND(i)/var_problem.Glob_jacobian(iquad)(i);
	      }
	    
	    var_helm.ref_mu(ref_domain).MltMatrix(var_problem, iquad, i, vec_v, vec_u);
	    
	    vec_u *= coef_v;
	    if (coef_ipp == Real_wp(1))
	      vec_u = -vec_u;
	    
	    for (int j = 0; j < Dimension::dim_N; j++)
	      x_sol(offset_u + (j+1)*Nvol + i) += vec_u(j);
	  }
      }
  }

  
    //! returns the memory used by the object in bytes
  template<class T, class Complexe, class TypeEquation>
  size_t CondensationBlockSolver_HelmDG<T, Complexe, TypeEquation>::GetMemorySize() const
  {
    size_t taille = CondensationBlockSolver_Fem<T>::GetMemorySize();
    
    taille += optim_elt.GetM()*sizeof(void*);
    for (int i = 0; i < optim_elt.GetM(); i++)
      taille += optim_elt(i).GetMemorySize();
    
    return taille;
  }

  
  /****************************
   * HelmholtzEquationDG_Base *
   ****************************/
  
  
  //! which derivatives to evaluate during matrix-vector product ?
  template<class T, class Dimension> template<class T0, class Vector1>
  void HelmholtzEquationDG_Base<T, Dimension>::
  GetNeededDerivative(const VarHelmholtz_Base<T, Dimension>& vars,
                      const GlobalGenericMatrix<T0>& nat_mat,
                      Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    unknown_to_derive.Fill(false);
    fct_test_to_derive.Fill(false);
    unknown_to_derive(0) = true;
    fct_test_to_derive(0) = true;
  }
  
  
  //! fills tensors D and E appearing in the variational formulation
  /*!
    \param[in] vars considered problem
    \param[in] num_elem element where D and E must be computed
    \param[in] jloc quadrature point where D and E must be computed
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the element
    \param[out] Ephi_grad tensor E
    \param[out] Dgrad_phi tensor D    
    The tensors D and E are appearing in the terms
    \int_K D \nabla u v + E u \nabla v dx
    of the variational formulation
   */
  template<class T, class Dimension> template<class T0, class MatStiff>
  void HelmholtzEquationDG_Base<T, Dimension>::
  GetGradPhiTensor(const VarHelmholtz_Base<T, Dimension>& vars,
                   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    for (int i = 0; i < Dimension::dim_N; i++)
      {
	Ephi_grad(i+1, 0)(i) = -nat_mat.GetCoefStiffness();
	Dgrad_phi(0, i+1)(i) = -nat_mat.GetCoefStiffness();
      }
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
  template<class T, class Dimension> template<class T0, class Vector1, class Vector2>
  void HelmholtzEquationDG_Base<T, Dimension>::
  ApplyGradientUnknown(const VarHelmholtz_Base<T, Dimension>& var,
                       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		       int ref, Vector1& Vn, Vector2& Un)
  {
    T0 zero; SetComplexZero(zero);
    Un.Fill(zero);
    for (int j = 0; j < Dimension::dim_N; j++)
      Un(j+1) = -Vn(0)(j)*nat_mat.GetCoefStiffness();
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
  template<class T, class Dimension> template<class T0, class Vector1, class Vector2>
  void HelmholtzEquationDG_Base<T, Dimension>::
  ApplyGradientFctTest(const VarHelmholtz_Base<T, Dimension>& var,
                       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		       int ref, Vector1& Un, Vector2& Vn)
  {
    FillZero(Vn);
    for (int j = 0; j < Dimension::dim_N; j++)
      Vn(0)(j) = -Un(j+1)*nat_mat.GetCoefStiffness();
  }
  

  /***********************
   * HelmholtzEquationDG *
   ***********************/

  
  //! returns the matrix M, in the integral \f$ \int_K M u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] mass matrix M
  */
  template<class Dimension> template<class TypeEquation, class MatMass>
  void HelmholtzEquationDG<Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int i, int j, const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref, MatMass& mass)
  {
    FillZero(mass);
    mass(0, 0) = -vars.GetOmega()*vars.ref_rho(ref).GetCoefficient(vars, i, j)
      *nat_mat.GetCoefMass();
    
    mass(0, 0) -= Iwp*vars.ref_sigma(ref).GetCoefficient(vars, i, j)*nat_mat.GetCoefMass();
    typename Dimension::MatrixN_Nsym_Complex_wp mu
      = vars.ref_invMu(ref).GetCoefficient(vars, i, j);
    
    Mlt(-vars.GetOmega()*nat_mat.GetCoefMass(), mu);
    for (int p = 0; p < Dimension::dim_N; p++)
      for (int q = 0; q < Dimension::dim_N; q++)
        mass(p+1, q+1) = mu(p, q);    
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
  template<class Dimension> template<class TypeEquation, class Vector1>
  void HelmholtzEquationDG<Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref, Vector1& Un, Vector1& Vn)
  {
    Vn.Fill(0);
    Vn(0) = -Un(0)*var.ref_rho(ref).GetCoefficient(var, i, j);
    
    typename Dimension::MatrixN_Nsym_Complex_wp invMu = var.ref_invMu(ref).GetCoefficient(var, i, j);
    for (int n = 0; n < Dimension::dim_N; n++)
      for (int q = 0; q < Dimension::dim_N; q++)
        Vn(n+1) -= Un(q+1)*invMu(n, q);
    
    Vn *= var.GetOmega()*nat_mat.GetCoefMass();
    Vn(0) -= Iwp*Un(0)*var.ref_sigma(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefMass();
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
  template<class Dimension>  template<class Matrix1, class TypeEquation>
  void HelmholtzEquationDG<Dimension>::
  GetNabc(Matrix1& Nabc, typename Dimension::R_N& normale,
	  int ref, int iquad, int k, const GlobalGenericMatrix<Complex_wp>& nat_mat,
	  int ref_d, const EllipticProblem<TypeEquation>& vars,
	  const ElementReference<Dimension, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    FillZero(Nabc);
    Real_wp coef = 1.0;
    if (vars.FirstOrderFormulation() && (!vars.GetSymmetrizationUse()))
      coef = -1.0;
    
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// condition Dirichlet
	// N = (0, n; -n^t, 0)
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Nabc(0, i+1) = -normale(i);
	    Nabc(i+1, 0) = -coef*normale(i);
	  }
        Nabc *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Nabc(0, i+1) = normale(i);
	    Nabc(i+1, 0) = coef*normale(i);
	  }
        Nabc *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
	Complex_wp mu = vars.ref_mu(ref_d).GetCoefficient(vars, iquad, k)(0, 0);
	// we multiply by the impedance present in the data file
	vars.MltParamCondition(ref, 0, mu);
	mu /= vars.GetOmega();

	// then we form the matrix
	if (vars.FirstOrderFormulation())
	  mu *= Iwp;
	
	Nabc(0, 0) = mu*nat_mat.GetCoefStiffness();	
	mu = nat_mat.GetCoefStiffness()/mu;
	if (vars.FirstOrderFormulation() && (!vars.GetSymmetrizationUse()))
	  mu = -mu;
	
	for (int i = 0; i < Dimension::dim_N; i++)
	  for (int j = 0; j < Dimension::dim_N; j++)
	    Nabc(i+1, j+1) = -mu*normale(i)*normale(j);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        // the velocity is equal to sqrt( n \cdot mu n / rho)
        Complex_wp rho = vars.ref_rho(ref_d).GetCoefficient(vars, iquad, k);
        Complex_wp sigma = vars.ref_sigma(ref_d).GetCoefficient(vars, iquad, k);
        rho += Iwp*sigma/vars.GetOmega();
        
        TinyVector<Complex_wp, Dimension::dim_N> mu_n;
        vars.ref_mu(ref_d).MltMatrix(vars, iquad, k, normale, mu_n);
        
        // coef = rho c where c is the velocity
        Complex_wp coef2 = sqrt(rho*DotProd(mu_n, normale));

        // first-order abc is given by the matrix
        // N = [ rho c, 0, 0; 0, nx^2/(rho c), nx ny/(rho c); 0, nx ny/(rho c), ny^2/(rho c)]
	Nabc(0, 0) = coef2; coef2 = Complex_wp(1, 0)/coef2;
	if (vars.FirstOrderFormulation() && vars.GetSymmetrizationUse())
	  coef2 = -coef2;
	
	for (int i = 0; i < Dimension::dim_N; i++)
	  for (int j = 0; j < Dimension::dim_N; j++)
	    Nabc(i+1, j+1) = coef2*normale(i)*normale(j);

	if (vars.FirstOrderFormulation())
	  Mlt(nat_mat.GetCoefStiffness(), Nabc);
	else
	  Mlt(-Iwp*nat_mat.GetCoefStiffness(), Nabc);
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
  template<class Dimension>  template<class TypeEquation, class Matrix1>
  void HelmholtzEquationDG<Dimension>
  ::GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale,
	       int iquad, int k, int num_face,
	       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref, int ref2,
	       const EllipticProblem<TypeEquation>& vars,
	       const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
    Complex_wp coef = -Iwp*nat_mat.GetCoefDamping();
    if (vars.FirstOrderFormulation())
      coef = nat_mat.GetCoefDamping();
    
    if (vars.upwind_fluxes)
      {        
        // upwind fluxes : see Hesthaven paper for Maxwell's equations
        Complex_wp Ym = vars.coefficient_impedance_absorbing(ref);
        Complex_wp Yp = vars.coefficient_impedance_absorbing(ref2);
        Complex_wp Zm = 1.0/Ym;
        Complex_wp Zp = 1.0/Yp;
        Complex_wp inv_Zbar = 1.0/(Zm + Zp); 
        Complex_wp inv_Ybar = 1.0/(Ym + Yp);
        Complex_wp coef_u = (Zp - Zm)*inv_Zbar*nat_mat.GetCoefStiffness();
        Complex_wp coef_v = (Yp - Ym)*inv_Ybar*nat_mat.GetCoefStiffness();

	// multiplication by two because of PoidsFlux which is divided by two
	coef *= Real_wp(2);

	if (vars.FirstOrderFormulation())
	  coef_u = -coef_u;
	
	Nabc(0, 0) = -inv_Zbar*coef;
	if (vars.FirstOrderFormulation() && vars.GetSymmetrizationUse())
	  {
	    coef = -coef;
	    coef_v = -coef_v;
	  }
	
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
            Nabc(0, i+1) = coef_u*normale(i);
            Nabc(i+1, 0) = -coef_v*normale(i);
            for (int j = 0; j < Dimension::dim_N; j++)
              Nabc(i+1, j+1) = -inv_Ybar*coef*normale(i)*normale(j);
          }

        return;
      }
    
    Nabc(0, 0) = coef*vars.alpha_penalization;
    if (vars.FirstOrderFormulation() && vars.GetSymmetrizationUse())
      coef = -coef;
    
    for (int i = 0; i < Dimension::dim_N; i++)
      for (int j = 0; j < Dimension::dim_N; j++)
	Nabc(i+1, j+1) = coef*vars.delta_penalization*normale(i)*normale(j);
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
  template<class Dimension>  template<class Vector1, class Vector2, class TypeEquation>
  void HelmholtzEquationDG<Dimension>::
  MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
	     int i, int k, int num_face, const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref, int ref2,
             const EllipticProblem<TypeEquation>& vars,
	     const ElementReference<Dimension, 1>& Fb)
  {
    // int ref_elem = vars.mesh.elements(i).GetReference();
    // condition absorbante sur les sauts
    Un(0) = Vn(0)*vars.alpha_penalization;
    Complex_wp vloc(0);
    for (int i = 0; i < Dimension::dim_N; i++)
      vloc += Vn(i+1)*normale(i);
    
    vloc *= vars.delta_penalization;
    for (int i = 0; i < Dimension::dim_N; i++)
      Un(i+1) = vloc*normale(i);
    
    if (vars.FirstOrderFormulation())
      Un *= nat_mat.GetCoefDamping();
    else
      Un *= -Iwp*nat_mat.GetCoefDamping();
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
  template<class Dimension>  template<class Vector1, class TypeEquation>
  void HelmholtzEquationDG<Dimension>::
  MltNabc(typename Dimension::R_N& normale, int ref,
          const Vector1& Vn, Vector1& Un, int num_elem1, int npoint,
          const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
          const EllipticProblem<TypeEquation>& vars,
	  const ElementReference<Dimension, 1>& Fb)
  {
    typename Vector1::value_type zero; SetComplexZero(zero);
    // typedef typename TypeElement::Dimension Dimension;
    // typedef typename Dimension::R_N R_N;
    int cond = vars.mesh.GetBoundaryCondition(ref); 
    Un.Fill(zero);
    Real_wp coef = 1.0;
    if (vars.FirstOrderFormulation())
      coef = -1.0;

    // DISP(cond); DISP(GlobalVariables::LINE_NEUMANN); exit(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// condition Dirichlet
	// N = (0, n; -n^t, 0)
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Un(0) += coef*Vn(i+1)*normale(i);
	    Un(i+1) = Vn(0)*normale(i);
	  }
        Un *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Un(0) -= coef*Vn(i+1)*normale(i);
	    Un(i+1) = -Vn(0)*normale(i);
	  }
        Un *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	// condition absorbante
	Un(0) = Vn(0)*vars.coefficient_impedance_absorbing(ref_d);;
	Complex_wp vloc(0);
	for (int i = 0; i < Dimension::dim_N; i++)
	  vloc += Vn(i+1)*normale(i);
	
        vloc /= vars.coefficient_impedance_absorbing(ref_d);
	for (int i = 0; i < Dimension::dim_N; i++)
	  Un(i+1) = vloc*normale(i);
	
        if (vars.FirstOrderFormulation())
          Mlt(nat_mat.GetCoefStiffness(), Un);
        else
          Mlt(-Iwp*nat_mat.GetCoefStiffness(), Un);
      }
    
    // DISP(Un); DISP(Vn);
  }


  /**************************
   * VarComputationRCS_Helm *
   **************************/
  

  template<class Dimension>
  void VarComputationRCS_Helm<Dimension>::ComputeRCS(const VectReal_wp& U0)
  {
    cout << "not possible" << endl;
    abort();
  }

  
  //! computation of radar cross section
  /*!
    \param[in] U0 solution vector
    The radar cross section is computed and written on output files 
   */
  template<class Dimension>
  void VarComputationRCS_Helm<Dimension>::ComputeRCS(const VectComplex_wp& U0)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    if ((rank_proc == 0)||(this->print_level >= 10))
      if (this->print_level >= 4)
	cout<<"Computation of RCS "<<endl;
    
    int nbPointsRCS = this->nb_angles_RCS;
    
    this->InitComputationRCS(true);
        
    VectComplex_wp trace_En;
    VectComplex_wp trace_Hn;
    var_problem.ComputeEnHnOnBoundary(this->var_mesh, U0, trace_En, trace_Hn);
    
    /*ofstream file_En("En.dat"), file_Hn("Hn.dat");
    file_En.precision(15);
    file_Hn.precision(15);
    for (int i = 0; i < trace_En.GetM(); i++)
      {
        file_En << Seldon::real(trace_En(i)) << " " << imag(trace_En(i)) << endl;
        file_Hn << Seldon::real(trace_Hn(i)) << " " << imag(trace_Hn(i)) << endl;
      }
    
    file_En.close();
    file_Hn.close();
    */
    
    R_N ur;
    Real_wp teta;
    Real_wp step_angle = (this->last_angle_RCS-this->first_angle_RCS)/(nbPointsRCS-1);
    Complex_wp rcs_teta;
    int nb_points_effective = nbPointsRCS;
    int num_angle0 = 0;
    GetParallelDistributionPoints(nb_proc, rank_proc, nbPointsRCS, nb_points_effective, num_angle0);
    
    VectComplex_wp FarField(nb_points_effective);
    VectReal_wp RCSField(nb_points_effective);
    VectReal_wp TetaField(nb_points_effective);
    
    // loop over angles
    for (int num_angle = num_angle0; num_angle < num_angle0+nb_points_effective; num_angle++)
      {
	rcs_teta = 0;
	teta = this->first_angle_RCS + step_angle*num_angle;
	ur(0) = cos(teta); ur(1) = sin(teta);
	
	for (int k = 0; k < this->var_mesh.GetNbAllQuadraturePoints(); k++)
          {
            Complex_wp arg = -Complex_wp(0, DotProd(ur, this->var_mesh.GetQuadraturePoint(k)))
              *var_problem.GetOmega();
            
            rcs_teta += this->var_mesh.GetQuadratureWeight(k)*exp(arg)
              *(Iwp*var_problem.GetOmega()*trace_En(k)
                *DotProd(ur, this->var_mesh.GetQuadratureNormale(k)) + trace_Hn(k));
          }
     	
	rcs_teta *= Complex_wp(1,1)/(4.0*sqrt(pi_wp*var_problem.GetOmega()));
        
        int m = num_angle - num_angle0;
	TetaField(m) = 180.0*teta/pi_wp;
        FarField(m) = rcs_teta;
	RCSField(m) = 10.0*log(2.0*pi_wp*abs(rcs_teta*conj(rcs_teta)))/log(10.0);
      }
    
    ofstream file_out;
    if (rank_proc == 0)
      {
        file_out.open(this->file_RCS.data());
        file_out.precision(cout.precision());
      }
    
    VectReal_wp AllRCS, AllTeta;

    if (nb_proc > 1)
      {

#ifdef SELDON_WITH_MPI
	Vector<int64_t> xtmp;
	IVect NbPointsPerProc, OffsetPointsProc;
	// nmax is the maximum number of points (among all proc)
	int nmax = nb_points_effective;
	NbPointsPerProc.Reallocate(nb_proc);
	OffsetPointsProc.Reallocate(nb_proc);
	
	MPI_Gather(&nb_points_effective, 1, MPI_INTEGER,
                   NbPointsPerProc.GetData(), 1, MPI_INTEGER, 0, var_problem.comm_group_mode);
	
	MPI_Gather(&num_angle0, 1, MPI_INTEGER,
                   OffsetPointsProc.GetData(), 1, MPI_INTEGER, 0, var_problem.comm_group_mode);
	
	nmax = NbPointsPerProc.GetNormInf();
	
	MPI_Bcast(&nmax, 1, MPI_INTEGER, 0, var_problem.comm_group_mode);
	
        if (nmax != nb_points_effective)
          {
            RCSField.Resize(nmax);
            TetaField.Resize(nmax);
          }
	
	int Nall = nmax*nb_proc;
        AllRCS.Reallocate(Nall);
        AllTeta.Reallocate(Nall);
        
	MpiGather(var_problem.comm_group_mode, RCSField, xtmp, AllRCS, nmax, 0);
	MpiGather(var_problem.comm_group_mode, TetaField, xtmp, AllTeta, nmax, 0);
	
        if (rank_proc == 0)
	  for (int p = 0; p < nb_proc; p++)
	    for (int k = 0; k < NbPointsPerProc(p); k++)
	      {
		int num_angle = OffsetPointsProc(p) + k;
		file_out<<AllTeta(num_angle)<<"  "<<AllRCS(num_angle)<<'\n';
	      }
#endif
        
      }
    else
      {
        for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
          file_out<<TetaField(num_angle)<<"  "<<RCSField(num_angle)<<endl;
      }
    
    if (rank_proc == 0)
      file_out.close();
        
    if ((rank_proc == 0)||(this->print_level >= 10))
      if (this->print_level >= 4)
	cout<<"Radar cross section OK"<<endl;
  } 


  template<class Dimension>
  void VarComputationRCS_Helm<Dimension>
  ::ComputeIntegralRepresentation(const VectReal_wp& trace_En, const VectReal_wp& trace_Hn,
                                  const MeshInterpolationFEM<Dimension>& mesh, const R_N& pointX,
                                  VectReal_wp& scal_u) const
  {
    cout << "not possible" << endl;
    abort();
  }
    

#ifdef MONTJOIE_WITH_TWO_DIM
  //! computation of u and du/dn with integral representation
  /*!
    This methods does the computation 
    \f$ u(x) = \int_\Gamma \frac{\partial \phi}{\partial n(y)}(x,y) u(y)
                         - \phi(x,y) \frac{\partial u}{\partial n(y)}(y)  dy , \quad x \in \Sigma 
    \ \frac{\partial u}{\partial n(x)}(x)  =  \int_\Gamma \frac{\partial^2 \phi}{\partial n(x)
             \partial n(y)} (x,y) u(y) - 
    \frac{\partial \phi}{n(x)}(x,y) \frac{\partial u}{\partial n(y)}(y)  dy ,\quad x \in Sigma \f$
    
    in 2D with \f$ \phi(x,y) = \frac{i}{4} H_0^{(1)}( k |x-y|) \f$
    
    \param[in] mesh_ description of \f$ \Gamma \f$ boundary 
    \param[in] pointX the point were u(x) and du/dn(x) will be computed
    \param[in] normaleX unit normale at this point
    \param[in] trace_En trace of u on Gamma
    \param[in] trace_Hn trace of du/dn on Gamma
    
    \param[out] scal_u u(x)
    \param[out] scal_du_dn \f$ \frac{\partial u(x)}{\partial n(x)} \f$
  */
  template<>
  void TransparencySolver_Helm<Dimension2>
  ::ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
				  const MeshInterpolationFEM<Dimension2>& mesh_,
				  const R2& pointX, const R2& normaleX,
				  TinyVector<Complex_wp, 1>& scal_u,
				  TinyVector<Complex_wp, 1>& scal_du_dn) const
  {
    
    Real_wp XminusY_dnX, XminusY_dnY, T;
    Complex_wp dudn, uj, H0_T, H1_T, h1st;
    R2 pointY, normaleY, XminusY;
    
    // initialization to zero
    FillZero(scal_u); FillZero(scal_du_dn);
    
    Real_wp omega = var_problem.GetOmega();
    Real_wp k_inf = omega*sqrt(abs(var_helm.rho0/var_helm.mu0));

    Real_wp k2 = k_inf*k_inf, poids;
    // loop on quadrature points
    //DISP(mesh_.GetNbAllQuadraturePoints());
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // we write n(y), y, du/dn(y) and u on local variables
        // we take opposite of n, because trace_normale(ind) is the inward normal
        // and the representation formula is written for an outward normal
        // we take also opposite of trace_dU_dN(ind) for the same reason
        normaleY(0) = -mesh_.GetQuadratureNormale(k)(0);
        normaleY(1) = -mesh_.GetQuadratureNormale(k)(1);
        pointY =  mesh_.GetQuadraturePoint(k);
        poids = mesh_.GetQuadratureWeight(k);
        dudn = -trace_Hn(k);
        uj = trace_En(k);

        // XminusY = x-y
        XminusY = pointX - pointY;
        // T = |x-y|
        T = Norm2(XminusY);
        // XminusY_dnX = (x-y) . n(x) 
        // XminusY_dnY = (x-y) . n(y) 
        XminusY_dnX = DotProd(XminusY, normaleX);
        XminusY_dnY = DotProd(XminusY, normaleY);
	
        // arg = k |x-y|
        Real_wp arg = k_inf*T;
        
        // H0_T = H_0^{(1)} (k|x-y|)
        // H1_T = H_1^{(1)} (k|x-y|)
        H0_T = Complex_wp(jn(0, arg), yn(0, arg));
        H1_T = Complex_wp(jn(1, arg), yn(1, arg));
        
        // invDist = 1 / |x-y|  h1_st = H_1^{(1)} (k|x-y|) / |x-y|
        Real_wp invDist = 1.0/T; h1st = H1_T*invDist;
        
        // XminusY_dnX = (x-y) . n(x) / |x-y|
        // XminusY_dnY = (x-y) . n(y) / |x-y|
        XminusY_dnX *= invDist; XminusY_dnY *= invDist;
	
        // u = \int_\Gamma d phi(x,y) / dn(y) u(y) - phi(x,y) du(y)/dn(y) dx
        // phi(x,y) = i/4 H_0^{(1)}( k |x-y|)  
        // so we multiply scal_u and scal_du_dn by i/4 at the end of the loop
        // d phi(x,y) / dn(y) = ik/4 H_1^{(1)} (k |x-y|) (x-y). n(y) / |x-y| 
	
        scal_u(0) += ( k_inf * H1_T * XminusY_dnY * uj  -  H0_T * dudn ) *poids;

	
        // du / dn(x)  =  \int_\Gamma d^2 phi(x,y) / dn(x) dn(y) u(y) 
        //                           - dphi(x,y) / dn(x) du(y)/dn(y)  dy 
        // d phi(x,y) / dn(x) = -ik/4 H_1^{(1)} (k |x-y|) (x-y). n(y) / |x-y| 
        // d^2 phi(x,y) / dn(x) dn(y) = i/4 * k^2 H_0^{(1)} (k |x-y|)
        //              (x-y) . n(x) (x-y) . n(y) / |x-y|^2
        //              - i/4 * 2 k H_1^{(1)} (k |x-y|) (x-y) . n(x) (x-y) . n(y) / |x-y|^3
        //                              + i/4 * k H_1^{(1)} (k |x-y|) n(x) . n(y) / |x-y|
	
        scal_du_dn(0) += ( ( XminusY_dnX * XminusY_dnY * ( k2 * H0_T - 2.0 * k_inf * h1st ) + \
                             k_inf * h1st * DotProd(normaleX,normaleY ) ) * uj \
                           + k_inf * H1_T * XminusY_dnX                   * dudn ) * poids;

      }
    
    scal_u(0) *= Real_wp(0.25)*Iwp; scal_du_dn(0) *=  Real_wp(0.25)*Iwp;
    
  } // end method ComputeIntegralRepresentation (2D case)


  //! computation of u with integral representation
  /*!
    This methods does the computation 
    \f$ u(x) = \int_\Gamma \frac{\partial \phi}{\partial n(y)}(x,y) u(y) 
     - \phi(x,y) \frac{\partial u}{\partial n(y)}(y)  dy , \quad x \in \Sigma \f$
    
    in 2D with \f$ \phi(x,y) = \frac{i}{4} H_0^{(1)}( k |x-y|) \f$
    
    \param[in] mesh_ description of \f$ \Gamma \f$ boundary 
    \param[in] pointX the point were u(x) and du/dn(x) will be computed
    \param[in] trace_En trace of u on Gamma
    \param[in] trace_Hn trace of du/dn on Gamma
    
    \param[out] scal_u u(x)
  */
  template<>
  void VarComputationRCS_Helm<Dimension2>::
  ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
				const MeshInterpolationFEM<Dimension2>& mesh_,
				const R2& pointX, VectComplex_wp& scal_u) const
  {
    Real_wp XminusY_dnY, T;
    Complex_wp dudn, uj, H0_T, H1_T, h1st;
    R2 pointY, normaleY, XminusY;
    
    // initialization to zero
    scal_u(0) = 0;
    
    Real_wp omega = var_problem.GetOmega();
    Real_wp mu0(1), rho0(1);
    var_source.GetCoefAB_Infinity(mu0, rho0);
    Real_wp k_inf = omega*sqrt(rho0/mu0);
    
    Real_wp poids;
    // loop on quadrature points
    //DISP(mesh_.GetNbAllQuadraturePoints());
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // we write n(y), y, du/dn(y) and u on local variables
        // we take opposite of n, because trace_normale(ind) is the inward normal
        // and the representation formula is written for an outward normal
        // we take also opposite of trace_dU_dN(ind) for the same reason
        normaleY(0) = -mesh_.GetQuadratureNormale(k)(0);
        normaleY(1) = -mesh_.GetQuadratureNormale(k)(1);
        pointY =  mesh_.GetQuadraturePoint(k);
        poids = mesh_.GetQuadratureWeight(k);
        dudn = -trace_Hn(k);
        uj = trace_En(k);
                
        // XminusY = x-y
        XminusY = pointX - pointY;
        // T = |x-y|
        T = Norm2(XminusY);
        // XminusY_dnY = (x-y) . n(y) 
        XminusY_dnY = DotProd(XminusY, normaleY);
	
        // arg = k |x-y|
        Real_wp arg = k_inf*T;
        
        // H0_T = H_0^{(1)} (k|x-y|)
        // H1_T = H_1^{(1)} (k|x-y|)
        H0_T = Complex_wp(jn(0, arg), yn(0, arg));
        H1_T = Complex_wp(jn(1, arg), yn(1, arg));
        
        // invDist = 1 / |x-y|  h1_st = H_1^{(1)} (k|x-y|) / |x-y|
        Real_wp invDist = 1.0/T; h1st = H1_T*invDist;
        
        // XminusY_dnY = (x-y) . n(y) / |x-y|
        XminusY_dnY *= invDist;
	
        // u = \int_\Gamma d phi(x,y) / dn(y) u(y) - phi(x,y) du(y)/dn(y) dx
        // phi(x,y) = i/4 H_0^{(1)}( k |x-y|)  
        // so we multiply scal_u and scal_du_dn by i/4 at the end of the loop
        // d phi(x,y) / dn(y) = ik/4 H_1^{(1)} (k |x-y|) (x-y). n(y) / |x-y| 
	
        scal_u(0) += ( k_inf * H1_T * XminusY_dnY * uj  -  H0_T * dudn ) *poids;
      }
    
    scal_u(0) *= Real_wp(0.25)*Iwp;
    
  } // end method ComputeIntegralRepresentation (2D case)
  
#endif
  
#ifdef MONTJOIE_WITH_THREE_DIM
  
  //! computation of u and du/dn from values defined on a closed surface, for Helmholtz equation
  /*!
    The representation formula used reads as
    \f$ u(x) = \int_\Gamma \frac{\partial \phi(x,y)}{\partial n(y)} u(y)
    - \phi(x,y) \frac{du(y)}{dn(y)}  dy  , x \in \Sigma \f$
    \f$ \frac{du}{dn(x)}  =  \int_\Gamma \frac{\partial^2 \phi(x,y) }{\partial n(x)
    \partial n(y)} u(y)
    - \frac{\partial \phi(x,y)}{\partial n(x)} \frac{du(y)}{dn(y)}  dy \; x \in \Sigma \f$
    in 3D with \f$ \phi(x,y) = \frac{exp( i k |x-y|)}{4 pi |x-y|} \f$
  
    INPUT
  
    trace_position : the positions y on Gamma
    trace_normale  : the normale vectors n(y) on Gamma
    trace_U        : values u(y) on Gamma
    trace_dU_dN    : values du(y) / dn(y) on Gamma
    gauss_int      : integration used to do the computation
    pointX         : the point were u(x) and du/dn(x) will be computed
    normaleX       : normale vector at this point
    info_section(0) : boundary condition associated with Gamma 
  
    OUTPUT
  
    scal_u     : u(x)
    scal_du_dn : du(x)/dn(x)
  */
  template<>
  void TransparencySolver_Helm<Dimension3>
  ::ComputeIntegralRepresentation(const VectComplex_wp& trace_En, const VectComplex_wp& trace_Hn,
				  const MeshInterpolationFEM<Dimension3>& mesh_,
				  const R3& pointX, const R3& normaleX,
				  TinyVector<Complex_wp, 1>  & scal_u,
				  TinyVector<Complex_wp, 1>& scal_du_dn) const
  {
    R3 pointY, normaleY;
    Complex_wp uj, du_dn, phi;
    R3_Complex_wp grad_phi, tmp_phi;
    Matrix3_3sym_Complex_wp hessian_phi;
    Real_wp poids;

    // initialization to zero
    scal_u(0) = 0.0; scal_du_dn(0) = 0.0;

    Real_wp omega = var_problem.GetOmega();
    Real_wp k_inf = omega*sqrt(abs(var_helm.rho0/var_helm.mu0));

    // loop on quadrature points
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // we write n(y), y, du/dn(y) and u on local variables
        // we take opposite of n, because trace_normale(ind) is the inward normal
        // and the representation formula is written for an outward normal
        // we take also opposite of trace_dU_dN(ind) for the same reason
        normaleY = mesh_.GetQuadratureNormale(k);
        Mlt(-1.0, normaleY);
        
        pointY =  mesh_.GetQuadraturePoint(k);
        poids = mesh_.GetQuadratureWeight(k);
        du_dn = -trace_Hn(k);
        uj = trace_En(k);
        
        this->ComputeGreenKernel(pointX, pointY, k_inf, phi, grad_phi, hessian_phi);
        
        scal_u(0) += (DotProd(grad_phi, normaleY)*uj - phi * du_dn)*poids;
	
        Mlt(hessian_phi, normaleX, tmp_phi);
        
        scal_du_dn(0) += (-DotProd(tmp_phi, normaleY)*uj
                          + DotProd(grad_phi, normaleX) * du_dn ) * poids;
      }
    
  } // end method ComputeIntegralRepresentation (3-D case)


  //! computation of u with integral representation
  /*!
    This methods does the computation 
    \f$ u(x) = \int_\Gamma \frac{\partial \phi}{\partial n(y)}(x,y) u(y)
    - \phi(x,y) \frac{\partial u}{\partial n(y)}(y)  dy , \quad x \in \Sigma \f$
    
    in 3D with \f$ \phi(x,y) = \frac{exp( i k |x-y|)}{4 \pi |x-y|} \f$
    
    \param[in] mesh_ description of \f$ \Gamma \f$ boundary 
    \param[in] pointX the point were u(x) and du/dn(x) will be computed
    \param[in] trace_En trace of u on Gamma
    \param[in] trace_Hn trace of du/dn on Gamma
    
    \param[out] scal_u u(x)
  */
  template<>
  void VarComputationRCS_Helm<Dimension3>::
  ComputeIntegralRepresentation(const VectComplex_wp& trace_En,	const VectComplex_wp& trace_Hn,
				const MeshInterpolationFEM<Dimension3>& mesh_, const R3& pointX,
				VectComplex_wp& scal_u) const
  {
    R3 pointY, normaleY;
    Complex_wp uj, du_dn, phi;
    R3_Complex_wp grad_phi, tmp_phi;
    Real_wp poids;
    R3 xMinusy; Real_wp T, inv_T; Complex_wp arg;
    Complex_wp phi_div_T, alpha;
    
    // initialization to zero
    scal_u(0) = 0.0;

    Real_wp omega = var_problem.GetOmega();
    Real_wp mu0(1), rho0(1);
    var_source.GetCoefAB_Infinity(mu0, rho0);
    Real_wp k_inf = omega*sqrt(rho0/mu0);

    // loop on quadrature points
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // we write n(y), y, du/dn(y) and u on local variables
        // we take opposite of n, because trace_normale(ind) is the inward normal
        // and the representation formula is written for an outward normal
        // we take also opposite of trace_dU_dN(ind) for the same reason
        normaleY = mesh_.GetQuadratureNormale(k);
        Mlt(-1.0, normaleY);
        
        pointY =  mesh_.GetQuadraturePoint(k);
        poids = mesh_.GetQuadratureWeight(k);
        du_dn = -trace_Hn(k);
        uj = trace_En(k);
        
        // in 3-D \phi(x,y) = exp(ik|x-y|) / (4 pi |x-y| )
        xMinusy = pointX - pointY;
        // T = |x-y|   T2 = |x-y|^2   arg = ik |x-y|
        T = Norm2(xMinusy); arg = Iwp*k_inf*T;
        inv_T = 1.0/T;
        
        // evaluation of green function
        phi = exp(arg) * inv_T * 0.25 / pi_wp;
        phi_div_T = phi * inv_T;
        
        // evaluation of its gradient
        grad_phi(0) = xMinusy(0);
        grad_phi(1) = xMinusy(1);
        grad_phi(2) = xMinusy(2);
    
        alpha = (-Iwp*k_inf + inv_T)*phi_div_T;
        Mlt( alpha , grad_phi);
        
        scal_u(0) += (DotProd(grad_phi, normaleY)*uj - phi * du_dn)*poids;
      }
    
  } // end method ComputeIntegralRepresentation (3-D case)
  
#endif
  
  
  //! computation of transparent source for Helmholtz equation
  template<class Dimension>
  void TransparencySolver_Helm<Dimension>
  ::GetSource(const VectComplex_wp& trace_u, const VectComplex_wp& trace_du_dn,
	      int n, const Real_wp& k_inf, const R_N& point, const R_N& normale,
	      Vector<VectComplex_wp>& scal_g, int j) const
  {
    Complex_wp scal_u = trace_u(n);
    Complex_wp scal_du_dn = trace_du_dn(n);
    if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	// Ap = \int_\Sigma (du_tld - ik u_tld) \phi_j ds
	scal_g(0)(j) = var_helm.mu0*(scal_du_dn - Iwp*k_inf*scal_u);
      }
    else
      {
	scal_g(0)(j) = 0.5*(scal_du_dn/var_problem.GetOmega()  - Iwp*scal_u);
	Complex_wp fpot = -0.5*(scal_u + Iwp/var_problem.GetOmega()*scal_du_dn);
        for (int i = 0; i < Dimension::dim_N; i++)
	  scal_g(i+1)(j) = fpot*normale(i);	
      }
  }  
  
}

#define MONTJOIE_FILE_VAR_HELMHOLTZ_CXX
#endif
