#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HCURL_3D_CXX

namespace Montjoie
{  
  
  //! generic matrix-vector product for edge finite elements in 3-D
  //! and Maxwell equations
  /*
    this algorithm separates the geometry \f$ DF_i^{*-1} \f$ from the 
    spatial derivatives. This allows low-storage and fast-algorithm
    in \f$ O(r^4) \f$ instead of \f$ O(r^6) \f$ if the full matrix is stored, where
    r is the order of approximation
  */
  template <class T, class T0, class TypeEquation>
  void MltAddHcurl3D(const T& alpha, const GlobalGenericMatrix<T0>& nat_mat,
                     const SeldonTranspose& trans, int level,
                     const FemMatrixFreeClass<T0, TypeEquation>& A,
                     const Vector<T>& B2, const T& beta, Vector<T>& C2, bool assemble)
  {
#ifdef SELDON_WITH_MPI
    Vector<T> B(B2), C(C2);
#else
    Vector<T> B, C;
#endif
    A.ApplyRightScaling(B2, C2, B, C);    

    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (beta == zero)
      C.Zero();
    else if (beta != one)
      Mlt(beta, C);

    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    typedef typename TypeEquation::Complexe Complexe;
    
    Complexe m_iomega; var.GetMiomega(m_iomega);
    Complexe m_omega2 = m_iomega*m_iomega;
    if (A.var.FirstOrderFormulation())
      {
	int Nvol = mesh_num.GetNbDof();
	int Nscal = var.offset_dof_unknown(1);

	T m = nat_mat.GetCoefMass()*alpha;
	T s = nat_mat.GetCoefStiffness()*alpha;
	T sig = nat_mat.GetCoefDamping()*alpha;
	
	bool sym = A.GetSymmetrizationUse();
	bool stiff = false, mass = true;
	if (s != zero)
	  stiff = true;
	
	if (m != zero)
	  mass = true;
	
	if (sig != zero)
	  mass = true;

	T coef_H = m_iomega*m;
	bool signE = true, signH = true;
	if (sym)
	  coef_H = -coef_H;
	
	if (!sym)
	  {
	    if (trans.Trans())
	      signE = false;
	    else
	      signH = false;
	  }

        Real_wp coef_sym(1);
        if (sym)
          coef_sym = Real_wp(-1);
        
	TinyVector<Vector<T>, 1> Eloc, ProdEloc, Eloc_star, ProdEloc_star;
	Vector<T> Equad, Equad_star, curlE_quad, ProdEquad, ProdEquad_star, contrib;

	// we form E* (for PML layers)
	Vector<T> Estar, Cstar;
	if (var.GetNbGlobalEltPML() > 0)
	  {
	    Estar.Reallocate(Nvol);
	    Cstar.Reallocate(Nvol);
	    Cstar.Zero();
	    for (int i = 0; i < Nvol; i++)
	      {
		Estar(i) = B(i);
		int npml = mesh_num.GetDofPML(i);
		if (npml >= 0)
		  Estar(i) = B(Nvol + npml);
	      }
	  }

	T coef_E = sig / m_iomega;
	for (int i0 = 0; i0 < mesh.GetNbElt(); i0++)
	  {
	    //int i = list_level.GetElementNumber(i0);
	    int i = i0;
	    
	    const ElementReference<Dimension3, 2>& Fb = var.GetReferenceElementHcurl(i);
	    const HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass_elt = var.Glob_matMass_elem(i);
	    
	    int nb_dof_elt = Fb.GetNbDof();
	    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	    bool diag_mass = mass_elt.IsDiagonalMass();
	    
	    bool affine = mesh.IsElementAffine(i);
	    int ref = mesh.Element(i).GetReference();
	    bool pml_elt = var.InsidePML(i);
	    bool ortho = (var.OrthogonalElement(i) == 0);
	    int i1 = i - mesh.GetNbElt() + var.GetNbEltPML();
	    if (!pml_elt)
	      ortho = true;
	    
	    Eloc(0).Reallocate(nb_dof_elt);
	    ProdEloc(0).Reallocate(nb_dof_elt);
	    Equad.Reallocate(3*nb_points_quad);
	    ProdEquad.Reallocate(3*nb_points_quad);
	    ProdEloc(0).Zero();
	    ProdEquad.Zero();
	    
	    TinyVector<T, 3> vec_u, vec_v, vec_w, vec_Hs, vec_H, vec_Q, vec_P, vec_Qv;
            TinyVector<T, 3> vec_Es;
	    R3 tau_r; TinyVector<int, 3> permut120(1, 2, 0), permut201(2, 0, 1);
	    TinyVector<Complexe, 3> tau;
	    Complexe epsilon = var.ref_epsilon(ref).GetConstant()(0, 0);
	    T coef_Es = sig / (epsilon*m_iomega);
	    T coefM_Es = m / epsilon;
	    
	    // we check that the element does not touch a PML layer
	    bool close_pml_elt = false;
	    if ( (!pml_elt) && (var.GetNbGlobalEltPML() > 0) && !mesh_num.drop_interface_pml_dof)
	      for (int j = 0; j < Fb.GetNbDof(); j++)
		{
		  int num_dof = mesh_num.Element(i).GetNumberDof(j);
		  if (num_dof >= 0)
		    {
		      int npml = mesh_num.GetDofPML(num_dof);
		      if (npml >= 0)
			close_pml_elt = true;
		    }
		}
	    
	    // we retrieve values of E on the element
	    var.GetLocalUnknownVector(B, i, Eloc);
            
	    // and values of E*
	    if (close_pml_elt || pml_elt)
	      {
		Eloc_star(0).Reallocate(nb_dof_elt);
		var.GetLocalUnknownVector(Estar, i, Eloc_star);
                
		ProdEloc_star(0).Reallocate(nb_dof_elt);
		ProdEloc_star(0).Zero();
	      }
            
	    if (Fb.UseQuadraturePointsForRh())
	      contrib.Reallocate(nb_dof_elt);
	    
	    if (pml_elt)
	      {
		if (Fb.UseQuadraturePointsForRh() || !ortho)
		  ProdEquad_star.Reallocate(3*nb_points_quad);		
	      }
	    
	    // offsets for values of H
	    int offset_H = Nscal + var.GetOffsetDofV(i);
	    int offset_Hs = offset_H;
	    if (pml_elt)
	      offset_Hs += 3*nb_points_quad;

	    int offset_Hs_rot_inv = offset_H;
	    int offset_Hs_rot = offset_Hs;
	    if (trans.Trans())
	      {
		offset_Hs_rot = offset_H;
		offset_Hs_rot_inv = offset_Hs;
	      }
            
	    // values of E on quadrature points
	    Fb.ApplyChTranspose(Eloc(0), Equad);

	    if (close_pml_elt || pml_elt)
	      {
		Equad_star.Reallocate(3*nb_points_quad);
		Fb.ApplyChTranspose(Eloc_star(0), Equad_star);
	      }
            
            bool drude_stiff = false;

            // stiffness terms
	    if (stiff)
	      {
		// curl of E*
		curlE_quad.Reallocate(3*nb_points_quad);
		if ((close_pml_elt || pml_elt) && (trans.NoTrans()))
		  {
		    if (Fb.UseQuadraturePointsForRh())
		      Fb.ApplyRhQuadratureTranspose(Equad_star, curlE_quad);
		    else
		      Fb.ApplyRhTranspose(Eloc_star(0), curlE_quad);
		  }
		else
		  {
		    if (Fb.UseQuadraturePointsForRh())
		      Fb.ApplyRhQuadratureTranspose(Equad, curlE_quad);
		    else
		      Fb.ApplyRhTranspose(Eloc(0), curlE_quad);
		  }
		
		if (affine)
		  {
		    Matrix3_3 dfj = var.Glob_DFj(i)(0);
		    for (int k = 0; k < nb_points_quad; k++)
		      {
			T poids = Fb.WeightsND(k)*s;

			// part curl E \cdot \psi
			CopyVector(curlE_quad, k, vec_u);
			Mlt(dfj, vec_u, vec_v);
			vec_v *= poids;
			if (signH)
			  vec_v = -vec_v;

			AddVector(vec_v, offset_Hs_rot_inv+3*k, C);
			
			// part H* \cdot curl phi
			ExtractVector(B, offset_Hs_rot + 3*k, vec_u);
			MltTrans(dfj, vec_u, vec_v);
			vec_v *= poids;
			if (signE)
			  vec_v = -vec_v;

			CopyVector(vec_v, k, curlE_quad);
		      }
		  }
		else
		  {
		    for (int k = 0; k < nb_points_quad; k++)
		      {
			T poids = Fb.WeightsND(k)*s;

			// part curl E \cdot \psi
			CopyVector(curlE_quad, k, vec_u);
			Mlt(var.Glob_DFj(i)(k), vec_u, vec_v);
			vec_v *= poids;
			if (signH)
			  vec_v = -vec_v;

			AddVector(vec_v, offset_Hs_rot_inv+3*k, C);

			// part H* \cdot curl phi
			ExtractVector(B, offset_Hs_rot + 3*k, vec_u);
			MltTrans(var.Glob_DFj(i)(k), vec_u, vec_v);
			vec_v *= poids;
			if (signE)
			  vec_v = -vec_v;
			
			CopyVector(vec_v, k, curlE_quad);
		      }
		  }

		if (pml_elt && trans.Trans())
		  {
		    if (Fb.UseQuadraturePointsForRh())
		      {
			Fb.ApplyRhQuadrature(curlE_quad, ProdEquad_star);
			Fb.ApplyCh(ProdEquad_star, ProdEloc_star(0));
		      }
		    else
		      Fb.ApplyRh(curlE_quad, ProdEloc_star(0));
		  }
		else
		  {
		    if (Fb.UseQuadraturePointsForRh())
		      Fb.ApplyRhQuadrature(curlE_quad, ProdEquad);
		    else
		      Fb.ApplyRh(curlE_quad, ProdEloc(0));

		    if (close_pml_elt && trans.Trans())
		      {
			if (Fb.UseQuadraturePointsForRh())
			  {
			    Fb.ApplyCh(ProdEquad, ProdEloc(0));
			    ProdEquad.Zero();
			  }

			for (int j = 0; j < Fb.GetNbDof(); j++)
			  {
			    int num_dof = mesh_num.Element(i).GetNumberDof(j);
			    if (num_dof >= 0)
			      {
				int npml = mesh_num.GetDofPML(num_dof);
				if (npml >= 0)
				  {
				    ProdEloc_star(0)(j) += ProdEloc(0)(j);
				    SetComplexZero(ProdEloc(0)(j));
				  }
			      }
			  }
		      }
		  }
	      }
            
            // Drude's model part
            if (var.ref_drude(ref).IsEnabled())
              {
                int nb_dof_drude_vec = var.GetNbVectorialDofDrude();
                int offset_Q = var.GetNbDof() - nb_dof_drude_vec;
                int offset_P = offset_Q - nb_dof_drude_vec;
                offset_P += var.OffsetDofDrudeV(i);
                offset_Q += var.OffsetDofDrudeV(i);

                int nPole = var.ref_drude(ref).gamma.GetM();
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp coefE_q(1), coefQ_e(1), coefME_q(0), coefMQ_e(0);
                    Real_wp coef1(1), coef2(1), coefP_q(1), coefQ_p(1);
                    
                    Real_wp gamma = var.ref_drude(ref).gamma(kp);
                    Real_wp coef_omega2 = var.ref_drude(ref).omega_02(kp);
                    Real_wp coef_eps_inf = var.ref_drude(ref).eps_omega_p2(kp);
                    Real_wp coef_sig = var.ref_drude(ref).eps_sigma(kp);
                
                    if (var.use_symm_drude)
                      {
                        coef1 = Real_wp(1)/coef_eps_inf;
                        coef2 = coef_omega2 * coef1;
                        coefP_q = coef2; coefQ_p = coef2;
                      }
                    else
                      {
                        coefQ_e = coef_eps_inf;
                        coefMQ_e = coef_sig;
                        coefQ_p = coef_omega2;
                      }
                
                    if (var.ref_drude(ref).IsModeTE())
                      {
                        coefE_q = coef_sym;
                        coef1 *= coef_sym;
                        coef2 *= coef_sym;
                        coefP_q *= -coef_sym;
                        coefMQ_e = -coef_sig;
                        if (var.use_symm_drude)
                          coefQ_e = -1.0;
                        else
                          coefQ_e = -coef_eps_inf;
                      }
                    else
                      {
                        coefQ_e *= -coef_sym; coefMQ_e *= -coef_sym;
                        coefP_q = -coefP_q; coefQ_p *= coef_sym;
                      }
                    
                    if (trans.Trans())
                      {
                        Real_wp coef_tmp = coefE_q;
                        coefE_q = coefQ_e;
                        coefQ_e = coef_tmp;

                        coef_tmp = coefME_q;
                        coefME_q = coefMQ_e;
                        coefMQ_e = coef_tmp;
                        
                        coef_tmp = coefQ_p;
                        coefQ_p = coefP_q;
                        coefP_q = coef_tmp;
                      }
                    
                    for (int k = 0; k < nb_points_quad; k++)
                      {
                        vec_Q(0) = B(offset_Q + 3*k);
                        vec_Q(1) = B(offset_Q + 3*k+1);
                        vec_Q(2) = B(offset_Q + 3*k+2);
                        
                        vec_P(0) = B(offset_P + 3*k);
                        vec_P(1) = B(offset_P + 3*k+1);
                        vec_P(2) = B(offset_P + 3*k+2);
                        
                        Real_wp jacobian(1);
                        CopyVector(Equad, k, vec_u);
                        if (affine)
                          {
                            if (var.ref_drude(ref).IsModeTM())
                              {
                                Mlt(var.Glob_DFjm1(i)(0), vec_Q, vec_Qv);
                                MltTrans(var.Glob_DFjm1(i)(0), vec_u, vec_v);
                              }
                            
                            jacobian = var.Glob_jacobian(i)(0)*Fb.WeightsND(k);
                          }
                        else
                          {
                            if (var.ref_drude(ref).IsModeTM())
                              {
                                Mlt(var.Glob_DFjm1(i)(k), vec_Q, vec_Qv);
                                MltTrans(var.Glob_DFjm1(i)(k), vec_u, vec_v);
                              }
                            
                            jacobian = var.Glob_jacobian(i)(k);
                          }

                        vec_H(0) = B(offset_H + 3*k);
                        vec_H(1) = B(offset_H + 3*k+1);
                        vec_H(2) = B(offset_H + 3*k+2);
                        
                        T poids = s*Fb.WeightsND(k);
                        if (stiff)
                          {
                            if (var.ref_drude(ref).IsModeTM())
                              {
                                // part +Q in E-equation
                                poids *= coefE_q;
                                ProdEquad(3*k) += poids*vec_Qv(0);
                                ProdEquad(3*k+1) += poids*vec_Qv(1);
                                ProdEquad(3*k+2) += poids*vec_Qv(2);
                                
                                // part +E in Q-equation                                
                                poids = s*Fb.WeightsND(k)*coefQ_e;
                                C(offset_Q+3*k) += poids*vec_v(0);
                                C(offset_Q+3*k+1) += poids*vec_v(1);
                                C(offset_Q+3*k+2) += poids*vec_v(2);
                              }
                            else
                              {
                                // part +Q in H-equation
                                poids = jacobian*s*coefE_q;
                                C(offset_H+3*k) += poids*vec_Q(0);
                                C(offset_H+3*k+1) += poids*vec_Q(1);
                                C(offset_H+3*k+2) += poids*vec_Q(2);
                                
                                // part -eps_inf omega_p^2 H in Q-equation
                                poids = jacobian*s*coefQ_e;
                                C(offset_Q+3*k) += poids*vec_H(0);
                                C(offset_Q+3*k+1) += poids*vec_H(1);
                                C(offset_Q+3*k+2) += poids*vec_H(2);
                              }
                            
                            // part -(omega_0^2) / (eps_inf omega_p^2) Q in P-equation
                            poids = jacobian*s*coefP_q;
                            C(offset_P+3*k) += poids*vec_Q(0);
                            C(offset_P+3*k+1) += poids*vec_Q(1);
                            C(offset_P+3*k+2) += poids*vec_Q(2);
                            
                            // symmetric part -(omega_0^2)/ (eps_inf omega_p^2) P in Q-equation
                            // and part - gamma / (eps_inf omega_p^2) Q in Q-equation
                            poids = jacobian*s*coefQ_p;
                            T poids2 = coef_sym*jacobian*coef1*gamma*sig;
                            C(offset_Q+3*k) += poids*vec_P(0) + poids2*vec_Q(0);
                            C(offset_Q+3*k+1) += poids*vec_P(1) + poids2*vec_Q(1);
                            C(offset_Q+3*k+2) += poids*vec_P(2) + poids2*vec_Q(2);
                          }
                        
                        if (mass)
                          {
                            // part -i omega (omega_0^2) / (eps_inf omega_p^2) P in P-equation
                            poids = jacobian*m_iomega*m*coef2;
                            C(offset_P+3*k) += poids*vec_P(0);
                            C(offset_P+3*k+1) += poids*vec_P(1);
                            C(offset_P+3*k+2) += poids*vec_P(2);
                            
                            // part i omega / (eps_inf omega_p^2) Q in Q-equation
                            poids = coef_sym*jacobian*coef1*m_iomega*m;
                            C(offset_Q+3*k) += poids*vec_Q(0);
                            C(offset_Q+3*k+1) += poids*vec_Q(1);
                            C(offset_Q+3*k+2) += poids*vec_Q(2);

                            // part i omega sigma eps_inf E in Q-equation
                            // and symmetric part (for transpose)
                            if (var.ref_drude(ref).IsModeTM())
                              {
                                poids = m*m_iomega*Fb.WeightsND(k)*coefME_q;
                                ProdEquad(3*k) += poids*vec_Qv(0);
                                ProdEquad(3*k+1) += poids*vec_Qv(1);
                                ProdEquad(3*k+2) += poids*vec_Qv(2);
                                
                                poids = m*m_iomega*Fb.WeightsND(k)*coefMQ_e;
                                C(offset_Q+3*k) += poids*vec_v(0);
                                C(offset_Q+3*k+1) += poids*vec_v(1);
                                C(offset_Q+3*k+2) += poids*vec_v(2);
                              }
                            else
                              {
                                poids = jacobian*m*m_iomega*coefME_q;
                                C(offset_H+3*k) += poids*vec_Q(0);
                                C(offset_H+3*k+1) += poids*vec_Q(1);
                                C(offset_H+3*k+2) += poids*vec_Q(2);
                                
                                // part -eps_inf omega_p^2 H in Q-equation
                                poids = jacobian*m*m_iomega*coefMQ_e;
                                C(offset_Q+3*k) += poids*vec_H(0);
                                C(offset_Q+3*k+1) += poids*vec_H(1);
                                C(offset_Q+3*k+2) += poids*vec_H(2);
                              }
                          }
                      }
                    
                    offset_P += 3*nb_points_quad; offset_Q += 3*nb_points_quad;
                  }
                
                if (var.ref_drude(ref).IsModeTE())
                  drude_stiff = true;
              }
            
            if (mass && !diag_mass)
              drude_stiff = false;
            
	    if (mass)
	      {
		T poids(zero); Real_wp jacobian;
		bool variable = var.UseNumericalIntegration(i);

		for (int k = 0; k < nb_points_quad; k++)
		  {
		    // mass term for E		    
		    if (!diag_mass)
		      {
			CopyVector(Equad, k, vec_u);
			if (!variable)
			  {
			    Mlt(var.Glob_matMass_elem(i).GetBh(0), vec_u, vec_v);
			    Mlt(var.Glob_matMass_elem(i).GetBhSigma(0), vec_u, vec_w);
			    
			    vec_u = m*vec_v + sig*vec_w;
			    vec_u *= Fb.WeightsND(k);
			  }
			else
			  {
			    Mlt(var.Glob_matMass_elem(i).GetBh(k), vec_u, vec_v);
			    Mlt(var.Glob_matMass_elem(i).GetBhSigma(k), vec_u, vec_w);
			    
			    vec_u = m*vec_v + sig*vec_w;
			  }

			ProdEquad(3*k) += vec_u(0);
			ProdEquad(3*k+1) += vec_u(1);
			ProdEquad(3*k+2) += vec_u(2);
		      }
		    
		    // mass term for H
		    if (affine)
		      jacobian = var.Glob_jacobian(i)(0)*Fb.WeightsND(k);
		    else
		      jacobian = var.Glob_jacobian(i)(k);
		    
		    ExtractVector(B, offset_H + 3*k, vec_u);

		    var.ref_mu(ref).MltMatrix(var, i, k, vec_u, vec_v);
		    poids = coef_H*jacobian;
                    T poids2 = m*m_iomega*jacobian;

		    vec_H(0) = vec_v(0)*poids;
		    vec_H(1) = vec_v(1)*poids;
		    vec_H(2) = vec_v(2)*poids;

		    if (pml_elt)
		      {
			tau = var.GetTauPML(i1, k);
			T poids_sig = sig*jacobian;
			ExtractVector(B, offset_Hs + 3*k, vec_w);
                        
			vec_H(0) += poids_sig*tau(1)*vec_v(0);
			vec_H(1) += poids_sig*tau(2)*vec_v(1);
			vec_H(2) += poids_sig*tau(0)*vec_v(2);
			
			if (trans.NoTrans())
			  {			    
			    vec_Hs(0) = poids2*(vec_w(0) - vec_u(0)) + poids_sig*(tau(2)*vec_w(0) - tau(0)*vec_u(0));
			    vec_Hs(1) = poids2*(vec_w(1) - vec_u(1)) + poids_sig*(tau(0)*vec_w(1) - tau(1)*vec_u(1));
			    vec_Hs(2) = poids2*(vec_w(2) - vec_u(2)) + poids_sig*(tau(1)*vec_w(2) - tau(2)*vec_u(2));
			  }
			else
			  {
			    vec_H(0) -= (poids2 + poids_sig*tau(0))*vec_w(0);
			    vec_H(1) -= (poids2 + poids_sig*tau(1))*vec_w(1);
			    vec_H(2) -= (poids2 + poids_sig*tau(2))*vec_w(2);
			    
			    vec_Hs(0) = (poids2 + poids_sig*tau(2))*vec_w(0);
			    vec_Hs(1) = (poids2 + poids_sig*tau(0))*vec_w(1);
			    vec_Hs(2) = (poids2 + poids_sig*tau(1))*vec_w(2); 
			  }
			
			AddVector(vec_Hs, offset_Hs + 3*k, C);
		      }

		    AddVector(vec_H, offset_H+3*k, C);
		  }
	      }
	    
            if (mass && !diag_mass && pml_elt)
              {
                T epsilon = var.ref_epsilon(ref).GetConstant()(0, 0)*sig;
		for (int k = 0; k < nb_points_quad; k++)
		  {
                    // we compute vec_v = DF_i^{*-1} Echap
                    CopyVector(Equad, k, vec_u);
                    CopyVector(Equad_star, k, vec_P);
                    if (affine)
                      {
                        Real_wp jacobian = var.Glob_jacobian(i)(0);
                        MltTrans(var.Glob_DFjm1(i)(0), vec_u, vec_v);
                        vec_v *= Real_wp(1)/jacobian;

                        MltTrans(var.Glob_DFjm1(i)(0), vec_P, vec_Q);
                        vec_Q *= Real_wp(1)/jacobian;
                      }
                    else
                      {
                        Real_wp jacobian = var.Glob_jacobian(i)(k) / Fb.WeightsND(k);
                        MltTrans(var.Glob_DFjm1(i)(k), vec_u, vec_v);
                        vec_v *= Real_wp(1)/jacobian;

                        MltTrans(var.Glob_DFjm1(i)(k), vec_P, vec_Q);
                        vec_Q *= Real_wp(1)/jacobian;
                      }

                    tau = var.GetTauPML(i1, k);
                    T poids = Fb.WeightsND(k)*epsilon;
                    vec_w(0) = tau(1)*vec_v(0)*poids;
                    vec_w(1) = tau(2)*vec_v(1)*poids;
                    vec_w(2) = tau(0)*vec_v(2)*poids;
                    
                    if (trans.NoTrans())
                      {
                        poids = m*m_iomega*Fb.WeightsND(k);
                        vec_Es(0) = poids*(vec_Q(0) - vec_v(0));
                        vec_Es(1) = poids*(vec_Q(1) - vec_v(1));
                        vec_Es(2) = poids*(vec_Q(2) - vec_v(2));
                        
                        poids = sig*Fb.WeightsND(k);
                        vec_Es(0) += poids*(tau(2)*vec_Q(0) - tau(0)*vec_v(0));
                        vec_Es(1) += poids*(tau(0)*vec_Q(1) - tau(1)*vec_v(1));
                        vec_Es(2) += poids*(tau(1)*vec_Q(2) - tau(2)*vec_v(2));
                      }
                    else
                      {
                        poids = m*m_iomega*Fb.WeightsND(k);
                        vec_Es(0) = poids*vec_Q(0); vec_w(0) -= poids*vec_Q(0);
                        vec_Es(1) = poids*vec_Q(1); vec_w(1) -= poids*vec_Q(1);
                        vec_Es(2) = poids*vec_Q(2); vec_w(2) -= poids*vec_Q(2);
                        
                        poids = sig*Fb.WeightsND(k);
                        vec_Es(0) += poids*tau(2)*vec_Q(0); vec_w(0) -= poids*tau(0)*vec_Q(0);
                        vec_Es(1) += poids*tau(0)*vec_Q(1); vec_w(1) -= poids*tau(1)*vec_Q(1);
                        vec_Es(2) += poids*tau(1)*vec_Q(2); vec_w(2) -= poids*tau(2)*vec_Q(2);
                      }
                    
                    if (affine)
                      {
                        Mlt(var.Glob_DFjm1(i)(0), vec_w, vec_u);
                        Mlt(var.Glob_DFjm1(i)(0), vec_Es, vec_Q);
                      }
                    else
                      {
                        Mlt(var.Glob_DFjm1(i)(k), vec_w, vec_u);
                        Mlt(var.Glob_DFjm1(i)(k), vec_Es, vec_Q);
                      }
                    
                    ProdEquad(3*k) += vec_u(0);
                    ProdEquad(3*k+1) += vec_u(1);
                    ProdEquad(3*k+2) += vec_u(2);

                    ProdEquad_star(3*k) = vec_Q(0);
                    ProdEquad_star(3*k+1) = vec_Q(1);
                    ProdEquad_star(3*k+2) = vec_Q(2);
                  }
                
                Fb.ApplyCh(ProdEquad_star, contrib);
                ProdEloc_star(0) += contrib;
              }
            
	    if (Fb.UseQuadraturePointsForRh() || drude_stiff)
	      {
		Fb.ApplyCh(ProdEquad, contrib);
		ProdEloc(0) += contrib;
	      }
	    else
	      {
		if (mass && !diag_mass)
		  {
		    Fb.ApplyCh(ProdEquad, Eloc(0));
		    ProdEloc(0) += Eloc(0);
		  }
	      }
	    
	    if (mass && diag_mass)
	      {
		for (int j = 0; j < Fb.GetNbDof(); j++)
		  {
		    T vloc = m*mass_elt.GetDh(j) + sig*mass_elt.GetDhSigma(j);
		    ProdEloc(0)(j) += vloc*Eloc(0)(j);
		  }

		if (pml_elt)
		  for (int j = 0; j < Fb.GetNbDof(); j++)
		    {
                      int num_dof = mesh_num.Element(i).GetNumberDof(j);
                      int npml = mesh_num.GetDofPML(num_dof);
                      if (npml >= 0)
                        {
                          tau_r = mass_elt.GetTauPML(j);
                          int n = Fb.GetCoordinateDof(j);
                          T vloc = mass_elt.GetDh(j) * coef_E * tau_r(permut120(n));
                          ProdEloc(0)(j) += vloc*Eloc(0)(j);
                          
                          if (trans.Trans())
                            {
                              vloc = mass_elt.GetDh(j)*coefM_Es;
                              ProdEloc_star(0)(j) += vloc*Eloc_star(0)(j);
                              ProdEloc(0)(j) -= vloc*Eloc_star(0)(j);
                              
                              vloc = mass_elt.GetDh(j)*coef_Es;
                              ProdEloc_star(0)(j) += vloc*Eloc_star(0)(j)*tau_r(permut201(n));
                              ProdEloc(0)(j) -= vloc*Eloc_star(0)(j)*tau_r(n);
                            }
                          else
                            {
                              vloc = mass_elt.GetDh(j)*coefM_Es;
                              ProdEloc_star(0)(j) += vloc*(Eloc_star(0)(j) - Eloc(0)(j));
                              
                              vloc = mass_elt.GetDh(j)*coef_Es;
                              ProdEloc_star(0)(j) += vloc*(Eloc_star(0)(j)*tau_r(permut201(n)) - Eloc(0)(j)*tau_r(n));
                            }
                        }
		    }
	      }
	    
	    // contributions are added to C
	    var.AddLocalUnknownVector(one, ProdEloc, i, C);

	    if ((pml_elt) || (close_pml_elt && trans.Trans()))
	      var.AddLocalUnknownVector(one, ProdEloc_star, i, Cstar);
	  }

	if (var.GetNbGlobalEltPML() > 0)
	  for (int i = 0; i < Nvol; i++)
	    {
	      int npml = mesh_num.GetDofPML(i);
	      if (npml >= 0)
		C(Nvol + npml) += Cstar(i);
              else
                C(i) += Cstar(i);
	    }
      }
    else
      {
	TinyVector<Vector<T>, 1> Eh_dof;
	Vector<T> Eh_quad, H_quad;
	TinyVector<Vector<T>, 1> Prod_Eh;

	bool stiff = true;
	if (nat_mat.GetCoefStiffness() == zero)
	  stiff = false;
	
	bool mass = true;
	if ((nat_mat.GetCoefDamping() == zero) && (nat_mat.GetCoefMass() == zero))
	  mass = false;

	Vector<T> Eh_mass, ProdEh_quad;
	
	const HarmonicMaxwell_3D<Complexe>& var_maxwell = static_cast<const HarmonicMaxwell_3D<Complexe>& >(var);
	
	// loop over all elements of the mesh
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    const ElementReference<Dimension3, 2>& Fb = var.GetReferenceElementHcurl(i);
	    const HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass_elt = var.Glob_matMass_elem(i);
 
	    // local E_h
	    int nb_dof_loc = Fb.GetNbDof();
	    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	    bool diag_mass = mass_elt.IsDiagonalMass();

	    Prod_Eh(0).Reallocate(nb_dof_loc);
	    
	    TinyMatrix<T0, Symmetric, 3, 3> Bmass, Astiff;
	    TinyVector<T, 3> vh1, vh2, Eh_loc, H_loc;

    	    if (stiff)
	      H_quad.Reallocate(3*nb_points_quad);
	    
	    if (stiff)
	      {
		if (Fb.UseQuadraturePointsForRh())
		  Eh_mass.Reallocate(3*nb_points_quad);
		else
		  Eh_mass.Reallocate(nb_dof_loc);
	      }
	    
	    // we retrieve local vector E
	    var.GetLocalUnknownVector(B, i, Eh_dof);
	    
	    // interpolation of E on quadrature points
	    if ((stiff && Fb.UseQuadraturePointsForRh()) || ((!diag_mass) && mass) || var.linearize_drude)
	      {
		Eh_quad.Reallocate(3*nb_points_quad);
		Fb.ApplyChTranspose(Eh_dof(0), Eh_quad);
	      }
	    
	    // curl(E) on quadrature points
	    if (stiff)
	      {
		if (Fb.UseQuadraturePointsForRh())
		  Fb.ApplyRhQuadratureTranspose(Eh_quad, H_quad);
		else
		  Fb.ApplyRhTranspose(Eh_dof(0), H_quad);
	      }
	    
            if (var.linearize_drude)
              {
                int ref_domain = mesh.Element(i).GetReference();
                Vector<T0> poidsP, poidsE, poidsDiag; 
                bool affine = mesh.IsElementAffine(i);            
                
                int offset_P = var.GetNbDof() - var.nb_dof_drude_vec;
                offset_P += var.OffsetDofDrudeV(i); 
                int nPole = var.ref_drude(ref_domain).gamma.GetM();
                poidsP.Reallocate(nPole);
                poidsE.Reallocate(nPole);
                poidsDiag.Reallocate(nPole);
                ProdEh_quad.Reallocate(3*nb_points_quad);
                for (int k = 0; k < nPole; k++)
                  {
                    poidsP(k) = -var.ref_drude(ref_domain).eps_omega_p2(k)*nat_mat.GetCoefStiffness()
                      - var.ref_drude(ref_domain).eps_sigma(k)*nat_mat.GetCoefDamping()*m_iomega;
                    if (trans.NoTrans())
                      poidsE(k) = m_omega2*nat_mat.GetCoefMass();
                    else
                      {
                        poidsE(k) = poidsP(k);
                        poidsP(k) = m_omega2*nat_mat.GetCoefMass();
                      }
                
                    poidsDiag(k) = m_omega2*nat_mat.GetCoefMass()
                      + m_iomega*nat_mat.GetCoefDamping()*var.ref_drude(ref_domain).gamma(k)
                      + var.ref_drude(ref_domain).omega_02(k)*nat_mat.GetCoefStiffness();
                  }

                TinyVector<T, 3> vec_u, vec_v, vec_P, vec_Q; Real_wp jacobian;
                for (int j = 0; j < nb_points_quad; j++)
                  {
                    vec_P.Zero();
                    for (int k = 0; k < nPole; k++)
                      {
                        int offset_P2 = offset_P + 3*(nb_points_quad*k + j);
                        vec_P(0) += B(offset_P2)*poidsE(k);
                        vec_P(1) += B(offset_P2+1)*poidsE(k);
                        vec_P(2) += B(offset_P2+2)*poidsE(k);
                      }

                    CopyVector(Eh_quad, j, vec_u);
                    if (affine)
                      {
                        Mlt(var.Glob_DFjm1(i)(0), vec_P, vec_Q);
                        MltTrans(var.Glob_DFjm1(i)(0), vec_u, vec_v);
                        jacobian = var.Glob_jacobian(i)(0)*Fb.WeightsND(j);
                      }
                    else
                      {
                        Mlt(var.Glob_DFjm1(i)(j), vec_P, vec_Q);
                        MltTrans(var.Glob_DFjm1(i)(j), vec_u, vec_v);
                        jacobian = var.Glob_jacobian(i)(j);
                      }
                    
                    ProdEh_quad(3*j) = Fb.WeightsND(j)*vec_Q(0);
                    ProdEh_quad(3*j+1) = Fb.WeightsND(j)*vec_Q(1);
                    ProdEh_quad(3*j+2) = Fb.WeightsND(j)*vec_Q(2);
                    
                    T coef1 = jacobian*alpha;
                    T coef2 = alpha*Fb.WeightsND(j);
                    for (int k = 0; k < nPole; k++)
                      {
                        int offset_P2 = offset_P + 3*(nb_points_quad*k + j);
                        C(offset_P2) = coef1*poidsDiag(k)*B(offset_P2) + coef2*poidsP(k)*vec_v(0);
                        C(offset_P2+1) = coef1*poidsDiag(k)*B(offset_P2+1) + coef2*poidsP(k)*vec_v(1);
                        C(offset_P2+2) = coef1*poidsDiag(k)*B(offset_P2+2) + coef2*poidsP(k)*vec_v(2);
                      }
                  }
              }
            
	    if (var.UseNumericalIntegration(i))
	      for (int j = 0; j < nb_points_quad; j++)
		{
		  if ((!diag_mass) && (mass))
		    {
		      var_maxwell.GetMassMatrix(i, j, nat_mat, Bmass);
		      
		      Eh_loc(0) = Eh_quad(3*j);
		      Eh_loc(1) = Eh_quad(3*j+1);
		      Eh_loc(2) = Eh_quad(3*j+2);
		      
		      Mlt(Bmass, Eh_loc, vh1);
		      
		      Eh_quad(3*j) = vh1(0);
		      Eh_quad(3*j+1) = vh1(1);
		      Eh_quad(3*j+2) = vh1(2);
		    }
		  
		  if (stiff)
		    {	
		      var_maxwell.GetStiffMatrix(i, j, nat_mat, Astiff);
		      
		      H_loc(0) = H_quad(3*j);
		      H_loc(1) = H_quad(3*j+1);
		      H_loc(2) = H_quad(3*j+2);
		      
		      Mlt(Astiff, H_loc, vh2);
		      
		      H_quad(3*j) = vh2(0);
		      H_quad(3*j+1) = vh2(1);
		      H_quad(3*j+2) = vh2(2);
		    }
		}
	    else
	      {
		if ((!diag_mass) && (mass))
		  var_maxwell.GetMassMatrix(i, 0, nat_mat, Bmass);    

		if (stiff)
		  var_maxwell.GetStiffMatrix(i, 0, nat_mat, Astiff);
		
		for (int j = 0; j < nb_points_quad; j++)
		  {
		    if ((!diag_mass) && (mass))
		      {
			Eh_loc(0) = Eh_quad(3*j);
			Eh_loc(1) = Eh_quad(3*j+1);
			Eh_loc(2) = Eh_quad(3*j+2);
			
			Mlt(Bmass, Eh_loc, vh1);
			Mlt(Fb.WeightsND(j), vh1);
			
			Eh_quad(3*j) = vh1(0);
			Eh_quad(3*j+1) = vh1(1);
			Eh_quad(3*j+2) = vh1(2);
		      }
		    
		    if (stiff)
		      {	
			H_loc(0) = H_quad(3*j);
			H_loc(1) = H_quad(3*j+1);
			H_loc(2) = H_quad(3*j+2);
			
			Mlt(Astiff, H_loc, vh2);
			Mlt(Fb.WeightsND(j), vh2);
			
			H_quad(3*j) = vh2(0);
			H_quad(3*j+1) = vh2(1);
			H_quad(3*j+2) = vh2(2);
		      }
		  }
	      }
    
	    // integration against curl(phi)
	    if (stiff)
	      {
		if (Fb.UseQuadraturePointsForRh())
		  {
		    Fb.ApplyRhQuadrature(H_quad, Eh_mass);
		    if ((!diag_mass) && (mass))
		      Eh_quad += Eh_mass;
		    else
		      Eh_quad = Eh_mass;
		  }
		else
		  Fb.ApplyRh(H_quad, Eh_mass);
	      }
            if (var.linearize_drude)
              {
                if ((stiff && Fb.UseQuadraturePointsForRh()) || ((!diag_mass) && mass))
                  Eh_quad += ProdEh_quad;
                else
                  Eh_quad = ProdEh_quad;
              }
              
	    if ((stiff && Fb.UseQuadraturePointsForRh()) || ((!diag_mass) && mass) || var.linearize_drude)
	      Fb.ApplyCh(Eh_quad, Prod_Eh(0));
	    else
	      Prod_Eh(0).Zero();
	    
	    if (stiff && (!Fb.UseQuadraturePointsForRh()))
	      Prod_Eh(0) += Eh_mass;

	    if (diag_mass && mass)
	      {
		for (int j = 0; j < Fb.GetNbDof(); j++)
		  {
		    T vloc = nat_mat.GetCoefMass()*mass_elt.GetDh(j) + nat_mat.GetCoefDamping()*mass_elt.GetDhSigma(j);
		    Prod_Eh(0)(j) += vloc*Eh_dof(0)(j);
		  }
	      }
	    
	    // we add contribution to C
	    var.AddLocalUnknownVector(alpha, Prod_Eh, i, C);
	  }
      }
    
    A.AddExtraBoundaryTerms(alpha, B, C);
    
    A.ApplyLeftScaling(B2, C2, B, C);    
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C2);
  }


  //! generic matrix-vector product for edge finite elements in 3-D
  //! and Maxwell equations (HDG formulation)
  /*
    this algorithm separates the geometry \f$ DF_i^{*-1} \f$ from the 
    spatial derivatives. This allows low-storage and fast-algorithm
    in \f$ O(r^4) \f$ instead of \f$ O(r^6) \f$ if the full matrix is stored, where
    r is the order of approximation
  */
  template <class T, class TypeEquation>
  void MltAddHcurlHdg3D(const T& alpha, const GlobalGenericMatrix<T>& nat_mat0,
			const SeldonTranspose& trans, int level,
			const FemMatrixFreeClass<T, TypeEquation>& A,
			const Vector<T>& B2, const T& beta, Vector<T>& C2, bool assemble)
  {
    
#ifdef SELDON_WITH_MPI
    Vector<T> B(B2), C(C2);
#else
    Vector<T> B, C;
#endif
    A.ApplyRightScaling(B2, C2, B, C);    

    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (beta == zero)
      C.Zero();
    else if (beta != one)
      Mlt(beta, C);

    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    typedef typename TypeEquation::Complexe Complexe;
    
    Complexe m_iomega; var.GetMiomega(m_iomega);
    int Nvol = var.GetNbMainUnknownDof();
    
    T m = nat_mat0.GetCoefMass()*alpha;
    T s = nat_mat0.GetCoefStiffness()*alpha;
    T sig = nat_mat0.GetCoefDamping()*alpha;
    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(m);
    nat_mat.SetCoefDamping(sig);
    nat_mat.SetCoefStiffness(s);
    
    bool sym = A.GetSymmetrizationUse();
    bool stiff = false, mass = true;
    if (s != zero)
      stiff = true;
    
    if (m != zero)
      mass = true;
    
    if (sig != zero)
      mass = true;
    
    T coef_H = m_iomega*m;
    bool signE = false, signH = false;
    if (sym)
      coef_H = -coef_H;
	
    if (!sym)
      {
	if (trans.Trans())
	  signE = true;
	else
	  signH = true;
      }
    
    TinyVector<Vector<T>, 1> Lloc, ProdLloc;
    Vector<T> Eloc, Hloc, Equad, Hquad, curlH_quad, ProdEquad, ProdEloc;
    TinyVector<T, 3> vec_u, vec_v, vec_w;
    Vector<T> Lambda, ProdLambda;
		
    int nb_dof_Lglob = var.GetOffsetDofUnknown(1);
    for (int i0 = 0; i0 < mesh.GetNbElt(); i0++)
      {
	//int i = list_level.GetElementNumber(i0);
	int i = i0;

	int offset_E = nb_dof_Lglob + var.GetOffsetDofV(i);
	int offset_H = nb_dof_Lglob + Nvol + var.GetOffsetDofV(i);
		
	const ElementReference<Dimension3, 2>& Fb = var.GetReferenceElementHcurl(i);
	
	int nb_dof_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	
	//bool affine = mesh.IsElementAffine(i);
	//int ref = mesh.Element(i).GetReference();

	int nb_dof_L = mesh_num.Element(i).GetNbDof();
	Lloc(0).Reallocate(nb_dof_L); Lloc(0).Zero();
	ProdLloc(0).Reallocate(nb_dof_L); ProdLloc(0).Zero();
	Eloc.Reallocate(nb_dof_elt);
	Hloc.Reallocate(nb_dof_elt);
	ProdEloc.Reallocate(nb_dof_elt); ProdEloc.Zero();
	curlH_quad.Reallocate(3*nb_points_quad);
	
	// we retrieve values of lambda on the element
	mesh_num.number_map.GetLocalUnknownVector(mesh_num, B, i, Lloc);

	// we retrieve values of E and H on the element
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    Eloc(j) = B(offset_E + j);
	    Hloc(j) = B(offset_H + j);
	  }
	
	if (Fb.LumpedMassMatrix())
	  {
	    const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	    const Matrix<Real_wp, General, ArrayRowSparse>& R = Fb_hex.GetPermutedRh();
	    
	    // stiffness terms
	    if (stiff)
	      {
		// curl of H
		//Fb.ApplyRhTranspose(Hloc, curlH_quad);
		Mlt(SeldonTrans, R, Hloc, curlH_quad);
		Mlt(R, Eloc, ProdEloc);

		T coefE = s, coefH = s;
		if (signE)
		  coefE = -coefE;
		
		if (signH)
		  coefH = -coefH;
		
		for (int j = 0; j < nb_dof_elt; j++)
		  {
		    C(offset_E+j) -= coefE*curlH_quad(j);
		    C(offset_H+j) -= coefH*ProdEloc(j);
		  }
	      }

	    if (mass)
	      {
		bool variable = var.UseNumericalIntegration(i);
		for (int j = 0; j < nb_points_quad; j++)
		  {		     
		    int jx, jy, jz;
		    Fb_hex.GetDofNumber_FromPointNode(j, jx, jy, jz);

		    vec_u.Init(Eloc(jx), Eloc(jy), Eloc(jz));

		    if (!variable)
		      {
			Mlt(var.Glob_matMass_elem(i).GetBh(0), vec_u, vec_v);
			Mlt(var.Glob_matMass_elem(i).GetBhSigma(0), vec_u, vec_w);
			    
			vec_u = m*vec_v + sig*vec_w;
			vec_u *= Fb.WeightsND(j);
		      }
		    else
		      {
			Mlt(var.Glob_matMass_elem(i).GetBh(j), vec_u, vec_v);
			Mlt(var.Glob_matMass_elem(i).GetBhSigma(j), vec_u, vec_w);
			
			vec_u = m*vec_v + sig*vec_w;
		      }
		    
		    C(offset_E+jx) += vec_u(0);
		    C(offset_E+jy) += vec_u(1);
		    C(offset_E+jz) += vec_u(2);

		    vec_u.Init(Hloc(jx), Hloc(jy), Hloc(jz));

		    if (!variable)
		      {
			Mlt(var.Glob_matMass_elem(i).GetAh(0), vec_u, vec_v);
			    
			vec_u = m*Fb.WeightsND(j)*vec_v;
		      }
		    else
		      {
			Mlt(var.Glob_matMass_elem(i).GetAh(j), vec_u, vec_v);
			vec_u = m*vec_v;
		      }

		    if (sym)
		      vec_u = -vec_u;
		    
		    C(offset_H+jx) += vec_u(0);
		    C(offset_H+jy) += vec_u(1);
		    C(offset_H+jz) += vec_u(2);
		  }
	      }

	    if (stiff)
	      {
		// loop over boundaries of the element
		int offsetS = 0, offset_L = 0; TinyVector<T, 2> vec_us, vec_vs;
		for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
		  {
		    int num_face, ref_boundary, rf, rot; bool new_face; Complexe phase, phase_conj;
		    var.GetGeometryPhaseData(i, num_loc,
					     num_face, ref_boundary, rf, new_face, rot,
					     phase, phase_conj);
		    
		    const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(var.GetSurfaceFiniteElement(num_face));
		    
		    int order = Fb_quad.GetOrder();
		    Lambda.Reallocate(Fb_quad.GetNbDof());
		    ProdLambda.Reallocate(Fb_quad.GetNbDof());
		    for (int j = 0; j < Fb_quad.GetNbDof(); j++)
		      {
			if ((j >= 2*(order+1)) && (j < 4*(order+1)))
			  Lambda(j) = -Lloc(0)(offset_L+j);
			else
			  Lambda(j) = Lloc(0)(offset_L+j);
		      }

		    bool abc_cond = false, dir_cond = false;
		    if (ref_boundary != 0)
		      {
			int cond_ref = mesh.GetBoundaryCondition(ref_boundary);
			if (cond_ref == BoundaryConditionEnum::LINE_DIRICHLET)
			  dir_cond = true;
			else if (cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
			  {
			    // treated correctly
			  }
			else if (cond_ref == BoundaryConditionEnum::LINE_ABSORBING)
			  abc_cond = true;
			else
			  {
			    cout << "Boundary condition not treated" << endl;
			    abort();
			  }
		      }
		    
		    bool signF = true;
		    if ((num_loc == 0) || (num_loc == 2) || (num_loc == 4))
		      signF = false;

		    if (!sym)
		      signF = !signF;
		    
		    for (int k = 0; k < Fb_quad.GetNbPointsQuadratureInside(); k++)
		      {
			int npoint = Fb.GetQuadNumber(num_loc, k);
			int kx, ky, kz;
			Fb_hex.GetDofNumber_FromPointNode(npoint, kx, ky, kz);
			int jx = kx, jy = ky;
			if ((num_loc == 0) || (num_loc == 5))
			  { jx = ky; jy = kz; }
			else if ((num_loc == 1) || (num_loc == 4))
			  { jx = kx; jy = kz; }

			kx = Fb_quad.GetXdofNumber(k);
			ky = Fb_quad.GetYdofNumber(k);

			if (!signE)
			  vec_us.Init(Eloc(jx)-Lambda(kx), Eloc(jy)-Lambda(ky));
			else
			  vec_us.Init(Eloc(jx)+Lambda(kx), Eloc(jy)+Lambda(ky));
			
			Mlt(var.Glob_matMass_elem(i).GetBhSurf(offsetS+k), vec_us, vec_vs);
			
			C(offset_E+jx) += s*vec_vs(0);
			C(offset_E+jy) += s*vec_vs(1);

			if (abc_cond)
			  {
			    if (!signE)
			      {
				vec_us(0) -= Lambda(kx); vec_us(1) -= Lambda(ky);
			      }
			    else
			      {
				vec_us(0) += Lambda(kx); vec_us(1) += Lambda(ky);
			      }
			    
			    Mlt(var.Glob_matMass_elem(i).GetBhSurf(offsetS+k), vec_us, vec_vs);
			  }
			
			if (signH)
			  vec_vs = -vec_vs;
			
			ProdLambda(kx) = -s*vec_vs(0);
			ProdLambda(ky) = -s*vec_vs(1);

			if (!dir_cond)
			  {
			    T coef = s*Fb_quad.WeightsND(k);
			    if (signF)
			      coef = -coef;

			    C(offset_H+jx) -= coef*Lambda(ky);
			    C(offset_H+jy) += coef*Lambda(kx);

			    ProdLambda(ky) -= coef*Hloc(jx);
			    ProdLambda(kx) += coef*Hloc(jy);
			  }
		      }
		    
		    for (int j = 0; j < Fb_quad.GetNbDof(); j++)
		      {
			if ((j >= 2*(order+1)) && (j < 4*(order+1)))
			  ProdLloc(0)(offset_L+j) = -ProdLambda(j);
			else
			  ProdLloc(0)(offset_L+j) = ProdLambda(j);
		      }
		    
		    offset_L += Fb_quad.GetNbDof();
		    offsetS += Fb_quad.GetNbPointsQuadratureInside();
		  }
	      }
	  }
	else
	  {
	    // values of E and H on quadrature points
	    Fb.ApplyChTranspose(Eloc, Equad);
	    Fb.ApplyChTranspose(Hloc, Hquad);
	  }
	    
	// contributions are added to C
	var.AddLocalUnknownVector(one, ProdLloc, i, C);	
      }
    
    A.AddExtraBoundaryTerms(alpha, B, C);
    
    A.ApplyLeftScaling(B2, C2, B, C);    
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C2);
  }

  
  //! computation of unknowns on the quadrature points of element
  template<class Complexe, class TypeEquation>
  void GetExtrapolationSecondOrder_Sipg(const EllipticProblem<TypeEquation>& var,
                                        int num_elem, const Vector<Complexe>& B,
                                        ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
                                        const GlobalGenericMatrix<Complexe>& nat_mat,
                                        const ElementReference<Dimension3, 2>& Fb)
  {
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    int offset_face = mesh_num.OffsetQuadElementNumber(num_elem);
    bool affine = mesh.IsElementAffine(num_elem);    

    // we get values of E on the element
    int nb_pts_quad = Fb.GetNbPointsQuadratureInside();  
    Vector<Complexe> Eloc(Fb.GetNbDof()), valE, curlE;
    Vector<Complexe>& Equad = var_extra.GetEnQuad(num_elem, Fb);
    Vector<Complexe>& rotEquad = var_extra.GetCurlEnQuad(num_elem, Fb);

    int offset_elt = mesh_num.Element(num_elem).GetNumberDof(0);
    for (int j = 0; j < Fb.GetNbDof(); j++)
      Eloc(j) = B(offset_elt + j);
    
    if (Fb.UseQuadraturePointsForSh())
      {
        // computation of E on quadrature points if needed
        Equad.Reallocate(3*nb_pts_quad);
        Fb.ApplyChTranspose(Eloc, Equad);
        
        // and curl of E
        rotEquad.Reallocate(3*nb_pts_quad);
        if (Fb.UseQuadraturePointsForRh())
          Fb.ApplyRhQuadratureTranspose(Equad, rotEquad);
        else
          Fb.ApplyRhTranspose(Eloc, rotEquad);
      }
    
    int ref_d = mesh.Element(num_elem).GetReference();
    int num_loc2, num_elem2, rot = 0;
    int ref, cond, offset_face2; Complexe phase, phase_conj; bool new_face;
    TinyVector<Complexe, 3> vec_u, vec_v, vec_E, vec_H;
    int num_point = nb_pts_quad; Real_wp dsj;
    
    // loop over boundaries of the element
    for (int num_loc = 0; num_loc < mesh.Element(num_elem).GetNbBoundary(); num_loc++)
      {        
	int num_face = mesh.Element(num_elem).numBoundary(num_loc);
        int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        if (var.FaceHasToBeConsideredForBoundaryIntegral(num_face))
          {            
            // the face is involved in the fluxes terms
            int rf = mesh_num.GetOrderQuadrature(num_face);        
            var_extra.GetFaceOrientationAndNumber(var, mesh, num_face,
                                                  num_elem, num_loc, offset_face,
                                                  ref, cond, new_face, num_elem2, num_loc2,
                                                  offset_face2, rot, phase, phase_conj);
            
            const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
              GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
            
            // computation of E on quadrature points of face
            valE.Reallocate(3*nb_points_face);            
            if (Fb.UseQuadraturePointsForSh())
              Fb.ApplyShQuadratureTranspose(num_loc, Equad, valE, rf);
            else
              Fb.ApplyShTranspose(num_loc, Eloc, valE, rf);
            
            // computation of curl E on quadrature points of face
            curlE.Reallocate(3*nb_points_face);
            if (Fb.UseQuadraturePointsForSh())
              {
                if (Fb.UseQuadraturePointsForRh())
                  Fb.ApplyShQuadratureTranspose(num_loc, rotEquad, curlE, rf);
                else
                  Fb.ApplyNablaShQuadratureTranspose(num_loc, Equad, curlE, rf);
              }
            else
              Fb.ApplyNablaShTranspose(num_loc, Eloc, curlE, rf);
            
            int offset_neighbor = 0;
            bool neighbor_face = false;
	    if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
	      {
                neighbor_face = true;
                offset_neighbor = var.GetOffsetNeighboringFace(num_face);
              }
            
            // loop over quadrature points
            for (int j = 0; j < nb_points_face; j++)
              {
                int jrot = FacesQuadRotation(rot, j);
                
                if (new_face)
                  dsj = var.Glob_dsj(num_face)(j);
                else
                  dsj = var.Glob_dsj(num_face)(jrot);
                
                int j2 = num_point + j - nb_pts_quad;
                if (affine)
                  j2 = 0;
                
                // storing u
                CopyVector(valE, j, vec_u);
                
                TimesProd(Fb.NormaleLoc(num_loc), vec_u, vec_v);
                Mlt(var.Glob_DFj(num_elem)(j2), vec_v, vec_E);
                vec_E *= 1.0/dsj;
                
                var_extra.extrapolU(0)(offset_face+j) = vec_E(0);
                var_extra.extrapolU(1)(offset_face+j) = vec_E(1);
                var_extra.extrapolU(2)(offset_face+j) = vec_E(2);
                
                // storing 1/mu rot u
                CopyVector(curlE, j, vec_u);
                
                Mlt(var.Glob_DFj(num_elem)(j2), vec_u, vec_v);
                vec_v *= var.Glob_invJacobian(num_elem)(j2);
                
                var.ref_invMu(ref_d).MltMatrix(var, num_elem, num_point+j, vec_v, vec_H);
                
                var_extra.extrapolRotU(0)(offset_face+j) = vec_H(0); 
                var_extra.extrapolRotU(1)(offset_face+j) = vec_H(1); 
                var_extra.extrapolRotU(2)(offset_face+j) = vec_H(2); 
                
                if (neighbor_face)
                  {
                    var_extra.Uneighbor(0)(offset_neighbor+j) = vec_E(0);
                    var_extra.Uneighbor(1)(offset_neighbor+j) = vec_E(1);
                    var_extra.Uneighbor(2)(offset_neighbor+j) = vec_E(2);
                    
                    var_extra.rotUneighbor(0)(offset_neighbor+j) = vec_H(0);
                    var_extra.rotUneighbor(1)(offset_neighbor+j) = vec_H(1);
                    var_extra.rotUneighbor(2)(offset_neighbor+j) = vec_H(2);
                  }
              }                      
          } 
        
        offset_face += nb_points_face;
        num_point += nb_points_face;
      }
  }
  
  
  //! computation of 1/mu rot u for SIPG formulation and Maxwell's equations
  //! values of 1/mu rot u on quadrature points are exchanged between processors
  template<class Complexe, class TypeEquation>
  void GetExtrapolationAndExchange_Sipg(const GlobalGenericMatrix<Complexe>& nat_mat,
                                        const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                        const Vector<Complexe>& B,
                                        ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra)
  {
    const Mesh<Dimension3>& mesh = A.var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = A.var.GetMeshNumbering(0);
    
    int size_extrapol = 0;
    int size_neighbor = 0;

    // we will store U on quadrature points of each boundary
    // and exchange U on interface (Uneighbor)
    size_extrapol = mesh_num.OffsetQuadElementNumber(mesh.GetNbElt());
    size_neighbor = A.var.GetNbPointsQuadratureNeighbor();
    
    for (int k = 0; k < 3; k++)
      {        
        var_extra.extrapolU(k).Reallocate(size_extrapol);
        var_extra.extrapolRotU(k).Reallocate(size_extrapol);        

        var_extra.Uneighbor(k).Reallocate(size_neighbor);
        var_extra.rotUneighbor(k).Reallocate(size_neighbor);
      }
        
    // loop over all elements of the mesh to compute u on quadrature points of faces
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        GetExtrapolationSecondOrder_Sipg(A.var, i, B, var_extra, nat_mat,
                                         A.var.GetReferenceElementHcurl(i));
      }
    
    
#ifdef SELDON_WITH_MPI    
    // exchanging datas between processors
    if (A.var.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
        var_extra.ReallocateExchangeVector(2);
        
        if (size_neighbor > 0)
          {
            // exchanging values of u and rot u
            for (int k = 0; k < 3; k++)
              {
                A.var.ExchangeUfaceDomains(var_extra.Uneighbor(k), var_extra.Xsend(0), var_extra.Xsend_tmp(0),
                                           var_extra.Xneighbor(0), var_extra.Xneighbor_tmp(0), var_extra.request(0), 90);
            
                // sending and receiving mu du/dn
                A.var.ExchangeUfaceDomains(var_extra.rotUneighbor(k), var_extra.Xsend(1), var_extra.Xsend_tmp(1),
                                           var_extra.Xneighbor(1), var_extra.Xneighbor_tmp(1), var_extra.request(1), 91);
            
                // finalizing transfers if needed
                A.var.GetUfaceDomains(var_extra.Uneighbor(k), var_extra.Xsend(0), var_extra.Xsend_tmp(0),
                                      var_extra.Xneighbor(0), var_extra.Xneighbor_tmp(0), var_extra.request(0), 90);
                
                A.var.GetUfaceDomains(var_extra.rotUneighbor(k), var_extra.Xsend(1), var_extra.Xsend_tmp(1),
                                      var_extra.Xneighbor(1), var_extra.Xneighbor_tmp(1), var_extra.request(1), 91);
              }
          }
      }
#endif
  }


  //! black-box matrix vector product for Maxwell equation with SIPG formulation
  template<class T0, class TypeEquation, class Complexe>
  void MltAdd_ElementHcurlSipg(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                               const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                               int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                               const ElementReference<Dimension3, 2>& Fb)
  {
    const EllipticProblem<TypeEquation>& vars = A.var;
    const Mesh<Dimension3>& mesh = vars.mesh;
    const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
    
    // local E_h
    int nb_dof_loc = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    
    Vector<Complexe> Eloc(nb_dof_loc), ProdE(nb_dof_loc), Etemp;
    if (Fb.UseQuadraturePointsForRh())
      Etemp.Reallocate(3*nb_points_quad);
    else
      Etemp.Reallocate(nb_dof_loc);
    
    int offset_elt = mesh_num.Element(i).GetNumberDof(0);
    ProdE.Fill(0);
    for (int j = 0; j < nb_dof_loc; j++)
      Eloc(j) = B(offset_elt + j);

    TinyMatrix<Complexe, Symmetric, 3, 3> Bmass, Astiff;
    TinyVector<Complexe,3> vec_u, vec_v, vec_E, vec_H;
    bool null_stiff = false;
    if (abs(nat_mat.GetCoefStiffness()) == Real_wp(0))
      null_stiff = true;
    
    bool affine = vars.mesh.IsElementAffine(i);
    bool variable = vars.UseNumericalIntegration(i);
    Complexe cone; SetComplexOne(cone);
    
    Complexe x_test;
    ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra
      = const_cast<FemMatrixFreeClass<Complexe, TypeEquation>& >(A)
      .GetExtrapolVariables(x_test);
    
    int offset_face = mesh_num.OffsetQuadElementNumber(i);
    int ref, cond, offset_face2; Complexe phase, phase_conj; bool new_face;
    int num_loc2, num_elem2, rot = 0; Real_wp dsj;
    Complexe coef_stiff = alpha*nat_mat.GetCoefStiffness();

    Vector<Complexe>& Equad = var_extra.GetEnQuad(i, Fb);
    Vector<Complexe>& rotEquad = var_extra.GetCurlEnQuad(i, Fb);
        
    // if true, Equad and rotEquad already contains the correct values
    if (!Fb.UseQuadraturePointsForSh())
      {
        Equad.Reallocate(3*nb_points_quad);
        rotEquad.Reallocate(3*nb_points_quad);
        
        // computes E on quadrature points
        Fb.ApplyChTranspose(Eloc, Equad);
        // curl(E) on quadrature points
        if (!null_stiff)
          {
            if (Fb.UseQuadraturePointsForRh())
              Fb.ApplyRhQuadratureTranspose(Equad, rotEquad);
            else
              Fb.ApplyRhTranspose(Eloc, rotEquad);
          }
      }
	      
    // geometry terms
    int ref_d = vars.mesh.Element(i).GetReference();
    if (vars.UseNumericalIntegration(i))
      for (int j = 0; j < nb_points_quad; j++)
	{
	  TypeEquation::GetMassMatrix(vars, i, j, nat_mat, ref_d, Bmass, 
                                      variable, affine, Fb);
	  
	  vec_u(0) = Equad(3*j);
	  vec_u(1) = Equad(3*j+1);
	  vec_u(2) = Equad(3*j+2);
	  
	  Mlt(Bmass, vec_u, vec_v);
	  
	  Equad(3*j) = alpha*vec_v(0);
	  Equad(3*j+1) = alpha*vec_v(1);
	  Equad(3*j+2) = alpha*vec_v(2);
	  
	  if (!null_stiff)
	    {	
	      TypeEquation::GetStiffMatrix(vars, i, j, nat_mat, ref_d, Astiff);
	      
	      vec_u(0) = rotEquad(3*j);
	      vec_u(1) = rotEquad(3*j+1);
	      vec_u(2) = rotEquad(3*j+2);
	      
	      Mlt(Astiff, vec_u, vec_v);
	      
	      rotEquad(3*j) = alpha*vec_v(0);
	      rotEquad(3*j+1) = alpha*vec_v(1);
	      rotEquad(3*j+2) = alpha*vec_v(2);
	    }
	}
    else
      {
	TypeEquation::GetMassMatrix(vars, i, 0, nat_mat, ref_d, Bmass,
                                    variable, affine, Fb);
        
	TypeEquation::GetStiffMatrix(vars, i, 0, nat_mat, ref_d, Astiff);
	
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    vec_u(0) = Equad(3*j);
	    vec_u(1) = Equad(3*j+1);
	    vec_u(2) = Equad(3*j+2);
	    
	    Mlt(Bmass, vec_u, vec_v);
	    Mlt(alpha*Fb.WeightsND(j), vec_v);
	    
	    Equad(3*j) = vec_v(0);
	    Equad(3*j+1) = vec_v(1);
	    Equad(3*j+2) = vec_v(2);
	    
	    if (!null_stiff)
	      {	
		vec_u(0) = rotEquad(3*j);
		vec_u(1) = rotEquad(3*j+1);
		vec_u(2) = rotEquad(3*j+2);
		
		Mlt(Astiff, vec_u, vec_v);
		Mlt(alpha*Fb.WeightsND(j), vec_v);
		
		rotEquad(3*j) = vec_v(0);
		rotEquad(3*j+1) = vec_v(1);
		rotEquad(3*j+2) = vec_v(2);
	      }
	  }
      }
    
    // treating fluxes terms
    int num_point = nb_points_quad;
    TinyVector<Complexe, 3> Eint, Eext, rotEint, rotEext, saut_E;
    if (!null_stiff)
      for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
        {
          int num_face = mesh.Element(i).numBoundary(num_loc);
          int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
          if (vars.FaceHasToBeConsideredForBoundaryIntegral(num_face))
            {            
              // the face is involved in the fluxes terms
              int rf = mesh_num.GetOrderQuadrature(num_face);        
              var_extra.GetFaceOrientationAndNumber(vars, mesh, num_face, i, num_loc, offset_face,
                                                    ref, cond, new_face, num_elem2, num_loc2,
                                                    offset_face2, rot, phase, phase_conj);
              
              const VectReal_wp& PoidsFlux
                = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
              
              const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
                GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
              
              bool face_on_gamma = true;
              int offset_neighbor = 0;
              bool neighbor_face = false;
              if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
                {
                  face_on_gamma = false;
                  neighbor_face = true;
                  offset_neighbor = vars.GetOffsetNeighboringFace(num_face);
                }
              else
                {
                  if (cond == BoundaryConditionEnum::LINE_INSIDE)
                    face_on_gamma = false;
                }
              
              if (!face_on_gamma)
                {
                  Complexe coef_penal = vars.alpha_penalization;
                  if (vars.automatic_choice_penalization)
                    coef_penal *= vars.Glob_CoefPenalDG(num_face);
                  
                  Vector<Complexe> feval_E(3*nb_points_face), feval_H(3*nb_points_face);
                  for (int j = 0; j < nb_points_face; j++)
                    {
                      int jrot = FacesQuadRotation(rot, j);
                      if (new_face)
                        dsj = vars.Glob_dsj(num_face)(j);
                      else
                        dsj = vars.Glob_dsj(num_face)(jrot);
                      
                      Real_wp poids = dsj*PoidsFlux(j);                  
                      Eint(0) = var_extra.extrapolU(0)(offset_face+j);
                      Eint(1) = var_extra.extrapolU(1)(offset_face+j);
                      Eint(2) = var_extra.extrapolU(2)(offset_face+j);
                      
                      rotEint(0) = var_extra.extrapolRotU(0)(offset_face+j);
                      rotEint(1) = var_extra.extrapolRotU(1)(offset_face+j);
                      rotEint(2) = var_extra.extrapolRotU(2)(offset_face+j);
                      
                      if (neighbor_face)
                        {
                          Eext(0) = var_extra.Uneighbor(0)(offset_neighbor+j);
                          Eext(1) = var_extra.Uneighbor(1)(offset_neighbor+j);
                          Eext(2) = var_extra.Uneighbor(2)(offset_neighbor+j);

                          rotEext(0) = var_extra.rotUneighbor(0)(offset_neighbor+j);
                          rotEext(1) = var_extra.rotUneighbor(1)(offset_neighbor+j);
                          rotEext(2) = var_extra.rotUneighbor(2)(offset_neighbor+j);
                        }
                      else
                        {
                          Eext(0) = var_extra.extrapolU(0)(offset_face2+jrot);
                          Eext(1) = var_extra.extrapolU(1)(offset_face2+jrot);
                          Eext(2) = var_extra.extrapolU(2)(offset_face2+jrot);
                          
                          rotEext(0) = var_extra.extrapolRotU(0)(offset_face2+jrot);
                          rotEext(1) = var_extra.extrapolRotU(1)(offset_face2+jrot);
                          rotEext(2) = var_extra.extrapolRotU(2)(offset_face2+jrot);
                        }
                      
                      // part {1/mu rot U} n \times phi
                      Add(rotEint, rotEext, vec_u);
                      Add(Eint, Eext, saut_E);
                      Add(coef_penal, saut_E, vec_u);
                      
                      int j2 = num_point + j - nb_points_quad;
                      if (affine)
                        j2 = 0;
                      
                      MltTrans(vars.Glob_DFj(i)(j2), vec_u, vec_v);
                      TimesProd(vec_v, Fb.NormaleLoc(num_loc), vec_H);
                      vec_H *= -1.0/dsj*coef_stiff*poids;
                      CopyVector(vec_H, j, feval_H);
                                            
                      // part n \times (u^+ - u^-) \cdot 1/mu rot phi
                      MltTrans(vars.Glob_DFj(i)(j2), saut_E, vec_v);
                      vars.ref_invMu(ref_d).MltMatrix(vars, i, num_point+j, vec_v, vec_E);
                      vec_E *= -vars.Glob_invJacobian(i)(j2)*coef_stiff*poids;
                      CopyVector(vec_E, j, feval_E);
                    }
                  
                  if (Fb.UseQuadraturePointsForSh())
                    {
                      Fb.ApplyShQuadrature(1.0, num_loc, feval_H, Equad, rf);
                      if (Fb.UseQuadraturePointsForRh())
                        Fb.ApplyShQuadrature(1.0, num_loc, feval_E, rotEquad, rf);
                      else
                        Fb.ApplyNablaShQuadrature(1.0, num_loc, feval_E, Equad, rf);
                    }
                  else
                    {
                      Fb.ApplySh(1.0, num_loc, feval_H, ProdE, rf);
                      Fb.ApplyNablaSh(1.0, num_loc, feval_E, ProdE, rf);
                    }
                }
            }
          
          num_point += nb_points_face;
          offset_face += nb_points_face;
        }    
    
    // integration against curl(phi)
    if (!null_stiff)
      {
	if (Fb.UseQuadraturePointsForRh())
	  {              
            Fb.ApplyRhQuadrature(rotEquad, Etemp);
	    Add(cone, Etemp, Equad);
            if (Fb.UseQuadraturePointsForSh())
              Fb.ApplyCh(Equad, ProdE);
            else
              {
                Fb.ApplyCh(Equad, Eloc);
                Add(cone, Eloc, ProdE);
              }
	  }
	else
	  {
            Fb.ApplyRh(rotEquad, Etemp);
            Add(cone, Etemp, ProdE);
            Fb.ApplyCh(Equad, Etemp);
            Add(cone, Etemp, ProdE);
          }
      }

    // adding contributions to the global vector C
    for (int j = 0; j < nb_dof_loc; j++)
      C(offset_elt+j) += ProdE(j);    
  }

  
  //! generic matrix-vector product for SIPG formulation
  //! and Maxwell equations
  /*
    this algorithm separates the geometry \f$ DF_i^{*-1} \f$ from the 
    spatial derivatives. This allows low-storage and fast-algorithm
    in \f$ O(r^4) \f$ instead of \f$ O(r^6) \f$ if the full matrix is stored, where
    r is the order of approximation
  */
  template <class T0, class Complexe, class TypeEquation, class Vector1>
  void MltAdd_SipgMaxwell(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                          const SeldonTranspose& trans, int level,
                          const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                          const Vector1& B2, const T0& beta, Vector1& C2, bool assemble)
  {
#ifdef SELDON_WITH_MPI
    Vector1 B(B2), C(C2);
#else
    Vector1 B, C;
#endif
    A.ApplyRightScaling(B2, C2, B, C);    

    if (beta == T0(0))
      C2.Fill(0);
    else
      Mlt(beta, C2);
    
    ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra =
      const_cast<FemMatrixFreeClass<Complexe, TypeEquation>& >(A).GetExtrapolVariables();
    
    var_extra.Reallocate(A.var.mesh.GetNbElt());
    
    // computes E and rot E on quadrature points of surfaces
    GetExtrapolationAndExchange_Sipg(nat_mat, A, B, var_extra);
        
    // loop over all elements of the mesh
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      {
        MltAdd_ElementHcurlSipg(alpha, nat_mat, A, i, B, C, A.var.GetReferenceElementHcurl(i));
      }
    
    A.AddExtraBoundaryTerms(Complexe(alpha), B, C);
    
    A.ApplyLeftScaling(B2, C2, B, C);    
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C2);
  }
    
} // namespace Montjoie

#define MONTJOIE_FILE_PROD_MAT_VECT_HCURL_3D_CXX
#endif
