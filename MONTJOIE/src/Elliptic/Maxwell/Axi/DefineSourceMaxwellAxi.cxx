#ifndef MONTJOIE_FILE_DEFINE_SOURCE_MAXWELL_AXI_CXX

namespace Montjoie
{

  /***********************
   * PlaneWaveMaxwellAxi *
   ***********************/
  

  void PlaneWaveMaxwellAxi::InitDefaultValues()
  {
#ifdef SELDON_WITH_MPI   
    int rank_proc; MPI_Comm_rank(var_maxwell.comm_group_mode, &rank_proc);
#endif

    epsilon0 = var_maxwell.epsilon0;
    mu0 = var_maxwell.mu0; 
    number_mode = var_maxwell.GetCurrentModeNumber();

    R3 k_wave = var_maxwell.GetWaveVector(); //DISP(k_wave);
    incidence_axial = true;
    // wave vector -> (k_r, k_theta, k_z) where k_z = kwave_z and (k_r,k_theta) = k_wave
    // (k_r, k_theta) is computed in form k_bot (cos theta, sin theta)
    k_bot = sqrt(k_wave(0)*k_wave(0) + k_wave(1)*k_wave(1));
    if ((k_bot/var_maxwell.GetOmega()) > 10.0*epsilon_machine)
      {
	incidence_axial = false; 
	teta0 = acos(k_wave(0)/k_bot);
	if (k_wave(1) < 0)
	  teta0 = -teta0;
      }
    
    kwave_z = k_wave(2); omega = var_maxwell.GetOmega();
    
    // polarization, eps and mu
    R3_Complex_wp polar; var_maxwell.GetPolarization(polar);
    //if (rank_proc == 0)
    //{ DISP(k_wave); DISP(polar); }
    
    if (abs(DotProd(k_wave, polar)) > 100*epsilon_machine)
      {
        cout << "Polarization and wave vector are not orthogonal" << endl;
        abort();
      }

    E0_x = polar(0);
    E0_y = polar(1);
    E0_z = polar(2); 
    
    // H^0 = k \times E^0 / || k ||
    H0_x = (k_wave(1)*E0_z - kwave_z*E0_y) / omega;
    H0_y = (kwave_z*E0_x - k_wave(0)*E0_z) / omega;
    H0_z = (k_wave(0)*E0_y - k_wave(1)*E0_x) / omega;
    
    threshold = 10.0*epsilon_machine;

    modified_form = false;
    if (var_maxwell.modified_formulation && (number_mode != 0))
      modified_form = true;
  }
  
  
  //! computation of incident wave
  /*!
    \param[in] x (r,z) coordinates of the point
    \param[out] E components of incident electric field (E_r, E_theta, E_z)
    \param[out] H components of incident megnetic field (H_r, H_theta, H_z)
    \param[in] nb quadrature point number
   */
  void PlaneWaveMaxwellAxi::Evaluate_E_H(const R2& x, R3_Complex_wp& E, R3_Complex_wp& H, int nb)
  {
    Real_wp arg = kwave_z*x(1);
    Complex_wp u_inc = exp(Complex_wp(0.0, arg)), coef_jm, coef_jm1, coef_jp1;
    if (incidence_axial)
      {
	// wave vector oriented along e_z
	if (number_mode == 1)
	  {
	    // E_{r,m} = 1/2 e^{ikz} ( E_{0,x} + I* E_{0,y} )
	    // E_{\theta,m} = 1/2*I* e^{ikz} ( E_{0,y} - I* E_{0,x} )
	    E(0) = u_inc * 0.5 * (E0_x + Iwp*E0_y);
	    E(2) = 0.0;
	    E(1) = Iwp * u_inc * 0.5 * (E0_y - Iwp*E0_x) ;
	    
	    H(0) = u_inc * 0.5 * (H0_x + Iwp*H0_y);
	    H(2) = 0.0;
	    H(1) = Iwp * u_inc * 0.5 * (H0_y - Iwp*H0_x) ;
	  }
	else if (number_mode == -1)
	  {
	    // E_{r,m} = 1/2 e^{ikz} ( E_{0,x} - I* E_{0,y} )
	    // E_{\theta,m} = 1/2*I* e^{ikz} ( E_{0,y} + I* E_{0,x} )
	    E(0) = u_inc * 0.5 *(E0_x - Iwp*E0_y);
	    E(2) = 0.0;
	    E(1) = -Iwp * u_inc * 0.5 * (E0_y + Iwp*E0_x);
	    
	    H(0) = -u_inc * 0.5 * (H0_x - Iwp*H0_y);
	    H(2) = 0.0;
	    H(1) = Iwp * u_inc * 0.5 * (H0_y + Iwp*H0_x);
	  }
	else
	  {
	    E.Fill(0);
	    H.Fill(0);
	  }
      }
    else
      {
	// general case (wave vector not oriented along e_z)
	VectReal_wp Jn(3);
	// we compute Jn(kr), J_{n-1}(kr) and J_{n+1}(kr)
	// These functions are pre-computed
	//Jn(0) = this->vars.GetBessel_Value(number_mode-1, nb);
	//Jn(1) = this->vars.GetBessel_Value(number_mode, nb);
	//Jn(2) = this->vars.GetBessel_Value(number_mode+1, nb);
	Real_wp kr = k_bot*x(0);
	ComputeBesselFunctions(number_mode-1, 3, kr, Jn);
	
	Complex_wp Im1, Im, Ip1; Im1 = ComputePowerI(number_mode-1);
	Im = Im1*Iwp; Ip1 = Im*Iwp;
	Complex_wp exp_m1 = exp(Complex_wp(number_mode-1)*Iwp*teta0);
	Complex_wp exp_m = exp(Complex_wp(number_mode)*Iwp*teta0);
	Complex_wp exp_p1 = exp(Complex_wp(number_mode+1)*Iwp*teta0);
	
	coef_jm1 = u_inc*0.5*Im1*Jn(0)*exp_m1;
	coef_jp1 = u_inc*0.5*Ip1*Jn(2)*exp_p1;
	coef_jm = u_inc*Im*Jn(1)*exp_m;
	E(0) = coef_jm1*(E0_x+Iwp*E0_y) + coef_jp1*(E0_x-Iwp*E0_y);
	E(2) = coef_jm*E0_z;
	E(1) = Iwp*(coef_jm1*(E0_y - Iwp*E0_x) + coef_jp1*(E0_y + Iwp*E0_x) );
	if (number_mode < 0)
	  E(1) = -E(1);
	
	H(0) = coef_jm1*(H0_x+Iwp*H0_y) + coef_jp1*(H0_x-Iwp*H0_y);
	H(2) = coef_jm*H0_z;
	H(1) = Iwp*(coef_jm1*(H0_y - Iwp*H0_x) + coef_jp1*(H0_y + Iwp*H0_x) );
	if (number_mode < 0)
	  {
	    H(0) = -H(0);
	    H(2) = -H(2);
	  }
      }
  }
  
  
  //! computation of incident electric field
  /*!
    \param[in] x (r,z) coordinates of the point
    \param[out] f components of incident electric field (E_r, E_theta, E_z)
  */
  void PlaneWaveMaxwellAxi::Evaluate_E(const R2& x, R3_Complex_wp& f)
  {
    f.Fill(0);
    Real_wp arg = kwave_z*x(1);
    Complex_wp u_inc = exp(Complex_wp(0.0,arg));
    if (incidence_axial)
      {
	if (number_mode == 1)
	  {
	    // E_{r,m} = e^{ikz} ( 1/2 E_{0,x} + i/2*I* E_{0,y} )
	    f(0) = u_inc * 0.5 * (E0_x + Iwp*E0_y);
	    f(1) = Iwp * u_inc * 0.5 * (E0_y - Iwp*E0_x) ;
	    f(2) = 0.0;
	  }
	else if (number_mode == -1)
	  {
	    // E_{r,m} = e^{ikz} ( 1/2 E_{0,x} - i/2*I* E_{0,y} )
	    f(0) = u_inc * 0.5 *(E0_x - Iwp*E0_y);
	    f(1) = -Iwp * u_inc * 0.5 * (E0_y + Iwp*E0_x) ;
	    f(2) = 0.0;
	  }
	else
	  f.Fill(0);
      }
    else
      {
	// general case, wave vector not oriented along e_z
	VectReal_wp Jn; Real_wp kr = k_bot*x(0);
	ComputeBesselFunctions(number_mode-1, 3, kr, Jn);
	// DISP(kr); DISP(Jn); DISP(k_bot); DISP(x);
	Complex_wp Im1, Im, Ip1; Im1 = ComputePowerI(number_mode-1);
	Im = Im1*Iwp; Ip1 = Im*Iwp; // DISP(number_mode); DISP(Im1); DISP(Im); DISP(Ip1);
	Complex_wp exp_m1 = exp(Complex_wp(number_mode-1)*Iwp*teta0);
	Complex_wp exp_m = exp(Complex_wp(number_mode)*Iwp*teta0);
	Complex_wp exp_p1 = exp(Complex_wp(number_mode+1)*Iwp*teta0);
	f(0) = u_inc*0.5* ( Im1*Jn(0)*exp_m1*(E0_x+Iwp*E0_y) + Ip1*Jn(2)*exp_p1*(E0_x-Iwp*E0_y) );
	f(1) = Iwp*u_inc*0.5* ( Im1*Jn(0)*exp_m1*(E0_y - Iwp*E0_x)
                                + Ip1*Jn(2)*exp_p1*(E0_y + Iwp*E0_x) );
	
        f(2) = u_inc*Im*Jn(1)*E0_z*exp_m;
	if (number_mode < 0)
	  f(1) = -f(1);

	// DISP(u_inc); DISP(E0_x); DISP(E0_y); DISP(f);
      }
  }
  
  
  //! computation of incident magnetic field
  /*!
    \param[in] x (r,z) coordinates of the point
    \param[out] f components of incident magnetic field (H_r, H_theta, H_z)
  */
  void PlaneWaveMaxwellAxi::Evaluate_H(const R2& x, R3_Complex_wp& f)
  {
    f.Fill(0);
    Real_wp arg = kwave_z*x(1);
    Complex_wp u_inc = exp(Complex_wp(0.0,arg));
    if (incidence_axial)
      {
	if (number_mode == 1)
	  {
	    // H_{r,m} = e^{ikz} ( 1/2 H_{0,x} + i/2*I* H_{0,y} )
	    f(0) = u_inc * 0.5 * (H0_x + Iwp*H0_y);
	    f(1) = Iwp * u_inc * 0.5 * (H0_y - Iwp*H0_x) ;
	    f(2) = 0.0;
	  }
	else if (number_mode == -1)
	  {
	    // H_{r,m} = e^{ikz} ( 1/2 H_{0,x} - i/2*I* H_{0,y} )
	    f(0) = -u_inc * 0.5 *(H0_x - Iwp*H0_y);
	    f(1) = Iwp * u_inc * 0.5 * (H0_y + Iwp*H0_x) ;
	    f(2) = 0.0;
	  }
	else
	  f.Fill(0);
      }
    else
      {
	// general case, wave vector not oriented along e_z
	VectReal_wp Jn; Real_wp kr = k_bot*x(0); 
	ComputeBesselFunctions(number_mode-1, 3, kr, Jn);
	// DISP(kr); DISP(Jn); DISP(k_bot); DISP(x);
	Complex_wp Im1,Im,Ip1; Im1 = ComputePowerI(number_mode-1);
	Im = Im1*Iwp; Ip1 = Im*Iwp; // DISP(number_mode); DISP(Im1); DISP(Im); DISP(Ip1);
	Complex_wp exp_m1 = exp(Complex_wp(number_mode-1)*Iwp*teta0);
	Complex_wp exp_m = exp(Complex_wp(number_mode)*Iwp*teta0);
	Complex_wp exp_p1 = exp(Complex_wp(number_mode+1)*Iwp*teta0);
	f(0) = u_inc*0.5* ( Im1*Jn(0)*exp_m1*(H0_x+Iwp*H0_y) + Ip1*Jn(2)*exp_p1*(H0_x-Iwp*H0_y) );
	f(1) = Iwp*u_inc*0.5* ( Im1*Jn(0)*exp_m1*(H0_y - Iwp*H0_x) 
                                + Ip1*Jn(2)*exp_p1*(H0_y + Iwp*H0_x) );
	f(2) = u_inc*Im*Jn(1)*H0_z*exp_m;

	if (number_mode < 0)
	  {
	    f(0) = -f(0);
	    f(2) = -f(2);
	  }

	// DISP(u_inc); DISP(H0_x); DISP(H0_y); DISP(f);
      }
  }


  void PlaneWaveMaxwellAxi::EvaluateIncidentField(int i, int j, const R2& x, R3_Complex_wp& f)
  {
    R3_Complex_wp E, H;
    if (modified_form)
      Evaluate_E_H(x, E, H, 0);
    else
      Evaluate_E(x, f);

    /* if (var_maxwell.InsidePML(i))
      {
      int i1 = i - var_maxwell.mesh.GetNbElt() + var_maxwell.GetNbEltPML();
	Complex_wp dr = 1.0/var_maxwell.GetTauPML(i1, j, 0);
	Complex_wp dz = 1.0/var_maxwell.GetTauPML(i1, j, 1);
        if (modified_form)
	  {
	    E(0) *= dr; H(0) *= dr;
	    E(2) *= dz; H(2) *= dz;
	  }
	else
	  {
	    f(0) *= dr;
	    f(2) *= dz;
	  }
	  } */
    
    if (modified_form)
      {
	// unknown U is equal to omega mu \tilde{H}
	f(0) = omega*H(2);
	f(2) = -omega*H(0);
	
	f(1) = E(1);
      }
  }


  void IncidentWaveProjector_MaxwellAxi
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    R3_Complex_wp fvec;
    this->EvaluateIncidentField(i, j, x, fvec);
    f(0) = fvec(0); f(1) = fvec(2); f(2) = fvec(1);
  }
  
  
  /***********************************
   * DiffractedWaveSource_MaxwellAxi *
   ***********************************/
  
  
  //! initialization before volumetric integration
  void DiffractedWaveSource_MaxwellAxi::InitElement(int num_elem, const VectR2& s)
  {
    pml_element = var_maxwell.InsidePML(num_elem);
    int ref = var_maxwell.mesh.Element(num_elem).GetReference();
    ref_domain = ref;
    coef_E.SetDiagonal(-epsilon0 + Iwp*var_maxwell.ref_sigma(ref).GetConstant()/omega);
    coef_E += var_maxwell.ref_epsilon(ref).GetConstant();
    coef_E *= omega;
    
    if (var_maxwell.modified_formulation)
      {
	coef_H.SetDiagonal(-1.0/mu0);
	coef_H += var_maxwell.ref_invMu(ref).GetConstant();
      }
    else
      {
	coef_H.SetDiagonal(-mu0);
	coef_H += var_maxwell.ref_mu(ref).GetConstant();
      }
    
    coef_H *= omega;

    int nb_pts_quad = var_maxwell.GetNbPointsQuadratureInside(num_elem);
    EvalF.Reallocate(nb_pts_quad);
    EvalH.Reallocate(nb_pts_quad);
  }
  

  //! evaluation of incident electric field
  void DiffractedWaveSource_MaxwellAxi
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    if (this->num_loc_ >= 0)
      {
	R3_Complex_wp fvec;
	this->EvaluateIncidentField(i, j, x, fvec);

	// Dirichlet condition is null for boundaries of the PML layer
	if (var_maxwell.InsidePML(i))
	  f.Fill(0);
	else
	  {
	    f(0) = -fvec(0);
	    f(1) = -fvec(2);
	    f(2) = -fvec(1);
	  }
      }
  }
  
  
  //! returns true if volumetric source
  bool DiffractedWaveSource_MaxwellAxi
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    if (var_maxwell.ref_epsilon(ref_domain).IsVarying())
      return true;

    if (var_maxwell.ref_sigma(ref_domain).IsVarying())
      return true;

    if (var_maxwell.ref_mu(ref_domain).IsVarying())
      return true;
    
    Real_wp norme_rhs = MaxAbs(coef_E);
    if (norme_rhs > threshold)
      return true;

    if (MaxAbs(coef_H) > threshold)
      return true;
    
    return false;
  }
  
  
  //! evaluation of volumetric source
  void DiffractedWaveSource_MaxwellAxi
  ::EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    R3_Complex_wp E, H;
    this->Evaluate_E_H(x, E, H, 0);
    
    if (var_maxwell.ref_epsilon(ref_domain).IsVarying() ||
	var_maxwell.ref_sigma(ref_domain).IsVarying() )
      {
	coef_E.SetDiagonal(-epsilon0 + Iwp*var_maxwell.ref_sigma(ref_domain)
			   .GetCoefficient(var_maxwell, i, j)/omega);
	
	coef_E += var_maxwell.ref_epsilon(ref_domain).GetCoefficient(var_maxwell, i, j);
	coef_E *= omega;
      }

    Mlt(coef_E, E, EvalF(j));
    
    if (var_maxwell.ref_mu(ref_domain).IsVarying())
      {
	if (var_maxwell.modified_formulation)
	  {
	    coef_H.SetDiagonal(-1.0/mu0);
	    coef_H += var_maxwell.ref_invMu(ref_domain).GetCoefficient(var_maxwell, i, j);
	  }
	else
	  {
	    coef_H.SetDiagonal(-mu0);
	    coef_H += var_maxwell.ref_mu(ref_domain).GetCoefficient(var_maxwell, i, j);
	  }
	
	coef_H *= omega;
      }

    Mlt(coef_H, H, EvalH(j));

    if (var_maxwell.modified_formulation)
      EvalF(j) *= omega;
    
    // part due to omega^2 (epsilon - eps0) E^inc
    if (this->modified_form)
      {
        f(0) = x(0)*EvalF(j)(0)*abs(this->number_mode);
        f(2) = EvalF(j)(1)*square(this->number_mode) + EvalF(j)(0)*abs(this->number_mode);
        f(1) = x(0)*EvalF(j)(2)*abs(this->number_mode);
	
	f(0) *= x(0); f(1) *= x(0); f(2) *= x(0);
      }
    else
      {
	f(0) = x(0)*EvalF(j)(0);  
	f(2) = x(0)*EvalF(j)(1);  
	f(1) = x(0)*EvalF(j)(2);  
      }

    // part due to (1/mu - 1/mu0) H^inc
    if (var_maxwell.modified_formulation)
      {
	if (this->number_mode == 0)
	  f(2) += EvalH(j)(2);
	else
	  {
	    Real_wp m2 = square(this->number_mode);
	    f(1) += x(0)*(EvalH(j)(0)*m2 - EvalH(j)(1));
	    f(0) -= x(0)*EvalH(j)(2)*m2;
	  }
      }
  }

  
  //! returns true if volumetric source against curl of basis functions
  bool DiffractedWaveSource_MaxwellAxi
  ::IsNonNullGradientSource(const VectR2& s)
  {
    if (this->modified_form)
      {
	if (var_maxwell.ref_epsilon(ref_domain).IsVarying())
	  return true;
	
	if (var_maxwell.ref_sigma(ref_domain).IsVarying())
	  return true;
	
	if (MaxAbs(coef_E) > threshold)
	  return true;
      }
    
    if (var_maxwell.ref_mu(ref_domain).IsVarying())
      return true;
	
    if (MaxAbs(coef_H) > threshold)
      return true;
    
    return false;
  }
    
  
  //! evalues f in the source term \int f curl(phi)
  void DiffractedWaveSource_MaxwellAxi
  ::EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    // part due to omega^2 (epsilon - epsilon0) E^{inc}
    if (this->modified_form)
      {
	f(1) = -x(0)*x(0)*EvalF(j)(0)*abs(this->number_mode);
	f(0) = x(0)*x(0)*EvalF(j)(2)*abs(this->number_mode);
      }
    
    // part due to -(1/mu - 1/mu0) H^inc
    if (this->number_mode == 0)
      {
	f(0) = -x(0)*EvalH(j)(0);
	f(2) = x(0)*EvalH(j)(1);
	f(1) = -x(0)*EvalH(j)(2);
      }
    else
      {
	f(2) = -x(0)*x(0)*EvalH(j)(1);
      }
  }
  
  
  //! returns true if there is source on volumetric dofs (Hrz, Hteta)
  bool DiffractedWaveSource_MaxwellAxi
  ::IsNonNull_SourceH(const VectR2& s)
  {
    if (var_maxwell.modified_formulation)
      return false;

    if (var_maxwell.ref_mu(ref_domain).IsVarying())
      return true;
	
    if (MaxAbs(coef_H) > threshold)
      return true;
    
    return false;
  }
  
  
  //! evaluates source for volumetric dofs (Hrz, Hteta)
  void DiffractedWaveSource_MaxwellAxi
  ::EvaluateSource_H(int i, int j, const R2& x, VectComplex_wp& f)
  {
    f(0) = EvalH(j)(0);
    f(2) = EvalH(j)(1);
    f(1) = EvalH(j)(2);
  }  
  
  
  //! surfacic source for Neumann boundary condition
  bool DiffractedWaveSource_MaxwellAxi
  ::IsNonNullSurfacicSource(int ref)
  {
    int ref_cond = var_maxwell.mesh.GetBoundaryCondition(ref);
    if ((ref_cond == BoundaryConditionEnum::LINE_NEUMANN)||(ref_cond == BoundaryConditionEnum::LINE_IMPEDANCE))
      return true;
    
    return false;
  }
  
  //! evaluation of surfacic source
  void DiffractedWaveSource_MaxwellAxi
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    f.Fill(0);
    int n1 = var_maxwell.mesh.Boundary(this->num_boundary_).numVertex(0);
    int n2 = var_maxwell.mesh.Boundary(this->num_boundary_).numVertex(1);
    if (var_maxwell.IsVertexOnAxis(n1) && var_maxwell.IsVertexOnAxis(n2))
      return;
    
    // Neumann condition -\int r H^{inc} \times n \lambda 
    R2 pt_glob = PointsElem.GetPointQuadratureBoundary(k);
    R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Real_wp radius = pt_glob(0);
    R3_Complex_wp E, H;
    this->Evaluate_E_H(pt_glob, E, H, 0);
    
    if (k == 0)
      {
	int Nquad = MatricesElem.GetNbPointsQuadratureBoundary();
	this->EvalH.Reallocate(Nquad);
      }
    
    this->EvalH(k)(0) = -H(1)*normale(1);
    this->EvalH(k)(2) = H(1)*normale(0);
    this->EvalH(k)(1) = normale(0)*H(2) - normale(1)*H(0);
    if (var_maxwell.modified_formulation)
      this->EvalH(k) *= omega;
   
    if (this->modified_form)
      {
	f(0) = radius*this->EvalH(k)(0)*abs(this->number_mode);
        f(2) = this->EvalH(k)(1)*square(this->number_mode)
	  + this->EvalH(k)(0)*abs(this->number_mode);
	
        f(1) = radius*this->EvalH(k)(2)*abs(this->number_mode);
      }
    else
      {
	f(0) = this->EvalH(k)(0);
	f(2) = this->EvalH(k)(1);
	f(1) = this->EvalH(k)(2);
      }

    f(0) *= radius; f(1) *= radius; f(2) *= radius;
  }

  
  bool DiffractedWaveSource_MaxwellAxi
  ::IsNonNullSurfacicSourceGradient(int ref)
  {
    int ref_cond = var_maxwell.mesh.GetBoundaryCondition(ref);
    if ((ref_cond == BoundaryConditionEnum::LINE_NEUMANN)||(ref_cond == BoundaryConditionEnum::LINE_IMPEDANCE))
      return this->modified_form;
    
    return false;
  }


  void DiffractedWaveSource_MaxwellAxi
  ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension2>& PointsElem,
				   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    f(1) = -x(0)*x(0)*this->EvalH(k)(0)*abs(this->number_mode);
    f(0) = x(0)*x(0)*this->EvalH(k)(2)*abs(this->number_mode);    
  }
  

  /******************************
   * TotalWaveSource_MaxwellAxi *
   ******************************/


  //! Dirichlet condition for total wave field source
  void TotalWaveSource_MaxwellAxi
  ::EvaluateFunction(int i, int j, const R2& x, VectComplex_wp& f)
  {
    if (this->num_loc_ >= 0)
      {
	R3_Complex_wp fvec;
	this->EvaluateIncidentField(i, j, x, fvec);

	// Dirichlet condition is null for all boundaries except PML outer boundary
	if (!var_maxwell.InsidePML(i))
	  f.Fill(0);
	else
	  {
            f(0) = fvec(0);
            f(2) = fvec(1);
            f(1) = fvec(2);
          }
      }
  }
  
  
  //! source for a PML element
  void TotalWaveSource_MaxwellAxi
  ::EvaluateVolumetricSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    cout << "PML not implemented with total source" << endl;
    abort();
    
    /*
    R3_Complex_wp E, H;
    this->Evaluate_E_H(x, E, H, 0);
    
    Complex_wp eps_teta, mu_teta;
    TinyMatrix<Complex_wp, Symmetric, 2, 2> eps_rz, mu_rz;
    this->vars.GetPhysicalCoefPML(i, x, eps_rz, eps_teta, mu_rz, mu_teta);
    
    eps_rz(0, 0) -= this->epsilon0; eps_rz(1, 1) -= this->epsilon0; eps_teta -= this->epsilon0;
    if (this->vars.modified_formulation)
      {
	GetInverse(mu_rz); mu_teta = 1.0/mu_teta;
	mu_rz(0, 0) -= 1.0/this->mu0; mu_rz(1, 1) -= 1.0/this->mu0; mu_teta -= 1.0/this->mu0;
      }
    else
      {
	mu_rz(0, 0) -= this->mu0; mu_rz(1, 1) -= this->mu0; mu_teta -= this->mu0;
      }
    
    this->EvalF(j)(0) = eps_rz(0, 0)*E(0) + eps_rz(0, 1)*E(2);
    this->EvalF(j)(1) = eps_teta*E(1);
    this->EvalF(j)(2) = eps_rz(1, 0)*E(0) + eps_rz(1, 1)*E(2);

    this->EvalH(j)(0) = mu_rz(0, 0)*H(0) + mu_rz(0, 1)*H(2);
    this->EvalH(j)(1) = mu_teta*H(1);
    this->EvalH(j)(2) = mu_rz(1, 0)*H(0) + mu_rz(1, 1)*H(2);

    this->EvalH(j) *= -this->omega;    
    if (this->vars.modified_formulation)
      this->EvalF(j) *= -this->omega*this->omega;
    else
      this->EvalF(j) *= -this->omega;
    
    if (this->modified_form)
      {
        f(0) = x(0)*this->EvalF(j)(0)*abs(this->number_mode);
        f(1) = this->EvalF(j)(1)*square(this->number_mode)
	  + this->EvalF(j)(0)*abs(this->number_mode);
        
	f(2) = x(0)*this->EvalF(j)(2)*abs(this->number_mode);
	
	f *= x(0);
      }
    else
      f = x(0)*this->EvalF(j);  
    
    if (this->vars.modified_formulation)
      {
	if (this->number_mode == 0)
	  f(1) += this->EvalH(j)(2);
	else
	  {
	    Real_wp m2 = square(this->number_mode);
	    f(2) += x(0)*(this->EvalH(j)(0)*m2 - this->EvalH(j)(1));
	    f(0) -= x(0)*this->EvalH(j)(2)*m2;
	  }
	  }    */
  }
  
  
  //! returns true for a PML element
  bool TotalWaveSource_MaxwellAxi
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    return this->pml_element;
  }
  
  
  void TotalWaveSource_MaxwellAxi
  ::EvaluateSource_H(int i, int j, const R2& x, VectComplex_wp& f)
  {
    f(0) = this->EvalH(j)(0);
    f(2) = this->EvalH(j)(1);
    f(1) = this->EvalH(j)(2);
  }
  
  
  bool TotalWaveSource_MaxwellAxi
  ::IsNonNull_SourceH(const VectR2& s)
  {
    if (var_maxwell.modified_formulation)
      return false;
    
    return this->pml_element;
  }
  
  
  bool TotalWaveSource_MaxwellAxi
  ::IsNonNullGradientSource(const VectR2& s)
  {
    if (!var_maxwell.modified_formulation)
      return false;
    
    return this->pml_element;
  }
    
  
  void TotalWaveSource_MaxwellAxi
  ::EvaluateGradientSource(int i, int j, const R2& x, VectComplex_wp& f)
  {
    f.Fill(0);
    
    // part due to omega^2 (epsilon - epsilon0) E^{inc}
    if (this->modified_form)
      {
	f(1) = -x(0)*x(0)*this->EvalF(j)(0)*abs(this->number_mode);
	f(0) = x(0)*x(0)*this->EvalF(j)(2)*abs(this->number_mode);
      }
    
    // part due to -(1/mu - 1/mu0) H^inc
    if (this->number_mode == 0)
      {
	f(0) = -x(0)*this->EvalH(j)(0);
	f(2) = x(0)*this->EvalH(j)(1);
	f(1) = -x(0)*this->EvalH(j)(2);
      }
    else
      {
	f(2) = -x(0)*x(0)*this->EvalH(j)(1);
      }
  }
  

  //! surfacic source for first order absorbing boundary condition
  bool TotalWaveSource_MaxwellAxi
  ::IsNonNullSurfacicSource(int ref)
  {
    int ref_cond = var_maxwell.mesh.GetBoundaryCondition(ref);
    if (ref_cond == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    return false;
  }
  
  
  //! evaluation of surfacic source
  void TotalWaveSource_MaxwellAxi
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension2>& PointsElem,
			   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    // source for E^\theta is equal to H^{inc} \times n + \ii E^{inc}_\theta
    R2 pt_glob = PointsElem.GetPointQuadratureBoundary(k);
    R2 normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    R3_Complex_wp E, H;
    
    this->Evaluate_E_H(pt_glob, E, H, 0);
    
    if (k == 0)
      {
	int Nquad = MatricesElem.GetNbPointsQuadratureBoundary();
	this->EvalF.Reallocate(Nquad);
      }

    this->EvalF(k)(0) = -normale(1)*(-H(1) + Iwp*(E(0)*normale(1)-E(2)*normale(0)));
    this->EvalF(k)(1) = (H(0)*normale(1)-H(2)*normale(0)) - Iwp*E(1);
    this->EvalF(k)(2) = normale(0)*(-H(1) + Iwp*(E(0)*normale(1)-E(2)*normale(0)));
    
    if (this->modified_form)
      {
	Real_wp radius = pt_glob(0);
        f(0) = radius*this->EvalF(k)(0)*abs(this->number_mode);
        f(2) = this->EvalF(k)(1)*square(this->number_mode)
	  + this->EvalF(k)(0)*abs(this->number_mode);
	
        f(1) = radius*this->EvalF(k)(2)*abs(this->number_mode);

      }
    else
      {
        f(0) = this->EvalF(k)(0);
        f(2) = this->EvalF(k)(1);
        f(1) = this->EvalF(k)(2);
      }

    if (var_maxwell.modified_formulation)
      Mlt(this->omega*pt_glob(0), f);
    else
      Mlt(pt_glob(0), f);
  }


  bool TotalWaveSource_MaxwellAxi
  ::IsNonNullSurfacicSourceGradient(int ref)
  {
    int ref_cond = var_maxwell.mesh.GetBoundaryCondition(ref);
    if (ref_cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	if (this->modified_form)
	  return true;
      }
    
    return false;
    
  }
  

  void TotalWaveSource_MaxwellAxi
  ::EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension2>& PointsElem,
				   const SetMatrices<Dimension2>& MatricesElem, VectComplex_wp& f)
  {
    R2 x = PointsElem.GetPointQuadratureBoundary(k);
    f.Fill(0);
    f(1) = -x(0)*x(0)*this->EvalF(k)(0)*abs(this->number_mode)*this->omega;
    f(0) = x(0)*x(0)*this->EvalF(k)(2)*abs(this->number_mode)*this->omega;    
  }


  /*******************************
   * VolumetricSource_MaxwellAxi *
   *******************************/
  
  
  //! default constructor with given problem
  template<class T>
  void VolumetricSource_MaxwellAxi<T>::InitDefaultValues(bool modified_formulation)
  { 
    fsrc = NULL;
    uniform_source = false;
    reference_source = -1;
    
    // omega coefficient in order to have the same solution with modified formulation
    // and mixed H(curl) formulation
    coef_vol = var_source.GetCoefficientVolumeSource();
    if (modified_formulation)
      coef_vol *= var_problem.GetOmega();
    
    number_mode = var_boundary.GetCurrentModeNumber();
    Real_wp kr = max(var_problem.GetOmega()*var_problem.GetXmax(),
		     Real_wp(abs(this->number_mode)));
    
    int N = toInteger(6*kr);
    gauss_teta.ConstructQuadrature(N);
    
    modified_form = false;
    if (modified_formulation && (this->number_mode != 0))
      modified_form = true;

    GaussianSourceField<T, Dimension3>* f_gauss = new GaussianSourceField<T, Dimension3>();
    Vector<T> polar;
    var_problem.GetPolarization(polar);
 
    for (int num = 0; num < var_source.GetNbParameterSource(0); num++)
      {
	const VectString& param = var_source.GetParameterSource(0, num);
	if ((param.GetM() > 1) && (param(0) == "SRC_VOLUME"))
	  {
	    int nb = 2;
	    if (param(1) == "GAUSSIAN")
	      {
		var_source.InitGaussianParameter(*f_gauss, param, nb);
		f_gauss->SetPolarization(polar);
		fsrc = f_gauss;
	      }
            else if (param(1) == "REF_UNIFORM")
              {
                reference_source = to_num<int>(param(2));
                uniform_source = true;
              }
            else if (param(1) == "UNIFORM")
              {
                Real_wp xmin = to_num<Real_wp>(param(nb++));
                Real_wp xmax = to_num<Real_wp>(param(nb++));
                Real_wp ymin = to_num<Real_wp>(param(nb++));
                Real_wp ymax = to_num<Real_wp>(param(nb++));
                T Jx = to_num<T>(param(nb++));
                T Jy = to_num<T>(param(nb++));
                left_corner.PushBack(R2(xmin, ymin));
                right_corner.PushBack(R2(xmax, ymax));
                liste_current.PushBack(TinyVector<T, 2>(Jx, Jy));
              }
	  }
      }
  }
  

  template<class T>
  VolumetricSource_MaxwellAxi<T>::~VolumetricSource_MaxwellAxi()
  {
    if (fsrc != NULL)
      delete fsrc;
  }
  
  
  //! volumetric source
  template<class T>
  bool VolumetricSource_MaxwellAxi<T>
  ::IsNonNullVolumetricSource(const VectR2& s)
  {
    if (uniform_source)
      {
        if (ref_domain == reference_source)
          return true;
        
        return false;
      }

    return true;
  }
   

  template<class T>
  void VolumetricSource_MaxwellAxi<T>::InitElement(int i, const VectR2& s)
  {
    int nb_pts_quad = var_problem.GetNbPointsQuadratureInside(i);
    EvalF.Reallocate(nb_pts_quad);
    ref_domain = var_problem.mesh.Element(i).GetReference();
  }

  template<>
  void VolumetricSource_MaxwellAxi<Real_wp>::EvaluateUniformSource(int i, int j, const R2& x, TinyVector<Real_wp, 3>& f)
  {    
  }


  template<>
  void VolumetricSource_MaxwellAxi<Complex_wp>::EvaluateUniformSource(int i, int j, const R2& x, TinyVector<Complex_wp, 3>& f)
  {
    //DISP(i); DISP(j); DISP(x);
    const Real_wp r(x(0));
    f.Zero();

    TinyVector<R2, 2> interX;
    TinyVector<R2, 3> interY;
    TinyVector<R2, 5> liste_inter_teta;
    int m = this->number_mode;
    //Real_wp err_cos(0), err_sin(0);
    for (int n = 0; n < left_corner.GetM(); n++)
      {
        Real_wp x1 = left_corner(n)(0);
        Real_wp x2 = right_corner(n)(0);
        Real_wp y1 = left_corner(n)(1);
        Real_wp y2 = right_corner(n)(1);

        Real_wp teta1(0), teta2(0);
        if (x1 <= -r)
          teta1 = pi_wp;
        else if (x1 >= r)
          teta1 = Real_wp(0);
        else
          teta1 = acos(x1/r);

        if (teta1 == Real_wp(0))
          continue;

        if (x2 <= -r)
          teta2 = pi_wp;
        else if (x2 >= r)
          teta2 = Real_wp(0);
        else
          teta2 = acos(x2/r);

        if (teta2 == pi_wp)
          continue;

        int nb_seg_x = 2;
        if (teta2 == Real_wp(0))
          {
            // cas d'un seul intervalle
            nb_seg_x = 1;
            interX(0).Init(-teta1, teta1);
          }
        else
          {
            interX(0).Init(teta2, teta1);
            interX(1).Init(-teta1, -teta2);
          }

        Real_wp tetaY1(0), tetaY2(0);
        if (y1 <= -r)
          tetaY1 = -pi_wp/2;
        else if (y1 >= r)
          tetaY1 = pi_wp/2;
        else
          tetaY1 = asin(y1/r);

        if (y1 >= r)
          continue;

        if (y2 <= -r)
          tetaY2 = -pi_wp/2;
        else if (y2 >= r)
          tetaY2 = pi_wp/2;
        else
          tetaY2 = asin(y2/r);

        if (y2 <= -r)
          continue;
        
        int nb_seg_y = 2;
        if (y2 >= r)
          {
            // cas d'un seul intervalle
            if (tetaY1 >= 0)
              {
                nb_seg_y = 1;
                interY(0).Init(tetaY1, pi_wp-tetaY1);
              }
            else
              {
                // on met deux intervalles pour rester entre -pi et pi
                interY(0).Init(tetaY1, pi_wp);
                interY(1).Init(-pi_wp, -tetaY1-pi_wp);
              }
          }
        else if (y1 <= -r)
          {
            if (tetaY2 <= 0)
              {
                nb_seg_y = 1;
                interY(0).Init(-pi_wp-tetaY2, tetaY2);
              }
            else
              {
                nb_seg_y = 2;
                interY(0).Init(pi_wp-tetaY2, pi_wp);
                interY(1).Init(-pi_wp, tetaY2);
              }
          }
        else
          {
            if (tetaY1*tetaY2 < 0)
              {
                // cas avec trois intervalles pour etre dans -pi, pi
                nb_seg_y = 3;
                interY(0).Init(pi_wp-tetaY2, pi_wp);
                interY(1).Init(-pi_wp, -pi_wp-tetaY1);
                interY(2).Init(tetaY1, tetaY2);
              }
            else if (tetaY2 <= 0)
              {
                interY(0).Init(-pi_wp-tetaY2, -pi_wp-tetaY1);
                interY(1).Init(tetaY1, tetaY2);
              }
            else
              {
                interY(0).Init(tetaY1, tetaY2);
                interY(1).Init(pi_wp-tetaY2, pi_wp-tetaY1);
              }
          }

        // intersection entre les intervalles interX et interY
        int nb_inter_teta = 0;
        for (int j = 0; j < nb_seg_x; j++)
          for (int k = 0; k < nb_seg_y; k++)
            {
              Real_wp a = interX(j)(0), b = interX(j)(1);
              Real_wp c = interY(k)(0), d = interY(k)(1);
              if (d <= a)
                continue;

              if (c <= a)
                {
                  if (d <= b)
                    liste_inter_teta(nb_inter_teta++).Init(a, d);
                  else
                    liste_inter_teta(nb_inter_teta++).Init(a, b);
                }
              else
                {
                  if (c >= b)
                    continue;

                  if (d <= b)
                    liste_inter_teta(nb_inter_teta++).Init(c, d);
                  else
                    liste_inter_teta(nb_inter_teta++).Init(c, b);
                }
            }

        //DISP(n); DISP(r); DISP(x1); DISP(x2); DISP(y1); DISP(y2);
        //DISP(interX); DISP(interY); DISP(nb_inter_teta); DISP(liste_inter_teta);

        // calcul de l'integrale de cos(teta) exp(i m theta) et sin(teta) exp(i m theta)
        Complex_wp Icos(0, 0), Isin(0, 0), phase1, phase2;
        for (int k = 0; k < nb_inter_teta; k++)
          {
            teta1 = liste_inter_teta(k)(0);
            teta2 = liste_inter_teta(k)(1);
            if (m == -1)
              {
                Icos += 0.5*(teta2-teta1);
                Isin -= 0.5*Iwp*(teta2-teta1);
              }
            else
              {
                phase1 = exp(Iwp*(m+1)*teta1);
                phase2 = exp(Iwp*(m+1)*teta2);
                Icos += (phase2 - phase1) / (2.0*Iwp*(m+1));
                Isin -= (phase2 - phase1) / (2.0*(m+1));
              }

            if (m == 1)
              {
                Icos += 0.5*(teta2-teta1);
                Isin += 0.5*Iwp*(teta2-teta1);
              }
            else
              {
                phase1 = exp(Iwp*(m-1)*teta1);
                phase2 = exp(Iwp*(m-1)*teta2);
                Icos += (phase2 - phase1) / (2.0*Iwp*(m-1));
                Isin += (phase2 - phase1) / (2.0*(m-1));
              }
          }
        
        Icos /= 2.0*pi_wp;
        Isin /= 2.0*pi_wp;
        
        // autre methode d'evaluation de l'integrale pour verification
        /*int N = 100000;
        Complex_wp Icos_ref(0), Isin_ref(0);
        for (int k = 0; k < N; k++)
          {
            Real_wp teta = 2.0*pi_wp*Real_wp(k)/N;
            Real_wp x = r*cos(teta), y = r*sin(teta);
            if ((x >= x1) && (x <= x2) && (y >= y1) && (y <= y2))
              {
                Icos_ref += cos(teta)*exp(Iwp*m*teta);
                Isin_ref += sin(teta)*exp(Iwp*m*teta);
              }
          }
        
        Icos_ref /= Real_wp(N);
        Isin_ref /= Real_wp(N); 
        
        err_cos = max(err_cos, abs(Icos_ref - Icos));
        err_sin = max(err_sin, abs(Isin_ref - Isin)); 

        Icos = Icos_ref; Isin = Isin_ref; */
        
        // on met a jour la source
        Complex_wp Jx = liste_current(n)(0);
        Complex_wp Jy = liste_current(n)(1);
        f(0) += (Jx*Icos + Jy*Isin);
        f(2) += Iwp*(-Jx*Isin + Jy*Icos);
      }
  }

  
  //! evaluation of gaussian
  template<>
  void VolumetricSource_MaxwellAxi<Real_wp>
  ::EvaluateVolumetricSource(int i, int j, const R2& x, Vector<Real_wp>& f)
  {
    // integration with respect to theta
    Vector<Real_wp> feval(f.GetM());
    feval.Fill(0);
    R3 pt3D; Real_wp teta, poids_cos, poids_sin, cos_teta, sin_teta;
    EvalF(j).Fill(0);
    for (int k = 0; k < gauss_teta.GetNbPointsQuad(); k++)
      {
        teta = 2*pi_wp*gauss_teta.Points(k);
        cos_teta = cos(teta); sin_teta = sin(teta);
        pt3D(0) = x(0)*cos_teta;
        pt3D(1) = x(0)*sin_teta;
        pt3D(2) = x(1);
        fsrc->EvaluateFunction(pt3D, feval);
        poids_cos = gauss_teta.Weights(k)*cos(this->number_mode*teta);
	poids_sin = gauss_teta.Weights(k)*sin(this->number_mode*teta);
	EvalF(j)(0) += (feval(0)*cos_teta + feval(1)*sin_teta)*poids_cos;
        EvalF(j)(1) += (-feval(0)*sin_teta + feval(1)*cos_teta)*poids_sin;
        EvalF(j)(2) += feval(2);
      }
    
    f(0) = EvalF(j)(0);
    f(2) = EvalF(j)(1);
    f(1) = EvalF(j)(2);
    if (this->number_mode > 0)
      Mlt(Real_wp(2), f);
    
    Mlt(coef_vol*x(0), f);
  }
  
  
  //! Evaluation of volumetric source f (term \int f \varphi)
  template<>
  void VolumetricSource_MaxwellAxi<Complex_wp>::
  EvaluateVolumetricSource(int i, int j, const R2& x, Vector<Complex_wp>& f)
  {
    Vector<Complex_wp> feval(f.GetM());
    feval.Zero();

    EvalF(j).Zero();    
    if (uniform_source)
      EvaluateUniformSource(i, j, x, EvalF(j));
    else
      {
        // integration in theta
        R3 pt3D; Real_wp teta, cos_teta, sin_teta; Complex_wp poids;
        for (int k = 0; k < gauss_teta.GetNbPointsQuad(); k++)
          {
            teta = 2*pi_wp*gauss_teta.Points(k);
            cos_teta = cos(teta); sin_teta = sin(teta);
            pt3D(0) = x(0)*cos_teta;
            pt3D(1) = x(0)*sin_teta;
            pt3D(2) = x(1);
            fsrc->EvaluateFunction(pt3D, feval);
            poids = gauss_teta.Weights(k)*exp(Iwp*Real_wp(this->number_mode)*teta);
            feval(0) *= poids; feval(1) *= poids; feval(2) *= poids;
            EvalF(j)(0) += feval(0)*cos_teta + feval(1)*sin_teta;
            EvalF(j)(1) += Iwp*(-feval(0)*sin_teta + feval(1)*cos_teta);
            EvalF(j)(2) += feval(2);
          }        
      }

    
    if (this->number_mode < 0)
      EvalF(j)(1) = -EvalF(j)(1);    
    
    if (this->modified_form)
      {
        f(0) = x(0)*EvalF(j)(0)*abs(this->number_mode);
        f(2) = EvalF(j)(1)*square(this->number_mode) + EvalF(j)(0)*abs(this->number_mode);
        f(1) = x(0)*EvalF(j)(2)*abs(this->number_mode);
      }
    else
      {
        f(0) = EvalF(j)(0);
        f(2) = EvalF(j)(1);
        f(1) = EvalF(j)(2);
      }
    
    Mlt(coef_vol*x(0), f);
  }


  template<class T>
  bool VolumetricSource_MaxwellAxi<T>
  ::IsNonNullGradientSource(const VectR2& s)
  {
    if (this->modified_form)
      return IsNonNullVolumetricSource(s);
    
    return false;
  }


  template<>
  void VolumetricSource_MaxwellAxi<Real_wp>
  ::EvaluateGradientSource(int i, int j, const R2& x, Vector<Real_wp>& f)
  {
  }


  template<>
  void VolumetricSource_MaxwellAxi<Complex_wp>
  ::EvaluateGradientSource(int i, int j, const R2& x, Vector<Complex_wp>& f)
  {
    f(1) = -x(0)*coef_vol*x(0)*EvalF(j)(0)*abs(this->number_mode);
    f(0) = x(0)*coef_vol*x(0)*EvalF(j)(2)*abs(this->number_mode);
  }
  
}

#define MONTJOIE_FILE_DEFINE_SOURCE_MAXWELL_AXI_CXX
#endif
