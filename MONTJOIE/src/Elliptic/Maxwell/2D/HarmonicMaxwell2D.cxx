#ifndef MONTJOIE_FILE_HARMONIC_MAXWELL_2D_CXX

namespace Montjoie
{

  
  /******************************
   * HarmonicMaxwellEquation_2D *
   ******************************/
  
  
  //! computation of mass coefficients
  /*!
    \param[in,out] var given problem
    \param[in] num_elem element number
    \param[in] N number of quadrature points
    \param[in] Points points after transformation F_i
    \param[in] mat jacobian matrices
   */
  template<class TypeEquation>
  void HarmonicMaxwellEquation_2D::
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

    int ref = var.mesh.Element(num_elem).GetReference();
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;
    
    var.Glob_matMass_Bh(num_elem).Reallocate(N);
    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    
    int i1 = num_elem - var.mesh.GetNbElt() + var.GetNbEltPML();
    Complex_wp mu, invMu;
    Matrix2_2sym_Complex_wp epsilon, sigma;
    
    for (int j = 0; j < N; j++)
      {
	    // computation of epsilon and mu on quadrature point if necessary
	epsilon = var.ref_epsilon(ref).GetCoefficient(var, num_elem, j);
	sigma = var.ref_sigma(ref).GetCoefficient(var, num_elem, j);
	mu = var.ref_mu(ref).GetCoefficient(var, num_elem, j);
	
	if (var.InsidePML(num_elem))
	  var.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, j);
	
        invMu = 1.0/mu;
	
	// mass coefficient -> -\omega^2 \omega_k J_i DF_i^{-1} epsilon DF_i^{*-1}
	Matrix2_2 dfj, dfjm1; Real_wp jacob; 
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
	
	Matrix2_2_Complex_wp A_tmp;
	Matrix2_2sym_Complex_wp Bmass; Complex_wp Astiff;
	
	MltTrans(epsilon, dfjm1, A_tmp);
	Real_wp poids_mass = Fb.WeightsND(j);
	if (N == 1)
	  poids_mass = 1.0;
	
	Mlt(dfjm1, A_tmp, Bmass);
	if (var.FirstOrderFormulation())
	  Mlt(-Iwp*var.GetOmega()*jacob*poids_mass, Bmass);	    
	else
	  Mlt(-var.GetSquareOmega()*jacob*poids_mass, Bmass);	    
	
	var.Glob_matMass_Bh(num_elem)(j) = Bmass;
	
	MltTrans(sigma, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Bmass);
        
        if (var.FirstOrderFormulation())
          var.Glob_matMass_Bh(num_elem)(j) += jacob*poids_mass*Bmass;
        else
          Add(-Iwp*var.GetOmega()*jacob*poids_mass, Bmass, var.Glob_matMass_Bh(num_elem)(j));
        
	// stiffness coefficient -> \omega_k/(mu * J_i)
	if (var.FirstOrderFormulation())
	  Astiff = -Iwp*var.GetOmega()*poids_mass*jacob*mu;
	else
	  Astiff = poids_mass/jacob*invMu;
	
	var.Glob_matMass_Dh(num_elem)(j) = Astiff;
      }
  }
  
#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! computation of mass coefficients for first family on quadrilaterals
  /*!
    \param[in,out] var given problem
    \param[in] num_elem element number
    \param[in] N number of quadrature points
    \param[in] Points points after transformation F_i
    \param[in] mat jacobian matrices
   */
  template<class TypeEquation>
  void HarmonicMaxwellEquation_2D::
  ComputeMassMatrixQuad(EllipticProblem<TypeEquation>& var,
			int num_elem, const QuadrangleHcurlFirstFamily& Fb)
  {
    int ref = var.mesh.Element(num_elem).GetReference();
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = (Fb.GetOrder()+1)*(Fb.GetOrder()+1);
    if (!variable)
      N = 1;
    
    if (var.IsVaryingMedia(ref))
      {
        cout << "variable case not working because of Gauss quadrature points" << endl;
        abort();
      }

    var.Glob_matMass_Bh(num_elem).Reallocate(N);

    int i1 = num_elem - var.mesh.GetNbElt() + var.GetNbEltPML();
    Complex_wp mu, invMu;
    Matrix2_2sym_Complex_wp epsilon, sigma;
    
    for (int j = 0; j < N; j++)
      {
	// computation of epsilon and mu on quadrature point if necessary
	epsilon = var.ref_epsilon(ref).GetCoefficient(var, num_elem, j);
	sigma = var.ref_sigma(ref).GetCoefficient(var, num_elem, j);
	mu = var.ref_mu(ref).GetCoefficient(var, num_elem, j);
	
	if (var.InsidePML(num_elem))
	  var.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, j);
	
        invMu = 1.0/mu;
	
	// mass coefficient -> -\omega^2 \omega_k J_i DF_i^{-1} epsilon DF_i^{*-1}
	Real_wp jacob;
	Matrix2_2 dfj, dfjm1;
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
	
	Matrix2_2_Complex_wp A_tmp;
	Matrix2_2sym_Complex_wp Bmass;
	    
	Real_wp poids_mass = Fb.WeightsND(j);
	if (N == 1)
	  poids_mass = 1.0;
	    
	// stiffness coefficient -> \omega_k/(mu * J_i)
	MltTrans(epsilon, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Bmass);
	if (var.FirstOrderFormulation())
	  Mlt(-Iwp*var.GetOmega()*jacob*poids_mass, Bmass);
	else
	  Mlt(-var.GetSquareOmega()*jacob*poids_mass, Bmass);
	
	var.Glob_matMass_Bh(num_elem)(j) = Bmass;
      
        MltTrans(sigma, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Bmass);
        
        if (var.FirstOrderFormulation())
          var.Glob_matMass_Bh(num_elem)(j) += jacob*poids_mass*Bmass;
        else
          Add(-Iwp*var.GetOmega()*jacob*poids_mass, Bmass, var.Glob_matMass_Bh(num_elem)(j));
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
    Complex_wp Astiff;
    for (int j = 0; j < N; j++)
      {
	mu = var.ref_mu(ref).GetCoefficient(var, num_elem, Ni+j);
	  
	// stiffness coefficient -> \omega_k/(mu * J_i)
	// on Gauss points
	if (var.InsidePML(num_elem))
          var.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, Ni+j);
	
        invMu = 1.0/mu;
	
	Real_wp jacob = Det(Mat.GetPointQuadrature(Ni+j));
	Real_wp poids_mass = Fb.WeightsGauss2D(j);
	if (N == 1)
	  poids_mass = 1.0;
	
	if (var.FirstOrderFormulation())
	  Astiff = poids_mass*jacob*mu;
	else
	  Astiff = poids_mass/jacob*invMu;
	
	var.Glob_matMass_Dh(num_elem)(j) = Astiff;
      }
  }
#endif
  

  /**********************
   * HarmonicMaxwell_2D *
   **********************/
  
  
  //! Copies input data from another similar object
  template<class Complexe> 
  void HarmonicMaxwell_2D<Complexe>
  ::CopyInputData(const HarmonicMaxwell_2D<Complexe>& var)
  {
    ref_mu = var.ref_mu;
    ref_invMu = var.ref_invMu;
    ref_epsilon = var.ref_epsilon;
    ref_invEpsilon = var.ref_invEpsilon;
    ref_sigma = var.ref_sigma;
    
    mu0 = var.mu0;
    invMu0 = var.invMu0;
    epsilon0 = var.epsilon0;    
  }
  

  //! initialization of physical indices
  template<class Complexe>
  void HarmonicMaxwell_2D<Complexe>::InitIndices(int n)
  {
    if (n <= 1)
      {
        cout << " ATTENTION VOUS N'AVEZ MIS QU'UN SEUL MILIEU 0!!!!!!!!!" << endl;
        cout << "Number of domains should be greater than 2" << endl; abort();
      }
    
    if (PhysicalConstant::adimensionalization != PhysicalConstant::ADIM_NO)
      {
	// non-physical units, epsilon_0 and mu_0 set to 1
	epsilon0 = 1.0;
	mu0 = 1.0; invMu0 = 1.0;
      }
    else
      {
	// use of physical units
	epsilon0 = PhysicalConstant::epsilon0_permittivity;
	mu0 = PhysicalConstant::mu0_permeability; invMu0 = 1.0/mu0;
      }
    
    ref_epsilon.Reallocate(n);
    ref_sigma.Reallocate(n);
    ref_invEpsilon.Reallocate(n);
    ref_mu.Reallocate(n);
    ref_invMu.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
	ref_epsilon(i).SetDiagonal(epsilon0);
	ref_invEpsilon(i).SetDiagonal(1.0/epsilon0);
	ref_mu(i).SetConstant(mu0);
	ref_invMu(i).SetConstant(invMu0);
        ref_sigma(i).Zero();
      }
  }
  
  
  //! returns the number of physical media
  template<class Complexe>
  int HarmonicMaxwell_2D<Complexe>::GetNbPhysicalIndices() const
  {
    return ref_epsilon.GetM();
  }


  //! reading of relative physical indices \f$ \varepsilon_r \; \mu_r \f$
  template<class Complexe>
  void HarmonicMaxwell_2D<Complexe>::SetIndices(int i, const VectString& parameters)
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
    
    // multiplication of relative physical indices by epsilon_0, mu_0
    ref_epsilon(i).Mlt(epsilon0);
    ref_mu(i).Mlt(mu0);
  }
  

  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Complexe>
  void HarmonicMaxwell_2D<Complexe>
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
      {
        ref_epsilon(i).SetInputData(nb, parameters, parameters(0));
        // multiplication of relative physical indices by epsilon_0, mu_0
        ref_epsilon(i).Mlt(epsilon0);
      }
    else if (name_media == "mu")
      {
        ref_mu(i).SetInputData(nb, parameters, parameters(0));
        ref_mu(i).Mlt(mu0);
      }
    else if (name_media == "sigma")
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }
  
  
  //! returns the name associated with the physical index num
  template<class Complexe>
  string HarmonicMaxwell_2D<Complexe>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("epsilon");
      case 1: return string("mu");
      case 2: return string("sigma");
      }
    
    return string();
  }  
    
  
  //! fills the array with varying indices present in the problem
  template<class Complexe>
  void HarmonicMaxwell_2D<Complexe>::
  GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complexe>* >& rho_complex,
		    IVect& num_ref, IVect& num_index, IVect& num_component,
		    Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    int nb = 0;
    for (int i = 0; i < ref_epsilon.GetM(); i++)
      {
        nb += ref_epsilon(i).GetNbVaryingMedia();
        nb += ref_mu(i).GetNbVaryingMedia();
        nb += ref_sigma(i).GetNbVaryingMedia();
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
    for (int i = 0; i < ref_epsilon.GetM(); i++)
      {
        int nb0 = nb;
        ref_epsilon(i).GetVaryingMedia(nb, rho_complex, num_component);
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
      }
  }
  
  
  //! return the velocity of waves in media of reference ref
  template<class Complexe>
  Real_wp HarmonicMaxwell_2D<Complexe>::GetVelocityOfMedia(int ref) const
  {
    return 1.0/sqrt(abs(this->ref_epsilon(ref)(0,0) * this->ref_mu(ref).GetConstant()));
  }
  

  //! return the velocity of waves in media of reference ref
  template<class Complexe>
  Real_wp HarmonicMaxwell_2D<Complexe>::GetVelocityOfInfinity() const
  {
    return 1.0/sqrt(abs(this->epsilon0 * this->mu0));
  }


  //! computation of invEpsilon and invMu from values of Epsilon and Mu
  template<class Complexe>
  void HarmonicMaxwell_2D<Complexe>::FinalizeComputationVaryingIndices()
  {
    // computation of inverses    
    for (int ref = 1; ref < ref_epsilon.GetM(); ref++)
      {
        ref_invEpsilon(ref) = ref_epsilon(ref);
        ref_invEpsilon(ref).GetInverse();
        ref_invMu(ref) = ref_mu(ref);
        ref_invMu(ref).GetInverse();        
      }
  }

    
  //! returns true if epsilon or mu is variable inside element i
  template<class Complexe>
  bool HarmonicMaxwell_2D<Complexe>::IsVaryingMedia(int i) const
  {
    if (ref_epsilon(i).IsVarying() || ref_mu(i).IsVarying() || ref_sigma(i).IsVarying())
      return true;
    
    return false;
  }
  

  //! updates epsilon and mu inside PML
  template<>
  void HarmonicMaxwell_2D<Complex_wp>
  ::ModifyPhysicalCoefPML(Matrix2_2sym_Complex_wp& eps, Complex_wp& mu, Matrix2_2sym_Complex_wp& sigma,
			  int i1, int i) const
  {
    Complex_wp one(1, 0);
    Complex_wp dx = 1.0/var_boundary.GetTauPML(i1, i, 0);
    Complex_wp dy = 1.0/var_boundary.GetTauPML(i1, i, 1);
    
    mu *= dx*dy;
    
    Complex_wp ratio = dy/dx;
    eps(0, 0) *= ratio;
    sigma(0, 0) *= ratio;
    
    ratio = one/ratio;
    eps(1, 1) *= ratio;
    sigma(1, 1) *= ratio;
  }


  template<>
  void HarmonicMaxwell_2D<Real_wp>
  ::ModifyPhysicalCoefPML(Matrix2_2sym& eps, Real_wp& mu, Matrix2_2sym& sigma,
			  int i1, int i) const
  {
  }

    
  //! allocation of arrays before computation of mass matrix
  template<class Complexe>
  void HarmonicMaxwell_2D<Complexe>::AllocateMassMatrices()
  {
    this->Glob_matMass_Bh.Reallocate(var_problem.mesh.GetNbElt());
    this->Glob_matMass_Dh.Reallocate(var_problem.mesh.GetNbElt());
    
    this->coefficient_impedance_absorbing.Reallocate(this->GetNbPhysicalIndices());
    this->coefficient_impedance_absorbing.Fill(1.0);
    for (int i = 0; i < this->coefficient_impedance_absorbing.GetM(); i++)
      {
	Complexe epsilon = this->ref_epsilon(i)(0,0), mu = this->ref_mu(i);
	if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	  this->coefficient_impedance_absorbing(i) = sqrt(mu/epsilon);
	else
	  this->coefficient_impedance_absorbing(i) = sqrt(epsilon/mu);
      }
    
    Mlt(var_boundary.GetImpedanceCoefficientABC(), this->coefficient_impedance_absorbing);
  }

  
  /****************
   * ImpedanceABC *
   ****************/
  

  //! impedance
  void ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_2D>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                         const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                         const SetPoints<Dimension2>& Pts, const SetMatrices<Dimension2>& Mat)
  {
    if (k == 0)
      {
        this->vec_normale.Reallocate(Mat.GetNbPointsQuadratureBoundary());
        for (int kp = 0; kp < Mat.GetNbPointsQuadratureBoundary(); kp++)
          this->vec_normale(kp) = Mat.GetNormaleQuadratureBoundary(kp);
      }

    this->read_param_condition = false;
    this->coef_grad = 0;
    if (var_maxwell.coefficient_impedance_absorbing.GetM() != 0)
      {
        Complex_wp m_iomega = -Iwp*var_problem.GetOmega();
	if (var_problem.FirstOrderFormulation())
	  SetComplexOne(m_iomega);
	
        this->coef_phi = m_iomega*var_maxwell.coefficient_impedance_absorbing(ref)*nat_mat.GetCoefDamping();
      }
    else
      {
        abort();
      }
  }
  
  
  /*************************
   * HarmonicMaxwell_2D_Eq *
   *************************/


  template<class TypeEquation>
  HarmonicMaxwell_2D_Eq<TypeEquation>::HarmonicMaxwell_2D_Eq()
    : HarmonicMaxwell_2D<typename TypeEquation::Complexe>(this->GetLeafClass())
  {
  }
  
  
  template<class TypeEquation>
  void HarmonicMaxwell_2D_Eq<TypeEquation>::InitIndices(int n)
  {
    HarmonicMaxwell_2D<Complexe>::InitIndices(n);
  }
  
  
  template<class TypeEquation>
  int HarmonicMaxwell_2D_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return HarmonicMaxwell_2D<Complexe>::GetNbPhysicalIndices();
  }

  
  template<class TypeEquation>
  void HarmonicMaxwell_2D_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    HarmonicMaxwell_2D<Complexe>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  void HarmonicMaxwell_2D_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    HarmonicMaxwell_2D<Complexe>::SetPhysicalIndex(name_media, i, parameters);
  }
    
  
  template<class TypeEquation>
  string HarmonicMaxwell_2D_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return HarmonicMaxwell_2D<Complexe>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  Real_wp HarmonicMaxwell_2D_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return HarmonicMaxwell_2D<Complexe>::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  Real_wp HarmonicMaxwell_2D_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return HarmonicMaxwell_2D<Complexe>::GetVelocityOfInfinity();
  }
  
  
  template<class TypeEquation>
  void HarmonicMaxwell_2D_Eq<TypeEquation>::FinalizeComputationVaryingIndices()
  {
    HarmonicMaxwell_2D<Complexe>::FinalizeComputationVaryingIndices();
  }
  
  
  template<class TypeEquation>
  bool HarmonicMaxwell_2D_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return HarmonicMaxwell_2D<Complexe>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool HarmonicMaxwell_2D_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return HarmonicMaxwell_2D<Complexe>::IsVaryingMedia(i);
  }
    
  
  template<class TypeEquation>
  void HarmonicMaxwell_2D_Eq<TypeEquation>::AllocateMassMatrices()
  {
    HarmonicMaxwell_2D<Complexe>::AllocateMassMatrices();
  }


  template<class TypeEquation>
  void HarmonicMaxwell_2D_Cplx<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
		      IVect& num_ref, IVect& num_index, IVect& num_component,
		      Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    HarmonicMaxwell_2D<Complexe>::GetVaryingIndices(rho_complex, num_ref, num_index, num_component,
						    compute_grad, compute_hess);
  }


  IncidentWaveProjector<Complex_wp, Dimension2>*
  EllipticProblem<HarmonicMaxwellEquation_2D>
  ::GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) const
  {
    return new IncidentWaveProjector_Maxwell2D(*this, u_inc);
  }


  void EllipticProblem<HarmonicMaxwellEquation_2D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicMaxwellEquation_2D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat,
				      *this, this->GetReferenceElementHcurl(i));
  }


  void FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_2D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Not possible" << endl;
    abort();
  }


  void FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_2D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddHcurl2D(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		   X, Complex_wp(1, 0), Y, false);
  }
  
} // end namespace

#define MONTJOIE_FILE_HARMONIC_MAXWELL_2D_CXX
#endif
