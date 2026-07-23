#ifndef MONTJOIE_FILE_HELIO_AXISYM_HELMHOLTZ_CXX

namespace Montjoie
{
  
  //! Default constructor.
  template<class TypeEquation>
  HELIO_VarHelmholtzAxisym<TypeEquation>::HELIO_VarHelmholtzAxisym()
  {
    alpha0 = 1.0;
    beta0  = 1.0;
    
    R_sun  = 696.8e6;
  }
  
  
  //! Initializes parameters from input file.
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    
    // Specific entries for HELIO_HelmAxi equation
    if (!description_field.compare("CalculEnveloppe"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of HELIO" << endl;
	    cout << "CalculEnveloppe needs 1 parameters, for instance :" << endl;
	    cout << "CalculEnveloppe = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        if (parameters(0) == "YES") 
          this->calcul_enveloppe = true;
        else
          this->calcul_enveloppe = false;
	
        if (this->calcul_enveloppe) 
          {
            cerr << "Enveloppe computation not implemented (yet?). Abort." << endl;
            abort();
          }	
      }
    else 
      VarHelmholtzAxi_Eq<TypeEquation>::SetInputData(description_field, parameters);
    
  }
  

  //! Computes FEM matrix terms of element k for the helioseismology axisymmetric Helmholtz equation
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>
  ::ComputeElementaryMatrix(int   ID_elt, IVect& num_dof, VirtualMatrix<Complexe>& mat,
			    const GlobalGenericMatrix<Complexe>& nat_mat)
  {
    
    const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElement(ID_elt);
    
    // Initializing matrix
    int nb_dof_elt = Fb.GetNbDof();
    mat.Reallocate(nb_dof_elt, nb_dof_elt);
    mat.Fill(0);
    
    // Retrieve degrees of freedom of element iquad
    num_dof = this->GetDofNumberOnElement(ID_elt); 
    int N   = Fb.GetNbPointsQuadratureInside();
    
    // The following matrices will be computed:
    // \int CoefMass  \phi_i\phi_j
    // \int CoefGrad  \nabla\phi_i\phi_j
    // \int CoefGradT \phi_i\nabla\phi_j
    // \int CoefStiff \nabla\phi_i\nabla\phi_j
    Vector<Complexe>                            CoefMass (N);
    Vector<TinyVector<Complexe, 2> >             CoefGrad (N);
    Vector<TinyVector<Complexe, 2> >             CoefGradT(N);
    Vector<TinyMatrix<Complexe, General, 2, 2> > CoefStiff(N);
    
    // Some shortcuts and initializations
    Complexe omega;     omega = this->GetOmega();
    Complexe m_iomega;  this->GetMiomega(m_iomega);
    Complexe m_omega2;  this->GetMomega2(m_omega2);
    
    int  ref    = this->mesh.Element(ID_elt).GetReference();
    bool affine = this->mesh.IsElementAffine(ID_elt);
    bool fr3    = this->formulation_r3 && (this->number_mode != 0);
    bool flow   = this->include_flow_term;

    Real_wp m = Real_wp(this->number_mode);
    Real_wp jacob;
    
    TinyMatrix<Complexe, Symmetric, 2, 2> A;
    TinyMatrix<Real_wp, General, 2, 2> dFm1;
    
    // Loop on all volumetric quadrature points of the element (\xi_j)
    for (int ID_pt = 0; ID_pt < N; ID_pt++)
      {
	// We get the transformation jacobian at point \xi_j
	if (affine)
	  {
	    jacob = this->Glob_jacobian(ID_elt)(0)*Fb.WeightsND(ID_pt);
	    dFm1  = this->Glob_DFjm1   (ID_elt)(0);                
	    Mlt(1.0/this->Glob_jacobian(ID_elt)(0),dFm1);
	  }
	else
	  {
	    jacob = this->Glob_jacobian(ID_elt)(ID_pt);
	    dFm1  = this->Glob_DFjm1   (ID_elt)(ID_pt);
	    Mlt(Fb.WeightsND(ID_pt)/jacob, dFm1);
	  }
	
	// Then the values of physical indices at point \xi_j
	Complexe r      = this->Glob_rtilde(ID_elt)(ID_pt);
	Complexe rho    = this->ref_rho  (ref).GetCoefficient(*this, ID_elt, ID_pt);
	Complexe mu     = this->ref_mu   (ref).GetCoefficient(*this, ID_elt, ID_pt);	    
	Complexe sigma  = this->ref_sigma(ref).GetCoefficient(*this, ID_elt, ID_pt);
	Complexe alpha;
	Complexe beta ;
	TinyVector<Complexe, 2> Galpha,dFGalpha;
	TinyVector<Complexe, 2> Gbeta ,dFGbeta ;
	TinyVector<Complexe, 2> dFGr;
	
	this->ref_alpha(ref).GetCoefGradient(*this, ID_elt, ID_pt, alpha, Galpha);
	this->ref_beta (ref).GetCoefGradient(*this, ID_elt, ID_pt, beta , Gbeta );
	
	Mlt(dFm1, Galpha, dFGalpha);
	Mlt(dFm1, Gbeta , dFGbeta );
	dFGr.Init(dFm1(0,0), dFm1(1,0));
	
	TinyVector<Complexe,3>   v0(0.0);
	TinyVector<Complexe,2> dFv0(0.0);
	if (flow)
	  {
	    v0 = this->ref_v0(ref).GetCoefficient(*this, ID_elt, ID_pt);
	    TinyVector<Complexe, 2> v0temp(v0(0),v0(2));
	    Mlt(dFm1, v0temp, dFv0);
	  }
	
	Complexe CM = nat_mat.GetCoefMass();
	Complexe CD = nat_mat.GetCoefDamping();
	Complexe CS = nat_mat.GetCoefStiffness();
	
	// Matrices computation
	if (!fr3) // "Normal" variational formulation
	  {
	    
	    // Mass matrix
            
	    CoefMass(ID_pt)  = r*(m_omega2*rho*CM + m_iomega*sigma*CD) -2.0*m*omega*v0(1)*CM;
	    //CoefMass(ID_pt)  = r*(m_omega2*rho*CM + m_iomega*sigma*CD + sigma*sigma*CD/4.e0) -2.0*m*omega*v0(1)*CM;
	    CoefMass(ID_pt) += ((abs(r)>epsilon_machine)? m*m/r*alpha*beta*mu*CS : Complexe(0.0));
	    CoefMass(ID_pt) += r*mu*DotProd(dFGalpha,dFGbeta)*CS;
	    CoefMass(ID_pt) *= jacob;
	    
	    // Gradient matrices
	    CoefGrad(ID_pt)  = mu*alpha*dFGbeta*CS;
	    CoefGrad(ID_pt) += m_iomega*dFv0*CD;
	    CoefGrad(ID_pt) *= r*jacob;
	    
	    CoefGradT(ID_pt)  = mu*beta*dFGalpha*CS;
	    CoefGradT(ID_pt) -= m_iomega*dFv0*CD;
	    CoefGradT(ID_pt) *= r*jacob;
	    
	    // Stiffness matrix
	    MltTrans(dFm1,dFm1,A);
	    CoefStiff(ID_pt)  = A;
	    CoefStiff(ID_pt) *= jacob*r*alpha*beta*mu*CS;
	    
	  }
	else // formulation_r3
	  {
	    
	    // Mass matrix
	    CoefMass(ID_pt)  = r*  (r* (r*(m_omega2*rho*CM + m_iomega*sigma*CD)-2.0*m*omega*v0(1)*CD) + m*m*alpha*beta*mu*CS);
	    CoefMass(ID_pt) += r*r*mu*CS*(r*DotProd(dFGalpha,dFGbeta) + alpha*dFGbeta(0) + dFGalpha(0)*beta);
	    CoefMass(ID_pt) *= jacob;
	    
	    // Gradient matrices
	    CoefGrad(ID_pt)  = m_iomega*dFv0*CD;
	    CoefGrad(ID_pt) += r*mu*alpha*(dFGbeta +beta*dFGr)*CS; 
	    CoefGrad(ID_pt) *= r*r*jacob; 
	    
	    CoefGradT(ID_pt)  = -m_iomega*dFv0*CD;
	    CoefGradT(ID_pt) += r*mu*beta*(dFGalpha +alpha*dFGr)*CS; 
	    CoefGradT(ID_pt) *= r*r*jacob; 
	    
	    // Stiffness matrix
	    MltTrans(dFm1,dFm1,A);
	    CoefStiff(ID_pt)  = A;
	    CoefStiff(ID_pt) *= jacob*r*r*r*alpha*beta*mu*CS;
	    
	  }
      }
    
    // Element matrices are added to the global matrix
    Fb.AddVariableMassMatrix             (0,0,CoefMass ,mat);    
    Fb.AddVariableGradientMatrix         (0,0,CoefGrad ,mat);
    Fb.AddVariableTransposeGradientMatrix(0,0,CoefGradT,mat);
    Fb.AddVariableStiffnessMatrix        (0,0,CoefStiff,mat);
    
    // Modification of signs if needed
    
    mesh_num.number_map.ModifyLocalRowMatrix   (this->mesh_num, mat, ID_elt);
    mesh_num.number_map.ModifyLocalColumnMatrix(this->mesh_num, mat, ID_elt);    
    
  }

 
  //! Allocation of arrays containing physical properties.
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>::InitIndices(int n)
  {
    
    VarHelmholtzAxi_Eq<TypeEquation>::InitIndices(n);
    
    ref_alpha.Reallocate(n+1);
    ref_beta .Reallocate(n+1);
    for (int i = 0; i <= n; i++)
      {
	ref_alpha(i).SetIdentity(); 
	ref_alpha(i).Mlt(alpha0);
	ref_beta (i).SetIdentity();
	ref_beta (i).Mlt(beta0 );
      }    
  }


  //! Reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>::
  SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    
    int nb = 1;
    if (i >= this->ref_rho.GetM())
      {
	cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
	cout << "Current reference i : " << i << " < " << this->ref_rho.GetM() << endl;
	abort();
      }
    
    if      (name_media == "rho")
      this->ref_rho(i)  .SetInputData(nb, parameters, parameters(0));
    else if (name_media == "mu")
      this->ref_mu(i)   .SetInputData(nb, parameters, parameters(0));
    else if (name_media == "alpha")
      this->ref_alpha(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "beta")
      this->ref_beta(i) .SetInputData(nb, parameters, parameters(0));
    else if (name_media == "sigma")
      this->ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "M")
      {
	if (this->include_flow_term)
	  this->ref_v0(i).SetInputData(nb, parameters, parameters(0));
      }
    else
      {
	cout << "Unknown media : " << name_media << endl;
	abort();
      }
    
  }
  

  //! returns the name associated with the physical index num
  template<class TypeEquation>
  string HELIO_VarHelmholzAxisym<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("rho");
      case 1: return string("mu");
      case 2: return string("alpha");
      case 2: return string("beta");
      case 2: return string("sigma");
      case 2: return string("M");
      }

    return string();
  }  
    

  //! Sets physical properties of media with reference i.
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {

    if (parameters.GetM() <= 2)
      {
	cout << "Not enough parameters in MateriauDielec = " << endl;
	cout << "Current parameters are : " << endl << parameters << endl;
	abort();
      }
    
    int nb = 1;
    this->ref_rho  (i).SetInputData(nb, parameters, parameters(0));
    this->ref_mu   (i).SetInputData(nb, parameters, parameters(0));
    this->ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    this->ref_alpha(i).SetInputData(nb, parameters, parameters(0));
    this->ref_beta (i).SetInputData(nb, parameters, parameters(0));
    if (this->include_flow_term)
       this->ref_v0(i).SetInputData(nb, parameters, parameters(0));

  }

  
  //! Fills arrays with the varying indices present in the current problem. 
  /*!
    \param[out] rho_complex will contain pointers to physical indices
    \param[out] num_ref corresponding references
    \param[out] num_index corresponding number defining the index nature
    \param[out] num_component 0 for scalars, 0,1 or 2 for vectorial indices
  */
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>::
  GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
		    Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
                    IVect& num_ref, IVect& num_index, IVect& num_component, 
                    Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {

    // We first compute the total number of varying media
    int nb = 0;
    for (int i = 0; i < this->ref_rho.GetM(); i++)
      {
        nb += this->ref_rho  (i).GetNbVaryingMedia();
        nb += this->ref_mu   (i).GetNbVaryingMedia();
        nb += this->ref_sigma(i).GetNbVaryingMedia();
        nb += this->ref_alpha(i).GetNbVaryingMedia();
        nb += this->ref_beta (i).GetNbVaryingMedia();
	if (this->include_flow_term)
	   nb += this->ref_v0(i).GetNbVaryingMedia();
      }

    rho_complex  .Reallocate(nb);
    num_ref      .Reallocate(nb);
    num_index    .Reallocate(nb);
    num_component.Reallocate(nb);
    compute_grad .Reallocate(nb);
    compute_hess .Reallocate(nb);

    compute_grad.Fill(false);
    compute_hess.Fill(false);
    nb = 0;

    // Listing all media properties
    for (int i = 0; i < this->ref_rho.GetM(); i++)
      {
        int nb0 = nb;
        this->ref_rho(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        this->ref_mu(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 1;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        this->ref_sigma(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 2;
            num_ref(j) = i;
          }

        nb0 = nb;
        this->ref_alpha(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 3;
            num_ref(j) = i;
            compute_grad(j) = true;
          }

        nb0 = nb;
        this->ref_beta(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 4;
            num_ref(j) = i;
            compute_grad(j) = true;
          }

	if (this->include_flow_term)
	  {
	    nb0 = nb;
	    this->ref_v0(i).GetVaryingMedia(nb, rho_complex, num_component);
	    for (int j = nb0; j < nb; j++)
	      {
		num_index(j) = 5;
		num_ref(j) = i;
	      }
	  }

      }

  }
  

  //! Finds the maximum values of physical constants
  //! among all references precised in RefUsed.
  template<class TypeEquation>
  void HELIO_VarHelmholtzAxisym<TypeEquation>::
  SetPhysicalIndexAtInfinity(const Vector<bool>& RefUsed)
  {
    
    // Taking the maximum among elements of the current processor
    Vector<Real_wp> coef(5);
    coef.Fill(-1);
    for (int ref = 0; ref < RefUsed.GetM(); ref++)
      if (RefUsed(ref))
	{
	  coef(0) = max(coef(0),abs(this->ref_rho  (ref).GetConstant()));
	  coef(1) = max(coef(1),abs(this->ref_mu   (ref).GetConstant()));
	  coef(2) = max(coef(2),abs(this->ref_sigma(ref).GetConstant()));
	  coef(3) = max(coef(3),abs(this->ref_alpha(ref).GetConstant()));
	  coef(4) = max(coef(4),abs(this->ref_beta (ref).GetConstant()));
	}
    
    // Then taking the maximum among processors
#ifdef SELDON_WITH_MPI
    Vector<Real_wp> max_coef0(5);
    Vector<int64_t> xtmp;
    MpiAllreduce(this->comm_group_mode, coef, xtmp, max_coef0, 5, MPI::MAX);    
    
    coef = max_coef0;
#endif
    
    this->rho0   = coef(0); 
    this->mu0    = coef(1); 
    this->sigma0 = coef(2);
    this->alpha0 = coef(3);
    this->beta0  = coef(4);
    
    if (coef(0) == -1.0) this->rho0   = 1.0;
    if (coef(1) == -1.0) this->mu0    = 1.0;
    if (coef(2) == -1.0) this->sigma0 = 1.0;
    if (coef(3) == -1.0) this->alpha0 = 1.0;
    if (coef(4) == -1.0) this->beta0  = 1.0;
    
    Complexe m_iomega; this->GetMiomega(m_iomega); 
    Complexe rho0_tilde = this->rho0 + this->sigma0/m_iomega;
    this->kwave_envelope = this->omega*sqrt(rho0_tilde/(this->alpha0*this->beta0*this->mu0));
    
  }
  

  //! evaluation of impedance coefficient
  void ImpedanceABC<Complex_wp, HELIO_HelmholtzAxisymEquation>::
  EvaluateImpedancePhi(int i, int ID_elt, int ID_edge, int ID_pt_loc, int k,
                       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                       const SetPoints<Dimension2>& Pts,
                       const SetMatrices<Dimension2>& Mat)
  {
    if (k == 0)
      {
	int nb_points = Mat.GetNbPointsQuadratureBoundary();
	this->stored_coef_phi.Reallocate(nb_points);
      }
    
    Real_wp  r    = Pts.GetPointQuadratureBoundary(k)(0);
    Complex_wp mu   = var_helm.ref_mu  (ref).GetConstant();
    Complex_wp beta = var_helm.ref_beta(ref).GetConstant();
    
    Complex_wp alpha, dnalpha;
    TinyVector<Complex_wp, 2> Galpha;
    var_helm.ref_alpha(ref).GetCoefGradient(var_helm, ID_elt, ID_pt_loc, alpha, Galpha);
    Galpha /= Complex_wp(var_helm.R_sun, 0);
    
    TinyVector<Real_wp, 2> RGa;
    RGa(0) = realpart(Galpha(0));
    RGa(1) = realpart(Galpha(1));
    dnalpha = DotProd(Galpha, var_helm.Glob_normale(ID_edge)(k));
    
    Complex_wp feval = -r*mu*dnalpha*beta -Iwp*var_helm.GetWaveNumberInfinity()*r*mu
      *alpha*beta*nat_mat.GetCoefDamping();
    //feval = -Iwp*r*var.kwave_envelope*mu*nat_mat.GetCoefDamping();
    
    if ((var_helm.UseFormulationR3()) && (var_helm.GetCurrentModeNumber() != 0))
      feval *= r*r;
    
    this->stored_coef_phi(k) = feval;    
  }
  
} 
#define MONTJOIE_FILE_HELIO_AXISYM_HELMHOLTZ_CXX
#endif


