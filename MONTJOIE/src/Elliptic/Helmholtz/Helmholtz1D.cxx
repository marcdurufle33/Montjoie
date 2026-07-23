#ifndef MONTJOIE_FILE_HELMHOLTZ_1D_CXX

#include "Helmholtz1D_Inline.cxx"

namespace Montjoie
{

  /*****************
   * VarHelmholtz_1D *
   *****************/


  //! default constructor
  template<class Complexe>
  VarHelmholtz_1D<Complexe>::VarHelmholtz_1D()
  {
    SetComplexOne(rho0);
    SetComplexOne(mu0);
    Complexe zero;
    SetComplexZero(zero);
    SetComplexZero(sigma0);

    kwave = this->omega;    
    this->InitIndices(PhysicalConstant::nb_max_indices);
    
    helio_scaling_model = false;
    
    beta_impedance_left = 0;
    beta_impedance_right = 0;   

    SetComplexZero(left_dirichlet_condition);
    SetComplexZero(right_dirichlet_condition);
  }

  
  //! reading the parameters contained in data file
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarProblem_1D::SetInputData(description_field, parameters);
    
    if (!description_field.compare("MateriauDielec"))
      {
        if (parameters.GetM() <= 3)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "MateriauDielec needs more parameters, for instance : " << endl;
            cout << "MateriauDielec = ref rho mu sigma" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
	
        int ref = to_num<int>(parameters(0));
        SetIndices(ref, parameters);
      }
    else if (!description_field.compare("PhysicalMedia"))
      {
        if (parameters.GetM() <= 2)
          {
            cout << "In SetInputData of VarHelmholtz_1D" << endl;
            cout << "PhysiacalMedia needs more parameters, for instance : " << endl;
            cout << "PhysicalMedia = rho ref valeur" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }

        SetPhysicalIndex(parameters(0), to_num<int>(parameters(1)), parameters);
      }
    else if (!description_field.compare("ConventionPhysicalIndex"))
      {
	if (parameters(0) == "Helio")
	  {
	    if (parameters.GetM() <= 2)
	      {
		cout << "In SetInputData of VarHelmholtz_1D" << endl;
		cout << "ConventionPhysicalIndex needs more parameters, for instance : " << endl;
		cout << "ConventionPhysicalIndex = rho ref valeur" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }
	    
	    helio_scaling_model = true;
	    int ref = to_num<int>(parameters(1));
	    log_rho_read(ref) = false;
            if (parameters(2) == "LOG")
	      {
		// in this case, we assume that the logarithm of rho is given
		// instead of rho'/rho such that rho and rho'/rho are directly available
		log_rho_read(ref) = true;
		ref_rho(ref).EnableComputationGradient();
	      }
	    else
	      {
		int nb = 2;
		ref_true_rho(ref).SetInputData(parameters, nb);
	      }
	  }
	else
	  helio_scaling_model = false;
      }
    else if (description_field == "TypeSource")
      {
	if (parameters(0) == "SRC_BOUNDARY")
	  {
	    left_dirichlet_condition = to_num<Complexe>(parameters(1));
	    right_dirichlet_condition = to_num<Complexe>(parameters(2));
	  }
      }
  }
  
  
  //! Evaluates indexes rho mu and sigma
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>
  ::GetPhysicalCoefficient(Complexe& rho, Complexe& mu, Complexe& sigma,
			   Complexe& v0, int i, int j, int ref)
  {
    Complexe drho;
    SetComplexZero(drho);
    if (log_rho_read(ref))
      ref_rho(ref).GetCoefGradient(*this, i, j, rho, drho);
    else
      rho = ref_rho(ref).GetCoefficient(*this, i, j);
    
    mu = ref_mu(ref).GetCoefficient(*this, i, j);
    sigma = ref_sigma(ref).GetCoefficient(*this, i, j);
    SetComplexZero(v0);

    if (helio_scaling_model)
      {
	Complexe ratio_rho = rho;

	// derivative of log(rho) = rho'/rho
	if (log_rho_read(ref))
	  ratio_rho = drho;
	
	Complexe c = mu;
	Complexe gamma = sigma;
	rho = 1.0/square(c);
	SetComplexOne(mu);
	sigma = 2.0*gamma*rho;
	v0 = ratio_rho;
      }
  }
  
  
  //! computes physical coefficients
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::ComputePhysicalCoefficients()
  {
    // first step : we update ElementRho
    int nb_elt = this->mesh.GetNbElt();
    this->ElementRho.Reallocate(nb_elt);
    this->ElementRho.Zero();
    IVect NbElement_ByRef(ref_rho.GetM()); NbElement_ByRef.Zero();
    for (int i = 0; i < nb_elt; i++)
      {
	int ref = this->mesh.Element(i).GetReference();
	this->ElementRho(i) = NbElement_ByRef(ref);
	NbElement_ByRef(ref)++;
      }
    
    // then allocating needed arrays in ref_rho, etc
    Vector< VectReal_wp > points(ref_rho.GetM());
    Vector< Vector<Complexe> > indices_rho(ref_rho.GetM());
    Vector< Vector<Complexe> > indices_mu(ref_rho.GetM());
    Vector< Vector<Complexe> > indices_sigma(ref_rho.GetM());
    for (int ref = 0; ref < ref_rho.GetM(); ref++)
      {
        ref_rho(ref).InitArray(ref, *this);
        ref_mu(ref).InitArray(ref, *this);
        ref_sigma(ref).InitArray(ref, *this);
        ref_true_rho(ref).InitArray(ref, *this);
	
        if (this->print_level >= 4)
	  {
	    points(ref).Reallocate(ref_rho(ref).GetNbPoints());
	    indices_rho(ref).Reallocate(ref_rho(ref).GetNbPoints());
	    indices_mu(ref).Reallocate(ref_rho(ref).GetNbPoints());
	    indices_sigma(ref).Reallocate(ref_rho(ref).GetNbPoints());
	  }
      }
    
    VectReal_wp s(2);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(i);
	// We get the two extremities of the edge
	this->mesh.GetVerticesElement(i, s);
        
	// number of integration points
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
	int ref = this->mesh.Element(i).GetReference();
        int offset = ref_rho(ref).GetOffsetPoint(*this, i);
	
	for (int j = 0; j <= nb_points_quadrature+1; j++)
	  {
            Real_wp x(0);
            if (j == 0)
              x = s(0);
            else if (j == nb_points_quadrature+1)
              x = s(1);
            else
              {
                // point inside the element
                x = Fb.Points(j-1)*s(1) + (1.0-Fb.Points(j-1))*s(0);
              }
                        
            ref_rho(ref).ComputeCoefficient(*this, i, j, x);
            ref_mu(ref).ComputeCoefficient(*this, i, j, x);
            ref_sigma(ref).ComputeCoefficient(*this, i, j, x);
            ref_true_rho(ref).ComputeCoefficient(*this, i, j, x);
           
	    if (this->print_level >= 4)
	      {	
                points(ref)(offset+j) = x;
		indices_rho(ref)(offset+j) = ref_rho(ref).GetCoefficient(*this, i, j);
		indices_mu(ref)(offset+j) = ref_mu(ref).GetCoefficient(*this, i, j);
		indices_sigma(ref)(offset+j) = ref_sigma(ref).GetCoefficient(*this, i, j);
	      }
	  }
      }
    for (int ref = 0; ref < ref_rho.GetM(); ref++)
      if ((points(ref).GetM() > 0) && (this->print_level >= 4))
        {
          points(ref).WriteText("nodal_points_ref" + to_str(ref)+".dat");
          indices_rho(ref).Write("indices_rho_ref" + to_str(ref)+".dat");
          indices_mu(ref).Write("indices_mu_ref" + to_str(ref)+".dat");
          indices_sigma(ref).Write("indices_sigma_ref" + to_str(ref)+".dat");
        }
  }
  
  
  //! computation of geometrical coefficients needed to compute finite element matrix
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::ComputeMassMatrix()
  {    
    int ref; Complexe rho, sigma, mu, invMu, sigma_v, v0, coef_mass, coef_damping;
    VectReal_wp s(2);
    mat_mass.Reallocate(this->mesh.GetNbElt()); 
    mat_sigma.Reallocate(this->mesh.GetNbElt()); 
    mat_massV.Reallocate(this->mesh.GetNbElt());
    mat_flow.Reallocate(this->mesh.GetNbElt());
    if (mixed_formulation)
      {
	mat_sigmaV.Reallocate(this->mesh.GetNbElt());
	if (this->type_coordinate != VarProblem_1D::CARTESIAN)
	  mat_mass_radius.Reallocate(this->mesh.GetNbElt());
	
	if (calcul_enveloppe)
	  {
	    cout << "case not implemented" << endl;
	    abort();
	  }
      }
    else
      mat_sigmaV.Reallocate(this->mesh.GetNbElt());
    
    // position of degrees of freedom
    CoorDofs.Reallocate(this->GetNbDof());
    
    // m_iomega = -i omega for complex numbers, 1 for real numbers
    Complexe m_iomega, m_omega2, zeta, ikwave, zero;
    SetComplexZero(zero);
    this->GetMiomega(m_iomega);
    this->GetIkwave(ikwave);
    m_omega2 = m_iomega*m_iomega;
    
    if ((this->type_coordinate == VarProblem_1D::SPHERICAL) ||
	(this->type_coordinate == VarProblem_1D::POLAR))
      mat_mass_mu.Reallocate(this->mesh.GetNbElt());
    
    // loop over elements
    Real_wp tau, tau_primitive;
    Complexe k2 = square(kwave), coef_k2, inv_zeta, coef_radius;
    SetComplexZero(coef_k2);
    int num_point = 0;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	SetComplexOne(inv_zeta);
	SetComplexOne(zeta);
	SetComplexZero(tau); SetComplexZero(tau_primitive);
	const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(i);
	
	// number of integration points
	int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
	mat_mass(i).Reallocate(nb_points_quadrature);
	mat_sigma(i).Reallocate(nb_points_quadrature);
	mat_flow(i).Reallocate(nb_points_quadrature);
        mat_massV(i).Reallocate(nb_points_quadrature);
	if (this->mixed_formulation)
	  {
	    mat_sigmaV(i).Reallocate(nb_points_quadrature);
	    if (this->type_coordinate != VarProblem_1D::CARTESIAN)
	      mat_mass_radius(i).Reallocate(nb_points_quadrature);
	  }
        else if (this->mesh.Element(i).IsPML())
	  {
	    mat_sigmaV(i).Reallocate(nb_points_quadrature);
	    
	    if (this->calcul_enveloppe)
	      mat_massV(i).Reallocate(nb_points_quadrature);
	  }
	
        if ((this->type_coordinate == VarProblem_1D::SPHERICAL)
	    || (this->type_coordinate == VarProblem_1D::POLAR))
          mat_mass_mu(i).Reallocate(nb_points_quadrature);
        
	// We get the two extremities of the edge
	this->mesh.GetVerticesElement(i, s);
	
	// length of the edge
	Real_wp Ji = abs(s(1)-s(0));
	
	// physical coefficients rho and mu
	ref = this->mesh.Element(i).GetReference();
	SetComplexZero(sigma_v);
        
	for (int j = 0; j < nb_points_quadrature; j++)
	  {
	    Real_wp x = Fb.Points(j)*s(1) + (1.0-Fb.Points(j))*s(0);
	    CoorDofs(this->mesh.GetNumberDof(i, j)) = x;
	    if (this->mixed_formulation)
	      CoorDofs(this->mesh.GetNbDof() + this->OffsetDofV(i) + j) = x;
	    
            SetComplexOne(coef_radius);
            
	    GetPhysicalCoefficient(rho, mu, sigma, v0, i, j+1, ref);
	    invMu = Real_wp(1)/mu;

	    Complexe rtilde; SetComplexReal(x, rtilde);
            //DISP(i); DISP(j); DISP(rho); DISP(mu); DISP(invMu); DISP(sigma);
	    //DISP(invMu); DISP(x); DISP(rtilde);

	    if (this->mesh.Element(i).IsPML())
	      {		
		this->pml_damping.GetDampingTau(x, tau, tau_primitive);
                
                if (sizeof(Complexe) == sizeof(Real_wp))
                  SetComplexOne(zeta);
                else
                  {
		    if (m_iomega == zero)
		      SetComplexOne(zeta);
		    else
		      {
			zeta = Real_wp(1) + tau/m_iomega;
			rtilde = x + tau_primitive / m_iomega;
			inv_zeta = Real_wp(1) / zeta;
		      }
		  }
		
		sigma_v = tau*invMu;
	      }

	    if (this->type_coordinate == VarProblem_1D::POLAR)
              coef_radius = rtilde;
            else if (this->type_coordinate == VarProblem_1D::SPHERICAL)
              coef_radius = rtilde*rtilde;
            	    
	    if (mixed_formulation)
              {
                coef_mass = m_iomega*rho;
                coef_damping = rho*tau + sigma*zeta;
		if (this->type_coordinate != VarProblem_1D::CARTESIAN)
		  mat_mass_radius(i)(j) = coef_radius;
              }
            else
              {
                coef_mass = m_omega2*rho*zeta;
                coef_damping = m_iomega*sigma*zeta;
		coef_k2 = k2*inv_zeta;
              }
            
            mat_mass(i)(j) = Ji*coef_mass*Fb.Weights(j)*coef_radius;
            mat_sigma(i)(j) = Ji*coef_damping*Fb.Weights(j)*coef_radius;
	    if (helio_scaling_model)
	      {
		// term r^2 (rho' / (rho zeta) ) u'
		mat_flow(i)(j) = v0*inv_zeta*Fb.Weights(j)*coef_radius;
	      }
	    else
	      SetComplexZero(mat_flow(i)(j));

	    if (mixed_formulation)
	      {
                coef_mass = m_iomega*invMu;
                coef_damping = sigma_v;
		mat_massV(i)(j) = Ji*coef_mass*Fb.Weights(j)*coef_radius;
		mat_sigmaV(i)(j) = Ji*coef_damping*Fb.Weights(j)*coef_radius;
	      }
            else if (this->mesh.Element(i).IsPML())
              {
		mat_sigmaV(i)(j) = Real_wp(1)/(invMu*Ji)*Fb.Weights(j)*inv_zeta*coef_radius;
                if (this->type_coordinate == VarProblem_1D::POLAR)
		  mat_mass_mu(i)(j) = Ji/(invMu*inv_zeta*Fb.Weights(j)*coef_radius);
		else if (this->type_coordinate == VarProblem_1D::SPHERICAL)
                  mat_mass_mu(i)(j) = Ji/invMu*inv_zeta*Fb.Weights(j);
	      }		
            else
              {
                mat_massV(i)(j) = 1.0/(invMu*Ji)*Fb.Weights(j)*coef_radius;
                if (this->type_coordinate == VarProblem_1D::POLAR)
                  mat_mass_mu(i)(j) = Ji/(invMu*Fb.Weights(j)*coef_radius);
		else if (this->type_coordinate == VarProblem_1D::SPHERICAL)
                  mat_mass_mu(i)(j) = Ji/invMu*Fb.Weights(j);
              }
	    
	    if (this->calcul_enveloppe)
	      {
		// mass part zeta (-\omega^2 - i \sigma \omega) rho + k(omega)^2  
                //  -i omega / (-i omega + tau)
		// should be equal to zero if k(omega) is well chosen and outside PML
		mat_mass(i)(j) += Ji/invMu*coef_k2*Fb.Weights(j);
		//DISP(i); DISP(j); DISP(m_omega2); DISP(k2); DISP(mat_mass(i)(j));
		
		if (this->mesh.Element(i).IsPML())
		  {
		    // gradient part i k(omega) m_iomega/ (m_iomega+tau)  stored in mat_massV
		    mat_massV(i)(j) = ikwave/invMu*inv_zeta*Fb.Weights(j);
		    
		    // stiffness part is not modified - i omega / (-i omega + tau)
		  }
	      }
	    
	    num_point++;
	  }
      }
  }
  
  
  //! computation of mass matrix (assuming that this matrix is diagonal)
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::GetMassMatrix(Vector<Complexe>& Dh)
  {
    Dh.Reallocate(this->nodl);
    Dh.Fill(0);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      for (int j = 0; j < mat_mass(i).GetM(); j++)
	Dh(this->mesh.GetNumberDof(i, j)) += mat_mass(i)(j);

    int Nvol = this->mesh.GetNbDof();
    if (this->mixed_formulation)
      for (int i = 0; i < this->mesh.GetNbElt(); i++)
	for (int j = 0; j < mat_massV(i).GetM(); j++)
	  Dh(Nvol + this->OffsetDofV(i) + j) += mat_massV(i)(j);
    
    this->AddDomains(Dh);
  }


  //! computation of damping matrix (assuming that this matrix is diagonal)
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::GetDampingMatrix(Vector<Complexe>& Dh)
  {
    Dh.Reallocate(this->nodl);
    Dh.Fill(0);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      for (int j = 0; j < mat_mass(i).GetM(); j++)
	Dh(this->mesh.GetNumberDof(i, j)) += mat_sigma(i)(j);

    int Nvol = this->mesh.GetNbDof();
    if (this->mixed_formulation)
      for (int i = 0; i < this->mesh.GetNbElt(); i++)
	for (int j = 0; j < mat_massV(i).GetM(); j++)
	  Dh(Nvol + this->OffsetDofV(i) + j) += mat_sigmaV(i)(j);   
    
    Complexe impedance = sqrt(mu0*rho0);
    Complexe m_iomega, ikwave;
    GetMiomega(m_iomega); GetIkwave(ikwave);
    if (this->mixed_formulation)
      SetComplexOne(m_iomega);
    
    if (this->condition_left_side == BoundaryConditionEnum::LINE_ABSORBING)
      {
	Dh(0) += m_iomega*impedance;
	//if (calcul_enveloppe)
	//Dh(0) += ikwave;
      }
    
    if (this->condition_right_side == BoundaryConditionEnum::LINE_ABSORBING)
      {
	Dh(Nvol-1) += m_iomega*impedance;
	//if (calcul_enveloppe)
	//Dh(Nvol-1) += ikwave;
      }

    this->AddDomains(Dh);
  }
  
  
  //! computation of elementary matrix on element iquad
  template<class Complexe> template<class T>
  void VarHelmholtz_1D<Complexe>
  ::ComputeElementaryMatrixGen(int iquad, IVect& num_dof,
			       VirtualMatrix<T>& mat_interac,
			       const GlobalGenericMatrix<T>& nat_mat)
  {
    const ElementReference<Dimension1, 1>& Fb = this->GetReferenceElementH1(iquad);
    int nb_dof_edge = Fb.GetNbDof();

    // extremities of the edge
    VectReal_wp s(2);
    this->mesh.GetVerticesElement(iquad, s);
    
    int nb_dof_loc = nb_dof_edge;
    if (this->mixed_formulation)
      nb_dof_loc *= 2;
    
    mat_interac.Reallocate(nb_dof_loc, nb_dof_loc);
    num_dof.Reallocate(nb_dof_loc);
    mat_interac.Zero();;
    int ref = this->mesh.Element(iquad).GetReference();
    Complexe mu, zero;
    SetComplexZero(zero);
    
    Complexe ikwave; this->GetIkwave(ikwave);
    
    if (this->mixed_formulation)
      {
	int Nvol = this->mesh.GetNbDof();
	// mass and damping part
	for (int i = 0; i < nb_dof_edge; i++)
	  {
	    num_dof(i) = this->mesh.GetNumberDof(iquad, i);
	    num_dof(nb_dof_edge+i) = Nvol + this->OffsetDofV(iquad) + i;
	    
	    mat_interac.SetEntry(i, i, this->mat_mass(iquad)(i)*nat_mat.GetCoefMass()
                                 + this->mat_sigma(iquad)(i)*nat_mat.GetCoefDamping());
	    
	    mat_interac.SetEntry(nb_dof_edge+i, nb_dof_edge+i,
                                 this->mat_massV(iquad)(i)*nat_mat.GetCoefMass()
                                 + this->mat_sigmaV(iquad)(i)*nat_mat.GetCoefDamping());	    
	  }
	
	// stiffness part
	if (type_coordinate == VarProblem_1D::CARTESIAN)
	  {
	    const Matrix<Real_wp>& gradient_matrix = Fb.GetGradientMatrix();
	    for (int i = 0; i < nb_dof_edge; i++)
	      for (int j = 0; j < nb_dof_edge; j++)
		{
		  mat_interac.SetEntry(i, nb_dof_edge+j, gradient_matrix(j, i)*nat_mat.GetCoefStiffness());
		  mat_interac.SetEntry(nb_dof_edge+j, i, -gradient_matrix(j, i)*nat_mat.GetCoefStiffness());
		}
	  }
	else if (Fb.LumpedMassMatrix())
	  {
	    const Matrix<Real_wp>& gradient_matrix = Fb.GetGradientMatrix();
	    //DISP(iquad); DISP(this->mat_mass_radius(iquad));
	    for (int i = 0; i < nb_dof_edge; i++)
	      for (int j = 0; j < nb_dof_edge; j++)
		{
		  mat_interac.SetEntry(i, nb_dof_edge+j, gradient_matrix(j, i)*this->mat_mass_radius(iquad)(j)*nat_mat.GetCoefStiffness());
		  mat_interac.SetEntry(nb_dof_edge+j, i, -gradient_matrix(j, i)*this->mat_mass_radius(iquad)(j)*nat_mat.GetCoefStiffness());
		}
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
      }
    else
      {
	Vector<T> coef_mass(nb_dof_edge), coef_D(nb_dof_edge), coef_E(nb_dof_edge),
	  coef_stiff(nb_dof_edge);
	
	TinyVector<bool, 4> null_term;
	
	null_term.Fill(false);
	coef_D.Zero(); coef_E.Zero(); coef_stiff.Zero();
	
	// dofs of the edge
	for (int i = 0; i < nb_dof_edge; i++)
	  {
	    num_dof(i) = this->mesh.GetNumberDof(iquad, i);
	    
	    // mass matrix coefficient
	    coef_mass(i) = this->mat_mass(iquad)(i)*nat_mat.GetCoefMass()
	      + this->mat_sigma(iquad)(i)*nat_mat.GetCoefDamping();
	  }
        
	if (nat_mat.GetCoefStiffness() != 0)
	  {
            // we add gradient contribution
	    for (int j = 0; j < coef_E.GetM(); j++)
	      coef_E(j) = this->mat_flow(iquad)(j)*nat_mat.GetCoefStiffness();
	    
            // and stiffness contribution
	    if (this->mesh.Element(iquad).IsPML())
              {
                for (int j = 0; j < coef_stiff.GetM(); j++)
                  coef_stiff(j) = this->mat_sigmaV(iquad)(j)*nat_mat.GetCoefStiffness();
              }
            else
              {
                if ((ref_mu(ref).IsVarying()) || (this->type_coordinate != VarProblem_1D::CARTESIAN))
                  {
                    for (int j = 0; j < coef_stiff.GetM(); j++)
                      coef_stiff(j) = this->mat_massV(iquad)(j)*nat_mat.GetCoefStiffness();
                  }
                else
                  {
                    T coef = ref_mu(ref).GetConstant()/(s(1)-s(0))*nat_mat.GetCoefStiffness();
		    for (int j = 0; j < coef_stiff.GetM(); j++)
		      coef_stiff(j) = coef*Fb.Weights(j);		    
                  }
              }

	    // gradient contribution for enveloppe
	    if (this->calcul_enveloppe)
	      {
		if (this->mesh.Element(iquad).IsPML())
		  {
		    for (int j = 0; j < coef_D.GetM(); j++)
		      {
			coef_D(j) = this->mat_massV(iquad)(j)*nat_mat.GetCoefStiffness();
			coef_E(j) -= coef_D(j);
		      }
		  }
		else
		  {
		    mu = ref_mu(ref).GetConstant();
		    T coef = mu*ikwave*nat_mat.GetCoefStiffness();
		    for (int j = 0; j < coef_D.GetM(); j++)
		      {
			coef_D(j) = coef*Fb.Weights(j);
			coef_E(j) -= coef_D(j);
		      }
		  }
	      }
	  }
	else
	  {
	    null_term(1) = true; null_term(2) = true; null_term(3) = true;
	  }

	Fb.AddVariableElemMatrix(0, 0, coef_mass, coef_stiff, coef_D, coef_E, null_term, mat_interac);
      }
  }
  
  
  //! initialization of physical indexes
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::InitIndices(int n)
  {
    ref_rho.Reallocate(n);
    ref_mu.Reallocate(n);
    ref_sigma.Reallocate(n);
    ref_true_rho.Reallocate(n);
    log_rho_read.Reallocate(n);
    log_rho_read.Fill(false);

    Complexe one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        ref_rho(i).SetConstant(rho0);
        ref_mu(i).SetConstant(mu0);
        ref_sigma(i).SetConstant(zero);
        ref_true_rho(i).SetConstant(one);
      }
  }
  
  
  //! modification of physical indexes according to parameters of data file
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::SetIndices(int i, const VectString& parameters)
  {
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

    int nb = 1;
    ref_rho(i).SetInputData(parameters, nb);
    ref_mu(i).SetInputData(parameters, nb);
    ref_sigma(i).SetInputData(parameters, nb);
  }
  
  
  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

    int nb = 2;
    if (name_media == "rho")
      ref_rho(i).SetInputData(parameters, nb);
    else if (name_media == "mu")
      ref_mu(i).SetInputData(parameters, nb);
    else if (name_media == "sigma")
      ref_sigma(i).SetInputData(parameters, nb);
    else if (name_media == "true_rho")
      ref_true_rho(i).SetInputData(parameters, nb);
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }
  
  
  //! returns the velocity of waves in the physical media of reference ref
  template<class Complexe>
  Real_wp VarHelmholtz_1D<Complexe>::GetVelocityOfMedia(int ref) const
  {
    if (helio_scaling_model)
      return abs(this->ref_mu(ref).GetConstant());
    
    return sqrt(abs(this->ref_mu(ref).GetConstant()/this->ref_rho(ref).GetConstant()));
  }


  //! computation of right hand side
  template<class T>
  void VarHelmholtz_1D<T>::ComputeRightHandSide(Vector<T>& bsrc)
  {
    bsrc.Reallocate(this->GetNbDof()); bsrc.Zero();
    
    // inhomogeneous Dirichlet
    //if (condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
    bsrc(0) = left_dirichlet_condition;
    
    //if (condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
    bsrc(this->mesh.GetNbDof()-1) = right_dirichlet_condition;
  }

  
  //! launching all the simulation
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::RunAll(const string& input_file, const string& name_elt)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    ReadInputFile(input_file, *this);
    
    if ((this->print_level >= 0) && (rank_proc == 0))
      cout << "The data file has been read" << endl;
    
    // initialisation calcul
    this->ComputeMeshAndFiniteElement(name_elt);    
    this->ComputePhysicalCoefficients();
    this->ComputeMassMatrix();

    GlobalGenericMatrix<Complexe> nat_mat;
    int N = this->GetNbDof();
    if ((this->print_level >= 0) && (rank_proc == 0))
      {
	cout << "Order = " << this->mesh.GetOrder() << endl;
	cout << "Number of dofs = " << N << endl;
      }

    // calcul et factorisation de la matrice element fini
    if (this->mixed_formulation)
      {
	Matrix<Complexe, General, ArrayRowSparse> mat_direct;	
	mat_direct.Reallocate(N, N);

	this->AddMatrixFEM(mat_direct, nat_mat);	
	this->AddBoundaryTerms(mat_direct, nat_mat);

	mat_direct.WriteText("mat_first.dat");
	
	// calcul du second membre
	Vector<Complexe> rhs(N), x_sol(N);
	this->ComputeRightHandSide(rhs);

	// calcul de la solution
	All_MatrixLU<Complexe> mat_lu;
	mat_lu.Factorize(mat_direct);
	
	x_sol = rhs;
	mat_lu.Solve(x_sol);
	
	// sortie
	this->WriteDatas(x_sol);

	// on sort aussi l'inconnue vectorielle
        DISP(mesh.GetNbElt());
	DISP(N); DISP(mesh.GetNbDof());
	int Nv = N - this->mesh.GetNbDof();
	
	Vector<Complexe> V(Nv);
	for (int i = 0; i < Nv; i++)
	  V(i) = x_sol(this->mesh.GetNbDof() + i);
	
	int Np = this->var_section.GetNbPointsGrid();
	Vector<Complexe> output_vector(Np);
	
	this->ComputeInterpolationU(V, this->var_section, output_vector, true, false);

	int double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
	bool ascii = false;
	
	GridInterpolationFull<Dimension1> grid;
	int type = GridInterpolationFull<Dimension1>::LINE;
	grid.Init(type, this->var_section.GetXmin(), this->var_section.GetXmax(), Np);
	
	string name_file = GetBaseString(this->output_section_param.GetTotalFieldFile()) + "_V.dat";
	WriteMatlab(output_vector, grid, name_file, double_prec, ascii);
      }
    else
      {
	Matrix<Complexe, General, BandedCol> mat_direct;
	
	mat_direct.Reallocate(N, N, this->mesh.GetOrder(), this->mesh.GetOrder());
	this->AddMatrixFEM(mat_direct, nat_mat);
	
	this->AddBoundaryTerms(mat_direct, nat_mat);
	mat_direct.WriteText("mat.dat");
    
	Vector<int> pivot;
	mat_direct.Factorize(pivot);
	
	// calcul du second membre
	Vector<Complexe> rhs(N), x_sol(N);
	this->ComputeRightHandSide(rhs);
	
	// calcul de la solution
	x_sol = rhs;
	mat_direct.Solve(pivot, x_sol);
	
	// sortie
	this->WriteDatas(x_sol);
      }
  }
  
  
  //! treatment of absorbing boundary condition
  template<class Complexe> template<class T>
  void VarHelmholtz_1D<Complexe>
  ::AddBoundaryTerms(VirtualMatrix<T>& mat_sp, GlobalGenericMatrix<T>& nat_mat)
  {
    // condition on the left extremity
    int ref = this->mesh.Element(0).GetReference();
    Complexe m_iomega; GetMiomega(m_iomega);
    Complexe rho, mu, sigma, v0;
    GetPhysicalCoefficient(rho, mu, sigma, v0, 0, 0, ref);
    rho += sigma/m_iomega;
        
    Complexe impedance = sqrt(mu*rho);
    if (this->mixed_formulation)
      SetComplexOne(m_iomega);
    
    Complexe coef_impedance = m_iomega*impedance;
    if (this->type_coordinate == VarProblem_1D::SPHERICAL)
      coef_impedance *= this->xmin_*this->xmin_;
    else if (this->type_coordinate == VarProblem_1D::POLAR)
      coef_impedance *= this->xmin_;
    
    if (this->condition_left_side == BoundaryConditionEnum::LINE_ABSORBING)
      mat_sp.AddInteraction(0, 0, coef_impedance*nat_mat.GetCoefDamping());
    
    coef_impedance = mu*beta_impedance_left;
    if (this->type_coordinate == VarProblem_1D::SPHERICAL)
      coef_impedance *= this->xmin_*this->xmin_;
    else if (this->type_coordinate == VarProblem_1D::POLAR)
      coef_impedance *= this->xmin_;
    
    if (this->condition_left_side == BoundaryConditionEnum::LINE_IMPEDANCE)
      mat_sp.AddInteraction(0, 0, coef_impedance*nat_mat.GetCoefDamping());

    // condition on the right extremity
    ref = this->mesh.Element(this->mesh.GetNbElt()-1).GetReference();
    int order = this->GetReferenceElementH1(this->mesh.GetNbElt()-1).GetNbPointsQuadratureInside();
    GetPhysicalCoefficient(rho, mu, sigma, v0, this->mesh.GetNbElt()-1, order+1, ref);

    rho += sigma/m_iomega;    
    int Nvol = this->mesh.GetNbDof();
    impedance = sqrt(mu*rho);
    coef_impedance = m_iomega*impedance;
    if (this->type_coordinate == VarProblem_1D::SPHERICAL)
      coef_impedance *= this->xmax_*this->xmax_;
    else if (this->type_coordinate == VarProblem_1D::POLAR)
      coef_impedance *= this->xmax_;

    if (this->condition_right_side == BoundaryConditionEnum::LINE_ABSORBING)
      mat_sp.AddInteraction(Nvol-1, Nvol-1, coef_impedance*nat_mat.GetCoefDamping());

    coef_impedance = mu*beta_impedance_right;
    if (this->type_coordinate == VarProblem_1D::SPHERICAL)
      coef_impedance *= this->xmax_*this->xmax_;
    else if (this->type_coordinate == VarProblem_1D::POLAR)
      coef_impedance *= this->xmax_;
    
    if (this->condition_right_side == BoundaryConditionEnum::LINE_IMPEDANCE)
      mat_sp.AddInteraction(Nvol-1, Nvol-1, coef_impedance*nat_mat.GetCoefDamping());

    // enforcing Dirichlet condition (matrix is assumed to be non-symmetric)
    T one; SetComplexOne(one);
    if (condition_left_side == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	mat_sp.ClearRow(0);
	mat_sp.AddInteraction(0, 0, one);
      }

    if (condition_right_side == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	mat_sp.ClearRow(Nvol-1);
	mat_sp.AddInteraction(Nvol-1, Nvol-1, one);
      }
  }
  

  //! reads the impedance from parameters
  template<class Complexe>
  void VarHelmholtz_1D<Complexe>::ReadImpedance(int& nb, const VectString& parameters, int side)
  {
    nb++;  
    if (side == 0)
      beta_impedance_left = to_num<Complexe>(parameters(nb));
    else
      beta_impedance_right = to_num<Complexe>(parameters(nb));
  }
  
} // end namespace

#define MONTJOIE_FILE_HELMHOLTZ_1D_CXX
#endif
