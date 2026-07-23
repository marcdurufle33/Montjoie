#ifndef MONTJOIE_FILE_MAXWELL_AXISYM_HARMONIC_CXX

namespace Montjoie
{

  //! additionnal input datas for axisymmetric computation
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {    
    VarAxisymProblem::SetInputData(description_field, parameters);
    
    if (!description_field.compare("OutputHy"))
      {
        int position_word = -1;
        for (int i = 0; i < parameters.GetM(); i++)
          if (isalpha(parameters(i)[0]))
	    {
              position_word = i;
              break;
            }
        
        if (position_word > 0)
          {
            ref_outputJ.Reallocate(position_word);
            for (int i = 0; i < position_word; i++)
              ref_outputJ(i) = to_num<int>(parameters(i));
	    
	    if (parameters.GetM() <= position_word+1)
	      {
		cout << "In SetInputData of HarmonicMaxwell_Axi" << endl;
		cout << "OutputHy needs more parameters, for instance :" << endl;
		cout << "OutputHy = ref PARAM nom_fichier subdiv" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            name_file_outputJ = parameters(position_word);
            nb_subdiv_outputJ = to_num<int>(parameters(position_word+1));
            if (parameters.GetM() > position_word+2)
              output_Hy_teta = false;
            else
              output_Hy_teta = true;
          }
      }    
  }
  
  
  //! allocation of arrays containing physical properties
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>::InitIndices(int n)
  {
    if (PhysicalConstant::adimensionalization != PhysicalConstant::ADIM_NO)
      {
        epsilon0 = 1.0;
        mu0 = 1.0;
        invMu0 = 1.0;
      }
    else
      {
	// use of physical units
	epsilon0 = PhysicalConstant::epsilon0_permittivity;
	
        mu0 = PhysicalConstant::mu0_permeability;
        invMu0 = 1.0/mu0;
      }
    
    ref_epsilon.Reallocate(n);
    ref_mu.Reallocate(n);
    ref_invMu.Reallocate(n);
    ref_invEpsilon.Reallocate(n);
    ref_sigma.Reallocate(n);
    Complexe zero; SetComplexZero(zero);
    for (int i = 0; i < n; i++)
      {
        ref_epsilon(i).SetDiagonal(epsilon0);
        ref_invEpsilon(i).SetDiagonal(1.0/epsilon0);
        ref_sigma(i).SetConstant(zero);
        ref_mu(i).SetDiagonal(mu0);
        ref_invMu(i).SetDiagonal(invMu0);
      }
  }
  
  
  //! set physical properties of media with reference i
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>
  ::SetIndices(int i, const VectString& parameters)
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
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      ref_sigma(i).Mlt(PhysicalConstant::impedance0);

    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
      {
        ref_epsilon(i).Mlt(PhysicalConstant::epsilon0_permittivity);
        ref_mu(i).Mlt(PhysicalConstant::mu0_permeability);
      }
    
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
  void HarmonicMaxwell_Axi<Complexe>
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
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
          ref_epsilon(i).Mlt(PhysicalConstant::epsilon0_permittivity);

        ref_epsilon(i).Mlt(epsilon0);
      }
    else if (name_media == "mu")
      {
        ref_mu(i).SetInputData(nb, parameters, parameters(0));
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
          ref_mu(i).Mlt(PhysicalConstant::mu0_permeability);

        ref_mu(i).Mlt(mu0);
      }
    else if (name_media == "sigma")
      {
        ref_sigma(i).SetInputData(nb, parameters, parameters(0));
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          ref_sigma(i).Mlt(PhysicalConstant::impedance0);
      }
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }


  //! returns the name associated with the physical index num
  template<class Complexe>
  string HarmonicMaxwell_Axi<Complexe>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("epsilon");
      case 1: return string("mu");
      case 2: return string("sigma");
      }

    return string();
  }  
    
  
  //! fills arrays with the varying indices present in the current problem
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>::
  GetVaryingIndex(Vector<PhysicalVaryingMedia<Dimension2, Complexe>* >& rho_complex, 
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
  

  //! computation of invEpsilon and invMu from values of Epsilon and Mu
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>::FinalizeComputationVaryingIndices()
  {
    // computation of inverses    
    for (int ref = 1; ref < ref_mu.GetM(); ref++)
      {
        ref_invEpsilon(ref) = ref_epsilon(ref);
        ref_invEpsilon(ref).GetInverse();

        ref_invMu(ref) = ref_mu(ref);
        ref_invMu(ref).GetInverse();
      }
  }

  
  //! mu and epsilon are modified in PML layers
  template<>
  void HarmonicMaxwell_Axi<Complex_wp>
  ::ModifyPMLCoefficient(TinyMatrix<Complex_wp, Symmetric, 2, 2>& epsilon_rz,
			 Complex_wp& epsilon_teta,
			 TinyMatrix<Complex_wp, Symmetric, 2, 2>& mu_rz,
			 Complex_wp& mu_teta, Complex_wp& dr, Complex_wp& dz, int i1, int j) const
  {
    if (var_problem.FirstOrderFormulation())
      {
	Complex_wp coef_om = Iwp/var_problem.GetOmega();
	dr = 1.0 + coef_om*var_boundary.GetTauPML(i1, j, 0);
	dz = 1.0 + coef_om*var_boundary.GetTauPML(i1, j, 1);
      }
    else
      {
	dr = 1.0/var_boundary.GetTauPML(i1, j, 0);
	dz = 1.0/var_boundary.GetTauPML(i1, j, 1);
      }
    
    epsilon_rz(0, 0) *= dz/dr;
    epsilon_rz(1, 1) *= dr/dz;

    mu_rz(0, 0) *= dz/dr;
    mu_rz(1, 1) *= dr/dz;
    
    epsilon_teta *= dr*dz;
    mu_teta *= dr*dz;
  }

  
  //! no change in real domain
  template<>
  void HarmonicMaxwell_Axi<Real_wp>
  ::ModifyPMLCoefficient(TinyMatrix<Real_wp, Symmetric, 2, 2>& epsilon_rz,
			 Real_wp& epsilon_teta,
			 TinyMatrix<Real_wp, Symmetric, 2, 2>& mu_rz,
			 Real_wp& mu_teta, Real_wp& dr, Real_wp& dz, int i1, int j) const
  {
  }
  
  
  //! allocation of arrays needed for the computation of elementary matrices
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>::AllocateMassMatrices()
  {
    Glob_rtilde.Reallocate(var_problem.mesh.GetNbElt());
  }
  
  
  //! computation of geometric quantities needed for the expression of elementary matrices
  template<class Complexe>
  void HarmonicMaxwell_Axi<Complexe>
  ::ComputeLocalMassMatrix(int i)
  {
    int N = var_problem.Glob_PointsQuadrature(i).GetM();
    Glob_rtilde(i).Reallocate(N);
    if (!var_problem.InsidePML(i))
      {
        for (int j = 0; j < N; j++)
          Glob_rtilde(i)(j) = var_problem.Glob_PointsQuadrature(i)(j)(0);
      }
    else
      {
	Complexe rtilde;
        int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
        for (int j = 0; j < N; j++)
          {
            rtilde = var_boundary.GetPrimitiveTauPML(i1, j, 0);
            Glob_rtilde(i)(j) = rtilde;
	  }
      }    
  }
  

  /****************
   * ImpedanceABC *
   ****************/
  
  
  //! evaluation of impedance coefficient
  void ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
                       const SetPoints<Dimension2>& Pts,
                       const SetMatrices<Dimension2>& Mat)
  {
    Real_wp radius = Pts.GetPointQuadratureBoundary(k)(0);
    Complex_wp epsilon = var_maxwell.ref_epsilon(ref)(0,0)
      + Iwp*var_maxwell.ref_sigma(ref).GetConstant()/var_maxwell.GetOmega();
    
    Complex_wp mu = var_maxwell.ref_mu(ref)(0, 0);
    Complex_wp feval = -Iwp*radius*sqrt(epsilon/mu);
    if (var_maxwell.modified_formulation)
      feval *= var_maxwell.GetOmega();
    
    if (k == 0)
      {
        VecNormale.Reallocate(Mat.GetNbPointsQuadratureBoundary());
	VecImpedance.Reallocate(VecNormale.GetM());
	VecRadius.Reallocate(VecNormale.GetM());
        for (int k2 = 0; k2 < Mat.GetNbPointsQuadratureBoundary(); k2++)
          {
	    VecNormale(k2) = Mat.GetNormaleQuadratureBoundary(k2);
	    VecRadius(k2) = Pts.GetPointQuadratureBoundary(k2)(0);
	  }
      }
    
    VecImpedance(k) = feval;
    modified_formulation = var_maxwell.modified_formulation;
    number_mode = var_maxwell.GetCurrentModeNumber();
  }
  
  
  //! application of impedance of first-order absorbing boundary condition 
  void ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  ApplyImpedancePhi_H1(int n, int j, int offset,
                       const TinyVector<Real_wp, 1>& phi_boundary,
                       const TinyVector<Real_wp, 2>& dphi_boundary, VectComplex_wp& feval_phi)
  {
    feval_phi.Zero();
    if ((modified_formulation) && (number_mode != 0))
      {
	Real_wp m2 = square(number_mode);
	Real_wp nr = VecNormale(j)(0), nz = VecNormale(j)(1);
	Real_wp r = VecRadius(j);
	Real_wp psi = phi_boundary(0);
	Real_wp rotPsi_times_n = nr*dphi_boundary(0) + nz*dphi_boundary(1);
	Complex_wp vloc = (psi*nz - r*rotPsi_times_n);	
	vloc *= VecImpedance(j);

	feval_phi(0) = r*nz*vloc;
        feval_phi(1) = -r*nr*vloc;
	feval_phi(2) = nz*vloc + VecImpedance(j)*m2*psi;
      }
    else
      {
	Complex_wp feval = VecImpedance(j);
	feval_phi(2) = feval*phi_boundary(0);
      }
  }


  //! application of impedance of first-order absorbing boundary condition 
  void ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  ApplyImpedancePhi_Hcurl(int n, int j, int offset,
                       const TinyVector<Real_wp, 2>& phi_boundary,
                          const TinyVector<Real_wp, 1>& dphi_boundary, VectComplex_wp& feval_phi)
  {
    feval_phi.Zero();
    Real_wp phir = phi_boundary(0);
    Real_wp phiz = phi_boundary(1);
    Real_wp phi_times_n = phir*VecNormale(j)(1) - phiz*VecNormale(j)(0);
    if ((modified_formulation) && (number_mode != 0))
      {
	//Real_wp m2 = square(number_mode);
	Real_wp nr = VecNormale(j)(0), nz = VecNormale(j)(1);
	Real_wp r = VecRadius(j);
	Complex_wp vloc = r*phi_times_n;	
	vloc *= VecImpedance(j);

	feval_phi(0) = r*nz*vloc;
	feval_phi(1) = -r*nr*vloc;
        feval_phi(2) = nz*vloc;
      }
    else
      {
	Complex_wp feval = VecImpedance(j);
	feval_phi(0) = feval*VecNormale(j)(1)*phi_times_n;
	feval_phi(1) = -feval*VecNormale(j)(0)*phi_times_n;
      }

  }
  
  
  //! terms for modified formulation
  void ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  ApplyImpedanceGrad(int n, int j, int, const TinyVector<Real_wp, 1>& phi_boundary,
		     const TinyVector<Real_wp, 2>& dphi_boundary, VectComplex_wp& feval_dphi)
  {
    feval_dphi.Zero();
    if ((modified_formulation) && (number_mode != 0))
      {
	Real_wp nr = VecNormale(j)(0), nz = VecNormale(j)(1);
	Real_wp psi = phi_boundary(0);
	Real_wp rotPsi_times_n = nr*dphi_boundary(0) + nz*dphi_boundary(1);
	Real_wp r = VecRadius(j);
	Complex_wp vloc = -(psi*nz - r*rotPsi_times_n);
	vloc *= VecImpedance(j);
	feval_dphi(1) = r*nr*vloc;
	feval_dphi(2) = r*nz*vloc;
      }
  }

  //! terms for modified formulation
  void ImpedanceABC<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  ApplyImpedanceCurl(int n, int j, int offset, const TinyVector<Real_wp, 2>& phi_boundary,
		     const TinyVector<Real_wp, 1>& dphi_boundary, VectComplex_wp& feval_dphi)
  {
    feval_dphi.Zero();
    if ((modified_formulation) && (number_mode != 0))
      {
	Real_wp nr = VecNormale(j)(0), nz = VecNormale(j)(1);
	Real_wp phir = phi_boundary(0);
	Real_wp phiz = phi_boundary(1);
	Real_wp phi_times_n = phir*nz - phiz*nr;
	Real_wp r = VecRadius(j);
	Complex_wp vloc = -r*phi_times_n;
	vloc *= VecImpedance(j);
	feval_dphi(1) = r*nr*vloc;
	feval_dphi(2) = r*nz*vloc;
      }
  }
  
  
  /********************
   * ImpedanceGeneric *
   ********************/


  //! evaluation of generic impedance
  void ImpedanceGeneric<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref_d,
                       const SetPoints<Dimension2>& Pts, const SetMatrices<Dimension2>& Mat)
  {
    abort();
  }


  //! multiplication by impedance
  void ImpedanceGeneric<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi_boundary,
                    const TinyVector<Real_wp, 2>& dphi_boundary, VectComplex_wp& feval_phi)
  {
    abort();
  }


  //! multiplication by impedance
  void ImpedanceGeneric<Complex_wp, HarmonicMaxwellEquation_HcurlAxi>::
  ApplyImpedancePhi_Hcurl(int n, int j, int, const TinyVector<Real_wp, 2>& phi_boundary,
                    const TinyVector<Real_wp, 1>& dphi_boundary, VectComplex_wp& feval_phi)
  {
    abort();
  }

  
  /*******************
   * LDG formulation *
   *******************/
  
  
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
  void HarmonicMaxwellEquationAxiDG::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    FillZero(mass);
    TinyMatrix<Complex_wp, Symmetric, 3, 3> Eps
      = vars.ref_epsilon(ref).GetCoefficient(vars, i, j);
    
    Complex_wp Sigma = vars.ref_sigma(ref).GetCoefficient(vars, i, j);
    Eps(0, 0) += Iwp*Sigma/vars.GetOmega();
    Eps(1, 1) += Iwp*Sigma/vars.GetOmega();
    Eps(2, 2) += Iwp*Sigma/vars.GetOmega();
    
    Complex_wp coef = -vars.GetOmega()*vars.Glob_rtilde(i)(j)*nat_mat.GetCoefMass();
    mass(0, 0) = coef*Eps(0, 0); mass(0, 1) = coef*Eps(0, 1); mass(0, 2) = coef*Eps(0, 2);
    mass(1, 0) = coef*Eps(1, 0); mass(1, 1) = coef*Eps(1, 1); mass(1, 2) = coef*Eps(1, 2);
    mass(2, 0) = coef*Eps(2, 0); mass(2, 1) = coef*Eps(2, 1); mass(2, 2) = coef*Eps(2, 2);
    
    Eps = vars.ref_mu(ref).GetCoefficient(vars, i, j);
    mass(3, 3) = coef*Eps(0, 0); mass(3, 4) = coef*Eps(0, 1); mass(3, 5) = coef*Eps(0, 2);
    mass(4, 3) = coef*Eps(1, 0); mass(4, 4) = coef*Eps(1, 1); mass(4, 5) = coef*Eps(1, 2);
    mass(5, 3) = coef*Eps(2, 0); mass(5, 4) = coef*Eps(2, 1); mass(5, 5) = coef*Eps(2, 2);
    
    mass(0, 4) = Real_wp(vars.GetCurrentModeNumber())*nat_mat.GetCoefStiffness();
    mass(1, 3) = -mass(0, 4);
    mass(4, 0) = mass(0, 4);
    mass(3, 1) = mass(1, 3);
    
    mass(2, 4) = -nat_mat.GetCoefStiffness();
    mass(4, 2) = mass(2, 4);
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
  void HarmonicMaxwellEquationAxiDG::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    abort();
  }
  
  
  //! which derivatives to evaluate during matrix-vector product ?
  template<class TypeEquation, class T0, class Vector1>
  void HarmonicMaxwellEquationAxiDG::
  GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                      const GlobalGenericMatrix<T0>& nat_mat,
                      Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    for (int i = 0; i < 3; i++)
      {
        unknown_to_derive(i) = true;
        fct_test_to_derive(i) = true;
        unknown_to_derive(i+3) = false;
        fct_test_to_derive(i+3) = false;
      }
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
  void HarmonicMaxwellEquationAxiDG::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                   int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                   MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    Complex_wp coef = vars.Glob_rtilde(num_elem)(jloc)*nat_mat.GetCoefStiffness();
    Dgrad_phi(0, 5)(1) = coef;
    Dgrad_phi(1, 5)(0) = -coef;
    Dgrad_phi(2, 3)(1) = coef;
    Dgrad_phi(2, 4)(0) = -coef;
    
    Ephi_grad(5, 0)(1) = Dgrad_phi(0, 5)(1);
    Ephi_grad(5, 1)(0) = Dgrad_phi(1, 5)(0);
    Ephi_grad(3, 2)(1) = Dgrad_phi(2, 3)(1);
    Ephi_grad(4, 2)(0) = Dgrad_phi(2, 4)(0);
    
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
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void HarmonicMaxwellEquationAxiDG::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& dU, Vector2& V)
  {
    abort();
  }

  
  //! Applying the tensor E to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result E u
  */
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void HarmonicMaxwellEquationAxiDG::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                       int ref, Vector1& U, Vector2& dV)
  {
    abort();
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
  void HarmonicMaxwellEquationAxiDG::
  GetNabc(Matrix1& Nabc, R2& normale,
          int ref, int iquad, int k,const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
          const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    FillZero(Nabc);
    int cond = vars.mesh.GetBoundaryCondition(ref);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
        Complex_wp coef = vars.Glob_rtilde(iquad)(k)*nat_mat.GetCoefStiffness();
        Nabc(0, 5) = -coef*normale(1);
        Nabc(1, 5) = coef*normale(0);
        Nabc(2, 3) = -coef*normale(1);
        Nabc(2, 4) = coef*normale(0);
        
        Nabc(5, 0) = Nabc(0, 5);
        Nabc(5, 1) = Nabc(1, 5);
        Nabc(3, 2) = Nabc(2, 3);
        Nabc(4, 2) = Nabc(2, 4);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        Complex_wp coef = -Iwp*vars.Glob_rtilde(iquad)(k)*nat_mat.GetCoefStiffness();
        Nabc(0, 0) = coef*normale(1)*normale(1);
        Nabc(1, 1) = coef*normale(0)*normale(0);
        Nabc(0, 1) = -coef*normale(0)*normale(1);
        Nabc(1, 0) = Nabc(0, 1);
        Nabc(2, 2) = coef;
        
        Nabc(3, 3) = Nabc(0, 0);
        Nabc(4, 4) = Nabc(1, 1);
        Nabc(3, 4) = Nabc(0, 1);
        Nabc(4, 3) = Nabc(1, 0);
        Nabc(5, 5) = Nabc(2, 2);
        
      }
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
    \param[in] refd reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Vector1, class TypeEquation, class T0>
  void HarmonicMaxwellEquationAxiDG::
  MltNabc(R2& normale, int ref,
          const Vector1& Vn, Vector1& Un, int num_elem1, int k,
          const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
          const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    abort();
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
  template<class Matrix1, class GenericPb, class T0>
  void HarmonicMaxwellEquationAxiDG::
  GetPenalDG(Matrix1& Nabc, R2& normale,
             int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
             const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    FillZero(Nabc);
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
  template<class Vector1, class Vector2, class GenericPb, class T0>
  void HarmonicMaxwellEquationAxiDG::
  MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
             int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
             const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    abort();
  }
  
  
  /********************************
   * VarComputationRCS_MaxwellAxi *
   ********************************/

  
  void VarComputationRCS_MaxwellAxi::ComputeRCS(const VectReal_wp& U0)
  {
    cout << "not implemented" << endl;
    abort();
  }


  void VarComputationRCS_MaxwellAxi::ComputeRCS(const VectComplex_wp& U0)
  {
    cout << "not implemented" << endl;
    abort();
  }

  
  //! computation of radar cross section
  /*!
    \param[in] U0 solution vector
    \param[out] RCS_Vector far field 
   */
  void VarComputationRCS_MaxwellAxi
  ::ComputeRCS(const Matrix<Complex_wp, General, ColMajor>& U0, VectR3_Complex_wp& RCS_Vector)
  { 
#ifdef SELDON_WITH_MPI   
    int rank_proc; MPI_Comm_rank(var_maxwell.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    int nbPointsRCS = this->nb_angles_RCS;
    if (nbPointsRCS <= 0)
      return;
    
    if (this->type_rcs == this->MONOSTATIC_RCS)
      nbPointsRCS *= 2;
    
    // allocation is performed only for the first mode encountered
    // otherwise we take old values and add contribution of the current mode
    if ((RCS_Vector.GetM() == 0) && (rank_proc == 0))
      {
	RCS_Vector.Reallocate(nbPointsRCS);
	for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
	  RCS_Vector(num_angle).Fill(0);
      }
    
    if (this->type_rcs == this->MONOSTATIC_RCS)
      {
	// monostatic rcs
	ComputeMonostaticRCS(U0, RCS_Vector);
	
	if (rank_proc == 0)
	  WriteMonostaticRCS(RCS_Vector);
	
	return;
      }

    // otherwise bistatic rcs
    Vector<Complex_wp> U0_vec;
    U0_vec.SetData(U0.GetM(), U0.GetData());

    // computes E and H on the quadrature points
    VectR3_Complex_wp trace_En, trace_Hn;
    ComputeEnHn(U0_vec, trace_En, trace_Hn);
    U0_vec.Nullify();
    
    R3 ur; Real_wp teta;
    Real_wp step_angle = (this->last_angle_RCS-this->first_angle_RCS) / (nbPointsRCS-1);

    // updating the far field
    // loop on each angle
    R3_Complex_wp rcs_teta; VectComplex_wp ContribRcs(3*nbPointsRCS);
    for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
      {
	rcs_teta.Fill(0);
	ur = this->incident_direction(num_angle);
	
	AddContributionMode(num_angle, ur, trace_En, trace_Hn, rcs_teta);

	ContribRcs(3*num_angle) = rcs_teta(0);
	ContribRcs(3*num_angle+1) = rcs_teta(1);
	ContribRcs(3*num_angle+2) = rcs_teta(2);
      }
    
    // summing value between processors
    VectComplex_wp AddRcs(ContribRcs);
#ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    MpiReduce(var_maxwell.comm_group_mode, ContribRcs, xtmp, AddRcs,
	      AddRcs.GetM(), MPI_SUM, 0);
#endif
    
    if (rank_proc == 0)
      {
	// updates RCS_Vector
	for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
	  {
	    rcs_teta(0) = AddRcs(3*num_angle);
	    rcs_teta(1) = AddRcs(3*num_angle+1);
	    rcs_teta(2) = AddRcs(3*num_angle+2);
	    RCS_Vector(num_angle) += rcs_teta;
	  }
	
	// writes in a file
	ofstream file_out(this->file_RCS.data()); file_out.precision(15);
	ofstream file_far;
	bool write_far_field = (this->file_far_field.size() > 0);
	if (write_far_field)
	  {
	    file_far.open(this->file_far_field.data());
	    file_far.precision(15);
	  }
	
	for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
	  {
	    teta = this->first_angle_RCS + step_angle*num_angle;
	    rcs_teta = RCS_Vector(num_angle);
	    
	    Real_wp module_h = Norm2(rcs_teta); module_h *= module_h; // DISP(module_h); 
	    Real_wp sigma = 10.0*log(1.0/(4.0*pi_wp)*module_h)/log(10.0);
	
	    file_out << (180.0*teta/pi_wp) << "  " << sigma << '\n';
	    
	    if (write_far_field)
	      {
		file_far << (180.0*teta/pi_wp) << "  "
			 << real(rcs_teta(0)) << " " << imag(rcs_teta(0)) << " "
			 << real(rcs_teta(1)) << " " << imag(rcs_teta(1)) << " "
			 << real(rcs_teta(2)) << " " << imag(rcs_teta(2)) << '\n';
	      }
	  }
	
	file_out.close();
	if (write_far_field)
	  file_far.close();
      }
    
    if ((rank_proc == 0)||(this->print_level >= 10))
      if (this->print_level >= 0)
	cout<<"Radar cross section OK"<<endl;
  }


  //! computes E and H on quadrature points of Gamma
  void VarComputationRCS_MaxwellAxi
  ::ComputeEnHn(const VectComplex_wp& U0, VectR3_Complex_wp& trace_En,
		VectR3_Complex_wp& trace_Hn)
  {
    Vector<VectComplex_wp> U0_vec(3);
    int N1 = var_maxwell.offset_dof_unknown(1);
    int N2 = var_maxwell.offset_dof_unknown(2) - N1;
    U0_vec(0).SetData(N1, const_cast<Complex_wp*>(&U0(0)));
    U0_vec(1).SetData(N2, const_cast<Complex_wp*>(&U0(N1)));
    if (var_maxwell.FirstOrderFormulation())
      U0_vec(2).SetData(U0.GetM()-(N1+N2), const_cast<Complex_wp*>(&U0(N1+N2)));
    
    // first we evaluate values of H on nodal points
    Vector<VectComplex_wp> EvalH;
    var_maxwell.EvaluateH_MixedFormulation(U0_vec, EvalH);
    
    // then E and H are evaluated on quadrature points of the boundary
    trace_En.Reallocate(this->var_mesh.GetNbAllQuadraturePoints());
    trace_Hn.Reallocate(this->var_mesh.GetNbAllQuadraturePoints());
    VectReal_wp points_surf = this->var_mesh.PointsReferenceSurface()(0);
    int nb = 0; R2 pt_glob, pt_loc; 
    Vector<TinyVector<Complex_wp, 3> > val_E(points_surf.GetM()), val_H(points_surf.GetM());
    for (int i = 0; i < this->var_mesh.GetNbBoundary(); i++)
      {
        int num_elem = this->var_mesh.GetElementNumberOfSurface(i);
        int num_loc = this->var_mesh.GetLocalPositionOfSurface(i);        

	const ElementReference<Dimension2, 2>& Fb_hcurl = var_maxwell.GetReferenceElementHcurl(num_elem, 0);
        const ElementReference<Dimension2, 1>& Fb_h1 = var_maxwell.GetReferenceElementH1(num_elem, 1);
	for (int j = 0; j < this->var_mesh.GetNbPointsQuadrature(i); j++)
	  {
	    pt_glob = this->var_mesh.PointsQuadrature(i, j);                  
	    Fb_h1.GetLocalCoordOnBoundary(num_loc, points_surf(j), pt_loc);
	    var_maxwell.ComputeEpolar(U0_vec, num_elem, pt_loc, pt_glob,
				      this->var_mesh.Dfjm1Quadrature(i, j), val_E(j), Fb_hcurl, Fb_h1);
	    
	    var_maxwell.ComputeHpolar(EvalH, num_elem, pt_loc, pt_glob,
				      this->var_mesh.Dfjm1Quadrature(i, j), val_H(j), Fb_hcurl, Fb_h1);
	  }

        //int ig = var_maxwell.mesh.Element(num_elem).numBoundary(num_loc);
        for (int j = 0; j < val_E.GetM(); j++)
          {
            trace_En(nb) = val_E(j);
            trace_Hn(nb) = val_H(j);
            nb++;
          }
      }

    U0_vec(0).Nullify();
    U0_vec(1).Nullify();
    if (var_maxwell.FirstOrderFormulation())
      U0_vec(2).Nullify();
    
    /* ofstream file_pt("Points.dat"), file_E, file_H, file_w("Weights.dat");
    file_pt.precision(15); file_E.precision(15); file_H.precision(15); file_w.precision(15);
    if (TypeEquation::number_mode == 1)
      {
	file_E.open("ValE1.dat"); file_H.open("ValH1.dat");
      }
    else
      {
	file_E.open("ValEm1.dat"); file_H.open("ValHm1.dat");
      }
    
    for (int i = 0; i < this->var_mesh.GetNbAllQuadraturePoints(); i++)
      {
        file_pt << this->var_mesh.GetQuadraturePoint(i)(0) 
                << " 0 " << this->var_mesh.GetQuadraturePoint(i)(1) << '\n';
	
	file_w << this->var_mesh.GetQuadratureWeight(i) << '\n';
        
        file_E << real(trace_En(i)(0)) << " " << imag(trace_En(i)(0)) << " "
               << real(trace_En(i)(1)) << " " << imag(trace_En(i)(1)) << " "
               << real(trace_En(i)(2)) << " " << imag(trace_En(i)(2)) << '\n';

        file_H << real(trace_Hn(i)(0)) << " " << imag(trace_Hn(i)(0)) << " "
               << real(trace_Hn(i)(1)) << " " << imag(trace_Hn(i)(1)) << " "
               << real(trace_Hn(i)(2)) << " " << imag(trace_Hn(i)(2)) << '\n';
      }
    
    file_pt.close(); file_E.close(); file_H.close(); file_w.close(); */
    
    //VectR3_Complex_wp trace_En, trace_Hn;
    //this->var_mesh.ComputeEnHnOnBoundary(var, U0, trace_En, trace_Hn);
  }
  
  
  //! adds contribution of the current mode to the far field
  void VarComputationRCS_MaxwellAxi
  ::AddContributionMode(int num_angle, const R3& ur, const VectR3_Complex_wp& trace_En,
			const VectR3_Complex_wp& trace_Hn, R3_Complex_wp& rcs_teta)
  {
    R3_Complex_wp Hn_u; 
    
    // integration over all quadrature points of the boundary
    R3_Complex_wp rcs_tmp; rcs_tmp.Fill(Complex_wp(0,0));
    Real_wp z, kortho, kz, teta0;
    
    // decomposition of exp(i k u \cdot x)
    // \, = \, exp( i k_\bot \, r cos(teta - teta_0)) * exp(i kz \, z) 
    kortho = sqrt(ur(0)*ur(0) + ur(1)*ur(1));
    if (ur(0) >= kortho)
      teta0 = 0.0;
    else if (ur(0) <= -kortho)
      teta0 = pi_wp;
    else
      teta0 = acos(ur(0)/kortho);
    
    if (ur(1) < 0)
      teta0 = -teta0;
	
    kortho *= var_maxwell.GetOmega();
    kz = var_maxwell.GetOmega()*ur(2);
    
    VectReal_wp Jn(3); Complex_wp exp_kz, Er, Eteta, Ez, Etimes_n, Hr, Hteta, Hz, Htimes_n;
    Complex_wp Im1, Im, Ip1, exp_m, exp_m1, exp_p1, int_cos, int_sin, int_one; 
    R3_Complex_wp Int_En, Int_Hn; Real_wp nx, nz, r;
    
    int number_mode = var_maxwell.GetCurrentModeNumber();
    Im1 = ComputePowerI(number_mode-1);
    Im = Im1*Iwp; Ip1 = Im*Iwp;
    exp_m1 = exp(-Complex_wp(number_mode-1)*Iwp*teta0);
    exp_m = exp(-Complex_wp(number_mode)*Iwp*teta0);
    exp_p1 = exp(-Complex_wp(number_mode+1)*Iwp*teta0);
    
    for (int ind = 0; ind < this->var_mesh.GetNbAllQuadraturePoints(); ind++)
      {
	r = this->var_mesh.GetQuadraturePoint(ind)(0);
	z = this->var_mesh.GetQuadraturePoint(ind)(1);
	//kr = kortho*r; 
	exp_kz = exp(Iwp*kz*z);
	int num_point = num_angle*this->var_mesh.GetNbAllQuadraturePoints()+ind;
	// Jn(0) = J_{m-1}(k_\bot r), Jn(1) = J_{m}(k_\bot r), Jn(2) = J_{m+1}(k_\bot r)
	Jn(0) = var_maxwell.GetBessel_Value(number_mode-1, num_point);
	Jn(1) = var_maxwell.GetBessel_Value(number_mode, num_point);
	Jn(2) = var_maxwell.GetBessel_Value(number_mode+1, num_point);
	
	// int_cos : integral of   cos(theta) exp(-i m theta) exp(i k u \cdot x)
	// int_sin : integral with weight sin(theta) and 1
	// use of Jacobi-Anger expansion exp(i k_\bot r cos(theta)) = \sum i^n J_n(k_\bot r) e^{i n theta}
	int_cos = pi_wp*(exp_m1*Jn(0)*Im1 + exp_p1*Jn(2)*Ip1)*exp_kz;
	int_sin = -pi_wp*Iwp*(exp_m1*Jn(0)*Im1 - exp_p1*Jn(2)*Ip1)*exp_kz;
	int_one = 2.0*pi_wp*exp_m*Jn(1)*Im*exp_kz;
	    
	// value of E and H in cylindrical coordinates
	Er = trace_En(ind)(0); Hr = trace_Hn(ind)(0);
	Eteta = trace_En(ind)(1); Hteta = trace_Hn(ind)(1);
	Ez = trace_En(ind)(2); Hz = trace_Hn(ind)(2);
	    
	// normale (nx, 0, nz)
	nx = this->var_mesh.GetQuadratureNormale(ind)(0);
	nz = this->var_mesh.GetQuadratureNormale(ind)(1);
	
	// Int_En contains the integral in theta of 
	// E \times n exp(-i m theta) exp(i k u \cdot x)
	Etimes_n = (Er*nz - Ez*nx);
	Int_En(0) = Etimes_n*int_sin + Eteta*nz*int_cos;
	Int_En(1) = -Etimes_n*int_cos + Eteta*nz*int_sin;
	Int_En(2) = -Eteta*nx*int_one;
	
	// Int_Hn contains the integral in theta of 
	// H \times n exp(-i m theta) exp(i k u \cdot x)
	Htimes_n = (Hr*nz - Hz*nx);
	Int_Hn(0) = Htimes_n*int_sin + Hteta*nz*int_cos;
	Int_Hn(1) = -Htimes_n*int_cos + Hteta*nz*int_sin;
	Int_Hn(2) = -Hteta*nx*int_one;
	
	// Hn_u will contain the integral in theta
	// of (E \times n) + u \times (H \times n) exp(i k u \cdot x)
	TimesProd(ur, Int_Hn, Hn_u);
	Hn_u += Int_En;
	
	// then performing the integration over the 2-D curve (r, z):
	Add(this->var_mesh.GetQuadratureWeight(ind)*r, Hn_u, rcs_tmp);
      }
    
    TimesProd(ur, rcs_tmp, Hn_u);
    Add(Iwp*var_maxwell.GetOmega(), Hn_u, rcs_teta);
  }
  
  
  //! updating monostatic RCS
  /*!
    \param[in] U0 solution vector for given mode
    \param[in,out] rcs_teta RCS vector to update
    \param[in] time_process timers
    \param[in] num_angle angle number
   */
  void VarComputationRCS_MaxwellAxi
  ::ComputeMonostaticRCS(const Matrix<Complex_wp, General, ColMajor>& U0, VectR3_Complex_wp& RCS_Vector)
  {        
#ifdef SELDON_WITH_MPI   
    int rank_proc; MPI_Comm_rank(var_maxwell.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    VectR3_Complex_wp trace_En, trace_Hn;

    R3_Complex_wp rcs_teta; VectComplex_wp ContribRcs(3*U0.GetN());
    VectComplex_wp U0_vec;
    for (int k = 0; k < U0.GetN(); k++)
      {
	rcs_teta.Fill(0);
	
	// computation of E and H on quadrature points
	U0_vec.SetData(U0.GetM(), const_cast<Complex_wp*>(&U0(0, k)));
	ComputeEnHn(U0_vec, trace_En, trace_Hn);
	U0_vec.Nullify();
	
	// then computation of rcs in the same direction than wave vector
	AddContributionMode(k/2, this->incident_direction(k/2), trace_En, trace_Hn, rcs_teta);	
	
	ContribRcs(3*k) = rcs_teta(0);
	ContribRcs(3*k+1) = rcs_teta(1);
	ContribRcs(3*k+2) = rcs_teta(2);
      }

    // summing value between processors
    VectComplex_wp AddRcs(ContribRcs);
#ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    MpiReduce(var_maxwell.comm_group_mode, ContribRcs, xtmp, AddRcs,
	      AddRcs.GetM(), MPI_SUM, 0);
#endif
    
    if (rank_proc == 0)
      {
	// updates RCS_Vector
	for (int k = 0; k < U0.GetN(); k++)
	  {
	    rcs_teta(0) = AddRcs(3*k);
	    rcs_teta(1) = AddRcs(3*k+1);
	    rcs_teta(2) = AddRcs(3*k+2);
	    RCS_Vector(k) += rcs_teta;
	  }
      }
  }
  
  
  //! writing radar cross section
  /*!
    \param[in] RCS_Vector 3-D far field vector
     The radar cross section is written on output files 
   */
  void VarComputationRCS_MaxwellAxi
  ::WriteMonostaticRCS(const VectR3_Complex_wp& RCS_Vector)
  {
    Real_wp rcs_hh, rcs_vv, rcs_hv, rcs_vh;
    Real_wp step_angle = (this->last_angle_RCS-this->first_angle_RCS) / (this->nb_angles_RCS-1);
    Real_wp teta;
    
    ofstream file_out(this->file_RCS.data()); file_out.precision(15);
    
    ofstream file_far;
    bool write_far_field = (this->file_far_field.size() > 0);
    if (write_far_field)
      {
	file_far.open(this->file_far_field.data());
	file_far.precision(15);
      }
    
    // loop on each angle
    Real_wp constante = 1.0/(4.0*pi_wp);
    for (int num_angle = 0; num_angle < this->nb_angles_RCS; num_angle++)
      {
	teta = this->first_angle_RCS + step_angle*num_angle;
	rcs_hh = absSquare(RCS_Vector(2*num_angle)(0)) + absSquare(RCS_Vector(2*num_angle)(2));
	rcs_hv = absSquare(RCS_Vector(2*num_angle)(1));

	rcs_vv = absSquare(RCS_Vector(2*num_angle+1)(1));
	rcs_vh = absSquare(RCS_Vector(2*num_angle+1)(0)) + absSquare(RCS_Vector(2*num_angle+1)(2));
	
	rcs_hh = 10.0*log(constante*rcs_hh)/log(10.0);
	rcs_vv = 10.0*log(constante*rcs_vv)/log(10.0); 
	rcs_hv = 10.0*log(constante*rcs_hv)/log(10.0);
	rcs_vh = 10.0*log(constante*rcs_vh)/log(10.0); 
	
	file_out << (180.0*teta/pi_wp) << " " << rcs_hh << " " << rcs_vv << " "
		 << rcs_hv << " " << rcs_vh << '\n';
	
	if (write_far_field)
	  {
	    file_far << (180.0*teta/pi_wp) << " " 
		     << real(RCS_Vector(2*num_angle)(0)) << " " << imag(RCS_Vector(2*num_angle)(0)) << " "
		     << real(RCS_Vector(2*num_angle)(1)) << " " << imag(RCS_Vector(2*num_angle)(1)) << " "
		     << real(RCS_Vector(2*num_angle)(2)) << " " << imag(RCS_Vector(2*num_angle)(2)) << " "
		     << real(RCS_Vector(2*num_angle+1)(0)) << " " << imag(RCS_Vector(2*num_angle+1)(0)) << " "
		     << real(RCS_Vector(2*num_angle+1)(1)) << " " << imag(RCS_Vector(2*num_angle+1)(1)) << " "
		     << real(RCS_Vector(2*num_angle+1)(2)) << " " << imag(RCS_Vector(2*num_angle+1)(2)) << '\n';
	  }
      }
    
    file_out.close();
    if (write_far_field)
      file_far.close();
  }    
  
}

#define MONTJOIE_FILE_MAXWELL_AXISYM_HARMONIC_CXX
#endif


