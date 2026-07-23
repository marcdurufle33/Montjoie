#ifndef MONTJOIE_FILE_HELMHOLTZ_HDIV_CXX

#include "HelmholtzHdivInline.cxx"

namespace Montjoie
{

  /*************************
   * HelmholtzEquationHdiv *
   *************************/

  
  //! computation of geometric quantities (omega/J_i and 1/J_i DF^* DF)
  template<class T, class Dimension>
  template<class TypeEquation>
  void HelmholtzEquationHdiv_Base<T, Dimension>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		    int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
    T momega2, m_iomega;
    TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> Hess_dfj;
    typename Dimension::MatrixN_N dfjm1, dfj, dfj_trans;    
    //int ref = var.mesh.Element(num_elem).GetReference();
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    
    // if constant DFi and constant physical constants, we compute
    // geometric quantities only for a quadrature point
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;
    
    Real_wp jacob, poids;    
    var.GetMiomega(m_iomega);
    var.GetMomega2(momega2);

    var.Glob_matMass_Dh(num_elem).Reallocate(N);
    var.Glob_matMass_Bh(num_elem).Reallocate(N);
    var.Glob_matMass_Ah(num_elem).Reallocate(N);
    var.Glob_matMass_AhSigma(num_elem).Reallocate(N);
    var.Glob_matMass_Vh(num_elem).Reallocate(N);
    
    // loop over quadrature points
    Real_wp rhoC2, rho, sigma;
    TinyVector<Real_wp, Dimension::dim_N> grad_P, vec_C;
    TinyMatrix<Real_wp, Symmetric, Dimension::dim_N, Dimension::dim_N> HessP;
    TinyMatrix<Real_wp, Symmetric, Dimension::dim_N, Dimension::dim_N> mat_mass, mat_mass_sigma;
    for (int j = 0; j < N; j++)
      {
        rhoC2 = var.eval_rhoC2(num_elem)(j);
	grad_P = var.grad_p0(num_elem)(j);
	HessP = var.hessian_p0(num_elem)(j);
        rho = var.eval_rho(num_elem)(j);
        sigma = var.eval_sigma(num_elem)(j);
        
	if (variable)
	  {
	    poids = Fb.WeightsND(j);
	    if (affine)
	      {
		jacob = var.Glob_jacobian(num_elem)(0);
		dfjm1 = var.Glob_DFjm1(num_elem)(0);
	      }
	    else
	      {
		jacob = var.Glob_jacobian(num_elem)(j) / Fb.WeightsND(j);
		dfjm1 = var.Glob_DFjm1(num_elem)(j);
	      }
	  }
	else
	  {
	    poids = 1.0;
	    jacob = var.Glob_jacobian(num_elem)(0);
	    dfjm1 = var.Glob_DFjm1(num_elem)(0);
	  }       
	
	Mlt(1.0/jacob, dfjm1);
	GetInverse(dfjm1, dfj);
	Transpose(dfj, dfj_trans);
	
	// geometrical term for -omega^2 rho u . v => -omega^2 (1/Ji DF_i^T rho DF_i) 
	Mlt(dfj_trans, dfj, mat_mass);
	Mlt(poids/jacob, mat_mass);
	var.Glob_matMass_Ah(num_elem)(j) = rho*momega2*mat_mass;
	
        // sigma term
        var.Glob_matMass_AhSigma(num_elem)(j) = 2.0*sigma*rho*m_iomega*mat_mass;
        mat_mass_sigma = mat_mass;
        
        // term with hessian of pressure and sigma^2 term
	Mlt(HessP, dfj, Hess_dfj);
	Mlt(dfj_trans, Hess_dfj, mat_mass);
	Mlt(poids/jacob, mat_mass);
        Add(sigma*sigma*rho, mat_mass_sigma, mat_mass);
	var.Glob_matMass_Bh(num_elem)(j) = mat_mass;
        
	// geometrical term for C . u  div v => 1/Ji DF_i^T C
	Mlt(dfj_trans, grad_P, vec_C);
	Mlt(poids/jacob, vec_C);
	var.Glob_matMass_Vh(num_elem)(j) = vec_C;
        
	// geometrical term for rho c^2 div u div v => (rho c^2)/Ji 
	var.Glob_matMass_Dh(num_elem)(j) = poids/jacob*rhoC2;
      }
  }

  
  /*************************
   * VarHelmholtzHdiv_Base *
   *************************/


  //! allocation of arrays containing physical indices
  template<class Dimension>
  void VarHelmholtzHdiv_Base<Dimension>::InitIndices(int n)
  {
    ref_rho0.Reallocate(n);
    ref_p0.Reallocate(n);
    ref_grad_p0.Reallocate(n);
    ref_c0.Reallocate(n);
    ref_sigma.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
        ref_sigma(i).SetConstant(0.0);
	ref_rho0(i).SetConstant(Real_wp(1));
	ref_c0(i).SetConstant(Real_wp(1));
	//ref_p0(i).SetConstant(Real_wp(1));
	typename Dimension::R_N zero;
        ref_grad_p0(i).SetConstant(zero);
      }
  }
  
  
  //! modifies physical indices with a line of the data file
  template<class Dimension>
  void VarHelmholtzHdiv_Base<Dimension>::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho0.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho0.GetM() << endl;
        abort();
      }

    ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    ref_rho0(i).SetInputData(nb, parameters, parameters(0));
    ref_c0(i).SetInputData(nb, parameters, parameters(0));
    //ref_p0(i).SetInputData(nb, parameters, parameters(0));
    ref_grad_p0(i).SetInputData(nb, parameters, parameters(0));
  }
  
  
  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Dimension>
  void VarHelmholtzHdiv_Base<Dimension>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho0.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho0.GetM() << endl;
        abort();
      }

    if (name_media == "sigma")
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "rho0")
      ref_rho0(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "c0")
      ref_c0(i).SetInputData(nb, parameters, parameters(0));
    //ref_p0(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "grad_p0")
      ref_grad_p0(i).SetInputData(nb, parameters, parameters(0));
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }


  //! returns the name associated with the physical index num
  template<class Dimension>
  string VarHelmholtzHdiv_Base<Dimension>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("sigma");
      case 1: return string("rho0");
      case 2: return string("c0");
      case 3: return string("grad_p0");
      }

    return string();
  }  
  
    
  //! fills varying indices present in the current problem
  template<class Dimension>
  void VarHelmholtzHdiv_Base<Dimension>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&  rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    int nb = 0;
    for (int i = 0; i < ref_rho0.GetM(); i++)
      {
        nb += ref_sigma(i).GetNbVaryingMedia();
        nb += ref_rho0(i).GetNbVaryingMedia();
        nb += ref_c0(i).GetNbVaryingMedia();
        //nb += ref_p0(i).GetNbVaryingMedia();
        nb += ref_grad_p0(i).GetNbVaryingMedia();
      }
    
    rho_real.Reallocate(nb);
    num_ref.Reallocate(nb);
    num_index.Reallocate(nb);
    num_component.Reallocate(nb);
    compute_grad.Reallocate(nb);
    compute_hess.Reallocate(nb);
    compute_grad.Fill(false);
    compute_hess.Fill(false);
    nb = 0;
    for (int i = 0; i < ref_rho0.GetM(); i++)
      {
        int nb0 = nb;
        ref_sigma(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }

	nb0 = nb;
        ref_rho0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 1;
            num_ref(j) = i;
          }

	nb0 = nb;
        ref_c0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 2;
            num_ref(j) = i;
          }

	nb0 = nb;
        //ref_p0(i).GetVaryingMedia(nb, rho_real, num_component);
        ref_grad_p0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 3;
            num_ref(j) = i;
            compute_grad(j) = true;
          }
      }
  }
  
  
  template<class Dimension>
  void VarHelmholtzHdiv_Base<Dimension>::ComputePhysicalCoefficients()
  {
    int nb_elt = var_problem.mesh.GetNbElt();
    eval_rho.Reallocate(nb_elt);
    eval_rhoC2.Reallocate(nb_elt);
    eval_sigma.Reallocate(nb_elt);
    
    grad_p0.Reallocate(nb_elt);
    hessian_p0.Reallocate(nb_elt);
    
    Real_wp c;
    TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> hessP;
    for (int i = 0; i < nb_elt; i++)
      {
        int ref = var_problem.mesh.Element(i).GetReference();
        int N = var_problem.Glob_PointsQuadrature(i).GetM();
        N = var_problem.GetNbPointsQuadratureInside(i);
        eval_sigma(i).Reallocate(N);
	eval_rho(i).Reallocate(N);
	eval_rhoC2(i).Reallocate(N);
        
        grad_p0(i).Reallocate(N);
        hessian_p0(i).Reallocate(N);
        
        for (int j = 0; j < N; j++)
          {
            eval_rho(i)(j) = this->ref_rho0(ref).GetCoefficient(var_problem, i, j);
            eval_sigma(i)(j) = this->ref_sigma(ref).GetCoefficient(var_problem, i, j);
            //this->ref_p0(ref).GetCoefGradient(var_problem, i, j, p0, grad_p0(i)(j));
            this->ref_grad_p0(ref).GetCoefGradient(var_problem, i, j, grad_p0(i)(j), hessP);
            hessian_p0(i)(j) = hessP;
	    c = this->ref_c0(ref).GetCoefficient(var_problem, i, j);
	    eval_rhoC2(i)(j) = eval_rho(i)(j)*c*c;
          }
      }
    
    for (int ref = 0; ref < ref_rho0.GetM(); ref++)
      {
	//this->ref_p0(ref).Clear();
        this->ref_grad_p0(ref).Clear();
	this->ref_c0(ref).Clear();
	this->ref_rho0(ref).Clear();
        this->ref_sigma(ref).Clear();
      }
  }


  /***********************
   * VarHelmholtzHdiv_Eq *
   ***********************/
  
  
  template<class TypeEquation>
  VarHelmholtzHdiv_Eq<TypeEquation>::VarHelmholtzHdiv_Eq()
    : VarHelmholtzHdiv_Base<Dimension>(this->GetLeafClass())
  {
  }


  template<class TypeEquation>
  void VarHelmholtzHdiv_Eq<TypeEquation>::InitIndices(int n)
  {
    VarHelmholtzHdiv_Base<Dimension>::InitIndices(n);
  }
  
  
  template<class TypeEquation>
  void VarHelmholtzHdiv_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    VarHelmholtzHdiv_Base<Dimension>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  void VarHelmholtzHdiv_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarHelmholtzHdiv_Base<Dimension>::SetPhysicalIndex(name_media, i, parameters);
  }
  

  template<class TypeEquation>
  string VarHelmholtzHdiv_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarHelmholtzHdiv_Base<Dimension>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  int VarHelmholtzHdiv_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarHelmholtzHdiv_Base<Dimension>::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  bool VarHelmholtzHdiv_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarHelmholtzHdiv_Base<Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool VarHelmholtzHdiv_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarHelmholtzHdiv_Base<Dimension>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  Real_wp VarHelmholtzHdiv_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return VarHelmholtzHdiv_Base<Dimension>::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  Real_wp VarHelmholtzHdiv_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return VarHelmholtzHdiv_Base<Dimension>::GetVelocityOfInfinity();
  }
    
  
  template<class TypeEquation>
  void VarHelmholtzHdiv_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >&  rho_cplx,
		      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >&  rho_real, IVect& num_ref,
 		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    VarHelmholtzHdiv_Base<Dimension>::GetVaryingIndices(rho_real, num_ref, num_index,
							num_component, compute_grad, compute_hess);
  }
  
  
  template<class TypeEquation>
  void VarHelmholtzHdiv_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();
    VarHelmholtzHdiv_Base<Dimension>::ComputePhysicalCoefficients();
  }

  template<class TypeEquation>
  int VarHelmholtzHdiv_Eq<TypeEquation>
  ::GetBoundaryConditionId(const IVect& ref, int pos, const VectString& parameters, bool& periodic)
  {
    if (parameters(pos) == "ATMOSPHERE")
      {
	periodic = false;
	if (this->Param_condition(ref(0)).GetM() <= 1)
	  {
	    cout << "You need to give the order of the atmosphere boundary condition and alpha parameter" << endl;
	    cout << "Current parameters are " << this->Param_condition(ref(0)) << endl;
	    abort();
	  }
		
	int order = toInteger(realpart(this->Param_condition(ref(0))(0)));
	this->order_ABC = 20 + order;
	this->gamma_cla_coef = 1.0/realpart(this->Param_condition(ref(0))(1));
	this->take_into_account_curvature_for_abc = true;
	return BoundaryConditionEnum::LINE_ABSORBING;
      }
    
    return VarHarmonic<TypeEquation>::GetBoundaryConditionId(ref, pos, parameters, periodic);
  }


  //! computation of elementary matrices for Helmholtz equation and H(div) formulation
  /*!
    Considered equation is equal to -\omega^2 v - grad( div v) = 0
    v is a vector that can be approximation with H(div) elements 
    Associated variational formulation is equal to
    -\omega^2 \int_\Omega u v + \int_\Omega div u div v
   */
  template<class TypeEquation> template<class T>
  void VarHelmholtzHdiv_Eq<TypeEquation>
  ::ComputeElementaryMatrixHelm(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
				const GlobalGenericMatrix<T>& nat_mat,
				const ElementReference<Dimension, 3>& Fb)
  {
    int nb_dof_elt = Fb.GetNbDof();
    int nb_dof_all = nb_dof_elt;
    
    const MeshNumbering<Dimension>& mesh_num = this->mesh_num;
        
    // number of integration points
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    //int ref_domain = vars.mesh.Element(iquad).GetReference();
    
    // dof numbers
    num_dof.Reallocate(nb_dof_all);
    IVect Nodle = this->GetDofNumberOnElement(iquad);
    for (int j = 0; j < nb_dof_elt; j++)
      num_dof(j) = Nodle(j);
    
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();

    bool variable = this->UseNumericalIntegration(iquad);
    //bool affine = vars.mesh.IsElementAffine(iquad);

    int N = nb_points_quad;
    if (!variable)
      N = 1;
    
    Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > MatMass(N);
    Vector<TinyVector<T, Dimension::dim_N> > VecGrad(N);
    Vector<T> CoefStiff(N);
    for (int j = 0; j < N; j++)
      {
        MatMass(j) = this->Glob_matMass_Bh(iquad)(j);
        Mlt(nat_mat.GetCoefStiffness(), MatMass(j));
        Add(nat_mat.GetCoefDamping(), this->Glob_matMass_AhSigma(iquad)(j), MatMass(j));
        Add(nat_mat.GetCoefMass(), this->Glob_matMass_Ah(iquad)(j), MatMass(j));
        
        VecGrad(j) = this->Glob_matMass_Vh(iquad)(j);
        Mlt(nat_mat.GetCoefStiffness(), VecGrad(j));
        
        CoefStiff(j) = this->Glob_matMass_Dh(iquad)(j)*nat_mat.GetCoefStiffness();
      }
    
    VectReal_wp Ones(nb_dof_elt);
    Vector<Real_wp> val_phi(Dimension::dim_N*nb_points_quad);
    Vector<Real_wp> div_phi(nb_points_quad);
    Vector<T> feval_phi(Dimension::dim_N*nb_points_quad);
    Vector<T> feval_div(nb_points_quad);
    Vector<T> contrib(nb_dof_elt), contrib_div(nb_dof_elt);
    val_phi.Fill(0); div_phi.Fill(0);  feval_phi.Fill(0); feval_div.Fill(0);
    contrib.Fill(0); contrib_div.Fill(0);
    TinyVector<Real_wp, Dimension::dim_N> vec_u;
    TinyVector<T, Dimension::dim_N> vec_v;
    // loop over columns
    for (int i = 0; i < nb_dof_elt; i++)
      {
	Ones.Fill(0); Ones(i) = 1.0;
	
	Fb.ApplyChTranspose(Ones, val_phi);
	if (Fb.UseQuadraturePointsForRh())
	  Fb.ApplyRhQuadratureTranspose(val_phi, div_phi);
	else
	  Fb.ApplyRhTranspose(Ones, div_phi);
	
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    CopyVector(val_phi, j, vec_u); 
	    if (variable)
              {
                Mlt(MatMass(j), vec_u, vec_v);
                Add(div_phi(j), VecGrad(j), vec_v);
              }
	    else
	      {
		Mlt(MatMass(0), vec_u, vec_v);
                Add(div_phi(j), VecGrad(0), vec_v);
		vec_v *= Fb.WeightsND(j);
	      }
	    
	    CopyVector(vec_v, j, feval_phi);
	    
	    if (variable)
	      {
                feval_div(j) = div_phi(j)*CoefStiff(j);
                feval_div(j) += DotProd(VecGrad(j), vec_u);
              }
	    else
	      {
                feval_div(j) = div_phi(j)*CoefStiff(0)*Fb.WeightsND(j);
                feval_div(j) += DotProd(VecGrad(0), vec_u)*Fb.WeightsND(j);
              }
	  }
	
	Fb.ApplyRh(feval_div, contrib_div);
	Fb.ApplyCh(feval_phi, contrib);
	
	for (int j = 0; j < nb_dof_elt; j++)
	  mat_interac.SetEntry(j, i, contrib(j) + contrib_div(j));
      }
    
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad, 1);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad, 1);
  }
  
  
  template<class TypeEquation>
  void VarHelmholtzHdiv_Eq<TypeEquation>::AllocateMassMatrices()
  {
    this->Glob_matMass_Dh.Reallocate(this->mesh.GetNbElt());
    this->Glob_matMass_Ah.Reallocate(this->mesh.GetNbElt());
    this->Glob_matMass_Bh.Reallocate(this->mesh.GetNbElt());
    this->Glob_matMass_AhSigma.Reallocate(this->mesh.GetNbElt());
    this->Glob_matMass_Vh.Reallocate(this->mesh.GetNbElt());
  }
  
  
  /*********************************
   * DiffractedWaveSource_HelmHdiv *
   *********************************/
  
  
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>::InitDefaultValues()
  {
    coef_vol = 1.0; invRho0 = 1.0; invRho = 1.0;
    
    // variables for gradient source
    coef_grad.SetIdentity(); invMu0 = 1.0;
    invMu.SetIdentity();
    
    omega2 = var_problem.GetSquareOmega();
    coef_invMu0 = 0.0; coef_invMu = 1.0;
    
    varying_rho = false; varying_mu = false;
  }
  
  
  //! initialization for each element
  /*!
    \param[in] num_elem element number
    if physical properties are constant inside the element
    coefficients \f$ -\omega^2(\rho - \rho_0) \f$ and (mu-mu0) are computed
   */
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>
  ::InitElement(int num_elem, const VectR_N& s)
  {
  }

  
  //! volumetric source if rho <> rho0
  template<class Dimension>
  bool DiffractedWaveSource_HelmHdiv<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    // if rho different from rho 0, we have a volumic source \int f \phi
    if ((coef_grad.IsZero())&&(!varying_mu))
      return false;
    
    return true;
  }

  
  //! evaluation of incident wave
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f)
  {        
    // for Dirichlet condition u = -u_inc
    if (this->num_loc_ >= 0)
      {
	Complex_wp u_inc; R_N_Complex_wp grad_uinc;
	this->incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);

        if (var_problem.InsidePML(i))
	  grad_uinc.Zero();
	else
	  grad_uinc *= -1.0;
	
	CopyVector(grad_uinc, 0, f);
      }
  }

  
  //! evaluation of volumetric source
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
  }
  

  //! volumetric source for \f$ \nabla phi \f$ if mu <> mu0
  template<class Dimension>
  bool DiffractedWaveSource_HelmHdiv<Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    // gradient source if mu different from mu_0
    if ((coef_vol != Complex_wp(0))||varying_rho)
      return true;
    
    return false;
  }
  
  //! evaluation of volumetric source g in \f$ \int g \nabla \varphi \f$
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectComplex_wp& f)
  {
  }
  
  //! initialization before evaluation of surfacic source
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Complex_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);
  }

  
  //! surfacic for Neumann and High-Conductivity boundary conditions
  template<class Dimension>
  bool DiffractedWaveSource_HelmHdiv<Dimension>::
  IsNonNullSurfacicSource(int ref)
  {    
    // for neumann condition and high conductivity condition
    int cond_ref = var_problem.mesh.GetBoundaryCondition(ref);
    reference_condition = cond_ref;
    if (cond_ref == BoundaryConditionEnum::LINE_NEUMANN)
      return true;
    
    return false;
  }

  
  //! evaluation of surfacic source
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>::
  EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                         const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    R_N x = PointsElem.GetPointQuadratureBoundary(k);
    Complex_wp u_inc;
    incident_wave.EvaluateFunction(x, u_inc);
    
    R_N_Complex_wp fvec;
    if (reference_condition == BoundaryConditionEnum::LINE_NEUMANN)
      {
	// neumann condition -\int \mu_0 du_inc/dn \phi 	
	fvec = MatricesElem.GetNormaleQuadratureBoundary(k);
        fvec *= -u_inc;
      }
    
    CopyVector(fvec, 0, f);
  }
  

  //! surfacic source \f$ \int f \frac{d\varphi}{ds} \f$  for high-conductivity boundary condition
  template<class Dimension>
  bool DiffractedWaveSource_HelmHdiv<Dimension>::
  IsNonNullSurfacicSourceGradient(int ref)
  {
    return false;
  }

  
  //! evaluation of surfacic source f in \f$ \int f \frac{d \varphi}{ds} \f$ 
  template<class Dimension>
  void DiffractedWaveSource_HelmHdiv<Dimension>::
  EvaluateSurfacicSourceGradient(int k, const SetPoints<Dimension>& PointsElem,
				  const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
  }
  
  
  template<class Dimension>
  void IncidentWaveProjector_HelmHdiv<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectComplex_wp& f)
  {
    Complex_wp u_inc; typename Dimension::R_N_Complex_wp grad_uinc;
    this->incident_wave.EvaluateFunctionGradient(x, u_inc, grad_uinc);
    
    grad_uinc *= -1.0/var_problem.GetSquareOmega();
    CopyVector(grad_uinc, 0, f);
  }
  

  /****************
   * ImpedanceABC *
   ****************/
  

  //! computation of impedance
  template<class Dimension>
  void ImpedanceABC<Complex_wp, HelmholtzEquationHdiv<Dimension> >::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
		       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
		       const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    // storing the normale
    if (k == 0)
      {
        this->vec_normale.Reallocate(Mat.GetNbPointsQuadratureBoundary());
        for (int kp = 0; kp < Mat.GetNbPointsQuadratureBoundary(); kp++)
          this->vec_normale(kp) = Mat.GetNormaleQuadratureBoundary(kp);
      }

    // atmospheric bc
    if (order == 21)
      {
        // condition -rho c^2 div \xi = rho (-i \omega + sigma)^2 delta (\xi \dot n)
        Real_wp rho = var_helm.ref_rho0(ref).GetConstant();
        Complex_wp m_iomega = -Iwp*var_problem.GetOmega() + var_helm.ref_sigma(ref).GetConstant(); 
        Complex_wp coef = rho*square(m_iomega)*delta;
        this->coef_phi = coef*nat_mat.GetCoefMass();
        return;     
      }
    
    if (order == 22)
      {
        // condition -rho c^2 div \xi = rho (-i \omega + sigma)^2 delta (1 - 2 \delta / r_t) (\xi \dot n)
        // courbure = 1 / (2 r_t)
        Real_wp rho = var_helm.ref_rho0(ref).GetConstant();
        Complex_wp m_iomega = -Iwp*var_problem.GetOmega() + var_helm.ref_sigma(ref).GetConstant(); 
        Real_wp k1 = Mat.GetK1QuadratureBoundary(k);
        Real_wp k2 = Mat.GetK2QuadratureBoundary(k);
        Real_wp courbure = 0.5*(k1+k2); // h = 1/R in 3-D, 1/(2R) in 2-D
        Complex_wp coef = rho*square(m_iomega)*delta*(1-4.0*delta*courbure);
        this->coef_phi = coef*nat_mat.GetCoefMass();
        return;
      }
    
    if (order ==23)
      {
        // condition -rho c^2 div \xi
        // = rho (-i \omega + sigma)^2 delta [1 - 2 \delta / r_t + 
        //             delta^2 (6/r_t^2 - (-i \omega + sigma)^2/c_0^2 + Laplace) ] (\xi \dot n)
        cout << " Atmospheric order 3 for Helmholtz H div not (yet) implemented" << endl;
        abort();
      }

    // absorbing bc   
    Real_wp mu = var_helm.ref_rho0(ref).GetConstant()*var_helm.ref_c0(ref).GetConstant();
    Complex_wp coef = -Iwp*var_problem.GetOmega()*mu*nat_mat.GetCoefDamping();
    this->coef_phi = coef;
  }


  /*******************
   * EllipticProblem *
   *******************/
  
  
  template<class Dimension>
  IncidentWaveProjector<Complex_wp, Dimension>*
  EllipticProblem<HelmholtzEquationHdiv<Dimension> >
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Complex_wp, Dimension>& u_inc) const
  {
    return new IncidentWaveProjector_HelmHdiv<Dimension>(*this, u_inc);
  }


  template<class Dimension>
  void EllipticProblem<HelmholtzEquationHdiv<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HelmholtzEquationHdiv<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelm(i, num_dof, mat_elem, nat_mat,
				      this->GetReferenceElementHdiv(i));
  }


  template<class Dimension>
  void EllipticProblem<LaplaceEquationHdiv<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelm(i, num_dof, mat_elem, nat_mat,
				      this->GetReferenceElementHdiv(i));
  }


  template<class Dimension>
  void EllipticProblem<LaplaceEquationHdiv<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
}

#define MONTJOIE_FILE_HELMHOLTZ_HDIV_CXX
#endif
