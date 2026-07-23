#ifndef MONTJOIE_FILE_VAR_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{

  /*********************
   * HelmholtzEquation *
   *********************/


  template<class T, class Dim>  
  inline bool HelmholtzEquation_Base<T, Dim>::SymmetricGlobalMatrix()
  {
    return true; 
  }
  

  template<class T, class Dim>
  inline bool HelmholtzEquation_Base<T, Dim>::SymmetricElementaryMatrix()
  {
    return true; 
  }
    
  
  template<class T, class Dim>
  inline bool HelmholtzEquation_Base<T, Dim>::ComputeDFjm1()
  {
    return store_dfjm1; 
  }


  template<class T, class Dimension> template<class TypeEquation>
  inline void HelmholtzEquation_Base<T, Dimension>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		      int i, const ElementReference_Dim<Dimension>& Fb)
  {
    ComputeMassMatrixGen(var, var, var, i,
                         dynamic_cast<const ElementReference<Dimension, 1>& >(Fb));
  }

  
  //! stiffness coefficients J_i DF_i^{-1} DF_i^{*-1} poids
  template<class T, class Dimension> template<class T0>
  inline TinyMatrix<T0, Symmetric,
		    Dimension::dim_N, Dimension::dim_N>
  HelmholtzEquation_Base<T, Dimension>::
  GetStiffCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
		      const GlobalGenericMatrix<T0>& nat_mat, int ref)
  {
    return TinyMatrix<T0, Symmetric, Dimension::dim_N, Dimension::dim_N>
      (nat_mat.GetCoefStiffness()*var.Glob_matMass_Bh(i)(j));
  } 
  
  
  //! mass coefficient \omega^2 \rho J_i poids
  template<class T, class Dimension> template<class T0>
  inline T0 HelmholtzEquation_Base<T, Dimension>
  ::GetMassCoefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
		       const GlobalGenericMatrix<T0>& nat_mat, int ref)
  {
    return var.Glob_matMass_Dh(i)(j)*nat_mat.GetCoefMass() 
      + var.Glob_matMass_DhSigma(i)(j)*nat_mat.GetCoefDamping();
  }
  

  template<class T, class Dimension> template<class T0>
  inline TinyVector<T0, Dimension::dim_N>
  HelmholtzEquation_Base<T, Dimension>
  ::GetNablaU_Coefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
			  const GlobalGenericMatrix<T0>& nat_mat, int ref)
  {
    TinyVector<T0, Dimension::dim_N> v;    
    if (var.include_flow_term)
      v = var.Glob_matDamp_Mh(i)(j)*nat_mat.GetCoefDamping();
    
    return v;
  }

  
  template<class T, class Dimension> template<class T0>
  inline TinyVector<T0, Dimension::dim_N>
  HelmholtzEquation_Base<T, Dimension>
  ::GetNablaPhi_Coefficient(const VarHelmholtz_Base<T, Dimension>& var, int i, int j,
			    const GlobalGenericMatrix<T0>& nat_mat, int ref)
  {
    TinyVector<T0, Dimension::dim_N> v;    
    if (var.include_flow_term)
      v = var.Glob_matDamp_Mh(i)(j)*-nat_mat.GetCoefDamping();
    
    return v;
  }
  
  
  //! mass coefficients
  template<class T, class Dimension> template<class T0, class MatStiff>
  inline void HelmholtzEquation_Base<T, Dimension>::
  GetMassPhiDFiTensor(const VarHelmholtz_Base<T, Dimension>& var, 
                      int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                      MatStiff& Amass, bool variable, bool affine,
		      const ElementReference_Dim<Dimension>& Fb)
  {
    Amass(0, 0) = var.Glob_matMass_Dh(i)(j)*nat_mat.GetCoefMass();
    Amass(0, 0) += var.Glob_matMass_DhSigma(i)(j)*nat_mat.GetCoefDamping();
  }
  
  
  //! fills stiffness tensor C
  template<class T, class Dimension> template<class TypeEquation, class T0, class MatStiff>
  inline void HelmholtzEquation_Base<T, Dimension>::
  GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                    int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                    int ref, MatStiff& Cgrad_grad)
  {
    TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>
      mu = vars.ref_mu(ref).GetCoefficient(vars, num_elem, jloc);
    
    if (vars.include_double_gradM_flow)
      {
        TinyVector<T, Dimension::dim_N>
          v0 = vars.ref_v0(ref).GetCoefficient(vars, num_elem, jloc);
        
        T beta = vars.ref_beta(ref).GetCoefficient(vars, num_elem, jloc);
        
        // we subtract v0 v0^T to mu
        Rank1Update(-beta, v0, mu);
      }
    
    Cgrad_grad(0, 0) = mu;
    Cgrad_grad *= nat_mat.GetCoefStiffness();
  }
  
  
  //! stiffness coefficients
  template<class T, class Dimension> template<class T0, class MatStiff>
  inline void HelmholtzEquation_Base<T, Dimension>::
  GetGradGradDFiTensor(const VarHelmholtz_Base<T, Dimension>& vars,
                       int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                       int ref, MatStiff& Cgrad_grad,
                       bool variable, bool affine,
		       const ElementReference_Dim<Dimension>& Fb)
  {
    Cgrad_grad(0, 0) = vars.Glob_matMass_Bh(num_elem)(jloc);
    Cgrad_grad(0, 0) *= nat_mat.GetCoefStiffness();
  }
  
  
  //! no terms in grad(u) v in Helmholtz equation
  template<class T, class Dimension> template<class T0, class MatStiff>
  void HelmholtzEquation_Base<T, Dimension>::
  GetGradPhiDFiTensor(const VarHelmholtz_Base<T, Dimension>& vars,
                      int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                      MatStiff& Dgrad_phi, MatStiff& Ephi_grad,
                      bool variable, bool affine,
		      const ElementReference_Dim<Dimension>& Fb)
  {
    if (vars.include_flow_term)
      {
	Dgrad_phi(0, 0) = vars.Glob_matDamp_Mh(num_elem)(jloc);
	Ephi_grad(0, 0) = vars.Glob_matDamp_Mh(num_elem)(jloc);
	Dgrad_phi(0, 0) *= -nat_mat.GetCoefDamping();
	Ephi_grad(0, 0) *= nat_mat.GetCoefDamping();
      }
    
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
  }
  
  
  template<class T, class Dim> template<class Matrix1, class TypeEquation, class T0>
  inline void HelmholtzEquation_Base<T, Dim>
  ::GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale,
               int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
               int ref, int ref2, const EllipticProblem<TypeEquation>& vars,
               const ElementReference<Dimension, 1>& Fb)
  {
    Nabc(0, 0) = vars.alpha_penalization*nat_mat.GetCoefStiffness();
  }
  
  
  template<class T, class Dim> template<class Vector1, class Vector2, class TypeEquation, class T0>
  inline void HelmholtzEquation_Base<T, Dim>
  ::MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector2& Un,
                         int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
                         int ref, int ref2,
                         const EllipticProblem<TypeEquation>& vars,
                         const ElementReference<Dimension, 1>& Fb)
  {
    Un(0) = Vn(0)*vars.alpha_penalization*nat_mat.GetCoefStiffness();
  }


  /*******************
   * VarHelmholtz_Base *
   *******************/


  //! default constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarHelmholtz_Base<Complexe, Dimension>
  ::VarHelmholtz_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_boundary(var)
  {
    InitDefaultValues();
  }


  template<class Complexe, class Dimension>
  inline TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N>
  VarHelmholtz_Base<Complexe, Dimension>::GetMu_QuadraturePoint(int ref, int i, int num_point) const
  {
    return ref_mu(ref).GetCoefficient(var_problem, i, num_point);
  }
  
  
  /*******************
   * VarHelmholtz_Eq *
   *******************/
  

  template<class TypeEquation>
  inline VarHelmholtz_Eq<TypeEquation>::VarHelmholtz_Eq()
    : VarHelmholtz_Base<Complexe, Dimension>(static_cast<EllipticProblem<TypeEquation>& >(*this))
  {
  }
    
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    VarHelmholtz_Base<Complexe, Dimension>::SetInputData(description_field, parameters);
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>
  ::CopyInputData(const VarProblem_Base& var_base)
  {
    const VarHelmholtz_Eq<TypeEquation>& var = static_cast<const VarHelmholtz_Eq<TypeEquation>& >(var_base);
    VarHarmonic<TypeEquation>::CopyInputData(var);
    VarHelmholtz_Base<Complexe, Dimension>::CopyInputData(var);
  }
    
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>::InitIndices(int n)
  {
    VarHelmholtz_Base<Complexe, Dimension>::InitIndices(n);
  }
  
  
  template<class TypeEquation>
  inline int VarHelmholtz_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarHelmholtz_Base<Complexe, Dimension>::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>
  ::SetIndices(int i, const VectString& parameters)
  {
    VarHelmholtz_Base<Complexe, Dimension>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarHelmholtz_Base<Complexe, Dimension>::SetPhysicalIndex(name_media, i, parameters);
  }
  

  template<class TypeEquation>
  inline string VarHelmholtz_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarHelmholtz_Base<Complexe, Dimension>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>::FinalizeComputationVaryingIndices()
  {
    VarHelmholtz_Base<Complexe, Dimension>::FinalizeComputationVaryingIndices();
  }
  
  
  template<class TypeEquation>
  inline bool VarHelmholtz_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarHelmholtz_Base<Complexe, Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  inline bool VarHelmholtz_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarHelmholtz_Base<Complexe, Dimension>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>::AllocateMassMatrices()
  {
    VarHelmholtz_Base<Complexe, Dimension>::AllocateMassMatrices();
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>::PutOtherGlobalDofs()
  {
    VarHelmholtz_Base<Complexe, Dimension>::PutOtherGlobalDofs();
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>
  ::MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Real_wp& coef) const
  { 
    coef *= realpart(this->ref_mu(ref).GetConstant()(0, 0));
  }
  
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>
  ::MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Complex_wp& coef) const
  { 
    coef *= this->ref_mu(ref).GetConstant()(0, 0);
  }
  
  
  //! if mix_ is true, we use a first-order formulation
  template<class TypeEquation>
  inline void VarHelmholtz_Eq<TypeEquation>::SetFirstOrderFormulation(bool mix_)
  {
    this->mixed_formulation = mix_;
    if (this->mixed_formulation)
      this->mesh_num.compute_dof_pml = true;
  }
  

  //! returns the velocity of waves in the physical media of reference ref
  template<class TypeEquation>
  inline Real_wp VarHelmholtz_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return sqrt(abs(this->ref_mu(ref)(0, 0)/this->ref_rho(ref).GetConstant()));
  }
  
  
  //! returns the velocity of waves at infinity
  template<class TypeEquation>
  inline Real_wp VarHelmholtz_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return sqrt(abs(this->mu0/this->rho0));
  }
  
  
  //! returns the coefficient |mu|
  template<class TypeEquation>
  inline Real_wp VarHelmholtz_Eq<TypeEquation>
  ::GetCoefficientPenaltyStiffness(int ref) const
  {
    return this->coefficient_mu_penalty(ref);
  }


  template<class TypeEquation>
  inline ModalSourceBoundary_Dim<typename TypeEquation::Complexe, typename TypeEquation::Dimension>* 
  VarHelmholtz_Eq<TypeEquation>::GetNewModalSourceEquation() const
  {
    return new ModalSourceLaplace<Complexe, Dimension>(this->GetLeafClass());
  }
  

  /*********************
   * VarHelmholtz_Cplx *
   *********************/
  

  template<class TypeEquation>
  inline void VarHelmholtz_Cplx<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    VarHelmholtz_Base<Complexe, Dimension>::GetVaryingIndices(rho_complex, num_ref, num_index,
                                                              num_component, compute_grad, compute_hess);
  }


  template<class TypeEquation>
  inline void VarHelmholtz_Cplx<TypeEquation>
  ::ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
                        const GridInterpolation<Dimension>& var_interp,
			int iquad, bool compute_grad) const
  {
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtz_Cplx<TypeEquation>
  ::ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
                        const GridInterpolation<Dimension>& var_interp,
			int iquad, bool compute_grad) const
  {
    if (TypeEquation::FirstOrderFormulation)
      {
        val_u(1) *= -this->GetOmega();
        val_u(2) *= -this->GetOmega();
        if (val_u.GetM() > 3)
          val_u(3) *= -this->GetOmega();
      }
  }


  template<class TypeEquation>
  inline void VarHelmholtz_Cplx<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			int i, bool u_component, bool grad_component) const 
  {
  }


  template<class TypeEquation>
  inline void VarHelmholtz_Cplx<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			int i, bool u_component, bool grad_component) const 
  {
    if (TypeEquation::FirstOrderFormulation)
      {
	if (u_component)
	  {
	    Unodal(1) *= -this->GetOmega();
	    Unodal(2) *= -this->GetOmega();
	    if (Unodal.GetM() > 3)
	      Unodal(3) *= -this->GetOmega();	
	  }
      }
  }
  
  template<class TypeEquation>
  Complex_wp VarHelmholtz_Cplx<TypeEquation>
  ::ComputeImpedanceCoefficient(const IVect& ref, VectComplex_wp& U)
  {
    return VarHelmholtz_Base<Complex_wp, Dimension>::ComputeImpedanceCoefficient(ref, U);
  }
  

  template<class T, class Dimension>
  inline CondensationBlockSolver<T, HelmholtzEquation<Dimension> >
  ::CondensationBlockSolver(EllipticProblem<HelmholtzEquation<Dimension> >& var)
    : CondensationBlockSolver_Helm<T, HelmholtzEquation<Dimension> >(var)
  {
  }
  
  

  /***********************
   * HelmholtzEquationDG *
   ***********************/

  
  template<class T, class Dim>
  inline bool HelmholtzEquationDG_Base<T, Dim>::SymmetricGlobalMatrix()
  {
    return true; 
  }
  
  
  template<class T, class Dim>
  inline bool HelmholtzEquationDG_Base<T, Dim>::SymmetricElementaryMatrix()
  {
    return false;
  }


  //! computation of mass coefficients
  template<class T, class Dimension> template<class TypeEquation>
  inline void HelmholtzEquationDG_Base<T, Dimension>::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                    int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }
  

  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline VarComputationRCS_Helm<Dimension>
  ::VarComputationRCS_Helm(EllipticProblem<TypeEquation>& var)
    : VarComputationRCS_Base<Dimension>(var), var_problem(var), var_source(var)
  {
  }


  template<class Dimension>
  inline VarComputationRCS<HelmholtzEquation<Dimension> >
  ::VarComputationRCS(EllipticProblem<HelmholtzEquation<Dimension> >& var)
    : VarComputationRCS_Helm<Dimension>(var) 
  {
  }
  
  
  template<class Dimension> template<class TypeEquation>    
  inline TransparencySolver_Helm<Dimension>
  ::TransparencySolver_Helm(EllipticProblem<TypeEquation>& var, All_LinearSolver& solver)
    : TransparencySolver_Fem<Dimension, 1, 1>(var, solver), var_problem(var), var_helm(var)
  {
  }  
  
  
  template<class Dimension>
  inline TransparencySolver<HelmholtzEquation<Dimension> >
  ::TransparencySolver(EllipticProblem<HelmholtzEquation<Dimension> >& var,
		       All_LinearSolver& solver) :
    TransparencySolver_Helm<Dimension>(var, solver) 
  {
  }
  
  
  template<class Dimension>
  inline VarComputationRCS<HelmholtzEquationDG<Dimension> >
  ::VarComputationRCS(EllipticProblem<HelmholtzEquationDG<Dimension> >& var)
    : VarComputationRCS_Helm<Dimension>(var)
  {
  }

  
  template<class Dimension>
  inline TransparencySolver<HelmholtzEquationDG<Dimension> >
  ::TransparencySolver(EllipticProblem<HelmholtzEquationDG<Dimension> >& var,
		       All_LinearSolver& solver) :
    TransparencySolver_Helm<Dimension>(var, solver) 
  {
  }
  

  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelm(i, num_dof, mat_elem, nat_mat, *this, *this,
					this->GetReferenceElementH1(i)); 
  }
  

  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquation<Dimension> >
  ::AddBoundaryConditionTerms(VirtualMatrix<Real_wp>& mat_sp,
                              const GlobalGenericMatrix<Real_wp>& nat_mat,
                              int offset_row, int offset_col)
  {
    cout << "not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquation<Dimension> >
  ::AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>& mat_sp,
                              const GlobalGenericMatrix<Complex_wp>& nat_mat,
                              int offset_row, int offset_col)
  {
    VarHarmonic<HelmholtzEquation<Dimension> >::
      AddBoundaryConditionTerms(mat_sp, nat_mat, offset_row, offset_col);
    
#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL
    this->TreatThinSlotCondition(mat_sp);
#endif      
  }


  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class T, class Complexe, class Dimension> template<class TypeEquation>
  inline CondensationBlockSolver_HelmDG<T, Complexe, Dimension>
  ::CondensationBlockSolver_HelmDG(EllipticProblem<TypeEquation>& var)
    : CondensationBlockSolver_Fem<T>(var), var_problem(var), var_helm(var)
  {    
  }
  
  
  template<class T, class Dimension>
  inline CondensationBlockSolver<T, HelmholtzEquationDG<Dimension> >
  ::CondensationBlockSolver(EllipticProblem<HelmholtzEquationDG<Dimension> >& var)
    : CondensationBlockSolver_HelmDG<T, Complex_wp, Dimension>(var)
  {
  }
  

  template<class Dimension>
  inline EllipticProblem<HelmholtzEquationDG<Dimension> >::EllipticProblem()
  {
  }
  
  
  template<class Dimension>
  inline IncidentWaveProjector<Complex_wp, Dimension>*
  EllipticProblem<HelmholtzEquationDG<Dimension> >
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Complex_wp, Dimension>& u_inc) const
  {
    return new IncidentWaveProjector_HelmDG<Dimension>(*this, u_inc);
  }
  
  
  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>& solver,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelmDG(i, num_dof, mat_elem, solver, nat_mat, *this, *this,
                                        this->GetReferenceElementH1(i)); 
  }
  

  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  inline void EllipticProblem<HelmholtzEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    this->AddElementaryFluxesHelmDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }

}

#define MONTJOIE_FILE_VAR_HELMHOLTZ_INLINE_CXX
#endif
