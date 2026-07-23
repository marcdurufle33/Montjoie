#ifndef MONTJOIE_FILE_AXISYM_HELMHOLTZ_INLINE_CXX

namespace Montjoie
{

  //! default constructor
  template<class Complexe> template<class TypeEquation>
  VarHelmholtz_Axi<Complexe>::VarHelmholtz_Axi(EllipticProblem<TypeEquation>& var)
    : VarAxisymProblem(var), var_problem(var), var_boundary(var), var_output(var), var_source(var),
      var_gibc_base(var.var_gibc)
  {
    InitDefaultValues();
  }


  template<class Complexe>
  inline TinyMatrix<Complexe, Symmetric, 2, 2>
  VarHelmholtz_Axi<Complexe>::GetMu_QuadraturePoint(int ref, int i, int num_point) const
  {
    const TinyMatrix<Complexe, Symmetric, 3, 3>& mu
      = this->ref_mu(ref).GetCoefficient(var_problem, i, num_point);
    
    TinyMatrix<Complexe, Symmetric, 2, 2> A;
    A(0, 0) = mu(0, 0); A(0, 1) = mu(0, 1);
    A(1, 1) = mu(1, 1);
    return A;
  }

  
  //! returns the number of physical indices
  template<class Complexe>
  inline int VarHelmholtz_Axi<Complexe>::GetNbPhysicalIndices() const
  {
    return ref_rho.GetM();
  }

  
  //! returns the velocity of waves in physical media whose reference is ref
  template<class Complexe>
  inline Real_wp VarHelmholtz_Axi<Complexe>::GetVelocityOfMedia(int ref) const
  {
    return sqrt(abs(this->ref_mu(ref).GetConstant()(0, 0)/this->ref_rho(ref).GetConstant()));
  }


  //! returns the velocity of waves at infinity
  template<class Complexe>
  inline Real_wp VarHelmholtz_Axi<Complexe>::GetVelocityOfInfinity() const
  {
    return 1.0;
  }


  //! returns rho + i sigma / omega on quadrature point
  template<>
  inline Complex_wp VarHelmholtz_Axi<Complex_wp>
  ::GetRhoTilde(int ref, int num_elem, int npoint) const
  {
    Complex_wp rho_tilde = this->ref_rho(ref).GetCoefficient(var_problem, num_elem, npoint);
    rho_tilde += Iwp*this->ref_sigma(ref).GetCoefficient(var_problem, num_elem, npoint)/var_problem.GetOmega();
    return rho_tilde;    
  }


  //! returns rho + i sigma / omega on quadrature point
  template<>
  inline Real_wp VarHelmholtz_Axi<Real_wp>
  ::GetRhoTilde(int ref, int num_elem, int npoint) const
  {
    return this->ref_rho(ref).GetCoefficient(var_problem, num_elem, npoint);
  }

  
  //! returns mu n \cdot n on quadrature point
  template<class Complexe>
  inline Complexe VarHelmholtz_Axi<Complexe>
  ::GetMuNormale(int ref, int num_elem, int npoint,
                 const R2& normale) const
  {
    Complexe mu = this->ref_mu(ref).GetCoefficient(var_problem, num_elem, npoint)(0, 0);
    return mu;
  }


  template<class Complexe>
  inline bool VarHelmholtz_Axi<Complexe>::UseFormulationR3() const
  {
    return formulation_r3; 
  }  

  
  /**********************
   * VarHelmholtzAxi_Eq *
   **********************/


  template<class TypeEquation>
  inline bool VarHelmholtzAxi_Eq<TypeEquation>
  ::UseNumericalIntegration(int i) const
  {
    return true; 
  }

  
  template<class TypeEquation>
  inline const R3& VarHelmholtzAxi_Eq<TypeEquation>::GetPhaseOrigin() const
  {
    return VarHelmholtz_Axi<Complexe>::GetPhaseOrigin();
  }
  
  
  template<class TypeEquation>
  inline int VarHelmholtzAxi_Eq<TypeEquation>::GetNbModesSource() const
  {
    return 1;
  }


  //! allocates an incident field
  template<class TypeEquation>
  inline IncidentWaveField<Complex_wp, Dimension2>* VarHelmholtzAxi_Eq<TypeEquation>
  ::GetNewIncidentField(int n, const Vector<VectString>&, const Complex_wp& val) const
  {
    return new PlaneWaveAxiIncidentField(this->kwave3D, this->rho0, this->sigma0,
					 this->mu0, this->GetOmega(), this->number_mode);
  }


  template<class TypeEquation>
  inline IncidentWaveField<Real_wp, Dimension2>* VarHelmholtzAxi_Eq<TypeEquation>
  ::GetNewIncidentField(int n, const Vector<VectString>&, const Real_wp& val) const
  {
    return NULL;
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    VarHelmholtz_Axi<Complexe>::SetInputData(description_field, parameters);
  }
  

  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>::CheckInputMesh()
  {
    VarHelmholtz_Axi<Complexe>::CheckSectionMeshAxi();
  }
  
  
  template<class TypeEquation>
  inline int VarHelmholtzAxi_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarHelmholtz_Axi<Complexe>::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  inline Real_wp VarHelmholtzAxi_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return VarHelmholtz_Axi<Complexe>::GetVelocityOfMedia(ref); 
  }
  
  
  template<class TypeEquation>
  inline Real_wp VarHelmholtzAxi_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return VarHelmholtz_Axi<Complexe>::GetVelocityOfInfinity();
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>::ConstructFiniteElement(const string& name_elt)
  {
    VarHarmonic<TypeEquation>::ConstructFiniteElement(name_elt);
    VarHelmholtz_Axi<Complexe>::ConstructFiniteElement(name_elt);
  }

  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>::InitIndices(int n)
  {
    VarHelmholtz_Axi<Complexe>::InitIndices(n);
  }

  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    VarHelmholtz_Axi<Complexe>::SetIndices(i, parameters);
  }
  

  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarHelmholtz_Axi<Complexe>::SetPhysicalIndex(name_media, i, parameters);
  }
  

  template<class TypeEquation>
  inline string VarHelmholtzAxi_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarHelmholtz_Axi<Complexe>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::SetPhysicalIndexAtInfinity(const Vector<bool>& RefUsed)
  {
    VarHelmholtz_Axi<Complexe>::SetPhysicalIndexAtInfinity(RefUsed);
  }

  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Real_wp& coef) const
  {
    coef *= realpart(this->ref_mu(ref).GetConstant()(0, 0));
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Complex_wp& coef) const
  { 
    coef *= this->ref_mu(ref).GetConstant()(0, 0);
  }

  
  template<class TypeEquation>
  inline bool VarHelmholtzAxi_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarHelmholtz_Axi<Complexe>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  inline bool VarHelmholtzAxi_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarHelmholtz_Axi<Complexe>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>::AllocateMassMatrices()
  {
    VarHelmholtz_Axi<Complexe>::AllocateMassMatrices();
  }
    
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
                           SetPoints<Dimension2>& PointsElem,
                           SetMatrices<Dimension2>& MatricesElem,
                           IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb)
  {
    VarHarmonic<TypeEquation>::ComputeLocalMassMatrix(i, N, linear_sparse, PointsElem,
                                                      MatricesElem, OrderFace, Fb);
    
    VarHelmholtz_Axi<Complexe>::ComputeLocalMassMatrix(i);
  }
  


  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::ConstructAll(const string& name_file, const string& name_element,
		 Vector<string>& lines_data_file)
  {
    VarHelmholtz_Axi<Complexe>::ConstructAll(name_file, name_element, lines_data_file);
  }


  template<class TypeEquation>
  inline void VarHelmholtzAxi_Eq<TypeEquation>
  ::ConstructAll(const Vector<string>& lines_data_file, const string& input_file,
                 const string& name_elt, const string& name_eq, All_LinearSolver*& solver,
                 bool compute_rho, bool delete_points, int num)
  {
    VarHelmholtz_Axi<Complexe>::ConstructAll(lines_data_file, input_file, name_elt, name_eq, solver,
                                             compute_rho, delete_points, num);
  }


  template<class TypeEquation>
  inline All_LinearSolver* VarHelmholtzAxi_Eq<TypeEquation>::GetNewLinearSolver()
  {
    return new All_LinearSolver(this->GetLeafClass());
  }


  /************************
   * VarHelmholtzAxi_Cplx *
   ************************/


  template<class TypeEquation>
  inline void VarHelmholtzAxi_Cplx<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    VarHelmholtz_Axi<Complexe>::GetVaryingIndices(rho_complex, num_ref, num_index,
                                                  num_component, compute_grad, compute_hess);
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Cplx<TypeEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Cplx<TypeEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelmAxi(i, num_dof, mat_elem, nat_mat,
					   this->GetReferenceElementH1(i));
  }

  
  template<class TypeEquation>
  inline IncidentWaveProjector<Complex_wp, Dimension2>*
  VarHelmholtzAxi_Cplx<TypeEquation>
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) const
  {
    return new IncidentWaveProjector_HelmAxi(this->GetLeafClass(), u_inc);
  }

  template<class TypeEquation>
  inline Complex_wp VarHelmholtzAxi_Cplx<TypeEquation>
  ::ComputeImpedanceCoefficient(const IVect& ref, VectComplex_wp& U)
  {
    return VarHelmholtz_Axi<Complex_wp>::ComputeImpedanceCoefficient(ref, U);
  }




  /************************
   * VarHelmholtzAxi_Real *
   ************************/

  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Real<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    VarHelmholtz_Axi<Complexe>::GetVaryingIndices(rho_real, num_ref, num_index,
                                                  num_component, compute_grad, compute_hess);
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Real<TypeEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
                            CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelmAxi(i, num_dof, mat_elem, nat_mat,
                                         this->GetReferenceElementH1(i));
  }
  
  
  template<class TypeEquation>
  inline void VarHelmholtzAxi_Real<TypeEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
                            CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHelmAxi(i, num_dof, mat_elem, nat_mat,
                                         this->GetReferenceElementH1(i));
  }
  

  /************************
   * HelmholtzEquationAxi *
   *************************/
  

  inline bool HelmholtzEquationAxi::SymmetricGlobalMatrix()
  {
    return true;
  }
  
  
  inline bool HelmholtzEquationAxi::SymmetricElementaryMatrix()
  {
    return false;
  }
  
  
  //! mass coefficient for non-linear Kerr term
  template<class TypeEquation, class T, class MatMass>
  inline void HelmholtzEquationAxi
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                  int i, int j, const GlobalGenericMatrix<T>& nat_mat, int ref, MatMass& mass)
  {
    if (vars.InsidePML(i))
      mass(0, 0) = vars.Glob_rtilde(i)(j)*vars.Glob_DrPML(i)(j)*vars.Glob_DzPML(i)(j);
    else
      mass(0, 0) = vars.Glob_rtilde(i)(j);
  }


  inline TinyMatrix<Complex_wp, Symmetric, 2, 2>
  HelmholtzEquationAxi
  ::GetStiffCoefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j, 
			const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref)
  {
    return TinyMatrix<Complex_wp, Symmetric, 2, 2>
      (nat_mat.GetCoefStiffness()*var.Glob_stiff(i)(j));
  }
  
  
  inline Complex_wp HelmholtzEquationAxi
  ::GetMassCoefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j,
		       const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref)
  {
    Real_wp m2 = square(var.GetCurrentModeNumber());
    Complex_wp coef = -var.GetSquareOmega()*var.Glob_mass_omega2(i)(j)
      - Complex_wp(0, var.GetOmega())*var.Glob_mass_iomega(i)(j)
      + m2*var.Glob_mass_m2(i)(j);

    if (var.include_flow_term)
      {
	Real_wp m = var.GetCurrentModeNumber();
	coef += m*var.GetOmega()*var.Glob_mass_m_omega(i)(j);
      }
    
    return coef;
  }
    
  
  inline TinyVector<Complex_wp, 2> HelmholtzEquationAxi
  ::GetNablaU_Coefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref)
  {
    TinyVector<Complex_wp, 2> v;
    if (var.include_flow_term)
      v = -Complex_wp(0, var.GetOmega())*var.Glob_grad_iomega(i)(j);

    if (var.include_double_gradM_flow)
      v -= Complex_wp(0, var.GetCurrentModeNumber())*var.Glob_grad_im(i)(j);
    
    return v;
  }

  
  inline TinyVector<Complex_wp, 2> HelmholtzEquationAxi
  ::GetNablaPhi_Coefficient(const EllipticProblem<HelmholtzEquationAxi>& var, int i, int j,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref)
  {
    TinyVector<Complex_wp, 2> v;
    if (var.include_flow_term)
      v = Complex_wp(0, var.GetOmega())*var.Glob_grad_iomega(i)(j);

    if (var.include_double_gradM_flow)
      v += Complex_wp(0, var.GetCurrentModeNumber())*var.Glob_grad_im(i)(j);
    
    return v;
  }
  
  
  /**********************
   * LaplaceEquationAxi *
   **********************/


  inline bool LaplaceEquationAxi::SymmetricGlobalMatrix()
  {
    return true;
  }
  
  
  inline bool LaplaceEquationAxi::SymmetricElementaryMatrix()
  {
    return false;
  }
  
  
  //! mass coefficient for non-linear Kerr term
  template<class TypeEquation, class T, class MatMass>
  inline void LaplaceEquationAxi
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                  int i, int j, const GlobalGenericMatrix<T>& nat_mat, int ref, MatMass& mass)
  {
    mass(0, 0) = vars.Glob_rtilde(i)(j);
  }


  /*************
   * Impedance *
   *************/
  
  
  template<class TypeEquation>
  inline ImpedanceABC<Complex_wp, HelmholtzEquationAxi>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceABC_Helm<Dimension2>(var), var_helm(var), var_boundary(var)
  {
  }


  template<class T> template<class TypeEquation>
  inline ImpedanceABC<T, LaplaceEquationAxi>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, Dimension2>(var), var_helm(var)
  {
  }

  
  template<class Complexe> template<class TypeEquation>
  inline ImpedanceGenericAxiHelm<Complexe>
  ::ImpedanceGenericAxiHelm(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<Complexe, Dimension2>(var), var_helm(var), var_boundary(var)    
  {
  }
  

  template<class TypeEquation>
  inline ImpedanceGeneric<Complex_wp, HelmholtzEquationAxi>
  ::ImpedanceGeneric(const EllipticProblem<TypeEquation>& var)
    : ImpedanceGenericAxiHelm<Complex_wp>(var)
  {
  }
  

  template<class TypeEquation>
  inline ImpedanceGeneric<Real_wp, LaplaceEquationAxi>
  ::ImpedanceGeneric(const EllipticProblem<TypeEquation>& var)
    : ImpedanceGenericAxiHelm<Real_wp>(var)
  {
  }


  /********************
   * Plane wave stuff *
   ********************/
  
  
  template<class TypeEquation>
  inline IncidentWaveProjector_HelmAxi
  ::IncidentWaveProjector_HelmAxi(const EllipticProblem<TypeEquation>& var,
                                  IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : IncidentWaveProjector<Complex_wp, Dimension2>(var, u_inc), var_helm(var)
  { 
    omega = var.GetOmega();
    number_mode = var.GetCurrentModeNumber(); 
    InitDefaultValues();
  }
  

  //! constructor with a given problem
  template<class TypeEquation>
  inline DiffractedWaveSource<HelmholtzEquationAxi>
  ::DiffractedWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
			 IncidentWaveField<Complex_wp, Dimension2>& fsrc_)
    : DiffractedWaveSource_HelmAxi(vars_helm, fsrc_)
  {
  }

  //! constructor with given problem
  template<class TypeEquation>
  inline DiffractedWaveSource_HelmAxi
  ::DiffractedWaveSource_HelmAxi(const EllipticProblem<TypeEquation>& var,
				 IncidentWaveField<Complex_wp, Dimension2>& u_inc)
    : VirtualSourceFEM<Complex_wp, Dimension2>(var), incident_wave(u_inc),
      var_problem(var), var_helm(var), var_boundary(var), fct_imped(var.GetGenericImpedanceFunction())
  {
    InitDefaultValues();
  }
  

  //! constructor with given problem
  template<class TypeEquation>
  inline TotalWaveSource_HelmAxi
  ::TotalWaveSource_HelmAxi(const EllipticProblem<TypeEquation>& var,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) 
    : VirtualSourceFEM<Complex_wp, Dimension2>(var), incident_wave(u_inc), var_problem(var),
      var_helm(var), var_boundary(var)
  {
    InitDefaultValues();
  }


  //! constructor with given problem
  template<class TypeEquation>
  inline TotalWaveSource<HelmholtzEquationAxi>
  ::TotalWaveSource(const EllipticProblem<TypeEquation>& vars_helm,
		    IncidentWaveField<Complex_wp, Dimension2>& fsrc_)
    : TotalWaveSource_HelmAxi(vars_helm, fsrc_) 
  {
  }


  template<class TypeEquation>
  inline VarGeneralizedImpedance<HelmholtzEquationAxi> 
  ::VarGeneralizedImpedance(EllipticProblem<TypeEquation>& var)
    : VarGeneralizedImpedance_Helm<Dimension2>(var) 
  {
  }


  //! constructor with given problem
  template<class TypeEquation>
  inline VolumetricSource<HelmholtzEquationAxi>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_AxiHelm<Complex_wp>(var, param) 
  {
  }
  
  
  //! constructor with given problem
  template<class TypeEquation>
  inline VolumetricSource<LaplaceEquationAxi>
  ::VolumetricSource(const EllipticProblem<TypeEquation>& var,
		     const Vector<VectString>& param)
    : VolumetricSource_AxiHelm<Real_wp>(var, param) 
  {
  }


  inline FemMatrixFreeClass<Complex_wp, HelmholtzEquationAxi>
  ::FemMatrixFreeClass(const EllipticProblem<HelmholtzEquationAxi>& var_) 
    : FemMatrixFreeClass_Eq<Complex_wp, HelmholtzEquationAxi>(var_)
  {
  }
  
}

#define MONTJOIE_FILE_AXISYM_HELMHOLTZ_INLINE_CXX
#endif
