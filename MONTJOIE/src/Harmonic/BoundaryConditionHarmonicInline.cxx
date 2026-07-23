#ifndef MONTJOIE_FILE_BOUNDARY_CONDITION_HARMONIC_INLINE_CXX

namespace Montjoie
{  

  /*****************************
   * VarBoundaryCondition_Base *
   *****************************/


  //! default constructor
  template<class TypeEquation>
  inline VarBoundaryCondition_Base
  ::VarBoundaryCondition_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(static_cast<VarProblem_Base&>(var)),
      var_computation(static_cast<VarComputationProblem&>(var)),
      var_comm(var), var_source(var)
  {
    InitDefaultValues();
  }
  
  
  //! returns the number of degrees of freedom with Dirichlet condition
  inline int VarBoundaryCondition_Base::GetNbDirichletDof() const
  {
    return nb_dof_dirichlet;
  }


  //! returns the number of degrees of freedom with Dirichlet condition
  inline int VarBoundaryCondition_Base::GetNbGlobalDirichletDof() const
  {
    return nb_dof_dirichlet_all;
  }
  
  
  //! returns dof number of the Dirichlet dof i
  inline int VarBoundaryCondition_Base::GetDirichletDofNumber(int i) const
  {
    return Dirichlet_dof(i);
  }
  
  
  //! returns true if the dof i is a dof with Dirichlet condition
  inline bool VarBoundaryCondition_Base::IsDofDirichlet(int i) const
  {
    return is_dof_dirichlet(i);
  }


  //! returns arrays is_dof_dirichlet
  inline const Vector<bool>& VarBoundaryCondition_Base::GetIsDofDirichlet() const
  {
    return is_dof_dirichlet;
  }
  
  
  //! returns true if the Dirichlet condition should be treated to obtain a symmetric matrix
  inline bool VarBoundaryCondition_Base::UseSymmetricDirichlet() const
  {
    return matrix_symmetric_dir;
  }


  //! forces symmetric treatment of Dirichlet condition (even for a non-symmetric matrix)
  inline void VarBoundaryCondition_Base::EnableSymmetricDirichlet(bool sym)
  {
    matrix_symmetric_dir = sym;
  }
  

  //! returns the coefficient alpha to put in the rows alpha u = alpha u_i (Dirichlet condition)
  inline Real_wp VarBoundaryCondition_Base::GetCoefficientDirichlet() const
  {
    return coef_dirichlet_matrix;
  }


  //! sets the coefficient alpha to put in the rows alpha u = alpha u_i (Dirichlet condition)
  inline void VarBoundaryCondition_Base::SetCoefficientDirichlet(const Real_wp& coef)
  {
    coef_dirichlet_matrix = coef;
  }

  
  //! informs that dof i is a Dirichlet dof or not
  inline void VarBoundaryCondition_Base::SetDirichletDof(int i, bool b)
  {
    is_dof_dirichlet(i) = b;
  }
  
  
  //! returns dof numbers with Dirichlet condition
  inline const IVect& VarBoundaryCondition_Base::GetDirichletDofNumber() const
  {
    return Dirichlet_dof;
  }
  
  
  //! returns the number of components set to 0 for supported condition
  inline int VarBoundaryCondition_Base::GetNbSupportedComponents(int ref) const
  {
    return supported_components_BC(ref).GetM();
  }
  
  
  //! returns component number for supported condition
  inline int VarBoundaryCondition_Base::GetSupportedComponent(int ref, int i) const
  {
    return supported_components_BC(ref)(i);  
  }


  //! sets the component numbers for supported condition
  inline void VarBoundaryCondition_Base::SetSupportedComponents(int ref, const IVect& comp)
  {
    supported_components_BC(ref) = comp;
  }

  
  //! sets homogeneous Dirichlet condition x_sol = 0 for Dirichlet dofs
  template<class T>
  inline void VarBoundaryCondition_Base::ImposeNullDirichletCondition(Vector<T>& x_sol) const
  {
    // Dirichlet dofs
    for (int i = 0; i < this->nb_dof_dirichlet; i++)
      x_sol(this->Dirichlet_dof(i)) = 0.0;
  }
  
  
  //! return order for high conductivity condition
  inline int VarBoundaryCondition_Base::GetHighConductivityOrder() const
  {
    return order_high_conductivity;
  }

  
  //! returns the number of modes
  /*!
    It returns the number of modes to compute in order to recover the solution
    for a source that has no symmetry (only the computational domain has symmetry).
    For example, in case of cyclic domains, it is the number of Fourier modes.
  */
  inline int VarBoundaryCondition_Base::GetNbModes() const
  {
    return list_number_mode.GetM();
  }
  

  //! returns the number of modes
  /*!
    It returns the number of modes to compute in order to recover the solution
    for a source that has no symmetry (only the computational domain has symmetry).
    For example, in case of cyclic domains, it is the number of Fourier modes.
  */
  inline int VarBoundaryCondition_Base::GetNbModesSource() const
  {
    return list_number_mode.GetM();
  }

  
  //! returns mode number of mode n
  inline int VarBoundaryCondition_Base::GetModeNumber(int n) const
  {
    return list_number_mode(n);
  }
  
  
  //! returns current mode number solved
  inline int VarBoundaryCondition_Base::GetCurrentModeNumber() const
  {
    return number_mode;
  }

  
  //! sets the current mode number
  inline void VarBoundaryCondition_Base::SetCurrentModeNumber(int n)
  {
    number_mode = n; 
  }
  
  
  //! returns true if the modes are not stored
  /*!
    When the modes are not stored, the final solution is modified at each computation.
    If the modes are stored, fft can be used to obtain quickly the final solution
  */
  inline bool VarBoundaryCondition_Base::ModesNotStored() const
  {
    return do_not_store_modes;
  }

  
  //! if true, the modes are stored (=> each mode is displayed)
  inline void VarBoundaryCondition_Base::ForceStorageModes(bool store)
  {
    do_not_store_modes = !store;
  }
  

  //! returns the type of symmetry for the domain
  /*!
    Several symmetries are possible (periodicity in x, y, z and/or theta)
  */
  inline int VarBoundaryCondition_Base::GetSymmetryType() const
  {
    if (plane_wave_quasi_periodic)
      return NO_SYMMETRY;
    
    return type_symmetry;
  }


  //! returns the number of periodic dofs
  inline int VarBoundaryCondition_Base::GetNbPeriodicDof(int n) const
  {
    return var_problem.GetMeshNumberingBase(n).GetNbPeriodicDof();
  }


  //! returns the dof number of periodic dof i
  inline int VarBoundaryCondition_Base::GetPeriodicDof(int i, int n) const
  {
    return var_problem.GetMeshNumberingBase(n).GetPeriodicDof(i);
  }
  
  
  //! returns the dof number of original periodic dof i
  inline int VarBoundaryCondition_Base::GetOriginalPeriodicDof(int i, int n) const
  {
    return var_problem.GetMeshNumberingBase(n).GetOriginalPeriodicDof(i);
  }


  //! returns the formulation used for periodic dofs
  inline int VarBoundaryCondition_Base::GetFormulationForPeriodicCondition() const
  {
    return var_problem.GetMeshNumberingBase().GetFormulationForPeriodicCondition();
  }
  

  //! sets the list of modes to compute
  inline void VarBoundaryCondition_Base::SetModesToCompute(const IVect& num)
  {
    this->list_number_mode = num;
  }
  
  
  //! adds mode n at the end of the list of modes
  inline void VarBoundaryCondition_Base::PushBackMode(int n)
  {
    this->list_number_mode.PushBack(-n);
    this->list_number_mode.PushBack(n);
  }

  
  //! returns the number of periodic modes in x-coordinate
  inline int VarBoundaryCondition_Base::GetNbPeriodicModesX() const
  {
    return nb_modes_periodic_xyz(0);
  }


  //! returns the number of periodic modes in y-coordinate
  inline int VarBoundaryCondition_Base::GetNbPeriodicModesY() const
  {
    return nb_modes_periodic_xyz(1);
  }

  
  //! returns the number of periodic modes in z-coordinate
  inline int VarBoundaryCondition_Base::GetNbPeriodicModesZ() const
  {
    return nb_modes_periodic_xyz(2);
  }
  

  //! returns order of absorbing boundary condition
  inline int VarBoundaryCondition_Base::GetOrderAbsorbingCondition() const
  {
    return order_ABC;
  }


  //! returns the number of elements contained in PML layers
  inline int VarBoundaryCondition_Base::GetNbEltPML() const
  {
    return this->nb_eltPML;
  }


  //! returns the global number of elements contained in PML layers
  inline int VarBoundaryCondition_Base::GetNbGlobalEltPML() const
  {
    return this->nb_eltPML_all;
  }

  
  /*inline Real_wp VarBoundaryCondition_Base::GetOmegaPML() const
  {
    return this->omega_pml;
  }


  inline void VarBoundaryCondition_Base::SetOmegaPML(const Real_wp& om)
  {
    this->omega_pml = om;
    }*/


  inline void VarBoundaryCondition_Base
  ::MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Real_wp& coef) const
  {
    // default value of mu is one, coef is not modified
  }
  
  
  inline void VarBoundaryCondition_Base
  ::MltMuIntegrationByParts(int ref, int ne, int num_loc, int k, Complex_wp& coef) const
  {
    // default value of mu is one, coef is not modified
  }

  
  //! ?
  inline void VarBoundaryCondition_Base::SetPhysicalIndexAtInfinity(const Vector<bool>&)
  {
  }    
  
  
  //! returns the maximal velocity in PML layers
  inline const Real_wp& VarBoundaryCondition_Base::GetMaximumVelocityPML() const
  {
    return max_velocity_pml;
  }
  
  
  /****************************
   * VarBoundaryCondition_Dim *
   ****************************/


  //! default constructor
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarBoundaryCondition_Dim<Complexe, Dimension>
  ::VarBoundaryCondition_Dim(EllipticProblem<TypeEquation>& var)
    : VarBoundaryCondition_Base(var), mesh_reference(var.mesh), var_comm(var)
  {
    InitDefaultValues();
  }

  
  //! phase for quasi-periodic condition
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetPeriodicPhase(int num_point, Real_wp& phase) const
  {
    // real numbers -> no phase
    SetComplexOne(phase);
  }


  //! phase for quasi-periodic condition
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetPeriodicPhase(int num_point, Complex_wp& phase) const
  {
    to_complex(periodic_dof_phase(num_point), phase);
  }
  
  
  //! retrieves the dof numbers associated with a periodic dof (original number and periodic number)
  template<class Complexe, class Dimension>
  void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetPeriodicDofNumbers(int i, int& num_per, int& num_orig, int n) const
  {
    num_per = var_comm.GetMeshNumbering(n).GetPeriodicDof(i);
    num_orig = var_comm.GetMeshNumbering(n).GetOriginalPeriodicDof(i);
  }
  

  //! returns modifiction of impedance coefficient for absorbing condition
  template<class Complexe, class Dimension> 
  inline const Complexe&
  VarBoundaryCondition_Dim<Complexe, Dimension>::GetImpedanceCoefficientABC() const
  {
    return Coef_Imped_ABC;
  }
  
  
  //! returns the damping coefficient in PML layer
  /*!
    \param[in] i element number in the PML
    \param[in] j quadrature point number
    This method returns (tau_x, tau_y, tau_z), the damping on each coordinate
  */
  template<class Complexe, class Dimension>
  inline const TinyVector<Complexe, Dimension::dim_N>&
  VarBoundaryCondition_Dim<Complexe, Dimension>::GetTauPML(int i, int j) const
  {
    return Glob_TauPML(i)(j);
  }

  
  //! returns the damping coefficient in PML layer
  /*!
    \param[in] i element number in the PML
    \param[in] j quadrature point number
    \param[in] m coordinate number
    This method returns tau_m, the damping on coordinate m
  */
  template<class Complexe, class Dimension>
  inline const Complexe&
  VarBoundaryCondition_Dim<Complexe, Dimension>::GetTauPML(int i, int j, int m) const
  {
    return Glob_TauPML(i)(j)(m);
  }


  //! returns the primitive of damping coefficient in PML layer
  /*!
    \param[in] i element number in the PML
    \param[in] j quadrature point number
    \param[in] m coordinate number
    This method returns \int tau_m, the primitive of damping on coordinate m
  */
  template<class Complexe, class Dimension>
  inline const Complexe&
  VarBoundaryCondition_Dim<Complexe, Dimension>::GetPrimitiveTauPML(int i, int j, int m) const
  {
    return Glob_TauPrimePML(i)(j)(m);
  }


  //! evaluates the damping coefficient in PML layer
  /*!
    \param[in] point point where the damping tau is evaluated
    \param[out] zeta damping tau
    \param[out] zeta primitive of damping tau
  */
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetDampingFactorPML(const R_N& point, int num_pml,
                        int type_pml, R_N& zeta, R_N& zeta_prime)
  {
    return VarBoundaryCondition_Base::
      GetDampingTauPML(mesh_reference, num_pml, type_pml, point, zeta, zeta_prime);
  }


  //! returns the parameters given in the data file for ConditionReference
  template<class Complexe, class Dimension>
  inline Vector<Vector<Complexe> >& VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetParamCondition()
  {
    return Param_condition;
  }
  
  
  //! returns the parameter j of reference i
  /*!
    Parameters of boundary conditions are parameters located after specifying
    boundary conditions in the data file. For instance, if you write in the data file
    ConditionReference = 1 DIRICHLET 3.0 4.0
    The parameters for reference 1 are 3.0 and 4.0
  */
  template<class Complexe, class Dimension>
  inline const Complexe&
  VarBoundaryCondition_Dim<Complexe, Dimension>::GetParamCondition(int i, int j) const
  {
    return Param_condition(i)(j);
  }


  //! Multiplies the coefficient coef by the impedance given in the data file
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::MltParamCondition(int ref, int k, Real_wp& coef) const
  {
    coef *= realpart(this->GetParamCondition(ref, k));
  }
  
  
  //! Multiplies the coefficient coef by the impedance given in the data file
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::MltParamCondition(int ref, int k, Complex_wp& coef) const
  {
    coef *= this->GetParamCondition(ref, k);
  }

  
  //! sets boundary condition for reference ref
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>::SetBoundaryConditionMesh(int ref, int type)
  {
    mesh_reference.SetBoundaryCondition(ref, type);
  }
  
  
  //! adds periodic boundary condition between two references
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::AddPeriodicConditionMesh(const TinyVector<int, 2>& ref, int type)
  {
    mesh_reference.AddPeriodicCondition(ref, type);
  }

  
  //! returns damping factor \zeta inside PML layers for x=point
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition_Dim<Complexe, Dimension>
  ::GetDampingFactorPML(const R_N& point, int num_pml, int type_pml,
                        R_N_Complex_wp& zeta, R_N_Complex_wp& zeta_prime)
  {
    return VarBoundaryCondition_Base::
      GetDampingFactorPML(mesh_reference, num_pml, type_pml, point, zeta, zeta_prime);
  }
  
    
  /************************
   * VarBoundaryCondition *
   ************************/
  
  
  //! constructor with a given problem
  template<class Complexe, class Dimension> template<class TypeEquation>
  inline VarBoundaryCondition<Complexe, Dimension>
  ::VarBoundaryCondition(EllipticProblem<TypeEquation>& var)
    : VarBoundaryCondition_Dim<Complexe, Dimension>(var), var_problem(var)
#ifdef MONTJOIE_WITH_TRANSMISSION
    ,var_transmission_base(var.var_transmission)
#endif
    , var_gibc_base(var.var_gibc)
  {
  }


#ifdef MONTJOIE_WITH_TRANSMISSION
  template<class Complexe, class Dimension>
  inline const VarTransmission_Base<Dimension>& 
  VarBoundaryCondition<Complexe, Dimension>::GetTransmissionProblem() const
  {
    return var_transmission_base; 
  }
#endif
  

  //! finds degrees of freedom associated with a set of references
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition<Complexe, Dimension>
  ::FindDofsOnReference(const Vector<int>& ref, Vector<int>& Dofs)
  {
    Vector<int> ref_cond(var_problem.mesh.GetNbReferences()+1);
    ref_cond.Zero();
    for (int i = 0; i < ref.GetM(); i++)
      ref_cond(ref(i)) = 1;
    
    VarBoundaryCondition_Base::FindDofsOnReference(var_problem, ref_cond, 1, Dofs);
  }

  
  //! finds degrees of freedom with Dirichlet condition
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition<Complexe, Dimension>::TreatDirichletCondition()
  {
    VarBoundaryCondition_Base::TreatDirichletCondition(var_problem);
  }
  
  
  //! Adds a term alpha \int f(phi_i, dphi_i) phi_j + g(phi_i, dphi_i) dphi_j  
  //! to the sparse matrix mat_sp
  /*!
    The definition of f and g is supplied by the user through the object fimpedance
    this object should derive from ImpedanceFunction_Base
    and overload virtual methods EvaluateImpedancePhi (evaluates f), EvaluateImpedanceGrad (evaluates g),
    ApplyImpedancePhi (computes f(phi_i, dphi_i)) and ApplyImpedanceGrad (computes g(phi_i, dphi_i))
  */
  template<class Complexe, class Dimension> template<class T>
  inline void VarBoundaryCondition<Complexe, Dimension>
  ::AddMatrixImpedanceBoundary(const T& alpha, const Vector<int>& ref_cond, int ref_target,
			       const GlobalGenericMatrix<T>& nat_mat,
			       VirtualMatrix<T>& mat_sp, int offset_row, int offset_col,
			       ImpedanceFunction_Base<T, Dimension>& fimpedance,
			       bool changeColonne, bool changeLigne)
  {
    VarBoundaryCondition_Dim<Complexe, Dimension>::
      AddMatrixImpedanceBoundary(alpha, ref_cond, ref_target, nat_mat, mat_sp,
				 offset_row, offset_col, fimpedance, changeColonne, changeLigne,
				 var_problem);
  }

  
  //! treats PML layers
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition<Complexe, Dimension>::FindElementsInsidePML()
  {
    VarBoundaryCondition_Base::FindElementsInsidePML(var_problem);
  }

  
  //! adds terms due to boundary conditions
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition<Complexe, Dimension>
  ::AddBoundaryConditionTerms(VirtualMatrix<Real_wp>& mat_sp,
			      const GlobalGenericMatrix<Real_wp>& nat_mat,
			      int offset_row, int offset_col)
  {
    AddBoundaryConditionGen(mat_sp, nat_mat, offset_row, offset_col);
  }


  //! adds terms due to boundary conditions  
  template<class Complexe, class Dimension>
  inline void VarBoundaryCondition<Complexe, Dimension>
  ::AddBoundaryConditionTerms(VirtualMatrix<Complex_wp>& mat_sp,
			      const GlobalGenericMatrix<Complex_wp>& nat_mat,
			      int offset_row, int offset_col)
  {
    AddBoundaryConditionGen(mat_sp, nat_mat, offset_row, offset_col);
  }

  
  /**************************
   * ImpedanceFunction_Base *
   **************************/
  

  //! default constructor
  template<class T, class Dimension>
  inline VirtualImpedanceFunction<T, Dimension>
  ::VirtualImpedanceFunction(const DistributedProblem<Dimension>& var)
    : var_problem(var), var_boundary(var.GetBoundaryConditionProblem())
  {
    InitDefaultValues();
  }


  //! destructor
  template<class T, class Dimension>
  inline VirtualImpedanceFunction<T, Dimension>::~VirtualImpedanceFunction()
  {
  }
  
  
  //! returns true if only dofs on the surface are involved
  template<class T, class Dimension>
  inline bool VirtualImpedanceFunction<T, Dimension>::InvolveOnlyTangentialDofs() const
  {
    return only_surface_dof_involved;
  }


  //! sets the impedance coefficients
  template<class T, class Dimension>
  inline void VirtualImpedanceFunction<T, Dimension>::SetCoefficient(const T& a, const T& b)
  {
    read_param_condition = false;
    coef_phi = a;
    coef_grad = b;
  }
  

  //! default constructor
  template<class T>
  inline ImpedanceFunction_Base<T, Dimension2>
  ::ImpedanceFunction_Base(const VarProblem<Dimension2>& var)
    : VirtualImpedanceFunction<T, Dimension2>(var),
      var_problem(var), var_boundary(var.GetBoundaryConditionProblem())
  {
  }


  //! default constructor
  template<class T>
  inline ImpedanceFunction_Base<T, Dimension3>
  ::ImpedanceFunction_Base(const VarProblem<Dimension3>& var)
    : VirtualImpedanceFunction<T, Dimension3>(var),
      var_problem(var), var_boundary(var.GetBoundaryConditionProblem())
  {
  }
  

  /**************************************
   * ImpedanceGeneric and other objects *
   **************************************/

  
  //! default constructor
  template<class T, class TypeEquation>
  ImpedanceGeneric<T, TypeEquation>
  ::ImpedanceGeneric(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, typename TypeEquation::Dimension>(var)
  {
  }
  
  
  //! default constructor
  template<class T, class TypeEquation>
  ImpedanceABC<T, TypeEquation>
  ::ImpedanceABC(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, typename TypeEquation::Dimension>(var)
  {
  }


#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
  //! default constructor
  template<class T, class TypeEquation>
  ImpedanceHighConductivity<T, TypeEquation>
  ::ImpedanceHighConductivity(const EllipticProblem<TypeEquation>& var)
    : ImpedanceFunction_Base<T, typename TypeEquation::Dimension>(var)
  {
  }
#endif
  
}

#define MONTJOIE_FILE_BOUNDARY_CONDITION_HARMONIC_INLINE_CXX
#endif
