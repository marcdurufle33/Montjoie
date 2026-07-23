#ifndef MONTJOIE_FILE_VAR_PROBLEM_BASE_INLINE_CXX

namespace Montjoie
{
  
  /*******************
   * VarProblem_Base *
   *******************/
  
  
  //! default constructor
  template<class TypeEquation>
  inline VarProblem_Base::VarProblem_Base(EllipticProblem<TypeEquation>& var)
    : nb_unknowns_scal(TypeEquation::nb_unknowns_scal),
      nb_unknowns_vec(TypeEquation::nb_unknowns_vec), nb_unknowns(TypeEquation::nb_unknowns),
      nb_components_en(TypeEquation::nb_components_en), nb_components_hn(TypeEquation::nb_components_hn),
      nb_unknowns_hdg(TypeEquation::nb_unknowns_hdg),
      type_element(TypeEquation::type_element), other_type_element(TypeEquation::GetOtherElementType()),
      first_order_formulation(TypeEquation::FirstOrderFormulation)
  {
    InitDefaultValues();
    
    this->compute_dfjm1 = TypeEquation::ComputeDFjm1();    
  }
  
  
  //! destructor
  inline VarProblem_Base::~VarProblem_Base()
  {
  }
  

  //! returns the number of degrees of freedom of the problem
  /*!
    In sequential, this number will represent the size of the 
    finite element matrix. It comprises all the degrees of freedom
    (mesh + boundary conditions + other models)
    
    In parallel, this is the size of the local finite element
    matrix (associated with the current processor)
   */
  inline int VarProblem_Base::GetNbDof() const
  {
    return this->nodl;
  }


  //! sets the number of degrees of freedom of the problem
  inline void VarProblem_Base::SetNbDof(int n)
  {
    this->nodl = n;
  }


  //! returns the number of mesh numberings
  inline int VarProblem_Base::GetNbMeshNumberings() const
  {
    return this->all_mesh_num.GetM();
  }  

  
  //! returns the number of degrees of freedom for numbering n
  inline int VarProblem_Base::GetNbMeshDof(int n) const
  {
    return this->all_mesh_num(n)->GetNbDof();
  }


  //! returns the mesh numbering n
  inline const MeshNumbering_Base<Real_wp>& VarProblem_Base::GetMeshNumberingBase(int n) const
  {
    return *this->all_mesh_num(n);
  }


  //! returns the offset associated with unknown m
  inline int VarProblem_Base::GetOffsetDofUnknown(int m) const
  {
    return offset_dof_unknown(m);
  }


  //! returns the offset associated with unknown m (with static condensation)
  inline int VarProblem_Base::GetOffsetDofCondensed(int m) const
  {
    return offset_dof_condensed(m);
  }


  //! returns the number of dofs for unknown m
  inline int VarProblem_Base::GetNbDofUnknown(int m) const
  {
    return offset_dof_unknown(m+1) - offset_dof_unknown(m);
  }
  

  //! returns the default order (used to subdivide the mesh for outputs)
  inline int VarProblem_Base::GetDefaultOrder() const
  {
    return default_order;
  }
  

  //! returns true if a Discontinuous Galerkin formulation is used
  inline int VarProblem_Base::FormulationDG() const
  {
    return dg_formulation;
  }

  
  //! returns true if inverse of jacobian matrix DF_i^-1 is computed and stored
  inline bool VarProblem_Base::ComputeDFjm1() const
  {
    return compute_dfjm1;
  }


  //! returns true if a first-order formulation (for continuous elements)
  /*!
    The first-order formulation of an equation is obtained by adding unknowns
    For example, Laplace equation Delta u = 0
    is transformed in div v = 0, v = grad u
    For Discontinuous Galerkin, if true the PML are rewritten to obtain
    a linear eigenvalue problem (in omega)
   */
  inline bool VarProblem_Base::FirstOrderFormulation() const
  {
    return this->mixed_formulation;
  }


  //! returns true if a first-order formulation (with only first-order derivatives) is used 
  /*!
    The first-order formulation of an equation is obtained by adding unknowns
    For example, Laplace equation Delta u = 0
    is transformed in div v = 0, v = grad u
   */
  inline bool VarProblem_Base::FirstOrderFormulationDG() const
  {
    return this->first_order_formulation;
  }
  

  //! tells to use a first-order formulation with only first-order derivatives
  inline void VarProblem_Base::SetFirstOrderFormulation(bool mix_)
  {
    this->mixed_formulation = mix_;
  }
  
  
  //! returns true if exact integration is performed
  /*!
    If true, some Discontinuous Galerkin formulation (e.g. in aeroacoustics)
    can be written without making the integration by parts. The
    induced stiffness matrix will not be explicitely skew-symmetric
    (i.e. skew-symmetric if the integrals are exact)
   */
  inline bool VarProblem_Base::UseExactIntegrationElement() const
  {
    return dg_exact_integration;
  }


  //! returns the over-integration used to compute integrals
  inline int VarProblem_Base::GetOverIntegration() const
  {
    return order_over_integration;
  }
  

  //! returns minimum of x-coordinates of the physical domain
  inline const Real_wp& VarProblem_Base::GetXmin() const
  {
    return xmin;
  }
  
  
  //! returns maximum of x-coordinates of the physical domain
  inline const Real_wp& VarProblem_Base::GetXmax() const
  {
    return xmax;
  }
  

  //! returns minimum of y-coordinates of the physical domain  
  inline const Real_wp& VarProblem_Base::GetYmin() const
  {
    return ymin;
  }
  
  
  //! returns maximum of y-coordinates of the physical domain
  inline const Real_wp& VarProblem_Base::GetYmax() const
  {
    return ymax;
  }


  //! returns minimum of z-coordinates of the physical domain
  inline const Real_wp& VarProblem_Base::GetZmin() const
  {
    return zmin;
  }
  
  
  //! returns maximum of z-coordinates of the physical domain
  inline const Real_wp& VarProblem_Base::GetZmax() const
  {
    return zmax;
  }

  
  //! sets boundaries of the computational domain
  inline void VarProblem_Base
  ::SetComputationalDomain(const Real_wp& x0, const Real_wp& x1, const Real_wp& y0,
			   const Real_wp& y1, const Real_wp& z0, const Real_wp& z1)
  {
    xmin = x0; xmax = x1; ymin = y0; ymax = y1; zmin = z0; zmax = z1;
  }


  //! returns the square of the pulsation
  inline const Real_wp& VarProblem_Base::GetSquareOmega() const
  {
    return this->omega2;
  }

 
  //! returns the pulsation omega
  inline const Real_wp& VarProblem_Base::GetOmega() const
  {
    return this->omega;
  }

  
  //! returns - i omega for complex numbers, 1 for real numbers
  inline void VarProblem_Base::GetMiomega(Real_wp& k) const
  {
    k = 1.0;
  }


  //! returns - i omega for complex numbers, 1 for real numbers
  inline void VarProblem_Base::GetMiomega(Complex_wp& k) const
  {
    k = -Iwp*this->omega;
  }


  //! returns - omega^2 for complex numbers, 1 for real numbers
  inline void VarProblem_Base::GetMomega2(Real_wp& k) const
  {
    k = 1.0;
  }


  //! returns - omega^2 for complex numbers, 1 for real numbers
  inline void VarProblem_Base::GetMomega2(Complex_wp& k) const
  {
    k = -this->omega2;
  }

  
  //! returns the frequency omega/ (2 pi)
  inline const Real_wp& VarProblem_Base::GetFrequency() const
  {
    return this->frequency;
  }
 

  //! returns the characteristical length (used for adimensionalization)
  inline const Real_wp& VarProblem_Base::GetWaveLengthAdim() const
  {
    return wavelength_adim;
  }
  
  
  
  /******************************
   * VarComputationProblem_Base *
   ******************************/


  //! Destructor
  inline VarComputationProblem_Base::~VarComputationProblem_Base()
  {
  }  


  //! returns the threshold used to drop values in finite element matrix
  inline Real_wp VarComputationProblem_Base::GetThresholdMatrix() const
  {
    return threshold_matrix;
  }


  //! sets the threshold used to drop values in finite element matrix
  inline void VarComputationProblem_Base::SetThresholdMatrix(const Real_wp& epsilon)
  {
    threshold_matrix = epsilon;
  }


  //! Fills nodes to be eliminated (static condensation)
  /*!
    \param[in] i element number
    \param[in] nb_dof_loc number of degrees of freedom for element i
    \param[out] nb_dof_edges number of degrees of freedom that cannot be eliminated
    \param[out] nb_dof_int number of degrees of freedom that can be eliminated
    \param[out] intern_node numbering of degrees freedom
    if intern_node(i) is positive, it is the local number of dofs that are kept
    if intern_node(i) is negative, it is -local number-1 of dofs that are eliminated
   */
  inline void VarComputationProblem_Base::
  GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
                          Vector<int>& intern_node) const
  {
    // no static condensation
    nb_dof_edges = nb_dof_loc;
    nb_dof_int = 0;
    intern_node.Reallocate(nb_dof_loc); intern_node.Fill();
  }

  
  //! returns a new condensed block solver
  inline CondensationBlockSolver_Base<Real_wp>* VarComputationProblem_Base
  ::GetNewCondensationSolver(const Real_wp&)
  {
    return NULL;
  }


  //! returns a new condensed block solver
  inline CondensationBlockSolver_Base<Complex_wp>* VarComputationProblem_Base
  ::GetNewCondensationSolver(const Complex_wp&)
  {
    return NULL;
  }
  

  /*************************
   * VarComputationProblem *
   *************************/


  //! default constructor
  template<class TypeEquation>
  inline VarComputationProblem::VarComputationProblem(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_boundary(var), var_comm(var)
  {
    InitDefaultValues();

    this->symmetric_elementary_matrix = TypeEquation::SymmetricElementaryMatrix();
    this->symmetric_global_matrix = TypeEquation::SymmetricGlobalMatrix();    
  }

    
  //! returns true if we should not store the matrix, and use a matrix-free algorithm
  inline bool VarComputationProblem::UseMatrixFreeAlgorithm() const
  {
    return (storage_finite_element_matrix == MATRIX_FREE);
  }


  //! returns true if the global matrix is symmetric
  inline bool VarComputationProblem::IsSymmetricGlobalMatrix() const
  {
    return symmetric_global_matrix;
  }


  //! returns the storage used for the finite element matrix
  inline int VarComputationProblem::GetStorageFiniteElementMatrix() const
  {
    return storage_finite_element_matrix;
  }
  
  
  //! sets which storage to use for the finite element matrix
  inline void VarComputationProblem::SetStorageFiniteElementMatrix(int type)
  {
    storage_finite_element_matrix = type;
  }


  //! informs that the elementary matrix is symmetric
  inline void VarComputationProblem::SetSymmetricElementaryMatrix(bool sym)
  {
    symmetric_elementary_matrix = sym;
  }


  //! temporary sets a static condensation to compute a condensed matrix
  inline void VarComputationProblem::SetLeafStaticCondensation(bool leaf) const
  {
    leaf_static_condensation = leaf;
  }


  //! returns true if the condensed matrix must be computed
  inline bool VarComputationProblem::GetLeafStaticCondensation() const
  {
    return leaf_static_condensation;
  }


  //! returns true if the static condensation is light (vectorial unknowns are eliminated)
  inline bool VarComputationProblem::LightStaticCondensation() const
  {
    return light_static_condensation;
  }


  //! returns the number of degrees of freedom of the problem
  inline int VarComputationProblem::GetNbRows() const
  {
    return var_problem.GetNbDof();
  }


  //! returns true if the matrix will by symmetrized
  inline bool VarComputationProblem::GetSymmetrizationUse() const
  {
    return use_symmetrization_when_possible;
  }


  //! sets if the matrix should be symmetrized
  inline void VarComputationProblem::SetSymmetrizationUse(bool sym)
  {
    use_symmetrization_when_possible = sym;
  }

  
  //! sets if homogeneous Dirichlet is used
  inline void VarComputationProblem::SetHomogeneousDirichlet(bool hg_dir)
  {
    // for homogeneous Dirichlet, no need to store columns of Dirichlet dofs
    erase_dirichlet_columns = hg_dir;
  }
  

  //! returns true if Dirichlet condition is homogeneous
  inline bool VarComputationProblem::IsHomogeneousDirichlet() const
  {
    return erase_dirichlet_columns;
  }


  //! returns the print level
  inline int VarComputationProblem::GetPrintLevel() const
  {
    return var_problem.print_level;
  }

  
  //! sets the verbosity level
  inline void VarComputationProblem::SetPrintLevel(int lvl)
  {
    var_problem.print_level = lvl;
  }
  

  //! returns the memory used by the object in bytes
  inline size_t VarComputationProblem::GetMemorySize() const
  {
    return 0;
  }
  
}
  
#define MONTJOIE_FILE_VAR_PROBLEM_BASE_INLINE_CXX
#endif
