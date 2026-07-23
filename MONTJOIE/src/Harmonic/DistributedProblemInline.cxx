#ifndef MONTJOIE_FILE_DISTRIBUTED_PROBLEM_INLINE_CXX

namespace Montjoie
{

  //! constructor with a given problem
  template<class TypeEquation>
  inline DistributedProblem_Base::DistributedProblem_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var), var_boundary(var)
  {
    InitDefaultValues();
  }
  

  //! returns a reference to the object VarComputationProblem
  inline VarProblem_Base& DistributedProblem_Base::GetVarProblemBase()
  {
    return var_problem;
  }


  //! returns a reference to the object VarComputationProblem
  inline const VarProblem_Base& DistributedProblem_Base::GetVarProblemBase() const
  {
    return var_problem;
  }

  
  //! Destructor
  inline DistributedProblem_Base::~DistributedProblem_Base()
  {
  }
  
  
  //! returns the global number of dofs contained in the mesh
  /*!
    This number of dofs is equal to the number of dofs in the mesh
    as if the simulation was completed in sequential
    \param[in] n unknown number
  */
  inline int DistributedProblem_Base::GetNbGlobalMeshDof(int n) const
  {
    return nodl_mesh(n);
  }


  //! returns the global number of dofs for an unknown in PML layers
  /*!
    This number of dofs is equal to the number of dofs in the PML layers
    as if the simulation was completed in sequential
  */
  inline int DistributedProblem_Base::GetNbGlobalDofPML(int n) const
  {
    return nodl_pml(n);
  }
  
  
  //! returns the global number of dofs
  /*!
    This number of dofs is equal to the size of the finite element matrix
    as if the simulation was completed in sequential
  */
  inline int DistributedProblem_Base::GetNbGlobalDof() const
  {
    return nodl_all;
  }


  //! returns the offset dof unknown (global dofs)
  inline int DistributedProblem_Base::GetOffsetGlobalUnknownDof(int n) const
  {
    return offset_dof_unknown_all(n);
  }

  
  //! returns the global number of dofs for a given unknown
  inline int DistributedProblem_Base::GetNbGlobalUnknownDof(int n) const
  {
    return offset_dof_unknown_all(n+1) - offset_dof_unknown_all(n);
  }
  
  
  //! returns the global number for an unknown (with static condensation)
  inline int DistributedProblem_Base::GetNbGlobalCondensedDof(int n) const
  {
    return offset_dof_condensed_all(n+1) - offset_dof_condensed_all(n);
  }
  
  
  //! returns the global dof number of local dof i
  inline int DistributedProblem_Base::GetGlobalDofNumber(int i) const
  {
    return GlobDofNumber_Subdomain(i);
  }


  //! returns the global dof numbers
  inline const IVect& DistributedProblem_Base::GetGlobalDofNumber() const
  {
    return GlobDofNumber_Subdomain;
  }
  

  //! returns the number of quadrature points on boundary condition NEIGHBOR
  inline int DistributedProblem_Base::GetNbPointsQuadratureNeighbor() const
  {
    return this->nb_points_quadrature_on_neighboring_interfaces;
  }
  
  
  //! returns the number of subdomain in interaction with the current one
  inline int DistributedProblem_Base::GetNbSubdomains() const
  {
    return MatchingDofOrig_Subdomain.GetM();
  }
  
  
  //! returns the number of processors involved for the computation of one mode
  inline int DistributedProblem_Base::GetNbProcPerMode() const
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(comm_group_mode, &nb_proc);
    return nb_proc;
#endif
    
    return 1;
  }


  //! returns the number of processors involved for the computation of one mode
  inline int DistributedProblem_Base::GetRankProcMode() const
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(comm_group_mode, &rank_proc);
    return rank_proc;
#endif
    
    return 0;
  }
  

  //! gives access to MatchingDofOrig_Subdomain
  inline Vector<IVect>& DistributedProblem_Base::GetOriginalMatchingDofNeighbor()
  {
    return MatchingDofOrig_Subdomain;
  }

  
  //! gives access to MatchingNumber_Subdomain
  inline IVect& DistributedProblem_Base::GetProcMatchingNeighbor()
  {
    return MatchingNumber_Subdomain;
  }
  

  //! returns the global element number adjacent to the neighboring face num_face
  inline int DistributedProblem_Base::GetElementNumberNeighboringFace(int num_face) const
  {
    return NeighboringConnectivity(num_face)(1);
  }


  //! returns the offset for quadrature points of the neighboring face num_face
  inline int DistributedProblem_Base::GetOffsetNeighboringFace(int num_face) const
  {
    return NeighboringConnectivity(num_face)(9);
  }
  

  //! returns the processor that owns the adjacent element of the neighboring face num_face
  inline int DistributedProblem_Base::GetProcessorNeighboringFace(int nf) const
  {
    return NeighboringConnectivity(nf)(2);
  }
   

  //! returns the local position of the face on the adjacent element (of the neighboring face nf)
  inline int DistributedProblem_Base::GetLocalPositionNeighboringFace(int nf) const
  {
    return NeighboringConnectivity(nf)(4); 
  }
  

  //! returns the rotation of the face between the two elements
  inline int DistributedProblem_Base::GetRotationNeighboringFace(int nf) const
  {
    return NeighboringConnectivity(nf)(5);
  }
  

  //! returns the order of the adjacent element of the neighboring face nf
  inline int DistributedProblem_Base::GetOrderEltNeighboringFace(int nf, int n) const
  {
    int p = NeighboringConnectivity(nf)(10 + n);
    return NeighboringConnectivity(nf)(p);
  }
  

  //! returns the hybrid type of the adjacent element (0, 1, 3 or 3) of the neighboring face nf
  inline int DistributedProblem_Base::GetTypeEltNeighboringFace(int nf) const
  {
    return NeighboringConnectivity(nf)(6);
  }
    

  //! returns the reference of the adjacent element of the neighboring face nf
  inline int DistributedProblem_Base::GetRefDomainNeighboringFace(int nf) const
  {
    return NeighboringConnectivity(nf)(7);
  }
    

  //! returns the size of array OffsetDofV (usually nb_elt + 1)
  inline int DistributedProblem_Base::GetSizeOffsetDofV() const
  {
    return OffsetDofV.GetM();
  }

  
  //! returns the offset for vectorial unknowns (element i)
  inline int DistributedProblem_Base::GetOffsetDofV(int i) const
  {
    return OffsetDofV(i);
  }


  //! sets the offset for vectorial unknowns (element i)
  inline void DistributedProblem_Base::SetOffsetDofV(int i, int num)
  {
    OffsetDofV(i) = num;
  }
    

  //! returns the number of dofs that belong to another processor
  inline int DistributedProblem_Base::GetNbOverlappedDof() const
  {
    return OverlapDofNumber_Subdomain.GetM();
  }
  

  //! returns the i-th dof number that belong to another processor
  inline int DistributedProblem_Base::GetOverlappedDofNumber(int i) const
  {
    return OverlapDofNumber_Subdomain(i);
  }


  //! returns the list of dof numbers that belong to another processor
  inline const IVect& DistributedProblem_Base::GetOverlappedDofNumber() const
  {
    return OverlapDofNumber_Subdomain;
  }


  //! returns the processor that owns overlapped dof i
  inline int DistributedProblem_Base::GetOverlappedProcNumber(int i) const
  {
    return OverlapProcNumber_Subdomain(i);
  }


  //! returns the list of processors that owns overlapped dofs
  inline const IVect& DistributedProblem_Base::GetOverlappedProcNumber() const
  {
    return OverlapProcNumber_Subdomain;
  }


#ifdef SELDON_WITH_MPI
  //! returns a new distributed vector with the correct references
  template<class T>
  inline DistributedVector<T>* DistributedProblem_Base
  ::AllocateDistributedVector(Vector<T>& x) const
  {
    DistributedVector<T>* d;
    d = new DistributedVector<T>(this->OverlapDofNumber_Subdomain, this->comm_group_mode);  
    d->SetData(x.GetM(), x.GetData());
    return d;
  }


  //! Nullifies x and deletes the pointer
  template<class T>
  inline void DistributedProblem_Base::NullifyDistributedVector(DistributedVector<T>* x) const
  {
    x->Nullify();
    delete x;
  }
#endif


  //! sets the array epart used to split the mesh between subdomains
  inline void DistributedProblem_Base::SetEpartSplitting(const IVect& epart)
  {
    splitting_algorithm = SUBDIV_DOMAIN_USER;
    split_mesh_epart_user = epart; 
  }
    

  //! sets the output file name for epart array
  inline void DistributedProblem_Base::SaveEpartSplitting(const string& name)
  {
    partition_subdomain_to_store = true;
    file_name_partition_subdomain = name;
  }
  

  /**********************
   * DistributedProblem *
   **********************/
  
  
  //! constructor with a given problem
  template<class Dimension> template<class TypeEquation>
  inline DistributedProblem<Dimension>::DistributedProblem(EllipticProblem<TypeEquation>& var)
    : DistributedProblem_Base(var), VarGeometryProblem<Dimension>(var),
      var_computation(var), var_boundary(var), var_output(var)
  {
  }

  
  //! returns a reference to the object VarComputationProblem
  template<class Dimension>
  inline VarComputationProblem& DistributedProblem<Dimension>::GetComputationProblem()
  {
    return var_computation;
  }


  //! returns a reference to the object VarComputationProblem
  template<class Dimension>
  inline const VarComputationProblem& DistributedProblem<Dimension>::GetComputationProblem() const
  {
    return var_computation;
  }

  
  //! returns a reference to the object VarBoundaryCondition_Base
  template<class Dimension>
  inline VarBoundaryCondition_Base& DistributedProblem<Dimension>::GetBoundaryConditionProblem()
  {
    return var_boundary;
  }


  //! returns a reference to the object VarOutputProblem_Dim
  template<class Dimension>
  inline VarOutputProblem_Dim<Dimension>& DistributedProblem<Dimension>::GetOutputProblem()
  {
    return var_output;
  }


  //! returns a reference to the object VarBoundaryCondition_Base
  template<class Dimension>
  inline const VarBoundaryCondition_Base& DistributedProblem<Dimension>
  ::GetBoundaryConditionProblem() const
  {
    return var_boundary;
  }

  
  template<class Dimension>
  inline int DistributedProblem<Dimension>::GetNbMainUnknownDof() const
  {
    if (this->FormulationDG() == ElementReference_Base::HDG)
      return this->OffsetDofV(this->mesh.GetNbElt());

    return this->GetNbDofUnknown(0);
  }
  
  
  //! Computes E \times n and H \times n on nodal points of a face
  /*!
    For Helmholtz equation, it would be u and du/dn
    the definition depends on the equation, that's why it is defined in derive classes
    
    \param[in] u_nodal solution on nodal points
    \param[in] grad_nodal gradient of solution on nodal points
    \param[in] num_elem element number
    \param[in] pts nodal points
    \param[in] normale normales associated with nodal points
    \param[out] En_nodal values of E \times n on nodal points
    \param[out] Hn_nodal values of H \times n on nodal points
    \warning ModifyOutputUnknown has been called to obtain u_nodal and grad_nodal
    
    u_nodal(m)(j) is a component of the solution (m is the component number,
    j the number of the nodal point), grad_nodal(m)(j) a component of the gradient
    (or curl or divergence depending on the type of the finite element)
    
    En_nodal(:)(j) will be the vector E \times j for the nodal point j
   */
  template<class Dimension>
  inline void DistributedProblem<Dimension>
  ::ComputeEnHnNodal(Vector<VectReal_wp>& u_nodal, Vector<VectReal_wp>& grad_nodal,
		     int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
		     Vector<VectReal_wp>& En_nodal, Vector<VectReal_wp>& Hn_nodal) const
  {
    cout << "This method is not available for this equation" << endl;
    abort();
  }
    
  
  template<class Dimension>
  inline void DistributedProblem<Dimension>
  ::ComputeEnHnNodal(Vector<VectComplex_wp>& u_nodal, Vector<VectComplex_wp>& grad_nodal,
		     int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
		     Vector<VectComplex_wp>& En_nodal, Vector<VectComplex_wp>& Hn_nodal) const 
  {
    cout << "This method is not available for this equation" << endl;
    abort();
  }
  
}

#define MONTJOIE_FILE_DISTRIBUTED_PROBLEM_INLINE_CXX
#endif

