#ifndef MONTJOIE_FILE_DISTRIBUTED_PROBLEM_HXX

namespace Montjoie
{
  
  //! base class for distributed problems
  class DistributedProblem_Base
  {
    friend class VarBoundaryCondition_Base;
    template<class Dim> friend class VarGeometryProblem;
    
  private :
    //! reference to the object var_problem
    VarProblem_Base& var_problem;
    VarBoundaryCondition_Base& var_boundary;
    
    void InitDefaultValues();
  protected:

    //! adjacent subdomains
    /*!
      This arrays contains the list of processors interacting with the current one
    */
    IVect MatchingNumber_Subdomain; 
    

    //! dof numbers of the current subdomain on the interface
    /*!
      MatchingDofOrig_Subdomain(i) contains all the dofs common with processor i
      these dofs are sorted so that the dofs contained in the corresponding array of processor i
      should match
    */
    Vector<IVect> MatchingDofOrig_Subdomain;

    //! for Dirichlet dofs only
    IVect DirichletNumber_Subdomain;
    Vector<IVect> DirichletDofOrig_Subdomain;
        
  public :    
#ifdef SELDON_WITH_MPI
    //! communicator for processors sharing the mesh
    /*!
      Usually the mesh is split into several parts, each part being associated with a processor
      this communicator regroups all the processor on which the mesh is distributed
    */
    MPI_Comm comm_group_mode;
#endif
    
    // available partitioning algorithms
    enum {SUBDIV_DOMAIN_SCOTCH, SUBDIV_DOMAIN_METIS, SUBDIV_DOMAIN_USER,
          SUBDIV_DOMAIN_CONCENTRIC, SUBDIV_DOMAIN_LAYERED};

  protected :
    //! which algorithm is used to split domains in small sub-domains ?
    int splitting_algorithm;
    int coef_split_default_domain, coef_split_default_pml;
        
    //! informations about neighbor elements for faces shared with other processors
    /*!
      0 : numero global de l'element
      1 : numero global de l'element en face
      2 : processeur de l'element en face
      3 : numero local de l'element en face (pas utilise)
      4 : position locale du bord dans l'element en faec
      5 : rotation entre les deux faces
      6 : ordre element en face
      7 : type element en face
      8 : nombre de ddl de l'element en face
      9 : offset du nombre de points de quadrature
      10 : nombre de ddl de l'element en face (pour PML)
      11 : reference de l'element en face
      12 -> 12+N : numero des ddls de l'element en face
      12+N -> 12 +2*N : numero des ddls PML de l'element en face
    */
    Vector<IVect> NeighboringConnectivity;

    //! quadrature points numbers of the current subdomain on the interface
    /*!
      MatchingQuadOrig_Subdomain(i) contains all the quadrature points common with processor i
      these points are sorted so that the points contained
      in the corresponding array of processor i
      should match, i.e. quadrature points should have same coordinates
    */
    Vector<IVect> MatchingQuadOrig_Subdomain;
    
    //! total number of quadrature points common with other processors
    int nb_points_quadrature_on_neighboring_interfaces;
    //! processors sharing the mesh
    IVect list_processor_on_group;
    
    //! dofs overlapping with other domains
    /*!
      This arrays lists all the dofs which are already contained on a previous processor
      The first processor to contain a dof is called "original processor"
    */
    IVect OverlapDofNumber_Subdomain;
    //! original processor where overlapped dofs are
    /*!
      For each dof shared with other processors,
      this array contains the number of the first processor
    */
    IVect OverlapProcNumber_Subdomain;
    
    Vector<int> nodl_mesh; //! total number of dofs associated with the meshes
    int nodl_all; //! total number of dofs (if all subdomains are joined)
    Vector<int> nodl_pml; //! total number of dofs in PML layers

    Vector<int> offset_dof_unknown_all, offset_dof_condensed_all;
    
    //! offsets for vectorial dofs (L2 unknowns)
    IVect OffsetDofV;
    //! global offsets for vectorial dofs (L2 unknowns)
    IVect GlobOffsetDofV_Subdomain;
    //! global row numbers
    IVect GlobDofNumber_Subdomain;
    
    // array "Epart" supplied by the user for SUBDIV_DOMAIN_USER
    IVect split_mesh_epart_user;
    bool partition_subdomain_to_store;
    string file_name_partition_subdomain;
    VectReal_wp radius_splitting_subdomain;
    Vector<IVect> ref_layered_subdivision;

  public :
    template<class TypeEquation>
    DistributedProblem_Base(EllipticProblem<TypeEquation>& var);

    virtual ~DistributedProblem_Base();
    virtual int GetNbMainUnknownDof() const = 0;
    virtual IVect GetDofNumberOnElement(int, int num=0) const = 0;
    virtual IVect GetScalarDofNumberOnElement(int) const = 0;
    
    // inline functions
    VarProblem_Base& GetVarProblemBase();
    const VarProblem_Base& GetVarProblemBase() const;
    
    int GetNbGlobalMeshDof(int n = 0) const;
    int GetNbGlobalDofPML(int n = 0) const;
    int GetNbGlobalDof() const;

    int GetOffsetGlobalUnknownDof(int n) const;
    int GetNbGlobalUnknownDof(int n = 0) const;
    int GetNbGlobalCondensedDof(int n = 0) const;
    
    int GetGlobalDofNumber(int i) const;
    const IVect& GetGlobalDofNumber() const;
    int GetNbPointsQuadratureNeighbor() const;
    int GetNbSubdomains() const;
    int GetNbProcPerMode() const;
    int GetRankProcMode() const;

    IVect& GetProcMatchingNeighbor();
    Vector<IVect>& GetOriginalMatchingDofNeighbor();
    int GetElementNumberNeighboringFace(int num_face) const;
    int GetOffsetNeighboringFace(int num_face) const;
    int GetProcessorNeighboringFace(int nf) const;
    int GetLocalPositionNeighboringFace(int nf) const;
    int GetRotationNeighboringFace(int nf) const;
    int GetOrderEltNeighboringFace(int nf, int n = 0) const;
    int GetTypeEltNeighboringFace(int nf) const;
    IVect GetNodleNeighboringFace(int nf, int n = 0) const;
    IVect GetNodlePmlNeighboringFace(int nf, int n = 0) const;

    int GetRefDomainNeighboringFace(int nf) const;
    
    int GetSizeOffsetDofV() const;
    int GetOffsetDofV(int i) const;
    void SetOffsetDofV(int i, int num);
    int GetNbOverlappedDof() const;
    int GetOverlappedDofNumber(int i) const;
    const IVect& GetOverlappedDofNumber() const;
    int GetOverlappedProcNumber(int i) const;
    const IVect& GetOverlappedProcNumber() const;    

#ifdef SELDON_WITH_MPI
    template<class T>
    DistributedVector<T>* AllocateDistributedVector(Vector<T>& x) const;

    template<class T>
    void NullifyDistributedVector(DistributedVector<T>* x) const;
#endif

    void SetEpartSplitting(const IVect& epart);
    void SaveEpartSplitting(const string& name);
    
    // other functions    
    virtual void GetMemoryUsed(map<string, size_t>& var) const;
    void DisplayMemoryUsed(map<string, size_t>& var) const;

    void SetInputData(const string& description_field,
                      const VectString& parameters);
    
    // purely virtual methods
    virtual void FindElementsInsidePML() = 0;
    virtual IVect GetNodleElement(int i, int n = 0) const = 0;
    
    // MPI methods 
    template<class Matrix1>
    void InitDistributedMatrix(Matrix1& mat);
    
    template<class T, class Prop, class Storage, class Allocator>
    void InitDistributedMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& mat);

    template<class T>
    void AddDomains(Vector<T>& X, int nb_u = -1) const;

    template<class T>
    void AssembleDirichlet(Vector<T>& X, int nb_u = -1, bool only_num = false) const;

#ifdef SELDON_WITH_MPI
    template<class T>
    void ReduceDirichlet(Vector<T>& X, const MPI_Op& oper, int nb_u = -1, bool only_num = false) const;
#endif

    template<class T>
    void ExchangeDomains(Vector<T>& X, int nb_u = -1) const;

    template<class T>
    void ExchangeRelaxDomains(Vector<T>& X, const Real_wp& omega,
			      int proc, int nb_u = -1) const;

    template<class T>
    void ExchangeQuadRelaxDomains(Vector<T>& X, const Real_wp& omega,
				  int proc, int nb_u = -1) const;
    
#ifdef SELDON_WITH_MPI
    void ConstructDirichletComm();
    
    template<class T>
    void ReduceDistributedVector(Vector<T>&, const MPI_Op& op,
                                 int nb_u = -1, bool only_num = false) const;
    
    template<class T>
    void ReduceDistributedVector(Vector<T>& X, const MPI_Op& oper, const IVect& ProcNumber,
				 const Vector<IVect>& DofNumber, int Nvol, int nb_u) const;

    template<class T>
    void ExchangeRelaxDistributedVector(Vector<T>& X, const Real_wp& omega,
					int proc, const IVect& ProcNumber,
					const Vector<IVect>& DofNumber, int Nvol, int nb_u) const;
    
    template<class T>
    void ExchangeDistributedVector(Vector<T>& X, const IVect& ProcNumber,
				   const Vector<IVect>& DofNumber, int Nvol, int nb_u) const;
    
    template<class T>
    void ExchangeUfaceDomains(const Vector<T>& X, Vector<Vector<T> >&, Vector<Vector<int64_t> >&,
                              Vector<Vector<T> >&, Vector<Vector<int64_t> >&,
                              Vector<MPI_Request>& request, int tag) const;
    
    template<class T>
    void GetUfaceDomains(Vector<T>& X, Vector<Vector<T> >& xdom, Vector<Vector<int64_t> >&,
                         Vector<Vector<T> >& xdom2, Vector<Vector<int64_t> >&,
                         Vector<MPI_Request>& request, int tag) const;    
#endif
    
    void CopyInputData(const DistributedProblem_Base& var);
    
  };
  
  
  //! class handling problem distributed on several processors
  /*!
    In this class, we put methods specific for parallel implementation
    we use the term "subdomain" to name the part of mesh affected on each processor
  */
  template<class Dimension>
  class DistributedProblem : public DistributedProblem_Base, public VarGeometryProblem<Dimension>    
  {
  private:
    VarComputationProblem& var_computation;
    VarBoundaryCondition_Base& var_boundary;
    VarOutputProblem_Dim<Dimension>& var_output;
    
  public :    
    typedef typename Dimension::R_N R_N;

    template<class TypeEquation>
    DistributedProblem(EllipticProblem<TypeEquation>& var);

    // Inline methods
    VarComputationProblem& GetComputationProblem();
    const VarComputationProblem& GetComputationProblem() const;
    VarBoundaryCondition_Base& GetBoundaryConditionProblem();
    VarOutputProblem_Dim<Dimension>& GetOutputProblem();
    const VarBoundaryCondition_Base& GetBoundaryConditionProblem() const;

    int GetNbMainUnknownDof() const;
    
    // other methods
    void SetInputData(const string& description_field, const VectString& parameters);
    IVect GetDofNumberOnElement(int, int num = 0) const;
    IVect GetScalarDofNumberOnElement(int) const;
    
    void GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
				 Vector<int>& intern_node) const;

    void GetStaticCondensedRows(IVect& IndexCondensedRows, IVect&, IVect&, IVect&,
				int& nb_scalar_dof, int&,
				IVect&, Vector<IVect>&) const;
    
    virtual void GetMemoryUsed(map<string, size_t>& var) const;    
    
    // main methods for distribute the mesh into several processors
    void SplitMeshForParallelComputation(const string& name_elt);    
    void SplitSubdomains(const string& name_elt);
    
#ifdef SELDON_WITH_MPI    
    template<class T>
    void ReduceDistributedVectorFace(Vector<T>& X, const MPI_Op& oper, int nb_u = 1) const;
#endif

    // Inline virtual methods
    virtual void ComputeLocalProlongation(FiniteElementInterpolator& proj,
					  DistributedProblem<Dimension>& var_coarse, int rc, 
					  TinyVector<Matrix<Real_wp>, 4>& ProlongationElement) const = 0;    
    
    virtual DistributedProblem<Dimension>* GetNewEllipticProblem() const = 0;
    virtual void CopyFiniteElement(const DistributedProblem<Dimension>& var) = 0;
    
    virtual void ComputeEnHnNodal(Vector<VectReal_wp>& u_nodal, Vector<VectReal_wp>& grad_nodal,
				  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
				  Vector<VectReal_wp>& En_nodal, Vector<VectReal_wp>& Hn_nodal) const;
    
    virtual void ComputeEnHnNodal(Vector<VectComplex_wp>& u_nodal, Vector<VectComplex_wp>& grad_nodal,
				  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
				  Vector<VectComplex_wp>& En_nodal, Vector<VectComplex_wp>& Hn_nodal) const;

    // other virtual methods
    bool IsSymmetricProblem(bool eigen = false) const;    
    bool IsSymmetricMassMatrix() const;
    void CopyInputData(const VarProblem_Base& var_base);
    IVect GetNodleElement(int i, int n = 0) const;
    
  };
  
}

#define MONTJOIE_FILE_DISTRIBUTED_PROBLEM_HXX
#endif

