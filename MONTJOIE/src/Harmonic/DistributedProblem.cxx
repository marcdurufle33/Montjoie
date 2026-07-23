#ifndef MONTJOIE_FILE_DISTRIBUTED_PROBLEM_CXX

namespace Montjoie
{
  
  /***************************
   * DistributedProblem_Base *
   ***************************/
  
  
  //! Default constructor
  void DistributedProblem_Base::InitDefaultValues()
  {    
#ifdef SELDON_WITH_MPI 
    comm_group_mode = MPI_COMM_WORLD;
    int nb_proc; MPI_Comm_size(comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif

    if (nb_proc == 1)
      {
	// running simulation on a single processor
	// no splitting of the mesh
	splitting_algorithm = -1;
      }
    else
      {
	// default splitting algorithm
	splitting_algorithm = SUBDIV_DOMAIN_METIS;

#ifdef MONTJOIE_WITH_SCOTCH
        // splitting_algorithm = SUBDIV_DOMAIN_SCOTCH;
#endif        
      }

    coef_split_default_domain = 2;
    coef_split_default_pml = 3;

    nb_points_quadrature_on_neighboring_interfaces = 0;
    
    nodl_all = 0;
    
    partition_subdomain_to_store = false;
  }
    
  
  //! adds memory used by the object
  void DistributedProblem_Base
  ::GetMemoryUsed(map<string, size_t>& var) const
  {
    size_t taille = sizeof(*this) + MatchingNumber_Subdomain.GetMemorySize()
      + Seldon::GetMemorySize(MatchingDofOrig_Subdomain) + Seldon::GetMemorySize(NeighboringConnectivity) 
      + Seldon::GetMemorySize(MatchingQuadOrig_Subdomain) + list_processor_on_group.GetMemorySize()
      + OverlapDofNumber_Subdomain.GetMemorySize() + OverlapProcNumber_Subdomain.GetMemorySize();
    
    taille += OffsetDofV.GetMemorySize() + GlobOffsetDofV_Subdomain.GetMemorySize()
      + split_mesh_epart_user.GetMemorySize() + radius_splitting_subdomain.GetMemorySize()
      + Seldon::GetMemorySize(ref_layered_subdivision);
    
    var["Parallel Number"] = taille;
  }
  

  //! displays memory used
  void DistributedProblem_Base::DisplayMemoryUsed(map<string, size_t>& var) const
  {
    // storing size of each element in the vector size_elt
    Vector<size_t> size_elt(var.size());
    int nb = 0;
    map<string, size_t>::iterator it;
    for (it = var.begin(); it != var.end(); ++it)
      size_elt(nb++) = it->second;

    // then suming contributions of other processors
#ifdef SELDON_WITH_MPI
    Vector<size_t> size_target(var.size());
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
    MPI_Reduce(size_elt.GetData(), size_target.GetData(), size_target.GetM(),
               MPI_UNSIGNED_LONG, MPI_SUM, 0, this->comm_group_mode);
    
    size_elt = size_target;
    if (rank_proc == 0)
      {
#endif        
        // displaying memory taken by each element
        size_t overall_size = 0;
        string size_chaine; nb = 0;
        for (it = var.begin(); it != var.end(); ++it)
          {
            size_t taille = size_elt(nb++);
            size_chaine = GetHumanReadableMemory(taille);
            overall_size += taille;
            cout << "Size taken by " << it->first << " : " << size_chaine << endl;
          }
        
        // memory taken by all the elements
        size_chaine = GetHumanReadableMemory(overall_size);
        cout << "Size taken by overall simulation : " << size_chaine << endl;

#ifdef SELDON_WITH_MPI
      }
#endif
    
  }

  
  //! returns the global dof numbers of the adjacent element of the neighboring face nf
  IVect DistributedProblem_Base::GetNodleNeighboringFace(int nf, int n) const
  {
    int p = NeighboringConnectivity(nf)(10 + n);
    int nb_ddl = NeighboringConnectivity(nf)(p+1);
    IVect tmp(nb_ddl);
    for (int i = 0; i < nb_ddl; i++)
      tmp(i) = NeighboringConnectivity(nf)(p+3+i);
    
    return tmp; 
  }
  

  //! returns the global pml dof numbers of the adjacent element of the neighboring face nf
  IVect DistributedProblem_Base::GetNodlePmlNeighboringFace(int nf, int n) const
  {
    int p = NeighboringConnectivity(nf)(10 + n);
    int nb_ddl = NeighboringConnectivity(nf)(p+2);
    IVect tmp;
    if (nb_ddl > 0)
      {
        int offset = p + 3 + NeighboringConnectivity(nf)(p+1);
	tmp.Reallocate(nb_ddl);
	for (int i = 0; i < nb_ddl; i++)
	  tmp(i) = NeighboringConnectivity(nf)(offset + i);
      }
    
    return tmp; 
  }

  
  //! setting parallel parameters from a data line
  void DistributedProblem_Base::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("NbProcessorsPerMode"))
      {
        
#ifdef SELDON_WITH_MPI
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of DistributedProblem" << endl;
	    cout << "NbProcessorsPerMode needs more parameters, for instance :" << endl;
	    cout << "NbProcessorsPerMode = num" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        int nb_proc = to_num<int>(parameters(0));
        if (nb_proc == 1)
          {
            comm_group_mode = MPI_COMM_SELF;
          }
        else
          {
            cout << "Not implemented" << endl;
            abort();
          }        
#endif
      }
    else if (!description_field.compare("SaveSplitDomain"))
      {
	if (parameters(0) == "YES")
	  {
	    this->partition_subdomain_to_store = true;
	    this->file_name_partition_subdomain = parameters(1);
	  }
      }
    else if (!description_field.compare("SplitDomain"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of DistributedProblem" << endl;
	    cout << "SplitDomain needs more parameters, for instance :" << endl;
	    cout << "SplitDomain = Scotch" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

#ifdef SELDON_WITH_MPI 
        int nb_proc; MPI_Comm_size(comm_group_mode, &nb_proc);
        int rank_proc; MPI_Comm_rank(comm_group_mode, &rank_proc);
#else
        int nb_proc(1), rank_proc(0);
#endif
        
        int nb = 1;
	if (!parameters(0).compare("Scotch"))
          this->splitting_algorithm = SUBDIV_DOMAIN_SCOTCH;
	else if (!parameters(0).compare("User"))
          {
	    if (parameters.GetM() <= 1)
	      {
		cout << "In SetInputData of DistributedProblem" << endl;
		cout << "SplitDomain needs more parameters, for instance :" << endl;
		cout << "SplitDomain = User file_name" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

            this->splitting_algorithm = SUBDIV_DOMAIN_USER;
            if (rank_proc == 0)
	      this->split_mesh_epart_user.Read(parameters(1));
	    
            nb = 2;
          }
        else if (!parameters(0).compare("Concentric"))
          {
            this->splitting_algorithm = SUBDIV_DOMAIN_CONCENTRIC;
            radius_splitting_subdomain.Reallocate(nb_proc-1);
            if (parameters(1) == "AUTO")
              {
                abort();
                nb = 2;
              }
            else
              {
                for (int i = 0; i < nb_proc-1; i++)
                  radius_splitting_subdomain(i) = to_num<Real_wp>(parameters(i+1));
                
                nb = nb_proc;
              }            
          }
        else if (!parameters(0).compare("Layered"))
          {
            this->splitting_algorithm = SUBDIV_DOMAIN_LAYERED;
            if (parameters(1) == "AUTO")
              {
                ref_layered_subdivision.Reallocate(nb_proc);
                for (int i = 0; i < ref_layered_subdivision.GetM(); i++)
                  {
                    ref_layered_subdivision(i).Reallocate(1);
                    ref_layered_subdivision(i)(0) = i+1;  
                  }
              }
            else
              {
                ref_layered_subdivision.Reallocate(nb_proc);
                int nb = 1;
                for (int i = 0; i < ref_layered_subdivision.GetM(); i++)
                  {
                    int nb_ref = to_num<int>(parameters(nb++));
                    ref_layered_subdivision(i).Reallocate(nb_ref);
                    for (int j = 0; j < nb_ref; j++)
                      ref_layered_subdivision(i)(j) = to_num<int>(parameters(nb++));
                  }
              }
          }
        else
          this->splitting_algorithm = SUBDIV_DOMAIN_METIS;
        
        if (parameters.GetM() > nb)
	  coef_split_default_domain = to_num<int>(parameters(nb++));
	
	if (parameters.GetM() > nb)
	  coef_split_default_pml = to_num<int>(parameters(nb++));
      }
    
  }
  

  //! initializes a distributed matrix
  template<class Matrix1>
  void DistributedProblem_Base::InitDistributedMatrix(Matrix1& mat)
  {
    // mat is not a distributed matrix, nothing to do
  }
  
  
  //! initializes a distributed matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedProblem_Base
  ::InitDistributedMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& mat)
  {
#ifdef SELDON_WITH_MPI
    if (this->GetNbProcPerMode() > 1)
      {
        int nb_u = 1; int Nvol = var_problem.GetNbDof();
        if (var_problem.GetNbMeshNumberings() == 1)
          {
            nb_u = var_problem.nb_unknowns_scal;
            Nvol = var_problem.GetNbDofUnknown(0);
            if (this->var_problem.FormulationDG() == ElementReference_Base::HDG)
              nb_u = var_problem.nb_unknowns_hdg;
          }
        
        mat.Init(this->nodl_all, &this->GlobDofNumber_Subdomain,
                 &this->OverlapDofNumber_Subdomain,
                 &this->OverlapProcNumber_Subdomain, Nvol, nb_u,
                 &this->MatchingNumber_Subdomain,
                 &this->MatchingDofOrig_Subdomain, this->comm_group_mode);
      }
    else
      mat.Init(mat.GetM(), NULL, NULL, NULL,
               0, 1, NULL, NULL, MPI_COMM_SELF);
#endif
  }


  //! adding contributions of other domains in vector X
  /*!
    \param[in,out] X vector to modify
    \param[in] nb_u number of unknowns in the vector X
  */
  template<class T>
  void DistributedProblem_Base::AddDomains(Vector<T>& X, int nb_u) const
  {
#ifdef SELDON_WITH_MPI
    // for DG formulation, no need to assemble solution
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    int Nvol = var_problem.GetNbDof();
    if (nb_u == -1)
      {
        if (var_problem.GetNbMeshNumberings() == 1)
          {
            Nvol = var_problem.GetNbDofUnknown(0);
            if (var_problem.FormulationDG() == ElementReference_Base::HDG)
              nb_u = var_problem.nb_unknowns_hdg;
            else
              nb_u = var_problem.nb_unknowns;
          }
        else
          nb_u = 1;
      }
    
    ReduceDistributedVector(X, MPI_SUM, this->MatchingNumber_Subdomain,
			    this->MatchingDofOrig_Subdomain, Nvol, nb_u);
#endif
  }


#ifdef SELDON_WITH_MPI
  //! reduction operation for Dirichlet dofs only
  template<class T>
  void DistributedProblem_Base::ReduceDirichlet(Vector<T>& X, const MPI_Op& oper,
                                                int nb_u, bool only_num) const
  {
    // for DG formulation, no need to assemble solution
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;
    
    int Nvol = var_problem.GetNbDof();
    if (only_num)
      {
        if (var_problem.GetNbMeshNumberings() == 1)
          ReduceDistributedVector(X, oper, this->DirichletNumber_Subdomain,
                                  this->DirichletDofOrig_Subdomain, Nvol, 1);
        else
          {
            /* int n0 = -1;
            for (int i = 0; i < var_problem.nb_unknowns; i++)
              if (var_problem.mesh_num_unknown(i) == nb_u)
                n0 = i;
            */
            
            int n0 = nb_u;
            if (n0 < 0)
              return;
            
            const Vector<int>& dir_number = this->DirichletNumber_Subdomain;
            Vector<IVect> dir_dof(dir_number.GetM());
            for (int p = 0; p < dir_number.GetM(); p++)
              {
                int nb_dof = 0;
                for (int j = 0; j < this->DirichletDofOrig_Subdomain(p).GetM(); j++)
                  {
                    int num_dof = this->DirichletDofOrig_Subdomain(p)(j);
                    if ((num_dof >= var_problem.offset_dof_unknown(n0))
                        && (num_dof < var_problem.offset_dof_unknown(n0+1)))
                      nb_dof++;
                  }

                dir_dof(p).Reallocate(nb_dof);
                nb_dof = 0;
                for (int j = 0; j < this->DirichletDofOrig_Subdomain(p).GetM(); j++)
                  {
                    int num_dof = this->DirichletDofOrig_Subdomain(p)(j);
                    if ((num_dof >= var_problem.offset_dof_unknown(n0))
                        && (num_dof < var_problem.offset_dof_unknown(n0+1)))
                      dir_dof(p)(nb_dof++) = num_dof - var_problem.offset_dof_unknown(n0);
                  }
              }
            
            ReduceDistributedVector(X, oper, dir_number, dir_dof, Nvol, 1);
          }
      }
    else
      {
        if (nb_u == -1)
          {
            if (var_problem.GetNbMeshNumberings() == 1)
              {
                Nvol = var_problem.GetNbDofUnknown(0);
                if (var_problem.FormulationDG() == ElementReference_Base::HDG)
                  nb_u = var_problem.nb_unknowns_hdg;
                else
                  nb_u = var_problem.nb_unknowns;
              }
            else
              nb_u = 1;
          }
        
        ReduceDistributedVector(X, oper, this->DirichletNumber_Subdomain,
                                this->DirichletDofOrig_Subdomain, Nvol, nb_u);
      }
  }
#endif

  
  //! adding contributions of other domains (only Dirichlet dofs) in vector X
  /*!
    \param[in,out] X vector to modify
    \param[in] nb_u number of unknowns in the vector X
  */
    template<class T>
  void DistributedProblem_Base::AssembleDirichlet(Vector<T>& X, int nb_u, bool only_num) const
  {
#ifdef SELDON_WITH_MPI
    ReduceDirichlet(X, MPI_SUM, nb_u, only_num);
#endif
  }
  

  //! exchanges values between neighboring domains
  /*!
    \param[in,out] X vector to modify
    \param[in] nb_u number of unknowns in the vector X
  */
  template<class T>
  void DistributedProblem_Base
  ::ExchangeDomains(Vector<T>& X, int nb_u) const
  {
#ifdef SELDON_WITH_MPI
    // for DG formulation, no exchanges are possible
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    int Nvol = var_problem.GetNbDof();
    if (nb_u == -1)
      {
        if (var_problem.GetNbMeshNumberings() == 1)
          {
            Nvol = var_problem.GetNbDofUnknown(0);
            if (var_problem.FormulationDG() == ElementReference_Base::HDG)
              nb_u = var_problem.nb_unknowns_hdg;
            else
              nb_u = var_problem.nb_unknowns;
          }
        else
          nb_u = 1;
      }
    
    ExchangeDistributedVector(X, this->MatchingNumber_Subdomain,
			      this->MatchingDofOrig_Subdomain, Nvol, nb_u);
#endif
  }

  
  //! processor proc sends the values to neighboring processors
  //! and neighboring processors perform a relaxation on these values
  /*!
    \param[in,out] X vector to modify
    \param[in] omega relaxation parameter
    \param[in] proc target processor
    \param[in] nb_u number of unknowns in the vector X
  */
  template<class T>
  void DistributedProblem_Base
  ::ExchangeRelaxDomains(Vector<T>& X,
			 const Real_wp& omega, int proc, int nb_u) const
  {
#ifdef SELDON_WITH_MPI
    // for DG formulation, no exchanges are possible
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    int Nvol = var_problem.GetNbDof();
    if (nb_u == -1)
      {
        if (var_problem.GetNbMeshNumberings() == 1)
          {
            Nvol = var_problem.GetNbDofUnknown(0);
            if (var_problem.FormulationDG() == ElementReference_Base::HDG)
              nb_u = var_problem.nb_unknowns_hdg;
            else
              nb_u = var_problem.nb_unknowns;
          }
        else
          nb_u = 1;
      }
    
    ExchangeRelaxDistributedVector(X, omega, proc, this->MatchingNumber_Subdomain,
				   this->MatchingDofOrig_Subdomain, Nvol, nb_u);
#endif
  }
  

  //! processor proc sends the values to neighboring processors
  //! and neighboring processors perform a relaxation on these values
  /*!
    \param[in,out] X vector to modify
    \param[in] omega relaxation parameter
    \param[in] proc target processor
    \param[in] nb_u number of unknowns in the vector X
  */
  template<class T>
  void DistributedProblem_Base
  ::ExchangeQuadRelaxDomains(Vector<T>& X,
			     const Real_wp& omega, int proc, int nb_u) const
  {
#ifdef SELDON_WITH_MPI
    // for DG formulation, no exchanges are possible
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    int Nvol = var_problem.GetNbDof();
    if (nb_u == -1)
      {
        if (var_problem.GetNbMeshNumberings() == 1)
          {
            Nvol = var_problem.GetNbDofUnknown(0);
            if (var_problem.FormulationDG() == ElementReference_Base::HDG)
              nb_u = var_problem.nb_unknowns_hdg;
            else
              nb_u = var_problem.nb_unknowns;
          }
        else
          nb_u = 1;
      }
    
    ExchangeRelaxDistributedVector(X, omega, proc, this->MatchingNumber_Subdomain,
				   this->MatchingQuadOrig_Subdomain, Nvol, nb_u);
#endif
  }


#ifdef SELDON_WITH_MPI
  void DistributedProblem_Base::ConstructDirichletComm()
  {
    if (MatchingDofOrig_Subdomain.GetM() <= 0)
      return;

    ExtractDistributedSharedNumbers(MatchingNumber_Subdomain,
				    MatchingDofOrig_Subdomain,
				    var_boundary.GetIsDofDirichlet(),
				    var_boundary.GetNbDirichletDof(),
				    DirichletNumber_Subdomain,
				    DirichletDofOrig_Subdomain);
  }
  
  
  //! reducing a distributed vector
  /*!
    \param[in,out] X vector to modify
    \param[in] oper operation to perform (MPI_SUM, MPI_MIN or MPI_MAX)
    \param[in] nb_u number of "vectors"
   */
  template<class T>
  void DistributedProblem_Base
  ::ReduceDistributedVector(Vector<T>& X, const MPI_Op& oper, int nb_u, bool only_num) const
  {

    // for DG formulation, no need to assemble vector
    if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      return;

    int Nvol = var_problem.GetNbDof();
    if (only_num)
      {
        if (var_problem.GetNbMeshNumberings() == 1)
          ReduceDistributedVector(X, oper, this->MatchingNumber_Subdomain,
                                  this->MatchingDofOrig_Subdomain, Nvol, 1);
        else
          {
            int n0 = -1;
            for (int i = 0; i < var_problem.nb_unknowns; i++)
              if (var_problem.mesh_num_unknown(i) == nb_u)
                n0 = i;

            if (n0 < 0)
              return;
            
            const Vector<int>& dir_number = this->MatchingNumber_Subdomain;
            Vector<IVect> dir_dof(dir_number.GetM());
            for (int p = 0; p < dir_number.GetM(); p++)
              {
                int nb_dof = 0;
                for (int j = 0; j < this->MatchingDofOrig_Subdomain(p).GetM(); j++)
                  {
                    int num_dof = this->MatchingDofOrig_Subdomain(p)(j);
                    if ((num_dof >= var_problem.offset_dof_unknown(n0))
                        && (num_dof < var_problem.offset_dof_unknown(n0+1)))
                      nb_dof++;
                  }

                dir_dof(p).Reallocate(nb_dof);
                nb_dof = 0;
                for (int j = 0; j < this->MatchingDofOrig_Subdomain(p).GetM(); j++)
                  {
                    int num_dof = this->MatchingDofOrig_Subdomain(p)(j);
                    if ((num_dof >= var_problem.offset_dof_unknown(n0))
                        && (num_dof < var_problem.offset_dof_unknown(n0+1)))
                      dir_dof(p)(nb_dof++) = num_dof - var_problem.offset_dof_unknown(n0);
                  }
              }
            
            ReduceDistributedVector(X, oper, dir_number, dir_dof, Nvol, 1);
          }
      }
    else
      {
        if (nb_u == -1)
          {
            if (var_problem.GetNbMeshNumberings() == 1)
              {
                Nvol = var_problem.GetNbDofUnknown(0);
                if (var_problem.FormulationDG() == ElementReference_Base::HDG)
                  nb_u = var_problem.nb_unknowns_hdg;
                else
                  nb_u = var_problem.nb_unknowns;
              }
            else
              nb_u = 1;
          }
        
        ReduceDistributedVector(X, oper, this->MatchingNumber_Subdomain,
                                this->MatchingDofOrig_Subdomain, Nvol, nb_u);
      }
  }
  

  //! reducing a vector
  template<class T>
  void DistributedProblem_Base::
  ReduceDistributedVector(Vector<T>& X, const MPI_Op& oper, const IVect& ProcNumber,
			  const Vector<IVect>& DofNumber, int Nvol, int nb_u) const
  {  
    // if only one processor, no communication needed, we leave method
    if (this->GetNbProcPerMode() <= 1)
      return;
    
    glob_chrono.Start(VirtualTimer::COMM);

    AssembleVector(X, oper, ProcNumber, DofNumber, comm_group_mode, Nvol, nb_u, 19);

    glob_chrono.Stop(VirtualTimer::COMM);
  }


  //! exchanges values between domains
  template<class T>
  void DistributedProblem_Base::
  ExchangeDistributedVector(Vector<T>& X, const IVect& ProcNumber,
			    const Vector<IVect>& DofNumber, int Nvol, int nb_u) const
  {  
    // if only one processor, no communication needed, we leave method
    if (this->GetNbProcPerMode() <= 1)
      return;
    
    glob_chrono.Start(VirtualTimer::COMM);

    ExchangeVector(X, ProcNumber, DofNumber, comm_group_mode, Nvol, nb_u, 26);
    
    glob_chrono.Stop(VirtualTimer::COMM);
  }


  //! exchanging a vector of doubles (or complex doubles)
  //! processor proc sends the values to neighboring processors
  //! and neighboring processors perform a relaxation on these values
  template<class T>
  void DistributedProblem_Base::
  ExchangeRelaxDistributedVector(Vector<T>& X, const Real_wp& omega,
				 int proc, const IVect& ProcNumber,
				 const Vector<IVect>& DofNumber, int Nvol, int nb_u) const
  {  
    // if only one processor, no communication needed, we leave method
    if (this->GetNbProcPerMode() <= 1)
      return;
    
    glob_chrono.Start(VirtualTimer::COMM);

    ExchangeRelaxVector(X, omega, proc, ProcNumber, DofNumber, comm_group_mode, Nvol, nb_u, 26);

    glob_chrono.Stop(VirtualTimer::COMM);
  }

  
  //! exchanging values of X between domains
  /*!
    \param[in] X values to send
    \param[out] xsend, xsend_tmp, xrecv, xrecv_tmp arrays used for communications
    \param[in,out] request list of request for all the domains
    \param[in] tag identifying number for messages
  */
  template<class T>
  void DistributedProblem_Base::
  ExchangeUfaceDomains(const Vector<T>& X,
                       Vector<Vector<T> >& xsend, Vector<Vector<int64_t> >& xsend_tmp,
                       Vector<Vector<T> >& xrecv, Vector<Vector<int64_t> >& xrecv_tmp,
                       Vector<MPI_Request>& request, int tag) const
  {
    glob_chrono.Start(VirtualTimer::COMM);
    // if only one processor, no communication needed, we leave method
    if (this->GetNbProcPerMode() <= 1)
      return;
    
    int nb_dom = this->MatchingNumber_Subdomain.GetM();
    xsend.Reallocate(nb_dom);
    xsend_tmp.Reallocate(nb_dom);
    request.Reallocate(2*nb_dom);
    // sending informations to other domains
    for (int i = 0; i < nb_dom; i++)
      {
	int j = this->MatchingNumber_Subdomain(i);
	int nb = this->MatchingQuadOrig_Subdomain(i).GetM();
	if (nb > 0)
	  {
	    xsend(i).Reallocate(nb);
	    for (int k = 0; k < nb; k++)
	      xsend(i)(k) = X(this->MatchingQuadOrig_Subdomain(i)(k));
	    
	    // sending the value to the corresponding processor
            request(i+nb_dom) = MpiIsend(comm_group_mode, xsend(i), xsend_tmp(i), nb, j, tag);
          }
      }
    
    // receiving the informations
    xrecv.Reallocate(nb_dom);
    xrecv_tmp.Reallocate(nb_dom);
    for (int i = 0; i < nb_dom; i++)
      {
	int j = this->MatchingNumber_Subdomain(i);
	int nb = this->MatchingQuadOrig_Subdomain(i).GetM();
	if (nb > 0)
	  {
	    xrecv(i).Reallocate(nb); xrecv(i).Fill(0);
	    
	    // receiving the values of domain j
            request(i) = MpiIrecv(comm_group_mode, xrecv(i), xrecv_tmp(i), nb, j, tag);
	  }
      }
    
    glob_chrono.Stop(VirtualTimer::COMM);
  }
  
  
  //! retrieving values of other domains in vector X
  /*!
    \param[out] X created vector
    \param[out] xdom values received on each domain
    \param[in,out] request list of request for all the domains
    \param[in] tag identifying number for messages
   */
  template<class T>
  void DistributedProblem_Base::
  GetUfaceDomains(Vector<T>& X, Vector<Vector<T> >& xsend, Vector<Vector<int64_t> >& xsend_tmp,
                  Vector<Vector<T> >& xdom, Vector<Vector<int64_t> >& xdom_tmp,
                  Vector<MPI_Request>& request, int tag) const
  {
    if (this->GetNbProcPerMode() <= 1)
      return;
    
    glob_chrono.Start(VirtualTimer::COMM);
    // now waiting all communications are effective
    MPI_Status status;
    int nb_dom = this->MatchingQuadOrig_Subdomain.GetM();
    for (int i = 0; i < nb_dom; i++)
      {
        int nb = this->MatchingQuadOrig_Subdomain(i).GetM();
	if (nb > 0)
          {
            MPI_Wait(&request(i), &status);
            MPI_Wait(&request(nb_dom+i), &status);
          }
      }
    
    for (int i = 0; i < nb_dom; i++)
      {
        if (this->MatchingQuadOrig_Subdomain(i).GetM() > 0)
          {
            MpiCompleteIrecv(xdom(i), xdom_tmp(i), xdom(i).GetM());
            // gathering all the values
            for (int k = 0; k < this->MatchingQuadOrig_Subdomain(i).GetM(); k++)
              X(this->MatchingQuadOrig_Subdomain(i)(k)) = xdom(i)(k);
          }
      }
    
    glob_chrono.Stop(VirtualTimer::COMM);
  }
#endif


  //! copies input parameters of another similar problem
  void DistributedProblem_Base::CopyInputData(const DistributedProblem_Base& var)
  {
    splitting_algorithm = var.splitting_algorithm;
    coef_split_default_domain = var.coef_split_default_domain;
    coef_split_default_pml = var.coef_split_default_pml;
    list_processor_on_group = var.list_processor_on_group;
    split_mesh_epart_user = var.split_mesh_epart_user;
    file_name_partition_subdomain = var.file_name_partition_subdomain;
    radius_splitting_subdomain = var.radius_splitting_subdomain;
    ref_layered_subdivision = var.ref_layered_subdivision;
#ifdef SELDON_WITH_MPI
    comm_group_mode = var.comm_group_mode;
#endif
  }
  
  
  /**********************
   * DistributedProblem *
   **********************/

  
  //! modifies the object with a line of the data file
  template<class Dimension>
  void DistributedProblem<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    DistributedProblem_Base::SetInputData(description_field, parameters);
    VarGeometryProblem<Dimension>::SetInputData(description_field, parameters);
  }
  

  //! return dof numbers of dofs associated with the element i
  template<class Dimension>
  IVect DistributedProblem<Dimension>::GetDofNumberOnElement(int i, int num) const
  {
    if (this->FormulationDG() == ElementReference_Base::HDG)
      {
	int nb_dof_elt = this->GetNbLocalDof(i);	
	IVect num(nb_dof_elt);
	for (int k = 0; k < nb_dof_elt; k++)
	  num(k) = this->OffsetDofV(i) + k;
	
	return num;
      }
    
    if (num == 0)
      return this->mesh_num.number_map.GetDofNumberOnElement(this->mesh_num, i);
    else
      return this->other_mesh_num(num-1)->number_map.
        GetDofNumberOnElement(*this->other_mesh_num(num-1), i);
  }


  //! return dof numbers of dofs associated with the element i
  template<class Dimension>
  IVect DistributedProblem<Dimension>::GetScalarDofNumberOnElement(int i) const
  {
    return GetDofNumberOnElement(i);
  }  


  //! Sets internal dofs (for static condensation)
  template<class Dimension>
  void DistributedProblem<Dimension>
  ::GetInternalNodesElement(int i, int nb_dof_loc,
			    int& nb_dof_edges, int& nb_dof_int,
			    Vector<int>& intern_node) const
  {
    // For HDG formulations
    int nb_u = this->nb_unknowns_hdg;
    nb_dof_edges = this->mesh_num.Element(i).GetNbDof()*nb_u;
    nb_dof_int = nb_dof_loc - nb_dof_edges;
    intern_node.Reallocate(nb_dof_loc);
    int nb = 0;
    for (int j = 0; j < nb_dof_edges; j++)
      intern_node(j) = nb++;
    
    nb = 1;
    for (int j = nb_dof_edges; j < nb_dof_loc; j++)
      intern_node(j) = -nb++;
  }
  

  //! Retrieves dofs after static condensation
  /*!
    \param[out] IndexCondensedRows(i) is equal to -1 if the dof is eliminated, and to index if the dof is kept 
    => The condensed matrix will be a smaller matrix using "index" as row/column numbers
    \param[out] global_row global row numbers of the condensed matrix
    \param[out] overlap_row dofs that are owned by another processor
    \param[out] overlap_proc processor that owns the overlapped dofs
    \param[out] nb_scalar_dof number of dofs per unknown for the condensed matrix
    \param[out] nb_global_dof global size of the condensed matrix
    \param[out] sharing_procs numbers of processors that share dofs with current processor
    \param[out] sharing_rows row numbers of shared dofs
   */
  template<class Dimension>
  inline void DistributedProblem<Dimension>
  ::GetStaticCondensedRows(IVect& IndexCondensedRows, IVect& global_row, IVect& overlap_row, IVect& overlap_proc,
			   int& nb_scalar_dof, int& nb_global_dof,
			   IVect& sharing_procs, Vector<IVect>& sharing_rows) const
  {
    int nb_proc = this->GetNbProcPerMode();
    if (this->FormulationDG() == ElementReference_Base::HDG)
      {
        // HDG formulation : only surface unknowns lambda are kept
	IndexCondensedRows.Reallocate(this->GetNbDof());
	IndexCondensedRows.Fill(-1);
	int offset = this->offset_dof_unknown(this->nb_unknowns_hdg);
	if (nb_proc > 1)
	  global_row.Reallocate(offset);
        
	for (int i = 0; i < offset; i++)
	  {
	    IndexCondensedRows(i) = i;
	    if (nb_proc > 1)
	      global_row(i) = this->GlobDofNumber_Subdomain(i);
	  }
        
	overlap_row = this->OverlapDofNumber_Subdomain;
	overlap_proc = this->OverlapProcNumber_Subdomain;
        if (nb_proc > 1)
          nb_global_dof = this->offset_dof_unknown_all(this->nb_unknowns_hdg);
        else
          nb_global_dof = this->offset_dof_unknown(this->nb_unknowns_hdg);
        
	if (this->GetNbMeshNumberings() == 1)
          nb_scalar_dof = this->offset_dof_unknown(1);
        else
          nb_scalar_dof = offset;
        
	sharing_procs = this->MatchingNumber_Subdomain;
	sharing_rows = this->MatchingDofOrig_Subdomain;
	
	return;
      }

    if (this->GetComputationProblem().LightStaticCondensation())
      {
        // case where only discontinuous unknowns (vectorial) are eliminated
        if (this->GetNbMeshNumberings() == 1)
          nb_scalar_dof = this->offset_dof_unknown(1);
        else
          nb_scalar_dof = this->offset_dof_unknown(this->nb_unknowns_scal);
        
	nb_global_dof = nb_scalar_dof;

        int Ns = this->offset_dof_unknown(this->nb_unknowns_scal);
	IndexCondensedRows.Reallocate(this->GetNbDof());
	IndexCondensedRows.Fill(-1);
	for (int i = 0; i < Ns; i++)
	  IndexCondensedRows(i) = i;
        
	if (nb_proc > 1)
	  {
	    global_row.Reallocate(Ns);
	    for (int i = 0; i < Ns; i++)
	      global_row(i) = GlobDofNumber_Subdomain(i);

	    overlap_row = this->OverlapDofNumber_Subdomain;
	    overlap_proc = this->OverlapProcNumber_Subdomain;

	    sharing_procs = this->MatchingNumber_Subdomain;
	    sharing_rows = this->MatchingDofOrig_Subdomain;

	    nb_global_dof = this->offset_dof_unknown_all(this->nb_unknowns_scal);
	  }
        
	return;
      }
    
    IndexCondensedRows.Reallocate(this->GetNbDof());
    IndexCondensedRows.Fill(-1);
    int nb = 0; // int Ncond = offset_dof_condensed_all(this->nb_unknowns_scal);
    if (nb_proc > 1)
      global_row.Reallocate(this->offset_dof_condensed(this->nb_unknowns_scal));

    for (int j = 0; j < this->nb_unknowns_scal; j++)
      for (int i0 = this->offset_dof_condensed(j); i0 < this->offset_dof_condensed(j+1); i0++)
	{
          int i = i0 - this->offset_dof_condensed(j);
	  if (nb_proc > 1)
            {
              int p = this->GlobDofNumber_Subdomain(i + this->offset_dof_unknown(j));
              p -= this->offset_dof_unknown_all(j);
              global_row(nb) = p + this->offset_dof_condensed_all(j);
            }
          
	  IndexCondensedRows(this->offset_dof_unknown(j) + i) = nb++;
	}
    
    nb_scalar_dof = this->offset_dof_condensed(1);
    if (this->GetNbMeshNumberings() > 1)
      nb_scalar_dof = this->offset_dof_condensed(this->nb_unknowns_scal);
    
    nb_global_dof = nb_scalar_dof;    
    if (nb_proc > 1)
      {
	int noverlap_old = this->OverlapDofNumber_Subdomain.GetM();
	int noverlap = 0;
	for (int i = 0; i < noverlap_old; i++)
	  if (IndexCondensedRows(this->OverlapDofNumber_Subdomain(i)) >= 0)
            noverlap++;
	
	overlap_row.Reallocate(noverlap);
	overlap_proc.Reallocate(noverlap);
	noverlap = 0;
        for (int i = 0; i < noverlap_old; i++)
	  if (IndexCondensedRows(this->OverlapDofNumber_Subdomain(i)) >= 0)
            {
              overlap_row(noverlap) = IndexCondensedRows(this->OverlapDofNumber_Subdomain(i));
              overlap_proc(noverlap) = this->OverlapProcNumber_Subdomain(i);
              noverlap++;
            }
        
	sharing_procs = this->MatchingNumber_Subdomain;
	sharing_rows.Reallocate(sharing_procs.GetM());
	for (int p = 0; p < sharing_procs.GetM(); p++)
	  {
	    int nb_dof = 0;
	    for (int i = 0; i < this->MatchingDofOrig_Subdomain(p).GetM(); i++)
	      if (IndexCondensedRows(this->MatchingDofOrig_Subdomain(p)(i)) >= 0)
                nb_dof++;

            sharing_rows(p).Reallocate(nb_dof);
            nb_dof = 0;
            for (int i = 0; i < this->MatchingDofOrig_Subdomain(p).GetM(); i++)
	      if (IndexCondensedRows(this->MatchingDofOrig_Subdomain(p)(i)) >= 0)
                {
                  sharing_rows(p)(nb_dof) = IndexCondensedRows(this->MatchingDofOrig_Subdomain(p)(i));
                  nb_dof++;
                }
	  }

        nb_global_dof = this->offset_dof_condensed_all(this->nb_unknowns_scal);
      }
  }

  
  //! adds memory used by the object
  template<class Dimension>
  void DistributedProblem<Dimension>::GetMemoryUsed(map<string, size_t>& var) const
  {
    DistributedProblem_Base::GetMemoryUsed(var);
    VarGeometryProblem<Dimension>::GetMemoryUsed(var);
  }

  
  //! splitting the mesh into several parts and distribute dofs between processors
  template<class Dimension>
  void DistributedProblem<Dimension>
  ::SplitMeshForParallelComputation(const string& name_elt)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc = this->GetNbProcPerMode();
    int rank_proc = this->GetRankProcMode();
    
    // if multi-processor, the mesh in split in small sub-meshes
    // for DG formulation, the algorithm is slightly different
    // because the numbers are not shared, it is discontinuous
    // so you need to know the neighbors
    if (nb_proc > 1)
      { 
        // splitting mesh into several parts
        this->SplitSubdomains(name_elt);    
		
	// degrees of freedom due to mesh numberings
        this->GlobDofNumber_Subdomain.Resize(this->nodl);
        this->GlobDofNumber_Subdomain.Fill(-1);

        // loop over unknowns
        int nb_u = this->nb_unknowns;
        if (this->FormulationDG() == ElementReference_Base::HDG)
          nb_u = this->nb_unknowns_hdg;
        
        for (int i = 0; i < nb_u; i++)
          {
            int n = this->mesh_num_unknown(i);
            MeshNumbering<Dimension>* mesh_num;
            if (n == 0)
              mesh_num = &this->mesh_num;
            else
              mesh_num = this->other_mesh_num(n-1);

            int offset = this->offset_dof_unknown(i);
            int offset_glob = this->offset_dof_unknown_all(i);
            int Nmax = mesh_num->GlobDofNumber_Subdomain.GetM();
            Nmax = min(Nmax, this->offset_dof_unknown(i+1) - offset);
            for (int j = 0; j < Nmax; j++)
              this->GlobDofNumber_Subdomain(offset + j) = offset_glob + mesh_num->GlobDofNumber_Subdomain(j);
          }
        
        if (this->FormulationDG() == ElementReference_Base::HDG)
          {
            int offset_uloc = this->offset_dof_unknown(this->nb_unknowns_hdg);
            int offset_ug = this->offset_dof_unknown_all(this->nb_unknowns_hdg);
            // vectorial unknowns for HDG formulation
            int Nvol = this->OffsetDofV(this->mesh.GetNbElt());
            int Nvol_all = (this->offset_dof_unknown_all(this->nb_unknowns_hdg + this->nb_unknowns) - offset_ug) / this->nb_unknowns;
            for (int i = 0; i < this->mesh.GetNbElt(); i++)
              for (int j = this->OffsetDofV(i); j < this->OffsetDofV(i+1); j++)
                for (int m = 0; m < this->nb_unknowns; m++)
                  this->GlobDofNumber_Subdomain(offset_uloc + j + m*Nvol)
                    = offset_ug + this->GlobOffsetDofV_Subdomain(i) + j - this->OffsetDofV(i) + m*Nvol_all;			
          }
        else if (this->FormulationDG() == ElementReference_Base::CONTINUOUS)
          {
            int offset_uloc = this->offset_dof_unknown(this->nb_unknowns);
            int offset_ug = this->offset_dof_unknown_all(this->nb_unknowns);
            // vectorial dofs for continuous formulation
            if (this->OffsetDofV.GetM() > 0)
              for (int i = 0; i < this->mesh.GetNbElt(); i++)
                for (int j = this->OffsetDofV(i); j < this->OffsetDofV(i+1); j++)
                  this->GlobDofNumber_Subdomain(offset_uloc + j)
                    = offset_ug + this->GlobOffsetDofV_Subdomain(i) + j - this->OffsetDofV(i);
          }
        
	// looking for PML
	this->FindElementsInsidePML();
        
        // filling missing dofs (assuming that there are no shared dofs)
        int ddl_max = 0; int nb_dof_m1 = 0;
        for (int i = 0; i < this->GlobDofNumber_Subdomain.GetM(); i++)
          {
            int n = this->GlobDofNumber_Subdomain(i);
            ddl_max = max(n, ddl_max);
            if (n == -1)
              nb_dof_m1++;
          }
        
        Vector<int> num_rcv, num_send;
        if (rank_proc == 0)
          num_rcv.Reallocate(2*nb_proc);

        TinyVector<int, 2> num_d;
        num_d(0) = ddl_max; num_d(1) = nb_dof_m1;
        MPI_Gather(&num_d(0), 2, MPI_INT, num_rcv.GetData(), 2, MPI_INT, 0,  this->comm_group_mode);
        
        int offset;
        if (rank_proc == 0)
          {
            int nb_all_m1 = 0;
            for (int i = 0; i < nb_proc; i++)
              {
                nb_all_m1 += num_rcv(2*i+1);
                ddl_max = max(ddl_max, num_rcv(2*i));
              }
            
            ddl_max++;
            if (nb_all_m1 != this->nodl_all - ddl_max)
              {
                cout << "Inconsistent numbering in parallel" << endl;
                DISP(this->nodl_all); DISP(ddl_max); DISP(nb_all_m1);
                abort();
              }

            num_send.Reallocate(nb_proc);
            offset = ddl_max;
            for (int i = 0; i < nb_proc; i++)
              {
                num_send(i) = offset;
                offset += num_rcv(2*i+1);
              }
          }

        MPI_Scatter(num_send.GetData(), 1, MPI_INT, &offset, 1, MPI_INT, 0, this->comm_group_mode);
        for (int i = 0; i < this->GlobDofNumber_Subdomain.GetM(); i++)
          {
            int n = this->GlobDofNumber_Subdomain(i);
            if (n == -1)
              this->GlobDofNumber_Subdomain(i) = offset++;
          }
        
        if (this->print_level >= 6)
          this->GlobDofNumber_Subdomain.WriteText("numGlob" + to_str(rank_proc) + ".dat");
      }
#endif
  }
  
    
  //! splitting of the mesh for parallel computation
  /*!
    \param[in] map informations for the numbering of the mesh
   */
  template<class Dimension>
  void DistributedProblem<Dimension>::SplitSubdomains(const string& name_elt)
  {
#ifdef SELDON_WITH_MPI
    
    Vector<IVect> NumElement_Subdomain;
    Vector<Mesh<Dimension> > mesh_subdomain;
    
    int rank; MPI_Comm_rank(this->comm_group_mode, &rank);
    int nb_proc; MPI_Comm_size(this->comm_group_mode, &nb_proc);
    
#ifdef MONTJOIE_WITH_TRANSMISSION
    int nb_surfaces_transmission, nb_ddl_transmission;
    IVect InfoSurfTransmission, InfoSurfMatchingTransmission;
#endif
    
    if (this->mesh.GetNbPeriodicReferences() > 0)
      if (this->mesh_num.GetFormulationForPeriodicCondition() == this->mesh_num.STRONG_PERIODIC)
	{
	  cout << "Formulation not handled in parallel" << endl;
	  cout << "Use the weak formulation, UseSameDofsForPeriodicCondition = WEAK" << endl;
	  cout << "If strictly periodic condition is sufficient, UseSameDofsForPeriodicCondition = YES" << endl;
	  abort();
	}

    Vector<IVect> offset_matching;
    if (rank == 0)
      {		
	//! splitting of the inital mesh into subdomains with Scotch or metis
	int nb_subdiv = nb_proc;
	
        // by default, each element has a weight equal to one
	IVect weight_elt(this->mesh.GetNbElt());
        weight_elt.Fill(this->coef_split_default_domain);
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          if (this->mesh.Element(i).IsPML())
            weight_elt(i) = this->coef_split_default_pml;
        
	switch(this->splitting_algorithm)
	  {
	  case DistributedProblem_Base::SUBDIV_DOMAIN_METIS :
	    this->mesh.SplitMetis(nb_subdiv, weight_elt, NumElement_Subdomain, mesh_subdomain,
				  0, false);
	    break;
	  case DistributedProblem_Base::SUBDIV_DOMAIN_SCOTCH :
	    this->mesh.SplitScotch(nb_subdiv, weight_elt, NumElement_Subdomain, mesh_subdomain,
				   0, false);
	    break;
	  case DistributedProblem_Base::SUBDIV_DOMAIN_CONCENTRIC :
	    this->mesh.SplitConcentric(nb_subdiv, this->radius_splitting_subdomain,
				       NumElement_Subdomain, mesh_subdomain, 0,
				       false);
	    break;
	  case DistributedProblem_Base::SUBDIV_DOMAIN_LAYERED :
	    this->mesh.SplitLayered(nb_subdiv, this->ref_layered_subdivision,
				    NumElement_Subdomain, mesh_subdomain, 0,
				    false);
	    break;
	  case DistributedProblem_Base::SUBDIV_DOMAIN_USER :
            // user-supplied epart
            this->mesh.PartMesh(nb_subdiv, this->split_mesh_epart_user,
				NumElement_Subdomain, mesh_subdomain, 0,
				false);
	    break;
	  default :
	    {
	      cout << "Unknown splitting algorithm" << endl;
	      abort();
	    }
          }
        
	// random generator is reinitialized
	srand(Seed);

        if (this->print_level >= 10)	
          this->mesh.print_level = 6;
        else
          this->mesh.print_level = 0;

        // Filling Epart(i) = processor number for each element
        // and NumLoc(i) = local number of each element in the processor
        weight_elt.Clear();
        IVect Epart(this->mesh.GetNbElt());
        IVect NumLoc(this->mesh.GetNbElt());
        for (int i = 0; i < NumElement_Subdomain.GetM(); i++)
          for (int j = 0; j < NumElement_Subdomain(i).GetM(); j++)
            {
              int ne = NumElement_Subdomain(i)(j);
              Epart(ne) = i;
              NumLoc(ne) = j;
            }
        
        if (this->partition_subdomain_to_store)
          Epart.Write(this->file_name_partition_subdomain);          
	
        if (nb_subdiv != nb_proc)
	  {
	    cout << "The number of subdomains is different "
                 << "from the number of processors \n. Exiting ..." << endl;
	    abort();
	  }
        
        if (this->mesh_num.IsOrderVariable())
	  this->ComputeVariableOrder();

        // numbering the global mesh
        this->ConstructFiniteElement(name_elt);        
        this->mesh_num.NumberMesh();
        for (int i = 0; i < this->other_mesh_num.GetM(); i++)
          this->other_mesh_num(i)->NumberMesh();
        
	// updating periodic references
	this->mesh.UpdatePeriodicReferenceSplit(Epart, NumLoc, mesh_subdomain);
	
#ifdef MONTJOIE_WITH_TRANSMISSION
        // treating global dofs associated with transmission conditions
        this->TreatTransmission(Epart);
        
        // sending informations to other processors
        this->SendTransmissionDofs(Epart, nb_surfaces_transmission, nb_ddl_transmission,
				   InfoSurfTransmission, InfoSurfMatchingTransmission);
#endif
        
        // treating global dofs associated with impedance conditions
        this->TreatGibc(Epart);
	
        // computation of dofs for the global problem
        this->ComputeNumberOfDofs();
        
        this->nodl_mesh.Reallocate(this->GetNbMeshNumberings());
        this->nodl_pml.Reallocate(this->GetNbMeshNumberings());
        this->nodl_mesh(0) = this->mesh_num.GetNbDof();
        this->nodl_pml(0) = this->mesh_num.GetNbDofPML();
        for (int i = 0; i < this->other_mesh_num.GetM(); i++)
          {
            this->nodl_mesh(i+1) = this->other_mesh_num(i)->GetNbDof();
            this->nodl_pml(i+1) = this->other_mesh_num(i)->GetNbDofPML();
          }
        
        this->nodl_all = this->nodl;
        this->offset_dof_unknown_all = this->offset_dof_unknown;
        this->offset_dof_condensed_all = this->offset_dof_condensed;
        
        // brodcasting the following variables :
        // nodl_mesh : global number of degrees on freedom on the mesh 
        // nodl_pml : global number of degrees on freedom associated with PML elements
        // nodl_all : global number of degrees of freedom for the considered problem
        // offset_dof* : dof offsets for each unknown
        int No = this->offset_dof_unknown.GetM();
        MPI_Bcast(this->nodl_mesh.GetData(), this->nodl_mesh.GetM(), MPI_INTEGER, 0, this->comm_group_mode);
        MPI_Bcast(this->nodl_pml.GetData(), this->nodl_mesh.GetM(), MPI_INTEGER, 0, this->comm_group_mode);
        MPI_Bcast(&this->nodl_all, 1, MPI_INTEGER, 0, this->comm_group_mode);
        MPI_Bcast(this->offset_dof_condensed_all.GetData(), No, MPI_INTEGER, 0, this->comm_group_mode);
	MPI_Bcast(this->offset_dof_unknown_all.GetData(), No, MPI_INTEGER, 0, this->comm_group_mode);
        
        // minimal processor(+1) for each dof
        ParamParallelMesh<Dimension> param;
        Vector<IVect>& MinimalProc = param.MinimalProc;
        MinimalProc.Reallocate(this->nodl_mesh.GetM());
        param.dg_form = this->FormulationDG();
        param.nodl_mesh = this->nodl_mesh;
        for (int n = 0; n < this->nodl_mesh.GetM(); n++)
          if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
            {
              MinimalProc(n).Reallocate(this->nodl_mesh(n));
              MinimalProc(n).Fill(nb_proc+1);
              for (int i = 0; i < this->mesh.GetNbElt(); i++)
                {
                  int nb_ddl = 0;
                  if (n == 0)
                    nb_ddl = this->mesh_num.Element(i).GetNbDof();
                  else
                    nb_ddl = this->other_mesh_num(n-1)->Element(i).GetNbDof();
                  
                  for (int j = 0; j < nb_ddl; j++)
                    {
                      int num_dof = 0;
                      if (n == 0)
                        num_dof = this->mesh_num.Element(i).GetNumberDof(j);
                      else
                        num_dof = this->other_mesh_num(n-1)->Element(i).GetNumberDof(j);
                      
                      if (num_dof >= 0)
                        MinimalProc(n)(num_dof) = min(MinimalProc(n)(num_dof), Epart(i)+1);
                    }
                }
            }
        
        Mesh<Dimension> mesh_root, sub_mesh;
        mesh_root.CopyInputData(this->mesh);
        Vector<MeshNumbering<Dimension>* > all_mesh_num(this->nodl_mesh.GetM());
	Vector<MeshNumbering<Dimension>* > mesh_root_num(this->nodl_mesh.GetM());
        for (int i = 0; i < this->nodl_mesh.GetM(); i++)
          {
            if (i == 0)
              all_mesh_num(i) = &this->mesh_num;
            else
              all_mesh_num(i) = this->other_mesh_num(i-1);
            
            mesh_root_num(i) = new MeshNumbering<Dimension>(mesh_root);        
            mesh_root_num(i)->SetOrder(all_mesh_num(i)->GetOrder());
            mesh_root_num(i)->compute_dof_pml = all_mesh_num(i)->compute_dof_pml;
          }
        
        if (this->print_level >= 8)
          this->mesh.Write(string("initial.mesh"));
        
        // loop on all the subdomains
	for (int i = 0; i < nb_proc; i++)
          {
            sub_mesh = mesh_subdomain(i);
            
            // checking that the number of elements is non-null
            if (mesh_subdomain(i).GetNbElt() == 0)
              {
                cout<<"The subdomain " << i << " has no elements " << endl;
		abort();
              }
            
            if (this->OffsetDofV.GetM() > 0)
              {
                int nb_elt = NumElement_Subdomain(i).GetM();
                param.OffsetDofV.Reallocate(nb_elt);          
                
                for (int j = 0; j < nb_elt; j++)
                  {
                    int ne = NumElement_Subdomain(i)(j);
                    param.OffsetDofV(j) = this->OffsetDofV(ne);
                  }
              }
            
            // sending mesh to the appropriate processor
            SendMeshToProcessor(this->mesh, all_mesh_num, NumElement_Subdomain(i), Epart, NumLoc,
                                sub_mesh, i, this->comm_group_mode, param);
            
            if (i == 0)
              {
                mesh_root_num(0)->number_map = this->mesh_num.number_map;
                for (int k = 0; k < this->other_mesh_num.GetM(); k++)
                  mesh_root_num(k+1)->number_map = this->other_mesh_num(k)->number_map;
                
                RecvMeshFromProcessor(mesh_root, mesh_root_num, NeighboringConnectivity,
                                      this->MatchingNumber_Subdomain,  this->MatchingDofOrig_Subdomain,
                                      param, 0, this->comm_group_mode);
                
                this->GlobOffsetDofV_Subdomain = param.OffsetDofV;
                offset_matching = param.offset_matching;
              }
            
            if (this->print_level >= 7)
              sub_mesh.Write(string("subdomain") + to_str(i) + ".mesh");
          }

	// global mesh is cleared
        this->mesh.ClearPeriodicCondition();
	this->mesh.Clear();
        this->mesh_num.Clear();
        
        // small meshes are cleared and mesh of root processor is kept
        mesh_subdomain.Clear();
        NumElement_Subdomain.Clear();
        this->mesh = mesh_root;
	this->mesh_num = *mesh_root_num(0);
        for (int i = 0; i < this->other_mesh_num.GetM(); i++)
          *this->other_mesh_num(i) = *mesh_root_num(i+1);
        
        for (int i = 0; i < this->other_mesh_num.GetM(); i++)
          delete mesh_root_num(i);
      }
    else
      {
	// random generator is reinitialized (to have the same seed for all procs)
	srand(Seed);
	
        if (this->print_level >= 10)	
          this->mesh.print_level = 7;
        else
          this->mesh.print_level = 0;
        
#ifdef MONTJOIE_WITH_TRANSMISSION	
        IVect EpartVoid;
        this->SendTransmissionDofs(EpartVoid, nb_surfaces_transmission, nb_ddl_transmission,
				   InfoSurfTransmission, InfoSurfMatchingTransmission);
#endif

	this->InitGibcReferences(this->mesh.GetNbReferences());
        
        // retrieving global number of degrees of freedom
        this->nodl_mesh.Reallocate(this->GetNbMeshNumberings());
        this->nodl_pml.Reallocate(this->GetNbMeshNumberings());
        int No = this->nb_unknowns + 1;
        if (this->dg_formulation == ElementReference_Base::HDG)
          No = this->nb_unknowns + this->nb_unknowns_hdg + 1;

        this->offset_dof_condensed_all.Reallocate(No);
        this->offset_dof_unknown_all.Reallocate(No);
        
        MPI_Bcast(this->nodl_mesh.GetData(), this->nodl_mesh.GetM(), MPI_INTEGER, 0, this->comm_group_mode);
        MPI_Bcast(this->nodl_pml.GetData(), this->nodl_mesh.GetM(), MPI_INTEGER, 0, this->comm_group_mode);
        MPI_Bcast(&this->nodl_all, 1, MPI_INTEGER, 0, this->comm_group_mode);
        MPI_Bcast(this->offset_dof_condensed_all.GetData(), No, MPI_INTEGER, 0, this->comm_group_mode);
	MPI_Bcast(this->offset_dof_unknown_all.GetData(), No, MPI_INTEGER, 0, this->comm_group_mode);
        
        ParamParallelMesh<Dimension> param;
        param.nodl_mesh = this->nodl_mesh;
        param.dg_form = this->FormulationDG();

        Vector<MeshNumbering<Dimension>* > all_mesh_num(this->nodl_mesh.GetM());
        for (int i = 0; i < this->nodl_mesh.GetM(); i++)
          {
            if (i == 0)
              all_mesh_num(i) = &this->mesh_num;
            else
              all_mesh_num(i) = this->other_mesh_num(i-1);
          }
        
        // receiving the mesh associated with the processor
        this->mesh.ClearPeriodicCondition();
        RecvMeshFromProcessor(this->mesh, all_mesh_num, NeighboringConnectivity,
                              this->MatchingNumber_Subdomain,  this->MatchingDofOrig_Subdomain,
                              param, 0, this->comm_group_mode);
        
        this->GlobOffsetDofV_Subdomain = param.OffsetDofV;
        offset_matching = param.offset_matching;
      }

    // checking the new mesh
    this->CheckInputMesh();
    
    // sorting boundaries in order to have the same order in two different processors
    IVect global_number, local_number;
    global_number.Reallocate(this->mesh.GetNbBoundaryRef());
    local_number.Reallocate(this->mesh.GetNbBoundaryRef());
    local_number.Fill();
    global_number.Fill(-1);
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
        if (this->NeighboringConnectivity(i).GetM() > 0)
          global_number(i) = this->NeighboringConnectivity(i)(0);
      }
    
    Sort(global_number, local_number);

    // rearranging dofs
    int Nm = this->other_mesh_num.GetM() + 1;    
    Vector<Vector<IVect> > Nodle(Nm); Vector<IVect> IndexPML(Nm);
    Vector<int> nb_dof_pml(Nm); nb_dof_pml.Zero();
    for (int n = 0; n < Nm; n++)
      if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        {
          MeshNumbering<Dimension>* mesh_num;
          if (n == 0)
            mesh_num = &this->mesh_num;
          else
            mesh_num = this->other_mesh_num(n-1);
          
          Nodle(n).Reallocate(this->mesh.GetNbElt());
          for (int i = 0; i < this->mesh.GetNbElt(); i++)
            {
              int nb_dof = mesh_num->Element(i).GetNbDof();
              Nodle(n)(i).Reallocate(nb_dof);
              for (int j = 0; j < nb_dof; j++)
                Nodle(n)(i)(j) = mesh_num->Element(i).GetNumberDof(j);
            }        
          
          nb_dof_pml(n) = mesh_num->GetNbDofPML();
          if (nb_dof_pml(n) > 0)
            {
              IndexPML(n).Reallocate(mesh_num->GetNbDof());
              for (int i = 0; i < mesh_num->GetNbDof(); i++)
                IndexPML(n)(i) = mesh_num->GetDofPML(i);
            }
        }
    
    this->ClearFiniteElement();
    this->ConstructFiniteElement(name_elt);
    this->mesh_num.treat_periodic_condition_during_number = false;
    this->mesh_num.NumberMesh();
    for (int n = 0; n < Nm-1; n++)
      {
        this->other_mesh_num(n)->treat_periodic_condition_during_number = false;
        this->other_mesh_num(n)->NumberMesh();
      }
    
    if (this->GetMeshNumbering(0).number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
        int nb_proc_interac = 0;
        Vector<bool> ProcUsed(nb_proc);
        ProcUsed.Fill(false);
        
        for (int j = 0; j < this->mesh.GetNbBoundaryRef(); j++)
          {
            int jloc = local_number(j);
            if (this->NeighboringConnectivity(jloc).GetM() > 0)
              {
                int proc = this->NeighboringConnectivity(jloc)(2);
                if (!ProcUsed(proc))
                  {
                    ProcUsed(proc) = true;
                    nb_proc_interac++;
                  }
              }
          }    
        
        // we list the processors interacting with the current one in 
        // array MatchingNumber_Subdomain
        this->MatchingNumber_Subdomain.Reallocate(nb_proc_interac);
        nb_proc_interac = 0;
        for (int i = 0; i < nb_proc; i++)
          if (ProcUsed(i))
            {
              this->MatchingNumber_Subdomain(nb_proc_interac) = i;
              nb_proc_interac++;
            }
      } 
    
    Vector<Vector<IVect> > matching_dof(Nm);
    for (int n = 0; n < Nm; n++)
      if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        {
          MeshNumbering<Dimension>* mesh_num;
          if (n == 0)
            mesh_num = &this->mesh_num;
          else
            mesh_num = this->other_mesh_num(n-1);
          
          int nb_dof = 0;
          for (int i = 0; i < this->mesh.GetNbElt(); i++)
            {
              mesh_num->Element(i).ReallocateDof(Nodle(n)(i).GetM());
              for (int j = 0; j < Nodle(n)(i).GetM(); j++)
                {
                  nb_dof = max(nb_dof, Nodle(n)(i)(j));
                  mesh_num->Element(i).SetNumberDof(j, Nodle(n)(i)(j));
                }
            }
          
          mesh_num->SetNbDof(nb_dof+1);
          mesh_num->ReconstructOffsetDofs();
          
          if (nb_dof_pml(n) > 0)
            {
              mesh_num->ReallocateDofPML(nb_dof_pml(n));
              for (int i = 0; i < mesh_num->GetNbDof(); i++)
                mesh_num->SetDofPML(i, IndexPML(n)(i));
            }
          
          // global dofs are sorted
          IVect GlobalDof(mesh_num->GlobDofNumber_Subdomain);
          IVect LocalDof(GlobalDof.GetM());
          LocalDof.Fill();
          
          Sort(GlobalDof.GetM(), GlobalDof, LocalDof);
          
          // array MatchingDofOrig_Subdomain is updated with local numbers
          // we don't use a Glob_to_local array (may be slow, but less expensive in memory)
          matching_dof(n).Reallocate(this->MatchingNumber_Subdomain.GetM());
          for (int p = 0; p < this->MatchingNumber_Subdomain.GetM(); p++)
            {
              // index in LocalDof array
              int k = 0; int jmin = offset_matching(p)(n);
              int jmax = offset_matching(p)(n+1);
              matching_dof(n)(p).Reallocate(jmax-jmin);
              for (int j = jmin; j < jmax; j++)
                {
                  while (GlobalDof(k) < this->MatchingDofOrig_Subdomain(p)(j))
                    k++;
                  
                  this->MatchingDofOrig_Subdomain(p)(j) = LocalDof(k);
                  matching_dof(n)(p)(j-jmin) = LocalDof(k);
                }
              //DISP(p); DISP(this->MatchingDofOrig_Subdomain(p));
            }        
        }
    
#ifdef MONTJOIE_WITH_TRANSMISSION	
    this->DistributeTransmissionDofs(nb_surfaces_transmission, nb_ddl_transmission,
                                     InfoSurfTransmission, InfoSurfMatchingTransmission);
#endif
        
    this->ComputeNumberOfDofs();          
    this->PutOtherGlobalDofs();
        
    if (this->print_level >= 7)
      cout <<"Trying to find the minimal processor for each dof" << endl;
    
    int nb_u = this->nb_unknowns;
    if (this->FormulationDG() == ElementReference_Base::HDG)
      nb_u = this->nb_unknowns_hdg;
    
    Vector<IVect> overlap_dof(Nm), overlap_proc(Nm);
    Vector<Vector<IVect> > overlap_number(Nm);
    for (int n = 0; n < Nm; n++)
      if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        {
          MeshNumbering<Dimension>* mesh_num;
          if (n == 0)
            mesh_num = &this->mesh_num;
          else
            mesh_num = this->other_mesh_num(n-1);
          
          int nb_dof_mesh = mesh_num->GetNbDof();
          for (int i = 0; i < nb_u; i++)
            if (this->mesh_num_unknown(i) == n)
              nb_dof_mesh = this->GetNbDofUnknown(i);
          
          // for each dof, we search the minimal processor
          IVect num_proc(nb_dof_mesh);
          num_proc.Fill(rank);
          this->ReduceDistributedVector(num_proc, MPI_MIN, this->MatchingNumber_Subdomain,
                                        matching_dof(n), nb_dof_mesh, 1);
          //DISP(num_proc);
          if (this->print_level >= 7)
            cout <<"minimal processor number retrieved" << endl;
          
          int nb_dof = 0;
          for (int i = 0; i < nb_dof_mesh; i++)
            if (num_proc(i) < rank)
              nb_dof++;
          
          // overlapped dofs
          overlap_dof(n).Reallocate(nb_dof);
          overlap_proc(n).Reallocate(nb_dof);
          
          nb_dof = 0;        
          for (int i = 0; i < nb_dof_mesh; i++)
            if (num_proc(i) < rank)
              {
                overlap_dof(n)(nb_dof) = i;
                overlap_proc(n)(nb_dof) = num_proc(i);
                nb_dof++;
              }
          
          /*overlap_number(n).Reallocate(nb_proc);
            for (int i = 0; i < matching_dof(n).GetM(); i++)
            {
            int proc2 = this->MatchingNumber_Subdomain(i);
            if (proc2 > rank)
            {  
            nb_dof = 0;
            for (int j = 0; j < matching_dof(n)(i).GetM(); j++)
                      if (num_proc(matching_dof(n)(i)(j)) == rank)
                        nb_dof++;
                    
                    overlap_number(n)(proc2).Reallocate(nb_dof);
                    nb_dof = 0;
                    for (int j = 0; j < matching_dof(n)(i).GetM(); j++)
                      if (num_proc(matching_dof(n)(i)(j)) == rank)
                        {
                          overlap_number(n)(proc2)(nb_dof) = matching_dof(n)(i)(j);
                          nb_dof++;
                        }
                  }
                  }*/
        }
    
    // constructing arrays OverlapDofNumber, OverlapProcNumber
    int nb_all_over = 0;
    for (int i = 0; i < nb_u; i++)
      nb_all_over += overlap_dof(this->mesh_num_unknown(i)).GetM();
    
    this->OverlapDofNumber_Subdomain.Reallocate(nb_all_over);
    this->OverlapProcNumber_Subdomain.Resize(nb_all_over);
    nb_all_over = 0;
    for (int i = 0; i < nb_u; i++)
      {
        int n = this->mesh_num_unknown(i);
        if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          for (int j = 0; j < overlap_dof(n).GetM(); j++)
            {
              this->OverlapDofNumber_Subdomain(nb_all_over) = overlap_dof(n)(j) + this->offset_dof_unknown(i);
              this->OverlapProcNumber_Subdomain(nb_all_over) = overlap_proc(n)(j);
              nb_all_over++;
            }
      }
    
    // updating MatchingDofOrig_Subdomain
    if (Nm > 1)
      {
        Vector<int> offset(this->MatchingNumber_Subdomain.GetM());
        offset.Zero();
        for (int i = 0; i < nb_u; i++)
          {
            int n = this->mesh_num_unknown(i);
            if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
              for (int p = 0; p < this->MatchingNumber_Subdomain.GetM(); p++)
                offset(p) += matching_dof(n)(p).GetM();
          }
        
        for (int p = 0; p < this->MatchingNumber_Subdomain.GetM(); p++)
          {
            this->MatchingDofOrig_Subdomain(p).Reallocate(offset(p));
            this->MatchingDofOrig_Subdomain(p).Fill(-1);
          }
        
        offset.Zero();
        for (int i = 0; i < nb_u; i++)
          {
            int n = this->mesh_num_unknown(i);
            if (this->GetMeshNumbering(n).number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
              for (int p = 0; p < this->MatchingNumber_Subdomain.GetM(); p++)
                {
                  for (int j = 0; j < matching_dof(n)(p).GetM(); j++)
                    MatchingDofOrig_Subdomain(p)(offset(p) + j) = matching_dof(n)(p)(j) + this->offset_dof_unknown(i);
                  
                  offset(p) += matching_dof(n)(p).GetM();
                }
          }
      }
    
    // constructing MatchingQuadOrig_Subdomain
    nb_points_quadrature_on_neighboring_interfaces = 0;
    for (int j = 0; j < this->mesh.GetNbBoundaryRef(); j++)
      {
        int jloc = local_number(j);
        int nb_points = this->mesh_num.GetNbPointsQuadratureBoundary(jloc);
        if (this->NeighboringConnectivity(jloc).GetM() > 0)
          {
            this->NeighboringConnectivity(jloc)(9) = this->nb_points_quadrature_on_neighboring_interfaces;
            this->nb_points_quadrature_on_neighboring_interfaces += nb_points;
          }
      }
    
    int nb_proc_interac = this->MatchingNumber_Subdomain.GetM();
    this->MatchingQuadOrig_Subdomain.Reallocate(nb_proc_interac);
    for (int i = 0; i < nb_proc_interac; i++)
      {
        int proc = this->MatchingNumber_Subdomain(i);
        int nb_dof = 0;
        for (int j = 0; j < this->mesh.GetNbBoundaryRef(); j++)
          {
            int jloc = local_number(j);
            if (this->NeighboringConnectivity(jloc).GetM() > 0)
              {
                int proc2 = this->NeighboringConnectivity(jloc)(2);
                if (proc2 == proc)
                  nb_dof += this->mesh_num.GetNbPointsQuadratureBoundary(jloc);
              }
          }
            
        this->MatchingQuadOrig_Subdomain(i).Reallocate(nb_dof);
        int offset = 0;
        for (int j = 0; j < this->mesh.GetNbBoundaryRef(); j++)
          {
            int jloc = local_number(j);
            if (this->NeighboringConnectivity(jloc).GetM() > 0)
              {
                int proc2 = this->NeighboringConnectivity(jloc)(2);
                if (proc2 == proc)
                  {
                    int rot = this->NeighboringConnectivity(jloc)(5);
                    int rot2 = this->mesh_num
		      .GetOppositeOrientationFace(rot, this->mesh.BoundaryRef(jloc));
                    
                    int nb_points_face = this->mesh_num.GetNbPointsQuadratureBoundary(jloc);
                    int rf = this->mesh_num.GetOrderQuadrature(jloc);
                    const Matrix<int>& FacesQuadRotation
                      = this->mesh_num.number_map
		      .GetRotationQuadraturePoints(rf, this->mesh.Boundary(jloc));
                    
                    int offset1 = this->NeighboringConnectivity(jloc)(9);
                    
                    for (int k = 0; k < nb_points_face; k++)
                      {
                        // matching numbers
                        int krot = FacesQuadRotation(rot2, k);
                        int num_quad1;                        
                        if (rank < proc2)
                          num_quad1 =  offset1 + k;
                        else
                          num_quad1 = offset1 + krot;
                        
                        this->MatchingQuadOrig_Subdomain(i)(offset + k) = num_quad1;
                      }
                    
                    offset += nb_points_face;
                  }
              }
          }
        
        // DISP(MatchingQuadOrig_Subdomain(i));
      }
    
      
    if (this->print_level >= 8)
      {
        /*sleep(rank);
        for (int i = 0; i < this->MatchingNumber_Subdomain.GetM(); i++)
          {
            DISP(i); DISP(this->MatchingNumber_Subdomain(i));
            DISP(this->MatchingQuadOrig_Subdomain(i));
            }
        */
        /* Mesh<Dimension> sub_mesh;
           VectReal_wp step(5); step.Fill(); Mlt(0.25, step);
           sub_mesh = mesh;
           sub_mesh.SubdivideMesh(step);
           sub_mesh.Write(string("subdiv") + to_str(rank) + ".mesh");    
        */
      }
    
#endif
  }
  
  

#ifdef SELDON_WITH_MPI
  template<class Dimension> template<class T>
  void DistributedProblem<Dimension>::
  ReduceDistributedVectorFace(Vector<T>& X, const MPI_Op& oper, int nb_u) const
  {    
    int Nvol = this->mesh.GetNbBoundary();

    // if only one processor, no communication needed, we leave method
    if (this->GetNbProcPerMode() <= 1)
      return;
    
    // sorting numbers in order to have the same order in all processors
    IVect global_number(this->mesh.GetNbBoundaryRef());
    IVect local_number(this->mesh.GetNbBoundaryRef());
    local_number.Fill();
    global_number.Fill(-1);
    int nb_process; MPI_Comm_size(this->comm_group_mode, &nb_process);
    IVect nb_face_per_proc(nb_process);
    nb_face_per_proc.Fill(0);
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
        if (this->NeighboringConnectivity(i).GetM() > 0)
          {
            global_number(i) = this->NeighboringConnectivity(i)(0);
            nb_face_per_proc(this->NeighboringConnectivity(i)(2))++;
          }
      }
    
    Sort(global_number, local_number);
    
    int nb_proc = 0; IVect IndexProc(nb_process);
    IndexProc.Fill(-1);
    for (int i = 0; i < nb_process; i++)
      if (nb_face_per_proc(i) > 0)
        IndexProc(i) = nb_proc++;
    
    // then filling the arrays
    IVect ProcNumber(nb_proc);
    Vector<IVect> DofNumber(nb_proc);
    nb_proc = 0;
    for (int i = 0; i < nb_process; i++)
      if (nb_face_per_proc(i) > 0)
        {
          ProcNumber(nb_proc) = i;
          DofNumber(nb_proc).Reallocate(nb_face_per_proc(i));
          nb_proc++;
        }
    
    nb_face_per_proc.Fill(0);
    for (int j = 0; j < this->mesh.GetNbBoundaryRef(); j++)
      {
        int jloc = local_number(j);
        if (this->NeighboringConnectivity(jloc).GetM() > 0)
          {
            int proc = this->NeighboringConnectivity(jloc)(2);
            int p = IndexProc(proc);
            DofNumber(p)(nb_face_per_proc(proc)) = jloc;
            nb_face_per_proc(proc)++;
          }
      }
        
    this->ReduceDistributedVector(X, oper, ProcNumber, DofNumber, Nvol, nb_u);
  }
#endif

  
  //! returns true if the finite element matrix is symmetric
  template<class Dimension>
  bool DistributedProblem<Dimension>::IsSymmetricProblem(bool eigen) const
  {
    if (this->mesh_num.GetFormulationForPeriodicCondition()
        == this->mesh_num.STRONG_PERIODIC)
      return false;
    
    for (int i = 0; i < this->mesh.GetNbPeriodicReferences(); i++)
      if (this->mesh.GetPeriodicityTypeReference(i) > 0)
        return false;
    
    if (!this->sipg_formulation)
      return false;
    
    return var_computation.IsSymmetricGlobalMatrix();
  }
  
  
  //! returns true if the mass matrix is symmetric (for eigenvalue problems)
  template<class Dimension>
  bool DistributedProblem<Dimension>::IsSymmetricMassMatrix() const
  {
    return true;
  }
  

  //! copies input parameters of another similar problem
  template<class Dimension>
  void DistributedProblem<Dimension>::CopyInputData(const VarProblem_Base& var_base)
  {
    const DistributedProblem<Dimension>& var
      = static_cast<const DistributedProblem<Dimension>& >(var_base);
    
    DistributedProblem_Base::CopyInputData(var);
    VarGeometryProblem<Dimension>::CopyInputData(var);
  }


  //! returns the dofs associated with element i
  template<class Dimension>
  IVect DistributedProblem<Dimension>::GetNodleElement(int i, int n) const
  {
    if (n == 0)
      return this->mesh_num.Element(i).GetNodle();

    return this->other_mesh_num(n-1)->Element(i).GetNodle();
  }
  
}

#define MONTJOIE_FILE_DISTRIBUTED_PROBLEM_CXX
#endif

