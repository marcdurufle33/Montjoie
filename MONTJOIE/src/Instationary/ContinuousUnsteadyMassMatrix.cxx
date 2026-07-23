#ifndef MONTJOIE_FILE_CONTINUOUS_UNSTEADY_MASS_MATRIX_CXX

namespace Montjoie
{
  
  /*********************
   * VirtualMassMatrix *
   *********************/
  
  
  VirtualMassMatrix::~VirtualMassMatrix()
  {
  }
  
  
  void VirtualMassMatrix::InitDefaultValues()
  {
    compute_Dh_default = false;
    compute_Sh_default = false;
    compute_invDh_default = true;
    compute_DhMinusdtSh_default = false;
    compute_invDhPlusdtSh_default = false;
    
    compute_Bh_default = false;
    compute_invBh_default = false;
    compute_BhMinusdtSh_default = false;
    compute_invBhPlusdtSh_default = false;
    compute_ShVec_default = false;    
  }


  //! which matrix do you need, Dh, Dh^{-1}, (Dh - dt/2 Sh), (Dh + dt/2 Sh)^{-1}, Sh ?
  void VirtualMassMatrix::
  FindMatricesToCompute(bool& compute_Dh, bool& compute_invDh,
			bool& compute_DhMinusdtSh, bool& compute_invDhPlusdtSh, bool& compute_Sh)
  {
    compute_Dh = compute_Dh_default;
    compute_invDh = true;
    compute_DhMinusdtSh = compute_DhMinusdtSh_default;
    compute_invDhPlusdtSh = compute_invDhPlusdtSh_default;
    compute_Sh = true;
    if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_ORDER2)
      {
	compute_Dh = true;
	compute_invDh = true;
	compute_DhMinusdtSh = true;
	compute_invDhPlusdtSh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_SYSTEM)
      {
	compute_Dh = compute_Dh_default;
        compute_invDh = compute_invDh_default;
        if (var_time.GetTimeSchemeOrder() == 4)
          {
            compute_Dh = true;
            compute_invDh = true;
            compute_Sh = true;
          }
	compute_DhMinusdtSh = true;
	compute_invDhPlusdtSh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::THETA_SCHEME)
      {
	compute_Dh = true;
	compute_invDh = compute_invDh_default;
	compute_DhMinusdtSh = true;
	compute_invDhPlusdtSh = compute_invDhPlusdtSh_default;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::IMPLICIT_SYMMETRIC_MULTISTEP)
      {
	compute_Dh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::PADE_SCHEME)
      {
	compute_Dh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LOCAL_IMPLICIT_RK)
      {
	compute_Dh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LOCAL_IMEX_RK)
      {
	compute_Dh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LINEAR_SDIRK)
      {
        compute_Dh = true;
      }
  }
  

  //! which matrix do you need, Bh, Bh^{-1}, (Bh - dt/2 ShV), (Bh + dt/2 ShV)^{-1}, ShV ?
  void VirtualMassMatrix
  ::FindMatricesToComputeVec(bool& compute_Bh, bool& compute_invBh,
			     bool& compute_BhMinusdtSh, bool& compute_invBhPlusdtSh,
			     bool& compute_ShVec)
  {
    compute_Bh = compute_Bh_default;
    compute_invBh = true;
    compute_BhMinusdtSh = compute_BhMinusdtSh_default;
    compute_invBhPlusdtSh = compute_invBhPlusdtSh_default;
    compute_ShVec = true;
    
    if (!var_time.FirstOrderScheme())
      {
        compute_Bh = false;
	compute_invBh = false;
        compute_BhMinusdtSh = false;
	compute_invBhPlusdtSh = false;
        compute_ShVec = false;
      }

    if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_SYSTEM)
      {
	compute_Bh = compute_Bh_default;
	compute_invBh = compute_invBh_default;
	compute_ShVec = compute_ShVec_default;
        if (var_time.GetTimeSchemeOrder() == 4)
          {
            compute_Bh = true;
            compute_invBh = true;
            compute_ShVec = true;
          }
	
	compute_invBh = true;
	compute_ShVec = true;
	compute_BhMinusdtSh = true;
	compute_invBhPlusdtSh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::PADE_SCHEME)
      {
	compute_Bh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LINEAR_SDIRK)
      {
        compute_Bh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LOCAL_IMPLICIT_RK)
      {
	compute_Bh = true;
      }
    else if (var_time.GetTimeSchemeType() == TimeSchemeEnum::LOCAL_IMEX_RK)
      {
	compute_Bh = true;
      }
  }


  //! overwrites Y by C_h^{-1} Y
  void VirtualMassMatrix::SolveOperatorCh(VectReal_wp& Y)
  {
    cout << "Not available for this matrix" << endl;
    abort();
  }


  /*************************************
   * ContinuousUnsteadyMassMatrix_Base *
   *************************************/
  
  
  void ContinuousUnsteadyMassMatrix_Base::InitDefaultValues()
  {
    Dh = NULL;
    DhMinusdtSh = NULL;
    DhPlusdtSh = NULL;
    Sh = NULL;
    Bh = NULL;
    BhMinusdtSh = NULL;
    BhPlusdtSh = NULL;
    ShVec = NULL;
  }

  
  ContinuousUnsteadyMassMatrix_Base::~ContinuousUnsteadyMassMatrix_Base()
  {
    Clear();
  }


  void ContinuousUnsteadyMassMatrix_Base::Clear()
  {
    if (Dh != NULL)
      delete Dh;
  
    if (DhMinusdtSh != NULL)
      delete DhMinusdtSh;

    if (DhPlusdtSh != NULL)
      delete DhPlusdtSh;

    if (Sh != NULL)
      delete Sh;

    if (Bh != NULL)
      delete Bh;
  
    if (BhMinusdtSh != NULL)
      delete BhMinusdtSh;

    if (BhPlusdtSh != NULL)
      delete BhPlusdtSh;

    if (ShVec != NULL)
      delete ShVec;
    
    InitDefaultValues();
  }

  
  //! applies both Dh and Bh
  void ContinuousUnsteadyMassMatrix_Base
  ::ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
		    const VectReal_wp& X, const Real_wp& beta, VectReal_wp& Y)
  {
    if (Dh == NULL)
      {
	cout << "Dh not allocated" << endl;
	abort();
      }
    
    int Ns = Dh->GetNumberOfRows();
    int Nv = Y.GetM() - Ns;
    VectReal_wp Uh, Vh, ProdUh, ProdVh;
    Uh.SetData(Ns, const_cast<Real_wp*>(&X(0)));
    ProdUh.SetData(Ns, &Y(0));
    if (Nv > 0)
      {
	Vh.SetData(Nv, const_cast<Real_wp*>(&X(Ns)));
	ProdVh.SetData(Nv, &Y(Ns));
      }
    
    ApplyOperatorDh(alpha, t, Uh, beta, ProdUh);
    if (Nv > 0)
      ApplyOperatorBh(alpha, t, Vh, beta, ProdVh);
    
    Uh.Nullify();
    Vh.Nullify();
    ProdUh.Nullify();
    ProdVh.Nullify();
  }
  

  //! solves both by Dh and Bh
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveMassMatrix(VectReal_wp& Y)
  {
    if (Dh == NULL)
      {
	cout << "Dh not allocated" << endl;
	abort();
      }
    
    int Ns = Dh->GetNumberOfRows();
    int Nv = Y.GetM() - Ns;
    VectReal_wp Uh, Vh;
    Uh.SetData(Ns, &Y(0));
    if (Nv > 0)
      Vh.SetData(Nv, &Y(Ns));
    
    SolveOperatorDh(Uh);
    if (Nv > 0)
      SolveOperatorBh(Vh);
    
    Uh.Nullify();
    Vh.Nullify();
  }
  
  
  void ContinuousUnsteadyMassMatrix_Base::InitDiagonalMass()
  { 
    if (Dh == NULL)
      Dh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (DhMinusdtSh == NULL)
      DhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (DhPlusdtSh == NULL)
      DhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (Sh == NULL)
      Sh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (Bh == NULL)
      Bh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (BhMinusdtSh == NULL)
      BhMinusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (BhPlusdtSh == NULL)
      BhPlusdtSh = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
    if (ShVec == NULL)
      ShVec = this->GetNewMassMatrix(FemMassMatrix::DIAGONAL);
    
  }
  
  
  //! computation of needed mass matrices Dh, Sh, Bh etc for the chosen time scheme
  void ContinuousUnsteadyMassMatrix_Base::Init(bool compute_time)
  {    
    int N = var_computation.GetNbRows();
    DistributedMatrix<Real_wp, General, ArrayRowSparse> mat_mass(N, N), mat_sigma(N, N);
    
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(1.0); nat_mat.SetCoefStiffness(0.0);
    nat_mat.SetCoefDamping(0.0);
    
    // computation of mass matrix
    var_boundary.SetCoefficientDirichlet(Real_wp(1));
    var_computation.AddMatrixWithBC(mat_mass, nat_mat);
    
    // computation of damping matrix
    nat_mat.SetCoefMass(0.0); nat_mat.SetCoefDamping(1.0);
    var_boundary.SetCoefficientDirichlet(Real_wp(0));
    var_computation.AddMatrixWithBC(mat_sigma, nat_mat);
    
    // mass matrix for main unknowns (continuous or in Hcurl)
    ComputeScalarMassMatrix(mat_mass, mat_sigma);
    
    // mass matrix for intermediary unknowns (in L2)
    if (var_problem.FirstOrderFormulation())
      ComputeVectorialMassMatrix(mat_mass, mat_sigma);
    
    // we clear matrices so that we have more space for factorization
    mat_mass.Clear(); mat_sigma.Clear();
  }
  

  //! returns the size used by the object in bytes
  size_t ContinuousUnsteadyMassMatrix_Base::GetMemorySize() const
  {
    size_t taille = NonLumpedElement.GetMemorySize() + NonLumpedDof.GetMemorySize() + InverseDof.GetMemorySize();
    
    if (Dh != NULL)
      taille += Dh->GetMemorySize();

    if (DhMinusdtSh != NULL)
      taille += DhMinusdtSh->GetMemorySize();
    
    if (DhPlusdtSh != NULL)
      taille += DhPlusdtSh->GetMemorySize();
    
    if (Sh != NULL)
      taille += Sh->GetMemorySize();

    if (Bh != NULL)
      taille += Bh->GetMemorySize();

    if (BhMinusdtSh != NULL)
      taille += BhMinusdtSh->GetMemorySize();
    
    if (BhPlusdtSh != NULL)
      taille += BhPlusdtSh->GetMemorySize();
    
    if (ShVec != NULL)
      taille += ShVec->GetMemorySize();

    return taille;
  }


  template<class Prop, class Storage>
  void ContinuousUnsteadyMassMatrix_Base
  ::ExtractScalarMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& A,
			int N, int nb_global_scalar_dof, Vector<int>& scalar_rows,
			DistributedMatrix<Real_wp, Prop, Storage>& B)
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);

    if (nb_proc > 1)
      {
	const IVect& all_rows = A.GetGlobalRowNumber();
	int nb_scalar_dof = A.GetNodlScalar();
	int nb_u = A.GetNbScalarUnknowns();
	scalar_rows.Reallocate(nb_scalar_dof);
	for (int i = 0; i < nb_scalar_dof; i++)
	  scalar_rows(i) = all_rows(i);

	B.Init(nb_global_scalar_dof*nb_u, &scalar_rows, const_cast<IVect*>(&A.GetOverlapRowNumber()),
	       const_cast<IVect*>(&A.GetOverlapProcNumber()), nb_scalar_dof, nb_u,
	       const_cast<IVect*>(&A.GetProcessorSharingRows()),
	       const_cast<Vector<IVect>* >(&A.GetSharingRowNumbers()), A.GetCommunicator());
      }

    GetSubMatrix(static_cast<const Matrix<Real_wp, General, ArrayRowSparse>& >(A),
		 0, N, 0, N, static_cast<Matrix<Real_wp, Prop, Storage>& >(B));
#else
    GetSubMatrix(A, 0, N, 0, N, B);
#endif
    
  }
  

  //! returns true if local part of mat_mass is symmetric
  bool ContinuousUnsteadyMassMatrix_Base
  ::IsLocalSymmetricMatrix(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
			   int N, Real_wp threshold)
  {
    Real_wp m_ref = MaxAbs(mat_mass);
    bool sym = true;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < mat_mass.GetRowSize(i); j++)
	if (abs(mat_mass.Value(i, j) - mat_mass(mat_mass.Index(i, j), i)) > threshold*m_ref)
	  sym = false;

#ifdef SELDON_WITH_MPI
    bool sym_tmp(sym);
    MPI_Allreduce(&sym_tmp, &sym, 1, MPI_CHAR, MPI_LAND, mat_mass.GetCommunicator());
#endif
	
    return sym;
  }


  //! constructs arrays InverseDof, NonLumpedElement and NonLumpedDof
  void ContinuousUnsteadyMassMatrix_Base::ConstructNonLumpedDofs(const Vector<bool>& diag_elt)
  {
    int Nvol = var_problem.offset_dof_unknown(1);
    Vector<int> Nodle;
    InverseDof.Reallocate(Nvol);
    InverseDof.Fill(-1);
    int ne = 0;
    for (int i = 0; i < diag_elt.GetM(); i++)
      if (!diag_elt(i))
	{
	  Nodle = var_comm.GetScalarDofNumberOnElement(i);
	  for (int j = 0; j < Nodle.GetM(); j++)
	    if (Nodle(j) >= 0)
	      if (InverseDof(Nodle(j)) == -1)
		InverseDof(Nodle(j)) = 0;
	  
	  ne++;
	}

#ifdef SELDON_WITH_MPI
    var_comm.ReduceDistributedVector(InverseDof, MPI_MAX, 1);
#endif

    int ndof = 0;
    for (int i = 0; i < Nvol; i++)
      if (InverseDof(i) == 0)
	ndof++;
    
    NonLumpedElement.Reallocate(ne);
    NonLumpedDof.Reallocate(ndof);
    ne = 0;
    for (int i = 0; i < diag_elt.GetM(); i++)
      if (!diag_elt(i))
	{
	  NonLumpedElement(ne) = i;
	  ne++;
	}
    
    ndof = 0;
    for (int i = 0; i < Nvol; i++)
      if (InverseDof(i) >= 0)
	{
	  NonLumpedDof(ndof) = i;
	  InverseDof(i) = ndof;
	  ndof++;
	}
  }


  //! construct lumped dofs (arrays LumpedDof and IndexDof)
  void ContinuousUnsteadyMassMatrix_Base::ConstructLumpedDofs(Vector<int>& IndexDof)
  {
    int Nvol = InverseDof.GetM();
    int nb_dof_sparse = NonLumpedDof.GetM();
    int nb_dof_diag = Nvol - nb_dof_sparse;
    LumpedDof.Reallocate(nb_dof_diag);
    IndexDof.Reallocate(Nvol); IndexDof.Fill(-1);
    nb_dof_diag = 0;
    for (int i = 0; i < Nvol; i++)
      if (InverseDof(i) == -1)
	{
	  LumpedDof(nb_dof_diag) = i;
	  IndexDof(i) = nb_dof_diag;
	  nb_dof_diag++;
	}
  }

  
  //! extracts and assembles diagonal part of mat_mass and mat_sigma
  void ContinuousUnsteadyMassMatrix_Base
  ::ExtractDiagonal(const DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
		    const DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_sigma,
		    int N, VectReal_wp& Dh_diagonal, VectReal_wp& Sh_diagonal,
		    VectReal_wp& Dh_ass, VectReal_wp& Sh_ass)
  {
    Dh_diagonal.Reallocate(N);
    Sh_diagonal.Reallocate(N);
    for (int i = 0; i < N; i++)
      {
	Dh_diagonal(i) = mat_mass(i, i);
	Sh_diagonal(i) = mat_sigma(i, i);
      }
    
    Dh_ass = Dh_diagonal;
    Sh_ass = Sh_diagonal;
    var_comm.AddDomains(Dh_ass);
    var_comm.AddDomains(Sh_ass);
  }
  

  //! fills diagonal matrices
  void ContinuousUnsteadyMassMatrix_Base
  ::FillDiagonal(const VectReal_wp& Dh_diagonal, const VectReal_wp& Sh_diagonal,
		 const VectReal_wp& Dh_ass, const VectReal_wp& Sh_ass,
		 bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
		 bool compute_invDhPlusdtSh, bool compute_Sh)
  {
    Real_wp dt = var_time.GetTimeStep();
    int N = Dh_diagonal.GetM();
    
    // assembling the needed matrices
    if (compute_invDh)
      {
	VectReal_wp& invDiagonal = Dh->GetInverseDiagonal();
	invDiagonal.Reallocate(N);
	for (int i = 0; i < N; i++)
	  invDiagonal(i) = Real_wp(1)/Dh_ass(i);
      }
    
    if (compute_DhMinusdtSh)
      {
	VectReal_wp& diagonal = DhMinusdtSh->GetDiagonal();
	diagonal.Reallocate(N);
	for (int i = 0; i < N; i++)
	  diagonal(i) = Dh_diagonal(i) - 0.5*dt*Sh_diagonal(i);
      }
    
    if (compute_invDhPlusdtSh)
      {
	VectReal_wp& diagonal = DhPlusdtSh->GetInverseDiagonal();
	diagonal.Reallocate(N);
	for (int i = 0; i < N; i++)
	  diagonal(i) = Real_wp(1)/(Dh_ass(i) + 0.5*dt*Sh_ass(i));
      }
    
    if (compute_Dh)
      Dh->GetDiagonal() = Dh_diagonal;
    
    if (compute_Sh)
      Sh->GetDiagonal() = Sh_diagonal;
  }


  //! fills block-diagonal matrices
  void ContinuousUnsteadyMassMatrix_Base::
  FillBlockDiagonal(const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& Dh_block,
		    const DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>& Sh_block,
		    bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
		    bool compute_invDhPlusdtSh, bool compute_Sh, bool invert)
  {
    typedef DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> TypeBlockDiag;
    Real_wp dt = var_time.GetTimeStep();
	
    if (compute_invDh)
      {
	TypeBlockDiag& block_diagonal = Dh->GetBlockDiagonal();
	block_diagonal = Dh_block;
	if (invert)
	  Dh->Invert(compute_Dh);	      
      }
    
    if (compute_DhMinusdtSh)
      {
	TypeBlockDiag& block_diagonal = DhMinusdtSh->GetBlockDiagonal();
	block_diagonal = Dh_block;
	Add(-0.5*dt, Sh_block, block_diagonal);
      }
    
    if (compute_invDhPlusdtSh)
      {
	TypeBlockDiag& block_diagonal = DhPlusdtSh->GetBlockDiagonal();
	block_diagonal = Dh_block;
	Add(0.5*dt, Sh_block, block_diagonal);
		
	if (invert)
	  DhPlusdtSh->Invert(false);
      }
    
    if (compute_Sh)
      Sh->GetBlockDiagonal() = Sh_block;
  }


  //! fills block-diagonal matrices
  void ContinuousUnsteadyMassMatrix_Base::
  FillBlockDiagonal(const DistributedMatrix<Real_wp, General, BlockDiagRow>& Dh_block,
		    const DistributedMatrix<Real_wp, General, BlockDiagRow>& Sh_block,
		    bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
		    bool compute_invDhPlusdtSh, bool compute_Sh, bool invert)
  {
    typedef DistributedMatrix<Real_wp, General, BlockDiagRow> TypeBlockDiag;
    Real_wp dt = var_time.GetTimeStep();
    
    // assembling the needed matrices
    if (compute_invDh)
      {
	TypeBlockDiag& block_diagonal = Dh->GetUnsymmetricBlockDiagonal();
	block_diagonal = Dh_block;
	if (invert)
	  Dh->Invert(compute_Dh);
      }
    
    if (compute_DhMinusdtSh)
      {
	TypeBlockDiag& block_diagonal = DhMinusdtSh->GetUnsymmetricBlockDiagonal();
	block_diagonal = Dh_block;
	Add(-0.5*dt, Sh_block, block_diagonal);
      }
    
    if (compute_invDhPlusdtSh)
      {
	TypeBlockDiag& block_diagonal = DhPlusdtSh->GetUnsymmetricBlockDiagonal();
	block_diagonal = Dh_block;
	Add(0.5*dt, Sh_block, block_diagonal);

	if (invert)
	  DhPlusdtSh->Invert(false);
      }
    
    if (compute_Sh)
      Sh->GetUnsymmetricBlockDiagonal() = Sh_block;
  }
  

  //! Fills sparse matrices
  void ContinuousUnsteadyMassMatrix_Base
  ::FillSparseMatrices(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Dh_s,
		       DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Sh_s,
		       bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
		       bool compute_invDhPlusdtSh, bool compute_Sh)
  {
    const MassMatrixParameter& param = var_time.GetMassMatrixParameter();
    bool it_solver = param.UseIterativeSolver();
    int nb_iter = param.GetNbMaxIterations();
    Real_wp eps = param.GetStoppingCriterion();
    int print_level = param.GetPrintLevel();
    Real_wp dt = var_time.GetTimeStep();    

    if (compute_DhMinusdtSh)
      {
	DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& M = DhMinusdtSh->GetSymSparseMatrix();
	M = Dh_s;
	Add(-Real_wp(0.5)*dt, Sh_s, M);
      }
    
    if (compute_invDhPlusdtSh)
      {
	DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& M = DhPlusdtSh->GetSymSparseMatrix();
	M = Dh_s;
	Add(Real_wp(0.5)*dt, Sh_s, M);
	
	DhPlusdtSh->SetIterativeSolver(it_solver, nb_iter, eps, print_level);
	DhPlusdtSh->Invert(false);
      }
    
    if (compute_Dh || compute_invDh)
      {
	if (compute_invDh)
	  {
	    Dh->SetIterativeSolver(it_solver, nb_iter, eps, print_level);
	    Dh->Invert(compute_Dh);
	  }
      }
    else
      Dh_s.Clear();
    
    if (!compute_Sh)
      Sh_s.Clear();
  }
  
  
  //! Fills sparse matrices
  void ContinuousUnsteadyMassMatrix_Base
  ::FillSparseMatrices(DistributedMatrix<Real_wp, General, ArrayRowSparse>& Dh_s,
		       DistributedMatrix<Real_wp, General, ArrayRowSparse>& Sh_s,
		       bool compute_Dh, bool compute_invDh, bool compute_DhMinusdtSh,
		       bool compute_invDhPlusdtSh, bool compute_Sh)
  {
    const MassMatrixParameter& param = var_time.GetMassMatrixParameter();
    bool it_solver = param.UseIterativeSolver();
    int nb_iter = param.GetNbMaxIterations();
    Real_wp eps = param.GetStoppingCriterion();
    int print_level = param.GetPrintLevel();
    Real_wp dt = var_time.GetTimeStep();
    
    if (compute_DhMinusdtSh)
      {
	DistributedMatrix<Real_wp, General, ArrayRowSparse>& M = DhMinusdtSh->GetSparseMatrix();
	M = Dh_s;
	Add(-Real_wp(0.5)*dt, Sh_s, M);
      }
    
    if (compute_invDhPlusdtSh)
      {
	DistributedMatrix<Real_wp, General, ArrayRowSparse>& M = DhPlusdtSh->GetSparseMatrix();
	M = Dh_s;
	Add(Real_wp(0.5)*dt, Sh_s, M);
	
	DhPlusdtSh->SetIterativeSolver(it_solver, nb_iter, eps, print_level);
	DhPlusdtSh->Invert(false);
      }
    
    if (compute_Dh || compute_invDh)
      {
	if (compute_invDh)
	  {
	    Dh->SetIterativeSolver(it_solver, nb_iter, eps, print_level);
	    Dh->Invert(compute_Dh);
	  }
      }
    else
      Dh_s.Clear();
    
    if (!compute_Sh)
      Sh_s.Clear();
  }


  //! symmetrizes the pattern of mat_mass
  void ContinuousUnsteadyMassMatrix_Base
  ::SymmetrizePattern(DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass)
  {
    for (int i = 0; i < mat_mass.GetM(); i++)
      for (int j = 0; j < mat_mass.GetRowSize(i); j++)
	if (mat_mass(mat_mass.Index(i, j), i) == 0)
	  mat_mass.Get(mat_mass.Index(i, j), i) = 0.0;
  }
  

  //! computation of needed scalar mass matrices (for continuous unknown)
  void ContinuousUnsteadyMassMatrix_Base::
  ComputeScalarMassMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
                          DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_sigma)
  {
    const MassMatrixParameter& param = var_time.GetMassMatrixParameter();
    bool store_mass = param.MassMatrixStored();
    
    // H1 or Hcurl formulation
    // finding the type of matrix
    Vector<bool> diag_elt;
    int type_matrix = var_computation.GetMassMatrixType(diag_elt);

    if ((type_matrix == FemMassMatrix::MATRIX_FREE) && store_mass)
      type_matrix = FemMassMatrix::MATRIX_SPARSE_UNSYM;
    
    if ((type_matrix == FemMassMatrix::DIAG_MATRIX_FREE) && store_mass)
      type_matrix = FemMassMatrix::DIAG_SPARSE_UNSYM;

    if ((type_matrix == FemMassMatrix::BLOCK_DIAG_MATRIX_FREE) && store_mass)
      type_matrix = FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM;

    int N = var_problem.offset_dof_unknown(var_problem.nb_unknowns_scal);
    int nb_u = var_problem.nb_unknowns_scal;
    bool matrix_sym = false;
    if (IsLocalSymmetricMatrix(mat_mass, N) && IsLocalSymmetricMatrix(mat_sigma, N))
      matrix_sym = true;

    //matrix_sym = false;
    if (matrix_sym)
      type_matrix = FemMassMatrix::GetSymmetryType(type_matrix);

    // do we have to compte Dh, Sh, Dh^{-1}, etc ?
    bool compute_Dh, compute_invDh, compute_DhMinusdtSh,
      compute_invDhPlusdtSh, compute_Sh;
    
    this->FindMatricesToCompute(compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh);
    
    // allocating matrices
    if (compute_invDh || compute_Dh)
      Dh = GetNewMassMatrix(type_matrix);
    
    if (compute_DhMinusdtSh)
      DhMinusdtSh = GetNewMassMatrix(type_matrix);

    if (compute_invDhPlusdtSh)
      DhPlusdtSh = GetNewMassMatrix(type_matrix);

    if (compute_Sh)
      Sh = GetNewMassMatrix(type_matrix);

    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // diagonal mass matrices
	  VectReal_wp Dh_diagonal, Sh_diagonal, Dh_ass, Sh_ass;
	  ExtractDiagonal(mat_mass, mat_sigma, N, Dh_diagonal, Sh_diagonal, Dh_ass, Sh_ass);
	  
	  FillDiagonal(Dh_diagonal, Sh_diagonal, Dh_ass, Sh_ass,
		       compute_Dh, compute_invDh, compute_DhMinusdtSh,
		       compute_invDhPlusdtSh, compute_Sh);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // symmetric block-diagonal mass matrices
	  typedef DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> TypeBlockDiag;
	  TypeBlockDiag Dh_block, Sh_block;
          
	  ConvertToBlockDiagonal(mat_mass, Dh_block, 0, N);
	  ConvertToBlockDiagonal(mat_sigma, Sh_block, 0, N);

	  FillBlockDiagonal(Dh_block, Sh_block,
			    compute_Dh, compute_invDh, compute_DhMinusdtSh,
			    compute_invDhPlusdtSh, compute_Sh);
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE :
	{
	  // symmetric sparse mass matrices
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Dh_s = Dh->GetSymSparseMatrix();
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Sh_s = Sh->GetSymSparseMatrix();

	  if (var_problem.FirstOrderFormulation())
	    {
              int Nglob = var_comm.GetOffsetGlobalUnknownDof(var_problem.nb_unknowns_scal);
	      ExtractScalarMatrix(mat_mass, N, Nglob, scalar_rows_, Dh_s);
	      ExtractScalarMatrix(mat_sigma, N, Nglob, scalar_rows_, Sh_s);
	    }
	  else
	    {
	      Copy(mat_mass, Dh_s);
	      Copy(mat_sigma, Sh_s);
	    }

	  FillSparseMatrices(Dh_s, Sh_s,
			     compute_Dh, compute_invDh, compute_DhMinusdtSh,
			     compute_invDhPlusdtSh, compute_Sh);
        }
	break;
      case FemMassMatrix::BLOCK_DIAGONAL_UNSYM :
	{
	  // block-diagonal mass matrices (unsymmetric)
	  typedef DistributedMatrix<Real_wp, General, BlockDiagRow> TypeBlockDiag;
	  TypeBlockDiag Dh_block, Sh_block;

          // the pattern of mat_mass and mat_sigma is symmetrized
	  SymmetrizePattern(mat_mass);
	  SymmetrizePattern(mat_sigma);

	  // matrices are converted to block-diagonal
	  ConvertToBlockDiagonal(mat_mass, Dh_block, 0, N);
	  ConvertToBlockDiagonal(mat_sigma, Sh_block, 0, N);

	  // then stored in Dh objects
	  FillBlockDiagonal(Dh_block, Sh_block,
			    compute_Dh, compute_invDh, compute_DhMinusdtSh,
			    compute_invDhPlusdtSh, compute_Sh);
	  
	}
	break;
      case FemMassMatrix::MATRIX_SPARSE_UNSYM :
	{
	  // sparse mass matrices (unsymmetric)
	  DistributedMatrix<Real_wp, General, ArrayRowSparse>& Dh_s = Dh->GetSparseMatrix();
	  DistributedMatrix<Real_wp, General, ArrayRowSparse>& Sh_s = Sh->GetSparseMatrix();

	  if (var_problem.FirstOrderFormulation())
	    {
              int Nglob = var_comm.GetOffsetGlobalUnknownDof(var_problem.nb_unknowns_scal);
	      ExtractScalarMatrix(mat_mass, N, Nglob, scalar_rows_, Dh_s);
	      ExtractScalarMatrix(mat_sigma, N, Nglob, scalar_rows_, Sh_s);
	    }
	  else
	    {
	      Dh_s = mat_mass;
	      Sh_s = mat_sigma;
	    }

	  FillSparseMatrices(Dh_s, Sh_s,
			     compute_Dh, compute_invDh, compute_DhMinusdtSh,
			     compute_invDhPlusdtSh, compute_Sh);
        }
	break;
      case FemMassMatrix::DIAG_SPARSE :
	{
	  ConstructNonLumpedDofs(diag_elt);

	  // diagonal part is extracted
	  VectReal_wp Dh_diagonal, Sh_diagonal, Dh_assembled, Sh_assembled;
	  ExtractDiagonal(mat_mass, mat_sigma, N, Dh_diagonal, Sh_diagonal, Dh_assembled, Sh_assembled);

	  VectReal_wp Dh_diag, Sh_diag, Dh_ass, Sh_ass;
          int Nvol = var_problem.offset_dof_unknown(1);
	  int nb_dof_sparse = NonLumpedDof.GetM();
	  int nb_dof_diag = Nvol - nb_dof_sparse;
	  Dh_diag.Reallocate(nb_u*nb_dof_diag);
	  Sh_diag.Reallocate(nb_u*nb_dof_diag);
	  Dh_ass.Reallocate(nb_u*nb_dof_diag);
	  Sh_ass.Reallocate(nb_u*nb_dof_diag);
	  int num = 0;
	  for (int i = 0; i < Nvol; i++)
	    if (InverseDof(i) == -1)
	      {
		for (int k = 0; k < nb_u; k++)
		  {
		    int k_glob = k*Nvol + i;
		    int k_loc = k*nb_dof_diag + num;
		    Dh_diag(k_loc) = Dh_diagonal(k_glob);
		    Sh_diag(k_loc) = Sh_diagonal(k_glob);
		    Dh_ass(k_loc) = Dh_assembled(k_glob);
		    Sh_ass(k_loc) = Sh_assembled(k_glob);
		  }
		
		num++;
	      }
	  
	  Dh_diagonal.Clear(); Sh_diagonal.Clear(); Dh_assembled.Clear(); Sh_assembled.Clear();
	  FillDiagonal(Dh_diag, Sh_diag, Dh_ass, Sh_ass,
		       compute_Dh, compute_invDh, compute_DhMinusdtSh,
		       compute_invDhPlusdtSh, compute_Sh);
	  
	  // sparse matrices are extracted
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Dh_s = Dh->GetSymSparseMatrix();
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Sh_s = Sh->GetSymSparseMatrix();

	  if (nb_u == 1)
	    {
#ifdef SELDON_WITH_MPI
	      ExtractDistributedInfo(mat_mass, NonLumpedDof, InverseDof, info_sparse);
#endif
	      
	      ExtractSubMatrix(mat_mass, NonLumpedDof, InverseDof,
			       info_sparse, Dh_s);
	      
	      ExtractSubMatrix(mat_sigma, NonLumpedDof, InverseDof,
			       info_sparse, Sh_s);
	    }
	  else
	    {
	      cout << "not implemented" << endl;
	      abort();
	    }

	  FillSparseMatrices(Dh_s, Sh_s,
			     compute_Dh, compute_invDh, compute_DhMinusdtSh,
			     compute_invDhPlusdtSh, compute_Sh);
        }
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE :
	{
	  ConstructNonLumpedDofs(diag_elt);
	  
	  Vector<int> IndexDof;
	  ConstructLumpedDofs(IndexDof);
	  int nb_dof_diag = LumpedDof.GetM();
	  
	  // block-diagonal part is extracted
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Dh_b;
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> Sh_b;
	  
	  if (nb_u == 1)
	    {
#ifdef SELDON_WITH_MPI
	      ExtractDistributedInfo(mat_mass, LumpedDof, IndexDof, info_block);
#endif
	      
	      ExtractSubMatrix(mat_mass, LumpedDof, IndexDof, info_block, Dh_b);
	      ExtractSubMatrix(mat_sigma, LumpedDof, IndexDof, info_block, Sh_b);
	    }
	  else
	    {
	      cout << "not implemented" << endl;
	      abort();
	    }
	  
	  // matrices are converted to block-diagonal
	  typedef DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> TypeBlockDiag;
	  TypeBlockDiag Dh_block, Sh_block;
	  
	  if (LumpedDof.GetM() > 0)
	    {
	      ConvertToBlockDiagonal(Dh_b, Dh_block, 0, nb_dof_diag);
	      ConvertToBlockDiagonal(Sh_b, Sh_block, 0, nb_dof_diag);
	      
	      // then stored in Dh objects
	      FillBlockDiagonal(Dh_block, Sh_block,
				compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh, false);
	    }
	  
	  Dh_block.Clear(); Sh_block.Clear();
	  Dh_b.Clear(); Sh_b.Clear();
	  
	  // sparse matrices are extracted
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Dh_s = Dh->GetSymSparseMatrix();
	  DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& Sh_s = Sh->GetSymSparseMatrix();
	  
	  if (nb_u == 1)
	    {
#ifdef SELDON_WITH_MPI
	      ExtractDistributedInfo(mat_mass, NonLumpedDof, InverseDof, info_sparse);
#endif
	      
	      ExtractSubMatrix(mat_mass, NonLumpedDof, InverseDof,
			       info_sparse, Dh_s);
	      
	      ExtractSubMatrix(mat_sigma, NonLumpedDof, InverseDof,
			       info_sparse, Sh_s);
	    }
	  else
	    {
	      cout << "not implemented" << endl;
	      abort();
	    }
	  
	  FillSparseMatrices(Dh_s, Sh_s,
			     compute_Dh, compute_invDh, compute_DhMinusdtSh,
			     compute_invDhPlusdtSh, compute_Sh);
        }
	break;
      case FemMassMatrix::DIAG_SPARSE_UNSYM :
	{
	  ConstructNonLumpedDofs(diag_elt);

	  // diagonal part is extracted
	  VectReal_wp Dh_diagonal, Sh_diagonal, Dh_assembled, Sh_assembled;
	  ExtractDiagonal(mat_mass, mat_sigma, N, Dh_diagonal, Sh_diagonal, Dh_assembled, Sh_assembled);

	  VectReal_wp Dh_diag, Sh_diag, Dh_ass, Sh_ass;
          int Nvol = var_problem.offset_dof_unknown(1);
	  int nb_dof_sparse = NonLumpedDof.GetM();
	  int nb_dof_diag = Nvol - nb_dof_sparse;
	  Dh_diag.Reallocate(nb_u*nb_dof_diag);
	  Sh_diag.Reallocate(nb_u*nb_dof_diag);
	  Dh_ass.Reallocate(nb_u*nb_dof_diag);
	  Sh_ass.Reallocate(nb_u*nb_dof_diag);
	  int num = 0;
	  for (int i = 0; i < Nvol; i++)
	    if (InverseDof(i) == -1)
	      {
		for (int k = 0; k < nb_u; k++)
		  {
		    int k_glob = k*Nvol + i;
		    int k_loc = k*nb_dof_diag + num;
		    Dh_diag(k_loc) = Dh_diagonal(k_glob);
		    Sh_diag(k_loc) = Sh_diagonal(k_glob);
		    Dh_ass(k_loc) = Dh_assembled(k_glob);
		    Sh_ass(k_loc) = Sh_assembled(k_glob);
		  }
		
		num++;
	      }

	  Dh_diagonal.Clear(); Sh_diagonal.Clear(); Dh_assembled.Clear(); Sh_assembled.Clear();
	  FillDiagonal(Dh_diag, Sh_diag, Dh_ass, Sh_ass,
		       compute_Dh, compute_invDh, compute_DhMinusdtSh,
		       compute_invDhPlusdtSh, compute_Sh);
	  
	  // sparse matrices are extracted
	  DistributedMatrix<Real_wp, General, ArrayRowSparse>& Dh_s = Dh->GetSparseMatrix();
	  DistributedMatrix<Real_wp, General, ArrayRowSparse>& Sh_s = Sh->GetSparseMatrix();
	  
	  if (nb_u == 1)
	    {
#ifdef SELDON_WITH_MPI
	      ExtractDistributedInfo(mat_mass, NonLumpedDof, InverseDof, info_sparse);
#endif
	      
	      ExtractSubMatrix(mat_mass, NonLumpedDof, InverseDof,
			       info_sparse, Dh_s);
	      
	      ExtractSubMatrix(mat_sigma, NonLumpedDof, InverseDof,
			       info_sparse, Sh_s);
	    }
	  else
	    {
	      cout << "not implemented" << endl;
	      abort();
	    }
	  
	  FillSparseMatrices(Dh_s, Sh_s,
			     compute_Dh, compute_invDh, compute_DhMinusdtSh,
			     compute_invDhPlusdtSh, compute_Sh);
        }
	break;
      case FemMassMatrix::BLOCK_DIAG_SPARSE_UNSYM :
	{
	  ConstructNonLumpedDofs(diag_elt);

	  Vector<int> IndexDof;
	  ConstructLumpedDofs(IndexDof);
	  int nb_dof_diag = LumpedDof.GetM();

	  // block-diagonal part is extracted
	  DistributedMatrix<Real_wp, General, ArrayRowSparse> Dh_b;
	  DistributedMatrix<Real_wp, General, ArrayRowSparse> Sh_b;

	  if (nb_u == 1)
	    {
#ifdef SELDON_WITH_MPI
	      ExtractDistributedInfo(mat_mass, LumpedDof, IndexDof, info_block);
#endif
	      
	      ExtractSubMatrix(mat_mass, LumpedDof, IndexDof, info_block, Dh_b);
	      ExtractSubMatrix(mat_sigma, LumpedDof, IndexDof, info_block, Sh_b);
	    }
	  else
	    {
	      cout << "not implemented" << endl;
	      abort();
	    }
	  
	  // the pattern of Dh_b and Sh_b is symmetrized
	  SymmetrizePattern(Dh_b);
	  SymmetrizePattern(Sh_b);
	  
	  // matrices are converted to block-diagonal
	  typedef DistributedMatrix<Real_wp, General, BlockDiagRow> TypeBlockDiag;
	  TypeBlockDiag Dh_block, Sh_block;
	  
	  if (nb_dof_diag > 0)
	    {
	      ConvertToBlockDiagonal(Dh_b, Dh_block, 0, nb_dof_diag);
	      ConvertToBlockDiagonal(Sh_b, Sh_block, 0, nb_dof_diag);

	      // then stored in Dh objects
	      FillBlockDiagonal(Dh_block, Sh_block,
				compute_Dh, compute_invDh, compute_DhMinusdtSh,
				compute_invDhPlusdtSh, compute_Sh, false);
	    }
	  
	  Dh_block.Clear(); Sh_block.Clear();
	  Dh_b.Clear(); Sh_b.Clear();
	  
	  // sparse matrices are extracted
	  DistributedMatrix<Real_wp, General, ArrayRowSparse>& Dh_s = Dh->GetSparseMatrix();
	  DistributedMatrix<Real_wp, General, ArrayRowSparse>& Sh_s = Sh->GetSparseMatrix();
	  
	  if (nb_u == 1)
	    {
#ifdef SELDON_WITH_MPI
	      ExtractDistributedInfo(mat_mass, NonLumpedDof, InverseDof, info_sparse);
#endif
	      
	      ExtractSubMatrix(mat_mass, NonLumpedDof, InverseDof,
			       info_sparse, Dh_s);
	      
	      ExtractSubMatrix(mat_sigma, NonLumpedDof, InverseDof,
			       info_sparse, Sh_s);
	    }
	  else
	    {
	      cout << "not implemented" << endl;
	      abort();
	    }

	  FillSparseMatrices(Dh_s, Sh_s,
			     compute_Dh, compute_invDh, compute_DhMinusdtSh,
			     compute_invDhPlusdtSh, compute_Sh);
        }
	break;
      default:
	{
	  cout << "Storage not implemented" << endl;
          abort();
        }
	break;
      }
  }
  

  //! computation of needed vectorial mass matrices (for discontinuous unknown)
  void ContinuousUnsteadyMassMatrix_Base
  ::ComputeVectorialMassMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_mass,
                               DistributedMatrix<Real_wp, General, ArrayRowSparse>& mat_sigma)
  {
    // H1 or Hcurl formulation
    // finding the type of matrix
    int type_matrix = FemMassMatrix::MATRIX_FREE;
    
    // mass matrices related to vectorial unknown have this type
    // do we have to compte Bh, Sh, Bh^{-1}, etc ?
    bool compute_Bh, compute_invBh, compute_BhMinusdtSh,
      compute_invBhPlusdtSh, compute_ShVec;
    
    this->FindMatricesToComputeVec(compute_Bh, compute_invBh, compute_BhMinusdtSh,
				   compute_invBhPlusdtSh, compute_ShVec);
    
    int Nscal = var_time.GetNbScalarUnknowns();
    int N = var_time.GetNbVectorialUnknowns();
    Real_wp dt = var_time.GetTimeStep();
    
    // allocating matrices
    if ((compute_invBh) || (compute_Bh))
      Bh = GetNewMassMatrix(type_matrix);
    
    if (compute_BhMinusdtSh)
      BhMinusdtSh = GetNewMassMatrix(type_matrix);

    if (compute_invBhPlusdtSh)
      BhPlusdtSh = GetNewMassMatrix(type_matrix);

    if (compute_ShVec)
      ShVec = GetNewMassMatrix(type_matrix);
    
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  VectReal_wp Bh_diagonal, ShVec_diagonal;
	  Bh_diagonal.Reallocate(N);
	  ShVec_diagonal.Reallocate(N);
	  for (int i = 0; i < N; i++)
	    {
	      Bh_diagonal(i) = -mat_mass(Nscal+i, Nscal+i);
	      ShVec_diagonal(i) = -mat_sigma(Nscal+i, Nscal+i);
	    }
	  
	  if (compute_invBh)
	    {
	      VectReal_wp& invDiagonal = Bh->GetInverseDiagonal();
	      invDiagonal = Bh_diagonal;
	      for (int i = 0; i < N; i++)
                invDiagonal(i) = 1.0/invDiagonal(i);
	    }
	  
	  if (compute_BhMinusdtSh)
	    {
	      VectReal_wp& diagonal = BhMinusdtSh->GetDiagonal();
	      diagonal.Reallocate(N);
	      for (int i = 0; i < N; i++)
		diagonal(i) = Bh_diagonal(i) - 0.5*dt*ShVec_diagonal(i);
	    }
	  
	  if (compute_invBhPlusdtSh)
	    {
	      VectReal_wp& diagonal = BhPlusdtSh->GetInverseDiagonal();
	      diagonal.Reallocate(N);
	      for (int i = 0; i < N; i++)
		diagonal(i) = Bh_diagonal(i) + 0.5*dt*ShVec_diagonal(i);
	      
	      for (int i = 0; i < N; i++)
		diagonal(i) = Real_wp(1)/diagonal(i);
	    }
	  
	  if (compute_Bh)
            Bh->GetDiagonal() = Bh_diagonal;
          
	  if (!compute_ShVec)
            ShVec->GetDiagonal() = ShVec_diagonal;
          
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  typedef DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> TypeBlockDiag;
	  TypeBlockDiag Bh_block, ShVec_block;

	  ConvertToBlockDiagonal(mat_mass, Bh_block, Nscal, Nscal+N);
	  ConvertToBlockDiagonal(mat_sigma, ShVec_block, Nscal, Nscal+N);
	  Mlt(Real_wp(-1), Bh_block);
	  Mlt(Real_wp(-1), ShVec_block);
	  
	  if (compute_invBh)
	    {
	      TypeBlockDiag& block_diagonal = Bh->GetBlockDiagonal();
	      block_diagonal = Bh_block;
	      Bh->Invert(compute_Bh);
	    }
	  
	  if (compute_BhMinusdtSh)
	    {
	      TypeBlockDiag& block_diagonal = BhMinusdtSh->GetBlockDiagonal();
	      block_diagonal = Bh_block;
	      Add(-0.5*dt, ShVec_block, block_diagonal);
            }
	  
	  if (compute_invBhPlusdtSh)
	    {
	      TypeBlockDiag& block_diagonal = BhPlusdtSh->GetBlockDiagonal();
	      block_diagonal = Bh_block;
	      Add(0.5*dt, ShVec_block, block_diagonal);
	      
	      BhPlusdtSh->Invert(false);
	    }
	  
	  if (compute_ShVec)
            ShVec->GetBlockDiagonal() = ShVec_block;
          
	}
	break;
      case FemMassMatrix::MATRIX_FREE:
	{
	  // matrix-free algorithm the computation is done in the leaf class HyperbolicProblem
	}
	break;
      default:
	{
	  cout << "Storage not implemented" << endl;
          abort();
        }
	break;
      }
  }

  
  //! applies operator Sh, Prod_Uh = beta Prod_Uh + alpha Sh Uh
  void ContinuousUnsteadyMassMatrix_Base::
  ApplyOperatorSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
		  const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    if (Sh == NULL)
      {
	cout << "Matrix Sh not computed" << endl;
	abort();
      }
    
    if (beta == Real_wp(0))
      Prod_Uh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);
    
    Sh->MltMass(alpha, Uh, Prod_Uh);
  }
  

  //! applies operator Dh, Prod_Uh = beta Prod_Uh + alpha Dh Uh
  void ContinuousUnsteadyMassMatrix_Base::
  ApplyOperatorDh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
		  const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    if (Dh == NULL)
      {
	cout << "Matrix Dh not computed" << endl;
	abort();
      }
    glob_chrono.Start(VirtualTimer::MASS);
    if (beta == Real_wp(0))
      {
        Real_wp zero; SetComplexZero(zero);
        Prod_Uh.Fill(zero);
      }
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);

    Dh->MltMass(alpha, Uh, Prod_Uh);
    glob_chrono.Stop(VirtualTimer::MASS);
  }

  
  
  //! applies operator Dh - dt/2 Sh, Prod_Uh = beta Prod_Uh + alpha (Dh - dt/2 Sh) Uh
  void ContinuousUnsteadyMassMatrix_Base::
  ApplyOperatorDhMinusdtSh(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			   const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    if (DhMinusdtSh == NULL)
      {
	cout << "Matrix DhMinusdtSh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    if (beta == Real_wp(0))
      Prod_Uh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);
    
    DhMinusdtSh->MltMass(alpha, Uh, Prod_Uh);
    glob_chrono.Stop(VirtualTimer::MASS);
  }
    

  //! solves by operator Dh + dt/2 Sh, Y = (Dh + dt/2)^{-1} Y
  void ContinuousUnsteadyMassMatrix_Base::
  SolveOperatorDhPlusdtSh(VectReal_wp& Y)
  {
    if (DhPlusdtSh == NULL)
      {
	cout << "Matrix DhPlusdtSh not computed" << endl;
	abort();
      }
    
    glob_chrono.Start(VirtualTimer::MASS);
    DhPlusdtSh->SolveMass(Y);
    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! solves by operator Dh + dt/2 Sh, Y = Y + alpha (Dh + dt/2)^{-1} X
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveOperatorDhPlusdtSh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    if (DhPlusdtSh == NULL)
      {
	cout << "Matrix DhPlusdtSh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    DhPlusdtSh->SolveMass(alpha, X, Y);
    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! solves by operator Dh, Y = Dh^{-1} Y
  void ContinuousUnsteadyMassMatrix_Base::
  SolveOperatorDh(VectReal_wp& Y)
  {
    if (Dh == NULL)
      {
	cout << "Matrix Dh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    Dh->SolveMass(Y);
    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! solves by Cholesky factor of Dh, Y = L^{-1} Y or Y = L^{-T} Y where Dh = L L^T
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveCholeskyDh(const SeldonTranspose& transA, VectReal_wp& Y)
  {
    abort();
  }
  

  //! solves by operator Dh, Y = Y + alpha Dh^{-1} X
  void ContinuousUnsteadyMassMatrix_Base::
  SolveOperatorDh(const Real_wp& alpha, const VectReal_wp& X, VectReal_wp& Y)
  {
    if (Dh == NULL)
      {
	cout << "Matrix Dh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    Dh->SolveMass(alpha, X, Y);
    glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! applies operator Bh, Vh = beta Vh + alpha Bh Uh    
  void ContinuousUnsteadyMassMatrix_Base
  ::ApplyOperatorBh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
		    const Real_wp& beta, VectReal_wp& Vh)
  {
    if (Bh == NULL)
      {
	cout << "Matrix Bh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    if (beta == Real_wp(0))
      Vh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Vh);
   
    if (Bh->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	Bh->MltMass(alpha, Uh, Vh);
	return;
      }
    
    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.AddElementStiffnessTensor(i, offset, alpha, Uh, Vh,
                                         Real_wp(1), Real_wp(0));
    
    glob_chrono.Stop(VirtualTimer::MASS);
  }


  //! applies operator (Bh - dt/2 ShV), Vh = beta Vh + alpha (Bh - dt/2 ShV) Uh      
  void ContinuousUnsteadyMassMatrix_Base::
  ApplyOperatorBhMinusdtSh(const Real_wp& alpha, const Real_wp& t, const VectReal_wp& Uh,
                           const Real_wp& beta, VectReal_wp& Vh)
  {
    if (BhMinusdtSh == NULL)
      {
	cout << "Matrix BhMinusdtSh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    if (beta == Real_wp(0))
      Vh.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, Vh);
   
    if (BhMinusdtSh->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	BhMinusdtSh->MltMass(alpha, Uh, Vh);
	return;
      }
    
    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.AddElementStiffnessTensor(i, offset, alpha, Uh, Vh,
                                         Real_wp(1), -0.5*var_time.GetTimeStep());

    glob_chrono.Stop(VirtualTimer::MASS);
  }
  

  //! solves operator (Bh + dt/2 ShV), Vh = (Bh + dt/2 ShV)^{-1} Vh
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveOperatorBhPlusdtSh(VectReal_wp& Vh)
  {
    if (BhPlusdtSh == NULL)
      {
	cout << "Matrix BhPlusdtSh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);

    if (BhPlusdtSh->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	BhPlusdtSh->SolveMass(Vh);
	return;
      }
    
    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.SolveElementStiffnessTensor(i, offset, Vh,
                                           0.5*var_time.GetTimeStep());

    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! solves operator Bh + dt/2 ShV, Vh = Vh + alpha (Bh + dt/2 ShV)^{-1} Uh
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveOperatorBhPlusdtSh(const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {
    if (BhPlusdtSh == NULL)
      {
	cout << "Matrix BhPlusdtSh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);

    if (BhPlusdtSh->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	BhPlusdtSh->SolveMass(alpha, Uh, Vh);
	return;
      }
    
    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.SolveAddElementStiffnessTensor(i, offset, alpha, Uh, Vh,
                                              0.5*var_time.GetTimeStep());	

    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! solves operator Bh, Vh = Bh^{-1} Vh
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveOperatorBh(VectReal_wp& Vh)
  {
    if (Bh == NULL)
      {
	cout << "Matrix Bh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);

    if (Bh->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	Bh->SolveMass(Vh);
	return;
      }

    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.SolveElementStiffnessTensor(i, offset, Vh, Real_wp(0));	
    
    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! solves operator Bh, Vh = Vh + alpha Bh^{-1} Uh
  void ContinuousUnsteadyMassMatrix_Base
  ::SolveOperatorBh(const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {
    if (Bh == NULL)
      {
	cout << "Matrix Bh not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);

    if (Bh->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	Bh->SolveMass(alpha, Uh, Vh);
	return;
      }

    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.SolveAddElementStiffnessTensor(i, offset, alpha, Uh, Vh,
                                              Real_wp(0));
    
    glob_chrono.Stop(VirtualTimer::MASS);
  }
  
  
  //! applies operator ShV, Prod_Uh = beta Prod_Uh + alpha ShV Uh    
  void ContinuousUnsteadyMassMatrix_Base
  ::ApplyOperatorShVectorial(const Real_wp& alpha, const Real_wp& tn, const VectReal_wp& Uh,
			     const Real_wp& beta, VectReal_wp& Prod_Uh)
  {
    if (ShVec == NULL)
      {
	cout << "Matrix ShVec not computed" << endl;
	abort();
      }

    glob_chrono.Start(VirtualTimer::MASS);
    if (beta == Real_wp(0))
      Prod_Uh.Fill(0);
    else if (beta != Real_wp(1))
      Mlt(beta, Prod_Uh);

    if (ShVec->GetMatrixType() != FemMassMatrix::MATRIX_FREE)
      {
	ShVec->MltMass(alpha, Uh, Prod_Uh);
	return;
      }

    int offset = 0;
    int Nelt = var_computation.GetNbElt();
    
    for (int i = 0; i < Nelt; i++)
      var_time.AddElementStiffnessTensor(i, offset, alpha, Uh, Prod_Uh,
                                         Real_wp(0), Real_wp(1));
    
    glob_chrono.Stop(VirtualTimer::MASS);
  }


  /********************************
   * ContinuousUnsteadyMassMatrix *
   ********************************/


  template<class Dimension>
  ContinuousMassMatrix_Base* ContinuousUnsteadyMassMatrix<Dimension>
  ::GetNewMassMatrix(int type_mat)
  {
    return new ContinuousMassMatrix<Dimension>(var_problem, type_mat,
                                               this->NonLumpedElement, this->NonLumpedDof,
                                               this->InverseDof, this->LumpedDof);
  }
  
} // end namespace

#define MONTJOIE_FILE_CONTINUOUS_UNSTEADY_MASS_MATRIX_CXX
#endif

