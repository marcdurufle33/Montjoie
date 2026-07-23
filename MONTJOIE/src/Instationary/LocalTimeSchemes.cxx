#ifndef MONTJOIE_FILE_LOCAL_TIME_SCHEMES_CXX

namespace Montjoie
{

  /************************
   * LocalTimeScheme_Base *
   ************************/
  

  //! default constructor
  LocalTimeScheme_Base::LocalTimeScheme_Base()
  {
    local_implicit = false;
    type_splitting = AUTO_SPLITTING;
    dt_splitting = Real_wp(0);
    p = 2;    
  }


  //! modifies how the domain is split into a coarse and fine region
  void LocalTimeScheme_Base::SetInputSplitting(const Vector<string>& param, int& nb)
  {
    if (param.GetM() <= nb+1)
      {
	cout << "Provide more parameters for LOCAL_IMPLICIT_RK" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }    

    if (param(nb) == "AUTO")
      {
	this->type_splitting = LocalTimeScheme_Base::AUTO_SPLITTING;
	this->dt_splitting = to_num<Real_wp>(param(nb+1));
      }
    else if (param(nb) == "FILE")
      {
	this->type_splitting = LocalTimeScheme_Base::FILE_SPLITTING;
	this->file_with_number_element = param(nb+1);
      }

    nb += 2;
  }


  /***********************
   * LocalOdeScheme_Base *
   ***********************/
  

  //! default constructor
  LocalOdeScheme_Base::LocalOdeScheme_Base(VectReal_wp& x, VectReal_wp& y, VectReal_wp& L, VectReal_wp& F)
    : x_glob(x), y_glob(y), L_glob(L), F_glob(F)
  {
    nb_dof_uv_local = 0; nb_dof_L_local = 0;
    nb_dof_uv_proche = 0;
    nb_dof_uv_extended = 0; nb_dof_L_extended = 0;
    nb_dof_uv_extended_nonloc = 0;
    matrix_stored = false;
    num_level_vol = -1; num_level_ch = -1;
  }


  //! internal method that creates a list of dof numbers from IndexDof
  void LocalOdeScheme_Base::CreateListIndex(int& nb_dof, IVect& liste_dof, IVect& IndexDof)
  {
    liste_dof.Reallocate(nb_dof);
    nb_dof = 0;
    for (int i = 0; i < IndexDof.GetM(); i++)
      if (IndexDof(i) == 0)
	{
	  liste_dof(nb_dof) = i;
	  IndexDof(i) = nb_dof;
	  nb_dof++;
	}
  }


  //! finds common dofs between dof_far and dof_loc
  int LocalOdeScheme_Base::FindCommonDofs(const Vector<int>& dof_far, const Vector<int>& dof_loc,
					  Vector<int>& dof_common_fine, Vector<int>& dof_common_coarse)
  {
    int Nfar_ext = dof_far.GetM();
    int Nimp_ext = dof_loc.GetM();
    
    Vector<int> Index(x_glob.GetM());
    Index.Fill(-1);
    
    for (int i = 0; i < Nfar_ext; i++)
      Index(dof_far(i)) = i;
    
    int nb_dof_common = 0;
    for (int i = 0; i < Nimp_ext; i++)
      if (Index(dof_loc(i)) >= 0)
	nb_dof_common++;
    
    dof_common_fine.Reallocate(nb_dof_common);
    dof_common_coarse.Reallocate(nb_dof_common);
    nb_dof_common = 0;
    for (int i = 0; i < Nimp_ext; i++)
      if (Index(dof_loc(i)) >= 0)
	{
	  dof_common_coarse(nb_dof_common) = Index(dof_loc(i));
	  dof_common_fine(nb_dof_common) = i;
	  nb_dof_common++;
	}

    return nb_dof_common;
  }
  

  //! creates arrays needed to compute y = A P x
  void LocalOdeScheme_Base::Init(const Vector<int>& num_element, VirtualOdeSystem<Real_wp>& sys,
				 Vector<IVect>& lvl_vol, Vector<IVect>& lvl_surface, bool compute_direct_info)
  {
    // cast of the object sys
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
    DistributedProblem_Base& var = var_time.GetDistributedProblem();
    VarProblem_Base& var_p = dynamic_cast<VarProblem_Base&>(var);

    Vector<int> IndexLocalUV, IndexLocalL, IndexProcheUV;
    Vector<int> IndexExtendedUV, IndexExtendedL;
    
    // creating list of dofs on the local zone (for u/v and lambda)
    int Nvol = var.GetNbMainUnknownDof();
    int nb_dof_L = var_p.GetNbMeshDof()*var_p.nb_unknowns_hdg;
    int nb_dof_uv = Nvol*var_p.nb_unknowns;

    // counting the dofs
    num_elt_local = num_element;
    IndexLocalL.Reallocate(nb_dof_L); IndexLocalUV.Reallocate(nb_dof_uv);
    IndexLocalL.Fill(-1); IndexLocalUV.Fill(-1);
    nb_dof_L_local = 0; nb_dof_uv_local = 0;
    for (int i2 = 0; i2 < num_element.GetM(); i2++)
      {
	int i = num_element(i2);
	IVect Nodle = var.GetNodleElement(i);
	for (int j = 0; j < Nodle.GetM(); j++)
	  if (IndexLocalL(Nodle(j)) == -1)
	    IndexLocalL(Nodle(j)) = 0;
	
	int v1 = var.GetOffsetDofV(i);
	int v2 = var.GetOffsetDofV(i+1);
	int nb_dof_v = v2 - v1;
	Vector<int> dof_v(var_p.nb_unknowns*nb_dof_v);
	for (int k = 0; k < var_p.nb_unknowns; k++)
	  for (int i = v1; i < v2; i++)
	    if (IndexLocalUV(k*Nvol + i) == -1)
	      {
		IndexLocalUV(k*Nvol + i) = 0;
		nb_dof_uv_local++;
	      }
      }

    // assembling IndexLocalL to take into account neighbor processors
#ifdef SELDON_WITH_MPI
    var.ReduceDistributedVector(IndexLocalL, MPI_MAX);
#endif

    for (int i = 0; i < IndexLocalL.GetM(); i++)
      if (IndexLocalL(i) == 0)
	nb_dof_L_local++;

    // listing the dofs
    CreateListIndex(nb_dof_L_local, dof_local_L, IndexLocalL);
    CreateListIndex(nb_dof_uv_local, dof_local_uv, IndexLocalUV);
    
    // finite element matrices involved
    DistributedMatrix<Real_wp, General, ArrayRowSparse>& A22 = var_time.Glob_mat_RhV;
    DistributedMatrix<Real_wp, General, ArrayRowSparse>& A21 = var_time.Glob_mat_Ch_Lambda;
    DistributedMatrix<Real_wp, General, ArrayRowSparse>& A12 = var_time.Glob_mat_RhS;

    DiscontinuousUnsteadyMassMatrix_Base& mat_mass
      = dynamic_cast<DiscontinuousUnsteadyMassMatrix_Base& >(*var_time.Glob_mat_Dh);
    
    DiscontinuousMassMatrix_Base& Ch = mat_mass.GetOperatorCh();
    if (Ch.GetMatrixType() != FemMassMatrix::BLOCK_ELT)
      {
	cout << "Storage of Ch not handled" << endl;
	abort();
      }

    matrix_stored = false;
    if (A22.GetM() > 0)
      matrix_stored = true;

    // levels used for the matrix vector product
    MatrixVectorProductLevel& level_Ch = Ch.GetLevelTime();
    MatrixVectorProductLevel& level_vol = var_time.GetTimeLevelDistribution();
    
    // calcul des ddls proches
    VectReal_wp prod_aleatoire(nb_dof_L), vec_aleatoire(nb_dof_uv);
    
    prod_aleatoire.Zero(); vec_aleatoire.Zero();
    for (int i = 0; i < nb_dof_L_local; i++)
      {
	int ig = dof_local_L(i);
	prod_aleatoire(ig) = Real_wp(rand()) / RAND_MAX;
      }

    for (int i = 0; i < nb_dof_uv_local; i++)
      {
	int ig = dof_local_uv(i);
	vec_aleatoire(ig) = Real_wp(rand()) / RAND_MAX;
      }
    
    y_glob.Reallocate(nb_dof_uv); y_glob.Zero();
    F_glob.Reallocate(nb_dof_uv + nb_dof_L); F_glob.Zero();

    x_glob.Nullify();
    x_glob.SetData(nb_dof_uv, &F_glob(nb_dof_L));

    L_glob.Nullify();
    L_glob.SetData(nb_dof_L, F_glob.GetData());
    
    if (matrix_stored)      
      MltAdd(Real_wp(1), A21, prod_aleatoire, Real_wp(1), vec_aleatoire); 
    else
      {
	var_time.GetExtrapolationLambdaHDG(-1, y_glob, L_glob);
	var.AddDomains(prod_aleatoire);
	var_time.MltAddStiffnessHDG(Real_wp(1), -1, -1, y_glob, prod_aleatoire, Real_wp(1), vec_aleatoire);
      }
    
    var_time.SolveMassMatrix(vec_aleatoire);
    
    nb_dof_uv_proche = 0;
    IndexProcheUV.Reallocate(nb_dof_uv);
    IndexProcheUV.Fill(-1);
    for (int i = 0; i < vec_aleatoire.GetM(); i++)
      if (vec_aleatoire(i) != Real_wp(0))
	{
	  IndexProcheUV(i) = 0;
	  nb_dof_uv_proche++;
	}
    
    CreateListIndex(nb_dof_uv_proche, dof_proche_uv, IndexProcheUV);
    
    // the elements that belong to the close dofs
    int nb_elt = level_vol.GetNbElt();
    IVect IndexElement(nb_elt);
    IndexElement.Fill(-1);
    int nb_elt_proche = 0;
    for (int i = 0; i < nb_elt; i++)
      for (int j = var.GetOffsetDofV(i); j < var.GetOffsetDofV(i+1); j++)
	if (IndexProcheUV(j) >= 0)
	  if (IndexElement(i) == -1)
	    {
	      IndexElement(i) = 0;
	      nb_elt_proche++;
	    }

    CreateListIndex(nb_elt_proche, num_elt_extended, IndexElement);
    //DISP(nb_elt_proche); DISP(num_elt_extended); DISP(dof_proche_uv);
    
    // we create an extended zone (with elements contained in num_elt_extended)
    IndexExtendedL.Reallocate(nb_dof_L); IndexExtendedUV.Reallocate(nb_dof_uv);
    IndexExtendedL.Fill(-1); IndexExtendedUV.Fill(-1);
    nb_dof_L_extended = 0; nb_dof_uv_extended = 0;
    for (int i2 = 0; i2 < num_elt_extended.GetM(); i2++)
      {
	int i = num_elt_extended(i2);
	IVect Nodle = var.GetNodleElement(i);
	for (int j = 0; j < Nodle.GetM(); j++)
	  if (IndexExtendedL(Nodle(j)) == -1)
	    IndexExtendedL(Nodle(j)) = 0;
	
	int v1 = var.GetOffsetDofV(i);
	int v2 = var.GetOffsetDofV(i+1);
	int nb_dof_v = v2 - v1;
	Vector<int> dof_v(var_p.nb_unknowns*nb_dof_v);
	for (int k = 0; k < var_p.nb_unknowns; k++)
	  for (int i = v1; i < v2; i++)
	    if (IndexExtendedUV(k*Nvol + i) == -1)
	      {
		IndexExtendedUV(k*Nvol + i) = 0;
		nb_dof_uv_extended++;
	      }
      }

#ifdef SELDON_WITH_MPI
    var.ReduceDistributedVector(IndexExtendedL, MPI_MAX);
#endif

    for (int i = 0; i < IndexExtendedL.GetM(); i++)
      if (IndexExtendedL(i) == 0)
	nb_dof_L_extended++;

    // listing the dofs
    CreateListIndex(nb_dof_L_extended, dof_extended_L, IndexExtendedL);
    CreateListIndex(nb_dof_uv_extended, dof_extended_uv, IndexExtendedUV);
    
    // matrices are extracted (if stored) for the close region
    int nb_proc = var.GetNbProcPerMode();
#ifdef SELDON_WITH_MPI
    int rank_proc = var.GetRankProcMode();
#endif
    if (matrix_stored)
      {
	ExtractSubMatrix(A12, dof_local_L, IndexLocalL, dof_local_uv, IndexLocalUV, A12_sub);
	ExtractSubMatrix(A21, dof_proche_uv, IndexProcheUV, dof_local_L, IndexLocalL, A21_sub);
	ExtractSubMatrix(A22, dof_proche_uv, IndexProcheUV, dof_local_uv, IndexLocalUV, A22_sub);

	if (nb_proc > 1)
	  {
	    cout << "Functionality not tested" << endl;
	    abort();
	  }
      }

    // levels for volume unknowns
    
    // level n : local zone
    // level n+1 : extended zone
    // level n+2 : only new elements
    num_level_vol = lvl_vol.GetM();
    lvl_vol.PushBack(num_elt_local);
    lvl_vol.PushBack(num_elt_extended);

    IndexElement.Fill(-1);
    for (int i = 0; i < num_elt_local.GetM(); i++)
      IndexElement(num_elt_local(i)) = 0;

    IVect num_elt_new;
    num_elt_new.Reallocate(num_elt_extended.GetM() - num_elt_local.GetM());
    nb_elt_proche = 0;
    for (int i = 0; i < num_elt_extended.GetM(); i++)
      if (IndexElement(num_elt_extended(i)) == -1)
	{
	  num_elt_new (nb_elt_proche) = num_elt_extended(i);
	  nb_elt_proche++;
	}

    lvl_vol.PushBack(num_elt_new);
	
    level_vol.SetLevelArray(lvl_vol);

    // setting level_vol in mat_mass
    DiscontinuousMassMatrix_Base& Dh = mat_mass.GetOperatorDh();
    DiscontinuousMassMatrix_Base& Bh = mat_mass.GetOperatorBh();

    Bh.GetLevelTime() = level_vol;
    Dh.GetLevelTime() = level_vol;
    
    // levels for unknown lambda
    IVect num_surface_proche;
    int nb_boundaries = level_Ch.GetNbElt();
    IVect IndexBoundary(nb_boundaries);
    IndexBoundary.Fill(-1);
    const IVect& block_size_L = Ch.GetOffsetSizeBlocks();
    int nblock = 0, nb_element_Ch_proche = 0;
    for (int i = 0; i < nb_dof_L_local; i++)
      {
	int j = dof_local_L(i);
	while ((nblock < nb_boundaries) && (j >= block_size_L(nblock+1)))
	  nblock++;
	
	if (nblock < nb_boundaries)
	  {
	    if (IndexBoundary(nblock) == -1)
	      {
		IndexBoundary(nblock) = nb_element_Ch_proche;
		nb_element_Ch_proche++;
	      }
	  }
      }
    
    num_surface_proche.Reallocate(nb_element_Ch_proche);
    nb_element_Ch_proche = 0;
    for (int i = 0; i < IndexBoundary.GetM(); i++)
      if (IndexBoundary(i) >= 0)
	{
	  num_surface_proche(nb_element_Ch_proche) = i;
	  nb_element_Ch_proche++;
	}
    
    num_level_ch = lvl_surface.GetM();
    lvl_surface.PushBack(num_surface_proche);
    
    level_Ch.SetLevelArray(lvl_surface);
    
    // dofs inside the extended region but non in the local region
    Vector<int> IndexDof(nb_dof_uv);
    IndexDof.Fill(-1);
    nb_dof_uv_extended_nonloc = 0;
    for (int i = 0; i < nb_dof_uv; i++)
      if ((IndexLocalUV(i) < 0) && (IndexExtendedUV(i) >= 0))
	{
	  IndexDof(i) = 0;
	  nb_dof_uv_extended_nonloc++;
	}
    
    CreateListIndex(nb_dof_uv_extended_nonloc, dof_extended_uv_nonloc, IndexDof);

    // creates arrays ProcNumberLocalL and RowNumberLocalL needed for AssembleLocal
#ifdef SELDON_WITH_MPI
    const Vector<IVect>& RowNumberGlobalL = var.GetOriginalMatchingDofNeighbor();
    const IVect& ProcNumberGlobalL = var.GetProcMatchingNeighbor();

    Vector<bool> is_local(nb_dof_L);
    is_local.Fill(false);
    for (int i = 0; i < nb_dof_L; i++)
      if (IndexLocalL(i) >= 0)
	is_local(i) = true;

    ExtractDistributedSharedNumbers(ProcNumberGlobalL, RowNumberGlobalL,
				    is_local, nb_dof_L_local,
				    ProcNumberLocalL, RowNumberLocalL);

    // creates arrays needed in mat_direct_info
    // OverlapRowNumbers, OverlapProcNumbers, GlobalRowNumbers, ProcSharingRows, SharingRowNumbers
    // nodl_scalar, nb_unknowns_scal, nglob, comm
    if (!compute_direct_info)
      return;
    
    const IVect& OverlapRowGlobalL = var.GetOverlappedDofNumber();
    const IVect& OverlapProcGlobalL = var.GetOverlappedProcNumber();
    Vector<int> IndexOverlap(nb_dof_L);
    IndexOverlap.Fill(-1);
    for (int i = 0; i < OverlapRowGlobalL.GetM(); i++)
      IndexOverlap(OverlapRowGlobalL(i)) = i;
    
    int nb_overlap = 0;
    for (int i = 0; i < nb_dof_L; i++)
      if ((IndexOverlap(i) >= 0) && (IndexLocalL(i) >= 0))
	nb_overlap++;

    IVect& OverlapRowLocalL = mat_direct_info.OverlapRowNumbers;
    IVect& OverlapProcLocalL = mat_direct_info.OverlapProcNumbers;
    OverlapRowLocalL.Reallocate(nb_overlap);
    OverlapProcLocalL.Reallocate(nb_overlap);
    nb_overlap = 0;
    for (int i = 0; i < nb_dof_L; i++)
      if ((IndexOverlap(i) >= 0) && (IndexLocalL(i) >= 0))
	{
	  OverlapRowLocalL(nb_overlap) = IndexLocalL(i);
	  OverlapProcLocalL(nb_overlap) = OverlapProcGlobalL(IndexOverlap(i));
	  nb_overlap++;
	}
    
    mat_direct_info.ProcSharingRows = ProcNumberLocalL;
    mat_direct_info.SharingRowNumbers.Reallocate(RowNumberLocalL.GetM());
    for (int i = 0; i < RowNumberLocalL.GetM(); i++)
      {
	int nb_dof = RowNumberLocalL(i).GetM();
	mat_direct_info.SharingRowNumbers(i).Reallocate(nb_dof);
	for (int j = 0; j < nb_dof; j++)
	  mat_direct_info.SharingRowNumbers(i)(j) = IndexLocalL(RowNumberLocalL(i)(j));

      }

    int nb_dof_L_non_overlap = nb_dof_L_local - nb_overlap;
    Vector<int> nb_rows_per_proc;
    if (rank_proc == 0)
      nb_rows_per_proc.Reallocate(nb_proc);
    
    MPI_Gather(&nb_dof_L_non_overlap, 1, MPI_INTEGER,
	       nb_rows_per_proc.GetData(), 1, MPI_INTEGER, 0, var.comm_group_mode);

    Vector<int> offset_row;
    int nglob;
    if (rank_proc == 0)
      {
	offset_row.Reallocate(nb_proc+1);
	offset_row(0) = 0;
	for (int i = 0; i < nb_proc; i++)
	  offset_row(i+1) = offset_row(i) + nb_rows_per_proc(i);

	nglob = offset_row(nb_proc);
      }

    int offset_row_loc = 0;
    MPI_Scatter(offset_row.GetData(), 1, MPI_INTEGER, &offset_row_loc,
		1, MPI_INTEGER, 0, var.comm_group_mode);

    MPI_Bcast(&nglob, 1, MPI_INTEGER, 0, var.comm_group_mode);

    mat_direct_info.nglob = nglob;
    IVect& global_row = mat_direct_info.GlobalRowNumbers;
    global_row.Reallocate(nb_dof_L_local);
    global_row.Fill(-1);
    int nb = 0;
    for (int i = 0; i < nb_dof_L_local; i++)
      if (IndexOverlap(dof_local_L(i)) == -1)
	global_row(i) = offset_row_loc + nb++;
    
    var.ReduceDistributedVector(global_row, MPI_MAX, mat_direct_info.ProcSharingRows,
				mat_direct_info.SharingRowNumbers, nb_dof_L_local, 1);

    mat_direct_info.nodl_scalar = nb_dof_L_local;
    mat_direct_info.nb_unknowns_scal = 1;
    mat_direct_info.nloc = nb_dof_L_local;
    mat_direct_info.comm = var.comm_group_mode;
#endif
  }


  //! adds contributions for Lambda (only local part)
  void LocalOdeScheme_Base::AssembleLocal(VirtualOdeSystem<Real_wp>& sys,
					  VectReal_wp& L)
  {
#ifdef SELDON_WITH_MPI
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
    DistributedProblem_Base& var = var_time.GetDistributedProblem();
    VarProblem_Base& var_p = dynamic_cast<VarProblem_Base&>(var);

    int nb_dof_L = var_p.GetNbMeshDof()*var_p.nb_unknowns_hdg;
    int nb_u = 1;
    
    var.ReduceDistributedVector(L, MPI_SUM, ProcNumberLocalL,
				RowNumberLocalL, nb_dof_L, nb_u);
#endif
  }


  //! adds contributions for Lambda (only local part)
  void LocalOdeScheme_Base::AssembleLocal(VirtualOdeSystem<Real_wp>& sys,
					  VectComplex_wp& L)
  {
#ifdef SELDON_WITH_MPI
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
    DistributedProblem_Base& var = var_time.GetDistributedProblem();
    VarProblem_Base& var_p = dynamic_cast<VarProblem_Base&>(var);

    int nb_dof_L = var_p.GetNbMeshDof()*var_p.nb_unknowns_hdg;
    int nb_u = 1;
    
    var.ReduceDistributedVector(L, MPI_SUM, ProcNumberLocalL,
				RowNumberLocalL, nb_dof_L, nb_u);
#endif
  }

  
  //! Computes Y = A P X
  void LocalOdeScheme_Base::ApplyOperatorAP(const VectReal_wp& X, VectReal_wp& Y,
					    VirtualOdeSystem<Real_wp>& sys, bool invert_mass)
  {
    // cast to retrieve var_time, Ch and level_Ch
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
    //DistributedProblem_Base& var = var_time.GetDistributedProblem();
 
    DiscontinuousUnsteadyMassMatrix_Base& mat_mass
      = dynamic_cast<DiscontinuousUnsteadyMassMatrix_Base& >(*var_time.Glob_mat_Dh);
    
    DiscontinuousMassMatrix_Base& Ch = mat_mass.GetOperatorCh();
    MatrixVectorProductLevel& level_Ch = Ch.GetLevelTime();
    
    VectReal_wp ProdLambda, ProdY;
    
    // multiplication by -A12
    if (matrix_stored)
      {
	ProdLambda.Reallocate(nb_dof_L_local);
	ProdLambda.Zero();
    
	MltAdd(Real_wp(-1), A12_sub, X, Real_wp(0), ProdLambda);

	for (int i = 0; i < nb_dof_L_local; i++)
	  L_glob(dof_local_L(i)) = ProdLambda(i);
      }
    else
      {
	for (int i = 0; i < nb_dof_uv_extended_nonloc; i++)
	  x_glob(dof_extended_uv_nonloc(i)) = 0;
	
	for (int i = 0; i < nb_dof_uv_local; i++)
	  x_glob(dof_local_uv(i)) = X(i);

	for (int i = 0; i < nb_dof_L_extended; i++)
	  L_glob(dof_extended_L(i)) = 0;

	var_time.GetExtrapolationLambdaHDG(num_level_vol, x_glob, L_glob);
      }

    //DISP(L_glob);
    
    // then we invert by A11 to recover lambda
    level_Ch.SetLevel(num_level_ch);

    // assembling locally L_glob (AddDomains should be slower)
    // var.AddDomains(L_glob);
    AssembleLocal(sys, L_glob);
    Ch.SolveMass(0, L_glob); 

    //DISP(L_glob);
    
    // then we compute - A21 Lambda - A22 X to obtain Y
    if (matrix_stored)
      {
	ProdY.Reallocate(nb_dof_uv_proche);
	ProdY.Zero();
	
	for (int i = 0; i < nb_dof_L_local; i++)
	  ProdLambda(i) = L_glob(dof_local_L(i));

	MltAdd(-Real_wp(1), A21_sub, ProdLambda, Real_wp(0), ProdY);
	MltAdd(-Real_wp(1), A22_sub, X, Real_wp(1), ProdY);

	// we invert by the mass matrix
	if (invert_mass)
	  this->SolveOperatorM(sys, ProdY, dof_proche_uv, false);

	Y = ProdY;
      }
    else
      {
	for (int i = 0; i < nb_dof_uv_extended; i++)
	  y_glob(dof_extended_uv(i)) = 0;
		
	var_time.MltAddStiffnessHDG(Real_wp(1), num_level_vol+1, num_level_vol,
				    x_glob, L_glob, Real_wp(1), y_glob);

	// we invert by the mass matrix
	if (invert_mass)
	  SolveOperatorM(sys, y_glob, dof_proche_uv, true);
	
	// then we extract the result in Y
	for (int i = 0; i < nb_dof_uv_proche; i++)
	  Y(i) = y_glob(dof_proche_uv(i));
      }

    // we go back to level -1 for Ch (for AddPrimitiveTimeSource)
    level_Ch.SetLevel(-1);
  }

  
  //! solves by mass matrix
  template<class T>
  void LocalOdeScheme_Base::SolveOperatorM(VirtualOdeSystem<Real_wp>& sys, Vector<T>& Y,
					   const Vector<int>& liste_dof, bool global_vec)
  {
    // cast to retrieve var_time, Ch and level_Ch
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
    DistributedProblem_Base& var = var_time.GetDistributedProblem();
    VarProblem_Base& var_p = dynamic_cast<VarProblem_Base&>(var);
 
    DiscontinuousUnsteadyMassMatrix_Base& mat_mass
      = dynamic_cast<DiscontinuousUnsteadyMassMatrix_Base& >(*var_time.Glob_mat_Dh);
    
    DiscontinuousMassMatrix_Base& Dh = mat_mass.GetOperatorDh();
    DiscontinuousMassMatrix_Base& Bh = mat_mass.GetOperatorBh();

    if ((Dh.GetMatrixType() == FemMassMatrix::DIAGONAL) && (Bh.GetMatrixType() == FemMassMatrix::DIAGONAL))
      {
	int Nvol = var.GetNbMainUnknownDof();
	int offset_u = Nvol*var_p.nb_unknowns_scal;
	const VectReal_wp& invDh = Dh.GetInverseDiagonal();
	const VectReal_wp& invBh = Bh.GetInverseDiagonal();

	if (global_vec)
	  for (int i = 0; i < liste_dof.GetM(); i++)
	    {
	      int ig = liste_dof(i);
	      if (ig < offset_u)
		Y(ig) *= invDh(ig);
	      else
		Y(ig) *= invBh(ig-offset_u);
	    }
	else
	  for (int i = 0; i < liste_dof.GetM(); i++)
	    {
	      int ig = liste_dof(i);
	      if (ig < offset_u)
		Y(i) *= invDh(ig);
	      else
		Y(i) *= invBh(ig-offset_u);
	    }	  
      }
    else if ((Dh.GetMatrixType() == FemMassMatrix::BLOCK_ELT) && (Bh.GetMatrixType() == FemMassMatrix::BLOCK_ELT))
      {
	Dh.GetLevelTime().SetLevel(num_level_vol+1);
	Bh.GetLevelTime().SetLevel(num_level_vol+1);
	
	if (!global_vec)
	  {
	    // for local vectors, we put values on y_glob
	    for (int i = 0; i < nb_dof_uv_extended_nonloc; i++)
	      {
		x_glob(dof_extended_uv_nonloc(i)) = 0;
		y_glob(dof_extended_uv_nonloc(i)) = 0;
	      }

	    for (int i = 0; i < liste_dof.GetM(); i++)
	      y_glob(liste_dof(i)) = realpart(Y(i));

	    if (sizeof(T) != sizeof(Real_wp))
	      {
		cout << "Impossible" << endl;
		abort();
	      }
	    
	    // then we solve y_glob (for extended zone)
	    mat_mass.SolveMassMatrix(y_glob);

	    // then we take back values on the local vector
	    for (int i = 0; i < liste_dof.GetM(); i++)
	      Y(i) = y_glob(liste_dof(i));
	  }
	else
	  {
	    for (int i = 0; i < nb_dof_uv_extended_nonloc; i++)
	      y_glob(dof_extended_uv_nonloc(i)) = 0;
	    
	    if (sizeof(T) == sizeof(Real_wp))
	      {
		for (int i = 0; i < liste_dof.GetM(); i++)
		  y_glob(liste_dof(i)) = realpart(Y(liste_dof(i)));
		
		mat_mass.SolveMassMatrix(y_glob);

		for (int i = 0; i < liste_dof.GetM(); i++)
		  Y(liste_dof(i)) = y_glob(liste_dof(i));
	      }
	    else
	      {
		for (int i = 0; i < nb_dof_uv_extended_nonloc; i++)
		  x_glob(dof_extended_uv_nonloc(i)) = 0;
		
		for (int i = 0; i < liste_dof.GetM(); i++)
		  {
		    x_glob(liste_dof(i)) = imagpart(Y(liste_dof(i)));
		    y_glob(liste_dof(i)) = realpart(Y(liste_dof(i)));
		  }
		
		mat_mass.SolveMassMatrix(x_glob);
		mat_mass.SolveMassMatrix(y_glob);

		for (int i = 0; i < liste_dof.GetM(); i++)
		  to_complex(Complex_wp(x_glob(liste_dof(i)), y_glob(liste_dof(i))), Y(liste_dof(i)));
	      }
	  }
      }
    else
      {	
	DISP(Dh.GetMatrixType()); DISP(Bh.GetMatrixType());
	cout << "Storage of Dh and Bh not implemented" << endl;
	abort();
      }

    //int test_input; cout << "waiting" << endl; cin >> test_input;
  }


  //! adds mass matrix
  void LocalOdeScheme_Base::AddOperatorM(VirtualOdeSystem<Real_wp>& sys, const Real_wp& alpha,
					 const VectReal_wp& X, VectReal_wp& Y,
					 const Vector<int>& liste_dof)
  {
    // cast to retrieve var_time, Ch and level_Ch
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
    DistributedProblem_Base& var = var_time.GetDistributedProblem();
    VarProblem_Base& var_p = dynamic_cast<VarProblem_Base&>(var);    
     
    DiscontinuousUnsteadyMassMatrix_Base& mat_mass
      = dynamic_cast<DiscontinuousUnsteadyMassMatrix_Base& >(*var_time.Glob_mat_Dh);
    
    DiscontinuousMassMatrix_Base& Dh = mat_mass.GetOperatorDh();
    DiscontinuousMassMatrix_Base& Bh = mat_mass.GetOperatorBh();

    if ((Dh.GetMatrixType() == FemMassMatrix::DIAGONAL) && (Bh.GetMatrixType() == FemMassMatrix::DIAGONAL))
      {
	int Nvol = var.GetNbMainUnknownDof();
	int offset_u = Nvol*var_p.nb_unknowns_scal;
	const VectReal_wp& diagDh = Dh.GetDiagonal();
	const VectReal_wp& diagBh = Bh.GetDiagonal();

	for (int i = 0; i < liste_dof.GetM(); i++)
	  {
	    int ig = liste_dof(i);
	    if (ig < offset_u)
	      Y(i) += alpha*diagDh(ig)*X(i);
	    else
	      Y(i) += alpha*diagBh(ig-offset_u)*X(i);
	  }	  
      }
    else if ((Dh.GetMatrixType() == FemMassMatrix::BLOCK_ELT) && (Bh.GetMatrixType() == FemMassMatrix::BLOCK_ELT))
      {
	Dh.GetLevelTime().SetLevel(num_level_vol+1);
	Bh.GetLevelTime().SetLevel(num_level_vol+1);
	
	// for local vectors, we put values on y_glob
	for (int i = 0; i < nb_dof_uv_extended_nonloc; i++)
	  {
	    x_glob(dof_extended_uv_nonloc(i)) = 0;
	    y_glob(dof_extended_uv_nonloc(i)) = 0;
	  }
	
	for (int i = 0; i < liste_dof.GetM(); i++)
	  {
	    x_glob(liste_dof(i)) = X(i);
	    y_glob(liste_dof(i)) = Y(i);
	  }

	// then we solve y_glob (for extended zone)
	mat_mass.ApplyMassMatrix(alpha, Real_wp(0), x_glob, Real_wp(1), y_glob);
	
	// then we take back values on the local vector
	for (int i = 0; i < liste_dof.GetM(); i++)
	  Y(i) = y_glob(liste_dof(i));
      }
    else
      {
	DISP(Dh.GetMatrixType()); DISP(Bh.GetMatrixType());
	cout << "Storage of Dh and Bh not implemented" << endl;
	abort();
      }
    
    //int test_input; cout << "waiting" << endl; cin >> test_input;
  }
  

  /******************************
   * VarComputationProblemLocal *
   ******************************/
  
  
  //! constructor
  VarComputationProblemLocal
  ::VarComputationProblemLocal(VarComputationProblem& var_global_,
			       LocalOdeScheme_Base& var_local_,
			       MatrixVectorProductLevel& level_vol_,
			       Vector<int>& IndexU, Vector<int>& IndexL)
    : var_global(var_global_), var_local(var_local_), level_vol(level_vol_),
      IndexProche(IndexU), IndexLambda(IndexL), level_inside(level_vol)
  {
    level_vol.SetLevel(var_local.GetVolumeLevel()+1);
    level_inside.SetLevel(var_local.GetVolumeLevel());
  }


  //! returns the number of elements needed to compute the local matrix
  int VarComputationProblemLocal::GetNbElt() const
  {
    return level_vol.GetNbElt();
  }
  

  //! returns the size of the local matrix
  int VarComputationProblemLocal::GetNbRows() const
  {
    return var_local.GetNbDofProche() + var_local.GetNbLocalDofLambda();
  }
  

  //! returns the print level
  int VarComputationProblemLocal::GetPrintLevel() const
  {
    return -1;
  }
    

  //! computes the elementary matrix
  void VarComputationProblemLocal
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem, CondensationBlockSolver_Base<Real_wp>& solver,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    // nelem : local element number, iquad : global element number
    solver.SetNbCondensedElt(level_inside.GetNbElt());
    int nelem = level_inside.GetLocalElementNumber();
    solver.SetElementNumber(nelem, i);
    int iquad = level_vol.GetElementNumber(i);
    bool treat_inside = level_inside.TreatElement(iquad);
    solver.SetTreatmentStiffnessInside(treat_inside);

    // we compute the elementary matrix of element iquad
    var_global.ComputeElementaryMatrix(iquad, num_dof, mat_elem, solver, nat_mat);

    // then we change the dof numbers (to local numbers)
    int nL = var_local.GetNbLocalDofLambda();
    int nb_dof_L = IndexLambda.GetM();
    
    for (int i = 0; i < num_dof.GetM(); i++)
      {
	if (num_dof(i) >= nb_dof_L)
	  {
	    int i0 = IndexProche(num_dof(i)-nb_dof_L);
	    if (i0 >= 0)
	      num_dof(i) = nL + i0;
	    else
	      num_dof(i) = -1;
	  }
	else
	  num_dof(i) = IndexLambda(num_dof(i));
      }
  }
    

  //! computes the elementary matrix
  void VarComputationProblemLocal
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem, CondensationBlockSolver_Base<Complex_wp>& solver,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    solver.SetNbCondensedElt(level_inside.GetNbElt());
    int nelem = level_inside.GetLocalElementNumber();
    solver.SetElementNumber(nelem, i);
    int iquad = level_vol.GetElementNumber(i);
    bool treat_inside = level_inside.TreatElement(iquad);
    solver.SetTreatmentStiffnessInside(treat_inside);
    
    var_global.ComputeElementaryMatrix(iquad, num_dof, mat_elem, solver, nat_mat);

    int nL = var_local.GetNbLocalDofLambda();
    int nb_dof_L = IndexLambda.GetM();
    for (int i = 0; i < num_dof.GetM(); i++)
      {
	if (num_dof(i) >= nb_dof_L)
	  {
	    int i0 = IndexProche(num_dof(i)-nb_dof_L);
	    if (i0 >= 0)
	      num_dof(i) = nL + i0;
	    else
	      num_dof(i) = -1;
	  }
	else
	  num_dof(i) = IndexLambda(num_dof(i));
      }
  }
  

  //! sets dof numbers to be condensed
  void VarComputationProblemLocal
  ::GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
			    Vector<int>& intern_node) const
  {
    nb_dof_int = 0; nb_dof_edges = nb_dof_loc;
    intern_node.Clear();
  }
  
  
  
  /*********************
   * LocalOdeSchemeHDG *
   *********************/


  //! constructor with interpolation points ci and time step dt0
  LocalOdeSchemeHDG::LocalOdeSchemeHDG(const VectReal_wp& ci, const Real_wp& dt0, VirtualOdeSystem<Real_wp>& glob_ode,
				       const IVect& dofL, VectReal_wp& x , VectReal_wp& y, VectReal_wp& L, VectReal_wp& Fg)
    : LocalOdeScheme_Base(x, y, L, Fg), global_ode(glob_ode), dof_local(dofL)
  {
    lob.AffectPoints(ci);
    dt = dt0;
    static_condensation = true;
    cond_real_solver = NULL;
  }


  //! destructor
  LocalOdeSchemeHDG::~LocalOdeSchemeHDG()
  {
    x_cplx.Nullify();
    L_cplx.Nullify();
    
    if (cond_real_solver != NULL)
      delete cond_real_solver;

    for (int i = 0; i < cond_cplx_solver.GetM(); i++)
      if (cond_cplx_solver(i) != NULL)
	delete cond_cplx_solver(i);
  }

  //! applies the function 
  void LocalOdeSchemeHDG::EvaluateFunction(const Real_wp& t, const VectReal_wp& y, VectReal_wp& f,
					   bool invert_mass, bool source)
  {
    y_imp.Reallocate(dof_local.GetM());
    for (int i = 0; i < dof_local.GetM(); i++)
      y_imp(i) = y(dof_local(i));
    
    this->ApplyOperatorAP(y_imp, f, global_ode, false);
    if (source)
      AddPrimitiveTimeSource(Real_wp(1), t, 0, f);

    if (invert_mass)
      this->SolveMassMatrix(f);
  }


  //! not used
  void LocalOdeSchemeHDG::SetDirichletCondition(const Real_wp& t, int n, VectReal_wp& Y, Real_wp alpha)
  {
  }


  //! there iw always a source because of the coarse contribution
  Real_wp LocalOdeSchemeHDG::GetFinalTimeSource() const
  {
    return 1e300;
  }


  //! sets the source (coarse contribution and local source)
  void LocalOdeSchemeHDG::SetSource(const Real_wp& t, const Vector<VectReal_wp>& wp_i, const Vector<int>& dof_wp,
				    bool presence_source_, const Vector<Vector<Real_wp, VectSparse> >& Fp_i)
  {
    tn = t;
    wp_i_ptr = &wp_i;
    dof_wp_ptr = &dof_wp;
    Fp_i_ptr = &Fp_i;
    presence_source = presence_source_;
  }
  

  //! solves by mass matrix
  void LocalOdeSchemeHDG::SolveMassMatrix(Vector<Real_wp>& X)
  {
    this->SolveOperatorM(global_ode, X, this->dof_proche_uv, false);
  }


  //! solves by mass matrix
  void LocalOdeSchemeHDG::ApplyMassMatrix(const Real_wp& alpha, const Real_wp& t,
					  const Vector<Real_wp>& x, const Real_wp& beta, Vector<Real_wp>& y)
  {
    if (beta != Real_wp(1))
      y *= beta;
    
    this->AddOperatorM(global_ode, alpha, x, y, this->dof_proche_uv);
  }
  

  //! adds alpha*F to b_src
  void LocalOdeSchemeHDG
  ::AddPrimitiveTimeSource(const Real_wp& alpha, const Real_wp& t,
			   int nb_deriv, Vector<Real_wp>& b_src)
  {
    Real_wp tau = t - tn;
    const Vector<VectReal_wp>& wp_i = *wp_i_ptr;
    const Vector<Vector<Real_wp, VectSparse> >& Fp_i = *Fp_i_ptr;
    const Vector<int>& dof_wp = *dof_wp_ptr;
    
    F = wp_i(wp_i.GetM()-1);
    for (int j = wp_i.GetM()-2; j >= 0; j--)
      F = tau*F + wp_i(j);
    
    if (presence_source)
      {
	Real_wp xi = tau / dt;
	VectReal_wp phi;
	lob.ComputeValuesPhiRef(xi, phi);
	for (int i = 0; i < phi.GetM(); i++)
          Add(alpha*phi(i), Fp_i(i), b_src);
      }
    
    for (int i = 0; i < F.GetM(); i++)
      {
	int ig = dof_wp(i);
	b_src(ig) = b_src(ig) + alpha*F(i);
      }
  }


  //! computes and factorizes a M + c K (real or complex)
  template<class T>
  void LocalOdeSchemeHDG::FactorizeGeneric(const T& a, const T& b, const T& c,
					   CondensationBlockSolver_Base<T>*& cond_solver,
					   GlobalGenericMatrix<T>& nat_mat, All_MatrixLU<T>& mat_lu)
  {
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(global_ode);
    DistributedProblem_Base& var = var_time.GetDistributedProblem();
    VarComputationProblem& var_computation = dynamic_cast<VarComputationProblem&>(var);
    MatrixVectorProductLevel& level_vol = var_time.GetTimeLevelDistribution();

    // global size of system
    int N = this->nb_dof_uv_local + this->nb_dof_L_local;

    // Global to local indexes
    Vector<int> IndexLocal(x_glob.GetM());
    IndexLocal.Fill(-1);
    for (int i = 0; i < this->nb_dof_uv_local; i++)
      IndexLocal(this->dof_local_uv(i)) = i;
    
    Vector<int> IndexLambda(L_glob.GetM());
    IndexLambda.Fill(-1);
    for (int i = 0; i < this->nb_dof_L_local; i++)
      IndexLambda(this->dof_local_L(i)) = i;
    
    VarComputationProblemLocal var_local(var_computation, *this, level_vol,
					 IndexLocal, IndexLambda);

    int rank_proc = var.GetRankProcMode();
    
    if (static_condensation)
      {
	// case where a static condensation is performed
	var_computation.SetLeafStaticCondensation(true);

	cond_solver = var_computation.GetNewCondensationSolver(T(0));

	CondensationBlockSolver_Fem<T>& solver_fem
	  = static_cast<CondensationBlockSolver_Fem<T>& >(*cond_solver);

	solver_fem.SetStoreBlock(true);

	// condensed matrix (to be solved)
	int nL = this->nb_dof_L_local;
	DistributedMatrix<T, General, ArrayRowSparse> mat_direct(nL, nL);
#ifdef SELDON_WITH_MPI
	mat_direct.Init(mat_direct_info);
#endif

        if (rank_proc == 0)
          cout << "The local finite element matrix is computed..." << endl;

	Matrix<T> mat_elt;
	nat_mat.SetCoefMass(a);
	nat_mat.SetCoefDamping(b);
	nat_mat.SetCoefStiffness(c);
	AssembleMatrix(mat_direct, mat_elt, nat_mat, var_local, *cond_solver, 0, 0);
        
        if (rank_proc == 0)
          cout << "The local finite element matrix is factorized..." << endl;

	// we factorized the condensed matrix
        mat_lu.SetPrintLevel(2);
        mat_lu.SetPivotThreshold(0.0);
	mat_lu.Factorize(mat_direct);

        mat_lu.SetPrintLevel(0);
	var_computation.SetLeafStaticCondensation(false);

        if (rank_proc == 0)
          cout << "Factorization completed" << endl;
      }    
    else
      {
	if (var.GetNbProcPerMode() > 1)
	  {
	    cout << "Case not implemented : use static condensation in parallel" << endl;
	    abort();
	  }
	
	// case without static condensation
	cond_solver = new CondensationBlockSolver_Base<T>();
	
        if (rank_proc == 0)
          cout << "The local finite element matrix is computed" << endl;

	// the finite element matrix is computed
	DistributedMatrix<T, General, ArrayRowSparse> mat_direct(N, N);    
	Matrix<T> mat_elt;
	nat_mat.SetCoefMass(a);
	nat_mat.SetCoefDamping(b);
	nat_mat.SetCoefStiffness(c);
	AssembleMatrix(mat_direct, mat_elt, nat_mat, var_local, *cond_solver, 0, 0);
	
        if (rank_proc == 0)
          cout << "The local finite element matrix is factorized" << endl;

	// and factorized
        mat_lu.SetPrintLevel(2);
	mat_lu.Factorize(mat_direct);

        mat_lu.SetPrintLevel(0);
      }    
  }

  
  //! factorizes the linear system a M + b K 
  void LocalOdeSchemeHDG
  ::FactorizeOperatorDhPlusGammaKh(const Real_wp& a, const Real_wp& b, const Real_wp& c)
  {
    FactorizeGeneric(a, b, c, cond_real_solver, nat_mat_real, mat_lu_real);    
  }


  //! solves the system (a M + b K) y = x (real or complex)
  template<class T>
  void LocalOdeSchemeHDG::
  SolveGeneric(CondensationBlockSolver_Base<T>*& cond_solver,
	       GlobalGenericMatrix<T>& nat_mat, All_MatrixLU<T>& mat_lu, Vector<T>& x,
	       Vector<T>& F_g, Vector<T>& L_g, Vector<T>& x_g)
  {
    VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(global_ode);
    MatrixVectorProductLevel& level_vol = var_time.GetTimeLevelDistribution();
    //DistributedProblem_Base& var = var_time.GetDistributedProblem();
    
    for (int i = 0; i < this->nb_dof_L_extended; i++)
      SetComplexZero(L_g(this->dof_extended_L(i)));
    
    if (static_condensation)
      {
	// modification de F
	int offset_u = L_g.GetM();
	for (int i = 0; i < this->nb_dof_uv_extended_nonloc; i++)
	  SetComplexZero(F_g(offset_u + this->dof_extended_uv_nonloc(i)));
	
	for (int i = 0; i < this->nb_dof_uv_proche; i++)
	  F_g(offset_u + this->dof_proche_uv(i)) = x(i);
	
	CondensationBlockSolver_Fem<T>& solver_fem
	  = static_cast<CondensationBlockSolver_Fem<T>& >(*cond_solver);
	
	solver_fem.GetVolumeLevel() = level_vol;
	solver_fem.GetVolumeLevel().SetLevel(this->num_level_vol);
	solver_fem.ModifyRhsStaticCondensation(SeldonNoTrans, F_g, nat_mat);

	this->AssembleLocal(global_ode, F_g);
	
	//this->dof_local_L.WriteText("dofL.dat");
	//this->F_glob.WriteText("F.dat");
	
	// calcul de Lambda
	Vector<T> Lambda(this->nb_dof_L_local);
	for (int i = 0; i < this->nb_dof_L_local; i++)
	  Lambda(i) = F_g(this->dof_local_L(i));
	
	mat_lu.Solve(Lambda);

        //Lambda.WriteText("L.dat");

	for (int i = 0; i < this->nb_dof_L_local; i++)
	  F_g(this->dof_local_L(i)) = Lambda(i);
	
	// recomposition de la solution x
	solver_fem.RecomposeSolution(SeldonNoTrans, F_g, nat_mat);
      }
    else
      {
	for (int i = 0; i < this->nb_dof_uv_extended_nonloc; i++)
	  SetComplexZero(x_g(this->dof_extended_uv_nonloc(i)));

	for (int i = 0; i < this->nb_dof_uv_proche; i++)
	  x_g(this->dof_proche_uv(i)) = x(i);

	int N = this->nb_dof_uv_local + this->nb_dof_L_local;
	Vector<T> F(N);
	F.Zero();
	
	for (int i = 0; i < this->nb_dof_uv_local; i++)
	  F(this->nb_dof_L_local + i) = x_g(this->dof_local_uv(i));
	
	mat_lu.Solve(F);
	
	for (int i = 0; i < this->nb_dof_uv_local; i++)
	  x_g(this->dof_local_uv(i)) = F(this->nb_dof_L_local + i);
	
	for (int i = 0; i < this->nb_dof_L_local; i++)
	  L_g(this->dof_local_L(i)) = F(i);	
      }
    
    // and close dofs
    var_time.MltAddLambdaHDG(nat_mat.GetCoefStiffness() / nat_mat.GetCoefMass(),
			     this->num_level_vol+2, L_g, x_g);
	
    this->SolveOperatorM(this->global_ode, x_g, this->dof_extended_uv_nonloc, true);
    
    for (int i = 0; i < this->nb_dof_uv_proche; i++)
      x(i) = x_g(this->dof_proche_uv(i));
    
    //x.WriteText("sol.dat");
  }


  //! Solves (a M + c K) y = f
  void LocalOdeSchemeHDG
  ::SolveOperatorDhPlusGammaKh(const Real_wp& t, const VectReal_wp& f, VectReal_wp& y)
  {
    y = f;
    SolveGeneric(cond_real_solver, nat_mat_real, mat_lu_real, y,
		 this->F_glob, this->L_glob, this->x_glob);
  }
  

  //! factorize alpha M + gamma K
  void LocalOdeSchemeHDG
  ::FactorizeOperatorComplex(const VectComplex_wp& alpha, const VectComplex_wp& beta,
			     const VectComplex_wp& gamma)
  {
    int N = alpha.GetM();
    cond_cplx_solver.Reallocate(N);
    mat_lu_cplx.Reallocate(N);
    nat_mat_cplx.Reallocate(N);
    for (int i = 0; i < N; i++)
      FactorizeGeneric(alpha(i), beta(i), gamma(i), cond_cplx_solver(i),
		       nat_mat_cplx(i), mat_lu_cplx(i));    
  }


  //! Solves (a M + c K) y = f (where a = alpha(n), c = gamma(n))
  void LocalOdeSchemeHDG
  ::SolveOperatorComplex(const Real_wp& t, const VectComplex_wp& X, VectComplex_wp& Y, int n)
  {
    if (F_cplx.GetM() != this->F_glob.GetM())
      {
	x_cplx.Nullify();
	L_cplx.Nullify();
	F_cplx.Reallocate(this->F_glob.GetM());
	x_cplx.SetData(this->x_glob.GetM(), &F_cplx(this->L_glob.GetM()));
	L_cplx.SetData(this->L_glob.GetM(), F_cplx.GetData());
      }
    
    Y = X;
    SolveGeneric(cond_cplx_solver(n), nat_mat_cplx(n), mat_lu_cplx(n), Y,
		 F_cplx, L_cplx, x_cplx);
  }
  

  void LocalOdeSchemeHDG::GetMemoryUsed(map<string, size_t>& var) const
  {
    if (cond_real_solver != NULL)
      {
        var["CondensationSolverFine"] = cond_real_solver->GetMemorySize();
        var["LuSolverFine"] = mat_lu_real.GetMemorySize();
      }
    
    for (int i = 0; i < cond_cplx_solver.GetM(); i++)
      if (cond_cplx_solver(i) != NULL)
        {
          var["CondensationSolverFine"] += cond_cplx_solver(i)->GetMemorySize();
          var["LuSolverFine"] += mat_lu_cplx(i).GetMemorySize();
        }
  }
  
  
  /***********************
   * LocalTimeSteppingRK *
   ***********************/


  //! default constructor
  template<class T>
  LocalTimeSteppingRK<T>::LocalTimeSteppingRK()
    : coarse_ode(x_glob, y_glob, L_glob, F_glob)
  {
    dt = Real_wp(0);
    fine_ode = NULL;
    implicit_scheme = NULL;
  }


  //! destructor
  template<class T>
  LocalTimeSteppingRK<T>::~LocalTimeSteppingRK()
  {
    x_glob.Nullify();
    L_glob.Nullify();
    if (fine_ode != NULL)
      delete fine_ode;
    
    if (implicit_scheme != NULL)
      delete implicit_scheme;
  }

  
  //! stability function = explicit scheme
  template<class T>
  UnivariatePolynomial<Real_wp> LocalTimeSteppingRK<T>::GetStabilityFunction() const
  {
    return explicit_scheme.GetStabilityFunction();
  }


  //! sets the time schemes to use (explicit and implicit part)
  template<class T>
  void LocalTimeSteppingRK<T>::SetInputData(int r, const Vector<string>& param)
  {
    if (param.GetM() <= 4)
      {
	cout << "Provide more parameters for LOCAL_IMPLICIT_RK" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }
    
    // explicit scheme is initialized
    int nb = 4;
    int order_exp = r;
    int s_exp = to_num<int>(param(1));
    explicit_scheme.SetOrder(order_exp, s_exp, false);

    // implicit scheme is initialized
    this->local_implicit = true;
    if (implicit_scheme != NULL)
      delete implicit_scheme;

    if (param(2) == "Pade")
      {
	int order_imp = to_num<int>(param(3));
	implicit_scheme = new PadeScheme_Iterator<T>(order_imp, true);
      }
    else if (param(2) == "LinearSdirk")
      {
	int order_imp = to_num<int>(param(3));
	int s_imp = to_num<int>(param(4));
	implicit_scheme =
	  new LinearSdirkScheme_Iterator<T>(order_imp, s_imp,
					    LinearSdirkScheme_Iterator<T>::STABLE_WEIGHTS);

	nb = 5;
      }
    else if (param(2) == "Explicit")
      {
	this->local_implicit = false;
	implicit_scheme = NULL;
	this->p = to_num<int>(param(3));
      }
    else
      {
	cout << "Scheme used on the fine region unknown" << endl;
	abort();
      }

    this->SetInputSplitting(param, nb);
  }
  

  //! sets the initial condition
  template<>
  void LocalTimeSteppingRK<Complex_wp>::SetInitialCondition(const Real_wp& t, const Real_wp& dt0,
							    Vector<Complex_wp>& U0, VirtualOdeSystem<Complex_wp>& sys)
  {
    cout << "Not available" << endl;
    abort();
  }
  

  //! sets the initial condition
  template<>
  void LocalTimeSteppingRK<Real_wp>::SetInitialCondition(const Real_wp& t, const Real_wp& dt0,
							 Vector<Real_wp>& U0, VirtualOdeSystem<Real_wp>& sys)
  {
    glob_chrono.SetMessage("Solver_Local", " solve the local ODE");
    glob_chrono.SetMessage("wtild_Local", " compute w-tilde local");
    glob_chrono.SetMessage("Source_Local", " compute the source");
    
    dt = dt0;
    
    const VectReal_wp& points_ci = explicit_scheme.GetPoints();
    DerMat = explicit_scheme.GetDerivativeMatrix();
    for (int l = 1; l < DerMat.GetM(); l++)
      {
	Real_wp coef = Real_wp(1) / pow(dt, l);
	for (int j = 0; j < DerMat.GetN(); j++)
	  DerMat(l, j) *= coef;
      }

    fine_ode = new LocalOdeSchemeHDG(points_ci, dt, sys, dof_fine_among_extended,
				     x_glob, y_glob, L_glob, F_glob);
    
    // we create num_explicit
    {
      VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
      DistributedProblem_Base& var = var_time.GetDistributedProblem();
      VarComputationProblem& var_p = dynamic_cast<VarComputationProblem&>(var);
      int nb_elt = var_p.GetNbElt();

      Vector<bool> ImplicitElement(nb_elt);
      ImplicitElement.Fill(false);
      for (int i = 0; i < this->num_implicit.GetM(); i++)
	ImplicitElement(this->num_implicit(i)) = true;

      this->num_explicit.Reallocate(nb_elt - this->num_implicit.GetM());
      int nb = 0;
      for (int i = 0; i < nb_elt; i++)
	if (!ImplicitElement(i))
	  {
	    this->num_explicit(nb) = i;
	    nb++;
	  }
    }

    //DISP(this->num_explicit);
    //DISP(this->num_implicit);
    
    // we initialize the coarse and fine region
    fine_ode->Init(this->num_implicit, sys, num_level_volume, num_level_surface, true);
    coarse_ode.Init(this->num_explicit, sys, num_level_volume, num_level_surface);
    
    int Nimp = fine_ode->GetNbLocalDof();
    int Nfar = coarse_ode.GetNbLocalDof();
    int Nfar_ext = coarse_ode.GetNbDofProche();
    
    {
      // we find dofs in common between the coarse and fine region (extended zones)
      const Vector<int>& dof_far = coarse_ode.GetDofProche();
      const Vector<int>& dof_loc = fine_ode->GetDofProche();
      
      coarse_ode.FindCommonDofs(dof_far, dof_loc, dof_common_fine, dof_common_coarse);

      Vector<bool> DofCommon(Nfar_ext);
      DofCommon.Fill(false);
      for (int i = 0; i < dof_common_coarse.GetM(); i++)
	DofCommon(dof_common_coarse(i)) = true;

      dof_original_coarse.Reallocate(Nfar_ext - dof_common_coarse.GetM());
      int nb = 0;
      for (int i = 0; i < Nfar_ext; i++)
	if (!DofCommon(i))
	  {
	    dof_original_coarse(nb) = i;
	    nb++;
	  }
      
      // and conversion between implicit dofs and extended dofs
      Vector<int> all_num;
      const Vector<int>& dof_imp = fine_ode->GetLocalDof();
      coarse_ode.FindCommonDofs(dof_loc, dof_imp, all_num, dof_fine_among_extended);
    }
    
    // we copy the initial condition
    Un = U0;
    U0.Clear();

    int Nproche = fine_ode->GetNbDofProche();
    int nb_dof_uv = Un.GetM();
    
    // vectors needed for the time scheme
    EvalF.Reallocate(points_ci.GetM());
    EvalFproche.Reallocate(points_ci.GetM());

    // for dense storage of EvalF/EvalFproche
    //for (int k = 0; k < points_ci.GetM(); k++)
    //  EvalF(k).Reallocate(nb_dof_uv);
    //for (int k = 0; k < points_ci.GetM(); k++)
    //  EvalFproche(k).Reallocate(Nproche);

    w_tild.Reallocate(nb_dof_uv);
    y_imp.Reallocate(Nimp); y_far.Reallocate(Nfar);
    AP_y.Reallocate(Nproche); AImP_y.Reallocate(Nfar_ext);
    AP_y.Zero(); AImP_y.Zero();

    if (!this->local_implicit)
      {
	explicit_scheme.SetInitialCondition(t, dt/this->p, AP_y, *fine_ode);
	AP_y.Reallocate(Nproche);
      }
    else
      {
	implicit_scheme->SetInitialCondition(t, dt, AP_y, *fine_ode);
	AP_y.Reallocate(Nproche);
      }

    const UnivariatePolynomial<Real_wp>& P = explicit_scheme.GetStabilityFunction();
    int orderN = P.GetOrder();

    w_n_i.Reallocate(orderN); wp_i.Reallocate(orderN);
    for (int i = 0; i < orderN; i++)
      {
	w_n_i(i).Reallocate(Nfar_ext);
	wp_i(i).Reallocate(dof_common_fine.GetM());
	w_n_i(i).Zero(); wp_i(i).Zero();
      }    
  }
  

  //! U^n+1 is computed from U^n and then replaces U^n
  template<>
  void LocalTimeSteppingRK<Complex_wp>::Advance(const Real_wp& t, int n, VirtualOdeSystem<Complex_wp>& sys)
  {
    cout << "Not implemented for complex numbers" << endl;
    abort();
  }


  //! U^n+1 is computed from U^n and then replaces U^n
  template<>
  void LocalTimeSteppingRK<Real_wp>::Advance(const Real_wp& t, int n, VirtualOdeSystem<Real_wp>& sys)
  {
    const VectReal_wp& points_ci = explicit_scheme.GetPoints();
    
    int Nproche = fine_ode->GetNbDofProche();
    int Nimp = fine_ode->GetNbLocalDof();
    int Nfar = coarse_ode.GetNbLocalDof();
    int Nfar_ext = coarse_ode.GetNbDofProche();

    /* #ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
    int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    string suffix;
    if (nb_proc > 1)
      suffix = "_" + to_str(rank_proc) + "N" + to_str(nb_proc);
    */
    
    const Vector<int>& dof_implicit = fine_ode->GetLocalDof();
    const Vector<int>& dof_explicit = coarse_ode.GetLocalDof();
    const Vector<int>& dof_proche = fine_ode->GetDofProche();
    const Vector<int>& dof_explicit_extended = coarse_ode.GetDofProche();

    const UnivariatePolynomial<Real_wp>& P = explicit_scheme.GetStabilityFunction();
    int orderN = P.GetOrder();

    glob_chrono.Start("Source_Local");
    // source is evaluated at quadrature points
    bool presence_source = false;
    if (t < sys.GetFinalTimeSource())
      {
	presence_source = true;
	const VectReal_wp& points_ci = explicit_scheme.GetPoints();
	for (int i = 0; i < points_ci.GetM(); i++)
	  {
            w_tild.Zero();
	    sys.AddPrimitiveTimeSource(Real_wp(1), t + points_ci(i)*dt, 0, w_tild);
            Copy(w_tild, EvalF(i));
	    
	    // local part is extracted
	    for (int j = 0; j < Nproche; j++)
              AP_y(j) = w_tild(dof_proche(j));	    
            
            Copy(AP_y, EvalFproche(i));
	  }    
      }
    glob_chrono.Stop("Source_Local");

    glob_chrono.Start("wtild_Local");
    //Un.WriteText("Un" + suffix + ".dat");
    // vectors w_n_j are computed
    w_tild = Un;
    VectReal_wp coef_alpha(orderN);
    for (int j = 0; j < orderN; j++)
      {
	// we compute P w and (I-P) w
	for (int i = 0; i < Nimp; i++)
	  y_imp(i) = w_tild(dof_implicit(i));
	
	for (int i = 0; i < Nfar; i++)
	  y_far(i) = w_tild(dof_explicit(i));

	// we compute AP w and A(I-P) w
	fine_ode->ApplyOperatorAP(y_imp, AP_y, sys, false);
	coarse_ode.ApplyOperatorAP(y_far, AImP_y, sys, false);

	// vectors w_n_j are stored
	w_n_i(j) = AImP_y;
	Real_wp coef = (j+1)*P(j+1);
	Mlt(coef, w_n_i(j));
	coef_alpha(j) = Real_wp(1)/(j+1);
	
	// we store wp_i
	for (int i = 0; i < dof_common_coarse.GetM(); i++)
	  wp_i(j)(i) = w_n_i(j)(dof_common_coarse(i));

	//wp_i(j).WriteText("wp" + to_str(j) + ".dat");
	
	// then w_tild = A w
	w_tild.Zero();
	for (int i = 0; i < Nfar_ext; i++)
	  w_tild(dof_explicit_extended(i)) = AImP_y(i);

	//DISP(j); w_tild.WriteText("AImP" + to_str(j) + suffix + ".dat");
	for(int i = 0; i < Nproche; i++)
	  w_tild(dof_proche(i)) += AP_y(i);

	//w_tild.WriteText("AP" + to_str(j) + suffix + ".dat");
	
	// the source is added
	if (presence_source)
	  if (j < orderN-1)
	    for (int i = 0; i < points_ci.GetM(); i++)
	      Add(DerMat(j, i), EvalF(i), w_tild);

	// we solve by the mass matrix
	sys.SolveMassMatrix(w_tild);
      }	
    glob_chrono.Stop("wtild_Local");
    
    // we extract the solution on the fine region
    for(int i = 0; i < Nproche; i++)
      AP_y(i) = Un(dof_proche(i));

    fine_ode->SetSource(t, wp_i, dof_common_fine, presence_source, EvalFproche);

    // then we chose which local time stepping is used for the fine region
    glob_chrono.Start("Solver_Local");
    if (this->local_implicit)
      {
	implicit_scheme->GetIterate() = AP_y;

	implicit_scheme->Advance(t, 0, *fine_ode);

	AP_y = implicit_scheme->GetIterate();
      }
    else
      {
	Real_wp dtau = dt/this->p;
	
	// the explicit scheme p times
	explicit_scheme.GetIterate() = AP_y;
	
	for(int m = 0; m < p; m++)
	  explicit_scheme.Advance(t + m*dtau, m, *fine_ode);
	
	AP_y = explicit_scheme.GetIterate();
      }

    // the solution in the fine zone is copied
    for (int i = 0; i < Nproche; i++)
      Un(dof_proche(i)) = AP_y(i);

    //Un.WriteText("UnProche" + suffix + ".dat");
    
    // the coarse zone (not common with extended fine zone) is advanced
    const Matrix<Real_wp>& bi = explicit_scheme.GetWeights();
    Vector<int> num(bi.GetM());
    num.Zero();
    for (int i = 0; i < Nfar_ext; i++)
      {
	int ig = dof_explicit_extended(i);
	Real_wp vloc = coef_alpha(orderN-1)*w_n_i(orderN-1)(i);
	for (int j = orderN-2; j >= 0; j--)
	  vloc = dt*vloc + coef_alpha(j)*w_n_i(j)(i);
	
	if (presence_source)
	  for (int j = 0; j < bi.GetN(); j++)            
	    {
              while ((num(j) < EvalF(j).GetM()) && (EvalF(j).Index(num(j)) < ig))
                num(j)++;
              
              if ((num(j) < EvalF(j).GetM()) && (EvalF(j).Index(num(j)) == ig))
                vloc += bi(0, j)*EvalF(j)(ig);
            }
        
	AImP_y(i) = vloc;
      }

    coarse_ode.SolveOperatorM(sys, AImP_y, coarse_ode.GetDofProche(), false);

    //dof_proche.WriteText("dof_proche.dat");

    for (int i = 0; i < dof_original_coarse.GetM(); i++)
      {
	int i0 = dof_original_coarse(i);
	int ig = dof_explicit_extended(i0);
	Un(ig) += dt*AImP_y(i0);	
      }

    glob_chrono.Stop("Solver_Local");
    
    //Un.WriteText("Unext" + suffix + ".dat");
    
    //AP_y.WriteText("yp.dat");
    //int test_input; cout << "Waiting final" << endl; cin >> test_input;
  }


  //! changes the time step
  template<class T>
  void LocalTimeSteppingRK<T>::ChangeTimeStep(const Real_wp& dt0)
  {
    dt = dt0;
  }
    

  //! returns the main iterate U^n
  template<class T>
  Vector<T>& LocalTimeSteppingRK<T>::GetIterate()
  {
    return Un;
  }


  //! returns the main iterate U^n
  template<class T>
  const Vector<T>& LocalTimeSteppingRK<T>::GetIterate() const
  {
    return Un;
  }


  //! returns the number of iterates
  template<class T>
  int LocalTimeSteppingRK<T>::GetNumberOfIterates() const
  {
    return 1;
  }


  //! returns the main iterate U^n
  template<class T>
  Vector<T>& LocalTimeSteppingRK<T>::GetIterate(int k)
  {
    return Un;
  }


  //! clears vectors stored
  template<class T>
  void LocalTimeSteppingRK<T>::Clear()
  {
    Un.Clear();
  }


  //! clears vectors stored and Uf = U^n
  template<class T>
  void LocalTimeSteppingRK<T>::ClearFirst(Vector<T>& Uf)
  {
    EvalF.Clear(); EvalFproche.Clear();
    Uf = Un;
    Clear();
  }
  
  
  template<class T>
  size_t LocalTimeSteppingRK<T>::GetMemorySize() const
  {
    size_t taille = 0;
    taille += Un.GetMemorySize() + Seldon::GetMemorySize(EvalF)
      + Seldon::GetMemorySize(EvalFproche) + Seldon::GetMemorySize(w_n_i)
      + Seldon::GetMemorySize(wp_i) + w_tild.GetMemorySize()
      + y_imp.GetMemorySize() + y_far.GetMemorySize() + AP_y.GetMemorySize()
      + AImP_y.GetMemorySize() + x_glob.GetMemorySize() + y_glob.GetMemorySize()
      + L_glob.GetMemorySize() + F_glob.GetMemorySize();
    
    return taille;
  }
  

  template<class T>
  void LocalTimeSteppingRK<T>::GetMemoryUsed(map<string, size_t>& var) const
  {
    fine_ode->GetMemoryUsed(var);
  }


  /*************************
   * LocalTimeSteppingImex *
   *************************/


  //! default constructor
  template<class T>
  LocalImexRK<T>::LocalImexRK()
    : coarse_ode(x_glob, y_glob, L_glob, F_glob)
  {
    dt = Real_wp(0);
    fine_ode = NULL;
  }


  //! destructor
  template<class T>
  LocalImexRK<T>::~LocalImexRK()
  {
    x_glob.Nullify();
    L_glob.Nullify();
    if (fine_ode != NULL)
      delete fine_ode;
  }

  
  //! stability function = explicit scheme
  template<class T>
  UnivariatePolynomial<Real_wp> LocalImexRK<T>::GetStabilityFunction() const
  {
    // PolK = polynomial associated with each ki
    Vector<UnivariatePolynomial<Real_wp> > PolK(bi.GetM());
    UnivariatePolynomial<Real_wp> z, one, P;
    one.SetOrder(0); one(0) = 1.0;
    z.SetOrder(1); z(0) = 0.0; z(1) = 1.0;
    
    P = one;
    for (int i = 0; i < bi.GetM(); i++)
      {
        PolK(i) = one;
        for (int j = 0; j < i; j++)
          Add(Aexp(i, j), PolK(j), PolK(i));
        
        PolK(i) *= z;
        Add(bi(i), PolK(i), P);
      }
    
    return P;
  }


  //! sets the time schemes to use (explicit and implicit part)
  template<class T>
  void LocalImexRK<T>::SetInputData(int r, const Vector<string>& param)
  {
    if (param.GetM() <= 1)
      {
	cout << "Provide more parameters for LOCAL_IMEX_RK" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }
    
    int nb = 1;    
    if (r == 3)
      {
        ci.Reallocate(4);
        ci(0) = Real_wp(0);
        ci(1) = Real_wp(1767732205903)/Real_wp(2027836641118);
        ci(2) = Real_wp(3)/Real_wp(5);
        ci(3) = Real_wp(1);
        
        Aexp.Reallocate(4, 4); Aexp.Zero();
        Aexp(1, 0) = ci(1);
        Aexp(2, 0) = Real_wp(5535828885825)/Real_wp(10492691773637);
        Aexp(2, 1) = ci(2) - Aexp(2, 0);
        Aexp(3, 0) = Real_wp(6485989280629)/Real_wp(16251701735622);
        Aexp(3, 1) = Real_wp(-4246266847089)/Real_wp(9704473918619);
        Aexp(3, 2) = ci(3) - Aexp(3, 0) - Aexp(3, 1);
        
        bi.Reallocate(4);
        bi(0) = Real_wp(1471266399579)/Real_wp(7840856788654);
        bi(1) = Real_wp(-4482444167858)/Real_wp(7529755066697);
        bi(2) = Real_wp(11266239266428)/Real_wp(11593286722821);
        bi(3) = Real_wp(1767732205903)/Real_wp(4055673282236);
        
        Aimp.Reallocate(4, 4); Aimp.Zero();
        Aimp(1, 0) = Aexp(1, 0)/2; Aimp(1, 1) = Aimp(1, 0);
        Aimp(2, 2) = Aimp(1, 1); Aimp(3, 3) = Aimp(1, 1);
        Aimp(2, 0) = Real_wp(2746238789719)/Real_wp(10658868560708);
        Aimp(2, 1) = ci(2) - Aimp(2, 2) - Aimp(2, 0);
        Aimp(3, 0) = Real_wp(1471266399579)/Real_wp(7840856788654);
        Aimp(3, 1) = Real_wp(-4482444167858)/Real_wp(7529755066697);
        Aimp(3, 2) = ci(3) - Aimp(3, 3) - Aimp(3, 0) - Aimp(3, 1);
      }
    else if (r == 4)
      {
        if (param(1) == "Asher")
          {
            ci.Reallocate(5);
            ci(0) = Real_wp(0);
            ci(1) = Real_wp(1)/Real_wp(2);
            ci(2) = Real_wp(2)/Real_wp(3);
            ci(3) = ci(1);
            ci(4) = Real_wp(1);
            
            Aexp.Reallocate(5, 5); Aexp.Zero();
            Aexp(1, 0) = ci(1);
            Aexp(2, 0) = Real_wp(11)/Real_wp(18); Aexp(2, 1) = Real_wp(1)/Real_wp(18);
            Aexp(3, 0) = Real_wp(5)/Real_wp(6); Aexp(3, 1) = -Aexp(3, 0); Aexp(3, 2) = ci(3);
            Aexp(4, 0) = Real_wp(1)/Real_wp(4); Aexp(4, 1) = Real_wp(7)/Real_wp(4);
            Aexp(4, 2) = Real_wp(3)/Real_wp(4); Aexp(4, 3) = -Aexp(4, 1);
            
            bi.Reallocate(5);
            for (int i = 0; i < 5; i++)
              bi(i) = Aexp(4, i);
            
            Aimp.Reallocate(5, 5); Aimp.Zero();
            Aimp(4, 4) = ci(1); Aimp(1, 1) = ci(1); Aimp(2, 2) = ci(1); Aimp(3, 3) = ci(1);
            Aimp(2, 1) = Real_wp(1)/Real_wp(6);
            Aimp(3, 1) = Real_wp(-1)/Real_wp(2); Aimp(3, 2) = Real_wp(1)/Real_wp(2);
            Aimp(4, 1) = Real_wp(3)/Real_wp(2); Aimp(4, 2) = -Aimp(4, 1);
            Aimp(4, 3) = Aimp(4, 4);
            
            nb++;
          }
        else if (param(1) == "Calvo")
          {
            ci.Reallocate(6);
            ci(0) = Real_wp(0);
            ci(1) = Real_wp(1)/Real_wp(4);
            ci(2) = Real_wp(3)/Real_wp(4);
            ci(3) = Real_wp(11)/Real_wp(20);
            ci(4) = Real_wp(1)/Real_wp(2);
            ci(5) = Real_wp(1);

            Aexp.Reallocate(6, 6); Aexp.Zero();
            Aexp(1, 0) = ci(1);
            Aexp(2, 0) = Real_wp(-1)/Real_wp(4); Aexp(2, 1) = Real_wp(1);
            Aexp(3, 0) = Real_wp(-13)/Real_wp(100); Aexp(3, 1) = Real_wp(43)/Real_wp(75);
            Aexp(3, 2) = ci(3) - Aexp(3, 0) - Aexp(3, 1);
            Aexp(4, 0) = Real_wp(-6)/Real_wp(85);
            Aexp(4, 1) = Real_wp(42)/Real_wp(85); Aexp(4, 2) = Real_wp(179)/Real_wp(1360);
            Aexp(4, 3) = ci(4) - Aexp(4, 0) - Aexp(4, 1) - Aexp(4, 2);
            Aexp(5, 1) = Real_wp(79)/Real_wp(24); Aexp(5, 2) = Real_wp(-5)/Real_wp(8);
            Aexp(5, 3) = Real_wp(25)/Real_wp(2); 
            Aexp(5, 4) = ci(5) - Aexp(5, 0) - Aexp(5, 1) - Aexp(5, 2) - Aexp(5, 3);

            bi.Reallocate(6);
            bi(0) = Real_wp(0); bi(1) = Real_wp(25)/Real_wp(24);
            bi(2) = Real_wp(-49)/Real_wp(48);
            bi(3) = Real_wp(125)/Real_wp(16);
            bi(4) = Real_wp(-85)/Real_wp(12);
            bi(5) = Real_wp(1) - bi(0) - bi(1) - bi(2) - bi(3) - bi(4);

            Aimp.Reallocate(6, 6); Aimp.Zero();
            Aimp(1, 1) = ci(1);
            Aimp(2, 2) = Aimp(1, 1);  Aimp(3, 3) = Aimp(1, 1); Aimp(4, 4) = Aimp(1, 1); Aimp(5, 5) = Aimp(1, 1);
            Aimp(2, 1) = Real_wp(1)/Real_wp(2);
            Aimp(3, 1) = Real_wp(17)/Real_wp(20);
            Aimp(3, 2) = ci(3) - Aimp(3, 1) - Aimp(3, 3);
            Aimp(4, 1) = Real_wp(371)/Real_wp(1360); Aimp(4, 2) = Real_wp(-137)/Real_wp(2720);
            Aimp(4, 3) = ci(4) - Aimp(4, 4) - Aimp(4, 1) - Aimp(4, 2);
            for (int i = 1; i < 6; i++)
              Aimp(5, i) = bi(i);
            
            nb++;
          }
        else
          {
            ci.Reallocate(6);
            ci(0) = Real_wp(0);
            ci(1) = Real_wp(1)/Real_wp(2);
            ci(2) = Real_wp(83)/Real_wp(250);
            ci(3) = Real_wp(31)/Real_wp(50);
            ci(4) = Real_wp(17)/Real_wp(20);
            ci(5) = Real_wp(1);
            
            Aexp.Reallocate(6, 6); Aexp.Zero();
            Aexp(1, 0) = ci(1);
            Aexp(2, 0) = Real_wp(13861)/Real_wp(62500);
            Aexp(2, 1) = ci(2) - Aexp(2, 0);
            Aexp(3, 0) = Real_wp(-116923316275)/Real_wp(2393684061468);
            Aexp(3, 1) = Real_wp(-2731218467317)/Real_wp(15368042101831);
            Aexp(3, 2) = ci(3) - Aexp(3, 0) - Aexp(3, 1);
            Aexp(4, 0) = Real_wp(-451086348788)/Real_wp(2902428689909);
            Aexp(4, 1) = Real_wp(-2682348792572)/Real_wp(7519795681897);
            Aexp(4, 2) = Real_wp(12662868775082)/Real_wp(11960479115383);
            Aexp(4, 3) = ci(4) - Aexp(4, 0) - Aexp(4, 1) - Aexp(4, 2);
            Aexp(5, 0) = Real_wp(647845179188)/Real_wp(3216320057751);
            Aexp(5, 1) = Real_wp(73281519250)/Real_wp(8382639484533);
            Aexp(5, 2) = Real_wp(552539513391)/Real_wp(3454668386233);
            Aexp(5, 3) = Real_wp(3354512671639)/Real_wp(8306763924573);
            Aexp(5, 4) = ci(5) - Aexp(5, 0) - Aexp(5, 1) - Aexp(5, 2) - Aexp(5, 3);
            
            bi.Reallocate(6);
            bi(0) = Real_wp(82889)/Real_wp(524892);
            bi(1) = Real_wp(0); bi(2) = Real_wp(15625)/Real_wp(83664);
            bi(3) = Real_wp(69875)/Real_wp(102672); bi(5) = Real_wp(1)/Real_wp(4);
            bi(4) = Real_wp(1) - bi(0) - bi(1) - bi(2) - bi(3) - bi(5);
            
            Aimp.Reallocate(6, 6); Aimp.Zero();
            Aimp(1, 0) = Aexp(1, 0)/2; Aimp(1, 1) = Aimp(1, 0);
            Aimp(2, 2) = Aimp(1, 1);  Aimp(3, 3) = Aimp(1, 1); Aimp(4, 4) = Aimp(1, 1); Aimp(5, 5) = Aimp(1, 1);
            Aimp(2, 0) = Real_wp(8611)/Real_wp(62500);
            Aimp(2, 1) = ci(2) - Aimp(2, 2) - Aimp(2, 0);
            Aimp(3, 0) = Real_wp(5012029)/Real_wp(34652500);
            Aimp(3, 1) = Real_wp(-654441)/Real_wp(2922500);
            Aimp(3, 2) = ci(3) - Aimp(3, 3) - Aimp(3, 0) - Aimp(3, 1);
            Aimp(4, 0) = Real_wp(15267082809)/Real_wp(155376265600);
            Aimp(4, 1) = Real_wp(-71443401)/Real_wp(120774400);
            Aimp(4, 2) = Real_wp(730878875)/Real_wp(902184768);
            Aimp(4, 3) = ci(4) - Aimp(4, 4) - Aimp(4, 0) - Aimp(4, 1) - Aimp(4, 2);
            Aimp(5, 0) = Real_wp(82889)/Real_wp(524892);
            Aimp(5, 1) = Real_wp(0); Aimp(5, 2) = Real_wp(15625)/Real_wp(83664);
            Aimp(5, 3) = Real_wp(69875)/Real_wp(102672);
            Aimp(5, 4) = ci(5) - Aimp(5, 5) - Aimp(5, 0) - Aimp(5, 1) - Aimp(5, 2) - Aimp(5, 3);
          }
      }
    else if (r == 5)
      {
        ci.Reallocate(8);
        ci(0) = Real_wp(0); ci(1) = Real_wp(41)/Real_wp(100);
        ci(2) = Real_wp(2935347310677)/Real_wp(11292855782102);
        ci(3) = Real_wp(1426016391358)/Real_wp(7196633302097);
        ci(4) = Real_wp(92)/Real_wp(100); ci(5) = Real_wp(24)/Real_wp(100);
        ci(6) = Real_wp(3)/Real_wp(5); ci(7) = Real_wp(1);

        Aexp.Reallocate(8, 8); Aexp.Zero();
        Aexp(1, 0) = ci(1); Aexp(2, 0) = Real_wp(367902744464)/Real_wp(2072280473677);
        Aexp(2, 1) = ci(2) - Aexp(2, 0);
        Aexp(3, 0) = Real_wp(1268023523408)/Real_wp(10340822734521); 
        Aexp(3, 1) = Real_wp(0); Aexp(3, 2) = ci(3) - Aexp(3, 0);
        Aexp(4, 0) = Real_wp(14463281900351)/Real_wp(6315353703477);
        Aexp(4, 1) = Real_wp(0); Aexp(4, 2) = Real_wp(66114435211212)/Real_wp(5879490589093);
        Aexp(4, 3) = ci(4) - Aexp(4, 2) - Aexp(4, 1) - Aexp(4, 0);
        Aexp(5, 0) = Real_wp(14090043504691)/Real_wp(34967701212078);
        Aexp(5, 1) = Real_wp(0); Aexp(5, 2) = Real_wp(15191511035443)/Real_wp(11219624916014);
        Aexp(5, 3) = Real_wp(-18461159152457)/Real_wp(12425892160975);
        Aexp(5, 4) = ci(5) - Aexp(5, 0) - Aexp(5, 1) - Aexp(5, 2) - Aexp(5, 3);
        Aexp(6, 0) = Real_wp(19230459214898)/Real_wp(13134317526959);
        Aexp(6, 1) = Real_wp(0); Aexp(6, 2) = Real_wp(21275331358303)/Real_wp(2942455364971);
        Aexp(6, 4) = Real_wp(-1)/Real_wp(8); Aexp(6, 5) = Real_wp(-1)/Real_wp(8);
        Aexp(6, 3) = ci(6) - Aexp(6, 0) - Aexp(6, 2) - Aexp(6, 4) - Aexp(6, 5);
        Aexp(7, 0) = Real_wp(-19977161125411)/Real_wp(11928030595625);
        Aexp(7, 1) = Real_wp(0); Aexp(7, 2) = Real_wp(-40795976796054)/Real_wp(6384907823539);
        Aexp(7, 3) = Real_wp(177454434618887)/Real_wp(12078138498510);
        Aexp(7, 4) = Real_wp(782672205425)/Real_wp(8267701900261);
        Aexp(7, 5) = Real_wp(-69563011059811)/Real_wp(9646580694205);
        Aexp(7, 6) = ci(7) - Aexp(7, 0) - Aexp(7, 2) - Aexp(7, 3) - Aexp(7, 4) - Aexp(7, 5);
        
        bi.Reallocate(8);
        bi(0) = Real_wp(-872700587467)/Real_wp(9133579230613);
        bi(1) = Real_wp(0); bi(2) = Real_wp(0);
        bi(3) = Real_wp(22348218063261)/Real_wp(9555858737531);
        bi(4) = Real_wp(-1143369518992)/Real_wp(8141816002931);
        bi(5) = Real_wp(-39379526789629)/Real_wp(19018526304540);
        bi(7) = Real_wp(41)/Real_wp(200);
        bi(6) = Real_wp(1) - bi(0) - bi(3) - bi(4) - bi(5) - bi(7);

        Aimp.Reallocate(8, 8); Aimp.Zero();
        Aimp(1, 0) = bi(7); Aimp(1, 1) = bi(7); Aimp(2, 2) = bi(7); Aimp(3, 3) = bi(7);
        Aimp(4, 4) = bi(7); Aimp(5, 5) = bi(7); Aimp(6, 6) = bi(7); Aimp(7, 7) = bi(7);
        Aimp(2, 0) = Real_wp(41)/Real_wp(400); Aimp(2, 1) = ci(2) - Aimp(2, 0) - Aimp(2, 2);
        Aimp(3, 0) = Real_wp(683785636431)/Real_wp(9252920307686);
        Aimp(3, 2) = ci(3) - Aimp(3, 0) - Aimp(3, 3);
        Aimp(4, 0) = Real_wp(3016520224154)/Real_wp(10081342136671);
        Aimp(4, 2) = Real_wp(30586259806659)/Real_wp(12414158314087);
        Aimp(4, 3) = ci(4) - Aimp(4, 0) - Aimp(4, 2) - Aimp(4, 4);
        Aimp(5, 0) = Real_wp(218866479029)/Real_wp(1489978393911);
        Aimp(5, 2) = Real_wp(638256894668)/Real_wp(5436446318841);
        Aimp(5, 3) = Real_wp(-1179710474555)/Real_wp(5321154724896);
        Aimp(5, 4) = ci(5) - Aimp(5, 0) - Aimp(5, 2) - Aimp(5, 3) - Aimp(5, 5);
        Aimp(6, 0) = Real_wp(1020004230633)/Real_wp(5715676835656);
        Aimp(6, 2) = Real_wp(25762820946817)/Real_wp(25263940353407);
        Aimp(6, 3) = Real_wp(-2161375909145)/Real_wp(9755907335909);
        Aimp(6, 4) = Real_wp(-211217309593)/Real_wp(5846859502534);
        Aimp(6, 5) = ci(6) - Aimp(6, 0) - Aimp(6, 2) - Aimp(6, 3) - Aimp(6, 4) - Aimp(6, 6);
        for (int i = 0; i < 8; i++)
          Aimp(7, i) = bi(i);
      }
    
    /*
    DISP(ci); DISP(bi); DISP(Aexp); DISP(Aimp);
    bi.WriteText("bi.dat"); ci.WriteText("ci.dat");
    Aexp.WriteText("Aexp.dat");
    Aimp.WriteText("Aimp.dat");
    
    // Cabane profile
    int N = 100;
    VectComplex_wp pts(6*N);
    for (int i = 0; i < N; i++)
      {
        Real_wp t = Real_wp(i+1)/N;
        pts(i) = Complex_wp(0, t);
        pts(6*N-1-i) = Complex_wp(0, -t);
        
        pts(N+i) = Complex_wp(-t, 1);
        pts(5*N-1-i) = Complex_wp(-t, -1);
        
        t = Real_wp(1) - t;
        Real_wp coef = t/10*(14.0-4.0*t);
        pts(2*N+i) = Complex_wp(t-2.0, coef);
        pts(4*N-1-i) = Complex_wp(t-2.0, -coef);
      }
    
    pts.Write("z.dat");
    
    // we compute the efficiency with bisection method
    UnivariatePolynomial<Real_wp> R = GetStabilityFunction();
    R.WriteText("R.dat");
    Real_wp a = 0, b = bi.GetM();
    while ((b-a) > 1e-12)
      {
        Real_wp c = (a+b)/2;
        bool stable_scheme = true;
        for (int j = 0; j < pts.GetM(); j++)
          {
            Real_wp G = abs(R.Evaluate(c*pts(j)));
            if (G > 1.0+1e-10)
              stable_scheme = false;
          }
        
        if (stable_scheme)
          a = c;
        else
          b = c;
      }

    DISP(a); DISP(a/bi.GetM());
    */
    
    this->SetInputSplitting(param, nb);
  }
  

  //! sets the initial condition
  template<>
  void LocalImexRK<Complex_wp>::SetInitialCondition(const Real_wp& t, const Real_wp& dt0,
							    Vector<Complex_wp>& U0, VirtualOdeSystem<Complex_wp>& sys)
  {
    cout << "Not available" << endl;
    abort();
  }
  

  //! sets the initial condition
  template<>
  void LocalImexRK<Real_wp>::SetInitialCondition(const Real_wp& t, const Real_wp& dt0,
                                                 Vector<Real_wp>& U0, VirtualOdeSystem<Real_wp>& sys)
  {
    //glob_chrono.SetMessage("Solver_Local", " solve the local ODE");
    //glob_chrono.SetMessage("Source_Local", " compute the source");
    
    dt = dt0;
    
    fine_ode = new LocalOdeSchemeHDG(ci, dt, sys, dof_fine_among_extended,
				     x_glob, y_glob, L_glob, F_glob);
    
    // we create num_explicit
    {
      VarInstationary_Base& var_time = dynamic_cast<VarInstationary_Base&>(sys);
      DistributedProblem_Base& var = var_time.GetDistributedProblem();
      VarComputationProblem& var_p = dynamic_cast<VarComputationProblem&>(var);
      int nb_elt = var_p.GetNbElt();

      Vector<bool> ImplicitElement(nb_elt);
      ImplicitElement.Fill(false);
      for (int i = 0; i < this->num_implicit.GetM(); i++)
	ImplicitElement(this->num_implicit(i)) = true;

      this->num_explicit.Reallocate(nb_elt - this->num_implicit.GetM());
      int nb = 0;
      for (int i = 0; i < nb_elt; i++)
	if (!ImplicitElement(i))
	  {
	    this->num_explicit(nb) = i;
	    nb++;
	  }
    }

    //DISP(this->num_explicit);
    //DISP(this->num_implicit);
    
    // we initialize the coarse and fine region
    fine_ode->Init(this->num_implicit, sys, num_level_volume, num_level_surface, true);
    coarse_ode.Init(this->num_explicit, sys, num_level_volume, num_level_surface);
    
    int Nimp = fine_ode->GetNbLocalDof();
    int Nfar = coarse_ode.GetNbLocalDof();
    int Nfar_ext = coarse_ode.GetNbDofProche();
    
    {
      // we find dofs in common between the coarse and fine region (extended zones)
      const Vector<int>& dof_far = coarse_ode.GetDofProche();
      const Vector<int>& dof_loc = fine_ode->GetDofProche();
      
      coarse_ode.FindCommonDofs(dof_far, dof_loc, dof_common_fine, dof_common_coarse);

      Vector<bool> DofCommon(Nfar_ext);
      DofCommon.Fill(false);
      for (int i = 0; i < dof_common_coarse.GetM(); i++)
	DofCommon(dof_common_coarse(i)) = true;

      dof_original_coarse.Reallocate(Nfar_ext - dof_common_coarse.GetM());
      int nb = 0;
      for (int i = 0; i < Nfar_ext; i++)
	if (!DofCommon(i))
	  {
	    dof_original_coarse(nb) = i;
	    nb++;
	  }
      
      // and conversion between implicit dofs and extended dofs
      Vector<int> all_num;
      const Vector<int>& dof_imp = fine_ode->GetLocalDof();
      coarse_ode.FindCommonDofs(dof_loc, dof_imp, all_num, dof_fine_among_extended);
    }
    
    // we copy the initial condition
    Yn = U0;
    U0.Clear();
    
    // factorization of matrix (M + dt gamma K)
    Real_wp coef = dt*Aimp(1, 1);
    fine_ode->FactorizeOperatorDhPlusGammaKh(Real_wp(1), coef, coef);
        
    int Nproche = fine_ode->GetNbDofProche();
    int nb_dof_uv = Yn.GetM();
    
    // vectors needed for the time scheme
    ui.Reallocate(nb_dof_uv); ui.Zero();
    Fi.Reallocate(nb_dof_uv); Fi.Zero();
    u_fine.Reallocate(Nproche); u_fine.Zero();
    y_imp.Reallocate(Nimp); y_imp.Zero();
    y_far.Reallocate(Nfar); y_far.Zero();
    AImP_y.Reallocate(ci.GetM());
    AP_y.Reallocate(ci.GetM());
    for (int i = 0; i < ci.GetM(); i++)
      {
        AImP_y(i).Reallocate(Nfar_ext); AImP_y(i).Zero();
        AP_y(i).Reallocate(Nproche); AP_y(i).Zero();
      }
  }
  

  //! U^n+1 is computed from U^n and then replaces U^n
  template<>
  void LocalImexRK<Complex_wp>::Advance(const Real_wp& t, int n, VirtualOdeSystem<Complex_wp>& sys)
  {
    cout << "Not implemented for complex numbers" << endl;
    abort();
  }


  //! U^n+1 is computed from U^n and then replaces U^n
  template<>
  void LocalImexRK<Real_wp>::Advance(const Real_wp& t, int n, VirtualOdeSystem<Real_wp>& sys)
  {
    int Nproche = fine_ode->GetNbDofProche();
    int Nimp = fine_ode->GetNbLocalDof();
    int Nfar = coarse_ode.GetNbLocalDof();
    int Nfar_ext = coarse_ode.GetNbDofProche();
    
    const Vector<int>& dof_implicit = fine_ode->GetLocalDof();
    const Vector<int>& dof_explicit = coarse_ode.GetLocalDof();
    const Vector<int>& dof_proche = fine_ode->GetDofProche();
    const Vector<int>& dof_explicit_extended = coarse_ode.GetDofProche();
    
    for (int i = 0; i < ci.GetM(); i++)
      {
        // source is evaluated at point c_i    
        bool presence_source = false;
        if (t < sys.GetFinalTimeSource())
          {
            presence_source = true;
            Fi.Zero();
            sys.AddPrimitiveTimeSource(Real_wp(1), t + ci(i)*dt, 0, Fi);
          }
        
        if (i == 0)
          {
            // first stage : explicit, nothing to do, u_i = y_n
            ui = Yn;

            // we replace Yn by M Yn
            sys.ApplyMassMatrix(Real_wp(1), t, ui, Real_wp(0), Yn);
          }
        else
          {
            //chrono.Reset(1);
            //chrono.Start(1);

            // we compute y_n + dt \sum a_ij k_i (with the explicit/implicit part)
            ui = Yn; Real_wp coef(1);
            for (int j = 0; j < i; j++)
              {
                coef = dt*Aexp(i, j); 
                for (int k = 0; k < Nfar_ext; k++)
                  ui(dof_explicit_extended(k)) += coef*AImP_y(j)(k);
                
                coef = dt*Aimp(i, j);
                for (int k = 0; k < Nproche; k++)
                  ui(dof_proche(k)) += coef*AP_y(j)(k);
              }
            
            // for implicit zone we add gamma dt PF_j
            if (presence_source)
              {
                coef = Aimp(i, i)*dt;
                for (int k = 0; k < Nproche; k++)
                  u_fine(k) = ui(dof_proche(k)) + coef*Fi(dof_proche(k));
              }
            else
              for (int k = 0; k < Nproche; k++)
                u_fine(k) = ui(dof_proche(k));
            
            // for explicit zone, we solve with the mass matrix            
            coarse_ode.SolveOperatorM(sys, ui, coarse_ode.GetDofProche(), true);
            
            // for implicit zone, we solve with (M - gamma dt AP)
            fine_ode->SolveOperatorDhPlusGammaKh(t, u_fine, u_fine);
            
            for (int k = 0; k < Nproche; k++)
              ui(dof_proche(k)) = u_fine(k);            
          }
        
	// we compute P u_i and (I-P) u_i
	for (int k = 0; k < Nimp; k++)
	  y_imp(k) = ui(dof_implicit(k));
	
	for (int k = 0; k < Nfar; k++)
	  y_far(k) = ui(dof_explicit(k));
        
        // we compute AP u_i and A(I-P) u_i
        fine_ode->ApplyOperatorAP(y_imp, AP_y(i), sys, false);
        coarse_ode.ApplyOperatorAP(y_far, AImP_y(i), sys, false);
        
	if (presence_source)
          {
            for (int k = 0; k < dof_original_coarse.GetM(); k++)
              {
                int i0 = dof_original_coarse(k);
                AImP_y(i)(i0) += Fi(dof_explicit_extended(i0));
              }
            
            for (int k = 0; k < Nproche; k++)
              AP_y(i)(k) += Fi(dof_proche(k));
          }
      }	
    
    // on calcule le nouveau itere
    for (int i = 0; i < ci.GetM(); i++)
      {
        Real_wp coef = dt*bi(i);
        for (int k = 0; k < Nproche; k++)
          Yn(dof_proche(k)) += coef*AP_y(i)(k);

        for (int k = 0; k < Nfar_ext; k++)
          Yn(dof_explicit_extended(k)) += coef*AImP_y(i)(k);
      }
        
    sys.SolveMassMatrix(Yn);
    
    //int test_input; cout << "Waiting final" << endl; cin >> test_input;
  }


  //! changes the time step
  template<class T>
  void LocalImexRK<T>::ChangeTimeStep(const Real_wp& dt0)
  {
    dt = dt0;
  }
    

  //! returns the main iterate U^n
  template<class T>
  Vector<T>& LocalImexRK<T>::GetIterate()
  {
    return Yn;
  }


  //! returns the main iterate U^n
  template<class T>
  const Vector<T>& LocalImexRK<T>::GetIterate() const
  {
    return Yn;
  }


  //! returns the number of iterates
  template<class T>
  int LocalImexRK<T>::GetNumberOfIterates() const
  {
    return 1;
  }


  //! returns the main iterate U^n
  template<class T>
  Vector<T>& LocalImexRK<T>::GetIterate(int k)
  {
    return Yn;
  }


  //! clears vectors stored
  template<class T>
  void LocalImexRK<T>::Clear()
  {
    Yn.Clear();
  }


  //! clears vectors stored and Uf = U^n
  template<class T>
  void LocalImexRK<T>::ClearFirst(Vector<T>& Uf)
  {
    Uf = Yn;
    Clear();
  }
  
  
  template<class T>
  size_t LocalImexRK<T>::GetMemorySize() const
  {
    size_t taille = 0;
    taille += Yn.GetMemorySize()
      + x_glob.GetMemorySize() + y_glob.GetMemorySize()
      + L_glob.GetMemorySize() + F_glob.GetMemorySize();
    
    return taille;
  }
  

  template<class T>
  void LocalImexRK<T>::GetMemoryUsed(map<string, size_t>& var) const
  {
    fine_ode->GetMemoryUsed(var);
  }

  
  /* template<class T>
  LocalImperialeScheme_Iterator<T>::LocalImperialeScheme_Iterator()
  {
    dt = 0;
    safety_coef = 1.0;
    two_level = false;
    new_dt = 0;
    change_dt = false;
    dt_opt = 0;
    unstable_algo = false;
    order = 2;
    alpha1 = 0;
    alpha2 = 0;
    rk_scheme = false;
  }
  

  template<class T>
  void LocalImperialeScheme_Iterator<T>::SetInputData(const VectString& param)
  {
    safety_coef = to_num<Real_wp>(param(1));
    rk_scheme = false;
    if (param.GetM() > 2)
      {
        if (!param(2).compare("TWO_LEVEL"))
          {
	    // two-level stable algorithm with modified equation
            two_level = true;
            unstable_algo = false;
          }
        else if (!param(2).compare("TWO_LEVEL_UNSTABLE"))
          {
	    // two-level unstable algorithm with modified equation
            two_level = true;
            unstable_algo = true;
          }
	else if (!param(2).compare("MULTILEVEL_RK"))
	  {
	    // multilevel unstable algorithm with RK scheme
	    rk_scheme = true;
	    unstable_algo = true;
	    two_level = false;
	  }
        else          
          {
	    // multilevel unstable algorithm with modified equation
            two_level = false;
            unstable_algo = true;
          }
      }
    
    if (param.GetM() > 3)
      {
	// Time step dt can be tuned manually to ensure stability
        change_dt = true;
        new_dt = to_num<Real_wp>(param(3));
      }
  }
  
  
  template<class T>
  void LocalImperialeScheme_Iterator<T>::
  SetInitialCondition(const T& t, const T& deltat, const Vector<T>& local_dt,
                      Vector<T>& Y0, VirtualOdeSystem<T>& sys)
  {
    // Yn : solution at current time
    // Ynm1 : solution at iteration n-1
    dt = deltat;
    order = sys.GetTimeSchemeOrder();
    dt_opt = deltat;
    Ynm1 = Y0;
    //Y0.Clear();
    Yn = Ynm1;
    Fy = Ynm1;
    Py = Ynm1;
    
    level.Reallocate(local_dt.GetM());
    level.Fill(0);
    Real_wp coef = safety_coef, one(1.0);
    int level_max = 0;
    // computation of a level for each element
    for (int i = 0; i < local_dt.GetM(); i++)
      {
        Real_wp c(0);
	if (rk_scheme)
	  {
	    // CFL = sqrt( (r-1)^2 -1)
	    // searching k such that c Delta t_e >= dt / sqrt( (k-1)^2 -1)
	    c = 1.0 + sqrt(1.0 + square(dt/(coef*local_dt(i))));
	  }
        else if (order == 2)
          {
            // CFL number of T_opt is 2k+1
	    // searching k such that c Delta t_e >=  dt / (2k+1)
	    // => k >=  1/2 dt / (c Delta t_e) - 1
            c = 0.5*(dt/(coef*local_dt(i)) - 1.0);
          }
        else if (order == 4)
          {
            // CFL >= 2 sqrt(r(r+1))
	    // searching k such that c Delta t_e >= dt / (2 sqrt(k (k+1)) )
	    // k(k+1) >= beta  with beta = (1/2 dt / (c dt_e))^2 
            Real_wp beta = square(0.5*dt/(coef*local_dt(i)));
            c = 0.5*(-1.0 + sqrt(1.0 + 4.0*beta));
          }
        else
          {
            cout << "Order not implemented " << endl;
            abort();
          }
        
        int k = toInteger(Real_wp(ceil(c)));
        if (order == 4)
          {
            if (k == 1)
              k--;
          }

	if (rk_scheme)
	  {
	    if (k <= 4)
	      k = 0;
	    else
	      k = (k-3)/2;
	  }
        
        level(i) = k;          
        level_max = max(level_max, k);
      }
        
    // two level
    if (two_level)
      {
        //level_max++;
        for (int i = 0; i < local_dt.GetM(); i++)
          if (level(i) > 0)
            level(i) = level_max;
      }
    //level.Fill(level_max);
    
    // computation of Chebyshev polynomials
    coefT.Reallocate(level_max+1);
    if (rk_scheme)
      {
	UnivariatePolynomial<Real_wp> Tn, Tnm1, Tnm2, Tnp1, X, ESin, ECos, x, pone;
	x.SetOrder(1); x(1) = 1.0;
	pone.SetOrder(0); pone(0) = 1.0;
	int lvl = 1;
	for (int r = 6; r <= 2*level_max+4; r += 2)
	  {
	    Real_wp beta = sqrt(Real_wp((r-1)*(r-1)-1));
	    X = x; X(1) /= beta;
	    Tnm2 = pone;
	    Tnm1 = X;
	    Tn = 2.0*X*Tnm1 - Tnm2;
	    for (int n = 3; n <= r; n++)
	      {
		Tnp1 = 2.0*X*Tn - Tnm1;
		Tnm2 = Tnm1;
		Tnm1 = Tn;
		Tn = Tnp1;
	      }
	    
	    Tn *= 0.5*(r-2)/sqrt(beta*beta+1);
	    Tnm2 *= -0.5*r/sqrt(beta*beta+1);
	    Tnm1 *= beta/sqrt(beta*beta+1);
	    
	    ESin.SetOrder(r-1);
	    ECos.SetOrder(r);
	    for (int k = 1; k <= r-1; k += 2)
	      ESin(k) = Tnm1(k);
	    
	    ECos(r) = Tn(r);
	    for (int k = 0; k <= r-2; k += 2)
	      ECos(k) = Tn(k)+Tnm2(k);
	    
	    if (r%4 == 2)
	      ECos *= -one;
	    else
	      ESin *= -one;
	    
	    DISP(ESin); DISP(ECos);
	    coefT(lvl).Reallocate(r-4);
	    for (int k = 5; k <= r; k++)
	      {
		if (k%2 == 1)
		  coefT(lvl)(k-5) = ESin(k);
		else
		  coefT(lvl)(k-5) = ECos(k);
		
		if (coefT(lvl)(k-5) < 0)
		  coefT(lvl)(k-5) *= -1.0;
	      }
	    
	    DISP(lvl); DISP(coefT(lvl));
	    lvl++;
	  }
	
	coefP.Reallocate(2*level_max);
	for (int r = 0; r < level_max; r++)
	  {
	    coefP(2*r).Reallocate(level_max+1);
	    coefP(2*r+1).Reallocate(level_max+1);
	    coefP(2*r).Fill(0);
	    coefP(2*r+1).Fill(0);
	    for (int k = 0; k <= r; k++)
	      {
		coefP(2*r)(k) = 0.0;
		coefP(2*r+1)(k) = 0.0;
	      }
	    
	    for (int k = r+1; k <= level_max; k++)
	      {
		coef = coefT(k)(2*r);
		for (int m = 0; m < 2*r; m++)
		  coef /= coefP(m)(k);
		
		coefP(2*r)(k) = coef;
		
		coef = coefT(k)(2*r+1);
		for (int m = 0; m <= 2*r; m++)
		  coef /= coefP(m)(k);
		
		coefP(2*r+1)(k) = coef;
	      }
	    
	    DISP(r); DISP(coefP(2*r)); DISP(coefP(2*r+1));
	  }	 
	
	Q0 = Yn;
	Ynm1.Clear();
      }
    else
      {
	for (int r = 1; r <= level_max; r++)
	  {
	    Real_wp xsi = 1.0;
	    if (order == 4)
	      xsi = 0.5*Real_wp(2.0*r+1)/sqrt(Real_wp(r*(r+1)));
	    
	    //UnivariatePolynomial<Real_wp> X, T0, T1, Tn;
	    UnivariatePolynomial<Real_wp> X, T0, T1, Tn;
	    X.SetOrder(1);
	    X(0) = 0;
	    X(1) = xsi/(2.0*r + 1);
	    if (r%2 == 1)
	      X(1) = -X(1);
	    
	    T0.SetOrder(0); T0(0) = 1.0;
	    T1 = X;
	    
	    for (int n = 1; n < 2*r+1; n++)
	      {
		Tn = 2.0*X*T1 - T0;
		T0 = T1;
		T1 = Tn;
	      }
	    
	    coefT(r).Reallocate(r);
	    for (int n = 1; n <= r; n++)
	      {
		coefT(r)(n-1) = Tn(2*n+1)/xsi;
		if (coefT(r)(n-1) < 0)
		  coefT(r)(n-1) = -coefT(r)(n-1);
	      }
	    
	    DISP(r);
	    DISP(coefT(r));
	  }
	
	// coefficients for operator P
	if (unstable_algo)
	  {
	    coefP.Reallocate(level_max+1);
	    for (int r = 1; r <= level_max; r++)
	      {
		coefP(r).Reallocate(level_max+1);
		for (int k = 0; k < r; k++)
		  coefP(r)(k) = 0.0;
		
		for (int k = r; k <= level_max; k++)
		  {
		    coef = sqrt(coefT(k)(r-1));
		    for (int m = 1; m < r; m++)
		      coef /= coefP(m)(k);
		    
		    coefP(r)(k) = coef;
		  }
		
		DISP(r); DISP(coefP(r));
	      }
	    
	    Qn.Reallocate(level_max+1);
	    for (int k = 1; k <= level_max; k++)
	      Qn(k) = Ynm1;
	    
	    alpha1 = 1.0/(6.0*coefP(1)(1));
	    alpha2 = 2.0/3.0;
	  }
	else
	  {
	    Q0 = Ynm1;
	    Q1 = Ynm1;
	    ImPy = Ynm1;
	    W = Ynm1;
	    
	    int r = level_max;
	    Real_wp xsi = 1.0;
	    if (order == 4)
	      xsi = 0.5*Real_wp(2.0*r+1)/sqrt(Real_wp(r*(r+1)));
	    
	    alpha1 = 1.0/(24.0*xsi)*square(2.0*r+1);
	    alpha2 = 1.0/3.0*(2.0*r+1);
	  }
      }
    
    // calcul de Sh
    typename SystemEqn::Nature_Matrix nat_mat;
    nat_mat.SetCoefMass(0.0);
    nat_mat.SetCoefStiffness(0.0);
    nat_mat.SetCoefDamping(1.0);
    sys.var_harmonic.AddMatrixWithBC(Glob_mat_Sh, nat_mat);

    nat_mat.SetCoefMass(0.0);
    nat_mat.SetCoefStiffness(1.0);
    nat_mat.SetCoefDamping(0.0);
    sys.var_harmonic.AddMatrixWithBC(sys.Glob_mat_Kh, nat_mat);

  }
  
  
  // main method
  template<class T>
  void LocalImperialeScheme_Iterator<T>::Advance(const T& t, int n, VirtualOdeSystem<T>& sys)
  {
    sys.GiveIterate(n, t, Yn);
    
    if (rk_scheme)
      {
	// fourth-order contribution common to all elements
	Real_wp one(1.0);
	Copy(Yn, Q0);
	sys.EvaluateDerivativeFunction(t, 0, Q0, Fy);
	Add(dt_opt, Fy, Yn);
	
	sys.EvaluateDerivativeFunction(t, 1, Fy, Q0);
	Add(dt_opt*dt_opt/2, Q0, Yn);
	
	sys.EvaluateDerivativeFunction(t, 2, Q0, Fy);
	Add(pow(dt_opt, 3.0)/6, Fy, Yn);
	
	sys.EvaluateDerivativeFunction(t, 3, Fy, Q0);
	Mlt(pow(dt_opt, 4.0), Q0); Add(one/24, Q0, Yn);
	
	//Copy(Yn, Fy); Q1 = Q0;
	// then additional terms stabilizing this scheme
	for (int k = 0; k < coefP.GetM(); k++)
	  {
	    Copy(Q0, Py);
	    ApplyOperatorP(k, Py, sys); //DISP(k); DISP(Norm2(Py));
	    sys.EvaluateFunction(t, dt_opt, -1, Py, Q0); //DISP(Norm2(Q0));
	    Add(one, Q0, Yn);
	    //Yn.Write("Yn_L"+to_str(k)+".dat");
	  }
	
	//int test_input; cout << "we wait" << endl; cin >> test_input;
	return;
      }
    
    if (unstable_algo)
      {
        // calcul de Fy = A Yn
        sys.EvaluateFunction(t, Yn, Fy);
        //Yn.Write("Yn.dat"); Fy.Write("Fy.dat");
        
        // calcul de 
        // Q1 = dt A P1 Fy
        // Q2 = dt A P2 Q1
        // Q3 = dt A P3 Q2
        // etc
        // Qn = dt A Pn Qn-1
        for (int r = 1; r < coefT.GetM(); r++)
          {
            if ((order == 4) &&(r == 1))
              {
                Copy(Fy, Py);
                Mlt(coefP(1)(1), Py);
                sys.EvaluateFunction(t, dt_opt, -1, Py, Qn(r));
              }
            else
              {
                if (r == 1)
                  Copy(Fy, Py);
                else
                  Copy(Qn(r-1), Py);
                
                this->ApplyOperatorP(r, Py, sys);
                sys.EvaluateFunction(t, dt_opt, r, Py, Qn(r));
              }
            //Qn(r).Write("Qn"+to_str(r)+".dat");
          }
        
        // calcul de 
        // Qn-1 = Qn-1 + dt A Pn Qn
        // Qn-2 = Qn-2 + dt A Pn-1 Qn-1
        // etc
        // Q1 = Q1 + dt A P2 Q2
        // Fy = Fy + dt A P1 Q1
        for (int r = coefT.GetM()-1; r > 0; r--)
          {
            if ((order == 4) && (r == 1))
              {
                sys.EvaluateFunction(t, coefP(1)(1)*dt_opt, -1, Qn(r), Py, false, false);
                sys.AddPrimitiveTimeSource(alpha2, t+0.5*dt, 0, Py);
                sys.AddPrimitiveTimeSource(alpha2, t-0.5*dt, 0, Py);
                sys.AddPrimitiveTimeSource(-2.0*alpha2, t, 0, Py);
                sys.SolveMassMatrix(Py);
                Add(1.0, Py, Fy);
                //Py.Write("Py.dat");
              }
            else if ((order == 4) && (r == 2))
              {
                this->ApplyOperatorP(r, Qn(r), sys);
                sys.EvaluateFunction(t, dt_opt, 1, Qn(r), Py, false, false);
                //sys.AddPrimitiveTimeSource(t, alpha1*dt, 1, Py);
                sys.AddPrimitiveTimeSource(alpha1, t+0.5*dt, 0, Py);
                sys.AddPrimitiveTimeSource(-alpha1, t-0.5*dt, 0, Py);
                sys.level_time_scheme.SetLevel(-1);
                sys.SolveMassMatrix(Py);
                Add(1.0, Py, Qn(r-1));
                //Qn(r-1).Write("Qn"+to_str(r-1)+".dat");
              }
            else
              {
                this->ApplyOperatorP(r, Qn(r), sys);
                sys.EvaluateFunction(t, dt_opt, r-1, Qn(r), Py);
                if (r == 1)
                  Add(1.0, Py, Fy);
                else
                  Add(1.0, Py, Qn(r-1));
                
                //Qn(r-1).Write("QnB"+to_str(r-1)+".dat");
              }
          }
        
        //Fy.Write("Fy.dat");
        Mlt(2.0*dt, Fy);
      }
    else
      {
        // computation of Fy = dt A P Yn and W = dt (A (I-P) Yn + Fn)
        Copy(Yn, Py);
        this->ApplyOperatorP(0, Py, sys);
                
        Copy(Yn, ImPy); Add(-1.0, Py, ImPy);
	
        // we add source to W
        sys.EvaluateFunction(t, dt_opt, -1, ImPy, W, false, true);
	// and Sh U^n-1 to W
	MltAdd(-dt_opt, Glob_mat_Sh, Ynm1, 1.0, W);
	sys.SolveMassMatrix(W);
        
        sys.EvaluateFunction(t, dt_opt, -1, Py, Fy);
        
        // then first step Q0 = Yn, Q1 = dt/(2r+1) A Yn
        int r = coefT.GetM() - 1;
        Real_wp xsi = 1.0;
        if (order == 4)
          xsi = 0.5*Real_wp(2.0*r+1)/sqrt(Real_wp(r*(r+1)));
        
        Real_wp coef = 2.0*xsi/(2.0*r+1);
        Copy(Yn, Q0);
        Copy(W, Q1); Add(1.0, Fy, Q1);
        Mlt(0.5*coef, Q1);
        
        int kmin = 1;
        if (order == 4)
          {
            kmin = 3;
            // Q2 = Q0 + 2 xsi/(2r+1) dt A Q1
            sys.EvaluateFunction(t, coef*dt_opt, -1, Q1, Fy, false, false);
            sys.AddPrimitiveTimeSource(coef*dt_opt*alpha1, t+0.5*dt, 0, Fy);
            sys.AddPrimitiveTimeSource(-coef*dt_opt*alpha1, t-0.5*dt, 0, Fy);
            sys.SolveMassMatrix(Fy);
            Add(1.0, Q0, Fy);
            Copy(Q1, Q0);
            Copy(Fy, Q1);
            
            // Q3 = Q1 + 2 xsi/(2r+1) dt A Q2
            sys.EvaluateFunction(t, coef*dt_opt, -1, Q1, Fy, false, true);
            sys.AddPrimitiveTimeSource(coef*dt_opt*alpha2, t+0.5*dt, 0, Fy);
            sys.AddPrimitiveTimeSource(coef*dt_opt*alpha2, t-0.5*dt, 0, Fy);
            sys.AddPrimitiveTimeSource(-2.0*coef*dt_opt*alpha2, t, 0, Fy);
            sys.SolveMassMatrix(Fy);
            Add(1.0, Q0, Fy);
            Copy(Q1, Q0);
            Copy(Fy, Q1);
          }
        
        // then small leap frog steps
        for (int k = kmin; k <= 2*r; k++)
          {
            // Q_k+1 = Q_k-1 + 2 dt/(2r+1) (AP Q_k + w delta_{k pair} )
            
            // Fy = dt AP Q_k
            Copy(Q1, Py);
            this->ApplyOperatorP(0, Py, sys);
                        
            sys.EvaluateFunction(t, dt_opt, 1, Py, Fy);
            
            if (k%2 == 0)
              Add(1.0, W, Fy);
            
            Mlt(coef, Fy);
            Add(1.0, Q0, Fy);
	    
            Copy(Q1, Q0);
            Copy(Fy, Q1);
          }
	
	if (r == 0)
	  Copy(Q1, Fy);
        
        Mlt(2.0*dt/(xsi*dt_opt), Fy);
      }    

    glob_chrono.Start(VirtualTimer::SCHEME);
    
    Add(1.0, Ynm1, Fy);
    
    Copy(Yn, Ynm1);
    Copy(Fy, Yn);
    
    glob_chrono.Stop(VirtualTimer::SCHEME);
    //int test_input; cout << "we wait " << endl; cin >> test_input;
  }

  
  template<class T>
  void LocalImperialeScheme_Iterator<T>::
  ApplyOperatorP(int k, VectReal_wp& Y, VirtualOdeSystem<T>& var)
  {
    abort();
    
    typedef typename TypeEquation::TypeEquationStationary TypeEquationStationary;
    typedef typename TypeElement::Dimension Dimension;
    
    const Mesh<Dimension>& mesh = var.var_harmonic.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.var_harmonic.mesh_num;
    int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int offset1 = mesh_num.Element(i).GetNumberDof(0);
        int offset2;
        if (i == mesh.GetNbElt()-1)
          offset2 = Nvol;
        else
          offset2 = mesh_num.Element(i+1).GetNumberDof(0);
           
        if (unstable_algo) 
          {
            if (level(i) < k)
              for (int j = offset1; j < offset2; j++)
                for (int m = 0; m < TypeEquationStationary::nb_unknowns; m++)
                  Y(j + m*Nvol) = 0;
            else
	    {
    	    Real_wp coef = coefP(k)(level(i));
	    for (int j = offset1; j < offset2; j++)
	      for (int m = 0; m < TypeEquationStationary::nb_unknowns; m++)
		Y(j + m*Nvol) *= coef;
	    // }
          }
        else
          {
            if (level(i) == 0)
              for (int j = offset1; j < offset2; j++)
                for (int m = 0; m < TypeEquationStationary::nb_unknowns; m++)
                  Y(j + m*Nvol) = 0;
            
          }
      }

  }

  
  template<class T>
  void RunTimeScheme(const T& t0, const T& tf, T& dt, Vector<T>& local_deltat,
		     Vector<T>& Y0, VirtualOdeSystem<T>& sys,
		     LocalImperialeScheme_Iterator<T>& local_time_scheme)
  {
    T t;
    int nb_max_iter;
    if (dt != T(0))
      nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
    else
      return;
    
    // if tf <= t0, we exit 
    if (nb_max_iter < 1)
      return;
    
    local_time_scheme.SetInitialCondition(t0, dt, local_deltat, Y0, sys);
        
    t = t0;
    int lmax = 0;
    for (int i = 0; i < local_time_scheme.level.GetM(); i++)
      lmax = max(lmax, local_time_scheme.level(i));
    
    IVect repartition(lmax+1); repartition.Fill(0);
    for (int i = 0; i < local_time_scheme.level.GetM(); i++)
      repartition(local_time_scheme.level(i))++;
    
    cout << "Repartition =  " << repartition << endl;
    
    if (local_time_scheme.unstable_algo)
      sys.GiveLevelTime(local_deltat, local_time_scheme.level, sys.LEVEL_COARSE_TO_FINE);
    else
      sys.GiveLevelTime(local_deltat, local_time_scheme.level, sys.LEVEL_I_MINUS_P_P);
    
    sys.GiveNumberIterations(dt, nb_max_iter);
    
    //local_time_scheme.Yn.Fill(1);
    //local_time_scheme.ApplyOperatorP(0, local_time_scheme.Yn, sys);
    //local_time_scheme.Yn.Write("PnVec.dat");
    
    //local_time_scheme.Yn.FillRand(); Mlt(1e-14,  local_time_scheme.Yn);
    //local_time_scheme.Advance(0, 0, sys);
    
    //exit(0);
    
    if (sys.print_level > 6)
      {
      }
    
    if (local_time_scheme.change_dt)
      {
        sys.deltat = dt;
        dt = local_time_scheme.new_dt;
        local_time_scheme.dt = dt;
        if (local_time_scheme.order == 4)
          local_time_scheme.dt_opt = dt;
      }
    
    // iterations in time
    // main loop
    for (int n_time = 0; n_time < nb_max_iter; n_time++)
      {
	t = t0 + n_time*dt;
	local_time_scheme.Advance(t, n_time, sys);
      }
    
  }


  template<class T>
  LocalPipernoScheme_Iterator<T>::LocalPipernoScheme_Iterator()
  {
    new_dt = 0;
    safety_coef = 1.0;
    change_dt = false;
    level_max = 0;
  }
  
  
  template<class T>
  void LocalPipernoScheme_Iterator<T>::SetInputData(const VectString& param)
  {
    safety_coef = to_num<Real_wp>(param(1));
    if (param.GetM() > 2)
      {
        // nothing, here for compatibility
      }
    
    if (param.GetM() > 3)
      {
        change_dt = true;
        new_dt = to_num<Real_wp>(param(3));
      } 
  }
    
  template<class T>
  void LocalPipernoScheme_Iterator<T>::
  SetInitialCondition(const T& t, const T& deltat, const Vector<T>& local_dt,
                      Vector<T>& E0, Vector<T>& H0, VirtualOdeSystem<T>& sys)
  {
    En = E0;
    Hn = H0;
    E0.Clear();
    H0.Clear();
    
    ProdEn = En;
    ProdHn = Hn;
    if (!change_dt)
      new_dt = deltat;

    level.Reallocate(local_dt.GetM());
    level.Fill(0);
    Real_wp coef = safety_coef;
    level_max = 0;
    // computation of a level for each element
    for (int i = 0; i < local_dt.GetM(); i++)
      {
        Real_wp c(0);
        c = log2(deltat/(coef*local_dt(i)));
        int k = toInteger(Real_wp(ceil(c)));
        k = max(k, 0);
        
        level(i) = k;          
        level_max = max(level_max, k);
      }
    
    // same maximum level for all processors
#ifdef SELDON_WITH_MPI
    int lvl = level_max;
    MPI_Allreduce(&lvl, &level_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
#endif
    
    NbEltPerLevel.Reallocate(level_max+1);
    NbEltPerLevel.Fill(0);
    for (int i = 0; i < local_dt.GetM(); i++)
      NbEltPerLevel(level(i))++;
    
  }
   
  
  template<class T>
  void LocalPipernoScheme_Iterator<T>::Advance(const T& t, int n, VirtualOdeSystem<T>& sys)
  {
    sys.GiveIterate(n, t, En);
    
    AdvanceLevel(t, new_dt, 0, sys);
    
  }
  
  
  template<class T>
  void LocalPipernoScheme_Iterator<T>
  ::AdvanceLevel(const T& tn, const Real_wp& dt, int level, VirtualOdeSystem<T>& sys)
  {
    // advancing fine level of delta t/2
    if (level_max > level)
      AdvanceLevel(tn, 0.5*dt, level+1, sys);
    
    // advancing coarse level of delta t with Verlet scheme
    if (NbEltPerLevel(level) > 0)
      {
        sys.SetLevel(level);
        sys.ApplyOperatorRhVectorial(1.0, tn+0.25*dt, En, 0, ProdHn, false);
        sys.AddVectorialTimeSource(1.0, tn+0.25*dt, 0, ProdHn);
        sys.SolveOperatorBh(ProdHn);

        sys.AddVector(0.5*dt, ProdHn, Hn);
      }
    
    sys.ComputeExtrapolationVectorial(En, Hn);
    
    if (NbEltPerLevel(level) > 0)
      {
        sys.ApplyOperatorRhScalar(1.0, tn+0.5*dt, Hn, 0, ProdEn, false);
        sys.AddScalarTimeSource(1.0, tn+0.5*dt, 0, ProdEn);
        sys.SolveOperatorDh(ProdEn);
        
        sys.AddVector(dt, ProdEn, En);
      }
    
    sys.ComputeExtrapolationScalar(En, Hn);
            
    if (NbEltPerLevel(level) > 0)
      {
        sys.ApplyOperatorRhVectorial(1.0, tn+0.75*dt, En, 0, ProdHn, false);
        sys.AddVectorialTimeSource(1.0, tn+0.75*dt, 0, ProdHn);
        sys.SolveOperatorBh(ProdHn);
        
        sys.AddVector(0.5*dt, ProdHn, Hn);
      }
    
    sys.ComputeExtrapolationVectorial(En, Hn);
    
    // advancing fine level of delta t/2
    if (level_max > level)
      AdvanceLevel(tn, 0.5*dt, level+1, sys);
  }


  template<class T>
  void RunTimeScheme(const T& t0, const T& tf, T& dt, Vector<T>& local_deltat,
		     Vector<T>& E0, Vector<T>& H0, VirtualOdeSystem<T>& sys,
		     LocalPipernoScheme_Iterator<T>& local_time_scheme)
  {
    T t;
    int nb_max_iter;
    if (dt != T(0))
      nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
    else
      return;
    
    // if tf <= t0, we exit 
    if (nb_max_iter < 1)
      return;
    
    local_time_scheme.SetInitialCondition(t0, dt, local_deltat, E0, H0, sys);
        
    t = t0;
    int lmax = 0;
    for (int i = 0; i < local_time_scheme.level.GetM(); i++)
      lmax = max(lmax, local_time_scheme.level(i));
    
    IVect repartition(lmax+1); repartition.Fill(0);
    for (int i = 0; i < local_time_scheme.level.GetM(); i++)
      repartition(local_time_scheme.level(i))++;
    
    cout << "Repartition =  " << repartition << endl;
    
    sys.GiveLevelTime(local_deltat, local_time_scheme.level, sys.LEVEL_DISJOINT);
    
    sys.GiveNumberIterations(dt, nb_max_iter);
    
    if (local_time_scheme.change_dt)
      {
        dt = local_time_scheme.new_dt;
        sys.deltat = dt;
        nb_max_iter = toInteger(ceil(abs(tf - t0)/dt ) );
      }
    
    // iterations in time
    // main loop
    for (int n_time = 0; n_time < nb_max_iter; n_time++)
      {
	t = t0 + n_time*dt;
	local_time_scheme.Advance(t, n_time, sys);
      }
    
  }
*/   
}

#define MONTJOIE_FILE_LOCAL_TIME_SCHEMES_CXX
#endif
