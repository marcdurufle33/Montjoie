#ifndef SELDON_FILE_DISTRIBUTED_BLOCK_DIAGONAL_MATRIX_CXX

namespace Seldon
{
  
  //! clears arrays used to perform a matrix vector product
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::EraseArrayForMltAdd()
  {
    local_number_distant_values = false;
    size_max_distant_col = 0;
    global_col_to_recv.Clear(); ptr_global_col_to_recv.Clear();
    local_col_to_send.Clear();
    proc_col_to_recv.Clear(); proc_col_to_send.Clear();
  }
  
    
  //! row numbers will contain global numbers instead of local numbers
  /*!
    Arrays allowing matrix-vector product are erased
    this method is usually called when interactions are added to the matrix
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::SwitchToGlobalNumbers()
  {    
    if (!local_number_distant_values)
      return;
    
    MPI_Comm& comm = comm_;
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    IVect row_num;
    row_num.SetData(this->GetNbRowNumbers(), this->GetRowNumData());
    for (int i = 0; i < row_num.GetM(); i++)
      if (proc_row(i) != rank)
        row_num(i) = global_col_to_recv(row_num(i));
    
    row_num.Nullify();
    
    EraseArrayForMltAdd();
  }
  

  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::CompleteWithAssociatedNumbers(IVect& num)
  {
    MPI_Comm& comm = this->GetCommunicator();
    
    // duplicate numbers are removed if present
    RemoveDuplicate(num);

    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    
    const IVect& global_row = this->GetGlobalRowNumber();
    //const Vector<IVect>& shared_rows = this->GetSharingRowNumbers();
    const IVect& shared_proc = this->GetProcessorSharingRows();
    
    // which processor are used
    Vector<bool> ProcUsed(nb_proc);
    ProcUsed.Fill(false);
    for (int i = 0; i < this->proc_row.GetM(); i++)
      if (this->proc_row(i) != rank)
        ProcUsed(this->proc_row(i)) = true;

    for (int i = 0; i < shared_proc.GetM(); i++)
      ProcUsed(shared_proc(i)) = true;
    
    // we send to all processors in contact global row numbers to clear

    IVect num_glob(num); int nb_row_to_send = num.GetM();
    for (int i = 0; i < num.GetM(); i++)
      num_glob(i) = global_row(num(i));

    // we send and receive data
    Vector<MPI_Request> request_send(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
	MPI_Isend(&nb_row_to_send, 1, MPI_INTEGER, i, 143, comm, &request_send(i));
    
    Vector<MPI_Request> request_recv(nb_proc);
    Vector<int> nb_row_to_recv(nb_proc);
    nb_row_to_recv.Zero();
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
        MPI_Irecv(&nb_row_to_recv(i), 1, MPI_INTEGER, i, 143, comm, &request_recv(i));

    MPI_Status status;
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
	{
	  MPI_Wait(&request_send(i), &status);
	  MPI_Wait(&request_recv(i), &status);
	}

    Vector<Vector<int> > num_to_recv(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      if ((ProcUsed(i)) && (nb_row_to_send > 0))
        MPI_Isend(num_glob.GetData(), nb_row_to_send, MPI_INTEGER, i, 144, comm, &request_send(i));
    
    for (int i = 0; i < nb_proc; i++)
      if (nb_row_to_recv(i) > 0)
	{
	  num_to_recv(i).Reallocate(nb_row_to_recv(i));
	  MPI_Irecv(num_to_recv(i).GetData(), nb_row_to_recv(i), MPI_INTEGER, i, 144, comm, &request_recv(i));
	}
    
    for (int i = 0; i < nb_proc; i++)
      {
	if ((ProcUsed(i)) && (nb_row_to_send > 0))
	  MPI_Wait(&request_send(i), &status);
	
	if (nb_row_to_recv(i) > 0)
	  MPI_Wait(&request_recv(i), &status);
      }

    // compressing array
    for (int i = 0; i < nb_proc; i++)
      if (nb_row_to_recv(i))
	{
	  num_glob.PushBack(num_to_recv(i));
	  num_to_recv(i).Clear();
	}

    RemoveDuplicate(num_glob);

    // generating global numbers for ind array
    IVect global_ind(this->size_ind_), permut(this->size_ind_);
    for (int i = 0; i < this->size_ind_; i++)
      {
	if (this->proc_row(i) == rank)
	  global_ind(i) = global_row(this->ind_[i]);
	else
	  {
	    int iglob = this->ind_[i];
	    if (this->local_number_distant_values)
	      iglob = this->global_col_to_recv(iglob);
	    
	    global_ind(i) = iglob;
	  }
      }

    permut.Fill();
    Sort(global_ind, permut);
    
    // now generating numbers to clear
    int k = 0;
    num.Clear();
    for (int i = 0; i < num_glob.GetM(); i++)
      {
	while ((k < global_ind.GetM()) && (global_ind(k) < num_glob(i)))
	  k++;

	if ((k < global_ind.GetM()) && (global_ind(k) == num_glob(i)))
	  {
	    // we found a row to clear
	    num.PushBack(permut(k));
	  }
      }

    Sort(num);

        // Alternative method to be improved (loop ?)
    // finding all the blocks
    /*    IVect num_block(num.GetM());
    for (int i = 0; i < num.GetM(); i++)
      {
	int nblock = this->GetBlockNumber(num(i));
	num_block(i) = nblock;
      }

    Sort(num_block, num);

    // loop on all the blocks
    int k = 0;
    Vector<Vector<int> > num_to_send(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
	{
	  num_to_send(i).Reallocate(2*this->size_ind_);
	  num_to_send(i).Fill(-1);
	}
    
    IVect offset(nb_proc);
    offset.Zero();
    while (k < num.GetM())
      {
	int k0 = k, nblock = num_block(k);
	while( (k < num.GetM()) && (num_block(k) == nblock))
	  k++;

	DISP(nblock);
	for (int k1 = k0; k1 < k; k1++)
	  DISP(num(k1));
	
	for (int j = row_ptr[nblock]; j < row_ptr[nblock+1]; j++)	  
	  {
	    int p = this->proc_row(j);
	    int jloc = this->ind_[j];
	    int jglob = jloc;
	    if (p != rank)
	      {
		if (this->local_number_distant_values)
		  jglob = this->global_col_to_recv(jloc);
		
		Check_And_ReallocateIfNecessary(num_to_send(p), offset(p)+k-k0+2);
		num_to_send(p)(offset(p)) = jglob;
		num_to_send(p)(offset(p)+1) = k-k0;
		for (int k1 = k0; k1 < k; k1++)
		  num_to_send(p)(offset(p)+k1-k0+2) = global_row(num(k1));
	      
		offset(p) += k-k0 + 2;
	      }
	    else
	      for (int p2 = 0; p2 < shared_proc.GetM(); p2++)
		{
		  // bisection to find jloc
		  int i0 = 0, i1 = shared_rows(p2).GetM()-1;
		  while (i1 > i0+1)
		    {
		      int i2 = (i0+i1)/2;
		      if (shared_rows(p2)(i2) < jloc)
			i0 = i2;
		      else 
			{
			  i1 = i2;
			  if (shared_rows(p2)(i2) == jloc)
			    break;
			}
		    }

		  if ((shared_rows(p2)(i0) == jloc) || (shared_rows(p2)(i1) == jloc))
		    {
		      p = shared_proc(p2);
		      jglob = global_row(jloc);
		      
		      Check_And_ReallocateIfNecessary(num_to_send(p), offset(p)+k-k0+2);
		      num_to_send(p)(offset(p)) = jglob;
		      num_to_send(p)(offset(p)+1) = k-k0;
		      for (int k1 = k0; k1 < k; k1++)
			num_to_send(p)(offset(p)+k1-k0+2) = global_row(num(k1));
		      
		      offset(p) += k-k0 + 2;		      
		    }		  
		}
	  }
      }
    
    DISP(offset); DISP(num_to_send);
    // we send data
    Vector<MPI_Request> request_send(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
	{
	  num_to_send(i).Resize(offset(i));	 	  
	  request_send(i) = comm.Isend(&offset(i), 1, MPI_INTEGER, i, 143);
	}
    
    // we receive the data
    Vector<MPI_Request> request_recv(nb_proc);
    Vector<int> offset_recv(nb_proc);
    offset_recv.Zero();
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
	request_recv(i) = comm.Irecv(&offset_recv(i), 1, MPI_INTEGER, i, 143);

    MPI_Status status;
    for (int i = 0; i < nb_proc; i++)
      if (ProcUsed(i))
	{
        MPI_Wait(&request_send(i), &status);
	  MPI_Wait(&request_recv(i), &status);
	}

    DISP(offset_recv);
    Vector<Vector<int> > num_to_recv(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      if (offset(i) > 0)
	request_send(i) = comm.Isend(num_to_send(i).GetData(), offset(i), MPI_INTEGER, i, 144);

    for (int i = 0; i < nb_proc; i++)
      if (offset_recv(i) > 0)
	{
	  num_to_recv(i).Reallocate(offset_recv(i));
	  request_recv(i) = comm.Irecv(num_to_recv(i).GetData(), offset_recv(i), MPI_INTEGER, i, 144);
	}

    for (int i = 0; i < nb_proc; i++)
      {
	if (offset(i) > 0)
	  request_send(i).Wait(status);
	
	if (offset_recv(i) > 0)
	  request_recv(i).Wait(status);
      }

    DISP(num_to_recv);
    // now we form all the numbers (local and distant) that need to be treated
    for (int i = 0; i < num.GetM(); i++)
      {
	int nblock = this->GetBlockNumber(num(i));
	num(i) = row_ptr[nblock] + this->rev_[num(i)];
      }

    DISP(num);
    IVect num_interac, proc_interac, pos_interac;
    for (int p = 0; p < nb_proc; p++)
      if (offset_recv(p) > 0)
	{
	  k = 0;
	  while (k < offset_recv(p))
	    {
	      num_interac.PushBack(num_to_recv(p)(k));
	      proc_interac.PushBack(p);
	      pos_interac.PushBack(k);
	      k += num_to_recv(p)(k+1) + 2;
	    }
	}

    Sort(num_interac, proc_interac, pos_interac);
    k = 0;
    DISP(num_interac); DISP(proc_interac); DISP(pos_interac);
    for (int i = 0; i < num_interac.GetM(); i++)
      {
	while ((k < global_row.GetM()) && (global_row(k) < num_interac(i)))
	  k++;

	if (global_row(k) != num_interac(i))
	  {
	    cout << "impossible case"  << endl;
	    abort();
	  }
	else
	  {
	    int nblock = this->GetBlockNumber(k);
	    int size_block = this->GetBlockSize(nblock);
	    int p = proc_interac(i);
	    int pos = pos_interac(i);
	    int nb = num_to_recv(p)(pos+1);
	    for (int j = 0; j < nb; j++)
	      {
		int jglob = num_to_recv(p)(pos+2+j);
		for (int j2 = 0; j2 < size_block; j2++)
		  if (this->proc_row(row_ptr[nblock] + j2) == p)
		    {
		      int jglob2 = this->Index(nblock, j2);
		      if (this->local_number_distant_values)
			jglob2 = this->global_col_to_recv(jglob2);
		      
		      if (jglob2 == jglob)
			num.PushBack(row_ptr[nblock]+j2);
		    }
	      }
	  }
      }

    RemoveDuplicate(num);
    */
  }
  

  //! changes global numbers in proc_row to local numbers
  /*!
    \param[in] proc_num processor rank for each row number stored in dist_num
    \param[in] dist_num list of row numbers to analyse
    \param[out] glob_num list of distant row numbers, sorted by processor rank
    it will contain row numbers for processor proc_glob(0), then processor proc_glob(1), etc
    \param[out] ptr_glob_num beginning of index for each processor in array glob_num
    the row numbers of processor proc_glob(i) in array glob_num are stored in
    glob_num(ptr_glob_num(i):ptr_glob_num(i+1))
    \param[out] proc_glob list of distant processors involved in array dist_num
    \param[out] local_num for each processor proc_local(i) local_num(i) will store
    the local row numbers that should be sent to processor proc_local(i) (since these local
    rows will be distant in distant processor proc_local(i)), somehow glob_num stores
    the distants rows that will be received, whereas local_num stores the local rows
    that will be sent
    \param[out] proc_local list of distant processors that will receive local row numbers
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::SortAndAssembleDistantInteractions(const IVect& proc_num, IVect& dist_num,
                                       IVect& glob_num, IVect& ptr_glob_num, IVect& proc_glob,
                                       Vector<IVect>& local_num, IVect& proc_local)
  {
    MPI_Comm& comm = comm_;
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    // counting the number of distant rows
    int nb_glob = 0;
    IVect nb_num_per_proc(nb_proc);
    nb_num_per_proc.Fill(0);
    for (int i = 0; i < proc_num.GetM(); i++)
      if (proc_num(i) != rank)
        {
          nb_num_per_proc(proc_num(i))++;
          nb_glob++;
        }
    
    // sorting numbers
    IVect permut(nb_glob); permut.Fill();
    glob_num.Reallocate(nb_glob);
    IVect all_proc(nb_glob);
    nb_glob = 0;
    for (int i = 0; i < proc_num.GetM(); i++)
      if (proc_num(i) != rank)
        {
          all_proc(nb_glob) = proc_num(i);
          glob_num(nb_glob) = dist_num(i);
          nb_glob++;
        }

    Sort(all_proc, glob_num, permut);
    
    // number of processors involved ?
    int nb_global_proc = 0;
    for (int i = 0 ; i < nb_num_per_proc.GetM(); i++)
      if (nb_num_per_proc(i) > 0)
        nb_global_proc++;
    
    proc_glob.Reallocate(nb_global_proc);
            
    // for each processor, sorting numbers
    all_proc.Clear();
    ptr_glob_num.Reallocate(nb_global_proc+1);
    nb_global_proc = 0; ptr_glob_num(0) = 0;
    for (int i = 0; i < nb_num_per_proc.GetM(); i++)
      if (nb_num_per_proc(i) > 0)
        {
          proc_glob(nb_global_proc) = i;
          ptr_glob_num(nb_global_proc+1) = ptr_glob_num(nb_global_proc) + nb_num_per_proc(i);
          Sort(ptr_glob_num(nb_global_proc),
               ptr_glob_num(nb_global_proc)+ nb_num_per_proc(i) -1, glob_num, permut);
          
          nb_global_proc++;
        }
    
    // changing numbers in proc_row
    IVect inv_permut(nb_glob);
    for (int i = 0; i < nb_glob; i++)
      inv_permut(permut(i)) = i;
    
    nb_glob = 0;
    for (int i = 0; i < proc_num.GetM(); i++)
      if (proc_num(i) != rank)
        {
          dist_num(i) = inv_permut(nb_glob);
          nb_glob++;
        }
    
    // exchanging nb_num_per_proc
    IVect nb_num_send(nb_proc);
    MPI_Alltoall(nb_num_per_proc.GetData(), 1, MPI_INTEGER,
                 nb_num_send.GetData(), 1, MPI_INTEGER, comm);
    
    // sending numbers
    Vector<MPI_Request> request_send(nb_proc), request_recv(nb_proc);
    nb_global_proc = 0;
    for (int p = 0; p < nb_proc; p++)
      if (nb_num_per_proc(p) > 0)
        {
          int size = nb_num_per_proc(p);
          MPI_Isend(&glob_num(ptr_glob_num(nb_global_proc)), size,
                    MPI_INTEGER, p, 17, comm, &request_send(p));
          
          nb_global_proc++;
        }
        
    int nb_local_proc = 0;
    for (int p = 0; p < nb_proc; p++)
      if (nb_num_send(p) > 0)
        nb_local_proc++;
    
    local_num.Reallocate(nb_local_proc);
    proc_local.Reallocate(nb_local_proc);
    
    // receiving numbers
    MPI_Status status; nb_local_proc = 0;
    for (int p = 0; p < nb_proc; p++)
      if (nb_num_send(p) > 0)
        {
          proc_local(nb_local_proc) = p;
          local_num(nb_local_proc).Reallocate(nb_num_send(p));
          MPI_Recv(local_num(nb_local_proc).GetData(), nb_num_send(p),
                   MPI_INTEGER, p, 17, comm, &status);
          
          nb_local_proc++;
        }

    for (int i = 0; i < request_send.GetM(); i++)
      if (nb_num_per_proc(i) > 0)
        MPI_Wait(&request_send(i), &status);
    
    // global to local conversion
    IVect Glob_to_local(this->GetGlobalM());
    const IVect& RowNumber = this->GetGlobalRowNumber();
    Glob_to_local.Fill(-1);
    for (int i = 0; i < RowNumber.GetM(); i++)
      Glob_to_local(RowNumber(i)) = i;
    
    // replacing global numbers with local numbers
    for (int i = 0; i < local_num.GetM(); i++)
      for (int j = 0; j < local_num(i).GetM(); j++)
        {
          local_num(i)(j) = Glob_to_local(local_num(i)(j));
          if (local_num(i)(j) < 0)
            {
              cout << "Impossible case" << endl;
              abort();
            }
        }
    
  }
  

  //! internal function
  template<class T, class Prop, class Storage, class Allocator> template<class T2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::ScatterValues(const Vector<T2>& X, const IVect& num_recv,
                  const IVect& ptr_num_recv, const IVect& proc_recv,
                  const Vector<IVect>& num_send,
                  const IVect& proc_send, Vector<T2>& Xcol) const
  {
    // sending datas
    const MPI_Comm& comm = comm_;
    Vector<Vector<T2> > xsend(proc_send.GetM()), xrecv(proc_recv.GetM());    
    Vector<Vector<int64_t> > xsend_tmp(proc_send.GetM()), xrecv_tmp(proc_recv.GetM());    
    int tag = 30;
    Vector<MPI_Request> request_send(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xsend(i).Reallocate(nb);
        for (int j = 0; j < nb; j++)
          xsend(i)(j) = X(num_send(i)(j));
        
        request_send(i) =
          MpiIsend(comm, xsend(i), xsend_tmp(i), nb, proc_send(i), tag);
      }
    
    // receiving datas
    Vector<MPI_Request> request_recv(proc_recv.GetM());
    int N = 0;
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i); N += nb;
        xrecv(i).Reallocate(nb);
        request_recv(i) =
          MpiIrecv(comm, xrecv(i), xrecv_tmp(i), nb, proc_recv(i), tag);
      }
    
    // waiting that transfers are effective
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      MPI_Wait(&request_recv(i), &status);

    xsend.Clear();
    // completing receives
    for (int i = 0; i < request_recv.GetM(); i++)
      MpiCompleteIrecv(xrecv(i), xrecv_tmp(i), xrecv(i).GetM());
    
    // values are stored in Xcol
    Xcol.Reallocate(N); N = 0;
    for (int i = 0; i < proc_recv.GetM(); i++)
      for (int j = 0; j < ptr_num_recv(i+1) - ptr_num_recv(i); j++)
        Xcol(N++) = xrecv(i)(j);
  }
  

  //! internal function  
  template<class T, class Prop, class Storage, class Allocator> template<class T2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::AssembleValues(const Vector<T2>& Xcol, const IVect& num_recv,
                   const IVect& ptr_num_recv, const IVect& proc_recv,
                   const Vector<IVect>& num_send,
                   const IVect& proc_send, Vector<T2>& X) const
  {
    // sending datas
    const MPI_Comm& comm = comm_;
    Vector<Vector<T2> > xsend(proc_recv.GetM()), xrecv(proc_send.GetM());    
    Vector<Vector<int64_t> > xsend_tmp(proc_recv.GetM()), xrecv_tmp(proc_send.GetM());
    int tag = 32, N = 0;
    Vector<MPI_Request> request_send(proc_recv.GetM());
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i);
        xsend(i).Reallocate(nb);
        for (int j = 0; j < nb; j++)
          xsend(i)(j) = Xcol(N++);
        
        request_send(i) = 
          MpiIsend(comm, xsend(i), xsend_tmp(i), nb, proc_recv(i), tag);
      }
    
    // receiving datas
    Vector<MPI_Request> request_recv(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xrecv(i).Reallocate(nb);
        request_recv(i) = 
          MpiIrecv(comm, xrecv(i), xrecv_tmp(i), nb, proc_send(i), tag);
      }
    
    // waiting that transfers are effective
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      MPI_Wait(&request_recv(i), &status);
    
    xsend.Clear();
    // completing receives
    for (int i = 0; i < request_recv.GetM(); i++)
      MpiCompleteIrecv(xrecv(i), xrecv_tmp(i), xrecv(i).GetM());
    
    // values are added to X
    for (int i = 0; i < num_send.GetM(); i++)
      for (int j = 0; j < num_send(i).GetM(); j++)
        X(num_send(i)(j)) += xrecv(i)(j);
  }
  

  //! selects minimal value instead of summing (as in AssembleRowValues)
  /*!
    The minimal Yproc is searched, and if equality the minimal Y is selected
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::
  AssembleValuesMin(const IVect& Xcol, const IVect& Xcol_proc,
                    const IVect& num_recv, const IVect& ptr_num_recv, const IVect& proc_recv,
                    const Vector<IVect>& num_send, const IVect& proc_send,
                    IVect& Y, IVect& Yproc) const
  {
    // sending datas
    const MPI_Comm& comm = comm_;
    Vector<Vector<int> > xsend(proc_recv.GetM()), xrecv(proc_send.GetM());    
    int cplx = GetRatioMpiDataType(Xcol);
    int tag = 35, N = 0;
    Vector<MPI_Request> request_send(proc_recv.GetM());
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i);
        xsend(i).Reallocate(2*nb);
        for (int j = 0; j < nb; j++)
          {
            xsend(i)(j) = Xcol(N);
            xsend(i)(nb+j) = Xcol_proc(N);
            N++;
          }
        
        MPI_Isend(xsend(i).GetDataVoid(), 2*cplx*nb,
                  GetMpiDataType(Xcol), proc_recv(i), tag,
                  comm, &request_send(i));
      }
    
    // receiving datas
    Vector<MPI_Request> request_recv(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xrecv(i).Reallocate(2*nb);
        MPI_Irecv(xrecv(i).GetDataVoid(), 2*cplx*nb,
                  GetMpiDataType(Xcol), proc_send(i), tag,
                  comm, &request_recv(i));
      }
    
    // waiting that transfers are effective
    MPI_Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      MPI_Wait(&request_send(i), &status);

    for (int i = 0; i < request_recv.GetM(); i++)
      MPI_Wait(&request_recv(i), &status);
    
    xsend.Clear();
    // values are assembled in X
    for (int i = 0; i < num_send.GetM(); i++)
      for (int j = 0; j < num_send(i).GetM(); j++)
        {
          int nb = num_send(i).GetM();
          int proc = xrecv(i)(nb+j);
          int col = xrecv(i)(j);
          if (proc < Yproc(num_send(i)(j)))
            {
              Yproc(num_send(i)(j)) = proc;
              Y(num_send(i)(j)) = col;
            }
          else if (proc == Yproc(num_send(i)(j)))
            {
              if (col < Y(num_send(i)(j)))
                Y(num_send(i)(j)) = col;
            }
        }
  }


  //! assembles the vector by selecting the minimum instead of summing values
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::AssembleVecMin(Vector<int>& X, Vector<int>& Xproc) const
  {
    AssembleVectorMin(X, Xproc, *ProcSharingRows, *SharingRowNumbers,
                      comm_, nodl_scalar_, nb_unknowns_scal_, 13);
  }
  
        
  //! default constructor
  template<class T, class Prop, class Storage, class Allocator>
  DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::DistributedMatrix_BlockDiag()
    : Matrix<T, Prop, Storage, Allocator>()
  {
    local_number_distant_values = false;
    size_max_distant_col = 0;
    
    OverlapRowNumbers = NULL;
    OverlapProcNumbers = NULL;
    GlobalRowNumbers = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nodl_scalar_ = 0;
    nb_unknowns_scal_ = 1;
    nglob_ = 0;

    // default case : the matrix is not distributed
    comm_ = MPI_COMM_SELF;
    
    Matrix<T, Prop, Storage, Allocator> mat;
    if (IsSymmetricMatrix(mat))
      sym_matrix = true;
    else
      sym_matrix = false;
  }
    
  
  //! clears matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::Clear()
  {
    proc_row.Clear();
    EraseArrayForMltAdd();
  }
    

  //! Initialisation of pointers
  /*!
    This method is mandatory, otherwise pointers are set to NULL
    \param[in] n global number of rows (as if it was on a single processor)
    \param[in] row_num local to global numbering
    \param[in] overlap_num rows already counted in another processor
    \param[in] proc_num original processor for overlapped rows
  */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::
  Init(int n, IVect* row_num, IVect* overlap_num, IVect* proc_num,
       int Nvol, int nb_u, IVect* MatchingProc, Vector<IVect>* MatchingDofNumber,
       const MPI_Comm& comm)
  {
    nglob_ = n;
    GlobalRowNumbers = row_num;
    OverlapRowNumbers = overlap_num;
    OverlapProcNumbers = proc_num;
    
    nodl_scalar_ = Nvol;
    nb_unknowns_scal_ = nb_u;
    ProcSharingRows = MatchingProc;
    SharingRowNumbers = MatchingDofNumber;
    
    comm_ = comm;    
  }

  
  //! inits pointers with those of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::Init(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    OverlapRowNumbers = const_cast<IVect*>(&A.GetOverlapRowNumber());
    OverlapProcNumbers = const_cast<IVect*>(&A.GetOverlapProcNumber());
    GlobalRowNumbers = const_cast<IVect*>(&A.GetGlobalRowNumber());
    ProcSharingRows = const_cast<IVect*>(&A.GetProcessorSharingRows());
    SharingRowNumbers = const_cast<Vector<IVect>* >(&A.GetSharingRowNumbers());
    nodl_scalar_ = A.GetNodlScalar();
    nb_unknowns_scal_ = A.GetNbScalarUnknowns();
    nglob_ = A.GetGlobalM();
    comm_ = A.GetCommunicator();
  }
  

  //! returns local to global numbering
  template<class T, class Prop, class Storage, class Allocator>
  const IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::
  GetGlobalRowNumber() const
  {
    if (this->GlobalRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *GlobalRowNumbers;
  }
  
  
  //! returns local to global numbering
  template<class T, class Prop, class Storage, class Allocator>
  IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetGlobalRowNumber()
  {
    if (this->GlobalRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *GlobalRowNumbers;
  }


  //! returns rows already counted in another processor
  template<class T, class Prop, class Storage, class Allocator>
  const IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::
  GetOverlapRowNumber() const
  {
    if (this->OverlapRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapRowNumbers;
  }
  

  //! returns rows already counted in another processor
  template<class T, class Prop, class Storage, class Allocator>
  IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetOverlapRowNumber()
  {
    if (this->OverlapRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapRowNumbers;
  }

  
  //! returns processor numbers of the original rows
  template<class T, class Prop, class Storage, class Allocator>
  const IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::
  GetOverlapProcNumber() const
  {
    if (this->OverlapProcNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapProcNumbers;
  }
  

  //! returns processor numbers of the original rows
  template<class T, class Prop, class Storage, class Allocator>
  IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetOverlapProcNumber()
  {
    if (this->OverlapProcNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapProcNumbers;
  }


  //! returns processor numbers for each set of shared rows
  template<class T, class Prop, class Storage, class Allocator>
  IVect& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::GetProcessorSharingRows()
  {
    if (this->ProcSharingRows == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *ProcSharingRows;
  }
  
  
  //! returns row numbers for each set of shared rows
  template<class T, class Prop, class Storage, class Allocator>
  Vector<IVect>& DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::GetSharingRowNumbers()
  {
    if (this->SharingRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *SharingRowNumbers;
  }


  //! sets blocks (row numbers and processors)
  template<class T, class Prop, class Storage, class Allocator> template<class Allocator1>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks)
  {    
    MPI_Comm& comm = comm_;
    Vector<IVect> ProcDofs;
    ProcDofs.Reallocate(NumDof_Blocks.GetM());
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    for (int i = 0; i < NumDof_Blocks.GetM(); i++)
      {
        ProcDofs(i).Reallocate(NumDof_Blocks(i).GetM());
        ProcDofs(i).Fill(rank_proc);
      }
    
    SetPattern(NumDof_Blocks, ProcDofs);
  }
  
  
  //! sets blocks (row numbers and processors)
  template<class T, class Prop, class Storage, class Allocator> template<class Allocator1>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks,
               const Vector<IVect>& ProcDofs)
  {
    Clear();
    
    MPI_Comm& comm = comm_;
    
    // calling SetPattern of the sequential matrix
    // distant rows are ignored
    int nblock = NumDof_Blocks.GetM();
    Vector<Vector<bool> > ignore_row(nblock);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    for (int i = 0; i < nblock; i++)
      {
        ignore_row(i).Reallocate(NumDof_Blocks(i).GetM());
        ignore_row(i).Fill(true);
        for (int j = 0; j < ProcDofs(i).GetM(); j++)
          if (ProcDofs(i)(j) ==  rank_proc)
            ignore_row(i)(j) = false;
      }
    
    IVect permut;
    Matrix<T, Prop, Storage, Allocator>::SetPattern(NumDof_Blocks,
                                                    ignore_row, permut);
    
    // storing proc_row
    IVect row_num, row_ptr;
    row_num.SetData(this->GetNbRowNumbers(), this->GetRowNumData());
    row_ptr.SetData(nblock+1, this->GetRowPtr());
    proc_row.Reallocate(row_num.GetM());
    for (int i = 0; i < nblock; i++)
      {
        int n = permut(i);
        for (int j = row_ptr(i); j < row_ptr(i+1); j++)
          proc_row(j) = ProcDofs(n)(j-row_ptr(i));
      }
    
    row_num.Nullify(); row_ptr.Nullify();
  }
  
    
  //! adding a distant interaction to the current matrix
  /*!
    If the distant interaction does not belong to the profile of the matrix,
    it is not added. No exception is raised.
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::AddDistantInteraction(int i, int jglob, int proc, const T& val)
  {
    // operation not implemented for local numbers
    SwitchToGlobalNumbers();
    
    const IVect& global_row = this->GetGlobalRowNumber();
    if (sym_matrix)
      {
        // for symmetric matrices, only adding the upper part
        if (global_row(i) > jglob)
          return;
      }
    
    int* row_num = this->GetRowNumData();
    int* row_ptr = this->GetRowPtr();
    int* rev_num = this->GetReverseNum();
    
    // block number associated with local dof
    int n1 = this->GetBlockNumber(i);
    int j1 = rev_num[i];

    // searching j2
    int j2 = -1;
    int rank; MPI_Comm_rank(comm_, &rank);
    for (int j = row_ptr[n1]; j < row_ptr[n1+1]; j++)
      if ((proc_row(j) != rank) && (row_num[j] == jglob))
        {
          j2 = j - row_ptr[n1];
          break;
        }
    
    // if the interaction is not found (j2 = -1)
    // the interaction is not added
    if (j2 >= 0)
      {
        // interaction in the pattern
        this->Value(n1, j1, j2)  += val;
      }
  }
  
  
  //! adding a distant interaction to the current matrix
  /*!
    If the distant interaction does not belong to the profile of the matrix,
    it is not added. No exception is raised.
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::AddRowDistantInteraction(int iglob, int j, int proc, const T& val)
  {
    SwitchToGlobalNumbers();
    
    const IVect& global_row = this->GetGlobalRowNumber();
    if (sym_matrix)
      {
        // for symmetric matrices, only adding the upper part
        if (global_row(j) < iglob)
          return;
      }

    int* row_num = this->GetRowNumData();
    int* row_ptr = this->GetRowPtr();
    int* rev_num = this->GetReverseNum();
        
    // block number associated with local dof
    int n1 = this->GetBlockNumber(j);
    int j2 = rev_num[j];
    
    // searching j1
    int j1 = -1;
    int rank; MPI_Comm_rank(comm_, &rank);
    for (int k = row_ptr[n1]; k < row_ptr[n1+1]; k++)
      if ((proc_row(k) != rank) && (row_num[k] == iglob))
        {
          j1 = k - row_ptr[n1];
          break;
        }
    
    if (j1 >= 0)
      {
        // interaction in the pattern
        this->Value(n1, j1, j2)  += val;
      }
  }
  
  
  //! returns global number
  template<class T, class Prop, class Storage, class Allocator>
  int DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::IndexGlobal(int i, int j) const
  {
    if (GlobalRowNumbers == NULL)
      return this->Index(i, j);
    
    const MPI_Comm& comm = comm_;
    int irow = this->ind_[this->ptr_[i] + j];
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    // if the row belongs to the processor
    // GlobalRowNumbers is used in order to return the global number
    if (proc_row(this->ptr_[i] + j) == rank_proc)
      return (*GlobalRowNumbers)(irow);
    
    // otherwise, the row does not belong to the processor
    if (local_number_distant_values)
      return global_col_to_recv(irow);
    
    return irow;
  }
  
  
  //! writes the matrix on a file
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::WriteText(const string& file_name) const
  {
    const MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::WriteText(file_name);
    
    // a different file name for each processor
    string name = GetBaseString(file_name) + "_P" + to_str(rank_proc)
      + "." + GetExtension(file_name);
    
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("DistributedMatrix_BlockDiag::WriteText(string FileName)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close(); 
  }
  
  
  //! writes the matrix on a file
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {
    const MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc == 1)
      return Matrix<T, Prop, Storage, Allocator>::WriteText(FileStream);
    
    int nblok = this->GetNbBlocks();
    const IVect& global_row = this->GetGlobalRowNumber();
    IVect row_ptr;
    row_ptr.SetData(nblok+1, this->GetRowPtr());

    for (int i = 0; i < nblok; i++)
      {
        int size = this->GetBlockSize(i);
        IVect row(size);
        for (int j = 0; j < size; j++)
          {
            if (proc_row(row_ptr(i) + j) == rank_proc)
              row(j) = global_row(this->Index(i, j));
            else
              {
                if (local_number_distant_values)
                  row(j) = global_col_to_recv(this->Index(i, j));
                else
                  row(j) = this->Index(i, j);
              }
          }

        for (int j = 0; j < size; j++)
          for (int k = 0; k < size; k++)
            FileStream << row(j) + 1 << ' ' << row(k) + 1 << ' ' << this->Value(i, j, k) << '\n';
      }
    
    row_ptr.Nullify();
  }
  
  
  //! conversion to a distributed sparse matrix
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Prop0, class Storage0>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::ConvertToSparse(DistributedMatrix<T0, Prop0, Storage0>& A) const
  {
    int nb_proc; MPI_Comm_size(comm_, &nb_proc);
    if (nb_proc == 1)
      return Seldon::
        ConvertToSparse(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), A);
    
    A.Clear();
    A.Reallocate(this->m_, this->m_);
    A.Init(nglob_, GlobalRowNumbers, OverlapRowNumbers, OverlapProcNumbers,
           nodl_scalar_, nb_unknowns_scal_, ProcSharingRows, SharingRowNumbers, comm_);
    
    int rank; MPI_Comm_rank(comm_, &rank);
    int size_max = this->GetSizeMaxBlock();
    Vector<int> num(size_max);
    Vector<T> val(size_max);
    int* ptr_ = this->GetRowPtr();
    T zero; SetComplexZero(zero);
    // loop over blocks
    for (int i = 0; i < this->GetNbBlocks(); i++)
      {
        int offset = ptr_[i];
        int size_block = this->GetBlockSize(i);
        int nb = 0;
        // local interactions will be added by calling AddInteractionRow
        // whereas distant interactions are added with 
        // AddDistantInteraction/AddRowDistantInteraction
        
        // local number are stored in num
        for (int j = 0; j < size_block; j++)
          if (proc_row(offset+j) == rank)
            {
              num(nb) = this->Index(i, j);
              nb++;
            }
        
        for (int j = 0; j < size_block; j++)
          {
            int irow = this->Index(i, j);
            if (proc_row(offset+j) != rank)
              {
                if (local_number_distant_values)
                  irow = global_col_to_recv(irow);
              }
            
            nb = 0;
            for (int k = 0; k < size_block; k++)
              {
                if (proc_row(offset+k) == rank)
                  {
                    if (proc_row(offset+j) == rank)
                      {
                        val(nb) = this->Value(i, j, k);
                        nb++;
                      }
                    else
                      A.AddRowDistantInteraction(irow, this->Index(i, k),
                                                 proc_row(offset+j), this->Value(i, j, k));
                  }
                else
                  {
                    int jcol = this->Index(i, k);
                    if (local_number_distant_values)
                      jcol = global_col_to_recv(jcol);

                    if (proc_row(offset+j) == rank)
                      A.AddDistantInteraction(irow, jcol,
                                              proc_row(offset+k), this->Value(i, j, k));
                    else
                      {
                        // case not treated
                        if (this->Value(i, j, k) != zero)
                          {
                            // case where row and column are distant
                            cout << "Case not treated" << endl;
                            abort();
                          }
                      }
                  }
              }
            
            if (nb > 0)
              A.AddInteractionRow(irow, nb, num, val);
          }
      }
  }


  //! returns true if the matrix is ready for a matrix-vector product
  template<class T, class Prop, class Storage, class Allocator>
  bool DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::IsReadyForMltAdd() const
  {
    return local_number_distant_values;
  }
  
  
  //! prepares matrix for matrix vector product
  /*!
    This method prepares the matrix-vector, such that a call to MltAdd
    with a distributed will perform exchanges of values between values in order
    to produce the correct result.
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::PrepareMltAdd()
  {
    if (local_number_distant_values)
      return;
      
    local_number_distant_values = true;
    
    IVect row_num;
    row_num.SetData(this->GetNbRowNumbers(), this->GetRowNumData());
    
    SortAndAssembleDistantInteractions(proc_row, row_num,
                                       global_col_to_recv, ptr_global_col_to_recv,
                                       proc_col_to_recv, local_col_to_send, proc_col_to_send);
    
    row_num.Nullify();
  }

  
  //! exchanges values on distant columns
  /*!
    this method exchanges values of X between processors so that Xcol will contain
    the values of X on distant columns dist_col
   */
  template<class T, class Prop, class Storage, class Allocator> template<class T2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::ScatterColValues(const Vector<T2>& X, Vector<T2>& Xcol) const
  {
    ScatterValues(X, global_col_to_recv, ptr_global_col_to_recv, proc_col_to_recv,
                  local_col_to_send, proc_col_to_send, Xcol);
  }
  
  
  //! sums results of distant rows
  /*!
    this method exchanges values of Xrow between processors, and adds them to X
    so that X will sum the results of distant rows dist_row
   */
  template<class T, class Prop, class Storage, class Allocator> template<class T2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::AssembleRowValues(const Vector<T2>& Xrow, Vector<T2>& X) const
  {
    AssembleValues(Xrow, global_col_to_recv, ptr_global_col_to_recv, proc_col_to_recv,
                   local_col_to_send, proc_col_to_send, X);
  }
  
  
  //! assembles the vector (when there are overlapped rows)
  template<class T, class Prop, class Storage, class Allocator> template<class T2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::AssembleVec(Vector<T2>& X) const
  {
    AssembleVector(X, MPI_SUM, *ProcSharingRows, *SharingRowNumbers,
                   comm_, nodl_scalar_, nb_unknowns_scal_, 14);
  }

  
  //! replaces X by A X
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Allocator2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::Mlt_(Vector<T2, VectFull, Allocator2>& X, bool assemble) const
  {
    const MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank; MPI_Comm_rank(comm, &rank);
    if (nb_proc == 1)
      {
        Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), X);
        return;
      }
        
    if (!IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>& >(*this)
          .PrepareMltAdd();
      }
    
    // scattering columns values
    Vector<T2> Xcol;
    this->ScatterColValues(X, Xcol);
    
    int N = global_col_to_recv.GetM();
    Vector<T2> Yrow(N);
    T2 zero; SetComplexZero(zero);
    Yrow.Fill(zero);
    
    // performing matrix vector product
    int* ptr_ = this->GetPtr();
    int* offset_ = this->GetOffset();
    int* ind_ = this->GetInd();
    T* data_ = this->GetData();
    T2* x_ptr = X.GetData(); T2 val, xj; 
    Vector<T2> x_loc(this->GetSizeMaxBlock());
    for (int i = 0; i < this->GetNbBlocks(); i++)
      {
	int nb = offset_[i];
	if (sym_matrix)
          {
            x_loc.Fill(zero);
            for (int j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                if (proc_row(j) == rank)
                  xj = x_ptr[ind_[j]];
                else
                  xj = Xcol(ind_[j]);
                
                val = data_[nb++]*xj;
                for (int k = j+1; k < ptr_[i+1]; k++)
                  {
                    x_loc(k-ptr_[i]) += data_[nb]*xj;
                    if (proc_row(k) == rank)
                      val += data_[nb++]*x_ptr[ind_[k]];
                    else
                      val += data_[nb++]*Xcol(ind_[k]);
                  }
              
                x_loc(j-ptr_[i]) += val;                
              }
          }
        else
          for (int j = ptr_[i]; j < ptr_[i+1]; j++)
            {
              val = zero;
              for (int k = ptr_[i]; k < ptr_[i+1]; k++)
                {
                  if (proc_row(k) == rank)
                    val += data_[nb++]*x_ptr[ind_[k]];
                  else
                    val += data_[nb++]*Xcol(ind_[k]);
                }
              
              x_loc(j-ptr_[i]) = val;
            }
        
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (proc_row(j) == rank)
              X(ind_[j]) = x_loc(j-ptr_[i]);
            else
              Yrow(ind_[j]) = x_loc(j-ptr_[i]);
          }
      }

    // then assembling
    if (assemble)
      {
        this->AssembleRowValues(Yrow, X);
        this->AssembleVec(X);
      }
  }
  
  
  //! computes Y = Y + alpha A X
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class T1, class Allocator1, class T3, class T2, class Allocator2>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::MltAdd_(const T0& alpha, const SeldonTranspose& trans,
            const Vector<T1, VectFull, Allocator1>& X, const T3& beta,
            Vector<T2, VectFull, Allocator2>& Yres, bool assemble) const
  {
    const MPI_Comm& comm = comm_;
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        MltAdd(alpha, trans,
               static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
               X, beta, Yres);
        
        return;
      }
        
    T3 zero; SetComplexZero(zero);
    Vector<T2, VectFull, Allocator2> Y;
    
    if (beta == zero)
      Y.SetData(Yres.GetM(), Yres.GetData());
    else
      Y.Reallocate(Yres.GetM());
    
    Y.Fill(zero);

    if (!IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>& >(*this)
          .PrepareMltAdd();
      }
    
    // scattering columns values
    Vector<T1> Xcol;
    this->ScatterColValues(X, Xcol);
    
    int N = global_col_to_recv.GetM();
    Vector<T2> Yrow(N);
    Yrow.Fill(zero);
    
    // performing matrix vector product
    int* ptr_ = this->GetPtr();
    int* offset_ = this->GetOffset();
    int* ind_ = this->GetInd();
    T* data_ = this->GetData();
    T1* x_ptr = X.GetData();
    T2* y_ptr = Y.GetData(); T2 val, xj; 
    Vector<T2> x_loc(this->GetSizeMaxBlock());
    for (int i = 0; i < this->GetNbBlocks(); i++)
      {
	int nb = offset_[i];
	if (sym_matrix)
          {
            x_loc.Fill(zero);
            if (trans.ConjTrans())
	      {
		for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		  {
		    if (proc_row(j) == rank)
		      xj = x_ptr[ind_[j]];
		    else
		      xj = Xcol(ind_[j]);
		    
		    val = conjugate(data_[nb++])*xj;
		    for (int k = j+1; k < ptr_[i+1]; k++)
		      {
			x_loc(k-ptr_[i]) += conjugate(data_[nb])*xj;
			if (proc_row(k) == rank)
			  val += conjugate(data_[nb++])*x_ptr[ind_[k]];
			else
			  val += conjugate(data_[nb++])*Xcol(ind_[k]);
		      }
		    
		    x_loc(j-ptr_[i]) += val;
		  }            
	      }
	    else
	      {
		for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		  {
		    if (proc_row(j) == rank)
		      xj = x_ptr[ind_[j]];
		    else
		      xj = Xcol(ind_[j]);
		    
		    val = data_[nb++]*xj;
		    for (int k = j+1; k < ptr_[i+1]; k++)
		      {
			x_loc(k-ptr_[i]) += data_[nb]*xj;
			if (proc_row(k) == rank)
			  val += data_[nb++]*x_ptr[ind_[k]];
			else
			  val += data_[nb++]*Xcol(ind_[k]);
		      }
		    
		    x_loc(j-ptr_[i]) += val;
		  }
	      }
	    
	    for (int j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		if (proc_row(j) == rank)
		  y_ptr[ind_[j]] = x_loc(j-ptr_[i]);
		else
		  Yrow(ind_[j]) = x_loc(j-ptr_[i]);
	      }
          }
        else
          {
	    if (trans.Trans())
	      {
		for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		  {
		    int row = ind_[j];
		    if (proc_row(j) == rank)
		      val = x_ptr[row];
		    else
		      val = Xcol(row);
		    
		    for (int k = ptr_[i]; k < ptr_[i+1]; k++)
		      {
			if (proc_row(k) == rank)
			  y_ptr[ind_[k]] += data_[nb++]*val;
			else
			  Yrow(ind_[k]) += data_[nb++]*val;
		      }
		  }		
	      }
	    else if (trans.ConjTrans())
	      {
		for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		  {
		    int row = ind_[j];
		    if (proc_row(j) == rank)
		      val = x_ptr[row];
		    else
		      val = Xcol(row);
		    
		    for (int k = ptr_[i]; k < ptr_[i+1]; k++)
		      {
			if (proc_row(k) == rank)
			  y_ptr[ind_[k]] += conjugate(data_[nb++])*val;
			else
			  Yrow(ind_[k]) += conjugate(data_[nb++])*val;
		      }
		  }
	      }
	    else
	      {
		for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		  {
		    int row = ind_[j];
		    val = zero;
		    for (int k = ptr_[i]; k < ptr_[i+1]; k++)
		      {
			if (proc_row(k) == rank)
			  val += data_[nb++]*x_ptr[ind_[k]];
			else
			  val += data_[nb++]*Xcol(ind_[k]);
		      }
		    
		    if (proc_row(j) == rank)
		      y_ptr[row] = val;
		    else
		      Yrow(row) = val;
		  }
	      }
          }
      }

    // then assembling
    if (assemble)
      {
        this->AssembleRowValues(Yrow, Y);
        this->AssembleVec(Y);
      }
    
    if (beta == zero)
      {
        Mlt(alpha, Y);
        Y.Nullify();
      }
    else
      {
        Mlt(beta, Yres);
        Add(alpha, Y, Yres);
      }
  }
    
  
  //! Adds alpha A to the current matrix
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>
  ::Add_(const T0& alpha, const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        Add(alpha, A,
            static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this));
        
        return;
      }
    
    this->SwitchToGlobalNumbers();
    
    // treating only the case where the profile of A
    // is included in the profile of the current matrix (*this)
    int* ptrA = A.GetRowPtr();
    int* indA = A.GetRowNumData();
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        // finding the block number in the current matrix (num)
        int num = -1;
        for (int j = ptrA[i]; j < ptrA[i+1]; j++)
          {
            int p = A.proc_row(j);
            if (p == rank)
              num = this->GetBlockNumber(indA[j]);
          }
        
        if (num == -1)
          {
            cout << "Empty block ? " << endl;
            abort();
          }
        
        // then checking that all rows of the block of A
        // are belonging to the block num of the current matrix
        for (int j = ptrA[i]; j < ptrA[i+1]; j++)
          {
            int p = A.proc_row(j);
            int b = -1;
            if (p == rank)
              b = this->GetBlockNumber(indA[j]);
            else
              {
                if (A.local_number_distant_values)
                  {
                    for (int k = this->ptr_[num]; k < this->ptr_[num+1]; k++)
                      if (this->proc_row(k) != rank)
                        if (this->ind_[k] == A.global_col_to_recv(indA[j]))
                          b = num;
                  }
                else
                  {
                    for (int k = this->ptr_[num]; k < this->ptr_[num+1]; k++)
                      if (this->proc_row(k) != rank)
                        if (this->ind_[k] == indA[j])
                          b = num;
                  }
              }
            
            if (b != num)
              {
                cout << "The matrix A is not included in the current matrix" << endl;
                cout << "The matrix A cannot be added (not implemented) " << endl;
                abort();
              }
          }
      }
    
    // adding contents of A with AddInteraction
    int size_max = A.GetSizeMaxBlock();
    Vector<int> num(size_max);
    Vector<T> val(size_max);
    T zero; SetComplexZero(zero);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        int offset = ptrA[i];
        int size_block = A.GetBlockSize(i);
        int nb = 0;
        for (int j = 0; j < size_block; j++)
          if (A.proc_row(offset+j) == rank)
            {
              num(nb) = A.Index(i, j);
              nb++;
            }
        
        for (int j = 0; j < size_block; j++)
          {
            int irow = A.Index(i, j);
            if (A.proc_row(offset+j) != rank)
              {
                if (A.local_number_distant_values)
                  irow = A.global_col_to_recv(irow);
              }
            
            nb = 0;
            for (int k = 0; k < size_block; k++)
              {
                if (A.proc_row(offset+k) == rank)
                  {
                    if (A.proc_row(offset+j) == rank)
                      {
                        val(nb) = alpha*A.Value(i, j, k);
                        nb++;
                      }
                    else
                      this->AddRowDistantInteraction(irow, A.Index(i, k),
                                                     A.proc_row(offset+j),
                                                     alpha*A.Value(i, j, k));
                  }
                else
                  {
                    int jcol = A.Index(i, k);
                    if (A.local_number_distant_values)
                      jcol = A.global_col_to_recv(jcol);

                    if (A.proc_row(offset+j) == rank)
                      this->AddDistantInteraction(irow, jcol,
                                                  A.proc_row(offset+k),
                                                  alpha*A.Value(i, j, k));
                    else
                      {
                        // case not treated
                        if (A.Value(i, j, k) != zero)
                          {
                            cout << "Case not treated" << endl;
                            abort();
                          }
                      }
                  }
              }
            
            if (nb > 0)
              this->AddInteractionRow(irow, nb, num, val);
          }
      }
  }

  
  //! assembles matrix
  /*!
    blocks stored in different processors are added, so that each processor
    contains the assembled block (equal to a sequential execution)
    As a result, the functions Mlt, MltAdd should be called with
    assemble = false in order to obtain the correct result
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix_BlockDiag<T, Prop, Storage, Allocator>::Assemble()
  {
    if (!local_number_distant_values)
      this->PrepareMltAdd();
    
    MPI_Comm& comm = comm_;
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return;
    
    // equivalent of MltMin to retrieve the minimal processor number and block number
    // for each row number
    // it tells to which processor each block should be send
    // the assembling procedure will be performed in this processor
    // and blocks will be sent back
    int m = this->GetM();
    Vector<int> Y(m), Yproc(m);
    Yproc.Fill(rank);
    for (int i = 0; i < m; i++)
      Y(i) = this->GetBlockNumber(i);
    
    int N = this->global_col_to_recv.GetM();
    Vector<int> Yrow(N), Yrow_proc(N);
    Yrow_proc.Fill(rank);
    for (int i = 0; i < this->GetNbBlocks(); i++)
      for (int j = this->ptr_[i]; j < this->ptr_[i+1]; j++)
        if (this->proc_row(j) != rank)
          Yrow(this->ind_[j]) = i;

    this->AssembleValuesMin(Yrow, Yrow_proc, global_col_to_recv,
                            ptr_global_col_to_recv, proc_col_to_recv,
                            local_col_to_send, proc_col_to_send, Y, Yproc);
    
    this->AssembleVecMin(Y, Yproc);
    
    // counting the number of blocks to send to each processor
    IVect proc_block(this->GetNbBlocks()), num_block_dist(this->GetNbBlocks());
    for (int i = 0; i < m; i++)
      {
        proc_block(this->GetBlockNumber(i)) = Yproc(i);
        num_block_dist(this->GetBlockNumber(i)) = Y(i);
      }
    
    Y.Clear(); Yproc.Clear(); Yrow.Clear(); Yrow_proc.Clear();
    
    IVect nb_block_send(nb_proc);
    nb_block_send.Fill(0);
    for (int i = 0; i < this->GetNbBlocks(); i++)
      nb_block_send(proc_block(i))++;
    
    nb_block_send(rank) = 0;
    
    // exchanging the number of blocks
    IVect nb_block_recv(nb_proc);
    MPI_Alltoall(nb_block_send.GetData(), 1, MPI_INTEGER,
                 nb_block_recv.GetData(), 1, MPI_INTEGER, comm);
    
    // sending blocks
    Vector<MPI_Request> request_send(nb_proc), request_send_int(nb_proc);
    Vector<Vector<int> > DataIntSend(nb_proc);
    Vector<Vector<T> > DataSend(nb_proc);
    Vector<Vector<int64_t> > DataSend_tmp(nb_proc);
    for (int p = 0; p < nb_block_send.GetM(); p++)
      if (nb_block_send(p) > 0)
        {
          int N = nb_block_send(p);
          int sizeS = 0;
          DataIntSend(p).Reallocate(3*N); N = 0;
          for (int k = 0; k < proc_block.GetM(); k++)
            if (proc_block(k) == p)
              {
                int block_size = this->GetBlockSize(k);
                int nb = 0;
                if (sym_matrix)
                  nb = block_size*(block_size+1)/2;
                else
                  nb = block_size*block_size;
                
                sizeS += nb;
                DataIntSend(p)(N++) = num_block_dist(k);
                DataIntSend(p)(N++) = k;
                DataIntSend(p)(N++) = nb;
              }
          
          MPI_Isend(DataIntSend(p).GetData(), N,
                    MPI_INTEGER, p, 26, comm, &request_send_int(p));
          
          DataSend(p).Reallocate(sizeS);
          sizeS = 0;
          for (int k = 0; k < proc_block.GetM(); k++)
            if (proc_block(k) == p)
              {
                int block_size = this->GetBlockSize(k);
                int nb = 0;
                if (sym_matrix)
                  nb = block_size*(block_size+1)/2;
                else
                  nb = block_size*block_size;
                
                for (int j = 0; j < nb; j++)
                  DataSend(p)(sizeS++) = this->data_[this->offset_[k]+j];
              }
          
          request_send(p) = 
            MpiIsend(comm, DataSend(p), DataSend_tmp(p), sizeS, p, 27);
        }
    
    // receiving blocks
    Vector<Vector<T> > DataRecv(nb_proc);
    Vector<Vector<int> > DataIntRecv(nb_proc);
    MPI_Status status;
    Vector<Vector<int64_t> > DataRecv_tmp(nb_proc);
    for (int p = 0; p < nb_block_recv.GetM(); p++)
      if (nb_block_recv(p) > 0)
        {
          DataIntRecv(p).Reallocate(3*nb_block_recv(p));
          MPI_Recv(DataIntRecv(p).GetData(), 3*nb_block_recv(p), MPI_INTEGER, p, 26, comm, &status);
          
          int sizeS = 0;
          for (int k = 0; k < nb_block_recv(p); k++)
            sizeS += DataIntRecv(p)(3*k+2);
          
          DataRecv(p).Reallocate(sizeS);
          MpiRecv(comm, DataRecv(p), DataRecv_tmp(p), sizeS, p, 27, status);
        }
    
    for (int i = 0; i < request_send.GetM(); i++)
      if (nb_block_send(i) > 0)
        {
          MPI_Wait(&request_send_int(i), &status);
          MPI_Wait(&request_send(i), &status);
        }
    
    // adding blocks
    for (int p = 0; p < nb_block_recv.GetM(); p++)
      {
        int sizeS = 0;
        for (int k = 0; k < nb_block_recv(p); k++)
          {
            int num = DataIntRecv(p)(3*k);
            int nb = DataIntRecv(p)(3*k+2);
            for (int j = 0; j < nb; j++)
              this->data_[this->offset_[num]+j] += DataRecv(p)(sizeS++);
          }
      }
    
    // sending back blocks to each processor
    for (int p = 0; p < nb_block_recv.GetM(); p++)
      if (nb_block_recv(p) > 0)
        {
          int sizeS = 0;
          for (int k = 0; k < nb_block_recv(p); k++)
            {
              int num = DataIntRecv(p)(3*k);
              int nb = DataIntRecv(p)(3*k+2);
              for (int j = 0; j < nb; j++)
                DataRecv(p)(sizeS++) = this->data_[this->offset_[num]+j];
            }
          
          request_send(p) =
            MpiIsend(comm, DataRecv(p), DataRecv_tmp(p), sizeS, p, 29);
        }
    
    for (int p = 0; p < nb_block_send.GetM(); p++)
      if (nb_block_send(p) > 0)
        MpiRecv(comm, DataSend(p), DataSend_tmp(p), DataSend(p).GetM(), p, 29, status);
    
    for (int i = 0; i < request_send.GetM(); i++)
      if (nb_block_recv(i) > 0)
        MPI_Wait(&request_send(i), &status);
    
    // updating blocks
    for (int p = 0; p < nb_block_send.GetM(); p++)
      if (nb_block_send(p) > 0)
        {
          int sizeS = 0;
          for (int k = 0; k < nb_block_send(p); k++)
            {
              int nb = DataIntSend(p)(3*k+2);
              int num = DataIntSend(p)(3*k+1);
              for (int j = 0; j < nb; j++)
                this->data_[this->offset_[num]+j] = DataSend(p)(sizeS++);
            }
        }
    
  }
  
  
  //! computes B = B + alpha A
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const DistributedMatrix<T1, General, BlockDiagRow, Allocator1>& A,
		 DistributedMatrix<T2, General, BlockDiagRow, Allocator2>& B)
  {
    B.Add_(alpha, A);
  }
  
  
  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, General, BlockDiagRow, Alloc>& A,
                       DistributedMatrix<T, General, ArrayRowSparse, Alloc2>& B)
  {
    A.ConvertToSparse(B);
  }


  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, General, BlockDiagRow, Alloc>& A,
                       DistributedMatrix<T, General, ArrayRowComplexSparse, Alloc2>& B)
  {
    A.ConvertToSparse(B);
  }
  
  
  //! computes B = B + alpha A
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const DistributedMatrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		 DistributedMatrix<T2, Symmetric, BlockDiagRowSym, Allocator2>& B)
  {
    B.Add_(alpha, A);
  }
  
  
  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       DistributedMatrix<T, Symmetric, ArrayRowSymSparse, Alloc2>& B)
  {
    A.ConvertToSparse(B);
  }


  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       DistributedMatrix<T, Symmetric, ArrayRowSymComplexSparse, Alloc2>& B)
  {
    A.ConvertToSparse(B);
  }
  
  
  //! converts a sparse matrix into a block-diagonal one
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
			      DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B,
			      int m1, int m2)
  {
    if (m2 == 0)
      m2 = A.GetM();
    
    const MPI_Comm& comm = A.GetCommunicator();
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return
        ConvertToBlockDiagonal(static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A),
                               static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(B), m1, m2);
    
    int m = A.GetM();
    // in order to find the blocks, we use the fact that by iterating products of
    // A (somehow calculating A^n), the blocks are filled
    // Each block is identified by the minimal processor storing the block
    // and the minimal global row number of the block in this processor
    // therefore Y will contain the minimal global row numbers and Yproc
    // the minimal processor
    Vector<int> Y(m), Yproc(m), Yprev(m), Yproc_prev(m);
    const IVect& global = A.GetGlobalRowNumber();

    // starting from a diagonal profile (blocks of size 1)
    Yproc.Fill(rank);
    for (int i = 0; i < m; i++)
      Y(i) = global(i);
    
    bool test_loop = true;
    while (test_loop)
      {
        Copy(Y, Yprev);
        Copy(Yproc, Yproc_prev);
        
        // iterating a matrix-vector product where we overwrite the minimal column
        // number and minimal processor number
        MltMin(A, Y, Yproc);
        
        // if an entry has been modified, the loop is continued
        bool modif = false;
        for (int i = 0; i < m; i++)
          if ((Y(i) != Yprev(i)) || (Yproc(i) != Yproc_prev(i)))
            modif = true;
        
        MPI_Allreduce(&modif, &test_loop, 1, MPI_CHAR, MPI_LOR, comm);
      }
    
    // now we have to construct blocks from Y and Yproc
    Yprev.Clear(); Yproc_prev.Clear();
    
    IVect IndexRow(A.GetGlobalM());
    IndexRow.Fill(-1);
    
    IVect nb_block_send(nb_proc);
    nb_block_send.Fill(0);
    
    int nb_block = 0;
    for (int i = m1; i < m2; i++)
      {
        if (IndexRow(Y(i)) == -1)
          {
            nb_block_send(Yproc(i))++;
            IndexRow(Y(i)) = nb_block++;
          }
      }
    
    nb_block_send(rank) = 0;
    
    // exchanging the number of blocks
    IVect nb_block_recv(nb_proc);
    MPI_Alltoall(nb_block_send.GetData(), 1, MPI_INTEGER,
                 nb_block_recv.GetData(), 1, MPI_INTEGER, comm);
    
    // size of blocks
    IVect size_block(nb_block);
    size_block.Fill(0);
    for (int i = m1; i < m2; i++)
      size_block(IndexRow(Y(i)))++;
    
    // first version of blocks
    Vector<IVect> pattern(nb_block), proc_row(nb_block);
    IVect proc_block(nb_block), root_dof(nb_block);
    for (int i = 0; i < nb_block; i++)
      {
        pattern(i).Reallocate(size_block(i));
        proc_row(i).Reallocate(size_block(i));
        pattern(i).Fill(-1);
        proc_row(i).Fill(-1);
      }
    
    size_block.Fill(0);
    for (int i = m1; i < m2; i++)
      {
        int n = IndexRow(Y(i));
        pattern(n)(size_block(n)) = global(i);
        proc_row(n)(size_block(n)) = rank;
        proc_block(n) = Yproc(i);
        root_dof(n) = Y(i);
        size_block(n)++;
      }    
    
    // sending blocks    
    Vector<MPI_Request> request_send(nb_proc), request_send_int(nb_proc);
    Vector<IVect> Data_send(nb_block_send.GetM());
    IVect size_data_send(nb_block_send.GetM());
    for (int p = 0; p < nb_block_send.GetM(); p++)
      if ((p != rank) && (nb_block_send(p) > 0))
        {
          int N = nb_block_send(p);
          int sizeS = 2*N;
          for (int k = 0; k < proc_block.GetM(); k++)
            if (proc_block(k) == p)
              sizeS += size_block(k);
          
          size_data_send(p) = sizeS;
          MPI_Isend(&size_data_send(p), 1, MPI_INTEGER, p, 21, comm, &request_send_int(p));
          Data_send(p).Reallocate(sizeS);
          sizeS = 0;
          for (int k = 0; k < proc_block.GetM(); k++)
            if (proc_block(k) == p)
              {
                Data_send(p)(sizeS++) = root_dof(k);
                Data_send(p)(sizeS++) = size_block(k);
                for (int n = 0; n < size_block(k); n++)
                  Data_send(p)(sizeS++) = pattern(k)(n);
              }
          
          MPI_Isend(Data_send(p).GetData(), sizeS, MPI_INTEGER, p, 22, comm, &request_send(p));
        }
    
    // receiving blocks
    MPI_Status status;
    Vector<IVect> DataRecv(nb_proc);
    for (int p = 0; p < nb_block_recv.GetM(); p++)
      if ((p != rank) && (nb_block_recv(p) > 0))
        {
          int sizeS = 0;
          MPI_Recv(&sizeS, 1, MPI_INTEGER, p, 21, comm, &status);
          DataRecv(p).Reallocate(sizeS);
          MPI_Recv(DataRecv(p).GetData(), sizeS, MPI_INTEGER, p, 22, comm, &status);
        }
    
    for (int i = 0; i < request_send.GetM(); i++)
      if ((i != rank) && (nb_block_send(i) > 0))
        {
          MPI_Wait(&request_send_int(i), &status);
          MPI_Wait(&request_send(i), &status);
        }
    
    // treating distant unknowns
    Vector<IVect> list_block_recv(nb_proc);
    for (int p = 0; p < nb_block_recv.GetM(); p++)
      if ((p != rank) && (nb_block_recv(p) > 0))
        {
          int sizeS = 0;
          list_block_recv(p).Reallocate(nb_block_recv(p));
          for (int n = 0; n < nb_block_recv(p); n++)
            {
              int col = DataRecv(p)(sizeS++);
              int sizeB = DataRecv(p)(sizeS++);
              IVect num(sizeB), proc(sizeB);
              for (int k = 0; k < sizeB; k++)
                {
                  num(k) = DataRecv(p)(sizeS++);
                  proc(k) = p;
                }
              
              // adding these rows to pattern
              int j = IndexRow(col);
              pattern(j).PushBack(num);
              proc_row(j).PushBack(proc);
              list_block_recv(p)(n) = j;
            }
        }
    
    // eliminating duplicates in blocks
    for (int i = 0; i < pattern.GetM(); i++)
      RemoveDuplicate(pattern(i), proc_row(i));
    
    // sending back blocks in processors that sent their blocks
    Data_send.Reallocate(nb_block_recv.GetM());
    size_data_send.Reallocate(nb_block_recv.GetM());
    for (int p = 0; p < nb_block_recv.GetM(); p++)
      if ((p != rank) && (nb_block_recv(p) > 0))
        {          
          int sizeS = 0;
          for (int i2 = 0; i2 < list_block_recv(p).GetM(); i2++)
            sizeS += 2 + 2*pattern(list_block_recv(p)(i2)).GetM();                  
          
          size_data_send(p) = sizeS;
          Data_send(p).Reallocate(sizeS);
          sizeS = 0;
          for (int i2 = 0; i2 < list_block_recv(p).GetM(); i2++)
            {
              int i = list_block_recv(p)(i2);
              Data_send(p)(sizeS++) = root_dof(i);
              Data_send(p)(sizeS++) = pattern(i).GetM();
              for (int k = 0; k < pattern(i).GetM(); k++)
                {
                  Data_send(p)(sizeS++) = pattern(i)(k);
                  Data_send(p)(sizeS++) = proc_row(i)(k);
                }
            }
          
          MPI_Isend(&size_data_send(p), 1, MPI_INTEGER, p, 23, comm, &request_send_int(p));
          MPI_Isend(Data_send(p).GetData(), sizeS, MPI_INTEGER, p, 24, comm, &request_send(p));
        }
    
    for (int p = 0; p < nb_block_send.GetM(); p++)
      if ((p != rank) && (nb_block_send(p) > 0))
        {
          int sizeS = 0;
          MPI_Recv(&sizeS, 1, MPI_INTEGER, p, 23, comm, &status);
          DataRecv(p).Reallocate(sizeS);
          MPI_Recv(DataRecv(p).GetData(), sizeS, MPI_INTEGER, p, 24, comm, &status);
        }
    
    for (int i = 0; i < request_send.GetM(); i++)
      if ((i != rank) && (nb_block_recv(i) > 0))
        {
          MPI_Wait(&request_send_int(i), &status);
          MPI_Wait(&request_send(i), &status);
        }
    
    // updating blocks
    for (int p = 0; p < nb_block_send.GetM(); p++)
      if ((p != rank) && (nb_block_send(p) > 0))
        {
          int sizeS = 0;
          for (int n = 0; n < nb_block_send(p); n++)
            {
              int col = DataRecv(p)(sizeS++);
              int sizeB = DataRecv(p)(sizeS++);
              IVect num(sizeB), proc(sizeB);
              for (int k = 0; k < sizeB; k++)
                {
                  num(k) = DataRecv(p)(sizeS++);
                  proc(k) = DataRecv(p)(sizeS++);
                }
              
              // adding these rows to pattern
              int j = IndexRow(col);
              pattern(j) = num;
              proc_row(j) = proc;
            }
        }
    
    // using local numbers instead of global numbers 
    // for rows belonging to the current processor
    IndexRow.Fill(-1);
    for (int i = 0; i < global.GetM(); i++)
      IndexRow(global(i)) = i;
    
    for (int i = 0; i < pattern.GetM(); i++)
      for (int j = 0; j < pattern(i).GetM(); j++)
        if (IndexRow(pattern(i)(j)) >= 0)
          {
            pattern(i)(j) = IndexRow(pattern(i)(j));
            proc_row(i)(j) = rank;
          }
    
    IndexRow.Clear();
    
    B.Init(A);    
    B.SetPattern(pattern, proc_row);

    // filling local values
    FillBlockDiagonal(static_cast<const Matrix<T0, Prop0, Storage0, Allocator0>& >(A),
                      static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(B), m1, m2);
    
    // then filling distant values
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetDistantColSize(i); j++)
        {
          int jglob = A.IndexGlobalCol(i, j);
          int proc = A.ProcessorDistantCol(i, j);
          B.AddDistantInteraction(i, jglob, proc, A.ValueDistantCol(i, j));
        }

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i < A.GetDistantRowSize(j); i++)
        {
          int iglob = A.IndexGlobalRow(j, i);
          int proc = A.ProcessorDistantRow(j, i);
          B.AddRowDistantInteraction(iglob, j, proc, A.ValueDistantRow(j, i));
        }
  }

  
  //! converts a sparse matrix into a block-diagonal one
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRowSym, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B);
  }


  //! converts a sparse matrix into a block-diagonal one
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRow, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B);
  }
  
  
  //! converts a sparse matrix into a block-diagonal one
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRowSym, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B);
  }

  
  //! converts a sparse matrix into a block-diagonal one
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	    DistributedMatrix<T, Prop1, BlockDiagRow, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B);
  }

  
  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRowSym, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }
  

  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRow, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }


  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRowSym, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowSymComplexSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }


  //! converts a block-diagonal matrix into a sparse matrix
  template<class T, class Prop0, class Allocator0,
	   class Prop1, class Allocator1>
  void Copy(const DistributedMatrix<T, Prop0, BlockDiagRow, Allocator0>& A,
	    DistributedMatrix<T, Prop1, ArrayRowComplexSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }
  
  
  //! overwrites X by A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const DistributedMatrix<T1, Prop, BlockDiagRow, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X, bool assemble)
  {
    A.Mlt_(X, assemble);
  }


  //! overwrites X by A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const DistributedMatrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X, bool assemble)
  {
    A.Mlt_(X, assemble);
  }


  //! computes B = beta B + alpha A X
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const DistributedMatrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble)
  {
    A.MltAdd_(alpha, SeldonNoTrans, X, beta, B, assemble);
  }
  

  //! computes B = beta B + alpha A X
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const DistributedMatrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble)
  {
    A.MltAdd_(alpha, SeldonNoTrans, X, beta, B, assemble);
  }


  //! computes B = beta B + alpha A^T X  
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const DistributedMatrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble)
  {
    A.MltAdd_(alpha, trans, X, beta, B, assemble);
  }


  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const DistributedMatrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B, bool assemble)
  {
    A.MltAdd_(alpha, trans, X, beta, B, assemble);
  }
  
  
  //! A is replaced by its inverse (in assembled form)
  template<class T, class Allocator>
  void GetInverse(DistributedMatrix<T, General, BlockDiagRow, Allocator>& A)
  {
    // first, we assemble the matrix
    A.Assemble();
    // then each block can be inverted as for a sequential matrix
    GetInverse(static_cast<Matrix<T, General, BlockDiagRow, Allocator>& >(A));    
  }
    

  //! A is replaced by its inverse (in assembled form)
  template<class T, class Allocator>
  void GetInverse(DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    // first, we assemble the matrix
    A.Assemble();
    // then each block can be inverted as for a sequential matrix
    GetInverse(static_cast<Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& >(A));    
  }


  //! A is replaced by its Cholesky factorisation (in assembled form)
  template<class T, class Allocator>
  void GetCholesky(DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    A.Assemble();
    GetCholesky(static_cast<Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& >(A));    
  }
  
  
  //! solves L x = b or L^T x = b
  template<class T, class Allocator, class Allocator2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                     Vector<T, VectFull, Allocator2>& x)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        SolveCholesky(trans,
                      static_cast<const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& >(A), x);
        
        return;
      }
    
    // scattering columns values
    Vector<T> Xcol;
    A.ScatterColValues(x, Xcol);
    
    int nblock = A.GetNbBlocks(), nb = 0;
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
    int* proc_row_ = A.GetProcNumber();
    T val;
    if (trans.Trans())
      {
        // we solve L^T x = x
        for (int k = 0; k < nblock; k++)
          {
            nb = offset_[k+1]-1;
            for (int i = ptr_[k+1]-1; i >= ptr_[k]; i--)
              {
                int irow = ind_[i];
                if (proc_row_[i] == rank)
                  val = x(irow);
                else
                  val = Xcol(irow);
                
                for (int j = ptr_[k+1]-1; j >= i+1; j--)
                  {
                    if (proc_row_[j] == rank)
                      val -= data_[nb--]*x(ind_[j]);
                    else
                      val -= data_[nb--]*Xcol(ind_[j]);
                  }
                
                if (proc_row_[i] == rank)
                  x(irow) = val/data_[nb--];
                else
                  Xcol(irow) = val/data_[nb--];
              }
          }
      }
    else
      {
        // we solve L x = x
        for (int k = 0; k < nblock; k++)
          {
            nb = offset_[k];
            for (int i = ptr_[k]; i < ptr_[k+1]; i++)
              {
                int irow = ind_[i];
                if (proc_row_[i] == rank)
                  {
                    x(irow) /= data_[nb++];
                    val = x(irow);
                  }
                else
                  {
                    Xcol(irow) /= data_[nb++];
                    val = Xcol(irow);
                  }
                
                for (int j = i+1; j < ptr_[k+1]; j++)
                  {
                    if (proc_row_[j] == rank)
                      x(ind_[j]) -= data_[nb++]*val;                
                    else
                      Xcol(ind_[j]) -= data_[nb++]*val;
                  }
              }
          }
      }
  }
  
  
  //! multiplying L x or L^T x
  template<class T, class Allocator, class Allocator2>
  void MltCholesky(const SeldonTranspose& trans,
                   const DistributedMatrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                   Vector<T, VectFull, Allocator2>& x)
  {
    const MPI_Comm& comm = A.GetCommunicator();
    int rank; MPI_Comm_rank(comm, &rank);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      {
        MltCholesky(trans,
                    static_cast<const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& >(A), x);
        
        return;
      }

    // scattering columns values
    Vector<T> Xcol;
    A.ScatterColValues(x, Xcol);
    
    int nblock = A.GetNbBlocks(), nb = 0;
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
    int* proc_row_ = A.GetProcNumber();
    T val;
    if (trans.Trans())
      {
        // we overwrite x by L^T x
        for (int k = 0; k < nblock; k++)
          {
            nb = offset_[k];
            for (int i = ptr_[k]; i < ptr_[k+1]; i++)
              {
                int irow = ind_[i];
                if (proc_row_[i] == rank)
                  val = x(irow)*data_[nb++];
                else
                  val = Xcol(irow)*data_[nb++];
                
                for (int j = i+1; j < ptr_[k+1]; j++)
                  {
                    if (proc_row_[j] == rank)
                      val += data_[nb++]*x(ind_[j]);
                    else
                      val += data_[nb++]*Xcol(ind_[j]);
                  }
                
                if (proc_row_[i] == rank)
                  x(irow) = val; 
                else
                  Xcol(irow) = val; 
              }
          }
      }
    else
      {
        // we overwrite x by L x
        for (int k = 0; k < nblock; k++)
          {
            nb = offset_[k+1]-1;
            for (int i = ptr_[k+1]-1; i >= ptr_[k]; i--)
              {
                int irow = ind_[i];
                if (proc_row_[i] == rank)
                  val = x(irow);
                else
                  val = Xcol(irow);
                
                for (int j = ptr_[k+1]-1; j >= i+1; j--)
                  {
                    if (proc_row_[j] == rank)
                      x(ind_[j]) += data_[nb--]*val;
                    else
                      Xcol(ind_[j]) += data_[nb--]*val;
                  }
                  
                if (proc_row_[i] == rank)
                  x(irow) *= data_[nb--];
                else
                  Xcol(irow) *= data_[nb--];
              }
          }
      }
  }


  template<class T1, class Prop1, class Allocator1>
  void EraseRow(const IVect& num_,
		DistributedMatrix<T1, Prop1, BlockDiagRow, Allocator1>& A)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return EraseRow(num_, static_cast<Matrix<T1, Prop1, BlockDiagRow, Allocator1>& >(A));
    
    // we retrieve the indexes in the array A.ind_
    IVect num(num_);
    A.CompleteWithAssociatedNumbers(num);
          
    // required entries are set to 0
    int k = 0;
    T1 zero; SetComplexZero(zero);
    for (int i = 0; i < num.GetM(); i++)
      {
	while ((k < A.GetNbBlocks()) && (num(i) >= A.ptr_[k+1]))
	  k++;

	int j = num(i) - A.ptr_[k];
	int size_row = A.ptr_[k+1] - A.ptr_[k];
	for (int q = 0; q < size_row; q++)
	  A.data_[A.offset_[k]  + size_row*j + q] = zero;
      }
  }


  template<class T1, class Prop1, class Allocator1>
  void EraseRow(const IVect& num_,
		DistributedMatrix<T1, Prop1, BlockDiagRowSym, Allocator1>& A)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return EraseRow(num_, static_cast<Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& >(A));
    
    // we retrieve the indexes in the array A.ind_
    IVect num(num_);
    A.CompleteWithAssociatedNumbers(num);
          
    // required entries are set to 0
    int k = 0;
    T1 zero; SetComplexZero(zero);
    for (int i = 0; i < num.GetM(); i++)
      {
	while ((k < A.GetNbBlocks()) && (num(i) >= A.ptr_[k+1]))
	  k++;

	int j = num(i) - A.ptr_[k];
	int size_row = A.ptr_[k+1] - A.ptr_[k];
	for (int q = 0; q <= j; q++)
	  A.data_[A.offset_[k]  + size_row*q - q*(q+1)/2 + j] = zero;
	
	for (int q = j+1; q < size_row; q++)
	  A.data_[A.offset_[k]  + size_row*j - j*(j+1)/2 + q] = zero;
      }
  }


  template<class T1, class Prop1, class Allocator1>
  void EraseCol(const IVect& num_,
		DistributedMatrix<T1, Prop1, BlockDiagRow, Allocator1>& A)
  {
    MPI_Comm& comm = A.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc == 1)
      return EraseCol(num_, static_cast<Matrix<T1, Prop1, BlockDiagRow, Allocator1>& >(A));
    
    // we retrieve the indexes in the array A.ind_
    IVect num(num_);
    A.CompleteWithAssociatedNumbers(num);
          
    // required entries are set to 0
    int k = 0;
    T1 zero; SetComplexZero(zero);
    for (int i = 0; i < num.GetM(); i++)
      {
	while ((k < A.GetNbBlocks()) && (num(i) >= A.ptr_[k+1]))
	  k++;

	int j = num(i) - A.ptr_[k];
	int size_row = A.ptr_[k+1] - A.ptr_[k];
	for (int q = 0; q < size_row; q++)
	  A.data_[A.offset_[k]  + size_row*q + j] = zero;
      }
  }


  template<class T1, class Prop1, class Allocator1>
  void EraseCol(const IVect& num,
		DistributedMatrix<T1, Prop1, BlockDiagRowSym, Allocator1>& A)
  {
    // by symmetry, rows and columns are both erased here
    EraseRow(num, A);
  }

  
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const DistributedMatrix<T, Prop, BlockDiagRow, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }

  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out,
                      const DistributedMatrix<T, Prop, BlockDiagRowSym, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }
  
}

#define SELDON_FILE_DISTRIBUTED_BLOCK_DIAGONAL_MATRIX_CXX
#endif
