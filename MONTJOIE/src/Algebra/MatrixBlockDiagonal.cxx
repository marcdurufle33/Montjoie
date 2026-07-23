#ifndef SELDON_FILE_MATRIX_BLOCK_DIAGONAL_CXX

#include "MatrixBlockDiagonal.hxx"

namespace Seldon
{

  /************************
   * Matrix_BlockDiagonal *
   ************************/
  
  
  //! Default constructor.
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  Matrix_BlockDiagonal() : Matrix_Base<T, Allocator>()
  {
    ptr_ = NULL;
    ind_ = NULL;
    rev_ = NULL;
    num_blok_ = NULL;
    offset_ = NULL;
    nz_ = 0;
    nblok_ = 0;
    size_ind_ = 0;
    size_max_blok_ = 0;
  }
  
  
  //! copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  Matrix_BlockDiagonal(const Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    ptr_ = NULL;
    ind_ = NULL;
    rev_ = NULL;
    num_blok_ = NULL;
    offset_ = NULL;
    nz_ = 0;
    nblok_ = 0;
    size_ind_ = 0;
    this->m_ = 0;
    this->n_ = 0;
    this->data_ = NULL;
    size_max_blok_ = 0;
    
    Copy(A);
  }
  
  
  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::~Matrix_BlockDiagonal()
  {
    Clear();
  }
  
  
  //! clears the matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Clear()
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    Allocator::deallocate(this->data_, nz_);
	    this->data_ = NULL;
	  }
	
	if (this->ind_ != NULL)
	  {
	    AllocatorInt::deallocate(this->ind_, this->size_ind_);
	    this->ind_ = NULL;
	  }
	
	if (this->rev_ != NULL)
	  {
	    AllocatorInt::deallocate(this->rev_, this->m_);
	    this->rev_ = NULL;
	  }
	
	if (this->ptr_ != NULL)
	  {
	    AllocatorInt::deallocate(this->ptr_, nblok_+1);
	    this->ptr_ = NULL;
	  }
	
	if (this->num_blok_ != NULL)
	  {
	    AllocatorInt::deallocate(this->num_blok_, this->m_);
	    this->num_blok_ = NULL;
	  }
	
	if (this->offset_ != NULL)
	  {
	    AllocatorInt::deallocate(this->offset_, this->nblok_+1);
	    this->offset_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
	ptr_ = NULL;
	ind_ = NULL;
	rev_ = NULL;
	num_blok_ = NULL;
      }
#endif
    
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0; this->nblok_ = 0;
    size_ind_ = 0;
    size_max_blok_ = 0;
  }
  

  //! sets to zero a row i
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::ClearRow(int i)
  {
    int n = this->num_blok_[i];
    int j = this->rev_[i];
    int size_row = ptr_[n+1] - ptr_[n];
    T zero; SetComplexZero(zero);
    for (int k = 0; k < size_row; k++)
      this->data_[offset_[n]  + size_row*j + k] = zero;
  }


  //! sets to zero a column i
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::ClearColumn(int i)
  {
    int n = this->num_blok_[i];
    int j = this->rev_[i];
    int size_row = ptr_[n+1] - ptr_[n];
    T zero; SetComplexZero(zero);
    for (int k = 0; k < size_row; k++)
      this->data_[offset_[n]  + size_row*k + j] = zero;
  }

  
  //! Initialisation of the profile of the matrix
  /*!
    \param[in] num_dof_blocks row numbers for all the blocks
    \param[in] size_blocks size of each block
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const IVect& num_dof_blocks, const IVect& size_blocks,
             const Vector<bool>& ignore_row, IVect& permut)
  {
    // clears previous entries
    Clear();
    
    // number of non-zero entries
    int nnz = 0;
    bool too_large = false;
    for (int i = 0; i < size_blocks.GetM(); i++)
      {
        int nnz_block = size_blocks(i)*size_blocks(i);
        nnz += nnz_block;
        if ((nnz_block < 0) || (nnz < 0))
          too_large = true;
      }
    
    if (too_large)
      {
	cout<<"Block-matrix too large"<<endl;
	abort();
      }
    
    // counting m
    this->size_ind_ = num_dof_blocks.GetM();
    this->m_ = 0;
    for (int i = 0; i < ignore_row.GetM(); i++)
      if (!ignore_row(i))
        this->m_++;
    
    this->n_ = this->m_;

    // allocating the arrays
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
        
        this->ind_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->ind_,
                                                                 this->size_ind_, this) );
        
        this->rev_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->rev_,
                                                                 this->m_, this) );
        
        this->num_blok_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->num_blok_,
                                                                 this->m_, this) );
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        this->m_ = 0;
        this->n_ = 0;
        this->size_ind_ = 0;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
    if ((this->ind_ == NULL)||(this->rev_ == NULL))
      {
        this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
#endif
    
    this->nblok_ = size_blocks.GetM();
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

        this->ptr_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->ptr_,
                                                                 this->nblok_+1, this) );
        
        this->offset_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->offset_,
                                                                 this->nblok_+1, this) );
        
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        this->m_ = 0; this->n_ = 0;
        this->nblok_ = 0;
        this->size_ind_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
    if (this->ptr_ == NULL)
      {
        this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
        this->nblok_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
#endif
    
    this->nz_ = nnz;
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
        
        this->data_ =
          reinterpret_cast<pointer>(Allocator::reallocate(this->data_,
                                                                this->nz_, this) );
        
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        this->m_ = 0; this->n_ = 0;
        this->nz_ = 0;
        this->size_ind_ = 0;
        this->nblok_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->data_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
    
    if (this->data_ == NULL)
      {
        this->m_ = 0; this->n_ = 0;
        this->nz_ = 0;
        this->size_ind_ = 0;
        this->nblok_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->data_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
#endif
    
    // now filling the arrays
    if ((this->nblok_ > 0)&&(this->m_ > 0))
      {
        size_max_blok_ = size_blocks.GetNormInf();
	this->ptr_[0] = 0;
        this->offset_[0] = 0;
	for (int i = 0; i < this->nblok_; i++)
	  {
	    this->ptr_[i+1] = this->ptr_[i] + size_blocks(i);
	    this->offset_[i+1] = this->offset_[i] + size_blocks(i)*size_blocks(i);
	  }
	
	for (int i = 0; i < this->size_ind_; i++)
	  this->ind_[i] = num_dof_blocks(i);
	
	for (int i = 0; i < this->nblok_; i++)
	  for (int j = this->ptr_[i]; j < this->ptr_[i+1]; j++)
	    {
	      if (!ignore_row(j))
                {
                  this->num_blok_[num_dof_blocks(j)] = i;
                  this->rev_[num_dof_blocks(j)] = j-this->ptr_[i];
                }
	    }
	
        permut.Reallocate(this->nblok_);
        permut.Fill();
      }
    
    // we set non-zero entries to 0
    T zero; SetComplexZero(zero);
    Fill(zero);
  }


  //! specifies all the blocks of the matrix
  /*!
    \param[in] NumDof_Blocks an array containing the row numbers for each block of the matrix
    this methods clears the matrix and affects a new structure, while setting to 0 all the entries
    of the matrix
   */
  template <class T, class Prop, class Storage, class Allocator> template<class Allocator1>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const Vector<IVect, Vect_Full, Allocator1>& NumDof_Blocks,
             const Vector<Vector<bool> >& ignore_row, IVect& permut)
  {
    int nb_blocks = NumDof_Blocks.GetM();
    IVect size_blocks(nb_blocks);
    int nodl = 0;
    for (int i = 0; i < nb_blocks; i++)
      {
	size_blocks(i) = NumDof_Blocks(i).GetM();
	nodl += size_blocks(i);
      }
    
    Vector<bool> ignore_row_flat(nodl);
    IVect num_dofs(nodl); int nb = 0;
    ignore_row_flat.Fill(false);
    for (int i = 0; i < nb_blocks; i++)
      {
	for (int j = 0; j < size_blocks(i); j++)
	  {
            if (ignore_row.GetM() > 0)
              ignore_row_flat(nb) = ignore_row(i)(j);
            
            num_dofs(nb) = NumDof_Blocks(i)(j);
            nb++;
          }
      }
    
    SetPattern(num_dofs, size_blocks, ignore_row_flat, permut);
  }

  
  //! specifies all the blocks of the matrix
  /*!
    \param[in] NumDof_Blocks an array containing the row numbers for each block of the matrix
    this methods clears the matrix and affects a new structure, while setting to 0 all the entries
    of the matrix
   */ 
  template <class T, class Prop, class Storage, class Allocator> template<class Allocator1>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const Vector<IVect, Vect_Full, Allocator1>& NumDof_Blocks)
  {
    Vector<Vector<bool> > ignore_row; IVect permut;
    SetPattern(NumDof_Blocks, ignore_row, permut);
  }
  
  
  //! present for compatibility purpose
  template<class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    // empty, allocation is made by SetPattern method
  }
  
  
  //! present for compatibility purpose
  template<class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    // empty, you can't change size of matrix
  }
  
  
  //! copies object
  template<class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  Copy(const Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& A)
  {
    Clear();
    
    this->m_ = A.GetM();
    this->n_ =  A.GetN();
    this->size_ind_ = A.size_ind_;
    
    if (this->m_ <= 0)
      return;
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	this->ind_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->ind_,
								 this->size_ind_, this) );
	    
	this->rev_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->rev_,
								 this->m_, this) );
	
	this->num_blok_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->num_blok_,
								 this->m_, this) );
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
        this->size_ind_ = 0;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
    if ((this->ind_ == NULL)||(this->rev_ == NULL))
      {
	this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
#endif

    this->nblok_ = A.GetNbBlocks();

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	this->ptr_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->ptr_,
							  this->nblok_+1, this) );
	
	this->offset_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->offset_,
							  this->nblok_+1, this) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0; this->n_ = 0;
	this->nblok_ = 0;
        this->size_ind_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
    if (this->ptr_ == NULL)
      {
	this->m_ = 0; this->n_ = 0;
	this->nblok_ = 0;
        this->size_ind_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
#endif
    
    this->nz_ = A.GetDataSize();

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	this->data_ =
	  reinterpret_cast<pointer>(Allocator::reallocate(this->data_,
								this->nz_, this) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0; this->n_ = 0;
	this->nz_ = 0;
        this->size_ind_ = 0;
	this->nblok_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->data_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0; this->n_ = 0;
	this->nz_ = 0;
        this->size_ind_ = 0;
	this->nblok_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->data_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
#endif

    // now filling the arrays    
    this->size_max_blok_ = A.GetSizeMaxBlock();
    AllocatorInt::memorycpy(this->ptr_, A.GetPtr(), this->nblok_+1);
    AllocatorInt::memorycpy(this->offset_, A.GetOffset(), this->nblok_+1);

    AllocatorInt::memorycpy(this->ind_, A.GetInd(), this->size_ind_);
    AllocatorInt::memorycpy(this->rev_, A.GetRev(), this->m_);
    AllocatorInt::memorycpy(this->num_blok_, A.GetNumBlock(), this->m_);

    Allocator::memorycpy(this->data_, A.GetData(), this->nz_);

  }
  

  //! returns A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  const typename Matrix_BlockDiagonal<T,Prop,Storage,Allocator>::value_type
  Matrix_BlockDiagonal<T,Prop,Storage,Allocator>::operator() (int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    T zero; SetComplexZero(zero);
    int n1 = this->num_blok_[i];
    int n2 = this->num_blok_[j];
    if (n1 == n2)
      {
	int row = this->rev_[i];
	int col = this->rev_[j];
	return this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1])+col];
      }
    
    return zero;
  }


  //! returns reference to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_BlockDiagonal<T,Prop,Storage,Allocator>::reference
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    int n1 = this->num_blok_[i];
    int n2 = this->num_blok_[j];
    if (n1 == n2)
      {
	int row = this->rev_[i];
	int col = this->rev_[j];
	return this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1])+col];
      }
    else
      {
        cout << "Entry not present in the matrix" << endl;
        abort();
      }
    
    return this->data_[0];
  }


  //! returns reference to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  typename Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::const_reference
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    int n1 = this->num_blok_[i];
    int n2 = this->num_blok_[j];
    if (n1 == n2)
      {
	int row = this->rev_[i];
	int col = this->rev_[j];
	return this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1])+col];
      }
    else
      {
        cout << "Entry not present in the matrix" << endl;
        abort();
      }
    
    return this->data_[0];
  }


  //! Sets A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Set(int i, int j, const T& val)
  {
    this->Get(i, j) = val;
  }


  //! Sets A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::SetEntry(int i, int j, const T& val)
  {
    this->Get(i, j) = val;
  }
  
    
  //! sets all non-zero entries to a given value
  template <class T, class Prop, class Storage, class Allocator> template<class T0>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (int i = 0; i < this->nz_; i++)
      this->data_[i] = x_;
  }
  
  
  //! sets randomly non-zero entries
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::FillRand()
  {
    for (int i = 0; i < this->nz_; i++)
      SetComplexReal(rand(), this->data_[i]);
  }
  
  
  //! sets non-zero entries to 0
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::Zero()
  {
    T zero; SetComplexZero(zero);
    for (int i = 0; i < this->nz_; i++)
      this->data_[i] = zero;
  }
  
  
  //! *this = A
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& 
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  operator= (const Matrix_BlockDiagonal<T, Prop, Storage, Allocator>& A)
  {
    Copy(A);
    return *this;
  }
  
  
  //! Multiplies non-zero entries of the matrix by a coefficient alpha
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0> Matrix_BlockDiagonal<T, Prop, Storage, Allocator>&
  Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->nz_; i++)
      this->data_[i] *= alpha;
    
    return *this;
  }
  
  
  //! adds terms to a row of the matrix
  /*!
    \param[in] i row where values are added
    \param[in] nb_interac number of values to add to the matrix
    \param[in] col_interac column numbers associated with values
    \param[in] val_itnterac values to add
    Equivalent Matlab operation A(i, col_interac) = A(i, col_interac) + val_interac
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::
  AddInteractionRow(int i, int nb_interac, const IVect& col_interac,
                    const Vector<T>& val_interac, bool sorted)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
#endif
    int n1 = this->num_blok_[i];
    int row = this->rev_[i];
    for (int j1 = 0; j1 < nb_interac; j1++)
      {
	int j = col_interac(j1);
#ifdef SELDON_CHECK_BOUNDS
	if (j < 0 || j >= this->m_)
	  throw WrongCol("Matrix_BlockDiagonal::operator()",
			 string("Index should be in [0, ")
			 + to_str(this->m_-1) + "], but is equal to "
			 + to_str(j) + ".");
#endif
	int n2 = this->num_blok_[j];
	if (n2 == n1)
	  {
	    int col = this->rev_[j];
	    this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1])+col] += val_interac(j1);
	  }
      }
  }
  
  
  //! Writes matrix in Matlab format (i j val)
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>
  ::WriteText(const string& file_name) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_BlockDiagonal::WriteText(string FileName)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! Writes matrix in Matlab format (i j val)
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_BlockDiagonal<T, Prop, Storage, Allocator>::WriteText(ostream& FileStream) const
  {
    int nb = 0, irow, icol;
    for (int i = 0; i < nblok_; i++)
      {
	int size_block = ptr_[i+1] - ptr_[i];
	for (int j = 0; j < size_block; j++)
	  {
	    irow = ind_[j + ptr_[i] ];
	    for (int k = 0; k < size_block; k++)
	      {
		icol = ind_[k + ptr_[i] ];
		FileStream << irow+1 << ' ' << icol+1 << ' ' << this->data_[nb++] << '\n';
	      }
	  }
      }
  }
  

  /***************************
   * MatrixSym_BlockDiagonal *
   ***************************/

  
  //! Default constructor.
  template <class T, class Prop, class Storage, class Allocator>
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  MatrixSym_BlockDiagonal() : Matrix_Base<T,Allocator>()
  {
    this->nz_ = 0;
    nblok_ = 0;
    size_ind_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
    rev_ = NULL;
    num_blok_ = NULL;
    offset_ = NULL;
    size_max_blok_ = 0;
  }
  
  
  //! copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  MatrixSym_BlockDiagonal(const MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0;
    this->data_ = NULL;
    size_ind_ = 0;
    nblok_ = 0;
    
    ptr_ = NULL;
    ind_ = NULL;
    rev_ = NULL;
    num_blok_ = NULL;
    offset_ = NULL;
    size_max_blok_ = 0;
    
    Copy(A);
  }
  
  
  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::~MatrixSym_BlockDiagonal()
  {
    this->Clear();
  }
  
  
  //! clears the matrix
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T,Prop,Storage,Allocator>::Clear()
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    Allocator::deallocate(this->data_, nz_);
	    this->data_ = NULL;
	  }
	
	if (this->ind_ != NULL)
	  {
	    AllocatorInt::deallocate(this->ind_, this->size_ind_);
	    this->ind_ = NULL;
	  }
	
	if (this->rev_ != NULL)
	  {
	    AllocatorInt::deallocate(this->rev_, this->m_);
	    this->rev_ = NULL;
	  }
	
	if (this->ptr_ != NULL)
	  {
	    AllocatorInt::deallocate(this->ptr_, nblok_+1);
	    this->ptr_ = NULL;
	  }
	
	if (this->num_blok_ != NULL)
	  {
	    AllocatorInt::deallocate(this->num_blok_, this->m_);
	    this->num_blok_ = NULL;
	  }
	
	if (this->offset_ != NULL)
	  {
	    AllocatorInt::deallocate(this->offset_, this->nblok_+1);
	    this->offset_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
	ptr_ = NULL;
	ind_ = NULL;
	rev_ = NULL;
	num_blok_ = NULL;
	offset_ = NULL;
      }
#endif
    
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0; this->nblok_ = 0;
    size_ind_ = 0;
    size_max_blok_ = 0;
  }
  

  //! sets to zero a row i
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::ClearRow(int i)
  {
    int n = this->num_blok_[i];
    int j = this->rev_[i];
    int size_row = ptr_[n+1] - ptr_[n];
    T zero; SetComplexZero(zero);
    for (int k = 0; k <= j; k++)
      this->data_[offset_[n]  + size_row*k - k*(k+1)/2 + j] = zero;
    
    for (int k = j+1; k < size_row; k++)
      this->data_[offset_[n]  + size_row*j - j*(j+1)/2 + k] = zero;
  }


  //! sets to zero a column i
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::ClearColumn(int i)
  {
    ClearRow(i);
  }

    
  //! Initialisation of the profile of the matrix
  /*!
    \param[in] num_dof_blocks row numbers for all the blocks
    \param[in] size_blocks size of each block
   */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const IVect& num_dof_blocks, const IVect& size_blocks,
             const Vector<bool>& ignore_row, IVect& permut)
  {
    // clears previous entries
    Clear();
    
    // number of non-zero entries
    bool too_large = false;
    int nnz = 0;
    for (int i = 0; i < size_blocks.GetM(); i++)
      {
	int nnz_block = size_blocks(i)*(size_blocks(i)+1)/2;
        nnz += nnz_block;
        if ((nnz_block < 0) || (nnz < 0))
          too_large = true;
      }
    
    if (too_large)
      {
	cout<<"Block-matrix too large"<<endl;
	abort();
      }
    // counting m
    this->size_ind_ = num_dof_blocks.GetM();
    this->m_ = 0;
    for (int i = 0; i < ignore_row.GetM(); i++)
      if (!ignore_row(i))
        this->m_++;
    
    this->n_ = this->m_;
    
    // allocating arrays
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
        
        this->ind_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->ind_,
                                                                 this->size_ind_, this) );
        
        this->rev_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->rev_,
                                                                 this->m_, this) );
        
        this->num_blok_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->num_blok_,
                                                                 this->m_, this) );
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
    if ((this->ind_ == NULL)||(this->rev_ == NULL))
      {
        this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
#endif
    
    this->nblok_ = size_blocks.GetM();
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
        
        this->ptr_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->ptr_,
                                                                 this->nblok_+1, this) );
        
        this->offset_ =
          reinterpret_cast<int*>(AllocatorInt::reallocate(this->offset_,
                                                                 this->nblok_+1, this) );
        
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        this->m_ = 0; this->n_ = 0;
        this->nblok_ = 0;
        this->size_ind_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
    if (this->ptr_ == NULL)
      {
        this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
        this->nblok_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
#endif
    
    this->nz_ = nnz;
    
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
        
        this->data_ =
          reinterpret_cast<pointer>(Allocator::reallocate(this->data_,
                                                                this->nz_, this) );
        
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        this->m_ = 0; this->n_ = 0;
        this->nz_ = 0;
        this->size_ind_ = 0;
        this->nblok_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->data_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
    if (this->data_ == NULL)
      {
        this->m_ = 0; this->n_ = 0;
        this->nz_ = 0;
        this->size_ind_ = 0;
        this->nblok_ = 0;
        this->ptr_ = NULL;
        this->offset_ = NULL;
        this->ind_ = NULL;
        this->rev_ = NULL;
        this->data_ = NULL;
        this->num_blok_ = NULL;
        return;
      }
#endif
    
    // now filling the arrays
    if ((this->nblok_ > 0)&&(this->m_ > 0))
      {
        this->size_max_blok_ = size_blocks.GetNormInf();
	this->ptr_[0] = 0; this->offset_[0] = 0;
	for (int i = 0; i < this->nblok_; i++)
	  {
	    this->ptr_[i+1] = this->ptr_[i] + size_blocks(i);
	    this->offset_[i+1] = this->offset_[i] + size_blocks(i)*(size_blocks(i)+1)/2;
	  }
	
	for (int i = 0; i < this->size_ind_; i++)
	  this->ind_[i] = num_dof_blocks(i);
	
	for (int i = 0; i < this->nblok_; i++)
	  for (int j = this->ptr_[i]; j < this->ptr_[i+1]; j++)
	    {
	      if (!ignore_row(j))
                {
                  this->num_blok_[num_dof_blocks(j)] = i;
                  this->rev_[num_dof_blocks(j)] = j-this->ptr_[i];
                }
	    }	

        permut.Reallocate(this->nblok_);
        permut.Fill();
      }
    
    // we set non-zero entries to 0
    T zero; SetComplexZero(zero);
    Fill(zero);
  }
  
  
  //! specifies all the blocks of the matrix
  /*!
    \param[in] NumDof_Blocks an array containing the row numbers for each block of the matrix
    this methods clears the matrix and affects a new structure, while setting to 0 all the entries
    of the matrix
   */
  template <class T, class Prop, class Storage, class Allocator> template<class Allocator1>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const Vector<IVect, Vect_Full, Allocator1>& NumDof_Blocks,
             const Vector<Vector<bool> >& ignore_row, IVect& permut)
  {
    int nb_blocks = NumDof_Blocks.GetM();
    IVect size_blocks(nb_blocks);
    int nodl = 0;
    for (int i = 0; i < nb_blocks; i++)
      {
	size_blocks(i) = NumDof_Blocks(i).GetM();
	nodl += size_blocks(i);
      }
    
    Vector<bool> ignore_row_flat(nodl);
    IVect num_dofs(nodl); int nb = 0;
    ignore_row_flat.Fill(false);
    for (int i = 0; i < nb_blocks; i++)
      {
	for (int j = 0; j < size_blocks(i); j++)
	  {
            if (ignore_row.GetM() > 0)
              ignore_row_flat(nb) = ignore_row(i)(j);
            
            num_dofs(nb) = NumDof_Blocks(i)(j);
            nb++;
          }
      }
    
    SetPattern(num_dofs, size_blocks, ignore_row_flat, permut);
  }

  
  //! specifies all the blocks of the matrix
  /*!
    \param[in] NumDof_Blocks an array containing the row numbers for each block of the matrix
    this methods clears the matrix and affects a new structure, while setting to 0 all the entries
    of the matrix
   */ 
  template <class T, class Prop, class Storage, class Allocator> template<class Allocator1>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const Vector<IVect, Vect_Full, Allocator1>& NumDof_Blocks)
  {
    Vector<Vector<bool> > ignore_row; IVect permut;
    SetPattern(NumDof_Blocks, ignore_row, permut);
  }
  

  //! present for compatibility purpose  
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    // empty, allocation is made by SetPattern method
  }
  
  
  //! present for compatibility purpose
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    // empty, you can't change size of matrix
  }
  
  
  //! copies object
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  Copy(const MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& A)
  {
    Clear();
    
    this->m_ = A.GetM();
    this->n_ =  A.GetN();
    this->size_ind_ = A.size_ind_;
    if (this->m_ <= 0)
      return;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	this->ind_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->ind_,
								 this->size_ind_, this) );
	    
	this->rev_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->rev_,
								 this->m_, this) );
	
	this->num_blok_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->num_blok_,
								 this->m_, this) );
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
        this->size_ind_ = 0;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
    if ((this->ind_ == NULL)||(this->rev_ == NULL))
      {
	this->m_ = 0; this->n_ = 0;
        this->size_ind_ = 0;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
#endif

    this->nblok_ = A.GetNbBlocks();

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	this->ptr_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->ptr_,
								 this->nblok_+1, this) );
	
	this->offset_ =
	  reinterpret_cast<int*>(AllocatorInt::reallocate(this->offset_,
								 this->nblok_+1, this) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0; this->n_ = 0;
	this->nblok_ = 0;
        this->size_ind_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
    if (this->ptr_ == NULL)
      {
	this->m_ = 0; this->n_ = 0;
	this->nblok_ = 0;
        this->size_ind_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
#endif
    
    this->nz_ = A.GetDataSize();

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	this->data_ =
	  reinterpret_cast<pointer>(Allocator::reallocate(this->data_,
								this->nz_, this) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0; this->n_ = 0;
	this->nz_ = 0;
        this->size_ind_ = 0;
	this->nblok_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->data_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0; this->n_ = 0;
	this->nz_ = 0;
        this->size_ind_ = 0;
	this->nblok_ = 0;
	this->ptr_ = NULL;
	this->offset_ = NULL;
	this->ind_ = NULL;
	this->rev_ = NULL;
	this->data_ = NULL;
	this->num_blok_ = NULL;
	return;
      }
#endif

    // now filling the arrays    
    this->size_max_blok_ = A.GetSizeMaxBlock();
    AllocatorInt::memorycpy(this->ptr_, A.GetPtr(), this->nblok_+1);
    AllocatorInt::memorycpy(this->offset_, A.GetOffset(), this->nblok_+1);

    AllocatorInt::memorycpy(this->ind_, A.GetInd(), this->size_ind_);
    AllocatorInt::memorycpy(this->rev_, A.GetRev(), this->m_);
    AllocatorInt::memorycpy(this->num_blok_, A.GetNumBlock(), this->m_);

    Allocator::memorycpy(this->data_, A.GetData(), this->nz_);

  }
  

  //! returns A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  typename MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::value_type
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::operator() (int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    T zero; SetComplexZero(zero);
    int n1 = this->num_blok_[i];
    int n2 = this->num_blok_[j];
    if (n1 == n2)
      {
	int row = this->rev_[i]; 
	int col = this->rev_[j];
	if (row > col)
	  {
	    int itmp = row;
	    row = col;
	    col = itmp;
	  }
	
	return this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1]) - row*(row+1)/2 + col];
      }
    
    return zero;
  }


  //! returns a reference to A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  typename MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::reference
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    T zero; SetComplexZero(zero);
    int n1 = this->num_blok_[i];
    int n2 = this->num_blok_[j];
    if (n1 == n2)
      {
	int row = this->rev_[i]; 
	int col = this->rev_[j];
	if (row > col)
	  {
	    int itmp = row;
	    row = col;
	    col = itmp;
	  }
	
	return this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1]) - row*(row+1)/2 + col];
      }
    else
      {
	cout << "Entry not present in the matrix" << endl;
        abort();
      }
    
    return this->data_[0];
  }


  //! returns A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  typename MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::const_reference
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    T zero; SetComplexZero(zero);
    int n1 = this->num_blok_[i];
    int n2 = this->num_blok_[j];
    if (n1 == n2)
      {
	int row = this->rev_[i]; 
	int col = this->rev_[j];
	if (row > col)
	  {
	    int itmp = row;
	    row = col;
	    col = itmp;
	  }
	
	return this->data_[offset_[n1]+row*(ptr_[n1+1]-ptr_[n1]) - row*(row+1)/2 + col];
      }
    else
      {
	cout << "Entry not present in the matrix" << endl;
        abort();
      }
    
    return this->data_[0];
  }


  //! Sets A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Set(int i, int j, const T& val)
  {
    this->Get(i, j) = val;
  }


  //! Sets A(i, j)
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::SetEntry(int i, int j, const T& val)
  {
    this->Get(i, j) = val;
  }
    
    
  //! sets all non-zero entries to a given value
  template <class T, class Prop, class Storage, class Allocator> template<class T0>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (int i = 0; i < this->nz_; i++)
      this->data_[i] = x_;
  }
  
  
  //! sets randomly non-zero entries
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::FillRand()
  {
    for (int i = 0; i < this->nz_; i++)
      SetComplexReal(rand(), this->data_[i]);
  }

  
  //! sets non-zero entries to 0
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Zero()
  {
    T zero; SetComplexZero(zero);
    for (int i = 0; i < this->nz_; i++)
      this->data_[i] = zero;
  }

  
  //! copies object
  template <class T, class Prop, class Storage, class Allocator>
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& 
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>
  ::operator= (const MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& A)
  {
    Copy(A);
    return *this;
  }

  
  //! Multiplies non-zero entries of the matrix by a coefficient alpha
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0> MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>& 
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->nz_; i++)
      this->data_[i] *= alpha;
    
    return *this;
  }
  

  //! adds terms to a row of the matrix
  /*!
    \param[in] i row where values are added
    \param[in] nb_interac number of values to add to the matrix
    \param[in] col_interac column numbers associated with values
    \param[in] val_itnterac values to add
    Equivalent Matlab operation A(i, col_interac) = A(i, col_interac) + val_interac
   */  
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T,Prop,Storage,Allocator>::
  AddInteractionRow(int i, int nb_interac, const IVect& col_interac,
                    const Vector<T>& val_interac, bool sorted)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("MatrixSym_BlockDiagonal::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
#endif
    
    int n1 = this->num_blok_[i];
    int row = this->rev_[i];
    int shift = offset_[n1] + row*(ptr_[n1+1]-ptr_[n1]) - row*(row+1)/2;
    for (int j1 = 0; j1 < nb_interac; j1++)
      {
	int j = col_interac(j1);
#ifdef SELDON_CHECK_BOUNDS
	if (j < 0 || j >= this->m_)
	  throw WrongCol("MatrixSym_BlockDiagonal::operator()",
			 string("Index should be in [0, ")
			 + to_str(this->m_-1) + "], but is equal to "
			 + to_str(j) + ".");
#endif
	int n2 = this->num_blok_[j];
	if (n2 == n1)
	  {
	    int col = this->rev_[j];
	    if (row <= col)
	      this->data_[shift + col] += val_interac(j1);
	    
	    // DISP(this->data_[shift + col]);
	  }
      }
  }
  
  
  //! writes matrix in Matlab sparse format
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>
  ::WriteText(const string& file_name) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("MatrixSym_BlockDiagonal::WriteText(string FileName)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! writes matrix in Matlab sparse format
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {
    int nb = 0, irow, icol;
    for (int i = 0; i < nblok_; i++)
      {
	int size_block = ptr_[i+1] - ptr_[i];
	for (int j = 0; j < size_block; j++)
	  {
	    irow = ind_[j + ptr_[i] ];
	    for (int k = j; k < size_block; k++)
	      {
		icol = ind_[k + ptr_[i] ];
		if (j != k)
		  FileStream << icol+1 << ' ' << irow+1 << ' ' << this->data_[nb] << '\n';
		
		FileStream << irow+1 << ' ' << icol+1 << ' ' << this->data_[nb++] << '\n';
	      }
	  }
      }
  }


  /*************
   * Functions *
   *************/

  
  //! adds alpha*A to the current matrix
  template<class T0, class T1, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		 Matrix<T1, General, BlockDiagRow, Allocator2>& B)
  {
    // treating only the case where the profile of A
    // is included in the profile of matrix B
    int* ptrA = A.GetRowPtr();
    int* indA = A.GetRowNumData();
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        // finding the block number in the current matrix (num)
        int num = B.GetBlockNumber(indA[ptrA[i]]);
        
        // then checking that all rows of the block of A
        // are belonging to the block num of B
        for (int j = ptrA[i]+1; j < ptrA[i+1]; j++)
          if (num != B.GetBlockNumber(indA[j]))
            {
              cout << "The matrix A is not included in the current matrix" << endl;
              cout << "The matrix A cannot be added (not implemented) " << endl;
              abort();
            }
      }
    
    // adding contents of A with AddInteraction
    int size_max = A.GetSizeMaxBlock();
    Vector<int> num(size_max);
    Vector<T1> val(size_max);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        int size_block = A.GetBlockSize(i);
        for (int j = 0; j < size_block; j++)
          num(j) = A.Index(i, j);
        
        for (int j = 0; j < size_block; j++)
          {
            int irow = A.Index(i, j);
            for (int k = 0; k < size_block; k++)
              val(k) = alpha*A.Value(i, j, k);
            
            B.AddInteractionRow(irow, size_block, num, val);
          }
      }
  }
  
  
  //! conversion from block-diagonal matrix to sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, General, BlockDiagRow, Alloc>& A,
                       Matrix<T, General, ArrayRowSparse, Alloc2>& B)
  {
    int n = A.GetM();
    B.Clear();
    B.Reallocate(n, n);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
	int nb_val = A.GetBlockSize(i);
	IVect col(nb_val);
	for (int j = 0; j < nb_val; j++)
	  col(j) = A.Index(i, j);
	
        // permut in case where col is not increasing
        IVect permut(nb_val); permut.Fill();
        Sort(col, permut);
        for (int j = 0; j < nb_val; j++)
	  {
            B.ReallocateRow(col(j), nb_val);
	    for (int k = 0; k < nb_val; k++)
	      {
                B.Index(col(j), k) = col(k);
                B.Value(col(j), k) = A.Value(i, permut(j), permut(k));
              }
          }        
      }
  }


  //! conversion from block-diagonal matrix to sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, General, BlockDiagRow, Alloc>& A,
                       Matrix<T, General, ArrayRowComplexSparse, Alloc2>& B)
  {
    int n = A.GetM();
    B.Clear();
    B.Reallocate(n, n);
    typename ClassComplexType<T>::Treal zero; SetComplexZero(zero);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
	int nb_val = A.GetBlockSize(i);
	IVect col(nb_val);
	for (int j = 0; j < nb_val; j++)
	  col(j) = A.Index(i, j);

        // permut in case where col is not increasing
        IVect permut(nb_val); permut.Fill();
        Sort(col, permut);	
        T val;
        for (int j = 0; j < nb_val; j++)
	  {
            int nb_val_real = 0, nb_val_imag = 0;
            for (int k = 0; k < nb_val; k++)
              {
                val = A.Value(i, permut(j), permut(k));
                if (real(val) != zero)
                  nb_val_real++;

                if (imag(val) != zero)
                  nb_val_imag++;
              }
            
            B.ReallocateRealRow(col(j), nb_val_real);
            B.ReallocateImagRow(col(j), nb_val_imag);
            nb_val_real = 0; nb_val_imag = 0;
	    for (int k = 0; k < nb_val; k++)
	      {
                val = A.Value(i, permut(j), permut(k));
                if (real(val) != zero)
                  {
                    B.IndexReal(col(j), nb_val_real) = col(k);
                    B.ValueReal(col(j), nb_val_real) = real(val);
                    nb_val_real++;
                  }
                
                if (imag(A.Value(i, permut(j), permut(k))) != zero)
                  {
                    B.IndexImag(col(j), nb_val_imag) = col(k);
                    B.ValueImag(col(j), nb_val_imag) = imag(val);
                    nb_val_imag++;
                  }
              }
          }
      }
  }
  
  
  //! conversion from block-diagonal matrix to dense matrix
  template<class T, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, General, BlockDiagRow, Allocator1>& A,
		  Matrix<T, General, RowMajor, Allocator2>& B)
  {
    int n = A.GetM();
    B.Reallocate(n, n);
    T zero; SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        int size_block = A.GetBlockSize(i);
        for (int j = 0; j < size_block; j++)
          for (int k = 0; k < size_block; k++)
            B(A.Index(i, j), A.Index(i, k)) = A.Value(i, j, k);
      }
  }
  
  
  //! adds alpha*A to the current matrix
  template<class T0, class T1, class Allocator1, class Allocator2>
  void AddMatrix(const T0& alpha, const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		 Matrix<T1, Symmetric, BlockDiagRowSym, Allocator2>& B)
  {
    // treating only the case where the profile of A
    // is included in the profile of matrix B
    int* ptrA = A.GetRowPtr();
    int* indA = A.GetRowNumData();
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        // finding the block number in the current matrix (num)
        int num = B.GetBlockNumber(indA[ptrA[i]]);
        
        // then checking that all rows of the block of A
        // are belonging to the block num of B
        for (int j = ptrA[i]+1; j < ptrA[i+1]; j++)
          if (num != B.GetBlockNumber(indA[j]))
            {
              cout << "The matrix A is not included in the current matrix" << endl;
              cout << "The matrix A cannot be added (not implemented) " << endl;
              abort();
            }
      }
    
    // adding contents of A with AddInteraction
    int size_max = A.GetSizeMaxBlock();
    Vector<int> num(size_max);
    Vector<T1> val(size_max);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        int size_block = A.GetBlockSize(i);
        for (int j = 0; j < size_block; j++)
          num(j) = A.Index(i, j);
        
        for (int j = 0; j < size_block; j++)
          {
            int irow = A.Index(i, j);
            for (int k = 0; k < size_block; k++)
              val(k) = alpha*A.Value(i, j, k);
            
            B.AddInteractionRow(irow, size_block, num, val);
          }
      }
  }


  /************************************
   * Conversions with sparse matrices *
   ************************************/
  
  
  //! convert a block-diagonal matrix A to a sparse matrix B
  template<class T, class Allocator1, class Allocator2>
  void ConvertToSparse(const Matrix<T, Symmetric, BlockDiagRowSym, Allocator1>& A,
                       Matrix<T, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int n = A.GetM();
    B.Clear();
    B.Reallocate(n, n);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
	int nb_val = A.GetBlockSize(i);
	IVect col(nb_val);
	for (int j = 0; j < nb_val; j++)
	  col(j) = A.Index(i, j);
	
        // permut in case where col is not increasing
        IVect permut(nb_val); permut.Fill();
        Sort(col, permut);
        for (int j = 0; j < nb_val; j++)
	  {
            B.ReallocateRow(col(j), nb_val-j);
	    for (int k = j; k < nb_val; k++)
	      {
                B.Index(col(j), k-j) = col(k);
                B.Value(col(j), k-j) = A.Value(i, permut(j), permut(k));
              }
          }
      }
  }
  

  //! conversion from block-diagonal matrix to sparse matrix
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       Matrix<T, Symmetric, ArrayRowSymComplexSparse, Alloc2>& B)
  {
    int n = A.GetM();
    B.Clear();
    B.Reallocate(n, n);
    typename ClassComplexType<T>::Treal zero; SetComplexZero(zero);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
	int nb_val = A.GetBlockSize(i);
	IVect col(nb_val);
	for (int j = 0; j < nb_val; j++)
	  col(j) = A.Index(i, j);
	
        // permut in case where col is not increasing
        IVect permut(nb_val); permut.Fill();
        Sort(col, permut);
        T val;
        for (int j = 0; j < nb_val; j++)
	  {
            int nb_val_real = 0, nb_val_imag = 0;
            for (int k = j; k < nb_val; k++)
              {
                val = A.Value(i, permut(j), permut(k));
                if (real(val) != zero)
                  nb_val_real++;

                if (imag(val) != zero)
                  nb_val_imag++;
              }
            
            B.ReallocateRealRow(col(j), nb_val_real);
            B.ReallocateImagRow(col(j), nb_val_imag);
            nb_val_real = 0; nb_val_imag = 0;
	    for (int k = j; k < nb_val; k++)
	      {
                val = A.Value(i, permut(j), permut(k));
                if (real(val) != zero)
                  {
                    B.IndexReal(col(j), nb_val_real) = col(k);
                    B.ValueReal(col(j), nb_val_real) = real(val);
                    nb_val_real++;
                  }
                
                if (imag(val) != zero)
                  {
                    B.IndexImag(col(j), nb_val_imag) = col(k);
                    B.ValueImag(col(j), nb_val_imag) = imag(val);
                    nb_val_imag++;
                  }
              }
          }
      }
  }

  
  //! conversion from block-diagonal matrix to dense matrix
  template<class T, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Symmetric, BlockDiagRowSym, Allocator1>& A,
		  Matrix<T, Symmetric, RowSymPacked, Allocator2>& B)
  {
    int n = A.GetM();
    B.Reallocate(n, n);
    T zero; SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        int size_block = A.GetBlockSize(i);
        for (int j = 0; j < size_block; j++)
          for (int k = j; k < size_block; k++)
            B(A.Index(i, j), A.Index(i, k)) = A.Value(i, j, k);
      }
  }

  
  void MarkUsedColumnsSparseMatrix(const Vector<int>& Ptr, const Vector<int>& Ind,
                                   int irow, int m, int m1,
                                   Vector<bool>& UsedRow, IVect& row_found)
  {
    for (int j = Ptr(irow); j < Ptr(irow+1); j++)
      {
        int icol = Ind(j);
        if ((icol >= 0)&&(icol < m))
          if (!UsedRow(icol))
            {
              // we found a new row to add to the block
              row_found.PushBack(icol);
              UsedRow(icol) = true;
            }
      }
  }

  
  void GetBlockProfile(const Vector<int>& Ptr, const Vector<int>& Ind,
                       Vector<IVect>& pattern, int m1, int m2)
  {
    if (m2 <= m1)
      return;
    
    int m = m2-m1;
    int nblock = 0, nrow = 0;
    int irow = 0, first_row = 0;
    Vector<bool> UsedRow(m); UsedRow.Fill(false);
    pattern.Reallocate(m);
    while (nrow < m)
      {
	// we look for all other rows, until we find an isolated block
	IVect new_row(1); new_row(0) = first_row;
	UsedRow(first_row) = true;
	IVect all_row;
	while (new_row.GetM() > 0)
	  {
	    // in all_row, we store all the row numbers of the given block
	    all_row.PushBack(new_row);
	    
	    // now finding new rows for the next iteration
	    IVect row_found;
	    for (int i = 0; i < new_row.GetM(); i++)
	      {
		irow = new_row(i) + m1;
                // loop on all entries of this row and marking them used
                MarkUsedColumnsSparseMatrix(Ptr, Ind, irow, m, m1, UsedRow, row_found);
              }
	    
	    new_row = row_found;
	  }
	
	// specifying the new block
	nrow += all_row.GetM();
	pattern(nblock) = all_row;
	nblock++;
	
	// finding an available row to find a new block
	if (nrow < m)
	  for (int k = first_row+1; k < m; k++)
	    if (!UsedRow(k))
	      {
		first_row = k;
		break;
	      }
	
      }

    pattern.Resize(nblock);
  }
  
  
  //! conversion of a sparse matrix into a block-diagonal matrix
  /*!
    \param[in] A sparse matrix
    \param[out] B block-diagonal matrix
    \param[in] m1 first subscript
    \param[in] m2 second subscript
    We convert a sub-matrix of A into block-diagonal matrix. 
    The submatrix has rows and columns numbers comprised between m1 and m2.
    Call ConverToBlockDiagonal(A, B, 0, A.GetM()); if you want to convert all the matrix
   */
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {
    Vector<IVect> pattern;

    {
      Vector<int> Ptr, Ind;
      GetSymmetricPattern(A, Ptr, Ind);
      GetBlockProfile(Ptr, Ind, pattern, m1, m2);
    }
    
    // initialization of the block-diagonal matrix
    B.Clear();
    B.SetPattern(pattern);
    pattern.Clear();
    
    FillBlockDiagonal(A, B, m1, m2);
  }


  //! conversion of a sparse matrix into a block-diagonal matrix
  /*!
    \param[in] A sparse matrix
    \param[out] B block-diagonal matrix
  */
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! Filling block-diagonal matrix B with values of A
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			 Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    FillBlockDiagonal(A, B, 0, A.GetM());
  }


  //! Filling block-diagonal matrix B with values of A
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			 Matrix<T1, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {  
    int m = m2 - m1;
    // filling the block-diagonal matrix with values of A
    for (int i = 0; i < m; i++)
      {
	int irow = i + m1;
	int nb_val  = A.GetRowSize(irow);
	IVect col_num(nb_val);
	Vector<T1> value(nb_val);
	nb_val = 0;
	for (int j = 0; j < A.GetRowSize(irow); j++)
	  {
	    int icol = A.Index(irow, j) - m1;
	    if ((icol >= 0)&&(icol < m))
	      {
		col_num(nb_val) = icol;
		value(nb_val) = A.Value(irow, j);
		nb_val++;
	      }
	  }
	
	B.AddInteractionRow(i, nb_val, col_num, value);
      }  
  }
  
  
  //! Filling block-diagonal matrix B with values of A
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {  
    int m = m2 - m1;
    // filling the block-diagonal matrix with values of A
    for (int i = 0; i < m; i++)
      {
	int irow = i + m1;
	int nb_r  = A.GetRealRowSize(irow);
	int nb_i  = A.GetImagRowSize(irow);
	IVect col_num(nb_r+nb_i);
	Vector<complex<T1> > value(nb_r+nb_i);
	int nb_val = 0;
	for (int j = 0; j < nb_r; j++)
	  {
	    int icol = A.IndexReal(irow, j) - m1;
	    if ((icol >= 0)&&(icol < m))
	      {
		col_num(nb_val) = icol;
		value(nb_val) = complex<T1>(A.ValueReal(irow, j), 0);
		nb_val++;
	      }
	  }
	
	for (int j = 0; j < nb_i; j++)
	  {
	    int icol = A.IndexImag(irow, j) - m1;
	    if ((icol >= 0)&&(icol < m))
	      {
		col_num(nb_val) = icol;
		value(nb_val) = complex<T1>(0, A.ValueImag(irow, j));
		nb_val++;
	      }
	  }
	
	B.AddInteractionRow(i, nb_val, col_num, value);
      }  
  }
  
  
  //! Filling block-diagonal matrix B with values of A
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {  
    int m = m2 - m1;
    // filling the block-diagonal matrix with values of A
    for (int i = 0; i < m; i++)
      {
	int irow = i + m1;
	int nb_r  = A.GetRealRowSize(irow);
	int nb_i  = A.GetImagRowSize(irow);
	IVect col_num(nb_r+nb_i);
	Vector<complex<T1> > value(nb_r+nb_i);
	int nb_val = 0;
	for (int j = 0; j < nb_r; j++)
	  {
	    int icol = A.IndexReal(irow, j) - m1;
	    if ((icol >= 0)&&(icol < m))
	      {
		col_num(nb_val) = icol;
		value(nb_val) = complex<T1>(A.ValueReal(irow, j), 0);
		nb_val++;
	      }
	  }
	
	for (int j = 0; j < nb_i; j++)
	  {
	    int icol = A.IndexImag(irow, j) - m1;
	    if ((icol >= 0)&&(icol < m))
	      {
		col_num(nb_val) = icol;
		value(nb_val) = complex<T1>(0, A.ValueImag(irow, j));
		nb_val++;
	      }
	  }
	
	B.AddInteractionRow(i, nb_val, col_num, value);
      }  
  }
  

  //! forbidden case
  template<class T0, class Prop0, class Allocator0,
	   class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<complex<T0>, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
			 Matrix<T0, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {
    throw WrongArgument("FillBlockDiagonal(A, B, m1, m2)", "B must be a complex matrix");
  }


  //! forbidden case
  template<class T0, class Prop0, class Allocator0,
	   class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<complex<T0>, Prop0, ArrayRowComplexSparse, Allocator0>& A,
			 Matrix<T0, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {
    throw WrongArgument("FillBlockDiagonal(A, B, m1, m2)", "B must be a complex matrix");
  }

  
  //! conversion of a sparse matrix into a block-diagonal matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! conversion of a sparse matrix into a block-diagonal matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! conversion of a sparse matrix into a block-diagonal matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! conversion of a sparse matrix into a block-diagonal matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
		  Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }

  
  //! conversion of a block-diagonal matrix into sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    ConvertToSparse(A, B);
  }


  //! conversion of a block-diagonal matrix into sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    ConvertToSparse(A, B);
  }


  //! conversion of a block-diagonal matrix into sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B)
  {
    ConvertToSparse(A, B);
  }


  //! conversion of a block-diagonal matrix into sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B)
  {
    ConvertToSparse(A, B);
  }

  
  /**************************
   * Matrix-vector products *
   **************************/


  //! replaces X by A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X)
  {
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    int* ind_ = A.GetInd();
    T1* data_ = A.GetData();
    T2* x_ptr = X.GetData();
    Vector<T2> bloc(A.GetSizeMaxBlock());
    T2 zero; SetComplexZero(zero);
    for (i = 0; i < nblock; i++)
      {
	nb = offset_[i];
        for (k = 0; k < ptr_[i+1]-ptr_[i]; k++)
          bloc(k) = zero;
        
	for (j = ptr_[i]; j < ptr_[i+1]; j++)
          for (k = ptr_[i]; k < ptr_[i+1]; k++)
            bloc(j-ptr_[i]) += data_[nb++]*x_ptr[ind_[k]];
        
        for (k = ptr_[i]; k < ptr_[i+1]; k++)
          x_ptr[ind_[k]] = bloc(k-ptr_[i]);
      }
  }
  

  //! replaces X by A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X)
  {
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k, row;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    int* ind_ = A.GetInd();
    T1* data_ = A.GetData();
    T2* x_ptr = X.GetData();
    Vector<T2> bloc(A.GetSizeMaxBlock());
    T2 val;
    T2 zero; SetComplexZero(zero);
    for (i = 0; i < nblock; i++)
      {
	nb = offset_[i];
        for (k = 0; k < ptr_[i+1]-ptr_[i]; k++)
          bloc(k) = zero;
        
        for (j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            row = ind_[j];
            // diagonal term
            val = data_[nb++]*x_ptr[row];
            // extra-diagonal contributions
            for (k = (j+1); k < ptr_[i+1]; k++)
              {
                val += data_[nb]*x_ptr[ind_[k]];
                bloc(k-ptr_[i]) += data_[nb]*x_ptr[row];
                nb++;
              }
            
            bloc(j-ptr_[i]) += val;
          }
        
        for (k = ptr_[i]; k < ptr_[i+1]; k++)
          x_ptr[ind_[k]] = bloc(k-ptr_[i]);
      }
  }

  
  //! Computes Y = A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y)
  {
    Copy(X, Y);
    Mlt(A, Y);
  }


  //! Computes Y = A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y)
  {
    Copy(X, Y);
    Mlt(A, Y);
  }


  //! Computes Y = A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const SeldonTranspose& trans,
		 const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y)
  {
    if (trans.NoTrans())
      MltVector(A, X, Y);
    else if (trans.Trans())
      {
	int nblock = A.GetNbBlocks(), nb = 0, i, j, k, row;
	int* ptr_ = A.GetPtr();
	int* offset_ = A.GetOffset();
	int* ind_ = A.GetInd();
	T1* data_ = A.GetData();
	T2* x_ptr = X.GetData();
	T3* y_ptr = Y.GetData();

	Y.Zero();
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		for (k = ptr_[i]; k < ptr_[i+1]; k++)
		  y_ptr[ind_[k]] += data_[nb++]*x_ptr[row];
	      }
	  }
      }
    else
      {
	cout << "Not implemented" << endl;
	abort();
      }
  }


  //! Computes Y = A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2,
	   class T3, class Allocator3>
  void MltVector(const SeldonTranspose&,
		 const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
		 const Vector<T2, VectFull, Allocator2>& X,
		 Vector<T3, VectFull, Allocator3>& Y)
  {
    Copy(X, Y);
    Mlt(A, Y);
  }

  
  //! computes C = A B
  template<class T0, class Allocator0, class T1, class Allocator1,
           class T2, class Allocator2, class T3, class T4>
  void MltAddMatrix(const T3& alpha,
		    const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& B,
		    const T4& beta,
		    Matrix<T2, Symmetric, BlockDiagRowSym, Allocator2>& C)
  {
    // we assume that A B and C are already constructed (with SetPattern)
    // and have the same pattern
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    T0* dataA = A.GetData();
    T1* dataB = B.GetData();
    T2* dataC = C.GetData();
    Matrix<T0, Symmetric, RowSymPacked, Allocator0> Aloc;
    Matrix<T1, Symmetric, RowSymPacked, Allocator1> Bloc;
    Matrix<T2, Symmetric, RowSymPacked, Allocator2> Cloc;
    for (int i = 0; i < A.GetNbBlocks(); i++)
      {
        int size_mat = ptr_[i+1] - ptr_[i];
        int nb = offset_[i];
        Aloc.SetData(size_mat, size_mat, &dataA[nb]);
        Bloc.SetData(size_mat, size_mat, &dataB[nb]);
        Cloc.SetData(size_mat, size_mat, &dataC[nb]);
        
        MltAdd(alpha, Aloc, Bloc, beta, Cloc);
        
        Aloc.Nullify();
        Bloc.Nullify();
        Cloc.Nullify();
      }
  }


  //! Matrix vector product with block-diagonal matrix
  /*!
    Performs operation : B = beta*B + alpha*A*X
   */
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T4 zero; SetComplexZero(zero);
    T4 one; SetComplexOne(one);
    
    if (beta == zero)
      B.Fill(zero);
    else if (beta != one)
      Mlt(beta, B);
    
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k, row;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    int* ind_ = A.GetInd();
    T1* data_ = A.GetData();
    T2* x_ptr = X.GetData();
    T4* b_ptr = B.GetData(); T4 val;     
    for (i = 0; i < nblock; i++)
      {
	nb = offset_[i];
	for (j = ptr_[i]; j < ptr_[i+1]; j++)
	  {
	    row = ind_[j];
	    val = zero;
	    for (k = ptr_[i]; k < ptr_[i+1]; k++)
	      val += data_[nb++]*x_ptr[ind_[k]];
	    
	    b_ptr[row] += alpha*val;
	  }
      }
  }


  //! Matrix vector product with block-diagonal matrix
  /*!
    Performs operation : B = beta*B + alpha*A*X
   */
  template<class T0, class T1, class Prop, class Allocator1,
	   class T2, class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, Prop, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T4 zero; SetComplexZero(zero);
    T4 one; SetComplexOne(one);

    if (beta == zero)
      B.Fill(zero);
    else if (beta != one)
      Mlt(beta, B);
    
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k, row;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    int* ind_ = A.GetInd();
    T1* data_ = A.GetData();
    T2* x_ptr = X.GetData();
    T4* b_ptr = B.GetData();
    T4 val; 
    if (alpha == one)
      {
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		// diagonal term
                val = data_[nb++]*x_ptr[row];
		// extra-diagonal contributions
		for (k = (j+1); k < ptr_[i+1]; k++)
		  {
		    val += data_[nb]*x_ptr[ind_[k]];
		    b_ptr[ind_[k]] += data_[nb]*x_ptr[row];
                    nb++;
		  }
		
		b_ptr[row] += val;
	      }
	  }
      }
    else
      {
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		// diagonal term
		val = data_[nb++]*x_ptr[row];
		// extra-diagonal contributions
		for (k = (j+1); k < ptr_[i+1]; k++)
		  {
		    val += data_[nb]*x_ptr[ind_[k]];
		    b_ptr[ind_[k]] += alpha*data_[nb]*x_ptr[row];
		    nb++;
		  }
		
		b_ptr[row] += alpha*val;
	      }
	  }
      }
  }
  
    
  //! Matrix vector product with block-diagonal matrix
  /*!
    Performs operation : B = beta*B + alpha*A*X
   */
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    if (trans.NoTrans())
      {
	MltAddVector(alpha, A, X, beta, B);
	return;
      }
    
    T4 zero; SetComplexZero(zero);
    T4 one; SetComplexOne(one);

    if (beta == zero)
      B.Fill(zero);
    else if (beta != one)
      Mlt(beta, B);
    
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k, row;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    int* ind_ = A.GetInd();
    T1* data_ = A.GetData();
    T2* x_ptr = X.GetData();
    T4* b_ptr = B.GetData();
    if (trans.Trans())
      {
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		for (k = ptr_[i]; k < ptr_[i+1]; k++)
		  b_ptr[ind_[k]] += alpha*data_[nb++]*x_ptr[row];
	      }
	  }
      }
    else
      {
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		for (k = ptr_[i]; k < ptr_[i+1]; k++)
		  b_ptr[ind_[k]] += alpha*conjugate(data_[nb++])*x_ptr[row];
	      }
	  }
      }
  }
  

  //! Matrix vector product with block-diagonal matrix
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const SeldonTranspose& trans,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    if (!trans.ConjTrans())
      {
	MltAddVector(alpha, A, X, beta, B);
	return;
      }
    
    T4 zero; SetComplexZero(zero);
    T4 one; SetComplexOne(one);

    if (beta == zero)
      B.Fill(zero);
    else if (beta != one)
      Mlt(beta, B);
    
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k, row;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    int* ind_ = A.GetInd();
    T1* data_ = A.GetData();
    T2* x_ptr = X.GetData();
    T4* b_ptr = B.GetData();
    T4 val; 
    if (alpha == one)
      {
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		// diagonal term
                val = conjugate(data_[nb++])*x_ptr[row];
		// extra-diagonal contributions
		for (k = (j+1); k < ptr_[i+1]; k++)
		  {
		    val += conjugate(data_[nb])*x_ptr[ind_[k]];
		    b_ptr[ind_[k]] += conjugate(data_[nb])*x_ptr[row];
                    nb++;
		  }
		
		b_ptr[row] += val;
	      }
	  }
      }
    else
      {
	for (i = 0; i < nblock; i++)
	  {
	    nb = offset_[i];
	    for (j = ptr_[i]; j < ptr_[i+1]; j++)
	      {
		row = ind_[j];
		// diagonal term
		val = conjugate(data_[nb++])*x_ptr[row];
		// extra-diagonal contributions
		for (k = (j+1); k < ptr_[i+1]; k++)
		  {
		    val += conjugate(data_[nb])*x_ptr[ind_[k]];
		    b_ptr[ind_[k]] += alpha*conjugate(data_[nb])*x_ptr[row];
		    nb++;
		  }
		
		b_ptr[row] += alpha*val;
	      }
	  }
      }
  }
  
  
  //! A is replaced by its inverse
  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, BlockDiagRow, Allocator>& A)
  {
    T zero; SetComplexZero(zero);
    T one; SetComplexOne(one);

    int nblock = A.GetNbBlocks(), nb = 0, i, j, k;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
    Matrix<T, General, RowMajor, Allocator> Asub;
    for (i = 0; i < nblock; i++)
      {
	int n = ptr_[i+1]-ptr_[i];
	Asub.Reallocate(n,n);
	nb = offset_[i];
	for (j = 0; j < n; j++)
	  for (k = 0; k < n; k++)
	    Asub(j,k) = data_[nb++];

        if (n == 1)
          {
            // pseudo-inverse, not inverting the block is equal to 0
            if (Asub(0, 0) != zero)
              Asub(0, 0) = one/Asub(0, 0);
          }
        else
          GetInverse(Asub);
        
	nb = offset_[i];
	for (j = 0; j < n; j++)
	  for (k = 0; k < n; k++)
	    data_[nb++] = Asub(j,k); 
      }
  }
  
  
  //! A is replaced by its inverse
  template<class T, class Allocator>
  void GetInverse(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
    T one; SetComplexOne(one);
    T zero; SetComplexZero(zero);
    Matrix<T, Symmetric, RowSymPacked, Allocator> Asub;
    for (i = 0; i < nblock; i++)
      {
	int n = ptr_[i+1]-ptr_[i];
	Asub.Reallocate(n,n);
	nb = offset_[i];
	for (j = 0; j < n; j++)
	  for (k = j; k < n; k++)
	    Asub(j,k) = data_[nb++];
        
        if (n == 1)
          {
            // pseudo-inverse, not inverting the block is equal to 0
            if (Asub(0, 0) != zero)
              Asub(0, 0) = one/Asub(0, 0);
          }
        else
          GetInverse(Asub);
        
	nb = offset_[i];
	for (j = 0; j < n; j++)
	  for (k = j; k < n; k++)
	    data_[nb++] = Asub(j,k); 
      }
  }

  
  //! A is replaced by its Cholesky factorisation
  template<class T, class Allocator>
  void GetCholesky(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    int nblock = A.GetNbBlocks(), nb = 0, i, j, k;
    int* ptr_ = A.GetPtr();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
    Matrix<T, Symmetric, RowSymPacked, Allocator> Asub;
    for (i = 0; i < nblock; i++)
      {
	int n = ptr_[i+1]-ptr_[i];
	Asub.Reallocate(n, n);
	nb = offset_[i];
	for (j = 0; j < n; j++)
	  for (k = j; k < n; k++)
	    Asub(j, k) = data_[nb++];
        
	GetCholesky(Asub);
        
	nb = offset_[i];
	for (j = 0; j < n; j++)
	  for (k = j; k < n; k++)
	    data_[nb++] = Asub(j,k); 
      }
  }
  
  
  //! solving L x = b or L^T x = b
  template<class T, class Allocator, class Allocator2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                     Vector<T, VectFull, Allocator2>& x)
  {
    int nblock = A.GetNbBlocks(), nb = 0;
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
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
                val = x(irow);
                for (int j = ptr_[k+1]-1; j >= i+1; j--)
                  val -= data_[nb--]*x(ind_[j]);
                
                x(irow) = val/data_[nb--];
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
                x(irow) /= data_[nb++];
                val = x(irow);
                for (int j = i+1; j < ptr_[k+1]; j++)
                  x(ind_[j]) -= data_[nb++]*val;                
              }
          }
      }
  }
  
  
  //! multiplying L x or L^T x
  template<class T, class Allocator, class Allocator2>
  void MltCholesky(const SeldonTranspose& trans,
                   const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                   Vector<T, VectFull, Allocator2>& x)
  {
    int nblock = A.GetNbBlocks(), nb = 0;
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    int* offset_ = A.GetOffset();
    T* data_ = A.GetData();
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
                val = x(irow)*data_[nb++];
                for (int j = i+1; j < ptr_[k+1]; j++)
                  val += data_[nb++]*x(ind_[j]);
                
                x(irow) = val; 
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
                val = x(irow);
                for (int j = ptr_[k+1]-1; j >= i+1; j--)
                  x(ind_[j]) += data_[nb--]*val;
                  
                x(irow) *= data_[nb--];
              }
          }
      }
  }
  
  
  //! A is multiplied by a scalar
  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, General, BlockDiagRow, Allocator>& A)
  {
    A *= alpha;
  }
  

  //! A is multiplied by a scalar
  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    A *= alpha;
  }
  
}

#define SELDON_FILE_MATRIX_BLOCK_DIAGONAL_CXX
#endif
