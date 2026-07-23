#ifndef SELDON_FILE_MATRIX_BLOCK_DIAGONAL_CXX

#include "OptMatrixBlockDiagonal.hxx"

namespace Seldon
{
  //! constructor for an empty matrix
  template<class T, class Prop, class Storage, class Allocator>
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::MatrixSym_BlockDiagonal()
  {
    m_ = 0;
    n_ = 0;
    size_max_blok_ = 0;
  }
  
  
  //! clears matrix
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Clear()
  {
    all_mat_.Clear(); blok_num_.Clear(); size_block_.Clear();
    mat99_.Clear(); mat88_.Clear(); mat77_.Clear(); mat66_.Clear(); mat55_.Clear();
    mat44_.Clear(); mat33_.Clear(); mat22_.Clear(); mat1010_.Clear(); mat1111_.Clear();
    diag_value_.Clear(); ptr_.Clear(); rev_num_.Clear(); row_num_.Clear();
    m_ = 0;
    n_ = 0;
    size_max_blok_ = 0;
    
    mat1212_.Clear(); mat1313_.Clear(); mat1414_.Clear(); //mat1515_.Clear();
    //mat1616_.Clear(); mat1717_.Clear(); mat1818_.Clear(); mat1919_.Clear();
    //mat2020_.Clear();
  }
  
  
  //! returns size of row i
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowSize(int i) const
  {
    return size_block_(i);
  }
  
  
  //! returns the maximal size among all the blocks
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetSizeMaxBlock() const
  {
    return size_max_blok_;
  }
  
  
  //! returns the block number associated with row i
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetBlockNumber(int i) const
  {
    int n = this->blok_num_(i);
    switch (size_block_(i))
      {
      case 1 :
        return n;
      case 2 :
        return n + offset_blok_(1);
      case 3 :
        return n + offset_blok_(2);
      case 4 :
        return n + offset_blok_(3);
      case 5 :
        return n + offset_blok_(4);
      case 6 :
        return n + offset_blok_(5);
      case 7 :
        return n + offset_blok_(6);
      case 8 :
        return n + offset_blok_(7);
      case 9 :
        return n + offset_blok_(8);
      case 10 :
        return n + offset_blok_(9);
      case 11 :
        return n + offset_blok_(10);
      case 12 :
        return n + offset_blok_(11);
      case 13 :
        return n + offset_blok_(12);
      case 14 :
        return n + offset_blok_(13);
      default :
        return n + offset_blok_(max_size_small-1);
      }
  }
  

  //! specifies all the blocks of the matrix
  /*!
    \param[in] NumDof_Blocks an array containing the row numbers for each block of the matrix
    this methods clears the matrix and affects a new structure, while setting to 0 all the entries
    of the matrix
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class Allocator1>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks)
  {
    Vector<Vector<bool> > ignore_row;
    IVect permut;
    SetPattern(NumDof_Blocks, ignore_row, permut);
  }
  
  
  //! specifies all the blocks of the matrix
  /*!
    \param[in] NumDof_Blocks an array containing the row numbers for each block of the matrix
    this methods clears the matrix and affects a new structure, while setting to 0 all the entries
    of the matrix
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class Allocator1>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  SetPattern(const Vector<IVect, VectFull, Allocator1>& NumDof_Blocks,
             const Vector<Vector<bool> >& ignore_row, IVect& permut)
  {
    // previous data is cleared
    Clear();
    
    // allocating arrays related to blocks
    int nblock = NumDof_Blocks.GetM();
    permut.Reallocate(nblock);
    ptr_.Reallocate(nblock+1); ptr_.Fill(0);
    
    int size_block;
    // counting blocks of size 1, 2x2, 3x3, etc
    IVect size_dd(max_size_small); size_dd.Fill(0);
    int nb_big_block = 0;
    int nb_row = 0; m_ = 0;
    for (int i = 0; i < NumDof_Blocks.GetM(); i++)
      {
	size_block = NumDof_Blocks(i).GetM();
	if (size_block >= max_size_small)
	  nb_big_block++;
	else
	  size_dd(size_block)++;
	
        if (ignore_row.GetM() > 0)
          {
            for (int j = 0; j < ignore_row(i).GetM(); j++)
              if (!ignore_row(i)(j))
                m_++;
          }
        else
          m_ += size_block;
        
	nb_row += size_block;
      }
    
    n_ = m_;
    // nb_row and m_ may be different for distributed matrices
    // allocating arrays related to row numbers
    T zero; SetComplexZero(zero);

    size_block_.Reallocate(m_);
    row_num_.Reallocate(nb_row);
    blok_num_.Reallocate(m_);
    rev_num_.Reallocate(m_);
    rev_num_.Fill(-1);
    blok_num_.Fill(-1);
    diag_value_.Reallocate(size_dd(1));
    diag_value_.Fill(zero);
    mat22_.Reallocate(size_dd(2));
    mat33_.Reallocate(size_dd(3));
    mat44_.Reallocate(size_dd(4));
    mat55_.Reallocate(size_dd(5));
    mat66_.Reallocate(size_dd(6));
    mat77_.Reallocate(size_dd(7));
    mat88_.Reallocate(size_dd(8));
    mat99_.Reallocate(size_dd(9));
    mat1010_.Reallocate(size_dd(10));
    mat1111_.Reallocate(size_dd(11));
    mat1212_.Reallocate(size_dd(12));
    mat1313_.Reallocate(size_dd(13));
    mat1414_.Reallocate(size_dd(14));
    /*mat1515_.Reallocate(size_dd(15));
    mat1616_.Reallocate(size_dd(16));
    mat1717_.Reallocate(size_dd(17));
    mat1818_.Reallocate(size_dd(18));
    mat1919_.Reallocate(size_dd(19));
    mat2020_.Reallocate(size_dd(20)); */
    all_mat_.Reallocate(nb_big_block);
    
    size_max_blok_ = 0;
    for (int i = 0; i < NumDof_Blocks.GetM(); i++)
      size_max_blok_ = max(size_max_blok_, NumDof_Blocks(i).GetM());
    
    int nb = 0; nblock = 0; int num = 0;
    // small blocks
    for (size_block = 1; size_block < max_size_small; size_block++)
      {
	num = 0;
	for (int i = 0; i < NumDof_Blocks.GetM(); i++)
	  if (NumDof_Blocks(i).GetM() == size_block)
	    {
	      for (int j = 0; j < size_block; j++)
		{
		  row_num_(nb) = NumDof_Blocks(i)(j);
		  bool treat_row = true;
                  if (ignore_row.GetM() > 0)
                    if (ignore_row(i)(j))
                      treat_row = false;
                  
                  if (treat_row)
                    {
                      rev_num_(row_num_(nb)) = j;
                      blok_num_(row_num_(nb)) = num;
                      size_block_(row_num_(nb)) = size_block;
                    }
		  nb++;
		}
	      
              permut(nblock) = i;
	      ptr_(nblock+1) = ptr_(nblock) + size_block;
	      nblock++; num++;
	    }
      }
    
    // big blocks
    num = 0;
    for (int i = 0; i < NumDof_Blocks.GetM(); i++)
      if (NumDof_Blocks(i).GetM() >= max_size_small)
	{
	  size_block = NumDof_Blocks(i).GetM();
	  for (int j = 0; j < size_block; j++)
	    {
	      row_num_(nb) = NumDof_Blocks(i)(j);
              bool treat_row = true;
              if (ignore_row.GetM() > 0)
                if (ignore_row(i)(j))
                  treat_row = false;
              
              if (treat_row)
                {
                  rev_num_(row_num_(nb)) = j;
                  blok_num_(row_num_(nb)) = num;
                  size_block_(row_num_(nb)) = size_block;
                }
	      nb++;
	    }
	  
          permut(nblock) = i;
	  all_mat_(num).Reallocate(size_block, size_block);
	  all_mat_(num).Fill(zero);
	  ptr_(nblock+1) = ptr_(nblock) + size_block;
	  nblock++; num++;
	}
    
    // offset for each type of block
    offset_blok_.Reallocate(max_size_small);
    offset_blok_(0) = 0;
    for (int i = 0; i < max_size_small-1; i++)
      offset_blok_(i+1) = offset_blok_(i) + size_dd(i+1);
  }
  

  //! changes the size of the matrix (previous elements may be lost)
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    // nothing, the profile must be changed with SetPattern
  }
  

  //! changes the size of the matrix (previous elements are kept)  
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    // nothing, the profile must be changed with SetPattern
  }
  
  
  //! returns the number of non-zero entries of the matrix
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    int size = diag_value_.GetM() + TinyMatrix<T, Prop, 2, 2>::size_*mat22_.GetM()
      + TinyMatrix<T, Prop, 3, 3>::size_*mat33_.GetM() 
      + TinyMatrix<T, Prop, 4, 4>::size_*mat44_.GetM()
      + TinyMatrix<T, Prop, 5, 5>::size_*mat55_.GetM()
      + TinyMatrix<T, Prop, 6, 6>::size_*mat66_.GetM()
      + TinyMatrix<T, Prop, 7, 7>::size_*mat77_.GetM()
      + TinyMatrix<T, Prop, 8, 8>::size_*mat88_.GetM()
      + TinyMatrix<T, Prop, 9, 9>::size_*mat99_.GetM()
      + TinyMatrix<T, Prop, 10, 10>::size_*mat1010_.GetM() 
      + TinyMatrix<T, Prop, 11, 11>::size_*mat1111_.GetM()
      + TinyMatrix<T, Prop, 12, 12>::size_*mat1212_.GetM()
      + TinyMatrix<T, Prop, 13, 13>::size_*mat1313_.GetM()
      + TinyMatrix<T, Prop, 14, 14>::size_*mat1414_.GetM();
      /*+ TinyMatrix<T, Prop, 15, 15>::size_*mat1515_.GetM()
      + TinyMatrix<T, Prop, 16, 16>::size_*mat1616_.GetM()
      + TinyMatrix<T, Prop, 17, 17>::size_*mat1717_.GetM()
      + TinyMatrix<T, Prop, 18, 18>::size_*mat1818_.GetM()
      + TinyMatrix<T, Prop, 19, 19>::size_*mat1919_.GetM()
      + TinyMatrix<T, Prop, 20, 20>::size_*mat2020_.GetM(); */
    
    for (int i = 0; i < all_mat_.GetM(); i++)
      size += all_mat_(i).GetDataSize();
    
    return size;
  }

  
  //! returns the numbers of rows 
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetM() const
  {
    return m_;
  }

  
  //! returns the numbers of columns
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetN() const
  {
    return n_;
  }

  
  //! return the number of blocks
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetNbBlocks() const
  {
    return ptr_.GetM()-1;
  }

  
  //! returns size of a block
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetBlockSize(int i) const
  {
    return ptr_(i+1) - ptr_(i);
  }
    
  
  //! return the number of different row numbers
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetNbRowNumbers() const
  {
    return row_num_.GetM();
  }
  
  
  //! returns the array containing row numbers
  template<class T, class Prop, class Storage, class Allocator>
  int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowNumData() const
  {
    return row_num_.GetData();
  }
  
  
  //! returns the array containing offsets for each block
  template<class T, class Prop, class Storage, class Allocator>
  int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetRowPtr() const
  {
    return ptr_.GetData();
  }
  
  
  //! returns the array containing the local position of each row in the block
  template<class T, class Prop, class Storage, class Allocator>
  int* MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetReverseNum() const
  {
    return rev_num_.GetData();
  }

  
  //! returns row number of row j of block i
  template<class T, class Prop, class Storage, class Allocator>
  int MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Index(int i, int j) const
  {
    return row_num_(ptr_(i) + j);
  }
    
  
  //! returns value (j, k) of block i
  template<class T, class Prop, class Storage, class Allocator>
  const T& MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>
  ::Value(int i, int j, int k) const
  {
    switch (ptr_(i+1) - ptr_(i))
      {
      case 1 :
        return diag_value_(i);
      case 2 :
        return mat22_(i-offset_blok_(1))(j, k);
      case 3 :
        return mat33_(i-offset_blok_(2))(j, k);
      case 4 :
        return mat44_(i-offset_blok_(3))(j, k);
      case 5 :
        return mat55_(i-offset_blok_(4))(j, k);
      case 6 :
        return mat66_(i-offset_blok_(5))(j, k);
      case 7 :
        return mat77_(i-offset_blok_(6))(j, k);
      case 8 :
        return mat88_(i-offset_blok_(7))(j, k);
      case 9 :
        return mat99_(i-offset_blok_(8))(j, k);
      case 10 :
        return mat1010_(i-offset_blok_(9))(j, k);
      case 11 :
        return mat1111_(i-offset_blok_(10))(j, k);
      case 12 :
        return mat1212_(i-offset_blok_(11))(j, k);
      case 13 :
        return mat1313_(i-offset_blok_(12))(j, k);
      case 14 :
        return mat1414_(i-offset_blok_(13))(j, k);
        /* case 15 :
        return mat1515_(i-offset_blok_(14))(j, k);
      case 16 :
        return mat1616_(i-offset_blok_(15))(j, k);
      case 17 :
        return mat1717_(i-offset_blok_(16))(j, k);
      case 18 :
        return mat1818_(i-offset_blok_(17))(j, k);
      case 19 :
        return mat1919_(i-offset_blok_(18))(j, k);
      case 20 :
        return mat2020_(i-offset_blok_(19))(j, k);
        */
      }
    
    return all_mat_(i-offset_blok_(max_size_small-1))(j, k);
  }
  
  
  //! returns access to value (j, k) of block i
  template<class T, class Prop, class Storage, class Allocator>
  T& MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Value(int i, int j, int k)
  {
    switch (ptr_(i+1) - ptr_(i))
      {
      case 1 :
        return diag_value_(i);
      case 2 :
        return mat22_(i-offset_blok_(1))(j, k);
      case 3 :
        return mat33_(i-offset_blok_(2))(j, k);
      case 4 :
        return mat44_(i-offset_blok_(3))(j, k);
      case 5 :
        return mat55_(i-offset_blok_(4))(j, k);
      case 6 :
        return mat66_(i-offset_blok_(5))(j, k);
      case 7 :
        return mat77_(i-offset_blok_(6))(j, k);
      case 8 :
        return mat88_(i-offset_blok_(7))(j, k);
      case 9 :
        return mat99_(i-offset_blok_(8))(j, k);
      case 10 :
        return mat1010_(i-offset_blok_(9))(j, k);
      case 11 :
        return mat1111_(i-offset_blok_(10))(j, k);
      case 12 :
        return mat1212_(i-offset_blok_(11))(j, k);
      case 13 :
        return mat1313_(i-offset_blok_(12))(j, k);
      case 14 :
        return mat1414_(i-offset_blok_(13))(j, k);
        /* case 15 :
        return mat1515_(i-offset_blok_(14))(j, k);
      case 16 :
        return mat1616_(i-offset_blok_(15))(j, k);
      case 17 :
        return mat1717_(i-offset_blok_(16))(j, k);
      case 18 :
        return mat1818_(i-offset_blok_(17))(j, k);
      case 19 :
        return mat1919_(i-offset_blok_(18))(j, k);
      case 20 :
        return mat2020_(i-offset_blok_(19))(j, k);
        */
      }
    
    return all_mat_(i-offset_blok_(max_size_small-1))(j, k);
  }
  
  
  //! all the non-zero entries of the matrix are set to the same value x
  template<class T, class Prop, class Storage, class Allocator> template<class T0>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    diag_value_.Fill(x);
    for (int i = 0; i < mat22_.GetM(); i++)
      mat22_(i).Fill(x);
    
    for (int i = 0; i < mat33_.GetM(); i++)
      mat33_(i).Fill(x);
    
    for (int i = 0; i < mat44_.GetM(); i++)
      mat44_(i).Fill(x);
    
    for (int i = 0; i < mat55_.GetM(); i++)
      mat55_(i).Fill(x);
    
    for (int i = 0; i < mat66_.GetM(); i++)
      mat66_(i).Fill(x);

    for (int i = 0; i < mat77_.GetM(); i++)
      mat77_(i).Fill(x);
    
    for (int i = 0; i < mat88_.GetM(); i++)
      mat88_(i).Fill(x);

    for (int i = 0; i < mat99_.GetM(); i++)
      mat99_(i).Fill(x);

    for (int i = 0; i < mat1010_.GetM(); i++)
      mat1010_(i).Fill(x);

    for (int i = 0; i < mat1111_.GetM(); i++)
      mat1111_(i).Fill(x);

    for (int i = 0; i < mat1212_.GetM(); i++)
      mat1212_(i).Fill(x);

    for (int i = 0; i < mat1313_.GetM(); i++)
      mat1313_(i).Fill(x);

    for (int i = 0; i < mat1414_.GetM(); i++)
      mat1414_(i).Fill(x);

    /* for (int i = 0; i < mat1515_.GetM(); i++)
      mat1515_(i).Fill(x);

    for (int i = 0; i < mat1616_.GetM(); i++)
      mat1616_(i).Fill(x);

    for (int i = 0; i < mat1717_.GetM(); i++)
      mat1717_(i).Fill(x);

    for (int i = 0; i < mat1818_.GetM(); i++)
      mat1818_(i).Fill(x);

    for (int i = 0; i < mat1919_.GetM(); i++)
      mat1919_(i).Fill(x);

    for (int i = 0; i < mat2020_.GetM(); i++)
    mat2020_(i).Fill(x); */
    
    for (int i = 0; i < all_mat_.GetM(); i++)
      all_mat_(i).Fill(x);
  }

  
  //! values of non-zero entries of the matrix are set to random values
  template<class T, class Prop, class Storage, class Allocator>  
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::FillRand()
  {
    diag_value_.FillRand();
    for (int i = 0; i < mat22_.GetM(); i++)
      mat22_(i).FillRand();
    
    for (int i = 0; i < mat33_.GetM(); i++)
      mat33_(i).FillRand();
    
    for (int i = 0; i < mat44_.GetM(); i++)
      mat44_(i).FillRand();
    
    for (int i = 0; i < mat55_.GetM(); i++)
      mat55_(i).FillRand();
    
    for (int i = 0; i < mat66_.GetM(); i++)
      mat66_(i).FillRand();

    for (int i = 0; i < mat77_.GetM(); i++)
      mat77_(i).FillRand();
    
    for (int i = 0; i < mat88_.GetM(); i++)
      mat88_(i).FillRand();

    for (int i = 0; i < mat99_.GetM(); i++)
      mat99_(i).FillRand();

    for (int i = 0; i < mat1010_.GetM(); i++)
      mat1010_(i).FillRand();
    
    for (int i = 0; i < mat1111_.GetM(); i++)
      mat1111_(i).FillRand();
    
    for (int i = 0; i < mat1212_.GetM(); i++)
      mat1212_(i).FillRand();
    
    for (int i = 0; i < mat1313_.GetM(); i++)
      mat1313_(i).FillRand();
    
    for (int i = 0; i < mat1414_.GetM(); i++)
      mat1414_(i).FillRand();

    /* for (int i = 0; i < mat1515_.GetM(); i++)
      mat1515_(i).FillRand();
    
    for (int i = 0; i < mat1616_.GetM(); i++)
      mat1616_(i).FillRand();

    for (int i = 0; i < mat1717_.GetM(); i++)
      mat1717_(i).FillRand();

    for (int i = 0; i < mat1818_.GetM(); i++)
      mat1818_(i).FillRand();
    
    for (int i = 0; i < mat1919_.GetM(); i++)
      mat1919_(i).FillRand();

    for (int i = 0; i < mat2020_.GetM(); i++)
    mat2020_(i).FillRand(); */
    
    for (int i = 0; i < all_mat_.GetM(); i++)
      all_mat_(i).FillRand();
  }
  
  
  //! values of non-zero entries are set to zero
  template<class T, class Prop, class Storage, class Allocator>  
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::Zero()
  {
    diag_value_.Zero();
    for (int i = 0; i < mat22_.GetM(); i++)
      mat22_(i).Zero();
    
    for (int i = 0; i < mat33_.GetM(); i++)
      mat33_(i).Zero();
    
    for (int i = 0; i < mat44_.GetM(); i++)
      mat44_(i).Zero();
    
    for (int i = 0; i < mat55_.GetM(); i++)
      mat55_(i).Zero();
    
    for (int i = 0; i < mat66_.GetM(); i++)
      mat66_(i).Zero();

    for (int i = 0; i < mat77_.GetM(); i++)
      mat77_(i).Zero();
    
    for (int i = 0; i < mat88_.GetM(); i++)
      mat88_(i).Zero();

    for (int i = 0; i < mat99_.GetM(); i++)
      mat99_(i).Zero();

    for (int i = 0; i < mat1010_.GetM(); i++)
      mat1010_(i).Zero();

    for (int i = 0; i < mat1111_.GetM(); i++)
      mat1111_(i).Zero();

    for (int i = 0; i < mat1212_.GetM(); i++)
      mat1212_(i).Zero();

    for (int i = 0; i < mat1313_.GetM(); i++)
      mat1313_(i).Zero();

    for (int i = 0; i < mat1414_.GetM(); i++)
      mat1414_(i).Zero();

    /* for (int i = 0; i < mat1515_.GetM(); i++)
      mat1515_(i).Zero();

    for (int i = 0; i < mat1616_.GetM(); i++)
      mat1616_(i).Zero();

    for (int i = 0; i < mat1717_.GetM(); i++)
      mat1717_(i).Zero();

    for (int i = 0; i < mat1818_.GetM(); i++)
      mat1818_(i).Zero();

    for (int i = 0; i < mat1919_.GetM(); i++)
      mat1919_(i).Zero();

    for (int i = 0; i < mat2020_.GetM(); i++)
    mat2020_(i).Zero(); */
    
    for (int i = 0; i < all_mat_.GetM(); i++)
      all_mat_(i).Zero();
  }
  

  //! returns A(i, j)
  template<class T, class Prop, class Storage, class Allocator>
  const T MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {
    int n1 = blok_num_(i);
    int n2 = blok_num_(j);
    if (size_block_(i) == size_block_(j))
      if (n1 == n2)
	{
	  int j1 = rev_num_(i);
	  int j2 = rev_num_(j);
	  switch (size_block_(i))
	    {
	    case 1 :
	      return diag_value_(n1);
	    case 2 :
	      return mat22_(n1)(j1, j2);
	    case 3 :
	      return mat33_(n1)(j1, j2);
	    case 4 :
	      return mat44_(n1)(j1, j2);
	    case 5 :
	      return mat55_(n1)(j1, j2);
	    case 6 :
	      return mat66_(n1)(j1, j2);
	    case 7 :
	      return mat77_(n1)(j1, j2);
	    case 8 :
	      return mat88_(n1)(j1, j2);
	    case 9 :
	      return mat99_(n1)(j1, j2);
	    case 10 :
	      return mat1010_(n1)(j1, j2);
	    case 11 :
	      return mat1111_(n1)(j1, j2);
            case 12 :
	      return mat1212_(n1)(j1, j2);
	    case 13 :
	      return mat1313_(n1)(j1, j2);
	    case 14 :
	      return mat1414_(n1)(j1, j2);
              /*case 15 :
	      return mat1515_(n1)(j1, j2);
	    case 16 :
	      return mat1616_(n1)(j1, j2);
	    case 17 :
	      return mat1717_(n1)(j1, j2);
	    case 18 :
	      return mat1818_(n1)(j1, j2);
	    case 19 :
	      return mat1919_(n1)(j1, j2);
	    case 20 :
            return mat2020_(n1)(j1, j2); */
	    default :
	      return all_mat_(n1)(j1, j2);
	    }
	}
    
    return T(0);
  }

  
  //! multiplies matrix by a scalar
  template <class T, class Prop, class Storage, class Allocator>
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>&
  MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::operator *= (const T& alpha)
  {
    Mlt(alpha, diag_value_);
    for (int i = 0; i < mat22_.GetM(); i++)
      Mlt(alpha, mat22_(i));

    for (int i = 0; i < mat33_.GetM(); i++)
      Mlt(alpha, mat33_(i));

    for (int i = 0; i < mat44_.GetM(); i++)
      Mlt(alpha, mat44_(i));

    for (int i = 0; i < mat55_.GetM(); i++)
      Mlt(alpha, mat55_(i));

    for (int i = 0; i < mat66_.GetM(); i++)
      Mlt(alpha, mat66_(i));

    for (int i = 0; i < mat77_.GetM(); i++)
      Mlt(alpha, mat77_(i));

    for (int i = 0; i < mat88_.GetM(); i++)
      Mlt(alpha, mat88_(i));

    for (int i = 0; i < mat99_.GetM(); i++)
      Mlt(alpha, mat99_(i));

    for (int i = 0; i < mat1010_.GetM(); i++)
      Mlt(alpha, mat1010_(i));

    for (int i = 0; i < mat1111_.GetM(); i++)
      Mlt(alpha, mat1111_(i));

    for (int i = 0; i < mat1212_.GetM(); i++)
      Mlt(alpha, mat1212_(i));

    for (int i = 0; i < mat1313_.GetM(); i++)
      Mlt(alpha, mat1313_(i));

    for (int i = 0; i < mat1414_.GetM(); i++)
      Mlt(alpha, mat1414_(i));

    /* for (int i = 0; i < mat1515_.GetM(); i++)
      Mlt(alpha, mat1515_(i));
    
    for (int i = 0; i < mat1616_.GetM(); i++)
      Mlt(alpha, mat1616_(i));

    for (int i = 0; i < mat1717_.GetM(); i++)
      Mlt(alpha, mat1717_(i));

    for (int i = 0; i < mat1818_.GetM(); i++)
      Mlt(alpha, mat1818_(i));

    for (int i = 0; i < mat1919_.GetM(); i++)
      Mlt(alpha, mat1919_(i));

    for (int i = 0; i < mat2020_.GetM(); i++)
    Mlt(alpha, mat2020_(i)); */

    for (int i = 0; i < all_mat_.GetM(); i++)
      Mlt(alpha, all_mat_(i));
    
    return *this;
  }
  
  
  //! conversion to a sparse matrix (CSR)
  /*!
    \param[out] A sparse matrix that will be equal to the current matrix
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Prop0, class Storage0>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  ConvertToSparse(Matrix<T0, Prop0, Storage0>& A) const
  {
    A.Clear(); A.Reallocate(m_, m_);
    for (int i = 0; i < diag_value_.GetM(); i++)
      A.AddInteraction(row_num_(i), row_num_(i), diag_value_(i));
    
    int nb = diag_value_.GetM();
    for (int i = 0; i < mat22_.GetM(); i++)
      {
	int nb_val = 2;
	IVect col(2); Vector<T> val(2);
	for (int k = 0; k < 2; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 2; j++)
	  {
	    for (int k = 0; k < 2; k++)
	      val(k) = mat22_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 2;
      }

    for (int i = 0; i < mat33_.GetM(); i++)
      {
	int nb_val = 3;
	IVect col(3); Vector<T> val(3);
	for (int k = 0; k < 3; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 3; j++)
	  {
	    for (int k = 0; k < 3; k++)
	      val(k) = mat33_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 3;
      }

    for (int i = 0; i < mat44_.GetM(); i++)
      {
	int nb_val = 4;
	IVect col(4); Vector<T> val(4);
	for (int k = 0; k < 4; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 4; j++)
	  {
	    for (int k = 0; k < 4; k++)
	      val(k) = mat44_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 4;
      }

    for (int i = 0; i < mat55_.GetM(); i++)
      {
	int nb_val = 5;
	IVect col(5); Vector<T> val(5);
	for (int k = 0; k < 5; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 5; j++)
	  {
	    for (int k = 0; k < 5; k++)
	      val(k) = mat55_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 5;
      }

    for (int i = 0; i < mat66_.GetM(); i++)
      {
	int nb_val = 6;
	IVect col(6); Vector<T> val(6);
	for (int k = 0; k < 6; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 6; j++)
	  {
	    for (int k = 0; k < 6; k++)
	      val(k) = mat66_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 6;
      }

    for (int i = 0; i < mat77_.GetM(); i++)
      {
	int nb_val = 7;
	IVect col(7); Vector<T> val(7);
	for (int k = 0; k < 7; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 7; j++)
	  {
	    for (int k = 0; k < 7; k++)
	      val(k) = mat77_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 7;
      }

    for (int i = 0; i < mat88_.GetM(); i++)
      {
	int nb_val = 8;
	IVect col(8); Vector<T> val(8);
	for (int k = 0; k < 8; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 8; j++)
	  {
	    for (int k = 0; k < 8; k++)
	      val(k) = mat88_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 8;
      }

    for (int i = 0; i < mat99_.GetM(); i++)
      {
	int nb_val = 9;
	IVect col(9); Vector<T> val(9);
	for (int k = 0; k < 9; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 9; j++)
	  {
	    for (int k = 0; k < 9; k++)
	      val(k) = mat99_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 9;
      }

    for (int i = 0; i < mat1010_.GetM(); i++)
      {
	int nb_val = 10;
	IVect col(10); Vector<T> val(10);
	for (int k = 0; k < 10; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 10; j++)
	  {
	    for (int k = 0; k < 10; k++)
	      val(k) = mat1010_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 10;
      }

    for (int i = 0; i < mat1111_.GetM(); i++)
      {
	int nb_val = 11;
	IVect col(11); Vector<T> val(11);
	for (int k = 0; k < 11; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 11; j++)
	  {
	    for (int k = 0; k < 11; k++)
	      val(k) = mat1111_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 11;
      }

    for (int i = 0; i < mat1212_.GetM(); i++)
      {
	int nb_val = 12;
	IVect col(12); Vector<T> val(12);
	for (int k = 0; k < 12; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 12; j++)
	  {
	    for (int k = 0; k < 12; k++)
	      val(k) = mat1212_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 12;
      }

    for (int i = 0; i < mat1313_.GetM(); i++)
      {
	int nb_val = 13;
	IVect col(13); Vector<T> val(13);
	for (int k = 0; k < 13; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 13; j++)
	  {
	    for (int k = 0; k < 13; k++)
	      val(k) = mat1313_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 13;
      }

    for (int i = 0; i < mat1414_.GetM(); i++)
      {
	int nb_val = 14;
	IVect col(14); Vector<T> val(14);
	for (int k = 0; k < 14; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 14; j++)
	  {
	    for (int k = 0; k < 14; k++)
	      val(k) = mat1414_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 14;
      }

    /*for (int i = 0; i < mat1515_.GetM(); i++)
      {
	int nb_val = 15;
	IVect col(15); Vector<T> val(15);
	for (int k = 0; k < 15; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 15; j++)
	  {
	    for (int k = 0; k < 15; k++)
	      val(k) = mat1515_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 15;
      }

    for (int i = 0; i < mat1616_.GetM(); i++)
      {
	int nb_val = 16;
	IVect col(16); Vector<T> val(16);
	for (int k = 0; k < 16; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 16; j++)
	  {
	    for (int k = 0; k < 16; k++)
	      val(k) = mat1616_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 16;
      }

    for (int i = 0; i < mat1717_.GetM(); i++)
      {
	int nb_val = 17;
	IVect col(17); Vector<T> val(17);
	for (int k = 0; k < 17; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 17; j++)
	  {
	    for (int k = 0; k < 17; k++)
	      val(k) = mat1717_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 17;
      }

    for (int i = 0; i < mat1818_.GetM(); i++)
      {
	int nb_val = 18;
	IVect col(18); Vector<T> val(18);
	for (int k = 0; k < 18; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 18; j++)
	  {
	    for (int k = 0; k < 18; k++)
	      val(k) = mat1818_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 18;
      }

    for (int i = 0; i < mat1919_.GetM(); i++)
      {
	int nb_val = 19;
	IVect col(19); Vector<T> val(19);
	for (int k = 0; k < 19; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 19; j++)
	  {
	    for (int k = 0; k < 19; k++)
	      val(k) = mat1919_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 19;
      }

    for (int i = 0; i < mat2020_.GetM(); i++)
      {
	int nb_val = 20;
	IVect col(20); Vector<T> val(20);
	for (int k = 0; k < 20; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < 20; j++)
	  {
	    for (int k = 0; k < 20; k++)
	      val(k) = mat2020_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += 20;
        } */

    for (int i = 0; i < all_mat_.GetM(); i++)
      {
	int nb_val = all_mat_(i).GetM();
	IVect col(nb_val); Vector<T> val(nb_val);
	for (int k = 0; k < nb_val; k++)
	  col(k) = row_num_(nb+k);
	
	for (int j = 0; j < nb_val; j++)
	  {
	    for (int k = 0; k < nb_val; k++)
	      val(k) = all_mat_(i)(j, k);
	    
	    A.AddInteractionRow(col(j), nb_val, col, val);
	  }
	
	nb += nb_val;
      }
  }
  
  
  //! writes the matrix on a file
  template<class T, class Prop, class Storage, class Allocator>
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
      throw IOError("Matrix_BlockDiagonal::WriteText(string FileName)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close(); 
  }
  
  
  //! writes the matrix on a file stream
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::WriteText(ostream& FileStream) const
  {
    int nb = 0;
    for (int i = 0; i < diag_value_.GetM(); i++)
      FileStream << row_num_(i)+1 << ' ' << row_num_(i)+1 << ' ' << diag_value_(i) << '\n';
    
    nb = diag_value_.GetM();
    for (int i = 0; i < mat22_.GetM(); i++)
      {
	for (int j = 0; j < 2; j++)
	  for (int k = 0; k < 2; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat22_(i)(j, k) << '\n';
	
	nb += 2;
      }

    for (int i = 0; i < mat33_.GetM(); i++)
      {
	for (int j = 0; j < 3; j++)
	  for (int k = 0; k < 3; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat33_(i)(j, k) << '\n';
	
	nb += 3;
      }

    for (int i = 0; i < mat44_.GetM(); i++)
      {
	for (int j = 0; j < 4; j++)
	  for (int k = 0; k < 4; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat44_(i)(j, k) << '\n';
	
	nb += 4;
      }

    for (int i = 0; i < mat55_.GetM(); i++)
      {
	for (int j = 0; j < 5; j++)
	  for (int k = 0; k < 5; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat55_(i)(j, k) << '\n';
	
	nb += 5;
      }

    for (int i = 0; i < mat66_.GetM(); i++)
      {
	for (int j = 0; j < 6; j++)
	  for (int k = 0; k < 6; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat66_(i)(j, k) << '\n';
	
	nb += 6;
      }

    for (int i = 0; i < mat77_.GetM(); i++)
      {
	for (int j = 0; j < 7; j++)
	  for (int k = 0; k < 7; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat77_(i)(j, k) << '\n';
	
	nb += 7;
      }

    for (int i = 0; i < mat88_.GetM(); i++)
      {
	for (int j = 0; j < 8; j++)
	  for (int k = 0; k < 8; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat88_(i)(j, k) << '\n';
	
	nb += 8;
      }

    for (int i = 0; i < mat99_.GetM(); i++)
      {
	for (int j = 0; j < 9; j++)
	  for (int k = 0; k < 9; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat99_(i)(j, k) << '\n';
	
	nb += 9;
      }

    for (int i = 0; i < mat1010_.GetM(); i++)
      {
	for (int j = 0; j < 10; j++)
	  for (int k = 0; k < 10; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat1010_(i)(j, k) << '\n';
	
	nb += 10;
      }

    for (int i = 0; i < mat1111_.GetM(); i++)
      {
	for (int j = 0; j < 11; j++)
	  for (int k = 0; k < 11; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat1111_(i)(j, k) << '\n';
	
	nb += 11;
      }

    for (int i = 0; i < mat1212_.GetM(); i++)
      {
	for (int j = 0; j < 12; j++)
	  for (int k = 0; k < 12; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat1212_(i)(j, k) << '\n';
	
	nb += 12;
      }

    for (int i = 0; i < mat1313_.GetM(); i++)
      {
	for (int j = 0; j < 13; j++)
	  for (int k = 0; k < 13; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat1313_(i)(j, k) << '\n';
	
	nb += 13;
      }

    for (int i = 0; i < mat1414_.GetM(); i++)
      {
	for (int j = 0; j < 14; j++)
	  for (int k = 0; k < 14; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << mat1414_(i)(j, k) << '\n';
	
	nb += 14;
      }

      /* for (int i = 0; i < mat1515_.GetM(); i++)
      {
	for (int j = 0; j < 15; j++)
	  for (int k = 0; k < 15; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
            << row_num_(nb+k)+1 << ' ' << mat1515_(i)(j, k) << '\n';
	
	nb += 15;
      }

    for (int i = 0; i < mat1616_.GetM(); i++)
      {
	for (int j = 0; j < 16; j++)
	  for (int k = 0; k < 16; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
            << row_num_(nb+k)+1 << ' ' << mat1616_(i)(j, k) << '\n';
	
	nb += 16;
      }

    for (int i = 0; i < mat1717_.GetM(); i++)
      {
	for (int j = 0; j < 17; j++)
	  for (int k = 0; k < 17; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
            << row_num_(nb+k)+1 << ' ' << mat1717_(i)(j, k) << '\n';
	
	nb += 17;
      }

    for (int i = 0; i < mat1818_.GetM(); i++)
      {
	for (int j = 0; j < 18; j++)
	  for (int k = 0; k < 18; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
            << row_num_(nb+k)+1 << ' ' << mat1818_(i)(j, k) << '\n';
	
	nb += 18;
      }

    for (int i = 0; i < mat1919_.GetM(); i++)
      {
	for (int j = 0; j < 19; j++)
	  for (int k = 0; k < 19; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
            << row_num_(nb+k)+1 << ' ' << mat1919_(i)(j, k) << '\n';
	
	nb += 19;
      }

    for (int i = 0; i < mat2020_.GetM(); i++)
      {
	for (int j = 0; j < 20; j++)
	  for (int k = 0; k < 20; k++)
	    FileStream << row_num_(nb+j)+1 << ' '
            << row_num_(nb+k)+1 << ' ' << mat2020_(i)(j, k) << '\n';
	
	nb += 20;
        } */

    for (int i = 0; i < all_mat_.GetM(); i++)
      {
	for (int j = 0; j < all_mat_(i).GetM(); j++)
	  for (int k = 0; k < all_mat_(i).GetM(); k++)
	    FileStream << row_num_(nb+j)+1 << ' '
                       << row_num_(nb+k)+1 << ' ' << all_mat_(i)(j, k) << '\n';
	
	nb += all_mat_(i).GetM();
      }
  }
  
  
  //! displays matrix
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, BlockDiagRow, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }


  //! displays matrix
  template<class T, class Prop, class Allocator>
  ostream& operator<<(ostream& out, const Matrix<T, Prop, BlockDiagRowSym, Allocator>& A)
  {
    A.WriteText(out);
    return out;
  }
  
  
  //! computes Y = Y + alpha*A*X
  template<class T, class Prop, class Storage, class Allocator>
  template< class T0, class T1, class Allocator1, class T2, class Allocator2>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  MltAdd_(const T0& alpha, const class_SeldonNoTrans& trans,
          const Vector<T1, VectFull, Allocator1>& X,
          Vector<T2, VectFull, Allocator2>& Y) const
  {
    for (int i = 0; i < diag_value_.GetM(); i++)
      Y(row_num_(i)) += alpha*diag_value_(i)*X(row_num_(i));
    
    int nb = diag_value_.GetM();
    {
      TinyVector<T2, 2> x, y;
      for (int i = 0; i < mat22_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  
	  Mlt(mat22_(i), x, y);

	  Y(row_num_(nb)) += alpha*y(0);
	  Y(row_num_(nb+1)) += alpha*y(1);
	  nb += 2;
	}
    }

    {
      TinyVector<T2, 3> x, y;
      for (int i = 0; i < mat33_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  x(2) = X(row_num_(nb+2));
	  
	  Mlt(mat33_(i), x, y);

	  Y(row_num_(nb)) += alpha*y(0);
	  Y(row_num_(nb+1)) += alpha*y(1);
	  Y(row_num_(nb+2)) += alpha*y(2);
	  nb += 3;
	}
    }

    {
      TinyVector<T2, 4> x, y;
      for (int i = 0; i < mat44_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  x(2) = X(row_num_(nb+2));
	  x(3) = X(row_num_(nb+3));
	  
	  Mlt(mat44_(i), x, y);

	  Y(row_num_(nb)) += alpha*y(0);
	  Y(row_num_(nb+1)) += alpha*y(1);
	  Y(row_num_(nb+2)) += alpha*y(2);
	  Y(row_num_(nb+3)) += alpha*y(3);
	  nb += 4;
	}
    }

    {
      TinyVector<T2, 5> x, y;
      for (int i = 0; i < mat55_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
          Mlt(mat55_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
	  nb += 5;
	}
    }

    {
      TinyVector<T2, 6> x, y;
      for (int i = 0; i < mat66_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
          Mlt(mat66_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
	  nb += 6;
	}
    }

    {
      TinyVector<T2, 7> x, y;
      for (int i = 0; i < mat77_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
          Mlt(mat77_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
	  nb += 7;
	}
    }

    {
      TinyVector<T2, 8> x, y;
      for (int i = 0; i < mat88_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
          Mlt(mat88_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
	  nb += 8;
	}
    }

    {
      TinyVector<T2, 9> x, y;
      for (int i = 0; i < mat99_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
          Mlt(mat99_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
	  nb += 9;
	}
    }

    {
      TinyVector<T2, 10> x, y;
      for (int i = 0; i < mat1010_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1010_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 10;
        }
    }


    {
      TinyVector<T2, 11> x, y;
      for (int i = 0; i < mat1111_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1111_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 11;
        }
    }


   {
      TinyVector<T2, 12> x, y;
      for (int i = 0; i < mat1212_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1212_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 12;
        }
    }


    {
      TinyVector<T2, 13> x, y;
      for (int i = 0; i < mat1313_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1313_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 13;
        }
    }


    {
      TinyVector<T2, 14> x, y;
      for (int i = 0; i < mat1414_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1414_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 14;
        }
    }

    /* {
      TinyVector<T2, 15> x, y;
      for (int i = 0; i < mat1515_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1515_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 15;
        }
    }

    {
      TinyVector<T2, 16> x, y;
      for (int i = 0; i < mat1616_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1616_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 16;
        }
    }

    {
      TinyVector<T2, 17> x, y;
      for (int i = 0; i < mat1717_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1717_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 17;
        }
    }

    {
      TinyVector<T2, 18> x, y;
      for (int i = 0; i < mat1818_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1818_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 18;
        }
    }

    {
      TinyVector<T2, 19> x, y;
      for (int i = 0; i < mat1919_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1919_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 19;
        }
    }

    {
      TinyVector<T2, 20> x, y;
      for (int i = 0; i < mat2020_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat2020_(i), x, y);
          AddVector(alpha, x, row_num_, nb, Y);
          nb += 20;
        }
        }*/          

    Vector<T> xs(size_max_blok_), ys(size_max_blok_), x, y;
    for (int i = 0; i < all_mat_.GetM(); i++)
      {
        T2 zero; SetComplexZero(zero);

	int size_mat = all_mat_(i).GetM();
        x.SetData(size_mat, xs.GetData());
        y.SetData(size_mat, ys.GetData());
	for (int j = 0; j < size_mat; j++)
	  {
            x(j) = X(row_num_(nb+j));
            y(j) = 0;
          }
	
	MltAdd(alpha, all_mat_(i), x, zero, y);
	
	for (int j = 0; j < size_mat; j++)
	  Y(row_num_(nb+j)) += y(j);
        
        x.Nullify(); y.Nullify();
	nb += size_mat;
      }    
  }
  

  //! computes Y = Y + alpha*A*X
  template<class T, class Prop, class Storage, class Allocator>
  template< class T0, class T1, class Allocator1, class T2, class Allocator2>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  MltAdd_(const T0& alpha, const class_SeldonTrans& trans,
          const Vector<T1, VectFull, Allocator1>& X,
          Vector<T2, VectFull, Allocator2>& Y) const
  {
    for (int i = 0; i < diag_value_.GetM(); i++)
      Y(row_num_(i)) += alpha*diag_value_(i)*X(row_num_(i));
    
    int nb = diag_value_.GetM();
    {
      TinyVector<T2, 2> x, y;
      for (int i = 0; i < mat22_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  
	  MltTrans(mat22_(i), x, y);

	  Y(row_num_(nb)) += alpha*y(0);
	  Y(row_num_(nb+1)) += alpha*y(1);
	  nb += 2;
	}
    }

    {
      TinyVector<T2, 3> x, y;
      for (int i = 0; i < mat33_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  x(2) = X(row_num_(nb+2));
	  
	  MltTrans(mat33_(i), x, y);

	  Y(row_num_(nb)) += alpha*y(0);
	  Y(row_num_(nb+1)) += alpha*y(1);
	  Y(row_num_(nb+2)) += alpha*y(2);
	  nb += 3;
	}
    }

    {
      TinyVector<T2, 4> x, y;
      for (int i = 0; i < mat44_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  x(2) = X(row_num_(nb+2));
	  x(3) = X(row_num_(nb+3));
	  
	  MltTrans(mat44_(i), x, y);

	  Y(row_num_(nb)) += alpha*y(0);
	  Y(row_num_(nb+1)) += alpha*y(1);
	  Y(row_num_(nb+2)) += alpha*y(2);
	  Y(row_num_(nb+3)) += alpha*y(3);
	  nb += 4;
	}
    }

    {
      TinyVector<T2, 5> x, y;
      for (int i = 0; i < mat55_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat55_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 5;
	}
    }

    {
      TinyVector<T2, 6> x, y;
      for (int i = 0; i < mat66_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat66_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 6;
	}
    }

    {
      TinyVector<T2, 7> x, y;
      for (int i = 0; i < mat77_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat77_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 7;
	}
    }

    {
      TinyVector<T2, 8> x, y;
      for (int i = 0; i < mat88_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat88_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 8;
	}
    }

    {
      TinyVector<T2, 9> x, y;
      for (int i = 0; i < mat99_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat99_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 9;
	}
    }

    {
      TinyVector<T2, 10> x, y;
      for (int i = 0; i < mat1010_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1010_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 10;
	}
    }

    {
      TinyVector<T2, 11> x, y;
      for (int i = 0; i < mat1111_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1111_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 11;
	}
    }

    {
      TinyVector<T2, 12> x, y;
      for (int i = 0; i < mat1212_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1212_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 12;
	}
    }

    {
      TinyVector<T2, 13> x, y;
      for (int i = 0; i < mat1313_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1313_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 13;
	}
    }

    {
      TinyVector<T2, 14> x, y;
      for (int i = 0; i < mat1414_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1414_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 14;
	}
    }

    /* {
      TinyVector<T2, 15> x, y;
      for (int i = 0; i < mat1515_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1515_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 15;
	}
    }

    {
      TinyVector<T2, 16> x, y;
      for (int i = 0; i < mat1616_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1616_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 16;
	}
    }

    {
      TinyVector<T2, 17> x, y;
      for (int i = 0; i < mat1717_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1717_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 17;
	}
    }

    {
      TinyVector<T2, 18> x, y;
      for (int i = 0; i < mat1818_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1818_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 18;
	}
    }

    {
      TinyVector<T2, 19> x, y;
      for (int i = 0; i < mat1919_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat1919_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 19;
	}
    }

    {
      TinyVector<T2, 20> x, y;
      for (int i = 0; i < mat2020_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);
	  MltTrans(mat2020_(i), x, y);
          AddVector(alpha, y, row_num_, nb, Y);
	  nb += 20;
	}
        } */

    Vector<T> xs(size_max_blok_), ys(size_max_blok_), x, y;
    for (int i = 0; i < all_mat_.GetM(); i++)
      {
        T2 zero; SetComplexZero(zero);
	int size_mat = all_mat_(i).GetM();
        x.SetData(size_mat, xs.GetData());
        y.SetData(size_mat, ys.GetData());
	for (int j = 0; j < size_mat; j++)
	  {
            x(j) = X(row_num_(nb+j));
            y(j) = 0;
          }
	
	MltAdd(alpha, SeldonTrans, all_mat_(i), x, zero, y);
	
	for (int j = 0; j < size_mat; j++)
	  Y(row_num_(nb+j)) += y(j);
	
        x.Nullify(); y.Nullify();
	nb += size_mat;
      }    
  }


  //! computes Y = Y + alpha*A*X
  template<class T, class Prop, class Storage, class Allocator>
  template< class T0, class T1, class Allocator1, class T2, class Allocator2>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  MltAdd_(const T0& alpha, const class_SeldonConjTrans& trans,
          const Vector<T1, VectFull, Allocator1>& X,
          Vector<T2, VectFull, Allocator2>& Y) const
  {
    for (int i = 0; i < diag_value_.GetM(); i++)
      Y(row_num_(i)) += alpha*conj(diag_value_(i))*X(row_num_(i));
    
    int nb = diag_value_.GetM();
    {
      TinyVector<T2, 2> x, y;
      for (int i = 0; i < mat22_.GetM(); i++)
	{
	  x(0) = conj(X(row_num_(nb)));
	  x(1) = conj(X(row_num_(nb+1)));
	  
	  MltTrans(mat22_(i), x, y);

	  Y(row_num_(nb)) += alpha*conj(y(0));
          Y(row_num_(nb+1)) += alpha*conj(y(1));
	  nb += 2;
	}
    }

    {
      TinyVector<T2, 3> x, y;
      for (int i = 0; i < mat33_.GetM(); i++)
	{
	  x(0) = conj(X(row_num_(nb)));
	  x(1) = conj(X(row_num_(nb+1)));
	  x(2) = conj(X(row_num_(nb+2)));
	  
	  MltTrans(mat33_(i), x, y);

	  Y(row_num_(nb)) += alpha*conj(y(0));
	  Y(row_num_(nb+1)) += alpha*conj(y(1));
	  Y(row_num_(nb+2)) += alpha*conj(y(2));
	  nb += 3;
	}
    }

    {
      TinyVector<T2, 4> x, y;
      for (int i = 0; i < mat44_.GetM(); i++)
	{
	  x(0) = conj(X(row_num_(nb)));
	  x(1) = conj(X(row_num_(nb+1)));
	  x(2) = conj(X(row_num_(nb+2)));
	  x(3) = conj(X(row_num_(nb+3)));
	  
	  MltTrans(mat44_(i), x, y);

	  Y(row_num_(nb)) += alpha*conj(y(0));
	  Y(row_num_(nb+1)) += alpha*conj(y(1));
	  Y(row_num_(nb+2)) += alpha*conj(y(2));
	  Y(row_num_(nb+3)) += alpha*conj(y(3));
	  nb += 4;
	}
    }

    {
      TinyVector<T2, 5> x, y;
      for (int i = 0; i < mat55_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat55_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 5;
	}
    }

    {
      TinyVector<T2, 6> x, y;
      for (int i = 0; i < mat66_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat66_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 6;
	}
    }

    {
      TinyVector<T2, 7> x, y;
      for (int i = 0; i < mat77_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat77_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 7;
	}
    }

    {
      TinyVector<T2, 8> x, y;
      for (int i = 0; i < mat88_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat88_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 8;
	}
    }

    {
      TinyVector<T2, 9> x, y;
      for (int i = 0; i < mat99_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat99_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 9;
	}
    }

    {
      TinyVector<T2, 10> x, y;
      for (int i = 0; i < mat1010_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1010_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 10;
	}
    }

    {
      TinyVector<T2, 11> x, y;
      for (int i = 0; i < mat1111_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1111_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 11;
	}
        }

    {
      TinyVector<T2, 12> x, y;
      for (int i = 0; i < mat1212_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1212_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 12;
	}
    }

    {
      TinyVector<T2, 13> x, y;
      for (int i = 0; i < mat1313_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1313_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 13;
	}
    }

    {
      TinyVector<T2, 14> x, y;
      for (int i = 0; i < mat1414_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1414_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 14;
	}
    }

    /* {
      TinyVector<T2, 15> x, y;
      for (int i = 0; i < mat1515_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1515_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 15;
	}
    }

    {
      TinyVector<T2, 16> x, y;
      for (int i = 0; i < mat1616_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1616_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 16;
	}
    }

    {
      TinyVector<T2, 17> x, y;
      for (int i = 0; i < mat1717_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1717_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 17;
	}
    }

    {
      TinyVector<T2, 18> x, y;
      for (int i = 0; i < mat1818_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1818_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 18;
	}
    }

    {
      TinyVector<T2, 19> x, y;
      for (int i = 0; i < mat1919_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat1919_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 19;
	}
    }

    {
      TinyVector<T2, 20> x, y;
      for (int i = 0; i < mat2020_.GetM(); i++)
	{
          ExtractVectorConj(X, row_num_, nb, x);
	  MltTrans(mat2020_(i), x, y);
          AddVectorConj(alpha, y, row_num_, nb, Y);
	  nb += 20;
	}
        }*/

    Vector<T> xs(size_max_blok_), ys(size_max_blok_), x, y;
    for (int i = 0; i < all_mat_.GetM(); i++)
      {
        T2 zero; SetComplexZero(zero);
	int size_mat = all_mat_(i).GetM();
	x.SetData(size_mat, xs.GetData());
        y.SetData(size_mat, ys.GetData());
	for (int j = 0; j < size_mat; j++)
	  {
            x(j) = X(row_num_(nb+j));
            y(j) = 0;
          }
	
	MltAdd(alpha, SeldonConjTrans, all_mat_(i), x, zero, y);
	
	for (int j = 0; j < size_mat; j++)
	  Y(row_num_(nb+j)) += y(j);
	
        x.Nullify(); y.Nullify();
	nb += size_mat;
      }
    
  }


  //! computes X = A*X (X is overwritten with the product of A by X)
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Allocator2>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  Mlt_(Vector<T2, VectFull, Allocator2>& X) const
  {
    for (int i = 0; i < diag_value_.GetM(); i++)
      X(row_num_(i)) *= diag_value_(i);
    
    int nb = diag_value_.GetM();
    
    {
      TinyVector<T2, 2> x, y;
      for (int i = 0; i < mat22_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  
	  Mlt(mat22_(i), x, y);

	  X(row_num_(nb)) = y(0);
	  X(row_num_(nb+1)) = y(1);
	  nb += 2;
	}
    }

    {
      TinyVector<T2, 3> x, y;
      for (int i = 0; i < mat33_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  x(2) = X(row_num_(nb+2));
	  
	  Mlt(mat33_(i), x, y);

	  X(row_num_(nb)) = y(0);
	  X(row_num_(nb+1)) = y(1);
	  X(row_num_(nb+2)) = y(2);
	  nb += 3;
	}
    }

    {
      TinyVector<T2, 4> x, y;
      for (int i = 0; i < mat44_.GetM(); i++)
	{
	  x(0) = X(row_num_(nb));
	  x(1) = X(row_num_(nb+1));
	  x(2) = X(row_num_(nb+2));
	  x(3) = X(row_num_(nb+3));
	  
	  Mlt(mat44_(i), x, y);

	  X(row_num_(nb)) = y(0);
	  X(row_num_(nb+1)) = y(1);
	  X(row_num_(nb+2)) = y(2);
	  X(row_num_(nb+3)) = y(3);
	  nb += 4;
	}
    }

    {
      TinyVector<T2, 5> x, y;
      for (int i = 0; i < mat55_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat55_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
	  nb += 5;
	}
    }

    {
      TinyVector<T2, 6> x, y;
      for (int i = 0; i < mat66_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat66_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
	  nb += 6;
	}
    }

    {
      TinyVector<T2, 7> x, y;
      for (int i = 0; i < mat77_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat77_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
	  nb += 7;
	}
    }

    {
      TinyVector<T2, 8> x, y;
      for (int i = 0; i < mat88_.GetM(); i++)
	{
	  ExtractVector(X, row_num_, nb, x);          
          Mlt(mat88_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 8;
	}
    }

    {
      TinyVector<T2, 9> x, y;
      for (int i = 0; i < mat99_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat99_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
	  nb += 9;
	}
    }

    {
      TinyVector<T2, 10> x, y;
      for (int i = 0; i < mat1010_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1010_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 10;
        }
    }

    {
      TinyVector<T2, 11> x, y;
      for (int i = 0; i < mat1111_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1111_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 11;
        }
        }

    {
      TinyVector<T2, 12> x, y;
      for (int i = 0; i < mat1212_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1212_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 12;
        }
    }

    {
      TinyVector<T2, 13> x, y;
      for (int i = 0; i < mat1313_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1313_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 13;
        }
    }

    {
      TinyVector<T2, 14> x, y;
      for (int i = 0; i < mat1414_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1414_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 14;
        }
    }

    /*{
      TinyVector<T2, 15> x, y;
      for (int i = 0; i < mat1515_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1515_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 15;
        }
    }

    {
      TinyVector<T2, 16> x, y;
      for (int i = 0; i < mat1616_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1616_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 16;
        }
    }

    {
      TinyVector<T2, 17> x, y;
      for (int i = 0; i < mat1717_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1717_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 17;
        }
    }

    {
      TinyVector<T2, 18> x, y;
      for (int i = 0; i < mat1818_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1818_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 18;
        }
    }

    {
      TinyVector<T2, 19> x, y;
      for (int i = 0; i < mat1919_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat1919_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 19;
        }
    }

    {
      TinyVector<T2, 20> x, y;
      for (int i = 0; i < mat2020_.GetM(); i++)
	{
          ExtractVector(X, row_num_, nb, x);          
          Mlt(mat2020_(i), x, y);
          ExtractVector(y, row_num_, nb, X);
          nb += 20;
        }
        } */
    
    Vector<T> xs(size_max_blok_), ys(size_max_blok_), x, y;
    for (int i = 0; i < all_mat_.GetM(); i++)
      {
	int size_mat = all_mat_(i).GetM();
        x.SetData(size_mat, xs.GetData());
        y.SetData(size_mat, ys.GetData());
	for (int j = 0; j < size_mat; j++)
	  {
            x(j) = X(row_num_(nb+j));
            y(j) = 0;
          }
	
	Mlt(all_mat_(i), x, y);
	
	for (int j = 0; j < size_mat; j++)
	  X(row_num_(nb+j)) = y(j);
	
        x.Nullify(); y.Nullify();
	nb += size_mat;
      }
    
  }

  
  //! computes C = this * B
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>
  ::Mlt_(const MatrixSym_BlockDiagonal<T0, Prop0, Storage0, Allocator0>& B,
         MatrixSym_BlockDiagonal<T1, Prop1, Storage1, Allocator1>& C) const
  {
    // we assume that patterns of this, B and C are the same
    for (int i = 0; i < diag_value_.GetM(); i++)
      C.diag_value_(i) = diag_value_(i) * B.diag_value_(i);
    
    for (int i = 0; i < mat22_.GetM(); i++)
      Mlt(mat22_(i), B.mat22_(i), C.mat22_(i));

    for (int i = 0; i < mat33_.GetM(); i++)
      Mlt(mat33_(i), B.mat33_(i), C.mat33_(i));

    for (int i = 0; i < mat44_.GetM(); i++)
      Mlt(mat44_(i), B.mat44_(i), C.mat44_(i));

    for (int i = 0; i < mat55_.GetM(); i++)
      Mlt(mat55_(i), B.mat55_(i), C.mat55_(i));

    for (int i = 0; i < mat66_.GetM(); i++)
      Mlt(mat66_(i), B.mat66_(i), C.mat66_(i));

    for (int i = 0; i < mat77_.GetM(); i++)
      Mlt(mat77_(i), B.mat77_(i), C.mat77_(i));
    
    for (int i = 0; i < mat88_.GetM(); i++)
      MltMatrix(mat88_(i), B.mat88_(i), C.mat88_(i));
    
    for (int i = 0; i < mat99_.GetM(); i++)
      MltMatrix(mat99_(i), B.mat99_(i), C.mat99_(i));

    for (int i = 0; i < mat1010_.GetM(); i++)
      MltMatrix(mat1010_(i), B.mat1010_(i), C.mat1010_(i));

    for (int i = 0; i < mat1111_.GetM(); i++)
      MltMatrix(mat1111_(i), B.mat1111_(i), C.mat1111_(i));
    
    for (int i = 0; i < mat1212_.GetM(); i++)
      MltMatrix(mat1212_(i), B.mat1212_(i), C.mat1212_(i));

    for (int i = 0; i < mat1313_.GetM(); i++)
      MltMatrix(mat1313_(i), B.mat1313_(i), C.mat1313_(i));

    for (int i = 0; i < mat1414_.GetM(); i++)
      MltMatrix(mat1414_(i), B.mat1414_(i), C.mat1414_(i));

    /* for (int i = 0; i < mat1515_.GetM(); i++)
      Mlt(mat1515_(i), B.mat1515_(i), C.mat1515_(i));

    for (int i = 0; i < mat1616_.GetM(); i++)
      Mlt(mat1616_(i), B.mat1616_(i), C.mat1616_(i));

    for (int i = 0; i < mat1717_.GetM(); i++)
      Mlt(mat1717_(i), B.mat1717_(i), C.mat1717_(i));

    for (int i = 0; i < mat1818_.GetM(); i++)
      Mlt(mat1818_(i), B.mat1818_(i), C.mat1818_(i));

    for (int i = 0; i < mat1919_.GetM(); i++)
      Mlt(mat1919_(i), B.mat1919_(i), C.mat1919_(i));

    for (int i = 0; i < mat2020_.GetM(); i++)
    Mlt(mat2020_(i), B.mat2020_(i), C.mat2020_(i)); */
    
    for (int i = 0; i < all_mat_.GetM(); i++)
      Mlt(all_mat_(i), B.all_mat_(i), C.all_mat_(i));
  }

  
  //! replaces block-diagonal matrix by its inverse
  template<class T, class Prop, class Storage, class Allocator>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::GetInverse_()
  {
    for (int i = 0; i < diag_value_.GetM(); i++)
      diag_value_(i) = 1.0/diag_value_(i);

    for (int i = 0; i < mat22_.GetM(); i++)
      GetInverse(mat22_(i));
    
    for (int i = 0; i < mat33_.GetM(); i++)
      GetInverse(mat33_(i));

    for (int i = 0; i < mat44_.GetM(); i++)
      GetInverse(mat44_(i));

    for (int i = 0; i < mat55_.GetM(); i++)
      GetInverse(mat55_(i));

    for (int i = 0; i < mat66_.GetM(); i++)
      GetInverse(mat66_(i));

    for (int i = 0; i < mat77_.GetM(); i++)
      GetInverse(mat77_(i));

    for (int i = 0; i < mat88_.GetM(); i++)
      GetInverseMatrix(mat88_(i));

    for (int i = 0; i < mat99_.GetM(); i++)
      GetInverseMatrix(mat99_(i));

    for (int i = 0; i < mat1010_.GetM(); i++)
      GetInverseMatrix(mat1010_(i));

    for (int i = 0; i < mat1111_.GetM(); i++)
      GetInverseMatrix(mat1111_(i));
    
    for (int i = 0; i < mat1212_.GetM(); i++)
      GetInverseMatrix(mat1212_(i));

    for (int i = 0; i < mat1313_.GetM(); i++)
      GetInverseMatrix(mat1313_(i));

    for (int i = 0; i < mat1414_.GetM(); i++)
      GetInverseMatrix(mat1414_(i));

    /*for (int i = 0; i < mat1515_.GetM(); i++)
      GetInverse(mat1515_(i));

    for (int i = 0; i < mat1616_.GetM(); i++)
      GetInverse(mat1616_(i));

    for (int i = 0; i < mat1717_.GetM(); i++)
      GetInverse(mat1717_(i));

    for (int i = 0; i < mat1818_.GetM(); i++)
      GetInverse(mat1818_(i));

    for (int i = 0; i < mat1919_.GetM(); i++)
      GetInverse(mat1919_(i));

    for (int i = 0; i < mat2020_.GetM(); i++)
    GetInverse(mat2020_(i)); */

    for (int i = 0; i < all_mat_.GetM(); i++)
      GetInverse(all_mat_(i));
  }
  
  
  //! adds alpha*A to the current matrix
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
  void MatrixSym_BlockDiagonal<T, Prop, Storage, Allocator>::
  Add_(const T0& alpha, const Matrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    bool same_pattern = true;
    for (int i = 0; i < size_block_.GetM(); i++)
      if (A.size_block_(i) != size_block_(i))
	same_pattern = false;
    
    if (same_pattern)
      for (int i = 0; i < row_num_.GetM(); i++)
	if (A.row_num_(i) != row_num_(i))
	  same_pattern = false;
    
    if (same_pattern)
      {
	// we add directly values if the sparsity pattern is the same
	Add(alpha, A.diag_value_, diag_value_);
	for (int i = 0; i < mat22_.GetM(); i++)
	  Add(alpha, A.mat22_(i), mat22_(i));
	
	for (int i = 0; i < mat33_.GetM(); i++)
	  Add(alpha, A.mat33_(i), mat33_(i));
	
	for (int i = 0; i < mat44_.GetM(); i++)
	  Add(alpha, A.mat44_(i), mat44_(i));
	
	for (int i = 0; i < mat55_.GetM(); i++)
	  Add(alpha, A.mat55_(i), mat55_(i));
	
	for (int i = 0; i < mat66_.GetM(); i++)
	  Add(alpha, A.mat66_(i), mat66_(i));
	
	for (int i = 0; i < mat77_.GetM(); i++)
	  Add(alpha, A.mat77_(i), mat77_(i));

	for (int i = 0; i < mat88_.GetM(); i++)
	  Add(alpha, A.mat88_(i), mat88_(i));
	
	for (int i = 0; i < mat99_.GetM(); i++)
	  Add(alpha, A.mat99_(i), mat99_(i));

	for (int i = 0; i < mat1010_.GetM(); i++)
	  Add(alpha, A.mat1010_(i), mat1010_(i));

	for (int i = 0; i < mat1111_.GetM(); i++)
	  Add(alpha, A.mat1111_(i), mat1111_(i));

	for (int i = 0; i < mat1212_.GetM(); i++)
	  Add(alpha, A.mat1212_(i), mat1212_(i));

	for (int i = 0; i < mat1313_.GetM(); i++)
	  Add(alpha, A.mat1313_(i), mat1313_(i));

	for (int i = 0; i < mat1414_.GetM(); i++)
	  Add(alpha, A.mat1414_(i), mat1414_(i));

	/* for (int i = 0; i < mat1515_.GetM(); i++)
	  Add(alpha, A.mat1515_(i), mat1515_(i));

	for (int i = 0; i < mat1616_.GetM(); i++)
	  Add(alpha, A.mat1616_(i), mat1616_(i));

	for (int i = 0; i < mat1717_.GetM(); i++)
	  Add(alpha, A.mat1717_(i), mat1717_(i));

	for (int i = 0; i < mat1818_.GetM(); i++)
	  Add(alpha, A.mat1818_(i), mat1818_(i));

	for (int i = 0; i < mat1919_.GetM(); i++)
	  Add(alpha, A.mat1919_(i), mat1919_(i));

	for (int i = 0; i < mat2020_.GetM(); i++)
        Add(alpha, A.mat2020_(i), mat2020_(i)); */
	
	for (int i = 0; i < all_mat_.GetM(); i++)
	  Add(alpha, A.all_mat_(i), all_mat_(i));
	
	return;
      }
    
    // otherwise we use conversion with sparse matrices
    Matrix<T, General, ArrayRowSparse> Asp, Bsp;
    this->ConvertToSparse(Bsp); this->Clear();
    A.ConvertToSparse(Asp); 
    Add(alpha, Asp, Bsp);
    Asp.Clear();
    
    ConvertToBlockDiagonal(Bsp, reinterpret_cast<Matrix<T, Prop, Storage, Allocator>& >(*this),
                           0, Bsp.GetM());
  }
  
  
  //! we add the value a to the entry (i, j) of the matrix
  /*!
    If the entry (i, j) is not contained in the sparsity pattern of the matrix,
    it will be ignored
   */
  template<class T, class Allocator>
  void Matrix<T, Symmetric, BlockDiagRowSym, Allocator>::AddInteraction(int i, int j, const T& a)
  {
    if (i > j)
      return;
    
    int n1 = this->blok_num_(i);
    int n2 = this->blok_num_(j);
    if (n1 == n2)
      {
	int j1 = this->rev_num_(i);
	int j2 = this->rev_num_(j);
	switch (this->size_block_(i))
	  {
	  case 1 :
	    this->diag_value_(n1) += a;
	    break;
	  case 2 :
	    this->mat22_(n1)(j1, j2) += a;
	    break;
	  case 3 :
	    this->mat33_(n1)(j1, j2) += a;
	    break;
	  case 4 :
	    this->mat44_(n1)(j1, j2) += a;
	    break;
	  case 5 :
	    this->mat55_(n1)(j1, j2) += a;
	    break;
	  case 6 :
	    this->mat66_(n1)(j1, j2) += a;
	    break;
	  case 7 :
	    this->mat77_(n1)(j1, j2) += a;
	    break;
	  case 8 :
	    this->mat88_(n1)(j1, j2) += a;
	    break;
	  case 9 :
	    this->mat99_(n1)(j1, j2) += a;
	    break;
	  case 10 :
	    this->mat1010_(n1)(j1, j2) += a;
	    break;
	  case 11 :
	    this->mat1111_(n1)(j1, j2) += a;
	    break;
          case 12 :
	    this->mat1212_(n1)(j1, j2) += a;
	    break;
	  case 13 :
	    this->mat1313_(n1)(j1, j2) += a;
	    break;
	  case 14 :
	    this->mat1414_(n1)(j1, j2) += a;
	    break;
            /* case 15 :
	    this->mat1515_(n1)(j1, j2) += a;
	    break;
	  case 16 :
	    this->mat1616_(n1)(j1, j2) += a;
	    break;
	  case 17 :
	    this->mat1717_(n1)(j1, j2) += a;
	    break;
	  case 18 :
	    this->mat1818_(n1)(j1, j2) += a;
	    break;
	  case 19 :
	    this->mat1919_(n1)(j1, j2) += a;
	    break;
	  case 20 :
	    this->mat2020_(n1)(j1, j2) += a;
	    break; */
	  default :
	    this->all_mat_(n1)(j1, j2) += a;
	    break;
	  }
      }
  }
  
  
  //! several values are added to the row i of the matrix
  /*!
    \param[in] i row number
    \param[in] nb_val number of values to add
    \param[in] col col numbers of values to add
    \param[in] val values to add
    If some entries don't belong to the sparsity pattern of the matrix,
    these entries will be ignored
   */
  template<class T, class Allocator>
  void Matrix<T, Symmetric, BlockDiagRowSym, Allocator>::
  AddInteractionRow(int i, int nb_val, const IVect& col,
		    const Vector<T>& val)
  {
    int n1 = this->blok_num_(i);
    int j1 = this->rev_num_(i);

    switch (this->size_block_(i))
      {
      case 1 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->diag_value_(n1) += val(j);
	}
	break;
      case 2 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)	      
		this->mat22_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 3 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat33_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 4 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat44_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 5 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat55_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 6 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat66_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 7 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat77_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 8 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat88_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 9 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat99_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 10 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1010_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 11 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1111_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 12 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1212_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 13 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1313_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 14 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1414_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
        /* case 15 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1515_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 16 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1616_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 17 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1717_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 18 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1818_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 19 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat1919_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 20 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->mat2020_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break; */
      default :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (i <= col(j))
	      if (this->blok_num_(col(j)) == n1)
		this->all_mat_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      }
  }
  

  //! we add the value a to the entry (i, j) of the matrix
  /*!
    If the entry (i, j) is not contained in the sparsity pattern of the matrix,
    it will be ignored
   */
  template<class T, class Allocator>
  void Matrix<T, General, BlockDiagRow, Allocator>::AddInteraction(int i, int j, const T& a)
  {
    int n1 = this->blok_num_(i);
    int n2 = this->blok_num_(j);
    if (n1 == n2)
      {
	int j1 = this->rev_num_(i);
	int j2 = this->rev_num_(j);
	switch (this->size_block_(i))
	  {
	  case 1 :
	    this->diag_value_(n1) += a;
	    break;
	  case 2 :
	    this->mat22_(n1)(j1, j2) += a;
	    break;
	  case 3 :
	    this->mat33_(n1)(j1, j2) += a;
	    break;
	  case 4 :
	    this->mat44_(n1)(j1, j2) += a;
	    break;
	  case 5 :
	    this->mat55_(n1)(j1, j2) += a;
	    break;
	  case 6 :
	    this->mat66_(n1)(j1, j2) += a;
	    break;
	  case 7 :
	    this->mat77_(n1)(j1, j2) += a;
	    break;
	  case 8 :
	    this->mat88_(n1)(j1, j2) += a;
	    break;
	  case 9 :
	    this->mat99_(n1)(j1, j2) += a;
	    break;
	  case 10 :
	    this->mat1010_(n1)(j1, j2) += a;
	    break;
	  case 11 :
	    this->mat1111_(n1)(j1, j2) += a;
	    break;
          case 12 :
	    this->mat1212_(n1)(j1, j2) += a;
	    break;
	  case 13 :
	    this->mat1313_(n1)(j1, j2) += a;
	    break;
	  case 14 :
	    this->mat1414_(n1)(j1, j2) += a;
	    break;
            /* case 15 :
            this->mat1515_(n1)(j1, j2) += a;
	    break;
	  case 16 :
	    this->mat1616_(n1)(j1, j2) += a;
	    break;
	  case 17 :
	    this->mat1717_(n1)(j1, j2) += a;
	    break;
	  case 18 :
	    this->mat1818_(n1)(j1, j2) += a;
	    break;
	  case 19 :
	    this->mat1919_(n1)(j1, j2) += a;
	    break;
	  case 20 :
	    this->mat2020_(n1)(j1, j2) += a;
	    break; */
	  default :
	    this->all_mat_(n1)(j1, j2) += a;
	    break;
	  }
      }
  }
  
  
  //! several values are added to the row i of the matrix
  /*!
    \param[in] i row number
    \param[in] nb_val number of values to add
    \param[in] col col numbers of values to add
    \param[in] val values to add
    If some entries don't belong to the sparsity pattern of the matrix,
    these entries will be ignored
   */
  template<class T, class Allocator>
  void Matrix<T, General, BlockDiagRow, Allocator>::
  AddInteractionRow(int i, int nb_val, const IVect& col,
		    const Vector<T>& val)
  {
    int n1 = this->blok_num_(i);
    int j1 = this->rev_num_(i);
    switch (this->size_block_(i))
      {
      case 1 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->diag_value_(n1) += val(j);
	}
	break;
      case 2 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)	      
	      this->mat22_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 3 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat33_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 4 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat44_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 5 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat55_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 6 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat66_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 7 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat77_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 8 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat88_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 9 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat99_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 10 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1010_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 11 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1111_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 12 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1212_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 13 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1313_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 14 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1414_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
        /* case 15 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1515_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 16 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1616_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 17 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1717_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 18 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1818_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 19 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat1919_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      case 20 :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->mat2020_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break; */
      default :
	{
	  for (int j = 0; j < nb_val; j++)
	    if (this->blok_num_(col(j)) == n1)
	      this->all_mat_(n1)(j1, this->rev_num_(col(j))) += val(j);
	}
	break;
      }
  }
  

  //! we add two block-diagonal matrices
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha, const Matrix<T1, General, BlockDiagRow, Allocator1>& B,
		 Matrix<T2, General, BlockDiagRow, Allocator2>& A)
  {
    A.Add_(alpha, B);
  }
  
  
  //! conversion of a block-diagonal matrix A to a sparse matrix B
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, General, BlockDiagRow, Alloc>& A,
                       Matrix<T, General, ArrayRowSparse, Alloc2>& B)
  {
    A.ConvertToSparse(B);
  }
  

  //! conversion of a block-diagonal matrix A to a dense matrix B
  template<class T, class Allocator1, class Allocator2>
  void Copy(const Matrix<T, General, BlockDiagRow, Allocator1>& A,
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
  
  
  //! we add two block-diagonal matrices
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void AddMatrix(const T0& alpha, const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& B,
		 Matrix<T2, Symmetric, BlockDiagRowSym, Allocator2>& A)
  {
    A.Add_(alpha, B);
  }
    
  
  //! conversion of a block-diagonal matrix A to a sparse matrix B
  template<class T, class Alloc, class Alloc2>
  void ConvertToSparse(const Matrix<T, Symmetric, BlockDiagRowSym, Alloc>& A,
                       Matrix<T, Symmetric, ArrayRowSymSparse, Alloc2>& B)
  {
    A.ConvertToSparse(B);
  }
  

  //! conversion of a block-diagonal matrix A to a dense matrix B
  template<class T, class Allocator1, class Allocator2>
  void Copy(const Matrix<T, Symmetric, BlockDiagRowSym, Allocator1>& A,
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

  
  template<class MatrixSparse>
  void MarkUsedColumnsSparseMatrix(const MatrixSparse& A, int irow, int m, int m1,
                                   Vector<bool>& UsedRow, IVect& row_found)
  {
    for (int j = 0; j < A.GetRowSize(irow); j++)
      {
        int icol = A.Index(irow, j) - m1;
        if ((icol >= 0)&&(icol < m))
          if (!UsedRow(icol))
            {
              // we found a new row to add to the block
              row_found.PushBack(icol);
              UsedRow(icol) = true;
            }
      }
  }

  
  template<class T, class Prop, class Allocator>
  void MarkUsedColumnsSparseMatrix(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
                                   int irow, int m, int m1,
                                   Vector<bool>& UsedRow, IVect& row_found)
  {
    for (int j = 0; j < A.GetRealRowSize(irow); j++)
      {
        int icol = A.IndexReal(irow, j) - m1;
        if ((icol >= 0)&&(icol < m))
          if (!UsedRow(icol))
            {
              // we found a new row to add to the block
              row_found.PushBack(icol);
              UsedRow(icol) = true;
            }
      }

    for (int j = 0; j < A.GetImagRowSize(irow); j++)
      {
        int icol = A.IndexImag(irow, j) - m1;
        if ((icol >= 0)&&(icol < m))
          if (!UsedRow(icol))
            {
              // we found a new row to add to the block
              row_found.PushBack(icol);
              UsedRow(icol) = true;
            }
      }
  }


  template<class T, class Prop, class Allocator>
  void MarkUsedColumnsSparseMatrix(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
                                   int irow, int m, int m1,
                                   Vector<bool>& UsedRow, IVect& row_found)
  {
    for (int j = 0; j < A.GetRealRowSize(irow); j++)
      {
        int icol = A.IndexReal(irow, j) - m1;
        if ((icol >= 0)&&(icol < m))
          if (!UsedRow(icol))
            {
              // we found a new row to add to the block
              row_found.PushBack(icol);
              UsedRow(icol) = true;
            }
      }

    for (int j = 0; j < A.GetImagRowSize(irow); j++)
      {
        int icol = A.IndexImag(irow, j) - m1;
        if ((icol >= 0)&&(icol < m))
          if (!UsedRow(icol))
            {
              // we found a new row to add to the block
              row_found.PushBack(icol);
              UsedRow(icol) = true;
            }
      }
  }
  
  
  //! conversion of a sparse matrix into a block-diagonal matrix
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {
    if (m2 <= m1)
      return;
    
    int m = m2 - m1;
    int nblock = 0, nrow = 0;
    int irow = 0, first_row = 0;
    Vector<bool> UsedRow(m); UsedRow.Fill(false);
    Vector<IVect, VectFull, NewAlloc<IVect> > pattern(m);
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
                MarkUsedColumnsSparseMatrix(A, irow, m, m1, UsedRow, row_found);
	      }
	    
	    new_row = row_found;
	  }
	
	// specifying the new block
	nrow += all_row.GetM();
	Sort(all_row.GetM(), all_row);
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
    
    // initialization of the block-diagonal matrix
    B.Clear();
    B.SetPattern(pattern);
    pattern.Clear(); UsedRow.Clear();
    
    FillBlockDiagonal(A, B, m1, m2);
  }
  

  //! conversion from any sparse matrix to a block-diagonal one
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void ConvertToBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			      Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  
  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			 Matrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    FillBlockDiagonal(A, B, 0, A.GetM());
  }

  template<class T0, class Prop0, class Storage0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
			 Matrix<T1, Prop1, Storage1, Allocator1>& B, int m1, int m2)
  {
    int m = m2-m1;
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
  
    
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B,
                         int m1, int m2)
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
  

  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Storage1, class Allocator1>
  void FillBlockDiagonal(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
			 Matrix<complex<T1>, Prop1, Storage1, Allocator1>& B,
                         int m1, int m2)
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


  //! conversion from symmetric sparse matrix to block-diagonal symmetric matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! conversion from unsymmetric sparse matrix to block-diagonal matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }


  //! conversion from symmetric sparse matrix to block-diagonal symmetric matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRowSym, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! conversion from unsymmetric sparse matrix to block-diagonal matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	    Matrix<T1, Prop1, BlockDiagRow, Allocator1>& B)
  {
    ConvertToBlockDiagonal(A, B, 0, A.GetM());
  }
  
  
  //! conversion from block-diagonal symmetric matrix to symmetric sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }


  //! conversion from block-diagonal matrix to unsymmetric sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }


  //! conversion from block-diagonal symmetric matrix to symmetric sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRowSym, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }


  //! conversion from block-diagonal matrix to unsymmetric sparse matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, BlockDiagRow, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B)
  {
    A.ConvertToSparse(B);
  }

  
  //! X is replaced by A X
  template<class T1, class Prop, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Prop, BlockDiagRow, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X)
  {
    A.Mlt_(X);
  }
  
  
  //! X is replaced by A X
  template<class T1, class Allocator1, class T2, class Allocator2>
  void Mlt(const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
           Vector<T2, VectFull, Allocator2>& X)
  {
    A.Mlt_(X);
  }
  
  
  //! computes C = A B
  template<class T0, class Allocator0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Mlt(const Matrix<T0, Symmetric, BlockDiagRowSym, Allocator0>& A,
           const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& B,
           Matrix<T2, Symmetric, BlockDiagRowSym, Allocator2>& C)
  {
    A.Mlt_(B, C);
  }
  
  
  //! computation of B = beta*B + alpha*A*X
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2, class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexZero(one);
    T4 czero; SetComplexZero(czero);
    if (beta == zero)
      B.Fill(czero);
    else if (beta != one)
      Mlt(beta, B);
    
    A.MltAdd_(alpha, SeldonNoTrans, X, B);
  }
  
  
  //! computation of B = beta*B + alpha*A*X
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexZero(one);
    T4 czero; SetComplexZero(czero);
    if (beta == zero)
      B.Fill(czero);
    else if (beta != one)
      Mlt(beta, B);
    
    A.MltAdd_(alpha, SeldonNoTrans, X, B);
  }

  
  //! computes B = beta B + alpha A^T X
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonTrans& trans,
		    const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexZero(one);
    T4 czero; SetComplexZero(czero);
    if (beta == zero)
      B.Fill(czero);
    else if (beta != one)
      Mlt(beta, B);
    
    A.MltAdd_(alpha, SeldonTrans, X, B);
  }


  //! computes B = beta B + alpha A^T X
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonTrans& trans,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexZero(one);
    T4 czero; SetComplexZero(czero);
    if (beta == zero)
      B.Fill(czero);
    else if (beta != one)
      Mlt(beta, B);
    
    A.MltAdd_(alpha, SeldonTrans, X, B);
  }
  

  //! computes B = beta B + alpha A^H X
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonConjTrans& trans,
		    const Matrix<T1, General, BlockDiagRow, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexZero(one);
    T4 czero; SetComplexZero(czero);
    if (beta == zero)
      B.Fill(czero);
    else if (beta != one)
      Mlt(beta, B);
    
    A.MltAdd_(alpha, SeldonConjTrans, X, B);
  }


  //! computes B = beta B + alpha A^H X
  template<class T0, class T1, class Allocator1, class T2, class Allocator2,
           class T3, class T4, class Allocator4>
  void MltAddVector(const T0& alpha, const class_SeldonConjTrans& trans,
		    const Matrix<T1, Symmetric, BlockDiagRowSym, Allocator1>& A,
		    const Vector<T2, VectFull, Allocator2>& X,
		    const T3& beta, Vector<T4, VectFull, Allocator4>& B)
  {
    T3 zero; SetComplexZero(zero);
    T3 one; SetComplexZero(one);
    T4 czero; SetComplexZero(czero);
    if (beta == zero)
      B.Fill(czero);
    else if (beta != one)
      Mlt(beta, B);
    
    A.MltAdd_(alpha, SeldonConjTrans, X, B);
  }

    
  //! we replace A by its inverse
  template<class T, class Allocator>
  void GetInverse(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    A.GetInverse_();
  }

  
  //! we replace A by its inverse
  template<class T, class Allocator>
  void GetInverse(Matrix<T, General, BlockDiagRow, Allocator>& A)
  {
    A.GetInverse_();
  }
  
  
  //! A is replaced by its Cholesky factorisation
  template<class T, class Allocator>
  void GetCholesky(Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    T one; SetComplexOne(one);
    for (int i = 0; i < A.diag_value_.GetM(); i++)
      A.diag_value_(i) = one/sqrt(A.diag_value_(i));

    for (int i = 0; i < A.mat22_.GetM(); i++)
      GetCholesky(A.mat22_(i));
    
    for (int i = 0; i < A.mat33_.GetM(); i++)
      GetCholesky(A.mat33_(i));

    for (int i = 0; i < A.mat44_.GetM(); i++)
      GetCholesky(A.mat44_(i));

    for (int i = 0; i < A.mat55_.GetM(); i++)
      GetCholesky(A.mat55_(i));

    for (int i = 0; i < A.mat66_.GetM(); i++)
      GetCholesky(A.mat66_(i));

    for (int i = 0; i < A.mat77_.GetM(); i++)
      GetCholesky(A.mat77_(i));

     for (int i = 0; i < A.mat88_.GetM(); i++)
      GetCholeskyMatrix(A.mat88_(i));

    for (int i = 0; i < A.mat99_.GetM(); i++)
      GetCholeskyMatrix(A.mat99_(i));

    for (int i = 0; i < A.mat1010_.GetM(); i++)
      GetCholeskyMatrix(A.mat1010_(i));

    for (int i = 0; i < A.mat1111_.GetM(); i++)
      GetCholeskyMatrix(A.mat1111_(i));
    
    for (int i = 0; i < A.mat1212_.GetM(); i++)
      GetCholeskyMatrix(A.mat1212_(i));

    for (int i = 0; i < A.mat1313_.GetM(); i++)
      GetCholeskyMatrix(A.mat1313_(i));

    for (int i = 0; i < A.mat1414_.GetM(); i++)
      GetCholeskyMatrix(A.mat1414_(i));

    /* for (int i = 0; i < A.mat1515_.GetM(); i++)
      GetCholesky(A.mat1515_(i));

    for (int i = 0; i < A.mat1616_.GetM(); i++)
      GetCholesky(A.mat1616_(i));

    for (int i = 0; i < A.mat1717_.GetM(); i++)
      GetCholesky(A.mat1717_(i));

    for (int i = 0; i < A.mat1818_.GetM(); i++)
      GetCholesky(A.mat1818_(i));

    for (int i = 0; i < A.mat1919_.GetM(); i++)
      GetCholesky(A.mat1919_(i));

    for (int i = 0; i < A.mat2020_.GetM(); i++)
    GetCholesky(A.mat2020_(i)); */

    for (int i = 0; i < A.all_mat_.GetM(); i++)
      GetCholesky(A.all_mat_(i));
  }
  
  
  //! solves L x = b or L^T x = b assuming that GetCholesky has been called
  template<class T, class Allocator, class T2, class Allocator2>
  void SolveCholesky(const SeldonTranspose& trans,
                     const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                     Vector<T2, VectFull, Allocator2>& X)
  {
    for (int i = 0; i < A.diag_value_.GetM(); i++)
      X(A.row_num_(i)) *= A.diag_value_(i);
    
    int nb = A.diag_value_.GetM();
    
    {
      TinyVector<T2, 2> x;
      for (int i = 0; i < A.mat22_.GetM(); i++)
	{
	  x(0) = X(A.row_num_(nb));
	  x(1) = X(A.row_num_(nb+1));
	  
	  SolveCholesky(trans, A.mat22_(i), x);

	  X(A.row_num_(nb)) = x(0);
	  X(A.row_num_(nb+1)) = x(1);
	  nb += 2;
	}
    }

    {
      TinyVector<T2, 3> x;
      for (int i = 0; i < A.mat33_.GetM(); i++)
	{
	  x(0) = X(A.row_num_(nb));
	  x(1) = X(A.row_num_(nb+1));
	  x(2) = X(A.row_num_(nb+2));
	  
	  SolveCholesky(trans, A.mat33_(i), x);

	  X(A.row_num_(nb)) = x(0);
	  X(A.row_num_(nb+1)) = x(1);
	  X(A.row_num_(nb+2)) = x(2);
	  nb += 3;
	}
    }

    {
      TinyVector<T2, 4> x;
      for (int i = 0; i < A.mat44_.GetM(); i++)
	{
	  x(0) = X(A.row_num_(nb));
	  x(1) = X(A.row_num_(nb+1));
	  x(2) = X(A.row_num_(nb+2));
	  x(3) = X(A.row_num_(nb+3));
	  
	  SolveCholesky(trans, A.mat44_(i), x);

	  X(A.row_num_(nb)) = x(0);
	  X(A.row_num_(nb+1)) = x(1);
	  X(A.row_num_(nb+2)) = x(2);
	  X(A.row_num_(nb+3)) = x(3);
	  nb += 4;
	}
    }

    {
      TinyVector<T2, 5> x;
      for (int i = 0; i < A.mat55_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat55_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 5;
	}
    }

    {
      TinyVector<T2, 6> x;
      for (int i = 0; i < A.mat66_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat66_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 6;
	}
    }

    {
      TinyVector<T2, 7> x;
      for (int i = 0; i < A.mat77_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat77_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 7;
	}
    }

    {
      TinyVector<T2, 8> x;
      for (int i = 0; i < A.mat88_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat88_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 8;
	}
    }

    {
      TinyVector<T2, 9> x;
      for (int i = 0; i < A.mat99_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat99_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 9;
	}
    }

    {
      TinyVector<T2, 10> x;
      for (int i = 0; i < A.mat1010_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1010_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 10;
	}
    }

    {
      TinyVector<T2, 11> x;
      for (int i = 0; i < A.mat1111_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1111_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 11;
	}
    }
    
    {
      TinyVector<T2, 12> x;
      for (int i = 0; i < A.mat1212_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1212_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 12;
	}
    }

    {
      TinyVector<T2, 13> x;
      for (int i = 0; i < A.mat1313_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1313_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 13;
	}
    }

    {
      TinyVector<T2, 14> x;
      for (int i = 0; i < A.mat1414_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1414_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 14;
	}
    }

    /*{
      TinyVector<T2, 15> x;
      for (int i = 0; i < A.mat1515_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1515_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 15;
	}
    }

    {
      TinyVector<T2, 16> x;
      for (int i = 0; i < A.mat1616_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1616_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 16;
	}
    }

    {
      TinyVector<T2, 17> x;
      for (int i = 0; i < A.mat1717_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1717_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 17;
	}
    }

    {
      TinyVector<T2, 18> x;
      for (int i = 0; i < A.mat1818_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1818_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 18;
	}
    }

    {
      TinyVector<T2, 19> x;
      for (int i = 0; i < A.mat1919_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat1919_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 19;
	}
    }

    {
      TinyVector<T2, 20> x;
      for (int i = 0; i < A.mat2020_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  SolveCholesky(trans, A.mat2020_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 20;
	}
    } */
    
    Vector<T> xs(A.size_max_blok_), x;
    for (int i = 0; i < A.all_mat_.GetM(); i++)
      {
	int size_mat = A.all_mat_(i).GetM();
	x.SetData(size_mat, xs.GetData());
        for (int j = 0; j < size_mat; j++)
          x(j) = X(A.row_num_(nb+j));
          	
	SolveCholesky(trans, A.all_mat_(i), x);
	
	for (int j = 0; j < size_mat; j++)
	  X(A.row_num_(nb+j)) = x(j);
	
        x.Nullify();
	nb += size_mat;
      }
  }


  //! computes y = L x or y = L^T x assuming that GetCholesky has been called
  template<class T, class Allocator, class T2, class Allocator2>
  void MltCholesky(const SeldonTranspose& trans,
                   const Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A,
                   Vector<T2, VectFull, Allocator2>& X)
  {
    for (int i = 0; i < A.diag_value_.GetM(); i++)
      X(A.row_num_(i)) /= A.diag_value_(i);
    
    int nb = A.diag_value_.GetM();
    
    {
      TinyVector<T2, 2> x;
      for (int i = 0; i < A.mat22_.GetM(); i++)
	{
	  x(0) = X(A.row_num_(nb));
	  x(1) = X(A.row_num_(nb+1));
	  
	  MltCholesky(trans, A.mat22_(i), x);

	  X(A.row_num_(nb)) = x(0);
	  X(A.row_num_(nb+1)) = x(1);
	  nb += 2;
	}
    }

    {
      TinyVector<T2, 3> x;
      for (int i = 0; i < A.mat33_.GetM(); i++)
	{
	  x(0) = X(A.row_num_(nb));
	  x(1) = X(A.row_num_(nb+1));
	  x(2) = X(A.row_num_(nb+2));
	  
	  MltCholesky(trans, A.mat33_(i), x);

	  X(A.row_num_(nb)) = x(0);
	  X(A.row_num_(nb+1)) = x(1);
	  X(A.row_num_(nb+2)) = x(2);
	  nb += 3;
	}
    }

    {
      TinyVector<T2, 4> x;
      for (int i = 0; i < A.mat44_.GetM(); i++)
	{
	  x(0) = X(A.row_num_(nb));
	  x(1) = X(A.row_num_(nb+1));
	  x(2) = X(A.row_num_(nb+2));
	  x(3) = X(A.row_num_(nb+3));
	  
	  MltCholesky(trans, A.mat44_(i), x);

	  X(A.row_num_(nb)) = x(0);
	  X(A.row_num_(nb+1)) = x(1);
	  X(A.row_num_(nb+2)) = x(2);
	  X(A.row_num_(nb+3)) = x(3);
	  nb += 4;
	}
    }


    {
      TinyVector<T2, 5> x;
      for (int i = 0; i < A.mat55_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat55_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 5;
	}
    }

    {
      TinyVector<T2, 6> x;
      for (int i = 0; i < A.mat66_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat66_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 6;
	}
    }

    {
      TinyVector<T2, 7> x;
      for (int i = 0; i < A.mat77_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat77_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 7;
	}
    }

    {
      TinyVector<T2, 8> x;
      for (int i = 0; i < A.mat88_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat88_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 8;
	}
    }

    {
      TinyVector<T2, 9> x;
      for (int i = 0; i < A.mat99_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat99_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 9;
	}
    }

    {
      TinyVector<T2, 10> x;
      for (int i = 0; i < A.mat1010_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1010_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 10;
	}
    }

    {
      TinyVector<T2, 11> x;
      for (int i = 0; i < A.mat1111_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1111_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 11;
	}
        }

    {
      TinyVector<T2, 12> x;
      for (int i = 0; i < A.mat1212_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1212_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 12;
	}
    }

    {
      TinyVector<T2, 13> x;
      for (int i = 0; i < A.mat1313_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1313_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 13;
	}
    }

    {
      TinyVector<T2, 14> x;
      for (int i = 0; i < A.mat1414_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1414_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 14;
	}
    }

    /* {
      TinyVector<T2, 15> x;
      for (int i = 0; i < A.mat1515_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1515_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 15;
	}
    }

    {
      TinyVector<T2, 16> x;
      for (int i = 0; i < A.mat1616_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1616_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 16;
	}
    }

    {
      TinyVector<T2, 17> x;
      for (int i = 0; i < A.mat1717_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1717_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 17;
	}
    }

    {
      TinyVector<T2, 18> x;
      for (int i = 0; i < A.mat1818_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1818_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 18;
	}
    }

    {
      TinyVector<T2, 19> x;
      for (int i = 0; i < A.mat1919_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat1919_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 19;
	}
    }

    {
      TinyVector<T2, 20> x;
      for (int i = 0; i < A.mat2020_.GetM(); i++)
	{
          ExtractVector(X, A.row_num_, nb, x);
	  MltCholesky(trans, A.mat2020_(i), x);
          ExtractVector(x, A.row_num_, nb, X);
	  nb += 20;
	}
        } */
    
    Vector<T> xs(A.size_max_blok_), x;
    for (int i = 0; i < A.all_mat_.GetM(); i++)
      {
	int size_mat = A.all_mat_(i).GetM();
        x.SetData(size_mat, xs.GetData());
	for (int j = 0; j < size_mat; j++)
          x(j) = X(A.row_num_(nb+j));
          	
	MltCholesky(trans, A.all_mat_(i), x);
	
	for (int j = 0; j < size_mat; j++)
	  X(A.row_num_(nb+j)) = x(j);
	
        x.Nullify();
	nb += size_mat;
      }
  }


  //! A = alpha*A
  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, General, BlockDiagRow, Allocator>& A)
  {
    A *= alpha;
  }
   
  
  //! A = alpha*A
  template<class T, class Allocator, class T0>
  void MltScalar(const T0& alpha, Matrix<T, Symmetric, BlockDiagRowSym, Allocator>& A)
  {
    A *= alpha;
  }
  
  
  //! computes C = A B with Seldon matrices (avoids too many template recursions)
  template<class T, int m>
  void MltMatrix(const TinyMatrix<T, Symmetric, m, m>& A, const TinyMatrix<T, Symmetric, m, m>& B,
                 TinyMatrix<T, Symmetric, m, m>& C)
  {
    Matrix<T, Symmetric, RowSymPacked> Ad, Bd, Cd;
    Ad.SetData(m, m, const_cast<T*>(&A(0, 0)));
    Bd.SetData(m, m, const_cast<T*>(&B(0, 0)));
    Cd.SetData(m, m, &C(0, 0));
    
    Mlt(Ad, Bd, Cd);
    
    Ad.Nullify(); Bd.Nullify(); Cd.Nullify();
  }


  //! replaces A by its inverse with Seldon matrix (avoids too many template recursions)
  template<class T, int m>
  void GetInverseMatrix(TinyMatrix<T, General, m, m>& A)
  {
    Matrix<T, General, RowMajor> Ad;
    Ad.SetData(m, m, &A(0, 0));
    
    GetInverse(Ad);
    
    Ad.Nullify();
  }

  
  //! replaces A by its inverse with Seldon matrix (avoids too many template recursions)
  template<class T, int m>
  void GetInverseMatrix(TinyMatrix<T, Symmetric, m, m>& A)
  {
    Matrix<T, Symmetric, RowSymPacked> Ad;
    Ad.SetData(m, m, &A(0, 0));
    
    GetInverse(Ad);
    
    Ad.Nullify();
  }


  //! replaces A by its Cholesky factorisation with Seldon matrix
  //! (avoids too many template recursions)
  template<class T, int m>
  void GetCholeskyMatrix(TinyMatrix<T, Symmetric, m, m>& A)
  {
    Matrix<T, Symmetric, RowSymPacked> Ad;
    Ad.SetData(m, m, &A(0, 0));
    
    GetCholesky(Ad);
    
    Ad.Nullify();
  }
  
}

#define SELDON_FILE_MATRIX_BLOCK_DIAGONAL_CXX
#endif
