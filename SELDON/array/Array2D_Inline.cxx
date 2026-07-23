// Copyright (C) 2022 Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_ARRAY2D_INLINE_CXX

#include "Array2D.hxx"

namespace Seldon
{
  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Array2D<T, Allocator>::~Array2D()
  {
    Clear();
  }


  /*****************
   * BASIC METHODS *
   *****************/


  //! Returns the length in dimension #1.
  /*!
    \return The length in dimension #1.
  */
  template <class T, class Allocator>
  inline int Array2D<T, Allocator>::GetM() const
  {
    return length1_;
  }


  //! Returns the length in dimension #2.
  /*!
    \return The length in dimension #2.
  */
  template <class T, class Allocator>
  inline int Array2D<T, Allocator>::GetN() const
  {
    return length2_;
  }


  //! Returns the number of elements in the 2D array.
  /*!
    Returns the number of elements stored by the 2D array, i.e.
    the product of the lengths in the two dimensions.
    \return The number of elements in the 2D array.
  */
  template <class T, class Allocator>
  inline long Array2D<T, Allocator>::GetSize() const
  {
    return length1_ * length2_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory by
    the array, i.e. the product of lengths in the two
    dimensions.
    \return The number of elements stored in the array.
  */
  template <class T, class Allocator>
  inline long Array2D<T, Allocator>::GetDataSize() const
  {
    return length1_ * length2_;
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_' which stores the
    values.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  inline typename Array2D<T, Allocator>::pointer Array2D<T, Allocator>
  ::GetData() const
  {
    return data_;
  }


  //! Returns a pointer to an element of data array.
  /*!
    Returns a pointer to an element of data array.
    \param i index along dimension #1.
    \param j index along dimension #2.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  inline typename Array2D<T, Allocator>::pointer Array2D<T, Allocator>
  ::GetDataPointer(int i, int j) const
  {
    return data_ + i * length2_ + j;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \return Element (i, j) of the 2D array.
  */
  template <class T, class Allocator>
  inline typename Array2D<T, Allocator>::reference
  Array2D<T, Allocator>::operator() (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, length1_, length2_, "Array2D");
#endif

    return data_[i * length2_ + j];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \return Element (i, j) of the 2D array.
  */
  template <class T, class Allocator>
  inline typename Array2D<T, Allocator>::const_reference
  Array2D<T, Allocator>::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, length1_, length2_, "Array2D");
#endif

    return data_[i*length2_ + j];
  }

  //! Duplicates a 2D array (assignment operator).
  /*!
    \param A 2D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Array2D<T, Allocator>& Array2D<T, Allocator>::operator=
  (const Array2D<T, Allocator>& A)
  {
    Copy(A);

    return *this;
  }

  //! Duplicates a 2D array.
  /*!
    \param A 2D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Array2D<T, Allocator>::Copy(const Array2D<T, Allocator>& A)
  {
    Reallocate(A.GetM(), A.GetN());

    Allocator::memorycpy(data_, A.GetData(), GetDataSize());
  }
  
} // namespace Seldon.

#define SELDON_FILE_ARRAY2D_INLINE_CXX
#endif
