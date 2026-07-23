#ifndef MONTJOIE_FILE_ASSEMBLE_MATRIX_INLINE_CXX

namespace Montjoie
{

  //! destructor
  template<class T>
  inline CondensationBlockSolver_Base<T>::~CondensationBlockSolver_Base()
  {
  }


  //! modifies the elementary matrix (applies static condensation if needed)
  template<class T>
  void CondensationBlockSolver_Base<T>
  ::ModifyElementaryMatrix(int i, IVect& num_ddl, VirtualMatrix<T>& mat_interac,
			   const GlobalGenericMatrix<T>& nat_mat)
  {
  }
  

  //! sets the treatment of inside stiffness
  template<class T>
  inline void CondensationBlockSolver_Base<T>::SetTreatmentStiffnessInside(bool t)
  {
    this->treatment_stiff_inside = t;
  }


  //! returns true if the internal stiffness matrix must be computed 
  template<class T>
  inline bool CondensationBlockSolver_Base<T>::TreatInsideStiffness() const
  {
    return this->treatment_stiff_inside;
  }


  //! changes the element number of the condensed element and global element
  template<class T>
  inline void CondensationBlockSolver_Base<T>::SetElementNumber(int n, int ng)
  {
    num_elem_condensed = n;
    num_elem_global = ng;
  }


  //! returns the element number of the condensed element
  template<class T>
  inline int CondensationBlockSolver_Base<T>::GetCondensedElementNumber() const
  {
    return num_elem_condensed;
  }


  //! returns the element number of the global element
  template<class T>
  inline int CondensationBlockSolver_Base<T>::GetGlobalElementNumber() const
  {
    return num_elem_global;
  }


  //! returns the number of condensed elements
  template<class T>
  inline int CondensationBlockSolver_Base<T>::GetNbCondensedElt() const
  {
    return nb_elt_condensed;
  }


  //! sets the number of condensed elements
  template<class T>
  inline void CondensationBlockSolver_Base<T>::SetNbCondensedElt(int n)
  {
    nb_elt_condensed = n;
  }

  
  //! returns mass coefficient
  template<class T>
  const T& GlobalGenericMatrix<T>::GetCoefMass() const
  {
    return coef_mass;
  }
  
  
  //! returns stiffness coefficient
  template<class T>
  const T& GlobalGenericMatrix<T>::GetCoefStiffness() const
  {
    return coef_stiff;
  }
  
  
  //! returns damping coefficient
  template<class T>
  const T& GlobalGenericMatrix<T>::GetCoefDamping() const
  {
    return coef_sigma;
  }
  
  
  //! changes mass coefficient
  template<class T>
  void GlobalGenericMatrix<T>::SetCoefMass(const T& alpha)
  {
    coef_mass = alpha;
  }
  
  
  //! changes stiffness coefficient
  template<class T>
  void GlobalGenericMatrix<T>::SetCoefStiffness(const T& alpha)
  {
    coef_stiff = alpha;
  }
  
  
  //! changes damping coefficient
  template<class T>
  void GlobalGenericMatrix<T>::SetCoefDamping(const T& alpha)
  {
    coef_sigma = alpha;
  }


  /****************************
   * MatrixVectorProductLevel *
   ****************************/
  
  
  //! sets to the current level
  inline void MatrixVectorProductLevel::SetLevel()
  {
    SetLevel(current_level);
  }
  
  
  //! returns the element numbers for each level
  inline Vector<IVect>& MatrixVectorProductLevel::GetLevelArray()
  {
    return *num_elem;
  }


  //! sets (shallow copy) the different levels
  inline void MatrixVectorProductLevel::SetLevelArray(Vector<IVect>& lvl)
  {
    num_elem = &lvl;
  }
  
  
  //! returns the number of elements associated with the current level
  inline int MatrixVectorProductLevel::GetNbElt() const
  {
    return GetNbElt(current_level);
  }


  //! returns the number of elements associated with a given level
  inline int MatrixVectorProductLevel::GetNbElt(int lvl) const
  {
    switch (lvl)
      {
      case ALL_LEVELS : return nb_elt;
      case LVL_PML: return nb_elt - nb_elt_outside_PML;
      case LVL_NOPML : return nb_elt_outside_PML;
      default: return (*num_elem)(lvl).GetM();
      }  
  }
  

  //! returns the number of the element i of the current level
  inline int MatrixVectorProductLevel::GetElementNumber(int i) const
  {
    switch (current_level)
      {
      case ALL_LEVELS : return i;
      case LVL_PML: return nb_elt_outside_PML + i;
      case LVL_NOPML : return i;
      default: return (*num_elem)(current_level)(i);
      }  
  }
  

  //! returns the local element number
  inline int MatrixVectorProductLevel::GetLocalElementNumber() const
  {
    return nelem;
  }
  
}

#define MONTJOIE_FILE_ASSEMBLE_MATRIX_INLINE_CXX
#endif

