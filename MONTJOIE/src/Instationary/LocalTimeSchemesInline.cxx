#ifndef MONTJOIE_FILE_LOCAL_TIME_SCHEMES_INLINE_CXX

namespace Montjoie
{

  //! strategy used to split the domain into a coarse region and a fine region
  inline int LocalTimeScheme_Base::GetSplittingType() const
  {
    return type_splitting;
  }


  //! returns the pivot time step used to split the domain into two regions
  inline Real_wp LocalTimeScheme_Base::GetSplittingDt() const
  {
    return dt_splitting;
  }


  //! return the file containing the element numbers
  inline string LocalTimeScheme_Base::GetNumberElementFileName() const
  {
    return file_with_number_element;
  }


  //! sets the list of elements defining the fine region
  inline void LocalTimeScheme_Base::SetImplicitElement(const IVect& num)
  {
    num_implicit = num;
  }


  inline int LocalOdeScheme_Base::GetNbLocalDof() const
  {
    return nb_dof_uv_local;
  }


  inline int LocalOdeScheme_Base::GetNbLocalDofLambda() const
  {
    return nb_dof_L_local;
  }


  inline int LocalOdeScheme_Base::GetNbDofProche() const
  {
    return nb_dof_uv_proche;
  }


  inline int LocalOdeScheme_Base::GetVolumeLevel() const
  {
    return num_level_vol;
  }
  

  inline const Vector<int>& LocalOdeScheme_Base::GetLocalDof() const
  {
    return dof_local_uv;
  }


  inline const Vector<int>& LocalOdeScheme_Base::GetDofProche() const
  {
    return dof_proche_uv;
  }
  
}

#define MONTJOIE_FILE_LOCAL_TIME_SCHEMES_INLINE_CXX
#endif

