#ifndef MONTJOIE_FILE_VIBRO_ACOUSTIC_SESSION_INLINE_CXX

namespace Montjoie
{
  //! returns sound speed
  template<class TypeEqPlate, class TypeEqVol>
  inline Real_wp VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetCelerity() const
  {
    return sqrt(c2);
  }
  
  
  //! returns frequency of the source
  template<class TypeEqPlate, class TypeEqVol>
  inline Real_wp VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetFrequency() const
  {
    return var_volume.var_harmonic.GetFrequency();
  }
  
  
  //! returns order of approximation for acoustic equation
  template<class TypeEqPlate, class TypeEqVol>
  inline int VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetOrder() const
  {
    return var_volume.var_harmonic.GetMeshNumbering(0).GetOrder();
  }
  
  
  //! returns mesh size used in the volume
  template<class TypeEqPlate, class TypeEqVol>
  inline Real_wp VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetMeshSize() const
  {
    return var_volume.var_harmonic.mesh.GetMeshSize();
  }

  
  //! returns time step
  template<class TypeEqPlate, class TypeEqVol>
  inline Real_wp VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetTimeStep() const
  {
    return var_plate.GetTimeStep();
  }
  
  
  //! returns initial time of the experiment
  template<class TypeEqPlate, class TypeEqVol>
  inline Real_wp VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetInitialTime() const
  {
    return var_plate.GetInitialTime();
  }
  
  
  //! returns final time of the experiment
  template<class TypeEqPlate, class TypeEqVol>
  inline Real_wp VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetFinalTime() const
  {
    return var_plate.GetFinalTime();
  }
  
  
  //! returns the number of dofs used for Reissner-Mindlin (mecanic) equations 
  //! (only for one unknown, u for instance)
  template<class TypeEqPlate, class TypeEqVol>
  inline int VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetNbMecanicDofSurface() const
  {
    return var_plate.var_harmonic.GetMeshNumbering(0).GetNbDof();
  }
  
  
  //! returns the number of dofs for acoustic equation that are on the plate
  template<class TypeEqPlate, class TypeEqVol>
  inline int VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetNbAcousticDofSurface() const
  {
    //return ddl_vol_minus.GetM();
    return nodl_acous;
  }

  
  //! returns the size of linear system to solve
  template<class TypeEqPlate, class TypeEqVol>
  inline int VibroAcousticSession<TypeEqPlate, TypeEqVol>::GetM() const
  {
    return GetNbMecanicDofSurface();
  }
  
  
  //! returns the number of volume dofs
  template<class TypeEqPlate, class TypeEqVol>
  inline int VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetNbDofVolume() const
  {
    //return var_volume.var_harmonic.mesh.GetNbDof();
    return var_volume.GetNbScalarUnknowns();
  }
  
  
  //! returns the diagonal mass matrix of acoustic equation
  template<class TypeEqPlate, class TypeEqVol>
  inline VectReal_wp& VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  GetAcousticMass()
  {
    return const_cast<VectReal_wp&>(var_volume.Glob_mat_Dh->GetDiagonalDh()); 
  }
  
  
  //! returns the diagonal damping matrix of acoustic equation
  template<class TypeEqPlate, class TypeEqVol>
  inline VectReal_wp& VibroAcousticSession<TypeEqPlate, TypeEqVol>::
  GetAcousticDamping()
  {
    return const_cast<VectReal_wp&>(var_volume.Glob_mat_Dh->GetDiagonalSh());
  }


  template<class TypeEqPlate, class TypeEqVol>
  inline const VectReal_wp& VibroAcousticSession<TypeEqPlate, TypeEqVol>
  ::GetMecanicMass() const
  {
    return MassLambda;
  }
  
}

#define MONTJOIE_FILE_VIBRO_ACOUSTIC_SESSION_INLINE_CXX
#endif
