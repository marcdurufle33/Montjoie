#ifndef MONTJOIE_FILE_FLUID_STRUCTURE_INTERACTION_INLINE_CXX

namespace Montjoie
{

  //! returns the number of unknowns of the evolution system
  //! (ie size of X when solving dX/dt = g(t, X))
  template<class Dimension>
  inline int FluidStructureInteraction_Base<Dimension>::GetNbDof() const
  {
    return this->nodl;
  }
  
  
  //! returns frequency of the time source
  template<class TypeEqFluid, class TypeEqSolid>
  inline Real_wp FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::GetFrequency() const
  {
    return var_fluid.var_harmonic.GetFrequency();
  }
  
  
  //! returns time step
  template<class TypeEqFluid, class TypeEqSolid>
  inline Real_wp FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::GetTimeStep() const
  {
    return var_fluid.GetTimeStep();
  }
  
  
  //! returns initial time
  template<class TypeEqFluid, class TypeEqSolid>
  inline Real_wp FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::GetInitialTime() const
  {
    return var_fluid.GetInitialTime();
  }
  
  
  //! returns final time
  template<class TypeEqFluid, class TypeEqSolid>
  inline Real_wp FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::GetFinalTime() const
  {
    return var_fluid.GetFinalTime();
  }
  
  
  //! not used
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::
  GiveIterate(int n, const Real_wp& t, VectReal_wp& Y)
  {
  }
    
  
  //! overwriting X by M_p^-1 X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::
  SolveMassPressure(VectReal_wp& X)
  {
    var_fluid.SolveOperatorDh(X);
  }

  
  //! overwriting X by M_v^-1 X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::SolveMassDisplacementFluid(VectReal_wp& X)
  {
    var_fluid.SolveOperatorBh(X);
  }
  
  
  //! overwriting X by M_u^-1 X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::SolveMassDisplacementSolid(VectReal_wp& X)
  {
    var_solid.SolveOperatorDh(X);
  }
  
  
  //! overwriting X by M_sigma^-1 X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::SolveMassStress(VectReal_wp& X)
  {
    var_solid.SolveOperatorBh(X);
  }
    
  
  //! computing Y = beta*Y + alpha* R_p X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::ApplyStiffnessPressure(const Real_wp& alpha, const VectReal_wp& X,
                           const Real_wp& beta, VectReal_wp& Y)
  {
    var_fluid.ApplyOperatorRhScalar(-alpha, 0.0, X, beta, Y);
  }
  
  
  //! computing Y = beta*Y + alpha* R_v X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::ApplyStiffnessDisplacementFluid(const Real_wp& alpha, const VectReal_wp& X,
                                    const Real_wp& beta, VectReal_wp& Y)
  {
    var_fluid.ApplyOperatorRhVectorial(alpha, 0.0, X, beta, Y);
  }
  
  
  //! computing Y = beta*Y + alpha* R_u X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::ApplyStiffnessDisplacementSolid(const Real_wp& alpha, const VectReal_wp& X,
                                    const Real_wp& beta, VectReal_wp& Y)
  {
    var_solid.ApplyOperatorRhScalar(-alpha, 0, X, beta, Y);
  }
  
  
  //! computing Y = beta*Y + alpha* R_sigma X
  template<class TypeEqFluid, class TypeEqSolid>
  inline void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::ApplyStiffnessStress(const Real_wp& alpha, const VectReal_wp& X,
                         const Real_wp& beta, VectReal_wp& Y)
  {
    var_solid.ApplyOperatorRhVectorial(alpha, 0, X, beta, Y);
  }
  
}

#define MONTJOIE_FILE_FLUID_STRUCTURE_INTERACTION_INLINE_CXX
#endif

