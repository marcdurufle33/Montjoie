#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/BoundaryConditionHarmonic.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class VirtualImpedanceFunction<Real_wp, Dimension2>;
  SELDON_EXTERN template class VirtualImpedanceFunction<Complex_wp, Dimension2>;

  SELDON_EXTERN template class VirtualImpedanceFunction<Real_wp, Dimension3>;
  SELDON_EXTERN template class VirtualImpedanceFunction<Complex_wp, Dimension3>;

  SELDON_EXTERN template class ImpedanceFunction_Base<Real_wp, Dimension2>;
  SELDON_EXTERN template class ImpedanceFunction_Base<Complex_wp, Dimension2>;

  SELDON_EXTERN template class ImpedanceFunction_Base<Real_wp, Dimension3>;
  SELDON_EXTERN template class ImpedanceFunction_Base<Complex_wp, Dimension3>;

  SELDON_EXTERN template void VarBoundaryCondition_Base::ImposeNullDirichletCondition(VectReal_wp&) const;
  SELDON_EXTERN template void VarBoundaryCondition_Base::ImposeNullDirichletCondition(VectComplex_wp&) const;

  SELDON_EXTERN template Real_wp VarBoundaryCondition_Base::ComputeDirichletCoef(VirtualMatrix<Real_wp>&);
  SELDON_EXTERN template Real_wp VarBoundaryCondition_Base::ComputeDirichletCoef(VirtualMatrix<Complex_wp>&);


  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Real_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Real_wp, Symmetric, DiagonalRow>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Complex_wp, Symmetric, DiagonalRow>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Complex_wp, General, BlockDiagRow>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(FemMatrixFreeClass_Base<Real_wp>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(FemMatrixFreeClass_Base<Complex_wp>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(Matrix<Complex_wp, General, ArrayRowSparse>&);

#ifdef SELDON_WITH_MPI    
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void VarBoundaryCondition_Base::SetPeriodicCondition(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
#endif


  SELDON_EXTERN template void VarBoundaryCondition_Base::ApplyPeriodicCondition(VectReal_wp&) const;
  SELDON_EXTERN template void VarBoundaryCondition_Base::ApplyPeriodicCondition(VectComplex_wp&) const;

  SELDON_EXTERN template void VarBoundaryCondition_Base::
  AddMatrixImpedanceBoundary(const Real_wp&, const Vector<int>&, int, const GlobalGenericMatrix<Real_wp>&,
                             VirtualMatrix<Real_wp>&, int, int, ImpedanceFunction_Base<Real_wp, Dimension2>&,
                             bool, bool, const VarProblem<Dimension2>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::
  AddMatrixImpedanceBoundary(const Complex_wp&, const Vector<int>&, int, const GlobalGenericMatrix<Complex_wp>&,
                             VirtualMatrix<Complex_wp>&, int, int, ImpedanceFunction_Base<Complex_wp, Dimension2>&,
                             bool, bool, const VarProblem<Dimension2>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::
  AddMatrixImpedanceBoundary(const Real_wp&, const Vector<int>&, int, const GlobalGenericMatrix<Real_wp>&,
                             VirtualMatrix<Real_wp>&, int, int, ImpedanceFunction_Base<Real_wp, Dimension3>&,
                             bool, bool, const VarProblem<Dimension3>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::
  AddMatrixImpedanceBoundary(const Complex_wp&, const Vector<int>&, int, const GlobalGenericMatrix<Complex_wp>&,
                             VirtualMatrix<Complex_wp>&, int, int, ImpedanceFunction_Base<Complex_wp, Dimension3>&,
                             bool, bool, const VarProblem<Dimension3>&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::GetDampingFactorPML(Mesh<Dimension2>&, int , int, const R2&, R2_Complex_wp&, R2_Complex_wp&);
  
  SELDON_EXTERN template void VarBoundaryCondition_Base::GetDampingFactorPML(Mesh<Dimension3>&, int , int, const R3&, R3_Complex_wp&, R3_Complex_wp&);

  SELDON_EXTERN template void VarBoundaryCondition_Base::TreatDirichletCondition(const VarProblem<Dimension2>& var);
  SELDON_EXTERN template void VarBoundaryCondition_Base::TreatDirichletCondition(const VarProblem<Dimension3>& var);
  
  SELDON_EXTERN template class VarBoundaryCondition_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarBoundaryCondition_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class VarBoundaryCondition_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarBoundaryCondition_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template class VarBoundaryCondition<Real_wp, Dimension2>;
  SELDON_EXTERN template class VarBoundaryCondition<Complex_wp, Dimension2>;

  SELDON_EXTERN template class VarBoundaryCondition<Real_wp, Dimension3>;
  SELDON_EXTERN template class VarBoundaryCondition<Complex_wp, Dimension3>;

  SELDON_EXTERN template void VarBoundaryCondition<Real_wp, Dimension2>::AddBoundaryConditionGen(VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, int, int);
  SELDON_EXTERN template void VarBoundaryCondition<Real_wp, Dimension2>::AddBoundaryConditionGen(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, int, int);
  SELDON_EXTERN template void VarBoundaryCondition<Complex_wp, Dimension2>::AddBoundaryConditionGen(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, int, int);
  SELDON_EXTERN template void VarBoundaryCondition<Complex_wp, Dimension2>::AddBoundaryConditionGen(VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, int, int);

  SELDON_EXTERN template void VarBoundaryCondition<Real_wp, Dimension3>::AddBoundaryConditionGen(VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, int, int);
  SELDON_EXTERN template void VarBoundaryCondition<Real_wp, Dimension3>::AddBoundaryConditionGen(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, int, int);
  SELDON_EXTERN template void VarBoundaryCondition<Complex_wp, Dimension3>::AddBoundaryConditionGen(VirtualMatrix<Complex_wp>&, const GlobalGenericMatrix<Complex_wp>&, int, int);
  SELDON_EXTERN template void VarBoundaryCondition<Complex_wp, Dimension3>::AddBoundaryConditionGen(VirtualMatrix<Real_wp>&, const GlobalGenericMatrix<Real_wp>&, int, int);
  
}
