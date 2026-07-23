#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Computation/FemMatrixFreeClass.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class ExtrapolVariablesProduct_Base<Real_wp>;
  SELDON_EXTERN template class ExtrapolVariablesProduct_Base<Complex_wp>;

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Real_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension2>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Real_wp&, Real_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Real_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension3>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Real_wp&, Real_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Real_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension2>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Complex_wp&, Complex_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Real_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension3>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Complex_wp&, Complex_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Complex_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension2>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Real_wp&, Real_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Complex_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension3>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Real_wp&, Real_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Complex_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension2>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Complex_wp&, Complex_wp&);

  SELDON_EXTERN template void ExtrapolVariablesProduct_Base<Complex_wp>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension3>&, int, int, int, int, int&, int&,
			      bool&, int&, int&, int&, int&, Complex_wp&, Complex_wp&);

  SELDON_EXTERN template class FemMatrixFreeClass_Base<Real_wp>;
  SELDON_EXTERN template class FemMatrixFreeClass_Base<Complex_wp>;

  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::ApplyRightScaling(const VectReal_wp&, VectReal_wp& C2, VectReal_wp& B, VectReal_wp&) const;
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::ApplyLeftScaling(const VectReal_wp&, VectReal_wp& C2, VectReal_wp& B, VectReal_wp&) const;

  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::ApplyRightScaling(const VectComplex_wp&, VectComplex_wp& C2, VectComplex_wp& B, VectComplex_wp&) const;
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::ApplyLeftScaling(const VectComplex_wp&, VectComplex_wp& C2, VectComplex_wp& B, VectComplex_wp&) const;

  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::ApplyRightScaling(const VectComplex_wp&, VectComplex_wp& C2, VectComplex_wp& B, VectComplex_wp&) const;
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::ApplyLeftScaling(const VectComplex_wp&, VectComplex_wp& C2, VectComplex_wp& B, VectComplex_wp&) const;
  
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::SetDirichletCondition(DistributedMatrix<Real_wp, General, ArrayRowSparse>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::SetDirichletCondition(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&, int, int, bool);

  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::SetDirichletCondition(Matrix<Real_wp, Symmetric, DiagonalRow>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(Matrix<Complex_wp, Symmetric, DiagonalRow>&, int, int, bool);

  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::SetDirichletCondition(DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&, int, int, bool);

  SELDON_EXTERN template void FemMatrixFreeClass_Base<Real_wp>::SetDirichletCondition(DistributedMatrix<Real_wp, General, BlockDiagRow>&, int, int, bool);
  SELDON_EXTERN template void FemMatrixFreeClass_Base<Complex_wp>::SetDirichletCondition(DistributedMatrix<Complex_wp, General, BlockDiagRow>&, int, int, bool);

  SELDON_EXTERN template void CopyReal(const FemMatrixFreeClass_Base<Complex_wp>&, FemMatrixFreeClass_Base<Real_wp>&);
  SELDON_EXTERN template void CopyReal(const FemMatrixFreeClass_Base<Real_wp>&, FemMatrixFreeClass_Base<Real_wp>&);
  
}

