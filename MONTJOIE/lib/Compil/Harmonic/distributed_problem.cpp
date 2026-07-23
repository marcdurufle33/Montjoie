#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Harmonic/DistributedProblem.cxx"
#endif

namespace Montjoie
{
  
  SELDON_EXTERN template class DistributedProblem<Dimension2>;
  SELDON_EXTERN template class DistributedProblem<Dimension3>;

  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Real_wp, General, RowSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Real_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Real_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, General, RowSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, Symmetric, RowSymSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, General, ArrayRowSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse>&);

  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, General, RowComplexSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse>&);

  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(Matrix<Real_wp, Symmetric, DiagonalRow>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(Matrix<Complex_wp, Symmetric, DiagonalRow>&);

  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Real_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, General, BlockDiagRow>&);
  SELDON_EXTERN template void DistributedProblem_Base::InitDistributedMatrix(DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym>&);

  SELDON_EXTERN template void DistributedProblem_Base::AddDomains(Vector<int>&, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::AddDomains(VectReal_wp&, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::AddDomains(VectComplex_wp&, int) const;

  SELDON_EXTERN template void DistributedProblem_Base::AssembleDirichlet(Vector<int>&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::AssembleDirichlet(VectReal_wp&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::AssembleDirichlet(VectComplex_wp&, int, bool) const;
  
  SELDON_EXTERN template void DistributedProblem_Base::ExchangeDomains(VectReal_wp&, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::ExchangeDomains(VectComplex_wp&, int) const;

  SELDON_EXTERN template void DistributedProblem_Base::ExchangeRelaxDomains(VectReal_wp&, const Real_wp&, int, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::ExchangeRelaxDomains(VectComplex_wp&, const Real_wp&, int, int) const;

  SELDON_EXTERN template void DistributedProblem_Base::ExchangeQuadRelaxDomains(VectReal_wp&, const Real_wp&, int, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::ExchangeQuadRelaxDomains(VectComplex_wp&, const Real_wp&, int, int) const;

#ifdef SELDON_WITH_MPI
  SELDON_EXTERN template void DistributedProblem_Base::ExchangeUfaceDomains(const Vector<Real_wp>&, Vector<Vector<Real_wp> >&, Vector<Vector<int64_t> >&, Vector<Vector<Real_wp> >&, Vector<Vector<int64_t> >&, Vector<MPI_Request>&, int) const;

  SELDON_EXTERN template void DistributedProblem_Base::ExchangeUfaceDomains(const Vector<Complex_wp>&, Vector<Vector<Complex_wp> >&, Vector<Vector<int64_t> >&, Vector<Vector<Complex_wp> >&, Vector<Vector<int64_t> >&, Vector<MPI_Request>&, int) const;
  
  SELDON_EXTERN template void DistributedProblem_Base::GetUfaceDomains(Vector<Real_wp>&, Vector<Vector<Real_wp> >&, Vector<Vector<int64_t> >&, Vector<Vector<Real_wp> >& xdom, Vector<Vector<int64_t> >& xdom_tmp, Vector<MPI_Request>&, int) const;

  SELDON_EXTERN template void DistributedProblem_Base::GetUfaceDomains(Vector<Complex_wp>&, Vector<Vector<Complex_wp> >&, Vector<Vector<int64_t> >&, Vector<Vector<Complex_wp> >& xdom, Vector<Vector<int64_t> >& xdom_tmp, Vector<MPI_Request>&, int) const;


  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<bool>&, const MPI_Op&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<int>&, const MPI_Op&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<Real_wp>&, const MPI_Op&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<Complex_wp>&, const MPI_Op&, int, bool) const;
  
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<bool>&, const MPI_Op&, const IVect&, const Vector<IVect>&, int, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<int>&, const MPI_Op&, const IVect&, const Vector<IVect>&, int, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<Real_wp>&, const MPI_Op&, const IVect&, const Vector<IVect>&, int, int) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDistributedVector(Vector<Complex_wp>&, const MPI_Op&, const IVect&, const Vector<IVect>&, int, int) const;

  SELDON_EXTERN template void DistributedProblem_Base::ReduceDirichlet(Vector<int>&, const MPI_Op&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDirichlet(VectReal_wp&, const MPI_Op&, int, bool) const;
  SELDON_EXTERN template void DistributedProblem_Base::ReduceDirichlet(VectComplex_wp&, const MPI_Op&, int, bool) const;

  SELDON_EXTERN template void DistributedProblem<Dimension2>::ReduceDistributedVectorFace(VectReal_wp&, const MPI_Op&, int) const;
  SELDON_EXTERN template void DistributedProblem<Dimension2>::ReduceDistributedVectorFace(VectComplex_wp&, const MPI_Op&, int) const;
  SELDON_EXTERN template void DistributedProblem<Dimension3>::ReduceDistributedVectorFace(VectReal_wp&, const MPI_Op&, int) const;
  SELDON_EXTERN template void DistributedProblem<Dimension3>::ReduceDistributedVectorFace(VectComplex_wp&, const MPI_Op&, int) const;

#endif
  
}
