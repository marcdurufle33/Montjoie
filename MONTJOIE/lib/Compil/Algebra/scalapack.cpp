#include "MontjoieFlag.hxx"

#include "Share/MontjoieCommonHeader.hxx"
#include "Share/MontjoieCommonInline.hxx"

#ifdef SELDON_WITH_MPI

#include "Algebra/DistributedDenseMatrix.cxx"

#ifdef SELDON_WITH_SCALAPACK

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Algebra/Scalapack.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template class DistributedMatrix_Pointers<double, General, ColMajor>;

  SELDON_EXTERN template class DistributedMatrix<double, General, ColMajor>;
  
  SELDON_EXTERN template void GetEigenvalues(DistributedMatrix<double, General, ColMajor>&, Vector<double>&, Vector<double>&);

  SELDON_EXTERN template void DistributeVector(Vector<Real_wp>&, DistributedMatrix<Real_wp, General, ColMajor>& A);
  SELDON_EXTERN template void DistributeVector(Vector<Complex_wp>&, DistributedMatrix<Complex_wp, General, ColMajor>& A);

  SELDON_EXTERN template void GatherVector(DistributedMatrix<Real_wp, General, ColMajor>&, Vector<Real_wp>&);
  SELDON_EXTERN template void GatherVector(DistributedMatrix<Complex_wp, General, ColMajor>&, Vector<Complex_wp>&);

  SELDON_EXTERN template void SolveCholesky(DistributedMatrix<Real_wp, General, ColMajor>&, DistributedMatrix<Real_wp, General, ColMajor>&, LapackInfo&);

  SELDON_EXTERN template void GetCholesky(DistributedMatrix<Real_wp, General, ColMajor>&, LapackInfo&);

}

#endif

namespace Seldon
{
  SELDON_EXTERN template class DistributedMatrix<Real_wp, General, RowMajor>;
  
  SELDON_EXTERN template void Mlt(const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, Vector<Real_wp>&, bool);
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowMajor>&, const Vector<Real_wp>&, const Real_wp&, Vector<Real_wp>&, bool); 
  SELDON_EXTERN template void MltAdd(const Real_wp&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowMajor>&, const SeldonTranspose&, const DistributedMatrix<Real_wp, General, RowMajor>&, const Real_wp&, Matrix<Real_wp, Symmetric, RowSymPacked>&);
  
}


#endif
