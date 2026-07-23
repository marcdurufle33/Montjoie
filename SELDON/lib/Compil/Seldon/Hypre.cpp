#include "SeldonFlag.hxx"

#include "SeldonSolverHeader.hxx"
#include "SeldonSolverInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/solver/preconditioner/Hypre.cxx"
#endif

namespace Seldon
{

  SELDON_EXTERN template class HyprePreconditioner<HYPRE_Complex>;
  
  SELDON_EXTERN template void HyprePreconditioner<HYPRE_Complex>::ConstructPreconditioner(DistributedMatrix<HYPRE_Complex, General, ArrayRowSparse>&, bool);
  SELDON_EXTERN template void HyprePreconditioner<HYPRE_Complex>::ConstructPreconditioner(DistributedMatrix<HYPRE_Complex, Symmetric, ArrayRowSymSparse>&, bool);

}
