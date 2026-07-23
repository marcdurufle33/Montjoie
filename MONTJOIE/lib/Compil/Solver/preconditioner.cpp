#include "MontjoieFlag.hxx"

#include "Harmonic/MontjoieHarmonicHeader.hxx"
#include "Harmonic/MontjoieHarmonicInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Solver/Preconditioner.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class JacobiPreconditioning<Real_wp>;
  SELDON_EXTERN template class JacobiPreconditioning<Complex_wp>;

  SELDON_EXTERN template class SubdomainPreconditioning_Base<Real_wp>;
  SELDON_EXTERN template class SubdomainPreconditioning_Base<Complex_wp>;

  SELDON_EXTERN template class SubdomainPreconditioning_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class SubdomainPreconditioning_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class SubdomainPreconditioning_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class SubdomainPreconditioning_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template class LowOrderPreconditioning_Base<Real_wp>;
  SELDON_EXTERN template class LowOrderPreconditioning_Base<Complex_wp>;

  SELDON_EXTERN template class LowOrderPreconditioning_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class LowOrderPreconditioning_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class LowOrderPreconditioning_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class LowOrderPreconditioning_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template class MultigridPreconditioning_Base<Real_wp>;
  SELDON_EXTERN template class MultigridPreconditioning_Base<Complex_wp>;

  SELDON_EXTERN template class MultigridPreconditioning_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class MultigridPreconditioning_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class MultigridPreconditioning_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class MultigridPreconditioning_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template class All_Preconditioner_Base<Real_wp>;
  SELDON_EXTERN template class All_Preconditioner_Base<Complex_wp>;

  SELDON_EXTERN template void All_Preconditioner_Base<Real_wp>::GetMaximumSize(Real_wp, Dimension2, int&, int&, int&, int&);
  SELDON_EXTERN template void All_Preconditioner_Base<Real_wp>::GetMaximumSize(Real_wp, Dimension3, int&, int&, int&, int&);
  SELDON_EXTERN template void All_Preconditioner_Base<Complex_wp>::GetMaximumSize(Complex_wp, Dimension2, int&, int&, int&, int&);
  SELDON_EXTERN template void All_Preconditioner_Base<Complex_wp>::GetMaximumSize(Complex_wp, Dimension3, int&, int&, int&, int&);

  SELDON_EXTERN template void All_Preconditioner_Base<Real_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void All_Preconditioner_Base<Real_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void All_Preconditioner_Base<Complex_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);

  SELDON_EXTERN template class LocalPreconditioning_Base<Real_wp>;
  SELDON_EXTERN template class LocalPreconditioning_Base<Complex_wp>;

  SELDON_EXTERN template class LocalPreconditioning_Dim<Real_wp, Dimension2>;
  SELDON_EXTERN template class LocalPreconditioning_Dim<Complex_wp, Dimension2>;
  SELDON_EXTERN template class LocalPreconditioning_Dim<Real_wp, Dimension3>;
  SELDON_EXTERN template class LocalPreconditioning_Dim<Complex_wp, Dimension3>;

  SELDON_EXTERN template void LocalPreconditioning_Base<Real_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, const Vector<Real_wp>&, Vector<Real_wp>&);
  SELDON_EXTERN template void LocalPreconditioning_Base<Real_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Real_wp>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  SELDON_EXTERN template void LocalPreconditioning_Base<Complex_wp>::SolveGen(const SeldonTranspose&, const VirtualMatrix<Complex_wp>&, const Vector<Complex_wp>&, Vector<Complex_wp>&);
  
}
