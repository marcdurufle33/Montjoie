#include "MontjoieFlag.hxx"

#include "Output/MontjoieOutputHeader.hxx"
#include "Output/MontjoieOutputInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Output/ParameterOutput.cxx"
#include "Output/ParameterOutputReprise.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template void ParameterOutputReprise::Read(Vector<TinyVector<Real_wp, 2> > &, bool);
  SELDON_EXTERN template void ParameterOutputReprise::Write(Vector<TinyVector<Real_wp, 2> > const &, bool);

  SELDON_EXTERN template void ParameterOutputReprise::Read(Vector<TinyVector<Real_wp, 3> > &, bool);
  SELDON_EXTERN template void ParameterOutputReprise::Write(Vector<TinyVector<Real_wp, 3> > const&, bool);

}
