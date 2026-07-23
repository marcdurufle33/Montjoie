#include "MontjoieFlag.hxx"

#include "Output/MontjoieOutputHeader.hxx"
#include "Output/MontjoieOutputInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Output/GridInterpolation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class GridInterpolation_Base<Dimension2>;
  SELDON_EXTERN template class GridInterpolation_Base<Dimension3>;

  SELDON_EXTERN template class GridInterpolationFull_Base<Dimension1>;
  SELDON_EXTERN template class GridInterpolationFull_Base<Dimension2>;
  SELDON_EXTERN template class GridInterpolationFull_Base<Dimension3>;

  SELDON_EXTERN template void WriteVtk(const VectReal_wp&, const GridInterpolationFull<Dimension2>&, const string&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const VectReal_wp&, const GridInterpolationFull<Dimension3>&, const string&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const VectComplex_wp&, const GridInterpolationFull<Dimension2>&, const string&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const VectComplex_wp&, const GridInterpolationFull<Dimension3>&, const string&, bool, bool);

  SELDON_EXTERN template void ReadVtk(VectReal_wp&, GridInterpolationFull<Dimension2>&, const string&, bool, bool);
  SELDON_EXTERN template void ReadVtk(VectReal_wp&, GridInterpolationFull<Dimension3>&, const string&, bool, bool);
  SELDON_EXTERN template void ReadVtk(VectComplex_wp&, GridInterpolationFull<Dimension2>&, const string&, bool, bool);
  SELDON_EXTERN template void ReadVtk(VectComplex_wp&, GridInterpolationFull<Dimension3>&, const string&, bool, bool);


  SELDON_EXTERN template void WriteMatlab(const VectReal_wp&, GridInterpolationFull<Dimension1>&, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectComplex_wp&, GridInterpolationFull<Dimension1>&, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectReal_wp&, GridInterpolationFull<Dimension2>&, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectReal_wp&, GridInterpolationFull<Dimension3>&, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectComplex_wp&, GridInterpolationFull<Dimension2>&, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectComplex_wp&, GridInterpolationFull<Dimension3>&, const string&, int, bool);

  SELDON_EXTERN template void WriteMatlab(const VectReal_wp&, int, const Vector<GridInterpolationFull<Dimension1> >&, int, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectComplex_wp&, int, const Vector<GridInterpolationFull<Dimension1> >&, int, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectReal_wp&, int, const Vector<GridInterpolationFull<Dimension2> >&, int, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectReal_wp&, int, const Vector<GridInterpolationFull<Dimension3> >&, int, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectComplex_wp&, int, const Vector<GridInterpolationFull<Dimension2> >&, int, const string&, int, bool);
  SELDON_EXTERN template void WriteMatlab(const VectComplex_wp&, int, const Vector<GridInterpolationFull<Dimension3> >&, int, const string&, int, bool);

  SELDON_EXTERN template void ReadMatlab(VectReal_wp&, GridInterpolationFull<Dimension2>&, const string&, bool);
  SELDON_EXTERN template void ReadMatlab(VectReal_wp&, GridInterpolationFull<Dimension3>&, const string&, bool);
  SELDON_EXTERN template void ReadMatlab(VectComplex_wp&, GridInterpolationFull<Dimension2>&, const string&, bool);
  SELDON_EXTERN template void ReadMatlab(VectComplex_wp&, GridInterpolationFull<Dimension3>&, const string&, bool);
  
}
