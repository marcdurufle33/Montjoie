#ifndef MONTJOIE_FILE_CROSS_CORRELATION_FUNCTION_HXX

namespace Montjoie
{

  void ComputeFirstOrderCrossCorrelation(const VectComplex_wp &, const VectComplex_wp &,
                                         VectComplex_wp &);
  void ComputeEnvelopFourier(const VectComplex_wp &, VectComplex_wp &);
  int ReturnArgmin(const VectReal_wp &);
  Real_wp ReturnArgmin(const VectReal_wp &, const VectReal_wp &);
}

#define MONTJOIE_FILE_CROSS_CORRELATION_FUNCTION_HXX
#endif
