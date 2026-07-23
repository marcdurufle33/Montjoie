#ifndef MONTJOIE_FILE_CROSS_CORRELATION_FUNCTION_CXX

#include "CrossCorrelation/CrossCorrelationFunction.hxx"

namespace Montjoie
{
  void ComputeFirstOrderCrossCorrelation(const VectComplex_wp & Entry1,
                                       const VectComplex_wp & Entry2,
                                       VectComplex_wp & Result)
  {
    const int Size = Entry1.GetM();
    Result.Resize(Size);
#pragma omp parallel for
    for(int Index = 0; Index < Size; Index++)
      {
        if ((Index & 1) == 0)
          {
            Result(Index) = Entry1(Index) * Entry2(Index);
          }
        else
          {
            Result(Index) = - Entry1(Index) * Entry2(Index);
          }
        Result(Index) *= 2 / Real_wp(Size);
      }
  }

  int ReturnArgmin(const VectReal_wp & Entry)
  {
    int ReturnValue = 0;
    const int Size = Entry.GetM();
    for(int Index = 0; Index < Size; Index++)
      {
        if ( abs(Entry(Index)) > abs(Entry(ReturnValue)) )
          {
            ReturnValue = Index;
          }
      }
    return ReturnValue;
  }

  Real_wp ReturnArgmin(const VectReal_wp & Time, const VectReal_wp & Entry)
  {
    const int ReturnValue = ReturnArgmin(Entry);
    cout << ReturnValue << " " << Time(ReturnValue) << endl;
    return Time(ReturnValue);
  }

  void ComputeEnvelopFourier(const VectComplex_wp & EntryFourier,
                             VectReal_wp & Envelop)
  {
    FftInterface<Complex_wp> Fourier;
    VectComplex_wp Hilbert;
    const int Size = EntryFourier.GetM();
    Fourier.Init(Size);
    Hilbert.Resize(Size);
    Hilbert.Copy(EntryFourier);
    Envelop.Resize(Size);
#pragma omp parallel for
    for(int Index=0; Index<Size; Index++)
      {
        if ( ( (Index == Size/2) && ( (Size & 1) == 0) ) || (Index == 0) )
          {
            // We keep a ratio *1
          }
        else if (Index < Size/2)
          {
            Hilbert(Index)*=2;
          }
        else
          {
            Hilbert(Index)=0;
          }
      }
    Fourier.ApplyInverse(Hilbert);
#pragma omp parallel for
    for(int Index=0; Index<Size; Index++)
      {
        Envelop(Index) = abs(Hilbert(Index));
      }
    Hilbert.Clear();
  }

}

#define MONTJOIE_FILE_CROSS_CORRELATION_FUNCTION_CXX
#endif
