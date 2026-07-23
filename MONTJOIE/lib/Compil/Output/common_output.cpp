#include "MontjoieFlag.hxx"

#include "Quadrature/MontjoieQuadratureHeader.hxx"
#include "Quadrature/MontjoieQuadratureInline.hxx"
#include "Output/CommonInputOutput.hxx"
#include "Output/SplineInterpolation.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Output/CommonInputOutput.cxx"
#include "Output/SplineInterpolation.cxx"
#endif

namespace Montjoie
{
  SELDON_EXTERN template class SplineInterpolation<Real_wp>;
  SELDON_EXTERN template class SplineInterpolation<Complex_wp>;
  
  SELDON_EXTERN template class BSplineInterpolation<Real_wp>;
  SELDON_EXTERN template class BSplineInterpolation<Complex_wp>;

  SELDON_EXTERN template void ReadMedit(Vector<Complex_wp>&, const string&);
  SELDON_EXTERN template void ReadMedit(Vector<Real_wp>&, const string&);
  
  SELDON_EXTERN template void WriteMedit(const Vector<Complex_wp>&, const string&, const Dimension2&, bool);
  SELDON_EXTERN template void WriteMedit(const Vector<Complex_wp>&, const string&, const Dimension3&, bool);
  SELDON_EXTERN template void WriteMedit(const Vector<Real_wp>&, const string&, const Dimension2&, bool);
  SELDON_EXTERN template void WriteMedit(const Vector<Real_wp>&, const string&, const Dimension3&, bool);

  SELDON_EXTERN template void ReadMeshData(Vector<Vector<Complex_wp> >&, const string&);
  SELDON_EXTERN template void ReadMeshData(Vector<Vector<Real_wp> >&, const string&);

  SELDON_EXTERN template void ReadMeshData(Vector<Vector<complex<float> > >&, const string&);
  SELDON_EXTERN template void ReadMeshData(Vector<Vector<float> >&, const string&);

  SELDON_EXTERN template void WriteMeshData(const Vector<Vector<Complex_wp> >&, const string&);
  SELDON_EXTERN template void WriteMeshData(const Vector<Vector<Real_wp> >&, const string&);

  SELDON_EXTERN template void WriteMeshData(const Vector<Vector<complex<float> > >&, const string&);
  SELDON_EXTERN template void WriteMeshData(const Vector<Vector<float> >&, const string&);
  
  SELDON_EXTERN template void WriteVtk(const Vector<Complex_wp>&, const string&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<Real_wp>&, const string&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<R3_Complex_wp>&, const string&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<R3>&, const string&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<TinyMatrix<Complex_wp, Symmetric, 3, 3> >&, const string&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> >&, const string&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<Vector<Complex_wp> >&, const IVect&, const Vector<IVect>&, const Vector<string>&, ostream&, bool, bool);
  SELDON_EXTERN template void WriteVtk(const Vector<Vector<Real_wp> >&, const IVect&, const Vector<IVect>&, const Vector<string>&, ostream&, bool, bool);
  
  SELDON_EXTERN template void WriteComplexVectorVtk(const Vector<Complex_wp>&, ostream&, bool, bool, bool);
  
  SELDON_EXTERN template void ReadVtk(Vector<Complex_wp>&, string&, istream&, bool ascii);
  SELDON_EXTERN template void ReadVtk(Vector<Real_wp>&, string&, istream&, bool ascii);
  SELDON_EXTERN template void ReadVtk(Vector<R3_Complex_wp>&, string&, istream&, bool ascii);
  SELDON_EXTERN template void ReadVtk(Vector<R3>&, string&, istream&, bool ascii);
  SELDON_EXTERN template void ReadVtk(Vector<TinyMatrix<Real_wp, Symmetric, 3, 3> >&, string&, istream&, bool ascii);
  SELDON_EXTERN template void WriteBinaryDoubleOrFloat(const Vector<float>&, ostream&, bool, bool, bool);
  SELDON_EXTERN template void WriteBinaryDoubleOrFloat(const Vector<double>&, ostream&, bool, bool, bool);
  SELDON_EXTERN template void WriteBinaryDoubleOrFloat(const Vector<complex<float> >&, ostream&, bool, bool, bool);
  SELDON_EXTERN template void WriteBinaryDoubleOrFloat(const Vector<complex<double> >&, ostream&, bool, bool, bool);

#ifdef MONTJOIE_WITH_MULTIPLE
  SELDON_EXTERN template void WriteBinaryDoubleOrFloat(const Vector<Real_wp>&, ostream&, bool, bool, bool);
  SELDON_EXTERN template void WriteBinaryDoubleOrFloat(const Vector<Complex_wp>&, ostream&, bool, bool, bool);
  SELDON_EXTERN template void ReadBinaryDoubleOrFloat(Vector<Real_wp>&, istream&, bool, bool, bool);
  SELDON_EXTERN template void ReadBinaryDoubleOrFloat(Vector<Complex_wp>&, istream&, bool, bool, bool);
#endif

  SELDON_EXTERN template void ReadBinaryDoubleOrFloat(Vector<float>&, istream&, bool, bool, bool);
  SELDON_EXTERN template void ReadBinaryDoubleOrFloat(Vector<double>&, istream&, bool, bool, bool);
  SELDON_EXTERN template void ReadBinaryDoubleOrFloat(Vector<complex<float> >&, istream&, bool, bool, bool);
  SELDON_EXTERN template void ReadBinaryDoubleOrFloat(Vector<complex<double> >&, istream&, bool, bool, bool);

  SELDON_EXTERN template class WriteOnTheGoWithBuffer<Real_wp>;
  SELDON_EXTERN template class WriteOnTheGoWithBuffer<Complex_wp>;

  SELDON_EXTERN template class WriteOnTheGoWithTinyBuffer<Real_wp, 1>;
  SELDON_EXTERN template class WriteOnTheGoWithTinyBuffer<Real_wp, 2>;
  SELDON_EXTERN template class WriteOnTheGoWithTinyBuffer<Real_wp, 3>;
  SELDON_EXTERN template class WriteOnTheGoWithTinyBuffer<Real_wp, 4>;
  SELDON_EXTERN template class WriteOnTheGoWithTinyBuffer<Real_wp, 5>;
  SELDON_EXTERN template class WriteOnTheGoWithTinyBuffer<Real_wp, 6>;
  
}
