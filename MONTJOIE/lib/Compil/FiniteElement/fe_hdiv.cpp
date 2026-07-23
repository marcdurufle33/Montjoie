#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/FiniteElementHdiv.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class FiniteElementHdiv<Dimension2>;
  SELDON_EXTERN template class FiniteElementHdiv<Dimension3>;

  SELDON_EXTERN template void FiniteElementHdiv<Dimension2>::GetCurlFromGradient(const Vector<VectReal_wp>&, Vector<VectReal_wp>&);
  SELDON_EXTERN template void FiniteElementHdiv<Dimension2>::GetCurlFromGradient(const Vector<VectComplex_wp>&, Vector<VectComplex_wp>&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeNodalValues(const SetMatrices<Dimension2>&,
                                                         const Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                                                         const Mesh<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeNodalValues(const SetMatrices<Dimension2>&,
                                                         const Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                                                         const Mesh<Dimension2>&, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeQuadratureValues(const SetMatrices<Dimension2>&,
                                                              const Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                              Vector<VectReal_wp >&, bool, bool,
                                                              const Mesh<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeQuadratureValues(const SetMatrices<Dimension2>&,
                                                              const Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                              Vector<VectComplex_wp >&, bool, bool,
                                                              const Mesh<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeProjectionDof(const SetMatrices<Dimension2>&,
                                                           Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                           const MeshNumbering<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeProjectionDof(const SetMatrices<Dimension2>&,
                                                           Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                           const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeProjectionSurfaceDof(const SetMatrices<Dimension2>& ,
                                                                  Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                  const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeProjectionSurfaceDof(const SetMatrices<Dimension2>& ,
                                                                  Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                  const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegral(const SetMatrices<Dimension2>&,
                                                      Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                      const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegral(const SetMatrices<Dimension2>&,
                                                      Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                      const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegralGradient(const SetMatrices<Dimension2>&,
                                                              Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                              const MeshNumbering<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegralGradient(const SetMatrices<Dimension2>&,
                                                              Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                              const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegralSurface(const SetMatrices<Dimension2>&,
                                                             Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                             const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegralSurface(const SetMatrices<Dimension2>&,
                                                             Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                             const MeshNumbering<Dimension2>&, int, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeGaussIntegralSurface(const Vector<Matrix2_2>&, const VectReal_wp&, const VectReal_wp&,
                                                                  Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                  const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeGaussIntegralSurface(const Vector<Matrix2_2>&, const VectReal_wp&, const VectReal_wp&,
                                                                  Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                  const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegralSurfaceGradient(const SetMatrices<Dimension2>&,
                                                                     Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                     const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeIntegralSurfaceGradient(const SetMatrices<Dimension2>&,
                                                                     Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                     const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeValueNodalBoundary(const VectReal_wp&, VectReal_wp&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeValueNodalBoundary(const VectComplex_wp&, VectComplex_wp&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeValueBoundary(const SetMatrices<Dimension2>&,
                                                           const Vector<VectReal_wp >&, Vector<Vector<TinyVector<Real_wp, 2> > >& res,
                                                           const Mesh<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeValueBoundary(const SetMatrices<Dimension2>&,
                                                           const Vector<VectComplex_wp >&, Vector<Vector<TinyVector<Complex_wp, 2> > >& res,
                                                           const Mesh<Dimension2>&, int, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeGradientBoundary(const SetMatrices<Dimension2>&,
                                                              const Vector<VectReal_wp >&, Vector<Vector<TinyVector<Real_wp, 1> > >&,
                                                              const Mesh<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ComputeGradientBoundary(const SetMatrices<Dimension2>&,
                                                              const Vector<VectComplex_wp >&, Vector<Vector<TinyVector<Complex_wp, 1> > >&,
                                                              const Mesh<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Real_wp, 2>&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Complex_wp, 2>&);
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Real_wp, 2>&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Complex_wp, 2>&);
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, VectReal_wp&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, VectComplex_wp&);


  // 3-D
  SELDON_EXTERN template void FiniteElementHdiv<Dimension3>::GetCurlFromGradient(const Vector<VectReal_wp>&, Vector<VectReal_wp>&);
  SELDON_EXTERN template void FiniteElementHdiv<Dimension3>::GetCurlFromGradient(const Vector<VectComplex_wp>&, Vector<VectComplex_wp>&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeNodalValues(const SetMatrices<Dimension3>&,
                                                         const Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                                                         const Mesh<Dimension3>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeNodalValues(const SetMatrices<Dimension3>&,
                                                         const Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                                                         const Mesh<Dimension3>&, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeQuadratureValues(const SetMatrices<Dimension3>&,
                                                              const Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                              Vector<VectReal_wp >&, bool, bool,
                                                              const Mesh<Dimension3>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeQuadratureValues(const SetMatrices<Dimension3>&,
                                                              const Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                              Vector<VectComplex_wp >&, bool, bool,
                                                              const Mesh<Dimension3>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeProjectionDof(const SetMatrices<Dimension3>&,
                                                           Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                           const MeshNumbering<Dimension3>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeProjectionDof(const SetMatrices<Dimension3>&,
                                                           Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                           const MeshNumbering<Dimension3>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeProjectionSurfaceDof(const SetMatrices<Dimension3>& ,
                                                                  Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                  const MeshNumbering<Dimension3>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeProjectionSurfaceDof(const SetMatrices<Dimension3>& ,
                                                                  Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                  const MeshNumbering<Dimension3>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegral(const SetMatrices<Dimension3>&,
                                                      Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                      const MeshNumbering<Dimension3>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegral(const SetMatrices<Dimension3>&,
                                                      Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                      const MeshNumbering<Dimension3>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegralGradient(const SetMatrices<Dimension3>&,
                                                              Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                              const MeshNumbering<Dimension3>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegralGradient(const SetMatrices<Dimension3>&,
                                                              Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                              const MeshNumbering<Dimension3>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegralSurface(const SetMatrices<Dimension3>&,
                                                             Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                             const MeshNumbering<Dimension3>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegralSurface(const SetMatrices<Dimension3>&,
                                                             Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                             const MeshNumbering<Dimension3>&, int, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeGaussIntegralSurface(const Vector<Matrix3_3>&, const VectReal_wp&, const VectReal_wp&,
                                                                  Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                  const MeshNumbering<Dimension3>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeGaussIntegralSurface(const Vector<Matrix3_3>&, const VectReal_wp&, const VectReal_wp&,
                                                                  Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                  const MeshNumbering<Dimension3>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegralSurfaceGradient(const SetMatrices<Dimension3>&,
                                                                     Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                     const MeshNumbering<Dimension3>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeIntegralSurfaceGradient(const SetMatrices<Dimension3>&,
                                                                     Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                     const MeshNumbering<Dimension3>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeValueNodalBoundary(const VectReal_wp&, VectReal_wp&, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeValueNodalBoundary(const VectComplex_wp&, VectComplex_wp&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeValueBoundary(const SetMatrices<Dimension3>&,
                                                           const Vector<VectReal_wp >&, Vector<Vector<TinyVector<Real_wp, 3> > >& res,
                                                           const Mesh<Dimension3>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeValueBoundary(const SetMatrices<Dimension3>&,
                                                           const Vector<VectComplex_wp >&, Vector<Vector<TinyVector<Complex_wp, 3> > >& res,
                                                           const Mesh<Dimension3>&, int, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeGradientBoundary(const SetMatrices<Dimension3>&,
                                                              const Vector<VectReal_wp >&, Vector<Vector<TinyVector<Real_wp, 1> > >&,
                                                              const Mesh<Dimension3>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ComputeGradientBoundary(const SetMatrices<Dimension3>&,
                                                              const Vector<VectComplex_wp >&, Vector<Vector<TinyVector<Complex_wp, 1> > >&,
                                                              const Mesh<Dimension3>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Real_wp, 3>&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Complex_wp, 3>&);
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Real_wp, 3>&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Complex_wp, 3>&);
  
  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, VectReal_wp&);

  SELDON_EXTERN template
  void FiniteElementHdiv<Dimension3>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, VectComplex_wp&);
  
}
