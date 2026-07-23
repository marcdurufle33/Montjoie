#include "MontjoieFlag.hxx"

#include "FiniteElement/MontjoieFiniteElementHeader.hxx"
#include "FiniteElement/MontjoieFiniteElementInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "FiniteElement/FiniteElementHcurl2D.cxx"
#endif

namespace Montjoie
{

  SELDON_EXTERN template class FiniteElementHcurl<Dimension2>;

  SELDON_EXTERN template void FiniteElementHcurl<Dimension2>::GetCurlFromGradient(const Vector<VectReal_wp>&, Vector<VectReal_wp>&);
  SELDON_EXTERN template void FiniteElementHcurl<Dimension2>::GetCurlFromGradient(const Vector<VectComplex_wp>&, Vector<VectComplex_wp>&);

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeNodalValues(const SetMatrices<Dimension2>&,
							  const Vector<VectReal_wp>&, Vector<VectReal_wp>&,
							  const Mesh<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeNodalValues(const SetMatrices<Dimension2>&,
                                                          const Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                                                          const Mesh<Dimension2>&, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeQuadratureValues(const SetMatrices<Dimension2>&,
                                                               const Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                               Vector<VectReal_wp >&, bool, bool,
                                                               const Mesh<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeQuadratureValues(const SetMatrices<Dimension2>&,
                                                               const Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                               Vector<VectComplex_wp >&, bool, bool,
                                                               const Mesh<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeProjectionDof(const SetMatrices<Dimension2>&,
                                                            Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                            const MeshNumbering<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeProjectionDof(const SetMatrices<Dimension2>&,
                                                            Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                            const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeProjectionSurfaceDof(const SetMatrices<Dimension2>& ,
                                                                   Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                   const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeProjectionSurfaceDof(const SetMatrices<Dimension2>& ,
                                                                   Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                   const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegral(const SetMatrices<Dimension2>&,
                                                       Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                       const MeshNumbering<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegral(const SetMatrices<Dimension2>&,
                                                       Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                       const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegralGradient(const SetMatrices<Dimension2>&,
                                                               Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                               const MeshNumbering<Dimension2>&, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegralGradient(const SetMatrices<Dimension2>&,
                                                               Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                               const MeshNumbering<Dimension2>&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegralSurface(const SetMatrices<Dimension2>&,
                                                              Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                              const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegralSurface(const SetMatrices<Dimension2>&,
                                                              Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                              const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeGaussIntegralSurface(const Vector<Matrix2_2>&, const VectReal_wp&, const VectReal_wp&,
                                                                   Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                   const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeGaussIntegralSurface(const Vector<Matrix2_2>&, const VectReal_wp&, const VectReal_wp&,
                                                                   Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                   const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegralSurfaceGradient(const SetMatrices<Dimension2>&,
                                                                      Vector<VectReal_wp >&, Vector<VectReal_wp >&,
                                                                      const MeshNumbering<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeIntegralSurfaceGradient(const SetMatrices<Dimension2>&,
                                                                      Vector<VectComplex_wp >&, Vector<VectComplex_wp >&,
                                                                      const MeshNumbering<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeValueNodalBoundary(const VectReal_wp&, VectReal_wp&, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeValueNodalBoundary(const VectComplex_wp&, VectComplex_wp&, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeValueBoundary(const SetMatrices<Dimension2>&,
                                                            const Vector<VectReal_wp >&, Vector<Vector<TinyVector<Real_wp, 2> > >& res,
                                                            const Mesh<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeValueBoundary(const SetMatrices<Dimension2>&,
                                                            const Vector<VectComplex_wp >&, Vector<Vector<TinyVector<Complex_wp, 2> > >& res,
                                                            const Mesh<Dimension2>&, int, int) const;
    
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeGradientBoundary(const SetMatrices<Dimension2>&,
                                                               const Vector<VectReal_wp >&, Vector<Vector<TinyVector<Real_wp, 1> > >&,
                                                               const Mesh<Dimension2>&, int, int) const;

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ComputeGradientBoundary(const SetMatrices<Dimension2>&,
                                                               const Vector<VectComplex_wp >&, Vector<Vector<TinyVector<Complex_wp, 1> > >&,
                                                               const Mesh<Dimension2>&, int, int) const;
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Real_wp, 2>&);

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ApplyRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Complex_wp, 2>&);
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Real_wp, 2>&);

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, TinyVector<Complex_wp, 2>&);
  
  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, VectReal_wp&);

  SELDON_EXTERN template
  void FiniteElementHcurl<Dimension2>::ApplyInverseRotationCyclic(int n, FftInterface<Complex_wp>&, VectComplex_wp&);

  
}
