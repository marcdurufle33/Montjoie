#ifndef MONTJOIE_FILE_POINTS_REFERENCE_INLINE_CXX

namespace Montjoie
{  
  
  /************************
   * SetPoints<Dimension> *
   ************************/
  
  
  //! constructor
  template<class Dimension>
  inline SetPoints<Dimension>::SetPoints()
  {
  }
  
  
  //! allocate the array containing dof points
  template<class Dimension>
  inline void SetPoints<Dimension>::ReallocatePointsDof(int N)
  {
    PointsDof.Reallocate(N);
  }
  
  
  //! allocate the array containing quadrature points
  template<class Dimension>
  inline void SetPoints<Dimension>::ReallocatePointsQuadrature(int N)
  {
    PointsQuadrature.Reallocate(N);
  }


  //! allocate the array containing nodal points
  template<class Dimension>
  inline void SetPoints<Dimension>::ReallocatePointsNodal(int N)
  {
    PointsNodal.Reallocate(N);
  }


  //! allocate the array containing quadrature points on the boundary
  template<class Dimension>
  inline void SetPoints<Dimension>::ReallocatePointsQuadratureBoundary(int N)
  {
    PointsBoundary.Reallocate(N);
  }

  

  //! allocate the array containing dof points on the boundary
  template<class Dimension>
  inline void SetPoints<Dimension>::ReallocatePointsDofBoundary(int N)
  {
    PointsDofBoundary.Reallocate(N);
  }

  
  //! modification of a dof point
  template<class Dimension>
  inline void SetPoints<Dimension>::SetPointDof(int i, const R_N& mat)
  {
    PointsDof(i) = mat;
  }
  
  
  //! modification of a quadrature point
  template<class Dimension>
  inline void SetPoints<Dimension>::SetPointQuadrature(int i, const R_N& mat)
  {
    PointsQuadrature(i) = mat;
  }
  
  
  //! modification of a nodal point
  template<class Dimension>
  inline void SetPoints<Dimension>::SetPointNodal(int i, const R_N& mat)
  {
    PointsNodal(i) = mat;
  }
  
  
  //! modification of a quadrature point on the boundary
  template<class Dimension>
  inline void SetPoints<Dimension>::SetPointQuadratureBoundary(int i, const R_N& mat)
  {
    PointsBoundary(i) = mat;
  }
  

  //! modification of a dof point on the boundary
  template<class Dimension>
  inline void SetPoints<Dimension>::SetPointDofBoundary(int i, const R_N& mat)
  {
    PointsDofBoundary(i) = mat;
  }


  //! returns the number of nodal points stored
  template<class Dimension>
  inline int SetPoints<Dimension>::GetNbPointsNodal() const
  {
    return PointsNodal.GetM();
  }
  
  
  //! returns the number of quadrature points
  template<class Dimension>
  inline int SetPoints<Dimension>::GetNbPointsQuadrature() const
  {
    return PointsQuadrature.GetM();
  }
  
  
  //! returns the number of dof points
  template<class Dimension>
  inline int SetPoints<Dimension>::GetNbPointsDof() const
  {
    return PointsDof.GetM();
  }
  
  
  //! returns the number of quadrature points on the boundary
  template<class Dimension>
  inline int SetPoints<Dimension>::GetNbPointsQuadratureBoundary() const
  {
    return PointsBoundary.GetM();
  }


  //! returns the number of dof points on the boundary
  template<class Dimension>
  inline int SetPoints<Dimension>::GetNbPointsDofBoundary() const
  {
    return PointsBoundary.GetM();
  }
  
  
  //! returns a dof point
  template<class Dimension>
  inline const typename Dimension::R_N& SetPoints<Dimension>::GetPointDof(int i) const
  {
    return PointsDof(i);
  }
  
  
  //! returns a quadrature point
  template<class Dimension>
  inline const typename Dimension::R_N& SetPoints<Dimension>::GetPointQuadrature(int i) const
  {
    return PointsQuadrature(i);
  }
  
  
  //! returns a nodal point
  template<class Dimension>
  inline const typename Dimension::R_N& SetPoints<Dimension>::GetPointNodal(int i) const
  {
    return PointsNodal(i);
  }
  
  
  //! returns nodal points
  template<class Dimension>
  inline const Vector<typename Dimension::R_N>& SetPoints<Dimension>::GetPointNodal() const
  {
    return PointsNodal;
  }
  
  
  //! returns a quadrature point on the boundary
  template<class Dimension>
  inline const typename Dimension::R_N& SetPoints<Dimension>
  ::GetPointQuadratureBoundary(int i) const
  {
    return PointsBoundary(i);
  }


  //! returns a dof point on the boundary
  template<class Dimension>
  inline const typename Dimension::R_N& SetPoints<Dimension>
  ::GetPointDofBoundary(int i) const
  {
    return PointsDofBoundary(i);
  }
  
  
  //! nodal points -> quadrature points
  template<class Dimension>
  inline void SetPoints<Dimension>::CopyNodalToQuadrature()
  {
    PointsQuadrature = PointsNodal;
  }
  
  
  //! nodal points -> dof points
  template<class Dimension>
  inline void SetPoints<Dimension>::CopyNodalToDof()
  {
    PointsDof = PointsNodal;
  }
  
  
  //! quadrature points -> dof points
  template<class Dimension>
  inline void SetPoints<Dimension>::CopyQuadratureToDof()
  {
    PointsDof = PointsQuadrature;
  }
    
  
  /**************************
   * SetMatrices<Dimension> * 
   **************************/
  
  
  //! constructor with given problem
  template<class Dimension>
  inline SetMatrices<Dimension>::SetMatrices()
  {
    points_elem = NULL;
  }
  

  //! allocate the array containing jacobian matrices on dof points  
  template<class Dimension>
  inline void SetMatrices<Dimension>::ReallocatePointsDof(int N)
  {
    MatricesDof.Reallocate(N);
  }
  
  
  //! allocate the array containing jacobian matrices on quadrature points  
  template<class Dimension>
  inline void SetMatrices<Dimension>::ReallocatePointsQuadrature(int N)
  {
    MatricesQuadrature.Reallocate(N);
    Weights_ND.Reallocate(N); Weights_ND.Fill(0);
  }


  //! allocate the array containing jacobian matrices on nodal points  
  template<class Dimension>
  inline void SetMatrices<Dimension>::ReallocatePointsNodal(int N)
  {
    MatricesNodal.Reallocate(N);
  }

  
  //! allocate the array containing jacobian matrices on boundary dof points 
  template<class Dimension>
  inline void SetMatrices<Dimension>::ReallocatePointsDofBoundary(int N)
  {
    MatricesDofBoundary.Reallocate(N);
  }


  //! returns the number of quadrature points on the boundary
  template<class Dimension>
  inline int SetMatrices<Dimension>::GetNbPointsQuadratureBoundary() const
  {
    return MatricesBoundary.GetM();
  }


  //! returns the number of dof points on the boundary
  template<class Dimension>
  inline int SetMatrices<Dimension>::GetNbPointsDofBoundary() const
  {
    return MatricesDofBoundary.GetM();
  }
  
  
  //! modification of DFi(x) where x is a dof point
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetPointDof(int i, const MatrixN_N& mat)
  {
    MatricesDof(i) = mat;
  }
  
  
  //! modification of DFi(x) where x is a quadrature point
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetPointQuadrature(int i, const MatrixN_N& mat)
  {
    MatricesQuadrature(i) = mat;
  }
  
  
  //! modification of DFi(x) where x is a nodal point
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetPointNodal(int i, const MatrixN_N& mat)
  {
    MatricesNodal(i) = mat;
  }
  
  
  //! modification of DFi(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetPointQuadratureBoundary(int i, const MatrixN_N& mat)
  {
    MatricesBoundary(i) = mat;
  }


  //! modification of DFi(x) where x is a dof point on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetPointDofBoundary(int i, const MatrixN_N& mat)
  {
    MatricesDofBoundary(i) = mat;
  }


  //! sets jacobian matrices on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetPointNodalBoundary(const Vector<MatrixN_N>& mat)
  {
    MatricesNodalBoundary = mat;
  }

  
  //! sets inverse of jacobian matrices on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetInversePointNodalBoundary(const Vector<MatrixN_N>& mat)
  {
    invMatricesNodalBoundary = mat;
  }
  
  
  //! modification of normale n(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetNormaleQuadratureBoundary(int i, const R_N& mat)
  {
    NormaleQuadrature(i) = mat;
  }
  
  
  //! modification of surfacic element ds(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetDsQuadratureBoundary(int i, const Real_wp& mat)
  {
    DsQuadrature(i) = mat;
  }
  
  
  //! modification of curvature K1(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetK1QuadratureBoundary(int i, const Real_wp& mat)
  {
    K1_curve(i) = mat;
  }
  
  
  //! modification of curvature K2(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline void SetMatrices<Dimension>::SetK2QuadratureBoundary(int i, const Real_wp& mat)
  {
    K2_curve(i) = mat;
  }
  
  
  //! constant matrix for all nodal points
  template<class Dimension>
  inline void SetMatrices<Dimension>::FillNodal(const MatrixN_N& mat)
  {
    MatricesNodal.Fill(mat);
  }
  
  
  //! constant matrix for all quadrature points
  template<class Dimension>
  inline void SetMatrices<Dimension>::FillQuadrature(const MatrixN_N& mat)
  {
    MatricesQuadrature.Fill(mat);
  }
  
  
  //! constant matrix for all dof points
  template<class Dimension>
  inline void SetMatrices<Dimension>::FillDof(const MatrixN_N& mat)
  {
    MatricesDof.Fill(mat);
  }
  
  
  //! returns DFi(x) where x is a quadrature point
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>
  ::GetPointQuadrature(int i) const
  {
    return MatricesQuadrature(i);
  }
  

  //! returns DFi(x) where x is a dof point
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>::GetPointDof(int i) const
  {
    return MatricesDof(i);
  }
  
    
  //! returns DFi(x) where x is a nodal point
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>::GetPointNodal(int i) const
  {
    return MatricesNodal(i);
  }
  
  
  //! returns DFi(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>
  ::GetPointQuadratureBoundary(int i) const
  {
    return MatricesBoundary(i);
  }


  //! returns DFi(x) where x is a dof point on the boundary
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>
  ::GetPointDofBoundary(int i) const
  {
    return MatricesDofBoundary(i);
  }
  
  
  //! returns DFi(x) where x is a nodal point on the boundary
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>
  ::GetPointNodalBoundary(int i) const
  {
    return MatricesNodalBoundary(i);
  }
  

  //! returns inverse of DFi(x) where x is a nodal point on the boundary  
  template<class Dimension>
  inline const typename Dimension::MatrixN_N& SetMatrices<Dimension>
  ::GetInversePointNodalBoundary(int i) const
  {
    return invMatricesNodalBoundary(i);
  }
  
  
  //! returns normale n(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline const typename Dimension::R_N& SetMatrices<Dimension>
  ::GetNormaleQuadratureBoundary(int i) const
  {
    return NormaleQuadrature(i);
  }
  
  
  //! returns normale n(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline const Real_wp& SetMatrices<Dimension>::GetDsQuadratureBoundary(int i) const
  {
    return DsQuadrature(i);
  }
  
  
  //! returns curvature K1(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline const Real_wp& SetMatrices<Dimension>::GetK1QuadratureBoundary(int i) const
  {
    return K1_curve(i);
  }
  
  
  //! returns curvature K2(x) where x is a quadrature point on the boundary
  template<class Dimension>
  inline const Real_wp& SetMatrices<Dimension>::GetK2QuadratureBoundary(int i) const
  {
    return K2_curve(i);
  }
  
  
  //! copy DFi
  template<class Dimension>
  inline void SetMatrices<Dimension>::CopyNodalToQuadrature()
  {
    MatricesQuadrature = MatricesNodal;
  }
  
  
  //! copy DFi
  template<class Dimension>
  inline void SetMatrices<Dimension>::CopyNodalToDof()
  {
    MatricesDof = MatricesNodal;
  }
  
  
  //! copy DFi
  template<class Dimension>
  inline void SetMatrices<Dimension>::CopyQuadratureToDof()
  {
    MatricesDof = MatricesQuadrature;
  }

  template<class Dimension>
  inline void SetMatrices<Dimension>::InitSetPoints(const SetPoints<Dimension>& pts)
  {
    points_elem = const_cast<SetPoints<Dimension>* >(&pts);
  }

  template<class Dimension>
  inline const SetPoints<Dimension>& SetMatrices<Dimension>::GetSetPoints() const
  {
    if (points_elem == NULL)
      abort();

    return *points_elem;
  }
  
    
  /********************
   * FjInverseProblem *
   ********************/
  

  //! returns vertices of the element
  template<class Dimension>
  inline typename Dimension::VectR_N& FjInverseProblem<Dimension>::GetVertices()
  {
    return Vertices;
  }
  
  
  //! returns nodal points for the element
  template<class Dimension>
  inline SetPoints<Dimension>& FjInverseProblem<Dimension>::GetSetPoints()
  {
    return PointsElem;
  }


  //! initial guess -> point at the center of the element
  template<class Dimension>
  inline void FjInverseProblem<Dimension>::FindInitGuess(R_N& x)
  {
    x = Fb.GetCenterReferenceElement();
  }

  
  template<class Dimension>
  inline void FjInverseProblem<Dimension>::ComputeScheme(const R_N& x, R_N& fvec)
  {
    EvaluateFunction(x, fvec);
  }


  template<class Dimension>
  inline void FjInverseProblem<Dimension>::ComputeAndFactoriseDiff(const R_N& x, const R_N&)
  {    
    EvaluateJacobian(x, invDF);
    GetInverse(invDF);
  }


  template<class Dimension>
  inline void FjInverseProblem<Dimension>::SolveDifferential(const R_N& x, R_N& y)
  {
    Mlt(invDF, x, y);
  }


  template<class Dimension>
  inline Real_wp FjInverseProblem<Dimension>::GetNorm2Vector(const R_N& x)
  {
    return Norm2(x);
  }


  template<class Dimension>
  inline void FjInverseProblem<Dimension>::EvaluateF(const R_N& x, R_N& feval)
  {
    EvaluateFunction(x, feval);
  }
  

  template<class Dimension>
  inline void FjInverseProblem<Dimension>::EvaluateJacobian(const R_N& x, R_N& feval, MatrixN_N& fjac)
  {
    EvaluateFunction(x, feval);
    EvaluateJacobian(x, fjac);
  }
  
} // end namespace

#define MONTJOIE_FILE_POINTS_REFERENCE_INLINE_CXX
#endif
