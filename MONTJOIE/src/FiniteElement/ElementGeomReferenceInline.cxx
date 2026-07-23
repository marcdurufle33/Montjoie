#ifndef MONTJOIE_FILE_ELEMENT_GEOM_REFERENCE_INLIN_CXX

namespace Montjoie
{
  //! default constructor
  template<class Dimension>
  inline ElementGeomReference_Base<Dimension>::ElementGeomReference_Base()
  {
    order_geom = -1;
    
    // we assume that all points are different
    quadrature_equal_nodal = false;
    dof_equal_nodal = false;
    dof_equal_quadrature = false;
  }
  
  
  //! destructor
  template<class Dimension>
  inline ElementGeomReference_Base<Dimension>::~ElementGeomReference_Base()
  {
  }  
 
  
  //! returns order used for geometry
  template<class Dimension>
  inline int ElementGeomReference_Base<Dimension>::GetGeometryOrder() const
  {
    return order_geom;
  }
  
  
  //! returns normale of the side n
  template<class Dimension>
  inline const typename Dimension::R_N& ElementGeomReference_Base<Dimension>::NormaleLoc(int n) const
  {
    return normale(n);
  }
  

  //! 1-D nodal points
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReference_Base<Dimension>::PointsNodal1D() const
  {
    return points_nodal1d;
  }
  
  
  //! 1-D nodal point
  template<class Dimension>
  inline const Real_wp& ElementGeomReference_Base<Dimension>::PointsNodal1D(int j) const
  {
    return points_nodal1d(j);
  }


  //! 2-D or 3-D nodal points
  template<class Dimension>
  inline const Vector<typename Dimension::R_N>&
  ElementGeomReference_Base<Dimension>::PointsNodalND() const
  {
    return points_nodal_nd;
  }

  
  //! 2-D or 3-D nodal points
  template<class Dimension>
  inline const typename Dimension::R_N&
  ElementGeomReference_Base<Dimension>::PointsNodalND(int j) const
  {
    return points_nodal_nd(j);
  }


  //! 2-D or 3-D nodal weights
  template<class Dimension>
  inline const Vector<Real_wp>&
  ElementGeomReference_Base<Dimension>::WeightsNodalND() const
  {
    return weights_nodal_nd;
  }

  
  //! 2-D or 3-D nodal weights
  template<class Dimension>
  inline const Real_wp&
  ElementGeomReference_Base<Dimension>::WeightsNodalND(int j) const
  {
    return weights_nodal_nd(j);
  }

  
  //! 1-D dof points
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReference_Base<Dimension>::PointsDof1D() const
  {
    return points_dof1d;
  }
  
  
  //! 1-D dof point
  template<class Dimension>
  inline const Real_wp& ElementGeomReference_Base<Dimension>::PointsDof1D(int j) const
  {
    return points_dof1d(j);
  }
  
  
  //! 1-D quadrature points
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReference_Base<Dimension>::Points1D() const
  {
    return points1d;
  }
  
  
  //! 1-D quadrature point
  template<class Dimension>
  inline const Real_wp& ElementGeomReference_Base<Dimension>::Points1D(int j) const
  {
    return points1d(j);
  }
  
  
  //! 1-D quadrature weights
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReference_Base<Dimension>::Weights1D() const
  {
    return weights1d;
  }
  
  
  //! 1-D quadrature weight
  template<class Dimension>
  inline const Real_wp& ElementGeomReference_Base<Dimension>::Weights1D(int j) const
  {
    return weights1d(j);
  }
    
  
  //! returns the hybrid type of the element
  template<class Dimension>
  inline int ElementGeomReference_Base<Dimension>::GetHybridType() const
  {
    return hybrid_type_elt;
  }
  

  //! return the number of vertices for the reference element
  template<class Dimension>
  inline int ElementGeomReference_Base<Dimension>::GetNbVertices() const
  { 
    return nb_vertices_elt;
  }


  //! return the number of edges (or faces in 3-D) for the reference element
  template<class Dimension>
  inline int ElementGeomReference_Base<Dimension>::GetNbBoundaries() const
  { 
    return nb_boundaries_elt;
  }

  
  //! return the number of edges for the reference element
  template<class Dimension>
  inline int ElementGeomReference_Base<Dimension>::GetNbEdges() const
  { 
    return nb_edges_elt;
  }


  template<class Dimension>
  inline void ElementGeomReference_Base<Dimension>::SetAxisymmetricGeometry()
  {
  }
  
  
  template<class Dimension>
  inline void ElementGeomReference_Base<Dimension>::ComputeCoefficientTransformation()
  {
    abort(); 
  }
  

  /**************************************
   * Container for ElementGeomReference *
   **************************************/
  
  template<class Dimension>
  inline ElementGeomReferenceContainer_Base<Dimension>
  ::ElementGeomReferenceContainer_Base(ElementGeomReference<Dimension>& elt)
    : elt_geom(elt)
  {
  }
    
  template<class Dimension>
  inline bool ElementGeomReferenceContainer_Base<Dimension>::QuadratureEqualNodal() const
  {
    return elt_geom.quadrature_equal_nodal;
  }
  
  template<class Dimension>
  inline bool ElementGeomReferenceContainer_Base<Dimension>::DofEqualNodal() const
  {
    return elt_geom.dof_equal_nodal; 
  }
  
  template<class Dimension>
  inline bool ElementGeomReferenceContainer_Base<Dimension>::DofEqualQuadrature() const
  {
    return elt_geom.dof_equal_quadrature;
  }

  template<class Dimension>
  inline int ElementGeomReferenceContainer_Base<Dimension>::GetGeometryOrder() const
  {
    return elt_geom.GetGeometryOrder();
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>::SetGeometryOrder(int r)
  {
    elt_geom.order_geom = r;
  }

  template<class Dimension>
  inline ElementGeomReference<Dimension>& ElementGeomReferenceContainer_Base<Dimension>
  ::GetGeometricElement()
  {
    return elt_geom;
  }
  
  template<class Dimension>
  inline const ElementGeomReference<Dimension>& 
  ElementGeomReferenceContainer_Base<Dimension>::GetGeometricElement() const
  {
    return elt_geom;
  }
  
  template<class Dimension>
  inline const typename Dimension::R_N&
  ElementGeomReferenceContainer_Base<Dimension>::NormaleLoc(int n) const
  {
    return elt_geom.NormaleLoc(n);
  }
    
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReferenceContainer_Base<Dimension>::PointsNodal1D() const
  {
    return elt_geom.PointsNodal1D();
  }
  
  template<class Dimension>
  inline const Real_wp& ElementGeomReferenceContainer_Base<Dimension>::PointsNodal1D(int i) const
  {
    return elt_geom.PointsNodal1D(i);
  }
    
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReferenceContainer_Base<Dimension>::PointsDof1D() const
  {
    return elt_geom.PointsDof1D();
  }
  
  template<class Dimension>
  inline const Real_wp& ElementGeomReferenceContainer_Base<Dimension>::PointsDof1D(int i) const
  {
    return elt_geom.PointsDof1D(i);
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>::SetPointsDof1D(const VectReal_wp& pts)
  {
    elt_geom.points_dof1d = pts;
  }
  
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReferenceContainer_Base<Dimension>::Points1D() const
  {
    return elt_geom.Points1D();
  }
  
  template<class Dimension>
  inline const Real_wp& ElementGeomReferenceContainer_Base<Dimension>::Points1D(int i) const
  {
    return elt_geom.Points1D(i);
  }
  
  template<class Dimension>
  inline const VectReal_wp& ElementGeomReferenceContainer_Base<Dimension>::Weights1D() const
  {
    return elt_geom.Weights1D();
  }
  
  template<class Dimension>
  inline const Real_wp& ElementGeomReferenceContainer_Base<Dimension>::Weights1D(int i) const
  {
    return elt_geom.Weights1D(i);
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>::SetPoints1D(const VectReal_wp& pts)
  {
    elt_geom.points1d = pts;
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>::SetWeights1D(const VectReal_wp& w)
  {
    elt_geom.weights1d = w;
  }
    
  template<class Dimension>
  inline int ElementGeomReferenceContainer_Base<Dimension>::GetHybridType() const
  {
    return elt_geom.GetHybridType();
  }
  
  template<class Dimension>
  inline int ElementGeomReferenceContainer_Base<Dimension>::GetNbVertices() const
  {
    return elt_geom.GetNbVertices();
  }
  
  template<class Dimension>
  inline int ElementGeomReferenceContainer_Base<Dimension>::GetNbEdges() const
  {
    return elt_geom.GetNbEdges();
  }
  
  template<class Dimension>
  inline int ElementGeomReferenceContainer_Base<Dimension>::GetNbBoundaries() const
  {
    return elt_geom.GetNbBoundaries();
  }
  
  template<class Dimension>
  inline FiniteElementProjector* ElementGeomReferenceContainer_Base<Dimension>
  ::GetNewNodalInterpolation() const
  {
    return elt_geom.GetNewNodalInterpolation();
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::SetMesh(Mesh<Dimension>& mesh)
  {
    elt_geom.SetMesh(mesh);
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::GetNormale(const MatrixN_N& dfjm1, R_N& normale_fj, Real_wp& dsj, int num_edge) const
  {
    elt_geom.GetNormale(dfjm1, normale_fj, dsj, num_edge);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::Fj(const VectR_N& s, const SetPoints<Dimension>& PTReel,
       const R_N& point, R_N& res, const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.Fj(s, PTReel, point, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::FjLinear(const VectR_N& s, const R_N& point, R_N& res) const
  {
    elt_geom.FjLinear(s, point, res);
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::DFjLinear(const VectR_N& s, const R_N& point, MatrixN_N& res) const
  {
    elt_geom.DFjLinear(s, point, res);
  }
    
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::DFj(const VectR_N& s, const SetPoints<Dimension>& PTReel,
	const R_N& point, MatrixN_N& res,
	const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.DFj(s, PTReel, point, res, mesh, nquad);
  }
    
  template<class Dimension>
  inline Real_wp ElementGeomReferenceContainer_Base<Dimension>
  ::GetMinimalSize(const VectR_N& s) const
  {
    return elt_geom.GetMinimalSize(s);
  }
    
  template<class Dimension>
  inline bool ElementGeomReferenceContainer_Base<Dimension>
  ::OutsideReferenceElement(const VectR_N& s, const R_N& Xn, const Real_wp& epsilon) const
  {
    return elt_geom.OutsideReferenceElement(s, Xn, epsilon);
  }
  
  template<class Dimension>
  inline Real_wp ElementGeomReferenceContainer_Base<Dimension>
  ::GetDistanceToBoundary(const R_N& pointloc) const
  {
    return elt_geom.GetDistanceToBoundary(pointloc);
  }
  
  template<class Dimension>
  inline int ElementGeomReferenceContainer_Base<Dimension>
  ::ProjectPointOnBoundary(R_N& pointloc) const
  {
    return elt_geom.ProjectPointOnBoundary(pointloc);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeCoefJacobian(const VectR_N& s, VectReal_wp& CoefJacobian) const
  {
    elt_geom.ComputeCoefJacobian(s, CoefJacobian);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::FjElem(const VectR_N& s, SetPoints<Dimension>& res,
	   const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.FjElem(s, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::FjElemNodal(const VectR_N& s, SetPoints<Dimension>& res,
		const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.FjElemNodal(s, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::FjElemQuadrature(const VectR_N& s, SetPoints<Dimension>& res,
		     const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.FjElemQuadrature(s, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::FjElemDof(const VectR_N& s, SetPoints<Dimension>& res,
	      const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.FjElemDof(s, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::DFjElem(const VectR_N& s, const SetPoints<Dimension>& PTReel,
	    SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
	    int nquad) const
  {
    elt_geom.DFjElem(s, PTReel, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::DFjElemNodal(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		 SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
		 int nquad) const
  {
    elt_geom.DFjElemNodal(s, PTReel, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::DFjElemQuadrature(const VectR_N& s, const SetPoints<Dimension>& PTReel,
		      SetMatrices<Dimension>& res, const Mesh<Dimension>& mesh,
		      int nquad) const
  {
    elt_geom.DFjElemQuadrature(s, PTReel, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::DFjElemDof(const VectR_N& s, const SetPoints<Dimension>& PTReel,
	       SetMatrices<Dimension>& res,
	       const Mesh<Dimension>& mesh, int nquad) const
  {
    elt_geom.DFjElemDof(s, PTReel, res, mesh, nquad);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeValuesNodalPhi1D(const Real_wp& x, VectReal_wp& phi) const
  {
    elt_geom.ComputeValuesNodalPhi1D(x, phi);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeValuesPhiNodalRef(const R_N& pointloc, VectReal_wp& phi) const
  {
    elt_geom.ComputeValuesPhiNodalRef(pointloc, phi);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeGradientPhiNodalRef(const R_N& pointloc, VectR_N& grad_phi) const
  {
    elt_geom.ComputeGradientPhiNodalRef(pointloc, grad_phi);
  }

  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeNodalGradientRef(const Vector<Real_wp>& Un, VectR_N& Unode) const
  {
    return this->elt_geom.ComputeNodalGradientRef(Un, Unode);
  }
  
  template<class Dimension>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeNodalGradientRef(const Vector<Complex_wp>& Un,
			    Vector<TinyVector<Complex_wp, Dimension::dim_N> >& Unode) const
  {
    return this->elt_geom.ComputeNodalGradientRef(Un, Unode);
  }
    
  template<class Dimension> template<class T>
  inline void ElementGeomReferenceContainer_Base<Dimension>
  ::ComputeNodalGradient(const SetMatrices<Dimension>& MatricesElem,
			 const Vector<Vector<T> >& U, Vector<Vector<T> >& gradU) const
  {
    return this->elt_geom.ComputeNodalGradient(MatricesElem, U, gradU);
  }

} // namespace Montjoie

#define MONTJOIE_FILE_ELEMENT_GEOM_REFERENCE_INLINE_CXX
#endif
