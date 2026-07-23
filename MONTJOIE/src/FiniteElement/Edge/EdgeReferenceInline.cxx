#ifndef MONTJOIE_FILE_EDGE_REFERENCE_INLINE_CXX

namespace Montjoie
{
  
  /************************************
   * ElementGeomReference<Dimension1> *
   ************************************/
  
  
  //! Destructor
  inline ElementGeomReference<Dimension1>::~ElementGeomReference<Dimension1>()
  {
  }
  
  
  //! returns true if the mass matrix is diagonal
  inline bool ElementGeomReference<Dimension1>::LumpedMassMatrix() const
  {
    return mass_lumping;
  }
  
  
  //! returns 0
  inline int ElementGeomReference<Dimension1>::GetHybridType() const
  {
    return 0;
  }
  
  
  //! returns the order of approximation
  inline int ElementGeomReference<Dimension1>::GetOrder() const
  {
    return order;
  }
  
  
  //! returns the order of approximation
  inline int ElementGeomReference<Dimension1>::GetGeometryOrder() const
  {
    return order_geom;
  }
  

  inline const VectReal_wp& ElementGeomReference<Dimension1>::PointsDof() const
  {
    return points_dof; 
  }

  
  //! returns the number of quadrature points
  inline int ElementGeomReference<Dimension1>::GetNbPointsQuadratureInside() const
  {
    return nb_points;
  }
  
  
  //! returns the number of degrees of freedom
  inline int ElementGeomReference<Dimension1>::GetNbDof() const
  {
    return nb_dof_loc;
  }
    
  
  //! returns the number of nodal points
  inline int ElementGeomReference<Dimension1>::GetNbPointsNodalElt() const
  {
    return nb_dof_loc;
  }
  
  
  //! allocates a projector (to project from nodal points to any set of points)
  inline FiniteElementProjector* ElementGeomReference<Dimension1>
  ::GetNewNodalInterpolation() const
  {
    return new DenseProjector<Dimension1>();
  }
  
  
  //! evaluates nodal shape functions at point x
  inline void ElementGeomReference<Dimension1>
  ::ComputeValuesPhiNodalRef(const Real_wp& x, VectReal_wp& res) const
  {
    ComputeValuesPhiRef(x, res);
  }
  
  
  //! returns stiffness matrix
  inline const Matrix<Real_wp, Symmetric, RowSymPacked>&
  ElementGeomReference<Dimension1>::GetStiffnessMatrix() const
  {
    return stiffness_matrix;
  }


  //! returns mass matrix
  inline const Matrix<Real_wp, Symmetric, RowSymPacked>&
  ElementGeomReference<Dimension1>::GetMassMatrix() const
  {
    return mass_matrix;
  }

  
  //! returns gradient matrix
  inline const Matrix<Real_wp>& ElementGeomReference<Dimension1>::GetGradientMatrix() const
  {
    return gradient_matrix;
  }

  
  /******************
   * EdgeHierarchic *
   ******************/
  

  inline Real_wp EdgeHierarchic::GetValuePhi1D(int i, const Real_wp& pointloc) const
  {
    abort();
    return Real_wp(0);
  }


  inline Real_wp EdgeHierarchic::GetGradientPhi1D(int i, const Real_wp& pointloc) const
  {
    abort();
    return Real_wp(0);
  }


  inline void EdgeHierarchic
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C,
			  const Real_wp& D, const Real_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
    
  
  inline void EdgeHierarchic
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C,
			  const Complex_wp& D, const Complex_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
    
  
  inline void EdgeHierarchic
  ::AddVariableElemMatrix(int off_row, int off_col, const VectReal_wp& mass,
			  const VectReal_wp& C, const VectReal_wp& D, const VectReal_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& mat) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
  

  inline void EdgeHierarchic
  ::AddVariableElemMatrix(int off_row, int off_col, const VectComplex_wp& mass,
			  const VectComplex_wp& C, const VectComplex_wp& D, const VectComplex_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& mat) const
  {
    cout << "Not implemented" << endl;
    abort();
  }

  
  /************************************
   * ElementReference_Dim<Dimension1> *
   ************************************/
  
  
  inline ElementReference_Dim<Dimension1>
  ::ElementReference_Dim(ElementGeomReference<Dimension1>& elt) : ElementReference_Base(elt.Weights), elt_geom(elt)
  {
  }
    
  
  inline ElementReference_Dim<Dimension1>::~ElementReference_Dim()
  {
  }
  
  
  inline int ElementReference_Dim<Dimension1>::GetGeometryOrder() const
  {
    return elt_geom.GetGeometryOrder();
  }

    
  inline size_t ElementReference_Dim<Dimension1>::GetMemorySize() const
  {
    return elt_geom.GetMemorySize() + sizeof(*this);
  }

  
  inline void ElementReference_Dim<Dimension1>
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad, int type_func)
  {
    elt_geom.ConstructFiniteElement(r, rgeom, rquad, type_quad, type_func);
    this->order = elt_geom.GetOrder();
    this->order_quad = elt_geom.GetOrder();
    this->nb_points_quadrature_inside = elt_geom.GetNbPointsQuadratureInside();
    this->nb_dof_loc = elt_geom.GetNbDof();
    this->mass_lumping = elt_geom.LumpedMassMatrix();
  }
  
  
  inline void ElementReference_Dim<Dimension1>::SetDofPoints(const VectReal_wp& pts)
  {
    elt_geom.SetDofPoints(pts);
  }
  
  
  inline const VectReal_wp& ElementReference_Dim<Dimension1>::PointsDof() const
  {
    return elt_geom.PointsDof(); 
  }

  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeValuesPhiRef(const Real_wp& pointloc, VectReal_wp& res) const
  {
    elt_geom.ComputeValuesPhiRef(pointloc, res);
  }
  
  
  inline Real_wp ElementReference_Dim<Dimension1>::GetValuePhi1D(int i, const Real_wp& pointloc) const
  {
    return elt_geom.GetValuePhi1D(i, pointloc);
  }


  inline Real_wp ElementReference_Dim<Dimension1>::GetGradientPhi1D(int i, const Real_wp& pointloc) const
  {
    return elt_geom.GetGradientPhi1D(i, pointloc);
  }
  
  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeGradientPhiRef(const Real_wp& pointloc, VectReal_wp& res) const
  {
    elt_geom.ComputeGradientPhiRef(pointloc, res);
  }

  
  inline void ElementReference_Dim<Dimension1>::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    elt_geom.GetValueSinglePhiQuadrature(i, phi);
  }
  
  inline const VectReal_wp& ElementReference_Dim<Dimension1>::Points() const
  {
    return elt_geom.Points;
  }
  
  
  inline const VectReal_wp& ElementReference_Dim<Dimension1>::Weights() const
  {
    return elt_geom.Weights;
  }

  
  inline const Real_wp& ElementReference_Dim<Dimension1>::Points(int i) const
  {
    return elt_geom.Points(i);
  }
  
  
  inline const Real_wp& ElementReference_Dim<Dimension1>::Weights(int i) const
  {
    return elt_geom.Weights(i);
  }

  
  inline const Real_wp& ElementReference_Dim<Dimension1>::WeightsND(int i) const
  {
    return elt_geom.Weights(i);
  }  


  inline const VectReal_wp& ElementReference_Dim<Dimension1>::WeightsND() const
  {
    return elt_geom.Weights;
  }  
  
  
  inline const Matrix<Real_wp, Symmetric, RowSymPacked>&
  ElementReference_Dim<Dimension1>::GetStiffnessMatrix() const
  {
    return elt_geom.GetStiffnessMatrix();
  }
  
  
  inline const Matrix<Real_wp>& ElementReference_Dim<Dimension1>::GetGradientMatrix() const
  {
    return elt_geom.GetGradientMatrix();
  }

  
  inline void ElementReference_Dim<Dimension1>
  ::GetPolynomialFunctions(Vector<UnivariatePolynomial<Real_wp> >& Phi,
			   Vector<UnivariatePolynomial<Real_wp> >& dPhi) const
  {
    elt_geom.GetPolynomialFunctions(Phi, dPhi);
  }
  

  inline void ElementReference_Dim<Dimension1>::FjElem(const VectReal_wp& s, VectReal_wp& points) const
  {
    return elt_geom.FjElem(s, points);
  }

  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeIntegralRef(const VectReal_wp& u, VectReal_wp& v) const
  {
    elt_geom.ComputeIntegralRef(u, v);
  }
  
  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeIntegralGradientRef(const VectReal_wp& u, VectReal_wp& v) const
  {
    elt_geom.ComputeIntegralGradientRef(u, v);
  }   
  
  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeIntegralRef(const VectComplex_wp& u, VectComplex_wp& v) const
  {
    elt_geom.ComputeIntegralRef(u, v);
  }
  
    
  inline void ElementReference_Dim<Dimension1>
  ::ComputeIntegralGradientRef(const VectComplex_wp& u, VectComplex_wp& v) const
  {
    elt_geom.ComputeIntegralGradientRef(u, v);
  }   


  inline void ElementReference_Dim<Dimension1>::ConstructNumberMap(NumberMap& nmap, int) const
  {
  }
  
  
  inline void ElementReference_Dim<Dimension1>
  ::ApplyCh(const VectReal_wp& u, VectReal_wp& v) const
  {
    elt_geom.ApplyCh(u, v);
  }


  inline void ElementReference_Dim<Dimension1>
  ::ApplyCh(const VectComplex_wp& u, VectComplex_wp& v) const
  {
    elt_geom.ApplyCh(u, v);
  }

  
  inline void ElementReference_Dim<Dimension1>
  ::ApplyChTranspose(const VectReal_wp& u, VectReal_wp& v) const
  {
    elt_geom.ApplyChTranspose(u, v);
  }


  inline void ElementReference_Dim<Dimension1>
  ::ApplyChTranspose(const VectComplex_wp& u, VectComplex_wp& v) const
  {
    elt_geom.ApplyChTranspose(u, v);
  }

  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    elt_geom.ComputeProjectionDofRef(feval, contrib);
  }
  
  
  inline void ElementReference_Dim<Dimension1>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    elt_geom.ComputeProjectionDofRef(feval, contrib);
  }


  inline void ElementReference_Dim<Dimension1>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int) const
  {
  }


  inline void ElementReference_Dim<Dimension1>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int) const
  {
  }
  

  inline void ElementReference_Dim<Dimension1>
  ::AddConstantElemMatrix(int m, int n, const Real_wp& mass, const Real_wp& C,
			  const Real_wp& D, const Real_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    elt_geom.AddConstantElemMatrix(m, n, mass, C, D, E, null_term, A);
  }
   
  
  inline void ElementReference_Dim<Dimension1>
  ::AddConstantElemMatrix(int m, int n, const Complex_wp& mass, const Complex_wp& C,
			  const Complex_wp& D, const Complex_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    elt_geom.AddConstantElemMatrix(m, n, mass, C, D, E, null_term, A);
  }
    
  
  inline void ElementReference_Dim<Dimension1>
  ::AddVariableElemMatrix(int m, int n, const VectReal_wp& mass,
			  const VectReal_wp& C, const VectReal_wp& D, const VectReal_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Real_wp>& A) const
  {
    elt_geom.AddVariableElemMatrix(m, n, mass, C, D, E, null_term, A);
  }
  

  inline void ElementReference_Dim<Dimension1>
  ::AddVariableElemMatrix(int m, int n, const VectComplex_wp& mass,
			  const VectComplex_wp& C, const VectComplex_wp& D, const VectComplex_wp& E,
			  const TinyVector<bool, 4>& null_term, VirtualMatrix<Complex_wp>& A) const
  {
    elt_geom.AddVariableElemMatrix(m, n, mass, C, D, E, null_term, A);
  }
  

  /***********************************************
   * EdgeLobattoReference and other constructors *
   ***********************************************/
  
  
  
  inline EdgeLobattoReference::EdgeLobattoReference() : ElementReference<Dimension1, 1>(edge)
  {
  }


  inline EdgeRadauReference::EdgeRadauReference() : ElementReference<Dimension1, 1>(edge)
  {
  }
  
  
  inline EdgeGaussReference::EdgeGaussReference() : ElementReference<Dimension1, 1>(edge)
  {
  }


  
  inline EdgeHierarchicReference::EdgeHierarchicReference() : ElementReference<Dimension1, 1>(edge)
  {
  }
  
}

#define MONTJOIE_FILE_EDGE_REFERENCE_INLINE_CXX
#endif
