#ifndef MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  QuadrangleHdivOptimalHpFirstFamily::QuadrangleHdivOptimalHpFirstFamily()
    : QuadrangleReference<3>()
  {
    // dof, nodal and quadrature points are all different
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number the mesh
  void QuadrangleHdivOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<3>::ConstructNumberMap(nmap, dg);

    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, nb_dof_loc - 4*order);
    
    nmap.SetEqualEdgesDofSymmetry(order, order+1);
  }
  

  size_t QuadrangleHdivOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<3>::GetMemorySize();
    taille += NumDofsX.GetMemorySize();
    taille += NumDofsY.GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += JacobiPolynom.GetMemorySize();
    return taille;
  }
  
      
  //! construction of finite element
  void QuadrangleHdivOptimalHpFirstFamily::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			 int rsurf, int type_surf)
  {
    QuadrangleReference<3>::
      ConstructFiniteElement(r, rgeom, r+1, Globatto<Real_wp>::QUADRATURE_GAUSS, r);
    
    ConstructFunctions();
    
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructHdivElementaryMatrix();
    
    FindDofsOnEdge();
    
    delete this->element_surface;
    EdgeHierarchicReference* edge = new EdgeHierarchicReference();
    edge->ConstructFiniteElement(r-1, r, r, type_surf, EdgeHierarchic::LEGENDRE);
    this->element_surface = edge;

  }
  
  
  //! construction of basis functions
  void QuadrangleHdivOptimalHpFirstFamily::ConstructFunctions()
  {
    // dofs on the three edges and inside dofs (two by inside nodes)
    nb_dof_loc = order*4 + 2*order*order;
    nb_dof_boundaries = 4*order;
    
    Globatto<Real_wp> lob;
    lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);

    VectReal_wp points_dof1d; VectR2 points_dof2d;    
    points_dof1d = this->Points1D();
    points_dof2d = this->PointsND();
    
    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    GetJacobiPolynomial(JacobiPolynom, order, Real_wp(1), Real_wp(1));
    
    NumDofsX.Reallocate(order, order);
    NumDofsY.Reallocate(order, order);
    NumDofsX.Fill(-1); NumDofsY.Fill(-1);
    int offset = 4*order;
    for (int diag = 0; diag <= order-1; diag++)
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  if (max(i, j) == diag)
	    {
	      NumDofsX(i, j) = offset++;
	      NumDofsY(i, j) = offset++;
	    }
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] res dof components
   */
  template<class Vector1, class Vector2>
  void QuadrangleHdivOptimalHpFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    ProjectQuadratureToDofRef(feval, contrib);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point_loc
  */
  void QuadrangleHdivOptimalHpFirstFamily::
  ComputeValuesPhiRef(const R2& point_loc, VectR2& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    FillZero(phi);
    
    Real_wp x, y;
    x = point_loc(0);
    y = point_loc(1);
    
    // evaluating P_i^{0, 0} and P_i^{1, 1}
    VectReal_wp Px, Py, Jx, Jy;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*y-1, Py);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*x-1, Jx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*y-1, Jy);
    
    int node = 0;
    // edge 0
    for (int i = 0; i < order; i++)
      {	
        phi(node)(1) = -(1-y)*Px(i);
        node++;
      }
    
    // edge 1
    for (int i = 0; i < order; i++)
      {	
        phi(node)(0) = x*Py(i);
        node++;
      }
    
    // edge 2
    for (int i = 0; i < order; i++)
      {	
        phi(node)(1) = y*Px(i);        
        if (i%2 == 1)
          phi(node)(1) = -phi(node)(1);
        
        node++;
      }
    
    // edge 3
    for (int i = 0; i < order; i++)
      {	
        phi(node)(0) = (x-1.0)*Py(i);
        if (i%2 == 1)
          phi(node)(0) = -phi(node)(0);
        
        node++;
      }
    
    // inside
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {	
	  node = NumDofsY(i, j);
          phi(node)(1) = y*(1.0-y)*Px(i)*Jy(j);

	  node = NumDofsX(i, j);
          phi(node)(0) = x*(1.0-x)*Jx(i)*Py(j);
	}
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void QuadrangleHdivOptimalHpFirstFamily
  ::ComputeDivPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    Real_wp x, y;
    x = point_loc(0);
    y = point_loc(1);
    
    // evaluating P_i^{0, 0} and P_i^{1, 1}
    VectReal_wp Px, Py, Jx, Jy, dJx, dJy;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*y-1, Py);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*x-1, Jx, dJx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*y-1, Jy, dJy);
    
    int node = 0;
    // edge 0
    for (int i = 0; i < order; i++)
      {	
        res(node) = Px(i);
        node++;
      }
    
    // edge 1
    for (int i = 0; i < order; i++)
      {	
        res(node) = Py(i);
        node++;
      }
    
    // edge 2
    for (int i = 0; i < order; i++)
      {	
        res(node) = Px(i);
        if (i%2 == 1)
          res(node) = -res(node);
        
        node++;
      }
    
    // edge 3
    for (int i = 0; i < order; i++)
      {	
        res(node) = Py(i);
        if (i%2 == 1)
          res(node) = -res(node);
        
        node++;
      }
    
    // inside
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
        {	
	  node = NumDofsY(i, j);
          res(node) = -Px(i)*( (2.0*y - 1.0)*Jy(j) - 2.0*y*(1.0-y)*dJy(j));
	  
	  node = NumDofsX(i, j);
          res(node) = Py(j)*( -(2.0*x - 1.0)*Jx(i) + 2.0*x*(1.0-x)*dJx(i));
        }
  }
  
  
  //! displays informations about class 
  ostream& operator <<(ostream& out, const QuadrangleHdivOptimalHpFirstFamily& e)
  {
    out<<static_cast<const QuadrangleReference<3>& >(e);
    out<<"Quadrangle of Nedelec's first family "<<endl;
    out<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    return out;
  }

}
  
#define MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
