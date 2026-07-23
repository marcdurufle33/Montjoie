#ifndef MONTJOIE_FILE_TRIANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TriangleHdivOptimalHpFirstFamily::TriangleHdivOptimalHpFirstFamily()
    : TriangleReference<3>()
  {
    // dof, nodal and quadrature points are all different
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number the mesh
  void TriangleHdivOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TriangleReference<3>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 0);
    nmap.SetNbDofTriangle(order, nb_dof_loc - 3*order);
    
    nmap.SetEqualEdgesDofSymmetry(order, order+1);
    
    // FindHcurlSignEdge(*this, nmap);
  }


  size_t TriangleHdivOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = TriangleReference<3>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += NumDofsX.GetMemorySize();
    taille += NumDofsY.GetMemorySize();
    return taille;
  }
  
  
  //! construction of finite element
  void TriangleHdivOptimalHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    TriangleReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
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
  void TriangleHdivOptimalHpFirstFamily::ConstructFunctions()
  {
    // dofs on the three edges and inside dofs (two by inside nodes)
    nb_dof_loc = order*3 + (order-1)*order;
    nb_dof_boundaries = 3*order;

    VectReal_wp points_dof1d; VectR2 points_dof2d;        
    Fb_geom.ConstructLobattoPoints(order, 0, points_dof1d, points_dof2d);

    points_dof1d = this->Points1D();
    points_dof2d = this->PointsND();
    
    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    NumDofsX.Reallocate(order-1, order-1);
    NumDofsY.Reallocate(order-1, order-1);
    NumDofsX.Fill(-1); NumDofsY.Fill(-1);
    int offset = 3*order;
    for (int diag = 0; diag <= order-2; diag++)
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag-i; j++)
	  if (i+j == diag)
	    {
	      NumDofsX(i, j) = offset++;
	      NumDofsY(i, j) = offset++;
	    }
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] contrib dof components
   */
  template<class Vector1, class Vector2>
  void TriangleHdivOptimalHpFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    ProjectQuadratureToDofRef(feval, contrib);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TriangleHdivOptimalHpFirstFamily::
  ComputeValuesPhiRef(const R2& point_loc, VectR2& phi) const
  {
    phi.Reallocate(nb_dof_loc);
	
    Real_wp x, y;
    x = point_loc(0);
    y = point_loc(1);
    
    VectReal_wp P1, P2, P3;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*x - 1.0 + y, P1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y-x, P2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, -2.0*y + 1.0 - x, P3);
    
    int node = 0;
    // edges
    // edge 0
    for (int i = 0; i < order; i++)
      {	
        phi(node)(0) = x*P1(i);
        phi(node)(1) = (y-1.0)*P1(i);
        node++;
      }
    
    // edge 2
    for (int i = 0; i < order; i++)
      {	
        phi(node)(0) = x*P2(i);
        phi(node)(1) = y*P2(i);
        node++;
      }
    
    // edge 3
    for (int i = 0; i < order; i++)
      {	
        phi(node)(0) = (x-1.0)*P3(i);
        phi(node)(1) = y*P3(i);
        node++;
      }
    
    // inside
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        {	
	  node = NumDofsX(i, j);
          phi(node)(0) = (x-1.0)*x*P1(i)*P3(j);
          phi(node)(1) = y*x*P1(i)*P3(j);
          if (j%2 == 1)
            Mlt(-1.0, phi(node));
          
	  node = NumDofsY(i, j);
          phi(node)(0) = x*y*P1(i)*P3(j);
          phi(node)(1) = (y-1.0)*y*P1(i)*P3(j);
          if (j%2 == 1)
            Mlt(-1.0, phi(node));
        }
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void TriangleHdivOptimalHpFirstFamily
  ::ComputeDivPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    Real_wp x, y;
    x = point_loc(0);
    y = point_loc(1);
    
    VectReal_wp P1, P2, P3, dP1, dP2, dP3;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*x - 1.0 + y, P1, dP1);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y-x, P2, dP2);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, -2.0*y + 1.0 - x, P3, dP3);
    
    int node = 0;
    // edges
    // edge 0
    for (int i = 0; i < order; i++)
      {	
        res(node) = 2.0*P1(i) + dP1(i)*(2.0*x - (1.0-y));
        node++;
      }
    
    // edge 2
    for (int i = 0; i < order; i++)
      {	
        res(node) = 2.0*P2(i) + dP2(i)*(-x + y);
        node++;
      }
    
    // edge 3
    for (int i = 0; i < order; i++)
      {	
        res(node) = 2.0*P3(i) + dP3(i)*(-(x-1.0) - 2.0*y);
        node++;
      }
    
    // inside
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        { 
	  node = NumDofsX(i, j);
          res(node) = 2.0*x*P1(i)*P3(j)
            + (x - 1.0)*( P1(i)*P3(j) + x*(2.0*dP1(i)*P3(j) - P1(i)*dP3(j)))
            + y*x*(dP1(i)*P3(j) - 2.0*P1(i)*dP3(j));
          
          if (j%2 == 1)
            res(node) = -res(node);
          
	  node = NumDofsY(i, j);
          res(node) = 2.0*y*P1(i)*P3(j) 
            + x*(y*(2.0*dP1(i)*P3(j) - P1(i)*dP3(j)))
            + (y - 1.0)*( P1(i)*P3(j) + y*(dP1(i)*P3(j) - 2.0*P1(i)*dP3(j)));
          
          if (j%2 == 1)
            res(node) = -res(node);
        }
    
  }


  //! displays informations about class TriangleHcurlOptimalHpFirstFamily
  ostream& operator <<(ostream& out, const TriangleHdivOptimalHpFirstFamily& e)
  {
    out << static_cast<const TriangleReference<3>& >(e);
    out << "Triangle of Nedelec's first family " << endl;
    out << "Number of local degrees of freedom " << e.GetNbDof() << endl;
    return out;
  }

}
  
#define MONTJOIE_FILE_TRIANGLE_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
