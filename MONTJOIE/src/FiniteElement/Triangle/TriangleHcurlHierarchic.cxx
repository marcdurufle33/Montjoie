#ifndef MONTJOIE_FILE_TRIANGLE_HCURL_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  TriangleHcurlHierarchic::TriangleHcurlHierarchic() : TriangleReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    nb_dof_based_edges = 0;
  }
  
  
  //! how to number a mesh
  void TriangleHcurlHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TriangleReference<2>::ConstructNumberMap(nmap, dg);

    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order+1);
    nmap.SetNbDofTriangle(order, nb_dof_loc-3*(order+1));
    nmap.SetNbDofQuadrangle(order, 0);
    
    nmap.SetEqualEdgesDofSymmetry(order, order+1);
    
  }


  size_t TriangleHcurlHierarchic::GetMemorySize() const
  {
    size_t taille = TriangleReference<2>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    return taille;
  }
  
  
  //! constructing finite element
  void TriangleHcurlHierarchic
  ::ConstructFiniteElement(int order_, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    TriangleReference<2>::ConstructFiniteElement(order_, rgeom, rquad, type_quad);
    
    order = order_;
    ConstructFunctions();

    this->Fb_geom.ComputeCoefficientTransformation();
    ConstructHcurlElementaryMatrix();

    delete this->element_surface;
    EdgeHierarchicReference* edge = new EdgeHierarchicReference();
    edge->ConstructFiniteElement(order, order, rquad, type_surf,
				 EdgeHierarchic::LEGENDRE_COMBINED);
    
    this->element_surface = edge;    
  }
  
  
  //! construction of basis functions
  void TriangleHcurlHierarchic::ConstructFunctions()
  {
    nb_dof_boundaries = 3*(order+1);
    int nb_dof_based_edges = nb_dof_boundaries;
    if (order >= 2)
      nb_dof_based_edges = nb_dof_boundaries + 3*(order-1);

    if (order == 0)
      nb_dof_loc = 3;
    else if (order == 1)
      nb_dof_loc = 6;
    else
      nb_dof_loc = (order+1)*(order+2);
    
    CoordinateDofs.Reallocate(nb_dof_loc,3);
    EdgesDof.Reallocate(order+1, 3);
    for (int num_dof = 0; num_dof < nb_dof_boundaries; num_dof++)
      {
	// c'est une fonction de base associée a une arete
        int num_loc_edge = (num_dof/(order+1));
	int k = num_dof - num_loc_edge*(order+1);
	EdgesDof(k, num_loc_edge) = num_dof;
	CoordinateDofs(num_dof,0) = num_loc_edge;
	CoordinateDofs(num_dof,1) = k;
	CoordinateDofs(num_dof,2) = -2;
      }
    
    for (int num_dof = nb_dof_boundaries; num_dof < nb_dof_based_edges; num_dof++)
      {
	// c'est une fonction de base intérieure associée à une arete
        int num_loc_edge = (num_dof-nb_dof_boundaries)/(order-1);
	int k = num_dof - nb_dof_boundaries - num_loc_edge*(order-1);
	CoordinateDofs(num_dof,0) = num_loc_edge;
	CoordinateDofs(num_dof,1) = k;
	CoordinateDofs(num_dof,2) = -1;
      }
    
    // on affecte les "coordonnées" aux fonctions intérieures
    int num_dof = nb_dof_based_edges;
    for (int k = 0; k <= order-3; k++)
      for (int l = 0; l <= order-3-k; l++)
	{
	  CoordinateDofs(num_dof, 0) = 0;
	  CoordinateDofs(num_dof, 1) = k;
	  CoordinateDofs(num_dof, 2) = l;
	  CoordinateDofs(num_dof+1, 0) = 1;
	  CoordinateDofs(num_dof+1, 1) = k;
	  CoordinateDofs(num_dof+1, 2) = l;
	  num_dof += 2;
	}

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    // on calcule les polynômes de Legendre
    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));

    this->SetPointsDofND(this->PointsND());    
    this->SetPointsDof1D(this->Points1D());

    this->nb_points_dof_inside = nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TriangleHcurlHierarchic::
  ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(R2(0,0));
    Real_wp x = point_loc(0); Real_wp y = point_loc(1);
    
    int num_dof = 0;
    // first edge
    res(num_dof)(0) = 1.0-y;
    res(num_dof)(1) = x;
    num_dof++;
    
    if (order >= 1)
      {
	res(num_dof)(0) = 2.0*x+y-1.0;
        res(num_dof)(1) = x;
        num_dof++;	
      }

        
    Real_wp coef_phi0, coef_phi1;    
    R2 phi0, phi1;
    Real_wp param_edge1 = 2.0*x+y-1.0; VectReal_wp Pn1;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, param_edge1, Pn1);
    
    // first edge
    phi0(0) = 1.0-y; phi0(1) = x;
    phi1(0) = 2.0*x + y -1.0; phi1(1) = x;
    for (int k = 2; k <= order; k++)
      {
	coef_phi0 = -(k-1)*Pn1(k-2)/k ;
	coef_phi1 = (2*k-1)*Pn1(k-1)/k;
	Add(coef_phi0, phi0, res(num_dof));
	Add(coef_phi1, phi1, res(num_dof));
	num_dof++;
      }    

    
    // second edge
    res(num_dof)(0) = -y;
    res(num_dof)(1) = x;
    num_dof++;
    if (order >= 1)
      {
	res(num_dof)(0) = -y;
        res(num_dof)(1) = -x;
        num_dof++;
      }

    phi0(0) = -y; phi0(1) = x;
    phi1(0) = -y; phi1(1) = -x;
    Real_wp param_edge2 = y-x; VectReal_wp Pn2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, param_edge2, Pn2);
    for (int k = 2; k <= order; k++)
      {
	coef_phi0 = -(k-1)*Pn2(k-2)/k ;
	coef_phi1 = (2*k-1)*Pn2(k-1)/k;
	Add(coef_phi0, phi0, res(num_dof));
	Add(coef_phi1, phi1, res(num_dof)); num_dof++;	
      }
    
    // third edge
    res(num_dof)(0) = -y;
    res(num_dof)(1) = x-1.0;
    num_dof++;
    if (order >= 1)
      {
	res(num_dof)(0) = y;
        res(num_dof)(1) = 2.0*y+x-1.0;
        num_dof++;
      }    

    phi0(0) = -y; phi0(1) = x-1.0;
    phi1(0) = y; phi1(1) = 2.0*y+x-1.0;
    Real_wp param_edge3 = 1.0-2.0*y-x; VectReal_wp Pn3;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, param_edge3, Pn3);
    for (int k = 2; k <= order; k++)
      {
	coef_phi0 = -(k-1)*Pn3(k-2)/k ;
	coef_phi1 = (2*k-1)*Pn3(k-1)/k;
	Add(coef_phi0, phi0, res(num_dof));
	Add(coef_phi1, phi1, res(num_dof)); num_dof++;	
      }
    
    // edge-based functions on first edge
    Real_wp coef = x*(1.0-x-y);
    for (int k = 0; k < order-1; k++)
      {
	res(num_dof) = this->NormaleLoc(0);
	Mlt(coef*Pn1(k), res(num_dof)); num_dof++;
      }
    
    // edge-based functions on second edge
    coef = x*y;
    for (int k = 0; k < order-1; k++)
      {
	res(num_dof) = this->NormaleLoc(1);
	Mlt(coef*Pn2(k), res(num_dof)); num_dof++;
      }
	
    // edge-based functions on third edge
    coef = y*(1.0-x-y);
    for (int k = 0; k < order-1; k++)
      {
	res(num_dof) = this->NormaleLoc(2);
	Mlt(coef*Pn3(k), res(num_dof)); num_dof++;
      }

    // interior functions
    // \phi_{k,l}^1 = x y (1-x-y) L_k ( 2x +y -1) L_l (1-x-2y) e_x
    // \phi_{k,l}^2 = x y (1-x-y) L_k ( 2x +y -1) L_l (1-x-2y) e_y
    coef = x*y*(1.0-x-y);
    for (int k = 0; k < order-2; k++)
      for (int l = 0; l < order-2-k; l++)
	{
	  res(num_dof)(0) = coef*Pn1(k)*Pn3(l);
	  res(num_dof+1)(1) = res(num_dof)(0);
	  num_dof += 2;
      }
      
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void TriangleHcurlHierarchic::ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    Real_wp x = point_loc(0); Real_wp y = point_loc(1);
    int num_dof = 0;
    
    // pour les ordres supérieurs, on utilise les formules suivantes :
    // \phi_k^1 = (2k-1)/k L_{k-1}(2.0*x+y-1.0) \phi_1^1
    //           - (k-1)/k L_{k-2}(2.0*x+y-1.0) \phi_0^1
    // \phi_k^2 = (2k-1)/k L_{k-1}(y-x) \phi_1^2
    //            - (k-1)/k L_{k-2}(y-x) \phi_0^2
    // \phi_k^3 = (2k-1)/k L_{k-1}(1.0-x-2.0*y) \phi_1^3
    //             - (k-1)/k L_{k-2}(1.0-x-2.0*y) \phi_0^3
    // où \phi_0^1 = [1-y; x]   \phi_0^2 = [-y; x]  \phi_0^3 = [-y; x-1]
    // et \phi_1^1 = [2*x+y-1; x]   \phi_1^2 = [-y; -x]  \phi_1^3 = [y; 2*y+x-1]
    // L_k est le polynome de Legendre d'ordre k
    
    // on sait que rot ( \lambda \vec(u) )
    //          = \nabla \lambda \times u + \lambda rot y
    // d'où
    // rot(\phi_k^1) = (2k-1)/k (1-y) L'_{k-1} (2x+y-1)
    //                 - (k-1)/k (2x+y-1) L'_{k-2} (2x+y-1)
    //                 -2(k-1)/k L_{k-2}(2x+y-1)
    // rot(\phi_k^2) = (2k-1)/k (x+y) L'_{k-1} (y-x)
    //                 - (k-1)/k (y-x) L'_{k-2} (y-x)
    //                 -2(k-1)/k L_{k-2}(y-x)
    // rot(\phi_k^3) = (2k-1)/k (1-x) L'_{k-1} (1-2y-x)
    //                 - (k-1)/k (1-x-2y) L'_{k-2} (1-x-2y)
    //                 -2(k-1)/k L_{k-2}(1-2y-x)
    
    // first edge
    res(num_dof) = 2.0; num_dof++;
    res(num_dof) = 0.0; num_dof++;
    
    Real_wp param_edge1 = 2.0*x+y-1.0; VectReal_wp Pn1, dPn1;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, param_edge1, Pn1, dPn1);
    for (int k = 2; k <= order; k++)
      {
	Real_wp coef_der_lkm1 = (1.0-y)*(2*k-1)/k;
	Real_wp coef_der_lkm2 = (2.0*x+y-1.0)*(k-1)/k;
	Real_wp coef_lkm2 = 2.0*(k-1)/k;
	coef_der_lkm1 *= dPn1(k-1);
	coef_der_lkm2 *= dPn1(k-2);
	coef_lkm2 *= Pn1(k-2);
	
	res(num_dof) = coef_der_lkm1 - coef_der_lkm2 - coef_lkm2; num_dof++;
      }
    
    // second edge
    res(num_dof) = 2.0; num_dof++;
    res(num_dof) = 0.0; num_dof++;
    Real_wp param_edge2 = y-x; VectReal_wp Pn2, dPn2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, param_edge2, Pn2, dPn2);
    for (int k = 2; k <= order; k++)
      {
	Real_wp coef_der_lkm1 = (x+y)*(2*k-1)/k;
	Real_wp coef_der_lkm2 = (y-x)*(k-1)/k;
	Real_wp coef_lkm2 = 2.0*(k-1)/k;
	coef_der_lkm1 *= dPn2(k-1);
	coef_der_lkm2 *= dPn2(k-2);
	coef_lkm2 *= Pn2(k-2);
	
	res(num_dof) = coef_der_lkm1 - coef_der_lkm2 - coef_lkm2; num_dof++;
      }
    
    // third edge
    res(num_dof) = 2.0; num_dof++;
    res(num_dof) = 0.0; num_dof++;
    
    Real_wp param_edge3 = 1.0-2.0*y-x; VectReal_wp Pn3, dPn3;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, param_edge3, Pn3, dPn3);
    for (int k = 2; k <= order; k++)
      {
	Real_wp coef_der_lkm1 = (1.0-x)*(2*k-1)/k;
	Real_wp coef_der_lkm2 = (1.0-2.0*y-x)*(k-1)/k;
	Real_wp coef_lkm2 = 2.0*(k-1)/k;
	coef_der_lkm1 *= dPn3(k-1);
	coef_der_lkm2 *= dPn3(k-2);
	coef_lkm2 *= Pn3(k-2);
	
	res(num_dof) = coef_der_lkm1 - coef_der_lkm2 - coef_lkm2; num_dof++;
      }

    // fonctions de base basées sur les aretes
    // \psi_k^1 = \lambda_1 \lambda_2 L_k(\lambda_2 - \lambda_1) \n_1
    //          = (1-x-y) x L_k (2*x+y-1) [0;-1]
    // \psi_k^2 = \lambda_2 \lambda_3 L_k(\lambda_3 - \lambda_2) \n_2
    //          = x y       L_k (y-x)     [1; 1]
    // \psi_k^3 = \lambda_3 \lambda_1 L_k(\lambda_1 - \lambda_3) \n_3
    //          = y (1-x-y) L_k(1-x-2*y)  [-1;0]
    // L_k est le polynôme de Legendre d'ordre k
    // le rotationnel s'écrit alors
    // rot(\psi_k^1) = (2x+y-1) L_k(2x+y-1) - 2x(1-x-y) L'_k (2x+y-1)
    // rot(\psi_k^2) = (y-x) L_k(y-x) - 2xy L'_k (y-x)
    // rot(\psi_k^3) = (1-x-2y) L_k(1-x-2y) - 2y(1-x-y) L'_k (1-x-2y)
    
    // first edge
    Real_wp coef_der;
    for (int k = 0; k < order-1; k++)
      {
	coef_der = -2.0*x*(1.0-x-y);
	res(num_dof) = param_edge1*Pn1(k) + coef_der*dPn1(k); num_dof++;
      }
    
    // second edge
    for (int k = 0; k < order-1; k++)
      {
	coef_der = -2.0*x*y;
	res(num_dof) = param_edge2*Pn2(k) + coef_der*dPn2(k); num_dof++;
      }
    
    // third edge
    for (int k = 0; k < order-1; k++)
      {
	coef_der = -2.0*y*(1.0-x-y);
	res(num_dof) = param_edge3*Pn3(k) + coef_der*dPn3(k); num_dof++;
      }
    
    // fonctions de base complètement intérieures
    // \phi_{k,l}^1 = x y (1-x-y) L_k ( 2x +y -1) L_l (1-x-2y) e_x
    // \phi_{k,l}^2 = x y (1-x-y) L_k ( 2x +y -1) L_l (1-x-2y) e_y
    Real_wp coef = x*y*(1.0-x-y);
    for (int k = 0; k < order-2; k++)
      for (int l = 0; l < order-2-k; l++)
	{
	  res(num_dof) = -x*(1.0-x-2.0*y)*Pn1(k)*Pn3(l)
            - coef * (-2.0*Pn1(k)*dPn3(l) + dPn1(k)*Pn3(l));
          
          num_dof++;
	  
          res(num_dof) = y*(1.0-2.0*x-y)*Pn1(k)*Pn3(l)
            + coef * (2.0*dPn1(k)*Pn3(l) - Pn1(k)*dPn3(l));
          
          num_dof++;
	}
  }
  
    
  //! displays informations about class TriangleHcurlHierarchic
  ostream& operator <<(ostream& out, const TriangleHcurlHierarchic& e)
  {
    out<<static_cast<const TriangleReference<2>&>(e);
    out<<"Triangle of nedelec's second family "<<endl;
    out<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    
    return out;
  }

}
  
#define MONTJOIE_FILE_TRIANGLE_HCURL_HIERARCHIC_CXX
#endif
