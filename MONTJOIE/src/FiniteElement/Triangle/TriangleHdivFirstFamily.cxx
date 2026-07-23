#ifndef MONTJOIE_FILE_TRIANGLE_HDIV_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TriangleHdivFirstFamily::TriangleHdivFirstFamily() : TriangleReference<3>()
  {
    // dof, nodal and quadrature points are all different
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    // nodal points used to express basis functions
    // NODAL_REGULAR -> regular points (Graglia expressions)
    // NODAL_GAUSS -> gauss points along edges (for conformity with quadrilaterals)
    type_nodal_basis = NODAL_GAUSS;
    //type_nodal_basis = NODAL_REGULAR;
  }
  
  
  //! how to number the mesh
  void TriangleHdivFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TriangleReference<3>::ConstructNumberMap(nmap, dg);

    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 0);
    nmap.SetNbDofTriangle(order, nb_dof_loc - 3*order);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t TriangleHdivFirstFamily::GetMemorySize() const
  {
    size_t taille = TriangleReference<3>::GetMemorySize();
    taille += InverseVDM_Edge.GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    return taille;
  }
  
    
  //! construction of finite element
  void TriangleHdivFirstFamily::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			 int rsurf, int type_surf)
  {
    TriangleReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    ConstructFunctions();    
    
    this->Fb_geom.ComputeCoefficientTransformation();
    ConstructHdivElementaryMatrix();    
  }
  
  
  //! construction of basis functions
  void TriangleHdivFirstFamily::ConstructFunctions()
  {
    // dofs on the three edges and inside dofs (two by inside nodes)
    nb_dof_loc = order*3 + (order-1)*order;
    nb_dof_boundaries = 3*order;
    
    EdgesDof.Reallocate(order, 3);
    
    // First Edge
    for (int i = 0; i < order; i++)
      EdgesDof(i, 0) = i;
          
    // Second Edge
    for (int i = 0; i < order; i++)
      EdgesDof(i, 1) = i + order;
          
    // Third Edge
    for (int i = 0; i < order; i++)
      EdgesDof(i, 2) = i + 2*order;	      

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);              
    nb_dof_boundaries = 3*order;
    VectReal_wp points_lob;
    VectReal_wp points_dof1d; VectR2 points_dof2d;
    points_dof1d.Reallocate(order);    
    
    // location of dofs on the edges
    VectReal_wp poids;
    ComputeGaussLegendre(points_dof1d, poids, order-1);
    if (type_nodal_basis == NODAL_LOBATTO)
      {
	ComputeGaussLobatto(points_lob, poids, order+1);
	points_dof1d.Reallocate(order);
	for (int i = 0; i < order; i++)
	  points_dof1d(i) = points_lob(i+1);
      }
    else if (type_nodal_basis == NODAL_REGULAR)
      {
	points_dof1d.Reallocate(order);
	for (int i = 0; i < order; i++)
	  points_dof1d(i) = Real_wp(i+1)/(order+1);
      }
    
    // weight for orthonormalization
    InvWeightBasisRr.Reallocate(nb_dof_loc);
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefficientLegendre();
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefficientOddJacobi();	
    int num = 0;
    for (int i = 0; i < order; i++)
      {
	for (int j = 0; j <= order-1-i; j++)
	  {
	    Real_wp coef = CoefLegendre(i)*CoefOddJacobi(i, j);
	    InvWeightBasisRr(num) = coef;
	    InvWeightBasisRr(num+1) = coef;
	    num += 2;
	  }
	
	InvWeightBasisRr(num) = CoefLegendre(i)*CoefOddJacobi(i, order-1-i);
	num++;
      }
        
    InvWeightBasisRr.Fill(1.0);
    VectReal_wp Scale(nb_dof_loc);
    VectR2 phi(nb_dof_loc); Scale.Fill(0);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
	ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
	for (int i = 0; i < nb_dof_loc; i++)
	  Scale(i) += DotProd(phi(i), phi(i))*this->WeightsND(k);
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      InvWeightBasisRr(i) = 1.0/sqrt(Scale(i));
    
    // points where dofs are defined, and tangent on these points
    VectR2 tangente_dof(nb_dof_loc), points_inside;
    Fb_geom.ConstructLobattoPoints(order+1, 0, points_lob, points_inside);
    points_dof2d.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      {
	points_dof2d(i).Init(points_dof1d(i), 0);
	tangente_dof(i).Init(0, -1.0);
            
	points_dof2d(order+i).Init(1.0-points_dof1d(i), points_dof1d(i));
	tangente_dof(order+i).Init(1.0, 1.0);
	
	points_dof2d(2*order+i).Init(0, 1.0-points_dof1d(i));
	tangente_dof(2*order+i).Init(-1.0, 0.0);
      }
    
    for (int i = 0; i < (order-1)*order/2; i++)
      {
	points_dof2d(3*order + 2*i) = points_inside(3*(order+1) + i);
	points_dof2d(3*order + 2*i+1) = points_inside(3*(order+1) + i);
	
	tangente_dof(3*order + 2*i).Init(1.0, 0);
	tangente_dof(3*order + 2*i+1).Init(0, 1.0);
      }
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    /* ofstream file_out("points.dat"); file_out.precision(15);
    for (int i = 0; i < nb_dof_loc; i++)
      {
        Real_wp x = points_dof2d(i)(0), y = points_dof2d(i)(1);
        Real_wp tx = tangente_dof(i)(0), ty = tangente_dof(i)(1);
        file_out << "\\draw[color=blue,line width=0.05cm,->] " << "(" << 10*x - 0.5*tx << ", " << 10*y - 0.5*ty << ") -- (" << 10*x + 0.5*tx  << ", " << 10*y + 0.5*ty << ");" << endl;
      }

      file_out.close(); */
    
    // computation of VDM = psi_j(xi_k) \cdot t_k
    InverseVDM_Edge.Reallocate(nb_dof_loc, nb_dof_loc);
    for (int k = 0; k < nb_dof_loc; k++)
      {
	ComputeValuesPhiOrthoRef(points_dof2d(k), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  InverseVDM_Edge(j, k) = DotProd(phi(j), tangente_dof(k));
      }
    
    //InverseVDM_Edge.WriteText("VdmTri"+to_str(order)+".dat");
    GetInverse(InverseVDM_Edge);
    
    // then checking phi_i(xi_k) \cdot t_k = delta_{i, k}
    for (int k = 0; k < nb_dof_loc; k++)
      {
	ComputeValuesPhiRef(points_dof2d(k), phi);            
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    Real_wp val = DotProd(phi(j), tangente_dof(k));
	    Real_wp val_exact = 0;
	    if (j == k)
	      val_exact = 1.0;
	    
	    if (abs(val - val_exact) > 1e4*epsilon_machine)
	      {
		DISP(val); DISP(val_exact);
		abort();
	      }
	  }
      }

    this->nb_points_dof_inside = nb_dof_loc;
    this->num_dof_points_surf.Reallocate(3);
    for (int n = 0; n < 3; n++)
      {
	this->num_dof_points_surf(n).Reallocate(order);
	for (int i = 0; i < order; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }    
  }
  
  
  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void TriangleHdivFirstFamily::
  ComputeValuesPhiOrthoRef(const R2& pointloc, VectR2& phi) const
  {
    Real_wp xtilde(0), ytilde(1);    
    if (abs(pointloc(1) - 1.0) > epsilon_machine)
      {
	xtilde = 2.0*pointloc(0)/(1.0-pointloc(1)) - 1.0;
	ytilde = 2.0*pointloc(1) - 1.0;
      }
    
    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    
    int num = 0;
    VectReal_wp Li, Pj; phi.Reallocate(nb_dof_loc);
    EvaluateJacobiPolynomial(LegendrePolynom, order, xtilde, Li);
    Real_wp pow_oneMinusytilde = 1.0, oneMinusy = 0.5*(1.0-ytilde), vali;
    for (int i = 0; i < order; i++)
      {
        if (i == 0)
          vali = 1.0;
        else if (i == 1)
          vali = 2.0*pointloc(0) + pointloc(1) -1.0;
        else
          vali = pow_oneMinusytilde*Li(i);
        
	EvaluateJacobiPolynomial(OddJacobiPolynom(i), order-i, ytilde, Pj);
        // part due to P_{r-1}^2
	for (int j = 0; j <= order-1-i; j++)
	  {
            phi(num).Init(Pj(j)*vali, 0);
            phi(num+1).Init(0, Pj(j)*vali);
            num += 2;
          }
        
        // part due to D_r
        phi(num++).Init(vali*Pj(order-1-i)*pointloc(0), vali*Pj(order-1-i)*pointloc(1));
        
	pow_oneMinusytilde *= oneMinusy;
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      phi(i) *= InvWeightBasisRr(i);
  }
  
  
  //! Evaluating divergence nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] div_phi divergence of nearly orthogonal functions on pointloc
  */
  void TriangleHdivFirstFamily::
  ComputeDivPhiOrthoRef(const R2& pointloc, VectReal_wp& div_phi) const
  {
    Real_wp xtilde(0), ytilde(1), dxtilde_dx(0), dxtilde_dy(0);
    if (abs(pointloc(1) - 1.0) > epsilon_machine)
      {
	xtilde = 2.0*pointloc(0)/(1.0-pointloc(1)) - 1.0;
        dxtilde_dx = 2.0/(1.0-pointloc(1));
        dxtilde_dy = dxtilde_dx*pointloc(0)/(1.0-pointloc(1));
        
        ytilde = 2.0*pointloc(1) - 1.0;
      }

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    
    int num = 0;
    VectReal_wp Li, Pj, dLi, dPj;
    div_phi.Reallocate(nb_dof_loc);
    EvaluateJacobiPolynomial(LegendrePolynom, order, xtilde, Li, dLi);
    Real_wp powm1_oneMinusytilde = 0.0, pow_oneMinusytilde = 1.0,
      oneMinusy = 0.5*(1.0-ytilde), vali;
    
    R2 gradi;
    for (int i = 0; i < order; i++)
      {
        if (i == 0)
          {
            vali = 1.0;
            gradi.Fill(0);
          }
        else if (i == 1)
          {
            vali = 2.0*pointloc(0) + pointloc(1) -1.0;
            gradi.Init(2.0, 1.0);
          }
        else
          {
            vali = pow_oneMinusytilde*Li(i);
            gradi(0) = pow_oneMinusytilde*dLi(i)*dxtilde_dx;
            gradi(1) = -i*powm1_oneMinusytilde*Li(i) + pow_oneMinusytilde*dLi(i)*dxtilde_dy;
          }
        
	EvaluateJacobiPolynomial(OddJacobiPolynom(i), order-i, ytilde, Pj, dPj);
        // part due to P_{r-1}^2
	for (int j = 0; j <= order-1-i; j++)
	  {
            div_phi(num) = Pj(j)*gradi(0);
            div_phi(num+1) = (2.0*vali*dPj(j) + gradi(1)*Pj(j));
            num += 2;
          }
        
        // part due to S_r
        div_phi(num) = vali*Pj(order-1-i) + pointloc(0)*Pj(order-1-i)*gradi(0);
        div_phi(num) += vali*Pj(order-1-i) + pointloc(1)*(2.0*vali*dPj(order-1-i)
                                                          + gradi(1)*Pj(order-1-i));
        num++;
                
        powm1_oneMinusytilde = pow_oneMinusytilde;
	pow_oneMinusytilde *= oneMinusy;
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      div_phi(i) *= InvWeightBasisRr(i);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TriangleHdivFirstFamily::
  ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectR2 psi;
    ComputeValuesPhiOrthoRef(point_loc, psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	Add(InverseVDM_Edge(i, j), psi(j), res(i));
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void TriangleHdivFirstFamily::ComputeDivPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    VectReal_wp div_psi;
    ComputeDivPhiOrthoRef(point_loc, div_psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	res(i) += InverseVDM_Edge(i, j)*div_psi(j);
  }
  
  
  //! displays informations about class TriangleHdivFirstFamily
  ostream& operator <<(ostream& out, const TriangleHdivFirstFamily& e)
  {
    out<<static_cast<const TriangleReference<3>& >(e);
    out<<"Triangle of Nedelec's first family "<<endl;
    out<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    return out;
  }


  //! default constructor
  TriangleHdivOptimalFirstFamily::TriangleHdivOptimalFirstFamily() : TriangleHdivFirstFamily()
  {
    type_nodal_basis = NODAL_LOBATTO;
  }

  //! displays informations about class TriangleHcurlOptimalFirstFamily
  ostream& operator <<(ostream& out, const TriangleHdivOptimalFirstFamily& e)
  {
    out<<static_cast<const TriangleReference<3>& >(e);
    out<<"Triangle of Nedelec's first family "<<endl;
    out<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    return out;
  }

}
  
#define MONTJOIE_FILE_TRIANGLE_HDIV_FIRST_FAMILY_CXX
#endif
