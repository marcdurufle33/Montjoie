#ifndef MONTJOIE_FILE_TRIANGLE_HIERARCHIC_CXX

namespace Montjoie
{
  
  //! default constructor
  TriangleHierarchic::TriangleHierarchic() : TriangleReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    //type_basis = TENSOR_BASIS;
    type_basis = INVARIANT_BASIS;
  }
  

  size_t TriangleHierarchic::GetMemorySize() const
  {
    size_t taille = TriangleReference<1>::GetMemorySize();
    taille += jacobi_11_pol.GetMemorySize();
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol);
    taille += CoefLeg11.GetMemorySize();
    taille += CoefJacobi.GetMemorySize();
    taille += NumDofsTri.GetMemorySize();
    return taille;
  }

  
  //! method useful to number mesh
  void TriangleHierarchic::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TriangleReference<1>::ConstructNumberMap(nmap, dg);

    if (dg == ElementReference_Base::HDG)
      {
	nmap.SetFormulationDG(ElementReference_Base::HDG);
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, this->order+1);
	nmap.SetNbDofTriangle(this->order, 0);

	nmap.SetEqualEdgesDofSymmetry(this->order, this->order+1);
	nmap.SetEvenEdgesDofToSkewSymmetric(this->order);

	nmap.SetEdgesDofSymmetry(this->order, 0, this->order, false);
	nmap.SetEdgesDofSymmetry(this->order, this->order, 0, false);
	
	return;
      }
    
    nmap.SetFormulationDG(ElementReference_Base::CONTINUOUS);
    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
    nmap.SetNbDofQuadrangle(order, (order-1)*(order-1));
    nmap.SetNbDofTriangle(order, (order-1)*(order-2)/2);
    
    nmap.SetEqualEdgesDofSymmetry(order, order-1);
    nmap.SetOddEdgesDofToSkewSymmetric(order);
  }

  
  //! constructing finite element
  void TriangleHierarchic::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						  int rsurf, int type_surf)
  {
    TriangleReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructElementaryMatrix(*this);
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    delete this->element_surface;
    EdgeHierarchicReference* edge = new EdgeHierarchicReference();
    edge->ConstructFiniteElement(r, r, rquad, type_surf);
    this->element_surface = edge;    
  }
  
  
  //! constructions of basis functions
  void TriangleHierarchic::ConstructFunctions()
  {
    nb_dof_boundaries = 3*order;
    
    EdgesDof.Reallocate(order+1, 3);
    // changing EdgesDof
    EdgesDof(0, 0) = 0; EdgesDof(order, 0) = 1;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 0) = 2+i;
    
    EdgesDof(0, 1) = 1; EdgesDof(order, 1) = 2;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 1) = 1+i+order;

    EdgesDof(0, 2) = 2; EdgesDof(order, 2) = 0;
    for (int i = 1; i < order; i++)
      EdgesDof(i, 2) = i+2*order;
    
    nb_dof_loc = (order+1)*(order+2)/2;
    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);    
    
    VectR2 points_dof2d; VectReal_wp points_dof1d;
    points_dof2d = this->PointsND();
    this->nb_points_dof_inside = nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    this->elt_geom.dof_equal_quadrature = true;
    
    points_dof1d = this->Points1D();
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);
    
    // computing P_m^{1,1}
    GetJacobiPolynomial(jacobi_11_pol, order, Real_wp(1), Real_wp(1));
    jacobi_2ip1_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip1_pol(i), order, Real_wp(2*(i+1)+1), Real_wp(1));
    
    CoefLeg11.Reallocate(order-1); CoefLeg11.Fill(0);
    CoefJacobi.Reallocate(order-1, order-1); CoefJacobi.Fill(0);
    VectReal_wp Pn;
    for (int i = 0; i <= order_quad; i++)
      {
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, 2*this->Points1D(i) - 1.0, Pn);
        for (int j = 1; j < order; j++)
          CoefLeg11(j-1) += this->Weights1D(i)*square(this->Points1D(i)*(1.0-this->Points1D(i))*Pn(j-1));

        for (int k = 0; k < order-1; k++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(k), order-2, 2*this->Points1D(i) - 1.0, Pn);
            for (int j = 1; j < order; j++)
              CoefJacobi(k, j-1) += this->Weights1D(i)
                *square(this->Points1D(i)*pow(1.0-this->Points1D(i), k+1)*Pn(j-1));
          }
      }
    
    for (int j = 1; j < order; j++)
      {
        CoefLeg11(j-1) = 1.0/sqrt(CoefLeg11(j-1));
        for (int k = 0; k < order - 1; k++)
          if (abs(CoefJacobi(k, j-1)) > 0)
            CoefJacobi(k, j-1) = 1.0/sqrt(CoefJacobi(k, j-1));
      }
    
    // hierarchic numbering of internal nodes
    NumDofsTri.Reallocate(order+1, order+1);
    NumDofsTri.Fill(-1);
    NumDofsTri(0, 0) = 0;
    NumDofsTri(order, 0) = 1;
    NumDofsTri(0, order) = 2;
    int nb = 3;
    for (int i = 1; i < order; i++)
      NumDofsTri(i, 0) = nb++;
    
    for (int i = 1; i < order; i++)
      NumDofsTri(order-i, i) = nb++;

    for (int i = 1; i < order; i++)
      NumDofsTri(0, i) = nb++;
    
    for (int diag = 1; diag < order; diag++)
      {
        for (int i = 1; i < diag; i++)
          NumDofsTri(diag-i, i) = nb++;
      }
  }
  
  
  //! computation of mass matrix
  void TriangleHierarchic::ConstructMassMatrix()
  {    
  }
  
  
  //! computation of stiffness matrix R_h
  void TriangleHierarchic::ConstructStiffnessMatrix()
  {
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point_loc
  */
  void TriangleHierarchic::ComputeValuesPhiRef(const R2& point_loc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    Real_wp L0 = 1.0-point_loc(0)-point_loc(1);
    Real_wp L1 = point_loc(0);
    Real_wp L2 = point_loc(1);
    phi(0) = L0;
    phi(1) = L1;
    phi(2) = L2;
    
    if (order <= 1)
      return;
    
    if (type_basis == TENSOR_BASIS)
      {
        // computing a and b
        Real_wp a(0);
        if (point_loc(1) != Real_wp(1))
          a = 2*point_loc(0)/(1-point_loc(1)) - 1.0;
        
        Real_wp b = 2*point_loc(1) - 1.0;
        VectReal_wp Pn1, Pn2;
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, a, Pn1);
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, b, Pn2);
        for (int k = 0; k < order-1; k++)
          {
            Pn1(k) *= CoefLeg11(k);
            Pn2(k) *= CoefLeg11(k);
          }
        
        // first edge
        int num_dof = 3;
        for (int i = 1; i < order; i++)
          phi(num_dof++) = 0.25*(1-a)*(1+a)*Pn1(i-1)*pow(0.5*(1-b), i+1);
        
        // second edge
        for (int i = 1; i < order; i++)
          phi(num_dof++) = 0.125*(1+a)*(1-b)*(1+b)*Pn2(i-1);
        
        // third edge
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, -b, Pn2);
        for (int k = 0; k < order-1; k++)
          Pn2(k) *= CoefLeg11(k);
        
        for (int i = 1; i < order; i++)
          phi(num_dof++) = 0.125*(1-a)*(1-b)*(1+b)*Pn2(i-1);
        
        Real_wp valx(1), valy(1);
        for (int i = 1; i < order; i++)
          {
            valx = 0.25*(1-a)*(1+a)*Pn1(i-1)*pow(0.5*(1-b), i+1);
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-i-2, b, Pn2);
            for (int k = 0; k < order-i-1; k++)
              Pn2(k) *= CoefJacobi(i-1, k);
 
            for (int j = 1; j < order-i; j++)
              {
                valy = 0.5*(1+b)*Pn2(j-1);
                phi(NumDofsTri(i, j)) = valx*valy;
              }
          }
      }
    else
      {
        VectReal_wp Pn1, Pn2, Pn3;
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L1-L0, Pn1);
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L1, Pn2);
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L0-L2, Pn3);
        for (int k = 0; k < order-1; k++)
          {
            Pn1(k) *= CoefLeg11(k);
            Pn2(k) *= CoefLeg11(k);            
            Pn3(k) *= CoefLeg11(k);
          }
        
        int num_dof = 3;
        for (int i = 1; i < order; i++)
          phi(num_dof++) = L0*L1*Pn1(i-1);
        
        for (int i = 1; i < order; i++)
          phi(num_dof++) = L1*L2*Pn2(i-1);
        
        for (int i = 1; i < order; i++)
          phi(num_dof++) = L0*L2*Pn3(i-1);
        
        Real_wp gamma = 0, b = 2.0*point_loc(1) - 1.0;
        Real_wp pow_OneMinusY = 1.0;
        if (abs(1.0-point_loc(1)) > epsilon_machine)
          gamma = 2.0*point_loc(0)/(1.0-point_loc(1)) - 1.0;
        
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma, Pn1);
        for (int k = 0; k < order-1; k++)
          Pn1(k) *= CoefLeg11(k);
                
        Real_wp vali = 1.0;
        for (int i = 1; i < order; i++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, b, Pn2);
            for (int k = 0; k < order-i-1; k++)
              Pn2(k) *= CoefJacobi(i-1, k);
            
            if (i == 1)
              vali = CoefLeg11(i-1);
            else if (i == 2)
              vali = (L1 - L0)*CoefLeg11(i-1);
            else
              vali = pow_OneMinusY*Pn1(i-1);
            
            for (int j = 1; j < order-i; j++)
              {
                phi(NumDofsTri(i, j)) = L0*L1*L2*vali*Pn2(j-1);
              }

            pow_OneMinusY *= 1.0-point_loc(1);
          }
        
      }
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res gradient of basis functions on point_loc
  */
  void TriangleHierarchic::ComputeGradientPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc); 
    res(0).Init(-1.0, -1.0);
    res(1).Init(1.0, 0.0);
    res(2).Init(0.0, 1.0);
    
    if (order <= 1)
      return;
    
    if (type_basis == TENSOR_BASIS)
      {
        // computing a and b
        Real_wp a(0), da_dx(0), da_dy(0);
        if (point_loc(1) != Real_wp(1))
          {
            a = 2.0*point_loc(0)/(1-point_loc(1)) - 1.0;
            da_dx = 2.0/(1-point_loc(1));
            da_dy = 2.0*point_loc(0)/((1-point_loc(1))*(1-point_loc(1)));
          }
        
        Real_wp b = 2.0*point_loc(1) - 1.0;
        
        VectReal_wp Pn1, dPn1, Pn2, dPn2; R2 grad;
        // first edge
        int num_dof = 3;
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, a, Pn1, dPn1);
        for (int k = 0; k < order-1; k++)
          {
            Pn1(k) *= CoefLeg11(k);
            dPn1(k) *= CoefLeg11(k);
          }
        
        for (int i = 1; i < order; i++)
          {
            grad(0) = (-0.5*a*Pn1(i-1) + 0.25*(1-a)*(1+a)*dPn1(i-1))*pow(0.5*(1-b), i+1);
            grad(1) = -0.125*(1-a)*(1+a)*Pn1(i-1)*(i+1)*pow(0.5*(1-b), i); 
            res(num_dof)(0) = da_dx*grad(0);
            res(num_dof)(1) = da_dy*grad(0) + 2.0*grad(1); num_dof++;
          }
        
        // second edge
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, b, Pn2, dPn2);
        for (int k = 0; k < order-1; k++)
          {
            Pn2(k) *= CoefLeg11(k);
            dPn2(k) *= CoefLeg11(k);
          }
 
        for (int i = 1; i < order; i++)
          {
            res(num_dof)(0) = point_loc(1)*Pn2(i-1);
            res(num_dof)(1) = point_loc(0)*Pn2(i-1) + 2.0*point_loc(0)*point_loc(1)*dPn2(i-1);
            num_dof++;
          }
        
        // third edge
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, -b, Pn2, dPn2);
        for (int k = 0; k < order-1; k++)
          {
            Pn2(k) *= CoefLeg11(k);
            dPn2(k) *= CoefLeg11(k);
          }
 
        for (int i = 1; i < order; i++)
          {
            res(num_dof)(0) = -point_loc(1)*Pn2(i-1);
            res(num_dof)(1) = (1.0-point_loc(0)-2.0*point_loc(1))*Pn2(i-1)
              - 2.0*(1.0-point_loc(0)-point_loc(1))*point_loc(1)*dPn2(i-1);
            num_dof++;
          }
        
        Real_wp val1, val2, dval1, dval2, val3, dval3;
        for (int i = 1; i < order; i++)
          {
            val1 = 0.25*(1-a)*(1+a)*Pn1(i-1);
            val2 = pow(0.5*(1-b), i+1);
            dval1 = -0.5*a*Pn1(i-1) + 0.25*(1-a)*(1+a)*dPn1(i-1);
            dval2 = -0.5*(i+1)*pow(0.5*(1-b), i);
            
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-i-2, b, Pn2, dPn2);
            for (int k = 0; k < order-1-i; k++)
              {
                Pn2(k) *= CoefJacobi(i-1, k);
                dPn2(k) *= CoefJacobi(i-1, k);
              }
 
            for (int j = 1; j < order-i; j++)
              {
                val3 = 0.5*(1+b)*Pn2(j-1);
                dval3 = 0.5*Pn2(j-1) + 0.5*(1+b)*dPn2(j-1);
                
                grad(0) = dval1*val2*val3;
                grad(1) = val1*(dval2*val3 + val2*dval3);
                num_dof = NumDofsTri(i, j);
                res(num_dof)(0) = da_dx*grad(0);
                res(num_dof)(1) = da_dy*grad(0) + 2.0*grad(1);
              }
          }
      }
    else
      {
        Real_wp L0 = 1.0-point_loc(0)-point_loc(1);
        Real_wp L1 = point_loc(0);
        Real_wp L2 = point_loc(1);
    
        VectReal_wp Pn1, Pn2, Pn3, Pn3b, dPn1, dPn2, dPn3;
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L1-L0, Pn1, dPn1);
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L2-L1, Pn2, dPn2);
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, L0-L2, Pn3, dPn3);
        for (int k = 0; k < order-1; k++)
          {
            Pn1(k) *= CoefLeg11(k);
            Pn2(k) *= CoefLeg11(k);
            Pn3(k) *= CoefLeg11(k);
            dPn1(k) *= CoefLeg11(k);
            dPn2(k) *= CoefLeg11(k);
            dPn3(k) *= CoefLeg11(k);
          }
        
        int num_dof = 3;
        for (int i = 1; i < order; i++)
          res(num_dof++).Init((L0-L1)*Pn1(i-1) + 2.0*L0*L1*dPn1(i-1),
                              -L1*Pn1(i-1) + L0*L1*dPn1(i-1));
        
        for (int i = 1; i < order; i++)
          res(num_dof++).Init(L2*Pn2(i-1) - L1*L2*dPn2(i-1), L1*Pn2(i-1) + L1*L2*dPn2(i-1));
        
        for (int i = 1; i < order; i++)
          res(num_dof++).Init(-L2*Pn3(i-1) - L0*L2*dPn3(i-1),
                              (L0-L2)*Pn3(i-1) - 2.0*L0*L2*dPn3(i-1));
        
        Real_wp vloc = L0*L1*L2;
        Real_wp dv_dx = (L0-L1)*L2;
        Real_wp dv_dy = (L0-L2)*L1;

        Real_wp gamma = 0, dgamma_dx = 0, dgamma_dy = 0;
        Real_wp b = 2.0*point_loc(1) - 1.0;
        Real_wp pow_OneMinusY = 1.0, powM1_OneMinusY = 0;
        if (abs(1.0-point_loc(1)) > epsilon_machine)
          {
            gamma = 2.0*point_loc(0)/(1.0-point_loc(1)) - 1.0;
            dgamma_dx = 2.0/(1.0-point_loc(1));
            dgamma_dy = 2.0*point_loc(0)/square(1.0-point_loc(1));
          }
        
        EvaluateJacobiPolynomial(jacobi_11_pol, order-2, gamma, Pn1, dPn1);
        for (int k = 0; k < order-1; k++)
          {
            Pn1(k) *= CoefLeg11(k);
            dPn1(k) *= CoefLeg11(k);
          }
        
        Real_wp vali = 1.0; R2 gradi;
        for (int i = 1; i < order; i++)
          {
            EvaluateJacobiPolynomial(jacobi_2ip1_pol(i-1), order-2-i, b, Pn2, dPn2);
            for (int k = 0; k < order-1-i; k++)
              {
                Pn2(k) *= CoefJacobi(i-1, k);
                dPn2(k) *= CoefJacobi(i-1, k);
              }
            
            if (i == 1)
              {
                vali = CoefLeg11(i-1);
                gradi.Fill(0);
              }
            else if (i == 2)
              {
                vali = (L1 - L0)*CoefLeg11(i-1);
                gradi.Init(2.0, 1.0);
                gradi *= CoefLeg11(i-1);
              }
            else
              {
                vali = pow_OneMinusY*Pn1(i-1);
                gradi(0) = dgamma_dx*dPn1(i-1)*pow_OneMinusY;
                gradi(1) = dgamma_dy*dPn1(i-1)*pow_OneMinusY
                  - Real_wp(i-1)*powM1_OneMinusY*Pn1(i-1);
              }
            
            for (int j = 1; j < order-i; j++)
              {
                Real_wp prod = vali*Pn2(j-1);
                Real_wp dprod_dx = gradi(0)*Pn2(j-1);
                Real_wp dprod_dy = gradi(1)*Pn2(j-1) + 2.0*vali*dPn2(j-1);
                num_dof = NumDofsTri(i, j);
                res(num_dof).Init(vloc*dprod_dx + dv_dx*prod, vloc*dprod_dy + dv_dy*prod);
              }
            
            powM1_OneMinusY = pow_OneMinusY;
            pow_OneMinusY *= 1.0-point_loc(1);
          }
      }
  }
  
  
  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] res dof components
   */
  template<class Vector1>
  void TriangleHierarchic::ComputeProjectionDofGen(const Vector1& feval, Vector1& contrib) const
  {
    Vector1 feval_weight = feval;
    for (int i = 0; i < this->nb_points_dof_inside; i++)
      feval_weight(i) *= this->WeightsND(i);
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  

  //! computation of prolongation operator inside an element
  /*!
    \param[in,out] proj prolongation operator
    \param[in,out] LocalProlongation prolongation operator
    \param[in] FaceCoarse coarse finite element
    \param[in] FaceFine fine finite element
   */
  void TriangleHierarchic::
  ComputeLocalProlongation(FiniteElementProjector& proj, Matrix<Real_wp>& LocalProlongation,
			   const ElementReference_Dim<Dimension2>& FaceCoarse,
			   const ElementReference_Dim<Dimension2>& FaceFine) const
  {
    int rc = FaceCoarse.GetOrder();
    int rf = FaceFine.GetOrder();
    IVect num_coarse(FaceCoarse.GetNbDof()), num_fine(FaceCoarse.GetNbDof());
    num_coarse.Fill();
    num_fine.Fill(-1);
    // vertices
    num_fine(0) = 0; num_fine(1) = 1;
    num_fine(2) = 2;
    
    // edges
    int node = 3;
    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 3+i;

    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 3+(rf-1)+i;

    for (int i = 0; i < rc-1; i++)
      num_fine(node++) = 3+2*(rf-1)+i;

    // interior
    const TriangleHierarchic& Fb_c = static_cast<const TriangleHierarchic&>(FaceCoarse);
    const TriangleHierarchic& Fb_f = static_cast<const TriangleHierarchic&>(FaceFine);
    for (int i = 1; i < rf; i++)
      for (int j = 1; j < rf-i; j++)
        {
          if ((i+j) < rc)
            num_fine(Fb_c.NumDofsTri(i, j)) = Fb_f.NumDofsTri(i, j); 
        }
    
    proj.SetIdentity(num_coarse, num_fine, FaceFine.GetNbDof());
  }


  //! displays informations about class TriangleReference
  ostream& operator <<(ostream& out, const TriangleHierarchic& e)
  {
    out<<static_cast<const TriangleReference<1>&>(e);
    return out;
  }
  
}
  
#define MONTJOIE_FILE_TRIANGLE_HIERARCHIC_CXX
#endif
