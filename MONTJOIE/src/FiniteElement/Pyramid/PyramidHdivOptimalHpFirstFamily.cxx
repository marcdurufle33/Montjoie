#ifndef MONTJOIE_FILE_PYRAMID_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidHdivOptimalHpFirstFamily::PyramidHdivOptimalHpFirstFamily() : PyramidReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! how to number mesh
  void PyramidHdivOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<3>::ConstructNumberMap(nmap, dg);
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, (order+1)*(order+1));
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
    
    // dofs inside the Pyramid
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ElementReference<Dimension2, 3>::FindHdivRotationTri(order, this->Points2D_tri(),
							 this->Weights2D_tri(),
							 ValuePhiTri2D, nmap);
    
    ElementReference<Dimension2, 3>::FindHdivRotationQuad(order, NumQuad2D,
							  ValuePhiQuad2D, nmap);
    
  }


  size_t PyramidHdivOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = PyramidReference<3>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize() + JacobiPolynom11.GetMemorySize();
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol);
    taille += ValLeg.GetMemorySize() + ValuePhiQuad2D.GetMemorySize() + ValuePhiTri2D.GetMemorySize();
    return taille;
  }
  
  
  //! construction of finite element
  void PyramidHdivOptimalHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                         int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    
    if (rquad <= r)
      rquad = r+1;
    
    //type_quad = PyramidQuadrature::QUADRATURE_GAUSS;
    type_quad = PyramidQuadrature::QUADRATURE_JACOBI1;    
    PyramidReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad, r, r);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Div_Phi and elementary matrices
    ConstructHdivElementaryMatrix();
    
    //this->mass_matrix.Write("MhPyrDiv"+to_str(this->order)+".dat");
    //this->const_div_matrix.Write("KhPyrDiv"+to_str(this->order)+".dat");
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace();

    TriangleDgOrtho* Fb_tri = new TriangleDgOrtho();
    element_tri_surf = Fb_tri;
    Fb_tri->ConstructFiniteElement(order-1, order, order, TriangleQuadrature::QUADRATURE_GAUSS);

    QuadrangleDgOrtho* Fb_quad = new QuadrangleDgOrtho();
    element_quad_surf = Fb_quad;
    Fb_quad->ConstructFiniteElement(order, order, order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
  }
  
  
  //! construction of basis functions
  void PyramidHdivOptimalHpFirstFamily::ConstructFunctions()
  {  		
    nb_dof_tri = order*(order+1)/2;
    nb_dof_quad = (order+1)*(order+1);
    nb_dof_boundaries = 4*nb_dof_tri + nb_dof_quad;
    nb_dof_loc = (order+1)*(2*order*order+7*order+2)/2;    
    
    VectR3 points_dof3d;
    points_dof3d = this->PointsND();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDof2D_tri(this->Points2D_tri());
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(points_dof3d);
    
    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    GetJacobiPolynomial(JacobiPolynom11, order, Real_wp(1), Real_wp(1));
    jacobi_2ip1_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip1_pol(i), order, Real_wp(2*i+1), Real_wp(0));
    
    ValLeg.Reallocate(order+1, order+1);
    ValLeg.Fill(0);
    VectReal_wp Px;
    for (int i = 0; i <= order; i++)
      {
        EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*this->Points1D(i)-1.0, Px);
	for (int j = 0; j <= order; j++)
          ValLeg(j, i) = Px(j);
      }
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ValuePhiQuad2D.Reallocate(nb_dof_quad, this->Points2D_quad().GetM());
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k1 = 0; k1 <= order; k1++)
	  for (int k2 = 0; k2 <= order; k2++)
	    ValuePhiQuad2D(i*(order+1) + j, NumQuad2D(k1, k2)) = ValLeg(i, k1)*ValLeg(j, k2);
    
    ValuePhiTri2D.Reallocate(nb_dof_tri, this->Points2D_tri().GetM());
    for (int k = 0; k < this->Points2D_tri().GetM(); k++)
      {
	Real_wp x = this->Points2D_tri()(k)(0), y = this->Points2D_tri()(k)(1);
	Real_wp a = 2.0*x/(1.0-y) - 1.0;
	
	VectReal_wp Px, Py;
	EvaluateJacobiPolynomial(LegendrePolynom, order-1, a, Px);
		
	Real_wp pow_one_minus_y = 1.0;
	int nb = 0;
	for (int i = 0; i < order; i++)
	  {
	    EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, 2.0*y-1.0, Py);
	    for (int j = 0; j < order-i; j++)
	      {
		Real_wp val = Px(i) * pow_one_minus_y * Py(j);
		ValuePhiTri2D(nb, k) = val;
		nb++;
	      }
	    
	    pow_one_minus_y *= 1.0-y;
	  }
      }
    //DISP(ValuePhiTri2D);
  }
  
  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void PyramidHdivOptimalHpFirstFamily
  ::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    for (int i = 0; i < feval_weight.GetM(); i+=3)
      {
        int j = i/3;
        feval_weight(i) *= this->WeightsND(j);
        feval_weight(i+1) *= this->WeightsND(j);
        feval_weight(i+2) *= this->WeightsND(j);
      }
    
    contrib.Reallocate(nb_dof_loc);
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void PyramidHdivOptimalHpFirstFamily
  ::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
		
    // on recupere x, y, z de la pyramide symetrique
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);   		
    Real_wp omz, invOmz(0);
    omz = 1.0 - z;
    
    int num = 0;
    
    // coordinates on the cube
    Real_wp a(0), b(0), c;
    if (abs(omz) > epsilon_machine)
      {
	invOmz = 1.0/(1.0-z);
	a = x*invOmz;
	b = y*invOmz;
      }
    
    c = 2.0*z - 1.0;
    
    VectReal_wp powOneMinusZ(order+1);
    powOneMinusZ(0) = 1.0;
    for (int i = 0; i <= order-1; i++)
      powOneMinusZ(i+1) = powOneMinusZ(i)*omz;
	  
    // values of Legendre polynomials for edges
    VectReal_wp Px, Py, Px2, Py2;
    EvaluateJacobiPolynomial(LegendrePolynom, order, x, Px2);
    EvaluateJacobiPolynomial(LegendrePolynom, order, y, Py2);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order, b, Py);
    
    Vector<VectReal_wp> Pz(order);
    for (int i = 0; i <= order-1; i++)
      EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, c, Pz(i));
    
    R3 v1( (2.0*(1.0+x-z) - x*z)/(1.0-z), -y*z/(1.0-z), z);
    R3 v2( (-2.0*(1.0-x-z) - x*z)/(1.0-z), -y*z/(1.0-z), z);
    R3 v3(-x*z/(1.0-z), (2.0*(1.0+y-z) - y*z)/(1.0-z), z);
    R3 v4(-x*z/(1.0-z), (-2.0*(1.0-y-z) - y*z)/(1.0-z), z);
    R3 v5(-x, -y, (1.0-z));
    
    R3 v6(x*z*y/(1.0-z), (1.0-z-y*y), -y*z);
    R3 v7(1.0-z-x*x, x*z*y/(1.0-z), -x*z);
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  Real_wp val = 0.25*Px(i)*Py(j)*powOneMinusZ(max(i, j))*invOmz;
	  res(num) = val*v5;
	  num++;
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{
	  Real_wp val = 0.125*Px(i)*Pz(i)(j)*powOneMinusZ(i);
	  res(num++) = val*v4;
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{	
	  Real_wp val = 0.125*Py(i)*Pz(i)(j)*powOneMinusZ(i);
	  res(num++) = val*v1;
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{	
	  Real_wp val = 0.125*Px(i)*Pz(i)(j)*powOneMinusZ(i);
	  res(num++) = -val*v3;
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{	
	  Real_wp val = 0.125*Py(i)*Pz(i)(j)*powOneMinusZ(i);
	  res(num++) = -val*v2;
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
	for (int k = 0; k <= order -1 - max(i, j); k++)
	  {
	    int m = max(i, j);
	    Real_wp val = Px(i)*Py(j)*powOneMinusZ(m)*Pz(m)(k)*invOmz;
	    res(num) = val*z*v5;
	    res(num+1) = val*v6;
	    res(num+2) = val*v7;
	    
	    num += 3;
	  }
    
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void PyramidHdivOptimalHpFirstFamily
  ::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    res.Fill(0);

    // on recupere x, y, z de la pyramide symetrique
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    Real_wp omz, invOmz(0);
    omz = 1.0 - z;
    
    int num = 0;
    
    // coordinates on the cube
    Real_wp a(0), b(0), c, da_dx(0), db_dy(0), da_dz(0), db_dz(0);
    if (abs(omz) > epsilon_machine)
      {
	invOmz = 1.0/(1.0-z);
	a = x*invOmz;
	b = y*invOmz;
	da_dx = invOmz;
	db_dy = invOmz;
	da_dz = x*invOmz*invOmz;
	db_dz = y*invOmz*invOmz;
      }
    
    c = 2.0*z - 1.0;
    
    VectReal_wp powOneMinusZ(order+2), dpowOneMinusZ(order+1);
    powOneMinusZ(0) = 1.0;
    dpowOneMinusZ(0) = 0.0;
    for (int i = 0; i <= order; i++)
      {
	powOneMinusZ(i+1) = powOneMinusZ(i)*omz;
	if (i > 0)
	  dpowOneMinusZ(i) = -Real_wp(i)*powOneMinusZ(i-1);
      }
    
    // values of Legendre polynomials for edges
    VectReal_wp Px, Py, Px2, Py2, dPx, dPy, dPx2, dPy2;
    EvaluateJacobiPolynomial(LegendrePolynom, order, x, Px2, dPx2);
    EvaluateJacobiPolynomial(LegendrePolynom, order, y, Py2, dPy2);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order, b, Py, dPy);
    
    Vector<VectReal_wp> Pz(order), dPz(order);
    for (int i = 0; i <= order-1; i++)
      EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, c, Pz(i), dPz(i));
    
    // divergence of v1, v2, v3, v4 is equal to 3
    R3 v1( (2.0*(1.0+x-z) - x*z)/(1.0-z), -y*z/(1.0-z), z);
    R3 v2( (-2.0*(1.0-x-z) - x*z)/(1.0-z), -y*z/(1.0-z), z);
    R3 v3(-x*z/(1.0-z), (2.0*(1.0+y-z) - y*z)/(1.0-z), z);
    R3 v4(-x*z/(1.0-z), (-2.0*(1.0-y-z) - y*z)/(1.0-z), z);

    // divergence of v5 is equal to -3
    R3 v5(-x, -y, (1.0-z));
    
    R3 v6(x*z*y/(1.0-z), (1.0-z-y*y), -y*z);
    R3 v7(1.0-z-x*x, x*z*y/(1.0-z), -x*z);
    
    // divergence of v6 and v7
    Real_wp div_v6, div_v7;
    div_v6 = y*(z/(1.0-z) - 3.0);
    div_v7 = x*(z/(1.0-z) - 3.0);
    
    Real_wp val; R3 grad;
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  int m = max(i, j);
	  val = Px(i)*Py(j)*powOneMinusZ(m)*invOmz;
	  grad(0) = da_dx*dPx(i)*Py(j)*powOneMinusZ(m)*invOmz;
	  grad(1) = db_dy*Px(i)*dPy(j)*powOneMinusZ(m)*invOmz;
	  grad(2) = powOneMinusZ(m)*invOmz*(da_dz*dPx(i)*Py(j) + db_dz*Px(i)*dPy(j))
	    + Px(i)*Py(j)*(powOneMinusZ(m)*invOmz*invOmz + dpowOneMinusZ(m)*invOmz);
	  
	  res(num) = 0.25*(-3.0*val + DotProd(grad, v5));
	  num++;
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{
	  val = Px(i)*Pz(i)(j)*powOneMinusZ(i);
	  grad(0) = dPx(i)*da_dx*Pz(i)(j)*powOneMinusZ(i);
	  grad(1) = 0.0;
	  grad(2) = dPx(i)*da_dz*Pz(i)(j)*powOneMinusZ(i)
            + 2.0*dPz(i)(j)*Px(i)*powOneMinusZ(i) + Px(i)*Pz(i)(j)*dpowOneMinusZ(i);
          
	  res(num++) = 0.125*(3.0*val + DotProd(grad, v4));
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{	
	  val = Py(i)*Pz(i)(j)*powOneMinusZ(i);
	  grad(0) = 0.0;
	  grad(1) = dPy(i)*db_dy*Pz(i)(j)*powOneMinusZ(i);	  
	  grad(2) = dPy(i)*db_dz*Pz(i)(j)*powOneMinusZ(i) 
            + 2.0*Py(i)*dPz(i)(j)*powOneMinusZ(i) + Py(i)*Pz(i)(j)*dpowOneMinusZ(i);

	  res(num++) = 0.125*(3.0*val + DotProd(grad, v1));
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{	
	  val = Px(i)*Pz(i)(j)*powOneMinusZ(i);
	  grad(0) = dPx(i)*da_dx*Pz(i)(j)*powOneMinusZ(i);
	  grad(1) = 0.0;
	  grad(2) = dPx(i)*da_dz*Pz(i)(j)*powOneMinusZ(i)
            + 2.0*dPz(i)(j)*Px(i)*powOneMinusZ(i) + Px(i)*Pz(i)(j)*dpowOneMinusZ(i);
          
	  res(num++) = 0.125*(-3.0*val - DotProd(grad, v3));
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
	{	
          val = Py(i)*Pz(i)(j)*powOneMinusZ(i);
	  grad(0) = 0.0;
	  grad(1) = dPy(i)*db_dy*Pz(i)(j)*powOneMinusZ(i);	  
	  grad(2) = dPy(i)*db_dz*Pz(i)(j)*powOneMinusZ(i) 
            + 2.0*Py(i)*dPz(i)(j)*powOneMinusZ(i) + Py(i)*Pz(i)(j)*dpowOneMinusZ(i);
          
	  res(num++) = 0.125*(-3.0*val - DotProd(grad, v2));
	}
    
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1; j++)
	for (int k = 0; k <= order -1 - max(i, j); k++)
	  {
	    int m = max(i, j);
	    val = Px(i)*Py(j)*powOneMinusZ(m)*Pz(m)(k)*invOmz;
	    grad(0) = da_dx*dPx(i)*Py(j)*powOneMinusZ(m)*Pz(m)(k)*invOmz;
	    grad(1) = db_dy*Px(i)*dPy(j)*powOneMinusZ(m)*Pz(m)(k)*invOmz;
	    grad(2) = powOneMinusZ(m)*Pz(m)(k)*invOmz*(da_dz*dPx(i)*Py(j) + db_dz*Px(i)*dPy(j))
	      + Px(i)*Py(j)*invOmz*(powOneMinusZ(m)*(2.0*dPz(m)(k) + Pz(m)(k)*invOmz)
				    + dpowOneMinusZ(m)*Pz(m)(k));
	    
	    res(num) = -3.0*val*z + z*(grad(0)*v5(0) + grad(1)*v5(1) + grad(2)*v5(2))
	      + val*v5(2);
	    
	    res(num+1) = val*div_v6 + DotProd(grad, v6);
	    res(num+2) = val*div_v7 + DotProd(grad, v7);
	    num += 3;
	  }
  }
  

  void PyramidHdivOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    const VectReal_wp& inv_weight
      = static_cast<const TriangleDgOrtho&>(*this->element_tri_surf).GetInverseWeightFunction();
    
    switch (num_loc)
      {
      case 0:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = -contrib(i);
	
	break;
      case 1:
      case 2:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) *= inv_weight(i);    
	
	break;
      case 3:
      case 4:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) *= -inv_weight(i);    

	break;
      }
  }


  void PyramidHdivOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {    
    const VectReal_wp& inv_weight
      = static_cast<const TriangleDgOrtho&>(*this->element_tri_surf).GetInverseWeightFunction();
    
    switch (num_loc)
      {
      case 0:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = -contrib(i);
	
	break;
      case 1:
      case 2:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) *= inv_weight(i);    
	
	break;
      case 3:
      case 4:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) *= -inv_weight(i);    

	break;
      }
  }

   
  //! displays details of class PyramidHdivOptimalHpFirstFamily
  ostream& operator <<(ostream& out, const PyramidHdivOptimalHpFirstFamily& e)
  {
    out<<static_cast<const PyramidReference<3>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
