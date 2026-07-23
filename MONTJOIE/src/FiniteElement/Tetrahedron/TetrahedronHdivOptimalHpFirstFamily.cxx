#ifndef MONTJOIE_FILE_TETRAHEDRON_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronHdivOptimalHpFirstFamily::TetrahedronHdivOptimalHpFirstFamily()
    : TetrahedronReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! how to number mesh
  void TetrahedronHdivOptimalHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<3>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
    
    // dofs inside the tetrahedron
    nmap.SetNbDofTetrahedron(order, nb_dof_loc-nb_dof_boundaries);
    nmap.SetNbDofHexahedron(order, 0);    

    ElementReference<Dimension2, 3>::FindHdivRotationTri(order, this->Points2D_tri(),
							 this->Weights2D_tri(),
							 ValuePhi2D, nmap);
    
    //this->FindHdivLinearCombinationRotation(*this, nmap, points2d_tri, points2d_quad,
    //					    ValuePhi2D, ValuePhi2D);
  }


  size_t TetrahedronHdivOptimalHpFirstFamily::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<3>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += NumDofsX.GetMemorySize() + NumDofsY.GetMemorySize() + NumDofsZ.GetMemorySize();
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol);
    taille += ValuePhi2D.GetMemorySize();
    return taille;
  }
  
    
  //! construction of finite element
  void TetrahedronHdivOptimalHpFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Div_Phi and elementary matrices
    ConstructHdivElementaryMatrix();

    //this->mass_matrix.Write("MhTetDiv"+to_str(this->order)+".dat");
    //this->const_div_matrix.Write("KhTetDiv"+to_str(this->order)+".dat");
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace();

    TriangleDgOrtho* Fb_tri = new TriangleDgOrtho();
    element_tri_surf = Fb_tri;
    Fb_tri->ConstructFiniteElement(order-1, order, order, TriangleQuadrature::QUADRATURE_GAUSS);
  }
  
  
  //! construction of basis functions
  void TetrahedronHdivOptimalHpFirstFamily::ConstructFunctions()
  {
    int nb_dof_inside_tri = order*(order+1)/2;	  
    nb_dof_boundaries = 4*nb_dof_inside_tri;
    nb_dof_loc = order*(order+1)*(order+3)/2;
    nb_dof_tri = nb_dof_inside_tri;
    nb_dof_quad = 0;

    VectR3 points_dof3d;    	  
    points_dof3d = this->PointsND();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDof2D_tri(this->Points2D_tri());
    this->SetPointsDofND(points_dof3d);
    
    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));

    jacobi_2ip1_pol.Reallocate(order);
    for (int i = 0; i < order; i++)
      GetJacobiPolynomial(jacobi_2ip1_pol(i), order, Real_wp(2*i+1), Real_wp(0));
    
    if (order >= 2)
      {
	NumDofsX.Reallocate(order-1, order-1, order-1);
	NumDofsY.Reallocate(order-1, order-1, order-1);
	NumDofsZ.Reallocate(order-1, order-1, order-1);
	NumDofsX.Fill(-1); NumDofsY.Fill(-1); NumDofsZ.Fill(-1);
	int offset = 3*order*(order+1);
	for (int diag = 0; diag <= order-2; diag++)
	  for (int i = 0; i <= diag; i++)
	    for (int j = 0; j <= diag-i; j++)
	      for (int k = 0; k <= diag-i-j; k++)
		if (i+j+k == diag)
		  {
		    NumDofsX(i, j, k) = offset++;
		    NumDofsY(i, j, k) = offset++;
		    NumDofsZ(i, j, k) = offset++;
		  }
      }
    
    ValuePhi2D.Reallocate(nb_dof_inside_tri, this->Points2D_tri().GetM());
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
		ValuePhi2D(nb, k) = val;
		nb++;
	      }
	    
	    pow_one_minus_y *= 1.0-y;
	  }
      }
    
    for (int i = 0; i < this->Points2D_tri().GetM(); i++)
      {
        Real_wp x = this->Points2D_tri()(i)(0);
        Real_wp y = this->Points2D_tri()(i)(1);
        Real_wp z = 0.0;
        R3 pt(x, y, z); VectR3 phi;
        ComputeValuesPhiRef(pt, phi);
        for (int k = 0; k < nb_dof_inside_tri; k++)
          {
            ValuePhi2D(k, i) = DotProd(phi(k), this->NormaleLoc(0));
          }
      }
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point_loc
  */
  void TetrahedronHdivOptimalHpFirstFamily
  ::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    VectReal_wp Pxy, Pxz, Pyz, Pxyz;
    Real_wp x_div_omz(0), x_div_omy(0), y_div_omz(0), x_div_omyz(0);

    if (abs(y+z - 1.0) > epsilon_machine)
      x_div_omyz = 2.0*x/(1.0-y-z) - 1.0;
      
    if (abs(y - 1.0) > epsilon_machine)
      x_div_omy = 2.0*x/(1.0-y) - 1.0;
    
    if (abs(z - 1.0) > epsilon_machine)
      {
	x_div_omz = 2.0*x/(1.0-z) - 1.0;
	y_div_omz = 2.0*y/(1.0-z) - 1.0;
      }    
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x_div_omy, Pxy);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x_div_omz, Pxz);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y_div_omz, Pyz);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x_div_omyz, Pxyz);
    
    Vector<VectReal_wp> Py(order), Pz(order), Pyz2(order);
    VectReal_wp pow_OneMinusY(order+1), pow_OneMinusZ(order+1), pow_OneMinusYZ(order+1);
    pow_OneMinusY(0) = 1.0; pow_OneMinusZ(0) = 1.0; pow_OneMinusYZ(0) = 1.0;
    for (int i = 0; i < order; i++)
      {
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, 2.0*y-1, Py(i));
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, 2.0*z-1, Pz(i));
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, y_div_omz, Pyz2(i));
	pow_OneMinusY(i+1) = pow_OneMinusY(i)*(1.0-y);
	pow_OneMinusZ(i+1) = pow_OneMinusZ(i)*(1.0-z);
	pow_OneMinusYZ(i+1) = pow_OneMinusYZ(i)*(1.0-y-z);
      }
    
    // faces triangulaires
    int node = 0, N = order*(order+1)/2; Real_wp val;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
        {		
	  val = Pxy(i)*Py(i)(j)*pow_OneMinusY(i);
          phi(node).Init(-val*x, -val*y, val*(1.0-z));
	  
	  val = Pxz(i)*Pz(i)(j)*pow_OneMinusZ(i);
	  phi(node + N).Init(val*x, val*(y-1.0), val*z);
	  
	  val = Pyz(i)*Pz(i)(j)*pow_OneMinusZ(i);
	  phi(node + 2*N).Init(val*(1.0-x), -val*y, -val*z);
	  phi(node + 3*N).Init(val*x, val*y, val*z);
	  
	  node++;
        }	
    
    node += 3*N;
    
    // interior
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        for (int k = 0; k <= order-2-i-j; k++)
          {
	    val = Pxyz(i)*pow_OneMinusYZ(i)*Pyz2(i)(j)*pow_OneMinusZ(j)*Pz(i+j+1)(k);
	    phi(node).Init(z*val*x, z*val*y, z*val*(z-1.0));
	    phi(node+1).Init(y*val*x, y*val*(y-1.0), y*val*z);
	    phi(node+2).Init(x*val*(x-1.0), x*val*y, x*val*z);
	    node += 3;
          }	
    
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void TetrahedronHdivOptimalHpFirstFamily
  ::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z; R3 grad; Real_wp coef;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    VectReal_wp Pxy, Pxz, Pyz, Pxyz, dPxy, dPxz, dPyz, dPxyz;
    Real_wp x_div_omz(0), x_div_omy(0), y_div_omz(0), x_div_omyz(0);
    Real_wp dxy_dx(0), dxy_dy(0), dxz_dx(0), dxz_dz(0), dyz_dy(0), dyz_dz(1);
    Real_wp dxyz_dx(0), dxyz_dz(0);
    
    if (abs(y+z - 1.0) > epsilon_machine)
      {
	x_div_omyz = 2.0*x/(1.0-y-z) - 1.0;
	dxyz_dx = 2.0/(1.0-y-z);
	dxyz_dz = 2.0*x/square(1.0-y-z);
      }
      
    if (abs(y - 1.0) > epsilon_machine)
      {
	x_div_omy = 2.0*x/(1.0-y) - 1.0;
	dxy_dx = 2.0/(1.0-y);
	dxy_dy = 2.0*x/square(1.0-y);
      }
    
    if (abs(z - 1.0) > epsilon_machine)
      {
	x_div_omz = 2.0*x/(1.0-z) - 1.0;
	y_div_omz = 2.0*y/(1.0-z) - 1.0;
	dxz_dx = 2.0/(1.0-z);
	dxz_dz = 2.0*x/square(1.0-z);
	dyz_dy = 2.0/(1.0-z);
	dyz_dz = 2.0*y/square(1.0-z);
      }    
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x_div_omy, Pxy, dPxy);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x_div_omz, Pxz, dPxz);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, y_div_omz, Pyz, dPyz);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, x_div_omyz, Pxyz, dPxyz);
    
    Vector<VectReal_wp> Py(order), Pz(order), Pyz2(order);
    Vector<VectReal_wp> dPy(order), dPz(order), dPyz2(order);
    VectReal_wp pow_OneMinusY(order+1), pow_OneMinusZ(order+1), pow_OneMinusYZ(order+1);
    pow_OneMinusY(0) = 1.0; pow_OneMinusZ(0) = 1.0; pow_OneMinusYZ(0) = 1.0;
    for (int i = 0; i < order; i++)
      {
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, 2.0*y-1, Py(i), dPy(i));
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, 2.0*z-1, Pz(i), dPz(i));
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, y_div_omz, Pyz2(i), dPyz2(i));
	pow_OneMinusY(i+1) = pow_OneMinusY(i)*(1.0-y);
	pow_OneMinusZ(i+1) = pow_OneMinusZ(i)*(1.0-z);
	pow_OneMinusYZ(i+1) = pow_OneMinusYZ(i)*(1.0-y-z);
      }
    
    // faces triangulaires
    int node = 0, N = order*(order+1)/2; Real_wp val;
    for (int i = 0; i <= order-1; i++)
      for (int j = 0; j <= order-1-i; j++)
        {
	  val = Pxy(i)*Py(i)(j)*pow_OneMinusY(i);
	  grad(0) = dxy_dx*dPxy(i)*Py(i)(j)*pow_OneMinusY(i);
	  grad(1) = pow_OneMinusY(i)*(dxy_dy*dPxy(i)*Py(i)(j) + Pxy(i)*2.0*dPy(i)(j));
	  if (i > 0)
	    grad(1) -= Real_wp(i)*pow_OneMinusY(i-1)*Pxy(i)*Py(i)(j);
	  
	  grad(2) = 0.0;
	  res(node) = -3.0*val - x*grad(0) - y*grad(1);
	  
	  val = Pxz(i)*Pz(i)(j)*pow_OneMinusZ(i);
	  grad(0) = dxz_dx*dPxz(i)*Pz(i)(j)*pow_OneMinusZ(i);
	  grad(1) = 0.0;
	  grad(2) = pow_OneMinusZ(i)*(dxz_dz*dPxz(i)*Pz(i)(j) + Pxz(i)*2.0*dPz(i)(j));
	  if (i > 0)
	    grad(2) -= Real_wp(i)*pow_OneMinusZ(i-1)*Pxz(i)*Pz(i)(j);
	  
	  res(node + N) = 3.0*val + x*grad(0) + z*grad(2);
	  
	  val = Pyz(i)*Pz(i)(j)*pow_OneMinusZ(i);
	  grad(0) = 0.0;
	  grad(1) = dyz_dy*dPyz(i)*Pz(i)(j)*pow_OneMinusZ(i);
	  grad(2) = pow_OneMinusZ(i)*(dyz_dz*dPyz(i)*Pz(i)(j) + Pyz(i)*2.0*dPz(i)(j));
	  if (i > 0)
	    grad(2) -= Real_wp(i)*pow_OneMinusZ(i-1)*Pyz(i)*Pz(i)(j);
	  
	  res(node + 2*N) = -3.0*val - y*grad(1) - z*grad(2);	  
	  res(node + 3*N) = 3.0*val + y*grad(1) + z*grad(2);
	  
	  node++;
        }	
    
    node += 3*N;
    
    // interior
    for (int i = 0; i <= order-2; i++)
      for (int j = 0; j <= order-2-i; j++)
        for (int k = 0; k <= order-2-i-j; k++)
          {
	    val = Pxyz(i)*pow_OneMinusYZ(i)*Pyz2(i)(j)*pow_OneMinusZ(j)*Pz(i+j+1)(k);
	    grad(0) = dxyz_dx*dPxyz(i)*pow_OneMinusYZ(i)*Pyz2(i)(j)
	      *pow_OneMinusZ(j)*Pz(i+j+1)(k);
	    
	    grad(1) = pow_OneMinusZ(j)*Pz(i+j+1)(k)*pow_OneMinusYZ(i)*
	      (dxyz_dz*dPxyz(i)*Pyz2(i)(j) + dyz_dy*Pxyz(i)*dPyz2(i)(j));
	    
	    grad(2) = pow_OneMinusZ(j)*Pz(i+j+1)(k)*pow_OneMinusYZ(i)*
	      (dxyz_dz*dPxyz(i)*Pyz2(i)(j) + dyz_dz*Pxyz(i)*dPyz2(i)(j));
	    
	    grad(2) += 2.0*Pxyz(i)*pow_OneMinusYZ(i)*Pyz2(i)(j)
	      *pow_OneMinusZ(j)*dPz(i+j+1)(k);
	    
	    if (i > 0)
	      {
		coef = Real_wp(i)*pow_OneMinusYZ(i-1)*Pxyz(i)*Pyz2(i)(j)
		  *pow_OneMinusZ(j)*Pz(i+j+1)(k);
		grad(1) -= coef;
		grad(2) -= coef;
	      }
	    
	    if (j > 0)
	      grad(2) -= Real_wp(j)*pow_OneMinusZ(j-1)*Pxyz(i)*pow_OneMinusYZ(i)
		*Pyz2(i)(j)*Pz(i+j+1)(k);
	    
	    res(node) = val*(4.0*z - 1.0) + grad(0)*z*x + grad(1)*z*y + grad(2)*z*(z-1.0);
	    res(node+1) = val*(4.0*y - 1.0) + grad(0)*y*x + grad(1)*(y-1.0)*y + grad(2)*y*z;
	    res(node+2) = val*(4.0*x - 1.0) + grad(0)*(x-1.0)*x + grad(1)*x*y + grad(2)*x*z;
	    node += 3;
          }    
  }


  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void TetrahedronHdivOptimalHpFirstFamily
  ::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    //TetrahedronReference<2>::ComputeProjectionDofRef(feval, contrib);
    
    const VectReal_wp& weights3d = this->WeightsND();    
    Vector1 feval_weight = feval;
    for (int i = 0; i < feval_weight.GetM()/3; i++)
      {
	feval_weight(3*i) *= weights3d(i);
	feval_weight(3*i+1) *= weights3d(i);
	feval_weight(3*i+2) *= weights3d(i);
      }
    
    ApplyCh(feval_weight, contrib);
    SolveMassMatrix(contrib);
  }
  
  
  void TetrahedronHdivOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    const VectReal_wp& inv_weight
      = static_cast<const TriangleDgOrtho&>(*this->element_tri_surf).GetInverseWeightFunction();
    
    for (int i = 0; i < contrib.GetM(); i++)
      contrib(i) *= inv_weight(i);
    
    if ((num_loc == 0) || (num_loc == 2))
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  void TetrahedronHdivOptimalHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    const VectReal_wp& inv_weight
      = static_cast<const TriangleDgOrtho&>(*this->element_tri_surf).GetInverseWeightFunction();
    
    for (int i = 0; i < contrib.GetM(); i++)
      contrib(i) *= inv_weight(i);

    if ((num_loc == 0) || (num_loc == 2))
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  void TetrahedronHdivOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  
  void TetrahedronHdivOptimalHpFirstFamily
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

	
  //! displays details of class TetrahedronHdivOptimalHpFirstFamily
  ostream& operator <<(ostream& out, const TetrahedronHdivOptimalHpFirstFamily& e)
  {
    out<<static_cast<const TetrahedronReference<3>& >(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_HDIV_OPTIMAL_HP_FIRST_FAMILY_CXX
#endif
