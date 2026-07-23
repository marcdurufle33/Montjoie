#ifndef MONTJOIE_FILE_WEDGE_HDIV_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  WedgeHdivHpFirstFamily::WedgeHdivHpFirstFamily() : WedgeReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! how to number mesh
  void WedgeHdivHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return WedgeReference<3>::ConstructNumberMap(nmap, dg);

    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, order*order);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
		    
    // dofs inside the wedge
    nmap.SetNbDofWedge(order, nb_dof_loc-nb_dof_boundaries);

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ElementReference<Dimension2, 3>::FindHdivRotationTri(order, this->Points2D_tri(),
							 this->Weights2D_tri(),
							 ValuePhiTri2D, nmap);

    ElementReference<Dimension2, 3>::FindHdivRotationQuad(order, NumQuad2D,
							  ValuePhiQuad2D, nmap);
  }


  size_t WedgeHdivHpFirstFamily::GetMemorySize() const
  {
    size_t taille = WedgeReference<3>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize() + JacobiPolynom11.GetMemorySize();    
    taille += Seldon::GetMemorySize(jacobi_2ip1_pol);
    taille += ValLeg.GetMemorySize() + ValuePhiQuad2D.GetMemorySize() + ValuePhiTri2D.GetMemorySize();
    return taille;
  }

    
  //! construction of finite element
  void WedgeHdivHpFirstFamily::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						    int rsurf_tri, int rsurf_quad, int type_surf_tri,
						    int type_surf_quad, int gauss_z)
  {    
    if (rquad <= r-1)
      rquad = r;
    
    WedgeReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad,
					      r, r, -1, -1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
        
    // construction of Value_PhiVec, Div_Phi and elementary matrices
    ConstructHdivElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace();

    TriangleDgOrtho* Fb_tri = new TriangleDgOrtho();
    element_tri_surf = Fb_tri;
    Fb_tri->ConstructFiniteElement(order-1, order, order, TriangleQuadrature::QUADRATURE_GAUSS);

    QuadrangleDgOrtho* Fb_quad = new QuadrangleDgOrtho();
    element_quad_surf = Fb_quad;
    Fb_quad->ConstructFiniteElement(order-1, order, order, Globatto<Real_wp>::QUADRATURE_GAUSS);
  }
  
  
  //! construction of basis functions
  void WedgeHdivHpFirstFamily::ConstructFunctions()
  {    
    nb_dof_tri = order*(order+1)/2;
    nb_dof_quad = order*order;
    nb_dof_boundaries = 2*nb_dof_tri + 3*nb_dof_quad;
    nb_dof_loc = order*(order+2)*order + (order+1)*order*(order+1)/2;

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

    ValLeg.Reallocate(order, order+1);
    ValLeg.Fill(0);
    VectReal_wp Px;
    for (int i = 0; i <= order; i++)
      {
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*this->Points1D(i)-1.0, Px);
	for (int j = 0; j < order; j++)
          ValLeg(j, i) = Px(j);
      }
    
    ValuePhiQuad2D.Reallocate(order*order, this->Points2D_quad().GetM());
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
	for (int k1 = 0; k1 <= order; k1++)
	  for (int k2 = 0; k2 <= order; k2++)
	    ValuePhiQuad2D(i*order + j, NumQuad2D(k1, k2)) = ValLeg(i, k1)*ValLeg(j, k2);    

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

  }
  
  
  //! projection of a function on the finite element space of approximation
  /*!
    \param[in] feval evaluation of f on points where dofs are located
    \param[out] contrib  result, dof components of f
  */
  template<class Vector1, class Vector2>
  void WedgeHdivHpFirstFamily
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
    \param[out] phi values of basis functions on point_loc
  */
  void WedgeHdivHpFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
		
    int r = order; Real_wp vloc, val;
    
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    Real_wp a = 2.0*x - 1.0, b = 2.0*y - 1.0, c = 2.0*z - 1.0, m = 2.0*x/(1.0-y)-1.0;
    
    VectReal_wp Px, Py, Pz, Px2;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, r, b, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, r, c, Pz);
    EvaluateJacobiPolynomial(LegendrePolynom, r, m, Px2);
    Vector<VectReal_wp> Py2(r);
    VectReal_wp pow_OneMinusY(order+1);
    pow_OneMinusY(0) = 1.0;
    for (int i = 0; i < r; i++)
      {
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, b, Py2(i));
	pow_OneMinusY(i+1) = pow_OneMinusY(i)*(1.0-y);
      }
    
    // Faces triangulaires
    int num = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	{          
	  vloc = Px2(i)*Py2(i)(j)*pow_OneMinusY(i);	  
	  phi(num).Init(0, 0, (1.0-z)*vloc);
	  num += 1;
	}
    
    // Faces quadrangulaires
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{          
	  vloc = Px(i)*Pz(j);
	  phi(num).Init(x*vloc, (y-1.0)*vloc, 0);
	  num++;
	}
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{          	  
	  vloc = Py(i)*Pz(j);
	  phi(num).Init(x*vloc, y*vloc, 0);
	  num++;
	}
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{          	     
	  vloc = Py(i)*Pz(j); 
	  phi(num).Init(-(x-1.0)*vloc, -y*vloc, 0);
	  num++;
	}
    
    // Faces triangulaires
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	{          
	  vloc = Px2(i)*Py2(i)(j)*pow_OneMinusY(i);	  
	  phi(num).Init(0, 0, z*vloc);
	  num++;
	}
    
    // Interieur 
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	for (int k = 0; k < r; k++)
	  {
	    val = Px2(i)*Py2(i)(j)*pow_OneMinusY(i)*Pz(k);
	    if (i+j < r-1)
	      {
		vloc = y*val;
		phi(num).Init(x*vloc, (y-1.0)*vloc, 0);
		
		vloc = x*val;
		phi(num+1).Init((x-1.0)*vloc, y*vloc, 0);
		num += 2;		
	      }
	    
	    if (k < r-1)
	      {		
		phi(num).Init(0, 0, z*(1.0-z)*val);
		num++;
	      }
	  }    
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void WedgeHdivHpFirstFamily::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    int r = order; Real_wp vloc, val; R3 grad;
    
    Real_wp x = point_loc(0), y = point_loc(1), z = point_loc(2);
    Real_wp a = 2.0*x - 1.0, b = 2.0*y - 1.0, c = 2.0*z - 1.0, m = 2.0*x/(1.0-y)-1.0;
    Real_wp dm_dx = 2.0/(1.0-y), dm_dy = 2.0*x/square(1.0-y);
    
    VectReal_wp Px, Py, Pz, Px2, dPx, dPy, dPz, dPx2;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, r, b, Py, dPy);
    EvaluateJacobiPolynomial(LegendrePolynom, r, c, Pz, dPz);
    EvaluateJacobiPolynomial(LegendrePolynom, r, m, Px2, dPx2);
    Vector<VectReal_wp> Py2(r), dPy2(r);
    VectReal_wp pow_OneMinusY(order+1), dpow_OneMinusY(order+1);
    pow_OneMinusY(0) = 1.0;
    dpow_OneMinusY(0) = 0.0;
    for (int i = 0; i < r; i++)
      {
	EvaluateJacobiPolynomial(jacobi_2ip1_pol(i), order-1-i, b, Py2(i), dPy2(i));
	pow_OneMinusY(i+1) = pow_OneMinusY(i)*(1.0-y);
	if (i > 0)
	  dpow_OneMinusY(i) = -Real_wp(i)*pow_OneMinusY(i-1);
      }
    
    // Faces triangulaires
    int num = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	{          
	  vloc = Px2(i)*Py2(i)(j)*pow_OneMinusY(i);	  
	  res(num) = -vloc;
	  num += 1;
	}
    
    // Faces quadrangulaires
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{          
	  vloc = Px(i)*Pz(j);
	  res(num) = 2.0*x*dPx(i)*Pz(j) + 2.0*vloc;
	  num++;
	}
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{          	  
	  vloc = Py(i)*Pz(j);
	  res(num) = 2.0*vloc + 2.0*y*dPy(i)*Pz(j);
	  num++;
	}
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	{          	     
	  vloc = Py(i)*Pz(j); 
	  res(num) = -2.0*vloc - 2.0*y*dPy(i)*Pz(j);   
	  num++;
	}
    
    // Faces triangulaires
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	{          
	  vloc = Px2(i)*Py2(i)(j)*pow_OneMinusY(i);
	  res(num) = vloc;
	  num++;
	}
    
    // Interieur 
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r-i; j++)
	for (int k = 0; k < r; k++)
	  {
	    val = Px2(i)*Py2(i)(j)*pow_OneMinusY(i)*Pz(k);
	    grad(0) = dm_dx*dPx2(i)*Py2(i)(j)*pow_OneMinusY(i)*Pz(k);
	    grad(1) = dm_dy*dPx2(i)*Py2(i)(j)*pow_OneMinusY(i)*Pz(k)
	      + Px2(i)*Pz(k)*( 2.0*dPy2(i)(j)*pow_OneMinusY(i) +
			       Py2(i)(j)*dpow_OneMinusY(i));
	    grad(2) = 2.0*dPz(k)*Px2(i)*Py2(i)(j)*pow_OneMinusY(i);
	    
	    if (i+j < r-1)
	      {
		res(num) = (3.0*y-1.0)*val + x*y*grad(0) + y*(y-1.0)*grad(1);
		res(num+1) = (3.0*x-1.0)*val + x*(x-1.0)*grad(0) + x*y*grad(1);
		num += 2;
	      }
	    
	    if (k < r-1)
	      {
		res(num) = (1.0 - 2.0*z)*val + z*(1.0-z)*grad(2);
		num++;
	      }
	  }
  }
  
  
  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
    (\omega, \xi_k) are here Gauss points and not Gauss-Lobatto points
  */  
  template<class Vector1, class Vector2>
  void WedgeHdivHpFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }


  void WedgeHdivHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    const VectReal_wp& inv_weight
      = static_cast<const TriangleDgOrtho&>(*this->element_tri_surf).GetInverseWeightFunction();

    switch(num_loc)
      {
      case 0:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = -inv_weight(i)*contrib(i);
	
	break;
      case 3:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = -contrib(i);
	
	break;
      case 4:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = inv_weight(i)*contrib(i);
	
	break;
      }
  }
  

  void WedgeHdivHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    const VectReal_wp& inv_weight
      = static_cast<const TriangleDgOrtho&>(*this->element_tri_surf).GetInverseWeightFunction();

    switch(num_loc)
      {
      case 0:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = -inv_weight(i)*contrib(i);
	
	break;
      case 3:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = -contrib(i);
	
	break;
      case 4:
	for (int i = 0; i < contrib.GetM(); i++)
	  contrib(i) = inv_weight(i)*contrib(i);
	
	break;
      }
  }
  
  
  //! displays details of class WedgeHdivHpFirstFamily
  ostream& operator <<(ostream& out, const WedgeHdivHpFirstFamily& e)
  {
    out << static_cast<const WedgeReference<3>& >(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_HDIV_HP_FIRST_FAMILY_CXX
#endif
