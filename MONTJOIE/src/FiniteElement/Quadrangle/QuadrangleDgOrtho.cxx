#ifndef MONTJOIE_FILE_QUADRANGLE_DG_ORTHO_CXX

namespace Montjoie
{
  
  QuadrangleDgOrtho::QuadrangleDgOrtho()
  {
  }


  size_t QuadrangleDgOrtho::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<1>::GetMemorySize();
    taille += LegendrePolynom.GetMemorySize();
    taille += InvWeightPolynomial.GetMemorySize();
    taille += NumFct2D.GetMemorySize();
    return taille;
  }

  
  // construction of finite element
  void QuadrangleDgOrtho
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    QuadrangleReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    order = r;
    nb_dof_loc = (r+1)*(r+1);

    this->ConstructFunctions();
    
    this->SetPointsDof1D(this->Points1D());
    this->SetPointsDofND(this->PointsND());

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;

    this->elt_geom.dof_equal_quadrature = true;

    this->Fb_geom.ComputeCoefficientTransformation();    
    
    this->ConstructElementaryMatrix(*this);
  }
    

  void QuadrangleDgOrtho::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    // no needed since this element is discontinuous
  }

  
  void QuadrangleDgOrtho::ConstructFunctions()
  {
    GetJacobiPolynomial(LegendrePolynom, order+1, Real_wp(0), Real_wp(0));

    int num = 0;
    NumFct2D.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	NumFct2D(i, j) = num++;
    
    VectReal_wp xi, omega, CoefLegendre;
    CoefLegendre.Reallocate(order+1); CoefLegendre.Fill(0);
    ComputeGaussLegendre(xi, omega, order);
    for (int q = 0; q <= order; q++)
      {
	VectReal_wp Pn;
	EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*xi(q) - 1.0, Pn);
	for (int i = 0; i <= order; i++)
	  CoefLegendre(i) += omega(q)*Pn(i)*Pn(i);
      }
    
    // we keep 1/sqrt(omega)
    for (int i = 0; i <= order; i++)
      CoefLegendre(i) = 1.0/sqrt(CoefLegendre(i));
    
    InvWeightPolynomial.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	InvWeightPolynomial(this->NumFct2D(i, j)) = CoefLegendre(i)*CoefLegendre(j);
  }
  

  template<class Vector1, class Vector2>
  void QuadrangleDgOrtho::ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    Vector1 feval_weight = feval;
    const VectReal_wp& weights2d = this->WeightsND();
    for (int i = 0; i < this->nb_points_dof_inside; i++)
      feval_weight(i) *= weights2d(i);

    ApplyCh(feval_weight, contrib);
    for (int i = 0; i < contrib.GetM(); i++)
      contrib(i) *= square(InvWeightPolynomial(i));
  }
  

  void QuadrangleDgOrtho
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }
  

  void QuadrangleDgOrtho
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    ComputeProjectionDofGen(feval, contrib);
  }

  
  void QuadrangleDgOrtho::ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const
  {
    VectReal_wp Px, Py;
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*pointloc(0) - 1.0, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*pointloc(1) - 1.0, Py);
    
    phi.Reallocate(nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	phi(this->NumFct2D(i, j)) = Px(i)*Py(j);

  }
  
  
  void QuadrangleDgOrtho::ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const
  {
    VectReal_wp Px, Py, dPx, dPy;
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*pointloc(0) - 1.0, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, order, 2.0*pointloc(1) - 1.0, Py, dPy);
    
    grad_phi.Reallocate(nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	grad_phi(this->NumFct2D(i, j)).Init(2.0*dPx(i)*Py(j), 2.0*Px(i)*dPy(j));
    
  }
  
}

#define MONTJOIE_FILE_QUADRANGLE_DG_ORTHO_CXX
#endif

