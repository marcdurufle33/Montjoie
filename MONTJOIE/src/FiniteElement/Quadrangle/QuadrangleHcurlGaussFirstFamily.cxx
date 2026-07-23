#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_GAUSS_FIRST_FAMILY_CXX

namespace Montjoie
{
  //! default constructor
  QuadrangleHcurlGaussFirstFamily::QuadrangleHcurlGaussFirstFamily() : QuadrangleReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void QuadrangleHcurlGaussFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<2>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 2*order*(order-1));
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t QuadrangleHcurlGaussFirstFamily::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<2>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize();
    taille += NumDofs_Y.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += normale_dof.GetMemorySize();
    taille += lob_ortho.GetMemorySize();
    return taille;
  }
  

  //! constructing finite element
  void QuadrangleHcurlGaussFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    QuadrangleReference<2>::ConstructFiniteElement(r, rgeom, r, type_quad);
    
    lob_basis.ConstructQuadrature(order-1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    lob_ortho.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    lob_ortho.ComputeGradPhi(1e3*epsilon_machine);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructHcurlElementaryMatrix();

    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r-1, r, r, type_quad, EdgeGauss::GAUSS);
    edge->SetDofPoints(this->PointsDof1D());
    this->element_surface = edge;    
  }


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  void QuadrangleHcurlGaussFirstFamily
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
  }
  
  
  //! construction of basis functions
  void QuadrangleHcurlGaussFirstFamily::ConstructFunctions()
  {
    this->nb_dof_H = order*order;
    EdgesDof.Reallocate(order, 4);
    NumDofs_X.Reallocate(order, order+1);
    NumDofs_Y.Reallocate(order+1, order);
    nb_dof_loc = order*(order+1)*2;
    normale_dof.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 0) = i;
	NumDofs_X(i, 0) = i;
	normale_dof(i).Init(0.0, -1.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 1) = i+order;
	NumDofs_Y(order, i) = i+order;
	normale_dof(order+i).Init(1.0, 0.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 2) = i+2*order;
	NumDofs_X(order-1-i, order) = EdgesDof(i, 2);
	normale_dof(i+2*order).Init(0.0, 1.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 3) = i + 3*order;
	NumDofs_Y(0, order-1-i) = EdgesDof(i, 3);
	normale_dof(i+3*order).Init(-1.0, 0.0);
      }

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);    
    int numero = 4*order;
    nb_dof_boundaries = numero;
    
    // horizontal dofs inside
    for (int i = 0; i < order; i++)
      for (int j = 1; j < order; j++)
	{
	  normale_dof(numero).Init(0.0, -1.0);
	  NumDofs_X(i,j) = numero++;
	}
    
    // vertical dofs inside
    for (int i = 1; i < order; i++)
      for (int j = 0; j < order; j++)
	{
	  normale_dof(numero).Init(1.0, 0.0);
	  NumDofs_Y(i,j) = numero++;
	}
    
    nb_dof_loc = numero;
    
    VectReal_wp points_dof1d; VectR2 points_dof2d;
    points_dof1d = lob_basis.Points();
    points_dof2d.Reallocate(nb_dof_loc);
    CoordinateDofs.Reallocate(nb_dof_loc, 3);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  points_dof2d(NumDofs_X(i,j)).Init(lob_basis.Points(i), lob_ortho.Points(j));
	  points_dof2d(NumDofs_Y(j,i)).Init(lob_ortho.Points(j), lob_basis.Points(i));
	  
	  // the two first indices are coordinates in space
	  // the third index is 0 if it is oriented by ex, 1 for ey
	  CoordinateDofs(NumDofs_X(i,j), 0) = i;
	  CoordinateDofs(NumDofs_X(i,j), 1) = j;
	  CoordinateDofs(NumDofs_X(i,j), 2) = 0;
	  CoordinateDofs(NumDofs_Y(j,i), 0) = j;
	  CoordinateDofs(NumDofs_Y(j,i), 1) = i;
	  CoordinateDofs(NumDofs_Y(j,i), 2) = 1;
	}
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    lob_quad.ComputeGradPhi(1e3*epsilon_machine);

    this->nb_points_dof_inside = points_dof2d.GetM();
    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(order);
	for (int i = 0; i < order; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] res dof components
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlGaussFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& res) const
  {
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  if (j == order)
	    res(NumDofs_X(i,j)) = -feval(2*NumDofs_X(i,j));
	  else
	    res(NumDofs_X(i,j)) = feval(2*NumDofs_X(i,j));
	  
	  if (j == 0)
	    res(NumDofs_Y(j,i)) = -feval(2*NumDofs_Y(j,i)+1);
	  else
	    res(NumDofs_Y(j,i)) = feval(2*NumDofs_Y(j,i)+1);
	}
  }

    
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void QuadrangleHcurlGaussFirstFamily
  ::ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc);
    FillZero(res);
    VectReal_wp phiGL_x(order+1), phiGL_y(order+1), phiG_x(order), phiG_y(order);
    for (int i = 0; i < order; i++)
      {
	phiG_x(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiG_y(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i <= order; i++)
      {
	phiGL_x(i) = lob_ortho.EvaluatePhi(i, point_loc(0));
	phiGL_y(i) = lob_ortho.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  int num_dof = NumDofs_X(i, j);
	  res(num_dof)(0) = phiG_x(i)*phiGL_y(j);
	  if (j == order)
	    res(num_dof)(0) = -res(num_dof)(0);
	  
	  num_dof = NumDofs_Y(j, i);
	  res(num_dof)(1) = phiGL_x(j)*phiG_y(i);
	  if (j == 0)
	    res(num_dof)(1) = -res(num_dof)(1);
	}
  }
  

  //! values of H
  void QuadrangleHcurlGaussFirstFamily
  ::ComputeValuesPhiHRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(order*order);
    FillZero(res);
    VectReal_wp phiG_x(order), phiG_y(order);
    for (int i = 0; i < order; i++)
      {
	phiG_x(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiG_y(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
        res(i*order+j) = phiG_x(i) * phiG_y(j);
  }
    
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void QuadrangleHcurlGaussFirstFamily
  ::ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    VectReal_wp dphiGL_x(order+1), dphiGL_y(order+1), phiG_x(order), phiG_y(order);
    for (int i = 0; i < order; i++)
      {
	phiG_x(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiG_y(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i <= order; i++)
      {
	dphiGL_x(i) = lob_ortho.EvaluatePhiGrad(i, point_loc(0));
	dphiGL_y(i) = lob_ortho.EvaluatePhiGrad(i, point_loc(1));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  int num_dof = NumDofs_X(i, j);
	  res(num_dof) = -dphiGL_y(j)*phiG_x(i);
	  if (j == order)
	    res(num_dof) = -res(num_dof);
	  
	  num_dof = NumDofs_Y(j, i);
	  res(num_dof) = dphiGL_x(j)*phiG_y(i);
	  if (j==0)
	    res(num_dof) = -res(num_dof);
      }
  }
  
  
  //! displays informations about class QuadrangleHcurlGaussFirstFamily
  ostream& operator <<(ostream& out, const QuadrangleHcurlGaussFirstFamily& e)
  {
    out<<static_cast<const QuadrangleReference<2>&>(e);    
    return out;
  }

}
  
#define MONTJOIE_FILE_QUADRANGLE_HCURL_GAUSS_FIRST_FAMILY_CXX
#endif
