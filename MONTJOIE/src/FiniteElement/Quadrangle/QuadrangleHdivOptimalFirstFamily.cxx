#ifndef MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_FIRST_FAMILY_CXX

namespace Montjoie
{
  //! default constructor
  QuadrangleHdivOptimalFirstFamily::QuadrangleHdivOptimalFirstFamily() : QuadrangleReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void QuadrangleHdivOptimalFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<3>::ConstructNumberMap(nmap, dg);

    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 2*order*order);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t QuadrangleHdivOptimalFirstFamily::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<3>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize();
    taille += NumDofs_Y.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += tangente_dof.GetMemorySize();
    taille += GL_Gexact.GetMemorySize();
    return taille;
  }

      
  //! constructing finite element
  void QuadrangleHdivOptimalFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    
    QuadrangleReference<3>::ConstructFiniteElement(r, rgeom, r+1, type_quad);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructHdivElementaryMatrix();
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();          
  }
  
  
  //! construction of basis functions
  void QuadrangleHdivOptimalFirstFamily::ConstructFunctions()
  {
    EdgesDof.Reallocate(order, 4);
    NumDofs_X.Reallocate(order+2, order);
    NumDofs_Y.Reallocate(order, order+2);
    nb_dof_loc = order*(order+2)*2;
    tangente_dof.Reallocate(nb_dof_loc);
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 0) = i;
	NumDofs_Y(i, 0) = i;
	tangente_dof(i).Init(0.0, -1.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 1) = i+order;
	NumDofs_X(order+1, i) = i+order;
	tangente_dof(order+i).Init(1.0, 0.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 2) = i+2*order;
	NumDofs_Y(order-1-i, order+1) = EdgesDof(i, 2);
	tangente_dof(i+2*order).Init(0.0, 1.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 3) = i + 3*order;
	NumDofs_X(0, order-1-i) = EdgesDof(i, 3);
	tangente_dof(i+3*order).Init(-1.0, 0.0);
      }

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    int numero = 4*order;
    nb_dof_boundaries = numero;
    
    // horizontal dofs inside
    for (int i = 0; i < order; i++)
      for (int j = 1; j <= order; j++)
	{
	  tangente_dof(numero).Init(0.0, -1.0);
	  NumDofs_Y(i, j) = numero++;
	}
    
    // vertical dofs inside
    for (int i = 1; i <= order; i++)
      for (int j = 0; j < order; j++)
	{
	  tangente_dof(numero).Init(1.0, 0.0);
	  NumDofs_X(i, j) = numero++;
	}
    
    nb_dof_loc = numero;

    VectReal_wp points_lob, weights_lob;
    ComputeGaussLobatto(points_lob, weights_lob, order+1);

    VectReal_wp points_dof1d; VectR2 points_dof2d;
    points_dof1d.Reallocate(order);
    for (int i = 0; i < order; i++)
      points_dof1d(i) = points_lob(i+1);
    
    lob_basis.AffectPoints(points_dof1d);
    
    points_dof2d.Reallocate(nb_dof_loc);
    CoordinateDofs.Reallocate(nb_dof_loc,3);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	{
	  points_dof2d(NumDofs_Y(i, j)).Init(lob_basis.Points(i), lob_quad.Points(j));
	  points_dof2d(NumDofs_X(j, i)).Init(lob_quad.Points(j), lob_basis.Points(i));
	  
	  // the two first indices are coordinates in space
	  // the third index is 0 if it is oriented by ex, 1 for ey
	  CoordinateDofs(NumDofs_Y(i, j), 0) = i;
	  CoordinateDofs(NumDofs_Y(i, j), 1) = j;
	  CoordinateDofs(NumDofs_Y(i, j), 2) = 1;
	  CoordinateDofs(NumDofs_X(j, i), 0) = j;
	  CoordinateDofs(NumDofs_X(j, i), 1) = i;
	  CoordinateDofs(NumDofs_X(j, i), 2) = 0;
	}

    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    this->nb_points_dof_inside = nb_dof_loc;
    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(order);
	for (int i = 0; i < order; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }            
  }
  
  
  //! computation of stiffness matrix
  void QuadrangleHdivOptimalFirstFamily::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
  }
  
  
  //! computation of mass matrix
  void QuadrangleHdivOptimalFirstFamily::ConstructMassMatrix()
  {
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
    GL_Gexact.Reallocate(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
	GL_Gexact(i, j) = lob_quad.EvaluatePhi(i, gauss.Points(j));
  }

  
  //! Integration against basis functions on an edge
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc edge number
  */  
  template<class Vector1, class Vector2>
  void QuadrangleHdivOptimalFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    Vector1 feval2(2*(order+2));
    feval2.Fill(0);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
        {
          feval2(2*i) += GL_Gexact(i, j)*feval(2*j);
          feval2(2*i+1) += GL_Gexact(i, j)*feval(2*j+1);
        }
    
    ComputeIntegralSurfaceRef(feval2, res, num_loc);
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void QuadrangleHdivOptimalFirstFamily
  ::ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc);
    FillZero(res);
    VectReal_wp phiGL_x(order+2), phiGL_y(order+2), phiG_x(order), phiG_y(order);
    for (int i = 0; i < order; i++)
      {
	phiG_x(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiG_y(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i <= order+1; i++)
      {
	phiGL_x(i) = lob_quad.EvaluatePhi(i, point_loc(0));
	phiGL_y(i) = lob_quad.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	{
	  int num_dof = NumDofs_Y(i, j);
	  res(num_dof)(1) = phiG_x(i)*phiGL_y(j);
	  if (j == 0)
	    res(num_dof)(1) = -res(num_dof)(1);
	  
	  num_dof = NumDofs_X(j, i);
	  res(num_dof)(0) = phiGL_x(j)*phiG_y(i);
	  if (j == 0)
	    res(num_dof)(0) = -res(num_dof)(0);
	}
  }
  
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void QuadrangleHdivOptimalFirstFamily
  ::ComputeDivPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    VectReal_wp dphiGL_x(order+2), dphiGL_y(order+2), phiG_x(order), phiG_y(order);
    for (int i = 0; i < order; i++)
      {
	phiG_x(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiG_y(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i <= order+1; i++)
      {
	dphiGL_x(i) = lob_quad.EvaluatePhiGrad(i, point_loc(0));
	dphiGL_y(i) = lob_quad.EvaluatePhiGrad(i, point_loc(1));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	{
	  int num_dof = NumDofs_Y(i, j);
	  res(num_dof) = dphiGL_y(j)*phiG_x(i);
	  if (j == 0)
	    res(num_dof) = -res(num_dof);
	  
	  num_dof = NumDofs_X(j, i);
	  res(num_dof) = dphiGL_x(j)*phiG_y(i);
	  if (j==0)
	    res(num_dof) = -res(num_dof);
      }
  }
  

  //! displays informations about class QuadrangleHdivFirstFamily
  ostream& operator <<(ostream& out, const QuadrangleHdivOptimalFirstFamily& e)
  {
    out<<static_cast<const QuadrangleReference<3>& >(e);    
    return out;
  }

}
  
#define MONTJOIE_FILE_QUADRANGLE_HDIV_OPTIMAL_FIRST_FAMILY_CXX
#endif
