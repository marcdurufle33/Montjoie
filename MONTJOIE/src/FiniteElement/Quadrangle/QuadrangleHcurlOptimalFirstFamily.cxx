#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_FIRST_FAMILY_CXX

namespace Montjoie
{
  //! default constructor
  QuadrangleHcurlOptimalFirstFamily::QuadrangleHcurlOptimalFirstFamily(bool proj) : QuadrangleReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    project_with_quadrature = proj;
  }
  
  
  //! how to number mesh
  void QuadrangleHcurlOptimalFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<2>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 2*order*order);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t QuadrangleHcurlOptimalFirstFamily::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<2>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize();
    taille += NumDofs_Y.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += normale_dof.GetMemorySize();
    taille += dGL_GL.GetMemorySize();
    taille += GL_Gexact.GetMemorySize();
    taille += ValGaussExtremity.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += ch_loc.GetMemorySize();
    return taille;
  }

    
  //! constructing finite element
  void QuadrangleHcurlOptimalFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    
    QuadrangleReference<2>::ConstructFiniteElement(r, rgeom, r+1, type_quad);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructHcurlElementaryMatrix();
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();      
  }


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  void QuadrangleHcurlOptimalFirstFamily
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r-1, r, r+1, Globatto<Real_wp>::QUADRATURE_LOBATTO,
				 EdgeGauss::LOBATTO_INT);
    
    this->element_surface = edge;
  }  
  
  //! construction of basis functions
  void QuadrangleHcurlOptimalFirstFamily::ConstructFunctions()
  {
    EdgesDof.Reallocate(order, 4);
    NumDofs_X.Reallocate(order, order+2);
    NumDofs_Y.Reallocate(order+2, order);
    nb_dof_loc = order*(order+2)*2;
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
	NumDofs_Y(order+1, i) = i+order;
	normale_dof(order+i).Init(1.0, 0.0);
      }
    
    for (int i = 0; i < order; i++)
      {
	EdgesDof(i, 2) = i+2*order;
	NumDofs_X(order-1-i, order+1) = EdgesDof(i, 2);
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
      for (int j = 1; j <= order; j++)
	{
	  normale_dof(numero).Init(0.0, -1.0);
	  NumDofs_X(i, j) = numero++;
	}
    
    // vertical dofs inside
    for (int i = 1; i <= order; i++)
      for (int j = 0; j < order; j++)
	{
	  normale_dof(numero).Init(1.0, 0.0);
	  NumDofs_Y(i, j) = numero++;
	}
    
    nb_dof_loc = numero;
    
    VectReal_wp points_dof1d; VectR2 points_dof2d;
    VectReal_wp points_lob, weights_lob;
    ComputeGaussLobatto(points_lob, weights_lob, order+1);
    points_dof1d.Reallocate(order);
    for (int i = 0; i < order; i++)
      points_dof1d(i) = points_lob(i+1);
    
    lob_basis.AffectPoints(points_dof1d);
    
    if (!project_with_quadrature)
      points_dof2d.Reallocate(nb_dof_loc);
    
    CoordinateDofs.Reallocate(nb_dof_loc, 3);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	{
	  if (!project_with_quadrature)
	    {
	      points_dof2d(NumDofs_X(i, j)).Init(lob_basis.Points(i), lob_quad.Points(j));
	      points_dof2d(NumDofs_Y(j, i)).Init(lob_quad.Points(j), lob_basis.Points(i));
	    }
	  
	  // the two first indices are coordinates in space
	  // the third index is 0 if it is oriented by ex, 1 for ey
	  CoordinateDofs(NumDofs_X(i, j), 0) = i;
	  CoordinateDofs(NumDofs_X(i, j), 1) = j;
	  CoordinateDofs(NumDofs_X(i, j), 2) = 0;
	  CoordinateDofs(NumDofs_Y(j, i), 0) = j;
	  CoordinateDofs(NumDofs_Y(j, i), 1) = i;
	  CoordinateDofs(NumDofs_Y(j, i), 2) = 1;
	}

    if (project_with_quadrature)
      {
	points_dof1d = points_lob;
	points_dof2d = this->PointsND();
	this->nb_points_dof_inside = nb_points_quadrature_inside;
      }
    else
      this->nb_points_dof_inside = points_dof2d.GetM();
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    if (project_with_quadrature)
      this->num_dof_points_surf = this->num_quad_points_surf;
    else
      {
	this->num_dof_points_surf.Reallocate(4);
	for (int n = 0; n < 4; n++)
	  {
	    this->num_dof_points_surf(n).Reallocate(order);
	    for (int i = 0; i < order; i++)
	      this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
	  }
      }
  }
  
  
  //! computation of stiffness matrix
  void QuadrangleHcurlOptimalFirstFamily::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    dGL_GL = lob_quad.GradPhi();
    
    int nb_points_fct_H1 = (order+2)*(order+2);
    Matrix<Real_wp, General, ArrayRowSparse> Ch(nb_dof_loc, 2*nb_points_fct_H1);    

    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order+1; j++)
	for (int k = 0; k < order; k++)
	  {
	    int node = NumQuad2D(i, j);
	    int num_dof = NumDofs_X(k, j);
            if (i == order+1)
              {
                Real_wp vloc = ValGaussExtremity(k, 1);
                if (j == order+1)
                  vloc = -vloc;
                
                Ch.AddInteraction(num_dof, 2*node, vloc);
              }
            else if (i == 0)
              {
                Real_wp vloc = ValGaussExtremity(k, 0);
                if (j == order+1)
                  vloc = -vloc;
                
                Ch.AddInteraction(num_dof, 2*node, vloc);
              }
            else
              {
                if (k+1 == i)
                  {
                    if (j == order+1)
                      Ch.AddInteraction(num_dof, 2*node, -1.0);
                    else
                      Ch.AddInteraction(num_dof, 2*node, 1.0);
                  }
              }
            
	    num_dof = NumDofs_Y(i, k);
            if (j == order+1)
              {
                Real_wp vloc = ValGaussExtremity(k, 1);
                if (i == 0)
                  vloc = -vloc;
                
                Ch.AddInteraction(num_dof, 2*node+1, vloc);
              }
            else if (j == 0)
              {
                Real_wp vloc = ValGaussExtremity(k, 0);
                if (i == 0)
                  vloc = -vloc;
                
                Ch.AddInteraction(num_dof, 2*node+1, vloc);
              }
            else
              {
                if (k+1 == j)
                  {
                    if (i == 0)
                      Ch.AddInteraction(num_dof, 2*node+1, -1.0);
                    else
                      Ch.AddInteraction(num_dof, 2*node+1, 1.0);
                  }
              }
          }
    
    Copy(Ch, ch_loc);
  }
  
  
  //! computation of mass matrix
  void QuadrangleHcurlOptimalFirstFamily::ConstructMassMatrix()
  {
    ValGaussExtremity.Reallocate(order, 2);
    for (int i = 0; i < order; i++)
      {
	ValGaussExtremity(i,0) = lob_basis.EvaluatePhi(i, 0.0);
	ValGaussExtremity(i,1) = lob_basis.EvaluatePhi(i, 1.0);
      }
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
    GL_Gexact.Reallocate(order+2, order+1);
    for (int i = 0; i <= order+1; i++)
      for (int j = 0; j <= order; j++)
	GL_Gexact(i, j) = lob_quad.EvaluatePhi(i, gauss.Points(j));
    
  }
  
  
  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] res dof components
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlOptimalFirstFamily::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& res) const
  {
    if (project_with_quadrature)
      for (int i = 0; i < order; i++)
	for (int j = 0; j <= order+1; j++)
	  {
	    if (j == order+1)
	      res(NumDofs_X(i,j)) = -feval(2*NumQuad2D(i+1,j));
	    else
	      res(NumDofs_X(i,j)) = feval(2*NumQuad2D(i+1,j));
	    
	    if (j == 0)
	      res(NumDofs_Y(j,i)) = -feval(2*NumQuad2D(j,i+1)+1);
	    else
	      res(NumDofs_Y(j,i)) = feval(2*NumQuad2D(j,i+1)+1);
	  }
    else
      for (int i = 0; i < order; i++)
	for (int j = 0; j <= order+1; j++)
	  {
	    if (j == order+1)
	      res(NumDofs_X(i,j)) = -feval(2*NumDofs_X(i,j));
	    else
	      res(NumDofs_X(i,j)) = feval(2*NumDofs_X(i,j));
	    
	    if (j == 0)
	      res(NumDofs_Y(j,i)) = -feval(2*NumDofs_Y(j,i)+1);
	    else
	      res(NumDofs_Y(j,i)) = feval(2*NumDofs_Y(j,i)+1);
	  }
  }
  
  
  //! projection from quadrature points (Gauss-Lobatto points) to dof components
  /*!
    \param[in] Equad values on quadrature points
    \param[out] Edof dof components
  */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlOptimalFirstFamily::
  ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const
  {
    FillZero(Edof);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order+1; j++)
	{
	  // dofs along x
	  int num_dof = NumDofs_X(i, j);
	  if (j == order)
            Edof(num_dof) = -Equad(2*NumQuad2D(i+1, j));
          else
            Edof(num_dof) = Equad(2*NumQuad2D(i+1, j));
	  
	  // dofs along y
	  num_dof = NumDofs_Y(j,i);
	  if (j == 0)
            Edof(num_dof) = -Equad(2*NumQuad2D(j, i+1)+1);
          else
            Edof(num_dof) = Equad(2*NumQuad2D(j, i+1)+1);
	}
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
  void QuadrangleHcurlOptimalFirstFamily::
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
  
  
  //! Integration against basis functions on an edge
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc edge number
  */  
  template<class T>
  void QuadrangleHcurlOptimalFirstFamily::
  ApplyShGen(const T& alpha, const Vector<T>& feval, Vector<T>& res, int num_loc, int r) const
  {
    int num_dof;
    if (num_loc == 0)
      {
	for (int j = 0; j <= order+1; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i, 0);
		if (j == 0)
                  res(num_dof) += alpha*ValGaussExtremity(i,0)*feval(2*j);
                else if (j == order+1)
                  res(num_dof) += alpha*ValGaussExtremity(i,1)*feval(2*j);
                else
                  {
                    if (j == i+1)
                      res(num_dof) += alpha*feval(2*j);
                  }
                
		num_dof = NumDofs_Y(j,i);
		if (j == 0)
		  res(num_dof) -= alpha*ValGaussExtremity(i,0)*feval(2*j+1);
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,0)*feval(2*j+1);
	      }
	  }
      }
    else if (num_loc == 1)
      {
	for (int j = 0; j <= order+1; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,j);
		if (j == order+1)
		  res(num_dof) -= alpha*ValGaussExtremity(i,1)*feval(2*j);
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,1)*feval(2*j);
		
                num_dof = NumDofs_Y(order+1, i);
		if (j == 0)
                  res(num_dof) += alpha*ValGaussExtremity(i,0)*feval(2*j+1);
                else if (j == order+1)
                  res(num_dof) += alpha*ValGaussExtremity(i,1)*feval(2*j+1);
                else
                  {
                    if (j == i+1)
                      res(num_dof) += alpha*feval(2*j+1);
                  }
	      }
	  }
      }
    else if (num_loc == 2)
      {
	for (int j = 0; j <= order+1; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,order+1);
		if (j == 0)
                  res(num_dof) -= alpha*ValGaussExtremity(i,0)*feval(2*(order+1-j));
                else if (j == order+1)
                  res(num_dof) -= alpha*ValGaussExtremity(i,1)*feval(2*(order+1-j));
                else
                  {
                    if (j == i+1)
                      res(num_dof) -= alpha*feval(2*(order+1-j));
                  }
		
                num_dof = NumDofs_Y(j,i);
		if (j == 0)
		  res(num_dof) -= alpha*ValGaussExtremity(i,1)*feval(2*(order+1-j)+1);
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,1)*feval(2*(order+1-j)+1);
	      }
	  }
      }
    else
      {
	for (int j = 0; j <= order+1; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,j);
		if (j == order+1)
		  res(num_dof) -= alpha*ValGaussExtremity(i,0)*feval(2*(order+1-j));
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,0)*feval(2*(order+1-j));
		
		num_dof = NumDofs_Y(0, i);
                if (j == 0)
                  res(num_dof) -= alpha*ValGaussExtremity(i,0)*feval(2*(order+1-j)+1);
                else if (j == order+1)
                  res(num_dof) -= alpha*ValGaussExtremity(i,1)*feval(2*(order+1-j)+1);
                else
                  {
                    if (j == i+1)
                      res(num_dof) -= alpha*feval(2*(order+1-j)+1);
                  }
	      }
	  }
      }
  }
  

  //! Integration against curl of basis functions on an edge
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \nabla \times \varphi_i dx
    \param[in] n edge number
  */  
  template<class T>
  void QuadrangleHcurlOptimalFirstFamily::
  ApplyNablaShGen(const T& alpha, const Vector<T>& feval, Vector<T>& res, int n, int r) const
  {
    for (int i = 0; i <= order+1; i++)
      {
	int j = this->num_quad_points_surf(n)(i);
        for (int k = 0; k < nb_dof_loc; k++)
          res(k) += alpha*Curl_Phi2D(k, j)*feval(i);
      }
  }
  
  
  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlOptimalFirstFamily
  ::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    //typename Vector1::value_type vloc;
    FillZero(Unode);
    
    QuadrangleReference<2>::ComputeNodalValuesRef(Un, Unode);
    return;
  }
  
  
  //! computation of curl u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode curl of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlOptimalFirstFamily
  ::ComputeNodalCurlGen(const Vector1& Un, Vector2 & Ucurl) const
  {
    abort();
    //typedef typename Vector2::value_type Complexe;
    //Vector<TinyVector<Complexe,2> > Unode(points_nodal2d.GetM());
    //Complexe dUy_dx(0), dUx_dy(0);
    
    //abort();
    //QuadrangleReference<2>::ComputeNodalCurlRef(Un, Ucurl);
    // first step, we compute nodal values
    /* ComputeNodalValuesRef(Un, Unode);
    
    // and now computing the curl
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  // du_y/dx
	  dUy_dx = 0;
	  for (int k = 0; k <= order; k++)
	    dUy_dx += dGL_GL(k,i)*Unode(2*NumNodes2D(k,j)+1);
	  
	  // du_x/dy
	  dUx_dy = 0;
	  for (int k = 0; k <= order; k++)
	    dUx_dy += dGL_GL(k,j)*Unode(2*NumNodes2D(i,k));
	  
	  Ucurl(NumNodes2D(i,j)) = dUy_dx - dUx_dy;
	}
    */
  }


  void QuadrangleHcurlOptimalFirstFamily
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    switch (num_loc)
      {
      case 0:
	for (int i = 0; i < order; i++)
	  contrib(i) = feval(2*(i+1));

	break;
      case 1:
	for (int i = 0; i < order; i++)
	  contrib(i) = feval(2*(i+1)+1);

	break;
      case 2:
	for (int i = 0; i < order; i++)
	  contrib(i) = -feval(2*(i+1));

	break;
      case 3:
	for (int i = 0; i < order; i++)
	  contrib(i) = -feval(2*(i+1)+1);

	break;
      }
  }
  

  void QuadrangleHcurlOptimalFirstFamily
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    switch (num_loc)
      {
      case 0:
	for (int i = 0; i < order; i++)
	  contrib(i) = feval(2*(i+1));

	break;
      case 1:
	for (int i = 0; i < order; i++)
	  contrib(i) = feval(2*(i+1)+1);

	break;
      case 2:
	for (int i = 0; i < order; i++)
	  contrib(i) = -feval(2*(i+1));

	break;
      case 3:
	for (int i = 0; i < order; i++)
	  contrib(i) = -feval(2*(i+1)+1);

	break;
      }
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void QuadrangleHcurlOptimalFirstFamily
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
	  int num_dof = NumDofs_X(i, j);
	  res(num_dof)(0) = phiG_x(i)*phiGL_y(j);
	  if (j == order+1)
	    res(num_dof)(0) = -res(num_dof)(0);
	  
	  num_dof = NumDofs_Y(j, i);
	  res(num_dof)(1) = phiGL_x(j)*phiG_y(i);
	  if (j == 0)
	    res(num_dof)(1) = -res(num_dof)(1);
	}
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void QuadrangleHcurlOptimalFirstFamily
  ::ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
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
	  int num_dof = NumDofs_X(i, j);
	  res(num_dof) = -dphiGL_y(j)*phiG_x(i);
	  if (j == order+1)
	    res(num_dof) = -res(num_dof);
	  
	  num_dof = NumDofs_Y(j, i);
	  res(num_dof) = dphiGL_x(j)*phiG_y(i);
	  if (j==0)
	    res(num_dof) = -res(num_dof);
      }
  }
  
  
  //! displays informations about class QuadrangleHcurlFirstFamily
  ostream& operator <<(ostream& out, const QuadrangleHcurlOptimalFirstFamily& e)
  {
    out<<static_cast<const QuadrangleReference<2>&>(e);    
    return out;
  }

}
  
#define MONTJOIE_FILE_QUADRANGLE_HCURL_OPTIMAL_FIRST_FAMILY_CXX
#endif
