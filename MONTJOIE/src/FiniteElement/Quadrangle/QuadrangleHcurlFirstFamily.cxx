#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_FIRST_FAMILY_CXX

namespace Montjoie
{
  //! default constructor
  QuadrangleHcurlFirstFamily::QuadrangleHcurlFirstFamily() : QuadrangleReference<2>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;

    this->optimized_elem_matrix = true;
  }
  
  
  //! how to number mesh
  void QuadrangleHcurlFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<2>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order);
    nmap.SetNbDofQuadrangle(order, 2*order*(order-1));
    
    nmap.SetOppositeEdgesDofSymmetry(order, order);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t QuadrangleHcurlFirstFamily::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<2>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize();
    taille += NumDofs_Y.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += NumNodes2D_gauss.GetMemorySize();
    taille += normale_dof.GetMemorySize();
    taille += dGL_G.GetMemorySize();
    taille += dGL_GL.GetMemorySize();
    taille += G_GL.GetMemorySize();
    taille += GL_Gexact.GetMemorySize();
    taille += GL_G.GetMemorySize();
    taille += ValGaussExtremity.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += ch_loc.GetMemorySize() + WeightsGauss.GetMemorySize();
    return taille;
  }

  
  //! constructing finite element
  void QuadrangleHcurlFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf, int type_surf)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    
    QuadrangleReference<2>::ConstructFiniteElement(r, rgeom, r, type_quad);
    
    lob_basis.ConstructQuadrature(order-1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    int nb_points_quadrature_gauss = order*order;
    VectR2 points2d = this->PointsND();
    points2d.Resize(nb_points_quadrature_inside + nb_points_quadrature_gauss);
    WeightsGauss.Reallocate(nb_points_quadrature_gauss);
    int nb = nb_points_quadrature_inside;
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
        {
          WeightsGauss(nb-nb_points_quadrature_inside) = lob_basis.Weights(i)*lob_basis.Weights(j);
          points2d(nb).Init(lob_basis.Points(i), lob_basis.Points(j));
          nb++;
        }
    
    this->SetPointsND(points2d);
    
    ConstructFunctions();
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructHcurlElementaryMatrix();
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();     

    if (this->element_surface != NULL)
      delete this->element_surface;

    EdgeGaussReference* edge = new EdgeGaussReference();
    edge->ConstructFiniteElement(r-1, r, r, Globatto<Real_wp>::QUADRATURE_LOBATTO,
				 EdgeGauss::GAUSS);

    edge->SetDofPoints(this->PointsDof1D());
    this->element_surface = edge;

  }


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  void QuadrangleHcurlFirstFamily
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {
  }  
  
  
  //! construction of basis functions
  void QuadrangleHcurlFirstFamily::ConstructFunctions()
  {
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
    CoordinateDofs.Reallocate(nb_dof_loc,3);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  points_dof2d(NumDofs_X(i,j)).Init(lob_basis.Points(i), lob_quad.Points(j));
	  points_dof2d(NumDofs_Y(j,i)).Init(lob_quad.Points(j), lob_basis.Points(i));
	  
	  // the two first indices are coordinates in space
	  // the third index is 0 if it is oriented by ex, 1 for ey
	  CoordinateDofs(NumDofs_X(i,j), 0) = i;
	  CoordinateDofs(NumDofs_X(i,j), 1) = j;
	  CoordinateDofs(NumDofs_X(i,j), 2) = 0;
	  CoordinateDofs(NumDofs_Y(j,i), 0) = j;
	  CoordinateDofs(NumDofs_Y(j,i), 1) = i;
	  CoordinateDofs(NumDofs_Y(j,i), 2) = 1;
	}
    
    NumNodes2D_gauss.Reallocate(order, order);
    WeightsGauss.Reallocate(order*order);
    int compt = 0;
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
	{
	  WeightsGauss(compt) = lob_basis.Weights(i)*lob_basis.Weights(j);
	  NumNodes2D_gauss(i,j) = compt++;
	}
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    this->nb_points_dof_inside = points_dof2d.GetM();
    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(order);
	for (int i = 0; i < order; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }
  }
  
  
  //! computation of stiffness matrix
  void QuadrangleHcurlFirstFamily::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    dGL_G.Reallocate(order+1, order);
    for (int j = 0; j <= order; j++)
      for (int i = 0; i < order; i++)
	dGL_G(j,i) = lob_quad.EvaluatePhiGrad(j, lob_basis.Points(i));
    
    dGL_GL = lob_quad.GradPhi();
    
    int nb_points_fct_H1 = (order+1)*(order+1);
    Matrix<Real_wp, General, ArrayRowSparse> Ch(nb_dof_loc, 2*nb_points_fct_H1);    
    Matrix<Real_wp, General, ArrayRowSparse> Rh(nb_dof_loc, order*order);

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k < order; k++)
	  {
	    int node = NumQuad2D(i, j);
	    int num_dof = NumDofs_X(k, j);
	    Real_wp vloc = G_GL(k,i);
	    if (j == order)
	      vloc = -vloc;
	    
	    Ch.AddInteraction(num_dof, 2*node, vloc);
	    
	    vloc = G_GL(k,j);
	    if (i == 0)
	      vloc = -vloc;
	    
	    num_dof = NumDofs_Y(i,k);
	    Ch.AddInteraction(num_dof, 2*node+1, vloc);
          }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumNodes2D_gauss(i,j);
	    int num_dof = NumDofs_X(i,k);
	    Real_wp vloc = -dGL_G(k,j);
	    if (k == order)
	      vloc = -vloc;
	    
	    Rh.AddInteraction(num_dof, node, vloc);
	    
	    num_dof = NumDofs_Y(k,j);
	    vloc = dGL_G(k,i);
	    if (k == 0)
	      vloc = -vloc;
	    
	    Rh.AddInteraction(num_dof, node, vloc);
          }
    
    Seldon::Copy(Rh, rh_loc);
    Seldon::Copy(Ch, ch_loc);

  }
  
  
  //! computation of mass matrix
  void QuadrangleHcurlFirstFamily::ConstructMassMatrix()
  {
    ValGaussExtremity.Reallocate(order, 2);
    for (int i = 0; i < order; i++)
      {
	ValGaussExtremity(i,0) = lob_basis.EvaluatePhi(i, 0.0);
	ValGaussExtremity(i,1) = lob_basis.EvaluatePhi(i, 1.0);
      }
    
    G_GL.Reallocate(order, order+1);
    GL_G.Reallocate(order+1, order);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  G_GL(i,j) = lob_basis.EvaluatePhi(i, lob_quad.Points(j));
	  GL_G(j,i) = lob_quad.EvaluatePhi(j, lob_basis.Points(i));
	}  
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
    GL_Gexact.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	GL_Gexact(i, j) = lob_quad.EvaluatePhi(i, gauss.Points(j));
    
  }
  

  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] res dof components
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlFirstFamily::
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
  
  
  //! projection from quadrature points (Gauss-Lobatto points) to dof components
  /*!
    \param[in] Equad values on quadrature points
    \param[out] Edof dof components
  */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlFirstFamily::
  ProjectQuadratureToDofGen(const Vector1& Equad, Vector2& Edof) const
  {
    FillZero(Edof);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  // dofs along x
	  int num_dof = NumDofs_X(i,j);
	  if (j == order)
            for (int k = 0; k <= order; k++)
              Edof(num_dof) -= GL_G(k,i)*Equad(2*NumQuad2D(k,j));
          else
            for (int k = 0; k <= order; k++)
              Edof(num_dof) += GL_G(k,i)*Equad(2*NumQuad2D(k,j));
	  
	  // dofs along y
	  num_dof = NumDofs_Y(j,i);
	  if (j == 0)
            for (int k = 0; k <= order; k++)
              Edof(num_dof) -= GL_G(k,i)*Equad(2*NumQuad2D(j,k)+1);
          else
            for (int k = 0; k <= order; k++)
              Edof(num_dof) += GL_G(k,i)*Equad(2*NumQuad2D(j,k)+1);
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
  template<class T>
  void QuadrangleHcurlFirstFamily::
  ApplyShGen(const T& alpha, int num_loc, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    int num_dof;
    if (num_loc == 0)
      {
	for (int j = 0; j <= order; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,0);
		res(num_dof) += alpha*G_GL(i,j)*feval(2*j);
                
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
	for (int j = 0; j <= order; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,j);
		if (j == order)
		  res(num_dof) -= alpha*ValGaussExtremity(i,1)*feval(2*j);
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,1)*feval(2*j);
		
                num_dof = NumDofs_Y(order,i);
		res(num_dof) += alpha*G_GL(i,j)*feval(2*j+1);
	      }
	  }
      }
    else if (num_loc == 2)
      {
	for (int j = 0; j <= order; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,order);
		res(num_dof) -= alpha*G_GL(i,j)*feval(2*(order-j));
		
                num_dof = NumDofs_Y(j,i);
		if (j == 0)
		  res(num_dof) -= alpha*ValGaussExtremity(i,1)*feval(2*(order-j)+1);
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,1)*feval(2*(order-j)+1);
	      }
	  }
      }
    else
      {
	for (int j = 0; j <= order; j++)
	  {
	    for (int i = 0; i < order; i++)
	      {
		num_dof = NumDofs_X(i,j);
		if (j == order)
		  res(num_dof) -= alpha*ValGaussExtremity(i,0)*feval(2*(order-j));
		else
		  res(num_dof) += alpha*ValGaussExtremity(i,0)*feval(2*(order-j));
		
		num_dof = NumDofs_Y(0, i);
		res(num_dof) -= alpha*G_GL(i,j)*feval(2*(order-j)+1);
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
  void QuadrangleHcurlFirstFamily::
  ApplyNablaShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    for (int i = 0; i <= order; i++)
      {
	int j = this->num_quad_points_surf(n)(i);
        for (int k = 0; k < nb_dof_loc; k++)
          res(k) += alpha*Curl_Phi2D(k, j)*feval(i);
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
  void QuadrangleHcurlFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    Vector1 feval2(2*(order+1));
    feval2.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        {
          feval2(2*i) += GL_Gexact(i, j)*feval(2*j);
          feval2(2*i+1) += GL_Gexact(i, j)*feval(2*j+1);
        }
    
    ComputeIntegralSurfaceRef(feval2, res, num_loc);
  }  
  
  
  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlFirstFamily::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    typename Vector1::value_type vloc;
    FillZero(Unode);
    
    const Matrix<int>& NumNodes2D = this->GetNumNodes2D();
    if (order != this->GetGeometryOrder())
      {
        QuadrangleReference<2>::ComputeNodalValuesRef(Un, Unode);
        return;
      }
    
    // loop on each nodal point
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  // dofs along x
	  vloc = 0;
	  for (int k = 0; k < order; k++)
	    vloc += G_GL(k,i)*Un(NumDofs_X(k,j));
	  
	  if (j == order)
	    vloc = -vloc;
	  
	  Unode(2*NumNodes2D(i,j)) = vloc;
	  
	  // dofs along y
	  vloc = 0;
	  for (int k = 0; k < order; k++)
	    vloc += G_GL(k,j)*Un(NumDofs_Y(i,k));
	  
	  if (i == 0)
	    vloc = -vloc;
	  
	  Unode(2*NumNodes2D(i,j)+1) = vloc;
	}
    
  }
  

  //! computation of curl u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode curl of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlFirstFamily::ComputeNodalCurlGen(const Vector1& Un, Vector2 & Ucurl) const
  {
    typedef typename Vector2::value_type Complexe;
    Vector<Complexe> Unode(2*this->GetNbPointsNodalElt());
    Complexe dUy_dx(0), dUx_dy(0);

    const Matrix<int>& NumNodes2D = this->GetNumNodes2D();    
    // first step, we compute nodal values
    ComputeNodalValuesRef(Un, Unode);
    
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
  }
  
  
  //! adds \int B varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] B tensor B on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int B varphi_j  varphi_i
  */
  template<class T, class Prop>
  void QuadrangleHcurlFirstFamily::
  AddVariableMassMatrixGen(int off_row, int off_col, const Vector<TinyMatrix<T, Prop, 2, 2> >& B,
			   VirtualMatrix<T>& mat) const
  {
    int num_dof_x, num_dof_y, num_dof2_x, num_dof2_y;
    T vloc_x, vloc_y;

    const Matrix<int>& NumNodes2D = this->NumQuad2D;
    // loop over unknowns
    bool sym = mat.IsSymmetric();
    
    // loop on dofs
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
        {
          // dofs oriented along e_x, e_y
          num_dof_x = NumDofs_X(i, j) + off_row;
          num_dof_y = NumDofs_Y(j, i) + off_row;
          
          // second loop on dofs
          for (int k = 0; k < order; k++)
            for (int l = 0; l <= order; l++)
              {
                num_dof2_x = NumDofs_X(k, l) + off_col;
                num_dof2_y = NumDofs_Y(l, k) + off_col;
                vloc_x = 0.0; vloc_y = 0.0;
                
                // interaction between num_dof_x and num_dof2_x
                // interaction between num_dof_y and num_dof2_y
                if (j==l)
                  {
                    for (int n = 0; n <= order; n++)
                      {
                        vloc_x += B(NumNodes2D(n, j))(0,0)*G_GL(i, n)*G_GL(k, n);
                        vloc_y += B(NumNodes2D(j, n))(1, 1)*G_GL(i, n)*G_GL(k, n);
                      }
                  }
                
                // now signs of local dofs are considered
                // negative for dofs oriented along e_y with i = 0
                // negative for dofs oriented along e_x with j = r
                if ((num_dof_x <= num_dof2_x) || (!sym) )
                  {
                    if (((j==order)&&(l!=order))||((j!=order)&&(l==order)))
                      vloc_x = -vloc_x;
                    
                    mat.AddInteraction(num_dof_x, num_dof2_x, vloc_x);
                  }
                
                if ((num_dof_y<=num_dof2_y) || (!sym) )
                  {
                    if (((j==0)&&(l!=0))||((j!=0)&&(l==0)))
                      vloc_y = -vloc_y;
                    
                    mat.AddInteraction(num_dof_y, num_dof2_y, vloc_y);
                  }
                
                // interaction between num_dof_x and num_dof2_y
                // interaction between num_dof_y and num_dof2_x
                if ((num_dof_x <= num_dof2_y) || (!sym) )
                  {
                    vloc_x = B(NumNodes2D(l, j))(0,1)*G_GL(i, l)*G_GL(k, j);
                    
                    if (((j==order)&&(l!=0))||((j!=order)&&(l==0)))
                      vloc_x=-vloc_x;
                    
                    mat.AddInteraction(num_dof_x, num_dof2_y, vloc_x);
                  }
                
                if ((num_dof_y<=num_dof2_x) || (!sym) )
                  {
                    vloc_y = B(NumNodes2D(j, l))(1,0)*G_GL(i, l)*G_GL(k, j);
                    
                    if (((j==0)&&(l!=order))||((j!=0)&&(l==order)))
                      vloc_y=-vloc_y;
                    
                    mat.AddInteraction(num_dof_y, num_dof2_x, vloc_y);
                  }
              }
        }
  }


  //! adds \int A \nabla \times varphi_j \nabla \times varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A \nabla \times varphi_j \nabla \times varphi_i
  */
  template<class T>
  void QuadrangleHcurlFirstFamily::
  AddVariableStiffnessMatrixGen(int off_row, int off_col, const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    int num_dof_x, num_dof_y, num_dof2_x, num_dof2_y;
    T vloc_x, vloc_y;
    
    // loop over unknowns
    bool sym = mat.IsSymmetric();
    
    // loop on dofs
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
        {
          // dofs oriented along e_x, e_y
          num_dof_x = NumDofs_X(i, j) + off_row;
          num_dof_y = NumDofs_Y(j, i) + off_row;
          
          // second loop on dofs
          for (int k = 0; k < order; k++)
            for (int l = 0; l <= order; l++)
              {
                num_dof2_x = NumDofs_X(k, l) + off_col;
                num_dof2_y = NumDofs_Y(l, k) + off_col;
                vloc_x = 0.0; vloc_y = 0.0;
                
                // interaction between num_dof_x and num_dof2_x
                // interaction between num_dof_y and num_dof2_y
                if (i==k)
                  {
                    for (int n = 0; n < order; n++)
                      {
                        vloc_x += A(NumNodes2D_gauss(i, n))*dGL_G(j, n)*dGL_G(l, n);
                        vloc_y += A(NumNodes2D_gauss(n, i))*dGL_G(j, n)*dGL_G(l, n);
                      }
                  }
                
                // now signs of local dofs are considered
                // negative for dofs oriented along e_y with i = 0
                // negative for dofs oriented along e_x with j = r
                if ((num_dof_x <= num_dof2_x) || (!sym) )
                  {
                    if (((j==order)&&(l!=order))||((j!=order)&&(l==order)))
                      vloc_x = -vloc_x;
                    
                    mat.AddInteraction(num_dof_x, num_dof2_x, vloc_x);
                  }
                
                if ((num_dof_y<=num_dof2_y) || (!sym) )
                  {
                    if (((j==0)&&(l!=0))||((j!=0)&&(l==0)))
                      vloc_y = -vloc_y;
                    
                    mat.AddInteraction(num_dof_y, num_dof2_y, vloc_y);
                  }
                
                // interaction between num_dof_x and num_dof2_y
                // interaction between num_dof_y and num_dof2_x
                if ((num_dof_x <= num_dof2_y) || (!sym) )
                  {
                    vloc_x = -A(NumNodes2D_gauss(i,k))*dGL_G(j, k)*dGL_G(l, i);
                    
                    if (((j==order)&&(l!=0))||((j!=order)&&(l==0)))
                      vloc_x=-vloc_x;
                    
                    mat.AddInteraction(num_dof_x, num_dof2_y, vloc_x);
                  }
                
                if ((num_dof_y<=num_dof2_x) || (!sym) )
                  {
                    vloc_y = -A(NumNodes2D_gauss(k,i))*dGL_G(j, k)*dGL_G(l, i);
                    
                    if (((j==0)&&(l!=order))||((j!=0)&&(l==order)))
                      vloc_y=-vloc_y;
                    
                    mat.AddInteraction(num_dof_y, num_dof2_x, vloc_y);
                  }
              }
        }
  }

      
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void QuadrangleHcurlFirstFamily::ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
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
	phiGL_x(i) = lob_quad.EvaluatePhi(i, point_loc(0));
	phiGL_y(i) = lob_quad.EvaluatePhi(i, point_loc(1));
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
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void QuadrangleHcurlFirstFamily::ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
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
	dphiGL_x(i) = lob_quad.EvaluatePhiGrad(i, point_loc(0));
	dphiGL_y(i) = lob_quad.EvaluatePhiGrad(i, point_loc(1));
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
  
  
  //! Evaluation of \f$ \psi_i(point\_loc) \f$
  /*!
    where \f$ \psi_i \f$  are basis functions used for H
    H is discretized on r^2 Gauss points
  */
  void QuadrangleHcurlFirstFamily::
  ComputeValuesPhiH(const R2& point_loc, const R2& point_glob,
		    Vector<TinyVector<Real_wp, 1> >& res, const Matrix2_2& dfjm1,
		    const Mesh<Dimension2>& mesh, int nquad) const
  {
    res.Reallocate(order*order);
    VectReal_wp phix(order), phiy(order); phix.Fill(0); phiy.Fill(0);
    for (int i = 0; i < order; i++)
      {
	phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
	{
	  int num_dof = NumNodes2D_gauss(i,j);
	  res(num_dof)(0) = phix(i)*phiy(j);
	}
  }

  
  //! projection from nodal points (Gauss-Lobatto) to Gauss points (used for H)
  template<class Vector1, class Vector2>
  void QuadrangleHcurlFirstFamily::ComputeNodalValuesH(const Vector1& Un, Vector2 & Unode) const
  {
    int r = this->order;
    typedef typename Vector2::value_type Complexe;
    Vector<Complexe> Uint(r*(r+1));
    // computing U on G - GL points
    for (int i = 0; i < r; i++)
      for (int j = 0; j <= r; j++)
	{
	  int node = j*r + i;
	  Uint(node) = 0;
	  for (int k = 0; k < r; k++)
	    Uint(node) += G_GL(k, j)*Un(NumNodes2D_gauss(i,k));
	}

    const Matrix<int>& NumNodes2D = this->GetNumNodes2D();    
    // computation on GL -GL now
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  int node = NumNodes2D(i,j);
	  Unode(node) = 0;
	  for (int k = 0; k < r; k++)
	    Unode(node) += G_GL(k,i)*Uint(j*r+k);
	}
  }

  
  //! projection from Gauss points to nodal points
  template<class Vector1>
  void QuadrangleHcurlFirstFamily::
  InterpolateHQuadrature(const SetPoints<Dimension2>& PointsElem,
                         const SetMatrices<Dimension2>& MatricesElem,
                         const Vector1& Hloc, Vector1& Hloc_node,
                         const ElementReference<Dimension2, 2>& FaceBasis,
                         const Mesh<Dimension2>& mesh, int nquad) const
  {
    const Matrix<int>& NumNodes2D = this->GetNumNodes2D();
    Vector1 Hint(order*(order+1)); Hint.Fill(0);
    for (int i = 0; i < order; i++)
      for (int j = 0; j <= order; j++)
	{
	  for (int k = 0; k < order; k++)
	    Hint(i*(order+1)+j) += G_GL(k,j)*Hloc(NumNodes2D_gauss(i,k));
	}
    
    Hloc_node.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  for (int k = 0; k < order; k++)
	    Hloc_node(NumNodes2D(i,j)) += G_GL(k,i)*Hint(k*(order+1)+j);
	}
    
  }
    

  //! displays informations about class QuadrangleHcurlFirstFamily
  ostream& operator <<(ostream& out, const QuadrangleHcurlFirstFamily& e)
  {
    out<<static_cast<const QuadrangleReference<2>&>(e);    
    return out;
  }

}
  
#define MONTJOIE_FILE_QUADRANGLE_HCURL_FIRST_FAMILY_CXX
#endif
