#ifndef MONTJOIE_FILE_QUADRANGLE_HCURL_LOBATTO_CXX

namespace Montjoie
{
  
  //! default constructor
  QuadrangleHcurlLobatto::QuadrangleHcurlLobatto() : QuadrangleReference<2>()
  {
    nb_points_quadrature_gauss = 0;
    mass_lumping = true;
  }
  
  
  //! how to number mesh
  void QuadrangleHcurlLobatto::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return QuadrangleReference<2>::ConstructNumberMap(nmap, dg);

    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, order+1);
    nmap.SetNbDofQuadrangle(order, 2*(order+1)*(order-1));
    
    nmap.SetOppositeEdgesDofSymmetry(order, order+1);
    nmap.SetAllEdgesDofToSkewSymmetric(order);
  }
  

  size_t QuadrangleHcurlLobatto::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<2>::GetMemorySize();
    taille += NumDofs_X.GetMemorySize();
    taille += NumDofs_Y.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += Node_to_DofX.GetMemorySize();
    taille += Node_to_DofY.GetMemorySize();
    taille += Sign_CrossMassMatrix.GetMemorySize();
    taille += WeightGauss.GetMemorySize();
    taille += GL_Gexact.GetMemorySize();
    return taille;
  }
  
      
  //! constructing finite element
  void QuadrangleHcurlLobatto::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						      int rsurf, int type_surf)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    
    QuadrangleReference<2>::ConstructFiniteElement(r, rgeom, r, type_quad);    

    ConstructFunctions();
    
    this->Fb_geom.ComputeCoefficientTransformation();
    ConstructHcurlElementaryMatrix();
    
    ConstructStiffnessMatrix();

  }


  //! constructs 1-D finite element (restriction of the 2-D finite element to an edge)
  void QuadrangleHcurlLobatto
  ::ConstructFiniteElement1D(int r, int rgeom, int rquad, int type_quad)
  {  
    EdgeLobattoReference* edge = new EdgeLobattoReference();
    edge->ConstructFiniteElement(order);
    this->element_surface = edge;
  }


  //! construction of basis functions
  void QuadrangleHcurlLobatto::ConstructFunctions()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    lob_basis = lob_quad;
    
    VectReal_wp points_dof1d; VectR2 points_dof2d;
    points_dof1d = this->Points1D();
    points_dof2d = this->PointsND();
    EdgesDof.Reallocate(order+1, 4);
    NumDofs_X.Reallocate(order+1, order+1);
    NumDofs_Y.Reallocate(order+1, order+1);
    nb_dof_loc = 2*(order+1)*(order+1);
    
    // first, tangential dofs on edges are numbered
    for (int i = 0; i < order+1; i++)
      {
	EdgesDof(i, 0) = i;
	NumDofs_X(i, 0) = i;
      }
    
    for (int i = 0; i < order+1; i++)
      {
	EdgesDof(i, 1) = i+order+1;
	NumDofs_Y(order, i) = i+order+1;	
      }
    
    for (int i = 0; i < order+1; i++)
      {
	EdgesDof(i, 2) = i+2*(order+1);
	NumDofs_X(order-i, order) = EdgesDof(i, 2);
      }
    
    for (int i = 0; i < order+1; i++)
      {
	EdgesDof(i, 3) = i+3*(order+1);
	NumDofs_Y(0, order-i) = EdgesDof(i, 3);
      }

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    int numero = 4*(order+1);
    nb_dof_boundaries = numero;
    
    // internal dofs
    for (int j = 1; j < order; j++)
      for (int i = 0; i < order+1; i++)
	{
	  NumDofs_X(i, j) = numero;
	  NumDofs_Y(j, i) = (order+1)*(order-1) + numero++;
	}
    
    CoordinateDofs.Reallocate(nb_dof_loc, 3);
    Node_to_DofX.Reallocate(points_dof2d.GetM());
    Node_to_DofY.Reallocate(points_dof2d.GetM());
    Sign_CrossMassMatrix.Reallocate(points_dof2d.GetM());
    for (int i = 0; i < order+1; i++)
      for (int j = 0; j < order+1; j++)
	{
	  // the two first indices are coordinates in space
	  // the third index is 0 if it is oriented by ex, 1 for ey
	  CoordinateDofs(NumDofs_X(i,j), 0) = i;
	  CoordinateDofs(NumDofs_X(i,j), 1) = j;
	  CoordinateDofs(NumDofs_X(i,j), 2) = 0;
	  
	  CoordinateDofs(NumDofs_Y(j,i), 0) = j;
	  CoordinateDofs(NumDofs_Y(j,i), 1) = i;
	  CoordinateDofs(NumDofs_Y(j,i), 2) = 1;
	  
	  Node_to_DofX(NumQuad2D(i,j)) = NumDofs_X(i,j);
	  Node_to_DofY(NumQuad2D(i,j)) = NumDofs_Y(i,j);
	  if ((i==0)^(j==order))
	    Sign_CrossMassMatrix(NumQuad2D(i,j)) = -1;
	  else
	    Sign_CrossMassMatrix(NumQuad2D(i,j)) = 1;
	}
    
    nb_points_quadrature_gauss = points_dof2d.GetM();
    WeightGauss.Reallocate(points_dof2d.GetM());
    for (int i = 0; i < points_dof2d.GetM(); i++)
      WeightGauss(i) = this->WeightsND(i);
    
    this->SetPointsDof1D(points_dof1d);
    this->SetPointsDofND(points_dof2d);

    this->nb_points_dof_inside = points_dof2d.GetM();
    this->num_dof_points_surf = this->num_quad_points_surf;
  }
  
  
  //! computation of stiffness matrix
  void QuadrangleHcurlLobatto::ConstructStiffnessMatrix()
  {
    int nb_points_nodal_elt = this->PointsDofND().GetM();
    // matrix R_h and transpose
    Matrix<Real_wp, General, ArrayRowSparse> Rh_curl(nb_dof_loc, nb_points_nodal_elt);
    
    // loop over nodes
    int numpsi; Real_wp vloc;
    for (int i = 0; i < order+1; i++)
      {
	for (int j = 0; j < order+1; j++)
	  {
	    // R_h {(k,l,1),(i,j)} = \int_\hat{K} rot( \phi{k,l} e_x) \phi{i,j}
	    //                     = \sum_{m,n} -\omega_{m,n} \phi_l^'(\xi_n)
	    //                       \phi_k(\xi_m) \phi_i(\xi_m) \phi_j(\xi_n)
	    // where \omega_{m,n} is the 2-D Gauss-Lobatto weight
	    // and \xi_m the Gauss-Lobatto point
	    // \phi_i is the basis function related to the point i
	    
	    // Finally, we get
	    // R_h {(k,l,1),(i,j)} = -\omega_{i,j} * \phi_l^'(\xi_j) \delta_{i,k}
	    // the scalar dofs (i,j) interact with vectorial dofs 
	    // oriented along e_x located on the same column (i=k)
	    
	    // R_h {(k,l,2),(i,j)} = \omega_{i,j} * \phi_k^'(\xi_i) \delta_{j,l}
	    // the scalar dofs (i,j) interact with vectorial dofs 
	    // oriented along e_y located on the same row (j=l)
	    
	    // node number
	    int num_loc = NumQuad2D(i,j);
	    
	    // interaction between vectorial dofs ex and scalar dofs
	    for (int l = 0; l < order+1; l++)
	      {
		//numero = NumQuad2D(i,l);
		numpsi = NumDofs_X(i,l); // vectorial dof number
		
		vloc = -lob_basis.GradPhi(l,j);
		if (l == (order+1-1))
		  vloc = -vloc;
		
		// we don't multiply by the weight
		// since this weight is rejected in mass matrix
		if (abs(vloc) >= epsilon_machine)
		  {
		    Rh_curl.AddInteraction(numpsi, num_loc, vloc);
                  }
	      }
	    
	    // interaction between vectorial dofs ey and scalar dofs
	    for (int k = 0; k < order+1; k++)
	      {
		//numero = NumQuad2D(k,j);
		numpsi = NumDofs_Y(k,j); // vectorial dof number
		
		vloc = lob_basis.GradPhi(k,i);
		if (k==0)
		  vloc = -vloc;
		
		// we don't multiply by the weight
		// since this weight is rejected in mass matrix
		if (abs(vloc) >= epsilon_machine)
		  {
		    Rh_curl.AddInteraction(numpsi, num_loc, vloc);
                  }
	      }
	  }
      }
    
    // conversion to Column Sparse Row format, to be more efficient
    Seldon::Copy(Rh_curl, rh_loc);
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
    GL_Gexact.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	GL_Gexact(i, j) = lob_quad.EvaluatePhi(i, gauss.Points(j));
  }
  
  
  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
  */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlLobatto::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Uh(Node_to_DofX(i)) = Vh(2*i);
	Uh(Node_to_DofY(i)) = Vh(2*i+1);
      }
    
    // opposite signs for horizontal dofs on y = 1, and vertical dofs on x = 0
    for (int i = 0; i <= order; i++)
      {
	Uh(NumDofs_X(i, order)) *= -1;
	Uh(NumDofs_Y(0, i)) *= -1;
      } 
  }
  
  
  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlLobatto::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Vh(2*i) = Uh(Node_to_DofX(i));
	Vh(2*i+1) = Uh(Node_to_DofY(i));
      }
    
    // opposite signs for horizontal dofs on y = 1, and vertical dofs on x = 0
    for (int i = 0; i <= order; i++)
      {
	Vh(2*NumQuad2D(i, order)) *= -1;
	Vh(2*NumQuad2D(0, i)+1) *= -1;
      }
  }
  

  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlLobatto::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    if (order != this->GetGeometryOrder())
      {
        QuadrangleReference<2>::ComputeNodalValuesRef(Un, Unode);
        return;
      }
    
    const Matrix<int>& NumNodes2D = this->GetNumNodes2D();
    // loop on each nodal point
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  Unode(2*NumNodes2D(i,j)) = Un(NumDofs_X(i,j));
	  Unode(2*NumNodes2D(i,j)+1) = Un(NumDofs_Y(i,j));
	  if (j == order)
	    Unode(2*NumNodes2D(i,j)) *= -1;
	  if (i == 0)
	    Unode(2*NumNodes2D(i,j)+1) *= -1;
	}
    
  }
  
  
  //! projection of feval on dofs
  /*!
    \param[in] feval evaluation of a function on dof points
    \param[out] contrib dof components
   */
  template<class Vector1, class Vector2>
  void QuadrangleHcurlLobatto::
  ComputeProjectionDofGen(const Vector1& feval, Vector2& contrib) const
  {
    typedef typename Vector2::value_type Complexe;
    TinyVector<Complexe,2> tmp;
    // E ( \xi_i) = \sum DF_i^{*-1} (E_j^1 \phi_j e_x + E_j^2 phi_j e_y )
    // [E_i^1; E_i^2] = DF_i^* E( \xi_i )
    contrib.Fill(0);
    // return;
    for (int node = 0; node < this->PointsDofND().GetM(); node++)
      {
	int i = CoordinateQuad(node, 0);
	int j = CoordinateQuad(node, 1);
	tmp(0) = feval(2*node); tmp(1) = feval(2*node+1);
	if (j == (order+1-1))
	  contrib(Node_to_DofX(node)) = -tmp(0);
	else
	  contrib(Node_to_DofX(node)) = tmp(0);
	
	if (i == 0)
	  contrib(Node_to_DofY(node)) = -tmp(1);
	else
	  contrib(Node_to_DofY(node)) = tmp(1);
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
  void QuadrangleHcurlLobatto::
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
  

  //! Integration against basis functions on an edge
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc edge number
  */  
  template<class T>
  void QuadrangleHcurlLobatto::
  ApplyShGen(const T& alpha, int num_loc, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    int num_dof;
    if (num_loc == 0)
      {
	for (int i = 0; i < order+1; i++)
	  {
	    num_dof = NumDofs_X(i, 0);
	    res(num_dof) += alpha*feval(2*i);
	    num_dof = NumDofs_Y(i, 0);
	    if (i == 0)
              res(num_dof) -= alpha*feval(2*i+1);
            else
              res(num_dof) += alpha*feval(2*i+1);
	  }
      }
    else if (num_loc == 1)
      {
	for (int i = 0; i < order+1; i++)
	  {
	    num_dof = NumDofs_X(order, i);
	    if (i == order)
              res(num_dof) -= alpha*feval(2*i);
            else
              res(num_dof) += alpha*feval(2*i);
            
	    num_dof = NumDofs_Y(order, i);
	    res(num_dof) += alpha*feval(2*i+1);
	  }
      }
    else if (num_loc == 2)
      {
	for (int i = 0; i < order+1; i++)
	  {
	    num_dof = NumDofs_X(order-i, order);
	    res(num_dof) -= alpha*feval(2*i);
	    
            num_dof = NumDofs_Y(order-i, order);
	    if (i == order)
              res(num_dof) -= alpha*feval(2*i+1);
            else
              res(num_dof) += alpha*feval(2*i+1);
	  }
      }
    else
      {
	for (int i = 0; i < order+1; i++)
	  {
	    num_dof = NumDofs_X(0, order-i);
	    if (i == 0)
              res(num_dof) -= alpha*feval(2*i);
            else
              res(num_dof) += alpha*feval(2*i);
            
	    num_dof = NumDofs_Y(0, order-i);
	    res(num_dof) -= alpha*feval(2*i+1);
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
  void QuadrangleHcurlLobatto::
  ApplyNablaShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    for (int i = 0; i <= order; i++)
      {
	int j = this->num_quad_points_surf(n)(i);
        for (int k = 0; k < nb_dof_loc; k++)
          res(k) += alpha*Curl_Phi2D(k, j)*feval(i);
      }
  }
  

  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void QuadrangleHcurlLobatto::ComputeValuesPhiRef(const R2& point_loc, VectR2& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(R2(0,0));
    VectReal_wp phix(order+1), phiy(order+1);
    for (int i = 0; i <= order; i++)
      {
	phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  res(NumDofs_X(i, j))(0) = phix(i)*phiy(j);
	  res(NumDofs_Y(i, j))(1) = phix(i)*phiy(j);

	  if (j == order)
	    res(NumDofs_X(i, j))(0) = -res(NumDofs_X(i, j))(0);

	  if (i == 0)
	    res(NumDofs_Y(i, j))(1) = -res(NumDofs_Y(i, j))(1);
	}
  }
  
  
  //! Evaluating curl of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res curl of basis functions on point_loc
  */
  void QuadrangleHcurlLobatto::ComputeCurlPhiRef(const R2& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(0);
    VectReal_wp phix(order+1), phiy(order+1), dphix(order+1), dphiy(order+1);
    for (int i = 0; i <= order; i++)
      {
	phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
	dphix(i) = lob_basis.EvaluatePhiGrad(i, point_loc(0));
	dphiy(i) = lob_basis.EvaluatePhiGrad(i, point_loc(1));
      }
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  res(NumDofs_X(i, j)) = -phix(i)*dphiy(j);
	  res(NumDofs_Y(i, j)) = dphix(i)*phiy(j);

	  if (j == order)
	    res(NumDofs_X(i, j)) = -res(NumDofs_X(i, j));

	  if (i == 0)
	    res(NumDofs_Y(i, j)) = -res(NumDofs_Y(i, j));
	}
  }


  //! displays informations about class QuadrangleHcurlLobatto
  ostream& operator <<(ostream& out, const QuadrangleHcurlLobatto& e)
  {
    out<<static_cast<const QuadrangleReference<2>&>(e);
    cout<<"Quadrilateral of Nedelec's second family "<<endl;
    cout<<"Number of local degrees of freedom "<<e.GetNbDof()<<endl;
    
    return out;
  }

}
  
#define MONTJOIE_FILE_QUADRANGLE_HCURL_LOBATTO_CXX
#endif
