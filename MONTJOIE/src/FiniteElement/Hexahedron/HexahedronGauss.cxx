#ifndef MONTJOIE_FILE_HEXAHEDRON_GAUSS_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronGauss::HexahedronGauss() : HexahedronReference<1>()
  {
    Fb_geom.quadrature_equal_nodal = false;
    Fb_geom.dof_equal_nodal = false;
    Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
    use_quadrature_for_sh = false;
  }
  

  size_t HexahedronGauss::GetMemorySize() const
  {
    size_t taille = HexahedronReference<1>::GetMemorySize();
    taille += dGL_G.GetMemorySize() + GL_G.GetMemorySize();
    taille += invGL_G.GetMemorySize() + Gquad_GL.GetMemorySize();
    taille += rh_loc.GetMemorySize() + rh_locX.GetMemorySize() +
      rh_locY.GetMemorySize() + rh_locZ.GetMemorySize();
    taille += rh_weight.GetMemorySize() + rh_weightX.GetMemorySize() +
      rh_weightY.GetMemorySize() + rh_weightZ.GetMemorySize();
    taille += ch1_node.GetMemorySize() + ch2_node.GetMemorySize() + ch3_node.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize();
    taille += NumDofs3D.GetMemorySize() + CoordinateDofs.GetMemorySize() +
      NumDofs2D.GetMemorySize() + CoordinateDofs_quad.GetMemorySize();
    taille += GL_Geom.GetMemorySize() + dGL_Geom.GetMemorySize();
    taille += invSqrtWeights3d.GetMemorySize() + invWeights3d.GetMemorySize();
    taille += stiffness_matrix1D.GetMemorySize() + mass_matrix1D.GetMemorySize()
      + gradient_matrix1D.GetMemorySize();

    taille += Seldon::GetMemorySize(NumProjOperator);
    return taille;
  }

  
  //! construction of finite element
  void HexahedronGauss::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					       int rsurf_tri, int rsurf_quad,
					       int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (type_quad == -1)
      type_quad = Globatto<Real_wp>::QUADRATURE_GAUSS;
    
    HexahedronReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    invWeights3d.Reallocate((order_quad+1)*(order_quad+1)*(order_quad+1));
    invSqrtWeights3d.Reallocate((order_quad+1)*(order_quad+1)*(order_quad+1));
    for (int k = 0; k < invWeights3d.GetM(); k++)
      {
        invWeights3d(k) = 1.0/this->WeightsND(k);
        invSqrtWeights3d(k) = 1.0/sqrt(this->WeightsND(k));
      }
    
    ConstructFunctions();
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();

    if (type_quad != Globatto<Real_wp>::QUADRATURE_LOBATTO)
      {
	//ConstructElementaryMatrix(*this);
	//Value_Phi.Clear(); Gradient_Phi.Clear();
	//mass_matrix.Clear(); mass_matrix_chol.Clear();
	//ValuePhi_Nodal.Clear();
      }

    int order_geom = rgeom;
    if (order_geom != order)
      {
        Matrix<Real_wp> Phi1D;
        Phi1D.Reallocate(order+1, order_geom+1);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order_geom; j++)
            Phi1D(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
        
        Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
        Ch1.Clear();
        Ch1.Reallocate((order+1)*(order+1)*(order_geom+1), (order+1)*(order+1)*(order+1));
        
        Ch2.Clear();
        Ch2.Reallocate((order+1)*(order_geom+1)*(order_geom+1),
                       (order+1)*(order_geom+1)*(order+1));

        Ch3.Clear();
        Ch3.Reallocate((order_geom+1)*(order_geom+1)*(order_geom+1),
                       (order_geom+1)*(order_geom+1)*(order+1));

        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            for (int k = 0; k <= order_geom; k++)
              for (int m = 0; m <= order; m++)
                Ch3.AddInteraction(Fb_geom.GetNumNodes3D(i, j, k),
                                   (order_geom+1)*(order+1)*i + (order+1)*j + m, Phi1D(m, k));
        
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            for (int k = 0; k <= order; k++)
              for (int m = 0; m <= order; m++)
                Ch2.AddInteraction((order_geom+1)*(order+1)*i + (order+1)*j + k,
                                   (order+1)*(order+1)*i + (order+1)*m + k, Phi1D(m, j));
        
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              for (int m = 0; m <= order; m++)
                Ch1.AddInteraction((order+1)*(order+1)*i + (order+1)*j + k,
                                   NumDofs3D(m, j, k), Phi1D(m, i));
        
        Copy(Ch1, ch1_node); Copy(Ch2, ch2_node); Copy(Ch3, ch3_node);
      }

    function_basis_quad = new QuadrangleGauss();
    function_basis_quad->ConstructFiniteElement(r, rgeom, rquad, type_quad);
    element_quad_surf = function_basis_quad;
  }
  
  
  //! construction of basis functions
  void HexahedronGauss::ConstructFunctions()
  {    
    int r = order;
    int order_geom = this->GetGeometryOrder();
    if (order == order_geom)
      {
	this->elt_geom.dof_equal_nodal = true;
        nb_dof_loc = this->GetNbPointsNodalElt();
        nb_dof_boundaries = 6*r*r + 2;
        
        lob_basis = this->GetNodalShapeFunctions1D();
        
        FacesDof = Fb_geom.GetNodalNumber();
        NumDofs3D = Fb_geom.GetNumNodes3D();
        CoordinateDofs = Fb_geom.GetCoordinateNodes3D();
        NumDofs2D = Fb_geom.GetNumNodes2D_quad();
        
        this->SetPointsDof1D(Fb_geom.PointsNodal1D());
        this->SetPointsDof2D_quad(this->PointsNodal2D_quad()); 
        this->SetPointsDofND(Fb_geom.PointsNodalND());
      }
    else
      {
	this->elt_geom.dof_equal_nodal = false;
        MeshNumbering<Dimension2>::
	  ConstructQuadrilateralNumbering(r, NumDofs2D, CoordinateDofs_quad);
        
	MeshNumbering<Dimension3>::
	  ConstructHexahedralNumbering(r, NumDofs3D, CoordinateDofs);
        
        lob_basis.ConstructQuadrature(r, lob_basis.QUADRATURE_LOBATTO);
        lob_basis.ComputeGradPhi(1e3*epsilon_machine);
        
        nb_dof_loc = (r+1)*(r+1)*(r+1);
        nb_dof_boundaries = 6*r*r + 2;
        
        this->SetPointsDof1D(lob_basis.Points());
	VectR2 points_dof2d_quad;
	VectR3 points_dof3d;
        points_dof2d_quad.Reallocate((r+1)*(r+1));
        points_dof3d.Reallocate((r+1)*(r+1)*(r+1));
    
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            points_dof2d_quad(NumDofs2D(i,j)) = R2(this->PointsDof1D(i), this->PointsDof1D(j));

        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            for (int k = 0; k <= r; k++)
              points_dof3d(NumDofs3D(i, j, k)) = R3(this->PointsDof1D(i),
                                                    this->PointsDof1D(j), this->PointsDof1D(k));

        this->SetPointsDof2D_quad(points_dof2d_quad); 
        this->SetPointsDofND(points_dof3d);
        
        FacesDof.Reallocate((r+1)*(r+1), 6);
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D(i, j), 0) = NumDofs3D(0, i, j);
        
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D(i, j), 1) = NumDofs3D(i, 0, j);
        
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D(i, j), 2) = NumDofs3D(i, j, 0);
    
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D(i, j), 3) = NumDofs3D(i, j, r);
 
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D(i, j), 4) = NumDofs3D(i, r, j);
    
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D(i, j), 5) = NumDofs3D(r, i, j);
      }

    this->nb_points_dof_inside = this->nb_dof_loc;
    this->num_dof_points_surf.Reallocate(6);
    for (int n = 0; n < 6; n++)
      {
	this->num_dof_points_surf(n).Reallocate((r+1)*(r+1));
	for (int i = 0; i < FacesDof.GetM(); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }

    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order_quad, gauss.QUADRATURE_GAUSS);
    GL_G.Reallocate(order+1, order_quad+1);
    Gquad_GL.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        {
          GL_G(i,j) = lob_basis.EvaluatePhi(i, gauss.Points(j));
          Gquad_GL(i, j) = gauss.EvaluatePhi(i, lob_basis.Points(j));
        }
    
    invGL_G = GL_G; GetInverse(invGL_G);


    GL_Geom.Reallocate(order + 1, order_geom + 1);
    dGL_Geom.Reallocate(order + 1, order_geom + 1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        {
          GL_Geom(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
          dGL_Geom(i, j) = lob_basis.EvaluatePhiGrad(i,this->PointsNodal1D(j));
        }
    
    Fb_geom.ComputeCoefficientTransformation();    
  }
  

  //! construction of stiffness matrix
  void HexahedronGauss::ConstructStiffnessMatrix()
  { 
    Matrix<Real_wp> Phi1D, dPhi1D;
    Phi1D.Reallocate(order+1, order_quad+1);
    dPhi1D.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int k = 0; k <= order_quad; k++)
        {
          Phi1D(i, k) = lob_basis.EvaluatePhi(i, this->Points1D(k));
          dPhi1D(i, k) = lob_basis.EvaluatePhiGrad(i, this->Points1D(k));
        }
    
    // 1-D stiffness matrix
    stiffness_matrix1D.Reallocate(order+1, order+1);
    gradient_matrix1D.Reallocate(order+1, order+1);
    mass_matrix1D.Reallocate(order+1, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  stiffness_matrix1D(i,j) = 0;
          gradient_matrix1D(i,j) = 0;
          mass_matrix1D(i,j) = 0;
	  for (int k = 0; k <= order_quad; k++)
	    {
              stiffness_matrix1D(i,j) += this->Weights1D(k)*dPhi1D(i,k)*dPhi1D(j,k);
              gradient_matrix1D(i,j) += this->Weights1D(k)*dPhi1D(i,k)*Phi1D(j,k);
              mass_matrix1D(i,j) += this->Weights1D(k)*Phi1D(i,k)*Phi1D(j,k);
            }
	}
   
    int numero,numpsi; Real_wp vloc;
    // matrix Rh and transpose
    Matrix<Real_wp, General, ArrayRowSparse>
      Rh(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_weight(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_x(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_y(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_z(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_weight_x(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_weight_y(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      Rh_weight_z(nb_points_quadrature_inside, 3*nb_points_quadrature_inside),
      RhB(nb_points_quadrature_inside, 3*this->nb_dof_boundaries);
    
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    // loop over quadrature points
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
	for (int k = 0; k <= order_quad; k++)
	  {
	    int num_loc = NumQuad3D(i,j,k); // scalar dof number
	    int l,m,n;
	    m = j; n = k;
	    
	    // part dU/dx v_x
	    for (l = 0; l <= order_quad; l++)
	      {
		numero = NumQuad3D(l, m, n);
		numpsi = 3*numero; // vectorial dof number 
		
		// entry = dphi_i/dx(\xi_l)
		vloc = lob_quad.GradPhi(i, l);
		if (abs(vloc) >= 1e4*epsilon_machine)
		  {
		    // dphi_{i,j,k} / dx (\xi_l,\xi_m,\xi_n) = \phi'_i(\xi_l) 
		    Rh.AddInteraction(num_loc, numpsi, vloc);
                    Rh_weight.AddInteraction(num_loc, numpsi, vloc*this->WeightsND(numero));
		    Rh_x.AddInteraction(num_loc, numpsi, vloc);
                    Rh_weight_x.AddInteraction(num_loc, numpsi, vloc*this->WeightsND(numero));
		    if (numpsi < 3*this->nb_dof_boundaries)
		      RhB.AddInteraction(num_loc, numpsi, vloc);
		  }
	      }
	    
	    l = i; n = k;
	    // part dU/dy v_y
	    for (m = 0; m <= order_quad; m++)
	      {
		numero = NumQuad3D(l,m,n);
		numpsi = 3*numero+1; // vectorial dof number
		
		vloc = lob_quad.GradPhi(j, m);
		if (abs(vloc) >= 1e4*epsilon_machine)
		  {
		    // dphi_{i,j,k} / dy (\xi_l,\xi_m,\xi_n) = \phi'_j(\xi_m) 
		    Rh.AddInteraction(num_loc, numpsi, vloc);
                    Rh_weight.AddInteraction(num_loc, numpsi, vloc*this->WeightsND(numero));
		    Rh_y.AddInteraction(num_loc, numpsi, vloc);
                    Rh_weight_y.AddInteraction(num_loc, numpsi, vloc*this->WeightsND(numero));
		    if (numpsi < 3*this->nb_dof_boundaries)
		      RhB.AddInteraction(num_loc, numpsi, vloc);
		  }
	      }
	    
	    l = i; m = j;
	    // part dU/dz v_z
	    for (n = 0; n <= order_quad; n++)
	      {
		numero = NumQuad3D(l,m,n);
		numpsi = 3*numero+2; // vectorial dof number
		
		vloc = lob_quad.GradPhi(k, n);
		if (abs(vloc) >= 1e4*epsilon_machine)
		  {
		    // dphi_{i,j,k} / dx (\xi_l,\xi_m,\xi_n) = \phi'_k(\xi_n) 
		    Rh.AddInteraction(num_loc, numpsi, vloc);
                    Rh_weight.AddInteraction(num_loc, numpsi, vloc*this->WeightsND(numero));
                    Rh_z.AddInteraction(num_loc, numpsi, vloc);
                    Rh_weight_z.AddInteraction(num_loc, numpsi, vloc*this->WeightsND(numero));
		    if (numpsi < 3*this->nb_dof_boundaries)
		      RhB.AddInteraction(num_loc, numpsi, vloc);
		  }
	      }
	  }
    
    // conversion to CSR format fore efficiency
    Seldon::Copy(Rh, rh_loc);
    Seldon::Copy(Rh_weight, rh_weight);
    Seldon::Copy(Rh_x, rh_locX);
    Seldon::Copy(Rh_y, rh_locY);
    Seldon::Copy(Rh_z, rh_locZ);
    Seldon::Copy(Rh_weight_x, rh_weightX);
    Seldon::Copy(Rh_weight_y, rh_weightY);
    Seldon::Copy(Rh_weight_z, rh_weightZ);
    Seldon::Copy(RhB, rh_boundary);
  }
  
  
  //! construction of mass matrix
  void HexahedronGauss::ConstructMassMatrix()
  {
    dGL_G.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        dGL_G(i,j) = lob_basis.EvaluatePhiGrad(i, this->Points1D(j));
    
    // computation of operator Ch
    int N = order_quad + 1;
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Ch3;
    Ch3.Reallocate((order+1)*(order+1)*N, (order+1)*(order+1)*(order+1));
    Ch2.Reallocate((order+1)*N*N, (order+1)*(order+1)*N);
    Ch1.Reallocate(N*N*N, (order+1)*N*N);
    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  for (int j = 0; j < N; j++)
	    Ch3.AddInteraction(i1*(order+1)*N+i2*N+j, NumDofs3D(i1,i2,i3), GL_G(i3,j));
    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 < N; i3++)
	  for (int j = 0; j < N; j++)
	    Ch2.AddInteraction(i1*N*N+j*N+i3, i1*(order+1)*N+i2*N+i3, GL_G(i2,j));
    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 < N; i2++)
	for (int i3 = 0; i3 < N; i3++)
	  for (int j = 0; j < N; j++)
	    Ch1.AddInteraction(NumQuad3D(j,i2,i3), i1*N*N+i2*N+i3, GL_G(i1,j));
    
    Seldon::Copy(Ch1, ch1_loc); Seldon::Copy(Ch2, ch2_loc); Seldon::Copy(Ch3, ch3_loc);
  }


  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
   */  
  template<class Vector1, class Vector2>
  void HexahedronGauss::ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
                                                       Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
    
  
  //! computation of U on nodal points
  /*!
    \param[in] Un components of U on dofs
    \param[out] Unode values of U on nodal points
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    int order_geom = this->GetGeometryOrder();
    if (order == order_geom)
      Unode = Un;
    else
      {
        Unode.Reallocate(this->GetNbPointsNodalElt());
        Vector2 z((order+1)*(order_geom+1)*(order_geom+1));
        Vector2 y((order+1)*(order+1)*(order_geom+1));
        Mlt(ch1_node, Un, y);
        Mlt(ch2_node, y, z);
        Mlt(ch3_node, z, Unode);
      }
  }

  
  //! computation of u on nodal points of a face
  /*!
    \param[in] Un dof components of u
    \param[out] Unode values of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1, class Vector2>
  void HexahedronGauss
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    int order_geom = this->GetGeometryOrder();
    Unode.Reallocate((order_geom+1)*(order_geom+1));
    const Matrix<int>& FacesNodal = Fb_geom.GetNodalNumber();
    const Matrix<int>& NumNodes2D_quad = Fb_geom.GetNumNodes2D_quad();
    if (order_geom == order)
      for (int j = 0; j < Unode.GetM(); j++)
        {
          int node = FacesNodal(j, num_loc);
          Unode(j) = Un(node);
        }
    else
      {
        Vector1 x((order_geom+1)*(order+1));
        x.Fill(0);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order_geom; j++)
            {
              for (int k = 0; k <= order; k++)
                {
                  int node = FacesDof(NumDofs2D(i, k), num_loc);
                  x(i*(order_geom+1) + j) += GL_Geom(k, j)*Un(node);
                }
            }
        
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            {
              int node = NumNodes2D_quad(i, j);
              Unode(node) = 0;
              for (int k = 0; k <= order; k++)
                Unode(node) += GL_Geom(k, i)*x(k*(order_geom+1) + j);
            }
      }
  }
  

  //! computation of grad u on nodal points of a face
  /*!
    \param[in] Un dof components of u
    \param[out] Unode gradient of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */  
  template<class Vector1, class Vector2>
  void HexahedronGauss
  ::ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Vector2 Udof(3*(order+1)*(order+1));    
    int ext = 0;
    if (num_loc > 2)
      ext = order;
    
    switch (num_loc)
      {
      case 0 :
      case 5 :
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            {
              int node = 3*NumDofs2D(i, j);
              Udof(node) = 0; Udof(node+1) = 0; Udof(node+2) = 0;
              for (int k = 0; k <= order; k++)
                {
                  Udof(node) += lob_basis.GradPhi(k, ext)*Un(NumDofs3D(k, i, j));
                  Udof(node+1) += lob_basis.GradPhi(k, i)*Un(NumDofs3D(ext, k, j));
                  Udof(node+2) += lob_basis.GradPhi(k, j)*Un(NumDofs3D(ext, i, k));
                }
            }
        break;
      case 1 :
      case 4 :
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            {
              int node = 3*NumDofs2D(i, j);
              Udof(node) = 0; Udof(node+1) = 0; Udof(node+2) = 0;
              for (int k = 0; k <= order; k++)
                {
                  Udof(node+1) += lob_basis.GradPhi(k, ext)*Un(NumDofs3D(i, k, j));
                  Udof(node) += lob_basis.GradPhi(k, i)*Un(NumDofs3D(k, ext, j));
                  Udof(node+2) += lob_basis.GradPhi(k, j)*Un(NumDofs3D(i, ext, k));
                }
            }
        break;
      case 2 :
      case 3 :
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            {
              int node = 3*NumDofs2D(i, j);
              Udof(node) = 0; Udof(node+1) = 0; Udof(node+2) = 0;
              for (int k = 0; k <= order; k++)
                {
                  Udof(node+2) += lob_basis.GradPhi(k, ext)*Un(NumDofs3D(i, j, k));
                  Udof(node) += lob_basis.GradPhi(k, i)*Un(NumDofs3D(k, j, ext));
                  Udof(node+1) += lob_basis.GradPhi(k, j)*Un(NumDofs3D(i, k, ext));
                }
            }
        break;
      }

    int order_geom = this->GetGeometryOrder();    
    Unode.Reallocate(3*(order_geom+1)*(order_geom+1));
    if (order_geom == order)
      {
        for (int j = 0; j < Unode.GetM(); j++)
          Unode(j) = Udof(j);
      }
    else
      {
        Vector2 Ux(3*(order+1)*(order_geom+1));
        Ux.Fill(0);
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            {
              int node = 3*(order+1)*i + 3*j;
              for (int k = 0; k <= order; k++)
                {
                  int p = 3*NumDofs2D(k, j);
                  Ux(node) += GL_Geom(k, i)*Udof(p);
                  Ux(node+1) += GL_Geom(k, i)*Udof(p+1);
                  Ux(node+2) += GL_Geom(k, i)*Udof(p+2);
                }
            }
        
        Unode.Fill(0);
	const Matrix<int>& NumNodes2D_quad = Fb_geom.GetNumNodes2D_quad();
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            {
              int node = 3*NumNodes2D_quad(i, j);
              for (int k = 0; k <= order; k++)
                {
                  int p = 3*((order+1)*i + k);
                  Unode(node) += GL_Geom(k, j)*Ux(p);
                  Unode(node+1) += GL_Geom(k, j)*Ux(p+1);
                  Unode(node+2) += GL_Geom(k, j)*Ux(p+2);
                }
            }
      }
  }
  
  
  //! x is overwritten by M^-1 x where M is the mass matrix
  template<class Vector1>
  void HexahedronGauss::SolveMassMatrixGen(Vector1& x) const
  {
    SolveCholesky(SeldonNoTrans, x);
    SolveCholesky(SeldonTrans, x);
  }
  
  
  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix
  template<class Vector1>
  void HexahedronGauss::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    if (order != order_quad)
      {
        abort();
      }
    
    Vector1 y(nb_dof_loc); typename Vector1::value_type vloc;
    if (TransA.Trans())
      {
        for (int i = 0; i < y.GetM(); i++)
          y(i) = x(i)*invSqrtWeights3d(i);
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                vloc = 0;
                for (int m = 0; m <= order; m++)
                  vloc += invGL_G(m, i)*y(NumDofs3D(m, j, k));
                
                x(NumDofs3D(i, j, k)) = vloc;
              }
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                vloc = 0;
                for (int m = 0; m <= order; m++)
                  vloc += invGL_G(m, j)*x(NumDofs3D(i, m, k));
                
                y(NumDofs3D(i, j, k)) = vloc;
              }
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                vloc = 0;
                for (int m = 0; m <= order; m++)
                  vloc += invGL_G(m, k)*y(NumDofs3D(i, j, m));
                
                x(NumDofs3D(i, j, k)) = vloc;
              }
      }
    else
      {
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                vloc = 0;
                for (int m = 0; m <= order; m++)
                  vloc += invGL_G(k, m)*x(NumDofs3D(i, j, m));
                
                y(NumDofs3D(i, j, k)) = vloc;
              }

        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                vloc = 0;
                for (int m = 0; m <= order; m++)
                  vloc += invGL_G(j, m)*y(NumDofs3D(i, m, k));
                
                x(NumDofs3D(i, j, k)) = vloc;
              }
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                vloc = 0;
                for (int m = 0; m <= order; m++)
                  vloc += invGL_G(i, m)*x(NumDofs3D(m, j, k));
                
                y(NumDofs3D(i, j, k)) = vloc;
              }
        
        for (int i = 0; i < y.GetM(); i++)
          x(i) = y(i)*invSqrtWeights3d(i);
      }
  }
    
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector2>
  void HexahedronGauss::MltMassMatrixGen(Vector2& Vh) const
  {
    Vector2 x(nb_points_quadrature_inside);
    ApplyChTranspose(Vh, x);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      x(i) *= this->WeightsND(i);
    
    ApplyCh(x, Vh);    
  }

  
  //! Integration against basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Ch Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 Vh1(ch2_loc.GetM()), Vh2(ch3_loc.GetM());
    Mlt(SeldonTrans, ch1_loc, Uh, Vh1);
    Mlt(SeldonTrans, ch2_loc, Vh1, Vh2);
    Mlt(SeldonTrans, ch3_loc, Vh2, Vh);
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
  void HexahedronGauss::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {    
    Vector2 Vh1(ch3_loc.GetM()), Vh2(ch2_loc.GetM());
    Mlt(ch3_loc, Uh, Vh1);
    Mlt(ch2_loc, Vh1, Vh2);
    Mlt(ch1_loc, Vh2, Vh);
  }
  
  
  //! Integration against gradient of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \varphi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \varphi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
    Mlt(rh_loc, Uh, Uquad);
    ApplyChGen(Uquad, Vh);    
  }
   
  
  //! computation of gradient of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j grad phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = grad phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
    ApplyChTransposeGen(Uh, Uquad);
    Mlt(SeldonTrans, rh_loc, Uquad, Vh);
  }
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 dUquad(3*nb_points_quadrature_inside), Uquad(nb_points_quadrature_inside);
    Vector1 Ux(nb_dof_loc), Uy(nb_dof_loc), Uz(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      {
        Ux(i) = Uh(3*i);
        Uy(i) = Uh(3*i+1);
        Uz(i) = Uh(3*i+2);
      }
    
    Vector1 dUx_quad(nb_points_quadrature_inside);
    Vector1 dUy_quad(nb_points_quadrature_inside);
    Vector1 dUz_quad(nb_points_quadrature_inside);

    ApplyChTransposeGen(Ux, dUx_quad);
    ApplyChTransposeGen(Uy, dUy_quad);
    ApplyChTransposeGen(Uz, dUz_quad);
    
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        dUquad(3*i) = dUx_quad(i);
        dUquad(3*i+1) = dUy_quad(i);
        dUquad(3*i+2) = dUz_quad(i);
      }
    
    Mlt(rh_weight, dUquad, Uquad);
    ApplyChGen(Uquad, Vh);
  }
   
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside), Vquad(3*nb_points_quadrature_inside);
    ApplyChTransposeGen(Uh, Uquad);
    Mlt(SeldonTrans, rh_weight, Uquad, Vquad);

    Vector1 dUx_quad(nb_points_quadrature_inside);
    Vector1 dUy_quad(nb_points_quadrature_inside);
    Vector1 dUz_quad(nb_points_quadrature_inside);

    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        dUx_quad(i) = Vquad(3*i);
        dUy_quad(i) = Vquad(3*i+1);
        dUz_quad(i) = Vquad(3*i+2);
      }

    Vector1 Ux(nb_dof_loc), Uy(nb_dof_loc), Uz(nb_dof_loc);    
    Ux.Fill(0); Uy.Fill(0); Uz.Fill(0);
    ApplyChGen(dUx_quad, Ux);
    ApplyChGen(dUy_quad, Uy);
    ApplyChGen(dUz_quad, Uz);
    
    for (int i = 0; i < nb_dof_loc; i++)
      {
        Vh(3*i) = Ux(i);
        Vh(3*i+1) = Uy(i);
        Vh(3*i+2) = Uz(i);
      }
  }

  
  //! Integration against derivatives of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Ux Ux_i = \int_K f d/dx(\varphi_i)  dx
    \param[out] Uy Uy_i = \int_K f d/dy(\varphi_i)  dx
    \param[out] Uz Uz_i = \int_K f d/dz(\varphi_i)  dx
    This operation is equivalent to a matrix vector product
    Ux = Rh^1 Uh, Uy = Rh^2 Uh, Uz = Rh^3 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
    where (Rh)^3_{i,j} = d/dz(\varphi_i(\xi_j) )
  */  
  template<class Vector1, class Vector2>
  void HexahedronGauss
  ::ApplyRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
    Mlt(rh_locX, Uh, Ux);
    Mlt(rh_locY, Uh, Uy);
    Mlt(rh_locZ, Uh, Uz);
  }
  

  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    (Rh^3)_{i, j} = \int_K \varphi_j d/dz(\varphi_i) dx
    then Ux = Rh^1 Uh, Uy = Rh^2 Uh, Uz = Rh^3 Uh
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss
  ::ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy, Vector2& Uz) const
  {
    Vector1 dUquad(3*nb_points_quadrature_inside), Uquad(nb_points_quadrature_inside);
    for (int i = 0; i < nb_dof_loc; i++)
      {
        Ux(i) = Uh(3*i);
        Uy(i) = Uh(3*i+1);
        Uz(i) = Uh(3*i+2);
      }
    
    Vector1 dUx_quad(nb_points_quadrature_inside);
    Vector1 dUy_quad(nb_points_quadrature_inside);
    Vector1 dUz_quad(nb_points_quadrature_inside);

    ApplyChTransposeGen(Ux, dUx_quad);
    ApplyChTransposeGen(Uy, dUy_quad);
    ApplyChTransposeGen(Uz, dUz_quad);
    
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        dUquad(3*i) = dUx_quad(i);
        dUquad(3*i+1) = dUy_quad(i);
        dUquad(3*i+2) = dUz_quad(i);
      }
    
    Mlt(rh_weightX, dUquad, dUx_quad);
    Mlt(rh_weightY, dUquad, dUy_quad);
    Mlt(rh_weightZ, dUquad, dUz_quad);
    
    ApplyChGen(dUx_quad, Ux);
    ApplyChGen(dUy_quad, Uy);
    ApplyChGen(dUz_quad, Uz);
  }
  
  
  //! Integration against gradient of basis functions associated with quadrature points
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh Vh_i = \int_K f \nabla \psi_i dx
    This operation is equivalent to a matrix vector product
    Vh = Rh Uh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(rh_loc, Uh, Vh);
  }
   
  
  //! Computation of gradient on quadrature points from values on quadrature points
  /*!
    \param[in] Uh values of u on quadrature points
    \param[out] Vh gradient of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Rh* Uh
    where (Rh)_{i,j} = \nabla \psi_i(\xi_j)
  */
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }
  
  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the face
    \param[in] r order of quadrature rule of the face
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((r == 0) || (r == order_quad))
      {
	typename Vector1::value_type vloc;
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order_quad; j++)
	    {
	      vloc = 0;
	      for (int k = 0; k <= order; k++)
		vloc += GL_G(k, j)*Uh(FacesDof(NumDofs2D(i, k), num_loc));
	      
	      for (int k = 0; k <= order_quad; k++)
		Vh(NumQuad2D(k, j)) += GL_G(i, k)*vloc;
	    }
      }
    else
      {
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= r; j++)
	    {
	      vloc = 0;
	      for (int k = 0; k <= order; k++)
		vloc += Ah(k, j)*Uh(FacesDof(NumDofs2D(i, k), num_loc));
	      
	      for (int k = 0; k <= r; k++)
		Vh(Num(k, j)) += Ah(i, k)*vloc;
	    }
      }
    
  }
  
  
  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
  */
  template<class T0, class Vector1, class Vector2>
  void HexahedronGauss
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((r == 0) || (r == order_quad))
      {
	typename Vector1::value_type vloc;
	for (int i = 0; i <= order_quad; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      vloc = 0;
	      for (int k = 0; k <= order_quad; k++)
		vloc += GL_G(j, k)*Uh(NumQuad2D(i, k));
	      
	      vloc *= alpha;
	      for (int k = 0; k <= order; k++)
		Vh(FacesDof(NumDofs2D(k, j), num_loc)) += GL_G(k, i)*vloc;
	    }
      }
    else
      {
	typename Vector1::value_type vloc;
	const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      vloc = 0;
	      for (int k = 0; k <= r; k++)
		vloc += Ah(j, k)*Uh(Num(i, k));
	      
	      vloc *= alpha;
	      for (int k = 0; k <= order; k++)
		Vh(FacesDof(NumDofs2D(k, j), num_loc)) += Ah(k, i)*vloc;
	    }
      }
  }

  
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] res components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* res
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void HexahedronGauss
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& res, Vector2& Vh, int r) const
  {
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();

    int Nquad = (order_quad+1)*(order_quad+1);
    Vector2 Udof(Nquad), feval(3*Nquad);
    Udof.Fill(0);
    feval.Fill(0);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    switch (num_loc)
      {
      case 0 :
      case 5 :
        {
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                Udof(p) = res(NumDofs3D(ext, j, k));
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(i, j, k);
                    feval(3*p) += dGL_Geom(i, ext)*res(node);
                  }
              }
          
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    feval(3*NumQuad2D(k, j) + 1) += lob_basis.GradPhi(i, k)*Udof(num);
                    feval(3*NumQuad2D(i, k) + 2) += lob_basis.GradPhi(j, k)*Udof(num);
                  }
              }
        }
        break;
      case 1 :
      case 4 :
        {
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                Udof(p) = res(NumDofs3D(j, ext, k));
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, i, k);
                    feval(3*p+1) += dGL_Geom(i, ext)*res(node);                    
                  }
              }
          
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    feval(3*NumQuad2D(k, j)) += lob_basis.GradPhi(i, k)*Udof(num);
                    feval(3*NumQuad2D(i, k) + 2) += lob_basis.GradPhi(j, k)*Udof(num);
                  }
              }
        }
        break;
      case 2 :
      case 3 :
        {
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                Udof(p) += res(NumDofs3D(j, k, ext));
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, k, i);
                    feval(3*p+2) += dGL_Geom(i, ext)*res(node);                    
                  }
              }
          
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    feval(3*NumQuad2D(k, j)) += lob_basis.GradPhi(i, k)*Udof(num);
                    feval(3*NumQuad2D(i, k) + 1) += lob_basis.GradPhi(j, k)*Udof(num);
                  }
              }
        }
        break;
      }

    if ((r == 0) || (r == order))
      {
        typename Vector1::value_type vx, vy, vz;
        Vh.Fill(0);
        for (int i = 0; i <= order; i++)
	  for (int j = 0; j <= order_quad; j++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int k = 0; k <= order; k++)
		{
                  vx += GL_G(k, j)*feval(3*NumDofs2D(i, k));
                  vy += GL_G(k, j)*feval(3*NumDofs2D(i, k)+1);
                  vz += GL_G(k, j)*feval(3*NumDofs2D(i, k)+2);
                }
	      
	      for (int k = 0; k <= order_quad; k++)
		{
                  Vh(3*NumQuad2D(k, j)) += GL_G(i, k)*vx;
                  Vh(3*NumQuad2D(k, j)+1) += GL_G(i, k)*vy;
                  Vh(3*NumQuad2D(k, j)+2) += GL_G(i, k)*vz;
                }
	    }
      }
    else
      {
        typename Vector1::value_type vx, vy, vz;
	Vh.Fill(0);
        const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
        for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int j = 0; j <= order; j++)
		{
                  vx += Ah(j, m)*feval(3*NumDofs2D(i, j));
                  vy += Ah(j, m)*feval(3*NumDofs2D(i, j)+1);
                  vz += Ah(j, m)*feval(3*NumDofs2D(i, j)+2);
                }
	      
	      for (int n = 0; n <= r; n++)
		{
                  Vh(3*Num(n, m)) += Ah(i, n)*vx;
                  Vh(3*Num(n, m)+1) += Ah(i, n)*vy;
                  Vh(3*Num(n, m)+2) += Ah(i, n)*vz;
                }
	    }
      }
  }
  
  
  //! integration against gradient of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh  vector containing \omega_k f(\xi_k) 
    \param[out] res res_i = res_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void HexahedronGauss::
  ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& res, int r) const
  {
    int Nquad = (order_quad+1)*(order_quad+1);
    Vector2 feval(3*Nquad);

    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    if ((r != 0) && (r != this->order)) 
      {
        typename Vector1::value_type vx, vy, vz;
        feval.Fill(0);
        const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	const Matrix<int>& Num = NumProjOperator(r);
	for (int i = 0; i <= order; i++)
	  for (int m = 0; m <= r; m++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int n = 0; n <= r; n++)
		{
                  vx += Ah(i, n)*Uh(3*Num(n, m));
                  vy += Ah(i, n)*Uh(3*Num(n, m)+1);
                  vz += Ah(i, n)*Uh(3*Num(n, m)+2);
                }
	      
	      for (int j = 0; j <= order; j++)
		{
                  feval(3*NumDofs2D(i, j)) += Ah(j, m)*vx;
                  feval(3*NumDofs2D(i, j)+1) += Ah(j, m)*vy;
                  feval(3*NumDofs2D(i, j)+2) += Ah(j, m)*vz;
                }
	    }
      }    
    else
      {
        typename Vector1::value_type vx, vy, vz;
        feval.Fill(0);
	for (int i = 0; i <= order_quad; i++)
	  for (int j = 0; j <= order; j++)
	    {
	      vx = 0; vy = 0; vz = 0;
	      for (int k = 0; k <= order_quad; k++)
		{
                  vx += GL_G(j, k)*Uh(3*NumQuad2D(i, k));
                  vy += GL_G(j, k)*Uh(3*NumQuad2D(i, k)+1);
                  vz += GL_G(j, k)*Uh(3*NumQuad2D(i, k)+2);
                }
	      
              for (int k = 0; k <= order; k++)
		{
                  feval(3*NumDofs2D(k, j)) += GL_G(k, i)*vx;
                  feval(3*NumDofs2D(k, j)+1) += GL_G(k, i)*vy;
                  feval(3*NumDofs2D(k, j)+2) += GL_G(k, i)*vz;
                }
	    }
      }
    
    Mlt(alpha, feval);
    
    int ext = 0;
    if (num_loc > 2)
      ext = this->GetGeometryOrder();

    Vector2 Udof(Nquad);
    Udof.Fill(0);
    switch (num_loc)
      {
      case 0 :
      case 5 :
        {
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    Udof(num) += lob_basis.GradPhi(i, k)*feval(3*NumQuad2D(k, j) + 1);
                    Udof(num) += lob_basis.GradPhi(j, k)*feval(3*NumQuad2D(i, k) + 2);
                  }
              }
          
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                res(NumDofs3D(ext, j, k)) += Udof(p);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(i, j, k);
                    res(node) += dGL_Geom(i, ext)*feval(3*p);
                  }
              }

        }
        break;
      case 1 :
      case 4 :
        {
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    Udof(num) += lob_basis.GradPhi(i, k)*feval(3*NumQuad2D(k, j));
                    Udof(num) += lob_basis.GradPhi(j, k)*feval(3*NumQuad2D(i, k) + 2);
                  }
              }
          
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                res(NumDofs3D(j, ext, k)) += Udof(p);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, i, k);
                    res(node) += dGL_Geom(i, ext)*feval(3*p+1);
                  }
              }
        }
        break;
      case 2 :
      case 3 :
        {
          for (int i = 0; i <= order; i++)
            for (int j = 0; j <= order; j++)
              {
                int num = NumQuad2D(i, j);
                for (int k = 0; k <= order; k++)
                  {
                    Udof(num) += lob_basis.GradPhi(i, k)*feval(3*NumQuad2D(k, j));
                    Udof(num) += lob_basis.GradPhi(j, k)*feval(3*NumQuad2D(i, k) + 1);
                  }
              }
          
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              {
                int p = NumQuad2D(j, k);
                res(NumDofs3D(j, k, ext)) += Udof(p);
                for (int i = 0; i <= order; i++)
                  {
                    int node = NumDofs3D(j, k, i);
                    res(node) += dGL_Geom(i, ext)*feval(3*p+2);
                  }
              }
        }
        break;
      }
  }
  
  
  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T>
  void HexahedronGauss::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    T zero; SetComplexZero(zero);
    if (mass == zero)
      return;
    
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;
    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  {
	    int i = NumDofs3D(i1, i2, i3);
	    for (int j1 = 0; j1 <= order; j1++)
	      for (int j2 = 0; j2 <= order; j2++)
		for (int j3 = 0; j3 <= order; j3++)
		  {
		    int j = NumDofs3D(j1, j2, j3);
		    val(j) = mass*mass_matrix1D(i1, j1)
		      *mass_matrix1D(i2, j2)*mass_matrix1D(i3, j3);
		  }
	    
	    A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
	  }    
  }
  
  
  //! we add elementary matrix for a constant jacobian
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass, C, D, E coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p \sum_{p,q}  C(p, q) S^{p,q}
    where S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
    where R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void HexahedronGauss
  ::AddConstantElemMatrixGen(int m, int n, const T& mass, const TinyMatrix<T, Prop, 3, 3>& C,
			     const TinyVector<T, 3>& D, const TinyVector<T, 3>& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;
    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
        for (int i3 = 0; i3 <= order; i3++)
          {
	    val.Zero();
            int i = NumDofs3D(i1, i2, i3);
	    
            for (int j1 = 0; j1 <= order; j1++)
              for (int j2 = 0; j2 <= order; j2++)
                for (int j3 = 0; j3 <= order; j3++)
                  {
                    int j = NumDofs3D(j1, j2, j3);
		    if (!null_term(0))
		      val(j) = mass*mass_matrix1D(i1, j1)
			*mass_matrix1D(i2, j2)*mass_matrix1D(i3, j3);
		    
		    if (!null_term(1))
		      {
			val(j) += C(0, 0)*stiffness_matrix1D(i1, j1)
                          *mass_matrix1D(i2, j2)*mass_matrix1D(i3, j3);
                        
                        val(j) += C(1, 1)*stiffness_matrix1D(i2, j2)
                          *mass_matrix1D(i1, j1)*mass_matrix1D(i3, j3);
                        
                        val(j) += C(2, 2)*stiffness_matrix1D(i3, j3)
                          *mass_matrix1D(i1, j1)*mass_matrix1D(i2, j2);
                        
                        val(j) += C(0, 1)*gradient_matrix1D(i1, j1)
                          *gradient_matrix1D(j2, i2)*mass_matrix1D(i3, j3);
                        
                        val(j) += C(1, 0)*gradient_matrix1D(j1, i1)
                          *gradient_matrix1D(i2, j2)*mass_matrix1D(i3, j3);
                        
                        val(j) += C(0, 2)*gradient_matrix1D(i1, j1)
                          *gradient_matrix1D(j3, i3)*mass_matrix1D(i2, j2);
                        
                        val(j) += C(2, 0)*gradient_matrix1D(j1, i1)
                          *gradient_matrix1D(i3, j3)*mass_matrix1D(i2, j2);
                        
                        val(j) += C(1, 2)*gradient_matrix1D(i2, j2)
                          *gradient_matrix1D(j3, i3)*mass_matrix1D(i1, j1);
                        
                        val(j) += C(2, 1)*gradient_matrix1D(j2, i2)
                          *gradient_matrix1D(i3, j3)*mass_matrix1D(i1, j1);
		      }
		    
		    if (!null_term(2))
		      {
                        val(j) += D(0)*gradient_matrix1D(i1, j1)
                          *mass_matrix1D(i2, j2)*mass_matrix1D(i3, j3);
                        
                        val(j) += D(1)*gradient_matrix1D(i2, j2)
                          *mass_matrix1D(i1, j1)*mass_matrix1D(i3, j3);
                        
                        val(j) += D(2)*gradient_matrix1D(i3, j3)
                          *mass_matrix1D(i2, j2)*mass_matrix1D(i1, j1);
		      }
		    
		    if (!null_term(3))
		      {
			val(j) += E(0)*gradient_matrix1D(j1, i1)
                          *mass_matrix1D(i2, j2)*mass_matrix1D(i3, j3);
                        
                        val(j) += E(1)*gradient_matrix1D(j2, i2)
                          *mass_matrix1D(i1, j1)*mass_matrix1D(i3, j3);
                        
                        val(j) += E(2)*gradient_matrix1D(j3, i3)
                          *mass_matrix1D(i1, j1)*mass_matrix1D(i2, j2);
                      }		    
                  }
	    
	    A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
          }
  }
  
  
  //! adds elementary matrix for a variable jacobian J_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A, C, D, E coefficients multiplied by omega_i J_i for each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + S + R + Rt
    where M_ij = \int A varphi_j varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
    where R_ij = \int E varphi_i  \nabla varphi_j
    where Rt_ij = \int D varphi_j  \nabla varphi_i
   */
  template<class T, class Prop>
  void HexahedronGauss
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& A,
			     const Vector<TinyMatrix<T, Prop, 3, 3> >& C,
			     const Vector<TinyVector<T, 3> >& D,
			     const Vector<TinyVector<T, 3> >& E,
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension3, 1>::
      AddVariableElemMatrixOpt(off_row, off_col, A, C, D, E, null_term, mat);
  }
    
  
  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
   */
  template<class T>
  void HexahedronGauss
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension3, 1>::AddVariableMassMatrixOpt(off_row, off_col, A, mat);
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] phi values of basis functions on point
    */
  void HexahedronGauss::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    int r = order;
    
    phi.Reallocate(nb_dof_loc);
    // computation of 1-D basis functions
    VectReal_wp phix(r+1), phiy(r+1), phiz(r+1);
    if (r >= 8)
      {
        lob_basis.ComputeValuesPhiRef(point_loc(0), phix);
        lob_basis.ComputeValuesPhiRef(point_loc(1), phiy);
        lob_basis.ComputeValuesPhiRef(point_loc(2), phiz);
      }
    else
      for (int i = 0; i <= r; i++)
        {
          phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
          phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
          phiz(i) = lob_basis.EvaluatePhi(i, point_loc(2));
        }
    
    // now tensorization
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        for (int k = 0; k <= r; k++)
          phi(NumDofs3D(i, j, k)) = phix(i)*phiy(j)*phiz(k);
    
  }

  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res gradient of basis functions on point
  */
  void HexahedronGauss::ComputeGradientPhiRef(const R3& point_loc, VectR3& res) const
  {
    int r = order;
    
    res.Reallocate(nb_dof_loc);
    // computation of 1-D basis functions
    VectReal_wp phix(r+1), phiy(r+1), phiz(r+1);
    VectReal_wp dphix(r+1), dphiy(r+1), dphiz(r+1);
    for (int i = 0; i <= r; i++)
      {
        phix(i) = lob_basis.EvaluatePhi(i, point_loc(0));
        dphix(i) = lob_basis.EvaluatePhiGrad(i, point_loc(0));
        phiy(i) = lob_basis.EvaluatePhi(i, point_loc(1));
        dphiy(i) = lob_basis.EvaluatePhiGrad(i, point_loc(1));
        phiz(i) = lob_basis.EvaluatePhi(i, point_loc(2));
        dphiz(i) = lob_basis.EvaluatePhiGrad(i, point_loc(2));
      }
    
    // now tensorization
    R3 grad;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        for (int k = 0; k <= r; k++)
          {
            grad(0) = dphix(i)*phiy(j)*phiz(k);
            grad(1) = phix(i)*dphiy(j)*phiz(k);
            grad(2) = phix(i)*phiy(j)*dphiz(k);
            res(NumDofs3D(i, j, k)) = grad;
          }
    
  }

  
  //! retrieves values of a single basis function on all quadrature points
  void HexahedronGauss::GetValueSinglePhiQuadrature(int i, VectReal_wp& phi) const
  {
    int i1 = CoordinateDofs(i, 0);
    int i2 = CoordinateDofs(i, 1);
    int i3 = CoordinateDofs(i, 2);
    for (int k1 = 0; k1 <= order; k1++)
      for (int k2 = 0; k2 <= order; k2++)
	for (int k3 = 0; k3 <= order; k3++)
	  phi(NumQuad3D(k1, k2, k3)) = GL_G(i1, k1)*GL_G(i2, k2)*GL_G(i3, k3);
  }


  //! retrieves values and gradients of a single basis function on all quadrature points
  void HexahedronGauss
  ::GetGradientSinglePhiQuadrature(int i, VectReal_wp& phi, VectR3& grad_phi) const
  {
    int i1 = CoordinateDofs(i, 0);
    int i2 = CoordinateDofs(i, 1);
    int i3 = CoordinateDofs(i, 2);
    for (int k1 = 0; k1 <= order; k1++)
      for (int k2 = 0; k2 <= order; k2++)
	for (int k3 = 0; k3 <= order; k3++)
	  {
	    int n = NumQuad3D(k1, k2, k3);
	    phi(n) = GL_G(i1, k1)*GL_G(i2, k2)*GL_G(i3, k3);
	    grad_phi(n)(0) = dGL_G(i1, k1)*GL_G(i2, k2)*GL_G(i3, k3);
	    grad_phi(n)(1) = GL_G(i1, k1)*dGL_G(i2, k2)*GL_G(i3, k3);
	    grad_phi(n)(2) = GL_G(i1, k1)*GL_G(i2, k2)*dGL_G(i3, k3);
	  }
  }


  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */  
  void HexahedronGauss::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad3D(k, 0);
        int k2 = CoordinateQuad3D(k, 1);
        int k3 = CoordinateQuad3D(k, 2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              phi(NumDofs3D(i, j, k)) = GL_G(i, k1)*GL_G(j, k2)*GL_G(k, k3);
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/((order_quad+1)*(order_quad+1));
        int ks = (k-nb_points_quadrature_inside)%((order_quad+1)*(order_quad+1));
        int k1 = CoordinateQuad2D(ks, 0);
        int k2 = CoordinateQuad2D(ks, 1);
        int ext = 0;
        if (num_loc > 2)
          ext = order;
        
        switch (num_loc)
          {
          case 0 :
          case 5 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  phi(NumDofs3D(ext, i, j)) = GL_G(i, k1)*GL_G(j, k2);
            }
            break;
          case 1:
          case 4 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  phi(NumDofs3D(i, ext, j)) = GL_G(i, k1)*GL_G(j, k2);
            }
            break;
          case 2:
          case 3 :
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  phi(NumDofs3D(i, j, ext)) = GL_G(i, k1)*GL_G(j, k2);
            }
            break;
          }
      }
  }
  
  
  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
   */
  void HexahedronGauss::GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    FillZero(grad_phi);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad3D(k, 0);
        int k2 = CoordinateQuad3D(k, 1);
        int k3 = CoordinateQuad3D(k, 2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              grad_phi(NumDofs3D(i, j, k)).Init(dGL_G(i, k1)*GL_G(j, k2)*GL_G(k, k3),
                                                GL_G(i, k1)*dGL_G(j, k2)*GL_G(k, k3),
                                                GL_G(i, k1)*GL_G(j, k2)*dGL_G(k, k3));
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/((order_quad+1)*(order_quad+1));
        int ks = (k-nb_points_quadrature_inside)%((order_quad+1)*(order_quad+1));
        int k1 = CoordinateQuad2D(ks, 0);
        int k2 = CoordinateQuad2D(ks, 1);
        int ext = 0;
        if (num_loc > 2)
          ext = order;
        
        switch (num_loc)
          {
          case 0:
          case 5:
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    grad_phi(NumDofs3D(ext, i, j))(1) = dGL_G(i, k1)*GL_G(j, k2);
                    grad_phi(NumDofs3D(ext, i, j))(2) = GL_G(i, k1)*dGL_G(j, k2);
                    
                    for (int m = 0; m <= order; m++)
                      grad_phi(NumDofs3D(m, i, j))(0)
                        = lob_basis.GradPhi(m, ext)*GL_G(i, k1)*GL_G(j, k2);
                  }
            }
            break;
          case 1:
          case 4:
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    grad_phi(NumDofs3D(i, ext, j))(0) = dGL_G(i, k1)*GL_G(j, k2);
                    grad_phi(NumDofs3D(i, ext, j))(2) = GL_G(i, k1)*dGL_G(j, k2);
                    
                    for (int m = 0; m <= order; m++)
                      grad_phi(NumDofs3D(i, m, j))(1)
                        = lob_basis.GradPhi(m, ext)*GL_G(i, k1)*GL_G(j, k2);
                  }
            }
            break;
          case 2:
          case 3:
            {
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    grad_phi(NumDofs3D(i, j, ext))(0) = dGL_G(i, k1)*GL_G(j, k2);
                    grad_phi(NumDofs3D(i, j, ext))(1) = GL_G(i, k1)*dGL_G(j, k2);
                    
                    for (int m = 0; m <= order; m++)
                      grad_phi(NumDofs3D(i, j, m))(2)
                        = lob_basis.GradPhi(m, ext)*GL_G(i, k1)*GL_G(j, k2);
                  }
            }
            break;
          }
      }
  }

  
  //! computation of projection between finite element of different orders    
  void HexahedronGauss::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi; VectR2 Points; IVect perm;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    NumProjOperator.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate(order+1, r+1);
	NumProjOperator(r).Reallocate(r+1, r+1);
	Points = Pts(i); perm.Reallocate(Points.GetM()); perm.Fill();
	Sort(Points, perm);
	int nb = 0;
	for (int j1 = 0; j1 <= r; j1++)
	  for (int j2 = 0; j2 <= r; j2++)
	    NumProjOperator(r)(j1, j2) = perm(nb++);
	
	for (int i1 = 0; i1 <= order; i1++)
	  for (int i2 = 0; i2 <= r; i2++)
	    ProjOperatorQuadOrder(r)(i1, i2) = lob_basis.EvaluatePhi(i1, Points(i2)(1));
	
      }
  }


  //! selects dofs near the current dof
  /*!
    \param[in] pos current dof
    \param[in] DofUsed some dofs are already in use and can't be selected
    \param[out] ListeDof selected dofs
    \param[in] nb_dof number of dofs to select
   */
  void HexahedronGauss::
  PickNearDofs(int pos, const VectBool& DofUsed, IVect& ListeDof, int nb_dof) const
  {
    if (nb_dof <= 0)
      return;
    
    const Matrix<int>& CoordinateNodes = Fb_geom.GetCoordinateNodes3D();
    int i0 = CoordinateNodes(pos,0);
    int j0 = CoordinateNodes(pos,1), k0 = CoordinateNodes(pos,2);
    ListeDof.Reallocate(nb_dof); ListeDof.Fill(-1);
    int k = 1, nb = 0, node;
    if (!DofUsed(pos))
      ListeDof(nb++) = pos;
    
    // DISP(DofUsed); DISP(i0); DISP(j0); DISP(k0); DISP(NumDofs3D);
    while (nb < nb_dof)
      {
	// loop on concentric cube at distance k
	if ((j0-k) >= 0)
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumDofs3D(m,j0-k,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((j0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,i0-k); m <= min(order,i0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumDofs3D(m,j0+k,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((i0-k) >= 0)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order+1,j0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumDofs3D(i0-k,m,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((i0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,j0-k+1); m < min(order+1,j0+k); m++)
	    for (int n = max(0,k0-k); n <= min(order,k0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumDofs3D(i0+k,m,n);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((k0-k) >= 0)&&(nb < nb_dof))
	  for (int m = max(0,i0-k+1); m < min(order+1,i0+k); m++)
	    for (int n = max(0,j0-k+1); n < min(order+1,j0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumDofs3D(m,n,k0-k);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	if (((k0+k) <= order)&&(nb < nb_dof))
	  for (int m = max(0,i0-k+1); m < min(order+1,i0+k); m++)
	    for (int n = max(0,j0-k+1); n < min(order+1,j0+k); n++)
	      if (nb < nb_dof)
		{
		  node = NumDofs3D(m,n,k0+k);
		  if (!DofUsed(node))
		    ListeDof(nb++) = node;
		}
	
	k++;
	if (k > order)
	  nb = nb_dof;
	
      }
	
  }

} // namespace Montjoie
  
#define MONTJOIE_FILE_HEXAHEDRON_GAUSS_CXX
#endif
