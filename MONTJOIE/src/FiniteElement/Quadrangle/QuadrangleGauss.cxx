#ifndef MONTJOIE_FILE_QUADRANGLE_GAUSS_CXX

namespace Montjoie
{

  int QuadrangleGauss::type_quadrature_default(0);
  

  //! default constructor
  QuadrangleGauss::QuadrangleGauss() : QuadrangleReference<1>()
  {
    Fb_geom.quadrature_equal_nodal = false;
    Fb_geom.dof_equal_nodal = false;
    Fb_geom.dof_equal_quadrature = false;
    this->use_quadrature_for_rh = true;
  }


  size_t QuadrangleGauss::GetMemorySize() const
  {
    size_t taille = QuadrangleReference<1>::GetMemorySize();
    taille += GL_G.GetMemorySize();
    taille += invGL_G.GetMemorySize();
    taille += G_GL.GetMemorySize();
    taille += dGL_G.GetMemorySize();
    taille += GL_Geom.GetMemorySize();
    taille += dGL_Geom.GetMemorySize();
    taille += NumFct2D.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    taille += invWeights2d.GetMemorySize();
    taille += invSqrtWeights2d.GetMemorySize();
    taille += stiffness_matrix1D.GetMemorySize();
    taille += mass_matrix1D.GetMemorySize();
    taille += gradient_matrix1D.GetMemorySize();
    taille += proj_dof_to_nodal.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += ch1_node.GetMemorySize() + ch2_node.GetMemorySize();
    taille += rh_locX.GetMemorySize() + rh_locY.GetMemorySize();
    taille += ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize();
    return taille;
  }
  
  
  //! constructing finite element
  void QuadrangleGauss::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					       int rsurf, int type_surf)
  {
    if (type_quad == -1)
      type_quad = type_quadrature_default;

    if (rquad == 0)
      rquad = r;

    if (rsurf == 0)
      rsurf = rquad;
    
    QuadrangleReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();
    
    if (type_quad == Globatto<Real_wp>::QUADRATURE_LOBATTO)
      {
        VectReal_wp points1d_gauss, omega;
        ComputeGaussLegendre(points1d_gauss, omega, order_quad);
        
        GL_Gquad.Reallocate(r+1, rsurf+1);
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= rsurf; j++)
            GL_Gquad(i, j) = lob_basis.EvaluatePhi(i, points1d_gauss(j));
      }
    else
      GL_Gquad.Clear();
    
    GL_G.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        GL_G(i,j) = lob_basis.EvaluatePhi(i, lob_quad.Points(j));
    
    int order_geom = Fb_geom.GetGeometryOrder();
    GL_Geom.Reallocate(order + 1, order_geom + 1);
    dGL_Geom.Reallocate(order + 1, order_geom + 1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        {
          GL_Geom(i, j) = lob_basis.EvaluatePhi(i, Fb_geom.PointsNodal1D(j));
          dGL_Geom(i, j) = lob_basis.EvaluatePhiGrad(i, Fb_geom.PointsNodal1D(j));
        }
    
    rquad = order_quad;
    if (order == order_quad)
      {
        invGL_G = GL_G;
        GetInverse(invGL_G);
      }
    
    G_GL.Reallocate(rquad+1, r+1);
    dGL_G.Reallocate(r+1, rquad+1);
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= rquad; j++)
	{
	  G_GL(j, i) = lob_quad.EvaluatePhi(j, lob_basis.Points(i));
	  dGL_G(i, j) = lob_basis.EvaluatePhiGrad(i, Fb_geom.Points1D(j));
	}
    
    if ((type_quad != Globatto<Real_wp>::QUADRATURE_LOBATTO) || (order != order_quad))
      {
        ConstructElementaryMatrix(*this);
	
        Value_Phi.Clear(); Gradient_Phi.Clear();
        mass_matrix.Clear();
	if (r == rquad)
	  mass_matrix_chol.Clear();
        
	ValuePhi_Nodal.Clear();
      }
  }
  
  
  //! constructions of basis functions
  void QuadrangleGauss::ConstructFunctions()
  {
    int r = order;
    int order_geom = Fb_geom.GetGeometryOrder();
    if (order == order_geom)
      {
	this->elt_geom.dof_equal_nodal = true;
        nb_dof_loc = this->GetNbPointsNodalElt();
        nb_dof_boundaries = 4*order;
        
        lob_basis = this->GetNodalShapeFunctions1D();
        
        EdgesDof = this->GetNodalNumber();
        NumFct2D = this->GetNumNodes2D();
        CoordinateDofs = this->GetCoordinateNodes2D();
        
        this->SetPointsDof1D(this->PointsNodal1D());
	this->SetPointsDofND(this->PointsNodalND()); 
      }
    else
      {
	this->elt_geom.dof_equal_nodal = false;
        MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumFct2D, CoordinateDofs);
        
        lob_basis.ConstructQuadrature(r, lob_basis.QUADRATURE_LOBATTO);
        
        nb_dof_loc = (r+1)*(r+1);
        nb_dof_boundaries = 4*r;
        
        this->SetPointsDof1D(lob_basis.Points());
        VectR2 points_dof2d((r+1)*(r+1));
	
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            points_dof2d(NumFct2D(i,j))
	      = R2(this->PointsDof1D(i), this->PointsDof1D(j));
        
	this->SetPointsDofND(points_dof2d);
	
        // the dof points  on edges are Gauss-Lobatto points
        EdgesDof.Reallocate(r+1, 4);
        // First Edge
        for (int i = 0; i <= r; i++)
          EdgesDof(i, 0) = NumFct2D(i,0);
        
        // Second Edge
        for (int i = 0; i <= r; i++)
          EdgesDof(i, 1) = NumFct2D(r,i);
        
        // Third Edge
        for (int i = 0; i <= r; i++)
          EdgesDof(i, 2) = NumFct2D(r-i,r);
        
        // Fourth Edge
        for (int i = 0; i <= r; i++)
          EdgesDof(i, 3) = NumFct2D(0,r-i);
      }

    this->nb_points_dof_inside = this->nb_dof_loc;
    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(r+1);
	for (int i = 0; i <= r; i++)
	  this->num_dof_points_surf(n)(i) = EdgesDof(i, n);
      }

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);        
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    
    invWeights2d.Reallocate((order_quad+1)*(order_quad+1));
    invSqrtWeights2d.Reallocate((order_quad+1)*(order_quad+1));
    for (int k = 0; k < this->WeightsND().GetM(); k++)
      {
        invWeights2d(k) = 1.0/this->WeightsND(k);
        invSqrtWeights2d(k) = 1.0/sqrt(this->WeightsND(k));
      }
    
    Matrix<Real_wp> Phi1D, dPhi1D;
    Phi1D.Reallocate(order+1, this->Points1D().GetM());
    dPhi1D.Reallocate(order+1, this->Points1D().GetM());
    for (int i = 0; i <= order; i++)
      for (int k = 0; k < this->Points1D().GetM(); k++)
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
    
    // ch1_loc, ch2_loc, rh_loc
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2, Rh, RhX, RhY;
    Ch1.Reallocate((order+1)*(order_quad+1), (order_quad+1)*(order_quad+1));
    Ch2.Reallocate((order+1)*(order+1), (order+1)*(order_quad+1));
    Rh.Reallocate((order_quad+1)*(order_quad+1), 2*(order_quad+1)*(order_quad+1));
    RhX.Reallocate((order_quad+1)*(order_quad+1), 2*(order_quad+1)*(order_quad+1));
    RhY.Reallocate((order_quad+1)*(order_quad+1), 2*(order_quad+1)*(order_quad+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order_quad; k++)
          Ch2.AddInteraction(NumFct2D(i, j), (order_quad+1)*i + k, Phi1D(j, k));
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        for (int k = 0; k <= order_quad; k++)
          Ch1.AddInteraction((order_quad+1)*i + j, NumQuad2D(k, j), Phi1D(i, k));
    
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
        for (int k = 0; k <= order_quad; k++)
          {
            Rh.AddInteraction(NumQuad2D(i, j), 2*NumQuad2D(k, j), lob_quad.GradPhi(i, k));
            Rh.AddInteraction(NumQuad2D(i, j), 2*NumQuad2D(i, k)+1, lob_quad.GradPhi(j, k));
            RhX.AddInteraction(NumQuad2D(i, j), 2*NumQuad2D(k, j), lob_quad.GradPhi(i, k));
            RhY.AddInteraction(NumQuad2D(i, j), 2*NumQuad2D(i, k)+1, lob_quad.GradPhi(j, k));
          }
    
    Copy(Ch1, ch1_loc); Copy(Ch2, ch2_loc); Copy(Rh, rh_loc);
    Copy(RhX, rh_locX); Copy(RhY, rh_locY);
    
    if (order_geom != order)
      {
        Phi1D.Reallocate(order+1, order_geom+1);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order_geom; j++)
            Phi1D(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));

        Ch1.Clear();
        Ch1.Reallocate((order+1)*(order_geom+1), (order+1)*(order+1));
        
        Ch2.Clear();
        Ch2.Reallocate((order_geom+1)*(order_geom+1), (order_geom+1)*(order+1));

        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            for (int k = 0; k <= order; k++)
              Ch2.AddInteraction(this->GetNumNodes2D(i, j), (order+1)*i + k, Phi1D(k, j));
        
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              Ch1.AddInteraction((order+1)*i + j, NumFct2D(k, j), Phi1D(k, i));
        
        Copy(Ch1, ch1_node); Copy(Ch2, ch2_node);
      }
    
    Fb_geom.ComputeCoefficientTransformation();
  }
  
  
  void QuadrangleGauss::ConstructMassMatrix()
  {
    Matrix<Real_wp> Phi1D, dPhi1D;
    Phi1D.Reallocate(order+1, this->Points1D().GetM());
    for (int i = 0; i <= order; i++)
      for (int k = 0; k < this->Points1D().GetM(); k++)
        Phi1D(i, k) = lob_basis.EvaluatePhi(i, this->Points1D(k));
       
    // ch1_loc, ch2_loc, rh_loc
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2;
    Ch1.Reallocate((order+1)*(order_quad+1), (order_quad+1)*(order_quad+1));
    Ch2.Reallocate((order+1)*(order+1), (order+1)*(order_quad+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order_quad; k++)
          Ch2.AddInteraction(NumFct2D(i, j), (order_quad+1)*i + k, Phi1D(j, k));
    
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        for (int k = 0; k <= order_quad; k++)
          Ch1.AddInteraction((order_quad+1)*i + j, NumQuad2D(k, j), Phi1D(i, k));
        
    Copy(Ch1, ch1_loc); Copy(Ch2, ch2_loc);
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
  void QuadrangleGauss::ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
						       Vector2& res, int num_loc) const
  {
    if (GL_Gquad.GetM() > 0)
      {
        res.Reallocate(nb_dof_loc);
        res.Fill(0.0);
        for (int j = 0; j <= order; j++)
          {
            int num_dof = EdgesDof(j, num_loc);
            res(num_dof) = 0;
            for (int i = 0; i <= order_quad; i++)
              res(num_dof) += feval(i)*GL_Gquad(j, i);
          }

        return;
      }
    
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }
  
  
  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleGauss::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    if (order == this->GetGeometryOrder())
      Unode = Un;
    else
      {
        Unode.Reallocate(this->GetNbPointsNodalElt());
        Vector2 y((order+1)*(this->GetGeometryOrder()+1));
        Mlt(ch1_node, Un, y);
        Mlt(ch2_node, y, Unode);
      }
  }

  
  //! computation of u on nodal points of an edge
  /*!
    \param[in] Un dof components of u
    \param[out] Unode values of u on quadrature points of the edge
    \param[in] num_loc local position of the edge in the element
  */
  template<class Vector1, class Vector2>
  void QuadrangleGauss
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(this->GetGeometryOrder()+1);
    if (this->GetGeometryOrder() == order)
      for (int j = 0; j <= this->GetGeometryOrder(); j++)
        {
          int node = this->GetNodalNumber(num_loc, j);
          Unode(j) = Un(node);
        }
    else
      {
        for (int j = 0; j <= this->GetGeometryOrder(); j++)
          {
            Unode(j) = 0;
            for (int k = 0; k <= order; k++)
              {
                int node = EdgesDof(k, num_loc);
                Unode(j) += GL_Geom(k, j)*Un(node);
              }
          }
      }
  }
  

  //! computation of grad u on nodal points of an edge
  /*!
    \param[in] Un dof components of u
    \param[out] Unode gradient of u on quadrature points of the edge
    \param[in] num_loc local position of the edge in the element
  */  
  template<class Vector1, class Vector2>
  void QuadrangleGauss
  ::ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Vector2 Udof(2*(order+1));    
    switch (num_loc)
      {
      case 0 :
        for (int i = 0; i <= order; i++)
          {
            Udof(2*i) = 0; Udof(2*i+1) = 0;
            for (int j = 0; j <= order; j++)
              {
                Udof(2*i) += lob_basis.GradPhi(j, i)*Un(NumFct2D(j, 0));
                Udof(2*i+1) += lob_basis.GradPhi(j, 0)*Un(NumFct2D(i, j));
              }
          }
        break;
      case 1 :
        for (int i = 0; i <= order; i++)
          {
            Udof(2*i) = 0; Udof(2*i+1) = 0;
            for (int j = 0; j <= order; j++)
              {
                Udof(2*i) += lob_basis.GradPhi(j, order)*Un(NumFct2D(j, i));
                Udof(2*i+1) += lob_basis.GradPhi(j, i)*Un(NumFct2D(order, j));
              }
          }
        break;
      case 2 :
        for (int i = 0; i <= order; i++)
          {
            Udof(2*(order-i)) = 0; Udof(2*(order-i)+1) = 0;
            for (int j = 0; j <= order; j++)
              {
                Udof(2*(order-i)) += lob_basis.GradPhi(j, i)*Un(NumFct2D(j, order));
                Udof(2*(order-i)+1) += lob_basis.GradPhi(j, order)*Un(NumFct2D(i, j));
              }
          }
        break;
      case 3 :
        for (int i = 0; i <= order; i++)
          {
            Udof(2*(order-i)) = 0; Udof(2*(order-i)+1) = 0;
            for (int j = 0; j <= order; j++)
              {
                Udof(2*(order-i)) += lob_basis.GradPhi(j, 0)*Un(NumFct2D(j, i));
                Udof(2*(order-i)+1) += lob_basis.GradPhi(j, i)*Un(NumFct2D(0, j));
              }
          }
        break;
      }
    
    Unode.Reallocate(2*(this->GetGeometryOrder()+1));
    if (this->GetGeometryOrder() == order)
      {
        for (int j = 0; j < 2*this->GetGeometryOrder()+2; j++)
          Unode(j) = Udof(j);
      }
    else
      {
        for (int j = 0; j <= this->GetGeometryOrder(); j++)
          {
            Unode(2*j) = 0; Unode(2*j+1) = 0;
            for (int k = 0; k <= order; k++)
              {
                Unode(2*j) += GL_Geom(k, j)*Udof(2*k);
                Unode(2*j+1) += GL_Geom(k, j)*Udof(2*k+1);
              }
          }
      }
  }
  

  //! x is overwritten by M^-1 x where M is the mass matrix
  template<class Vector1>
  void QuadrangleGauss::SolveMassMatrixGen(Vector1& x) const
  {
    SolveCholeskyGen(SeldonNoTrans, x);
    SolveCholeskyGen(SeldonTrans, x);
  }
  

  //! x is overwritten by L^-1 x or L^-T where M = L L^t is the mass matrix  
  template<class Vector1>
  void QuadrangleGauss::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& x) const
  {
    if (order != order_quad)
      return ElementReference<Dimension2, 1>::SolveCholesky(TransA, x);
    
    Vector1 y(nb_dof_loc); typename Vector1::value_type vloc;
    if (TransA.Trans())
      {
        for (int i = 0; i < y.GetM(); i++)
          x(i) *= invSqrtWeights2d(i);
        
        for (int i = 0; i <= order; i++)
          {
            for (int j = 0; j <= order; j++)
              {
                vloc = 0;
                for (int k = 0; k <= order; k++)
                  vloc += invGL_G(k, i)*x(NumFct2D(k, j));
                
                y(NumFct2D(i, j)) = vloc;
              }
          }
        
        for (int i = 0; i <= order; i++)
          {
            for (int j = 0; j <= order; j++)
              {
                vloc = 0;
                for (int k = 0; k <= order; k++)
                  vloc += invGL_G(k, i)*y(NumFct2D(j, k));
                
                x(NumFct2D(j, i)) = vloc;
              }
          }
      }
    else
      {
        for (int i = 0; i <= order; i++)
          {
            for (int j = 0; j <= order; j++)
              {
                vloc = 0;
                for (int k = 0; k <= order; k++)
                  vloc += invGL_G(i, k)*x(NumFct2D(j, k));
                
                y(NumFct2D(j, i)) = vloc;
              }
          }
                
        for (int i = 0; i <= order; i++)
          {
            for (int j = 0; j <= order; j++)
              {
                vloc = 0;
                for (int k = 0; k <= order; k++)
                  vloc += invGL_G(i, k)*y(NumFct2D(k, j));
                
                x(NumFct2D(i, j)) = vloc;
              }
          } 
        
        for (int i = 0; i < y.GetM(); i++)
          x(i) *= invSqrtWeights2d(i);
      }
  }
  
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector2>
  void QuadrangleGauss::MltMassMatrixGen(Vector2& Vh) const
  {
    Vector2 x(nb_points_quadrature_inside);
    ApplyChTransposeGen(Vh, x);
    for (int i = 0; i < this->WeightsND().GetM(); i++)
      x(i) *= this->WeightsND(i);
    
    ApplyChGen(x, Vh);    
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
  void QuadrangleGauss::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 x(ch1_loc.GetM());
    Mlt(ch1_loc, Uh, x);
    Mlt(ch2_loc, x, Vh);
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
  void QuadrangleGauss::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector2 x(ch1_loc.GetM());
    Mlt(SeldonTrans, ch2_loc, Uh, x);
    Mlt(SeldonTrans, ch1_loc, x, Vh);
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
  void QuadrangleGauss::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
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
  void QuadrangleGauss::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Vector1 Uquad(nb_points_quadrature_inside);
    ApplyChTransposeGen(Uh, Uquad);
    Mlt(SeldonTrans, rh_loc, Uquad, Vh);
  }

  
  //! Integration against derivatives of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Ux Ux_i = \int_K f d/dx(\varphi_i)  dx
    \param[out] Uy Uy_i = \int_K f d/dy(\varphi_i)  dx
    This operation is equivalent to a matrix vector product
    Ux = Rh^1 Uh, Uy = Rh^2 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
  */  
  template<class Vector1, class Vector2>
  void QuadrangleGauss::ApplyRhSplitGen(const Vector1& Uh, Vector2& Ux, Vector2& Uy) const
  {
    Vector2 U_quad(nb_points_quadrature_inside);
    Mlt(rh_locX, Uh, U_quad); ApplyChGen(U_quad, Ux);
    Mlt(rh_locY, Uh, U_quad); ApplyChGen(U_quad, Uy);
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
  void QuadrangleGauss::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
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
  void QuadrangleGauss::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  }

  
  //! computation of u on quadrature points of an edge
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>
  void QuadrangleGauss::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    int num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= order_quad; j++)
	      Vh(j) += GL_G(i, j)*Uh(num_dof);
	  }
      }
    else
      {
	const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	for (int i = 0; i <= order; i++)
	  {
	    int num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= r; j++)
	      Vh(j) += Ah(i, j)*Uh(num_dof);
	  }
      }
  }
  
  
  //! integration against basis functions on a face
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc face number
    \param[in] Uquad vector containing \omega_k f(\xi_k) 
    \param[out] Vh res_i = res_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the edge
  */
  template<class T0, class Vector1, class Vector2>
  void QuadrangleGauss
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    int num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= order_quad; j++)
	      Vh(num_dof) += alpha*GL_G(i, j)*Uh(j);
	  }
      }
    else
      {
	const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	for (int i = 0; i <= order; i++)
	  {
	    int num_dof = EdgesDof(i, num_loc);
	    for (int j = 0; j <= r; j++)
	      Vh(num_dof) += alpha*Ah(i, j)*Uh(j);
	  }
      }
  }
  
  
  //! computation of gradient of u on quadrature points of an edge
  /*!
    \param[in] num_loc face number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the edge
    \param[in] r order of quadrature rule of the edge
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
  */
  template<class Vector1, class Vector2>
  void QuadrangleGauss
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vector2 Uquad(2*(order+1));
    Uquad.Fill(0);
    switch (num_loc)
      {
      case 0 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*i) += Uh(NumFct2D(k, 0)) * lob_basis.GradPhi(k, i);
                Uquad(2*i+1) += Uh(NumFct2D(i, k)) * lob_basis.GradPhi(k, 0);
              }            
          }
        break;
      case 1 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*i) += Uh(NumFct2D(k, i)) * lob_basis.GradPhi(k, order);
                Uquad(2*i+1) += Uh(NumFct2D(order, k)) * lob_basis.GradPhi(k, i);
              }
          }
        break;
      case 2 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*(order-i)) += Uh(NumFct2D(k, order)) * lob_basis.GradPhi(k, i);
                Uquad(2*(order-i)+1) += Uh(NumFct2D(i, k)) * lob_basis.GradPhi(k, order);
              }
          }
        break;
      case 3 :
        for (int i = 0; i <= order; i++)
          {
            for (int k = 0; k <= order; k++)
              {
                Uquad(2*(order-i)) += Uh(NumFct2D(k, i)) * lob_basis.GradPhi(k, 0);
                Uquad(2*(order-i)+1) += Uh(NumFct2D(0, k)) * lob_basis.GradPhi(k, i);
              }
          }
        break;
      }
    
    Vh.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
	  {
	    for (int j = 0; j <= order_quad; j++)
	      {
                Vh(2*j) += GL_G(i, j)*Uquad(2*i);
                Vh(2*j+1) += GL_G(i, j)*Uquad(2*i+1);
              }
	  }
      }
    else
      {
	const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	for (int i = 0; i <= order; i++)
	  {
	    for (int j = 0; j <= r; j++)
	      {
                Vh(2*j) += Ah(i, j)*Uquad(2*i);
                Vh(2*j+1) += Ah(i, j)*Uquad(2*i+1);
              }
	  }
      }
  }
  
  
  //! integration against gradient of basis functions on an edge
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc edge number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is the quadrature rule of edge num_loc
   */
  template<class T0, class Vector1, class Vector2>
  void QuadrangleGauss
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    Vector2 Uquad(2*(order+1));
    Uquad.Fill(0);
    if ((r == 0) || (r == order_quad))
      {
	for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order_quad; j++)
            {
              Uquad(2*i) += GL_G(i, j)*Uh(2*j);
              Uquad(2*i+1) += GL_G(i, j)*Uh(2*j+1);
            }
      }
    else
      {
	const Matrix<Real_wp>& Ah = ProjOperatorOrder(r);
	for (int i = 0; i <= order; i++)
          for (int j = 0; j <= r; j++)
            {
              Uquad(2*i) += Ah(i, j)*Uh(2*j);
              Uquad(2*i+1) += Ah(i, j)*Uh(2*j+1);
            }
      }

    Mlt(alpha, Uquad);
    
    typedef typename Vector2::value_type Complexe;
    switch (num_loc)
      {
      case 0 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*i);
            Complexe fy = Uquad(2*i+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, 0)) += fx*lob_basis.GradPhi(k, i);
                Vh(NumFct2D(i, k)) += fy*lob_basis.GradPhi(k, 0);
              }
          }
        break;
      case 1 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*i);
            Complexe fy = Uquad(2*i+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, i)) += fx*lob_basis.GradPhi(k, order);
                Vh(NumFct2D(order, k)) += fy*lob_basis.GradPhi(k, i);
              }
          }
        break;
      case 2 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*(order-i));
            Complexe fy = Uquad(2*(order-i)+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, order)) += fx*lob_basis.GradPhi(k, i);
                Vh(NumFct2D(i, k)) += fy*lob_basis.GradPhi(k, order);
              }
          }
        break;
      case 3 :
        for (int i = 0; i <= order; i++)
          {
            Complexe fx = Uquad(2*(order-i));
            Complexe fy = Uquad(2*(order-i)+1);
            for (int k = 0; k <= order; k++)
              {
                Vh(NumFct2D(k, i)) += fx*lob_basis.GradPhi(k, 0);
                Vh(NumFct2D(0, k)) += fy*lob_basis.GradPhi(k, i);
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
  void QuadrangleGauss::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    T zero; SetComplexZero(zero);
    if (mass == zero)
      return;
    
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	{
	  int i = NumFct2D(i1, i2);
	  for (int j1 = 0; j1 <= order; j1++)
	    for (int j2 = 0; j2 <= order; j2++)
	      {
		int j = NumFct2D(j1, j2);
		A.AddInteraction(m+i, n+j, mass*mass_matrix1D(i1, j1)*mass_matrix1D(i2, j2));
	      }
	}
  }
  
  
  //! we add constant elementary matrix to the matrix A
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + mass M + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p + \sum_p  C(p, q) S^{p,q}
    where M is the mass matrix equal to
    M_{i, j} = \int_K \phi_j \phi_i dx
    R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
    S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void QuadrangleGauss::
  AddConstantElemMatrixGen(int m, int n, const T& mass,
			   const TinyMatrix<T, Prop, 2, 2>& C,
			   const TinyVector<T, 2>& D, 
			   const TinyVector<T, 2>& E,
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    T vloc;
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
        {
          val.Zero();
	  int i = NumFct2D(i1, i2);
	  
          for (int j1 = 0; j1 <= order; j1++)
            for (int j2 = 0; j2 <= order; j2++)
              {
                int j = NumFct2D(j1, j2);
		SetComplexZero(vloc);
		
		if (!null_term(0))
		  vloc += mass*mass_matrix1D(i1, j1)*mass_matrix1D(i2, j2);
		
		if (!null_term(1))
		  vloc += C(0, 0)*stiffness_matrix1D(i1, j1)*mass_matrix1D(i2, j2)
                    + C(1, 1)*stiffness_matrix1D(i2, j2)*mass_matrix1D(i1, j1)
                    + C(0, 1)*gradient_matrix1D(i1, j1)*gradient_matrix1D(j2, i2)
                    + C(1, 0)*gradient_matrix1D(j1, i1)*gradient_matrix1D(i2, j2);
		
		if (!null_term(2))
		  vloc += D(0)*gradient_matrix1D(i1, j1)*mass_matrix1D(i2, j2)
		    + D(1)*gradient_matrix1D(i2, j2)*mass_matrix1D(i1, j1);

		if (!null_term(3))
		  vloc += E(0)*gradient_matrix1D(j1, i1)*mass_matrix1D(i2, j2)
		    + E(1)*gradient_matrix1D(j2, i2)*mass_matrix1D(i1, j1);
		
		val(j) = vloc;
              }
	  
	  A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
        }
  }


  //! adds variable elementary matrix to mat
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + R + Rt + S
    M_ij = \int mass \varphi_j \varphi_i
    R_ij = \int D \cdot \nabla \varphi_i \varphi_j
    Rt_ij = \int E \cdot \nabla \varphi_j \varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
   */
  template<class T, class Prop>
  void QuadrangleGauss
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
			     const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			     const Vector<TinyVector<T, 2> >& D,
			     const Vector<TinyVector<T, 2> >& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension2, 1>::
      AddVariableElemMatrixOpt(off_row, off_col, mass, C, D, E, null_term, mat);
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
  void QuadrangleGauss::AddVariableMassMatrixGen(int off_row, int off_col,
						 const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension2, 1>::AddVariableMassMatrixOpt(off_row, off_col, A, mat);
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */
  void QuadrangleGauss::ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const
  {
    int r = order;
    phi.Reallocate((r+1)*(r+1));
    VectReal_wp phix(r+1), phiy(r+1);
    lob_basis.ComputeValuesPhiRef(pointloc(0), phix);
    lob_basis.ComputeValuesPhiRef(pointloc(1), phiy);
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	phi(NumFct2D(i, j)) = phix(i)*phiy(j);
    
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on pointloc
  */
  void QuadrangleGauss::ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const
  {
    int r = order;
    grad_phi.Reallocate((r+1)*(r+1));
    VectReal_wp phix(r+1), phiy(r+1);
    VectReal_wp dphix(r+1), dphiy(r+1);
    for (int i = 0; i <= r; i++)
      {
	phix(i) = lob_basis.EvaluatePhi(i, pointloc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, pointloc(1));
	dphix(i) = lob_basis.EvaluatePhiGrad(i, pointloc(0));
	dphiy(i) = lob_basis.EvaluatePhiGrad(i, pointloc(1));
      }
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	grad_phi(NumFct2D(i, j)).Init(dphix(i)*phiy(j), phix(i)*dphiy(j));

  }
  
  
  //! retrieves values of a single basis function on all quadrature points
  void QuadrangleGauss::GetValueSinglePhiQuadrature(int k, VectReal_wp& phi) const
  {
    int k1 = CoordinateDofs(k, 0);
    int k2 = CoordinateDofs(k, 1);
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
	phi(NumQuad2D(i, j)) = GL_G(k1, i)*GL_G(k2, j);
  }


  //! retrieves values and gradients of a single basis function on all quadrature points
  void QuadrangleGauss
  ::GetGradientSinglePhiQuadrature(int k, VectReal_wp& phi, VectR2& grad_phi) const
  {
    int k1 = CoordinateDofs(k, 0);
    int k2 = CoordinateDofs(k, 1);
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
	{
	  int n = NumQuad2D(i, j);
	  phi(n) = GL_G(k1, i)*GL_G(k2, j);
	  grad_phi(n)(0) = dGL_G(k1, i)*GL_G(k2, j);
	  grad_phi(n)(1) = GL_G(k1, i)*dGL_G(k2, j);
	}
  }
  
  
  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void QuadrangleGauss::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad(k, 0);
        int k2 = CoordinateQuad(k, 1);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            phi(NumFct2D(i, j)) = GL_G(i, k1)*GL_G(j, k2);
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/(order_quad+1);
        int k2 = (k-nb_points_quadrature_inside)%(order_quad+1);
        switch (num_loc)
          {
          case 0:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(i, 0)) = GL_G(i, k2);
            }
            break;
          case 1:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(order, i)) = GL_G(i, k2);
            }
            break;
          case 2:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(order-i, order)) = GL_G(i, k2);
            }
            break;
          case 3:
            {
              for (int i = 0; i <= order; i++)
                phi(NumFct2D(0, order-i)) = GL_G(i, k2);
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
  void QuadrangleGauss::GetGradientPhiOnQuadraturePoint(int k, VectR2& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    FillZero(grad_phi);
    if (k < nb_points_quadrature_inside)
      {
        int k1 = CoordinateQuad(k, 0);
        int k2 = CoordinateQuad(k, 1);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            grad_phi(NumFct2D(i, j)).Init(dGL_G(i, k1)*GL_G(j, k2), GL_G(i, k1)*dGL_G(j, k2));
      }
    else
      {
        int num_loc = (k-nb_points_quadrature_inside)/(order_quad+1);
        int k2 = (k-nb_points_quadrature_inside)%(order_quad+1);
        switch (num_loc)
          {
          case 0:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumFct2D(i, 0))(0) = dGL_G(i, k2);
                  for (int j = 0; j <= order; j++)
                    grad_phi(NumFct2D(i, j))(1) = lob_basis.GradPhi(j, 0)*GL_G(i, k2);
                }
            }
            break;
          case 1:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumFct2D(order, i))(1) = dGL_G(i, k2);
                  for (int j = 0; j <= order; j++)
                    grad_phi(NumFct2D(j, i))(0) = lob_basis.GradPhi(j, order)*GL_G(i, k2);
                }
            }
            break;
          case 2:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumFct2D(i, order))(0) = dGL_G(i, order_quad-k2);
                  for (int j = 0; j <= order; j++)
                    grad_phi(NumFct2D(i, j))(1)
                      = lob_basis.GradPhi(j, order)*GL_G(i, order_quad-k2);
                }
            }
            break;
          case 3:
            {
              for (int i = 0; i <= order; i++)
                {
                  grad_phi(NumFct2D(0, i))(1) = dGL_G(i, order_quad-k2);
                  for (int j = 0; j <= order; j++)
                    grad_phi(NumFct2D(j, i))(0) = lob_basis.GradPhi(j, 0)*GL_G(i, order_quad-k2);
                }
            }
            break;
          }
      }
  }
  
}
  
#define MONTJOIE_FILE_QUADRANGLE_GAUSS_CXX
#endif
