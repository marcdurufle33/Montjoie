#ifndef MONTJOIE_FILE_QUADRANGLE_RADAU_CXX

namespace Montjoie
{
  
  //! default constructor
  QuadrangleRadau::QuadrangleRadau() : QuadrangleLobatto()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    radau_points = false;
    mass_lumping = true;
    diagonal_mass = true;
  }
  

  size_t QuadrangleRadau::GetMemorySize() const
  {
    size_t taille = QuadrangleLobatto::GetMemorySize();
    taille += lob_x.GetMemorySize();
    return taille;
  }

  
  //! constructing finite element
  void QuadrangleRadau::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
					       int rsurf, int type_surf)
  {
    if (rquad == 0)
      rquad = r;

    if (r != rquad)
      {
        cout << "A quadrature order different from order not implemented for QuadrangleRadau : " << r << " not equal to " << rquad << endl;
        abort();
      }

    QuadrangleLobatto::ConstructFiniteElement(r, rgeom, r, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    if (radau_points)
      {
        VectReal_wp points_x, weights_x;
        VectReal_wp points_gauss, weights_gauss;
        bool regular_element = false;
        if (regular_element)
          ComputeGaussRadauJacobi(points_x, weights_x, r,
                                  Real_wp(0), Real_wp(1), true);
        else
          ComputeGaussRadauJacobi(points_x, weights_x, r,
                                  Real_wp(0), Real_wp(0), true);
        
        lob_x.AffectPoints(points_x);
	
        //DISP(this->points1d);
        //DISP(points_x); DISP(weights_x);
        //DISP(NumQuad2D); DISP(NumFct2D);
        //DISP(this->points_dof2d);
        VectR2 points2d; VectReal_wp weights2d;
        points2d.Reallocate((r+1)*(r+1) + 4*(r+1));
        weights2d.Reallocate((r+1)*(r+1) + 4*(r+1));
        weights2d.Fill(0);
        
        // cas avec points de Gauss
        /* ComputeGaussLegendre(points_gauss, weights_gauss, r);
         for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            {
              // this->points2d(NumQuad2D(i, j)).Init(points_gauss(i), points_gauss(j));
              // this->weights2d(NumQuad2D(i, j)) = weights_gauss(i)*weights_gauss(j);
              } */
        
        // cas avec points de Gauss-Radau (=> mass lumping)
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            {
              points2d(NumQuad2D(i, j)).Init(points_x(i), this->Points1D(j));
              if (regular_element)
                weights2d(NumQuad2D(i, j)) = weights_x(i)*this->Weights1D(j)/points_x(i);
              else
                weights2d(NumQuad2D(i, j)) = weights_x(i)*this->Weights1D(j);
            }
        
        //DISP(this->points1d);
        int nb = (r+1)*(r+1);
        int rsurf = r;
	this->num_quad_points_surf.Reallocate(4);
	for (int n = 0; n < 4; n++)
	  this->num_quad_points_surf(n).Reallocate(rsurf+1);
	
        // Points of integration for the first edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(0)(i) = nb;
            points2d(nb++).Init(this->Points1D(i), 0.0);
          }
        
        // Points of integration for the second edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(1)(i) = nb;
            points2d(nb++).Init(1.0, this->Points1D(i));
          }
        
        // Points of integration for the third edge
        for (int i = 0; i <= rsurf; i++)
          {
	    this->num_quad_points_surf(2)(i) = nb;
            points2d(nb++).Init(1.0-this->Points1D(i), 1.0);
          }
        
        // Points of integration for the fourth edge
        for (int i = 0; i <= rsurf; i++)
          {
            this->num_quad_points_surf(3)(i) = nb;
	    points2d(nb++).Init(0.0, 1.0-this->Points1D(i));
          }
        
        this->points2d_all = points2d;
	this->SetPointsND(points2d);
	this->SetWeightsND(weights2d);
        //DISP(this->points2d); DISP(this->weights2d);
        
	this->num_dof_points_surf = this->num_quad_points_surf;
	this->SetPointsDofND(points2d_all);
	this->elt_geom.quadrature_equal_nodal = false;
	this->elt_geom.dof_equal_nodal = false;
	this->elt_geom.dof_equal_quadrature = true;
	
        this->Fb_geom.ComputeCoefficientTransformation();
	
        ConstructElementaryMatrix(*this);
        ConstructMassMatrix();
      }
  }
  
  
  //! constructing mass matrix
  void QuadrangleRadau::ConstructMassMatrix()
  {
    int order_geom = this->GetGeometryOrder();
    Matrix<Real_wp> Phi1D_x, Phi1D_y;
    Matrix<Real_wp, General, ArrayRowSparse> Ch1, Ch2;
    Phi1D_x.Reallocate(order+1, order_geom+1);
    Phi1D_y.Reallocate(order+1, order_geom+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        {
          Phi1D_x(i, j) = lob_x.EvaluatePhi(i, this->PointsNodal1D(j));
          Phi1D_y(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
        }
    
    Ch1.Clear();
    Ch1.Reallocate((order+1)*(order_geom+1), (order+1)*(order+1));
    
    Ch2.Clear();
    if (order != order_geom)
      {
        Ch2.Reallocate((order_geom+1)*(order_geom+1), (order_geom+1)*(order+1));
        
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order_geom; j++)
            for (int k = 0; k <= order; k++)
              Ch2.AddInteraction(this->GetNumNodes2D(i, j), (order+1)*i + k, Phi1D_y(k, j));
        
        Copy(Ch2, ch2_node);
        for (int i = 0; i <= order_geom; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
              Ch1.AddInteraction((order+1)*i + j, NumFct2D(k, j), Phi1D_x(k, i));
      }
    else
      for (int i = 0; i <= order_geom; i++)
        for (int j = 0; j <= order; j++)
          for (int k = 0; k <= order; k++)
            Ch1.AddInteraction(this->GetNumNodes2D(i, j), NumFct2D(k, j), Phi1D_x(k, i));
    
    Copy(Ch1, ch1_node); 
    //DISP(Ch1); DISP(Ch2);
  }    
  
  
  //! computation of stiffness matrix R_h
  /*!
    this matrix is used to compute gradients
   */
  void QuadrangleRadau::ConstructStiffnessMatrix()
  {    
  }
  

  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<class Vector1, class Vector2>
  void QuadrangleRadau::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    if (!radau_points)
      return QuadrangleLobatto::ComputeNodalValuesRef(Un, Unode);

    int order_geom = this->GetGeometryOrder();
    if (order == order_geom)
      {
        Mlt(ch1_node, Un, Unode);
      }
    else
      {
        Unode.Reallocate(this->GetNbPointsNodalElt());
        Vector2 y((order+1)*(order_geom+1));
        Mlt(ch1_node, Un, y);
        Mlt(ch2_node, y, Unode);
      }
  }
  

  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral against each basis function
    \param[in] num_loc local edge number
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1, class Vector2>
  void QuadrangleRadau::ComputeGaussIntegralSurfaceGen(const Vector1 & feval,
						       Vector2& res, int num_loc) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
    else
      QuadrangleLobatto::ComputeGaussIntegralSurfaceGen(feval, res, num_loc);
  }
  
  
  //! Uh is overwritten by M^-1 Uh where M is the mass matrix
  template<class Vector1>
  void QuadrangleRadau::SolveMassMatrixGen(Vector1& Uh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::SolveMassMatrix(Uh);
        
    for (int i = 0; i < nb_dof_loc; i++)
      Uh(i) *= invWeights2d(i);
  }
   
  
  //! Uh is overwritten by L^-1 Uh or L^-T Uh where M = L L^t is the mass matrix  
  template<class Vector1>
  void QuadrangleRadau::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& Uh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::SolveCholesky(TransA, Uh);

    for (int i = 0; i < nb_dof_loc; i++)
      Uh(i) *= invSqrtWeights2d(i);
  }
  
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  template<class Vector2>
  void QuadrangleRadau::MltMassMatrixGen(Vector2& Vh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::MltMassMatrix(Vh);
      
    const VectReal_wp& weights2d = this->WeightsND();
    for (int i = 0; i < nb_dof_loc; i++)
      Vh(i) *= weights2d(i);
  }
  
  
  //! Integration against basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Vh res_i = \int_K f \varphi_i dx
   */
  template<class Vector1, class Vector2>
  void QuadrangleRadau::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::ApplyChGen(Uh, Vh);

    // Ch = Id (quadrature points = interpolation points)
    Vh = Uh;
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
  void QuadrangleRadau::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::ApplyChTransposeGen(Uh, Vh);
 
    // Ch = Id (quadrature points = interpolation points)
    Vh = Uh;
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
  void QuadrangleRadau::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::ApplyRhGen(Uh, Vh);

    Mlt(rh_loc, Uh, Vh);
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
  void QuadrangleRadau::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (radau_points)
      return ElementReference<Dimension2, 1>::ApplyRhTransposeGen(Uh, Vh);
    
    Mlt(SeldonTrans, rh_loc, Uh, Vh);
  } 
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void QuadrangleRadau::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyConstantRhGen(Uh, Vh);
    else
      Mlt(rh_weight, Uh, Vh);
  }
  
  
  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    then Vx = Rh^1 Uh, Vy = Rh^2 Uh
   */
  template<class Vector1, class Vector2>
  void QuadrangleRadau::ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyConstantRhSplitGen(Uh, Vx, Vy);
    else
      {
        Mlt(rh_weightX, Uh, Vx);
        Mlt(rh_weightY, Uh, Vy);
      }
  }
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void QuadrangleRadau::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyConstantRhTransposeGen(Uh, Vh);
    else
      Mlt(SeldonTrans, rh_weight, Uh, Vh);
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
  void QuadrangleRadau
  ::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyShTransposeGen(num_loc, Uh, Vh);
    else
      QuadrangleLobatto::ApplyShTranspose(num_loc, Uh, Vh);
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
  void QuadrangleRadau
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyShGen(alpha, num_loc, Uh, Vh, r);
    else
      QuadrangleLobatto::ApplyShGen(alpha, num_loc, Uh, Vh, r);
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
  void QuadrangleRadau
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
    else
      QuadrangleLobatto::ApplyNablaShTransposeGen(num_loc, Uh, Vh, r);
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
  void QuadrangleRadau
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
    else
      QuadrangleLobatto::ApplyNablaShGen(alpha, num_loc, Uh, Vh, r);
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
  void QuadrangleRadau::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::AddConstantMassMatrixGen(m, n, mass, A);
    else
      QuadrangleLobatto::AddConstantMassMatrixGen(m, n, mass, A);
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
  void QuadrangleRadau::
  AddConstantElemMatrixGen(int m, int n, const T& mass,
			   const TinyMatrix<T, Prop, 2, 2>& C,
			   const TinyVector<T, 2>& D, 
			   const TinyVector<T, 2>& E,
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
    else
      QuadrangleLobatto::AddConstantElemMatrixGen(m, n, mass, C, D, E, null_term, A);
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
  void QuadrangleRadau
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
			     const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			     const Vector<TinyVector<T, 2> >& D,
			     const Vector<TinyVector<T, 2> >& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::
	AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, mat);
    else
      QuadrangleLobatto::
	AddVariableElemMatrixGen(off_row, off_col, mass, C, D, E, null_term, mat);
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
  void QuadrangleRadau::
  AddVariableMassMatrixGen(int off_row, int off_col,
			   const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::AddVariableMassMatrixGen(off_row, off_col, A, mat);
    else
      QuadrangleLobatto::AddVariableMassMatrixGen(off_row, off_col, A, mat);
  }


  //! computation of gradient on quadrature points of the edge of the square
  /*!
    \param[in] pts reference points
    \param[in] mat jacobian matrices
    \param[in] feval components of u (by dofs)
    \param[out] res value of gradient of u on quadrature points of the edge
    \param[in] vars given problem
    \param[in] FaceBasis given finite element 
    \param[in] n element number
    \param[in] num_loc local number of the edge
   */
  template<class Vector1, class Vector2>
  void QuadrangleRadau::
  ComputeGradientBoundaryGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::ComputeGradientBoundaryGen(feval, res, num_loc);
    else
      QuadrangleLobatto::ComputeGradientBoundaryGen(feval, res, num_loc);
  }
    

  //! retrieves values of a single basis functions on all quadrature points
  void QuadrangleRadau::GetValueSinglePhiQuadrature(int k, VectReal_wp& phi) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::GetValueSinglePhiQuadrature(k, phi);
    else
      QuadrangleLobatto::GetValueSinglePhiQuadrature(k, phi);
  }


  //! retrieves values of a single basis functions on all quadrature points
  void QuadrangleRadau
  ::GetGradientSinglePhiQuadrature(int k, VectReal_wp& phi, VectR2& grad_phi) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::GetGradientSinglePhiQuadrature(k, phi, grad_phi);
    else
      QuadrangleLobatto::GetGradientSinglePhiQuadrature(k, phi, grad_phi);
  }

      
  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */
  void QuadrangleRadau::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::GetValuePhiOnQuadraturePoint(k, phi);
    else
      QuadrangleLobatto::GetValuePhiOnQuadraturePoint(k, phi);
  }
  
  
  //! Retrieving gradient of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] grad_phi gradient of basis functions on quadrature point k
  */
  void QuadrangleRadau::GetGradientPhiOnQuadraturePoint(int k, VectR2& grad_phi) const
  {
    if (radau_points)
      ElementReference<Dimension2, 1>::GetGradientPhiOnQuadraturePoint(k, grad_phi);
    else
      QuadrangleLobatto::GetGradientPhiOnQuadraturePoint(k, grad_phi);    
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */
  void QuadrangleRadau::ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const
  {
    if (!radau_points)
      return QuadrangleLobatto::ComputeValuesPhiRef(pointloc, phi);
    
    int r = order;
    phi.Reallocate((r+1)*(r+1));
    VectReal_wp phix(r+1), phiy(r+1);
    for (int i = 0; i <= r; i++)
      {
	phix(i) = lob_x.EvaluatePhi(i, pointloc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, pointloc(1));
      }
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	phi(NumFct2D(i, j)) = phix(i)*phiy(j);

  }
  

  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on pointloc
  */
  void QuadrangleRadau::ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const
  {
    if (!radau_points)
      return QuadrangleLobatto::ComputeGradientPhiRef(pointloc, grad_phi);

    int r = order;
    grad_phi.Reallocate((r+1)*(r+1));
    VectReal_wp phix(r+1), phiy(r+1);
    VectReal_wp dphix(r+1), dphiy(r+1);
    for (int i = 0; i <= r; i++)
      {
	phix(i) = lob_x.EvaluatePhi(i, pointloc(0));
	phiy(i) = lob_basis.EvaluatePhi(i, pointloc(1));
	dphix(i) = lob_x.EvaluatePhiGrad(i, pointloc(0));
	dphiy(i) = lob_basis.EvaluatePhiGrad(i, pointloc(1));
      }
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	grad_phi(NumFct2D(i, j)).Init(dphix(i)*phiy(j), phix(i)*dphiy(j));

  }


  //! displays informations about class QuadrangleRadau
  ostream& operator <<(ostream& out, const QuadrangleRadau& e)
  {
    out<<"Number of degrees of freedom on a face "<<e.nb_dof_loc<<endl;
    return out;
  }
    
}
  
#define MONTJOIE_FILE_QUADRANGLE_RADAU_CXX
#endif
