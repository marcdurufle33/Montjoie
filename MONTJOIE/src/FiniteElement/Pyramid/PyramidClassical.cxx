#ifndef MONTJOIE_FILE_PYRAMID_CLASSICAL_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidClassical::PyramidClassical() : PyramidReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    use_quadrature_for_rh = true;
    
    type_quad_basis = Globatto<Real_wp>::QUADRATURE_GAUSS;
  }
  

  size_t PyramidClassical::GetMemorySize() const
  {
    size_t taille = PyramidReference<1>::GetMemorySize();
    taille += InvWeightFct.GetMemorySize();
    taille += InverseBasisVDM.GetMemorySize();
    taille += NumDofs2D_quad.GetMemorySize();
    taille += NumFct3D.GetMemorySize();
    for (int k = 0; k < 5; k++)
      taille += GradientPhi_Boundary(k).GetMemorySize();
    
    taille += lob_basis.GetMemorySize() + gauss_basis.GetMemorySize();
    taille += Seldon::GetMemorySize(ProjOperatorTriSurfOrder);
    taille += DerivDxtildeDx.GetMemorySize() + DerivDxtildeDz.GetMemorySize() +
      DerivDytildeDy.GetMemorySize() + DerivDytildeDz.GetMemorySize();
    taille += rh_loc.GetMemorySize();
    taille += GL_G.GetMemorySize();
    return taille;
  }
  
  
  //! constructing finite element
  void PyramidClassical::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                                                int rtri, int r_quad,
                                                int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    PyramidReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                                rtri, r_quad, type_surf_tri, type_surf_quad);
    
    ConstructFunctions();
    
    ConstructElementaryMatrix(*this);
    this->Fb_geom.ComputeCoefficientTransformation();
    
    lob_basis.ConstructQuadrature(order, lob_basis.QUADRATURE_LOBATTO);
    gauss_basis.ConstructQuadrature(order, lob_basis.QUADRATURE_GAUSS);
        
    ConstructStiffnessMatrix();
    
    GradientPhi_Boundary(0).Reallocate(nb_dof_loc, nb_dof_quad);
    R3 pt_loc; VectR3 grad_phi;
    const VectR2& points_dof2d_quad = this->PointsDof2D_quad();
    for (int k = 0; k < nb_dof_quad; k++)
      {
        this->GetLocalCoordOnBoundary(0, points_dof2d_quad(k), pt_loc);
        this->ComputeGradientPhiRef(pt_loc, grad_phi);
        for (int i = 0; i < nb_dof_loc; i++)
          GradientPhi_Boundary(0)(i, k) = grad_phi(i);            
      }
    
    //this->mass_matrix.WriteText("mat_mass_reg"+to_str(r)+".dat");
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    /*switch (r)
      {
      case 1:
	elt_r1.ConstructFiniteElement(*this);
	break;
      case 2:
	elt_r2.ConstructFiniteElement(*this);
	break;
      case 3:
	elt_r3.ConstructFiniteElement(*this);
	break;
      case 4:
	elt_r4.ConstructFiniteElement(*this);
	break;
      case 5:
	elt_r5.ConstructFiniteElement(*this);
	break;
        }*/
#endif
    
  }
  
#ifdef MONTJOIE_WITH_META_PROGRAMMING
  bool PyramidClassical::IsOrderImplemented(int r)
  {
    return (r <= 5);
  }
#endif
  
  
  //! constructing stiffness matrix
  void PyramidClassical::ConstructStiffnessMatrix()
  {
    Globatto<Real_wp> lob_quad, lob_z;
    lob_quad.AffectPoints(gauss_basis.Points());
    lob_z.AffectPoints(points1d_z);
    
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
    lob_z.ComputeGradPhi(1e3*epsilon_machine);

    GL_G.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        GL_G(i,j) = lob_basis.EvaluatePhi(i, gauss_basis.Points(j));
        
    DerivDxtildeDx.Reallocate(this->PointsND().GetM());
    DerivDxtildeDz.Reallocate(this->PointsND().GetM());
    DerivDytildeDy.Reallocate(this->PointsND().GetM());
    DerivDytildeDz.Reallocate(this->PointsND().GetM());
    
    Matrix<Real_wp, General, ArrayRowSparse> Rh;
    Rh.Reallocate(nb_points_quadrature_inside, 3*nb_points_quadrature_inside);
    for (int i = 0; i <= order_quad; i++)
      for (int j = 0; j <= order_quad; j++)
        for (int k = 0; k <= order_quad; k++)
          {
            int node = (order_quad+1)*((order_quad+1)*i + j) + k;
            for (int m = 0; m <= order_quad; m++)
              {
                int nx = (order_quad+1)*((order_quad+1)*m + j) + k;
                Rh.AddInteraction(node, 3*nx, lob_quad.GradPhi(i, m));
                
                int ny = (order_quad+1)*((order_quad+1)*i + m) + k;
                Rh.AddInteraction(node, 3*ny+1, lob_quad.GradPhi(j, m));
                
                int nz = (order_quad+1)*((order_quad+1)*i + j) + m;
                Rh.AddInteraction(node, 3*nz+2, lob_z.GradPhi(k, m));
              }
            
	    Real_wp x = this->PointsND()(node)(0), y = this->PointsND()(node)(1), z = this->PointsND()(node)(2);
	    DerivDxtildeDx(node) = 0.5/(1.0-z);
	    DerivDxtildeDz(node) = 0.5*x/square(1.0-z);
	    DerivDytildeDy(node) = 0.5/(1.0-z);
	    DerivDytildeDz(node) = 0.5*y/square(1.0-z);
          }
    
    Copy(Rh, rh_loc);
  }
  
  
  //! construction of basis functions
  void PyramidClassical::ConstructFunctions()
  {
    nb_dof_loc = this->Fb_geom.GetNbPointsNodalElt();
    nb_dof_tri = (order+1)*(order+2)/2;
    nb_dof_quad = (order+1)*(order+1);
    nb_dof_boundaries = 3*order*order+2;

    if (order != this->GetGeometryOrder())
      {
	this->elt_geom.dof_equal_nodal = false;
        int r = order;
        Matrix<int> NumDofs2D_tri, CoordinateDofs; Array3D<int> NumDofs3D;
        // construction of nodes numbering
        MeshNumbering<Dimension2>::
          ConstructTriangularNumbering(r, NumDofs2D_tri, CoordinateDofs);
        
        MeshNumbering<Dimension2>::
          ConstructQuadrilateralNumbering(r, NumDofs2D_quad, CoordinateDofs);

        MeshNumbering<Dimension3>::ConstructPyramidalNumbering(r, NumDofs3D, CoordinateDofs);
        
        FacesDof.Reallocate((r+1)*(r+1), 5);
	FacesDof.Fill(-1);
        // dofs on faces
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            FacesDof(NumDofs2D_quad(i, j), 0) = NumDofs3D(i, j, 0);
        
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            {
              int node = NumDofs2D_tri(i,j);
              FacesDof(node,1) = NumDofs3D(i, 0, j);
              FacesDof(node,2) = NumDofs3D(r-j, i, j);
              FacesDof(node,3) = NumDofs3D(i, r-j, j);
              FacesDof(node,4) = NumDofs3D(0, i, j);
            }
        
	VectReal_wp points_dof1d;
	VectR2 points_dof2d_tri, points_dof2d_quad;
	VectR3 points_dof3d;
	PyramidGeomReference::
	  ConstructLobattoPoints(r, points_dof1d, points_dof2d_tri,
				 points_dof2d_quad, NumDofs2D_quad, points_dof3d);
        
	this->SetPointsDof1D(points_dof1d);
	this->SetPointsDof2D_tri(points_dof2d_tri);
	this->SetPointsDof2D_quad(points_dof2d_quad);
	this->SetPointsDofND(points_dof3d);
	
	const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
	const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
	
        nb_dof_loc = points_dof3d.GetM();
        
        InvWeightFct.Reallocate(nb_dof_loc);
        NumFct3D.Reallocate(r+1, r+1, r+1); NumFct3D.Fill(-1);
        int nb = 0;
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            for (int k = 0; k <= r-max(i,j); k++)
              {
                NumFct3D(i, j, k) = nb;
                InvWeightFct(nb) = CoefLegendre(i)*CoefLegendre(j)*CoefEvenJacobi(max(i,j), k);
                nb++;
              }
        
        InverseBasisVDM.Reallocate(nb_dof_loc, nb_dof_loc);
        InverseBasisVDM.Fill(Real_wp(0));
            
        // evaluation of orthogonal polynomials
        VectReal_wp phi;
        for (int k = 0; k < nb_dof_loc; k++)
          {
            Fb_geom.ComputeValuesPhiOrthoRef(r, NumFct3D, InvWeightFct, points_dof3d(k), phi);
            for (int j = 0; j < nb_dof_loc; j++)
              InverseBasisVDM(j, k) = phi(j);
          }
        
        GetInverse(InverseBasisVDM);
      }
    else
      {
	this->elt_geom.dof_equal_nodal = true;
        FacesDof = this->Fb_geom.GetNodalNumber();
        this->SetPointsDof2D_tri(this->PointsNodal2D_tri());
        this->SetPointsDof2D_quad(this->PointsNodal2D_quad());
        this->SetPointsDofND(this->PointsNodalND());
	this->SetPointsDof1D(this->PointsNodal1D());
        
        Matrix<int> CoordinateDofs;
        MeshNumbering<Dimension2>::
          ConstructQuadrilateralNumbering(order, NumDofs2D_quad, CoordinateDofs);

      }

    this->nb_points_dof_inside = nb_dof_loc;
    this->num_dof_points_surf.Reallocate(5);
    for (int n = 0; n < 5; n++)
      {
	if (this->elt_geom.IsLocalFaceQuadrilateral(n))
	  this->num_dof_points_surf(n).Reallocate(nb_dof_quad);
	else
	  this->num_dof_points_surf(n).Reallocate(nb_dof_tri);
	
	for (int i = 0; i < this->num_dof_points_surf(n).GetM(); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }
    
    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    int order_geom = this->GetGeometryOrder();
    function_basis_tri = new TriangleClassical();
    function_basis_quad = new QuadrangleGauss();
    function_basis_tri->ConstructFiniteElement(order, order, order_quad);    
    function_basis_quad->ConstructFiniteElement(order, order_geom, order_quad, type_quad_basis);
    element_tri_surf = function_basis_tri;
    element_quad_surf = function_basis_quad;
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
  void PyramidClassical::ApplyRhQuadratureGen(const Vector1& Uh, Vector2& Vh) const
  {
    // applying chaine rule to get derivatives on symmetric pyramid
    Vector1 dUh(Uh.GetM());
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	dUh(3*i) = DerivDxtildeDx(i)*Uh(3*i) + DerivDxtildeDz(i)*Uh(3*i+2);
	dUh(3*i+1) = Uh(3*i+1)*DerivDytildeDy(i) + Uh(3*i+2)*DerivDytildeDz(i);
	dUh(3*i+2) = Uh(3*i+2);
      }
    
    Mlt(rh_loc, dUh, Vh);
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
  void PyramidClassical::ApplyRhQuadratureTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, rh_loc, Uh, Vh);

    // applying chaine rule to get derivatives on symmetric pyramid
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	Vh(3*i+2) += Vh(3*i+1)*DerivDytildeDz(i) + Vh(3*i)*DerivDxtildeDz(i);
	Vh(3*i+1) *= DerivDytildeDy(i);
	Vh(3*i) *= DerivDxtildeDx(i);
      }
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
  void PyramidClassical::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order))
      {
	if (num_loc == 0)
	  {
	    Vector1 Udof(nb_dof_quad);
	    for (int i = 0; i < nb_dof_quad; i++)
	      Udof(i) = Uh(FacesDof(i, 0));
	    
	    function_basis_quad->ApplyChTranspose(Udof, Vh);
	  }
	else
	  {
	    Vector1 Udof(nb_dof_tri);
	    for (int i = 0; i < nb_dof_tri; i++)
	      {
                int num_dof = FacesDof(i, num_loc);
                Udof(i) = Uh(num_dof);
              }
            
	    function_basis_tri->ApplyChTranspose(Udof, Vh);
	  }
      }
    else
      {
	Vh.Fill(0);
	typename Vector1::value_type vloc;
	if (num_loc == 0)
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	    for (int i = 0; i < nb_dof_quad; i++)
	      {
		vloc = Uh(FacesDof(i, 0));
		for (int j = 0; j < Ah.GetN(); j++)
		  Vh(j) += Ah(i, j)*vloc;
	      }
	  }
	else
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorTriSurfOrder(r);
	    for (int i = 0; i < nb_dof_tri; i++)
	      {
                int num_dof = FacesDof(i, num_loc);
                vloc = Uh(num_dof);
                
		for (int j = 0; j < Ah.GetN(); j++)
		  Vh(j) += Ah(i, j)*vloc;
	      }
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
  void PyramidClassical
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order))
      {
	if (num_loc == 0)
	  {
	    Vector1 Udof(nb_dof_quad);
	    function_basis_quad->ApplyCh(Uh, Udof);
	    
	    for (int i = 0; i < nb_dof_quad; i++)
	      Vh(FacesDof(i, num_loc)) += alpha*Udof(i);
	  }
	else
	  {
	    Vector1 Udof(nb_dof_tri);
	    function_basis_tri->ApplyCh(Uh, Udof);
	    
	    for (int i = 0; i < nb_dof_tri; i++)
	      {
                int num_dof = FacesDof(i, num_loc);
                Vh(num_dof) += alpha*Udof(i);
              }
	  }
      }
    else
      {
	typename Vector1::value_type vloc;
	
	if (num_loc == 0)
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	    for (int i = 0; i < nb_dof_quad; i++)
	      {
		vloc = 0;
		for (int j = 0; j < Ah.GetN(); j++)
		  vloc += Ah(i, j)*Uh(j);
		
		Vh(FacesDof(i, num_loc)) += alpha*vloc;
	      }
	  }
	else
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorTriSurfOrder(r);
	    for (int i = 0; i < nb_dof_tri; i++)
	      {
		vloc = 0;
		for (int j = 0; j < Ah.GetN(); j++)
		  vloc += Ah(i, j)*Uh(j);
		
		int num_dof = FacesDof(i, num_loc);
                Vh(num_dof) += alpha*vloc;
	      }
	  }
      }	
  }

  
  //! returns the projector to finite element of order r for the face num_loc
  const Matrix<Real_wp>& PyramidClassical::GetProjectionOperatorFace(int num_loc, int r) const
  {
    if (num_loc == 0)
      return ProjOperatorQuadOrder(r);
    
    return ProjOperatorTriOrder(r);
  }
  
  
  //! computation of gradient of u on quadrature points of a boundary
  /*!
    \param[in] num_loc boundary number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = grad \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class Vector1, class Vector2>
  void PyramidClassical
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order_quad))
      ElementReference<Dimension3, 1>::ApplyNablaShTranspose(num_loc, Uh, Vh, r);
    else
      {
	Vh.Fill(0);
        typedef typename Vector2::value_type T;
        TinyVector<T, 3> vloc; Real_wp val;
	const Matrix<Real_wp>& Ah = GetProjectionOperatorFace(num_loc, r);
	if (num_loc == 0)
          {
            for (int i = 0; i < Ah.GetM(); i++)
              {
                vloc.Fill(0);
                for (int j = 0; j < nb_dof_loc; j++)
                  Add(Uh(j), GradientPhi_Boundary(num_loc)(j, i), vloc);
                
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    val = Ah(i, j);
                    Vh(3*j) += val*vloc(0);
                    Vh(3*j+1) += val*vloc(1);
                    Vh(3*j+2) += val*vloc(2);
                  }
              }
          }
        else
          {
            
            for (int i = 0; i < nb_dof_loc; i++)
              {
                int offset = 3*(num_loc-1)*nb_dof_loc + 3*i;
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    Vh(3*j) += Ah(offset, j)*Uh(i);
                    Vh(3*j+1) += Ah(offset+1, j)*Uh(i);
                    Vh(3*j+2) += Ah(offset+2, j)*Uh(i);
                  }
              }
          }
      }
  }
   
    
  //! integration against gradient of basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k) 
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f grad(\psi_i) ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule on the face
   */
  template<class T0, class Vector1, class Vector2>
  void PyramidClassical
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    if ((r == 0)||(r == order_quad))
      ElementReference<Dimension3, 1>::ApplyNablaSh(alpha, num_loc, Uh, Vh, r);
    else
      {
        typedef typename Vector2::value_type T;
        TinyVector<T, 3> vloc; Real_wp val;
	const Matrix<Real_wp>& Ah = GetProjectionOperatorFace(num_loc, r);
        if (num_loc == 0)
          {
            for (int i = 0; i < Ah.GetM(); i++)
              {
                vloc.Fill(0);
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    val = Ah(i, j);
                    vloc(0) += val*Uh(3*j);
                    vloc(1) += val*Uh(3*j+1);
                    vloc(2) += val*Uh(3*j+2);
                  }
                
                Mlt(alpha, vloc);
                for (int j = 0; j < nb_dof_loc; j++)
                  Vh(j) += DotProd(vloc, GradientPhi_Boundary(num_loc)(j, i));
                
              }
          }
        else
          {
            
            for (int i = 0; i < nb_dof_loc; i++)
              {
                int offset = 3*(num_loc-1)*nb_dof_loc + 3*i;
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    Vh(i) += alpha*Ah(offset, j)*Uh(3*j);
                    Vh(i) += alpha*Ah(offset+1, j)*Uh(3*j+1);
                    Vh(i) += alpha*Ah(offset+2, j)*Uh(3*j+2);
                  }
              }
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
   */  
  template<class Vector1,class Vector2>
  void PyramidClassical::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    if (num_loc == 0)
      {
	const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
        res.Fill(0);
        Vector2 w1((order_quad+1)*(order_quad+1)); FillZero(w1);
        Vector2 w2(w1.GetM()); FillZero(w2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order_quad; j++)
            for (int k = 0; k <= order_quad; k++)
              w1(i*(order_quad+1) + j) += GL_G(i, k)*feval(NumQuad2D(k, j));
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order_quad; k++)
              w2(NumDofs2D_quad(i, j)) += GL_G(j, k)*w1(i*(order_quad+1) + k);
        
        for (int j = 0; j < (order+1)*(order+1); j++)
          {
            int num_dof = FacesDof(j, num_loc);
            res(num_dof) = w2(j);
          }
      }
    else
      ComputeIntegralSurfaceRef(feval, res, num_loc);
  }


  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] x local point where functions are evaluated
    \param[out] phi values of basis functions on x
  */
  void PyramidClassical::ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const
  {
    int order_geom = this->GetGeometryOrder();
    if (order != order_geom)
      {
        VectReal_wp psi; phi.Reallocate(nb_dof_loc); phi.Fill(0);
        Fb_geom.ComputeValuesPhiOrthoRef(order, NumFct3D, InvWeightFct, x, psi);
        Mlt(InverseBasisVDM, psi, phi);
      }
    else
      ComputeValuesPhiNodalRef(x, phi);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] x local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on x
  */
  void PyramidClassical::ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const
  {
    int order_geom = this->GetGeometryOrder();
    if (order != order_geom)
      {
        VectReal_wp psi, phi; VectR3 grad_psi; 
        psi.Reallocate(nb_dof_loc); psi.Fill(0);
        phi.Reallocate(nb_dof_loc); phi.Fill(0);
        grad_phi.Reallocate(nb_dof_loc);
        Fb_geom.ComputeGradientPhiOrthoRef(order, NumFct3D, InvWeightFct, x, grad_psi);
        for (int i = 0; i < psi.GetM(); i++)
          psi(i) = grad_psi(i)(0);
	
        Mlt(InverseBasisVDM, psi, phi);
        
        for (int i = 0; i < psi.GetM(); i++)
          {
            grad_phi(i)(0) = phi(i);
            psi(i) = grad_psi(i)(1);
          }
        
        Mlt(InverseBasisVDM, psi, phi);
        for (int i = 0; i < psi.GetM(); i++)
          {
            grad_phi(i)(1) = phi(i);
            psi(i) = grad_psi(i)(2);
          }
        
        Mlt(InverseBasisVDM, psi, phi);
        for (int i = 0; i < psi.GetM(); i++)
          grad_phi(i)(2) = phi(i);
        
      }
    else
      ComputeGradientPhiNodalRef(x, grad_phi);
  }
  

  //! computation of projection between finite element of different orders  
  void PyramidClassical::
  ComputeTriangularInterpolationProjector(const IVect& order_elt,
                                          const ElementReference<Dimension2, 1>& Fb,
                                          const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi; VectR3 grad_phi; R3 pt_loc;
    ProjOperatorTriOrder.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorTriOrder(r).Reallocate(12*nb_dof_loc, Pts(i).GetM());
	for (int num_loc = 1; num_loc < 5; num_loc++)
          for (int j = 0; j < Pts(i).GetM(); j++)
            {
              GetLocalCoordOnBoundary(num_loc, Pts(i)(j), pt_loc);
              ComputeGradientPhiRef(pt_loc, grad_phi);
              for (int k = 0; k < nb_dof_loc; k++)
                for (int m = 0; m < 3; m++)
                  ProjOperatorTriOrder(r)(3*k + m + 3*(num_loc-1)*nb_dof_loc, j) = grad_phi(k)(m);
            }
      }

    ProjOperatorTriSurfOrder.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorTriSurfOrder(r).Reallocate((order+1)*(order+2)/2, Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
            Fb.ComputeValuesPhiRef(Pts(i)(j), phi);
	    for (int k = 0; k < Fb.GetNbDof(); k++)
              ProjOperatorTriSurfOrder(r)(k, j) = phi(k);
          }
      }
  }
  
  
  //! computation of projection between finite element of different orders  
  void PyramidClassical::
  ComputeQuadrangularInterpolationProjector(const IVect& order_elt,
                                            const ElementReference<Dimension2, 1>& Fb,
                                            const Vector<VectR2>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi;
    ProjOperatorQuadOrder.Reallocate(rmax+1);
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();    
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate((order+1)*(order+1), Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumQuad2D(i1, i2);
		  Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0);
		  ProjOperatorQuadOrder(r)(k, j) = lob_basis.EvaluatePhi(i1, xt)*lob_basis.EvaluatePhi(i2, yt);
		}	    
	  }
      }
  }

  
  //! displays details about class PyramidClassical
  ostream& operator <<(ostream& out, const PyramidClassical& e)
  {
    out<<static_cast<const PyramidReference<1>&>(e);
    out<<"Number dof over the pyramid "<<e.nb_dof_loc<<endl;
    return out;
  }


  //! default constructor
  PyramidLobatto::PyramidLobatto() : PyramidClassical()
  {
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  
  
  //! constructing finite element
  void PyramidLobatto::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                         int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    type_quad_basis = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    PyramidClassical::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                             r, r, TriangleQuadrature::QUADRATURE_GAUSS,
                                             Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
  }
  
} // end namespace
  
#define MONTJOIE_FILE_PYRAMID_CLASSICAL_CXX
#endif
